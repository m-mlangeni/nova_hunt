"""
nova_lightcurves.py -- pull BlackGEM/MeerLICHT and ATLAS light curves for the new nova
candidates and plot them together, grouped by host galaxy.

Input : new_nova_candidates.csv  (rank, kind, galaxy, host_D, ra, dec, M_peak, first_mjd, last_mjd)
Output: lc_plots/<galaxy>_#<rank>.png, nova_lightcurves.html (all plots grouped by galaxy),
        lc_data/bg_detections.csv, lc_data/bg_limits.csv, lc_data/atlas_#<rank>.txt

BlackGEM/MeerLICHT (BigQuery, same tables as query_lightcurve.py):
  * every detection within 3" of each candidate (any CLASS_REAL), and
  * every image whose footprint covers the position, with its transient depth T_LMAG
    (non-detections are drawn as upper limits).
ATLAS forced photometry (https://fallingstar-data.com/forcedphot/):
  * needs a (free) ATLAS forced-photometry account. Give either
      export ATLAS_TOKEN=...                         (token from the website), or
      export ATLAS_USER=...  ATLAS_PASSWORD=...      (the script fetches a token)
  * full ATLAS history (from MJD 57000) on difference images; results cached in lc_data/.

Usage
  python nova_lightcurves.py --dry-run      # BigQuery cost estimate only
  python nova_lightcurves.py                # BG/ML + ATLAS + plots
  python nova_lightcurves.py --no-atlas     # BG/ML only
"""

import argparse, io, os, re, sys, time, base64, html
import numpy as np
import pandas as pd

PROJECT = "blackgem-databases-access"
FALLBACK_PROJECTS = ["blackgem-transient-db", "blackgem"]
DET_TABLE = "`blackgem-transient-db.blackgem_transients_raw.historic_detections`"
IMG_TABLE = "`blackgem-transient-db.blackgem_transients_raw.images`"
ARCSEC_M = 30.887
MATCH_ARCSEC = 3.0
FIELD_HALF_DEG = 0.64          # BG/ML field half-width
ATLAS_URL = "https://fallingstar-data.com/forcedphot"
ATLAS_MJD_MIN = 57000.0
OUT_PLOTS, OUT_DATA = "lc_plots", "lc_data"


# ---------------------------------------------------------------- BigQuery
def bq_client(project, sql):
    from google.cloud import bigquery
    from google.api_core.exceptions import Forbidden
    for proj in [project] + [p for p in FALLBACK_PROJECTS if p != project]:
        client = bigquery.Client(project=proj)
        try:
            dj = client.query(sql, job_config=bigquery.QueryJobConfig(dry_run=True, use_query_cache=False))
            return client, dj.total_bytes_processed / 1e9, proj
        except Forbidden:
            print(f"  no permission in project '{proj}', trying the next one")
    sys.exit("No BigQuery project accepted the job.")


def run_bq(sql, project, dry_run, label):
    client, gb, proj = bq_client(project, sql)
    print(f"{label}: ~{gb:.2f} GB (~${gb / 1e3 * 5:.2f}) in project '{proj}'")
    if dry_run:
        return None
    return client.query(sql).result().to_dataframe()


def targets_sql(C):
    return "[" + ",".join(f"STRUCT({int(r['rank'])} AS rank, {r.ra:.7f} AS ra, {r.dec:.7f} AS dec)"
                          for _, r in C.iterrows()) + "]"


def bg_detections_sql(C):
    return f"""
    SELECT t.rank, i.MJD_OBS AS mjd, i.FILTER AS filt, d.IMAGE_ID AS image_id,
           d.RA_PSF_D AS ra, d.DEC_PSF_D AS dec, d.MAG_PSF_D AS mag, d.MAGERR_PSF_D AS magerr,
           d.FNU_PSF_D AS fnu, d.FNUERR_PSF_D AS fnuerr, d.SNR_ZOGY AS snr, d.CLASS_REAL AS rb,
           CASE WHEN SUBSTR(REGEXP_EXTRACT(i.FILENAME, r'[^/]+$'), 1, 3) IN ('BG2','BG3','BG4')
                  THEN SUBSTR(REGEXP_EXTRACT(i.FILENAME, r'[^/]+$'), 1, 3)
                WHEN SUBSTR(REGEXP_EXTRACT(i.FILENAME, r'[^/]+$'), 1, 2) = 'ML' THEN 'ML1' END AS telescop
    FROM UNNEST({targets_sql(C)}) AS t
    JOIN {DET_TABLE} d ON ST_DWITHIN(d.GEOG_POINT, ST_GEOGPOINT(t.ra, t.dec), {MATCH_ARCSEC * ARCSEC_M:.1f})
    JOIN {IMG_TABLE} i ON d.IMAGE_ID = i.IMAGE_ID
    """


def bg_images_sql(C):
    return f"""
    SELECT t.rank, i.IMAGE_ID AS image_id, i.MJD_OBS AS mjd, i.FILTER AS filt,
           i.T_LMAG AS t_lmag, i.TQC_FLAG AS tqc, i.RA_CNTR AS ra_c, i.DEC_CNTR AS dec_c
    FROM UNNEST({targets_sql(C)}) AS t
    JOIN {IMG_TABLE} i
      ON ST_DWITHIN(ST_GEOGPOINT(i.RA_CNTR, i.DEC_CNTR), ST_GEOGPOINT(t.ra, t.dec), {0.9 * 3600 * ARCSEC_M:.0f})
    """


def get_bg(C, project, dry_run):
    os.makedirs(OUT_DATA, exist_ok=True)
    fd, fi = f"{OUT_DATA}/bg_detections.csv", f"{OUT_DATA}/bg_images.csv"
    if os.path.exists(fd) and os.path.exists(fi) and os.path.exists(f"{OUT_DATA}/bg_targets.csv") and not dry_run:
        # reuse only if the cache was made for exactly these positions
        cached = pd.read_csv(f"{OUT_DATA}/bg_targets.csv")
        key = lambda t: set(zip(t["rank"], t.ra.round(5), t.dec.round(5)))
        if key(cached) == key(C):
            print(f"Using cached {fd} and {fi} (delete them to re-query)")
            return pd.read_csv(fd), pd.read_csv(fi)
        print("Candidate list changed since the BG/ML data were pulled: re-querying")
    det = run_bq(bg_detections_sql(C), project, dry_run, "BG/ML detections query")
    img = run_bq(bg_images_sql(C), project, dry_run, "BG/ML images query")
    if dry_run:
        return None, None
    det.to_csv(fd, index=False)
    img.to_csv(fi, index=False)
    C[["rank", "ra", "dec"]].to_csv(f"{OUT_DATA}/bg_targets.csv", index=False)
    return det, img


def bg_limits(C, det, img):
    """Images covering each position in which nothing was detected there -> upper limits."""
    img = img.merge(C[["rank", "ra", "dec"]].rename(columns={"ra": "t_ra", "dec": "t_dec"}), on="rank")
    dra = ((img.ra_c - img.t_ra + 180) % 360 - 180) * np.cos(np.radians(img.t_dec))
    img = img[(dra.abs() < FIELD_HALF_DEG) & ((img.dec_c - img.t_dec).abs() < FIELD_HALF_DEG)]
    img = img[~img.tqc.astype(str).str.lower().str.startswith("red")]
    hit = set(zip(det["rank"], det.image_id.astype(str)))
    keep = [(r, str(i)) not in hit for r, i in zip(img["rank"], img.image_id)]
    lim = img[keep].dropna(subset=["t_lmag"])
    lim.to_csv(f"{OUT_DATA}/bg_limits.csv", index=False)
    return lim


# ---------------------------------------------------------------- ATLAS
def atlas_headers():
    import requests
    tok = os.environ.get("ATLAS_TOKEN")
    if not tok:
        user, pw = os.environ.get("ATLAS_USER"), os.environ.get("ATLAS_PASSWORD")
        if not (user and pw):
            return None
        r = requests.post(f"{ATLAS_URL}/api-token-auth/", data={"username": user, "password": pw})
        if r.status_code != 200:
            print(f"ATLAS login failed: HTTP {r.status_code} {r.text[:200]}")
            return None
        tok = r.json()["token"]
        print("ATLAS: got a token (store it as ATLAS_TOKEN to skip the login next time)")
    return {"Authorization": f"Token {tok}", "Accept": "application/json"}


def atlas_one(ra, dec, headers, path):
    """Queue one forced-photometry task, wait for it, save the text result to path."""
    import requests
    task_url = None
    while not task_url:
        r = requests.post(f"{ATLAS_URL}/queue/", headers=headers,
                          data={"ra": ra, "dec": dec, "mjd_min": ATLAS_MJD_MIN, "send_email": False})
        if r.status_code == 201:
            task_url = r.json()["url"]
        elif r.status_code == 429:                       # throttled: wait as told
            msg = r.json().get("detail", "")
            s = re.findall(r"available in (\d+) seconds", msg)
            m = re.findall(r"available in (\d+) minutes", msg)
            wait = int(s[0]) if s else int(m[0]) * 60 if m else 10
            print(f"    throttled, waiting {wait}s")
            time.sleep(wait)
        else:
            print(f"    could not queue: HTTP {r.status_code} {r.text[:200]}")
            return False
    t0 = time.time()
    while True:
        r = requests.get(task_url, headers=headers)
        if r.status_code != 200:
            print(f"    status check failed: HTTP {r.status_code}")
            return False
        j = r.json()
        if j.get("finishtimestamp"):
            if j.get("error_msg"):
                print(f"    ATLAS: {j['error_msg']}")
                return False
            res = requests.get(j["result_url"], headers=headers)
            with open(path, "w") as fh:
                fh.write(res.text)
            requests.delete(task_url, headers=headers)     # tidy the user's queue
            return True
        time.sleep(5 if j.get("starttimestamp") else 10)
        if time.time() - t0 > 3600:
            print("    gave up after 1 h")
            return False


def read_atlas(path):
    """ATLAS forced phot text -> DataFrame (MJD, filter, uJy, duJy, mag5sig, ...)."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    with open(path) as fh:
        lines = [l for l in fh if l.strip()]
    if not lines or not lines[0].lstrip().startswith("#"):
        return None
    cols = lines[0].lstrip("#").split()
    rows = [l.split() for l in lines[1:]]
    a = pd.DataFrame(rows, columns=cols)
    for c in ["MJD", "m", "dm", "uJy", "duJy", "err", "chi/N", "mag5sig"]:
        if c in a:
            a[c] = pd.to_numeric(a[c], errors="coerce")
    return a.rename(columns={"MJD": "mjd", "F": "filt"})


def atlas_nightly(a):
    """Clean (err==0, sane chi/N) and combine each night's exposures per filter
    (weighted mean flux), then 3-sigma detections / 3-sigma upper limits in AB mag."""
    if a is None or a.empty:
        return None
    a = a[(a.get("err", 0) == 0) & (a.duJy > 0) & (a.duJy < 1e4)]
    if "chi/N" in a:
        a = a[a["chi/N"] < 10]
    a = a.assign(night=np.floor(a.mjd - 0.5))
    out = []
    for (n, f), g in a.groupby(["night", "filt"]):
        w = 1 / g.duJy ** 2
        flux, eflux = (g.uJy * w).sum() / w.sum(), 1 / np.sqrt(w.sum())
        out.append(dict(mjd=g.mjd.mean(), filt=f, uJy=flux, duJy=eflux, n=len(g)))
    N = pd.DataFrame(out)
    N["snr"] = N.uJy / N.duJy
    N["det"] = N.snr >= 3
    N["mag"] = np.where(N.det, 23.9 - 2.5 * np.log10(N.uJy.clip(lower=1e-3)), np.nan)
    N["magerr"] = np.where(N.det, 1.0857 / N.snr, np.nan)
    N["lim"] = 23.9 - 2.5 * np.log10(3 * N.duJy)
    return N


def get_atlas(C, enabled):
    os.makedirs(OUT_DATA, exist_ok=True)
    res = {}
    headers = atlas_headers() if enabled else None
    if enabled and headers is None:
        print("ATLAS: no ATLAS_TOKEN or ATLAS_USER/ATLAS_PASSWORD set -- using only cached ATLAS files")
    for _, r in C.iterrows():
        path = f"{OUT_DATA}/atlas_#{int(r['rank'])}.txt"
        if not os.path.exists(path) and headers is not None:
            print(f"ATLAS #{int(r['rank'])} ({r.galaxy}) RA {r.ra:.5f} Dec {r.dec:.5f} ...")
            atlas_one(r.ra, r.dec, headers, path)
        res[int(r["rank"])] = atlas_nightly(read_atlas(path))
    return res


# ---------------------------------------------------------------- plots
BGCOL = dict(u="#7b5ea7", g="#2e8b57", q="#e07b00", r="#c0392b", i="#7a4a2a", z="#333333")
ATCOL = dict(o="#d62728", c="#17becf")


def plot_one(r, det, lim, atl, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rk = int(r["rank"])
    d = det[det["rank"] == rk]
    l = lim[lim["rank"] == rk]
    a = atl.get(rk)
    t0, t1 = r.first_mjd, r.last_mjd
    fig, axs = plt.subplots(2, 1, figsize=(10, 7.2), gridspec_kw=dict(height_ratios=[1, 1.15]))
    windows = [(57200, 61500, "full history (BG/ML + ATLAS)"),
               (t0 - 80, t1 + 160, f"zoom on the outburst (MJD {t0:.0f}-{t1:.0f})")]
    for ax, (lo, hi, title) in zip(axs, windows):
        full = hi - lo > 2000
        if a is not None and len(a):
            aa = a[(a.mjd > lo) & (a.mjd < hi)]
            for f, g in aa.groupby("filt"):
                u = g[~g.det]
                ax.scatter(u.mjd, u.lim, marker="v", s=8 if full else 14, color=ATCOL.get(f, "grey"),
                           alpha=0.18, lw=0)
                dd = g[g.det]
                ax.errorbar(dd.mjd, dd.mag, dd.magerr, fmt="s", ms=3.5, color=ATCOL.get(f, "grey"),
                            alpha=0.9, lw=0.7, label=f"ATLAS {f} (nightly, ≥3σ)")
        ll = l[(l.mjd > lo) & (l.mjd < hi)]
        for f, g in ll.groupby("filt"):
            ax.scatter(g.mjd, g.t_lmag, marker="v", s=6 if full else 12, color=BGCOL.get(f, "grey"),
                       alpha=0.25, lw=0)
        dd = d[(d.mjd > lo) & (d.mjd < hi)]
        for f, g in dd.groupby("filt"):
            good, weak = g[g.rb >= 0.5], g[g.rb < 0.5]
            ax.errorbar(good.mjd, good.mag, good.magerr, fmt="o", ms=4.5, color=BGCOL.get(f, "grey"),
                        lw=0.8, label=f"BG/ML {f}")
            ax.scatter(weak.mjd, weak.mag, marker="o", s=22, facecolors="none",
                       edgecolors=BGCOL.get(f, "grey"), lw=0.8)
        ax.axvspan(t0 - 0.5, t1 + 0.5, color="gold", alpha=0.15, lw=0)
        ax.set_xlim(lo, hi)
        ys = list(d.mag.dropna().values) + [r.M_peak + 5 * np.log10(r.host_D * 1e5)]
        if a is not None and len(a):
            ad = a[a.det & (a.mjd > lo) & (a.mjd < hi) & (a.magerr < 0.3)]
            ys += list(ad.mag.values)
        ys = np.array(ys, float)
        ax.set_ylim(min(np.nanmax(ys) + 1.5, 22.5), max(np.nanmin(ys) - 0.7, 12.0))
        ax.set_xlabel("MJD", fontsize=8)
        ax.set_ylabel("AB mag", fontsize=8)
        ax.set_title(title, fontsize=9, loc="left")
        ax.tick_params(labelsize=7)
        if not full:
            h, lab = ax.get_legend_handles_labels()
            seen = dict(zip(lab, h))
            ax.legend(seen.values(), seen.keys(), fontsize=6.5, loc="lower left", ncol=3)
    dist_mod = 5 * np.log10(r.host_D * 1e6) - 5
    fig.suptitle(f"#{rk}  {r.kind}  ·  {r.galaxy} ({r.host_D:.1f} Mpc, μ={dist_mod:.2f})  ·  "
                 f"RA {r.ra:.5f}  Dec {r.dec:.5f}  ·  M_peak {r.M_peak:.1f}", fontsize=10)
    fig.text(0.01, 0.005, "filled: BG/ML CLASS_REAL≥0.5 · open: <0.5 · triangles: non-detection limits "
             "(BG/ML T_LMAG; ATLAS 3σ nightly) · shaded: candidate's detection window", fontsize=6.5)
    fig.tight_layout(rect=(0, 0.02, 1, 0.96))
    fig.savefig(path, dpi=110)
    plt.close(fig)


def summarise(r, det, atl):
    rk = int(r["rank"])
    d = det[det["rank"] == rk]
    a = atl.get(rk)
    s = {"bg_det": len(d), "bg_det_rb05": int((d.rb >= 0.5).sum()),
         "bg_nights": int(np.floor(d.mjd - 0.5).nunique()) if len(d) else 0}
    if a is not None and len(a):
        win = a[(a.mjd > r.first_mjd - 5) & (a.mjd < r.last_mjd + 30)]
        s.update(atlas_nights=len(a), atlas_nights_window=len(win),
                 atlas_det_window=int(win.det.sum()),
                 atlas_det_elsewhere=int(a[(a.mjd < r.first_mjd - 60) | (a.mjd > r.last_mjd + 200)].det.sum()))
    else:
        s.update(atlas_nights=0, atlas_nights_window=0, atlas_det_window=0, atlas_det_elsewhere=0)
    return s


def write_page(C, S):
    parts = []
    for gal, G in C.groupby("galaxy", sort=False):
        parts.append(f"<h2>{html.escape(gal)} <span class=sub>({len(G)} candidate{'s' if len(G) > 1 else ''})</span></h2>")
        for _, r in G.iterrows():
            rk = int(r["rank"])
            p = f"{OUT_PLOTS}/{gal.split(' (')[0].replace(' ', '')}_#{rk}.png"
            b64 = base64.b64encode(open(p, "rb").read()).decode() if os.path.exists(p) else ""
            s = S[rk]
            parts.append(
                f"<div class=card><div class=meta><b>#{rk} · {html.escape(r.kind)}</b><br>"
                f"RA {r.ra:.5f} · Dec {r.dec:.5f}<br>M<sub>peak</sub> {r.M_peak:.1f}<br>"
                f"BG/ML: {s['bg_det']} detections ({s['bg_det_rb05']} with RB≥0.5) on {s['bg_nights']} nights<br>"
                f"ATLAS: {s['atlas_nights']} nights; in the outburst window {s['atlas_det_window']} of "
                f"{s['atlas_nights_window']} nights ≥3σ; other epochs ≥3σ: {s['atlas_det_elsewhere']}<br>"
                f"<a target=_blank href='https://www.legacysurvey.org/viewer?ra={r.ra:.6f}&dec={r.dec:.6f}&zoom=15&mark={r.ra:.6f},{r.dec:.6f}'>Legacy viewer</a>"
                f"</div><img src='data:image/png;base64,{b64}'></div>")
    counts = C.groupby("galaxy", sort=False).size()
    table = "".join(f"<tr><td>{html.escape(g)}</td><td>{n}</td></tr>" for g, n in counts.items())
    doc = ("<!DOCTYPE html><html><head><meta charset=utf-8><title>New nova candidates: BG/ML + ATLAS light curves</title>"
           "<style>body{font-family:system-ui,sans-serif;font-size:13px;margin:24px;max-width:1300px;color:#1f2933}"
           "h2{border-bottom:1px solid #ddd;margin-top:28px}.sub{color:#6b7280;font-weight:normal}"
           ".card{display:flex;gap:14px;border-bottom:1px solid #eee;padding:10px 0}.meta{width:300px;flex:none;line-height:1.5}"
           "img{max-width:950px}table{border-collapse:collapse}td{padding:2px 12px;border-bottom:1px solid #eee}</style></head><body>"
           "<h1>New nova candidates within 5 Mpc: BlackGEM/MeerLICHT + ATLAS light curves</h1>"
           f"<table><tr><th align=left>Galaxy</th><th>Candidates</th></tr>{table}"
           f"<tr><td><b>Total</b></td><td><b>{len(C)}</b></td></tr></table>{''.join(parts)}</body></html>")
    with open("nova_lightcurves.html", "w") as fh:
        fh.write(doc)
    print("wrote nova_lightcurves.html")


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidates", default="new_nova_candidates.csv")
    ap.add_argument("--project", default=PROJECT)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--no-atlas", action="store_true")
    a = ap.parse_args()

    C = pd.read_csv(a.candidates)
    order = C.groupby("galaxy").size().sort_values(ascending=False).index
    C["galaxy"] = pd.Categorical(C.galaxy, categories=order, ordered=True)
    C = C.sort_values(["galaxy", "rank"]).reset_index(drop=True)
    C["galaxy"] = C.galaxy.astype(str)
    print("Candidates per galaxy:")
    print(C.groupby("galaxy", sort=False).size().to_string())

    det, img = get_bg(C, a.project, a.dry_run)
    if a.dry_run:
        return
    lim = bg_limits(C, det, img)
    print(f"BG/ML: {len(det)} detections, {len(lim)} non-detection limits")
    atl = get_atlas(C, not a.no_atlas)

    os.makedirs(OUT_PLOTS, exist_ok=True)
    S = {}
    for _, r in C.iterrows():
        p = f"{OUT_PLOTS}/{r.galaxy.split(' (')[0].replace(' ', '')}_#{int(r['rank'])}.png"
        plot_one(r, det, lim, atl, p)
        S[int(r["rank"])] = summarise(r, det, atl)
    pd.DataFrame(S).T.to_csv(f"{OUT_DATA}/summary.csv")
    write_page(C, S)


if __name__ == "__main__":
    main()
