"""
nova_hunt.py -- systematic search for missed novae in nearby galaxies
using the MeerLICHT/BlackGEM transient database (BigQuery) + REGALADE.

Pipeline
  1. REGALADE galaxies within DMAX Mpc in the BG/ML sky (VizieR or a local file)
  2. ONE BigQuery spatial join: all real/bogus-passing detections inside them
  3. Elliptical host association (REGALADE R1/R2/PA)
  4. Group detections into sources, build light curves, apply nova cuts
  5. "Was it there before?" check against image depth limits (T_LMAG)
  6. Ranked HTML page for visual inspection

Table / column names are taken from Pull_BG_ML/BG_Trans_Query_ID/query_lightcurve.py.

Usage
  python nova_hunt.py --dry-run                 # cost estimate only, nothing pulled
  python nova_hunt.py                           # run (galaxies fetched from VizieR)
  python nova_hunt.py --regalade my_cut.fits    # use a local REGALADE file instead
  python nova_hunt.py --dmax 12 --rb 0.6        # tweak cuts
"""

import argparse
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

# ---------------------------------------------------------------- CONFIG
CFG = dict(
    dmax_mpc=20.0,           # novae (M ~ -6..-10) are hard to see beyond ~15-20 Mpc
    dmin_mpc=0.1,            # skip Local Group dwarfs / bad distances below this
    dec_max=30.0,            # BG/ML footprint (southern)
    ellipse_scale=2.0,       # keep detections within 2x the REGALADE ellipse
    default_r1_arcsec=30.0,  # used when a galaxy has no ellipse fit
    max_r1_arcsec=1800.0,    # cap huge galaxies (e.g. LMC/SMC-like entries)
    rb_min=0.5,              # source-level cut: MEAN CLASS_REAL over its detections must exceed this
    rb_floor=0.1,            # per-detection floor in the SQL (keeps weak detections so the mean is honest)
    min_det=2,               # at least this many detections
    det_frac_min=0.5,        # detected in >= this fraction of exposures on its detected nights (0.25 = looser)
    field_half_deg=0.64,     # BG/ML field half-width (1.65 deg^2 square) for the exposure count
    m_lim=21.5,              # typical single-epoch depth, for galaxy ranking
    Mabs_bright=-10.5,       # brighter -> SN, not a nova
    Mabs_faint=-4.5,         # fainter -> unlikely nova
    quiet_days=30.0,         # must be undetected this long before first detection
    min_decline_mag=0.5,     # must fade at least this much after peak
    max_duration_days=200.0, # detected over longer -> variable / AGN
    min_nights=2,            # use --min-nights 1 to also allow single-night sources
    fade_or_within_days=30., # must fade >= min_decline_mag, OR all detections within this window
    ref_arcsec=1.5,          # a reference-image source this close ...
    ref_brighter_by=1.0,     # ... and brighter than peak+this -> quiescent counterpart = variable
    tns_arcsec=3.0,          # TNS crossmatch radius
    nucleus_arcsec=2.0,      # reject sources this close to the galaxy centre
    gaia_star_arcsec=1.5,    # reject sources this close to a Gaia star
    # (G brighter than, within arcsec) pairs for a bright-star cut. OFF by default: BG's "bright"
    # Gaia stars are G<~14, and real novae in Cen A sit 17-90" from G 11-13 stars, so any useful
    # radius also removes real transients. Artefact groups are handled by the cluster flag instead.
    bright_star=[],
    cluster_arcsec=90.0,     # >= cluster_n other candidates this close -> likely artefact group
    cluster_n=3,
    link_arcsec=2.0,         # detections closer than this = same source
    fov_radius_deg=0.9,      # BG/ML half-diagonal, for the image-depth query
)

PROJECT = "blackgem-databases-access"   # project the query jobs are billed to
FALLBACK_PROJECTS = ["blackgem-transient-db", "blackgem"]
DET_TABLE = "`blackgem-transient-db.blackgem_transients_raw.historic_detections`"
IMG_TABLE = "`blackgem-transient-db.blackgem_transients_raw.images`"
ARCSEC_M = 30.887      # metres on the BigQuery sphere per arcsec (as in query_lightcurve.py)

TELESCOP_EXPR = """CASE
        WHEN SUBSTR(REGEXP_EXTRACT(i.FILENAME, r'[^/]+$'), 1, 3) IN ('BG2','BG3','BG4')
            THEN SUBSTR(REGEXP_EXTRACT(i.FILENAME, r'[^/]+$'), 1, 3)
        WHEN SUBSTR(REGEXP_EXTRACT(i.FILENAME, r'[^/]+$'), 1, 2) = 'ML' THEN 'ML1'
        ELSE NULL END"""


# ---------------------------------------------------------------- 1. galaxies
VIZIER_MIRRORS = [
    "vizier.cds.unistra.fr",     # Strasbourg (main)
    "vizier.cfa.harvard.edu",    # Harvard
    "vizier.nao.ac.jp",          # Japan
    "vizier.hia.nrc.ca",         # Canada
    "vizier.ast.cam.ac.uk",      # Cambridge
    "vizier.china-vo.org",       # China
]
REGALADE_VIZIER_ID = "J/A+A/706/A284"   # A&A 706, A284 (2026) -- best guess, verified at runtime
def regalade_cache():
    return f"regalade_D{CFG['dmax_mpc']:g}Mpc.csv"


def _pick(cols, *cands):
    low = {c.lower(): c for c in cols}
    for c in cands:
        if c.lower() in low:
            return low[c.lower()]
    return None


def _download_strips(server, cat_id, cache):
    """Discover column names, then download only the needed columns in 15-deg RA strips
    with the distance/Dec cuts applied on the server. Retries each strip; keeps partial
    progress in <cache>.part<i>.csv so a rerun resumes where it stopped."""
    import os, time
    from astroquery.vizier import Vizier
    head = Vizier(columns=["**"], row_limit=1, vizier_server=server, timeout=60).get_catalogs(cat_id)
    tab = max(head, key=lambda t: len(t.colnames))      # the main table
    cols = tab.colnames
    print("  columns:", ", ".join(cols))
    c_ra = _pick(cols, "RAJ2000", "RA_ICRS", "ra", "_RA")
    c_de = _pick(cols, "DEJ2000", "DE_ICRS", "dec", "_DE")
    c_d = _pick(cols, "D", "Dist", "Dmpc")
    if not (c_ra and c_de and c_d):
        raise RuntimeError(f"could not identify ra/dec/D columns in {cols}")
    want = [c for c in [_pick(cols, "name", "Name"), c_ra, c_de, _pick(cols, "R1"), _pick(cols, "R2"),
                        _pick(cols, "PA"), c_d, _pick(cols, "logM", "logMstar")] if c]
    table_id = tab.meta.get("name", cat_id)
    print(f"  table {table_id}; filtering {c_d} in {CFG['dmin_mpc']}..{CFG['dmax_mpc']}, {c_de} < {CFG['dec_max']}")
    parts = []
    for k, ra0 in enumerate(range(0, 360, 15)):
        part = f"{cache}.part{k:02d}.csv"
        if os.path.exists(part):
            parts.append(pd.read_csv(part)); continue
        filt = {c_d: f"{CFG['dmin_mpc']}..{CFG['dmax_mpc']}",
                c_ra: f"{ra0}..{ra0 + 15}", c_de: f"<{CFG['dec_max']}"}
        for attempt in range(3):
            try:
                v = Vizier(columns=want, column_filters=filt, row_limit=-1,
                           vizier_server=server, timeout=300)
                res = v.query_constraints(catalog=table_id)
                df = res[0].to_pandas() if len(res) else pd.DataFrame(columns=want)
                break
            except Exception as e:
                print(f"    RA {ra0}-{ra0 + 15}: attempt {attempt + 1} failed ({type(e).__name__})")
                time.sleep(5)
        else:
            raise RuntimeError(f"RA strip {ra0}-{ra0 + 15} failed 3 times (rerun to resume)")
        df.to_csv(part, index=False)
        parts.append(df)
        print(f"    RA {ra0:3d}-{ra0 + 15:3d}: {len(df)} galaxies")
    t = pd.concat(parts, ignore_index=True)
    t = t.rename(columns={c_ra: "ra", c_de: "dec", c_d: "D"})
    t = t.drop_duplicates(subset=["ra", "dec"]).reset_index(drop=True)
    t.to_csv(cache, index=False)
    for k in range(24):
        try:
            os.remove(f"{cache}.part{k:02d}.csv")
        except OSError:
            pass
    print(f"  {len(t)} galaxies saved to {cache}")
    return t


def fetch_regalade_vizier():
    """Pull only the nearby southern galaxies from VizieR (server-side cut), trying
    several mirrors with a short timeout. Result is cached per distance limit (regalade_D<dmax>Mpc.csv)."""
    import os
    cache = regalade_cache()
    if os.path.exists(cache):
        print(f"Using cached {cache} (delete it to re-download)")
        return pd.read_csv(cache)
    from astroquery.vizier import Vizier
    errors = []
    for server in VIZIER_MIRRORS:
        print(f"Trying VizieR mirror {server} ...")
        kw = dict(vizier_server=server, timeout=20)
        try:
            cat_id = REGALADE_VIZIER_ID
            try:
                Vizier(row_limit=1, **kw).get_catalogs(cat_id)[0]
            except Exception:
                cats = Vizier(**kw).find_catalogs("REGALADE")
                if not cats:
                    raise RuntimeError("REGALADE not listed on this mirror")
                cat_id = list(cats.keys())[0]
            print("  found catalogue", cat_id)
            t = _download_strips(server, cat_id, cache)
            return t
        except Exception as e:
            errors.append(f"{server}: {type(e).__name__}: {str(e)[:120]}")
            print("  failed:", errors[-1])
    raise SystemExit(
        "\nCould not get REGALADE from any VizieR mirror:\n  " + "\n  ".join(errors) +
        "\n\nOptions:\n"
        "  1. Download it by hand (browser/VPN), e.g. from CDS 'J/A+A/706/A284', or ask the\n"
        "     BlackGEM team (Hugo Tranin / blackpearl.blackgem.org) for the D<20 Mpc subset,\n"
        "     then run:  python nova_hunt.py --regalade that_file.fits\n"
        "  2. Check whether your network blocks VizieR:  curl -I https://vizier.cfa.harvard.edu\n")


def load_galaxies(path=None):
    if path is None:
        g = fetch_regalade_vizier()
    elif path.endswith(".parquet"):
        g = pd.read_parquet(path)
    elif path.endswith((".fits", ".fit", ".fits.gz")):
        from astropy.table import Table
        g = Table.read(path).to_pandas()
    else:
        g = pd.read_csv(path)
    # normalise column names (VizieR / local files may differ in case or naming)
    alias = {"raj2000": "ra", "dej2000": "dec", "ra_icrs": "ra", "de_icrs": "dec", "de": "dec",
             "d": "D", "dist": "D", "r1": "R1", "r2": "R2", "pa": "PA",
             "logm": "logM", "logmstar": "logM", "name": "name"}
    g = g.rename(columns={c: alias[c.lower()] for c in g.columns if c.lower() in alias})
    missing = [c for c in ["ra", "dec", "D"] if c not in g]
    if missing:
        raise SystemExit(f"REGALADE table lacks columns {missing}; columns are: {list(g.columns)}")
    for c, v in [("name", ""), ("R1", np.nan), ("R2", np.nan), ("PA", np.nan), ("logM", np.nan)]:
        if c not in g:
            g[c] = v
    if "name" in g and g["name"].dtype == object:
        g["name"] = g["name"].apply(lambda s: s.decode() if isinstance(s, bytes) else s).str.strip()
    g = g[(g.D > CFG["dmin_mpc"]) & (g.D < CFG["dmax_mpc"]) & (g.dec < CFG["dec_max"])].copy()
    g["R1"] = g["R1"].fillna(CFG["default_r1_arcsec"]).clip(upper=CFG["max_r1_arcsec"])
    g["R2"] = g["R2"].fillna(g["R1"]).clip(upper=CFG["max_r1_arcsec"])
    g["PA"] = g["PA"].fillna(0.0)
    g["logM"] = g["logM"].fillna(g["logM"].median() if g["logM"].notna().any() else 9.0)
    g["dm"] = 5 * np.log10(g.D * 1e6) - 5
    # relative nova yield: ~ stellar mass x fraction of nova LF above the depth
    g["priority"] = 10 ** (g.logM - 10) * np.clip((CFG["m_lim"] - g.dm + 10.5) / 6, 0, 1)
    g["name"] = g["name"].astype(str)
    g = g.sort_values(["priority", "name", "ra"], ascending=[False, True, True],
                      kind="mergesort").reset_index(drop=True)   # deterministic order
    g["gid"] = np.arange(len(g))
    return g


# ---------------------------------------------------------------- 2. BigQuery
def _galaxy_array_sql(gal):
    """Inline the galaxy list as an UNNEST'ed array of STRUCTs (fine for <~20k galaxies)."""
    rows = [f"STRUCT({int(r.gid)} AS gid, {r.ra:.7f} AS ra, {r.dec:.7f} AS dec, "
            f"{CFG['ellipse_scale'] * r.R1 * ARCSEC_M:.1f} AS rad_m)"
            for r in gal.itertuples()]
    return "[" + ",\n".join(rows) + "]"


def detections_query(gal):
    return f"""
    SELECT
        g.gid,
        d.DETECTION_ID, d.IMAGE_ID,
        i.MJD_OBS AS mjd, i.FILTER AS filt, i.OBJECT AS field_id,
        {TELESCOP_EXPR} AS telescop,
        d.RA_PSF_D AS ra, d.DEC_PSF_D AS dec,
        d.FNU_PSF_D AS fnu, d.MAG_PSF_D AS mag, d.MAGERR_PSF_D AS magerr,
        d.SNR_ZOGY AS snr, d.CLASS_REAL AS rb,
        d.sep_near_gaia, d.mag_g_near_gaia, d.sep_near_ref, d.mag_auto_near_ref,
        d.sep_bright_gaia, d.mag_g_bright_gaia
    FROM UNNEST({_galaxy_array_sql(gal)}) AS g
    JOIN {DET_TABLE} d
      ON ST_DWITHIN(d.GEOG_POINT, ST_GEOGPOINT(g.ra, g.dec), g.rad_m)
    JOIN {IMG_TABLE} i
      ON d.IMAGE_ID = i.IMAGE_ID
    WHERE d.CLASS_REAL >= {CFG['rb_floor']}
      AND d.FNU_PSF_D > 0
    """


def run_query(sql, project, dry_run=False):
    from google.cloud import bigquery
    from google.api_core.exceptions import Forbidden
    client = bigquery.Client(project=project)
    dj = None
    for proj in [project] + [p for p in FALLBACK_PROJECTS if p != project]:
        client = bigquery.Client(project=proj)
        try:
            dj = client.query(sql, job_config=bigquery.QueryJobConfig(dry_run=True, use_query_cache=False))
            print(f"  running in project '{proj}'")
            break
        except Forbidden:
            print(f"  no permission to run jobs in project '{proj}', trying the next one")
    if dj is None:
        raise SystemExit("No BigQuery project accepted the job -- ask for bigquery.jobs.create permission.")
    gb = dj.total_bytes_processed / 1e9
    print(f"  query will scan ~{gb:.2f} GB (~${gb / 1e3 * 5:.2f} on-demand)")
    if dry_run:
        return None
    return client.query(sql).result().to_dataframe()


def images_query(sources):
    """All images whose centre is near each candidate: depth (T_LMAG), quality, centre."""
    rows = [f"STRUCT({int(s)} AS src, {r.ra:.7f} AS ra, {r.dec:.7f} AS dec)"
            for s, r in sources.iterrows()]
    return f"""
    SELECT c.src, i.IMAGE_ID AS image_id, i.MJD_OBS AS mjd, i.FILTER AS filt,
           i.T_LMAG AS t_lmag, i.TQC_FLAG AS tqc, i.RA_CNTR AS ra_c, i.DEC_CNTR AS dec_c
    FROM UNNEST([{",".join(rows)}]) AS c
    JOIN {IMG_TABLE} i
      ON ST_DWITHIN(ST_GEOGPOINT(i.RA_CNTR, i.DEC_CNTR), ST_GEOGPOINT(c.ra, c.dec),
                    {CFG['fov_radius_deg'] * 3600 * ARCSEC_M:.0f})
    """


def fetch_images(F, project, cache, chunk=3000):
    """Images covering each candidate, cached; only missing candidates are re-queried.
    Queried in chunks so the inline array stays under BigQuery's query-length limit."""
    import os
    have = pd.read_parquet(cache) if os.path.exists(cache) else None
    done = set(have.src.unique()) if have is not None else set()
    todo = F[~F.index.isin(done)]
    parts = [have] if have is not None else []
    for k in range(0, len(todo), chunk):
        print(f"  image query {k // chunk + 1}/{(len(todo) - 1) // chunk + 1} ({len(todo[k:k + chunk])} sources)")
        res = run_query(images_query(todo.iloc[k:k + chunk]), project)
        if res is not None:
            parts.append(res)
    if not parts:
        return None
    images = pd.concat(parts, ignore_index=True)
    # record queried sources with no images too, so they are not re-queried
    empty = set(todo.index) - set(images.src.unique())
    if empty:
        images = pd.concat([images, pd.DataFrame({"src": list(empty)})], ignore_index=True)
    images.to_parquet(cache)
    images = images.dropna(subset=["mjd"])
    # keep only images whose (square) footprint really contains the source
    pos = F[["ra", "dec"]].rename(columns={"ra": "s_ra", "dec": "s_dec"})
    images = images.join(pos, on="src", how="inner")
    dra = ((images.ra_c - images.s_ra + 180) % 360 - 180) * np.cos(np.radians(images.s_dec))
    inside = (dra.abs() < CFG["field_half_deg"]) & ((images.dec_c - images.s_dec).abs() < CFG["field_half_deg"])
    good_q = ~images.tqc.astype(str).str.lower().str.startswith("red")
    return images[inside & good_q]


def add_detection_fraction(F, det, images):
    """Colleague's criterion: detected in >= 1/2 (or 1/4) of the exposures taken on the nights
    it was detected. Counted PER FILTER (novae are found mostly in q, less in i/u), and the
    best filter with >= 2 exposures is used."""
    F = F.copy()
    F["det_frac"] = np.nan
    F["det_frac_filt"] = ""
    if images is None or images.empty:
        return F
    d = det[det.src.isin(F.index)][["src", "IMAGE_ID", "mjd", "filt"]].copy()
    d["IMAGE_ID"] = d.IMAGE_ID.astype(str)
    d["night"] = np.floor(d.mjd - 0.5)
    im = images[["src", "image_id", "mjd", "filt"]].copy()
    im["image_id"] = im.image_id.astype(str)
    im["night"] = np.floor(im.mjd - 0.5)
    det_nights = d[["src", "night"]].drop_duplicates()
    im = im.merge(det_nights, on=["src", "night"])
    exp_sets = im.groupby(["src", "filt"]).image_id.agg(set)
    det_sets = d.groupby(["src", "filt"]).IMAGE_ID.agg(set)
    best, bestf = {}, {}
    for (s, f), ids in det_sets.items():
        n_exp = len(ids | exp_sets.get((s, f), set()))   # an exposure it was detected in counts
        if n_exp < 2:
            continue
        frac = len(ids) / n_exp
        if frac > best.get(s, -1):
            best[s], bestf[s] = frac, f
    F["det_frac"] = pd.Series(best).reindex(F.index)
    F["det_frac_filt"] = pd.Series(bestf).reindex(F.index).fillna("")
    return F


def flag_clusters(F):
    """Candidates with several other candidates very close by are usually artefacts
    (bright-star bleeds, ghosts, crosstalk); score them down and flag them."""
    F = F.copy()
    F["n_cand_near"] = 0
    if len(F) > 1:
        c = SkyCoord(F.ra.values * u.deg, F.dec.values * u.deg)
        i, _, _, _ = c.search_around_sky(c, CFG["cluster_arcsec"] * u.arcsec)
        F["n_cand_near"] = np.bincount(i, minlength=len(F)) - 1
    F["score"] *= np.where(F.n_cand_near >= CFG["cluster_n"], 0.3, 1.0)
    return F.sort_values("score", ascending=False)


# ---------------------------------------------------------------- 3. host association
def associate(det, gal):
    """Elliptical radius of each detection w.r.t. the galaxy it was matched to; keep best host."""
    G = gal.set_index("gid").loc[det.gid]
    dra = (det.ra.values - G.ra.values)
    dra = (dra + 180) % 360 - 180
    dra *= np.cos(np.radians(G.dec.values)) * 3600
    ddec = (det.dec.values - G.dec.values) * 3600
    pa = np.radians(G.PA.values)
    x = dra * np.sin(pa) + ddec * np.cos(pa)      # along major axis (PA east of north)
    y = -dra * np.cos(pa) + ddec * np.sin(pa)     # along minor axis
    det = det.copy()
    det["r_ell"] = np.sqrt((x / G.R1.values) ** 2 + (y / np.maximum(G.R2.values, 1)) ** 2)
    det["r_nuc"] = np.hypot(dra, ddec)
    det = det[det.r_ell < CFG["ellipse_scale"]]
    det = det.sort_values("r_ell").drop_duplicates("DETECTION_ID")
    for col in ["name", "D", "dm", "logM", "priority"]:
        det["host_" + col] = gal.set_index("gid").loc[det.gid, col].values
    return det.reset_index(drop=True)


# ---------------------------------------------------------------- 4. sources + nova features
def group_sources(det):
    """Friends-of-friends linking of detections into sources (connected components)."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    c = SkyCoord(det.ra.values * u.deg, det.dec.values * u.deg)
    i, j, _, _ = c.search_around_sky(c, CFG["link_arcsec"] * u.arcsec)
    n = len(det)
    _, labels = connected_components(coo_matrix((np.ones(len(i)), (i, j)), shape=(n, n)),
                                     directed=False)
    det = det.copy()
    det["src"] = labels
    return det


def nova_features(lc):
    lc = lc.sort_values("mjd")
    good = lc[np.isfinite(lc.mag)]
    if good.empty:
        return None
    peak = good.loc[good.mag.idxmin()]
    after = good[good.mjd > peak.mjd + 0.5]
    f = dict(
        ra=lc.ra.median(), dec=lc.dec.median(), n_det=len(lc),
        n_nights=lc.mjd.round().nunique(), telescopes=",".join(sorted(lc.telescop.dropna().unique())),
        first_mjd=lc.mjd.min(), last_mjd=lc.mjd.max(),
        peak_mag=peak.mag, peak_filt=peak.filt, peak_mjd=peak.mjd,
        M_peak=peak.mag - lc.host_dm.iloc[0],
        decline=(after.mag.max() - peak.mag) if len(after) else np.nan,
        rb_med=lc.rb.median(), rb_mean=lc.rb.mean(),
        host=lc.host_name.iloc[0], host_D=lc.host_D.iloc[0], host_logM=lc.host_logM.iloc[0],
        host_priority=lc.host_priority.iloc[0], r_ell=lc.r_ell.min(), r_nuc=lc.r_nuc.median(),
        sep_gaia=lc.sep_near_gaia.median(), gaia_g=lc.mag_g_near_gaia.median(),
        sep_ref=lc.sep_near_ref.median(), ref_mag=lc.mag_auto_near_ref.median(),
        sep_bright=lc.sep_bright_gaia.median() if "sep_bright_gaia" in lc else np.nan,
        bright_g=lc.mag_g_bright_gaia.median() if "mag_g_bright_gaia" in lc else np.nan,
        t2_days=np.nan,
    )
    if len(after) >= 2 and np.ptp(after.mjd) > 0.5:
        k = np.polyfit(after.mjd - peak.mjd, after.mag - peak.mag, 1)[0]
        f["t2_days"] = 2.0 / k if k > 0 else np.nan
    return pd.Series(f)


def apply_cuts(F):
    F = F.copy()
    dur = F.last_mjd - F.first_mjd
    c = pd.DataFrame(index=F.index)
    c["Mabs"] = F.M_peak.between(CFG["Mabs_bright"], CFG["Mabs_faint"])
    c["off_nucleus"] = F.r_nuc > CFG["nucleus_arcsec"]
    c["n_det"] = F.n_det >= CFG["min_det"]
    c["mean_rb"] = F.rb_mean > CFG["rb_min"]
    c["nights"] = F.n_nights >= CFG["min_nights"]
    c["duration"] = dur < CFG["max_duration_days"]
    c["fades"] = (F.decline >= CFG["min_decline_mag"]) | (dur < CFG["fade_or_within_days"])
    c["no_gaia_star"] = ~(F.sep_gaia < CFG["gaia_star_arcsec"])
    c["no_ref_counterpart"] = ~((F.sep_ref < CFG["ref_arcsec"])
                                & (F.ref_mag < F.peak_mag + CFG["ref_brighter_by"]))
    near_bright = np.zeros(len(F), bool)
    for gmax, rad in CFG["bright_star"]:
        near_bright |= ((F.bright_g < gmax) & (F.sep_bright < rad)).values
    c["no_bright_star"] = ~near_bright
    keep = np.ones(len(F), bool)
    print("  cut                 survivors")
    for col in c:
        keep &= c[col].values
        print(f"  {col:<20s}{keep.sum():>8d}")
    s = np.ones(len(F))
    s *= np.where(F.decline >= CFG["min_decline_mag"], 1.0, 0.6)
    s *= np.where(F.t2_days.between(2, 150), 1.5, 1.0)
    s *= np.clip(F.n_nights / 4, 0.5, 1.0)
    s *= np.clip((F.rb_med - 0.4) / 0.5, 0.2, 1.0)
    s *= np.clip(F.host_priority / F.host_priority.max(), 0.05, 1) ** 0.3
    F["score"] = s
    return F[keep].sort_values("score", ascending=False)


def flag_tns(F, path):
    """Attach the nearest TNS object (name, type, discovery date) within tns_arcsec."""
    import os
    F = F.copy()
    F["tns_name"] = ""
    F["tns_type"] = ""
    F["tns_date"] = ""
    if not path or not os.path.exists(path) or F.empty:
        print("  (no TNS file found -- skipping known-transient check)")
        return F
    t = pd.read_csv(path, usecols=["name_prefix", "name", "ra", "declination", "type", "discoverydate"])
    t = t.dropna(subset=["ra", "declination"])
    ct = SkyCoord(t.ra.values * u.deg, t.declination.values * u.deg)
    cf = SkyCoord(F.ra.values * u.deg, F.dec.values * u.deg)
    idx, sep, _ = cf.match_to_catalog_sky(ct)
    ok = sep.arcsec < CFG["tns_arcsec"]
    T = t.iloc[idx]
    F.loc[ok, "tns_name"] = (T.name_prefix.fillna("").values + " " + T.name.astype(str).values)[ok]
    F.loc[ok, "tns_type"] = T.type.fillna("").values[ok]
    F.loc[ok, "tns_date"] = T.discoverydate.fillna("").astype(str).str[:10].values[ok]
    print(f"  {ok.sum()} of {len(F)} candidates match a TNS object within {CFG['tns_arcsec']}\"")
    return F


# ---------------------------------------------------------------- 5. quiet-before check
def add_quiet_flag(F, images):
    """quiet=True if a good image >= quiet_days before the first detection was deep enough
    to have seen the source at (peak - 0.5 mag) and it wasn't detected."""
    F = F.copy()
    F["quiet"] = pd.Series([None] * len(F), index=F.index, dtype=object)
    F["n_pre_images"] = 0
    if images is None or images.empty:
        return F
    for s, r in F.iterrows():
        im = images[(images.src == s) & (images.mjd < r.first_mjd - 1)
                    & (images.mjd > r.first_mjd - 365)]
        if "tqc" in im:
            im = im[im.tqc.astype(str).str.lower().isin(["green", "yellow", "0", "1", "none", "nan"])]
        deep = im[im.t_lmag > r.peak_mag + 0.5]
        F.loc[s, "n_pre_images"] = len(deep)
        F.at[s, "quiet"] = bool(len(deep) > 0)
    # no deep pre-image -> can't confirm it's new; mild penalty only
    F["score"] *= np.where(F.quiet.eq(False), 0.5, 1.0)
    return F.sort_values("score", ascending=False)


# ---------------------------------------------------------------- 6. inspection page
def write_html(F, det, path="nova_candidates.html", n=300):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import base64, io
    colors = dict(u="tab:purple", g="tab:green", q="tab:orange", r="tab:red", i="tab:brown", z="k")
    rows = []
    for rank, (src, s) in enumerate(F.head(n).iterrows(), 1):
        lc = det[det.src == src].sort_values("mjd")
        fig, ax = plt.subplots(figsize=(4.2, 2.6))
        for flt, d in lc.groupby("filt"):
            ax.errorbar(d.mjd, d.mag, d.magerr, fmt="o", ms=3, label=flt, color=colors.get(flt))
        ax.invert_yaxis(); ax.set_xlabel("MJD"); ax.set_ylabel("mag"); ax.legend(fontsize=6)
        buf = io.BytesIO(); fig.tight_layout(); fig.savefig(buf, format="png", dpi=90); plt.close(fig)
        img = base64.b64encode(buf.getvalue()).decode()
        cut = (f"https://www.legacysurvey.org/viewer/cutout.jpg?ra={s.ra:.6f}&dec={s.dec:.6f}"
               f"&layer=ls-dr10&pixscale=0.4&size=200")
        view = (f"https://www.legacysurvey.org/viewer?ra={s.ra:.6f}&dec={s.dec:.6f}"
                f"&zoom=15&mark={s.ra:.6f},{s.dec:.6f}")
        quiet = {True: "yes", False: "<b style='color:red'>no</b>"}.get(s.quiet, "?")
        tn = s.get("tns_name", "")
        tns = (f"<b style='color:#b06000'>TNS: {tn} {s.get('tns_type', '')} ({s.get('tns_date', '')})</b>"
               if isinstance(tn, str) and tn.strip() else "<b style='color:green'>not in TNS</b>")
        rows.append(f"""<tr><td>{rank}<br>src {src}</td>
          <td><b>{s.host}</b> ({s.host_D:.1f} Mpc)<br>RA {s.ra:.5f} Dec {s.dec:.5f}<br>
          peak {s.peak_mag:.2f} ({s.peak_filt}) &rarr; M = {s.M_peak:.1f}<br>
          t<sub>2</sub> ~ {s.t2_days:.0f} d &middot; fade {s.decline:.1f} mag<br>
          {s.n_det} det / {s.n_nights} nights &middot; {s.telescopes}<br>
          quiet before: {quiet} &middot; r/R<sub>ell</sub> = {s.r_ell:.2f}<br>
          {tns}<br>
          Gaia sep {s.sep_gaia:.1f}" &middot; score {s.score:.2f}<br>
          mean RB {s.rb_mean:.2f} &middot; detected in {s.get('det_frac', float('nan')):.0%} of {s.get('det_frac_filt', '')} exposures<br>
          {"<b style='color:red'>in a group of " + str(int(s.n_cand_near)) + " nearby candidates (artefact?)</b><br>" if s.get('n_cand_near', 0) >= CFG['cluster_n'] else ""}
          <a href="{view}" target=_blank>Legacy viewer</a></td>
          <td><img src="{cut}" width=200 height=200></td>
          <td><img src="data:image/png;base64,{img}"></td>
          <td><select><option>?</option><option>nova</option><option>SN</option>
          <option>variable</option><option>bogus</option></select></td></tr>""")
    html = ("<html><head><meta charset=utf-8><title>Nova candidates</title><style>"
            "body{font-family:sans-serif;font-size:13px}td{vertical-align:top;padding:6px;"
            "border-bottom:1px solid #ccc}</style></head><body>"
            f"<h2>Nova candidates: {len(F)} passing cuts (top {min(n, len(F))} shown)</h2>"
            f"<table>{''.join(rows)}</table></body></html>")
    with open(path, "w") as fh:
        fh.write(html)
    print("wrote", path)


# ---------------------------------------------------------------- main
def run(regalade=None, project=PROJECT, dry_run=False, with_images=True,
        tns="../BG_Trans_Query_ID/tns_public_objects.csv"):
    gal = load_galaxies(regalade)
    print(f"{len(gal)} galaxies with {CFG['dmin_mpc']} < D < {CFG['dmax_mpc']} Mpc, Dec < {CFG['dec_max']}")

    import os
    raw = f"nova_hunt_detections_raw_D{CFG['dmax_mpc']:g}_floor{CFG['rb_floor']:g}_v3.parquet"
    gal_file = raw.replace(".parquet", "_galaxies.csv")   # galaxy numbering used by that query
    if os.path.exists(raw) and not dry_run:
        print(f"Using cached detections {raw} (delete it to re-query BigQuery)")
        det = pd.read_parquet(raw)
        # galaxy ids (gid) in the cached detections refer to the list used at query time
        for f in [gal_file, "nova_hunt_galaxies.csv"]:
            if os.path.exists(f):
                gal = pd.read_csv(f)
                gal["name"] = gal["name"].astype(str)
                break
        gal.to_csv(gal_file, index=False)
    else:
        print("Detections query:")
        det = run_query(detections_query(gal), project, dry_run=dry_run)
        if det is None:
            return None
        det.to_parquet(raw)
        gal.to_csv(gal_file, index=False)
    gal.to_csv("nova_hunt_galaxies.csv", index=False)
    print(f"{len(det)} detections with CLASS_REAL >= {CFG['rb_floor']}")

    det = associate(det, gal)
    det = group_sources(det)
    F = det.groupby("src")[[c for c in det.columns if c != "src"]].apply(nova_features).dropna(how="all")
    print(f"{len(F)} distinct sources")
    F = apply_cuts(F)
    print(f"{len(F)} pass the hard nova cuts")

    if with_images and len(F):
        print("Image query (exposure counts + 'quiet before' check):")
        img_cache = raw.replace("detections_raw", "images")
        images = fetch_images(F, project, img_cache)
        F = add_detection_fraction(F, det, images)
        frac_ok = ~(F.det_frac < CFG["det_frac_min"])      # NaN (no image info) is kept
        print(f"  det_frac >= {CFG['det_frac_min']:<6}{frac_ok.sum():>8d}")
        F = F[frac_ok]
        F = add_quiet_flag(F, images)
    F = flag_clusters(F)

    F = flag_tns(F, tns)
    F.to_csv("nova_candidates.csv")
    det.to_csv("nova_hunt_detections.csv", index=False)
    write_html(F, det)
    return F, det


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regalade", help="local REGALADE file (fits/csv/parquet); default: VizieR")
    ap.add_argument("--project", default=PROJECT)
    ap.add_argument("--dmax", type=float, default=CFG["dmax_mpc"])
    ap.add_argument("--rb", type=float, default=CFG["rb_min"], help="min MEAN CLASS_REAL per source")
    ap.add_argument("--rb-floor", type=float, default=CFG["rb_floor"], help="per-detection floor in SQL")
    ap.add_argument("--det-frac", type=float, default=CFG["det_frac_min"],
                    help="min fraction of exposures detected on detected nights (0.5; 0.25 looser)")
    ap.add_argument("--dry-run", action="store_true", help="only print the BigQuery cost estimate")
    ap.add_argument("--tns", default="../BG_Trans_Query_ID/tns_public_objects.csv",
                    help="local TNS public-objects CSV for flagging known transients")
    ap.add_argument("--min-nights", type=int, default=CFG["min_nights"],
                    help="set to 1 to include single-night sources (many more, mostly junk)")
    ap.add_argument("--no-images", action="store_true", help="skip the quiet-before image query")
    a = ap.parse_args()
    CFG["dmax_mpc"], CFG["rb_min"], CFG["min_nights"] = a.dmax, a.rb, a.min_nights
    CFG["rb_floor"], CFG["det_frac_min"] = a.rb_floor, a.det_frac
    run(a.regalade, a.project, a.dry_run, not a.no_images, a.tns)
