# nova_hunt

Systematic search for **missed novae in nearby galaxies** in the MeerLICHT/BlackGEM transient database.

Every MeerLICHT/BlackGEM transient detection inside the elliptical footprint of a nearby REGALADE galaxy is pulled with a single BigQuery spatial join. The detections are grouped into sources and filtered with real/bogus and nova light-curve cuts. The script then writes a ranked HTML page for visual inspection, with known TNS objects flagged.

## Requirements

- Read access to the BlackGEM transient database on Google BigQuery (`blackgem-transient-db`), and permission to run query jobs in a billing project (default `blackgem-databases-access`; the script falls back to `blackgem-transient-db` and `blackgem`).
- Google Cloud credentials on your machine (`gcloud auth application-default login`).
- Python 3.10+ and the packages in `requirements.txt`:

```bash
pip install -r requirements.txt
```

## Usage

```bash
python nova_hunt.py --dmax 5 --dry-run           # galaxy list + BigQuery cost estimate, nothing pulled
python nova_hunt.py --dmax 5 --det-frac 0.25     # full run
```

| Option | Default | Meaning |
|---|---|---|
| `--dmax` | 20 | Maximum galaxy distance (Mpc) |
| `--rb` | 0.5 | Minimum **mean** CLASS_REAL per source |
| `--rb-floor` | 0.1 | Per-detection CLASS_REAL floor in the SQL query |
| `--det-frac` | 0.5 | Minimum fraction of exposures detected on the nights a source was seen (0.25 = looser) |
| `--min-nights` | 2 | Minimum number of nights for multi-night sources |
| `--no-single-night` | off | Drop single-night sources entirely (see Method, step 6) |
| `--no-name-lookup` | off | Don't look up host-galaxy names online (`galaxy_names.csv` is still used) |
| `--rank-by-host` | off | Multiply the score by the host-galaxy prior (old behaviour) |
| `--allow-empty-strips` | off | Accept REGALADE RA strips that VizieR returns empty (normally a VizieR error; see Notes) |
| `--regalade FILE` | VizieR | Read galaxies from a local REGALADE file instead of VizieR, e.g. the REGALADE v2 FITS (4 GB) linked from [ellipse_xmatch](https://github.com/htranin/ellipse_xmatch). Only the needed columns and rows are read; the subset is saved as `regalade_D<dmax>Mpc.csv`, so later runs need neither the file nor VizieR. Recommended: VizieR has become unreliable |
| `--tns FILE` | `../BG_Trans_Query_ID/tns_public_objects.csv` | Local TNS public-objects CSV used to flag known transients |
| `--project` | `blackgem-databases-access` | Google Cloud project the BigQuery jobs are billed to |
| `--no-images` | off | Skip the image query (no exposure fraction or "quiet before" check) |

**Outputs** (named by distance limit, so runs at different distances do not overwrite each other)
- `nova_candidates_D<dmax>.html`: ranked inspection page with Legacy Survey cutouts, light curves, host galaxy, TNS match and flags, plus a dropdown to mark each candidate (nova, SN, variable or bogus).
- `nova_candidates_D<dmax>.csv`: the same candidates with all computed features.

**Host names:** REGALADE v2 'thin' has no name column. Candidate hosts are named from `galaxy_names.csv` (public names and positions; pre-filled with the 175 galaxies within 5 Mpc from REGALADE on VizieR, and tracked in git), otherwise by a NED lookup by position (SIMBAD as fallback). Lookups are added to the file, so each galaxy is looked up once. PGC-only names are shown with their common name, e.g. `PGC2789 = NGC 253`.

**Caching:** the galaxy list, detections and image query results are saved locally, so rerunning with different cuts costs nothing. A new query is triggered when `--dmax` or `--rb-floor` changes, or when the galaxy list differs from the one the saved detections were made with. Saved image results store each source's position and are re-queried if the source numbering changes. A detections query scans about 100 GB (about $0.50 on-demand); image queries scan under 0.1 GB.

## Method

1. **Galaxies:** [REGALADE](https://doi.org/10.1051/0004-6361/202556896) (Tranin et al. 2026, A&A 706, A284; VizieR `J/A+A/706/A284`) with 0.1 < D < `dmax` Mpc and Dec < +30°. Galaxies are ranked by expected nova yield (stellar mass × fraction of the nova luminosity function above the depth limit). The download tries several VizieR mirrors in 15° RA strips and saves the result locally.
2. **Detections:** one BigQuery join of `historic_detections` and `images` against all galaxies, keeping detections within 2× each galaxy's ellipse with CLASS_REAL ≥ `rb-floor`.
3. **Host matching:** normalised elliptical separation (R1, R2, PA east of north), the same convention as [ellipse_xmatch](https://github.com/htranin/ellipse_xmatch).
4. **Sources:** detections within 2″ are linked into one source, and light-curve features are computed (peak magnitude, absolute magnitude, fade, t₂, number of nights, duration).
5. **Real/bogus:** at least 2 detections, mean CLASS_REAL > `rb`, and detected in at least `det-frac` of the exposures covering the position on the nights it was seen. The fraction is **counted per filter** (best filter with ≥ 2 exposures), because novae are detected mostly in q and less in i/u. Peak brightness and fade are measured from confident detections (CLASS_REAL ≥ `rb`) when there are any, so a stray low-score point cannot set the peak.
6. **Nights:** detected on ≥ 2 nights, **or** on a single night with independent confirmation, which means either:
   - ≥ 4 detections from ≥ 2 telescopes or in ≥ 2 filters, seen in ≥ 50% of that night's exposures; or
   - 2–3 detections that all have CLASS_REAL ≥ 0.7 and span ≥ 30 minutes. Junk pairs are back-to-back exposures minutes apart; AT 2024pid spans 1.25 h.

   Single-night candidates are dropped if ≥ 3 of them share a host galaxy and night (a bad epoch).
7. **Nova cuts:**
   - −10.5 < M_peak < −4.5.
   - More than 2″ from the galaxy nucleus.
   - Duration < 200 days.
   - Fades ≥ 0.5 mag, or all detections within 30 days.
   - **Gaia:** reject if a Gaia source lies within 0.5″ (the transient's own star). Between 0.5″ and 1.5″, reject only if the Gaia source is brighter than peak + 1 mag; faint Gaia sources near novae in nearby galaxies are clusters or blends (AT 2025andm).
   - No reference-image source within 1.5″ brighter than peak + 1 mag.
8. **Flags and ranking:**
   - "Quiet before": a deep earlier image (from T_LMAG) with no detection.
   - Artefact groups: candidates with ≥ 3 other candidates within 90″ are scored down.
   - TNS: match to a known object within 3″.
   - **Score:** an ad hoc product of fade, t₂, number of nights, median CLASS_REAL, the group flag and the quiet flag. It orders the list only. The host-galaxy prior (stellar mass × visible fraction of the nova luminosity function) is shown as a separate column.

All thresholds are in the `CFG` dictionary at the top of `nova_hunt.py`.

## Validation (5 Mpc, November 2026)

At 5 Mpc there are 175 galaxies, 473,613 detections and 140,326 sources, of which **156 candidates** pass with `--det-frac 0.25` (95 seen on several nights, 61 on a single night). The candidates include **17 known TNS transients**, among them **all four classified novae** in the sample: AT 2025nok, AT 2024aeql, AT 2025cfa and AT 2024pid. The BlackGEM nova AT 2025andm in NGC 55 is also recovered.

A recall test against all 69 TNS transients in these galaxies since 2020 found:
- 42 were never detected by BlackGEM/MeerLICHT (coverage or depth).
- Of the 27 that were detected, 16 are recovered.
- The rest are single detections or long-lived variables.

`tests/test_rules.py` builds one synthetic source per rule (21 cases) and checks that each is kept or rejected as intended.

## Refreshing the TNS list: `tns_update.py`

`nova_hunt.py` flags known transients from a local TNS export (`../BG_Trans_Query_ID/tns_public_objects.csv`). Keep it current, because an old export misses recent objects (e.g. AT 2025andm):

```bash
python tns_update.py            # uses the MeerLICHT_Bot1 TNS bot; key from api_key_MeerLICHT_Bot1_TNS.txt or $TNS_API_KEY
```

The previous file is kept as `.bak`. The key file is a credential and is excluded by `.gitignore`.

## Light curves: `nova_lightcurves.py`

For a short list of candidates (`new_nova_candidates.csv`: rank, kind, galaxy, host_D, ra, dec, M_peak, first_mjd, last_mjd), this script pulls:
- every BlackGEM/MeerLICHT detection within 3″ at any CLASS_REAL, plus the depth limit of every image covering the position;
- ATLAS forced photometry for the full ATLAS history (needs a free account at fallingstar-data.com/forcedphot; set `ATLAS_TOKEN`, or `ATLAS_USER` and `ATLAS_PASSWORD`).

It plots both together, grouped by host galaxy, in `nova_lightcurves.html`.

```bash
python nova_lightcurves.py --dry-run     # BigQuery cost estimate
python nova_lightcurves.py               # BG/ML + ATLAS + plots
python nova_lightcurves.py --no-atlas    # BG/ML only
```

## Notes and known limitations

- A bright-star distance cut was tested and disabled. BlackGEM's "bright Gaia" stars are G < ~14, and real novae in Cen A lie 17–90″ from G = 11–13 stars.
- Galaxies without an ellipse fit get a 30″ radius; very large galaxies are capped at 30′.
- Novae beyond about 15–20 Mpc are close to the single-epoch detection limit.
- **VizieR can return an empty result instead of an error.** The galaxy download treats an empty RA strip as a failure: it retries, never saves the strip, and stops if it stays empty. Every new list is also checked to contain all galaxies of any smaller-distance list already downloaded.
- Data products are not committed (see `.gitignore`); they contain unpublished survey data.

## Acknowledgements

The method follows the approach and real/bogus criteria of H. Tranin (ellipse_xmatch, REGALADE). It uses MeerLICHT/BlackGEM data, REGALADE, Gaia, the Legacy Surveys viewer and the Transient Name Server.
