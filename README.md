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
| `--min-nights` | 2 | Minimum number of nights with a detection |
| `--regalade FILE` | VizieR | Use a local REGALADE file (FITS/CSV/parquet) instead of downloading it |
| `--tns FILE` | `../BG_Trans_Query_ID/tns_public_objects.csv` | Local TNS public-objects CSV used to flag known transients |
| `--project` | `blackgem-databases-access` | Google Cloud project the BigQuery jobs are billed to |
| `--no-images` | off | Skip the image query (no exposure fraction or "quiet before" check) |

**Outputs**
- `nova_candidates.html`: ranked inspection page with Legacy Survey cutouts, light curves, host galaxy, TNS match and flags, plus a dropdown to mark each candidate (nova, SN, variable or bogus).
- `nova_candidates.csv`: the same candidates with all computed features.

**Caching:** the galaxy list, detections and image query results are saved locally. Rerunning with different cuts costs nothing; changing `--dmax` or `--rb-floor` triggers a new query. A detections query scans about 100 GB (about $0.50 on-demand); image queries scan under 0.1 GB.

## Method

1. **Galaxies:** [REGALADE](https://doi.org/10.1051/0004-6361/202556896) (Tranin et al. 2026, A&A 706, A284; VizieR `J/A+A/706/A284`) with 0.1 < D < `dmax` Mpc and Dec < +30°. Galaxies are ranked by expected nova yield (stellar mass × fraction of the nova luminosity function above the depth limit). The download tries several VizieR mirrors in 15° RA strips and saves the result locally.
2. **Detections:** one BigQuery join of `historic_detections` and `images` against all galaxies, keeping detections within 2× each galaxy's ellipse with CLASS_REAL ≥ `rb-floor`.
3. **Host matching:** normalised elliptical separation (R1, R2, PA east of north), the same convention as [ellipse_xmatch](https://github.com/htranin/ellipse_xmatch).
4. **Sources:** detections within 2″ are linked into one source, and light-curve features are computed (peak magnitude, absolute magnitude, fade, t₂, number of nights, duration).
5. **Real/bogus:** at least 2 detections, mean CLASS_REAL > `rb`, and detected in at least `det-frac` of the exposures covering the position on the nights it was seen. The fraction is **counted per filter** (best filter with ≥ 2 exposures), because novae are detected mostly in q and less in i/u.
6. **Nova cuts:** −10.5 < M_peak < −4.5; more than 2″ from the galaxy nucleus; detected on ≥ 2 nights; duration < 200 days; fades ≥ 0.5 mag or all detections within 30 days; no Gaia star within 1.5″; no reference-image source within 1.5″ brighter than peak + 1 mag.
7. **Flags and ranking:**
   - "Quiet before": a deep earlier image (from T_LMAG) with no detection.
   - Artefact groups: candidates with ≥ 3 other candidates within 90″, typical of bright-star bleeds and ghosts, are scored down.
   - TNS: match to a known object within 3″.

All thresholds are in the `CFG` dictionary at the top of `nova_hunt.py`.

## Validation (5 Mpc, November 2026)

At 5 Mpc there are 175 galaxies, 473,613 detections and 140,326 sources, of which **105 candidates** pass with `--det-frac 0.25`. The candidates include **13 of the 14** known TNS transients in the sample, among them all three classified novae (AT 2025nok, AT 2024aeql, AT 2025cfa). The top 7 ranked candidates are all known TNS transients in Centaurus A.

## Notes and known limitations

- A bright-star distance cut was tested and disabled. BlackGEM's "bright Gaia" stars are G < ~14, and real novae in Cen A lie 17–90″ from G = 11–13 stars.
- Galaxies without an ellipse fit get a 30″ radius; very large galaxies are capped at 30′.
- Novae beyond about 15–20 Mpc are close to the single-epoch detection limit.
- Data products are not committed (see `.gitignore`); they contain unpublished survey data.

## Acknowledgements

The method follows the approach and real/bogus criteria of H. Tranin (ellipse_xmatch, REGALADE). It uses MeerLICHT/BlackGEM data, REGALADE, Gaia, the Legacy Surveys viewer and the Transient Name Server.
