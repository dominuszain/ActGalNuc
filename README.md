# INFERNUS [www.doi.org/10.13140/RG.2.2.18744.71680]

# Classification of Emission-Line Galaxies with SDSS and WISE

A fully reproducible astrophysical data analysis pipeline that classifies narrow emission-line galaxies using optical (BPT, WHAN) and mid-infrared (WISE) diagnostics, and derives physical properties including electron density and supermassive black hole mass via the M-sigma relation.

**Author:** Zain Ul Abideen — Erasmus Mundus MASS Programme

## Scientific Summary

The project consists of three linked tasks:

1. **Optical Classification (BPT & WHAN)** — Query SDSS DR18 for ~10,000 narrow emission-line galaxies (z < 0.35). Classify each using the Baldwin–Phillips–Terlevich (BPT) diagnostic diagram with both Kewley+01 (3-class) and Kauffmann+03 (2-class) boundaries, and the WHAN diagram (Cid Fernandes+11, 5-class scheme).

2. **WISE Mid-IR Cross-Matching** — Spatially cross-match the BPT AGN subsample with WISE all-sky photometry. Apply the Mateos+12 AGN wedge in WISE colour-colour space to identify mid-IR AGN and compare with optical classification.

3. **Electron Density & Black Hole Mass** — Compute electron density from the [S II] 6717/6731 doublet using PyNeb. Map [O III] velocity dispersion across the BPT plane, and estimate supermassive black hole masses via the M-sigma relation (Tremaine+02, Gebhardt+00, Wu 2009).

## Dependencies

Python 3 with:

```
numpy>=1.24
scipy>=1.10
matplotlib>=3.7
astropy>=5.2
pyneb>=1.1
requests>=2.28
```

Install with:

```bash
pip install -r requirements.txt
```

To compile the report, you also need a LaTeX distribution with `pdflatex` and `bibtex`.

## Quick Start

```bash
# Full pipeline: download data, run all tasks, compile the report
make all

# After any data change, update the report statistics:
python3 scripts/07_generate_numbers.py
make report

# Clean all generated files
make clean
```

## Pipeline

The pipeline is a linear sequence of 7 scripts, orchestrated by the `Makefile`:

```
SDSS DR18 SkyServer (internet)
      │
      ▼
01_fetch_sdss.py ────────► data/sdss_sample.fits
      │
      ▼
02_bpt_whan.py ──────────► data/sdss_sample_classified.fits
                           output/bpt_diagram.pdf
                           output/whan_diagram.pdf
      │
      ├──► 03_fetch_wise.py ──► data/wise_matched.fits
      │         │
      │         └──► 04_wise_diagram.py ──► output/wise_diagram.pdf
      │
      ├──► 05_filter_sii.py ──► data/sdss_sample_sii.fits
      │         │
      │         └──► 06_density_mass.py ──► 8 output PDFs
      │
      └──► 07_generate_numbers.py ──► report/numbers.tex
                                          │
                                          ▼
                                    make report ──► report/report.pdf
```

### Script Summary

| Script | Input | Output | Requires Internet |
|--------|-------|--------|-------------------|
| `01_fetch_sdss.py` | — | `data/sdss_sample.fits` | Yes |
| `02_bpt_whan.py` | `data/sdss_sample.fits` | `data/sdss_sample_classified.fits` + 2 PDFs | No |
| `03_fetch_wise.py` | `data/sdss_sample_classified.fits` | `data/wise_matched.fits` | Yes |
| `04_wise_diagram.py` | `data/wise_matched.fits` | `output/wise_diagram.pdf` | No |
| `05_filter_sii.py` | `data/sdss_sample_classified.fits` | `data/sdss_sample_sii.fits` | No |
| `06_density_mass.py` | `data/sdss_sample_sii.fits` | 8 output PDFs | No |
| `07_generate_numbers.py` | All 4 FITS files | `report/numbers.tex` | No |

## Project Structure

```
.
├── agn_project/           # Python library (7 modules)
│   ├── __init__.py        # Public API
│   ├── config.py          # Constants, classification curves, M-sigma params
│   ├── sdss.py            # SDSS DR18 SkyServer paginated SQL query engine
│   ├── classification.py  # BPT, WHAN, WISE classification functions
│   ├── plots.py           # All matplotlib figure generation
│   ├── wise.py            # WISE cross-matching and colour computations
│   ├── physics.py         # Electron density (PyNeb) and MBH (M-sigma)
│   └── utils.py           # Cosmology helper (FlatLambdaCDM)
├── scripts/               # Executable pipeline scripts (7 scripts)
├── report/                # LaTeX source and bibliography
│   ├── report.tex         # Full report source
│   ├── references.bib     # 11 BibTeX entries
│   └── numbers.tex        # 74 auto-generated LaTeX \newcommand stubs
├── Makefile               # Build automation
├── requirements.txt       # Python dependencies
└── .gitignore
```

## Package Design

The `agn_project/` library is a reusable Python package. Its `__init__.py` exposes all public functions so scripts import cleanly:

```python
from agn_project import query_sdss, classify_bpt_kewley, electron_density, plot_bpt, ...
```

All constants and analytic curves live in a single `config.py` — the single source of truth. Every other module imports from it, so changing a parameter once updates the entire pipeline.

## Figures Generated

| Figure | Description |
|--------|-------------|
| `bpt_diagram.pdf` | BPT diagram with Kewley+01 (solid) and Kauffmann+03 (dashed) curves |
| `whan_diagram.pdf` | WHAN diagram colored by BPT class |
| `wise_diagram.pdf` | WISE colour-colour diagram with Mateos+12 AGN wedge |
| `bpt_density_map.pdf` | Electron density 2D heatmap on the BPT plane |
| `bpt_sigma_map.pdf` | [O III] velocity dispersion 2D map on the BPT plane |
| `mbh_histogram.pdf` | SMBH mass histogram (Tremaine+02 calibration) |
| `mbh_by_class_histogram.pdf` | SMBH mass distributions split by BPT class |
| `bpt_sigma_mbh_*.pdf` | Sigma maps binned by MBH (4 mass bins: 5-6, 6-7, 7-8, 8-9 dex) |

## Key Design Decisions

- **Cursor-based pagination**: DR18 SkyServer does not support `OFFSET/FETCH`, so pagination uses `WHERE specObjID > last_id ORDER BY specObjID`.
- **SSH/retry logic**: Up to 10 retries with exponential backoff for SkyServer timeouts and SSL errors.
- **Agg backend**: Matplotlib is configured with the `Agg` backend for headless PDF generation.
- **Auto-synced report statistics**: `07_generate_numbers.py` reads all FITS data files and writes 74 `\newcommand` definitions to `report/numbers.tex`, which the LaTeX report includes via `\input{numbers.tex}`. Regenerate this file whenever any upstream data changes, then recompile.

## Caveats

- FITS data files are **not committed**. They are binary, can be large, and are fully reproducible by running `make data` (requires internet access to SDSS SkyServer).
- Output PDFs and LaTeX build artifacts are also excluded.
- The `report/numbers.tex` file in the repository contains placeholder stubs; it must be regenerated before compiling the report.
- WISE cross-matching queries the SDSS SkyServer's `wise_allsky` table rather than IRSA — kept to a single query endpoint for simplicity.

## References

- Baldwin, Phillips & Terlevich (1981)
- Kewley et al. (2001)
- Kauffmann et al. (2003)
- Cid Fernandes et al. (2011)
- Mateos et al. (2012)
- Wright et al. (2010) — WISE
- Tremaine et al. (2002) — M-sigma
- Gebhardt et al. (2000) — M-sigma (alternate)
- Wu (2009) — M-sigma (alternate)
