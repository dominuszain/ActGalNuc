#!/usr/bin/env python3
"""Task 1a: Download SDSS narrow emission-line galaxies (z < 0.35)."""

import os
import sys
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT_DIR)

import numpy as np
from astropy.table import Table
from agn_project.sdss import query_sdss, sdss_sql_narrow_line_galaxies
from agn_project.config import SDSS_MAX_ROWS

DATA_DIR = os.path.join(ROOT_DIR, "data")
OUT_FITS = os.path.join(DATA_DIR, "sdss_sample.fits")


def main():
    sql = sdss_sql_narrow_line_galaxies()
    print("Querying SDSS DR18 SkyServer (paginated, up to 10000 rows) ...")
    tbl = query_sdss(sql, max_rows=SDSS_MAX_ROWS, chunk_size=400, delay=8.0)
    print(f"Downloaded {len(tbl)} rows.")

    # Remove rows with missing essential data
    good = (
        (tbl["oiii_5007_flux"] > 0)
        & (tbl["h_beta_flux"] > 0)
        & (tbl["nii_6584_flux"] > 0)
        & (tbl["h_alpha_flux"] > 0)
    )
    tbl = tbl[good]
    print(f"After flux > 0 filter: {len(tbl)} rows.")

    # Remove objects with physically impossible fluxes (fiber contamination)
    bad_ids = [309661264211109888, 315351506763147264]
    tbl = tbl[~np.isin(tbl["specObjID"], bad_ids)]
    print(f"After removing {len(bad_ids)} objects with impossible fluxes: {len(tbl)} rows.")

    tbl.write(OUT_FITS, overwrite=True)
    print(f"Saved to {OUT_FITS}")


if __name__ == "__main__":
    main()
