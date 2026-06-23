#!/usr/bin/env python3
"""Task 3: Filter sample to objects with good [SII] and sigma measurements.

Since DR18 galSpecLine contains all needed columns (fluxes, EQW, sigma),
this step filters the classified sample to those with S/N >= 3 on
[SII] 6717, 6731 and valid sigma_forbidden.
"""

import os
import sys
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT_DIR)

import numpy as np
from astropy.table import Table

from agn_project.config import SN_CUT, VELDISP_COL, VELDISP_COL_ALT

DATA_DIR = os.path.join(ROOT_DIR, "data")
IN_FITS = os.path.join(DATA_DIR, "sdss_sample_classified.fits")
OUT_FITS = os.path.join(DATA_DIR, "sdss_sample_sii.fits")


def main():
    tbl = Table.read(IN_FITS)
    print(f"Loaded {len(tbl)} classified galaxies.")

    # Apply [SII] S/N cuts
    for col in ["sii_6717", "sii_6731"]:
        fcol = f"{col}_flux"
        ecol = f"{col}_flux_err"
        if fcol in tbl.colnames and ecol in tbl.colnames:
            sn = tbl[fcol] / tbl[ecol]
            good = (tbl[fcol] > 0) & (tbl[ecol] > 0) & (sn >= SN_CUT)
            tbl = tbl[good]
            print(f"  After {col} S/N >= {SN_CUT}: {len(tbl)}")
        else:
            print(f"  WARNING: column {fcol} not found, keeping all")

    # Require valid sigma
    sigma_col = VELDISP_COL if VELDISP_COL in tbl.colnames else None
    if sigma_col is None and VELDISP_COL_ALT in tbl.colnames:
        sigma_col = VELDISP_COL_ALT

    if sigma_col:
        good = np.isfinite(tbl[sigma_col]) & (tbl[sigma_col] > 0)
        tbl = tbl[good]
        print(f"  After sigma > 0: {len(tbl)}")
    else:
        print("  WARNING: no sigma column found")

    tbl.write(OUT_FITS, overwrite=True)
    print(f"Saved {len(tbl)} objects to {OUT_FITS}")


if __name__ == "__main__":
    main()
