#!/usr/bin/env python3
"""Task 2: Cross-match SDSS AGN sample with WISE All-Sky catalogue."""

import sys
sys.path.insert(0, "/home/zain/Desktop/test")

from astropy.table import Table
from agn_project.wise import crossmatch_wise

DATA_DIR = "/home/zain/Desktop/test/data"
IN_FITS = f"{DATA_DIR}/sdss_sample_classified.fits"
OUT_FITS = f"{DATA_DIR}/wise_matched.fits"


def main():
    tbl = Table.read(IN_FITS)
    print(f"Loaded {len(tbl)} galaxies from {IN_FITS}")

    # Cross-match only AGN subsample to keep load light
    agn_mask = tbl["bpt_kewley"] == "AGN/Seyfert"
    agn_tbl = tbl[agn_mask]
    print(f"Cross-matching {len(agn_tbl)} BPT AGN with WISE (chunks of 200, 12s delay) ...")
    matched = crossmatch_wise(agn_tbl, radius_arcsec=3.0, chunk_size=200, delay=12.0)
    print(f"Matched objects: {len(matched)}")

    if len(matched) > 0:
        matched.write(OUT_FITS, overwrite=True)
        print(f"Saved to {OUT_FITS}")
    else:
        print("WARNING: no WISE matches found; writing empty table.")
        Table().write(OUT_FITS, overwrite=True)


if __name__ == "__main__":
    main()
