#!/usr/bin/env python3
"""Task 2: WISE colour-colour diagram + AGN classification stats."""

import sys
sys.path.insert(0, "/home/zain/Desktop/test")

import numpy as np
from astropy.table import Table
from agn_project.wise import wise_colors
from agn_project.classification import classify_wise_agn
from agn_project.plots import plot_wise

DATA_DIR = "/home/zain/Desktop/test/data"
OUT_DIR = "/home/zain/Desktop/test/output"
IN_MATCHED = f"{DATA_DIR}/wise_matched.fits"
IN_SDSS = f"{DATA_DIR}/sdss_sample_classified.fits"


def main():
    matched = Table.read(IN_MATCHED)
    sdss_all = Table.read(IN_SDSS)
    print(f"WISE-matched sample: {len(matched)}")
    print(f"Full SDSS sample:    {len(sdss_all)}")

    # Check required WISE columns
    for col in ["w_w1mpro", "w_w2mpro", "w_w3mpro"]:
        if col not in matched.colnames:
            print(f"ERROR: missing WISE column '{col}' in matched table.")
            return

    # Remove objects with bad photometry (null / upper limits)
    good = (
        (matched["w_w1mpro"] > 0)
        & (matched["w_w2mpro"] > 0)
        & (matched["w_w3mpro"] > 0)
    )
    matched = matched[good]
    print(f"After quality cuts:   {len(matched)}")

    if len(matched) == 0:
        print("No usable WISE data; aborting.")
        return

    x, y = wise_colors(matched)
    wise_agn = classify_wise_agn(x, y)

    print(f"\n{'='*50}")
    print(f"  WISE AGN Classification (Mateos+12 wedge)")
    print(f"{'='*50}")
    print(f"  WISE AGN:           {wise_agn.sum():5d}  ({100*wise_agn.sum()/len(matched):.1f}%)")
    print(f"  WISE non-AGN:       {(~wise_agn).sum():5d}  ({100*(~wise_agn).sum()/len(matched):.1f}%)")

    plot_wise(x, y, wise_agn, title="WISE Colour-Colour – Mateos+12 AGN Wedge",
              fname=f"{OUT_DIR}/wise_diagram.pdf")

    # Cross-comparison with SDSS BPT
    bpt_agn = matched["bpt_kewley"] == "AGN/Seyfert"
    both_agn = bpt_agn & wise_agn
    print(f"\n  SDSS BPT AGN:        {bpt_agn.sum():5d}  ({100*bpt_agn.sum()/len(matched):.1f}%)")
    print(f"  AGN by both surveys: {both_agn.sum():5d}  ({100*both_agn.sum()/len(matched):.1f}%)")
    print(f"  BPT-only AGN:        {(bpt_agn & ~wise_agn).sum():5d}")
    print(f"  WISE-only AGN:       {((~bpt_agn) & wise_agn).sum():5d}")

    # Overall stats percentages
    pct_matched = 100 * len(matched) / len(sdss_all)
    pct_wagn = 100 * wise_agn.sum() / len(matched)
    pct_both = 100 * both_agn.sum() / len(matched)
    print(f"\n  % SDSS objects with WISE match: {pct_matched:.1f}%")
    print(f"  % matched objects that are WISE AGN: {pct_wagn:.1f}%")
    print(f"  % matched objects AGN in both SDSS+WISE: {pct_both:.1f}%")


if __name__ == "__main__":
    main()
