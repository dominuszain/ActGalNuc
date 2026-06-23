#!/usr/bin/env python3
"""Task 1: Produce BPT and WHAN diagrams + classification statistics."""

import sys
sys.path.insert(0, "/home/zain/Desktop/test")

import numpy as np
from astropy.table import Table
from agn_project.classification import (
    bpt_flux_ratios,
    classify_bpt_kewley,
    classify_bpt_kauffmann,
    whan_ratios,
    classify_whan,
)
from agn_project.plots import plot_bpt, plot_whan

DATA_DIR = "/home/zain/Desktop/test/data"
IN_FITS = f"{DATA_DIR}/sdss_sample.fits"


def stats_table(labels, name):
    """Print label distribution."""
    unique, counts = np.unique(labels, return_counts=True)
    total = len(labels)
    print(f"\n{'='*50}")
    print(f"  {name}  (N = {total})")
    print(f"{'='*50}")
    for u, c in sorted(zip(unique, counts), key=lambda x: -x[1]):
        print(f"  {u:20s}  {c:6d}  ({100*c/total:.1f}%)")


def main():
    tbl = Table.read(IN_FITS)
    print(f"Loaded {len(tbl)} galaxies from {IN_FITS}")

    # ---- BPT ----
    log_o3hb, log_n2ha = bpt_flux_ratios(tbl)
    labels_kewley = classify_bpt_kewley(log_n2ha, log_o3hb)
    labels_kauff = classify_bpt_kauffmann(log_n2ha, log_o3hb)

    stats_table(labels_kewley, "BPT – Kewley+01 (3-class)")
    stats_table(labels_kauff, "BPT – Kauffmann+03 (2-class)")

    plot_bpt(log_n2ha, log_o3hb, labels_kewley,
             title="BPT Diagram – Kewley+01 & Kauffmann+03",
             fname="/home/zain/Desktop/test/output/bpt_diagram.pdf")

    # ---- WHAN ----
    log_n2ha_w, ew_ha = whan_ratios(tbl)
    whan_labels = classify_whan(log_n2ha_w, ew_ha)
    stats_table(whan_labels, "WHAN (Cid Fernandes+11)")

    plot_whan(log_n2ha_w, ew_ha, labels_kewley,  # BPT labels per task 1b
              title="WHAN Diagram — coloured by BPT class",
              fname="/home/zain/Desktop/test/output/whan_diagram.pdf")

    # Cross-comparison: BPT(Kewley) vs WHAN AGN
    bpt_agn = labels_kewley == "AGN/Seyfert"
    whan_agn = (whan_labels == "Strong AGN") | (whan_labels == "wAGN")
    both_agn = bpt_agn & whan_agn
    bpt_only = bpt_agn & ~whan_agn
    whan_only = ~bpt_agn & whan_agn

    print(f"\n{'='*50}")
    print(f"  Cross-comparison: BPT(Kewley) vs WHAN AGN")
    print(f"{'='*50}")
    print(f"  AGN in BPT:             {bpt_agn.sum():5d}  ({100*bpt_agn.sum()/len(tbl):.1f}%)")
    print(f"  AGN in WHAN:            {whan_agn.sum():5d}  ({100*whan_agn.sum()/len(tbl):.1f}%)")
    print(f"  AGN by both:            {both_agn.sum():5d}  ({100*both_agn.sum()/len(tbl):.1f}%)")
    print(f"  BPT-only AGN:           {bpt_only.sum():5d}  ({100*bpt_only.sum()/len(tbl):.1f}%)")
    print(f"  WHAN-only AGN:          {whan_only.sum():5d}  ({100*whan_only.sum()/len(tbl):.1f}%)")

    print(f"\n  WHAN-only AGN breakdown:")
    whan_only_labs = whan_labels[whan_only]
    for lab in ["Strong AGN", "wAGN"]:
        n = np.sum(whan_only_labs == lab)
        print(f"    {lab:15s}: {n:5d}  ({100*n/whan_only.sum():.1f}%)")

    print(f"\n  BPT-only AGN breakdown (by WHAN class):")
    bpt_only_labs = whan_labels[bpt_only]
    for lab in np.unique(bpt_only_labs):
        n = np.sum(bpt_only_labs == lab)
        print(f"    {lab:15s}: {n:5d}  ({100*n/bpt_only.sum():.1f}%)")

    # Merge WHAN labels into data for later use — ensure consistent string types
    tbl["bpt_kewley"] = np.array(labels_kewley, dtype="U20")
    tbl["bpt_kauff"] = np.array(labels_kauff, dtype="U20")
    tbl["whan_label"] = np.array(whan_labels, dtype="U20")
    tbl["log_o3hb"] = log_o3hb
    tbl["log_n2ha"] = log_n2ha
    tbl.write(IN_FITS.replace(".fits", "_classified.fits"), overwrite=True)
    print(f"\nSaved classified table to {IN_FITS.replace('.fits', '_classified.fits')}")


if __name__ == "__main__":
    main()
