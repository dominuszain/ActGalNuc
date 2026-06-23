#!/usr/bin/env python3
"""Task 3: Electron density, velocity-dispersion maps, and black-hole mass.

Produces:
  1. 2D BPT map coloured by electron density (from [SII] ratio)
  2. 2D BPT map coloured by [OIII] velocity dispersion
  3. Supermassive black hole mass histograms (M-sigma)
  4. Binned (by MBH) [OIII] velocity-dispersion maps on the BPT plane
"""

import sys
sys.path.insert(0, "/home/zain/Desktop/test")

import numpy as np
from astropy.table import Table
from agn_project.physics import electron_density, mbh_from_sigma
from agn_project.plots import (
    plot_bpt_2d_map,
    plot_histogram,
    plot_histogram_by_class,
)
from agn_project.classification import bpt_flux_ratios

DATA_DIR = "/home/zain/Desktop/test/data"
OUT_DIR = "/home/zain/Desktop/test/output"
IN_FITS = f"{DATA_DIR}/sdss_sample_sii.fits"


def main():
    tbl = Table.read(IN_FITS)
    print(f"Loaded {len(tbl)} objects from {IN_FITS}")

    # Ensure we have BPT ratios
    if "log_n2ha" not in tbl.colnames or "log_o3hb" not in tbl.colnames:
        log_o3hb, log_n2ha = bpt_flux_ratios(tbl)
    else:
        log_n2ha = tbl["log_n2ha"]
        log_o3hb = tbl["log_o3hb"]

    # ---- Electron density ----
    print("\nComputing electron density from [SII] 6717/6731 ...")
    sii_ratio = tbl["sii_6717_flux"] / tbl["sii_6731_flux"]
    ne = electron_density(sii_ratio)
    mask_ne = np.isfinite(ne)
    print(f"  Objects with valid ne: {mask_ne.sum()} / {len(ne)}")

    if mask_ne.sum() > 20:
        plot_bpt_2d_map(
            log_n2ha[mask_ne], log_o3hb[mask_ne], ne[mask_ne],
            cbar_label=r"$n_e\ [{\rm cm}^{-3}]$",
            title="Electron Density on the BPT Diagram",
            fname=f"{OUT_DIR}/bpt_density_map.pdf",
            log_norm=True,
        )

    # ---- [OIII] velocity dispersion ----
    from agn_project.config import VELDISP_COL, VELDISP_COL_ALT
    sigma_col = None
    if VELDISP_COL in tbl.colnames:
        sigma_col = VELDISP_COL
    elif VELDISP_COL_ALT in tbl.colnames:
        sigma_col = VELDISP_COL_ALT
    if sigma_col is None:
        print("WARNING: no velocity dispersion column found; skipping kinematics.")
        return
    sigma_oiii = tbl[sigma_col]

    mask_sig = np.isfinite(sigma_oiii) & (sigma_oiii > 0) & (sigma_oiii < 500)
    print(f"\nObjects with valid sigma: {mask_sig.sum()} / {len(sigma_oiii)}")

    if mask_sig.sum() > 20:
        plot_bpt_2d_map(
            log_n2ha[mask_sig], log_o3hb[mask_sig], sigma_oiii[mask_sig],
            cbar_label=r"$\sigma_{[{\rm O\,III}]}\ [{\rm km\,s}^{-1}]$",
            title=r"[O III] Velocity Dispersion on the BPT Diagram",
            fname=f"{OUT_DIR}/bpt_sigma_map.pdf",
            log_norm=False,
        )

    # ---- Supermassive black hole mass (M-sigma) ----
    print("\nComputing SMBH mass (Tremaine+02 M-sigma) ...")
    log_mbh = mbh_from_sigma(sigma_oiii, relation="tremaine2002")
    mask_mbh = np.isfinite(log_mbh) & (log_mbh > 4) & (log_mbh < 10)
    print(f"  Objects with valid M_BH: {mask_mbh.sum()}")

    if mask_mbh.sum() > 10:
        plot_histogram(
            log_mbh[mask_mbh], bins=30,
            xlabel=r"$\log_{10}\,(M_{\rm BH} / M_\odot)$",
            title="Supermassive Black Hole Mass (Tremaine+02 M-sigma)",
            fname=f"{OUT_DIR}/mbh_histogram.pdf",
        )

        # Histogram by BPT class
        if "bpt_kewley" in tbl.colnames:
            classes = ["Star-forming", "Composite", "AGN/Seyfert"]
            labels_str = tbl["bpt_kewley"].astype(str)  # FITS stores as bytes
            mask_all = mask_mbh & np.isin(labels_str, classes)
            plot_histogram_by_class(
                log_mbh[mask_all], tbl["bpt_kewley"][mask_all], classes,
                bins=25,
                xlabel=r"$\log_{10}\,(M_{\rm BH} / M_\odot)$",
                title="SMBH Mass by BPT Class (Kewley+01)",
                fname=f"{OUT_DIR}/mbh_by_class_histogram.pdf",
            )

        # Compare three M-sigma relations
        from agn_project.physics import mbh_all_relations
        print("\nMBH comparison (Tremaine+02 vs Gebhardt+00 vs Wu 2009):")
        all_mbh = mbh_all_relations(sigma_oiii)
        for rel, vals in all_mbh.items():
            mask = np.isfinite(vals) & (vals > 4) & (vals < 10)
            print(f"  {rel}: median log(MBH) = {np.median(vals[mask]):.2f}, "
                  f"N = {mask.sum()}")

        # ---- Binned analysis (1 dex bins in MBH) ----
        print("\nBinned analysis by M_BH (1 dex bins) ...")
        try:
            bins = np.arange(5, 10, 1)  # 5-6, 6-7, 7-8, 8-9
            bin_mask = mask_mbh & mask_sig
            for i in range(len(bins) - 1):
                lo, hi = bins[i], bins[i + 1]
                in_bin = bin_mask & (log_mbh >= lo) & (log_mbh < hi)
                if in_bin.sum() < 10:
                    print(f"  Bin {lo}-{hi}: too few objects ({in_bin.sum()}), skipping")
                    continue
                print(f"  Bin {lo}-{hi}: {in_bin.sum()} objects")
                plot_bpt_2d_map(
                    log_n2ha[in_bin], log_o3hb[in_bin], sigma_oiii[in_bin],
                    cbar_label=r"$\sigma_{[{\rm O\,III}]}\ [{\rm km\,s}^{-1}]$",
                    title=rf"[O III] $\sigma$ on BPT — $M_{{\rm BH}}$ {lo}-{hi} dex",
                    fname=f"{OUT_DIR}/bpt_sigma_mbh_{lo}_{hi}.pdf",
                    log_norm=False,
                )
        except Exception as e:
            print(f"  Binned analysis failed: {e}")

    print("\nTask 3 complete. All figures saved to output/.")


if __name__ == "__main__":
    main()
