"""BPT, WHAN, and WISE classification functions."""

import numpy as np
from .config import (
    _kewley_line,
    _kauffmann_line,
    WHAN_NII_HA_THRESHOLD,
    WHAN_EWHA_SF_THRESHOLD,
    WHAN_EWHA_WAGN_THRESHOLD,
    WHAN_EWHA_RETIRED_THRESHOLD,
    WHAN_EW_PASSIVE_THRESHOLD,
    is_wise_agn,
)


def bpt_flux_ratios(table):
    """Compute log([OIII]/Hb) and log([NII]/Ha) from flux columns."""
    o3hb = table["oiii_5007_flux"] / table["h_beta_flux"]
    n2ha = table["nii_6584_flux"] / table["h_alpha_flux"]
    log_o3hb = np.log10(o3hb)
    log_n2ha = np.log10(n2ha)
    return log_o3hb, log_n2ha


def classify_bpt_kewley(log_n2ha, log_o3hb):
    """Classify using Kewley+2001 theoretical line (3 classes).

    Star-forming: below Kauffmann line
    Composite:    between Kauffmann and Kewley
    AGN/Seyfert:  above Kewley line, or log([NII]/Ha) >= 0.47
    """
    labels = np.full(len(log_n2ha), "Unclassified", dtype=object)
    kew = _kewley_line(log_n2ha)
    kau = _kauffmann_line(log_n2ha)

    kau_valid = log_n2ha < 0.05
    sf_mask = kau_valid & np.isfinite(kau) & (log_o3hb < kau)
    labels[sf_mask] = "Star-forming"

    kew_valid = log_n2ha < 0.47
    comp_mask = kew_valid & kau_valid & (log_o3hb >= kau) & (log_o3hb < kew)
    labels[comp_mask] = "Composite"

    agn1 = kew_valid & (log_o3hb >= kew)
    labels[agn1] = "AGN/Seyfert"

    # Objects with log >= 0.05 (beyond Kauffmann) but within Kewley range
    comp2 = (~kau_valid) & kew_valid & (log_o3hb < kew)
    labels[comp2] = "Composite"

    agn2 = ~kew_valid  # log >= 0.47 -> Seyfert/LINER territory
    labels[agn2] = "AGN/Seyfert"

    return labels


def classify_bpt_kauffmann(log_n2ha, log_o3hb):
    """Classify using Kauffmann+2003 empirical line (2-class split).

    Objects with log([NII]/Ha) >= 0.05 are classified as AGN.
    """
    labels = np.full(len(log_n2ha), "Unclassified", dtype=object)
    kau = _kauffmann_line(log_n2ha)

    valid = log_n2ha < 0.05
    sf = valid & np.isfinite(kau) & (log_o3hb < kau)
    agn_in_range = valid & np.isfinite(kau) & (log_o3hb >= kau)
    agn_beyond = ~valid

    labels[sf] = "Star-forming"
    labels[agn_in_range] = "AGN"
    labels[agn_beyond] = "AGN"
    return labels


def whan_ratios(table):
    """Return log([NII]/Ha) and H-alpha equivalent width (Ang).

    DR18 galSpecLine stores EQWs with emission as negative values;
    we take the absolute value for the WHAN diagram.
    """
    n2ha = table["nii_6584_flux"] / table["h_alpha_flux"]
    log_n2ha = np.log10(np.clip(n2ha, 1e-30, None))
    ew_ha = np.abs(table["h_alpha_eqw"])
    return log_n2ha, ew_ha


def classify_whan(log_n2ha, ew_ha):
    """Classify galaxies via the WHAN diagram (Cid Fernandes+11).

    Returns labels: Passive, Retired, wAGN, Strong AGN, Star-forming,
    or Unclassified.
    """
    labels = np.full(len(log_n2ha), "Unclassified", dtype=object)

    passive = ew_ha < WHAN_EW_PASSIVE_THRESHOLD
    retired = (ew_ha >= WHAN_EW_PASSIVE_THRESHOLD) & (ew_ha < WHAN_EWHA_RETIRED_THRESHOLD)
    low_n2 = log_n2ha <= WHAN_NII_HA_THRESHOLD
    high_n2 = log_n2ha > WHAN_NII_HA_THRESHOLD

    sf = low_n2 & (ew_ha >= WHAN_EWHA_SF_THRESHOLD)
    strong_agn = high_n2 & (ew_ha >= WHAN_EWHA_WAGN_THRESHOLD)
    weak_agn = high_n2 & (ew_ha >= WHAN_EWHA_SF_THRESHOLD) & (ew_ha < WHAN_EWHA_WAGN_THRESHOLD)

    labels[passive] = "Passive"
    labels[retired] = "Retired"
    labels[sf] = "Star-forming"
    labels[strong_agn] = "Strong AGN"
    labels[weak_agn] = "wAGN"
    return labels


def classify_wise_agn(x, y):
    """Wrapper for WISE AGN wedge classification."""
    return is_wise_agn(x, y)
