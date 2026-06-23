"""Physics computations: electron density, velocity dispersion, M-sigma."""

import numpy as np
import pyneb as pn

from .config import TE_NLR, MSIGMA_PARAMS


# ---------------------------------------------------------------------------
# Electron density from [S II] 6717 / 6731
# ---------------------------------------------------------------------------

def electron_density(sii_ratio, te=TE_NLR):
    """Compute electron density from [S II] 6717/6731 flux ratio.

    Uses PyNeb's 5-level atom model for [S II].

    Parameters
    ----------
    sii_ratio : array
        F(6717) / F(6731) flux ratio.
    te : float
        Electron temperature in K (default 10 000 K for NLR).

    Returns
    -------
    ne : array
        Electron density in cm^-3.  Values clipped to [1, 1e6].
        NaN where ratio out of theoretical bounds.
    """
    # PyNeb's temden.getTemDen expects emissivity ratio, which for [SII]
    # (6717/6731) is the flux ratio.
    sii = pn.Atom("S", 2)
    ne = np.full_like(sii_ratio, np.nan, dtype=float)

    for i, r in enumerate(sii_ratio):
        if np.isnan(r) or r <= 0:
            continue
        try:
            dens = sii.getTemDen(r, tem=te, wave1=6717, wave2=6731)
            ne[i] = dens
        except Exception:
            # ratio out of theoretical range -> leave as NaN
            pass

    # Clip to physically meaningful range for NLR
    ne = np.clip(ne, 1, 1e6)
    return ne


# ---------------------------------------------------------------------------
# Supermassive black hole mass from M-sigma relation
# ---------------------------------------------------------------------------

def mbh_from_sigma(sigma_vel, relation="tremaine2002"):
    """Compute SMBH mass from stellar / gas velocity dispersion.

    Parameters
    ----------
    sigma_vel : array
        Velocity dispersion in km/s.
    relation : str
        One of 'tremaine2002', 'gebhardt2000', 'wu2009'.

    Returns
    -------
    mbh : array
        log10(M_BH / M_sun).
    """
    params = MSIGMA_PARAMS[relation]
    alpha = params["alpha"]
    beta = params["beta"]
    sigma0 = params["sigma0"]

    sigma_ratio = np.clip(sigma_vel, 10, 500) / sigma0
    log_mbh = alpha + beta * np.log10(sigma_ratio)
    return log_mbh


def mbh_all_relations(sigma_vel):
    """Return a dict of log(MBH) estimates for all three relations."""
    return {key: mbh_from_sigma(sigma_vel, key) for key in MSIGMA_PARAMS}
