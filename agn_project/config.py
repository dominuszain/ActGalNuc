"""Configuration, constants, and classification curves for the AGN project."""

import numpy as np

# ---------------------------------------------------------------------------
# SDSS DR18 Skyserver
# ---------------------------------------------------------------------------
SDSS_SQL_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"
SDSS_CHUNK_SIZE = 400
SDSS_MAX_ROWS = 10000
Z_MAX = 0.35
SN_CUT = 3.0

# ---------------------------------------------------------------------------
# Cosmology (Planck 2018, flat LambdaCDM)
# ---------------------------------------------------------------------------
H0 = 67.4          # km / s / Mpc
OMEGA_M = 0.315

# ---------------------------------------------------------------------------
# WISE
# ---------------------------------------------------------------------------
WISE_MATCH_RADIUS = 3.0   # arcsec

# ---------------------------------------------------------------------------
# DR18 galSpecLine columns for velocity dispersion
# ---------------------------------------------------------------------------
VELDISP_COL = "sigma_forbidden"
VELDISP_COL_ALT = "sigma_balmer"

# ---------------------------------------------------------------------------
# BPT -- Kewley et al. (2001) theoretical maximum starburst line
#      log([OIII]/Hb) = 0.61 / (log([NII]/Ha) - 0.47) + 1.19
# ---------------------------------------------------------------------------
def _kewley_line(log_n2ha):
    """Kewley+01: valid where log([NII]/Ha) < 0.47."""
    valid = (log_n2ha < 0.47)
    result = np.full_like(log_n2ha, np.nan)
    result[valid] = 0.61 / (log_n2ha[valid] - 0.47) + 1.19
    return result


# ---------------------------------------------------------------------------
# BPT -- Kauffmann et al. (2003) empirical division
#      log([OIII]/Hb) = 0.61 / (log([NII]/Ha) - 0.05) + 1.30
# ---------------------------------------------------------------------------
def _kauffmann_line(log_n2ha):
    """Kauffmann+03: valid where log([NII]/Ha) < 0.05."""
    valid = (log_n2ha < 0.05)
    result = np.full_like(log_n2ha, np.nan)
    result[valid] = 0.61 / (log_n2ha[valid] - 0.05) + 1.30
    return result


# Convenience: smooth x-grid for plotting
BPT_X_GRID = np.linspace(-2.0, 0.47, 500)
BPT_Y_KEWLEY = _kewley_line(BPT_X_GRID)
BPT_Y_KAUFFMANN = _kauffmann_line(BPT_X_GRID)


# ---------------------------------------------------------------------------
# WHAN diagram boundaries (Cid Fernandes et al. 2011)
# ---------------------------------------------------------------------------
WHAN_NII_HA_THRESHOLD = -0.4
WHAN_EWHA_SF_THRESHOLD = 3.0
WHAN_EWHA_WAGN_THRESHOLD = 6.0
WHAN_EWHA_RETIRED_THRESHOLD = 3.0
WHAN_EW_PASSIVE_THRESHOLD = 0.5


# ---------------------------------------------------------------------------
# WISE AGN wedge (Mateos et al. 2012, 2013)
#   y = log10(f4.6/f3.4),   x = log10(f12/f4.6)
# ---------------------------------------------------------------------------
WISE_WEDGE_SLOPE = 0.315
WISE_WEDGE_TOP_OFFSET = 0.297
WISE_WEDGE_BOT_OFFSET = -0.110
WISE_POWERLAW_SLOPE = -3.172
WISE_POWERLAW_INTERCEPT = 0.436


def is_wise_agn(x, y):
    """Return bool mask for objects inside the Mateos+2012 AGN wedge."""
    y_center = WISE_WEDGE_SLOPE * x
    y_top = y_center + WISE_WEDGE_TOP_OFFSET
    y_bot = y_center + WISE_WEDGE_BOT_OFFSET
    y_pl = WISE_POWERLAW_SLOPE * x + WISE_POWERLAW_INTERCEPT
    inside = (y >= y_bot) & (y <= y_top) & (y >= y_pl)
    return inside


# ---------------------------------------------------------------------------
# Electron density (PyNeb) -- default NLR temperature
# ---------------------------------------------------------------------------
TE_NLR = 10000.0  # K


# ---------------------------------------------------------------------------
# M-sigma relation parameters
# log(M_BH / M_sun) = alpha + beta * log10(sigma / 200 km/s)
# ---------------------------------------------------------------------------
MSIGMA_PARAMS = {
    "tremaine2002": {"alpha": 8.13, "beta": 4.02, "sigma0": 200., "label": "Tremaine+02"},
    "gebhardt2000": {"alpha": 8.08, "beta": 3.75, "sigma0": 200., "label": "Gebhardt+00"},
    "wu2009":       {"alpha": 8.09, "beta": 4.08, "sigma0": 200., "label": "Wu 09"},
}
