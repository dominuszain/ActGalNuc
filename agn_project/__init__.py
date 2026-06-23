"""AGN Final Project -- SDSS + WISE emission-line galaxy classification."""

from .config import (
    H0, OMEGA_M,
    SDSS_SQL_URL, SDSS_CHUNK_SIZE, SDSS_MAX_ROWS, Z_MAX, SN_CUT,
    WISE_MATCH_RADIUS,
    VELDISP_COL, VELDISP_COL_ALT,
    TE_NLR,
    MSIGMA_PARAMS,
    _kewley_line, _kauffmann_line,
    BPT_X_GRID, BPT_Y_KEWLEY, BPT_Y_KAUFFMANN,
    WHAN_NII_HA_THRESHOLD, WHAN_EWHA_SF_THRESHOLD,
    WHAN_EWHA_WAGN_THRESHOLD, WHAN_EWHA_RETIRED_THRESHOLD,
    WHAN_EW_PASSIVE_THRESHOLD,
    WISE_WEDGE_SLOPE, WISE_WEDGE_TOP_OFFSET, WISE_WEDGE_BOT_OFFSET,
    WISE_POWERLAW_SLOPE, WISE_POWERLAW_INTERCEPT,
    is_wise_agn,
)

from .utils import cosmo

from .sdss import query_sdss, sdss_sql_narrow_line_galaxies

from .classification import (
    bpt_flux_ratios,
    classify_bpt_kewley, classify_bpt_kauffmann,
    whan_ratios, classify_whan,
    classify_wise_agn,
)

from .wise import crossmatch_wise, wise_colors

from .physics import electron_density, mbh_from_sigma, mbh_all_relations
