"""SDSS DR18 SkyServer paginated SQL queries.

DR18 column names in galSpecLine are *lowercase* and include all needed
quantities (fluxes, equivalent widths, velocity dispersions) in one table.
"""

import time
import numpy as np
import requests
from io import BytesIO
from astropy.table import Table
from .config import SDSS_SQL_URL, SDSS_CHUNK_SIZE, SDSS_MAX_ROWS, Z_MAX, SN_CUT


def _fetch_chunk(sql, fmt="csv", timeout=90, retries=10):
    """Run a single SQL query against SDSS SkyServer, return raw CSV text.

    Retries on SSL / connection errors (SkyServer is intermittent).
    """
    session = requests.Session()
    adapter = requests.adapters.HTTPAdapter(
        max_retries=3,
        pool_connections=1,
        pool_maxsize=1,
    )
    session.mount("https://", adapter)
    session.mount("http://", adapter)

    params = {"cmd": sql, "format": fmt}
    for attempt in range(retries):
        try:
            resp = session.get(SDSS_SQL_URL, params=params,
                             timeout=timeout, verify=True)
            resp.raise_for_status()
            return resp.text
        except Exception as e:
            wait = min(8 * (attempt + 1), 60)
            if attempt < retries - 1:
                time.sleep(wait)
            else:
                raise RuntimeError(
                    f"SDSS query failed after {retries} attempts: {e}"
                ) from e


def query_sdss(base_sql, max_rows=SDSS_MAX_ROWS, chunk_size=SDSS_CHUNK_SIZE,
               delay=6.0):
    """Cursor-based paginated query to SDSS SkyServer.

    DR18 does NOT support OFFSET/FETCH.  We use cursor pagination:
    ``SELECT TOP N ... WHERE ... AND s.specObjID > last``.

    The base_sql must:
      - Start with ``SELECT <cols>`` (we replace 'SELECT' with
        'SELECT TOP N')
      - Alias SpecObj as ``s``
      - NOT end with ORDER BY (added internally)
      - Put all filtering in a single WHERE clause

    Parameters
    ----------
    base_sql : str
    max_rows : int
    chunk_size : int
    delay : float
        Seconds between requests (keep >= 6 to avoid 500s).

    Returns
    -------
    astropy.table.Table
    """
    # Build a template: strip leading SELECT, we'll prepend SELECT TOP N
    inner = base_sql.strip()
    assert inner.upper().startswith("SELECT"), \
        "base_sql must start with SELECT"
    # Remove the SELECT keyword but keep the column list
    col_start = inner[6:].strip()  # after 'SELECT'

    tables = []
    fetched = 0
    last_id = 0

    while fetched < max_rows:
        chunk = min(chunk_size, max_rows - fetched)
        cursor = f" AND s.specObjID > {last_id}" if last_id > 0 else ""

        pag_sql = (
            f"SELECT TOP {chunk} {col_start}"
            f"{cursor}"
            f" ORDER BY s.specObjID"
        )

        raw = _fetch_chunk(pag_sql)

        if raw.strip() == "":
            break

        tbl = Table.read(BytesIO(raw.encode() if isinstance(raw, str) else raw),
                         format="ascii.csv", header_start=1, data_start=2)
        if len(tbl) == 0:
            break

        tables.append(tbl)
        fetched += len(tbl)
        last_id = int(tbl["specObjID"].max())
        print(f"  fetched {fetched} rows (last specObjID={last_id})")
        time.sleep(delay)

    if not tables:
        raise RuntimeError("SDSS query returned zero rows.")

    # Stack tables --- all should have same columns
    full = Table(np.hstack([t.as_array() for t in tables]),
                 names=tables[0].colnames)
    return full


def sdss_sql_narrow_line_galaxies():
    """SQL for narrow emission-line galaxies (z<0.35, S/N>=3 on 4 lines).

    DR18 galSpecLine contains everything: fluxes, EQW, sigma_balmer,
    sigma_forbidden.  No joins to galSpecIndx/galSpecInfo needed.
    """
    sql = f"""
    SELECT
        s.specObjID, s.ra, s.dec, s.z,
        g.oiii_5007_flux,  g.oiii_5007_flux_err,
        g.h_beta_flux,     g.h_beta_flux_err,
        g.nii_6584_flux,   g.nii_6584_flux_err,
        g.h_alpha_flux,    g.h_alpha_flux_err,
        g.sii_6717_flux,   g.sii_6717_flux_err,
        g.sii_6731_flux,   g.sii_6731_flux_err,
        g.oi_6300_flux,    g.oi_6300_flux_err,
        g.h_alpha_eqw,     g.h_alpha_eqw_err,
        g.sigma_balmer,    g.sigma_forbidden
    FROM SpecObj s
    JOIN galSpecLine g ON s.specObjID = g.specObjID
    WHERE
        s.z < {Z_MAX} AND s.z > 0
        AND s.class = 'GALAXY'
        AND g.h_alpha_flux > 0   AND g.h_alpha_flux_err > 0
        AND g.h_beta_flux > 0    AND g.h_beta_flux_err > 0
        AND g.oiii_5007_flux > 0 AND g.oiii_5007_flux_err > 0
        AND g.nii_6584_flux > 0  AND g.nii_6584_flux_err > 0
        AND (g.h_alpha_flux / g.h_alpha_flux_err)   >= {SN_CUT}
        AND (g.h_beta_flux / g.h_beta_flux_err)     >= {SN_CUT}
        AND (g.oiii_5007_flux / g.oiii_5007_flux_err) >= {SN_CUT}
        AND (g.nii_6584_flux / g.nii_6584_flux_err) >= {SN_CUT}
    """
    return sql
