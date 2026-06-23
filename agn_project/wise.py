"""WISE cross-matching utilities.

Queries the SDSS DR18 SkyServer's ``wise_allsky`` table and does
positional cross-matching in Python (avoiding expensive SQL JOINs).
"""

import time
import numpy as np
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from io import BytesIO
import requests

from .config import WISE_MATCH_RADIUS, SDSS_SQL_URL


def crossmatch_wise(sdss_table, radius_arcsec=WISE_MATCH_RADIUS,
                    chunk_size=500, delay=6.0):
    """Cross-match SDSS positions with WISE via SDSS SkyServer.

    Queries ``wise_allsky`` by sky region (no SQL JOIN needed) and
    does the positional cross-match locally with astropy.

    Parameters
    ----------
    sdss_table : astropy.table.Table
        Must contain columns: 'specObjID', 'ra', 'dec'.
    radius_arcsec : float
        Matching radius in arcsec.
    chunk_size : int
        SDSS objects processed per chunk.
    delay : float
        Seconds between requests.

    Returns
    -------
    matched : astropy.table.Table
        SDSS + WISE columns.
    """
    radius_deg = radius_arcsec / 3600.0
    radius_q = radius_deg * u.deg

    all_matched = []
    n = len(sdss_table)

    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        chunk = sdss_table[start:end]

        # Bounding box around this chunk (+ small margin)
        margin = 5.0 / 3600.0  # 5 arcsec extra margin
        ra_min = float(np.min(chunk["ra"])) - radius_deg - margin
        ra_max = float(np.max(chunk["ra"])) + radius_deg + margin
        dec_min = float(np.min(chunk["dec"])) - radius_deg - margin
        dec_max = float(np.max(chunk["dec"])) + radius_deg + margin

        sql = f"""
        SELECT cntr, ra, dec, w1mpro, w2mpro, w3mpro, w4mpro
        FROM wise_allsky
        WHERE
            ra BETWEEN {ra_min} AND {ra_max}
            AND dec BETWEEN {dec_min} AND {dec_max}
        """

        for attempt in range(5):
            try:
                resp = requests.get(
                    SDSS_SQL_URL,
                    params={"cmd": sql, "format": "csv"},
                    timeout=120,
                )
                if resp.status_code == 200:
                    break
            except Exception:
                if attempt < 4:
                    time.sleep(5 * (attempt + 1))

        if resp.status_code != 200:
            print(f"  WISE query failed for chunk {start}-{end}")
            continue

        raw = resp.text.strip()
        # Response is "#Table1\ncol1,col2,...\ndata1,data2,..."
        # If only "#Table1" (no columns), skip
        if raw == "#Table1" or raw == "":
            continue
        lines = raw.split("\n")
        if len(lines) < 3:  # header + column names but no data rows
            continue

        wise_tbl = Table.read(
            BytesIO(resp.content), format="ascii.csv",
            header_start=1, data_start=2,
        )
        if len(wise_tbl) == 0:
            continue

        # Local positional cross-match
        sdss_coords = SkyCoord(ra=np.array(chunk["ra"]),
                               dec=np.array(chunk["dec"]), unit="deg")
        wise_coords = SkyCoord(ra=np.array(wise_tbl["ra"]),
                               dec=np.array(wise_tbl["dec"]), unit="deg")

        idx_wise, idx_sdss, sep2d, _ = sdss_coords.search_around_sky(
            wise_coords, radius_q
        )

        if len(idx_sdss) == 0:
            continue

        s_part = chunk[idx_sdss]
        for col in wise_tbl.colnames:
            s_part[f"w_{col}"] = wise_tbl[col][idx_wise]
        s_part["w_sep_arcsec"] = sep2d.to("arcsec").value
        all_matched.append(s_part)

        print(f"  chunk {start}-{end}: {len(idx_sdss)} matches")
        time.sleep(delay)

    if not all_matched:
        print("  No WISE matches found.")
        return Table()

    full = Table(np.hstack([t.as_array() for t in all_matched]),
                 names=all_matched[0].colnames)

    # Keep closest match per SDSS object
    if "w_sep_arcsec" in full.colnames and "specObjID" in full.colnames:
        ids = full["specObjID"]
        seps = full["w_sep_arcsec"]
        order = np.argsort(seps)
        full = full[order]
        _, keep = np.unique(ids[order], return_index=True)
        full = full[keep]

    return full


def wise_colors(table):
    """Compute WISE colour indices.

    Expects columns: w_w1mpro, w_w2mpro, w_w3mpro.

    Returns x = log10(f12/f4.6), y = log10(f4.6/f3.4).
    """
    w1 = table["w_w1mpro"]
    w2 = table["w_w2mpro"]
    w3 = table["w_w3mpro"]

    f1 = 10 ** (-w1 / 2.5)
    f2 = 10 ** (-w2 / 2.5)
    f3 = 10 ** (-w3 / 2.5)

    x = np.log10(f3 / f2)
    y = np.log10(f2 / f1)
    return x, y
