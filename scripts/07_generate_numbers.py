#!/usr/bin/env python3
"""Generate numbers.tex -- all report statistics from FITS data.

Run this whenever the FITS data changes. The report uses
\\input{numbers.tex} to pull in all values automatically.
All command names are letter-only (no digits) for TeX compatibility.
"""

import sys
sys.path.insert(0, "/home/zain/Desktop/test")

import numpy as np
from astropy.table import Table
from agn_project.physics import mbh_from_sigma

DATA = "/home/zain/Desktop/test/data"
OUT = "/home/zain/Desktop/test/report/numbers.tex"

def pct(n, d):
    return f"{100*n/d:.1f}"

def main():
    tbl = Table.read(f"{DATA}/sdss_sample_classified.fits")
    tbl_s = Table.read(f"{DATA}/sdss_sample_sii.fits")
    tbl_w = Table.read(f"{DATA}/wise_matched.fits")
    N = len(tbl)

    # ---- helpers ----
    bpt = tbl["bpt_kewley"]
    wh  = tbl["whan_label"]
    bpt_agn = bpt == "AGN/Seyfert"
    wh_agn  = (wh == "Strong AGN") | (wh == "wAGN")

    cmds = {}

    # Sample
    cmds["samplesize"] = N

    # BPT
    cmds["bptSfCount"]      = int(np.sum(bpt == "Star-forming"))
    cmds["bptCompCount"]    = int(np.sum(bpt == "Composite"))
    cmds["bptAgnCount"]     = int(np.sum(bpt == "AGN/Seyfert"))
    cmds["bptSfPct"]        = pct(cmds["bptSfCount"], N)
    cmds["bptCompPct"]      = pct(cmds["bptCompCount"], N)
    cmds["bptAgnPct"]       = pct(cmds["bptAgnCount"], N)
    cmds["kauffSfCount"]    = int(np.sum(tbl["bpt_kauff"] == "Star-forming"))
    cmds["kauffAgnCount"]   = int(np.sum(tbl["bpt_kauff"] == "AGN"))
    cmds["kauffAgnPct"]     = pct(cmds["kauffAgnCount"], N)

    # WHAN
    for label, key in [("Star-forming","whanSf"),("Strong AGN","whanStrong"),
                       ("wAGN","whanWagn"),("Retired","whanRet"),("Passive","whanPass")]:
        cmds[key+"Count"] = int(np.sum(wh == label))
        cmds[key+"Pct"]   = pct(cmds[key+"Count"], N)
    cmds["whanAgnTotCount"] = cmds["whanStrongCount"] + cmds["whanWagnCount"]
    cmds["whanAgnTotPct"]   = pct(cmds["whanAgnTotCount"], N)

    # Cross-comparison
    cmds["bothAgnCount"]    = int(np.sum(bpt_agn & wh_agn))
    cmds["bptOnlyCount"]    = int(np.sum(bpt_agn & ~wh_agn))
    cmds["whanOnlyCount"]   = int(np.sum(~bpt_agn & wh_agn))
    cmds["bothAgnPct"]      = pct(cmds["bothAgnCount"], N)
    cmds["bptOnlyPct"]      = pct(cmds["bptOnlyCount"], N)
    cmds["whanOnlyPct"]     = pct(cmds["whanOnlyCount"], N)

    # WHAN-only breakdown
    wonly = ~bpt_agn & wh_agn
    cmds["woStrongCount"] = int(np.sum(wh[wonly] == "Strong AGN"))
    cmds["woWagnCount"]   = int(np.sum(wh[wonly] == "wAGN"))
    cmds["woStrongPct"]   = pct(cmds["woStrongCount"], cmds["whanOnlyCount"])
    cmds["woWagnPct"]     = pct(cmds["woWagnCount"], cmds["whanOnlyCount"])
    cmds["woBptSfCount"]  = int(np.sum((bpt=="Star-forming") & wonly))
    cmds["woBptCompCount"]= int(np.sum((bpt=="Composite") & wonly))
    cmds["woBptSfPct"]    = pct(cmds["woBptSfCount"], cmds["whanOnlyCount"])
    cmds["woBptCompPct"]  = pct(cmds["woBptCompCount"], cmds["whanOnlyCount"])

    # BPT-only breakdown
    bonly = bpt_agn & ~wh_agn
    cmds["boRetCount"] = int(np.sum(wh[bonly] == "Retired"))
    cmds["boSfCount"]  = int(np.sum(wh[bonly] == "Star-forming"))
    cmds["boPassCount"]= int(np.sum(wh[bonly] == "Passive"))
    cmds["boRetPct"]   = pct(cmds["boRetCount"], cmds["bptOnlyCount"])
    cmds["boSfPct"]    = pct(cmds["boSfCount"], cmds["bptOnlyCount"])
    cmds["boPassPct"]  = pct(cmds["boPassCount"], cmds["bptOnlyCount"])

    # WISE
    cmds["wiseAgnSample"]  = 1436
    cmds["wiseMatched"]    = len(tbl_w)
    cmds["wiseMatchRate"]  = pct(len(tbl_w), 1436)
    from agn_project.config import is_wise_agn
    from agn_project.wise import wise_colors
    x, y = wise_colors(tbl_w)
    n_wedge = int(np.sum(is_wise_agn(x, y)))
    cmds["wiseWedge"]      = n_wedge
    cmds["wiseWedgePct"]   = pct(n_wedge, len(tbl_w))
    cmds["wiseOutside"]    = len(tbl_w) - n_wedge
    cmds["wiseOutsidePct"] = pct(cmds["wiseOutside"], len(tbl_w))

    # SII / sigma
    sii_r = tbl_s["sii_6717_flux"] / tbl_s["sii_6731_flux"]
    cmds["siiSample"]      = len(tbl_s)
    cmds["siiAboveLimit"]  = int(np.sum(sii_r > 1.45))
    cmds["siiAboveLimitPct"] = pct(cmds["siiAboveLimit"], len(tbl_s))
    from agn_project.physics import electron_density
    ne = electron_density(sii_r)
    cmds["neValid"]        = int(np.sum(np.isfinite(ne)))
    sigma = tbl_s["sigma_forbidden"]
    cmds["sigmaMedian"]    = int(round(np.median(sigma)))
    cmds["sigmaMin"]       = int(round(sigma.min()))
    cmds["sigmaMax"]       = int(round(sigma.max()))
    cmds["sigmaSubTen"]    = int(np.sum(sigma < 10))

    # MBH
    logT = mbh_from_sigma(sigma, "tremaine2002")
    logG = mbh_from_sigma(sigma, "gebhardt2000")
    logW = mbh_from_sigma(sigma, "wu2009")
    mask = lambda x: np.isfinite(x) & (x > 4) & (x < 10)
    cmds["mbhNtrem"]   = int(np.sum(mask(logT)))
    cmds["mbhNgeb"]    = int(np.sum(mask(logG)))
    cmds["mbhNwu"]     = int(np.sum(mask(logW)))
    cmds["mbhMedTrem"] = f"{np.median(logT[mask(logT)]):.2f}"
    cmds["mbhMedGeb"]  = f"{np.median(logG[mask(logG)]):.2f}"
    cmds["mbhMedWu"]   = f"{np.median(logW[mask(logW)]):.2f}"
    cmds["mbhBelowFour"] = int(np.sum(logT <= 4))
    for lo, hi, k in [(5,6,"FiveSix"),(6,7,"SixSeven"),(7,8,"SevenEight"),(8,9,"EightNine")]:
        cmds[f"mbhBin{k}Count"] = int(np.sum((logT>=lo)&(logT<hi)))
        cmds[f"mbhBin{k}Pct"]   = pct(cmds[f"mbhBin{k}Count"], cmds["mbhNtrem"])
    other = int(np.sum(mask(logT) & ((logT<5)|(logT>=9))))
    cmds["mbhBinOtherCount"] = other
    cmds["mbhBinOtherPct"]   = pct(other, cmds["mbhNtrem"])

    # ---- write ----
    lines = ["% Auto-generated by scripts/07_generate_numbers.py -- do not edit"]
    for name, val in cmds.items():
        lines.append(f"\\newcommand{{\\{name}}}{{{val}}}")
    with open(OUT, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Written {len(cmds)} LaTeX commands to {OUT}")

if __name__ == "__main__":
    main()
