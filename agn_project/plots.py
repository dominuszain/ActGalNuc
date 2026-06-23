"""Plotting functions for BPT, WHAN, WISE, density, and kinematics diagrams."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm
from scipy.stats import binned_statistic_2d

from .config import (
    BPT_X_GRID,
    BPT_Y_KEWLEY,
    BPT_Y_KAUFFMANN,
    WISE_WEDGE_SLOPE,
    WISE_WEDGE_TOP_OFFSET,
    WISE_WEDGE_BOT_OFFSET,
    WISE_POWERLAW_SLOPE,
    WISE_POWERLAW_INTERCEPT,
    WHAN_NII_HA_THRESHOLD,
    WHAN_EWHA_SF_THRESHOLD,
    WHAN_EWHA_WAGN_THRESHOLD,
    WHAN_EW_PASSIVE_THRESHOLD,
)

# Consistent colour palette across all plots
COLORS = {
    "Star-forming": "#1f77b4",  # blue
    "Composite":    "#7f7f7f",  # grey
    "AGN/Seyfert":  "#d62728",  # red
    "AGN":          "#d62728",  # red
    "Unclassified": "#bcbd22",  # olive
    "Strong AGN":   "#d62728",
    "wAGN":         "#ff7f0e",  # orange
    "Retired":      "#2ca02c",  # green
    "Passive":      "#9467bd",  # purple
    "HII":          "#1f77b4",
}

MARKER = "o"
MS = 2.5
ALPHA = 0.6


def _init_figure(figsize=(8, 7)):
    fig, ax = plt.subplots(figsize=figsize)
    ax.tick_params(direction="in", top=True, right=True)
    return fig, ax


def plot_bpt(log_n2ha, log_o3hb, labels, title="BPT Diagram",
             fname="output/bpt_diagram.pdf",
             kewley_label="Kewley+01", kauff_label="Kauffmann+03"):
    """NII BPT diagram with classification curves.

    Parameters
    ----------
    log_n2ha, log_o3hb : array
        log10 flux ratios.
    labels : array of str
        Classification label per object (e.g. 'Star-forming', 'Composite',
        'AGN/Seyfert').
    """
    fig, ax = _init_figure()

    # Plot data points class by class for a clean legend
    unique_labels = sorted(set(labels), key=lambda x: list(COLORS).index(x)
                            if x in COLORS else 99)
    for lab in unique_labels:
        mask = labels == lab
        ax.scatter(log_n2ha[mask], log_o3hb[mask],
                   c=COLORS.get(lab, COLORS["Unclassified"]),
                   marker=MARKER, s=MS, alpha=ALPHA, label=lab,
                   edgecolors="none", rasterized=True)

    # Classification curves
    ax.plot(BPT_X_GRID, BPT_Y_KEWLEY, "k-", lw=1.5, label=kewley_label)
    ax.plot(BPT_X_GRID, BPT_Y_KAUFFMANN, "k--", lw=1.5, label=kauff_label)

    ax.set_xlabel(r"$\log_{10}\,([{\rm N\,II]}\,\lambda 6583 / {\rm H}\alpha)$")
    ax.set_ylabel(r"$\log_{10}\,([{\rm O\,III]}\,\lambda 5007 / {\rm H}\beta)$")
    ax.set_title(title)
    ax.legend(loc="upper left", fontsize=7, markerscale=3)
    ax.set_xlim(-2.0, 0.8)
    ax.set_ylim(-1.5, 1.5)
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f"  saved {fname}")


# =========================================================================
# WHAN diagram
# =========================================================================

def plot_whan(log_n2ha, ew_ha, labels, title="WHAN Diagram",
              fname="output/whan_diagram.pdf"):
    """WHAN diagram: EW(Ha) vs log([NII]/Ha) with log y-axis."""
    fig, ax = _init_figure()

    # Log scale: clip very small/negative EWs to a floor
    ew_plot = np.clip(ew_ha, 0.05, None)

    unique_labels = sorted(set(labels), key=lambda x: list(COLORS).index(x)
                            if x in COLORS else 99)
    for lab in unique_labels:
        mask = labels == lab
        ax.scatter(log_n2ha[mask], ew_plot[mask],
                   c=COLORS.get(lab, COLORS["Unclassified"]),
                   marker=MARKER, s=MS, alpha=ALPHA, label=lab,
                   edgecolors="none", rasterized=True)

    ax.set_yscale("log")

    # Classification boundaries (horizontal lines on log scale)
    xlim = ax.get_xlim()
    ax.axvline(WHAN_NII_HA_THRESHOLD, color="k", ls="--", lw=1.2,
               label=r"$\log([{\rm N\,II}]/{\rm H}\alpha) = -0.4$")
    ax.axhline(WHAN_EWHA_SF_THRESHOLD, color="grey", ls=":", lw=1.2)
    ax.axhline(WHAN_EWHA_WAGN_THRESHOLD, color="grey", ls="-.", lw=1.2)
    ax.axhline(WHAN_EW_PASSIVE_THRESHOLD, color="grey", ls="--", lw=1.2)

    ax.set_xlabel(r"$\log_{10}\,([{\rm N\,II]}\,\lambda 6583 / {\rm H}\alpha)$")
    ax.set_ylabel(r"$W_{{\rm H}\alpha}$ (\AA, log scale)")
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=7, markerscale=3)
    ax.set_xlim(xlim)
    ax.set_ylim(bottom=0.03, top=3e3)
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f"  saved {fname}")


# =========================================================================
# WISE colour-colour diagram
# =========================================================================

def plot_wise(x, y, agn_mask=None, title="WISE Colour-Colour Diagram",
              fname="output/wise_diagram.pdf"):
    """WISE log(f12/f4.6) vs log(f4.6/f3.4) with Mateos+2012 AGN wedge.

    Parameters
    ----------
    x : array --- log10(f12/f4.6)
    y : array --- log10(f4.6/f3.4)
    agn_mask : bool array or None
        If given, colour points by AGN / non-AGN.
    """
    fig, ax = _init_figure()

    if agn_mask is None:
        agn_mask = np.zeros(len(x), dtype=bool)

    x_grid = np.linspace(-0.5, 1.8, 200)
    y_center = WISE_WEDGE_SLOPE * x_grid
    y_top = y_center + WISE_WEDGE_TOP_OFFSET
    y_bot = y_center + WISE_WEDGE_BOT_OFFSET
    y_pl = WISE_POWERLAW_SLOPE * x_grid + WISE_POWERLAW_INTERCEPT

    # Mateos+12 AGN selection wedge
    ax.fill_between(x_grid, y_bot, y_top, color="grey", alpha=0.12,
                    label="Mateos+12 wedge")
    ax.plot(x_grid, y_center, "-", color="grey", lw=1.2, alpha=0.8,
            label=r"Wedge centre ($y=0.315x$)")
    ax.plot(x_grid, y_top, "--", color="red", lw=0.8, alpha=0.7,
            label=r"Top boundary ($+0.297$)")
    ax.plot(x_grid, y_bot, "--", color="blue", lw=0.8, alpha=0.7,
            label=r"Bottom boundary ($-0.110$)")
    ax.plot(x_grid, y_pl, "k-.", lw=1, label=r"$\alpha=-0.3$ limit")

    # Data
    ax.scatter(x[~agn_mask], y[~agn_mask], c="blue", s=MS, alpha=0.4,
               edgecolors="none", rasterized=True, label="non-AGN")
    ax.scatter(x[agn_mask], y[agn_mask], c="red", s=MS, alpha=0.6,
               edgecolors="none", rasterized=True, label="WISE AGN")

    ax.set_xlabel(r"$\log_{10}\,(f_{12\mu{\rm m}} / f_{4.6\mu{\rm m}})$")
    ax.set_ylabel(r"$\log_{10}\,(f_{4.6\mu{\rm m}} / f_{3.4\mu{\rm m}})$")
    ax.set_title(title)
    ax.legend(fontsize=7, markerscale=3)
    ax.set_xlim(-0.6, 1.8)
    ax.set_ylim(-0.6, 0.6)
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f"  saved {fname}")


# =========================================================================
# 2D BPT density/kinesis maps
# =========================================================================

def plot_bpt_2d_map(log_n2ha, log_o3hb, values, cbar_label, title,
                    fname, nbins=30, cmap="viridis", log_norm=False):
    """Plot a 2D histogram or mean-value map on the BPT plane.

    Parameters
    ----------
    values : array
        Quantity to colour-code (e.g. electron density or vel dispersion).
    """
    fig, ax = _init_figure()

    norm = LogNorm() if log_norm else None

    stat, x_edges, y_edges, _ = binned_statistic_2d(
        log_n2ha, log_o3hb, values, statistic="median",
        bins=[np.linspace(-2, 0.8, nbins), np.linspace(-1.5, 1.5, nbins)]
    )

    xi = 0.5 * (x_edges[1:] + x_edges[:-1])
    yi = 0.5 * (y_edges[1:] + y_edges[:-1])
    mesh = ax.pcolormesh(x_edges, y_edges, stat.T, cmap=cmap, norm=norm,
                         shading="flat", rasterized=True)

    # Overlay classification curves
    ax.plot(BPT_X_GRID, BPT_Y_KEWLEY, "k-", lw=1)
    ax.plot(BPT_X_GRID, BPT_Y_KAUFFMANN, "k--", lw=1)

    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(cbar_label)

    ax.set_xlabel(r"$\log_{10}\,([{\rm N\,II}] / {\rm H}\alpha)$")
    ax.set_ylabel(r"$\log_{10}\,([{\rm O\,III}] / {\rm H}\beta)$")
    ax.set_title(title)
    ax.set_xlim(-2, 0.8)
    ax.set_ylim(-1.5, 1.5)
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f"  saved {fname}")


# =========================================================================
# Histogram helper
# =========================================================================

def plot_histogram(values, bins=30, xlabel="", title="",
                   fname="output/histogram.pdf", log=False):
    fig, ax = _init_figure((7, 4))
    if log:
        ax.set_yscale("log")
        ax.hist(values, bins=bins, color="#1f77b4", edgecolor="k", lw=0.3)
        ax.set_ylabel("Count (log scale)")
    else:
        ax.hist(values, bins=bins, color="#1f77b4", edgecolor="k", lw=0.3)
        ax.set_ylabel("Count")
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f"  saved {fname}")


def plot_histogram_by_class(values, labels, classes, bins=30, xlabel="",
                            title="", fname="output/histogram.pdf", log=False):
    """Overlaid histograms for a list of class labels."""
    fig, ax = _init_figure((8, 5))
    for cls in classes:
        mask = labels == cls
        if mask.sum() == 0:
            continue
        ax.hist(values[mask], bins=bins, alpha=0.5,
                color=COLORS.get(cls, "#bcbd22"), label=cls,
                edgecolor="k", lw=0.3)
    if log:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    ax.set_title(title)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f"  saved {fname}")
