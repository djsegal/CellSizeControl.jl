#!/usr/bin/env python3
"""The SHAPE of the daughter-RLS-vs-maternal-age prediction is convex, not linear.

The two-bucket Kennedy fold only compares endpoints. The full curve carries a sharper,
falsifiable claim: daughter replicative lifespan declines with maternal age along a *convex*
(decelerating) path -- steepest at young/mid maternal age, flattening toward old age -- and so
lies below the straight chord joining its endpoints. This is the shape directly extractable from
modern single-mother-resolved microfluidic RLS studies (Lee 2012; Jo 2015), a stronger test than
the pooled two-bucket fold.

Reads daughter_rls_convexity.csv / daughter_rls_convexity_summary.csv (from
gen_daughter_rls_convexity.jl, which builds the curve from the CellSizeControl package with
nothing refit). Okabe-Ito palette. Run via a venv with matplotlib + numpy.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _pubstyle import apply_style, halo, opaque_legend, pub_audit, BLUE, VERM

HERE = Path(__file__).resolve().parent


def _read_curve():
    x, y, chord = [], [], []
    with open(HERE / "daughter_rls_convexity.csv") as f:
        for row in csv.DictReader(f):
            x.append(float(row["frac_mid"]))
            y.append(float(row["daughter_rls_mean"]))
            chord.append(float(row["chord"]))
    return np.array(x), np.array(y), np.array(chord)


def _read_summary():
    out = {}
    with open(HERE / "daughter_rls_convexity_summary.csv") as f:
        for row in csv.DictReader(f):
            out[row["quantity"]] = float(row["value"])
    return out


def main():
    x, y, chord = _read_curve()
    s = _read_summary()

    apply_style()
    fig, ax = plt.subplots(figsize=(7.6, 5.2))

    # the convexity: shade the sag between the curve and its endpoint chord
    ax.fill_between(x, y, chord, color=BLUE, alpha=0.12, lw=0,
                    label="convexity (curve below chord)")
    # the straight-line null (a linear decline would trace this)
    ax.plot([x[0], x[-1]], [y[0], y[-1]], "--", color="0.45", lw=1.6,
            label="linear-decline null (endpoint chord)")
    # the predicted curve
    ax.plot(x, y, "-o", color=BLUE, lw=2.3, ms=5.5, mfc=BLUE, mec="white", mew=0.6,
            label="Model: daughter RLS vs maternal age\n(CellSizeControl, nothing refit)")

    # mark the peak sag (depth of convexity)
    isag = int(np.argmin(y - chord))
    ax.annotate("", xy=(x[isag], y[isag]), xytext=(x[isag], chord[isag]),
                arrowprops=dict(arrowstyle="<->", color="0.30", lw=1.3))
    halo(ax.text(x[isag] + 0.02, (y[isag] + chord[isag]) / 2,
                 f"peak sag\n{s['peak_sag_divisions']:.1f} div",
                 color="0.25", fontsize=11, ha="left", va="center"))

    ax.set(xlabel="Maternal replicative age (fraction of mother's lifespan)",
           ylabel="Daughter replicative lifespan (divisions)",
           title="Daughter lifespan declines convexly with maternal age")
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, max(y) * 1.18)
    ax.grid(axis="y", which="major", color="0.9", lw=0.7)
    ax.set_axisbelow(True)

    # honest characterization box
    txt = (f"convex (decelerating), not linear\n"
           f"$x^2$ coeff $= {s['quad_x2_coeff']:+.1f}$ (>0)\n"
           f"curvature removes "
           f"{100 * (1 - s['rms_quadratic'] / s['rms_linear']):.0f}% of the linear residual\n"
           f"steepest decline at frac $\\approx$ {s['steepest_decline_frac']:.2f} "
           f"(early/mid, not a late cliff)")
    ax.text(0.97, 0.97, txt, transform=ax.transAxes, ha="right", va="top", fontsize=10.5,
            linespacing=1.5,
            bbox=dict(boxstyle="round,pad=0.45", facecolor="white", edgecolor="0.6", alpha=0.95))
    opaque_legend(ax, loc="lower left", fontsize=11)

    fig.tight_layout()
    issues = pub_audit(fig)
    assert not issues, "convexity pub_audit: " + "; ".join(issues)
    out = HERE / "daughter_rls_convexity.png"
    fig.savefig(out, bbox_inches="tight")
    print(f"wrote {out}")
    print(f"  young {s['young_rls']:.1f} -> old {s['old_rls']:.1f} div (fold {s['fold']:.2f}x); "
          f"convex below chord, peak sag {s['peak_sag_divisions']:.2f} div; audit clean")


if __name__ == "__main__":
    main()
