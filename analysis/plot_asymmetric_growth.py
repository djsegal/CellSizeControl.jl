#!/usr/bin/env python3
"""Publication-grade asymmetric-growth figure for the aging-unification paper.

(a) The volume-versus-time lineage over a replicative lifespan -- the mother body is
monotonic (never shrinks; only the bud detaches) and enlarges with age, while each
bud grows and pinches off as a daughter. (b) Per generation, the model set-point
volumes (markers) with the closed-form laws the mechanism follows overlaid: the
mother enlargement is a single saturating exponential V*(g) and the daughter birth
size the product of two saturating processes (asymmetry erosion x mother
enlargement). No goodness-of-fit statistic is shown: the markers are noiseless
deterministic output, so a smooth-curve R^2 is tautologically ~1 and uninformative;
the panel reports the functional form, not a fit quality.

Okabe-Ito palette, markers with the closed-form overlay, 300 dpi. Data:
lineage_timecourse.csv + cs_da_lineage.csv (rigorous_cell_size.jl). Run via a venv
with matplotlib + scipy.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit

from _pubstyle import (apply_style, pub_arrow, halo, opaque_legend, pub_audit,
                       BLUE, VERM, GREEN)

HERE = Path(__file__).resolve().parent
# Okabe-Ito colorblind-safe palette (science-space standard): mother = blue,
# daughter = vermillion.
RLS = 30  # replicative-lifespan horizon (divisions); budding yeast ~25-30 (Schnitzer 2022)


def _r2(y, yhat):
    y = np.asarray(y, float)
    return 1.0 - np.sum((y - yhat) ** 2) / np.sum((y - y.mean()) ** 2)


def single_exp_fit(x, y):
    """Mother enlargement V*(g) = a + b*exp(c*g) -- the model's set-point law."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    f = lambda g, a, b, c: a + b * np.exp(c * g)
    p, _ = curve_fit(f, x, y, p0=[y[-1], y[0] - y[-1], -0.1], maxfev=40000)
    xs = np.linspace(x.min(), x.max(), 300)
    return xs, f(xs, *p), p, _r2(y, f(x, *p))


def product_fit(x, y):
    """Daughter V_d(g) = v0 + A*(1-exp(-g/t1))*(1+k*(1-exp(-g/t2))) -- the product of
    the asymmetry erosion beta(g) and the mother enlargement V*(g), so it carries two
    timescales rather than one."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    f = lambda g, v0, A, t1, k, t2: v0 + A * (1 - np.exp(-g / t1)) * (1 + k * (1 - np.exp(-g / t2)))
    p, _ = curve_fit(f, x, y, p0=[20, 15, 12, 0.45, 8], maxfev=200000)
    xs = np.linspace(x.min(), x.max(), 300)
    return xs, f(xs, *p), p, _r2(y, f(x, *p))


def panel_timecourse(ax):
    t, vm, vt = [], [], []
    with open(HERE / "lineage_timecourse.csv") as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            m = float(row[1]) * 1e15  # L -> fL
            b = float(row[2]) * 1e15
            t.append(float(row[0]) / 60.0)  # s -> min
            vm.append(m)
            vt.append(m + b)
    t = np.asarray(t)
    vm = np.asarray(vm)
    vt = np.asarray(vt)
    ax.fill_between(t, 0, vm, color=BLUE, alpha=0.22, linewidth=0)
    ax.fill_between(t, vm, vt, color=VERM, alpha=0.32, linewidth=0)
    ax.plot(t, vm, color=BLUE, lw=1.3, label="Mother body")
    ax.plot(t, vt, color=VERM, lw=1.0, label="Mother + bud")
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Volume (fL)")
    ax.set_title("(a) Lineage growth over the replicative lifespan")
    ax.set_xlim(0, t.max())
    ax.set_ylim(0, vt.max() * 1.14)
    opaque_legend(ax, loc="upper left", fontsize=11, handlelength=1.4)
    # the monotonic-mother invariant, stated where it reads (haloed so the fill can't cut it)
    # near-black for contrast on the light-blue fill (printer-grayscale / CVD / dim-screen safe)
    halo(ax.annotate("Mother body never shrinks\n(only the bud detaches)",
                     xy=(t.max() * 0.34, vm[int(len(vm) * 0.34)] * 0.55),
                     fontsize=11, color="0.12", ha="center", va="center"))


def panel_maternal(ax):
    gen, dau, mom = [], [], []
    with open(HERE / "cs_da_lineage.csv") as f:
        for row in csv.DictReader(f):
            g = int(row["gen"])
            if g > RLS:
                continue
            gen.append(g)
            dau.append(float(row["Vdaughter"]))
            mom.append(float(row["Vmother"]))
    xm, ym, pm, r2m = single_exp_fit(gen, mom)
    xd, yd, pd, r2d = product_fit(gen, dau)
    # data markers (the model output) with the mechanistic fit overlaid -- show both
    ax.scatter(gen, mom, s=26, facecolors="white", edgecolors=BLUE, linewidths=1.1,
               zorder=4, label="Mother size at Start")
    ax.scatter(gen, dau, s=22, facecolors=VERM, edgecolors="white", linewidths=0.5,
               zorder=4, label="Daughter birth size")
    ax.plot(xm, ym, "-", color=BLUE, lw=1.8, zorder=3, solid_capstyle="round")
    ax.plot(xd, yd, "-", color=VERM, lw=1.8, zorder=3, solid_capstyle="round")
    lo = min(min(dau), min(mom)) - 2
    hi = max(max(dau), max(mom)) + 18  # headroom so the upper-left annotation clears the curves
    ax.set_ylim(lo, hi)
    ax.set_xlim(0, RLS + 1)
    # the markers are deterministic model output; the lines are the closed-form laws the
    # mechanism follows (a saturating exponential; a product of two saturating processes).
    # No R^2 is reported -- fitting a smooth curve to noiseless output is tautologically ~1.
    a, b, c = pm
    # Mother eq lifted into the clear band above the mother curve; daughter eq moved left, off the
    # lower-right legend. Each label is colour-matched to its curve (no arrow needed).
    halo(ax.text(0.045, 0.646,
                 f"Mother: $V_m = {a:.0f} - {abs(b):.0f}\\,e^{{{c:.2g}\\,g}}$",
                 transform=ax.transAxes, color=BLUE, fontsize=11, va="top", ha="left"))
    halo(ax.text(0.745, 0.06,
                 "Daughter: $V_d = r(g)\\,V_m,\\ r{:}\\,0.7{\\to}0.9$",
                 transform=ax.transAxes, color=VERM, fontsize=11, va="bottom", ha="right"))
    # Published DIRECTION of the trend: daughters of older mothers are born larger
    # (Johnston 1966, Antonie van Leeuwenhoek 32:94; Yang et al. 2011, Cell Cycle 10:144).
    # We annotate the DIRECTION with the model's OWN honest fold-changes -- NOT Yang's ~4x
    # MOTHER-enlargement number, which belongs to the mother and would be mislabelled if put
    # on the daughter. In this minimal model the daughter and mother set-points are coupled
    # (V_dau = r(g) V*(g)), so it cannot reach a ~4x mother without unphysical daughters; over
    # the lifespan the daughter birth size rises ~1.8x and the mother ~1.4x.
    g_early, g_late = min(gen), max(gen)
    d_early, d_late = dau[gen.index(g_early)], dau[gen.index(g_late)]
    m_early, m_late = mom[gen.index(g_early)], mom[gen.index(g_late)]
    fold_d, fold_m = d_late / d_early, m_late / m_early
    # text in the open upper band (above both saturating curves, which top out ~52 fL),
    # arrow down to the late daughter marker; haloed so the y-gridlines can't cut the glyphs.
    # Annotate the DIRECTION (text in the clear band above the curves, each line its own short
    # haloed label so no single bbox spans the rising mother curve) + a bare arrow to the late
    # daughter marker.
    # One boxed comment in the clear upper band (no arrow: it states a general direction the rising
    # daughter curve already shows, and a cross-panel arrow only crowded the data).
    ax.text(0.629, 0.944,
            "Old mothers make larger daughters\n(direction: Johnston 1966; Yang 2011)\n"
            f"Model: {fold_d:.1f}× daughter, {fold_m:.1f}× mother",
            transform=ax.transAxes, color="0.15", fontsize=12, va="top", ha="center",
            linespacing=1.35,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="0.6", alpha=0.95))
    ax.set_xlabel("Maternal replicative age (generations)")
    ax.set_ylabel("Volume (fL)")
    ax.set_title("(b) Maternal-age asymmetry, to the replicative lifespan")
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.grid(axis="y", which="major", color="0.9", lw=0.7)
    ax.set_axisbelow(True)
    opaque_legend(ax, loc="lower right", fontsize=12, handletextpad=0.4,
                  bbox_to_anchor=(1.0, 0.16))
    return r2m, r2d


def main():
    apply_style()
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.0, 4.8))
    panel_timecourse(axL)
    r2m, r2d = panel_maternal(axR)
    fig.tight_layout(w_pad=2.0)
    issues = pub_audit(fig)
    assert not issues, "asymmetric_growth pub_audit: " + "; ".join(issues)
    out = HERE / "asymmetric_growth.png"
    fig.savefig(out, bbox_inches="tight")
    print(f"wrote {out}  (mother R2={r2m:.5f}, daughter R2={r2d:.5f}) | audit clean")


if __name__ == "__main__":
    main()
