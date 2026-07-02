#!/usr/bin/env python3
"""Asymmetric cell growth from our energetic cell-size model, run to the budding-yeast
replicative lifespan (~30 divisions), in the class email's visual style:

  (a) the Volume-vs-time lineage: the per-compartment growth fix keeps the bud (daughter)
      viable while the mother body enlarges with replicative age (the sawtooth of division);
  (b) the same lineage per generation: mother size and daughter birth size both rise with
      maternal age, each fit by its mechanistic form (mother = a single saturating
      exponential; daughter = bud + beta(g)*mother, a product of two saturating processes).

Data: lineage_timecourse.csv + cs_da_lineage.csv (rigorous_cell_size.jl). Run via rr_env.
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

from _pubstyle import apply_style, BLUE, VERM, GREEN, OKABE

HERE = Path(__file__).resolve().parent
# Okabe-Ito colorblind-safe palette (PubPlots / science-space standard), chosen for high
# contrast + legibility: mother = blue, daughter = VERMILLION (red-orange, not the hard-to-
# read yellow-orange). Fills are medium-saturation tints of the same hues (not washed out).
MOTHER_FILL, MOTHER_EDGE = "#7fb3d5", BLUE
BUD_FILL, BUD_EDGE = "#e8895a", VERM
RLS = 30  # replicative lifespan horizon (divisions); yeast ~25-30 (Schnitzer 2022)


def _r2(y, yhat):
    y = np.asarray(y, float)
    return 1.0 - np.sum((y - yhat) ** 2) / np.sum((y - y.mean()) ** 2)


def single_exp_fit(x, y):
    """Mother enlargement V*(g) = a + b*exp(c*g) — exactly the model's set-point law."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    f = lambda g, a, b, c: a + b * np.exp(c * g)
    p, _ = curve_fit(f, x, y, p0=[y[-1], y[0] - y[-1], -0.1], maxfev=40000)
    xs = np.linspace(x.min(), x.max(), 300)
    return xs, f(xs, *p), p, _r2(y, f(x, *p))


def product_fit(x, y):
    """Daughter V_d(g) = v0 + A*(1-exp(-g/t1))*(1+k*(1-exp(-g/t2))) — the exact product of
    the asymmetry erosion beta(g) and the mother enlargement V*(g) (two saturating
    processes), so it carries two timescales rather than one."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    f = lambda g, v0, A, t1, k, t2: v0 + A * (1 - np.exp(-g / t1)) * (1 + k * (1 - np.exp(-g / t2)))
    p, _ = curve_fit(f, x, y, p0=[20, 15, 12, 0.45, 8], maxfev=200000)
    xs = np.linspace(x.min(), x.max(), 300)
    return xs, f(xs, *p), p, _r2(y, f(x, *p))


def mother_eq(p):
    a, b, c = p
    return f"$V_m = {a:.0f} - {abs(b):.0f}\\,\\mathrm{{exp}}({c:.2g}\\,g)$"


def panel_timecourse(ax):
    t, vm, vt = [], [], []
    with open(HERE / "lineage_timecourse.csv") as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            t.append(float(row[0]) / 60.0)          # s -> min
            m = float(row[1]) * 1e15                 # L -> fL
            b = float(row[2]) * 1e15
            vm.append(m)
            vt.append(m + b)
    ax.fill_between(t, 0, vm, color=MOTHER_FILL, label="Mother", linewidth=0)
    ax.fill_between(t, vm, vt, color=BUD_FILL, label="Bud (daughter)", linewidth=0)
    ax.plot(t, vm, color=MOTHER_EDGE, lw=0.9)
    ax.plot(t, vt, color=BUD_EDGE, lw=0.9)
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Volume (fL)")
    ax.set_title("(a) Lineage growth over the replicative lifespan", fontsize=11)
    ax.set_xlim(0, max(t))
    ax.set_ylim(0, max(vt) * 1.12)
    ax.legend(loc="upper left", frameon=False, fontsize=12)


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
    ax.plot(xm, ym, "-", color=MOTHER_EDGE, lw=2.0, solid_capstyle="round")
    ax.plot(xd, yd, "-", color=BUD_EDGE, lw=2.0, solid_capstyle="round")
    lo = min(min(dau), min(mom)) - 2
    hi = max(max(dau), max(mom)) + 3
    ax.set_ylim(lo, hi)
    ax.set_xlim(0, RLS + 1)
    # direct, color-matched labels in clear space (each curve fits its mechanism to R^2=1)
    ax.text(0.045, 0.95, "Mother size at Start\n" + mother_eq(pm), transform=ax.transAxes,
            color=MOTHER_EDGE, fontsize=8.5, va="top", ha="left")
    ax.text(0.97, 0.06,
            "Daughter birth\n$V_d = r(g)\\,V_m,\\ r{:}\\,0.7{\\to}0.9$",
            transform=ax.transAxes, color=BUD_EDGE, fontsize=8.5, va="bottom", ha="right")
    ax.set_xlabel("Maternal replicative age (generations)")
    ax.set_ylabel("Volume (fL)")
    ax.set_title("(b) Maternal-age asymmetry, to the replicative lifespan", fontsize=11)
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.grid(axis="y", which="major", color="0.9", lw=0.7)
    ax.set_axisbelow(True)
    return r2m, r2d


def main():
    apply_style()
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.5, 4.6))
    fig.suptitle("Asymmetric cell growth across the replicative lifespan "
                 "(energetic cell-size model)", y=0.99, fontsize=12.5)
    panel_timecourse(axL)
    r2m, r2d = panel_maternal(axR)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = HERE / "asymmetric_growth_course.png"
    fig.savefig(out, bbox_inches="tight")
    print(f"wrote {out}  (mother R2={r2m:.5f}, daughter R2={r2d:.5f})")


if __name__ == "__main__":
    main()
