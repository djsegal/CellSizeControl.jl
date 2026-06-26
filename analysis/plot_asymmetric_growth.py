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

from _pubstyle import apply_style, BLUE, VERM, GREEN, OKABE

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
    ax.set_title("(a) Lineage growth over the replicative lifespan", fontsize=11)
    ax.set_xlim(0, t.max())
    ax.set_ylim(0, vt.max() * 1.14)
    ax.legend(loc="upper left", frameon=False, fontsize=9, handlelength=1.4)
    # the monotonic-mother invariant, stated where it reads
    ax.annotate("Mother body never shrinks\n(only the bud detaches)",
                xy=(t.max() * 0.34, vm[int(len(vm) * 0.34)] * 0.55),
                fontsize=8.0, color=BLUE, ha="center", va="center")


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
    hi = max(max(dau), max(mom)) + 4
    ax.set_ylim(lo, hi)
    ax.set_xlim(0, RLS + 1)
    # the markers are deterministic model output; the lines are the closed-form laws the
    # mechanism follows (a saturating exponential; a product of two saturating processes).
    # No R^2 is reported -- fitting a smooth curve to noiseless output is tautologically ~1.
    a, b, c = pm
    ax.text(0.045, 0.96,
            f"Mother: $V_m = {a:.0f} - {abs(b):.0f}\\,e^{{{c:.2g}\\,g}}$",
            transform=ax.transAxes, color=BLUE, fontsize=8.5, va="top", ha="left")
    ax.text(0.955, 0.07,
            "Daughter: $V_d = r(g)\\,V_m,\\ r{:}\\,0.7{\\to}0.9$",
            transform=ax.transAxes, color=VERM, fontsize=8.5, va="bottom", ha="right")
    # Published reference for the DIRECTION + magnitude of the trend: daughters of older
    # mothers are born larger (Johnston 1966, Antonie van Leeuwenhoek 32:94 -- directional).
    # Yang et al. 2011 (Cell Cycle 10:144, Table 1 / Fig 2A) quantify it: virgin daughters
    # ~6.9 um diameter -> terminal ~11.2 um, ~4x volume over the lifespan. Diameter-based and
    # not directly overlayable on this fL axis, so it is stated as a reference trend, with an
    # arrow marking the published direction (small early -> large late), not a fake curve.
    g_late = max(gen)
    d_early, d_late = dau[gen.index(min(gen))], dau[gen.index(g_late)]
    ax.annotate(
        "Daughters of old mothers are\nborn larger (Johnston 1966;\n"
        "Yang 2011, ~6.9 → 11 µm, ~4× vol)",
        xy=(g_late * 0.86, d_late * 0.965), xycoords="data",
        xytext=(g_late * 0.40, (d_early + d_late) * 0.5), textcoords="data",
        fontsize=7.6, color="0.30", ha="left", va="center",
        arrowprops=dict(arrowstyle="->", color="0.45", lw=1.0,
                        connectionstyle="arc3,rad=-0.15"))
    ax.set_xlabel("Maternal replicative age (generations)")
    ax.set_ylabel("Volume (fL)")
    ax.set_title("(b) Maternal-age asymmetry, to the replicative lifespan", fontsize=11)
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.grid(axis="y", which="major", color="0.9", lw=0.7)
    ax.set_axisbelow(True)
    leg = ax.legend(loc="lower right", frameon=True, fontsize=8.5, handletextpad=0.4,
                    bbox_to_anchor=(1.0, 0.16))
    leg.get_frame().set_facecolor("white")
    leg.get_frame().set_edgecolor("0.15")
    leg.get_frame().set_linewidth(0.8)
    leg.get_frame().set_alpha(1.0)
    return r2m, r2d


def main():
    apply_style()
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.0, 4.4))
    panel_timecourse(axL)
    r2m, r2d = panel_maternal(axR)
    fig.tight_layout(w_pad=2.0)
    out = HERE / "asymmetric_growth.png"
    fig.savefig(out, bbox_inches="tight")
    print(f"wrote {out}  (mother R2={r2m:.5f}, daughter R2={r2d:.5f})")


if __name__ == "__main__":
    main()
