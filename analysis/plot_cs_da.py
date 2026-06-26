#!/usr/bin/env python3
"""Maternal-age phenomenology from the energetic cell-size model (rigorous_cell_size.jl ->
cs_da_lineage.csv), to the budding-yeast replicative lifespan:

  (a) the mother enlarges with replicative age and the daughters she produces grow with her
      age (Kennedy 1994), each fit by its mechanistic form;
  (b) the cell cycle slows with replicative age (accumulated damage lengthens the cycle).

Okabe-Ito palette (colorblind-safe): mother = blue, daughter = vermillion, cycle = green.
Run via a venv with matplotlib + scipy.
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
MOTHER, DAUGHTER, CYCLE = BLUE, VERM, GREEN  # Okabe-Ito


def _r2(y, yhat):
    y = np.asarray(y, float)
    return 1.0 - np.sum((y - yhat) ** 2) / np.sum((y - y.mean()) ** 2)


def single_exp(x, y):
    f = lambda g, a, b, c: a + b * np.exp(c * g)
    p, _ = curve_fit(f, x, y, p0=[y[-1], y[0] - y[-1], -0.1], maxfev=40000)
    xs = np.linspace(x.min(), x.max(), 300)
    return xs, f(xs, *p), p, _r2(y, f(x, *p))


def main():
    gen, dau, mom, cyc = [], [], [], []
    with open(HERE / "cs_da_lineage.csv") as f:
        for row in csv.DictReader(f):
            gen.append(int(row["gen"]))
            dau.append(float(row["Vdaughter"]))
            mom.append(float(row["Vmother"]))
            cyc.append(float(row["cycle"]))
    gen = np.array(gen, float)

    apply_style()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.6))
    fig.suptitle("Maternal-age phenomenology across the replicative lifespan "
                 "(energetic cell-size model)", y=0.99, fontsize=12.5)

    # (a) size: mother + daughter, curves only (fit-by-mechanism), Okabe-Ito
    xm, ym, _, _ = single_exp(gen, mom)
    xd, yd, _, _ = single_exp(gen, dau)
    axA.plot(xm, ym, "-", color=MOTHER, lw=2.0, solid_capstyle="round")
    axA.plot(xd, yd, "-", color=DAUGHTER, lw=2.0, solid_capstyle="round")
    axA.text(0.045, 0.95, "Mother size at Start", transform=axA.transAxes,
             color=MOTHER, fontsize=9, va="top", ha="left")
    axA.text(0.97, 0.06, "Daughter birth", transform=axA.transAxes,
             color=DAUGHTER, fontsize=9, va="bottom", ha="right")
    axA.set(xlabel="Maternal replicative age (generations)", ylabel="Volume (fL)",
            title="(a) The mother enlarges; daughters grow with her age")
    axA.set_xlim(0, max(gen) + 1)
    axA.yaxis.set_major_locator(MultipleLocator(5))
    axA.grid(axis="y", which="major", color="0.9", lw=0.7)
    axA.set_axisbelow(True)

    # (b) cycle time slows with replicative age (accumulated damage lengthens it)
    axB.plot(gen, cyc, "-", color=CYCLE, lw=2.0, solid_capstyle="round")
    axB.set(xlabel="Maternal replicative age (generations)", ylabel="Cycle time (min)",
            title="(b) The cell cycle slows with replicative age")
    axB.set_xlim(0, max(gen) + 1)
    axB.grid(axis="y", which="major", color="0.9", lw=0.7)
    axB.set_axisbelow(True)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = HERE / "cs_da_maternal_age.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
