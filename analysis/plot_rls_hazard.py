#!/usr/bin/env python3
"""Does the autocatalytic-damage mechanism produce GOMPERTZ mortality — the canonical aging
law (hazard rising exponentially with age)? From the emergent RLS survival data
(rls_survival.csv) we compute the discrete hazard h(a) = [S(a) - S(a+1)] / S(a) and fit the
Gompertz form ln h(a) = ln h0 + gamma * a. An exponentially-rising hazard (gamma > 0, good
log-linear fit) is the Gompertz signature; the slope gamma is the actuarial aging rate. We
compare to the budding-yeast literature (mean RLS ~24-30, max ~67; e.g. the 29.3-generation
microfluidic WT mean). Okabe-Ito. Run via a venv with matplotlib + numpy.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
BLUE, VERM = "#0072b2", "#d55e00"  # Okabe-Ito


def main():
    age, surv = [], []
    with open(HERE / "rls_survival.csv") as f:
        for row in csv.DictReader(f):
            age.append(int(row["age"]))
            surv.append(float(row["survival"]))
    age, surv = np.array(age, float), np.array(surv)

    # discrete hazard h(a) = (S(a) - S(a+1)) / S(a), over ages with enough cells
    a_h, h = [], []
    for i in range(len(age) - 1):
        if surv[i] > 0.02 and surv[i + 1] >= 0:   # drop the noisy far tail
            haz = (surv[i] - surv[i + 1]) / surv[i]
            if haz > 0:
                a_h.append(age[i])
                h.append(haz)
    a_h, h = np.array(a_h, float), np.array(h)

    # Gompertz fit: ln h = ln h0 + gamma * a   (mortality-rate doubling time = ln2/gamma)
    lo = a_h >= 8                                  # fit the aging phase (skip the flat young plateau)
    gamma, lnh0 = np.polyfit(a_h[lo], np.log(h[lo]), 1)
    fit_h = np.exp(lnh0 + gamma * a_h)
    # R^2 of the log-linear fit on the aging phase
    resid = np.log(h[lo]) - (lnh0 + gamma * a_h[lo])
    r2 = 1 - np.sum(resid**2) / np.sum((np.log(h[lo]) - np.mean(np.log(h[lo]))) ** 2)
    mrdt = np.log(2) / gamma                        # mortality-rate doubling time (divisions)

    plt.rcParams.update({"font.size": 11, "figure.dpi": 150, "savefig.dpi": 150,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    ax.semilogy(a_h, h, "o", ms=5, color=BLUE, alpha=0.8, label="Emergent Hazard $h(a)$")
    ax.semilogy(a_h, fit_h, "-", lw=2.2, color=VERM,
                label=f"Gompertz Fit ($\\gamma$ = {gamma:.3f}, $R^2$ = {r2:.2f})")
    ax.set(xlabel="Replicative Age (Divisions)",
           ylabel="Mortality Hazard $h(a)$ (log scale)",
           title="Gompertz-Like Mortality with Late-Life Deceleration")
    ax.annotate(f"Aging-phase doubling\ntime $\\approx$ {mrdt:.1f} divisions\n"
                "(+ young ramp, late plateau)",
                xy=(0.04, 0.96), xycoords="axes fraction", ha="left", va="top", fontsize=9,
                bbox=dict(boxstyle="round", fc="white", ec="0.8", alpha=0.9))
    ax.legend(loc="lower right", frameon=False, fontsize=9.5)

    fig.tight_layout()
    out = HERE / "rls_hazard.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out, f"| Gompertz gamma={gamma:.3f}, R2={r2:.2f}, MRDT={mrdt:.1f} div")


if __name__ == "__main__":
    main()
