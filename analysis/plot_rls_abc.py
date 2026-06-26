#!/usr/bin/env python3
"""Gallery (CC-4): ABC calibration of the emergent-RLS model to real wild-type budding-yeast
data (from gen_rls_abc.jl + gen_rls_abc_predictive.jl). Target = McCormick et al. 2015 (Cell
Metab 22:895-906, n=29,383 WT mother cells): mean RLS 26.6, SD 9.7. (a) the joint posterior over
the two damage parameters is a tight ridge (they trade off -- the RLS mean alone cannot separate
threshold from autocatalysis); (b) the pooled posterior-predictive RLS distribution matches the
measured mean and spread. Okabe-Ito. Run via a venv with matplotlib + scipy.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

BLUE, VERM, GREEN = "#0072b2", "#d55e00", "#009e73"  # Okabe-Ito
HERE = Path(__file__).resolve().parent
MEAN_T, SD_T = 26.6, 9.7  # McCormick 2015 WT pooled


def main() -> None:
    D, k = [], []
    with open(HERE / "rls_abc_posterior.csv") as fh:
        for r in csv.DictReader(fh):
            D.append(float(r["D_crit"]))
            k.append(float(r["kappa"]))
    D, k = np.array(D), np.array(k)
    rho = np.corrcoef(D, k)[0, 1]

    pp = np.array([int(r["rls"]) for r in csv.DictReader(open(HERE / "rls_abc_predictive.csv"))])

    plt.rcParams.update({"font.size": 11, "figure.dpi": 150, "savefig.dpi": 150,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.4))
    fig.suptitle("Calibrating the Emergent Replicative Lifespan to McCormick 2015 (ABC)",
                 y=0.99, fontsize=12)

    # (a) joint posterior (D_crit, kappa): the identifiability ridge
    hb = axA.hexbin(D, k, gridsize=34, cmap="Blues", mincnt=1)
    fig.colorbar(hb, ax=axA, pad=0.02, label="posterior density")
    axA.text(0.04, 0.93, rf"ridge corr $\rho={rho:.2f}$", transform=axA.transAxes,
             fontsize=10, color=VERM, fontweight="bold")
    axA.set(xlabel=r"Viability Threshold $D_{\rm crit}$",
            ylabel=r"Autocatalysis $\kappa$",
            title="(a) Joint Posterior: a Trade-Off Ridge")

    # (b) posterior-predictive vs the data
    bins = np.arange(0, pp.max() + 2) - 0.5
    axB.hist(pp, bins=bins, density=True, color=BLUE, alpha=0.55,
             label=f"model predictive\n(mean {pp.mean():.1f}, SD {pp.std():.1f})")
    xs = np.linspace(0, pp.max(), 300)
    axB.plot(xs, norm.pdf(xs, MEAN_T, SD_T), "-", lw=2.6, color=VERM,
             label=f"McCormick 2015 WT\n(mean {MEAN_T}, SD {SD_T})")
    axB.axvline(pp.mean(), color=BLUE, lw=1.4, ls="--")
    axB.axvline(MEAN_T, color=VERM, lw=1.4, ls=":")
    axB.set(xlabel="Replicative Lifespan (divisions)", ylabel="Probability Density",
            title="(b) Posterior-Predictive Matches the Data", xlim=(0, 60))
    axB.legend(loc="upper right", frameon=False, fontsize=8.5)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = HERE / "rls_abc.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
