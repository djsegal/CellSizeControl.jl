#!/usr/bin/env python3
"""CC-N: the population newborn-size law (from gen_newborn_size_law.jl). (a) the predicted
geometric mixture — mothers of replicative age a are a fraction 2^{-(a+1)} of the culture and
bud a daughter of size frac(a)*V*(a), so the age comb (stems, weight = marker area) lands on the
simulated newborn-size histogram; the distribution is right-skewed above the young-mother floor
0.32*V*. (b) the scale-free signature: rescaled by V*, the newborn distribution collapses onto one
curve for every set-point, so its CV/skew/ratio are pure numbers. Okabe-Ito. Run via a venv with
matplotlib + numpy.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _pubstyle import apply_style, BLUE, VERM, GREEN, halo, pub_audit

HERE = Path(__file__).resolve().parent
A0, VSTAR = 0.32, 60.0


def read_csv(name):
    with open(HERE / name) as f:
        return list(csv.DictReader(f))


def main():
    comb = read_csv("newborn_size_law_comb.csv")
    ages = np.array([int(r["age"]) for r in comb])
    weights = np.array([float(r["weight"]) for r in comb])
    csizes = np.array([float(r["newborn_size"]) for r in comb])

    hist = read_csv("newborn_size_hist.csv")
    xr = np.array([float(r["size_over_Vstar"]) for r in hist])
    hc = np.array([int(r["count"]) for r in hist], dtype=float)
    hc /= hc.sum()

    rows = read_csv("newborn_size_law.csv")
    law = {r["source"]: r for r in rows if r["case"] == "calibrated"}
    a_mean = float(law["analytic"]["mean"])
    a_cv = float(law["analytic"]["cv"])
    a_skew = float(law["analytic"]["skew"])
    a_ratio = float(law["analytic"]["ratio"])

    apply_style()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.4, 4.5))
    fig.suptitle(
        "Newborn-size law: the age-eroding asymmetry sampled by the geometric age structure",
        y=0.99, fontsize=12.5,
    )

    # (a) simulated newborn-size histogram (in fL) + the predicted geometric-mixture comb
    xr_fl = xr * VSTAR
    w_bar = xr_fl[1] - xr_fl[0]
    axA.bar(xr_fl, hc, width=w_bar * 0.95, color=GREEN, alpha=0.55,
            edgecolor="white", linewidth=0.2, label="Simulated newborns", zorder=1)
    # comb: a stem at each predicted size frac(a)*V*(a), marker area ~ geometric weight
    axA.vlines(csizes, 0, weights, color=BLUE, lw=1.4, alpha=0.9, zorder=3)
    axA.scatter(csizes, weights, s=20 + 900 * weights, color=BLUE, zorder=4,
                edgecolor="white", linewidth=0.6,
                label=r"Predicted comb $2^{-(a+1)}$")
    for a in (0, 1, 2, 3):
        axA.annotate(f"a={a}", (csizes[a], weights[a]), textcoords="offset points",
                     xytext=(6, 4), fontsize=10, color=BLUE)
    axA.axvline(A0 * VSTAR, color="0.45", lw=1.6, ls="--", zorder=2)
    halo(axA.text(A0 * VSTAR + 0.4, 0.44, r"floor $\alpha_0 V^*$", fontsize=10, color="0.3"))
    axA.axvline(a_mean, color=VERM, lw=2.0, zorder=2, label=f"Mean {a_mean:.1f} fL")
    axA.set(xlabel="Newborn birth volume (fL)", ylabel="Fraction of daughters",
            title="(a) Right-skewed geometric mixture", xlim=(16, 34))
    axA.legend(loc="upper right", frameon=True, fontsize=11)

    # (b) scale-free collapse: rescaled by V* the distribution is one curve; annotate the invariants
    axB.bar(xr, hc, width=(xr[1] - xr[0]) * 0.95, color=GREEN, alpha=0.7,
            edgecolor="white", linewidth=0.2)
    axB.axvline(a_ratio * A0, color=VERM, lw=2.0)
    txt = (f"scale-free (any $V^*$):\n"
           f"  ratio $= \\overline{{V_b}}/(\\alpha_0 V^*) = {a_ratio:.3f}$\n"
           f"  CV $= {a_cv:.3f}$\n"
           f"  skew $= {a_skew:.2f}$")
    halo(axB.text(0.97, 0.95, txt, transform=axB.transAxes, ha="right", va="top",
                  fontsize=11.5, color="0.15"))
    axB.set(xlabel=r"Newborn volume / $V^*$", ylabel="Fraction of daughters",
            title="(b) Scale-free signature")

    for viol in pub_audit(fig):
        print("AUDIT:", viol)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = HERE / "newborn_size_law.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
