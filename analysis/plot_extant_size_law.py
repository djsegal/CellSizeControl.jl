#!/usr/bin/env python3
"""CC-X (from gen_extant_size_law.jl): the extant-vs-newborn size divergence + the senescence
age-law correction. (a) In a balanced exponentially-growing culture a snapshot samples every cell
at its last division: age-0 cells are the small buds (newborn distribution), age a>=1 cells are
mothers carrying their full retained body. The standing population over-represents the larger,
older mother bodies, so the MEAN EXTANT cell is D~1.97x the mean newborn -- a scale-free signature.
(b) The geometric age law 2^{-(a+1)} carries a senescence correction: at a short lifespan the
dividing population's age law is the truncated geometric lambda^{-(a+1)}, lambda<2 solving the
Euler-Lotka equation; lambda->2 (recovering 2^{-(a+1)}) as the lifespan grows, and lambda=phi at
rls=2. Okabe-Ito. Run via a venv with matplotlib + numpy.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _pubstyle import apply_style, BLUE, VERM, GREEN, ORANGE, halo, pub_audit

HERE = Path(__file__).resolve().parent
A0, VSTAR = 0.32, 60.0


def read_csv(name):
    with open(HERE / name) as f:
        return list(csv.DictReader(f))


def main():
    law = {r["source"]: r for r in read_csv("extant_size_law.csv") if r["case"] == "calibrated"}
    nb_mean = float(law["analytic"]["newborn_mean"])
    ext_mean = float(law["analytic"]["extant_mean"])
    D = float(law["analytic"]["divergence"])

    hist = read_csv("extant_size_hist.csv")
    ext = [(float(r["size_over_Vstar"]), float(r["fraction"]))
           for r in hist if r["population"] == "extant"]
    nb = [(float(r["size_over_Vstar"]), float(r["fraction"]))
          for r in hist if r["population"] == "newborn"]
    ext_x = np.array([p[0] for p in ext]) * VSTAR
    ext_y = np.array([p[1] for p in ext])
    nb_x = np.array([p[0] for p in nb]) * VSTAR
    nb_y = np.array([p[1] for p in nb])

    lam = read_csv("senescence_age_law.csv")
    lam_m = np.array([int(r["rls"]) for r in lam])
    lam_l = np.array([float(r["lambda"]) for r in lam])

    short = read_csv("senescence_short_rls.csv")
    ages = np.array([int(r["age"]) for r in short])
    p_lotka = np.array([float(r["lotka"]) for r in short])
    p_naive = np.array([float(r["naive_geometric"]) for r in short])
    p_sim = np.array([float(r["simulation"]) for r in short])
    m_short = len(ages)
    lam_short = float(lam[[int(r["rls"]) for r in lam].index(m_short)]["lambda"])

    apply_style()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.6, 4.6))
    fig.suptitle(
        "Size-structure signatures of balanced exponential growth", y=0.99, fontsize=12.5,
    )

    # (a) extant vs newborn size distributions
    wA = ext_x[1] - ext_x[0]
    axA.bar(nb_x, nb_y, width=wA * 0.95, color=VERM, alpha=0.55, edgecolor="white",
            linewidth=0.2, label="Newborns (age 0)", zorder=2)
    axA.bar(ext_x, ext_y, width=wA * 0.95, color=BLUE, alpha=0.40, edgecolor="white",
            linewidth=0.2, label="All extant cells", zorder=1)
    axA.set_ylim(0, 0.58)
    top = 0.58
    axA.axvline(nb_mean, color=VERM, lw=2.0, zorder=4)
    axA.axvline(ext_mean, color=BLUE, lw=2.0, zorder=4)
    halo(axA.text(nb_mean - 1.0, top * 0.40, f"newborn\n{nb_mean:.1f} fL",
                  color=VERM, fontsize=10.5, ha="right", va="top"))
    halo(axA.text(ext_mean + 1.2, top * 0.40, f"extant\n{ext_mean:.1f} fL",
                  color=BLUE, fontsize=10.5, ha="left", va="top"))
    halo(axA.annotate("", xy=(ext_mean, top * 0.22),
                      xytext=(nb_mean, top * 0.22),
                      arrowprops=dict(arrowstyle="<->", color="0.25", lw=1.6)))
    halo(axA.text(0.5 * (nb_mean + ext_mean), top * 0.26,
                  rf"$D=\overline{{V}}_{{\rm ext}}/\overline{{V}}_{{\rm nb}}={D:.2f}$",
                  ha="center", fontsize=11.5, color="0.15"))
    axA.set(xlabel="Cell volume (fL)", ylabel="Fraction of cells",
            title="(a) Extant cells outsize newborns (scale-free)", xlim=(0, 80))
    axA.legend(loc="upper right", frameon=True, fontsize=11)

    # (b) senescence correction to the geometric age law at a short lifespan
    axB.bar(ages, p_sim, width=0.62, color=GREEN, alpha=0.6, edgecolor="white",
            linewidth=0.2, label=f"Simulated dividers (rls={m_short})", zorder=1)
    axB.plot(ages, p_lotka, "o-", color=BLUE, lw=1.8, ms=6, zorder=3,
             label=r"Euler-Lotka $\lambda^{-(a+1)}$")
    axB.plot(ages, p_naive, "s--", color=VERM, lw=1.4, ms=5, zorder=2, alpha=0.9,
             label=r"Naive $2^{-(a+1)}$")
    axB.set_ylim(0, 0.58)
    halo(axB.text(0.97, 0.50,
                  rf"$\lambda={lam_short:.3f}<2$" "\n"
                  r"(solves $\lambda=\sum_{a=0}^{rls-1}\lambda^{-a}$)",
                  transform=axB.transAxes, ha="right", va="top", fontsize=11, color="0.15"))
    axB.set(xlabel="Replicative age $a$", ylabel="Fraction of dividing cells",
            title="(b) Senescence flattens the age law at short RLS")
    axB.set_xticks(ages)
    axB.legend(loc="upper right", frameon=True, fontsize=10.5)

    # inset: lambda(rls) -> 2, golden ratio at rls=2 (upper-centre, clear of the decaying age law)
    axins = axB.inset_axes((0.30, 0.52, 0.30, 0.34))
    axins.plot(lam_m, lam_l, "-", color="0.25", lw=1.6)
    axins.axhline(2.0, color=VERM, lw=1.0, ls=":")
    axins.plot([2], [lam_l[0]], "*", color=ORANGE, ms=11, zorder=5)
    axins.annotate(r"$\varphi$", (2, lam_l[0]), textcoords="offset points",
                   xytext=(7, -2), fontsize=11, color=ORANGE)
    axins.set(xlabel="mean RLS", ylabel=r"$\lambda$", xlim=(1, 40), ylim=(1.55, 2.05))
    axins.tick_params(labelsize=8)
    axins.xaxis.label.set_size(8.5)
    axins.yaxis.label.set_size(8.5)

    for viol in pub_audit(fig):
        print("AUDIT:", viol)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = HERE / "extant_size_law.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
