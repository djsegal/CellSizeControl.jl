#!/usr/bin/env python3
"""CC-2: the stochastic lineage ensemble (1e6 lineages, from gen_lineage_ensemble.jl).
(a) the daughter birth-size distribution, (b) the inherited-damage distribution, and
(c) the mother(at-division)→daughter(birth) size relationship. The tight coupling
(r≈0.98) is mechanistic, not homeostatic: the daughter is a fixed age-dependent fraction
r(a) of the enlarging mother, so the two rise together along the maternal-age trajectory.
Okabe-Ito. Run via a venv with matplotlib + numpy.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
BLUE, VERM, GREEN = "#0072b2", "#d55e00", "#009e73"  # Okabe-Ito


def main():
    summ = {}
    with open(HERE / "lineage_ensemble_summary.csv") as f:
        for row in csv.DictReader(f):
            summ[row["metric"]] = float(row["value"])
    vm, vd, dmg = [], [], []
    with open(HERE / "lineage_ensemble_pairs.csv") as f:
        for row in csv.DictReader(f):
            vm.append(float(row["Vmother_div"]))
            vd.append(float(row["Vdaughter"]))
            dmg.append(float(row["Ddaughter"]))
    vm, vd, dmg = np.array(vm), np.array(vd), np.array(dmg)
    r = summ["mother_daughter_size_corr"]
    N = int(summ["N_lineages"])

    plt.rcParams.update({"font.size": 10, "figure.dpi": 150, "savefig.dpi": 150,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(13.5, 4.2))
    fig.suptitle(f"Stochastic Lineage Ensemble (N = {N:,} Lineages, Division + Damage Noise)",
                 y=0.99, fontsize=12)

    axA.hist(vd, bins=50, color=BLUE, alpha=0.85, edgecolor="white", linewidth=0.3)
    axA.axvline(summ["daughter_size_mean"], color=VERM, lw=2.0,
                label=f"Mean {summ['daughter_size_mean']:.1f} fL  (CV {summ['daughter_size_cv']:.2f})")
    axA.set(xlabel="Daughter Birth Volume (fL)", ylabel="Daughters (Subsample)",
            title="(a) Daughter Birth-Size Distribution")
    axA.legend(loc="upper left", frameon=False, fontsize=9)

    axB.hist(dmg, bins=50, color=GREEN, alpha=0.85, edgecolor="white", linewidth=0.3)
    axB.axvline(summ["daughter_damage_mean"], color=VERM, lw=2.0,
                label=f"Mean {summ['daughter_damage_mean']:.2f}  (CV {summ['daughter_damage_cv']:.2f})")
    axB.set(xlabel="Inherited Damage (a.u.)", ylabel="Daughters (Subsample)",
            title="(b) Inherited-Damage Distribution")
    axB.legend(loc="upper right", frameon=False, fontsize=9)

    axC.scatter(vm, vd, s=5, color=BLUE, alpha=0.35, edgecolor="none")
    axC.set(xlabel="Mother Volume at Division (fL)", ylabel="Daughter Birth Volume (fL)",
            title=f"(c) Mother $\\to$ Daughter Size Coupling (r = {r:.2f})")

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = HERE / "lineage_ensemble.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
