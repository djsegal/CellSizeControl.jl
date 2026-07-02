#!/usr/bin/env python3
"""CC-P: steady-state structure of an exponentially growing population (~10^6 cells, from
gen_population.jl). (a) the replicative-age distribution vs the geometric law P(age=a)=2^{-(a+1)}
(a straight line on a log axis — half the cells are virgin daughters, each older class half the
previous; Hartwell & Unger 1977); (b) the newborn (age-0) birth-size distribution, right-skewed
because the rare, enlarged old mothers bud the largest daughters. Okabe-Ito. Run via a venv with
matplotlib + numpy.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _pubstyle import apply_style, BLUE, VERM, GREEN

HERE = Path(__file__).resolve().parent


def main():
    ages, counts, fracs, geo = [], [], [], []
    with open(HERE / "population_age_structure.csv") as f:
        for row in csv.DictReader(f):
            ages.append(int(row["age"]))
            counts.append(int(row["count"]))
            fracs.append(float(row["fraction"]))
            geo.append(float(row["geometric"]))
    ages = np.array(ages)
    fracs = np.array(fracs)
    geo = np.array(geo)

    centers, ncount = [], []
    with open(HERE / "population_newborn_size.csv") as f:
        for row in csv.DictReader(f):
            centers.append(float(row["bin_center"]))
            ncount.append(int(row["count"]))
    centers = np.array(centers)
    ncount = np.array(ncount)

    summ = {}
    with open(HERE / "population_summary.csv") as f:
        for row in csv.DictReader(f):
            summ[row["metric"]] = row["value"]
    N = int(summ["N_cells"])
    ngen = int(summ["generations"])
    nb_mean = float(summ["newborn_size_mean"])

    apply_style()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.0, 4.4))
    fig.suptitle(
        f"Steady-state population structure (N = {N:,} cells, {ngen} generations)",
        y=0.99, fontsize=12,
    )

    # (a) replicative-age structure vs the geometric law, on a log axis
    keep = fracs > 0
    axA.semilogy(ages[keep], fracs[keep], "o", color=BLUE, markersize=8,
                 label="Simulated population")
    axA.semilogy(ages, geo, "-", color=VERM, lw=2.2,
                 label=r"Geometric law $2^{-(a+1)}$")
    axA.set(xlabel="Replicative age $a$ (buds produced)",
            ylabel="Fraction of population",
            title="(a) Replicative-age structure")
    axA.set_ylim(max(fracs[keep].min() / 2, 1e-7), 1.0)
    axA.legend(loc="upper right", frameon=False, fontsize=12)

    # (b) newborn birth-size distribution
    w = centers[1] - centers[0]
    axB.bar(centers, ncount / ncount.sum(), width=w * 0.95, color=GREEN,
            alpha=0.85, edgecolor="white", linewidth=0.2)
    axB.axvline(nb_mean, color=VERM, lw=2.0,
                label=f"Mean {nb_mean:.1f} fL")
    axB.axvline(0.32 * 60.0, color=BLUE, lw=2.0, ls="--",
                label="Young-mother floor 19.2 fL")
    axB.set(xlabel="Newborn birth volume (fL)", ylabel="Fraction of daughters",
            title="(b) Newborn-size distribution")
    axB.legend(loc="upper right", frameon=False, fontsize=12)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = HERE / "population_structure.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
