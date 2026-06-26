#!/usr/bin/env python3
"""The fitness face of replicative aging (from cs_da_lineage.csv): the same age-eroding
division asymmetry that grows the daughters (the size face, shown elsewhere) also
  (a) loads the daughter with inherited damage, and
  (b) slows the mother's cell cycle,
both rising across the replicative lifespan. Okabe-Ito palette. Run via a venv with
matplotlib + numpy.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
DAMAGE, CYCLE = "#cc79a7", "#009e73"  # Okabe-Ito reddish-purple + bluish-green


def main():
    gen, dmg, cyc = [], [], []
    with open(HERE / "cs_da_lineage.csv") as f:
        for row in csv.DictReader(f):
            gen.append(int(row["gen"]))
            dmg.append(float(row["Ddaughter"]))
            cyc.append(float(row["cycle"]))
    gen = np.array(gen, float)

    plt.rcParams.update({"font.size": 10, "figure.dpi": 150, "savefig.dpi": 150,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.4))
    fig.suptitle("The Fitness Face of Replicative Aging (Same Age-Eroding Asymmetry)",
                 y=0.99, fontsize=12)

    axA.plot(gen, dmg, "-", lw=2.0, color=DAMAGE, solid_capstyle="round")
    axA.set(xlabel="Maternal Replicative Age (Generation)",
            ylabel="Daughter Inherited Damage (a.u.)",
            title="(a) Daughters of Old Mothers Inherit More Damage")
    axA.set_xlim(0, max(gen) + 1)
    axA.grid(axis="y", which="major", color="0.9", lw=0.7)
    axA.set_axisbelow(True)

    axB.plot(gen, cyc, "-", lw=2.0, color=CYCLE, solid_capstyle="round")
    axB.set(xlabel="Maternal Replicative Age (Generation)", ylabel="Cycle Time (min)",
            title="(b) The Cell Cycle Slows With Replicative Age")
    axB.set_xlim(0, max(gen) + 1)
    axB.grid(axis="y", which="major", color="0.9", lw=0.7)
    axB.set_axisbelow(True)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = HERE / "fitness_face.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
