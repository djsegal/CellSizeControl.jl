#!/usr/bin/env python3
"""Gallery: the emergent replicative-lifespan survival curve (from gen_rls_survival.jl).
S(a) = fraction of mothers still dividing after a divisions — the classic budding-yeast
mortality curve — from the autocatalytic-damage + viability-threshold model. Okabe-Ito.
Run via a venv with matplotlib + numpy.
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
    age, surv = np.array(age), np.array(surv)
    # median survival age (where S crosses 0.5)
    med = int(np.argmin(np.abs(surv - 0.5)))

    plt.rcParams.update({"font.size": 11, "figure.dpi": 150, "savefig.dpi": 150,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.step(age, surv, where="post", lw=2.4, color=BLUE)
    ax.fill_between(age, surv, step="post", alpha=0.12, color=BLUE)
    ax.axvline(med, color=VERM, lw=1.8, ls="--",
               label=f"Median Lifespan {med} Divisions")
    ax.axhline(0.5, color="0.7", lw=0.8)
    ax.set(xlabel="Replicative Age (Divisions)", ylabel="Surviving Fraction $S(a)$",
           title="Emergent Replicative-Lifespan Survival Curve",
           xlim=(0, age.max()), ylim=(0, 1.02))
    ax.legend(loc="upper right", frameon=False, fontsize=10)

    fig.tight_layout()
    out = HERE / "rls_survival.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out, f"| median {med} divisions")


if __name__ == "__main__":
    main()
