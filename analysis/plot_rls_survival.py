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

from _pubstyle import apply_style, BLUE, VERM, GREEN, OKABE

HERE = Path(__file__).resolve().parent


def main():
    age, surv = [], []
    with open(HERE / "rls_survival.csv") as f:
        for row in csv.DictReader(f):
            age.append(int(row["age"]))
            surv.append(float(row["survival"]))
    age, surv = np.array(age), np.array(surv)
    # median survival age (where S crosses 0.5)
    med = int(np.argmin(np.abs(surv - 0.5)))

    apply_style()
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.step(age, surv, where="post", lw=2.4, color=BLUE)
    ax.fill_between(age, surv, step="post", alpha=0.12, color=BLUE)
    ax.axvline(med, color=VERM, lw=1.8, ls="--",
               label=f"Median lifespan {med} divisions")
    ax.axhline(0.5, color="0.7", lw=0.8)
    ax.set(xlabel="Replicative age (divisions)", ylabel="Surviving fraction $S(a)$",
           title="Emergent replicative-lifespan survival curve",
           xlim=(0, age.max()), ylim=(0, 1.02))
    ax.legend(loc="upper right", frameon=False, fontsize=12)

    fig.tight_layout()
    out = HERE / "rls_survival.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out, f"| median {med} divisions")


if __name__ == "__main__":
    main()
