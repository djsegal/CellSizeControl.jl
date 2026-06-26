#!/usr/bin/env python3
"""Size-control discriminator figure (from the package, via gen_discriminator.jl):
  (a) the Soifer-Amir Vd-vs-Vb slope recovers timer (~2), adder (~1), sizer (~0);
  (b) a sub-doubling timer collapses a lineage while the inhibitor-dilution sizer holds it.
Okabe-Ito palette. Run via a venv with matplotlib + numpy.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
COL = {"timer": "#0072b2", "adder": "#d55e00", "sizer": "#009e73"}  # Okabe-Ito


def main():
    by_rule = defaultdict(lambda: ([], []))
    with open(HERE / "discriminator.csv") as f:
        for row in csv.DictReader(f):
            by_rule[row["rule"]][0].append(float(row["Vb"]))
            by_rule[row["rule"]][1].append(float(row["Vd"]))
    gen, tvb, svb = [], [], []
    with open(HERE / "collapse.csv") as f:
        for row in csv.DictReader(f):
            gen.append(int(row["gen"]))
            tvb.append(float(row["timer_Vb"]))
            svb.append(float(row["sizer_Vb"]))

    plt.rcParams.update({"font.size": 10, "figure.dpi": 150, "savefig.dpi": 150,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.4))

    for rule in ("timer", "adder", "sizer"):
        vb, vd = by_rule[rule]
        s = np.polyfit(vb, vd, 1)[0]
        axA.scatter(vb, vd, s=6, color=COL[rule], alpha=0.5,
                    label=f"{rule} (slope {s:.2f})")
    axA.set(xlabel=r"Birth Volume $V_b$", ylabel=r"Division Volume $V_d$",
            title="(a) The Size-Control Slope Discriminator")
    axA.legend(loc="upper left", frameon=False, fontsize=9)

    axB.plot(gen, tvb, "-", lw=2.0, color="#d55e00", label="Sub-Doubling Timer (Collapses)")
    axB.plot(gen, svb, "-", lw=2.0, color="#0072b2", label="Inhibitor-Dilution Sizer (Stable)")
    axB.set(xlabel="Generation", ylabel="Birth Volume",
            title="(b) The Sizer Stabilizes What the Timer Collapses")
    axB.legend(loc="center right", frameon=False, fontsize=9)

    fig.tight_layout()
    out = HERE / "discriminator.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
