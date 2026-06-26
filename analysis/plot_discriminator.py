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

from _pubstyle import apply_style, BLUE, VERM, GREEN, OKABE

HERE = Path(__file__).resolve().parent
COL = {"timer": BLUE, "adder": VERM, "sizer": GREEN}  # Okabe-Ito


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

    apply_style()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.4))

    for rule in ("timer", "adder", "sizer"):
        vb, vd = by_rule[rule]
        s = np.polyfit(vb, vd, 1)[0]
        axA.scatter(vb, vd, s=6, color=COL[rule], alpha=0.5,
                    label=f"{rule} (slope {s:.2f})")
    axA.set(xlabel=r"Birth volume $V_b$", ylabel=r"Division volume $V_d$",
            title="(a) The size-control slope discriminator")
    axA.legend(loc="upper left", frameon=False, fontsize=9)

    axB.plot(gen, tvb, "-", lw=2.0, color=VERM, label="Sub-doubling timer (collapses)")
    axB.plot(gen, svb, "-", lw=2.0, color=BLUE, label="Inhibitor-dilution sizer (stable)")
    axB.set(xlabel="Generation", ylabel="Birth volume",
            title="(b) The sizer stabilizes what the timer collapses")
    axB.legend(loc="center right", frameon=False, fontsize=9)

    fig.tight_layout()
    out = HERE / "discriminator.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
