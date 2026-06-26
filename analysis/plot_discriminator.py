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
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 4.8))

    # Soifer & Amir 2016 (Curr Biol 26:356) reference slope targets for the three regimes.
    # Drawn as faint reference lines through the data cloud's centroid so the model slope
    # can be read against the published target (sizer 0 / adder 1 / timer 2).
    SOIFER = {"sizer": 0.0, "adder": 1.0, "timer": 2.0}
    for rule in ("timer", "adder", "sizer"):
        vb, vd = by_rule[rule]
        vb = np.asarray(vb, float)
        vd = np.asarray(vd, float)
        s = np.polyfit(vb, vd, 1)[0]
        axA.scatter(vb, vd, s=6, color=COL[rule], alpha=0.5,
                    label=f"{rule} (model slope {s:.2f})")
        # reference target line: slope = Soifer-Amir value, anchored at the cloud centroid
        ref = SOIFER[rule]
        xc, yc = vb.mean(), vd.mean()
        xr = np.array([vb.min(), vb.max()])
        axA.plot(xr, yc + ref * (xr - xc), ls=(0, (4, 3)), lw=1.3, color=COL[rule],
                 alpha=0.9, zorder=2)
    # one legend entry naming the reference, color-neutral
    axA.plot([], [], ls=(0, (4, 3)), lw=1.3, color="0.35",
             label="Soifer-Amir 2016 targets:\nsizer 0 / adder 1 / timer 2")
    axA.set(xlabel=r"Birth volume $V_b$", ylabel=r"Division volume $V_d$",
            title="(a) The size-control slope discriminator")
    axA.legend(loc="upper left", frameon=False, fontsize=10.5)

    axB.plot(gen, tvb, "-", lw=2.0, color=VERM, label="Sub-doubling timer (collapses)")
    axB.plot(gen, svb, "-", lw=2.0, color=BLUE, label="Inhibitor-dilution sizer (stable)")
    axB.set(xlabel="Generation", ylabel="Birth volume",
            title="(b) The sizer stabilizes what the timer collapses")
    axB.legend(loc="center right", frameon=False, fontsize=11)

    fig.tight_layout()
    out = HERE / "discriminator.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
