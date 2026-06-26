#!/usr/bin/env python3
"""Gallery: the inhibitor-dilution sizer mechanism (from gen_whi5_dilution.jl). Whi5 is made
in a fixed dose W; [Whi5] = W/V dilutes as the cell grows, and Start fires when it crosses the
threshold θ (critical size V* = W/θ; Schmoller 2015). A daughter born small dilutes over a
long G1; a mother born ≥ V* is already past threshold and fires at once — the Di Talia
mother/daughter G1 asymmetry, mechanistically. Okabe-Ito. Run via a venv with matplotlib.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
BLUE, VERM, GREEN = "#0072b2", "#d55e00", "#009e73"  # Okabe-Ito daughter / threshold / mother
THRESH, VSTAR = 1.5, 40.0


def main():
    t = defaultdict(list)
    V = defaultdict(list)
    C = defaultdict(list)
    with open(HERE / "whi5_dilution.csv") as f:
        for row in csv.DictReader(f):
            c = row["cell"]
            t[c].append(float(row["t"]))
            V[c].append(float(row["V"]))
            C[c].append(float(row["whi5_conc"]))
    col = {"daughter": BLUE, "mother": GREEN}

    plt.rcParams.update({"font.size": 11, "figure.dpi": 150, "savefig.dpi": 150,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.4))
    fig.suptitle("Inhibitor (Whi5) Dilution Sets Start at a Critical Size", y=0.99,
                 fontsize=12)

    for c in ("daughter", "mother"):
        g1 = t[c][-1]
        axA.plot(t[c], C[c], "-", lw=2.4, color=col[c],
                 label=f"{c.capitalize()} (G1 sizer step {g1:.0f} min)")
        axA.plot([t[c][-1]], [C[c][-1]], "o", ms=7, color=col[c], zorder=5)
    axA.axhline(THRESH, color=VERM, lw=1.8, ls="--", label=r"Start Threshold $\theta$")
    axA.set(xlabel="Time in G1 (min)", ylabel=r"Inhibitor Concentration $[W]=W/V$",
            title="(a) Whi5 Dilutes to the Start Threshold")
    axA.legend(loc="upper right", frameon=False, fontsize=9)

    for c in ("daughter", "mother"):
        axB.plot(t[c], V[c], "-", lw=2.4, color=col[c])
        axB.plot([t[c][-1]], [V[c][-1]], "o", ms=7, color=col[c], zorder=5)
    axB.axhline(VSTAR, color=VERM, lw=1.8, ls="--", label=r"Critical Size $V^\ast = W/\theta$")
    axB.set(xlabel="Time in G1 (min)", ylabel="Cell Volume (fL)",
            title="(b) Growth to the Critical Size $V^\\ast$")
    axB.legend(loc="lower right", frameon=False, fontsize=9)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = HERE / "whi5_dilution.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
