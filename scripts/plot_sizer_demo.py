#!/usr/bin/env python3
"""V4 — cell-size-control demo figure from docs/figdata/{lineages,scatter}.csv.
(left) a sub-doubling timer collapses a lineage while the inhibitor-dilution sizer
holds birth size steady; (right) the Soifer-Amir Vd-vs-Vb discriminator —
timer→2, adder→1, sizer→0. Run: <scratch>/rr_env/bin/python scripts/plot_sizer_demo.py
Writes docs/figures/cell_size_control_demo.png."""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FD = ROOT / "docs" / "figdata"
OUT = ROOT / "docs" / "figures" / "cell_size_control_demo.png"


def main():
    gen, tV, sV = [], [], []
    with open(FD / "lineages.csv") as f:
        for r in csv.DictReader(f):
            gen.append(int(r["gen"]))
            tV.append(float(r["timer_Vb"]))
            sV.append(float(r["sizer_Vb"]))

    pts = defaultdict(lambda: ([], []))
    slope = {}
    with open(FD / "scatter.csv") as f:
        for r in csv.DictReader(f):
            pts[r["regime"]][0].append(float(r["Vb"]))
            pts[r["regime"]][1].append(float(r["Vd"]))
            slope[r["regime"]] = float(r["slope"])

    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 130})
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Cell-size control: the sizer stabilizes what the timer collapses", y=1.0)

    axL.plot(gen, tV, "o-", color="#d62728", lw=2, ms=4,
             label="sub-doubling timer (collapses → 0)")
    axL.plot(gen, sV, "o-", color="#2ca02c", lw=2, ms=4,
             label="inhibitor-dilution sizer (stable at V*/2)")
    axL.axhline(20.0, color="#2ca02c", ls=":", lw=1, alpha=0.6)
    axL.set(xlabel="generation", ylabel="birth volume (a.u.)",
            title="(a) lineage birth size over generations")
    axL.legend(frameon=False, fontsize=9, loc="center right")

    colors = {"timer": "#1f77b4", "adder": "#ff7f0e", "sizer": "#2ca02c"}
    for name in ("timer", "adder", "sizer"):
        x, y = pts[name]
        axR.scatter(x, y, s=10, alpha=0.45, color=colors[name],
                    label=f"{name} (slope {slope[name]:.2f})")
    axR.set(xlabel="birth volume $V_b$", ylabel="division volume $V_d$",
            title="(b) Soifer–Amir discriminator: timer 2 / adder 1 / sizer 0")
    axR.legend(frameon=False, fontsize=9, loc="upper left")

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
