#!/usr/bin/env python3
"""Cross-validation summary (from gen_crossval.jl): the model reproduces seven independent
published budding-yeast benchmarks from ONE parameterization, with no per-target refitting. Each
row shows the model's deviation from the published target in units of the literature tolerance
band (so different units -- slopes, minutes, divisions -- are comparable on one axis); points
inside the shaded band agree with the data. Okabe-Ito. Run via a venv with matplotlib.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _pubstyle import apply_style, BLUE, VERM, GREEN, OKABE

HERE = Path(__file__).resolve().parent


def main() -> None:
    rows = list(csv.DictReader(open(HERE / "crossval.csv")))
    rows.reverse()  # first metric on top
    labels, z, annot = [], [], []
    for r in rows:
        mod, ref = float(r["model"]), float(r["reference"])
        lo, hi = float(r["ref_lo"]), float(r["ref_hi"])
        half = (hi - lo) / 2
        z.append((mod - ref) / half)
        u = (" " + r["unit"]) if r["unit"] else ""
        labels.append(r["metric"])
        annot.append(f"model {mod:.2f} vs {ref:.2f}{u}  [{r['source']}]")

    apply_style()
    fig, ax = plt.subplots(figsize=(9.2, 4.6))
    y = range(len(labels))
    ax.axvspan(-1, 1, color=GREEN, alpha=0.13, label="Within published tolerance")
    ax.axvline(0, color="0.4", lw=1.4, ls="--", label="Published target")
    ax.scatter(z, y, s=90, color=BLUE, zorder=5)
    for yi, (zi, a) in enumerate(zip(z, annot)):
        ax.text(2.55, yi, a, va="center", fontsize=8.5, color="0.25")
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=12)
    ax.set_xlim(-2.5, 2.5)
    ax.set_xlabel("Model deviation from the published value (in units of the tolerance band)")
    ax.set_title("Cross-validation: one parameterization reproduces seven benchmarks",
                 fontsize=12)
    ax.legend(loc="lower left", frameon=False, fontsize=8.5)

    fig.tight_layout()
    out = HERE / "crossval.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
