#!/usr/bin/env python3
"""AGE-2 figure (Okabe-Ito): the replicative lifespan EMERGES from autocatalytic damage and a
viability threshold, rather than being a hard-coded generation cap.
  (a) the emergent RLS distribution (mean ~25, CV ~0.3 — calibrated to Schnitzer 2022), and
  (b) example damage trajectories accelerating to each cell's threshold (where it senesces).
Run via a venv with matplotlib + numpy.
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


def main():
    rls = []
    with open(HERE / "lifespan_samples.csv") as f:
        for row in csv.DictReader(f):
            rls.append(int(row["rls"]))
    rls = np.array(rls)
    m, sd = rls.mean(), rls.std(ddof=1)

    traces = defaultdict(lambda: ([], []))
    thr = {}
    with open(HERE / "damage_traces.csv") as f:
        for row in csv.DictReader(f):
            c = int(row["cell"])
            traces[c][0].append(float(row["age"]))
            traces[c][1].append(float(row["damage"]))
            thr[c] = float(row["threshold"])

    apply_style()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.4))
    fig.suptitle("The replicative lifespan emerges from autocatalytic damage", y=0.99,
                 fontsize=12)

    # (a) RLS distribution
    axA.hist(rls, bins=np.arange(0, rls.max() + 2) - 0.5, color=BLUE, alpha=0.85,
             edgecolor="white", linewidth=0.3)
    axA.axvline(m, color=VERM, lw=2.0, label=f"Mean {m:.1f} divisions")
    axA.axvspan(24, 26, color="0.6", alpha=0.18, label="Schnitzer 2022 (~24–26)")
    axA.set(xlabel="Replicative lifespan (divisions)", ylabel="Cells",
            title=f"(a) Emergent RLS distribution (CV = {sd / m:.2f})")
    axA.legend(loc="upper right", frameon=False, fontsize=9)
    axA.set_xlim(0, np.percentile(rls, 99.5) + 3)

    # (b) damage trajectories accelerating to each cell's threshold
    for i, (c, (ages, dmg)) in enumerate(sorted(traces.items())):
        axB.plot(ages, dmg, "-", lw=1.8, color=GREEN, alpha=0.85,
                 solid_capstyle="round")
        axB.plot([ages[-1]], [dmg[-1]], "o", ms=5, color=VERM, zorder=5)
        axB.hlines(thr[c], 0, ages[-1], color="0.75", lw=0.8, linestyles="dotted")
    axB.plot([], [], "-", color=GREEN, label="Accumulated damage $D(a)$")
    axB.plot([], [], "o", color=VERM, label="Senescence (threshold crossed)")
    axB.set(xlabel="Replicative age (divisions)", ylabel="Mother accumulated damage (a.u.)",
            title="(b) Autocatalytic damage reaches the viability threshold")
    axB.legend(loc="upper left", frameon=False, fontsize=9)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = HERE / "emergent_lifespan.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out, f"| mean={m:.1f} CV={sd / m:.2f}")


if __name__ == "__main__":
    main()
