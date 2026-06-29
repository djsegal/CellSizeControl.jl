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

from _pubstyle import (apply_style, pub_arrow, halo, opaque_legend, pub_audit,
                       BLUE, VERM, GREEN, REDPURPLE)

HERE = Path(__file__).resolve().parent
DAMAGE, CYCLE = REDPURPLE, GREEN  # Okabe-Ito reddish-purple + bluish-green


def main():
    gen, dmg, cyc = [], [], []
    with open(HERE / "cs_da_lineage.csv") as f:
        for row in csv.DictReader(f):
            gen.append(int(row["gen"]))
            dmg.append(float(row["Ddaughter"]))
            cyc.append(float(row["cycle"]))
    gen = np.array(gen, float)

    apply_style()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 4.8))
    fig.suptitle("The fitness face of replicative aging (same age-eroding asymmetry)",
                 y=0.99, fontsize=14)

    # Panel (a): inherited damage is illustrative (arbitrary units) -- there is no published
    # per-generation damage curve to target, so it is honestly labelled as such, no fake data.
    axA.plot(gen, dmg, "-", lw=2.0, color=DAMAGE, solid_capstyle="round")
    axA.set(xlabel="Maternal replicative age (generations)",
            ylabel="Daughter inherited damage (illustrative, a.u.)",
            title="(a) Daughters of old mothers inherit more damage")
    axA.set_xlim(0, max(gen) + 1)
    axA.grid(axis="y", which="major", color="0.9", lw=0.7)
    axA.set_axisbelow(True)

    # Panel (b): the cell-cycle-slows-with-age trend IS published. Reference targets (the
    # DIRECTION + magnitude, not a per-cell overlay -- model magnitudes are illustrative):
    #   Egilmez & Jazwinski 1989 (J Bacteriol 171:37): generation time rises ~5-6x by end of life;
    #   Fehrmann/Charvin 2013 (Cell Rep 5:1589): 78.3 min in young cells, then abrupt senescence entry;
    #   Moreno et al. 2019 (eLife 8:e48240): the lengthening is G1-specific (Whi5 ~3x in final cycles).
    cyc_arr = np.asarray(cyc, float)
    fold = cyc_arr.max() / cyc_arr.min()  # model end-of-life slowing fold (now ~5.2x)
    axB.plot(gen, cyc, "-", lw=2.0, color=CYCLE, solid_capstyle="round",
             label=f"Model cycle time ({fold:.1f}x over life)")
    # mark the published young-cell anchor (Fehrmann/Charvin 2013, 78.3 min) as a reference line
    axB.axhline(78.3, color="0.45", lw=1.2, ls="--",
                label="Fehrmann/Charvin 2013: 78.3 min (young)")
    axB.set(xlabel="Maternal replicative age (generations)", ylabel="Cycle time (min)",
            title="(b) The cell cycle slows with replicative age")
    axB.set_xlim(0, max(gen) + 1)
    axB.set_ylim(min(cyc) - 8, max(cyc) * 1.06)
    axB.grid(axis="y", which="major", color="0.9", lw=0.7)
    axB.set_axisbelow(True)
    # No arrow: the cycle curve is already green and named in the legend, so this is just a boxed
    # reference note set cleanly inside the 250-300 gridline band (checklist: drop the arrow when the
    # series is already colour-coded + legended).
    axB.text(max(gen) * 0.5, 135,
             "Cycle lengthens with age (Egilmez &\n"
             f"Jazwinski 1989, ~5-6×; G1-specific,\nMoreno 2019). Model: {fold:.1f}×",
             color="0.20", fontsize=10, ha="left", va="center", linespacing=1.35,
             bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="0.6", alpha=0.95))
    opaque_legend(axB, loc="upper left", fontsize=11)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    issues = pub_audit(fig)
    assert not issues, "fitness_face pub_audit: " + "; ".join(issues)
    out = HERE / "fitness_face.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out, f"| cycle {cyc_arr.min():.0f}->{cyc_arr.max():.0f} min ({fold:.2f}x) | audit clean")


if __name__ == "__main__":
    main()
