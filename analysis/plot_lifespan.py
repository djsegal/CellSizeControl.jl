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

from _pubstyle import (apply_style, opaque_legend, halo, pub_audit,
                       BLUE, VERM, GREEN)

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
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 4.8))
    fig.suptitle("The replicative lifespan emerges from autocatalytic damage", y=0.99,
                 fontsize=14)

    # (a) RLS distribution. Published reference target: McCormick et al. 2015 (Cell Metab
    # 22:895) wild-type RLS, mean 26.6 with SD 9.7 (CV 0.365). These are not the main-text
    # headline numbers (the 780,000-daughter dissection) -- they are pooled/inferred from the
    # two WT controls in the supplemental data, and are shown here only as a summary target
    # band, not a per-cell overlay (the per-cell RLS distribution is not publicly deposited).
    MC_MEAN, MC_SD = 26.6, 9.7
    axA.hist(rls, bins=np.arange(0, rls.max() + 2) - 0.5, color=BLUE, alpha=0.85,
             edgecolor="white", linewidth=0.3)
    axA.axvspan(MC_MEAN - MC_SD, MC_MEAN + MC_SD, color="0.6", alpha=0.16,
                label=f"McCormick 2015 WT, suppl. data\n(mean {MC_MEAN}, SD {MC_SD}, CV "
                      f"{MC_SD / MC_MEAN:.2f})")
    axA.axvline(MC_MEAN, color="0.35", lw=1.6, ls=":")
    axA.axvline(m, color=VERM, lw=2.0,
                label=f"Model mean {m:.1f} (CV {sd / m:.2f})")
    axA.set(xlabel="Replicative lifespan (divisions)", ylabel="Cells",
            title="(a) Emergent RLS distribution vs published target")
    axA.set_xlim(0, np.percentile(rls, 99.5) + 3)
    axA.set_ylim(0, 400)  # headroom so the upper-right legend clears the histogram bars
    opaque_legend(axA, loc="upper right", fontsize=12)

    # (b) each cell's damage trajectory accelerates to ITS OWN viability threshold D_crit (a
    # per-cell lognormal draw -- the source of RLS spread). The senescence marker sits where
    # the green damage curve crosses that cell's dashed grey D_crit line. Each dashed line is
    # drawn at the cell's threshold and a short tick to the right of the crossing makes clear
    # the dot is ON its own threshold, not a single shared line.
    for i, (c, (ages, dmg)) in enumerate(sorted(traces.items())):
        ages = np.asarray(ages, float); dmg = np.asarray(dmg, float); T = thr[c]
        # Senescence is the FIRST age D >= D_crit; the raw last point overshoots the threshold,
        # which made the red dot float above its dashed line. Interpolate the exact crossing so the
        # green curve, the dashed threshold, and the senescence dot all meet at one point.
        k = int(np.argmax(dmg >= T)) if np.any(dmg >= T) else len(dmg) - 1
        if 0 < k < len(dmg) and dmg[k] >= T > dmg[k - 1]:
            f = (T - dmg[k - 1]) / (dmg[k] - dmg[k - 1])
            a_cross = ages[k - 1] + f * (ages[k] - ages[k - 1])
            ax_a = np.append(ages[:k], a_cross); ax_d = np.append(dmg[:k], T)
        else:
            a_cross = ages[-1]; ax_a, ax_d = ages, dmg
        axB.plot(ax_a, ax_d, "-", lw=1.8, color=GREEN, alpha=0.85, solid_capstyle="round")
        axB.hlines(T, 0, a_cross, color="0.55", lw=1.0, linestyles=(0, (4, 3)), zorder=2)
        axB.plot([a_cross], [T], "o", ms=6, color=VERM, zorder=5,
                 markeredgecolor="white", markeredgewidth=0.6)
    # legend proxies: curve, per-cell threshold line, senescence marker
    axB.plot([], [], "-", color=GREEN, lw=1.8, label=r"Accumulated damage $D(a)$")
    axB.plot([], [], ls=(0, (4, 3)), color="0.55", lw=1.0,
             label=r"Per-cell viability threshold $D_{\rm crit}$")
    axB.plot([], [], "o", color=VERM, markeredgecolor="white", markeredgewidth=0.6,
             label=r"Senescence ($D$ crosses $D_{\rm crit}$)")
    axB.set(xlabel="Replicative age (divisions)", ylabel="Mother accumulated damage (a.u.)",
            title="(b) Autocatalytic damage reaches the viability threshold")
    opaque_legend(axB, loc="upper left", bbox_to_anchor=(0.02, 0.92), fontsize=12)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    issues = pub_audit(fig)
    assert not issues, "emergent_lifespan pub_audit: " + "; ".join(issues)
    out = HERE / "emergent_lifespan.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out, f"| mean={m:.1f} CV={sd / m:.2f} | audit clean")


if __name__ == "__main__":
    main()
