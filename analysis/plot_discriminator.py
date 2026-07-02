#!/usr/bin/env python3
"""Size-control discriminator figure (from the package, via gen_discriminator.jl):
  (a) the Soifer-Amir Vd-vs-Vb slope recovers timer (~2), adder (~1), sizer (~0);
  (b) a sub-doubling timer collapses a lineage while the inhibitor-dilution sizer holds it.
Okabe-Ito palette, redundant colour+marker encoding, opaque big-marker legends.
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
# Okabe-Ito, paired with a glyph so each regime reads in grayscale / for CVD.
COL = {"timer": BLUE, "adder": VERM, "sizer": GREEN}
MRK = {"timer": "^", "adder": "s", "sizer": "o"}


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
    # Drawn as a reference line through the data cloud's centroid so the model slope can be
    # read against the published target (sizer 0 / adder 1 / timer 2).
    SOIFER = {"sizer": 0.0, "adder": 1.0, "timer": 2.0}
    for rule in ("timer", "adder", "sizer"):
        vb, vd = by_rule[rule]
        vb = np.asarray(vb, float)
        vd = np.asarray(vd, float)
        s = np.polyfit(vb, vd, 1)[0]
        axA.scatter(vb, vd, s=20, color=COL[rule], marker=MRK[rule], alpha=0.45,
                    edgecolors="none", zorder=2,
                    label=f"{rule}  (model slope {s:.2f}, target {SOIFER[rule]:.0f})")
        # reference target line: slope = Soifer-Amir value, anchored at the cloud centroid
        ref = SOIFER[rule]
        xc, yc = vb.mean(), vd.mean()
        xr = np.array([vb.min(), vb.max()])
        axA.plot(xr, yc + ref * (xr - xc), ls=(0, (5, 3)), lw=2.0, color=COL[rule],
                 alpha=0.95, zorder=3, solid_capstyle="round")
    axA.set(xlabel=r"Birth volume $V_b$ (normalized)",
            ylabel=r"Division volume $V_d$ (normalized)",
            title="(a) The size-control slope discriminator")
    axA.set_ylim(bottom=0)  # honest baseline + opens room under the lower-right legend
    axA.set_xlim(right=1.4)  # cap the birth-volume axis
    leg = opaque_legend(axA, loc="lower right", markerscale=2.6, fontsize=12,
                        title="Soifer-Amir 2016 regimes\n(dashed = published target)",
                        title_fontsize=12)
    leg._legend_box.align = "left"

    # Dual y-axis: the sizer holds near 20 fL while the sub-doubling timer collapses toward 0; on a
    # shared axis the collapse is invisible (squished at the bottom), so each series gets its own scale.
    lS, = axB.plot(gen, svb, "-", lw=2.4, color=BLUE, marker="o", ms=5, markevery=4,
                   markeredgecolor="white", markeredgewidth=0.5,
                   label="Inhibitor-dilution sizer (stable)")
    axB.set_xlabel("Generation")
    axB.set_ylabel("Sizer birth volume (fL)", color=BLUE)
    axB.tick_params(axis="y", labelcolor=BLUE)
    axB.set_ylim(bottom=0)
    axB.set_title("(b) The sizer stabilizes what the timer collapses")
    axB.set_xlim(min(gen), max(gen))
    axBt = axB.twinx()
    axBt.spines["top"].set_visible(False)
    lT, = axBt.plot(gen, tvb, "-", lw=2.4, color=VERM, marker="s", ms=5, markevery=4,
                    markeredgecolor="white", markeredgewidth=0.5,
                    label="Sub-doubling timer (collapses)")
    axBt.set_ylabel("Timer birth volume (fL)", color=VERM)
    axBt.tick_params(axis="y", labelcolor=VERM)
    axBt.set_ylim(bottom=0)
    opaque_legend(axB, loc="center right", fontsize=11, handles=[lS, lT])

    fig.tight_layout()
    issues = pub_audit(fig)
    assert not issues, "discriminator pub_audit: " + "; ".join(issues)
    out = HERE / "discriminator.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out, "| audit clean")


if __name__ == "__main__":
    main()
