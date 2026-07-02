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

from _pubstyle import (apply_style, pub_arrow, halo, opaque_legend, pub_audit,
                       BLUE, VERM, GREEN)

HERE = Path(__file__).resolve().parent
# Okabe-Ito: daughter / threshold / mother. The Start threshold here is the SAME calibrated
# setpoint as the bistable-switch figure: theta = c* = 0.449, with total dose W = 18, giving
# V* = W/theta = 40 fL. (Earlier this schematic used W = 60, theta = 1.5 for the same V* = 40,
# which read as inconsistent with the switch figure's c* ~ 0.45 even though both describe the
# identical setpoint -- only the Whi5 dose, and hence the concentration scale, differed.)
THRESH, VSTAR = 0.449, 40.0


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

    apply_style()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 4.8))
    fig.suptitle("Inhibitor (Whi5) dilution sets Start at a critical size", y=0.99,
                 fontsize=14)

    # Each cell's sizer step is its Whi5-dilution time; the fixed 19-min CLN2 timer
    # follows Start, so the total two-step G1 is the sizer step + 19 min.
    T_CLN2 = 19.0
    lbl = {
        "daughter": "Daughter: ~{step:.0f} min sizer step",
        "mother": "Mother: ~{step:.0f} min sizer step (born past $V^\\ast$)",
    }
    # Daughter: a genuine dilution trajectory (the blue curve) whose [Whi5] falls to the Start
    # threshold over its ~22-min "sizer step"; the marker is where she fires Start.
    axA.plot(t["daughter"], C["daughter"], "-", lw=2.4, color=col["daughter"],
             label=lbl["daughter"].format(step=t["daughter"][-1]))
    axA.plot([t["daughter"][-1]], [C["daughter"][-1]], "o", ms=7, color=col["daughter"], zorder=5)
    # Mother: born already past V*, so her [Whi5] is below threshold at birth -- she fires Start at
    # once and has NO dilution curve, only a point. Drawn (and shown in the legend) as a MARKER,
    # not a line, so the legend matches what is actually on the axes.
    axA.plot([t["mother"][-1]], [C["mother"][-1]], "o", ms=8.5, color=col["mother"], zorder=5,
             markeredgecolor="white", markeredgewidth=0.7,
             label="Mother: fires at once (born past $V^\\ast$)")
    axA.axhline(THRESH, color=VERM, lw=1.8, ls="--",
                label=r"Start threshold $\theta=c^\ast\approx0.45$")
    axA.set(xlabel="Time in G1 (min)", ylabel=r"Inhibitor concentration $[W]=W/V$",
            title="(a) Whi5 dilutes to the Start threshold")
    # headroom below the threshold so the annotations sit in the clear band under the curve
    cmin = min(min(C[c]) for c in ("daughter", "mother"))
    cmax = max(max(C[c]) for c in ("daughter", "mother"))
    axA.set_ylim(0.35, 0.65)  # headroom so the upper-right legend clears the blue daughter curve
    opaque_legend(axA, loc="upper right", fontsize=11)

    # Spell out the two-step G1 for each cell: sizer step + 19-min CLN2 timer. Each arrow is a
    # thick filled triangle (pub_arrow) and TARGETS ITS OWN COLOURED MARKER (the cell's Start
    # firing point), NOT the red threshold line. shrinkB keeps the triangle head off the marker
    # so it isn't buried under it.
    # No arrows: each curve is already colour-coded (blue daughter, green mother), so the
    # two-step-G1 note is just a colour-matched label near each curve (checklist: drop the arrow
    # when the series is already colour-coded).
    d_step, m_step = t["daughter"][-1], t["mother"][-1]
    ylo = cmin - 0.06
    halo(axA.text(d_step * 0.62 + 5.0, ylo + 0.018,
                  f"+ 19 min CLN2 timer\n~{d_step + T_CLN2:.0f} min total G1",
                  color=col["daughter"], fontsize=12, ha="center", va="bottom"))
    halo(axA.text(d_step * 0.22 - 0.33, ylo + 0.018,
                  f"+ 19 min CLN2 timer\n~{m_step + T_CLN2:.0f} min G1 (timer only)",
                  color=col["mother"], fontsize=12, ha="center", va="bottom"))

    # Same logic in volume: the daughter grows from her small birth size up to V* (blue curve);
    # the mother is born past V* and fires at once (green marker, no curve).
    axB.plot(t["daughter"], V["daughter"], "-", lw=2.4, color=col["daughter"])
    axB.plot([t["daughter"][-1]], [V["daughter"][-1]], "o", ms=7, color=col["daughter"], zorder=5)
    axB.plot([t["mother"][-1]], [V["mother"][-1]], "o", ms=8.5, color=col["mother"], zorder=5,
             markeredgecolor="white", markeredgewidth=0.7)
    axB.axhline(VSTAR, color=VERM, lw=1.8, ls="--", label=r"Critical size $V^\ast = W/\theta$")
    axB.set(xlabel="Time in G1 (min)", ylabel="Cell volume (fL)",
            title="(b) Growth to the critical size $V^\\ast$")
    opaque_legend(axB, loc="lower right", fontsize=11)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    issues = pub_audit(fig)
    assert not issues, "whi5_dilution pub_audit: " + "; ".join(issues)
    out = HERE / "whi5_dilution.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out, "| audit clean")


if __name__ == "__main__":
    main()
