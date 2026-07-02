#!/usr/bin/env python3
"""Size-control variability amplification and the aging homeostasis boundary
(from the package, via gen_size_cv_amplification.jl):
  (a) birth-size CV amplifies toward the homeostasis boundary along the sizer->adder->timer
      axis; measured points track CV = cv / sqrt(1 - (alpha*f)^2). The critical slope
      alpha_c = 1/f moves left as division symmetrizes (f: 0.32 -> 0.50), reaching the timer
      slope alpha = 2 at f = 0.5.
  (b) with age-eroding asymmetry (f: 0.32 -> 0.50), the CV amplification of the three canonical
      modes vs maternal replicative age: the sizer stays flat, the adder rises modestly, the
      timer runs to the boundary and diverges.
Okabe-Ito palette, redundant colour+marker encoding, opaque big-marker legends.
Run via a venv (or an env) with matplotlib + numpy.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _pubstyle import apply_style, opaque_legend, pub_audit, BLUE, VERM, GREEN, ORANGE, SKY

HERE = Path(__file__).resolve().parent
CV = 0.06

# panel (a): one colour+marker per division asymmetry f (young -> aged)
FCOL = {0.32: SKY, 0.40: ORANGE, 0.50: VERM}
FMRK = {0.32: "o", 0.40: "s", 0.50: "^"}
# panel (b): the three canonical control modes (match the discriminator figure)
MCOL = {"sizer": GREEN, "adder": ORANGE, "timer": BLUE}
MMRK = {"sizer": "o", "adder": "s", "timer": "^"}


def main():
    # ---- (a) CV vs control slope alpha, per asymmetry ----
    by_f = defaultdict(lambda: ([], [], []))   # f -> (alpha, cv_measured, cv_pred)
    with open(HERE / "size_cv_amplification.csv") as fh:
        for row in csv.DictReader(fh):
            f = float(row["f"])
            by_f[f][0].append(float(row["alpha"]))
            by_f[f][1].append(float(row["cv_measured"]))
            by_f[f][2].append(float(row["cv_pred"]))

    # ---- (b) CV amplification vs replicative age, per mode ----
    by_mode = defaultdict(lambda: ([], [], []))   # mode -> (age, amp_measured, amp_pred)
    with open(HERE / "size_cv_aging.csv") as fh:
        for row in csv.DictReader(fh):
            m = row["mode"]
            if row["amp_measured"] == "NaN":
                continue
            by_mode[m][0].append(int(row["age"]))
            by_mode[m][1].append(float(row["amp_measured"]))
            by_mode[m][2].append(float(row["amp_pred"]))

    apply_style()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 4.8))

    # (a)
    for f in (0.32, 0.40, 0.50):
        a, cm, cp = (np.asarray(v, float) for v in by_f[f])
        order = np.argsort(a)
        a, cm, cp = a[order], cm[order], cp[order]
        axA.plot(a, cp, ls=(0, (5, 3)), lw=2.0, color=FCOL[f], alpha=0.9, zorder=2,
                 solid_capstyle="round")
        axA.scatter(a, cm, s=26, color=FCOL[f], marker=FMRK[f], alpha=0.9, edgecolors="white",
                    linewidths=0.4, zorder=3,
                    label=fr"$f={f:.2f}$   ($\alpha_c=1/f={1/f:.2f}$)")
    axA.axhline(CV, ls=":", lw=1.4, color="0.5", zorder=1)
    for x, name in ((0.0, "sizer"), (1.0, "adder"), (2.0, "timer")):
        axA.axvline(x, ls=(0, (1, 3)), lw=1.1, color="0.7", zorder=0)
    axA.set(xlabel=r"Control slope $\alpha$  (sizer 0 $\to$ adder 1 $\to$ timer 2)",
            ylabel=r"Birth-size CV$(V_b)$",
            title="(a) Variability amplifies toward the homeostasis boundary")
    axA.set_xlim(-0.05, 2.0)
    axA.set_ylim(bottom=0.055)
    leg = opaque_legend(axA, loc="upper left", fontsize=9.5,
                        title=r"$CV = cv/\sqrt{1-(\alpha f)^2}$" + "\n(dashed = theory)",
                        title_fontsize=9.5)
    leg._legend_box.align = "left"

    # (b)
    for name in ("sizer", "adder", "timer"):
        ag, am, ap = (np.asarray(v, float) for v in by_mode[name])
        order = np.argsort(ag)
        ag, am, ap = ag[order], am[order], ap[order]
        axB.plot(ag, ap, ls=(0, (5, 3)), lw=2.0, color=MCOL[name], alpha=0.9, zorder=2,
                 solid_capstyle="round")
        axB.scatter(ag, am, s=26, color=MCOL[name], marker=MMRK[name], alpha=0.9,
                    edgecolors="white", linewidths=0.4, zorder=3,
                    label=f"{name}" + (r"  ($\alpha=%g$)" % {"sizer": 0, "adder": 1, "timer": 2}[name]))
    axB.axhline(1.0, ls=":", lw=1.4, color="0.5", zorder=1)
    axB.set(xlabel="Maternal replicative age (divisions)",
            ylabel=r"Birth-size CV amplification  (relative to sizer)",
            title="(b) Asymmetry erosion drives the timer to the boundary")
    axB.set_xlim(0, 60)
    axB.set_ylim(0.8, 8.5)
    opaque_legend(axB, loc="upper left", fontsize=10.5)

    fig.tight_layout()
    issues = pub_audit(fig)
    assert not issues, "size_cv pub_audit: " + "; ".join(issues)
    out = HERE / "size_cv_amplification.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out, "| audit clean")


if __name__ == "__main__":
    main()
