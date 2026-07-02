#!/usr/bin/env python3
"""Coupling the size-CV amplification to the damage/RLS recursion — the size-noise -> aging channel
(from the package, via gen_size_noise_aging.jl). Damage production scales with cell size, so a
timer's amplified birth-size CV (CV = cv_size / sqrt(1-(alpha f)^2)) passes into the per-division
damage-production noise. Panels:
  (a) the emergent RLS distribution (crit_cv=0, so mode enters ONLY through the size-noise channel):
      the timer's is BROADER than the sizer's while their means coincide (Delta ~ 0.06 divisions).
  (b) RLS CV vs the control slope alpha at fixed mean asymmetry f=0.40: broadening rises along
      sizer -> adder -> timer, tracking the size-CV amplification A(alpha,f).
  (c) the aging axis: as division symmetrizes with replicative age (f: 0.32 -> 0.50) the timer's
      RLS-CV/sizer ratio diverges (A_timer -> inf at f=0.5), the sizer stays flat.
Okabe-Ito palette, redundant colour+marker, opaque legends.
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

from _pubstyle import apply_style, opaque_legend, pub_audit, GREEN, ORANGE, BLUE

HERE = Path(__file__).resolve().parent

MCOL = {"sizer": GREEN, "adder": ORANGE, "timer": BLUE}
MMRK = {"sizer": "o", "adder": "s", "timer": "^"}
MA = {"sizer": 0, "adder": 1, "timer": 2}


def main():
    # ---- (a) RLS pmfs (sizer vs timer, isolated channel) ----
    rls, sizer_f, timer_f = [], [], []
    with open(HERE / "size_noise_rls_hist.csv") as fh:
        for row in csv.DictReader(fh):
            rls.append(int(row["rls"]))
            sizer_f.append(float(row["sizer_frac"]))
            timer_f.append(float(row["timer_frac"]))
    rls = np.asarray(rls); sizer_f = np.asarray(sizer_f); timer_f = np.asarray(timer_f)

    # ---- (b) RLS CV vs alpha at f=0.40, isolated channel ----
    modeA = {}       # mode -> (alpha, rls_cv)
    with open(HERE / "size_noise_rls.csv") as fh:
        for row in csv.DictReader(fh):
            if float(row["crit_cv"]) != 0.0:
                continue
            modeA[row["mode"]] = (float(row["alpha"]), float(row["rls_cv"]))

    # ---- (c) aging axis: RLS-CV ratio vs f, per mode ----
    by_mode = defaultdict(lambda: ([], []))    # mode -> (f, ratio)
    with open(HERE / "size_noise_rls_aging.csv") as fh:
        for row in csv.DictReader(fh):
            by_mode[row["mode"]][0].append(float(row["f"]))
            by_mode[row["mode"]][1].append(float(row["rls_cv_ratio_vs_sizer"]))

    apply_style()
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(16.5, 4.7))

    # (a) distributions — step outlines + light fill, means marked
    for name in ("sizer", "timer"):
        f = sizer_f if name == "sizer" else timer_f
        m = float(np.sum(rls * f))
        axA.step(rls, f, where="mid", lw=2.0, color=MCOL[name], alpha=0.95, zorder=3,
                 label=fr"{name} ($\alpha={MA[name]}$), mean$=${m:.2f}")
        axA.fill_between(rls, f, step="mid", color=MCOL[name], alpha=0.18, zorder=1)
        axA.axvline(m, ls=(0, (1, 3)), lw=1.3, color=MCOL[name], alpha=0.8, zorder=2)
    axA.set(xlabel="Replicative lifespan (divisions)", ylabel="Fraction of cells",
            title="(a) Timer's amplified size noise broadens RLS\n(equal means; crit$_{cv}=0$)")
    axA.set_xlim(rls.min() - 0.5, rls.max() + 0.5)
    opaque_legend(axA, loc="upper right", fontsize=10.5)

    # (b) RLS CV vs control slope
    order = ("sizer", "adder", "timer")
    xs = [modeA[n][0] for n in order]
    ys = [modeA[n][1] for n in order]
    axB.plot(xs, ys, ls="-", lw=1.8, color="0.55", alpha=0.7, zorder=1)
    for n in order:
        a, cv = modeA[n]
        axB.scatter([a], [cv], s=90, color=MCOL[n], marker=MMRK[n], edgecolors="white",
                    linewidths=0.6, zorder=3, label=f"{n} " + fr"($\alpha={MA[n]}$)")
    axB.set(xlabel=r"Control slope $\alpha$  (sizer 0 $\to$ adder 1 $\to$ timer 2)",
            ylabel=r"RLS CV  (emergent, crit$_{cv}=0$)",
            title=r"(b) RLS spread rises with the size-CV amplification" + "\n" +
                  r"($f=0.40$; cv$_{\mathrm{damage}}=A(\alpha,f)\,$cv$_{\mathrm{size}}$)")
    axB.set_xlim(-0.15, 2.15)
    axB.set_ylim(bottom=0.0)
    opaque_legend(axB, loc="upper left", fontsize=10.5)

    # (c) aging axis — RLS-CV ratio vs f
    for name in order:
        fs, rr = (np.asarray(v, float) for v in by_mode[name])
        o = np.argsort(fs)
        axC.plot(fs[o], rr[o], ls="-", lw=2.0, color=MCOL[name], alpha=0.9, zorder=2,
                 solid_capstyle="round")
        axC.scatter(fs[o], rr[o], s=34, color=MCOL[name], marker=MMRK[name], alpha=0.95,
                    edgecolors="white", linewidths=0.4, zorder=3,
                    label=f"{name} " + fr"($\alpha={MA[name]}$)")
    axC.axhline(1.0, ls=":", lw=1.4, color="0.5", zorder=1)
    axC.set(xlabel=r"Division asymmetry $f$  (young 0.32 $\to$ aged 0.50)",
            ylabel="RLS CV  (relative to sizer)",
            title="(c) Aging symmetrizes division:\nthe timer's RLS broadening diverges")
    axC.set_xlim(0.31, 0.49)
    axC.set_ylim(0.9, 3.0)
    opaque_legend(axC, loc="upper left", fontsize=10.5)

    fig.tight_layout()
    issues = pub_audit(fig)
    assert not issues, "size_noise pub_audit: " + "; ".join(issues)
    out = HERE / "size_noise_aging.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out, "| audit clean")


if __name__ == "__main__":
    main()
