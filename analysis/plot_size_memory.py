#!/usr/bin/env python3
"""Lineage birth-size MEMORY set by the single return-map pole r = alpha*f
(from the package, via gen_size_memory.jl):
  (a) lag-k lineage autocorrelation rho_k = r^k -- geometric decay at one rate r. The
      sizer (r=0) is memoryless; the timer (r=0.8) carries the longest memory.
  (b) the single-lineage invariant cv = sqrt(CV(Vb)^2 (1 - rho1^2)) recovers the intrinsic
      per-division noise (0.06) across every control mode -- mode- and set-point-free.
  (c) nutrient-shift step response (set-point doubling): mean birth size relaxes
      geometrically at rate r; the sizer absorbs the shift in one division, the timer takes
      several generations (memory -1/ln r).
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

from _pubstyle import apply_style, opaque_legend, pub_audit, halo, GREEN, ORANGE, BLUE

HERE = Path(__file__).resolve().parent
CV = 0.06

# the three canonical control modes (match the discriminator + CV figures)
MCOL = {"sizer": GREEN, "adder": ORANGE, "timer": BLUE}
MMRK = {"sizer": "o", "adder": "s", "timer": "^"}
# panel (a): one colour+marker per return slope r along the sizer->timer axis
ACOL = {0.00: GREEN, 0.25: "#7f7f7f", 0.50: ORANGE, 0.75: "#8064b0", 0.80: BLUE}
AMRK = {0.00: "o", 0.25: "v", 0.50: "s", 0.75: "D", 0.80: "^"}


def main():
    # ---- (a) autocorrelation rho_k vs lag, per return slope r ----
    by_r = defaultdict(lambda: ([], [], []))   # r -> (lag, rho_measured, rho_pred)
    with open(HERE / "size_memory_autocorr.csv") as fh:
        for row in csv.DictReader(fh):
            r = round(float(row["r"]), 2)
            by_r[r][0].append(int(row["lag"]))
            by_r[r][1].append(float(row["rho_measured"]))
            by_r[r][2].append(float(row["rho_pred"]))

    # ---- (b) invariant: recovered cv per mode (one point per mode, averaged over beta/f) ----
    by_mode_inv = defaultdict(list)            # mode -> [cv_recovered...]
    r_of_mode = {}
    with open(HERE / "size_memory_invariant.csv") as fh:
        for row in csv.DictReader(fh):
            by_mode_inv[row["mode"]].append(float(row["cv_recovered"]))
            r_of_mode[row["mode"]] = float(row["r"])

    # ---- (c) nutrient-step relaxation, per mode ----
    by_mode_step = defaultdict(lambda: ([], [], []))   # mode -> (gen, vb_measured, vb_geom)
    mem_of_mode = {}
    with open(HERE / "size_memory_step.csv") as fh:
        for row in csv.DictReader(fh):
            m = row["mode"]
            by_mode_step[m][0].append(int(row["gen"]))
            by_mode_step[m][1].append(float(row["vb_measured"]))
            by_mode_step[m][2].append(float(row["vb_geom"]))
            mem_of_mode[m] = float(row["memory_gen"])

    apply_style()
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(15, 4.6))

    # (a) geometric autocorrelation decay
    for r in sorted(by_r):
        lag, rm, rp = (np.asarray(v, float) for v in by_r[r])
        order = np.argsort(lag)
        lag, rm, rp = lag[order], rm[order], rp[order]
        col = ACOL.get(round(r, 2), "0.4")
        mrk = AMRK.get(round(r, 2), "o")
        # theory line only where the geometric law is above the noise floor
        axA.plot(lag, np.clip(rp, 1e-4, None), ls=(0, (5, 3)), lw=1.8, color=col, alpha=0.9,
                 zorder=2, solid_capstyle="round")
        axA.scatter(lag, np.clip(np.abs(rm), 1e-4, None), s=30, color=col, marker=mrk, alpha=0.9,
                    edgecolors="white", linewidths=0.4, zorder=3, label=fr"$r={r:.2f}$")
    axA.set_yscale("log")
    axA.set(xlabel="Lineage lag $k$ (generations)",
            ylabel=r"Birth-size autocorrelation $\rho_k$",
            title=r"(a) $\rho_k = r^k$: geometric memory decay")
    axA.set_xlim(0.7, 6.3)
    axA.set_ylim(1e-3, 1.3)
    leg = opaque_legend(axA, loc="lower left", fontsize=9.5,
                        title=r"return slope $r=\alpha f$" + "\n(dashed = $r^k$)",
                        title_fontsize=9.5)
    leg._legend_box.align = "left"

    # (b) the single-lineage invariant recovers the intrinsic noise
    order = ["sizer", "adder", "timer"]
    xs = np.arange(len(order))
    for i, m in enumerate(order):
        vals = np.asarray(by_mode_inv[m], float)
        axB.scatter(np.full_like(vals, i), vals, s=46, color=MCOL[m], marker=MMRK[m], alpha=0.9,
                    edgecolors="white", linewidths=0.5, zorder=3,
                    label=fr"{m}  ($r={r_of_mode[m]:.2f}$)")
    axB.axhline(CV, ls="--", lw=1.6, color="0.35", zorder=1)
    halo(axB.text(1.0, CV + 0.0016, r"intrinsic $cv=0.06$", ha="center", va="bottom",
                  fontsize=10.5, color="0.2"))
    axB.set_xticks(xs)
    axB.set_xticklabels([m for m in order])
    axB.set(xlabel="Control mode",
            ylabel=r"Recovered $\sqrt{CV(V_b)^2\,(1-\rho_1^2)}$",
            title=r"(b) $CV^2(1-\rho_1^2)=cv^2$: one-lineage invariant")
    axB.set_xlim(-0.5, 2.5)
    axB.set_ylim(0.055, 0.065)
    opaque_legend(axB, loc="upper right", fontsize=9.5)

    # (c) nutrient-shift step response
    for m in order:
        g, vm, vg = (np.asarray(v, float) for v in by_mode_step[m])
        o = np.argsort(g)
        g, vm, vg = g[o], vm[o], vg[o]
        axC.plot(g, vg, ls=(0, (5, 3)), lw=1.8, color=MCOL[m], alpha=0.9, zorder=2,
                 solid_capstyle="round")
        axC.scatter(g, vm, s=24, color=MCOL[m], marker=MMRK[m], alpha=0.9, edgecolors="white",
                    linewidths=0.4, zorder=3,
                    label=fr"{m}  ($\tau={mem_of_mode[m]:.2f}$ gen)")
    axC.axvline(0, ls=(0, (1, 3)), lw=1.1, color="0.7", zorder=0)
    halo(axC.text(10.5, 44.0, "set-point $V^*$ doubles", ha="center", va="bottom",
                  fontsize=10.5, color="0.25"))
    axC.set(xlabel="Generations after set-point step",
            ylabel=r"Mean birth size $\langle V_b \rangle$",
            title=r"(c) Nutrient-shift relaxation at rate $r$")
    axC.set_xlim(-0.5, 20)
    axC.set_ylim(top=118)
    opaque_legend(axC, loc="upper left", fontsize=9.5,
                  title=r"memory $\tau=-1/\ln r$", title_fontsize=9.5)

    fig.tight_layout()
    issues = pub_audit(fig)
    assert not issues, "size_memory pub_audit: " + "; ".join(issues)
    out = HERE / "size_memory.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out, "| audit clean")


if __name__ == "__main__":
    main()
