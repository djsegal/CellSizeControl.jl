#!/usr/bin/env python3
"""Gallery (CC-5): the inhibitor-dilution sizer from a bistable mechanism (from
gen_whi5_sbf_switch.jl). A Whi5:SBF double-negative feedback makes SBF activity bistable;
growth dilutes Whi5 (c = W/V), and the OFF/G1 state disappears at the saddle-node c*, firing
Start. The emergent set-point V* = W/c* is exactly linear in W — the phenomenological
inhibitor-dilution law V* = W/theta, with theta = c* now derived, not imposed. Okabe-Ito.
Run via a venv with matplotlib.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _pubstyle import apply_style, BLUE, VERM, GREEN, OKABE

GREY = "#999999"
HERE = Path(__file__).resolve().parent


def main() -> None:
    branches: dict[str, tuple[list[float], list[float]]] = {
        "off": ([], []),
        "unstable": ([], []),
        "on": ([], []),
    }
    with open(HERE / "whi5_sbf_bifurcation.csv") as fh:
        for row in csv.DictReader(fh):
            c, x = float(row["c"]), float(row["x"])
            branches[row["branch"]][0].append(c)
            branches[row["branch"]][1].append(x)
    cstar = min(branches["off"][0])  # lower fold = OFF saddle-node

    W, Vmech, Vlaw, thetas = [], [], [], []
    with open(HERE / "whi5_sbf_setpoint.csv") as fh:
        for row in csv.DictReader(fh):
            W.append(float(row["W"]))
            Vmech.append(float(row["Vstar_mech"]))
            Vlaw.append(float(row["Vstar_law"]))
            thetas.append(float(row["theta"]))
    theta = sum(thetas) / len(thetas)  # representative c* (constant to <1e-3 across W)

    apply_style()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 4.8))
    fig.suptitle("A bistable Whi5:SBF switch implements the inhibitor-dilution sizer", y=0.99,
                 fontsize=14)

    # (a) bifurcation / hysteresis
    axA.plot(branches["off"][0], branches["off"][1], "-", lw=2.6, color=BLUE,
             label="OFF / G1 (stable)")
    axA.plot(branches["on"][0], branches["on"][1], "-", lw=2.6, color=GREEN,
             label="ON / Start fired (stable)")
    axA.plot(branches["unstable"][0], branches["unstable"][1], "--", lw=1.8, color=GREY,
             label="Unstable threshold")
    axA.axvline(cstar, color=VERM, lw=1.8, ls=":")
    axA.text(cstar + 0.07, 0.72, r"Start at $c^\ast=W/V^\ast$", color=VERM, fontsize=11,
             ha="left", va="center")
    # point the dilution arrow AT the OFF-branch saddle-node (where the blue OFF curve ends and
    # turns into the grey unstable branch, i.e. where Start fires), not into empty space above it
    isn = branches["off"][0].index(cstar)
    x_sn = branches["off"][1][isn]
    axA.annotate("", xy=(cstar + 0.02, x_sn), xytext=(2.5, x_sn + 0.16),
                 arrowprops=dict(arrowstyle="-|>", color="0.35", lw=1.7, mutation_scale=14,
                                 shrinkA=2, shrinkB=4))
    axA.text(1.9, x_sn + 0.22, "Growth dilutes Whi5", fontsize=11, color="0.35", ha="center")
    axA.set(xlabel=r"Whi5 concentration $c=W/V$", ylabel="SBF activity (Start commitment)",
            title="(a) Whi5 dilution drives a bistable switch", xlim=(0, 3.0), ylim=(-0.03, 1.05))
    axA.legend(loc="lower right", frameon=False, fontsize=11)

    # (b) emergent sizer law V* = W/theta
    grid = [0] + W
    axB.plot(grid, [g / theta for g in grid], "-", lw=2.2, color=GREY,
             label=r"Law $V^\ast=W/\theta$  ($\theta=c^\ast$)")
    axB.plot(W, Vmech, "o", ms=8, color=BLUE, zorder=5, label="Mechanistic switch")
    axB.set(xlabel="Total Whi5 per cycle $W$", ylabel=r"Emergent set-point $V^\ast$ (fL)",
            title=r"(b) $V^\ast$ is exactly linear in $W$", xlim=(0, max(W) * 1.05),
            ylim=(0, max(Vmech) * 1.05))
    axB.legend(loc="upper left", frameon=False, fontsize=11)
    axB.text(0.97, 0.05, rf"$\theta=c^\ast={theta:.3f}$", transform=axB.transAxes,
             ha="right", fontsize=11, color="0.35")

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = HERE / "whi5_sbf_switch.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
