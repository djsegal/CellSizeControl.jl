#!/usr/bin/env python3
"""Gallery (CC-3): the size-control phase diagram (from gen_phase_diagram.jl). (a) the lineage
phase field over (control strength alpha, division asymmetry f): homeostatic vs runaway, with
the analytic boundary alpha*f = 1 and the Soifer-Amir sizer/adder/timer slope contours. (b)
where the discriminator misclassifies under measurement noise: the fraction of finite-sample
lineages whose recovered slope lands in the wrong bin, over (alpha, noise cv) — reliable except
near the bin edges, which blur as noise grows. Okabe-Ito accents. Run via a venv with matplotlib.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _pubstyle import apply_style, BLUE, VERM, GREEN, OKABE

HERE = Path(__file__).resolve().parent


def load_grid(fname, zkey):
    xs, ys, z = {}, {}, {}
    with open(HERE / fname) as fh:
        rows = list(csv.DictReader(fh))
    X = sorted({float(r[list(r)[0]]) for r in rows})
    Y = sorted({float(r[list(r)[1]]) for r in rows})
    xi = {v: i for i, v in enumerate(X)}
    yi = {v: i for i, v in enumerate(Y)}
    M = np.full((len(Y), len(X)), np.nan)
    k0, k1 = list(rows[0])[0], list(rows[0])[1]
    for r in rows:
        M[yi[float(r[k1])], xi[float(r[k0])]] = float(r[zkey])
    return np.array(X), np.array(Y), M


def main() -> None:
    apply_style()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.4, 4.9))
    fig.suptitle("Size-control phase diagram: homeostasis and discriminator reliability",
                 y=0.99, fontsize=14)

    # (a) phase field over (alpha, f): logratio = homeostatic (~0) vs runaway (large +)
    A, F, LR = load_grid("phase_alpha_f.csv", "logratio")
    _, _, SL = load_grid("phase_alpha_f.csv", "slope")
    pcm = axA.pcolormesh(A, F, LR, cmap="RdBu_r", vmin=-1.5, vmax=1.5, shading="auto")
    cb = fig.colorbar(pcm, ax=axA, pad=0.02)
    cb.set_label(r"$\log_{10}(V_{\rm end}/V_0)$  (runaway $\to$)", fontsize=11)
    # sizer/adder/timer bin edges as slope contours
    cs = axA.contour(A, F, SL, levels=[0.5, 1.5], colors="k", linewidths=1.0, linestyles="-")
    axA.clabel(cs, fmt={0.5: "sizer | adder", 1.5: "adder | timer"}, fontsize=10)
    # analytic homeostasis boundary alpha*f = 1  ->  f = 1/alpha
    aa = np.linspace(1 / F.max(), A.max(), 200)
    axA.plot(aa, 1 / aa, "--", color=VERM, lw=2.0, label=r"Homeostasis bound $\alpha f=1$")
    axA.set(xlabel=r"Control strength $\alpha$ (sizer 0 / adder 1 / timer 2)",
            ylabel=r"Division asymmetry $f$ (daughter fraction)",
            title="(a) Homeostatic vs runaway lineages", xlim=(A.min(), A.max()),
            ylim=(F.min(), F.max()))
    axA.legend(loc="lower left", frameon=True, framealpha=0.9, fontsize=10.5)

    # (b) misclassification over (alpha, cv)
    A2, CV, MC = load_grid("phase_misclass.csv", "misclass")
    pcm2 = axB.pcolormesh(A2, CV, MC, cmap="magma", vmin=0, vmax=1, shading="auto")
    cb2 = fig.colorbar(pcm2, ax=axB, pad=0.02)
    cb2.set_label("P(misclassified)", fontsize=11)
    for edge in (0.35, 0.65, 1.35, 1.65):
        axB.axvline(edge, color="w", lw=0.6, ls=":", alpha=0.5)
    axB.set(xlabel=r"Control strength $\alpha$",
            ylabel=r"Measurement noise $cv$",
            title="(b) Where the discriminator becomes unreliable", xlim=(A2.min(), A2.max()),
            ylim=(CV.min(), CV.max()))
    axB.text(0.04, 0.275, "n=80 cells, 300 replicates", fontsize=10, color="w")
    axB.text(0.04, 0.245, "hard: strong sizers + bin edges", fontsize=10, color="w")

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = HERE / "phase_diagram.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
