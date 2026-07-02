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

from matplotlib.colors import TwoSlopeNorm

from _pubstyle import (apply_style, opaque_legend, halo, pub_audit,
                       SEQ_CMAP, DIV_CMAP, BLUE, VERM, GREEN)

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
    # signed log-ratio (homeostatic ~0 vs runaway large +): a DIVERGING quantity, so use the
    # CB-safe diverging colormap (vik) centred at 0 via TwoSlopeNorm. No RdBu.
    # vik (CB-safe + print-safe, via cmcrameri) centred at 0. Clip to the homeostatic range so the
    # near-0 structure (90% of cells lie in [-0.5, 1.1]) is resolved; the runaway tail (top ~5%, up
    # to +41) saturates and is marked by the colorbar's extend arrow. A full -|max|..|max| range
    # would give the tiny negative side half the bar (the misleading "-40" look).
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-0.6, vmax=1.2)
    pcm = axA.pcolormesh(A, F, LR, cmap=DIV_CMAP, norm=norm, shading="auto")
    cb = fig.colorbar(pcm, ax=axA, pad=0.02, extend="max", ticks=[-0.5, 0.0, 0.5, 1.0])
    cb.set_label(r"$\log_{10}(V_{\rm end}/V_0)$  (runaway $\to$)", fontsize=11)
    # sizer/adder/timer bin edges as slope contours; labels haloed so the field can't cut them
    cs = axA.contour(A, F, SL, levels=[0.5, 1.5], colors="k", linewidths=1.1, linestyles="-")
    lbls = axA.clabel(cs, fmt={0.5: "sizer | adder", 1.5: "adder | timer"}, fontsize=12)
    halo(lbls)
    # analytic homeostasis boundary alpha*f = 1  ->  f = 1/alpha
    aa = np.linspace(1 / F.max(), A.max(), 200)
    axA.plot(aa, 1 / aa, "--", color=VERM, lw=2.2, label=r"Homeostasis bound $\alpha f=1$")
    axA.set(xlabel=r"Control strength $\alpha$ (sizer 0 / adder 1 / timer 2)",
            ylabel=r"Division asymmetry $f$ (daughter fraction)",
            title="(a) Homeostatic vs runaway lineages", xlim=(A.min(), A.max()),
            ylim=(F.min(), F.max()))
    opaque_legend(axA, loc="lower left", fontsize=12)

    # (b) misclassification over (alpha, cv)
    A2, CV, MC = load_grid("phase_misclass.csv", "misclass")
    # misclassification probability is a SEQUENTIAL quantity (0..1): use viridis. No magma.
    pcm2 = axB.pcolormesh(A2, CV, MC, cmap=SEQ_CMAP, vmin=0, vmax=1, shading="auto")
    cb2 = fig.colorbar(pcm2, ax=axB, pad=0.02)
    cb2.set_label("P(misclassified)", fontsize=11)
    for edge in (0.35, 0.65, 1.35, 1.65):
        axB.axvline(edge, color="w", lw=0.8, ls=":", alpha=0.6)
    axB.set(xlabel=r"Control strength $\alpha$",
            ylabel=r"Measurement noise $cv$",
            title="(b) Where the discriminator becomes unreliable", xlim=(A2.min(), A2.max()),
            ylim=(CV.min(), CV.max()))
    # white text on the dark (low-misclass) viridis region, haloed dark so it survives lighter cells
    axB.text(0.04, 0.30, "n=80 cells, 300 replicates\nhard: strong sizers + bin edges",
             transform=axB.transAxes, fontsize=12, color="0.15", va="top", linespacing=1.3,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="0.6", alpha=0.92))

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    issues = pub_audit(fig)
    assert not issues, "phase_diagram pub_audit: " + "; ".join(issues)
    out = HERE / "phase_diagram.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out, "| audit clean")


if __name__ == "__main__":
    main()
