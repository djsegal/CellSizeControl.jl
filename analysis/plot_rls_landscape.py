#!/usr/bin/env python3
"""CC-1: the emergent-RLS parameter landscape (from gen_rls_landscape.jl). The mean and CV of
the replicative lifespan across the autocatalytic-damage parameters (viability threshold
D_crit x autocatalysis kappa). The Schnitzer-2022 calibration (mean ~25, CV ~0.3) sits in a
BROAD basin (contoured), not on a fine-tuned point — the emergent RLS is robust. Sequential
colorblind-safe colormap (cividis) for the scalar fields; contour marks the matching region.
Run via a venv with matplotlib + numpy.
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


def main():
    D, K, M, C = [], [], [], []
    with open(HERE / "rls_landscape.csv") as f:
        for row in csv.DictReader(f):
            D.append(float(row["D_crit"]))
            K.append(float(row["kappa"]))
            M.append(float(row["mean_rls"]))
            C.append(float(row["cv_rls"]))
    Dv, Kv = sorted(set(D)), sorted(set(K))
    nd, nk = len(Dv), len(Kv)
    di = {v: i for i, v in enumerate(Dv)}
    ki = {v: i for i, v in enumerate(Kv)}
    mean = np.full((nk, nd), np.nan)
    cv = np.full((nk, nd), np.nan)
    for d, k, m, c in zip(D, K, M, C):
        mean[ki[k], di[d]] = m
        cv[ki[k], di[d]] = c

    ext = [Dv[0], Dv[-1], Kv[0], Kv[-1]]

    apply_style()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.4))
    fig.suptitle("Emergent replicative-lifespan landscape (robust, not fine-tuned)",
                 y=0.99, fontsize=12)

    for ax, Z, label, title in (
        (axA, mean, "Mean RLS (divisions)", "(a) Mean lifespan"),
        (axB, cv, "Lifespan CV", "(b) Lifespan variability"),
    ):
        im = ax.imshow(Z, origin="lower", aspect="auto", extent=ext, cmap="cividis")
        fig.colorbar(im, ax=ax, label=label)
        ax.set(xlabel=r"Viability threshold $D_\mathrm{crit}$",
               ylabel=r"Autocatalysis $\kappa$", title=title)

    # contour the Schnitzer-matching region (mean 23-27 AND CV 0.25-0.35) on both panels
    match = ((mean >= 23) & (mean <= 27) & (cv >= 0.25) & (cv <= 0.35)).astype(float)
    for ax in (axA, axB):
        ax.contour(np.linspace(Dv[0], Dv[-1], match.shape[1]),
                   np.linspace(Kv[0], Kv[-1], match.shape[0]),
                   match, levels=[0.5], colors=VERM, linewidths=2.0)
    axA.plot([], [], color=VERM, lw=2.0, label="Schnitzer 2022 region")
    axA.legend(loc="upper right", frameon=True, framealpha=0.85, fontsize=8.5)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = HERE / "rls_landscape.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
