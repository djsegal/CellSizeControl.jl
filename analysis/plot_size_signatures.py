#!/usr/bin/env python3
"""Gallery: four dynamical signatures of size control (from gen_size_signatures.jl), at matched
mean birth size. (a) birth-size distributions: a sizer is narrow, a timer broad; (b) the
consecutive birth-size inheritance map, whose slope is the size-control memory alpha*f; (c) the
step-response, relaxation of mean birth size after a 2x perturbation (a sizer forgets in one
generation); (d) the birth-size autocorrelation: a sizer is memoryless, a timer carries memory
across generations. Okabe-Ito. Run via a venv with matplotlib.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _pubstyle import apply_style, BLUE, VERM, GREEN, OKABE

HERE = Path(__file__).resolve().parent
COL = {"sizer": BLUE, "adder": GREEN, "timer": VERM}  # Okabe-Ito
ORDER = ["sizer", "adder", "timer"]
LBL = {"sizer": r"sizer ($\alpha$=0)", "adder": r"adder ($\alpha$=1)",
       "timer": r"timer ($\alpha$=1.5)"}


def load(fname):
    rows = list(csv.DictReader(open(HERE / fname)))
    return rows


def main() -> None:
    apply_style()
    plt.rcParams.update({"font.size": 10.5})
    fig, ax = plt.subplots(2, 2, figsize=(10.5, 8.2))
    fig.suptitle("Dynamical signatures of size control (matched mean birth size)",
                 y=0.995, fontsize=13)

    # (a) birth-size distributions
    H = defaultdict(list)
    for r in load("sig_hist.csv"):
        H[r["regime"]].append(float(r["Vb"]))
    bins = np.linspace(8, 36, 60)
    for reg in ORDER:
        v = np.array(H[reg])
        ax[0, 0].hist(v, bins=bins, density=True, histtype="step", lw=2.2, color=COL[reg],
                      label=f"{LBL[reg]}, CV={v.std() / v.mean():.2f}")
    ax[0, 0].set(xlabel="Birth volume $V_b$", ylabel="Probability density",
                 title="(a) Birth-size distributions")
    ax[0, 0].legend(frameon=False, fontsize=8.5)

    # (b) inheritance map Vb_{n+1} vs Vb_n
    P = defaultdict(lambda: ([], []))
    for r in load("sig_pairs.csv"):
        P[r["regime"]][0].append(float(r["Vb_n"]))
        P[r["regime"]][1].append(float(r["Vb_next"]))
    for reg in ORDER:
        x, y = np.array(P[reg][0]), np.array(P[reg][1])
        ax[0, 1].scatter(x, y, s=4, alpha=0.12, color=COL[reg], edgecolors="none")
        m = np.polyfit(x, y, 1)[0]
        xs = np.array([x.min(), x.max()])
        ax[0, 1].plot(xs, np.polyval(np.polyfit(x, y, 1), xs), "-", lw=2.2, color=COL[reg],
                      label=f"{LBL[reg]}, slope={m:.2f}")
    ax[0, 1].set(xlabel="Birth volume $V_b$ (gen $n$)",
                 ylabel="Birth volume $V_b$ (gen $n{+}1$)",
                 title=r"(b) Inheritance map (slope = memory $\alpha f$)")
    ax[0, 1].legend(frameon=False, fontsize=8.5, loc="upper left")

    # (c) step-response
    S = defaultdict(lambda: ([], []))
    for r in load("sig_step.csv"):
        S[r["regime"]][0].append(int(r["gen"]))
        S[r["regime"]][1].append(float(r["meanVb"]))
    for reg in ORDER:
        ax[1, 0].plot(S[reg][0], S[reg][1], "o-", ms=4, lw=2.0, color=COL[reg], label=LBL[reg])
    ax[1, 0].axhline(20.0, color="0.5", lw=1.0, ls="--", label="Set-point")
    ax[1, 0].set(xlabel="Generation after perturbation", ylabel="Mean birth volume $V_b$",
                 title="(c) Step-response (relaxation from $2\\times$)")
    ax[1, 0].legend(frameon=False, fontsize=8.5)

    # (d) autocorrelation
    A = defaultdict(lambda: ([], []))
    for r in load("sig_acf.csv"):
        A[r["regime"]][0].append(int(r["lag"]))
        A[r["regime"]][1].append(float(r["acf"]))
    for reg in ORDER:
        ax[1, 1].plot(A[reg][0], A[reg][1], "o-", ms=4, lw=2.0, color=COL[reg], label=LBL[reg])
    ax[1, 1].axhline(0.0, color="0.5", lw=0.8)
    ax[1, 1].set(xlabel="Generation lag", ylabel="Birth-size autocorrelation",
                 title="(d) Memory: a sizer forgets immediately")
    ax[1, 1].legend(frameon=False, fontsize=8.5)

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = HERE / "size_signatures.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
