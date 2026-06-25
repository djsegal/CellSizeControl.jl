#!/usr/bin/env python3
"""CS-DA figure: the maternal-age phenomenology from rigorous_cell_size.jl (the
literature-faithful cell-size-control model). From cs_da_lineage.csv:
  (a) daughter birth size RISES with maternal generation while the mother enlarges
      (asymmetry erosion; Kennedy 1994: old mothers -> larger daughters);
  (b) the cell cycle slows with replicative age (aged mothers divide progressively
      more slowly).
Run via rr_env. Writes cs_da_maternal_age.png."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent


def sat_fit(x, y):
    """Fit the saturating exponential V(a) = A + B·(1 − e^(−a/τ)) — the model's form for
    the maternal enlargement V*(a), the asymmetry erosion β(a), and the age-slowed cycle.
    Returns a dense (xs, ys) trend line plus the asymptote A+B and time constant τ."""
    import numpy as np
    from scipy.optimize import curve_fit

    def f(a, A, B, tau):
        return A + B * (1.0 - np.exp(-a / tau))

    x = np.asarray(x, float)
    y = np.asarray(y, float)
    p0 = [y[0], y[-1] - y[0], max(1.0, (x[-1] - x[0]) / 3.0)]
    try:
        popt, _ = curve_fit(f, x, y, p0=p0, maxfev=20000,
                            bounds=([-np.inf, -np.inf, 1e-3], [np.inf, np.inf, np.inf]))
    except Exception:
        popt = p0
    resid = y - f(x, *popt)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1.0 - np.sum(resid**2) / ss_tot if ss_tot > 0 else float("nan")
    xs = np.linspace(x.min(), x.max(), 200)
    return xs, f(xs, *popt), tuple(popt), r2


def _fit_caption(label, popt, r2):
    A, B, tau = popt
    return (f"{label}: $V\\!=\\!{A:.1f}\\!+\\!{B:.1f}(1\\!-\\!e^{{-a/{tau:.1f}}})$"
            f"  $R^2\\!=\\!{r2:.3f}$")


def main():
    gen, dau, mom, cyc = [], [], [], []
    with open(HERE / "cs_da_lineage.csv") as f:
        r = csv.DictReader(f)
        for row in r:
            gen.append(int(row["gen"]))
            dau.append(float(row["Vdaughter"]))
            mom.append(float(row["Vmother"]))
            cyc.append(float(row["cycle"]))

    plt.rcParams.update({"font.size": 10, "figure.dpi": 130,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 4.6))
    fig.suptitle("Maternal-age phenomenology from the energetic cell-size model "
                 "(Whi5 sizer + Di Talia G1 + asymmetry erosion)", y=0.99, fontsize=12)

    # saturating-exponential trend fits V(a)=A+B(1−e^(−a/τ)) (the model's form)
    xm, ym, pm, r2m = sat_fit(gen, mom)
    xd, yd, pd, r2d = sat_fit(gen, dau)
    axA.plot(gen, mom, "s", color="#3b3bbf", ms=5, label="mother size at Start")
    axA.plot(gen, dau, "o", color="#e08a1e", ms=5, label="daughter birth size")
    axA.plot(xm, ym, "-", color="#3b3bbf", lw=1.3, zorder=1)
    axA.plot(xd, yd, "-", color="#e08a1e", lw=1.3, zorder=1)
    axA.text(0.97, 0.03,
             _fit_caption("mother", pm, r2m) + "\n" + _fit_caption("daughter", pd, r2d),
             transform=axA.transAxes, va="bottom", ha="right", fontsize=7.5,
             bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))
    axA.set(xlabel="maternal replicative age (generation)", ylabel="volume (fL)",
            title="(a) mother enlarges; daughters grow with maternal age")
    axA.legend(loc="upper left", frameon=False, fontsize=9)
    axA.annotate("1st daughter\nsmall, pristine", xy=(gen[0], dau[0]),
                 xytext=(4.5, dau[0] + 4.5), fontsize=8, color="#b06010",
                 arrowprops=dict(arrowstyle="->", color="#b06010", lw=0.8))
    axA.annotate("late daughters\nlarge (asymmetry lost)", xy=(gen[-1], dau[-1]),
                 xytext=(gen[-1] - 10.5, dau[-1] - 11), fontsize=8, color="#b06010",
                 ha="center", arrowprops=dict(arrowstyle="->", color="#b06010", lw=0.8))

    xc, yc, pc, r2c = sat_fit(gen, cyc)
    axB.plot(gen, cyc, "^", color="#c44e52", ms=5, label="cycle time")
    axB.plot(xc, yc, "-", color="#c44e52", lw=1.3, zorder=1)
    axB.text(0.97, 0.03, _fit_caption("cycle", pc, r2c), transform=axB.transAxes,
             va="bottom", ha="right", fontsize=7.5,
             bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))
    axB.legend(loc="upper left", frameon=False, fontsize=9)
    axB.set(xlabel="maternal replicative age (generation)", ylabel="cycle time (min)",
            title="(b) the cell cycle slows with replicative age")

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = HERE / "cs_da_maternal_age.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
