#!/usr/bin/env python3
"""Out-of-sample test of the single-asymmetry model: daughter replicative lifespan vs maternal
age (from daughter_rls_fraction.csv / daughter_rls_kennedy.csv).

The SAME age-eroding division asymmetry r(a) that enlarges the daughters of old mothers (the
size face) also -- under passive volume-proportional damage segregation -- loads those daughters
with a larger share of the mother's accumulated damage, so a daughter is born partway up the
autocatalytic damage trajectory and her own emergent lifespan is shortened. With NO parameter
beyond the McCormick-2015 wild-type damage calibration and the size-face asymmetry, the model
predicts the maternal-age daughter-lifespan deficit and is tested against Kennedy, Austriaco &
Guarente 1994 (J Cell Biol 127:1985), which was used nowhere in the fit.

Okabe-Ito palette (model = blue, independent data = vermillion). Run via a venv with
matplotlib + numpy.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _pubstyle import (apply_style, halo, opaque_legend, pub_audit, BLUE, VERM)

HERE = Path(__file__).resolve().parent


def _read_fraction():
    mid, mean, sd = [], [], []
    with open(HERE / "daughter_rls_fraction.csv") as f:
        for row in csv.DictReader(f):
            mid.append(float(row["frac_mid"]))
            mean.append(float(row["daughter_rls_mean"]))
            sd.append(float(row["daughter_rls_sd"]))
    return np.array(mid), np.array(mean), np.array(sd)


def _read_buckets():
    out = {}
    with open(HERE / "daughter_rls_kennedy.csv") as f:
        for row in csv.DictReader(f):
            out[row["bucket"]] = (float(row["frac_lo"]), float(row["frac_hi"]),
                                  float(row["model_rls"]), float(row["kennedy_rls"]))
    return out


def _read_posterior_summary():
    """Posterior-predictive median + 95% credible interval per quantity (from propagating the
    full McCormick-2015 ABC posterior through the Kennedy daughter-RLS prediction)."""
    out = {}
    with open(HERE / "daughter_rls_posterior_summary.csv") as f:
        for row in csv.DictReader(f):
            out[row["quantity"]] = (float(row["median"]), float(row["lo95"]),
                                    float(row["hi95"]), float(row["kennedy"]))
    return out


def main():
    mid, mean, sd = _read_fraction()
    bk = _read_buckets()
    ps = _read_posterior_summary()
    f70 = bk["first70"]  # (lo, hi, model, kennedy)
    l10 = bk["last10"]
    yng = ps["young_rls"]   # (median, lo95, hi95, kennedy)
    old = ps["old_rls"]
    fold_ci = ps["fold"]
    model_fold = fold_ci[0]
    ken_fold = f70[3] / l10[3]

    apply_style()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 4.9),
                                   gridspec_kw={"width_ratios": [1.55, 1.0]})
    fig.suptitle("Held-out prediction: daughters of old mothers are short-lived "
                 "(one age-eroding asymmetry)", y=0.995, fontsize=13.5)

    # ---- Panel (a): the full predicted curve vs Kennedy's two buckets -----------------------
    axA.fill_between(mid, mean - sd, mean + sd, color=BLUE, alpha=0.16, lw=0,
                     label="Model spread (±1 SD)")
    axA.plot(mid, mean, "-o", color=BLUE, lw=2.2, ms=5, mfc=BLUE, mec="white", mew=0.6,
             label="Model: inherited damage scaling\n(McCormick-2015 + size-face $r(a)$)")
    # Kennedy bucket bands (independent data): horizontal segments over their fraction ranges
    axA.hlines(f70[3], f70[0], f70[1], color=VERM, lw=3.0, zorder=5)
    axA.hlines(l10[3], l10[0], l10[1], color=VERM, lw=3.0, zorder=5,
               label="Kennedy 1994 bucket mean")
    axA.plot([(f70[0] + f70[1]) / 2, (l10[0] + l10[1]) / 2], [f70[3], l10[3]],
             "D", color=VERM, ms=9, mec="white", mew=0.8, zorder=6)
    # shade the two bucket ranges faintly so the reader sees the averaging windows
    axA.axvspan(f70[0], f70[1], color=VERM, alpha=0.05, lw=0)
    axA.axvspan(l10[0], l10[1], color=VERM, alpha=0.10, lw=0)
    halo(axA.text(0.35, 28.0, "first 70%", color=VERM, fontsize=12, ha="center", va="bottom"))
    halo(axA.text(0.95, 9.6, "last 10%", color=VERM, fontsize=12, ha="center", va="bottom"))
    axA.set(xlabel="Maternal replicative age (fraction of mother's lifespan)",
            ylabel="Daughter replicative lifespan (divisions)",
            title="(a) Daughter lifespan declines with maternal age")
    axA.set_xlim(0, 1.0)
    axA.set_ylim(0, 33)
    axA.grid(axis="y", which="major", color="0.9", lw=0.7)
    axA.set_axisbelow(True)
    opaque_legend(axA, loc="lower left", fontsize=12)

    # ---- Panel (b): the headline number -- the old-mother deficit FOLD, model vs data -------
    x = np.arange(2)
    w = 0.38
    model_vals = [yng[0], old[0]]          # posterior-predictive medians
    err_lo = [yng[0] - yng[1], old[0] - old[1]]
    err_hi = [yng[2] - yng[0], old[2] - old[0]]
    ken_vals = [f70[3], l10[3]]
    bM = axB.bar(x - w / 2, model_vals, w, color=BLUE, label="Model (posterior median)",
                 edgecolor="white", yerr=[err_lo, err_hi],
                 error_kw=dict(ecolor="0.25", elinewidth=1.3, capsize=4, capthick=1.3))
    bK = axB.bar(x + w / 2, ken_vals, w, color=VERM, label="Kennedy 1994", edgecolor="white")
    for i, r in enumerate(bM):     # label above the upper credible-interval cap
        halo(axB.text(r.get_x() + r.get_width() / 2, r.get_height() + err_hi[i] + 0.5,
                      f"{r.get_height():.1f}", ha="center", va="bottom", fontsize=12))
    for r in bK:
        halo(axB.text(r.get_x() + r.get_width() / 2, r.get_height() + 0.5,
                      f"{r.get_height():.1f}", ha="center", va="bottom", fontsize=12))
    axB.set_xticks(x)
    axB.set_xticklabels(["daughters of\nyoung mothers\n(first 70%)",
                         "daughters of\nold mothers\n(last 10%)"], fontsize=12)
    axB.set_ylabel("Daughter replicative lifespan (divisions)")
    axB.set_title("(b) The rejuvenation deficit", fontsize=13)
    axB.set_ylim(0, 40)
    axB.grid(axis="y", which="major", color="0.9", lw=0.7)
    axB.set_axisbelow(True)
    # Centered over the old-mother bar pair (data x=1.0, the divider between its blue/orange bars),
    # in the clear space above the short old-mother bars. Boxed for readability (not a halo hack).
    axB.text(1.0, 16.0,
             f"fold-drop\nmodel {model_fold:.1f}× ({fold_ci[1]:.1f}–{fold_ci[2]:.1f})"
             f"\nvs data {ken_fold:.1f}×",
             ha="center", va="center", fontsize=12, fontweight="bold", linespacing=1.4,
             bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="0.6", alpha=0.95))
    # the absolute-offset caveat (model runs ~1/3 low; the held-out claim is the fold, not the
    # level) lives in the caption, so the panel stays uncluttered.
    opaque_legend(axB, loc="upper right", bbox_to_anchor=(1.0, 0.99), fontsize=10.5)

    fig.tight_layout(rect=(0, 0, 1, 0.94))
    issues = pub_audit(fig)
    assert not issues, "daughter_rls pub_audit: " + "; ".join(issues)
    out = HERE / "daughter_rls.png"
    fig.savefig(out, bbox_inches="tight")
    print(f"wrote {out}")
    print(f"  model fold {model_fold:.2f}x vs Kennedy {ken_fold:.2f}x | "
          f"first70 {f70[2]:.1f} vs {f70[3]} | last10 {l10[2]:.1f} vs {l10[3]} | audit clean")


if __name__ == "__main__":
    main()
