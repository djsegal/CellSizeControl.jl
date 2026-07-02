#!/usr/bin/env python3
"""Gallery (CC-6): the model's central falsifiable prediction against real data. The bistable
Whi5:SBF switch predicts the commitment size V* is proportional to the total Whi5 dose W
(V* = W/c*; from gen_whi5_dosage_test.jl). We test this against the Schmoller-2015 WHI5-dosage
series (cell size vs WHI5 gene copy number and ploidy), digitized from Heldt 2018 Fig 3C and
normalized to haploid 1xWHI5 (whi5_dosage_data.csv).

(a) Cell size vs Whi5 dose, with the model proportionality line. Within a fixed ploidy the
    measured size rises with Whi5 dose, as predicted; the two ploidies do not collapse onto one
    line, because ploidy adds a Whi5-independent size effect (genomic SBF-site titration) that
    pure inhibitor dilution does not capture.
(b) Fold-change in size for a doubling of WHI5 dose at fixed ploidy: the proportional law
    predicts 2.0x; the data give a sub-proportional ~1.2-1.4x. Direction confirmed, strict
    proportionality only approximate. Okabe-Ito. Run via a venv with matplotlib.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _pubstyle import apply_style, halo, opaque_legend, pub_audit, BLUE, VERM, GREEN

GREY = "#999999"
HERE = Path(__file__).resolve().parent


def load_model():
    a, vn = [], []
    with open(HERE / "whi5_dosage_model.csv") as fh:
        for row in csv.DictReader(fh):
            a.append(float(row["whi5_amount"]))
            vn.append(float(row["Vstar_norm"]))
    return a, vn


def load_data():
    rows = []
    with open(HERE / "whi5_dosage_data.csv") as fh:
        for row in csv.DictReader((r for r in fh if not r.startswith("#"))):
            rows.append(row)
    return rows


def main() -> None:
    amt, vnorm = load_model()
    rows = load_data()
    hap = {int(r["whi5_amount"][0]): float(r["size_rel"]) for r in rows if r["ploidy"] == "1"}
    dip = {int(r["whi5_amount"][0]): float(r["size_rel"]) for r in rows if r["ploidy"] == "2"}

    # honest head-to-head numbers
    fold_model = 2.0                       # proportional law: doubling W doubles V*
    fold_hap = hap[2] / hap[1]             # 1.40
    fold_dip = dip[2] / dip[1]             # 1.18
    ploidy_effect = dip[1] / hap[1]        # 1.71 at fixed (1x) Whi5

    apply_style()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 4.8))
    fig.suptitle(r"The predicted $V^\ast\propto W$ against the Schmoller 2015 Whi5-dosage data",
                 y=0.99, fontsize=14)

    # (a) size vs Whi5 dose: model proportionality line + data by ploidy
    axA.plot(amt, vnorm, "-", lw=2.4, color=GREY, zorder=2,
             label=r"Model $V^\ast=W/c^\ast$ (proportional)")
    # within-ploidy connectors (thin) to read the measured slope
    axA.plot([1, 2], [hap[1], hap[2]], ":", lw=1.6, color=BLUE, zorder=3)
    axA.plot([1, 2], [dip[1], dip[2]], ":", lw=1.6, color=VERM, zorder=3)
    axA.plot([1, 2], [hap[1], hap[2]], "o", ms=9, color=BLUE, zorder=5,
             label="Data: haploid (Schmoller 2015)")
    axA.plot([1, 2], [dip[1], dip[2]], "s", ms=8.5, color=VERM, zorder=5,
             label="Data: diploid (Schmoller 2015)")
    halo(axA.text(2.52, 0.16,
                  "At fixed Whi5 dose, diploids are larger:\n"
                  "ploidy adds a Whi5-independent effect\n"
                  "that pure dilution does not capture",
                  fontsize=12, color="0.35", ha="right", va="bottom"))
    axA.set(xlabel=r"Whi5 dose $W$ (relative to $1\times$ WHI5)",
            ylabel=r"Cell size $V^\ast$ (relative to $1\times$ WHI5)",
            title="(a) Critical size rises with Whi5 dose", xlim=(0, 2.6), ylim=(0, 3.0))
    opaque_legend(axA, loc="upper left", fontsize=12, markerscale=1.0, labelspacing=0.7)

    # (b) fold-change for a 2x Whi5 dose: predicted vs observed
    labels = ["Model\n(proportional)", "Data\nhaploid", "Data\ndiploid"]
    folds = [fold_model, fold_hap, fold_dip]
    colors = [GREY, BLUE, VERM]
    xpos = [0, 1, 2]
    axB.bar(xpos, folds, width=0.62, color=colors, zorder=2,
            edgecolor="white", linewidth=0.6)
    # reference lines drawn ABOVE the bars (z-order) so the gray bar's white border can't chop the
    # y=2.0 dashed line to half-thickness; labels boxed (not halo) with a clean unicode x
    axB.axhline(2.0, color="0.45", lw=1.3, ls="--", zorder=3)
    axB.axhline(1.0, color="0.1", lw=1.3, ls=":", zorder=3)  # near-black so it reads over the orange + gray bars
    _rbox = dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="0.7", alpha=0.95)
    axB.text(2.46, 2.08, "proportional (2×)", fontsize=12, color="0.25",
             ha="right", va="bottom", zorder=6, bbox=_rbox)
    axB.text(2.46, 0.92, "no change (1×)", fontsize=12, color="0.25",
             ha="right", va="top", zorder=6, bbox=_rbox)
    for x, f in zip(xpos, folds):
        halo(axB.text(x, f + 0.05, f"{f:.2f}$\\times$", ha="center", va="bottom",
                      fontsize=11.5, color="0.15"))
    # overlay the mixed dilution+titration extension (V* = g*V0 + W/c*, one extra param rho):
    # it turns the pure model's 2x into the observed sub-proportional folds, and haploid > diploid.
    mixed = {}
    mpath = HERE / "whi5_dosage_mixed.csv"
    if mpath.exists():
        with open(mpath) as fh:
            for row in csv.DictReader(fh):
                mixed[int(row["ploidy"])] = float(row["fold"])
    if mixed:
        mx = [1, 2]; my = [mixed[1], mixed[2]]
        axB.plot(mx, my, "D", ms=11, mfc="none", mec=GREEN, mew=2.2, zorder=7,
                 label="Mixed dilution+titration")
        opaque_legend(axB, loc="upper right", fontsize=11, markerscale=1.0)
    axB.set_xticks(xpos)
    axB.set_xticklabels(labels, fontsize=12)
    axB.set(ylabel=r"Size fold-change for $1\times\!\to\!2\times$ WHI5",
            title="(b) Doubling Whi5 dose (fixed ploidy)", ylim=(0, 2.45))

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    issues = pub_audit(fig)
    assert not issues, "whi5_dosage_test pub_audit: " + "; ".join(issues)
    out = HERE / "whi5_dosage_test.png"
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out, "| audit clean")
    print(f"fold-change 2x Whi5: model {fold_model:.2f}, haploid {fold_hap:.2f}, "
          f"diploid {fold_dip:.2f}; ploidy effect at 1x Whi5 = {ploidy_effect:.2f}")


if __name__ == "__main__":
    main()
