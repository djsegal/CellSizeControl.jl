# Three-layer testing — CellSizeControl

Mirrors the science-space three-layer discipline (analytic / reference / cross-code),
in `test/runtests.jl`.

- **L1 — analytic limits.** The slope discriminator recovers each control law's
  closed-form slope: `TimerRule(2)` → 2, `AdderRule` → 1, `SizerRule`/
  `InhibitorDilutionSizer` → 0. The inhibitor-dilution setpoint equals `W/thresh`
  exactly and is birth-size-independent.
- **L2 — reference reproduction.** Reproduces the documented qualitative result
  (and the parent yeast-wcm `CellSize` finding): a **sub-doubling timer collapses**
  a symmetric-division lineage toward 0, while the **inhibitor-dilution sizer holds**
  births at `V*/2`. Grounded in Schmoller 2015 (the dilution sizer) and the
  Soifer–Amir 2016 phenomenology.
- **L3 — cross-source consistency.** The regression slope is monotone in timer
  strength (fold), a cross-check that the discriminator orders control regimes as
  the theory requires.

Next: add a **L2 quantitative** reproduction against the Chandler-Brown/Schmoller
2018 mechanistic Whi5:SBF titration ODEs (BioModels MODEL1803220001) and a
two-phase (G1-sizer + post-Start-timer) composition reproducing the whole-cycle
adder slope ≈ 1.
