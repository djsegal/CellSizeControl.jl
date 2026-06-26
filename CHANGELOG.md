# Changelog

All notable changes to the `CellSizeControl` package are documented here. This
project loosely follows [Keep a Changelog](https://keepachangelog.com/) and
[Semantic Versioning](https://semver.org/).

## [0.1.0] - 2026-06-26

### Added

- **Size-control rules and the slope discriminator.** The sizer / adder / timer
  rules as a birth-to-division volume map, the budding-yeast inhibitor-dilution
  (Whi5) sizer with set-point `V* = W/theta`, the saturating open-loop timer
  failure mode, and the model-agnostic Soifer--Amir slope discriminator
  (`size_control_slope`, `classify_control`) that recovers timer (2), adder (1),
  and sizer (0).
- **A mechanistic Whi5:SBF Start switch** (`Whi5SBFSwitch`). A bistable
  double-negative feedback whose saddle-node threshold reproduces the
  inhibitor-dilution set-point `V* = W/theta` from first principles, with the
  threshold emergent rather than imposed.
- **The continuous linear size-control map** (`LinearSizeControl`,
  `V_d = alpha V_b + beta`) unifying sizer / adder / timer in one parameter, with
  the homeostasis condition `alpha f < 1`.
- **A two-step G1 cell cycle** (`cell_cycle`, `lineage_timecourse`) reproducing
  the Di Talia mother and daughter G1 durations as an emergent asymmetry, on top
  of surface-area-limited and exponential single-cell growth.
- **A maternal-age aging layer.** One age-eroding division asymmetry
  (`aging_daughter_fraction`) drives both the size face (larger daughters) and the
  fitness face (inherited damage), with a biologically correct monotonic-mother
  division accounting (`simulate_aging_lineage`).
- **An emergent replicative lifespan** (`replicative_lifespan`,
  `lifespan_distribution`) from autocatalytic damage crossing a viability
  threshold, calibrated to the measured budding-yeast lifespan distribution.
- **Three-layer test suite** (analytic / published-reference / cross-source) plus
  Aqua, ExplicitImports, and JET quality gates; CI on Julia 1.10 and 1.11.
