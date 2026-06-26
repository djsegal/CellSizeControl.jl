# CellSizeControl.jl

Generic cell-size-control primitives in Julia: the **sizer / adder / timer** rules, the
model-agnostic **slope discriminator**, the budding-yeast **inhibitor-dilution sizer** (and a
mechanistic **bistable Whi5:SBF switch** that reproduces it from first principles), a two-step
G1 **cell cycle**, and a **maternal-age aging** layer in which one age-eroding division
asymmetry ties daughter size, inherited damage, and the **replicative lifespan** to a single
mechanism.

```julia
using CellSizeControl

# Which control law does a lineage obey?  (timer→2, adder→1, sizer→0)
s = simulate_lineage(AdderRule(1.0); n=600)
classify_control(size_control_slope(s.Vb, s.Vd))      # :adder

# Whi5-style inhibitor dilution: Start at V* = W/[Whi5]* — a mechanistic sizer
ids = InhibitorDilutionSizer(60.0, 1.5)               # V* = 40
setpoint_volume(ids)                                  # 40.0

# ...and the SAME sizer law emerges from a bistable Whi5:SBF Start switch (θ derived,
# not imposed): growth dilutes Whi5 until the OFF/G1 state vanishes at a saddle-node.
sw = Whi5SBFSwitch(18.0)
setpoint_volume(sw), whi5_sbf_threshold(sw)           # (≈ 40, ≈ 0.45)  →  V* = W/θ

# Two-step G1 (Di Talia 2007): a mother (born ≥ V*) runs ≈ the Cln2 timer; a daughter
# (born small) spends extra time growing to V* — the G1 asymmetry emerges, not imposed.
cell_cycle(40.0; Vstar=36.0).G1                       # ≈ 19 min (mother)

# The replicative lifespan EMERGES from autocatalytic damage + a viability threshold.
# At the default damage params it is ≈ 25 divisions; ABC-calibrating them to real data
# (McCormick 2015 WT) gives mean 26.6, CV 0.37 — see analysis/gen_rls_abc.jl.
using Statistics: mean
mean(lifespan_distribution(2000))                     # ≈ 25 divisions (defaults)
```

## What it is

The size-control question — *what division rule keeps a growing, dividing population's size
distribution stable?* — is generic and reusable, so it lives in its own package. It holds the
three control laws (`TimerRule`, `AdderRule`, `SizerRule`) and the continuous linear map
unifying them (`LinearSizeControl`, `Vd = αVb + β`); the mechanistic `InhibitorDilutionSizer`
and the bistable `Whi5SBFSwitch` that reproduces its `V* = W/θ` law from a double-negative
feedback; the Soifer–Amir 2016 `Vd`-vs-`Vb` slope classifier (timer 2 / adder 1 / sizer 0); an
energetic two-step-G1 `cell_cycle`; and a maternal-age layer (`simulate_aging_lineage`,
`replicative_lifespan`, `lifespan_distribution`) in which a single age-eroding asymmetry `r(a)`
is both the *size* face (larger daughters) and the *fitness* face (more inherited damage, finite
lifespan). The mother keeps her cell body at division (monotonic, never shrinks); only the bud
leaves.

## Science grounding

One parameterization reproduces seven independent published benchmarks (see the cross-validation
in `analysis/`):

- **Slope discriminator:** Soifer, Robert & Amir 2016, *Curr Biol* 26:356 (timer 2 / adder 1 /
  sizer 0). **Linear size-control map:** Amir 2014, *Phys Rev Lett* 112:208102.
- **Inhibitor-dilution sizer:** Schmoller, Turner, Kõivomägi & Skotheim 2015, *Nature* 526:268
  (V\* = W/[Whi5]\*). **Bistable Start switch:** Skotheim et al. 2008, *Nature* 454:291; Charvin
  et al. 2010, *PLoS Biol* 8:e1000284 — the `Whi5SBFSwitch` realizes the sizer law from this.
- **Two-step G1:** Di Talia, Skotheim, Bean, Siggia & Cross 2007, *Nature* 448:947 (mother ≈ 19
  min, daughter ≈ 45 min). **Budded timer:** Allard, Decker, Weiner, Toettcher & Graziano 2018,
  *PLoS One* 13:e0209301.
- **Maternal-age daughter size:** Kennedy, Austriaco & Guarente 1994, *J Cell Biol* 127:1985.
  **Replicative lifespan:** ABC-calibrated to McCormick et al. 2015, *Cell Metab* 22:895
  (wild-type mean ≈ 26.6, CV ≈ 0.37).
- Full source notes: the parent workshop's `docs/literature/AUDIT_NOTES.md` §4.

## Testing (three-layer)

`julia --project=. -e 'using Pkg; Pkg.test()'` — L1 analytic (the slope recovers each regime;
V\*=W/thresh; the bistable switch is genuinely hysteretic and gives V\*=W/c\*), L2 reference (a
sub-doubling timer collapses a lineage while the inhibitor-dilution sizer holds it; the Di Talia
G1 split), L3 cross-source (slope monotone in timer strength), plus the maternal-age lineage
(monotonic mother, Kennedy magnitude), the emergent-lifespan calibration, and the
Aqua / ExplicitImports / JET quality gates. Self-contained — no domain dependencies.

## Status

Spun out of the `yeast-wcm-workshop` incubator. **Private pending the public-release gate**;
when released it becomes its own public, citable repo with CI and a Documenter site, destined
for the General registry (see `RELEASE.md` and `CHANGELOG.md`).
