# CellSizeControl.jl

Generic cell-size-control primitives in Julia: the **sizer / adder / timer** rules, the
model-agnostic **slope discriminator**, the budding-yeast **inhibitor-dilution sizer**, and
a **maternal-age aging** layer in which one age-eroding division asymmetry ties daughter
size, inherited damage, and the **replicative lifespan** to a single mechanism.

```julia
using CellSizeControl

# Which control law does a lineage obey?  (timer→2, adder→1, sizer→0)
s = simulate_lineage(AdderRule(1.0); n=600)
classify_control(size_control_slope(s.Vb, s.Vd))      # :adder

# Whi5-style inhibitor dilution: Start at V* = W/[Whi5]* — a mechanistic sizer
ids = InhibitorDilutionSizer(60.0, 1.5)               # V* = 40
setpoint_volume(ids)                                  # 40.0

# Two-step G1 (Di Talia 2007): a mother (born ≥ V*) runs ≈ the Cln2 timer; a daughter
# (born small) spends extra time growing to V* — the G1 asymmetry emerges, not imposed.
cell_cycle(40.0; Vstar=36.0).G1                       # ≈ 19 min (mother)

# Maternal-age asymmetry erosion: one r(a) drives daughter SIZE and inherited DAMAGE
L = simulate_aging_lineage(SizerRule(60.0); n=30, enlarge_max=0.45)
L.Vdaughter[end] > L.Vdaughter[1]                     # old mothers → larger daughters (Kennedy)

# The replicative lifespan EMERGES from autocatalytic damage + a viability threshold
using Statistics: mean
mean(lifespan_distribution(2000))                     # ≈ 25 divisions (Schnitzer 2022)
```

## What it is

The size-control question — *what division rule keeps a growing, dividing population's size
distribution stable?* — is generic and reusable, so it lives in its own package. It holds
the three control laws (`TimerRule`, `AdderRule`, `SizerRule`), the mechanistic
`InhibitorDilutionSizer`, the Soifer–Amir 2016 `Vd`-vs-`Vb` slope classifier (timer 2 /
adder 1 / sizer 0), an energetic two-step-G1 `cell_cycle`, and a maternal-age layer
(`simulate_aging_lineage`, `replicative_lifespan`, `lifespan_distribution`) in which a
single age-eroding asymmetry `r(a)` is both the *size* face (larger daughters) and the
*fitness* face (more inherited damage, finite lifespan). The mother keeps her cell body at
division (monotonic, never shrinks); only the bud leaves.

## Science grounding
- **Slope discriminator:** Soifer, Robert & Amir 2016, *Curr Biol* (arXiv:1410.4771).
- **Inhibitor-dilution sizer:** Schmoller, Turner, Kõivomägi & Skotheim 2015, *Nature*
  (V\* = W/[Whi5]\*); mechanistic ODEs in Chandler-Brown/Schmoller 2018 (BioModels
  MODEL1803220001) are the planned heavyweight variant.
- **Two-step G1:** Di Talia, Skotheim, Bean, Siggia & Cross 2007, *Nature* (mother ≈ 19 min,
  daughter ≈ 45 min). **Budded timer:** Leitão & Lucena 2018.
- **Maternal-age daughter size:** Kennedy, Austriaco & Guarente 1994, *J Cell Biol*.
  **Replicative lifespan:** ≈ 25 divisions, CV ≈ 0.3 (Schnitzer et al. 2022).
- Full source notes: the parent workshop's `docs/literature/AUDIT_NOTES.md` §4.

## Testing (three-layer)
`julia --project=. -e 'using Pkg; Pkg.test()'` — L1 analytic (slope recovers each regime;
V\*=W/thresh), L2 reference (a sub-doubling timer collapses a lineage to 0 while the
inhibitor-dilution sizer holds it at V\*/2; the Di Talia G1 split), L3 cross-source (slope
monotone in timer strength), plus the maternal-age lineage (monotonic mother, Kennedy
daughter-size magnitude), the emergent-lifespan calibration (mean + CV vs Schnitzer), and
Aqua / ExplicitImports / JET quality gates. Self-contained — no domain dependencies.

## Status
Spun out of the `yeast-wcm-workshop` incubator. **Private pending the public-pop gate**;
when popped it becomes its own public, citable repo with CI, destined for the General
registry (see `RELEASE.md`). Roadmap: the Chandler-Brown/Schmoller 2018 mechanistic Whi5:SBF
titration ODEs + a two-phase (G1-sizer + post-Start-timer) composition reproducing the
whole-cycle adder.
