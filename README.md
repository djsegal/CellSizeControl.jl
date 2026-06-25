# CellSizeControl.jl

Generic cell-size-control primitives in Julia: the **sizer / adder / timer** rules,
the model-agnostic **slope discriminator**, and the budding-yeast **inhibitor-dilution
sizer**.

```julia
using CellSizeControl

# Which control law does a lineage obey?  (timer→2, adder→1, sizer→0)
s = simulate_lineage(AdderRule(1.0); n=600)
classify_control(size_control_slope(s.Vb, s.Vd))      # :adder

# Whi5-style inhibitor dilution: Start at V* = W/[Whi5]* — a mechanistic sizer
ids = InhibitorDilutionSizer(60.0, 1.5)               # V* = 40
setpoint_volume(ids)                                  # 40.0
```

## What it is

The size-control question — *what division rule keeps a growing, dividing
population's size distribution stable?* — is generic and reusable, so it lives in
its own package. Three control laws (`TimerRule`, `AdderRule`, `SizerRule`), the
mechanistic `InhibitorDilutionSizer`, a lineage simulator, and the Soifer–Amir
2016 `Vd`-vs-`Vb` slope classifier (timer 2 / adder 1 / sizer 0).

## Science grounding
- **Slope discriminator:** Soifer, Robert & Amir 2016, *Curr Biol* (arXiv:1410.4771).
- **Inhibitor-dilution sizer:** Schmoller, Turner, Kõivomägi & Skotheim 2015,
  *Nature* (V\* = W/[Whi5]\*); mechanistic ODEs in Chandler-Brown/Schmoller 2018
  (BioModels MODEL1803220001) are the planned heavyweight variant.
- Full source notes: the parent workshop's `docs/literature/AUDIT_NOTES.md` §4.

## Testing (three-layer)
`julia --project=. -e 'using Pkg; Pkg.test()'` — L1 analytic (slope recovers each
regime; V\*=W/thresh), L2 reference (a sub-doubling timer collapses a lineage to 0
while the inhibitor-dilution sizer holds it at V\*/2), L3 cross-source (slope
monotone in timer strength). Self-contained — no domain dependencies.

## Status
Spun out of the `yeast-wcm-workshop` incubator (Track-A class artifact + Track-B
package). **Private pending the public-pop gate** (see the workshop's
`docs/decisions/0003-public-pop-strategy.md`); when popped it becomes its own
public, citable repo with CI. Roadmap: add the Chandler-Brown/Schmoller 2018
mechanistic Whi5:SBF titration ODEs + a two-phase (G1-sizer + post-Start-timer)
composition reproducing the whole-cycle adder.
