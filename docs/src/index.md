```@meta
CurrentModule = CellSizeControl
```

# CellSizeControl.jl

Generic, dependency-free primitives for **budding-yeast cell-size control,
asymmetric division, and replicative aging**. The package implements the
sizer / adder / timer rules, the model-agnostic Soifer--Amir slope discriminator,
the Whi5 inhibitor-dilution sizer (and a mechanistic bistable Whi5:SBF switch that
reproduces it from first principles), a two-step-G1 cell cycle, and a maternal-age
aging layer in which one age-eroding division asymmetry ties daughter size,
inherited damage, and the replicative lifespan to a single mechanism.

It depends only on `Random` and `Statistics`.

```@docs
CellSizeControl
```

## Quick start

```julia
using CellSizeControl

# the slope discriminator: timer -> 2, adder -> 1, sizer -> 0
s = simulate_lineage(SizerRule(40.0); n = 600)
classify_control(size_control_slope(s.Vb, s.Vd))      # :sizer

# the inhibitor-dilution sizer: V* = W / theta
setpoint_volume(InhibitorDilutionSizer(60.0, 1.5))    # 40.0

# the mechanistic Whi5:SBF switch reproduces it (theta emergent)
sw = Whi5SBFSwitch(18.0)
setpoint_volume(sw), whi5_sbf_threshold(sw)           # (~40, ~0.45)

# the emergent replicative lifespan
using Statistics
mean(lifespan_distribution(2000))                     # ~25 divisions
```

## Size-control rules

```@docs
SizeControlRule
TimerRule
AdderRule
SizerRule
InhibitorDilutionSizer
SaturatingTimerRule
LinearSizeControl
Whi5SBFSwitch
division_volume
setpoint_volume
whi5_sbf_steady
whi5_sbf_threshold
saturating_timer_buds
```

## The slope discriminator

```@docs
size_control_slope
classify_control
```

## Lineages, asymmetry, and aging

```@docs
simulate_lineage
aging_daughter_fraction
simulate_aging_lineage
replicative_lifespan
damage_trajectory
lifespan_distribution
simulate_population
```

## Single-cell growth and the cell cycle

```@docs
qss_growth_rate
exponential_growth_rate
grow_to
grow_for
cell_cycle
lineage_timecourse
```
