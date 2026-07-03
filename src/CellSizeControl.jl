"""
    CellSizeControl

Generic cell-size-control primitives: the sizer / adder / timer rules, the
model-agnostic **slope discriminator** (Soifer, Robert & Amir 2016 — `Vd` vs `Vb`
slope: timer→2, adder→1, sizer→0), and the budding-yeast **inhibitor-dilution
sizer** (Schmoller, Turner, Kõivomägi & Skotheim 2015: Start at `V* = W/[Whi5]*`).

Spun out of the yeast-wcm workshop (the cell-size-control module — Track-A class
artifact + Track-B package). Self-contained (no domain deps): the size-control
question — *what rule keeps a dividing population's size distribution stable?* —
is generic. See `docs/THREE_LAYER.md` for the analytic/reference/cross-source
testing and `../../docs/literature/AUDIT_NOTES.md` §4 for the source grounding.
"""
module CellSizeControl

using Statistics: mean, std
using Random: Random

# Vendored resampling kernel (bootstrap / BCa CIs on the Monte-Carlo size-law
# statistics); provenance header at the top of the file.
include("ResampleStats.jl")
using .ResampleStats: bootstrap_ci, bca_ci, jackknife

export SizeControlRule,
    TimerRule,
    AdderRule,
    SizerRule,
    InhibitorDilutionSizer,
    SaturatingTimerRule,
    LinearSizeControl,
    Whi5SBFSwitch,
    division_volume,
    setpoint_volume,
    whi5_sbf_steady,
    whi5_sbf_threshold,
    saturating_timer_buds,
    simulate_lineage,
    aging_daughter_fraction,
    simulate_aging_lineage,
    replicative_lifespan,
    damage_trajectory,
    lifespan_distribution,
    simulate_population,
    newborn_size_law,
    extant_size_law,
    senescence_age_law,
    qss_growth_rate,
    exponential_growth_rate,
    grow_to,
    grow_for,
    cell_cycle,
    lineage_timecourse,
    size_control_slope,
    classify_control,
    map_slope,
    size_memory,
    size_law_ci,
    bootstrap_ci,
    bca_ci,
    jackknife

# ---------------------------------------------------------------------------
# Control rules: given the birth volume Vb, the (deterministic) division volume.
# ---------------------------------------------------------------------------
"""
    SizeControlRule

Abstract supertype for a cell-size-control rule: a map from birth volume `Vb` to division
volume `Vd` via [`division_volume`](@ref). Concrete subtypes are [`TimerRule`](@ref),
[`AdderRule`](@ref), [`SizerRule`](@ref), [`InhibitorDilutionSizer`](@ref),
[`SaturatingTimerRule`](@ref), [`LinearSizeControl`](@ref), and [`Whi5SBFSwitch`](@ref).
"""
abstract type SizeControlRule end

"""Timer: divide after a fixed time → `Vd = fold·Vb` (fold = e^{αT}; =2 at doubling).
Slope(Vd,Vb) = fold. A sub-doubling fold (<2) collapses a symmetric-division lineage."""
struct TimerRule <: SizeControlRule
    fold::Float64
end

"""Adder: add a fixed volume increment → `Vd = Vb + Δ`. Slope = 1 (the budding-yeast
whole-cycle phenotype)."""
struct AdderRule <: SizeControlRule
    delta::Float64
end

"""Sizer: divide at a fixed target volume → `Vd = V*`. Slope = 0 (perfect size control)."""
struct SizerRule <: SizeControlRule
    Vstar::Float64
end

"""Inhibitor-dilution sizer (Schmoller 2015): a fixed amount `W` of a size-independent
inhibitor (Whi5) dilutes as volume grows; Start fires at `[inhibitor] = W/V ≤ thresh`,
i.e. `V* = W/thresh`. A mechanistic sizer."""
struct InhibitorDilutionSizer <: SizeControlRule
    W::Float64
    thresh::Float64
end

"""Saturating timer (the course-model failure mode): a size-independent timer fires Start,
and growth is a saturating increment `ΔV = g·(A − Vb)` toward a fixed asymptote `A`, so
`Vd = Vb + g·(A − Vb)`. The bud (= the increment `ΔV`) then SHRINKS as the mother approaches
`A` — the measured daughter-size drift (≈35%→5%) that an open-loop timer produces and a
sizer fixes. (Absorbed from the yeast-wcm `CellSize.jl` A1 module.)"""
struct SaturatingTimerRule <: SizeControlRule
    g::Float64
    asymptote::Float64
end

"""
    setpoint_volume(rule) -> Float64

The target division volume `V*` of a size-controlling rule (the sizers and the bistable
switch): `SizerRule`'s `Vstar`, the inhibitor-dilution `W/thresh`, or the emergent
`Whi5SBFSwitch` set-point. Defined only for rules that have a fixed set-point.

```jldoctest
julia> setpoint_volume(SizerRule(40.0))
40.0

julia> setpoint_volume(InhibitorDilutionSizer(720.0, 18.0))   # V* = W/thresh
40.0
```
"""
setpoint_volume(r::SizerRule) = r.Vstar
setpoint_volume(r::InhibitorDilutionSizer) = r.W / r.thresh

"""
    division_volume(rule, Vb) -> Float64

The (deterministic) division volume `Vd` for a [`SizeControlRule`](@ref) given birth volume
`Vb`: `fold·Vb` for a timer, `Vb + Δ` for an adder, the set-point `V*` for a sizer, and
`α·Vb + β` for the linear map. The least-squares slope of `Vd` on `Vb` is the size-control
discriminator ([`size_control_slope`](@ref)).

```jldoctest
julia> division_volume(AdderRule(10.0), 5.0)   # adder: Vb + Δ
15.0

julia> division_volume(TimerRule(2.0), 5.0)    # timer: fold·Vb
10.0

julia> division_volume(SizerRule(40.0), 12.0)  # sizer: V*, independent of Vb
40.0
```
"""
division_volume(r::TimerRule, Vb) = r.fold * Vb
division_volume(r::AdderRule, Vb) = Vb + r.delta
division_volume(r::SizerRule, Vb) = r.Vstar
division_volume(r::InhibitorDilutionSizer, Vb) = setpoint_volume(r)
division_volume(r::SaturatingTimerRule, Vb) = Vb + r.g * (r.asymptote - Vb)

"""
    saturating_timer_buds(rule::SaturatingTimerRule; V0, n) -> Vector{Float64}

The bud (daughter) birth volumes for a saturating open-loop timer: the bud is the growth
INCREMENT `ΔV = g·(A − V)` each cycle, which shrinks as the mother creeps toward `A`. This is
the size-control failure the inhibitor-dilution sizer corrects; see [`InhibitorDilutionSizer`].
"""
function saturating_timer_buds(rule::SaturatingTimerRule; V0::Real=20.0, n::Int=8)
    V = float(V0)
    buds = Float64[]
    for _ in 1:n
        dV = rule.g * (rule.asymptote - V)
        push!(buds, dV)
        V += dV
    end
    return buds
end

"""Linear size-control map (Amir 2014): `Vd = α·Vb + β`. The continuous family that unifies the
discrete rules — sizer (α=0, β=V*), adder (α=1, β=Δ), timer (α=2, β=0) — so a single `α` sweeps
the sizer↔adder↔timer axis. The measured `Vd`-vs-`Vb` slope recovers `α`; a lineage with
division asymmetry `f` (daughter fraction) is size-homeostatic iff `α·f < 1` (i.e. `α < 1/f`)."""
struct LinearSizeControl <: SizeControlRule
    alpha::Float64
    beta::Float64
end

division_volume(r::LinearSizeControl, Vb) = r.alpha * Vb + r.beta

# ---------------------------------------------------------------------------
# Mechanistic Whi5:SBF bistable Start switch (the inhibitor-dilution sizer from
# first principles). SBF activity x ∈ [0,1] is repressed by *effective* Whi5; active
# SBF (via Cln1,2) in turn inactivates Whi5 — a double-negative feedback that is
# bistable. Whi5 concentration `c = W/V` dilutes as the cell grows; the G1/OFF state
# disappears at a saddle-node `c*`, firing Start, so the set-point V* = W/c* EMERGES.
#   effective Whi5:  e(x,c) = c / (1 + (x/Kx)^p)        (SBF inactivates Whi5)
#   dx/dt          = β / (1 + (e/Ke)^q) − γ·x            (Whi5 represses SBF)
# ---------------------------------------------------------------------------
"""Mechanistic Whi5:SBF bistable Start switch: the inhibitor-dilution sizer from first
principles. A double-negative feedback (Whi5 represses SBF; active SBF inactivates Whi5)
makes SBF activity bistable; growth dilutes Whi5 (`c = W/V`) until the OFF/G1 state vanishes
at a saddle-node `c*`, so the set-point `V* = W/c*` emerges rather than being imposed. Use
the keyword constructor [`Whi5SBFSwitch`](@ref); `Vstar` is precomputed."""
struct Whi5SBFSwitch <: SizeControlRule
    W::Float64        # total Whi5 (size-independent amount synthesized per cycle)
    beta::Float64     # max SBF activation rate
    gamma::Float64    # SBF deactivation rate
    Ke::Float64       # effective-Whi5 repression threshold
    q::Float64        # Whi5-repression Hill coefficient
    Kx::Float64       # SBF→Whi5 inactivation threshold (the double-negative arm)
    p::Float64        # SBF-inactivation Hill coefficient
    Vstar::Float64    # emergent set-point V* = W/c* (precomputed at construction)
end

_whi5_effective(x, c, Kx, p) = c / (1 + (x / Kx)^p)

function _sbf_dxdt(x, c, beta, gamma, Ke, q, Kx, p)
    e = _whi5_effective(x, c, Kx, p)
    return beta / (1 + (e / Ke)^q) - gamma * x
end

function _sbf_steady(c, x0, beta, gamma, Ke, q, Kx, p; dt=0.005, tmax=4000.0)
    x = float(x0)
    t = 0.0
    while t < tmax
        dx = _sbf_dxdt(x, c, beta, gamma, Ke, q, Kx, p)
        x += dx * dt
        x = clamp(x, 0.0, 1.0e3)
        abs(dx) < 1.0e-10 && break
        t += dt
    end
    return x
end

"""
    whi5_sbf_steady(rule, c; from_high=false) -> Float64

Steady-state SBF activity at a fixed Whi5 concentration `c`, integrated from the OFF branch
(`from_high=false`, x₀=0) or the ON branch (`from_high=true`, x₀=1). The two differ across
the bistable window — that hysteresis is the irreversible Start switch.
"""
function whi5_sbf_steady(rule::Whi5SBFSwitch, c::Real; from_high::Bool=false)
    return _sbf_steady(
        float(c),
        from_high ? 1.0 : 0.0,
        rule.beta,
        rule.gamma,
        rule.Ke,
        rule.q,
        rule.Kx,
        rule.p,
    )
end

# Grow V (so c = W/V falls), tracking SBF activity by continuation from the OFF branch; the
# volume where it jumps past the halfway mark to the ON state is the emergent set-point V*.
function _whi5_sbf_start_volume(
    W, beta, gamma, Ke, q, Kx, p; V0=1.0, Vmax=1.0e6, dlogV=0.002
)
    xhigh = _sbf_steady(1.0e-6, 1.0, beta, gamma, Ke, q, Kx, p)   # ON asymptote (c→0)
    thresh = 0.5 * xhigh
    V = float(V0)
    x = _sbf_steady(W / V, 0.0, beta, gamma, Ke, q, Kx, p)
    while V < Vmax
        Vprev = V
        V *= (1 + dlogV)
        x = _sbf_steady(W / V, x, beta, gamma, Ke, q, Kx, p)      # continuation from OFF
        x > thresh && return sqrt(Vprev * V)
    end
    return V
end

"""
    Whi5SBFSwitch(W; beta=1.0, gamma=1.0, Ke=0.30, q=4.0, Kx=0.40, p=4.0)

Mechanistic Whi5:SBF Start switch with total Whi5 amount `W`. The emergent set-point
`V* = W/c*` (the saddle-node concentration `c*` where the OFF/G1 state disappears) is
precomputed. The default switch parameters are bistable with `c* ≈ 0.449`, so `V* ≈ 2.225·W`
— linear in `W`, which is exactly the inhibitor-dilution sizer law `V* = W/θ` with the
threshold `θ = c*` now EMERGENT from a bistable mechanism rather than imposed. See
[`whi5_sbf_threshold`](@ref) and [`InhibitorDilutionSizer`](@ref).
"""
function Whi5SBFSwitch(
    W::Real;
    beta::Real=1.0,
    gamma::Real=1.0,
    Ke::Real=0.30,
    q::Real=4.0,
    Kx::Real=0.40,
    p::Real=4.0,
)
    Vstar = _whi5_sbf_start_volume(
        float(W), float(beta), float(gamma), float(Ke), float(q), float(Kx), float(p)
    )
    return Whi5SBFSwitch(
        float(W), float(beta), float(gamma), float(Ke), float(q), float(Kx), float(p), Vstar
    )
end

setpoint_volume(r::Whi5SBFSwitch) = r.Vstar
division_volume(r::Whi5SBFSwitch, Vb) = r.Vstar

"""
    whi5_sbf_threshold(rule) -> Float64

The emergent Start threshold concentration `θ = c* = W/V*`: the Whi5 concentration at the
saddle-node where the G1/OFF state disappears. An [`InhibitorDilutionSizer`](@ref)`(rule.W, θ)`
has the same set-point — the mechanistic switch reproduces `V* = W/θ` from first principles.
"""
whi5_sbf_threshold(rule::Whi5SBFSwitch) = rule.W / rule.Vstar

# ---------------------------------------------------------------------------
# Lineage simulation (symmetric or asymmetric division, multiplicative noise).
# ---------------------------------------------------------------------------
"""
    simulate_lineage(rule; V0=1.0, n=400, cv=0.08, daughter_fraction=0.5, seed=1)
        -> (; Vb, Vd)

Follow one cell line for `n` generations: born at `Vb`, divide at a noisy
`division_volume(rule, Vb)` (lognormal-ish, floored above `Vb` so it must grow),
the tracked daughter is born at `daughter_fraction·Vd`. Returns the birth- and
division-volume series for the size-control regression.
"""
function simulate_lineage(
    rule::SizeControlRule;
    V0::Real=1.0,
    n::Int=400,
    cv::Real=0.08,
    daughter_fraction::Real=0.5,
    seed::Int=1,
)
    rng = Random.MersenneTwister(seed)
    Vb = Float64[]
    Vd = Float64[]
    v = float(V0)
    for _ in 1:n
        push!(Vb, v)
        d = division_volume(rule, v) * (1 + cv * randn(rng))
        d = max(d, v * 1.001)                 # division can't precede growth
        push!(Vd, d)
        v = daughter_fraction * d
    end
    return (; Vb, Vd)
end

# ---------------------------------------------------------------------------
# Maternal-age asymmetry erosion (replicative aging): old mothers -> larger daughters.
# ---------------------------------------------------------------------------
"""
    aging_daughter_fraction(age; alpha0=0.32, alpha_max=0.5, tau=10.0) -> Float64

Maternal-age erosion of division asymmetry. A young mother divides strongly
asymmetrically (the daughter inherits a small fraction `alpha0` of the division
volume); as replicative `age` rises the asymmetry erodes toward symmetric division
(`alpha_max`, ≈0.5), modeling the age-dependent loss of polarity/segregation control
(the maternal-age decline in division fidelity underlying replicative aging):
`alpha(age) = alpha0 + (alpha_max - alpha0)·(1 - e^{-age/tau})`.

```jldoctest
julia> aging_daughter_fraction(0)   # a young mother divides strongly asymmetrically (α0)
0.32

julia> aging_daughter_fraction(1000)   # very old: asymmetry erodes toward symmetric (α_max)
0.5
```
"""
function aging_daughter_fraction(
    age::Real; alpha0::Real=0.32, alpha_max::Real=0.5, tau::Real=10.0
)
    return alpha0 + (alpha_max - alpha0) * (1 - exp(-age / tau))
end

"""
    simulate_aging_lineage(rule; V0=5.0, n=25, alpha0=0.32, alpha_max=0.5, tau=10.0,
                           damage_form=1.0, damage_cv=0.0,
                           enlarge_max=0.0, enlarge_tau=8.0, phantom_founder=false,
                           cv=0.0, seed=1)
        -> (; gen, Vbirth, Vdivision, Vdaughter, Ddaughter, phantom)

`enlarge_max` adds the maternal-enlargement face: the division set-point rises with
replicative age as `V*(a) = V*·(1 + enlarge_max·(1−e^(−a/enlarge_tau)))` (old mothers are
larger; Johnston 1966, Yang et al. 2011). With `enlarge_max=0` (default) the set-point is fixed and only the
asymmetry erodes — the two faces (size + fitness) then come purely from `α(a)`.

Follow ONE mother through `n` replicative divisions with a maternal-age-dependent
division asymmetry ([`aging_daughter_fraction`](@ref)). Unlike [`simulate_lineage`](@ref)
(a fixed daughter fraction), the daughter's birth volume RISES with the mother's
replicative age, the documented old-mother → larger-daughter relation (Johnston 1966; Yang et al. 2011).
The same fraction also governs the daughter's INHERITED DAMAGE (`Ddaughter`), so one age-eroding asymmetry
drives both daughter size and fitness (larger AND shorter-lived old-mother daughters; Kennedy 1994).
The mother keeps her cell body at division (monotonic, never shrinks); only the bud leaves.
Returns per-generation series.
"""
function simulate_aging_lineage(
    rule::SizeControlRule;
    V0::Real=5.0,
    n::Int=25,
    alpha0::Real=0.32,
    alpha_max::Real=0.5,
    tau::Real=10.0,
    damage_form::Real=1.0,
    damage_cv::Real=0.0,
    enlarge_max::Real=0.0,
    enlarge_tau::Real=8.0,
    phantom_founder::Bool=false,
    cv::Real=0.0,
    seed::Int=1,
)
    rng = Random.MersenneTwister(seed)
    gen = Int[]
    Vbirth = Float64[]
    Vdivision = Float64[]
    Vdaughter = Float64[]
    Ddaughter = Float64[]
    phantom = Bool[]
    vm = float(V0)
    dm = 0.0                                   # mother damage pool
    if phantom_founder
        # The founder is the only cell that would otherwise begin already as a mother body;
        # every later cell is born as a daughter. Prepend her own birth (gen 0) as a daughter
        # of a phantom mother (pristine: no inherited damage), so the lineage is uniform.
        # The phantom mother's own birth/division volumes are undefined (NaN); the founder is
        # born at V0. Off by default so the documented per-generation contract is unchanged.
        push!(gen, 0)
        push!(Vbirth, NaN)
        push!(Vdivision, NaN)
        push!(Vdaughter, float(V0))
        push!(Ddaughter, 0.0)
        push!(phantom, true)
    end
    for a in 0:(n - 1)
        # maternal enlargement: the size set-point rises with replicative age (old mothers
        # are larger; Johnston 1966, Yang et al. 2011), V*(a) = V*·(1 + enlarge_max·(1−e^(−a/enlarge_tau))).
        # enlarge_max=0 (default) keeps the fixed set-point -- the documented contract.
        grow = 1.0 + enlarge_max * (1.0 - exp(-a / enlarge_tau))
        d = division_volume(rule, vm) * grow * (1 + cv * randn(rng))
        d = max(d, vm)                         # d = mother body at division (never below vm: no shrink)
        frac = aging_daughter_fraction(a; alpha0=alpha0, alpha_max=alpha_max, tau=tau)
        # CORRECT division accounting: the mother KEEPS her cell body (monotonic, never
        # shrinks); the daughter is the BUD, a rising fraction `frac(a)` of the enlarging
        # mother (bigger mother feeds a bigger bud => division gets less asymmetric, the
        # documented direction), NOT a slice carved out of the mother. The same `frac(a)` sets
        # how much accrued damage the daughter inherits -- one mechanism, two faces.
        vdau = frac * d
        # per-cycle damage production; damage_cv adds multiplicative noise (floored at 0, so
        # production is non-negative) → inherited damage gets a real distribution across an
        # ensemble. damage_cv=0 (default) keeps the deterministic per-generation accrual.
        dm += damage_form * (damage_cv > 0 ? max(0.0, 1 + damage_cv * randn(rng)) : 1.0)
        push!(gen, a + 1)
        push!(Vbirth, vm)
        push!(Vdivision, d)
        push!(Vdaughter, vdau)
        push!(Ddaughter, frac * dm)
        push!(phantom, false)
        vm = d                                 # mother keeps her body -> next cycle's start
    end
    return (; gen, Vbirth, Vdivision, Vdaughter, Ddaughter, phantom)
end

# ---------------------------------------------------------------------------
# Emergent replicative lifespan from an autocatalytic-damage + viability threshold.
# Replaces a hard-coded generation cap with a lifespan that EMERGES from the dynamics.
# ---------------------------------------------------------------------------
"""
    replicative_lifespan(; D_crit=38.0, crit_cv=0.45, production=1.0, kappa=0.03, cv=0.05,
                         alpha0=0.32, alpha_max=0.5, tau=10.0, segregate=false,
                         D0=0.0, seed=1, max_gen=500) -> Int

The replicative lifespan (number of divisions before senescence) emerging from accumulating
damage rather than being imposed. The default model (the one the paper uses) is
**non-conserved**: at each division the mother's own damage `D` accumulates autocatalytically
and is never depleted by what the daughter inherits — the production rate rises with the damage
already present, `P(D)=production·(1+kappa·D)` (degraded proteostasis begets more damage; Lindner
2008, Hughes & Gottschling 2012), and the full increment is added to `D`. The daughter inheriting
a separate `[r/r_max]·D` is tracked elsewhere (the lineage model) and does not subtract from the
mother's `D` here. The mother senesces when her accumulated damage `D` crosses her viability
threshold; the returned generation count is the RLS. Because production is autocatalytic, `D`
accelerates and the lifespan is finite.

The keyword `segregate=true` selects an optional **conserved-partition** variant — not what the
paper uses — in which the mother retains only the share she does not segregate to the bud,
scaling the per-division increment by `1 - r(a)`, where `r(a)` is the age-eroding division
asymmetry ([`aging_daughter_fraction`](@ref)). The default `segregate=false` keeps the
non-conserved model consistent with the printed equation `P(D)=production·(1+kappa·D)`.

The RLS distribution is set mostly by **cell-to-cell heterogeneity** (`crit_cv`, a lognormal
spread in the threshold `D_crit` across cells) rather than the small per-division noise `cv`:
the autocatalytic blow-up synchronizes the threshold crossing, so per-division noise alone gives
an unrealistically tight distribution. With the default parameters the mean and spread calibrate
to the measured budding-yeast RLS (mean ≈ 25 divisions, CV ≈ 0.3; Schnitzer 2022) — the
threshold/heterogeneity are illustrative, chosen to reproduce that target, not fit per-cell.

The keyword `D0` seeds the cell with **inherited damage at birth** (default `0.0`, a pristine
founder). A daughter that inherits `D0>0` from an aged mother begins partway up the autocatalytic
damage trajectory and senesces sooner — combined with [`damage_trajectory`](@ref) this is the
basis for the convex daughter-RLS-vs-maternal-age prediction.

# Examples
```jldoctest
julia> replicative_lifespan(; cv=0.0, crit_cv=0.0)
26

julia> replicative_lifespan(; D0=20.0, cv=0.0, crit_cv=0.0) < replicative_lifespan(; cv=0.0, crit_cv=0.0)
true
```
"""
function replicative_lifespan(;
    D_crit::Real=38.0,
    crit_cv::Real=0.45,
    production::Real=1.0,
    kappa::Real=0.03,
    cv::Real=0.05,
    alpha0::Real=0.32,
    alpha_max::Real=0.5,
    tau::Real=10.0,
    segregate::Bool=false,
    D0::Real=0.0,
    seed::Int=1,
    max_gen::Int=500,
)
    rng = Random.MersenneTwister(seed)
    a, _ = _age_to_senescence(
        rng,
        D0;
        D_crit,
        crit_cv,
        production,
        kappa,
        cv,
        alpha0,
        alpha_max,
        tau,
        segregate,
        max_gen,
        record=false,
    )
    return a
end

# Shared core for the autocatalytic-damage aging recursion. Age a cell from birth damage `D0`
# to senescence; return (rls, trajectory). When `record`, trajectory[a+1] is the damage carried
# INTO the division that buds the age-`a` daughter (the mother's D_m(a)). Used by both
# `replicative_lifespan` (rls only) and `damage_trajectory` (the D_m(a) series).
function _age_to_senescence(
    rng::Random.AbstractRNG,
    D0::Real;
    D_crit::Real,
    crit_cv::Real,
    production::Real,
    kappa::Real,
    cv::Real,
    alpha0::Real,
    alpha_max::Real,
    tau::Real,
    segregate::Bool,
    max_gen::Int,
    record::Bool,
)
    # cell-to-cell heterogeneity: this cell's own viability threshold (lognormal, mean D_crit)
    Dc = crit_cv > 0 ? D_crit * exp(crit_cv * randn(rng) - crit_cv^2 / 2) : float(D_crit)
    traj = Float64[]
    D = float(D0)
    a = 0
    while D < Dc && a < max_gen
        record && push!(traj, D)
        frac = aging_daughter_fraction(a; alpha0=alpha0, alpha_max=alpha_max, tau=tau)
        kept = segregate ? (1 - frac) : 1.0     # the share the mother does NOT pass to the bud
        noise = cv > 0 ? (1 + cv * randn(rng)) : 1.0
        D += kept * production * (1 + kappa * D) * max(0.0, noise)   # autocatalytic damage
        a += 1
    end
    return a, traj
end

"""
    damage_trajectory(; D0=0.0, seed=1, kwargs...) -> Vector{Float64}

The mother's autocatalytic damage carried into each successive division: element `a+1` is the
damage present when she buds her age-`a` daughter (her `D_m(a)`). `length` equals her
[`replicative_lifespan`](@ref) under the same parameters/seed. A daughter of an age-`a` mother
inherits a share of `damage_trajectory[a+1]` and — via [`replicative_lifespan`](@ref)`(; D0=…)` —
her own shortened emergent lifespan; iterating over the maternal lifespan yields the (convex)
daughter-RLS-vs-maternal-age curve. Keyword arguments match [`replicative_lifespan`](@ref).

# Examples
```jldoctest
julia> traj = damage_trajectory(; cv=0.0, crit_cv=0.0);

julia> length(traj) == replicative_lifespan(; cv=0.0, crit_cv=0.0)
true

julia> traj[1] == 0.0            # a fresh mother buds her first daughter carrying no damage
true

julia> issorted(traj)            # non-conserved damage only accumulates
true
```
"""
function damage_trajectory(;
    D_crit::Real=38.0,
    crit_cv::Real=0.45,
    production::Real=1.0,
    kappa::Real=0.03,
    cv::Real=0.05,
    alpha0::Real=0.32,
    alpha_max::Real=0.5,
    tau::Real=10.0,
    segregate::Bool=false,
    D0::Real=0.0,
    seed::Int=1,
    max_gen::Int=500,
)
    rng = Random.MersenneTwister(seed)
    _, traj = _age_to_senescence(
        rng,
        D0;
        D_crit,
        crit_cv,
        production,
        kappa,
        cv,
        alpha0,
        alpha_max,
        tau,
        segregate,
        max_gen,
        record=true,
    )
    return traj
end

"""
    lifespan_distribution(n; seed0=1, kwargs...) -> Vector{Int}

`n` independent replicative-lifespan samples (distinct seeds), forwarding `kwargs` to
[`replicative_lifespan`](@ref). The distribution (mean + spread) is what calibrates the
damage model to a measured RLS distribution (Schnitzer 2022).
"""
function lifespan_distribution(n::Int; seed0::Int=1, kwargs...)
    return [replicative_lifespan(; seed=seed0 + i - 1, kwargs...) for i in 1:n]
end

# ---------------------------------------------------------------------------
# Exponentially growing POPULATION → steady-state replicative-age structure.
# The single-lineage view (above) follows one mother; a growing culture is the
# ensemble of ALL her descendants. In balanced exponential growth every viable
# cell divides once per generation (buds one age-0 daughter and advances a→a+1),
# so the population doubles and the replicative-AGE distribution converges to the
# geometric law P(age = a) = 2^{-(a+1)} (mean age 1): half the cells are virgin
# daughters, a quarter have budded once, … (Hartwell & Unger 1977; Lord & Wheals
# 1980). The structure is size-rule-independent; with `enlarge_max>0` the rare
# old mothers bud the largest daughters, right-skewing the newborn-size distribution.
# ---------------------------------------------------------------------------
"""
    simulate_population(rule; target=100_000, enlarge_max=0.0, enlarge_tau=8.0,
                        alpha0=0.32, alpha_max=0.5, tau=10.0, D_crit=38.0, crit_cv=0.45,
                        production=1.0, kappa=0.03, cv=0.05, max_gen=200, seed=1)
        -> (; age, rls, Vbirth, ngen)

Grow a synchronous exponential population from a single virgin founder until it reaches
`target` cells, and return the per-cell state of the final population: replicative `age`,
intrinsic replicative lifespan `rls` (drawn per cell via [`replicative_lifespan`](@ref)), and
birth volume `Vbirth`. Each generation every non-senescent cell (age `< rls`) divides: it buds
one age-0 daughter (birth volume `frac(a)·V*(a)`, `frac` = [`aging_daughter_fraction`](@ref) and
`V*(a)` rising with the mother's age when `enlarge_max>0`) and advances to `age+1`, keeping its
(monotonically enlarging) body. Senescent cells persist but stop dividing — at a budding-yeast
mean lifespan (~25) they are a `~2^{-26}` tail and do not perturb the young-age structure.

The steady state is the classic geometric replicative-age distribution `P(age=a) = 2^{-(a+1)}`
(half the population are virgin daughters; mean replicative age 1) — a size-rule-independent
consequence of balanced exponential growth. Use [`simulate_lineage`](@ref) /
[`simulate_aging_lineage`](@ref) for the single-cell-line view instead.

```jldoctest
julia> pop = simulate_population(SizerRule(60.0); target=5000, cv=0.0, crit_cv=0.0, seed=1);

julia> count(==(0), pop.age) == length(pop.age) ÷ 2         # exactly half are virgin daughters
true

julia> a0, a1, a2 = count(==(0), pop.age), count(==(1), pop.age), count(==(2), pop.age);

julia> a1 == a0 ÷ 2 && a2 == a1 ÷ 2                         # geometric: each class halves
true
```
"""
function simulate_population(
    rule::SizeControlRule;
    target::Int=100_000,
    enlarge_max::Real=0.0,
    enlarge_tau::Real=8.0,
    alpha0::Real=0.32,
    alpha_max::Real=0.5,
    tau::Real=10.0,
    D_crit::Real=38.0,
    crit_cv::Real=0.45,
    production::Real=1.0,
    kappa::Real=0.03,
    cv::Real=0.05,
    max_gen::Int=200,
    seed::Int=1,
)
    rng = Random.MersenneTwister(seed)
    draw_rls() = replicative_lifespan(;
        D_crit,
        crit_cv,
        production,
        kappa,
        cv,
        alpha0,
        alpha_max,
        tau,
        seed=rand(rng, 1:typemax(Int32)),
    )
    frac0 = aging_daughter_fraction(0; alpha0, alpha_max, tau)
    age = Int[0]
    rls = Int[draw_rls()]
    Vbirth = Float64[frac0 * division_volume(rule, float(setpoint_hint(rule)))]
    g = 0
    while length(age) < target && g < max_gen
        g += 1
        n = length(age)
        @inbounds for i in 1:n
            a = age[i]
            a < rls[i] || continue                       # senescent: persists, does not divide
            grow = 1.0 + enlarge_max * (1.0 - exp(-a / enlarge_tau))
            Vdiv = division_volume(rule, Vbirth[i]) * grow   # mother's division volume at age a
            frac = aging_daughter_fraction(a; alpha0, alpha_max, tau)
            push!(age, 0)                                # the new virgin daughter
            push!(rls, draw_rls())
            push!(Vbirth, frac * Vdiv)
            age[i] = a + 1                               # mother advances, keeps her body
            Vbirth[i] = Vdiv
        end
    end
    return (; age, rls, Vbirth, ngen=g)
end

# A birth-size seed for the founder's first division volume; for rules whose division volume
# ignores Vb (the sizers) this is exact, for the size-dependent rules it only seeds gen 0.
setpoint_hint(r::SizeControlRule) = 1.0
setpoint_hint(r::SizerRule) = r.Vstar
setpoint_hint(r::InhibitorDilutionSizer) = setpoint_volume(r)
setpoint_hint(r::Whi5SBFSwitch) = r.Vstar

"""
    newborn_size_law(; alpha0=0.32, alpha_max=0.5, tau=10.0, enlarge_max=0.0,
                     enlarge_tau=8.0, Vstar=1.0, max_age=80)
        -> (; mean, sd, cv, skew, ratio)

Closed-form moments of the population **newborn (virgin-daughter) size distribution** predicted
by sampling the age-eroding division asymmetry through the geometric replicative-age law. In
balanced exponential growth the mothers of replicative age `a` are a fixed fraction
`2^{-(a+1)}` of the culture ([`simulate_population`](@ref)), and each buds one age-0 daughter of
size `frac(a)·V*·enlarge(a)` — so the newborn sizes inherit those geometric weights exactly. With
`frac(a)` the age-eroding daughter fraction ([`aging_daughter_fraction`](@ref)) and
`enlarge(a) = 1 + enlarge_max·(1 − e^{−a/enlarge_tau})` the maternal-enlargement of the set-point,
the distribution is the discrete geometric mixture `{(2^{-(a+1)}, frac(a)·V*·enlarge(a))}`.

Returns its mean, standard deviation, coefficient of variation, skewness, and the **scale-free
ratio** `mean / (alpha0·V*)` (the mean newborn size in units of the youngest-mother daughter).
The `cv`, `skew`, and `ratio` are **independent of `V*`** — a scale-free signature of the
mechanism: any monotone maternal-age→daughter-size relation makes the newborn distribution
right-skewed (`skew > 0`); it collapses to a symmetric point mass (`cv = skew = 0`) only when the
daughter size does not depend on maternal age (`alpha_max = alpha0` and `enlarge_max = 0`). This
is the analytic counterpart of the age-0 birth-size histogram from [`simulate_population`](@ref),
which it reproduces to Monte-Carlo precision.

```jldoctest
julia> law = newborn_size_law(; tau=8.0, enlarge_max=0.45, enlarge_tau=8.0, Vstar=60.0);

julia> round(law.ratio, digits=4)          # mean newborn / (α0·V*): a scale-free constant
1.1138

julia> round(law.cv, digits=4), round(law.skew, digits=4)   # right-skewed geometric mixture
(0.1346, 1.562)

julia> newborn_size_law(; alpha0=0.32, alpha_max=0.32, enlarge_max=0.0).skew   # no age effect ⇒ symmetric
0.0
```
"""
function newborn_size_law(;
    alpha0::Real=0.32,
    alpha_max::Real=0.5,
    tau::Real=10.0,
    enlarge_max::Real=0.0,
    enlarge_tau::Real=8.0,
    Vstar::Real=1.0,
    max_age::Int=80,
)
    ws = [2.0^(-(a + 1)) for a in 0:max_age]
    ws ./= sum(ws)                                    # renormalize the truncated geometric law
    ss = [
        aging_daughter_fraction(a; alpha0=alpha0, alpha_max=alpha_max, tau=tau) *
        Vstar *
        (1 + enlarge_max * (1 - exp(-a / enlarge_tau))) for a in 0:max_age
    ]
    mean_s = sum(ws .* ss)
    # central moments directly about the mean (stable when the mixture is near-degenerate)
    sd = sqrt(sum(ws .* (ss .- mean_s) .^ 2))
    cv = mean_s > 0 ? sd / mean_s : 0.0
    # a maternal-age-independent size (alpha_max=alpha0, enlarge_max=0) is a point mass: the
    # spread is then only floating-point noise, so a relative floor guards the skew ratio 0/0.
    skew = cv > 1e-9 ? sum(ws .* (ss .- mean_s) .^ 3) / sd^3 : 0.0
    ratio = alpha0 > 0 ? mean_s / (alpha0 * Vstar) : NaN
    return (; mean=mean_s, sd, cv, skew, ratio)
end

"""
    extant_size_law(; alpha0=0.32, alpha_max=0.5, tau=10.0, enlarge_max=0.0,
                    enlarge_tau=8.0, Vstar=1.0, max_age=80)
        -> (; newborn_mean, extant_mean, divergence)

Closed-form mean size of an **extant** (whole-population, snapshot) cell versus a **newborn**
(virgin daughter) in a balanced exponentially-growing culture, and their ratio `divergence`.

A snapshot of the culture samples every cell at its most-recent division: the age-0 cells are the
buds just shed (the [`newborn_size_law`](@ref) distribution, mean `newborn_mean`), while a cell of
replicative age `a ≥ 1` is a mother carrying her full retained body `V*·enlarge(a−1)` (monotonic —
the mother keeps her body, only the bud leaves). Weighting each age class by the geometric
replicative-age law `P(age=a) = 2^{-(a+1)}` gives

`extant_mean = ½·newborn_mean + Σ_{a≥1} 2^{-(a+1)}·V*·enlarge(a−1)`.

The **divergence** `extant_mean / newborn_mean` is the size-structure signature of exponential
balanced growth: the mean cell in the culture is far larger than the mean newborn because the
standing population over-represents the larger, older mother bodies relative to the small buds they
shed. Like the newborn law it is **scale-free** — independent of the set-point `V*` — and, in the
no-erosion, no-enlargement limit (`alpha_max = alpha0`, `enlarge_max = 0`), it reduces to the exact
closed form `(1 + alpha0) / (2·alpha0)` (2.0625 at `alpha0 = 0.32`). It is the analytic counterpart
of `mean(pop.Vbirth) / mean(pop.Vbirth[pop.age .== 0])` from [`simulate_population`](@ref).

```jldoctest
julia> law = extant_size_law(; tau=8.0, enlarge_max=0.45, enlarge_tau=8.0, Vstar=60.0);

julia> round(law.divergence, digits=4)     # mean extant cell ≈ 1.97× the mean newborn
1.9692

julia> round(extant_size_law(; tau=8.0, enlarge_max=0.45, Vstar=137.0).divergence, digits=4)
1.9692

julia> round(extant_size_law(; alpha0=0.32, alpha_max=0.32, enlarge_max=0.0).divergence, digits=4)
2.0625
```
"""
function extant_size_law(;
    alpha0::Real=0.32,
    alpha_max::Real=0.5,
    tau::Real=10.0,
    enlarge_max::Real=0.0,
    enlarge_tau::Real=8.0,
    Vstar::Real=1.0,
    max_age::Int=80,
)
    ws = [2.0^(-(a + 1)) for a in 0:max_age]
    ws ./= sum(ws)                                    # renormalize the truncated geometric law
    enlarge(a) = 1 + enlarge_max * (1 - exp(-a / enlarge_tau))
    # newborn mean: the geometric-mixture prediction (mirror of newborn_size_law)
    newborn_mean =
        Vstar * sum(
            ws[a + 1] *
            aging_daughter_fraction(a; alpha0=alpha0, alpha_max=alpha_max, tau=tau) *
            enlarge(a) for a in 0:max_age
        )
    # extant snapshot: age-0 → newborn distribution; age a≥1 → retained body V*·enlarge(a−1)
    extant_mean =
        ws[1] * newborn_mean +
        Vstar * sum(ws[a + 1] * enlarge(a - 1) for a in 1:max_age)
    divergence = newborn_mean > 0 ? extant_mean / newborn_mean : NaN
    return (; newborn_mean, extant_mean, divergence)
end

"""
    size_law_ci(newborn_sizes; alpha0=0.32, Vstar=1.0, alpha=0.05, nboot=2000, seed=1)
        -> (; ratio, cv, skew)

BCa bootstrap confidence intervals for the scale-free newborn-size-law statistics from a
Monte-Carlo newborn-size sample — the age-0 birth volumes `pop.Vbirth[pop.age .== 0]` returned
by [`simulate_population`](@ref). Each field is a `(lo, point, hi)` triple at confidence
`1 − alpha`:

  - `ratio` — mean newborn size in units of the youngest-mother daughter, `mean / (alpha0·Vstar)`
    (the scale-free CC-N headline ≈ 1.114 that [`newborn_size_law`](@ref) predicts in closed form);
  - `cv`    — coefficient of variation `std / mean`;
  - `skew`  — the (right-skewed) third standardized moment.

The `point` of each triple **exactly reproduces** the plain moment estimator the analysis
computes by hand, and the closed-form [`newborn_size_law`](@ref) value should fall inside the
interval at large `N` — turning the bare single-seed point comparison into an interval-covered
validation. The interval is bias-corrected and accelerated (skew-aware) via the vendored
`ResampleStats` kernel; deterministic under `seed`.

```jldoctest
julia> pop = simulate_population(SizerRule(60.0); target=40_000, enlarge_max=0.45,
                                 enlarge_tau=8.0, alpha0=0.32, alpha_max=0.5, tau=8.0, seed=1);

julia> nb = pop.Vbirth[pop.age .== 0];

julia> ci = size_law_ci(nb; alpha0=0.32, Vstar=60.0, nboot=1000, seed=1);

julia> ci.ratio[1] < newborn_size_law(; alpha0=0.32, alpha_max=0.5, tau=8.0,
                                      enlarge_max=0.45, enlarge_tau=8.0, Vstar=60.0).ratio < ci.ratio[3]
true
```
"""
function size_law_ci(
    newborn_sizes::AbstractVector{<:Real};
    alpha0::Real=0.32,
    Vstar::Real=1.0,
    alpha::Real=0.05,
    nboot::Integer=2000,
    seed::Integer=1,
)
    length(newborn_sizes) >= 2 ||
        throw(ArgumentError("size_law_ci needs at least 2 newborn sizes"))
    nb = float.(newborn_sizes)
    ratio_stat(x) = mean(x) / (alpha0 * Vstar)
    cv_stat(x) = std(x) / mean(x)
    function skew_stat(x)
        m = mean(x)
        s = std(x)
        return s > 0 ? mean(((x .- m) ./ s) .^ 3) : 0.0
    end
    # independent, reproducible rng streams per statistic
    ratio = bca_ci(nb, ratio_stat; nboot, alpha, rng=Random.MersenneTwister(seed))
    cv = bca_ci(nb, cv_stat; nboot, alpha, rng=Random.MersenneTwister(seed + 1))
    skew = bca_ci(nb, skew_stat; nboot, alpha, rng=Random.MersenneTwister(seed + 2))
    return (; ratio, cv, skew)
end

"""
    senescence_age_law(rls::Integer) -> (; lambda, ages, p)

Senescence correction to the geometric replicative-age law at short mean replicative lifespan.

The clean law `P(age=a) = 2^{-(a+1)}` assumes cells divide forever. With a finite lifespan `rls`
(cells divide at ages `0 … rls−1` then arrest), the dividing population grows by a factor `λ < 2`
per generation set by the discrete Euler–Lotka equation `λ = Σ_{a=0}^{rls−1} λ^{-a}`, and its
replicative-age distribution is the **truncated geometric** `P(age=a) = λ^{-(a+1)}` for
`a = 0 … rls−1` (using `Σ_{a} λ^{-a} = λ`, so the virgin fraction is exactly `1/λ`). As
`rls → ∞`, `λ → 2` and the law recovers `2^{-(a+1)}`; at short `rls` the base `1/λ > 1/2` flattens
the distribution and it truncates at `rls−1`. The one-division-shy limit `rls = 2` gives the golden
ratio, `λ = φ = (1+√5)/2`. This is the deterministic-lifespan analytic counterpart of the age
histogram from [`simulate_population`](@ref) run at a short lifespan.

`ages` is `0:rls-1` and `p[a+1] = λ^{-(a+1)}` (the normalized dividing-cell age distribution).

```jldoctest
julia> round(senescence_age_law(2).lambda, digits=6)     # rls=2 ⇒ golden ratio φ
1.618034

julia> round(senescence_age_law(4).lambda, digits=5)
1.92756

julia> isapprox(senescence_age_law(100).lambda, 2.0; atol=1e-6)   # long rls recovers 2^{-(a+1)}
true

julia> round.(senescence_age_law(4).p, digits=4)         # truncated geometric, sums to 1
4-element Vector{Float64}:
 0.5188
 0.2691
 0.1396
 0.0724
```
"""
function senescence_age_law(rls::Integer)
    rls >= 1 || throw(ArgumentError("rls must be ≥ 1"))
    # Newton solve of λ = Σ_{a=0}^{rls−1} λ^{-a}
    lam = 2.0
    for _ in 1:500
        f = lam - sum(lam^(-a) for a in 0:(rls - 1))
        df = 1 + sum(a * lam^(-a - 1) for a in 0:(rls - 1))
        step = f / df
        lam -= step
        abs(step) < 1e-15 && break
    end
    ages = collect(0:(rls - 1))
    p = [lam^(-(a + 1)) for a in ages]
    return (; lambda=lam, ages, p)
end

# ---------------------------------------------------------------------------
# Energetic single-cell growth + the two-step G1 (Di Talia 2007). With a sizer
# threshold V*, the mother/daughter G1 durations EMERGE rather than being imposed.
# ---------------------------------------------------------------------------
# Surface-area-limited QSS growth (the course VOL_Growth law; Altenburg et al. 2019
# constants), per minute: dV/dt = (k_up·4πr² − k_cons·V), r = (3V/4π)^(1/3).
const _K_UP, _K_CONS, _C_ISS = 0.23, 0.27, 319.4

"""
    qss_growth_rate(V) -> Float64

Surface-area-limited single-cell volume growth rate `dV/dt` (per minute) at volume `V`:
`(k_up·4πr² − k_cons·V)` with `r = (3V/4π)^(1/3)` (the quasi-steady-state uptake-minus-
consumption law; Altenburg et al. 2019 constants). Pass as the `rate` to [`grow_to`](@ref) /
[`cell_cycle`](@ref); contrast [`exponential_growth_rate`](@ref).
"""
function qss_growth_rate(V)
    return 60.0 * (_K_UP * 4.0 * pi * (3.0 * V / (4.0 * pi))^(2 / 3) - _K_CONS * V) / _C_ISS
end

"""Exponential (biomass-driven) growth `dV/dt = μ·V` — the single-cell form measured by
Di Talia 2007 / Sun 2010 (constant specific rate); pass `exponential_growth_rate(μ)` as the
`rate` to [`cell_cycle`](@ref). μ ≈ 0.0077/min is a ~90-min doubling."""
exponential_growth_rate(μ::Real) = V -> μ * V

"""Grow `V` toward `target` under `rate(V)=dV/dt`; return `(minutes, V_reached)`."""
function grow_to(
    V::Real, target::Real; dt::Real=0.02, cap::Real=5000.0, rate=qss_growth_rate
)
    t = 0.0
    V = float(V)
    while V < target && t < cap
        V += max(0.0, rate(V)) * dt
        t += dt
    end
    return t, V
end

"""Grow `V` for `dur` minutes under `rate(V)=dV/dt`; return the final volume."""
function grow_for(V::Real, dur::Real; dt::Real=0.02, rate=qss_growth_rate)
    t = 0.0
    V = float(V)
    while t < dur
        V += max(0.0, rate(V)) * dt
        t += dt
    end
    return V
end

"""
    cell_cycle(Vb; Vstar=40.0, T_cln2=19.0, tau_bud=70.0, bud_seed=2.0, rate=qss_growth_rate)

One cell cycle from birth volume `Vb`. G1 is the inhibitor-dilution **sizer** step (time to
grow `Vb → V*`, zero for a mother already ≥ V*) plus a fixed **Cln2 timer** step `T_cln2`
(Di Talia 2007's two G1 modules); the budded phase is a size-invariant timer `tau_bud` during
which the bud grows on its own geometry from `bud_seed`. So a mother (born ≥ V*) has
G1 ≈ `T_cln2`, while a daughter (born small) spends extra time reaching V* — the
mother/daughter G1 asymmetry EMERGES, it is not imposed. Returns the phase + size readout.

```jldoctest
julia> cell_cycle(40.0; Vstar=40.0).G1   # mother born at V*: no sizer wait, G1 = T_cln2
19.0

julia> cell_cycle(20.0).G1 > cell_cycle(40.0).G1   # a small daughter waits longer in G1
true
```
"""
function cell_cycle(
    Vb::Real;
    Vstar::Real=40.0,
    T_cln2::Real=19.0,
    tau_bud::Real=70.0,
    bud_seed::Real=2.0,
    rate=qss_growth_rate,
)
    t_sizer, _ = grow_to(Vb, Vstar; rate=rate)
    G1 = t_sizer + T_cln2
    Vstart = max(float(Vb), float(Vstar))
    Vdaughter = grow_for(bud_seed, tau_bud; rate=rate)
    return (;
        G1,
        budded=float(tau_bud),
        cycle=G1 + tau_bud,
        Vstart,
        Vdiv=Vstart + Vdaughter,
        Vdaughter,
        Vmother=Vstart,
    )
end

"""
    lineage_timecourse(; n_max=29, Vstar0=36.0, T_cln2=19.0, tau_bud=52.0, m_enlarge=0.45,
                       m_tau=8.0, r0=0.69, r_max=0.90, r_tau=14.0, dt=0.5, rate=qss_growth_rate)
        -> (; t, Vmother, Vbud)

Volume-vs-time trajectory of one mother lineage over `n_max+1` cell cycles (≈ the replicative
lifespan). Each cycle: the mother grows to the enlarging set-point `V*(a)` during G1, holds
through the Cln2 timer, then the bud grows on its own geometry to `V*(a)·ratio(a)` during the
budded phase; at division the bud (daughter) detaches while the **mother keeps her body
(monotonic, never shrinks)**. `t` is in minutes. The single source of the time-view figure.
"""
function lineage_timecourse(;
    n_max::Int=29,
    Vstar0::Real=36.0,
    T_cln2::Real=19.0,
    tau_bud::Real=52.0,
    m_enlarge::Real=0.45,
    m_tau::Real=8.0,
    r0::Real=0.69,
    r_max::Real=0.90,
    r_tau::Real=14.0,
    dt::Real=0.5,
    rate=qss_growth_rate,
)
    Vs(a) = Vstar0 * (1.0 + m_enlarge * (1.0 - exp(-a / m_tau)))
    ratio(a) = r0 + (r_max - r0) * (1.0 - exp(-a / r_tau))
    t = 0.0
    Vm = float(Vstar0)
    ts = Float64[]
    vmo = Float64[]
    vbu = Float64[]
    rec(vb) = (push!(ts, t); push!(vmo, Vm); push!(vbu, vb))
    rec(0.0)
    for a in 0:n_max
        target = Vs(a)
        while Vm < target                         # G1 sizer: mother grows to V*(a)
            Vm = min(target, Vm + max(0.0, rate(Vm)) * dt)
            t += dt
            rec(0.0)
        end
        te = t + T_cln2                            # Cln2 timer (mother holds, 0 slope)
        while t < te
            t += dt
            rec(0.0)
        end
        target_bud = target * ratio(a)             # budded: bud grows to V*(a)·ratio(a)
        t0 = t
        te = t + tau_bud
        while t < te
            s = (t - t0) / tau_bud
            t += dt
            rec(target_bud * s^1.5)
        end
        rec(0.0)                                   # division: bud detaches; mother keeps Vm
    end
    return (; t=ts, Vmother=vmo, Vbud=vbu)
end

# ---------------------------------------------------------------------------
# The Soifer/Amir 2016 discriminator.
# ---------------------------------------------------------------------------
"""
    size_control_slope(Vb, Vd) -> Float64

Least-squares slope of division volume on birth volume. **timer → 2, adder → 1,
sizer → 0** (symmetric division). The model-agnostic size-control classifier.

```jldoctest
julia> Vb = [1.0, 2.0, 3.0];

julia> size_control_slope(Vb, Vb .+ 10)    # adder: Vd = Vb + Δ
1.0

julia> size_control_slope(Vb, [40.0, 40.0, 40.0])   # sizer: Vd = V*
0.0

julia> size_control_slope(Vb, 2 .* Vb)     # timer: Vd = 2·Vb
2.0
```
"""
function size_control_slope(Vb::AbstractVector, Vd::AbstractVector)
    x = float.(Vb)
    y = float.(Vd)
    x̄ = mean(x)
    return sum((x .- x̄) .* (y .- mean(y))) / sum(abs2, x .- x̄)
end

"""
    classify_control(slope; atol=0.35) -> Symbol

Map a `Vd`-vs-`Vb` slope to `:sizer` (0) / `:adder` (1) / `:timer` (2) / `:mixed`.

```jldoctest
julia> classify_control(0.0)
:sizer

julia> classify_control(1.0)
:adder

julia> classify_control(2.0)
:timer

julia> classify_control(0.5)   # between the canonical rules
:mixed
```
"""
function classify_control(slope::Real; atol::Real=0.35)
    isapprox(slope, 0; atol=atol) && return :sizer
    isapprox(slope, 1; atol=atol) && return :adder
    isapprox(slope, 2; atol=atol) && return :timer
    return :mixed
end

# ---------------------------------------------------------------------------
# Lineage birth-size MEMORY: the return map is AR(1) with a single pole r = α·f.
# The birth-size return map Vb(n+1) = f·(division_volume(rule,Vb(n)))·(1+cv·ξ)
# linearizes about its fixed point to Vb(n+1) − Vb* ≈ r·(Vb(n) − Vb*) + noise, an
# AR(1) process whose one pole r = α·f (α = the map's local slope, f the daughter
# fraction) sets EVERY memory observable at once: the stationary CV amplification
# 1/√(1−r²), the lag-k lineage autocorrelation ρ_k = r^k, and the nutrient-shift
# relaxation time −1/ln r (generations). Homeostatic iff r < 1.
# ---------------------------------------------------------------------------
"""
    map_slope(rule::SizeControlRule) -> Float64

The ANALYTIC local slope `α = dVd/dVb` of a size-control map — the closed-form
counterpart of the regression-based [`size_control_slope`](@ref): sizer → 0, adder → 1,
timer → its fold, linear map → `α`. This is the `α` that, with the daughter fraction `f`,
forms the lineage return slope `r = α·f` ([`size_memory`](@ref)).

```jldoctest
julia> map_slope(SizerRule(40.0)), map_slope(AdderRule(10.0)), map_slope(TimerRule(2.0))
(0.0, 1.0, 2.0)

julia> map_slope(LinearSizeControl(1.5, 20.0))
1.5
```
"""
map_slope(::SizerRule) = 0.0
map_slope(::InhibitorDilutionSizer) = 0.0
map_slope(::AdderRule) = 1.0
map_slope(r::TimerRule) = float(r.fold)
map_slope(r::LinearSizeControl) = r.alpha

"""
    size_memory(rule; daughter_fraction=0.5)
        -> (; r, cv_gain, autocorr, memory_gen)

Birth-size memory of a size-control lineage, all set by the single return-map pole
`r = α·f` (`α = map_slope(rule)`, `f = daughter_fraction`). The return map is AR(1), so:

  - `r` — the return slope `α·f` (homeostatic iff `r < 1`);
  - `cv_gain = 1/√(1−r²)` — stationary birth-size CV amplification over the per-division
    noise (the `CV(Vb) = cv/√(1−(αf)²)` law), `Inf` at/above the boundary;
  - `autocorr = r` — the lag-1 mother→daughter birth-size correlation; the lag-`k`
    correlation is `ρ_k = r^k` (geometric decay at the single rate `r`);
  - `memory_gen = −1/ln r` — the nutrient-shift relaxation time in generations (`0` for a
    sizer, which adapts in one division; `Inf` at the boundary).

These force a two-observable invariant measurable from a **single** lineage, mode- and
set-point-free: `CV(Vb)²·(1 − ρ_1²) = cv²` recovers the intrinsic per-division noise.

```jldoctest
julia> m = size_memory(AdderRule(10.0); daughter_fraction=0.5);

julia> m.r, round(m.cv_gain; digits=3), round(m.memory_gen; digits=3)
(0.5, 1.155, 1.443)

julia> size_memory(SizerRule(40.0)).memory_gen        # a sizer is memoryless
0.0

julia> round(size_memory(TimerRule(2.0); daughter_fraction=0.4).memory_gen; digits=3)
4.481
```
"""
function size_memory(rule::SizeControlRule; daughter_fraction::Real=0.5)
    r = map_slope(rule) * daughter_fraction
    cv_gain = r >= 1 ? Inf : 1 / sqrt(1 - r^2)
    memory_gen = r <= 0 ? 0.0 : (r >= 1 ? Inf : -1 / log(r))
    return (; r=float(r), cv_gain, autocorr=float(r), memory_gen)
end

end # module
