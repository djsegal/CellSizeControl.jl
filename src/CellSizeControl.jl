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

using Statistics: mean
using Random: Random

export SizeControlRule,
    TimerRule,
    AdderRule,
    SizerRule,
    InhibitorDilutionSizer,
    SaturatingTimerRule,
    division_volume,
    setpoint_volume,
    saturating_timer_buds,
    simulate_lineage,
    aging_daughter_fraction,
    simulate_aging_lineage,
    replicative_lifespan,
    lifespan_distribution,
    qss_growth_rate,
    exponential_growth_rate,
    grow_to,
    grow_for,
    cell_cycle,
    lineage_timecourse,
    size_control_slope,
    classify_control

# ---------------------------------------------------------------------------
# Control rules: given the birth volume Vb, the (deterministic) division volume.
# ---------------------------------------------------------------------------
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

setpoint_volume(r::SizerRule) = r.Vstar
setpoint_volume(r::InhibitorDilutionSizer) = r.W / r.thresh

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
(Kennedy, Austriaco & Guarente 1994):
`alpha(age) = alpha0 + (alpha_max - alpha0)·(1 - e^{-age/tau})`.
"""
function aging_daughter_fraction(
    age::Real; alpha0::Real=0.32, alpha_max::Real=0.5, tau::Real=10.0
)
    return alpha0 + (alpha_max - alpha0) * (1 - exp(-age / tau))
end

"""
    simulate_aging_lineage(rule; V0=5.0, n=25, alpha0=0.32, alpha_max=0.5, tau=10.0,
                           enlarge_max=0.0, enlarge_tau=8.0, phantom_founder=false,
                           cv=0.0, seed=1)
        -> (; gen, Vbirth, Vdivision, Vdaughter, Ddaughter, phantom)

`enlarge_max` adds the maternal-enlargement face: the division set-point rises with
replicative age as `V*(a) = V*·(1 + enlarge_max·(1−e^(−a/enlarge_tau)))` (old mothers are
larger; Kennedy 1994). With `enlarge_max=0` (default) the set-point is fixed and only the
asymmetry erodes — the two faces (size + fitness) then come purely from `α(a)`.

Follow ONE mother through `n` replicative divisions with a maternal-age-dependent
division asymmetry ([`aging_daughter_fraction`](@ref)). Unlike [`simulate_lineage`](@ref)
(a fixed daughter fraction), the daughter's birth volume RISES with the mother's
replicative age, the documented old-mother → larger-daughter relation (Kennedy 1994). The same fraction
also governs the daughter's INHERITED DAMAGE (`Ddaughter`), so one age-eroding asymmetry
drives both daughter size and fitness (larger AND shorter-lived old-mother daughters).
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
        # are larger; Kennedy 1994), V*(a) = V*·(1 + enlarge_max·(1−e^(−a/enlarge_tau))).
        # enlarge_max=0 (default) keeps the fixed set-point -- the documented contract.
        grow = 1.0 + enlarge_max * (1.0 - exp(-a / enlarge_tau))
        d = division_volume(rule, vm) * grow * (1 + cv * randn(rng))
        d = max(d, vm)                         # d = mother body at division (never below vm: no shrink)
        frac = aging_daughter_fraction(a; alpha0=alpha0, alpha_max=alpha_max, tau=tau)
        # CORRECT division accounting: the mother KEEPS her cell body (monotonic, never
        # shrinks); the daughter is the BUD, a rising fraction `frac(a)` of the enlarging
        # mother (bigger mother feeds a bigger bud => division gets less asymmetric, the
        # Kennedy direction), NOT a slice carved out of the mother. The same `frac(a)` sets
        # how much accrued damage the daughter inherits -- one mechanism, two faces.
        vdau = frac * d
        dm += damage_form
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
    replicative_lifespan(; D_crit=32.0, crit_cv=0.55, production=1.0, kappa=0.10, cv=0.05,
                         alpha0=0.32, alpha_max=0.5, tau=10.0, segregate=true,
                         seed=1, max_gen=500) -> Int

The replicative lifespan (number of divisions before senescence) emerging from accumulating
damage rather than being imposed. Each division adds damage that is **autocatalytic** — the
production rate rises with the damage already present, `P(D)=production·(1+kappa·D)` (degraded
proteostasis begets more damage; Lindner 2008, Hughes & Gottschling 2012) — and the mother
**retains the share she does not segregate** to the bud, `1 - r(a)`, where `r(a)` is the same
age-eroding division asymmetry that drives daughter size and inherited damage
([`aging_daughter_fraction`](@ref); `segregate=false` makes the mother keep all of it). The
mother senesces when her retained damage `D` crosses her viability threshold; the returned
generation count is the RLS. Because production is autocatalytic, `D` accelerates and the
lifespan is finite.

The RLS distribution is set mostly by **cell-to-cell heterogeneity** (`crit_cv`, a lognormal
spread in the threshold `D_crit` across cells) rather than the small per-division noise `cv`:
the autocatalytic blow-up synchronizes the threshold crossing, so per-division noise alone gives
an unrealistically tight distribution. With the default parameters the mean and spread calibrate
to the measured budding-yeast RLS (mean ≈ 25 divisions, CV ≈ 0.3; Schnitzer 2022) — the
threshold/heterogeneity are illustrative, chosen to reproduce that target, not fit per-cell.
"""
function replicative_lifespan(;
    D_crit::Real=32.0,
    crit_cv::Real=0.55,
    production::Real=1.0,
    kappa::Real=0.10,
    cv::Real=0.05,
    alpha0::Real=0.32,
    alpha_max::Real=0.5,
    tau::Real=10.0,
    segregate::Bool=true,
    seed::Int=1,
    max_gen::Int=500,
)
    rng = Random.MersenneTwister(seed)
    # cell-to-cell heterogeneity: this cell's own viability threshold (lognormal, mean D_crit)
    Dc = crit_cv > 0 ? D_crit * exp(crit_cv * randn(rng) - crit_cv^2 / 2) : float(D_crit)
    D = 0.0
    a = 0
    while D < Dc && a < max_gen
        frac = aging_daughter_fraction(a; alpha0=alpha0, alpha_max=alpha_max, tau=tau)
        kept = segregate ? (1 - frac) : 1.0     # the share the mother does NOT pass to the bud
        noise = cv > 0 ? (1 + cv * randn(rng)) : 1.0
        D += kept * production * (1 + kappa * D) * max(0.0, noise)   # autocatalytic damage
        a += 1
    end
    return a
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
# Energetic single-cell growth + the two-step G1 (Di Talia 2007). With a sizer
# threshold V*, the mother/daughter G1 durations EMERGE rather than being imposed.
# ---------------------------------------------------------------------------
# Surface-area-limited QSS growth (the course VOL_Growth law; Altenburg et al. 2019
# constants), per minute: dV/dt = (k_up·4πr² − k_cons·V), r = (3V/4π)^(1/3).
const _K_UP, _K_CONS, _C_ISS = 0.23, 0.27, 319.4
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
"""
function classify_control(slope::Real; atol::Real=0.35)
    isapprox(slope, 0; atol=atol) && return :sizer
    isapprox(slope, 1; atol=atol) && return :adder
    isapprox(slope, 2; atol=atol) && return :timer
    return :mixed
end

end # module
