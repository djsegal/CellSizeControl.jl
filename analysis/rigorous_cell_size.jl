#!/usr/bin/env julia
# Rigorous, literature-faithful budding-yeast cell-size-control lineage model (Track A,
# firewall-clean: generic textbook cell-cycle biology, NOT our TDE/WCM IP). Built to take
# the class cell-size question to a scientifically complete state, validated against the
# published phenomenology that the course whole-cell sim cannot resolve (it follows one
# mother on a fixed clock).
#
# Physics:
#   - Growth: the course VOL_Growth QSS law, PER CELL on its own geometry (the validated
#     per-compartment fix): dV/dt = (k_up*G(V) - k_cons*V)/c_i_ss, G = 4*pi*r^2.
#   - G1 = Whi5-dilution SIZER (Schmoller 2015: Start at the critical size V* = W0/thresh)
#     followed by a fixed Cln2 TIMER step (Di Talia 2007, two-module G1). A cell born small
#     spends real time growing to V* (long G1); a cell already >= V* starts on the timer
#     alone (short G1). => mother short-G1, daughter long-G1 EMERGES, not imposed.
#   - Budded (S/G2/M) = a size-invariant TIMER tau_bud (Allard 2018).
#
# Step 1 (this file): validate the CORE against Di Talia 2007 — mother G1 ~= 19 min,
# daughter G1 ~= 45 min — and the Soifer-Amir size-control slope (sizer -> ~0 in small
# daughters, timer -> >0 in large mothers). Maternal-age asymmetry (CS-DA) + RLS come next.

using Pkg: Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)   # the CellSizeControl package
using CellSizeControl

# Use the package's CANONICAL model — no local duplication. The energetic growth laws
# (qss_growth_rate / exponential_growth_rate), grow_to / grow_for, the two-step-G1
# `cell_cycle`, and the `lineage_timecourse` all live in CellSizeControl now; the aliases
# below keep the maternal-age lineage + validation code readable.
const dVdt = qss_growth_rate
const exp_rate = exponential_growth_rate
const slope = size_control_slope

# CS-DA: a mother lineage with maternal-age asymmetry erosion. The Start size V*(a) rises
# (the mother enlarges with replicative age) and the division asymmetry erodes (beta(a) up,
# so the daughter takes a growing slice of the mother's body, toward symmetric division
# near death; Cdc42 polarity loss). ONE age function then makes late
# daughters LARGER (Johnston 1966; Yang 2011) -- the size face of the same asymmetry erosion as
# AGE-3's damage face.
function mother_lineage(;
    n_max=29,                         # 30 generations = the budding-yeast replicative lifespan (~25-30; Schnitzer 2022)
    Vstar0=36.0,
    T_cln2=19.0,
    tau_bud=52.0,  # n_max+1 gens ≈ yeast RLS (~25-30 divisions; Schnitzer 2022)
    m_enlarge=0.45,
    m_tau=8.0,        # V*(a) = Vstar0*(1 + m_enlarge*(1-e^{-a/m_tau})): the mother enlarges
    r0=0.69,
    r_max=0.90,
    r_tau=14.0,  # daughter/mother size ratio rises with age (asymmetry erosion)
    damage_form=1.0,                  # damage formed per cycle (arb. units)
    dmg_slow=12.0,                    # cycle lengthening per unit accumulated damage (min): calibrated
    # so the cycle reaches ~5.2x its young value by end of life (Egilmez & Jazwinski 1989 ~5-6x;
    # young ~83 min ~ Fehrmann/Charvin 2013 78.3 min)
    phantom_founder=false,            # dormant opt-in (default OFF; see notes below)
    rate=dVdt,
)
    Vstar(a) = Vstar0 * (1.0 + m_enlarge * (1.0 - exp(-a / m_tau)))
    ratio(a) = r0 + (r_max - r0) * (1.0 - exp(-a / r_tau))
    Vm = Vstar0
    Dm = 0.0                                           # mother damage pool (accumulates with age)
    gens = NamedTuple[]
    if phantom_founder
        push!(
            gens,
            (;
                gen=0,
                Vdaughter=Vstar0 * r0,
                Vmother=Vstar0,
                Ddaughter=0.0,
                Dmother=0.0,
                G1=NaN,
                cycle=NaN,
                phantom=true,
            ),
        )
    end
    for a in 0:n_max
        Vs = Vstar(a)
        # The mother is born at/near V* every cycle and is timer-dominated (Di Talia 2007):
        # her G1 is essentially the fixed Cln2 timer. (We deliberately do NOT add the tiny
        # "catch up to the slowly-rising V*" sizer step here: it is a per-generation
        # discretization artifact that puts a spurious knee at gen 2, because only the founder
        # starts exactly at V*(0) with a zero sizer step.)
        G1 = T_cln2
        # CORRECT division accounting: the mother KEEPS her cell body (monotonic, never
        # shrinks); only the BUD pinches off as the daughter. "Old mothers -> larger
        # daughters" (Johnston 1966; Yang 2011) comes from the bud growing as a rising fraction ratio(a)
        # of the ENLARGING mother (a bigger mother feeds a bigger bud => division gets less
        # asymmetric), NOT from the mother giving up a slice of herself.
        Vdaughter = Vs * ratio(a)                      # the bud (daughter); mother stays at Vs
        Vm = Vs                                         # mother kept her body -> next cycle's start
        Dm += damage_form                              # mother accrues damage with replicative age
        # the SAME age-eroding asymmetry that grows the daughter also sets how much of the
        # mother's accrued damage she inherits: young -> rejuvenated, old -> larger (Johnston
        # 1966; Yang 2011) AND more damaged (Kennedy 1994: old-mother daughters reduced
        # lifespan). One mechanism, two faces.
        Ddaughter = (ratio(a) / r_max) * Dm
        cycle = G1 + tau_bud + dmg_slow * Dm           # timer G1 + budded + damage-driven slowing
        push!(
            gens,
            (;
                gen=a + 1,
                Vdaughter,
                Vmother=Vs,
                Ddaughter,
                Dmother=Dm,
                G1,
                cycle,
                phantom=false,
            ),
        )
    end
    return gens
end

# Write the Volume-vs-time lineage trajectory (from the package's `lineage_timecourse`) to a
# CSV in the course growth_*.csv schema so the plotter consumes it identically.
function emit_timecourse(; kwargs...)
    tc = lineage_timecourse(; kwargs...)
    open(joinpath(@__DIR__, "lineage_timecourse.csv"), "w") do io
        println(io, "Time_s,Volume_Mother_L,Volume_Bud_L")
        for i in eachindex(tc.t)
            println(
                io,
                round(tc.t[i] * 60; digits=1),
                ",",
                tc.Vmother[i] * 1e-15,
                ",",
                tc.Vbud[i] * 1e-15,
            )
        end
    end
    println("wrote lineage_timecourse.csv")
    return nothing
end

function main()
    Vstar, T_cln2, tau_bud = 36.0, 19.0, 52.0   # calibrated: mother G1 19, daughter G1 ~45, daughter ~25 fL
    println("=== CORE validation vs Di Talia 2007 (mother 19 / daughter 45 min G1) ===")
    # a representative young-mother daughter birth size (bud of a cell at V*)
    Vdau = cell_cycle(Vstar; Vstar, T_cln2, tau_bud).Vdaughter
    md = cell_cycle(Vstar + 5.0; Vstar, T_cln2, tau_bud)   # a mother (already > V*)
    dd = cell_cycle(Vdau; Vstar, T_cln2, tau_bud)          # a daughter (born small)
    println("  daughter birth size = $(round(Vdau; digits=1)) fL")
    println("  mother   G1 = $(round(md.G1; digits=1)) min   (lit 19.0; mother born >= V*)")
    println(
        "  daughter G1 = $(round(dd.G1; digits=1)) min   (lit 45.5; born small, grows to V*)",
    )
    println(
        "  budded(S/G2/M) = $(round(dd.budded; digits=1)) min ; daughter cycle = $(round(dd.cycle; digits=1)) min",
    )
    println("  division size V_div = $(round(dd.Vdiv; digits=1)) fL")

    # Soifer-Amir slope: V_div vs V_birth. Sizer -> ~0 (small daughters), timer -> >0 (mothers).
    small = collect(range(20.0, Vstar; length=25))         # daughter-range births (<= V*)
    large = collect(range(Vstar, 65.0; length=25))         # mother-range births (>= V*)
    sd = slope(small, [cell_cycle(v; Vstar, T_cln2, tau_bud).Vdiv for v in small])
    sl = slope(large, [cell_cycle(v; Vstar, T_cln2, tau_bud).Vdiv for v in large])
    println("=== Soifer-Amir size-control slope (V_div vs V_birth) ===")
    println(
        "  small/daughter range slope = $(round(sd; digits=2))  (sizer -> ~0; Di Talia daughter sizer)",
    )
    println(
        "  large/mother  range slope = $(round(sl; digits=2))  (timer -> >0; Di Talia mother near-timer)",
    )

    # CS-DA: daughter size vs maternal generation (Johnston 1966; Yang 2011: old mothers -> larger daughters)
    println("=== CS-DA: daughter size vs maternal generation (Johnston 1966; Yang 2011) ===")
    L = mother_lineage()
    real = filter(r -> !r.phantom, L)          # T8: validated generations (1..N), phantom excluded
    ph = filter(r -> r.phantom, L)
    if !isempty(ph)
        println(
            "  gen  0 (PHANTOM founder-bud): daughter = $(round(ph[1].Vdaughter; digits=1)) fL " *
            "-- the founder's OWN birth-as-a-bud (age 0, pristine); makes the lineage uniform (T8)",
        )
    end
    for g in (1, 4, 8, 12, length(real))
        r = real[g]
        println(
            "  gen $(lpad(r.gen, 2)): daughter = $(round(r.Vdaughter; digits=1)) fL, " *
            "mother(Start) = $(round(r.Vmother; digits=1)) fL, " *
            "daughter damage = $(round(r.Ddaughter; digits=2)), cycle = $(round(r.cycle; digits=0)) min",
        )
    end
    d1, dN = real[1].Vdaughter, real[end].Vdaughter
    println(
        "  -> daughter birth size RISES $(round(d1; digits=1)) -> $(round(dN; digits=1)) fL " *
        "across the lineage (1st small/pristine -> late large; asymmetry erosion). " *
        "mother enlarges $(round(L[1].Vmother; digits=1)) -> $(round(L[end].Vmother; digits=1)) fL.",
    )
    # T7: the SAME beta(a) makes daughter DAMAGE rise too (the fitness face; Kennedy 1994)
    println(
        "  -> T7 UNIFY: daughter inherited DAMAGE RISES $(round(L[1].Ddaughter; digits=2)) -> " *
        "$(round(L[end].Ddaughter; digits=2)) over the lineage (old-mother daughters shorter-lived) " *
        "-- driven by the SAME asymmetry erosion beta(a) as the size rise. One function, two faces (size + fitness).",
    )

    # T3: exponential (biomass-driven) growth -- the data-matching law (Di Talia 2007; Sun
    # 2010). Under dV/dt=mu*V the G1 sizer is EXACTLY Di Talia's ln(V*/Vb)/mu + T_cln2, and
    # the bud = the mother's exponential growth over the budded timer, V_start*(e^{mu*tau}-1),
    # which grows WITH cell size -- so daughters self-stabilize with no surface-area
    # starvation (exponential growth fixes the collapse without the per-compartment trick).
    println(
        "=== T3: exponential growth mu=0.0077/min (~90-min doubling) -- Di Talia-consistent ===",
    )
    mu, Vse, taub = 0.0077, 30.0, 78.0
    g1e(Vb) = (Vb < Vse ? log(Vse / Vb) / mu : 0.0) + T_cln2
    Vd0 = Vse * (exp(mu * taub) - 1.0)                 # daughter = budded-phase exp growth
    vdive(Vb) = max(Vb, Vse) * exp(mu * taub)
    se = slope(collect(20.0:0.5:Vse), [vdive(v) for v in 20.0:0.5:Vse])
    le = slope(collect(Vse:0.5:65.0), [vdive(v) for v in Vse:0.5:65.0])
    println(
        "  mother G1 = $(round(g1e(Vse + 5); digits=1)) min (lit 19); " *
        "daughter (born $(round(Vd0; digits=1)) fL) G1 = $(round(g1e(Vd0); digits=1)) min (lit 45)",
    )
    println(
        "  slope: daughters $(round(se; digits=2)) (sizer->0), mothers $(round(le; digits=2)) " *
        "(near-timer; e^{mu*tau}=$(round(exp(mu * taub); digits=2)))",
    )
    println(
        "  daughters STABLE at $(round(Vd0; digits=1)) fL (= V*(e^{mu*tau}-1), constant): " *
        "exponential growth self-stabilizes -- no collapse, no per-compartment trick needed.",
    )
    open(joinpath(@__DIR__, "cs_da_lineage.csv"), "w") do io
        println(io, "gen,Vdaughter,Vmother,G1,cycle,Ddaughter")
        for r in L
            println(
                io,
                r.gen,
                ",",
                round(r.Vdaughter; digits=3),
                ",",
                round(r.Vmother; digits=3),
                ",",
                round(r.G1; digits=3),
                ",",
                round(r.cycle; digits=3),
                ",",
                round(r.Ddaughter; digits=3),
            )
        end
    end
    println("wrote cs_da_lineage.csv")
    emit_timecourse()
    return nothing
end

main()
