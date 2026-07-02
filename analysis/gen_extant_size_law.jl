#!/usr/bin/env julia
# CC-X (extant-vs-newborn size divergence + the senescence age-law correction).
#
# A snapshot of a balanced exponentially-growing culture samples every cell at its most-recent
# division: age-0 cells are the small buds (the newborn-size law), age a≥1 cells are mothers
# carrying their full retained body V*·enlarge(a−1). Weighting by the geometric replicative-age
# law 2^{-(a+1)} makes the MEAN EXTANT cell far larger than the mean newborn — the size-structure
# signature of exponential growth, where the standing population over-represents the larger, older
# mother bodies relative to the small buds they shed. CellSizeControl.extant_size_law returns the
# closed-form divergence D = extant_mean/newborn_mean; it is scale-free (V*-independent) and, with
# no erosion/enlargement, reduces to the exact (1+α0)/(2α0).
#
# The geometric age law itself carries a senescence correction: at a finite lifespan the dividing
# population's age distribution is the truncated geometric λ^{-(a+1)}, λ<2 solving the discrete
# Euler–Lotka λ=Σ_{a=0}^{rls−1}λ^{-a} (senescence_age_law); λ→2 recovers 2^{-(a+1)}, rls=2 gives φ.
#
# Writes extant_size_law.csv, extant_age_decomp.csv, extant_size_hist.csv, senescence_age_law.csv.
# Run: julia --project=.. gen_extant_size_law.jl [TARGET]

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)
using CellSizeControl
using Statistics: mean, std

const HERE = @__DIR__
const A0, AMAX, TAU = 0.32, 0.5, 8.0
const EM, ET = 0.45, 8.0

function popmeans(Vstar; target, seed=1, kwargs...)
    pop = simulate_population(
        SizerRule(Vstar); target=target, enlarge_max=EM, enlarge_tau=ET,
        alpha0=A0, alpha_max=AMAX, tau=TAU, seed=seed, kwargs...,
    )
    ext = mean(pop.Vbirth)
    nb = mean(pop.Vbirth[pop.age .== 0])
    return (pop=pop, ext=ext, nb=nb, D=ext / nb)
end

function main()
    target = isempty(ARGS) ? 300_000 : parse(Int, ARGS[1])
    Vstar = 60.0

    law = extant_size_law(; alpha0=A0, alpha_max=AMAX, tau=TAU, enlarge_max=EM,
        enlarge_tau=ET, Vstar=Vstar)
    # deterministic MC isolates the demographic law (no per-cell noise / senescence tail)
    sim = popmeans(Vstar; target=target, seed=1, cv=0.0, crit_cv=0.0)

    # ---- (1) analytic vs simulation + scale-free across set-points + dependence ----
    open(joinpath(HERE, "extant_size_law.csv"), "w") do io
        println(io, "case,Vstar,source,newborn_mean,extant_mean,divergence")
        println(io, "calibrated,", Vstar, ",analytic,",
            round(law.newborn_mean; digits=4), ",", round(law.extant_mean; digits=4), ",",
            round(law.divergence; digits=5))
        println(io, "calibrated,", Vstar, ",simulation,",
            round(sim.nb; digits=4), ",", round(sim.ext; digits=4), ",",
            round(sim.D; digits=5))
        for V in (30.0, 120.0, 240.0)
            l = extant_size_law(; alpha0=A0, alpha_max=AMAX, tau=TAU, enlarge_max=EM,
                enlarge_tau=ET, Vstar=V)
            s = popmeans(V; target=target, seed=2, cv=0.0, crit_cv=0.0)
            println(io, "scalefree,", V, ",analytic,",
                round(l.newborn_mean; digits=4), ",", round(l.extant_mean; digits=4), ",",
                round(l.divergence; digits=5))
            println(io, "scalefree,", V, ",simulation,",
                round(s.nb; digits=4), ",", round(s.ext; digits=4), ",",
                round(s.D; digits=5))
        end
        # dependence: how the demographic parameters move D (V* never does)
        deps = (
            ("full", EM, AMAX),
            ("erosion_only", 0.0, AMAX),
            ("enlarge_only", EM, A0),
            ("no_erosion_analytic", 0.0, A0),   # D reduces to (1+α0)/(2α0) = 2.0625
        )
        for (name, em, amax) in deps
            l = extant_size_law(; alpha0=A0, alpha_max=amax, tau=TAU, enlarge_max=em,
                enlarge_tau=ET, Vstar=Vstar)
            println(io, name, ",", Vstar, ",analytic,",
                round(l.newborn_mean; digits=4), ",", round(l.extant_mean; digits=4), ",",
                round(l.divergence; digits=5))
        end
    end

    # ---- (2) age decomposition of the extant mean: newborn bud vs retained mother body ----
    open(joinpath(HERE, "extant_age_decomp.csv"), "w") do io
        println(io, "age,weight,newborn_size,mother_body")
        for a in 0:20
            w = 2.0^(-(a + 1))
            frac = aging_daughter_fraction(a; alpha0=A0, alpha_max=AMAX, tau=TAU)
            nbsize = frac * Vstar * (1 + EM * (1 - exp(-a / ET)))
            body = a == 0 ? NaN : Vstar * (1 + EM * (1 - exp(-(a - 1) / ET)))
            println(io, a, ",", round(w; digits=8), ",", round(nbsize; digits=4), ",",
                a == 0 ? "" : string(round(body; digits=4)))
        end
    end

    # ---- (3) histograms (rescaled by V*) of the extant vs newborn populations ----
    open(joinpath(HERE, "extant_size_hist.csv"), "w") do io
        println(io, "population,size_over_Vstar,fraction")
        allv = sim.pop.Vbirth ./ Vstar
        nbv = sim.pop.Vbirth[sim.pop.age .== 0] ./ Vstar
        lo, hi = 0.0, maximum(allv) * 1.02
        nbins = 70
        for (tag, data) in (("extant", allv), ("newborn", nbv))
            h = zeros(Int, nbins)
            for v in data
                b = clamp(floor(Int, (v - lo) / (hi - lo) * nbins) + 1, 1, nbins)
                h[b] += 1
            end
            tot = sum(h)
            for b in 1:nbins
                c = lo + (b - 0.5) * (hi - lo) / nbins
                println(io, tag, ",", round(c; digits=5), ",", round(h[b] / tot; digits=6))
            end
        end
    end

    # ---- (4) senescence correction: Euler–Lotka λ vs mean lifespan, + a short-RLS age law ----
    open(joinpath(HERE, "senescence_age_law.csv"), "w") do io
        println(io, "rls,lambda,virgin_fraction")
        for m in 2:40
            sl = senescence_age_law(m)
            println(io, m, ",", round(sl.lambda; digits=8), ",", round(1 / sl.lambda; digits=6))
        end
    end
    # per-age departure at a short deterministic lifespan (analytic λ^{-(a+1)} + naive 2^{-(a+1)}
    # + Monte-Carlo dividing-cell fractions) — the falsifiable senescence signature
    Dc_short = 3.1
    m_short = replicative_lifespan(; D_crit=Dc_short, cv=0.0, crit_cv=0.0)
    sl = senescence_age_law(m_short)
    sp = simulate_population(SizerRule(Vstar); target=target, D_crit=Dc_short, cv=0.0,
        crit_cv=0.0, max_gen=400, seed=1)
    divmask = sp.age .< sp.rls
    Nd = count(divmask)
    open(joinpath(HERE, "senescence_short_rls.csv"), "w") do io
        println(io, "age,lotka,naive_geometric,simulation")
        for a in 0:(m_short - 1)
            emp = count(i -> divmask[i] && sp.age[i] == a, eachindex(sp.age)) / Nd
            println(io, a, ",", round(sl.p[a + 1]; digits=6), ",",
                round(2.0^(-(a + 1)); digits=6), ",", round(emp; digits=6))
        end
    end

    println("CC-X done (target=$target)")
    println("  extant/newborn divergence  analytic D=$(round(law.divergence;digits=5)) ",
        "simulation D=$(round(sim.D;digits=5))  (extant $(round(law.extant_mean;digits=2)) fL vs ",
        "newborn $(round(law.newborn_mean;digits=2)) fL at V*=$Vstar)")
    println("  senescence  short rls m=$m_short  λ=$(round(sl.lambda;digits=5)) ",
        "(virgin frac $(round(1/sl.lambda;digits=4)) vs geometric 0.5); λ(2)=φ, λ→2 as rls→∞")
    println("wrote extant_size_law.csv, extant_age_decomp.csv, extant_size_hist.csv, ",
        "senescence_age_law.csv, senescence_short_rls.csv")
    return nothing
end

main()
