#!/usr/bin/env julia
# CC-P (population campaign): grow a synchronous exponential culture to ~10^6 cells and record
# its STEADY-STATE structure — the replicative-age distribution (which converges to the geometric
# law P(age=a)=2^{-(a+1)}, half the cells virgin daughters; Hartwell & Unger 1977) and the
# newborn (age-0) birth-size distribution (right-skewed by the rare, enlarged old mothers that
# bud the biggest daughters). Writes population_age_structure.csv + population_newborn_size.csv +
# population_summary.csv. Run: julia --project=.. gen_population.jl [TARGET]  (default 1_000_000)

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)
using CellSizeControl
using Statistics: mean, std

const HERE = @__DIR__

function main()
    target = isempty(ARGS) ? 1_000_000 : parse(Int, ARGS[1])
    Vstar, frac0 = 60.0, 0.32
    pop = simulate_population(
        SizerRule(Vstar);
        target=target,
        enlarge_max=0.45,
        enlarge_tau=8.0,
        alpha0=frac0,
        alpha_max=0.5,
        tau=8.0,
        seed=1,
    )
    N = length(pop.age)
    amax = maximum(pop.age)

    # ---- (a) replicative-age structure vs the geometric prediction 2^{-(a+1)} ----
    counts = zeros(Int, amax + 1)
    for a in pop.age
        counts[a + 1] += 1
    end
    open(joinpath(HERE, "population_age_structure.csv"), "w") do io
        println(io, "age,count,fraction,geometric")
        for a in 0:amax
            frac = counts[a + 1] / N
            geo = 2.0^(-(a + 1))
            println(
                io,
                a,
                ",",
                counts[a + 1],
                ",",
                round(frac; digits=8),
                ",",
                round(geo; digits=8),
            )
        end
    end

    # ---- (b) newborn (age-0) birth-size distribution ----
    nb = pop.Vbirth[pop.age .== 0]
    lo, hi = extrema(nb)
    nbins = 60
    edges = range(lo, hi; length=nbins + 1)
    h = zeros(Int, nbins)
    for v in nb
        b = clamp(floor(Int, (v - lo) / (hi - lo) * nbins) + 1, 1, nbins)
        h[b] += 1
    end
    open(joinpath(HERE, "population_newborn_size.csv"), "w") do io
        println(io, "bin_center,count")
        for b in 1:nbins
            c = 0.5 * (edges[b] + edges[b + 1])
            println(io, round(c; digits=4), ",", h[b])
        end
    end

    # ---- (c) summary statistics (the gated headline numbers) ----
    senescent = count(i -> pop.age[i] >= pop.rls[i], 1:N)
    open(joinpath(HERE, "population_summary.csv"), "w") do io
        println(io, "metric,value")
        println(io, "N_cells,", N)
        println(io, "generations,", pop.ngen)
        println(io, "virgin_fraction,", round(counts[1] / N; digits=6))
        println(io, "age1_fraction,", round(counts[2] / N; digits=6))
        println(io, "ratio_age1_age0,", round(counts[2] / counts[1]; digits=6))
        println(io, "mean_replicative_age,", round(mean(pop.age); digits=6))
        println(io, "max_replicative_age,", amax)
        println(io, "senescent_fraction,", senescent / N)
        println(io, "newborn_size_mean,", round(mean(nb); digits=4))
        println(io, "newborn_size_sd,", round(std(nb); digits=4))
        println(io, "newborn_size_min,", round(lo; digits=4))
        println(io, "newborn_size_max,", round(hi; digits=4))
    end

    println("CC-P done: N=$N cells over $(pop.ngen) generations")
    println(
        "  virgin frac=$(round(counts[1]/N;digits=4)) (geom 0.5); ",
        "age1/age0=$(round(counts[2]/counts[1];digits=4)) (geom 0.5); ",
        "mean age=$(round(mean(pop.age);digits=4)) (geom 1.0)",
    )
    println(
        "  senescent frac=$(round(senescent/N;sigdigits=3)); ",
        "newborn size mean=$(round(mean(nb);digits=2)) sd=$(round(std(nb);digits=2))",
    )
    println(
        "wrote population_age_structure.csv + population_newborn_size.csv + population_summary.csv",
    )
    return nothing
end

main()
