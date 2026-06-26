#!/usr/bin/env julia
# AGE-2 emergent-replicative-lifespan data for the figure: (1) a large sample of RLS values
# from the autocatalytic-damage model (CellSizeControl.replicative_lifespan), and (2) a few
# example damage trajectories D(a) showing the acceleration to each cell's viability threshold.
# Run: julia --project=.. gen_lifespan.jl   (writes lifespan_samples.csv + damage_traces.csv)

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)
using CellSizeControl
using Random: MersenneTwister
using Statistics: mean, std

const HERE = @__DIR__

# (1) the RLS distribution (default-calibrated to Schnitzer 2022: mean ~25, CV ~0.3)
samples = lifespan_distribution(5000)
open(joinpath(HERE, "lifespan_samples.csv"), "w") do io
    println(io, "rls")
    for r in samples
        println(io, r)
    end
end
m, sd = mean(samples), std(samples)
println("RLS samples: n=", length(samples), "  mean=", round(m; digits=2),
    "  sd=", round(sd; digits=2), "  CV=", round(sd / m; digits=3))

# (2) example damage trajectories: re-derive D(a) under the same defaults for a few cells,
# so the figure can show the autocatalytic acceleration crossing each cell's threshold.
function damage_trace(; seed, D_crit=32.0, crit_cv=0.55, production=1.0, kappa=0.10,
    cv=0.05, alpha0=0.32, alpha_max=0.5, tau=10.0, max_gen=500)
    rng = MersenneTwister(seed)
    Dc = crit_cv > 0 ? D_crit * exp(crit_cv * randn(rng) - crit_cv^2 / 2) : float(D_crit)
    D = 0.0
    ages = Float64[0.0]
    dvals = Float64[0.0]
    a = 0
    while D < Dc && a < max_gen
        frac = aging_daughter_fraction(a; alpha0=alpha0, alpha_max=alpha_max, tau=tau)
        kept = 1 - frac
        noise = cv > 0 ? (1 + cv * randn(rng)) : 1.0
        D += kept * production * (1 + kappa * D) * max(0.0, noise)
        a += 1
        push!(ages, a)
        push!(dvals, D)
    end
    return ages, dvals, Dc
end

open(joinpath(HERE, "damage_traces.csv"), "w") do io
    println(io, "cell,age,damage,threshold")
    for (i, seed) in enumerate((3, 11, 27, 42, 88))
        ages, dvals, Dc = damage_trace(; seed=seed)
        for (a, d) in zip(ages, dvals)
            println(io, i, ",", a, ",", round(d; digits=4), ",", round(Dc; digits=4))
        end
    end
end
println("wrote lifespan_samples.csv + damage_traces.csv")
