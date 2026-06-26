#!/usr/bin/env julia
# CC-1: the emergent-RLS parameter landscape. Sweep the autocatalytic-damage params
# (viability threshold D_crit × autocatalysis kappa) and, at each grid point, the
# Monte-Carlo mean and CV of the replicative lifespan (CellSizeControl.lifespan_distribution).
# Tests whether the Schnitzer-2022 calibration (mean ~25, CV ~0.3) sits on a fine-tuned point
# or a broad basin. Run: julia --project=.. gen_rls_landscape.jl  (writes rls_landscape.csv)

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)
using CellSizeControl
using Statistics: mean, std

const HERE = @__DIR__
const N = 1500                                   # lineages per grid point
const DCRIT = collect(20.0:2.0:48.0)             # viability threshold
const KAPPA = collect(0.04:0.01:0.18)            # autocatalysis strength

open(joinpath(HERE, "rls_landscape.csv"), "w") do io
    println(io, "D_crit,kappa,mean_rls,cv_rls")
    for Dc in DCRIT, k in KAPPA
        s = lifespan_distribution(N; D_crit=Dc, kappa=k)
        m = mean(s)
        c = std(s) / m
        println(io, Dc, ",", k, ",", round(m; digits=3), ",", round(c; digits=4))
    end
end

# report the Schnitzer-matching region (mean 23-27 AND CV 0.25-0.35)
hits = 0
total = length(DCRIT) * length(KAPPA)
for Dc in DCRIT, k in KAPPA
    s = lifespan_distribution(N; D_crit=Dc, kappa=k)
    m = mean(s)
    c = std(s) / m
    (23 <= m <= 27 && 0.25 <= c <= 0.35) && (global hits += 1)
end
println("wrote rls_landscape.csv ($(length(DCRIT))x$(length(KAPPA)) grid, N=$N)")
println("Schnitzer-matching grid points (mean 23-27 & CV 0.25-0.35): $hits / $total")
