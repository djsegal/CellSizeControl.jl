#!/usr/bin/env julia
# Gallery: replicative-lifespan survival data. A large sample of emergent RLS values
# (CellSizeControl.lifespan_distribution) → the survival function S(a) = fraction of mothers
# still dividing after a divisions (a Kaplan-Meier-style mortality curve, the classic aging
# readout), plus the raw distribution. Run: julia --project=.. gen_rls_survival.jl
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)
using CellSizeControl
using Statistics: mean, median

const HERE = @__DIR__

function main()
    rls = lifespan_distribution(50_000)
    amax = maximum(rls)
    n = length(rls)
    # survival S(a) = P(RLS > a): fraction still dividing after a divisions
    open(joinpath(HERE, "rls_survival.csv"), "w") do io
        println(io, "age,survival,n_dividing")
        for a in 0:amax
            alive = count(>(a), rls)
            println(io, a, ",", round(alive / n; digits=5), ",", alive)
        end
    end
    println(
        "RLS survival: N=$n, mean=$(round(mean(rls);digits=2)), ",
        "median=$(median(rls)), max=$amax",
    )
    println("wrote rls_survival.csv")
    return nothing
end

main()
