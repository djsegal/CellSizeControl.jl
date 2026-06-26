#!/usr/bin/env julia
# CC-2: large stochastic lineage ensemble. Many noisy aging lineages
# (simulate_aging_lineage, cv>0) → the joint distributions of daughter birth SIZE and
# inherited DAMAGE, and the mother(at-division)→daughter(birth) correlation structure (the
# Soifer-Amir-style inheritance plot, here under the maternal-age asymmetry). Writes pooled
# summary stats + a subsample of the (Vmother_div, Vdaughter) pairs for plotting.
# Run: julia --project=.. gen_lineage_ensemble.jl [N]   (N lineages; default 200_000)

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)
using CellSizeControl
using Random: MersenneTwister

const HERE = @__DIR__
const NGEN = 30
const SUBSAMPLE = 6000          # pairs written for the scatter

function run_ensemble(N::Int)
    xm = Float64[]              # mother volume at division (subsample)
    xd = Float64[]              # daughter birth volume (subsample)
    xdmg = Float64[]            # daughter inherited damage (subsample)
    sum_v = 0.0;
    sumsq_v = 0.0  # daughter birth size moments
    sum_d = 0.0;
    sumsq_d = 0.0  # inherited damage moments
    sxy = 0.0;
    sx = 0.0;
    sy = 0.0;
    sxx = 0.0;
    syy = 0.0
    npair = 0
    # RANDOM (unbiased) subsample for the scatter — NOT a deterministic every-k stride, which
    # would be periodic in generation and band the plot. Keep each pair with prob p.
    rng = MersenneTwister(20260626)
    p = SUBSAMPLE / (N * NGEN)
    for i in 1:N
        L = simulate_aging_lineage(
            SizerRule(60.0);
            n=NGEN,
            alpha0=0.32,
            alpha_max=0.5,
            tau=8.0,
            enlarge_max=0.45,
            enlarge_tau=8.0,
            cv=0.10,
            seed=i,
        )
        for a in 1:NGEN
            vm = L.Vdivision[a]    # mother body at division
            vd = L.Vdaughter[a]    # her bud / daughter birth size
            dd = L.Ddaughter[a]
            sum_v += vd;
            sumsq_v += vd^2
            sum_d += dd;
            sumsq_d += dd^2
            sx += vm;
            sy += vd;
            sxx += vm^2;
            syy += vd^2;
            sxy += vm * vd
            npair += 1
            if length(xm) < SUBSAMPLE && rand(rng) < p
                push!(xm, vm)
                push!(xd, vd)
                push!(xdmg, dd)
            end
        end
    end
    mv = sum_v / npair
    sdv = sqrt(max(0.0, sumsq_v / npair - mv^2))
    md = sum_d / npair
    sdd = sqrt(max(0.0, sumsq_d / npair - md^2))
    r = (npair * sxy - sx * sy) / sqrt((npair * sxx - sx^2) * (npair * syy - sy^2))
    return (; N, npair, mv, sdv, md, sdd, r, xm, xd, xdmg)
end

function main()
    N = isempty(ARGS) ? 200_000 : parse(Int, ARGS[1])
    e = run_ensemble(N)
    open(joinpath(HERE, "lineage_ensemble_summary.csv"), "w") do io
        println(io, "metric,value")
        println(io, "N_lineages,", e.N)
        println(io, "n_pairs,", e.npair)
        println(io, "daughter_size_mean,", round(e.mv; digits=4))
        println(io, "daughter_size_cv,", round(e.sdv / e.mv; digits=4))
        println(io, "daughter_damage_mean,", round(e.md; digits=4))
        println(io, "daughter_damage_cv,", round(e.sdd / e.md; digits=4))
        return println(io, "mother_daughter_size_corr,", round(e.r; digits=4))
    end
    open(joinpath(HERE, "lineage_ensemble_pairs.csv"), "w") do io
        println(io, "Vmother_div,Vdaughter,Ddaughter")
        for (a, b, c) in zip(e.xm, e.xd, e.xdmg)
            println(
                io, round(a; digits=4), ",", round(b; digits=4), ",", round(c; digits=4)
            )
        end
    end
    println("CC-2 done: N=$(e.N) lineages, $(e.npair) pairs")
    println(
        "  daughter size mean=$(round(e.mv;digits=2)) CV=$(round(e.sdv/e.mv;digits=3)); ",
        "mother-daughter size corr r=$(round(e.r;digits=3))",
    )
    println("wrote lineage_ensemble_summary.csv + lineage_ensemble_pairs.csv")
    return nothing
end

main()
