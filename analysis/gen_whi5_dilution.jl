#!/usr/bin/env julia
# Gallery: the inhibitor-dilution sizer mechanism. Whi5 is made in a fixed dose W per cell;
# its concentration [Whi5] = W/V DILUTES as the cell grows, and Start fires when it crosses
# the threshold θ — i.e. at the critical size V* = W/θ (Schmoller 2015). A daughter (born
# small) must dilute Whi5 over a long G1; a mother (born ≥ V*) is already past threshold and
# fires almost immediately — the Di Talia mother/daughter G1 asymmetry, mechanistically.
# Integrates dV/dt = qss_growth_rate(V) recording (t, V, [Whi5]). Run: julia --project=.. gen_whi5_dilution.jl
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)
using CellSizeControl

const HERE = @__DIR__
const W, THRESH = 60.0, 1.5          # InhibitorDilutionSizer(W, θ) → V* = W/θ = 40
const VSTAR = W / THRESH

function dilute_to_start(Vb; dt=0.1, tmax=400.0)
    t, V = 0.0, float(Vb)
    ts, vs, cs = [0.0], [V], [W / V]
    while W / V > THRESH && t < tmax          # grow until [Whi5] diluted to threshold (Start)
        V += max(0.0, qss_growth_rate(V)) * dt
        t += dt
        push!(ts, t)
        push!(vs, V)
        push!(cs, W / V)
    end
    return ts, vs, cs
end

function main()
    open(joinpath(HERE, "whi5_dilution.csv"), "w") do io
        println(io, "cell,t,V,whi5_conc")
        # daughter born at Vb=30 dilutes to V*=40 in ~22 min (the sizer step; the
        # 19-min CLN2 timer follows → ~41 min total daughter G1). mother born at Vb=41
        # is already past V*, so ~0 sizer step and her G1 is the 19-min timer only.
        for (label, Vb) in (("daughter", 30.0), ("mother", 41.0))
            ts, vs, cs = dilute_to_start(Vb)
            for (t, v, c) in zip(ts, vs, cs)
                println(
                    io,
                    label,
                    ",",
                    round(t; digits=2),
                    ",",
                    round(v; digits=3),
                    ",",
                    round(c; digits=4),
                )
            end
            println(
                stderr,
                "  $label born V=$Vb → Start at t=$(round(ts[end];digits=1)) min, ",
                "V=$(round(vs[end];digits=1))",
            )
        end
    end
    println(
        "V* = ",
        VSTAR,
        " (threshold θ = ",
        THRESH,
        ", W = ",
        W,
        "); wrote whi5_dilution.csv",
    )
    return nothing
end

main()
