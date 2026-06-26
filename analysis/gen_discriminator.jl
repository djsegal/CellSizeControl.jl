#!/usr/bin/env julia
# Generate the size-control discriminator data from the CellSizeControl package: the
# Vb -> Vd slope for a timer / adder / sizer (Soifer-Amir), and a sub-doubling-timer-
# collapses vs inhibitor-dilution-sizer-stable lineage. Run: julia gen_discriminator.jl
using Pkg: Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)
using CellSizeControl
using Printf

const OUT = @__DIR__

open(joinpath(OUT, "discriminator.csv"), "w") do io
    println(io, "rule,Vb,Vd")
    for (nm, r) in
        (("timer", TimerRule(2.0)), ("adder", AdderRule(1.0)), ("sizer", SizerRule(2.0)))
        s = simulate_lineage(r; n=400, seed=1)
        for (b, d) in zip(s.Vb, s.Vd)
            @printf(io, "%s,%.5f,%.5f\n", nm, b, d)
        end
    end
end

open(joinpath(OUT, "collapse.csv"), "w") do io
    println(io, "gen,timer_Vb,sizer_Vb")
    t = simulate_lineage(TimerRule(1.6); n=40, cv=0.0, seed=5)
    z = simulate_lineage(InhibitorDilutionSizer(60.0, 1.5); n=40, cv=0.0, seed=6)
    for i in 1:40
        @printf(io, "%d,%.5f,%.5f\n", i, t.Vb[i], z.Vb[i])
    end
end

println("wrote discriminator.csv, collapse.csv")
