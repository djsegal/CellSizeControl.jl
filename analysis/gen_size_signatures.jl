# Gallery: four dynamical signatures of size control, computed from the linear map Vd=αVb+β
# at matched mean birth size (Vb* = 20 for all α, via β = 40(1−α/2) at f=0.5). Sizer (α=0),
# adder (α=1), weak timer (α=1.5, still homeostatic since α·f=0.75<1). Emits: birth-size
# distributions; the consecutive birth-size inheritance map; the step-response relaxation; and
# the birth-size autocorrelation (memory). Run: julia --project=. .
using CellSizeControl
using Statistics: cor, mean
using Printf

here = @__DIR__
regimes = [(name="sizer", α=0.0), (name="adder", α=1.0), (name="timer", α=1.5)]
β(α) = 40.0 * (1 - α / 2)        # keeps the fixed-point birth size at 20 for every α (f=0.5)

# ---- (a) birth-size distributions + (d) autocorrelation, from one long lineage each ----
open(joinpath(here, "sig_hist.csv"), "w") do io
    println(io, "regime,Vb")
    for r in regimes
        s = simulate_lineage(LinearSizeControl(r.α, β(r.α)); V0=20.0, n=6000, cv=0.12, seed=3)
        for v in s.Vb[1001:end]               # drop transient
            @printf(io, "%s,%.4f\n", r.name, v)
        end
    end
end

open(joinpath(here, "sig_acf.csv"), "w") do io
    println(io, "regime,lag,acf")
    for r in regimes
        s = simulate_lineage(LinearSizeControl(r.α, β(r.α)); V0=20.0, n=20000, cv=0.12, seed=5)
        x = s.Vb[1001:end]
        for k in 0:8
            a = k == 0 ? 1.0 : cor(@view(x[1:(end - k)]), @view(x[(1 + k):end]))
            @printf(io, "%s,%d,%.5f\n", r.name, k, a)
        end
    end
end

# ---- (b) inheritance map: consecutive birth sizes Vb_{n+1} vs Vb_n (slope = α·f) ----
open(joinpath(here, "sig_pairs.csv"), "w") do io
    println(io, "regime,Vb_n,Vb_next")
    for r in regimes
        s = simulate_lineage(LinearSizeControl(r.α, β(r.α)); V0=20.0, n=2200, cv=0.12, seed=7)
        x = s.Vb[1001:end]
        for i in 1:(length(x) - 1)
            @printf(io, "%s,%.4f,%.4f\n", r.name, x[i], x[i + 1])
        end
    end
end

# ---- (c) step-response: relaxation of mean birth size after a 2x perturbation ----
open(joinpath(here, "sig_step.csv"), "w") do io
    println(io, "regime,gen,meanVb")
    for r in regimes
        reps = [simulate_lineage(LinearSizeControl(r.α, β(r.α)); V0=40.0, n=14, cv=0.12,
                                 seed=100 + j).Vb for j in 1:400]
        for g in 1:14
            @printf(io, "%s,%d,%.4f\n", r.name, g, mean(rep[g] for rep in reps))
        end
    end
end

println("wrote sig_hist.csv, sig_acf.csv, sig_pairs.csv, sig_step.csv")
