# Lineage birth-size MEMORY: the return map is AR(1) with a single pole r = α·f.
#
# The birth-size return map Vb(n+1) = f·(α·Vb(n)+β)·(1+cv·ξ) linearizes about its fixed
# point Vb* = fβ/(1−αf) to an AR(1) process whose one pole r = α·f (α = map_slope, f =
# daughter fraction) sets EVERY memory observable at once:
#
#   • lag-k lineage autocorrelation   ρ_k = r^k             (geometric decay, rate r)
#   • stationary CV amplification      CV(Vb) = cv/√(1−r²)
#   • nutrient-shift relaxation        ⟨Vb(n)⟩ → Vb*' at rate r, memory −1/ln r generations
#
# Two of these force a single-lineage, mode- and set-point-free invariant that recovers the
# intrinsic per-division noise:  CV(Vb)²·(1 − ρ_1²) = cv².
#
# Prediction (mode-diagnostic memory): a SIZER (α=0, r=0) is memoryless — every birth is
# independent (ρ_k≡0) and a nutrient shift is absorbed in a single division; an ADDER (α=1,
# r=f) carries ~1.4 generations of size memory; a TIMER (α=2, r=2f) carries the longest
# memory and diverges as division symmetrizes with age (r=2f→1). Falsification: nonzero
# mother→daughter birth-size correlation in a sizer, a correlation that does not equal the
# return slope αf, a nutrient-shift transient whose relaxation rate is mode-independent, or
# CV²(1−ρ1²) that fails to collapse onto the same intrinsic noise across modes/set-points.
#
# Emits: size_memory_autocorr.csv (measured ρ_k vs r^k across the sizer→timer axis),
# size_memory_invariant.csv (recovered cv = √(CV²(1−ρ1²)) vs true cv, per mode),
# size_memory_step.csv (nutrient-step mean-birth-size relaxation, measured vs geometric).
# Run: julia --project=. analysis/gen_size_memory.jl
using CellSizeControl
using Statistics: mean, std, cor
using Printf

here = @__DIR__

const CV = 0.06      # per-division multiplicative noise
const BETA = 20.0    # set-point scale (β)

fixedpoint(α, β, f) = f * β / (1 - α * f)

# Sample lag-k autocorrelation of the birth-size series (burn-in discarded).
function autocorr(Vb, k; burn=2000)
    x = @view Vb[(burn + 1):end]
    cor(@view(x[1:(end - k)]), @view(x[(k + 1):end]))
end

# ---- (a) lag-k autocorrelation ρ_k vs the closed form r^k, along the sizer→timer axis ----
# One long lineage per α (fully mixed); f = 0.5 (symmetric division).
alphas = (0.0, 0.5, 1.0, 1.5, 1.6)
f = 0.5
open(joinpath(here, "size_memory_autocorr.csv"), "w") do io
    println(io, "alpha,f,r,lag,rho_measured,rho_pred")
    for α in alphas
        r = α * f
        s = simulate_lineage(
            LinearSizeControl(α, BETA); V0=BETA, n=80_000, cv=CV, daughter_fraction=f, seed=7
        )
        for k in 1:6
            @printf(io, "%.3f,%.3f,%.4f,%d,%.6f,%.6f\n", α, f, r, k, autocorr(s.Vb, k), r^k)
        end
    end
end

# ---- (b) the single-lineage invariant  cv_recovered = √(CV(Vb)²(1−ρ1²)) ≈ cv ----
# Mode- and set-point-free: the same intrinsic noise falls out for every (α, f, β).
modes = (("sizer", 0.0), ("adder", 1.0), ("timer", 1.6))
open(joinpath(here, "size_memory_invariant.csv"), "w") do io
    println(io, "mode,alpha,f,beta,r,cv_measured,rho1,cv_recovered,cv_true")
    for (name, α) in modes, β in (20.0, 60.0), fr in (0.4, 0.5)
        r = α * fr
        s = simulate_lineage(
            LinearSizeControl(α, β); V0=β, n=80_000, cv=CV, daughter_fraction=fr, seed=11
        )
        x = @view s.Vb[2001:end]
        CVb = std(x) / mean(x)
        ρ1 = autocorr(s.Vb, 1)
        cvrec = sqrt(CVb^2 * (1 - ρ1^2))
        @printf(
            io, "%s,%.2f,%.2f,%.1f,%.4f,%.6f,%.6f,%.6f,%.6f\n",
            name, α, fr, β, r, CVb, ρ1, cvrec, CV
        )
    end
end

# ---- (c) nutrient-shift step response: β: 20 → 40 (a set-point doubling) ----
# Burn in each of R lineages under β1, step the set-point to β2, then average the birth-size
# trajectory across seeds. The mean relaxes geometrically: ⟨Vb(n)⟩ = Vb*₂ + (Vb*₁−Vb*₂)·r^n.
function step_response(α, fr; β1=20.0, β2=40.0, burn=400, M=20, R=6000)
    traj = zeros(M + 1)
    for s in 1:R
        v0 = last(
            simulate_lineage(
                LinearSizeControl(α, β1); V0=β1, n=burn, cv=CV, daughter_fraction=fr, seed=s
            ).Vb,
        )
        vb = simulate_lineage(
            LinearSizeControl(α, β2); V0=v0, n=M + 1, cv=CV, daughter_fraction=fr, seed=s + 10^6
        ).Vb
        traj .+= vb
    end
    return traj ./ R
end

open(joinpath(here, "size_memory_step.csv"), "w") do io
    println(io, "mode,alpha,r,memory_gen,gen,vb_measured,vb_geom,vstar1,vstar2")
    fr = 0.5
    for (name, α) in modes
        r = α * fr
        mem = size_memory(LinearSizeControl(α, BETA); daughter_fraction=fr).memory_gen
        v1 = fixedpoint(α, 20.0, fr)
        v2 = fixedpoint(α, 40.0, fr)
        tr = step_response(α, fr)
        for n in 0:20
            geom = v2 + (v1 - v2) * r^n
            @printf(
                io, "%s,%.2f,%.4f,%.4f,%d,%.6f,%.6f,%.4f,%.4f\n",
                name, α, r, mem, n, tr[n + 1], geom, v1, v2
            )
        end
    end
end

println("wrote size_memory_autocorr.csv, size_memory_invariant.csv, size_memory_step.csv")
