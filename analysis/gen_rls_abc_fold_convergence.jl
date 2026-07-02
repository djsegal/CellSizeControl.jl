# Predictive convergence of the HEADLINE quantity — the old/young daughter-RLS fold deficit —
# across independent ABC-MCMC chains. Companion to gen_rls_abc_diagnostics.jl (which shows the
# ridge parameters D_crit/kappa mix slowly, as expected for a non-identifiable ridge) and to
# gen_daughter_rls_posterior.jl (which reports the fold CI from the released single chain).
#
# The point (reviewer ask): the reported claim is the fold, not the ridge parameters. Predictions
# along the ridge are identifiable even though the marginals are not, so the fold should be stable
# across dispersed chains. This script quantifies that: 4 dispersed chains, the fold posterior
# per chain, Gelman-Rubin R-hat and ESS on the fold, and the pooled 95% CI.
# Run: julia --project=. analysis/gen_rls_abc_fold_convergence.jl
using CellSizeControl
using Statistics: mean, std, var, quantile
using Random: MersenneTwister
using Printf

# ---- ABC sampler (identical to gen_rls_abc_diagnostics.jl) ----------------------------------
const TARGET_MEAN, TARGET_SD = 26.6, 9.7
const M = 1200
const ABC_MEAN_BW, ABC_SD_BW = 0.8, 0.8
summ(D, k, c, seed) = (ls = lifespan_distribution(M; seed0=seed, D_crit=D, kappa=k, crit_cv=c,
                        production=1.0, cv=0.05, max_gen=400, segregate=false); (mean(ls), std(ls)))
inbox(D, k, c) = (10.0 ≤ D ≤ 80.0) && (0.0 ≤ k ≤ 0.30) && (0.05 ≤ c ≤ 1.00)
function logpost(D, k, c, seed)
    inbox(D, k, c) || return -Inf
    m, s = summ(D, k, c, seed)
    return -0.5 * ((m - TARGET_MEAN) / ABC_MEAN_BW)^2 - 0.5 * ((s - TARGET_SD) / ABC_SD_BW)^2
end
function run_chain(D0, k0, c0, seed; nsteps=6000, burn=1500)
    rng = MersenneTwister(seed)
    D, k, c = D0, k0, c0; lp = logpost(D, k, c, seed)
    stepD, stepK, stepC = 3.0, 0.02, 0.06
    Ds = Float64[]; ks = Float64[]; cs = Float64[]; nacc = 0
    for it in 1:nsteps
        Dp = D + stepD*randn(rng); kp = k + stepK*randn(rng); cp = c + stepC*randn(rng)
        lpp = logpost(Dp, kp, cp, seed*10^6 + it)
        if log(rand(rng)) < lpp - lp; D,k,c,lp = Dp,kp,cp,lpp; nacc += 1 end
        if it > burn; push!(Ds,D); push!(ks,k); push!(cs,c) end
    end
    return Ds, ks, cs, nacc/nsteps
end
function rhat(chains)
    m = length(chains); n = minimum(length.(chains))
    ch = [c[1:n] for c in chains]
    means = mean.(ch); B = n*var(means); W = mean(var.(ch))
    return sqrt(((n-1)/n*W + B/n) / W)
end
function ess(chains)
    n = minimum(length.(chains)); ch = [c[1:n] for c in chains]
    x = vcat(ch...); N = length(x); mu = mean(x); v = var(x); v == 0 && return float(N)
    ac(k) = sum((x[1:N-k] .- mu).*(x[1+k:N] .- mu)) / ((N-k)*v)
    s = 0.0; k = 1
    while k < min(1000, N-1); r = ac(k); r <= 0 && break; s += r; k += 1 end
    return N / (1 + 2s)
end

# ---- fold machinery (identical to gen_daughter_rls_posterior.jl) -----------------------------
const R0, R_MAX, R_TAU = 0.69, 0.90, 14.0
const PROD, CV, MAX_GEN = 1.0, 0.05, 400
r_of(a, r0, rmax, rtau) = r0 + (rmax - r0) * (1.0 - exp(-a / rtau))
@inline step_damage(D, kappa, rng) = D + PROD * (1.0 + kappa * D) * max(0.0, 1.0 + CV * randn(rng))
function age_cell(D0, D_crit, kappa, crit_cv, rng)
    Dc = crit_cv > 0 ? D_crit * exp(crit_cv * randn(rng) - crit_cv^2 / 2) : D_crit
    traj = Float64[]; D = float(D0); a = 0
    while D < Dc && a < MAX_GEN; push!(traj, D); D = step_damage(D, kappa, rng); a += 1 end
    return a, traj
end
function fold_at(D_crit, kappa, crit_cv; N::Int, rng)
    young_sum = 0.0; young_n = 0; old_sum = 0.0; old_n = 0
    phi(a) = r_of(a, R0, R_MAX, R_TAU) / R_MAX
    for _ in 1:N
        L, traj = age_cell(0.0, D_crit, kappa, crit_cv, rng)
        L == 0 && continue
        for a in 0:(L - 1)
            Ld, _ = age_cell(phi(a) * traj[a + 1], D_crit, kappa, crit_cv, rng)
            frac = (L > 1) ? a / (L - 1) : 0.0
            frac <= 0.7 ? (young_sum += Ld; young_n += 1) : (frac >= 0.9 ? (old_sum += Ld; old_n += 1) : nothing)
        end
    end
    ym = young_sum / max(young_n,1); om = old_sum / max(old_n,1)
    return ym / max(om, 1e-9)
end

# ---- 4 dispersed chains -> per-chain fold posterior -----------------------------------------
inits = [(20.0,0.05,0.3),(40.0,0.15,0.7),(55.0,0.08,0.5),(30.0,0.20,0.9)]
const DRAWS_PER_CHAIN = 50   # thinned draws per chain fed through the fold simulator
const N_CELLS = 3000         # lineage sample per fold evaluation
here = @__DIR__
foldchains = Vector{Vector{Float64}}()
t0 = time()
for (i,(D0,k0,c0)) in enumerate(inits)
    Ds,ks,cs,acc = run_chain(D0,k0,c0, 20260701+i)
    idx = round.(Int, range(1, length(Ds); length=DRAWS_PER_CHAIN))
    rng = MersenneTwister(90000+i)
    fs = [fold_at(Ds[j], ks[j], cs[j]; N=N_CELLS, rng=rng) for j in idx]
    push!(foldchains, fs)
    @printf("chain %d: accept=%.2f  fold median=%.2f  [%.2f, %.2f]  (%.0fs)\n",
            i, acc, quantile(fs,.5), quantile(fs,.025), quantile(fs,.975), time()-t0)
end
allf = vcat(foldchains...)
@printf("\n=== daughter-fold predictive convergence across 4 dispersed chains ===\n")
@printf("per-chain fold medians: %s\n", join((@sprintf("%.2f",quantile(c,.5)) for c in foldchains), ", "))
@printf("Gelman-Rubin R-hat (fold): %.3f   (want <1.05)\n", rhat(foldchains))
@printf("ESS (fold, pooled): %.0f of %d draws\n", ess(foldchains), length(allf))
@printf("pooled fold: median %.2f, 95%% CI [%.2f, %.2f]   (Kennedy 3.35)\n",
        quantile(allf,.5), quantile(allf,.025), quantile(allf,.975))
open(joinpath(here, "rls_abc_fold_convergence.csv"), "w") do io
    println(io, "chain,fold_median,fold_lo,fold_hi")
    for (i,c) in enumerate(foldchains)
        @printf(io, "%d,%.4f,%.4f,%.4f\n", i, quantile(c,.5), quantile(c,.025), quantile(c,.975))
    end
end
println("wrote rls_abc_fold_convergence.csv")
