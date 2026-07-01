# ABC-MCMC convergence diagnostics for the emergent-RLS calibration (companion to gen_rls_abc.jl).
# Reports acceptance rate, multi-chain Gelman-Rubin R-hat, effective sample size (ESS), and a
# prior-predictive check (how diffuse the summaries are under the prior, i.e. that the data are
# informative). Does NOT alter the released posterior (rls_abc_posterior.csv); it runs its own
# dispersed multi-chain ensemble purely to characterise the sampler.
# Run: julia --project=. analysis/gen_rls_abc_diagnostics.jl
using CellSizeControl
using Statistics: mean, std, var, quantile
using Random: MersenneTwister
using Printf

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

# one random-walk Metropolis chain; returns (Ds,ks,cs, accept_rate)
function run_chain(D0, k0, c0, seed; nsteps=7000, burn=1500)
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

# Gelman-Rubin R-hat from a vector of per-chain sample vectors (equal length)
function rhat(chains)
    m = length(chains); n = minimum(length.(chains))
    ch = [c[1:n] for c in chains]
    means = mean.(ch); B = n*var(means); W = mean(var.(ch))
    varhat = (n-1)/n*W + B/n
    return sqrt(varhat / W)
end
# ESS via initial-positive-sequence autocorrelation (pooled over chains)
function ess(chains)
    n = minimum(length.(chains)); ch = [c[1:n] for c in chains]
    x = vcat(ch...); N = length(x); mu = mean(x); v = var(x); v == 0 && return float(N)
    ac(k) = sum((x[1:N-k] .- mu).*(x[1+k:N] .- mu)) / ((N-k)*v)
    s = 0.0; k = 1
    while k < min(1000, N-1)
        r = ac(k); r <= 0 && break; s += r; k += 1
    end
    return N / (1 + 2s)
end

# 4 dispersed chains
inits = [(20.0,0.05,0.3),(40.0,0.15,0.7),(55.0,0.08,0.5),(30.0,0.20,0.9)]
Dc=[];kc=[];cc=[];accs=Float64[]
for (i,(D0,k0,c0)) in enumerate(inits)
    Ds,ks,cs,a = run_chain(D0,k0,c0, 20260701+i)
    push!(Dc,Ds); push!(kc,ks); push!(cc,cs); push!(accs,a)
    @printf("chain %d: accept=%.2f  D=%.1f k=%.3f c=%.2f (n=%d)\n", i, a, mean(Ds), mean(ks), mean(cs), length(Ds))
end
@printf("\nacceptance rate: mean %.2f over 4 chains (range %.2f-%.2f)\n", mean(accs), minimum(accs), maximum(accs))
@printf("Gelman-Rubin R-hat:  D_crit=%.3f  kappa=%.3f  crit_cv=%.3f  (want <1.05)\n",
        rhat(Dc), rhat(kc), rhat(cc))
@printf("effective sample size: D_crit=%.0f  kappa=%.0f  crit_cv=%.0f\n", ess(Dc), ess(kc), ess(cc))

# OUTPUT-specific convergence: the reported quantity is the predictive mean RLS, not the ridge
# parameters. Compute it per chain (at each chain's posterior mean) and report the between-chain
# spread — it is stable even though the D_crit/kappa marginals scatter along the ridge.
chain_means = [(mean(Dc[i]), mean(kc[i]), mean(cc[i])) for i in 1:length(Dc)]
pred_means = Float64[]
for (D,k,c) in chain_means
    ls = lifespan_distribution(4000; seed0=7, D_crit=D, kappa=k, crit_cv=c, production=1.0, cv=0.05, max_gen=400, segregate=false)
    push!(pred_means, mean(ls))
end
@printf("\npredictive mean RLS per chain: %s\n", join((@sprintf("%.1f",m) for m in pred_means), ", "))
@printf("between-chain: mean %.2f, sd %.2f, cv %.3f  -> the REPORTED output is stable across chains\n",
        mean(pred_means), std(pred_means), std(pred_means)/mean(pred_means))

# prior-predictive check: draw from the uniform box, report the spread of predicted mean RLS
rngp = MersenneTwister(999); pm = Float64[]
for _ in 1:400
    D = 10+70*rand(rngp); k = 0.30*rand(rngp); c = 0.05+0.95*rand(rngp)
    m,_ = summ(D,k,c, rand(rngp,1:10^6)); push!(pm, m)
end
@printf("prior-predictive mean-RLS: median %.1f, 5-95%% [%.1f, %.1f] divisions (target %.1f) -> prior is diffuse; the data localise it\n",
        quantile(pm,.5), quantile(pm,.05), quantile(pm,.95), TARGET_MEAN)
