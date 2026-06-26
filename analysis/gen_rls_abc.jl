# CC-4: ABC-MCMC calibration of the emergent-RLS model to a REAL wild-type budding-yeast
# replicative-lifespan distribution. Target = McCormick et al. 2015 (Cell Metab 22:895-906),
# the largest WT microdissection dataset (pooled BY4741/BY4742, n=29,383): mean RLS = 26.6
# divisions, SD = 9.7 (CV = 0.365). We infer the posterior over the three damage parameters
# (D_crit, kappa, crit_cv) of `replicative_lifespan` via ABC on the summary stats (mean, SD).
# Run: julia --project=. gen_rls_abc.jl
using CellSizeControl
using Statistics: mean, std, quantile
using Random: MersenneTwister
using Printf

here = @__DIR__
const TARGET_MEAN, TARGET_SD = 26.6, 9.7      # McCormick 2015 WT pooled
const M = 1200                                  # lifespan samples per summary evaluation
const ABC_MEAN_BW, ABC_SD_BW = 0.8, 0.8         # ABC kernel bandwidths (divisions)

# summary: simulate M lifespans at θ=(D_crit,kappa,crit_cv) → (mean, sd). seed varies per call.
function summary(D_crit, kappa, crit_cv, seed)
    ls = lifespan_distribution(M; seed0=seed, D_crit=D_crit, kappa=kappa, crit_cv=crit_cv,
                               production=1.0, cv=0.05, max_gen=400, segregate=false)
    return mean(ls), std(ls)
end

# log ABC pseudo-likelihood (Gaussian kernel on the two summaries) + uniform-box prior
inbox(D, k, c) = (10.0 ≤ D ≤ 80.0) && (0.0 ≤ k ≤ 0.30) && (0.05 ≤ c ≤ 1.00)
function logpost(D, k, c, seed)
    inbox(D, k, c) || return -Inf
    m, s = summary(D, k, c, seed)
    return -0.5 * ((m - TARGET_MEAN) / ABC_MEAN_BW)^2 - 0.5 * ((s - TARGET_SD) / ABC_SD_BW)^2
end

# anchor: where do the package defaults land?
let (m, s) = summary(32.0, 0.10, 0.55, 1)
    @printf("defaults (D=32,k=0.10,c=0.55): mean=%.2f sd=%.2f  (target %.1f / %.1f)\n",
            m, s, TARGET_MEAN, TARGET_SD)
end

# ABC-MCMC (random-walk Metropolis); re-draw the simulation seed each step so the kernel
# integrates over simulation noise (a pseudo-marginal-style move).
rng = MersenneTwister(20260626)
nsteps, burn = 16000, 3000
D, k, c = 32.0, 0.10, 0.55
lp = logpost(D, k, c, 1)
stepD, stepK, stepC = 3.0, 0.02, 0.06
chain = Vector{NTuple{5,Float64}}()
nacc = 0
for it in 1:nsteps
    global D, k, c, lp, nacc
    Dp = D + stepD * randn(rng)
    kp = k + stepK * randn(rng)
    cp = c + stepC * randn(rng)
    lpp = logpost(Dp, kp, cp, 1000 + it)        # fresh seed each proposal
    if log(rand(rng)) < lpp - lp
        D, k, c, lp = Dp, kp, cp, lpp
        nacc += 1
    end
    it > burn && push!(chain, (D, k, c, lp, 0.0))
    it % 2000 == 0 && @printf("  step %d/%d  acc=%.2f  D=%.1f k=%.3f c=%.2f\n",
                              it, nsteps, nacc / it, D, k, c)
end

open(joinpath(here, "rls_abc_posterior.csv"), "w") do io
    println(io, "D_crit,kappa,crit_cv")
    for (D, k, c, _, _) in chain
        @printf(io, "%.4f,%.5f,%.5f\n", D, k, c)
    end
end

Ds = [x[1] for x in chain]; ks = [x[2] for x in chain]; cs = [x[3] for x in chain]
pe = (mean(Ds), mean(ks), mean(cs))           # posterior means
@printf("\nposterior mean: D_crit=%.2f [%.2f,%.2f]  kappa=%.4f [%.4f,%.4f]  crit_cv=%.3f [%.3f,%.3f]\n",
        pe[1], quantile(Ds, .025), quantile(Ds, .975),
        pe[2], quantile(ks, .025), quantile(ks, .975),
        pe[3], quantile(cs, .025), quantile(cs, .975))

# posterior-predictive RLS distribution at the posterior mean (a big sample for the figure)
pp = lifespan_distribution(20000; seed0=77, D_crit=pe[1], kappa=pe[2], crit_cv=pe[3],
                           production=1.0, cv=0.05, max_gen=400, segregate=false)
open(joinpath(here, "rls_abc_predictive.csv"), "w") do io
    println(io, "rls")
    for v in pp
        println(io, v)
    end
end
@printf("posterior-predictive RLS: mean=%.2f sd=%.2f cv=%.3f  (target %.1f / %.1f / %.3f)\n",
        mean(pp), std(pp), std(pp) / mean(pp), TARGET_MEAN, TARGET_SD, TARGET_SD / TARGET_MEAN)
println("wrote rls_abc_posterior.csv (", length(chain), " draws) + rls_abc_predictive.csv")
