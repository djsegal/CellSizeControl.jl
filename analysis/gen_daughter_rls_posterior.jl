# Posterior-predictive uncertainty for the out-of-sample daughter-RLS deficit.
#
# gen_daughter_rls.jl reports the Kennedy buckets at the ABC posterior MEAN only. This script
# propagates the full McCormick-2015 ABC posterior: it draws (D_crit, kappa, crit_cv) triples
# from rls_abc_posterior.csv and, for each draw, re-runs the daughter-RLS simulation, recording
# the young-mother bucket (first 70% of mother life), the old-mother bucket (last 10%), and the
# fold deficit. The spread across draws is the posterior-predictive credible interval on the
# prediction -- the uncertainty the point estimate hides.
#
# It also does a one-at-a-time sensitivity of the fold and absolute buckets to the size-face
# asymmetry parameters (r0, r_max, tau_r), which are illustrative / literature-anchored rather
# than regression-fit, to show how much the Kennedy comparison leans on them.
#
# NOTHING is fitted to Kennedy: (D_crit, kappa, crit_cv) are the WT lifespan posterior and
# (r0, r_max, tau_r) are the size-face asymmetry fixed by the daughter-size increase.
#
# Run: julia --project=. gen_daughter_rls_posterior.jl

using Statistics: mean, std, median, quantile
using Random: MersenneTwister
using DelimitedFiles: readdlm
using Printf

here = @__DIR__

# --- size-face asymmetry r(a): the bud volume fraction, fixed by the daughter-size data -------
const R0, R_MAX, R_TAU = 0.69, 0.90, 14.0
const PROD = 1.0       # damage formed per cycle (a.u.)
const CV   = 0.05      # per-division multiplicative noise (as in the ABC calibration)
const MAX_GEN = 400

r_of(a, r0, rmax, rtau) = r0 + (rmax - r0) * (1.0 - exp(-a / rtau))

@inline step_damage(D, kappa, rng) =
    D + PROD * (1.0 + kappa * D) * max(0.0, 1.0 + CV * randn(rng))

# Age a cell from seed damage D0 with its own lognormal viability threshold; return
# (lifespan, damage trajectory carried into each division).
function age_cell(D0, D_crit, kappa, crit_cv, rng)
    Dc = crit_cv > 0 ? D_crit * exp(crit_cv * randn(rng) - crit_cv^2 / 2) : D_crit
    traj = Float64[]
    D = float(D0)
    a = 0
    while D < Dc && a < MAX_GEN
        push!(traj, D)
        D = step_damage(D, kappa, rng)
        a += 1
    end
    return a, traj
end

# One full daughter-RLS run at a parameter set; returns (young_mean, old_mean, fold, mother_mean).
function run_buckets(D_crit, kappa, crit_cv, r0, rmax, rtau; N::Int, rng)
    young_sum = 0.0; young_n = 0
    old_sum = 0.0; old_n = 0
    mother_sum = 0.0; mother_n = 0
    phi(a) = r_of(a, r0, rmax, rtau) / rmax
    for _ in 1:N
        L, traj = age_cell(0.0, D_crit, kappa, crit_cv, rng)
        L == 0 && continue
        mother_sum += L; mother_n += 1
        for a in 0:(L - 1)
            Dm = traj[a + 1]
            D0 = phi(a) * Dm
            Ld, _ = age_cell(D0, D_crit, kappa, crit_cv, rng)
            frac = (L > 1) ? a / (L - 1) : 0.0
            if frac <= 0.7
                young_sum += Ld; young_n += 1
            elseif frac >= 0.9
                old_sum += Ld; old_n += 1
            end
        end
    end
    ym = young_sum / young_n
    om = old_sum / old_n
    return ym, om, ym / om, mother_sum / mother_n
end

function main(; ndraw::Int=400, N::Int=6000, seed::Int=20260627)
    post = readdlm(joinpath(here, "rls_abc_posterior.csv"), ',', Float64; skipstart=1)
    ndraw = min(ndraw, size(post, 1))
    # evenly thin the chain to ndraw draws (decorrelates better than the first ndraw rows)
    idx = round.(Int, range(1, size(post, 1); length=ndraw))
    rng = MersenneTwister(seed)

    young = Float64[]; old = Float64[]; fold = Float64[]; mom = Float64[]
    t0 = time()
    for (k, i) in enumerate(idx)
        Dc, ka, cv = post[i, 1], post[i, 2], post[i, 3]
        ym, om, fd, mm = run_buckets(Dc, ka, cv, R0, R_MAX, R_TAU; N=N, rng=rng)
        push!(young, ym); push!(old, om); push!(fold, fd); push!(mom, mm)
        if k % 50 == 0
            @printf("  draw %d/%d  (%.0fs)  young=%.1f old=%.1f fold=%.2f\n",
                    k, ndraw, time() - t0, ym, om, fd)
        end
    end

    qs = (0.025, 0.5, 0.975)
    yq = quantile(young, qs); oq = quantile(old, qs); fq = quantile(fold, qs)
    mq = quantile(mom, qs)

    println("\n=== Posterior-predictive Kennedy buckets ($ndraw draws x $N mothers) ===")
    @printf("  young-mother daughters (first 70%%): median %.1f  95%% CI [%.1f, %.1f]   (Kennedy 26.5)\n",
            yq[2], yq[1], yq[3])
    @printf("  old-mother   daughters (last 10%%) : median %.1f  95%% CI [%.1f, %.1f]   (Kennedy  7.9)\n",
            oq[2], oq[1], oq[3])
    @printf("  fold deficit                      : median %.2f 95%% CI [%.2f, %.2f]   (Kennedy 3.35)\n",
            fq[2], fq[1], fq[3])
    @printf("  pooled mother RLS                 : median %.1f  95%% CI [%.1f, %.1f]   (McCormick 26.6)\n",
            mq[2], mq[1], mq[3])

    open(joinpath(here, "daughter_rls_posterior.csv"), "w") do io
        println(io, "draw,young_rls,old_rls,fold,mother_rls")
        for k in 1:length(young)
            @printf(io, "%d,%.4f,%.4f,%.4f,%.4f\n", k, young[k], old[k], fold[k], mom[k])
        end
    end
    open(joinpath(here, "daughter_rls_posterior_summary.csv"), "w") do io
        println(io, "quantity,median,lo95,hi95,kennedy")
        @printf(io, "young_rls,%.4f,%.4f,%.4f,26.5\n", yq[2], yq[1], yq[3])
        @printf(io, "old_rls,%.4f,%.4f,%.4f,7.9\n", oq[2], oq[1], oq[3])
        @printf(io, "fold,%.4f,%.4f,%.4f,3.354\n", fq[2], fq[1], fq[3])
        @printf(io, "mother_rls,%.4f,%.4f,%.4f,26.6\n", mq[2], mq[1], mq[3])
    end

    # --- one-at-a-time sensitivity to the size-face asymmetry, at the posterior mean ----------
    Dc0, ka0, cv0 = mean(post[:, 1]), mean(post[:, 2]), mean(post[:, 3])
    rngs = MersenneTwister(seed + 1)
    base = run_buckets(Dc0, ka0, cv0, R0, R_MAX, R_TAU; N=20_000, rng=rngs)
    variants = [
        ("base", R0, R_MAX, R_TAU),
        ("r0=0.66", 0.66, R_MAX, R_TAU),
        ("r0=0.72", 0.72, R_MAX, R_TAU),
        ("rmax=0.86", R0, 0.86, R_TAU),
        ("rmax=0.94", R0, 0.94, R_TAU),
        ("tau_r=10", R0, R_MAX, 10.0),
        ("tau_r=18", R0, R_MAX, 18.0),
    ]
    println("\n=== Sensitivity of the Kennedy buckets to the (illustrative) asymmetry params ===")
    open(joinpath(here, "daughter_rls_sensitivity.csv"), "w") do io
        println(io, "variant,r0,r_max,tau_r,young_rls,old_rls,fold")
        for (name, r0, rmax, rtau) in variants
            ym, om, fd, _ = run_buckets(Dc0, ka0, cv0, r0, rmax, rtau; N=20_000, rng=rngs)
            @printf("  %-10s young=%.1f old=%.1f fold=%.2f\n", name, ym, om, fd)
            @printf(io, "%s,%.3f,%.3f,%.1f,%.4f,%.4f,%.4f\n", name, r0, rmax, rtau, ym, om, fd)
        end
    end
    @printf("\n(base posterior-mean run: young=%.1f old=%.1f fold=%.2f)\n", base[1], base[2], base[3])
    println("\nwrote daughter_rls_posterior.csv, daughter_rls_posterior_summary.csv, daughter_rls_sensitivity.csv")
    return nothing
end

main()
