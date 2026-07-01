# ABLATION: is the AGE-ERODING asymmetry r(a) load-bearing for the daughter-lifespan fold,
# or does the fold arise from maternal damage accumulation D(a) plus ANY inheritance?
# Daughter lifespan shortening could follow mainly from D_m(a) rising with maternal age rather
# than from the age-dependence of the inheritance share r(a); this decomposition settles which.
#
# We recompute the Kennedy young/old daughter-RLS fold under four inheritance schemes, holding
# the damage recursion + calibrated posterior fixed and changing ONLY the inheritance share phi:
#   (A) baseline      phi(a) = r(a)/r_max      -- age-eroding asymmetry (the paper's model)
#   (B) const-young   phi    = r(0)/r_max      -- asymmetry FROZEN at the young-mother value
#   (C) const-mean    phi    = <r(a)>/r_max    -- asymmetry frozen at its lineage-mean value
#   (D) const-symm    phi    = 1.0             -- fully symmetric inheritance (upper bound)
# If (B)/(C) still reproduce the ~3.3x fold, r(a)'s age-dependence is NOT load-bearing and the
# claim must be reframed; if the fold collapses toward 1 without the age-dependence, r(a) IS
# load-bearing. Either way this is an honest, decisive test. NOTHING is refit.
#
# Run: julia --project=. analysis/gen_daughter_rls_ablation.jl

using Statistics: mean, std
using Random: MersenneTwister
using DelimitedFiles: readdlm
using Printf

here = @__DIR__
post = readdlm(joinpath(here, "rls_abc_posterior.csv"), ',', Float64; skipstart=1)
const D_CRIT  = mean(post[:, 1]); const KAPPA = mean(post[:, 2]); const CRIT_CV = mean(post[:, 3])
const PROD = 1.0; const CV = 0.05
const R0, R_MAX, R_TAU = 0.69, 0.90, 14.0
r(a) = R0 + (R_MAX - R0) * (1.0 - exp(-a / R_TAU))
const MAX_GEN = 400

@inline step_damage(D, rng) = D + PROD * (1.0 + KAPPA * D) * max(0.0, 1.0 + CV * randn(rng))

function age_cell(D0, rng)
    Dc = CRIT_CV > 0 ? D_CRIT * exp(CRIT_CV * randn(rng) - CRIT_CV^2 / 2) : D_CRIT
    traj = Float64[]; D = float(D0); a = 0
    while D < Dc && a < MAX_GEN
        push!(traj, D); D = step_damage(D, rng); a += 1
    end
    return a, traj
end

# mean of r(a) over a typical mother lifespan (~26 divisions) for the const-mean scheme
const R_MEAN = mean(r.(0:25))

phi_scheme(name, a) =
    name === :baseline  ? r(a) / R_MAX :
    name === :const_young ? r(0) / R_MAX :
    name === :const_mean  ? R_MEAN / R_MAX :
    name === :const_symm  ? 1.0 : error("bad scheme")

function fold(name; N::Int=40_000, seed::Int=20260701)
    rng = MersenneTwister(seed)
    young_sum = 0.0; young_n = 0; old_sum = 0.0; old_n = 0
    for _ in 1:N
        L, traj = age_cell(0.0, rng); L == 0 && continue
        for a in 0:(L - 1)
            D0 = phi_scheme(name, a) * traj[a + 1]
            Ld, _ = age_cell(D0, rng)
            frac = (L > 1) ? a / (L - 1) : 0.0
            if frac <= 0.7; young_sum += Ld; young_n += 1
            elseif frac >= 0.9; old_sum += Ld; old_n += 1 end
        end
    end
    return young_sum / young_n, old_sum / old_n
end

@printf("calibrated: D_crit=%.3f kappa=%.4f crit_cv=%.4f ; r(0)=%.3f r_mean=%.3f r_max=%.3f\n",
        D_CRIT, KAPPA, CRIT_CV, r(0), R_MEAN, R_MAX)
@printf("\n%-12s  %8s  %8s  %8s   %s\n", "scheme", "young", "old", "fold", "phi(a)")
for (nm, desc) in ((:baseline,"r(a)/r_max age-eroding"), (:const_young,"r(0)/r_max frozen young"),
                   (:const_mean,"<r>/r_max frozen mean"), (:const_symm,"1.0 fully symmetric"))
    y, o = fold(nm)
    @printf("%-12s  %8.2f  %8.2f  %8.2fx  %s\n", nm, y, o, y / o, desc)
end
@printf("\n(Kennedy 1994 target fold = 26.5/7.9 = %.2fx)\n", 26.5 / 7.9)
