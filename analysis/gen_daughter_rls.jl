# Out-of-sample prediction: daughter replicative lifespan vs maternal age.
#
# The single age-eroding division asymmetry r(a) ties the SIZE face (old mothers make larger
# daughters) to the FITNESS face. Under passive volume-proportional damage segregation the bud
# of an old, nearly-symmetric division carries the fraction phi(a) = r(a)/r_max of the mother's
# accumulated damage, so a daughter is BORN partway up the autocatalytic damage trajectory and
# her own emergent replicative lifespan is shortened. NOTHING here is fitted: the damage
# parameters (D_crit, kappa, crit_cv) are the McCormick-2015 WT ABC posterior used for the
# mother lifespan, and (r0, r_max, r_tau) are the size-face asymmetry already fixed by the
# Johnston-1966 / Yang-2011 daughter-size increase. The Kennedy-1994 daughter-vs-mother-age
# numbers (26.5 div from young mothers, 7.9 div from old mothers) are NOT used to build this --
# they are the independent test target.
#
# Self-consistency note: the lineage figure uses a LINEAR damage accrual for the illustrative
# damage axis, but a replicative LIFESPAN is by construction the crossing time of the
# AUTOCATALYTIC damage state, so this prediction uses the one autocatalytic damage variable
# throughout (mother and daughter follow the same recursion; the daughter is seeded by passive
# inheritance). This is the only coherent way to turn inherited damage into an inherited lifespan.
#
# Run: julia --project=. gen_daughter_rls.jl

using Statistics: mean, std
using Random: MersenneTwister
using DelimitedFiles: readdlm
using Printf

here = @__DIR__

# --- calibrated damage parameters: read the McCormick-2015 ABC posterior MEAN (no re-fit) ----
post = readdlm(joinpath(here, "rls_abc_posterior.csv"), ',', Float64; skipstart=1)
const D_CRIT  = mean(post[:, 1])
const KAPPA   = mean(post[:, 2])
const CRIT_CV = mean(post[:, 3])
const PROD    = 1.0     # damage formed per cycle (a.u.), as in replicative_lifespan
const CV      = 0.05    # per-division multiplicative noise, as used in the ABC calibration
@printf("calibrated damage params (McCormick-2015 ABC posterior mean): D_crit=%.3f kappa=%.4f crit_cv=%.4f\n",
        D_CRIT, KAPPA, CRIT_CV)

# --- size-face asymmetry r(a): the bud volume fraction, fixed by the daughter-size data -------
const R0, R_MAX, R_TAU = 0.69, 0.90, 14.0
r(a)   = R0 + (R_MAX - R0) * (1.0 - exp(-a / R_TAU))   # rising bud fraction (asymmetry erosion)
phi(a) = r(a) / R_MAX                                   # passive volume-proportional inheritance share

const MAX_GEN = 400

# Autocatalytic damage step (identical recursion to CellSizeControl.replicative_lifespan with
# segregate=false): D <- D + production*(1+kappa*D)*max(0, 1+cv*randn).
@inline step_damage(D, rng) = D + PROD * (1.0 + KAPPA * D) * max(0.0, 1.0 + CV * randn(rng))

# Age a cell from a seed damage D0 with its OWN fresh lognormal viability threshold; return
# (replicative lifespan, full damage trajectory D[0..L-1] = damage present when it buds age a).
function age_cell(D0, rng)
    Dc = CRIT_CV > 0 ? D_CRIT * exp(CRIT_CV * randn(rng) - CRIT_CV^2 / 2) : D_CRIT
    traj = Float64[]
    D = float(D0)
    a = 0
    while D < Dc && a < MAX_GEN
        push!(traj, D)          # damage carried into the division that buds the age-a daughter
        D = step_damage(D, rng)
        a += 1
    end
    return a, traj
end

function main(; N::Int=40_000, seed::Int=20260627)
    rng = MersenneTwister(seed)

    # accumulate daughter RLS by maternal-age fraction bin and by absolute maternal age
    nbin = 20
    binsum = zeros(nbin); binsq = zeros(nbin); binN = zeros(Int, nbin)
    # Kennedy 1994 buckets: first 70% of mother life, and last 10% of mother life
    young_sum = 0.0; young_n = 0          # a/L <= 0.7
    old_sum   = 0.0; old_n   = 0          # a/L >= 0.9
    # absolute maternal age (generations) -> mean daughter RLS
    amax = 60
    absum = zeros(amax); abN = zeros(Int, amax)
    mother_rls = Int[]

    for _ in 1:N
        L, traj = age_cell(0.0, rng)      # a fresh mother: seed damage 0
        L == 0 && continue
        push!(mother_rls, L)
        for a in 0:(L - 1)
            Dm = traj[a + 1]              # mother damage when budding the age-a daughter
            D0 = phi(a) * Dm             # passive volume-proportional inheritance
            Ld, _ = age_cell(D0, rng)    # the daughter's own emergent lifespan
            frac = (L > 1) ? a / (L - 1) : 0.0   # fraction of the mother life completed
            b = clamp(floor(Int, frac * nbin) + 1, 1, nbin)
            binsum[b] += Ld; binsq[b] += Ld^2; binN[b] += 1
            if frac <= 0.7
                young_sum += Ld; young_n += 1
            elseif frac >= 0.9
                old_sum += Ld; old_n += 1
            end
            if a + 1 <= amax
                absum[a + 1] += Ld; abN[a + 1] += 1
            end
        end
    end

    @printf("\nmother RLS ensemble: mean=%.2f sd=%.2f cv=%.3f  (McCormick target 26.6/9.7/0.365)\n",
            mean(mother_rls), std(mother_rls), std(mother_rls) / mean(mother_rls))

    young_mean = young_sum / young_n
    old_mean   = old_sum / old_n
    @printf("\n=== Kennedy-1994 buckets (INDEPENDENT test; not used in any fit) ===\n")
    @printf("  daughters from mothers in FIRST 70%% of life : model %.1f div   (Kennedy 26.5)\n", young_mean)
    @printf("  daughters from mothers in LAST  10%% of life : model %.1f div   (Kennedy  7.9)\n", old_mean)
    @printf("  model fold-drop %.2fx  (Kennedy 26.5/7.9 = %.2fx)\n", young_mean / old_mean, 26.5 / 7.9)

    # write the maternal-age-fraction curve (the figure backbone)
    open(joinpath(here, "daughter_rls_fraction.csv"), "w") do io
        println(io, "frac_lo,frac_mid,frac_hi,daughter_rls_mean,daughter_rls_sd,n")
        for b in 1:nbin
            binN[b] == 0 && continue
            m = binsum[b] / binN[b]
            sd = sqrt(max(0.0, binsq[b] / binN[b] - m^2))
            @printf(io, "%.4f,%.4f,%.4f,%.4f,%.4f,%d\n",
                    (b - 1) / nbin, (b - 0.5) / nbin, b / nbin, m, sd, binN[b])
        end
    end

    # write the absolute-maternal-age curve (auxiliary)
    open(joinpath(here, "daughter_rls_abs_age.csv"), "w") do io
        println(io, "maternal_age,daughter_rls_mean,n")
        for a in 1:amax
            abN[a] == 0 && continue
            @printf(io, "%d,%.4f,%d\n", a - 1, absum[a] / abN[a], abN[a])
        end
    end

    # write the two-bucket summary + Kennedy targets for the plotter
    open(joinpath(here, "daughter_rls_kennedy.csv"), "w") do io
        println(io, "bucket,frac_lo,frac_hi,model_rls,kennedy_rls")
        @printf(io, "first70,0.0,0.7,%.4f,26.5\n", young_mean)
        @printf(io, "last10,0.9,1.0,%.4f,7.9\n", old_mean)
    end

    println("\nwrote daughter_rls_fraction.csv, daughter_rls_abs_age.csv, daughter_rls_kennedy.csv")
    return nothing
end

main()
