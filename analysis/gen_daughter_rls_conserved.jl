# Robustness of the daughter-lifespan fold to a CONSERVED damage partition.
# The main model tracks inherited damage as a non-conserved state (the daughter's share is NOT
# subtracted from the mother). Here we test the alternative: at each division the daughter inherits
# phi(a)*D_m AND that share is removed from the mother (strict conservation of a damage pool).
# We report the mother RLS and the young/old daughter fold under both, so the manuscript can state
# whether the qualitative result is an artifact of non-conservation.
# Run: julia --project=. analysis/gen_daughter_rls_conserved.jl
using Statistics: mean, std
using Random: MersenneTwister
using DelimitedFiles: readdlm
using Printf

here = @__DIR__
post = readdlm(joinpath(here, "rls_abc_posterior.csv"), ',', Float64; skipstart=1)
const D_CRIT = mean(post[:,1]); const KAPPA = mean(post[:,2]); const CRIT_CV = mean(post[:,3])
const PROD = 1.0; const CV = 0.05
const R0, R_MAX, R_TAU = 0.69, 0.90, 14.0
r(a) = R0 + (R_MAX - R0)*(1.0 - exp(-a/R_TAU)); phi(a) = r(a)/R_MAX
const MAX_GEN = 400
@inline incr(D, rng) = PROD*(1.0 + KAPPA*D)*max(0.0, 1.0 + CV*randn(rng))

# Age a mother from seed D0. conserved=true removes phi(a)*D from the mother when she buds the
# age-a daughter (and records the removed amount as that daughter's seed); conserved=false keeps
# the released non-conserved recursion. Returns (lifespan, daughter-seed per age).
function age_mother(D0, rng; conserved::Bool)
    Dc = CRIT_CV > 0 ? D_CRIT*exp(CRIT_CV*randn(rng) - CRIT_CV^2/2) : D_CRIT
    seeds = Float64[]; D = float(D0); a = 0
    while D < Dc && a < MAX_GEN
        push!(seeds, (conserved ? phi(a)*D : phi(a)*D))   # daughter's inherited seed (same share)
        if conserved
            D -= phi(a)*D                                  # strict partition: mother loses the share
        end
        D += incr(D, rng)                                  # autocatalytic production this cycle
        a += 1
    end
    return a, seeds
end

function fold(conserved; N=40_000, seed=20260701)
    rng = MersenneTwister(seed)
    ml = Int[]; ys=0.0; yn=0; os=0.0; on=0
    for _ in 1:N
        L, seeds = age_mother(0.0, rng; conserved=conserved); L == 0 && continue
        push!(ml, L)
        for a in 0:(L-1)
            Ld, _ = age_mother(seeds[a+1], rng; conserved=conserved)
            fr = (L>1) ? a/(L-1) : 0.0
            if fr <= 0.7; ys += Ld; yn += 1 elseif fr >= 0.9; os += Ld; on += 1 end
        end
    end
    return mean(ml), std(ml)/mean(ml), ys/yn, os/on
end

for (nm, cons) in (("non-conserved (released)", false), ("conserved partition", true))
    m, cv, y, o = fold(cons)
    @printf("%-24s : mother RLS mean=%.1f cv=%.2f | daughter young=%.1f old=%.1f fold=%.2fx\n",
            nm, m, cv, y, o, y/o)
end
@printf("(McCormick mother target mean 26.6 cv 0.365 ; Kennedy daughter fold 3.35x)\n")
