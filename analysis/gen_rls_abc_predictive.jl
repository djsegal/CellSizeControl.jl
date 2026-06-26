# CC-4 (cont.): the proper Bayesian posterior-PREDICTIVE RLS distribution — pool lifespan
# samples over draws FROM the posterior chain (integrating parameter uncertainty), rather than
# simulating only at the posterior mean (which sits off the curved (D_crit,kappa) ridge). Reads
# rls_abc_posterior.csv; writes rls_abc_predictive.csv. Run: julia --project=. .
using CellSizeControl
using Statistics: mean, std, cor
using Random: MersenneTwister
using Printf

here = @__DIR__

# read the posterior draws
Ds = Float64[]; ks = Float64[]; cs = Float64[]
for (i, line) in enumerate(eachline(joinpath(here, "rls_abc_posterior.csv")))
    i == 1 && continue
    p = split(line, ',')
    push!(Ds, parse(Float64, p[1])); push!(ks, parse(Float64, p[2])); push!(cs, parse(Float64, p[3]))
end
n = length(Ds)
@printf("posterior: %d draws.  (D_crit,kappa) ridge corr = %.3f\n", n, cor(Ds, ks))

# pooled posterior-predictive: for each of NP random posterior draws, simulate K lifespans
rng = MersenneTwister(424242)
NP, K = 4000, 12
pp = Int[]
for _ in 1:NP
    j = rand(rng, 1:n)
    ls = lifespan_distribution(K; seed0=rand(rng, 1:10^7), D_crit=Ds[j], kappa=ks[j],
                               crit_cv=cs[j], production=1.0, cv=0.05, max_gen=400)
    append!(pp, ls)
end
open(joinpath(here, "rls_abc_predictive.csv"), "w") do io
    println(io, "rls")
    for v in pp
        println(io, v)
    end
end
@printf("pooled posterior-predictive RLS (%d samples): mean=%.2f sd=%.2f cv=%.3f  (McCormick 26.6 / 9.7 / 0.365)\n",
        length(pp), mean(pp), std(pp), std(pp) / mean(pp))
println("wrote rls_abc_predictive.csv")
