#!/usr/bin/env julia
# CC-N (newborn-size law): the population newborn (virgin-daughter) size distribution is the
# age-eroding division asymmetry sampled through the fixed geometric replicative-age law
# P(age=a)=2^{-(a+1)} — a discrete geometric MIXTURE {(2^{-(a+1)}, frac(a)·V*·enlarge(a))} whose
# moments CellSizeControl.newborn_size_law predicts in closed form. This campaign (1) writes the
# predicted comb, (2) checks it against a Monte-Carlo population, (3) shows the cv/skew/ratio are
# V*-independent (scale-free), and (4) decomposes the skew into its erosion / enlargement parts.
# Writes newborn_size_law_comb.csv + newborn_size_law.csv. Run: julia --project=.. gen_newborn_size_law.jl [TARGET]

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)
using CellSizeControl
using Statistics: mean, std

const HERE = @__DIR__

# calibrated demographic parameters (as in the CC-P population campaign)
const A0, AMAX, TAU = 0.32, 0.5, 8.0
const EM, ET = 0.45, 8.0

skewof(x) = (m = mean(x); s = std(x); mean(((x .- m) ./ s) .^ 3))

function simstat(Vstar; target, seed=1)
    pop = simulate_population(
        SizerRule(Vstar);
        target=target,
        enlarge_max=EM,
        enlarge_tau=ET,
        alpha0=A0,
        alpha_max=AMAX,
        tau=TAU,
        seed=seed,
    )
    nb = pop.Vbirth[pop.age .== 0]
    return (nb=nb, mean=mean(nb), cv=std(nb) / mean(nb), skew=skewof(nb), n0=length(nb))
end

function main()
    target = isempty(ARGS) ? 300_000 : parse(Int, ARGS[1])
    Vstar = 60.0

    # ---- (1) the predicted geometric-mixture comb (age → weight, newborn size) ----
    open(joinpath(HERE, "newborn_size_law_comb.csv"), "w") do io
        println(io, "age,weight,newborn_size")
        for a in 0:25
            w = 2.0^(-(a + 1))
            frac = aging_daughter_fraction(a; alpha0=A0, alpha_max=AMAX, tau=TAU)
            s = frac * Vstar * (1 + EM * (1 - exp(-a / ET)))
            println(io, a, ",", round(w; digits=8), ",", round(s; digits=4))
        end
    end

    law = newborn_size_law(; alpha0=A0, alpha_max=AMAX, tau=TAU, enlarge_max=EM,
        enlarge_tau=ET, Vstar=Vstar)
    sim = simstat(Vstar; target=target, seed=1)

    # ---- (2)+(3) analytic vs simulation, and V*-invariance across setpoints ----
    open(joinpath(HERE, "newborn_size_law.csv"), "w") do io
        println(io, "case,Vstar,source,mean,cv,skew,ratio")
        println(io, "calibrated,", Vstar, ",analytic,",
            round(law.mean; digits=4), ",", round(law.cv; digits=4), ",",
            round(law.skew; digits=4), ",", round(law.ratio; digits=5))
        println(io, "calibrated,", Vstar, ",simulation,",
            round(sim.mean; digits=4), ",", round(sim.cv; digits=4), ",",
            round(sim.skew; digits=4), ",", round(sim.mean / (A0 * Vstar); digits=5))
        # scale-free: same demographic params, different set-point → identical cv/skew/ratio
        for V in (30.0, 120.0)
            l = newborn_size_law(; alpha0=A0, alpha_max=AMAX, tau=TAU, enlarge_max=EM,
                enlarge_tau=ET, Vstar=V)
            s = simstat(V; target=target, seed=2)
            println(io, "scalefree,", V, ",analytic,",
                round(l.mean; digits=4), ",", round(l.cv; digits=4), ",",
                round(l.skew; digits=4), ",", round(l.ratio; digits=5))
            println(io, "scalefree,", V, ",simulation,",
                round(s.mean; digits=4), ",", round(s.cv; digits=4), ",",
                round(s.skew; digits=4), ",", round(s.mean / (A0 * V); digits=5))
        end
        # ---- (4) mechanism decomposition of the skew (analytic) ----
        decomp = (
            ("full", EM, AMAX),
            ("erosion_only", 0.0, AMAX),
            ("enlarge_only", EM, A0),
            ("neither", 0.0, A0),
        )
        for (name, em, amax) in decomp
            l = newborn_size_law(; alpha0=A0, alpha_max=amax, tau=TAU, enlarge_max=em,
                enlarge_tau=ET, Vstar=Vstar)
            println(io, name, ",", Vstar, ",analytic,",
                round(l.mean; digits=4), ",", round(l.cv; digits=4), ",",
                round(l.skew; digits=4), ",", round(l.ratio; digits=5))
        end
    end

    # ---- histogram of the simulated newborn sizes rescaled to V*=1 (for the collapse figure) ----
    open(joinpath(HERE, "newborn_size_hist.csv"), "w") do io
        println(io, "Vstar,size_over_Vstar,count")
        for (V, sd) in ((Vstar, sim.nb),)
            r = sd ./ V
            lo, hi = extrema(r)
            nb, edges = 60, nothing
            edges = range(lo, hi; length=nb + 1)
            h = zeros(Int, nb)
            for v in r
                b = clamp(floor(Int, (v - lo) / (hi - lo) * nb) + 1, 1, nb)
                h[b] += 1
            end
            for b in 1:nb
                c = 0.5 * (edges[b] + edges[b + 1])
                println(io, V, ",", round(c; digits=5), ",", h[b])
            end
        end
    end

    # ---- BCa bootstrap CIs on the headline scale-free statistics (ratio/cv/skew) ----
    # The simulation row above is a single-seed point estimate; a bias-corrected-accelerated
    # bootstrap over the newborn sample puts a 95% interval on each, and confirms the analytic
    # closed form lands inside it. Uses CellSizeControl.size_law_ci (vendored ResampleStats).
    ci = size_law_ci(sim.nb; alpha0=A0, Vstar=Vstar, alpha=0.05, nboot=4000, seed=1)
    inside(c, v) = c[1] <= v <= c[3]
    open(joinpath(HERE, "newborn_size_law_ci.csv"), "w") do io
        println(io, "stat,analytic,sim_point,ci_lo,ci_hi,analytic_inside")
        for (name, an, c) in
            (("ratio", law.ratio, ci.ratio), ("cv", law.cv, ci.cv), ("skew", law.skew, ci.skew))
            println(io, name, ",", round(an; digits=5), ",", round(c[2]; digits=5), ",",
                round(c[1]; digits=5), ",", round(c[3]; digits=5), ",", inside(c, an))
        end
    end

    println("CC-N done (target=$target)")
    println("  analytic  : mean=$(round(law.mean;digits=3)) cv=$(round(law.cv;digits=4)) ",
        "skew=$(round(law.skew;digits=4)) ratio=$(round(law.ratio;digits=5))")
    println("  simulation: mean=$(round(sim.mean;digits=3)) cv=$(round(sim.cv;digits=4)) ",
        "skew=$(round(sim.skew;digits=4)) N0=$(sim.n0)")
    println("  BCa 95% CI: ratio ", round.((ci.ratio[1], ci.ratio[3]); digits=4),
        " (analytic ", round(law.ratio; digits=4), inside(ci.ratio, law.ratio) ? " ✓in) " : " ✗out) ",
        "| skew ", round.((ci.skew[1], ci.skew[3]); digits=3))
    println("wrote newborn_size_law_comb.csv + newborn_size_law.csv + newborn_size_hist.csv ",
        "+ newborn_size_law_ci.csv")
    return nothing
end

main()
