# CC-3 data: the size-control phase diagram. (a) the recovered Vd-vs-Vb slope landscape over
# (control strength α, division asymmetry f) with the homeostasis boundary α·f = 1; (b) where
# the Soifer-Amir discriminator MISCLASSIFIES under measurement noise — the fraction of
# finite-sample lineages whose recovered slope lands in the wrong sizer/adder/timer bin, over
# (α, noise cv). Built on the linear size-control map Vd = αVb + β. Run: julia --project=. .
using CellSizeControl
using Statistics: mean
using Printf

here = @__DIR__

# ---- (a) slope landscape over (α, f) + homeostasis indicator ----
alphas_a = collect(0.0:0.05:2.2)
fracs = collect(0.30:0.01:0.62)
open(joinpath(here, "phase_alpha_f.csv"), "w") do io
    println(io, "alpha,f,slope,logratio")
    for f in fracs, α in alphas_a
        s = simulate_lineage(
            LinearSizeControl(α, 20.0); V0=20.0, n=300, cv=0.04, daughter_fraction=f, seed=7
        )
        slope = size_control_slope(s.Vb, s.Vd)
        logratio = log10(max(last(s.Vb), 1e-12) / first(s.Vb))   # >0 grows, <0 collapses
        @printf(io, "%.3f,%.3f,%.4f,%.4f\n", α, f, slope, logratio)
    end
end

# ---- (b) discriminator misclassification under noise over (α, cv) ----
# finite-sample experiment: n lineage points, R independent replicates; a replicate is
# "correctly classified" iff classify_control(recovered slope) == classify_control(true α).
# cv starts at 0.02: a noiseless perfect sizer makes identical cells (zero Vb variance → the
# regression slope is undefined), so a slope can only be measured with real measurement noise.
alphas_b = collect(0.0:0.1:2.0)
cvs = collect(0.02:0.02:0.30)
n_obs, R = 80, 300
open(joinpath(here, "phase_misclass.csv"), "w") do io
    println(io, "alpha,cv,misclass")
    for cv in cvs, α in alphas_b
        truth = classify_control(α)
        wrong = 0
        for r in 1:R
            # V0=10 is a realistic newborn (below the β=20 set-point), so cells grow INTO
            # division without tripping the "division can't precede growth" floor.
            s = simulate_lineage(
                LinearSizeControl(α, 20.0); V0=10.0, n=n_obs, cv=cv, daughter_fraction=0.5,
                seed=1000 * r + 1,
            )
            classify_control(size_control_slope(s.Vb, s.Vd)) === truth || (wrong += 1)
        end
        @printf(io, "%.3f,%.3f,%.4f\n", α, cv, wrong / R)
    end
end

println("wrote phase_alpha_f.csv (", length(alphas_a) * length(fracs), " cells) + ",
        "phase_misclass.csv (", length(alphas_b) * length(cvs), " cells, R=$R each)")
