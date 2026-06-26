# Cross-validation: recompute the model's value for each independent published benchmark from
# the package primitives (one parameterization, no per-target refitting) and tabulate it against
# the literature value + source. Emits crossval.csv. Run: julia --project=. gen_crossval.jl
using CellSizeControl
using Statistics: mean, std
using Printf

here = @__DIR__

# --- size-control slopes (Soifer-Amir 2016: timer 2 / adder 1 / sizer 0) ---
slope(rule; kw...) = (s = simulate_lineage(rule; n=600, kw...); size_control_slope(s.Vb, s.Vd))
sizer_slope = slope(SizerRule(2.0); seed=3)
adder_slope = slope(AdderRule(1.0); seed=2)
timer_slope = slope(TimerRule(2.0); seed=1)

# --- two-step G1 mother/daughter (Di Talia 2007: ~19 / ~45 min) ---
Vstar, T_cln2, tau_bud = 36.0, 19.0, 52.0
mother_G1 = cell_cycle(Vstar + 5; Vstar=Vstar, T_cln2=T_cln2, tau_bud=tau_bud).G1
Vdau = cell_cycle(Vstar; Vstar=Vstar, T_cln2=T_cln2, tau_bud=tau_bud).Vdaughter
daughter_G1 = cell_cycle(Vdau; Vstar=Vstar, T_cln2=T_cln2, tau_bud=tau_bud).G1

# --- the inhibitor-dilution law from the bistable switch (V* = W/theta) ---
sw = Whi5SBFSwitch(18.0)
vstar_ratio = setpoint_volume(sw) / (sw.W / whi5_sbf_threshold(sw))   # = 1 exactly (consistency)

# --- replicative lifespan (data-calibrated ABC posterior-predictive vs McCormick 2015) ---
pp = [parse(Int, l) for (i, l) in enumerate(eachline(joinpath(here, "rls_abc_predictive.csv"))) if i > 1]
rls_mean = mean(pp)
rls_cv = std(pp) / mean(pp)

# metric, model, reference, ref_lo, ref_hi, unit, source
rows = [
    ("Sizer slope", sizer_slope, 0.0, -0.2, 0.2, "", "Soifer-Amir 2016"),
    ("Adder slope", adder_slope, 1.0, 0.8, 1.2, "", "Soifer-Amir 2016"),
    ("Timer slope", timer_slope, 2.0, 1.8, 2.2, "", "Soifer-Amir 2016"),
    ("Mother G1", mother_G1, 19.0, 17.0, 21.0, "min", "Di Talia 2007"),
    ("Daughter G1", daughter_G1, 45.5, 40.0, 50.0, "min", "Di Talia 2007"),
    ("RLS mean", rls_mean, 26.6, 25.0, 28.0, "div", "McCormick 2015"),
    ("RLS CV", rls_cv, 0.365, 0.34, 0.39, "", "McCormick 2015"),
]

open(joinpath(here, "crossval.csv"), "w") do io
    println(io, "metric,model,reference,ref_lo,ref_hi,unit,source")
    for (m, mod, ref, lo, hi, u, src) in rows
        @printf(io, "%s,%.4f,%.4f,%.4f,%.4f,%s,%s\n", m, mod, ref, lo, hi, u, src)
    end
end

println("cross-validation (model vs published, one parameterization):")
for (m, mod, ref, lo, hi, u, src) in rows
    @printf("  %-13s model %7.3f  vs  ref %6.3f %-3s  [%s]\n", m, mod, ref, u, src)
end
@printf("Whi5:SBF switch reproduces V*=W/theta to ratio %.6f (consistency check)\n", vstar_ratio)
println("wrote crossval.csv")
