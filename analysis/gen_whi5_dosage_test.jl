# CC-6 data: pilot the mechanistic Whi5:SBF switch across a matched range of total Whi5 W and
# emit the emergent commitment size V* = W/c*, NORMALIZED to the 1xWHI5 reference dose. This is
# the model curve for the head-to-head against the Schmoller-2015 Whi5-dosage series (cell size
# vs WHI5 copy number), digitized from Heldt 2018 Fig 3C (see whi5_dosage_data.csv). Because the
# saddle-node threshold c* is a pure number set by the switch shape (independent of W), V* is
# exactly linear in W, so the normalized model curve is the proportionality line through the
# origin. Run: julia --project=. analysis/gen_whi5_dosage_test.jl
using CellSizeControl
using Printf

here = @__DIR__

# Reference dose = 1xWHI5 (the manuscript's canonical W = 18, giving V* ~ 40 fL). Doubling WHI5
# copy number doubles the size-independent Whi5 dose, so whi5_amount (relative to 1x) maps to
# W = W_ref * amount. We sweep a fine grid spanning the experimental range (whi5Delta ~ 0 up to
# 2xWHI5) so the emitted curve can be drawn as a continuous proportionality line.
const W_REF = 18.0
Vstar(W) = setpoint_volume(Whi5SBFSwitch(Float64(W)))
const VSTAR_REF = Vstar(W_REF)

amounts = collect(range(0.0, 2.6, length=53))   # Whi5 dose relative to 1xWHI5
open(joinpath(here, "whi5_dosage_model.csv"), "w") do io
    println(io, "whi5_amount,W,Vstar_fL,Vstar_norm")
    for a in amounts
        W = W_REF * a
        # V* = 0 at W = 0 (no inhibitor, the switch has no OFF branch to lose); use the exact
        # linear law there to avoid the degenerate constructor, and the package elsewhere.
        v = a == 0.0 ? 0.0 : Vstar(W)
        @printf(io, "%.4f,%.4f,%.5f,%.6f\n", a, W, v, v / VSTAR_REF)
    end
end

# Report the proportionality check at the experimental doses (1x and 2x).
v1, v2 = Vstar(W_REF), Vstar(2 * W_REF)
@printf("model V*(1xWHI5) = %.3f fL, V*(2xWHI5) = %.3f fL, fold = %.4f (proportional law -> 2.0)\n",
        v1, v2, v2 / v1)
println("wrote whi5_dosage_model.csv")
