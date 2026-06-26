# CC-5 data: the mechanistic Whi5:SBF bistable Start switch. Emits (a) the full bifurcation
# S-curve SBF-activity x vs Whi5 concentration c — including the UNSTABLE middle branch, from
# the exact parametric steady state — and (b) the emergent set-point V* vs total Whi5 W, which
# falls exactly on the phenomenological inhibitor-dilution line V* = W/θ. Run: julia --project=. .
using CellSizeControl
using Printf

# default switch parameters (match the package constructor)
const BETA, GAMMA, KE, Q, KX, P = 1.0, 1.0, 0.30, 4.0, 0.40, 4.0

# exact steady state, parameterized by x: γx = β/(1+(e/Ke)^q) with e = c/(1+(x/Kx)^p).
# Solve for c(x):  e = Ke·(β/(γx) − 1)^(1/q),  c = e·(1 + (x/Kx)^p),  for x ∈ (0, β/γ).
c_of_x(x) = KE * (BETA / (GAMMA * x) - 1)^(1 / Q) * (1 + (x / KX)^P)

here = @__DIR__

# ---- (a) bifurcation curve + branch classification (off / unstable / on) ----
# x grid with log-refined tails: near x→0 the OFF branch runs to high c, near x→(β/γ) the ON
# branch runs to c→0, so both stable branches must be sampled into their asymptotes.
xmax = BETA / GAMMA
lo_tail = exp10.(range(-9, -1.2, length=300))                 # x ≈ 1e-9 … 0.063
mid = collect(range(0.05, stop=0.95 * xmax, length=900))
hi_tail = xmax .- exp10.(range(-1.2, -9, length=300))         # x ≈ (β/γ − 0.063) … −1e-9
xs = sort(unique(vcat(lo_tail, mid, hi_tail)))
cs = [c_of_x(x) for x in xs]
# the two folds = local extrema of c(x); below the lower (c_F1) the OFF state is gone (Start).
dc = diff(cs)
fold_idx = [i for i in 2:length(dc) if sign(dc[i]) != sign(dc[i - 1])]
i_min = fold_idx[1]                       # lower fold: local min of c  → c* (OFF disappears)
i_max = length(fold_idx) >= 2 ? fold_idx[2] : length(xs)
cstar = cs[i_min]

open(joinpath(here, "whi5_sbf_bifurcation.csv"), "w") do io
    println(io, "x,c,branch")
    for (i, (x, c)) in enumerate(zip(xs, cs))
        branch = i <= i_min ? "off" : (i <= i_max ? "unstable" : "on")
        @printf(io, "%.6f,%.6f,%s\n", x, c, branch)
    end
end

# ---- (b) emergent set-point V* vs W, mechanistic vs the law V* = W/θ ----
θ = whi5_sbf_threshold(Whi5SBFSwitch(18.0))   # emergent threshold (W-independent)
open(joinpath(here, "whi5_sbf_setpoint.csv"), "w") do io
    println(io, "W,Vstar_mech,Vstar_law,theta")
    for W in (9.0, 18.0, 36.0, 72.0, 144.0)
        sw = Whi5SBFSwitch(W)
        @printf(io, "%.3f,%.5f,%.5f,%.6f\n", W, setpoint_volume(sw), W / θ, whi5_sbf_threshold(sw))
    end
end

@printf("c* (OFF saddle-node) = %.5f   θ from package = %.5f   match=%s\n",
        cstar, θ, isapprox(cstar, θ; atol=3e-3))
@printf("hysteresis window: c in [%.4f, %.4f]\n", cs[i_min], cs[i_max])
println("wrote whi5_sbf_bifurcation.csv + whi5_sbf_setpoint.csv")
