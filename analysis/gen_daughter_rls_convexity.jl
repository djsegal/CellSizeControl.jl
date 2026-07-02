# SHAPE of the daughter-RLS-vs-maternal-age prediction: is it a straight line, or convex?
#
# The two-bucket Kennedy fold (daughters of young vs old mothers) only tests the endpoints.
# A sharper, falsifiable claim is the SHAPE of the full curve. Under the single age-eroding
# asymmetry, a daughter of an age-a mother inherits a share phi(a) = alpha(a)/alpha_max of the
# mother's autocatalytic damage D_m(a) and begins partway up her own damage trajectory, so her
# emergent lifespan is shortened. This script builds the population fraction-binned curve from
# the PACKAGE primitives (damage_trajectory for D_m(a); replicative_lifespan(; D0=...) for the
# daughter's seeded lifespan -- NOTHING is refit) and characterizes its convexity:
#   * chord test        -- a convex curve lies BELOW the straight line joining its endpoints
#   * quadratic fit      -- convex <=> positive x^2 coefficient; linear-vs-quadratic residual
#   * steepest-decline location -- where the slope is most negative (early/mid vs a late cliff)
#
# Result (reproducible): the population curve is monotone-decreasing and CONVEX -- a decelerating
# decline, steepest at young/mid maternal age and flattening toward old age -- NOT a linear
# decline and NOT a late cliff. Mechanistically the convexity is emergent: a single noiseless
# mother's daughter curve is ~linear; averaging over the lifespan + threshold heterogeneity
# (crit_cv) bends the population shape convex, because daughter RLS is a concave (~log) function
# of inherited damage and so is most sensitive to D_m at low damage (young mothers).
#
# Run: julia --project=. analysis/gen_daughter_rls_convexity.jl

using CellSizeControl
using Statistics: mean
using Printf

here = @__DIR__

# passive volume-proportional inheritance share phi(a) = alpha(a)/alpha_max, alpha the size-face
# asymmetry fixed by the Johnston-1966 / Yang-2011 daughter-size increase (r0,r_max,r_tau).
const R0, R_MAX, R_TAU = 0.69, 0.90, 14.0
phi(a) = aging_daughter_fraction(a; alpha0=R0, alpha_max=R_MAX, tau=R_TAU) / R_MAX

# population fraction-binned daughter-RLS curve, built entirely from package functions
function daughter_curve(; N::Int=40_000, nbin::Int=20, mseed0::Int=0, dseed0::Int=10_000_000)
    binsum = zeros(nbin); binN = zeros(Int, nbin); dc = dseed0
    for m in 1:N
        traj = damage_trajectory(; seed=mseed0 + m)   # the mother's D_m(a) series
        L = length(traj); L < 2 && continue
        for a in 0:(L - 1)
            Ld = replicative_lifespan(; D0=phi(a) * traj[a + 1], seed=(dc += 1))
            b = clamp(floor(Int, (a / (L - 1)) * nbin) + 1, 1, nbin)
            binsum[b] += Ld; binN[b] += 1
        end
    end
    keep = binN .> 0
    x = [(b - 0.5) / nbin for b in 1:nbin][keep]
    y = (binsum ./ max.(binN, 1))[keep]
    return x, y, binN[keep]
end

# least-squares polynomial fit (degree deg) -> coefficients highest-power first
function polyfit(x, y, deg)
    A = reduce(hcat, [x .^ p for p in deg:-1:0])
    return A \ y
end
polyval(c, x) = sum(c[i] * x .^ (length(c) - i) for i in eachindex(c))
rms(r) = sqrt(mean(abs2, r))

x, y, n = daughter_curve()

# --- convexity diagnostics -------------------------------------------------------------------
chord = y[1] .+ (y[end] - y[1]) .* (x .- x[1]) ./ (x[end] - x[1])
dev = y .- chord                                   # convex => dev <= 0 everywhere
below = all(dev .<= 1e-9)
max_gap = -minimum(dev)                            # peak sag below the chord (convexity depth)
frac_at_max_gap = x[argmin(dev)]

qc = polyfit(x, y, 2)                              # [a, b, c]  (a>0 => convex)
lc = polyfit(x, y, 1)
rms_lin = rms(y .- polyval(lc, x))
rms_quad = rms(y .- polyval(qc, x))

# discrete slope: steepest-decline location (early/mid vs a late cliff)
slope = diff(y) ./ diff(x)
xmid = (x[1:(end - 1)] .+ x[2:end]) ./ 2
steep_frac = xmid[argmin(slope)]                   # most negative slope
shallow_frac = xmid[argmax(slope)]

# second differences (discrete curvature) sign count
d2 = diff(diff(y))

@printf("\n=== daughter-RLS-vs-maternal-age SHAPE (package-native, nothing refit) ===\n")
@printf("bins=%d  young (frac %.2f) RLS=%.2f   old (frac %.2f) RLS=%.2f   fold=%.2fx\n",
        length(x), x[1], y[1], x[end], y[end], y[1] / y[end])
@printf("monotone decreasing         : %s\n", all(diff(y) .< 0))
@printf("convex (curve below chord)  : %s   (peak sag %.2f div at maternal-age frac %.2f)\n",
        below, max_gap, frac_at_max_gap)
@printf("quadratic x^2 coefficient   : %+.2f   (>0 => convex)\n", qc[1])
@printf("2nd-difference signs        : #pos=%d  #neg=%d   (mostly + => convex)\n",
        count(>(0), d2), count(<(0), d2))
@printf("linear vs quadratic RMS     : %.3f  ->  %.3f div  (%.0f%% residual removed by curvature)\n",
        rms_lin, rms_quad, 100 * (1 - rms_quad / rms_lin))
@printf("steepest decline at frac    : %.2f   shallowest at frac : %.2f   (early/mid, not a late cliff)\n",
        steep_frac, shallow_frac)

open(joinpath(here, "daughter_rls_convexity.csv"), "w") do io
    println(io, "frac_mid,daughter_rls_mean,chord,dev_below_chord,n")
    for i in eachindex(x)
        @printf(io, "%.4f,%.4f,%.4f,%.4f,%d\n", x[i], y[i], chord[i], dev[i], n[i])
    end
end
open(joinpath(here, "daughter_rls_convexity_summary.csv"), "w") do io
    println(io, "quantity,value")
    @printf(io, "young_rls,%.4f\n", y[1])
    @printf(io, "old_rls,%.4f\n", y[end])
    @printf(io, "fold,%.4f\n", y[1] / y[end])
    @printf(io, "convex_below_chord,%d\n", below ? 1 : 0)
    @printf(io, "peak_sag_divisions,%.4f\n", max_gap)
    @printf(io, "peak_sag_frac,%.4f\n", frac_at_max_gap)
    @printf(io, "quad_x2_coeff,%.4f\n", qc[1])
    @printf(io, "rms_linear,%.4f\n", rms_lin)
    @printf(io, "rms_quadratic,%.4f\n", rms_quad)
    @printf(io, "steepest_decline_frac,%.4f\n", steep_frac)
end

println("\nwrote daughter_rls_convexity.csv, daughter_rls_convexity_summary.csv")
