# ===========================================================================
# VENDORED KERNEL — provenance
# ---------------------------------------------------------------------------
# Source : science-space/tools/numerics/ResampleStats.jl
# Commit : 07d7e7f304c47158a16eb0109645ddf4fa37d670 (2026-07-03)
# Vendored into CellSizeControl to put bootstrap/BCa confidence intervals on
# the package's headline Monte-Carlo size-law statistics (newborn-size ratio,
# CV, skew), which previously shipped as bare single-seed point estimates.
#
# Only cosmetic adaptation vs the upstream file: the two blanket `using Random`
# / `using Statistics` lines are narrowed to explicit `using X: name` imports
# and `Random.default_rng()` is spelled `default_rng()`, so the module passes
# CellSizeControl's ExplicitImports / Aqua release gates. The numerics are
# byte-identical to upstream; keep this copy in sync when the kernel changes.
# ===========================================================================
"""
    ResampleStats

Dependency-light resampling and weighted-summary primitives.

Self-contained, stdlib-only (`Random`, `Statistics`), usable by `include`:

    include("ResampleStats.jl")
    using .ResampleStats

Everything is deterministic: any routine that resamples takes an explicit
`rng::AbstractRNG`, so results reproduce bit-for-bit under a fixed seed.

Public interface:
  - `weighted_quantile(x, w, q) -> Float64`
      The `q`-quantile (`q ∈ [0,1]`) of values `x` with nonnegative weights
      `w`, by linear interpolation on the centered cumulative weight. Reduces
      **exactly** to the plain type-7 quantile (`Statistics.quantile` default)
      when all weights are equal.
  - `weighted_mean(x, w) -> Float64`
      `Σ wᵢ xᵢ / Σ wᵢ`.
  - `weighted_std(x, w; corrected=false) -> Float64`
      Weighted standard deviation. Default is the *uncorrected* (population)
      form `sqrt(Σ wᵢ (xᵢ-μ)² / Σ wᵢ)`; with equal weights this equals
      `Statistics.std(x; corrected=false)`. Pass `corrected=true` for the
      reliability-weight (frequency-style) `Σw/(Σw-1)` correction, which with
      equal weights equals `Statistics.std(x)`.
  - `bootstrap_ci(data, statistic; nboot=2000, alpha=0.05, rng) -> (lo, point, hi)`
      Percentile bootstrap CI of `statistic(sample)::Float64` over `nboot`
      resamples-with-replacement; `point = statistic(data)`.
  - `bca_ci(data, statistic; nboot=2000, alpha=0.05, rng) -> (lo, point, hi)`
      Bias-corrected and accelerated (BCa) bootstrap CI. Same bootstrap
      distribution as `bootstrap_ci`, but the percentiles are shifted by a
      bias-correction `z0` (from the fraction of `θ*` below `point`) and an
      acceleration `a` (from the jackknife skewness), giving a
      transformation-respecting, skew-aware interval that reduces to the plain
      percentile interval when `z0 ≈ 0` and `a ≈ 0`.
  - `jackknife(data, statistic) -> (estimate, bias, se)`
      Leave-one-out jackknife: bias-corrected `estimate`, `bias`, and standard
      error `se`.
  - `ecdf(x) -> Function`
      Empirical CDF as a callable `F(t)` = fraction of `x` that is `≤ t`.
      O(log n) per query via a sorted copy.
  - `kde(x; bandwidth=nothing, npoints=256) -> (grid, density)`
      Gaussian kernel density estimate on an evenly spaced `grid`; default
      bandwidth by Silverman's rule. `density` integrates to ≈ 1 over `grid`.
"""
module ResampleStats

using Random: AbstractRNG, default_rng
using Statistics: std

export weighted_quantile, weighted_mean, weighted_std,
    bootstrap_ci, bca_ci, jackknife, ecdf, kde

# ---------------------------------------------------------------------------
# Weighted summaries
# ---------------------------------------------------------------------------

"""
    weighted_quantile(x, w, q) -> Float64

`q`-quantile of `x` weighted by nonnegative `w`, `q ∈ [0,1]`.

Order statistics are placed at the *centered* cumulative-weight positions
`gᵢ = (Σ_{j≤i} wⱼ) - wᵢ/2`, rescaled affinely so the minimum sits at 0 and the
maximum at 1; the quantile is a linear interpolation of `x` against that
position. With equal weights this collapses to `(i-1)/(n-1)`, i.e. the type-7
quantile used by `Statistics.quantile`. Up-weighting a point drags the
quantile toward it, as expected.
"""
function weighted_quantile(x::AbstractVector{<:Real}, w::AbstractVector{<:Real},
                           q::Real)
    n = length(x)
    n == length(w) || throw(ArgumentError("x and w must have equal length"))
    n == 0 && throw(ArgumentError("empty input"))
    (0.0 <= q <= 1.0) || throw(ArgumentError("q must be in [0,1]"))
    any(<(0), w) && throw(ArgumentError("weights must be nonnegative"))
    W = sum(w)
    W > 0 || throw(ArgumentError("weights must sum to a positive value"))

    p = sortperm(x)
    xs = float.(x[p])
    ws = float.(w[p])
    n == 1 && return xs[1]

    # centered cumulative weight of each order statistic
    g = similar(ws)
    acc = 0.0
    @inbounds for i in 1:n
        acc += ws[i]
        g[i] = acc - ws[i] / 2
    end
    g1 = g[1]
    gn = g[n]
    denom = gn - g1
    # degenerate: all weight on coincident endpoints
    denom > 0 || return xs[n]

    target = q
    # clamp outside the interpolation range
    P1 = 0.0
    Pn = 1.0
    Pfirst = (g[1] - g1) / denom            # == 0
    Plast = (g[n] - g1) / denom             # == 1
    target <= Pfirst && return xs[1]
    target >= Plast && return xs[n]

    @inbounds for i in 1:(n - 1)
        Pi = (g[i] - g1) / denom
        Pj = (g[i + 1] - g1) / denom
        if target <= Pj
            Pj == Pi && return xs[i]
            t = (target - Pi) / (Pj - Pi)
            return xs[i] + t * (xs[i + 1] - xs[i])
        end
    end
    return xs[n]
end

"""
    weighted_mean(x, w) -> Float64

Weighted arithmetic mean `Σ wᵢ xᵢ / Σ wᵢ` (nonnegative weights).
"""
function weighted_mean(x::AbstractVector{<:Real}, w::AbstractVector{<:Real})
    length(x) == length(w) || throw(ArgumentError("x and w must have equal length"))
    !isempty(x) || throw(ArgumentError("empty input"))
    any(<(0), w) && throw(ArgumentError("weights must be nonnegative"))
    W = sum(w)
    W > 0 || throw(ArgumentError("weights must sum to a positive value"))
    s = 0.0
    @inbounds for i in eachindex(x)
        s += w[i] * x[i]
    end
    return s / W
end

"""
    weighted_std(x, w; corrected=false) -> Float64

Weighted standard deviation about the weighted mean.

Uncorrected (default): `sqrt(Σ wᵢ (xᵢ-μ)² / Σ wᵢ)` — equals
`Statistics.std(x; corrected=false)` at equal weights. With `corrected=true`
the frequency-weight correction `Σw/(Σw−1)` is applied to the denominator,
which at equal weights reproduces `Statistics.std(x)` (the n−1 correction).
"""
function weighted_std(x::AbstractVector{<:Real}, w::AbstractVector{<:Real};
                      corrected::Bool=false)
    length(x) == length(w) || throw(ArgumentError("x and w must have equal length"))
    !isempty(x) || throw(ArgumentError("empty input"))
    any(<(0), w) && throw(ArgumentError("weights must be nonnegative"))
    W = sum(w)
    W > 0 || throw(ArgumentError("weights must sum to a positive value"))
    μ = weighted_mean(x, w)
    ss = 0.0
    @inbounds for i in eachindex(x)
        d = x[i] - μ
        ss += w[i] * d * d
    end
    if corrected
        W > 1 || throw(ArgumentError("corrected std needs sum(w) > 1"))
        return sqrt(ss / (W - 1))
    else
        return sqrt(ss / W)
    end
end

# ---------------------------------------------------------------------------
# Plain type-7 quantile (internal; used by the bootstrap)
# ---------------------------------------------------------------------------

# type-7 quantile of an unsorted vector v at probability q
function _quantile7(v::AbstractVector{<:Real}, q::Real)
    n = length(v)
    n == 0 && throw(ArgumentError("empty input"))
    s = sort(float.(v))
    n == 1 && return s[1]
    h = (n - 1) * q
    lo = floor(Int, h)
    frac = h - lo
    lo0 = lo + 1                    # 1-based lower index
    lo0 >= n && return s[n]
    return s[lo0] + frac * (s[lo0 + 1] - s[lo0])
end

# ---------------------------------------------------------------------------
# Bootstrap
# ---------------------------------------------------------------------------

"""
    bootstrap_ci(data, statistic; nboot=2000, alpha=0.05, rng) -> (lo, point, hi)

Percentile bootstrap confidence interval. Draws `nboot` resamples of size
`length(data)` with replacement from `data` (using `rng`), evaluates
`statistic(sample)::Float64` on each, and returns the `alpha/2` and
`1-alpha/2` percentiles of that bootstrap distribution as `(lo, hi)`, together
with `point = statistic(data)` computed on the full sample.
"""
function bootstrap_ci(data::AbstractVector, statistic;
                      nboot::Integer=2000, alpha::Real=0.05,
                      rng::AbstractRNG=default_rng())
    n = length(data)
    n > 0 || throw(ArgumentError("empty data"))
    nboot > 0 || throw(ArgumentError("nboot must be positive"))
    (0.0 < alpha < 1.0) || throw(ArgumentError("alpha must be in (0,1)"))

    point = Float64(statistic(data))
    boots = Vector{Float64}(undef, nboot)
    sample = similar(data, n)
    @inbounds for b in 1:nboot
        for i in 1:n
            sample[i] = data[rand(rng, 1:n)]
        end
        boots[b] = Float64(statistic(sample))
    end
    lo = _quantile7(boots, alpha / 2)
    hi = _quantile7(boots, 1 - alpha / 2)
    return (lo, point, hi)
end

# ---------------------------------------------------------------------------
# Jackknife
# ---------------------------------------------------------------------------

"""
    jackknife(data, statistic) -> (estimate, bias, se)

Leave-one-out jackknife. With `θ̂ = statistic(data)` and `θ̂₍ᵢ₎` the statistic on
`data` with element `i` removed:

  - `estimate` = bias-corrected `n·θ̂ − (n−1)·mean(θ̂₍ᵢ₎)`,
  - `bias`     = `(n−1)·(mean(θ̂₍ᵢ₎) − θ̂)`,
  - `se`       = `sqrt((n−1)/n · Σ (θ̂₍ᵢ₎ − mean(θ̂₍ᵢ₎))²)`.

For a linear statistic like the mean the bias is exactly 0 and `se = std(data)/√n`.
"""
function jackknife(data::AbstractVector, statistic)
    n = length(data)
    n >= 2 || throw(ArgumentError("jackknife needs at least 2 observations"))
    full = Float64(statistic(data))
    loo = Vector{Float64}(undef, n)
    keep = Vector{eltype(data)}(undef, n - 1)
    @inbounds for i in 1:n
        k = 1
        for j in 1:n
            j == i && continue
            keep[k] = data[j]
            k += 1
        end
        loo[i] = Float64(statistic(keep))
    end
    m = sum(loo) / n
    bias = (n - 1) * (m - full)
    estimate = n * full - (n - 1) * m
    ss = 0.0
    @inbounds for i in 1:n
        d = loo[i] - m
        ss += d * d
    end
    se = sqrt((n - 1) / n * ss)
    return (estimate, bias, se)
end

# ---------------------------------------------------------------------------
# BCa bootstrap  (bias-corrected and accelerated)
# ---------------------------------------------------------------------------
#
# Self-contained standard-normal CDF Φ and quantile Φ⁻¹ (probit), so the module
# stays stdlib-only (no dependency on SpecialFns.jl). The BCa construction needs
# both: Φ⁻¹ to turn the bias fraction into z0 and to map the target coverage
# levels, and Φ to turn the adjusted z-scores back into percentiles.

# erf via Abramowitz & Stegun 7.1.26 (|error| ≤ 1.5e-7); odd extension for x<0.
function _erf_as(x::Float64)
    s = sign(x)
    z = abs(x)
    t = 1.0 / (1.0 + 0.3275911 * z)
    y = 1.0 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t -
                0.284496736) * t + 0.254829592) * t * exp(-z * z)
    return s * y
end

# Standard-normal CDF Φ(x) = ½(1 + erf(x/√2)).
_norm_cdf(x::Real) = 0.5 * (1.0 + _erf_as(float(Float64(x)) / sqrt(2.0)))

# Standard-normal quantile Φ⁻¹(p) (probit) via Acklam's rational approximation
# (relative error ≈ 1.15e-9), refined by one Halley step against `_norm_cdf`.
# p = 0 / 1 map to ∓Inf.
function _norm_quantile(p::Real)
    pp = float(Float64(p))
    (0.0 <= pp <= 1.0) || throw(DomainError(pp, "_norm_quantile needs 0 ≤ p ≤ 1"))
    pp == 0.0 && return -Inf
    pp == 1.0 && return Inf

    a = (-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00)
    b = (-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01)
    c = (-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00)
    d = (7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00)

    plow = 0.02425
    phigh = 1.0 - plow
    local z::Float64
    if pp < plow
        q = sqrt(-2.0 * log(pp))
        z = (((((c[1]*q + c[2])*q + c[3])*q + c[4])*q + c[5])*q + c[6]) /
            ((((d[1]*q + d[2])*q + d[3])*q + d[4])*q + 1.0)
    elseif pp <= phigh
        q = pp - 0.5
        r = q * q
        z = (((((a[1]*r + a[2])*r + a[3])*r + a[4])*r + a[5])*r + a[6]) * q /
            (((((b[1]*r + b[2])*r + b[3])*r + b[4])*r + b[5])*r + 1.0)
    else
        q = sqrt(-2.0 * log(1.0 - pp))
        z = -(((((c[1]*q + c[2])*q + c[3])*q + c[4])*q + c[5])*q + c[6]) /
              ((((d[1]*q + d[2])*q + d[3])*q + d[4])*q + 1.0)
    end
    # one Halley refinement: e = Φ(z) - p, u = e·√(2π)·e^{z²/2}
    e = _norm_cdf(z) - pp
    u = e * sqrt(2.0 * pi) * exp(0.5 * z * z)
    z = z - u / (1.0 + 0.5 * z * u)
    return z
end

# Acceleration `a` from the leave-one-out jackknife skewness:
#   a = Σ(θ̄ − θ₋ᵢ)³ / (6 · (Σ(θ̄ − θ₋ᵢ)²)^{3/2}).
# Deterministic (no rng); a = 0 exactly for a symmetric jackknife (e.g. the
# mean of a symmetric dataset), and takes the sign of the statistic's skew.
function _bca_acceleration(data::AbstractVector, statistic)
    n = length(data)
    n >= 2 || throw(ArgumentError("BCa acceleration needs at least 2 observations"))
    loo = Vector{Float64}(undef, n)
    keep = Vector{eltype(data)}(undef, n - 1)
    @inbounds for i in 1:n
        k = 1
        for j in 1:n
            j == i && continue
            keep[k] = data[j]
            k += 1
        end
        loo[i] = Float64(statistic(keep))
    end
    m = sum(loo) / n
    num = 0.0
    den = 0.0
    @inbounds for i in 1:n
        d = m - loo[i]
        num += d * d * d
        den += d * d
    end
    den == 0.0 && return 0.0        # constant statistic: no acceleration
    return num / (6.0 * den^1.5)
end

# Bias correction z0 = Φ⁻¹(fraction of θ* strictly below `point`). The fraction
# is clamped away from 0/1 (so z0 stays finite when point sits at/near the tail
# of the bootstrap distribution).
function _bca_z0(boots::AbstractVector{<:Real}, point::Real)
    nb = length(boots)
    nb > 0 || throw(ArgumentError("empty bootstrap distribution"))
    cnt = 0
    @inbounds for b in boots
        b < point && (cnt += 1)
    end
    frac = cnt / nb
    lo = 1.0 / (2.0 * nb)           # half a resample's worth of probability
    frac = clamp(frac, lo, 1.0 - lo)
    return _norm_quantile(frac)
end

"""
    bca_ci(data, statistic; nboot=2000, alpha=0.05, rng) -> (lo, point, hi)

Bias-corrected and accelerated (BCa) bootstrap confidence interval.

Draws `nboot` resamples of size `length(data)` with replacement (using `rng`)
and evaluates `θ*ᵦ = statistic(sample)` on each, exactly as `bootstrap_ci` does.
Instead of reading the raw `alpha/2` and `1-alpha/2` percentiles, BCa reads
*adjusted* percentiles that correct for median bias and skew:

  - **bias correction** `z0 = Φ⁻¹(#{θ* < point} / nboot)` — zero when the
    bootstrap distribution is centered on `point`;
  - **acceleration** `a = Σ(θ̄ − θ₋ᵢ)³ / (6 (Σ(θ̄ − θ₋ᵢ)²)^{3/2})` from the
    leave-one-out jackknife values `θ₋ᵢ` (`θ̄` their mean) — zero for a
    symmetric statistic;
  - the two coverage endpoints `z_α ∈ {Φ⁻¹(alpha/2), Φ⁻¹(1-alpha/2)}` are mapped
    to `α = Φ(z0 + (z0 + z_α) / (1 − a (z0 + z_α)))`, and `(lo, hi)` are those
    (type-7) percentiles of the bootstrap distribution.

`point = statistic(data)` on the full sample. When `z0 ≈ 0` and `a ≈ 0` the
adjusted percentiles collapse back to `alpha/2` and `1-alpha/2`, so BCa
reproduces the percentile interval; for a skewed statistic the interval becomes
asymmetric about `point` in the transformation-respecting direction.
Deterministic under a fixed `rng`.
"""
function bca_ci(data::AbstractVector, statistic;
                nboot::Integer=2000, alpha::Real=0.05,
                rng::AbstractRNG=default_rng())
    n = length(data)
    n >= 2 || throw(ArgumentError("BCa needs at least 2 observations"))
    nboot > 0 || throw(ArgumentError("nboot must be positive"))
    (0.0 < alpha < 1.0) || throw(ArgumentError("alpha must be in (0,1)"))

    point = Float64(statistic(data))

    # bootstrap distribution (same construction as bootstrap_ci)
    boots = Vector{Float64}(undef, nboot)
    sample = similar(data, n)
    @inbounds for b in 1:nboot
        for i in 1:n
            sample[i] = data[rand(rng, 1:n)]
        end
        boots[b] = Float64(statistic(sample))
    end

    z0 = _bca_z0(boots, point)
    a = _bca_acceleration(data, statistic)

    function adjusted_level(ptail)
        zt = _norm_quantile(ptail)
        s = z0 + zt
        denom = 1.0 - a * s
        # guard the (rare) singular/near-singular denominator
        adj = denom == 0.0 ? z0 + s * 1e12 : z0 + s / denom
        return clamp(_norm_cdf(adj), 0.0, 1.0)
    end

    a1 = adjusted_level(alpha / 2)
    a2 = adjusted_level(1.0 - alpha / 2)
    lo = _quantile7(boots, a1)
    hi = _quantile7(boots, a2)
    return (lo, point, hi)
end

# ---------------------------------------------------------------------------
# Empirical CDF
# ---------------------------------------------------------------------------

"""
    ecdf(x) -> Function

Return the empirical cumulative distribution function of `x` as a callable
`F(t)` giving the fraction of samples `≤ t`. Each evaluation is O(log n) via a
sorted copy captured in the closure.
"""
function ecdf(x::AbstractVector{<:Real})
    !isempty(x) || throw(ArgumentError("empty input"))
    s = sort(float.(x))
    n = length(s)
    return t -> searchsortedlast(s, float(t)) / n
end

# ---------------------------------------------------------------------------
# Gaussian KDE
# ---------------------------------------------------------------------------

# Silverman's rule-of-thumb bandwidth
function _silverman(x::AbstractVector{<:Real})
    n = length(x)
    σ = std(x)                                   # corrected
    iqr = _quantile7(x, 0.75) - _quantile7(x, 0.25)
    spread = iqr > 0 ? min(σ, iqr / 1.349) : σ
    spread > 0 || (spread = σ > 0 ? σ : 1.0)
    return 0.9 * spread * n^(-1 / 5)
end

"""
    kde(x; bandwidth=nothing, npoints=256) -> (grid, density)

Gaussian kernel density estimate. Evaluates
`f̂(g) = (1/nh) Σᵢ φ((g − xᵢ)/h)` on `npoints` evenly spaced points spanning
`[min(x) − 3h, max(x) + 3h]`, where `φ` is the standard normal pdf. When
`bandwidth` is `nothing` it defaults to Silverman's rule. The returned
`density` integrates to ≈ 1 over `grid` (exact up to boundary truncation).
"""
function kde(x::AbstractVector{<:Real}; bandwidth=nothing, npoints::Integer=256)
    !isempty(x) || throw(ArgumentError("empty input"))
    npoints >= 2 || throw(ArgumentError("npoints must be ≥ 2"))
    h = bandwidth === nothing ? _silverman(x) : float(bandwidth)
    h > 0 || throw(ArgumentError("bandwidth must be positive"))
    n = length(x)
    lo = minimum(x) - 3h
    hi = maximum(x) + 3h
    grid = collect(range(lo, hi; length=npoints))
    density = Vector{Float64}(undef, npoints)
    c = 1 / (n * h * sqrt(2π))
    @inbounds for j in 1:npoints
        acc = 0.0
        g = grid[j]
        for i in 1:n
            u = (g - x[i]) / h
            acc += exp(-0.5 * u * u)
        end
        density[j] = c * acc
    end
    return (grid, density)
end

end # module
