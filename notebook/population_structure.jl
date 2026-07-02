### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #= mock-bind for running outside Pluto =#
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 5b1f0a20-3d11-4efc-bb61-000000000001
md"""
# Three scale-free laws of the growing culture

A budding-yeast culture in balanced exponential growth is not a bag of identical cells. Because
every mother buds one daughter per division and the *daughter fraction erodes with maternal age*,
the population carries a fixed demographic structure: the mothers of replicative age `a` are a
constant fraction `2^{-(a+1)}` of the culture (half virgin daughters, a quarter age-1 mothers, an
eighth age-2 …). Sampling the age-eroding division asymmetry through that geometric age law yields
three closed-form, **scale-free** predictions — signatures of the mechanism that survive without
knowing the set-point `V^\ast`:

1. **The newborn-size law** — the virgin-daughter size distribution is a right-skewed geometric
   *mixture*, not a symmetric bell.
2. **The extant-vs-newborn divergence** — the mean cell you pull from an unsorted culture is
   `D ≈ 1.97×` larger than the mean virgin daughter.
3. **The senescence age-law correction** — at short lifespan the geometric law flattens to a
   *truncated* geometric set by the discrete Euler–Lotka equation, `λ < 2`.

Everything below reuses the exported `newborn_size_law`, `extant_size_law`, and
`senescence_age_law`. Move the sliders: the set-point knob shows scale-freeness (the moments do
not move); the lifespan knob shows the geometric→truncated-geometric flattening.
"""

# ╔═╡ 5b1f0a20-3d11-4efc-bb61-000000000002
using CellSizeControl, PlutoUI, Plots, Statistics

# ╔═╡ 5b1f0a20-3d11-4efc-bb61-000000000003
md"""
### The calibrated demographic knobs

These are the model's fitted asymmetry and maternal-enlargement parameters (the same ones that
reproduce the Soifer–Amir slopes and the Johnston/Yang old-mother→large-daughter trend). The
age-eroding daughter fraction runs `α₀ = 0.32 → α_max = 0.5` with timescale `τ = 8` divisions, and
the set-point enlarges with maternal age by up to `45%` (`τ_e = 8`).
"""

# ╔═╡ 5b1f0a20-3d11-4efc-bb61-000000000004
begin
    A0, AMAX, TAU = 0.32, 0.5, 8.0
    EM, ET = 0.45, 8.0
    # Okabe–Ito (the repo's palette): blue = mother/extant, vermillion = daughter/newborn.
    BLUE, VERM, ORANGE, GREEN = "#0072B2", "#D55E00", "#E69F00", "#009E73"
end

# ╔═╡ 5b1f0a20-3d11-4efc-bb61-000000000005
md"""
## 1 · The newborn-size law — a right-skewed geometric mixture

Each age-`a` mother is a fraction `2^{-(a+1)}` of the culture and buds one daughter of size
`frac(a)·V^\ast·enlarge(a)`. So the newborn (virgin-daughter) sizes are the discrete mixture
`{(2^{-(a+1)}, frac(a)·V^\ast·enlarge(a))}` — the comb below. Old, rare mothers make large, rare
daughters: the tail is to the **right**.

**Set-point** `V^\ast` (fL): $(@bind Vstar Slider(20:5:240; default = 60, show_value = true))
"""

# ╔═╡ 5b1f0a20-3d11-4efc-bb61-000000000006
nb = newborn_size_law(; alpha0 = A0, alpha_max = AMAX, tau = TAU,
                      enlarge_max = EM, enlarge_tau = ET, Vstar = Float64(Vstar))

# ╔═╡ 5b1f0a20-3d11-4efc-bb61-000000000007
let
    ages = 0:20
    ws = [2.0^(-(a + 1)) for a in ages]
    ss = [aging_daughter_fraction(a; alpha0 = A0, alpha_max = AMAX, tau = TAU) *
          Float64(Vstar) * (1 + EM * (1 - exp(-a / ET))) for a in ages]
    plot(ss, ws; seriestype = :sticks, lw = 2, marker = :circle, ms = 4,
         color = VERM, legend = false,
         xlabel = "newborn size  (fL)", ylabel = "population weight  2^-(a+1)",
         title = "newborn-size comb at V* = $(Vstar) fL   (mean = $(round(nb.mean, digits = 2)) fL)")
    vline!([nb.mean]; color = BLUE, ls = :dash, lw = 2)
end

# ╔═╡ 5b1f0a20-3d11-4efc-bb61-000000000008
md"""
At this set-point the mean newborn size is **$(round(nb.mean, digits = 2)) fL**, a fixed multiple
**R = $(round(nb.ratio, digits = 4))** of the youngest-mother daughter `α₀·V^\ast`, with
**CV = $(round(nb.cv, digits = 4))** and **skewness = $(round(nb.skew, digits = 3)) > 0**
(right-skewed). Drag `V^\ast` above: the *mean* tracks the set-point, but **R, CV, and skew do not
move** — they are scale-free. The distribution collapses to a symmetric point mass (CV = skew = 0)
only if the daughter size stops depending on maternal age.
"""

# ╔═╡ 5b1f0a20-3d11-4efc-bb61-000000000009
let
    Vs = [30.0, 60.0, 120.0, 240.0]
    laws = [newborn_size_law(; alpha0 = A0, alpha_max = AMAX, tau = TAU,
                             enlarge_max = EM, enlarge_tau = ET, Vstar = V) for V in Vs]
    p1 = plot(Vs, [l.ratio for l in laws]; marker = :circle, ms = 6, lw = 2, color = ORANGE,
              legend = false, xlabel = "set-point V*  (fL)", ylabel = "ratio R",
              title = "R is flat in V* (scale-free)", ylim = (0, 2))
    p2 = plot(Vs, [l.cv for l in laws]; marker = :square, ms = 6, lw = 2, color = GREEN,
              legend = false, xlabel = "set-point V*  (fL)", ylabel = "CV",
              title = "CV is flat in V*", ylim = (0, 0.3))
    p3 = plot(Vs, [l.skew for l in laws]; marker = :diamond, ms = 6, lw = 2, color = VERM,
              legend = false, xlabel = "set-point V*  (fL)", ylabel = "skewness",
              title = "skew is flat in V*", ylim = (0, 2))
    plot(p1, p2, p3; layout = (1, 3), size = (900, 280))
end

# ╔═╡ 5b1f0a20-3d11-4efc-bb61-00000000000a
md"""
## 2 · Extant vs newborn — the standing culture is bigger than its daughters

Snapshot the culture and you sample every cell at its most-recent division. The age-0 cells are the
small buds (the newborn law), but a cell of age `a ≥ 1` is a **mother carrying her full retained
body** `V^\ast·enlarge(a−1)` — the mother keeps her body; only the bud leaves. Weighting each age
class by `2^{-(a+1)}` makes the mean *extant* cell far larger than the mean *newborn*.
"""

# ╔═╡ 5b1f0a20-3d11-4efc-bb61-00000000000b
ex = extant_size_law(; alpha0 = A0, alpha_max = AMAX, tau = TAU,
                     enlarge_max = EM, enlarge_tau = ET, Vstar = Float64(Vstar))

# ╔═╡ 5b1f0a20-3d11-4efc-bb61-00000000000c
let
    bar(["mean newborn", "mean extant"], [ex.newborn_mean, ex.extant_mean];
        color = [VERM, BLUE], legend = false, ylabel = "mean size  (fL)",
        title = "divergence D = extant / newborn = $(round(ex.divergence, digits = 3))   (V* = $(Vstar) fL)")
end

# ╔═╡ 5b1f0a20-3d11-4efc-bb61-00000000000d
md"""
The mean extant cell is **D = $(round(ex.divergence, digits = 3))×** the mean newborn
($(round(ex.extant_mean, digits = 1)) fL vs $(round(ex.newborn_mean, digits = 1)) fL). Slide
`V^\ast`: both means scale, but **D never moves** — another scale-free signature. In the clean
no-erosion, no-enlargement limit it collapses to the exact closed form

`D = (1 + α₀) / (2 α₀) = $(round((1 + A0) / (2 * A0), digits = 4))`   (at α₀ = $(A0)),

which the standing population's over-representation of large old mothers lifts to ≈1.97 once the
age-eroding asymmetry is switched on.
"""

# ╔═╡ 5b1f0a20-3d11-4efc-bb61-00000000000e
md"""
## 3 · The senescence correction — Euler–Lotka flattens the geometric law

The clean `2^{-(a+1)}` law assumes cells divide forever. Real mothers senesce: with a finite
replicative lifespan `rls` the dividing population grows by `λ < 2` per generation, set by the
discrete **Euler–Lotka** equation `λ = Σ_{a=0}^{rls−1} λ^{-a}`, and its age structure is the
**truncated geometric** `P(age=a) = λ^{-(a+1)}` for `a = 0 … rls−1`. As `rls → ∞`, `λ → 2` and the
naive law returns; at short `rls` the base `1/λ > 1/2` **flattens** the histogram.

**Replicative lifespan** `rls`: $(@bind rls Slider(2:1:30; default = 4, show_value = true)) divisions
"""

# ╔═╡ 5b1f0a20-3d11-4efc-bb61-00000000000f
sl = senescence_age_law(rls)

# ╔═╡ 5b1f0a20-3d11-4efc-bb61-000000000010
let
    ages = sl.ages
    naive = [2.0^(-(a + 1)) for a in ages]
    plot(ages, naive; seriestype = :sticks, lw = 2, marker = :circle, ms = 4, color = BLUE,
         label = "naive geometric 2^-(a+1)", xlabel = "replicative age a",
         ylabel = "fraction of dividing cells", legend = :topright,
         title = "rls = $(rls):  λ = $(round(sl.lambda, digits = 4)),  virgin fraction = $(round(1 / sl.lambda, digits = 4))")
    plot!(ages .+ 0.15, sl.p; seriestype = :sticks, lw = 2, marker = :diamond, ms = 4,
          color = VERM, label = "Euler–Lotka λ^-(a+1)")
end

# ╔═╡ 5b1f0a20-3d11-4efc-bb61-000000000011
md"""
Here `λ = $(round(sl.lambda, digits = 4))` and the virgin fraction is exactly `1/λ =
$(round(1 / sl.lambda, digits = 4))` — above `1/2`, and it truncates at age `rls−1`. Two anchors:
the one-division-shy limit `rls = 2` gives the **golden ratio** `λ = φ = $(round(senescence_age_law(2).lambda, digits = 5))`,
and a long lifespan recovers the clean law
(`λ($(100)) = $(round(senescence_age_law(100).lambda, digits = 5)) ≈ 2`). Drag `rls` up and watch
the vermillion truncated-geometric relax onto the blue `2^{-(a+1)}` comb.
"""

# ╔═╡ 5b1f0a20-3d11-4efc-bb61-000000000012
md"""
## Headline numbers & how to falsify them

At the calibrated parameters (`α₀ = 0.32 → α_max = 0.5`, `τ = 8`; enlarge `0.45`, `τ_e = 8`):

| Law | Headline | Scale-free invariant |
|---|---|---|
| Newborn size | mean `R ≈ 1.114·α₀V^\ast`, **CV ≈ 0.135**, **skew ≈ 1.562** (right-skewed) | R, CV, skew are independent of `V^\ast` |
| Extant vs newborn | **D ≈ 1.97** (no-erosion closed form `(1+α₀)/(2α₀) = 2.0625`) | D is independent of `V^\ast` |
| Senescence | `λ < 2`, virgin fraction `1/λ`; `λ(2) = φ`, `λ(4) = 1.9276`, `λ→2` as `rls→∞` | the age law is `λ^{-(a+1)}` |

**Falsification conditions.**
- A measured **newborn-size distribution that is symmetric or left-skewed** (skew ≤ 0), or whose
  CV/skew **fail to collapse onto the same numbers across strains of different mean size** (different
  `V^\ast`), refutes the geometric-age-sampling account. Genetically flattening the
  maternal-age→daughter-size relation should drive the newborn skew toward zero; if it does not, the
  single-asymmetry mechanism is wrong.
- Measuring the mean size of an **unsorted culture vs isolated virgin daughters** (Coulter/imaging)
  and finding a ratio **far from ≈1.97** — or one that is not the same number across set-points —
  refutes the balanced-growth, mother-keeps-body account.
- In a **short-lived (low-RLS) strain**, a dividing-cell age structure that stays `2^{-(a+1)}`
  rather than **flattening toward `λ^{-(a+1)}`** (`λ < 2`) refutes the Euler–Lotka correction.

Run it live: `julia --project=.` in the package, `using Pluto; Pluto.run()`, open this file.
"""

# ╔═╡ Cell order:
# ╟─5b1f0a20-3d11-4efc-bb61-000000000001
# ╠═5b1f0a20-3d11-4efc-bb61-000000000002
# ╟─5b1f0a20-3d11-4efc-bb61-000000000003
# ╠═5b1f0a20-3d11-4efc-bb61-000000000004
# ╟─5b1f0a20-3d11-4efc-bb61-000000000005
# ╠═5b1f0a20-3d11-4efc-bb61-000000000006
# ╠═5b1f0a20-3d11-4efc-bb61-000000000007
# ╟─5b1f0a20-3d11-4efc-bb61-000000000008
# ╠═5b1f0a20-3d11-4efc-bb61-000000000009
# ╟─5b1f0a20-3d11-4efc-bb61-00000000000a
# ╠═5b1f0a20-3d11-4efc-bb61-00000000000b
# ╠═5b1f0a20-3d11-4efc-bb61-00000000000c
# ╟─5b1f0a20-3d11-4efc-bb61-00000000000d
# ╟─5b1f0a20-3d11-4efc-bb61-00000000000e
# ╠═5b1f0a20-3d11-4efc-bb61-00000000000f
# ╠═5b1f0a20-3d11-4efc-bb61-000000000010
# ╟─5b1f0a20-3d11-4efc-bb61-000000000011
# ╟─5b1f0a20-3d11-4efc-bb61-000000000012
