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

# ╔═╡ e4e5ed1e-2d11-4efc-bb61-5906d0cc28d2
md"""
# Watch the cell-size sizer at work

A cell can control its size three ways: a **timer** (divide after a fixed time), an **adder** (add a fixed volume), or a **sizer** (divide at a fixed target size). Budding yeast uses an **inhibitor-dilution sizer** — Whi5 is diluted as the cell grows and Start fires at a threshold size.

The tell-tale is the **slope of division size on birth size**: ~0 for a sizer, ~1 for an adder, ~2 for a timer. Pick a rule, perturb the starting size, and watch how fast (or whether) the lineage forgets it.
"""

# ╔═╡ 66f09ecb-3970-4c81-9706-55cd537f5f1b
using CellSizeControl, PlutoUI, Plots, Statistics

# ╔═╡ 312c301c-0ad9-47b6-9b4a-810b206f21fc
md"""
control rule: $(@bind ruletype Select(["sizer (inhibitor dilution)", "adder", "timer"]))

start size V0: $(@bind V0 Slider(20:5:120; default = 80, show_value = true)) fL    noise CV: $(@bind cv Slider(0:0.01:0.2; default = 0.05, show_value = true))    generations: $(@bind n Slider(5:1:30; default = 15, show_value = true))
"""

# ╔═╡ df4303d2-27fd-4440-80aa-67d1b3198ae5
begin
	rule = ruletype == "adder" ? AdderRule(25.0) :
	       ruletype == "timer" ? TimerRule(2.0)   :
	                              SizerRule(50.0)
	o = simulate_lineage(rule; V0 = Float64(V0), n = n, cv = cv, daughter_fraction = 0.5, seed = 1)
end

# ╔═╡ 560803a3-e21c-4132-8261-67c655e2b356
plot(1:length(o.Vb), o.Vb; lw = 3, marker = :circle, legend = false,
	 xlabel = "generation", ylabel = "birth volume V_b  (fL)",
	 title = "$(ruletype): birth size across the lineage (started at V0 = $(V0) fL)")

# ╔═╡ 9c2b19b3-7d7d-485f-b212-353861aac28b
let
	b, d  = o.Vb, o.Vd
	slope = cov(b, d) / var(b)
	scatter(b, d; legend = false, xlabel = "V_b  (birth, fL)", ylabel = "V_d  (division, fL)",
	        title = "discriminator slope = $(round(slope, digits = 2))   (sizer~0, adder~1, timer~2)")
	xs = range(minimum(b), maximum(b); length = 2)
	plot!(xs, mean(d) .+ slope .* (xs .- mean(b)); lw = 2)
end

# ╔═╡ b0b210f2-2951-4908-8a0b-e8af6cff53ac
md"""
A **sizer** drives the slope toward **0**: division size is pinned near the target $V^\ast$ regardless of birth size, so a cell born too big or too small is corrected within a generation. An **adder** gives slope **1**, a **timer** slope **2** (size errors amplify). The inhibitor-dilution sizer is what makes budding yeast's size so tightly controlled.
"""

# ╔═╡ Cell order:
# ╟─e4e5ed1e-2d11-4efc-bb61-5906d0cc28d2
# ╠═66f09ecb-3970-4c81-9706-55cd537f5f1b
# ╟─312c301c-0ad9-47b6-9b4a-810b206f21fc
# ╠═df4303d2-27fd-4440-80aa-67d1b3198ae5
# ╠═560803a3-e21c-4132-8261-67c655e2b356
# ╠═9c2b19b3-7d7d-485f-b212-353861aac28b
# ╟─b0b210f2-2951-4908-8a0b-e8af6cff53ac
