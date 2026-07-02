using CellSizeControl
using Test
using Aqua
using ExplicitImports
using JET
using Statistics: mean, std

@testset "CellSizeControl" begin
    # ---- Q: package-quality gates (release-readiness) ----
    @testset "Q — Aqua + ExplicitImports" begin
        Aqua.test_all(CellSizeControl; ambiguities=false)
        @test check_no_implicit_imports(CellSizeControl) === nothing
        @test check_all_explicit_imports_via_owners(CellSizeControl) === nothing
    end

    # ---- Q: JET static analysis (no inference errors / undefined-var / no-method) ----
    # The `rate` callback is an untyped keyword (`rate=qss_growth_rate`), so grow_to/
    # grow_for/cell_cycle/lineage_timecourse dispatch on it dynamically. That is a known,
    # accepted runtime-dispatch cost (a deliberate extension point, not a bug) and is NOT
    # an error JET's default analyzer reports; if a future JET flags it, scope it out with
    # target_modules rather than monomorphizing the API. This gate catches real inference
    # errors (typos, undefined vars, guaranteed no-method) only.
    @testset "Q — JET package analysis" begin
        JET.test_package(CellSizeControl; target_defined_modules=true)
    end

    # ---- L1: analytic limits — the slope discriminator recovers each regime ----
    @testset "L1 — slope discriminator (timer 2 / adder 1 / sizer 0)" begin
        s_timer = simulate_lineage(TimerRule(2.0); n=600, seed=1)
        s_adder = simulate_lineage(AdderRule(1.0); n=600, seed=2)
        s_sizer = simulate_lineage(SizerRule(2.0); n=600, seed=3)

        @test isapprox(size_control_slope(s_timer.Vb, s_timer.Vd), 2.0; atol=0.3)
        @test isapprox(size_control_slope(s_adder.Vb, s_adder.Vd), 1.0; atol=0.3)
        @test isapprox(size_control_slope(s_sizer.Vb, s_sizer.Vd), 0.0; atol=0.3)

        @test classify_control(size_control_slope(s_timer.Vb, s_timer.Vd)) === :timer
        @test classify_control(size_control_slope(s_adder.Vb, s_adder.Vd)) === :adder
        @test classify_control(size_control_slope(s_sizer.Vb, s_sizer.Vd)) === :sizer

        # inhibitor-dilution sizer: setpoint is exactly W/thresh, and it's a sizer
        ids = InhibitorDilutionSizer(60.0, 1.5)        # V* = 40
        @test setpoint_volume(ids) == 40.0
        @test division_volume(ids, 7.0) == 40.0        # birth-size-independent
        s_ids = simulate_lineage(ids; V0=5.0, n=600, seed=4)
        @test isapprox(size_control_slope(s_ids.Vb, s_ids.Vd), 0.0; atol=0.3)
    end

    # ---- CC-5: the mechanistic Whi5:SBF bistable switch reproduces V* = W/θ ----
    # The phenomenological inhibitor-dilution sizer's law emerges from a double-negative
    # feedback: the OFF/G1 state disappears at a saddle-node c*, and V* = W/c* is linear in W.
    @testset "CC-5 — mechanistic Whi5:SBF switch (V* = W/θ from a bistable mechanism)" begin
        sw = Whi5SBFSwitch(18.0)

        # (1) genuinely bistable: across the switching window the OFF and ON branches differ
        @test whi5_sbf_steady(sw, 5.0; from_high=false) < 0.05          # small cell → OFF
        @test whi5_sbf_steady(sw, 1.0e-4; from_high=false) > 0.9        # big cell → ON
        lo = whi5_sbf_steady(sw, 0.8; from_high=false)
        hi = whi5_sbf_steady(sw, 0.8; from_high=true)
        @test hi - lo > 0.5                                             # hysteresis (two stable states)

        # (2) the emergent set-point is EXACTLY linear in W (c* is W-independent → a true sizer)
        @test setpoint_volume(Whi5SBFSwitch(36.0)) / setpoint_volume(Whi5SBFSwitch(18.0)) ≈ 2 atol =
            1e-2
        # the emergent threshold θ = c* does not depend on W
        @test whi5_sbf_threshold(Whi5SBFSwitch(36.0)) ≈ whi5_sbf_threshold(Whi5SBFSwitch(18.0)) atol =
            1e-3

        # (3) it reproduces the inhibitor-dilution law: an InhibitorDilutionSizer with the
        #     emergent θ has the SAME set-point as the mechanistic switch
        θ = whi5_sbf_threshold(sw)
        @test setpoint_volume(InhibitorDilutionSizer(sw.W, θ)) ≈ setpoint_volume(sw) rtol = 1e-6

        # (4) and it behaves as a sizer in a lineage (Vd–Vb slope → 0)
        s_sw = simulate_lineage(sw; V0=5.0, n=400, seed=7)
        @test isapprox(size_control_slope(s_sw.Vb, s_sw.Vd), 0.0; atol=0.3)
    end

    # ---- the continuous linear size-control map Vd = αVb + β (Amir 2014) ----
    # one knob α sweeps sizer(0)→adder(1)→timer(2); the recovered slope must track α, and a
    # lineage with daughter fraction f is homeostatic iff α·f < 1.
    @testset "LinearSizeControl — slope tracks α; homeostasis iff α·f<1" begin
        @test division_volume(LinearSizeControl(0.7, 12.0), 20.0) == 0.7 * 20.0 + 12.0
        for α in (0.0, 0.5, 1.0, 1.5)
            s = simulate_lineage(LinearSizeControl(α, 20.0); V0=20.0, n=800, cv=0.06, seed=11)
            @test isapprox(size_control_slope(s.Vb, s.Vd), α; atol=0.2)
        end
        # homeostasis boundary at α = 1/f: α·f<1 stays bounded; α·f>1 runs away
        stable = simulate_lineage(LinearSizeControl(1.5, 20.0); daughter_fraction=0.5, n=200,
                                  cv=0.0, seed=12)                       # α·f = 0.75 < 1
        @test maximum(stable.Vb) < 1e3
        runaway = simulate_lineage(LinearSizeControl(2.5, 20.0); daughter_fraction=0.5, n=200,
                                   cv=0.0, seed=13)                      # α·f = 1.25 > 1
        @test last(runaway.Vb) > 10 * first(runaway.Vb)
    end

    # ---- L1: analytic fixed point of the linear size-control map ----
    # The deterministic (noiseless) lineage Vb_{n+1} = f·(α·Vb_n + β) is an affine contraction
    # when α·f < 1, with the CLOSED-FORM fixed point Vb* = f·β/(1−α·f) (division Vd* = α·Vb*+β)
    # reached at geometric rate (α·f)^n. This pins the actual homeostatic set-point the
    # sizer/adder/timer family converges to — sharper than "stays bounded" — and specialises to
    # the textbook sizer (α=0 ⇒ Vb*=f·V*) and adder (α=1 ⇒ Vb*=f·Δ/(1−f)) limits.
    @testset "L1 — linear-map fixed point Vb* = fβ/(1−αf) + scale invariance" begin
        for (α, β, f) in ((0.0, 40.0, 0.5), (1.0, 20.0, 0.5), (1.5, 20.0, 0.5),
                          (0.7, 12.0, 0.6), (1.2, 8.0, 0.4))
            s = simulate_lineage(LinearSizeControl(α, β); V0=5.0, n=400, cv=0.0,
                                 daughter_fraction=f)
            Vb_star = f * β / (1 - α * f)
            @test isapprox(last(s.Vb), Vb_star; rtol=1e-6)              # birth set-point
            @test isapprox(last(s.Vd), α * Vb_star + β; rtol=1e-6)      # division set-point
        end
        # sizer/adder specialisations of the closed form (the textbook limits)
        sz = simulate_lineage(SizerRule(60.0); V0=5.0, n=300, cv=0.0, daughter_fraction=0.5)
        @test isapprox(last(sz.Vb), 0.5 * 60.0; rtol=1e-6)             # α=0 → Vb* = f·V*
        ad = simulate_lineage(AdderRule(20.0); V0=5.0, n=400, cv=0.0, daughter_fraction=0.5)
        @test isapprox(last(ad.Vb), 0.5 * 20.0 / (1 - 0.5); rtol=1e-6) # α=1 → Vb* = fΔ/(1−f)

        # geometric convergence at the contraction rate α·f: |Vb_n − Vb*| ~ (αf)^n exactly
        α, β, f = 1.5, 20.0, 0.5
        s = simulate_lineage(LinearSizeControl(α, β); V0=5.0, n=60, cv=0.0, daughter_fraction=f)
        Vb_star = f * β / (1 - α * f)
        err = abs.(s.Vb .- Vb_star)
        @test isapprox(err[20] / err[10], (α * f)^10; rtol=1e-6)        # rate = (αf)^n

        # dimensional invariance: rescaling the volume UNIT (Vb, Vd → k·) leaves the
        # sizer/adder/timer slope unchanged — it is a dimensionless classifier
        base = simulate_lineage(AdderRule(10.0); n=300, seed=3)
        s0 = size_control_slope(base.Vb, base.Vd)
        for k in (1e-3, 1e3)
            @test isapprox(size_control_slope(k .* base.Vb, k .* base.Vd), s0; rtol=1e-9)
        end
    end

    # ---- L2: reference behaviour — the timer collapses, the sizer is stable ----
    # (reproduces the yeast-wcm CellSize finding: a sub-doubling timer drives
    #  daughters toward 0; the inhibitor-dilution sizer holds them at V*/2.)
    @testset "L2 — sizer stabilizes a lineage the timer collapses" begin
        timer = simulate_lineage(TimerRule(1.6); n=40, cv=0.0, seed=5)   # sub-doubling
        @test last(timer.Vb) < 0.2 * first(timer.Vb)                     # collapse → 0

        sizer = simulate_lineage(InhibitorDilutionSizer(60.0, 1.5); n=40, cv=0.0, seed=6)
        born = sizer.Vb[10:end]                                          # after transient
        @test all(b -> isapprox(b, 20.0; atol=1e-6), born)              # stable at V*/2 = 20

        # saturating open-loop timer (the course-model failure mode borged from CellSize.jl):
        # the bud = the growth increment g(A−V), which SHRINKS as the mother creeps toward A
        sat = SaturatingTimerRule(0.6, 69.0)
        @test division_volume(sat, 20.0) ≈ 20.0 + 0.6 * (69.0 - 20.0)
        buds = saturating_timer_buds(sat; V0=20.0, n=8)
        @test all(diff(buds) .< 0)                                       # daughters shrink (the bug)
        @test buds[end] / buds[1] < 0.1                                  # ~35%→5% drift the sizer fixes
    end

    # ---- maternal-age asymmetry erosion: old mothers → larger daughters ----
    @testset "maternal-age asymmetry erosion (CS-DA)" begin
        # the erosion function: alpha0 at age 0, -> alpha_max as age grows
        @test aging_daughter_fraction(0; alpha0=0.3, alpha_max=0.5, tau=10) == 0.3
        @test isapprox(
            aging_daughter_fraction(1e6; alpha0=0.3, alpha_max=0.5, tau=10), 0.5; atol=1e-6
        )
        # alpha0 == alpha_max reduces to a fixed fraction (no aging)
        @test aging_daughter_fraction(7; alpha0=0.4, alpha_max=0.4, tau=5) == 0.4

        # on a sizer, the daughter birth size RISES with maternal generation (Kennedy 1994)
        L = simulate_aging_lineage(
            SizerRule(60.0); n=25, alpha0=0.32, alpha_max=0.5, tau=8.0, cv=0.0
        )
        @test L.Vdaughter[end] > L.Vdaughter[1]                 # late daughters larger
        @test issorted(L.Vdaughter)                             # monotone rise (deterministic)
        @test isapprox(L.Vdaughter[1], 0.32 * 60.0; atol=1e-6)  # first daughter = alpha0·V*
        # the SAME erosion drives inherited DAMAGE: daughters of old mothers more damaged
        @test L.Ddaughter[end] > L.Ddaughter[1]                 # damage rises with maternal age
        @test issorted(L.Ddaughter)                             # monotone (size + fitness, one fn)

        # stochastic-AGE refinement: damage_cv adds noise to the per-cycle damage production,
        # so the inherited damage varies cell-to-cell (a real distribution); damage_cv=0 is
        # deterministic (seed-independent). Check at a fixed late generation across seeds.
        det1 =
            simulate_aging_lineage(SizerRule(60.0); n=20, damage_cv=0.0, seed=1).Ddaughter
        det2 =
            simulate_aging_lineage(SizerRule(60.0); n=20, damage_cv=0.0, seed=7).Ddaughter
        @test det1 == det2                                      # no damage noise → deterministic
        noisy = [
            simulate_aging_lineage(SizerRule(60.0); n=20, damage_cv=0.3, seed=s).Ddaughter[end]
            for s in 1:200
        ]
        @test std(noisy) > 0.05 * mean(noisy)                   # damage_cv>0 → real spread
        @test isapprox(mean(noisy), det1[end]; rtol=0.1)        # noise is mean-preserving-ish

        # no erosion (alpha0 == alpha_max) → constant daughter size
        flat = simulate_aging_lineage(
            SizerRule(60.0); n=15, alpha0=0.4, alpha_max=0.4, cv=0.0
        )
        @test all(d -> isapprox(d, 0.4 * 60.0; atol=1e-6), flat.Vdaughter)

        # maternal enlargement: the division set-point rises with age (old mothers larger),
        # so Vdivision grows from V* toward V*·(1+enlarge_max); default enlarge_max=0 is flat
        enl = simulate_aging_lineage(
            SizerRule(60.0);
            n=25,
            alpha0=0.32,
            alpha_max=0.5,
            tau=8.0,
            enlarge_max=0.45,
            enlarge_tau=8.0,
            cv=0.0,
        )
        @test issorted(enl.Vdivision)                          # mother at division enlarges
        @test isapprox(enl.Vdivision[1], 60.0; atol=1e-6)      # gen 1 at the base set-point
        @test enl.Vdivision[end] > 1.3 * enl.Vdivision[1]      # rises toward V*(1+0.45)
        @test all(isapprox.(flat.Vdivision, 60.0; atol=1e-6))  # enlarge_max=0 stays fixed

        # the MOTHER NEVER SHRINKS at division (monotonic body) — the corrected accounting:
        # she keeps her body, only the bud leaves, so the next cycle's birth = this division
        @test issorted(enl.Vbirth)                             # mother at Start monotonic non-decreasing
        @test enl.Vbirth[2:end] == enl.Vdivision[1:(end - 1)]  # division size carried over, no slice lost

        # Kennedy 1994 (magnitude): across the replicative lifespan (~25–30 divisions;
        # Schnitzer 2022), an old mother's daughters are substantially larger than a young
        # mother's — here the late daughter is ≳1.4× the first.
        ken = simulate_aging_lineage(
            SizerRule(60.0);
            n=30,
            alpha0=0.32,
            alpha_max=0.5,
            tau=8.0,
            enlarge_max=0.45,
            cv=0.0,
        )
        @test ken.Vdaughter[end] > 1.4 * ken.Vdaughter[1]
        @test length(ken.gen) == 30                            # spans the replicative lifespan

        # phantom founder: prepend gen 0 = the founder's OWN birth (born at V0), so the
        # lineage is uniform (every cell born a daughter) without shifting the real gens
        ph = simulate_aging_lineage(
            SizerRule(60.0);
            V0=12.0,
            n=25,
            alpha0=0.32,
            alpha_max=0.5,
            tau=8.0,
            phantom_founder=true,
            cv=0.0,
        )
        @test ph.gen[1] == 0                                    # gen-0 founder row prepended
        @test ph.phantom[1] && !any(ph.phantom[2:end])          # exactly one phantom, at front
        @test ph.Vdaughter[1] == 12.0                           # founder born at V0
        @test isnan(ph.Vbirth[1]) && isnan(ph.Vdivision[1])     # phantom mother undefined
        @test ph.Ddaughter[1] == 0.0                            # founder pristine (no inherited damage)
        @test length(ph.gen) == length(L.gen) + 1               # exactly one extra row vs default
        @test ph.Vdaughter[2:end] == L.Vdaughter                # real gens 1..N unchanged

        # Mechanism decomposition (ablation): the daughter-lifespan deficit is driven principally
        # by the maternal-damage trajectory, NOT by the age-erosion of the asymmetry. Inherited
        # damage Ddaughter rises with maternal age even when the
        # asymmetry fraction is held CONSTANT (alpha0 == alpha_max), because the mother's damage
        # pool D_m accumulates. This is the honest basis for the "r(a) couples, D_m(a) drives"
        # framing in the manuscript (gen_daughter_rls_ablation.jl quantifies the fold).
        const_share = simulate_aging_lineage(
            SizerRule(60.0); n=25, alpha0=0.4, alpha_max=0.4, tau=8.0, damage_form=1.0, cv=0.0,
        )
        @test const_share.Ddaughter[end] > const_share.Ddaughter[2]     # rises w/ age at fixed share
        # and the age-eroding share amplifies (does not generate) the late-age inherited burden:
        eroding_share = simulate_aging_lineage(
            SizerRule(60.0); n=25, alpha0=0.4, alpha_max=0.9, tau=8.0, damage_form=1.0, cv=0.0,
        )
        @test eroding_share.Ddaughter[end] > const_share.Ddaughter[end]  # erosion amplifies late burden
    end

    # ---- AGE-2: the replicative lifespan EMERGES from autocatalytic damage + a threshold ----
    @testset "emergent replicative lifespan (AGE-2)" begin
        # deterministic limit (no per-division noise, no cell-to-cell spread) → a fixed RLS
        r1 = replicative_lifespan(; cv=0.0, crit_cv=0.0, seed=1)
        r2 = replicative_lifespan(; cv=0.0, crit_cv=0.0, seed=99)
        @test r1 == r2 && r1 > 0                                 # seed-independent + finite

        # monotonic in the mechanism: higher viability threshold → longer life; faster
        # autocatalysis (kappa) → shorter life (the senescence accelerator)
        lo = replicative_lifespan(; D_crit=20.0, cv=0.0, crit_cv=0.0)
        hi = replicative_lifespan(; D_crit=40.0, cv=0.0, crit_cv=0.0)
        @test hi > lo
        slow = replicative_lifespan(; kappa=0.05, cv=0.0, crit_cv=0.0)
        fast = replicative_lifespan(; kappa=0.30, cv=0.0, crit_cv=0.0)
        @test slow > fast

        # a mother that segregates damage to her bud outlives one that keeps all of it
        @test replicative_lifespan(; segregate=true, cv=0.0, crit_cv=0.0) >
            replicative_lifespan(; segregate=false, cv=0.0, crit_cv=0.0)

        # Schnitzer 2022 calibration: the DEFAULT distribution has mean ≈ 25 and CV ≈ 0.3
        s = lifespan_distribution(3000)
        m = sum(s) / length(s)
        sd = sqrt(sum((x - m)^2 for x in s) / (length(s) - 1))
        @test 22.0 < m < 28.0                                   # mean RLS in the budding-yeast range
        @test 0.22 < sd / m < 0.40                              # realistic spread (cell-to-cell)
        @test all(>(0), s)                                      # every cell divides at least once
        @test length(lifespan_distribution(50)) == 50
    end

    # ---- AGE-3: daughter RLS vs maternal age — the CONVEX shape prediction ----
    # One age-eroding asymmetry ties the size face to the fitness face: a daughter born to an
    # age-`a` mother inherits a share φ(a) of the mother's accumulated damage D_m(a) and so
    # begins partway up her own autocatalytic damage trajectory. This is an out-of-sample
    # prediction (the McCormick mother-RLS calibration + the size-face asymmetry, nothing refit)
    # against Kennedy 1994. The population curve (averaged over lifespan + threshold
    # heterogeneity — what a microfluidic RLS study actually measures) is monotone-decreasing
    # AND convex: the decline is steepest at young/mid maternal age and flattens toward old age,
    # because daughter RLS is a concave (≈ logarithmic) function of inherited damage, so it is
    # most sensitive to D_m at low damage. This is a sharper, falsifiable shape claim than the
    # two-bucket Kennedy fold.
    @testset "daughter RLS vs maternal age — convex shape (AGE-3)" begin
        # inherited damage seed shortens life, deterministically (the D0 primitive)
        @test replicative_lifespan(; D0=20.0, cv=0.0, crit_cv=0.0) <
            replicative_lifespan(; cv=0.0, crit_cv=0.0)
        # damage_trajectory contract: one entry per division, monotone, born pristine
        tr = damage_trajectory(; cv=0.0, crit_cv=0.0)
        @test length(tr) == replicative_lifespan(; cv=0.0, crit_cv=0.0)
        @test tr[1] == 0.0 && issorted(tr)

        # passive volume-proportional inheritance share φ(a) = α(a)/α_max, α the size-face
        # asymmetry fixed by the daughter-size increase (Johnston 1966 / Yang 2011).
        φ(a) = aging_daughter_fraction(a; alpha0=0.69, alpha_max=0.90, tau=14.0) / 0.90

        # build the population fraction-binned daughter-RLS curve from the package model
        N, nbin, dseed = 4000, 10, 1_000_000
        binsum = zeros(nbin); binN = zeros(Int, nbin)
        for m in 1:N
            traj = damage_trajectory(; seed=m)              # a fresh mother: her D_m(a) series
            L = length(traj); L < 2 && continue
            for a in 0:(L - 1)
                Ld = replicative_lifespan(; D0=φ(a) * traj[a + 1], seed=(dseed += 1))
                b = clamp(floor(Int, (a / (L - 1)) * nbin) + 1, 1, nbin)
                binsum[b] += Ld; binN[b] += 1
            end
        end
        @test all(>(0), binN)
        x = [(b - 0.5) / nbin for b in 1:nbin]
        y = binsum ./ binN

        # (1) monotone decreasing in maternal age
        @test all(<(0), diff(y))
        # (2) substantial deficit young → old (young daughters near the full ≈25-div lifespan;
        #     old-mother daughters strongly shortened)
        @test y[1] > 20.0 && y[end] < 10.0
        # (3) CONVEX: the curve lies below the straight chord joining its endpoints, with a real
        #     margin (a linear decline would sit ON the chord). Decelerating decline.
        chord = y[1] .+ (y[end] - y[1]) .* (x .- x[1]) ./ (x[end] - x[1])
        @test all(y .<= chord .+ 1e-9)
        @test minimum(y .- chord) < -0.5           # observed margin ≈ -1.6; robust across seeds
        # (4) positive quadratic curvature (least-squares aᵢx²+bx+c ⇒ a>0 for convex)
        A = hcat(x .^ 2, x, ones(length(x)))
        a2 = (A \ y)[1]
        @test a2 > 1.0                             # observed ≈ 7

        # mechanism check: the convexity is an EMERGENT POPULATION effect. A single noiseless
        # mother's daughter curve is essentially LINEAR (a near-constant ≈ −1 div/gen decline),
        # NOT convex — it stays close to its own endpoint chord. Averaging over the lifespan +
        # threshold heterogeneity is what bends the population shape convex.
        det_traj = damage_trajectory(; cv=0.0, crit_cv=0.0)
        Ldet = length(det_traj)
        yd = Float64[replicative_lifespan(; D0=φ(a) * det_traj[a + 1], cv=0.0, crit_cv=0.0)
                     for a in 0:(Ldet - 1)]
        xd = collect(0:(Ldet - 1))
        chord_d = yd[1] .+ (yd[end] - yd[1]) .* (xd .- xd[1]) ./ (xd[end] - xd[1])
        @test all(<=(0), diff(yd))                 # single mother: monotone non-increasing
        @test minimum(yd .- chord_d) > -0.5        # …and ≈ linear (hugs its chord), unlike (3)
    end

    # ---- lineage timecourse: the mother is monotonic (never shrinks) over the lifespan ----
    @testset "lineage_timecourse — monotonic mother, buds detach" begin
        tc = lineage_timecourse(; n_max=29)
        @test issorted(tc.Vmother)                       # mother volume never decreases
        @test tc.Vmother[end] > tc.Vmother[1]            # mother enlarges with age
        @test minimum(tc.Vbud) ≈ 0.0                     # bud returns to 0 at each division
        @test maximum(tc.Vbud) > 0.0                     # buds grow during the budded phase
    end

    # ---- L2 reference: the two-step G1 reproduces Di Talia 2007 mother/daughter G1 ----
    @testset "L2 — Di Talia two-step G1 (mother ~19, daughter ~45 min)" begin
        Vstar, T_cln2, tau_bud = 36.0, 19.0, 52.0
        # a mother is born already ≥ V* → only the Cln2 timer runs → G1 ≈ T_cln2
        mom = cell_cycle(Vstar + 5; Vstar=Vstar, T_cln2=T_cln2, tau_bud=tau_bud)
        @test isapprox(mom.G1, T_cln2; atol=1e-6)
        # a daughter is born small → spends real time growing to V* → a longer G1 (~45 min)
        Vdau = cell_cycle(Vstar; Vstar=Vstar, T_cln2=T_cln2, tau_bud=tau_bud).Vdaughter
        dau = cell_cycle(Vdau; Vstar=Vstar, T_cln2=T_cln2, tau_bud=tau_bud)
        @test 40.0 < dau.G1 < 50.0                       # Di Talia daughter G1 ≈ 45.5
        @test dau.G1 > mom.G1                            # the asymmetry emerges, not imposed
        # exponential growth gives the same qualitative split (a known-good rate function)
        e = cell_cycle(20.0; Vstar=Vstar, rate=exponential_growth_rate(0.0077))
        @test e.G1 > T_cln2
    end

    # ---- L3: cross-source consistency — slope ordering is monotone in fold ----
    @testset "L3 — monotone slope vs control strength" begin
        slopes = [
            size_control_slope(
                simulate_lineage(TimerRule(f); n=500, seed=7).Vb,
                simulate_lineage(TimerRule(f); n=500, seed=7).Vd,
            ) for f in (1.6, 1.8, 2.0, 2.2)
        ]
        @test issorted(slopes)                                          # stronger timer → larger slope
        @test all(>(0), slopes)
    end

    # ---- perf guard: the compute-campaign hot paths stay O(n), not O(n²) ----
    # A regression tripwire ahead of the CC-1..CC-5 population campaigns (10^5–10^6 cells):
    # for each hot path, measure the bytes allocated at a base workload and at 10× it, and
    # assert the growth is near-linear — far below the ~100× a quadratic blow-up would show —
    # plus a generous absolute per-element ceiling. Uses Base `@allocated` (deterministic, no
    # extra dependency); the *scaling ratio* — not the raw byte count — is the portable guard,
    # so it holds across Julia versions where the absolute allocation sizes drift.
    @testset "perf — hot paths allocate O(n) (no quadratic blow-up)" begin
        alloc(f) = (f(); @allocated f())        # warm up, then bytes for one call

        # 10× the workload must not cost ~100× (quadratic); linear ⇒ ratio ≲ 10. A bound of
        # 25 leaves headroom for measurement noise / fixed overhead yet stays 4× below O(n²).
        # ceil_per is a loose per-element byte ceiling catching a constant-factor blow-up.
        cases = (
            ("simulate_lineage",
                n -> simulate_lineage(SizerRule(2.0); n=n, seed=1), 400, 4000, 1_000),
            ("lifespan_distribution",
                n -> lifespan_distribution(n; seed0=1), 50, 500, 100_000),
            ("simulate_aging_lineage",
                n -> simulate_aging_lineage(InhibitorDilutionSizer(1.0, 0.025); n=n, seed=1),
                25, 250, 8_000),
            ("lineage_timecourse",
                n -> lineage_timecourse(; n_max=n), 29, 290, 100_000),
        )
        for (name, f, base, big, ceil_per) in cases
            a_base = alloc(() -> f(base))
            a_big = alloc(() -> f(big))
            @test a_big / a_base < 25            # near-linear scaling, not O(n²)
            @test a_big / big < ceil_per         # loose absolute per-element ceiling
        end
    end
end
