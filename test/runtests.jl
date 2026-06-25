using CellSizeControl
using Test

@testset "CellSizeControl" begin
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
end
