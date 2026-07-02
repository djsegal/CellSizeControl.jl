# Coupling the size-CV amplification to the damage/RLS recursion (size-noise → aging channel).
#
# In the base model size and damage are INDEPENDENT: the replicative-lifespan recursion
#   D_{a+1} = D_a + production·(1 + κ·D_a)·(1 + cv·ξ),   senesce when D crosses D_crit,
# never references the size-control mode (with the default non-conserved damage, segregate=false,
# the per-division increment does not depend on α or f). The only free stochastic input is the
# per-division production noise `cv`.
#
# Minimal, parameter-free coupling — SIZE-DEPENDENT DAMAGE PRODUCTION. Per-cycle damage production
# scales with cell size (metabolic / biosynthetic / proteostatic load ∝ volume), so a fractional
# birth-size fluctuation δV/V passes 1:1 into a fractional production fluctuation δP/P. The
# stationary birth-size CV is mode-dependent (the CV-amplification frontier):
#   CV(Vb) = cv_size · A(α,f),   A(α,f) = 1 / sqrt(1 − (α·f)^2),
# so the damage-production noise the RLS recursion sees is cv_damage = A(α,f)·cv_size. A timer
# (α=2) amplifies size noise (A>1); a sizer (α=0) does not (A=1). This is the whole coupling: no
# new biology beyond "production ∝ volume," and mode enters RLS ONLY through cv_damage.
#
# Prediction. Mean RLS is set by the autocatalytic threshold crossing and is essentially invariant
# to the (mean-1) production noise; the RLS DISTRIBUTION broadens with cv_damage. So at matched
# mean division asymmetry f a timer has the same mean RLS as a sizer but a broader one, and the
# broadening GROWS as replicative aging symmetrizes division (f: 0.32→0.5, A_timer→∞). Isolating
# the channel (crit_cv=0, no cell-to-cell threshold spread) exposes it cleanly; at the realistic
# threshold heterogeneity (crit_cv=0.45) it is a subdominant correction swamped by the threshold
# spread — the falsification handle.
#
# Deterministic cross-check: at cv=0 the RLS is a single fixed number, IDENTICAL across modes —
# the mode effect lives entirely in the noise-driven spread. Monte-Carlo cross-check: the ensemble
# RLS mean/CV at fixed seeds.
#
# Emits size_noise_rls.csv (mode comparison, isolated + realistic threshold), size_noise_rls_aging.csv
# (RLS broadening vs the aging asymmetry axis), size_noise_rls_hist.csv (the RLS pmfs for the figure).
# Run: julia --project=. analysis/gen_size_noise_aging.jl
using CellSizeControl
using Statistics: mean, std
using Printf

here = @__DIR__

A(α, f) = 1 / sqrt(1 - (α * f)^2)          # size-CV amplification (return slope r = α·f)

const CV_SIZE = 0.06                        # per-division size noise (matches the CV-amplification run)
const N = 40_000                            # Monte-Carlo ensemble
const SEED0 = 1
const F0 = 0.40                             # matched mean division asymmetry (young 0.32 → aged 0.5)
const MODES = (("sizer", 0.0), ("adder", 1.0), ("timer", 2.0))

# Monte-Carlo RLS ensemble at a given production-noise cv (fixed seed bank).
function rls_samples(cv; crit_cv=0.0, n=N, seed0=SEED0, kwargs...)
    return [replicative_lifespan(; cv=cv, crit_cv=crit_cv, seed=seed0 + i, kwargs...) for i in 1:n]
end
rls_stats(cv; kwargs...) = (s = rls_samples(cv; kwargs...); (mean=mean(s), cv=std(s) / mean(s)))

# Deterministic cross-check: cv=0 → a single fixed RLS, mode-independent.
rls_det = replicative_lifespan(; cv=0.0, crit_cv=0.0)

# ---- (1) mode comparison at f = F0, isolated (crit_cv=0) and realistic (crit_cv=0.45) ----
cv_sizer = A(0.0, F0) * CV_SIZE             # sizer baseline damage noise (= cv_size)
open(joinpath(here, "size_noise_rls.csv"), "w") do io
    println(io, "mode,alpha,f,A,cv_size,cv_damage,crit_cv,rls_mean,rls_cv,rls_cv_ratio_vs_sizer")
    for crit_cv in (0.0, 0.45)
        base = rls_stats(cv_sizer; crit_cv=crit_cv).cv
        for (name, α) in MODES
            amp = A(α, F0)
            st = rls_stats(amp * CV_SIZE; crit_cv=crit_cv)
            @printf(io, "%s,%.1f,%.3f,%.5f,%.4f,%.5f,%.2f,%.4f,%.6f,%.5f\n",
                name, α, F0, amp, CV_SIZE, amp * CV_SIZE, crit_cv, st.mean, st.cv, st.cv / base)
        end
    end
end

# ---- (2) aging axis: RLS broadening vs division asymmetry f (isolated channel) ----
# As replicative aging symmetrizes division (f rises 0.32→0.5) the timer's size-CV amplifies
# (A→∞) so the sizer/timer RLS-CV gap widens; the sizer stays flat.
open(joinpath(here, "size_noise_rls_aging.csv"), "w") do io
    println(io, "f,alpha,mode,A,cv_damage,rls_mean,rls_cv,rls_cv_ratio_vs_sizer")
    for f in 0.32:0.02:0.48
        base = rls_stats(A(0.0, f) * CV_SIZE).cv
        for (name, α) in MODES
            amp = A(α, f)
            st = rls_stats(amp * CV_SIZE)
            @printf(io, "%.2f,%.1f,%s,%.5f,%.5f,%.4f,%.6f,%.5f\n",
                f, α, name, amp, amp * CV_SIZE, st.mean, st.cv, st.cv / base)
        end
    end
end

# ---- (3) RLS pmfs for the figure (sizer vs timer at f = F0, isolated channel) ----
s_sizer = rls_samples(A(0.0, F0) * CV_SIZE)
s_timer = rls_samples(A(2.0, F0) * CV_SIZE)
lo, hi = extrema(vcat(s_sizer, s_timer))
open(joinpath(here, "size_noise_rls_hist.csv"), "w") do io
    println(io, "rls,sizer_frac,timer_frac")
    for r in lo:hi
        @printf(io, "%d,%.6f,%.6f\n",
            r, count(==(r), s_sizer) / length(s_sizer), count(==(r), s_timer) / length(s_timer))
    end
end

# ---- headline receipt ----
st_s = rls_stats(A(0.0, F0) * CV_SIZE)
st_t = rls_stats(A(2.0, F0) * CV_SIZE)
@printf("deterministic RLS (cv=0) = %d divisions, mode-independent\n", rls_det)
@printf("f=%.2f  A_timer=%.4f (=5/3=%.4f)\n", F0, A(2.0, F0), 5 / 3)
@printf("isolated channel: RLS_CV sizer=%.4f timer=%.4f  ratio=%.3f  Δmean=%.3f div\n",
    st_s.cv, st_t.cv, st_t.cv / st_s.cv, st_t.mean - st_s.mean)
@printf("wrote size_noise_rls.csv + size_noise_rls_aging.csv + size_noise_rls_hist.csv\n")
