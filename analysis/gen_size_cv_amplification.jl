# Size-control variability amplification and the aging homeostasis boundary.
#
# In the linear size-control map Vd = αVb + β, a lineage with division asymmetry f (daughter
# fraction) has return slope r = α·f and is size-homeostatic iff r < 1 (Amir 2014). With
# per-division multiplicative noise cv, the steady-state birth-size variability obeys
#
#     CV(Vb) = cv / sqrt(1 - (α·f)^2),
#
# so cell-size CV amplifies toward the boundary r → 1 and diverges there. Because replicative
# aging erodes division asymmetry from f0 ≈ 0.32 toward symmetric f = 0.5
# (aging_daughter_fraction), the return slope r(a) = α·f(a) RISES with the mother's replicative
# age at fixed control mode. The critical control slope α_c(f) = 1/f equals exactly 2 — the
# timer slope — at the aging endpoint f = 0.5. Prediction: sizer (α=0) and adder (α=1) lineages
# stay homeostatic at every replicative age, but a timer-controlled lineage (α → 2) is driven to
# marginal loss of size homeostasis precisely as division symmetrizes with age; its birth-size CV
# grows without bound while the sizer's stays flat.
#
# Emits: size_cv_amplification.csv (measured vs predicted CV across the sizer→timer axis at
# three fixed asymmetries) and size_cv_aging.csv (the CV-amplification-vs-replicative-age curve
# for the three canonical modes). Run: julia --project=. analysis/gen_size_cv_amplification.jl
using CellSizeControl
using Statistics: mean, std
using Printf

here = @__DIR__

const CV = 0.06      # per-division multiplicative noise
const BETA = 20.0    # set-point scale (β)

# Steady-state CV of birth volume, estimated across the final birth of R independent lineages
# (each fully burned in after n divisions: |Vb_n − Vb*| ~ (α·f)^n → 0).
function cv_birth(α, f; cv=CV, β=BETA, n=200, R=4000, seed0=1)
    finals = [
        last(
            simulate_lineage(
                LinearSizeControl(α, β); V0=β, n=n, cv=cv, daughter_fraction=f, seed=seed0 + r
            ).Vb,
        ) for r in 1:R
    ]
    return std(finals) / mean(finals)
end

cv_pred(r; cv=CV) = cv / sqrt(1 - r^2)   # CV(Vb) = cv / sqrt(1 - (α·f)^2)

# ---- (a) CV amplification across the sizer→adder→timer axis, at three asymmetries ----
alphas = collect(0.0:0.1:1.9)
fracs = (0.32, 0.40, 0.50)          # young, mid, aged (symmetric) division
open(joinpath(here, "size_cv_amplification.csv"), "w") do io
    println(io, "f,alpha,r,cv_measured,cv_pred")
    for f in fracs, α in alphas
        r = α * f
        r >= 0.995 && continue      # skip the divergence itself (r → 1)
        cm = cv_birth(α, f)
        @printf(io, "%.3f,%.3f,%.4f,%.6f,%.6f\n", f, α, r, cm, cv_pred(r))
    end
end

# ---- (b) CV amplification vs replicative age for the three canonical control modes ----
# At maternal age a the division asymmetry is f(a) = aging_daughter_fraction(a) (0.32 → 0.5), so
# the return slope a cell experiences is r(a) = α·f(a). Amplification is reported relative to the
# sizer baseline (CV = cv); the sizer is flat (r ≡ 0), the timer runs to the boundary r → 1.
modes = (("sizer", 0.0), ("adder", 1.0), ("timer", 2.0))
ages = 0:2:60
base = cv_birth(0.0, 0.5)           # sizer baseline CV (≈ cv), age-independent
open(joinpath(here, "size_cv_aging.csv"), "w") do io
    println(io, "mode,alpha,age,f,r,amp_measured,amp_pred,homeostatic")
    for (name, α) in modes, a in ages
        f = aging_daughter_fraction(a)
        r = α * f
        if r >= 0.995
            @printf(io, "%s,%.1f,%d,%.4f,%.4f,%s,%s,%d\n", name, α, a, f, r, "NaN", "Inf", 0)
        else
            amp_m = cv_birth(α, f; R=2500) / base
            @printf(
                io, "%s,%.1f,%d,%.4f,%.4f,%.4f,%.4f,%d\n",
                name, α, a, f, r, amp_m, cv_pred(r) / CV, 1
            )
        end
    end
end

# headline receipt: the amplification factor at (α=1.6, f=0.5) → r=0.8 → theory 5/3
A = cv_birth(1.6, 0.5) / cv_birth(0.0, 0.5)
@printf(
    "wrote size_cv_amplification.csv + size_cv_aging.csv | A(α=1.6,f=0.5)=%.4f vs 5/3=%.4f | α_c(f=0.5)=%.1f (timer)\n",
    A, 5 / 3, 1 / 0.5
)
