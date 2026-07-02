# The return-map pole r = α·f — one eigenvalue, three signatures

The birth-size return map of a size-control lineage,

    Vb(n+1) = f·(α·Vb(n) + β)·(1 + cv·ξ),

linearises about its fixed point `Vb* = fβ/(1−αf)` to an **AR(1) process with a single
pole** `r = α·f` (`α = map_slope(rule)` the control slope, `f` the daughter fraction).
That one eigenvalue simultaneously sets the **stationary variance**, the **temporal
memory**, and the **transient relaxation** of the lineage — so three observables that look
independent are all reports of the same number. The lineage is size-homeostatic iff
`r < 1` (Amir 2014); at `r → 1` every signature diverges together.

This note unifies three results in the package around that pole. Each has an equation, a
headline, a gating test in `test/runtests.jl`, and a falsification.

## The three signatures

### 1. CV amplification — stationary variance
`CV(Vb) = cv / √(1 − r²)` (`analysis/gen_size_cv_amplification.jl`). Per-division
multiplicative noise `cv` is amplified toward the homeostasis boundary; the gain is
`1/√(1−r²)` and diverges at `r → 1`.

- **Headline:** at `(α=1.6, f=0.5) ⇒ r=0.8` the amplification is exactly `5/3`; the sizer
  (`r=0`) stays flat at `CV = cv`.
- **Gate:** `prediction — CV(Vb)=cv/√(1−(αf)²) + timer-critical aging boundary`.
- **Falsify:** birth-size CV flat across the sizer→timer axis, or independent of `r` as
  asymmetry erodes.

### 2. Size-noise → RLS broadening — a downstream coupling
Under size-dependent damage production (production ∝ volume) the RLS recursion sees noise
`cv_damage = A(α,f)·cv_size`, `A = 1/√(1−r²)` (`analysis/gen_size_noise_aging.jl`). Mean
RLS is invariant to this mean-1 noise; the RLS **distribution** broadens with `cv_damage`.

- **Headline:** at `f=0.4` the timer's RLS-CV is `≈1.42×` the sizer's; the ratio grows past
  `2×` as division symmetrizes (`f→0.5`, `A_timer→∞`). Honest detectability bound: at
  realistic threshold heterogeneity (`crit_cv=0.45`) the channel is swamped (ratio ≈ 1).
- **Gate:** `AGE-4 — size-noise → RLS broadening (timer vs sizer)`.
- **Falsify:** RLS-CV independent of control mode at matched `f` with the threshold spread
  controlled down.

### 3. Birth-size memory — temporal correlation + transient relaxation
The AR(1) pole sets both memory faces at once (`analysis/gen_size_memory.jl`):
lag-`k` lineage autocorrelation `ρ_k = r^k`, and nutrient-shift mean relaxation
`⟨Vb(n)⟩ = Vb*₂ + (Vb*₁−Vb*₂)·rⁿ` with memory time `τ = −1/ln r` generations.

- **Headline:** a sizer (`r=0`) is memoryless — every birth independent, a set-point shift
  absorbed in one division; an adder carries `≈1.44` generations; a timer carries the
  longest and diverges as `r=2f→1`. Two observables force the mode- and set-point-free
  **invariant** `CV(Vb)²·(1−ρ₁²) = cv²`, recovering the intrinsic per-division noise.
- **Gate:** `prediction — birth-size memory ρ_k=r^k + CV²(1−ρ1²)=cv² + nutrient relaxation`.
- **Falsify:** nonzero mother→daughter birth-size correlation in a sizer, `ρ₁ ≠ αf`, a
  mode-independent relaxation rate, or `CV²(1−ρ₁²)` failing to collapse across modes/set-points.

## The unifying table — {mode → r → CV-gain → ρ₁ → memory τ}

One pole `r = α·f`; everything below is a function of it (`f = 0.5`). `size_memory(rule; f)`
returns `(r, cv_gain, autocorr=ρ₁, memory_gen=τ)`.

| Mode  | α | r = α·f | CV-gain 1/√(1−r²) | ρ₁ = r | memory τ = −1/ln r (gen) |
|-------|---|---------|-------------------|--------|--------------------------|
| Sizer | 0 | 0.0     | 1.000             | 0.0    | 0 (memoryless)           |
| Adder | 1 | 0.5     | 1.155             | 0.5    | 1.443                    |
| Timer | 2 | 1.0     | ∞ (marginal)      | 1.0    | ∞                        |

At `f = 0.5` the timer sits exactly on the boundary `r = 1`. A young timer (`f = 0.4`,
`r = 0.8`) is still homeostatic: CV-gain `5/3`, `ρ₁ = 0.8`, `τ = 4.48` gen.

## The aging axis — r → 1 for the timer

Replicative aging erodes division asymmetry from `f₀ ≈ 0.32` toward symmetric `f = 0.5`
(`aging_daughter_fraction`), so the pole `r(a) = α·f(a)` **rises with maternal age** at
fixed mode. The critical control slope `α_c(f) = 1/f` equals exactly `2` — the timer slope —
at the aging endpoint `f = 0.5`. So sizer (`α=0`) and adder (`α=1`) stay homeostatic at
every age, but a **timer is driven to marginal loss of homeostasis precisely as division
symmetrizes** (`r = 2f → 1`): its birth-size CV, its RLS-CV broadening, and its size memory
all grow without bound together, because they are one eigenvalue.

## The combined guard

One testset (`combined — pole r ties CV-amplification, autocorrelation, and relaxation`)
runs a single lineage/config `(α=1.6, β=20, f=0.5, cv=0.06 ⇒ r=0.8)` and checks that the
three observables report the **same** pole to tolerance:

- `r_autocorr = ρ₁`;
- `r_cv = √(1 − (cv/CV(Vb))²)` from the CV-amplification law (intrinsic `cv` known);
- `r_step` = the geometric relaxation rate of a `β: 20→40` nutrient shift;

all equal `α·f = 0.8`, their relaxation memories `−1/ln r` agree with `size_memory`'s
`memory_gen`, and the invariant `√(CV(Vb)²(1−ρ₁²))` recovers the intrinsic noise `cv`. It is
one guard that the thread's observables are mutually consistent — if any drifts from the
shared pole, the map is not the AR(1) claimed.
