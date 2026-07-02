# Mixed dilution + titration extension of the inhibitor-dilution sizer.
# The pure-dilution law V* = W/c* is exactly proportional in Whi5 dose W (fold = 2 for 2x W).
# Schmoller 2015 instead reports a SUB-proportional dose response (2x WHI5 -> ~1.40x size in
# haploids, ~1.18x in diploids). A minimal combined model adds a titration baseline V0 (a
# Cln3-vs-fixed-SBF-site term, ploidy-scaled: V0 propto genome copies g): V*(W,g) = g*V0 + W/c*.
# This is closed-form (analytically tractable) and makes the dose response sub-proportional, with a
# larger baseline (more sites) giving a smaller fold -- the observed ploidy direction. We report the
# fold for a doubling of Whi5 dose across the baseline ratio rho = c* V0 / W_ref, and the rho that
# reproduces the haploid 1.40 fold, then the diploid fold that the same c* V0 (doubled sites) predicts.
# Run: julia --project=. analysis/gen_mixed_dilution_titration.jl
using Printf

# V*(W,g) = g*V0 + W/c* ; in units of W_ref/c*: v(w,g) = g*rho + w   (rho = c* V0 / W_ref, w = W/W_ref)
fold(rho, g) = (g*rho + 2.0) / (g*rho + 1.0)   # fold in size for doubling the dilution term (w: 1 -> 2)

@printf("pure dilution (rho=0): haploid fold = %.2f (proportional, the minimal-switch prediction)\n", fold(0.0, 1))
@printf("\ntitration baseline sweep (haploid, g=1):\n  rho   fold(2x W)\n")
for rho in (0.0, 0.5, 1.0, 1.5, 2.0, 3.0)
    @printf("  %.1f    %.2f\n", rho, fold(rho, 1))
end

# solve haploid fold(rho,1) = 1.40  ->  (rho+2)/(rho+1) = 1.40  ->  rho = (2 - 1.40)/(1.40 - 1)
rho_h = (2 - 1.40) / (1.40 - 1)
@printf("\nrho that reproduces the haploid 1.40 fold: rho = %.2f\n", rho_h)
@printf("  -> haploid fold = %.2f (target 1.40)\n", fold(rho_h, 1))

# if the titration baseline scales with genome copies (diploid g=2, twice the SBF sites), the SAME
# per-copy rho predicts the diploid fold: a larger baseline -> smaller fold, the observed direction.
@printf("  -> diploid fold (g=2, same per-copy rho) = %.2f (Schmoller ~1.18)\n", fold(rho_h, 2))
@printf("\nDirection captured: baseline makes the response sub-proportional (fold<2) and the larger\n")
@printf("diploid baseline lowers the fold further, matching haploid>diploid. One extra parameter (rho),\n")
@printf("closed-form. (Precise two-fold values are not jointly fit: only %d digitized points exist.)\n", 4)

# emit the mixed-model folds for the Fig 3 overlay (reproducible; read by plot_whi5_dosage_test.py)
here = @__DIR__
open(joinpath(here, "whi5_dosage_mixed.csv"), "w") do io
    println(io, "series,ploidy,fold,rho")
    @printf(io, "mixed,1,%.4f,%.4f\n", fold(rho_h, 1), rho_h)
    @printf(io, "mixed,2,%.4f,%.4f\n", fold(rho_h, 2), rho_h)
end
println("wrote whi5_dosage_mixed.csv")
