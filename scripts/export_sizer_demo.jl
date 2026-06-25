# V4 — export data for the cell-size-control demo figure:
#  (1) lineages: a sub-doubling TIMER collapses while the inhibitor-dilution SIZER
#      holds birth size steady (the core size-control result / the course bug fix);
#  (2) the Soifer-Amir Vd-vs-Vb scatter for timer/adder/sizer regimes (slope 2/1/0).
# Run: julia --project=. scripts/export_sizer_demo.jl   (writes docs/figdata/*.csv)
using CellSizeControl

root = normpath(joinpath(@__DIR__, ".."))
fd = joinpath(root, "docs", "figdata")
mkpath(fd)

# (1) lineages over generations (deterministic, cv=0)
ngen = 25
timer = simulate_lineage(TimerRule(1.6); V0=20.0, n=ngen, cv=0.0, seed=1)        # collapses
sizer = simulate_lineage(InhibitorDilutionSizer(60.0, 1.5); V0=5.0, n=ngen, cv=0.0, seed=2)  # V*=40 → Vb→20
open(joinpath(fd, "lineages.csv"), "w") do io
    println(io, "gen,timer_Vb,sizer_Vb")
    for i in 1:ngen
        println(io, "$i,$(timer.Vb[i]),$(sizer.Vb[i])")
    end
end

# (2) Vd-vs-Vb scatter + fitted slope per regime
regimes = (
    ("timer", TimerRule(2.0), 1),
    ("adder", AdderRule(20.0), 2),
    ("sizer", SizerRule(40.0), 3),
)
open(joinpath(fd, "scatter.csv"), "w") do io
    println(io, "regime,slope,Vb,Vd")
    for (name, rule, seed) in regimes
        s = simulate_lineage(rule; V0=20.0, n=400, cv=0.12, seed=seed)
        sl = size_control_slope(s.Vb, s.Vd)
        for i in 1:length(s.Vb)
            println(io, "$name,$sl,$(s.Vb[i]),$(s.Vd[i])")
        end
    end
end

println("wrote ", joinpath(fd, "lineages.csv"), " and scatter.csv")
println(
    "  timer last/first Vb = ",
    round(timer.Vb[end] / timer.Vb[1]; digits=4),
    " (collapse) ; sizer last Vb = ",
    round(sizer.Vb[end]; digits=3),
    " (→ V*/2=20)",
)
