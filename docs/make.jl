#!/usr/bin/env julia
# Minimal Documenter.jl scaffold for CellSizeControl.
#
# Build (standalone) from the package directory:
#
#   julia --project=docs -e 'using Pkg; Pkg.instantiate()'   # once
#   julia --project=docs docs/make.jl
#
# Output goes to docs/build/. Doctests embedded in the docstrings run as part of
# the build (doctest = true). The HTML is not deployed; this scaffold exists so
# the exported API is documented and checkable.
#
# ADDITIVE: lives in its own environment (docs/Project.toml) and loads the
# package from the parent directory via LOAD_PATH; touches no existing files.

push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))

using Documenter
using CellSizeControl

DocMeta.setdocmeta!(CellSizeControl, :DocTestSetup,
                    :(using CellSizeControl); recursive = true)

makedocs(
    sitename = "CellSizeControl.jl",
    modules  = [CellSizeControl],
    authors  = "Daniel J. Segal",
    pages    = ["Home" => "index.md"],
    # Run doctests during the build; fail on broken/missing cross-references.
    doctest  = true,
    checkdocs = :exports,
    format   = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true"),
)
