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
    pages    = ["Home" => "index.md", "Showcase" => ["showcase/discriminator.md", "showcase/sizer_forgets.md", "showcase/whi5_dose.md"]],
    # Run doctests during the build; fail on broken/missing cross-references.
    doctest  = true,
    checkdocs = :exports,
    format   = Documenter.HTML(; prettyurls = false),
)

# Publish to GitHub Pages (gh-pages branch). A no-op when run locally; on CI it
# deploys via GITHUB_TOKEN / DOCUMENTER_KEY (active once the repo is public).
deploydocs(;
    repo      = "github.com/djsegal/CellSizeControl.jl",
    devbranch = "main",
)
