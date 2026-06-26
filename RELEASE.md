# Release checklist — Julia General registry

`CellSizeControl.jl` is intended for the Julia **General** registry (a citable, reusable
package backing the cell-size-control paper). Size is not a barrier — General has no minimum;
the criteria are a coherent reusable API + a clean, loadable, tested package, which this is.

## AutoMerge prerequisites (status)

- [x] **License** — MIT (`LICENSE`).
- [x] **`[compat]`** for all `[deps]` (Random, Statistics) + `julia`, and for the test
  extras (Aqua, ExplicitImports, Test).
- [x] **Loads + tests pass** — three-layer suite + Aqua + ExplicitImports green.
- [x] **Specific name** — "CellSizeControl" is descriptive (not a single generic word); check
  it's free in the registry at registration (`Pkg` name-clash check) — expected clear.
- [ ] **Public repo** — currently PRIVATE; flip to public on release (Daniel-gated).
- [ ] **Real version** — bump `0.1.0-dev` → `0.1.0` in `Project.toml` at release (drop the
  prerelease suffix; AutoMerge registers `0.1.0`).
- [ ] **Git tag** + run the JuliaRegistries/Registrator bot (`@JuliaRegistrator register`)
  on the tagged commit; let AutoMerge run (≈3-day waiting period for a new 0.x package).
- [ ] **(Recommended) docs + CI badge** — the README + `.github/workflows/ci.yml` are in
  place; a Documenter site is optional for a small package.

## Order of operations (all Daniel-gated)

1. Bump version to `0.1.0`, commit, push.
2. Make `github.com/djsegal/cell-size-control` public.
3. Tag `v0.1.0`; comment `@JuliaRegistrator register` on the release commit.
4. Mint the Zenodo DOI (the `.zenodo.json` deposit); cite the DOI + the registered package
   in the paper.

## Sibling note

`TranscriptionMultiplier` (in `tf-coupling-fit`) is the same-shape, already-public package
and is an equally good General candidate — same checklist (verify its `[compat]` + a real
tagged version). Registering both gives the two papers proper citable software artifacts.
