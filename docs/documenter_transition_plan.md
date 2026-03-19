# Documenter Transition Plan

This note records the first incremental step toward a more standard Julia
package documentation layout.

## 1. Current layout

The repository already has strong documentation content, but the current docs
system is still a flat-note system:

- many markdown notes directly under `docs/`
- no `docs/src/`
- no `docs/make.jl`
- no `docs/Project.toml`
- no Documenter-based API reference

That is workable for active research notes, but it is not the usual Julia
package documentation pattern.

## 2. Target layout

The target is the standard Documenter-based structure used by many Julia
packages:

- `docs/Project.toml`
- `docs/make.jl`
- `docs/src/`
- a small curated docs home page
- a later API reference built from docstrings

Within that layout, the current material should eventually separate into:

- tutorials / getting started
- how-to pages
- explanations
- API reference from docstrings
- development and history notes

## 3. Why this should be incremental

The repository already has a large amount of scientifically useful text. A
single rewrite would risk mixing three jobs at once:

- rebuilding the docs toolchain
- rewriting the teaching flow
- rewriting the scientific content

That is not necessary.

The better first transition step is:

1. add the standard Documenter skeleton
2. create a small curated site from the strongest current pages
3. keep the flat note tree in place for now
4. migrate or merge the supporting notes only later

This keeps the package readable while avoiding a disruptive docs rewrite.

## 4. First curated site

The first Documenter site should stay intentionally small. A good first page
set is:

- home / index
- first radial workflow
- recommended atomic setup
- example guide
- current atomic branch
- current ordinary branch
- architecture

That is enough to prove the package can build a standard Julia docs site
without pretending that the full note tree has already been migrated.

## 5. Later split of the flat notes

The eventual split should look roughly like this:

- tutorials:
  - first radial workflow
  - first ordinary-branch workflow when that becomes stable enough
- how-to:
  - recommended atomic setup
  - example running guides
  - export usage notes
- explanations:
  - branch-status pages
  - architecture
  - terminology
- API reference:
  - pages built from real docstrings on exported symbols
- development/history notes:
  - supporting note chains that explain how current interpretations were reached

## 6. Exported symbols that most need docstrings later

The exported symbols that most urgently need real API docstrings for a later
reference pass are:

- basis construction and mappings:
  - `build_basis`
  - `UniformBasisSpec`
  - `MappedUniformBasisSpec`
  - `HalfLineBasisSpec`
  - `RadialBasisSpec`
  - `AsinhMapping`
  - `fit_asinh_mapping_for_extent`
  - `fit_asinh_mapping_for_strength`
- basic operator and diagnostics layer:
  - `basis_diagnostics`
  - `radial_quadrature`
  - `overlap_matrix`
  - `kinetic_matrix`
  - `nuclear_matrix`
  - `centrifugal_matrix`
  - `multipole_matrix`
- atomic line:
  - `atomic_operators`
  - `atomic_one_body_operators`
  - `atomic_ida_operators`
  - `direct_matrix`
  - `exchange_matrix`
  - `fock_matrix_alpha`
  - `fock_matrix_beta`
  - `uhf_scf`
- ordinary mapped line:
  - `hybrid_mapped_ordinary_basis`
  - `mapped_ordinary_one_body_operators`
  - `mapped_cartesian_hydrogen_energy`
  - `ordinary_sho_hamiltonian`
  - `ordinary_sho_spectrum`
  - `ordinary_cartesian_ida_operators`
- export layer:
  - `write_fullida_dense_jld2`
  - `write_sliced_ham_jld2`

Those are the exported entry points most likely to matter once the package has
a real API reference page.
