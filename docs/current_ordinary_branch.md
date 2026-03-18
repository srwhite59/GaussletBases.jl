# Current Ordinary Branch

This page is the shortest current status read for the ordinary Cartesian line.

## 1. What the ordinary branch is today

The current ordinary story is:

- Coulomb-expansion first, not 3D grid first
- one-dimensional mapped ordinary bases on each Cartesian axis
- an explicit backend split between validation and experimental PGDG-style
  analytic construction
- the friendlier hybrid/core-supported regime as the practical target

This line is still experimental, but it now has a clear current
interpretation.

## 2. Start here for the ordinary line

Read these first:

- [`docs/ordinary_coulomb_expansion_path.md`](ordinary_coulomb_expansion_path.md)
- [`docs/mapped_ordinary_basis.md`](mapped_ordinary_basis.md)
- [`docs/ordinary_pgdg_hybrid_regime.md`](ordinary_pgdg_hybrid_regime.md)
- [`docs/ordinary_sho_spectral_test.md`](ordinary_sho_spectral_test.md)
- [`docs/ordinary_pgdg_hybrid_consolidation.md`](ordinary_pgdg_hybrid_consolidation.md)

Then run:

1. `examples/23_cartesian_hydrogen_coulomb_expansion.jl`
2. `examples/24_mapped_cartesian_hydrogen.jl`
3. `examples/25_mapped_cartesian_hydrogen_backends.jl`
4. `examples/29_hybrid_mapped_cartesian_hydrogen.jl`
5. `examples/30_ordinary_sho_spectra.jl`

## 3. What counts as current workflow documentation

These are the current-status pages for the ordinary branch:

- [`docs/ordinary_coulomb_expansion_path.md`](ordinary_coulomb_expansion_path.md)
- [`docs/mapped_ordinary_basis.md`](mapped_ordinary_basis.md)
- [`docs/ordinary_pgdg_hybrid_regime.md`](ordinary_pgdg_hybrid_regime.md)
- [`docs/ordinary_sho_spectral_test.md`](ordinary_sho_spectral_test.md)
- [`docs/ordinary_pgdg_hybrid_consolidation.md`](ordinary_pgdg_hybrid_consolidation.md)

These pages describe the present recommended interpretation of the ordinary
branch.

## 4. Supporting notes for the ordinary line

These notes are supporting historical or development notes:

- [`docs/ordinary_pgdg_decision.md`](ordinary_pgdg_decision.md)
- [`docs/ordinary_pgdg_comx.md`](ordinary_pgdg_comx.md)
- [`docs/ordinary_pgdg_proxy_refinement.md`](ordinary_pgdg_proxy_refinement.md)
- [`docs/ordinary_pgdg_distortion_regime.md`](ordinary_pgdg_distortion_regime.md)
- [`docs/ordinary_pgdg_backend_pivot.md`](ordinary_pgdg_backend_pivot.md)
- [`docs/ordinary_cartesian_ida.md`](ordinary_cartesian_ida.md)
- [`docs/ordinary_pgdg_localized_backend.md`](ordinary_pgdg_localized_backend.md)
- [`docs/ordinary_pgdg_one_body_fidelity.md`](ordinary_pgdg_one_body_fidelity.md)

They record why the branch is interpreted the way it is now, but they are not
the first pages a new ordinary-branch reader should open.

## 5. Current interpretation

The current wording discipline for the ordinary line is:

- `:numerical_reference` remains the validation route
- the PGDG-style analytic route is good enough in the friendly hybrid regime
- hard pure mapped small-`c` cases remain stress tests
- `AsinhMapping` is the current working map, not final truth
- current `c,s` heuristics are provisional
- the radial branch remains numerical rather than PGDG-driven

That is much closer to the practical White–Lindsey-style hybrid picture than
to the harsher pure mapped stress-test cases.

If you want the broader package context after that, go back to:

- [`docs/index.md`](index.md)
- [`docs/example_guide.md`](example_guide.md)
