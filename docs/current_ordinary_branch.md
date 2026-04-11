# Current Ordinary Branch

This page is the shortest current status read for the newer ordinary Cartesian
line.

## 1. What the ordinary branch is today

The current ordinary story is:

- Coulomb-expansion first, not 3D grid first
- one-dimensional mapped ordinary bases on each Cartesian axis
- an explicit backend split between validation and experimental PGDG-style
  analytic construction
- a separate paper-faithful Qiu-White residual-Gaussian reference path
- the old COMX-cleaned 1D hybrid route quarantined as legacy/internal code,
  not current workflow

This line is still experimental, but it now has a clear current
interpretation and should be read as a real second workflow in the package.

## 2. Start here for the ordinary line

Read these first:

- [`docs/ordinary_coulomb_expansion_path.md`](ordinary_coulomb_expansion_path.md)
- [`docs/mapped_ordinary_basis.md`](mapped_ordinary_basis.md)
- [`docs/ordinary_cartesian_vee_validation.md`](ordinary_cartesian_vee_validation.md)
- [`docs/ordinary_cartesian_qiu_white_reference.md`](ordinary_cartesian_qiu_white_reference.md)
- [`docs/ordinary_cartesian_qiu_white_crossblock_correction.md`](ordinary_cartesian_qiu_white_crossblock_correction.md)
- [`docs/qw_pgdg_base_milestone_note.md`](qw_pgdg_base_milestone_note.md)
- [`docs/qw_pgdg_fixed_a_mapping_note.md`](qw_pgdg_fixed_a_mapping_note.md)

Then run:

1. `examples/23_cartesian_hydrogen_coulomb_expansion.jl`
2. `examples/24_mapped_cartesian_hydrogen.jl`
3. `examples/25_mapped_cartesian_hydrogen_backends.jl`
4. `examples/33_ordinary_cartesian_1s2_vee.jl`
5. `examples/38_qiu_white_reference_vee.jl`

The older 1D COMX-cleaned hybrid examples remain in `examples/` only as
legacy/internal experimental regressions:

- `29_hybrid_mapped_cartesian_hydrogen.jl`
- `30_ordinary_sho_spectra.jl`
- `34_hybrid_cartesian_1s2_vee.jl`
- `35_hybrid_cartesian_residual_vee.jl`
- `36_hybrid_cartesian_legacy_he_s_vee.jl`
- `37_hybrid_cartesian_mwg_vee.jl`

## 3. What counts as current workflow documentation

These are the current-status pages for the ordinary branch:

- [`docs/ordinary_coulomb_expansion_path.md`](ordinary_coulomb_expansion_path.md)
- [`docs/mapped_ordinary_basis.md`](mapped_ordinary_basis.md)
- [`docs/ordinary_cartesian_vee_validation.md`](ordinary_cartesian_vee_validation.md)
- [`docs/ordinary_cartesian_qiu_white_reference.md`](ordinary_cartesian_qiu_white_reference.md)
- [`docs/ordinary_cartesian_qiu_white_crossblock_correction.md`](ordinary_cartesian_qiu_white_crossblock_correction.md)

These pages describe the present recommended interpretation of the ordinary
branch.

The older 1D hybrid-surrogate notes remain in the flat `docs/` tree as
legacy/internal supporting notes, not as current workflow pages:

- `ordinary_cartesian_hybrid_vee_validation.md`
- `ordinary_cartesian_residual_gaussian_interaction.md`
- `ordinary_cartesian_mwg_interaction.md`
- `ordinary_cartesian_legacy_he_s_adapter.md`
- `ordinary_pgdg_hybrid_regime.md`
- `ordinary_sho_spectral_test.md`
- `ordinary_pgdg_hybrid_consolidation.md`

## 4. Supporting notes for the ordinary line

The main supporting-note stack is now grouped here:

- [`docs/ordinary_pgdg_supporting_notes.md`](ordinary_pgdg_supporting_notes.md)

One additional narrow supporting note remains useful on top of that:

- [`docs/ordinary_cartesian_ida.md`](ordinary_cartesian_ida.md)

These notes record why the branch is interpreted the way it is now, but they
are not the first pages a new ordinary-branch reader should open.

## 5. Current interpretation

The current wording discipline for the ordinary line is:

- `:numerical_reference` remains the validation route
- the PGDG-style analytic route is good enough on the mapped ordinary backbone
- Gaussian-supplement comparisons should go through the separate paper-faithful
  3D Qiu-White route, not the old 1D COMX-cleaned surrogate path
- hard pure mapped small-`c` cases remain stress tests
- `AsinhMapping` is the current working map, not final truth
- for pre-nesting Cartesian convergence tests, the current default family is
  fixed `a = 1/(2Z)` with `s` solved from `count` and `xmax` rather than a
  fixed-`s` scan
- the radial branch remains numerical rather than PGDG-driven

That keeps the public ordinary workflow aligned with the corrected 3D
gausslet-plus-Gaussian interpretation, without treating the older 1D
COMX-cleaned hybrid surrogate as a supported route.

If you want the broader package context after that, go back to:

- [`docs/index.md`](index.md)
- [`docs/example_guide.md`](example_guide.md)
