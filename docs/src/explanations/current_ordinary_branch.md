# Current ordinary branch

This page is the shortest current user-facing status read for the newer
ordinary Cartesian line.

## What the ordinary branch is today

The current ordinary story is:

- Coulomb-expansion first, not 3D grid first
- one-dimensional mapped ordinary bases on each Cartesian axis
- an explicit backend split between validation and experimental PGDG-style
  analytic construction
- a separate paper-faithful Qiu-White residual-Gaussian reference path
- the old COMX-cleaned 1D hybrid route quarantined as legacy/internal code,
  not current workflow

The Gaussian-supplement motivation is still the same as in the hybrid
gausslet/Gaussian papers: improve the near-nuclear/core region without asking
the gausslet backbone alone to resolve the whole core. But the supported repo
route for that comparison is now the separate full-3D Qiu-White-style path,
not the older 1D COMX-cleaned surrogate branch.

This line is still experimental, but it now has a clear current
interpretation and should be read as a real second workflow in the package.

For terms involving residual Gaussians, the intended interaction stays in the
same two-index integral-diagonal-approximation (IDA) representation used for
the gausslet channel, rather than being treated as a four-index tensor with
only diagonal index patterns retained.

## Start here

Within the manual, the best supporting pages are:

- [Examples](../examples/index.md)
- [Example guide](../howto/example_guide.md)
- [First radial workflow](../tutorials/first_radial_workflow.md)

The radial tutorial is still useful context because the package’s operator and
quadrature story was established there first.

## Notes that are not yet migrated into the first site

The current ordinary-branch workflow and supporting notes still live in the
flat `docs/` tree for now. The most important filenames are:

- `ordinary_coulomb_expansion_path.md`
- `mapped_ordinary_basis.md`
- `ordinary_cartesian_vee_validation.md`
- `ordinary_cartesian_qiu_white_reference.md`
- `ordinary_cartesian_qiu_white_crossblock_correction.md`
- `qw_pgdg_base_milestone_note.md`
- `qw_pgdg_fixed_a_mapping_note.md`
- `ordinary_pgdg_supporting_notes.md`
- `ordinary_homonuclear_chain_experimental_note.md`

The older 1D hybrid-surrogate notes remain in the flat `docs/` tree only as
legacy/internal supporting notes:

- `ordinary_cartesian_hybrid_vee_validation.md`
- `ordinary_cartesian_residual_gaussian_interaction.md`
- `ordinary_cartesian_mwg_interaction.md`
- `ordinary_cartesian_legacy_he_s_adapter.md`
- `ordinary_pgdg_hybrid_regime.md`
- `ordinary_sho_spectral_test.md`
- `ordinary_pgdg_hybrid_consolidation.md`

## Current interpretation

The wording discipline for the ordinary line remains:

- `:numerical_reference` is the validation route
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

The ordinary branch also now includes one narrow experimental homonuclear
chain line. Its current producer-side export surface is:

- `experimental_homonuclear_chain_nested_dense_payload`
- `write_experimental_homonuclear_chain_nested_dense_jld2`

That line should still be read as explicitly experimental:

- homonuclear chains only
- nested-chain policy still explicit rather than frozen
- producer-side continuation/export only
- not a settled general-purpose chain workflow

For the current status and milestone framing of that line, use:

- `docs/ordinary_homonuclear_chain_experimental_note.md`

If you want the lower-priority architecture and note-history context after
that, continue with:

- [Developer Notes](../developer/index.md)
