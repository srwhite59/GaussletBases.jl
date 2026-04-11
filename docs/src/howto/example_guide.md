# Example guide

This page is the running guide for the examples.

For the underlying API entry points used by these examples, also see:

- [Bases and mappings](../reference/bases_and_mappings.md)
- [Operators and diagnostics](../reference/operators_and_diagnostics.md)
- [Atomic and ordinary workflows](../reference/atomic_and_ordinary.md)
- [Export layer](../reference/export.md)

From the repository root, each example runs as:

```bash
julia --project=. examples/NAME.jl
```

From a fresh checkout, instantiate once first:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Core starting sequence

If you are new, start here:

1. `01_first_gausslet.jl`
2. `02_radial_basis.jl`
3. `03_radial_operators.jl`
4. `04_hydrogen_ground_state.jl`

Those four examples are still the clearest public entry path.

For the radial and atomic examples below, first runs can still be compile-heavy.
Some of the heavier examples use cheaper quadrature settings where that is
enough to keep the workflow clear. Reruns are much faster.

## Radial and atomic sequence

After the core starting sequence, the atomic line is:

5. `15_atomic_hydrogen_ylm.jl`
6. `16_atomic_ida_ingredients.jl`
7. `19_atomic_ida_direct.jl`
8. `20_atomic_ida_exchange.jl`
9. `21_atomic_ida_fock.jl`
10. `22_atomic_ida_uhf.jl`

The key conceptual jump happens between:

- `04_hydrogen_ground_state.jl`, which is still a single fixed-`l` radial
  block
- and `15_atomic_hydrogen_ylm.jl`, which introduces the explicit `(l,m)`
  channel list and the full one-electron atomic block structure

The tiny exact-interacting checks are:

11. `17_atomic_ida_two_electron.jl`
12. `18_atomic_ida_two_electron_lanczos.jl`

The current export examples are:

13. `31_atomic_fullida_dense_export.jl`
14. `32_atomic_sliced_export.jl`

These are now explicit export demonstrations, not just raw smoke tests: they
check a few structural invariants after writing the file. They are still
heavier than the earlier atomic examples and are best treated as advanced
workflow examples. They also use a deliberately simpler radial setup than the
default two-`xgaussian` atomic front door so the export layer itself stays the
focus.

## Ordinary Cartesian sequence

For the ordinary Cartesian line, use this sequence:

15. `23_cartesian_hydrogen_coulomb_expansion.jl`
16. `24_mapped_cartesian_hydrogen.jl`
17. `25_mapped_cartesian_hydrogen_backends.jl`
18. `33_ordinary_cartesian_1s2_vee.jl`
19. `38_qiu_white_reference_vee.jl`

`38_qiu_white_reference_vee.jl` is a slow reference example. Its nearest/GGT
path is part of the public ordinary workflow; the MWG branch is still
experimental and is skipped unless you opt in with
`GAUSSLETBASES_RUN_EXPERIMENTAL_MWG=1`.

The older 1D COMX-cleaned hybrid examples remain in `examples/` only as
legacy/internal experimental regressions and are intentionally omitted from
the supported ordinary sequence:

- `29_hybrid_mapped_cartesian_hydrogen.jl`
- `30_ordinary_sho_spectra.jl`
- `34_hybrid_cartesian_1s2_vee.jl`
- `35_hybrid_cartesian_residual_vee.jl`
- `36_hybrid_cartesian_legacy_he_s_vee.jl`
- `37_hybrid_cartesian_mwg_vee.jl`

The more diagnostic ordinary examples are:

20. `26_ordinary_cartesian_ida.jl`
21. `27_ordinary_cartesian_ida_localized_backends.jl`
22. `28_ordinary_one_body_fidelity.jl`

There is also one narrow experimental homonuclear-chain export line. It does
not yet have a full public example script, but the current shape is:

```julia
basis = bond_aligned_homonuclear_chain_qw_basis(...)
path = experimental_bond_aligned_homonuclear_chain_nested_qw_operators(basis)
write_experimental_homonuclear_chain_nested_dense_jld2("chain.jld2", path)
```

Treat that line as:

- homonuclear-chain-specific
- producer-side only
- explicitly experimental

For the current milestone/status framing, see
`docs/ordinary_homonuclear_chain_experimental_note.md`.

## Primitive and hierarchy sequence

If your focus is the primitive/contraction architecture, continue with:

29. `05_primitive_sets.jl`
30. `06_basis_contraction.jl`
31. `07_position_contraction.jl`
32. `08_basis_representation.jl`
33. `14_radial_primitive_operators.jl`
34. `09_basis_partition.jl`
35. `10_hierarchical_partition.jl`
36. `13_global_leaf_contraction.jl`

The prototype side branch remains:

37. `11_leaf_pgdg.jl`
38. `12_leaf_pgdg_augmentation.jl`

## Read these pages alongside the examples

Use these pages to keep the examples in context:

- [First radial workflow](../tutorials/first_radial_workflow.md)
- [Recommended atomic setup](recommended_atomic_setup.md)
- [Visualization utilities](visualization.md)
- [Current atomic branch](../explanations/current_atomic_branch.md)
- [Current ordinary branch](../explanations/current_ordinary_branch.md)
- [Developer Notes](../developer/index.md)
