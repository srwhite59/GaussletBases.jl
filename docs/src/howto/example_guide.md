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

## Core starting sequence

If you are new, start here:

1. `01_first_gausslet.jl`
2. `02_radial_basis.jl`
3. `03_radial_operators.jl`
4. `04_hydrogen_ground_state.jl`

Those four examples are still the clearest public entry path.

## Radial and atomic sequence

After the core starting sequence, the atomic line is:

5. `15_atomic_hydrogen_ylm.jl`
6. `16_atomic_ida_ingredients.jl`
7. `19_atomic_ida_direct.jl`
8. `20_atomic_ida_exchange.jl`
9. `21_atomic_ida_fock.jl`
10. `22_atomic_ida_uhf.jl`

The tiny exact-interacting checks are:

11. `17_atomic_ida_two_electron.jl`
12. `18_atomic_ida_two_electron_lanczos.jl`

The current export examples are:

13. `31_atomic_fullida_dense_export.jl`
14. `32_atomic_sliced_export.jl`

## Ordinary Cartesian sequence

For the ordinary Cartesian line, use this sequence:

15. `23_cartesian_hydrogen_coulomb_expansion.jl`
16. `24_mapped_cartesian_hydrogen.jl`
17. `25_mapped_cartesian_hydrogen_backends.jl`
18. `33_ordinary_cartesian_1s2_vee.jl`
19. `29_hybrid_mapped_cartesian_hydrogen.jl`
20. `30_ordinary_sho_spectra.jl`

The more diagnostic ordinary examples are:

21. `26_ordinary_cartesian_ida.jl`
22. `27_ordinary_cartesian_ida_localized_backends.jl`
23. `28_ordinary_one_body_fidelity.jl`

## Primitive and hierarchy sequence

If your focus is the primitive/contraction architecture, continue with:

24. `05_primitive_sets.jl`
25. `06_basis_contraction.jl`
26. `07_position_contraction.jl`
27. `08_basis_representation.jl`
28. `14_radial_primitive_operators.jl`
29. `09_basis_partition.jl`
30. `10_hierarchical_partition.jl`
31. `13_global_leaf_contraction.jl`

The prototype side branch remains:

32. `11_leaf_pgdg.jl`
33. `12_leaf_pgdg_augmentation.jl`

## Read these pages alongside the examples

Use these pages to keep the examples in context:

- [First radial workflow](../tutorials/first_radial_workflow.md)
- [Recommended atomic setup](recommended_atomic_setup.md)
- [Current atomic branch](../explanations/current_atomic_branch.md)
- [Current ordinary branch](../explanations/current_ordinary_branch.md)
- [Developer Notes](../developer/index.md)
