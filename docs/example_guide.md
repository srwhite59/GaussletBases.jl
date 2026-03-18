# Reading and running the examples

This page is the running guide for the examples.

It is not the full package manual. Its job is simpler:

- tell you where to start
- group the examples by purpose
- point you to the right status pages and supporting notes

From the repository root, each example runs as:

```bash
julia --project=. examples/NAME.jl
```

If you are not sure where to begin, first read:

- [`docs/index.md`](index.md)
- [`docs/first_radial_workflow.md`](first_radial_workflow.md)

## Core starting sequence

These are the best first examples for a new reader.

1. `01_first_gausslet.jl`  
   Build one ordinary gausslet and inspect its exact Gaussian expansion.
2. `02_radial_basis.jl`  
   Build a small demo radial basis and inspect the basic diagnostics.
3. `03_radial_operators.jl`  
   Build the basic radial one-body operators on a small demo basis.
4. `04_hydrogen_ground_state.jl`  
   Solve the radial hydrogen ground-state problem.

Those four examples are still the clearest public entry path.

## Radial and atomic sequence

If your main interest is radial atoms, continue here after the core starting
sequence.

Read first:

- [`docs/current_atomic_branch.md`](current_atomic_branch.md)
- [`docs/recommended_atomic_setup.md`](recommended_atomic_setup.md)

Then run:

5. `15_atomic_hydrogen_ylm.jl`  
   Add explicit `(l,m)` channels on top of the radial substrate.
6. `16_atomic_ida_ingredients.jl`  
   Build the static He / IDA-style atomic ingredients.
7. `19_atomic_ida_direct.jl`  
   Build the direct/Hartree term from a trial density.
8. `20_atomic_ida_exchange.jl`  
   Build the matching exchange term.
9. `21_atomic_ida_fock.jl`  
   Combine the one-body, direct, and exchange pieces into the current small
   Fock-style helper.
10. `22_atomic_ida_uhf.jl`  
    Run the present minimal UHF fixed-point kernel.

If you want the tiny exact-interacting checks after that:

11. `17_atomic_ida_two_electron.jl`
12. `18_atomic_ida_two_electron_lanczos.jl`

These are still intentionally narrow examples, not a broad atomic workflow.

## Ordinary Cartesian sequence

If your main interest is the ordinary Cartesian branch, use this sequence.

Read first:

- [`docs/current_ordinary_branch.md`](current_ordinary_branch.md)

Then run:

13. `23_cartesian_hydrogen_coulomb_expansion.jl`  
    First ordinary Cartesian hydrogen path: Coulomb expansion first, not grid
    first.
14. `24_mapped_cartesian_hydrogen.jl`  
    Move the same hydrogen problem onto the public mapped full-line basis
    route.
15. `25_mapped_cartesian_hydrogen_backends.jl`  
    Compare the numerical-reference and experimental PGDG-style one-body
    backends.
16. `29_hybrid_mapped_cartesian_hydrogen.jl`  
    Move into the friendlier hybrid/core-supported regime.
17. `30_ordinary_sho_spectra.jl`  
    Compare low-energy harmonic-oscillator spectra instead of only comparing
    raw matrix norms.

The examples below are supporting ordinary-branch diagnostics rather than the
first pages a new ordinary-branch reader should open:

18. `26_ordinary_cartesian_ida.jl`  
    Build the first static He-style ordinary Cartesian IDA ingredients.
19. `27_ordinary_cartesian_ida_localized_backends.jl`  
    Compare the backend split for the ordinary static ingredients.
20. `28_ordinary_one_body_fidelity.jl`  
    Decompose the remaining mild-mapped one-body differences.

This keeps the ordinary branch readable as:

- hydrogen first
- hybrid practical regime next
- backend diagnostics after that

## Primitive-layer and contraction sequence

If your main interest is the shared primitive/contraction architecture, use
this sequence after the core starting examples.

Read first:

- [`docs/intermediate_primitive_layer.md`](intermediate_primitive_layer.md)

Then run:

21. `05_primitive_sets.jl`
22. `06_basis_contraction.jl`
23. `07_position_contraction.jl`
24. `08_basis_representation.jl`
25. `14_radial_primitive_operators.jl`

That is the shortest path from “basis functions as final objects” to “basis
functions as contractions of a shared primitive layer.”

## Hierarchy and current corrected direction

If your main interest is the current hierarchy/contraction research direction,
read this only after the primitive-layer line is comfortable.

Read first:

- [`docs/architecture.md`](architecture.md)
- [`docs/global_map_local_contraction.md`](global_map_local_contraction.md)

Then run:

26. `09_basis_partition.jl`
27. `10_hierarchical_partition.jl`
28. `13_global_leaf_contraction.jl`

This is the clearest current research direction:

- one global mapped common layer
- optional local contraction on top of it

## Prototype side branch

These examples are still useful, but they are prototype studies rather than
the main public story:

29. `11_leaf_pgdg.jl`
30. `12_leaf_pgdg_augmentation.jl`

Read them only after the corrected-direction examples above.

## Supporting notes

For the current interpreted state of each main line, use:

- [`docs/current_atomic_branch.md`](current_atomic_branch.md)
- [`docs/current_ordinary_branch.md`](current_ordinary_branch.md)

For narrower supporting and historical decision notes, start from those branch
pages rather than opening the whole note stack at random.
