# Reading and running the examples

The examples in this repository are meant to be read in sequence.

They do not all serve the same purpose. Some are beginner examples, some are the main radial scientific path, and some belong to the more experimental contraction and hierarchy line.

From the repository root, each example runs as:

```bash
julia --project=. examples/NAME.jl
```

## Stage 1: first contact

These are the best starting examples for a new user.

### `01_first_gausslet.jl`

Build one ordinary gausslet, evaluate it, and inspect its exact Gaussian expansion.

This is the smallest possible introduction to the package.

### `02_radial_basis.jl`

Build a **small demo radial basis**, inspect one basis function, and evaluate the basic diagnostics.

This example is deliberately tiny and fast. It is a smoke-test-style introduction, not the package’s recommended atom setup.

### `03_radial_operators.jl`

Build the basic radial one-body operators and the current radial operator bundle on a **small demo basis**.

Again, this is a fast introductory example rather than a production recommendation.

### `04_hydrogen_ground_state.jl`

Solve the hydrogen ground-state problem in the radial basis.

This is the first full scientific check in the repository, and it is the most important example after the first three.

### `23_cartesian_hydrogen_coulomb_expansion.jl`

Build the first small ordinary Cartesian hydrogen Hamiltonian by:

- taking a short Gaussian expansion of `1/r`
- building one-dimensional Gaussianized operator factors
- assembling the 3D nuclear-attraction operator from separable pieces

This is the right first ordinary-gausslet hydrogen path. It should be read
before any future 3D grid-based debug or comparison path.

### `24_mapped_cartesian_hydrogen.jl`

Build the same ordinary Cartesian hydrogen Hamiltonian, but now on the public
globally mapped full-line basis route.

This is the packaged form of the successful Asinh-mapped hydrogen study:

- one global map on each Cartesian axis
- one-dimensional mapped ordinary bases
- Coulomb expansion first
- no 3D grid-first design

### `25_mapped_cartesian_hydrogen_backends.jl`

Build the same mapped Cartesian hydrogen problem with the explicit backend
split now used on the one-body mapped ordinary branch:

- `:pgdg_experimental` as the preferred analytic route in the mild/moderate
  regime
- `:numerical_reference` as the validation route

This is the example that makes the present mapped ordinary status explicit.

### `15_atomic_hydrogen_ylm.jl`

Take the same hydrogen problem one step further by adding explicit angular channels `(l,m)` on top of the radial operator substrate.

This is the first atomic example in the repository in the usual radial-plus-angular sense.

### `16_atomic_ida_ingredients.jl`

Build the static He-style IDA ingredients on top of the same radial-plus-angular structure.

This example does **not** solve the many-electron problem. It makes the one-body blocks, radial multipoles, angular kernels, and orbital indexing explicit.

### `19_atomic_ida_direct.jl`

Build the first direct/Hartree one-body term from a trial spatial density matrix.

This is the first small physical consumer of the sectorized angular preparation.

### `20_atomic_ida_exchange.jl`

Build the matching exchange/Fock-style one-body term from a trial spatial density matrix.

This completes the first small direct-plus-exchange mean-field-style pair, still without a full SCF workflow.

### `21_atomic_ida_fock.jl`

Combine the one-body atomic term with the direct and exchange pieces into the first tiny Fock-style helper.

This is still not an SCF driver. It is only the effective one-body matrix assembly step.

### `22_atomic_ida_uhf.jl`

Run the smallest self-consistent He-like UHF iteration in the current atomic IDA model.

This is still intentionally narrow:

- fixed occupations
- simple damping
- no DIIS
- no broad workflow layer
- no claim of a general HF implementation

### `17_atomic_ida_two_electron.jl`

Use those same IDA ingredients in the smallest real interacting application: one spin-up and one spin-down electron in a tiny atomic basis.

This is still intentionally small and explicit. It is not yet a general many-electron workflow.

### `18_atomic_ida_two_electron_lanczos.jl`

Take a somewhat larger He-like model and solve it with a standard Hermitian Lanczos iteration.

This is the first step beyond the tiny dense reference problem, but it is still a narrow demonstration rather than a broad many-electron framework.

## Stage 2: primitive layers and contraction

These examples explain the common Gaussian primitive layer that lies behind the basis functions.

### `05_primitive_sets.jl`

Build primitive sets directly and form simple matrices on them.

### `06_basis_contraction.jl`

Start from a basis, recover its primitive layer and contraction matrix, and rebuild basis-level matrices by contracting primitive matrices upward.

### `07_position_contraction.jl`

Do the same thing for the position matrix.

### `08_basis_representation.jl`

Build a compact in-memory representation of a basis together with its primitive layer and selected matrices.

These examples are for users who want to understand the package at a deeper level than “build a basis and call operators.”

### `14_radial_primitive_operators.jl`

Build the primitive-space radial one-body operators on an explicit quadrature grid, contract them upward to the radial basis, and compare them with the current direct radial operator path.

This is the clean bridge between the radial workflow and the shared primitive/contraction architecture.

## Stage 3: partitions and hierarchy

These examples introduce local grouping in physical space.

### `09_basis_partition.jl`

Partition basis functions into interval boxes using their physical-space centers, and inspect local blocks and couplings.

### `10_hierarchical_partition.jl`

Refine one box and inspect the resulting parent-child hierarchy.

These are advanced organizational examples. They do not yet build a new basis.

## Stage 4: corrected current direction

These are the examples to read only after the earlier structure is clear.
They express the current corrected nested/contraction direction more faithfully than the earlier leaf-local prototype line.

### `13_global_leaf_contraction.jl`

Build one globally mapped common basis over a region, then contract locally on the leaf boxes of a hierarchy.

This is the example that best matches the current corrected research direction:

- one global mapped common layer
- optional local contraction as a refinement

## Stage 5: prototype side branch

These examples are still useful, but they should be read as prototypes rather than as the main conceptual direction.

### `11_leaf_pgdg.jl`

Generate a simple leaf-local basis on a hierarchy.

This is a useful prototype, but it is not the best conceptual picture of the long-term direction.

### `12_leaf_pgdg_augmentation.jl`

Add extra user-supplied Gaussian primitives in selected leaves.

This is another prototype example showing local enrichment.

## Recommended reading orders

### If your main interest is radial atomic work

Read:

1. `01_first_gausslet.jl`
2. `02_radial_basis.jl`
3. `03_radial_operators.jl`
4. `04_hydrogen_ground_state.jl`
5. `15_atomic_hydrogen_ylm.jl`

Then, if you want the first ordinary Cartesian hydrogen path, read:

6. `23_cartesian_hydrogen_coulomb_expansion.jl`
7. `24_mapped_cartesian_hydrogen.jl`
8. `25_mapped_cartesian_hydrogen_backends.jl`

Then read:

- `docs/first_radial_workflow.md`
- `docs/recommended_atomic_setup.md`
- `docs/atomic_ylm_layer.md`
- `docs/ordinary_coulomb_expansion_path.md`
- `docs/mapped_ordinary_basis.md`
- `docs/ordinary_pgdg_backend_pivot.md`

### If your main interest is primitive layers and contraction

After the radial examples, continue with:

6. `05_primitive_sets.jl`
7. `06_basis_contraction.jl`
8. `07_position_contraction.jl`
9. `08_basis_representation.jl`
10. `14_radial_primitive_operators.jl`

### If your main interest is the current nested/hierarchy research direction

Only after the earlier stages, continue with:

11. `09_basis_partition.jl`
12. `10_hierarchical_partition.jl`
13. `13_global_leaf_contraction.jl`

Then, only if you want the prototype side branch, treat:

14. `11_leaf_pgdg.jl`
15. `12_leaf_pgdg_augmentation.jl`

as prototype studies rather than as the main public story of the package.

### If your main interest is the current small atomic IDA / HF line

After `15_atomic_hydrogen_ylm.jl`, continue with:

16. `16_atomic_ida_ingredients.jl`
17. `19_atomic_ida_direct.jl`
18. `20_atomic_ida_exchange.jl`
19. `21_atomic_ida_fock.jl`
20. `22_atomic_ida_uhf.jl`
21. `17_atomic_ida_two_electron.jl`
22. `18_atomic_ida_two_electron_lanczos.jl`

Then read:

- `docs/atomic_ida_layer.md`
- `docs/atomic_ida_direct.md`
- `docs/atomic_ida_exchange.md`
- `docs/atomic_ida_fock.md`
- `docs/atomic_ida_spin_fock.md`
- `docs/atomic_ida_uhf.md`
- `docs/atomic_ida_two_electron.md`

The important interpretation here is:

- `16` gives the static interacting ingredients
- `19` and `20` give the direct and exchange pieces
- `21` gives the algebraic Fock-style combination
- `22` gives the first minimal UHF fixed-point kernel in the current atomic IDA model

This is already a meaningful atomic mean-field line, but it should still be read
as a small current-model workflow rather than as a general HF package.
