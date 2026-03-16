# GaussletBases.jl

GaussletBases.jl is a Julia package for working with gausslet basis functions.

Its broad foundation is the ordinary one-dimensional gausslet story:

- explicit Gaussian constructions
- mapped and distorted coordinates
- basis functions built from a visible primitive layer
- matrix construction through that shared primitive layer

Today its most mature public-facing path is **atom-centered radial calculations**: building radial bases, constructing matching quadrature grids, checking basis diagnostics, and forming the basic radial operator matrices needed for one-electron and mean-field-style work.

The package now also has a first explicit atomic angular layer built from spherical-harmonic channels `(l,m)` on top of that radial substrate. That makes it possible to assemble and diagonalize the one-electron hydrogen Hamiltonian in an explicit radial-plus-angular basis without yet adding the full two-electron atomic story.

The repository also contains a newer one-dimensional line of work on **shared primitive layers, contraction, partitions, and hierarchy**. That line is scientifically important because it is the bridge from ordinary gausslets toward nested constructions, but it is still the more experimental side of the package. If you are new, start with the radial path, while keeping in mind that ordinary gausslets are the broader foundation underneath it.

## Who this package is for

GaussletBases is most useful today for people who want to:

- explore ordinary gausslets as localized basis functions built from explicit Gaussian primitive layers
- explore radial gausslet bases for atoms and related model problems
- inspect the underlying Gaussian layer behind a basis
- study how primitive Gaussian layers can be contracted into more useful localized functions
- experiment with simple one-dimensional hierarchy and contraction ideas without yet committing to a full molecular workflow

It is **not** yet a complete electronic-structure workflow package. It is a basis, quadrature, and operator package, with a newer experimental line for contraction and hierarchy in 1D.

## Installation

At present, install the package directly from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/srwhite59/GaussletBases.jl")
```

Then load it with:

```julia
using GaussletBases
```

## What is solid today

For a new user, the package currently has three visible layers:

- the broad foundation of ordinary 1D gausslets and explicit primitive-layer constructions
- the mature radial line
- the newer experimental contraction and hierarchy line

The most established public-facing path is the radial line:

- ordinary 1D, half-line, and radial gausslet bases
- explicit coordinate mappings
- explicit radial quadrature, separate from the basis
- basis diagnostics
- radial one-body matrices
- a first explicit one-electron atomic `(l,m)` layer built on those radial matrices
- the current two-index IDA-style radial multipole matrices
- a high-level `RadialAtomicOperators` bundle

The newer 1D contraction line is also real code, but it should be treated as **advanced** rather than as the default starting point.

## Best first path through the repository

If you are new to the package, a good first path is:

1. run the first four examples
2. read the radial quickstart
3. read the recommended atomic setup note
4. only then move on to primitive layers, contraction, and hierarchy

From the repository root:

```bash
julia --project=. examples/01_first_gausslet.jl
julia --project=. examples/02_radial_basis.jl
julia --project=. examples/03_radial_operators.jl
julia --project=. examples/04_hydrogen_ground_state.jl
```

Those four examples take you from one gausslet, to a radial basis, to radial operators, to a real hydrogen ground-state calculation.

The next natural step after those four is:

```bash
julia --project=. examples/15_atomic_hydrogen_ylm.jl
```

which adds the explicit `(l,m)` angular channels on top of the same radial substrate.

## A good first atom-centered setup

For an atom of nuclear charge `Z`, the package documentation recommends starting with:

- `s = 0.2`
- `c = s / (2Z)`
- `tails = 6`
- `odd_even_kmax = 6`
- `rmax = 30.0` bohr as a good first-row starting point

In code:

```julia
using GaussletBases

Z = 2.0
s = 0.2

map = AsinhMapping(c = s / (2Z), s = s)

rb = build_basis(RadialBasisSpec(:G10;
    rmax = 30.0,
    mapping = map,
    reference_spacing = 1.0,
    tails = 6,
    odd_even_kmax = 6,
    xgaussians = XGaussian[],
))
```

This is the recommended package starting point, not the only possible one. In practice, somewhat larger values of `s`, such as `0.3` to `0.5`, can still work surprisingly well in some situations. But if you want one clean starting recipe, use `s = 0.2`.

## What you usually do next

Once you have a radial basis, the normal next steps are:

```julia
diag = basis_diagnostics(rb)
grid = radial_quadrature(rb)
ops = atomic_operators(rb, grid; Z = Z, lmax = 2)
```

Here:
- `diag` checks whether the basis is behaving well numerically
- `grid` is the separate quadrature grid used for integrals
- `ops` bundles the basic radial matrices for atomic-style work

That basis/quadrature separation is one of the central ideas in the package.

## Hydrogen is the first real scientific check

The most useful first physics example is hydrogen.

The repository includes:

- `examples/04_hydrogen_ground_state.jl`
- `examples/15_atomic_hydrogen_ylm.jl`

The first of these builds a recommended radial basis, forms the one-electron hydrogen Hamiltonian for `l = 0`, and checks the lowest eigenvalue against the exact nonrelativistic ground-state energy `-0.5 Ha`.

The second adds explicit angular channels `(l,m)` up to a chosen `lmax`, assembles the corresponding block-diagonal one-electron atomic Hamiltonian, and reports the low-lying hydrogen energies.

If you are trying to decide whether the package is working for you, this is the first example that really answers the question.

## What about the primitive layer and contraction line?

Every basis function in the package is built from a lower-level Gaussian-type primitive layer. For most users, this matters only later. But it becomes important if you want to:

- inspect the shared Gaussian building blocks behind a basis
- contract primitive-space matrices up to basis-space matrices
- study partitions and local contraction in 1D
- explore the newer nested/hierarchy direction in the repository

That line is advanced rather than beginner-facing. It is worth learning, but only after the radial path is comfortable.

## Documentation map

The main documentation surface is:

- [`docs/first_radial_workflow.md`](docs/first_radial_workflow.md)  
  The best first read after this README.
- [`docs/recommended_atomic_setup.md`](docs/recommended_atomic_setup.md)  
  Practical parameter choices for atom-centered radial work.
- [`docs/atomic_ylm_layer.md`](docs/atomic_ylm_layer.md)
  The first explicit angular `(l,m)` layer built on top of the radial substrate.
- [`docs/atomic_ida_layer.md`](docs/atomic_ida_layer.md)
  The first static interacting He / IDA-style ingredients built on top of the radial-plus-angular layer.
- [`docs/atomic_ida_uhf.md`](docs/atomic_ida_uhf.md)
  The current minimal UHF kernel built on top of the present atomic IDA model.
- [`docs/mapped_ordinary_basis.md`](docs/mapped_ordinary_basis.md)
  The first public globally mapped full-line basis route for the ordinary Cartesian hydrogen branch.
- [`docs/example_guide.md`](docs/example_guide.md)  
  The best order in which to read and run the examples.
- [`docs/architecture.md`](docs/architecture.md)  
  A compact map of the package structure: ordinary gausslets, radial work, primitive layers, hierarchy, and the current corrected nested direction.
- [`docs/terminology.md`](docs/terminology.md)  
  Plain-language explanations of the few package-specific terms.
- [`docs/intermediate_primitive_layer.md`](docs/intermediate_primitive_layer.md)  
  The first advanced note about primitive layers, contraction, and why they matter.

There are also narrower supporting notes about the current research direction in 1D contraction and hierarchy. Those are useful, but they are **not** the best starting point for a new user.

## Example order

A good reading/running order is:

### Start here
- `01_first_gausslet.jl`
- `02_radial_basis.jl`
- `03_radial_operators.jl`
- `04_hydrogen_ground_state.jl`
- `15_atomic_hydrogen_ylm.jl`
- `16_atomic_ida_ingredients.jl`

### Then move to the advanced primitive-layer line
- `05_primitive_sets.jl`
- `06_basis_contraction.jl`
- `07_position_contraction.jl`
- `08_basis_representation.jl`

### Then, if you want the ordinary Cartesian hydrogen branch
- `23_cartesian_hydrogen_coulomb_expansion.jl`
- `24_mapped_cartesian_hydrogen.jl`

### Then, if you want the current corrected nested/contraction direction
- `09_basis_partition.jl`
- `10_hierarchical_partition.jl`
- `13_global_leaf_contraction.jl`

### Prototype side branch
- `11_leaf_pgdg.jl`
- `12_leaf_pgdg_augmentation.jl`

The first four examples are the public-facing heart of the package. The later ones are there to explain and test the more experimental contraction and hierarchy ideas. In particular, `13_global_leaf_contraction.jl` is now the clearest example of the corrected current direction, while `11` and `12` should be read as useful prototypes rather than as the main conceptual path.

## Current scope and limits

What is already useful today:

- ordinary 1D gausslet objects and explicit Gaussian constructions
- the first public globally mapped full-line ordinary basis route
- radial gausslet bases and radial one-body operators
- the first explicit one-electron `(l,m)` atomic layer
- the first static He / IDA-style interacting atomic ingredients
- the first direct / exchange / Fock / minimal-UHF atomic line in the current IDA model
- explicit quadrature and diagnostics
- primitive-layer matrix construction and contraction
- the first 1D partition/hierarchy/contraction experiments

What is not yet here:

- exact non-diagonal electron-electron operators
- the fuller interacting spherical-angular atomic layer beyond the present one-electron `(l,m)` path
- a broad general HF workflow beyond the present minimal atomic IDA UHF kernel
- 2D or 3D nested workflows
- named Gaussian chemistry basis libraries
- larger solver layers such as DMRG or related workflows
- Python and Fortran interoperability layers

That is deliberate. The package is still settling its scientific structure before it grows into those directions.
