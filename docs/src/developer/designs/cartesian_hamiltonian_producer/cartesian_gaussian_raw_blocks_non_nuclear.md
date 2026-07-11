# Cartesian Gaussian Raw Blocks - Non-Nuclear Slice

This page owns the implemented neutral `G-A` and `A-A` overlap, kinetic,
coordinate-moment, and second-moment raw blocks. It does not own residual
selection, final-basis transforms, terminal `G-G` products, or Qiu-White
provider semantics.

## Lifecycle

| ID | Lifecycle | Current boundary |
| --- | --- | --- |
| `HP-CGRB-NN-FILE-01` | Implemented | Non-nuclear owner file and module include |
| `HP-CGRB-NN-FN-01` | Implemented | Full and overlap-only neutral raw blocks |
| `HP-CGRB-NN-WIRE-01` | Implemented | Main diatomic Residual-Gaussian/Qiu-White callers use the neutral owner |
| `HP-CGRB-NN-TEST-01` | Validation completed; tracked coverage is indirect | H2 endpoint plus accepted QW/Be2/Cr2 parity evidence |

Implementation commits are `00d052e29` for extraction, `9fa0cc16d` for the
overlap-only path, `806b37e32` for `A-A` family reuse, and `71a89433c` for
`G-A` family reuse. Manager-log Passes 084B-086B record acceptance.

## Ownership And Callers

The owner is:

```text
src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl
```

It is loaded by
`src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl`. Shared
private axis families come from the same internal module, and analytic
overlap/kinetic/position/x2 tables come from
`src/cartesian_gaussian_axis_integrals.jl`.

Implemented entry points are:

```julia
gaussian_non_nuclear_raw_blocks(proxy, supplement, expansion)
gaussian_non_nuclear_overlap_blocks(proxy, supplement)
```

Direct live callers are in:

- `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`;
- `src/ordinary_qw_raw_blocks.jl`.

The full helper supplies exact augmented operators and the main diatomic
Qiu-White wrapper. The overlap-only helper supplies the residual setup mixed
overlap without constructing unused kinetic, moment, or `A-A` blocks.

## Returned Blocks

Let `nG = proxy.ncart` and `nA = length(supplement.orbitals)`. The full helper
returns exactly:

```text
ga.overlap       :: nG x nA
ga.kinetic       :: nG x nA
ga.position.x/y/z:: nG x nA
ga.x2.x/y/z      :: nG x nA

aa.overlap       :: nA x nA
aa.kinetic       :: nA x nA
aa.position.x/y/z:: nA x nA
aa.x2.x/y/z      :: nA x nA
```

The overlap-only helper returns `(; ga = (; overlap))` with one `nG x nA`
matrix. Columns and `A-A` rows follow `supplement.orbitals` order exactly.
There is no sorting or owner regrouping. All `A-A` outputs are explicitly
symmetric; `G-A` outputs are rectangular.

The full helper retains an `expansion` argument for caller-signature parity,
but non-nuclear values do not depend on Coulomb expansion terms.

## Numerical Convention

One-dimensional tables use the shared private axis helper terms:

```text
:overlap
:kinetic
:position
:x2
```

The existing Gaussian polynomial prefactors and analytic-integral sign and
normalization conventions are unchanged. Three-dimensional primitive products
are coupled across axes. Kinetic is the sum of the three terms with one
kinetic axis and two overlap axes; each coordinate or second moment replaces
the corresponding axis overlap table only.

The full helper builds one supplement axis-family inventory keyed by exponents,
axis center, Cartesian power, and axis prefactors. Its `G-A` tables are reused
by family, while `A-A` tables use canonical family-pair IDs and reversal flags.
The overlap-only helper skips every unused operator and all `A-A` work. For an
off-diagonal `A-A` orbital pair, both requested orientations are evaluated and
averaged before the value is mirrored; this preserves the prior Qiu-White
roundoff/symmetry convention.

All tables, maps, and scratch matrices are function-local. The owner returns
raw matrices only, with no status, route, provider, cache, or artifact object.

Accepted Cr2 q4 evidence reduced full non-nuclear raw-block allocation from
about `10.9 GiB` to about `0.86 GiB`, with block parity at roundoff. The
overlap-only residual setup fell from about `10.99 GiB` to about `199 MiB`.
These are compact implementation evidence, not Cr2-specific behavior.

## Consumer Boundaries

`CartesianResidualGaussians` owns residual selection, basis orientation,
augmented transformations, moment-matched descriptors, and MWG interaction.
It consumes these raw matrices but does not own their analytic construction.

The main diatomic Qiu-White path consumes the neutral owner. The following
Qiu-White helpers remain live and must not be deleted under `HP-CGRB-NN-*`:

```text
_qwrg_cartesian_shell_cross_moment_blocks_3d
_qwrg_cartesian_shell_self_moment_blocks_3d
_qwrg_atomic_cartesian_blocks_3d
```

They still serve atomic Qiu-White reference/operator assembly, hybrid
representation sidecars, dense-parent GTO probes, and CPB/provider callers in
`ordinary_qw_operator_assembly.jl`, `cartesian_qw_hybrid_representation.jl`,
`cartesian_gto_probes.jl`, and `CartesianCPBBlockProviders.jl`. Their
`factor_ga`/`factor_aa` outputs are outside this overlap/kinetic/moment
contract. Rewiring or deleting them requires a separate caller/ownership audit.

Terminal final-basis `G-G` product matrices are separately owned by the R3
terminal optimization contract. Mixed-Hartree `GG/GA/AA` reference blocks are
separately owned by the reference-Hartree contract. Neither belongs here.

## Validation And Failure Boundary

Tracked endpoint coverage is:

```text
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
```

It exercises overlap setup and full raw blocks through residual construction
and exact augmented kinetic/position/x2 operators. Accepted ignored
Qiu-White, Be2, and Cr2 parity evidence remains in manager-log Passes 084B-086B;
no dedicated committed
`cartesian_gaussian_raw_blocks_non_nuclear_runtests.jl` file exists.

The internal helpers assume validated proxy/supplement data. Unsupported axis
terms, inconsistent primitive arrays, and linear-algebra dimension errors
throw or propagate; there is no fallback result.

This owner does not authorize nuclear or mixed-Hartree changes, final-basis
`G-G` optimization, residual selection/transforms, Qiu-White semantic changes,
factor blocks, persistent providers/caches, artifacts, drivers, solvers, or
Cr2 workflow.
