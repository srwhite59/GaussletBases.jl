# Cartesian Gaussian Raw Blocks - Nuclear Slice

This page owns the implemented neutral kernel for exact uncharged Cartesian
Gaussian nuclear `G-A` and `A-A` raw blocks. It is internal numerical
infrastructure, not a route, cache, artifact, or public API.

## Lifecycle

| ID | Lifecycle | Current boundary |
| --- | --- | --- |
| `HP-CGRB-FILE-01` | Implemented | Neutral module and nuclear owner file |
| `HP-CGRB-FN-01` | Implemented | Exact uncharged by-center `G-A`/`A-A` blocks |
| `HP-CGRB-FN-02` | Implemented | Term-first one-dimensional axis-family reuse |
| `HP-CGAI-FN-01` | Unused optional authority; superseded as an endpoint | Proposed nonallocating axis-table helper never landed |
| `HP-CGRB-WIRE-01` | Implemented | Residual-Gaussian and Qiu-White callers use the neutral owner |
| `HP-CGRB-TEST-01` | Validation completed; tracked coverage is indirect | H2/core tests plus accepted QW/Be2/Cr2 parity evidence |

Implementation commits are `5da4c8a6e` for extraction, `47d9b2a3e` for
streamed `A-A`, `82b3f697f` for axis-family reuse, and `7c9afa8bd` for the
stable high-exponent formula. Manager-log Passes 078-082B record acceptance.

## Ownership And Callers

The internal module is loaded from:

```text
src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl
src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl
```

`src/GaussletBases.jl` includes the module but exports no raw-block API. The
module also includes non-nuclear and mixed-Hartree files; those have separate
contracts. In particular, `mixed_hartree_blocks.jl` is not part of this slice.

The implemented entry point is:

```julia
CartesianGaussianRawBlocks.gaussian_nuclear_raw_blocks_by_center(
    proxy, supplement, expansion, atom_locations)
```

Its direct live callers are:

- `_r3a_qw_blocks(...)` in
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`;
- `_qwrg_diatomic_cartesian_shell_blocks_3d(...)` in
  `src/ordinary_qw_raw_blocks.jl`.

The latter remains part of the live Qiu-White operator path in
`src/ordinary_qw_operator_assembly.jl`. Residual selection, exact augmented
transforms, and Qiu-White route semantics remain with those consumers.

## Returned Blocks

Let `nG = proxy.ncart`, `nA = length(supplement.orbitals)`, and `nC` be the
number of nuclear centers. The return value is exactly:

```text
ga :: Vector{Matrix{Float64}}   # nC matrices, each nG x nA
aa :: Vector{Matrix{Float64}}   # nC matrices, each nA x nA
```

The fields are returned as `(; ga, aa)`. Center order is exactly
`atom_locations` order. Supplement columns/rows follow
`supplement.orbitals` order. `A-A` matrices are explicitly mirrored and
symmetrized; `G-A` is rectangular.

Each center matrix is the unit attractive operator

```text
U_A = -1 / |r - R_A|.
```

No physical nuclear charge is accepted or applied. Hamiltonian/one-body
consumers form `sum_A Z_A U_A`; the current Qiu-White caller does this
explicitly after raw-block construction.

## Numerical Convention

The caller supplies one resolved `CoulombGaussianExpansion`. The kernel uses
that exact coefficient/exponent sequence for every center and both block
families; it has no local compact/high selector or fallback.

For left, right, and Coulomb-factor exponents `alpha_l`, `alpha_r`, and
`alpha_f`, the one-dimensional factor integral uses

```text
gamma = alpha_l + alpha_r + alpha_f

Q = (alpha_l*alpha_r*(x_l-x_r)^2
   + alpha_l*alpha_f*(x_l-x_f)^2
   + alpha_r*alpha_f*(x_r-x_f)^2) / gamma.
```

The polynomial moment sum is evaluated about the weighted center, multiplied
by the existing left/right axis prefactors and `exp(-Q)`. The pairwise form is
algebraically equivalent to the old weighted-variance expression but avoids
high-exponent cancellation. Polynomial powers, moments, primitive
normalization, Coulomb coefficients, and the negative nuclear sign are
unchanged.

Three-dimensional contracted elements retain the coupled primitive form:

```text
sum_pq c_p c_q I_x[p,q] I_y[p,q] I_z[p,q].
```

The x/y/z tables must not be contracted independently into three unrelated
shell scalars.

## Reuse And Ordering

The implementation builds one function-local supplement axis-family inventory.
An axis family is keyed by primitive exponents, axis center, Cartesian power,
and axis prefactors. Orbital-to-family integer maps preserve the original
flattened supplement order.

For each Coulomb term, the kernel:

1. fills unique `G-A` tables by axis family and unique nuclear coordinate;
2. fills canonical `A-A` family-pair tables with reversal flags;
3. reuses those tables across all 3D orbitals and upper-triangle orbital pairs;
4. accumulates the coupled primitive products immediately.

Repeated x, y, or z coordinates among different nuclei reuse the same axis
table while output center order remains unchanged. All inventories and scratch
matrices are construction-local; no persistent cache or provider object is
created.

The accepted Cr2 q4 audit reduced raw nuclear allocation from roughly
`44.5 GiB` to tens of MiB while preserving `G-A`/`A-A` values at about
`1e-14`. This is evidence for the implemented organization, not a runtime or
Cr2-specific contract.

## Axis-Helper Disposition

`HP-CGAI-FN-01` proposed an in-place
`_cartesian_gaussian_axis_integral_table!(...)` or specialized nonallocating
factor helper in `src/cartesian_gaussian_axis_integrals.jl`. That symbol and
consumer never landed. `HP-CGRB-FN-02` succeeded with the private
`_factor_axis_integral(...)` and `_fill_axis_factor_table!(...)` owned directly
by `nuclear_blocks.jl`; therefore `HP-CGAI-FN-01` remains unused optional
authority and is superseded as a performance endpoint.

The later allocating `_cartesian_gaussian_axis_integral(...)` and
`_cartesian_gaussian_axis_integral_table(...)` helpers are committed and live
for Qiu-White and non-nuclear callers. They were introduced under later neutral
axis-kernel consolidation and do not make the unimplemented in-place helper
source-backed.

## Validation And Failure Boundary

Tracked coverage is:

```text
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
test/core/runtests.jl
```

The H2 test exercises the neutral blocks through exact augmented operators.
The core BigFloat oracle checks the stable factor-axis formula through compact
and tight high-exponent cases. Accepted ignored Qiu-White, Be2, and Cr2 parity
evidence remains in manager-log Passes 078-082B; no dedicated committed
`cartesian_gaussian_raw_blocks_nuclear_runtests.jl` file exists.

This internal kernel assumes validated proxy, supplement, expansion, and center
objects. Invalid axes, inconsistent primitive shapes, and linear-algebra
dimension errors throw or propagate; there is no status-bearing fallback.

The nuclear owner does not own or authorize overlap/kinetic/moment blocks,
terminal projection, residual selection or transforms, Qiu-White route objects,
pair/MWG interactions, final-basis `G-G` products, caches, metadata, artifacts,
drivers, solvers, or element/Cr2 policy.
