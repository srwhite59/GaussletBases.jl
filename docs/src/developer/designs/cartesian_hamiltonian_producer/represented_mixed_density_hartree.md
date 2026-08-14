# Represented Mixed-Density Hartree Producer

Status: approved internal implementation and bounded validation under
`HP-REP-MIXDENS-HARTREE-FN-01` and
`HP-REP-MIXDENS-HARTREE-TEST-01`. No public API, artifact, solver, or Cr2
endpoint is authorized.

## Purpose

Construct the complete direct Hartree field of a represented multicenter
molecular determinant whose occupied orbitals contain both terminal-gausslet
and Gaussian-supplement content. The facility owns three ordered operations:

1. reconstruct the represented spin densities in the mixed source basis;
2. construct their direct Coulomb potential with one explicit Gaussian Coulomb
   expansion;
3. evaluate complete terminal/supplement `GG`, `GA`, and `AA` one-body field
   blocks and transform them into the native final basis.

This is a reusable internal Hartree producer. It is not screened-Hartree
algebra, an occupied-only image, a general four-index engine, or a molecular
workflow.

## Existing Boundaries

The implementation must reuse, not duplicate:

- `CartesianResidualGaussianBasis.T_G/T_A` for the final-to-mixed
  reconstruction;
- terminal Gaussian-sum and neutral mixed-Hartree target kernels in
  `cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl`;
- polynomial-Gaussian product, convolution, and pair-factor algebra in
  `GaussianAnalyticIntegrals.jl`, with the dense Gaussian Coulomb reference
  serving only as an extraction donor and bounded oracle;
- the explicit `CoulombGaussianExpansion` object supplied by the caller;
- `transform_augmented_operator` for the exact one-body transformation from
  raw `GG/GA/AA` blocks to native `[G,R]` order;
- `gaussian_coulomb_pair_matrix` only as a bounded test oracle.

The existing same-center atomic wrappers remain the fast path governed by
[Reference Hartree numerics](reference_hartree_numerics.md). This lane must not
change their accepted outputs or turn them into multicenter workflow objects.

## Layer 1: Represented Density

Let `G` be the orthonormal terminal basis, `A` the ordered supplement, and

```text
X    = <G|A>
S_AA = <A|A>.
```

The native final basis is `B = [G,R]`, with each residual column represented by
`T_G` and `T_A`. For each spin and final occupied coefficient block

```text
C_final = [C_G_final; C_R],
```

reconstruct the mixed-source coefficients exactly as

```text
C_G = C_G_final + T_G*C_R
C_A = T_A*C_R.
```

The mixed metric and spin densities are

```text
S_mixed = [I   X;
           X'  S_AA]

P_sigma = C_sigma * Diagonal(n_sigma) * C_sigma'
P        = P_alpha + P_beta.
```

Hartree construction uses the spin sum, but alpha and beta coefficient blocks,
occupations, traces, and reconstruction errors remain separately visible in
diagnostics. Low-rank occupied coefficients are the primary representation;
the source must not materialize a dense terminal density merely to carry the
object between layers.

The private `RepresentedMixedDensity` object owns:

- the ordered `C_G` and `C_A` blocks by spin and their occupations;
- terminal, supplement, residual-transform, and final-order fingerprints;
- the explicit `X` and `S_AA` identities used for validation;
- spin and total mixed-metric charge diagnostics.

For block densities `P_GG`, `P_GA`, and `P_AA`, total represented charge is

```text
N = Tr(P_GG) + 2*Tr(X' * P_GA) + Tr(S_AA * P_AA).
```

Construction must fail on final-order, dimension, fingerprint, metric,
occupation, or reconstruction mismatch. It must verify that mixed-metric
orbital Grams and spin traces reproduce the supplied orthonormal final-basis
state to the declared numerical tolerance. Labels alone are not identity.

## Layer 2: Direct Coulomb Potential

`DirectRepresentedHartreePotential` is a private concrete route constructed
from one `RepresentedMixedDensity` and one explicit
`CoulombGaussianExpansion`. It forms the complete pair density, including

```text
G-G source products
G-A and A-G source products
A-A source products.
```

Dropping mixed source products, projecting back to an external AO span, or
using only occupied-space field images defines a different object and is
forbidden.

The direct route may stream, combine algebraically identical Gaussian product
terms, cache one-dimensional factors, and retain a factorized potential. Any
compression must be algebraically exact for the represented density and the
supplied expansion. It may not fit, truncate, screen, or silently substitute a
different density.

### Contracted Pair-Product Action

The production source operation missing at Pass 459 is a nonmaterializing
contracted PGDG/Gaussian pair-product action. For each occupied spin orbital
and each Coulomb expansion exponent, it must contract

```text
rho_GG = sum_s sum_k n[s,k] sum_ij C_G[i,k,s] C_G[j,k,s] G_i G_j
rho_GA = sum_s sum_k n[s,k] sum_ia C_G[i,k,s] C_A[a,k,s]
         (G_i A_a + A_a G_i)
rho_AA = sum_s sum_k n[s,k] sum_ab C_A[a,k,s] C_A[b,k,s] A_a A_b.
```

All diagonal and off-diagonal products are required. The real-orbital
`GA/AG` multiplicity must be represented exactly once with its factor of two;
symmetry must not drop transition products or count them twice.

The action consumes block-local terminal parent states and contraction
coefficients together with explicit supplement primitives. For the separable
kernel term

```text
exp(-gamma * |r-r'|^2)
  = product_axis exp(-gamma * (x_axis-x'_axis)^2),
```

it constructs or applies the exact one-dimensional polynomial-Gaussian
source-pair convolution, contracts it immediately with the low-rank occupied
coefficients, and accumulates the target potential factors. It may tile source
block pairs and cache repeated one-dimensional state/family factors.

It must not construct or retain:

- an ordered or compact source-pair Coulomb matrix;
- a `pair_count x pair_count` kernel;
- all terminal/supplement primitive-pair terms at once;
- dense `P_GG` merely to drive the contraction.

The required dense `GG/GA/AA` output matrices are not forbidden by this rule.
Workspace and caches must scale with occupied rank, source/target tiles, and
one-dimensional state/family inventories, not with the square of the global
source-pair count.

Existing parent/PGDG IDA factors contain charge-density approximations and
cannot replace `G_i G_j` transition products. Calling the dense
`gaussian_coulomb_pair_matrix` from production is also forbidden.

### Shared Analytic Kernel

The implementation may extract one narrow allocation-conscious
polynomial-Gaussian pair/convolution primitive into
`GaussianAnalyticIntegrals.jl`. It may then replace duplicated private algebra
in `gaussian_coulomb_reference.jl` with that shared primitive, provided the
dense oracle remains numerically unchanged and remains an oracle only.

This extraction does not authorize a general ERI API, shell-quartet engine,
provider hierarchy, or new root export. Three-dimensional represented-density
orchestration remains in the approved raw/reference-density owners.

“Direct” means exact contraction of the represented density under the recorded
Gaussian Coulomb expansion. It is not a claim that a finite expansion is exact
continuum `1/r`. Diagnostics must retain the expansion policy, parameters,
term count, and fingerprint separately from density and basis fingerprints.

## Direct And Fitted Routes

The first source pass is direct-only. It must not add a fitted potential type,
fit-only scaffolding, or a symbolic `source=:direct|:fit` switch.

Any later auxiliary fit requires a separate docs-only amendment and a distinct
concrete type, such as `FittedRepresentedHartreePotential`. That type must
retain the source `RepresentedMixedDensity` and direct-potential fingerprints,
fit controls, rank, and certified errors. It may never dispatch through the
direct/oracle method or satisfy an oracle parameter by shared fields alone.

A fitted route can be accepted only after the direct route passes. Its complete
field must agree with the direct route within the governing matrix and energy
gates, or carry a rigorous certified bound over the complete final space.
Occupied columns, selected entries, and random probes are diagnostics, not a
full-space certificate.

## Layer 3: Complete Mixed-Basis Field

For a represented density `P`, evaluate the complete raw one-body field

```text
J_GG = <G | v[P] | G>
J_GA = <G | v[P] | A>
J_AA = <A | v[P] | A>.
```

The private `RepresentedHartreeField` result owns the three materialized raw
blocks, the concrete potential route, ordering and kernel fingerprints,
Hermiticity diagnostics, and the explicit no-half direct scalar. The first
accepted implementation must materialize all three blocks. A factorized apply
may coexist later only if it reproduces those blocks; it cannot replace the
complete-matrix acceptance gate.

Transform the raw triple with the existing one-body operator transform:

```text
J_B = transform_augmented_operator(J_GG, J_GA, J_AA, residual),
```

where `B = [G,R]` is native final order. This is an exact one-body field
transformation. It does not authorize `C' V C`, a transformed interaction, or
a generalized final-basis overlap workflow.

The producer reports separately

```text
E0 = Tr(P_B * J_B) = (rho|rho)
EH = 0.5*E0.
```

No hidden factor of one half is permitted. For a trace-preserving variation
`delta P`, the independent functional check is

```text
d[0.5*(rho|rho)] = Tr(delta P * J[P]).
```

The field and scalar must come from the same represented density and Coulomb
route. A density/potential/field fingerprint mismatch is a hard failure.

## Initial Internal Surface

The first implementation may add only these private names inside
`CartesianReferenceDensity`:

- `RepresentedMixedDensity`;
- `ContractedMixedDensityPairAction`;
- `DirectRepresentedHartreePotential`;
- `RepresentedHartreeField`;
- narrow constructors/evaluators implementing the three layers above.

Exact private spelling of the constructors may follow local Julia style, but no
name is root-exported. Do not add an abstract provider hierarchy, persistent
artifact shape, metadata carrier, public keyword, or alternate correction
object.

Approved source ownership is limited to:

- private `src/cartesian_reference_density/represented_molecular_hartree.jl`;
- `src/cartesian_reference_density/CartesianReferenceDensity.jl` for one
  include;
- `src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl` for narrow reuse
  or generalization of neutral pair-density/potential target kernels;
- `src/GaussianAnalyticIntegrals.jl` for the narrow shared one-dimensional
  pair/convolution primitive;
- `src/gaussian_coulomb_reference.jl` only to remove duplicated private
  polynomial-kernel algebra and preserve oracle parity.

Existing residual transforms and screened-Hartree assembly are consumers, not
edit surfaces in this pass.

## Bounded Repository Validation

The committed test owner is
`test/nested/cartesian_represented_molecular_hartree_runtests.jl`. Use a small
multicenter H2-style or synthetic terminal-plus-supplement fixture, small enough
for an independent dense primitive-pair oracle built with the same explicit
Coulomb expansion.

The test must cover:

- alpha, beta, and total mixed-metric charge and final/mixed Gram recovery;
- nonzero `G-A` source-density terms, so omitting mixed products fails visibly;
- isolated `GG`, `GA/AG`, and `AA` source sectors plus their combined result,
  including translated centers and off-diagonal signs/multiplicities;
- complete `GG`, `GA`, and `AA` agreement with an independently contracted
  small four-index oracle;
- native `[G,R]` transformed-field parity through
  `transform_augmented_operator`;
- Hermiticity, finite values, no-half `E0`, and `0.5*E0` derivative identity;
- predeclared trace-preserving terminal, supplement, and mixed-sector density
  variations rather than occupied-only probes;
- rejection of malformed ordering, transforms, metrics, occupations, and
  density/potential fingerprints;
- exact unchanged output from the existing same-center atomic path.

Add a bounded resource-scaling gate over two fixture sizes. It must report
occupied rank, source/target dimensions, source block-pair count, retained
one-dimensional cache entries, output bytes, peak additional workspace, and
the largest retained workspace shapes. The test and source diff must establish
structurally that no global pair-index matrix exists; do not add a persistent
status flag for that claim. The larger fixture must demonstrate that additional
workspace follows the tiled/action representation rather than global
source-pair squared or `N^4` storage. Do not turn wall time into a fragile
committed threshold.

These are numerical-object tests, not a solver or Cr2 fixture. The test may use
the existing dense `N^4` Gaussian Coulomb helper only at its bounded oracle
size; production source must not call it as a large-system backend.

## Separate Cr2 Consumer Acceptance

After repository tests and manager review pass, CR2 may run a separate ignored
or external consumer acceptance against the frozen represented determinant.
That pass must retain the precommitted REQ-084 gates:

- charge error at most `1e-10 e`;
- pre-affine direct-energy error at most `1e-5 Ha`;
- complete and terminal spectral-norm field errors at most `1e-5 Ha`;
- complete external-row/column operator-norm error at most `1e-5 Ha`;
- Hermiticity infinity error at most `1e-10 Ha`;
- independent derivative error at most `1e-8 Ha`.

The complete-field errors require an independent complete represented-density
oracle or a rigorous certified full-space bound. An external-AO projected `J`,
occupied action, selected rows, or random vectors cannot satisfy that gate.
If the direct producer cannot scale and no rigorous full-space certificate is
available, stop before molecular-full energies. Do not admit a fit as the
oracle or weaken the gates.

Passing this consumer gate validates the reusable producer for that state. It
does not create a Cr2 default, production endpoint, artifact, or solver path.

## Failure Rule And Budget

The preferred implementation budget is at most `350` added source lines and
the hard limit is `450` added source lines across all approved source files,
including module wiring and any shared-kernel extraction. The existing test
limit remains `180` added lines in the one planned test file. No other source
or test file is approved.

If a coherent direct implementation exceeds the hard limit, makes the atomic
path regress, or cannot provide complete raw blocks, make no source commit.
Report the smallest justified reusable extraction or user-approved budget
increase. Do not split partial scaffolding across commits.

## Explicit Non-Goals

This authority does not permit:

- a general four-index ERI engine or production call to the dense `N^4` oracle;
- dense global source-pair matrices, IDA transition-product substitution,
  quadrature, product truncation, or pair screening;
- a fitted molecular potential in the first source pass;
- an occupied-only, external-AO, diagonal, row-gauge, or random-probe substitute;
- changes to screened-Hartree formulas, `q0/P0`, residual selection, MWG, IDA,
  exchange, EGOI, or one-body physics;
- public exports, facade/driver controls, artifacts, sidecars, solvers, SCF/HF,
  or Cr2-specific source behavior;
- molecular refitting, endpoint interpretation, or a production Cr2 claim.
