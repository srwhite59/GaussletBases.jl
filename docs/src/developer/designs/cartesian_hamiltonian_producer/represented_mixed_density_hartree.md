# Represented Mixed-Density Hartree Producer

Status: the bounded exact implementation at `a77ceed5d` is accepted as a small
oracle, but its production contraction and residual revalidation are not
Cr2-usable. Their replacement remains approved under
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

### Residual Validity Versus State Accuracy

The residual basis and the represented determinant have different numerical
contracts. One call-local tolerance must not control both.

An existing `CartesianResidualGaussianBasis` is revalidated with the current
Residual Gaussian owner rules:

```text
G_R = T_G + X*T_A
S_RR = residual_gaussian_overlap(T_G, T_A, X, S_AA)

norm(G_R, Inf) <= 1e-10

identity_error = maximum(abs, S_RR - I)
identity_scale = maximum(abs, S_RR)
identity_error <= 5e-8 * (1 + max(1, identity_scale)).
```

These are the authoritative cross and scale-aware final-identity checks from
the residual construction contract. Revalidation also requires exact
dimensions, finite transforms, a supported non-injected native `[G,R]`
ordering, and matching candidate identity. It does not rerun selection,
reorthogonalize the basis, change its rank or orientation, or infer validity
from represented charge.

The represented occupied state is then checked independently at `1e-10` for
final/mixed Gram recovery and spin-resolved charge recovery. Field symmetry,
energy closure, and derivatives retain their own later tolerances. No consumer
may pass a looser residual override to weaken charge, state, field, or
derivative gates. Diagnostics report the residual cross error, residual
identity maximum and infinity norms, scale-aware bound, state Gram errors, and
spin charges separately.

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

### Occupied-Contracted Separable Block Action

The component-pair implementation at `a77ceed5d` is exact on bounded fixtures,
but it expands the Cr2 state to `116091` source components, then visits
`6738618186` unique source pairs and `909713455110` high-135 pair terms. That
ordering is a bounded oracle, not the production action.

For each occupied spin orbital and each Coulomb expansion exponent, the
production action must still contract

```text
rho_GG = sum_s sum_k n[s,k] sum_ij C_G[i,k,s] C_G[j,k,s] G_i G_j
rho_GA = sum_s sum_k n[s,k] sum_ia C_G[i,k,s] C_A[a,k,s]
         (G_i A_a + A_a G_i)
rho_AA = sum_s sum_k n[s,k] sum_ab C_A[a,k,s] C_A[b,k,s] A_a A_b.
```

All diagonal and off-diagonal products are required. The real-orbital
`GA/AG` multiplicity must be represented exactly once with its factor of two;
symmetry must not drop transition products or count them twice.

The action first maps each terminal occupied block into its own parent support.
For terminal block `b`, support map `B_b`, and spin occupations `n_sigma`, form
only block-local occupied arrays and density tiles:

```text
D_b_sigma       = B_b * C_G_sigma[b.columns, :]
P_GG[b,c]       = sum_sigma D_b_sigma * Diagonal(n_sigma) * D_c_sigma'
P_GA[b,:]       = sum_sigma D_b_sigma * Diagonal(n_sigma) * C_A_sigma'
P_AA            = sum_sigma C_A_sigma * Diagonal(n_sigma) * C_A_sigma'.
```

`B_b` is the identity on an identity block and the committed terminal
coefficient matrix otherwise. The implementation may retain the occupied
block arrays, but not a global parent density or global parent-pair inventory.
Supplement contractions remain in the contracted orbital/family
representation; primitive coefficients are expanded only inside the relevant
axis contraction.

For the separable kernel term

```text
exp(-gamma * |r-r'|^2)
  = product_axis exp(-gamma * (x_axis-x'_axis)^2),
```

the implementation builds reusable one-dimensional maps keyed by the explicit
coupling exponent and the actual source/target state or Gaussian-family pair.
Each `GG`, `GA/AG`, or `AA` density tile is transformed by a deterministic
sequence of axis-wise matrix/tensor contractions and accumulated directly into
one target terminal-block pair, terminal/supplement block, or supplement block.
It must not implement the six-index operation as scalar source-target quartets.

This changes contraction order only. Every target block receives every source
sector, and terminal support locality is an execution partition rather than a
screening rule. No block, coefficient, primitive, or Coulomb term may be
dropped because it is small or distant.

It must not construct or retain:

- an ordered or compact source-pair Coulomb matrix;
- a `pair_count x pair_count` kernel;
- a flattened global list of parent sites plus supplement primitives for the
  production field path;
- a global parent coefficient or density matrix when block-local occupied
  arrays suffice;
- all terminal/supplement primitive-pair terms at once;
- dense `P_GG` merely to drive the contraction.

The required dense `GG/GA/AA` output matrices are not forbidden by this rule.
Apart from those outputs, storage must scale with occupied rank, the sum of
block-local support arrays, the largest live source/target tensor tile, and
unique one-dimensional state/family tables. It must not scale as the square of
the global parent or primitive-component count. Runtime must be stated as the
sum of the actually scheduled block/tensor contractions under all expansion
terms, not as a nominal global component-pair loop.

Existing parent/PGDG IDA factors contain charge-density approximations and
cannot replace `G_i G_j` transition products. Calling the dense
`gaussian_coulomb_pair_matrix` from production is also forbidden.

The current component-pair loop must be removed from production dispatch. It
may survive only as an explicitly bounded small-fixture oracle with a hard
dimension/component guard if that materially simplifies parity testing;
otherwise the independent Gaussian oracle replaces it. Do not retain two
unbounded production algorithms.

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

The implementation owns only these private names inside
`CartesianReferenceDensity`:

- `RepresentedMixedDensity`;
- `ContractedMixedDensityPairAction`;
- `DirectRepresentedHartreePotential`;
- `RepresentedHartreeField`;
- narrow constructors/evaluators implementing the three layers above.

Exact private spelling of the constructors may follow local Julia style, but no
name is root-exported. Do not add an abstract provider hierarchy, persistent
artifact shape, metadata carrier, public keyword, or alternate correction
object. `ContractedMixedDensityPairAction` should be repurposed as the compact
occupied-block contraction plan rather than supplemented by a parallel plan
type. Its variable-size inventories must be vector-backed.

Approved source ownership is limited to:

- private `src/cartesian_reference_density/represented_molecular_hartree.jl`;
- one private
  `src/cartesian_reference_density/represented_hartree_contractions.jl` for
  the occupied-block/separable contraction kernel;
- `src/cartesian_reference_density/CartesianReferenceDensity.jl` for one
  include;
- `src/cartesian_residual_gaussians/residual_basis.jl` only to expose one
  private reusable residual-validity calculation matching its existing
  construction contract;
- `src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl` for narrow reuse
  or generalization of neutral pair-density/potential target kernels;
- `src/GaussianAnalyticIntegrals.jl` for the narrow shared one-dimensional
  pair/convolution primitive;
- `src/gaussian_coulomb_reference.jl` only to remove duplicated private
  polynomial-kernel algebra and preserve oracle parity.

Residual selection/construction and screened-Hartree assembly remain unchanged.

## Bounded Repository Validation

The committed test owner is
`test/nested/cartesian_represented_molecular_hartree_runtests.jl`. Use a small
multicenter H2-style or synthetic terminal-plus-supplement fixture, small enough
for an independent dense primitive-pair oracle built with the same explicit
Coulomb expansion.

The existing bounded exact checks remain and the replacement test must cover:

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

The fixture must contain at least two terminal blocks, nontrivial terminal
contraction coefficients, multiple occupied alpha/beta columns, and translated
multi-primitive supplement functions. Compare each isolated source sector and
their sum against both the bounded component-pair result, while retained, and
the independent dense Gaussian oracle. This must distinguish occupied-first
block contraction from flattened global pair enumeration.

Add residual-policy checks in which a finite unused residual direction has an
identity error above `1e-10` but within the scale-aware `5e-8` contract and is
accepted while strict state charge/Gram checks still pass. A cross error above
`1e-10`, a scale-aware identity failure, or a state recovery error above
`1e-10` must fail independently. Do not add an `allow_invalid_residual` or
consumer-selected tolerance knob.

Add a bounded resource-scaling gate over two fixture sizes. It must report
occupied rank; terminal block/support and supplement orbital/primitive counts;
evaluated source-sector and target-block contraction counts; unique
one-dimensional table counts and bytes; output bytes; operation estimates by
sector; peak additional workspace; and the largest live tile/intermediate
shapes. Counts must come from the actual deterministic contraction plan. The
test and source diff must establish structurally that no global component-pair
inventory or pair-index matrix exists; do not add a persistent status flag for
that claim. The larger fixture must demonstrate block/occupied scaling rather
than global source-pair-squared or `N^4` storage. Do not turn wall time into a
fragile committed threshold.

These are numerical-object tests, not a solver or Cr2 fixture. The test may use
the existing dense `N^4` Gaussian Coulomb helper only at its bounded oracle
size; production source must not call it as a large-system backend.

## Separate Cr2 Consumer Acceptance

After repository tests and manager review pass, CR2 may run a separate ignored
or external consumer acceptance against the frozen represented determinant.
Before entering complete field construction, that consumer must build the
actual contraction plan and report its exact block/tensor counts, cached
one-dimensional tables, operation estimate, output bytes, peak workspace, and
representative warmed sector timings. The producing host gate is provisionally
`24` hours and `48 GiB`; this is an external workflow stop rule, not a source
default or runtime promise. If the model cannot defend both bounds, stop before
the full field rather than committing preflight-only source.

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

The scaling correction may add at most `650` preferred and `800` hard source
lines across the approved files, including the one new private contraction file
and wiring. It must delete or hard-bound the current flattened component-pair
production loop in the same commit; adding a second unbounded path is not
accepted. The existing test owner may add `120` preferred and `180` hard lines,
for a final file no longer than `360` lines. No other source or test file is
approved.

If a coherent exact contraction exceeds the hard limits, cannot defend the
external `24`-hour/`48`-GiB preflight bounds, makes the atomic or screened path
regress, or cannot provide complete raw blocks, make no source commit. Report
the smallest justified reusable extraction or budget correction. Do not split
partial planning/scaffolding from the working numerical replacement.

## Explicit Non-Goals

This authority does not permit:

- a general four-index ERI engine or production call to the dense `N^4` oracle;
- dense global source-pair matrices, IDA transition-product substitution,
  quadrature, product truncation, or pair screening;
- a fitted molecular potential in the first source pass;
- an occupied-only, external-AO, diagonal, row-gauge, or random-probe substitute;
- changes to screened-Hartree formulas, `q0/P0`, residual selection, MWG, IDA,
  exchange, EGOI, or one-body physics;
- residual reorthogonalization, rank/orientation changes, or a call-local
  validity override;
- public exports, facade/driver controls, artifacts, sidecars, solvers, SCF/HF,
  or Cr2-specific source behavior;
- molecular refitting, endpoint interpretation, or a production Cr2 claim.
