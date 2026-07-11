# Residual Gaussian Domain Module

This page is the canonical contract for the implemented internal
`CartesianResidualGaussians` domain module. It owns the ordinary owner-local
residual basis, exact augmented one-body and moment transforms, final-basis
moment-matched Gaussian descriptors, and residual-containing MWG/IDA blocks.

It does not own supplement loading, terminal construction, Gaussian raw
integrals, artifacts, facade parsing, or solver behavior.

## Lifecycle

| ID | Lifecycle | Current boundary |
| --- | --- | --- |
| `HP-RG-FILE-01` | Implemented | Internal module and three owner files |
| `HP-RG-OBJ-01` | Implemented | `CartesianResidualGaussianBasis` |
| `HP-RG-FN-01` | Implemented | Owner-local basis construction and one final merge |
| `HP-RG-FN-02` | Implemented | Exact augmented operator transformation |
| `HP-RG-FN-03` | Implemented | Final-residual moment-matched descriptors |
| `HP-RG-FN-04` | Implemented | Residual-containing MWG/IDA interaction |
| `HP-RG-WIRE-01` | Implemented | Live terminal, facade, and driver composition |
| `HP-RG-TEST-01` | Validation completed; active maintenance | Existing H2 endpoint |

The module extraction landed in `93ac83ab5`, `7c0651be5`, and `7e45d5da1`.
Manager-log Passes 065-073 record the migration and closeout. Numerical
robustness and cutoff precedence are canonical in
[Residual Gaussian orthogonality and cutoff policy](residual_gaussian_orthogonality_robustness.md).

## Source Ownership

The internal owner is:

```text
src/cartesian_residual_gaussians/CartesianResidualGaussians.jl
src/cartesian_residual_gaussians/residual_basis.jl
src/cartesian_residual_gaussians/augmented_operators.jl
src/cartesian_residual_gaussians/mwg_interaction.jl
```

`src/GaussletBases.jl` includes the module but exports no RG-domain API.

Current composition remains in:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/cartesian_base_hamiltonian.jl
src/cartesian_protected_ladder_bundle.jl
bin/cartesian_ham_builder.jl
```

The terminal file constructs exact raw `GG/GA/AA` inputs, validates dimensions
and Coulomb-expansion parity, delegates domain mathematics, assembles the
existing `CartesianIDAHamiltonian`, and owns the compact artifact hook. The
base facade and canonical driver remain callers. The protected-ladder owner
uses the same RG builder and exact transform for separate opt-in compositions.

These are live compatibility and composition surfaces, not duplicate owners
of residual selection, exact transforms, or MWG mathematics.

## Basis Object

`CartesianResidualGaussianBasis` is a numerical construction object with the
following live fields:

```text
base_dimension
candidate_count
residual_dimension

candidate_labels
candidate_owner_indices
candidate_centers
compact_source_candidate_indices

residual_source_owner_indices
residual_occupations
owner_retained_counts
residual_labels

T_G
T_A
injected_G
injected_A

occupation_cutoff
residual_injection_cutoff
tau_neg_abs
tau_neg_rel
tau_merge_abs
tau_merge_rel

selection_rule
orientation
sign_rule
```

For `nG` orthonormal terminal functions, `nA` supplement functions, and `nR`
retained residuals:

```text
T_G :: nG x nR
T_A :: nA x nR
R   = G*T_G + A*T_A.
```

The ordinary path has `injected_G = injected_A = nothing` and
`residual_injection_cutoff = 0`. Compact-source and injection fields support
separately documented opt-in consumers; their policy does not belong to this
core contract.

The object is not a route status, report payload, artifact record, or public
input. Candidate metadata records the supplied supplement ordering and exact
physical owner assignment. Residual metadata records selected owner-local
source modes and the policy used to construct them.

## Ordinary Basis Construction

Inputs are:

```text
G       orthonormal terminal basis
A       explicit Cartesian Gaussian supplement
X       = <G|A>, size nG x nA
S_AA    = <A|A>, size nA x nA
owner   one physical nucleus for every supplement column
```

Every candidate center must exactly match one and only one supplied nucleus.
Selection is performed separately for each physical owner. For owner `a`, the
residual metric is

```text
M_a = S_AA[a,a] - X[:,a]' * X[:,a].
```

The metric is explicitly symmetrized. Its negative-eigenvalue tolerance is

```text
tau = max(tau_neg_abs,
          tau_neg_rel * max(maximum(eigenvalues), 1)).
```

An eigenvalue below `-tau` is a construction failure. Retention uses the
strict physical rule

```text
lambda > residual_occupation_cutoff.
```

No eigenvalue is floored or clamped upward. If an owner block retains all
directions, the implementation uses the symmetric inverse square root of its
metric. If rank is physically lost at the cutoff, it uses the retained natural
modes ordered by decreasing occupation. If no owner retains a direction, the
construction fails.

The retained owner blocks are concatenated once. Before the final merge,

```text
T_G0 = -X*T_A0,
```

so the residuals are orthogonal to `G` in the represented metric. The single
inter-owner merge forms

```text
S_merge = <R0|R0>
```

and applies one symmetric inverse square root. The merge uses the same
absolute/relative negative-metric check and additionally requires

```text
minimum(eigenvalues(S_merge)) > tau_merge.
```

A zero, negative, or numerically near-singular final merge is a hard failure.
There is no merge-eigenvalue flooring. Finally, signs are deterministic: the
largest-magnitude entry in each `T_A` column is positive.

The ordinary result must satisfy

```text
norm(T_G + X*T_A, Inf) <= orthogonality_atol
<R|R> approximately I under the current identity policy.
```

Current production defaults and the exact scale-aware identity check are in
the [orthogonality and cutoff contract](residual_gaussian_orthogonality_robustness.md).
The production cutoff is `1e-6`. The separate numerical-complete opt-in uses
an explicit `1e-10`; it does not alter this default.

## Exact Augmented Operators

`transform_augmented_operator(O_GG, O_GA, O_AA, residual)` owns the exact
ordinary `[G,A] -> [G,R]` transformation. For the ordinary non-injected path:

```text
O_GR = O_GG*T_G + O_GA*T_A

O_RR = T_G'*O_GG*T_G
     + T_G'*O_GA*T_A
     + T_A'*O_GA'*T_G
     + T_A'*O_AA*T_A.
```

The returned matrix is assembled in native `[G,R]` order and explicitly
symmetrized. This transform is used for:

- kinetic energy;
- uncharged nuclear attraction by center;
- `x`, `y`, and `z` moments;
- `x2`, `y2`, and `z2` moments.

The ordinary `GG` block remains the supplied exact terminal block. Physical
nuclear charges are applied by the Hamiltonian composition owner after the
uncharged matrices are transformed. Protected fixed-sector transforms and
packet/reference transforms are separate contracts in the same source file;
they do not change this ordinary formula.

## Moment-Matched Residual Descriptors

`moment_matched_gaussians(operators, residual)` reads diagonal first and
second moments of the final merged residual functions. For residual `r` and
axis `u`:

```text
center[r,u]   = <r|u|r>
variance[r,u] = <r|u2|r> - center[r,u]^2
width[r,u]    = sqrt(2*variance[r,u]).
```

Centers and variances must be finite, and every variance and width must be
strictly positive. Descriptors are computed after the final merge because MWG
is not invariant under arbitrary residual rotations.

These matched Gaussians describe only residual-containing two-index density
interaction blocks. They are not exact residual-GTO four-index integrals and
are not replacement basis functions.

## Residual MWG/IDA Interaction

`assemble_residual_ida_interaction(...)` combines:

```text
V_GG    existing terminal IDA/MWG interaction, unchanged
V_GM    terminal-to-matched-residual density interaction
V_MM    matched-residual density interaction.
```

One producer-owned `CoulombGaussianExpansion` must match the PGDG exponent
sequence on every axis. Residual pair terms use that exact expansion.

For a direct terminal block, `V_GM` copies support values directly. For a
support-local compact block with coefficient matrix `C`, support weights `w_s`,
and final density weights

```text
w_f = C' * w_s,
```

the density action uses

```text
C_density = C .* w_s ./ w_f.
```

Every final weight must be finite and greater than `1e-12`. This weight-aware
normalization is required; plain wavefunction projection is not an equivalent
density transform. `V_MM` is assembled symmetrically from the three analytic
matched-Gaussian axis factors.

The resulting two-index interaction is an IDA/MWG approximation for density
interactions. Exact augmented one-body transformation and approximate MWG
interaction must not be conflated, and the interaction is never rotated as
`C' V C`.

## Validation

The tracked endpoint is:

```text
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
```

It currently checks:

- owner-local candidate and retained-owner metadata;
- the production `occupation_cutoff == 1e-6`;
- base/residual/final dimensions `487 / 18 / 505`;
- `G-R` orthogonality and residual identity;
- exact kinetic, unit-nuclear, coordinate, and second-moment transforms;
- finite symmetric final interaction;
- an independent weight-aware `V_GM` reconstruction;
- compact H2 self-Coulomb `0.4574161883692301` within `1e-10`;
- facade, artifact readback, and provenance compatibility.

`test/misc/runtests.jl` contains the compact numerical-complete metric
contract at explicit cutoff `1e-10`. That test protects the separate opt-in
policy and must not be cited as the production default.

Validation IDs authorize maintenance of these existing contracts. They do not
grant new source, fixture, endpoint, or Cr2 workflow authority.

## Related Contracts And Non-Goals

- Numerical tolerances and cutoff precedence:
  [orthogonality and cutoff policy](residual_gaussian_orthogonality_robustness.md).
- Explicit numerical-complete `1e-10` composition:
  [numerical-complete residual basis](numerical_complete_residual_basis.md).
- Default-off direct injection:
  [residual injection compatibility](residual_gaussian_injection_hybrid.md).
- Protected replacement and additive-reference consumers:
  [protected-localized basis](protected_localized_basis.md) and
  [protected additive reference](protected_additive_reference_correction.md).
- Raw exact Gaussian matrices:
  [nuclear raw blocks](cartesian_gaussian_raw_blocks_nuclear.md) and
  [non-nuclear raw blocks](cartesian_gaussian_raw_blocks_non_nuclear.md).
- Persistence:
  [ordinary artifact manifest](cartesian_hamiltonian_artifact_manifest.md).

This core owner does not authorize global residual selection, eigenvalue
flooring, automatic width/zeta filtering, residual localization, Gaussian-array
enrichment, protected/injection policy, EGOI, screened Hartree, artifact schema,
public exports, solver workflow, or Cr2-specific behavior.
