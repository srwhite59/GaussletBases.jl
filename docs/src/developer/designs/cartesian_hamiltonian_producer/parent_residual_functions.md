# Parent Residual Functions And Parent-Backed Gaussian Direct Interaction

Status: implemented internal facility under the four authority IDs below.

Authority IDs:

- `HP-PQS-PRF-FN-01`
- `HP-PQS-PRF-TEST-01`
- `HP-PQS-PARENT-GDIRECT-FN-01`
- `HP-PQS-PARENT-GDIRECT-TEST-01`

This page is the canonical contract for two adjacent facilities:

1. parent residual functions (PRFs), which preserve selected omitted content
   from the mapped Cartesian parent basis while remaining orthogonal to the
   current terminal basis; and
2. an onsite-calibrated Gaussian direct-Coulomb resource for functions with
   explicit coefficients in that same parent basis.

The implemented facility is internal and opt-in. It does not choose shell source
orders, targets, mode counts, cutoffs, or physical states, and it does not
assemble a complete Hamiltonian.

## Purpose And Evidence

Fixed-parent Cr2 measurements showed that lowering one atom-local shell from
source order `q = 7` to `q = 6` and restoring eight selected omitted modes per
owner can recover useful one-body content with fewer functions than the raw
`q = 7` basis. In the measured ladder, the `q = 6` basis plus 16 PRFs was 116
functions smaller than raw `q = 7` and lay about `1.070 mHa` above the
`q = 9` one-body plateau for the imported determinant.

The corresponding direct-interaction audits compared the complete 16 by 6783
PRF-to-terminal block. Relative to full parent IDA, the onsite-calibrated
Gaussian model changed the projected determinant's density-density PRF-G
contribution by `-0.01320 mHa`; current MWG changed it by `-1.07335 mHa`.
Bounded continuum checks also found the Gaussian model materially closer than
MWG for the tested direct integrals.

This evidence supports reusable PRF mechanics and a compact direct resource.
It does not validate transition-density exchange, PRF-to-GTO-residual
interactions, a relaxed HF endpoint, or automatic PRF selection.

## Definitions

Let `P` be the current mapped Cartesian parent product basis, with parent
metric `S_P`. Let `G` be the current orthonormal terminal realization with
implicit parent coefficients `C_G`. A parent-backed function has explicit
coefficients in `P` and uses the same parent axes, ordering, mapping, and
Coulomb expansion as `G`.

A PRF block is a set of parent-backed columns `R` that:

- is confined to one caller-selected terminal shell support;
- is orthogonal to the complete current terminal basis;
- is orthonormal in the parent metric; and
- retains its explicit parent coefficients for exact operators and direct
  interaction evaluation.

The required identities are

```text
C_G' S_P R = 0
R'   S_P R = I.
```

PRFs are not terminal IDA sites, Gaussian supplement orbitals, Residual
Gaussians, protected-original replacements, or fitted density/potential
terms. Their near-zero linear parent/IDA weight is expected and must not be
used as a rejection rule.

## Parent Residual Function Construction

### Consumer-selected targets

The consumer owns the physical selection policy. It supplies:

- the terminal block or exact shell support to inspect;
- one or more selected parent target columns; and
- the final requested column count and orientation.

Targets may come from an external occupied-state analysis, omitted source
modes, a consumer RDM, or another measurement. The repo helper receives
columns, not an RDM-selection policy. It must not choose a shell, `source_q`,
target, occupation threshold, mode count, or top-`k` rule.

For selected targets `T`, first restrict them to the exact chosen support and
project out the complete terminal basis:

```text
X_GT   = C_G' S_P T
T_perp = T - C_G X_GT
M_T    = T_perp' S_P T_perp.
```

The production helper may use the terminal realization block structure and
factorized parent metric; it must not form a dense parent-by-terminal matrix.
It must fail if the selected columns are nonfinite, have material support
leakage, or lose rank after projection. It must not silently discard a column,
floor an eigenvalue, weaken a metric threshold, or replace the consumer's
selection.

When `M_T` is full rank, use one symmetric Lowdin normalization:

```text
R = T_perp * inv(sqrt(Symmetric(M_T))).
```

Column phases use a deterministic largest-parent-coefficient convention.
Rotations, ordering, and selection supplied by the consumer otherwise remain
part of the input convention. No stable artifact identity for individual PRF
columns is established in this lane.

### Internal object

One compact internal vector-backed PRF block type is approved. It may carry
only:

- the source terminal block identity;
- exact support indices/states;
- support-local parent coefficients;
- a compact construction/validation summary.

The object must not store a dense parent-by-terminal transform, an RDM, shell
selection policy, Hamiltonian, solver state, artifact payload, or variable
cloud of per-stage metadata. Multiple owner/shell blocks remain explicit
objects; combining them must validate global PRF orthogonality and terminal
cross overlap.

Required diagnostics are support identity and bounds, parent and PRF counts,
the selected-target metric spectrum, projection loss, maximum support leakage,
`C_G' S_P R` error, `R' S_P R` error, deterministic-phase pivots, parent-charge
sums, and position/spread moments when position operators are available.

### Exact one-body blocks

For each approved parent one-body operator `O_P`, build only the new native
`[G,R]` blocks:

```text
O_GR = C_G' O_P R
O_RR = R'   O_P R.
```

The current `O_GG` block is unchanged. The initial facility covers kinetic,
per-center unit-nuclear attraction, assembled physical `H1` blocks, and
position/second-moment diagnostics through existing factorized parent-axis and
terminal block machinery. It must use the same resolved parent construction
and Coulomb expansion as the terminal basis.

No source-Hamiltonian transform, generalized final overlap, `C' V C`
interaction rotation, or dense full-parent operator is approved.

## Parent-Backed Gaussian Direct Interaction

### Parent-site resource

For each current mapped parent site `p`, use its mapped PGDG center `r_p` and
the onsite value `U_p` from the same parent IDA Coulomb expansion used by the
Hamiltonian construction. `U_p` must be finite and positive. Define a
normalized spherical Gaussian charge exponent

```text
alpha_p = pi * U_p^2 / 2.
```

The direct interaction between parent-site charges is

```text
beta_pq = alpha_p * alpha_q / (alpha_p + alpha_q)

K_pq = erf(sqrt(beta_pq) * R_pq) / R_pq,
R_pq = |r_p - r_q|,
```

with the analytic coincident-center limit

```text
K_pq = 2 * sqrt(beta_pq) / sqrt(pi).
```

Thus `K_pp = U_p`. The implementation must use the analytic limit rather than
divide by zero. It must reject nonfinite/nonpositive onsite data and must not
substitute bare `1/R`, final-row center metadata, MWG widths, or independently
resolved Coulomb parameters.

One compact internal parent-Gaussian direct resource is approved. It may carry
the mapped parent centers, onsite values/exponents, the Hamiltonian-wide
Coulomb-expansion fingerprint, and compact validation facts. It is in-memory
only and is not an artifact contract.

### Function charges and direct blocks

For a normalized parent-backed function with parent coefficients `A[:,a]`,
the model charge is

```text
q_a[p] = abs2(A[p,a]).
```

This convention is valid only for the current orthonormal parent-backed
representation. The implementation must report `sum(q_a)` and reject material
deviation from unity; it must not silently renormalize a malformed column.

For parent-backed functions `a` and `b`, the direct model is

```text
V_ab = q_a' K q_b.
```

The initial composition surface returns only PRF-PRF and PRF-G blocks. Native
terminal-G parent coefficients must be consumed from terminal blocks, with
identity columns and contracted columns handled exactly. The implementation
must tile or stream parent-site kernels and terminal columns; it must not form
a dense parent-by-parent kernel or dense parent-by-terminal coefficient matrix.

Existing G-G `Vee` remains byte-for-byte and semantically unchanged. The new
resource is a direct block evaluator, not authority to splice those blocks
into `CartesianIDAHamiltonian`, MWG, a solver Hamiltonian, or an artifact.

### Oracle boundary

A bounded validation helper may contract the existing full parent IDA direct
kernel against the same parent-backed charges. Full parent IDA is the direct
interaction oracle for this lane; it is not a continuum-exact exchange
oracle. The comparator must remain tiled/bounded and must not become a second
producer interaction policy.

## Exchange Boundary

The Gaussian resource approximates direct charge-charge integrals such as
`(aa|ii)`. It does not construct the transition product required by true
exchange,

```text
(ai|ia).
```

A density-density Hamiltonian may use a direct matrix element in an
exchange-like contraction, but that does not make the result continuum
exchange. Selected Cr2 pairs showed that the direct and transition integrals
can differ substantially. Individual PRF matrix elements are not a physical
error estimate because PRF occupation is small and their values depend on the
chosen PRF orientation.

The next exchange measurement belongs to the consumer: contract the actual
occupied orbitals in `[G,PRF]`, compare their total PRF-induced `Delta J` and
`Delta K` with the G-only projected determinant, and test stability under PRF
rotations and selection changes. No sparse pair correction, transition-density
kernel, exact-exchange helper, or exchange uncertainty policy is approved here.

PRF-to-true-GTO-residual interactions are also unresolved. GTO residuals are
not parent-backed, so parent centers and onsite values do not define that
block.

## Implemented Source Boundary

Implementation is limited to existing owners:

- `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
  for parent/terminal projection, validation, and the compact PRF block;
- `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl` for exact
  factorized PRF one-body blocks;
- `src/cartesian_final_basis_realization/pqs_terminal_ida.jl` for streamed
  parent-backed direct contractions and the bounded parent-IDA comparator;
- `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`
  only for narrow internal module wiring;
- `src/ordinary_coulomb.jl` only for the scalar onsite-calibrated Gaussian
  Coulomb formula and its compact resource.

No root export is approved. A narrow export from the existing internal
`CartesianFinalBasisRealization` module is allowed only if the qualified
consumer/test surface requires it. No new file or module is approved.

Commit `5b46ae073` implemented:

```text
CartesianParentResidualFunctionBlock
build_parent_residual_function_block(...)
parent_residual_one_body_blocks(...)
parent_gaussian_direct_resource(...)
parent_gaussian_direct_blocks(...)
parent_ida_direct_blocks(...)
```

The implementation added no source file or module and leaves ordinary `G-G`
one-body and interaction matrices unchanged.

## Validation Contract

Use existing bounded tests only:

- `test/core/runtests.jl` for the scalar Gaussian direct formula;
- `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` for PRF
  projection, exact one-body blocks, and bounded direct contractions.

Required gates:

1. A compact synthetic metric fixture distinguishes a valid selected target
   from support leakage, terminal contamination, nonfinite input, and rank
   loss. The helper never drops a supplied column.
2. A bounded H2 PQS case uses consumer-selected omitted parent targets and
   verifies support identity, deterministic phases, global `G-R` and `R-R`
   metric identities, and unchanged G construction.
3. PRF kinetic, unit-nuclear, and `H1` blocks agree with a bounded dense parent
   oracle and remain finite/symmetric where applicable.
4. The Gaussian kernel reproduces every tested onsite, is symmetric, finite,
   positive semidefinite within numerical tolerance, and approaches `1/R` in
   the far field.
5. Tiled PRF-PRF and PRF-G results match bounded explicit contractions. Parent
   charge sums remain near one without renormalization, and G-G is unchanged.
6. A bounded parent-IDA comparison reports direct residuals without asserting
   continuum exchange accuracy.
7. Existing shell-q, terminal realization, exact one-body, IDA, and public
   facade gates remain green. An ignored Be2 smoke may be used if H2 does not
   exercise both identity and contracted terminal columns.

No Cr2 fixture, endpoint energy, HF convergence, or exchange assertion belongs
in committed tests.

Accepted validation for `5b46ae073` passed package load, core `440/440`, the
focused H2 gate `358/358`, and the supplemented facade `69/69`. Exact one-body
oracle errors remained near `1e-15`; ordinary `H1` and `Vee` were unchanged.

## Failure And Stop Rules

Stop and report if:

- selected targets cannot be represented as full-rank support-local residuals;
- PRFs are not orthonormal or terminal-orthogonal at the existing numerical
  tolerances;
- exact operators require a dense parent operator or alternate Hamiltonian;
- parent onsite values are nonpositive or disagree with the resolved Coulomb
  expansion;
- parent-backed charge sums are materially nonunit;
- tiled and explicit bounded direct contractions disagree;
- G-G construction or ordinary omitted-PRF behavior changes; or
- use under these IDs alone requires interaction assembly, GTO residual
  blocks, artifacts, solver wiring, or exchange support.

Do not repair a failure by selecting fewer modes, changing a consumer target,
renormalizing charges, flooring metrics, switching Coulomb policy, or adding a
compatibility representation.

## Explicit Non-Goals

These IDs do not approve:

- shell, `source_q`, RDM, mode-count, cutoff, or automatic selection policy;
- public API, facade, driver, default, artifact, sidecar, or solver workflow;
- a combined PRF Hamiltonian or mutation of current `H1`/`Vee`;
- PRF-to-GTO-residual interactions;
- transition-density or exact exchange;
- screened-Hartree, EGOI, residual-GTO, protected-basis, or MWG changes;
- parent-grid, shell-support, terminal-contraction, or Coulomb-policy changes;
- dense parent-by-parent, parent-by-terminal, or final-by-final matrices;
- Cr2-specific source behavior, relaxed-HF interpretation, production energy,
  or paper claims.

The separately approved
[parent-backed injected composition](parent_backed_injected_composition.md)
may consume these implemented objects under its own IDs. It does not broaden
the PRF/direct IDs themselves.
