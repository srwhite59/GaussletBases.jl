# Cartesian Hamiltonian Producer Design

Status: **Slice A, Slice B, Slice C1, and Slice C2 implementation authority; Slice D remains candidate**

This document is the implementation authority for Slice A, Slice B, Slice C1,
and the Slice C2 Hamiltonian construction boundary of the Cartesian/PQS
Hamiltonian producer. It remains the candidate design record for Slice D
driver/artifact simplification. It is intentionally more detailed than a target
card and less historical than the manager running log.

The normative PQS mathematics remains in
`docs/src/algorithms/pqs_shell_construction.md`. This document controls the
implementation shape, allowed objects, pseudocode, line budgets, deletion
requirements, and merge gates.

`AGENTS.md` points here for Cartesian Hamiltonian producer work. Registry
entries explicitly marked **approved** below are binding implementation
authority. Registry entries marked candidate do not permit source
implementation without a later docs-only approval.

## 1. Why this document exists

The Cartesian lane repeatedly accumulated preflights, payload wrappers,
metadata-carried numerical objects, status mirrors, and route-specific duplicate
physics while trying to expose the next blocker. Cleanup through Pass 339 removed
most of that accumulated machinery. The active route is now narrow enough that
new construction should be designed before it is coded.

The process rule is:

> Design may discover the next implementation boundary. Merged source code must
> cross a real numerical or physics boundary.

No source commit may be accepted merely because it adds a more precise blocker,
preflight, status, summary, or planned object.

## 2. Current code state

Current generic terminal/PQS stages already provide:

1. terminal shellification geometry;
2. selected terminal lowering contracts;
3. retained-unit plans;
4. retained-unit transform contracts;
5. ordered terminal support records;
6. boundary-product retained rules and retained counts.

The current source-plan path stops because terminal PQS records lack the three
numerical shell-realization ingredients:

```text
shell_projection
shell_overlap
lowdin_cleanup
```

These ingredients are **not** terminal contract fields today and should not be
added as metadata or staged summaries. Slice A must construct them inside the
terminal realizer from the typed support, retained, transform, and parent metric
objects.

Existing reusable seams are:

```julia
_nested_projected_q_shell_layer(...)
_nested_projected_q_shell_staged_unit_descriptor(...)
CartesianContractedParentMetrics._pqs_shell_realization_plan(descriptor, metrics)
CartesianFinalBasisRealization.pqs_source_shell_realization_final_basis(...)
CartesianIDAHamiltonian(...)
write_cartesian_ida_hamiltonian(...)
```

`_pqs_shell_realization_plan` already computes the shell-support overlap matrix,
but currently does not return it.

The active H2 and Cr2 routes must remain one generic terminal-topology route.
No H2 numerical compatibility route may be reintroduced.

One-center atomic and bond-aligned diatomic PQS routes may differ in
shellification geometry only. Once terminal support, retained, and transform
records exist, both routes must use the same terminal-basis realizer and produce
`CartesianTerminalBasisRealization`.

## 3. Scope

This design ultimately covers the base all-electron PQS Hamiltonian producer:

```text
terminal support and retained contracts
-> terminal localized final basis
-> final-basis kinetic matrix
-> final-basis unit nuclear-attraction matrices by center
-> localized IDA electron-electron matrix
-> CartesianIDAHamiltonian
-> existing minimal artifact writer
```

The following are outside this design and require separate approval:

- residual-GTO or MWG supplement augmentation;
- EGOI or other Hamiltonian corrections;
- RHF/DMRG solver implementation;
- White-Lindsey pair-framework completion;
- distorted-product COMX realization;
- a public change to the visible driver stage sequence.

A supplement request may remain visible, but supplement construction must not
block or complicate the base Hamiltonian producer.

Current authority covers the realized terminal basis, blockwise one-body
operators, localized IDA assembly, and the narrow final Hamiltonian construction
boundary. Driver/materialization simplification, artifact routing, and solver
integration remain future candidate work.

## 4. Binding invariants

### 4.0 Atomic/diatomic unification boundary

Geometry-specific code may produce terminal regions differently:

```text
one-center shellification
bond-aligned diatomic shellification
```

After that boundary, these steps are shared:

```text
terminal lowering
retained rules
terminal basis realization
one-body assembly
IDA assembly
CartesianIDAHamiltonian construction
```

The terminal basis realizer must not inspect system classification, atom count,
route kind, bond axis, or terminal role vocabulary to choose an atomic or
diatomic algorithm. It may dispatch only on terminal lowering/transform kind,
such as direct identity, PQS shell, or unsupported distorted COMX.

Slice A is not complete unless at least one one-center atomic terminal plan and
the H2/Cr2 terminal plans pass through the same realization entry point. No
atomic-specific final-basis object, adapter, or overload is permitted.

If the current one-center route can reach terminal records only through the old
route-skeleton shape input, Slice A must connect that one-center terminal plan
to the typed terminal records before implementation proceeds. It must not
introduce an atomic adapter around the terminal basis realizer.

### 4.1 PQS basis construction

For terminal unit `i`, let `C_i` map unit support rows to its retained localized
columns.

- Direct core/slab/boundary units use an implicit identity map.
- PQS shell units first project the shell candidate out of every already
  retained earlier terminal block in terminal order, then apply shell-local
  Lowdin only within the new shell's retained columns.
- Previous blocks are never rotated by a later shell.
- No global Lowdin is permitted.
- Concatenation is accepted only after explicit cross-block overlap checks.
- Projecting a new shell against previous blocks may put nonzero coefficients
  on previous terminal rows. A block's `support_indices` are therefore its
  effective parent-row coefficient support after projection, not necessarily
  its original terminal region rows.
- A later direct block must be checked against all previously accepted PQS and
  direct blocks. Direct identity does not imply automatic orthogonality to
  earlier blocks when the terminal order interleaves direct and PQS regions.

For every realized PQS shell:

```text
C_i' * S_ii * C_i ~= I
```

For every distinct pair of completed terminal blocks:

```text
C_i' * S_ij * C_j ~= 0
```

Failure of the second condition is
`:terminal_pqs_cross_block_projection_required`; it is not repaired globally.
It is a construction failure that prevents merge of a production slice, not a
new accepted endpoint.

### 4.2 Operator construction

- Base one-body and IDA operators are assembled directly in the final localized
  basis by terminal block pairs.
- A global dense parent/support operator matrix is forbidden in the production
  route.
- A global dense final coefficient matrix is diagnostic/export-only, not the
  working representation.
- Direct identity blocks remain implicit; do not allocate large identity
  matrices.
- Nuclear-attraction matrices remain unit-charge and separated by center until
  consumed by `CartesianIDAHamiltonian`.

### 4.3 Localized IDA convention

The producer stores a one-basis localized IDA Hamiltonian. In the produced final
basis, `electron_electron_ida[a,b]` is the density-density interaction between
localized final basis functions `a` and `b`.

For a later orbital rotation or consumer coefficient matrix `X`, the represented
two-body tensor is:

```text
(pq|rs)_IDA =
    sum_ab (X[a,p] * X[a,q]) *
           electron_electron_ida[a,b] *
           (X[b,r] * X[b,s])
```

Equivalently, the density map is the row-wise pair product
`D[a,p,q] = X[a,p] * X[a,q]`. The base all-electron producer does not export a
separate orbital-to-density transform.

During construction, each final basis column is sign-canonicalized so its
localized IDA weight is positive. One-body matrices and IDA contractions must
use the same canonicalized final basis. A near-zero or nonfinite final IDA
weight is a construction error.

`X` is real in the current producer contract. A future complex-valued consumer
must define the corresponding conjugation convention before it can reuse this
artifact.

### 4.4 Data contracts

- Algorithmic matrices, transforms, rules, source plans, dimensions, and
  coefficients may not be stored in `.metadata`.
- Summaries are disposable views and may not be read as numerical inputs.
- No runtime-generated `NamedTuple{unit_keys}` is allowed in the new producer.
- No result repeats a child object together with child-status and
  child-materialized mirrors.
- The final product is the existing `CartesianIDAHamiltonian`; no broad artifact
  payload is allowed around it.

## 5. Design-change process

### 5.1 Separate design from implementation

A design amendment must be a documentation-only commit or PR. It may not include
`src`, `test`, `bin`, or `tools` changes.

Implementation may begin only from an approved design revision.

### 5.2 Object-first gate

Any proposed new production item in the following categories must first appear
in the registry in Section 6 and be explicitly approved by a frozen design
revision:

- `struct`, abstract type, enum, or exception type;
- module or source file;
- persistent NamedTuple/result shape;
- stage return field;
- public or module-owned helper;
- status or blocker symbol;
- metadata key;
- report or artifact field;
- committed test or probe.

If coding reveals a need for an unlisted item, the agent stops coding and opens a
design-only amendment. The implementation patch must not smuggle the item in and
ask for approval afterward.

### 5.3 Review roles

Before freezing a design revision, obtain three read-only reviews:

1. **Numerical review** — checks PQS projection, Lowdin locality, IDA gauge, and
   center conventions.
2. **Performance review** — checks asymptotic work, peak allocations, dense
   matrix counts, and Cr2-scale viability.
3. **Deletion review** — identifies the exact preflight, payload, adapter,
   summary, and test code that the implementation makes unnecessary.

Reviewers edit or comment on this document only. They do not add scaffolding to
`src`.

### 5.4 Implementation branch rule

Agents may use several exploratory commits on a non-main implementation branch.
Before merge, the branch must be reduced or squashed so the mainline change:

- implements an approved numerical slice;
- switches its live consumer;
- deletes the replaced path in the same merge;
- passes a real physics endpoint;
- stays within the approved added-line budget.

Prototype-only or blocker-only commits do not merge.

## 6. Object and surface registry

This registry records approved surfaces, candidate future surfaces, and
rejected surfaces. Only items explicitly marked **approved** are permitted in
implementation.

### HP-OBJ-01 — `CartesianTerminalBasisBlock` — approved Slice A

Owner:
`CartesianFinalBasisRealization`

Purpose:
One successfully realized terminal basis block. Direct identity blocks are
represented without an allocated identity matrix.

Exact fields:

```julia
struct CartesianTerminalBasisBlock
    unit_key::Symbol
    support_indices::Vector{Int}
    support_states::Vector{NTuple{3,Int}}
    coefficients::Union{Nothing,Matrix{Float64}}
    column_range::UnitRange{Int}
end
```

Meaning:

- `coefficients === nothing` means direct identity on the listed support rows.
- Otherwise `coefficients` is support-local rows by retained columns.
- A block exists only after successful realization; it has no status or blocker
  fields.
- `support_states` may remain only because blockwise Cartesian operator
  assembly uses `(ix, iy, iz)` directly. It is not a second source of support
  authority; it must be derived from `support_indices` and parent axis
  dimensions at construction.
- `support_indices` are the effective coefficient support after projection and
  sign canonicalization. For a PQS shell, they may include rows from previous
  terminal regions.
- `coefficients` for a PQS block are already sign-canonicalized to positive
  localized IDA weights.
- Direct `coefficients === nothing` blocks are valid only after the direct
  support overlap is identity within tolerance and the direct support IDA
  weights are finite and positive.

### HP-OBJ-02 — `CartesianTerminalBasisRealization` — approved Slice A

Owner:
`CartesianFinalBasisRealization`

Purpose:
The block-structured localized final basis consumed by operator assembly.

Exact fields:

```julia
struct CartesianTerminalBasisRealization
    blocks::Vector{CartesianTerminalBasisBlock}
    final_dimension::Int
    max_cross_overlap::Float64
end
```

The object does not store parent bundles, stage objects, summaries, metadata,
global coefficients, or global self-overlap matrices.
`max_cross_overlap` is a diagnostic only and must not become an algorithmic
input.

### HP-RES-01 — terminal basis build result — rejected

Do not introduce a persistent terminal-basis result wrapper. The Slice A
realizer should return `CartesianTerminalBasisRealization` on success.

Expected unsupported topology such as distorted COMX should be rejected before
calling the Slice A PQS realizer. Missing shell-projection/overlap/Lowdin inputs
are implementation defects after Slice A and must not survive as route-state
blockers. Rank loss or failed cross projection is a numerical construction
error, not a mergeable endpoint.

### HP-FILE-01 — terminal realization implementation file — approved Slice A

Path:

```text
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
```

It would be included by the existing `CartesianFinalBasisRealization` module. No
new module is proposed. All source additions in this file count against the
single Slice A budget in Section 9.

### HP-FN-00 — projected terminal shell realization — approved Slice A

Candidate internal helper. It may be file-local inside `HP-FILE-01`.

Purpose:
Realize one PQS terminal shell by explicitly projecting its raw boundary seed
against all previously accepted blocks, forming the projected Gram matrix,
performing shell-local Lowdin, and returning effective support plus final
coefficients. It does not own final IDA sign-gauge canonicalization.

Conceptual signature:

```julia
realize_projected_terminal_shell(
    seed_support_indices,
    seed_coefficients,
    previous_blocks::Vector{CartesianTerminalBasisBlock},
    metrics;
    rank_atol = 1.0e-10,
    identity_atol = 1.0e-8,
    projection_atol = 1.0e-12,
)
```

Normative projection:

```text
B = all previous accepted blocks embedded on their effective support
X = current boundary-mode candidate on its seed support

X_projected = X - B * (B' S B)^(-1) * (B' S X)
G = X_projected' S X_projected
L = G^(-1/2)
C_new = X_projected * L
```

Because previous blocks are intended to be orthonormal, an implementation may
use sequential orthogonal projection:

```text
for previous block C_j:
    residual = C_j' * S * X
    if norm(residual, Inf) <= projection_atol:
        do not subtract and do not enlarge effective support
    else:
        X <- X - C_j * residual
```

It must validate previous-block orthonormality and run a second
reorthogonalization pass or fail if the projection residual remains above
tolerance. The helper returns only data consumed immediately by `HP-FN-01`; it
does not create a public result object.

The default `projection_atol = 1.0e-12` is deliberately tighter than the final
cross-overlap acceptance tolerance. It is the threshold for deciding whether a
candidate is already orthogonal enough to avoid a roundoff-only subtraction
that would artificially grow effective support.

Production Slice A must use the recursively projected coefficients of previous
PQS blocks and their effective supports. It must not approximate previous PQS
blocks by their original shell-local coefficients when projecting later shells.

### HP-FN-01 — terminal basis realizer — approved Slice A

Proposed internal signature:

```julia
pqs_terminal_basis_realization(
    support_records,
    retained_records,
    transform_contracts,
    bundles;
    identity_atol = 1.0e-8,
    cross_atol = 1.0e-8,
)
```

Return:
`CartesianTerminalBasisRealization` on success.

The helper may consume typed objects and vectors. It may not read numerical data
from summaries. Until transform contracts have typed numerical fields, it may
read the already-existing raw source plan and retained rule from the current
contract metadata, but this exception is temporary and must not be extended to
new fields. Shell projection, overlap, and Lowdin must never be added to
metadata.

`HP-FN-01` owns final basis-block construction. It derives support weights from
bundles, validates direct-sector weights, sign-canonicalizes completed block
columns to positive localized IDA weights, and constructs
`CartesianTerminalBasisBlock`.

### HP-FN-02 — cross-block overlap audit — approved Slice A

Purpose:
Compute `C_i' * S_ij * C_j` one block pair at a time and return the largest
infinity norm.

Requirements:

- no global self-overlap matrix;
- no global coefficient matrix;
- no global repair;
- one temporary support cross block at a time;
- symmetric pair traversal only.

### HP-WIRE-01 — generic terminal-basis stage integration — approved Slice A

Owner:
`cartesian_transforms`

Purpose:
Connect Slice A terminal-basis realization to every supported
`:pqs_source_box` terminal plan from typed terminal support, retained, and
transform records.

Contract:

- `cartesian_transforms` owns terminal basis realization for supported PQS
  terminal plans.
- It invokes `pqs_terminal_basis_realization(...)` from typed terminal support,
  retained, and transform records.
- It must not dispatch on system classification, atom count, route kind, bond
  axis, or terminal role names.
- It may dispatch only on terminal lowering/transform kind:
  direct identity, PQS shell, or unsupported distorted COMX.
- It must serve one-center atomic, contact-core diatomic, and separated
  diatomic terminal plans through the same entry point.
- It must not create atomic-specific or diatomic-specific final-basis adapters.

### HP-CHANGE-01 — return shell overlap from existing shell plan — rejected/deferred

Returning `shell_overlap_matrix` from
`CartesianContractedParentMetrics._pqs_shell_realization_plan` may still be a
useful local edit, but it is not sufficient Slice A authority. The required
previous-block projection can enlarge effective support beyond the original
descriptor support, so the projected Gram and overlap diagnostics must be
computed on the effective support used by `HP-FN-00`.

If implementation uses this existing return, it is a helper detail under
`HP-FN-00`, not a standalone approved production surface.

### HP-FN-03 — blockwise final one-body assembly — approved Slice B

Owner:
`CartesianFinalBasisRealization`, in the exact source file
`src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`, included by
`src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`.
This file path is part of the approved Slice B surface. The implementation must
not revive the retired pair-materialization framework.

Purpose:
Fill a dense final-basis matrix from terminal block pairs and caller-supplied 1D
axis matrices without constructing a global support matrix.

Required pre-coding navigation:
Use `docs/src/developer/algorithm_implementation_index.md` before Slice B
implementation work, especially its term-first Coulomb Gaussian contraction and
Gaussian factor reuse entries. Inspect these source anchors before coding:

- `src/ordinary_mapped_backends.jl`
- `src/ordinary_coulomb.jl`
- `src/ordinary_cartesian_ida.jl`
- `src/ordinary_qw_raw_blocks.jl`
- `src/ordinary_qw_operator_assembly.jl`

Direct reuse is preferred when the ordinary helpers match terminal-basis layout.
If the existing ordinary helpers cannot be called directly because their basis,
indexing, or matrix layout is different, Slice B must still reuse their
organization: build and reuse Gaussian factor packets once, then contract
Coulomb terms term-first with the Gaussian expansion index as the short inner
reduction.

Proposed conceptual signature:

```julia
assemble_terminal_product_operator!(
    destination,
    basis::CartesianTerminalBasisRealization,
    axis_x,
    axis_y,
    axis_z;
    scale::Float64 = 1.0,
)
```

The helper handles a single separable Cartesian product term. For terminal
blocks `L` and `R`, it computes local block actions:

```text
C_L' * P_lr * C_R
P_lr[a,b] = axis_x[ix_a, ix_b] * axis_y[iy_a, iy_b] * axis_z[iz_a, iz_b]
```

Direct block coefficients remain implicit identity. Kinetic and
Gaussian-expanded unit nuclear attraction must not locally rebuild one dense
final-basis matrix per Coulomb term as the long-lived implementation pattern.
The approved implementation shape is term-first: for each terminal support-pair
tile, reuse the centered 1D Gaussian factor data and reduce over the Gaussian
expansion index before accumulating the final block. The helper always
accumulates into `destination`; it does not allocate or return a result matrix.
For Gaussian-expanded nuclear attraction, `HP-FN-03` may use private file-local
helpers in the approved Slice B file to perform term-first support-tile
contraction from reusable factor packets. Those helpers are not persistent
production surfaces and must not introduce new result, cache, stage, metadata,
or status objects.

Input exception for the efficient nuclear path:
`assemble_terminal_product_operator!` retains its single-term public signature.
One unexported helper in the same approved file may additionally consume the
Coulomb coefficient vector and three term-first factor arrays:

```julia
_accumulate_terminal_gaussian_sum!(
    destination,
    basis::CartesianTerminalBasisRealization,
    coefficients,
    factors_x,
    factors_y,
    factors_z;
    scale::Float64 = -1.0,
)
```

This helper may be called by validation or future approved orchestration code,
but it must not be stored in a stage, returned as a result, exported, or treated
as route orchestration authority.

Implementation target: **90 source lines**.

No result object is proposed; the destination matrix is the result. `HP-FN-03`
does not approve a `K`/`U_A` payload, stage-return field, report object,
persistent one-body orchestration API, or status vocabulary. If source-level
K/U orchestration proves necessary, implementation stops for a docs-only design
amendment.

Stop rule:
If efficient Slice B implementation requires a new persistent factor-cache or
result type, stage field, metadata key, or orchestration API, stop and request a
docs-only design amendment before coding it. HP-FN-03 remains the only approved
Slice B source surface unless a later amendment explicitly approves another
surface.

Rules:

- consume only `CartesianTerminalBasisRealization` plus the caller-supplied 1D
  factor data: three matrices for the single-product public helper, or the
  coefficient vector plus three term-first factor arrays for the private
  Gaussian-sum helper;
- no atomic/diatomic branches;
- no global support-space operator matrix;
- no global dense coefficient matrix;
- no IDA, Hamiltonian payload/artifact work, or retired
  CartesianPairBlockMaterialization route revival;
- direct blocks are implicit row selectors; do not allocate dense identity
  matrices for them;
- for Slice B use, validate symmetric axis factors before upper-triangular
  block traversal;
- at most one terminal support-pair contraction workspace live;
- tile or stream any simultaneously live local contraction workspace above the
  `64 MiB` cap;
- preserve symmetry by filling both final-basis block triangles from one
  computed upper-triangular block. Do not compute a full nonsymmetric result and
  hide it with final averaging.

### HP-FN-04 — blockwise localized IDA assembly — Slice C1 implemented

HP-FN-04 was implemented as the narrow Slice C1 localized IDA kernel. The
approved surface remains limited to final-basis `electron_electron_ida`
assembly; it does not authorize Hamiltonian construction, route wiring, or
artifact work.

Purpose:
Construct positive-gauge final IDA weights and the final localized
`electron_electron_ida` matrix blockwise from raw PGDG pair-factor terms.

Owner:
`CartesianFinalBasisRealization`, in
`src/cartesian_final_basis_realization/pqs_terminal_ida.jl`, included by
`src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`.
This surface must not be reimplemented in the retired pair-materialization
framework.

Required pre-coding navigation:
Use `docs/src/developer/algorithm_implementation_index.md` before Slice C
implementation work. Inspect these anchors before coding:

- `src/cartesian_nested_faces.jl`:
  `_nested_factorized_weight_aware_pair_terms`,
  `_nested_weight_aware_pair_terms`
- `src/ordinary_qw_raw_blocks.jl`:
  `_qwrg_fixed_block_interaction_matrix`
- `src/cartesian_contracted_parent_metrics/core.jl`:
  `_pqs_source_box_ida_factor_provenance`
- `src/ordinary_mapped_backends.jl`:
  `_mapped_coulomb_expanded_symmetric_matrix`
- `src/ordinary_cartesian_ida.jl`:
  `_ordinary_cartesian_ida_from_pair_factors`,
  `_ordinary_cartesian_ida_from_gausslet_bundle`

The current conclusion is that these are donor/oracle patterns, not directly
callable Slice C kernels for terminal block layout. Slice C must reuse this
organization and convention:

- raw numerator factors come from PGDG `pair_factor_terms_raw`, or from
  `_pqs_source_box_ida_factor_provenance` as a validation/provenance wrapper
  over those same arrays;
- density-normalized `pair_factor_terms` are donor/reference data only for the
  final-density convention and must not replace the raw numerator in the base
  terminal contraction;
- `_mapped_coulomb_expanded_symmetric_matrix` is the fast loop-shape donor:
  terminal Slice C keeps the Gaussian expansion index as the short inner
  reduction, but does not call it directly because it constructs a full ordinary
  Cartesian product matrix;
- `_nested_weight_aware_pair_terms` and `_nested_support_reference_pair_sum` are
  support-set oracle patterns, not terminal block-pair kernels;
- `_ordinary_cartesian_ida_from_pair_factors` is a consumer/convention example;
  its term-by-term `kron` construction is not an approved implementation pattern
  when term-first mapped/terminal contraction is available.

Implementation target: **120 source lines**.

Proposed internal function surface:

```julia
assemble_terminal_ida_interaction!(
    destination,
    basis::CartesianTerminalBasisRealization,
    coefficients,
    raw_pair_terms_x,
    raw_pair_terms_y,
    raw_pair_terms_z,
    weights_x,
    weights_y,
    weights_z;
    weight_atol = 1.0e-12,
    symmetry_atol = 1.0e-10,
)
```

This function is internal to `CartesianFinalBasisRealization`; it is not
exported and does not introduce a result object.

Expected inputs:

- `coefficients` are the Coulomb Gaussian expansion coefficients;
- `raw_pair_terms_x/y/z` are term-first PGDG raw pair-factor tensors with shape
  `(nterms, n_axis, n_axis)`;
- `weights_x/y/z` are ordinary parent axis integral-weight vectors using the
  same convention as Slice A sign canonicalization;
- the helper may accept arrays obtained directly from
  `pgdg.pair_factor_terms_raw` or from the raw-axis fields of
  `_pqs_source_box_ida_factor_provenance`;
- implementation may use three explicit vectors or a fixed three-axis tuple. It
  must not introduce a runtime-keyed axis `NamedTuple` or persistent weight
  payload;
- the helper must not build those tensors per terminal block pair and must not
  allocate any `(nterms, tile_rows, tile_cols)` workspace.

Input validation:

- coefficient count must match all three raw pair tensors;
- all raw tensors must have matching term count and adequate axis dimensions for
  every terminal support state;
- exponent ordering must match the Coulomb expansion used to build the raw pair
  tensors;
- every raw term matrix must be finite and symmetric before upper-triangle
  mirroring can hide bad inputs;
- every axis weight vector must be finite and long enough for the corresponding
  terminal support indices.

Rules:

- do not form a global support raw-pair matrix;
- do not form a global dense coefficient matrix;
- do not introduce all-pairs inventory/status plumbing;
- consume the sign-canonicalized basis produced by Slice A;
- compute final localized IDA weights from the canonicalized terminal basis and
  parent 1D weights at the contraction boundary;
- for each terminal block, compute support product weights from
  `weights_x/y/z`;
- compute final localized IDA weights as

  ```text
  d = C' * support_weights
  ```

  with direct blocks treated as implicit identity;
- a near-zero, nonpositive, or nonfinite final IDA weight is a construction
  error. Residual-Gaussian near-zero weights are not part of this base Slice C;
- do not repeat Slice A sign flips in Slice C. Negative recomputed weights mean
  the supplied basis and weights are inconsistent;
- tile a terminal block pair when one full support-pair workspace would exceed
  the reviewed memory limit;
- use the same **64 MiB** simultaneous local-workspace cap as Slices A/B unless
  a later design amendment changes it;
- reduce over the Gaussian expansion term as the short inner loop for each
  terminal support-pair tile:

  ```text
  for each final support pair (a,b) in the current tile:
      ix, iy, iz = left_state[a]
      jx, jy, jz = right_state[b]
      value = sum(coeff[k] *
                  raw_x[k, ix, jx] *
                  raw_y[k, iy, jy] *
                  raw_z[k, iz, jz])
  ```

- apply the final IDA density gauge only after the raw numerator tile has been
  contracted through the terminal block coefficients and positive final weights;
- for each block pair, first accumulate the unnormalized final block:

  ```text
  B_ij = C_i' * N_ij * C_j
  ```

  after all support tiles have contributed, normalize:

  ```text
  V_ij[p,q] = B_ij[p,q] / (d_i[p] * d_j[q])
  ```

- do not materialize weight-scaled coefficient copies as the primary design;
- do not use `_ordinary_cartesian_ida_from_pair_factors` or any term-by-term
  `kron` construction as production authority;
- no separate IDA payload object is proposed.

Stop rule:
If efficient implementation requires a persistent pair-factor cache/result type,
stage field, metadata key, orchestration API, or status framework, stop and
request a docs-only design amendment before coding it.

### HP-FN-05 — final Hamiltonian producer — Slice C2 approved

C2 approval:
`HP-FN-05` is approved only for constructing the existing
`CartesianIDAHamiltonian` from already assembled Slice B/C1 matrices and
physical center/electron data. It does not approve driver/materialization,
artifact writing, route-stage plumbing, or a Hamiltonian wrapper payload.

Approved construction boundary:

```julia
build_cartesian_ida_hamiltonian(
    kinetic,
    nuclear_attraction_unit_by_center,
    electron_electron_ida,
    nup,
    ndn;
    nuclear_charges,
    nuclear_positions,
)::CartesianIDAHamiltonian{Float64}
```

Implementation may call the existing `CartesianIDAHamiltonian(...)`
constructor directly if no extra source helper is needed. If a helper is added,
it must be internal to `CartesianFinalBasisRealization` and must only validate
the C2 boundary before delegating to the existing constructor.

Inputs:

- `kinetic` is the final-basis kinetic matrix assembled by Slice B;
- `nuclear_attraction_unit_by_center` is a vector of uncharged unit nuclear
  attraction matrices from Slice B, one per physical center;
- `electron_electron_ida` is the final-basis localized IDA matrix from C1;
- `nuclear_charges`, `nuclear_positions`, `nup`, and `ndn` are physical system
  data, not route-status data.

Caller-side provenance requirements:

- before invoking C1, the caller must verify that Coulomb expansion
  coefficients and every raw pair-factor tensor came from the same exponent
  vector in the same order;
- before invoking C2, the caller must verify that `K`, every unit `U_A`, and
  `electron_electron_ida` were built from the same
  `CartesianTerminalBasisRealization`;
- C2 must not recover provenance from metadata or infer it from route labels.

Rules:

- do not apply nuclear charges to `U_A` before construction; the existing
  `CartesianIDAHamiltonian` stores unit attraction matrices and applies charges
  through `one_body_hamiltonian`;
- do not allocate a global support operator or global dense coefficient matrix;
- do not add a result wrapper, route payload, stage field, report field,
  metadata key, status framework, artifact path, or persistent factor cache;
- do not add center-specific or atom/diatomic branches. Center count is data;
- do not revive the retired CartesianPairBlockMaterialization route.

Wrapper/orchestration target: **60 source lines**. If C2 cannot stay as a thin
constructor boundary, stop for a docs-only amendment.

### HP-WIRE-02 — direct materialization Hamiltonian handoff — Slice D candidate

Proposed purpose:
Wire the already realized terminal basis and already approved A/B/C kernels into
the canonical materialization stage so a real `CartesianIDAHamiltonian` is
returned and optionally written by the existing artifact writer.

Proposed boundary:

```julia
cartesian_materialization(
    parent,
    recipe,
    terminal_basis_realization::CartesianTerminalBasisRealization,
    materialization_inputs,
)::Union{Nothing,CartesianIDAHamiltonian{Float64}}
```

The exact implementation may use an equivalent positional boundary, but it must
pass the existing `transforms.terminal_basis_realization` directly from the
driver/harness call site. It may use `parent` and `recipe` as non-summary input
objects for parent bundles, center records, physical charges/positions,
electron counts, and approved recipe controls. If `report` remains in the
signature for compatibility with nearby printing/saving code, it must not be
used as a computational payload except for already non-summary fields that are
also available from `parent`/`recipe`; it must not read target/source-plan or
supplement summaries as construction data.

Return contract:

- return `nothing` when no base Hamiltonian materialization is requested;
- return the `CartesianIDAHamiltonian{Float64}` object itself on success;
- if `save_ham_artifact = true`, write the Hamiltonian with
  `write_cartesian_ida_hamiltonian` and still return the same Hamiltonian;
- never return a `NamedTuple` containing `ida_hamiltonian`, `result_kind`,
  `materialized`, blocker/status mirrors, or durable materialization fields;
- for a requested PQS Hamiltonian, missing terminal basis or missing provenance
  is an error, not a blocked-result wrapper.

Non-PQS behavior:

- preserve the existing White-Lindsey materialization behavior as a separate
  route if it still has live callers, or keep a two-argument compatibility
  method only for that route;
- do not treat `terminal_basis_realization === nothing` as a PQS fallback;
- if a unified signature accepts `Union{Nothing,CartesianTerminalBasisRealization}`,
  it must dispatch cleanly by route family and error on requested PQS
  Hamiltonian construction without a terminal basis.

It must not carry the basis through `cartesian_report`, embed `transforms`,
reconstruct terminal realization from report summaries, or introduce a new
build-input payload.

Candidate rules:

- update live driver/harness/probe callers to pass
  `transforms.terminal_basis_realization` directly;
- materialization owns A/B/C composition for the base Hamiltonian:
  build Slice B `K`, by-center unit `U_A`, Slice C1 `V`, then construct the
  existing `CartesianIDAHamiltonian`;
- enforce C1 coefficient/raw tensor exponent provenance before IDA assembly;
- use existing public `write_cartesian_ida_hamiltonian` for production artifact
  writes. Use `read_cartesian_ida_hamiltonian` only in validation/readback
  checks, not automatically after every production write. Do not invent a new
  artifact shape;
- delete the false blocked-source-plan materialization/report contract in the
  same merge, or report exact live callers that prevent deletion;
- delete or simplify the old materialization wrapper plumbing, including
  `_pqs_source_box_route_driver_durable_materialization`, generic
  materialization serialization, and result-kind/materialized print paths;
- no new report field may duplicate the terminal basis, Hamiltonian matrices,
  raw pair tensors, factor packets, or route stages.

Stop rule:
If Slice D needs to move basis-sized numerical data through report metadata, a
new route payload, recursive stage embedding, a new artifact wrapper, or a
compatibility adapter for the old blocked source-plan path, stop for review.

### HP-OBJ-03 — generic build-result wrapper — rejected

Do not introduce `CartesianHamiltonianBuildResult`, another `Payload`, or a
broad status wrapper around `CartesianIDAHamiltonian`.

At the route boundary, use either:

```text
hamiltonian::CartesianIDAHamiltonian
```

or the two-field result pattern `(value, blocker)` if a blocked route must be
represented in a future approved slice.

### HP-TEST-01 — new committed terminal smoke — rejected

No new committed smoke or preflight test is proposed.

Existing temporary/internal-vocabulary smoke assertions must be deleted or
reduced when the real Hamiltonian endpoint is restored.

## 7. End-to-end pseudocode

### 7.1 Realize the terminal localized basis

```text
function realize_terminal_basis(...):
    verify ordered support, retained, and transform records agree by key
    metrics = axis metrics from parent bundles
    blocks = empty vector
    next_column = 1

    for terminal record in authoritative terminal order:
        if direct identity sector:
            require direct support overlap ~= I
            require direct support IDA weights finite and positive
            r = support_count
            append block with:
                support rows/states
                coefficients = nothing
                column_range = next_column : next_column+r-1

        else if PQS filled-source sector:
            raw_plan = existing transform-contract raw source plan
            retained_rule = existing transform-contract boundary rule

            seed candidate = retained boundary product modes from raw_plan and
                             retained_rule on the terminal shell seed support

            shell_basis = realize_projected_terminal_shell(
                seed_support_indices,
                seed_coefficients,
                previous_blocks = blocks,
                metrics,
            )

            r = shell_basis.retained_count
            append block with:
                support rows/states = shell_basis.effective support
                coefficients = shell_basis.sign-canonicalized coefficients

        else if distorted COMX sector:
            fail before calling this Slice A PQS realizer

        else:
            fail; no compatibility fallback

        next_column += r

    max_cross = audit every off-diagonal block overlap
    if max_cross > cross_atol:
        fail; previous-block projection is not correct enough for production

    return CartesianTerminalBasisRealization(...)
```

The current terminal source-realization preflight and its summary are deleted
when this real object is connected. The realizer does not coexist with the old
preflight mirror.

### 7.2 Assemble final-basis kinetic energy

```text
K = zeros(final_dimension, final_dimension)

for each Cartesian active axis a in x,y,z:
    axis factors = kinetic on a, overlap on other axes
    assemble_terminal_product_operator!(K, basis, factors...)

symmetry-check K
```

Only the final `K` and one block workspace are live.

### 7.3 Assemble unit nuclear attraction by center

```text
U_by_center = one final matrix per nucleus

for center A:
    U_A = zeros(final_dimension, final_dimension)
    build/reuse centered 1D Gaussian factor packet for A

    for each upper-triangular terminal block pair and support tile:
        reduce over Gaussian Coulomb term k as the short inner loop:
            sum -coefficient[k] * Fx[k] * Fy[k] * Fz[k]
        accumulate the reduced tile into U_A

    symmetry-check U_A
```

`U_A` is uncharged `-1/r_A`. Physical charges are passed separately to
`CartesianIDAHamiltonian`.

### 7.4 Localized IDA gauge already fixed by Slice A

```text
inside terminal basis finalization for each block i:
    support_weights_i = product of parent 1D weights on its support states
    final_weights_i = C_i' * support_weights_i

    for each retained column c:
        if abs(weight[c]) <= tolerance or nonfinite:
            fail
        if weight[c] < 0:
            multiply that basis column by -1
            weight[c] = -weight[c]
```

This is column-sign canonicalization, not signed-final-weight division. Slice B
and Slice C consume the canonicalized terminal basis unchanged.

### 7.5 Assemble localized IDA electron-electron matrix

```text
V = zeros(final_dimension, final_dimension)

for upper-triangular terminal block pair (i,j):
    compute d_i and d_j from parent axis product weights:
        d = C' * support_weights
    reject nonfinite, nonpositive, or <= weight_atol entries

    for each tile of terminal support block pair (i,j):
        raw_support_pair_tile = zero tile

        reduce over Gaussian Coulomb term k as the short inner loop:
            raw_support_pair_tile += coefficient[k] *
                raw_pair_x[k] * raw_pair_y[k] * raw_pair_z[k]

        local C_i/C_j rows are block-local coefficients or implicit identity
        accumulate unnormalized B_ij += C_i' * raw_support_pair_tile * C_j

    normalize after all support tiles contribute:
        V_ij[p,q] = B_ij[p,q] / (d_i[p] * d_j[q])
    place V_ij and transpose partner

symmetry-check V
```

For direct blocks, local `C` is implicit identity and must not be allocated. The
weights consumed here are the positive localized final weights fixed by Slice A
canonicalization and recomputed from terminal block coefficients plus parent
axis weights at the Slice C boundary; they are not metadata-carried numerical
data.

### 7.6 Produce the Hamiltonian

```text
ham = CartesianIDAHamiltonian(
    K,
    U_by_center,
    V,
    system.nup,
    system.ndn;
    nuclear_charges = parent.nuclear_charges,
    nuclear_positions = parent.atom_locations,
)
```

The driver may return `ham` or write it with the existing minimal writer. It does
not construct a second artifact payload containing duplicate matrices.

## 8. Performance model and hard limits

Let `N` be the final retained dimension and `b_i` the support size of terminal
block `i`.

Required production-memory shape:

```text
owned final matrices: O(N^2)
block workspace:      O(tile_rows * tile_cols)
working basis:        sum_i O(b_i*r_i) for PQS blocks
```

Projection and audit compute terminal block cross actions:

```text
C_left' * S_lr * C_right
```

They must not construct a global parent overlap matrix or a global final
overlap matrix. At most one terminal support-pair workspace may be live. The
Slice A support-pair workspace cap is **64 MiB** unless a later design amendment
approves a different cap. If a local action exceeds the cap, it must be tiled or
streamed. A dense local pair block is acceptable when bounded; the forbidden
operations are global overlap construction and simultaneous retention of many
pair blocks.

Forbidden shape:

```text
global support matrix: O(M^2), where M is all terminal support rows
global dense C:        O(M*N) as the normal working representation
all pair workspaces retained simultaneously
```

For the current Cr2 preflight estimate `N = 4291`, one dense `Float64` matrix is
about 140.5 MiB. `CartesianIDAHamiltonian` with kinetic, two center matrices, and
IDA interaction owns about 562 MiB before object overhead. For the base
two-center Cr2 fixture:

```text
target peak RSS increase: <= 1.2 GiB
absolute merge cap:       <= 1.5 GiB
raw pair tile workspace:  <= 64 MiB
```

Measure peak RSS of an isolated producer process relative to a package-load
baseline. Report cumulative allocations separately; do not confuse cumulative
allocated bytes with peak live memory.

If `max_ij b_i*b_j` exceeds the reviewed memory budget for one workspace, IDA
assembly must tile that terminal block pair. A full support-pair workspace is
allowed only when the implementation report shows the peak allocation remains
below the reviewed target and cap on the representative Cr2 fixture.

A full dense eigendecomposition is not part of producer construction. H2 may use
it as a bounded physics validation; Cr2 validation should not require it.

## 9. Implementation slices and budgets

Exploratory commits may be separate on a branch. The mergeable slices below must
produce a real consumer-visible result and delete the replaced path.

Slice A, Slice B, and Slice C1 have been implemented under this design. Slice
C2 is approved as a narrow Hamiltonian construction boundary. Slice D remains a
future candidate.

### Slice A — terminal basis realization

Target:
Real terminal basis realization on generic one-center and bond-aligned diatomic
PQS terminal topology, validated on one-center atomic, contact-core H2, and
separated Cr2 records. This slice does not assemble one-body operators, IDA, or
a Hamiltonian artifact.

Added-source target: **150 lines**.
Redesign threshold: **225 added source lines**.
Requirement: net source decrease after deleting the terminal preflight path.

Crossing the redesign threshold returns the work to design. Projection
correctness and numerical validation take priority over meeting the target by
omitting checks.

Must delete in the same merge:

- terminal source-realization preflight implementation;
- terminal preflight summary mirror;
- source-plan fields whose only purpose was reporting that preflight;
- internal-vocabulary assertions made obsolete in the H2 smoke.

Merge validation:

- package load;
- one-center atomic terminal records pass through the same realization entry
  point as diatomics;
- one-center direct sectors pass identity and positive-weight checks;
- one-center PQS shells use previous-block projection and shell-local Lowdin;
- generic H2 terminal basis realizes without the retired compatibility route;
- completed H2 final overlap is near identity;
- the reviewed Cr2 fixture produces a real terminal basis;
- distorted-COMX rejection is allowed only for an input whose typed transform
  inventory actually contains distorted COMX;
- implementation report includes raw cross overlaps, projected cross overlaps,
  recursive effective support sizes, shell ranks, coefficient memory, and
  one-center/H2/Cr2 terminal basis facts.

A branch commit that only moves Cr2 to another blocker does not merge.

The existing atomic multilayer-plan/common-H1 path is migration-only. Slice A
must not delete it and must not add features to it. It should be used only as a
numerical oracle until Slice B or C reproduces the reviewed atomic one-body
endpoint, including hydrogen energy near `-0.5`; then the migration path should
be deleted rather than preserved through adapters.

### Slice B — blockwise final one-body operators

Status: approved Slice B authority.

Target:
Final-basis kinetic `K` and by-center unit nuclear attraction matrices `U_A`
without global support matrices. Slice B does not assemble IDA, does not build a
`CartesianIDAHamiltonian`, and does not write artifacts.

`assemble_terminal_product_operator!` consumes only
`CartesianTerminalBasisRealization`. Atomic and diatomic operator branches are
not permitted.

Added-source target: **150 lines**.
Redesign threshold: **225 added source lines**.
Requirement: net source decrease or a documented net deletion path in the same
producer lane. Crossing the redesign threshold returns the work to design.

Must delete or stop calling:

- any terminal-route dense support-matrix adapter replaced by the blockwise
  kernel;
- duplicate route-local kinetic/nuclear contractions.
- any Slice A/B report or smoke field whose only purpose is to preserve the old
  blocked/operator-preflight vocabulary.

Merge validation:

- source work must cite `docs/src/developer/algorithm_implementation_index.md`
  and report whether it directly reused the ordinary Gaussian-factor/term-first
  helpers or followed their organization because terminal-basis layout prevented
  direct calls;
- before source coding, the implementation target card must name the exact
  one-center/H oracle baseline and tolerance to be used for the unified
  terminal basis check. If no reviewed one-center/H baseline is available,
  establish it first in ignored `tmp/work` code and do not commit source;
- H2 reviewed one-body lowest energy remains
  `-0.7946037173365863` within `1e-10`;
- H2 old/common dense result and new blockwise result agree as a temporary
  cross-check if such an oracle is still live, after which the cross-check
  adapter is deleted;
- N2 light separated-diatomic one-body validation reduces cumulative
  allocations by at least `10x` from the recorded `~18,063 MiB` baseline while
  preserving the H/H2 energies and the `64 MiB` simultaneous local-workspace
  cap;
- Cr2 one-body construction is a later whole-producer stress/performance gate,
  not required to close the first Slice B correctness/performance merge unless
  explicitly requested;
- no IDA, Hamiltonian artifact, residual-GTO, or driver simplification work is
  included;
- no global support-space operator, global dense coefficient matrix, or retired
  CartesianPairBlockMaterialization route is used as production authority;
- no new factor-cache/result type, stage field, metadata key, or orchestration
  API is introduced without a docs-only design amendment;
- implementation report includes final dimensions, matrix symmetry errors,
  finite checks, largest local workspace, elapsed time, and allocation/peak
  memory observations where practical.

### Slice C — localized IDA and `CartesianIDAHamiltonian`

Status: C1 has been implemented under `HP-FN-04`. C2 is approved under
`HP-FN-05` as the narrow construction boundary for the existing
`CartesianIDAHamiltonian`. Driver/materialization/artifact wiring remains Slice
D candidate work.

Target:
Complete base in-memory Hamiltonian producer using the existing
`CartesianIDAHamiltonian` type. Driver/materialization simplification and
artifact routing are deferred to Slice D unless a later Slice C amendment
explicitly changes that boundary.

Slice C1:
localized IDA matrix assembly only. It consumes only
`CartesianTerminalBasisRealization`, Coulomb expansion coefficients, raw axis
pair-factor term tensors derived from the parent bundles, and explicit parent
axis weight vectors. It produces a real final-basis `electron_electron_ida`
matrix. It does not construct `CartesianIDAHamiltonian`.

Slice C2:
approved construction of the existing `CartesianIDAHamiltonian` from Slice B
`K`, by-center unit `U_A`, and Slice C1 `V`. The number of centers is data, not
dispatch. C2 may use the existing constructor directly or a thin internal helper
that only validates the C2 boundary and delegates to that constructor.

Added-source target: **150 lines**. Exceeding the target requires manager
review and a measured memory/validation justification.

Must delete or account for in the same merge:

- any duplicate IDA weight/raw-pair implementation exposed during C2 wiring;
- any new validation scratch code must remain ignored under `tmp/work`;
- blocked source-plan payload surfaces, H2 terminal-stage smoke cleanup, and
  route/materialization report cleanup are deferred to Slice D unless C2 source
  directly touches or supersedes them.

Must not add:

- a new Hamiltonian payload object instead of `CartesianIDAHamiltonian`;
- a global support matrix or global dense coefficient matrix;
- all-pairs inventory/status plumbing;
- metadata-carried numerical pair data;
- CPBM route authority;
- a persistent pair-factor cache/result type, stage field, report field, or
  orchestration API without docs-only amendment.

Merge validation:

- C1 first correctness gate: H2 one-body energy plus reviewed H2
  self-Coulomb/localized IDA parity. The reviewed H2 self-Coulomb value is
  `0.4569117646737212`;
- an ignored H2 oracle comparison against the still-present complete-core/shell
  IDA formula is allowed for validation, but not production authority;
- Cr2 is not the first Slice C correctness gate;
- a light separated diatomic may be used as the first topology/performance smoke
  after H2 parity, but not as a replacement for H2 physics parity;
- C2 validation checks that constructed `CartesianIDAHamiltonian` has physically
  correct center charges, positions, electron counts, finite/symmetric `K`, all
  unit `U_A`, and finite/symmetric `electron_electron_ida`;
- C2 validation checks `one_body_hamiltonian(ham)` reproduces the Slice B H2
  H1 matrix and lowest eigenvalue;
- C2 validation checks the H2 self-Coulomb computed from
  `ham.electron_electron_ida` remains at the reviewed value
  `0.4569117646737212`;
- `electron_electron_ida` is finite and symmetric;
- final localized IDA weights are positive and finite;
- implementation report states which IDA anchors were directly reused versus
  used only as donor patterns;
- implementation report states whether raw pair-factor terms came directly from
  `pgdg.pair_factor_terms_raw` or through
  `_pqs_source_box_ida_factor_provenance`, and confirms they were not rebuilt
  inside terminal block loops;
- implementation report states how C1 caller-side exponent provenance was
  checked before IDA assembly;
- Cr2 base Hamiltonian construction is a later stress/performance gate, not a
  blocker for the first H2 Slice C correctness merge.

### Slice D — driver simplification

Status: future candidate, not approved by the current Slice A/B/C authority.

Target:
Make the canonical materialization stage return/write the real Hamiltonian with
no duplicate payload/report graph.

The final producer receives the same `CartesianTerminalBasisRealization` for
atoms and diatomics by direct handoff from `cartesian_transforms` to
`cartesian_materialization`.

Required handoff:

```julia
transforms = cartesian_transforms(units, recipe)
...
materialization = cartesian_materialization(
    parent,
    recipe,
    transforms.terminal_basis_realization,
    materialization_inputs,
)
```

The signature may be adjusted only to preserve this dataflow. Slice D must not:

- embed `terminal_basis_realization` in `cartesian_report`;
- pass the full `transforms` stage recursively into materialization;
- reconstruct terminal realization from report summaries;
- introduce a new Hamiltonian build-input payload;
- place terminal basis, matrices, raw pair tensors, or factor packets in
  metadata.

`cartesian_materialization` return contract:

- `nothing` when no base Hamiltonian is requested;
- `CartesianIDAHamiltonian{Float64}` on successful base Hamiltonian
  materialization;
- the same Hamiltonian object when artifact writing is requested;
- no materialization `NamedTuple`, `result_kind`, `materialized`,
  `ida_hamiltonian` wrapper field, or blocked-status mirror.

Production work:

1. Materialization constructs Slice B `K` and by-center uncharged unit `U_A`
   from the handed-off terminal basis.
2. Materialization verifies C1 coefficient/raw-tensor exponent provenance, then
   constructs Slice C1 `electron_electron_ida`.
3. Materialization constructs the existing `CartesianIDAHamiltonian` directly.
4. If artifact output is requested, materialization uses the existing
   `write_cartesian_ida_hamiltonian` contract. Validation may read the artifact
   back with `read_cartesian_ida_hamiltonian`; production should not
   automatically perform readback. It must not add a new artifact shape.

Deletion target:

- delete the obsolete blocked source-plan payload/report contract whose only
  remaining message is that `source_plan = nothing`;
- delete or stop publishing report fields such as
  `physical_gausslet_source_plan_summary` when they only preserve the old
  blocked source-plan story;
- audit `physical_gausslet_target_*` and supplement-preflight report fields and
  either delete them in the same merge or report exact live callers that still
  require them;
- delete or simplify `_pqs_source_box_route_driver_durable_materialization`,
  generic report-artifact serialization of materialization data, and
  old `result_kind`/`materialized` printing;
- remove H2 smoke/probe assertions that inspect obsolete blocker/status
  vocabulary after the real Hamiltonian endpoint is available.

Added-source target: **80 lines**, with a net source decrease required.

No new driver stage, public type, report field cloud, or artifact wrapper is
proposed.

Validation:

- package load;
- H2 terminal smoke updated to validate real materialization rather than the old
  blocked source-plan contract;
- H2 materialization constructs `CartesianIDAHamiltonian{Float64}` with final
  dimension `471`;
- `one_body_hamiltonian(ham)` reproduces H2 H1 lowest
  `-0.79460371733658908`;
- `ham.electron_electron_ida` reproduces H2 self-Coulomb
  `0.4569117646737212`;
- requested artifact output writes with existing `write_cartesian_ida_hamiltonian`;
- validation-only artifact readback uses existing
  `read_cartesian_ida_hamiltonian`;
- light separated diatomic may be used as a topology/performance smoke after H2
  passes. Cr2 remains a later stress/performance gate unless explicitly
  requested.

## 10. Test policy

No new committed test file is proposed by this design.

Acceptance should reuse or replace existing tests so the durable checks are
physical:

- one-center terminal-basis realization through the same entry point as H2/Cr2;
- hydrogen energy;
- H2 one-body energy and localized IDA quantity;
- symmetry and finiteness of physical matrices;
- real `CartesianIDAHamiltonian` construction;
- artifact consumption only when Slice D is approved.

Do not preserve tests whose main purpose is asserting internal blocker symbols,
preflight names, field presence, terminal role vocabulary, or metadata flags.
Temporary debugging belongs in ignored `tmp/work` scripts.

### 10.1 Optional uncommitted numerical spike

During Slice A implementation or future design amendments, an ignored
`tmp/work` script may be used to measure the projection design. It must not be
committed or treated as implementation authority.

The spike should report:

- one-center, H2, and Cr2 terminal records used;
- raw `B_previous' S X` cross overlaps;
- projected cross overlaps after previous-block projection;
- effective support sizes for each shell;
- projected shell ranks and retained counts;
- coefficient memory estimate;
- allocation/RSS observations when practical.
- whether the calculation used true recursively projected previous blocks or a
  shell-local approximation.

If the spike shows that previous-block projection produces unexpectedly dense
effective supports or unstable ranks, Slice A returns to design before source
implementation.

## 11. Mechanical implementation gate

Every implementation PR must list the approved design IDs it uses. For the
current Slice A/B/C authority, the approved ID set is:

```text
Design IDs: HP-OBJ-01, HP-OBJ-02, HP-FILE-01, HP-FN-00, HP-FN-01, HP-FN-02, HP-WIRE-01, HP-FN-03, HP-FN-04, HP-FN-05
```

Review added lines before scientific review:

```bash
git diff --check

git diff --numstat BASE..HEAD -- src bin tools test docs

git diff -U0 BASE..HEAD -- src bin tools test |
  grep '^+' |
  grep -Ev '^\+\+\+' |
  grep -nE '(^\+.*(struct|abstract type|@enum|module)|NamedTuple\{|\.metadata|get\(.*metadata|haskey\(.*metadata|_materialized|status.*=|blocker.*=|catch$|catch err|::Any|Payload|summary.*=)' || true
```

For every suspicious addition, the PR must cite an approved design ID from the
frozen revision. An uncited new object or surface is an automatic rejection.

The PR must also report:

```text
added src lines:
deleted src lines:
objects added:
objects deleted:
metadata keys added:
status/blocker symbols added:
tests added:
peak allocation/memory:
physics endpoints:
```

## 12. Design amendment template

```text
Design amendment ID:
Problem demonstrated:
Why existing approved objects/functions are insufficient:
Exact proposed object/function/file:
Exact fields or signature:
Owner module/file:
Added-source budget:
Code deleted or simplified by the addition:
CPU scaling:
Memory scaling:
Physics endpoint that consumes it:
Why this is not metadata, a status mirror, or a compatibility adapter:
```

An amendment with no concrete deletion/simplification target and no physics
consumer should normally be rejected.

## 13. Slice A implementation watchpoints

The final recursive spike answered the pre-freeze numerical questions well
enough to approve Slice A. Implementation must still report:

1. whether one-center atomic, contact-core H2, and separated Cr2 terminal plans
   realize through the same entry point within the 150-line target and below the
   225-line redesign threshold;
2. which current H2 smoke assertions were deleted when the real terminal basis
   returned;
3. effective support sizes, shell ranks, cross-overlap residuals, and workspace
   sizes for one-center, H2, and Cr2;
4. whether one-center public staging was connected to typed terminal records
   without freezing the old route-skeleton shape as the intended contract.

No source implementation may add unlisted production surfaces while answering
these questions. Any new production surface requires a prior docs-only
amendment.
