# Cartesian Hamiltonian Producer Design

Status: **draft v2 candidate for milestone review; not yet an implementation authority**

This document is the proposed implementation authority for completing the
Cartesian/PQS Hamiltonian producer from the current terminal-topology route to a
`CartesianIDAHamiltonian`. It is intentionally more detailed than a target card
and less historical than the manager running log.

The normative PQS mathematics remains in
`docs/src/algorithms/pqs_shell_construction.md`. This document controls the
implementation shape, allowed objects, pseudocode, line budgets, deletion
requirements, and merge gates.

After review and explicit manager/user approval, `AGENTS.md` should point here
and the object registry below becomes binding. Until then, agents may edit this
document but must not treat it as permission to change `src`. Registry entries
are candidates unless this document is later frozen as implementation authority.

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

## 3. Scope

This design covers the base all-electron PQS Hamiltonian producer:

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

## 4. Binding invariants

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

## 6. Candidate object and surface registry

This draft records candidate surfaces and rejected surfaces. It is not yet
implementation authority. Only items later marked **approved** in a frozen
design revision are permitted in implementation. Items marked **candidate**
require another design review before coding.

### HP-OBJ-01 — `CartesianTerminalBasisBlock` — candidate

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
  authority; it must be derived from `support_indices` and `parent_dims` at
  construction.

Definition and validation target: **25 added source lines**.

### HP-OBJ-02 — `CartesianTerminalBasisRealization` — candidate

Owner:
`CartesianFinalBasisRealization`

Purpose:
The block-structured localized final basis consumed by operator assembly.

Exact fields:

```julia
struct CartesianTerminalBasisRealization
    blocks::Vector{CartesianTerminalBasisBlock}
    parent_dims::NTuple{3,Int}
    final_dimension::Int
    max_cross_overlap::Float64
end
```

The object does not store parent bundles, stage objects, summaries, metadata,
global coefficients, or global self-overlap matrices.
`max_cross_overlap` is a diagnostic only and must not become an algorithmic
input.

Definition and validation target: **20 added source lines**.

### HP-RES-01 — terminal basis build result — candidate

Persistent result shape:

```julia
(
    basis::Union{Nothing,CartesianTerminalBasisRealization},
    blocker::Union{Nothing,Symbol},
)
```

No `status`, `materialized`, `summary`, `metadata`, `next_blocker`, or child
mirrors are permitted.
Invariant: `basis === nothing` if and only if `blocker !== nothing`.

Candidate blocker vocabulary:

```text
:missing_terminal_shell_projection
:missing_terminal_shell_overlap
:missing_terminal_shell_lowdin_cleanup
:terminal_pqs_cross_block_projection_required
:distorted_product_realization_missing
```

No additional blocker symbol may be added without a design amendment.

### HP-FILE-01 — terminal realization implementation file — candidate

Path:

```text
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
```

It would be included by the existing `CartesianFinalBasisRealization` module. No
new module is proposed.

Target file size in the first accepted implementation: **150 source lines**,
including object definitions and all helpers. Exceeding the target requires
manager review; numerical projection correctness takes priority over meeting
the target by omitting validation.

### HP-FN-01 — terminal basis realizer — candidate

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
`HP-RES-01`.

The helper may consume typed objects and vectors. It may not read numerical data
from summaries. Until transform contracts have typed numerical fields, it may
read the already-existing raw source plan and retained rule from the current
contract metadata, but this exception is temporary and must not be extended to
new fields. Shell projection, overlap, and Lowdin must never be added to
metadata.

Implementation target within `HP-FILE-01`: **85 source lines**.

### HP-FN-02 — cross-block overlap audit — candidate

Purpose:
Compute `C_i' * S_ij * C_j` one block pair at a time and return the largest
infinity norm.

Requirements:

- no global self-overlap matrix;
- no global coefficient matrix;
- no global repair;
- one temporary support cross block at a time;
- symmetric pair traversal only.

Implementation target within `HP-FILE-01`: **35 source lines**.

### HP-CHANGE-01 — return shell overlap from existing shell plan — candidate

Allowed modification:
Add `shell_overlap_matrix` to the existing return value of
`CartesianContractedParentMetrics._pqs_shell_realization_plan`.

Target budget: **2 added source lines**.

It must be a normal returned field, not metadata.

### HP-FN-03 — blockwise final one-body assembly — candidate

Candidate owner:
`src/pqs_multilayer_complete_core_shell_h1.jl` or a smaller existing
`CartesianFinalBasisRealization` operator file. This must be resolved before
approval and must not revive the retired pair-materialization framework.

Purpose:
Fill a dense final-basis matrix from terminal block pairs and caller-supplied 1D
axis matrices without constructing a global support matrix.

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

The helper handles a single separable Cartesian product term. Kinetic and
Gaussian-expanded nuclear attraction call it repeatedly. The helper always
accumulates into `destination`; it does not allocate or return a result matrix.

Implementation target: **90 source lines**.

No result object is proposed; the destination matrix is the result.

### HP-FN-04 — blockwise IDA assembly — candidate

Purpose:
Construct positive-gauge final IDA weights and the final localized
`electron_electron_ida` matrix blockwise from raw PGDG pair-factor terms.

Candidate owner:
unresolved. The owner must be chosen before approval and must not be the retired
pair-materialization framework.

Implementation target: **120 source lines**.

Rules:

- do not form a global support raw-pair matrix;
- canonicalize a nonzero column sign before positive-weight division rather than
  using signed-final-weight division;
- a near-zero or nonfinite final weight is a construction error;
- tile a terminal block pair when one full support-pair workspace would exceed
  the reviewed memory limit;
- no separate IDA payload object is proposed.

### HP-FN-05 — final Hamiltonian producer — candidate

Conceptual signature:

```julia
build_cartesian_ida_hamiltonian(
    basis::CartesianTerminalBasisRealization,
    parent,
    system,
    coulomb_expansion,
)::CartesianIDAHamiltonian{Float64}
```

It allocates only the final owned matrices required by
`CartesianIDAHamiltonian` plus bounded block workspaces.

Wrapper/orchestration target: **60 source lines**.

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
            r = support_count
            append block with:
                support rows/states
                coefficients = nothing
                column_range = next_column : next_column+r-1

        else if PQS filled-source sector:
            raw_plan = existing transform-contract raw source plan
            retained_rule = existing transform-contract boundary rule
            projection_basis = all previously accepted terminal blocks

            descriptor = projected q-shell descriptor from:
                outer_box
                inner_exclusion_box
                source_mode_shape
                bond_axis
                parent bundles
                projection_basis

            descriptor/projection must remove the candidate shell from the span
            of every previous terminal block before shell-local Gram/Lowdin

            shell_plan = existing _pqs_shell_realization_plan(descriptor, metrics)

            shell_basis = existing pqs_source_shell_realization_final_basis(
                raw_plan,
                retained_rule,
                shell_support_indices = record.support_indices,
                shell_overlap = shell_plan.shell_overlap_matrix,
                shell_projection = shell_plan.shell_projection_matrix,
                lowdin_cleanup = shell_plan.lowdin_cleanup,
            )

            require shell_basis available
            r = shell_basis.final_retained_count
            append block with shell_basis.final_shell_coefficients

        else if distorted COMX sector:
            return blocker :distorted_product_realization_missing

        else:
            fail; no compatibility fallback

        next_column += r

    max_cross = audit every off-diagonal block overlap
    if max_cross > cross_atol:
        return blocker :terminal_pqs_cross_block_projection_required

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
    build centered 1D Gaussian factor matrices for A

    for Gaussian Coulomb term k:
        assemble_terminal_product_operator!(
            U_A,
            basis,
            factor_x[k], factor_y[k], factor_z[k],
            scale = -coefficient[k],
        )

    symmetry-check U_A
```

`U_A` is uncharged `-1/r_A`. Physical charges are passed separately to
`CartesianIDAHamiltonian`.

### 7.4 Fix the localized IDA gauge

```text
for each terminal block i:
    support_weights_i = product of parent 1D weights on its support states
    final_weights_i = C_i' * support_weights_i

    for each retained column c:
        if abs(weight[c]) <= tolerance or nonfinite:
            fail
        if weight[c] < 0:
            multiply that basis column by -1
            weight[c] = -weight[c]
```

This is column-sign canonicalization, not signed-final-weight division.

The sign choice must be applied consistently to later one-body and IDA
contractions.

### 7.5 Assemble localized IDA electron-electron matrix

```text
V = zeros(final_dimension, final_dimension)

for upper-triangular terminal block pair (i,j):
    for each tile of terminal support block pair (i,j):
        raw_support_pair_tile = zero tile

        for Gaussian Coulomb term k:
            raw_support_pair_tile += coefficient[k] *
                product_pair_factor_tile(i, j, k)

        W_i_tile = local C_i rows divided by positive final weights_i
        W_j_tile = local C_j rows divided by positive final weights_j

        accumulate V_ij += W_i_tile' * raw_support_pair_tile * W_j_tile
    place V_ij and transpose partner

symmetry-check V
```

For direct blocks, local `C` is implicit identity and must not be allocated.

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
block workspace:      O(max_ij b_i*b_j)
working basis:        sum_i O(b_i*r_i) for PQS blocks
```

Forbidden shape:

```text
global support matrix: O(M^2), where M is all terminal support rows
global dense C:        O(M*N) as the normal working representation
all pair workspaces retained simultaneously
```

For the current Cr2 preflight estimate `N = 4291`, one dense `Float64` matrix is
about 140.5 MiB. `CartesianIDAHamiltonian` with kinetic, two center matrices, and
IDA interaction owns about 562 MiB before object overhead. The producer should
target peak memory below **1.2 GiB** for the base two-center Cr2 fixture and must
report measured peak allocation before being called production-ready.

If `max_ij b_i*b_j` exceeds the reviewed memory budget for one workspace, IDA
assembly must tile that terminal block pair. A full support-pair workspace is
allowed only when the implementation report shows the peak allocation remains
below the reviewed limit on the representative Cr2 fixture.

A full dense eigendecomposition is not part of producer construction. H2 may use
it as a bounded physics validation; Cr2 validation should not require it.

## 9. Implementation slices and budgets

Exploratory commits may be separate on a branch. The mergeable slices below must
produce a real consumer-visible result and delete the replaced path.

### Slice A — terminal basis and H2 physics restoration

Target:
Real terminal basis realization on the generic H2/Cr2 topology plus restoration
of the H2 one-body physics endpoint through existing common kernels where safe.

Added-source target: **150 lines**. Exceeding the target requires manager
review, but projection correctness and numerical validation take priority over
meeting the target by omitting checks.

Must delete in the same merge:

- terminal source-realization preflight implementation;
- terminal preflight summary mirror;
- source-plan fields whose only purpose was reporting that preflight;
- internal-vocabulary assertions made obsolete in the H2 smoke.

Merge validation:

- package load;
- H atom one-body energy remains near `-0.5`;
- generic H2 route restores its reviewed one-body lowest energy;
- completed H2 final overlap is near identity;
- Cr2 produces a real terminal basis or stops only at a reviewed
  distorted-product blocker.

A branch commit that only moves Cr2 to another blocker does not merge.

### Slice B — blockwise final one-body operators

Target:
Final-basis kinetic and by-center unit nuclear attraction without global support
matrices.

Added-source target: **150 lines**. Exceeding the target requires manager
review and a deletion/simplification explanation.

Must delete or stop calling:

- any terminal-route dense support-matrix adapter replaced by the blockwise
  kernel;
- duplicate route-local kinetic/nuclear contractions.

Merge validation:

- H and H2 reviewed one-body energies unchanged within tolerance;
- H2 old/common dense result and new blockwise result agree as a temporary
  cross-check, after which the cross-check adapter is deleted;
- Cr2 final K and center matrices are finite and symmetric;
- measured allocation and peak-memory report.

### Slice C — localized IDA and `CartesianIDAHamiltonian`

Target:
Complete base Hamiltonian producer and existing minimal artifact output.

Added-source target: **150 lines**. Exceeding the target requires manager
review and a measured memory/validation justification.

Must delete in the same merge:

- blocked source-plan payload surfaces superseded by the real Hamiltonian;
- temporary H2 terminal-stage smoke assertions that inspect blocker/status
  vocabulary;
- any duplicate IDA weight/raw-pair implementation exposed during wiring.

Merge validation:

- H2 one-body energy and reviewed self-Coulomb/IDA endpoint;
- constructed `CartesianIDAHamiltonian` has physically correct center charges,
  positions, and electron counts;
- artifact write/read round trip preserves all matrices;
- Cr2 base Hamiltonian construction completes within the reviewed memory limit,
  or the implementation does not merge.

### Slice D — driver simplification

Target:
Make the canonical materialization stage return/write the real Hamiltonian with
no duplicate payload/report graph.

Added-source target: **80 lines**, with a net source decrease required.

No new driver stage, public type, report field cloud, or artifact wrapper is
proposed.

## 10. Test policy

No new committed test file is proposed by this design.

Acceptance should reuse or replace existing tests so the durable checks are
physical:

- hydrogen energy;
- H2 one-body energy and localized IDA quantity;
- symmetry and finiteness of physical matrices;
- real `CartesianIDAHamiltonian` construction and artifact consumption.

Do not preserve tests whose main purpose is asserting internal blocker symbols,
preflight names, field presence, terminal role vocabulary, or metadata flags.
Temporary debugging belongs in ignored `tmp/work` scripts.

## 11. Mechanical implementation gate

After a future frozen design revision, every implementation PR must list the
approved design IDs it uses, for example:

```text
Design IDs: HP-OBJ-01, HP-OBJ-02, HP-FN-01, HP-FN-02, HP-CHANGE-01
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

## 13. Open review questions

The three design reviewers should resolve these before implementation begins:

1. Can Slice A restore the H2 physics endpoint within the 150-line target while
   keeping the terminal basis representation generic for Cr2?
2. Is the one-basis localized IDA convention in Section 4.3 the correct public
   mathematical contract for `CartesianIDAHamiltonian`?
3. Which existing operator file should own `HP-FN-03` and `HP-FN-04` without
   reviving the inactive pair-materialization framework?
4. Can previous-block projection be implemented through the current
   projected-shell descriptor path, or is a small explicit projection operation
   needed before approval?
5. Which current H2 smoke assertions can be deleted immediately when the real
   basis and Hamiltonian return?
6. Is the 1.2 GiB Cr2 peak-memory ceiling appropriate for the intended producer
   machine, or should v2 use a smaller target plus an absolute hard cap?

No source implementation begins until these questions are answered in this
document or explicitly deferred by the user/manager.
