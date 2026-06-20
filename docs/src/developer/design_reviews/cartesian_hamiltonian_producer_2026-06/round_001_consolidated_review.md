# Round 001 Consolidated Review

Reviewer: repo-doer@macmini

Date: 2026-06-20

## Verdict

needs design revision before ChatGPT-Pro review

The draft is directionally right: it rejects global Lowdin cleanup, keeps the
public product as `CartesianIDAHamiltonian`, and treats residual-GTO/MWG work as
future generic-producer work rather than resurrecting the deleted H2 private
lane. It is not ready to become implementation authority because several
approved registry items are broader than the reviewed numerical boundary, and
the shell-projection and IDA-gauge contracts are still underspecified at the
point where incorrect physics could be implemented cleanly.

## Blocking Issues

1. The shell-realization algorithm does not state how each PQS shell candidate
   is projected out of all already-retained earlier terminal blocks. It says to
   construct a projected q-shell descriptor and then audit cross-block overlaps.
   Auditing after the fact is necessary but not sufficient; the design must
   specify the projection authority that should make the audit pass.

2. `HP-FN-03`, `HP-FN-04`, and `HP-FN-05` are marked approved even though file
   ownership, exact inputs, tiling strategy, and numerical conventions are not
   settled. They should be candidate surfaces until v2 resolves those details.

3. The localized IDA gauge section needs a precise convention. The document
   states `final_weights_i = C_i' * support_weights_i` and then divides columns
   by positive final weights. That may be the intended localized IDA convention,
   but it must be explicitly tied to the electron-electron matrix convention
   consumed by `CartesianIDAHamiltonian`; otherwise a future implementation can
   produce correctly shaped but physically wrong `V`.

4. The performance model forbids global support matrices, but `HP-FN-04` still
   allows a full `raw_support_pair` for the largest terminal block pair. For the
   Cr2 estimate, a single 4802 x 4802 block is about 184 MiB before factor
   matrices and output matrices. The design needs an explicit block/tiling rule
   or a measured reason why one full support-pair workspace remains acceptable.

5. The registry uses "approved" language inside a document that is still
   explicitly draft-only. For round 001, every positive HP-* item should be
   treated as "candidate for approval in v2" unless the review explicitly keeps
   it as approved.

## Numerical Review

### Findings

- The no-global-Lowdin invariant is correct and should remain binding.
- Direct terminal sectors as identity blocks are correct, but the identity
  convention must say whether direct slabs and boundary slabs are accepted as
  final orthonormal rows only because parent rows are orthonormal in the working
  parent basis.
- The shell-local Lowdin sequence is incomplete without an explicit "project
  candidate against accumulated prior blocks" step. The current pseudocode can
  be read as using only the raw shell descriptor's inner exclusion geometry.
- `:terminal_pqs_cross_block_projection_required` is the right failure mode,
  but the design should say it is a construction failure, not a new mergeable
  blocker endpoint.
- The by-center unit nuclear attraction convention is correct: `U_A` is
  uncharged, and charges remain in `CartesianIDAHamiltonian`.
- The IDA convention is the weakest numerical section. The document should add
  the orbital/density relation represented by `electron_electron_ida` and how
  column signs affect both one-body and IDA contractions.

### Proposed Document Edits

In Section 4.1, replace:

```text
- PQS shell units use shell projection followed by shell-local Lowdin.
```

with:

```text
- PQS shell units first project the shell candidate out of every already
  retained earlier terminal block in terminal order, then apply shell-local
  Lowdin only within the new shell's retained columns. Earlier blocks are never
  changed by this projection.
```

In Section 7.1, insert before `shell_plan = existing ...`:

```text
            projection_basis = all previously accepted terminal blocks
            descriptor/projection must remove the candidate shell from
            projection_basis before shell-local Gram/Lowdin
```

In Section 7.5, add the convention:

```text
The stored `electron_electron_ida` is the localized IDA density-density matrix
in the same basis as `K` and `U_A`. For an orbital rotation X applied later by a
consumer, the induced two-body tensor must be documented as
sum_ab (D_a,ij) V_ab (D_b,kl), where D is the density map implied by the
localized IDA basis. This producer does not export a separate T map for the
base all-electron basis.
```

If that relation is not the intended convention, the design must replace it
with the exact intended relation before implementation.

## Performance Review

### Findings

- The owned final matrices are correctly identified as the dominant durable
  memory: about 562 MiB for Cr2 with `N = 4291` and two centers.
- The 1.2 GiB peak-memory ceiling is plausible only if support-pair workspaces,
  one-body factor matrices, and temporary contractions are not retained across
  block pairs.
- `HP-FN-03` needs an explicit accumulation interface. The pseudocode later
  calls it with `scale = -coefficient[k]`, but the approved signature does not
  include `scale`, `alpha`, or additive mode.
- `HP-FN-04` likely needs tiling or row/column chunking for large support block
  pairs. One full `raw_support_pair` at a time may be acceptable for H2 and
  small Cr2, but the design should make that an explicit measured allowance,
  not an accidental implementation choice.
- The line budgets for Slice A and the terminal realizer look too tight if the
  implementation must include key joins, validation, previous-block projection,
  direct/PQS/distorted sectors, and cross-block auditing. Tight budgets are
  useful review pressure, but they must not force shortcuts around projection.

### Proposed Document Edits

In `HP-FN-03`, replace the signature with:

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

Add:

```text
The helper always accumulates into `destination`; it does not allocate or return
a result matrix.
```

In Section 8, add:

```text
If `max_ij b_i*b_j` exceeds the reviewed memory budget for one workspace, IDA
assembly must tile that terminal block pair. A full support-pair workspace is
allowed only when the implementation report shows the peak allocation remains
below the approved limit on the representative Cr2 fixture.
```

Change the Slice A line-budget language to:

```text
Added-source target: 150 lines. Exceeding the target requires manager review,
but numerical projection correctness takes priority over meeting the budget by
omitting validation.
```

## Minimality And Deletion Review

### Findings

- The deletion bias is good: each slice names old preflight/payload/test
  surfaces that should be removed.
- `HP-OBJ-01` may store redundant `support_indices` and `support_states`.
  Keeping both is acceptable only if operator assembly needs states without
  repeated index decoding. The design should say that explicitly.
- `HP-OBJ-02.max_cross_overlap` is a compact construction diagnostic and can
  stay, but it should not become an algorithmic input.
- `HP-FILE-01` is acceptable as a candidate file, but adding a new file should
  be tied to deleting the current terminal source-realization preflight, not
  merely moving code.
- `HP-TEST-01` as a rejected item is useful. Keep it as a negative registry
  guard.
- The H2 smoke should be deleted or reduced again once Slice A restores the real
  H2 terminal basis. The draft already says this, and v2 should name the exact
  smoke assertions to remove.

### Proposed Document Edits

Under `HP-OBJ-01`, add:

```text
`support_states` may remain only because blockwise Cartesian operator assembly
uses `(ix, iy, iz)` directly. It is not a second source of support authority;
it must be derived from `support_indices` and `parent_dims` at construction.
```

Under Slice A "Must delete", add:

```text
- the current H2 terminal-stage smoke assertions that inspect
  `source_plan_status`, `source_plan_blocker`, and terminal preflight names,
  once a real basis object is returned.
```

## Registry Review

- `HP-OBJ-01`: keep but demote to candidate until `support_states` derivation
  and direct-sector orthonormality assumptions are stated.
- `HP-OBJ-02`: keep as approved candidate; fields are compact enough.
- `HP-RES-01`: keep as approved candidate, but add invariant
  `basis === nothing` iff `blocker !== nothing`.
- `HP-FILE-01`: keep but demote to candidate until Slice A deletion target is
  tied to this file.
- `HP-FN-01`: keep but demote to candidate until previous-block projection is
  specified.
- `HP-FN-02`: keep as approved candidate.
- `HP-CHANGE-01`: keep as approved candidate if it is truly only two lines and
  returns a normal field, not metadata.
- `HP-FN-03`: keep but demote to candidate; signature and owner are
  underspecified.
- `HP-FN-04`: keep but demote to candidate; IDA convention and tiling strategy
  are underspecified.
- `HP-FN-05`: keep but demote to candidate; it depends on unresolved operator
  and IDA contracts.
- `HP-OBJ-03`: keep as rejected.
- `HP-TEST-01`: keep as rejected.

## Missing Design IDs Or Overbroad IDs

Missing or underspecified production surfaces:

- A projection operation that explicitly projects a PQS shell candidate against
  accumulated previous terminal blocks.
- A support-index to support-state derivation/check helper, unless folded into
  `CartesianTerminalBasisBlock` construction.
- A block-pair workspace policy for IDA assembly, including when tiling is
  required.
- The exact owner file for `HP-FN-03` and `HP-FN-04`.

Overbroad IDs:

- `HP-FN-04` is too broad as a single approved helper. It covers positive-gauge
  weights, raw pair-factor construction, final contraction, sign handling, and
  symmetry placement. Consider splitting the design text into subrequirements
  without necessarily approving multiple helpers.
- `HP-FN-05` is too broad until the producer's inputs are fixed. It should not
  hide operator construction choices inside one wrapper.

## Slice Review

Slice A:
- Crosses a real numerical boundary only if it returns an actual terminal basis
  and restores H2 one-body physics.
- Must delete terminal source-realization preflight and blocker-vocabulary smoke
  assertions.
- Validation is sufficient if it includes package load, H one-body, H2
  one-body, final overlap, and Cr2 terminal basis/blocker behavior.
- The 150-line budget is aggressive and should be a target, not an excuse to
  skip projection correctness.

Slice B:
- Crosses a real numerical boundary by replacing final one-body assembly.
- Must delete or stop calling any dense support-matrix adapter it replaces.
- Validation should include H/H2 energy parity and a temporary old/new numeric
  comparison that is deleted before merge.
- Budget is plausible if the separable product helper is narrow.

Slice C:
- Crosses the main physics boundary by producing `CartesianIDAHamiltonian`.
- Must delete remaining blocked source-plan payload/report surfaces and smoke
  checks.
- Validation should include H2 IDA/self-Coulomb, artifact round trip, and Cr2
  memory report.
- Budget is plausible only after `HP-FN-04` is tightened.

Slice D:
- Crosses a usability boundary, not a new physics boundary.
- Must delete duplicate materialization/report graph rather than add a wrapper.
- Net-deletion requirement is appropriate.

## Questions For Manager Or ChatGPT-Pro

1. What is the exact mathematical IDA convention for `electron_electron_ida`
   after localized PQS column-sign canonicalization?
2. Should Slice A implement projection against accumulated previous terminal
   blocks immediately, or is a documented stop at
   `:terminal_pqs_cross_block_projection_required` acceptable for one local
   design iteration only?
3. Which existing source file owns blockwise one-body and IDA assembly without
   reviving the retired pair-materialization framework?
4. Is the 1.2 GiB peak-memory ceiling firm for Cr2, or should v2 specify a
   smaller target plus an absolute hard cap?

## Recommended Next Design Revision

Before ChatGPT-Pro review, revise the design to:

1. Demote positive HP-* items to candidates except the compact negative
   registry guards and the clearly bounded `HP-CHANGE-01`/`HP-FN-02`.
2. Add the accumulated-previous-block projection rule to Sections 4 and 7.
3. Define the IDA matrix convention precisely enough for an external consumer
   to build the same Hamiltonian.
4. Add an IDA block workspace/tiling rule tied to the Cr2 memory ceiling.
5. Fix `HP-FN-03` signature to include additive scaling and name the owner file
   for both one-body and IDA assembly.
6. Add the support-state derivation invariant to `HP-OBJ-01`.
7. Convert line budgets from hard caps to review targets wherever the hard cap
   could encourage omitting numerical validation.
