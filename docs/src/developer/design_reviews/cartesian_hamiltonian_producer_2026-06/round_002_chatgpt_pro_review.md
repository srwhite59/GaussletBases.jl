# Round 002 ChatGPT-Pro Milestone Review

Reviewer: ChatGPT-Pro milestone review, pasted by user

Date: 2026-06-20

## Verdict

Do not freeze v2 yet. Run one focused docs-only Design Pass 003 to resolve the
Slice A mathematics and object semantics, then freeze Slice A only. Operator and
IDA slices should remain candidates until the terminal basis exists and reports
support growth, cross-overlaps, and memory behavior.

## Main Findings

- The design-first process is working and caught a real mathematical issue
  before implementation.
- The current blocker is numerical: terminal PQS shell realization lacks real
  projection/overlap/Lowdin construction.
- A new PQS shell must be projected against all previously accepted blocks
  before shell-local Lowdin.
- Existing `_CartesianNestedProjectedQShellLayer3D` explicitly does not project
  against previously locked spans, so v2's `projection_basis` wording was not
  implementable.
- After projection, a shell's effective support can include previous terminal
  rows. `support_indices` must mean effective coefficient support, not original
  terminal region support.
- Column sign canonicalization belongs in Slice A terminal basis finalization,
  before one-body or IDA assembly.
- Direct identity sectors need overlap and positive-weight validation.
- `HP-CHANGE-01` is not sufficient because returning shell overlap from the
  existing shell plan does not handle effective support growth after
  previous-block projection.
- The two-field terminal result wrapper and transitional missing-shell-input
  blockers should be rejected or deferred; missing shell inputs become
  implementation defects once the realizer exists.
- Slice A should use target plus redesign threshold: target 150 added source
  lines, redesign threshold 225, and net source decrease after preflight
  deletion.
- Permit an ignored `tmp/work` numerical spike to measure raw/projected
  cross-overlaps, effective supports, ranks, coefficient memory, and H2/Cr2
  behavior before implementation.

## Design Pass 003 Required Revisions

1. Replace the fictitious descriptor `projection_basis` input with an explicit
   previous-block projection formula.
2. Define block `support_indices` as effective support after projection.
3. Move positive IDA sign canonicalization into terminal basis realization.
4. Require direct sectors to validate direct overlap and positive IDA weights.
5. Reject or defer `HP-CHANGE-01`.
6. Remove transitional missing-shell-input blockers from the final realizer
   contract and reconsider `HP-RES-01`.
7. Remove `parent_dims` from `HP-OBJ-02` unless an identified consumer requires
   it.
8. Approve or freeze only Slice A surfaces; keep `HP-FN-03`, `HP-FN-04`, and
   `HP-FN-05` as future candidates.
9. Replace conflicting micro-budgets with Slice A target/redesign threshold and
   net-deletion requirement.
10. Define the uncommitted numerical spike and its required measurements.

## Deferred Questions

- Whether the Slice A freeze should mark the remaining Slice A candidates as
  implementation authority now or after one more user/manager signoff.
- Whether `HP-FN-00` should be file-local or a cross-file helper after the spike
  clarifies support growth.
