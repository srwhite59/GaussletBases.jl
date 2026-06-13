Pass 155 - no-edit test-retirement audit for PQS/RHF/Be2 route work

Role: repo-doer@macmini

Task type: audit only. Do not edit files, do not delete tests, do not commit.

Governing policy:

- `AGENTS.md` says tests are code with runtime, maintenance, and conceptual
  cost.
- Development scaffolding tests should be deleted or quarantined once the
  transition they supported is complete.
- Stale helpers/tests should be deleted rather than preserved through adapters
  when a stronger downstream semantic or scientific check supersedes them.
- Do not add tests that mainly preserve old helper names, route-shadow
  vocabulary, all-pairs inventory details, transitional metadata flags, or
  repeated nonclaim flags.

Physics/workflow target:

Make the Be2/PQS Hamiltonian-constructor route usable enough for downstream CR2
inspection and eventual WL-vs-PQS comparison, while keeping the test suite
focused on numerical correctness and live route contracts. RHF remains a
private Hamiltonian validator only; HFDMRG is the serious HF/DMRG package.

Read first:

- `AGENTS.md`, especially the test scope/deletion policy
- `docs/src/developer/fixture_role_policy_2026-06-12.md`
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.154.md`
- `test/runtests.jl`
- the focused PQS/RHF/Be2 nested tests named below

Audit surfaces:

1. RHF private diagnostic tests:

```text
test/nested/pqs_multilayer_complete_core_shell_rhf_input_contract_runtests.jl
test/nested/pqs_multilayer_complete_core_shell_rhf_initial_density_runtests.jl
test/nested/pqs_multilayer_complete_core_shell_rhf_one_step_runtests.jl
test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl
```

2. Be2/PQS route fingerprint:

```text
test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```

3. Any default-runner includes in `test/runtests.jl` that preserve transitional
   payload vocabulary rather than a live numerical or route contract.

Audit questions:

For each candidate test or assertion group, decide:

- What non-obvious bug would this catch?
- Is that bug still live after the current Be2 source-plan/final-basis/H1/Ham
  handoff route exists?
- Is it superseded by a stronger downstream semantic or scientific check?
- Does it mainly preserve a private transitional payload field, old blocker,
  route-shadow vocabulary, or repeated nonclaim flag?
- Should the next pass delete it, keep it, or temporarily quarantine it as
  debug/manual?

Pay special attention to deletion, not just moving:

- RHF is not the product. If RHF tests remain, they should protect one compact
  Hamiltonian-validator convention, not four standing seam files by default.
- The Be2 fingerprint file should not permanently assert every construction
  payload field. The mature shape should be closer to:
  - one source-plan semantic smoke;
  - one final-basis semantic smoke;
  - one H1/Ham numerical check or handoff-readiness check.
- Repeated checks that public/export/artifact/RHF/WL/H1-J flags are false should
  be consolidated to one or two meaningful boundary tests, not repeated at every
  private seam.

Decision rules:

- If a test is only protecting an obsolete or superseded helper/blocker, mark it
  as `delete`.
- If a test is still useful only while one active transition is incomplete, mark
  it as `temporary_keep` with the exact blocker that should delete it later.
- If a test protects a live numerical convention or route contract that a
  downstream endpoint would not localize, mark it as `keep`.
- Do not propose adding new tests in this audit unless deleting multiple
  scaffolding tests requires one smaller replacement smoke.

Report back with:

- default-runner inventory for the RHF/PQS/Be2 tests inspected;
- deletion candidates, grouped by file/assertion region;
- tests or assertion groups to keep and why;
- any temporary keeps and exact deletion blockers;
- a proposed pass 156 deletion plan with expected net test/code shrinkage;
- whether any source helper deletion should be paired with test deletion;
- git status.

Validation:

Read-only audit only. Run:

```text
git status --short --branch
```

Do not run long tests for this audit unless you find an ambiguity that cannot be
resolved from source inspection. If any command is expected to take more than 60
seconds, explain why before running it.

Deletion/shrinkage report:

- deleted: none in this audit
- simplified: none in this audit
- quarantined: none in this audit
- not deleted because: pass 155 is inventory/decision only
- exact remaining caller/blocker: report the blockers for any tests you propose
  to keep temporarily

-- repo-manager@macmini
