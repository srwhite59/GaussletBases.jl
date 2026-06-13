Pass 154 - refactor diatomic route payloads out of helper file

Purpose:

Do a behavior-preserving file split. The recent Be2/PQS route objects are
conceptually useful, but `src/pqs_source_box_route_driver_helpers.jl` has become
too large and now hides the subsystem structure.

This pass should move code, not change behavior.

Target:

Start extracting the private diatomic complete-core/shell route payloads from:

```text
src/pqs_source_box_route_driver_helpers.jl
```

into a new private include file, suggested:

```text
src/pqs_source_box_diatomic_complete_core_shell.jl
```

Keep `cartesian_assembly(...)` and broad route-driver glue in
`pqs_source_box_route_driver_helpers.jl`.

Move only coherent private diatomic definitions/helpers that are already used by
`cartesian_assembly`, for example:

- `_PQSDiatomicCompleteCoreShellSupportWindowPayload`
- `_PQSDiatomicRawBoxRoutePayload`
- `_PQSDiatomicCompleteCoreShellSourceRealizationPayload`
- `_PQSDiatomicCompleteCoreShellSourcePlan`
- `_PQSDiatomicCompleteCoreShellSourcePlanPayload`
- `_PQSDiatomicCompleteCoreShellFinalBasisPayload`
- `_PQSDiatomicCompleteCoreShellH1Payload`
- `_PQSDiatomicCompleteCoreShellHamInputPayload`
- `_PQSDiatomicCompleteCoreShellHamiltonianHandoffPayload`
- their directly supporting private helpers

Decision rules:

- This is a move-only/refactor-only pass.
- Preserve function names, struct names, statuses, blockers, summaries, and
  behavior.
- Do not change test assertions except for any unavoidable include/load fallout.
- Do not add new payloads.
- Do not add H1-J, dense `Vee`, export, CR2/HFDMRG readiness, RHF, or WL work.
- Do not split the Be2 test file yet unless the code move makes a tiny import
  fix unavoidable.
- If a helper is shared by one-center and diatomic paths, leave it in the
  original helper file unless moving it is obviously behavior-preserving and
  reduces dependency friction.
- If dependency order becomes ambiguous, stop and report rather than doing a
  broad include reshuffle.

Include wiring:

- Add the new file to the existing source include chain in the smallest local
  place that preserves current method availability.
- Prefer one new include, not a module boundary change.
- Keep all code private/internal; no exports.

Validation:

- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- If validation passes, commit with a clear message such as:
  `Split diatomic PQS route helpers`
- If the focused test takes more than 60 seconds, report elapsed time and
  whether precompilation dominated.

Report back:

- Commit SHA if committed.
- New file path(s).
- What moved and what intentionally stayed.
- Confirmation behavior/statuses did not change.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
