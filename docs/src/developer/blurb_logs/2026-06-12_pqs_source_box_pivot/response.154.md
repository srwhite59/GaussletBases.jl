Pass 154 response - split diatomic PQS route helpers

Implemented and committed the behavior-preserving file split.

Commit:

```text
b9f5fdae Split diatomic PQS route helpers
```

New file:

```text
src/pqs_source_box_diatomic_complete_core_shell.jl
```

Include wiring:

- Added one include in `src/GaussletBases.jl`, immediately after
  `include("pqs_source_box_route_driver_helpers.jl")`.

Moved:

- `_PQSDiatomicCompleteCoreShellHamReadinessPayload`
- `_PQSDiatomicCompleteCoreShellSupportWindowPayload`
- `_PQSDiatomicRawBoxRoutePayload`
- `_PQSDiatomicCompleteCoreShellSourceRealizationPayload`
- `_PQSDiatomicCompleteCoreShellSourcePlan`
- `_PQSDiatomicCompleteCoreShellSourcePlanPayload`
- `_PQSDiatomicCompleteCoreShellFinalBasisPayload`
- `_PQSDiatomicCompleteCoreShellH1Payload`
- `_PQSDiatomicCompleteCoreShellHamInputPayload`
- `_PQSDiatomicCompleteCoreShellHamiltonianHandoffPayload`
- their directly supporting private diatomic complete-core/shell helpers.

Intentionally stayed in `src/pqs_source_box_route_driver_helpers.jl`:

- `cartesian_assembly(...)`
- broad route-driver glue
- one-center/old complete core-shell diagnostic payload helpers
- shared descriptor/center/density helper surfaces used outside the diatomic
  complete-core/shell block

Behavior/status confirmation:

- This pass was a move-only refactor.
- Function names, struct names, statuses, blockers, summaries, and test
  assertions were not changed.
- No dense `Vee`, public export/HamV6, CR2/HFDMRG readiness, RHF/SCF, WL work,
  or H1-J materialization was added.

Validation:

```text
git diff --check
```

Result: passed.

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Result: passed. `GaussletBases` precompiled in about 56.6 s after the include
split, then printed `load ok`.

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```

Result: passed.

```text
Be2 PQS Ham payload readiness fingerprint: 225 passed in 14.1 s
Be2 PQS probe-enabled Ham readiness fingerprint: 350 passed in 43.0 s
```

Git status after commit:

```text
## main...origin/main [ahead 1]
```

Deletion/shrinkage report:

- deleted: no behavior deleted; 3,465 lines were removed from
  `pqs_source_box_route_driver_helpers.jl` by moving them to the new private
  file
- simplified: `pqs_source_box_route_driver_helpers.jl` now keeps
  `cartesian_assembly(...)` and broad route-driver glue instead of also owning
  the full diatomic complete-core/shell payload subsystem
- quarantined: diatomic payloads remain private/internal in the new file; dense
  `Vee`, export/HamV6, CR2/HFDMRG readiness, RHF/SCF, WL, and H1-J remain
  nonclaimed
- not deleted because: this pass was move-only and compatibility/report
  surfaces remain live
- exact remaining caller/blocker: `cartesian_assembly(...)` still calls the
  moved private helpers; route readiness remains blocked at
  `:missing_diatomic_hamiltonian_consumer_contract`

-- repo-doer@macmini
