Pass 263 - retire current-route metadata export stack if caller audit permits

Context:
- Current HEAD should include
  `275bd81b Shrink terminal assembly low-order flat mirrors`.
- Pass 261 added a valuable support-partition payload with a `+435` scoped
  implementation exception. Pass 262 paid down `-45`; more cleanup is still
  needed.
- A deletion audit identified the old current-route metadata export stack as
  the largest remaining mature candidate:
  - `src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl`
    is about 5.5k lines;
  - the main tracked test pressure appears to be
    `test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl` and
    `test/nested/pqs_source_metadata_real_artifact_acceptance_support.jl`;
  - `test/nested/integration_runtests.jl` includes that old acceptance test.
- This stack is old flat Cartesian/PQS migration scaffolding around
  `_pqs_current_route_*` inventories, source-shell/source-mode metadata, fixed
  column/source relations, safe-term authority comparisons, and source-metadata
  export tables. It is not the independent H2 PQS route authority and should
  not become public architecture by inertia.

Task:
Do a caller-audited retirement pass for the current-route metadata export stack.

Decision rule:
1. First audit callers of:
   - `src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl`;
   - `_pqs_current_route_*`;
   - `_be2_pqs_q5_source_metadata_*`;
   - `pqs_source_metadata_real_artifact_acceptance_*`.
2. If the only live non-doc callers are the old source-metadata acceptance
   support/test and the slow integration include, retire the stack:
   - delete `src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl`;
   - remove its include from `src/CartesianContractedParentMetrics.jl`;
   - delete
     `test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl`;
   - delete
     `test/nested/pqs_source_metadata_real_artifact_acceptance_support.jl`;
   - remove the acceptance test include from
     `test/nested/integration_runtests.jl`.
3. If any active source path or non-old-flat test still depends on these
   helpers, do not add adapters. Stop with a deletion card listing exact
   blockers and the smallest safe subset to delete next.

Strict exclusions:
- Do not touch independent H2 PQS source-plan/final-basis/H1/H1-J/RHF/
  support-partition/provider-block code.
- Do not touch provider blocks, mixed/GTO matrices, residual MWG, supplemented
  values, CR2/export, HamV6, or public API.
- Do not revive old source-metadata acceptance via a compatibility wrapper.
- Do not run the slow integration runner as validation.
- Do not chase every historical doc/log reference. If deletion succeeds, add at
  most a short note to an existing retirement/audit doc if needed; old pass
  logs may continue to mention retired helpers as history.

Validation:
- `git diff --check`.
- Package load:
  `julia --project=. -e 'using GaussletBases; println("load ok")'`
- Caller proof:
  run `rg` checks showing no remaining non-doc source/test callers of deleted
  `_pqs_current_route_*` and source-metadata acceptance helpers.
- Parse the edited integration runner if touched; do not execute the slow
  integration runner.

Line budget:
- If deletion proceeds, scoped `src + test + bin` should be strongly
  net-negative, likely several thousand lines.
- If deletion is blocked by live callers, scoped impact may be `0`, but the
  response must include deletion cards with exact blockers.

Report:
- Caller audit results, separated into live source callers, tests-only callers,
  docs/log historical references, and deleted callers.
- Files deleted/edited.
- Scoped line count for `src + test + bin`.
- Validation commands.
- Deletion/shrinkage result:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
