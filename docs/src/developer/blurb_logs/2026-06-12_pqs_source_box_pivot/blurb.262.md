Pass 262 - retire old current-route metadata export stack

Context:
- Current HEAD should include
  `4aa03d8e Materialize independent H2 PQS support partition`.
- Pass 261 was an accepted implementation exception: scoped `src + test + bin`
  was `+435 / -0`. The next cleanup-capable pass should pay that down with a
  mature deletion candidate.
- Read-only deletion audit found the strongest mature candidate:
  `src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl`
  plus the old Be2 source-metadata acceptance test/support.
- The candidate appears to have no live source callers outside its own include.
  Known remaining non-doc caller pressure is:
  `test/nested/pqs_source_metadata_real_artifact_acceptance_support.jl`, which
  is included only from `test/nested/integration_runtests.jl`.

Deletion candidate card:
- files/functions/tests:
  - `src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl`
  - include line in `src/CartesianContractedParentMetrics.jl`
  - `test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl`
  - `test/nested/pqs_source_metadata_real_artifact_acceptance_support.jl`
  - include line in `test/nested/integration_runtests.jl`
  - private entry points such as:
    `_pqs_current_route_retained_unit_inventory`,
    `_pqs_current_route_source_shell_mode_inventory`,
    `_pqs_source_metadata_export_contract`,
    `_pqs_current_route_safe_term_matrices`,
    `_pqs_current_route_safe_term_authority_comparison`
- why obsolete/stale:
  - This is an old current-route/source-metadata export and oracle comparison
    layer from the flat route-shadow period.
  - The active independent H2 PQS line now uses route-owned target/source-plan,
    final-basis, H1/H1-J/RHF, supplement preflight, and support-partition
    artifacts, not this current-route metadata TSV export.
  - It preserves old product/doside/support-local/current-route vocabulary that
    should not become public architecture.
- expected line savings:
  - about 5.5k source lines plus about 490 test/support lines, before small
    include/doc cleanup.
- risk class:
  - yellow: large deletion, but the audited live callers look tests-only.
- replacement/current authority:
  - independent H2 PQS driver artifacts and source-plan/final-basis payloads;
  - compact source-plan/support-partition summaries;
  - lower-level CPB/source-box provider tests for active local contracts.

Task:
1. Re-run a local caller audit before deleting:
   - search for the private entry points listed above;
   - distinguish live `src/` callers from docs/test-only callers.
2. If there are no live source callers outside
   `current_route_metadata_export.jl`, delete the old stack:
   - remove the include from `src/CartesianContractedParentMetrics.jl`;
   - delete `src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl`;
   - delete the two `pqs_source_metadata_real_artifact_acceptance*` test files;
   - remove their include from `test/nested/integration_runtests.jl`.
3. Update the retirement ledger or nearby docs only as needed to avoid claiming
   this stack remains live. Do not rewrite broad historical docs in this pass.
4. If a real live source caller is found, do not add a compatibility bridge.
   Stop and report the exact caller/blocker instead.

Strict exclusions:
- Do not touch independent H2 PQS supplement provider-block implementation.
- Do not modify pass-261 support partition behavior except if package load
  reveals a direct dependency issue.
- Do not run or preserve the old source-metadata acceptance test; it is the
  stale test pressure being retired.
- Do not run broad `test/nested/integration_runtests.jl`,
  `test/nested/cartesian_pair_stage_low_order_policy_runtests.jl`, or other
  slow old-flat integration gates.
- Do not create new public API, CR2/export, HamV6, residual MWG, supplemented
  values, H1/H1-J/RHF, or route-global matrices.

Minimum validation:
- `git diff --check`.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`.
- A small active-contract test only if package load points to one. Prefer
  `test/nested/cartesian_contracted_parent_metric_packet_runtests.jl` if the
  include deletion could affect the metric core.
- Do not run the deleted old acceptance test.

Line budget:
- This should be strongly net-negative in scoped `src + test + bin`.
- Report exact `git diff --numstat -- src test bin`.

Report:
- Confirm caller audit result:
  - live source callers:
  - tests-only callers:
  - docs-only references:
- Files deleted/edited.
- Validation commands and timings if any test exceeds 60 seconds.
- Scoped line count for `src + test + bin`.
- Deletion/shrinkage result:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
