Pass 263 response - retire current-route metadata export stack

Caller audit result:
- Live source callers:
  - none outside the candidate source file and its include.
  - proof command before deletion:
    `rg -n "_pqs_current_route_" src test bin --glob '!src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl' --glob '!test/nested/pqs_source_metadata_real_artifact_acceptance_support.jl' --glob '!test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl'`
    returned no matches.
- Tests-only callers:
  - `test/nested/pqs_source_metadata_real_artifact_acceptance_support.jl`
  - `test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl`
  - `test/nested/integration_runtests.jl` included the old acceptance test.
- Docs/log historical references:
  - many references remain in developer docs, baton logs, archived handoffs, and the current pass blurb.
  - I did not chase those historical references per the blurb.
- Deleted callers:
  - source include in `src/CartesianContractedParentMetrics.jl`;
  - slow integration include in `test/nested/integration_runtests.jl`;
  - the old source-metadata acceptance support/test files.

Files deleted/edited:
- deleted `src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl`
- edited `src/CartesianContractedParentMetrics.jl` to remove its include
- deleted `test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl`
- deleted `test/nested/pqs_source_metadata_real_artifact_acceptance_support.jl`
- edited `test/nested/integration_runtests.jl` to remove the acceptance test include

Validation:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed after precompile: `load ok`
- `julia --project=. -e 'Meta.parseall(read("test/nested/integration_runtests.jl", String)); println("integration runner parse ok")'`
  - passed: `integration runner parse ok`
- `git diff --check`
  - passed.
- Caller proof after deletion:
  - `rg -n "_pqs_current_route_|_be2_pqs_q5_source_metadata_|pqs_source_metadata_real_artifact_acceptance_" src test bin --glob '!*.tmp'`
    returned no matches.
  - `rg -n "current_route_metadata_export" src test bin --glob '!*.tmp'`
    returned no matches.

I did not run the slow integration runner.

Scoped line count for `src + test + bin`:
- `src/CartesianContractedParentMetrics.jl`: +0 / -1
- `src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl`: +0 / -5489
- `test/nested/integration_runtests.jl`: +0 / -1
- `test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl`: +0 / -53
- `test/nested/pqs_source_metadata_real_artifact_acceptance_support.jl`: +0 / -437
- total: +0 / -5981

Git status:
```text
## main...origin/main
 M src/CartesianContractedParentMetrics.jl
 D src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl
 M test/nested/integration_runtests.jl
 D test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl
 D test/nested/pqs_source_metadata_real_artifact_acceptance_support.jl
```

Deletion/shrinkage result:
- deleted: 5,981 scoped lines from the old current-route metadata export stack and its acceptance test pressure.
- simplified: `CartesianContractedParentMetrics` no longer includes the old `_pqs_current_route_*` metadata/export helper stack, and the slow integration runner no longer pulls in the retired source-metadata acceptance test.
- quarantined: historical docs/log references remain as history; no independent H2 PQS source-plan/final-basis/H1/H1-J/RHF/support-partition/provider-block code was touched.
- not deleted because: developer docs and baton/archive references were not chased by design; they are historical references outside the live source/test caller set.
- exact remaining caller/blocker: no remaining `src/test/bin` caller found for `_pqs_current_route_*`, `_be2_pqs_q5_source_metadata_*`, `pqs_source_metadata_real_artifact_acceptance_*`, or `current_route_metadata_export`.

-- repo-doer@macmini
