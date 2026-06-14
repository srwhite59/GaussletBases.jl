# Pass 242 response - shrink PQS source metadata acceptance scaffold

Implemented the cleanup/shrink pass for the Be2 PQS source metadata acceptance
scaffold. This pass did not touch H2 independent PQS source-plan/final-basis/H1
work or fake-PQS endpoint coverage.

## Files changed

- `test/nested/pqs_source_metadata_real_artifact_acceptance_support.jl`
- `test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl`

Unrelated pre-existing dirty file left untouched:

- `AGENTS.md`

## Deleted / Simplified

Deleted:

- `_be2_pqs_q5_source_metadata_export_tables(...)`
- the explicit export-wrapper test block
- path/privacy assertions for missing-artifact export paths
- `source_shells_table_path` / `source_modes_table_path` plumbing
- path-writing branch in `_pqs_source_metadata_acceptance_table_records(...)`
- fixed-column source-relation inventory construction/checks from this
  source-metadata acceptance scaffold
- no-op rows for operator assembly, QW Hamiltonian construction, and full
  postprocess export

Simplified:

- `_pqs_source_metadata_acceptance_table_records(...)` now uses only the live
  IOBuffer writer path.
- source metadata no-go assertions now check the source-shell/source-mode
  inventory directly, without also carrying fixed-column relation diagnostics.

## Live Contract Preserved

Preserved source-shell/source-mode metadata acceptance for:

- source shell and source mode TSV presence;
- source shell and source mode headers matching the current export contract;
- source-local axis labels staying within each source shell's contracted
  dimensions;
- explicit parent-lattice axis fields/statuses;
- no ray/shell/radial labels in source metadata;
- no inference from centers, nearest grid, support order, support indices, or
  raw-to-final support;
- no retained-weight/IDA division, route-construction change, packet adoption,
  QW Hamiltonian change, or public-default consumption.

## Validation

Focused test, without `BE2_PQS_Q5_ARTIFACT_DIR`:

```text
julia --project=. test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl
```

Result:

```text
Test Summary:                                       | Broken  Total  Time
Be2 strict-PQS q5 source metadata acceptance opt-in |      1      1  0.0s
```

Whitespace:

```text
git diff --check
```

passed.

## Scoped Line Budget

`src + test + bin`:

```text
0	41	test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl
11	102	test/nested/pqs_source_metadata_real_artifact_acceptance_support.jl
```

Total:

```text
11 added / 143 deleted, net -132
```

## Git Status

```text
## main...origin/main
 M AGENTS.md
 M test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl
 M test/nested/pqs_source_metadata_real_artifact_acceptance_support.jl
?? docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/blurb.242.md
```

`AGENTS.md` was already dirty and was not edited in this pass.

## Deletion/Shrinkage Accounting

deleted:
- explicit missing-artifact export wrapper and test;
- fixed-column source-relation route-shadow checks in this acceptance scaffold;
- path-writing/privacy assertion scaffold.

simplified:
- table parsing now uses the in-memory writer path only;
- no-go diagnostics are source-metadata scoped instead of duplicated through
  fixed-column relation metadata.

quarantined:
- live source-shell/source-mode metadata acceptance remains opt-in and skipped
  when no real Be2 artifact directory is configured.

not deleted because:
- the support file still protects live source metadata export contracts;
- the real-artifact branch remains useful when `BE2_PQS_Q5_ARTIFACT_DIR` is
  explicitly configured;
- fake-PQS endpoint coverage and independent H2 PQS work were out of scope.

exact remaining caller/blocker:
- `test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl`
  remains the only caller of
  `_be2_pqs_q5_source_metadata_acceptance(...)`.
- Further shrinkage is blocked on deciding which exact source count/category
  assertions are still acceptance-level versus migration-era inventory detail.

-- repo-doer@macmini
