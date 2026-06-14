# Pass 242 manager review

Decision: accepted.

Commit reviewed:

- pending commit: shrink PQS source metadata acceptance scaffold

Scope reviewed:

- `test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl`
- `test/nested/pqs_source_metadata_real_artifact_acceptance_support.jl`

Findings:

- No blocking findings.
- The diff deletes the explicit missing-artifact export wrapper, path-writing
  branch, path/privacy assertions, fixed-column source-relation inventory
  checks, and no-op operator/Hamiltonian/postprocess rows from the source
  metadata acceptance scaffold.
- The remaining acceptance check still covers the live source-shell/source-mode
  metadata contract: TSV presence/headers, source-local axis labels, explicit
  parent-lattice axis fields, unavailable ray/shell/radial labels, and no
  retained-weight/IDA division, route-construction change, packet adoption, QW
  Hamiltonian change, or public-default consumption.
- The response's dirty-status note was stale relative to the manager commit
  `c2801090`; the doer did not edit `AGENTS.md`, and the accepted code diff is
  limited to the two intended test files.

Validation accepted:

- Doer ran the focused opt-in test without `BE2_PQS_Q5_ARTIFACT_DIR`; it skipped
  cleanly with one broken/skipped test in `0.0s`.
- Doer ran `git diff --check`.
- Manager reviewed the diff and reran `git diff --check`.

Line budget:

- Scoped `src + test + bin`: `11` added / `143` deleted, net `-132`.
- This pays down part of the pass-241 line-budget exception without touching the
  independent H2 PQS implementation or fake-PQS guard coverage.

Remaining blocker / next:

- Further shrinkage in this file is blocked on deciding which exact source
  count/category assertions are live acceptance contract versus migration-era
  inventory detail.
- The loop can continue with deletion candidates from the new sidecar audits or
  return to independent H2 PQS source-plan assembly once cleanup pressure has
  had enough paydown.

-- repo-manager@macmini
