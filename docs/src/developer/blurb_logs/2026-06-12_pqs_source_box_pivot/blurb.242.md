# Pass 242 blurb - shrink PQS source metadata acceptance scaffold

Role: repo-doer.

Read before editing:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_manager_running_log.md`
- `docs/src/developer/old_flat_cartesian_retirement_audit_2026-06-14.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.241.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.241.md`

Task type: cleanup/shrink after pass-241 line-budget exception.

Purpose:

Pay down the pass-241 `+278` line-budget exception with a careful shrink of a
known yellow-risk source metadata acceptance scaffold. This is not a physics
implementation pass.

Deletion candidate:

```text
file:
  test/nested/pqs_source_metadata_real_artifact_acceptance_support.jl

related test:
  test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl
```

Subagent dossier:

```text
files/functions/tests:
  support file is 528 lines.
  main callable surface:
    _be2_pqs_q5_source_metadata_acceptance(...)
    _be2_pqs_q5_source_metadata_export_tables(...)
  test caller:
    test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl

why obsolete:
  Do not delete the whole support file. Source-shell/source-mode metadata is
  still a live private acceptance lane.

  Obsolete pressure is the historical scaffolding around it:
    explicit export wrapper for a missing artifact;
    path/privacy assertions;
    fixed-column source-relation assertions;
    no-op operator/Hamiltonian/postprocess rows;
    exact route-shadow count clouds preserving migration-era inventory details.

live source callers:
  None for the support-file functions themselves.

replacement/current authority:
  src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl
    _pqs_current_route_source_shell_mode_inventory
    _pqs_source_metadata_export_contract
    _write_pqs_source_shells_table
    _write_pqs_source_modes_table

risk class:
  yellow

expected line savings:
  easy target: about 100-140 lines
  broader compact shrink: about 150-180 lines if safe
```

Required cleanup target:

Shrink, do not delete wholesale.

Prioritize:

1. Delete `_be2_pqs_q5_source_metadata_export_tables(...)` if it is only
   testing an explicit export wrapper for a missing artifact.
2. Remove the explicit-export wrapper test block in
   `pqs_source_metadata_real_artifact_acceptance_runtests.jl`.
3. Simplify `_pqs_source_metadata_acceptance_table_records(...)` to the
   IOBuffer-only path if the path-writing branch becomes unused.
4. Drop `path` keyword plumbing in `_be2_pqs_q5_source_metadata_acceptance(...)`
   if it becomes unused.
5. Remove fixed-column source-relation setup/checks and no-op row assertions if
   they are route-shadow metadata pressure rather than live source metadata
   contract.
6. Optionally collapse exact category-count/no-go row clouds to compact
   live-contract assertions.

Preserve:

- source-shell/source-mode metadata acceptance;
- local axes are source-shell-local labels;
- parent-lattice axes are explicit separate fields;
- no inference of ray IDs, shell/ray/radial labels, relation weights, or
  source-to-final maps from this metadata;
- fake-PQS guardrails and independent H2 PQS work.

Forbidden:

- do not delete the whole support file;
- do not touch H2 independent source-plan/final-basis/H1 work;
- do not touch fake-PQS endpoint coverage;
- do not edit source code unless a trivial test-support helper removal requires
  no source behavior change;
- do not run CR2, Hamiltonian assembly, H1, H1-J, RHF, or full Be2 artifact
  acceptance.

Line budget:

- This cleanup pass must be net-negative in `src + test + bin`.
- Target at least `100` net deleted lines if the dossier is accurate.
- If you cannot safely delete at least `80` net lines, write `ATTENTION.md`
  explaining the blocker and exact remaining candidates.

Minimum validation:

Default focused validation:

```text
julia --project=. test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl
git diff --check
```

Run without `BE2_PQS_Q5_ARTIFACT_DIR`; it should load and skip the real-artifact
branch cleanly if no artifact is configured.

Validation judgment:

If you delete or substantially shrink the old scaffold that this test mostly
protected, do not spend minutes forcing the obsolete test shape to keep passing.
The validation question is whether active authority still loads and whether the
remaining live source-shell/source-mode metadata contract is covered. In that
case, run package load plus the smallest surviving metadata smoke/check you
leave behind, record why the deleted branch was stale, and do not revive old
export-wrapper or route-shadow assertions just to satisfy the former test.

Do not run:

- `test/nested/integration_runtests.jl`;
- `test/nested/cartesian_pair_stage_low_order_policy_runtests.jl`;
- slow old-flat or metadata-scaffold suites whose only purpose is to protect the
  branch being deleted;
- full Be2/CR2/Hamiltonian/HF/RHF routes.

Report:

- exact lines added/deleted;
- functions/tests deleted or simplified;
- live metadata contract preserved;
- validation and timings;
- deleted/simplified/quarantined/not-deleted accounting;
- exact remaining caller/blocker.

-- repo-manager@macmini
