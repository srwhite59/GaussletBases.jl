Pass 160 response - no-edit audit for diatomic handoff/consumer/readiness duplication

Role: repo-doer@macmini

Task type: no-edit source shrink audit. I edited no source, tests, or docs and
ran no Julia commands.

Duplicated field/flag inventory:

- Payload type fields:
  - Handoff owns real route objects and metadata at
    `src/pqs_source_box_diatomic_complete_core_shell.jl:322-358`.
  - Consumer repeats many handoff-derived scalars at
    `src/pqs_source_box_diatomic_complete_core_shell.jl:360-386`:
    `final_dimension`, one-body reference/status, representation kind/status,
    density gauge, raw pair convention, support count, pre-final/final shapes,
    nuclear coordinates/repulsion, electron count, and spin sector.
- Handoff builder owns the actual references and first summaries:
  - one-body and two-body references: lines `3000-3048`;
  - center/nuclear/electron/spin facts: lines `3050-3058`;
  - scalar summaries copied from those references: lines `3082-3180`;
  - nonclaim false flags repeated in conventions/summary/metadata:
    lines `3131-3135`, `3185-3193`, `3207-3218`.
- Consumer builder repeats handoff summaries:
  - copied scalar extraction: lines `3304-3360`;
  - downstream false flags in `readiness`: lines `3364-3378`;
  - the same flags again in `summary`: lines `3380-3408`;
  - the same flags again in `metadata`: lines `3410-3428`;
  - repeated constructor field list: lines `3430-3456`.
- Readiness repeats downstream flags again at lines `3677-3685`, after already
  making the downstream blocker/missing-object decision at lines `3566-3644`.
- The compact test now repeats consumer-derived field assertions at
  `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl:194-210`.
  The essential blocker-replacement assertions are only lines `212-217`.

What should stay:

- Handoff fields that carry actual route-owned objects should stay for now:
  `one_body_hamiltonian`, `density_interaction`, `pre_final_pair_matrix`,
  `final_to_pre_final_coefficients`, `pre_final_weights`, support/raw numerator
  data, center metadata, and nuclear/electron/spin facts. These are the current
  private inspectable route payload.
- Handoff summary fields for final dimension, density gauge, raw pair
  convention, support count, pair shapes, and nuclear/electron/spin facts can
  stay as the primary compact inspection summary.
- Readiness should keep only route-stage statuses, available/missing objects,
  and the active blocker. It is the right owner for
  `:missing_hfdmrg_density_density_contract`.
- The consumer contract should keep `status`, `blocker`, `route_family`,
  `source_handoff`, `source_handoff_status`, `readiness`, `available_objects`,
  `missing_objects`, `summary`, and `metadata`.

What should be consolidated or deleted:

- Delete the consumer payload's copied scalar fields unless a downstream caller
  actually needs them. The consumer can point at the handoff and expose a
  compact readiness summary instead.
- Delete consumer-builder scalar-copy code at lines `3304-3360`; use the
  handoff summary directly where summary fields are still needed.
- Consolidate downstream false flags into one tiny private helper or one local
  named tuple, used by consumer/readiness summaries. Current repeated false
  fields across handoff, consumer, and readiness are transitional convenience,
  not separate contracts.
- Trim the test to prove only:
  consumer payload available, `private_inspector_ready`, source handoff identity,
  downstream flags false from one readiness object, and readiness blocker
  replacement. It should not assert copied handoff scalars through the consumer.

Recommended pass 161 shrink plan:

1. Thin `_PQSDiatomicCompleteCoreShellHamiltonianConsumerContractPayload` to
   source-handoff reference plus compact readiness/summary/metadata.
2. Replace the repeated downstream-ready false bundle with one private helper,
   for example `_pqs_source_box_route_driver_diatomic_downstream_ham_nonclaims()`,
   returning the false flags and downstream missing objects.
3. Update readiness to consume the helper and keep
   `:missing_hfdmrg_density_density_contract`.
4. Shrink the existing Be2 fingerprint assertions to blocker replacement and
   nonpromotion facts only.

Proposed source/test line budget for pass 161:

- Target: at least 70 deletions and no more than 30 additions across tracked
  `src/` + `test/`, net reduction of at least 40 lines.
- Better target if the helper stays tiny: 90-120 deletions, 15-25 additions,
  net reduction of 65+ lines.
- Stale surfaces that can pay for the edit:
  consumer struct copied fields at lines `366-380`;
  consumer scalar extraction at lines `3304-3360`;
  duplicated downstream false flags at lines `3364-3378`, `3380-3408`,
  `3410-3428`, and readiness lines `3677-3685`;
  test assertions at lines `197-210` that only preserve copied scalar fields.

Exact blockers if no shrink is safe:

- No current downstream caller requires the consumer's copied scalar fields.
  Tests alone are not enough to keep them.
- The only safe reason to avoid shrinking would be a manager decision that the
  consumer summary itself is the review surface for downstream format selection.
  I did not find that requirement in pass 160.

Validation:

```text
git status --short --branch
## main...origin/main
```

Recent log context:

```text
eaa97c90 Publish PQS blurb 160
10108c7d Add diatomic Hamiltonian consumer contract
41a5dd56 Publish PQS blurb 159
```

Deletion/shrinkage report:

- deleted: none in this audit
- simplified: none in this audit
- quarantined: none in this audit
- not deleted because: audit only
- exact remaining caller/blocker:
  active readiness remains blocked on
  `:missing_hfdmrg_density_density_contract`; next source/test pass should
  shrink the consumer/readiness duplication before adding any HFDMRG, CR2,
  HamV6, dense `Vee`, or export format work

-- repo-doer@macmini
