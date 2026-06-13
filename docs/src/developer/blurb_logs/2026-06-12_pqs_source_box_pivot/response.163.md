Pass 163 response - add CR2 read-only Hamiltonian inspection view

Role: repo-doer@macmini

Implemented and committed:

```text
51aace02 Add CR2 read-only Hamiltonian inspection view
```

Where the CR2 read-only view lives:

- The view lives on the existing consumer contract readiness tuple:
  `consumer_contract_payload.readiness`.
- No new payload type, struct, test file, CR2/HFDMRG run, export, artifact, or
  public API was added.

Status/blocker labels added:

- `cr2_read_only_inspector_ready = true`
- `cr2_solver_ready = false`
- `cr2_export_ready = false`
- `cr2_handoff_blocker = :missing_cr2_solver_handoff_format`
- `two_body_representation_kind = :pre_final_density_interaction`
- `density_gauge = :pre_final_localized_positive_weight`
- `raw_pair_factor_convention = :raw_numerator`

Confirmation:

- CR2 solver/export readiness remains false.
- HFDMRG density-density and sliced readiness remain false.
- HamV6/export/artifact/public readiness remains false.
- Overall readiness blocker remains
  `:missing_hfdmrg_density_density_contract`.

Line-budget gate:

```text
git diff --numstat -- src test
15	27	src/pqs_source_box_diatomic_complete_core_shell.jl
8	4	test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```

Totals:

- src/test added: 23
- src/test deleted: 31
- net: -8

Validation:

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
Test Summary:                                   | Pass  Total   Time
Be2 PQS probe-enabled Ham readiness fingerprint |   58     58  49.4s
```

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
load ok
```

```text
git diff --check
```

passed with no output.

Git status after commit, before writing this response:

```text
## main...origin/main [ahead 1]
```

Deletion/shrinkage report:

- deleted: remaining uncalled duplicate handoff nonclaim/materialization fields
  in handoff conventions, summary, and metadata; redundant test assertions for
  materialized/private-inspect flags already covered by statuses and the
  consumer readiness view
- simplified: CR2 read-only inspectability is represented on the existing
  consumer readiness tuple instead of a new payload layer
- quarantined: CR2 solver/export, HFDMRG density-density/sliced, HamV6/export,
  artifacts, dense `Vee`, final-space `V`, `V6`/`Vblocks`, H1/J, RHF/SCF, and
  public API remain unavailable
- not deleted because: handoff still owns real H1, pre-final density
  interaction, convention labels, and Be2 nuclear/electron/spin facts
- exact remaining caller/blocker:
  `diatomic_complete_core_shell_ham_readiness_payload` still blocks on
  `:missing_hfdmrg_density_density_contract`

-- repo-doer@macmini
