Pass 151 response - add private diatomic Ham-input payload

Implemented and committed the private Be2/PQS Ham-input payload boundary.

Commit:

```text
8c2186ef Add diatomic Ham input payload
```

Files changed:

- `src/pqs_source_box_route_driver_helpers.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

Payload/helper names:

- `_PQSDiatomicCompleteCoreShellHamInputPayload`
- `_pqs_source_box_route_driver_diatomic_complete_core_shell_ham_input_payload`
- assembly field: `diatomic_complete_core_shell_ham_input_payload`

Payload fields:

```text
status
blocker
route_family
source_plan
source_plan_status
final_basis
final_basis_status
h1_payload
h1_payload_status
density_provenance
density_provenance_status
support_weights
support_weights_status
raw_pair_factor_terms
raw_pair_factor_status
support_pair_raw_numerator
support_pair_raw_numerator_status
density_interaction
density_interaction_status
available_objects
missing_objects
summary
metadata
```

Default behavior:

- `diatomic_complete_core_shell_ham_input_payload.status == :blocked_diatomic_complete_core_shell_ham_input_payload`
- blocker: `:missing_diatomic_complete_core_shell_source_plan`
- readiness remains blocked at the upstream source-plan producer:
  `:missing_diatomic_complete_core_shell_source_plan_producer`

Probe-enabled behavior:

- `diatomic_complete_core_shell_ham_input_payload.status == :available_diatomic_complete_core_shell_ham_input_payload`
- blocker: `nothing`
- readiness blocker moved to `:missing_diatomic_ham_consumer_contract`
- readiness no longer reports `:missing_diatomic_complete_core_shell_h1_j_consumer`
  once the Ham-input payload is available.

Observed compact Ham-input facts from the focused test:

- density gauge: `:pre_final_localized_positive_weight`
- raw-pair convention: `:raw_numerator`
- final dimension: `221`
- pre-final dimension: `221`
- support weight count: `275`
- support raw-pair shape: `(275, 275)`
- pre-final pair matrix shape: `(221, 221)`

Boundary confirmations:

- no old one-center H1-J payload call was added
- no full Ham payload/export/handoff object was added
- no RHF/SCF/Fock behavior was added
- no public API, exports, artifacts, WL payload, CR2 execution, or hfdmrg
  execution was added
- no density inputs were synthesized from retained diagnostic/self-integral
  weights
- no old `:pqs_multilayer_shell_source_plan` density helper guards were
  broadened

Validation:

```text
git diff --check
```

Result: passed.

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Result: passed. `GaussletBases` precompiled in about 56.8 s, then printed
`load ok`.

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```

Result: passed.

```text
Be2 PQS Ham payload readiness fingerprint: 202 passed in 13.5 s
Be2 PQS probe-enabled Ham readiness fingerprint: 296 passed in 40.4 s
```

Git status after commit:

```text
## main...origin/main [ahead 1]
```

Deletion/shrinkage report:

- deleted: none
- simplified: readiness now depends on the compact private Ham-input payload
  instead of treating H1 availability as directly blocked on H1-J consumer work
- quarantined: H1-J/full-Ham/RHF/public/export/artifact behavior remains
  private/unmaterialized or absent
- not deleted because: existing report-facing fields and old complete
  core/shell diagnostic route aliases are still compatibility surfaces
- exact remaining caller/blocker: next blocked seam is
  `:missing_diatomic_ham_consumer_contract`

-- repo-doer@macmini
