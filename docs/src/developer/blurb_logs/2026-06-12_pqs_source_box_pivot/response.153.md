Pass 153 response - add private diatomic Hamiltonian handoff payload

Implemented and committed the private factor-first Be2/PQS Hamiltonian handoff
payload.

Commit:

```text
6cd55566 Add diatomic Hamiltonian handoff payload
```

Files changed:

- `src/pqs_source_box_route_driver_helpers.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

Payload/helper names:

- `_PQSDiatomicCompleteCoreShellHamiltonianHandoffPayload`
- `_pqs_source_box_route_driver_diatomic_complete_core_shell_hamiltonian_handoff_payload`
- assembly field:
  `diatomic_complete_core_shell_hamiltonian_handoff_payload`

Default behavior:

- `diatomic_complete_core_shell_hamiltonian_handoff_payload.status == :blocked_diatomic_complete_core_shell_hamiltonian_handoff_payload`
- blocker: `:missing_diatomic_complete_core_shell_source_plan`
- readiness remains blocked upstream on
  `:missing_diatomic_complete_core_shell_source_plan_producer`

Probe-enabled behavior:

- `diatomic_complete_core_shell_hamiltonian_handoff_payload.status == :available_diatomic_complete_core_shell_hamiltonian_handoff_payload`
- blocker: `nothing`
- readiness blocker moved to
  `:missing_diatomic_hamiltonian_consumer_contract`

Observed compact handoff facts:

- final dimension: `221`
- pre-final dimension: `221`
- density gauge: `:pre_final_localized_positive_weight`
- raw-pair convention: `:raw_numerator`
- support weight count: `275`
- pre-final pair matrix shape: `(221, 221)`
- final-to-pre-final coefficient shape: `(221, 221)`
- center count: `2`
- nuclear charges: `(4.0, 4.0)`
- nuclear coordinates: `((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0))`
- nuclear repulsion status/value:
  `:available_diatomic_nuclear_repulsion`, `4.0`
- electron count status/value:
  `:available_diatomic_electron_count_convention`, `8`
- spin sector status/value:
  `:available_diatomic_spin_sector_convention`,
  `:closed_shell_singlet`

Payload carries references to:

- source plan
- final basis
- H1 payload
- Ham-input payload
- one-body Hamiltonian matrix from the H1 payload
- pre-final density interaction
- pre-final pair matrix
- final-to-pre-final coefficients
- pre-final weights
- support weights
- support raw-pair numerator
- raw pair factor terms
- center records / compact center metadata

Boundary confirmations:

- no dense `Vee` was built
- no public export or HamV6 adapter was added
- no CR2/HFDMRG readiness was claimed
- no RHF/SCF/DIIS behavior was added
- no H1-J scalar diagnostic or H1-J materialization was added
- no WL payload implementation was added
- no retained diagnostic/self-integral weights were reinterpreted as
  IDA/quadrature weights

Validation:

```text
git diff --check
```

Result: passed.

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Result: passed. After the final source edit, `GaussletBases` precompiled in
about 56.9 s, then printed `load ok`.

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```

Result: passed.

```text
Be2 PQS Ham payload readiness fingerprint: 225 passed in 14.2 s
Be2 PQS probe-enabled Ham readiness fingerprint: 350 passed in 43.2 s
```

Git status after commit:

```text
## main...origin/main [ahead 1]
```

Deletion/shrinkage report:

- deleted: none
- simplified: readiness now advances from Ham-input availability to a compact
  private Hamiltonian handoff payload and one sharper downstream consumer
  blocker
- quarantined: dense `Vee`, public export/HamV6, CR2/HFDMRG readiness,
  RHF/SCF, WL payload, and H1-J remain unmaterialized/nonclaimed
- not deleted because: existing report-facing fields and old complete
  core/shell diagnostic-route compatibility surfaces remain live until manager
  assigns a retirement pass
- exact remaining caller/blocker: `cartesian_assembly(...)` now exposes
  `diatomic_complete_core_shell_hamiltonian_handoff_payload`; the remaining
  blocker is `:missing_diatomic_hamiltonian_consumer_contract`

-- repo-doer@macmini
