Pass 159 response - add private Be2/PQS Hamiltonian consumer contract

Role: repo-doer@macmini

Implemented and committed:

```text
10108c7d Add diatomic Hamiltonian consumer contract
```

Files changed:

- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

New private payload:

```julia
_PQSDiatomicCompleteCoreShellHamiltonianConsumerContractPayload
```

Builder:

```julia
_pqs_source_box_route_driver_diatomic_complete_core_shell_hamiltonian_consumer_contract_payload
```

Status/blocker vocabulary:

- Available inspect-only contract:
  `:available_diatomic_complete_core_shell_hamiltonian_consumer_contract_payload`
  with `blocker = nothing`.
- Blocked contract when the handoff is absent:
  `:blocked_diatomic_complete_core_shell_hamiltonian_consumer_contract_payload`
  with `:missing_diatomic_complete_core_shell_hamiltonian_handoff_payload`.
- Non-PQS route:
  `:not_applicable_diatomic_complete_core_shell_hamiltonian_consumer_contract_non_pqs_route`.

Readiness blocker replacement:

- Old active blocker:
  `:missing_diatomic_hamiltonian_consumer_contract`.
- New active blocker when the consumer contract is available:
  `:missing_hfdmrg_density_density_contract`.
- Readiness now puts `:diatomic_hamiltonian_consumer_contract` in available
  objects and keeps it out of missing objects on the active Be2 assembly path.
- Remaining compact missing downstream objects:
  `:hfdmrg_density_density_contract`, `:hfdmrg_sliced_integrals`,
  `:hamv6_export_contract`, `:cr2_handoff_format`.

Payload behavior:

- Carries the existing handoff reference, final dimension, one-body H status
  and reference, pre-final density-interaction representation/status, density
  gauge, raw pair-factor convention, support/pair shapes, nuclear metadata,
  electron count, and spin sector.
- Sets `private_inspector_ready = true` only for the available handoff case.
- Keeps downstream readiness false:
  `hfdmrg_density_density_ready`, `hfdmrg_sliced_ready`,
  `hamv6_export_ready`, `cr2_ready`, `public_api`,
  `exports_materialized`, and `artifacts_materialized`.
- Metadata records the expiration condition:
  `:replace_when_real_downstream_hamiltonian_consumer_contract_is_chosen`.

Net source/test line change:

```text
src/pqs_source_box_diatomic_complete_core_shell.jl | 255 +++++++++++++++++++++
src/pqs_source_box_route_driver_helpers.jl         |   8 +
test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl | 27 ++-
3 files changed, 288 insertions(+), 2 deletions(-)
```

Validation:

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
Test Summary:                                   | Pass  Total   Time
Be2 PQS probe-enabled Ham readiness fingerprint |   63     63  49.0s
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

- deleted: old active readiness pressure on
  `:missing_diatomic_hamiltonian_consumer_contract` for the active Be2 assembly
  path
- simplified: readiness now distinguishes private inspectability from
  downstream `H,V`/export readiness
- quarantined: the new contract is private inspect-only and temporary; HFDMRG,
  CR2, HamV6/export, artifacts, dense `Vee`, sliced integrals, H1/J, and
  RHF/SCF remain unavailable
- not deleted because: the existing handoff remains the route-owned source of
  H1, pre-final density interaction, convention labels, and Be2 metadata
- exact remaining caller/blocker:
  `diatomic_complete_core_shell_ham_readiness_payload` now blocks on
  `:missing_hfdmrg_density_density_contract`; real downstream work still needs
  a reviewed HFDMRG density-density `H,V` contract, sliced integral/HamV6
  contract, or CR2 handoff format

-- repo-doer@macmini
