Pass 161 response - shrink diatomic Hamiltonian consumer duplication

Role: repo-doer@macmini

Implemented and committed:

```text
1274c356 Thin diatomic Hamiltonian consumer contract
```

Files changed:

- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

Fields removed from the consumer payload:

- `final_dimension`
- `one_body_hamiltonian`
- `one_body_hamiltonian_status`
- `two_body_representation_kind`
- `two_body_representation_status`
- `density_gauge`
- `raw_pair_factor_convention`
- `support_weight_count`
- `pre_final_pair_matrix_shape`
- `final_to_pre_final_coefficient_shape`
- `nuclear_charges`
- `nuclear_coordinates`
- `nuclear_repulsion`
- `electron_count`
- `spin_sector`

The consumer now keeps only the live contract fields: status/blocker,
route family, source handoff reference/status, readiness, available/missing
objects, summary, and metadata.

Downstream nonclaim consolidation:

- Added one private helper:
  `_pqs_source_box_route_driver_diatomic_downstream_ham_nonclaims()`.
- The helper centralizes false downstream readiness flags and the compact
  downstream missing-object tuple for HFDMRG density-density, sliced integrals,
  HamV6 export, and CR2 handoff format.
- Consumer and readiness summaries now use that bundle instead of maintaining
  separate repeated false flag lists.

Readiness behavior preserved:

- Consumer contract remains available on the Be2 probe-enabled path.
- `private_inspector_ready` remains true.
- `:diatomic_hamiltonian_consumer_contract` remains available and not missing.
- Readiness blocker remains `:missing_hfdmrg_density_density_contract`.
- No downstream readiness flag becomes true.

Old and new Be2 assertion counts:

- Before pass 161: 63 assertions.
- After pass 161: 55 assertions.

Line-budget gate:

```text
git diff --numstat -- src test
45	142	src/pqs_source_box_diatomic_complete_core_shell.jl
9	12	test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```

Totals:

- src/test added: 54
- src/test deleted: 154
- net: -100

Validation:

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
Test Summary:                                   | Pass  Total   Time
Be2 PQS probe-enabled Ham readiness fingerprint |   55     55  49.2s
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

- deleted: copied consumer scalar fields and scalar-copy builder code; repeated
  consumer/readiness downstream false-flag lists; test assertions that only
  preserved copied consumer fields
- simplified: consumer contract now points at the handoff as authority and
  carries a compact readiness/nonclaim summary
- quarantined: consumer remains private inspect-only; HFDMRG density-density,
  sliced integrals, HamV6/export, CR2 format, dense `Vee`, H1/J, and RHF/SCF
  remain unavailable
- not deleted because: the handoff still owns the real route objects and
  compact scalar inspection summary; readiness still owns the active downstream
  blocker
- exact remaining caller/blocker:
  `diatomic_complete_core_shell_ham_readiness_payload` still blocks on
  `:missing_hfdmrg_density_density_contract`

-- repo-doer@macmini
