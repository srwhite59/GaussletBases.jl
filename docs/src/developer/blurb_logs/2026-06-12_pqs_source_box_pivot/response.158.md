Pass 158 response - audit Be2/PQS Hamiltonian consumer contract

Role: repo-doer@macmini

Task type: no-edit audit. I edited no source, tests, fixtures, or downstream
files. I ran only read-only inspection commands plus `git status --short
--branch`.

Current Be2/PQS handoff inventory:

- Source-side framing remains source-box-first PQS. The handoff is private
  inspect-only and carries `source_plan`, `final_basis`, `h1_payload`, and
  `ham_input_payload` references.
- One-body: `one_body_hamiltonian` is the final H1 matrix reference from
  `h1_payload.final_hamiltonian.hamiltonian_matrix`; the focused Be2 test
  checks final dimension `221`, finite/symmetric H1, and lowest eigenvalue
  `-0.27746109235228694`.
- Two-body representation: `ham_input_payload` materializes
  `density_interaction` with representation
  `:pre_final_density_interaction`; the handoff exposes
  `pre_final_pair_matrix`, `final_to_pre_final_coefficients`,
  `pre_final_weights`, `support_weights`, `support_pair_raw_numerator`, and
  `raw_pair_factor_terms`.
- Conventions: density gauge is `:pre_final_localized_positive_weight`;
  raw pair-factor convention is `:raw_numerator`; support row order is
  `:core_then_shell`; metadata explicitly says source-box-first and that
  retained diagnostic weights are not IDA weights.
- Nuclear/electron/spin metadata: Be2 compact handoff carries nuclear charges
  `(4.0, 4.0)`, coordinates `((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0))`, nuclear
  repulsion `4.0`, electron count `8`, and spin sector
  `:closed_shell_singlet`.
- Explicit nonclaims already present: `private_inspect_only = true`,
  `dense_vee_materialized = false`, `h1_j_materialized = false`,
  `rhf_materialized = false`, `public_api = false`,
  `exports_materialized = false`, `artifacts_materialized = false`,
  `hamv6_materialized = false`, `cr2_ready = false`, and
  `hfdmrg_ready = false`.

Downstream expectations:

- HFDMRG public entrypoints expect `H, V` as `N x N` real/symmetric matrices
  plus orthonormal seed orbitals (`work/hfdmrg/src/HFDMRG.jl:13-24`,
  `work/hfdmrg/src/core.jl:773-810`). Its density-density backend wraps the
  supplied `V::AbstractMatrix` directly (`work/hfdmrg/src/backends/density_density.jl:1-30`)
  and the backend API documents density-density `Vee` as an `N x N` matrix
  with matched Fock/energy conventions (`work/hfdmrg/src/backend_api.jl:60-89`).
- HFDMRG also has sliced-integral backends, but those require either fixed
  `V6[a,b,c,d,n,m]` or ragged `Vblocks[n][m]` four-index slice blocks
  (`work/hfdmrg/src/backends/sliced_basis.jl:1-57`). The current PQS handoff
  does not provide that representation.
- CR2's current driver uses branch records with nuclear charges, `nup`, `ndn`,
  and `hf_mode`, constructs sorted one-body matrices from a basis bundle, builds
  H1 seed orbitals by diagonalizing `H`, and calls `HFDMRG.solve_hfdmrg` with
  sorted `H`, sorted `V`, and seed orbitals
  (`work/cr2/scripts/he2_diatomic_cartesian_cp_pec_driver.jl:156-233`,
  `447-455`, `582-623`, `1865-1900`).
- CR2's build path writes a `*_basis_bundle.jld2` with `include_ham = true`
  after building ordinary Cartesian Qiu-White operators and an interaction
  matrix (`work/cr2/scripts/he2_diatomic_cartesian_cp_pec_driver.jl:1644-1745`).
  Existing run summaries record `bundle_path`, `basis_dim`, H1 spectrum,
  branch `nup/ndn`, `hf_mode`, final energy, traces, and `hf_state_path`.

Gap analysis:

- The current handoff is enough for a private read-only inspector/comparison
  probe that checks dimensions, references, convention labels, nuclear/electron
  metadata, and nonclaim flags.
- It is not enough to mark HFDMRG ready. HFDMRG's density-density path expects
  a final-space `V::N x N` with solver/Fock conventions already reviewed. The
  PQS handoff currently exposes a pre-final density-interaction representation
  plus transform/provenance objects; that may be inspectable, but it is not yet
  a reviewed HFDMRG `V` contract.
- It is not enough for HFDMRG sliced backends or HamV6 export because no
  `V6`/`Vblocks`/dense two-electron integral object is materialized.
- It is not enough for CR2 production handoff because there is no agreed CR2
  bundle/JLD2/HamV6 artifact format, no seed-orbital contract, and no export
  behavior.

Recommended next implementation pass:

Add a private route-owned consumer contract payload, consuming the existing
`diatomic_complete_core_shell_hamiltonian_handoff_payload` without building new
physics. Suggested name:

`_PQSDiatomicCompleteCoreShellHamiltonianConsumerContractPayload`

Smallest fields:

- `status`
- `blocker`
- `route_family`
- `source_handoff`
- `source_handoff_status`
- `final_dimension`
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
- `readiness`
- `available_objects`
- `missing_objects`
- `summary`
- `metadata`

The `readiness` summary should say `private_inspector_ready = true` only when
the current handoff is available, while keeping
`hfdmrg_density_density_ready = false`, `hfdmrg_sliced_ready = false`,
`hamv6_export_ready = false`, `cr2_ready = false`, `public_api = false`, and
`artifacts_materialized = false`.

Proposed status/blocker vocabulary:

- Before this object exists, keep the current readiness blocker:
  `:missing_diatomic_hamiltonian_consumer_contract`.
- For the new private object when the current handoff is available:
  `status = :available_diatomic_complete_core_shell_hamiltonian_consumer_contract_payload`
  and `blocker = :missing_downstream_hamiltonian_consumer_contract` or another
  downstream-only blocker in its readiness summary, not in the inspect-only
  object availability.
- Downstream blockers should remain explicit:
  `:missing_hfdmrg_density_density_contract` for the final-space `H,V` solver
  contract; `:missing_hfdmrg_sliced_integrals` or
  `:missing_dense_two_body_integrals` for sliced/four-index integrals;
  `:missing_hamv6_export_contract` for HamV6; and
  `:missing_cr2_handoff_format` for CR2 bundle/artifact integration.

Deletion/shrinkage forecast:

- deleted: none in this audit
- simplified: none in this audit
- quarantined: no new quarantine; existing handoff remains private inspect-only
- not deleted because: audit only, and the current handoff/test still records
  the active Be2/PQS route-spine blocker
- exact remaining caller/blocker:
  `diatomic_complete_core_shell_ham_readiness_payload` still blocks on
  `:missing_diatomic_hamiltonian_consumer_contract`; after a private consumer
  contract object exists, the remaining blockers should move to downstream
  `H,V`, sliced-integral/HamV6, and CR2 handoff format contracts.

Validation:

```text
git status --short --branch
## main...origin/main
```

No Julia tests or downstream scripts were run, per audit instructions.

-- repo-doer@macmini
