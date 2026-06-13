Pass 153 - add private diatomic Hamiltonian handoff payload

Purpose:

Add the private factor-first Be2/PQS Hamiltonian handoff payload recommended by
pass 152. This should move the route from:

```text
Ham-input/electron-electron density interaction available
blocker = :missing_diatomic_ham_consumer_contract
```

to:

```text
private Hamiltonian handoff available
next blocker = missing downstream consumer/export/inspection contract
```

This is the object that should make the Be2/PQS constructor output inspectable.
It is not public export and not CR2/HFDMRG ready yet.

Implementation target:

- Work in `src/pqs_source_box_route_driver_helpers.jl`.
- Add one private payload and helper. Suggested names:
  - `_PQSDiatomicCompleteCoreShellHamiltonianHandoffPayload`
  - `_pqs_source_box_route_driver_diatomic_complete_core_shell_hamiltonian_handoff_payload`
- Wire it into `cartesian_assembly(...)` as one private field, for example:
  - `diatomic_complete_core_shell_hamiltonian_handoff_payload`
- Feed it from:
  - `diatomic_complete_core_shell_source_plan_payload`
  - `diatomic_complete_core_shell_final_basis_payload`
  - `diatomic_complete_core_shell_h1_payload`
  - `diatomic_complete_core_shell_ham_input_payload`
  - center records / parent center table if needed for nuclear metadata

Payload shape:

Carry references rather than copying large matrices:

- source plan
- final basis
- H1 payload
- Ham-input payload
- one-body Hamiltonian reference from H1 payload
- density interaction reference from Ham-input payload
- pre-final pair matrix reference from density interaction, if already carried
- final-to-pre-final coefficients reference, if already carried
- pre-final weights reference, if already carried
- support weights reference
- support raw-pair numerator reference
- raw pair factor terms reference
- center records or compact center metadata

Carry compact inspect metadata:

- status/blocker
- route family
- final dimension `221`
- pre-final dimension
- support weight count
- pre-final pair matrix shape
- final-to-pre-final coefficient shape
- density gauge
- raw-pair convention
- support/pre-final/final ordering labels
- center count
- nuclear charges
- nuclear coordinates
- nuclear repulsion status and value, if computable from center records
- electron count status/value
- spin sector status/value
- `private_inspect_only = true`
- public/export/HamV6/CR2/HFDMRG/RHF flags false

Nuclear/electron-count convention:

- If center records clearly give two Be centers with nuclear charges `(4, 4)`,
  record:
  - `electron_count = 8`
  - `electron_count_source = :neutral_sum_nuclear_charges_private_route_smoke`
  - `spin_sector = :closed_shell_singlet`
  - `spin_sector_source = :private_be2_route_smoke_default`
- If that cannot be established from the current parent/center data, keep the
  handoff blocked with:
  - `:missing_diatomic_electron_count_convention`
  or
  - `:missing_diatomic_spin_sector_convention`
- Do not make CR2 infer electron count.

Nuclear repulsion:

- If center locations and charges are available, compute the two-center nuclear
  repulsion as `Z1 * Z2 / norm(R1 - R2)` and label:
  - `nuclear_repulsion_source = :center_record_charge_distance`
- If units or coordinates are ambiguous, block with
  `:missing_diatomic_nuclear_repulsion`.

Expected probe-enabled behavior:

- `diatomic_complete_core_shell_hamiltonian_handoff_payload.status ==
  :available_diatomic_complete_core_shell_hamiltonian_handoff_payload`
  if all required inspect facts are available.
- It should report:
  - final dimension `221`
  - density gauge `:pre_final_localized_positive_weight`
  - raw-pair convention `:raw_numerator`
  - electron count `8` and spin sector `:closed_shell_singlet`, if established
  - nuclear repulsion value, if established
  - private inspect-only true
  - public/export/HamV6/CR2/HFDMRG/RHF false
- Readiness should move to:
  - `:missing_diatomic_hamiltonian_consumer_contract`
  or
  - `:missing_diatomic_cr2_inspection_contract`
  depending on the naming that fits the current helper.

Default behavior:

- Remains blocked upstream without source plan/final basis/H1/Ham-input.

Decision rules:

- Do not build dense `Vee`.
- Do not create a HamV6/export adapter.
- Do not run CR2 or hfdmrg.
- Do not add public API.
- Do not copy large matrices into summary metadata.
- If a required convention is missing, block with a precise blocker rather than
  guessing.
- Keep H1-J materialized false.

Test target:

- Update only:
  `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- Assert default blocked behavior remains behavior-preserving.
- Assert probe-enabled private handoff availability, or the precise blocker if
  electron count / nuclear repulsion cannot be established.
- Assert compact fields and nonclaim flags.
- Assert no dense `Vee`, public export, CR2/HFDMRG readiness, RHF, or H1-J
  behavior is added.

Trust boundary:

- Private/internal inspect-only handoff.
- No dense `Vee`.
- No public export/HamV6 adapter.
- No H1-J scalar diagnostic.
- No RHF/SCF/DIIS work.
- No WL payload implementation.
- No CR2 or hfdmrg execution.
- No edits outside this repo.
- Do not promote shell/support-row contraction, raw product-box probes, or old
  WL adapter paths to route authority.
- Do not reinterpret retained diagnostic/self-integral weights as
  IDA/quadrature weights.
- Do not request interactive command approval during unattended baton work. If
  approval would be required, write `.agent_handoffs/ATTENTION.md` and stop.

Validation:

- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- If validation passes, commit with a clear message such as:
  `Add diatomic Hamiltonian handoff payload`
- If the focused test takes more than 60 seconds, report elapsed time and
  whether precompilation dominated.

Report back:

- Commit SHA if committed.
- Exact payload/helper names.
- Default/probe-enabled handoff status/blocker.
- Observed final dimension, density gauge, raw-pair convention, electron count,
  spin sector, and nuclear repulsion status/value.
- Confirmation that no dense `Vee`, public export, CR2/HFDMRG readiness, RHF, or
  H1-J behavior was added.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
