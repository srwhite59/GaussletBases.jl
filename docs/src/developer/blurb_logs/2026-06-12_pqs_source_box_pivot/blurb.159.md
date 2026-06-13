Pass 159 - add private Be2/PQS Hamiltonian consumer contract

Role: repo-doer@macmini

Task type: narrow implementation with blocker replacement.

Purpose:

Replace the current Be2/PQS readiness blocker

```text
:missing_diatomic_hamiltonian_consumer_contract
```

with a real private inspect-only Hamiltonian consumer contract payload. This
must not make the route HFDMRG-ready, CR2-ready, HamV6-ready, export-ready, or
physics-endpoint-ready. It only makes the current handoff inspectable through an
explicit contract and moves the remaining blocker to downstream-specific
formats/contracts.

Hard bloat gate:

This pass is invalid if it merely adds another payload while the active
probe-enabled readiness still reports
`:missing_diatomic_hamiltonian_consumer_contract`. The new object must replace
that missing contract in readiness. Do not add a new test file. Update only the
existing compact Be2 fingerprint test.

Expiration condition:

The inspect-only consumer contract is temporary. When a real downstream
consumer is chosen, this object should be replaced by or collapsed into the
format-specific contract for HFDMRG density-density `H,V`, sliced integrals,
HamV6, or CR2 bundle handoff. Record that expiration condition in metadata or
summary, not as a long doc block.

Governing policy:

- `AGENTS.md`: every line has carrying cost; new code must replace a live
  blocker or enable a current workflow.
- `docs/src/developer/pqs_source_box_operator_framework.md`: source-box-first
  PQS remains the algorithmic frame; support-row contraction is oracle/debug.
- `docs/src/developer/pqs_source_box_fixture_policy.md`: compact Be2/PQS route
  facts are route-smoke/convention diagnostics, not physics endpoints.
- Pass 158 audit: HFDMRG expects reviewed final-space `H,V` or sliced integral
  contracts; CR2 expects an agreed bundle/export handoff. The current handoff
  is not those.

Exact source surfaces:

- `src/pqs_source_box_diatomic_complete_core_shell.jl`
  - add a compact private payload type near
    `_PQSDiatomicCompleteCoreShellHamiltonianHandoffPayload`;
  - add a builder near
    `_pqs_source_box_route_driver_diatomic_complete_core_shell_hamiltonian_handoff_payload`;
  - update
    `_pqs_source_box_route_driver_diatomic_complete_core_shell_ham_readiness_payload`
    to accept/use the new consumer contract payload.
- `src/pqs_source_box_route_driver_helpers.jl`
  - wire the new builder in `cartesian_assembly(...)` immediately after the
    Hamiltonian handoff payload and before readiness.
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  - update the existing compact test only.

Suggested object name:

```julia
_PQSDiatomicCompleteCoreShellHamiltonianConsumerContractPayload
```

Suggested builder name:

```julia
_pqs_source_box_route_driver_diatomic_complete_core_shell_hamiltonian_consumer_contract_payload
```

Minimum behavior:

When the existing Hamiltonian handoff is available, the new payload should be
available as a private inspect-only contract:

```text
status = :available_diatomic_complete_core_shell_hamiltonian_consumer_contract_payload
blocker = nothing
private_inspector_ready = true
```

It should carry compact references/statuses/summaries for:

- source handoff status;
- final dimension;
- one-body Hamiltonian status/reference;
- two-body representation kind/status, using the current pre-final density
  interaction representation;
- density gauge;
- raw pair-factor convention;
- support weight count;
- pre-final pair matrix shape;
- final-to-pre-final coefficient shape;
- nuclear charges/coordinates/repulsion;
- electron count and spin sector.

It must explicitly keep downstream readiness false:

```text
hfdmrg_density_density_ready = false
hfdmrg_sliced_ready = false
hamv6_export_ready = false
cr2_ready = false
public_api = false
exports_materialized = false
artifacts_materialized = false
```

Readiness replacement:

Update diatomic Ham readiness so that when the consumer contract payload is
available:

- `:diatomic_hamiltonian_consumer_contract` is in available objects;
- `:diatomic_hamiltonian_consumer_contract` is not in missing objects;
- readiness no longer blocks on
  `:missing_diatomic_hamiltonian_consumer_contract`;
- readiness blocks on a downstream-specific blocker, preferably:

```text
:missing_hfdmrg_density_density_contract
```

and missing objects may include:

```text
:hfdmrg_density_density_contract
:hfdmrg_sliced_integrals
:hamv6_export_contract
:cr2_handoff_format
```

Keep this compact. Do not add a long missing-object inventory.

Test updates:

In the existing Be2 compact test, add only the assertions needed to prove the
blocker replacement:

- `assembly.diatomic_complete_core_shell_hamiltonian_consumer_contract_payload`
  exists and is available;
- `private_inspector_ready` is true;
- final dimension is `221`;
- density gauge and raw pair convention match the handoff;
- downstream readiness flags are false;
- readiness blocker is now the downstream blocker, not
  `:missing_diatomic_hamiltonian_consumer_contract`.

Do not grow the test back into a field inventory.

Forbidden:

- No dense `Vee`.
- No final-space `V` construction.
- No `V6`/`Vblocks`.
- No H1/J materialization.
- No RHF/SCF/DIIS.
- No WL comparison.
- No CR2 or HFDMRG run.
- No public API.
- No export/HamV6/JLD2/artifact writing.
- No new test file.
- No large report-field cloud.

Validation:

Run:

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
```

Do not run full nested default suite unless a concrete ambiguity appears and
you explain why first.

Decision rules:

- If the new consumer contract cannot be wired without reintroducing a field
  cloud or growing the Be2 test wall, stop and report the blocker.
- If the readiness blocker remains
  `:missing_diatomic_hamiltonian_consumer_contract`, do not commit; the pass has
  not replaced the active blocker.
- If downstream inspection requires choosing a real export or solver format,
  stop at that design boundary.

Report back:

- new payload status/blocker vocabulary;
- old readiness blocker and new readiness blocker;
- net source/test line change if easy to report;
- validation commands/results;
- git status;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Commit if validation passes, with a message like:

```text
Add diatomic Hamiltonian consumer contract
```

-- repo-manager@macmini
