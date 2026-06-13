Pass 147 - add private diatomic final-basis payload

Purpose:

Use the new private Be2/PQS source-plan object to materialize only the private
diatomic complete-core/shell final basis. This should move the probe-enabled
Be2/PQS route from:

```text
source plan available
blocker = :missing_diatomic_complete_core_shell_final_basis_consumer
```

to:

```text
diatomic final basis available
next blocker = missing diatomic H1/Ham consumer
```

Physics target:

- Medium-term target: make the PQS driver usable enough for downstream Be2 WL
  versus PQS Hamiltonian comparisons.
- This pass is still a private route-seam pass, not a physics endpoint.

Implementation target:

- Work in `src/pqs_source_box_route_driver_helpers.jl`.
- Add one private final-basis payload and helper. Suggested names:
  - `_PQSDiatomicCompleteCoreShellFinalBasisPayload`
  - `_pqs_source_box_route_driver_diatomic_complete_core_shell_final_basis_payload`
- Consume:
  - support-window payload
  - raw-box route payload
  - source-realization payload
  - private diatomic source-plan payload / source-plan object
- Wire the payload into `cartesian_assembly(...)` as one private field, for
  example `diatomic_complete_core_shell_final_basis_payload`.

Implementation shape:

- Do not call `pqs_multilayer_complete_core_shell_final_basis(plan)`, because it
  correctly hard-gates on the old one-center source-plan object kind.
- Reuse the lower complete-core/shell final-basis realization helper directly,
  using the same overlap-block construction strategy as the existing one-center
  wrapper:
  - core support indices/states from the private diatomic source plan
  - shell support indices/states from the private diatomic source plan
  - `shell_final_coefficients` from the private diatomic source plan
  - overlap blocks from source-plan metrics and support states
- Preserve the source-plan ordering contract:
  - support row order is product/core first, then left/right PQS shell rows
  - shell coefficient block structure is left/right block diagonal
  - retained/pre-final map is carried as metadata, not expanded into scalar
    report fields.

Expected probe-enabled behavior:

- `diatomic_complete_core_shell_final_basis_payload.status` is available.
- It carries a real `final_basis`.
- It reports compact facts such as:
  - source-plan status
  - precleanup retained dimension `221`
  - final dimension, observed from the built final basis
  - core support count `25`
  - shell support count `250`
  - shell retained count `196`
  - support row order
  - old source-plan object kind false
  - final basis materialized true
- It keeps:
  - H1 materialized false
  - H1-J materialized false
  - Ham materialized false
  - route/public/export/artifact flags false

Default Be2/PQS behavior:

- Remains blocked upstream on missing parent axis-bundle/raw-box/source-plan
  data, preserving the existing default blocked path.

Downstream blocker:

- Update the Be2 readiness/Ham fingerprint so the probe-enabled path no longer
  reports missing final-basis consumer if the final basis is available.
- The next blocker should become a sharper label such as:

```text
:missing_diatomic_complete_core_shell_h1_consumer
```

or:

```text
:missing_diatomic_complete_core_shell_ham_consumer
```

depending on the existing readiness helper vocabulary. Prefer H1-consumer if no
one-body support/H1 payload is built in this pass.

Decision rules:

- If the lower final-basis helper requires a convention/fact that the private
  diatomic source plan does not carry, stop and report the exact missing fact.
- If final-basis construction fails because of overlap rank/conditioning,
  report the status and blocker; do not paper over it by changing tolerances
  unless the existing helper has a clearly appropriate tolerance knob.
- Do not add old one-center fields to the diatomic source plan.
- Do not add scalar report-field clouds.

Test target:

- Update only:
  `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- Assert default blocked behavior remains behavior-preserving.
- Assert probe-enabled final-basis payload availability.
- Assert the final-basis payload references the private diatomic source-plan
  object.
- Assert compact dimensions/counts and nonclaim flags.
- Assert no H1/H1-J/Ham/public/export/artifact behavior is added.

Trust boundary:

- Private/internal only.
- No H1/H1-J/Ham materialization.
- No support-one-body or support-density consumers yet, except overlap blocks
  needed by the final-basis helper.
- No RHF/SCF/DIIS work.
- No WL payload, public API, exports, artifacts, hfdmrg, or CR2 execution.
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
- If the focused test takes more than 60 seconds, report elapsed time and
  whether precompilation dominated.

Report back:

- Commit SHA if committed.
- Exact payload/helper names.
- Default/probe-enabled final-basis status/blocker.
- Observed final dimension and relevant counts.
- Confirmation that no H1/H1-J/Ham/public/export/artifact behavior was added.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
