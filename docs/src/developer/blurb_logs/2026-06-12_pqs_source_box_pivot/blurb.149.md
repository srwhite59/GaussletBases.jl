Pass 149 - add private diatomic H1 payload

Purpose:

Implement the private Be2/PQS H1 seam identified in pass 148. This should move
the probe-enabled Be2/PQS route from:

```text
final basis available
blocker = :missing_diatomic_complete_core_shell_h1_consumer
```

to:

```text
private diatomic H1 available
next blocker = missing H1-J / electron-electron Ham consumer
```

Physics target:

- Medium-term target: make the PQS driver usable enough for downstream Be2 WL
  versus PQS Hamiltonian comparisons.
- This pass adds only a private one-electron Hamiltonian diagnostic seam. It is
  not a full Hamiltonian payload and not a physics endpoint.

Implementation target:

- Work in `src/pqs_source_box_route_driver_helpers.jl`.
- Add one private H1 payload and helper. Suggested names:
  - `_PQSDiatomicCompleteCoreShellH1Payload`
  - `_pqs_source_box_route_driver_diatomic_complete_core_shell_h1_payload`
- Wire it into `cartesian_assembly(...)` as one private field, for example:
  - `diatomic_complete_core_shell_h1_payload`
- Feed it from:
  - `diatomic_complete_core_shell_source_plan_payload`
  - `diatomic_complete_core_shell_final_basis_payload`
  - parent center records
  - axis layers from `parent.parent_axis_bundle_object`
  - Coulomb Gaussian expansion

Implementation shape:

- Do not call `pqs_multilayer_complete_core_shell_h1_payload(...)`, because it
  correctly hard-gates on the old one-center source-plan object kind.
- Keep old `pqs_multilayer_support_*` guards unchanged in this pass.
- In the new private helper, assert:
  - source plan object kind is `:pqs_diatomic_complete_core_shell_source_plan`
  - source plan status is available
  - final basis status is available
  - final basis support row order is `:core_then_shell`
  - diatomic support rows are core/product first, then left/right PQS shell rows
- Build support one-body data with private/local code or tiny private wrappers:
  - support kinetic from `source_plan.metrics` and
    `vcat(source_plan.core_support_states, source_plan.shell_support_states)`
  - support electron-nuclear by-center from center records, axis layers, Coulomb
    expansion coefficients, and existing lower centered-Gaussian factor
    primitives
- Transfer to the final basis using the existing lower final-basis transfer:
  - `CartesianFinalBasisRealization.pqs_complete_core_shell_final_one_body_matrix`
- Assemble and solve H1 using existing lower `CartesianFinalBasisRealization`
  helpers where possible.

Expected probe-enabled behavior:

- `diatomic_complete_core_shell_h1_payload.status` is available.
- It carries compact fields/statuses for:
  - source-plan status
  - final-basis status
  - support kinetic status
  - support electron-nuclear status
  - final kinetic status
  - final electron-nuclear status
  - final H1 status
  - final dimension `221`
  - lowest H1 energy, finite and reported as observed
- It should mark:
  - H1 materialized true
  - H1-J materialized false
  - density-density materialized false
  - Ham payload materialized false
  - RHF materialized false
  - public/export/artifact flags false

Default Be2/PQS behavior:

- Remains blocked upstream without parent-axis bundle/source plan/final basis.

Readiness/Ham fingerprint:

- Update the Be2 readiness/Ham fingerprint so the probe-enabled path no longer
  reports missing H1 consumer once H1 is available.
- The next blocker should become a sharper label such as:

```text
:missing_diatomic_complete_core_shell_h1_j_consumer
```

or:

```text
:missing_diatomic_complete_core_shell_electron_electron_consumer
```

Use the label that best matches the existing readiness vocabulary. Do not claim
a full Ham payload.

Decision rules:

- If centered Gaussian nuclear support construction needs a fact not present in
  parent/source-plan/final-basis payloads, stop and report the exact missing
  fact.
- If H1 assembly gives nonsymmetric/nonfinite matrices, return a blocked payload
  with the precise blocker rather than loosening checks.
- Do not broaden old one-center object-kind guards unless it is a tiny internal
  extraction with no behavior change for existing callers.
- Do not add scalar report-field clouds.

Test target:

- Update only:
  `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- Assert default blocked behavior remains behavior-preserving.
- Assert probe-enabled private H1 payload availability.
- Assert final dimension `221`.
- Assert lowest H1 energy is finite and report the observed value in the
  response.
- Assert no H1-J/Ham/RHF/public/export/artifact behavior is added.

Trust boundary:

- Private/internal only.
- H1 only; no H1-J, density-density, electron-electron Ham payload, or full Ham
  payload materialization.
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
- Default/probe-enabled H1 status/blocker.
- Observed final dimension and lowest H1 energy.
- Confirmation that no H1-J/Ham/RHF/public/export/artifact behavior was added.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
