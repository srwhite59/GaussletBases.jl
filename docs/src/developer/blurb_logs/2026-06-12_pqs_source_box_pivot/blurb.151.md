Pass 151 - add private diatomic Ham-input payload

Purpose:

Add the private Be2/PQS electron-electron/Ham-input payload recommended by pass
150. This should move the route from:

```text
H1 available
blocker = :missing_diatomic_complete_core_shell_h1_j_consumer
```

to:

```text
H1 + private Ham-input/electron-electron density interaction available
next blocker = missing Ham handoff/consumer contract
```

This is not a full Hamiltonian export. It is the next private constructor seam
needed before CR2 can inspect a PQS Hamiltonian handoff.

Physics target:

- Medium-term target: make the PQS driver usable enough for downstream Be2 WL
  versus PQS Hamiltonian comparisons.
- This pass provides private two-electron input data and conventions. It is not
  a physics endpoint.

Implementation target:

- Work in `src/pqs_source_box_route_driver_helpers.jl`.
- Add one private payload and helper. Suggested names:
  - `_PQSDiatomicCompleteCoreShellHamInputPayload`
  - `_pqs_source_box_route_driver_diatomic_complete_core_shell_ham_input_payload`
- Wire it into `cartesian_assembly(...)` as one private field, for example:
  - `diatomic_complete_core_shell_ham_input_payload`
- Feed it from:
  - `diatomic_complete_core_shell_source_plan_payload`
  - `diatomic_complete_core_shell_final_basis_payload`
  - `diatomic_complete_core_shell_h1_payload`
  - source-plan bundles / parent axis bundle
  - Coulomb Gaussian expansion

Implementation shape:

- Do not call the old one-center H1-J payload.
- Do not broaden old `:pqs_multilayer_shell_source_plan` density helper guards.
- Reuse lower density/final-basis primitives where possible.
- If an old helper is only blocked by object-kind but its internals are
  structurally compatible, add a private diatomic wrapper rather than changing
  the old guard.

Build only:

- support weights from the source-plan bundles / axis-factor provenance;
- raw pair factor terms from the same provenance;
- support raw pair numerator matrix in the diatomic support row order;
- pre-final density interaction via
  `CartesianFinalBasisRealization.pqs_complete_core_shell_pre_final_density_interaction(...)`;
- compact ordering/convention metadata;
- references to the already materialized H1/final Hamiltonian from the private
  H1 payload.

Expected probe-enabled behavior:

- `diatomic_complete_core_shell_ham_input_payload.status ==
  :available_diatomic_complete_core_shell_ham_input_payload`
- It should report compact facts/statuses such as:
  - source-plan status
  - final-basis status
  - H1 payload status
  - density provenance status
  - support weights status
  - raw pair factor status
  - support raw pair numerator status
  - pre-final density interaction status
  - density gauge
  - raw pair factor convention
  - final dimension `221`
  - pre-final dimension, observed from the density interaction
  - support row order / source-plan support row order
- It should carry, by reference where appropriate:
  - source plan
  - final basis
  - H1 payload or final H1 matrix reference
  - density interaction object
- It should mark:
  - Ham-input materialized true
  - H1-J materialized false
  - full Ham payload materialized false
  - RHF materialized false
  - public/export/artifact flags false

Default Be2/PQS behavior:

- Remains blocked upstream without source plan/final basis/H1.

Readiness/Ham fingerprint:

- Update the Be2 readiness/Ham fingerprint so the probe-enabled path no longer
  reports missing H1-J as the primary next blocker if the Ham-input payload is
  available.
- The next blocker should become a sharper label such as:

```text
:missing_diatomic_ham_consumer_contract
```

or:

```text
:missing_diatomic_complete_core_shell_ham_handoff_consumer
```

Use the name that best fits the current code vocabulary.

Decision rules:

- If source-plan bundles do not expose the needed axis-weight/raw-pair
  provenance, block with `:missing_diatomic_support_density_provenance`.
- If support weights are missing, block with `:missing_diatomic_support_weights`.
- If raw pair factors are missing, block with
  `:missing_diatomic_raw_pair_factor_terms`.
- If pre-final density interaction construction fails, block with
  `:blocked_diatomic_pre_final_density_interaction` and preserve the precise
  missing reason.
- Do not synthesize density inputs from retained diagnostic/self-integral
  weights.
- Do not compute H1-J self-Coulomb in this pass.
- Do not add scalar report-field clouds.

Test target:

- Update only:
  `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- Assert default blocked behavior remains behavior-preserving.
- Assert probe-enabled private Ham-input payload availability.
- Assert final dimension `221`.
- Assert density gauge and raw-pair convention labels are present.
- Assert no H1-J/full-Ham/RHF/public/export/artifact behavior is added.

Trust boundary:

- Private/internal only.
- No H1-J scalar diagnostic.
- No full Ham payload/export/handoff object yet.
- No RHF/SCF/DIIS work.
- No WL payload, public API, exports, artifacts, hfdmrg execution, or CR2
  execution.
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
  `Add diatomic Ham input payload`
- If the focused test takes more than 60 seconds, report elapsed time and
  whether precompilation dominated.

Report back:

- Commit SHA if committed.
- Exact payload/helper names.
- Default/probe-enabled Ham-input status/blocker.
- Observed density gauge, raw-pair convention, final dimension, and pre-final
  dimension.
- Confirmation that no H1-J/full-Ham/RHF/public/export/artifact behavior was
  added.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
