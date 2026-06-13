Pass 139 - audit diatomic source-realization inputs

Purpose:

Do a no-edit audit of the actual source-realization inputs needed for Be2/PQS.
The goal is to determine whether existing structured route data can produce the
support states, support indices, retained/source coefficients, and factor
provenance required by the source-plan consumer shape, or whether a new
diatomic source-realization object must be designed.

Why now:

Pass 138 added a blocked private source-plan payload. In the probe-enabled Be2
path the parent axis bundle is available, and the remaining blocker is:

```text
:missing_diatomic_complete_core_shell_source_realization_contract
```

We should not build a fake `:pqs_multilayer_shell_source_plan`. The next
implementation has to be grounded in the actual Be2/PQS route objects.

Audit surfaces:

- `src/CartesianContractedParentMetrics.jl`
  - `_pqs_pqs_product_source_box_route_skeleton`
  - source boxes, retained units, pair entries
  - raw product box probe helpers around
    `_pqs_explicit_core_spacing_route_raw_product_box_plan_probe`
- `src/pqs_source_box_route_driver_helpers.jl`
  - pass-138 source-plan payload
  - Be2 readiness/source-plan fingerprints
  - probe option handling for parent axis construction and raw product boxes
- `src/pqs_multilayer_shell_source_plan.jl`
  - required source-plan fields and how support states are represented
- `src/pqs_multilayer_support_one_body.jl`
- `src/pqs_multilayer_support_density.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- Existing route-driver report tests that mention Be2 raw product box/materializer
  blockers.

Questions to answer:

1. Do Be2 retained-unit records contain actual retained/source coefficient
   matrices, or only retained counts/ranges/rules?
2. Do Be2 source boxes contain enough interval/window information to derive
   support states and support indices against the parent axis bundle?
3. Does existing raw product-box probe machinery produce structured support
   states/indices or only debug/materializer probes?
4. Are source-box density-density helper records relevant to source realization,
   or only pair/operator diagnostics?
5. What would the Be2 support order be if we tried to satisfy the existing
   `core_then_shell` source-plan consumer shape?
6. Is the honest next coding pass:
   - a source-realization readiness payload;
   - a raw product-box probe fingerprint;
   - a small source-box support-window extractor;
   - or a real diatomic source-realization materializer?

Deliverable:

Recommend the next coding pass with:

- proposed object/helper names;
- exact input objects and fields;
- output fields or blockers;
- whether raw product-box probe data should remain probe-only or become a
  route-owned structured input;
- one focused validation command;
- deletion/shrinkage forecast.

Trust boundary:

- No source edits.
- No test edits.
- No commits except the handoff response.
- No source-plan/final-basis/H1/H1-J/Ham materialization.
- No RHF/SCF work.
- No WL payload, public API, exports, artifacts, hfdmrg, or CR2 execution.
- No scalar report-field clouds.
- Do not promote raw product-box probes, shell/support-row contraction, or WL
  adapter paths to route authority during this audit.
- Do not reinterpret retained diagnostic/self-integral weights as
  IDA/quadrature weights.
- Do not request interactive command approval during unattended baton work. If
  approval would be required, write `.agent_handoffs/ATTENTION.md` and stop.

Validation:

- Read-only inspection should be enough.
- Use `rg`, `sed`, `git show`, and similar read-only commands.
- Do not run tests unless a short focused probe is necessary to answer the
  contract question.
- If any command is expected to exceed 60 seconds, report why instead of
  running it by default.

Report back:

- Files/helpers inspected.
- Whether retained units carry coefficients or only metadata.
- Whether source boxes can produce support windows/states.
- Status of raw product-box probe machinery and whether it is suitable route
  authority.
- Recommended next coding pass.
- Git status.
- Deletion/shrinkage forecast:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
