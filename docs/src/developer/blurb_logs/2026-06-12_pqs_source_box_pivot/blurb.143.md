Pass 143 - audit raw-box to source-plan realization mapping

Purpose:

Do a no-edit audit of how the private Be2/PQS raw-box route payload could be
converted into the source-plan consumer shape required by existing
complete-core/shell final-basis and H1 helpers.

Why now:

Pass 142 made the raw objects route-owned private candidate data:

```text
product/doside unit support count: 25
left raw PQS boundary selected count: 98
right raw PQS boundary selected count: 98
total: 221
descriptor retained dimension: 221
```

That strongly suggests a possible source realization:

```text
core/body sector: product unit
shell/source sector: pqs_left then pqs_right
```

But the existing retained route order is `(:pqs_left, :pqs_right, :product)`,
so the mapping needs an explicit ordering/permutation contract before any
materializer can safely return an object that looks like
`:pqs_multilayer_shell_source_plan`.

Audit surfaces:

- `src/pqs_source_box_route_driver_helpers.jl`
  - `diatomic_complete_core_shell_support_window_payload`
  - `diatomic_raw_box_route_payload`
  - `diatomic_complete_core_shell_source_plan_payload`
- `src/pqs_multilayer_shell_source_plan.jl`
- `src/pqs_multilayer_complete_core_shell_h1.jl`
- `src/pqs_multilayer_support_one_body.jl`
- `src/pqs_multilayer_support_density.jl`
- `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`
- `src/CartesianContractedParentMetrics.jl`
  - raw product-box/PQS plan structures and product unit shape
- Focused Be2 fingerprint test.

Questions to answer:

1. Can `producer.product_unit.support_states` and `support_indices` honestly
   serve as `core_support_states` and `core_support_indices`?
2. Can left/right raw PQS plan views provide the shell/source support states,
   support indices, and coefficient blocks required for
   `shell_support_states`, `shell_support_indices`, and
   `shell_final_coefficients`?
3. Does `shell_final_coefficients` need to be block diagonal over left/right
   PQS sectors, and what is the expected shape?
4. Where should the retained-order to support-order permutation live?
5. Does the existing `bundles` field expected by density-input helpers map to
   the parent axis bundle, raw PQS plan bundles, or a new diatomic bundle object?
6. Are `metrics`, `core_box`, `outer_box`, `bond_axis`, `metadata`, and
   `layer_count` meaningful for Be2/PQS, or does returning
   `object_kind = :pqs_multilayer_shell_source_plan` require an adapter/new
   object first?
7. What is the smallest implementation pass that can safely move from raw-box
   payload to a source-realization payload or source-plan materializer?

Deliverable:

Recommend one next coding pass:

- either a source-realization payload that still does not claim
  `:pqs_multilayer_shell_source_plan`;
- or a minimal materializer that returns the existing source-plan consumer
  shape, if and only if the audit finds the mapping is exact.

Include:

- proposed object/helper names;
- required input objects;
- output fields;
- blockers;
- focused validation;
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
  adapter paths to route authority.
- Do not reinterpret retained diagnostic/self-integral weights as
  IDA/quadrature weights.
- Do not request interactive command approval during unattended baton work. If
  approval would be required, write `.agent_handoffs/ATTENTION.md` and stop.

Validation:

- Read-only inspection should be enough.
- Use `rg`, `sed`, `git show`, and similar read-only commands.
- Do not run tests unless a short focused probe is necessary to answer the
  mapping question.
- If any command is expected to exceed 60 seconds, report why instead of
  running it by default.

Report back:

- Files/helpers inspected.
- Proposed mapping from raw-box payload to source-plan fields.
- Whether returning `:pqs_multilayer_shell_source_plan` is honest yet.
- Recommended next coding pass.
- Git status.
- Deletion/shrinkage forecast:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
