Pass 141 - fingerprint Be2 raw-box route producer

Purpose:

Add a focused test-only fingerprint of the existing private raw product-box
route producer for the Be2/PQS route. The goal is to prove exactly what raw
plan/coefficient objects are already constructible before we wrap any of them
as route-owned source-realization inputs.

Why now:

Pass 140 established compact support-window/order facts. The remaining source
realization blockers are:

```text
:raw_product_box_plan_objects
:pqs_axis_local_coefficients
:diatomic_complete_core_shell_source_plan_materializer
```

Pass 139 found that `_pqs_pqs_product_raw_box_route_producer` appears to return
raw product-box plans, raw PQS plans, a product unit, a descriptor, and pair
inventory, but it is currently private/shadow-only. We need a fingerprint before
promoting any of that path into route authority.

Task:

- Prefer a test-only change.
- Work in the focused Be2 fingerprint test or a new standalone focused test if
  that is cleaner:
  `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- First inspect the private producer signature in
  `src/CartesianContractedParentMetrics.jl`.
- Using the same Be2/PQS fixture ingredients as the focused route-driver test,
  call the existing private raw-box route producer if it is directly callable
  without production changes.
- Assert only compact facts:
  - status/object kind if present;
  - route shape or source keys;
  - raw product-box plan objects are present or exactly why unavailable;
  - raw PQS/source coefficient objects are present or exactly why unavailable;
  - product/doside unit object is present or exactly why unavailable;
  - pair inventory count/families if present;
  - private/shadow/probe-only flags remain true;
  - no source-plan/final-basis/H1/H1-J/Ham/public/export/artifact claims.

Decision rules:

- If the producer is not safely callable from a focused test without adding
  production glue, stop and report the signature/blocker. Do not add adapters.
- If calling it requires old WL materializer authority, stop and report.
- If it materializes large objects or takes more than 60 seconds, stop or report
  timing and whether precompilation dominated.
- Do not make this producer route authority in this pass.
- Do not wire it into `cartesian_assembly(...)`.

Trust boundary:

- Test/fingerprint only unless a tiny typo is exposed.
- No route-owned raw-box payload yet.
- No source-plan materializer.
- No final-basis/H1/H1-J/Ham materialization.
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

- Run the focused test you edited directly.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- If validation passes, commit with a clear message such as:
  `Fingerprint Be2 raw box producer`

Report back:

- Commit SHA if committed.
- Exact producer/helper called.
- Whether raw product-box plans, PQS coefficient objects, product unit, and pair
  inventory were available.
- Whether the producer is still private/shadow/probe-only.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
