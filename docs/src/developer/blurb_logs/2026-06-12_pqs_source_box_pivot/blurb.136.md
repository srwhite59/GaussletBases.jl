Pass 136 - fingerprint probe-enabled Be2 readiness

Purpose:

Add one focused test-only fingerprint that separates two Be2/PQS readiness
blockers:

```text
parent axis-bundle object availability
diatomic complete-core/shell source-plan producer availability
```

Why now:

Pass 135 added the private diatomic/PQS Ham readiness payload. In the default
Be2 fixture it reports both:

- missing `:parent_axis_bundle_object`
- missing `:diatomic_complete_core_shell_source_plan_producer`

Before implementing a diatomic source-plan producer, we need to know whether the
parent-axis-bundle blocker is already solvable by existing structured route
options such as `probe_parent_axis_construction = :auto`, or whether it is part
of the producer work.

Task:

- Work only in
  `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  unless the test exposes a tiny typo.
- Add a second Be2/PQS assembly fixture or helper variant that is identical to
  the existing focused Be2 Ham fingerprint except for the minimum existing route
  option needed to populate a parent axis-bundle object, likely:

```julia
probe_parent_axis_construction = :auto
```

- Call `cartesian_assembly(...)`.
- Inspect `assembly.diatomic_complete_core_shell_ham_readiness_payload`.

Expected acceptable outcomes:

1. Probe-enabled route has `parent_axis_bundle_object_available == true`.
   Then assert:
   - readiness remains
     `:blocked_diatomic_complete_core_shell_ham_readiness`;
   - blocker remains
     `:missing_diatomic_complete_core_shell_source_plan_producer`;
   - `:parent_axis_bundle_object in readiness.available_objects`;
   - `:parent_axis_bundle_object` is not in `readiness.missing_objects`;
   - current private Ham payload remains blocked as before.

2. Probe-enabled route still cannot expose a parent axis-bundle object.
   Then assert the exact current status/blocker and report the missing
   structured carry. Do not invent scalar aliases.

Decision rule:

- If enabling the probe changes route semantics in a way that makes the focused
  Be2 fixture no longer comparable, stop and report the observed changed facts
  instead of forcing assertions.
- Do not make probe-enabled behavior the default.
- Do not add production code unless there is a tiny typo in the readiness
  payload.

Trust boundary:

- Test/fingerprint only.
- No diatomic source-plan producer.
- No final-basis, H1, H1/J, density interaction, RHF, SCF, WL payload, exports,
  artifacts, hfdmrg, CR2 execution, public API, or fixture promotion.
- No scalar report-field clouds.
- Do not make shell/support-row contraction route authority.
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
- Whether probe-enabled Be2 readiness has a parent axis-bundle object.
- Readiness status/blocker/missing objects for both default and probe-enabled
  fixture paths.
- Confirmation that the private Ham payload remains blocked.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
