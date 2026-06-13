Pass 194 - shrink H2 readiness report aliases and audit the source-plan producer blocker.

Purpose:

Pass 193 correctly added an honest H2 gausslet-only readiness artifact, but it
also added a fairly broad set of diatomic readiness scalar report aliases. Do
not build H2 H1/H1-J on top of that surface by inertia.

This pass should pay down that carrying cost and identify the exact next
producer seam for the live blocker:

```text
:missing_diatomic_complete_core_shell_source_plan_producer
```

Task type:

Refactor/shrink plus audit. No H2 materialization.

Current state:

- H2 driver input exists:
  `test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl`
- H2 readiness test exists:
  `test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
- Artifact currently records H2 geometry, `q = n_s = 5`,
  `supplement_policy = :none`, `comparison_ready = false`, and the live
  readiness blocker.
- Current blocker:
  `:missing_diatomic_complete_core_shell_source_plan_producer`

Read/inspect:

```text
bin/cartesian_ham_builder.jl
src/pqs_source_box_route_driver_helpers.jl
src/pqs_source_box_route_driver_reporting.jl
src/pqs_source_box_diatomic_complete_core_shell.jl
test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.193.md
```

Exact task:

1. Shrink the report alias surface from pass 193.

   Inspect `_pqs_source_box_route_driver_diatomic_complete_core_shell_report_fields`
   and the diatomic readiness writer. Prefer a compact report field such as:

   ```text
   diatomic_complete_core_shell_readiness_summary
   ```

   and have the artifact writer read route statuses/materialization booleans
   from that summary instead of requiring many duplicated scalar aliases.

   Keep only scalar report aliases that are genuinely needed by existing live
   callers. If only the new artifact writer consumes them, remove them and
   consume the compact summary instead.

   Preserve the on-disk artifact keys from pass 193. The artifact consumer
   should not need to know whether the report held compact summary data or many
   scalar aliases.

2. Audit the source-plan producer blocker.

   Answer:

   ```text
   What exact current object/status produces
   :missing_diatomic_complete_core_shell_source_plan_producer?

   Which existing private payloads already exist before this blocker?

   Which inputs are missing for a real diatomic source-plan producer?

   Is the blocker primarily:
     parent basis/axis bundle construction,
     route-configured diatomic materializer inputs,
     raw-box/source-realization payload wiring,
     fixed-q/source-mode policy,
     or a guard that rejects the H2 route_kind?
   ```

   Do not implement the producer in this pass. The goal is to make pass 195
   precise.

3. Keep the H2 readiness artifact behavior unchanged.

   Re-run the explicit readiness test and confirm it still passes.

Trust boundary:

- No H2 source-plan producer implementation.
- No H2 final basis or H1 materialization.
- No H1-J/density interaction.
- No private RHF for H2.
- No supplemented WL/QW comparison.
- No supplement support.
- No Be2/Cr2 artifact work.
- No HFDMRG, DMRG, ECP, exports, public solver behavior, or Qiu-White
  correction work.
- Do not add this test to default runners.
- Preserve the visible staged driver style.
- Do not request interactive approval or sandbox escalation. If approval is
  genuinely required, write `.agent_handoffs/ATTENTION.md` and stop.

Line-budget rule:

This pass must be net-negative for:

```text
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Expected source of deletions: shrink duplicated scalar aliases and related
artifact-writer boilerplate. If that is not enough, find a safe stale scaffold
deletion; do not delete endpoints or reference tests.

Validation:

```text
julia --project=. test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
git status --short --branch
```

Report back:

- what report aliases were removed or retained and why;
- confirmation that artifact keys stayed stable;
- exact source-plan producer blocker diagnosis;
- recommended pass-195 implementation seam;
- scoped line-budget arithmetic;
- validation results;
- deletion/shrinkage report with exact remaining caller/blocker.

-- repo-manager@macmini
