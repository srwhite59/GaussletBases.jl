Review 077: accepted.

The H1 fixture now uses the shellification/lowering-backed region-plan source
plan as the active `fixture.plan`. Final-basis construction, support states,
support overlap/kinetic/nuclear matrices, and the H1 solve all consume that
region-backed plan. The explicit-box source plan remains only as
`explicit_box_plan` for compact bridge comparison.

The comparison stayed appropriately small: status, layer count, core/shell
support counts, shell retained count, final retained count, and final overlap
identity. No new broad test file or metadata vocabulary sweep was added.

Manager validation:

- `julia --project=. -e 't = @elapsed include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); println("elapsed_s=", t)'`
  passed, 44/44, elapsed about 5.94s.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage:

- explicit-box planning is no longer the active H1 route authority in the
  durable gate;
- the explicit-box entry point remains as a bridge/oracle comparison because
  descriptor/Lowdin source realization still delegates through it;
- the test did not shrink in line count, but it now exercises the intended
  authority boundary.

Follow-up:

- The explicit-box entry point still owns duplicate box-depth/layer-box
  arithmetic internally. A later cleanup can quarantine it more clearly or
  reduce it to a compatibility wrapper once source realization consumes the
  region plan directly.
- Keep dense support-space one-body helpers scoped to the H1 seam; do not grow
  them into the general PQS operator algorithm.

-- repo-manager@macmini
