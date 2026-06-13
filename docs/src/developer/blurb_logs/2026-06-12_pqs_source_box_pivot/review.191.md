Pass 191 review

Accepted.

The audit did what this pass needed: it recovered the old WL/QW H2 reference
without implementing H2, and it kept the line-budget pressure real by deleting
one stale scaffold test rather than adding another readiness layer.

Key accepted H2 facts:

- The old H2 reference is the WL/QW `R = 4.0` bohr fixture with H centers at
  `(0,0,-2)` and `(0,0,2)`, nuclear charges `(1,1)`, `:G10`,
  `core_spacing = 0.5`, `xmax_parallel = 6.0`, `xmax_transverse = 4.0`, and
  bond axis `:z`.
- The documented default complete-rectangular WL/QW row has final dimension
  `481`, residual count `18`, HF total `-0.910938264352`, and ED total
  `-1.015613837691`.
- The documented endcap/panel row has final dimension `461`, residual count
  `18`, HF total `-0.910977315003`, and ED total `-1.015663743783`.
- The documented WL/QW HF totals include the H/cc-pVTZ S/P residual supplement,
  so a future PQS gausslet-only driver target must not compare directly to
  those HF totals. Either the same supplement/residual policy must be carried,
  or the first comparison must be restricted to comparable inventory/H1/H1-J
  diagnostics.

The proposed first target is correct: start from the PQS source-box analog of
the default complete-rectangular WL/QW H2 route, not the endcap/panel route and
not a mismatched scalar `q = n_s` comparison. The equivalence condition is the
parent mapping, shellification/retained-basis policy, supplement policy, and
density-interaction convention.

Deletion accepted:

- Deleted `test/nested/cartesian_pair_block_one_body_block_set_consumption_skeleton_runtests.jl`.
- This was an unreferenced metadata-only mixed one-body scaffold.
- The smaller live smoke `cartesian_pair_block_one_body_consumer_smoke_runtests.jl`
  remains, and no physics endpoint/reference test was removed.

Validation reviewed:

- deleted-file search had no live `src/test/bin` callers;
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed;
- `git diff --check` passed;
- scoped line budget is `0 added / 205 deleted`.

Next manager decision:

Do not implement H2 until the next pass audits the driver support needed for a
bond-aligned diatomic fixed-q complete-rectangular PQS endpoint, especially the
residual supplement question. The next pass should keep the net-negative
source/test/bin/generator line budget active.

-- repo-manager@macmini
