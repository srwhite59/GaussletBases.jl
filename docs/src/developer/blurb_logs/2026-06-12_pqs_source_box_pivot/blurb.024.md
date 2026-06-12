Purpose:

Audit the existing PQS retained centered electron-nuclear by-center path on the
real projected-q-shell fixture before extending final-basis one-body
materialization to nuclear attraction or H1.

Why now:

Pass 023 made overlap and kinetic transformable from proven retained boundary
operators into the shell-realized final basis. The next blocker is:

```text
:missing_pqs_shell_boundary_electron_nuclear_operator_source
```

Existing CPBM helpers already exist:

```text
pqs_source_pair_centered_electron_nuclear_by_center_block(...)
pqs_source_pair_retained_centered_electron_nuclear_by_center_block(...)
```

but they must be checked on the real projected-q-shell fixture against a
shell-support boundary oracle before they can be treated like overlap/kinetic.

Exact task:

Create a `tmp/work` probe. Do not change production code or tests unless the
audit exposes a tiny missing introspection seam.

Use the same real fixture pattern as passes 020 and 022:

- `current_box = (1:5, 1:5, 1:5)`
- `inner_box = (2:4, 2:4, 2:4)`
- `q = 5`, `L = 5`
- old `_pqs_shell_realization_plan(...)` only as oracle/input for `P` and `L`
- CPBM raw source plan, PQS boundary retained rule, and final-basis object

Build a one-center test center record, preferably at the origin first:

```text
center_key = :origin
center_index = 1
location = (0.0, 0.0, 0.0)
charge recorded but not applied
```

Use the existing centered PQS source helpers to build:

```text
retained_source_electron_nuclear_by_center
```

Then build an independent shell-support oracle for the same uncharged
by-center operator from support states, axis layers, and the same Coulomb
Gaussian expansion:

```text
V_shell_support(center)
V_boundary_oracle = P' * V_shell_support(center) * P
```

Compare:

```text
retained_source_nuclear_block vs V_boundary_oracle
```

and record:

- max error;
- term label;
- center key/index/location;
- nuclear charge recorded;
- nuclear charge applied is false;
- centers summed is false;
- uncharged by-center convention is true.

If the origin center passes cleanly and the probe remains simple, optionally
try one off-origin center to catch axis-centering mistakes. Do not expand into a
large test matrix.

Trust boundary:

This is an audit/probe only. Do not extend
`pqs_source_shell_final_one_body_from_boundary_matrix(...)` to nuclear in this
pass. Do not assemble H1. Do not sum centers or apply charge. Do not run IDA,
density-density, RHF, drivers, exports, or artifacts. Do not call
`_pqs_current_route_safe_term_matrices(...)`.

Questions to answer:

- Does retained-source electron-nuclear by-center equal the shell-projected
  boundary nuclear operator on the real fixture?
- Does it preserve charge-recorded/not-applied and centers-not-summed metadata?
- Is an off-origin center clean, if checked?
- What exact implementation should pass 025 do?

Test policy:

No permanent tests in this pass. Use `tmp/work` artifacts only.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete yet;
- whether any new probe artifact was added and why it earned its temporary
  carrying cost;
- any remaining stale or duplicate surfaces to retire next.

Validation:

- `julia --project=. <tmp/work probe>`;
- `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- `git diff --check`.

Report back:

- write `.agent_handoffs/response.024.md.tmp`, then atomically rename to
  `.agent_handoffs/response.024.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.024.md`;
- include probe artifact path;
- include nuclear comparison numbers or exact blocker;
- include recommended pass-025 implementation target;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
