Pass 205 - extract low-order/WL materialization from route-driver helpers

Role:
You are `repo-doer@macmini` implementing one bounded cleanup pass for
GaussletBases. Follow `AGENTS.md`, `JuliaStyle.md`, and `BlurbStyle.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should be:
  `f4da4cfb Split source box route shadow`
- `src/pqs_source_box_route_driver_helpers.jl` is still the main flatness
  hotspot.
- The private-global-overlap hook is already retired.
- The next extraction target is the route-configured White-Lindsey / low-order
  materialization subtree.

Goal:
Move the route-configured low-order/WL materialization implementation out of:

```text
src/pqs_source_box_route_driver_helpers.jl
```

into a new private file, preferably:

```text
src/pqs_source_box_low_order_materialization.jl
```

or, if the ownership reads better after inspection:

```text
src/cartesian_low_order_route_materialization.jl
```

This is an extraction-only pass. The driver helper may keep
`cartesian_materialization(...)` and/or a tiny
`_pqs_source_box_route_driver_materialization(...)` wrapper, but it should stop
owning the low-order materializer implementation.

Move these families if they are in `pqs_source_box_route_driver_helpers.jl` and
can be moved without behavior changes:
- `_pqs_source_box_route_driver_white_lindsey_preflight_fixed_block`
- `_pqs_source_box_route_driver_white_lindsey_ham_preflight`
- `_pqs_source_box_route_driver_one_center_materializer_probe`
- `_pqs_source_box_route_driver_route_configured_one_center_report`
- `_pqs_source_box_route_driver_diatomic_materializer_probe`
- `_pqs_source_box_route_driver_diatomic_atom_growth_materializer_probe`
- `_pqs_source_box_route_driver_route_configured_diatomic_basis_adapter`
- `_pqs_source_box_route_driver_route_configured_diatomic_basis_adapter_summary`
- `_pqs_source_box_route_driver_diatomic_atom_growth_basis_adapter`
- `_pqs_source_box_route_driver_route_configured_diatomic_ham_adapter`
- `_pqs_source_box_route_driver_route_configured_diatomic_ham_adapter_summary`
- `_pqs_source_box_route_driver_low_order_shellization_policy`
- the route-configured low-order materialization branch, if it can be moved
  cleanly.

Also inspect existing:

```text
src/cartesian_atom_growth_route_driver_helpers.jl
```

Do not duplicate code already owned there. If a helper already belongs there,
leave it there and only move the remaining low-order materialization glue out
of `pqs_source_box_route_driver_helpers.jl`.

Do not:
- change behavior;
- change return shapes or report field names;
- add PQS H1/H1-J/RHF logic;
- add artifact formats;
- delete `:legacy_diatomic_source` in this pass;
- split parent construction or general shell/unit/transform code in this pass;
- add tests.

Deletion audit only:
Audit whether `:legacy_diatomic_source` is still a live materialization path.
Report live callers/tests/driver inputs if found. Do not delete it in this pass.

Line-count rule:
The active line-count rule still applies. For edits under `src`, `test`, `bin`,
and the CR2 generator script, final tracked diff must be net-negative:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Acceptance condition:
`sum(deleted) > sum(added)` for those scoped files.

If this extraction cannot be made net-negative without behavior changes or
deleting live endpoint/reference tests, write `.agent_handoffs/ATTENTION.md`
and stop.

Validation:
Run:

```sh
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Then run the smallest focused driver/materialization smoke that covers
`materialize_route=false` or the current driver default path. Prefer an existing
short test or checked-in driver input. If the only WL materialization smoke is
expected to run longer than 60 seconds, do not run it by default; report the
exact test/command and why it was skipped.

Always run:

```sh
git diff --check
git diff --cached --check
```

after staging if you stage.

Response file:
Write `.agent_handoffs/response.205.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.205.md
```

Report:
- new file name chosen and why;
- functions/families moved;
- functions deliberately left behind and why;
- line counts before/after for `pqs_source_box_route_driver_helpers.jl` and the
  new file;
- source/test/bin scoped line budget added/deleted/net;
- `:legacy_diatomic_source` live-path audit;
- validation commands and results;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Stop after writing the response. Manager will review later.

-- repo-manager@macmini
