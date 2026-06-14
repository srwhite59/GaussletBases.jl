# Pass 235 blurb - delete route-shadow density-density fixture pressure

Role: repo-doer.

Read before editing:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_manager_running_log.md`
- `docs/src/developer/old_flat_cartesian_retirement_audit_2026-06-14.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.234.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.234.md`

Task:

Pay down the pass-234 line-budget exception by deleting old flat/route-shadow
density-density fixture pressure. This is a cleanup pass, not a physics or
independent-H2 implementation pass.

Primary target:

```text
test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl
```

and corresponding route-shadow fixture helpers in:

```text
src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl
```

Focus on the PQS/PQS/product route-shaped density-density fixture path and raw
box density-density route producer/consumer sections called only by the slow
integration test.

Goal:

- Remove or sharply shrink route-shaped density-density fixture assertions,
  nonclaim metadata checks, helper-name assertions, and old route-shadow
  vocabulary.
- Delete source helpers only if caller search proves they are only used by the
  deleted/trimmed slow test sections.
- Preserve compact live math checks if no newer test covers them:
  - boundary retained count convention;
  - density-normalized vs raw-weighted pair convention;
  - nuclear charge/sign convention.

Forbidden:

- no changes to independent H2 PQS support-plan implementation;
- no H1, H1-J, RHF, supplements, CR2, export, public API;
- no deletion of fake-PQS guard fields or fake-PQS golden regression;
- no deletion of active layered route spine modules;
- no broad rewrite of current-route metadata export beyond the identified
  fixture pressure.

Decision rules:

- If a helper has non-test source callers, do not delete it. Report exact
  callers and stop at test shrinkage.
- If only tests call a helper and the test section is route-shadow pressure,
  delete both helper and test pressure.
- If a scientific/math convention would be lost, preserve it in the smallest
  existing or new compact test. Do not keep a 6000-line integration section just
  to preserve metadata vocabulary.

Line-budget target:

- This pass should be strongly net-negative in `src + test + bin`, ideally more
  than the `+108` debt from pass 234.
- Measure with `git diff --numstat -- src test bin`.

Validation:

- Run package load.
- Run the smallest relevant focused test if one remains after shrinkage.
- Run `git diff --check`.
- Do not run expensive physics.

Report:

- exact helpers/tests deleted or shrunk;
- caller-search evidence for deleted source helpers;
- preserved math/scientific checks, if any;
- validation and timing;
- scoped line-budget totals;
- deleted/simplified/quarantined/not-deleted accounting.

-- repo-manager@macmini
