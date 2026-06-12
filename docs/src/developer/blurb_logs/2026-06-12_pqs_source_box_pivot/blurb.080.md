Purpose:
  Prevent the explicit-box PQS multi-layer source entry point from reading as
  production route authority. It is now compatibility/bridge machinery; the
  active H1 fixture uses the shellification/lowering-backed region-plan path.

Task:
  Relabel the explicit-box entry point clearly as bridge-only in code comments,
  docstring text, and near-term docs.

Scope:
  - Update the `pqs_multilayer_shell_source_plan(bundles, core_box, outer_box; ...)`
    docstring to say explicit boxes are a compatibility/probe bridge and should
    not be used as shellification authority for new route work.
  - Add a compact metadata flag or summary field such as
    `explicit_box_bridge = true` for the explicit-box path if it fits cleanly.
  - Keep the region-plan path marked as shellification/lowering-backed route
    authority.
  - Update `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
    if it still describes the explicit-box entry point as route-owned source
    planning without the bridge caveat.
  - Only change `source_kind` if an `rg` scan shows no live tests/callers depend
    on the old symbol. If uncertain, leave `source_kind` stable and use metadata
    / doc wording instead.

Do not:
  - change numerical behavior;
  - change H1/RHF/IDA/density-density/driver behavior;
  - add tests;
  - add broad metadata assertions;
  - remove the explicit-box entry point;
  - request UI escalation. In unattended baton mode, write
    `.agent_handoffs/ATTENTION.md` and stop if permission is genuinely needed.

Validation:
  - `git diff --check`;
  - load check if source text changes beyond comments/docstrings/metadata.

Report:
  - what labels/docs changed;
  - whether `source_kind` was left stable or changed, and why;
  - validation run;
  - deletion/shrinkage report:
      - what conceptual authority was downgraded;
      - what explicit-box bridge remains;
      - whether any test coverage changed.

-- repo-manager@macmini
