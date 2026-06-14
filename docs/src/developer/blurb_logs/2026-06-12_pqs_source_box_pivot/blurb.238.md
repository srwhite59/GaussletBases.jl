# Pass 238 blurb - independent H2 PQS source-plan materializer audit

Role: repo-doer.

Read before working:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_manager_running_log.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.237.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.237.md`
- `docs/src/developer/old_flat_cartesian_retirement_audit_2026-06-14.md`

Task type: no-edit audit.

Purpose:

The independent H2 PQS route now has generated support regions and a
route-owned retained-rule readiness plan:

```text
support counts:  (275, 578, 362)
retained counts: (275, 98, 98)
expected final dimension: 471
```

The remaining blocker is:

```text
:missing_independent_pqs_physical_source_plan_materializer
```

Before implementing that materializer, audit the exact seam. The risk is
accidentally importing fake-PQS/WL fixed-source coefficient data while trying to
make source-plan progress.

Audit questions:

1. What should the independent H2 PQS physical source-plan object contain at
   this stage?
   - support-region plan;
   - retained-rule plan;
   - per-unit source-box/source-mode descriptors;
   - per-unit retained-rule descriptors;
   - missing materialization objects.

2. For `:atom_contact_core`, can the source-plan represent direct source modes
   without coefficient matrices? Identify whether identity source modes are
   enough for the next source-plan pass, and what exact object should own them.

3. For `:shared_shell_1` and `:shared_shell_2`, identify which existing
   source-box/PQS machinery should supply the q=5 filled source boxes,
   one-dimensional source transforms, and boundary COMX product-mode retained
   rule descriptors.

4. Identify which existing modules/files are the right owners. Start from:
   - `src/pqs_source_box_diatomic_complete_core_shell.jl`
   - `src/pqs_multilayer_shell_source_plan.jl`
   - `src/CartesianRawProductSources.jl`
   - `src/CartesianShellification.jl`
   - `src/CartesianTerminalLowering.jl`
   - `src/CartesianRetainedUnits.jl`

5. Identify the exact old/flat paths that must not be used:
   - fake-PQS H2 463 source-backed WL/QW route;
   - `bond_aligned_diatomic_nested_fixed_source(...)`;
   - old fixed-source coefficient matrices;
   - dense parent or shell-row oracle paths as route authority.

6. Identify same-surface deletion candidates for the eventual implementation
   pass. Prefer stale blocker aliases, obsolete source-plan candidate report
   fields, or old route-shadow tests directly superseded by the materializer.

Do not edit files in this pass.

Do not implement:

- source coefficient matrices;
- final basis;
- H1, H1-J, RHF;
- supplements;
- CR2/export/public API;
- fake-PQS/WL comparison logic.

Report:

- recommended source-plan object shape;
- exact implementation seam for the next pass;
- existing functions/modules to reuse;
- blockers or missing primitives;
- forbidden paths confirmed avoided;
- deletion candidates for keeping the next implementation net-negative;
- whether the next pass should implement a source-plan payload or stop earlier.

-- repo-manager@macmini
