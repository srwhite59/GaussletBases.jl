Pass 209 - audit physical H2 source-plan construction inputs

Role:
You are `repo-doer@macmini` doing a no-edit audit for GaussletBases. Follow
`AGENTS.md`, `JuliaStyle.md`, `BlurbStyle.md`, and the unattended baton rules
in `.agent_handoffs/RUN.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should be:
  `5de14d69 Add blocked H2 physical source plan seam`
- The H2 463 physical target source-plan seam exists but is blocked on:
  `:missing_atom_contact_core_support_rows`
- The target inventory remains:
  support counts `(275, 578, 362)`, retained counts `(251, 98, 114)`,
  expected final dimension `463`, no supplement.

Task type:
No-edit audit. Do not modify source, tests, docs, inputs, or generated files.

Purpose:
Identify the exact live construction source for the H2 physical source plan:

```text
atom-contact core support rows / coefficients
shared shell 1 support rows / coefficients
shared shell 2 support rows / coefficients
```

The pass should answer whether the existing atom-growth / complete-rectangular
materialization path already constructs those objects in a form the PQS physical
source-plan seam can consume, or whether a new route-owned producer is still
missing.

Known surfaces to inspect:
- `src/cartesian_shellification_plan.jl`
  - `_cartesian_materialize_atom_growth_complete_rectangular_sequence_low_order`
  - `_cartesian_materialize_atom_growth_complete_rectangular_shellification_low_order`
  - `_cartesian_materialize_atom_local_child_shellification_low_order`
  - `_cartesian_materialize_shared_complete_shell_region`
  - `_cartesian_materialize_direct_box_region`
- `src/bond_aligned_diatomic_geometry.jl`
  - `_bond_aligned_nested_shell_provenance`
  - `_bond_aligned_source_region_points`
- `src/pqs_source_box_route_driver_skeletons.jl`
  - `_pqs_source_box_route_driver_physical_gausslet_core_shell_skeleton`
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
  - `_pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_payload`
  - existing diagnostic `_PQSDiatomicCompleteCoreShellSourcePlan`
- Any current route-configured atom-growth materializer path moved to:
  `src/pqs_source_box_low_order_materialization.jl`

Questions to answer:
1. Where do the atom-contact-core support rows and coefficient blocks currently
   exist, if anywhere?
   - Are they in `core_support_indices`, `core_support_blocks`,
     `core_coefficient_blocks`, `sequence.core_column_range`, child sequences,
     `contact_cap_data`, or another object?
   - Do their counts match support `275` and retained `251` for H2 R=4 q=5?

2. Where do shared shell 1 and shared shell 2 rows/coefficients exist?
   - Are they available as `shared_shell_layers`, `shared_shell_column_ranges`,
     or shell provenance?
   - Do they match support counts `578`, `362` and retained counts `98`, `114`?

3. Is the existing materialization path low-order/WL-only, or can it be treated
   as the physical gausslet source-plan construction authority?
   - Pay close attention to flags such as `active_source_authority = false`,
     `private_development_only = true`, and `route_behavior_changed = false`.
   - Do not recommend adopting it as authority unless the convention and role
     are coherent.

4. If the needed objects exist, what is the smallest implementation seam for
   pass 210?
   - A wrapper from atom-growth materialization into a physical PQS source-plan
     object?
   - A new route-owned physical source-plan producer?
   - A blocked payload with a sharper missing object?

5. If the objects do not exist, what is the exact missing object?
   Use precise blockers, for example:

```text
:missing_atom_contact_core_coefficients
:missing_atom_contact_core_support_rows
:missing_shared_shell_layer_coefficients
:missing_physical_gausslet_source_plan_materializer
:low_order_materializer_not_source_authority
```

Do not:
- edit files;
- run or write H1/H1-J/RHF/HFDMRG/CR2;
- materialize final basis;
- compare to supplemented WL/QW H2 references;
- use the H2 221 source-box diagnostic source plan as the H2 463 source plan;
- turn low-order/WL diagnostic materialization into PQS source authority by
  wording alone.

Validation:
No Julia test is required for this audit. Use read-only `rg`, `sed`, `git grep`,
or short Julia introspection only if necessary. If you run any Julia command,
report it and keep it read-only.

Response file:
Write `.agent_handoffs/response.209.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.209.md
```

Report:
- current git status;
- exact files/functions inspected;
- whether atom-contact-core rows/coefs exist and where;
- whether shared-shell rows/coefs exist and where;
- whether existing materialization can be source authority or remains
  low-order/private diagnostic;
- recommended pass 210 implementation seam;
- exact blocker labels to use if blocked;
- deletion/shrinkage opportunity noticed during audit, if any.

Continue polling after writing the response. Manager will review, commit/push,
and publish the next blurb unless a real blocker appears.

-- repo-manager@macmini
