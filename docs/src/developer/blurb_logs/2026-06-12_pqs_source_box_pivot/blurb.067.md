Purpose:
  Move the repeated probe-local "multi-layer shell source plan -> complete
  core/shell final basis" assembly step into a narrow route-owned helper.

Context:
  Pass 066 identified `cartesian_assembly` as the first missing driver seam.
  Current side13 probes repeatedly do this locally:

  ```julia
  plan = pqs_multilayer_shell_source_plan(...)
  core_overlap = support_product(plan.core_support_states, plan.core_support_states, ...)
  core_shell_overlap = support_product(plan.core_support_states, plan.shell_support_states, ...)
  shell_overlap = support_product(plan.shell_support_states, plan.shell_support_states, ...)
  final_basis = CartesianFinalBasisRealization.pqs_complete_core_shell_final_basis(...)
  ```

  That is stable enough to own as a small helper, but not broad enough to be a
  full driver integration.

Task:
  Add a narrow helper, likely near `pqs_multilayer_shell_source_plan(...)`, with
  a name like:

  ```julia
  pqs_multilayer_complete_core_shell_final_basis(plan; metadata = (;), ...)
  ```

  It should:

  - require `plan.object_kind == :pqs_multilayer_shell_source_plan`;
  - require `plan.status == :available_pqs_multilayer_shell_source_plan`;
  - compute core/core, core/shell, and shell/shell overlap blocks from
    `plan.core_support_states`, `plan.shell_support_states`, and
    `plan.metrics.x/y/z.overlap`;
  - call
    `CartesianFinalBasisRealization.pqs_complete_core_shell_final_basis(...)`;
  - return that final-basis result or a compact blocked result if the plan is
    blocked;
  - preserve nonclaims: no one-body operators, H1, IDA, density-density, RHF,
    driver wiring, exports, artifacts, or fixture acceptance.

  Reuse existing support product helpers if they are available in production
  code. Do not copy a large probe-only support-matrix implementation if a
  production helper already exists.

Validation:
  Add or update one compact test only if needed. Preferred check:

  - build a small multi-layer plan;
  - call the new helper;
  - compare its final overlap/error/final retained count to the existing manual
    call path;
  - avoid broad metadata assertions.

  Run:

  - the focused test covering the helper;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Docs:
  Update `docs/src/developer/pqs_near_term_final_basis_realization_plan.md` to
  mark this assembly seam as implemented only if the helper lands. Keep the note
  concise and keep side13 RHF non-acceptance language unchanged.

Do not:
  - implement H1/RHF/IDA/density-density in this pass;
  - add driver wiring;
  - add exports unless an existing local convention absolutely requires it;
  - study or codify `Z,d,s,ns`;
  - add a new physics fixture or gate;
  - move tmp/work probe scripts into tracked tests.

Deletion/shrinkage report required:
  - what old code, test, metadata, or compatibility path became unnecessary;
  - what was deleted or simplified;
  - if nothing was deleted, why no existing surface was made obsolete;
  - whether any new test replaces/shrinks older coverage or is genuinely new
    live-contract coverage;
  - any remaining stale or duplicate surfaces to retire next.

-- repo-manager@macmini
