Pass 200 - define the physical H2 gausslet-only target before implementation.

Purpose:

The current H2 221-dimensional route is now explicitly labeled
`source_box_diagnostic`. Do not continue it to H1-J/RHF.

The next useful H2 target is a physical gausslet-only driver endpoint with full
retained atom-core interiors plus shell layers. Before implementing it, recover
the old WL gausslet-only inventory and define the PQS analog precisely enough
that implementation does not guess the basis shape.

Task type:

Audit/design plus deletion. No route implementation.

Read/inspect:

```text
docs/src/developer/high_order_endcap_panel_h2_chemistry_reproduction_2026-05-16.md
docs/src/developer/high_order_mainline_import_readiness_2026-05-15.md
docs/src/developer/blurb_logs/2026-06-11_mwg_decomposed_wl_cleanup/blurb.009.md
test/diatomic/runtests.jl
test/nested/*diatomic*
src/pqs_source_box_diatomic_complete_core_shell.jl
test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl
```

Audit questions:

1. Old WL gausslet-only H2 inventory.

   The documented default row says:

   ```text
   fixed block size = (1215, 463)
   final dimension with supplement = 481
   residual count = 18
   shared layer columns = (98, 114)
   ```

   Confirm what the `463` gausslet-only columns represent:

   ```text
   atom core columns
   shell/layer columns
   midpoint/bridge columns, if any
   retained unit/order convention
   whether full 5^3 atom cores are included
   ```

   Do not infer from arithmetic alone if a tracked test/doc/object exposes the
   inventory.

2. Physical H2 PQS target.

   Define the first gausslet-only PQS target that should replace the current
   221 diagnostic route. It should be comparable in parent geometry and
   shellification policy to the old WL gausslet-only H2 route, but without the
   H/cc-pVTZ residual supplement.

   State:

   ```text
   source/final units
   support row order
   retained order
   expected retained counts if known
   expected final dimension if known
   whether midpoint/product slab remains
   how full atom-core interiors are retained
   how shell layers are represented
   ```

3. Implementation seam.

   Recommend the smallest implementation seam for pass 201. Examples:

   ```text
   new route_kind = :bond_aligned_diatomic_atom_core_shell_pqs
   or extend existing source-realization payload with full-core retained units
   or reuse an existing complete-rectangular shellification object
   ```

   Be explicit about which existing helpers should be reused and which current
   diagnostic route objects should not be mutated.

4. Artifact status.

   State what the future physical H2 artifact should say before H1-J/RHF:

   ```text
   artifact_role
   physics_endpoint_ready
   retained_atom_core_interiors
   source_plan_role
   comparison_ready
   supplement_policy
   ```

5. Delete one safe stale scaffold.

   Requirements:

   - no include in default/integration/top-level runners;
   - no source caller;
   - not an explicit He/H2 driver endpoint;
   - not a WL H/H2 reference or scientific endpoint;
   - not needed by this audit.

   If no safe deletion exists, write `.agent_handoffs/ATTENTION.md` and stop.

Trust boundary:

- No physical H2 route implementation.
- No H2 source-plan/final-basis changes.
- No H1-J/density interaction.
- No private RHF for H2.
- No supplemented WL/QW comparison.
- No supplement support.
- No Be2/Cr2 artifact work.
- No HFDMRG, DMRG, ECP, exports, public solver behavior, or Qiu-White
  correction work.
- Preserve the visible staged driver style.
- Do not request interactive approval or sandbox escalation. If approval is
  genuinely required, write `.agent_handoffs/ATTENTION.md` and stop.

Line-budget rule:

This pass must be net-negative for:

```text
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Expected implementation diff is deletion only. The audit response/log does not
count toward this source/test/bin budget.

Validation:

Use read-only inspection and focused commands only. Do not run broad suites.
At minimum:

```text
rg -n "<deleted_file_basename>|<deleted_testset_name>" test src bin
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
git status --short --branch
```

Report back:

- old WL gausslet-only inventory findings;
- proposed physical H2 PQS target;
- recommended pass-201 implementation seam;
- stale scaffold deleted and why safe;
- validation results;
- deletion/shrinkage report with exact remaining caller/blocker.

-- repo-manager@macmini
