Pass 192 - audit H2 driver support and supplement policy, with one stale-test deletion.

Purpose:

He q=5/n_s=5 PQS is now driver-owned for H1/H1-J and optional private RHF. The
next physics target is H2, but pass 191 found an important blocker: the old
WL/QW H2 reference HF totals include the H/cc-pVTZ S/P residual supplement.

Do not implement H2 yet. This pass should answer whether the current visible
driver can express the intended H2 PQS target with the same supplement/residual
policy, or whether the first H2 driver endpoint must be a gausslet-only
diagnostic with WL HF comparison explicitly disabled.

This pass also keeps the source/test/bin/generator line budget negative by
deleting one stale development scaffold if safe.

Task type:

Audit plus deletion. No H2 implementation.

Physics target:

Future driver-owned H2 at `R = 4.0` bohr, bond axis `:z`, centered atoms
`(0,0,-2)` and `(0,0,2)`, nuclear charges `(1,1)`, `:G10`,
`core_spacing = 0.5`, `xmax_parallel = 6.0`, `xmax_transverse = 4.0`,
`q = 5`, `n_s = 5`, default complete-rectangular WL/QW analog first.

Important pass-191 facts:

- default complete-rectangular WL/QW reference:
  - final dimension `481`
  - residual count `18`
  - HF total `-0.910938264352`
  - ED total `-1.015613837691`
- endcap/panel reference:
  - final dimension `461`
  - HF total `-0.910977315003`
  - ED total `-1.015663743783`
- These documented HF/ED references include the H/cc-pVTZ S/P residual
  supplement. Do not compare a gausslet-only PQS route directly to these
  numbers.

Read/inspect:

```text
bin/cartesian_ham_builder.jl
test/driver_inputs/he_pqs_q5_wlmap.jl
src/pqs_source_box_diatomic_complete_core_shell.jl
src/pqs_source_box_route_driver_helpers.jl
docs/src/developer/high_order_endcap_panel_h2_chemistry_reproduction_2026-05-16.md
docs/src/developer/high_order_mainline_import_readiness_2026-05-15.md
docs/src/developer/blurb_logs/2026-06-11_mwg_decomposed_wl_cleanup/blurb.009.md
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.191.md
```

Use `rg` first. Do not run broad tests.

Audit questions:

1. Driver expressibility:

   Can `bin/cartesian_ham_builder.jl`, as a visible staged laboratory script,
   currently express a bond-aligned diatomic PQS source-box complete-rectangular
   target with the pass-191 H2 parameters?

   Report exact current surfaces and blockers. For example:

   ```text
   system/geometry inputs:
   parent mapping/extents inputs:
   route_family/route_kind support:
   diatomic shellification policy support:
   fixed-q/source-mode support:
   final-basis support:
   H1 support:
   H1-J/density-interaction support:
   private RHF support:
   save/artifact support:
   ```

2. Supplement/residual policy:

   Where, if anywhere, is the H/cc-pVTZ S/P residual supplement represented in
   the current driver or route-owned code?

   Distinguish:

   ```text
   old WL/QW reference path supports supplement
   current visible driver supports supplement
   PQS route-owned source-box path supports supplement
   artifact/save path can record supplement status
   ```

   If supplement support is missing from the PQS driver path, say that plainly.
   Do not synthesize a supplement adapter in this pass.

3. First H2 implementation decision:

   Recommend one of these, with reasons:

   ```text
   A. Implement supplemented H2 PQS first, because the driver can already carry
      the same H/cc-pVTZ S/P residual policy.

   B. Implement gausslet-only H2 PQS first, with comparison_ready=false and no
      direct HF comparison to the old WL/QW supplemented total.

   C. Stop before H2 implementation because a lower-level route/support
      convention is still ambiguous.
   ```

   The manager will choose based on this audit; do not implement the choice in
   this pass.

4. Future input/artifact shape:

   If implementation is safe soon, sketch the exact future driver input keys
   for `test/driver_inputs/h2_pqs_q5_wlmap.jl`.

   The artifact should keep final-basis self-overlap as a scalar identity-error
   diagnostic only. Do not request downstream storage or consumption of an `S`
   matrix for an orthonormal final working basis.

5. Delete one stale scaffold:

   Find one development scaffold test/probe whose role is superseded and whose
   deletion is safe. Requirements:

   - no include in `test/nested/runtests.jl`;
   - no include in `test/nested/integration_runtests.jl`;
   - no include in top-level/group runners;
   - no source caller;
   - not one of the explicit He driver endpoint tests;
   - not a WL He/H/H2 acceptance/reference test;
   - not an active scientific endpoint;
   - not needed for the current H2 audit.

   If no safe deletion exists, write `.agent_handoffs/ATTENTION.md` and stop.

Trust boundary:

- No H2 implementation.
- No new driver input file.
- No new source code.
- No new tests.
- No RHF/SCF algorithm changes.
- No Be2/Cr2 artifact work.
- No HFDMRG, DMRG, exports, public solver behavior, ECP, or Qiu-White
  correction implementation.
- Preserve the visible staged driver style. Do not propose hiding construction
  behind an opaque `run_driver(config)` wrapper.
- Do not request interactive approval or sandbox escalation. If approval is
  genuinely required, write `.agent_handoffs/ATTENTION.md` and stop.

Line-budget rule:

This pass must be net-negative for:

```text
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Expected:

```text
added: 0
deleted: >0
net: negative
```

The audit response and curated log do not count toward source/test/bin budget.

Validation:

```text
rg -n "<deleted_file_basename>|<deleted_helper_or_testset_name>" test src bin
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
git status --short --branch
```

Report back:

- supplement support status and exact blockers;
- recommended first H2 implementation option A/B/C;
- future driver input key sketch, if safe;
- stale scaffold deleted and why it is safe;
- validation results;
- deletion/shrinkage report with exact remaining caller/blocker.

-- repo-manager@macmini
