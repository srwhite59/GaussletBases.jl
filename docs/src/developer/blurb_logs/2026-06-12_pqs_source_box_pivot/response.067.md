Implemented pass 067.

Changes:
- Added `pqs_multilayer_complete_core_shell_final_basis(plan; ...)` in
  `src/pqs_multilayer_shell_source_plan.jl`.
  - Requires `plan.object_kind == :pqs_multilayer_shell_source_plan`.
  - Returns a compact blocked result if the source plan is not available.
  - Builds core/core, core/shell, and shell/shell overlap blocks from
    `plan.core_support_states`, `plan.shell_support_states`, and
    `plan.metrics.x/y/z.overlap`.
  - Delegates to
    `CartesianFinalBasisRealization.pqs_complete_core_shell_final_basis(...)`.
  - Preserves nonclaims: no one-body operators, H1, IDA, density-density, RHF,
    driver wiring, exports, or artifacts.
- Added a narrow private rectangular support-product helper in the same file.
  Existing production `_nested_support_product_matrix(...)` only covers square
  same-support blocks, while this seam needs the rectangular core/shell overlap.
- Added one compact check to
  `test/nested/pqs_direct_retained_final_h1_runtests.jl` comparing the new
  helper result with the existing manual final-basis construction on the
  established small one-shell fixture.
- Updated
  `docs/src/developer/pqs_near_term_final_basis_realization_plan.md` to mark
  this narrow assembly seam implemented, while keeping side13 RHF and fixture
  acceptance nonclaims unchanged.
- No exports were added.

Validation:
- `julia --project=. -e 't = @elapsed include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); println("elapsed_s=", t)'`
  - Passed: `PQS complete core-shell final H1 gate | 33`.
  - Elapsed: `62.99s`, including precompilation/startup.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - Passed.
- `git diff --check`
  - Passed.

Notes:
- The first focused test run exposed a keyword-only call bug when delegating to
  `pqs_complete_core_shell_final_basis(...)`; that was fixed before final
  validation.

Deletion/shrinkage report:
- No production code was deleted in this pass because the helper is the first
  route-owned replacement for probe-local final-basis assembly.
- No compatibility shim or new public export was added.
- The new test is genuinely new live-contract coverage for the helper and uses
  the existing H1 gate fixture rather than adding a standalone test file.
- Probe-local repeated blocks in `tmp/work` are now less necessary as assembly
  authority; they should be simplified or rewritten to call the new helper when
  those probes are next touched.
- Remaining stale/duplicate surface to retire next: probe-local manual
  `plan -> overlap blocks -> pqs_complete_core_shell_final_basis` assembly in
  side13 H1/J/RHF scripts.

-- repo-doer@macmini
