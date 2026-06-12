Review 086:

Accepted. This was another clean same-module file split. Dense support-space
one-body machinery is now visibly scoped to the H1 seam/reference layer rather
than living at the top of the PQS source-realization file.

What I checked:

- `src/pqs_multilayer_support_one_body.jl` contains the support product,
  kinetic, and electron-nuclear by-center helpers.
- `src/pqs_multilayer_shell_source_plan.jl` now starts with explicit-box and
  region-plan source realization, not support operator construction.
- `src/GaussletBases.jl` includes the support one-body file after the region
  plan and before the source-plan file.
- The extraction preserves function names and return shapes.
- The support helpers use `_pqs_multilayer_property(...)` from the region-plan
  file through include order, without creating circular include pressure.
- No tests, physics behavior, driver wiring, density, RHF, exports, artifacts,
  result-shape redesign, or optimization work were added.

Validation rerun by manager:

- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  passed, 46 tests.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion / shrinkage:

- Conceptual responsibility separated: dense support-space one-body helpers are
  no longer mixed into source-plan realization.
- `src/pqs_multilayer_shell_source_plan.jl` shrank again and is less of a
  private route island.
- No obsolete production behavior was deleted because this was intentionally a
  behavior-preserving extraction.
- No test was added; the existing H1 gate remains the behavior check.

Remaining concern:

- `src/pqs_multilayer_shell_source_plan.jl` still owns complete core/shell
  final-basis adapter and H1 payload helpers. The next cleanup should extract
  that final-basis/H1 seam into a same-module file, leaving the source-plan file
  focused on source realization plus the explicit-box bridge.

-- repo-manager@macmini
