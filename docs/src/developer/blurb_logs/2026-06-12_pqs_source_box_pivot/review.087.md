Review 087:

Accepted. This completes the immediate file-boundary cleanup around the
multi-layer PQS H1 seam. The source-plan file now owns source realization and
the explicit-box compatibility bridge, while the final-basis/H1 payload lives
in its own same-module file.

What I checked:

- `src/pqs_multilayer_complete_core_shell_h1.jl` contains the complete
  core/shell final-basis adapter, H1 payload, and blocked final-basis helper.
- `src/pqs_multilayer_shell_source_plan.jl` now contains only layer-spec
  helpers, common source realization, and the public source-plan entry points.
- `src/GaussletBases.jl` includes files in dependency order: region plan,
  support one-body, source plan, then complete core/shell H1.
- Function names and return shapes were preserved.
- No tests, physics behavior, driver wiring, density, RHF, exports, artifacts,
  result-shape redesign, or optimization work were added.

Validation rerun by manager:

- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  passed, 46 tests.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion / shrinkage:

- Conceptual responsibility separated: complete core/shell final-basis
  adaptation and H1 payload assembly are no longer embedded in source
  realization.
- `src/pqs_multilayer_shell_source_plan.jl` is now substantially less of a
  private route island.
- No obsolete behavior was deleted because this was intentionally a mechanical
  extraction.
- No test was added; the existing H1 gate remains the behavior check.

Remaining concern:

- The explicit-box bridge remains compatibility/probe machinery.
- The dense support one-body helpers and final-basis/H1 seam are now separate,
  but future H1/J work must avoid growing those files into a private density or
  RHF route island.
- Any move of the final-basis/H1 seam into `CartesianFinalBasisRealization`
  should be a later design decision, not an accidental mechanical cleanup.

-- repo-manager@macmini
