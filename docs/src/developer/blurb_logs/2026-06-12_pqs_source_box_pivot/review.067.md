Review 067: accepted.

`pqs_multilayer_complete_core_shell_final_basis(plan; ...)` is the right narrow
assembly seam. It consumes an available multi-layer source plan, builds the
core/core, core/shell, and shell/shell overlap blocks from the plan support
states, and delegates to `CartesianFinalBasisRealization`. It does not build
H1, IDA, density-density, RHF, driver wiring, exports, or artifacts.

One remaining cleanup is now visible: the tracked H1 fixture still builds
`fixture.final_basis` manually, then compares the new helper result against it.
That was acceptable for the first validation, but it leaves duplicate
test-side assembly as a parallel authority. The next pass should flip the H1
gate to use the helper as the normal final-basis construction path and delete
the manual `core_overlap/core_shell_overlap/shell_overlap -> pqs_complete...`
fixture assembly where possible.

-- repo-manager@macmini
