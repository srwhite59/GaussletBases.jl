What changed:

- Shrank `test/nested/cartesian_wl_gausslet_he_atom_acceptance_runtests.jl` from 1020 lines to 744 lines.
- Removed the test-local phase timing harness and broad timing fields from the He acceptance audit. Dedicated probes and production `TimeG` scopes remain the performance surfaces.
- Trimmed the He acceptance report to the active scientific contract instead of printing a wide report field cloud.
- Reduced the assertion set to the live endpoint checks:
  - fixture identity, retained dimension, decomposed unit count, and pair count;
  - shellification-backed decomposed WL inventory and anti-seed check;
  - route-global overlap, kinetic, electron-nuclear, Hamiltonian, and density-density materialization;
  - overlap positive-definiteness/rank;
  - H1 energy and hydrogenic `1s` self-Coulomb diagnostics;
  - key IDA retained-weight boundary convention;
  - RHF convergence, energy decomposition, density trace/electron count, and accepted small-He energy baseline;
  - anti-fallback checks for full-parent CPB, direct Cartesian product assembly, and ordinary Cartesian IDA.

Validation:

- `julia --project=. test/nested/cartesian_wl_gausslet_he_atom_acceptance_runtests.jl`
  - pass: 79 / 79
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - pass: `load ok`
- `git diff --check`
  - pass

Audits:

- `test/nested/runtests.jl` was inspected. The inline owned-unit coverage block remains because it is compact module-contract coverage for dense/sparse retained-unit construction and coverage failure modes, not an obsolete acceptance-route test.
- `test/nested/integration_runtests.jl` was not edited; no runner boundary changes required an integration gate.

Deletion/shrinkage:

- Deleted 391 lines and added 115 lines in the He acceptance file, for net 276 lines removed.
- Removed stale/low-value assertion pressure around:
  - broad report-shape fields;
  - route/status duplicates;
  - detailed pair-factor metadata vocabulary;
  - retained-density min/max/count metadata;
  - repeated nonclaim flags;
  - timing nonnegativity and timing relationship assertions;
  - density peak/core-fraction diagnostics that were not part of the active scientific contract.
- Simplified the test-local audit helper by removing phase `@elapsed` wrappers and pre-initialization scaffolding.
- No new tests were added. This pass shrank an existing scientific acceptance test instead of adding coverage.
- No production code changed, so no production compatibility path became obsolete in this pass.
- Remaining stale/duplicate surface to consider later: the He acceptance helper still builds a low-order seed report as parent-axis setup input before using shellification-backed inventory. That is outside this deletion-only test pass, but it remains worth retiring when the standard driver-facing parent/mapping setup fully replaces the seed bridge.

Trust boundary:

- No GTO/MWG implementation path was touched.
- No physics baseline changed.
- No full-parent CPB, direct Cartesian product, or ordinary Cartesian IDA fallback was introduced.
- No broad slow integration gate was run.

-- repo-doer@macmini
