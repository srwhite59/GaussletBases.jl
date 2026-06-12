Pass 032 completed.

Files changed:

- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.032.md`

What changed:

- Updated the `CartesianRawProductSources` ownership wording so it owns raw source CPB/source-mode facts plus the narrow PQS source-mode boundary selector metadata tied to raw source ordering.
- Clarified that `CartesianRawProductSources` does not own general retained-rule policy, shell projection, Lowdin cleanup, final retained units, IDA weights, pair blocks, Hamiltonians, exports, or artifacts.
- Added a concise `2026-06-12 PQS Final-Basis H1 Status` section recording that the cubic `q=5/L=5` explicit final-basis H1 probe succeeds against a shell-support oracle.
- Marked that seam as oracle-backed and validated, not a fully production-owned PQS route, because shell projection and Lowdin inputs still come from the shell-realization/oracle layer.
- Kept IDA, density-density, RHF, driver adoption, exports, and artifacts explicitly out of scope.

Validation run:

- `git diff --check -- docs/src/developer/pqs_source_box_operator_framework.md` passed.
- `LC_ALL=C rg -n "[^ -~\\t]" docs/src/developer/pqs_source_box_operator_framework.md` found no non-ASCII bytes.
- `git status --short --branch` showed only the intended doc edit before response files were written.
- No Julia tests were run because this was a documentation-only pass.

Deletion/shrinkage report:

- Corrected stale ownership wording that implied `CartesianRawProductSources` owned only source CPBs/source-mode facts and generic axis transform facts, without naming the narrow PQS boundary selector carveout.
- Removed the stale under-specified "does not own retained rules ... or pair blocks" wording and replaced it with the narrower current contract: no general retained-rule policy, shell projection, Lowdin cleanup, final retained units, IDA weights, pair blocks, Hamiltonians, exports, or artifacts.
- The new status note earns its carrying cost by preventing the validated H1 seam from being mistaken for production route ownership.
- No source code, tests, route behavior, or compatibility paths became obsolete in this doc-only pass.
- Remaining stale/duplicate surfaces for the next pass: final-basis realization/operator-transfer functions still live outside a dedicated final-basis realization module, and the raw-source selector path remains separate from the planned direct retained-boundary overlap/kinetic kernels.

-- repo-doer@macmini
