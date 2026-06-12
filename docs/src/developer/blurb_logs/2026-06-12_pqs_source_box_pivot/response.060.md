Pass 060 complete: implemented the narrow multi-layer PQS shell source plan and ran the side13 final-basis/H1 smoke.

Files edited:

- `src/GaussletBases.jl`
- `src/pqs_multilayer_shell_source_plan.jl`
- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`

Ignored artifacts written:

- `tmp/work/pqs_multilayer_shell_side13_h1_smoke.jl`
- `tmp/work/pqs_multilayer_shell_side13_h1_smoke_summary.txt`

Helper/object names:

- Added `pqs_multilayer_shell_source_plan(...)`.
- Added internal helpers:
  - `_pqs_multilayer_axis_metrics`
  - `_pqs_multilayer_box_depth`
  - `_pqs_multilayer_core_box_at_depth`
  - `_pqs_multilayer_block_concatenate_shell_coefficients`
  - `_pqs_multilayer_duplicate_count`

Implementation shape:

- The helper accepts a parent axis bundle plus ordered `core_box` and `outer_box`.
- It builds repeated legal one-cell projected-q-shell descriptors from the core outward.
- It calls the existing `CartesianContractedParentMetrics._pqs_shell_realization_plan` for each layer.
- It validates shell support disjointness, core/shell disjointness, and full intended outer-box coverage.
- It concatenates shell support rows and block-diagonal shell final coefficients into one collapsed shell sector.
- It returns support indices/states, shell records, metrics, collapsed shell coefficients, and compact summary fields for `pqs_complete_core_shell_final_basis(...)`.
- It does not build final-basis overlap, H1, IDA, RHF, driver wiring, exports, or artifacts.

Final-basis/H1 smoke result:

- Probe: `tmp/work/pqs_multilayer_shell_side13_h1_smoke.jl`
- Fixture:
  - parent count: 13
  - mapping: `AsinhMapping(c=0.1,s=1.0,tail_spacing=10.0)`
  - endpoints: `(-8.565228460168399, 8.565228460168399)`
  - core: `(4:10)^3`
  - outer box: `(1:13)^3`
  - shell layers: 3
- Plan:
  - status: `available_pqs_multilayer_shell_source_plan`
  - source kind: `repeated_one_cell_projected_q_shell_layers`
  - core support count: 343
  - shell support count: 1854
  - shell final retained count: 1206
  - combined support count: 2197
  - intended support count: 2197
  - support coverage: true
  - shell duplicate count: 0
  - core/shell duplicate count: 0
- Shell records:
  - shell1 `(3:11)^3 / (4:10)^3`, raw dims `(7,7,7)`, support 386, retained 218
  - shell2 `(2:12)^3 / (3:11)^3`, raw dims `(9,9,9)`, support 602, retained 386
  - shell3 `(1:13)^3 / (2:12)^3`, raw dims `(11,11,11)`, support 866, retained 602
- Final basis:
  - status: `available_pqs_complete_core_shell_final_basis`
  - final dimension: 1549
  - final overlap identity error: `5.5061337433937041e-13`
- H1:
  - Z = 1 lowest energy: `-0.494223730383033`
  - Z = 1 error vs `-0.5`: `0.005776269616967`
  - Z = 2 lowest energy: `-1.975561823201342`
  - Z = 2 error vs `-2.0`: `0.024438176798658`
- Timing:
  - multilayer plan: `2.248074` s
  - final basis: `0.901341` s
  - support one-body build: `1.083695` s
  - H1: `0.873846` s
  - total: `5.106957` s

Final-basis helper status:

- `pqs_complete_core_shell_final_basis(...)` stayed unchanged.
- The new helper supplies the collapsed shell sector needed by that existing helper.

Driver-stage shape:

- The seam is shaped as a shells/transforms/final-basis input: parent axis bundles plus core/outer boxes become shell records, support rows, shell isometries, and summary coverage data.
- Lower modules do not call upward into driver helpers.
- The helper is route-owned planning data, not a private H1/RHF lane.

Docs update:

- Updated the PQS near-term plan with the implemented seam and the side13 H1 smoke result.
- The note explicitly says this is not an acceptance gate or RHF fixture and records the spacing caveat.

Validation:

- `julia --project=. tmp/work/pqs_multilayer_shell_side13_h1_smoke.jl` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage report:

- No old production code was deleted in this pass because this was the first route-owned replacement for the missing multi-layer shell source seam.
- No permanent test was added; the scoped `tmp/work` smoke was enough for this implementation pass and avoids carrying a new long-running gate.
- The failed pass-058 single-descriptor multi-shell shape is superseded by `pqs_multilayer_shell_source_plan(...)`.
- The one-shell q ladder is less necessary as future physical-fixture pressure; it remains route/scaling evidence only.
- Remaining stale/duplicate surfaces to retire next: keep fixed-block matrices as oracle/reference only, and avoid promoting q-only one-shell PQS probes as He acceptance gates.

-- repo-doer@macmini
