Pass 075 complete.

Boundary recommendation:

`pqs_multilayer_shell_source_plan(...)` should stay a tactical bridge and should not become the owner of shellification geometry. The smallest next boundary is a shellification/lowering-backed PQS region plan that consumes `CartesianShellification.ShellificationPlan` plus `CartesianTerminalLowering.lower_terminal_regions(..., PQSLowering(...))`.

That intermediate object should carry:

- the direct core terminal region, its `owned_support`, and its raw/core box;
- ordered complete-shell layer records from shellification, including raw outer/inner boxes and owned support;
- the selected PQS lowering contracts/source CPBs for those shell regions;
- outer support coverage and disjointness/duplicate-support fingerprints from shellification/lowering provenance;
- explicit provenance that geometry came from `CartesianShellification` and `CartesianTerminalLowering`.

It should not carry PQS descriptors, shell projection/Lowdin matrices, support operator blocks, H1, IDA, RHF, or final-basis transfer data. A follow-up `pqs_multilayer_shell_source_plan(bundles, region_plan; ...)` can then consume that region plan and own only the PQS-specific source realization: projected-q shell descriptors, shell realization plans, collapsed shell-sector coefficients, retained counts, and support ordering.

Inspected files/functions:

- `src/cartesian_shellification/CartesianShellification.jl`
  - module contract/docstring;
  - `TerminalRegion`;
  - `ShellificationPlan`;
  - `_owned_support_for_region`;
  - `shellify`;
  - `terminal_regions`, `coverage`, `raw_plan`.
- `src/cartesian_shellification/terminal_geometry.jl`
  - `raw_terminal_geometry`;
  - `push_region!`;
  - single-atom shell growth;
  - coverage audit;
  - emitted raw region fields such as `outer_box`, `inner_exclusion_box`, `role`, `region_kind`, `shell_index`, and `support_count`.
- `src/cartesian_terminal_lowering/CartesianTerminalLowering.jl`
  - module ownership contract.
- `src/cartesian_terminal_lowering/contracts.jl`
  - `TerminalLoweringContract`;
  - `TerminalLoweringPlan`;
  - `source_cpbs`, `selected_contracts`.
- `src/cartesian_terminal_lowering/region_contracts.jl`
  - `_direct_terminal_contract`;
  - `_pqs_complete_shell_contract`;
  - `available_contracts`.
- `src/cartesian_terminal_lowering/selection.jl`
  - `lower_terminal_regions`;
  - selected-vs-available contract enumeration.
- `src/pqs_multilayer_shell_source_plan.jl`
  - `_pqs_multilayer_box_depth`;
  - `_pqs_multilayer_core_box_at_depth`;
  - `pqs_multilayer_shell_source_plan`.
- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`.

Fact split:

- Shellification/geometry ownership facts currently duplicated in `pqs_multilayer_shell_source_plan(...)`:
  - explicit `core_box`/`outer_box` interpretation;
  - symmetric equal shell-depth validation;
  - layer box growth via `_pqs_multilayer_core_box_at_depth`;
  - direct core support indices/states derived from the box;
  - intended outer support coverage;
  - core/shell and shell/shell duplicate checks;
  - ordered shell layer identity.
- PQS lowering/source-plan facts that belong in the PQS helper:
  - building one-cell projected-q shell layers from shellification/lowering-owned regions;
  - descriptor extraction;
  - `_pqs_shell_realization_plan`;
  - shell projection/Lowdin cleanup;
  - collapsed shell-sector coefficients and retained counts;
  - support ordering consumed by `pqs_multilayer_complete_core_shell_final_basis(...)`.

Can a small implementation pass follow?

Yes, for the one-center atom-outward case. `CartesianShellification` already exposes typed terminal regions with raw region boxes, owned support, ordering, roles/kinds, and coverage. `CartesianTerminalLowering` already selects PQS complete-shell contracts with source CPBs and carries the terminal region key/kind/owned support. The likely next pass is mechanical: build a compact `PQS` region-plan helper from shellification/lowering output, then add a new `pqs_multilayer_shell_source_plan` entry point that consumes it. I did not find a missing required geometry fact for the current one-center multi-layer use case.

Docs edited:

- Updated `docs/src/developer/pqs_near_term_final_basis_realization_plan.md` to record that `pqs_multilayer_shell_source_plan(...)` is a bridge, not shellification authority, and to define the shellification/lowering-backed boundary.

Validation:

- `git diff --check` passed.

Deletion/shrinkage report:

- No code was deleted in this docs/design pass because the blurb asked to prefer no source edits and the current helper remains the active bridge until the shellification/lowering-backed region plan exists.
- Responsibilities that should move out of `pqs_multilayer_shell_source_plan(...)` next:
  - box-depth arithmetic;
  - layer-box construction;
  - core/outer support coverage authority;
  - duplicate/disjointness checks;
  - geometry provenance.
- Probe/test glue that should become unnecessary after that move:
  - explicit test/probe construction of `core_box`/`outer_box` as PQS source-plan authority;
  - route-shadow checks that only restate private PQS box coverage;
  - any helper vocabulary asserting that the PQS source helper owns shell growth.
- No new tests were added; this is docs-only live-contract clarification.
- Remaining stale surface to retire next: the explicit-box `pqs_multilayer_shell_source_plan(bundles, core_box, outer_box; ...)` entry point should be quarantined as a bridge once the shellification/lowering-backed entry point exists.

-- repo-doer@macmini
