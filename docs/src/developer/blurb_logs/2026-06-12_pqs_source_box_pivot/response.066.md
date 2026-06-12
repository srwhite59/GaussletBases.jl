Implemented pass 066 as an audit-only docs pass.

Changed:
- Added `Driver-Spine Integration Audit` to
  `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`.

Audit result:
- Mapped the current PQS objects/helpers onto the canonical driver lifecycle:
  - `cartesian_system` / `cartesian_recipe`: centers, charges, requested terms,
    source-box policy, and provisional fixture inputs.
  - `cartesian_parent`: mapped parent axes, PGDG/intermediate axis data, center
    tables, parent support counts, and raw parent/source facts.
  - `cartesian_shells`: shellification policy and core/outer box decomposition.
  - `cartesian_units`: direct core support, PQS shell source CPBs, raw product
    source plans, and retained boundary source-mode rules.
  - `cartesian_transforms`: shell projection/Lowdin data and the route-owned
    equivalent of `pqs_multilayer_shell_source_plan(...)`.
  - `cartesian_pairs`: retained source one-body blocks and pair/operator
    inventories.
  - `cartesian_assembly`: complete core/shell final basis, final one-body/H1,
    final IDA weights, and pre-final density interaction.
  - `cartesian_report` / `cartesian_materialization`: compact status, blocker,
    timing, and nonclaim summaries.
- Identified the first missing integration seam:
  `cartesian_assembly` does not yet consume a route-owned multi-layer PQS
  source plan and publish a complete core/shell final-basis payload for
  H1/J/RHF.
- Marked ready future driver inputs:
  `pqs_multilayer_shell_source_plan(...)`, direct retained PQS one-body blocks
  plus the generic retained selector, and `CartesianFinalBasisRealization`
  complete core/shell final-basis helpers.
- Marked private/oracle surfaces:
  `tmp/work` side-13 H1/J/RHF probes, shell-support projected operator oracles,
  old fixed-block matrices, and WL side-13 comparisons.
- Recorded that RHF acceptance, exports, artifacts, GTO supplement work, and a
  permanent side-13 fixture gate should not be integrated until the physical
  fixture rule is reviewed.

Validation:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - Passed.
- `git diff --check`
  - Passed.

Deletion/shrinkage report:
- No source, tests, metadata, or compatibility paths were deleted because this
  blurb was an audit/planning pass and explicitly forbade implementation work.
- No tests were added.
- The note identifies future shrink targets: probe-local support/operator
  assembly code and report-stage route-shadow fields that only restate private
  probe status, once the missing assembly seam exists.
- Remaining stale/duplicate surfaces to retire next: private side-13 probe
  assembly helpers should stay private until route-owned assembly consumes the
  multi-layer plan; then they can be deleted or reduced to oracle-only probes.

-- repo-doer@macmini
