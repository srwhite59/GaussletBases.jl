Pass 203 - simplify `CartesianContractedParentMetrics.jl` by splitting the true core from legacy route-shadow code.

Target card:

- Cleanup target: `src/CartesianContractedParentMetrics.jl`, currently about
  19k lines.
- Physics target protected: physical H2 PQS and He/WL endpoint work should not
  depend on a giant mixed metrics/route-shadow file.
- Allowed implementation surface: mechanical file split, include ordering,
  caller audit, and deletion/shrinkage of stale scaffolding proven not to be a
  live endpoint/reference.
- Forbidden surfaces: numerical behavior changes, H2 source-plan
  implementation, H1/H1-J/RHF work, WL comparison, CR2 solver/export work, new
  public APIs, new tests.
- Success condition: the exported metric-packet API lives in a visibly small
  core include, and non-core PQS/product/source-box/reporting/sidecar code is
  explicitly quarantined in legacy/private include files or deleted if proven
  stale.

Governing analysis:

- The public export surface is small:

  ```julia
  CartesianContractedParentMetricPacket3D
  cartesian_contracted_parent_metric_packet
  cartesian_contracted_parent_metric_packet_dense_reference
  contracted_parent_metric_packet_parent
  contracted_parent_metric_packet_overlap
  contracted_parent_metric_packet_weights
  contracted_parent_metric_packet_centers
  contracted_parent_metric_packet_diagnostics
  ```

- The coherent metric core stores retained overlap, retained weights, retained
  position matrices/centers, first moments, and diagnostics.
- The current 19k-line file also owns product/PQS source-box shadows, PQS
  current-route inventory/export, component smoke reports, CR2 sidecar schemas,
  density-density fixtures, local Gaussian helpers, and old route descriptors.
- That non-core code may still be called by integration tests, so do not delete
  blindly.

Read first:

- `AGENTS.md`, especially code cleanup, test deletion bias, and line-budget
  discipline.
- `JuliaStyle.md`.
- `docs/code_bloat_and_wrong_contract_cleanup_note.md`.
- Current file:
  `src/CartesianContractedParentMetrics.jl`.
- Live tests/callers:
  - `test/nested/cartesian_contracted_parent_metric_packet_runtests.jl`
  - `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`
  - `test/nested/pqs_component_route_report_adapter_runtests.jl`
  - `test/nested/pqs_source_metadata_real_artifact_acceptance_runtests.jl`
  - `test/nested/pqs_source_metadata_real_artifact_acceptance_support.jl`
  - `test/nested/integration_runtests.jl`

Implementation task:

1. Create a directory for includes:

   ```text
   src/cartesian_contracted_parent_metrics/
   ```

2. Split `src/CartesianContractedParentMetrics.jl` mechanically into a thin
   module wrapper plus named includes. Suggested first structure:

   ```text
   src/CartesianContractedParentMetrics.jl
     module, imports, exports, includes only

   src/cartesian_contracted_parent_metrics/core.jl
     CartesianContractedParentMetricPacket3D
     _AxisMetricData1D
     _ParentCoefficientEntry3D
     accessors
     axis metric extraction/validation
     coefficient column entries
     cartesian_contracted_parent_metric_packet
     cartesian_contracted_parent_metric_packet_dense_reference

   src/cartesian_contracted_parent_metrics/source_box_route_shadow.jl
     product/product, PQS/product, PQS/PQS source-box shadows
     product slab/source-box fixtures
     density-density source-box blocks
     local Gaussian/nuclear source-box helpers

   src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl
     _pqs_current_route_* inventory
     source-shell/source-mode inventory
     metadata export contracts/writers

   src/cartesian_contracted_parent_metrics/component_smoke_sidecars.jl
     component smoke reports
     CR2 sidecar schemas
     component-report writers/adapters
   ```

   If the exact boundaries differ after reading the code, choose better names,
   but the split must make the metric core obvious and keep legacy/private
   route-shadow code out of `core.jl`.

3. Keep behavior unchanged.

   - Do not change public exports.
   - Do not rename private functions in this pass unless the old names are kept
     and all callers are updated.
   - Do not move code into a different top-level module.
   - Do not add new algorithms.

4. Build a caller inventory while splitting:

   ```text
   exported core functions and their tests
   private source-box route-shadow entry points and callers
   private current-route metadata/export entry points and callers
   private component-smoke/CR2 sidecar entry points and callers
   ```

5. Deletion/shrinkage:

   The line budget remains mandatory. A pure split is not enough if it is
   line-neutral or line-positive.

   You may delete or shrink code/tests only after proving they are stale. Good
   audit questions:

   - Are any component-smoke/CR2 sidecar helpers only called by
     `pqs_component_route_report_adapter_runtests.jl` and not by current driver
     artifacts or CR2 handoff work?
   - Are any source metadata export compatibility columns/tests preserving old
     route-shadow vocabulary rather than a live artifact contract?
   - Are any PQS product/source-box shadow fixtures superseded by
     `CartesianCPBBlockProviders` or current driver endpoints?

   Candidate files to inspect for shrinkage, not blind deletion:

   ```text
   test/nested/pqs_component_route_report_adapter_runtests.jl
   test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl
   test/nested/pqs_source_metadata_real_artifact_acceptance_support.jl
   ```

   Do not delete:

   - `test/nested/cartesian_contracted_parent_metric_packet_runtests.jl`
   - He/H2 physics endpoint tests
   - WL/H/H2 reference/acceptance tests
   - current physical H2 target inventory tests

6. If you cannot keep the source/test/bin budget line-negative without deleting
   live endpoint/reference coverage, stop and write `ATTENTION.md`. Do not
   weaken the line budget.

Line budget:

```text
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
sum(deleted) > sum(added)
```

Docs/blurb logs do not count toward the budget.

Validation:

- Run `julia --project=. -e 'using GaussletBases; println("load ok")'`.
- Run `julia --project=. test/nested/cartesian_contracted_parent_metric_packet_runtests.jl`.
- If you move current-route metadata/export code, run the smallest focused
  source-metadata acceptance test that covers the moved code.
- If you move component-smoke/sidecar code and keep it live, run the smallest
  focused test that covers it.
- Run `git diff --check`.
- Run deleted-helper/test-name searches for any deletion.

Reporting:

Report:

```text
core include line count:
source_box_route_shadow include line count:
current_route_metadata_export include line count:
component_smoke_sidecars include line count:
remaining wrapper line count:
source/test/bin lines added:
source/test/bin lines deleted:
net:
behavior changes: yes/no
tests run:
```

Include deletion/shrinkage report:

```text
deleted:
simplified:
quarantined:
not deleted because:
exact remaining caller/blocker:
```

Stop conditions:

- If the split would require numerical edits, stop and write `ATTENTION.md`.
- If the file’s private functions are too entangled to split mechanically in one
  pass, stop and write `ATTENTION.md` with the smallest safe split.
- If line-negative cannot be achieved safely, stop and write `ATTENTION.md`.
- Do not request UI escalation; write `ATTENTION.md` if approval would be
  required.

When complete:

- This pass is intended to be a stopping point. After response/review/commit,
  the manager should stop the loop for discussion rather than automatically
  issuing the next H2 or cleanup blurb.

-- repo-manager@macmini
