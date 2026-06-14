Pass 202 - shrink the driver helper by removing or extracting private global overlap.

Target card:

- Cleanup target: `src/pqs_source_box_route_driver_helpers.jl`.
- Physics target protected: physical H2 PQS route development should not be
  blocked by old private overlap scaffolding in the driver helper.
- Allowed implementation surface: audit, deletion or extraction of the
  `private_global_overlap` subtree, minimal include/wrapper changes, and
  deletion/shrinkage of tests that only preserve this stale private surface.
- Forbidden surfaces: H2 source-plan implementation, final basis, H1, H1-J,
  density interaction, RHF, WL comparison, CR2/export/artifact behavior, public
  APIs, broad report alias expansion, unrelated file splitting.
- Success condition: `pqs_source_box_route_driver_helpers.jl` no longer owns the
  private global-overlap implementation. It either calls one compact extracted
  stage function, or the path is deleted if no live target/caller justifies it.

Why this pass:

- `pqs_source_box_route_driver_helpers.jl` is still about 13k lines.
- The `private_global_overlap` path is private, optional, default-off, and
  currently lives in the human-facing driver helper even though it owns overlap
  facts, local CPB collection adapters, placement requirements, placement
  candidates, placement-plan skeletons, reviewed placement plans, and result
  summaries.
- That is a placement/global-assembly sub-pipeline, not top-level driver
  orchestration.

Read first:

- `AGENTS.md`, especially code cleanup, test deletion bias, and line-budget
  discipline.
- `JuliaStyle.md` for module/file ownership and compact structured concepts.
- `docs/code_bloat_and_wrong_contract_cleanup_note.md`.
- Recent pass 201 review so you do not accidentally resume H2 implementation.

Audit task:

1. Inventory all live `private_global_overlap` references:

   ```text
   rg -n "private_global_overlap" src test bin examples docs
   ```

2. Classify the path:

   ```text
   live current target:
     yes / no

   production/public/default behavior:
     yes / no

   driver option default:
     on / off

   tests that are endpoint/reference:
     list

   tests that are private scaffold:
     list

   exact source callers:
     list
   ```

3. Decide:

   - If the path is not a live target and only private scaffold tests preserve
     it, delete it.
   - If it is still live, extract it intact to a narrower ownership file and
     leave only a compact wrapper call in the driver helper.

Preferred deletion path if stale:

- Remove the private global-overlap stage and helper family from
  `src/pqs_source_box_route_driver_helpers.jl`.
- Remove default-off driver input knobs from `bin/cartesian_ham_builder.jl` only
  if they have no live external use:

  ```text
  private_global_overlap_requested
  private_global_overlap_global_dimension
  private_global_overlap_inputs
  ```

- Remove print/report/materialization fields that only exist for the deleted
  private path.
- Delete private-global-overlap tests that only preserve this path. Candidate
  files to audit carefully:

  ```text
  test/nested/cartesian_cpb_local_overlap_fingerprint_runtests.jl
  test/nested/cartesian_pair_block_driver_global_overlap_runtests.jl
  examples/private_global_overlap_option.jl
  ```

  Do not delete a real endpoint/reference test. If one of these files also
  carries still-live CPB local-overlap provider coverage, split or preserve the
  live part only if that is smaller than retaining the stale surface.

Preferred extraction path if live:

- Create one new private file, for example:

  ```text
  src/cartesian_route_private_global_overlap.jl
  ```

- Move all private-global-overlap implementation helpers there.
- Include it from `src/GaussletBases.jl` before
  `pqs_source_box_route_driver_helpers.jl` if needed.
- Leave only one compact call in `pqs_source_box_route_driver_helpers.jl`, e.g.

  ```julia
  private_global_overlap_stage =
      _pqs_source_box_route_driver_private_global_overlap_stage(report; ...)
  ```

- Do not add new behavior or new report fields. This should be behavior
  preserving except for file ownership.

Line budget:

- The source/test/bin/generator line budget is still mandatory:

  ```text
  git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
  sum(deleted) > sum(added)
  ```

- If extraction is necessary and is line-neutral or line-positive, offset it by
  deleting stale private-global-overlap tests or stale mixed one-body scaffolds.
  Do not use He/H2 physics endpoint tests or WL/H2 reference tests as budget.

Validation:

- Run `julia --project=. -e 'using GaussletBases; println("load ok")'`.
- Run `git diff --check`.
- If deleting the path, run live searches for deleted helper names and test file
  names.
- If extracting but preserving behavior, run the smallest private-global-overlap
  focused test that proves the wrapper still works. Do not run a broad route
  suite unless the audit shows it is necessary.
- If deleting driver option knobs, run a focused driver smoke that proves
  ordinary driver defaults still load/run.

Reporting:

- State whether the path was deleted or extracted.
- Report:

  ```text
  private_global_overlap references before/after:
  lines removed from pqs_source_box_route_driver_helpers.jl:
  source/test/bin lines added:
  source/test/bin lines deleted:
  net:
  tests deleted:
  tests preserved and why:
  validation run:
  ```

- Include deletion/shrinkage report:

  ```text
  deleted:
  simplified:
  quarantined:
  not deleted because:
  exact remaining caller/blocker:
  ```

Stop conditions:

- If a live current physics target depends on private global overlap, do not
  delete it. Extract it or write `ATTENTION.md` if extraction is unsafe.
- If deleting it would alter the H2 physical target inventory path from pass
  201, stop and write `ATTENTION.md`.
- If this cleanup cannot be kept line-negative without deleting accepted
  endpoint/reference tests, stop and write `ATTENTION.md`.
- Do not request UI escalation; write `ATTENTION.md` if a required command would
  need approval.

-- repo-manager@macmini
