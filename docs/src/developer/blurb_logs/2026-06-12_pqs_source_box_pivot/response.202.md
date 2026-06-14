Pass 202 response - private global-overlap driver cleanup

Path decision:
- Deleted, not extracted.
- Classification:
  - live current target: no
  - production/public/default behavior: no
  - driver option default: off before deletion
  - endpoint/reference tests: none for the private driver hook
  - private scaffold tests: `test/nested/cartesian_cpb_local_overlap_fingerprint_runtests.jl`, `test/nested/cartesian_pair_block_driver_global_overlap_runtests.jl`
  - exact source callers before deletion: `bin/cartesian_ham_builder.jl` materialization input knobs and `_pqs_source_box_route_driver_materialization` in `src/pqs_source_box_route_driver_helpers.jl`

Implementation:
- Removed the private global-overlap helper subtree from `src/pqs_source_box_route_driver_helpers.jl`.
- Removed private global-overlap keyword plumbing and result/summary spreading from `_pqs_source_box_route_driver_materialization`.
- Removed default-off driver knobs from `bin/cartesian_ham_builder.jl`:
  - `private_global_overlap_requested`
  - `private_global_overlap_global_dimension`
  - `private_global_overlap_inputs`
- Deleted `examples/private_global_overlap_option.jl`.
- Deleted private scaffold tests:
  - `test/nested/cartesian_cpb_local_overlap_fingerprint_runtests.jl`
  - `test/nested/cartesian_pair_block_driver_global_overlap_runtests.jl`
- Updated two older docs notes so they describe the private driver bridge as retired rather than current.

private_global_overlap references before/after:
- Before: `rg -n "private_global_overlap" src test bin examples docs` found hits in `src`, `test`, `bin`, `examples`, two older docs notes, and the active pass-202 blurb log.
- After: no hits in `src test bin examples`; only the active pass-202 blurb log still contains the exact `private_global_overlap` string.

Line budget:
- lines removed from `pqs_source_box_route_driver_helpers.jl`: 1514
- source/test/bin lines added: 0
- source/test/bin lines deleted: 4510
- net source/test/bin: -4510
- overall diff stat: 23 insertions, 4564 deletions

Tests preserved and why:
- Preserved route-global one-body oracle/reference tests such as the existing global overlap/kinetic/position/x2 tests; they exercise the lower route-global one-body materialization layer, not the deleted private driver hook.
- Preserved CPB/provider and retained-unit tests that do not call the deleted private driver helper names.

Validation run:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed; precompiled `GaussletBases`, then printed `load ok`
- `julia --project=. bin/cartesian_ham_builder.jl save_artifact=false save_tsv=false`
  - passed; printed `driver complete`
- `git diff --check`
  - passed
- `rg -n "private_global_overlap" src test bin examples`
  - no matches, exit 1 expected
- `rg -n "cartesian_cpb_local_overlap_fingerprint_runtests|cartesian_pair_block_driver_global_overlap_runtests|private_global_overlap_option|_pqs_source_box_route_driver_private_global_overlap" src test bin examples`
  - no matches, exit 1 expected
- `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - `0 5 bin/cartesian_ham_builder.jl`
  - `0 1514 src/pqs_source_box_route_driver_helpers.jl`
  - `0 705 test/nested/cartesian_cpb_local_overlap_fingerprint_runtests.jl`
  - `0 2286 test/nested/cartesian_pair_block_driver_global_overlap_runtests.jl`

Git status:
```text
## main...origin/main
 M bin/cartesian_ham_builder.jl
 M docs/src/developer/cartesian_parent_factors_and_cpb_kernels.md
 M docs/src/developer/cartesian_route_retirement_ledger.md
 D examples/private_global_overlap_option.jl
 M src/pqs_source_box_route_driver_helpers.jl
 D test/nested/cartesian_cpb_local_overlap_fingerprint_runtests.jl
 D test/nested/cartesian_pair_block_driver_global_overlap_runtests.jl
```

Deletion/shrinkage report:
- deleted: private global-overlap driver helper subtree; driver option knobs; option example; two private scaffold tests.
- simplified: `_pqs_source_box_route_driver_materialization` no longer coordinates private global-overlap result/summary fields.
- quarantined: lower route-global one-body oracle/reference paths remain outside the driver helper; older docs now label the private driver bridge as retired.
- not deleted because: lower route-global one-body materialization and oracle/reference tests are still live replacement pressure below the driver layer.
- exact remaining caller/blocker: no remaining `src/test/bin/examples` caller of the deleted private driver hook; H2 physical target remains blocked on its separate missing source-plan/final-basis/H1/H1-J work, not on private global overlap.

-- repo-doer@macmini
