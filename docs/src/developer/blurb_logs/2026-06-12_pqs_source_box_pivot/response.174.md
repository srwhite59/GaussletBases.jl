Pass 174 generator timing and scaffold deletion complete.

Commit:

- `c84b7a03 Instrument CR2 generator phases`

Files changed/deleted:

- modified `tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
- deleted `test/nested/cartesian_cpb_overlap_placement_facts_runtests.jl`

Tracked implementation line budget:

```text
0	711	test/nested/cartesian_cpb_overlap_placement_facts_runtests.jl
24	11	tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

- Added: 24
- Deleted: 722
- Net: -698
- Requirement satisfied: deleted > added.

Generator timing output:

```text
generator.phase.pqs_route_build=44.31484975
generator.phase.pqs_payload=0.380829417
generator.phase.provenance=0.308937
generator.phase.wl_route_build=186.2394435
generator.phase.wl_population=0.187260042
generator.phase.bundle_write=0.786514
generator.phase.total=232.239269167
```

Interpretation:

- The bottleneck is route construction, especially `wl_route_build`.
- Payload extraction, provenance labeling, WL payload population, and JLD2 write
  are all subsecond.
- Optimizing bundle writing or JLD2 storage will not fix the main runtime.

Artifact readback from the pre-commit validation run:

```text
producer_commit=fe4a943a6e79c98e95c82e9e19e6a77599cea1ee
producer_dirty=true
pqs_h1=(221, 221)
pqs_vee=(221, 221)
wl_h1=(2287, 2287)
wl_vee=(2287, 2287)
pqs_solver_ready=false
wl_solver_ready=false
```

The artifact was not regenerated after committing this pass because another
generator optimization pass is expected before the next CR2 handoff.

Validation commands/results:

- `rg -n "cartesian_cpb_overlap_placement_facts_runtests" test src docs/src/developer`
  - Found only old blurb-log references, no live source/default-runner caller.
- `julia --project=. tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - Passed; timing output above.
- JLD2 readback command passed; values summarized above.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - Passed: `load ok`.
- `git diff --check`
  - Passed before commit.
- `git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - Passed line-budget rule; details above.

Generated outputs ignored status:

```text
!! tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_fingerprint.tsv
!! tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_inspection_bundle.jld2
```

Deletion/shrinkage report:

- deleted: stale metadata-only `cartesian_cpb_overlap_placement_facts_runtests.jl`
- simplified: generator now reports top-level route/payload/write timing without a separate timing artifact
- quarantined: remaining CPB overlap placement pilot test was left for a later deletion decision
- not deleted because: the pilot test exercises materialized placement behavior, while this pass only needed to retire the facts-only scaffold
- exact remaining caller/blocker: generator runtime is dominated by `wl_route_build`; next pass should reduce or bypass duplicated WL route construction work, not tune JLD2 writing

-- repo-manager@macmini
