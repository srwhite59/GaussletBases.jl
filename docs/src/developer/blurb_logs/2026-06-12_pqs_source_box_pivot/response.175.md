Pass 175 WL HF-minimal one-body handoff shrink complete.

Commit:

- `6c251ea2 Shrink WL CR2 one-body handoff`

Files changed:

- `tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`

Tracked implementation line budget:

```text
3	8	tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

- Added: 3
- Deleted: 8
- Net: -5
- Requirement satisfied: deleted > added.

Artifact size:

- Before: `202M`
- After: `82M`
- Shrink: about `120M`

Generator timing output:

```text
generator.phase.pqs_route_build=45.027609625
generator.phase.pqs_payload=0.3788265
generator.phase.provenance=0.307479583
generator.phase.wl_route_build=186.269930834
generator.phase.wl_population=0.169682791
generator.phase.bundle_write=0.687189083
generator.phase.total=232.861505792
```

Interpretation:

- This pass improved artifact size/read cost, not route build time.
- Runtime remains dominated by route construction, especially
  `wl_route_build`.
- Payload population and JLD2 writing remain subsecond.

Artifact readback from the pre-commit validation run:

```text
producer_commit=59ed67d253c31d7648f42a5f55d919ba2895e751
producer_dirty=true
wl_h1=(2287, 2287)
wl_vee=(2287, 2287)
has_wl_kinetic=false
has_wl_nuclear_by_center=false
wl_kinetic_stored=false
wl_nuclear_by_center_stored=false
wl_solver_ready=false
```

The artifact was not regenerated after committing this pass because another
generator optimization pass is expected before the next CR2 handoff.

Validation commands/results:

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

- deleted: dense WL `kinetic_one_body` and `nuclear_one_body_by_center` sidecars from the generated artifact payload
- simplified: WL route now hands CR2 only final-basis `H1` plus final-basis density-density `Vee` for the first HF-facing matrix contract
- quarantined: one-body split data remains route-internal/debug data, not downstream CR2 working data
- not deleted because: WL route still constructs the full operator payload internally; generator-only shrink does not change route construction
- exact remaining caller/blocker: generator runtime is still dominated by `wl_route_build`; next pass should target avoiding unnecessary WL construction work or adding a faster HF-minimal WL route build, not JLD2 writes

-- repo-manager@macmini
