Pass 176 response - WL generator subphase timing and stale pilot deletion

Files changed/deleted:
- Modified `tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
- Deleted `test/nested/cartesian_cpb_overlap_placement_pilot_runtests.jl`

What changed:
- Added compact Julia `@elapsed` prints inside `_be2_white_lindsey_atom_growth_route()` for the existing sequential WL construction calls:
  - `generator.phase.wl.system`
  - `generator.phase.wl.recipe`
  - `generator.phase.wl.parent`
  - `generator.phase.wl.shells`
  - `generator.phase.wl.units`
  - `generator.phase.wl.transforms`
  - `generator.phase.wl.pairs`
  - `generator.phase.wl.assembly`
  - `generator.phase.wl.report`
  - `generator.phase.wl.materialization`
- Did not change fixture size, route semantics, H1/Vee values, overlap storage policy, solver flags, or artifact schema.
- Deleted the stale overlap-placement pilot test after caller inspection found no live source/default-runner reference.

Line-budget result:

```text
git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
0       719     test/nested/cartesian_cpb_overlap_placement_pilot_runtests.jl
25      10      tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Total added: 25
Total deleted: 729
Net: -704 lines

WL subphase timing output:

```text
generator.phase.wl.system=0.0
generator.phase.wl.recipe=0.0
generator.phase.wl.parent=0.559229875
generator.phase.wl.shells=8.628097625
generator.phase.wl.units=44.913618875
generator.phase.wl.transforms=49.239544334
generator.phase.wl.pairs=13.697098334
generator.phase.wl.assembly=15.020014208
generator.phase.wl.report=12.690385917
generator.phase.wl.materialization=39.619359542
generator.phase.wl_route_build=184.368656083
generator.phase.wl_population=0.17125825
generator.phase.bundle_write=0.673834833
generator.phase.total=227.951998625
```

Largest repeated WL route-build costs in this run:
- `transforms`: 49.24 s
- `units`: 44.91 s
- `materialization`: 39.62 s
- then `assembly`, `pairs`, and `report` at 12-15 s each.

Artifact size:

```text
-rw-r--r--@ 1 srw  staff    82M Jun 13 11:30 tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_inspection_bundle.jld2
-rw-r--r--@ 1 srw  staff   884B Jun 13 11:30 tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_fingerprint.tsv
```

Readback:

```text
producer_commit=87d22e380c9e5bcb00858beb23b9afdb5e3981af
producer_dirty=true
pqs_h1=(221, 221)
pqs_vee=(221, 221)
wl_h1=(2287, 2287)
wl_vee=(2287, 2287)
wl_solver_ready=false
```

Generated outputs ignored status:

```text
## main...origin/main
 M tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
!! tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_fingerprint.tsv
!! tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_inspection_bundle.jld2
```

Validation commands/results:

```text
rg -n "cartesian_cpb_overlap_placement_pilot_runtests" test src docs/src/developer
```

Result: only blurb-log mentions; no live source/default-runner caller.

```text
julia --project=. tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Result: passed; bundle regenerated; `cr2_solver_ready=false`; `white_lindsey_status=available_route_configured_diatomic_atom_growth_ham_payload`.

```text
julia --project=. -e 'using JLD2; p="tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_inspection_bundle.jld2"; jldopen(p,"r") do f; println("producer_commit=", f["producer/repo_commit"]); println("producer_dirty=", f["producer/dirty"]); println("pqs_h1=", size(f["routes/pqs_source_box/one_body/hamiltonian"])); println("pqs_vee=", size(f["routes/pqs_source_box/two_body/interaction_matrix"])); println("wl_h1=", size(f["routes/white_lindsey/one_body/hamiltonian"])); println("wl_vee=", size(f["routes/white_lindsey/two_body/interaction_matrix"])); println("wl_solver_ready=", f["routes/white_lindsey/readiness/cr2_solver_ready"]); end'
```

Result: passed with readback shown above.

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Result:

```text
load ok
```

```text
git diff --check
```

Result: passed with no output.

```text
git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Result: net-negative by 704 lines, shown above.

```text
git status --short --branch --ignored tmp/work/be2_wl_pqs_cr2_inspection_artifact
```

Result: generated JLD2/TSV outputs remain ignored, shown above.

Git status:

```text
## main...origin/main
 D test/nested/cartesian_cpb_overlap_placement_pilot_runtests.jl
 M tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Deletion/shrinkage report:
- deleted: `test/nested/cartesian_cpb_overlap_placement_pilot_runtests.jl`
- simplified: WL route-build observability now lives directly at the existing generator seam without a new report file or dependency.
- quarantined: none.
- not deleted because: the active tracked generator remains needed to regenerate the Be2 WL/PQS CR2 inspection bundle.
- exact remaining caller/blocker: no live caller for the deleted pilot test; next optimization blocker is measured WL route-build cost concentrated in `cartesian_transforms`, `cartesian_units`, and materialization inside `_be2_white_lindsey_atom_growth_route()`.

-- repo-doer@macmini
