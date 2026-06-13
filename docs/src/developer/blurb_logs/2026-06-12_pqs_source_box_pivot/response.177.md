Pass 177 response - WL/PQS comparison-readiness gate and WL Z audit labels

Files changed:
- Modified `tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`

What changed:
- Added a compact generated artifact group at `comparison/wl_pqs`:
  - `pqs_final_dimension`
  - `wl_final_dimension`
  - `final_dimension_match`
  - `comparison_ready`
  - `comparison_blocker`
  - `comparison_role`
- Added WL route-local metadata fields:
  - `routes/white_lindsey/metadata/white_lindsey_Z`
  - `routes/white_lindsey/metadata/nuclear_charge_tuple`
  - `routes/white_lindsey/metadata/white_lindsey_Z_audit_status`
  - `routes/white_lindsey/metadata/white_lindsey_Z_audit_blocker`
- Kept primary route readiness keys intact:
  - `routes/pqs_source_box/readiness/cr2_read_only_inspector_ready`
  - `routes/white_lindsey/readiness/cr2_read_only_inspector_ready`
  - `routes/pqs_source_box/readiness/cr2_solver_ready`
  - `routes/white_lindsey/readiness/cr2_solver_ready`
- Removed stale WL generator-only readiness aliases and development-era unavailable metadata labels to keep the generator-only diff net-negative.
- Did not change fixture size, route construction, H1/Vee values, physics, solver/export readiness, or overlap storage.

Line-budget result:

```text
git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
13      17      tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Total added: 13
Total deleted: 17
Net: -4 lines

Comparison-readiness readback:

```text
pqs_dim=221
wl_dim=2287
dim_match=false
comparison_ready=false
comparison_blocker=wl_pqs_final_dimension_mismatch
comparison_role=separate_route_inspection_not_same_basis_comparison
```

WL Z audit readback:

```text
wl_Z=2.0
wl_Z_audit_status=requires_review
wl_Z_audit_blocker=white_lindsey_Z_differs_from_nuclear_charge
```

H1/Vee dimensions and route readiness flags:

```text
pqs_h1=(221, 221)
pqs_vee=(221, 221)
wl_h1=(2287, 2287)
wl_vee=(2287, 2287)
pqs_readonly=true
wl_readonly=true
pqs_solver_ready=false
wl_solver_ready=false
```

Artifact size:

```text
-rw-r--r--@ 1 srw  staff   884B Jun 13 11:41 tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_fingerprint.tsv
-rw-r--r--@ 1 srw  staff    82M Jun 13 11:41 tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_inspection_bundle.jld2
```

Generator timing:

```text
generator.phase.pqs_route_build=41.803258209
generator.phase.pqs_payload=0.36876925
generator.phase.provenance=0.304728917
generator.phase.wl.system=4.1e-8
generator.phase.wl.recipe=4.2e-8
generator.phase.wl.parent=0.552521459
generator.phase.wl.shells=8.500714166
generator.phase.wl.units=44.930928417
generator.phase.wl.transforms=49.467097875
generator.phase.wl.pairs=13.883270667
generator.phase.wl.assembly=14.682398167
generator.phase.wl.report=12.770052125
generator.phase.wl.materialization=39.329469583
generator.phase.wl_route_build=184.117750042
generator.phase.wl_population=0.277398375
generator.phase.bundle_write=0.681244916
generator.phase.total=227.573432417
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
julia --project=. tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Result: passed; bundle regenerated; `cr2_solver_ready=false`; `white_lindsey_status=available_route_configured_diatomic_atom_growth_ham_payload`.

```text
julia --project=. -e 'using JLD2; p="tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_inspection_bundle.jld2"; jldopen(p,"r") do f; println("pqs_dim=", f["comparison/wl_pqs/pqs_final_dimension"]); println("wl_dim=", f["comparison/wl_pqs/wl_final_dimension"]); println("dim_match=", f["comparison/wl_pqs/final_dimension_match"]); println("comparison_ready=", f["comparison/wl_pqs/comparison_ready"]); println("comparison_blocker=", f["comparison/wl_pqs/comparison_blocker"]); println("comparison_role=", f["comparison/wl_pqs/comparison_role"]); println("wl_Z=", f["routes/white_lindsey/metadata/white_lindsey_Z"]); println("wl_Z_audit_status=", f["routes/white_lindsey/metadata/white_lindsey_Z_audit_status"]); println("wl_Z_audit_blocker=", f["routes/white_lindsey/metadata/white_lindsey_Z_audit_blocker"]); println("pqs_readonly=", f["routes/pqs_source_box/readiness/cr2_read_only_inspector_ready"]); println("wl_readonly=", f["routes/white_lindsey/readiness/cr2_read_only_inspector_ready"]); println("pqs_solver_ready=", f["routes/pqs_source_box/readiness/cr2_solver_ready"]); println("wl_solver_ready=", f["routes/white_lindsey/readiness/cr2_solver_ready"]); end'
```

Result: passed with requested readback values shown above.

```text
julia --project=. -e 'using JLD2; p="tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_inspection_bundle.jld2"; jldopen(p,"r") do f; println("pqs_h1=", size(f["routes/pqs_source_box/one_body/hamiltonian"])); println("pqs_vee=", size(f["routes/pqs_source_box/two_body/interaction_matrix"])); println("wl_h1=", size(f["routes/white_lindsey/one_body/hamiltonian"])); println("wl_vee=", size(f["routes/white_lindsey/two_body/interaction_matrix"])); end'
```

Result: passed with dimensions shown above.

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

Result: net-negative by 4 lines, shown above.

```text
git status --short --branch --ignored tmp/work/be2_wl_pqs_cr2_inspection_artifact
```

Result: generated JLD2/TSV outputs remain ignored, shown above.

Git status:

```text
## main...origin/main
 M tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Deletion/shrinkage report:
- deleted: stale generator-only WL readiness aliases (`solver_ready`, `hfdmrg_ready`, `rhf_ready`, `dmrg_ready`, `hamv6_export_ready`) and development-era unavailable metadata labels (`supplement_residual_gto_status`, `qiu_white_atom_local_hf_inputs_status`, `correction_egoi_stationary_cusp_status`, `mwg_ida_route_configured_diatomic_ham_status`, `old_seed_one_center_promoted`) plus three consumed-flag metadata fields.
- simplified: WL/PQS comparison status is now one compact `comparison/wl_pqs` group rather than implied by unrelated route-local dimensions.
- quarantined: WL `white_lindsey_Z=2.0` versus nuclear charges `(4, 4)` is explicitly labeled `requires_review` with blocker `white_lindsey_Z_differs_from_nuclear_charge`.
- not deleted because: primary `cr2_*` readiness keys remain the artifact contract for CR2 read-only inspection and solver/export false claims.
- exact remaining caller/blocker: same-basis WL/PQS comparison remains blocked by `wl_pqs_final_dimension_mismatch` (`221` vs `2287`); WL Z meaning remains blocked by `white_lindsey_Z_differs_from_nuclear_charge`.

-- repo-doer@macmini
