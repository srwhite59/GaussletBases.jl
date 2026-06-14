Pass 220 response - H2 WL gausslet-only 463 reference value probe

Probe path:
- `tmp/work/h2_wl_gausslet_only_463_reference_probe.jl`
- `tmp/work/h2_wl_gausslet_only_463_reference_probe_summary.txt`
- Both are ignored local `tmp/work` files.

Exact commands run:
```sh
julia --project=. tmp/work/h2_wl_gausslet_only_463_reference_probe.jl
git status --short --branch
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Notes on local probe corrections:
- First run failed before heavy work because the local `timed!` helper was missing the `do`-block method.
- Second run built the fixed block/operators, then failed while formatting inventory because `BondAlignedDiatomicQWBasis3D` has no `size(...)` method.
- I fixed only those local ignored probe-script issues and reran successfully.

Final dimension and inventory:
- final dimension 463 reproduced: yes
- support/retained inventory matched candidate: yes
- geometry: H2 R=4 with nuclei `[(0.0, 0.0, -2.0), (0.0, 0.0, 2.0)]`
- parent axis counts: `(9, 9, 15)`
- support counts: `(275, 578, 362)`
- retained counts: `(251, 98, 114)`
- retained order: `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`
- fixed block size: `(1215, 463)`
- residual count: `0`
- supplement policy: `:none`
- retained transform kind: `:white_lindsey_old_qw_gausslet_retained_transform`
- shared layer types: `(:_CartesianNestedCompleteShell3D, :_CartesianNestedCompleteShell3D)`
- final overlap identity error: `1.3816121305989328e-14`

WL gausslet-only values obtained:
- H1 lowest energy: `-0.7946609179724673`
- H1 orbital self-Coulomb: `0.45696639804337047`
- RHF converged: `true`
- RHF iterations: `9`
- RHF one-electron energy: `-1.5611571934181985`
- RHF electron-electron energy: `0.40220533775308426`
- RHF electronic energy: `-1.1589518556651142`
- RHF nuclear repulsion: `0.25`
- RHF total with nuclear repulsion: `-0.9089518556651142`
- RHF density trace: `0.9999999999999952`
- RHF electron count: `1.9999999999999905`
- RHF orbital residual: `1.5840620684318253e-14`

Provenance:
- constructor: `bond_aligned_homonuclear_qw_basis(R=4.0, core_spacing=0.5, xmax_parallel=6.0, xmax_transverse=4.0, bond_axis=:z)`
- fixed block source: old nested/QW complete-rectangular source via `_nested_bond_aligned_diatomic_source(...)` and `_nested_fixed_block(...)`
- operator source: `ordinary_cartesian_qiu_white_operators(fixed_block; supplement omitted)`
- interaction treatment: `:ggt_nearest`
- gausslet backend: `:numerical_reference`
- no GTO/MWG supplement was constructed or passed.

Timing by phase:
```text
coulomb expansion                         0.018949916 s
bond-aligned basis                        0.022052125 s
complete-rectangular source and fixed block 10.554786666 s
gausslet-only ordinary QW operators       4.793927000 s
H1 lowest solve                           0.035003208 s
H1 orbital self-Coulomb                   0.000598916 s
restricted closed-shell HF                0.463980791 s
```

Blocker:
- Values were obtained, so the previous blocker `:missing_wl_h2_gausslet_only_reference_values` can be retired after review.
- Remaining review/implementation blocker before driver adoption: `:wl_h2_gausslet_only_reference_values_not_yet_wired_to_driver_artifact`.

Recommended next pass:
- Add the reviewed WL no-supplement H2 463 scalar values to the H2 physical driver artifact/report as compact comparison fields.
- Keep the visible driver shape.
- Do not compare to supplemented WL/QW scalar values.
- Assert the PQS-vs-WL deltas using:
  - H1 lowest energy,
  - H1 orbital self-Coulomb,
  - RHF electronic energy and/or a clearly labeled nuclear-repulsion-included total.

Tracked source/test/bin line budget:
- No tracked `src`, `test`, or `bin` files changed.
- Scoped `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl` produced no output.

Git status:
```text
## main...origin/main
```

Deletion/shrinkage report:
- deleted: none
- simplified: none
- quarantined: old supplemented WL/QW H2 scalar references remain quarantined as non-comparable for the no-supplement endpoint
- not deleted because: audit/probe-only pass; no tracked source/test/bin changes requested
- exact remaining caller/blocker: `:wl_h2_gausslet_only_reference_values_not_yet_wired_to_driver_artifact`

-- repo-doer@macmini
