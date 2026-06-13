Pass 192 complete.

No H2 implementation was done. I audited the current visible driver/PQS route
surfaces, deleted one stale development scaffold test, and ran only the focused
validation requested by the blurb.

Driver expressibility audit:

- system/geometry inputs: current `bin/cartesian_ham_builder.jl` can set
  `atom_symbols`, `nuclear_charges`, and `atom_locations`, so H2 centered at
  `(0,0,-2)` and `(0,0,2)` is expressible at the raw coordinate level.
- parent mapping/extents inputs: current driver has `radius`,
  `parent_axis_counts`, `parent_axis_family`, `parent_mapping_rule`,
  `parent_mapping_Z`, `parent_mapping_d`, and `parent_mapping_tail_spacing`.
  It does not expose the old H2 reference controls as first-class inputs:
  `bond_axis`, `bond_length`, `xmax_parallel`, or `xmax_transverse`. Matching
  the documented H2 WL/QW reference geometry would therefore require either a
  new visible input translation layer or manual axis-count/geometry reasoning.
- route_family/route_kind support: `route_family = :pqs_source_box` is already
  visible. The current He driver input uses
  `route_kind = :one_center_fixed_q_complete_core_shell`. I found diatomic
  complete core/shell route-owned payloads in
  `src/pqs_source_box_diatomic_complete_core_shell.jl`, but no existing H2
  driver input file and no supplemented H2 route-kind contract ready for a
  physics endpoint.
- diatomic shellification policy support: route-configured bond-aligned
  diatomic shellification/materializer metadata exists, and the assembly builds
  `diatomic_complete_core_shell_*` payloads. The driver does not yet expose an
  old-standard complete-rectangular H2 policy knob that cleanly says
  "match the pass-191 WL/QW complete-rectangular H2 reference".
- fixed-q/source-mode support: `q`, `n_s`, and fixed PQS source-mode style are
  present for the current PQS path. The diatomic source-plan payloads can
  materialize a source-box-first plan when the parent/route skeleton is ready.
- final-basis support: diatomic complete core/shell final-basis payload support
  exists and records final dimension/support summaries.
- H1 support: diatomic complete core/shell H1 payload support exists. It builds
  support kinetic, support electron-nuclear-by-center matrices, final one-body
  matrices, final H1, and lowest-energy summary when inputs are available.
- H1-J/density-interaction support: diatomic ham-input/handoff payload support
  exists for private diagnostic density-interaction inputs
  (`pre_final_density_interaction`, density gauge, raw numerator convention).
  It is not the same report/artifact-facing He H1-J diagnostic scalar path, and
  it is marked private/non-public/non-export in metadata/readiness.
- private RHF support: current private RHF wiring consumes the one-center
  `complete_core_shell_diagnostic_route_payload`, not the diatomic
  complete-core/shell handoff payload. For H2, explicit electron count could be
  passed, but the route-owned diatomic private-RHF endpoint is not wired as the
  driver artifact authority.
- save/artifact support: `cartesian_save` writes the durable report and, when
  `complete_core_shell_h1_j_driver_route_materialized`, writes the current He
  artifact groups (`config`, `basis`, `physics`, `density_interaction`,
  `comparison`, optional `private_rhf`). It does not write an H2-specific
  artifact, supplement status, residual count, or diatomic private-RHF group.

Supplement/residual policy status:

- old WL/QW reference path supports supplement: yes. The pass-191 references
  use `legacy_bond_aligned_diatomic_gaussian_supplement("H", "cc-pVTZ", nuclei;
  lmax = 1)` and `ordinary_cartesian_qiu_white_operators(...)`, with
  `interaction_treatment = :ggt_nearest`.
- current visible driver supports supplement: no. I found no visible driver
  input keys for `supplement`, `supplement_basis`, `supplement_lmax`,
  `legacy_bond_aligned_diatomic_gaussian_supplement`, residual owner metadata,
  or residual count.
- PQS route-owned source-box path supports supplement: no supported route-owned
  supplement/residual integration is visible in the PQS source-box complete
  core/shell path. The diatomic PQS path is gausslet/source-box-first and uses
  support/final-basis density-interaction diagnostics, not the old WL/QW
  residual Gaussian supplement policy.
- artifact/save path can record supplement status: not currently. The He
  artifact writer records route/mapping/q/n_s, final-basis summary, H1/H1-J,
  density gauge, WL comparison constants, and optional private RHF, but no
  supplement or residual fields.

Recommendation:

Option B: implement gausslet-only H2 PQS first, with `comparison_ready=false`
and no direct HF comparison to the old WL/QW supplemented total.

Reason: the current driver can express a two-center atom list and the PQS
diatomic route has private source-plan/final-basis/H1/ham-handoff pieces, but
the documented H2 HF/ED references include the H/cc-pVTZ S/P residual
supplement and the current visible PQS driver path cannot carry that same
supplement/residual policy. Option A is blocked by missing supplement support.
Option C is stronger than needed because a gausslet-only H2 diagnostic can be
made honest by disabling supplemented WL HF comparison.

Future gausslet-only H2 driver input sketch, if manager approves:

```julia
route_family = :pqs_source_box
route_kind = :bond_aligned_diatomic_fixed_q_complete_core_shell

atom_symbols = ("H", "H")
nuclear_charges = (1, 1)
atom_locations = ((0.0, 0.0, -2.0), (0.0, 0.0, 2.0))
bond_axis = :z
bond_length = 4.0

parent_axis_family = :G10
core_spacing = 0.5
xmax_parallel = 6.0
xmax_transverse = 4.0

q = 5
n_s = 5
fixed_source_mode_shape = true

supplement_policy = :none
comparison_ready = false
comparison_reference_label = "gausslet-only H2 PQS R=4.0"
wl_rhf_total = nothing

run_h1 = true
run_h1_j = true
run_private_rhf = false
save_artifact = true
save_tsv = true
```

For a future supplemented target, the input would need explicit supplement
keys before comparing to the pass-191 HF total, for example
`supplement_policy = :legacy_bond_aligned_diatomic_gaussian`,
`supplement_basis = "cc-pVTZ"`, and `supplement_lmax = 1`. Those keys are not
supported by the current visible driver/PQS path.

Artifact shape note:

- Future H2 artifact should record final dimension, support/source/final-basis
  summaries, H1 lowest energy when available, density-interaction gauge/raw
  convention, nuclear repulsion/electron count/spin sector, and
  `comparison_ready = false` for gausslet-only.
- Final-basis self-overlap should remain a scalar identity-error diagnostic
  only. Do not save or consume an `S` matrix as downstream working data for an
  orthonormal final basis.
- Supplement fields should be explicit before any WL/QW HF total comparison:
  supplement policy, basis, lmax, residual count, residual owner set, and
  interaction treatment.

Stale scaffold deleted:

- Deleted `test/nested/cartesian_pair_block_one_body_block_set_summary_runtests.jl`.
- Safety: it is not included in `test/nested/runtests.jl`,
  `test/nested/integration_runtests.jl`, or top-level runners; it is not a He
  endpoint, H/H2 WL acceptance/reference, or scientific endpoint; it is not
  needed for this H2 audit. Its role was a metadata-only block-set summary shape
  scaffold for future mixed one-body term-set consumption, and the preferred
  mixed one-body consumer smoke/accessor tests already cover the active compact
  summary path.

Validation results:

- `rg -n "cartesian_pair_block_one_body_block_set_summary_runtests|CartesianPairBlockMaterialization one-body block-set summary shape" test src bin`
  - passed: no hits
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check`
  - passed
- `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - `0	237	test/nested/cartesian_pair_block_one_body_block_set_summary_runtests.jl`
- `git status --short --branch`
  - `## main...origin/main`
  - ` D test/nested/cartesian_pair_block_one_body_block_set_summary_runtests.jl`

Deletion/shrinkage report:

- deleted: `test/nested/cartesian_pair_block_one_body_block_set_summary_runtests.jl`
- simplified: mixed one-body block-set summary carrying cost reduced by 237
  development-scaffold lines
- quarantined: nothing
- not deleted because: active He driver endpoint tests, WL H/H2 reference tests,
  and included runner tests remain live
- exact remaining caller/blocker: H2 PQS endpoint implementation remains blocked
  on a manager-approved choice to either add explicit driver/PQS supplement
  support or proceed gausslet-only with `comparison_ready=false`; old
  supplemented WL/QW HF totals must not be used as direct comparisons until the
  H/cc-pVTZ S/P residual policy is represented in the driver artifact.

-- repo-doer@macmini
