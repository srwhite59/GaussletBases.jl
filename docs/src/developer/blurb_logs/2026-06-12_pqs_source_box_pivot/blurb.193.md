Pass 193 - add driver-visible H2 gausslet-only readiness artifact.

Purpose:

Start the H2 path in the same driver-owned style as the He 419 endpoint, but do
not turn it into a physics comparison yet. Pass 192 established that the old
WL/QW H2 HF/ED references include an H/cc-pVTZ S/P residual supplement that the
current PQS driver path cannot represent. Therefore the first H2 driver target
must be gausslet-only and explicitly not comparable to those supplemented
references.

This pass should make the visible driver honestly express the intended H2
gausslet-only target and write an artifact/readiness record even if downstream
diatomic H1/H1-J/RHF stages are still blocked.

Task type:

Implementation, driver/artifact readiness only.

Physics target:

Future H2 driver-owned gausslet-only PQS diagnostic:

```text
H2 at R = 4.0 bohr
bond axis = :z
atoms at (0,0,-2) and (0,0,2)
nuclear charges = (1,1)
route_family = :pqs_source_box
q = 5
n_s = 5
supplement_policy = :none
comparison_ready = false
```

Important guardrail:

Do not compare this artifact to the old WL/QW H2 HF totals
`-0.910938264352` or `-0.910977315003`. Those are supplemented references.
This pass is about honest driver expression and readiness/blocker reporting,
not physics acceptance.

Read first:

```text
bin/cartesian_ham_builder.jl
test/driver_inputs/he_pqs_q5_wlmap.jl
test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl
test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_rhf_runtests.jl
src/pqs_source_box_diatomic_complete_core_shell.jl
src/pqs_source_box_route_driver_helpers.jl
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.192.md
```

Preserve the visible laboratory-driver style:

- defaults and include/override flow stay visible;
- staged construction in `bin/cartesian_ham_builder.jl` stays visible;
- do not hide this behind an opaque `run_driver(config)` wrapper.

Exact task:

1. Add a small driver input file:

   ```text
   test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl
   ```

   It should set the H2 geometry and target semantics clearly. Suggested keys:

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
   comparison_blocker = :supplemented_reference_not_comparable_to_gausslet_only
   comparison_reference_label = "WL/QW H2 R=4 supplemented reference not used"

   run_h1 = false
   run_h1_j = false
   run_private_rhf = false
   save_artifact = true
   save_tsv = true
   ```

   Adjust names to match existing driver conventions, but keep the meaning.

2. Extend the visible driver/artifact path just enough to record these fields:

   ```text
   system/atom_symbols
   system/nuclear_charges
   system/atom_locations
   system/bond_axis
   system/bond_length

   config/route_family
   config/route_kind
   config/q
   config/n_s
   config/supplement_policy = :none
   config/comparison_ready = false

   comparison/ready = false
   comparison/blocker = :supplemented_reference_not_comparable_to_gausslet_only
   comparison/reference_label
   ```

   Also record whatever diatomic PQS readiness/materialization status the
   current route can honestly expose. If source-plan/final-basis/H1/H1-J/RHF is
   blocked, preserve the exact blocker. Do not invent availability.

3. Add one thin explicit/manual test that runs the driver and reads the artifact.

   The test should assert artifact honesty:

   ```text
   H2 geometry recorded
   q == 5
   n_s == 5
   supplement_policy == :none
   comparison_ready == false
   comparison blocker explains supplemented reference mismatch
   no WL HF total is treated as an accepted comparison
   no private RHF materialized
   route/global/export/public flags remain false if such fields are present
   ```

   Keep this out of default `test/nested/runtests.jl`. It is an explicit driver
   endpoint/readiness test.

4. If the current driver cannot write an artifact when later diatomic stages
   are blocked, implement the smallest blocker-aware save path needed. Do not
   force H1 just to make saving easy.

5. Keep final-basis self-overlap out of downstream working data. If an overlap
   check is available, save at most a scalar identity-error diagnostic.

Trust boundary:

- Do not require H1 to materialize.
- Do not require H1-J/density interaction to materialize.
- Do not run or wire private RHF for H2.
- Do not compare to supplemented WL/QW HF/ED references.
- Do not add supplement support.
- Do not revive Be2/Cr2 artifact work.
- Do not add HFDMRG, DMRG, ECP, exports, public solver behavior, or Qiu-White
  correction work.
- Do not add this test to default runners.
- Do not request interactive approval or sandbox escalation. If approval is
  genuinely required, write `.agent_handoffs/ATTENTION.md` and stop.

Line-budget rule:

This pass must be net-negative for:

```text
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Because this pass adds a driver input and a test, pay for it in the same pass by
deleting or shrinking stale development scaffolding. Do not delete scientific
endpoint/reference tests, the explicit He driver endpoints, or live runner
coverage. If you cannot keep the scoped budget net-negative without unsafe
deletion, write `.agent_handoffs/ATTENTION.md` and stop.

Useful deletion-search starting point, not authority:

- Look for remaining mixed one-body/pair-block metadata scaffold tests that are
  not included in default/integration runners and have no source callers.
- Verify every deletion with `rg`; do not delete by pattern alone.

Validation:

Run only focused checks:

```text
julia --project=. test/nested/<new_h2_driver_readiness_test>.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
git status --short --branch
```

If the driver endpoint takes more than 60 seconds, report the elapsed time and
the phase if available, but do not broaden the suite.

Report back:

- artifact path produced by the test;
- H2 geometry/config/comparison fields verified;
- route/materialization readiness status and exact blocker, if any;
- scoped line-budget arithmetic;
- stale code/test deleted or shrunk;
- validation results;
- deletion/shrinkage report with exact remaining caller/blocker.

-- repo-manager@macmini
