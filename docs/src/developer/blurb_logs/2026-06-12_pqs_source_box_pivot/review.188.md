Pass 188 review: accepted.

This pass successfully moved the accepted He q=5/n_s=5 PQS H1/H1-J gate from a
hand-built direct test into the visible driver path:

```text
test/driver_inputs/he_pqs_q5_wlmap.jl
-> bin/cartesian_ham_builder.jl
-> saved JLD2 artifact
-> thin readback test
```

Key checks:

- The visible driver spine remains intact. There is no opaque
  `run_driver(config)` wrapper.
- The one-center parent builder now accepts
  `parent_mapping_rule = :white_lindsey_atomic_mapping` and preserves
  `:identity_mapping` as the default.
- The driver input uses the intended WL-mapped He setup:
  `Z = 2`, `d = 0.3`, `tail_spacing = 10.0`, `q = n_s = 5`,
  parent axis counts `(11, 11, 11)`.
- The artifact writes the requested compact groups:
  `config`, `basis`, `physics`, `density_interaction`, and `comparison`.
- The new test checks the live physics facts:
  final dimension `419`, three shell layers, retained shell counts
  `(98, 98, 98)`, final-overlap identity noise, H1 lowest, H1-J self-Coulomb,
  density gauge, and WL deltas.
- The direct hand-built test
  `test/nested/pqs_direct_retained_final_h1_runtests.jl` was deleted, and its
  default-runner include was replaced.

Validation I reran:

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
julia --project=. test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl
git diff --check
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Results:

```text
load ok
cartesian_ham_builder He PQS q5 WL-map artifact: 10/10 pass, 3m04.7s
git diff --check: passed
line budget: 176 added, 178 deleted, net -2
```

Runtime note:

The focused driver test is now a multi-minute endpoint when cold. The output
shows the dominant phases are first-call compilation and assembly:
transform around 12s, pair terms around 25s, assembly around 28s, plus driver
report/save output. This is acceptable for this replacement pass, but before
adding private RHF or more driver endpoint tests, the next pass should audit the
driver test lane and runtime policy. We should not let the default nested suite
quietly grow into a cold multi-minute physics endpoint lane without an explicit
decision.

Deletion/shrinkage:

- deleted: `test/nested/pqs_direct_retained_final_h1_runtests.jl`
- simplified: the He PQS gate is now driver-owned and artifact-readback based.
- quarantined: none.
- not deleted because: the new compact artifact writer and driver input are
  live surfaces for the driver-owned physics gate.
- exact remaining caller/blocker: `test/nested/runtests.jl` now includes the
  driver-owned He PQS artifact test; its default-runner placement should be
  audited before adding private RHF.

Recommended next pass:

No-edit audit of the private RHF driver seam and the driver endpoint test lane:
decide whether the He driver endpoint remains in the default nested runner,
what artifact keys are needed for optional private RHF, and how to avoid
adding another multi-minute default test.

-- repo-manager@macmini
