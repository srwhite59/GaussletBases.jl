Pass 190 review: accepted.

This pass added optional private RHF as a driver-owned, default-off diagnostic
for the He 419 PQS endpoint. It did not make RHF public, default, exported, or
downstream-solver-facing.

Key accepted changes:

- Added `run_private_rhf = false` and compact private RHF controls to
  `bin/cartesian_ham_builder.jl`.
- Routed the controls through the existing visible driver spine.
- Added a route-owned private RHF payload built from
  `complete_core_shell_diagnostic_route_payload` using the existing private RHF
  input contract and SCF payload.
- Wrote only compact RHF artifact scalars under `private_rhf/*` and
  `comparison/delta_rhf`; no matrices, densities, orbitals, Fock matrices, or
  iteration histories were written.
- Added an explicit RHF endpoint test, not included by default runners.
- Deleted stale probe pressure:
  `test/nested/pqs_route_driver_one_center_materializer_probe_runtests.jl` and
  `test/nested/cartesian_white_lindsey_adapter_timing_probe.jl`.

Validation I reran:

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
julia --project=. test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_rhf_runtests.jl
git diff --check
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Results:

```text
load ok
cartesian_ham_builder_he_pqs_q5_wlmap_rhf artifact: 5/5 pass, 3m15.2s
git diff --check: passed
tracked numstat excluding the untracked test: 187 added, 335 deleted
new RHF endpoint test: 40 lines
true source/test/bin/generator budget: 227 added, 335 deleted, net -108
```

The endpoint reproduces the expected route-smoke values:

```text
PQS RHF total ≈ -2.850817886618113
WL  RHF total ≈ -2.85080350301779
delta ≈ -1.4383600322798173e-5
```

Runtime note:

The explicit RHF endpoint is still multi-minute cold and compilation-heavy.
Keeping it out of default nested/integration runners is correct. It is a
physics endpoint, not a per-edit module-contract gate.

Deletion/shrinkage:

- deleted:
  `test/nested/pqs_route_driver_one_center_materializer_probe_runtests.jl`;
  `test/nested/cartesian_white_lindsey_adapter_timing_probe.jl`
- simplified: private RHF now runs through the driver artifact path instead of a
  hand-built test fixture.
- quarantined: RHF remains default-off, private, and explicit-endpoint only.
- not deleted because: the explicit H1/H1-J and RHF driver endpoint tests are
  current He PQS physics gates.
- exact remaining caller/blocker: no default runner includes the RHF endpoint;
  H2/PQS should not start until the next pass defines the H2 driver-owned target
  and WL reference path clearly.

Recommended next pass:

Audit the H2 driver-owned PQS target against the old WL/QW H2 reference before
implementation. The audit should identify the exact driver input shape, WL
baseline values, shellification/parent mapping policy, artifact keys, and stale
test/probe deletion target that will pay for implementation.

-- repo-manager@macmini
