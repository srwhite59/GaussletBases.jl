Pass 189 review: accepted.

This pass did the right cleanup after pass 188:

- removed the multi-minute He driver endpoint from `test/nested/runtests.jl`;
- kept `test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl` as an
  explicit physics endpoint test;
- changed only the testset label so grep validation can distinguish the
  explicit endpoint from the default nested runner;
- audited the optional private RHF driver seam without implementing RHF.

Validation I reran:

```text
rg -n "cartesian_ham_builder_he_pqs_q5_wlmap" test/nested/runtests.jl test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Results:

```text
rg: endpoint appears only in the explicit endpoint test
load ok
git diff --check: passed
line budget: 1 added, 2 deleted, net -1
```

Accepted RHF audit conclusions:

- A future driver flag should call the existing private complete-core/shell RHF
  contract and SCF payload, not expose the initial-density or one-step helpers
  as separate driver stages.
- The needed route data are already available through the
  `complete_core_shell_diagnostic_route_payload`: source plan, Coulomb expansion,
  final basis, H1 payload, density inputs, H1-J payload, and density interaction.
- Suggested driver controls are reasonable:
  `run_private_rhf`, `private_rhf_max_iterations`, density/energy/residual
  tolerances, `private_rhf_mixing_kind`, `private_rhf_fixture_role`, and
  `wl_rhf_total`.
- Future artifact keys should be compact under `private_rhf/*` plus
  `comparison/wl_rhf_total` and `comparison/delta_rhf`.
- The future RHF test should be explicit, not default nested, and should read
  the driver artifact rather than reconstructing SCF internals.

Deletion/shrinkage:

- deleted: default nested include for the multi-minute endpoint.
- simplified: default nested runner no longer pays the explicit physics
  endpoint cost.
- quarantined: He driver endpoint remains explicit/manual.
- not deleted because: the endpoint is the current driver-owned He PQS physics
  acceptance gate.
- exact remaining caller/blocker: optional private RHF driver behavior is not
  implemented; next pass can add it behind `run_private_rhf=false` defaults,
  but must keep endpoint testing explicit and source/test/bin line budget
  negative.

Recommended next pass:

Implement optional private RHF in the driver artifact path, using the existing
private RHF helpers, with the endpoint test kept explicit and no default-runner
include.

-- repo-manager@macmini
