Pass 189 - keep the He driver endpoint explicit, not a default nested gate, and audit private RHF seam.

Purpose:

Pass 188 did the important architectural replacement:

```text
hand-built He PQS test
-> visible driver input
-> bin/cartesian_ham_builder.jl
-> saved artifact
-> thin readback test
```

That was the right direction. But the new driver endpoint took several minutes
when run cold, mostly due to compilation and assembly. It should remain as a
real driver-owned physics endpoint test, but it should not quietly become a
default per-edit nested runner cost.

This pass has two parts:

1. Make the He driver endpoint explicit/manual rather than default nested.
2. Audit the optional private RHF driver seam, but do not implement RHF yet.

Task type:

Small implementation cleanup plus read-only audit. This pass must still be
source/test/bin/generator net-negative.

Exact implementation task:

Remove the include for:

```text
test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl
```

from:

```text
test/nested/runtests.jl
```

Do not delete the driver endpoint test itself. It is now the explicit physics
endpoint to run directly:

```text
julia --project=. test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl
```

Reason:

- The endpoint is valuable and should be preserved.
- It is not a compact module-contract test.
- Scientific endpoint tests are high value as explicit acceptance gates, but
  should not automatically become the default per-edit nested lane when they
  cost multiple minutes cold.

Do not add a replacement default-runner test in this pass. The old hand-built
test is already deleted. Do not reintroduce it or add a metadata-only substitute.

Exact audit task:

Read the current private RHF surfaces and the new driver artifact path enough
to answer what a future optional private RHF driver flag would require.

Inspect at least:

```text
bin/cartesian_ham_builder.jl
src/pqs_source_box_route_driver_helpers.jl
src/pqs_source_box_route_driver_reporting.jl
src/pqs_multilayer_complete_core_shell_rhf.jl
test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl
```

Answer:

1. What current private RHF helper or payload should the driver call, if any?
2. What existing data is already available from the complete-core/shell
   diagnostic route payload or artifact?
3. What small input flag names would be appropriate?

   Candidate names:

   ```julia
   run_private_rhf = false
   private_rhf_max_iter = ...
   private_rhf_residual_atol = ...
   private_rhf_fixture_role = :route_smoke
   ```

4. What artifact keys would be needed for optional RHF?

   Keep this compact. Likely keys:

   ```text
   private_rhf/status
   private_rhf/blocker
   private_rhf/total_energy
   private_rhf/iteration_count
   private_rhf/converged
   private_rhf/residual
   comparison/wl_rhf_total
   comparison/delta_rhf
   ```

5. What test shape should cover RHF later?

   It should run the same driver input with `run_private_rhf=true`, read the
   artifact, and check the known route-smoke values. It should not reconstruct
   RHF internals in the test.

Known RHF route-smoke comparison values from the WL-mapped He 419 diagnostics:

```text
PQS RHF total ≈ -2.850817886618113
WL  RHF total ≈ -2.85080350301779
delta ≈ -1.4383600322798173e-5
```

Trust boundary:

- Do not implement `run_private_rhf` yet.
- Do not add RHF artifact keys yet.
- Do not add a new RHF test yet.
- Do not run broad nested/integration suites.
- Do not delete the new He driver endpoint test.
- Do not restore the deleted hand-built direct He test.
- Do not add H2, Be2, Cr2, HFDMRG, DMRG, exports, public solver behavior, or
  public APIs.
- Do not request interactive approval or sandbox escalation. If approval is
  genuinely required, write `.agent_handoffs/ATTENTION.md` and stop.

Line-budget rule:

This pass edits `test/nested/runtests.jl`, so it must be net-negative for:

```text
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Expected source/test/bin/generator line budget:

```text
added: 0
deleted: at least 1
net: negative
```

Do not add source/test/bin lines in this pass unless you delete more stale
source/test/bin lines in the same pass.

Validation:

```text
rg -n "cartesian_ham_builder_he_pqs_q5_wlmap" test/nested/runtests.jl test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
git status --short --branch
```

Do not rerun the multi-minute driver endpoint unless you changed the endpoint
test itself, which this pass should not do.

Report:

- confirm the default nested runner no longer includes the multi-minute driver
  endpoint;
- confirm the explicit endpoint test remains available;
- summarize the private RHF driver seam audit;
- recommended next implementation pass for RHF, or blocker if not ready;
- line-budget arithmetic;
- validation results;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Write the result to `.agent_handoffs/response.189.md` and copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.189.md
```

-- repo-manager@macmini
