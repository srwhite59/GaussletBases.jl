Pass 258 response - verify independent H2 PQS supplement preflight artifact

Result:
- The new independent H2 PQS supplement-preflight input produced the expected artifact facts.
- No source, test, bin, docs, or driver-input fix was needed.
- Provider blocks and supplemented values remain blocked.

Focused route command:
```sh
julia --project=. tmp/work/h2_independent_pqs_supplement_preflight_probe.jl
```

Temporary local probe:
- `tmp/work/h2_independent_pqs_supplement_preflight_probe.jl`
- Writes artifact: `/private/tmp/h2_independent_pqs_supplement_preflight_258.jld2`
- Captures driver stdout: `/private/tmp/h2_independent_pqs_supplement_preflight_258_stdout.txt`
- First attempt failed before route execution because the local probe included `bin/cartesian_ham_builder.jl` relative to `tmp/work/`; fixed only the local probe path and reran.

Elapsed time:
- Focused route/artifact probe elapsed: `78.408378416` seconds.

Artifact facts observed:
```text
artifact_role=independent_h2_pqs_supplement_preflight_diagnostic
fake_pqs_enabled=false
source_backed_fixed_source_oracle_used=false
retained_transform_authority=pqs_source_box_construction
supplement_policy=mwg_residual_gto
private_rhf_requested=false
supplement_request_status=available_pqs_physical_gausslet_supplement_request
supplement_representation_status=available_pqs_physical_gausslet_gto_supplement_representation
supplement_representation_orbital_count=18
supplement_preflight_status=blocked_pqs_physical_gausslet_mwg_residual_gto_preflight
supplement_preflight_blocker=missing_provider_gto_supplement_blocks
supplement_preflight_support_counts=(275, 578, 362)
supplement_preflight_retained_counts=(275, 98, 98)
supplement_preflight_gausslet_final_dimension=471
supplement_preflight_missing_fact_labels=(:missing_provider_gto_supplement_blocks, :missing_mixed_gausslet_gto_blocks, :missing_gto_gto_blocks, :missing_combined_raw_moment_matrices, :missing_residual_mwg_representation, :missing_combined_density_density_readiness)
```

Additional asserted facts:
- `config/comparison_ready == false`
- `physics/endpoint_ready == false`
- missing fact labels do not include fake-PQS/source-backed adapter facts
- missing fact labels do not include WL/QW scalar-reference blockers

Validation:
- `julia --project=. tmp/work/h2_independent_pqs_supplement_preflight_probe.jl`
  - passed on rerun; artifact facts above
- `git diff --check`
  - passed
- No package load was run because no source was touched in this pass.
- No broad stale integration gates were run.

Scoped line count:
- `src + test + bin`: 0 tracked change.
- The local probe is ignored under `tmp/work/` and not part of tracked scope.

Git status:
```text
## main...origin/main
```

Deletion/shrinkage result:
- deleted: none.
- simplified: none.
- quarantined: provider blocks, residual MWG representation, combined density-density readiness, supplemented values, CR2/export/HamV6/public API remain blocked.
- not deleted because: artifact was already correct, so no cleanup offset was needed.
- exact remaining caller/blocker: `cartesian_assembly` still builds supplement request/representation/preflight metadata; route advancement remains blocked at `:missing_provider_gto_supplement_blocks`.

-- repo-doer@macmini
