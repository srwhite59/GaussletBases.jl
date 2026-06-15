Pass 258 - verify independent H2 PQS supplement preflight artifact

Context:
- Current HEAD should include
  `51ad985f Add independent H2 PQS supplement preflight input`.
- Pass 257 added:
  `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_supplement_preflight.jl`.
- The new input is preflight-only for MWG/GTO supplements. It must not imply
  provider-block availability, supplemented values, CR2/export readiness, or
  public API readiness.

Task:
Run the new independent H2 PQS supplement-preflight input once through the
driver, inspect the JLD2 artifact, and repair only narrow preflight fact plumbing
if the artifact is wrong.

Allowed route run:
- This focused driver run is allowed to exceed 60 seconds because it builds the
  current independent H2 final basis/H1/H1-J context before writing supplement
  preflight facts.

Expected artifact facts:
- `route/artifact_role == :independent_h2_pqs_supplement_preflight_diagnostic`.
- `route/fake_pqs_enabled == false`.
- `route/source_backed_fixed_source_oracle_used == false`.
- `route/retained_transform_authority == :pqs_source_box_construction`.
- `config/supplement_policy == :mwg_residual_gto`.
- `config/comparison_ready == false`.
- `physics/endpoint_ready == false`.
- `private_rhf/requested == false`.
- `supplement_request/status ==
  :available_pqs_physical_gausslet_supplement_request`.
- `supplement_representation/status ==
  :available_pqs_physical_gausslet_gto_supplement_representation`.
- `supplement_representation/orbital_count == 18`.
- `supplement_preflight/status ==
  :blocked_pqs_physical_gausslet_mwg_residual_gto_preflight`.
- `supplement_preflight/blocker == :missing_provider_gto_supplement_blocks`.
- `supplement_preflight/support_counts == (275, 578, 362)`.
- `supplement_preflight/retained_counts == (275, 98, 98)`.
- `supplement_preflight/gausslet_final_dimension == 471`.
- Missing fact labels include provider/mixed/GTO/MWG/density blockers, not fake
  PQS or WL/QW scalar-reference blockers.

If the artifact is already correct:
- Prefer no source edits.
- A small tests-only cleanup offset is allowed if it is clearly stale, but do
  not delay the pass hunting for one.

If the artifact is wrong:
- Repair only the narrow preflight/reporting/classifier fact plumbing.
- Keep provider blocks and supplemented values blocked.
- Keep the fix scoped and line-negative if source/test/bin are touched.

Strict exclusions:
- Do not implement provider blocks.
- Do not build mixed gausslet/GTO matrices.
- Do not build residual MWG representation.
- Do not add supplemented H1/H1-J/RHF values.
- Do not compare to supplemented WL/QW scalar references.
- Do not use fake-PQS evidence as independent-PQS evidence.
- Do not add CR2/export/HamV6/public API readiness.

Validation:
- Run the focused driver/artifact probe for the new input.
- Run package load if source is touched.
- Run `git diff --check`.
- Do not run broad stale integration gates.

Report:
- Artifact facts observed, especially supplement request/representation/
  preflight status, counts, dimensions, missing facts, and blockers.
- Whether any source fix was needed.
- Validation command(s) and elapsed time for the focused route run.
- Scoped line count for `src + test + bin`.
- Deletion/shrinkage result:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
