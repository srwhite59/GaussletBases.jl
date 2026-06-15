Pass 256 - independent H2 PQS supplement staging audit

Context:
- Current HEAD should include `e5167d7a Shrink atom-growth report-stage RouteCore mirrors`.
- Independent H2 PQS now has fake-free diagnostics through private RHF:
  readiness -> final basis -> H1 -> H1-J -> private RHF.
- The route remains private/diagnostic-only. Supplements, CR2/export, public
  API, and public solver readiness remain blocked.
- MT4 says MWG/GTO supplements may only be staged after retained-transform
  authority is clear. That condition is now much better, but provider-block
  work should not start until the exact seam is reviewed.

Relevant supplement surfaces:
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
  - `_pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_request_payload`
  - `_pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_representation_payload`
  - `_pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_preflight_payload`
- `src/pqs_source_box_route_driver_reporting.jl`
  - artifact groups `supplement_request`, `supplement_representation`,
    `supplement_preflight`
- Existing provider/kernel tests that may be relevant as references only:
  - `test/nested/cartesian_cpb_mixed_gto_overlap_block_runtests.jl`
  - `test/nested/cartesian_cpb_gto_supplement_one_body_block_runtests.jl`
  - `test/nested/cartesian_cpb_gto_nuclear_by_center_block_runtests.jl`
  - `test/nested/cartesian_combined_gto_density_density_readiness_runtests.jl`
- Existing fake-PQS supplement preflight test is schema/history only:
  - `test/nested/cartesian_ham_builder_h2_fake_pqs_wl_source_backed_r4_runtests.jl`

Task:
Do a focused audit of the next independent-H2-PQS MWG/GTO supplement staging
seam. This is primarily a no-physics audit pass.

Questions to answer:
1. Can the existing supplement request/representation/preflight payloads already
   attach to the independent H2 PQS route kind and artifact roles, or are they
   still effectively fake/old-physical-route shaped?
2. What should the first independent-H2 supplement input be called if created
   later? For example, should it include from private RHF or from H1-J and set
   `supplement_policy = :mwg_residual_gto` with all physics flags off?
3. What exact artifact facts should the first independent supplement preflight
   expose?
   - `fake_pqs/enabled = false`
   - `source_backed_fixed_source_oracle_used = false`
   - `retained_transform_authority = :pqs_source_box_construction`
   - support counts `(275, 578, 362)`
   - retained counts `(275, 98, 98)`
   - gausslet final dimension `471`
   - supplement basis `H/cc-pVTZ`, `lmax = 1`, orbital count `18`
   - blocker `:missing_provider_gto_supplement_blocks` unless a narrower
     blocker is more honest
4. Which existing provider/kernel objects can be reused, and what is still
   missing for independent PQS?
   - mixed gausslet/GTO overlap and one-body blocks;
   - GTO/GTO self blocks;
   - mixed nuclear by-center blocks;
   - combined raw moment matrices;
   - residual MWG representation;
   - combined density-density readiness.
5. What must not be copied from the fake-PQS/WL source-backed supplement
   preflight?
6. What is the smallest safe pass after this audit: input/preflight only,
   provider-block object, or another cleanup first?

Allowed edits:
- Prefer no source edits.
- If the audit can be reported entirely in `response.256.md`, do that.
- To keep cleanup pressure, you may make a small tests-only deletion offset in
  already-classified low-order/report-stage stale mirror assertions. Do not
  touch supplement tests unless the audit proves the assertions are stale and
  not protecting a current contract.

Strict exclusions:
- Do not implement provider blocks.
- Do not build mixed gausslet/GTO matrices.
- Do not build residual MWG representation.
- Do not add supplemented H1/H1-J/RHF values.
- Do not compare to supplemented WL/QW scalar references.
- Do not mutate fake-PQS as independent-PQS evidence.
- Do not add public/export/CR2/HamV6 readiness.

Validation:
- If no source is touched, `git diff --check` is enough plus any parse smoke for
  a test file you cleaned.
- If source is touched, run package load.
- Do not run slow H2 route tests or broad stale integration gates.

Line budget:
- If source/test/bin are touched, keep scoped `src + test + bin` net-negative.
- If this remains truly no-edit, report scoped line impact as `0`.

Report:
- Direct answers to the audit questions above.
- Recommended next pass with exact allowed/forbidden surfaces.
- Any cleanup offset performed.
- Validation commands.
- Scoped line count for `src + test + bin`.
- Deletion/shrinkage result:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
