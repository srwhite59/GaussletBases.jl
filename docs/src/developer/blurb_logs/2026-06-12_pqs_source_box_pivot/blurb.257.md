Pass 257 - independent H2 PQS supplement preflight input

Context:
- Current HEAD should include `579b4684 Record independent H2 PQS supplement audit`.
- Pass 256 audited the MWG/GTO supplement staging seam. The conclusion was:
  the request/representation/preflight helpers can mechanically attach to the
  independent target, but the first safe pass should be input/preflight only.
- Provider blocks, mixed matrices, residual MWG representation, combined
  density readiness, supplemented values, CR2/export, and public API remain
  blocked.

Task:
Add the explicit independent-H2-PQS supplement-preflight input/artifact role.

Required changes:
1. Add a tiny include/override input:
   `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_supplement_preflight.jl`
   based on `h2_pqs_q5_independent_source_box_r4_h1_j.jl`.
2. The input should set:
   - `artifact_role = :independent_h2_pqs_supplement_preflight_diagnostic`
   - `supplement_policy = :mwg_residual_gto`
   - `comparison_ready = false`
   - `comparison_blocker = :independent_pqs_supplement_preflight_only`
   - `physics_endpoint_ready = false`
   - `physics_endpoint_blocker = :missing_provider_gto_supplement_blocks`
   - `run_final_basis = true`
   - `run_h1 = true`
   - `run_h1_j = true`
   - `run_private_rhf = false`
3. Update the compact independent-H2 artifact-role classifier so this new role
   still receives:
   - `source_backed_fixed_source_oracle_used = false`
   - `retained_transform_authority = :pqs_source_box_construction`
   - fake-PQS guard behavior unchanged.
4. Add a manifest row stating this is independent H2 PQS supplement preflight,
   not provider-block/supplemented-value/public/export ready.

Do not:
- Implement provider blocks.
- Build mixed gausslet/GTO matrices.
- Build GTO/GTO matrices.
- Build residual MWG representation.
- Build combined density-density readiness.
- Add supplemented H1/H1-J/RHF values.
- Compare to supplemented WL/QW scalar references.
- Use fake-PQS evidence as independent-PQS evidence.
- Add public/export/CR2/HamV6 readiness.

Validation:
- Include/flag smoke for the new input. Assert at least:
  - independent H2 route family/kind;
  - artifact role;
  - `supplement_policy === :mwg_residual_gto`;
  - `comparison_ready == false`;
  - blocker is the supplement-preflight-only blocker;
  - `(run_final_basis, run_h1, run_h1_j, run_private_rhf) ==
    (true, true, true, false)`.
- Run package load if source is touched.
- Run `git diff --check`.
- Do not run the slow H2 route unless the classifier/report change cannot be
  validated by smoke/load.

Line budget:
- Keep scoped `src + test + bin` net-negative.
- Offset new input/classifier lines with a small stale low-order/report-stage
  mirror deletion if needed.

Report:
- New input path and flag values.
- Classifier/manifest updates.
- Confirmation that provider blocks and supplemented values remain blocked.
- Validation commands.
- Scoped line count for `src + test + bin`.
- Deletion/shrinkage result:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
