Pass 259 - independent H2 PQS supplement provider-block seam audit

Context:
- Current HEAD should include
  `1f1e31e2 Record independent H2 PQS supplement preflight`.
- Pass 258 verified the independent supplement-preflight artifact:
  support counts `(275, 578, 362)`, retained counts `(275, 98, 98)`, final
  dimension `471`, H/cc-pVTZ lmax-1 representation with 18 orbitals, and
  blocker `:missing_provider_gto_supplement_blocks`.
- The next implementation must not jump directly to supplemented values. First
  we need the provider-block seam clearly specified.
- After this pass is accepted, manager will add the required medium-term
  checkpoint for passes 255-259.

Task:
Do a no-implementation audit of the first route-owned provider-block payload for
independent H2 PQS MWG/GTO supplement staging.

Questions to answer:
1. What exact route-owned payload should be implemented first?
   - Suggested name if appropriate:
     `_pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_provider_blocks_payload`
     or a better independent-PQS name.
2. What are the required inputs?
   - independent target/source-plan/final-basis payloads;
   - supplement request/representation payloads;
   - parent axis/factor data;
   - source/final coefficient maps;
   - CPB/local row ownership or equivalent;
   - nuclear centers/charges for by-center blocks.
3. Which existing provider functions should it call or wrap?
   - mixed gausslet/GTO overlap/position/x2/kinetic;
   - GTO/GTO supplement self overlap/position/x2/kinetic;
   - mixed and self nuclear by-center blocks.
4. What should the payload produce in the first implementation pass?
   - Prefer a compact readiness/materialization summary and small set of block
     handles/fingerprints over broad scalar report fields.
   - State whether matrices should remain local/provider-level, route-global, or
     absent in the first pass.
5. What exact artifact fields/status changes should result after provider
   blocks exist, and what blockers should remain?
6. What must remain explicitly blocked after provider blocks?
   - combined raw moment matrices;
   - residual MWG representation;
   - combined density-density readiness;
   - supplemented H1/H1-J/RHF values;
   - CR2/export/HamV6/public API readiness.
7. What is the smallest validation fixture for the first provider-block
   implementation?
   - Avoid the slow full H2 preflight route unless needed.
   - Prefer a compact route-owned smoke over copying large CPB provider tests.

Strict exclusions:
- Do not implement provider blocks in this pass.
- Do not build mixed matrices, residual MWG representation, combined density
  readiness, or supplemented values.
- Do not edit source unless you find and delete an unrelated stale test mirror
  as cleanup offset.
- Do not use fake-PQS/WL source-backed data as independent-PQS evidence.
- Do not compare to supplemented WL/QW scalar references.

Validation:
- No Julia run is required for a no-edit audit.
- Run `git diff --check`.
- If you do a small cleanup offset, run parse smoke for the touched test file.

Line budget:
- Prefer no source/test/bin edits; scoped impact `0` is acceptable for this
  design audit.
- If source/test/bin are touched, keep scoped net-negative.

Report:
- Direct answers to the audit questions.
- Recommended exact next implementation pass.
- Any cleanup offset performed.
- Validation command(s).
- Scoped line count for `src + test + bin`.
- Deletion/shrinkage result:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
