Pass 070 response: PQS support one-body assembly seam audit

Summary:
- Added a docs-only audit note to `docs/src/developer/pqs_near_term_final_basis_realization_plan.md` under the driver/assembly section.
- The note records that support kinetic assembly is safe to promote as a narrow helper over the multi-layer plan support states and `plan.metrics`.
- The note records that support electron-nuclear assembly should wait for an explicit by-center convention before promotion.

Audit answers:
- Support kinetic:
  - Safe to promote now.
  - Proposed helper shape: `pqs_multilayer_support_kinetic_matrix(plan)` or an equivalent narrow support one-body helper.
  - Inputs: route-owned multi-layer PQS plan, using plan support states and axis overlap/kinetic metrics.
  - Output: support-space kinetic matrix plus compact provenance.
  - Nonclaims: no final-basis transfer, H1, IDA, RHF, driver wiring, exports, or artifacts.
- Support electron-nuclear:
  - Not safe to promote as a generic helper until the by-center convention is explicit.
  - The promoted helper should produce separated by-center support matrices, keep them uncharged, and leave charge application and center summation to Hamiltonian assembly.
  - Sign convention should be explicit: the uncharged by-center matrix represents the electron-nuclear potential `-1/r_center`; Hamiltonian assembly can then add `Z_center * V_center` if `V_center` already carries the negative sign.
  - Current H1 gate is centered/origin only through `pgdg_intermediate.gaussian_factor_terms`.
  - Off-origin centers need the centered Gaussian factor-term source already used by retained by-center PQS kernels before route-owned support assembly is promoted.
  - Old fixed-block and WL matrices remain oracle comparisons only.
- If implemented next, the helper should shrink the remaining test-local `_pqs_h1_support_kinetic_matrix` and `_pqs_h1_support_nuclear_matrix` code.

Validation:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - PASS: `load ok`
- `git diff --check`
  - PASS

Diff:
- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
  - 30 insertions

Deletion/shrinkage report:
- No source, test, metadata, or compatibility path became unnecessary in this docs-only audit pass.
- Nothing was deleted because the blurb explicitly asked not to implement the helper or change tests; the remaining local H1 support helpers still serve the current gate.
- No new tests were added. This was not new live-contract coverage; it was an audit note to prevent promoting unclear nuclear support semantics.
- Remaining stale/duplicate surfaces to retire next:
  - `_pqs_h1_support_kinetic_matrix` after a route-owned support kinetic helper exists;
  - `_pqs_h1_support_nuclear_matrix` after the uncharged separated by-center support nuclear helper has a reviewed centered/off-origin convention.

-- repo-doer@macmini
