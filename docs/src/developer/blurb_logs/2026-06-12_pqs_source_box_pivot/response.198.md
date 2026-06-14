Pass 198 response - H2 source-plan/final-basis inventory audit before H1-J.

Task scope:

- No H2 implementation changes.
- No H1-J, density interaction, private RHF, supplemented comparison, export, or public solver work.
- Deleted one safe stale scaffold test.

Inventory findings from live route/source-plan objects:

- Probe command:
  - `julia --project=. -e '<focused staged-driver audit probe with run_h1=false, run_h1_j=false, run_private_rhf=false, save_artifact=false, save_tsv=false>'`
  - `AUDIT|elapsed_s=84.579863833`
- Parent axis counts:
  - `(x = 9, y = 9, z = 15)`
  - parent product size is `9 * 9 * 15 = 1215`
- Source plan status:
  - `:available_pqs_diatomic_complete_core_shell_source_plan`
- Source boxes from the live driver report:
  - `pqs_left = (x = 1:5, y = 1:5, z = 1:5)`
  - `product = (x = 1:5, y = 1:5, z = 8:8)`
  - `pqs_right = (x = 1:5, y = 1:5, z = 11:15)`
- Retained unit keys/order:
  - `(:pqs_left, :pqs_right, :product)`
- Retained unit kinds:
  - `(:pqs, :pqs, :product_doside)`
- Support order:
  - `(:product, :pqs_left, :pqs_right)`
- Final support row order:
  - `:core_then_shell`
- Source dimensions:
  - `pqs_left = (5, 5, 5)`
  - `pqs_right = (5, 5, 5)`
  - `product = (5, 5, 1)`
- Retained counts:
  - `pqs_left = 98`
  - `pqs_right = 98`
  - `product = 25`
- Source-plan core/shell split:
  - `core_unit_key = :product`
  - `shell_unit_keys = (:pqs_left, :pqs_right)`
  - `core_support_count = 25`
  - `shell_support_count = 250`
  - `shell_support_counts = (pqs_left = 125, pqs_right = 125)`
  - `shell_retained_count = 196`
  - `shell_retained_counts = (pqs_left = 98, pqs_right = 98)`
- Shell coefficients:
  - `shell_final_coefficients_shape = (250, 196)`
  - `shell_coefficient_block_structure = :block_diagonal_left_right_pqs`
- Final basis:
  - `final_basis_status = :available_pqs_complete_core_shell_final_basis`
  - `final_dimension = 221`
  - `final_overlap_identity_error = 2.510171831978689e-13`

Why the dimension is 221:

- Each atom-side PQS raw source box is `5^3 = 125` source modes.
- The retained atom-side PQS count is boundary-only:
  - `5^3 - 3^3 = 125 - 27 = 98`
  - left atom box: `98`
  - right atom box: `98`
- The midpoint/product slab is `5 * 5 * 1 = 25` and is fully retained as `:product_doside`.
- Total retained dimension:
  - `98 + 98 + 25 = 221`
- The source plan also reports `precleanup_retained_dimension = 221`, and the final basis keeps `final_dimension = 221`.

Whether full atom-core interiors are present:

- Full `5^3` atom-side raw source boxes are present as source/support objects.
- Full `5^3` atom-core interiors are not present in the retained/final H2 basis.
- The atom-side PQS units retain only the 98 boundary COM/product modes per atom box.
- The 27 interior modes per `5^3` atom box are not retained as atom-core bound-state content.
- The only fully retained interior-like block is the midpoint/product slab with 25 functions.
- There is no multishell atom-centered physics basis here; there are two atom-side PQS boundary units plus one midpoint/product slab.

Comparison to intended H2 physics target:

- The current parent row count `1215` matches the old WL/QW complete-rectangular parent row size.
- The current final dimension `221` does not match the old WL/QW default final dimension `481`.
- The old documented WL/QW H2 physics reference also carries an H/cc-pVTZ S/P supplement, so its HF/ED totals are not directly comparable to this gausslet-only diagnostic.
- More importantly, this 221-dimensional basis is not the user-described physical minimal atom-centered target because it does not retain full `5^3` atom cores before adding shells.

Interpretation of the positive H1:

- The positive lowest H1 value from pass 197 is consistent with the current basis being a source-box/boundary diagnostic rather than a physical atom-centered H2 basis.
- With atom interiors removed from the retained atom-side units, the basis lacks the near-nuclear bound-state content expected for a one-electron H2/H2+ style check.
- The current H1 result should therefore be treated as a route-shape warning, not as evidence that H1-J/RHF should proceed on this target.

Recommended next target:

- B. Replace/extend the H2 driver target to a physical gausslet-only basis with full atom-core interiors plus shell layers before continuing to H1-J/RHF.
- Exact next seam:
  - introduce or select a diatomic source-plan/final-basis route variant where the atom-side full `5^3` interiors are retained as atom-core units;
  - layer PQS/shell retained modes around those atom cores;
  - keep the midpoint/product or bridge region explicit if it is still part of the intended H2 route;
  - label the current 221-dimensional route as route-smoke/source-box diagnostic only before any future H1-J/RHF work.

Stale scaffold deleted:

- `test/nested/cartesian_pair_block_one_body_factor_inputs_runtests.jl`
- Why safe:
  - not included by default/integration/top-level runners;
  - no remaining basename/testset references in `test src bin`;
  - not an explicit He/H2 driver endpoint;
  - not a WL H/H2 reference or scientific endpoint;
  - not needed by this H2 source-plan/final-basis audit;
  - it was a private mixed one-body input-convention scaffold, while the live routine path remains the smaller mixed one-body consumer smoke/materialization contract coverage.

Validation results:

- `rg -n "cartesian_pair_block_one_body_factor_inputs_runtests|CartesianPairBlockMaterialization one-body named-tuple factor inputs|CartesianPairBlockMaterialization one-body provider factor inputs|CartesianPairBlockMaterialization one-body term-set factor inputs" test src bin`
  - no output; exit 1 from no matches
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check`
  - passed with no output
- `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - `0	247	test/nested/cartesian_pair_block_one_body_factor_inputs_runtests.jl`
- `git status --short --branch`
  - `## main...origin/main`
  - ` D test/nested/cartesian_pair_block_one_body_factor_inputs_runtests.jl`

Deletion/shrinkage report:

- deleted:
  - `test/nested/cartesian_pair_block_one_body_factor_inputs_runtests.jl`
- simplified:
  - removed stale private mixed one-body factor-input scaffold pressure from nested tests.
- quarantined:
  - none
- not deleted because:
  - the active H2 readiness driver endpoint remains the current audit target;
  - H2 source-plan/final-basis code was intentionally not changed in this audit pass;
  - scientific endpoints, WL/H2 references, and explicit He/H2 driver endpoints were preserved.
- exact remaining caller/blocker:
  - blocker before H1-J/RHF is the route target mismatch: the current 221-dimensional H2 basis retains atom-side boundary modes only and lacks full `5^3` atom-core retained interiors.

-- repo-doer@macmini
