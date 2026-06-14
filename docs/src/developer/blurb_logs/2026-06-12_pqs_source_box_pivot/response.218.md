Pass 218 response - H2 gausslet-only WL reference comparison path audit

Role/task boundary:
- repo-doer@macmini, read-only comparison-path audit.
- No `src`, `test`, or `bin` edits.
- No Julia probes or broad tests run.

Files/docs inspected:
- `.agent_handoffs/blurb.218.md`
- `.agent_handoffs/state.md`
- `.agent_handoffs/RUN.md`
- `.agent_handoffs/DOER_STARTUP.md`
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/high_order_endcap_panel_h2_chemistry_reproduction_2026-05-16.md`
- `docs/diatomic_h2_completion_state.md`
- `docs/src/reference/atomic_and_ordinary.md`
- `docs/src/developer/high_order_mainline_import_readiness_2026-05-15.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.191.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.200.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.210.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.192.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.193.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.208.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.209.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.210.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.212.md`
- `test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
- `bin/cartesian_ham_builder.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `src/pqs_source_box_low_order_materialization.jl`

Old WL/QW gausslet-only reference source and values found:
- The old gausslet-only H2 R=4 fixed block is the 463-column block underlying the documented WL/QW H2 runs.
- The inventory accepted in prior review is:
  - parent axis lengths: `(9, 9, 15)`, parent size `1215`;
  - fixed gausslet-only block shape: `(1215, 463)`;
  - retained order: `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`;
  - retained counts: `(251, 98, 114)`;
  - support counts: `(275, 578, 362)`;
  - child/core support: `5 x 5 x 11`, i.e. two `5^3` atom cores plus the 25-row contact plane;
  - no separate midpoint/product retained unit in the 463-column gausslet-only block.
- The 463 fixed block itself is supplement-free as a gausslet-only block.
- The old published/scalar WL/QW chemistry rows are not supplement-free endpoints: the default complete-rectangular row used final dimension `481 = 463 + 18` after adding H/cc-pVTZ S/P residual supplement columns. The endcap/panel row likewise carried 18 residual columns with final dimension 461.

Supplement status of old references:
- The old scalar WL/QW HF/ED values are supplemented and are therefore forbidden as direct comparison values for the current H2 gausslet-only PQS endpoint.
- In particular, the documented HF/ED totals from the high-order endcap-panel reproduction note include residual H/cc-pVTZ molecular Gaussian supplement columns. They may remain historical chemistry references, but they are not an honest no-supplement gausslet-only comparison target.

Current visible-driver/materializer path:
- The current H2 PQS visible-driver input already expresses the matching physical target geometry and inventory:
  - H2 centers at `(0, 0, -2)` and `(0, 0, 2)`;
  - bond length `4.0`;
  - `q = n_s = 5`;
  - `supplement_policy = :none`;
  - support counts `(275, 578, 362)`;
  - retained counts `(251, 98, 114)`;
  - final dimension `463`.
- That path is the current PQS physical endpoint, not a matching WL reference route.
- The visible driver has low-order route knobs and a private route-configured diatomic atom-growth materializer seam under `route_family = :white_lindsey_low_order`.
- I did not find a current driver input/test/artifact path proving that `:white_lindsey_low_order` can be constrained to the exact H2 R=4, no-supplement, final-dimension-463, common-support inventory target.
- The private low-order materializer reports aggregate retained/support dimensions and readiness flags. As currently audited, it is not enough by itself for endpoint comparison because it does not establish the exact support split, retained split/order, retained-transform kind, and no-supplement comparison label for the H2 463 target.

Required comparability conditions:
- same geometry: H2 R=4 with centers equivalent to `(0, 0, -2)` and `(0, 0, 2)`;
- same physical shell/support plan: support counts `(275, 578, 362)` and retained order `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`;
- no supplement: no GTO/MWG residual columns or supplement residual policy;
- final dimension: exactly `463`;
- explicit retained-transform kind: WL/old-QW gausslet retained transform must be named separately from the PQS retained transform;
- endpoint comparison label: something explicit such as `WL/QW H2 R=4 gausslet-only 463`.

Can current route-configured WL materialization be constrained to the match?
- Not yet proven.
- Shared `q=5` and `n_s=5` are not sufficient.
- The current low-order materializer seam may be useful, but until it produces or records the exact H2 463 common support/retained inventory, no supplement policy, and retained-transform kind, it should remain a candidate seam rather than the accepted WL comparison route.

Recommended next pass:
- Do not compare scalar energies yet.
- First add a small no-matrix WL gausslet-only reference-candidate/readiness seam for the visible driver.
- The seam should be route-owned enough to report:
  - geometry fingerprint;
  - parent axis counts `(9, 9, 15)`;
  - support counts `(275, 578, 362)`;
  - retained counts `(251, 98, 114)`;
  - retained order `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`;
  - final dimension `463`;
  - supplement policy `:none`;
  - explicit retained-transform kind for WL/old-QW;
  - endpoint comparison label `WL/QW H2 R=4 gausslet-only 463`;
  - blocker on any mismatch.
- If that readiness seam can be attached to `cartesian_report`/materialization without building matrices, the following pass can decide whether the actual WL one-body/two-body/RHF comparison is available.

Exact blocker:
- `:missing_route_configured_wl_h2_gausslet_only_463_reference_candidate`
- Secondary forbidden-reference blocker for the old scalar values:
  `:supplemented_wl_qw_h2_reference_not_gausslet_only`

Source/test/bin scoped line budget:
- No `src`, `test`, or `bin` files changed.
- Line-count rule not applicable.

Validation:
- `git status --short --branch`
  - `## main...origin/main`
- No Julia command run; this was a read-only audit.

Deletion/shrinkage report:
- deleted: none
- simplified: none
- quarantined: old supplemented WL/QW H2 HF/ED scalar references are explicitly quarantined as non-comparable for the no-supplement endpoint
- not deleted because: audit-only pass; no source/test/bin edits requested
- exact remaining caller/blocker: `:missing_route_configured_wl_h2_gausslet_only_463_reference_candidate`

-- repo-doer@macmini
