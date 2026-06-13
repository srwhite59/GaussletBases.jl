Pass 169 read-only audit complete.

Purpose:
Audit the White-Lindsey population seam for the Be2 CR2 inspection artifact without editing source, tests, docs, generator, or generated artifacts.

Inventory:
- Existing PQS artifact generator: `tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl` currently builds only PQS route payloads and marks `routes/white_lindsey/status = unavailable`.
- Legacy/private WL seed surface: `_white_lindsey_low_order_materialized_seed_fixture`, `_white_lindsey_low_order_materialized_seed_report`, `_WhiteLindseyLowOrderHamBundleAdapter`, and related one-center route-driver materialization. This can produce 223 x 223 matrices but is a one-center/seed/debug route, not an honest Be2 comparator.
- Route-configured diatomic WL surface: `_cartesian_shellization_route_materialize_bond_aligned_diatomic`, `_pqs_source_box_route_driver_route_configured_diatomic_basis_adapter`, `_pqs_source_box_route_driver_route_configured_diatomic_ham_adapter`, and the builder path tested by `test/nested/cartesian_ham_builder_diatomic_config_smoke_runtests.jl`. This can produce Be2-relevant arrays.
- Route-configured diatomic atom-growth WL surface: `_cartesian_shellization_route_materialize_bond_aligned_diatomic` plus `low_order_shellization_policy = :atom_growth_complete_rectangular`, `_pqs_source_box_route_driver_diatomic_atom_growth_materialization`, and `_pqs_source_box_route_driver_route_configured_diatomic_ham_adapter`. This is the smallest honest seam for the next artifact population pass because the report explicitly identifies the atom-growth route, basis adapter, ham adapter, shellification materialization kind, and blocked/default-behavior facts.

Answers:
1. The WL path that can produce Be2-relevant arrays is the route-configured diatomic atom-growth route, with `low_order_shellization_policy = :atom_growth_complete_rectangular`, consumed through the route driver materialization. It produces a basis bundle and a ham bundle when `save_ham_artifact = true`, `white_lindsey_expansion` is supplied, and `route_configured_diatomic_ham_interaction_treatment = :ggt_nearest`.
2. That path is diatomic Be2-comparable in the limited sense of final-basis ordinary Cartesian/Qiu-White arrays for the same Be2 nuclear charges and locations. It is not PQS source-box comparable at the pre-final/source-box level. The old materialized seed and one-center route-configured WL bundle are not Be2-comparable; they are one-center/seed/debug fixtures even when the surrounding metadata records a Be2 shellization request.
3. Schema population forecast:
   - `system`: populate from the route driver config and ham export metadata: atom symbols, nuclear charges, atom locations, bond axis, electron count if the generator already owns that convention, units/status fields, and nuclear term storage.
   - `final_basis`: populate from `basis/final_dimension`, `basis/final_integral_weights`, `basis/basis_labels`, `basis/basis_centers`, basis kind, materialized report kind, shellization source/authority, and atom-growth consumption flags.
   - `one_body`: populate `ham/overlap`, `ham/one_body_hamiltonian`, optional `ham/kinetic_one_body`, `ham/nuclear_one_body_by_center/*`, and status fields from `_pqs_source_box_route_driver_route_configured_diatomic_ham_adapter_summary`.
   - `two_body`: populate `ham/interaction_matrix` only as a final-basis density-density interaction matrix with `interaction_model = density_density` and `interaction_treatment = ggt_nearest`. Do not invent PQS-style pre-final/source-box WL pair fields.
   - `validation`: finite/symmetry/shape checks, final dimension consistency, integral-weight positivity, nuclear metadata availability, interaction treatment consumed, and route/ham adapter status checks.
   - `metadata`: generator provenance, route family `white_lindsey_low_order`, `private_development_only`, materializer backend/nside/d/reference/tail spacing, low-order shellization policy requested/resolved/source/status, shellification materialization kind, route default behavior unchanged, and explicit unsupported placeholders.
4. Two-body representation should be labeled as a final-basis density-density interaction matrix from `OrdinaryCartesianOperators3D`, not as PQS pre-final pair data.
5. If population cannot be run honestly in the generator, keep WL unavailable with a blocker such as `:missing_route_configured_diatomic_atom_growth_ham_payload` or the more specific route status returned by the materialization, for example `:pending_route_configured_diatomic_mwg_operator_support` if MWG/IDA is requested.
6. Richer provenance can come from compact route objects already present: materialization status, materialized report kind, atom-growth report, basis/ham adapter summaries, shellization source/authority, materializer config/status, low-order shellization policy, and ham export statuses. Avoid copying the existing field cloud wholesale into the artifact; carry a compact WL route summary/fingerprint.
7. Supplement/residual-GTO, correction/EGOI, and QW atom-local HF should remain explicit unavailable metadata. `src/ordinary_qw_raw_blocks.jl` still throws for true active `l > 0` atomic supplements on the ordinary QW path, MWG/IDA route-configured diatomic ham support is explicitly pending, and the audited path does not provide correction/EGOI or QW atom-local HF solver inputs.
8. Smallest implementation seam: add a local generator helper that invokes the existing route-driver materialization for Be2 WL with `low_order_shellization_policy = :atom_growth_complete_rectangular`, captures the in-memory materialization and/or its ham bundle values, and maps only the existing final-basis arrays plus compact adapter/materialization summaries into `routes/white_lindsey`.
9. Deletion/shrink candidate for the next line-budget implementation: keep the generator mapping small enough that old seed/one-center WL helper pressure does not grow. Do not adapt the seed report into the artifact. A later cleanup can remove or quarantine artifact-facing references to seed/one-center WL paths once the diatomic atom-growth population is in place.

Recommended first implementation pass:
- In the local artifact generator only, add a private WL population helper for `routes/white_lindsey` that calls the existing Be2 route-driver materialization with atom-growth policy and `:ggt_nearest`.
- Store only final-basis arrays and compact route summaries. Mark unavailable fields explicitly instead of creating placeholder array conventions.
- Validate with a focused generator run and JLD2 readback. Do not touch source, tests, production route behavior, lattice size, checkpointing, caching, RHF/SCF/Fock, GTO, IDA/MWG, exports, artifacts outside the local generator outputs, or broad tests.

Read-only commands run:
- `git status --short --branch`
- `sed -n ...` on `AGENTS.md`, `BlurbStyle.md`, the CR2 review, the generator script, `src/cartesian_bundle_export.jl`, `src/cartesian_shellization_route.jl`, `src/cartesian_atom_growth_route_driver_helpers.jl`, `src/pqs_source_box_route_driver_helpers.jl`, `src/ordinary_qw_raw_blocks.jl`, and the named focused tests.
- `rg -n ...` on the named WL source/test surfaces.
- One read-only `rg` command included a nonexistent `src/cartesian_ham_builder.jl` path and exited 2 after returning useful matches from other files; no files were changed.

Validation:
- No Julia commands run.
- Read-only inspection only.
- Final git status: `## main...origin/main`

Deletion/shrinkage report:
deleted: none
simplified: none; no-edit audit only
quarantined: recommended keeping old materialized seed, one-center WL, supplement/residual-GTO, correction/EGOI, QW atom-local HF, and MWG/IDA as unavailable/not-comparable in the artifact
not deleted because: this pass was read-only and the live Be2-comparable WL route still needs a small generator population seam first
exact remaining caller/blocker: `tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl` has no WL materialization call yet; honest population should use route-configured diatomic atom-growth materialization plus `_pqs_source_box_route_driver_route_configured_diatomic_ham_adapter` with `:ggt_nearest`

-- repo-doer@macmini
