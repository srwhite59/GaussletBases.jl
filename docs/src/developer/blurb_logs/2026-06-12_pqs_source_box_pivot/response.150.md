Pass 150 response - audit Be2 PQS post-H1 Hamiltonian seam

Scope:
- No source edits.
- No commits.
- No Julia commands.
- Read-only inspection only.

H1-J diagnostic contract summary:
- Current one-center H1-J entrypoint is `pqs_multilayer_complete_core_shell_h1_j_payload(plan; ...)`.
  - It hard-requires `object_kind === :pqs_multilayer_shell_source_plan`.
  - It hard-requires `status === :available_pqs_multilayer_shell_source_plan`.
  - It hard-requires an available complete core/shell final basis.
  - It hard-requires an old `:pqs_multilayer_complete_core_shell_h1_payload`.
  - References: `src/pqs_multilayer_complete_core_shell_h1.jl:233`.
- Required density inputs are:
  - `axis_weights`;
  - `raw_pair_factor_terms`;
  - `coulomb_expansion`;
  - support states from the old source plan.
  - References: `src/pqs_multilayer_complete_core_shell_h1.jl:256`, `src/pqs_multilayer_support_density.jl:23`, and `src/pqs_multilayer_support_density.jl:66`.
- The support-density helpers also hard-gate on `:pqs_multilayer_shell_source_plan`.
  - Reference: `src/pqs_multilayer_support_density.jl:3`.
- Density gauge:
  - `pqs_complete_core_shell_pre_final_density_interaction(...)` uses `:pre_final_localized_positive_weight`;
  - it computes `pre_final_weights = pre_final_coefficients' * support_weights`;
  - it applies the old fixed-block-style positive-weight boundary as `weighted_coefficients = pre_final_coefficients ./ pre_final_weights`;
  - it builds a pre-final pair matrix from a raw support-row numerator.
  - References: `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl:527` and `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl:547`.
- Axis weights/raw pair factors:
  - driver density inputs come from `CartesianContractedParentMetrics._pqs_source_box_ida_factor_provenance(bundles; expected_term_count)`;
  - they are explicitly marked private diagnostic-only and do not use retained PQS weights as authority.
  - Reference: `src/pqs_source_box_route_driver_helpers.jl:10841`.
- What scalar it reports:
  - H1-J diagonalizes the final H1 matrix, takes the lowest H1 orbital, maps it into the pre-final gauge, and reports one restricted-orbital direct-minus-exchange self-Coulomb scalar.
  - References: `src/pqs_multilayer_complete_core_shell_h1.jl:311` and `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl:698`.
- Why it is diagnostic/private:
  - It reports a one-orbital scalar and convention facts, not a complete solver-facing Hamiltonian.
  - It does not run RHF/GTO/driver/export/artifact code.
  - It is tied to a pre-final density gauge; current fixture policy explicitly classifies compact H1/J as convention diagnostic, not endpoint validation.
  - References: `src/pqs_multilayer_complete_core_shell_h1.jl:327` and `docs/src/developer/pqs_source_box_fixture_policy.md:13`.

Downstream Ham-constructor needs visible from current code:
- A downstream solver-facing object needs at least:
  - final/retained basis description and ordering;
  - one-body matrix `H1`;
  - an electron-electron representation, usually a density-density `Vee`-like object or an explicit factor payload;
  - nuclear metadata and charge convention;
  - electron count/occupation metadata when a solver handoff is intended;
  - ordering and producer/source metadata.
- Existing export docs describe current solver-facing exports as producer-side payloads with explicit `H1`, `Vee`, ordering metadata, and manifest data.
  - References: `docs/src/reference/export.md:9`, `docs/src/reference/export.md:20`, and `docs/src/reference/export.md:36`.
- The dense IDA export shape is concrete: payload keys include `H1`, `Vee`, `dims_per_shell`, `orders`, and `basis_centers`; metadata records interaction model and ordering.
  - Reference: `src/fullida_dense_export.jl:192`.
- Existing private complete-core/shell Ham payload code expects:
  - available final basis;
  - materialized H1 payload and one-body Hamiltonian;
  - available density inputs;
  - materialized H1-J diagnostic payload;
  - materialized pre-final density interaction.
  - References: `src/pqs_source_box_route_driver_helpers.jl:11263` and `src/pqs_source_box_route_driver_helpers.jl:11310`.
- When materialized, that Ham payload labels the electron-electron representation as `:pre_final_density_interaction`, records `density_gauge`, `raw_pair_factor_convention`, support row order, dimension summary, ordering summary, and keeps public/export/artifact/RHF claims false.
  - Reference: `src/pqs_source_box_route_driver_helpers.jl:11372`.

What Be2/PQS already has:
- Private diatomic source plan:
  - carries parent axis bundles, metrics, core/shell support indices/states, shell coefficients, support order, retained/pre-final map, source summaries, and convention labels.
  - Reference: `src/pqs_source_box_route_driver_helpers.jl:11741`.
- Private diatomic final basis:
  - built directly from private diatomic source-plan support states/metrics and the lower final-basis helper;
  - records final support row order `:core_then_shell`.
  - Reference: `src/pqs_source_box_route_driver_helpers.jl:13092`.
- Private diatomic H1 payload:
  - carries source-plan status, final-basis status, support kinetic/nuclear statuses, final one-body statuses, final H1 Hamiltonian, H1 solve status, final dimension, energy, and nonclaim flags.
  - Reference: `src/pqs_source_box_route_driver_helpers.jl:11777`.
- Parent axis bundle and Coulomb expansion access:
  - H1 already uses `source_plan.bundles`, parent center records, axis layers, and `coulomb_gaussian_expansion(doacc = false)`.
  - Reference: `src/pqs_source_box_route_driver_helpers.jl:13464`.
- Current readiness:
  - H1 available advances Be2/PQS readiness to `:missing_diatomic_complete_core_shell_h1_j_consumer`.
  - Reference: `src/pqs_source_box_route_driver_helpers.jl:13759`.

What remains missing:
- Densities/two-electron inputs have not been made diatomic-private:
  - no diatomic support weights payload;
  - no diatomic raw pair-factor provenance payload;
  - no diatomic support raw pair numerator payload;
  - no diatomic pre-final density interaction payload;
  - no solver-facing Ham-input/handoff object.
- Existing density-input helpers are still old one-center object-kind gated.
  - References: `src/pqs_source_box_route_driver_helpers.jl:10872` and `src/pqs_multilayer_support_density.jl:3`.
- No downstream compatibility decision exists for whether the first Be2/PQS handoff should be dense `H1`/`Vee`, factorized density inputs, or a current-route comparison-only object.
- No H1-J scalar is needed to define a Hamiltonian constructor; it is useful as a convention diagnostic once the density interaction exists.

Recommendation:
- Do not implement private H1-J diagnostic first.
- Implement a private Be2/PQS Ham-input/electron-electron payload first.
- Suggested name:
  - `_PQSDiatomicCompleteCoreShellHamInputPayload`
  - helper: `_pqs_source_box_route_driver_diatomic_complete_core_shell_ham_input_payload`
- Rationale:
  - H1-J would mostly add the self-Coulomb scalar after building the same density interaction that a Hamiltonian handoff actually needs.
  - The downstream WL/PQS comparison needs a structured electron-electron representation plus H1 and ordering/convention metadata, not just a scalar.
  - A Ham-input payload can carry the reusable density-interaction object without claiming public export, HamV6 compatibility, RHF readiness, or endpoint physics.

Smallest safe next implementation pass:
- Add a private Ham-input/electron-electron payload that consumes:
  - private diatomic source plan;
  - private final basis;
  - private H1 payload;
  - Coulomb expansion;
  - source-plan bundles for axis weights/raw pair factors.
- Build only:
  - support weights;
  - support raw pair numerator;
  - pre-final density interaction through `CartesianFinalBasisRealization.pqs_complete_core_shell_pre_final_density_interaction(...)`;
  - compact ordering/convention metadata;
  - references to the H1 matrix/final Hamiltonian already materialized by the private H1 payload.
- Do not build H1-J self-Coulomb in that first pass.
- Do not claim full Ham payload/export/public/RHF/CR2/hfdmrg readiness.
- After that, a tiny H1-J diagnostic pass can be optional and should consume the Ham-input payload's density interaction rather than rebuild density data.

Proposed status/blocker labels:
- Payload statuses:
  - `:blocked_diatomic_complete_core_shell_ham_input_payload`
  - `:available_diatomic_complete_core_shell_ham_input_payload`
- Missing/blocker labels:
  - `:missing_diatomic_complete_core_shell_source_plan`
  - `:missing_diatomic_complete_core_shell_final_basis`
  - `:missing_diatomic_complete_core_shell_h1_payload`
  - `:missing_diatomic_complete_core_shell_density_inputs`
  - `:missing_diatomic_support_density_provenance`
  - `:missing_diatomic_support_weights`
  - `:missing_diatomic_raw_pair_factor_terms`
  - `:blocked_diatomic_pre_final_density_interaction`
  - `:missing_diatomic_ham_consumer_contract`
- Readiness after available Ham-input payload should advance to something like:
  - `:missing_diatomic_ham_consumer_contract`
  - or `:missing_diatomic_complete_core_shell_ham_handoff_consumer`
  depending on whether the manager wants a comparison-only object or a downstream handoff boundary next.

Git status:
```text
## main...origin/main
```

Deletion/shrinkage forecast:
- deleted: none in this no-edit audit.
- simplified: a Ham-input payload would make the current H1-J blocker sharper by separating reusable electron-electron/Ham inputs from the optional H1-J scalar diagnostic.
- quarantined: H1-J self-Coulomb scalar as diagnostic-only; old `:pqs_multilayer_shell_source_plan` density helpers; public API/export/artifact/HamV6/hfdmrg/CR2/RHF claims.
- not deleted because: the one-center H1-J path remains the validated convention diagnostic, and no diatomic electron-electron payload has been approved yet.
- exact remaining caller/blocker: Be2/PQS readiness currently blocks at `:missing_diatomic_complete_core_shell_h1_j_consumer`; the next more useful blocker should be a private Ham-input/electron-electron payload, then a consumer-contract blocker.

-- repo-doer@macmini
