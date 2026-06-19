# Complete core/shell final-basis and H1 payload helpers for multi-layer PQS.

"""
    pqs_multilayer_complete_core_shell_final_basis(plan; metadata = (;), ...)

Build the complete core/shell final-basis realization for a route-owned
multi-layer PQS shell source plan. This helper performs only the overlap
assembly needed by `CartesianFinalBasisRealization`; it does not materialize
one-body operators, H1, IDA, density-density, RHF, driver wiring, exports, or
artifacts.
"""
function pqs_multilayer_complete_core_shell_final_basis(
    plan;
    identity_atol::Real = 1.0e-8,
    rank_atol::Real = 1.0e-10,
    metadata = (;),
)
    _pqs_multilayer_property(plan, :object_kind) ===
        :pqs_multilayer_shell_source_plan ||
        throw(ArgumentError("PQS multi-layer final basis requires a pqs_multilayer_shell_source_plan"))
    if plan.status !== :available_pqs_multilayer_shell_source_plan
        return _blocked_pqs_multilayer_complete_core_shell_final_basis(
            plan;
            blocker = isnothing(plan.blocker) ?
                :pqs_multilayer_shell_source_plan_not_available :
                plan.blocker,
            metadata,
        )
    end

    metrics = plan.metrics
    core_overlap = _pqs_multilayer_support_product_matrix(
        plan.core_support_states,
        plan.core_support_states,
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    core_shell_overlap = _pqs_multilayer_support_product_matrix(
        plan.core_support_states,
        plan.shell_support_states,
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    shell_overlap = _pqs_multilayer_support_product_matrix(
        plan.shell_support_states,
        plan.shell_support_states,
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    merged_metadata = merge(
        NamedTuple(plan.metadata),
        NamedTuple(metadata),
        (;
            source = :pqs_multilayer_complete_core_shell_final_basis,
            input_source_plan = :pqs_multilayer_shell_source_plan,
            overlap_blocks_built_from_plan_support_states = true,
            h1_materialized = false,
            ida_data_materialized = false,
            density_density_materialized = false,
            rhf_materialized = false,
            driver_route_materialized = false,
            exports_materialized = false,
            artifacts_materialized = false,
        ),
    )
    final_basis =
        CartesianFinalBasisRealization.pqs_complete_core_shell_final_basis(
            core_support_indices = plan.core_support_indices,
            shell_support_indices = plan.shell_support_indices,
            core_overlap = core_overlap,
            core_shell_overlap = core_shell_overlap,
            shell_overlap = shell_overlap,
            shell_final_coefficients = plan.shell_final_coefficients,
            identity_atol = identity_atol,
            rank_atol = rank_atol,
            metadata = merged_metadata,
        )
    return merge(
        final_basis,
        (;
            source_plan_object_kind = plan.object_kind,
            source_plan_status = plan.status,
            source_plan_layer_count = plan.layer_count,
            source_plan_core_support_count = length(plan.core_support_indices),
            source_plan_shell_support_count = length(plan.shell_support_indices),
            source_plan_final_basis_helper =
                :pqs_multilayer_complete_core_shell_final_basis,
            h1_materialized = false,
            ida_data_materialized = false,
            density_density_materialized = false,
            rhf_materialized = false,
            driver_route_materialized = false,
            exports_materialized = false,
            artifacts_materialized = false,
        ),
    )
end

"""
    pqs_multilayer_complete_core_shell_h1_payload(plan; final_basis, ...)

Build the narrow complete core/shell PQS H1 assembly payload for an available
multi-layer source plan and final basis. This helper assembles only kinetic,
separated by-center electron-nuclear one-body matrices, the final H1
Hamiltonian, and the ordinary H1 solve. It does not materialize IDA,
density-density, RHF, GTO, driver wiring, exports, artifacts, or fixture-rule
policy.
"""
function pqs_multilayer_complete_core_shell_h1_payload(
    plan;
    final_basis,
    coulomb_expansion,
    center_records,
    axis_layers = nothing,
    gaussian_factor_terms_by_center = nothing,
    metadata = (;),
)
    _pqs_multilayer_property(plan, :object_kind) ===
        :pqs_multilayer_shell_source_plan ||
        throw(ArgumentError("PQS multi-layer H1 payload requires a pqs_multilayer_shell_source_plan"))
    plan.status === :available_pqs_multilayer_shell_source_plan ||
        throw(ArgumentError("PQS multi-layer H1 payload requires an available source plan"))
    get(final_basis, :status, nothing) === :available_pqs_complete_core_shell_final_basis ||
        throw(ArgumentError("PQS multi-layer H1 payload requires an available complete core/shell final basis"))

    support_kinetic = pqs_multilayer_support_kinetic_matrix(plan)
    support_nuclear_by_center =
        pqs_multilayer_support_electron_nuclear_by_center_matrices(
            plan;
            coulomb_expansion,
            center_records,
            axis_layers,
            gaussian_factor_terms_by_center,
        )
    final_kinetic =
        CartesianFinalBasisRealization.pqs_complete_core_shell_final_one_body_matrix(
            final_basis,
            support_kinetic;
            term = :kinetic,
        )
    raw_base_layer_gaussian_factor_matrices_used =
        isnothing(gaussian_factor_terms_by_center) && !isnothing(axis_layers)
    nuclear_factor_source =
        isnothing(gaussian_factor_terms_by_center) ?
        :centered_axis_layer_gaussian_factor_terms :
        :pgdg_intermediate_gaussian_factor_terms
    final_nuclear_by_center = map(support_nuclear_by_center.records) do record
        CartesianFinalBasisRealization.pqs_complete_core_shell_final_one_body_matrix(
            final_basis,
            record.support_operator;
            term = :electron_nuclear_by_center,
            center_record = record,
            metadata = merge(
                record.metadata,
                (;
                    nuclear_factor_source,
                    support_gaussian_factor_terms_source =
                        record.gaussian_factor_terms_source,
                    raw_base_layer_gaussian_factor_matrices_used,
                ),
            ),
        )
    end
    final_hamiltonian =
        CartesianFinalBasisRealization.pqs_complete_core_shell_final_one_electron_hamiltonian(
            final_kinetic,
            final_nuclear_by_center,
        )
    h1 = CartesianFinalBasisRealization.pqs_complete_core_shell_final_h1_solve(
        final_hamiltonian,
    )
    summary = (;
        status = :materialized_pqs_multilayer_complete_core_shell_h1_payload,
        blocker = nothing,
        final_dimension = h1.final_dimension,
        lowest_energy = h1.lowest_energy,
        solve_kind = h1.solve_kind,
        center_count = support_nuclear_by_center.center_count,
        final_basis_status = final_basis.status,
        source_plan_status = plan.status,
        support_kinetic_materialized = true,
        support_nuclear_by_center_materialized = true,
        final_one_body_transfer_materialized = true,
        hamiltonian_data_materialized = true,
        h1_solve_materialized = true,
        ida_data_materialized = false,
        density_density_materialized = false,
        rhf_materialized = false,
        gto_materialized = false,
        driver_route_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
        fixture_rule_policy_materialized = false,
    )
    return (;
        object_kind = :pqs_multilayer_complete_core_shell_h1_payload,
        status = summary.status,
        blocker = nothing,
        support_kinetic,
        support_nuclear_by_center,
        final_kinetic,
        final_nuclear_by_center,
        final_hamiltonian,
        h1,
        summary,
        metadata = merge(
            NamedTuple(metadata),
            (;
                source = :pqs_multilayer_complete_core_shell_h1_payload,
                support_space_one_body_helpers_scoped_to_h1 = true,
                ida_data_materialized = false,
                density_density_materialized = false,
                rhf_materialized = false,
                gto_materialized = false,
                driver_route_materialized = false,
                exports_materialized = false,
                artifacts_materialized = false,
                fixture_rule_policy_materialized = false,
            ),
        ),
    )
end

"""
    pqs_multilayer_complete_core_shell_h1_j_payload(plan; ...)

Build the narrow complete core/shell PQS H1/J diagnostic payload from an
available H1 payload and route-owned support-density inputs. This helper uses
the reviewed localized IDA density gauge. It does not materialize
RHF, GTO, driver wiring, exports, artifacts, or fixture-rule policy.
"""
function pqs_multilayer_complete_core_shell_h1_j_payload(
    plan;
    final_basis,
    h1_payload,
    axis_weights,
    raw_pair_factor_terms,
    coulomb_expansion,
    metadata = (;),
)
    _pqs_multilayer_property(plan, :object_kind) ===
        :pqs_multilayer_shell_source_plan ||
        throw(ArgumentError("PQS multi-layer H1/J payload requires a pqs_multilayer_shell_source_plan"))
    plan.status === :available_pqs_multilayer_shell_source_plan ||
        throw(ArgumentError("PQS multi-layer H1/J payload requires an available source plan"))
    get(final_basis, :status, nothing) === :available_pqs_complete_core_shell_final_basis ||
        throw(ArgumentError("PQS multi-layer H1/J payload requires an available complete core/shell final basis"))
    get(h1_payload, :object_kind, nothing) ===
        :pqs_multilayer_complete_core_shell_h1_payload ||
        throw(ArgumentError("PQS multi-layer H1/J payload requires a complete core/shell H1 payload"))
    get(h1_payload, :status, nothing) ===
        :materialized_pqs_multilayer_complete_core_shell_h1_payload ||
        throw(ArgumentError("PQS multi-layer H1/J payload requires a materialized H1 payload"))

    support_weights = pqs_multilayer_support_weights(plan; axis_weights)
    support_pair_raw = pqs_multilayer_support_pair_raw_numerator_matrix(
        plan;
        raw_pair_factor_terms,
        coulomb_expansion,
    )
    density_interaction =
        CartesianFinalBasisRealization.pqs_complete_core_shell_ida_density_interaction(
            final_basis,
            support_pair_raw,
            support_weights;
            metadata = merge(
                NamedTuple(metadata),
                (;
                    source = :pqs_multilayer_complete_core_shell_h1_j_payload,
                    support_density_input_source =
                        :pqs_multilayer_support_density_helpers,
                ),
            ),
        )
    if density_interaction.status !==
       :materialized_pqs_complete_core_shell_ida_density_interaction
        summary = (;
            status = :blocked_pqs_multilayer_complete_core_shell_h1_j_payload,
            blocker = density_interaction.blocker,
            final_dimension = get(final_basis, :final_retained_count, 0),
            h1_energy = h1_payload.h1.lowest_energy,
            density_interaction_status = density_interaction.status,
            self_coulomb_materialized = false,
            ida_data_materialized = true,
            density_density_materialized = true,
            rhf_materialized = false,
            gto_materialized = false,
            driver_route_materialized = false,
            exports_materialized = false,
            artifacts_materialized = false,
        )
        return (;
            object_kind = :pqs_multilayer_complete_core_shell_h1_j_payload,
            status = summary.status,
            blocker = summary.blocker,
            density_interaction,
            self_coulomb = nothing,
            summary,
            metadata = merge(
                NamedTuple(metadata),
                (;
                    source = :pqs_multilayer_complete_core_shell_h1_j_payload,
                    support_density_input_source =
                        :pqs_multilayer_support_density_helpers,
                ),
            ),
        )
    end

    hamiltonian = Matrix{Float64}(h1_payload.final_hamiltonian.hamiltonian_matrix)
    decomposition = eigen(Symmetric((hamiltonian + transpose(hamiltonian)) ./ 2))
    lowest_energy = first(decomposition.values)
    lowest_orbital = decomposition.vectors[:, 1]
    self_coulomb =
        CartesianFinalBasisRealization.pqs_complete_core_shell_ida_orbital_self_coulomb(
            density_interaction,
            lowest_orbital;
            metadata = merge(
                NamedTuple(metadata),
                (;
                    source = :pqs_multilayer_complete_core_shell_h1_j_payload,
                    h1_orbital_source = :h1_payload_final_hamiltonian_lowest_eigenvector,
                ),
            ),
        )
    summary = (;
        status = :materialized_pqs_multilayer_complete_core_shell_h1_j_payload,
        blocker = nothing,
        final_dimension = h1_payload.h1.final_dimension,
        h1_energy = h1_payload.h1.lowest_energy,
        h1_orbital_energy = lowest_energy,
        h1_energy_reconstruction_error = abs(lowest_energy - h1_payload.h1.lowest_energy),
        density_interaction_status = density_interaction.status,
        density_gauge = density_interaction.density_gauge,
        ida_weights_all_positive =
            density_interaction.ida_weights_all_positive,
        electron_electron_ida_finite =
            density_interaction.electron_electron_ida_finite,
        self_coulomb = self_coulomb.self_coulomb,
        support_density_input_source = :pqs_multilayer_support_density_helpers,
        h1_orbital_source = :h1_payload_final_hamiltonian_lowest_eigenvector,
        signed_final_weight_division_used = false,
        raw_no_division_used = false,
        density_normalized_pair_terms_used_as_authority = false,
        ida_data_materialized = true,
        density_density_materialized = true,
        rhf_materialized = false,
        gto_materialized = false,
        driver_route_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
    )
    return (;
        object_kind = :pqs_multilayer_complete_core_shell_h1_j_payload,
        status = summary.status,
        blocker = nothing,
        density_interaction,
        self_coulomb,
        summary,
        metadata = merge(
            NamedTuple(metadata),
            (;
                source = :pqs_multilayer_complete_core_shell_h1_j_payload,
                support_density_input_source =
                    :pqs_multilayer_support_density_helpers,
                density_gauge = :localized_ida,
                signed_final_weight_division_used = false,
                raw_no_division_used = false,
                density_normalized_pair_terms_used_as_authority = false,
            ),
        ),
    )
end

function _blocked_pqs_multilayer_complete_core_shell_final_basis(
    plan;
    blocker,
    metadata = (;),
)
    return (;
        object_kind = :pqs_multilayer_complete_core_shell_final_basis,
        status = :blocked_pqs_multilayer_complete_core_shell_final_basis,
        blocker,
        source_plan_object_kind = _pqs_multilayer_property(plan, :object_kind),
        source_plan_status = _pqs_multilayer_property(plan, :status),
        source_plan_blocker = _pqs_multilayer_property(plan, :blocker),
        source_plan_layer_count = _pqs_multilayer_property(plan, :layer_count, 0),
        source_plan_core_support_count =
            hasproperty(plan, :core_support_indices) ?
            length(plan.core_support_indices) : 0,
        source_plan_shell_support_count =
            hasproperty(plan, :shell_support_indices) ?
            length(plan.shell_support_indices) : 0,
        final_basis_materialized = false,
        one_body_operator_materialized = false,
        h1_materialized = false,
        ida_data_materialized = false,
        density_density_materialized = false,
        rhf_materialized = false,
        driver_route_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
        metadata = merge(
            hasproperty(plan, :metadata) ? NamedTuple(plan.metadata) : (;),
            NamedTuple(metadata),
            (;
                source = :pqs_multilayer_complete_core_shell_final_basis,
                input_source_plan = :pqs_multilayer_shell_source_plan,
                final_basis_helper = :pqs_complete_core_shell_final_basis,
            ),
        ),
    )
end
