# Compact old-seed oracle summary for White--Lindsey adapter validation.

const _WHITE_LINDSEY_ORACLE_ONE_BODY_TERMS = (
    :overlap,
    :position_x,
    :position_y,
    :position_z,
    :x2_x,
    :x2_y,
    :x2_z,
    :kinetic,
)

"""
    white_lindsey_materialized_seed_oracle_summary([report]; kwargs...)

Return compact metadata from the old White--Lindsey materialized seed for use
as a validation oracle. The zero-argument form builds the private seed report
and then discards large fixture/matrix objects from the returned summary.

This summary is oracle/reference metadata only. It is not route authority for
the new adapter and does not build new coefficient maps, adapter pair blocks,
Hamiltonian data, exports, artifacts, IDA/MWG data, or Coulomb.
"""
function white_lindsey_materialized_seed_oracle_summary(; kwargs...)
    seed_report_builder = getproperty(
        parentmodule(@__MODULE__),
        :_white_lindsey_low_order_materialized_seed_report,
    )
    return white_lindsey_materialized_seed_oracle_summary(
        seed_report_builder(; kwargs...),
    )
end

function white_lindsey_materialized_seed_oracle_summary(report)
    inventory = report.inventory
    route_units = report.route_units
    operator_inventory = report.operator_inventory
    retained_units = route_units.retained_units
    operator_terms = Tuple(operator_inventory.terms)

    return (;
        object_kind = :white_lindsey_materialized_seed_oracle_summary,
        status = :available_white_lindsey_materialized_seed_oracle_summary,
        oracle_role = :validation_oracle_only,
        route_authority = false,
        adapter_authority = false,
        private_development_only = true,
        seed_report_kind = report.object_kind,
        seed_report_status = report.status,
        route_family = report.route_family,
        shellization_source = report.shellization_source,
        packet_kernel = report.packet_kernel,
        packet_inventory_available = !isnothing(report.packet_kernel),
        retained_dimension = report.retained_dimension,
        retained_unit_count = length(retained_units),
        unit_keys = Tuple(unit.unit_key for unit in retained_units),
        unit_roles = Tuple(unit.unit_role for unit in retained_units),
        retained_unit_kinds =
            Tuple(unit.retained_unit_kind for unit in retained_units),
        retained_counts = route_units.unit_inventory.retained_counts,
        retained_ranges = route_units.unit_inventory.ranges,
        piece_counts = inventory.piece_counts,
        support_counts = inventory.support_counts,
        seed_retained_counts = inventory.retained_counts,
        seed_retained_ranges = inventory.retained_ranges,
        fixed_block_ready = inventory.fixed_block_ready,
        overlap_ready = inventory.overlap_ready,
        retained_basis_integral_weights_ready =
            inventory.retained_basis_integral_weights_ready,
        weight_semantics = report.weight_semantics,
        operator_inventory_available = true,
        operator_source = operator_inventory.operator_source,
        operator_terms,
        one_body_operator_matrix_available =
            _white_lindsey_seed_oracle_term_availability(operator_terms),
        fixed_block_operator_matrix_sizes = operator_inventory.matrix_sizes,
        operator_finite_ready = operator_inventory.finite_ready,
        all_operator_matrices_finite = operator_inventory.all_finite,
        operator_symmetric_ready = operator_inventory.symmetric_ready,
        overlap_identity_ready = operator_inventory.overlap_identity_ready,
        operator_pairs_materialized = report.operator_pairs_materialized,
        electron_electron_materialized = report.electron_electron_materialized,
        new_adapter_coefficient_maps_materialized = false,
        new_adapter_pair_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

function _white_lindsey_seed_oracle_term_availability(operator_terms)
    return NamedTuple{_WHITE_LINDSEY_ORACLE_ONE_BODY_TERMS}(
        Tuple(term in operator_terms for term in _WHITE_LINDSEY_ORACLE_ONE_BODY_TERMS),
    )
end
