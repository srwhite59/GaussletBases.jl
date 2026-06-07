# Narrow White--Lindsey boundary-stratum unit coefficient adapters.

"""
    white_lindsey_boundary_stratum_unit_coefficients(unit_or_descriptor)

Materialize the narrow unit coefficient map supported by the current LW adapter
checkpoint. This pass supports only corner strata as a support-local one-row,
one-column map with value `1.0`.

Facet and edge strata block explicitly. This helper does not call old
White--Lindsey kernels, does not build one-body pair blocks, and does not build
Hamiltonian data, exports, artifacts, IDA/MWG data, or Coulomb.
"""
function white_lindsey_boundary_stratum_unit_coefficients(unit_or_descriptor)
    descriptor =
        _white_lindsey_is_unit_adapter_descriptor(unit_or_descriptor) ?
        unit_or_descriptor :
        _white_lindsey_has_retained_unit_fields(unit_or_descriptor) ?
        white_lindsey_boundary_stratum_unit_adapter_descriptor(unit_or_descriptor) :
        nothing

    isnothing(descriptor) && return _white_lindsey_blocked_unit_coefficients(
        unit_or_descriptor,
        :not_white_lindsey_boundary_stratum_unit_or_descriptor,
    )

    return _white_lindsey_boundary_stratum_unit_coefficients_from_descriptor(
        descriptor,
    )
end

function _white_lindsey_boundary_stratum_unit_coefficients_from_descriptor(
    descriptor,
)
    status, blocker =
        _white_lindsey_unit_coefficients_status(descriptor)
    status === :materialized_white_lindsey_corner_unit_coefficients ||
        return _white_lindsey_unit_coefficients_result(
            descriptor,
            status,
            blocker,
            nothing,
        )

    coefficient_matrix = reshape([1.0], 1, 1)
    return _white_lindsey_unit_coefficients_result(
        descriptor,
        status,
        blocker,
        coefficient_matrix,
    )
end

function _white_lindsey_unit_coefficients_status(descriptor)
    _white_lindsey_descriptor_property(descriptor, :status) ===
        :available_metadata_only_white_lindsey_unit_adapter_descriptor || return (
        :blocked_white_lindsey_boundary_stratum_unit_coefficients,
        _white_lindsey_descriptor_property(
            descriptor,
            :blocker,
            :unavailable_white_lindsey_unit_adapter_descriptor,
        ),
    )
    _white_lindsey_descriptor_property(descriptor, :stratum_kind) ===
        :corner_cpb || return (
        :blocked_white_lindsey_boundary_stratum_unit_coefficients,
        :white_lindsey_unit_coefficients_not_implemented_for_stratum,
    )
    _white_lindsey_descriptor_property(descriptor, :planned_old_kernel) ===
        :_nested_corner_piece || return (
        :blocked_white_lindsey_boundary_stratum_unit_coefficients,
        :white_lindsey_corner_kernel_plan_not_available,
    )
    _white_lindsey_descriptor_property(descriptor, :source_cpb_shape) == (1, 1, 1) ||
        return (
            :blocked_white_lindsey_boundary_stratum_unit_coefficients,
            :white_lindsey_corner_source_cpb_not_support_local,
        )
    fixed_coordinates =
        _white_lindsey_descriptor_property(descriptor, :fixed_axis_coordinates)
    (!isnothing(fixed_coordinates) && length(fixed_coordinates) == 3) || return (
        :blocked_white_lindsey_boundary_stratum_unit_coefficients,
        :missing_white_lindsey_corner_fixed_coordinates,
    )
    return :materialized_white_lindsey_corner_unit_coefficients, nothing
end

function _white_lindsey_unit_coefficients_result(
    descriptor,
    status::Symbol,
    blocker,
    coefficient_matrix,
)
    materialized = !isnothing(coefficient_matrix)
    return (;
        object_kind = :white_lindsey_boundary_stratum_unit_coefficients,
        status,
        blocker,
        unit_key = _white_lindsey_descriptor_property(descriptor, :unit_key),
        unit_index = _white_lindsey_descriptor_property(descriptor, :unit_index),
        unit_kind = _white_lindsey_descriptor_property(descriptor, :unit_kind),
        stratum_kind =
            _white_lindsey_descriptor_property(descriptor, :stratum_kind),
        source_cpb_role =
            _white_lindsey_descriptor_property(descriptor, :source_cpb_role),
        source_cpb_shape =
            _white_lindsey_descriptor_property(descriptor, :source_cpb_shape),
        fixed_axis_coordinates =
            _white_lindsey_descriptor_property(
                descriptor,
                :fixed_axis_coordinates,
            ),
        planned_old_kernel =
            _white_lindsey_descriptor_property(descriptor, :planned_old_kernel),
        coefficient_space = :source_cpb_support_local,
        parent_row_indices_available = false,
        source_support_row_count = materialized ? 1 : nothing,
        retained_column_count = materialized ? 1 : nothing,
        source_support_row_index = materialized ? 1 : nothing,
        retained_column_index = materialized ? 1 : nothing,
        coefficient_matrix,
        nonzero_count = materialized ? 1 : 0,
        nonzero_values = materialized ? (1.0,) : (),
        coefficient_maps_materialized = materialized,
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

function _white_lindsey_blocked_unit_coefficients(input, blocker::Symbol)
    descriptor = (;
        unit_key = _white_lindsey_descriptor_property(input, :unit_key),
        unit_index = _white_lindsey_descriptor_property(input, :unit_index),
        unit_kind = _white_lindsey_descriptor_property(input, :unit_kind),
        stratum_kind = _white_lindsey_descriptor_property(input, :stratum_kind),
        source_cpb_role =
            _white_lindsey_descriptor_property(input, :source_cpb_role),
        source_cpb_shape =
            _white_lindsey_descriptor_property(input, :source_cpb_shape),
        fixed_axis_coordinates =
            _white_lindsey_descriptor_property(input, :fixed_axis_coordinates),
        planned_old_kernel =
            _white_lindsey_descriptor_property(input, :planned_old_kernel),
    )
    return _white_lindsey_unit_coefficients_result(
        descriptor,
        :blocked_white_lindsey_boundary_stratum_unit_coefficients,
        blocker,
        nothing,
    )
end

function _white_lindsey_is_unit_adapter_descriptor(value)
    return _white_lindsey_descriptor_property(value, :object_kind) ===
           :white_lindsey_boundary_stratum_unit_adapter_descriptor
end

function _white_lindsey_has_retained_unit_fields(value)
    return hasproperty(value, :unit_key) &&
           hasproperty(value, :unit_kind) &&
           hasproperty(value, :source_cpbs) &&
           hasproperty(value, :metadata)
end

function _white_lindsey_descriptor_property(value, key::Symbol, default = nothing)
    hasproperty(value, key) && return getproperty(value, key)
    return default
end
