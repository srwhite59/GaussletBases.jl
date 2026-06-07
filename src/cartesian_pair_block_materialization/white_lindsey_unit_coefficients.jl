# Narrow White--Lindsey boundary-stratum unit coefficient adapters.

"""
    white_lindsey_boundary_stratum_unit_coefficients(unit_or_descriptor)

Materialize the narrow unit coefficient map supported by the current LW adapter
checkpoint. This pass supports corner strata as a support-local one-row,
one-column map with value `1.0`, plus edge and facet strata only when the
metadata context has real old-kernel inputs.

This helper does not build one-body pair blocks, Hamiltonian data, exports,
artifacts, IDA/MWG data, or Coulomb.
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

"""
    white_lindsey_boundary_stratum_unit_coefficient_context(unit_or_descriptor)

Return metadata-only old-kernel input context for one White--Lindsey
boundary-stratum unit. This reports the facts needed by future edge/facet
coefficient adapters without calling old kernels or materializing those
coefficient maps.
"""
function white_lindsey_boundary_stratum_unit_coefficient_context(unit_or_descriptor)
    descriptor =
        _white_lindsey_is_unit_adapter_descriptor(unit_or_descriptor) ?
        unit_or_descriptor :
        _white_lindsey_has_retained_unit_fields(unit_or_descriptor) ?
        white_lindsey_boundary_stratum_unit_adapter_descriptor(unit_or_descriptor) :
        nothing

    isnothing(descriptor) && return _white_lindsey_unit_coefficient_context(
        (;),
        :blocked_unknown_white_lindsey_unit_coefficient_context,
        :not_white_lindsey_boundary_stratum_unit_or_descriptor,
        (),
    )

    status, blocker, planned_old_calls =
        _white_lindsey_unit_coefficient_context_status(descriptor)
    return _white_lindsey_unit_coefficient_context(
        descriptor,
        status,
        blocker,
        planned_old_calls,
    )
end

function _white_lindsey_boundary_stratum_unit_coefficients_from_descriptor(
    descriptor,
)
    status, blocker =
        _white_lindsey_unit_coefficients_status(descriptor)
    status === :materialized_white_lindsey_facet_unit_coefficients &&
        return _white_lindsey_facet_unit_coefficients_result(descriptor)
    status === :materialized_white_lindsey_edge_unit_coefficients &&
        return _white_lindsey_edge_unit_coefficients_result(descriptor)
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
    stratum_kind = _white_lindsey_descriptor_property(descriptor, :stratum_kind)
    if stratum_kind in (:facet_cpb, :face_cpb)
        context =
            white_lindsey_boundary_stratum_unit_coefficient_context(descriptor)
        if context.status ===
           :ready_white_lindsey_facet_kernel_context_not_materialized
            _white_lindsey_facet_context_sources_materializable(context) &&
                return :materialized_white_lindsey_facet_unit_coefficients,
                       nothing
            return (
                :blocked_white_lindsey_boundary_stratum_unit_coefficients,
                :white_lindsey_facet_doside_source_not_materializable,
            )
        end
        return (
            :blocked_white_lindsey_boundary_stratum_unit_coefficients,
            :incomplete_white_lindsey_edge_facet_kernel_input_context,
        )
    elseif stratum_kind === :edge_cpb
        context =
            white_lindsey_boundary_stratum_unit_coefficient_context(descriptor)
        if context.status === :ready_white_lindsey_edge_kernel_context_not_materialized
            _white_lindsey_edge_context_source_materializable(context) &&
                return :materialized_white_lindsey_edge_unit_coefficients, nothing
            return (
                :blocked_white_lindsey_boundary_stratum_unit_coefficients,
                :white_lindsey_edge_doside_source_not_materializable,
            )
        end
        return (
            :blocked_white_lindsey_boundary_stratum_unit_coefficients,
            :incomplete_white_lindsey_edge_facet_kernel_input_context,
        )
    end
    stratum_kind === :corner_cpb || return (
        :blocked_white_lindsey_boundary_stratum_unit_coefficients,
        :unknown_white_lindsey_unit_coefficients_stratum,
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

function _white_lindsey_facet_unit_coefficients_result(descriptor)
    context =
        white_lindsey_boundary_stratum_unit_coefficient_context(descriptor)
    context.status === :ready_white_lindsey_facet_kernel_context_not_materialized ||
        return _white_lindsey_unit_coefficients_result(
            descriptor,
            :blocked_white_lindsey_boundary_stratum_unit_coefficients,
            :incomplete_white_lindsey_edge_facet_kernel_input_context,
            nothing,
        )

    active_entries = _white_lindsey_facet_active_axis_entries(context)
    first_entry, second_entry = active_entries
    first_source, first_source_policy = _white_lindsey_doside_source_for_axis(
        context.doside_source_1d,
        first_entry.axis,
    )
    second_source, second_source_policy = _white_lindsey_doside_source_for_axis(
        context.doside_source_1d,
        second_entry.axis,
    )
    first_retained, retained_policy =
        _white_lindsey_retained_count_for_axis(context, first_entry.axis)
    second_retained, second_retained_policy =
        _white_lindsey_retained_count_for_axis(context, second_entry.axis)
    retained_policy =
        retained_policy == second_retained_policy ?
        retained_policy :
        :mixed_axis_retained_counts
    fixed_side = Symbol(context.fixed_side)
    fixed_index = Int(context.fixed_index)
    parent_dims = _white_lindsey_checked_context_tuple(
        context.parent_dims,
        Val(3),
        Int,
        :parent_dims,
    )
    side_first = ParentGaussletBases._nested_doside_1d(
        first_source,
        first_entry.interval,
        first_retained;
        enforce_symmetric_odd = true,
    )
    side_second = ParentGaussletBases._nested_doside_1d(
        second_source,
        second_entry.interval,
        second_retained;
        enforce_symmetric_odd = true,
    )
    face = ParentGaussletBases._nested_face_product(
        context.face_kind,
        fixed_side,
        side_first,
        side_second,
        fixed_index,
        parent_dims,
    )
    coefficient_matrix = face.coefficient_matrix
    return _white_lindsey_unit_coefficients_result(
        descriptor,
        :materialized_white_lindsey_facet_unit_coefficients,
        nothing,
        coefficient_matrix;
        coefficient_space = :parent_cartesian_sparse_adapter,
        parent_row_indices_available = true,
        support_indices = face.support_indices,
        source_support_row_count = length(face.support_indices),
        retained_column_count = size(coefficient_matrix, 2),
        source_support_row_index = nothing,
        retained_column_index = nothing,
        old_kernels_used = context.planned_old_calls,
        nonzero_values = (),
        active_axes = (first_entry.axis, second_entry.axis),
        active_axis_intervals = (first_entry.interval, second_entry.interval),
        active_axis_retained_counts = (first_retained, second_retained),
        retained_count_policy = retained_policy,
        doside_source_policy = (first_source_policy, second_source_policy),
    )
end

function _white_lindsey_edge_unit_coefficients_result(descriptor)
    context =
        white_lindsey_boundary_stratum_unit_coefficient_context(descriptor)
    context.status === :ready_white_lindsey_edge_kernel_context_not_materialized ||
        return _white_lindsey_unit_coefficients_result(
            descriptor,
            :blocked_white_lindsey_boundary_stratum_unit_coefficients,
            :incomplete_white_lindsey_edge_facet_kernel_input_context,
            nothing,
        )

    fixed_sides = _white_lindsey_checked_context_tuple(
        context.fixed_sides,
        Val(2),
        Symbol,
        :fixed_sides,
    )
    fixed_indices = _white_lindsey_checked_context_tuple(
        context.fixed_indices,
        Val(2),
        Int,
        :fixed_indices,
    )
    parent_dims = _white_lindsey_checked_context_tuple(
        context.parent_dims,
        Val(3),
        Int,
        :parent_dims,
    )
    doside_source, doside_source_policy =
        _white_lindsey_doside_source_for_axis(
            context.doside_source_1d,
            context.free_axis,
        )
    side = ParentGaussletBases._nested_doside_1d(
        doside_source,
        context.free_axis_interval,
        context.retained_count;
        enforce_symmetric_odd = true,
    )
    edge = ParentGaussletBases._nested_edge_product(
        context.free_axis,
        fixed_sides,
        side,
        fixed_indices,
        parent_dims,
    )
    coefficient_matrix = edge.coefficient_matrix
    return _white_lindsey_unit_coefficients_result(
        descriptor,
        :materialized_white_lindsey_edge_unit_coefficients,
        nothing,
        coefficient_matrix;
        coefficient_space = :parent_cartesian_sparse_adapter,
        parent_row_indices_available = true,
        support_indices = edge.support_indices,
        source_support_row_count = length(edge.support_indices),
        retained_column_count = size(coefficient_matrix, 2),
        source_support_row_index = nothing,
        retained_column_index = nothing,
        old_kernels_used = context.planned_old_calls,
        nonzero_values = (),
        active_axes = (context.free_axis,),
        active_axis_intervals = (context.free_axis_interval,),
        active_axis_retained_counts = (context.retained_count,),
        retained_count_policy = :scalar_retained_count,
        doside_source_policy = doside_source_policy,
    )
end

function _white_lindsey_unit_coefficients_result(
    descriptor,
    status::Symbol,
    blocker,
    coefficient_matrix;
    coefficient_space = :source_cpb_support_local,
    parent_row_indices_available = false,
    support_indices = nothing,
    source_support_row_count = nothing,
    retained_column_count = nothing,
    source_support_row_index = nothing,
    retained_column_index = nothing,
    old_kernels_used = (),
    nonzero_values = nothing,
    active_axes = nothing,
    active_axis_intervals = nothing,
    active_axis_retained_counts = nothing,
    retained_count_policy = nothing,
    doside_source_policy = nothing,
)
    materialized = !isnothing(coefficient_matrix)
    source_support_row_count_resolved =
        isnothing(source_support_row_count) && materialized ? 1 :
        source_support_row_count
    retained_column_count_resolved =
        isnothing(retained_column_count) && materialized ? 1 :
        retained_column_count
    source_support_row_index_resolved =
        isnothing(source_support_row_index) && materialized &&
        source_support_row_count_resolved == 1 ? 1 : source_support_row_index
    retained_column_index_resolved =
        isnothing(retained_column_index) && materialized &&
        retained_column_count_resolved == 1 ? 1 : retained_column_index
    nonzero_count = materialized ? _white_lindsey_nonzero_count(
        coefficient_matrix,
    ) : 0
    nonzero_values_resolved = isnothing(nonzero_values) ?
                              (materialized ? (1.0,) : ()) :
                              nonzero_values
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
        source_cpb_intervals =
            _white_lindsey_descriptor_property(descriptor, :source_cpb_intervals),
        source_axis_intervals =
            _white_lindsey_descriptor_property(descriptor, :source_axis_intervals),
        active_product_axis_intervals =
            _white_lindsey_descriptor_property(
                descriptor,
                :active_product_axis_intervals,
            ),
        free_axis = _white_lindsey_descriptor_property(descriptor, :free_axis),
        free_axis_interval =
            _white_lindsey_descriptor_property(descriptor, :free_axis_interval),
        fixed_axis_coordinates =
            _white_lindsey_descriptor_property(
                descriptor,
                :fixed_axis_coordinates,
            ),
        fixed_side_metadata =
            _white_lindsey_descriptor_property(descriptor, :fixed_side_metadata),
        planned_old_kernel =
            _white_lindsey_descriptor_property(descriptor, :planned_old_kernel),
        coefficient_input_requirements =
            _white_lindsey_unit_coefficient_input_requirements(descriptor),
        missing_coefficient_inputs =
            _white_lindsey_missing_unit_coefficient_inputs(descriptor),
        coefficient_space,
        parent_row_indices_available,
        support_indices,
        source_support_row_count = source_support_row_count_resolved,
        retained_column_count = retained_column_count_resolved,
        source_support_row_index = source_support_row_index_resolved,
        retained_column_index = retained_column_index_resolved,
        coefficient_matrix,
        nonzero_count,
        nonzero_values = nonzero_values_resolved,
        old_kernels_used,
        active_axes,
        active_axis_intervals,
        active_axis_retained_counts,
        retained_count_policy,
        doside_source_policy,
        coefficient_maps_materialized = materialized,
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

function _white_lindsey_checked_context_tuple(
    value,
    ::Val{N},
    ::Type{T},
    field::Symbol,
) where {N,T}
    value isa Tuple && length(value) == N || throw(
        ArgumentError("White--Lindsey coefficient context requires $(field) length $(N)"),
    )
    return ntuple(index -> T(value[index]), N)
end

function _white_lindsey_facet_active_axis_entries(context)
    active_entries = context.active_product_axis_intervals
    active_entries isa Tuple && length(active_entries) == 2 || throw(
        ArgumentError("White--Lindsey facet context requires two active axis intervals"),
    )
    return active_entries
end

function _white_lindsey_doside_source_for_axis(source, axis::Symbol)
    isnothing(source) && return nothing, :missing_doside_source_1d
    hasproperty(source, axis) &&
        return getproperty(source, axis), :axis_keyed_doside_source_1d
    source isa AbstractDict && haskey(source, axis) &&
        return source[axis], :axis_keyed_doside_source_1d
    return source, :shared_doside_source_1d
end

function _white_lindsey_retained_count_for_axis(context, axis::Symbol)
    retained_counts =
        _white_lindsey_descriptor_property(context, :retained_counts)
    if !isnothing(retained_counts)
        hasproperty(retained_counts, axis) &&
            return Int(getproperty(retained_counts, axis)),
                   :axis_keyed_retained_counts
        retained_counts isa AbstractDict && haskey(retained_counts, axis) &&
            return Int(retained_counts[axis]), :axis_keyed_retained_counts
    end
    return Int(context.retained_count),
           :scalar_retained_count_reused_for_active_axes
end

function _white_lindsey_doside_source_materializable(source)
    isnothing(source) && return false
    return hasmethod(
        ParentGaussletBases._nested_doside_1d,
        Tuple{typeof(source),UnitRange{Int},Int},
    )
end

function _white_lindsey_edge_context_source_materializable(context)
    doside_source, _ =
        _white_lindsey_doside_source_for_axis(
            context.doside_source_1d,
            context.free_axis,
        )
    return _white_lindsey_doside_source_materializable(doside_source)
end

function _white_lindsey_facet_context_sources_materializable(context)
    active_entries = _white_lindsey_facet_active_axis_entries(context)
    return all(active_entries) do entry
        doside_source, _ =
            _white_lindsey_doside_source_for_axis(
                context.doside_source_1d,
                entry.axis,
            )
        _white_lindsey_doside_source_materializable(doside_source)
    end
end

function _white_lindsey_nonzero_count(coefficient_matrix)
    SparseArrays.issparse(coefficient_matrix) &&
        return SparseArrays.nnz(coefficient_matrix)
    return count(!iszero, coefficient_matrix)
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
        source_cpb_intervals =
            _white_lindsey_descriptor_property(input, :source_cpb_intervals),
        source_axis_intervals =
            _white_lindsey_descriptor_property(input, :source_axis_intervals),
        active_product_axis_intervals =
            _white_lindsey_descriptor_property(
                input,
                :active_product_axis_intervals,
            ),
        free_axis = _white_lindsey_descriptor_property(input, :free_axis),
        free_axis_interval =
            _white_lindsey_descriptor_property(input, :free_axis_interval),
        fixed_axis_coordinates =
            _white_lindsey_descriptor_property(input, :fixed_axis_coordinates),
        fixed_side_metadata =
            _white_lindsey_descriptor_property(input, :fixed_side_metadata),
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

function _white_lindsey_unit_coefficient_context_status(descriptor)
    descriptor_status = _white_lindsey_descriptor_property(descriptor, :status)
    descriptor_status ===
        :available_metadata_only_white_lindsey_unit_adapter_descriptor || return (
        :blocked_unavailable_white_lindsey_unit_coefficient_context,
        _white_lindsey_descriptor_property(
            descriptor,
            :blocker,
            :unavailable_white_lindsey_unit_adapter_descriptor,
        ),
        (),
    )

    stratum_kind = _white_lindsey_descriptor_property(descriptor, :stratum_kind)
    missing = _white_lindsey_context_missing_inputs(descriptor)
    if stratum_kind in (:facet_cpb, :face_cpb)
        isempty(missing) && return (
            :ready_white_lindsey_facet_kernel_context_not_materialized,
            nothing,
            (:_nested_doside_1d, :_nested_face_product),
        )
        return (
            :blocked_missing_white_lindsey_facet_kernel_context,
            first(missing),
            (:_nested_doside_1d, :_nested_face_product),
        )
    elseif stratum_kind === :edge_cpb
        isempty(missing) && return (
            :ready_white_lindsey_edge_kernel_context_not_materialized,
            nothing,
            (:_nested_doside_1d, :_nested_edge_product),
        )
        return (
            :blocked_missing_white_lindsey_edge_kernel_context,
            first(missing),
            (:_nested_doside_1d, :_nested_edge_product),
        )
    elseif stratum_kind === :corner_cpb
        isempty(missing) && return (
            :ready_white_lindsey_corner_kernel_context_not_materialized,
            nothing,
            (:_nested_corner_piece,),
        )
        return (
            :blocked_missing_white_lindsey_corner_kernel_context,
            first(missing),
            (:_nested_corner_piece,),
        )
    end
    return (
        :blocked_unknown_white_lindsey_unit_coefficient_context,
        :unknown_white_lindsey_unit_coefficients_stratum,
        (),
    )
end

function _white_lindsey_unit_coefficient_context(
    descriptor,
    status::Symbol,
    blocker,
    planned_old_calls,
)
    stratum_kind = _white_lindsey_descriptor_property(descriptor, :stratum_kind)
    return (;
        object_kind = :white_lindsey_boundary_stratum_unit_coefficient_context,
        status,
        blocker,
        unit_key = _white_lindsey_descriptor_property(descriptor, :unit_key),
        unit_index = _white_lindsey_descriptor_property(descriptor, :unit_index),
        unit_kind = _white_lindsey_descriptor_property(descriptor, :unit_kind),
        stratum_kind,
        planned_old_calls,
        missing_inputs = _white_lindsey_context_missing_inputs(descriptor),
        source_cpb_role =
            _white_lindsey_descriptor_property(descriptor, :source_cpb_role),
        source_cpb_intervals =
            _white_lindsey_descriptor_property(descriptor, :source_cpb_intervals),
        source_axis_intervals =
            _white_lindsey_descriptor_property(descriptor, :source_axis_intervals),
        active_product_axis_intervals =
            _white_lindsey_descriptor_property(
                descriptor,
                :active_product_axis_intervals,
            ),
        face_kind = _white_lindsey_context_face_kind(descriptor),
        fixed_side = _white_lindsey_context_fixed_side(descriptor),
        fixed_index = _white_lindsey_context_fixed_index(descriptor),
        free_axis = _white_lindsey_descriptor_property(descriptor, :free_axis),
        free_axis_interval =
            _white_lindsey_descriptor_property(descriptor, :free_axis_interval),
        fixed_sides = _white_lindsey_context_fixed_sides(descriptor),
        fixed_indices = _white_lindsey_context_fixed_indices(descriptor),
        retained_count =
            _white_lindsey_descriptor_property(descriptor, :retained_count),
        retained_counts =
            _white_lindsey_descriptor_property(descriptor, :retained_counts),
        parent_dims = _white_lindsey_descriptor_property(descriptor, :parent_dims),
        doside_source_1d =
            _white_lindsey_descriptor_property(descriptor, :doside_source_1d),
        coefficient_maps_materialized = false,
        edge_facet_coefficient_maps_materialized = false,
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

function _white_lindsey_context_missing_inputs(descriptor)
    stratum_kind = _white_lindsey_descriptor_property(descriptor, :stratum_kind)
    if stratum_kind in (:facet_cpb, :face_cpb, :edge_cpb)
        return _white_lindsey_missing_unit_coefficient_inputs(descriptor)
    elseif stratum_kind === :corner_cpb
        missing = Symbol[]
        isnothing(_white_lindsey_context_fixed_sides(descriptor)) &&
            push!(missing, :missing_white_lindsey_fixed_side_metadata)
        isnothing(_white_lindsey_context_fixed_indices(descriptor)) &&
            push!(missing, :missing_white_lindsey_corner_fixed_coordinates)
        return Tuple(missing)
    end
    return ()
end

function _white_lindsey_context_fixed_side_metadata(descriptor)
    return _white_lindsey_descriptor_property(descriptor, :fixed_side_metadata)
end

function _white_lindsey_context_fixed_sides(descriptor)
    fixed_side_metadata = _white_lindsey_context_fixed_side_metadata(descriptor)
    isnothing(fixed_side_metadata) && return nothing
    return Tuple(entry.side for entry in fixed_side_metadata)
end

function _white_lindsey_context_fixed_indices(descriptor)
    fixed_coordinates =
        _white_lindsey_descriptor_property(descriptor, :fixed_axis_coordinates)
    isnothing(fixed_coordinates) && return nothing
    return Tuple(entry.coordinate for entry in fixed_coordinates)
end

function _white_lindsey_context_fixed_side(descriptor)
    fixed_sides = _white_lindsey_context_fixed_sides(descriptor)
    isnothing(fixed_sides) && return nothing
    length(fixed_sides) == 1 || return nothing
    return only(fixed_sides)
end

function _white_lindsey_context_fixed_index(descriptor)
    fixed_indices = _white_lindsey_context_fixed_indices(descriptor)
    isnothing(fixed_indices) && return nothing
    length(fixed_indices) == 1 || return nothing
    return only(fixed_indices)
end

function _white_lindsey_context_face_kind(descriptor)
    active_axes = Tuple(
        entry.axis for entry in _white_lindsey_descriptor_property(
            descriptor,
            :active_product_axis_intervals,
            (),
        )
    )
    length(active_axes) == 2 || return nothing
    active_axes == (:x, :y) && return :xy
    active_axes == (:x, :z) && return :xz
    active_axes == (:y, :z) && return :yz
    return nothing
end

function _white_lindsey_unit_coefficient_input_requirements(descriptor)
    stratum_kind = _white_lindsey_descriptor_property(descriptor, :stratum_kind)
    if stratum_kind in (:facet_cpb, :face_cpb)
        missing = _white_lindsey_missing_unit_coefficient_inputs(descriptor)
        return (;
            object_kind = :white_lindsey_unit_coefficient_input_requirements,
            status = isempty(missing) ?
                     :available_white_lindsey_facet_kernel_context_inputs :
                     :blocked_missing_white_lindsey_facet_kernel_inputs,
            required_old_kernel = :_nested_face_product,
            required_1d_helper = :_nested_doside_1d,
            required_inputs = (
                :doside_source_1d,
                :active_product_axis_intervals,
                :fixed_side_metadata,
                :retained_count,
                :parent_dims,
            ),
            missing_inputs = missing,
        )
    elseif stratum_kind === :edge_cpb
        missing = _white_lindsey_missing_unit_coefficient_inputs(descriptor)
        return (;
            object_kind = :white_lindsey_unit_coefficient_input_requirements,
            status = isempty(missing) ?
                     :available_white_lindsey_edge_kernel_context_inputs :
                     :blocked_missing_white_lindsey_edge_kernel_inputs,
            required_old_kernel = :_nested_edge_product,
            required_1d_helper = :_nested_doside_1d,
            required_inputs = (
                :doside_source_1d,
                :free_axis_interval,
                :fixed_side_metadata,
                :retained_count,
                :parent_dims,
            ),
            missing_inputs = missing,
        )
    end
    stratum_kind === :corner_cpb && return (;
        object_kind = :white_lindsey_unit_coefficient_input_requirements,
        status = :available_corner_support_local_coefficients,
        required_old_kernel = :_nested_corner_piece,
        required_1d_helper = nothing,
        required_inputs = (:fixed_axis_coordinates,),
        missing_inputs =
            _white_lindsey_missing_unit_coefficient_inputs(descriptor),
    )
    return (;
        object_kind = :white_lindsey_unit_coefficient_input_requirements,
        status = :blocked_unknown_white_lindsey_unit_coefficient_inputs,
        required_old_kernel = nothing,
        required_1d_helper = nothing,
        required_inputs = (),
        missing_inputs =
            _white_lindsey_missing_unit_coefficient_inputs(descriptor),
    )
end

function _white_lindsey_missing_unit_coefficient_inputs(descriptor)
    stratum_kind = _white_lindsey_descriptor_property(descriptor, :stratum_kind)
    missing = Symbol[]
    if stratum_kind in (:facet_cpb, :face_cpb, :edge_cpb)
        isnothing(_white_lindsey_descriptor_property(descriptor, :doside_source_1d)) &&
            push!(missing, :missing_white_lindsey_doside_source_1d)
        isnothing(_white_lindsey_descriptor_property(descriptor, :retained_count)) &&
            push!(missing, :missing_white_lindsey_retained_count)
        isnothing(_white_lindsey_descriptor_property(descriptor, :parent_dims)) &&
            push!(missing, :missing_white_lindsey_parent_dims)
        isnothing(
            _white_lindsey_descriptor_property(descriptor, :fixed_side_metadata),
        ) && push!(missing, :missing_white_lindsey_fixed_side_metadata)
        if stratum_kind === :edge_cpb
            isnothing(
                _white_lindsey_descriptor_property(
                    descriptor,
                    :free_axis_interval,
                ),
            ) && push!(missing, :missing_white_lindsey_free_axis_interval)
        else
            isempty(
                _white_lindsey_descriptor_property(
                    descriptor,
                    :active_product_axis_intervals,
                    (),
                ),
            ) && push!(
                missing,
                :missing_white_lindsey_active_product_axis_intervals,
            )
        end
    elseif stratum_kind === :corner_cpb
        fixed_coordinates =
            _white_lindsey_descriptor_property(descriptor, :fixed_axis_coordinates)
        (!isnothing(fixed_coordinates) && length(fixed_coordinates) == 3) ||
            push!(missing, :missing_white_lindsey_corner_fixed_coordinates)
    end
    return Tuple(missing)
end
