# Overlap-only White--Lindsey pair-block pilot.

"""
    white_lindsey_boundary_stratum_overlap_block(pair_unit_coefficients; parent_axis_counts, overlap_1d)
    white_lindsey_boundary_stratum_overlap_block(unit_pair; parent_axis_counts, overlap_1d)

Materialize one White--Lindsey boundary-stratum overlap block from prepared
left/right unit coefficient maps. This contracts only over the left/right
support rows and does not build dense parent-parent matrices, one-body terms
beyond overlap, Hamiltonian data, exports, artifacts, IDA/MWG data, or Coulomb.
"""
function white_lindsey_boundary_stratum_overlap_block(
    pair_unit_coefficients;
    parent_axis_counts,
    overlap_1d,
)
    axis_counts = _axis_counts_tuple(parent_axis_counts)
    overlap_axes = _overlap_1d_tuple(overlap_1d)
    _assert_overlap_axis_sizes(overlap_axes, axis_counts)
    return _white_lindsey_boundary_stratum_product_block(
        pair_unit_coefficients,
        :overlap,
        axis_counts,
        overlap_axes,
        :white_lindsey_boundary_stratum_overlap_adapter,
        (;),
    )
end

function white_lindsey_boundary_stratum_overlap_block(
    unit_pair::CUP.UnitPairRecord;
    parent_axis_counts,
    overlap_1d,
)
    return white_lindsey_boundary_stratum_overlap_block(
        white_lindsey_boundary_stratum_pair_unit_coefficients(unit_pair);
        parent_axis_counts,
        overlap_1d,
    )
end

function _white_lindsey_boundary_stratum_product_block(
    pair_unit_coefficients,
    term::Symbol,
    axis_counts::NTuple{3,Int},
    operator_axes,
    materialization_path::Symbol,
    metadata,
)
    _assert_white_lindsey_pair_unit_coefficients_ready(pair_unit_coefficients)

    left_support = pair_unit_coefficients.left_support_indices
    right_support = pair_unit_coefficients.right_support_indices
    _assert_white_lindsey_overlap_support(left_support, :left)
    _assert_white_lindsey_overlap_support(right_support, :right)

    left_states =
        Tuple(ParentGaussletBases._cartesian_unflat_index(index, axis_counts)
              for index in left_support)
    right_states =
        Tuple(ParentGaussletBases._cartesian_unflat_index(index, axis_counts)
              for index in right_support)
    support_overlap =
        Matrix{Float64}(undef, length(left_states), length(right_states))
    _fill_source_mode_product_block!(
        support_overlap,
        left_states,
        right_states,
        operator_axes[1],
        operator_axes[2],
        operator_axes[3],
    )

    left_support_coefficients = Matrix{Float64}(
        pair_unit_coefficients.left_coefficient_matrix[left_support, :],
    )
    right_support_coefficients = Matrix{Float64}(
        pair_unit_coefficients.right_coefficient_matrix[right_support, :],
    )
    block = Matrix{Float64}(
        transpose(left_support_coefficients) *
        support_overlap *
        right_support_coefficients,
    )

    return PairBlockMaterializationResult(
        term,
        pair_unit_coefficients.pair_key,
        block,
        true,
        true,
        true,
        false,
        false,
        false,
        merge(
            (;
                materialization_path,
                pair_family = pair_unit_coefficients.pair_family,
                left_stratum_kind = pair_unit_coefficients.left_stratum_kind,
                right_stratum_kind = pair_unit_coefficients.right_stratum_kind,
                left_support_count = length(left_support),
                right_support_count = length(right_support),
                left_retained_column_count =
                    pair_unit_coefficients.left_retained_column_count,
                right_retained_column_count =
                    pair_unit_coefficients.right_retained_column_count,
                parent_axis_counts = axis_counts,
                support_overlap_shape = size(support_overlap),
                left_support_coefficient_shape =
                    size(left_support_coefficients),
                right_support_coefficient_shape =
                    size(right_support_coefficients),
                coefficient_source =
                    :white_lindsey_boundary_stratum_pair_unit_coefficients,
                local_pair_block_materialized = true,
                source_operator_blocks_materialized = true,
                final_pair_blocks_materialized = true,
                operator_blocks_materialized = false,
                hamiltonian_data_materialized = false,
                artifacts_materialized = false,
                ida_mwg_data_materialized = false,
                dense_parent_parent_overlap_materialized = false,
            ),
            NamedTuple(metadata),
        ),
    )
end

function _assert_white_lindsey_pair_unit_coefficients_ready(pair_unit_coefficients)
    _white_lindsey_descriptor_property(pair_unit_coefficients, :object_kind) ===
        :white_lindsey_boundary_stratum_pair_unit_coefficients || throw(
        ArgumentError("White--Lindsey overlap requires prepared pair unit coefficients"),
    )
    pair_unit_coefficients.status ===
        :materialized_white_lindsey_pair_unit_coefficients || throw(
        ArgumentError("White--Lindsey overlap requires materialized left/right unit coefficients"),
    )
    pair_unit_coefficients.pair_unit_coefficient_maps_materialized || throw(
        ArgumentError("White--Lindsey overlap requires materialized pair unit coefficient maps"),
    )
    return nothing
end

function _assert_white_lindsey_overlap_support(support_indices, side::Symbol)
    isnothing(support_indices) && throw(
        ArgumentError("White--Lindsey overlap requires $(side) parent support indices"),
    )
    isempty(support_indices) && throw(
        ArgumentError("White--Lindsey overlap requires nonempty $(side) support indices"),
    )
    return nothing
end
