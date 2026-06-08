# Kinetic-only White--Lindsey pair-block pilot.

"""
    white_lindsey_boundary_stratum_kinetic_block(pair_unit_coefficients; parent_axis_counts, overlap_1d, kinetic_1d)
    white_lindsey_boundary_stratum_kinetic_block(unit_pair; parent_axis_counts, overlap_1d, kinetic_1d)

Materialize one White--Lindsey boundary-stratum kinetic block from prepared
left/right unit coefficient maps using the factorized one-body form
`(K,S,S) + (S,K,S) + (S,S,K)`. Caller-supplied 1D kinetic factors own signs
and prefactors. This does not build Hamiltonian data, exports, artifacts,
IDA/MWG data, or Coulomb.
"""
function white_lindsey_boundary_stratum_kinetic_block(
    pair_unit_coefficients;
    parent_axis_counts,
    overlap_1d,
    kinetic_1d,
)
    axis_counts = _axis_counts_tuple(parent_axis_counts)
    overlap_axes = _overlap_1d_tuple(overlap_1d)
    kinetic_axes = _operator_1d_tuple(kinetic_1d, "kinetic_1d")
    _assert_overlap_axis_sizes(overlap_axes, axis_counts)
    _assert_operator_axis_sizes(kinetic_axes, axis_counts, "kinetic_1d")

    component_x = _white_lindsey_boundary_stratum_product_block(
        pair_unit_coefficients,
        :kinetic_x_component,
        axis_counts,
        (kinetic_axes[1], overlap_axes[2], overlap_axes[3]),
        :white_lindsey_boundary_stratum_kinetic_component_adapter,
        (; kinetic_component_axis = :x),
    )
    component_y = _white_lindsey_boundary_stratum_product_block(
        pair_unit_coefficients,
        :kinetic_y_component,
        axis_counts,
        (overlap_axes[1], kinetic_axes[2], overlap_axes[3]),
        :white_lindsey_boundary_stratum_kinetic_component_adapter,
        (; kinetic_component_axis = :y),
    )
    component_z = _white_lindsey_boundary_stratum_product_block(
        pair_unit_coefficients,
        :kinetic_z_component,
        axis_counts,
        (overlap_axes[1], overlap_axes[2], kinetic_axes[3]),
        :white_lindsey_boundary_stratum_kinetic_component_adapter,
        (; kinetic_component_axis = :z),
    )
    block = component_x.block + component_y.block + component_z.block

    return PairBlockMaterializationResult(
        :kinetic,
        pair_unit_coefficients.pair_key,
        block,
        true,
        true,
        true,
        false,
        false,
        false,
        merge(
            component_x.metadata,
            (;
                materialization_path =
                    :white_lindsey_boundary_stratum_kinetic_adapter,
                kinetic_factor_form =
                    :factorized_cartesian_sum_kss_sks_ssk,
                kinetic_component_axes = (:x, :y, :z),
                kinetic_component_terms = (
                    component_x.term,
                    component_y.term,
                    component_z.term,
                ),
                kinetic_component_block_shapes = (
                    size(component_x.block),
                    size(component_y.block),
                    size(component_z.block),
                ),
                operator_blocks_materialized = false,
                hamiltonian_data_materialized = false,
                artifacts_materialized = false,
                ida_mwg_data_materialized = false,
                dense_parent_parent_overlap_materialized = false,
            ),
        ),
    )
end

function white_lindsey_boundary_stratum_kinetic_block(
    unit_pair::CUP.UnitPairRecord;
    parent_axis_counts,
    overlap_1d,
    kinetic_1d,
)
    return white_lindsey_boundary_stratum_kinetic_block(
        white_lindsey_boundary_stratum_pair_unit_coefficients(unit_pair);
        parent_axis_counts,
        overlap_1d,
        kinetic_1d,
    )
end
