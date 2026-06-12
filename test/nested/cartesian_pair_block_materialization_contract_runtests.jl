using Test
using GaussletBases
using LinearAlgebra: I

const CPBM = GaussletBases.CartesianPairBlockMaterialization
const CPOPForPairBlocks = GaussletBases.CartesianPairOperatorPlans
const CUPForPairBlocks = GaussletBases.CartesianUnitPairs
const CRTCForPairBlocks = GaussletBases.CartesianRetainedUnitTransformContracts
const CRUForPairBlocks = GaussletBases.CartesianRetainedUnits
const CTLForPairBlocks = GaussletBases.CartesianTerminalLowering
const CRCForPairBlocks = GaussletBases.CartesianRouteCore
const CPBForPairBlocks = GaussletBases.CartesianCPB
const CRPSForPairBlocks = GaussletBases.CartesianRawProductSources

function _pair_block_count(counts, field::Symbol, value)
    matches = Tuple(entry for entry in counts if getproperty(entry, field) == value)
    isempty(matches) && return 0
    return only(matches).count
end

function _pair_block_minimal_lowering_plan()
    return CTLForPairBlocks.TerminalLoweringPlan(
        CTLForPairBlocks.PQSLowering(q = 3),
        (),
        (),
        (;
            object_kind = :synthetic_terminal_lowering_plan_summary,
            status = :available_terminal_lowering_plan,
            policy_kind = :pair_block_materialization_synthetic_lowering,
            terminal_region_count = 0,
            available_contract_count = 0,
            selected_contract_count = 0,
            selected_contract_kinds = (),
            all_terminal_regions_have_selected_contract = true,
            materialized = false,
            retained_spaces_materialized = false,
            final_retained_units_materialized = false,
            pair_inventory_materialized = false,
            operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
        (; fixture = :cartesian_pair_block_materialization_contract),
    )
end

function _pair_block_retained_unit(
    unit_key::Symbol,
    unit_index::Int,
    unit_kind::Symbol,
    lowering_kind::Symbol,
    retained_rule::Symbol,
    realization_rule;
    owned_support = nothing,
    source_cpbs = (),
    source_cpb_index = nothing,
    metadata = (; route_core_sidecar_status = :not_materialized),
)
    return CRUForPairBlocks.RetainedUnitRecord(
        unit_key,
        unit_index,
        unit_kind,
        Symbol(unit_key, "_contract"),
        Symbol(unit_key, "_terminal_region"),
        :synthetic_terminal_region,
        :synthetic_terminal_region,
        lowering_kind,
        retained_rule,
        realization_rule,
        owned_support,
        Tuple(source_cpbs),
        source_cpb_index,
        :not_materialized,
        nothing,
        :not_materialized,
        nothing,
        nothing,
        false,
        NamedTuple(metadata),
    )
end

function _pair_block_retained_plan()
    direct_source = CPBForPairBlocks.filled_cpb(
        1:2,
        2:3,
        1:2;
        role = :pair_block_direct_source_cpb,
    )
    direct = _pair_block_retained_unit(
        :pair_block_direct_unit,
        1,
        :direct_cpb_retained_unit,
        :direct_core_identity_cpb,
        :direct_source_modes,
        :direct_or_trivial_embedding;
        owned_support = CRCForPairBlocks.owned_cpb(direct_source),
        source_cpbs = (direct_source,),
    )
    pqs = _pair_block_retained_unit(
        :pair_block_pqs_unit,
        2,
        :pqs_shell_retained_unit,
        :pqs_filled_source_cpb,
        :pqs_boundary_comx_product_modes,
        :shell_projection_lowdin,
    )
    units = (direct, pqs)
    return CRUForPairBlocks.RetainedUnitPlan(
        CRUForPairBlocks.MetadataOnlyRetainedUnits(),
        _pair_block_minimal_lowering_plan(),
        units,
        (;
            object_kind = :synthetic_retained_unit_plan_summary,
            status = :available_retained_unit_plan,
            retained_unit_count = length(units),
            materialized = false,
            transforms_materialized = false,
            coefficient_maps_materialized = false,
            pair_inventory_materialized = false,
            operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
        (; fixture = :cartesian_pair_block_materialization_contract),
    )
end

function _pair_block_rectangular_retained_plan()
    left_source = CPBForPairBlocks.cpb(
        1:2,
        1:2,
        1:1;
        role = :pair_block_direct_left_source_cpb,
    )
    right_source = CPBForPairBlocks.cpb(
        3:3,
        2:3,
        1:3;
        role = :pair_block_direct_right_source_cpb,
    )
    left_direct = _pair_block_retained_unit(
        :pair_block_direct_left_unit,
        1,
        :direct_cpb_retained_unit,
        :direct_core_identity_cpb,
        :direct_source_modes,
        :direct_or_trivial_embedding;
        owned_support = CRCForPairBlocks.owned_cpb(left_source),
        source_cpbs = (left_source,),
    )
    right_direct = _pair_block_retained_unit(
        :pair_block_direct_right_unit,
        2,
        :direct_cpb_retained_unit,
        :direct_core_identity_cpb,
        :direct_source_modes,
        :direct_or_trivial_embedding;
        owned_support = CRCForPairBlocks.owned_cpb(right_source),
        source_cpbs = (right_source,),
    )
    pqs = _pair_block_retained_unit(
        :pair_block_rectangular_pqs_unit,
        3,
        :pqs_shell_retained_unit,
        :pqs_filled_source_cpb,
        :pqs_boundary_comx_product_modes,
        :shell_projection_lowdin,
    )
    units = (left_direct, right_direct, pqs)
    return CRUForPairBlocks.RetainedUnitPlan(
        CRUForPairBlocks.MetadataOnlyRetainedUnits(),
        _pair_block_minimal_lowering_plan(),
        units,
        (;
            object_kind = :synthetic_retained_unit_plan_summary,
            status = :available_retained_unit_plan,
            retained_unit_count = length(units),
            materialized = false,
            transforms_materialized = false,
            coefficient_maps_materialized = false,
            pair_inventory_materialized = false,
            operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
        (; fixture = :cartesian_pair_block_rectangular_direct_direct),
    )
end

function _pair_block_pqs_retained_plan()
    left_source = CPBForPairBlocks.filled_cpb(
        1:3,
        1:3,
        1:3;
        role = :pair_block_pqs_left_source_cpb,
    )
    right_source = CPBForPairBlocks.filled_cpb(
        4:6,
        1:3,
        1:3;
        role = :pair_block_pqs_right_source_cpb,
    )
    left_pqs = _pair_block_retained_unit(
        :pair_block_pqs_left_unit,
        1,
        :pqs_shell_retained_unit,
        :pqs_filled_source_cpb,
        :pqs_boundary_comx_product_modes,
        :shell_projection_lowdin;
        source_cpbs = (left_source,),
        metadata = (;
            q = 3,
            source_mode_shape = nothing,
            raw_product_source_axis_transform_matrices = (;
                x = [
                    1.0 0.0 0.0
                    0.0 1.0 0.0
                    0.0 0.0 1.0
                ],
                y = [
                    1.0 0.0 0.0
                    0.0 1.0 0.0
                    0.0 0.0 1.0
                ],
                z = [
                    1.0 0.0 0.0
                    0.0 1.0 0.0
                    0.0 0.0 1.0
                ],
            ),
        ),
    )
    right_pqs = _pair_block_retained_unit(
        :pair_block_pqs_right_unit,
        2,
        :pqs_shell_retained_unit,
        :pqs_filled_source_cpb,
        :pqs_boundary_comx_product_modes,
        :shell_projection_lowdin;
        source_cpbs = (right_source,),
        metadata = (;
            q = 9,
            source_mode_shape = (5, 4, 3),
            raw_product_source_axis_transform_matrices = (;
                x = ones(3, 5),
                y = ones(3, 4),
                z = ones(3, 3),
            ),
        ),
    )
    units = (left_pqs, right_pqs)
    return CRUForPairBlocks.RetainedUnitPlan(
        CRUForPairBlocks.MetadataOnlyRetainedUnits(),
        _pair_block_minimal_lowering_plan(),
        units,
        (;
            object_kind = :synthetic_retained_unit_plan_summary,
            status = :available_retained_unit_plan,
            retained_unit_count = length(units),
            materialized = false,
            transforms_materialized = false,
            coefficient_maps_materialized = false,
            pair_inventory_materialized = false,
            operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
        (; fixture = :cartesian_pair_block_pqs_source_pair_preflight),
    )
end

function _pair_block_mixed_pqs_retained_plan()
    direct_source = CPBForPairBlocks.filled_cpb(
        1:2,
        1:2,
        1:2;
        role = :pair_block_selector_direct_source_cpb,
    )
    pqs_source = CPBForPairBlocks.filled_cpb(
        1:3,
        1:3,
        1:3;
        role = :pair_block_selector_pqs_source_cpb,
    )
    direct = _pair_block_retained_unit(
        :pair_block_selector_direct_unit,
        1,
        :direct_cpb_retained_unit,
        :direct_core_identity_cpb,
        :direct_source_modes,
        :direct_or_trivial_embedding;
        owned_support = CRCForPairBlocks.owned_cpb(direct_source),
        source_cpbs = (direct_source,),
    )
    pqs = _pair_block_retained_unit(
        :pair_block_selector_pqs_unit,
        2,
        :pqs_shell_retained_unit,
        :pqs_filled_source_cpb,
        :pqs_boundary_comx_product_modes,
        :shell_projection_lowdin;
        source_cpbs = (pqs_source,),
        metadata = (; q = 3, source_mode_shape = nothing),
    )
    units = (direct, pqs)
    return CRUForPairBlocks.RetainedUnitPlan(
        CRUForPairBlocks.MetadataOnlyRetainedUnits(),
        _pair_block_minimal_lowering_plan(),
        units,
        (;
            object_kind = :synthetic_retained_unit_plan_summary,
            status = :available_retained_unit_plan,
            retained_unit_count = length(units),
            materialized = false,
            transforms_materialized = false,
            coefficient_maps_materialized = false,
            pair_inventory_materialized = false,
            operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
        (; fixture = :cartesian_pair_block_pqs_source_selector),
    )
end

function _pair_block_white_lindsey_retained_plan()
    left_source = CPBForPairBlocks.slab_cpb(
        1:1,
        1:3,
        1:3;
        role = :pair_block_lw_left_facet_source_cpb,
        metadata = (; stratum_kind = :facet_cpb, source_cpb_index = 1),
    )
    right_source = CPBForPairBlocks.cpb(
        4:4,
        2:2,
        1:3;
        role = :pair_block_lw_right_edge_source_cpb,
        metadata = (; stratum_kind = :edge_cpb, source_cpb_index = 2),
    )
    corner_source = CPBForPairBlocks.cpb(
        4:4,
        3:3,
        3:3;
        role = :pair_block_lw_corner_source_cpb,
        metadata = (; stratum_kind = :corner_cpb, source_cpb_index = 3),
    )
    left_lw = _pair_block_retained_unit(
        :pair_block_lw_left_unit,
        1,
        :white_lindsey_boundary_stratum_retained_unit,
        :white_lindsey_boundary_strata,
        :white_lindsey_boundary_stratum_product,
        :direct_or_trivial_embedding;
        owned_support = CRCForPairBlocks.owned_cpb(left_source),
        source_cpbs = (left_source,),
        source_cpb_index = 1,
        metadata = (; stratum_kind = :facet_cpb, source_cpb_index = 1),
    )
    right_lw = _pair_block_retained_unit(
        :pair_block_lw_right_unit,
        2,
        :white_lindsey_boundary_stratum_retained_unit,
        :white_lindsey_boundary_strata,
        :white_lindsey_boundary_stratum_product,
        :direct_or_trivial_embedding;
        owned_support = CRCForPairBlocks.owned_cpb(right_source),
        source_cpbs = (right_source,),
        source_cpb_index = 2,
        metadata = (; stratum_kind = :edge_cpb, source_cpb_index = 2),
    )
    corner_lw = _pair_block_retained_unit(
        :pair_block_lw_corner_unit,
        3,
        :white_lindsey_boundary_stratum_retained_unit,
        :white_lindsey_boundary_strata,
        :white_lindsey_boundary_stratum_product,
        :direct_or_trivial_embedding;
        owned_support = CRCForPairBlocks.owned_cpb(corner_source),
        source_cpbs = (corner_source,),
        source_cpb_index = 3,
        metadata = (; stratum_kind = :corner_cpb, source_cpb_index = 3),
    )
    units = (left_lw, right_lw, corner_lw)
    return CRUForPairBlocks.RetainedUnitPlan(
        CRUForPairBlocks.MetadataOnlyRetainedUnits(),
        _pair_block_minimal_lowering_plan(),
        units,
        (;
            object_kind = :synthetic_retained_unit_plan_summary,
            status = :available_retained_unit_plan,
            retained_unit_count = length(units),
            materialized = false,
            transforms_materialized = false,
            coefficient_maps_materialized = false,
            pair_inventory_materialized = false,
            operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
        (; fixture = :cartesian_pair_block_white_lindsey_boundary_strata),
    )
end

function _pair_block_transform_contract_with_metadata(contract, metadata)
    return CRTCForPairBlocks.RetainedUnitTransformContract(
        contract.unit_key,
        contract.unit_index,
        contract.unit_kind,
        contract.lowering_kind,
        contract.retained_rule,
        contract.realization_rule,
        contract.source_cpbs,
        contract.transform_path,
        contract.realization_path,
        contract.dimension_status,
        contract.column_range_status,
        contract.materialized,
        contract.blocker,
        NamedTuple(metadata),
    )
end

function _pair_block_transform_plan_with_contracts(transform_plan, contracts)
    return CRTCForPairBlocks.RetainedUnitTransformContractPlan(
        transform_plan.policy,
        transform_plan.retained_unit_plan,
        Tuple(contracts),
        transform_plan.summary,
        transform_plan.metadata,
    )
end

function _pair_block_record_for(plan, left_key::Symbol, right_key::Symbol)
    return only(
        record for record in CPBM.pair_block_materialization_records(plan)
        if record.pair_key == (left_key, right_key)
    )
end

function _pair_block_unit_pair_for(plan, left_key::Symbol, right_key::Symbol)
    return only(
        pair for pair in CUPForPairBlocks.unit_pairs(plan)
        if pair.pair_key == (left_key, right_key)
    )
end

function _pair_block_batch_result_for(batch_result, left_key::Symbol, right_key::Symbol)
    return only(
        result for result in batch_result.materialized_results
        if result.pair_key == (left_key, right_key)
    )
end

function _pair_block_result_with_metadata(result, metadata)
    return CPBM.PairBlockMaterializationResult(
        result.term,
        result.pair_key,
        result.block,
        result.materialized,
        result.source_operator_blocks_materialized,
        result.final_pair_blocks_materialized,
        result.operator_blocks_materialized,
        result.hamiltonian_data_materialized,
        result.artifacts_materialized,
        NamedTuple(metadata),
    )
end

function _pair_block_record_with_metadata(record, metadata)
    return CPBM.PairBlockMaterializationRecord(
        record.pair_key,
        record.pair_index,
        record.pair_family,
        record.source_operator_path,
        record.transform_path,
        record.realization_path,
        record.final_block_path,
        record.supported_terms,
        record.materialization_path,
        record.readiness_status,
        record.blocker,
        record.materialized,
        NamedTuple(metadata),
    )
end

function _pair_block_overlap_axes()
    overlap_x = [
        1.00 0.11 0.13 0.17
        0.11 1.20 0.19 0.23
        0.13 0.19 1.40 0.29
        0.17 0.23 0.29 1.60
    ]
    overlap_y = [
        1.70 0.31 0.37 0.41
        0.31 1.90 0.43 0.47
        0.37 0.43 2.10 0.53
        0.41 0.47 0.53 2.30
    ]
    overlap_z = [
        2.50 0.59 0.61 0.67
        0.59 2.70 0.71 0.73
        0.61 0.71 2.90 0.79
        0.67 0.73 0.79 3.10
    ]
    return overlap_x, overlap_y, overlap_z
end

function _pair_block_position_axes()
    position_x = [
        0.20 1.10 1.30 1.70
        1.10 0.40 1.90 2.30
        1.30 1.90 0.60 2.90
        1.70 2.30 2.90 0.80
    ]
    position_y = [
        3.10 0.21 0.23 0.27
        0.21 3.30 0.29 0.31
        0.23 0.29 3.50 0.37
        0.27 0.31 0.37 3.70
    ]
    position_z = [
        4.10 0.41 0.43 0.47
        0.41 4.30 0.49 0.51
        0.43 0.49 4.50 0.57
        0.47 0.51 0.57 4.70
    ]
    return position_x, position_y, position_z
end

function _pair_block_x2_axes()
    x2_x = [
        5.10 0.61 0.63 0.67
        0.61 5.30 0.69 0.71
        0.63 0.69 5.50 0.77
        0.67 0.71 0.77 5.70
    ]
    x2_y = [
        6.10 0.81 0.83 0.87
        0.81 6.30 0.89 0.91
        0.83 0.89 6.50 0.97
        0.87 0.91 0.97 6.70
    ]
    x2_z = [
        7.10 1.01 1.03 1.07
        1.01 7.30 1.09 1.11
        1.03 1.09 7.50 1.17
        1.07 1.11 1.17 7.70
    ]
    return x2_x, x2_y, x2_z
end

function _pair_block_kinetic_axes()
    kinetic_x = [
        -1.10 0.12 0.14 0.18
        0.12 -1.30 0.20 0.24
        0.14 0.20 -1.50 0.30
        0.18 0.24 0.30 -1.70
    ]
    kinetic_y = [
        -2.10 0.32 0.34 0.38
        0.32 -2.30 0.40 0.44
        0.34 0.40 -2.50 0.50
        0.38 0.44 0.50 -2.70
    ]
    kinetic_z = [
        -3.10 0.52 0.54 0.58
        0.52 -3.30 0.60 0.64
        0.54 0.60 -3.50 0.70
        0.58 0.64 0.70 -3.70
    ]
    return kinetic_x, kinetic_y, kinetic_z
end

function _pair_block_cpb_states(source_cpb)
    ix, iy, iz = CPBForPairBlocks.intervals(source_cpb)
    return Tuple((x, y, z) for x in ix for y in iy for z in iz)
end

function _pair_block_expected_product(left_cpb, right_cpb, operator_x, operator_y, operator_z)
    left_states = _pair_block_cpb_states(left_cpb)
    right_states = _pair_block_cpb_states(right_cpb)
    return [
        operator_x[left[1], right[1]] *
        operator_y[left[2], right[2]] *
        operator_z[left[3], right[3]]
        for left in left_states, right in right_states
    ]
end

function _pair_block_expected_overlap(left_cpb, right_cpb, overlap_x, overlap_y, overlap_z)
    return _pair_block_expected_product(left_cpb, right_cpb, overlap_x, overlap_y, overlap_z)
end

function _pair_block_pqs_source_overlap_axes(left_dims, right_dims)
    overlap_x = [
        Float64(10 * lx + rx)
        for lx in 1:left_dims[1], rx in 1:right_dims[1]
    ]
    overlap_y = [
        Float64(100 * ly + 2 * ry)
        for ly in 1:left_dims[2], ry in 1:right_dims[2]
    ]
    overlap_z = [
        Float64(1000 * lz + 3 * rz)
        for lz in 1:left_dims[3], rz in 1:right_dims[3]
    ]
    return overlap_x, overlap_y, overlap_z
end

function _pair_block_pqs_source_position_axes(left_dims, right_dims)
    position_x = [
        Float64(20 * lx + 5 * rx)
        for lx in 1:left_dims[1], rx in 1:right_dims[1]
    ]
    position_y = [
        Float64(200 * ly + 7 * ry)
        for ly in 1:left_dims[2], ry in 1:right_dims[2]
    ]
    position_z = [
        Float64(2000 * lz + 11 * rz)
        for lz in 1:left_dims[3], rz in 1:right_dims[3]
    ]
    return position_x, position_y, position_z
end

function _pair_block_pqs_source_x2_axes(left_dims, right_dims)
    x2_x = [
        Float64(30 * lx + 13 * rx)
        for lx in 1:left_dims[1], rx in 1:right_dims[1]
    ]
    x2_y = [
        Float64(300 * ly + 17 * ry)
        for ly in 1:left_dims[2], ry in 1:right_dims[2]
    ]
    x2_z = [
        Float64(3000 * lz + 19 * rz)
        for lz in 1:left_dims[3], rz in 1:right_dims[3]
    ]
    return x2_x, x2_y, x2_z
end

function _pair_block_pqs_source_kinetic_axes(left_dims, right_dims)
    kinetic_x = [
        Float64(-40 * lx + 23 * rx)
        for lx in 1:left_dims[1], rx in 1:right_dims[1]
    ]
    kinetic_y = [
        Float64(-400 * ly + 29 * ry)
        for ly in 1:left_dims[2], ry in 1:right_dims[2]
    ]
    kinetic_z = [
        Float64(-4000 * lz + 31 * rz)
        for lz in 1:left_dims[3], rz in 1:right_dims[3]
    ]
    return kinetic_x, kinetic_y, kinetic_z
end

function _pair_block_pqs_source_gaussian_factor_terms(left_dims, right_dims)
    gx = Array{Float64,3}(undef, 2, left_dims[1], right_dims[1])
    gy = Array{Float64,3}(undef, 2, left_dims[2], right_dims[2])
    gz = Array{Float64,3}(undef, 2, left_dims[3], right_dims[3])
    for term in 1:2
        for lx in 1:left_dims[1], rx in 1:right_dims[1]
            gx[term, lx, rx] = Float64(100 * term + 10 * lx + rx)
        end
        for ly in 1:left_dims[2], ry in 1:right_dims[2]
            gy[term, ly, ry] = Float64(200 * term + 20 * ly + 2 * ry)
        end
        for lz in 1:left_dims[3], rz in 1:right_dims[3]
            gz[term, lz, rz] = Float64(300 * term + 30 * lz + 3 * rz)
        end
    end
    return gx, gy, gz
end

function _pair_block_pqs_source_support_gaussian_factor_terms()
    gx = Array{Float64,3}(undef, 2, 3, 3)
    gy = Array{Float64,3}(undef, 2, 3, 3)
    gz = Array{Float64,3}(undef, 2, 3, 3)
    for term in 1:2
        for left in 1:3, right in 1:3
            gx[term, left, right] = Float64(10 * term + left + right / 10)
            gy[term, left, right] = Float64(20 * term + 2 * left + right / 7)
            gz[term, left, right] = Float64(30 * term + 3 * left + right / 5)
        end
    end
    return gx, gy, gz
end

function _pair_block_project_axis_terms(axis_terms, left_transform, right_transform)
    projected = Array{Float64,3}(
        undef,
        size(axis_terms, 1),
        size(left_transform, 2),
        size(right_transform, 2),
    )
    for term in axes(axis_terms, 1)
        projected[term, :, :] .=
            transpose(left_transform) * axis_terms[term, :, :] * right_transform
    end
    return projected
end

function _pair_block_fixed_a_nested_test_basis(
    count::Int;
    a::Float64 = 0.25,
    xmax::Float64 = 10.0,
    tail_spacing::Float64 = 10.0,
)
    endpoint = (count - 1) / 2
    s = asinh(xmax / a) / (endpoint - xmax / tail_spacing)
    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count,
            mapping = AsinhMapping(; a, s, tail_spacing),
            reference_spacing = 1.0,
        ),
    )
    return basis, s
end

function _pair_block_centered_axis_support_terms(
    layer,
    exponents,
    center,
    left_interval,
    right_interval,
)
    matrices = gaussian_factor_matrices(layer; exponents, center)
    terms = Array{Float64,3}(
        undef,
        length(matrices),
        length(left_interval),
        length(right_interval),
    )
    for (term, matrix) in pairs(matrices)
        terms[term, :, :] .= matrix[left_interval, right_interval]
    end
    return terms
end

function _pair_block_expected_source_product(
    left_dims,
    right_dims,
    source_mode_ordering,
    operator_x,
    operator_y,
    operator_z,
)
    left_modes = CRPSForPairBlocks.source_mode_indices(
        left_dims;
        source_mode_ordering,
    )
    right_modes = CRPSForPairBlocks.source_mode_indices(
        right_dims;
        source_mode_ordering,
    )
    return [
        operator_x[left[1], right[1]] *
        operator_y[left[2], right[2]] *
        operator_z[left[3], right[3]]
        for left in left_modes, right in right_modes
    ]
end

function _pair_block_expected_source_overlap(
    left_dims,
    right_dims,
    source_mode_ordering,
    overlap_x,
    overlap_y,
    overlap_z,
)
    return _pair_block_expected_source_product(
        left_dims,
        right_dims,
        source_mode_ordering,
        overlap_x,
        overlap_y,
        overlap_z,
    )
end

function _pair_block_expected_source_position(
    left_dims,
    right_dims,
    source_mode_ordering,
    axis,
    overlap_x,
    overlap_y,
    overlap_z,
    position_x,
    position_y,
    position_z,
)
    operator_x, operator_y, operator_z =
        axis === :x ? (position_x, overlap_y, overlap_z) :
        axis === :y ? (overlap_x, position_y, overlap_z) :
        (overlap_x, overlap_y, position_z)
    return _pair_block_expected_source_product(
        left_dims,
        right_dims,
        source_mode_ordering,
        operator_x,
        operator_y,
        operator_z,
    )
end

function _pair_block_expected_source_x2(
    left_dims,
    right_dims,
    source_mode_ordering,
    axis,
    overlap_x,
    overlap_y,
    overlap_z,
    x2_x,
    x2_y,
    x2_z,
)
    operator_x, operator_y, operator_z =
        axis === :x ? (x2_x, overlap_y, overlap_z) :
        axis === :y ? (overlap_x, x2_y, overlap_z) :
        (overlap_x, overlap_y, x2_z)
    return _pair_block_expected_source_product(
        left_dims,
        right_dims,
        source_mode_ordering,
        operator_x,
        operator_y,
        operator_z,
    )
end

function _pair_block_expected_source_kinetic(
    left_dims,
    right_dims,
    source_mode_ordering,
    overlap_x,
    overlap_y,
    overlap_z,
    kinetic_x,
    kinetic_y,
    kinetic_z,
)
    return _pair_block_expected_source_product(
        left_dims,
        right_dims,
        source_mode_ordering,
        kinetic_x,
        overlap_y,
        overlap_z,
    ) +
           _pair_block_expected_source_product(
        left_dims,
        right_dims,
        source_mode_ordering,
        overlap_x,
        kinetic_y,
        overlap_z,
    ) +
           _pair_block_expected_source_product(
        left_dims,
        right_dims,
        source_mode_ordering,
        overlap_x,
        overlap_y,
        kinetic_z,
    )
end

function _pair_block_expected_source_electron_nuclear(
    left_dims,
    right_dims,
    source_mode_ordering,
    coefficients,
    gaussian_x,
    gaussian_y,
    gaussian_z,
)
    left_modes = CRPSForPairBlocks.source_mode_indices(
        left_dims;
        source_mode_ordering,
    )
    right_modes = CRPSForPairBlocks.source_mode_indices(
        right_dims;
        source_mode_ordering,
    )
    return [
        sum(
            -Float64(coefficients[term]) *
            gaussian_x[term, left[1], right[1]] *
            gaussian_y[term, left[2], right[2]] *
            gaussian_z[term, left[3], right[3]]
            for term in eachindex(coefficients)
        )
        for left in left_modes, right in right_modes
    ]
end

function _pair_block_expected_position(
    left_cpb,
    right_cpb,
    axis,
    overlap_x,
    overlap_y,
    overlap_z,
    position_x,
    position_y,
    position_z,
)
    operator_x, operator_y, operator_z =
        axis === :x ? (position_x, overlap_y, overlap_z) :
        axis === :y ? (overlap_x, position_y, overlap_z) :
        (overlap_x, overlap_y, position_z)
    return _pair_block_expected_product(left_cpb, right_cpb, operator_x, operator_y, operator_z)
end

function _pair_block_expected_x2(
    left_cpb,
    right_cpb,
    axis,
    overlap_x,
    overlap_y,
    overlap_z,
    x2_x,
    x2_y,
    x2_z,
)
    operator_x, operator_y, operator_z =
        axis === :x ? (x2_x, overlap_y, overlap_z) :
        axis === :y ? (overlap_x, x2_y, overlap_z) :
        (overlap_x, overlap_y, x2_z)
    return _pair_block_expected_product(left_cpb, right_cpb, operator_x, operator_y, operator_z)
end

function _pair_block_expected_kinetic(
    left_cpb,
    right_cpb,
    overlap_x,
    overlap_y,
    overlap_z,
    kinetic_x,
    kinetic_y,
    kinetic_z,
)
    return _pair_block_expected_product(
        left_cpb,
        right_cpb,
        kinetic_x,
        overlap_y,
        overlap_z,
    ) +
           _pair_block_expected_product(
        left_cpb,
        right_cpb,
        overlap_x,
        kinetic_y,
        overlap_z,
    ) +
           _pair_block_expected_product(
        left_cpb,
        right_cpb,
        overlap_x,
        overlap_y,
        kinetic_z,
    )
end

@testset "CartesianPairBlockMaterialization unavailable summary" begin
    summary = CPBM.unavailable_summary(:not_selected, :not_selected_route)

    @test summary.object_kind == :cartesian_pair_block_materialization_plan_summary
    @test summary.status == :not_selected
    @test summary.blocker == :not_selected_route
    @test summary.pair_operator_plan_count == 0
    @test summary.pair_block_record_count == 0
    @test summary.ready_record_count == 0
    @test summary.blocked_record_count == 0
    @test summary.materialization_path_counts == ()
    @test summary.readiness_status_counts == ()
    @test summary.blocker_counts == ()
    @test !summary.materialized
    @test !summary.source_operator_blocks_materialized
    @test !summary.final_pair_blocks_materialized
    @test !summary.operator_blocks_materialized
    @test !summary.hamiltonian_data_materialized
    @test !summary.artifacts_materialized
end

@testset "CartesianPairBlockMaterialization direct/direct preflight" begin
    retained_plan = _pair_block_retained_plan()
    unit_pair_plan = CUPForPairBlocks.unit_pair_plan(retained_plan)
    transform_plan =
        CRTCForPairBlocks.retained_unit_transform_contract_plan(retained_plan)
    pair_operator_plan =
        CPOPForPairBlocks.pair_operator_plan(
            unit_pair_plan,
            transform_plan;
            route_core_sidecars = false,
        )
    materialization_plan =
        CPBM.pair_block_materialization_plan(pair_operator_plan)
    materialization_summary = CPBM.summary(materialization_plan)

    @test materialization_plan isa CPBM.PairBlockMaterializationPlan
    @test length(CPBM.pair_block_materialization_records(materialization_plan)) == 3
    @test materialization_summary.status == :blocked_pair_block_materialization_plan
    @test materialization_summary.pair_operator_plan_count == 3
    @test materialization_summary.pair_block_record_count == 3
    @test materialization_summary.ready_record_count == 1
    @test materialization_summary.blocked_record_count == 2
    @test _pair_block_count(
        materialization_summary.materialization_path_counts,
        :materialization_path,
        :direct_direct_pair_block_materialization_pilot,
    ) == 1
    @test _pair_block_count(
        materialization_summary.materialization_path_counts,
        :materialization_path,
        :deferred_pair_block_materialization_path,
    ) == 1
    @test _pair_block_count(
        materialization_summary.materialization_path_counts,
        :materialization_path,
        :pqs_source_pair_preflight,
    ) == 1
    @test _pair_block_count(
        materialization_summary.readiness_status_counts,
        :readiness_status,
        :ready_metadata_only_not_materialized,
    ) == 1
    @test _pair_block_count(
        materialization_summary.readiness_status_counts,
        :readiness_status,
        :blocked_pair_block_materialization_not_implemented,
    ) == 1
    @test _pair_block_count(
        materialization_summary.readiness_status_counts,
        :readiness_status,
        :blocked_missing_raw_product_source_plan,
    ) == 1
    @test _pair_block_count(
        materialization_summary.blocker_counts,
        :blocker,
        :non_direct_direct_pair_block_materialization_not_implemented,
    ) == 1
    @test _pair_block_count(
        materialization_summary.blocker_counts,
        :blocker,
        :missing_left_raw_product_source_plan,
    ) == 1

    direct_record = _pair_block_record_for(
        materialization_plan,
        :pair_block_direct_unit,
        :pair_block_direct_unit,
    )
    direct_pqs_record = _pair_block_record_for(
        materialization_plan,
        :pair_block_direct_unit,
        :pair_block_pqs_unit,
    )
    pqs_pqs_record = _pair_block_record_for(
        materialization_plan,
        :pair_block_pqs_unit,
        :pair_block_pqs_unit,
    )

    @test direct_record.materialization_path ==
          :direct_direct_pair_block_materialization_pilot
    @test direct_record.readiness_status == :ready_metadata_only_not_materialized
    @test direct_record.blocker === nothing
    @test !direct_record.materialized

    overlap_x, overlap_y, overlap_z = _pair_block_overlap_axes()
    overlap_result = CPBM.direct_direct_overlap_block(
        direct_record;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
    )
    direct_source = only(direct_record.metadata.left_source_cpbs)
    expected_overlap =
        _pair_block_expected_overlap(
            direct_source,
            direct_source,
            overlap_x,
            overlap_y,
            overlap_z,
        )

    @test overlap_result isa CPBM.PairBlockMaterializationResult
    @test overlap_result.term == :overlap
    @test overlap_result.pair_key == (:pair_block_direct_unit, :pair_block_direct_unit)
    @test size(overlap_result.block) == (8, 8)
    @test overlap_result.block ≈ expected_overlap
    @test overlap_result.materialized
    @test overlap_result.source_operator_blocks_materialized
    @test overlap_result.final_pair_blocks_materialized
    @test !overlap_result.operator_blocks_materialized
    @test !overlap_result.hamiltonian_data_materialized
    @test !overlap_result.artifacts_materialized

    @test direct_pqs_record.readiness_status ==
          :blocked_pair_block_materialization_not_implemented
    @test direct_pqs_record.blocker ==
          :non_direct_direct_pair_block_materialization_not_implemented
    @test_throws ArgumentError CPBM.direct_direct_overlap_block(
        direct_pqs_record;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
    )
    position_x, position_y, position_z = _pair_block_position_axes()
    @test_throws ArgumentError CPBM.direct_direct_position_block(
        direct_pqs_record;
        axis = :x,
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        position_1d = (; x = position_x, y = position_y, z = position_z),
    )
    x2_x, x2_y, x2_z = _pair_block_x2_axes()
    @test_throws ArgumentError CPBM.direct_direct_x2_block(
        direct_pqs_record;
        axis = :x,
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        x2_1d = (; x = x2_x, y = x2_y, z = x2_z),
    )
    kinetic_x, kinetic_y, kinetic_z = _pair_block_kinetic_axes()
    @test_throws ArgumentError CPBM.direct_direct_kinetic_block(
        direct_pqs_record;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        kinetic_1d = (; x = kinetic_x, y = kinetic_y, z = kinetic_z),
    )
    @test pqs_pqs_record.materialization_path == :pqs_source_pair_preflight
    @test pqs_pqs_record.readiness_status ==
          :blocked_missing_raw_product_source_plan
    @test pqs_pqs_record.blocker == :missing_left_raw_product_source_plan
    @test isnothing(pqs_pqs_record.metadata.left_source_mode_dims)
    @test isnothing(pqs_pqs_record.metadata.right_source_mode_dims)
    @test !materialization_summary.materialized
    @test !materialization_summary.source_operator_blocks_materialized
    @test !materialization_summary.final_pair_blocks_materialized
    @test !materialization_summary.operator_blocks_materialized
    @test !materialization_summary.hamiltonian_data_materialized
    @test !materialization_summary.artifacts_materialized
end

@testset "CartesianPairBlockMaterialization White-Lindsey boundary-stratum preflight" begin
    retained_plan = _pair_block_white_lindsey_retained_plan()
    lw_units = CRUForPairBlocks.units(retained_plan)
    @test length(lw_units) == 3

    facet_descriptor =
        CPBM.white_lindsey_boundary_stratum_unit_adapter_descriptor(lw_units[1])
    @test facet_descriptor.object_kind ==
          :white_lindsey_boundary_stratum_unit_adapter_descriptor
    @test facet_descriptor.status ==
          :available_metadata_only_white_lindsey_unit_adapter_descriptor
    @test isnothing(facet_descriptor.blocker)
    @test facet_descriptor.unit_key == :pair_block_lw_left_unit
    @test facet_descriptor.unit_index == 1
    @test facet_descriptor.unit_kind ==
          :white_lindsey_boundary_stratum_retained_unit
    @test facet_descriptor.stratum_kind == :facet_cpb
    @test facet_descriptor.source_cpb_count == 1
    @test facet_descriptor.source_cpb_role ==
          :pair_block_lw_left_facet_source_cpb
    @test facet_descriptor.source_cpb_codimension == 1
    @test facet_descriptor.source_cpb_shape == (1, 3, 3)
    @test facet_descriptor.active_product_axes == (:y, :z)
    @test facet_descriptor.fixed_axes == (:x,)
    @test facet_descriptor.fixed_axis_coordinates ==
          ((; axis = :x, coordinate = 1),)
    @test facet_descriptor.planned_old_kernel == :_nested_face_product
    @test facet_descriptor.planned_1d_helper == :_nested_doside_1d
    @test !facet_descriptor.coefficient_maps_materialized
    @test !facet_descriptor.final_pair_blocks_materialized
    @test !facet_descriptor.hamiltonian_data_materialized

    edge_descriptor =
        CPBM.white_lindsey_boundary_stratum_unit_adapter_descriptor(lw_units[2])
    @test edge_descriptor.status ==
          :available_metadata_only_white_lindsey_unit_adapter_descriptor
    @test edge_descriptor.unit_key == :pair_block_lw_right_unit
    @test edge_descriptor.stratum_kind == :edge_cpb
    @test edge_descriptor.source_cpb_role ==
          :pair_block_lw_right_edge_source_cpb
    @test edge_descriptor.source_cpb_codimension == 2
    @test edge_descriptor.source_cpb_shape == (1, 1, 3)
    @test edge_descriptor.active_product_axes == (:z,)
    @test edge_descriptor.fixed_axes == (:x, :y)
    @test edge_descriptor.fixed_axis_coordinates ==
          ((; axis = :x, coordinate = 4), (; axis = :y, coordinate = 2))
    @test edge_descriptor.planned_old_kernel == :_nested_edge_product
    @test edge_descriptor.planned_1d_helper == :_nested_doside_1d
    @test !edge_descriptor.coefficient_maps_materialized

    corner_descriptor =
        CPBM.white_lindsey_boundary_stratum_unit_adapter_descriptor(lw_units[3])
    @test corner_descriptor.status ==
          :available_metadata_only_white_lindsey_unit_adapter_descriptor
    @test corner_descriptor.unit_key == :pair_block_lw_corner_unit
    @test corner_descriptor.stratum_kind == :corner_cpb
    @test corner_descriptor.source_cpb_role == :pair_block_lw_corner_source_cpb
    @test corner_descriptor.source_cpb_codimension == 3
    @test corner_descriptor.source_cpb_shape == (1, 1, 1)
    @test corner_descriptor.active_product_axes == ()
    @test corner_descriptor.fixed_axes == (:x, :y, :z)
    @test corner_descriptor.fixed_axis_coordinates == (
        (; axis = :x, coordinate = 4),
        (; axis = :y, coordinate = 3),
        (; axis = :z, coordinate = 3),
    )
    @test corner_descriptor.planned_old_kernel == :_nested_corner_piece
    @test isnothing(corner_descriptor.planned_1d_helper)
    @test !corner_descriptor.coefficient_maps_materialized
    @test !corner_descriptor.source_operator_blocks_materialized
    @test !corner_descriptor.final_pair_blocks_materialized
    @test !corner_descriptor.operator_blocks_materialized
    @test !corner_descriptor.hamiltonian_data_materialized
    @test !corner_descriptor.artifacts_materialized

    missing_source_unit = _pair_block_retained_unit(
        :pair_block_lw_missing_source_unit,
        4,
        :white_lindsey_boundary_stratum_retained_unit,
        :white_lindsey_boundary_strata,
        :white_lindsey_boundary_stratum_product,
        :direct_or_trivial_embedding;
        metadata = (; stratum_kind = :facet_cpb, source_cpb_index = 4),
    )
    missing_source_descriptor =
        CPBM.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            missing_source_unit,
        )
    @test missing_source_descriptor.status ==
          :blocked_white_lindsey_boundary_stratum_unit_adapter_descriptor
    @test missing_source_descriptor.blocker ==
          :white_lindsey_unit_source_cpb_count_not_one
    @test missing_source_descriptor.source_cpb_count == 0
    @test isnothing(missing_source_descriptor.source_cpb_role)
    @test isnothing(missing_source_descriptor.active_product_axes)
    @test missing_source_descriptor.planned_old_kernel == :_nested_face_product
    @test !missing_source_descriptor.coefficient_maps_materialized

    unit_pair_plan = CUPForPairBlocks.unit_pair_plan(retained_plan)
    transform_plan =
        CRTCForPairBlocks.retained_unit_transform_contract_plan(retained_plan)
    pair_operator_plan =
        CPOPForPairBlocks.pair_operator_plan(
            unit_pair_plan,
            transform_plan;
            route_core_sidecars = false,
        )
    materialization_plan =
        CPBM.pair_block_materialization_plan(pair_operator_plan)
    materialization_summary = CPBM.summary(materialization_plan)

    @test length(CPBM.pair_block_materialization_records(materialization_plan)) == 6
    @test materialization_summary.status == :blocked_pair_block_materialization_plan
    @test materialization_summary.ready_record_count == 0
    @test materialization_summary.blocked_record_count == 6
    @test _pair_block_count(
        materialization_summary.materialization_path_counts,
        :materialization_path,
        :white_lindsey_boundary_stratum_adapter_preflight,
    ) == 6
    @test _pair_block_count(
        materialization_summary.readiness_status_counts,
        :readiness_status,
        :blocked_white_lindsey_boundary_stratum_adapter_not_available,
    ) == 6
    @test _pair_block_count(
        materialization_summary.blocker_counts,
        :blocker,
        :white_lindsey_boundary_stratum_pair_block_adapter_not_materialized,
    ) == 6

    lw_records = CPBM.pair_block_materialization_records(materialization_plan)
    lw_pair_descriptors = Tuple(
        CPBM.white_lindsey_boundary_stratum_pair_adapter_descriptor(
            record,
            _pair_block_unit_pair_for(unit_pair_plan, record.pair_key...),
        )
        for record in lw_records
    )
    @test length(lw_pair_descriptors) == 6
    @test count(
        descriptor ->
            descriptor.pair_family_classification === :facet_facet,
        lw_pair_descriptors,
    ) == 1
    @test count(
        descriptor -> descriptor.pair_family_classification === :facet_edge,
        lw_pair_descriptors,
    ) == 1
    @test count(
        descriptor -> descriptor.pair_family_classification === :facet_corner,
        lw_pair_descriptors,
    ) == 1
    @test count(
        descriptor -> descriptor.pair_family_classification === :edge_edge,
        lw_pair_descriptors,
    ) == 1
    @test count(
        descriptor -> descriptor.pair_family_classification === :edge_corner,
        lw_pair_descriptors,
    ) == 1
    @test count(
        descriptor -> descriptor.pair_family_classification === :corner_corner,
        lw_pair_descriptors,
    ) == 1

    lw_facet_edge_descriptor = only(
        descriptor for descriptor in lw_pair_descriptors
        if descriptor.pair_key ==
           (:pair_block_lw_left_unit, :pair_block_lw_right_unit)
    )
    @test lw_facet_edge_descriptor.object_kind ==
          :white_lindsey_boundary_stratum_pair_adapter_descriptor
    @test lw_facet_edge_descriptor.status ==
          :available_metadata_only_white_lindsey_pair_adapter_descriptor
    @test isnothing(lw_facet_edge_descriptor.blocker)
    @test lw_facet_edge_descriptor.descriptor_source ==
          :unit_pair_retained_unit_descriptors
    @test lw_facet_edge_descriptor.pair_family_classification == :facet_edge
    @test lw_facet_edge_descriptor.left_unit_key == :pair_block_lw_left_unit
    @test lw_facet_edge_descriptor.right_unit_key == :pair_block_lw_right_unit
    @test lw_facet_edge_descriptor.left_unit_kind ==
          :white_lindsey_boundary_stratum_retained_unit
    @test lw_facet_edge_descriptor.right_unit_kind ==
          :white_lindsey_boundary_stratum_retained_unit
    @test lw_facet_edge_descriptor.left_stratum_kind == :facet_cpb
    @test lw_facet_edge_descriptor.right_stratum_kind == :edge_cpb
    @test lw_facet_edge_descriptor.left_planned_old_kernel ==
          :_nested_face_product
    @test lw_facet_edge_descriptor.right_planned_old_kernel ==
          :_nested_edge_product
    @test lw_facet_edge_descriptor.left_planned_1d_helper ==
          :_nested_doside_1d
    @test lw_facet_edge_descriptor.right_planned_1d_helper ==
          :_nested_doside_1d
    @test lw_facet_edge_descriptor.shared_1d_helper == :_nested_doside_1d
    @test lw_facet_edge_descriptor.materialization_path ==
          :white_lindsey_boundary_stratum_adapter_preflight
    @test lw_facet_edge_descriptor.pair_block_readiness_status ==
          :blocked_white_lindsey_boundary_stratum_adapter_not_available
    @test lw_facet_edge_descriptor.pair_block_blocker ==
          :white_lindsey_boundary_stratum_pair_block_adapter_not_materialized
    @test !lw_facet_edge_descriptor.coefficient_maps_materialized
    @test !lw_facet_edge_descriptor.source_operator_blocks_materialized
    @test !lw_facet_edge_descriptor.final_pair_blocks_materialized
    @test !lw_facet_edge_descriptor.operator_blocks_materialized
    @test !lw_facet_edge_descriptor.hamiltonian_data_materialized
    @test !lw_facet_edge_descriptor.artifacts_materialized

    lw_cross_record = _pair_block_record_for(
        materialization_plan,
        :pair_block_lw_left_unit,
        :pair_block_lw_right_unit,
    )
    @test lw_cross_record.source_operator_path ==
          :white_lindsey_boundary_stratum_adapter_path
    @test lw_cross_record.materialization_path ==
          :white_lindsey_boundary_stratum_adapter_preflight
    @test lw_cross_record.readiness_status ==
          :blocked_white_lindsey_boundary_stratum_adapter_not_available
    @test lw_cross_record.blocker ==
          :white_lindsey_boundary_stratum_pair_block_adapter_not_materialized
    @test lw_cross_record.transform_path.left ==
          :white_lindsey_boundary_stratum_product_contract
    @test lw_cross_record.transform_path.right ==
          :white_lindsey_boundary_stratum_product_contract
    @test lw_cross_record.realization_path.left == :identity_or_trivial_embedding
    @test lw_cross_record.realization_path.right == :identity_or_trivial_embedding
    @test lw_cross_record.final_block_path == :source_block_direct_to_final_block
    @test lw_cross_record.metadata.transform_contract_keys.left ==
          :pair_block_lw_left_unit
    @test lw_cross_record.metadata.transform_contract_keys.right ==
          :pair_block_lw_right_unit
    @test lw_cross_record.metadata.source_contract_keys.left ==
          :pair_block_lw_left_unit_contract
    @test lw_cross_record.metadata.source_contract_keys.right ==
          :pair_block_lw_right_unit_contract
    @test lw_cross_record.metadata.left_unit_kind ==
          :white_lindsey_boundary_stratum_retained_unit
    @test lw_cross_record.metadata.right_unit_kind ==
          :white_lindsey_boundary_stratum_retained_unit
    @test lw_cross_record.metadata.left_stratum_kind == :facet_cpb
    @test lw_cross_record.metadata.right_stratum_kind == :edge_cpb
    @test lw_cross_record.metadata.left_source_cpb_count == 1
    @test lw_cross_record.metadata.right_source_cpb_count == 1
    @test lw_cross_record.metadata.left_source_cpb_roles ==
          (:pair_block_lw_left_facet_source_cpb,)
    @test lw_cross_record.metadata.right_source_cpb_roles ==
          (:pair_block_lw_right_edge_source_cpb,)
    @test !lw_cross_record.materialized
    @test !lw_cross_record.metadata.source_operator_blocks_materialized
    @test !lw_cross_record.metadata.final_pair_blocks_materialized
    @test !lw_cross_record.metadata.operator_blocks_materialized
    @test !lw_cross_record.metadata.hamiltonian_data_materialized
    @test !lw_cross_record.metadata.artifacts_materialized

    lw_adapter_summary =
        CPBM.white_lindsey_boundary_stratum_adapter_summary(lw_cross_record)
    @test lw_adapter_summary.object_kind ==
          :white_lindsey_boundary_stratum_adapter_summary
    @test lw_adapter_summary.status == lw_cross_record.readiness_status
    @test lw_adapter_summary.blocker == lw_cross_record.blocker
    @test lw_adapter_summary.pair_key == lw_cross_record.pair_key
    @test lw_adapter_summary.pair_index == lw_cross_record.pair_index
    @test lw_adapter_summary.materialization_path ==
          :white_lindsey_boundary_stratum_adapter_preflight
    @test lw_adapter_summary.left_stratum_kind == :facet_cpb
    @test lw_adapter_summary.right_stratum_kind == :edge_cpb
    @test lw_adapter_summary.left_planned_kernel == :_nested_face_product
    @test lw_adapter_summary.right_planned_kernel == :_nested_edge_product
    @test lw_adapter_summary.left_side_1d_helper == :_nested_doside_1d
    @test lw_adapter_summary.right_side_1d_helper == :_nested_doside_1d
    @test lw_adapter_summary.shared_1d_helper == :_nested_doside_1d
    @test lw_adapter_summary.legacy_kernel_reuse_map.corner ==
          :_nested_corner_piece
    @test !lw_adapter_summary.source_operator_blocks_materialized
    @test !lw_adapter_summary.final_pair_blocks_materialized
    @test !lw_adapter_summary.operator_blocks_materialized
    @test !lw_adapter_summary.hamiltonian_data_materialized
    @test !lw_adapter_summary.artifacts_materialized

    record_only_pair_descriptor =
        CPBM.white_lindsey_boundary_stratum_pair_adapter_descriptor(
            lw_cross_record,
        )
    @test record_only_pair_descriptor.status ==
          :available_metadata_only_white_lindsey_pair_adapter_descriptor
    @test record_only_pair_descriptor.descriptor_source ==
          :pair_block_materialization_record_metadata
    @test record_only_pair_descriptor.pair_family_classification == :facet_edge
    @test record_only_pair_descriptor.left_planned_old_kernel ==
          :_nested_face_product
    @test record_only_pair_descriptor.right_planned_old_kernel ==
          :_nested_edge_product
    @test isnothing(record_only_pair_descriptor.left_unit_descriptor_status)
    @test !record_only_pair_descriptor.coefficient_maps_materialized

    missing_stratum_record =
        _pair_block_record_with_metadata(
            lw_cross_record,
            merge(lw_cross_record.metadata, (; left_stratum_kind = nothing)),
        )
    missing_stratum_summary =
        CPBM.white_lindsey_boundary_stratum_adapter_summary(
            missing_stratum_record,
        )
    @test missing_stratum_summary.status ==
          :blocked_white_lindsey_boundary_stratum_adapter_summary
    @test missing_stratum_summary.blocker ==
          :missing_white_lindsey_stratum_kind
    @test isnothing(missing_stratum_summary.left_planned_kernel)
    @test !missing_stratum_summary.final_pair_blocks_materialized

    missing_stratum_pair_descriptor =
        CPBM.white_lindsey_boundary_stratum_pair_adapter_descriptor(
            missing_stratum_record,
        )
    @test missing_stratum_pair_descriptor.status ==
          :blocked_white_lindsey_boundary_stratum_pair_adapter_descriptor
    @test missing_stratum_pair_descriptor.blocker ==
          :missing_white_lindsey_stratum_kind
    @test isnothing(missing_stratum_pair_descriptor.pair_family_classification)
    @test !missing_stratum_pair_descriptor.final_pair_blocks_materialized

    unknown_stratum_record =
        _pair_block_record_with_metadata(
            lw_cross_record,
            merge(lw_cross_record.metadata, (; right_stratum_kind = :ridge_cpb)),
        )
    unknown_stratum_summary =
        CPBM.white_lindsey_boundary_stratum_adapter_summary(
            unknown_stratum_record,
        )
    @test unknown_stratum_summary.status ==
          :blocked_white_lindsey_boundary_stratum_adapter_summary
    @test unknown_stratum_summary.blocker ==
          :unknown_white_lindsey_stratum_kind
    @test isnothing(unknown_stratum_summary.right_planned_kernel)
    @test !unknown_stratum_summary.final_pair_blocks_materialized

    unknown_stratum_pair_descriptor =
        CPBM.white_lindsey_boundary_stratum_pair_adapter_descriptor(
            unknown_stratum_record,
        )
    @test unknown_stratum_pair_descriptor.status ==
          :blocked_white_lindsey_boundary_stratum_pair_adapter_descriptor
    @test unknown_stratum_pair_descriptor.blocker ==
          :unknown_white_lindsey_stratum_kind
    @test isnothing(unknown_stratum_pair_descriptor.pair_family_classification)
    @test !unknown_stratum_pair_descriptor.coefficient_maps_materialized

    lw_batch_summary =
        CPBM.white_lindsey_boundary_stratum_adapter_summary(materialization_plan)
    @test lw_batch_summary.object_kind ==
          :white_lindsey_boundary_stratum_adapter_batch_summary
    @test lw_batch_summary.status ==
          :available_metadata_only_white_lindsey_adapter_reuse_batch
    @test isnothing(lw_batch_summary.blocker)
    @test lw_batch_summary.input_record_count == 6
    @test lw_batch_summary.record_count == 6
    @test lw_batch_summary.available_count == 6
    @test lw_batch_summary.blocked_count == 0
    @test lw_batch_summary.reuse_metadata_available_count == 6
    @test lw_batch_summary.reuse_metadata_blocked_count == 0
    @test lw_batch_summary.skipped_record_count == 0
    @test _pair_block_count(
        lw_batch_summary.status_counts,
        :status,
        :blocked_white_lindsey_boundary_stratum_adapter_not_available,
    ) == 6
    @test _pair_block_count(
        lw_batch_summary.blocker_counts,
        :blocker,
        :white_lindsey_boundary_stratum_pair_block_adapter_not_materialized,
    ) == 6
    @test _pair_block_count(
        lw_batch_summary.left_planned_kernel_counts,
        :planned_kernel,
        :_nested_face_product,
    ) == 3
    @test _pair_block_count(
        lw_batch_summary.left_planned_kernel_counts,
        :planned_kernel,
        :_nested_edge_product,
    ) == 2
    @test _pair_block_count(
        lw_batch_summary.left_planned_kernel_counts,
        :planned_kernel,
        :_nested_corner_piece,
    ) == 1
    @test _pair_block_count(
        lw_batch_summary.right_planned_kernel_counts,
        :planned_kernel,
        :_nested_face_product,
    ) == 1
    @test _pair_block_count(
        lw_batch_summary.right_planned_kernel_counts,
        :planned_kernel,
        :_nested_edge_product,
    ) == 2
    @test _pair_block_count(
        lw_batch_summary.right_planned_kernel_counts,
        :planned_kernel,
        :_nested_corner_piece,
    ) == 3
    @test _pair_block_count(
        lw_batch_summary.all_planned_kernel_counts,
        :planned_kernel,
        :_nested_face_product,
    ) == 4
    @test _pair_block_count(
        lw_batch_summary.all_planned_kernel_counts,
        :planned_kernel,
        :_nested_edge_product,
    ) == 4
    @test _pair_block_count(
        lw_batch_summary.all_planned_kernel_counts,
        :planned_kernel,
        :_nested_corner_piece,
    ) == 4
    @test _pair_block_count(
        lw_batch_summary.stratum_kind_counts,
        :stratum_kind,
        :facet_cpb,
    ) == 4
    @test _pair_block_count(
        lw_batch_summary.stratum_kind_counts,
        :stratum_kind,
        :edge_cpb,
    ) == 4
    @test _pair_block_count(
        lw_batch_summary.stratum_kind_counts,
        :stratum_kind,
        :corner_cpb,
    ) == 4
    @test _pair_block_count(
        lw_batch_summary.shared_1d_helper_counts,
        :helper,
        :_nested_doside_1d,
    ) == 5
    @test !lw_batch_summary.source_operator_blocks_materialized
    @test !lw_batch_summary.final_pair_blocks_materialized
    @test !lw_batch_summary.operator_blocks_materialized
    @test !lw_batch_summary.hamiltonian_data_materialized
    @test !lw_batch_summary.artifacts_materialized

    non_lw_retained_plan = _pair_block_retained_plan()
    non_lw_plan =
        CPBM.pair_block_materialization_plan(
            CPOPForPairBlocks.pair_operator_plan(
                CUPForPairBlocks.unit_pair_plan(non_lw_retained_plan),
                CRTCForPairBlocks.retained_unit_transform_contract_plan(
                    non_lw_retained_plan,
                );
                route_core_sidecars = false,
            ),
        )
    non_lw_record = first(CPBM.pair_block_materialization_records(non_lw_plan))
    mixed_batch_summary =
        CPBM.white_lindsey_boundary_stratum_adapter_summary((
            lw_cross_record,
            non_lw_record,
        ))
    @test mixed_batch_summary.input_record_count == 2
    @test mixed_batch_summary.record_count == 1
    @test mixed_batch_summary.available_count == 1
    @test mixed_batch_summary.blocked_count == 0
    @test mixed_batch_summary.reuse_metadata_available_count == 1
    @test mixed_batch_summary.reuse_metadata_blocked_count == 0
    @test mixed_batch_summary.skipped_record_count == 1
    @test _pair_block_count(
        mixed_batch_summary.skipped_blocker_counts,
        :blocker,
        :not_white_lindsey_boundary_stratum_adapter_preflight,
    ) == 1
    @test !mixed_batch_summary.final_pair_blocks_materialized
end

@testset "CartesianPairBlockMaterialization PQS source-axis transform facts" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    basis, _ = _pair_block_fixed_a_nested_test_basis(13)
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :pgdg_localized_experimental,
        refinement_levels = 0,
    )
    interval = 2:(length(basis) - 1)
    old_plan = GaussletBases._cartesian_source_box_axis_transform_plan(
        GaussletBases._CartesianNestedAxisBundles3D(bundle, bundle, bundle),
        (interval, interval, interval),
        (5, 5, 5);
        enforce_symmetric_odd = false,
    )
    transform_result =
        CPBM.pqs_source_axis_transform_facts_from_pgdg_axes(
            (; x = bundle, y = bundle, z = bundle);
            source_intervals = (interval, interval, interval),
            source_mode_dims = (5, 5, 5),
            enforce_symmetric_odd = false,
        )

    @test transform_result.object_kind == :pqs_source_axis_transform_fact_set
    @test transform_result.status == :available_pqs_source_axis_transform_facts
    @test transform_result.transform_source ==
          :repo_owned_pgdg_doside_source_axis_transform
    @test transform_result.source_mode_dims_requested == (5, 5, 5)
    @test transform_result.source_mode_dims == (5, 5, 5)
    @test !transform_result.source_mode_dims_adjusted
    @test transform_result.source_mode_count == 125
    @test transform_result.axis_coefficient_shapes ==
          ((11, 5), (11, 5), (11, 5))
    @test transform_result.max_axis_overlap_error < 1.0e-10
    @test transform_result.source_product_modes_orthogonal
    @test !transform_result.shell_realization_materialized
    @test !transform_result.lowdin_cleanup_used
    @test !transform_result.ida_data_materialized
    @test !transform_result.hamiltonian_data_materialized
    @test !transform_result.driver_route_materialized
    @test !transform_result.artifacts_materialized

    facts = transform_result.axis_transform_facts
    @test length(facts) == 3
    for axis in 1:3
        @test facts[axis] isa CRPSForPairBlocks.AxisSourceTransformFact
        @test facts[axis].axis == axis
        @test facts[axis].source_interval == interval
        @test facts[axis].source_mode_dim == 5
        @test facts[axis].coefficient_status == :materialized
        @test size(facts[axis].coefficient_matrix) == (11, 5)
        @test facts[axis].coefficient_matrix ≈
              old_plan.axes[axis].local_coefficients
        @test facts[axis].metadata.transform_source ==
              :repo_owned_pgdg_doside_source_axis_transform
        @test facts[axis].metadata.coefficient_overlap_error < 1.0e-10
        @test !facts[axis].metadata.shell_realization_materialized
        @test !facts[axis].metadata.lowdin_cleanup_used
        @test !facts[axis].metadata.ida_data_materialized
        @test !facts[axis].metadata.hamiltonian_data_materialized
    end
end

@testset "CartesianPairBlockMaterialization PQS source-pair preflight" begin
    retained_plan = _pair_block_pqs_retained_plan()
    unit_pair_plan = CUPForPairBlocks.unit_pair_plan(retained_plan)
    transform_plan =
        CRTCForPairBlocks.retained_unit_transform_contract_plan(retained_plan)
    pair_operator_plan =
        CPOPForPairBlocks.pair_operator_plan(
            unit_pair_plan,
            transform_plan;
            route_core_sidecars = false,
        )
    materialization_plan =
        CPBM.pair_block_materialization_plan(pair_operator_plan)
    materialization_summary = CPBM.summary(materialization_plan)
    pqs_cross_record = _pair_block_record_for(
        materialization_plan,
        :pair_block_pqs_left_unit,
        :pair_block_pqs_right_unit,
    )

    @test materialization_summary.ready_record_count == 3
    @test materialization_summary.blocked_record_count == 0
    @test _pair_block_count(
        materialization_summary.materialization_path_counts,
        :materialization_path,
        :pqs_source_pair_preflight,
    ) == 3
    @test pqs_cross_record.materialization_path == :pqs_source_pair_preflight
    @test pqs_cross_record.readiness_status == :ready_metadata_only_not_materialized
    @test isnothing(pqs_cross_record.blocker)
    @test !pqs_cross_record.materialized
    @test pqs_cross_record.metadata.left_raw_product_source_plan_status ==
          :available_raw_product_box_plan
    @test pqs_cross_record.metadata.right_raw_product_source_plan_status ==
          :available_raw_product_box_plan
    @test pqs_cross_record.metadata.transform_contract_keys.left ==
          :pair_block_pqs_left_unit
    @test pqs_cross_record.metadata.transform_contract_keys.right ==
          :pair_block_pqs_right_unit
    @test pqs_cross_record.metadata.source_contract_keys.left ==
          :pair_block_pqs_left_unit_contract
    @test pqs_cross_record.metadata.source_contract_keys.right ==
          :pair_block_pqs_right_unit_contract
    @test pqs_cross_record.metadata.left_source_mode_dims == (3, 3, 3)
    @test pqs_cross_record.metadata.right_source_mode_dims == (5, 4, 3)
    @test pqs_cross_record.metadata.left_source_mode_count == 27
    @test pqs_cross_record.metadata.right_source_mode_count == 60
    @test Tuple(
        fact.coefficient_status
        for fact in pqs_cross_record.metadata.left_raw_product_source_axis_transform_facts
    ) == (:materialized, :materialized, :materialized)
    @test Tuple(
        fact.coefficient_status
        for fact in pqs_cross_record.metadata.right_raw_product_source_axis_transform_facts
    ) == (:materialized, :materialized, :materialized)
    @test Tuple(
        size(fact.coefficient_matrix)
        for fact in pqs_cross_record.metadata.left_raw_product_source_axis_transform_facts
    ) == ((3, 3), (3, 3), (3, 3))
    @test Tuple(
        size(fact.coefficient_matrix)
        for fact in pqs_cross_record.metadata.right_raw_product_source_axis_transform_facts
    ) == ((3, 5), (3, 4), (3, 3))
    support_gaussian_x, support_gaussian_y, support_gaussian_z =
        _pair_block_pqs_source_support_gaussian_factor_terms()
    projected_gaussian =
        CPBM.pqs_source_pair_gaussian_factor_terms_1d(
            pqs_cross_record;
            gaussian_factor_terms_axis = (;
                x = support_gaussian_x,
                y = support_gaussian_y,
                z = support_gaussian_z,
            ),
        )
    left_facts =
        pqs_cross_record.metadata.left_raw_product_source_axis_transform_facts
    right_facts =
        pqs_cross_record.metadata.right_raw_product_source_axis_transform_facts
    @test size(projected_gaussian.x) == (2, 3, 5)
    @test size(projected_gaussian.y) == (2, 3, 4)
    @test size(projected_gaussian.z) == (2, 3, 3)
    @test projected_gaussian.x ≈ _pair_block_project_axis_terms(
        support_gaussian_x,
        left_facts[1].coefficient_matrix,
        right_facts[1].coefficient_matrix,
    )
    @test projected_gaussian.y ≈ _pair_block_project_axis_terms(
        support_gaussian_y,
        left_facts[2].coefficient_matrix,
        right_facts[2].coefficient_matrix,
    )
    @test projected_gaussian.z ≈ _pair_block_project_axis_terms(
        support_gaussian_z,
        left_facts[3].coefficient_matrix,
        right_facts[3].coefficient_matrix,
    )
    analytic_layers = (;
        x = build_basis(
            UniformBasisSpec(
                :G10;
                xmin = 1.0,
                xmax = 6.0,
                spacing = 1.0,
            ),
        ),
        y = build_basis(
            UniformBasisSpec(
                :G10;
                xmin = 1.0,
                xmax = 3.0,
                spacing = 1.0,
            ),
        ),
        z = build_basis(
            UniformBasisSpec(
                :G10;
                xmin = 1.0,
                xmax = 3.0,
                spacing = 1.0,
            ),
        ),
    )
    analytic_expansion = CoulombGaussianExpansion(
        [1.25, 0.5],
        [0.2, 0.7];
        del = 0.1,
        s = 1.0,
        c = 0.25,
        maxu = 4.0,
    )
    analytic_center = (;
        center_key = :pair_block_pqs_analytic_center,
        center_index = 1,
        location = (0.15, -0.2, 0.25),
        charge = 1.0,
    )
    centered_projected =
        CPBM.pqs_source_pair_centered_gaussian_factor_terms_1d(
            pqs_cross_record;
            axis_layers = analytic_layers,
            coulomb_expansion = analytic_expansion,
            center_record = analytic_center,
        )
    centered_support = (;
        x = _pair_block_centered_axis_support_terms(
            analytic_layers.x,
            analytic_expansion.exponents,
            analytic_center.location[1],
            left_facts[1].source_interval,
            right_facts[1].source_interval,
        ),
        y = _pair_block_centered_axis_support_terms(
            analytic_layers.y,
            analytic_expansion.exponents,
            analytic_center.location[2],
            left_facts[2].source_interval,
            right_facts[2].source_interval,
        ),
        z = _pair_block_centered_axis_support_terms(
            analytic_layers.z,
            analytic_expansion.exponents,
            analytic_center.location[3],
            left_facts[3].source_interval,
            right_facts[3].source_interval,
        ),
    )
    expected_centered =
        CPBM.pqs_source_pair_gaussian_factor_terms_1d(
            pqs_cross_record;
            gaussian_factor_terms_axis = centered_support,
        )
    @test size(centered_projected.x) == (2, 3, 5)
    @test size(centered_projected.y) == (2, 3, 4)
    @test size(centered_projected.z) == (2, 3, 3)
    @test centered_projected.x ≈ expected_centered.x
    @test centered_projected.y ≈ expected_centered.y
    @test centered_projected.z ≈ expected_centered.z
    supplied_centered_nuclear =
        CPBM.pqs_source_pair_electron_nuclear_by_center_block(
            pqs_cross_record;
            coulomb_expansion = analytic_expansion,
            center_record = analytic_center,
            gaussian_factor_terms_1d = centered_projected,
        )
    centered_nuclear =
        CPBM.pqs_source_pair_centered_electron_nuclear_by_center_block(
            pqs_cross_record;
            axis_layers = analytic_layers,
            coulomb_expansion = analytic_expansion,
            center_record = analytic_center,
        )
    @test centered_nuclear.term == :source_electron_nuclear_by_center
    @test size(centered_nuclear.block) == (27, 60)
    @test centered_nuclear.block ≈ supplied_centered_nuclear.block
    @test centered_nuclear.metadata.by_center
    @test centered_nuclear.metadata.nuclear_charge_recorded
    @test !centered_nuclear.metadata.nuclear_charge_applied
    @test !centered_nuclear.metadata.centers_summed
    @test centered_nuclear.metadata.uncharged_by_center_convention
    @test !centered_nuclear.metadata.shell_realization_materialized
    @test !centered_nuclear.metadata.lowdin_cleanup_used
    @test !centered_nuclear.metadata.ida_data_materialized
    @test !centered_nuclear.metadata.hamiltonian_data_materialized
    centered_retained =
        CPBM.pqs_source_pair_retained_centered_electron_nuclear_by_center_block(
            pqs_cross_record;
            axis_layers = analytic_layers,
            coulomb_expansion = analytic_expansion,
            center_record = analytic_center,
        )
    expected_centered_retained =
        CPBM.pqs_source_pair_retained_one_body_block(centered_nuclear)
    @test centered_retained.term == :retained_source_electron_nuclear_by_center
    @test size(centered_retained.block) == (26, 54)
    @test centered_retained.block ≈ expected_centered_retained.block
    @test centered_retained.metadata.by_center
    @test centered_retained.metadata.nuclear_charge_recorded
    @test !centered_retained.metadata.nuclear_charge_applied
    @test !centered_retained.metadata.centers_summed
    @test centered_retained.metadata.retained_source_operator_block_materialized
    @test !centered_retained.metadata.shell_realization_materialized
    @test !centered_retained.metadata.lowdin_cleanup_used

    not_materialized_retained_plan = _pair_block_mixed_pqs_retained_plan()
    not_materialized_unit_pair_plan =
        CUPForPairBlocks.unit_pair_plan(not_materialized_retained_plan)
    not_materialized_transform_plan =
        CRTCForPairBlocks.retained_unit_transform_contract_plan(
            not_materialized_retained_plan,
        )
    not_materialized_pair_operator_plan =
        CPOPForPairBlocks.pair_operator_plan(
            not_materialized_unit_pair_plan,
            not_materialized_transform_plan;
            route_core_sidecars = false,
        )
    not_materialized_plan =
        CPBM.pair_block_materialization_plan(not_materialized_pair_operator_plan)
    not_materialized_record = _pair_block_record_for(
        not_materialized_plan,
        :pair_block_selector_pqs_unit,
        :pair_block_selector_pqs_unit,
    )
    @test_throws ArgumentError CPBM.pqs_source_pair_gaussian_factor_terms_1d(
        not_materialized_record;
        gaussian_factor_terms_axis = (;
            x = support_gaussian_x,
            y = support_gaussian_y,
            z = support_gaussian_z,
        ),
    )
    @test pqs_cross_record.metadata.source_mode_ordering ==
          :x_major_y_major_z_fast
    @test pqs_cross_record.metadata.left_source_mode_ordering ==
          :x_major_y_major_z_fast
    @test pqs_cross_record.metadata.right_source_mode_ordering ==
          :x_major_y_major_z_fast
    @test pqs_cross_record.metadata.left_raw_product_source_summary.source_mode_dims ==
          (3, 3, 3)
    @test pqs_cross_record.metadata.right_raw_product_source_summary.source_mode_dims ==
          (5, 4, 3)
    @test pqs_cross_record.metadata.left_raw_product_source_retained_rule isa
          CRPSForPairBlocks.PQSBoundaryProductModeRetainedRule
    @test pqs_cross_record.metadata.right_raw_product_source_retained_rule isa
          CRPSForPairBlocks.PQSBoundaryProductModeRetainedRule
    @test pqs_cross_record.metadata.left_raw_product_source_retained_rule.retained_count ==
          26
    @test pqs_cross_record.metadata.right_raw_product_source_retained_rule.retained_count ==
          54
    @test pqs_cross_record.metadata.left_raw_product_source_retained_rule_summary.retained_count ==
          26
    @test pqs_cross_record.metadata.right_raw_product_source_retained_rule_summary.retained_count ==
          54
    @test !pqs_cross_record.metadata.source_operator_blocks_materialized
    @test !pqs_cross_record.metadata.final_pair_blocks_materialized
    @test !materialization_summary.materialized
    @test !materialization_summary.source_operator_blocks_materialized
    @test !materialization_summary.final_pair_blocks_materialized
    @test !materialization_summary.operator_blocks_materialized
    @test !materialization_summary.hamiltonian_data_materialized
    @test !materialization_summary.artifacts_materialized

    left_source_dims = pqs_cross_record.metadata.left_source_mode_dims
    right_source_dims = pqs_cross_record.metadata.right_source_mode_dims
    pqs_overlap_x, pqs_overlap_y, pqs_overlap_z =
        _pair_block_pqs_source_overlap_axes(left_source_dims, right_source_dims)
    pqs_overlap_result = CPBM.pqs_source_pair_overlap_block(
        pqs_cross_record;
        overlap_1d = (;
            x = pqs_overlap_x,
            y = pqs_overlap_y,
            z = pqs_overlap_z,
        ),
    )
    expected_pqs_overlap = _pair_block_expected_source_overlap(
        left_source_dims,
        right_source_dims,
        pqs_cross_record.metadata.source_mode_ordering,
        pqs_overlap_x,
        pqs_overlap_y,
        pqs_overlap_z,
    )
    @test pqs_overlap_result isa CPBM.PairBlockMaterializationResult
    @test pqs_overlap_result.term == :source_overlap
    @test pqs_overlap_result.pair_key ==
          (:pair_block_pqs_left_unit, :pair_block_pqs_right_unit)
    @test size(pqs_overlap_result.block) == (27, 60)
    @test pqs_overlap_result.block ≈ expected_pqs_overlap
    @test pqs_overlap_result.materialized
    @test pqs_overlap_result.source_operator_blocks_materialized
    @test !pqs_overlap_result.final_pair_blocks_materialized
    @test !pqs_overlap_result.operator_blocks_materialized
    @test !pqs_overlap_result.hamiltonian_data_materialized
    @test !pqs_overlap_result.artifacts_materialized
    @test pqs_overlap_result.metadata.block_space == :raw_product_source_modes
    @test pqs_overlap_result.metadata.left_source_mode_dims == left_source_dims
    @test pqs_overlap_result.metadata.right_source_mode_dims == right_source_dims
    @test pqs_overlap_result.metadata.source_mode_ordering ==
          :x_major_y_major_z_fast
    @test pqs_overlap_result.metadata.transform_contract_keys.left ==
          :pair_block_pqs_left_unit
    @test pqs_overlap_result.metadata.transform_contract_keys.right ==
          :pair_block_pqs_right_unit
    @test pqs_overlap_result.metadata.transform_paths.left ==
          :pqs_source_modes_boundary_selection_shell_realization_contract
    @test pqs_overlap_result.metadata.transform_paths.right ==
          :pqs_source_modes_boundary_selection_shell_realization_contract
    @test pqs_overlap_result.metadata.realization_paths.left ==
          :shell_projection_lowdin_planned
    @test pqs_overlap_result.metadata.realization_paths.right ==
          :shell_projection_lowdin_planned
    @test pqs_overlap_result.metadata.source_operator_blocks_materialized
    @test !pqs_overlap_result.metadata.final_pair_blocks_materialized
    @test !pqs_overlap_result.metadata.shell_realization_materialized
    @test pqs_overlap_result.metadata.left_raw_product_source_retained_rule ===
          pqs_cross_record.metadata.left_raw_product_source_retained_rule
    @test pqs_overlap_result.metadata.right_raw_product_source_retained_rule ===
          pqs_cross_record.metadata.right_raw_product_source_retained_rule

    left_retained_columns = CRPSForPairBlocks.retained_column_indices(
        pqs_cross_record.metadata.left_raw_product_source_retained_rule,
    )
    right_retained_columns = CRPSForPairBlocks.retained_column_indices(
        pqs_cross_record.metadata.right_raw_product_source_retained_rule,
    )
    retained_overlap_result =
        CPBM.pqs_source_pair_retained_one_body_block(pqs_overlap_result)
    wrapper_retained_overlap =
        CPBM.pqs_source_pair_retained_overlap_block(
            pqs_cross_record;
            overlap_1d = (;
                x = pqs_overlap_x,
                y = pqs_overlap_y,
                z = pqs_overlap_z,
            ),
        )
    @test retained_overlap_result.term == :retained_source_overlap
    @test size(retained_overlap_result.block) == (26, 54)
    @test retained_overlap_result.block ≈
          pqs_overlap_result.block[left_retained_columns, right_retained_columns]
    overlap_direct_oracle_maxdiff =
        maximum(abs.(wrapper_retained_overlap.block .- retained_overlap_result.block))
    @test wrapper_retained_overlap.block ≈ retained_overlap_result.block
    @test overlap_direct_oracle_maxdiff == 0.0
    @test retained_overlap_result.materialized
    @test retained_overlap_result.source_operator_blocks_materialized
    @test !retained_overlap_result.final_pair_blocks_materialized
    @test !retained_overlap_result.operator_blocks_materialized
    @test !retained_overlap_result.hamiltonian_data_materialized
    @test !retained_overlap_result.artifacts_materialized
    @test retained_overlap_result.metadata.block_space == :retained_pqs_source_modes
    @test retained_overlap_result.metadata.source_block_space ==
          :raw_product_source_modes
    @test retained_overlap_result.metadata.source_space_input_used
    @test wrapper_retained_overlap.metadata.retained_direct_boundary_product_used
    @test !wrapper_retained_overlap.metadata.source_space_input_used
    @test !wrapper_retained_overlap.metadata.raw_source_operator_block_materialized
    @test retained_overlap_result.metadata.retained_transform_kind ==
          :source_mode_column_selector
    @test wrapper_retained_overlap.metadata.retained_transform_kind ==
          :source_mode_column_selector
    @test retained_overlap_result.metadata.left_retained_count == 26
    @test retained_overlap_result.metadata.right_retained_count == 54
    @test retained_overlap_result.metadata.retained_source_operator_block_materialized
    @test !retained_overlap_result.metadata.shell_realization_materialized
    @test !retained_overlap_result.metadata.lowdin_cleanup_used

    shell_bridge =
        CPBM.pqs_source_pair_shell_realization_bridge_summary(pqs_overlap_result)
    @test shell_bridge.object_kind ==
          :pqs_source_pair_shell_realization_bridge_summary
    @test shell_bridge.status ==
          :available_metadata_only_shell_realization_bridge
    @test isnothing(shell_bridge.blocker)
    @test shell_bridge.source_block_term == :source_overlap
    @test shell_bridge.source_block_status == :source_operator_block_materialized
    @test shell_bridge.block_space == :raw_product_source_modes
    @test shell_bridge.left_source_mode_dims == left_source_dims
    @test shell_bridge.right_source_mode_dims == right_source_dims
    @test shell_bridge.left_source_mode_count == 27
    @test shell_bridge.right_source_mode_count == 60
    @test shell_bridge.source_mode_ordering == :x_major_y_major_z_fast
    @test shell_bridge.transform_contract_keys.left == :pair_block_pqs_left_unit
    @test shell_bridge.transform_contract_keys.right == :pair_block_pqs_right_unit
    @test shell_bridge.source_contract_keys.left ==
          :pair_block_pqs_left_unit_contract
    @test shell_bridge.source_contract_keys.right ==
          :pair_block_pqs_right_unit_contract
    @test shell_bridge.realization_paths.left == :shell_projection_lowdin_planned
    @test shell_bridge.realization_paths.right == :shell_projection_lowdin_planned
    @test shell_bridge.source_operator_blocks_materialized
    @test !shell_bridge.final_pair_blocks_materialized
    @test !shell_bridge.shell_realization_materialized
    @test !shell_bridge.operator_blocks_materialized
    @test !shell_bridge.hamiltonian_data_materialized
    @test !shell_bridge.artifacts_materialized

    final_readiness =
        CPBM.pqs_source_pair_final_block_readiness_summary(shell_bridge)
    @test final_readiness.object_kind ==
          :pqs_source_pair_final_block_readiness_summary
    @test final_readiness.status == :blocked_final_pqs_pair_block_not_ready
    @test final_readiness.blocker == :shell_realization_not_materialized
    @test final_readiness.bridge_status ==
          :available_metadata_only_shell_realization_bridge
    @test final_readiness.pair_key == shell_bridge.pair_key
    @test final_readiness.source_block_term == :source_overlap
    @test final_readiness.left_source_mode_dims == left_source_dims
    @test final_readiness.right_source_mode_dims == right_source_dims
    @test final_readiness.transform_contract_keys.left ==
          :pair_block_pqs_left_unit
    @test final_readiness.realization_paths.left ==
          :shell_projection_lowdin_planned
    @test final_readiness.source_operator_blocks_materialized
    @test !final_readiness.shell_realization_materialized
    @test !final_readiness.final_pair_blocks_materialized
    @test !final_readiness.operator_blocks_materialized
    @test !final_readiness.hamiltonian_data_materialized
    @test !final_readiness.artifacts_materialized

    non_bridge_readiness =
        CPBM.pqs_source_pair_final_block_readiness_summary(
            (; object_kind = :not_a_pqs_source_bridge_summary),
        )
    @test non_bridge_readiness.status ==
          :blocked_final_pqs_pair_block_not_ready
    @test non_bridge_readiness.blocker ==
          :not_pqs_source_shell_realization_bridge_summary
    @test !non_bridge_readiness.final_pair_blocks_materialized
    @test !non_bridge_readiness.shell_realization_materialized

    missing_ordering_bridge = merge(shell_bridge, (; source_mode_ordering = nothing))
    missing_ordering_readiness =
        CPBM.pqs_source_pair_final_block_readiness_summary(
            missing_ordering_bridge,
        )
    @test missing_ordering_readiness.status ==
          :blocked_final_pqs_pair_block_not_ready
    @test missing_ordering_readiness.blocker == :missing_source_mode_ordering
    @test !missing_ordering_readiness.final_pair_blocks_materialized

    synthetic_blocked_bridge =
        merge(
            shell_bridge,
            (;
                status = :blocked_missing_shell_realization_facts,
                blocker = :synthetic_bridge_blocker,
            ),
        )
    synthetic_blocked_readiness =
        CPBM.pqs_source_pair_final_block_readiness_summary(
            synthetic_blocked_bridge,
        )
    @test synthetic_blocked_readiness.status ==
          :blocked_final_pqs_pair_block_not_ready
    @test synthetic_blocked_readiness.blocker == :synthetic_bridge_blocker

    missing_bridge_realization_result =
        _pair_block_result_with_metadata(
            pqs_overlap_result,
            merge(pqs_overlap_result.metadata, (; realization_paths = nothing)),
        )
    missing_bridge_realization =
        CPBM.pqs_source_pair_shell_realization_bridge_summary(
            missing_bridge_realization_result,
        )
    @test missing_bridge_realization.status ==
          :blocked_missing_shell_realization_facts
    @test missing_bridge_realization.blocker == :missing_shell_realization_path
    @test !missing_bridge_realization.final_pair_blocks_materialized
    @test !missing_bridge_realization.shell_realization_materialized
    missing_realization_readiness =
        CPBM.pqs_source_pair_final_block_readiness_summary(
            missing_bridge_realization,
        )
    @test missing_realization_readiness.status ==
          :blocked_final_pqs_pair_block_not_ready
    @test missing_realization_readiness.blocker == :missing_shell_realization_path
    @test !missing_realization_readiness.final_pair_blocks_materialized

    missing_bridge_source_metadata_result =
        _pair_block_result_with_metadata(
            pqs_overlap_result,
            merge(pqs_overlap_result.metadata, (; left_source_mode_dims = nothing)),
        )
    missing_bridge_source_metadata =
        CPBM.pqs_source_pair_shell_realization_bridge_summary(
            missing_bridge_source_metadata_result,
        )
    missing_source_metadata_readiness =
        CPBM.pqs_source_pair_final_block_readiness_summary(
            missing_bridge_source_metadata,
        )
    @test missing_source_metadata_readiness.status ==
          :blocked_final_pqs_pair_block_not_ready
    @test missing_source_metadata_readiness.blocker ==
          :missing_left_source_mode_dims
    @test !missing_source_metadata_readiness.final_pair_blocks_materialized

    missing_bridge_contract_result =
        _pair_block_result_with_metadata(
            pqs_overlap_result,
            merge(pqs_overlap_result.metadata, (; transform_contract_keys = nothing)),
        )
    missing_bridge_contract =
        CPBM.pqs_source_pair_shell_realization_bridge_summary(
            missing_bridge_contract_result,
        )
    @test missing_bridge_contract.status ==
          :blocked_missing_shell_realization_facts
    @test missing_bridge_contract.blocker == :missing_transform_contract_keys
    @test !missing_bridge_contract.final_pair_blocks_materialized
    @test !missing_bridge_contract.shell_realization_materialized
    missing_contract_readiness =
        CPBM.pqs_source_pair_final_block_readiness_summary(missing_bridge_contract)
    @test missing_contract_readiness.status ==
          :blocked_final_pqs_pair_block_not_ready
    @test missing_contract_readiness.blocker == :missing_transform_contract_keys

    pqs_position_x, pqs_position_y, pqs_position_z =
        _pair_block_pqs_source_position_axes(left_source_dims, right_source_dims)
    pqs_x2_x, pqs_x2_y, pqs_x2_z =
        _pair_block_pqs_source_x2_axes(left_source_dims, right_source_dims)
    pqs_kinetic_x, pqs_kinetic_y, pqs_kinetic_z =
        _pair_block_pqs_source_kinetic_axes(left_source_dims, right_source_dims)

    for axis in (:x, :y, :z)
        position_overlap_1d =
            axis === :y ?
            (pqs_overlap_x, pqs_overlap_y, pqs_overlap_z) :
            (; x = pqs_overlap_x, y = pqs_overlap_y, z = pqs_overlap_z)
        position_1d =
            axis === :y ?
            (pqs_position_x, pqs_position_y, pqs_position_z) :
            (; x = pqs_position_x, y = pqs_position_y, z = pqs_position_z)
        position_result = CPBM.pqs_source_pair_position_block(
            pqs_cross_record;
            axis,
            overlap_1d = position_overlap_1d,
            position_1d,
        )
        expected_position = _pair_block_expected_source_position(
            left_source_dims,
            right_source_dims,
            pqs_cross_record.metadata.source_mode_ordering,
            axis,
            pqs_overlap_x,
            pqs_overlap_y,
            pqs_overlap_z,
            pqs_position_x,
            pqs_position_y,
            pqs_position_z,
        )

        @test position_result.term == Symbol(:source_position_, axis)
        @test size(position_result.block) == (27, 60)
        @test position_result.block ≈ expected_position
        @test position_result.materialized
        @test position_result.source_operator_blocks_materialized
        @test !position_result.final_pair_blocks_materialized
        @test !position_result.operator_blocks_materialized
        @test !position_result.hamiltonian_data_materialized
        @test !position_result.artifacts_materialized
        @test position_result.metadata.block_space == :raw_product_source_modes
        @test position_result.metadata.position_axis == axis
        @test position_result.metadata.left_source_mode_dims == left_source_dims
        @test position_result.metadata.right_source_mode_dims == right_source_dims
        @test position_result.metadata.source_operator_blocks_materialized
        @test !position_result.metadata.final_pair_blocks_materialized

        x2_1d =
            axis === :z ?
            (pqs_x2_x, pqs_x2_y, pqs_x2_z) :
            (; x = pqs_x2_x, y = pqs_x2_y, z = pqs_x2_z)
        x2_result = CPBM.pqs_source_pair_x2_block(
            pqs_cross_record;
            axis,
            overlap_1d = (; x = pqs_overlap_x, y = pqs_overlap_y, z = pqs_overlap_z),
            x2_1d,
        )
        expected_x2 = _pair_block_expected_source_x2(
            left_source_dims,
            right_source_dims,
            pqs_cross_record.metadata.source_mode_ordering,
            axis,
            pqs_overlap_x,
            pqs_overlap_y,
            pqs_overlap_z,
            pqs_x2_x,
            pqs_x2_y,
            pqs_x2_z,
        )

        @test x2_result.term == Symbol(:source_x2_, axis)
        @test size(x2_result.block) == (27, 60)
        @test x2_result.block ≈ expected_x2
        @test x2_result.materialized
        @test x2_result.source_operator_blocks_materialized
        @test !x2_result.final_pair_blocks_materialized
        @test !x2_result.operator_blocks_materialized
        @test !x2_result.hamiltonian_data_materialized
        @test !x2_result.artifacts_materialized
        @test x2_result.metadata.block_space == :raw_product_source_modes
        @test x2_result.metadata.x2_axis == axis
        @test x2_result.metadata.left_source_mode_dims == left_source_dims
        @test x2_result.metadata.right_source_mode_dims == right_source_dims
        @test x2_result.metadata.source_operator_blocks_materialized
        @test !x2_result.metadata.final_pair_blocks_materialized
    end

    kinetic_result = CPBM.pqs_source_pair_kinetic_block(
        pqs_cross_record;
        overlap_1d = (; x = pqs_overlap_x, y = pqs_overlap_y, z = pqs_overlap_z),
        kinetic_1d = (pqs_kinetic_x, pqs_kinetic_y, pqs_kinetic_z),
    )
    expected_kinetic = _pair_block_expected_source_kinetic(
        left_source_dims,
        right_source_dims,
        pqs_cross_record.metadata.source_mode_ordering,
        pqs_overlap_x,
        pqs_overlap_y,
        pqs_overlap_z,
        pqs_kinetic_x,
        pqs_kinetic_y,
        pqs_kinetic_z,
    )
    @test kinetic_result.term == :source_kinetic
    @test size(kinetic_result.block) == (27, 60)
    @test kinetic_result.block ≈ expected_kinetic
    @test kinetic_result.materialized
    @test kinetic_result.source_operator_blocks_materialized
    @test !kinetic_result.final_pair_blocks_materialized
    @test !kinetic_result.operator_blocks_materialized
    @test !kinetic_result.hamiltonian_data_materialized
    @test !kinetic_result.artifacts_materialized
    @test kinetic_result.metadata.block_space == :raw_product_source_modes
    @test kinetic_result.metadata.left_source_mode_dims == left_source_dims
    @test kinetic_result.metadata.right_source_mode_dims == right_source_dims
    @test kinetic_result.metadata.kinetic_factor_form == (
        (:kinetic, :overlap, :overlap),
        (:overlap, :kinetic, :overlap),
        (:overlap, :overlap, :kinetic),
    )
    @test kinetic_result.metadata.source_operator_blocks_materialized
    @test !kinetic_result.metadata.final_pair_blocks_materialized

    retained_kinetic_result =
        CPBM.pqs_source_pair_retained_one_body_block(
            kinetic_result,
            pqs_cross_record.metadata.left_raw_product_source_retained_rule,
            pqs_cross_record.metadata.right_raw_product_source_retained_rule,
        )
    wrapper_retained_kinetic =
        CPBM.pqs_source_pair_retained_kinetic_block(
            pqs_cross_record;
            overlap_1d = (; x = pqs_overlap_x, y = pqs_overlap_y, z = pqs_overlap_z),
            kinetic_1d = (pqs_kinetic_x, pqs_kinetic_y, pqs_kinetic_z),
        )
    @test retained_kinetic_result.term == :retained_source_kinetic
    @test size(retained_kinetic_result.block) == (26, 54)
    @test retained_kinetic_result.block ≈
          kinetic_result.block[left_retained_columns, right_retained_columns]
    kinetic_direct_oracle_maxdiff =
        maximum(abs.(wrapper_retained_kinetic.block .- retained_kinetic_result.block))
    @test wrapper_retained_kinetic.block ≈ retained_kinetic_result.block
    @test kinetic_direct_oracle_maxdiff == 0.0
    @test retained_kinetic_result.metadata.block_space == :retained_pqs_source_modes
    @test retained_kinetic_result.metadata.source_block_term == :source_kinetic
    @test retained_kinetic_result.metadata.source_space_input_used
    @test wrapper_retained_kinetic.metadata.retained_direct_boundary_product_used
    @test !wrapper_retained_kinetic.metadata.source_space_input_used
    @test !wrapper_retained_kinetic.metadata.raw_source_operator_block_materialized
    @test retained_kinetic_result.metadata.retained_source_operator_block_materialized
    @test !retained_kinetic_result.metadata.shell_realization_materialized
    @test !retained_kinetic_result.metadata.lowdin_cleanup_used

    gaussian_x, gaussian_y, gaussian_z =
        _pair_block_pqs_source_gaussian_factor_terms(
            left_source_dims,
            right_source_dims,
        )
    expansion = CoulombGaussianExpansion(
        [0.7, 0.2],
        [0.5, 1.5];
        del = 0.6,
        s = 0.5,
        c = 0.03,
        maxu = 1.2,
    )
    center_record = (;
        center_key = :synthetic_origin,
        center_index = 1,
        charge = 2.0,
        location = (0.0, 0.0, 0.0),
    )
    nuclear_result = CPBM.pqs_source_pair_electron_nuclear_by_center_block(
        pqs_cross_record;
        coulomb_expansion = expansion,
        center_record,
        gaussian_factor_terms_1d = (;
            x = gaussian_x,
            y = gaussian_y,
            z = gaussian_z,
        ),
    )
    expected_nuclear = _pair_block_expected_source_electron_nuclear(
        left_source_dims,
        right_source_dims,
        pqs_cross_record.metadata.source_mode_ordering,
        expansion.coefficients,
        gaussian_x,
        gaussian_y,
        gaussian_z,
    )
    @test nuclear_result.term == :source_electron_nuclear_by_center
    @test size(nuclear_result.block) == (27, 60)
    @test nuclear_result.block ≈ expected_nuclear
    @test nuclear_result.metadata.block_space == :raw_product_source_modes
    @test nuclear_result.metadata.by_center
    @test nuclear_result.metadata.center_key == :synthetic_origin
    @test nuclear_result.metadata.center_index == 1
    @test nuclear_result.metadata.center_location == (0.0, 0.0, 0.0)
    @test nuclear_result.metadata.nuclear_charge == 2.0
    @test nuclear_result.metadata.nuclear_charge_recorded
    @test !nuclear_result.metadata.nuclear_charge_applied
    @test !nuclear_result.metadata.centers_summed
    @test nuclear_result.metadata.uncharged_by_center_convention
    @test nuclear_result.metadata.gaussian_term_count == 2
    @test nuclear_result.metadata.axis_factor_term_shapes ==
          (x = size(gaussian_x), y = size(gaussian_y), z = size(gaussian_z))
    @test nuclear_result.source_operator_blocks_materialized
    @test !nuclear_result.final_pair_blocks_materialized
    @test !nuclear_result.operator_blocks_materialized
    @test !nuclear_result.hamiltonian_data_materialized
    @test !nuclear_result.artifacts_materialized
    @test !nuclear_result.metadata.shell_realization_materialized
    @test !nuclear_result.metadata.lowdin_cleanup_used
    @test !nuclear_result.metadata.ida_data_materialized
    @test !nuclear_result.metadata.driver_route_materialized

    retained_nuclear_result =
        CPBM.pqs_source_pair_retained_one_body_block(nuclear_result)
    wrapper_retained_nuclear =
        CPBM.pqs_source_pair_retained_electron_nuclear_by_center_block(
            pqs_cross_record;
            coulomb_expansion = expansion,
            center_record,
            gaussian_factor_terms_1d = (gaussian_x, gaussian_y, gaussian_z),
        )
    @test retained_nuclear_result.term ==
          :retained_source_electron_nuclear_by_center
    @test size(retained_nuclear_result.block) == (26, 54)
    @test retained_nuclear_result.block ≈
          nuclear_result.block[left_retained_columns, right_retained_columns]
    @test wrapper_retained_nuclear.block ≈ retained_nuclear_result.block
    @test retained_nuclear_result.metadata.block_space ==
          :retained_pqs_source_modes
    @test retained_nuclear_result.metadata.source_block_term ==
          :source_electron_nuclear_by_center
    @test retained_nuclear_result.metadata.by_center
    @test retained_nuclear_result.metadata.nuclear_charge_recorded
    @test !retained_nuclear_result.metadata.nuclear_charge_applied
    @test !retained_nuclear_result.metadata.centers_summed
    @test retained_nuclear_result.metadata.retained_source_operator_block_materialized
    @test !retained_nuclear_result.metadata.shell_realization_materialized
    @test !retained_nuclear_result.metadata.lowdin_cleanup_used

    bad_pqs_overlap_x = zeros(Float64, left_source_dims[1] + 1, right_source_dims[1])
    @test_throws ArgumentError CPBM.pqs_source_pair_overlap_block(
        pqs_cross_record;
        overlap_1d = (;
            x = bad_pqs_overlap_x,
            y = pqs_overlap_y,
            z = pqs_overlap_z,
        ),
    )
    bad_pqs_position_y =
        zeros(Float64, left_source_dims[2], right_source_dims[2] + 1)
    @test_throws ArgumentError CPBM.pqs_source_pair_position_block(
        pqs_cross_record;
        axis = :x,
        overlap_1d = (; x = pqs_overlap_x, y = pqs_overlap_y, z = pqs_overlap_z),
        position_1d = (;
            x = pqs_position_x,
            y = bad_pqs_position_y,
            z = pqs_position_z,
        ),
    )
    bad_pqs_x2_x = zeros(Float64, left_source_dims[1] + 1, right_source_dims[1])
    @test_throws ArgumentError CPBM.pqs_source_pair_x2_block(
        pqs_cross_record;
        axis = :x,
        overlap_1d = (; x = pqs_overlap_x, y = pqs_overlap_y, z = pqs_overlap_z),
        x2_1d = (; x = bad_pqs_x2_x, y = pqs_x2_y, z = pqs_x2_z),
    )
    bad_pqs_kinetic_z =
        zeros(Float64, left_source_dims[3], right_source_dims[3] + 1)
    @test_throws ArgumentError CPBM.pqs_source_pair_kinetic_block(
        pqs_cross_record;
        overlap_1d = (; x = pqs_overlap_x, y = pqs_overlap_y, z = pqs_overlap_z),
        kinetic_1d = (;
            x = pqs_kinetic_x,
            y = pqs_kinetic_y,
            z = bad_pqs_kinetic_z,
        ),
    )

    left_contract, right_contract =
        CRTCForPairBlocks.transform_contracts(transform_plan)
    missing_left_contract =
        _pair_block_transform_contract_with_metadata(
            left_contract,
            merge(
                left_contract.metadata,
                (;
                    raw_product_source_plan = nothing,
                    raw_product_source_summary =
                        CRPSForPairBlocks.unavailable_summary(
                            :blocked_missing_source_mode_dims,
                            :missing_pqs_source_mode_dims,
                        ),
                    raw_product_source_plan_status =
                        :blocked_missing_source_mode_dims,
                ),
            ),
        )
    missing_transform_plan =
        _pair_block_transform_plan_with_contracts(
            transform_plan,
            (missing_left_contract, right_contract),
        )
    missing_pair_operator_plan =
        CPOPForPairBlocks.pair_operator_plan(
            unit_pair_plan,
            missing_transform_plan;
            route_core_sidecars = false,
        )
    missing_materialization_plan =
        CPBM.pair_block_materialization_plan(missing_pair_operator_plan)
    missing_cross_record = _pair_block_record_for(
        missing_materialization_plan,
        :pair_block_pqs_left_unit,
        :pair_block_pqs_right_unit,
    )

    @test missing_cross_record.materialization_path == :pqs_source_pair_preflight
    @test missing_cross_record.readiness_status ==
          :blocked_missing_raw_product_source_plan
    @test missing_cross_record.blocker == :missing_left_raw_product_source_plan
    @test missing_cross_record.metadata.left_raw_product_source_plan_status ==
          :blocked_missing_source_mode_dims
    @test missing_cross_record.metadata.right_raw_product_source_plan_status ==
          :available_raw_product_box_plan
    @test_throws ArgumentError CPBM.pqs_source_pair_overlap_block(
        missing_cross_record;
        overlap_1d = (;
            x = pqs_overlap_x,
            y = pqs_overlap_y,
            z = pqs_overlap_z,
        ),
    )
    @test_throws ArgumentError CPBM.pqs_source_pair_position_block(
        missing_cross_record;
        axis = :x,
        overlap_1d = (; x = pqs_overlap_x, y = pqs_overlap_y, z = pqs_overlap_z),
        position_1d = (;
            x = pqs_position_x,
            y = pqs_position_y,
            z = pqs_position_z,
        ),
    )
    @test_throws ArgumentError CPBM.pqs_source_pair_kinetic_block(
        missing_cross_record;
        overlap_1d = (; x = pqs_overlap_x, y = pqs_overlap_y, z = pqs_overlap_z),
        kinetic_1d = (;
            x = pqs_kinetic_x,
            y = pqs_kinetic_y,
            z = pqs_kinetic_z,
        ),
    )

    incompatible_right_summary =
        merge(
            right_contract.metadata.raw_product_source_summary,
            (; source_mode_ordering = :synthetic_other_ordering),
        )
    incompatible_right_contract =
        _pair_block_transform_contract_with_metadata(
            right_contract,
            merge(
                right_contract.metadata,
                (; raw_product_source_summary = incompatible_right_summary),
            ),
        )
    incompatible_transform_plan =
        _pair_block_transform_plan_with_contracts(
            transform_plan,
            (left_contract, incompatible_right_contract),
        )
    incompatible_pair_operator_plan =
        CPOPForPairBlocks.pair_operator_plan(
            unit_pair_plan,
            incompatible_transform_plan;
            route_core_sidecars = false,
        )
    incompatible_materialization_plan =
        CPBM.pair_block_materialization_plan(incompatible_pair_operator_plan)
    incompatible_cross_record = _pair_block_record_for(
        incompatible_materialization_plan,
        :pair_block_pqs_left_unit,
        :pair_block_pqs_right_unit,
    )

    @test incompatible_cross_record.materialization_path ==
          :pqs_source_pair_preflight
    @test incompatible_cross_record.readiness_status ==
          :blocked_incompatible_raw_product_source_ordering
    @test incompatible_cross_record.blocker ==
          :incompatible_raw_product_source_ordering
    @test incompatible_cross_record.metadata.left_source_mode_ordering ==
          :x_major_y_major_z_fast
    @test incompatible_cross_record.metadata.right_source_mode_ordering ==
          :synthetic_other_ordering
    @test isnothing(incompatible_cross_record.metadata.source_mode_ordering)
    @test_throws ArgumentError CPBM.pqs_source_pair_overlap_block(
        incompatible_cross_record;
        overlap_1d = (;
            x = pqs_overlap_x,
            y = pqs_overlap_y,
            z = pqs_overlap_z,
        ),
    )
    @test_throws ArgumentError CPBM.pqs_source_pair_x2_block(
        incompatible_cross_record;
        axis = :x,
        overlap_1d = (; x = pqs_overlap_x, y = pqs_overlap_y, z = pqs_overlap_z),
        x2_1d = (; x = pqs_x2_x, y = pqs_x2_y, z = pqs_x2_z),
    )
    @test_throws ArgumentError CPBM.pqs_source_pair_kinetic_block(
        incompatible_cross_record;
        overlap_1d = (; x = pqs_overlap_x, y = pqs_overlap_y, z = pqs_overlap_z),
        kinetic_1d = (;
            x = pqs_kinetic_x,
            y = pqs_kinetic_y,
            z = pqs_kinetic_z,
        ),
    )
end

@testset "CartesianPairBlockMaterialization PQS source-pair batch selectors" begin
    retained_plan = _pair_block_mixed_pqs_retained_plan()
    unit_pair_plan = CUPForPairBlocks.unit_pair_plan(retained_plan)
    transform_plan =
        CRTCForPairBlocks.retained_unit_transform_contract_plan(retained_plan)
    pair_operator_plan =
        CPOPForPairBlocks.pair_operator_plan(
            unit_pair_plan,
            transform_plan;
            route_core_sidecars = false,
        )
    materialization_plan =
        CPBM.pair_block_materialization_plan(pair_operator_plan)
    pqs_record = _pair_block_record_for(
        materialization_plan,
        :pair_block_selector_pqs_unit,
        :pair_block_selector_pqs_unit,
    )

    @test pqs_record.materialization_path == :pqs_source_pair_preflight
    @test pqs_record.readiness_status == :ready_metadata_only_not_materialized
    @test pqs_record.metadata.left_source_mode_dims == (3, 3, 3)
    @test pqs_record.metadata.right_source_mode_dims == (3, 3, 3)

    overlap_descriptor = CPBM._pqs_source_safe_term_descriptor(:overlap)
    position_descriptor = CPBM._pqs_source_safe_term_descriptor(:position_y)
    x2_descriptor = CPBM._pqs_source_safe_term_descriptor(:x2_z)
    kinetic_descriptor = CPBM._pqs_source_safe_term_descriptor(:kinetic)
    unsupported_descriptor =
        CPBM._pqs_source_safe_term_descriptor(:unsupported_safe_term)

    @test overlap_descriptor.status == :available_pqs_source_safe_term
    @test overlap_descriptor.source_term == :source_overlap
    @test overlap_descriptor.family == :overlap
    @test overlap_descriptor.required_factor_name == "overlap_1d"
    @test overlap_descriptor.batch_materialization_path ==
          :ready_pqs_source_overlap_blocks_only
    @test position_descriptor.status == :available_pqs_source_safe_term
    @test position_descriptor.source_term == :source_position_y
    @test position_descriptor.family == :position
    @test position_descriptor.axis == :y
    @test position_descriptor.required_factor_name == "position_1d"
    @test position_descriptor.unsupported_record_blocker ==
          :unsupported_pqs_source_position_materialization_record
    @test x2_descriptor.status == :available_pqs_source_safe_term
    @test x2_descriptor.source_term == :source_x2_z
    @test x2_descriptor.family == :x2
    @test x2_descriptor.axis == :z
    @test x2_descriptor.required_factor_name == "x2_1d"
    @test kinetic_descriptor.status == :available_pqs_source_safe_term
    @test kinetic_descriptor.source_term == :source_kinetic
    @test kinetic_descriptor.family == :kinetic
    @test kinetic_descriptor.required_factor_name == "kinetic_1d"
    @test unsupported_descriptor.status ==
          :blocked_unsupported_pqs_source_safe_term
    @test unsupported_descriptor.blocker == :unsupported_pqs_source_one_body_term
    @test unsupported_descriptor.source_term === nothing

    source_dims = pqs_record.metadata.left_source_mode_dims
    selector_overlap_x, selector_overlap_y, selector_overlap_z =
        _pair_block_pqs_source_overlap_axes(source_dims, source_dims)
    selector_position_x, selector_position_y, selector_position_z =
        _pair_block_pqs_source_position_axes(source_dims, source_dims)
    selector_x2_x, selector_x2_y, selector_x2_z =
        _pair_block_pqs_source_x2_axes(source_dims, source_dims)
    selector_kinetic_x, selector_kinetic_y, selector_kinetic_z =
        _pair_block_pqs_source_kinetic_axes(source_dims, source_dims)
    selector_overlap_1d =
        (; x = selector_overlap_x, y = selector_overlap_y, z = selector_overlap_z)

    record_overlap = CPBM.pqs_source_pair_overlap_block(
        pqs_record;
        overlap_1d = selector_overlap_1d,
    )
    batch_overlap = CPBM.pqs_source_pair_overlap_blocks(
        materialization_plan;
        overlap_1d = selector_overlap_1d,
    )
    batch_overlap_record = _pair_block_batch_result_for(
        batch_overlap,
        :pair_block_selector_pqs_unit,
        :pair_block_selector_pqs_unit,
    )

    @test batch_overlap.term == :source_overlap
    @test batch_overlap.materialized_count == 1
    @test batch_overlap.skipped_count == 2
    @test batch_overlap.materialized
    @test batch_overlap.source_operator_blocks_materialized
    @test !batch_overlap.final_pair_blocks_materialized
    @test !batch_overlap.operator_blocks_materialized
    @test !batch_overlap.hamiltonian_data_materialized
    @test !batch_overlap.artifacts_materialized
    @test batch_overlap.metadata.materialization_path ==
          :ready_pqs_source_overlap_blocks_only
    @test batch_overlap.metadata.pair_block_record_count == 3
    @test batch_overlap_record.block ≈ record_overlap.block
    @test :unsupported_pqs_source_overlap_materialization_record in
          Tuple(skipped.blocker for skipped in batch_overlap.skipped_records)

    record_position = CPBM.pqs_source_pair_position_block(
        pqs_record;
        axis = :y,
        overlap_1d = selector_overlap_1d,
        position_1d = (
            selector_position_x,
            selector_position_y,
            selector_position_z,
        ),
    )
    batch_position = CPBM.pqs_source_pair_position_blocks(
        materialization_plan;
        axis = :y,
        overlap_1d = selector_overlap_1d,
        position_1d = (;
            x = selector_position_x,
            y = selector_position_y,
            z = selector_position_z,
        ),
    )
    batch_position_record = _pair_block_batch_result_for(
        batch_position,
        :pair_block_selector_pqs_unit,
        :pair_block_selector_pqs_unit,
    )

    @test batch_position.term == :source_position_y
    @test batch_position.materialized_count == 1
    @test batch_position.skipped_count == 2
    @test batch_position.materialized
    @test batch_position.source_operator_blocks_materialized
    @test !batch_position.final_pair_blocks_materialized
    @test !batch_position.operator_blocks_materialized
    @test !batch_position.hamiltonian_data_materialized
    @test !batch_position.artifacts_materialized
    @test batch_position.metadata.materialization_path ==
          :ready_pqs_source_position_blocks_only
    @test batch_position.metadata.position_axis == :y
    @test batch_position_record.block ≈ record_position.block

    record_x2 = CPBM.pqs_source_pair_x2_block(
        pqs_record;
        axis = :z,
        overlap_1d = selector_overlap_1d,
        x2_1d = (; x = selector_x2_x, y = selector_x2_y, z = selector_x2_z),
    )
    batch_x2 = CPBM.pqs_source_pair_x2_blocks(
        materialization_plan;
        axis = :z,
        overlap_1d = selector_overlap_1d,
        x2_1d = (selector_x2_x, selector_x2_y, selector_x2_z),
    )
    batch_x2_record = _pair_block_batch_result_for(
        batch_x2,
        :pair_block_selector_pqs_unit,
        :pair_block_selector_pqs_unit,
    )

    @test batch_x2.term == :source_x2_z
    @test batch_x2.materialized_count == 1
    @test batch_x2.skipped_count == 2
    @test batch_x2.materialized
    @test batch_x2.source_operator_blocks_materialized
    @test !batch_x2.final_pair_blocks_materialized
    @test !batch_x2.operator_blocks_materialized
    @test !batch_x2.hamiltonian_data_materialized
    @test !batch_x2.artifacts_materialized
    @test batch_x2.metadata.materialization_path ==
          :ready_pqs_source_x2_blocks_only
    @test batch_x2.metadata.x2_axis == :z
    @test batch_x2_record.block ≈ record_x2.block

    record_kinetic = CPBM.pqs_source_pair_kinetic_block(
        pqs_record;
        overlap_1d = selector_overlap_1d,
        kinetic_1d = (;
            x = selector_kinetic_x,
            y = selector_kinetic_y,
            z = selector_kinetic_z,
        ),
    )
    batch_kinetic = CPBM.pqs_source_pair_kinetic_blocks(
        materialization_plan;
        overlap_1d = selector_overlap_1d,
        kinetic_1d = (
            selector_kinetic_x,
            selector_kinetic_y,
            selector_kinetic_z,
        ),
    )
    batch_kinetic_record = _pair_block_batch_result_for(
        batch_kinetic,
        :pair_block_selector_pqs_unit,
        :pair_block_selector_pqs_unit,
    )

    @test batch_kinetic.term == :source_kinetic
    @test batch_kinetic.materialized_count == 1
    @test batch_kinetic.skipped_count == 2
    @test batch_kinetic.materialized
    @test batch_kinetic.source_operator_blocks_materialized
    @test !batch_kinetic.final_pair_blocks_materialized
    @test !batch_kinetic.operator_blocks_materialized
    @test !batch_kinetic.hamiltonian_data_materialized
    @test !batch_kinetic.artifacts_materialized
    @test batch_kinetic.metadata.materialization_path ==
          :ready_pqs_source_kinetic_blocks_only
    @test batch_kinetic.metadata.kinetic_factor_form == (
        (:kinetic, :overlap, :overlap),
        (:overlap, :kinetic, :overlap),
        (:overlap, :overlap, :kinetic),
    )
    @test batch_kinetic_record.block ≈ record_kinetic.block

    selector_record_overlap = CPBM.pqs_source_pair_one_body_block(
        pqs_record,
        :overlap;
        overlap_1d = selector_overlap_1d,
    )
    selector_record_position = CPBM.pqs_source_pair_one_body_block(
        pqs_record,
        :position_y;
        overlap_1d = selector_overlap_1d,
        position_1d = (;
            x = selector_position_x,
            y = selector_position_y,
            z = selector_position_z,
        ),
    )
    selector_record_x2 = CPBM.pqs_source_pair_one_body_block(
        pqs_record,
        :x2_z;
        overlap_1d = selector_overlap_1d,
        x2_1d = (selector_x2_x, selector_x2_y, selector_x2_z),
    )
    selector_record_kinetic = CPBM.pqs_source_pair_one_body_block(
        pqs_record,
        :kinetic;
        overlap_1d = selector_overlap_1d,
        kinetic_1d = (;
            x = selector_kinetic_x,
            y = selector_kinetic_y,
            z = selector_kinetic_z,
        ),
    )

    @test selector_record_overlap.term == :source_overlap
    @test selector_record_overlap.block ≈ record_overlap.block
    @test selector_record_position.term == :source_position_y
    @test selector_record_position.block ≈ record_position.block
    @test selector_record_x2.term == :source_x2_z
    @test selector_record_x2.block ≈ record_x2.block
    @test selector_record_kinetic.term == :source_kinetic
    @test selector_record_kinetic.block ≈ record_kinetic.block
    @test selector_record_kinetic.source_operator_blocks_materialized
    @test !selector_record_kinetic.final_pair_blocks_materialized
    @test !selector_record_kinetic.operator_blocks_materialized
    @test !selector_record_kinetic.hamiltonian_data_materialized
    @test !selector_record_kinetic.artifacts_materialized

    retained_selector_record_overlap =
        CPBM.pqs_source_pair_retained_one_body_block(
            pqs_record,
            :overlap;
            overlap_1d = selector_overlap_1d,
        )
    retained_helper_record_overlap =
        CPBM.pqs_source_pair_retained_overlap_block(
            pqs_record;
            overlap_1d = selector_overlap_1d,
        )
    retained_selector_record_kinetic =
        CPBM.pqs_source_pair_retained_one_body_block(
            pqs_record,
            :kinetic;
            overlap_1d = selector_overlap_1d,
            kinetic_1d = (;
                x = selector_kinetic_x,
                y = selector_kinetic_y,
                z = selector_kinetic_z,
            ),
        )
    retained_helper_record_kinetic =
        CPBM.pqs_source_pair_retained_kinetic_block(
            pqs_record;
            overlap_1d = selector_overlap_1d,
            kinetic_1d = (
                selector_kinetic_x,
                selector_kinetic_y,
                selector_kinetic_z,
            ),
        )

    @test retained_selector_record_overlap.block ≈
          retained_helper_record_overlap.block
    @test retained_selector_record_kinetic.block ≈
          retained_helper_record_kinetic.block
    @test retained_selector_record_overlap.metadata.block_space ==
          :retained_pqs_source_modes
    @test retained_selector_record_kinetic.metadata.block_space ==
          :retained_pqs_source_modes
    @test !retained_selector_record_overlap.metadata.shell_realization_materialized
    @test !retained_selector_record_overlap.metadata.lowdin_cleanup_used

    selector_batch_overlap = CPBM.pqs_source_pair_one_body_blocks(
        materialization_plan,
        :overlap;
        overlap_1d = selector_overlap_1d,
    )
    selector_batch_position = CPBM.pqs_source_pair_one_body_blocks(
        materialization_plan,
        :position_y;
        overlap_1d = selector_overlap_1d,
        position_1d = (
            selector_position_x,
            selector_position_y,
            selector_position_z,
        ),
    )
    selector_batch_x2 = CPBM.pqs_source_pair_one_body_blocks(
        materialization_plan,
        :x2_z;
        overlap_1d = selector_overlap_1d,
        x2_1d = (; x = selector_x2_x, y = selector_x2_y, z = selector_x2_z),
    )
    selector_batch_kinetic = CPBM.pqs_source_pair_one_body_blocks(
        materialization_plan,
        :kinetic;
        overlap_1d = selector_overlap_1d,
        kinetic_1d = (
            selector_kinetic_x,
            selector_kinetic_y,
            selector_kinetic_z,
        ),
    )

    @test selector_batch_overlap.term == batch_overlap.term
    @test _pair_block_batch_result_for(
        selector_batch_overlap,
        :pair_block_selector_pqs_unit,
        :pair_block_selector_pqs_unit,
    ).block ≈ batch_overlap_record.block
    @test selector_batch_position.term == batch_position.term
    @test _pair_block_batch_result_for(
        selector_batch_position,
        :pair_block_selector_pqs_unit,
        :pair_block_selector_pqs_unit,
    ).block ≈ batch_position_record.block
    @test selector_batch_x2.term == batch_x2.term
    @test _pair_block_batch_result_for(
        selector_batch_x2,
        :pair_block_selector_pqs_unit,
        :pair_block_selector_pqs_unit,
    ).block ≈ batch_x2_record.block
    @test selector_batch_kinetic.term == batch_kinetic.term
    @test _pair_block_batch_result_for(
        selector_batch_kinetic,
        :pair_block_selector_pqs_unit,
        :pair_block_selector_pqs_unit,
    ).block ≈ batch_kinetic_record.block
    @test selector_batch_kinetic.source_operator_blocks_materialized
    @test !selector_batch_kinetic.final_pair_blocks_materialized
    @test !selector_batch_kinetic.operator_blocks_materialized
    @test !selector_batch_kinetic.hamiltonian_data_materialized
    @test !selector_batch_kinetic.artifacts_materialized

    retained_selector_batch_overlap =
        CPBM.pqs_source_pair_retained_one_body_blocks(
            materialization_plan,
            :overlap;
            overlap_1d = selector_overlap_1d,
        )
    retained_selector_batch_kinetic =
        CPBM.pqs_source_pair_retained_one_body_blocks(
            materialization_plan,
            :kinetic;
            overlap_1d = selector_overlap_1d,
            kinetic_1d = (
                selector_kinetic_x,
                selector_kinetic_y,
                selector_kinetic_z,
            ),
        )
    retained_batch_overlap_record = _pair_block_batch_result_for(
        retained_selector_batch_overlap,
        :pair_block_selector_pqs_unit,
        :pair_block_selector_pqs_unit,
    )
    retained_batch_kinetic_record = _pair_block_batch_result_for(
        retained_selector_batch_kinetic,
        :pair_block_selector_pqs_unit,
        :pair_block_selector_pqs_unit,
    )

    @test retained_selector_batch_overlap.materialized_count ==
          selector_batch_overlap.materialized_count
    @test retained_selector_batch_overlap.skipped_count ==
          selector_batch_overlap.skipped_count
    @test retained_selector_batch_kinetic.materialized_count ==
          selector_batch_kinetic.materialized_count
    @test retained_batch_overlap_record.block ≈ retained_selector_record_overlap.block
    @test retained_batch_kinetic_record.block ≈ retained_selector_record_kinetic.block
    @test retained_batch_overlap_record.metadata.block_space ==
          :retained_pqs_source_modes
    @test retained_batch_kinetic_record.metadata.block_space ==
          :retained_pqs_source_modes
    @test !retained_batch_overlap_record.metadata.shell_realization_materialized
    @test !retained_batch_overlap_record.metadata.lowdin_cleanup_used

    symmetric_overlap_axis =
        [i == j ? 1.0 : 0.0 for i in 1:source_dims[1], j in 1:source_dims[1]]
    symmetric_overlap_1d =
        (; x = symmetric_overlap_axis, y = symmetric_overlap_axis, z = symmetric_overlap_axis)
    symmetric_kinetic_axis =
        [i == j ? -Float64(i) : 0.01 * (i + j) for i in 1:source_dims[1], j in 1:source_dims[1]]
    symmetric_retained_batch_overlap =
        CPBM.pqs_source_pair_retained_one_body_blocks(
            materialization_plan,
            :overlap;
            overlap_1d = symmetric_overlap_1d,
        )
    symmetric_retained_batch_kinetic =
        CPBM.pqs_source_pair_retained_one_body_blocks(
            materialization_plan,
            :kinetic;
            overlap_1d = symmetric_overlap_1d,
            kinetic_1d = (;
                x = symmetric_kinetic_axis,
                y = symmetric_kinetic_axis,
                z = symmetric_kinetic_axis,
            ),
        )
    symmetric_retained_batch_overlap_record = _pair_block_batch_result_for(
        symmetric_retained_batch_overlap,
        :pair_block_selector_pqs_unit,
        :pair_block_selector_pqs_unit,
    )
    symmetric_retained_batch_kinetic_record = _pair_block_batch_result_for(
        symmetric_retained_batch_kinetic,
        :pair_block_selector_pqs_unit,
        :pair_block_selector_pqs_unit,
    )
    retained_overlap_matrix =
        CPBM.pqs_retained_source_one_body_matrix(symmetric_retained_batch_overlap)
    retained_kinetic_matrix =
        CPBM.pqs_retained_source_one_body_matrix(symmetric_retained_batch_kinetic)

    @test retained_overlap_matrix.status ==
          :materialized_pqs_retained_source_one_body_matrix
    @test retained_kinetic_matrix.status ==
          :materialized_pqs_retained_source_one_body_matrix
    @test retained_overlap_matrix.matrix ===
          symmetric_retained_batch_overlap_record.block
    @test retained_kinetic_matrix.matrix ===
          symmetric_retained_batch_kinetic_record.block
    @test retained_overlap_matrix.retained_dimension ==
          size(symmetric_retained_batch_overlap_record.block, 1)
    @test retained_kinetic_matrix.retained_dimension ==
          size(symmetric_retained_batch_kinetic_record.block, 1)
    @test retained_overlap_matrix.matrix_space == :retained_pqs_source_modes
    @test retained_kinetic_matrix.matrix_space == :retained_pqs_source_modes
    @test all(isfinite, retained_overlap_matrix.matrix)
    @test all(isfinite, retained_kinetic_matrix.matrix)
    @test retained_overlap_matrix.matrix ≈ transpose(retained_overlap_matrix.matrix)
    @test retained_kinetic_matrix.matrix ≈ transpose(retained_kinetic_matrix.matrix)
    @test !retained_overlap_matrix.shell_realization_materialized
    @test !retained_overlap_matrix.lowdin_cleanup_used
    @test !retained_kinetic_matrix.shell_realization_materialized
    @test !retained_kinetic_matrix.lowdin_cleanup_used

    CFBRForPairBlocks = GaussletBases.CartesianFinalBasisRealization
    @test CPBM.pqs_source_shell_realization_final_basis ===
          CFBRForPairBlocks.pqs_source_shell_realization_final_basis
    @test CPBM.pqs_source_shell_projected_one_body_matrix ===
          CFBRForPairBlocks.pqs_source_shell_projected_one_body_matrix
    @test CPBM.pqs_source_shell_final_one_body_from_boundary_matrix ===
          CFBRForPairBlocks.pqs_source_shell_final_one_body_from_boundary_matrix

    final_basis_source_box = CPBForPairBlocks.filled_cpb(
        1:5,
        1:5,
        1:5;
        role = :pqs_final_basis_contract_source,
    )
    final_basis_raw_plan = CRPSForPairBlocks.raw_product_box_plan(
        final_basis_source_box;
        source_mode_dims = (5, 5, 5),
        source_key = :pqs_final_basis_contract_source,
    )
    final_basis_retained_rule =
        CRPSForPairBlocks.pqs_boundary_product_mode_retained_rule(
            final_basis_raw_plan,
        )
    final_basis_boundary_count = final_basis_retained_rule.retained_count
    final_basis_identity = Matrix{Float64}(
        I,
        final_basis_boundary_count,
        final_basis_boundary_count,
    )
    final_basis = CPBM.pqs_source_shell_realization_final_basis(
        final_basis_raw_plan,
        final_basis_retained_rule;
        shell_support_indices = collect(1:final_basis_boundary_count),
        shell_overlap = final_basis_identity,
        shell_projection = final_basis_identity,
        lowdin_cleanup = final_basis_identity,
    )
    shell_operator = copy(final_basis_identity)
    shell_operator[1, 1] = 2.0
    shell_operator[2, 2] = 3.0
    shell_operator[1, 2] = 0.25
    shell_operator[2, 1] = 0.25
    boundary_kinetic_result = (;
        object_kind = :pqs_retained_source_one_body_matrix,
        status = :materialized_pqs_retained_source_one_body_matrix,
        blocker = nothing,
        term = :retained_source_kinetic,
        matrix = shell_operator,
        matrix_space = :retained_pqs_source_modes,
        retained_dimension = final_basis_boundary_count,
        matrix_materialized = true,
    )
    final_kinetic_from_boundary =
        CPBM.pqs_source_shell_final_one_body_from_boundary_matrix(
            final_basis,
            boundary_kinetic_result;
            term = :kinetic,
        )

    nuclear_boundary = copy(final_basis_identity)
    nuclear_boundary[1, 1] = -2.0
    nuclear_boundary[2, 2] = -3.0
    nuclear_boundary[1, 2] = -0.125
    nuclear_boundary[2, 1] = -0.125
    retained_nuclear_boundary = CPBM.PairBlockMaterializationResult(
        :retained_source_electron_nuclear_by_center,
        (:pqs_source_shell, :pqs_source_shell),
        nuclear_boundary,
        true,
        true,
        false,
        false,
        false,
        false,
        (;
            block_space = :retained_pqs_source_modes,
            by_center = true,
            center_key = :synthetic_origin,
            center_index = 1,
            center_location = (0.0, 0.0, 0.0),
            nuclear_charge = 2.0,
            nuclear_charge_recorded = true,
            nuclear_charge_applied = false,
            centers_summed = false,
            uncharged_by_center_convention = true,
            retained_source_operator_block_materialized = true,
            shell_realization_materialized = false,
            lowdin_cleanup_used = false,
        ),
    )
    final_nuclear_from_boundary =
        CPBM.pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block(
            final_basis,
            retained_nuclear_boundary,
        )
    @test final_nuclear_from_boundary.object_kind ==
          :pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block
    @test final_nuclear_from_boundary.status ==
          :materialized_pqs_shell_final_electron_nuclear_by_center_from_boundary_block
    @test final_nuclear_from_boundary.blocker === nothing
    @test final_nuclear_from_boundary.term == :electron_nuclear_by_center
    @test final_nuclear_from_boundary.input_term ==
          :retained_source_electron_nuclear_by_center
    @test final_nuclear_from_boundary.input_block_space ==
          :retained_pqs_source_modes
    @test final_nuclear_from_boundary.boundary_operator == nuclear_boundary
    @test final_nuclear_from_boundary.final_operator == nuclear_boundary
    @test final_nuclear_from_boundary.center_key == :synthetic_origin
    @test final_nuclear_from_boundary.center_index == 1
    @test final_nuclear_from_boundary.center_location == (0.0, 0.0, 0.0)
    @test final_nuclear_from_boundary.nuclear_charge == 2.0
    @test final_nuclear_from_boundary.nuclear_charge_recorded
    @test !final_nuclear_from_boundary.nuclear_charge_applied
    @test !final_nuclear_from_boundary.centers_summed
    @test final_nuclear_from_boundary.uncharged_by_center_convention
    @test final_nuclear_from_boundary.retained_boundary_operator_input_used
    @test !final_nuclear_from_boundary.raw_source_operator_input_used
    @test !final_nuclear_from_boundary.shell_support_operator_generated
    @test final_nuclear_from_boundary.one_body_operator_materialized
    @test final_nuclear_from_boundary.electron_nuclear_materialized
    @test !final_nuclear_from_boundary.charge_summing_materialized
    @test !final_nuclear_from_boundary.h1_solve_materialized
    @test !final_nuclear_from_boundary.hamiltonian_data_materialized
    @test !final_nuclear_from_boundary.ida_data_materialized
    @test !final_nuclear_from_boundary.density_density_materialized
    @test !final_nuclear_from_boundary.rhf_materialized
    @test !final_nuclear_from_boundary.driver_route_materialized
    @test !final_nuclear_from_boundary.exports_materialized
    @test !final_nuclear_from_boundary.artifacts_materialized
    @test final_nuclear_from_boundary.next_blocker ==
          :missing_pqs_final_one_electron_hamiltonian_assembly
    @test !final_nuclear_from_boundary.metadata.current_route_safe_term_matrices_used
    @test !final_nuclear_from_boundary.metadata.old_fixed_block_matrix_authority_used

    second_nuclear_boundary = copy(final_basis_identity)
    second_nuclear_boundary[1, 1] = -0.5
    second_nuclear_boundary[3, 3] = -0.25
    retained_second_nuclear_boundary = CPBM.PairBlockMaterializationResult(
        :retained_source_electron_nuclear_by_center,
        (:pqs_source_shell, :pqs_source_shell),
        second_nuclear_boundary,
        true,
        true,
        false,
        false,
        false,
        false,
        merge(
            retained_nuclear_boundary.metadata,
            (;
                center_key = :synthetic_shifted,
                center_index = 2,
                center_location = (0.25, 0.0, -0.25),
                nuclear_charge = 3.0,
            ),
        ),
    )
    final_second_nuclear_from_boundary =
        CPBM.pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block(
            final_basis,
            retained_second_nuclear_boundary,
        )
    final_hamiltonian =
        CPBM.pqs_source_shell_final_one_electron_hamiltonian(
            final_kinetic_from_boundary,
            (final_nuclear_from_boundary, final_second_nuclear_from_boundary),
        )
    expected_charged_nuclear =
        2.0 .* nuclear_boundary .+ 3.0 .* second_nuclear_boundary
    expected_hamiltonian = shell_operator .+ expected_charged_nuclear
    @test final_hamiltonian.object_kind ==
          :pqs_source_shell_final_one_electron_hamiltonian
    @test final_hamiltonian.status ==
          :materialized_pqs_shell_final_one_electron_hamiltonian
    @test final_hamiltonian.blocker === nothing
    @test final_hamiltonian.final_dimension == final_basis.final_retained_count
    @test final_hamiltonian.kinetic_matrix == shell_operator
    @test final_hamiltonian.charged_nuclear_matrix == expected_charged_nuclear
    @test final_hamiltonian.hamiltonian_matrix == expected_hamiltonian
    @test final_hamiltonian.center_count == 2
    @test final_hamiltonian.applied_nuclear_charges == (2.0, 3.0)
    @test final_hamiltonian.centers_summed
    @test final_hamiltonian.nuclear_charges_applied
    @test !final_hamiltonian.input_nuclear_charges_applied
    @test !final_hamiltonian.input_centers_summed
    @test final_hamiltonian.hamiltonian_data_materialized
    @test !final_hamiltonian.h1_solve_materialized
    @test !final_hamiltonian.eigensolve_materialized
    @test !final_hamiltonian.generalized_overlap_solve_materialized
    @test !final_hamiltonian.ida_data_materialized
    @test !final_hamiltonian.density_density_materialized
    @test !final_hamiltonian.rhf_materialized
    @test !final_hamiltonian.driver_route_materialized
    @test !final_hamiltonian.exports_materialized
    @test !final_hamiltonian.artifacts_materialized
    @test final_hamiltonian.next_blocker ==
          :missing_pqs_final_one_electron_solve
    @test !final_hamiltonian.metadata.current_route_safe_term_matrices_used
    @test !final_hamiltonian.metadata.old_fixed_block_matrix_authority_used

    overlap_bridge_batch =
        CPBM.pqs_source_pair_shell_realization_bridge_summary(selector_batch_overlap)
    @test overlap_bridge_batch.object_kind ==
          :pqs_source_pair_shell_realization_bridge_batch_summary
    @test overlap_bridge_batch.status ==
          :available_metadata_only_shell_realization_bridge_batch
    @test isnothing(overlap_bridge_batch.blocker)
    @test overlap_bridge_batch.term == :source_overlap
    @test _pair_block_count(overlap_bridge_batch.term_counts, :term, :source_overlap) == 1
    @test overlap_bridge_batch.result_count == 1
    @test overlap_bridge_batch.available_count == 1
    @test overlap_bridge_batch.blocked_count == 0
    @test overlap_bridge_batch.skipped_record_count == 2
    @test _pair_block_count(
        overlap_bridge_batch.skipped_blocker_counts,
        :blocker,
        :unsupported_pqs_source_overlap_materialization_record,
    ) == 1
    @test overlap_bridge_batch.source_mode_ordering_status ==
          :uniform_source_mode_ordering
    @test overlap_bridge_batch.source_mode_ordering == :x_major_y_major_z_fast
    @test overlap_bridge_batch.source_operator_blocks_materialized
    @test !overlap_bridge_batch.final_pair_blocks_materialized
    @test !overlap_bridge_batch.shell_realization_materialized
    @test !overlap_bridge_batch.operator_blocks_materialized
    @test !overlap_bridge_batch.hamiltonian_data_materialized
    @test !overlap_bridge_batch.artifacts_materialized

    overlap_final_readiness =
        CPBM.pqs_source_pair_final_block_readiness_summary(overlap_bridge_batch)
    @test overlap_final_readiness.object_kind ==
          :pqs_source_pair_final_block_readiness_batch_summary
    @test overlap_final_readiness.status ==
          :blocked_final_pqs_pair_block_not_ready
    @test overlap_final_readiness.blocker == :shell_realization_not_materialized
    @test overlap_final_readiness.bridge_status ==
          :available_metadata_only_shell_realization_bridge_batch
    @test overlap_final_readiness.term == :source_overlap
    @test overlap_final_readiness.result_count == 1
    @test overlap_final_readiness.available_count == 1
    @test overlap_final_readiness.blocked_count == 0
    @test overlap_final_readiness.source_operator_blocks_materialized
    @test !overlap_final_readiness.shell_realization_materialized
    @test !overlap_final_readiness.final_pair_blocks_materialized
    @test !overlap_final_readiness.operator_blocks_materialized

    empty_bridge_batch =
        merge(
            overlap_bridge_batch,
            (;
                status = :available_metadata_only_shell_realization_bridge_batch,
                blocker = nothing,
                term_counts = (),
                result_count = 0,
                available_count = 0,
                blocked_count = 0,
                bridge_status_counts = (),
                blocker_counts = (),
            ),
        )
    empty_batch_readiness =
        CPBM.pqs_source_pair_final_block_readiness_summary(empty_bridge_batch)
    @test empty_batch_readiness.status ==
          :blocked_final_pqs_pair_block_not_ready
    @test empty_batch_readiness.blocker == :no_pqs_source_block_results
    @test empty_batch_readiness.result_count == 0
    @test empty_batch_readiness.available_count == 0
    @test empty_batch_readiness.blocked_count == 0
    @test !empty_batch_readiness.final_pair_blocks_materialized
    @test !empty_batch_readiness.shell_realization_materialized

    bridge_tuple =
        CPBM.pqs_source_pair_shell_realization_bridge_summary((
            selector_record_overlap,
            selector_record_position,
            selector_record_kinetic,
        ))
    @test bridge_tuple.status ==
          :available_metadata_only_shell_realization_bridge_batch
    @test bridge_tuple.term == :mixed_source_terms
    @test bridge_tuple.result_count == 3
    @test bridge_tuple.available_count == 3
    @test bridge_tuple.blocked_count == 0
    @test _pair_block_count(bridge_tuple.term_counts, :term, :source_overlap) == 1
    @test _pair_block_count(bridge_tuple.term_counts, :term, :source_position_y) == 1
    @test _pair_block_count(bridge_tuple.term_counts, :term, :source_kinetic) == 1
    @test bridge_tuple.source_mode_ordering_status == :uniform_source_mode_ordering
    @test bridge_tuple.source_operator_blocks_materialized
    @test !bridge_tuple.final_pair_blocks_materialized
    @test !bridge_tuple.shell_realization_materialized

    bridge_vector =
        CPBM.pqs_source_pair_shell_realization_bridge_summary([
            selector_record_overlap,
            selector_record_position,
        ])
    @test bridge_vector.term == :mixed_source_terms
    @test bridge_vector.result_count == 2
    @test bridge_vector.available_count == 2
    @test bridge_vector.blocked_count == 0

    non_pqs_source_result =
        _pair_block_result_with_metadata(
            selector_record_overlap,
            merge(
                selector_record_overlap.metadata,
                (; block_space = :direct_source_cpb_modes),
            ),
        )
    missing_source_metadata_result =
        _pair_block_result_with_metadata(
            selector_record_overlap,
            merge(selector_record_overlap.metadata, (; left_source_mode_dims = nothing)),
        )
    missing_transform_keys_result =
        _pair_block_result_with_metadata(
            selector_record_overlap,
            merge(
                selector_record_overlap.metadata,
                (; transform_contract_keys = nothing),
            ),
        )
    missing_realization_paths_result =
        _pair_block_result_with_metadata(
            selector_record_overlap,
            merge(selector_record_overlap.metadata, (; realization_paths = nothing)),
        )
    blocked_bridge_batch =
        CPBM.pqs_source_pair_shell_realization_bridge_summary((
            non_pqs_source_result,
            missing_source_metadata_result,
            missing_transform_keys_result,
            missing_realization_paths_result,
        ))
    @test blocked_bridge_batch.status ==
          :blocked_pqs_source_shell_realization_bridge_batch
    @test blocked_bridge_batch.blocker == :not_raw_product_source_modes
    @test blocked_bridge_batch.result_count == 4
    @test blocked_bridge_batch.available_count == 0
    @test blocked_bridge_batch.blocked_count == 4
    @test _pair_block_count(
        blocked_bridge_batch.bridge_status_counts,
        :status,
        :blocked_not_pqs_source_space_block,
    ) == 1
    @test _pair_block_count(
        blocked_bridge_batch.bridge_status_counts,
        :status,
        :blocked_missing_source_block_metadata,
    ) == 1
    @test _pair_block_count(
        blocked_bridge_batch.bridge_status_counts,
        :status,
        :blocked_missing_shell_realization_facts,
    ) == 2
    @test _pair_block_count(
        blocked_bridge_batch.blocker_counts,
        :blocker,
        :not_raw_product_source_modes,
    ) == 1
    @test _pair_block_count(
        blocked_bridge_batch.blocker_counts,
        :blocker,
        :missing_left_source_mode_dims,
    ) == 1
    @test _pair_block_count(
        blocked_bridge_batch.blocker_counts,
        :blocker,
        :missing_transform_contract_keys,
    ) == 1
    @test _pair_block_count(
        blocked_bridge_batch.blocker_counts,
        :blocker,
        :missing_shell_realization_path,
    ) == 1
    @test blocked_bridge_batch.source_mode_ordering_status ==
          :uniform_source_mode_ordering
    @test blocked_bridge_batch.source_operator_blocks_materialized
    @test !blocked_bridge_batch.final_pair_blocks_materialized
    @test !blocked_bridge_batch.shell_realization_materialized
    blocked_final_readiness =
        CPBM.pqs_source_pair_final_block_readiness_summary(blocked_bridge_batch)
    @test blocked_final_readiness.status ==
          :blocked_final_pqs_pair_block_not_ready
    @test blocked_final_readiness.blocker == :not_raw_product_source_modes
    @test blocked_final_readiness.result_count == 4
    @test blocked_final_readiness.blocked_count == 4
    @test !blocked_final_readiness.final_pair_blocks_materialized
    @test !blocked_final_readiness.shell_realization_materialized

    @test_throws ArgumentError CPBM.pqs_source_pair_one_body_block(
        pqs_record,
        :position_x;
        overlap_1d = selector_overlap_1d,
    )
    @test_throws ArgumentError CPBM.pqs_source_pair_one_body_block(
        pqs_record,
        :x2_z;
        overlap_1d = selector_overlap_1d,
    )
    @test_throws ArgumentError CPBM.pqs_source_pair_one_body_block(
        pqs_record,
        :kinetic;
        overlap_1d = selector_overlap_1d,
    )
    @test_throws ArgumentError CPBM.pqs_source_pair_one_body_block(
        pqs_record,
        :unsupported_safe_term;
        overlap_1d = selector_overlap_1d,
    )
    @test_throws ArgumentError CPBM.pqs_source_pair_one_body_blocks(
        materialization_plan,
        :position_x;
        overlap_1d = selector_overlap_1d,
    )
    @test_throws ArgumentError CPBM.pqs_source_pair_one_body_blocks(
        materialization_plan,
        :unsupported_safe_term;
        overlap_1d = selector_overlap_1d,
    )
end

@testset "CartesianPairBlockMaterialization direct/direct overlap selector" begin
    retained_plan = _pair_block_rectangular_retained_plan()
    unit_pair_plan = CUPForPairBlocks.unit_pair_plan(retained_plan)
    transform_plan =
        CRTCForPairBlocks.retained_unit_transform_contract_plan(retained_plan)
    pair_operator_plan =
        CPOPForPairBlocks.pair_operator_plan(
            unit_pair_plan,
            transform_plan;
            route_core_sidecars = false,
        )
    materialization_plan =
        CPBM.pair_block_materialization_plan(pair_operator_plan)
    materialization_summary = CPBM.summary(materialization_plan)

    @test length(CPBM.pair_block_materialization_records(materialization_plan)) == 6
    @test materialization_summary.ready_record_count == 3
    @test materialization_summary.blocked_record_count == 3

    cross_record = _pair_block_record_for(
        materialization_plan,
        :pair_block_direct_left_unit,
        :pair_block_direct_right_unit,
    )
    overlap_x, overlap_y, overlap_z = _pair_block_overlap_axes()
    cross_result = CPBM.direct_direct_overlap_block(
        cross_record;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
    )
    left_source = only(cross_record.metadata.left_source_cpbs)
    right_source = only(cross_record.metadata.right_source_cpbs)
    expected_cross =
        _pair_block_expected_overlap(
            left_source,
            right_source,
            overlap_x,
            overlap_y,
            overlap_z,
        )

    @test size(cross_result.block) == (4, 6)
    @test cross_result.block ≈ expected_cross

    position_x, position_y, position_z = _pair_block_position_axes()
    expected_position_terms = (;
        x = :position_x,
        y = :position_y,
        z = :position_z,
    )
    for axis in (:x, :y, :z)
        position_result = CPBM.direct_direct_position_block(
            cross_record;
            axis,
            parent_axis_counts = (4, 4, 4),
            overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
            position_1d = (; x = position_x, y = position_y, z = position_z),
        )
        expected_position =
            _pair_block_expected_position(
                left_source,
                right_source,
                axis,
                overlap_x,
                overlap_y,
                overlap_z,
                position_x,
                position_y,
                position_z,
            )

        @test position_result.term == getproperty(expected_position_terms, axis)
        @test size(position_result.block) == (4, 6)
        @test position_result.block ≈ expected_position
        @test position_result.materialized
        @test position_result.source_operator_blocks_materialized
        @test position_result.final_pair_blocks_materialized
        @test !position_result.operator_blocks_materialized
        @test !position_result.hamiltonian_data_materialized
        @test !position_result.artifacts_materialized
    end

    x2_x, x2_y, x2_z = _pair_block_x2_axes()
    expected_x2_terms = (;
        x = :x2_x,
        y = :x2_y,
        z = :x2_z,
    )
    for axis in (:x, :y, :z)
        x2_result = CPBM.direct_direct_x2_block(
            cross_record;
            axis,
            parent_axis_counts = (4, 4, 4),
            overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
            x2_1d = (; x = x2_x, y = x2_y, z = x2_z),
        )
        expected_x2 =
            _pair_block_expected_x2(
                left_source,
                right_source,
                axis,
                overlap_x,
                overlap_y,
                overlap_z,
                x2_x,
                x2_y,
                x2_z,
            )

        @test x2_result.term == getproperty(expected_x2_terms, axis)
        @test size(x2_result.block) == (4, 6)
        @test x2_result.block ≈ expected_x2
        @test x2_result.materialized
        @test x2_result.source_operator_blocks_materialized
        @test x2_result.final_pair_blocks_materialized
        @test !x2_result.operator_blocks_materialized
        @test !x2_result.hamiltonian_data_materialized
        @test !x2_result.artifacts_materialized
    end

    kinetic_x, kinetic_y, kinetic_z = _pair_block_kinetic_axes()
    kinetic_result = CPBM.direct_direct_kinetic_block(
        cross_record;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        kinetic_1d = (; x = kinetic_x, y = kinetic_y, z = kinetic_z),
    )
    expected_kinetic =
        _pair_block_expected_kinetic(
            left_source,
            right_source,
            overlap_x,
            overlap_y,
            overlap_z,
            kinetic_x,
            kinetic_y,
            kinetic_z,
        )

    @test kinetic_result.term == :kinetic
    @test size(kinetic_result.block) == (4, 6)
    @test kinetic_result.block ≈ expected_kinetic
    @test kinetic_result.materialized
    @test kinetic_result.source_operator_blocks_materialized
    @test kinetic_result.final_pair_blocks_materialized
    @test !kinetic_result.operator_blocks_materialized
    @test !kinetic_result.hamiltonian_data_materialized
    @test !kinetic_result.artifacts_materialized

    selector_overlap = CPBM.direct_direct_one_body_block(
        cross_record,
        :overlap;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
    )
    position_y_result = CPBM.direct_direct_position_block(
        cross_record;
        axis = :y,
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        position_1d = (; x = position_x, y = position_y, z = position_z),
    )
    selector_position_y = CPBM.direct_direct_one_body_block(
        cross_record,
        :position_y;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        position_1d = (; x = position_x, y = position_y, z = position_z),
    )
    x2_z_result = CPBM.direct_direct_x2_block(
        cross_record;
        axis = :z,
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        x2_1d = (; x = x2_x, y = x2_y, z = x2_z),
    )
    selector_x2_z = CPBM.direct_direct_one_body_block(
        cross_record,
        :x2_z;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        x2_1d = (; x = x2_x, y = x2_y, z = x2_z),
    )
    selector_kinetic = CPBM.direct_direct_one_body_block(
        cross_record,
        :kinetic;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        kinetic_1d = (; x = kinetic_x, y = kinetic_y, z = kinetic_z),
    )

    @test selector_overlap.block ≈ cross_result.block
    @test selector_position_y.block ≈ position_y_result.block
    @test selector_x2_z.block ≈ x2_z_result.block
    @test selector_kinetic.block ≈ kinetic_result.block
    @test_throws ArgumentError CPBM.direct_direct_one_body_block(
        cross_record,
        :coulomb;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
    )
    @test_throws ArgumentError CPBM.direct_direct_one_body_block(
        cross_record,
        :position_x;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
    )
    @test_throws ArgumentError CPBM.direct_direct_one_body_block(
        cross_record,
        :x2_x;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
    )
    @test_throws ArgumentError CPBM.direct_direct_one_body_block(
        cross_record,
        :kinetic;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
    )

    batch_result = CPBM.direct_direct_overlap_blocks(
        materialization_plan;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
    )
    position_batch_result = CPBM.direct_direct_position_blocks(
        materialization_plan;
        axis = :z,
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        position_1d = (; x = position_x, y = position_y, z = position_z),
    )
    x2_batch_result = CPBM.direct_direct_x2_blocks(
        materialization_plan;
        axis = :z,
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        x2_1d = (; x = x2_x, y = x2_y, z = x2_z),
    )
    kinetic_batch_result = CPBM.direct_direct_kinetic_blocks(
        materialization_plan;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        kinetic_1d = (; x = kinetic_x, y = kinetic_y, z = kinetic_z),
    )
    selector_kinetic_batch = CPBM.direct_direct_one_body_blocks(
        materialization_plan,
        :kinetic;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        kinetic_1d = (; x = kinetic_x, y = kinetic_y, z = kinetic_z),
    )
    batch_cross = only(
        result for result in batch_result.materialized_results
        if result.pair_key == (:pair_block_direct_left_unit, :pair_block_direct_right_unit)
    )
    materialized_keys = Set(result.pair_key for result in batch_result.materialized_results)

    @test batch_result isa CPBM.PairBlockMaterializationBatchResult
    @test batch_result.term == :overlap
    @test batch_result.materialized_count == 3
    @test batch_result.skipped_count == 3
    @test materialized_keys == Set((
        (:pair_block_direct_left_unit, :pair_block_direct_left_unit),
        (:pair_block_direct_left_unit, :pair_block_direct_right_unit),
        (:pair_block_direct_right_unit, :pair_block_direct_right_unit),
    ))
    skipped_blockers = Tuple(skipped.blocker for skipped in batch_result.skipped_records)
    @test count(
        ==(:non_direct_direct_pair_block_materialization_not_implemented),
        skipped_blockers,
    ) == 2
    @test count(==(:missing_left_raw_product_source_plan), skipped_blockers) == 1
    @test batch_cross.block ≈ expected_cross
    @test batch_result.materialized
    @test batch_result.source_operator_blocks_materialized
    @test batch_result.final_pair_blocks_materialized
    @test !batch_result.operator_blocks_materialized
    @test !batch_result.hamiltonian_data_materialized
    @test !batch_result.artifacts_materialized

    @test position_batch_result.term == :position_z
    @test position_batch_result.materialized_count == 3
    @test position_batch_result.skipped_count == 3
    @test position_batch_result.materialized
    @test position_batch_result.source_operator_blocks_materialized
    @test position_batch_result.final_pair_blocks_materialized
    @test !position_batch_result.operator_blocks_materialized
    @test !position_batch_result.hamiltonian_data_materialized
    @test !position_batch_result.artifacts_materialized

    @test x2_batch_result.term == :x2_z
    @test x2_batch_result.materialized_count == 3
    @test x2_batch_result.skipped_count == 3
    @test x2_batch_result.materialized
    @test x2_batch_result.source_operator_blocks_materialized
    @test x2_batch_result.final_pair_blocks_materialized
    @test !x2_batch_result.operator_blocks_materialized
    @test !x2_batch_result.hamiltonian_data_materialized
    @test !x2_batch_result.artifacts_materialized

    @test kinetic_batch_result.term == :kinetic
    @test kinetic_batch_result.materialized_count == 3
    @test kinetic_batch_result.skipped_count == 3
    @test kinetic_batch_result.materialized
    @test kinetic_batch_result.source_operator_blocks_materialized
    @test kinetic_batch_result.final_pair_blocks_materialized
    @test !kinetic_batch_result.operator_blocks_materialized
    @test !kinetic_batch_result.hamiltonian_data_materialized
    @test !kinetic_batch_result.artifacts_materialized

    @test selector_kinetic_batch.materialized_count ==
          kinetic_batch_result.materialized_count
    @test selector_kinetic_batch.skipped_count == kinetic_batch_result.skipped_count
    @test selector_kinetic_batch.term == kinetic_batch_result.term
end
