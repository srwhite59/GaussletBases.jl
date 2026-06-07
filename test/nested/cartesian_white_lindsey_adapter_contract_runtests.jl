using Test
using GaussletBases

const CPBMForLWAdapter = GaussletBases.CartesianPairBlockMaterialization
const CRUForLWAdapter = GaussletBases.CartesianRetainedUnits
const CRCForLWAdapter = GaussletBases.CartesianRouteCore
const CUPForLWAdapter = GaussletBases.CartesianUnitPairs
const CPBForLWAdapter = GaussletBases.CartesianCPB

function _lw_adapter_retained_unit(
    unit_key::Symbol,
    unit_index::Int,
    source_cpb,
    stratum_kind::Symbol,
    source_cpb_index::Int;
    dimension_status::Symbol = :not_materialized,
    dimension = nothing,
    extra_metadata = (;),
)
    return CRUForLWAdapter.RetainedUnitRecord(
        unit_key,
        unit_index,
        :white_lindsey_boundary_stratum_retained_unit,
        Symbol(unit_key, "_contract"),
        Symbol(unit_key, "_terminal_region"),
        :synthetic_terminal_region,
        :synthetic_terminal_region,
        :white_lindsey_boundary_strata,
        :white_lindsey_boundary_stratum_product,
        :direct_or_trivial_embedding,
        CRCForLWAdapter.owned_cpb(source_cpb),
        (source_cpb,),
        source_cpb_index,
        dimension_status,
        dimension,
        :not_materialized,
        nothing,
        nothing,
        false,
        merge((; stratum_kind, source_cpb_index), extra_metadata),
    )
end

function _lw_adapter_unit_pair(left_unit, right_unit, pair_index::Int)
    return CUPForLWAdapter.UnitPairRecord(
        (left_unit.unit_key, right_unit.unit_key),
        pair_index,
        Symbol(String(left_unit.unit_kind), "__", String(right_unit.unit_kind)),
        left_unit,
        right_unit,
        left_unit.unit_index,
        right_unit.unit_index,
        left_unit.unit_key,
        right_unit.unit_key,
        left_unit.unit_kind,
        right_unit.unit_kind,
        nothing,
        false,
        (;),
    )
end

function _lw_adapter_descriptor_units()
    facet_source = CPBForLWAdapter.slab_cpb(
        1:1,
        1:3,
        1:3;
        role = :lw_adapter_test_facet_source_cpb,
        metadata = (; stratum_kind = :facet_cpb, source_cpb_index = 1),
    )
    edge_source = CPBForLWAdapter.cpb(
        4:4,
        2:2,
        1:3;
        role = :lw_adapter_test_edge_source_cpb,
        metadata = (; stratum_kind = :edge_cpb, source_cpb_index = 2),
    )
    corner_source = CPBForLWAdapter.cpb(
        4:4,
        3:3,
        3:3;
        role = :lw_adapter_test_corner_source_cpb,
        metadata = (; stratum_kind = :corner_cpb, source_cpb_index = 3),
    )
    return (
        _lw_adapter_retained_unit(
            :lw_adapter_test_facet_unit,
            1,
            facet_source,
            :facet_cpb,
            1,
        ),
        _lw_adapter_retained_unit(
            :lw_adapter_test_edge_unit,
            2,
            edge_source,
            :edge_cpb,
            2,
        ),
        _lw_adapter_retained_unit(
            :lw_adapter_test_corner_unit,
            3,
            corner_source,
            :corner_cpb,
            3,
        ),
    )
end

function _lw_adapter_doside_source_1d()
    count = 7
    endpoint = (count - 1) / 2
    a = 0.25
    xmax = 10.0
    tail_spacing = 10.0
    s = asinh(xmax / a) / (endpoint - xmax / tail_spacing)
    basis = build_basis(MappedUniformBasisSpec(:G10;
        count,
        mapping = AsinhMapping(; a, s, tail_spacing),
        reference_spacing = 1.0,
    ))
    expansion = coulomb_gaussian_expansion(doacc = false)
    return GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
end

@testset "CartesianPairBlockMaterialization White-Lindsey corner coefficients" begin
    facet_unit, edge_unit, corner_unit = _lw_adapter_descriptor_units()
    facet_descriptor =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            facet_unit,
        )
    edge_descriptor =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            edge_unit,
        )
    corner_descriptor =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            corner_unit,
        )
    @test facet_descriptor.source_cpb_intervals == (1:1, 1:3, 1:3)
    @test facet_descriptor.source_axis_intervals ==
          (x = 1:1, y = 1:3, z = 1:3)
    @test facet_descriptor.active_product_axis_intervals == (
        (; axis = :y, interval = 1:3),
        (; axis = :z, interval = 1:3),
    )
    @test isnothing(facet_descriptor.free_axis)
    @test isnothing(facet_descriptor.free_axis_interval)
    @test isnothing(facet_descriptor.fixed_side_metadata)

    @test edge_descriptor.source_cpb_intervals == (4:4, 2:2, 1:3)
    @test edge_descriptor.source_axis_intervals ==
          (x = 4:4, y = 2:2, z = 1:3)
    @test edge_descriptor.active_product_axis_intervals ==
          ((; axis = :z, interval = 1:3),)
    @test edge_descriptor.free_axis == :z
    @test edge_descriptor.free_axis_interval == 1:3
    @test isnothing(edge_descriptor.fixed_side_metadata)

    @test corner_descriptor.source_cpb_intervals == (4:4, 3:3, 3:3)
    @test corner_descriptor.source_axis_intervals ==
          (x = 4:4, y = 3:3, z = 3:3)
    @test corner_descriptor.active_product_axis_intervals == ()
    @test isnothing(corner_descriptor.free_axis)
    @test isnothing(corner_descriptor.free_axis_interval)

    corner_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficients(
            corner_descriptor,
        )
    @test corner_coefficients.object_kind ==
          :white_lindsey_boundary_stratum_unit_coefficients
    @test corner_coefficients.status ==
          :materialized_white_lindsey_corner_unit_coefficients
    @test isnothing(corner_coefficients.blocker)
    @test corner_coefficients.unit_key == :lw_adapter_test_corner_unit
    @test corner_coefficients.stratum_kind == :corner_cpb
    @test corner_coefficients.source_cpb_role ==
          :lw_adapter_test_corner_source_cpb
    @test corner_coefficients.source_cpb_shape == (1, 1, 1)
    @test corner_coefficients.fixed_axis_coordinates == (
        (; axis = :x, coordinate = 4),
        (; axis = :y, coordinate = 3),
        (; axis = :z, coordinate = 3),
    )
    @test corner_coefficients.planned_old_kernel == :_nested_corner_piece
    @test corner_coefficients.coefficient_space == :source_cpb_support_local
    @test corner_coefficients.source_cpb_intervals == (4:4, 3:3, 3:3)
    @test corner_coefficients.source_axis_intervals ==
          (x = 4:4, y = 3:3, z = 3:3)
    @test corner_coefficients.missing_coefficient_inputs == ()
    @test corner_coefficients.coefficient_input_requirements.status ==
          :available_corner_support_local_coefficients
    @test !corner_coefficients.parent_row_indices_available
    @test size(corner_coefficients.coefficient_matrix) == (1, 1)
    @test corner_coefficients.coefficient_matrix[1, 1] == 1.0
    @test corner_coefficients.source_support_row_count == 1
    @test corner_coefficients.retained_column_count == 1
    @test corner_coefficients.source_support_row_index == 1
    @test corner_coefficients.retained_column_index == 1
    @test corner_coefficients.nonzero_count == 1
    @test corner_coefficients.nonzero_values == (1.0,)
    @test corner_coefficients.coefficient_maps_materialized
    @test !corner_coefficients.source_operator_blocks_materialized
    @test !corner_coefficients.final_pair_blocks_materialized
    @test !corner_coefficients.operator_blocks_materialized
    @test !corner_coefficients.hamiltonian_data_materialized
    @test !corner_coefficients.artifacts_materialized

    corner_unit_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficients(
            corner_unit,
        )
    @test corner_unit_coefficients.status ==
          :materialized_white_lindsey_corner_unit_coefficients
    @test corner_unit_coefficients.coefficient_matrix == [1.0;;]

    facet_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficients(
            facet_descriptor,
        )
    @test facet_coefficients.status ==
          :blocked_white_lindsey_boundary_stratum_unit_coefficients
    @test facet_coefficients.blocker ==
          :incomplete_white_lindsey_edge_facet_kernel_input_context
    @test facet_coefficients.stratum_kind == :facet_cpb
    @test facet_coefficients.coefficient_input_requirements.status ==
          :blocked_missing_white_lindsey_facet_kernel_inputs
    @test facet_coefficients.coefficient_input_requirements.required_old_kernel ==
          :_nested_face_product
    @test facet_coefficients.coefficient_input_requirements.required_1d_helper ==
          :_nested_doside_1d
    @test facet_coefficients.missing_coefficient_inputs == (
        :missing_white_lindsey_doside_source_1d,
        :missing_white_lindsey_retained_count,
        :missing_white_lindsey_parent_dims,
        :missing_white_lindsey_fixed_side_metadata,
    )
    @test facet_coefficients.active_product_axis_intervals == (
        (; axis = :y, interval = 1:3),
        (; axis = :z, interval = 1:3),
    )
    @test isnothing(facet_coefficients.coefficient_matrix)
    @test !facet_coefficients.coefficient_maps_materialized

    edge_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficients(
            edge_descriptor,
        )
    @test edge_coefficients.status ==
          :blocked_white_lindsey_boundary_stratum_unit_coefficients
    @test edge_coefficients.blocker ==
          :incomplete_white_lindsey_edge_facet_kernel_input_context
    @test edge_coefficients.stratum_kind == :edge_cpb
    @test edge_coefficients.coefficient_input_requirements.status ==
          :blocked_missing_white_lindsey_edge_kernel_inputs
    @test edge_coefficients.coefficient_input_requirements.required_old_kernel ==
          :_nested_edge_product
    @test edge_coefficients.coefficient_input_requirements.required_1d_helper ==
          :_nested_doside_1d
    @test edge_coefficients.missing_coefficient_inputs == (
        :missing_white_lindsey_doside_source_1d,
        :missing_white_lindsey_retained_count,
        :missing_white_lindsey_parent_dims,
        :missing_white_lindsey_fixed_side_metadata,
    )
    @test edge_coefficients.free_axis == :z
    @test edge_coefficients.free_axis_interval == 1:3
    @test isnothing(edge_coefficients.coefficient_matrix)
    @test !edge_coefficients.coefficient_maps_materialized

    facet_context =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficient_context(
            facet_descriptor,
        )
    @test facet_context.object_kind ==
          :white_lindsey_boundary_stratum_unit_coefficient_context
    @test facet_context.status ==
          :blocked_missing_white_lindsey_facet_kernel_context
    @test facet_context.blocker ==
          :missing_white_lindsey_doside_source_1d
    @test facet_context.planned_old_calls ==
          (:_nested_doside_1d, :_nested_face_product)
    @test facet_context.face_kind == :yz
    @test facet_context.fixed_index == 1
    @test isnothing(facet_context.fixed_side)
    @test facet_context.missing_inputs == (
        :missing_white_lindsey_doside_source_1d,
        :missing_white_lindsey_retained_count,
        :missing_white_lindsey_parent_dims,
        :missing_white_lindsey_fixed_side_metadata,
    )
    @test !facet_context.coefficient_maps_materialized
    @test !facet_context.source_operator_blocks_materialized
    @test !facet_context.final_pair_blocks_materialized

    edge_context =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficient_context(
            edge_descriptor,
        )
    @test edge_context.status ==
          :blocked_missing_white_lindsey_edge_kernel_context
    @test edge_context.blocker ==
          :missing_white_lindsey_doside_source_1d
    @test edge_context.planned_old_calls ==
          (:_nested_doside_1d, :_nested_edge_product)
    @test edge_context.free_axis == :z
    @test edge_context.free_axis_interval == 1:3
    @test edge_context.fixed_indices == (4, 2)
    @test isnothing(edge_context.fixed_sides)
    @test edge_context.missing_inputs == (
        :missing_white_lindsey_doside_source_1d,
        :missing_white_lindsey_retained_count,
        :missing_white_lindsey_parent_dims,
        :missing_white_lindsey_fixed_side_metadata,
    )
    @test !edge_context.coefficient_maps_materialized

    ready_facet_descriptor = merge(
        facet_descriptor,
        (;
            fixed_side_metadata = ((; axis = :x, side = :low),),
            retained_count = 9,
            parent_dims = (7, 7, 7),
            doside_source_1d = :synthetic_doside_source_1d,
        ),
    )
    ready_facet_context =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficient_context(
            ready_facet_descriptor,
        )
    @test ready_facet_context.status ==
          :ready_white_lindsey_facet_kernel_context_not_materialized
    @test isnothing(ready_facet_context.blocker)
    @test ready_facet_context.face_kind == :yz
    @test ready_facet_context.fixed_side == :low
    @test ready_facet_context.fixed_index == 1
    @test ready_facet_context.active_product_axis_intervals == (
        (; axis = :y, interval = 1:3),
        (; axis = :z, interval = 1:3),
    )
    @test ready_facet_context.retained_count == 9
    @test ready_facet_context.parent_dims == (7, 7, 7)
    @test ready_facet_context.doside_source_1d == :synthetic_doside_source_1d
    @test ready_facet_context.planned_old_calls ==
          (:_nested_doside_1d, :_nested_face_product)
    @test ready_facet_context.missing_inputs == ()
    @test !ready_facet_context.coefficient_maps_materialized
    ready_facet_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficients(
            ready_facet_descriptor,
        )
    @test ready_facet_coefficients.status ==
          :blocked_white_lindsey_boundary_stratum_unit_coefficients
    @test ready_facet_coefficients.blocker ==
          :white_lindsey_facet_doside_source_not_materializable
    @test ready_facet_coefficients.coefficient_input_requirements.status ==
          :available_white_lindsey_facet_kernel_context_inputs
    @test ready_facet_coefficients.missing_coefficient_inputs == ()
    @test isnothing(ready_facet_coefficients.coefficient_matrix)
    @test !ready_facet_coefficients.coefficient_maps_materialized

    ready_edge_descriptor = merge(
        edge_descriptor,
        (;
            fixed_side_metadata = (
                (; axis = :x, side = :high),
                (; axis = :y, side = :low),
            ),
            retained_count = 3,
            parent_dims = (7, 7, 7),
            doside_source_1d = :synthetic_doside_source_1d,
        ),
    )
    ready_edge_context =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficient_context(
            ready_edge_descriptor,
        )
    @test ready_edge_context.status ==
          :ready_white_lindsey_edge_kernel_context_not_materialized
    @test isnothing(ready_edge_context.blocker)
    @test ready_edge_context.free_axis == :z
    @test ready_edge_context.free_axis_interval == 1:3
    @test ready_edge_context.fixed_sides == (:high, :low)
    @test ready_edge_context.fixed_indices == (4, 2)
    @test ready_edge_context.retained_count == 3
    @test ready_edge_context.parent_dims == (7, 7, 7)
    @test ready_edge_context.doside_source_1d == :synthetic_doside_source_1d
    @test ready_edge_context.planned_old_calls ==
          (:_nested_doside_1d, :_nested_edge_product)
    @test ready_edge_context.missing_inputs == ()
    @test !ready_edge_context.edge_facet_coefficient_maps_materialized

    real_doside_source_1d = _lw_adapter_doside_source_1d()

    real_facet_source = CPBForLWAdapter.slab_cpb(
        1:1,
        2:6,
        2:6;
        role = :lw_adapter_test_real_facet_source_cpb,
        metadata = (;
            stratum_kind = :facet_cpb,
            source_cpb_index = 5,
            fixed_axes = (:x,),
            sides = (:low,),
        ),
    )
    real_facet_unit = _lw_adapter_retained_unit(
        :lw_adapter_test_real_facet_unit,
        5,
        real_facet_source,
        :facet_cpb,
        5;
        dimension_status = :available,
        dimension = 3,
        extra_metadata = (;
            parent_dims = (7, 7, 7),
            doside_source_1d = real_doside_source_1d,
        ),
    )
    real_facet_descriptor = merge(
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            real_facet_unit,
        ),
        (;
            retained_count = 3,
            parent_dims = (7, 7, 7),
            doside_source_1d = real_doside_source_1d,
        ),
    )
    real_facet_context =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficient_context(
            real_facet_descriptor,
        )
    @test real_facet_context.status ==
          :ready_white_lindsey_facet_kernel_context_not_materialized
    @test real_facet_context.face_kind == :yz
    @test real_facet_context.fixed_side == :low
    @test real_facet_context.fixed_index == 1
    @test real_facet_context.active_product_axis_intervals == (
        (; axis = :y, interval = 2:6),
        (; axis = :z, interval = 2:6),
    )
    real_facet_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficients(
            real_facet_descriptor,
        )
    @test real_facet_coefficients.status ==
          :materialized_white_lindsey_facet_unit_coefficients
    @test isnothing(real_facet_coefficients.blocker)
    @test real_facet_coefficients.coefficient_space ==
          :parent_cartesian_sparse_adapter
    @test size(real_facet_coefficients.coefficient_matrix) == (7^3, 9)
    @test real_facet_coefficients.source_support_row_count == 25
    @test real_facet_coefficients.retained_column_count == 9
    @test length(real_facet_coefficients.support_indices) == 25
    @test issorted(real_facet_coefficients.support_indices)
    @test real_facet_coefficients.nonzero_count > 0
    @test real_facet_coefficients.old_kernels_used ==
          (:_nested_doside_1d, :_nested_face_product)
    @test real_facet_coefficients.active_axes == (:y, :z)
    @test real_facet_coefficients.active_axis_intervals == (2:6, 2:6)
    @test real_facet_coefficients.active_axis_retained_counts == (3, 3)
    @test real_facet_coefficients.retained_count_policy ==
          :scalar_retained_count_reused_for_active_axes
    @test real_facet_coefficients.doside_source_policy ==
          (:shared_doside_source_1d, :shared_doside_source_1d)
    @test real_facet_coefficients.coefficient_input_requirements.status ==
          :available_white_lindsey_facet_kernel_context_inputs
    @test real_facet_coefficients.missing_coefficient_inputs == ()
    @test real_facet_coefficients.coefficient_maps_materialized
    @test real_facet_coefficients.parent_row_indices_available
    @test !real_facet_coefficients.source_operator_blocks_materialized
    @test !real_facet_coefficients.final_pair_blocks_materialized
    @test !real_facet_coefficients.operator_blocks_materialized
    @test !real_facet_coefficients.hamiltonian_data_materialized
    @test !real_facet_coefficients.artifacts_materialized

    real_edge_source = CPBForLWAdapter.cpb(
        7:7,
        1:1,
        2:6;
        role = :lw_adapter_test_real_edge_source_cpb,
        metadata = (;
            stratum_kind = :edge_cpb,
            source_cpb_index = 4,
            fixed_axes = (:x, :y),
            sides = (:high, :low),
        ),
    )
    real_edge_unit = _lw_adapter_retained_unit(
        :lw_adapter_test_real_edge_unit,
        4,
        real_edge_source,
        :edge_cpb,
        4;
        dimension_status = :available,
        dimension = 3,
        extra_metadata = (;
            parent_dims = (7, 7, 7),
            doside_source_1d = real_doside_source_1d,
        ),
    )
    real_edge_descriptor = merge(
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            real_edge_unit,
        ),
        (;
            retained_count = 3,
            parent_dims = (7, 7, 7),
            doside_source_1d = real_doside_source_1d,
        ),
    )
    real_edge_context =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficient_context(
            real_edge_descriptor,
        )
    @test real_edge_context.status ==
          :ready_white_lindsey_edge_kernel_context_not_materialized
    @test real_edge_context.free_axis == :z
    @test real_edge_context.free_axis_interval == 2:6
    @test real_edge_context.fixed_sides == (:high, :low)
    @test real_edge_context.fixed_indices == (7, 1)
    @test real_edge_context.parent_dims == (7, 7, 7)
    real_edge_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficients(
            real_edge_descriptor,
        )
    @test real_edge_coefficients.status ==
          :materialized_white_lindsey_edge_unit_coefficients
    @test isnothing(real_edge_coefficients.blocker)
    @test real_edge_coefficients.coefficient_space ==
          :parent_cartesian_sparse_adapter
    @test size(real_edge_coefficients.coefficient_matrix) == (7^3, 3)
    @test real_edge_coefficients.source_support_row_count == 5
    @test real_edge_coefficients.retained_column_count == 3
    @test length(real_edge_coefficients.support_indices) == 5
    @test issorted(real_edge_coefficients.support_indices)
    @test real_edge_coefficients.nonzero_count > 0
    @test real_edge_coefficients.old_kernels_used ==
          (:_nested_doside_1d, :_nested_edge_product)
    @test real_edge_coefficients.coefficient_input_requirements.status ==
          :available_white_lindsey_edge_kernel_context_inputs
    @test real_edge_coefficients.missing_coefficient_inputs == ()
    @test real_edge_coefficients.coefficient_maps_materialized
    @test real_edge_coefficients.parent_row_indices_available
    @test !real_edge_coefficients.source_operator_blocks_materialized
    @test !real_edge_coefficients.final_pair_blocks_materialized
    @test !real_edge_coefficients.operator_blocks_materialized
    @test !real_edge_coefficients.hamiltonian_data_materialized
    @test !real_edge_coefficients.artifacts_materialized

    real_facet_edge_pair = _lw_adapter_unit_pair(real_facet_unit, real_edge_unit, 1)
    real_pair_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_pair_unit_coefficients(
            real_facet_edge_pair,
        )
    @test real_pair_coefficients.object_kind ==
          :white_lindsey_boundary_stratum_pair_unit_coefficients
    @test real_pair_coefficients.status ==
          :materialized_white_lindsey_pair_unit_coefficients
    @test isnothing(real_pair_coefficients.blocker)
    @test real_pair_coefficients.pair_key == (
        :lw_adapter_test_real_facet_unit,
        :lw_adapter_test_real_edge_unit,
    )
    @test real_pair_coefficients.left_stratum_kind == :facet_cpb
    @test real_pair_coefficients.right_stratum_kind == :edge_cpb
    @test real_pair_coefficients.left_coefficient_status ==
          :materialized_white_lindsey_facet_unit_coefficients
    @test real_pair_coefficients.right_coefficient_status ==
          :materialized_white_lindsey_edge_unit_coefficients
    @test size(real_pair_coefficients.left_coefficient_matrix) == (7^3, 9)
    @test size(real_pair_coefficients.right_coefficient_matrix) == (7^3, 3)
    @test length(real_pair_coefficients.left_support_indices) == 25
    @test length(real_pair_coefficients.right_support_indices) == 5
    @test real_pair_coefficients.left_retained_column_count == 9
    @test real_pair_coefficients.right_retained_column_count == 3
    @test real_pair_coefficients.unit_coefficient_cache_scope ==
          :local_pair_or_batch_call
    @test real_pair_coefficients.unit_coefficient_cache_entry_count == 2
    @test real_pair_coefficients.coefficient_maps_materialized
    @test real_pair_coefficients.pair_unit_coefficient_maps_materialized
    @test !real_pair_coefficients.pair_blocks_materialized
    @test !real_pair_coefficients.source_operator_blocks_materialized
    @test !real_pair_coefficients.final_pair_blocks_materialized
    @test !real_pair_coefficients.operator_blocks_materialized
    @test !real_pair_coefficients.hamiltonian_data_materialized
    @test !real_pair_coefficients.artifacts_materialized

    overlap_1d = (;
        x = ones(Float64, 7, 7),
        y = ones(Float64, 7, 7),
        z = ones(Float64, 7, 7),
    )
    overlap_result =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_overlap_block(
            real_pair_coefficients;
            parent_axis_counts = (7, 7, 7),
            overlap_1d,
        )
    @test overlap_result.term == :overlap
    @test overlap_result.pair_key == real_pair_coefficients.pair_key
    @test size(overlap_result.block) == (9, 3)
    @test all(isfinite, overlap_result.block)
    @test sum(abs, overlap_result.block) > 0.0
    @test overlap_result.materialized
    @test overlap_result.source_operator_blocks_materialized
    @test overlap_result.final_pair_blocks_materialized
    @test !overlap_result.operator_blocks_materialized
    @test !overlap_result.hamiltonian_data_materialized
    @test !overlap_result.artifacts_materialized
    @test overlap_result.metadata.materialization_path ==
          :white_lindsey_boundary_stratum_overlap_adapter
    @test overlap_result.metadata.left_stratum_kind == :facet_cpb
    @test overlap_result.metadata.right_stratum_kind == :edge_cpb
    @test overlap_result.metadata.left_support_count == 25
    @test overlap_result.metadata.right_support_count == 5
    @test overlap_result.metadata.left_retained_column_count == 9
    @test overlap_result.metadata.right_retained_column_count == 3
    @test overlap_result.metadata.parent_axis_counts == (7, 7, 7)
    @test overlap_result.metadata.support_overlap_shape == (25, 5)
    @test overlap_result.metadata.local_pair_block_materialized
    @test !overlap_result.metadata.operator_blocks_materialized
    @test !overlap_result.metadata.hamiltonian_data_materialized
    @test !overlap_result.metadata.artifacts_materialized
    @test !overlap_result.metadata.dense_parent_parent_overlap_materialized

    blocked_pair_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_pair_unit_coefficients(
            ready_facet_descriptor,
            ready_edge_descriptor;
            pair_index = 2,
        )
    @test blocked_pair_coefficients.status ==
          :blocked_white_lindsey_pair_unit_coefficients
    @test blocked_pair_coefficients.blocker ==
          :left_white_lindsey_unit_coefficients_not_materialized
    @test blocked_pair_coefficients.left_coefficient_status ==
          :blocked_white_lindsey_boundary_stratum_unit_coefficients
    @test blocked_pair_coefficients.right_coefficient_status ==
          :blocked_white_lindsey_boundary_stratum_unit_coefficients
    @test blocked_pair_coefficients.left_blocker ==
          :white_lindsey_facet_doside_source_not_materializable
    @test blocked_pair_coefficients.right_blocker ==
          :white_lindsey_edge_doside_source_not_materializable
    @test isnothing(blocked_pair_coefficients.left_coefficient_matrix)
    @test isnothing(blocked_pair_coefficients.right_coefficient_matrix)
    @test !blocked_pair_coefficients.coefficient_maps_materialized
    @test !blocked_pair_coefficients.pair_blocks_materialized

    edge_corner_pair = _lw_adapter_unit_pair(real_edge_unit, corner_unit, 3)
    edge_corner_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_pair_unit_coefficients(
            edge_corner_pair,
        )
    @test edge_corner_coefficients.status ==
          :materialized_white_lindsey_pair_unit_coefficients
    @test edge_corner_coefficients.left_coefficient_status ==
          :materialized_white_lindsey_edge_unit_coefficients
    @test edge_corner_coefficients.right_coefficient_status ==
          :materialized_white_lindsey_corner_unit_coefficients
    @test size(edge_corner_coefficients.left_coefficient_matrix) == (7^3, 3)
    @test size(edge_corner_coefficients.right_coefficient_matrix) == (1, 1)
    @test edge_corner_coefficients.right_retained_column_count == 1
    @test edge_corner_coefficients.unit_coefficient_cache_entry_count == 2
    @test edge_corner_coefficients.pair_unit_coefficient_maps_materialized
    @test !edge_corner_coefficients.pair_blocks_materialized

    ready_corner_descriptor = merge(
        corner_descriptor,
        (;
            fixed_side_metadata = (
                (; axis = :x, side = :high),
                (; axis = :y, side = :high),
                (; axis = :z, side = :high),
            ),
            parent_dims = (7, 7, 7),
        ),
    )
    ready_corner_context =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficient_context(
            ready_corner_descriptor,
        )
    @test ready_corner_context.status ==
          :ready_white_lindsey_corner_kernel_context_not_materialized
    @test isnothing(ready_corner_context.blocker)
    @test ready_corner_context.fixed_sides == (:high, :high, :high)
    @test ready_corner_context.fixed_indices == (4, 3, 3)
    @test ready_corner_context.parent_dims == (7, 7, 7)
    @test ready_corner_context.planned_old_calls == (:_nested_corner_piece,)
    @test ready_corner_context.missing_inputs == ()
    @test !ready_corner_context.coefficient_maps_materialized

    ready_corner_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficients(
            ready_corner_descriptor,
        )
    @test ready_corner_coefficients.status ==
          :materialized_white_lindsey_corner_unit_coefficients
    @test ready_corner_coefficients.coefficient_matrix == [1.0;;]

    bad_descriptor = (;
        object_kind = :white_lindsey_boundary_stratum_unit_adapter_descriptor,
        status = :available_metadata_only_white_lindsey_unit_adapter_descriptor,
        unit_key = :lw_adapter_test_bad_corner_unit,
        stratum_kind = :corner_cpb,
        planned_old_kernel = :_nested_corner_piece,
    )
    bad_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficients(
            bad_descriptor,
        )
    @test bad_coefficients.status ==
          :blocked_white_lindsey_boundary_stratum_unit_coefficients
    @test bad_coefficients.blocker ==
          :white_lindsey_corner_source_cpb_not_support_local
    @test !bad_coefficients.coefficient_maps_materialized
end

@testset "CartesianPairBlockMaterialization White-Lindsey seed oracle summary" begin
    oracle_summary =
        CPBMForLWAdapter.white_lindsey_materialized_seed_oracle_summary()

    @test oracle_summary.object_kind ==
          :white_lindsey_materialized_seed_oracle_summary
    @test oracle_summary.status ==
          :available_white_lindsey_materialized_seed_oracle_summary
    @test oracle_summary.oracle_role == :validation_oracle_only
    @test !oracle_summary.route_authority
    @test !oracle_summary.adapter_authority
    @test oracle_summary.private_development_only
    @test oracle_summary.seed_report_kind ==
          :white_lindsey_low_order_materialized_seed_report
    @test oracle_summary.seed_report_status == :private_development_seed
    @test oracle_summary.route_family == :white_lindsey_low_order
    @test oracle_summary.shellization_source == :white_lindsey_one_center_seed
    @test oracle_summary.packet_kernel == :factorized_direct
    @test oracle_summary.packet_inventory_available

    @test oracle_summary.retained_dimension == 223
    @test oracle_summary.retained_unit_count == 4
    @test oracle_summary.unit_keys == (
        :low_order_core_direct,
        :low_order_face_interiors,
        :low_order_edges,
        :low_order_corners,
    )
    @test oracle_summary.unit_roles ==
          (:direct_core, :face_interiors, :edges, :corners)
    @test oracle_summary.retained_unit_kinds == (
        :white_lindsey_direct_core,
        :white_lindsey_face_interior_2d_products,
        :white_lindsey_edge_1d_side_functions,
        :white_lindsey_corner_direct_single_sites,
    )
    @test oracle_summary.retained_counts == (
        low_order_core_direct = 125,
        low_order_face_interiors = 54,
        low_order_edges = 36,
        low_order_corners = 8,
    )
    @test oracle_summary.retained_ranges == (
        low_order_core_direct = 1:125,
        low_order_face_interiors = 126:179,
        low_order_edges = 180:215,
        low_order_corners = 216:223,
    )
    @test oracle_summary.piece_counts ==
          (core = 1, faces = 6, edges = 12, corners = 8)
    @test oracle_summary.support_counts ==
          (core = 125, shell = 218, total_source = 343)
    @test oracle_summary.seed_retained_counts == (
        core = 125,
        faces = 54,
        edges = 36,
        corners = 8,
        shell = 98,
        total = 223,
    )
    @test oracle_summary.fixed_block_ready
    @test oracle_summary.overlap_ready
    @test oracle_summary.retained_basis_integral_weights_ready
    @test oracle_summary.weight_semantics == :retained_basis_integral_weights

    @test oracle_summary.operator_inventory_available
    @test oracle_summary.operator_source == :nested_fixed_block
    @test oracle_summary.operator_terms == (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    @test all(values(oracle_summary.one_body_operator_matrix_available))
    @test oracle_summary.fixed_block_operator_matrix_sizes.overlap == (223, 223)
    @test oracle_summary.fixed_block_operator_matrix_sizes.kinetic == (223, 223)
    @test oracle_summary.all_operator_matrices_finite
    @test oracle_summary.operator_symmetric_ready.overlap
    @test oracle_summary.operator_symmetric_ready.kinetic
    @test oracle_summary.overlap_identity_ready
    @test !oracle_summary.operator_pairs_materialized
    @test !oracle_summary.electron_electron_materialized
    @test !oracle_summary.new_adapter_coefficient_maps_materialized
    @test !oracle_summary.new_adapter_pair_blocks_materialized
    @test !oracle_summary.hamiltonian_data_materialized
    @test !oracle_summary.artifacts_materialized

    @test !haskey(oracle_summary, :fixture)
    @test !haskey(oracle_summary, :fixed_block)
    @test !haskey(oracle_summary, :operator_matrices)
end
