# Runtime role: CPB-local axis-product operator block contract.
#
# This test validates the generic CPB-local axis-product dense block primitive
# and the overlap wrapper around it. It does not add WL/PQS realization,
# retained transforms, route/global placement, driver wiring, kinetic,
# position, x2, Coulomb, Hamiltonian, IDA/MWG, PQS Lowdin/projection, exports,
# or artifacts.

using Test
using GaussletBases

const CPBOperatorBlock = GaussletBases.CartesianCPB
const CPGBOperatorBlock = GaussletBases.CartesianParentGaussletBases
const CBPOperatorBlock = GaussletBases.CartesianCPBBlockProviders

function _axis_product_expected_dense(axis_ops)
    nx_left, nx_right = size(axis_ops.x)
    ny_left, ny_right = size(axis_ops.y)
    nz_left, nz_right = size(axis_ops.z)
    element_type = promote_type(
        eltype(axis_ops.x),
        eltype(axis_ops.y),
        eltype(axis_ops.z),
    )
    dense = Matrix{element_type}(
        undef,
        nx_left * ny_left * nz_left,
        nx_right * ny_right * nz_right,
    )
    for ix_left in 1:nx_left, iy_left in 1:ny_left, iz_left in 1:nz_left
        left_index =
            (ix_left - 1) * ny_left * nz_left +
            (iy_left - 1) * nz_left +
            iz_left
        for ix_right in 1:nx_right, iy_right in 1:ny_right, iz_right in 1:nz_right
            right_index =
                (ix_right - 1) * ny_right * nz_right +
                (iy_right - 1) * nz_right +
                iz_right
            dense[left_index, right_index] =
                axis_ops.x[ix_left, ix_right] *
                axis_ops.y[iy_left, iy_right] *
                axis_ops.z[iz_left, iz_right]
        end
    end
    return dense
end

function _axis_product_case(axis_ops; expected_eltype = Float64)
    block = CBPOperatorBlock.cpb_axis_product_operator_block(
        axis_ops;
        term = :test_axis_product,
        factor_space = :test_factor_space,
        factor_convention = :test_factor_convention,
        index_domain = :parent_axis_indices,
    )
    block_summary = CBPOperatorBlock.summary(block)
    @test block_summary.status === :materialized_cpb_axis_product_operator_block
    @test block_summary.blocker === nothing
    @test block_summary.term === :test_axis_product
    @test block_summary.representation === :dense_local_cpb_product_space
    @test block_summary.local_ordering ===
          :parent_compatible_x_slowest_z_fastest
    @test block_summary.left_shape == (
        x = size(axis_ops.x, 1),
        y = size(axis_ops.y, 1),
        z = size(axis_ops.z, 1),
    )
    @test block_summary.right_shape == (
        x = size(axis_ops.x, 2),
        y = size(axis_ops.y, 2),
        z = size(axis_ops.z, 2),
    )
    @test block_summary.dense_block_shape ==
          (block_summary.left_support_count, block_summary.right_support_count)
    @test block_summary.dense_block_eltype === expected_eltype
    @test block_summary.provider_level_local_matrix_materialized === true
    @test block_summary.realization_status === :unrealized
    @test block_summary.route_global_status === :unassigned
    @test block_summary.route_driver_wiring === false
    @test block_summary.route_global_matrix_materialized === false
    @test block_summary.global_matrix_materialized === false
    @test !hasproperty(block_summary, :dense_block)
    @test !hasproperty(block_summary, :axis_ops)
    @test !hasproperty(block_summary, :global_overlap_matrix)
    @test !hasproperty(block_summary, :retained_blocks)
    @test block.dense_block == _axis_product_expected_dense(axis_ops)
    return block
end

function _test_blocked_axis_product(axis_ops, expected_blocker)
    block = CBPOperatorBlock.cpb_axis_product_operator_block(axis_ops)
    block_summary = CBPOperatorBlock.summary(block)
    @test block_summary.status === :blocked_cpb_axis_product_operator_block
    @test block_summary.blocker === expected_blocker
    @test isnothing(block.dense_block)
    @test block_summary.dense_block_available === false
    @test block_summary.dense_block_shape === :unavailable
    @test block_summary.dense_block_eltype === :unavailable
    @test block_summary.provider_level_local_matrix_materialized === false
    @test block_summary.realization_status === :unrealized
    @test block_summary.route_global_status === :unassigned
    @test block_summary.route_driver_wiring === false
    @test block_summary.route_global_matrix_materialized === false
    @test block_summary.global_matrix_materialized === false
    @test !hasproperty(block_summary, :dense_block)
    @test !hasproperty(block_summary, :axis_ops)
    @test !hasproperty(block_summary, :global_overlap_matrix)
    @test !hasproperty(block_summary, :retained_blocks)
    return nothing
end

function _sum_axis_products_expected_dense(product_terms)
    first_dense = _axis_product_expected_dense(first(product_terms).axis_ops)
    dense_eltype = Union{}
    for product_term in product_terms
        dense_eltype = promote_type(
            dense_eltype,
            typeof(product_term.coefficient),
            eltype(product_term.axis_ops.x),
            eltype(product_term.axis_ops.y),
            eltype(product_term.axis_ops.z),
        )
    end
    dense = zeros(dense_eltype, size(first_dense))
    for product_term in product_terms
        dense .+= product_term.coefficient .*
                  _axis_product_expected_dense(product_term.axis_ops)
    end
    return dense
end

function _sum_axis_products_case(product_terms; expected_dense)
    block = CBPOperatorBlock.cpb_sum_of_axis_products_operator_block(
        product_terms;
        term = :test_sum_of_axis_products,
        factor_space = :test_factor_space,
        factor_convention = :test_factor_convention,
        index_domain = :parent_axis_indices,
    )
    block_summary = CBPOperatorBlock.summary(block)
    @test block_summary.status ===
          :materialized_cpb_sum_of_axis_products_operator_block
    @test block_summary.blocker === nothing
    @test block_summary.term === :test_sum_of_axis_products
    @test block_summary.representation === :dense_local_cpb_sum_of_axis_products
    @test block_summary.product_term_count == length(product_terms)
    @test block_summary.product_term_labels ==
          Tuple(product_term.label for product_term in product_terms)
    @test block_summary.dense_block_available === true
    @test block_summary.dense_block_shape == size(expected_dense)
    @test block_summary.dense_block_eltype === eltype(expected_dense)
    @test block_summary.provider_level_local_matrix_materialized === true
    @test block_summary.realization_status === :unrealized
    @test block_summary.route_global_status === :unassigned
    @test block_summary.route_driver_wiring === false
    @test block_summary.route_global_matrix_materialized === false
    @test block_summary.global_matrix_materialized === false
    @test !hasproperty(block_summary, :dense_block)
    @test !hasproperty(block_summary, :axis_ops)
    @test !hasproperty(block_summary, :global_overlap_matrix)
    @test !hasproperty(block_summary, :retained_blocks)
    @test block.dense_block == expected_dense
    return block
end

function _test_blocked_sum_axis_products(product_terms, expected_blocker)
    block = CBPOperatorBlock.cpb_sum_of_axis_products_operator_block(product_terms)
    block_summary = CBPOperatorBlock.summary(block)
    @test block_summary.status ===
          :blocked_cpb_sum_of_axis_products_operator_block
    @test block_summary.blocker === expected_blocker
    @test isnothing(block.dense_block)
    @test block_summary.dense_block_available === false
    @test block_summary.dense_block_shape === :unavailable
    @test block_summary.dense_block_eltype === :unavailable
    @test block_summary.provider_level_local_matrix_materialized === false
    @test block_summary.realization_status === :unrealized
    @test block_summary.route_global_status === :unassigned
    @test block_summary.route_driver_wiring === false
    @test block_summary.route_global_matrix_materialized === false
    @test block_summary.global_matrix_materialized === false
    @test !hasproperty(block_summary, :dense_block)
    @test !hasproperty(block_summary, :axis_ops)
    @test !hasproperty(block_summary, :global_overlap_matrix)
    @test !hasproperty(block_summary, :retained_blocks)
    return block
end

function _operator_block_parent(; count = 3)
    axis = build_basis(MappedUniformBasisSpec(
        :G10;
        count,
        mapping = IdentityMapping(),
        reference_spacing = 1.0,
    ))
    return CPGBOperatorBlock.CartesianParentGaussletBasis3D(axis)
end

function _operator_block_overlap_1d()
    return (;
        x = [
            1.0 0.1 0.2
            0.3 1.1 0.4
            0.5 0.6 1.2
        ],
        y = [
            2.0 0.7 0.8
            0.9 2.1 1.0
            1.1 1.2 2.2
        ],
        z = [
            3.0 1.3 1.4
            1.5 3.1 1.6
            1.7 1.8 3.2
        ],
    )
end

function _operator_block_kinetic_1d()
    return (;
        x = [
            2.0 -0.1 -0.2
            -0.3 2.1 -0.4
            -0.5 -0.6 2.2
        ],
        y = [
            2.3 -0.7 -0.8
            -0.9 2.4 -1.0
            -1.1 -1.2 2.5
        ],
        z = [
            2.6 -1.3 -1.4
            -1.5 2.7 -1.6
            -1.7 -1.8 2.8
        ],
    )
end

function _operator_block_axis_bundle(overlap_1d; kinetic_1d = nothing)
    axis_bundle(axis, overlap_matrix) =
        isnothing(kinetic_1d) ?
        (; pgdg_intermediate = (; overlap = overlap_matrix)) :
        (; pgdg_intermediate = (;
            overlap = overlap_matrix,
            kinetic = getproperty(kinetic_1d, axis),
        ))
    return (;
        x = axis_bundle(:x, overlap_1d.x),
        y = axis_bundle(:y, overlap_1d.y),
        z = axis_bundle(:z, overlap_1d.z),
    )
end

@testset "CPB axis-product operator block" begin
    point_ops = (;
        x = fill(2.0, 1, 1),
        y = fill(3.0, 1, 1),
        z = fill(5.0, 1, 1),
    )
    point_block = _axis_product_case(point_ops)
    @test point_block.dense_block == fill(30.0, 1, 1)

    edge_ops = (;
        x = [1.0 2.0; 3.0 4.0],
        y = fill(5.0, 1, 1),
        z = fill(7.0, 1, 1),
    )
    edge_block = _axis_product_case(edge_ops)
    @test edge_block.dense_block == 35.0 .* edge_ops.x

    face_ops = (;
        x = reshape([1.0, 2.0], 2, 1),
        y = [3.0 4.0; 5.0 6.0],
        z = fill(7.0, 1, 1),
    )
    _axis_product_case(face_ops)

    cube_ops = (;
        x = [1.0 2.0; 3.0 4.0],
        y = reshape([5.0, 6.0], 2, 1),
        z = [7.0 8.0; 9.0 10.0],
    )
    cube_block = _axis_product_case(cube_ops)
    @test CBPOperatorBlock.summary(cube_block).left_shape == (x = 2, y = 2, z = 2)
    @test CBPOperatorBlock.summary(cube_block).right_shape == (x = 2, y = 1, z = 2)

    mixed_ops = (;
        x = Int16[1 2],
        y = reshape(Float32[3.0, 4.0], 2, 1),
        z = [5.0 6.0],
    )
    mixed_eltype = promote_type(Int16, Float32, Float64)
    mixed_block = _axis_product_case(mixed_ops; expected_eltype = mixed_eltype)
    @test eltype(mixed_block.dense_block) === mixed_eltype

    _test_blocked_axis_product(
        (; y = fill(1.0, 1, 1), z = fill(1.0, 1, 1)),
        :missing_x_axis_operator,
    )
    _test_blocked_axis_product(
        (; x = fill(1.0, 1, 1), z = fill(1.0, 1, 1)),
        :missing_y_axis_operator,
    )
    _test_blocked_axis_product(
        (; x = fill(1.0, 1, 1), y = fill(1.0, 1, 1)),
        :missing_z_axis_operator,
    )
    _test_blocked_axis_product(
        (; x = [1.0, 2.0], y = fill(1.0, 1, 1), z = fill(1.0, 1, 1)),
        :x_axis_operator_not_matrix,
    )
    _test_blocked_axis_product(
        (;
            x = Matrix{Float64}(undef, 0, 1),
            y = fill(1.0, 1, 1),
            z = fill(1.0, 1, 1),
        ),
        :x_axis_operator_empty,
    )
    _test_blocked_axis_product(
        (;
            x = fill(1.0, 1, 1),
            y = Matrix{Float64}(undef, 1, 0),
            z = fill(1.0, 1, 1),
        ),
        :y_axis_operator_empty,
    )
    _test_blocked_axis_product(
        (;
            x = fill(1.0, 1, 1),
            y = fill(1.0, 1, 1),
            z = Matrix{Float64}(undef, 0, 0),
        ),
        :z_axis_operator_empty,
    )

    one_term = (
        (;
            coefficient = 2.0,
            axis_ops = edge_ops,
            label = :scaled_edge_product,
        ),
    )
    one_term_sum = _sum_axis_products_case(
        one_term;
        expected_dense = 2.0 .* edge_block.dense_block,
    )
    @test one_term_sum.dense_block == 2.0 .* edge_block.dense_block

    two_term_ops = (;
        x = [0.5 1.0; 1.5 2.0],
        y = fill(3.0, 1, 1),
        z = fill(4.0, 1, 1),
    )
    two_terms = (
        (; coefficient = 1.5, axis_ops = edge_ops, label = :left_product),
        (; coefficient = -0.25, axis_ops = two_term_ops, label = :right_product),
    )
    _sum_axis_products_case(
        two_terms;
        expected_dense = _sum_axis_products_expected_dense(two_terms),
    )

    sx = [1.0 0.1; 0.2 1.1]
    sy = [1.2 0.3; 0.4 1.3]
    sz = [1.4 0.5; 0.6 1.5]
    kx = [2.0 -0.7; -0.8 2.1]
    ky = [2.2 -0.9; -1.0 2.3]
    kz = [2.4 -1.1; -1.2 2.5]
    kinetic_shaped_terms = (
        (; coefficient = 1.0, axis_ops = (x = kx, y = sy, z = sz), label = :x_piece),
        (; coefficient = 1.0, axis_ops = (x = sx, y = ky, z = sz), label = :y_piece),
        (; coefficient = 1.0, axis_ops = (x = sx, y = sy, z = kz), label = :z_piece),
    )
    kinetic_shaped_block = _sum_axis_products_case(
        kinetic_shaped_terms;
        expected_dense = _sum_axis_products_expected_dense(kinetic_shaped_terms),
    )
    @test CBPOperatorBlock.summary(kinetic_shaped_block).product_term_labels ==
          (:x_piece, :y_piece, :z_piece)

    _test_blocked_sum_axis_products((), :empty_axis_product_term_list)

    shape_mismatch_terms = (
        (; coefficient = 1.0, axis_ops = edge_ops, label = :base_shape),
        (; coefficient = 1.0, axis_ops = face_ops, label = :other_shape),
    )
    _test_blocked_sum_axis_products(
        shape_mismatch_terms,
        :axis_product_term_shape_mismatch,
    )

    blocked_product_terms = (
        (;
            coefficient = 1.0,
            axis_ops = (x = [1.0, 2.0], y = fill(1.0, 1, 1), z = fill(1.0, 1, 1)),
            label = :bad_axis_product,
        ),
    )
    blocked_product_sum = _test_blocked_sum_axis_products(
        blocked_product_terms,
        :axis_product_term_blocked,
    )
    blocked_product_summary = CBPOperatorBlock.summary(blocked_product_sum)
    @test blocked_product_summary.product_term_summaries[1].axis_product_blocker ===
          :x_axis_operator_not_matrix

    coefficient_blocked_sum = _test_blocked_sum_axis_products(
        ((; axis_ops = edge_ops, label = :missing_coefficient),),
        :axis_product_term_coefficient_unavailable,
    )
    coefficient_blocked_summary =
        CBPOperatorBlock.summary(coefficient_blocked_sum)
    @test coefficient_blocked_summary.product_term_summaries[1].coefficient_available ===
          false

    parent = _operator_block_parent()
    overlap_1d = _operator_block_overlap_1d()
    packet = CPGBOperatorBlock.parent_overlap_axis_factor_packet(
        parent,
        _operator_block_axis_bundle(overlap_1d),
    )
    left = CPBOperatorBlock.cpb(1:2, 2:3, 1:3; role = :operator_left)
    right = CPBOperatorBlock.cpb(2:3, 1:2, 2:2; role = :operator_right)
    interval_pair = CBPOperatorBlock.cpb_interval_pair(parent, left, right)
    overlap_operator =
        CBPOperatorBlock.cpb_overlap_operator_block(packet, interval_pair)
    overlap_summary = CBPOperatorBlock.summary(overlap_operator)
    axis_block_set = CBPOperatorBlock.cpb_overlap_axis_blocks(packet, interval_pair)
    existing_dense = CBPOperatorBlock.cpb_overlap_dense_block(axis_block_set)

    @test overlap_summary.status === :materialized_cpb_overlap_operator_block
    @test overlap_summary.blocker === nothing
    @test overlap_summary.term === :overlap
    @test overlap_summary.representation === :dense_local_cpb_product_space
    @test overlap_summary.realization_status === :unrealized
    @test overlap_summary.route_global_status === :unassigned
    @test overlap_summary.route_driver_wiring === false
    @test overlap_summary.route_global_matrix_materialized === false
    @test overlap_summary.provider_level_local_matrix_materialized === true
    @test overlap_summary.factor_convention === :axis_bundle_one_body_overlap
    @test !isnothing(overlap_operator.axis_product_block)
    @test overlap_operator.axis_product_block.dense_block == existing_dense.dense_block
    @test existing_dense.dense_block == overlap_operator.axis_product_block.dense_block
    @test CBPOperatorBlock.summary(overlap_operator.axis_product_block).term === :overlap
    @test !hasproperty(overlap_summary, :dense_block)
    @test !hasproperty(overlap_summary, :global_overlap_matrix)
    @test !hasproperty(overlap_summary, :axis_ops)
    @test !hasproperty(overlap_summary, :retained_blocks)

    missing_kinetic_operator =
        CBPOperatorBlock.cpb_kinetic_operator_block(packet, interval_pair)
    missing_kinetic_summary =
        CBPOperatorBlock.summary(missing_kinetic_operator)

    @test missing_kinetic_summary.status === :blocked_cpb_kinetic_operator_block
    @test missing_kinetic_summary.blocker ===
          :missing_parent_axis_bundle_kinetic_factors
    @test isnothing(missing_kinetic_operator.sum_axis_product_block)
    @test missing_kinetic_summary.dense_block_available === false
    @test missing_kinetic_summary.route_driver_wiring === false
    @test missing_kinetic_summary.global_matrix_materialized === false

    kinetic_1d = _operator_block_kinetic_1d()
    kinetic_packet = CPGBOperatorBlock.parent_overlap_axis_factor_packet(
        parent,
        _operator_block_axis_bundle(overlap_1d; kinetic_1d),
    )
    kinetic_operator =
        CBPOperatorBlock.cpb_kinetic_operator_block(kinetic_packet, interval_pair)
    kinetic_summary = CBPOperatorBlock.summary(kinetic_operator)
    left_intervals = CBPOperatorBlock.summary(interval_pair).left_intervals
    right_intervals = CBPOperatorBlock.summary(interval_pair).right_intervals
    overlap_blocks = (;
        x = view(overlap_1d.x, left_intervals.x, right_intervals.x),
        y = view(overlap_1d.y, left_intervals.y, right_intervals.y),
        z = view(overlap_1d.z, left_intervals.z, right_intervals.z),
    )
    kinetic_blocks = (;
        x = view(kinetic_1d.x, left_intervals.x, right_intervals.x),
        y = view(kinetic_1d.y, left_intervals.y, right_intervals.y),
        z = view(kinetic_1d.z, left_intervals.z, right_intervals.z),
    )
    expected_kinetic_sum =
        CBPOperatorBlock.cpb_sum_of_axis_products_operator_block((
            (;
                coefficient = 1.0,
                axis_ops = (x = kinetic_blocks.x, y = overlap_blocks.y, z = overlap_blocks.z),
                label = :kinetic_x_component,
            ),
            (;
                coefficient = 1.0,
                axis_ops = (x = overlap_blocks.x, y = kinetic_blocks.y, z = overlap_blocks.z),
                label = :kinetic_y_component,
            ),
            (;
                coefficient = 1.0,
                axis_ops = (x = overlap_blocks.x, y = overlap_blocks.y, z = kinetic_blocks.z),
                label = :kinetic_z_component,
            ),
        ))

    @test kinetic_summary.status === :materialized_cpb_kinetic_operator_block
    @test kinetic_summary.blocker === nothing
    @test kinetic_summary.term === :kinetic
    @test kinetic_summary.representation === :dense_local_cpb_sum_of_axis_products
    @test kinetic_summary.kinetic_factor_form === :sum_of_axis_products
    @test kinetic_summary.kinetic_component_axes == (:x, :y, :z)
    @test kinetic_summary.product_term_labels ==
          (:kinetic_x_component, :kinetic_y_component, :kinetic_z_component)
    @test kinetic_summary.factor_space === :parent_axis_bundle_pgdg_intermediate
    @test kinetic_summary.factor_convention === :axis_bundle_one_body_kinetic_sum
    @test kinetic_summary.index_domain === :parent_axis_indices
    @test kinetic_summary.provider_level_local_matrix_materialized === true
    @test kinetic_summary.realization_status === :unrealized
    @test kinetic_summary.route_global_status === :unassigned
    @test kinetic_summary.route_driver_wiring === false
    @test kinetic_summary.route_global_matrix_materialized === false
    @test kinetic_summary.global_matrix_materialized === false
    @test !hasproperty(kinetic_summary, :dense_block)
    @test !hasproperty(kinetic_summary, :axis_ops)
    @test !hasproperty(kinetic_summary, :global_overlap_matrix)
    @test !hasproperty(kinetic_summary, :retained_blocks)
    @test !isnothing(kinetic_operator.sum_axis_product_block)
    @test kinetic_operator.sum_axis_product_block.dense_block ==
          expected_kinetic_sum.dense_block
end
