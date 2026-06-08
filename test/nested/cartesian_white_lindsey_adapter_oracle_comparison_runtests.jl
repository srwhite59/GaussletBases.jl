using Test

include("cartesian_white_lindsey_adapter_fixture_helpers.jl")

function _lw_adapter_oracle_pair_family(adapter_block)
    metadata = adapter_block.metadata
    return Symbol(
        String(metadata.left_stratum_kind),
        "__",
        String(metadata.right_stratum_kind),
    )
end

function _lw_adapter_seed_compatible_doside_source_1d(seed; expansion)
    return GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        seed.basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
end

function _lw_adapter_seed_local_facet_edge_facts(seed_report)
    seed = seed_report.fixture
    shell = only(seed.sequence.shell_layers)
    face_index = findfirst(
        face -> (
            face.face_kind === :yz &&
            face.fixed_axis === :x &&
            face.fixed_side === :low
        ),
        shell.faces,
    )
    edge_index = findfirst(
        edge -> (
            edge.free_axis === :z &&
            edge.fixed_axes == (:x, :y) &&
            edge.fixed_sides == (:high, :low)
        ),
        shell.edges,
    )
    corner_index = findfirst(
        corner -> corner.fixed_sides == (:high, :low, :high),
        shell.corners,
    )
    if isnothing(face_index) || isnothing(edge_index) || isnothing(corner_index)
        return (;
            object_kind = :white_lindsey_seed_local_pair_facts,
            status = :blocked_seed_local_pair_slice_not_found,
            blocker = :old_seed_local_pair_slice_not_exposed,
            route_unit_range_status = :grouped_ranges_available,
            per_piece_range_status = :missing_representative_piece_range,
        )
    end

    face_range = seed.inventory.retained_ranges.faces[face_index]
    edge_range = seed.inventory.retained_ranges.edges[edge_index]
    corner_range = seed.inventory.retained_ranges.corners[corner_index]
    face = shell.faces[face_index]
    edge = shell.edges[edge_index]
    corner = shell.corners[corner_index]
    route_ranges = seed_report.route_units.unit_inventory.ranges

    return (;
        object_kind = :white_lindsey_seed_local_pair_facts,
        status = :available_seed_local_pair_slice,
        blocker = nothing,
        route_unit_range_status =
            :grouped_route_ranges_available_not_local_slice_authority,
        per_piece_range_status = :per_piece_seed_inventory_ranges_available,
        per_piece_range_source = :seed_inventory_retained_ranges,
        face_index,
        edge_index,
        face_signature = (;
            face_kind = face.face_kind,
            fixed_axis = face.fixed_axis,
            fixed_side = face.fixed_side,
            fixed_index = face.fixed_index,
        ),
        edge_signature = (;
            free_axis = edge.free_axis,
            fixed_axes = edge.fixed_axes,
            fixed_sides = edge.fixed_sides,
            fixed_indices = edge.fixed_indices,
        ),
        corner_index,
        corner_signature = (;
            fixed_sides = corner.fixed_sides,
            fixed_indices = corner.fixed_indices,
        ),
        face_global_range = face_range,
        edge_global_range = edge_range,
        corner_global_range = corner_range,
        face_local_range =
            seed.inventory.materialized_shell_local_ranges.faces[face_index],
        edge_local_range =
            seed.inventory.materialized_shell_local_ranges.edges[edge_index],
        corner_local_range =
            seed.inventory.materialized_shell_local_ranges.corners[corner_index],
        grouped_face_range = route_ranges.low_order_face_interiors,
        grouped_edge_range = route_ranges.low_order_edges,
        grouped_corner_range = route_ranges.low_order_corners,
        slice_shape = (length(face_range), length(edge_range)),
        edge_corner_slice_shape = (length(edge_range), length(corner_range)),
        overlap_slice_available = true,
        validation_oracle_only = true,
        route_authority = false,
        adapter_authority = false,
    )
end

function _lw_adapter_real_corner_unit(
    doside_source_1d;
    prefix::AbstractString = "lw_oracle_comparison",
)
    real_corner_source = CPBForLWAdapter.cpb(
        7:7,
        1:1,
        7:7;
        role = Symbol(prefix, "_real_corner_source_cpb"),
        metadata = (;
            stratum_kind = :corner_cpb,
            source_cpb_index = 6,
            fixed_axes = (:x, :y, :z),
            sides = (:high, :low, :high),
        ),
    )
    return _lw_adapter_retained_unit(
        Symbol(prefix, "_real_corner_unit"),
        6,
        real_corner_source,
        :corner_cpb,
        6;
        dimension_status = :available,
        dimension = 1,
        extra_metadata = (;
            parent_dims = (7, 7, 7),
            doside_source_1d,
        ),
    )
end

function _lw_adapter_seed_one_body_factors(doside_source_1d)
    pgdg = doside_source_1d.pgdg_intermediate
    return (;
        overlap_1d = (;
            x = pgdg.overlap,
            y = pgdg.overlap,
            z = pgdg.overlap,
        ),
        position_1d = (;
            x = pgdg.position,
            y = pgdg.position,
            z = pgdg.position,
        ),
        x2_1d = (;
            x = pgdg.x2,
            y = pgdg.x2,
            z = pgdg.x2,
        ),
        kinetic_1d = (;
            x = pgdg.kinetic,
            y = pgdg.kinetic,
            z = pgdg.kinetic,
        ),
    )
end

function _lw_adapter_seed_one_body_blocks(pair_unit_coefficients, factors)
    return (;
        overlap = CPBMForLWAdapter.white_lindsey_boundary_stratum_overlap_block(
            pair_unit_coefficients;
            parent_axis_counts = (7, 7, 7),
            overlap_1d = factors.overlap_1d,
        ),
        position_x = CPBMForLWAdapter.white_lindsey_boundary_stratum_position_block(
            pair_unit_coefficients;
            axis = :x,
            parent_axis_counts = (7, 7, 7),
            overlap_1d = factors.overlap_1d,
            position_1d = factors.position_1d,
        ),
        position_y = CPBMForLWAdapter.white_lindsey_boundary_stratum_position_block(
            pair_unit_coefficients;
            axis = :y,
            parent_axis_counts = (7, 7, 7),
            overlap_1d = factors.overlap_1d,
            position_1d = factors.position_1d,
        ),
        position_z = CPBMForLWAdapter.white_lindsey_boundary_stratum_position_block(
            pair_unit_coefficients;
            axis = :z,
            parent_axis_counts = (7, 7, 7),
            overlap_1d = factors.overlap_1d,
            position_1d = factors.position_1d,
        ),
        x2_x = CPBMForLWAdapter.white_lindsey_boundary_stratum_x2_block(
            pair_unit_coefficients;
            axis = :x,
            parent_axis_counts = (7, 7, 7),
            overlap_1d = factors.overlap_1d,
            x2_1d = factors.x2_1d,
        ),
        x2_y = CPBMForLWAdapter.white_lindsey_boundary_stratum_x2_block(
            pair_unit_coefficients;
            axis = :y,
            parent_axis_counts = (7, 7, 7),
            overlap_1d = factors.overlap_1d,
            x2_1d = factors.x2_1d,
        ),
        x2_z = CPBMForLWAdapter.white_lindsey_boundary_stratum_x2_block(
            pair_unit_coefficients;
            axis = :z,
            parent_axis_counts = (7, 7, 7),
            overlap_1d = factors.overlap_1d,
            x2_1d = factors.x2_1d,
        ),
        kinetic = CPBMForLWAdapter.white_lindsey_boundary_stratum_kinetic_block(
            pair_unit_coefficients;
            parent_axis_counts = (7, 7, 7),
            overlap_1d = factors.overlap_1d,
            kinetic_1d = factors.kinetic_1d,
        ),
    )
end

function _lw_adapter_seed_operator_slice(
    seed,
    left_range::UnitRange{Int},
    right_range::UnitRange{Int},
    term::Symbol,
)
    matrix = getproperty(seed.fixed_block, term)
    return matrix[left_range, right_range]
end

function _lw_adapter_seed_operator_slice(seed, seed_local_facts, term::Symbol)
    seed_local_facts.status === :available_seed_local_pair_slice || throw(
        ArgumentError("White-Lindsey seed local facts must expose a local pair slice"),
    )
    return _lw_adapter_seed_operator_slice(
        seed,
        seed_local_facts.face_global_range,
        seed_local_facts.edge_global_range,
        term,
    )
end

function _lw_adapter_oracle_value_comparison(
    adapter_block,
    oracle_summary;
    seed_local_facts = nothing,
    seed_operator_slice = nothing,
    tolerance::Float64 = 1.0e-10,
)
    term = adapter_block.term
    expected_adapter_shape = (
        adapter_block.metadata.left_retained_column_count,
        adapter_block.metadata.right_retained_column_count,
    )
    adapter_shape = size(adapter_block.block)
    oracle_shape_available =
        hasproperty(oracle_summary.fixed_block_operator_matrix_sizes, term)
    oracle_term_available =
        hasproperty(oracle_summary.one_body_operator_matrix_available, term) &&
        getproperty(oracle_summary.one_body_operator_matrix_available, term)
    oracle_shape = oracle_shape_available ?
                   getproperty(oracle_summary.fixed_block_operator_matrix_sizes, term) :
                   nothing
    oracle_shape_status =
        oracle_summary.overlap_ready && oracle_term_available && oracle_shape_available ?
        Symbol("available_global_retained_", String(term), "_shape") :
        Symbol("blocked_missing_global_retained_", String(term), "_shape")

    shape_status =
        adapter_shape == expected_adapter_shape &&
        oracle_term_available &&
        oracle_shape_available ?
        :local_pair_shape_matches_adapter_metadata_global_oracle_shape_available :
        :blocked_shape_metadata_mismatch
    seed_slice_status =
        isnothing(seed_local_facts) ?
        :not_requested :
        seed_local_facts.status
    seed_slice_shape =
        isnothing(seed_operator_slice) ? nothing : size(seed_operator_slice)
    max_abs_error =
        isnothing(seed_operator_slice) ?
        nothing :
        maximum(abs.(adapter_block.block .- seed_operator_slice))
    max_abs_error_status =
        isnothing(max_abs_error) ?
        :not_compared_local_seed_pair_block_not_available :
        max_abs_error <= tolerance ?
        :compared_local_seed_pair_slice_within_tolerance :
        :compared_local_seed_pair_slice_exceeds_tolerance
    status =
        shape_status !==
        :local_pair_shape_matches_adapter_metadata_global_oracle_shape_available ?
        :blocked_oracle_shape_comparison :
        max_abs_error_status ===
        :compared_local_seed_pair_slice_within_tolerance ?
        :value_compared :
        isnothing(max_abs_error) ?
        :metadata_shape_only :
        :blocked_oracle_value_comparison
    blocker =
        status === :value_compared || status === :metadata_shape_only ?
        nothing :
        status === :blocked_oracle_value_comparison ?
        :local_seed_pair_slice_mismatch :
        :shape_metadata_mismatch
    symmetry_error =
        adapter_shape[1] == adapter_shape[2] ?
        maximum(abs.(adapter_block.block .- transpose(adapter_block.block))) :
        nothing
    symmetry_error_status =
        isnothing(symmetry_error) ?
        :not_applicable_rectangular_local_pair_block :
        symmetry_error <= tolerance ?
        :square_local_pair_symmetric_within_tolerance :
        :square_local_pair_symmetry_exceeds_tolerance

    return (;
        object_kind = :white_lindsey_adapter_oracle_comparison,
        term,
        pair_key = adapter_block.pair_key,
        pair_family = _lw_adapter_oracle_pair_family(adapter_block),
        adapter_materialization_path = adapter_block.metadata.materialization_path,
        adapter_shape,
        expected_adapter_shape,
        oracle_shape,
        oracle_shape_status,
        shape_status,
        seed_slice_status,
        seed_slice_shape,
        tolerance,
        max_abs_error,
        max_abs_error_status,
        symmetry_error,
        symmetry_error_status,
        status,
        blocker,
        oracle_role = oracle_summary.oracle_role,
        route_authority = false,
        adapter_authority = false,
        local_pair_block_materialized =
            adapter_block.metadata.local_pair_block_materialized,
        source_operator_blocks_materialized =
            adapter_block.metadata.source_operator_blocks_materialized,
        final_pair_blocks_materialized =
            adapter_block.metadata.final_pair_blocks_materialized,
        operator_blocks_materialized =
            adapter_block.metadata.operator_blocks_materialized,
        hamiltonian_data_materialized =
            adapter_block.metadata.hamiltonian_data_materialized,
        artifacts_materialized = adapter_block.metadata.artifacts_materialized,
        dense_parent_parent_overlap_materialized =
            adapter_block.metadata.dense_parent_parent_overlap_materialized,
        kinetic_factor_form =
            get(adapter_block.metadata, :kinetic_factor_form, nothing),
        kinetic_component_axes =
            get(adapter_block.metadata, :kinetic_component_axes, nothing),
        kinetic_component_terms =
            get(adapter_block.metadata, :kinetic_component_terms, nothing),
    )
end

function _lw_adapter_oracle_value_comparisons(
    adapter_blocks,
    oracle_summary,
    seed,
    seed_local_facts,
    comparison_terms,
    left_range::UnitRange{Int},
    right_range::UnitRange{Int},
)
    return Tuple(
        _lw_adapter_oracle_value_comparison(
            getproperty(adapter_blocks, term),
            oracle_summary;
            seed_local_facts,
            seed_operator_slice =
                _lw_adapter_seed_operator_slice(seed, left_range, right_range, term),
        )
        for term in comparison_terms
    )
end

function _lw_adapter_oracle_comparison_summary(comparisons)
    comparison_tuple = Tuple(comparisons)
    return (;
        object_kind = :white_lindsey_adapter_oracle_comparison_summary,
        comparison_count = length(comparison_tuple),
        terms = Tuple(comparison.term for comparison in comparison_tuple),
        pair_families =
            Tuple(comparison.pair_family for comparison in comparison_tuple),
        statuses = Tuple(comparison.status for comparison in comparison_tuple),
        metadata_shape_only_count =
            count(comparison -> comparison.status === :metadata_shape_only,
                comparison_tuple),
        value_comparison_count =
            count(comparison -> comparison.max_abs_error !== nothing,
                comparison_tuple),
        blocked_count =
            count(comparison -> !isnothing(comparison.blocker),
                comparison_tuple),
    )
end

function _lw_adapter_unique_tuple(values)
    unique_values = Any[]
    for value in values
        any(existing -> existing == value, unique_values) && continue
        push!(unique_values, value)
    end
    return Tuple(unique_values)
end

function _lw_adapter_oracle_validation_coverage_summary(comparisons)
    comparison_tuple = Tuple(comparisons)
    pair_families = _lw_adapter_unique_tuple(
        Tuple(comparison.pair_family for comparison in comparison_tuple),
    )
    terms = _lw_adapter_unique_tuple(
        Tuple(comparison.term for comparison in comparison_tuple),
    )
    max_abs_errors = Tuple(
        comparison.max_abs_error for comparison in comparison_tuple
        if !isnothing(comparison.max_abs_error)
    )
    return (;
        object_kind =
            :white_lindsey_adapter_oracle_validation_coverage_summary,
        pair_family_count = length(pair_families),
        pair_families,
        term_count = length(terms),
        terms,
        comparison_count = length(comparison_tuple),
        value_comparison_count =
            count(comparison -> !isnothing(comparison.max_abs_error),
                comparison_tuple),
        blocked_count =
            count(comparison -> !isnothing(comparison.blocker),
                comparison_tuple),
        metadata_shape_only_count =
            count(comparison -> comparison.status === :metadata_shape_only,
                comparison_tuple),
        max_abs_error =
            isempty(max_abs_errors) ? nothing : maximum(max_abs_errors),
        all_within_tolerance = all(
            comparison -> (
                comparison.status === :value_compared &&
                !isnothing(comparison.max_abs_error) &&
                comparison.max_abs_error <= comparison.tolerance
            ),
            comparison_tuple,
        ),
        old_seed_validation_oracle_only =
            all(comparison -> comparison.oracle_role === :validation_oracle_only,
                comparison_tuple),
        route_authority = any(comparison -> comparison.route_authority,
            comparison_tuple),
        adapter_authority = any(comparison -> comparison.adapter_authority,
            comparison_tuple),
        local_pair_blocks_materialized =
            all(comparison -> comparison.local_pair_block_materialized,
                comparison_tuple),
        source_operator_blocks_materialized =
            all(comparison -> comparison.source_operator_blocks_materialized,
                comparison_tuple),
        final_pair_blocks_materialized =
            all(comparison -> comparison.final_pair_blocks_materialized,
                comparison_tuple),
        operator_blocks_materialized =
            any(comparison -> comparison.operator_blocks_materialized,
                comparison_tuple),
        hamiltonian_data_materialized =
            any(comparison -> comparison.hamiltonian_data_materialized,
                comparison_tuple),
        exports_materialized = false,
        artifacts_materialized =
            any(comparison -> comparison.artifacts_materialized,
                comparison_tuple),
        dense_parent_parent_overlap_materialized =
            any(
                comparison -> comparison.dense_parent_parent_overlap_materialized,
                comparison_tuple,
            ),
    )
end

@testset "CartesianPairBlockMaterialization White-Lindsey oracle comparison scaffold" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    seed_report =
        GaussletBases._white_lindsey_low_order_materialized_seed_report(;
            expansion,
        )
    seed = seed_report.fixture
    oracle_summary =
        CPBMForLWAdapter.white_lindsey_materialized_seed_oracle_summary(
            seed_report,
    )
    seed_local_facts = _lw_adapter_seed_local_facet_edge_facts(seed_report)
    doside_source_1d =
        _lw_adapter_seed_compatible_doside_source_1d(seed; expansion)
    factors = _lw_adapter_seed_one_body_factors(doside_source_1d)
    real_units = _lw_adapter_real_units(
        doside_source_1d;
        prefix = "lw_oracle_comparison",
    )
    real_corner_unit = _lw_adapter_real_corner_unit(
        doside_source_1d;
        prefix = "lw_oracle_comparison",
    )
    real_pair_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_pair_unit_coefficients(
            _lw_adapter_unit_pair(
                real_units.real_facet_unit,
                real_units.real_edge_unit,
                1,
            ),
        )
    facet_edge_adapter_blocks =
        _lw_adapter_seed_one_body_blocks(real_pair_coefficients, factors)
    facet_facet_pair_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_pair_unit_coefficients(
            _lw_adapter_unit_pair(
                real_units.real_facet_unit,
                real_units.real_facet_unit,
                2,
            ),
        )
    facet_facet_adapter_blocks =
        _lw_adapter_seed_one_body_blocks(facet_facet_pair_coefficients, factors)
    edge_edge_pair_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_pair_unit_coefficients(
            _lw_adapter_unit_pair(
                real_units.real_edge_unit,
                real_units.real_edge_unit,
                3,
            ),
        )
    edge_edge_adapter_blocks =
        _lw_adapter_seed_one_body_blocks(edge_edge_pair_coefficients, factors)
    edge_corner_pair_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_pair_unit_coefficients(
            _lw_adapter_unit_pair(
                real_units.real_edge_unit,
                real_corner_unit,
                4,
            ),
        )
    edge_corner_adapter_blocks =
        _lw_adapter_seed_one_body_blocks(edge_corner_pair_coefficients, factors)
    corner_corner_pair_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_pair_unit_coefficients(
            _lw_adapter_unit_pair(
                real_corner_unit,
                real_corner_unit,
                5,
            ),
        )
    corner_corner_adapter_blocks =
        _lw_adapter_seed_one_body_blocks(corner_corner_pair_coefficients, factors)

    @test seed_local_facts.object_kind ==
          :white_lindsey_seed_local_pair_facts
    @test seed_local_facts.status == :available_seed_local_pair_slice
    @test isnothing(seed_local_facts.blocker)
    @test seed_local_facts.route_unit_range_status ==
          :grouped_route_ranges_available_not_local_slice_authority
    @test seed_local_facts.per_piece_range_status ==
          :per_piece_seed_inventory_ranges_available
    @test seed_local_facts.per_piece_range_source ==
          :seed_inventory_retained_ranges
    @test seed_local_facts.face_index == 5
    @test seed_local_facts.edge_index == 11
    @test seed_local_facts.face_signature == (;
        face_kind = :yz,
        fixed_axis = :x,
        fixed_side = :low,
        fixed_index = 1,
    )
    @test seed_local_facts.edge_signature == (;
        free_axis = :z,
        fixed_axes = (:x, :y),
        fixed_sides = (:high, :low),
        fixed_indices = (7, 1),
    )
    @test seed_local_facts.corner_index == 6
    @test seed_local_facts.corner_signature == (;
        fixed_sides = (:high, :low, :high),
        fixed_indices = (7, 1, 7),
    )
    @test seed_local_facts.face_global_range == 162:170
    @test seed_local_facts.edge_global_range == 210:212
    @test seed_local_facts.corner_global_range == 221:221
    @test seed_local_facts.face_local_range == 37:45
    @test seed_local_facts.edge_local_range == 85:87
    @test seed_local_facts.corner_local_range == 96:96
    @test seed_local_facts.grouped_face_range == 126:179
    @test seed_local_facts.grouped_edge_range == 180:215
    @test seed_local_facts.grouped_corner_range == 216:223
    @test seed_local_facts.slice_shape == (9, 3)
    @test seed_local_facts.edge_corner_slice_shape == (3, 1)
    @test seed_local_facts.overlap_slice_available
    @test seed_local_facts.validation_oracle_only
    @test !seed_local_facts.route_authority
    @test !seed_local_facts.adapter_authority

    comparison_terms = (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    adapter_paths = (;
        overlap = :white_lindsey_boundary_stratum_overlap_adapter,
        position_x = :white_lindsey_boundary_stratum_position_adapter,
        position_y = :white_lindsey_boundary_stratum_position_adapter,
        position_z = :white_lindsey_boundary_stratum_position_adapter,
        x2_x = :white_lindsey_boundary_stratum_x2_adapter,
        x2_y = :white_lindsey_boundary_stratum_x2_adapter,
        x2_z = :white_lindsey_boundary_stratum_x2_adapter,
        kinetic = :white_lindsey_boundary_stratum_kinetic_adapter,
    )
    facet_edge_comparisons = _lw_adapter_oracle_value_comparisons(
        facet_edge_adapter_blocks,
        oracle_summary,
        seed,
        seed_local_facts,
        comparison_terms,
        seed_local_facts.face_global_range,
        seed_local_facts.edge_global_range,
    )
    facet_facet_comparisons = _lw_adapter_oracle_value_comparisons(
        facet_facet_adapter_blocks,
        oracle_summary,
        seed,
        seed_local_facts,
        comparison_terms,
        seed_local_facts.face_global_range,
        seed_local_facts.face_global_range,
    )
    edge_edge_comparisons = _lw_adapter_oracle_value_comparisons(
        edge_edge_adapter_blocks,
        oracle_summary,
        seed,
        seed_local_facts,
        comparison_terms,
        seed_local_facts.edge_global_range,
        seed_local_facts.edge_global_range,
    )
    @test edge_corner_pair_coefficients.status ==
          :materialized_white_lindsey_pair_unit_coefficients
    @test edge_corner_pair_coefficients.left_stratum_kind == :edge_cpb
    @test edge_corner_pair_coefficients.right_stratum_kind == :corner_cpb
    @test edge_corner_pair_coefficients.left_retained_column_count == 3
    @test edge_corner_pair_coefficients.right_retained_column_count == 1
    @test edge_corner_pair_coefficients.left_support_indices !== nothing
    @test edge_corner_pair_coefficients.right_support_indices ==
          [GaussletBases._cartesian_flat_index(7, 1, 7, (7, 7, 7))]
    edge_corner_comparisons = _lw_adapter_oracle_value_comparisons(
        edge_corner_adapter_blocks,
        oracle_summary,
        seed,
        seed_local_facts,
        comparison_terms,
        seed_local_facts.edge_global_range,
        seed_local_facts.corner_global_range,
    )
    corner_corner_comparisons = _lw_adapter_oracle_value_comparisons(
        corner_corner_adapter_blocks,
        oracle_summary,
        seed,
        seed_local_facts,
        comparison_terms,
        seed_local_facts.corner_global_range,
        seed_local_facts.corner_global_range,
    )
    comparison_batches = (
        (;
            comparisons = facet_edge_comparisons,
            expected_pair_key = (
                :lw_oracle_comparison_real_facet_unit,
                :lw_oracle_comparison_real_edge_unit,
            ),
            expected_pair_family = :facet_cpb__edge_cpb,
            expected_shape = (9, 3),
            expected_symmetry_status =
                :not_applicable_rectangular_local_pair_block,
        ),
        (;
            comparisons = facet_facet_comparisons,
            expected_pair_key = (
                :lw_oracle_comparison_real_facet_unit,
                :lw_oracle_comparison_real_facet_unit,
            ),
            expected_pair_family = :facet_cpb__facet_cpb,
            expected_shape = (9, 9),
            expected_symmetry_status =
                :square_local_pair_symmetric_within_tolerance,
        ),
        (;
            comparisons = edge_edge_comparisons,
            expected_pair_key = (
                :lw_oracle_comparison_real_edge_unit,
                :lw_oracle_comparison_real_edge_unit,
            ),
            expected_pair_family = :edge_cpb__edge_cpb,
            expected_shape = (3, 3),
            expected_symmetry_status =
                :square_local_pair_symmetric_within_tolerance,
        ),
        (;
            comparisons = edge_corner_comparisons,
            expected_pair_key = (
                :lw_oracle_comparison_real_edge_unit,
                :lw_oracle_comparison_real_corner_unit,
            ),
            expected_pair_family = :edge_cpb__corner_cpb,
            expected_shape = (3, 1),
            expected_symmetry_status =
                :not_applicable_rectangular_local_pair_block,
        ),
        (;
            comparisons = corner_corner_comparisons,
            expected_pair_key = (
                :lw_oracle_comparison_real_corner_unit,
                :lw_oracle_comparison_real_corner_unit,
            ),
            expected_pair_family = :corner_cpb__corner_cpb,
            expected_shape = (1, 1),
            expected_symmetry_status =
                :square_local_pair_symmetric_within_tolerance,
        ),
    )

    coverage_comparisons = Tuple(
        comparison for batch in comparison_batches
        for comparison in batch.comparisons
    )
    coverage_summary = _lw_adapter_oracle_validation_coverage_summary(
        coverage_comparisons,
    )
    @test coverage_summary.object_kind ==
          :white_lindsey_adapter_oracle_validation_coverage_summary
    @test coverage_summary.pair_family_count == 5
    @test coverage_summary.pair_families == (
        :facet_cpb__edge_cpb,
        :facet_cpb__facet_cpb,
        :edge_cpb__edge_cpb,
        :edge_cpb__corner_cpb,
        :corner_cpb__corner_cpb,
    )
    @test coverage_summary.term_count == length(comparison_terms)
    @test coverage_summary.terms == comparison_terms
    @test coverage_summary.comparison_count == 5 * length(comparison_terms)
    @test coverage_summary.value_comparison_count ==
          5 * length(comparison_terms)
    @test coverage_summary.blocked_count == 0
    @test coverage_summary.metadata_shape_only_count == 0
    @test coverage_summary.max_abs_error <= 1.0e-10
    @test coverage_summary.all_within_tolerance
    @test coverage_summary.old_seed_validation_oracle_only
    @test !coverage_summary.route_authority
    @test !coverage_summary.adapter_authority
    @test coverage_summary.local_pair_blocks_materialized
    @test coverage_summary.source_operator_blocks_materialized
    @test coverage_summary.final_pair_blocks_materialized
    @test !coverage_summary.operator_blocks_materialized
    @test !coverage_summary.hamiltonian_data_materialized
    @test !coverage_summary.exports_materialized
    @test !coverage_summary.artifacts_materialized
    @test !coverage_summary.dense_parent_parent_overlap_materialized

    for batch in comparison_batches, comparison in batch.comparisons
        @test comparison.object_kind ==
              :white_lindsey_adapter_oracle_comparison
        @test comparison.pair_key == batch.expected_pair_key
        @test comparison.pair_family == batch.expected_pair_family
        @test comparison.adapter_materialization_path ==
              getproperty(adapter_paths, comparison.term)
        @test comparison.adapter_shape == batch.expected_shape
        @test comparison.expected_adapter_shape == batch.expected_shape
        @test comparison.oracle_shape == (223, 223)
        @test comparison.oracle_shape_status ==
              Symbol("available_global_retained_", String(comparison.term), "_shape")
        @test comparison.shape_status ==
              :local_pair_shape_matches_adapter_metadata_global_oracle_shape_available
        @test comparison.seed_slice_status == :available_seed_local_pair_slice
        @test comparison.seed_slice_shape == batch.expected_shape
        @test comparison.status == :value_compared
        @test isnothing(comparison.blocker)
        @test comparison.max_abs_error <= comparison.tolerance
        @test comparison.max_abs_error_status ==
              :compared_local_seed_pair_slice_within_tolerance
        @test comparison.symmetry_error_status == batch.expected_symmetry_status
        if batch.expected_shape[1] == batch.expected_shape[2]
            @test comparison.symmetry_error <= comparison.tolerance
        else
            @test isnothing(comparison.symmetry_error)
        end
        @test comparison.oracle_role == :validation_oracle_only
        @test !comparison.route_authority
        @test !comparison.adapter_authority
        @test comparison.local_pair_block_materialized
        @test comparison.source_operator_blocks_materialized
        @test comparison.final_pair_blocks_materialized
        @test !comparison.operator_blocks_materialized
        @test !comparison.hamiltonian_data_materialized
        @test !comparison.artifacts_materialized
        @test !comparison.dense_parent_parent_overlap_materialized
    end
    kinetic_comparison = facet_edge_comparisons[end]
    @test kinetic_comparison.term == :kinetic
    @test kinetic_comparison.kinetic_factor_form ==
          :factorized_cartesian_sum_kss_sks_ssk
    @test kinetic_comparison.kinetic_component_axes == (:x, :y, :z)
    @test kinetic_comparison.kinetic_component_terms == (
        :kinetic_x_component,
        :kinetic_y_component,
        :kinetic_z_component,
    )

    facet_edge_summary = _lw_adapter_oracle_comparison_summary(
        facet_edge_comparisons,
    )
    @test facet_edge_summary.object_kind ==
          :white_lindsey_adapter_oracle_comparison_summary
    @test facet_edge_summary.comparison_count == length(comparison_terms)
    @test facet_edge_summary.terms == comparison_terms
    @test facet_edge_summary.pair_families ==
          ntuple(_ -> :facet_cpb__edge_cpb, length(comparison_terms))
    @test facet_edge_summary.statuses ==
          ntuple(_ -> :value_compared, length(comparison_terms))
    @test facet_edge_summary.metadata_shape_only_count == 0
    @test facet_edge_summary.value_comparison_count == length(comparison_terms)
    @test facet_edge_summary.blocked_count == 0

    facet_facet_summary = _lw_adapter_oracle_comparison_summary(
        facet_facet_comparisons,
    )
    @test facet_facet_summary.object_kind ==
          :white_lindsey_adapter_oracle_comparison_summary
    @test facet_facet_summary.comparison_count == length(comparison_terms)
    @test facet_facet_summary.terms == comparison_terms
    @test facet_facet_summary.pair_families ==
          ntuple(_ -> :facet_cpb__facet_cpb, length(comparison_terms))
    @test facet_facet_summary.statuses ==
          ntuple(_ -> :value_compared, length(comparison_terms))
    @test facet_facet_summary.metadata_shape_only_count == 0
    @test facet_facet_summary.value_comparison_count == length(comparison_terms)
    @test facet_facet_summary.blocked_count == 0

    edge_edge_summary = _lw_adapter_oracle_comparison_summary(
        edge_edge_comparisons,
    )
    @test edge_edge_summary.object_kind ==
          :white_lindsey_adapter_oracle_comparison_summary
    @test edge_edge_summary.comparison_count == length(comparison_terms)
    @test edge_edge_summary.terms == comparison_terms
    @test edge_edge_summary.pair_families ==
          ntuple(_ -> :edge_cpb__edge_cpb, length(comparison_terms))
    @test edge_edge_summary.statuses ==
          ntuple(_ -> :value_compared, length(comparison_terms))
    @test edge_edge_summary.metadata_shape_only_count == 0
    @test edge_edge_summary.value_comparison_count == length(comparison_terms)
    @test edge_edge_summary.blocked_count == 0

    edge_corner_summary = _lw_adapter_oracle_comparison_summary(
        edge_corner_comparisons,
    )
    @test edge_corner_summary.object_kind ==
          :white_lindsey_adapter_oracle_comparison_summary
    @test edge_corner_summary.comparison_count == length(comparison_terms)
    @test edge_corner_summary.terms == comparison_terms
    @test edge_corner_summary.pair_families ==
          ntuple(_ -> :edge_cpb__corner_cpb, length(comparison_terms))
    @test edge_corner_summary.statuses ==
          ntuple(_ -> :value_compared, length(comparison_terms))
    @test edge_corner_summary.metadata_shape_only_count == 0
    @test edge_corner_summary.value_comparison_count == length(comparison_terms)
    @test edge_corner_summary.blocked_count == 0
    for comparison in edge_corner_comparisons
        @test comparison.seed_slice_shape == seed_local_facts.edge_corner_slice_shape
        @test comparison.max_abs_error <= comparison.tolerance
    end

    corner_corner_summary = _lw_adapter_oracle_comparison_summary(
        corner_corner_comparisons,
    )
    @test corner_corner_summary.object_kind ==
          :white_lindsey_adapter_oracle_comparison_summary
    @test corner_corner_summary.comparison_count == length(comparison_terms)
    @test corner_corner_summary.terms == comparison_terms
    @test corner_corner_summary.pair_families ==
          ntuple(_ -> :corner_cpb__corner_cpb, length(comparison_terms))
    @test corner_corner_summary.statuses ==
          ntuple(_ -> :value_compared, length(comparison_terms))
    @test corner_corner_summary.metadata_shape_only_count == 0
    @test corner_corner_summary.value_comparison_count == length(comparison_terms)
    @test corner_corner_summary.blocked_count == 0
end
