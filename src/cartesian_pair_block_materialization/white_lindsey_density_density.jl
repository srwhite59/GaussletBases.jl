# Decomposed White--Lindsey density-density electron-electron route.
#
# This assembles the full retained two-index interaction matrix needed by the
# current IDA/HF convention from decomposed WL unit pairs. The Gaussian expansion
# index stays inside each support-block contraction. This path does not build a
# full parent CPB/window, call ordinary Cartesian IDA operators, assemble a
# Hamiltonian, run SCF, add GTO supplements, or apply PQS transforms.

function route_global_decomposed_wl_density_density_matrix(
    source;
    parent_axis_counts = nothing,
    parent_axis_bundle_object = nothing,
    coulomb_expansion = nothing,
    metadata = (;),
)
    inventory = _route_global_electron_nuclear_by_center_inventory(
        source;
        parent_axis_counts,
        parent_axis_bundle_object,
    )
    if inventory.status !== :available_white_lindsey_decomposed_unit_pair_inventory
        return _route_global_density_density_result(
            :blocked_route_global_density_density_interaction_matrix,
            inventory.blocker,
            inventory,
            nothing,
            0;
            metadata,
        )
    end

    blocker = _white_lindsey_density_density_input_blocker(
        parent_axis_counts,
        parent_axis_bundle_object,
        coulomb_expansion,
    )
    if !isnothing(blocker)
        return _route_global_density_density_result(
            :blocked_route_global_density_density_interaction_matrix,
            blocker,
            inventory,
            nothing,
            0;
            metadata,
        )
    end

    axis_counts = _axis_counts_tuple(parent_axis_counts)
    pair_factor_terms =
        _white_lindsey_density_density_pair_factor_terms(parent_axis_bundle_object)
    coefficients =
        _white_lindsey_density_density_expansion_coefficients(coulomb_expansion)
    retained_dimension = Int(inventory.retained_dimension)
    matrix = zeros(Float64, retained_dimension, retained_dimension)
    if inventory.unit_pairs isa CUP.UnitPairIndexTable
        stream_state = @timeg "decomposed_wl.density_density.local_pair_stream" begin
            _route_global_density_density_streaming_fill!(
                matrix,
                inventory.unit_pairs,
                axis_counts,
                pair_factor_terms,
                coefficients,
            )
        end
        if !isnothing(stream_state.blocker)
            return _route_global_density_density_result(
                :blocked_route_global_density_density_interaction_matrix,
                stream_state.blocker,
                inventory,
                nothing,
                stream_state.placed_block_count;
                metadata = merge(
                    NamedTuple(metadata),
                    (;
                        blocked_pair_key = stream_state.blocked_pair_key,
                        materialization_path = stream_state.materialization_path,
                    ),
                ),
            )
        end
        return _route_global_density_density_result(
            :materialized_route_global_density_density_interaction_matrix,
            nothing,
            inventory,
            matrix,
            stream_state.placed_block_count;
            gaussian_term_count = length(coefficients),
            unit_coefficient_cache_entry_count =
                stream_state.unit_coefficient_cache_entry_count,
            pair_factor_term_shapes =
                _white_lindsey_density_density_pair_factor_term_shapes(
                    pair_factor_terms,
                ),
            metadata = merge(
                NamedTuple(metadata),
                (;
                    materialization_path = stream_state.materialization_path,
                    pair_input_kind = :unit_pair_index_table,
                    local_pair_blocks_collected = false,
                    pair_coefficients_batch_materialized = false,
                    streaming_retained_matrix_insertion = true,
                ),
            ),
        )
    end

    placed_block_count = 0
    pair_coefficients_batch =
        white_lindsey_boundary_stratum_pair_unit_coefficients(inventory.unit_pairs)
    if pair_coefficients_batch.status !==
       :materialized_white_lindsey_pair_unit_coefficients_batch
        return _route_global_density_density_result(
            :blocked_route_global_density_density_interaction_matrix,
            pair_coefficients_batch.blocker,
            inventory,
            nothing,
            placed_block_count;
            metadata,
        )
    end

    for (unit_pair, pair_coefficients) in
        zip(inventory.unit_pairs, pair_coefficients_batch.results)
        block = _white_lindsey_density_density_pair_block(
            pair_coefficients,
            axis_counts,
            pair_factor_terms,
            coefficients,
        )
        rows = unit_pair.left_unit.column_range
        columns = unit_pair.right_unit.column_range
        size(block) == (length(rows), length(columns)) ||
            return _route_global_density_density_result(
                :blocked_route_global_density_density_interaction_matrix,
                :density_density_pair_block_shape_mismatch,
                inventory,
                nothing,
                placed_block_count;
                metadata = merge(
                    NamedTuple(metadata),
                    (; blocked_pair_key = unit_pair.pair_key),
                ),
            )

        matrix[rows, columns] .+= block
        placed_block_count += 1
        rows == columns && continue
        matrix[columns, rows] .+= transpose(block)
        placed_block_count += 1
    end

    return _route_global_density_density_result(
        :materialized_route_global_density_density_interaction_matrix,
        nothing,
        inventory,
        matrix,
        placed_block_count;
        gaussian_term_count = length(coefficients),
        unit_coefficient_cache_entry_count =
            pair_coefficients_batch.unit_coefficient_cache_entry_count,
        pair_factor_term_shapes =
            _white_lindsey_density_density_pair_factor_term_shapes(pair_factor_terms),
        metadata,
    )
end

function _route_global_density_density_streaming_fill!(
    matrix::AbstractMatrix{Float64},
    unit_pairs::CUP.UnitPairIndexTable,
    axis_counts::NTuple{3,Int},
    pair_factor_terms,
    coefficients::Vector{Float64},
)
    cache = Dict{Symbol,Any}()
    materialized_count = 0
    skipped_count = 0
    placed_block_count = 0
    blocker = nothing
    blocked_pair_key = nothing
    dimension = size(matrix, 1)
    for unit_pair in unit_pairs
        pair_coefficients = _white_lindsey_pair_unit_coefficients_result(
            unit_pair,
            _white_lindsey_unit_coefficients_from_local_cache(
                cache,
                unit_pair.left_unit,
            ),
            _white_lindsey_unit_coefficients_from_local_cache(
                cache,
                unit_pair.right_unit,
            ),
            length(cache),
        )
        if !_is_ready_white_lindsey_pair_unit_coefficients(pair_coefficients)
            skipped_count += 1
            if isnothing(blocker)
                blocker = pair_coefficients.blocker
                blocked_pair_key = pair_coefficients.pair_key
            end
            continue
        end

        block = _white_lindsey_density_density_pair_block(
            pair_coefficients,
            axis_counts,
            pair_factor_terms,
            coefficients,
        )
        rows = unit_pair.left_unit.column_range
        columns = unit_pair.right_unit.column_range
        placement_blocker =
            _route_global_density_density_streaming_placement_blocker(
                block,
                rows,
                columns,
                dimension,
            )
        if !isnothing(placement_blocker)
            skipped_count += 1
            if isnothing(blocker)
                blocker = placement_blocker
                blocked_pair_key = unit_pair.pair_key
            end
            continue
        end

        matrix[rows, columns] .+= block
        placed_block_count += 1
        if rows != columns
            matrix[columns, rows] .+= transpose(block)
            placed_block_count += 1
        end
        materialized_count += 1
    end
    return (;
        materialized_count,
        skipped_count,
        placed_block_count,
        blocker,
        blocked_pair_key,
        unit_coefficient_cache_entry_count = length(cache),
        pair_count = length(unit_pairs),
        materialization_path =
            :white_lindsey_decomposed_wl_streaming_density_density_global_assembly,
    )
end

function _route_global_density_density_streaming_placement_blocker(
    block,
    rows,
    columns,
    dimension::Int,
)
    block isa AbstractMatrix{<:Real} || return :missing_density_density_pair_block
    (
        _one_body_placement_valid_column_range(rows) &&
        _one_body_placement_valid_column_range(columns)
    ) || return :missing_column_ranges
    size(block) == (length(rows), length(columns)) ||
        return :density_density_pair_block_shape_mismatch
    _one_body_placement_ranges_inside_global_dimension(rows, columns, dimension) ||
        return :target_ranges_outside_global_dimension
    return nothing
end

function _white_lindsey_density_density_input_blocker(
    parent_axis_counts,
    parent_axis_bundle_object,
    coulomb_expansion,
)
    isnothing(parent_axis_counts) && return :missing_parent_axis_counts
    isnothing(parent_axis_bundle_object) && return :missing_parent_axis_bundle_object
    isnothing(coulomb_expansion) && return :missing_coulomb_gaussian_expansion

    coefficients = _white_lindsey_descriptor_property(
        coulomb_expansion,
        :coefficients,
    )
    coefficients isa AbstractVector ||
        return :missing_coulomb_expansion_coefficients

    axis_counts = _axis_counts_tuple(parent_axis_counts)
    pair_factor_terms =
        _white_lindsey_density_density_pair_factor_terms(parent_axis_bundle_object)
    for (axis_index, axis) in pairs((:x, :y, :z))
        terms = getproperty(pair_factor_terms, axis)
        terms isa AbstractArray{<:Real,3} ||
            return Symbol("missing_$(axis)_pair_factor_terms")
        size(terms, 1) == length(coefficients) ||
            return :pair_factor_term_count_mismatch
        size(terms, 2) == axis_counts[axis_index] ||
            return Symbol("$(axis)_pair_factor_terms_size_mismatch")
        size(terms, 3) == axis_counts[axis_index] ||
            return Symbol("$(axis)_pair_factor_terms_size_mismatch")
    end
    return nothing
end

function _white_lindsey_density_density_pair_factor_terms(parent_axis_bundle_object)
    return (;
        x = _white_lindsey_density_density_axis_pair_factor_terms(
            parent_axis_bundle_object,
            :x,
        ),
        y = _white_lindsey_density_density_axis_pair_factor_terms(
            parent_axis_bundle_object,
            :y,
        ),
        z = _white_lindsey_density_density_axis_pair_factor_terms(
            parent_axis_bundle_object,
            :z,
        ),
    )
end

function _white_lindsey_density_density_axis_pair_factor_terms(
    parent_axis_bundle_object,
    axis::Symbol,
)
    axis_bundle = _white_lindsey_descriptor_property(parent_axis_bundle_object, axis)
    pgdg_intermediate =
        _white_lindsey_descriptor_property(axis_bundle, :pgdg_intermediate)
    return _white_lindsey_descriptor_property(
        pgdg_intermediate,
        :pair_factor_terms,
    )
end

function _white_lindsey_density_density_expansion_coefficients(coulomb_expansion)
    coefficients = _white_lindsey_descriptor_property(
        coulomb_expansion,
        :coefficients,
    )
    return Float64[Float64(value) for value in coefficients]
end

function _white_lindsey_density_density_pair_block(
    pair_unit_coefficients,
    axis_counts::NTuple{3,Int},
    pair_factor_terms,
    coefficients::Vector{Float64},
)
    left_support = pair_unit_coefficients.left_support_indices
    right_support = pair_unit_coefficients.right_support_indices
    _assert_white_lindsey_overlap_support(left_support, :left)
    _assert_white_lindsey_overlap_support(right_support, :right)

    support_block = _white_lindsey_density_density_support_block(
        left_support,
        right_support,
        axis_counts,
        pair_factor_terms,
        coefficients,
    )
    left_support_coefficients = _white_lindsey_support_coefficient_matrix(
        pair_unit_coefficients.left_coefficient_matrix,
        left_support,
        pair_unit_coefficients.left_coefficient_space,
        :left,
    )
    right_support_coefficients = _white_lindsey_support_coefficient_matrix(
        pair_unit_coefficients.right_coefficient_matrix,
        right_support,
        pair_unit_coefficients.right_coefficient_space,
        :right,
    )
    return Matrix{Float64}(
        transpose(left_support_coefficients) *
        support_block *
        right_support_coefficients,
    )
end

function _white_lindsey_density_density_support_block(
    left_support,
    right_support,
    axis_counts::NTuple{3,Int},
    pair_factor_terms,
    coefficients::Vector{Float64},
)
    block = zeros(Float64, length(left_support), length(right_support))
    left_states =
        Tuple(
            ParentGaussletBases._cartesian_unflat_index(index, axis_counts)
            for index in left_support
        )
    right_states =
        Tuple(
            ParentGaussletBases._cartesian_unflat_index(index, axis_counts)
            for index in right_support
        )
    # The term-first PGDG pair-factor tables already match the WL interaction
    # convention. Do not apply integral weights here; IDA/HF consumes the full
    # retained interaction matrix at the density interpretation stage.
    for (row, left_state) in pairs(left_states)
        ix_left, iy_left, iz_left = left_state
        for (column, right_state) in pairs(right_states)
            ix_right, iy_right, iz_right = right_state
            value = 0.0
            @inbounds for term in eachindex(coefficients)
                value +=
                    coefficients[term] *
                    pair_factor_terms.x[term, ix_left, ix_right] *
                    pair_factor_terms.y[term, iy_left, iy_right] *
                    pair_factor_terms.z[term, iz_left, iz_right]
            end
            block[row, column] = value
        end
    end
    return block
end

function _white_lindsey_density_density_pair_factor_term_shapes(pair_factor_terms)
    return (;
        x = size(pair_factor_terms.x),
        y = size(pair_factor_terms.y),
        z = size(pair_factor_terms.z),
    )
end

function _route_global_density_density_result(
    status::Symbol,
    blocker,
    inventory,
    matrix,
    placed_block_count::Int;
    gaussian_term_count = 0,
    unit_coefficient_cache_entry_count = 0,
    pair_factor_term_shapes = :unavailable,
    metadata = (;),
)
    materialized =
        status === :materialized_route_global_density_density_interaction_matrix
    return (;
        object_kind =
            :white_lindsey_decomposed_route_global_density_density_interaction_matrix,
        term = :electron_electron_density_density,
        status,
        blocker,
        matrix,
        matrix_shape = materialized ? size(matrix) : :unavailable,
        retained_dimension =
            _route_global_one_body_value(inventory, :retained_dimension, nothing),
        retained_dimension_status =
            _route_global_one_body_value(
                inventory,
                :retained_dimension_status,
                :not_available,
            ),
        inventory_source_kind =
            _route_global_one_body_value(inventory, :source_kind, nothing),
        decomposed_unit_pair_inventory_status =
            _route_global_one_body_value(inventory, :status, nothing),
        unit_count = _route_global_one_body_value(inventory, :unit_count, 0),
        pair_count = _route_global_one_body_value(inventory, :pair_count, 0),
        placed_block_count,
        unit_coefficient_cache_entry_count,
        gaussian_term_count,
        pair_factor_source = :axis_pgdg_intermediate_pair_factor_terms,
        pair_factor_weighting = :wl_pair_factor_terms_existing_convention,
        pair_factor_term_shapes,
        axis_integral_weights_applied = false,
        axis_integral_weights_deferred = true,
        weight_application_stage = :ida_density_interpretation,
        ida_density_density_convention =
            :full_retained_two_index_interaction_matrix,
        route_global_density_density_matrix_materialized = materialized,
        electron_electron_density_density_matrix_materialized = materialized,
        local_operator_blocks_materialized = materialized,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        hamiltonian_data_materialized = false,
        gto_supplement_used = false,
        pqs_transforms_materialized = false,
        exports_or_artifacts = false,
        metadata = NamedTuple(metadata),
    )
end
