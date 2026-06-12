# Narrow route-global combined gausslet+GTO matrix assembly.
#
# This consumes existing decomposed WL gausslet route-global matrices and
# provider-level mixed/GTO blocks. Nuclear by-center blocks remain separated
# until this Hamiltonian assembly stage, where recorded nuclear charges are
# applied and centers are summed. No direct Cartesian fallback, full-window CPB
# route, PQS transform, export, or artifact is introduced here.

function route_global_combined_gto_one_electron_matrices(
    layout;
    overlap_result = nothing,
    kinetic_result = nothing,
    electron_nuclear_by_center_results = nothing,
    gto_bundle = nothing,
    mixed_gausslet_row_range = nothing,
    metadata = (;),
)
    blocker = _route_global_combined_gto_matrix_blocker(
        layout,
        overlap_result,
        kinetic_result,
        electron_nuclear_by_center_results,
        gto_bundle,
        mixed_gausslet_row_range,
    )
    if !isnothing(blocker)
        return _route_global_combined_gto_matrix_result(
            :blocked_route_global_combined_gto_one_electron_matrices,
            blocker,
            layout,
            nothing,
            nothing,
            ();
            gto_bundle,
            mixed_gausslet_row_range,
            metadata,
        )
    end

    overlap_matrix = _route_global_combined_gto_overlap_matrix(
        layout,
        _route_global_combined_gto_property(overlap_result, :matrix, nothing),
        _route_global_combined_gto_nested_property(
            gto_bundle,
            (:mixed_blocks, :overlap, :dense_block),
            nothing,
        ),
        _route_global_combined_gto_nested_property(
            gto_bundle,
            (:gto_blocks, :overlap, :dense_block),
            nothing,
        ),
    )
    kinetic_matrix = _route_global_combined_gto_overlap_matrix(
        layout,
        _route_global_combined_gto_property(kinetic_result, :matrix, nothing),
        _route_global_combined_gto_nested_property(
            gto_bundle,
            (:mixed_blocks, :kinetic, :dense_block),
            nothing,
        ),
        _route_global_combined_gto_nested_property(
            gto_bundle,
            (:gto_blocks, :kinetic, :dense_block),
            nothing,
        ),
    )
    nuclear_matrices =
        _route_global_combined_gto_nuclear_by_center_matrices(
            layout,
            electron_nuclear_by_center_results,
            gto_bundle,
        )
    hamiltonian = Matrix{Float64}(kinetic_matrix)
    for nuclear in nuclear_matrices
        hamiltonian .+= Float64(nuclear.nuclear_charge) .* nuclear.matrix
    end

    return _route_global_combined_gto_matrix_result(
        :materialized_route_global_combined_gto_one_electron_matrices,
        nothing,
        layout,
        overlap_matrix,
        hamiltonian,
        nuclear_matrices;
        gto_bundle,
        mixed_gausslet_row_range,
        metadata,
    )
end

function route_global_combined_gto_residual_moment_matrices(
    layout;
    moment_matrix_set = nothing,
    position_x_result = nothing,
    position_y_result = nothing,
    position_z_result = nothing,
    x2_x_result = nothing,
    x2_y_result = nothing,
    x2_z_result = nothing,
    gto_bundle = nothing,
    mixed_gausslet_row_range = nothing,
    metadata = (;),
)
    fields = (:position_x, :position_y, :position_z, :x2_x, :x2_y, :x2_z)
    gg_results = isnothing(moment_matrix_set) ?
                 (;
                     position_x = position_x_result,
                     position_y = position_y_result,
                     position_z = position_z_result,
                     x2_x = x2_x_result,
                     x2_y = x2_y_result,
                     x2_z = x2_z_result,
                 ) :
                 moment_matrix_set
    blocker = _route_global_combined_gto_moment_matrix_blocker(
        layout,
        gg_results,
        gto_bundle,
        mixed_gausslet_row_range,
        fields,
    )
    if !isnothing(blocker)
        return _route_global_combined_gto_moment_matrix_result(
            :blocked_route_global_combined_gto_residual_moment_matrices,
            blocker,
            layout,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing;
            gto_bundle,
            mixed_gausslet_row_range,
            metadata,
        )
    end

    position_x = _route_global_combined_gto_moment_matrix(layout, gg_results, gto_bundle, :position_x)
    position_y = _route_global_combined_gto_moment_matrix(layout, gg_results, gto_bundle, :position_y)
    position_z = _route_global_combined_gto_moment_matrix(layout, gg_results, gto_bundle, :position_z)
    x2_x = _route_global_combined_gto_moment_matrix(layout, gg_results, gto_bundle, :x2_x)
    x2_y = _route_global_combined_gto_moment_matrix(layout, gg_results, gto_bundle, :x2_y)
    x2_z = _route_global_combined_gto_moment_matrix(layout, gg_results, gto_bundle, :x2_z)

    return _route_global_combined_gto_moment_matrix_result(
        :materialized_route_global_combined_gto_residual_moment_matrices,
        nothing,
        layout,
        position_x,
        position_y,
        position_z,
        x2_x,
        x2_y,
        x2_z;
        gto_bundle,
        mixed_gausslet_row_range,
        metadata,
    )
end

function _route_global_combined_gto_moment_matrix_blocker(
    layout,
    gg_results,
    gto_bundle,
    mixed_gausslet_row_range,
    fields,
)
    _route_global_combined_gto_property(layout, :status, nothing) ===
        :available_route_global_combined_gto_basis_layout ||
        return :missing_route_global_combined_gto_basis_layout
    isnothing(gto_bundle) && return :missing_provider_level_gto_bundle

    gausslet_dimension = _route_global_combined_gto_property(
        layout,
        :gausslet_retained_dimension,
        nothing,
    )
    gausslet_dimension isa Integer && gausslet_dimension > 0 ||
        return :missing_gausslet_retained_dimension
    gto_count = _route_global_combined_gto_property(
        layout,
        :gto_supplement_orbital_count,
        nothing,
    )
    gto_count isa Integer && gto_count > 0 ||
        return :missing_gto_supplement_orbital_count

    first_mixed = _route_global_combined_gto_moment_bundle_block(
        gto_bundle,
        :mixed_blocks,
        first(fields),
    )
    first_mixed isa AbstractMatrix ||
        return Symbol("missing_mixed_gto_", String(first(fields)), "_block")
    mixed_range = _route_global_combined_gto_mixed_row_range(
        first_mixed,
        mixed_gausslet_row_range,
    )
    if mixed_range != _route_global_combined_gto_property(
        layout,
        :gausslet_retained_range,
        nothing,
    )
        return :missing_mixed_gto_route_global_row_coverage
    end

    for field in fields
        gg_matrix = _route_global_combined_gto_moment_gg_matrix(gg_results, field)
        gg_matrix isa AbstractMatrix ||
            return Symbol("missing_route_global_", String(field), "_matrix")
        size(gg_matrix) == (gausslet_dimension, gausslet_dimension) ||
            return Symbol("route_global_", String(field), "_matrix_dimension_mismatch")

        mixed_matrix =
            _route_global_combined_gto_moment_bundle_block(gto_bundle, :mixed_blocks, field)
        mixed_matrix isa AbstractMatrix ||
            return Symbol("missing_mixed_gto_", String(field), "_block")
        size(mixed_matrix) == (gausslet_dimension, gto_count) ||
            return Symbol("mixed_gto_", String(field), "_dimension_mismatch")

        gto_matrix =
            _route_global_combined_gto_moment_bundle_block(gto_bundle, :gto_blocks, field)
        gto_matrix isa AbstractMatrix ||
            return Symbol("missing_gto_gto_", String(field), "_block")
        size(gto_matrix) == (gto_count, gto_count) ||
            return Symbol("gto_gto_", String(field), "_dimension_mismatch")
    end
    return nothing
end

function _route_global_combined_gto_moment_matrix(
    layout,
    gg_results,
    gto_bundle,
    field::Symbol,
)
    return _route_global_combined_gto_overlap_matrix(
        layout,
        _route_global_combined_gto_moment_gg_matrix(gg_results, field),
        _route_global_combined_gto_moment_bundle_block(gto_bundle, :mixed_blocks, field),
        _route_global_combined_gto_moment_bundle_block(gto_bundle, :gto_blocks, field),
    )
end

function _route_global_combined_gto_moment_gg_matrix(gg_results, field::Symbol)
    result = getproperty(gg_results, field)
    result isa AbstractMatrix && return result
    return _route_global_combined_gto_property(result, :matrix, nothing)
end

function _route_global_combined_gto_moment_bundle_block(
    gto_bundle,
    block_group::Symbol,
    field::Symbol,
)
    return _route_global_combined_gto_nested_property(
        gto_bundle,
        (block_group, field, :dense_block),
        nothing,
    )
end

function _route_global_combined_gto_moment_matrix_result(
    status::Symbol,
    blocker,
    layout,
    position_x,
    position_y,
    position_z,
    x2_x,
    x2_y,
    x2_z;
    gto_bundle,
    mixed_gausslet_row_range,
    metadata,
)
    materialized =
        status === :materialized_route_global_combined_gto_residual_moment_matrices
    return (;
        object_kind = :route_global_combined_gto_residual_moment_matrices,
        status,
        blocker,
        layout_status =
            _route_global_combined_gto_property(layout, :status, :unavailable),
        layout_kind =
            _route_global_combined_gto_property(layout, :layout_kind, :unavailable),
        gausslet_retained_range =
            _route_global_combined_gto_property(layout, :gausslet_retained_range, :unavailable),
        gausslet_retained_dimension =
            _route_global_combined_gto_property(layout, :gausslet_retained_dimension, :unavailable),
        gto_supplement_range =
            _route_global_combined_gto_property(layout, :gto_supplement_range, :unavailable),
        gto_supplement_orbital_count =
            _route_global_combined_gto_property(layout, :gto_supplement_orbital_count, :unavailable),
        total_combined_dimension =
            _route_global_combined_gto_property(layout, :total_combined_dimension, :unavailable),
        mixed_gausslet_row_range =
            _route_global_combined_gto_effective_mixed_row_range(
                gto_bundle,
                mixed_gausslet_row_range,
            ),
        position_x,
        position_y,
        position_z,
        x2_x,
        x2_y,
        x2_z,
        raw_moment_matrix_fields =
            materialized ? (:position_x, :position_y, :position_z, :x2_x, :x2_y, :x2_z) : (),
        moment_matrix_shapes =
            materialized ?
            (;
                position_x = size(position_x),
                position_y = size(position_y),
                position_z = size(position_z),
                x2_x = size(x2_x),
                x2_y = size(x2_y),
                x2_z = size(x2_z),
            ) :
            (;),
        route_global_combined_gto_residual_moment_matrices_materialized =
            materialized,
        raw_gto_density_density_used_as_final_operator = false,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        generalized_overlap_final_solve = false,
        pqs_transforms_materialized = false,
        exports_or_artifacts = false,
        metadata = NamedTuple(metadata),
    )
end

function _route_global_combined_gto_matrix_blocker(
    layout,
    overlap_result,
    kinetic_result,
    electron_nuclear_by_center_results,
    gto_bundle,
    mixed_gausslet_row_range,
)
    _route_global_combined_gto_property(layout, :status, nothing) ===
        :available_route_global_combined_gto_basis_layout ||
        return :missing_route_global_combined_gto_basis_layout

    isnothing(gto_bundle) && return :missing_provider_level_gto_bundle
    mixed_overlap = _route_global_combined_gto_nested_property(
        gto_bundle,
        (:mixed_blocks, :overlap, :dense_block),
        nothing,
    )
    mixed_overlap isa AbstractMatrix ||
        return :missing_mixed_gto_overlap_block

    gausslet_dimension = _route_global_combined_gto_property(
        layout,
        :gausslet_retained_dimension,
        nothing,
    )
    gausslet_dimension isa Integer && gausslet_dimension > 0 ||
        return :missing_gausslet_retained_dimension
    mixed_range = _route_global_combined_gto_mixed_row_range(
        mixed_overlap,
        mixed_gausslet_row_range,
    )
    if mixed_range != _route_global_combined_gto_property(
        layout,
        :gausslet_retained_range,
        nothing,
    )
        return :missing_mixed_gto_route_global_row_coverage
    end

    overlap_matrix = _route_global_combined_gto_property(
        overlap_result,
        :matrix,
        nothing,
    )
    overlap_matrix isa AbstractMatrix ||
        return :missing_route_global_overlap_matrix
    size(overlap_matrix) == (gausslet_dimension, gausslet_dimension) ||
        return :route_global_overlap_matrix_dimension_mismatch

    kinetic_matrix = _route_global_combined_gto_property(
        kinetic_result,
        :matrix,
        nothing,
    )
    kinetic_matrix isa AbstractMatrix ||
        return :missing_route_global_kinetic_matrix
    size(kinetic_matrix) == (gausslet_dimension, gausslet_dimension) ||
        return :route_global_kinetic_matrix_dimension_mismatch

    mixed_kinetic = _route_global_combined_gto_nested_property(
        gto_bundle,
        (:mixed_blocks, :kinetic, :dense_block),
        nothing,
    )
    mixed_kinetic isa AbstractMatrix ||
        return :missing_mixed_gto_kinetic_block

    gto_count = _route_global_combined_gto_property(
        layout,
        :gto_supplement_orbital_count,
        nothing,
    )
    gto_count isa Integer && gto_count > 0 ||
        return :missing_gto_supplement_orbital_count
    size(mixed_overlap) == (gausslet_dimension, gto_count) ||
        return :mixed_gto_overlap_dimension_mismatch
    size(mixed_kinetic) == (gausslet_dimension, gto_count) ||
        return :mixed_gto_kinetic_dimension_mismatch

    gto_overlap = _route_global_combined_gto_nested_property(
        gto_bundle,
        (:gto_blocks, :overlap, :dense_block),
        nothing,
    )
    gto_overlap isa AbstractMatrix ||
        return :missing_gto_gto_overlap_block
    size(gto_overlap) == (gto_count, gto_count) ||
        return :gto_gto_overlap_dimension_mismatch

    gto_kinetic = _route_global_combined_gto_nested_property(
        gto_bundle,
        (:gto_blocks, :kinetic, :dense_block),
        nothing,
    )
    gto_kinetic isa AbstractMatrix ||
        return :missing_gto_gto_kinetic_block
    size(gto_kinetic) == (gto_count, gto_count) ||
        return :gto_gto_kinetic_dimension_mismatch

    return _route_global_combined_gto_nuclear_blocker(
        layout,
        electron_nuclear_by_center_results,
        gto_bundle,
    )
end

function _route_global_combined_gto_nuclear_blocker(
    layout,
    electron_nuclear_by_center_results,
    gto_bundle,
)
    matrix_results = _route_global_combined_gto_property(
        electron_nuclear_by_center_results,
        :matrix_results,
        nothing,
    )
    matrix_results isa Tuple ||
        return :missing_route_global_electron_nuclear_by_center_matrices
    mixed_blocks = _route_global_combined_gto_property(
        gto_bundle,
        :mixed_nuclear_by_center_blocks,
        (),
    )
    gto_blocks = _route_global_combined_gto_property(
        gto_bundle,
        :gto_nuclear_by_center_blocks,
        (),
    )
    length(matrix_results) == length(mixed_blocks) == length(gto_blocks) ||
        return :combined_gto_nuclear_center_count_mismatch

    gausslet_dimension = layout.gausslet_retained_dimension
    gto_count = layout.gto_supplement_orbital_count
    for index in eachindex(matrix_results)
        gg_matrix = _route_global_combined_gto_property(
            matrix_results[index],
            :matrix,
            nothing,
        )
        gg_matrix isa AbstractMatrix ||
            return :missing_route_global_electron_nuclear_by_center_matrix
        size(gg_matrix) == (gausslet_dimension, gausslet_dimension) ||
            return :route_global_electron_nuclear_matrix_dimension_mismatch

        mixed_matrix = _route_global_combined_gto_property(
            mixed_blocks[index],
            :dense_block,
            nothing,
        )
        mixed_matrix isa AbstractMatrix ||
            return :missing_mixed_gto_electron_nuclear_by_center_block
        size(mixed_matrix) == (gausslet_dimension, gto_count) ||
            return :mixed_gto_electron_nuclear_dimension_mismatch

        gto_matrix = _route_global_combined_gto_property(
            gto_blocks[index],
            :dense_block,
            nothing,
        )
        gto_matrix isa AbstractMatrix ||
            return :missing_gto_gto_electron_nuclear_by_center_block
        size(gto_matrix) == (gto_count, gto_count) ||
            return :gto_gto_electron_nuclear_dimension_mismatch

        _route_global_combined_gto_property(
            matrix_results[index],
            :nuclear_charge_applied,
            false,
        ) === false ||
            return :nuclear_charge_already_applied_to_route_global_by_center_matrix
    end
    return nothing
end

function _route_global_combined_gto_overlap_matrix(layout, gg_matrix, gG_matrix, GG_matrix)
    matrix = zeros(
        Float64,
        layout.total_combined_dimension,
        layout.total_combined_dimension,
    )
    gg = layout.gausslet_retained_range
    gto = layout.gto_supplement_range
    matrix[gg, gg] .= gg_matrix
    matrix[gg, gto] .= gG_matrix
    matrix[gto, gg] .= transpose(gG_matrix)
    matrix[gto, gto] .= GG_matrix
    return matrix
end

function _route_global_combined_gto_nuclear_by_center_matrices(
    layout,
    electron_nuclear_by_center_results,
    gto_bundle,
)
    gg_results = electron_nuclear_by_center_results.matrix_results
    mixed_blocks = gto_bundle.mixed_nuclear_by_center_blocks
    gto_blocks = gto_bundle.gto_nuclear_by_center_blocks
    return Tuple(
        _route_global_combined_gto_nuclear_by_center_matrix(
            layout,
            gg_results[index],
            mixed_blocks[index],
            gto_blocks[index],
        ) for index in eachindex(gg_results)
    )
end

function _route_global_combined_gto_nuclear_by_center_matrix(
    layout,
    gg_result,
    mixed_block,
    gto_block,
)
    matrix = _route_global_combined_gto_overlap_matrix(
        layout,
        gg_result.matrix,
        mixed_block.dense_block,
        gto_block.dense_block,
    )
    return (;
        center_index = _route_global_combined_gto_property(
            gg_result,
            :center_index,
            nothing,
        ),
        center_key = _route_global_combined_gto_property(
            gg_result,
            :center_key,
            nothing,
        ),
        center_location = _route_global_combined_gto_property(
            gg_result,
            :center_location,
            nothing,
        ),
        nuclear_charge = _route_global_combined_gto_property(
            gg_result,
            :nuclear_charge,
            nothing,
        ),
        matrix,
        matrix_shape = size(matrix),
        by_center = true,
        nuclear_charge_applied = false,
        centers_summed = false,
    )
end

function _route_global_combined_gto_matrix_result(
    status::Symbol,
    blocker,
    layout,
    overlap_matrix,
    hamiltonian_matrix,
    nuclear_by_center_matrices;
    gto_bundle,
    mixed_gausslet_row_range,
    metadata,
)
    materialized =
        status === :materialized_route_global_combined_gto_one_electron_matrices
    overlap_diagnostics =
        materialized ?
        _route_global_combined_gto_overlap_diagnostics(overlap_matrix) :
        _route_global_combined_gto_unavailable_overlap_diagnostics()
    return (;
        object_kind = :route_global_combined_gto_one_electron_matrix_assembly,
        status,
        blocker,
        layout_status =
            _route_global_combined_gto_property(layout, :status, :unavailable),
        layout_kind =
            _route_global_combined_gto_property(layout, :layout_kind, :unavailable),
        gausslet_retained_range =
            _route_global_combined_gto_property(layout, :gausslet_retained_range, :unavailable),
        gausslet_retained_dimension =
            _route_global_combined_gto_property(layout, :gausslet_retained_dimension, :unavailable),
        gto_supplement_range =
            _route_global_combined_gto_property(layout, :gto_supplement_range, :unavailable),
        gto_supplement_orbital_count =
            _route_global_combined_gto_property(layout, :gto_supplement_orbital_count, :unavailable),
        total_combined_dimension =
            _route_global_combined_gto_property(layout, :total_combined_dimension, :unavailable),
        block_layout_keys =
            _route_global_combined_gto_property(layout, :block_layout_keys, ()),
        mixed_gausslet_row_range =
            _route_global_combined_gto_effective_mixed_row_range(
                gto_bundle,
                mixed_gausslet_row_range,
            ),
        mixed_gausslet_row_count =
            _route_global_combined_gto_effective_mixed_row_count(gto_bundle),
        mixed_gausslet_row_coverage_status =
            _route_global_combined_gto_mixed_row_coverage_status(
                layout,
                gto_bundle,
                mixed_gausslet_row_range,
            ),
        overlap_matrix,
        hamiltonian_matrix,
        nuclear_by_center_matrices,
        overlap_matrix_shape =
            materialized ? size(overlap_matrix) : :unavailable,
        hamiltonian_matrix_shape =
            materialized ? size(hamiltonian_matrix) : :unavailable,
        by_center_matrix_count = length(nuclear_by_center_matrices),
        nuclear_charge_application_stage =
            :route_global_combined_gto_hamiltonian_assembly,
        nuclear_charge_applied_at_hamiltonian_assembly = materialized,
        center_summation_stage = :route_global_combined_gto_hamiltonian_assembly,
        centers_summed_at_hamiltonian_assembly = materialized,
        by_center_nuclear_matrices_remain_available = materialized,
        overlap_symmetry_error = overlap_diagnostics.symmetry_error,
        overlap_minimum_eigenvalue = overlap_diagnostics.minimum_eigenvalue,
        overlap_maximum_eigenvalue = overlap_diagnostics.maximum_eigenvalue,
        overlap_condition_estimate = overlap_diagnostics.condition_estimate,
        overlap_near_zero_eigenvalue_count =
            overlap_diagnostics.near_zero_eigenvalue_count,
        overlap_negative_eigenvalue_count =
            overlap_diagnostics.negative_eigenvalue_count,
        overlap_positive_definite =
            overlap_diagnostics.positive_definite,
        solve_ready = materialized && overlap_diagnostics.positive_definite,
        combined_overlap_matrix_materialized = materialized,
        combined_hamiltonian_matrix_materialized = materialized,
        route_global_combined_matrix_materialized = materialized,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        pqs_transforms_materialized = false,
        exports_or_artifacts = false,
        metadata = NamedTuple(metadata),
    )
end

function _route_global_combined_gto_overlap_diagnostics(matrix)
    sym = Symmetric((matrix + transpose(matrix)) ./ 2)
    eigenvalues = eigvals(sym)
    minimum_eigenvalue = minimum(eigenvalues)
    maximum_eigenvalue = maximum(eigenvalues)
    return (;
        symmetry_error = maximum(abs.(matrix - transpose(matrix))),
        minimum_eigenvalue,
        maximum_eigenvalue,
        condition_estimate = maximum_eigenvalue / minimum_eigenvalue,
        near_zero_eigenvalue_count = count(abs.(eigenvalues) .<= 1.0e-10),
        negative_eigenvalue_count = count(eigenvalues .< -1.0e-10),
        positive_definite = minimum_eigenvalue > 1.0e-10,
    )
end

function _route_global_combined_gto_unavailable_overlap_diagnostics()
    return (;
        symmetry_error = :unavailable,
        minimum_eigenvalue = :unavailable,
        maximum_eigenvalue = :unavailable,
        condition_estimate = :unavailable,
        near_zero_eigenvalue_count = :unavailable,
        negative_eigenvalue_count = :unavailable,
        positive_definite = false,
    )
end

function _route_global_combined_gto_nested_property(object, fields::Tuple, default)
    current = object
    for field in fields
        hasproperty(current, field) || return default
        current = getproperty(current, field)
    end
    return current
end

function _route_global_combined_gto_mixed_row_range(mixed_matrix, provided_range)
    !isnothing(provided_range) && return provided_range
    return 1:size(mixed_matrix, 1)
end

function _route_global_combined_gto_effective_mixed_row_range(
    gto_bundle,
    provided_range,
)
    mixed_overlap = _route_global_combined_gto_nested_property(
        gto_bundle,
        (:mixed_blocks, :overlap, :dense_block),
        nothing,
    )
    mixed_overlap isa AbstractMatrix || return :unavailable
    return _route_global_combined_gto_mixed_row_range(mixed_overlap, provided_range)
end

function _route_global_combined_gto_effective_mixed_row_count(gto_bundle)
    mixed_overlap = _route_global_combined_gto_nested_property(
        gto_bundle,
        (:mixed_blocks, :overlap, :dense_block),
        nothing,
    )
    mixed_overlap isa AbstractMatrix || return :unavailable
    return size(mixed_overlap, 1)
end

function _route_global_combined_gto_mixed_row_coverage_status(
    layout,
    gto_bundle,
    provided_range,
)
    row_range =
        _route_global_combined_gto_effective_mixed_row_range(gto_bundle, provided_range)
    row_range isa AbstractUnitRange ||
        return :unavailable_mixed_gto_route_global_row_coverage
    return row_range == _route_global_combined_gto_property(
        layout,
        :gausslet_retained_range,
        :unavailable,
    ) ?
           :full_mixed_gto_route_global_row_coverage :
           :partial_mixed_gto_route_global_row_coverage
end
