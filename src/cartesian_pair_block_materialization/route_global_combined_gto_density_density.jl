# Final-basis density-density readiness for a decomposed WL gausslet sector plus
# residualized GTO supplement directions. This does not build raw GTO/GTO
# electron-electron blocks and does not treat raw GTO columns as final MWG data.

function route_global_combined_gto_final_basis_density_density_matrix(
    combined_matrices,
    final_basis_projection;
    gausslet_density_density_result,
    residual_mwg_representation = nothing,
    metadata = (;),
)
    blocker = _route_global_combined_gto_density_density_blocker(
        combined_matrices,
        final_basis_projection,
        gausslet_density_density_result,
    )
    if !isnothing(blocker)
        return _route_global_combined_gto_density_density_result(
            :blocked_route_global_combined_gto_final_basis_density_density_matrix,
            blocker,
            combined_matrices,
            final_basis_projection,
            gausslet_density_density_result,
            nothing;
            residual_mwg_representation,
            metadata,
        )
    end

    if isnothing(residual_mwg_representation)
        return _route_global_combined_gto_density_density_result(
            :blocked_route_global_combined_gto_final_basis_density_density_matrix,
            :missing_residual_gto_to_mwg_effective_gaussian_representation,
            combined_matrices,
            final_basis_projection,
            gausslet_density_density_result,
            nothing;
            residual_mwg_representation,
            metadata,
        )
    end

    return _route_global_combined_gto_density_density_result(
        :blocked_route_global_combined_gto_final_basis_density_density_matrix,
        :missing_residual_mwg_density_density_kernel_for_final_basis_projection,
        combined_matrices,
        final_basis_projection,
        gausslet_density_density_result,
        nothing;
        residual_mwg_representation,
        metadata,
    )
end

function _route_global_combined_gto_density_density_blocker(
    combined_matrices,
    final_basis_projection,
    gausslet_density_density_result,
)
    _route_global_combined_gto_property(combined_matrices, :status, nothing) ===
        :materialized_route_global_combined_gto_one_electron_matrices ||
        return :missing_route_global_combined_gto_matrices
    _route_global_combined_gto_property(final_basis_projection, :status, nothing) ===
        :materialized_route_global_combined_gto_final_basis_projection ||
        return :missing_route_global_combined_gto_final_basis_projection
    _route_global_combined_gto_property(
        gausslet_density_density_result,
        :status,
        nothing,
    ) === :materialized_route_global_density_density_interaction_matrix ||
        return :missing_gausslet_density_density_interaction_matrix

    gausslet_dimension = _route_global_combined_gto_property(
        combined_matrices,
        :gausslet_retained_dimension,
        nothing,
    )
    gausslet_dimension isa Integer || return :missing_gausslet_retained_dimension
    final_dimension = _route_global_combined_gto_property(
        final_basis_projection,
        :final_dimension,
        nothing,
    )
    final_dimension isa Integer || return :missing_final_basis_dimension
    final_dimension >= gausslet_dimension ||
        return :final_basis_dimension_smaller_than_gausslet_sector

    gausslet_matrix =
        _route_global_combined_gto_property(gausslet_density_density_result, :matrix, nothing)
    gausslet_matrix isa AbstractMatrix || return :missing_gausslet_density_density_matrix
    size(gausslet_matrix) == (gausslet_dimension, gausslet_dimension) ||
        return :gausslet_density_density_matrix_shape_mismatch

    transform = _route_global_combined_gto_property(
        final_basis_projection,
        :transformation_matrix,
        nothing,
    )
    transform isa AbstractMatrix || return :missing_final_basis_transformation_matrix
    raw_dimension = _route_global_combined_gto_property(
        combined_matrices,
        :total_combined_dimension,
        nothing,
    )
    raw_dimension isa Integer || return :missing_raw_combined_dimension
    size(transform) == (raw_dimension, final_dimension) ||
        return :final_basis_transformation_matrix_shape_mismatch

    return nothing
end

function _route_global_combined_gto_density_density_result(
    status::Symbol,
    blocker,
    combined_matrices,
    final_basis_projection,
    gausslet_density_density_result,
    final_density_density_matrix;
    residual_mwg_representation,
    metadata,
)
    materialized =
        status ===
        :materialized_route_global_combined_gto_final_basis_density_density_matrix
    gausslet_matrix =
        _route_global_combined_gto_property(gausslet_density_density_result, :matrix, nothing)
    gausslet_dimension = _route_global_combined_gto_property(
        combined_matrices,
        :gausslet_retained_dimension,
        :unavailable,
    )
    raw_supplement_count = _route_global_combined_gto_property(
        final_basis_projection,
        :raw_supplement_count,
        :unavailable,
    )
    retained_supplement_count = _route_global_combined_gto_property(
        final_basis_projection,
        :retained_supplement_count,
        :unavailable,
    )
    final_dimension = _route_global_combined_gto_property(
        final_basis_projection,
        :final_dimension,
        :unavailable,
    )
    gausslet_shape =
        gausslet_matrix isa AbstractMatrix ? size(gausslet_matrix) : :unavailable
    final_shape =
        final_density_density_matrix isa AbstractMatrix ?
        size(final_density_density_matrix) :
        :unavailable
    expected_gausslet_shape =
        gausslet_dimension isa Integer ?
        (gausslet_dimension, gausslet_dimension) :
        :unavailable
    missing_blocks = materialized ? () : (
        :mixed_gausslet_residual_density_density,
        :residual_residual_density_density,
    )
    return (;
        object_kind =
            :route_global_combined_gto_final_basis_density_density_matrix,
        status,
        blocker,
        density_density_block_structure =
            (:gausslet_gausslet, :gausslet_residual, :residual_gausslet, :residual_residual),
        gausslet_density_density_matrix =
            materialized ? gausslet_matrix : nothing,
        gausslet_density_density_matrix_shape = gausslet_shape,
        gausslet_density_density_matrix_matches_existing_wl_block =
            gausslet_shape == expected_gausslet_shape,
        final_density_density_matrix,
        final_density_density_matrix_shape = final_shape,
        gausslet_retained_dimension = gausslet_dimension,
        raw_supplement_count,
        retained_supplement_count,
        final_dimension,
        residual_supplement_kind = :orthonormalized_gto_residual_against_wl,
        residual_mwg_representation_available =
            !isnothing(residual_mwg_representation),
        missing_density_density_blocks = missing_blocks,
        raw_gto_density_density_used_as_final_operator = false,
        generalized_overlap_final_solve = false,
        final_density_density_matrix_materialized = materialized,
        gausslet_gausslet_density_density_available =
            gausslet_matrix isa AbstractMatrix,
        gausslet_gausslet_dense_payload_retained_in_blocked_summary =
            materialized,
        mixed_gausslet_residual_density_density_available = false,
        residual_residual_density_density_available = false,
        gausslet_ida_weight_division_stage =
            :after_final_density_projection_at_density_interaction_boundary,
        gausslet_raw_pgdg_pair_numerator_projected_before_weight_division = true,
        residual_mwg_density_weight_convention =
            :blocked_pending_effective_residual_representation,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        pqs_transforms_materialized = false,
        exports_or_artifacts = false,
        metadata = NamedTuple(metadata),
    )
end
