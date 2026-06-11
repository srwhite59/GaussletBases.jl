# Final-basis density-density readiness for a decomposed WL gausslet sector plus
# residualized GTO supplement directions. This does not build raw GTO/GTO
# electron-electron blocks and does not treat raw GTO columns as final MWG data.

function route_global_residual_gto_mwg_representation(
    final_basis_projection,
    combined_matrices,
    supplement;
    raw_moment_matrices = nothing,
    residual_accept_tolerance = 1.0e-8,
    metadata = (;),
)
    blocker = _route_global_residual_gto_mwg_input_blocker(
        final_basis_projection,
        combined_matrices,
        supplement,
        residual_accept_tolerance,
    )
    if !isnothing(blocker)
        return _route_global_residual_gto_mwg_result(
            :blocked_route_global_residual_gto_mwg_representation,
            blocker,
            final_basis_projection,
            combined_matrices,
            supplement,
            nothing,
            nothing,
            :unavailable;
            raw_moment_matrices,
            residual_accept_tolerance,
            metadata,
        )
    end
    if isnothing(raw_moment_matrices)
        return _route_global_residual_gto_mwg_result(
            :blocked_route_global_residual_gto_mwg_representation,
            :missing_route_global_combined_gto_residual_moment_matrices,
            final_basis_projection,
            combined_matrices,
            supplement,
            nothing,
            nothing,
            :unavailable;
            raw_moment_matrices,
            residual_accept_tolerance,
            metadata,
        )
    end

    moment_blocker = _route_global_residual_gto_mwg_moment_blocker(
        combined_matrices,
        raw_moment_matrices,
    )
    if !isnothing(moment_blocker)
        return _route_global_residual_gto_mwg_result(
            :blocked_route_global_residual_gto_mwg_representation,
            moment_blocker,
            final_basis_projection,
            combined_matrices,
            supplement,
            nothing,
            nothing,
            :unavailable;
            raw_moment_matrices,
            residual_accept_tolerance,
            metadata,
        )
    end

    raw_overlap = _route_global_combined_gto_property(
        combined_matrices,
        :overlap_matrix,
        nothing,
    )
    transform = _route_global_combined_gto_property(
        final_basis_projection,
        :transformation_matrix,
        nothing,
    )
    gausslet_dimension = Int(
        _route_global_combined_gto_property(
            combined_matrices,
            :gausslet_retained_dimension,
            0,
        ),
    )
    moment_data = try
        center_data = ParentGaussletBases._qwrg_residual_center_data(
            raw_overlap,
            _route_global_residual_gto_moment_matrix(raw_moment_matrices, :position_x),
            _route_global_residual_gto_moment_matrix(raw_moment_matrices, :position_y),
            _route_global_residual_gto_moment_matrix(raw_moment_matrices, :position_z),
            transform,
            gausslet_dimension,
        )
        ParentGaussletBases._qwrg_residual_moment_data(
            raw_overlap,
            _route_global_residual_gto_moment_matrix(raw_moment_matrices, :position_x),
            _route_global_residual_gto_moment_matrix(raw_moment_matrices, :x2_x),
            _route_global_residual_gto_moment_matrix(raw_moment_matrices, :x2_y),
            _route_global_residual_gto_moment_matrix(raw_moment_matrices, :x2_z),
            center_data,
        )
    catch error
        return _route_global_residual_gto_mwg_result(
            :blocked_route_global_residual_gto_mwg_representation,
            :invalid_residual_gto_mwg_moment_data,
            final_basis_projection,
            combined_matrices,
            supplement,
            nothing,
            nothing,
            :unavailable;
            raw_moment_matrices,
            residual_accept_tolerance,
            metadata = merge(
                NamedTuple(metadata),
                (;
                    residual_mwg_moment_error_type = typeof(error),
                    residual_mwg_moment_error = sprint(showerror, error),
                ),
            ),
        )
    end
    moment_data.overlap_error <= residual_accept_tolerance ||
        return _route_global_residual_gto_mwg_result(
            :blocked_route_global_residual_gto_mwg_representation,
            :residual_gto_mwg_final_basis_not_orthonormal,
            final_basis_projection,
            combined_matrices,
            supplement,
            nothing,
            nothing,
            moment_data.overlap_error;
            raw_moment_matrices,
            residual_accept_tolerance,
            metadata,
        )

    return _route_global_residual_gto_mwg_result(
        :materialized_route_global_residual_gto_mwg_representation,
        nothing,
        final_basis_projection,
        combined_matrices,
        supplement,
        Matrix{Float64}(moment_data.centers),
        Matrix{Float64}(moment_data.widths),
        moment_data.overlap_error;
        raw_moment_matrices,
        residual_accept_tolerance,
        metadata,
    )
end

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
    _route_global_combined_gto_property(
        residual_mwg_representation,
        :status,
        nothing,
    ) === :materialized_route_global_residual_gto_mwg_representation ||
        return _route_global_combined_gto_density_density_result(
            :blocked_route_global_combined_gto_final_basis_density_density_matrix,
            _route_global_combined_gto_property(
                residual_mwg_representation,
                :blocker,
                :blocked_residual_gto_mwg_representation,
            ),
            combined_matrices,
            final_basis_projection,
            gausslet_density_density_result,
            nothing;
            residual_mwg_representation,
            metadata,
        )

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

function _route_global_residual_gto_mwg_input_blocker(
    final_basis_projection,
    combined_matrices,
    supplement,
    residual_accept_tolerance,
)
    residual_accept_tolerance isa Real && residual_accept_tolerance > 0 ||
        return :invalid_residual_gto_mwg_accept_tolerance
    _route_global_combined_gto_property(final_basis_projection, :status, nothing) ===
        :materialized_route_global_combined_gto_final_basis_projection ||
        return :missing_route_global_combined_gto_final_basis_projection
    _route_global_combined_gto_property(combined_matrices, :status, nothing) ===
        :materialized_route_global_combined_gto_one_electron_matrices ||
        return :missing_route_global_combined_gto_matrices
    _route_global_combined_gto_has_orbitals(supplement) ||
        return :missing_gto_supplement_orbitals
    raw_supplement_count = _route_global_combined_gto_property(
        final_basis_projection,
        :raw_supplement_count,
        nothing,
    )
    raw_supplement_count isa Integer || return :missing_raw_supplement_count
    length(supplement.orbitals) == raw_supplement_count ||
        return :residual_gto_raw_supplement_count_mismatch
    retained_supplement_count = _route_global_combined_gto_property(
        final_basis_projection,
        :retained_supplement_count,
        nothing,
    )
    retained_supplement_count isa Integer && retained_supplement_count > 0 ||
        return :missing_retained_residual_supplement_count
    final_dimension = _route_global_combined_gto_property(
        final_basis_projection,
        :final_dimension,
        nothing,
    )
    final_dimension isa Integer || return :missing_final_basis_dimension
    gausslet_dimension = _route_global_combined_gto_property(
        combined_matrices,
        :gausslet_retained_dimension,
        nothing,
    )
    gausslet_dimension isa Integer || return :missing_gausslet_retained_dimension
    final_dimension == gausslet_dimension + retained_supplement_count ||
        return :residual_gto_final_dimension_mismatch
    raw_dimension = _route_global_combined_gto_property(
        combined_matrices,
        :total_combined_dimension,
        nothing,
    )
    raw_dimension isa Integer || return :missing_raw_combined_dimension
    raw_dimension == gausslet_dimension + raw_supplement_count ||
        return :residual_gto_raw_dimension_mismatch
    transform = _route_global_combined_gto_property(
        final_basis_projection,
        :transformation_matrix,
        nothing,
    )
    transform isa AbstractMatrix || return :missing_final_basis_transformation_matrix
    size(transform) == (raw_dimension, final_dimension) ||
        return :final_basis_transformation_matrix_shape_mismatch
    overlap =
        _route_global_combined_gto_property(combined_matrices, :overlap_matrix, nothing)
    overlap isa AbstractMatrix || return :missing_combined_overlap_matrix
    size(overlap) == (raw_dimension, raw_dimension) ||
        return :combined_overlap_matrix_shape_mismatch
    return nothing
end

function _route_global_residual_gto_mwg_moment_blocker(
    combined_matrices,
    raw_moment_matrices,
)
    raw_dimension = _route_global_combined_gto_property(
        combined_matrices,
        :total_combined_dimension,
        nothing,
    )
    raw_dimension isa Integer || return :missing_raw_combined_dimension
    for field in (:position_x, :position_y, :position_z, :x2_x, :x2_y, :x2_z)
        matrix =
            _route_global_residual_gto_moment_matrix(raw_moment_matrices, field)
        matrix isa AbstractMatrix ||
            return Symbol("missing_route_global_combined_gto_", String(field), "_matrix")
        size(matrix) == (raw_dimension, raw_dimension) ||
            return Symbol("combined_gto_", String(field), "_matrix_shape_mismatch")
    end
    return nothing
end

function _route_global_residual_gto_moment_matrix(raw_moment_matrices, field::Symbol)
    return _route_global_combined_gto_property(raw_moment_matrices, field, nothing)
end

function _route_global_residual_gto_mwg_result(
    status::Symbol,
    blocker,
    final_basis_projection,
    combined_matrices,
    supplement,
    residual_centers,
    residual_widths,
    residual_overlap_error;
    raw_moment_matrices,
    residual_accept_tolerance,
    metadata,
)
    materialized =
        status === :materialized_route_global_residual_gto_mwg_representation
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
    gausslet_dimension = _route_global_combined_gto_property(
        combined_matrices,
        :gausslet_retained_dimension,
        :unavailable,
    )
    final_dimension = _route_global_combined_gto_property(
        final_basis_projection,
        :final_dimension,
        :unavailable,
    )
    moment_matrix_fields =
        isnothing(raw_moment_matrices) ? () :
        Tuple(
            field for field in (:position_x, :position_y, :position_z, :x2_x, :x2_y, :x2_z)
            if _route_global_residual_gto_moment_matrix(raw_moment_matrices, field) isa
               AbstractMatrix
        )
    return (;
        object_kind = :route_global_residual_gto_mwg_representation,
        status,
        blocker,
        representation_kind = :matched_width_gaussian_from_final_residual_moments,
        residual_supplement_kind = :orthonormalized_gto_residual_against_wl,
        raw_supplement_count,
        retained_supplement_count,
        gausslet_retained_dimension = gausslet_dimension,
        final_dimension,
        supplement_orbital_count = _route_global_combined_gto_orbital_count(supplement),
        residual_centers = materialized ? residual_centers : nothing,
        residual_widths = materialized ? residual_widths : nothing,
        residual_centers_shape =
            materialized ? size(residual_centers) : :unavailable,
        residual_widths_shape =
            materialized ? size(residual_widths) : :unavailable,
        residual_widths_finite =
            materialized ? all(isfinite, residual_widths) : false,
        residual_widths_positive =
            materialized ? all(>(0.0), residual_widths) : false,
        residual_overlap_error,
        residual_accept_tolerance,
        residual_moment_source =
            :raw_combined_overlap_position_and_x2_matrices,
        raw_moment_matrix_fields = moment_matrix_fields,
        raw_gto_density_density_used_as_final_operator = false,
        density_weight_division_stage =
            :after_final_density_projection_at_density_interaction_boundary,
        pqs_weight_division_future_stage = :after_projection_or_lowdin,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        generalized_overlap_final_solve = false,
        pqs_transforms_materialized = false,
        exports_or_artifacts = false,
        metadata = NamedTuple(metadata),
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
            _route_global_combined_gto_property(
                residual_mwg_representation,
                :status,
                nothing,
            ) === :materialized_route_global_residual_gto_mwg_representation,
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
