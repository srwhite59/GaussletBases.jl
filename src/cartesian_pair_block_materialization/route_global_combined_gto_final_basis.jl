# Final-basis projection for a decomposed WL gausslet sector plus a small GTO
# supplement. This consumes already assembled raw combined overlap/Hamiltonian
# matrices and constructs an orthonormal residual supplement basis against the
# final gausslet sector. It does not add generalized-overlap transfer logic,
# direct Cartesian fallback, PQS transforms, exports, or artifacts.

function route_global_combined_gto_final_basis_projection(
    combined_matrices;
    residual_drop_tolerance = 1.0e-10,
    identity_tolerance = 1.0e-9,
    metadata = (;),
)
    blocker = _route_global_combined_gto_final_basis_blocker(
        combined_matrices,
        residual_drop_tolerance,
        identity_tolerance,
    )
    if !isnothing(blocker)
        return _route_global_combined_gto_final_basis_result(
            :blocked_route_global_combined_gto_final_basis_projection,
            blocker,
            combined_matrices,
            nothing,
            nothing,
            nothing,
            ();
            residual_drop_tolerance,
            identity_tolerance,
            metadata,
        )
    end

    overlap = Matrix{Float64}(combined_matrices.overlap_matrix)
    hamiltonian = Matrix{Float64}(combined_matrices.hamiltonian_matrix)
    gausslet_range = combined_matrices.gausslet_retained_range
    supplement_range = combined_matrices.gto_supplement_range
    gausslet_dimension = length(gausslet_range)
    raw_supplement_count = length(supplement_range)

    s_ga = overlap[gausslet_range, supplement_range]
    s_ag = transpose(s_ga)
    s_aa = overlap[supplement_range, supplement_range]
    residual_overlap = Symmetric((s_aa - s_ag * s_ga + transpose(s_aa - s_ag * s_ga)) ./ 2)
    residual_eigen = eigen(residual_overlap)
    keep_indices = findall(>(residual_drop_tolerance), residual_eigen.values)
    retained_supplement_count = length(keep_indices)
    retained_supplement_count > 0 ||
        return _route_global_combined_gto_final_basis_result(
            :blocked_route_global_combined_gto_final_basis_projection,
            :empty_gto_supplement_residual_basis,
            combined_matrices,
            nothing,
            residual_overlap,
            nothing,
            residual_eigen.values;
            residual_drop_tolerance,
            identity_tolerance,
            retained_supplement_count,
            metadata,
        )

    kept_vectors = residual_eigen.vectors[:, keep_indices]
    kept_values = residual_eigen.values[keep_indices]
    residual_normalizer = kept_vectors * Diagonal(1.0 ./ sqrt.(kept_values))
    transform = zeros(
        Float64,
        gausslet_dimension + raw_supplement_count,
        gausslet_dimension + retained_supplement_count,
    )
    transform[gausslet_range, 1:gausslet_dimension] .=
        Matrix{Float64}(I, gausslet_dimension, gausslet_dimension)
    transform[gausslet_range, (gausslet_dimension + 1):end] .=
        -s_ga * residual_normalizer
    transform[supplement_range, (gausslet_dimension + 1):end] .=
        residual_normalizer

    final_overlap = Matrix{Float64}(transpose(transform) * overlap * transform)
    final_hamiltonian =
        Matrix{Float64}(transpose(transform) * hamiltonian * transform)
    final_identity_error =
        _route_global_combined_gto_final_identity_error(final_overlap)
    if final_identity_error > identity_tolerance
        return _route_global_combined_gto_final_basis_result(
            :blocked_route_global_combined_gto_final_basis_projection,
            :final_basis_overlap_not_identity,
            combined_matrices,
            transform,
            residual_overlap,
            final_overlap,
            residual_eigen.values;
            residual_drop_tolerance,
            identity_tolerance,
            retained_supplement_count,
            final_hamiltonian,
            metadata,
        )
    end

    return _route_global_combined_gto_final_basis_result(
        :materialized_route_global_combined_gto_final_basis_projection,
        nothing,
        combined_matrices,
        transform,
        residual_overlap,
        final_overlap,
        residual_eigen.values;
        residual_drop_tolerance,
        identity_tolerance,
        retained_supplement_count,
        final_hamiltonian,
        metadata,
    )
end

function _route_global_combined_gto_final_basis_blocker(
    combined_matrices,
    residual_drop_tolerance,
    identity_tolerance,
)
    residual_drop_tolerance isa Real && residual_drop_tolerance > 0 ||
        return :invalid_gto_residual_drop_tolerance
    identity_tolerance isa Real && identity_tolerance > 0 ||
        return :invalid_final_basis_identity_tolerance
    _route_global_combined_gto_property(combined_matrices, :status, nothing) ===
        :materialized_route_global_combined_gto_one_electron_matrices ||
        return :missing_route_global_combined_gto_matrices
    overlap =
        _route_global_combined_gto_property(combined_matrices, :overlap_matrix, nothing)
    hamiltonian =
        _route_global_combined_gto_property(
            combined_matrices,
            :hamiltonian_matrix,
            nothing,
        )
    overlap isa AbstractMatrix || return :missing_combined_overlap_matrix
    hamiltonian isa AbstractMatrix || return :missing_combined_hamiltonian_matrix
    size(overlap) == size(hamiltonian) || return :combined_matrix_shape_mismatch
    gausslet_range =
        _route_global_combined_gto_property(
            combined_matrices,
            :gausslet_retained_range,
            nothing,
        )
    supplement_range =
        _route_global_combined_gto_property(
            combined_matrices,
            :gto_supplement_range,
            nothing,
        )
    gausslet_range isa AbstractUnitRange{<:Integer} ||
        return :missing_gausslet_retained_range
    supplement_range isa AbstractUnitRange{<:Integer} ||
        return :missing_gto_supplement_range
    last(supplement_range) <= size(overlap, 1) ||
        return :combined_basis_range_dimension_mismatch
    s_gg = overlap[gausslet_range, gausslet_range]
    identity_error = maximum(
        abs.(s_gg - Matrix{Float64}(I, size(s_gg, 1), size(s_gg, 2))),
    )
    identity_error <= identity_tolerance ||
        return :gausslet_final_overlap_not_identity
    return nothing
end

function _route_global_combined_gto_final_basis_result(
    status::Symbol,
    blocker,
    combined_matrices,
    transform,
    residual_overlap,
    final_overlap,
    residual_eigenvalues;
    residual_drop_tolerance,
    identity_tolerance,
    retained_supplement_count = 0,
    final_hamiltonian = nothing,
    metadata = (;),
)
    materialized =
        status === :materialized_route_global_combined_gto_final_basis_projection
    raw_supplement_count =
        _route_global_combined_gto_property(
            combined_matrices,
            :gto_supplement_orbital_count,
            :unavailable,
        )
    gausslet_dimension =
        _route_global_combined_gto_property(
            combined_matrices,
            :gausslet_retained_dimension,
            :unavailable,
        )
    final_dimension =
        gausslet_dimension isa Integer ? gausslet_dimension + retained_supplement_count :
        :unavailable
    final_identity_error =
        materialized ? _route_global_combined_gto_final_identity_error(final_overlap) :
        :unavailable
    final_overlap_diagnostics =
        materialized ?
        _route_global_combined_gto_overlap_diagnostics(final_overlap) :
        _route_global_combined_gto_unavailable_overlap_diagnostics()
    return (;
        object_kind = :route_global_combined_gto_final_basis_projection,
        status,
        blocker,
        projection_kind = :gto_residual_orthogonalization_against_decomposed_wl,
        raw_supplement_count,
        retained_supplement_count,
        dropped_supplement_count =
            raw_supplement_count isa Integer ?
            raw_supplement_count - retained_supplement_count :
            :unavailable,
        residual_overlap_eigenvalues = Tuple(residual_eigenvalues),
        residual_drop_tolerance,
        final_dimension,
        final_overlap_identity_error = final_identity_error,
        final_overlap_identity_tolerance = identity_tolerance,
        final_overlap_identity =
            materialized && final_identity_error <= identity_tolerance,
        transformation_provenance =
            :residual_supplement_eigendecomposition_from_combined_overlap,
        transformation_matrix = transform,
        residual_overlap,
        final_overlap_matrix = final_overlap,
        final_hamiltonian_matrix = final_hamiltonian,
        final_overlap_matrix_shape =
            materialized ? size(final_overlap) : :unavailable,
        final_hamiltonian_matrix_shape =
            materialized ? size(final_hamiltonian) : :unavailable,
        final_overlap_symmetry_error = final_overlap_diagnostics.symmetry_error,
        final_overlap_minimum_eigenvalue =
            final_overlap_diagnostics.minimum_eigenvalue,
        final_overlap_maximum_eigenvalue =
            final_overlap_diagnostics.maximum_eigenvalue,
        final_overlap_condition_estimate =
            final_overlap_diagnostics.condition_estimate,
        final_overlap_near_zero_eigenvalue_count =
            final_overlap_diagnostics.near_zero_eigenvalue_count,
        final_overlap_negative_eigenvalue_count =
            final_overlap_diagnostics.negative_eigenvalue_count,
        final_basis_overlap_identity =
            materialized && final_identity_error <= identity_tolerance,
        final_basis_hamiltonian_available = materialized,
        final_orthogonalized_basis_available = materialized,
        generalized_overlap_final_solve = false,
        ordinary_hermitian_final_solve_ready =
            materialized && final_identity_error <= identity_tolerance,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        pqs_transforms_materialized = false,
        exports_or_artifacts = false,
        metadata = NamedTuple(metadata),
    )
end

function _route_global_combined_gto_final_identity_error(final_overlap)
    identity = Matrix{Float64}(I, size(final_overlap, 1), size(final_overlap, 2))
    return maximum(abs.(final_overlap - identity))
end
