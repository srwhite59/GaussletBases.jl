# PQS source shell-realization final-basis object.

"""
    pqs_source_shell_realization_final_basis(raw_source_plan, retained_rule; ...)

Materialize the first PQS shell-realization final-basis object from an already
selected raw source and PQS boundary source-mode rule. This builds only the
final shell coefficients and overlap/isometry diagnostics:

```text
R = shell_projection * lowdin_cleanup
S_final = R' * shell_overlap * R
```

It deliberately does not transform one-body operators. A later route must
materialize shell-projected one-body operators before H1 or Hamiltonian data can
be claimed.
"""
function pqs_source_shell_realization_final_basis(
    raw_source_plan::CRPS.RawProductBoxPlan,
    retained_rule::CRPS.PQSBoundaryProductModeRetainedRule;
    shell_support_indices,
    shell_overlap,
    shell_projection,
    lowdin_cleanup,
    identity_atol::Real = 1.0e-8,
    rank_atol::Real = 1.0e-10,
    metadata = (;),
)
    _pqs_shell_final_basis_validate_source(raw_source_plan, retained_rule)
    support_indices = Int[index for index in shell_support_indices]
    support_count = length(support_indices)
    boundary_count = retained_rule.retained_count

    projection = Matrix{Float64}(shell_projection)
    cleanup = Matrix{Float64}(lowdin_cleanup)
    overlap = Matrix{Float64}(shell_overlap)
    size(projection) == (support_count, boundary_count) ||
        throw(DimensionMismatch("PQS shell projection must be shell support rows x boundary source modes"))
    size(cleanup, 1) == boundary_count ||
        throw(DimensionMismatch("PQS Lowdin cleanup row count must match boundary source modes"))
    size(cleanup, 2) > 0 ||
        throw(ArgumentError("PQS Lowdin cleanup must retain at least one final column"))
    size(overlap) == (support_count, support_count) ||
        throw(DimensionMismatch("PQS shell overlap must be shell support rows x shell support rows"))

    final_shell_coefficients = projection * cleanup
    projected_boundary_overlap = transpose(projection) * overlap * projection
    final_overlap =
        transpose(final_shell_coefficients) * overlap * final_shell_coefficients

    projected_diagnostics =
        _pqs_shell_realization_overlap_diagnostics(
            projected_boundary_overlap,
            boundary_count;
            identity_target = false,
            rank_atol,
            identity_atol,
        )
    final_diagnostics =
        _pqs_shell_realization_overlap_diagnostics(
            final_overlap,
            size(cleanup, 2);
            identity_target = true,
            rank_atol,
            identity_atol,
        )
    rank_blocker =
        projected_diagnostics.rank == boundary_count &&
        final_diagnostics.rank == size(cleanup, 2) ? nothing :
        :shell_realization_rank_deficient
    identity_blocker =
        isnothing(rank_blocker) && final_diagnostics.identity_error <= identity_atol ?
        nothing : :final_overlap_not_identity
    blocker = isnothing(rank_blocker) ? identity_blocker : rank_blocker
    status = isnothing(blocker) ?
        :available_pqs_shell_realization_final_basis :
        :blocked_pqs_shell_realization_final_basis

    return (;
        object_kind = :pqs_source_shell_realization_final_basis,
        status,
        blocker,
        source_key = raw_source_plan.source_key,
        source_mode_dims = raw_source_plan.source_mode_dims,
        source_mode_count = raw_source_plan.source_mode_count,
        source_mode_ordering = raw_source_plan.source_mode_ordering,
        retained_rule_kind = retained_rule.retained_rule_kind,
        boundary_source_mode_count = boundary_count,
        boundary_source_column_indices =
            copy(retained_rule.retained_column_indices),
        shell_support_indices = support_indices,
        shell_support_count = support_count,
        final_retained_count = size(cleanup, 2),
        shell_projection = projection,
        lowdin_cleanup = cleanup,
        final_shell_coefficients,
        projected_boundary_overlap,
        final_overlap,
        projected_boundary_overlap_diagnostics = projected_diagnostics,
        final_overlap_diagnostics = final_diagnostics,
        final_overlap_identity_error = final_diagnostics.identity_error,
        final_overlap_is_identity =
            final_diagnostics.identity_error <= identity_atol,
        final_basis_materialized = isnothing(blocker),
        shell_realization_materialized = true,
        lowdin_cleanup_used = true,
        one_body_operator_materialized = false,
        operator_blocks_materialized = false,
        one_body_operator_blocker =
            :missing_pqs_shell_projected_one_body_operator_materialization,
        h1_solve_materialized = false,
        ida_data_materialized = false,
        density_density_materialized = false,
        rhf_materialized = false,
        hamiltonian_data_materialized = false,
        driver_route_materialized = false,
        artifacts_materialized = false,
        metadata = merge(
            NamedTuple(metadata),
            (;
                source = :pqs_source_shell_realization_final_basis,
                raw_source_operator_blocks_transformed = false,
                lowdin_alone_used_as_raw_to_final_transform = false,
                current_route_safe_term_matrices_used = false,
                support_local_operator_contraction_authority = false,
                old_fixed_block_matrix_authority_used = false,
            ),
        ),
    )
end

function _pqs_shell_final_basis_validate_source(
    raw_source_plan::CRPS.RawProductBoxPlan,
    retained_rule::CRPS.PQSBoundaryProductModeRetainedRule,
)
    retained_rule.retained_rule_kind === :boundary_comx_product_mode_selection ||
        throw(ArgumentError("PQS final basis requires the boundary source-mode retained rule"))
    retained_rule.transform_kind === :source_mode_column_selector ||
        throw(ArgumentError("PQS final basis requires source-mode column selector retained rule"))
    retained_rule.source_key == raw_source_plan.source_key ||
        throw(ArgumentError("PQS final basis source key mismatch"))
    retained_rule.source_mode_dims == raw_source_plan.source_mode_dims ||
        throw(ArgumentError("PQS final basis source-mode dimension mismatch"))
    retained_rule.source_mode_ordering == raw_source_plan.source_mode_ordering ||
        throw(ArgumentError("PQS final basis source-mode ordering mismatch"))
    retained_rule.retained_count == length(retained_rule.retained_column_indices) ||
        throw(ArgumentError("PQS final basis retained-column count mismatch"))
    retained_rule.shell_realization_materialized &&
        throw(ArgumentError("PQS final basis expects a raw boundary retained rule without shell realization"))
    retained_rule.lowdin_cleanup_used &&
        throw(ArgumentError("PQS final basis expects a raw boundary retained rule without Lowdin cleanup"))
    return nothing
end

function _pqs_shell_realization_overlap_diagnostics(
    overlap::AbstractMatrix{<:Real},
    expected_dimension::Int;
    identity_target::Bool,
    rank_atol::Real,
    identity_atol::Real,
)
    size(overlap) == (expected_dimension, expected_dimension) ||
        throw(DimensionMismatch("PQS shell-realization overlap shape mismatch"))
    all(isfinite, overlap) ||
        throw(ArgumentError("PQS shell-realization overlap contains non-finite entries"))
    overlap_matrix = Matrix{Float64}(overlap)
    symmetric_overlap = Symmetric((overlap_matrix + transpose(overlap_matrix)) ./ 2)
    symmetry_error = norm(overlap_matrix - transpose(overlap_matrix), Inf)
    eigenvalues = eigvals(symmetric_overlap)
    min_eigenvalue = minimum(eigenvalues)
    max_eigenvalue = maximum(eigenvalues)
    threshold = max(Float64(rank_atol), eps(Float64) * max(expected_dimension, 1))
    rank = count(value -> value > threshold, eigenvalues)
    identity_error = identity_target ?
        norm(
            overlap_matrix -
            Matrix{Float64}(I, expected_dimension, expected_dimension),
            Inf,
        ) :
        nothing
    return (;
        shape = size(overlap_matrix),
        symmetric = symmetry_error <= Float64(identity_atol),
        symmetry_error,
        eigenvalue_min = min_eigenvalue,
        eigenvalue_max = max_eigenvalue,
        rank,
        expected_dimension,
        rank_atol = Float64(rank_atol),
        identity_target,
        identity_error,
        full_rank = rank == expected_dimension,
    )
end
