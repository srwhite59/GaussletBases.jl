# Complete one-center direct-core plus surrounding-shell final-basis realization.

"""
    pqs_complete_core_shell_final_basis(; ...)

Build the narrow route-owned final-basis realization for a direct core plus an
independent surrounding PQS shell sector.

The caller supplies the core support overlap, shell support overlap,
core/shell cross overlap, and already Lowdin-cleaned shell coefficients. This
helper only builds the combined final basis and overlap diagnostics; it does
not materialize one-body operators, Hamiltonians, IDA data, RHF data, driver
wiring, exports, or artifacts.
"""
function pqs_complete_core_shell_final_basis(;
    core_support_indices,
    shell_support_indices,
    core_overlap,
    core_shell_overlap,
    shell_overlap,
    shell_final_coefficients,
    identity_atol::Real = 1.0e-8,
    rank_atol::Real = 1.0e-10,
    metadata = (;),
)
    core_indices = Int[index for index in core_support_indices]
    shell_indices = Int[index for index in shell_support_indices]
    core_count = length(core_indices)
    shell_support_count = length(shell_indices)
    isempty(core_indices) &&
        throw(ArgumentError("complete core/shell final basis requires core support"))
    isempty(shell_indices) &&
        throw(ArgumentError("complete core/shell final basis requires shell support"))
    disjoint_support = isempty(intersect(core_indices, shell_indices))
    disjoint_support ||
        throw(ArgumentError("complete core/shell final basis requires disjoint core and shell support"))

    core_matrix = Matrix{Float64}(core_overlap)
    cross_matrix = Matrix{Float64}(core_shell_overlap)
    shell_matrix = Matrix{Float64}(shell_overlap)
    shell_coefficients = Matrix{Float64}(shell_final_coefficients)
    size(core_matrix) == (core_count, core_count) ||
        throw(DimensionMismatch("core overlap shape must match core support"))
    size(cross_matrix) == (core_count, shell_support_count) ||
        throw(DimensionMismatch("core/shell overlap shape must be core support x shell support"))
    size(shell_matrix) == (shell_support_count, shell_support_count) ||
        throw(DimensionMismatch("shell overlap shape must match shell support"))
    size(shell_coefficients, 1) == shell_support_count ||
        throw(DimensionMismatch("shell final coefficients row count must match shell support"))
    size(shell_coefficients, 2) > 0 ||
        throw(ArgumentError("shell final coefficients must retain at least one shell column"))
    all(isfinite, core_matrix) ||
        throw(ArgumentError("core overlap contains non-finite entries"))
    all(isfinite, cross_matrix) ||
        throw(ArgumentError("core/shell overlap contains non-finite entries"))
    all(isfinite, shell_matrix) ||
        throw(ArgumentError("shell overlap contains non-finite entries"))
    all(isfinite, shell_coefficients) ||
        throw(ArgumentError("shell final coefficients contain non-finite entries"))

    shell_count = size(shell_coefficients, 2)
    total_count = core_count + shell_count
    core_range = 1:core_count
    shell_range = (core_count + 1):total_count

    shell_final_overlap =
        transpose(shell_coefficients) * shell_matrix * shell_coefficients
    pre_final_overlap = Matrix{Float64}(undef, total_count, total_count)
    pre_final_overlap[core_range, core_range] .= core_matrix
    pre_final_overlap[core_range, shell_range] .= cross_matrix * shell_coefficients
    pre_final_overlap[shell_range, core_range] .=
        transpose(shell_coefficients) * transpose(cross_matrix)
    pre_final_overlap[shell_range, shell_range] .= shell_final_overlap

    diagnostics = _pqs_complete_core_shell_overlap_diagnostics(
        pre_final_overlap;
        identity_atol,
        rank_atol,
    )
    if !diagnostics.full_rank
        return _pqs_complete_core_shell_final_basis_result(
            :blocked_pqs_complete_core_shell_final_basis,
            :combined_core_shell_overlap_rank_deficient,
            core_indices,
            shell_indices,
            core_count,
            shell_support_count,
            shell_count,
            core_range,
            shell_range,
            disjoint_support,
            shell_final_overlap,
            pre_final_overlap,
            diagnostics.summary,
            nothing,
            nothing,
            nothing,
            metadata,
        )
    end

    cleanup =
        diagnostics.eigenvectors * Diagonal(1.0 ./ sqrt.(diagnostics.eigenvalues))
    final_overlap = transpose(cleanup) * pre_final_overlap * cleanup
    final_identity_error = norm(
        final_overlap -
        Matrix{Float64}(I, size(final_overlap, 1), size(final_overlap, 2)),
        Inf,
    )
    blocker = final_identity_error <= Float64(identity_atol) ?
        nothing : :combined_core_shell_final_overlap_not_identity
    status = isnothing(blocker) ?
        :available_pqs_complete_core_shell_final_basis :
        :blocked_pqs_complete_core_shell_final_basis

    pre_final_coefficients = zeros(
        Float64,
        core_count + shell_support_count,
        total_count,
    )
    pre_final_coefficients[1:core_count, core_range] .=
        Matrix{Float64}(I, core_count, core_count)
    shell_support_rows = (core_count + 1):(core_count + shell_support_count)
    pre_final_coefficients[shell_support_rows, shell_range] .= shell_coefficients
    final_coefficients = pre_final_coefficients * cleanup

    return _pqs_complete_core_shell_final_basis_result(
        status,
        blocker,
        core_indices,
        shell_indices,
        core_count,
        shell_support_count,
        shell_count,
        core_range,
        shell_range,
        disjoint_support,
        shell_final_overlap,
        pre_final_overlap,
        diagnostics.summary,
        cleanup,
        final_overlap,
        (pre_final_coefficients = pre_final_coefficients, final_coefficients = final_coefficients),
        metadata,
        final_identity_error = final_identity_error,
    )
end

function _pqs_complete_core_shell_final_basis_result(
    status::Symbol,
    blocker,
    core_indices,
    shell_indices,
    core_count::Int,
    shell_support_count::Int,
    shell_count::Int,
    core_range,
    shell_range,
    disjoint_support::Bool,
    shell_final_overlap,
    pre_final_overlap,
    diagnostics,
    cleanup,
    final_overlap,
    coefficients,
    metadata;
    final_identity_error = nothing,
)
    materialized = status === :available_pqs_complete_core_shell_final_basis
    return (;
        object_kind = :pqs_complete_core_shell_final_basis,
        status,
        blocker,
        core_support_indices = core_indices,
        shell_support_indices = shell_indices,
        support_row_order = :core_then_shell,
        core_support_count = core_count,
        shell_support_count,
        shell_final_retained_count = shell_count,
        final_retained_count = core_count + shell_count,
        core_range,
        shell_range,
        core_shell_support_disjoint = disjoint_support,
        shell_final_overlap,
        shell_final_overlap_identity_error = norm(
            shell_final_overlap -
            Matrix{Float64}(I, size(shell_final_overlap, 1), size(shell_final_overlap, 2)),
            Inf,
        ),
        pre_final_overlap,
        pre_final_overlap_diagnostics = diagnostics,
        combined_lowdin_cleanup = cleanup,
        final_overlap,
        final_overlap_identity_error = final_identity_error,
        final_overlap_is_identity =
            !isnothing(final_identity_error) && final_identity_error <= diagnostics.identity_atol,
        pre_final_coefficients = isnothing(coefficients) ? nothing : coefficients.pre_final_coefficients,
        final_coefficients = isnothing(coefficients) ? nothing : coefficients.final_coefficients,
        final_basis_materialized = materialized,
        shell_realization_materialized = true,
        combined_lowdin_cleanup_used = materialized,
        one_body_operator_materialized = false,
        operator_blocks_materialized = false,
        h1_solve_materialized = false,
        generalized_overlap_solve_materialized = false,
        hamiltonian_data_materialized = false,
        ida_data_materialized = false,
        density_density_materialized = false,
        rhf_materialized = false,
        driver_route_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
        metadata = merge(
            NamedTuple(metadata),
            (;
                source = :pqs_complete_core_shell_final_basis,
                core_sector = :direct_core_support_identity,
                shell_sector = :surrounding_shell_lowdin_coefficients,
                old_fixed_block_matrix_authority_used = false,
                current_route_safe_term_matrices_used = false,
            ),
        ),
    )
end

function _pqs_complete_core_shell_overlap_diagnostics(
    overlap::AbstractMatrix{<:Real};
    identity_atol::Real,
    rank_atol::Real,
)
    matrix = Matrix{Float64}(overlap)
    size(matrix, 1) == size(matrix, 2) ||
        throw(DimensionMismatch("complete core/shell overlap must be square"))
    all(isfinite, matrix) ||
        throw(ArgumentError("complete core/shell overlap contains non-finite entries"))
    symmetric_matrix = Symmetric((matrix + transpose(matrix)) ./ 2)
    symmetry_error = norm(matrix - transpose(matrix), Inf)
    decomposition = eigen(symmetric_matrix)
    eigenvalues = decomposition.values
    threshold = max(Float64(rank_atol), eps(Float64) * max(size(matrix, 1), 1))
    rank = count(value -> value > threshold, eigenvalues)
    summary = (;
        shape = size(matrix),
        symmetry_error,
        eigenvalue_min = minimum(eigenvalues),
        eigenvalue_max = maximum(eigenvalues),
        rank,
        expected_dimension = size(matrix, 1),
        rank_atol = Float64(rank_atol),
        identity_atol = Float64(identity_atol),
        full_rank = rank == size(matrix, 1),
    )
    return (;
        summary,
        eigenvalues,
        eigenvectors = decomposition.vectors,
        full_rank = rank == size(matrix, 1),
    )
end
