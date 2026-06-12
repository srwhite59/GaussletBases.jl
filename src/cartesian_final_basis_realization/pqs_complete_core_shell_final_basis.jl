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

"""
    pqs_complete_core_shell_final_one_body_matrix(final_basis, support_operator; term)

Transform a caller-supplied one-body operator on the combined core/shell support
rows into the complete final basis. The support row order must match
`final_basis.support_row_order == :core_then_shell`.
"""
function pqs_complete_core_shell_final_one_body_matrix(
    final_basis::NamedTuple,
    support_operator;
    term::Symbol,
    center_record = nothing,
    symmetry_atol::Real = 1.0e-8,
    metadata = (;),
)
    _pqs_complete_core_shell_validate_final_basis(final_basis)
    operator = Matrix{Float64}(support_operator)
    support_count = final_basis.core_support_count + final_basis.shell_support_count
    size(operator) == (support_count, support_count) ||
        throw(DimensionMismatch("complete core/shell support operator shape mismatch"))
    all(isfinite, operator) ||
        throw(ArgumentError("complete core/shell support operator contains non-finite entries"))
    support_symmetry_error = norm(operator - transpose(operator), Inf)
    support_symmetry_error <= Float64(symmetry_atol) ||
        throw(ArgumentError("complete core/shell support operator must be symmetric"))

    coefficients = Matrix{Float64}(final_basis.final_coefficients)
    size(coefficients, 1) == support_count ||
        throw(DimensionMismatch("complete core/shell final coefficient row mismatch"))
    final_operator = transpose(coefficients) * operator * coefficients
    final_symmetry_error = norm(final_operator - transpose(final_operator), Inf)
    term_metadata = _pqs_complete_core_shell_one_body_metadata(
        term,
        center_record,
        metadata,
    )
    return (;
        object_kind = :pqs_complete_core_shell_final_one_body_matrix,
        status = :materialized_pqs_complete_core_shell_final_one_body_matrix,
        blocker = nothing,
        term,
        final_basis_object_kind = final_basis.object_kind,
        final_basis_status = final_basis.status,
        support_row_order = final_basis.support_row_order,
        support_operator = operator,
        support_operator_shape = size(operator),
        support_operator_symmetry_error = support_symmetry_error,
        final_operator,
        final_operator_shape = size(final_operator),
        final_operator_finite = all(isfinite, final_operator),
        final_operator_symmetry_error = final_symmetry_error,
        final_retained_count = final_basis.final_retained_count,
        route_owned_product_operator_used = true,
        old_fixed_block_matrix_authority_used = false,
        current_route_safe_term_matrices_used = false,
        generalized_overlap_solve_materialized = false,
        one_body_operator_materialized = true,
        hamiltonian_data_materialized = false,
        h1_solve_materialized = false,
        ida_data_materialized = false,
        density_density_materialized = false,
        rhf_materialized = false,
        driver_route_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
        metadata = term_metadata,
    )
end

function pqs_complete_core_shell_final_one_electron_hamiltonian(
    final_kinetic_result::NamedTuple,
    final_nuclear_by_center_results;
    symmetry_atol::Real = 1.0e-8,
)
    kinetic = _pqs_complete_core_shell_validate_final_one_body(
        final_kinetic_result,
        :kinetic;
        symmetry_atol,
    )
    nuclear_results = Tuple(final_nuclear_by_center_results)
    isempty(nuclear_results) &&
        throw(ArgumentError("complete core/shell H1 requires at least one nuclear center"))

    dimension = size(kinetic, 1)
    charged_nuclear_sum = zeros(Float64, dimension, dimension)
    center_summaries = map(nuclear_results) do result
        nuclear = _pqs_complete_core_shell_validate_final_one_body(
            result,
            :electron_nuclear_by_center;
            symmetry_atol,
        )
        metadata = result.metadata
        get(metadata, :nuclear_charge_recorded, false) ||
            throw(ArgumentError("complete core/shell H1 requires recorded nuclear charge"))
        get(metadata, :nuclear_charge_applied, true) &&
            throw(ArgumentError("complete core/shell H1 requires uncharged nuclear input"))
        get(metadata, :centers_summed, true) &&
            throw(ArgumentError("complete core/shell H1 requires separated nuclear centers"))
        charge = Float64(metadata.nuclear_charge)
        charged_nuclear_sum .+= charge .* nuclear
        (;
            center_key = get(metadata, :center_key, :unknown),
            center_index = get(metadata, :center_index, :unknown),
            center_location = get(metadata, :center_location, nothing),
            nuclear_charge = charge,
            nuclear_charge_recorded = true,
            input_nuclear_charge_applied = false,
            input_centers_summed = false,
            uncharged_by_center_convention =
                get(metadata, :uncharged_by_center_convention, false),
        )
    end

    hamiltonian_matrix = kinetic + charged_nuclear_sum
    hamiltonian_symmetry_error =
        norm(hamiltonian_matrix - transpose(hamiltonian_matrix), Inf)
    hamiltonian_symmetry_error <= Float64(symmetry_atol) ||
        throw(ArgumentError("complete core/shell H1 Hamiltonian must be symmetric"))
    return (;
        object_kind = :pqs_complete_core_shell_final_one_electron_hamiltonian,
        status = :materialized_pqs_complete_core_shell_final_one_electron_hamiltonian,
        blocker = nothing,
        final_dimension = dimension,
        kinetic_matrix = kinetic,
        charged_nuclear_matrix = charged_nuclear_sum,
        hamiltonian_matrix,
        hamiltonian_matrix_shape = size(hamiltonian_matrix),
        hamiltonian_matrix_finite = all(isfinite, hamiltonian_matrix),
        hamiltonian_matrix_symmetry_error = hamiltonian_symmetry_error,
        center_count = length(center_summaries),
        center_summaries,
        nuclear_charges_applied = true,
        centers_summed = true,
        hamiltonian_data_materialized = true,
        h1_solve_materialized = false,
        generalized_overlap_solve_materialized = false,
        ida_data_materialized = false,
        density_density_materialized = false,
        rhf_materialized = false,
        driver_route_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
        metadata = (;
            source = :pqs_complete_core_shell_final_one_electron_hamiltonian,
            nuclear_charge_application_stage = :hamiltonian_assembly,
            center_summation_stage = :hamiltonian_assembly,
            old_fixed_block_matrix_authority_used = false,
            current_route_safe_term_matrices_used = false,
        ),
    )
end

function pqs_complete_core_shell_final_h1_solve(
    hamiltonian_result::NamedTuple;
    reference_energy::Real = -0.5,
)
    get(hamiltonian_result, :object_kind, nothing) ===
        :pqs_complete_core_shell_final_one_electron_hamiltonian ||
        throw(ArgumentError("complete core/shell H1 solve requires complete H1 Hamiltonian"))
    matrix = Matrix{Float64}(hamiltonian_result.hamiltonian_matrix)
    all(isfinite, matrix) ||
        throw(ArgumentError("complete core/shell H1 matrix contains non-finite entries"))
    eigenvalues = eigvals(Symmetric((matrix + transpose(matrix)) ./ 2))
    energy = first(eigenvalues)
    return (;
        object_kind = :pqs_complete_core_shell_final_h1_solve,
        status = :materialized_pqs_complete_core_shell_final_h1_solve,
        blocker = nothing,
        solve_kind = :ordinary_symmetric,
        final_dimension = size(matrix, 1),
        lowest_energy = energy,
        reference_energy = Float64(reference_energy),
        reference_error = energy - Float64(reference_energy),
        eigenvalue_min = minimum(eigenvalues),
        eigenvalue_max = maximum(eigenvalues),
        h1_solve_materialized = true,
        generalized_overlap_solve_materialized = false,
        old_fixed_block_matrix_authority_used = false,
        current_route_safe_term_matrices_used = false,
        ida_data_materialized = false,
        density_density_materialized = false,
        rhf_materialized = false,
        driver_route_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
        metadata = (;
            source = :pqs_complete_core_shell_final_h1_solve,
            ordinary_symmetric_final_solve = true,
            old_fixed_block_matrix_authority_used = false,
            current_route_safe_term_matrices_used = false,
        ),
    )
end

"""
    pqs_complete_core_shell_final_ida_weights(final_basis, support_weights)

Project support-row integral weights into the complete core/shell final basis.
The support row order must match
`final_basis.support_row_order == :core_then_shell`.
"""
function pqs_complete_core_shell_final_ida_weights(
    final_basis::NamedTuple,
    support_weights;
    near_zero_atol::Real = 1.0e-12,
    metadata = (;),
)
    _pqs_complete_core_shell_validate_final_basis(final_basis)
    weights = Float64[Float64(weight) for weight in support_weights]
    support_count = final_basis.core_support_count + final_basis.shell_support_count
    length(weights) == support_count ||
        throw(DimensionMismatch("complete core/shell support weight length mismatch"))
    all(isfinite, weights) ||
        throw(ArgumentError("complete core/shell support weights contain non-finite entries"))
    coefficients = Matrix{Float64}(final_basis.final_coefficients)
    size(coefficients, 1) == support_count ||
        throw(DimensionMismatch("complete core/shell final coefficient row mismatch"))
    final_weights = vec(transpose(coefficients) * weights)
    all(isfinite, final_weights) ||
        throw(ArgumentError("complete core/shell final IDA weights contain non-finite entries"))
    near_zero_threshold = Float64(near_zero_atol)
    near_zero_count = count(weight -> abs(weight) <= near_zero_threshold, final_weights)
    negative_count = count(<(0.0), final_weights)
    positive_count = count(>(0.0), final_weights)
    return (;
        object_kind = :pqs_complete_core_shell_final_ida_weights,
        status = :materialized_pqs_complete_core_shell_final_ida_weights,
        blocker = nothing,
        final_basis_object_kind = final_basis.object_kind,
        final_basis_status = final_basis.status,
        support_row_order = final_basis.support_row_order,
        support_weight_count = length(weights),
        final_ida_weights = final_weights,
        final_ida_weight_count = length(final_weights),
        final_retained_count = final_basis.final_retained_count,
        weight_min = minimum(final_weights),
        weight_max = maximum(final_weights),
        weight_sum = sum(final_weights),
        support_weight_sum = sum(weights),
        near_zero_atol = near_zero_threshold,
        near_zero_count,
        negative_count,
        positive_count,
        all_finite = true,
        final_ida_weights_materialized = true,
        old_fixed_block_weight_authority_used = false,
        raw_source_weights_used_as_final_weights = false,
        boundary_shell_diagnostic_weights_used_as_final_weights = false,
        density_density_materialized = false,
        rhf_materialized = false,
        driver_route_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
        metadata = merge(
            NamedTuple(metadata),
            (;
                source = :pqs_complete_core_shell_final_ida_weights,
                support_weight_projection =
                    :transpose_final_coefficients_times_support_weights,
                support_row_order = final_basis.support_row_order,
                old_fixed_block_weight_authority_used = false,
            ),
        ),
    )
end

function _pqs_complete_core_shell_validate_final_basis(final_basis::NamedTuple)
    get(final_basis, :object_kind, nothing) === :pqs_complete_core_shell_final_basis ||
        throw(ArgumentError("complete core/shell one-body transfer requires complete final basis"))
    get(final_basis, :status, nothing) ===
        :available_pqs_complete_core_shell_final_basis ||
        throw(ArgumentError("complete core/shell one-body transfer requires available final basis"))
    get(final_basis, :final_basis_materialized, false) ||
        throw(ArgumentError("complete core/shell final basis is not materialized"))
    get(final_basis, :support_row_order, nothing) === :core_then_shell ||
        throw(ArgumentError("complete core/shell final basis support row order mismatch"))
    return nothing
end

function _pqs_complete_core_shell_validate_final_one_body(
    result::NamedTuple,
    term::Symbol;
    symmetry_atol::Real,
)
    get(result, :object_kind, nothing) ===
        :pqs_complete_core_shell_final_one_body_matrix ||
        throw(ArgumentError("complete core/shell Hamiltonian requires final one-body input"))
    get(result, :status, nothing) ===
        :materialized_pqs_complete_core_shell_final_one_body_matrix ||
        throw(ArgumentError("complete core/shell Hamiltonian requires materialized one-body input"))
    get(result, :term, nothing) === term ||
        throw(ArgumentError("complete core/shell Hamiltonian term mismatch"))
    matrix = Matrix{Float64}(result.final_operator)
    all(isfinite, matrix) ||
        throw(ArgumentError("complete core/shell final one-body matrix contains non-finite entries"))
    norm(matrix - transpose(matrix), Inf) <= Float64(symmetry_atol) ||
        throw(ArgumentError("complete core/shell final one-body matrix must be symmetric"))
    return matrix
end

function _pqs_complete_core_shell_one_body_metadata(
    term::Symbol,
    center_record,
    metadata,
)
    base = merge(
        NamedTuple(metadata),
        (;
            source = :pqs_complete_core_shell_final_one_body_matrix,
            route_owned_product_operator_used = true,
            old_fixed_block_matrix_authority_used = false,
            current_route_safe_term_matrices_used = false,
        ),
    )
    term !== :electron_nuclear_by_center && return base
    isnothing(center_record) &&
        throw(ArgumentError("electron-nuclear by-center transfer requires a center record"))
    charge = _pqs_complete_core_shell_descriptor_property(center_record, :charge)
    isnothing(charge) &&
        (charge = _pqs_complete_core_shell_descriptor_property(center_record, :nuclear_charge))
    location = _pqs_complete_core_shell_descriptor_property(center_record, :location)
    center_key = _pqs_complete_core_shell_descriptor_property(center_record, :center_key)
    center_index = _pqs_complete_core_shell_descriptor_property(center_record, :center_index)
    isnothing(charge) &&
        throw(ArgumentError("electron-nuclear by-center transfer requires recorded charge"))
    return merge(
        base,
        (;
            by_center = true,
            center_key = isnothing(center_key) ? :unknown : center_key,
            center_index = isnothing(center_index) ? :unknown : center_index,
            center_location = location,
            nuclear_charge = Float64(charge),
            nuclear_charge_recorded = true,
            nuclear_charge_applied = false,
            centers_summed = false,
            center_summation = false,
            uncharged_by_center_convention = true,
        ),
    )
end

function _pqs_complete_core_shell_descriptor_property(object, name::Symbol)
    if object isa NamedTuple && haskey(object, name)
        return getfield(object, name)
    end
    hasproperty(object, name) && return getproperty(object, name)
    return nothing
end
