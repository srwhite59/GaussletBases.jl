# CPBM-owned PQS final-basis helpers that still depend on CPBM result types.
#
# Shell-realization final-basis construction, oracle shell-support projection,
# and retained-boundary overlap/kinetic transfer are owned by
# CartesianFinalBasisRealization and reexported through CPBM aliases.

function _pqs_shell_final_basis_validate_available(final_basis::NamedTuple)
    get(final_basis, :object_kind, nothing) ===
        :pqs_source_shell_realization_final_basis ||
        throw(ArgumentError("PQS final helper requires a PQS shell-realization final basis"))
    get(final_basis, :status, nothing) ===
        :available_pqs_shell_realization_final_basis ||
        throw(ArgumentError("PQS final helper requires an available final basis"))
    get(final_basis, :final_basis_materialized, false) ||
        throw(ArgumentError("PQS final helper requires a materialized final basis"))
    return nothing
end

"""
    pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block(
        final_basis,
        retained_boundary_result,
    )

Transform a proven PQS retained-source electron-nuclear by-center boundary
block into the shell-realized final basis:

```text
V_final(center) = L' * V_boundary(center) * L
```

This is a by-center, uncharged one-body materialization seam only. Nuclear
charge application and center summation remain Hamiltonian-stage work.
"""
function pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block(
    final_basis::NamedTuple,
    retained_boundary_result::PairBlockMaterializationResult;
    symmetry_atol::Real = 1.0e-8,
)
    _pqs_shell_final_basis_validate_available(final_basis)
    nuclear_metadata =
        _pqs_shell_final_nuclear_validate_input(retained_boundary_result)

    boundary_operator = Matrix{Float64}(retained_boundary_result.block)
    boundary_count = final_basis.boundary_source_mode_count
    size(boundary_operator) == (boundary_count, boundary_count) ||
        throw(DimensionMismatch("PQS boundary electron-nuclear block shape must match boundary source modes"))
    all(isfinite, boundary_operator) ||
        throw(ArgumentError("PQS boundary electron-nuclear block contains non-finite entries"))
    boundary_symmetry_error =
        norm(boundary_operator - transpose(boundary_operator), Inf)
    boundary_symmetry_error <= Float64(symmetry_atol) ||
        throw(ArgumentError("PQS boundary electron-nuclear block must be symmetric"))

    cleanup = final_basis.lowdin_cleanup
    final_operator = transpose(cleanup) * boundary_operator * cleanup
    final_symmetry_error = norm(final_operator - transpose(final_operator), Inf)

    return (;
        object_kind =
            :pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block,
        status =
            :materialized_pqs_shell_final_electron_nuclear_by_center_from_boundary_block,
        blocker = nothing,
        term = :electron_nuclear_by_center,
        input_term = retained_boundary_result.term,
        input_pair_key = retained_boundary_result.pair_key,
        input_block_space = nuclear_metadata.block_space,
        boundary_operator,
        boundary_operator_shape = size(boundary_operator),
        boundary_operator_finite = true,
        boundary_operator_symmetry_error = boundary_symmetry_error,
        final_operator,
        final_operator_shape = size(final_operator),
        final_operator_finite = all(isfinite, final_operator),
        final_operator_symmetry_error = final_symmetry_error,
        final_basis_object_kind = final_basis.object_kind,
        final_basis_status = final_basis.status,
        boundary_source_mode_count = final_basis.boundary_source_mode_count,
        final_retained_count = final_basis.final_retained_count,
        center_key = nuclear_metadata.center_key,
        center_index = nuclear_metadata.center_index,
        center_location = nuclear_metadata.center_location,
        nuclear_charge = nuclear_metadata.nuclear_charge,
        nuclear_charge_recorded = nuclear_metadata.nuclear_charge_recorded,
        nuclear_charge_applied = false,
        centers_summed = false,
        uncharged_by_center_convention = true,
        retained_boundary_operator_input_used = true,
        raw_source_operator_input_used = false,
        shell_support_operator_generated = false,
        shell_realization_materialized = true,
        lowdin_cleanup_used = true,
        one_body_operator_materialized = true,
        electron_nuclear_materialized = true,
        charge_summing_materialized = false,
        h1_solve_materialized = false,
        hamiltonian_data_materialized = false,
        ida_data_materialized = false,
        density_density_materialized = false,
        rhf_materialized = false,
        driver_route_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
        next_blocker = :missing_pqs_final_one_electron_hamiltonian_assembly,
        metadata = (;
            source =
                :pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block,
            boundary_operator_contract =
                :retained_pqs_source_modes_equal_shell_projected_boundary_operator,
            by_center = true,
            nuclear_charge_recorded = nuclear_metadata.nuclear_charge_recorded,
            nuclear_charge_applied = false,
            centers_summed = false,
            uncharged_by_center_convention = true,
            current_route_safe_term_matrices_used = false,
            old_fixed_block_matrix_authority_used = false,
        ),
    )
end

function _pqs_shell_final_nuclear_validate_input(
    result::PairBlockMaterializationResult,
)
    result.term === :retained_source_electron_nuclear_by_center ||
        throw(ArgumentError("PQS final nuclear helper requires a retained-source electron-nuclear by-center result"))
    result.materialized ||
        throw(ArgumentError("PQS final nuclear helper requires a materialized retained-source block"))
    metadata = result.metadata
    get(metadata, :block_space, nothing) === :retained_pqs_source_modes ||
        throw(ArgumentError("PQS final nuclear helper requires retained PQS source-mode block space"))
    get(metadata, :by_center, false) ||
        throw(ArgumentError("PQS final nuclear helper requires by-center metadata"))
    get(metadata, :nuclear_charge_recorded, false) ||
        throw(ArgumentError("PQS final nuclear helper requires recorded nuclear charge metadata"))
    get(metadata, :nuclear_charge_applied, false) &&
        throw(ArgumentError("PQS final nuclear helper requires uncharged by-center input"))
    get(metadata, :centers_summed, false) &&
        throw(ArgumentError("PQS final nuclear helper requires separated by-center input"))
    get(metadata, :uncharged_by_center_convention, false) ||
        throw(ArgumentError("PQS final nuclear helper requires uncharged by-center convention metadata"))
    return (;
        block_space = metadata.block_space,
        center_key = metadata.center_key,
        center_index = metadata.center_index,
        center_location = metadata.center_location,
        nuclear_charge = metadata.nuclear_charge,
        nuclear_charge_recorded = metadata.nuclear_charge_recorded,
        nuclear_charge_applied = metadata.nuclear_charge_applied,
        centers_summed = metadata.centers_summed,
        uncharged_by_center_convention = metadata.uncharged_by_center_convention,
    )
end

"""
    pqs_source_shell_final_one_electron_hamiltonian(
        final_kinetic_result,
        final_nuclear_by_center_results,
    )

Assemble the PQS shell-realized final one-electron Hamiltonian:

```text
H = T_final + sum_center Z_center * V_final_unit_charge(center)
```

The nuclear inputs must be separated, uncharged by-center matrices. This helper
is the first Hamiltonian-stage consumer: it applies recorded charges and sums
centers, but it does not solve H1 or materialize downstream density/driver data.
"""
function pqs_source_shell_final_one_electron_hamiltonian(
    final_kinetic_result::NamedTuple,
    final_nuclear_by_center_results;
    symmetry_atol::Real = 1.0e-8,
)
    kinetic = _pqs_shell_final_hamiltonian_validate_kinetic(
        final_kinetic_result;
        symmetry_atol,
    )
    nuclear_results =
        _pqs_shell_final_hamiltonian_nuclear_tuple(final_nuclear_by_center_results)
    isempty(nuclear_results) &&
        throw(ArgumentError("PQS final one-electron Hamiltonian requires at least one nuclear center"))

    dimension = size(kinetic, 1)
    charged_nuclear_sum = zeros(Float64, dimension, dimension)
    center_summaries = map(nuclear_results) do result
        nuclear_matrix, summary =
            _pqs_shell_final_hamiltonian_validate_nuclear(
                result,
                dimension;
                symmetry_atol,
            )
        charged_nuclear_sum .+= summary.nuclear_charge .* nuclear_matrix
        summary
    end

    hamiltonian_matrix = kinetic + charged_nuclear_sum
    hamiltonian_symmetry_error =
        norm(hamiltonian_matrix - transpose(hamiltonian_matrix), Inf)
    hamiltonian_symmetry_error <= Float64(symmetry_atol) ||
        throw(ArgumentError("PQS final one-electron Hamiltonian must be symmetric"))

    return (;
        object_kind = :pqs_source_shell_final_one_electron_hamiltonian,
        status = :materialized_pqs_shell_final_one_electron_hamiltonian,
        blocker = nothing,
        term = :one_electron_hamiltonian,
        final_dimension = dimension,
        kinetic_matrix = kinetic,
        kinetic_matrix_shape = size(kinetic),
        kinetic_matrix_finite = true,
        charged_nuclear_matrix = charged_nuclear_sum,
        charged_nuclear_matrix_shape = size(charged_nuclear_sum),
        charged_nuclear_matrix_finite = all(isfinite, charged_nuclear_sum),
        hamiltonian_matrix,
        hamiltonian_matrix_shape = size(hamiltonian_matrix),
        hamiltonian_matrix_finite = all(isfinite, hamiltonian_matrix),
        hamiltonian_matrix_symmetry_error = hamiltonian_symmetry_error,
        center_count = length(center_summaries),
        center_summaries,
        applied_nuclear_charges =
            Tuple(summary.nuclear_charge for summary in center_summaries),
        centers_summed = true,
        nuclear_charges_applied = true,
        input_nuclear_charges_applied = false,
        input_centers_summed = false,
        one_body_operator_materialized = true,
        electron_nuclear_materialized = true,
        hamiltonian_data_materialized = true,
        h1_solve_materialized = false,
        eigensolve_materialized = false,
        generalized_overlap_solve_materialized = false,
        ida_data_materialized = false,
        density_density_materialized = false,
        rhf_materialized = false,
        driver_route_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
        next_blocker = :missing_pqs_final_one_electron_solve,
        metadata = (;
            source = :pqs_source_shell_final_one_electron_hamiltonian,
            nuclear_charge_application_stage = :hamiltonian_assembly,
            center_summation_stage = :hamiltonian_assembly,
            current_route_safe_term_matrices_used = false,
            old_fixed_block_matrix_authority_used = false,
        ),
    )
end

function _pqs_shell_final_hamiltonian_validate_kinetic(
    result::NamedTuple;
    symmetry_atol::Real,
)
    get(result, :object_kind, nothing) ===
        :pqs_source_shell_final_one_body_from_boundary_matrix ||
        throw(ArgumentError("PQS final Hamiltonian requires final-boundary one-body kinetic input"))
    get(result, :status, nothing) ===
        :materialized_pqs_shell_final_one_body_from_boundary_matrix ||
        throw(ArgumentError("PQS final Hamiltonian requires materialized kinetic input"))
    get(result, :term, nothing) === :kinetic ||
        throw(ArgumentError("PQS final Hamiltonian requires kinetic term input"))
    get(result, :one_body_operator_materialized, false) ||
        throw(ArgumentError("PQS final Hamiltonian requires a materialized kinetic matrix"))
    haskey(result, :final_operator) ||
        throw(ArgumentError("PQS final Hamiltonian kinetic input requires final_operator"))
    kinetic = Matrix{Float64}(result.final_operator)
    size(kinetic, 1) == size(kinetic, 2) ||
        throw(DimensionMismatch("PQS final Hamiltonian kinetic matrix must be square"))
    all(isfinite, kinetic) ||
        throw(ArgumentError("PQS final Hamiltonian kinetic matrix contains non-finite entries"))
    norm(kinetic - transpose(kinetic), Inf) <= Float64(symmetry_atol) ||
        throw(ArgumentError("PQS final Hamiltonian kinetic matrix must be symmetric"))
    return kinetic
end

function _pqs_shell_final_hamiltonian_nuclear_tuple(result::NamedTuple)
    return (result,)
end

function _pqs_shell_final_hamiltonian_nuclear_tuple(results)
    return Tuple(results)
end

function _pqs_shell_final_hamiltonian_validate_nuclear(
    result::NamedTuple,
    dimension::Int;
    symmetry_atol::Real,
)
    get(result, :object_kind, nothing) ===
        :pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block ||
        throw(ArgumentError("PQS final Hamiltonian requires final by-center nuclear input"))
    get(result, :status, nothing) ===
        :materialized_pqs_shell_final_electron_nuclear_by_center_from_boundary_block ||
        throw(ArgumentError("PQS final Hamiltonian requires materialized nuclear input"))
    get(result, :term, nothing) === :electron_nuclear_by_center ||
        throw(ArgumentError("PQS final Hamiltonian requires electron-nuclear by-center input"))
    get(result, :nuclear_charge_recorded, false) ||
        throw(ArgumentError("PQS final Hamiltonian requires recorded nuclear charge"))
    get(result, :nuclear_charge_applied, true) &&
        throw(ArgumentError("PQS final Hamiltonian requires uncharged nuclear input"))
    get(result, :centers_summed, true) &&
        throw(ArgumentError("PQS final Hamiltonian requires separated by-center nuclear input"))
    get(result, :uncharged_by_center_convention, false) ||
        throw(ArgumentError("PQS final Hamiltonian requires uncharged by-center convention"))
    nuclear = Matrix{Float64}(result.final_operator)
    size(nuclear) == (dimension, dimension) ||
        throw(DimensionMismatch("PQS final Hamiltonian nuclear matrix shape mismatch"))
    all(isfinite, nuclear) ||
        throw(ArgumentError("PQS final Hamiltonian nuclear matrix contains non-finite entries"))
    symmetry_error = norm(nuclear - transpose(nuclear), Inf)
    symmetry_error <= Float64(symmetry_atol) ||
        throw(ArgumentError("PQS final Hamiltonian nuclear matrix must be symmetric"))
    summary = (;
        center_key = result.center_key,
        center_index = result.center_index,
        center_location = result.center_location,
        nuclear_charge = Float64(result.nuclear_charge),
        nuclear_charge_recorded = true,
        input_nuclear_charge_applied = false,
        input_centers_summed = false,
        uncharged_by_center_convention = true,
        matrix_shape = size(nuclear),
        matrix_symmetry_error = symmetry_error,
    )
    return nuclear, summary
end
