function _pqs_source_box_route_driver_atomic_bundle(
    axis_basis;
    expansion,
    backend::Symbol,
)
    return _mapped_ordinary_gausslet_1d_bundle(
        axis_basis;
        exponents = expansion.exponents,
        center = 0.0,
        backend,
        refinement_levels = 0,
    )
end

function _pqs_source_box_route_driver_atomic_multilayer_shell_records(sequence)
    records = NamedTuple[]
    for (layer_index, shell) in enumerate(sequence.shell_layers)
        support_indices = Vector{Int}(shell.support_indices)
        local_coefficients =
            _nested_support_coefficient_slice(
                shell.coefficient_matrix,
                support_indices,
            )
        push!(
            records,
            (;
                layer_index,
                shell_support_indices = support_indices,
                shell_support_states = shell.support_states,
                shell_final_coefficients = local_coefficients,
                support_count = length(support_indices),
                retained_count = size(local_coefficients, 2),
                provenance = :one_center_atomic_full_parent_shell_sequence,
            ),
        )
    end
    return Tuple(records)
end

function _pqs_source_box_route_driver_atomic_multilayer_plan(sequence, bundles)
    shell_records =
        _pqs_source_box_route_driver_atomic_multilayer_shell_records(sequence)
    shell_support_indices =
        reduce(vcat, (record.shell_support_indices for record in shell_records); init = Int[])
    shell_support_states =
        reduce(
            vcat,
            (record.shell_support_states for record in shell_records);
            init = Tuple{Int,Int,Int}[],
        )
    shell_final_coefficients =
        _pqs_multilayer_block_concatenate_shell_coefficients(shell_records)
    shell_duplicate_count =
        length(shell_support_indices) - length(unique(shell_support_indices))
    core_shell_duplicate_count =
        length(intersect(sequence.core_indices, shell_support_indices))
    shell_duplicate_count == 0 ||
        throw(ArgumentError("atomic common operator plan requires disjoint shell support"))
    core_shell_duplicate_count == 0 ||
        throw(ArgumentError("atomic common operator plan requires disjoint core/shell support"))
    support_indices = vcat(sequence.core_indices, shell_support_indices)
    if !isnothing(sequence.support_indices)
        sort!(copy(support_indices)) == sort!(copy(sequence.support_indices)) ||
            throw(ArgumentError("atomic common operator plan support indices do not match sequence support"))
    end
    core_count = length(sequence.core_indices)
    shell_retained_count = size(shell_final_coefficients, 2)
    core_count + shell_retained_count == size(sequence.coefficient_matrix, 2) ||
        throw(ArgumentError("atomic common operator plan retained count mismatch"))
    size(shell_final_coefficients, 1) == length(shell_support_indices) ||
        throw(ArgumentError("atomic common operator plan shell coefficient row mismatch"))
    metrics = _pqs_multilayer_axis_metrics(bundles)
    summary = (;
        status = :available_pqs_multilayer_shell_source_plan,
        blocker = nothing,
        source_plan_family = :one_center_atomic_complete_core_shell,
        layer_count = length(shell_records),
        core_support_count = core_count,
        shell_support_count = length(shell_support_indices),
        shell_final_retained_count = shell_retained_count,
        support_counts = (core_count, map(record -> record.support_count, shell_records)...),
        retained_counts = (core_count, map(record -> record.retained_count, shell_records)...),
        collapsed_shell_sector = true,
        final_basis_helper = :pqs_complete_core_shell_final_basis,
    )
    return (;
        object_kind = :pqs_multilayer_shell_source_plan,
        status = :available_pqs_multilayer_shell_source_plan,
        blocker = nothing,
        source_kind = :one_center_atomic_complete_core_shell_adapter,
        bundles,
        metrics,
        core_box = nothing,
        outer_box = nothing,
        bond_axis = nothing,
        layer_count = length(shell_records),
        shell_records,
        core_support_indices = Vector{Int}(sequence.core_indices),
        core_support_states = sequence.core_states,
        shell_support_indices,
        shell_support_states,
        shell_final_coefficients,
        summary,
        metadata = (;
            source =
                :pqs_source_box_route_driver_atomic_multilayer_plan,
            input_source = :one_center_atomic_full_parent_shell_sequence,
            route_owned_authority = true,
            collapsed_shell_sector = true,
        ),
    )
end

function _pqs_source_box_route_driver_atomic_axis_layers(bundle)
    auxiliary_layer = bundle.pgdg_intermediate.auxiliary_layer
    return (x = auxiliary_layer, y = auxiliary_layer, z = auxiliary_layer)
end

function _pqs_source_box_route_driver_atomic_center_records(system, z::Real)
    location =
        hasproperty(system, :atom_locations) && !isempty(system.atom_locations) ?
        system.atom_locations[1] :
        (0.0, 0.0, 0.0)
    return ((;
        center_key = :atom_1,
        center_index = 1,
        location,
        charge = Float64(z),
        nuclear_charge = Float64(z),
    ),)
end

function _pqs_source_box_route_driver_centered_factor_terms(pgdg, expansion, center)
    center == pgdg.center && return pgdg.gaussian_factor_terms
    ops = mapped_ordinary_one_body_operators(
        pgdg.basis; exponents = expansion.exponents, center, backend = pgdg.backend)
    return ops.gaussian_factors
end

function _pqs_source_box_route_driver_validate_pgdg_expansion(pgdg, expansion)
    length(pgdg) == 3 || throw(DimensionMismatch("PQS IDA requires three PGDG axes"))
    for (axis_name, axis) in zip((:x, :y, :z), pgdg)
        length(axis.exponents) == length(expansion) &&
            axis.exponents == expansion.exponents ||
            throw(ArgumentError("PQS $(axis_name)-axis Coulomb exponent sequence mismatch"))
    end
    return nothing
end

function _pqs_source_box_route_driver_terminal_products(terminal_basis_realization, pgdg)
    S = Tuple(axis.overlap for axis in pgdg)
    T = Tuple(axis.kinetic for axis in pgdg)
    n = terminal_basis_realization.final_dimension
    K = zeros(Float64, n, n)
    C = CartesianFinalBasisRealization
    C.assemble_terminal_product_operator!(K, terminal_basis_realization, T[1], S[2], S[3])
    C.assemble_terminal_product_operator!(K, terminal_basis_realization, S[1], T[2], S[3])
    C.assemble_terminal_product_operator!(K, terminal_basis_realization, S[1], S[2], T[3])
    return (; kinetic = K)
end

function _pqs_source_box_route_driver_terminal_unit_nuclear(terminal_basis_realization, expansion, pgdg, atom_locations::Vector{NTuple{3,Float64}})
    _pqs_source_box_route_driver_validate_pgdg_expansion(pgdg, expansion)
    n = terminal_basis_realization.final_dimension
    C = CartesianFinalBasisRealization
    U = Matrix{Float64}[]
    for location in atom_locations
        matrix = zeros(Float64, n, n)
        factors = ntuple(axis -> _pqs_source_box_route_driver_centered_factor_terms(
            pgdg[axis], expansion, location[axis]), 3)
        C._accumulate_terminal_gaussian_sum!(
            matrix, terminal_basis_realization, expansion.coefficients,
            factors[1], factors[2], factors[3])
        push!(U, matrix)
    end
    return U
end

function _pqs_source_box_route_driver_terminal_vee(terminal_basis_realization, expansion, pgdg)
    _pqs_source_box_route_driver_validate_pgdg_expansion(pgdg, expansion)
    V = zeros(Float64, terminal_basis_realization.final_dimension,
        terminal_basis_realization.final_dimension)
    CartesianFinalBasisRealization.assemble_terminal_ida_interaction!(
        V, terminal_basis_realization, expansion.coefficients,
        pgdg[1].pair_factor_terms_raw, pgdg[2].pair_factor_terms_raw,
        pgdg[3].pair_factor_terms_raw,
        pgdg[1].weights, pgdg[2].weights, pgdg[3].weights)
    return V
end
