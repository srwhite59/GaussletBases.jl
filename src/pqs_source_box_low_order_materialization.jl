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

function _pqs_source_box_route_driver_atomic_common_h1(
    axis_basis,
    system;
    expansion,
    backend::Symbol,
    nside::Int,
    z::Real,
)
    bundle = _pqs_source_box_route_driver_atomic_bundle(
        axis_basis;
        expansion,
        backend,
    )
    term_coefficients = _one_center_atomic_term_coefficients(bundle, expansion)
    n = length(bundle.basis)
    sequence = _build_one_center_atomic_shell_sequence(
        bundle,
        (1:n, 1:n, 1:n);
        nside,
        term_coefficients,
    )
    bundles = _CartesianNestedAxisBundles3D(bundle, bundle, bundle)
    plan = _pqs_source_box_route_driver_atomic_multilayer_plan(sequence, bundles)
    final_basis =
        pqs_multilayer_complete_core_shell_final_basis(
            plan;
            metadata = (;
                source =
                    :pqs_source_box_route_driver_atomic_common_h1,
                route_kind = :one_center_atomic_complete_core_shell,
            ),
        )
    final_basis.status === :available_pqs_complete_core_shell_final_basis ||
        throw(ArgumentError("atomic common final basis did not materialize"))
    h1_payload = pqs_multilayer_complete_core_shell_h1_payload(
        plan;
        final_basis,
        coulomb_expansion = expansion,
        center_records =
            _pqs_source_box_route_driver_atomic_center_records(system, z),
        axis_layers = _pqs_source_box_route_driver_atomic_axis_layers(bundle),
        metadata = (;
            source = :pqs_source_box_route_driver_atomic_common_h1,
            h1_operator_authority =
                :pqs_multilayer_complete_core_shell_h1_payload,
        ),
    )
    return (; bundle, sequence, plan, final_basis, h1_payload)
end

function _pqs_source_box_route_driver_atomic_artifact_sidecar(common)
    packet = common.sequence.packet
    isnothing(packet) && throw(
        ArgumentError("atomic artifact sidecar requires an assembled shell sequence packet"),
    )
    final_basis = common.final_basis
    final_coefficients = Matrix{Float64}(final_basis.final_coefficients)
    support_count, final_count = size(final_coefficients)
    length(packet.weights) == support_count ||
        throw(DimensionMismatch("atomic artifact sidecar weight dimension mismatch"))

    function project_matrix(matrix, name)
        local_matrix = Matrix{Float64}(matrix)
        size(local_matrix) == (support_count, support_count) ||
            throw(DimensionMismatch("atomic artifact sidecar $(name) dimension mismatch"))
        return transpose(final_coefficients) * local_matrix * final_coefficients
    end

    weights = vec(transpose(final_coefficients) * Float64.(packet.weights))
    length(weights) == final_count ||
        throw(DimensionMismatch("atomic artifact sidecar final weight dimension mismatch"))
    position_x = project_matrix(packet.position_x, "position_x")
    position_y = project_matrix(packet.position_y, "position_y")
    position_z = project_matrix(packet.position_z, "position_z")
    pair_sum =
        isnothing(packet.pair_sum) ? nothing :
        project_matrix(packet.pair_sum, "pair_sum")
    return (;
        weights,
        fixed_centers = hcat(
            diag(position_x),
            diag(position_y),
            diag(position_z),
        ),
        pair_sum,
    )
end

function _pqs_source_box_route_driver_wl_atomic_pure_gausslet_materialization(
    report;
    save_basis_artifact::Bool,
    save_ham_artifact::Bool,
    basisfile,
    hamfile,
    materializer_backend,
    materializer_nside,
    white_lindsey_expansion,
    white_lindsey_Z,
)
    basis =
        hasproperty(report, :parent_basis_object) ? report.parent_basis_object : nothing
    isnothing(basis) && return nothing
    parent_axes = CartesianParentGaussletBases.parent_axes(basis)
    axis_basis = parent_axes.x

    system = report.system_metadata
    recipe = report.recipe_metadata
    atom_count =
        hasproperty(system, :atom_symbols) ? length(system.atom_symbols) : nothing
    atom_count == 1 || return nothing
    supplement_policy =
        hasproperty(recipe, :supplement_policy) ? recipe.supplement_policy : nothing
    (isnothing(supplement_policy) || supplement_policy == :none) || return nothing

    backend =
        isnothing(materializer_backend) ?
        (
            hasproperty(system, :map_backend) ?
            system.map_backend :
            :pgdg_localized_experimental
        ) :
        materializer_backend
    nside =
        isnothing(materializer_nside) ?
        (
            hasproperty(recipe, :n_s) && !isnothing(recipe.n_s) ?
            recipe.n_s :
            recipe.q
        ) :
        materializer_nside
    z =
        isnothing(white_lindsey_Z) ?
        Float64(first(system.nuclear_charges)) :
        Float64(white_lindsey_Z)
    expansion =
        isnothing(white_lindsey_expansion) ?
        coulomb_gaussian_expansion(doacc = false) :
        white_lindsey_expansion

    common =
        _pqs_source_box_route_driver_atomic_common_h1(
            axis_basis,
            system;
            expansion,
            backend,
            nside,
            z,
        )
    final_basis = common.final_basis
    h1_payload = common.h1_payload
    h1_hamiltonian = h1_payload.final_hamiltonian
    nuclear_one_body = h1_hamiltonian.charged_nuclear_matrix
    h1 = h1_hamiltonian.hamiltonian_matrix
    dim = h1_payload.h1.final_dimension
    overlap_identity_error = final_basis.final_overlap_identity_error
    h1_symmetry_error = h1_hamiltonian.hamiltonian_matrix_symmetry_error
    h1_lowest = h1_payload.h1.lowest_energy
    artifact_sidecar =
        (save_basis_artifact || save_ham_artifact) ?
        _pqs_source_box_route_driver_atomic_artifact_sidecar(common) :
        nothing
    summary = (;
        route_family = report.route_family,
        geometry = :atomic,
        supplement_policy = :none,
        source = :pqs_multilayer_complete_core_shell_h1_payload,
        backend,
        nside,
        nuclear_charge = z,
        support_dimension =
            final_basis.core_support_count + final_basis.shell_support_count,
        retained_dimension = dim,
        coefficient_matrix_size = size(final_basis.final_coefficients),
        overlap_identity_error,
        h1_lowest,
        h1_finite = all(isfinite, h1),
        h1_symmetry_error,
        density_density_pair_sum_present =
            !isnothing(artifact_sidecar) && !isnothing(artifact_sidecar.pair_sum),
        h1_operator_authority = :pqs_multilayer_complete_core_shell_h1_payload,
    )
    basis_artifact_path = nothing
    if save_basis_artifact
        basis_artifact_path = String(basisfile)
        jldopen(basis_artifact_path, "w") do file
            file["artifact_kind"] = :white_lindsey_atomic_pure_gausslet_basis
            file["summary"] = summary
            file["coefficient_matrix"] = final_basis.final_coefficients
            file["support_indices"] =
                vcat(final_basis.core_support_indices, final_basis.shell_support_indices)
            file["overlap"] = final_basis.final_overlap
            file["weights"] = artifact_sidecar.weights
            file["fixed_centers"] = artifact_sidecar.fixed_centers
        end
    end
    ham_artifact_path = nothing
    if save_ham_artifact
        ham_artifact_path = String(hamfile)
        jldopen(ham_artifact_path, "w") do file
            file["artifact_kind"] = :white_lindsey_atomic_pure_gausslet_hamiltonian
            file["summary"] = summary
            file["overlap"] = final_basis.final_overlap
            file["kinetic"] = h1_hamiltonian.kinetic_matrix
            file["nuclear_one_body"] = nuclear_one_body
            file["one_body_hamiltonian"] = h1
            file["density_density_pair_sum"] = artifact_sidecar.pair_sum
            file["weights"] = artifact_sidecar.weights
            file["fixed_centers"] = artifact_sidecar.fixed_centers
        end
    end

    return (;
        route_family = report.route_family,
        result_kind = :white_lindsey_atomic_pure_gausslet,
        requested = true,
        materialized = true,
        summary,
        retained_dimension = dim,
        support_dimension = summary.support_dimension,
        final_dimension = dim,
        h1_lowest,
        h1_finite = summary.h1_finite,
        h1_symmetry_error,
        overlap_identity_error,
        save_basis_artifact_requested = save_basis_artifact,
        save_ham_artifact_requested = save_ham_artifact,
        basisfile = basis_artifact_path,
        hamfile = ham_artifact_path,
        basis_artifact_written = save_basis_artifact,
        ham_artifact_written = save_ham_artifact,
    )
end

function _pqs_source_box_route_driver_pqs_h2_residual_gto_materialization(
    report;
    save_basis_artifact::Bool,
    save_ham_artifact::Bool,
    basisfile,
    hamfile,
    residual_gto_provider_blocks::Symbol = :none,
)
    report.route_family === :pqs_source_box || return nothing
    report.route_kind === :bond_aligned_diatomic_independent_pqs_source_box_core_shell ||
        return nothing
    get(report.recipe_metadata, :supplement_policy, nothing) === :mwg_residual_gto ||
        return nothing
    inputs =
        hasproperty(report, :pqs_gto_sidecar_inputs) ?
        report.pqs_gto_sidecar_inputs :
        nothing
    isnothing(inputs) && throw(
        ArgumentError(
            "H2 PQS residual-GTO sidecar materialization requires source plan, final basis, H1, density interaction, and supplement representation",
        ),
    )

    provider_block_mode =
        _pqs_source_box_route_driver_pqs_h2_provider_block_mode(
            residual_gto_provider_blocks,
        )
    sidecar = _pqs_source_box_route_driver_pqs_gto_sidecar(inputs)
    route_metadata = _pqs_source_box_route_driver_pqs_h2_route_metadata(
        report,
        merge(inputs, (; provider_block_mode,)),
    )
    inputs =
        merge(inputs, (; route_metadata, provider_block_mode,))
    residual =
        provider_block_mode !== :none ?
        _pqs_source_box_route_driver_pqs_gto_residual_transform(sidecar) :
        nothing
    provider_packet =
        provider_block_mode !== :none ?
        _pqs_source_box_route_driver_pqs_h2_residual_gto_provider_packet(
            inputs,
            sidecar,
            residual,
        ) :
        nothing
    density_descriptor =
        provider_block_mode !== :none ?
        _pqs_source_box_route_driver_pqs_h2_residual_gto_density_descriptor(
            inputs,
            sidecar,
            residual,
        ) :
        nothing
    one_body_blocks =
        provider_block_mode !== :none ?
        _pqs_source_box_route_driver_pqs_gto_one_body_blocks(provider_packet) :
        nothing
    density_blocks =
        provider_block_mode === :one_body_and_density_provider ?
        _pqs_source_box_route_driver_pqs_h2_residual_gto_density_blocks(
            provider_packet,
            density_descriptor,
        ) :
        nothing
    ida_hamiltonian =
        provider_block_mode === :one_body_and_density_provider ?
        _pqs_source_box_route_driver_pqs_h2_residual_gto_ida_hamiltonian(
            one_body_blocks,
            density_blocks,
            route_metadata,
        ) :
        nothing
    ida_hamiltonian_summary =
        _pqs_source_box_route_driver_ida_hamiltonian_summary(ida_hamiltonian)
    final_basis = inputs.final_basis
    h1_hamiltonian = inputs.h1_hamiltonian
    h1 = inputs.h1
    h1_matrix = Matrix{Float64}(h1_hamiltonian.hamiltonian_matrix)
    h1_symmetry_error = norm(h1_matrix - transpose(h1_matrix), Inf)
    h1_finite = all(isfinite, h1_matrix)
    optional = _pqs_source_box_route_driver_optional_property
    augmented_dimension = optional(one_body_blocks, :augmented_dimension)
    summary = (;
        route_family = report.route_family,
        route_kind = report.route_kind,
        result_kind = :h2_pqs_residual_gto_ida_hamiltonian,
        final_dimension = final_basis.final_dimension,
        overlap_identity_error = final_basis.final_overlap_identity_error,
        h1_lowest = h1.lowest_energy,
        h1_finite,
        h1_symmetry_error,
        supplement_policy = report.recipe_metadata.supplement_policy,
        provider_block_mode,
        residual_rank = optional(residual, :residual_rank),
        residual_overlap_identity_error =
            optional(residual, :residual_overlap_identity_error),
        residual_overlap_cutoff = optional(residual, :residual_overlap_cutoff),
        augmented_dimension,
        augmented_density_space =
            optional(density_descriptor, :augmented_density_space),
        augmented_density_gauge =
            optional(density_descriptor, :density_gauge),
        ida_hamiltonian_summary...,
    )

    ham_artifact_path = nothing
    ham_artifact_reloaded = false
    if save_ham_artifact
        ham_artifact_path = String(hamfile)
        isnothing(ida_hamiltonian) &&
            throw(ArgumentError("H2 residual-GTO Ham artifact requires IDA Hamiltonian provider blocks"))
        write_cartesian_ida_hamiltonian(ham_artifact_path, ida_hamiltonian)
        _pqs_source_box_route_driver_ida_hamiltonian_smoke(
            read_cartesian_ida_hamiltonian(ham_artifact_path),
        )
        ham_artifact_reloaded = true
    end

    return (;
        route_family = report.route_family,
        result_kind = :h2_pqs_residual_gto_ida_hamiltonian,
        requested = true,
        materialized = true,
        summary,
        final_dimension = final_basis.final_dimension,
        retained_dimension = report.retained_dimension,
        h1_lowest = h1.lowest_energy,
        h1_finite,
        h1_symmetry_error,
        overlap_identity_error = final_basis.final_overlap_identity_error,
        provider_block_mode,
        residual_rank = summary.residual_rank,
        residual_overlap_identity_error = summary.residual_overlap_identity_error,
        ida_hamiltonian,
        ida_hamiltonian_summary...,
        augmented_dimension,
        save_ham_artifact_requested = save_ham_artifact,
        hamfile = ham_artifact_path,
        ham_artifact_written = save_ham_artifact,
        ham_artifact_reloaded,
    )
end

function _pqs_source_box_route_driver_materialization(
    report;
    materialize_route::Bool = false,
    save_basis_artifact::Bool = false,
    save_ham_artifact::Bool = false,
    basisfile = nothing,
    hamfile = nothing,
    materializer_backend = nothing,
    materializer_nside = nothing,
    white_lindsey_expansion = nothing,
    white_lindsey_Z = nothing,
    residual_gto_provider_blocks::Symbol = :none,
)
    requested = materialize_route || save_basis_artifact || save_ham_artifact
    retained_dimension =
        hasproperty(report, :retained_dimension) ? report.retained_dimension : nothing
    if requested && hasproperty(report, :route_family) &&
       report.route_family == :white_lindsey_low_order
        wl_materialization =
            _pqs_source_box_route_driver_wl_atomic_pure_gausslet_materialization(
                report;
                save_basis_artifact,
                save_ham_artifact,
                basisfile,
                hamfile,
                materializer_backend,
                materializer_nside,
                white_lindsey_expansion,
                white_lindsey_Z,
            )
        !isnothing(wl_materialization) && return wl_materialization
    end
    if requested && hasproperty(report, :route_family) &&
       report.route_family == :pqs_source_box
        pqs_materialization =
            _pqs_source_box_route_driver_pqs_h2_residual_gto_materialization(
                report;
                save_basis_artifact,
                save_ham_artifact,
                basisfile,
                hamfile,
                residual_gto_provider_blocks,
            )
        !isnothing(pqs_materialization) && return pqs_materialization
    end
    return (;
        route_family = hasproperty(report, :route_family) ? report.route_family : nothing,
        result_kind = requested ? :not_materialized : :not_requested,
        requested,
        materialized = false,
        retained_dimension,
        save_basis_artifact_requested = save_basis_artifact,
        save_ham_artifact_requested = save_ham_artifact,
        basisfile,
        hamfile,
    )
end
