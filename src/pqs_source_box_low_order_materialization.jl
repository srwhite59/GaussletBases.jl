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
    cleanup = Matrix{Float64}(final_basis.combined_lowdin_cleanup)
    pre_final_count, final_count = size(cleanup)
    length(packet.weights) == pre_final_count ||
        throw(DimensionMismatch("atomic artifact sidecar weight dimension mismatch"))

    function project_matrix(matrix, name)
        local_matrix = Matrix{Float64}(matrix)
        size(local_matrix) == (pre_final_count, pre_final_count) ||
            throw(DimensionMismatch("atomic artifact sidecar $(name) dimension mismatch"))
        return transpose(cleanup) * local_matrix * cleanup
    end

    weights = vec(transpose(cleanup) * Float64.(packet.weights))
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

function _pqs_source_box_route_driver_axis_representation_for_gto(axis_bundles, axis::Symbol)
    pgdg = _nested_axis_pgdg(axis_bundles, axis)
    hasproperty(pgdg, :auxiliary_layer) ||
        throw(ArgumentError("PQS/GTO sidecar requires PGDG auxiliary layers"))
    auxiliary_layer = pgdg.auxiliary_layer
    if auxiliary_layer isa MappedPGDGLocalized1D
        metadata = BasisMetadata1D(
            :mapped_pgdg_localized,
            _basis_family_name(auxiliary_layer),
            mapping(auxiliary_layer),
            Float64[Float64(value) for value in centers(auxiliary_layer)],
            Float64[Float64(value) for value in reference_centers(auxiliary_layer)],
            Float64[Float64(value) for value in integral_weights(auxiliary_layer)],
            _basis_label_vector(length(auxiliary_layer)),
            primitive_set(auxiliary_layer),
            Matrix{Float64}(stencil_matrix(auxiliary_layer)),
        )
        return _basis_representation_from_metadata(metadata; operators = (:overlap,))
    end
    return basis_representation(auxiliary_layer; operators = (:overlap,))
end

function _pqs_source_box_route_driver_pqs_gto_sidecar(inputs)
    final_basis = inputs.final_basis
    source_plan = inputs.source_plan
    supplement = inputs.supplement_representation
    support_states =
        _pqs_source_box_route_driver_physical_gausslet_support_states(source_plan)
    axis_representations = (;
        x = _pqs_source_box_route_driver_axis_representation_for_gto(
            source_plan.axis_bundles,
            :x,
        ),
        y = _pqs_source_box_route_driver_axis_representation_for_gto(
            source_plan.axis_bundles,
            :y,
        ),
        z = _pqs_source_box_route_driver_axis_representation_for_gto(
            source_plan.axis_bundles,
            :z,
        ),
    )
    support_gto_cross =
        zeros(Float64, length(support_states), length(supplement.orbitals))
    for (column, orbital) in pairs(supplement.orbitals)
        overlap_x =
            _cartesian_basis_supplement_axis_primitive_cross(
                axis_representations.x,
                orbital,
                :x,
            )
        overlap_y =
            _cartesian_basis_supplement_axis_primitive_cross(
                axis_representations.y,
                orbital,
                :y,
            )
        overlap_z =
            _cartesian_basis_supplement_axis_primitive_cross(
                axis_representations.z,
                orbital,
                :z,
            )
        for (row, (ix, iy, iz)) in pairs(support_states)
            value = 0.0
            for primitive in eachindex(orbital.coefficients)
                value +=
                    Float64(orbital.coefficients[primitive]) *
                    overlap_x[ix, primitive] *
                    overlap_y[iy, primitive] *
                    overlap_z[iz, primitive]
            end
            support_gto_cross[row, column] = value
        end
    end
    final_coefficients = Matrix{Float64}(final_basis.final_coefficients)
    final_gto_cross_overlap =
        Matrix{Float64}(transpose(final_coefficients) * support_gto_cross)
    gto_self_overlap =
        Matrix{Float64}(_cartesian_supplement_cross_overlap(supplement, supplement))
    gto_residual_overlap =
        Matrix{Float64}(
            gto_self_overlap -
            transpose(final_gto_cross_overlap) * final_gto_cross_overlap,
        )
    residual_symmetry_error =
        norm(gto_residual_overlap - transpose(gto_residual_overlap), Inf)
    residual_eigenvalues =
        isempty(gto_residual_overlap) ?
        Float64[] :
        eigvals(Symmetric(0.5 .* (gto_residual_overlap + transpose(gto_residual_overlap))))
    diagnostics = (;
        sidecar_kind = :pqs_h2_residual_gto_sidecar,
        final_dimension = size(final_gto_cross_overlap, 1),
        gto_dimension = size(final_gto_cross_overlap, 2),
        support_dimension = length(support_states),
        final_gto_cross_overlap_finite = all(isfinite, final_gto_cross_overlap),
        gto_self_overlap_finite = all(isfinite, gto_self_overlap),
        gto_residual_overlap_finite = all(isfinite, gto_residual_overlap),
        gto_residual_overlap_symmetry_error = residual_symmetry_error,
        gto_residual_overlap_eigenvalue_min =
            isempty(residual_eigenvalues) ? nothing : minimum(residual_eigenvalues),
        gto_residual_overlap_eigenvalue_max =
            isempty(residual_eigenvalues) ? nothing : maximum(residual_eigenvalues),
        donor_kernels = (
            :basis_representation,
            :_cartesian_basis_supplement_axis_primitive_cross,
            :_cartesian_supplement_cross_overlap,
        ),
        provider_blocks_included = false,
    )
    return (;
        final_gto_cross_overlap,
        support_gto_cross,
        support_states,
        axis_representations,
        gto_self_overlap,
        gto_residual_overlap,
        diagnostics,
    )
end

function _pqs_source_box_route_driver_pqs_gto_residual_transform(
    sidecar;
    residual_overlap_cutoff::Real = 1.0e-10,
)
    residual_overlap = Matrix{Float64}(sidecar.gto_residual_overlap)
    residual_overlap_sym =
        Matrix{Float64}(0.5 .* (residual_overlap .+ transpose(residual_overlap)))
    decomposition = eigen(Symmetric(residual_overlap_sym))
    eigenvalues = Float64[decomposition.values...]
    cutoff = Float64(residual_overlap_cutoff)
    any(<(-cutoff), eigenvalues) && throw(
        ArgumentError(
            "H2 PQS residual-GTO overlap has eigenvalues below -$(cutoff)",
        ),
    )
    keep = findall(>(cutoff), eigenvalues)
    residual_rank = length(keep)
    residual_transform =
        isempty(keep) ?
        zeros(Float64, size(residual_overlap, 1), 0) :
        Matrix{Float64}(
            decomposition.vectors[:, keep] *
            Diagonal(1.0 ./ sqrt.(eigenvalues[keep])),
        )
    residual_identity =
        transpose(residual_transform) *
        residual_overlap_sym *
        residual_transform
    residual_overlap_identity_error =
        residual_rank == 0 ?
        0.0 :
        norm(residual_identity - Matrix{Float64}(I, residual_rank, residual_rank), Inf)
    return (;
        residual_transform,
        residual_rank,
        residual_overlap_eigenvalues = eigenvalues,
        residual_overlap_identity_error,
        residual_overlap_eigenvalue_min =
            isempty(eigenvalues) ? nothing : minimum(eigenvalues),
        residual_overlap_eigenvalue_max =
            isempty(eigenvalues) ? nothing : maximum(eigenvalues),
        residual_overlap_cutoff = cutoff,
    )
end

function _pqs_source_box_route_driver_pqs_h2_residual_gto_provider_packet(
    inputs,
    sidecar,
    residual,
)
    expansion = coulomb_gaussian_expansion(doacc = false)
    route_metadata = inputs.route_metadata
    gto_primitive_arrays = map(
        _pqs_source_box_route_driver_gto_primitive_arrays,
        inputs.supplement_representation.orbitals,
    )
    return (;
        final_basis = inputs.final_basis,
        source_plan = inputs.source_plan,
        support_states = sidecar.support_states,
        axis_representations = sidecar.axis_representations,
        supplement_representation = inputs.supplement_representation,
        gto_primitive_arrays,
        final_coefficients = Matrix{Float64}(inputs.final_basis.final_coefficients),
        support_gto_cross = Matrix{Float64}(sidecar.support_gto_cross),
        h_ff = Matrix{Float64}(inputs.h1_hamiltonian.hamiltonian_matrix),
        density_interaction = inputs.density_interaction,
        s_fg = Matrix{Float64}(sidecar.final_gto_cross_overlap),
        s_gg = Matrix{Float64}(sidecar.gto_self_overlap),
        s_r = Matrix{Float64}(sidecar.gto_residual_overlap),
        residual_transform = Matrix{Float64}(residual.residual_transform),
        residual_rank = residual.residual_rank,
        route_metadata,
        nuclear_charges = Float64.(route_metadata.nuclear_charges),
        atom_locations = route_metadata.atom_locations,
        coulomb_coefficients = Float64.(expansion.coefficients),
        coulomb_exponents = Float64.(expansion.exponents),
    )
end

function _pqs_source_box_route_driver_pqs_h2_residual_gto_density_descriptor(
    inputs,
    sidecar,
    residual,
)
    final_basis = inputs.final_basis
    density_interaction = inputs.density_interaction
    cleanup = Matrix{Float64}(final_basis.combined_lowdin_cleanup)
    s_fg = Matrix{Float64}(sidecar.final_gto_cross_overlap)
    residual_transform = Matrix{Float64}(residual.residual_transform)
    p_dimension, f_dimension = size(cleanup)
    size(s_fg, 1) == f_dimension ||
        throw(DimensionMismatch("density descriptor S_FG row count must match final dimension"))
    g_dimension = size(s_fg, 2)
    size(residual_transform, 1) == g_dimension ||
        throw(DimensionMismatch("density descriptor residual transform row count must match GTO dimension"))
    residual_rank = size(residual_transform, 2)
    residual_rank == residual.residual_rank ||
        throw(DimensionMismatch("density descriptor residual rank mismatch"))
    p_projection_of_g = Matrix{Float64}(cleanup * s_fg)
    residual_carrier =
        Matrix{Float64}(
            vcat(
                -p_projection_of_g,
                Matrix{Float64}(I, g_dimension, g_dimension),
            ) * residual_transform,
        )
    pair_matrix = Matrix{Float64}(density_interaction.pre_final_pair_matrix)
    size(pair_matrix) == (p_dimension, p_dimension) ||
        throw(DimensionMismatch("density descriptor P-P pair matrix dimension mismatch"))
    get(density_interaction, :density_gauge, nothing) ===
        :pre_final_localized_positive_weight ||
        throw(ArgumentError("density descriptor requires the pre-final localized positive-weight gauge"))
    return (;
        density_gauge = :pre_final_localized_positive_weight,
        augmented_density_space = (:pre_final_pqs, :residual_gto),
        p_dimension,
        f_dimension,
        g_dimension,
        residual_rank,
        p_projection_of_g,
        p_projection_of_g_size = size(p_projection_of_g),
        residual_orbital_coefficients_in_density_carrier = residual_carrier,
        residual_orbital_coefficients_in_density_carrier_size =
            size(residual_carrier),
    )
end

function _pqs_source_box_route_driver_pqs_gto_support_moment_matrix(
    packet,
    axis_index::Int,
    term::Symbol,
)
    states = packet.support_states
    axes = packet.axis_representations
    result = zeros(Float64, length(states), length(packet.gto_primitive_arrays))
    scratch = zeros(Float64, length(states))
    for (column, orbital_arrays) in pairs(packet.gto_primitive_arrays)
        overlap_tables =
            _pqs_source_box_route_driver_gto_axis_cross_tables(
                axes,
                orbital_arrays,
                :overlap,
            )
        moment_tables =
            _pqs_source_box_route_driver_gto_axis_cross_tables(
                axes,
                orbital_arrays,
                term,
            )
        axis_tables = ntuple(
            axis -> axis == axis_index ? moment_tables[axis] : overlap_tables[axis],
            3,
        )
        _pqs_source_box_route_driver_gto_support_column!(
            scratch,
            states,
            orbital_arrays.coefficients,
            axis_tables,
        )
        result[:, column] .= scratch
    end
    return result
end

function _pqs_source_box_route_driver_pqs_gto_self_moment_matrix(
    packet,
    axis_index::Int,
    term::Symbol,
)
    norbital = length(packet.gto_primitive_arrays)
    result = zeros(Float64, norbital, norbital)
    for column in eachindex(packet.gto_primitive_arrays)
        right_arrays = packet.gto_primitive_arrays[column]
        for row in eachindex(packet.gto_primitive_arrays)
            left_arrays = packet.gto_primitive_arrays[row]
            overlap_tables =
                _pqs_source_box_route_driver_gto_axis_self_tables(
                    left_arrays,
                    right_arrays,
                    :overlap,
                )
            moment_tables =
                _pqs_source_box_route_driver_gto_axis_self_tables(
                    left_arrays,
                    right_arrays,
                    term,
                )
            axis_tables = ntuple(
                axis -> axis == axis_index ? moment_tables[axis] : overlap_tables[axis],
                3,
            )
            result[row, column] =
                _pqs_source_box_route_driver_gto_weighted_hadamard(
                    left_arrays.coefficients,
                    axis_tables[1],
                    axis_tables[2],
                    axis_tables[3],
                    right_arrays.coefficients,
                )
        end
    end
    return Matrix{Float64}(0.5 .* (result .+ transpose(result)))
end

function _pqs_source_box_route_driver_pqs_support_moment_matrix(
    packet,
    axis_index::Int,
    term::Symbol,
)
    axis_index in 1:3 || throw(ArgumentError("axis index must be 1, 2, or 3"))
    states = packet.support_states
    bundles = packet.source_plan.axis_bundles
    pgdgs = ntuple(axis -> _nested_axis_pgdg(
        bundles,
        _pqs_source_box_route_driver_axis_symbol(axis),
    ), 3)
    operators = ntuple(3) do axis
        if axis == axis_index
            term === :position && return pgdgs[axis].position
            term === :x2 && return pgdgs[axis].x2
            throw(ArgumentError("support moment term must be :position or :x2"))
        end
        return pgdgs[axis].overlap
    end
    return _pqs_multilayer_support_product_matrix(
        states,
        states,
        operators[1],
        operators[2],
        operators[3],
    )
end

function _pqs_source_box_route_driver_pqs_support_overlap_matrix(packet)
    states = packet.support_states
    bundles = packet.source_plan.axis_bundles
    pgdgs = ntuple(axis -> _nested_axis_pgdg(
        bundles,
        _pqs_source_box_route_driver_axis_symbol(axis),
    ), 3)
    return _pqs_multilayer_support_product_matrix(
        states,
        states,
        pgdgs[1].overlap,
        pgdgs[2].overlap,
        pgdgs[3].overlap,
    )
end

function _pqs_source_box_route_driver_pqs_support_weights(packet)
    states = packet.support_states
    bundles = packet.source_plan.axis_bundles
    pgdgs = ntuple(axis -> _nested_axis_pgdg(
        bundles,
        _pqs_source_box_route_driver_axis_symbol(axis),
    ), 3)
    weights = Vector{Float64}(undef, length(states))
    @inbounds for (index, (ix, iy, iz)) in pairs(states)
        weights[index] =
            Float64(pgdgs[1].weights[ix]) *
            Float64(pgdgs[2].weights[iy]) *
            Float64(pgdgs[3].weights[iz])
    end
    all(isfinite, weights) ||
        throw(ArgumentError("support weights for residual-GTO density provider contain non-finite entries"))
    all(>(0.0), weights) ||
        throw(ArgumentError("support weights for residual-GTO density provider must be positive"))
    return weights
end

function _pqs_source_box_route_driver_density_carrier_moment_matrix(
    packet,
    density_descriptor,
    axis_index::Int,
    term::Symbol,
)
    final_basis = packet.final_basis
    pre_final_coefficients = Matrix{Float64}(final_basis.pre_final_coefficients)
    support_pp =
        _pqs_source_box_route_driver_pqs_support_moment_matrix(
            packet,
            axis_index,
            term,
        )
    support_pg =
        _pqs_source_box_route_driver_pqs_gto_support_moment_matrix(
            packet,
            axis_index,
            term,
        )
    pp = Matrix{Float64}(transpose(pre_final_coefficients) * support_pp * pre_final_coefficients)
    pg = Matrix{Float64}(transpose(pre_final_coefficients) * support_pg)
    gg =
        _pqs_source_box_route_driver_pqs_gto_self_moment_matrix(
            packet,
            axis_index,
            term,
        )
    size(pp) == (density_descriptor.p_dimension, density_descriptor.p_dimension) ||
        throw(DimensionMismatch("density carrier P-P moment dimension mismatch"))
    size(pg) == (density_descriptor.p_dimension, density_descriptor.g_dimension) ||
        throw(DimensionMismatch("density carrier P-G moment dimension mismatch"))
    size(gg) == (density_descriptor.g_dimension, density_descriptor.g_dimension) ||
        throw(DimensionMismatch("density carrier G-G moment dimension mismatch"))
    return [pp pg; transpose(pg) gg]
end

function _pqs_source_box_route_driver_pqs_h2_residual_gto_moments(
    packet,
    density_descriptor,
)
    carrier = Matrix{Float64}(
        density_descriptor.residual_orbital_coefficients_in_density_carrier,
    )
    pre_final_coefficients =
        Matrix{Float64}(packet.final_basis.pre_final_coefficients)
    support_overlap =
        _pqs_source_box_route_driver_pqs_support_overlap_matrix(packet)
    pp_overlap =
        Matrix{Float64}(
            transpose(pre_final_coefficients) *
            support_overlap *
            pre_final_coefficients,
        )
    pg_overlap =
        Matrix{Float64}(transpose(pre_final_coefficients) * packet.support_gto_cross)
    raw_overlap = [
        pp_overlap pg_overlap
        transpose(pg_overlap) packet.s_gg
    ]
    residual_overlap = Matrix{Float64}(transpose(carrier) * raw_overlap * carrier)
    residual_overlap_error =
        norm(
            residual_overlap -
            Matrix{Float64}(I, size(residual_overlap, 1), size(residual_overlap, 2)),
            Inf,
        )
    residual_overlap_error <= 1.0e-8 ||
        throw(ArgumentError("residual density moments require an orthonormal residual carrier"))
    nresidual = size(carrier, 2)
    centers = zeros(Float64, nresidual, 3)
    widths = zeros(Float64, nresidual, 3)
    norms = zeros(Float64, nresidual)
    position_matrices = ntuple(
        axis_index ->
            _pqs_source_box_route_driver_density_carrier_moment_matrix(
                packet,
                density_descriptor,
                axis_index,
                :position,
            ),
        3,
    )
    second_moment_matrices = ntuple(
        axis_index ->
            _pqs_source_box_route_driver_density_carrier_moment_matrix(
                packet,
                density_descriptor,
                axis_index,
                :x2,
            ),
        3,
    )
    for residual in 1:nresidual
        vector = view(carrier, :, residual)
        norm_value = Float64(dot(vector, raw_overlap * vector))
        norm_value > 1.0e-12 ||
            throw(ArgumentError("residual density moment extraction requires nonzero residual norm"))
        norms[residual] = norm_value
        for axis_index in 1:3
            position = position_matrices[axis_index]
            second_moment = second_moment_matrices[axis_index]
            center = Float64(dot(vector, position * vector) / norm_value)
            x2 = Float64(dot(vector, second_moment * vector) / norm_value)
            variance = x2 - center^2
            variance > 1.0e-12 ||
                throw(ArgumentError("MWG residual density moment extraction requires positive variances"))
            centers[residual, axis_index] = center
            widths[residual, axis_index] = sqrt(2.0 * variance)
        end
    end
    return (;
        residual_centers = centers,
        residual_widths = widths,
        residual_norms = norms,
        residual_overlap_error,
        residual_center_min = minimum(centers),
        residual_center_max = maximum(centers),
        residual_width_min = minimum(widths),
        residual_width_max = maximum(widths),
        residual_widths_positive = all(>(0.0), widths),
    )
end

function _pqs_source_box_route_driver_pqs_h2_residual_gto_density_blocks(
    packet,
    density_descriptor,
)
    moments =
        _pqs_source_box_route_driver_pqs_h2_residual_gto_moments(
            packet,
            density_descriptor,
        )
    bundles = packet.source_plan.axis_bundles
    components = _qwrg_mwg_interaction_components(
        _nested_axis_bundle(bundles, :x),
        _nested_axis_bundle(bundles, :y),
        _nested_axis_bundle(bundles, :z),
        coulomb_gaussian_expansion(doacc = false),
        moments.residual_centers,
        moments.residual_widths,
    )
    dims = (
        size(_nested_axis_pgdg(bundles, :x).overlap, 1),
        size(_nested_axis_pgdg(bundles, :y).overlap, 1),
        size(_nested_axis_pgdg(bundles, :z).overlap, 1),
    )
    support_states = packet.support_states
    support_rows =
        Int[_cartesian_flat_index(ix, iy, iz, dims) for (ix, iy, iz) in support_states]
    support_residual_density_normalized =
        Matrix{Float64}(components.gausslet_residual[support_rows, :])
    support_weights =
        _pqs_source_box_route_driver_pqs_support_weights(packet)
    support_residual_raw_numerator =
        support_residual_density_normalized .* reshape(support_weights, :, 1)
    density_interaction = packet.density_interaction
    pre_final_coefficients = Matrix{Float64}(packet.final_basis.pre_final_coefficients)
    pre_final_weights = Float64.(density_interaction.pre_final_weights)
    weighted_coefficients =
        pre_final_coefficients .* reshape(1.0 ./ pre_final_weights, 1, :)
    v_pr =
        Matrix{Float64}(
            transpose(weighted_coefficients) * support_residual_raw_numerator,
        )
    v_rr = Matrix{Float64}(components.residual_residual)
    v_rr = Matrix{Float64}(0.5 .* (v_rr .+ transpose(v_rr)))
    v_pp = Matrix{Float64}(density_interaction.pre_final_pair_matrix)
    v_dd = [
        v_pp v_pr
        transpose(v_pr) v_rr
    ]
    v_dd_symmetry_error = norm(v_dd - transpose(v_dd), Inf)
    return (;
        augmented_density_gauge = :pre_final_localized_positive_weight,
        augmented_density_space = (:pre_final_pqs, :residual_gto),
        v_pp_pair_matrix = v_pp,
        v_pr_pair_matrix = v_pr,
        v_rr_pair_matrix = v_rr,
        augmented_pair_matrix = v_dd,
        augmented_density_dimension = size(v_dd, 1),
        v_pr_pair_matrix_size = size(v_pr),
        v_rr_pair_matrix_size = size(v_rr),
        augmented_pair_matrix_size = size(v_dd),
        augmented_pair_matrix_finite = all(isfinite, v_dd),
        augmented_pair_matrix_symmetry_error = v_dd_symmetry_error,
        residual_centers = moments.residual_centers,
        residual_widths = moments.residual_widths,
        residual_overlap_error = moments.residual_overlap_error,
        residual_width_min = moments.residual_width_min,
        residual_width_max = moments.residual_width_max,
        residual_widths_positive = moments.residual_widths_positive,
        support_weight_min = minimum(support_weights),
        support_weight_max = maximum(support_weights),
    )
end

function _pqs_source_box_route_driver_restricted_one_orbital_self_coulomb(
    pair_matrix,
    orbital_coefficients,
)
    matrix = Matrix{Float64}(pair_matrix)
    vector = Float64[Float64(value) for value in orbital_coefficients]
    size(matrix) == (length(vector), length(vector)) ||
        throw(DimensionMismatch("one-orbital pair matrix and coefficient dimensions mismatch"))
    all(isfinite, matrix) ||
        throw(ArgumentError("one-orbital pair matrix contains non-finite entries"))
    all(isfinite, vector) ||
        throw(ArgumentError("one-orbital coefficients contain non-finite entries"))
    v = 0.5 .* (matrix .+ transpose(matrix))
    density = vector * transpose(vector)
    rho = 0.5 .* (density .+ transpose(density))
    occupations = vec(diag(rho))
    direct = 2.0 * dot(occupations, v * occupations)
    exchange = dot(vec(rho), vec(v .* rho))
    return direct - exchange
end

function _pqs_source_box_route_driver_pqs_h2_augmented_h1_j_diagnostic(
    one_body_blocks,
    density_blocks,
    density_interaction,
)
    isnothing(one_body_blocks) &&
        throw(ArgumentError("augmented H1-J diagnostic requires one-body provider blocks"))
    isnothing(density_blocks) &&
        throw(ArgumentError("augmented H1-J diagnostic requires density provider blocks"))
    augmented_h1 =
        Matrix{Float64}(one_body_blocks.augmented_one_body_hamiltonian)
    augmented_pair = Matrix{Float64}(density_blocks.augmented_pair_matrix)
    augmented_dimension = one_body_blocks.augmented_dimension
    density_dimension = density_blocks.augmented_density_dimension
    size(augmented_h1) == (augmented_dimension, augmented_dimension) ||
        throw(DimensionMismatch("augmented H1 matrix dimension mismatch"))
    size(augmented_pair) == (density_dimension, density_dimension) ||
        throw(DimensionMismatch("augmented pair matrix dimension mismatch"))
    final_to_pre_final =
        Matrix{Float64}(density_interaction.final_to_pre_final_coefficients)
    p_dimension, final_dimension = size(final_to_pre_final)
    residual_rank = density_dimension - p_dimension
    augmented_dimension == final_dimension + residual_rank ||
        throw(DimensionMismatch("augmented H1 and density gauge dimensions mismatch"))
    eigensystem = eigen(Symmetric(0.5 .* (augmented_h1 .+ transpose(augmented_h1))))
    orbital = Vector{Float64}(eigensystem.vectors[:, 1])
    f_coefficients = @view orbital[1:final_dimension]
    r_coefficients = @view orbital[(final_dimension + 1):augmented_dimension]
    density_coefficients =
        vcat(final_to_pre_final * Vector{Float64}(f_coefficients),
             Vector{Float64}(r_coefficients))
    self_coulomb =
        _pqs_source_box_route_driver_restricted_one_orbital_self_coulomb(
            augmented_pair,
            density_coefficients,
        )
    isfinite(self_coulomb) && self_coulomb > 0 ||
        throw(ArgumentError("augmented H1-J self Coulomb diagnostic must be finite and positive"))
    return (;
        density_gauge = :pre_final_pqs_plus_residual_gto,
        final_coefficients_length = final_dimension,
        residual_coefficients_length = residual_rank,
        density_coefficients_length = length(density_coefficients),
        augmented_h1_lowest = Float64(eigensystem.values[1]),
        self_coulomb,
        h1_j_self_coulomb = self_coulomb,
        finite = isfinite(self_coulomb),
        positive = self_coulomb > 0,
    )
end

function _pqs_source_box_route_driver_nuclear_repulsion(route_metadata)
    charges = Float64.(route_metadata.nuclear_charges)
    locations = route_metadata.atom_locations
    length(charges) == length(locations) ||
        throw(DimensionMismatch("nuclear charge and location counts differ"))
    repulsion = 0.0
    for right in 2:length(charges), left in 1:(right - 1)
        distance = norm(Float64.(locations[right] .- locations[left]))
        distance > 0.0 ||
            throw(ArgumentError("nuclear repulsion requires distinct centers"))
        repulsion += charges[left] * charges[right] / distance
    end
    return repulsion
end

function _pqs_source_box_route_driver_augmented_density_transform(
    density_interaction,
    one_body_blocks,
    density_blocks,
)
    final_to_pre_final =
        Matrix{Float64}(density_interaction.final_to_pre_final_coefficients)
    p_dimension, final_dimension = size(final_to_pre_final)
    augmented_dimension = one_body_blocks.augmented_dimension
    density_dimension = density_blocks.augmented_density_dimension
    residual_rank = density_dimension - p_dimension
    augmented_dimension == final_dimension + residual_rank ||
        throw(DimensionMismatch("augmented RHF density transform dimension mismatch"))
    transform = zeros(Float64, density_dimension, augmented_dimension)
    transform[1:p_dimension, 1:final_dimension] .= final_to_pre_final
    transform[(p_dimension + 1):density_dimension,
        (final_dimension + 1):augmented_dimension] .= I(residual_rank)
    return transform
end

function _pqs_source_box_route_driver_augmented_rhf_fock_components(
    hamiltonian,
    pair_matrix,
    density_transform,
    spin_summed_density,
)
    h = Matrix{Float64}(hamiltonian)
    v = Matrix{Float64}(pair_matrix)
    transform = Matrix{Float64}(density_transform)
    density = Matrix{Float64}(spin_summed_density)
    size(h, 1) == size(h, 2) ||
        throw(DimensionMismatch("augmented RHF Hamiltonian must be square"))
    size(density) == size(h) ||
        throw(DimensionMismatch("augmented RHF density dimension mismatch"))
    size(v, 1) == size(v, 2) ||
        throw(DimensionMismatch("augmented RHF pair matrix must be square"))
    size(transform, 2) == size(h, 1) ||
        throw(DimensionMismatch("augmented RHF density transform column mismatch"))
    size(transform, 1) == size(v, 1) ||
        throw(DimensionMismatch("augmented RHF density transform row mismatch"))
    all(isfinite, h) && all(isfinite, v) && all(isfinite, transform) &&
        all(isfinite, density) ||
        throw(ArgumentError("augmented RHF inputs contain non-finite entries"))

    orbital_density = 0.5 .* density
    density_gauge =
        transform * orbital_density * transpose(transform)
    density_gauge =
        0.5 .* (density_gauge .+ transpose(density_gauge))
    interaction = 0.5 .* (v .+ transpose(v))
    occupations = vec(diag(density_gauge))
    density_fock =
        2.0 .* Diagonal(interaction * occupations) .-
        density_gauge .* interaction
    coulomb_exchange =
        transpose(transform) * Matrix(density_fock) * transform
    fock = h + coulomb_exchange
    fock = 0.5 .* (fock .+ transpose(fock))
    one_body_energy = tr(density * h)
    two_body_energy = 0.5 * tr(density * coulomb_exchange)
    return (;
        fock,
        coulomb_exchange,
        one_body_energy,
        two_body_energy,
        electronic_energy = one_body_energy + two_body_energy,
        density_gauge,
    )
end

function _pqs_source_box_route_driver_closed_shell_density(orbital)
    vector = Vector{Float64}(orbital)
    all(isfinite, vector) ||
        throw(ArgumentError("closed-shell RHF orbital contains non-finite entries"))
    return 2.0 .* (vector * transpose(vector))
end

function _pqs_source_box_route_driver_pqs_h2_private_augmented_rhf_smoke(
    one_body_blocks,
    density_blocks,
    density_interaction,
    route_metadata;
    max_iterations::Int = 25,
    energy_atol::Real = 1.0e-10,
    residual_atol::Real = 1.0e-8,
    density_atol::Real = 1.0e-8,
)
    isnothing(one_body_blocks) &&
        throw(ArgumentError("private augmented RHF requires one-body provider blocks"))
    isnothing(density_blocks) &&
        throw(ArgumentError("private augmented RHF requires density provider blocks"))
    max_iterations > 0 ||
        throw(ArgumentError("private augmented RHF requires a positive iteration limit"))
    h = Matrix{Float64}(one_body_blocks.augmented_one_body_hamiltonian)
    pair_matrix = Matrix{Float64}(density_blocks.augmented_pair_matrix)
    transform =
        _pqs_source_box_route_driver_augmented_density_transform(
            density_interaction,
            one_body_blocks,
            density_blocks,
        )
    initial_eigensystem = eigen(Symmetric(0.5 .* (h .+ transpose(h))))
    density =
        _pqs_source_box_route_driver_closed_shell_density(
            initial_eigensystem.vectors[:, 1],
        )
    previous_energy = Inf
    converged = false
    iteration_count = 0
    final_energy_delta = Inf
    final_density_delta = Inf
    for iteration in 1:max_iterations
        components =
            _pqs_source_box_route_driver_augmented_rhf_fock_components(
                h,
                pair_matrix,
                transform,
                density,
            )
        fock_eigensystem = eigen(Symmetric(components.fock))
        updated_density =
            _pqs_source_box_route_driver_closed_shell_density(
                fock_eigensystem.vectors[:, 1],
            )
        mixed_density =
            iteration < 4 ?
            0.35 .* updated_density .+ 0.65 .* density :
            updated_density
        energy_delta = abs(components.electronic_energy - previous_energy)
        density_delta = norm(mixed_density - density, Inf)
        commutator = norm(components.fock * density - density * components.fock, Inf)
        previous_energy = components.electronic_energy
        density = mixed_density
        iteration_count = iteration
        final_energy_delta = energy_delta
        final_density_delta = density_delta
        if energy_delta <= Float64(energy_atol) &&
           density_delta <= Float64(density_atol) &&
           commutator <= Float64(residual_atol)
            density = updated_density
            converged = true
            break
        end
    end
    final_components =
        _pqs_source_box_route_driver_augmented_rhf_fock_components(
            h,
            pair_matrix,
            transform,
            density,
        )
    density_trace = tr(density)
    idempotency_error = norm(density * density - 2.0 .* density, Inf)
    commutator_residual =
        norm(final_components.fock * density - density * final_components.fock, Inf)
    trace_error = abs(density_trace - 2.0)
    converged = converged ||
        (
            final_energy_delta <= Float64(energy_atol) &&
            final_density_delta <= Float64(density_atol) &&
            commutator_residual <= Float64(residual_atol)
        )
    converged ||
        throw(ArgumentError("private augmented RHF smoke did not converge"))
    trace_error <= Float64(density_atol) ||
        throw(ArgumentError("private augmented RHF density trace mismatch"))
    idempotency_error <= 10.0 * Float64(density_atol) ||
        throw(ArgumentError("private augmented RHF density idempotency mismatch"))
    nuclear_repulsion =
        _pqs_source_box_route_driver_nuclear_repulsion(route_metadata)
    return (;
        rhf_kind = :private_augmented_residual_gto_rhf_smoke,
        converged,
        iteration_count,
        density_trace,
        idempotency_error,
        commutator_residual,
        one_body_energy = final_components.one_body_energy,
        two_body_energy = final_components.two_body_energy,
        electronic_energy = final_components.electronic_energy,
        nuclear_repulsion,
        total_with_nuclear_repulsion =
            final_components.electronic_energy + nuclear_repulsion,
    )
end

function _pqs_source_box_route_driver_pqs_h2_residual_gto_ham_handoff(
    one_body_blocks,
    density_blocks,
    density_interaction,
    augmented_h1_j,
    private_augmented_rhf,
)
    isnothing(one_body_blocks) &&
        throw(ArgumentError("H2 residual-GTO Ham handoff requires one-body blocks"))
    isnothing(density_blocks) &&
        throw(ArgumentError("H2 residual-GTO Ham handoff requires density blocks"))
    isnothing(private_augmented_rhf) &&
        throw(ArgumentError("H2 residual-GTO Ham handoff requires private RHF facts"))
    one_body = Matrix{Float64}(one_body_blocks.augmented_one_body_hamiltonian)
    pair_matrix = Matrix{Float64}(density_blocks.augmented_pair_matrix)
    orbital_to_density =
        _pqs_source_box_route_driver_augmented_density_transform(
            density_interaction,
            one_body_blocks,
            density_blocks,
        )
    orbital_dimension = size(one_body, 1)
    density_dimension = size(pair_matrix, 1)
    size(one_body, 2) == orbital_dimension ||
        throw(DimensionMismatch("Ham handoff one-body matrix must be square"))
    size(pair_matrix, 2) == density_dimension ||
        throw(DimensionMismatch("Ham handoff density pair matrix must be square"))
    size(orbital_to_density) == (density_dimension, orbital_dimension) ||
        throw(DimensionMismatch("Ham handoff orbital-to-density transform shape mismatch"))
    all(isfinite, one_body) && all(isfinite, pair_matrix) &&
        all(isfinite, orbital_to_density) ||
        throw(ArgumentError("Ham handoff matrices contain non-finite entries"))
    one_body_symmetry_error = norm(one_body - transpose(one_body), Inf)
    pair_symmetry_error = norm(pair_matrix - transpose(pair_matrix), Inf)
    return (;
        handoff_kind = :pqs_h2_residual_gto_ham_handoff,
        visibility = :private_experimental,
        model = :density_density,
        orbital_basis = (:final_pqs, :residual_gto),
        density_basis = (:pre_final_pqs, :residual_gto),
        one_body_hamiltonian = one_body,
        density_pair_matrix = pair_matrix,
        orbital_to_density,
        electron_count = 2,
        spin_sectors = (; nup = 1, ndn = 1),
        nuclear_repulsion = private_augmented_rhf.nuclear_repulsion,
        diagnostics = (;
            orbital_dimension,
            density_dimension,
            residual_rank = density_blocks.augmented_density_dimension -
                size(density_interaction.final_to_pre_final_coefficients, 1),
            one_body_symmetry_error,
            pair_symmetry_error,
            h1_lowest = one_body_blocks.augmented_h1_lowest,
            h1_j_self_coulomb = augmented_h1_j.self_coulomb,
            private_rhf_converged = private_augmented_rhf.converged,
            private_rhf_iterations = private_augmented_rhf.iteration_count,
            private_rhf_total_with_nuclear_repulsion =
                private_augmented_rhf.total_with_nuclear_repulsion,
            private_rhf_commutator_residual =
                private_augmented_rhf.commutator_residual,
        ),
    )
end

@inline function _pqs_source_box_route_driver_axis_index(axis::Symbol)
    axis === :x && return 1
    axis === :y && return 2
    axis === :z && return 3
    throw(ArgumentError("axis must be :x, :y, or :z"))
end

@inline function _pqs_source_box_route_driver_axis_symbol(index::Int)
    1 <= index <= 3 || throw(ArgumentError("axis index must be 1, 2, or 3"))
    return (:x, :y, :z)[index]
end

function _pqs_source_box_route_driver_gto_primitive_arrays(orbital)
    if hasproperty(orbital, :primitive_normalization)
        orbital.primitive_normalization == :axiswise_normalized_cartesian_gaussian ||
            throw(ArgumentError("unsupported GTO primitive normalization"))
    end
    length(orbital.coefficients) == length(orbital.exponents) ||
        throw(ArgumentError("GTO coefficients and exponents must have the same length"))
    length(orbital.center) == 3 ||
        throw(ArgumentError("GTO orbital center must have exactly three coordinates"))
    length(orbital.angular_powers) == 3 ||
        throw(ArgumentError("GTO orbital angular_powers must have exactly three entries"))
    # Axis-last layout keeps the axis-wise primitive views contiguous:
    # centers[:, axis_index] and angular_powers[:, axis_index].
    nprimitive = length(orbital.exponents)
    centers = zeros(Float64, nprimitive, 3)
    angular_powers = zeros(Int, nprimitive, 3)
    for primitive in 1:nprimitive, axis_index in 1:3
        centers[primitive, axis_index] = Float64(orbital.center[axis_index])
        angular_powers[primitive, axis_index] =
            Int(orbital.angular_powers[axis_index])
    end
    return (;
        exponents = Float64.(orbital.exponents),
        coefficients = Float64.(orbital.coefficients),
        centers,
        angular_powers,
    )
end

function _pqs_source_box_route_driver_gto_axis_cross(
    basis::BasisRepresentation1D,
    orbital_arrays,
    axis::Symbol,
    term::Symbol;
    factor_center::Union{Nothing,Real} = nothing,
    factor_exponent::Union{Nothing,Real} = nothing,
)
    basis_primitives = collect(primitives(primitive_set(basis)))
    all(primitive -> primitive isa Gaussian, basis_primitives) || throw(
        ArgumentError(
            "PQS/GTO one-body provider blocks require Gaussian 1D primitives on the Cartesian axis representation",
        ),
    )
    axis_index = _pqs_source_box_route_driver_axis_index(axis)
    centers_axis = @view orbital_arrays.centers[:, axis_index]
    powers_axis = @view orbital_arrays.angular_powers[:, axis_index]
    primitive_cross = zeros(
        Float64,
        length(basis_primitives),
        length(orbital_arrays.exponents),
    )
    for (row, primitive) in pairs(basis_primitives)
        alpha_basis = _qwrg_gaussian_exponent(primitive::Gaussian)
        for column in eachindex(orbital_arrays.exponents)
            exponent = orbital_arrays.exponents[column]
            center_value = centers_axis[column]
            power = powers_axis[column]
            prefactor = _qwrg_atomic_shell_prefactor(exponent, power)
            primitive_cross[row, column] =
                term === :overlap ?
                _qwrg_atomic_basic_integral(
                    alpha_basis,
                    primitive.center_value,
                    0,
                    1.0,
                    exponent,
                    center_value,
                    power,
                    prefactor,
                ) :
                term === :kinetic ?
                _qwrg_atomic_kinetic_integral(
                    alpha_basis,
                    primitive.center_value,
                    0,
                    1.0,
                    exponent,
                    center_value,
                    power,
                    prefactor,
                ) :
                term === :position ?
                _qwrg_atomic_basic_integral(
                    alpha_basis,
                    primitive.center_value,
                    0,
                    1.0,
                    exponent,
                    center_value,
                    power,
                    prefactor;
                    xpower = 1,
                ) :
                term === :x2 ?
                _qwrg_atomic_basic_integral(
                    alpha_basis,
                    primitive.center_value,
                    0,
                    1.0,
                    exponent,
                    center_value,
                    power,
                    prefactor;
                    xpower = 2,
                ) :
                term === :factor ?
                _qwrg_atomic_basic_integral(
                    alpha_basis,
                    primitive.center_value,
                    0,
                    1.0,
                    exponent,
                    center_value,
                    power,
                    prefactor;
                    extra_exponent = Float64(factor_exponent),
                    extra_center = Float64(factor_center),
                ) :
                throw(ArgumentError("unsupported PQS/GTO axis cross term :$(term)"))
        end
    end
    return Matrix{Float64}(transpose(basis.coefficient_matrix) * primitive_cross)
end

function _pqs_source_box_route_driver_gto_axis_cross_tables(
    axes,
    orbital_arrays,
    term::Symbol;
    factor_center = nothing,
    factor_exponent = nothing,
)
    return ntuple(3) do index
        axis = _pqs_source_box_route_driver_axis_symbol(index)
        _pqs_source_box_route_driver_gto_axis_cross(
            getproperty(axes, axis),
            orbital_arrays,
            axis,
            term;
            factor_center =
                isnothing(factor_center) ? nothing : factor_center[index],
            factor_exponent,
        )
    end
end

function _pqs_source_box_route_driver_gto_axis_self(
    left_arrays,
    right_arrays,
    axis::Symbol,
    term::Symbol;
    factor_center::Union{Nothing,Real} = nothing,
    factor_exponent::Union{Nothing,Real} = nothing,
)
    axis_index = _pqs_source_box_route_driver_axis_index(axis)
    center_left_axis = @view left_arrays.centers[:, axis_index]
    center_right_axis = @view right_arrays.centers[:, axis_index]
    power_left_axis = @view left_arrays.angular_powers[:, axis_index]
    power_right_axis = @view right_arrays.angular_powers[:, axis_index]
    matrix = zeros(
        Float64,
        length(left_arrays.exponents),
        length(right_arrays.exponents),
    )
    for j in eachindex(right_arrays.exponents)
        exponent_right = right_arrays.exponents[j]
        center_right = center_right_axis[j]
        power_right = power_right_axis[j]
        prefactor_right = _qwrg_atomic_shell_prefactor(exponent_right, power_right)
        for i in eachindex(left_arrays.exponents)
            exponent_left = left_arrays.exponents[i]
            center_left = center_left_axis[i]
            power_left = power_left_axis[i]
            prefactor_left = _qwrg_atomic_shell_prefactor(exponent_left, power_left)
            matrix[i, j] =
                term === :overlap ?
                _qwrg_atomic_basic_integral(
                    exponent_left,
                    center_left,
                    power_left,
                    prefactor_left,
                    exponent_right,
                    center_right,
                    power_right,
                    prefactor_right,
                ) :
                term === :kinetic ?
                _qwrg_atomic_kinetic_integral(
                    exponent_left,
                    center_left,
                    power_left,
                    prefactor_left,
                    exponent_right,
                    center_right,
                    power_right,
                    prefactor_right,
                ) :
                term === :position ?
                _qwrg_atomic_basic_integral(
                    exponent_left,
                    center_left,
                    power_left,
                    prefactor_left,
                    exponent_right,
                    center_right,
                    power_right,
                    prefactor_right;
                    xpower = 1,
                ) :
                term === :x2 ?
                _qwrg_atomic_basic_integral(
                    exponent_left,
                    center_left,
                    power_left,
                    prefactor_left,
                    exponent_right,
                    center_right,
                    power_right,
                    prefactor_right;
                    xpower = 2,
                ) :
                term === :factor ?
                _qwrg_atomic_basic_integral(
                    exponent_left,
                    center_left,
                    power_left,
                    prefactor_left,
                    exponent_right,
                    center_right,
                    power_right,
                    prefactor_right;
                    extra_exponent = Float64(factor_exponent),
                    extra_center = Float64(factor_center),
                ) :
                throw(ArgumentError("unsupported GTO/GTO axis term :$(term)"))
        end
    end
    return matrix
end

function _pqs_source_box_route_driver_gto_axis_self_tables(
    left_arrays,
    right_arrays,
    term::Symbol;
    factor_center = nothing,
    factor_exponent = nothing,
)
    return ntuple(3) do index
        axis = _pqs_source_box_route_driver_axis_symbol(index)
        _pqs_source_box_route_driver_gto_axis_self(
            left_arrays,
            right_arrays,
            axis,
            term;
            factor_center =
                isnothing(factor_center) ? nothing : factor_center[index],
            factor_exponent,
        )
    end
end

function _pqs_source_box_route_driver_gto_support_column!(
    destination::AbstractVector{<:Real},
    states,
    coefficients::AbstractVector{<:Real},
    axis_tables,
)
    fill!(destination, 0.0)
    @inbounds for (row, (ix, iy, iz)) in pairs(states)
        value = 0.0
        for primitive in eachindex(coefficients)
            value +=
                coefficients[primitive] *
                axis_tables[1][ix, primitive] *
                axis_tables[2][iy, primitive] *
                axis_tables[3][iz, primitive]
        end
        destination[row] = value
    end
    return destination
end

function _pqs_source_box_route_driver_gto_weighted_hadamard(
    left_coefficients::AbstractVector{<:Real},
    x::AbstractMatrix{<:Real},
    y::AbstractMatrix{<:Real},
    z::AbstractMatrix{<:Real},
    right_coefficients::AbstractVector{<:Real},
)
    matrix = Matrix{Float64}(x) .* Matrix{Float64}(y) .* Matrix{Float64}(z)
    return Float64(dot(left_coefficients, matrix * right_coefficients))
end

function _pqs_source_box_route_driver_pqs_gto_support_one_body(packet)
    states = packet.support_states
    axes = packet.axis_representations
    coefficients = packet.coulomb_coefficients
    exponents = packet.coulomb_exponents
    nuclear_charges = packet.nuclear_charges
    atom_locations = packet.atom_locations
    support_gto_kinetic =
        zeros(Float64, length(states), length(packet.gto_primitive_arrays))
    support_gto_charged_nuclear =
        zeros(Float64, length(states), length(packet.gto_primitive_arrays))
    scratch = zeros(Float64, length(states))

    for (column, orbital_arrays) in pairs(packet.gto_primitive_arrays)
        overlap_tables =
            _pqs_source_box_route_driver_gto_axis_cross_tables(
                axes,
                orbital_arrays,
                :overlap,
            )
        kinetic_tables =
            _pqs_source_box_route_driver_gto_axis_cross_tables(
                axes,
                orbital_arrays,
                :kinetic,
            )
        for kinetic_axis in 1:3
            axis_tables = ntuple(
                axis -> axis == kinetic_axis ?
                    kinetic_tables[axis] :
                    overlap_tables[axis],
                3,
            )
            _pqs_source_box_route_driver_gto_support_column!(
                scratch,
                states,
                orbital_arrays.coefficients,
                axis_tables,
            )
            support_gto_kinetic[:, column] .+= scratch
        end
        for (center_index, location) in pairs(atom_locations)
            charge = nuclear_charges[center_index]
            for term_index in eachindex(coefficients)
                factor_tables =
                    _pqs_source_box_route_driver_gto_axis_cross_tables(
                        axes,
                        orbital_arrays,
                        :factor;
                        factor_center = location,
                        factor_exponent = exponents[term_index],
                    )
                _pqs_source_box_route_driver_gto_support_column!(
                    scratch,
                    states,
                    orbital_arrays.coefficients,
                    factor_tables,
                )
                support_gto_charged_nuclear[:, column] .-=
                    charge * coefficients[term_index] .* scratch
            end
        end
    end

    h_fg_kinetic =
        Matrix{Float64}(transpose(packet.final_coefficients) * support_gto_kinetic)
    h_fg_charged_nuclear =
        Matrix{Float64}(
            transpose(packet.final_coefficients) * support_gto_charged_nuclear,
        )
    return (;
        h_fg_kinetic,
        h_fg_charged_nuclear,
        h_fg_one_body = h_fg_kinetic + h_fg_charged_nuclear,
    )
end

function _pqs_source_box_route_driver_pqs_gto_self_one_body(packet)
    coefficients = packet.coulomb_coefficients
    exponents = packet.coulomb_exponents
    nuclear_charges = packet.nuclear_charges
    atom_locations = packet.atom_locations
    norbital = length(packet.gto_primitive_arrays)
    h_gg_kinetic = zeros(Float64, norbital, norbital)
    h_gg_charged_nuclear = zeros(Float64, norbital, norbital)
    for column in eachindex(packet.gto_primitive_arrays)
        right_arrays = packet.gto_primitive_arrays[column]
        for row in eachindex(packet.gto_primitive_arrays)
            left_arrays = packet.gto_primitive_arrays[row]
            overlap_tables =
                _pqs_source_box_route_driver_gto_axis_self_tables(
                    left_arrays,
                    right_arrays,
                    :overlap,
                )
            kinetic_tables =
                _pqs_source_box_route_driver_gto_axis_self_tables(
                    left_arrays,
                    right_arrays,
                    :kinetic,
                )
            for kinetic_axis in 1:3
                axis_tables = ntuple(
                    axis -> axis == kinetic_axis ?
                        kinetic_tables[axis] :
                        overlap_tables[axis],
                    3,
                )
                h_gg_kinetic[row, column] +=
                    _pqs_source_box_route_driver_gto_weighted_hadamard(
                        left_arrays.coefficients,
                        axis_tables[1],
                        axis_tables[2],
                        axis_tables[3],
                        right_arrays.coefficients,
                    )
            end
            for (center_index, location) in pairs(atom_locations)
                charge = nuclear_charges[center_index]
                for term_index in eachindex(coefficients)
                    factor_tables =
                        _pqs_source_box_route_driver_gto_axis_self_tables(
                            left_arrays,
                            right_arrays,
                            :factor;
                            factor_center = location,
                            factor_exponent = exponents[term_index],
                        )
                    h_gg_charged_nuclear[row, column] -=
                        charge *
                        coefficients[term_index] *
                        _pqs_source_box_route_driver_gto_weighted_hadamard(
                            left_arrays.coefficients,
                            factor_tables[1],
                            factor_tables[2],
                            factor_tables[3],
                            right_arrays.coefficients,
                        )
                end
            end
        end
    end
    h_gg_kinetic = Matrix{Float64}(0.5 .* (h_gg_kinetic .+ transpose(h_gg_kinetic)))
    h_gg_charged_nuclear =
        Matrix{Float64}(0.5 .* (h_gg_charged_nuclear .+ transpose(h_gg_charged_nuclear)))
    return (;
        h_gg_kinetic,
        h_gg_charged_nuclear,
        h_gg_one_body = h_gg_kinetic + h_gg_charged_nuclear,
    )
end

function _pqs_source_box_route_driver_pqs_gto_one_body_blocks(packet)
    h_ff = packet.h_ff
    s_fg = packet.s_fg
    l = packet.residual_transform
    mixed = _pqs_source_box_route_driver_pqs_gto_support_one_body(packet)
    self = _pqs_source_box_route_driver_pqs_gto_self_one_body(packet)
    h_fg = Matrix{Float64}(mixed.h_fg_one_body)
    h_gg = Matrix{Float64}(self.h_gg_one_body)
    h_fr = Matrix{Float64}((h_fg - h_ff * s_fg) * l)
    residual_core = Matrix{Float64}(
        h_gg -
        transpose(s_fg) * h_fg -
        transpose(h_fg) * s_fg +
        transpose(s_fg) * h_ff * s_fg,
    )
    h_rr = Matrix{Float64}(transpose(l) * residual_core * l)
    h_rr = Matrix{Float64}(0.5 .* (h_rr .+ transpose(h_rr)))
    augmented_one_body_hamiltonian = [
        h_ff h_fr
        transpose(h_fr) h_rr
    ]
    augmented_symmetry_error =
        norm(augmented_one_body_hamiltonian - transpose(augmented_one_body_hamiltonian), Inf)
    augmented_h1_lowest =
        minimum(eigvals(Symmetric(0.5 .* (
            augmented_one_body_hamiltonian + transpose(augmented_one_body_hamiltonian)
        ))))
    return (;
        h_fg_kinetic = mixed.h_fg_kinetic,
        h_fg_charged_nuclear = mixed.h_fg_charged_nuclear,
        h_gg_kinetic = self.h_gg_kinetic,
        h_gg_charged_nuclear = self.h_gg_charged_nuclear,
        h_fg_one_body = h_fg,
        h_gg_one_body = h_gg,
        h_fr_one_body = h_fr,
        h_rr_one_body = h_rr,
        augmented_one_body_hamiltonian,
        augmented_dimension = size(augmented_one_body_hamiltonian, 1),
        augmented_h1_lowest,
        augmented_h1_symmetry_error = augmented_symmetry_error,
        nuclear_mixed_block_convention = :charged_nuclear,
    )
end

function _pqs_source_box_route_driver_pqs_h2_route_metadata(report, inputs)
    final_basis = inputs.final_basis
    source_plan = inputs.source_plan
    recipe = report.recipe_metadata
    system = report.system_metadata
    return (;
        route_family = report.route_family,
        route_kind = report.route_kind,
        atom_symbols = system.atom_symbols,
        nuclear_charges = system.nuclear_charges,
        atom_locations = system.atom_locations,
        bond_axis = system.bond_axis,
        bond_length = system.bond_length,
        q = recipe.q,
        n_s = recipe.n_s,
        core_spacing = recipe.core_spacing,
        supplement_policy = recipe.supplement_policy,
        support_order = source_plan.support_order,
        retained_order = source_plan.retained_order,
        support_counts = final_basis.support_counts,
        retained_counts = final_basis.retained_counts,
        retained_ranges = final_basis.retained_ranges,
        final_dimension = final_basis.final_dimension,
        artifact_scope = :pqs_ham_basis_plus_residual_gto_sidecar,
        provider_blocks_included =
            hasproperty(inputs, :provider_blocks_included) ?
            inputs.provider_blocks_included :
            false,
    )
end

function _pqs_source_box_route_driver_require_jld2_keys(file, keys)
    missing = Symbol[key for key in keys if !haskey(file, String(key))]
    isempty(missing) || throw(ArgumentError("H2 PQS sidecar artifact missing keys $(missing)"))
    return nothing
end

function _pqs_source_box_route_driver_square_matrix(name::AbstractString, matrix)
    size(matrix, 1) == size(matrix, 2) ||
        throw(ArgumentError("$name must be square; got $(size(matrix))"))
    return nothing
end

function _pqs_source_box_route_driver_pqs_h2_provider_block_mode(mode)
    mode === false && return false
    mode === :one_body_only && return :one_body_only
    mode === :one_body_and_density_provider && return :one_body_and_density_provider
    throw(ArgumentError("provider_blocks_included must be false, :one_body_only, or :one_body_and_density_provider; got $(mode)"))
end

function _pqs_source_box_route_driver_finite_matrix(
    name::AbstractString,
    matrix,
    expected_size::Tuple{Int,Int},
)
    size(matrix) == expected_size ||
        throw(ArgumentError("$name shape must be $(expected_size); got $(size(matrix))"))
    all(isfinite, matrix) ||
        throw(ArgumentError("$name contains non-finite entries"))
    return nothing
end

function _pqs_source_box_route_driver_symmetric_matrix(
    name::AbstractString,
    matrix,
    expected_size::Tuple{Int,Int},
    symmetry_atol::Real,
)
    _pqs_source_box_route_driver_finite_matrix(name, matrix, expected_size)
    size(matrix, 1) == size(matrix, 2) ||
        throw(ArgumentError("$name must be square; got $(size(matrix))"))
    symmetry_error = norm(matrix - transpose(matrix), Inf)
    symmetry_error <= Float64(symmetry_atol) ||
        throw(ArgumentError("$name is not symmetric"))
    return symmetry_error
end

"""
    _pqs_source_box_route_driver_pqs_h2_residual_gto_sidecar_artifact_roundtrip(
        basisfile,
        hamfile,
    )

Read the narrow H2 independent-PQS Ham/Basis plus residual-GTO sidecar
artifacts and return compact consumer facts. This is intentionally specific to
the current P1 artifact; it is not a general artifact registry.
"""
function _pqs_source_box_route_driver_pqs_h2_residual_gto_sidecar_artifact_roundtrip(
    basisfile,
    hamfile;
    symmetry_atol::Real = 1.0e-8,
)
    basis_facts = jldopen(String(basisfile), "r") do file
        _pqs_source_box_route_driver_require_jld2_keys(
            file,
            (
                :artifact_kind,
                :final_dimension,
                :final_coefficients,
                :final_gto_cross_overlap,
                :gto_self_overlap,
                :gto_residual_overlap,
                :provider_blocks_included,
            ),
        )
        artifact_kind = file["artifact_kind"]
        artifact_kind === :pqs_h2_residual_gto_sidecar_basis_bundle ||
            throw(ArgumentError("unexpected H2 PQS sidecar basis artifact kind $(artifact_kind)"))
        final_dimension = Int(file["final_dimension"])
        final_coefficients = file["final_coefficients"]
        final_gto_cross_overlap = file["final_gto_cross_overlap"]
        gto_self_overlap = file["gto_self_overlap"]
        gto_residual_overlap = file["gto_residual_overlap"]
        provider_blocks_included =
            _pqs_source_box_route_driver_pqs_h2_provider_block_mode(
                file["provider_blocks_included"],
            )
        size(final_coefficients, 2) == final_dimension ||
            throw(ArgumentError("basis final_coefficients column count must equal final_dimension"))
        size(final_gto_cross_overlap, 1) == final_dimension ||
            throw(ArgumentError("basis final_gto_cross_overlap row count must equal final_dimension"))
        _pqs_source_box_route_driver_square_matrix(
            "basis gto_self_overlap",
            gto_self_overlap,
        )
        size(gto_residual_overlap) == size(gto_self_overlap) ||
            throw(ArgumentError("basis gto_residual_overlap shape must match gto_self_overlap"))
        all(isfinite, gto_residual_overlap) ||
            throw(ArgumentError("basis gto_residual_overlap contains non-finite entries"))
        residual_symmetry_error =
            norm(gto_residual_overlap - transpose(gto_residual_overlap), Inf)
        residual_symmetry_error <= Float64(symmetry_atol) ||
            throw(ArgumentError("basis gto_residual_overlap is not symmetric"))
        residual_facts =
            if provider_blocks_included !== false
                _pqs_source_box_route_driver_require_jld2_keys(
                    file,
                    (
                        :residual_transform,
                        :residual_rank,
                        :residual_overlap_eigenvalues,
                        :residual_overlap_identity_error,
                        :residual_overlap_cutoff,
                        :augmented_density_gauge,
                        :augmented_density_space,
                        :p_dimension,
                        :f_dimension,
                        :g_dimension,
                        :p_projection_of_g,
                        :residual_orbital_coefficients_in_density_carrier,
                    ),
                )
                residual_transform = file["residual_transform"]
                residual_rank = Int(file["residual_rank"])
                residual_overlap_identity_error =
                    Float64(file["residual_overlap_identity_error"])
                residual_rank > 0 ||
                    throw(ArgumentError("one-body provider artifact requires positive residual_rank"))
                size(residual_transform) == (size(gto_self_overlap, 1), residual_rank) ||
                    throw(ArgumentError("residual_transform shape must be gto_dimension x residual_rank"))
                isfinite(residual_overlap_identity_error) ||
                    throw(ArgumentError("residual_overlap_identity_error must be finite"))
                residual_overlap_identity_error <= Float64(symmetry_atol) ||
                    throw(ArgumentError("residual overlap identity error is too large"))
                file["augmented_density_gauge"] ===
                    :pre_final_localized_positive_weight ||
                    throw(ArgumentError("augmented density gauge must be pre-final PQS"))
                Tuple(file["augmented_density_space"]) ===
                    (:pre_final_pqs, :residual_gto) ||
                    throw(ArgumentError("augmented density space must be (:pre_final_pqs, :residual_gto)"))
                p_dimension = Int(file["p_dimension"])
                f_dimension = Int(file["f_dimension"])
                g_dimension = Int(file["g_dimension"])
                f_dimension == final_dimension ||
                    throw(ArgumentError("density descriptor F dimension must match final_dimension"))
                p_projection_of_g = file["p_projection_of_g"]
                residual_carrier =
                    file["residual_orbital_coefficients_in_density_carrier"]
                size(p_projection_of_g) == (p_dimension, g_dimension) ||
                    throw(ArgumentError("p_projection_of_g shape mismatch"))
                all(isfinite, p_projection_of_g) ||
                    throw(ArgumentError("p_projection_of_g contains non-finite entries"))
                size(residual_carrier) == (p_dimension + g_dimension, residual_rank) ||
                    throw(ArgumentError("residual density carrier coefficient shape mismatch"))
                all(isfinite, residual_carrier) ||
                    throw(ArgumentError("residual density carrier contains non-finite entries"))
                density_moment_facts =
                    if provider_blocks_included === :one_body_and_density_provider
                        _pqs_source_box_route_driver_require_jld2_keys(
                            file,
                            (
                                :residual_centers,
                                :residual_widths,
                                :residual_moment_overlap_error,
                                :residual_width_min,
                                :residual_width_max,
                                :residual_widths_positive,
                            ),
                        )
                        residual_centers = file["residual_centers"]
                        residual_widths = file["residual_widths"]
                        size(residual_centers) == (residual_rank, 3) ||
                            throw(ArgumentError("residual_centers shape mismatch"))
                        size(residual_widths) == (residual_rank, 3) ||
                            throw(ArgumentError("residual_widths shape mismatch"))
                        all(isfinite, residual_centers) ||
                            throw(ArgumentError("residual_centers contain non-finite entries"))
                        all(isfinite, residual_widths) ||
                            throw(ArgumentError("residual_widths contain non-finite entries"))
                        all(>(0.0), residual_widths) ||
                            throw(ArgumentError("residual_widths must be positive"))
                        residual_moment_overlap_error =
                            Float64(file["residual_moment_overlap_error"])
                        residual_moment_overlap_error <= Float64(symmetry_atol) ||
                            throw(ArgumentError("residual moment overlap error is too large"))
                        file["residual_widths_positive"] === true ||
                            throw(ArgumentError("residual_widths_positive must be true"))
                        (;
                            residual_centers_size = size(residual_centers),
                            residual_widths_size = size(residual_widths),
                            residual_moment_overlap_error,
                            residual_width_min = Float64(file["residual_width_min"]),
                            residual_width_max = Float64(file["residual_width_max"]),
                        )
                    else
                        (;
                            residual_centers_size = nothing,
                            residual_widths_size = nothing,
                            residual_moment_overlap_error = nothing,
                            residual_width_min = nothing,
                            residual_width_max = nothing,
                        )
                    end
                (;
                    residual_transform_size = size(residual_transform),
                    residual_rank,
                    residual_overlap_identity_error,
                    residual_overlap_cutoff = Float64(file["residual_overlap_cutoff"]),
                    provider_blocks_included,
                    augmented_density_gauge = file["augmented_density_gauge"],
                    augmented_density_space = Tuple(file["augmented_density_space"]),
                    p_dimension,
                    f_dimension,
                    g_dimension,
                    p_projection_of_g_size = size(p_projection_of_g),
                    residual_orbital_coefficients_in_density_carrier_size =
                        size(residual_carrier),
                    density_moment_facts...,
                )
            else
                (;
                    residual_transform_size = nothing,
                    residual_rank = nothing,
                    residual_overlap_identity_error = nothing,
                    residual_overlap_cutoff = nothing,
                    provider_blocks_included,
                    augmented_density_gauge = nothing,
                    augmented_density_space = nothing,
                    p_dimension = nothing,
                    f_dimension = nothing,
                    g_dimension = nothing,
                    p_projection_of_g_size = nothing,
                    residual_orbital_coefficients_in_density_carrier_size = nothing,
                )
            end
        return (;
            basis_artifact_kind = artifact_kind,
            final_dimension,
            provider_blocks_included,
            final_coefficients_size = size(final_coefficients),
            final_gto_cross_overlap_size = size(final_gto_cross_overlap),
            gto_self_overlap_size = size(gto_self_overlap),
            gto_residual_overlap_size = size(gto_residual_overlap),
            gto_residual_overlap_finite = true,
            gto_residual_overlap_symmetry_error = residual_symmetry_error,
            residual_facts...,
        )
    end

    ham_facts = jldopen(String(hamfile), "r") do file
        _pqs_source_box_route_driver_require_jld2_keys(
            file,
            (
                :artifact_kind,
                :final_dimension,
                :one_body_hamiltonian,
                :pre_final_pair_matrix,
                :h1_lowest,
                :h1_j_self_coulomb,
                :final_gto_cross_overlap,
                :gto_self_overlap,
                :gto_residual_overlap,
                :provider_blocks_included,
            ),
        )
        artifact_kind = file["artifact_kind"]
        artifact_kind === :pqs_h2_residual_gto_sidecar_ham_bundle ||
            throw(ArgumentError("unexpected H2 PQS sidecar ham artifact kind $(artifact_kind)"))
        final_dimension = Int(file["final_dimension"])
        final_dimension == basis_facts.final_dimension ||
            throw(ArgumentError("basis and ham final dimensions differ"))
        one_body_hamiltonian = file["one_body_hamiltonian"]
        pre_final_pair_matrix = file["pre_final_pair_matrix"]
        h1_lowest = Float64(file["h1_lowest"])
        h1_j_self_coulomb = Float64(file["h1_j_self_coulomb"])
        final_gto_cross_overlap = file["final_gto_cross_overlap"]
        gto_self_overlap = file["gto_self_overlap"]
        gto_residual_overlap = file["gto_residual_overlap"]
        provider_blocks_included =
            _pqs_source_box_route_driver_pqs_h2_provider_block_mode(
                file["provider_blocks_included"],
            )
        provider_blocks_included === basis_facts.provider_blocks_included ||
            throw(ArgumentError("basis and ham provider block modes differ"))
        size(one_body_hamiltonian) == (final_dimension, final_dimension) ||
            throw(ArgumentError("one_body_hamiltonian shape must be final_dimension x final_dimension"))
        all(isfinite, one_body_hamiltonian) ||
            throw(ArgumentError("one_body_hamiltonian contains non-finite entries"))
        h1_symmetry_error =
            norm(one_body_hamiltonian - transpose(one_body_hamiltonian), Inf)
        h1_symmetry_error <= Float64(symmetry_atol) ||
            throw(ArgumentError("one_body_hamiltonian is not symmetric"))
        isfinite(h1_lowest) ||
            throw(ArgumentError("h1_lowest must be finite"))
        _pqs_source_box_route_driver_square_matrix(
            "pre_final_pair_matrix",
            pre_final_pair_matrix,
        )
        isfinite(h1_j_self_coulomb) && h1_j_self_coulomb > 0 ||
            throw(ArgumentError("h1_j_self_coulomb must be finite and positive"))
        size(final_gto_cross_overlap, 1) == final_dimension ||
            throw(ArgumentError("ham final_gto_cross_overlap row count must equal final_dimension"))
        _pqs_source_box_route_driver_square_matrix(
            "ham gto_self_overlap",
            gto_self_overlap,
        )
        size(gto_residual_overlap) == size(gto_self_overlap) ||
            throw(ArgumentError("ham gto_residual_overlap shape must match gto_self_overlap"))
        all(isfinite, gto_residual_overlap) ||
            throw(ArgumentError("ham gto_residual_overlap contains non-finite entries"))
        residual_symmetry_error =
            norm(gto_residual_overlap - transpose(gto_residual_overlap), Inf)
        residual_symmetry_error <= Float64(symmetry_atol) ||
            throw(ArgumentError("ham gto_residual_overlap is not symmetric"))
        augmented_facts =
            if provider_blocks_included !== false
                _pqs_source_box_route_driver_require_jld2_keys(
                    file,
                    (
                        :augmented_dimension,
                        :augmented_one_body_hamiltonian,
                        :augmented_h1_lowest,
                        :augmented_h1_symmetry_error,
                        :h_fg_kinetic,
                        :h_fg_charged_nuclear,
                        :h_fg_one_body,
                        :h_gg_kinetic,
                        :h_gg_charged_nuclear,
                        :h_gg_one_body,
                        :h_fr_one_body,
                        :h_rr_one_body,
                        :nuclear_mixed_block_convention,
                        :augmented_density_gauge,
                        :augmented_density_space,
                        :p_dimension,
                        :f_dimension,
                        :g_dimension,
                        :p_projection_of_g_size,
                        :residual_orbital_coefficients_in_density_carrier_size,
                    ),
                )
                residual_rank = basis_facts.residual_rank
                gto_dimension = size(gto_self_overlap, 1)
                augmented_dimension = Int(file["augmented_dimension"])
                augmented_dimension == final_dimension + residual_rank ||
                    throw(ArgumentError("augmented_dimension must equal final_dimension + residual_rank"))
                augmented_one_body_hamiltonian =
                    file["augmented_one_body_hamiltonian"]
                augmented_h1_symmetry_error =
                    _pqs_source_box_route_driver_symmetric_matrix(
                        "augmented one-body Hamiltonian",
                        augmented_one_body_hamiltonian,
                        (augmented_dimension, augmented_dimension),
                        symmetry_atol,
                    )
                augmented_h1_lowest = Float64(file["augmented_h1_lowest"])
                isfinite(augmented_h1_lowest) ||
                    throw(ArgumentError("augmented_h1_lowest must be finite"))
                fg_size = (final_dimension, gto_dimension)
                gg_size = (gto_dimension, gto_dimension)
                fr_size = (final_dimension, residual_rank)
                rr_size = (residual_rank, residual_rank)
                h_fg_kinetic = file["h_fg_kinetic"]
                h_fg_charged_nuclear = file["h_fg_charged_nuclear"]
                h_fg_one_body = file["h_fg_one_body"]
                h_gg_kinetic = file["h_gg_kinetic"]
                h_gg_charged_nuclear = file["h_gg_charged_nuclear"]
                h_gg_one_body = file["h_gg_one_body"]
                h_fr_one_body = file["h_fr_one_body"]
                h_rr_one_body = file["h_rr_one_body"]
                _pqs_source_box_route_driver_finite_matrix(
                    "h_fg_kinetic",
                    h_fg_kinetic,
                    fg_size,
                )
                _pqs_source_box_route_driver_finite_matrix(
                    "h_fg_charged_nuclear",
                    h_fg_charged_nuclear,
                    fg_size,
                )
                _pqs_source_box_route_driver_finite_matrix(
                    "h_fg_one_body",
                    h_fg_one_body,
                    fg_size,
                )
                _pqs_source_box_route_driver_finite_matrix(
                    "h_gg_kinetic",
                    h_gg_kinetic,
                    gg_size,
                )
                _pqs_source_box_route_driver_finite_matrix(
                    "h_gg_charged_nuclear",
                    h_gg_charged_nuclear,
                    gg_size,
                )
                h_gg_one_body_symmetry_error =
                    _pqs_source_box_route_driver_symmetric_matrix(
                        "h_gg_one_body",
                        h_gg_one_body,
                        gg_size,
                        symmetry_atol,
                    )
                _pqs_source_box_route_driver_finite_matrix(
                    "h_fr_one_body",
                    h_fr_one_body,
                    fr_size,
                )
                h_rr_one_body_symmetry_error =
                    _pqs_source_box_route_driver_symmetric_matrix(
                        "h_rr_one_body",
                        h_rr_one_body,
                        rr_size,
                        symmetry_atol,
                    )
                file["nuclear_mixed_block_convention"] === :charged_nuclear ||
                    throw(
                        ArgumentError(
                            "nuclear_mixed_block_convention must be :charged_nuclear",
                        ),
                    )
                file["augmented_density_gauge"] ===
                    :pre_final_localized_positive_weight ||
                    throw(ArgumentError("ham augmented density gauge must be pre-final PQS"))
                Tuple(file["augmented_density_space"]) ===
                    (:pre_final_pqs, :residual_gto) ||
                    throw(ArgumentError("ham augmented density space must be (:pre_final_pqs, :residual_gto)"))
                Int(file["p_dimension"]) == basis_facts.p_dimension ||
                    throw(ArgumentError("basis and ham P dimensions differ"))
                Int(file["f_dimension"]) == final_dimension ||
                    throw(ArgumentError("ham F dimension must match final_dimension"))
                Int(file["g_dimension"]) == gto_dimension ||
                    throw(ArgumentError("ham G dimension must match GTO dimension"))
                Tuple(file["p_projection_of_g_size"]) ==
                    basis_facts.p_projection_of_g_size ||
                    throw(ArgumentError("basis and ham P-projection sizes differ"))
                Tuple(file["residual_orbital_coefficients_in_density_carrier_size"]) ==
                    basis_facts.residual_orbital_coefficients_in_density_carrier_size ||
                    throw(ArgumentError("basis and ham residual carrier sizes differ"))
                density_matrix_facts =
                    if provider_blocks_included === :one_body_and_density_provider
                        _pqs_source_box_route_driver_require_jld2_keys(
                            file,
                            (
                                :v_pr_pair_matrix,
                                :v_rr_pair_matrix,
                                :augmented_pair_matrix,
                                :augmented_density_dimension,
                                :augmented_pair_matrix_symmetry_error,
                                :residual_centers,
                                :residual_widths,
                                :residual_moment_overlap_error,
                                :augmented_h1_j_self_coulomb,
                                :augmented_h1_j_density_coefficients_length,
                            ),
                        )
                        augmented_density_dimension =
                            Int(file["augmented_density_dimension"])
                        augmented_density_dimension ==
                            basis_facts.p_dimension + residual_rank ||
                            throw(
                                ArgumentError(
                                    "augmented_density_dimension must equal p_dimension + residual_rank",
                                ),
                            )
                        v_pr_pair_matrix = file["v_pr_pair_matrix"]
                        v_rr_pair_matrix = file["v_rr_pair_matrix"]
                        augmented_pair_matrix = file["augmented_pair_matrix"]
                        _pqs_source_box_route_driver_finite_matrix(
                            "v_pr_pair_matrix",
                            v_pr_pair_matrix,
                            (basis_facts.p_dimension, residual_rank),
                        )
                        _pqs_source_box_route_driver_symmetric_matrix(
                            "v_rr_pair_matrix",
                            v_rr_pair_matrix,
                            (residual_rank, residual_rank),
                            symmetry_atol,
                        )
                        augmented_pair_matrix_symmetry_error =
                            _pqs_source_box_route_driver_symmetric_matrix(
                                "augmented_pair_matrix",
                                augmented_pair_matrix,
                                (
                                    augmented_density_dimension,
                                    augmented_density_dimension,
                                ),
                                symmetry_atol,
                            )
                        augmented_h1_j_self_coulomb =
                            Float64(file["augmented_h1_j_self_coulomb"])
                        isfinite(augmented_h1_j_self_coulomb) &&
                            augmented_h1_j_self_coulomb > 0 ||
                            throw(ArgumentError("augmented H1-J self Coulomb must be finite and positive"))
                        Int(file["augmented_h1_j_density_coefficients_length"]) ==
                            augmented_density_dimension ||
                            throw(ArgumentError("augmented H1-J density coefficient length mismatch"))
                        (;
                            augmented_density_dimension,
                            v_pr_pair_matrix_size = size(v_pr_pair_matrix),
                            v_rr_pair_matrix_size = size(v_rr_pair_matrix),
                            augmented_pair_matrix_size = size(augmented_pair_matrix),
                            augmented_pair_matrix_symmetry_error,
                            augmented_h1_j_self_coulomb,
                        )
                    else
                        (;
                            augmented_density_dimension = nothing,
                            v_pr_pair_matrix_size = nothing,
                            v_rr_pair_matrix_size = nothing,
                            augmented_pair_matrix_size = nothing,
                            augmented_pair_matrix_symmetry_error = nothing,
                            augmented_h1_j_self_coulomb = nothing,
                        )
                    end
                (;
                    augmented_dimension,
                    augmented_one_body_hamiltonian_size =
                        size(augmented_one_body_hamiltonian),
                    augmented_h1_lowest,
                    augmented_h1_symmetry_error,
                    h_fg_one_body_size = size(h_fg_one_body),
                    h_gg_one_body_size = size(h_gg_one_body),
                    h_gg_one_body_symmetry_error,
                    h_fr_one_body_size = size(h_fr_one_body),
                    h_rr_one_body_size = size(h_rr_one_body),
                    h_rr_one_body_symmetry_error,
                    nuclear_mixed_block_convention = :charged_nuclear,
                    augmented_density_gauge = file["augmented_density_gauge"],
                    augmented_density_space =
                        Tuple(file["augmented_density_space"]),
                    p_projection_of_g_size =
                        Tuple(file["p_projection_of_g_size"]),
                    residual_orbital_coefficients_in_density_carrier_size =
                        Tuple(file["residual_orbital_coefficients_in_density_carrier_size"]),
                    density_matrix_facts...,
                )
            else
                (;
                    augmented_dimension = nothing,
                    augmented_one_body_hamiltonian_size = nothing,
                    augmented_h1_lowest = nothing,
                    augmented_h1_symmetry_error = nothing,
                    augmented_h1_j_self_coulomb = nothing,
                    h_fg_one_body_size = nothing,
                    h_gg_one_body_size = nothing,
                    h_gg_one_body_symmetry_error = nothing,
                    h_fr_one_body_size = nothing,
                    h_rr_one_body_size = nothing,
                    h_rr_one_body_symmetry_error = nothing,
                    nuclear_mixed_block_convention = nothing,
                    augmented_density_gauge = nothing,
                    augmented_density_space = nothing,
                    p_projection_of_g_size = nothing,
                    residual_orbital_coefficients_in_density_carrier_size = nothing,
                )
            end
        return (;
            ham_artifact_kind = artifact_kind,
            provider_blocks_included,
            one_body_hamiltonian_size = size(one_body_hamiltonian),
            h1_finite = true,
            h1_symmetry_error,
            h1_lowest,
            pre_final_pair_matrix_size = size(pre_final_pair_matrix),
            h1_j_self_coulomb,
            h1_j_self_coulomb_positive = true,
            final_gto_cross_overlap_size = size(final_gto_cross_overlap),
            gto_self_overlap_size = size(gto_self_overlap),
            gto_residual_overlap_size = size(gto_residual_overlap),
            gto_residual_overlap_finite = true,
            gto_residual_overlap_symmetry_error = residual_symmetry_error,
            augmented_facts...,
        )
    end

    return (;
        basisfile = String(basisfile),
        hamfile = String(hamfile),
        final_dimension = basis_facts.final_dimension,
        basis_artifact_kind = basis_facts.basis_artifact_kind,
        ham_artifact_kind = ham_facts.ham_artifact_kind,
        provider_blocks_included = basis_facts.provider_blocks_included,
        final_coefficients_size = basis_facts.final_coefficients_size,
        final_gto_cross_overlap_size = basis_facts.final_gto_cross_overlap_size,
        gto_self_overlap_size = basis_facts.gto_self_overlap_size,
        gto_residual_overlap_size = basis_facts.gto_residual_overlap_size,
        h1_finite = ham_facts.h1_finite,
        h1_symmetry_error = ham_facts.h1_symmetry_error,
        h1_lowest = ham_facts.h1_lowest,
        h1_j_self_coulomb = ham_facts.h1_j_self_coulomb,
        h1_j_self_coulomb_positive = ham_facts.h1_j_self_coulomb_positive,
        residual_rank = basis_facts.residual_rank,
        residual_overlap_identity_error =
            basis_facts.residual_overlap_identity_error,
        augmented_density_gauge = basis_facts.augmented_density_gauge,
        augmented_density_space = basis_facts.augmented_density_space,
        p_dimension = basis_facts.p_dimension,
        f_dimension = basis_facts.f_dimension,
        g_dimension = basis_facts.g_dimension,
        p_projection_of_g_size = basis_facts.p_projection_of_g_size,
        residual_orbital_coefficients_in_density_carrier_size =
            basis_facts.residual_orbital_coefficients_in_density_carrier_size,
        augmented_dimension = ham_facts.augmented_dimension,
        augmented_one_body_hamiltonian_size =
            ham_facts.augmented_one_body_hamiltonian_size,
        augmented_h1_lowest = ham_facts.augmented_h1_lowest,
        augmented_h1_symmetry_error = ham_facts.augmented_h1_symmetry_error,
        augmented_h1_j_self_coulomb =
            ham_facts.augmented_h1_j_self_coulomb,
        gto_residual_overlap_finite =
            basis_facts.gto_residual_overlap_finite &&
            ham_facts.gto_residual_overlap_finite,
        basis_gto_residual_overlap_symmetry_error =
            basis_facts.gto_residual_overlap_symmetry_error,
        ham_gto_residual_overlap_symmetry_error =
            ham_facts.gto_residual_overlap_symmetry_error,
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
            "H2 PQS residual-GTO sidecar materialization requires source plan, final basis, H1, H1-J, and supplement representation",
        ),
    )

    residual_gto_provider_blocks in
        (:none, :one_body_only, :one_body_and_density_provider) || throw(
        ArgumentError(
            "residual_gto_provider_blocks must be :none, :one_body_only, or :one_body_and_density_provider",
        ),
    )
    provider_blocks_included =
        residual_gto_provider_blocks === :none ? false : residual_gto_provider_blocks
    sidecar = _pqs_source_box_route_driver_pqs_gto_sidecar(inputs)
    route_metadata = _pqs_source_box_route_driver_pqs_h2_route_metadata(
        report,
        merge(inputs, (; provider_blocks_included,)),
    )
    inputs = merge(inputs, (; route_metadata, provider_blocks_included,))
    residual =
        residual_gto_provider_blocks !== :none ?
        _pqs_source_box_route_driver_pqs_gto_residual_transform(sidecar) :
        nothing
    provider_packet =
        residual_gto_provider_blocks !== :none ?
        _pqs_source_box_route_driver_pqs_h2_residual_gto_provider_packet(
            inputs,
            sidecar,
            residual,
        ) :
        nothing
    density_descriptor =
        residual_gto_provider_blocks !== :none ?
        _pqs_source_box_route_driver_pqs_h2_residual_gto_density_descriptor(
            inputs,
            sidecar,
            residual,
        ) :
        nothing
    one_body_blocks =
        residual_gto_provider_blocks !== :none ?
        _pqs_source_box_route_driver_pqs_gto_one_body_blocks(provider_packet) :
        nothing
    density_blocks =
        residual_gto_provider_blocks === :one_body_and_density_provider ?
        _pqs_source_box_route_driver_pqs_h2_residual_gto_density_blocks(
            provider_packet,
            density_descriptor,
        ) :
        nothing
    augmented_h1_j =
        residual_gto_provider_blocks === :one_body_and_density_provider ?
        _pqs_source_box_route_driver_pqs_h2_augmented_h1_j_diagnostic(
            one_body_blocks,
            density_blocks,
            provider_packet.density_interaction,
        ) :
        nothing
    private_augmented_rhf =
        residual_gto_provider_blocks === :one_body_and_density_provider ?
        _pqs_source_box_route_driver_pqs_h2_private_augmented_rhf_smoke(
            one_body_blocks,
            density_blocks,
            provider_packet.density_interaction,
            route_metadata,
        ) :
        nothing
    ham_handoff =
        residual_gto_provider_blocks === :one_body_and_density_provider ?
        _pqs_source_box_route_driver_pqs_h2_residual_gto_ham_handoff(
            one_body_blocks,
            density_blocks,
            provider_packet.density_interaction,
            augmented_h1_j,
            private_augmented_rhf,
        ) :
        nothing
    final_basis = inputs.final_basis
    h1_hamiltonian = inputs.h1_hamiltonian
    h1 = inputs.h1
    density_interaction = inputs.density_interaction
    h1_j = inputs.h1_j
    h1_matrix = Matrix{Float64}(h1_hamiltonian.hamiltonian_matrix)
    h1_symmetry_error = norm(h1_matrix - transpose(h1_matrix), Inf)
    h1_finite = all(isfinite, h1_matrix)
    h1_j_self_coulomb = h1_j.self_coulomb
    final_gto_cross_overlap_size = size(sidecar.final_gto_cross_overlap)
    gto_self_overlap_size = size(sidecar.gto_self_overlap)
    gto_residual_overlap_size = size(sidecar.gto_residual_overlap)
    augmented_dimension =
        isnothing(one_body_blocks) ? nothing : one_body_blocks.augmented_dimension
    augmented_h1_lowest =
        isnothing(one_body_blocks) ? nothing : one_body_blocks.augmented_h1_lowest
    augmented_h1_symmetry_error =
        isnothing(one_body_blocks) ? nothing :
        one_body_blocks.augmented_h1_symmetry_error
    summary = (;
        route_family = report.route_family,
        route_kind = report.route_kind,
        result_kind = :h2_pqs_ham_with_residual_gto_sidecar,
        final_dimension = final_basis.final_dimension,
        overlap_identity_error = final_basis.final_overlap_identity_error,
        h1_lowest = h1.lowest_energy,
        h1_finite,
        h1_symmetry_error,
        h1_j_self_coulomb,
        h1_j_self_coulomb_positive =
            isfinite(h1_j_self_coulomb) && h1_j_self_coulomb > 0,
        supplement_policy = report.recipe_metadata.supplement_policy,
        gto_dimension = sidecar.diagnostics.gto_dimension,
        final_gto_cross_overlap_size,
        gto_self_overlap_size,
        gto_residual_overlap_size,
        gto_residual_overlap_eigenvalue_min =
            sidecar.diagnostics.gto_residual_overlap_eigenvalue_min,
        gto_residual_overlap_eigenvalue_max =
            sidecar.diagnostics.gto_residual_overlap_eigenvalue_max,
        provider_blocks_included,
        residual_rank = isnothing(residual) ? nothing : residual.residual_rank,
        residual_overlap_identity_error =
            isnothing(residual) ? nothing : residual.residual_overlap_identity_error,
        residual_overlap_cutoff =
            isnothing(residual) ? nothing : residual.residual_overlap_cutoff,
        augmented_dimension,
        augmented_h1_lowest,
        augmented_h1_symmetry_error,
        nuclear_mixed_block_convention =
            isnothing(one_body_blocks) ?
            nothing :
            one_body_blocks.nuclear_mixed_block_convention,
        h_fg_kinetic_size =
            isnothing(one_body_blocks) ? nothing : size(one_body_blocks.h_fg_kinetic),
        h_fg_charged_nuclear_size =
            isnothing(one_body_blocks) ?
            nothing :
            size(one_body_blocks.h_fg_charged_nuclear),
        h_gg_kinetic_size =
            isnothing(one_body_blocks) ? nothing : size(one_body_blocks.h_gg_kinetic),
        h_gg_charged_nuclear_size =
            isnothing(one_body_blocks) ?
            nothing :
            size(one_body_blocks.h_gg_charged_nuclear),
        augmented_density_space =
            isnothing(density_descriptor) ?
            nothing :
            density_descriptor.augmented_density_space,
        augmented_density_gauge =
            isnothing(density_descriptor) ? nothing : density_descriptor.density_gauge,
        p_projection_of_g_size =
            isnothing(density_descriptor) ?
            nothing :
            density_descriptor.p_projection_of_g_size,
        residual_orbital_coefficients_in_density_carrier_size =
            isnothing(density_descriptor) ?
            nothing :
            density_descriptor.residual_orbital_coefficients_in_density_carrier_size,
        augmented_density_dimension =
            isnothing(density_blocks) ?
            nothing :
            density_blocks.augmented_density_dimension,
        v_pr_pair_matrix_size =
            isnothing(density_blocks) ? nothing : density_blocks.v_pr_pair_matrix_size,
        v_rr_pair_matrix_size =
            isnothing(density_blocks) ? nothing : density_blocks.v_rr_pair_matrix_size,
        augmented_pair_matrix_size =
            isnothing(density_blocks) ?
            nothing :
            density_blocks.augmented_pair_matrix_size,
        augmented_pair_matrix_symmetry_error =
            isnothing(density_blocks) ?
            nothing :
            density_blocks.augmented_pair_matrix_symmetry_error,
        augmented_h1_j_self_coulomb =
            isnothing(augmented_h1_j) ? nothing : augmented_h1_j.self_coulomb,
        residual_width_min =
            isnothing(density_blocks) ? nothing : density_blocks.residual_width_min,
        residual_width_max =
            isnothing(density_blocks) ? nothing : density_blocks.residual_width_max,
        ham_handoff_kind =
            isnothing(ham_handoff) ? nothing : ham_handoff.handoff_kind,
        ham_handoff_orbital_basis =
            isnothing(ham_handoff) ? nothing : ham_handoff.orbital_basis,
        ham_handoff_density_basis =
            isnothing(ham_handoff) ? nothing : ham_handoff.density_basis,
        ham_handoff_orbital_dimension =
            isnothing(ham_handoff) ?
            nothing :
            ham_handoff.diagnostics.orbital_dimension,
        ham_handoff_density_dimension =
            isnothing(ham_handoff) ?
            nothing :
            ham_handoff.diagnostics.density_dimension,
    )

    basis_artifact_path = nothing
    if save_basis_artifact
        basis_artifact_path = String(basisfile)
        jldopen(basis_artifact_path, "w") do file
            file["artifact_kind"] = :pqs_h2_residual_gto_sidecar_basis_bundle
            file["summary"] = summary
            file["route_metadata"] = route_metadata
            file["support_counts"] = final_basis.support_counts
            file["retained_counts"] = final_basis.retained_counts
            file["retained_ranges"] = final_basis.retained_ranges
            file["final_dimension"] = final_basis.final_dimension
            file["final_coefficients"] = final_basis.final_coefficients
            file["pre_final_coefficients"] = final_basis.pre_final_coefficients
            file["final_overlap_identity_error"] =
                final_basis.final_overlap_identity_error
            file["gto_supplement_metadata"] =
                inputs.supplement_representation.metadata
            file["final_gto_cross_overlap"] = sidecar.final_gto_cross_overlap
            file["gto_self_overlap"] = sidecar.gto_self_overlap
            file["gto_residual_overlap"] = sidecar.gto_residual_overlap
            file["gto_sidecar_diagnostics"] = sidecar.diagnostics
            file["provider_blocks_included"] = provider_blocks_included
            if !isnothing(residual)
                file["residual_transform"] = residual.residual_transform
                file["residual_rank"] = residual.residual_rank
                file["residual_overlap_eigenvalues"] =
                    residual.residual_overlap_eigenvalues
                file["residual_overlap_identity_error"] =
                    residual.residual_overlap_identity_error
                file["residual_overlap_eigenvalue_min"] =
                    residual.residual_overlap_eigenvalue_min
                file["residual_overlap_eigenvalue_max"] =
                    residual.residual_overlap_eigenvalue_max
                file["residual_overlap_cutoff"] =
                    residual.residual_overlap_cutoff
            end
            if !isnothing(density_descriptor)
                file["augmented_density_gauge"] =
                    density_descriptor.density_gauge
                file["augmented_density_space"] =
                    density_descriptor.augmented_density_space
                file["p_dimension"] = density_descriptor.p_dimension
                file["f_dimension"] = density_descriptor.f_dimension
                file["g_dimension"] = density_descriptor.g_dimension
                file["p_projection_of_g"] =
                    density_descriptor.p_projection_of_g
                file["residual_orbital_coefficients_in_density_carrier"] =
                    density_descriptor.residual_orbital_coefficients_in_density_carrier
                if !isnothing(density_blocks)
                    file["residual_centers"] = density_blocks.residual_centers
                    file["residual_widths"] = density_blocks.residual_widths
                    file["residual_moment_overlap_error"] =
                        density_blocks.residual_overlap_error
                    file["residual_width_min"] = density_blocks.residual_width_min
                    file["residual_width_max"] = density_blocks.residual_width_max
                    file["residual_widths_positive"] =
                        density_blocks.residual_widths_positive
                end
            end
        end
    end

    ham_artifact_path = nothing
    if save_ham_artifact
        ham_artifact_path = String(hamfile)
        jldopen(ham_artifact_path, "w") do file
            file["artifact_kind"] = :pqs_h2_residual_gto_sidecar_ham_bundle
            file["summary"] = summary
            file["route_metadata"] = route_metadata
            file["final_dimension"] = final_basis.final_dimension
            file["one_body_hamiltonian"] = h1_hamiltonian.hamiltonian_matrix
            file["kinetic"] = h1_hamiltonian.kinetic_matrix
            file["charged_nuclear"] = h1_hamiltonian.charged_nuclear_matrix
            file["pre_final_pair_matrix"] =
                density_interaction.pre_final_pair_matrix
            file["final_to_pre_final_coefficients"] =
                density_interaction.final_to_pre_final_coefficients
            file["h1_lowest"] = h1.lowest_energy
            file["h1_j_self_coulomb"] = h1_j_self_coulomb
            file["gto_supplement_metadata"] =
                inputs.supplement_representation.metadata
            file["final_gto_cross_overlap"] = sidecar.final_gto_cross_overlap
            file["gto_self_overlap"] = sidecar.gto_self_overlap
            file["gto_residual_overlap"] = sidecar.gto_residual_overlap
            file["gto_sidecar_diagnostics"] = sidecar.diagnostics
            file["provider_blocks_included"] = provider_blocks_included
            if !isnothing(residual)
                file["residual_transform"] = residual.residual_transform
                file["residual_rank"] = residual.residual_rank
                file["residual_overlap_eigenvalues"] =
                    residual.residual_overlap_eigenvalues
                file["residual_overlap_identity_error"] =
                    residual.residual_overlap_identity_error
                file["residual_overlap_eigenvalue_min"] =
                    residual.residual_overlap_eigenvalue_min
                file["residual_overlap_eigenvalue_max"] =
                    residual.residual_overlap_eigenvalue_max
                file["residual_overlap_cutoff"] =
                    residual.residual_overlap_cutoff
            end
            if !isnothing(density_descriptor)
                file["augmented_density_gauge"] =
                    density_descriptor.density_gauge
                file["augmented_density_space"] =
                    density_descriptor.augmented_density_space
                file["p_dimension"] = density_descriptor.p_dimension
                file["f_dimension"] = density_descriptor.f_dimension
                file["g_dimension"] = density_descriptor.g_dimension
                file["p_projection_of_g_size"] =
                    density_descriptor.p_projection_of_g_size
                file["residual_orbital_coefficients_in_density_carrier_size"] =
                    density_descriptor.residual_orbital_coefficients_in_density_carrier_size
                if !isnothing(density_blocks)
                    file["v_pr_pair_matrix"] = density_blocks.v_pr_pair_matrix
                    file["v_rr_pair_matrix"] = density_blocks.v_rr_pair_matrix
                    file["augmented_pair_matrix"] =
                        density_blocks.augmented_pair_matrix
                    file["augmented_density_dimension"] =
                        density_blocks.augmented_density_dimension
                    file["augmented_pair_matrix_symmetry_error"] =
                        density_blocks.augmented_pair_matrix_symmetry_error
                    file["residual_centers"] = density_blocks.residual_centers
                    file["residual_widths"] = density_blocks.residual_widths
                    file["residual_moment_overlap_error"] =
                        density_blocks.residual_overlap_error
                    file["augmented_h1_j_self_coulomb"] =
                        augmented_h1_j.self_coulomb
                    file["augmented_h1_j_density_coefficients_length"] =
                        augmented_h1_j.density_coefficients_length
                    file["private_augmented_rhf_kind"] =
                        private_augmented_rhf.rhf_kind
                    file["private_augmented_rhf_converged"] =
                        private_augmented_rhf.converged
                    file["private_augmented_rhf_iterations"] =
                        private_augmented_rhf.iteration_count
                    file["private_augmented_rhf_density_trace"] =
                        private_augmented_rhf.density_trace
                    file["private_augmented_rhf_idempotency_error"] =
                        private_augmented_rhf.idempotency_error
                    file["private_augmented_rhf_commutator_residual"] =
                        private_augmented_rhf.commutator_residual
                    file["private_augmented_rhf_one_body_energy"] =
                        private_augmented_rhf.one_body_energy
                    file["private_augmented_rhf_two_body_energy"] =
                        private_augmented_rhf.two_body_energy
                    file["private_augmented_rhf_electronic_energy"] =
                        private_augmented_rhf.electronic_energy
                    file["private_augmented_rhf_nuclear_repulsion"] =
                        private_augmented_rhf.nuclear_repulsion
                    file["private_augmented_rhf_total_with_nuclear_repulsion"] =
                        private_augmented_rhf.total_with_nuclear_repulsion
                end
            end
            if !isnothing(one_body_blocks)
                file["augmented_dimension"] = one_body_blocks.augmented_dimension
                file["h_fg_kinetic"] = one_body_blocks.h_fg_kinetic
                file["h_fg_charged_nuclear"] =
                    one_body_blocks.h_fg_charged_nuclear
                file["h_gg_kinetic"] = one_body_blocks.h_gg_kinetic
                file["h_gg_charged_nuclear"] =
                    one_body_blocks.h_gg_charged_nuclear
                file["h_fg_one_body"] = one_body_blocks.h_fg_one_body
                file["h_gg_one_body"] = one_body_blocks.h_gg_one_body
                file["h_fr_one_body"] = one_body_blocks.h_fr_one_body
                file["h_rr_one_body"] = one_body_blocks.h_rr_one_body
                file["augmented_one_body_hamiltonian"] =
                    one_body_blocks.augmented_one_body_hamiltonian
                file["augmented_h1_lowest"] = one_body_blocks.augmented_h1_lowest
                file["augmented_h1_symmetry_error"] =
                    one_body_blocks.augmented_h1_symmetry_error
                file["nuclear_mixed_block_convention"] =
                    one_body_blocks.nuclear_mixed_block_convention
            end
        end
    end

    artifact_roundtrip =
        save_basis_artifact && save_ham_artifact ?
        _pqs_source_box_route_driver_pqs_h2_residual_gto_sidecar_artifact_roundtrip(
            basis_artifact_path,
            ham_artifact_path,
        ) :
        nothing
    basis_artifact_reloaded = !isnothing(artifact_roundtrip)
    ham_artifact_reloaded = !isnothing(artifact_roundtrip)

    return (;
        route_family = report.route_family,
        result_kind = :h2_pqs_ham_with_residual_gto_sidecar,
        requested = true,
        materialized = true,
        summary,
        final_dimension = final_basis.final_dimension,
        retained_dimension = report.retained_dimension,
        h1_lowest = h1.lowest_energy,
        h1_finite,
        h1_symmetry_error,
        h1_j_self_coulomb,
        h1_j_self_coulomb_positive = summary.h1_j_self_coulomb_positive,
        overlap_identity_error = final_basis.final_overlap_identity_error,
        gto_dimension = sidecar.diagnostics.gto_dimension,
        final_gto_cross_overlap_size,
        gto_self_overlap_size,
        gto_residual_overlap_size,
        gto_residual_overlap_symmetry_error =
            sidecar.diagnostics.gto_residual_overlap_symmetry_error,
        gto_residual_overlap_eigenvalue_min =
            sidecar.diagnostics.gto_residual_overlap_eigenvalue_min,
        gto_residual_overlap_eigenvalue_max =
            sidecar.diagnostics.gto_residual_overlap_eigenvalue_max,
        provider_blocks_included,
        residual_rank = isnothing(residual) ? nothing : residual.residual_rank,
        residual_overlap_identity_error =
            isnothing(residual) ? nothing : residual.residual_overlap_identity_error,
        augmented_density_gauge =
            isnothing(density_descriptor) ? nothing : density_descriptor.density_gauge,
        augmented_density_dimension =
            isnothing(density_blocks) ?
            nothing :
            density_blocks.augmented_density_dimension,
        augmented_pair_matrix_symmetry_error =
            isnothing(density_blocks) ?
            nothing :
            density_blocks.augmented_pair_matrix_symmetry_error,
        augmented_h1_j_self_coulomb =
            isnothing(augmented_h1_j) ? nothing : augmented_h1_j.self_coulomb,
        private_augmented_rhf_kind =
            isnothing(private_augmented_rhf) ?
            nothing :
            private_augmented_rhf.rhf_kind,
        private_augmented_rhf_converged =
            isnothing(private_augmented_rhf) ?
            nothing :
            private_augmented_rhf.converged,
        private_augmented_rhf_iterations =
            isnothing(private_augmented_rhf) ?
            nothing :
            private_augmented_rhf.iteration_count,
        private_augmented_rhf_electronic_energy =
            isnothing(private_augmented_rhf) ?
            nothing :
            private_augmented_rhf.electronic_energy,
        private_augmented_rhf_total_with_nuclear_repulsion =
            isnothing(private_augmented_rhf) ?
            nothing :
            private_augmented_rhf.total_with_nuclear_repulsion,
        private_augmented_rhf_commutator_residual =
            isnothing(private_augmented_rhf) ?
            nothing :
            private_augmented_rhf.commutator_residual,
        ham_handoff,
        ham_handoff_kind =
            isnothing(ham_handoff) ? nothing : ham_handoff.handoff_kind,
        ham_handoff_orbital_basis =
            isnothing(ham_handoff) ? nothing : ham_handoff.orbital_basis,
        ham_handoff_density_basis =
            isnothing(ham_handoff) ? nothing : ham_handoff.density_basis,
        ham_handoff_orbital_dimension =
            isnothing(ham_handoff) ?
            nothing :
            ham_handoff.diagnostics.orbital_dimension,
        ham_handoff_density_dimension =
            isnothing(ham_handoff) ?
            nothing :
            ham_handoff.diagnostics.density_dimension,
        augmented_dimension,
        augmented_h1_lowest,
        augmented_h1_symmetry_error,
        save_basis_artifact_requested = save_basis_artifact,
        save_ham_artifact_requested = save_ham_artifact,
        basisfile = basis_artifact_path,
        hamfile = ham_artifact_path,
        basis_artifact_written = save_basis_artifact,
        ham_artifact_written = save_ham_artifact,
        basis_artifact_reloaded,
        ham_artifact_reloaded,
        artifact_roundtrip,
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
