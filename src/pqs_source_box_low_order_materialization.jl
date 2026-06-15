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

    fixed_block =
        one_center_atomic_full_parent_fixed_block(
            axis_basis;
            expansion,
            nside,
            gausslet_backend = backend,
            timing = false,
        )
    nuclear_one_body = -z .* fixed_block.gaussian_sum
    h1 = fixed_block.kinetic .+ nuclear_one_body
    h1_symmetric = 0.5 .* (h1 .+ transpose(h1))
    dim = size(h1, 1)
    overlap_identity_error =
        maximum(abs.(fixed_block.overlap .- Matrix{Float64}(I, dim, dim)))
    h1_symmetry_error = maximum(abs.(h1 .- transpose(h1)))
    h1_lowest = minimum(eigvals(Symmetric(h1_symmetric)))
    summary = (;
        route_family = report.route_family,
        geometry = :atomic,
        supplement_policy = :none,
        source = :one_center_atomic_full_parent_fixed_block,
        backend,
        nside,
        nuclear_charge = z,
        support_dimension = size(fixed_block.coefficient_matrix, 1),
        retained_dimension = dim,
        coefficient_matrix_size = size(fixed_block.coefficient_matrix),
        overlap_identity_error,
        h1_lowest,
        h1_finite = all(isfinite, h1),
        h1_symmetry_error,
        density_density_pair_sum_present = hasproperty(fixed_block, :pair_sum),
    )
    basis_artifact_path = nothing
    if save_basis_artifact
        basis_artifact_path = String(basisfile)
        jldopen(basis_artifact_path, "w") do file
            file["artifact_kind"] = :white_lindsey_atomic_pure_gausslet_basis
            file["summary"] = summary
            file["coefficient_matrix"] = fixed_block.coefficient_matrix
            file["support_indices"] = fixed_block.support_indices
            file["overlap"] = fixed_block.overlap
            file["weights"] = fixed_block.weights
            file["fixed_centers"] = fixed_block.fixed_centers
        end
    end
    ham_artifact_path = nothing
    if save_ham_artifact
        ham_artifact_path = String(hamfile)
        jldopen(ham_artifact_path, "w") do file
            file["artifact_kind"] = :white_lindsey_atomic_pure_gausslet_hamiltonian
            file["summary"] = summary
            file["overlap"] = fixed_block.overlap
            file["kinetic"] = fixed_block.kinetic
            file["nuclear_one_body"] = nuclear_one_body
            file["one_body_hamiltonian"] = h1
            file["density_density_pair_sum"] = fixed_block.pair_sum
            file["weights"] = fixed_block.weights
            file["fixed_centers"] = fixed_block.fixed_centers
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
        gto_self_overlap,
        gto_residual_overlap,
        diagnostics,
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
        provider_blocks_included = false,
    )
end

function _pqs_source_box_route_driver_jld2_contains(file, required_keys)
    return all(key -> haskey(file, String(key)), required_keys)
end

function _pqs_source_box_route_driver_h2_pqs_basis_artifact_reload_ok(path)
    jldopen(String(path), "r") do file
        _pqs_source_box_route_driver_jld2_contains(
            file,
            (
                :artifact_kind,
                :final_dimension,
                :final_coefficients,
                :final_gto_cross_overlap,
                :gto_self_overlap,
                :gto_residual_overlap,
            ),
        ) || return false
        final_dimension = Int(file["final_dimension"])
        final_dimension == 471 || return false
        final_coefficients = file["final_coefficients"]
        final_gto_cross_overlap = file["final_gto_cross_overlap"]
        gto_self_overlap = file["gto_self_overlap"]
        gto_residual_overlap = file["gto_residual_overlap"]
        size(final_coefficients, 2) == final_dimension || return false
        size(final_gto_cross_overlap, 1) == final_dimension || return false
        size(gto_self_overlap, 1) == size(gto_self_overlap, 2) || return false
        size(gto_residual_overlap) == size(gto_self_overlap) || return false
        return true
    end
end

function _pqs_source_box_route_driver_h2_pqs_ham_artifact_reload_ok(
    path,
    final_dimension::Int,
)
    jldopen(String(path), "r") do file
        _pqs_source_box_route_driver_jld2_contains(
            file,
            (
                :artifact_kind,
                :one_body_hamiltonian,
                :pre_final_pair_matrix,
                :h1_j_self_coulomb,
            ),
        ) || return false
        one_body_hamiltonian = file["one_body_hamiltonian"]
        pre_final_pair_matrix = file["pre_final_pair_matrix"]
        h1_j_self_coulomb = Float64(file["h1_j_self_coulomb"])
        size(one_body_hamiltonian) == (final_dimension, final_dimension) || return false
        size(pre_final_pair_matrix, 1) == size(pre_final_pair_matrix, 2) ||
            return false
        isfinite(h1_j_self_coulomb) || return false
        return true
    end
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

"""
    pqs_h2_residual_gto_sidecar_artifact_roundtrip(basisfile, hamfile)

Read the narrow H2 independent-PQS Ham/Basis plus residual-GTO sidecar
artifacts and return compact consumer facts. This is intentionally specific to
the current P1 artifact; it is not a general artifact registry.
"""
function pqs_h2_residual_gto_sidecar_artifact_roundtrip(
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
        return (;
            basis_artifact_kind = artifact_kind,
            final_dimension,
            final_coefficients_size = size(final_coefficients),
            final_gto_cross_overlap_size = size(final_gto_cross_overlap),
            gto_self_overlap_size = size(gto_self_overlap),
            gto_residual_overlap_size = size(gto_residual_overlap),
            gto_residual_overlap_finite = true,
            gto_residual_overlap_symmetry_error = residual_symmetry_error,
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
        return (;
            ham_artifact_kind = artifact_kind,
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
        )
    end

    return (;
        basisfile = String(basisfile),
        hamfile = String(hamfile),
        final_dimension = basis_facts.final_dimension,
        basis_artifact_kind = basis_facts.basis_artifact_kind,
        ham_artifact_kind = ham_facts.ham_artifact_kind,
        final_coefficients_size = basis_facts.final_coefficients_size,
        final_gto_cross_overlap_size = basis_facts.final_gto_cross_overlap_size,
        gto_self_overlap_size = basis_facts.gto_self_overlap_size,
        gto_residual_overlap_size = basis_facts.gto_residual_overlap_size,
        h1_finite = ham_facts.h1_finite,
        h1_symmetry_error = ham_facts.h1_symmetry_error,
        h1_lowest = ham_facts.h1_lowest,
        h1_j_self_coulomb = ham_facts.h1_j_self_coulomb,
        h1_j_self_coulomb_positive = ham_facts.h1_j_self_coulomb_positive,
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

    sidecar = _pqs_source_box_route_driver_pqs_gto_sidecar(inputs)
    route_metadata =
        _pqs_source_box_route_driver_pqs_h2_route_metadata(report, inputs)
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
        provider_blocks_included = false,
    )

    basis_artifact_path = nothing
    basis_artifact_reloaded = false
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
        end
        basis_artifact_reloaded =
            _pqs_source_box_route_driver_h2_pqs_basis_artifact_reload_ok(
                basis_artifact_path,
            )
    end

    ham_artifact_path = nothing
    ham_artifact_reloaded = false
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
        end
        ham_artifact_reloaded =
            _pqs_source_box_route_driver_h2_pqs_ham_artifact_reload_ok(
                ham_artifact_path,
                final_basis.final_dimension,
            )
    end

    artifact_roundtrip =
        save_basis_artifact && save_ham_artifact ?
        pqs_h2_residual_gto_sidecar_artifact_roundtrip(
            basis_artifact_path,
            ham_artifact_path,
        ) :
        nothing

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
        provider_blocks_included = false,
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
