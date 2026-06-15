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
    payload = report.route_materializer_payload
    basis =
        hasproperty(payload, :parent_basis_object) ? payload.parent_basis_object : nothing
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
        object_kind = :white_lindsey_atomic_pure_gausslet_materialization_summary,
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
        one_body_hamiltonian_available = true,
        density_density_matrix_available = hasproperty(fixed_block, :pair_sum),
        jld2_artifact_status =
            (save_basis_artifact || save_ham_artifact) ?
            :written_if_requested :
            :not_requested,
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
        object_kind = :cartesian_nesting_route_driver_materialization,
        route_family = report.route_family,
        status = :materialized_white_lindsey_atomic_pure_gausslet,
        blocker = nothing,
        materialized_report = summary,
        materialized_report_kind = summary.object_kind,
        retained_dimension = dim,
        support_dimension = summary.support_dimension,
        final_dimension = dim,
        h1_lowest,
        h1_finite = summary.h1_finite,
        h1_symmetry_error,
        overlap_identity_error,
        pqs_materialization_status = :not_applicable,
        save_basis_artifact_requested = save_basis_artifact,
        save_ham_artifact_requested = save_ham_artifact,
        basisfile = basis_artifact_path,
        hamfile = ham_artifact_path,
        jld2_basis_artifact_status =
            save_basis_artifact ? :written : :not_requested,
        jld2_ham_artifact_status =
            save_ham_artifact ? :written : :not_requested,
        ignored_legacy_keyword_count = 0,
    )
end

function _pqs_source_box_route_driver_materialization(
    report;
    materialize_route::Bool = false,
    save_basis_artifact::Bool = false,
    save_ham_artifact::Bool = false,
    basisfile = nothing,
    hamfile = nothing,
    kwargs...,
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
                materializer_backend = get(kwargs, :materializer_backend, nothing),
                materializer_nside = get(kwargs, :materializer_nside, nothing),
                white_lindsey_expansion =
                    get(kwargs, :white_lindsey_expansion, nothing),
                white_lindsey_Z = get(kwargs, :white_lindsey_Z, nothing),
            )
        !isnothing(wl_materialization) && return wl_materialization
    end
    return (;
        object_kind = :cartesian_nesting_route_driver_materialization,
        route_family = hasproperty(report, :route_family) ? report.route_family : nothing,
        status =
            requested ?
            :blocked_materialization_after_route_scaffold_demolition :
            :not_requested,
        blocker =
            requested ?
            :route_configured_low_order_materializer_removed :
            nothing,
        materialized_report = nothing,
        materialized_report_kind = nothing,
        retained_dimension,
        pqs_materialization_status =
            requested ?
            :blocked_materialization_after_route_scaffold_demolition :
            :not_requested,
        save_basis_artifact_requested = save_basis_artifact,
        save_ham_artifact_requested = save_ham_artifact,
        basisfile,
        hamfile,
        ignored_legacy_keyword_count = length(kwargs),
    )
end
