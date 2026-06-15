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

function _pqs_source_box_route_driver_pqs_gto_support_one_body(
    inputs,
    sidecar,
)
    source_plan = inputs.source_plan
    supplement = inputs.supplement_representation
    states = sidecar.support_states
    axes = sidecar.axis_representations
    coefficients = Float64.(coulomb_gaussian_expansion(doacc = false).coefficients)
    exponents = Float64.(coulomb_gaussian_expansion(doacc = false).exponents)
    route_metadata = inputs.route_metadata
    nuclear_charges = Float64.(route_metadata.nuclear_charges)
    atom_locations = route_metadata.atom_locations
    support_gto_kinetic = zeros(Float64, length(states), length(supplement.orbitals))
    support_gto_charged_nuclear =
        zeros(Float64, length(states), length(supplement.orbitals))
    scratch = zeros(Float64, length(states))

    for (column, orbital) in pairs(supplement.orbitals)
        orbital_arrays = _pqs_source_box_route_driver_gto_primitive_arrays(orbital)
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

    final_coefficients = Matrix{Float64}(inputs.final_basis.final_coefficients)
    h_fg_kinetic = Matrix{Float64}(transpose(final_coefficients) * support_gto_kinetic)
    h_fg_charged_nuclear =
        Matrix{Float64}(transpose(final_coefficients) * support_gto_charged_nuclear)
    return (;
        h_fg_kinetic,
        h_fg_charged_nuclear,
        h_fg_one_body = h_fg_kinetic + h_fg_charged_nuclear,
    )
end

function _pqs_source_box_route_driver_pqs_gto_self_one_body(inputs)
    supplement = inputs.supplement_representation
    expansion = coulomb_gaussian_expansion(doacc = false)
    coefficients = Float64.(expansion.coefficients)
    exponents = Float64.(expansion.exponents)
    route_metadata = inputs.route_metadata
    nuclear_charges = Float64.(route_metadata.nuclear_charges)
    atom_locations = route_metadata.atom_locations
    norbital = length(supplement.orbitals)
    h_gg_kinetic = zeros(Float64, norbital, norbital)
    h_gg_charged_nuclear = zeros(Float64, norbital, norbital)
    orbital_arrays = map(
        _pqs_source_box_route_driver_gto_primitive_arrays,
        supplement.orbitals,
    )
    for row in eachindex(orbital_arrays), column in eachindex(orbital_arrays)
        left_arrays = orbital_arrays[row]
        right_arrays = orbital_arrays[column]
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
    h_gg_kinetic = Matrix{Float64}(0.5 .* (h_gg_kinetic .+ transpose(h_gg_kinetic)))
    h_gg_charged_nuclear =
        Matrix{Float64}(0.5 .* (h_gg_charged_nuclear .+ transpose(h_gg_charged_nuclear)))
    return (;
        h_gg_kinetic,
        h_gg_charged_nuclear,
        h_gg_one_body = h_gg_kinetic + h_gg_charged_nuclear,
    )
end

function _pqs_source_box_route_driver_pqs_gto_one_body_blocks(
    inputs,
    sidecar,
    residual,
)
    h_ff = Matrix{Float64}(inputs.h1_hamiltonian.hamiltonian_matrix)
    s_fg = Matrix{Float64}(sidecar.final_gto_cross_overlap)
    l = Matrix{Float64}(residual.residual_transform)
    mixed = _pqs_source_box_route_driver_pqs_gto_support_one_body(inputs, sidecar)
    self = _pqs_source_box_route_driver_pqs_gto_self_one_body(inputs)
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
        provider_blocks_included = file["provider_blocks_included"]
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
            if provider_blocks_included === :one_body_only
                _pqs_source_box_route_driver_require_jld2_keys(
                    file,
                    (
                        :residual_transform,
                        :residual_rank,
                        :residual_overlap_eigenvalues,
                        :residual_overlap_identity_error,
                        :residual_overlap_cutoff,
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
                (;
                    residual_transform_size = size(residual_transform),
                    residual_rank,
                    residual_overlap_identity_error,
                    residual_overlap_cutoff = Float64(file["residual_overlap_cutoff"]),
                    provider_blocks_included,
                )
            else
                (;
                    residual_transform_size = nothing,
                    residual_rank = nothing,
                    residual_overlap_identity_error = nothing,
                    residual_overlap_cutoff = nothing,
                    provider_blocks_included,
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
        provider_blocks_included = file["provider_blocks_included"]
        provider_blocks_included == basis_facts.provider_blocks_included ||
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
            if provider_blocks_included === :one_body_only
                _pqs_source_box_route_driver_require_jld2_keys(
                    file,
                    (
                        :augmented_dimension,
                        :augmented_one_body_hamiltonian,
                        :augmented_h1_lowest,
                        :augmented_h1_symmetry_error,
                        :h_fg_kinetic,
                        :h_fg_charged_nuclear,
                        :h_gg_kinetic,
                        :h_gg_charged_nuclear,
                    ),
                )
                residual_rank = basis_facts.residual_rank
                augmented_dimension = Int(file["augmented_dimension"])
                augmented_dimension == final_dimension + residual_rank ||
                    throw(ArgumentError("augmented_dimension must equal final_dimension + residual_rank"))
                augmented_one_body_hamiltonian =
                    file["augmented_one_body_hamiltonian"]
                size(augmented_one_body_hamiltonian) ==
                    (augmented_dimension, augmented_dimension) ||
                    throw(ArgumentError("augmented one-body Hamiltonian shape mismatch"))
                all(isfinite, augmented_one_body_hamiltonian) ||
                    throw(ArgumentError("augmented one-body Hamiltonian contains non-finite entries"))
                augmented_h1_symmetry_error =
                    norm(
                        augmented_one_body_hamiltonian -
                        transpose(augmented_one_body_hamiltonian),
                        Inf,
                    )
                augmented_h1_symmetry_error <= Float64(symmetry_atol) ||
                    throw(ArgumentError("augmented one-body Hamiltonian is not symmetric"))
                augmented_h1_lowest = Float64(file["augmented_h1_lowest"])
                isfinite(augmented_h1_lowest) ||
                    throw(ArgumentError("augmented_h1_lowest must be finite"))
                size(file["h_fg_kinetic"]) == (final_dimension, size(gto_self_overlap, 1)) ||
                    throw(ArgumentError("h_fg_kinetic shape mismatch"))
                size(file["h_fg_charged_nuclear"]) ==
                    (final_dimension, size(gto_self_overlap, 1)) ||
                    throw(ArgumentError("h_fg_charged_nuclear shape mismatch"))
                _pqs_source_box_route_driver_square_matrix(
                    "h_gg_kinetic",
                    file["h_gg_kinetic"],
                )
                _pqs_source_box_route_driver_square_matrix(
                    "h_gg_charged_nuclear",
                    file["h_gg_charged_nuclear"],
                )
                (;
                    augmented_dimension,
                    augmented_one_body_hamiltonian_size =
                        size(augmented_one_body_hamiltonian),
                    augmented_h1_lowest,
                    augmented_h1_symmetry_error,
                )
            else
                (;
                    augmented_dimension = nothing,
                    augmented_one_body_hamiltonian_size = nothing,
                    augmented_h1_lowest = nothing,
                    augmented_h1_symmetry_error = nothing,
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
        augmented_dimension = ham_facts.augmented_dimension,
        augmented_one_body_hamiltonian_size =
            ham_facts.augmented_one_body_hamiltonian_size,
        augmented_h1_lowest = ham_facts.augmented_h1_lowest,
        augmented_h1_symmetry_error = ham_facts.augmented_h1_symmetry_error,
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

    residual_gto_provider_blocks in (:none, :one_body_only) || throw(
        ArgumentError(
            "residual_gto_provider_blocks must be :none or :one_body_only",
        ),
    )
    provider_blocks_included =
        residual_gto_provider_blocks === :one_body_only ? :one_body_only : false
    sidecar = _pqs_source_box_route_driver_pqs_gto_sidecar(inputs)
    route_metadata = _pqs_source_box_route_driver_pqs_h2_route_metadata(
        report,
        merge(inputs, (; provider_blocks_included,)),
    )
    inputs = merge(inputs, (; route_metadata, provider_blocks_included,))
    residual =
        residual_gto_provider_blocks === :one_body_only ?
        _pqs_source_box_route_driver_pqs_gto_residual_transform(sidecar) :
        nothing
    one_body_blocks =
        residual_gto_provider_blocks === :one_body_only ?
        _pqs_source_box_route_driver_pqs_gto_one_body_blocks(
            inputs,
            sidecar,
            residual,
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
