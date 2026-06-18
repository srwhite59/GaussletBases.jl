# Private H2 PQS residual-GTO sidecar, provider-block, and Ham handoff helpers.

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
    diagnostics = (;
        sidecar_kind = :pqs_h2_residual_gto_sidecar,
        gto_residual_overlap_symmetry_error = residual_symmetry_error,
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
        gto_primitive_arrays,
        final_coefficients = Matrix{Float64}(inputs.final_basis.final_coefficients),
        support_gto_cross = Matrix{Float64}(sidecar.support_gto_cross),
        h_ff = Matrix{Float64}(inputs.h1_hamiltonian.hamiltonian_matrix),
        density_interaction = inputs.density_interaction,
        s_fg = Matrix{Float64}(sidecar.final_gto_cross_overlap),
        s_gg = Matrix{Float64}(sidecar.gto_self_overlap),
        residual_transform = Matrix{Float64}(residual.residual_transform),
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
        residual_overlap_error,
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
    v_dd = [
        Matrix{Float64}(density_interaction.pre_final_pair_matrix) v_pr
        transpose(v_pr) v_rr
    ]
    v_dd_symmetry_error = norm(v_dd - transpose(v_dd), Inf)
    return (;
        augmented_density_gauge = :pre_final_localized_positive_weight,
        augmented_density_space = (:pre_final_pqs, :residual_gto),
        augmented_pair_matrix = v_dd,
        augmented_density_dimension = size(v_dd, 1),
        augmented_pair_matrix_size = size(v_dd),
        augmented_pair_matrix_symmetry_error = v_dd_symmetry_error,
        residual_centers = moments.residual_centers,
        residual_widths = moments.residual_widths,
        residual_overlap_error = moments.residual_overlap_error,
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

function _pqs_source_box_route_driver_pqs_h2_residual_gto_ham_handoff(
    one_body_blocks,
    density_blocks,
    density_interaction,
    nuclear_repulsion,
)
    isnothing(one_body_blocks) &&
        throw(ArgumentError("H2 residual-GTO Ham handoff requires one-body blocks"))
    isnothing(density_blocks) &&
        throw(ArgumentError("H2 residual-GTO Ham handoff requires density blocks"))
    isfinite(Float64(nuclear_repulsion)) ||
        throw(ArgumentError("H2 residual-GTO Ham handoff requires finite nuclear repulsion"))
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
        orbital_basis = (:final_pqs, :residual_gto),
        density_basis = (:pre_final_pqs, :residual_gto),
        one_body_hamiltonian = one_body,
        density_pair_matrix = pair_matrix,
        orbital_to_density,
        electron_count = 2,
        spin_sectors = (; nup = 1, ndn = 1),
        nuclear_repulsion = Float64(nuclear_repulsion),
        diagnostics = (;
            orbital_dimension,
            density_dimension,
            residual_rank = density_blocks.augmented_density_dimension -
                size(density_interaction.final_to_pre_final_coefficients, 1),
            one_body_symmetry_error,
            pair_symmetry_error,
        ),
    )
end

function _pqs_source_box_route_driver_ham_handoff_summary(ham_handoff)
    isnothing(ham_handoff) &&
        return (;
            ham_handoff_orbital_basis = nothing,
            ham_handoff_density_basis = nothing,
            ham_handoff_orbital_dimension = nothing,
            ham_handoff_density_dimension = nothing,
        )

    return (;
        ham_handoff_orbital_basis = ham_handoff.orbital_basis,
        ham_handoff_density_basis = ham_handoff.density_basis,
        ham_handoff_orbital_dimension = ham_handoff.diagnostics.orbital_dimension,
        ham_handoff_density_dimension = ham_handoff.diagnostics.density_dimension,
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
    return h_fg_kinetic + h_fg_charged_nuclear
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
    return h_gg_kinetic + h_gg_charged_nuclear
end

function _pqs_source_box_route_driver_pqs_gto_one_body_blocks(packet)
    h_ff = packet.h_ff
    s_fg = packet.s_fg
    l = packet.residual_transform
    h_fg = Matrix{Float64}(
        _pqs_source_box_route_driver_pqs_gto_support_one_body(packet),
    )
    h_gg = Matrix{Float64}(
        _pqs_source_box_route_driver_pqs_gto_self_one_body(packet),
    )
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
    provider_block_mode = _pqs_source_box_route_driver_pqs_h2_provider_block_mode(
        get(inputs, :provider_block_mode, false),
    )
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
        provider_block_mode,
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
    mode === false && return :none
    mode === :none && return :none
    mode === :one_body_only && return :one_body_only
    mode === :one_body_and_density_provider && return :one_body_and_density_provider
    throw(ArgumentError("provider block mode must be :none, :one_body_only, or :one_body_and_density_provider; got $(mode)"))
end

function _pqs_source_box_route_driver_read_pqs_h2_provider_block_mode(file)
    haskey(file, "provider_block_mode") ||
        throw(ArgumentError("H2 PQS sidecar artifact missing provider_block_mode"))
    return _pqs_source_box_route_driver_pqs_h2_provider_block_mode(
        file["provider_block_mode"],
    )
end

function _pqs_source_box_route_driver_optional_property(source, property::Symbol)
    return isnothing(source) ? nothing : getproperty(source, property)
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

function _pqs_source_box_route_driver_ham_handoff_consumer_invariant(
    one_body_hamiltonian,
    density_pair_matrix,
    orbital_to_density,
)
    h1_orbital = eigen(Symmetric(one_body_hamiltonian)).vectors[:, 1]
    density_coefficients = orbital_to_density * h1_orbital
    self_coulomb =
        _pqs_source_box_route_driver_restricted_one_orbital_self_coulomb(
            density_pair_matrix,
            density_coefficients,
        )
    isfinite(self_coulomb) && self_coulomb > 0 ||
        throw(ArgumentError("Ham handoff consumer H/V/T scalar must be finite and positive"))
    return (; self_coulomb,)
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
                :provider_block_mode,
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
        provider_block_mode =
            _pqs_source_box_route_driver_read_pqs_h2_provider_block_mode(file)
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
            if provider_block_mode !== :none
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
                    if provider_block_mode === :one_body_and_density_provider
                        _pqs_source_box_route_driver_require_jld2_keys(
                            file,
                            (
                                :residual_centers,
                                :residual_widths,
                                :residual_moment_overlap_error,
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
                        (;
                            residual_centers_size = size(residual_centers),
                            residual_widths_size = size(residual_widths),
                            residual_moment_overlap_error,
                        )
                    else
                        (;
                            residual_centers_size = nothing,
                            residual_widths_size = nothing,
                            residual_moment_overlap_error = nothing,
                        )
                    end
                (;
                    residual_transform_size = size(residual_transform),
                    residual_rank,
                    residual_overlap_identity_error,
                    residual_overlap_cutoff = Float64(file["residual_overlap_cutoff"]),
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
            provider_block_mode,
            final_coefficients_size = size(final_coefficients),
            final_gto_cross_overlap_size = size(final_gto_cross_overlap),
            gto_self_overlap_size = size(gto_self_overlap),
            gto_residual_overlap_size = size(gto_residual_overlap),
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
                :final_gto_cross_overlap,
                :gto_self_overlap,
                :gto_residual_overlap,
                :provider_block_mode,
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
        final_gto_cross_overlap = file["final_gto_cross_overlap"]
        gto_self_overlap = file["gto_self_overlap"]
        gto_residual_overlap = file["gto_residual_overlap"]
        provider_block_mode =
            _pqs_source_box_route_driver_read_pqs_h2_provider_block_mode(file)
        provider_block_mode === basis_facts.provider_block_mode ||
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
            if provider_block_mode !== :none
                _pqs_source_box_route_driver_require_jld2_keys(
                    file,
                    (
                        :augmented_dimension,
                        :augmented_one_body_hamiltonian,
                        :augmented_h1_lowest,
                        :augmented_h1_symmetry_error,
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
                    if provider_block_mode === :one_body_and_density_provider
                        _pqs_source_box_route_driver_require_jld2_keys(
                            file,
                            (
                                :augmented_pair_matrix,
                                :augmented_density_dimension,
                                :augmented_pair_matrix_symmetry_error,
                                :residual_centers,
                                :residual_widths,
                                :residual_moment_overlap_error,
                                :ham_handoff_orbital_basis,
                                :ham_handoff_density_basis,
                                :ham_handoff_orbital_to_density,
                                :ham_handoff_electron_count,
                                :ham_handoff_spin_nup,
                                :ham_handoff_spin_ndn,
                                :ham_handoff_nuclear_repulsion,
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
                        augmented_pair_matrix = file["augmented_pair_matrix"]
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
                        Tuple(file["ham_handoff_orbital_basis"]) ===
                            (:final_pqs, :residual_gto) ||
                            throw(ArgumentError("unexpected Ham handoff orbital basis"))
                        Tuple(file["ham_handoff_density_basis"]) ===
                            (:pre_final_pqs, :residual_gto) ||
                            throw(ArgumentError("unexpected Ham handoff density basis"))
                        orbital_to_density = file["ham_handoff_orbital_to_density"]
                        _pqs_source_box_route_driver_finite_matrix(
                            "ham_handoff_orbital_to_density",
                            orbital_to_density,
                            (augmented_density_dimension, augmented_dimension),
                        )
                        Int(file["ham_handoff_electron_count"]) == 2 ||
                            throw(ArgumentError("Ham handoff electron count must be 2"))
                        Int(file["ham_handoff_spin_nup"]) == 1 ||
                            throw(ArgumentError("Ham handoff spin nup must be 1"))
                        Int(file["ham_handoff_spin_ndn"]) == 1 ||
                            throw(ArgumentError("Ham handoff spin ndn must be 1"))
                        isfinite(Float64(file["ham_handoff_nuclear_repulsion"])) ||
                            throw(ArgumentError("Ham handoff nuclear repulsion must be finite"))
                        consumer_invariant =
                            _pqs_source_box_route_driver_ham_handoff_consumer_invariant(
                                augmented_one_body_hamiltonian,
                                augmented_pair_matrix,
                                orbital_to_density,
                            )
                        (;
                            augmented_density_dimension,
                            augmented_pair_matrix_size = size(augmented_pair_matrix),
                            augmented_pair_matrix_symmetry_error,
                            ham_handoff_orbital_to_density_size =
                                size(orbital_to_density),
                            ham_handoff_consumer_self_coulomb =
                                consumer_invariant.self_coulomb,
                        )
                    else
                        (;
                            augmented_density_dimension = nothing,
                            augmented_pair_matrix_size = nothing,
                            augmented_pair_matrix_symmetry_error = nothing,
                            ham_handoff_orbital_to_density_size = nothing,
                            ham_handoff_consumer_self_coulomb = nothing,
                        )
                    end
                (;
                    augmented_dimension,
                    augmented_one_body_hamiltonian_size =
                        size(augmented_one_body_hamiltonian),
                    augmented_h1_lowest,
                    augmented_h1_symmetry_error,
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
                    nuclear_mixed_block_convention = nothing,
                    augmented_density_gauge = nothing,
                    augmented_density_space = nothing,
                    p_projection_of_g_size = nothing,
                    residual_orbital_coefficients_in_density_carrier_size = nothing,
                    ham_handoff_orbital_to_density_size = nothing,
                )
            end
        return (;
            ham_artifact_kind = artifact_kind,
            provider_block_mode,
            one_body_hamiltonian_size = size(one_body_hamiltonian),
            h1_symmetry_error,
            h1_lowest,
            pre_final_pair_matrix_size = size(pre_final_pair_matrix),
            final_gto_cross_overlap_size = size(final_gto_cross_overlap),
            gto_self_overlap_size = size(gto_self_overlap),
            gto_residual_overlap_size = size(gto_residual_overlap),
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
        provider_block_mode = basis_facts.provider_block_mode,
        final_coefficients_size = basis_facts.final_coefficients_size,
        final_gto_cross_overlap_size = basis_facts.final_gto_cross_overlap_size,
        gto_self_overlap_size = basis_facts.gto_self_overlap_size,
        gto_residual_overlap_size = basis_facts.gto_residual_overlap_size,
        h1_symmetry_error = ham_facts.h1_symmetry_error,
        h1_lowest = ham_facts.h1_lowest,
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
        ham_handoff_orbital_to_density_size =
            ham_facts.ham_handoff_orbital_to_density_size,
        ham_handoff_consumer_self_coulomb =
            ham_facts.ham_handoff_consumer_self_coulomb,
        basis_gto_residual_overlap_symmetry_error =
            basis_facts.gto_residual_overlap_symmetry_error,
        ham_gto_residual_overlap_symmetry_error =
            ham_facts.gto_residual_overlap_symmetry_error,
    )
end

function _pqs_source_box_route_driver_write_pqs_h2_sidecar_common!(
    file,
    inputs,
    sidecar,
    provider_block_mode,
)
    file["gto_supplement_metadata"] = inputs.supplement_representation.metadata
    file["final_gto_cross_overlap"] = sidecar.final_gto_cross_overlap
    file["gto_self_overlap"] = sidecar.gto_self_overlap
    file["gto_residual_overlap"] = sidecar.gto_residual_overlap
    file["gto_sidecar_diagnostics"] = sidecar.diagnostics
    file["provider_block_mode"] = provider_block_mode
    return nothing
end

function _pqs_source_box_route_driver_write_jld2_fields!(file, source, fields)
    isnothing(source) && return nothing
    for field in fields
        key = field isa Pair ? first(field) : field
        property = field isa Pair ? last(field) : field
        file[String(key)] = getproperty(source, property)
    end
    return nothing
end

const _PQS_H2_RESIDUAL_TRANSFORM_ARTIFACT_FIELDS = (
    :residual_transform,
    :residual_rank,
    :residual_overlap_eigenvalues,
    :residual_overlap_identity_error,
    :residual_overlap_cutoff,
)

const _PQS_H2_DENSITY_MOMENT_ARTIFACT_FIELDS = (
    :residual_centers,
    :residual_widths,
    :residual_moment_overlap_error => :residual_overlap_error,
)
