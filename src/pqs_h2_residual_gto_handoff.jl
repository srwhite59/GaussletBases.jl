# Private H2 PQS residual-GTO provider-block and Ham artifact helpers.

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
        kinetic_ff = Matrix{Float64}(inputs.final_kinetic.final_operator),
        nuclear_unit_by_center_ff =
            Matrix{Float64}[Matrix{Float64}(record.final_operator)
                            for record in inputs.final_nuclear_by_center],
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
    s_fg = Matrix{Float64}(sidecar.final_gto_cross_overlap)
    residual_transform = Matrix{Float64}(residual.residual_transform)
    f_dimension = final_basis.final_dimension
    p_dimension = f_dimension
    size(s_fg, 1) == f_dimension ||
        throw(DimensionMismatch("density descriptor S_FG row count must match final dimension"))
    g_dimension = size(s_fg, 2)
    size(residual_transform, 1) == g_dimension ||
        throw(DimensionMismatch("density descriptor residual transform row count must match GTO dimension"))
    residual_rank = size(residual_transform, 2)
    residual_rank == residual.residual_rank ||
        throw(DimensionMismatch("density descriptor residual rank mismatch"))
    p_projection_of_g = s_fg
    residual_carrier =
        Matrix{Float64}(
            vcat(
                -p_projection_of_g,
                Matrix{Float64}(I, g_dimension, g_dimension),
            ) * residual_transform,
        )
    pair_matrix = Matrix{Float64}(density_interaction.electron_electron_ida)
    size(pair_matrix) == (p_dimension, p_dimension) ||
        throw(DimensionMismatch("density descriptor P-P pair matrix dimension mismatch"))
    get(density_interaction, :density_gauge, nothing) === :localized_ida ||
        throw(ArgumentError("density descriptor requires the localized IDA gauge"))
    return (;
        density_gauge = :localized_ida,
        augmented_density_space = (:localized_ida_pqs, :residual_gto),
        p_dimension,
        f_dimension,
        g_dimension,
        residual_rank,
        p_projection_of_g,
        residual_orbital_coefficients_in_density_carrier = residual_carrier,
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
                _cartesian_weighted_hadamard3(
                    left_arrays.coefficients,
                    right_arrays.coefficients,
                    axis_tables[1],
                    axis_tables[2],
                    axis_tables[3],
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
    final_coefficients = Matrix{Float64}(final_basis.final_coefficients)
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
    pp = Matrix{Float64}(transpose(final_coefficients) * support_pp * final_coefficients)
    pg = Matrix{Float64}(transpose(final_coefficients) * support_pg)
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
    final_coefficients =
        Matrix{Float64}(packet.final_basis.final_coefficients)
    support_overlap =
        _pqs_source_box_route_driver_pqs_support_overlap_matrix(packet)
    pp_overlap =
        Matrix{Float64}(
            transpose(final_coefficients) *
            support_overlap *
            final_coefficients,
        )
    pg_overlap =
        Matrix{Float64}(transpose(final_coefficients) * packet.support_gto_cross)
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
    final_coefficients = Matrix{Float64}(packet.final_basis.final_coefficients)
    ida_weights = Float64.(density_interaction.ida_weights)
    weighted_coefficients =
        final_coefficients .* reshape(1.0 ./ ida_weights, 1, :)
    v_pr =
        Matrix{Float64}(
            transpose(weighted_coefficients) * support_residual_raw_numerator,
        )
    v_rr = Matrix{Float64}(components.residual_residual)
    v_rr = Matrix{Float64}(0.5 .* (v_rr .+ transpose(v_rr)))
    v_dd = [
        Matrix{Float64}(density_interaction.electron_electron_ida) v_pr
        transpose(v_pr) v_rr
    ]
    return (;
        augmented_density_gauge = :localized_ida,
        augmented_density_space = (:localized_ida_pqs, :residual_gto),
        augmented_pair_matrix = v_dd,
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

function _pqs_source_box_route_driver_pqs_h2_residual_gto_ida_hamiltonian(
    one_body_blocks,
    density_blocks,
    route_metadata,
)
    isnothing(one_body_blocks) &&
        throw(ArgumentError("H2 residual-GTO IDA Hamiltonian requires one-body blocks"))
    isnothing(density_blocks) &&
        throw(ArgumentError("H2 residual-GTO IDA Hamiltonian requires density blocks"))
    isnothing(route_metadata.nup) &&
        throw(ArgumentError("H2 residual-GTO IDA Hamiltonian requires explicit nup system input"))
    isnothing(route_metadata.ndn) &&
        throw(ArgumentError("H2 residual-GTO IDA Hamiltonian requires explicit ndn system input"))
    return CartesianIDAHamiltonian(
        one_body_blocks.augmented_kinetic,
        one_body_blocks.augmented_nuclear_attraction_unit_by_center,
        density_blocks.augmented_pair_matrix,
        route_metadata.nup,
        route_metadata.ndn;
        nuclear_charges = route_metadata.nuclear_charges,
        nuclear_positions = route_metadata.atom_locations,
    )
end

function _pqs_source_box_route_driver_ida_hamiltonian_summary(ida_hamiltonian)
    isnothing(ida_hamiltonian) &&
        return (;
            ida_orbital_dimension = nothing,
            ida_center_count = nothing,
        )

    return (;
        ida_orbital_dimension = size(ida_hamiltonian.kinetic, 1),
        ida_center_count = length(ida_hamiltonian.nuclear_charges),
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
    right_centers = @view orbital_arrays.centers[:, axis_index]
    right_powers = @view orbital_arrays.angular_powers[:, axis_index]
    primitive_cross = _cartesian_gaussian_axis_integral_table(
        Float64[
            GaussianAnalyticIntegrals.gaussian_exponent(primitive)
            for primitive in basis_primitives
        ],
        Float64[primitive.center_value for primitive in basis_primitives],
        zeros(Int, length(basis_primitives)),
        ones(Float64, length(basis_primitives)),
        orbital_arrays.exponents,
        right_centers,
        right_powers,
        _cartesian_gaussian_axis_prefactors(orbital_arrays.exponents, right_powers),
        term;
        factor_exponent = isnothing(factor_exponent) ? 0.0 : Float64(factor_exponent),
        factor_center = isnothing(factor_center) ? 0.0 : Float64(factor_center),
    )
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
    return _cartesian_gaussian_axis_integral_table(
        left_arrays.exponents,
        center_left_axis,
        power_left_axis,
        _cartesian_gaussian_axis_prefactors(left_arrays.exponents, power_left_axis),
        right_arrays.exponents,
        center_right_axis,
        power_right_axis,
        _cartesian_gaussian_axis_prefactors(right_arrays.exponents, power_right_axis),
        term;
        factor_exponent = isnothing(factor_exponent) ? 0.0 : Float64(factor_exponent),
        factor_center = isnothing(factor_center) ? 0.0 : Float64(factor_center),
    )
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

function _pqs_source_box_route_driver_pqs_gto_support_one_body(packet)
    states = packet.support_states
    axes = packet.axis_representations
    coefficients = packet.coulomb_coefficients
    exponents = packet.coulomb_exponents
    atom_locations = packet.atom_locations
    support_gto_kinetic =
        zeros(Float64, length(states), length(packet.gto_primitive_arrays))
    support_gto_nuclear_by_center =
        [zeros(Float64, length(states), length(packet.gto_primitive_arrays))
         for _ in atom_locations]
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
                support_gto_nuclear_by_center[center_index][:, column] .-=
                    coefficients[term_index] .* scratch
            end
        end
    end

    kinetic_fg =
        Matrix{Float64}(transpose(packet.final_coefficients) * support_gto_kinetic)
    nuclear_unit_by_center_fg =
        Matrix{Float64}[
            Matrix{Float64}(transpose(packet.final_coefficients) * nuclear)
            for nuclear in support_gto_nuclear_by_center
        ]
    return (; kinetic_fg, nuclear_unit_by_center_fg)
end

function _pqs_source_box_route_driver_pqs_gto_self_one_body(packet)
    coefficients = packet.coulomb_coefficients
    exponents = packet.coulomb_exponents
    atom_locations = packet.atom_locations
    norbital = length(packet.gto_primitive_arrays)
    kinetic_gg = zeros(Float64, norbital, norbital)
    nuclear_unit_by_center_gg =
        [zeros(Float64, norbital, norbital) for _ in atom_locations]
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
                kinetic_gg[row, column] +=
                    _cartesian_weighted_hadamard3(
                        left_arrays.coefficients,
                        right_arrays.coefficients,
                        axis_tables[1],
                        axis_tables[2],
                        axis_tables[3],
                    )
            end
            for (center_index, location) in pairs(atom_locations)
                for term_index in eachindex(coefficients)
                    factor_tables =
                        _pqs_source_box_route_driver_gto_axis_self_tables(
                            left_arrays,
                            right_arrays,
                            :factor;
                            factor_center = location,
                            factor_exponent = exponents[term_index],
                        )
                    nuclear_unit_by_center_gg[center_index][row, column] -=
                        coefficients[term_index] *
                        _cartesian_weighted_hadamard3(
                            left_arrays.coefficients,
                            right_arrays.coefficients,
                            factor_tables[1],
                            factor_tables[2],
                            factor_tables[3],
                        )
                end
            end
        end
    end
    kinetic_gg = Matrix{Float64}(0.5 .* (kinetic_gg .+ transpose(kinetic_gg)))
    nuclear_unit_by_center_gg =
        Matrix{Float64}[
            Matrix{Float64}(0.5 .* (nuclear .+ transpose(nuclear)))
            for nuclear in nuclear_unit_by_center_gg
        ]
    return (; kinetic_gg, nuclear_unit_by_center_gg)
end

function _pqs_source_box_route_driver_pqs_gto_augmented_one_body(
    o_ff,
    o_fg,
    o_gg,
    s_fg,
    residual_transform,
)
    ff = Matrix{Float64}(o_ff)
    fg = Matrix{Float64}(o_fg)
    gg = Matrix{Float64}(o_gg)
    s = Matrix{Float64}(s_fg)
    l = Matrix{Float64}(residual_transform)
    o_fr = Matrix{Float64}((fg - ff * s) * l)
    residual_core = Matrix{Float64}(
        gg -
        transpose(s) * fg -
        transpose(fg) * s +
        transpose(s) * ff * s,
    )
    o_rr = Matrix{Float64}(transpose(l) * residual_core * l)
    o_rr = Matrix{Float64}(0.5 .* (o_rr .+ transpose(o_rr)))
    augmented = [
        ff o_fr
        transpose(o_fr) o_rr
    ]
    return Matrix{Float64}(0.5 .* (augmented .+ transpose(augmented)))
end

function _pqs_source_box_route_driver_pqs_gto_one_body_blocks(packet)
    s_fg = packet.s_fg
    l = packet.residual_transform
    support = _pqs_source_box_route_driver_pqs_gto_support_one_body(packet)
    self = _pqs_source_box_route_driver_pqs_gto_self_one_body(packet)
    augmented_kinetic =
        _pqs_source_box_route_driver_pqs_gto_augmented_one_body(
            packet.kinetic_ff,
            support.kinetic_fg,
            self.kinetic_gg,
            s_fg,
            l,
        )
    center_count = length(packet.nuclear_charges)
    length(packet.nuclear_unit_by_center_ff) == center_count ||
        throw(DimensionMismatch("PQS/GTO one-body blocks require one PQS unit-nuclear block per center"))
    length(support.nuclear_unit_by_center_fg) == center_count ||
        throw(DimensionMismatch("PQS/GTO one-body blocks require one support/GTO unit-nuclear block per center"))
    length(self.nuclear_unit_by_center_gg) == center_count ||
        throw(DimensionMismatch("PQS/GTO one-body blocks require one GTO/GTO unit-nuclear block per center"))
    augmented_nuclear_attraction_unit_by_center =
        Matrix{Float64}[
            _pqs_source_box_route_driver_pqs_gto_augmented_one_body(
                packet.nuclear_unit_by_center_ff[index],
                support.nuclear_unit_by_center_fg[index],
                self.nuclear_unit_by_center_gg[index],
                s_fg,
                l,
            ) for index in 1:center_count
        ]
    base_reconstructed =
        Matrix{Float64}(packet.kinetic_ff)
    for (charge, nuclear) in zip(packet.nuclear_charges, packet.nuclear_unit_by_center_ff)
        base_reconstructed .+= Float64(charge) .* Matrix{Float64}(nuclear)
    end
    base_reconstruction_error = norm(base_reconstructed - packet.h_ff, Inf)
    base_reconstruction_error <= 1.0e-8 ||
        throw(ArgumentError("PQS one-body H1 reconstruction from K + Z*U centers failed"))
    return (;
        augmented_kinetic,
        augmented_nuclear_attraction_unit_by_center,
        final_dimension = size(packet.h_ff, 1),
        augmented_dimension = size(augmented_kinetic, 1),
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
        nup = system.nup,
        ndn = system.ndn,
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
        artifact_scope = :pqs_h2_residual_gto_ida_hamiltonian,
        provider_block_mode,
    )
end

function _pqs_source_box_route_driver_pqs_h2_provider_block_mode(mode)
    mode === false && return :none
    mode === :none && return :none
    mode === :one_body_and_density_provider && return :one_body_and_density_provider
    throw(ArgumentError("provider block mode must be :none or :one_body_and_density_provider; got $(mode)"))
end

function _pqs_source_box_route_driver_optional_property(source, property::Symbol)
    return isnothing(source) ? nothing : getproperty(source, property)
end

function _pqs_source_box_route_driver_ida_hamiltonian_smoke(ida_hamiltonian)
    full_h1 = one_body_hamiltonian(ida_hamiltonian)
    full_h1_symmetry_error =
        norm(full_h1 - transpose(full_h1), Inf)
    full_h1_symmetry_error <= 1.0e-8 ||
        throw(ArgumentError("IDA full H1 reconstruction is not symmetric"))
    h1_eigen = eigen(Symmetric(full_h1))
    h1_orbital = h1_eigen.vectors[:, 1]
    self_coulomb =
        _pqs_source_box_route_driver_restricted_one_orbital_self_coulomb(
            ida_hamiltonian.electron_electron_ida,
            h1_orbital,
        )
    isfinite(self_coulomb) && self_coulomb > 0 ||
        throw(ArgumentError("IDA self-Coulomb scalar must be finite and positive"))

    center_count = length(ida_hamiltonian.nuclear_charges)
    center_count == 2 ||
        throw(ArgumentError("H2 IDA counterpoise smoke expects exactly two centers"))
    for active_center in 1:center_count
        center_weights = zeros(Float64, center_count)
        center_weights[active_center] = 1.0
        branch_h1 =
            one_body_hamiltonian(ida_hamiltonian; center_weights)
        norm(branch_h1 - transpose(branch_h1), Inf) <= 1.0e-8 ||
            throw(ArgumentError("IDA counterpoise branch H1 is not symmetric"))
        isfinite(eigen(Symmetric(branch_h1)).values[1]) ||
            throw(ArgumentError("IDA counterpoise branch H1 lowest value is not finite"))
        nuclear_repulsion(
            ida_hamiltonian;
            center_weights,
        ) == 0.0 ||
            throw(ArgumentError("IDA ghost branch nuclear repulsion must be zero"))
    end
    full_nuclear_repulsion = nuclear_repulsion(ida_hamiltonian)
    abs(full_nuclear_repulsion - ida_hamiltonian.nuclear_repulsion) <= 1.0e-8 ||
        throw(ArgumentError("IDA nuclear repulsion reconstruction mismatch"))
    return (;
        full_self_coulomb = self_coulomb,
        counterpoise_branch_count = center_count,
    )
end
