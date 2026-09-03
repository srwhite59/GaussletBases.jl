struct _SlicedLongitudinal
    centers::Vector{Float64}
    lattice_indices::Vector{Int}
    class_indices::Vector{Int}
    spacing::Float64
    coefficients::Vector{Float64}
    correlation::Vector{Float64}
    product_terms::Matrix{Float64}
    product_bounds::Vector{Float64}
    integral_weight::Float64
    ida_moments::NTuple{4,Float64}
    prototype_errors::NTuple{3,Float64}
end

struct _SlicedTransverse
    exponents::Vector{Float64}
    coefficients::Matrix{Float64}
    overlap::Matrix{Float64}
    kinetic::Matrix{Float64}
    radial_moments::Array{Float64,3}
    factor_terms::Array{Float64,3}
    pair_exponents::Vector{Float64}
    density_coefficients::Matrix{Float64}
    ranks::Vector{Int}
    gaps::Vector{Float64}
    phase_pivots::Vector{Int}
end

struct _SlicedCoulomb
    expansion::CoulombGaussianExpansion
    longitudinal_terms::Matrix{Float64}
    transverse_terms::Array{Float64,3}
    periodic_blocks::Array{Float64,3}
    h1_band::Int
    h1_offband_bound::Float64
    vee_tail_distance::Float64
    vee_tail_bound::Float64
    vee_transition_error::Float64
    h1_tail_distance::Float64
    h1_tail_bound::Float64
    h1_transition_error::Float64
end

struct _SlicedDiagnostics
    basis_fingerprint::String
    longitudinal_overlap_error::Float64
    transverse_normalization_error::Float64
    periodic_class_error::Float64
    omitted_projected_norm::Float64
    source_normalization_error::Float64
    actual_left_padding::Float64
    actual_right_padding::Float64
    maximum_separation::Float64
    retained_bytes::Int
end

"""Compact, in-memory minimal sliced hydrogen-chain operator owner."""
struct SlicedHydrogenChain
    nuclei::Vector{Float64}
    atom_spacing::Float64
    sites_per_atom::Int
    mode::Symbol
    longitudinal::_SlicedLongitudinal
    transverse::_SlicedTransverse
    coulomb::_SlicedCoulomb
    h1_bands::Matrix{Float64}
    nuclear_repulsion::Float64
    one_body_tolerance::Float64
    interaction_tolerance::Float64
    diagnostics::_SlicedDiagnostics
end

struct _SlicedH1View <: AbstractMatrix{Float64}
    chain::SlicedHydrogenChain
end

struct _SlicedVeeView <: AbstractMatrix{Float64}
    chain::SlicedHydrogenChain
end

@inline function _sliced_even_moment(x::Float64, variance::Float64, order::Int)
    order == 0 && return 1.0
    order == 2 && return x^2 + variance
    order == 4 && return x^4 + 6x^2 * variance + 3variance^2
    order == 6 && return x^6 + 15x^4 * variance + 45x^2 * variance^2 + 15variance^3
    order == 8 && return x^8 + 28x^6 * variance + 210x^4 * variance^2 +
        420x^2 * variance^3 + 105variance^4
    throw(ArgumentError("unsupported sliced even moment"))
end

function _sliced_longitudinal(spacing::Float64, centers::Vector{Float64},
    lattice::Vector{Int}, classes::Vector{Int}, tolerance::Float64)
    prototype = build_basis(UniformBasisSpec(:G10; xmin = -2spacing, xmax = 2spacing, spacing))
    representation = basis_representation(prototype; operators = (:overlap, :kinetic))
    midpoint = 3
    local_stencil = stencil(prototype[midpoint])
    coefficients_local = Float64.(coefficients(local_stencil))
    primitive_width = spacing / 3
    ncoeff = length(coefficients_local)
    stencil_error = 0.0
    for column in 1:length(prototype)
        translated = stencil(prototype[column])
        translated_coefficients = coefficients(translated)
        length(translated_coefficients) == ncoeff ||
            throw(ArgumentError("bounded G10 donor stencil sizes disagree"))
        for k in eachindex(translated_coefficients)
            expected_center = center(prototype[column]) + (k - (ncoeff + 1) / 2) * primitive_width
            stencil_error = max(stencil_error, abs(translated_coefficients[k] - coefficients_local[k]),
                abs(center(primitives(translated)[k]) - expected_center))
        end
    end
    stencil_error <= tolerance || throw(ArgumentError("bounded G10 translated stencil parity failed"))
    correlation = zeros(Float64, 2ncoeff - 1)
    @inbounds for a in 1:ncoeff, b in 1:ncoeff
        correlation[a - b + ncoeff] += coefficients_local[a] * coefficients_local[b]
    end
    pair_value(delta, kinetic) = begin
        value_out = 0.0
        @inbounds for k in eachindex(correlation)
            distance = delta + (k - ncoeff) * primitive_width
            overlap = sqrt(pi) * primitive_width * exp(-distance^2 / (4primitive_width^2))
            factor = kinetic ? (1 - distance^2 / (2primitive_width^2)) / (4primitive_width^2) : 1.0
            value_out += correlation[k] * overlap * factor
        end
        value_out
    end
    donor_overlap = representation.basis_matrices.overlap
    donor_kinetic = representation.basis_matrices.kinetic
    overlap_error = maximum(abs(pair_value((a - b) * spacing, false) - donor_overlap[a, b])
        for a in axes(donor_overlap, 1), b in axes(donor_overlap, 2))
    kinetic_error = maximum(abs(pair_value((a - b) * spacing, true) - donor_kinetic[a, b])
        for a in axes(donor_kinetic, 1), b in axes(donor_kinetic, 2))
    weights = integral_weights(prototype)
    weight = Float64(weights[midpoint])
    weight_error = maximum(abs(value - weight) for value in weights)
    maximum((overlap_error, kinetic_error, weight_error)) <= 50tolerance ||
        throw(ArgumentError("compact G10 realization disagrees with bounded donor"))
    isfinite(weight) && abs(weight) > tolerance ||
        throw(ArgumentError("sliced G10 integral weight is nonfinite or zero"))
    max_product_distance = 2cld((ncoeff - 1) ÷ 2, 3) + 24
    product_terms = zeros(Float64, 2ncoeff - 1, max_product_distance + 1)
    product_bounds = zeros(Float64, max_product_distance + 1)
    for distance in 0:max_product_distance, a in 1:ncoeff, b in 1:ncoeff
        separation = distance * spacing + (a - b) * primitive_width
        term = coefficients_local[a] * coefficients_local[b] * sqrt(pi) * primitive_width *
            exp(-separation^2 / (4primitive_width^2))
        product_terms[a + b - 1, distance + 1] += term
        product_bounds[distance + 1] += abs(term)
    end
    primitive_weight = sqrt(2pi) * primitive_width
    ida_moments = ntuple(4) do moment_index
        order = 2moment_index
        sum(coefficients_local[a] * primitive_weight * _sliced_even_moment(
            (a - (ncoeff + 1) / 2) * primitive_width, primitive_width^2, order)
            for a in 1:ncoeff) / weight
    end
    return _SlicedLongitudinal(centers, lattice, classes, spacing, coefficients_local,
        correlation, product_terms, product_bounds, weight, ida_moments,
        (overlap_error, kinetic_error, weight_error))
end

@inline function _sliced_longitudinal_pair(longitudinal::_SlicedLongitudinal,
    delta::Float64, kinetic::Bool)
    width = longitudinal.spacing / 3
    midpoint = (length(longitudinal.correlation) + 1) ÷ 2
    value_out = 0.0
    @inbounds for k in eachindex(longitudinal.correlation)
        distance = delta + (k - midpoint) * width
        overlap = sqrt(pi) * width * exp(-distance^2 / (4width^2))
        factor = kinetic ? (1 - distance^2 / (2width^2)) / (4width^2) : 1.0
        value_out += longitudinal.correlation[k] * overlap * factor
    end
    return value_out
end

function _sliced_longitudinal_pair_factor(longitudinal::_SlicedLongitudinal,
    delta::Float64, exponent::Float64)
    width = longitudinal.spacing / 3
    midpoint = (length(longitudinal.correlation) + 1) ÷ 2
    alpha = inv(width^2)
    determinant = alpha^2 + 4exponent * alpha
    prefactor = 2pi / sqrt(determinant)
    decay = exponent * alpha^2 / determinant
    value_out = 0.0
    @inbounds for k in eachindex(longitudinal.correlation)
        distance = delta + (k - midpoint) * width
        value_out += longitudinal.correlation[k] * prefactor * exp(-decay * distance^2)
    end
    return value_out
end

function _sliced_longitudinal_product_moment(longitudinal::_SlicedLongitudinal,
    distance::Int, order::Int)
    width = longitudinal.spacing / 3
    value_out = 0.0
    @inbounds for sum_index in axes(longitudinal.product_terms, 1)
        offset = (sum_index - length(longitudinal.coefficients)) * width / 2
        value_out += longitudinal.product_terms[sum_index, distance + 1] *
            _sliced_even_moment(offset, width^2 / 2, order)
    end
    return value_out
end

function _sliced_transverse_vector(z::Float64, atom_positions, exponents, contraction,
    prefactors, overlap, sqrt_overlap, invsqrt_overlap, tolerance::Float64)
    rho = zeros(Float64, length(exponents), length(exponents))
    for atom_z in atom_positions
        primitive = contraction .* prefactors .* exp.(-exponents .* (z - atom_z)^2)
        eta = sqrt_overlap * primitive
        rho .+= eta * eta'
    end
    decomposition = eigen(Symmetric(rho))
    lambda1 = decomposition.values[end]
    lambda2 = decomposition.values[end - 1]
    isfinite(lambda1) && lambda1 > 0 ||
        throw(ArgumentError("sliced transverse density lost numerical rank"))
    minimum(decomposition.values) >= -tolerance * lambda1 ||
        throw(ArgumentError("sliced transverse density is materially indefinite"))
    gap = (lambda1 - lambda2) / lambda1
    gap > tolerance || throw(ArgumentError("sliced transverse leading eigenspace is unresolved"))
    vector = invsqrt_overlap * decomposition.vectors[:, end]
    pivot = argmax(abs.(vector))
    vector[pivot] < 0 && (vector .*= -1)
    norm_error = abs(dot(vector, overlap * vector) - 1)
    norm_error <= 20tolerance || throw(ArgumentError("sliced transverse normalization failed"))
    rank = count(value -> value > tolerance * lambda1, decomposition.values)
    return vector, rank, gap, pivot, norm_error
end

function _sliced_transverse(nuclei, site_positions, mode::Symbol, sites_per_atom::Int,
    expansion::CoulombGaussianExpansion, tolerance::Float64, atom_spacing::Float64,
    phase_offset::Float64)
    source = legacy_s_gaussian_data("H", "STO-6G")
    length(source.shells) == 1 || throw(ArgumentError("sliced H/STO-6G requires one s shell"))
    shell = only(source.shells)
    exponents = copy(shell.exponents)
    contraction = copy(shell.coefficients)
    prefactors = Float64[GaussianAnalyticIntegrals.polynomial_gaussian_shell_prefactor(a, 0)
        for a in exponents]
    nprimitive = length(exponents)
    overlap1 = [GaussianAnalyticIntegrals.polynomial_gaussian_basic_integral(exponents[p],
        0.0, 0, prefactors[p], exponents[q], 0.0, 0, prefactors[q])
        for p in 1:nprimitive, q in 1:nprimitive]
    kinetic1 = [GaussianAnalyticIntegrals.polynomial_gaussian_kinetic_integral(exponents[p],
        0.0, 0, prefactors[p], exponents[q], 0.0, 0, prefactors[q])
        for p in 1:nprimitive, q in 1:nprimitive]
    overlap2 = overlap1 .^ 2
    kinetic2 = 2 .* overlap1 .* kinetic1
    source_norm = dot(contraction, (overlap1 .^ 3) * contraction)
    isfinite(source_norm) && source_norm > tolerance ||
        throw(ArgumentError("sliced H/STO-6G contracted source has invalid norm"))
    source_norm_error = abs(source_norm - 1)
    contraction ./= sqrt(source_norm)
    metric = eigen(Symmetric(overlap2))
    minimum(metric.values) > tolerance ||
        throw(ArgumentError("sliced transverse primitive metric lost rank"))
    sqrt_overlap = metric.vectors * Diagonal(sqrt.(metric.values)) * metric.vectors'
    invsqrt_overlap = metric.vectors * Diagonal(inv.(sqrt.(metric.values))) * metric.vectors'
    spacing = atom_spacing / sites_per_atom
    positions = mode == :periodic_template ?
        phase_offset .+ collect(0:(sites_per_atom - 1)) .* spacing : site_positions
    nclass = length(positions)
    vectors = zeros(Float64, nprimitive, nclass)
    ranks = zeros(Int, nclass)
    gaps = zeros(Float64, nclass)
    pivots = zeros(Int, nclass)
    max_norm_error = 0.0
    omitted_max = 0.0
    periodic_error = 0.0
    for site in eachindex(positions)
        atoms = nuclei
        shifted_atoms = nuclei
        if mode == :periodic_template
            window = max(2, ceil(Int, sqrt(-log(tolerance)) /
                (sqrt(2minimum(exponents)) * atom_spacing)) + 1)
            atoms = collect((-window):window) .* atom_spacing
            shifted_atoms = collect((-window + 1):(window + 1)) .* atom_spacing
            omitted = 0.0
            for cell in (window + 1):(window + 64), sign in (-1, 1)
                primitive = contraction .* prefactors .*
                    exp.(-exponents .* (positions[site] - sign * cell * atom_spacing)^2)
                omitted += sqrt(abs(dot(primitive, overlap2 * primitive)))
            end
            omitted_max = max(omitted_max, omitted)
            omitted_max <= 10tolerance ||
                throw(ArgumentError("periodic transverse projection window did not converge"))
        end
        vector, ranks[site], gaps[site], pivots[site], norm_error =
            _sliced_transverse_vector(positions[site], atoms, exponents, contraction,
                prefactors, overlap2, sqrt_overlap, invsqrt_overlap, tolerance)
        vectors[:, site] .= vector
        max_norm_error = max(max_norm_error, norm_error)
        if mode == :periodic_template
            shifted, _, _, _, _ = _sliced_transverse_vector(positions[site] + atom_spacing,
                shifted_atoms, exponents, contraction, prefactors, overlap2,
                sqrt_overlap, invsqrt_overlap, tolerance)
            delta = vector - shifted
            periodic_error = max(periodic_error, sqrt(abs(dot(delta, overlap2 * delta))))
        end
    end
    periodic_error <= 20tolerance ||
        throw(ArgumentError("periodic sliced transverse classes failed translation parity"))
    radial_moments = zeros(Float64, 4, nprimitive, nprimitive)
    for order_index in 1:4, p in 1:nprimitive, q in 1:nprimitive
        radial_moments[order_index, p, q] = overlap2[p, q] * factorial(order_index) /
            (exponents[p] + exponents[q])^order_index
    end
    factor_terms = zeros(Float64, length(expansion), nprimitive, nprimitive)
    for k in eachindex(expansion.exponents), p in 1:nprimitive, q in 1:nprimitive
        value = GaussianAnalyticIntegrals.polynomial_gaussian_basic_integral(exponents[p],
            0.0, 0, prefactors[p], exponents[q], 0.0, 0, prefactors[q];
            extra_exponent = expansion.exponents[k])
        factor_terms[k, p, q] = value^2
    end
    pair_exponents = Float64[exponents[p] + exponents[q]
        for p in 1:nprimitive for q in 1:nprimitive]
    density = zeros(Float64, nprimitive^2, nclass)
    for site in 1:nclass, p in 1:nprimitive, q in 1:nprimitive
        density[(p - 1) * nprimitive + q, site] = vectors[p, site] * vectors[q, site] *
            (prefactors[p] * prefactors[q])^2
    end
    transverse = _SlicedTransverse(exponents, vectors, overlap2, kinetic2,
        radial_moments, factor_terms, pair_exponents, density, ranks, gaps, pivots)
    fingerprint = bytes2hex(sha256(reinterpret(UInt8, vcat(exponents, contraction))))
    return transverse, max_norm_error, omitted_max, periodic_error, source_norm_error, fingerprint
end

@inline function _sliced_transverse_bilinear(transverse::_SlicedTransverse,
    matrix::Matrix{Float64}, a::Int, b::Int)
    value_out = 0.0
    @inbounds for p in axes(matrix, 1), q in axes(matrix, 2)
        value_out += transverse.coefficients[p, a] * transverse.coefficients[q, b] * matrix[p, q]
    end
    return value_out
end

@inline function _sliced_transverse_moment(transverse::_SlicedTransverse,
    order::Int, a::Int, b::Int)
    value_out = 0.0
    moment_index = order ÷ 2
    @inbounds for p in axes(transverse.coefficients, 1), q in axes(transverse.coefficients, 1)
        value_out += transverse.coefficients[p, a] * transverse.coefficients[q, b] *
            transverse.radial_moments[moment_index, p, q]
    end
    return value_out
end

@inline function _sliced_transverse_factor(transverse::_SlicedTransverse,
    k::Int, a::Int, b::Int)
    value_out = 0.0
    @inbounds for p in axes(transverse.coefficients, 1), q in axes(transverse.coefficients, 1)
        value_out += transverse.coefficients[p, a] * transverse.coefficients[q, b] *
            transverse.factor_terms[k, p, q]
    end
    return value_out
end

function _sliced_vee_transverse(transverse::_SlicedTransverse,
    exponent::Float64, a::Int, b::Int)
    value_out = 0.0
    @inbounds for p in eachindex(transverse.pair_exponents), q in eachindex(transverse.pair_exponents)
        ep = transverse.pair_exponents[p]
        eq = transverse.pair_exponents[q]
        denominator = ep * eq + exponent * (ep + eq)
        value_out += transverse.density_coefficients[p, a] *
            transverse.density_coefficients[q, b] * pi^2 / denominator
    end
    return value_out
end

function _sliced_h1_multipoles(longitudinal::_SlicedLongitudinal,
    transverse::_SlicedTransverse, distance::Int, a::Int, b::Int)
    z0 = _sliced_longitudinal_product_moment(longitudinal, distance, 0)
    z2 = _sliced_longitudinal_product_moment(longitudinal, distance, 2)
    z4 = _sliced_longitudinal_product_moment(longitudinal, distance, 4)
    z6 = _sliced_longitudinal_product_moment(longitudinal, distance, 6)
    z8 = _sliced_longitudinal_product_moment(longitudinal, distance, 8)
    r0 = _sliced_transverse_bilinear(transverse, transverse.overlap, a, b)
    r2 = _sliced_transverse_moment(transverse, 2, a, b)
    r4 = _sliced_transverse_moment(transverse, 4, a, b)
    r6 = _sliced_transverse_moment(transverse, 6, a, b)
    r8 = _sliced_transverse_moment(transverse, 8, a, b)
    return (z0 * r0,
        (2z2 * r0 - z0 * r2) / 2,
        z4 * r0 - 3z2 * r2 + 3z0 * r4 / 8,
        z6 * r0 - 15z4 * r2 / 2 + 45z2 * r4 / 8 - 5z0 * r6 / 16,
        z8 * r0 - 14z6 * r2 + 105z4 * r4 / 4 - 35z2 * r6 / 4 + 35z0 * r8 / 128)
end

@inline function _sliced_tail_value(multipoles, distance::Float64)
    inverse2 = inv(distance^2)
    return inv(distance) * (multipoles[1] + inverse2 * (multipoles[2] +
        inverse2 * (multipoles[3] + inverse2 * multipoles[4])))
end

function _sliced_vee_multipoles(longitudinal::_SlicedLongitudinal,
    transverse::_SlicedTransverse, a::Int, b::Int)
    z2, z4, z6 = longitudinal.ida_moments[1:3]
    b2 = _sliced_transverse_moment(transverse, 2, a, a)
    c2 = _sliced_transverse_moment(transverse, 2, b, b)
    b4 = _sliced_transverse_moment(transverse, 4, a, a)
    c4 = _sliced_transverse_moment(transverse, 4, b, b)
    b6 = _sliced_transverse_moment(transverse, 6, a, a)
    c6 = _sliced_transverse_moment(transverse, 6, b, b)
    axial2 = 2z2
    radial2 = b2 + c2
    axial4 = 2z4 + 6z2^2
    radial4 = b4 + c4 + 4b2 * c2
    axial6 = 2z6 + 30z4 * z2
    radial6 = b6 + c6 + 9(b4 * c2 + b2 * c4)
    q2 = (2axial2 - radial2) / 2
    q4 = axial4 - 3axial2 * radial2 + 3radial4 / 8
    q6 = axial6 - 15axial4 * radial2 / 2 + 45axial2 * radial4 / 8 - 5radial6 / 16
    return q2, q4, q6
end

function _sliced_longitudinal_nuclear(longitudinal::_SlicedLongitudinal,
    distance::Int, midpoint::Float64, exponent::Float64, nucleus::Float64)
    width = longitudinal.spacing / 3
    denominator = 1 + exponent * width^2
    value_out = 0.0
    @inbounds for sum_index in axes(longitudinal.product_terms, 1)
        center_value = midpoint + (sum_index - length(longitudinal.coefficients)) * width / 2
        value_out += longitudinal.product_terms[sum_index, distance + 1] /
            sqrt(denominator) * exp(-exponent * (center_value - nucleus)^2 / denominator)
    end
    return value_out
end

function _sliced_h1_band(longitudinal, transverse, expansion, natoms::Int, tolerance::Float64)
    nclass = size(transverse.coefficients, 2)
    kinetic_max = maximum(_sliced_transverse_bilinear(transverse, transverse.kinetic, a, a)
        for a in 1:nclass)
    factor_max = [_sliced_transverse_factor(transverse, k, a, a)
        for k in eachindex(expansion.coefficients), a in 1:nclass]
    bounds = zeros(Float64, length(longitudinal.product_bounds))
    for distance in eachindex(bounds)
        d = distance - 1
        bounds[distance] = abs(_sliced_longitudinal_pair(longitudinal,
            d * longitudinal.spacing, true)) +
            abs(_sliced_longitudinal_pair(longitudinal, d * longitudinal.spacing, false)) * kinetic_max
        for k in eachindex(expansion.coefficients)
            bounds[distance] += natoms * expansion.coefficients[k] * maximum(view(factor_max, k, :)) *
                longitudinal.product_bounds[distance] /
                sqrt(1 + expansion.exponents[k] * (longitudinal.spacing / 3)^2)
        end
    end
    band = findlast(value -> value > tolerance, bounds)
    band = isnothing(band) ? 0 : band - 1
    offband = band + 2 <= length(bounds) ? maximum(view(bounds, (band + 2):length(bounds))) : 0.0
    bounds[end] <= tolerance * 1e-6 ||
        throw(ArgumentError("sliced H1 product envelope did not reach its requested bound"))
    return band, offband
end

function _sliced_build_coulomb(longitudinal, transverse, expansion, mode, sites_per_atom,
    natoms, atom_spacing, one_body_tolerance, interaction_tolerance)
    h1_band, h1_offband = _sliced_h1_band(longitudinal, transverse, expansion,
        natoms, one_body_tolerance / 4)
    nclass = size(transverse.coefficients, 2)
    radial_bounds = ntuple(order_index -> maximum(abs(_sliced_transverse_moment(
        transverse, 2order_index, a, a)) for a in 1:nclass), 3)
    z2, z4, z6 = abs.(longitudinal.ida_moments[1:3])
    axial2 = 2z2
    axial4 = 2z4 + 6z2^2
    axial6 = 2z6 + 30z4 * z2
    radial2 = 2radial_bounds[1]
    radial4 = 2radial_bounds[2] + 4radial_bounds[1]^2
    radial6 = 2radial_bounds[3] + 18radial_bounds[2] * radial_bounds[1]
    vee_q6_bound = axial6 + 15axial4 * radial2 / 2 +
        45axial2 * radial4 / 8 + 5radial6 / 16
    vee_tail_distance = 20.0
    while 4vee_q6_bound / vee_tail_distance^7 > interaction_tolerance / 2
        vee_tail_distance += longitudinal.spacing
    end
    h1_q8_bound = 0.0
    transverse_bounds = ntuple(order_index -> maximum(abs(_sliced_transverse_moment(
        transverse, 2order_index, a, a)) for a in 1:nclass), 4)
    for distance in 0:h1_band
        z = ntuple(order_index -> abs(_sliced_longitudinal_product_moment(
            longitudinal, distance, 2order_index)), 4)
        z0 = abs(_sliced_longitudinal_product_moment(longitudinal, distance, 0))
        h1_q8_bound = max(h1_q8_bound, z[4] + 14z[3] * transverse_bounds[1] +
            105z[2] * transverse_bounds[2] / 4 + 35z[1] * transverse_bounds[3] / 4 +
            35z0 * transverse_bounds[4] / 128)
    end
    h1_tail_distance = 20.0
    cumulative_bound(distance) = 4h1_q8_bound *
        (2 / distance^9 + 1 / (4atom_spacing * distance^8))
    while cumulative_bound(h1_tail_distance) > one_body_tolerance / 4
        h1_tail_distance += longitudinal.spacing
    end
    maximum((vee_tail_distance, h1_tail_distance)) < 400 ||
        throw(ArgumentError("sliced Coulomb tail cannot meet requested tolerances in the validated range"))
    near_count = min(length(longitudinal.centers) - 1,
        ceil(Int, vee_tail_distance / longitudinal.spacing))
    longitudinal_terms = zeros(Float64, length(expansion), near_count + 1)
    for distance in 0:near_count, k in eachindex(expansion.exponents)
        longitudinal_terms[k, distance + 1] = _sliced_longitudinal_pair_factor(
            longitudinal, distance * longitudinal.spacing, expansion.exponents[k])
    end
    transverse_terms = mode == :periodic_template ?
        zeros(Float64, length(expansion), nclass, nclass) : zeros(Float64, 0, 0, 0)
    if mode == :periodic_template
        for k in eachindex(expansion.exponents), a in 1:nclass, b in 1:nclass
            transverse_terms[k, a, b] = _sliced_vee_transverse(
                transverse, expansion.exponents[k], a, b)
        end
    end
    dmax = cld(last(longitudinal.lattice_indices) - first(longitudinal.lattice_indices),
        sites_per_atom) + 1
    blocks = mode == :periodic_template ?
        zeros(Float64, sites_per_atom, sites_per_atom, dmax + 1) : zeros(Float64, 0, 0, 0)
    vee_transition_error = 0.0
    if mode == :periodic_template
        weight2 = longitudinal.integral_weight^2
        for cell_distance in 0:dmax, a in 1:sites_per_atom, b in 1:sites_per_atom
            site_distance = abs(cell_distance * sites_per_atom + b - a)
            physical_distance = site_distance * longitudinal.spacing
            if site_distance <= near_count
                value_out = 0.0
                @inbounds for k in eachindex(expansion.coefficients)
                    value_out += expansion.coefficients[k] * transverse_terms[k, a, b] *
                        longitudinal_terms[k, site_distance + 1]
                end
                blocks[a, b, cell_distance + 1] = value_out / weight2
                if physical_distance >= vee_tail_distance - longitudinal.spacing && physical_distance > 0
                    q2, q4, _ = _sliced_vee_multipoles(longitudinal, transverse, a, b)
                    tail = inv(physical_distance) + q2 / physical_distance^3 + q4 / physical_distance^5
                    vee_transition_error = max(vee_transition_error, abs(value_out / weight2 - tail))
                end
            else
                q2, q4, _ = _sliced_vee_multipoles(longitudinal, transverse, a, b)
                blocks[a, b, cell_distance + 1] = inv(physical_distance) +
                    q2 / physical_distance^3 + q4 / physical_distance^5
            end
        end
    else
        transition_site_distance = ceil(Int, vee_tail_distance / longitudinal.spacing)
        transition_distance = transition_site_distance * longitudinal.spacing
        classes = unique((1, cld(nclass, 2), nclass))
        weight2 = longitudinal.integral_weight^2
        for a in classes, b in classes
            exact_value = 0.0
            @inbounds for k in eachindex(expansion.coefficients)
                exact_value += expansion.coefficients[k] *
                    _sliced_vee_transverse(transverse, expansion.exponents[k], a, b) *
                    _sliced_longitudinal_pair_factor(longitudinal, transition_distance,
                        expansion.exponents[k])
            end
            q2, q4, _ = _sliced_vee_multipoles(longitudinal, transverse, a, b)
            tail_value = inv(transition_distance) + q2 / transition_distance^3 +
                q4 / transition_distance^5
            vee_transition_error = max(vee_transition_error,
                abs(exact_value / weight2 - tail_value))
        end
    end
    h1_transition_error = 0.0
    for distance in 0:h1_band
        pair_count = mode == :periodic_template ? sites_per_atom :
            min(nclass - distance, max(1, sites_per_atom))
        for a in 1:pair_count
            b = mode == :periodic_template ? mod1(a + distance, sites_per_atom) : a + distance
            multipoles = _sliced_h1_multipoles(longitudinal, transverse, distance, a, b)
            exact_value = 0.0
            @inbounds for k in eachindex(expansion.coefficients)
                exact_value += expansion.coefficients[k] * _sliced_transverse_factor(transverse, k, a, b) *
                    _sliced_longitudinal_nuclear(longitudinal, distance, 0.0,
                        expansion.exponents[k], h1_tail_distance)
            end
            h1_transition_error = max(h1_transition_error,
                abs(exact_value - _sliced_tail_value(multipoles, h1_tail_distance)))
        end
    end
    vee_tail_bound = 4vee_q6_bound / vee_tail_distance^7
    h1_tail_bound = cumulative_bound(h1_tail_distance)
    vee_tail_bound + vee_transition_error <= interaction_tolerance ||
        throw(ArgumentError("sliced Vee transition exceeds requested tolerance"))
    h1_offband + h1_tail_bound + h1_transition_error <= one_body_tolerance ||
        throw(ArgumentError("sliced H1 transition exceeds requested tolerance"))
    return _SlicedCoulomb(expansion, longitudinal_terms, transverse_terms, blocks,
        h1_band, h1_offband, vee_tail_distance, vee_tail_bound, vee_transition_error,
        h1_tail_distance, h1_tail_bound, h1_transition_error)
end

function _sliced_inverse_sums(longitudinal, first_atom_lattice::Int, natoms::Int,
    sites_per_atom::Int)
    mmin = 2first(longitudinal.lattice_indices) - 2first_atom_lattice
    mmax = 2last(longitudinal.lattice_indices) - 2first_atom_lattice
    sums = zeros(Float64, 4, mmax - mmin + 1)
    period = 2sites_per_atom
    for residue in 0:(period - 1)
        first_m = mmin + mod(residue - mmin, period)
        first_m > mmax && continue
        first_q = fld(first_m, period)
        last_q = fld(mmax - mod(mmax - residue, period), period)
        nmin = first_q - natoms + 1
        nmax = last_q
        prefix = zeros(Float64, 4, nmax - nmin + 2)
        for n in nmin:nmax
            index = n - nmin + 1
            @inbounds prefix[:, index + 1] .= prefix[:, index]
            coordinate = period * n + residue
            if coordinate != 0
                inverse = 2 / (longitudinal.spacing * abs(coordinate))
                @inbounds for order_index in 1:4
                    prefix[order_index, index + 1] += inverse^(2order_index - 1)
                end
            end
        end
        for midpoint in first_m:period:mmax
            q = fld(midpoint, period)
            low = q - natoms + 1
            @inbounds for order_index in 1:4
                sums[order_index, midpoint - mmin + 1] =
                    prefix[order_index, q - nmin + 2] -
                    prefix[order_index, low - nmin + 1]
            end
        end
    end
    return sums, mmin
end

function _sliced_h1_bands(longitudinal, transverse, coulomb, nuclei,
    sites_per_atom::Int)
    nsites = length(longitudinal.centers)
    band = coulomb.h1_band
    first_atom_lattice = first(longitudinal.lattice_indices) +
        round(Int, (first(nuclei) - first(longitudinal.centers)) / longitudinal.spacing)
    inverse_sums, midpoint_min = _sliced_inverse_sums(longitudinal,
        first_atom_lattice, length(nuclei), sites_per_atom)
    max_relative_atom = ceil(Int, (coulomb.h1_tail_distance +
        (sites_per_atom + band) * longitudinal.spacing) /
        (sites_per_atom * longitudinal.spacing)) + 1
    nprimitive = size(transverse.coefficients, 1)
    exact_terms = zeros(Float64, nprimitive, nprimitive, band + 1,
        sites_per_atom, 2max_relative_atom + 1)
    for distance in 0:band, residue in 1:sites_per_atom,
        relative_atom in -max_relative_atom:max_relative_atom
        separation = longitudinal.spacing *
            (residue - 1 + distance / 2 - sites_per_atom * relative_atom)
        abs(separation) <= coulomb.h1_tail_distance || continue
        @inbounds for k in eachindex(coulomb.expansion.coefficients)
            longitudinal_value = _sliced_longitudinal_nuclear(longitudinal,
                distance, separation, coulomb.expansion.exponents[k], 0.0)
            scale = coulomb.expansion.coefficients[k] * longitudinal_value
            for p in 1:nprimitive, q in 1:nprimitive
                exact_terms[p, q, distance + 1, residue,
                    relative_atom + max_relative_atom + 1] +=
                    scale * transverse.factor_terms[k, p, q]
            end
        end
    end
    bands = zeros(Float64, band + 1, nsites)
    for site in 1:nsites
        lattice = longitudinal.lattice_indices[site]
        cell = fld(lattice - first_atom_lattice, sites_per_atom)
        residue = mod(lattice - first_atom_lattice, sites_per_atom) + 1
        for distance in 0:min(band, nsites - site)
            other = site + distance
            a = longitudinal.class_indices[site]
            b = longitudinal.class_indices[other]
            sz = _sliced_longitudinal_pair(longitudinal,
                distance * longitudinal.spacing, false)
            tz = _sliced_longitudinal_pair(longitudinal,
                distance * longitudinal.spacing, true)
            local_value = tz * _sliced_transverse_bilinear(transverse,
                transverse.overlap, a, b) + sz * _sliced_transverse_bilinear(
                transverse, transverse.kinetic, a, b)
            multipoles = _sliced_h1_multipoles(longitudinal, transverse,
                distance, a, b)
            midpoint = 2lattice + distance - 2first_atom_lattice
            potential = 0.0
            @inbounds for order_index in 1:4
                potential += multipoles[order_index] *
                    inverse_sums[order_index, midpoint - midpoint_min + 1]
            end
            for relative_atom in -max_relative_atom:max_relative_atom
                atom = cell + relative_atom
                0 <= atom < length(nuclei) || continue
                separation = longitudinal.spacing *
                    (residue - 1 + distance / 2 - sites_per_atom * relative_atom)
                abs(separation) <= coulomb.h1_tail_distance || continue
                exact_value = 0.0
                @inbounds for p in 1:nprimitive, q in 1:nprimitive
                    exact_value += transverse.coefficients[p, a] *
                        transverse.coefficients[q, b] * exact_terms[p, q,
                            distance + 1, residue,
                            relative_atom + max_relative_atom + 1]
                end
                tail_value = iszero(separation) ? 0.0 :
                    _sliced_tail_value(multipoles, abs(separation))
                potential += exact_value - tail_value
            end
            bands[distance + 1, site] = local_value - potential
        end
    end
    return bands
end

"""
    sliced_hydrogen_chain(natoms; R, spacing, sites_per_atom, ...)

Build the bare minimal sliced-chain basis and compact H1/IntegralDiagonal Vee resources.
"""
function sliced_hydrogen_chain(natoms::Integer; R::Real, spacing = nothing,
    sites_per_atom = nothing, phase::Real = 0.0, left_padding::Real,
    right_padding::Real = left_padding, alignment::Symbol = :atom_centered,
    boundary::Symbol = :open, transverse_mode::Symbol = :finite,
    basis_name::AbstractString = "STO-6G", one_body_tolerance::Real = 1e-8,
    interaction_tolerance::Real = 1e-5)
    natoms > 0 || throw(ArgumentError("sliced hydrogen chain requires natoms > 0"))
    rv = Float64(R)
    phase_value = Float64(phase)
    left = Float64(left_padding)
    right = Float64(right_padding)
    all(isfinite, (rv, phase_value, left, right)) && rv > 0 && left >= 0 && right >= 0 ||
        throw(ArgumentError("sliced-chain geometry must be finite with positive R and padding"))
    boundary == :open || throw(ArgumentError("only open sliced-chain boundaries are supported"))
    alignment in (:atom_centered, :explicit_phase) ||
        throw(ArgumentError("alignment must be :atom_centered or :explicit_phase"))
    transverse_mode in (:finite, :periodic_template) ||
        throw(ArgumentError("unknown sliced transverse mode"))
    basis_name == "STO-6G" || throw(ArgumentError("initial sliced producer supports H/STO-6G"))
    h1_tolerance = Float64(one_body_tolerance)
    vee_tolerance = Float64(interaction_tolerance)
    isfinite(h1_tolerance) && 0 < h1_tolerance < 1e-2 &&
        isfinite(vee_tolerance) && 0 < vee_tolerance < 1e-2 ||
        throw(ArgumentError("invalid sliced-chain tolerance"))
    (spacing === nothing) != (sites_per_atom === nothing) ||
        throw(ArgumentError("supply exactly one of spacing or sites_per_atom"))
    b = sites_per_atom === nothing ? round(Int, rv / Float64(spacing)) : Int(sites_per_atom)
    b > 0 || throw(ArgumentError("sites_per_atom must be positive"))
    h = spacing === nothing ? rv / b : Float64(spacing)
    isfinite(h) && h > 0 && abs(rv / h - b) <= h1_tolerance ||
        throw(ArgumentError("R and longitudinal spacing must be commensurate"))
    nuclei = Float64[(atom - (natoms + 1) / 2) * rv for atom in 1:natoms]
    reference_index = floor(Int, (first(nuclei) - phase_value) / h + h1_tolerance)
    phase_offset = phase_value + reference_index * h - first(nuclei)
    alignment == :atom_centered && abs(phase_offset) > h1_tolerance &&
        throw(ArgumentError("atom-centered alignment requires nuclei on gausslet centers"))
    kmin = floor(Int, (first(nuclei) - left - phase_value) / h)
    kmax = ceil(Int, (last(nuclei) + right - phase_value) / h)
    lattice = collect(kmin:kmax)
    centers = phase_value .+ h .* lattice
    all(isfinite, centers) && all(diff(centers) .> 0) ||
        throw(ArgumentError("sliced longitudinal centers must be finite and distinct"))
    classes = transverse_mode == :periodic_template ?
        mod.(lattice .- reference_index, b) .+ 1 : collect(eachindex(centers))
    longitudinal = _sliced_longitudinal(h, centers, lattice, classes, h1_tolerance)
    overlap_error = maximum(abs(_sliced_longitudinal_pair(longitudinal, distance * h, false) -
        (distance == 0)) for distance in 0:min(length(centers) - 1, 80))
    overlap_error <= 50h1_tolerance ||
        throw(ArgumentError("uniform G10 longitudinal overlap check failed"))
    expansion = coulomb_gaussian_expansion(doacc = true)
    transverse, norm_error, omitted_norm, periodic_error, source_norm_error, fingerprint =
        _sliced_transverse(nuclei, centers, transverse_mode, b, expansion,
            h1_tolerance, rv, phase_offset)
    coulomb = _sliced_build_coulomb(longitudinal, transverse, expansion, transverse_mode,
        b, natoms, rv, h1_tolerance, vee_tolerance)
    h1_bands = _sliced_h1_bands(longitudinal, transverse, coulomb, nuclei, b)
    nuclear_repulsion = sum((natoms - distance) / (distance * rv)
        for distance in 1:(natoms - 1); init = 0.0)
    retained_bytes = Base.summarysize(longitudinal) + Base.summarysize(transverse) +
        Base.summarysize(coulomb) + Base.summarysize(h1_bands)
    diagnostics = _SlicedDiagnostics(fingerprint, overlap_error, norm_error, periodic_error,
        omitted_norm, source_norm_error, first(nuclei) - first(centers),
        last(centers) - last(nuclei), last(centers) - first(centers), retained_bytes)
    return SlicedHydrogenChain(nuclei, rv, b, transverse_mode, longitudinal, transverse,
        coulomb, h1_bands, nuclear_repulsion, h1_tolerance, vee_tolerance, diagnostics)
end

Base.length(chain::SlicedHydrogenChain) = length(chain.longitudinal.centers)

"""
    sliced_h1(chain::SlicedHydrogenChain)

Return a read-only, lazily evaluated `AbstractMatrix{Float64}` view of the
electronic one-body operator. The view retains `chain` and supports matrix
size, scalar indexing, and `mul!` without owning an `N x N` dense matrix;
`chain.nuclear_repulsion` remains a separate scalar.

Use `sliced_h1_bandwidth(chain)` for the represented structural
half-bandwidth. Entries outside it are exactly zero in this representation,
not necessarily in the unrepresented continuum operator. The concrete view
type is private; materialize `Matrix(sliced_h1(chain))` only for bounded
diagnostics.
"""
sliced_h1(chain::SlicedHydrogenChain) = _SlicedH1View(chain)

"""
    sliced_vee(chain::SlicedHydrogenChain)

Return a read-only, lazily evaluated two-index density-density interaction
view, not a four-index electron-repulsion tensor. The view retains `chain` and
supports matrix size, scalar indexing, and `mul!` without owning an `N x N`
dense matrix. Its concrete type is private; explicit `Matrix` materialization
is intended only for bounded diagnostics.
"""
sliced_vee(chain::SlicedHydrogenChain) = _SlicedVeeView(chain)

"""
    sliced_h1_bandwidth(chain::SlicedHydrogenChain)::Int

Return the structural half-bandwidth of `sliced_h1(chain)`. Entries farther
than this distance are exactly zero. The query reads compact construction
state and does not inspect H1 coefficients.
"""
@inline sliced_h1_bandwidth(chain::SlicedHydrogenChain)::Int =
    min(chain.coulomb.h1_band, length(chain) - 1)

Base.size(view::Union{_SlicedH1View,_SlicedVeeView}) = (length(view.chain), length(view.chain))
Base.eltype(::Type{<:_SlicedH1View}) = Float64
Base.eltype(::Type{<:_SlicedVeeView}) = Float64

function Base.getindex(view::_SlicedH1View, i::Int, j::Int)
    checkbounds(view, i, j)
    distance = abs(i - j)
    distance > view.chain.coulomb.h1_band && return 0.0
    return view.chain.h1_bands[distance + 1, min(i, j)]
end

function Base.getindex(view::_SlicedVeeView, i::Int, j::Int)
    checkbounds(view, i, j)
    chain = view.chain
    longitudinal = chain.longitudinal
    coulomb = chain.coulomb
    ai = longitudinal.class_indices[i]
    aj = longitudinal.class_indices[j]
    if chain.mode == :periodic_template
        cell_i = fld(longitudinal.lattice_indices[i] - (ai - 1), chain.sites_per_atom)
        cell_j = fld(longitudinal.lattice_indices[j] - (aj - 1), chain.sites_per_atom)
        return cell_j >= cell_i ? coulomb.periodic_blocks[ai, aj, cell_j - cell_i + 1] :
            coulomb.periodic_blocks[aj, ai, cell_i - cell_j + 1]
    end
    site_distance = abs(longitudinal.lattice_indices[i] - longitudinal.lattice_indices[j])
    physical_distance = site_distance * longitudinal.spacing
    if physical_distance > coulomb.vee_tail_distance
        q2, q4, _ = _sliced_vee_multipoles(longitudinal, chain.transverse, ai, aj)
        return inv(physical_distance) + q2 / physical_distance^3 + q4 / physical_distance^5
    end
    value_out = 0.0
    @inbounds for k in eachindex(coulomb.expansion.coefficients)
        value_out += coulomb.expansion.coefficients[k] *
            _sliced_vee_transverse(chain.transverse, coulomb.expansion.exponents[k], ai, aj) *
            coulomb.longitudinal_terms[k, site_distance + 1]
    end
    return value_out / longitudinal.integral_weight^2
end

"""
    sliced_row!(destination, view, row)

Validate `destination` length and `row`, completely overwrite the
caller-owned destination, and return it. A compatible buffer follows the
zero-allocation path. H1 extraction evaluates only its represented band and
writes zeros elsewhere; Vee extraction evaluates every column.

Obtain supported views from `sliced_h1` or `sliced_vee`; their concrete types
are private.
"""
function sliced_row!(destination::AbstractVector, view::_SlicedH1View, row::Integer)
    length(destination) == size(view, 2) ||
        throw(DimensionMismatch("sliced row destination has wrong length"))
    row_index = Int(row)
    checkbounds(view, row_index, 1)
    fill!(destination, zero(eltype(destination)))
    band = view.chain.coulomb.h1_band
    @inbounds for column in max(1, row_index - band):min(length(destination), row_index + band)
        destination[column] = view[row_index, column]
    end
    return destination
end

function sliced_row!(destination::AbstractVector, view::_SlicedVeeView, row::Integer)
    length(destination) == size(view, 2) ||
        throw(DimensionMismatch("sliced row destination has wrong length"))
    row_index = Int(row)
    checkbounds(view, row_index, 1)
    @inbounds for column in eachindex(destination)
        destination[column] = view[row_index, column]
    end
    return destination
end

function LinearAlgebra.mul!(destination::AbstractVector, view::_SlicedH1View,
    source::AbstractVector)
    length(source) == size(view, 2) && length(destination) == size(view, 1) ||
        throw(DimensionMismatch("sliced H1 action dimensions do not match"))
    band = view.chain.coulomb.h1_band
    @inbounds for row in eachindex(destination)
        value_out = zero(eltype(destination))
        for column in max(1, row - band):min(length(source), row + band)
            value_out += view[row, column] * source[column]
        end
        destination[row] = value_out
    end
    return destination
end

function LinearAlgebra.mul!(destination::AbstractVector, view::_SlicedVeeView,
    source::AbstractVector)
    length(source) == size(view, 2) && length(destination) == size(view, 1) ||
        throw(DimensionMismatch("sliced Vee action dimensions do not match"))
    @inbounds for row in eachindex(destination)
        value_out = zero(eltype(destination))
        for column in eachindex(source)
            value_out += view[row, column] * source[column]
        end
        destination[row] = value_out
    end
    return destination
end
