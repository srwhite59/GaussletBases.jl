struct _AtomicReferencePairDensityTerm
    coefficient::Float64
    exponent::Float64
    center::NTuple{3,Float64}
    powers::NTuple{3,Int}
end

struct _AtomicReferenceHartreeGGDiagnostics
    reference_orbital_count::Int
    reference_pair_count::Int
    primitive_pair_term_count::Int
    compressed_pair_term_count::Int
    coulomb_expansion_term_count::Int
    one_body_packet_term_count::Int
    reference_electron_count::Float64
    reference_density_trace::Float64
    density_symmetry_error::Float64
    density_finite::Bool
    gg_dimension::Int
    gg_trace::Float64
    gg_symmetry_error::Float64
    gg_finite::Bool
end

function _atomic_reference_center(supplement)
    isempty(supplement.orbitals) &&
        throw(ArgumentError("atomic reference Hartree GG requires at least one reference orbital"))
    center = _float_center(first(supplement.orbitals).center)
    for orbital in supplement.orbitals
        other = _float_center(orbital.center)
        maximum(abs.(Tuple(other[axis] - center[axis] for axis in 1:3))) <= 1.0e-12 ||
            throw(ArgumentError("atomic reference Hartree GG currently requires one same-center atomic reference"))
        orbital.primitive_normalization == :axiswise_normalized_cartesian_gaussian ||
            throw(ArgumentError("reference orbitals must use :axiswise_normalized_cartesian_gaussian primitives"))
    end
    return center
end

function _validate_atomic_reference_density(supplement, density)
    n = length(supplement.orbitals)
    size(density) == (n, n) ||
        throw(DimensionMismatch("atomic reference density matrix size must match reference orbital count"))
    all(isfinite, density) ||
        throw(ArgumentError("atomic reference density matrix must be finite"))
    symmetry_error = norm(Matrix{Float64}(density) - transpose(Matrix{Float64}(density)), Inf)
    symmetry_error <= 1.0e-10 ||
        throw(ArgumentError("atomic reference density matrix must be symmetric"))
    overlap = Matrix{Float64}(
        getfield(_GB_PARENT, :_cartesian_supplement_cross_overlap)(supplement, supplement))
    electron_count = tr(Matrix{Float64}(density) * overlap)
    return (; density = Matrix{Float64}(density), overlap, symmetry_error, electron_count)
end

function _reference_orbital_axis_prefactor(orbital, primitive, axis)
    exponent = Float64(orbital.exponents[primitive])
    power = orbital.angular_powers[axis]
    return getfield(_GB_PARENT.GaussianAnalyticIntegrals,
        :polynomial_gaussian_shell_prefactor)(exponent, power)
end

function _atomic_reference_pair_density_terms(supplement, density)
    center = _atomic_reference_center(supplement)
    raw_count = 0
    compressed = Dict{Tuple{Float64,NTuple{3,Int},NTuple{3,Float64}},Float64}()
    for b in eachindex(supplement.orbitals), a in eachindex(supplement.orbitals)
        scale = Float64(density[a, b])
        scale == 0.0 && continue
        left = supplement.orbitals[a]
        right = supplement.orbitals[b]
        for rb in eachindex(right.exponents), la in eachindex(left.exponents)
            raw_count += 1
            exponent = Float64(left.exponents[la]) + Float64(right.exponents[rb])
            powers = ntuple(axis -> left.angular_powers[axis] + right.angular_powers[axis], 3)
            coefficient = scale * Float64(left.coefficients[la]) *
                Float64(right.coefficients[rb])
            for axis in 1:3
                coefficient *= _reference_orbital_axis_prefactor(left, la, axis) *
                    _reference_orbital_axis_prefactor(right, rb, axis)
            end
            key = (exponent, powers, center)
            compressed[key] = get(compressed, key, 0.0) + coefficient
        end
    end
    terms = _AtomicReferencePairDensityTerm[]
    for ((exponent, powers, term_center), coefficient) in compressed
        coefficient == 0.0 && continue
        push!(terms, _AtomicReferencePairDensityTerm(
            Float64(coefficient), Float64(exponent), term_center, powers))
    end
    sort!(terms; by = term -> (term.exponent, term.powers, term.center))
    return terms, raw_count
end

function _mixed_hartree_even_moment(gamma::Float64, power::Int)
    isodd(power) && return 0.0
    value = sqrt(pi / gamma)
    for k in 1:div(power, 2)
        value *= (2k - 1) / (2gamma)
    end
    return value
end

function _mixed_hartree_axis_convolution_components(power::Int, density_exponent::Float64,
    coupling_exponent::Float64)
    total = density_exponent + coupling_exponent
    ratio = coupling_exponent / total
    components = Pair{Int,Float64}[]
    for moment_power in 0:power
        isodd(moment_power) && continue
        factor_power = power - moment_power
        coefficient = binomial(power, moment_power) * ratio^factor_power *
            _mixed_hartree_even_moment(total, moment_power)
        push!(components, factor_power => coefficient)
    end
    return components
end

function _mixed_hartree_primitive_polynomial_factor(a, b, exponent::Float64,
    center::Float64, power::Int)
    gaussian_exponent = getfield(_GB_PARENT.GaussianAnalyticIntegrals, :gaussian_exponent)
    alpha_a = gaussian_exponent(a)
    alpha_b = gaussian_exponent(b)
    gamma = alpha_a + alpha_b + exponent
    weighted_center = (alpha_a * a.center_value + alpha_b * b.center_value +
        exponent * center) / gamma
    constant = alpha_a * a.center_value^2 + alpha_b * b.center_value^2 +
        exponent * center^2 - gamma * weighted_center^2
    shift = weighted_center - center
    value = 0.0
    for moment_power in 0:power
        value += binomial(power, moment_power) * shift^(power - moment_power) *
            _mixed_hartree_even_moment(gamma, moment_power)
    end
    return exp(-constant) * value
end

function _mixed_hartree_axis_factor_matrix(layer, exponent::Float64, center::Float64,
    power::Int)
    primitive_layer = getfield(_GB_PARENT, :primitive_set)(layer)
    primitives = getfield(_GB_PARENT, :primitives)(primitive_layer)
    matrix = zeros(Float64, length(primitive_layer), length(primitive_layer))
    for column in axes(matrix, 2), row in 1:column
        value = _mixed_hartree_primitive_polynomial_factor(
            primitives[row], primitives[column], exponent, center, power)
        matrix[row, column] = value
        matrix[column, row] = value
    end
    return Matrix{Float64}(getfield(_GB_PARENT, :contract_primitive_matrix)(layer, matrix))
end

function _mixed_hartree_cached_axis_factor!(cache, axis, exponent, center, power)
    key = (objectid(axis.base_layer), Float64(exponent), Float64(center), Int(power))
    return get!(cache, key) do
        _mixed_hartree_axis_factor_matrix(
            axis.base_layer, Float64(exponent), Float64(center), Int(power))
    end
end

function _atomic_reference_hartree_factor_packets(bundles, terms, expansion)
    axes = Tuple(getfield(_GB_PARENT, :_nested_axis_pgdg)(bundles, axis) for axis in (:x, :y, :z))
    coefficients = Float64[]
    factors = (Matrix{Float64}[], Matrix{Float64}[], Matrix{Float64}[])
    cache = Dict{Tuple{UInt,Float64,Float64,Int},Matrix{Float64}}()
    for term in terms, expansion_index in eachindex(expansion.coefficients)
        coupling = Float64(expansion.exponents[expansion_index])
        factor_exponent = term.exponent * coupling / (term.exponent + coupling)
        axis_components = ntuple(axis -> _mixed_hartree_axis_convolution_components(
            term.powers[axis], term.exponent, coupling), 3)
        for xcomp in axis_components[1], ycomp in axis_components[2], zcomp in axis_components[3]
            coefficient = term.coefficient * Float64(expansion.coefficients[expansion_index]) *
                xcomp.second * ycomp.second * zcomp.second
            coefficient == 0.0 && continue
            push!(coefficients, coefficient)
            push!(factors[1], _mixed_hartree_cached_axis_factor!(
                cache, axes[1], factor_exponent, term.center[1], xcomp.first))
            push!(factors[2], _mixed_hartree_cached_axis_factor!(
                cache, axes[2], factor_exponent, term.center[2], ycomp.first))
            push!(factors[3], _mixed_hartree_cached_axis_factor!(
                cache, axes[3], factor_exponent, term.center[3], zcomp.first))
        end
    end
    return (; coefficients, x = factors[1], y = factors[2], z = factors[3])
end

function atomic_reference_hartree_gg_block(basis, bundles, reference_supplement,
    reference_density; expansion = getfield(_GB_PARENT, :coulomb_gaussian_expansion)(doacc = false))
    density = _validate_atomic_reference_density(reference_supplement, reference_density)
    terms, primitive_count = _atomic_reference_pair_density_terms(
        reference_supplement, density.density)
    packets = _atomic_reference_hartree_factor_packets(bundles, terms, expansion)
    GG = zeros(Float64, basis.final_dimension, basis.final_dimension)
    getfield(_GB_PARENT.CartesianFinalBasisRealization, :_accumulate_terminal_gaussian_sum!)(
        GG, basis, packets.coefficients, packets.x, packets.y, packets.z; scale = 1.0)
    GG = _symmetrize_raw_block(GG)
    diagnostics = _AtomicReferenceHartreeGGDiagnostics(
        length(reference_supplement.orbitals),
        length(reference_supplement.orbitals)^2,
        primitive_count,
        length(terms),
        length(expansion.coefficients),
        length(packets.coefficients),
        density.electron_count,
        tr(density.density * density.overlap),
        density.symmetry_error,
        all(isfinite, density.density),
        size(GG, 1),
        tr(GG),
        norm(GG - transpose(GG), Inf),
        all(isfinite, GG),
    )
    return (; GG, diagnostics)
end
