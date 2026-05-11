struct _GaussianCoulombAxisPairTerm
    alpha_left::Float64
    center_left::Float64
    power_left::Int
    prefactor_left::Float64
    alpha_right::Float64
    center_right::Float64
    power_right::Int
    prefactor_right::Float64
end

struct _GaussianCoulombAxisKernelTerm
    alpha_sum::Float64
    center_weight::Float64
    constant::Float64
    prefactor::Float64
    polynomial_coefficients::Vector{Float64}
end

struct _GaussianCoulombPairTerm3D
    coefficient::Float64
    x::_GaussianCoulombAxisPairTerm
    y::_GaussianCoulombAxisPairTerm
    z::_GaussianCoulombAxisPairTerm
end

struct _GaussianCoulombPairTermDescriptor3D
    x::_GaussianCoulombAxisPairTerm
    y::_GaussianCoulombAxisPairTerm
    z::_GaussianCoulombAxisPairTerm
end

struct _GaussianCoulombPairCoefficient
    term_index::Int
    coefficient::Float64
end

struct _GaussianCenteredPairTerm3D
    coefficient::Float64
    beta::Float64
    powers::NTuple{3,Int}
end

struct _GaussianCenteredTermDescriptor3D
    beta::Float64
    powers::NTuple{3,Int}
end

struct _GaussianCenteredPairCoefficient
    term_index::Int
    coefficient::Float64
end

"""
    gaussian_coulomb_pair_index(p, q, n)

Return the dense pair-index used by [`gaussian_coulomb_pair_matrix`](@ref).
The flattened pair `(p, q)` is stored at `(p - 1) * n + q`.
"""
function gaussian_coulomb_pair_index(p::Integer, q::Integer, n::Integer)
    n >= 0 || throw(ArgumentError("Gaussian Coulomb pair indexing requires n >= 0"))
    1 <= p <= n ||
        throw(ArgumentError("Gaussian Coulomb pair index p=$(p) is outside 1:$(n)"))
    1 <= q <= n ||
        throw(ArgumentError("Gaussian Coulomb pair index q=$(q) is outside 1:$(n)"))
    return (p - 1) * n + q
end

"""
    gaussian_coulomb_pair_matrix(orbitals_or_supplement; expansion = coulomb_gaussian_expansion(doacc = false), max_orbitals = 64)

Build a dense pure-Gaussian Coulomb reference matrix `G` with entries
`G[pq, rs] = (pq|rs)`, using the repo `CoulombGaussianExpansion`
approximation to `1/r`.

This utility is analytic and quadrature-free. It is intended for small
reference checks, Qiu-White/MWG residual-Gaussian validation, and
stationary-Fock/EGOI diagnostics. It is not a production large-system ERI
backend: both work and memory scale as `N^4`.

Accepted inputs include a vector of
`CartesianGaussianShellOrbitalRepresentation3D`, a
`CartesianGaussianShellSupplementRepresentation3D`, and
`LegacyAtomicGaussianSupplement`. The pair index convention is
`gaussian_coulomb_pair_index(p, q, n) == (p - 1) * n + q`.
"""
function gaussian_coulomb_pair_matrix(
    orbitals::AbstractVector{<:CartesianGaussianShellOrbitalRepresentation3D};
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    max_orbitals = 64,
)
    orbital_count = length(orbitals)
    _validate_gaussian_coulomb_orbital_count(orbital_count, max_orbitals)
    for orbital in orbitals
        _validate_gaussian_coulomb_orbital(orbital)
    end
    _gaussian_coulomb_orbitals_share_center(orbitals) &&
        return _gaussian_coulomb_pair_matrix_same_center(orbitals, expansion)
    return _gaussian_coulomb_pair_matrix_compressed_checked(orbitals, expansion)
end

function _gaussian_coulomb_pair_matrix_general(
    orbitals::AbstractVector{<:CartesianGaussianShellOrbitalRepresentation3D};
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    max_orbitals = 64,
)
    orbital_count = length(orbitals)
    _validate_gaussian_coulomb_orbital_count(orbital_count, max_orbitals)
    for orbital in orbitals
        _validate_gaussian_coulomb_orbital(orbital)
    end
    return _gaussian_coulomb_pair_matrix_general_checked(orbitals, expansion)
end

function _gaussian_coulomb_pair_matrix_general_checked(
    orbitals::AbstractVector{<:CartesianGaussianShellOrbitalRepresentation3D},
    expansion::CoulombGaussianExpansion,
)
    orbital_count = length(orbitals)
    pair_count = orbital_count^2
    pair_terms = Vector{Vector{_GaussianCoulombPairTerm3D}}(undef, pair_count)
    for p in 1:orbital_count, q in 1:orbital_count
        pair_terms[gaussian_coulomb_pair_index(p, q, orbital_count)] =
            _gaussian_coulomb_pair_terms(orbitals[p], orbitals[q])
    end

    matrix = zeros(Float64, pair_count, pair_count)
    for column in 1:pair_count, row in 1:column
        value = _gaussian_coulomb_pair_integral(pair_terms[row], pair_terms[column], expansion)
        matrix[row, column] = value
        matrix[column, row] = value
    end
    return matrix
end

function gaussian_coulomb_pair_matrix(
    supplement::CartesianGaussianShellSupplementRepresentation3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    max_orbitals = 64,
)
    return gaussian_coulomb_pair_matrix(
        supplement.orbitals;
        expansion,
        max_orbitals,
    )
end

function gaussian_coulomb_pair_matrix(
    supplement::LegacyAtomicGaussianSupplement;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    max_orbitals = 64,
)
    return gaussian_coulomb_pair_matrix(
        basis_representation(supplement);
        expansion,
        max_orbitals,
    )
end

function _validate_gaussian_coulomb_orbital_count(orbital_count::Int, max_orbitals)
    max_orbitals === nothing && return nothing
    max_count = Int(max_orbitals)
    max_count >= 0 || throw(ArgumentError("max_orbitals must be nonnegative or nothing"))
    orbital_count <= max_count ||
        throw(
            ArgumentError(
                "gaussian_coulomb_pair_matrix refuses $(orbital_count) orbitals with max_orbitals=$(max_count); this dense reference scales as N^4, so pass max_orbitals=nothing only for deliberate small-system diagnostics",
            ),
        )
    return nothing
end

function _gaussian_coulomb_orbitals_share_center(
    orbitals::AbstractVector{<:CartesianGaussianShellOrbitalRepresentation3D},
)
    isempty(orbitals) && return true
    reference_center = first(orbitals).center
    return all(orbital -> orbital.center == reference_center, orbitals)
end

function _gaussian_coulomb_pair_terms(
    left::CartesianGaussianShellOrbitalRepresentation3D,
    right::CartesianGaussianShellOrbitalRepresentation3D,
)
    terms = _GaussianCoulombPairTerm3D[]
    for right_primitive in eachindex(right.exponents), left_primitive in eachindex(left.exponents)
        coefficient =
            Float64(left.coefficients[left_primitive]) *
            Float64(right.coefficients[right_primitive])
        push!(
            terms,
            _GaussianCoulombPairTerm3D(
                coefficient,
                _gaussian_coulomb_axis_pair_term(left, right, left_primitive, right_primitive, 1),
                _gaussian_coulomb_axis_pair_term(left, right, left_primitive, right_primitive, 2),
                _gaussian_coulomb_axis_pair_term(left, right, left_primitive, right_primitive, 3),
            ),
        )
    end
    return terms
end

function _validate_gaussian_coulomb_orbital(
    orbital::CartesianGaussianShellOrbitalRepresentation3D,
)
    orbital.primitive_normalization == :axiswise_normalized_cartesian_gaussian ||
        throw(
            ArgumentError(
                "gaussian_coulomb_pair_matrix requires :axiswise_normalized_cartesian_gaussian primitive normalization",
            ),
        )
    length(orbital.exponents) == length(orbital.coefficients) ||
        throw(ArgumentError("Gaussian shell orbital requires matching exponent and coefficient counts"))
    all(exponent -> exponent > 0.0, orbital.exponents) ||
        throw(ArgumentError("Gaussian shell orbital exponents must be positive"))
    all(power -> power >= 0, orbital.angular_powers) ||
        throw(ArgumentError("Gaussian shell orbital angular powers must be nonnegative"))
    return nothing
end

function _gaussian_coulomb_axis_pair_term(
    left::CartesianGaussianShellOrbitalRepresentation3D,
    right::CartesianGaussianShellOrbitalRepresentation3D,
    left_primitive::Int,
    right_primitive::Int,
    axis::Int,
)
    alpha_left = Float64(left.exponents[left_primitive])
    alpha_right = Float64(right.exponents[right_primitive])
    power_left = left.angular_powers[axis]
    power_right = right.angular_powers[axis]
    return _GaussianCoulombAxisPairTerm(
        alpha_left,
        Float64(left.center[axis]),
        power_left,
        GaussianAnalyticIntegrals.polynomial_gaussian_shell_prefactor(alpha_left, power_left),
        alpha_right,
        Float64(right.center[axis]),
        power_right,
        GaussianAnalyticIntegrals.polynomial_gaussian_shell_prefactor(alpha_right, power_right),
    )
end

function _gaussian_coulomb_pair_integral(
    left_terms::Vector{_GaussianCoulombPairTerm3D},
    right_terms::Vector{_GaussianCoulombPairTerm3D},
    expansion::CoulombGaussianExpansion,
)
    value = 0.0
    for term_index in eachindex(expansion.coefficients)
        coupling_exponent = Float64(expansion.exponents[term_index])
        coefficient = Float64(expansion.coefficients[term_index])
        term_value = 0.0
        for left in left_terms, right in right_terms
            term_value +=
                left.coefficient *
                right.coefficient *
                _gaussian_coulomb_axis_integral(left.x, right.x, coupling_exponent) *
                _gaussian_coulomb_axis_integral(left.y, right.y, coupling_exponent) *
                _gaussian_coulomb_axis_integral(left.z, right.z, coupling_exponent)
        end
        value += coefficient * term_value
    end
    return value
end

function _gaussian_coulomb_pair_matrix_compressed_checked(
    orbitals::AbstractVector{<:CartesianGaussianShellOrbitalRepresentation3D},
    expansion::CoulombGaussianExpansion,
)
    orbital_count = length(orbitals)
    pair_count = orbital_count^2
    compact_pairs, compact_pair_index = _gaussian_coulomb_compact_pair_index(orbital_count)
    compact_count = length(compact_pairs)
    compact_terms = Vector{Vector{_GaussianCoulombPairTerm3D}}(undef, compact_count)
    for (compact_index, (p, q)) in pairs(compact_pairs)
        compact_terms[compact_index] = _gaussian_coulomb_pair_terms(orbitals[p], orbitals[q])
    end
    compact_coefficients, term_descriptors =
        _gaussian_coulomb_global_term_coefficients(compact_terms)
    term_kernel = _gaussian_coulomb_global_term_kernel(term_descriptors, expansion)
    compact_matrix = _gaussian_coulomb_compact_pair_matrix(
        compact_coefficients,
        term_kernel,
    )

    ordered_to_compact = Vector{Int}(undef, pair_count)
    for p in 1:orbital_count, q in 1:orbital_count
        ordered_to_compact[gaussian_coulomb_pair_index(p, q, orbital_count)] =
            compact_pair_index[p, q]
    end
    matrix = Matrix{Float64}(undef, pair_count, pair_count)
    Base.Threads.@threads :static for column in 1:pair_count
        compact_column = ordered_to_compact[column]
        for row in 1:pair_count
            matrix[row, column] = compact_matrix[ordered_to_compact[row], compact_column]
        end
    end
    return matrix
end

function _gaussian_coulomb_compact_pair_matrix(
    compact_coefficients::Vector{Vector{_GaussianCoulombPairCoefficient}},
    term_kernel::AbstractMatrix{Float64},
)
    compact_count = length(compact_coefficients)
    matrix = zeros(Float64, compact_count, compact_count)
    Base.Threads.@threads :static for column in 1:compact_count
        column_coefficients = compact_coefficients[column]
        for row in 1:column
            value = 0.0
            for right in column_coefficients
                right_index = right.term_index
                right_coefficient = right.coefficient
                for left in compact_coefficients[row]
                    value +=
                        left.coefficient *
                        right_coefficient *
                        term_kernel[left.term_index, right_index]
                end
            end
            matrix[row, column] = value
            matrix[column, row] = value
        end
    end
    return matrix
end

function _gaussian_coulomb_term_descriptor(term::_GaussianCoulombPairTerm3D)
    return _GaussianCoulombPairTermDescriptor3D(term.x, term.y, term.z)
end

function _gaussian_coulomb_global_term_coefficients(
    compact_terms::Vector{Vector{_GaussianCoulombPairTerm3D}},
)
    term_index_by_descriptor = Dict{_GaussianCoulombPairTermDescriptor3D,Int}()
    term_descriptors = _GaussianCoulombPairTermDescriptor3D[]
    compact_coefficients = Vector{Vector{_GaussianCoulombPairCoefficient}}(
        undef,
        length(compact_terms),
    )
    for compact_index in eachindex(compact_terms)
        local_coefficients = Dict{Int,Float64}()
        for term in compact_terms[compact_index]
            descriptor = _gaussian_coulomb_term_descriptor(term)
            term_index = get(term_index_by_descriptor, descriptor, 0)
            if term_index == 0
                push!(term_descriptors, descriptor)
                term_index = length(term_descriptors)
                term_index_by_descriptor[descriptor] = term_index
            end
            local_coefficients[term_index] =
                get(local_coefficients, term_index, 0.0) + term.coefficient
        end
        entries = collect(local_coefficients)
        sort!(entries; by = first)
        compact_coefficients[compact_index] = [
            _GaussianCoulombPairCoefficient(term_index, coefficient)
            for (term_index, coefficient) in entries if coefficient != 0.0
        ]
    end
    return compact_coefficients, term_descriptors
end

function _gaussian_coulomb_global_term_kernel(
    term_descriptors::Vector{_GaussianCoulombPairTermDescriptor3D},
    expansion::CoulombGaussianExpansion,
)
    term_count = length(term_descriptors)
    kernel = zeros(Float64, term_count, term_count)
    axis_terms, descriptor_axis_indices =
        _gaussian_coulomb_axis_term_indices(term_descriptors)
    axis_kernel_terms = _gaussian_coulomb_axis_kernel_terms(axis_terms)
    for expansion_index in eachindex(expansion.coefficients)
        coupling_exponent = Float64(expansion.exponents[expansion_index])
        coefficient = Float64(expansion.coefficients[expansion_index])
        axis_kernel = _gaussian_coulomb_axis_term_kernel(
            axis_kernel_terms,
            coupling_exponent,
        )
        Base.Threads.@threads :static for column in 1:term_count
            right = descriptor_axis_indices[column]
            for row in 1:column
                left = descriptor_axis_indices[row]
                value =
                    coefficient *
                    axis_kernel[left[1], right[1]] *
                    axis_kernel[left[2], right[2]] *
                    axis_kernel[left[3], right[3]]
                kernel[row, column] += value
                row == column || (kernel[column, row] += value)
            end
        end
    end
    return kernel
end

function _gaussian_coulomb_axis_term_indices(
    term_descriptors::Vector{_GaussianCoulombPairTermDescriptor3D},
)
    axis_index_by_term = Dict{_GaussianCoulombAxisPairTerm,Int}()
    axis_terms = _GaussianCoulombAxisPairTerm[]
    descriptor_axis_indices = Vector{NTuple{3,Int}}(undef, length(term_descriptors))
    for descriptor_index in eachindex(term_descriptors)
        descriptor = term_descriptors[descriptor_index]
        descriptor_axis_indices[descriptor_index] = (
            _gaussian_coulomb_axis_term_index!(axis_index_by_term, axis_terms, descriptor.x),
            _gaussian_coulomb_axis_term_index!(axis_index_by_term, axis_terms, descriptor.y),
            _gaussian_coulomb_axis_term_index!(axis_index_by_term, axis_terms, descriptor.z),
        )
    end
    return axis_terms, descriptor_axis_indices
end

function _gaussian_coulomb_axis_term_index!(
    axis_index_by_term::Dict{_GaussianCoulombAxisPairTerm,Int},
    axis_terms::Vector{_GaussianCoulombAxisPairTerm},
    axis_term::_GaussianCoulombAxisPairTerm,
)
    axis_index = get(axis_index_by_term, axis_term, 0)
    axis_index != 0 && return axis_index
    push!(axis_terms, axis_term)
    axis_index = length(axis_terms)
    axis_index_by_term[axis_term] = axis_index
    return axis_index
end

function _gaussian_coulomb_axis_kernel_terms(
    axis_terms::Vector{_GaussianCoulombAxisPairTerm},
)
    return [_gaussian_coulomb_axis_kernel_term(axis_term) for axis_term in axis_terms]
end

function _gaussian_coulomb_axis_kernel_term(axis_term::_GaussianCoulombAxisPairTerm)
    polynomial_coefficients = Float64[1.0]
    polynomial_coefficients = GaussianAnalyticIntegrals.polynomial_shift_multiply(
        polynomial_coefficients,
        -axis_term.center_left,
        axis_term.power_left,
    )
    polynomial_coefficients = GaussianAnalyticIntegrals.polynomial_shift_multiply(
        polynomial_coefficients,
        -axis_term.center_right,
        axis_term.power_right,
    )
    return _GaussianCoulombAxisKernelTerm(
        axis_term.alpha_left + axis_term.alpha_right,
        axis_term.alpha_left * axis_term.center_left +
        axis_term.alpha_right * axis_term.center_right,
        axis_term.alpha_left * axis_term.center_left^2 +
        axis_term.alpha_right * axis_term.center_right^2,
        axis_term.prefactor_left * axis_term.prefactor_right,
        polynomial_coefficients,
    )
end

function _gaussian_coulomb_axis_term_kernel(
    axis_terms::Vector{_GaussianCoulombAxisKernelTerm},
    coupling_exponent::Float64,
)
    axis_count = length(axis_terms)
    kernel = zeros(Float64, axis_count, axis_count)
    Base.Threads.@threads :static for column in 1:axis_count
        for row in 1:column
            value = _gaussian_coulomb_axis_integral(
                axis_terms[row],
                axis_terms[column],
                coupling_exponent,
            )
            kernel[row, column] = value
            kernel[column, row] = value
        end
    end
    return kernel
end

function _gaussian_coulomb_axis_integral(
    left::_GaussianCoulombAxisKernelTerm,
    right::_GaussianCoulombAxisKernelTerm,
    coupling_exponent::Float64,
)
    a11 = left.alpha_sum + coupling_exponent
    a22 = right.alpha_sum + coupling_exponent
    a12 = -coupling_exponent
    determinant = a11 * a22 - a12^2
    determinant > 0.0 ||
        throw(ArgumentError("polynomial Gaussian pair quadratic form must be positive definite"))

    d1 = left.center_weight
    d2 = right.center_weight
    constant = left.constant + right.constant
    quadratic_term = (a22 * d1^2 - 2.0 * a12 * d1 * d2 + a11 * d2^2) / determinant
    mean_x = (a22 * d1 - a12 * d2) / determinant
    mean_y = (-a12 * d1 + a11 * d2) / determinant
    sigma_xx = 0.5 * a22 / determinant
    sigma_yy = 0.5 * a11 / determinant
    sigma_xy = -0.5 * a12 / determinant

    moment_value = 0.0
    left_coefficients = left.polynomial_coefficients
    right_coefficients = right.polynomial_coefficients
    @inbounds for x_index in eachindex(left_coefficients)
        left_coefficient = left_coefficients[x_index]
        left_coefficient == 0.0 && continue
        for y_index in eachindex(right_coefficients)
            right_coefficient = right_coefficients[y_index]
            right_coefficient == 0.0 && continue
            moment_value +=
                left_coefficient *
                right_coefficient *
                GaussianAnalyticIntegrals.shifted_wick2_moment(
                    x_index - 1,
                    y_index - 1,
                    mean_x,
                    mean_y,
                    sigma_xx,
                    sigma_yy,
                    sigma_xy,
                )
        end
    end

    return left.prefactor *
           right.prefactor *
           (pi / sqrt(determinant)) *
           exp(-constant + quadratic_term) *
           moment_value
end

function _gaussian_coulomb_axis_integral(
    left::_GaussianCoulombAxisPairTerm,
    right::_GaussianCoulombAxisPairTerm,
    coupling_exponent::Float64,
)
    return _gaussian_coulomb_axis_integral(
        _gaussian_coulomb_axis_kernel_term(left),
        _gaussian_coulomb_axis_kernel_term(right),
        coupling_exponent,
    )
end

function _gaussian_coulomb_pair_matrix_same_center(
    orbitals::AbstractVector{<:CartesianGaussianShellOrbitalRepresentation3D},
    expansion::CoulombGaussianExpansion,
)
    orbital_count = length(orbitals)
    pair_count = orbital_count^2
    compact_pairs, compact_pair_index = _gaussian_coulomb_compact_pair_index(orbital_count)
    compact_count = length(compact_pairs)
    compact_terms = Vector{Vector{_GaussianCenteredPairTerm3D}}(undef, compact_count)
    for (compact_index, (p, q)) in pairs(compact_pairs)
        compact_terms[compact_index] = _gaussian_centered_pair_terms(orbitals[p], orbitals[q])
    end
    compact_coefficients, term_descriptors =
        _gaussian_centered_global_term_coefficients(compact_terms)
    term_kernel = _gaussian_centered_term_kernel(term_descriptors, expansion)
    compact_matrix = zeros(Float64, compact_count, compact_count)
    for column in 1:compact_count, row in 1:column
        value = _gaussian_centered_pair_integral(
            compact_coefficients[row],
            compact_coefficients[column],
            term_kernel,
        )
        compact_matrix[row, column] = value
        compact_matrix[column, row] = value
    end

    matrix = zeros(Float64, pair_count, pair_count)
    for p in 1:orbital_count, q in 1:orbital_count
        row = gaussian_coulomb_pair_index(p, q, orbital_count)
        compact_row = compact_pair_index[p, q]
        for r in 1:orbital_count, s in 1:orbital_count
            column = gaussian_coulomb_pair_index(r, s, orbital_count)
            matrix[row, column] = compact_matrix[compact_row, compact_pair_index[r, s]]
        end
    end
    return matrix
end

function _gaussian_coulomb_compact_pair_index(orbital_count::Int)
    compact_pairs = Vector{NTuple{2,Int}}()
    sizehint!(compact_pairs, div(orbital_count * (orbital_count + 1), 2))
    compact_pair_index = zeros(Int, orbital_count, orbital_count)
    compact_index = 0
    for p in 1:orbital_count, q in p:orbital_count
        compact_index += 1
        push!(compact_pairs, (p, q))
        compact_pair_index[p, q] = compact_index
        compact_pair_index[q, p] = compact_index
    end
    return compact_pairs, compact_pair_index
end

function _gaussian_centered_pair_terms(
    left::CartesianGaussianShellOrbitalRepresentation3D,
    right::CartesianGaussianShellOrbitalRepresentation3D,
)
    term_map = Dict{Tuple{Float64,Int,Int,Int},Float64}()
    for right_primitive in eachindex(right.exponents), left_primitive in eachindex(left.exponents)
        alpha_left = Float64(left.exponents[left_primitive])
        alpha_right = Float64(right.exponents[right_primitive])
        beta = alpha_left + alpha_right
        powers = (
            left.angular_powers[1] + right.angular_powers[1],
            left.angular_powers[2] + right.angular_powers[2],
            left.angular_powers[3] + right.angular_powers[3],
        )
        coefficient =
            Float64(left.coefficients[left_primitive]) *
            Float64(right.coefficients[right_primitive])
        for axis in 1:3
            coefficient *= GaussianAnalyticIntegrals.polynomial_gaussian_shell_prefactor(
                alpha_left,
                left.angular_powers[axis],
            )
            coefficient *= GaussianAnalyticIntegrals.polynomial_gaussian_shell_prefactor(
                alpha_right,
                right.angular_powers[axis],
            )
        end
        key = (beta, powers[1], powers[2], powers[3])
        term_map[key] = get(term_map, key, 0.0) + coefficient
    end

    entries = collect(term_map)
    sort!(
        entries;
        by = entry -> (entry.first[1], entry.first[2], entry.first[3], entry.first[4]),
    )
    terms = _GaussianCenteredPairTerm3D[]
    sizehint!(terms, length(entries))
    for (key, coefficient) in entries
        coefficient == 0.0 && continue
        push!(
            terms,
            _GaussianCenteredPairTerm3D(
                coefficient,
                key[1],
                (key[2], key[3], key[4]),
            ),
        )
    end
    return terms
end

function _gaussian_centered_global_term_coefficients(
    compact_terms::Vector{Vector{_GaussianCenteredPairTerm3D}},
)
    term_index_by_key = Dict{Tuple{Float64,Int,Int,Int},Int}()
    term_descriptors = _GaussianCenteredTermDescriptor3D[]
    compact_coefficients = Vector{Vector{_GaussianCenteredPairCoefficient}}(
        undef,
        length(compact_terms),
    )
    for compact_index in eachindex(compact_terms)
        coefficients = _GaussianCenteredPairCoefficient[]
        sizehint!(coefficients, length(compact_terms[compact_index]))
        for term in compact_terms[compact_index]
            key = (term.beta, term.powers[1], term.powers[2], term.powers[3])
            term_index = get(term_index_by_key, key, 0)
            if term_index == 0
                push!(
                    term_descriptors,
                    _GaussianCenteredTermDescriptor3D(term.beta, term.powers),
                )
                term_index = length(term_descriptors)
                term_index_by_key[key] = term_index
            end
            push!(
                coefficients,
                _GaussianCenteredPairCoefficient(term_index, term.coefficient),
            )
        end
        compact_coefficients[compact_index] = coefficients
    end
    return compact_coefficients, term_descriptors
end

function _gaussian_centered_term_kernel(
    term_descriptors::Vector{_GaussianCenteredTermDescriptor3D},
    expansion::CoulombGaussianExpansion,
)
    axis_descriptors, term_axis_indices =
        _gaussian_centered_axis_descriptors(term_descriptors)
    axis_kernel = zeros(Float64, length(axis_descriptors), length(axis_descriptors))
    term_count = length(term_descriptors)
    kernel = zeros(Float64, term_count, term_count)
    for expansion_index in eachindex(expansion.coefficients)
        coupling_exponent = Float64(expansion.exponents[expansion_index])
        coefficient = Float64(expansion.coefficients[expansion_index])
        for column in eachindex(axis_descriptors), row in 1:column
            left = axis_descriptors[row]
            right = axis_descriptors[column]
            value = GaussianAnalyticIntegrals.centered_polynomial_gaussian_pair_factor_integral(
                left[1],
                left[2],
                right[1],
                right[2],
                coupling_exponent,
            )
            axis_kernel[row, column] = value
            axis_kernel[column, row] = value
        end
        for column in 1:term_count, row in 1:column
            value =
                coefficient *
                axis_kernel[term_axis_indices[row, 1], term_axis_indices[column, 1]] *
                axis_kernel[term_axis_indices[row, 2], term_axis_indices[column, 2]] *
                axis_kernel[term_axis_indices[row, 3], term_axis_indices[column, 3]]
            kernel[row, column] += value
            row == column || (kernel[column, row] += value)
        end
    end
    return kernel
end

function _gaussian_centered_axis_descriptors(
    term_descriptors::Vector{_GaussianCenteredTermDescriptor3D},
)
    axis_index_by_key = Dict{Tuple{Float64,Int},Int}()
    axis_descriptors = Tuple{Float64,Int}[]
    term_axis_indices = Matrix{Int}(undef, length(term_descriptors), 3)
    for term_index in eachindex(term_descriptors)
        term = term_descriptors[term_index]
        for axis in 1:3
            key = (term.beta, term.powers[axis])
            axis_index = get(axis_index_by_key, key, 0)
            if axis_index == 0
                push!(axis_descriptors, key)
                axis_index = length(axis_descriptors)
                axis_index_by_key[key] = axis_index
            end
            term_axis_indices[term_index, axis] = axis_index
        end
    end
    return axis_descriptors, term_axis_indices
end

function _gaussian_centered_pair_integral(
    left_terms::Vector{_GaussianCenteredPairCoefficient},
    right_terms::Vector{_GaussianCenteredPairCoefficient},
    term_kernel::Matrix{Float64},
)
    value = 0.0
    for left in left_terms, right in right_terms
        value +=
            left.coefficient *
            right.coefficient *
            term_kernel[left.term_index, right.term_index]
    end
    return value
end
