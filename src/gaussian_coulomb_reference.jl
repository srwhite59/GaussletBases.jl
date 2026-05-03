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

struct _GaussianCoulombPairTerm3D
    coefficient::Float64
    x::_GaussianCoulombAxisPairTerm
    y::_GaussianCoulombAxisPairTerm
    z::_GaussianCoulombAxisPairTerm
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

function _gaussian_coulomb_pair_terms(
    left::CartesianGaussianShellOrbitalRepresentation3D,
    right::CartesianGaussianShellOrbitalRepresentation3D,
)
    _validate_gaussian_coulomb_orbital(left)
    _validate_gaussian_coulomb_orbital(right)
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

function _gaussian_coulomb_axis_integral(
    left::_GaussianCoulombAxisPairTerm,
    right::_GaussianCoulombAxisPairTerm,
    coupling_exponent::Float64,
)
    return GaussianAnalyticIntegrals.polynomial_gaussian_pair_factor_integral(
        left.alpha_left,
        left.center_left,
        left.power_left,
        left.prefactor_left,
        left.alpha_right,
        left.center_right,
        left.power_right,
        left.prefactor_right,
        right.alpha_left,
        right.center_left,
        right.power_left,
        right.prefactor_left,
        right.alpha_right,
        right.center_right,
        right.power_right,
        right.prefactor_right,
        coupling_exponent,
    )
end
