struct CartesianProductOrbital3D
    index::Int
    ix::Int
    iy::Int
    iz::Int
    x::Float64
    y::Float64
    z::Float64
end

function Base.show(io::IO, orbital::CartesianProductOrbital3D)
    print(
        io,
        "CartesianProductOrbital3D(index=",
        orbital.index,
        ", ix=",
        orbital.ix,
        ", iy=",
        orbital.iy,
        ", iz=",
        orbital.iz,
        ", center=(",
        orbital.x,
        ", ",
        orbital.y,
        ", ",
        orbital.z,
        "))",
    )
end

struct OrdinaryCartesianIDAOperators
    basis::MappedUniformBasis
    backend::Symbol
    expansion::CoulombGaussianExpansion
    one_body_1d::MappedOrdinaryOneBody1D
    overlap_3d::Matrix{Float64}
    one_body_hamiltonian::Matrix{Float64}
    pair_factors_1d::Vector{Matrix{Float64}}
    interaction_matrix::Matrix{Float64}
    orbital_data::Vector{CartesianProductOrbital3D}
    weight_1d::Vector{Float64}
    weight_3d::Vector{Float64}
end

function Base.show(io::IO, operators::OrdinaryCartesianIDAOperators)
    print(
        io,
        "OrdinaryCartesianIDAOperators(backend=:",
        operators.backend,
        ", n1d=",
        length(operators.weight_1d),
        ", norbitals=",
        length(operators.orbital_data),
        ", nterms=",
        length(operators.pair_factors_1d),
    )
    if operators.backend != :numerical_reference
        print(io, ", experimental=true")
    end
    print(io, ")")
end

orbitals(operators::OrdinaryCartesianIDAOperators) = operators.orbital_data

function _gaussian_pair_factor(a::Gaussian, b::Gaussian, exponent::Float64)
    exponent >= 0.0 || throw(ArgumentError("pair-factor exponent must be >= 0"))

    alpha_a = inv(a.width^2)
    alpha_b = inv(b.width^2)
    a11 = alpha_a + 2.0 * exponent
    a22 = alpha_b + 2.0 * exponent
    a12 = -2.0 * exponent
    determinant = a11 * a22 - a12^2
    determinant > 0.0 || throw(ArgumentError("pair-factor quadratic form must be positive definite"))

    d1 = alpha_a * a.center_value
    d2 = alpha_b * b.center_value
    constant_term = alpha_a * a.center_value^2 + alpha_b * b.center_value^2
    quadratic_term = (a22 * d1^2 - 2.0 * a12 * d1 * d2 + a11 * d2^2) / determinant
    return (2.0 * pi / sqrt(determinant)) * exp(-0.5 * (constant_term - quadratic_term))
end

function _primitive_pair_gaussian_factor_matrix(
    set::PrimitiveSet1D,
    ::_AnalyticPrimitiveMatrixBackend;
    exponent::Float64,
)
    matrix = zeros(Float64, length(set), length(set))
    for a in 1:length(set)
        pa = primitives(set)[a]
        for b in a:length(set)
            pb = primitives(set)[b]
            value_ab = _gaussian_pair_factor(pa, pb, exponent)
            matrix[a, b] = value_ab
            matrix[b, a] = value_ab
        end
    end
    return matrix
end

function _primitive_pair_gaussian_factor_matrix(
    set::PrimitiveSet1D,
    ::_NumericalPrimitiveMatrixBackend;
    exponent::Float64,
    h = nothing,
)
    exponent >= 0.0 || throw(ArgumentError("pair-factor exponent must be >= 0"))

    xlo, xhi = _primitive_set_bounds(set)
    h_try = h === nothing ? _primitive_matrix_start_h(set) : Float64(h)
    h_try > 0.0 || throw(ArgumentError("numerical pair-factor matrix requires h > 0"))

    previous = nothing
    current = nothing
    for _ in 1:_PRIMITIVE_MATRIX_MAXITER
        points, weights = _make_midpoint_grid(xlo, xhi, h_try)
        values = _primitive_sample_matrix(set, points)
        weighted_values = weights .* values
        kernel = exp.(-exponent .* ((points .- transpose(points)) .^ 2))
        current = _symmetrize_primitive_matrix(transpose(weighted_values) * (kernel * weighted_values))
        if previous !== nothing && norm(current - previous, Inf) <= _PRIMITIVE_MATRIX_TOL
            return current
        end
        previous = current
        h_try /= 2.0
    end
    return current
end

function _primitive_pair_gaussian_factor_matrices(
    set::PrimitiveSet1D,
    backend;
    exponents::AbstractVector{<:Real},
)
    return [
        _primitive_pair_gaussian_factor_matrix(
            set,
            backend;
            exponent = Float64(exponent),
        ) for exponent in exponents
    ]
end

function _primitive_pair_gaussian_factor_matrices(
    set::PrimitiveSet1D,
    ::_NumericalPrimitiveMatrixBackend;
    exponents::AbstractVector{<:Real},
    h = nothing,
)
    xlo, xhi = _primitive_set_bounds(set)
    h_try = h === nothing ? _primitive_matrix_start_h(set) : Float64(h)
    h_try > 0.0 || throw(ArgumentError("numerical pair-factor matrices require h > 0"))

    previous = nothing
    current = nothing
    exponent_values = Float64[Float64(exponent) for exponent in exponents]
    for _ in 1:_PRIMITIVE_MATRIX_MAXITER
        points, weights = _make_midpoint_grid(xlo, xhi, h_try)
        values = _primitive_sample_matrix(set, points)
        weighted_values = weights .* values
        current = Matrix{Float64}[]
        for exponent in exponent_values
            kernel = exp.(-exponent .* ((points .- transpose(points)) .^ 2))
            matrix = _symmetrize_primitive_matrix(transpose(weighted_values) * (kernel * weighted_values))
            push!(current, matrix)
        end
        if previous !== nothing
            maxdiff = maximum(norm(current[index] - previous[index], Inf) for index in eachindex(current))
            maxdiff <= _PRIMITIVE_MATRIX_TOL && return current
        end
        previous = current
        h_try /= 2.0
    end
    return current
end

function _pair_gaussian_factor_matrices(
    layer;
    exponents::AbstractVector{<:Real},
)
    primitive_layer = primitive_set(layer)
    primitive_matrices = _primitive_pair_gaussian_factor_matrices(
        primitive_layer,
        _select_primitive_matrix_backend(primitive_layer);
        exponents = exponents,
    )
    return [Matrix{Float64}(contract_primitive_matrix(layer, matrix)) for matrix in primitive_matrices]
end

function _mapped_cartesian_orbitals(centers_1d::AbstractVector{<:Real})
    orbitals_out = CartesianProductOrbital3D[]
    index = 0
    for ix in eachindex(centers_1d), iy in eachindex(centers_1d), iz in eachindex(centers_1d)
        index += 1
        push!(
            orbitals_out,
            CartesianProductOrbital3D(
                index,
                ix,
                iy,
                iz,
                Float64(centers_1d[ix]),
                Float64(centers_1d[iy]),
                Float64(centers_1d[iz]),
            ),
        )
    end
    return orbitals_out
end

function _mapped_cartesian_weights(weight_1d::AbstractVector{<:Real})
    weights_out = Float64[]
    for wx in weight_1d, wy in weight_1d, wz in weight_1d
        push!(weights_out, Float64(wx) * Float64(wy) * Float64(wz))
    end
    return weights_out
end

function _mapped_cartesian_one_body_matrix(
    one_body::MappedOrdinaryOneBody1D,
    expansion::CoulombGaussianExpansion;
    Z::Real,
)
    length(one_body.gaussian_factors) == length(expansion) ||
        throw(ArgumentError("one-body operator set must carry one Gaussian factor per Coulomb term"))

    overlap = one_body.overlap
    kinetic = one_body.kinetic
    overlap_3d = kron(overlap, kron(overlap, overlap))
    hamiltonian =
        kron(kinetic, kron(overlap, overlap)) +
        kron(overlap, kron(kinetic, overlap)) +
        kron(overlap, kron(overlap, kinetic))

    for term in eachindex(expansion.coefficients)
        factor = one_body.gaussian_factors[term]
        hamiltonian .-= Float64(Z) * expansion.coefficients[term] .* kron(factor, kron(factor, factor))
    end

    return overlap_3d, hamiltonian
end

"""
    ordinary_cartesian_ida_operators(
        basis::MappedUniformBasis;
        expansion = coulomb_gaussian_expansion(doacc = false),
        Z = 2.0,
        backend = :pgdg_experimental,
    )

Build the current ordinary Cartesian static IDA ingredients on the raw
Cartesian product basis.

This is the next ordinary-branch milestone after the one-body hydrogen path:

- one global map on each Cartesian axis
- one-body ingredients from the mapped ordinary backend split
- separable electron-electron IDA factors built from the same Coulomb
  expansion

The result is a static object, not a He solver:

- `overlap_3d`
- `one_body_hamiltonian`
- `interaction_matrix`
- explicit product-orbital indexing

`backend = :pgdg_localized_experimental` is the candidate solver-ready
implementation route in the mild-to-moderate mapped regime. In the current
experimental implementation, that route uses the cleaned/localized PGDG-style
one-dimensional basis together with a more derivative-aware analytic primitive
proxy than the pre-COMX path.
`:pgdg_experimental` retains the pre-COMX refined proxy path.
`:numerical_reference` remains the validation route.
"""
function ordinary_cartesian_ida_operators(
    basis::MappedUniformBasis;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
    backend::Symbol = :pgdg_experimental,
)
    one_body = mapped_ordinary_one_body_operators(
        basis;
        exponents = expansion.exponents,
        center = 0.0,
        backend = backend,
    )
    layer = _mapped_ordinary_backend_layer(basis, backend)
    pair_factors_basis = _pair_gaussian_factor_matrices(layer; exponents = expansion.exponents)

    weight_1d = Float64[Float64(weight) for weight in integral_weights(layer)]
    any(weight -> abs(weight) <= 1.0e-12, weight_1d) &&
        throw(ArgumentError("ordinary_cartesian_ida_operators requires nonzero 1D basis weights"))

    weight_outer = weight_1d * transpose(weight_1d)
    pair_factors_1d = [factor ./ weight_outer for factor in pair_factors_basis]
    overlap_3d, one_body_hamiltonian = _mapped_cartesian_one_body_matrix(one_body, expansion; Z = Z)

    interaction_matrix = zeros(Float64, size(overlap_3d))
    for term in eachindex(expansion.coefficients)
        factor = pair_factors_1d[term]
        interaction_matrix .+= expansion.coefficients[term] .* kron(factor, kron(factor, factor))
    end

    return OrdinaryCartesianIDAOperators(
        basis,
        backend,
        expansion,
        one_body,
        overlap_3d,
        one_body_hamiltonian,
        pair_factors_1d,
        interaction_matrix,
        _mapped_cartesian_orbitals(centers(basis)),
        weight_1d,
        _mapped_cartesian_weights(weight_1d),
    )
end
