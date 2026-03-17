struct MappedOrdinaryOneBody1D
    basis::MappedUniformBasis
    backend::Symbol
    overlap::Matrix{Float64}
    kinetic::Matrix{Float64}
    gaussian_factors::Vector{Matrix{Float64}}
    exponents::Vector{Float64}
    center::Float64
end

function Base.show(io::IO, operators::MappedOrdinaryOneBody1D)
    print(
        io,
        "MappedOrdinaryOneBody1D(backend=:",
        operators.backend,
        ", nbasis=",
        size(operators.overlap, 1),
        ", nfactors=",
        length(operators.gaussian_factors),
    )
    if operators.backend == :pgdg_experimental
        print(io, ", experimental=true")
    end
    print(io, ")")
end

function _mapped_ordinary_backend_layer(
    basis::MappedUniformBasis,
    backend::Symbol,
)
    if backend == :numerical_reference
        return basis
    elseif backend == :pgdg_experimental
        return mapped_pgdg_logfit_prototype(basis)
    else
        throw(ArgumentError("mapped ordinary backend must be :numerical_reference or :pgdg_experimental"))
    end
end

"""
    mapped_ordinary_one_body_operators(
        basis::MappedUniformBasis;
        exponents = Real[],
        center = 0.0,
        backend = :pgdg_experimental,
    )

Build the one-dimensional mapped ordinary one-body ingredients for a full-line
mapped basis.

The backend choice is explicit:

- `:numerical_reference` keeps the trusted mapped numerical path
- `:pgdg_experimental` uses the experimental analytic PGDG-style proxy

The experimental backend is intended for mild-to-moderate distortion regimes,
with the numerical path retained as the reference route.
"""
function mapped_ordinary_one_body_operators(
    basis::MappedUniformBasis;
    exponents::AbstractVector{<:Real} = Float64[],
    center::Real = 0.0,
    backend::Symbol = :pgdg_experimental,
)
    exponents_value = Float64[Float64(exponent) for exponent in exponents]
    center_value = Float64(center)
    layer = _mapped_ordinary_backend_layer(basis, backend)

    if backend == :numerical_reference
        representation = basis_representation(basis; operators = (:overlap, :kinetic))
        overlap = Matrix{Float64}(representation.basis_matrices.overlap)
        kinetic = Matrix{Float64}(representation.basis_matrices.kinetic)
    else
        overlap = Matrix{Float64}(overlap_matrix(layer))
        kinetic = Matrix{Float64}(kinetic_matrix(layer))
    end

    gaussian_factors = Matrix{Float64}[
        Matrix{Float64}(factor) for factor in gaussian_factor_matrices(
            layer;
            exponents = exponents_value,
            center = center_value,
        )
    ]

    return MappedOrdinaryOneBody1D(
        basis,
        backend,
        overlap,
        kinetic,
        gaussian_factors,
        exponents_value,
        center_value,
    )
end

function _orthonormalize_cartesian_1d(
    overlap::AbstractMatrix{<:Real},
    operators::AbstractVector{<:AbstractMatrix{<:Real}},
)
    decomposition = eigen(Symmetric(Matrix{Float64}(overlap)))
    invhalf =
        decomposition.vectors *
        Diagonal(1.0 ./ sqrt.(decomposition.values)) *
        transpose(decomposition.vectors)
    transformed = [invhalf * Matrix{Float64}(operator) * invhalf for operator in operators]
    return transformed
end

function _mapped_cartesian_hydrogen_energy(
    operators::MappedOrdinaryOneBody1D,
    expansion::CoulombGaussianExpansion;
    Z::Real = 1.0,
)
    length(operators.gaussian_factors) == length(expansion) ||
        throw(ArgumentError("mapped_cartesian_hydrogen_energy requires one Gaussian factor matrix per Coulomb-expansion term"))

    transformed = _orthonormalize_cartesian_1d(
        operators.overlap,
        [operators.kinetic, operators.gaussian_factors...],
    )
    kinetic_orth = first(transformed)
    gaussian_orth = transformed[2:end]
    n1d = size(operators.overlap, 1)
    identity_1d = Matrix{Float64}(I, n1d, n1d)
    hamiltonian =
        kron(kinetic_orth, kron(identity_1d, identity_1d)) +
        kron(identity_1d, kron(kinetic_orth, identity_1d)) +
        kron(identity_1d, kron(identity_1d, kinetic_orth))

    for term in eachindex(expansion.coefficients)
        factor = gaussian_orth[term]
        hamiltonian .-= Float64(Z) * expansion.coefficients[term] .* kron(factor, kron(factor, factor))
    end

    return minimum(eigen(Hermitian(hamiltonian)).values)
end

"""
    mapped_cartesian_hydrogen_energy(
        basis::MappedUniformBasis;
        expansion = coulomb_gaussian_expansion(doacc = false),
        Z = 1.0,
        backend = :pgdg_experimental,
    )

Build the current mapped Cartesian hydrogen one-body Hamiltonian with the
chosen backend and return its ground-state energy.
"""
function mapped_cartesian_hydrogen_energy(
    basis::MappedUniformBasis;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 1.0,
    backend::Symbol = :pgdg_experimental,
)
    operators = mapped_ordinary_one_body_operators(
        basis;
        exponents = expansion.exponents,
        center = 0.0,
        backend = backend,
    )
    return _mapped_cartesian_hydrogen_energy(operators, expansion; Z = Z)
end

function mapped_cartesian_hydrogen_energy(
    operators::MappedOrdinaryOneBody1D,
    expansion::CoulombGaussianExpansion;
    Z::Real = 1.0,
)
    return _mapped_cartesian_hydrogen_energy(operators, expansion; Z = Z)
end
