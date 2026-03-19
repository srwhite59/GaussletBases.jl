"""
    MappedOrdinaryOneBody1D

Bundle of one-dimensional one-body ingredients for the mapped ordinary branch.

The object stores:

- the source basis or hybrid basis
- the backend label
- the overlap matrix
- the kinetic matrix
- Gaussianized one-body factor matrices for the requested exponents

It is the common input object for the current mapped Cartesian hydrogen,
harmonic-oscillator, and ordinary Cartesian IDA helpers.
"""
struct MappedOrdinaryOneBody1D{B}
    basis::B
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
    if operators.backend != :numerical_reference
        print(io, ", experimental=true")
    end
    print(io, ")")
end

function _basis_sample_matrix(
    basis_like,
    points::AbstractVector{<:Real},
)
    primitive_values = _primitive_sample_matrix(
        primitive_set(basis_like),
        Float64[Float64(point) for point in points],
    )
    return primitive_values * Matrix{Float64}(stencil_matrix(basis_like))
end

function _localized_reference_data(basis::MappedUniformBasis)
    representation = basis_representation(basis; operators = (:overlap, :position, :kinetic))
    transform, center_values = _cleanup_comx_transform(
        representation.basis_matrices.overlap,
        representation.basis_matrices.position,
        integral_weights(basis),
    )
    kinetic = Matrix{Float64}(transpose(transform) * representation.basis_matrices.kinetic * transform)
    return (
        transform = transform,
        centers = center_values,
        kinetic = kinetic,
    )
end

function _localized_alignment_transform(
    basis::MappedUniformBasis,
    reference_transform::AbstractMatrix{<:Real},
    localized::MappedPGDGLocalized1D;
    h::Real = 0.02,
)
    xlo, xhi = _primitive_set_bounds(primitive_set(basis))
    points, weights = _make_midpoint_grid(xlo, xhi, Float64(h))
    reference_values = _basis_sample_matrix(basis, points) * Matrix{Float64}(reference_transform)
    localized_values = _basis_sample_matrix(localized, points)
    cross_overlap = transpose(reference_values) * (weights .* localized_values)
    decomposition = svd(cross_overlap)
    return Matrix{Float64}(decomposition.U * decomposition.Vt)
end

function _localized_corrected_kinetic(
    basis::MappedUniformBasis,
    localized::MappedPGDGLocalized1D,
)
    reference = _localized_reference_data(basis)
    alignment = _localized_alignment_transform(basis, reference.transform, localized)
    return Matrix{Float64}(transpose(alignment) * reference.kinetic * alignment)
end

function _mapped_ordinary_backend_layer(
    basis::MappedUniformBasis,
    backend::Symbol,
)
    if backend == :numerical_reference
        return basis
    elseif backend == :pgdg_experimental
        return mapped_pgdg_logfit_prototype(basis)
    elseif backend == :pgdg_localized_experimental
        return mapped_pgdg_localized(mapped_pgdg_derivativefit_prototype(basis))
    else
        throw(ArgumentError("mapped ordinary backend must be :numerical_reference, :pgdg_experimental, or :pgdg_localized_experimental"))
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
- `:pgdg_experimental` uses the refined pre-COMX analytic PGDG-style proxy
- `:pgdg_localized_experimental` uses the refined proxy followed by overlap
  cleanup and COMX-style localization, using the current one-Gaussian
  derivative-aware proxy

The experimental backends are intended for mild-to-moderate distortion
regimes, with the numerical path retained as the reference route.
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

function mapped_ordinary_one_body_operators(
    basis::HybridMappedOrdinaryBasis1D;
    exponents::AbstractVector{<:Real} = Float64[],
    center::Real = 0.0,
)
    exponents_value = Float64[Float64(exponent) for exponent in exponents]
    center_value = Float64(center)
    overlap = Matrix{Float64}(overlap_matrix(basis))
    kinetic = Matrix{Float64}(kinetic_matrix(basis))
    gaussian_factors = Matrix{Float64}[
        Matrix{Float64}(factor) for factor in gaussian_factor_matrices(
            basis;
            exponents = exponents_value,
            center = center_value,
        )
    ]
    return MappedOrdinaryOneBody1D(
        basis,
        basis.backend,
        overlap,
        kinetic,
        gaussian_factors,
        exponents_value,
        center_value,
    )
end

function _mapped_ordinary_localized_oracle_operators(
    basis::MappedUniformBasis;
    exponents::AbstractVector{<:Real} = Float64[],
    center::Real = 0.0,
)
    localized = _mapped_ordinary_backend_layer(basis, :pgdg_localized_experimental)
    exponents_value = Float64[Float64(exponent) for exponent in exponents]
    center_value = Float64(center)
    overlap = Matrix{Float64}(overlap_matrix(localized))
    kinetic = _localized_corrected_kinetic(basis, localized)
    gaussian_factors = Matrix{Float64}[
        Matrix{Float64}(factor) for factor in gaussian_factor_matrices(
            localized;
            exponents = exponents_value,
            center = center_value,
        )
    ]
    return MappedOrdinaryOneBody1D(
        basis,
        :pgdg_localized_oracle,
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
    basis::HybridMappedOrdinaryBasis1D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 1.0,
)
    operators = mapped_ordinary_one_body_operators(
        basis;
        exponents = expansion.exponents,
        center = 0.0,
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
