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

struct _MappedOrdinaryGausslet1DBundle{B,L}
    basis::B
    layer::L
    backend::Symbol
    overlap::Matrix{Float64}
    kinetic::Matrix{Float64}
    position::Matrix{Float64}
    x2::Matrix{Float64}
    gaussian_factors::Vector{Matrix{Float64}}
    gaussian_factor_terms::Array{Float64,3}
    pair_factors::Vector{Matrix{Float64}}
    pair_factor_terms::Array{Float64,3}
    exponents::Vector{Float64}
    center::Float64
    weights::Vector{Float64}
    centers::Vector{Float64}
end

struct _MappedLegacyProxyLayer1D
    primitive_layer::PrimitiveSet1D
    coefficient_matrix::Matrix{Float64}
    center_data::Vector{Float64}
    integral_weight_data::Vector{Float64}
end

primitive_set(layer::_MappedLegacyProxyLayer1D) = layer.primitive_layer
stencil_matrix(layer::_MappedLegacyProxyLayer1D) = layer.coefficient_matrix
centers(layer::_MappedLegacyProxyLayer1D) = layer.center_data
integral_weights(layer::_MappedLegacyProxyLayer1D) = layer.integral_weight_data

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

function _term_tensor(matrices::AbstractVector{<:AbstractMatrix{<:Real}})
    isempty(matrices) && return zeros(Float64, 0, 0, 0)
    nrow, ncol = size(first(matrices))
    tensor = zeros(Float64, length(matrices), nrow, ncol)
    for term in eachindex(matrices)
        matrix = matrices[term]
        size(matrix) == (nrow, ncol) || throw(
            ArgumentError("term-first tensor conversion requires matrices of uniform size"),
        )
        @views tensor[term, :, :] .= Matrix{Float64}(matrix)
    end
    return tensor
end

function _tensor_to_matrix_vector(tensor::AbstractArray{<:Real,3})
    nterms = size(tensor, 1)
    matrices = Matrix{Float64}[]
    for term in 1:nterms
        push!(matrices, Matrix{Float64}(tensor[term, :, :]))
    end
    return matrices
end

function _basis_space_scalar_factor_data(
    basis_like;
    exponents::AbstractVector{<:Real},
    center::Real = 0.0,
    h = nothing,
)
    primitive_layer = primitive_set(basis_like)
    xlo, xhi = _primitive_set_bounds(primitive_layer)
    h_try = h === nothing ? _primitive_matrix_start_h(primitive_layer) : Float64(h)
    h_try > 0.0 || throw(ArgumentError("basis-space scalar factor construction requires h > 0"))

    exponent_values = Float64[Float64(exponent) for exponent in exponents]
    center_value = Float64(center)
    previous_x2 = nothing
    previous_terms = nothing
    current_x2 = nothing
    current_terms = zeros(Float64, length(exponent_values), length(basis_like), length(basis_like))

    for _ in 1:_PRIMITIVE_MATRIX_MAXITER
        points, weights = _make_midpoint_grid(xlo, xhi, h_try)
        basis_values = _basis_sample_matrix(basis_like, points)
        weighted_basis = weights .* basis_values
        current_x2 = Matrix{Float64}(transpose(basis_values) * ((weights .* (points .^ 2)) .* basis_values))
        current_x2 = 0.5 .* (current_x2 .+ transpose(current_x2))

        point_distance2 = (points .- center_value) .^ 2
        current_terms = zeros(Float64, length(exponent_values), size(current_x2, 1), size(current_x2, 2))
        for term in eachindex(exponent_values)
            factor_values = exp.(-exponent_values[term] .* point_distance2)
            matrix = Matrix{Float64}(transpose(basis_values) * ((weights .* factor_values) .* basis_values))
            @views current_terms[term, :, :] .= 0.5 .* (matrix .+ transpose(matrix))
        end

        if previous_x2 !== nothing
            maxdiff = maximum(
                vcat(
                    norm(current_x2 - previous_x2, Inf),
                    [norm(current_terms[term, :, :] - previous_terms[term, :, :], Inf) for term in eachindex(exponent_values)]...,
                ),
            )
            maxdiff <= _PRIMITIVE_MATRIX_TOL && return (
                x2 = current_x2,
                factor_terms = current_terms,
            )
        end

        previous_x2 = current_x2
        previous_terms = current_terms
        h_try /= 2.0
    end

    return (
        x2 = current_x2,
        factor_terms = current_terms,
    )
end

function _basis_space_pair_factor_terms(
    basis_like;
    exponents::AbstractVector{<:Real},
    h = nothing,
)
    primitive_layer = primitive_set(basis_like)
    xlo, xhi = _primitive_set_bounds(primitive_layer)
    h_try = h === nothing ? _primitive_matrix_start_h(primitive_layer) : Float64(h)
    h_try > 0.0 || throw(ArgumentError("basis-space pair-factor construction requires h > 0"))

    exponent_values = Float64[Float64(exponent) for exponent in exponents]
    previous_terms = nothing
    current_terms = zeros(Float64, length(exponent_values), length(basis_like), length(basis_like))

    for _ in 1:_PRIMITIVE_MATRIX_MAXITER
        points, weights = _make_midpoint_grid(xlo, xhi, h_try)
        basis_values = _basis_sample_matrix(basis_like, points)
        weighted_basis = weights .* basis_values
        distance2 = (points .- transpose(points)) .^ 2
        current_terms = zeros(Float64, length(exponent_values), size(basis_values, 2), size(basis_values, 2))
        for term in eachindex(exponent_values)
            kernel = exp.(-exponent_values[term] .* distance2)
            matrix = Matrix{Float64}(transpose(weighted_basis) * (kernel * weighted_basis))
            @views current_terms[term, :, :] .= 0.5 .* (matrix .+ transpose(matrix))
        end

        if previous_terms !== nothing
            maxdiff = maximum(
                norm(current_terms[term, :, :] - previous_terms[term, :, :], Inf) for term in eachindex(exponent_values)
            )
            maxdiff <= _PRIMITIVE_MATRIX_TOL && return current_terms
        end

        previous_terms = current_terms
        h_try /= 2.0
    end

    return current_terms
end

function _mapped_legacy_proxy_layer(
    basis::MappedUniformBasis,
)
    coefficient_matrix = Matrix{Float64}(stencil_matrix(basis))
    proxy_coefficients = similar(coefficient_matrix)
    proxy_primitives = AbstractPrimitiveFunction1D[]

    for (row, primitive) in enumerate(primitives(basis))
        if primitive isa Distorted{<:Gaussian}
            x0 = center(primitive)
            du = dudx(primitive.mapping, x0)
            du > 0.0 || throw(
                ArgumentError("legacy-style mapped proxy layer requires positive local Jacobian"),
            )
            width = primitive.primitive.width / du
            push!(proxy_primitives, Gaussian(center = x0, width = width))
            @views proxy_coefficients[row, :] .= sqrt(du) .* coefficient_matrix[row, :]
        elseif primitive isa Gaussian
            push!(proxy_primitives, primitive)
            @views proxy_coefficients[row, :] .= coefficient_matrix[row, :]
        else
            throw(
                ArgumentError(
                    "legacy-style mapped proxy layer currently requires Gaussian or Distorted{Gaussian} primitives",
                ),
            )
        end
    end

    primitive_layer = PrimitiveSet1D(
        proxy_primitives;
        name = :mapped_legacy_proxy_gaussians,
    )
    primitive_weights = Float64[integral_weight(primitive) for primitive in proxy_primitives]
    integral_weight_data = vec(transpose(proxy_coefficients) * primitive_weights)
    return _MappedLegacyProxyLayer1D(
        primitive_layer,
        proxy_coefficients,
        Float64[Float64(center_value) for center_value in centers(basis)],
        integral_weight_data,
    )
end

function _mapped_coulomb_expanded_symmetric_matrix(
    coefficients::AbstractVector{<:Real},
    x_terms::AbstractArray{<:Real,3},
    y_terms::AbstractArray{<:Real,3},
    z_terms::AbstractArray{<:Real,3},
)
    nterms = length(coefficients)
    size(x_terms, 1) == nterms == size(y_terms, 1) == size(z_terms, 1) || throw(
        ArgumentError("term-first Coulomb-expanded assembly requires matching term dimensions"),
    )
    n1x = size(x_terms, 2)
    n1y = size(y_terms, 2)
    n1z = size(z_terms, 2)
    orbitals = _mapped_cartesian_orbitals(1:n1x)
    length(orbitals) == n1x * n1y * n1z || throw(
        ArgumentError("term-first Coulomb-expanded assembly currently expects equal 1D dimensions on each axis"),
    )
    matrix = zeros(Float64, length(orbitals), length(orbitals))
    coeffs = Float64[Float64(value) for value in coefficients]

    for i in eachindex(orbitals)
        oi = orbitals[i]
        for j in 1:i
            oj = orbitals[j]
            value = 0.0
            @inbounds for term in 1:nterms
                value += coeffs[term] *
                    x_terms[term, oi.ix, oj.ix] *
                    y_terms[term, oi.iy, oj.iy] *
                    z_terms[term, oi.iz, oj.iz]
            end
            matrix[i, j] = value
            matrix[j, i] = value
        end
    end

    return matrix
end

function _mapped_ordinary_gausslet_1d_bundle(
    basis::MappedUniformBasis;
    exponents::AbstractVector{<:Real} = Float64[],
    center::Real = 0.0,
    backend::Symbol = :pgdg_experimental,
)
    exponents_value = Float64[Float64(exponent) for exponent in exponents]
    center_value = Float64(center)
    layer = _mapped_ordinary_backend_layer(basis, backend)

    if backend == :numerical_reference
        representation = basis_representation(basis; operators = (:overlap, :position, :kinetic))
        overlap = Matrix{Float64}(representation.basis_matrices.overlap)
        position = Matrix{Float64}(representation.basis_matrices.position)
        kinetic = Matrix{Float64}(representation.basis_matrices.kinetic)
        scalar_data = _basis_space_scalar_factor_data(
            basis;
            exponents = exponents_value,
            center = center_value,
        )
        x2 = scalar_data.x2
        gaussian_factor_terms = scalar_data.factor_terms
        # Follow the legacy ordinary pattern more closely here: build the
        # Coulomb-expansion pair terms from a mapped pure-Gaussian proxy raw
        # layer, then contract to the working 1D gausslet basis.
        proxy_layer = _mapped_legacy_proxy_layer(basis)
        pair_factors = Matrix{Float64}[
            Matrix{Float64}(factor) for factor in _pair_gaussian_factor_matrices(
                proxy_layer;
                exponents = exponents_value,
            )
        ]
        pair_factor_terms_raw = _term_tensor(pair_factors)
    else
        overlap = Matrix{Float64}(overlap_matrix(layer))
        position = Matrix{Float64}(position_matrix(layer))
        kinetic = Matrix{Float64}(kinetic_matrix(layer))
        x2 = Matrix{Float64}(_x2_matrix(layer))
        gaussian_factors = Matrix{Float64}[
            Matrix{Float64}(factor) for factor in gaussian_factor_matrices(
                layer;
                exponents = exponents_value,
                center = center_value,
            )
        ]
        gaussian_factor_terms = _term_tensor(gaussian_factors)
        pair_factors = Matrix{Float64}[
            Matrix{Float64}(factor) for factor in _pair_gaussian_factor_matrices(
                layer;
                exponents = exponents_value,
            )
        ]
        pair_factor_terms_raw = _term_tensor(pair_factors)
    end

    weight_values = Float64[Float64(weight) for weight in integral_weights(layer)]
    center_values = Float64[Float64(point) for point in centers(layer)]
    gaussian_factors = backend == :numerical_reference ? _tensor_to_matrix_vector(gaussian_factor_terms) : gaussian_factors
    weight_outer = weight_values * transpose(weight_values)
    pair_factors = [
        Matrix{Float64}(pair_factor_terms_raw[term, :, :] ./ weight_outer) for term in eachindex(exponents_value)
    ]
    pair_factor_terms = _term_tensor(pair_factors)

    return _MappedOrdinaryGausslet1DBundle(
        basis,
        layer,
        backend,
        overlap,
        kinetic,
        position,
        x2,
        gaussian_factors,
        gaussian_factor_terms,
        pair_factors,
        pair_factor_terms,
        exponents_value,
        center_value,
        weight_values,
        center_values,
    )
end

function _mapped_ordinary_one_body_from_bundle(
    bundle::_MappedOrdinaryGausslet1DBundle,
)
    return MappedOrdinaryOneBody1D(
        bundle.basis,
        bundle.backend,
        bundle.overlap,
        bundle.kinetic,
        bundle.gaussian_factors,
        bundle.exponents,
        bundle.center,
    )
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
