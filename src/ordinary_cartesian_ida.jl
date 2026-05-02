"""
    CartesianProductOrbital3D

One orbital index in the current ordinary Cartesian product basis.

The object records the Cartesian product indices `(ix, iy, iz)` together with
the corresponding product-basis center `(x, y, z)`.
"""
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

"""
    OrdinaryCartesianIDAOperators

Static ordinary Cartesian one-body and density-density / IDA interaction data
built on the current mapped or hybrid one-dimensional ordinary basis.

The object bundles:

- the one-dimensional one-body ingredients used to build it
- the full Cartesian-product overlap and one-body matrices
- the separable pair factors used in the current Coulomb-expansion assembly
- the dense two-index interaction matrix
- the interaction treatment label used for the current hybrid ordinary branch
- explicit product-orbital indexing metadata

It is a solver-facing static Hamiltonian object, not a solver itself.
"""
struct OrdinaryCartesianIDAOperators{B}
    basis::B
    backend::Symbol
    interaction_treatment::Symbol
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
        ", interaction=:",
        operators.interaction_treatment,
    )
    if operators.backend != :numerical_reference
        print(io, ", experimental=true")
    end
    print(io, ")")
end

orbitals(operators::OrdinaryCartesianIDAOperators) = operators.orbital_data

"""
    ordinary_cartesian_vee_expectation(
        operators::OrdinaryCartesianIDAOperators,
        orbital::AbstractVector;
        overlap_tol = 1.0e-8,
    )

Return the current density-density / IDA interaction expectation value for a
single spatial orbital occupied once with spin up and once with spin down in
the ordinary Cartesian branch.

This helper is intentionally narrow. It assumes the Cartesian product basis is
already orthonormal to within `overlap_tol`, which is the intended pure
ordinary validation regime. It is therefore suitable for the current identity-
mapped ordinary-gausslet `1s^2` check, but it is not a general nonorthogonal-
basis two-electron expectation evaluator.
"""
function ordinary_cartesian_vee_expectation(
    operators::OrdinaryCartesianIDAOperators,
    orbital::AbstractVector;
    overlap_tol::Real = 1.0e-8,
)
    length(orbital) == length(operators.orbital_data) ||
        throw(ArgumentError("orbital length must match the ordinary Cartesian orbital dimension"))
    overlap_error = norm(operators.overlap_3d - I, Inf)
    overlap_error <= Float64(overlap_tol) || throw(
        ArgumentError(
            "ordinary_cartesian_vee_expectation currently requires an orthonormal Cartesian product basis; got overlap error $(overlap_error)",
        ),
    )
    weights = Float64[abs2(coefficient) for coefficient in orbital]
    norm2 = sum(weights)
    norm2 > 0.0 || throw(ArgumentError("orbital must have nonzero norm"))
    weights ./= norm2
    return Float64(real(dot(weights, operators.interaction_matrix * weights)))
end

"""
    ordinary_cartesian_1s2_check(
        operators::OrdinaryCartesianIDAOperators;
        overlap_tol = 1.0e-8,
    )

Return a small validation package for the doubly occupied noninteracting
`1s`-style reference state on the current ordinary Cartesian IDA object.

The helper diagonalizes the one-body Hamiltonian in the current basis, takes
the lowest orbital, and evaluates its doubly occupied density-density / IDA
interaction expectation through [`ordinary_cartesian_vee_expectation`](@ref).
Like that lower-level helper, it is intended for the orthonormal pure-
ordinary and localized-hybrid validation routes, not for a generic
nonorthogonal-basis two-electron solve.
"""
function ordinary_cartesian_1s2_check(
    operators::OrdinaryCartesianIDAOperators;
    overlap_tol::Real = 1.0e-8,
)
    decomposition = eigen(Hermitian(operators.one_body_hamiltonian))
    orbital = decomposition.vectors[:, 1]
    return (
        orbital_energy = Float64(decomposition.values[1]),
        orbital = orbital,
        vee_expectation = ordinary_cartesian_vee_expectation(
            operators,
            orbital;
            overlap_tol = overlap_tol,
        ),
        overlap_error = norm(operators.overlap_3d - I, Inf),
    )
end

struct _AuxiliaryContractedLayer1D
    primitive_layer::PrimitiveSet1D
    coefficient_matrix::Matrix{Float64}
    center_data::Vector{Float64}
    integral_weight_data::Vector{Float64}
end

primitive_set(layer::_AuxiliaryContractedLayer1D) = layer.primitive_layer
stencil_matrix(layer::_AuxiliaryContractedLayer1D) = layer.coefficient_matrix
centers(layer::_AuxiliaryContractedLayer1D) = layer.center_data
integral_weights(layer::_AuxiliaryContractedLayer1D) = layer.integral_weight_data

_gaussian_pair_factor(a::Gaussian, b::Gaussian, exponent::Float64) =
    GaussianAnalyticIntegrals.gaussian_pair_factor(a, b, exponent)

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

function _mapped_cartesian_orbitals(
    x_centers::AbstractVector{<:Real},
    y_centers::AbstractVector{<:Real},
    z_centers::AbstractVector{<:Real},
)
    orbitals_out = CartesianProductOrbital3D[]
    index = 0
    for ix in eachindex(x_centers), iy in eachindex(y_centers), iz in eachindex(z_centers)
        index += 1
        push!(
            orbitals_out,
            CartesianProductOrbital3D(
                index,
                ix,
                iy,
                iz,
                Float64(x_centers[ix]),
                Float64(y_centers[iy]),
                Float64(z_centers[iz]),
            ),
        )
    end
    return orbitals_out
end

function _mapped_cartesian_orbitals(centers_1d::AbstractVector{<:Real})
    return _mapped_cartesian_orbitals(centers_1d, centers_1d, centers_1d)
end

function _mapped_cartesian_weights(
    weight_x::AbstractVector{<:Real},
    weight_y::AbstractVector{<:Real},
    weight_z::AbstractVector{<:Real},
)
    weights_out = Float64[]
    for wx in weight_x, wy in weight_y, wz in weight_z
        push!(weights_out, Float64(wx) * Float64(wy) * Float64(wz))
    end
    return weights_out
end

function _mapped_cartesian_weights(weight_1d::AbstractVector{<:Real})
    return _mapped_cartesian_weights(weight_1d, weight_1d, weight_1d)
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

function _symmetrize_ida_matrix(matrix::AbstractMatrix{<:Real})
    matrix_value = Matrix{Float64}(matrix)
    return 0.5 .* (matrix_value .+ transpose(matrix_value))
end

function _normalized_density_transfer(overlaps::AbstractMatrix{<:Real})
    transfer = abs2.(Matrix{Float64}(overlaps))
    for column in axes(transfer, 2)
        total = sum(view(transfer, :, column))
        total > 0.0 || throw(ArgumentError("density transfer requires nonzero overlap into every target column"))
        transfer[:, column] ./= total
    end
    return transfer
end

function _nearest_density_transfer(
    source_centers::AbstractVector{<:Real},
    target_centers::AbstractVector{<:Real},
)
    transfer = zeros(Float64, length(source_centers), length(target_centers))
    for (column, target_center) in enumerate(target_centers)
        index = argmin(abs.(Float64.(source_centers) .- Float64(target_center)))
        transfer[index, column] = 1.0
    end
    return transfer
end

function _ordinary_cartesian_ida_from_pair_factors(
    basis,
    backend::Symbol,
    interaction_treatment::Symbol,
    expansion::CoulombGaussianExpansion,
    one_body::MappedOrdinaryOneBody1D,
    pair_factors_1d::Vector{Matrix{Float64}},
    one_body_hamiltonian::AbstractMatrix{<:Real},
    weight_1d::AbstractVector{<:Real},
)
    overlap_3d, _ = _mapped_cartesian_one_body_matrix(one_body, expansion; Z = 0.0)
    interaction_matrix = zeros(Float64, size(overlap_3d))
    for term in eachindex(expansion.coefficients)
        factor = pair_factors_1d[term]
        interaction_matrix .+= expansion.coefficients[term] .* kron(factor, kron(factor, factor))
    end

    weight_values = Float64[Float64(weight) for weight in weight_1d]

    return OrdinaryCartesianIDAOperators(
        basis,
        backend,
        interaction_treatment,
        expansion,
        one_body,
        overlap_3d,
        Matrix{Float64}(one_body_hamiltonian),
        pair_factors_1d,
        interaction_matrix,
        _mapped_cartesian_orbitals(centers(basis)),
        weight_values,
        _mapped_cartesian_weights(weight_values),
    )
end

function _ordinary_cartesian_ida_from_layer(
    basis,
    backend::Symbol,
    expansion::CoulombGaussianExpansion,
    one_body::MappedOrdinaryOneBody1D,
    layer;
    Z::Real,
)
    pair_factors_basis = _pair_gaussian_factor_matrices(layer; exponents = expansion.exponents)

    weight_1d = Float64[Float64(weight) for weight in integral_weights(layer)]
    any(weight -> abs(weight) <= 1.0e-12, weight_1d) &&
        throw(ArgumentError("ordinary Cartesian IDA layer requires nonzero 1D basis weights"))

    weight_outer = weight_1d * transpose(weight_1d)
    pair_factors_1d = [factor ./ weight_outer for factor in pair_factors_basis]
    _, one_body_hamiltonian = _mapped_cartesian_one_body_matrix(one_body, expansion; Z = Z)
    return _ordinary_cartesian_ida_from_pair_factors(
        basis,
        backend,
        :combined_basis,
        expansion,
        one_body,
        pair_factors_1d,
        one_body_hamiltonian,
        weight_1d,
    )
end

function _ordinary_cartesian_ida_from_gausslet_bundle(
    bundle::_MappedOrdinaryGausslet1DBundle,
    expansion::CoulombGaussianExpansion;
    Z::Real,
)
    one_body = _mapped_ordinary_one_body_from_bundle(bundle)
    pgdg_intermediate = bundle.pgdg_intermediate
    overlap_3d, _ = _mapped_cartesian_one_body_matrix(one_body, expansion; Z = 0.0)
    interaction_matrix = _mapped_coulomb_expanded_symmetric_matrix(
        expansion.coefficients,
        pgdg_intermediate.pair_factor_terms,
        pgdg_intermediate.pair_factor_terms,
        pgdg_intermediate.pair_factor_terms,
    )
    one_body_hamiltonian =
        _mapped_coulomb_expanded_symmetric_matrix(
            -Float64(Z) .* expansion.coefficients,
            pgdg_intermediate.gaussian_factor_terms,
            pgdg_intermediate.gaussian_factor_terms,
            pgdg_intermediate.gaussian_factor_terms,
        )
    n1d = size(pgdg_intermediate.overlap, 1)
    identity_1d = Matrix{Float64}(I, n1d, n1d)
    kinetic = pgdg_intermediate.kinetic
    one_body_hamiltonian .+= kron(kinetic, kron(identity_1d, identity_1d))
    one_body_hamiltonian .+= kron(identity_1d, kron(kinetic, identity_1d))
    one_body_hamiltonian .+= kron(identity_1d, kron(identity_1d, kinetic))

    return OrdinaryCartesianIDAOperators(
        bundle.basis,
        bundle.backend,
        :combined_basis,
        expansion,
        one_body,
        overlap_3d,
        one_body_hamiltonian,
        pgdg_intermediate.pair_factors,
        interaction_matrix,
        _mapped_cartesian_orbitals(pgdg_intermediate.centers),
        pgdg_intermediate.weights,
        _mapped_cartesian_weights(pgdg_intermediate.weights),
    )
end

# Alg QW-RG step 4: Define residual Gaussians by orthogonalizing the added Gaussian channel.
# See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
function _hybrid_residual_gaussian_seed(hybrid::HybridMappedOrdinaryBasis1D)
    backbone = _mapped_ordinary_backend_layer(hybrid.source_basis, hybrid.backend)
    nbackbone = length(backbone)
    ncore = length(hybrid.core_gaussians)
    ncore > 0 || throw(ArgumentError("residual Gaussian interaction treatment requires at least one added core Gaussian"))

    backbone_coefficients_small = Matrix{Float64}(stencil_matrix(backbone))
    hybrid_coefficients = Matrix{Float64}(stencil_matrix(hybrid))
    nbackbone_primitives = size(backbone_coefficients_small, 1)
    core_coefficients_small = Matrix{Float64}(hybrid.core_coefficient_matrix)
    ncore_primitives = size(core_coefficients_small, 1)
    nprimitives = size(hybrid_coefficients, 1)
    nprimitives == nbackbone_primitives + ncore_primitives || throw(
        ArgumentError("hybrid primitive layout must be backbone primitives followed by contracted core Gaussian primitives"),
    )
    size(core_coefficients_small, 2) == ncore || throw(
        ArgumentError("hybrid residual seed requires core contraction columns to match core supplement count"),
    )

    backbone_coefficients = zeros(Float64, nprimitives, nbackbone)
    backbone_coefficients[1:nbackbone_primitives, :] .= backbone_coefficients_small
    raw_core_coefficients = zeros(Float64, nprimitives, ncore)
    raw_core_coefficients[(nbackbone_primitives + 1):end, :] .= core_coefficients_small

    primitive_overlap = overlap_matrix(primitive_set(hybrid))
    overlap_backbone = Matrix{Float64}(backbone_coefficients' * primitive_overlap * backbone_coefficients)
    overlap_cross = Matrix{Float64}(backbone_coefficients' * primitive_overlap * raw_core_coefficients)
    residual_projector = vcat(
        -(overlap_backbone \ overlap_cross),
        Matrix{Float64}(I, ncore, ncore),
    )

    seed_coefficients = hcat(backbone_coefficients, raw_core_coefficients)
    residual_overlap = Matrix{Float64}(residual_projector' * (seed_coefficients' * primitive_overlap * seed_coefficients) * residual_projector)
    residual_decomposition = eigen(Symmetric(residual_overlap))
    keep = findall(>(1.0e-10), residual_decomposition.values)
    isempty(keep) && throw(ArgumentError("residual Gaussian interaction treatment produced no nontrivial residual Gaussian directions"))
    residual_transform =
        residual_decomposition.vectors[:, keep] *
        Diagonal(1.0 ./ sqrt.(residual_decomposition.values[keep]))
    residual_coefficients = seed_coefficients * residual_projector * residual_transform
    orthogonal_seed = hcat(backbone_coefficients, residual_coefficients)

    transfer = _normalized_density_transfer(orthogonal_seed' * primitive_overlap * hybrid_coefficients)
    raw_centers = Float64[gaussian.center_value for gaussian in hybrid.core_gaussians]
    residual_mixing = abs2.(Matrix{Float64}(residual_transform))
    residual_centers = Float64[]
    for column in axes(residual_mixing, 2)
        column_view = view(residual_mixing, :, column)
        total = sum(column_view)
        total > 0.0 || throw(ArgumentError("residual Gaussian center reconstruction requires nonzero Gaussian weight"))
        push!(residual_centers, dot(column_view, raw_centers) / total)
    end

    return (
        backbone = backbone,
        backbone_coefficients = backbone_coefficients,
        residual_coefficients = residual_coefficients,
        transfer = transfer,
        residual_centers = residual_centers,
    )
end

# Alg QW-RG step 8b: Compute residual-Gaussian first and second moments for MWG.
# See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
function _hybrid_residual_gaussian_mwg_data(hybrid::HybridMappedOrdinaryBasis1D)
    seed = _hybrid_residual_gaussian_seed(hybrid)
    primitive_layer = primitive_set(hybrid)
    primitive_overlap = Matrix{Float64}(overlap_matrix(primitive_layer))
    primitive_position = Matrix{Float64}(position_matrix(primitive_layer))
    primitive_x2 = Matrix{Float64}(_x2_matrix(primitive_layer))

    residual_overlap = Matrix{Float64}(seed.residual_coefficients' * primitive_overlap * seed.residual_coefficients)
    residual_position = Matrix{Float64}(seed.residual_coefficients' * primitive_position * seed.residual_coefficients)
    residual_x2 = Matrix{Float64}(seed.residual_coefficients' * primitive_x2 * seed.residual_coefficients)

    residual_centers = Float64[]
    residual_widths = Float64[]
    effective_gaussians = Gaussian[]
    for index in axes(residual_overlap, 1)
        norm_value = Float64(residual_overlap[index, index])
        norm_value > 1.0e-12 || throw(
            ArgumentError("MWG residual Gaussian reconstruction requires nonzero residual norm"),
        )
        center_value = Float64(residual_position[index, index] / norm_value)
        second_moment = Float64(residual_x2[index, index] / norm_value)
        variance = second_moment - center_value^2
        variance > 1.0e-12 || throw(
            ArgumentError("MWG residual Gaussian reconstruction requires positive residual variance"),
        )
        width_value = sqrt(2.0 * variance)
        push!(residual_centers, center_value)
        push!(residual_widths, width_value)
        push!(effective_gaussians, Gaussian(center = center_value, width = width_value))
    end

    return (
        backbone = seed.backbone,
        backbone_coefficients = seed.backbone_coefficients,
        residual_coefficients = seed.residual_coefficients,
        transfer = seed.transfer,
        residual_centers = residual_centers,
        residual_widths = residual_widths,
        effective_gaussians = effective_gaussians,
    )
end

function _hybrid_residual_gaussian_mwg_seed_layer(
    hybrid::HybridMappedOrdinaryBasis1D,
    mwg_data,
)
    hybrid_primitives = primitives(primitive_set(hybrid))
    effective_primitives = mwg_data.effective_gaussians
    primitive_layer = PrimitiveSet1D(
        AbstractPrimitiveFunction1D[vcat(hybrid_primitives, effective_primitives)...];
        name = :hybrid_residual_mwg_primitives,
    )

    nprimitive_hybrid = length(hybrid_primitives)
    nresidual = length(effective_primitives)
    nbackbone = size(mwg_data.backbone_coefficients, 2)
    coefficient_matrix = zeros(Float64, nprimitive_hybrid + nresidual, nbackbone + nresidual)
    coefficient_matrix[1:nprimitive_hybrid, 1:nbackbone] .= mwg_data.backbone_coefficients
    for index in 1:nresidual
        coefficient_matrix[nprimitive_hybrid + index, nbackbone + index] = 1.0
    end

    primitive_weights = Float64[integral_weight(primitive) for primitive in primitives(primitive_layer)]
    integral_weight_data = vec(transpose(coefficient_matrix) * primitive_weights)
    center_data = vcat(Float64.(centers(mwg_data.backbone)), mwg_data.residual_centers)
    return _AuxiliaryContractedLayer1D(
        primitive_layer,
        coefficient_matrix,
        center_data,
        integral_weight_data,
    )
end

# Alg QW-RG step 8a: Build nearest-center / GGT residual-Gaussian interaction data.
# See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
function _hybrid_residual_gaussian_pair_factors_nearest(
    hybrid::HybridMappedOrdinaryBasis1D,
    expansion::CoulombGaussianExpansion,
)
    seed = _hybrid_residual_gaussian_seed(hybrid)
    backbone = seed.backbone
    backbone_pair_factors = _pair_gaussian_factor_matrices(backbone; exponents = expansion.exponents)
    weight_backbone = Float64[Float64(weight) for weight in integral_weights(backbone)]
    any(weight -> abs(weight) <= 1.0e-12, weight_backbone) &&
        throw(ArgumentError("residual Gaussian interaction treatment requires nonzero backbone integral weights"))
    weight_outer = weight_backbone * transpose(weight_backbone)
    pair_backbone = [factor ./ weight_outer for factor in backbone_pair_factors]

    transfer_rg = _nearest_density_transfer(centers(backbone), seed.residual_centers)
    pair_seed = Matrix{Float64}[]
    for factor in pair_backbone
        factor_bg = factor * transfer_rg
        factor_gg = transpose(transfer_rg) * factor * transfer_rg
        push!(pair_seed, _symmetrize_ida_matrix([factor factor_bg; transpose(factor_bg) factor_gg]))
    end

    return [Matrix{Float64}(_symmetrize_ida_matrix(transpose(seed.transfer) * factor * seed.transfer)) for factor in pair_seed]
end

# Alg QW-RG step 8c and 8d: Build MWG effective Gaussian data and diagonal RG interaction factors.
# See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
function _hybrid_residual_gaussian_pair_factors_mwg(
    hybrid::HybridMappedOrdinaryBasis1D,
    expansion::CoulombGaussianExpansion,
)
    mwg_data = _hybrid_residual_gaussian_mwg_data(hybrid)
    seed_layer = _hybrid_residual_gaussian_mwg_seed_layer(hybrid, mwg_data)
    pair_seed_basis = _pair_gaussian_factor_matrices(seed_layer; exponents = expansion.exponents)
    weight_seed = Float64[Float64(weight) for weight in integral_weights(seed_layer)]
    any(weight -> abs(weight) <= 1.0e-12, weight_seed) &&
        throw(ArgumentError("MWG residual Gaussian interaction treatment requires nonzero seed integral weights"))
    weight_outer = weight_seed * transpose(weight_seed)
    pair_seed = [factor ./ weight_outer for factor in pair_seed_basis]
    return [Matrix{Float64}(_symmetrize_ida_matrix(transpose(mwg_data.transfer) * factor * mwg_data.transfer)) for factor in pair_seed]
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
one-dimensional basis on top of the quadrature-free local-linear analytic
Gaussian proxy.
`:pgdg_experimental` uses the same unlocalized analytic proxy.
`:numerical_reference` remains the validation route.
"""
function ordinary_cartesian_ida_operators(
    basis::MappedUniformBasis;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
    backend::Symbol = :pgdg_experimental,
)
    gausslet_bundle = _mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        center = 0.0,
        backend = backend,
    )
    return _ordinary_cartesian_ida_from_gausslet_bundle(
        gausslet_bundle,
        expansion,
        Z = Z,
    )
end

"""
    ordinary_cartesian_ida_operators(
        basis::HybridMappedOrdinaryBasis1D;
        expansion = coulomb_gaussian_expansion(doacc = false),
        Z = 2.0,
        interaction_treatment = :combined_basis,
    )

Build the current ordinary Cartesian static IDA ingredients for the hybrid
ordinary branch.

The hybrid path keeps the existing one-body construction and offers two narrow
interaction treatments:

- `:combined_basis` keeps the present combined hybrid-basis density-density
  construction
- `:residual_gaussian_nearest` uses a separate residual-Gaussian interaction
  ansatz that assigns the residual Gaussian channel to the nearest ordinary
  backbone centers before transferring the interaction back into the final
  localized hybrid basis
- `:residual_gaussian_mwg` uses the paper-style matched-width-Gaussian (MWG)
  approximation: residual Gaussians are orthogonalized against the backbone,
  their exact `x` and `x^2` moments are matched by effective Gaussian
  orbitals, and those effective orbitals are used to build the residual
  interaction seed before transferring it back into the final localized hybrid
  basis

The residual-Gaussian treatments are experimental and validation-oriented.
They are intended for the current `1s^2` scalar checks, not as a claim that
the hybrid ordinary branch is already solver-ready.
"""
function ordinary_cartesian_ida_operators(
    basis::HybridMappedOrdinaryBasis1D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
    interaction_treatment::Symbol = :combined_basis,
)
    one_body = mapped_ordinary_one_body_operators(
        basis;
        exponents = expansion.exponents,
        center = 0.0,
    )
    if interaction_treatment == :combined_basis
        return _ordinary_cartesian_ida_from_layer(
            basis,
            basis.backend,
            expansion,
            one_body,
            basis;
            Z = Z,
        )
    elseif interaction_treatment == :residual_gaussian_nearest
        pair_factors_1d = _hybrid_residual_gaussian_pair_factors_nearest(basis, expansion)
        _, one_body_hamiltonian = _mapped_cartesian_one_body_matrix(one_body, expansion; Z = Z)
        return _ordinary_cartesian_ida_from_pair_factors(
            basis,
            basis.backend,
            interaction_treatment,
            expansion,
            one_body,
            pair_factors_1d,
            one_body_hamiltonian,
            integral_weights(basis),
        )
    elseif interaction_treatment == :residual_gaussian_mwg
        pair_factors_1d = _hybrid_residual_gaussian_pair_factors_mwg(basis, expansion)
        _, one_body_hamiltonian = _mapped_cartesian_one_body_matrix(one_body, expansion; Z = Z)
        return _ordinary_cartesian_ida_from_pair_factors(
            basis,
            basis.backend,
            interaction_treatment,
            expansion,
            one_body,
            pair_factors_1d,
            one_body_hamiltonian,
            integral_weights(basis),
        )
    else
        throw(ArgumentError("hybrid ordinary Cartesian interaction_treatment must be :combined_basis, :residual_gaussian_nearest, or :residual_gaussian_mwg"))
    end
end
