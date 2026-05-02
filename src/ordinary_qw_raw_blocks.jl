function _qwrg_gausslet_1d_blocks(bundle::_MappedOrdinaryGausslet1DBundle)
    pgdg_intermediate = bundle.pgdg_intermediate
    return (
        overlap_gg = pgdg_intermediate.overlap,
        kinetic_gg = pgdg_intermediate.kinetic,
        position_gg = pgdg_intermediate.position,
        x2_gg = pgdg_intermediate.x2,
        factor_gg = pgdg_intermediate.gaussian_factors,
        factor_gg_terms = pgdg_intermediate.gaussian_factor_terms,
        pair_gg = pgdg_intermediate.pair_factors,
        pair_gg_terms = pgdg_intermediate.pair_factor_terms,
        weight_gg = pgdg_intermediate.weights,
        center_gg = pgdg_intermediate.centers,
    )
end

function _qwrg_basis_derivative_matrix(
    basis::MappedUniformBasis,
    points::AbstractVector{Float64},
)
    primitive_derivatives = _primitive_sample_matrix(
        primitive_set(basis),
        points;
        derivative_order = 1,
    )
    return primitive_derivatives * Matrix{Float64}(stencil_matrix(basis))
end

function _qwrg_gaussian_value_matrix(
    gaussians::AbstractVector{<:Gaussian},
    points::AbstractVector{Float64},
)
    return hcat([
        Float64[value(gaussian, point) for point in points] for gaussian in gaussians
    ]...)
end

function _qwrg_gaussian_derivative_matrix(
    gaussians::AbstractVector{<:Gaussian},
    points::AbstractVector{Float64},
)
    return hcat([
        Float64[derivative(gaussian, point; order = 1) for point in points] for gaussian in gaussians
    ]...)
end

function _qwrg_support_bounds(
    basis::MappedUniformBasis,
    gaussians::AbstractVector{<:Gaussian},
)
    basis_lo, basis_hi = _primitive_set_bounds(primitive_set(basis))
    if isempty(gaussians)
        return basis_lo, basis_hi
    end
    gaussian_set = PrimitiveSet1D(
        AbstractPrimitiveFunction1D[gaussian for gaussian in gaussians];
        name = :qiu_white_support_gaussians,
    )
    gaussian_lo, gaussian_hi = _primitive_set_bounds(gaussian_set)
    return min(basis_lo, gaussian_lo), max(basis_hi, gaussian_hi)
end

function _qwrg_cross_1d_blocks_midpoint(
    basis::MappedUniformBasis,
    gaussians::AbstractVector{<:Gaussian},
    expansion::CoulombGaussianExpansion;
    h = nothing,
    include_pair::Bool = true,
)
    xlo, xhi = _qwrg_support_bounds(basis, gaussians)
    h_try = h === nothing ? _primitive_matrix_start_h(primitive_set(basis)) : Float64(h)
    h_try > 0.0 || throw(ArgumentError("Qiu-White split cross-block construction requires h > 0"))

    previous = nothing
    current = nothing
    exponent_values = Float64[Float64(exponent) for exponent in expansion.exponents]
    for _ in 1:_PRIMITIVE_MATRIX_MAXITER
        points, weights = _make_midpoint_grid(xlo, xhi, h_try)
        basis_values = _basis_sample_matrix(basis, points)
        basis_derivatives = _qwrg_basis_derivative_matrix(basis, points)
        gaussian_values = _qwrg_gaussian_value_matrix(gaussians, points)
        gaussian_derivatives = _qwrg_gaussian_derivative_matrix(gaussians, points)
        weighted_basis = weights .* basis_values
        weighted_gaussians = weights .* gaussian_values
        overlap = Matrix{Float64}(transpose(basis_values) * weighted_gaussians)
        kinetic = Matrix{Float64}(0.5 .* transpose(basis_derivatives) * (weights .* gaussian_derivatives))
        position = Matrix{Float64}(transpose(basis_values) * ((weights .* points) .* gaussian_values))
        x2 = Matrix{Float64}(transpose(basis_values) * ((weights .* (points .^ 2)) .* gaussian_values))
        factors = Matrix{Float64}[
            Matrix{Float64}(transpose(basis_values) * ((weights .* exp.(-exponent .* (points .^ 2))) .* gaussian_values)) for exponent in exponent_values
        ]
        pair_factors = Matrix{Float64}[]
        if include_pair
            distance2 = (points .- transpose(points)) .^ 2
            for exponent in exponent_values
                kernel = exp.(-exponent .* distance2)
                push!(pair_factors, Matrix{Float64}(transpose(weighted_basis) * (kernel * weighted_gaussians)))
            end
        end

        current = (
            overlap_ga = overlap,
            kinetic_ga = kinetic,
            position_ga = position,
            x2_ga = x2,
            factor_ga = factors,
            pair_ga = pair_factors,
        )
        if previous !== nothing
            differences = Float64[
                norm(current.overlap_ga - previous.overlap_ga, Inf),
                norm(current.kinetic_ga - previous.kinetic_ga, Inf),
                norm(current.position_ga - previous.position_ga, Inf),
                norm(current.x2_ga - previous.x2_ga, Inf),
            ]
            append!(
                differences,
                [norm(current.factor_ga[index] - previous.factor_ga[index], Inf) for index in eachindex(current.factor_ga)],
            )
            if include_pair
                append!(
                    differences,
                    [norm(current.pair_ga[index] - previous.pair_ga[index], Inf) for index in eachindex(current.pair_ga)],
                )
            end
            maxdiff = maximum(differences)
            maxdiff <= _PRIMITIVE_MATRIX_TOL && return current
        end
        previous = current
        h_try /= 2.0
    end
    return current
end

function _qwrg_proxy_gaussian_primitives(
    proxy_layer::_MappedLegacyProxyLayer1D,
)
    return Gaussian[primitive for primitive in primitives(primitive_set(proxy_layer))]
end

function _qwrg_supplement_primitives_and_contraction(
    gaussians::AbstractVector{<:Gaussian},
)
    primitive_gaussians = Gaussian[gaussian for gaussian in gaussians]
    contraction_matrix = Matrix{Float64}(I, length(primitive_gaussians), length(primitive_gaussians))
    return primitive_gaussians, contraction_matrix
end

function _qwrg_supplement_primitives_and_contraction(
    gaussian_data::LegacyAtomicGaussianSupplement,
)
    _legacy_atomic_has_nonseparable_shells(gaussian_data) && throw(
        ArgumentError(
            "ordinary_cartesian_qiu_white_operators does not yet support true active l > 0 atomic supplements; the current QW path still assumes one centered separable 1D supplement channel and needs an explicit 3D Cartesian shell supplement route before lmax = 1 can be used honestly",
        ),
    )
    return (
        Gaussian[gaussian for gaussian in gaussian_data.primitive_gaussians],
        Matrix{Float64}(gaussian_data.contraction_matrix),
    )
end

_qwrg_gaussian_exponent(gaussian::Gaussian) =
    GaussianAnalyticIntegrals.gaussian_exponent(gaussian)

function _qwrg_atomic_shell_prefactor(exponent::Float64, power::Int)
    power >= 0 || throw(ArgumentError("atomic shell prefactor requires power >= 0"))
    return GaussianAnalyticIntegrals.polynomial_gaussian_shell_prefactor(exponent, power)
end

_qwrg_doublefactorial(n::Int) =
    GaussianAnalyticIntegrals.doublefactorial(n)

_qwrg_shifted_gaussian_moment(gamma::Float64, power::Int) =
    GaussianAnalyticIntegrals.shifted_gaussian_moment(gamma, power)

_qwrg_poly_shift_multiply(coefficients::Vector{Float64}, shift::Float64, power::Int) =
    GaussianAnalyticIntegrals.polynomial_shift_multiply(coefficients, shift, power)

function _qwrg_atomic_basic_integral(
    alpha_left::Float64,
    center_left::Float64,
    power_left::Int,
    prefactor_left::Float64,
    alpha_right::Float64,
    center_right::Float64,
    power_right::Int,
    prefactor_right::Float64;
    xpower::Int = 0,
    extra_exponent::Float64 = 0.0,
    extra_center::Float64 = 0.0,
)
    gamma = alpha_left + alpha_right + extra_exponent
    gamma > 0.0 ||
        throw(ArgumentError("atomic shell integral requires positive total Gaussian exponent"))
    return GaussianAnalyticIntegrals.polynomial_gaussian_basic_integral(
        alpha_left,
        center_left,
        power_left,
        prefactor_left,
        alpha_right,
        center_right,
        power_right,
        prefactor_right;
        xpower,
        extra_exponent,
        extra_center,
    )
end

function _qwrg_atomic_derivative_terms(power::Int, exponent::Float64)
    if !(power in (0, 1, 2))
        throw(
            ArgumentError(
                "atomic shell derivative terms currently support only powers 0, 1, and 2",
            ),
        )
    end
    return GaussianAnalyticIntegrals.polynomial_gaussian_derivative_terms(power, exponent)
end

function _qwrg_atomic_kinetic_integral(
    alpha_left::Float64,
    center_left::Float64,
    power_left::Int,
    prefactor_left::Float64,
    alpha_right::Float64,
    center_right::Float64,
    power_right::Int,
    prefactor_right::Float64,
)
    if !(power_left in (0, 1, 2)) || !(power_right in (0, 1, 2))
        throw(
            ArgumentError(
                "atomic shell derivative terms currently support only powers 0, 1, and 2",
            ),
        )
    end
    return GaussianAnalyticIntegrals.polynomial_gaussian_kinetic_integral(
        alpha_left,
        center_left,
        power_left,
        prefactor_left,
        alpha_right,
        center_right,
        power_right,
        prefactor_right,
    )
end

function _qwrg_atomic_orbital_axis_power(
    orbital::_AtomicCartesianShellOrbital3D,
    axis::Symbol,
)
    axis == :x && return orbital.lx
    axis == :y && return orbital.ly
    axis == :z && return orbital.lz
    throw(ArgumentError("axis must be :x, :y, or :z"))
end

function _qwrg_atomic_axis_cross_data(
    proxy_layer::_MappedLegacyProxyLayer1D,
    orbital::_AtomicCartesianShellOrbital3D,
    axis::Symbol,
    expansion::CoulombGaussianExpansion,
)
    proxy_gaussians = _qwrg_proxy_gaussian_primitives(proxy_layer)
    center_value = axis == :x ? orbital.center[1] : axis == :y ? orbital.center[2] : orbital.center[3]
    power = _qwrg_atomic_orbital_axis_power(orbital, axis)
    nproxy = length(proxy_gaussians)
    nprimitive = length(orbital.exponents)

    overlap = zeros(Float64, nproxy, nprimitive)
    kinetic = zeros(Float64, nproxy, nprimitive)
    position = zeros(Float64, nproxy, nprimitive)
    x2 = zeros(Float64, nproxy, nprimitive)
    factor = [zeros(Float64, nproxy, nprimitive) for _ in expansion.exponents]

    for primitive in 1:nprimitive
        exponent = Float64(orbital.exponents[primitive])
        prefactor = _qwrg_atomic_shell_prefactor(exponent, power)
        for proxy in 1:nproxy
            gaussian = proxy_gaussians[proxy]
            alpha_proxy = _qwrg_gaussian_exponent(gaussian)
            overlap[proxy, primitive] = _qwrg_atomic_basic_integral(
                alpha_proxy,
                gaussian.center_value,
                0,
                1.0,
                exponent,
                center_value,
                power,
                prefactor,
            )
            kinetic[proxy, primitive] = _qwrg_atomic_kinetic_integral(
                alpha_proxy,
                gaussian.center_value,
                0,
                1.0,
                exponent,
                center_value,
                power,
                prefactor,
            )
            position[proxy, primitive] = _qwrg_atomic_basic_integral(
                alpha_proxy,
                gaussian.center_value,
                0,
                1.0,
                exponent,
                center_value,
                power,
                prefactor;
                xpower = 1,
            )
            x2[proxy, primitive] = _qwrg_atomic_basic_integral(
                alpha_proxy,
                gaussian.center_value,
                0,
                1.0,
                exponent,
                center_value,
                power,
                prefactor;
                xpower = 2,
            )
            for term in eachindex(expansion.exponents)
                factor[term][proxy, primitive] = _qwrg_atomic_basic_integral(
                    alpha_proxy,
                    gaussian.center_value,
                    0,
                    1.0,
                    exponent,
                    center_value,
                    power,
                    prefactor;
                    extra_exponent = Float64(expansion.exponents[term]),
                    extra_center = 0.0,
                )
            end
        end
    end

    return (
        overlap = _qwrg_left_contract_cross_matrix(proxy_layer, overlap),
        kinetic = _qwrg_left_contract_cross_matrix(proxy_layer, kinetic),
        position = _qwrg_left_contract_cross_matrix(proxy_layer, position),
        x2 = _qwrg_left_contract_cross_matrix(proxy_layer, x2),
        factor = [_qwrg_left_contract_cross_matrix(proxy_layer, matrix) for matrix in factor],
    )
end

function _qwrg_atomic_axis_aa_data(
    left::_AtomicCartesianShellOrbital3D,
    right::_AtomicCartesianShellOrbital3D,
    axis::Symbol,
    expansion::CoulombGaussianExpansion,
)
    center_left = axis == :x ? left.center[1] : axis == :y ? left.center[2] : left.center[3]
    center_right = axis == :x ? right.center[1] : axis == :y ? right.center[2] : right.center[3]
    power_left = _qwrg_atomic_orbital_axis_power(left, axis)
    power_right = _qwrg_atomic_orbital_axis_power(right, axis)
    nleft = length(left.exponents)
    nright = length(right.exponents)

    overlap = zeros(Float64, nleft, nright)
    kinetic = zeros(Float64, nleft, nright)
    position = zeros(Float64, nleft, nright)
    x2 = zeros(Float64, nleft, nright)
    factor = [zeros(Float64, nleft, nright) for _ in expansion.exponents]

    for j in 1:nright
        exponent_right = Float64(right.exponents[j])
        prefactor_right = _qwrg_atomic_shell_prefactor(exponent_right, power_right)
        for i in 1:nleft
            exponent_left = Float64(left.exponents[i])
            prefactor_left = _qwrg_atomic_shell_prefactor(exponent_left, power_left)
            overlap[i, j] = _qwrg_atomic_basic_integral(
                exponent_left,
                center_left,
                power_left,
                prefactor_left,
                exponent_right,
                center_right,
                power_right,
                prefactor_right,
            )
            kinetic[i, j] = _qwrg_atomic_kinetic_integral(
                exponent_left,
                center_left,
                power_left,
                prefactor_left,
                exponent_right,
                center_right,
                power_right,
                prefactor_right,
            )
            position[i, j] = _qwrg_atomic_basic_integral(
                exponent_left,
                center_left,
                power_left,
                prefactor_left,
                exponent_right,
                center_right,
                power_right,
                prefactor_right;
                xpower = 1,
            )
            x2[i, j] = _qwrg_atomic_basic_integral(
                exponent_left,
                center_left,
                power_left,
                prefactor_left,
                exponent_right,
                center_right,
                power_right,
                prefactor_right;
                xpower = 2,
            )
            for term in eachindex(expansion.exponents)
                factor[term][i, j] = _qwrg_atomic_basic_integral(
                    exponent_left,
                    center_left,
                    power_left,
                    prefactor_left,
                    exponent_right,
                    center_right,
                    power_right,
                    prefactor_right;
                    extra_exponent = Float64(expansion.exponents[term]),
                    extra_center = 0.0,
                )
            end
        end
    end

    return (
        overlap = overlap,
        kinetic = kinetic,
        position = position,
        x2 = x2,
        factor = factor,
    )
end

function _qwrg_atomic_axis_aa_one_body_data(
    left::_AtomicCartesianShellOrbital3D,
    right::_AtomicCartesianShellOrbital3D,
    axis::Symbol,
    expansion::CoulombGaussianExpansion,
)
    center_left = axis == :x ? left.center[1] : axis == :y ? left.center[2] : left.center[3]
    center_right = axis == :x ? right.center[1] : axis == :y ? right.center[2] : right.center[3]
    power_left = _qwrg_atomic_orbital_axis_power(left, axis)
    power_right = _qwrg_atomic_orbital_axis_power(right, axis)
    nleft = length(left.exponents)
    nright = length(right.exponents)

    overlap = zeros(Float64, nleft, nright)
    kinetic = zeros(Float64, nleft, nright)
    position = zeros(Float64, nleft, nright)
    x2 = zeros(Float64, nleft, nright)

    for j in 1:nright
        exponent_right = Float64(right.exponents[j])
        prefactor_right = _qwrg_atomic_shell_prefactor(exponent_right, power_right)
        for i in 1:nleft
            exponent_left = Float64(left.exponents[i])
            prefactor_left = _qwrg_atomic_shell_prefactor(exponent_left, power_left)
            overlap[i, j] = _qwrg_atomic_basic_integral(
                exponent_left,
                center_left,
                power_left,
                prefactor_left,
                exponent_right,
                center_right,
                power_right,
                prefactor_right,
            )
            kinetic[i, j] = _qwrg_atomic_kinetic_integral(
                exponent_left,
                center_left,
                power_left,
                prefactor_left,
                exponent_right,
                center_right,
                power_right,
                prefactor_right,
            )
            position[i, j] = _qwrg_atomic_basic_integral(
                exponent_left,
                center_left,
                power_left,
                prefactor_left,
                exponent_right,
                center_right,
                power_right,
                prefactor_right;
                xpower = 1,
            )
            x2[i, j] = _qwrg_atomic_basic_integral(
                exponent_left,
                center_left,
                power_left,
                prefactor_left,
                exponent_right,
                center_right,
                power_right,
                prefactor_right;
                xpower = 2,
            )
        end
    end

    return (
        overlap = overlap,
        kinetic = kinetic,
        position = position,
        x2 = x2,
    )
end

function _qwrg_atomic_axis_factor_cross_data(
    proxy_layer::_MappedLegacyProxyLayer1D,
    orbital::_AtomicCartesianShellOrbital3D,
    axis::Symbol,
    expansion::CoulombGaussianExpansion,
    factor_center::Float64,
)
    proxy_gaussians = _qwrg_proxy_gaussian_primitives(proxy_layer)
    center_value = axis == :x ? orbital.center[1] : axis == :y ? orbital.center[2] : orbital.center[3]
    power = _qwrg_atomic_orbital_axis_power(orbital, axis)
    nproxy = length(proxy_gaussians)
    nprimitive = length(orbital.exponents)
    factor = [zeros(Float64, nproxy, nprimitive) for _ in expansion.exponents]

    for primitive in 1:nprimitive
        exponent = Float64(orbital.exponents[primitive])
        prefactor = _qwrg_atomic_shell_prefactor(exponent, power)
        for proxy in 1:nproxy
            gaussian = proxy_gaussians[proxy]
            alpha_proxy = _qwrg_gaussian_exponent(gaussian)
            for term in eachindex(expansion.exponents)
                factor[term][proxy, primitive] = _qwrg_atomic_basic_integral(
                    alpha_proxy,
                    gaussian.center_value,
                    0,
                    1.0,
                    exponent,
                    center_value,
                    power,
                    prefactor;
                    extra_exponent = Float64(expansion.exponents[term]),
                    extra_center = factor_center,
                )
            end
        end
    end

    return [_qwrg_left_contract_cross_matrix(proxy_layer, matrix) for matrix in factor]
end

function _qwrg_atomic_axis_factor_aa_data(
    left::_AtomicCartesianShellOrbital3D,
    right::_AtomicCartesianShellOrbital3D,
    axis::Symbol,
    expansion::CoulombGaussianExpansion,
    factor_center::Float64,
)
    center_left = axis == :x ? left.center[1] : axis == :y ? left.center[2] : left.center[3]
    center_right = axis == :x ? right.center[1] : axis == :y ? right.center[2] : right.center[3]
    power_left = _qwrg_atomic_orbital_axis_power(left, axis)
    power_right = _qwrg_atomic_orbital_axis_power(right, axis)
    nleft = length(left.exponents)
    nright = length(right.exponents)
    factor = [zeros(Float64, nleft, nright) for _ in expansion.exponents]

    for j in 1:nright
        exponent_right = Float64(right.exponents[j])
        prefactor_right = _qwrg_atomic_shell_prefactor(exponent_right, power_right)
        for i in 1:nleft
            exponent_left = Float64(left.exponents[i])
            prefactor_left = _qwrg_atomic_shell_prefactor(exponent_left, power_left)
            for term in eachindex(expansion.exponents)
                factor[term][i, j] = _qwrg_atomic_basic_integral(
                    exponent_left,
                    center_left,
                    power_left,
                    prefactor_left,
                    exponent_right,
                    center_right,
                    power_right,
                    prefactor_right;
                    extra_exponent = Float64(expansion.exponents[term]),
                    extra_center = factor_center,
                )
            end
        end
    end

    return factor
end

function _qwrg_atomic_weighted_hadamard(
    left_coefficients::AbstractVector{<:Real},
    x::AbstractMatrix{<:Real},
    y::AbstractMatrix{<:Real},
    z::AbstractMatrix{<:Real},
    right_coefficients::AbstractVector{<:Real},
)
    matrix = Matrix{Float64}(x) .* Matrix{Float64}(y) .* Matrix{Float64}(z)
    return Float64(dot(left_coefficients, matrix * right_coefficients))
end

function _qwrg_gausslet_axis_matrix(data, axis::Symbol; squared::Bool = false)
    gg_axis = squared ? data.x2_gg : data.position_gg
    overlap_gg = data.overlap_gg
    ngausslet3d = length(_mapped_cartesian_orbitals(1:size(overlap_gg, 1)))
    gg_block = zeros(Float64, ngausslet3d, ngausslet3d)
    if axis == :x
        _qwrg_fill_product_matrix!(gg_block, gg_axis, overlap_gg, overlap_gg)
    elseif axis == :y
        _qwrg_fill_product_matrix!(gg_block, overlap_gg, gg_axis, overlap_gg)
    elseif axis == :z
        _qwrg_fill_product_matrix!(gg_block, overlap_gg, overlap_gg, gg_axis)
    else
        throw(ArgumentError("axis must be :x, :y, or :z"))
    end
    return gg_block
end

function _qwrg_atomic_cartesian_one_body_aa(
    blocks,
    expansion::CoulombGaussianExpansion;
    Z::Real,
)
    matrix = Matrix{Float64}(blocks.kinetic_aa)
    z_value = Float64(Z)
    for term in eachindex(expansion.coefficients)
        matrix .-= z_value * expansion.coefficients[term] .* blocks.factor_aa[term]
    end
    return Matrix{Float64}(0.5 .* (matrix .+ transpose(matrix)))
end

function _qwrg_atomic_cartesian_blocks_3d(
    gausslet_bundle::_MappedOrdinaryGausslet1DBundle,
    supplement::_AtomicCartesianShellSupplement3D,
    expansion::CoulombGaussianExpansion,
)
    proxy_layer = gausslet_bundle.pgdg_intermediate.auxiliary_layer
    proxy_layer isa _MappedLegacyProxyLayer1D || throw(
        ArgumentError("explicit atomic Cartesian shell supplement currently requires the base refinement_levels = 0 PGDG proxy line on the gausslet side"),
    )

    n1d = size(gausslet_bundle.pgdg_intermediate.overlap, 1)
    ngausslet3d = n1d^3
    norbital = length(supplement.orbitals)
    overlap_ga = zeros(Float64, ngausslet3d, norbital)
    kinetic_ga = zeros(Float64, ngausslet3d, norbital)
    position_x_ga = zeros(Float64, ngausslet3d, norbital)
    position_y_ga = zeros(Float64, ngausslet3d, norbital)
    position_z_ga = zeros(Float64, ngausslet3d, norbital)
    x2_x_ga = zeros(Float64, ngausslet3d, norbital)
    x2_y_ga = zeros(Float64, ngausslet3d, norbital)
    x2_z_ga = zeros(Float64, ngausslet3d, norbital)
    factor_ga = [zeros(Float64, ngausslet3d, norbital) for _ in expansion.exponents]

    overlap_aa = zeros(Float64, norbital, norbital)
    kinetic_aa = zeros(Float64, norbital, norbital)
    position_x_aa = zeros(Float64, norbital, norbital)
    position_y_aa = zeros(Float64, norbital, norbital)
    position_z_aa = zeros(Float64, norbital, norbital)
    x2_x_aa = zeros(Float64, norbital, norbital)
    x2_y_aa = zeros(Float64, norbital, norbital)
    x2_z_aa = zeros(Float64, norbital, norbital)
    factor_aa = [zeros(Float64, norbital, norbital) for _ in expansion.exponents]

    cross_cache = Dict{Tuple{Int,Symbol},Any}()
    aa_cache = Dict{Tuple{Int,Int,Symbol},Any}()

    for (orbital_index, orbital) in pairs(supplement.orbitals)
        x_data = get!(cross_cache, (orbital_index, :x)) do
            _qwrg_atomic_axis_cross_data(proxy_layer, orbital, :x, expansion)
        end
        y_data = get!(cross_cache, (orbital_index, :y)) do
            _qwrg_atomic_axis_cross_data(proxy_layer, orbital, :y, expansion)
        end
        z_data = get!(cross_cache, (orbital_index, :z)) do
            _qwrg_atomic_axis_cross_data(proxy_layer, orbital, :z, expansion)
        end

        scratch = zeros(Float64, ngausslet3d)
        for primitive in eachindex(orbital.coefficients)
            coefficient = Float64(orbital.coefficients[primitive])
            _qwrg_fill_product_column!(
                scratch,
                view(x_data.overlap, :, primitive),
                view(y_data.overlap, :, primitive),
                view(z_data.overlap, :, primitive),
            )
            overlap_ga[:, orbital_index] .+= coefficient .* scratch

            _qwrg_fill_product_column!(
                scratch,
                view(x_data.kinetic, :, primitive),
                view(y_data.overlap, :, primitive),
                view(z_data.overlap, :, primitive),
            )
            kinetic_ga[:, orbital_index] .+= coefficient .* scratch
            _qwrg_fill_product_column!(
                scratch,
                view(x_data.overlap, :, primitive),
                view(y_data.kinetic, :, primitive),
                view(z_data.overlap, :, primitive),
            )
            kinetic_ga[:, orbital_index] .+= coefficient .* scratch
            _qwrg_fill_product_column!(
                scratch,
                view(x_data.overlap, :, primitive),
                view(y_data.overlap, :, primitive),
                view(z_data.kinetic, :, primitive),
            )
            kinetic_ga[:, orbital_index] .+= coefficient .* scratch

            _qwrg_fill_product_column!(
                scratch,
                view(x_data.position, :, primitive),
                view(y_data.overlap, :, primitive),
                view(z_data.overlap, :, primitive),
            )
            position_x_ga[:, orbital_index] .+= coefficient .* scratch
            _qwrg_fill_product_column!(
                scratch,
                view(x_data.overlap, :, primitive),
                view(y_data.position, :, primitive),
                view(z_data.overlap, :, primitive),
            )
            position_y_ga[:, orbital_index] .+= coefficient .* scratch
            _qwrg_fill_product_column!(
                scratch,
                view(x_data.overlap, :, primitive),
                view(y_data.overlap, :, primitive),
                view(z_data.position, :, primitive),
            )
            position_z_ga[:, orbital_index] .+= coefficient .* scratch

            _qwrg_fill_product_column!(
                scratch,
                view(x_data.x2, :, primitive),
                view(y_data.overlap, :, primitive),
                view(z_data.overlap, :, primitive),
            )
            x2_x_ga[:, orbital_index] .+= coefficient .* scratch
            _qwrg_fill_product_column!(
                scratch,
                view(x_data.overlap, :, primitive),
                view(y_data.x2, :, primitive),
                view(z_data.overlap, :, primitive),
            )
            x2_y_ga[:, orbital_index] .+= coefficient .* scratch
            _qwrg_fill_product_column!(
                scratch,
                view(x_data.overlap, :, primitive),
                view(y_data.overlap, :, primitive),
                view(z_data.x2, :, primitive),
            )
            x2_z_ga[:, orbital_index] .+= coefficient .* scratch

            for term in eachindex(expansion.exponents)
                _qwrg_fill_product_column!(
                    scratch,
                    view(x_data.factor[term], :, primitive),
                    view(y_data.factor[term], :, primitive),
                    view(z_data.factor[term], :, primitive),
                )
                factor_ga[term][:, orbital_index] .+= coefficient .* scratch
            end
        end
    end

    for (left_index, left) in pairs(supplement.orbitals), (right_index, right) in pairs(supplement.orbitals)
        x_data = get!(aa_cache, (left_index, right_index, :x)) do
            _qwrg_atomic_axis_aa_data(left, right, :x, expansion)
        end
        y_data = get!(aa_cache, (left_index, right_index, :y)) do
            _qwrg_atomic_axis_aa_data(left, right, :y, expansion)
        end
        z_data = get!(aa_cache, (left_index, right_index, :z)) do
            _qwrg_atomic_axis_aa_data(left, right, :z, expansion)
        end
        overlap_aa[left_index, right_index] = _qwrg_atomic_weighted_hadamard(
            left.coefficients,
            x_data.overlap,
            y_data.overlap,
            z_data.overlap,
            right.coefficients,
        )
        kinetic_aa[left_index, right_index] =
            _qwrg_atomic_weighted_hadamard(
                left.coefficients,
                x_data.kinetic,
                y_data.overlap,
                z_data.overlap,
                right.coefficients,
            ) +
            _qwrg_atomic_weighted_hadamard(
                left.coefficients,
                x_data.overlap,
                y_data.kinetic,
                z_data.overlap,
                right.coefficients,
            ) +
            _qwrg_atomic_weighted_hadamard(
                left.coefficients,
                x_data.overlap,
                y_data.overlap,
                z_data.kinetic,
                right.coefficients,
            )
        position_x_aa[left_index, right_index] = _qwrg_atomic_weighted_hadamard(
            left.coefficients,
            x_data.position,
            y_data.overlap,
            z_data.overlap,
            right.coefficients,
        )
        position_y_aa[left_index, right_index] = _qwrg_atomic_weighted_hadamard(
            left.coefficients,
            x_data.overlap,
            y_data.position,
            z_data.overlap,
            right.coefficients,
        )
        position_z_aa[left_index, right_index] = _qwrg_atomic_weighted_hadamard(
            left.coefficients,
            x_data.overlap,
            y_data.overlap,
            z_data.position,
            right.coefficients,
        )
        x2_x_aa[left_index, right_index] = _qwrg_atomic_weighted_hadamard(
            left.coefficients,
            x_data.x2,
            y_data.overlap,
            z_data.overlap,
            right.coefficients,
        )
        x2_y_aa[left_index, right_index] = _qwrg_atomic_weighted_hadamard(
            left.coefficients,
            x_data.overlap,
            y_data.x2,
            z_data.overlap,
            right.coefficients,
        )
        x2_z_aa[left_index, right_index] = _qwrg_atomic_weighted_hadamard(
            left.coefficients,
            x_data.overlap,
            y_data.overlap,
            z_data.x2,
            right.coefficients,
        )
        for term in eachindex(expansion.exponents)
            factor_aa[term][left_index, right_index] = _qwrg_atomic_weighted_hadamard(
                left.coefficients,
                x_data.factor[term],
                y_data.factor[term],
                z_data.factor[term],
                right.coefficients,
            )
        end
    end

    return (
        overlap_ga = overlap_ga,
        overlap_aa = overlap_aa,
        kinetic_ga = kinetic_ga,
        kinetic_aa = kinetic_aa,
        position_x_ga = position_x_ga,
        position_y_ga = position_y_ga,
        position_z_ga = position_z_ga,
        position_x_aa = position_x_aa,
        position_y_aa = position_y_aa,
        position_z_aa = position_z_aa,
        x2_x_ga = x2_x_ga,
        x2_y_ga = x2_y_ga,
        x2_z_ga = x2_z_ga,
        x2_x_aa = x2_x_aa,
        x2_y_aa = x2_y_aa,
        x2_z_aa = x2_z_aa,
        factor_ga = factor_ga,
        factor_aa = factor_aa,
    )
end

function _qwrg_left_contract_cross_matrix(
    proxy_layer::_MappedLegacyProxyLayer1D,
    primitive_cross::AbstractMatrix{<:Real},
)
    coefficient_matrix = Matrix{Float64}(stencil_matrix(proxy_layer))
    size(primitive_cross, 1) == size(coefficient_matrix, 1) || throw(
        DimensionMismatch("Qiu-White proxy cross contraction requires primitive row count to match proxy coefficients"),
    )
    return Matrix{Float64}(transpose(coefficient_matrix) * primitive_cross)
end

function _qwrg_diatomic_supplement_proxy_layer(
    axis_basis::MappedUniformBasis,
    bundle::_MappedOrdinaryGausslet1DBundle,
    axis::Symbol,
)
    if bundle.backend == :pgdg_localized_experimental
        return _mapped_legacy_proxy_localized(_mapped_legacy_proxy_layer(axis_basis)).layer
    end

    proxy_layer = bundle.pgdg_intermediate.auxiliary_layer
    proxy_layer isa _MappedLegacyProxyLayer1D || throw(
        ArgumentError("bond-aligned diatomic molecular supplement currently requires the base refinement_levels = 0 PGDG proxy line on the $(axis) axis"),
    )
    return proxy_layer
end

function _qwrg_right_contract_cross_matrix(
    primitive_cross::AbstractMatrix{<:Real},
    contraction_matrix::AbstractMatrix{<:Real},
)
    size(primitive_cross, 2) == size(contraction_matrix, 1) || throw(
        DimensionMismatch("Qiu-White proxy cross contraction requires primitive column count to match supplement contraction rows"),
    )
    return Matrix{Float64}(Matrix{Float64}(primitive_cross) * Matrix{Float64}(contraction_matrix))
end

function _qwrg_gaussian_cross_matrix(
    left::AbstractVector{<:Gaussian},
    right::AbstractVector{<:Gaussian},
    builder,
)
    matrix = zeros(Float64, length(left), length(right))
    @inbounds for j in eachindex(right), i in eachindex(left)
        matrix[i, j] = builder(left[i], right[j])
    end
    return matrix
end

function _qwrg_gaussian_cross_matrices(
    left::AbstractVector{<:Gaussian},
    right::AbstractVector{<:Gaussian},
    builders::NamedTuple,
)
    return (
        overlap = _qwrg_gaussian_cross_matrix(left, right, builders.overlap),
        kinetic = _qwrg_gaussian_cross_matrix(left, right, builders.kinetic),
        position = _qwrg_gaussian_cross_matrix(left, right, builders.position),
        x2 = _qwrg_gaussian_cross_matrix(left, right, builders.x2),
    )
end

function _qwrg_gaussian_cross_matrices(
    left::AbstractVector{<:Gaussian},
    right::AbstractVector{<:Gaussian},
    exponents::AbstractVector{<:Real},
    builder,
)
    matrices = Matrix{Float64}[]
    for exponent in exponents
        push!(matrices, _qwrg_gaussian_cross_matrix(left, right, (a, b) -> builder(a, b, Float64(exponent))))
    end
    return matrices
end

function _qwrg_cross_1d_blocks_proxy(
    proxy_layer::_MappedLegacyProxyLayer1D,
    supplement,
    expansion::CoulombGaussianExpansion,
)
    gaussians, contraction_matrix = _qwrg_supplement_primitives_and_contraction(supplement)
    proxy_gaussians = _qwrg_proxy_gaussian_primitives(proxy_layer)
    scalar_raw = _qwrg_gaussian_cross_matrices(
        proxy_gaussians,
        gaussians,
        (
            overlap = _gaussian_overlap,
            kinetic = _gaussian_kinetic,
            position = _gaussian_position,
            x2 = _gaussian_x2,
        ),
    )
    factor_raw = _qwrg_gaussian_cross_matrices(
        proxy_gaussians,
        gaussians,
        expansion.exponents,
        (a, b, exponent) -> _gaussian_factor(a, b, exponent, 0.0),
    )
    pair_raw = _qwrg_gaussian_cross_matrices(
        proxy_gaussians,
        gaussians,
        expansion.exponents,
        _gaussian_pair_factor,
    )

    return (
        overlap_ga = _qwrg_right_contract_cross_matrix(
            _qwrg_left_contract_cross_matrix(proxy_layer, scalar_raw.overlap),
            contraction_matrix,
        ),
        kinetic_ga = _qwrg_right_contract_cross_matrix(
            _qwrg_left_contract_cross_matrix(proxy_layer, scalar_raw.kinetic),
            contraction_matrix,
        ),
        position_ga = _qwrg_right_contract_cross_matrix(
            _qwrg_left_contract_cross_matrix(proxy_layer, scalar_raw.position),
            contraction_matrix,
        ),
        x2_ga = _qwrg_right_contract_cross_matrix(
            _qwrg_left_contract_cross_matrix(proxy_layer, scalar_raw.x2),
            contraction_matrix,
        ),
        factor_ga = [
            _qwrg_right_contract_cross_matrix(
                _qwrg_left_contract_cross_matrix(proxy_layer, matrix),
                contraction_matrix,
            ) for matrix in factor_raw
        ],
        pair_ga = [
            _qwrg_right_contract_cross_matrix(
                _qwrg_left_contract_cross_matrix(proxy_layer, matrix),
                contraction_matrix,
            ) for matrix in pair_raw
        ],
    )
end

function _qwrg_cross_1d_blocks(
    gausslet_bundle::_MappedOrdinaryGausslet1DBundle,
    supplement,
    expansion::CoulombGaussianExpansion,
)
    # Alg QW-RG step 6: Build the gausslet-Gaussian raw 1D cross blocks on the
    # same contracted proxy route used successfully for the pair terms:
    # analytic Gaussian primitive cross blocks plus contraction only on the
    # gausslet side. This keeps the base PGDG one-body side off midpoint
    # sampling as well.
    # See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
    proxy_layer = gausslet_bundle.pgdg_intermediate.auxiliary_layer
    proxy_layer isa _MappedLegacyProxyLayer1D || throw(
        ArgumentError("Qiu-White cross-block construction currently requires the base refinement_levels = 0 PGDG proxy line on the gausslet side"),
    )
    return _qwrg_cross_1d_blocks_proxy(proxy_layer, supplement, expansion)
end

function _qwrg_split_block_matrices(
    gausslet_bundle::_MappedOrdinaryGausslet1DBundle,
    supplement,
    expansion::CoulombGaussianExpansion,
)
    # Alg QW-RG step 6: Build split 1D raw-space blocks in the legacy style:
    # gausslet-gausslet from the existing gausslet side, Gaussian-Gaussian
    # analytically, and gausslet-Gaussian by a dedicated cross-block route.
    # See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
    gausslet_blocks = _qwrg_gausslet_1d_blocks(gausslet_bundle)
    gaussian_blocks = _qwrg_gaussian_analytic_blocks(supplement, expansion)
    cross_blocks = _qwrg_cross_1d_blocks(gausslet_bundle, supplement, expansion)
    return merge(gausslet_blocks, cross_blocks, gaussian_blocks)
end

function _qwrg_fill_product_column!(
    destination::AbstractVector{<:Real},
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    z::AbstractVector{<:Real},
)
    length(destination) == length(x) * length(y) * length(z) || throw(
        ArgumentError("Qiu-White separable 3D column assembly requires a matching destination length"),
    )
    index = 0
    @inbounds for ix in eachindex(x), iy in eachindex(y), iz in eachindex(z)
        index += 1
        destination[index] = x[ix] * y[iy] * z[iz]
    end
    return destination
end

function _qwrg_fill_product_matrix!(
    destination::AbstractMatrix{<:Real},
    x::AbstractMatrix{<:Real},
    y::AbstractMatrix{<:Real},
    z::AbstractMatrix{<:Real},
)
    size(destination, 1) == size(x, 1) * size(y, 1) * size(z, 1) || throw(
        ArgumentError("Qiu-White separable 3D matrix assembly requires matching row dimensions"),
    )
    size(destination, 2) == size(x, 2) * size(y, 2) * size(z, 2) || throw(
        ArgumentError("Qiu-White separable 3D matrix assembly requires matching column dimensions"),
    )
    row = 0
    @inbounds for ix in axes(x, 1), iy in axes(y, 1), iz in axes(z, 1)
        row += 1
        column = 0
        for jx in axes(x, 2), jy in axes(y, 2), jz in axes(z, 2)
            column += 1
            destination[row, column] =
                x[ix, jx] * y[iy, jy] * z[iz, jz]
        end
    end
    return destination
end

function _qwrg_raw_overlap_blocks(data)
    ngaussian = size(data.overlap_ga, 2)
    ngausslet3d = length(_mapped_cartesian_orbitals(1:size(data.overlap_gg, 1)))
    overlap_ga = zeros(Float64, ngausslet3d, ngaussian)
    overlap_aa = zeros(Float64, ngaussian, ngaussian)
    for index in 1:ngaussian
        overlap_vector = view(data.overlap_ga, :, index)
        _qwrg_fill_product_column!(view(overlap_ga, :, index), overlap_vector, overlap_vector, overlap_vector)
    end
    for i in 1:ngaussian, j in 1:ngaussian
        overlap_aa[i, j] = data.overlap_aa[i, j]^3
    end
    return overlap_ga, overlap_aa
end

function _qwrg_raw_kinetic_cross_block(data)
    ngaussian = size(data.overlap_ga, 2)
    ngausslet3d = size(data.overlap_ga, 1)^3
    kinetic_ga = zeros(Float64, ngausslet3d, ngaussian)
    scratch = zeros(Float64, ngausslet3d)
    for index in 1:ngaussian
        overlap_vector = view(data.overlap_ga, :, index)
        kinetic_vector = view(data.kinetic_ga, :, index)
        column = view(kinetic_ga, :, index)
        _qwrg_fill_product_column!(column, kinetic_vector, overlap_vector, overlap_vector)
        _qwrg_fill_product_column!(scratch, overlap_vector, kinetic_vector, overlap_vector)
        column .+= scratch
        _qwrg_fill_product_column!(scratch, overlap_vector, overlap_vector, kinetic_vector)
        column .+= scratch
    end
    return kinetic_ga
end

function _qwrg_raw_factor_cross_blocks(data, expansion::CoulombGaussianExpansion)
    ngaussian = size(data.overlap_ga, 2)
    ngausslet3d = size(data.overlap_ga, 1)^3
    factor_ga = Matrix{Float64}[]
    scratch = zeros(Float64, ngausslet3d)
    for term in eachindex(expansion.coefficients)
        matrix = zeros(Float64, ngausslet3d, ngaussian)
        for index in 1:ngaussian
            factor_vector = view(data.factor_ga[term], :, index)
            _qwrg_fill_product_column!(scratch, factor_vector, factor_vector, factor_vector)
            view(matrix, :, index) .= scratch
        end
        push!(factor_ga, matrix)
    end
    return factor_ga
end

function _qwrg_raw_one_body_aa_block(
    data,
    expansion::CoulombGaussianExpansion;
    Z::Real,
)
    ngaussian = size(data.overlap_aa, 1)
    one_body_aa = zeros(Float64, ngaussian, ngaussian)
    z_value = Float64(Z)
    for i in 1:ngaussian, j in i:ngaussian
        value = 3.0 * data.kinetic_aa[i, j] * data.overlap_aa[i, j]^2
        for term in eachindex(expansion.coefficients)
            value -= z_value * expansion.coefficients[term] * data.factor_aa[term][i, j]^3
        end
        one_body_aa[i, j] = value
        one_body_aa[j, i] = value
    end
    return one_body_aa
end

function _qwrg_raw_one_body_blocks(
    data,
    expansion::CoulombGaussianExpansion;
    Z::Real,
)
    ngaussian = size(data.overlap_ga, 2)
    ngausslet3d = size(data.overlap_ga, 1)^3
    one_body_ga = zeros(Float64, ngausslet3d, ngaussian)
    one_body_aa = zeros(Float64, ngaussian, ngaussian)
    z_value = Float64(Z)
    scratch = zeros(Float64, ngausslet3d)

    for index in 1:ngaussian
        overlap_vector = view(data.overlap_ga, :, index)
        kinetic_vector = view(data.kinetic_ga, :, index)
        column = view(one_body_ga, :, index)
        _qwrg_fill_product_column!(column, kinetic_vector, overlap_vector, overlap_vector)
        _qwrg_fill_product_column!(scratch, overlap_vector, kinetic_vector, overlap_vector)
        column .+= scratch
        _qwrg_fill_product_column!(scratch, overlap_vector, overlap_vector, kinetic_vector)
        column .+= scratch
        for term in eachindex(expansion.coefficients)
            factor_vector = view(data.factor_ga[term], :, index)
            _qwrg_fill_product_column!(scratch, factor_vector, factor_vector, factor_vector)
            column .-= z_value * expansion.coefficients[term] .* scratch
        end
    end

    one_body_aa .= _qwrg_raw_one_body_aa_block(data, expansion; Z = Z)

    return one_body_ga, one_body_aa
end

function _qwrg_gausslet_one_body_matrix(
    data,
    expansion::CoulombGaussianExpansion;
    Z::Real,
)
    overlap = data.overlap_gg
    kinetic = data.kinetic_gg
    n3 = size(overlap, 1)^3
    hamiltonian = zeros(Float64, n3, n3)
    scratch = zeros(Float64, n3, n3)
    _qwrg_fill_product_matrix!(hamiltonian, kinetic, overlap, overlap)
    _qwrg_fill_product_matrix!(scratch, overlap, kinetic, overlap)
    hamiltonian .+= scratch
    _qwrg_fill_product_matrix!(scratch, overlap, overlap, kinetic)
    hamiltonian .+= scratch
    hamiltonian .+= _mapped_coulomb_expanded_symmetric_matrix(
        -Float64(Z) .* expansion.coefficients,
        data.factor_gg_terms,
        data.factor_gg_terms,
        data.factor_gg_terms,
    )
    return 0.5 .* (hamiltonian .+ transpose(hamiltonian))
end

function _qwrg_gausslet_interaction_matrix(data, expansion::CoulombGaussianExpansion)
    return _mapped_coulomb_expanded_symmetric_matrix(
        expansion.coefficients,
        data.pair_gg_terms,
        data.pair_gg_terms,
        data.pair_gg_terms,
    )
end

function _qwrg_raw_axis_blocks(data, axis::Symbol; squared::Bool = false)
    gg_axis = squared ? data.x2_gg : data.position_gg
    ga_axis = squared ? data.x2_ga : data.position_ga
    aa_axis = squared ? data.x2_aa : data.position_aa
    overlap_gg = data.overlap_gg
    overlap_ga = data.overlap_ga
    overlap_aa = data.overlap_aa
    ngaussian = size(overlap_ga, 2)
    ngausslet3d = length(_mapped_cartesian_orbitals(1:size(overlap_gg, 1)))
    ga_block = zeros(Float64, ngausslet3d, ngaussian)
    aa_block = zeros(Float64, ngaussian, ngaussian)

    if axis == :x
        gg_block = zeros(Float64, ngausslet3d, ngausslet3d)
        _qwrg_fill_product_matrix!(gg_block, gg_axis, overlap_gg, overlap_gg)
        for index in 1:ngaussian
            _qwrg_fill_product_column!(
                view(ga_block, :, index),
                view(ga_axis, :, index),
                view(overlap_ga, :, index),
                view(overlap_ga, :, index),
            )
        end
    elseif axis == :y
        gg_block = zeros(Float64, ngausslet3d, ngausslet3d)
        _qwrg_fill_product_matrix!(gg_block, overlap_gg, gg_axis, overlap_gg)
        for index in 1:ngaussian
            _qwrg_fill_product_column!(
                view(ga_block, :, index),
                view(overlap_ga, :, index),
                view(ga_axis, :, index),
                view(overlap_ga, :, index),
            )
        end
    elseif axis == :z
        gg_block = zeros(Float64, ngausslet3d, ngausslet3d)
        _qwrg_fill_product_matrix!(gg_block, overlap_gg, overlap_gg, gg_axis)
        for index in 1:ngaussian
            _qwrg_fill_product_column!(
                view(ga_block, :, index),
                view(overlap_ga, :, index),
                view(overlap_ga, :, index),
                view(ga_axis, :, index),
            )
        end
    else
        throw(ArgumentError("axis must be :x, :y, or :z"))
    end

    for i in 1:ngaussian, j in 1:ngaussian
        aa_block[i, j] = aa_axis[i, j] * overlap_aa[i, j]^2
    end

    return gg_block, ga_block, aa_block
end

function _qwrg_gaussian_analytic_blocks(
    supplement,
    expansion::CoulombGaussianExpansion,
)
    gaussians, contraction_matrix = _qwrg_supplement_primitives_and_contraction(supplement)
    gaussian_set = PrimitiveSet1D(
        AbstractPrimitiveFunction1D[gaussian for gaussian in gaussians];
        name = :qiu_white_analytic_gaussians,
    )
    primitive_factors = Matrix{Float64}[
        Matrix{Float64}(gaussian_factor_matrix(
            gaussian_set;
            exponent = exponent,
            center = 0.0,
        )) for exponent in expansion.exponents
    ]
    primitive_pair_factors = _primitive_pair_gaussian_factor_matrices(
        gaussian_set,
        _select_primitive_matrix_backend(gaussian_set);
        exponents = expansion.exponents,
    )
    factor_terms = _congruence_transform_terms(_term_tensor(primitive_factors), contraction_matrix)
    pair_terms = _congruence_transform_terms(_term_tensor(primitive_pair_factors), contraction_matrix)
    return (
        overlap_aa = _congruence_transform(Matrix{Float64}(overlap_matrix(gaussian_set)), contraction_matrix),
        kinetic_aa = _congruence_transform(Matrix{Float64}(kinetic_matrix(gaussian_set)), contraction_matrix),
        position_aa = _congruence_transform(Matrix{Float64}(position_matrix(gaussian_set)), contraction_matrix),
        x2_aa = _congruence_transform(Matrix{Float64}(_x2_matrix(gaussian_set)), contraction_matrix),
        factor_aa = _tensor_to_matrix_vector(factor_terms),
        pair_aa = _tensor_to_matrix_vector(pair_terms),
    )
end

function _qwrg_diatomic_overlap_matrix(
    bundle_x::_MappedOrdinaryGausslet1DBundle,
    bundle_y::_MappedOrdinaryGausslet1DBundle,
    bundle_z::_MappedOrdinaryGausslet1DBundle,
)
    overlap_x = bundle_x.pgdg_intermediate.overlap
    overlap_y = bundle_y.pgdg_intermediate.overlap
    overlap_z = bundle_z.pgdg_intermediate.overlap
    matrix = zeros(
        Float64,
        size(overlap_x, 1) * size(overlap_y, 1) * size(overlap_z, 1),
        size(overlap_x, 2) * size(overlap_y, 2) * size(overlap_z, 2),
    )
    _qwrg_fill_product_matrix!(matrix, overlap_x, overlap_y, overlap_z)
    return matrix
end

function _qwrg_diatomic_interaction_matrix(
    bundle_x::_MappedOrdinaryGausslet1DBundle,
    bundle_y::_MappedOrdinaryGausslet1DBundle,
    bundle_z::_MappedOrdinaryGausslet1DBundle,
    expansion::CoulombGaussianExpansion,
)
    return _mapped_coulomb_expanded_symmetric_matrix(
        expansion.coefficients,
        bundle_x.pgdg_intermediate.pair_factor_terms,
        bundle_y.pgdg_intermediate.pair_factor_terms,
        bundle_z.pgdg_intermediate.pair_factor_terms,
    )
end

function _qwrg_diatomic_factor_term_cache(
    basis::MappedUniformBasis,
    centers_1d::AbstractVector{<:Real},
    expansion::CoulombGaussianExpansion,
    gausslet_backend::Symbol,
)
    cache = Dict{Float64,Array{Float64,3}}()
    for center_value in unique(Float64[Float64(value) for value in centers_1d])
        bundle = _mapped_ordinary_gausslet_1d_bundle(
            basis;
            exponents = expansion.exponents,
            center = center_value,
            backend = gausslet_backend,
        )
        cache[center_value] = bundle.pgdg_intermediate.gaussian_factor_terms
    end
    return cache
end

function _qwrg_full_cartesian_product_factorized_basis(dims::NTuple{3,Int})
    nx, ny, nz = dims
    x_functions = Matrix{Float64}(I, nx, nx)
    y_functions = Matrix{Float64}(I, ny, ny)
    z_functions = Matrix{Float64}(I, nz, nz)
    basis_triplets = NTuple{3,Int}[]
    sizehint!(basis_triplets, nx * ny * nz)
    for ix in 1:nx, iy in 1:ny, iz in 1:nz
        push!(basis_triplets, (ix, iy, iz))
    end
    return _CartesianNestedFactorizedBasis3D(
        dims,
        x_functions,
        y_functions,
        z_functions,
        basis_triplets,
        ones(Float64, length(basis_triplets)),
        0.0,
    )
end

function _qwrg_contracted_nuclear_axis_term_tables(
    axis_functions::AbstractMatrix{<:Real},
    basis::MappedUniformBasis,
    centers_1d::AbstractVector{<:Real},
    expansion::CoulombGaussianExpansion,
    gausslet_backend::Symbol,
)
    axis_matrix = Matrix{Float64}(axis_functions)
    center_values = Float64[Float64(value) for value in centers_1d]
    term_table_cache = Dict{Float64,Array{Float64,3}}()
    ordered_term_tables = Vector{Array{Float64,3}}(undef, length(center_values))
    @inbounds for center_index in eachindex(center_values)
        center_value = center_values[center_index]
        ordered_term_tables[center_index] = get!(term_table_cache, center_value) do
            bundle = _mapped_ordinary_gausslet_1d_bundle(
                basis;
                exponents = expansion.exponents,
                center = center_value,
                backend = gausslet_backend,
            )
            _nested_factorized_axis_term_tables(
                bundle.pgdg_intermediate.gaussian_factor_terms,
                axis_matrix,
            )
        end
    end
    return ordered_term_tables
end

function _qwrg_bond_aligned_nuclear_centers_by_axis(
    nuclei::AbstractVector{<:NTuple{3,<:Real}},
)
    nnuclei = length(nuclei)
    centers_x = Vector{Float64}(undef, nnuclei)
    centers_y = Vector{Float64}(undef, nnuclei)
    centers_z = Vector{Float64}(undef, nnuclei)
    @inbounds for nucleus_index in eachindex(nuclei)
        nucleus = nuclei[nucleus_index]
        centers_x[nucleus_index] = Float64(nucleus[1])
        centers_y[nucleus_index] = Float64(nucleus[2])
        centers_z[nucleus_index] = Float64(nucleus[3])
    end
    return centers_x, centers_y, centers_z
end

function _qwrg_factorized_basis_axis_indices(
    factorized_basis::_CartesianNestedFactorizedBasis3D,
)
    triplets = factorized_basis.basis_triplets
    nbasis = length(triplets)
    x_indices = Vector{Int}(undef, nbasis)
    y_indices = Vector{Int}(undef, nbasis)
    z_indices = Vector{Int}(undef, nbasis)
    @inbounds for basis_index in 1:nbasis
        ix, iy, iz = triplets[basis_index]
        x_indices[basis_index] = ix
        y_indices[basis_index] = iy
        z_indices[basis_index] = iz
    end
    return x_indices, y_indices, z_indices
end

function _qwrg_fill_direct_contracted_nuclear_matrix!(
    destination::Matrix{Float64},
    x_indices::Vector{Int},
    y_indices::Vector{Int},
    z_indices::Vector{Int},
    amplitudes::Vector{Float64},
    term_coefficients::Vector{Float64},
    operator_terms_x::Array{Float64,3},
    operator_terms_y::Array{Float64,3},
    operator_terms_z::Array{Float64,3},
)
    nbasis = length(x_indices)
    nterms = length(term_coefficients)
    size(destination) == (nbasis, nbasis) || throw(
        ArgumentError("direct contracted nuclear fill requires square output sized to the retained fixed basis"),
    )
    x_terms = vec(operator_terms_x)
    y_terms = vec(operator_terms_y)
    z_terms = vec(operator_terms_z)
    x_row_stride = stride(operator_terms_x, 2)
    y_row_stride = stride(operator_terms_y, 2)
    z_row_stride = stride(operator_terms_z, 2)
    x_column_stride = stride(operator_terms_x, 3)
    y_column_stride = stride(operator_terms_y, 3)
    z_column_stride = stride(operator_terms_z, 3)
    @inbounds for column in 1:nbasis
        xj = x_indices[column]
        yj = y_indices[column]
        zj = z_indices[column]
        amplitude_j = amplitudes[column]
        x_column_offset = (xj - 1) * x_column_stride
        y_column_offset = (yj - 1) * y_column_stride
        z_column_offset = (zj - 1) * z_column_stride
        for row in 1:column
            xi = x_indices[row]
            yi = y_indices[row]
            zi = z_indices[row]
            x_offset = (xi - 1) * x_row_stride + x_column_offset
            y_offset = (yi - 1) * y_row_stride + y_column_offset
            z_offset = (zi - 1) * z_row_stride + z_column_offset
            value = 0.0
            @simd for term in 1:nterms
                value +=
                    term_coefficients[term] *
                    x_terms[x_offset + term] *
                    y_terms[y_offset + term] *
                    z_terms[z_offset + term]
            end
            value *= amplitudes[row] * amplitude_j
            destination[row, column] = value
            destination[column, row] = value
        end
    end
    return destination
end

function _qwrg_bond_aligned_direct_contracted_nuclear_one_body_by_center(
    basis::AbstractBondAlignedOrdinaryQWBasis3D,
    factorized_basis::_CartesianNestedFactorizedBasis3D,
    bundle_x::_MappedOrdinaryGausslet1DBundle,
    bundle_y::_MappedOrdinaryGausslet1DBundle,
    bundle_z::_MappedOrdinaryGausslet1DBundle,
    expansion::CoulombGaussianExpansion,
    ;
    timing_setup_label::AbstractString = "qwrg.nuclear.direct_contracted.setup",
    timing_contract_label::AbstractString = "qwrg.nuclear.direct_contracted.contract",
)
    factorized_basis.dims == (
        size(bundle_x.pgdg_intermediate.overlap, 1),
        size(bundle_y.pgdg_intermediate.overlap, 1),
        size(bundle_z.pgdg_intermediate.overlap, 1),
    ) || throw(
        ArgumentError("direct contracted molecular nuclear assembly requires factorized-basis dimensions to match the parent Cartesian product basis"),
    )

    term_coefficients = Float64[-Float64(value) for value in expansion.coefficients]
    nuclei = basis.nuclei
    nnuclei = length(nuclei)
    nbasis = length(factorized_basis.basis_triplets)
    centers_x, centers_y, centers_z = _qwrg_bond_aligned_nuclear_centers_by_axis(nuclei)
    axis_term_tables_x, axis_term_tables_y, axis_term_tables_z = @timeg timing_setup_label begin
        axis_term_tables_x = _qwrg_contracted_nuclear_axis_term_tables(
            factorized_basis.x_functions,
            basis.basis_x,
            centers_x,
            expansion,
            bundle_x.backend,
        )
        axis_term_tables_y = _qwrg_contracted_nuclear_axis_term_tables(
            factorized_basis.y_functions,
            basis.basis_y,
            centers_y,
            expansion,
            bundle_y.backend,
        )
        axis_term_tables_z = _qwrg_contracted_nuclear_axis_term_tables(
            factorized_basis.z_functions,
            basis.basis_z,
            centers_z,
            expansion,
            bundle_z.backend,
        )
        (axis_term_tables_x, axis_term_tables_y, axis_term_tables_z)
    end

    return @timeg timing_contract_label begin
        x_indices, y_indices, z_indices = _qwrg_factorized_basis_axis_indices(factorized_basis)
        amplitudes = factorized_basis.basis_amplitudes
        matrices = Vector{Matrix{Float64}}(undef, nnuclei)
        for nucleus_index in 1:nnuclei
            matrix = Matrix{Float64}(undef, nbasis, nbasis)
            _qwrg_fill_direct_contracted_nuclear_matrix!(
                matrix,
                x_indices,
                y_indices,
                z_indices,
                amplitudes,
                term_coefficients,
                axis_term_tables_x[nucleus_index],
                axis_term_tables_y[nucleus_index],
                axis_term_tables_z[nucleus_index],
            )
            matrices[nucleus_index] = matrix
        end
        matrices
    end
end

function _qwrg_diatomic_kinetic_matrix(
    bundle_x::_MappedOrdinaryGausslet1DBundle,
    bundle_y::_MappedOrdinaryGausslet1DBundle,
    bundle_z::_MappedOrdinaryGausslet1DBundle,
)
    overlap_x = bundle_x.pgdg_intermediate.overlap
    overlap_y = bundle_y.pgdg_intermediate.overlap
    overlap_z = bundle_z.pgdg_intermediate.overlap
    kinetic_x = bundle_x.pgdg_intermediate.kinetic
    kinetic_y = bundle_y.pgdg_intermediate.kinetic
    kinetic_z = bundle_z.pgdg_intermediate.kinetic

    matrix = zeros(
        Float64,
        size(overlap_x, 1) * size(overlap_y, 1) * size(overlap_z, 1),
        size(overlap_x, 2) * size(overlap_y, 2) * size(overlap_z, 2),
    )
    scratch = similar(matrix)
    _qwrg_fill_product_matrix!(matrix, kinetic_x, overlap_y, overlap_z)
    _qwrg_fill_product_matrix!(scratch, overlap_x, kinetic_y, overlap_z)
    matrix .+= scratch
    _qwrg_fill_product_matrix!(scratch, overlap_x, overlap_y, kinetic_z)
    matrix .+= scratch
    return 0.5 .* (matrix .+ transpose(matrix))
end

function _qwrg_diatomic_nuclear_one_body_by_center(
    basis::AbstractBondAlignedOrdinaryQWBasis3D,
    bundle_x::_MappedOrdinaryGausslet1DBundle,
    bundle_y::_MappedOrdinaryGausslet1DBundle,
    bundle_z::_MappedOrdinaryGausslet1DBundle,
    expansion::CoulombGaussianExpansion,
)
    factor_x = _qwrg_diatomic_factor_term_cache(
        basis.basis_x,
        [nucleus[1] for nucleus in basis.nuclei],
        expansion,
        bundle_x.backend,
    )
    factor_y = _qwrg_diatomic_factor_term_cache(
        basis.basis_y,
        [nucleus[2] for nucleus in basis.nuclei],
        expansion,
        bundle_y.backend,
    )
    factor_z = _qwrg_diatomic_factor_term_cache(
        basis.basis_z,
        [nucleus[3] for nucleus in basis.nuclei],
        expansion,
        bundle_z.backend,
    )

    matrices = Matrix{Float64}[]
    for nucleus in basis.nuclei
        matrix = _mapped_coulomb_expanded_symmetric_matrix(
            -expansion.coefficients,
            factor_x[nucleus[1]],
            factor_y[nucleus[2]],
            factor_z[nucleus[3]],
        )
        push!(matrices, matrix)
    end
    return matrices
end

function _qwrg_same_nuclei(
    left::AbstractVector{<:NTuple{3,<:Real}},
    right::AbstractVector{<:NTuple{3,<:Real}};
    tol::Real = 1.0e-12,
)
    length(left) == length(right) || return false
    for index in eachindex(left, right)
        for axis in 1:3
            abs(Float64(left[index][axis]) - Float64(right[index][axis])) <= tol || return false
        end
    end
    return true
end

function _qwrg_diatomic_gausslet_axis_matrix(
    bundle_x::_MappedOrdinaryGausslet1DBundle,
    bundle_y::_MappedOrdinaryGausslet1DBundle,
    bundle_z::_MappedOrdinaryGausslet1DBundle,
    axis::Symbol;
    squared::Bool = false,
)
    overlap_x = bundle_x.pgdg_intermediate.overlap
    overlap_y = bundle_y.pgdg_intermediate.overlap
    overlap_z = bundle_z.pgdg_intermediate.overlap
    axis_x = squared ? bundle_x.pgdg_intermediate.x2 : bundle_x.pgdg_intermediate.position
    axis_y = squared ? bundle_y.pgdg_intermediate.x2 : bundle_y.pgdg_intermediate.position
    axis_z = squared ? bundle_z.pgdg_intermediate.x2 : bundle_z.pgdg_intermediate.position
    matrix = zeros(
        Float64,
        size(overlap_x, 1) * size(overlap_y, 1) * size(overlap_z, 1),
        size(overlap_x, 2) * size(overlap_y, 2) * size(overlap_z, 2),
    )
    if axis == :x
        _qwrg_fill_product_matrix!(matrix, axis_x, overlap_y, overlap_z)
    elseif axis == :y
        _qwrg_fill_product_matrix!(matrix, overlap_x, axis_y, overlap_z)
    elseif axis == :z
        _qwrg_fill_product_matrix!(matrix, overlap_x, overlap_y, axis_z)
    else
        throw(ArgumentError("axis must be :x, :y, or :z"))
    end
    return matrix
end

function _qwrg_diatomic_cartesian_shell_context(
    bundles::_CartesianNestedAxisBundles3D,
    supplement::_BondAlignedDiatomicCartesianShellSupplement3D,
    basis::BondAlignedDiatomicQWBasis3D,
)
    _qwrg_same_nuclei(supplement.source.nuclei, basis.nuclei) || throw(
        ArgumentError("bond-aligned diatomic molecular supplement nuclei must match the basis nuclei"),
    )

    return (
        proxy_x = _qwrg_diatomic_supplement_proxy_layer(basis.basis_x, bundles.bundle_x, :x),
        proxy_y = _qwrg_diatomic_supplement_proxy_layer(basis.basis_y, bundles.bundle_y, :y),
        proxy_z = _qwrg_diatomic_supplement_proxy_layer(basis.basis_z, bundles.bundle_z, :z),
        ngausslet3d =
            size(bundles.bundle_x.pgdg_intermediate.overlap, 1) *
            size(bundles.bundle_y.pgdg_intermediate.overlap, 1) *
            size(bundles.bundle_z.pgdg_intermediate.overlap, 1),
        norbital = length(supplement.orbitals),
    )
end

function _qwrg_diatomic_cartesian_shell_overlap_blocks_3d(
    bundles::_CartesianNestedAxisBundles3D,
    supplement::_BondAlignedDiatomicCartesianShellSupplement3D,
    basis::BondAlignedDiatomicQWBasis3D,
    expansion::CoulombGaussianExpansion,
)
    context = _qwrg_diatomic_cartesian_shell_context(bundles, supplement, basis)
    overlap_ga = zeros(Float64, context.ngausslet3d, context.norbital)
    overlap_aa = zeros(Float64, context.norbital, context.norbital)

    cross_cache = Dict{Tuple{Int,Symbol},Any}()
    aa_cache = Dict{Tuple{Int,Int,Symbol},Any}()
    scratch = zeros(Float64, context.ngausslet3d)

    for (orbital_index, orbital) in pairs(supplement.orbitals)
        x_data = get!(cross_cache, (orbital_index, :x)) do
            _qwrg_atomic_axis_cross_data(context.proxy_x, orbital, :x, expansion)
        end
        y_data = get!(cross_cache, (orbital_index, :y)) do
            _qwrg_atomic_axis_cross_data(context.proxy_y, orbital, :y, expansion)
        end
        z_data = get!(cross_cache, (orbital_index, :z)) do
            _qwrg_atomic_axis_cross_data(context.proxy_z, orbital, :z, expansion)
        end

        for primitive in eachindex(orbital.coefficients)
            coefficient = Float64(orbital.coefficients[primitive])
            _qwrg_fill_product_column!(
                scratch,
                view(x_data.overlap, :, primitive),
                view(y_data.overlap, :, primitive),
                view(z_data.overlap, :, primitive),
            )
            overlap_ga[:, orbital_index] .+= coefficient .* scratch
        end
    end

    for (left_index, left) in pairs(supplement.orbitals), (right_index, right) in pairs(supplement.orbitals)
        x_data = get!(aa_cache, (left_index, right_index, :x)) do
            _qwrg_atomic_axis_aa_one_body_data(left, right, :x, expansion)
        end
        y_data = get!(aa_cache, (left_index, right_index, :y)) do
            _qwrg_atomic_axis_aa_one_body_data(left, right, :y, expansion)
        end
        z_data = get!(aa_cache, (left_index, right_index, :z)) do
            _qwrg_atomic_axis_aa_one_body_data(left, right, :z, expansion)
        end
        overlap_aa[left_index, right_index] = _qwrg_atomic_weighted_hadamard(
            left.coefficients,
            x_data.overlap,
            y_data.overlap,
            z_data.overlap,
            right.coefficients,
        )
    end

    return (
        overlap_ga = overlap_ga,
        overlap_aa = Matrix{Float64}(0.5 .* (overlap_aa .+ transpose(overlap_aa))),
    )
end

function _qwrg_diatomic_cartesian_shell_blocks_3d(
    bundles::_CartesianNestedAxisBundles3D,
    supplement::_BondAlignedDiatomicCartesianShellSupplement3D,
    basis::BondAlignedDiatomicQWBasis3D,
    expansion::CoulombGaussianExpansion,
    nuclear_charges::AbstractVector{<:Real},
)
    length(nuclear_charges) == length(basis.nuclei) || throw(
        ArgumentError("bond-aligned diatomic molecular supplement route requires one nuclear charge per nucleus"),
    )
    context = _qwrg_diatomic_cartesian_shell_context(bundles, supplement, basis)
    proxy_x = context.proxy_x
    proxy_y = context.proxy_y
    proxy_z = context.proxy_z
    ngausslet3d = context.ngausslet3d
    norbital = context.norbital

    overlap_ga = zeros(Float64, ngausslet3d, norbital)
    kinetic_ga = zeros(Float64, ngausslet3d, norbital)
    one_body_ga = zeros(Float64, ngausslet3d, norbital)
    position_x_ga = zeros(Float64, ngausslet3d, norbital)
    position_y_ga = zeros(Float64, ngausslet3d, norbital)
    position_z_ga = zeros(Float64, ngausslet3d, norbital)
    x2_x_ga = zeros(Float64, ngausslet3d, norbital)
    x2_y_ga = zeros(Float64, ngausslet3d, norbital)
    x2_z_ga = zeros(Float64, ngausslet3d, norbital)

    overlap_aa = zeros(Float64, norbital, norbital)
    kinetic_aa = zeros(Float64, norbital, norbital)
    one_body_aa = zeros(Float64, norbital, norbital)
    position_x_aa = zeros(Float64, norbital, norbital)
    position_y_aa = zeros(Float64, norbital, norbital)
    position_z_aa = zeros(Float64, norbital, norbital)
    x2_x_aa = zeros(Float64, norbital, norbital)
    x2_y_aa = zeros(Float64, norbital, norbital)
    x2_z_aa = zeros(Float64, norbital, norbital)
    nuclear_ga_by_center = [zeros(Float64, ngausslet3d, norbital) for _ in basis.nuclei]
    nuclear_aa_by_center = [zeros(Float64, norbital, norbital) for _ in basis.nuclei]

    cross_cache = Dict{Tuple{Int,Symbol},Any}()
    factor_cross_cache = Dict{Tuple{Int,Symbol,Int},Any}()
    aa_cache = Dict{Tuple{Int,Int,Symbol},Any}()
    scratch = zeros(Float64, ngausslet3d)

    # Alg QW-RG step 3: build the added 3D molecular Gaussian orbitals as one
    # explicit two-center Cartesian shell set on the bond-aligned nuclei.
    # See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
    for (orbital_index, orbital) in pairs(supplement.orbitals)
        x_data = get!(cross_cache, (orbital_index, :x)) do
            _qwrg_atomic_axis_cross_data(proxy_x, orbital, :x, expansion)
        end
        y_data = get!(cross_cache, (orbital_index, :y)) do
            _qwrg_atomic_axis_cross_data(proxy_y, orbital, :y, expansion)
        end
        z_data = get!(cross_cache, (orbital_index, :z)) do
            _qwrg_atomic_axis_cross_data(proxy_z, orbital, :z, expansion)
        end

        for primitive in eachindex(orbital.coefficients)
            coefficient = Float64(orbital.coefficients[primitive])
            _qwrg_fill_product_column!(
                scratch,
                view(x_data.overlap, :, primitive),
                view(y_data.overlap, :, primitive),
                view(z_data.overlap, :, primitive),
            )
            overlap_ga[:, orbital_index] .+= coefficient .* scratch

            _qwrg_fill_product_column!(
                scratch,
                view(x_data.kinetic, :, primitive),
                view(y_data.overlap, :, primitive),
                view(z_data.overlap, :, primitive),
            )
            kinetic_ga[:, orbital_index] .+= coefficient .* scratch
            _qwrg_fill_product_column!(
                scratch,
                view(x_data.overlap, :, primitive),
                view(y_data.kinetic, :, primitive),
                view(z_data.overlap, :, primitive),
            )
            kinetic_ga[:, orbital_index] .+= coefficient .* scratch
            _qwrg_fill_product_column!(
                scratch,
                view(x_data.overlap, :, primitive),
                view(y_data.overlap, :, primitive),
                view(z_data.kinetic, :, primitive),
            )
            kinetic_ga[:, orbital_index] .+= coefficient .* scratch

            _qwrg_fill_product_column!(
                scratch,
                view(x_data.position, :, primitive),
                view(y_data.overlap, :, primitive),
                view(z_data.overlap, :, primitive),
            )
            position_x_ga[:, orbital_index] .+= coefficient .* scratch
            _qwrg_fill_product_column!(
                scratch,
                view(x_data.overlap, :, primitive),
                view(y_data.position, :, primitive),
                view(z_data.overlap, :, primitive),
            )
            position_y_ga[:, orbital_index] .+= coefficient .* scratch
            _qwrg_fill_product_column!(
                scratch,
                view(x_data.overlap, :, primitive),
                view(y_data.overlap, :, primitive),
                view(z_data.position, :, primitive),
            )
            position_z_ga[:, orbital_index] .+= coefficient .* scratch

            _qwrg_fill_product_column!(
                scratch,
                view(x_data.x2, :, primitive),
                view(y_data.overlap, :, primitive),
                view(z_data.overlap, :, primitive),
            )
            x2_x_ga[:, orbital_index] .+= coefficient .* scratch
            _qwrg_fill_product_column!(
                scratch,
                view(x_data.overlap, :, primitive),
                view(y_data.x2, :, primitive),
                view(z_data.overlap, :, primitive),
            )
            x2_y_ga[:, orbital_index] .+= coefficient .* scratch
            _qwrg_fill_product_column!(
                scratch,
                view(x_data.overlap, :, primitive),
                view(y_data.overlap, :, primitive),
                view(z_data.x2, :, primitive),
            )
            x2_z_ga[:, orbital_index] .+= coefficient .* scratch
        end

        one_body_ga[:, orbital_index] .= kinetic_ga[:, orbital_index]
        for (nucleus_index, nucleus) in pairs(basis.nuclei)
            charge_value = Float64(nuclear_charges[nucleus_index])
            x_factor = get!(factor_cross_cache, (orbital_index, :x, nucleus_index)) do
                _qwrg_atomic_axis_factor_cross_data(
                    proxy_x,
                    orbital,
                    :x,
                    expansion,
                    nucleus[1],
                )
            end
            y_factor = get!(factor_cross_cache, (orbital_index, :y, nucleus_index)) do
                _qwrg_atomic_axis_factor_cross_data(
                    proxy_y,
                    orbital,
                    :y,
                    expansion,
                    nucleus[2],
                )
            end
            z_factor = get!(factor_cross_cache, (orbital_index, :z, nucleus_index)) do
                _qwrg_atomic_axis_factor_cross_data(
                    proxy_z,
                    orbital,
                    :z,
                    expansion,
                    nucleus[3],
                )
            end
            for term in eachindex(expansion.coefficients), primitive in eachindex(orbital.coefficients)
                coefficient = Float64(orbital.coefficients[primitive])
                _qwrg_fill_product_column!(
                    scratch,
                    view(x_factor[term], :, primitive),
                    view(y_factor[term], :, primitive),
                    view(z_factor[term], :, primitive),
                )
                nuclear_ga_by_center[nucleus_index][:, orbital_index] .-=
                    expansion.coefficients[term] .* coefficient .* scratch
            end
            one_body_ga[:, orbital_index] .+=
                charge_value .* nuclear_ga_by_center[nucleus_index][:, orbital_index]
        end
    end

    for (left_index, left) in pairs(supplement.orbitals), (right_index, right) in pairs(supplement.orbitals)
        x_data = get!(aa_cache, (left_index, right_index, :x)) do
            _qwrg_atomic_axis_aa_one_body_data(left, right, :x, expansion)
        end
        y_data = get!(aa_cache, (left_index, right_index, :y)) do
            _qwrg_atomic_axis_aa_one_body_data(left, right, :y, expansion)
        end
        z_data = get!(aa_cache, (left_index, right_index, :z)) do
            _qwrg_atomic_axis_aa_one_body_data(left, right, :z, expansion)
        end
        overlap_aa[left_index, right_index] = _qwrg_atomic_weighted_hadamard(
            left.coefficients,
            x_data.overlap,
            y_data.overlap,
            z_data.overlap,
            right.coefficients,
        )
        kinetic_aa[left_index, right_index] =
            _qwrg_atomic_weighted_hadamard(
                left.coefficients,
                x_data.kinetic,
                y_data.overlap,
                z_data.overlap,
                right.coefficients,
            ) +
            _qwrg_atomic_weighted_hadamard(
                left.coefficients,
                x_data.overlap,
                y_data.kinetic,
                z_data.overlap,
                right.coefficients,
            ) +
            _qwrg_atomic_weighted_hadamard(
                left.coefficients,
                x_data.overlap,
                y_data.overlap,
                z_data.kinetic,
                right.coefficients,
            )
        position_x_aa[left_index, right_index] = _qwrg_atomic_weighted_hadamard(
            left.coefficients,
            x_data.position,
            y_data.overlap,
            z_data.overlap,
            right.coefficients,
        )
        position_y_aa[left_index, right_index] = _qwrg_atomic_weighted_hadamard(
            left.coefficients,
            x_data.overlap,
            y_data.position,
            z_data.overlap,
            right.coefficients,
        )
        position_z_aa[left_index, right_index] = _qwrg_atomic_weighted_hadamard(
            left.coefficients,
            x_data.overlap,
            y_data.overlap,
            z_data.position,
            right.coefficients,
        )
        x2_x_aa[left_index, right_index] = _qwrg_atomic_weighted_hadamard(
            left.coefficients,
            x_data.x2,
            y_data.overlap,
            z_data.overlap,
            right.coefficients,
        )
        x2_y_aa[left_index, right_index] = _qwrg_atomic_weighted_hadamard(
            left.coefficients,
            x_data.overlap,
            y_data.x2,
            z_data.overlap,
            right.coefficients,
        )
        x2_z_aa[left_index, right_index] = _qwrg_atomic_weighted_hadamard(
            left.coefficients,
            x_data.overlap,
            y_data.overlap,
            z_data.x2,
            right.coefficients,
        )

        one_body_value = kinetic_aa[left_index, right_index]
        factor_aa_local = Dict{Tuple{Symbol,Float64},Any}()
        for (nucleus_index, nucleus) in pairs(basis.nuclei)
            charge_value = Float64(nuclear_charges[nucleus_index])
            x_center = Float64(nucleus[1])
            x_factor = get!(factor_aa_local, (:x, x_center)) do
                _qwrg_atomic_axis_factor_aa_data(left, right, :x, expansion, x_center)
            end
            y_center = Float64(nucleus[2])
            y_factor = get!(factor_aa_local, (:y, y_center)) do
                _qwrg_atomic_axis_factor_aa_data(left, right, :y, expansion, y_center)
            end
            z_center = Float64(nucleus[3])
            z_factor = get!(factor_aa_local, (:z, z_center)) do
                _qwrg_atomic_axis_factor_aa_data(left, right, :z, expansion, z_center)
            end
            center_value = 0.0
            for term in eachindex(expansion.coefficients)
                center_value -= expansion.coefficients[term] * _qwrg_atomic_weighted_hadamard(
                    left.coefficients,
                    x_factor[term],
                    y_factor[term],
                    z_factor[term],
                    right.coefficients,
                )
            end
            nuclear_aa_by_center[nucleus_index][left_index, right_index] = center_value
            one_body_value += charge_value * center_value
        end
        one_body_aa[left_index, right_index] = one_body_value
    end

    return (
        overlap_ga = overlap_ga,
        overlap_aa = Matrix{Float64}(0.5 .* (overlap_aa .+ transpose(overlap_aa))),
        kinetic_ga = kinetic_ga,
        kinetic_aa = Matrix{Float64}(0.5 .* (kinetic_aa .+ transpose(kinetic_aa))),
        one_body_ga = one_body_ga,
        one_body_aa = Matrix{Float64}(0.5 .* (one_body_aa .+ transpose(one_body_aa))),
        nuclear_ga_by_center = nuclear_ga_by_center,
        nuclear_aa_by_center =
            [Matrix{Float64}(0.5 .* (matrix .+ transpose(matrix))) for matrix in nuclear_aa_by_center],
        position_x_ga = position_x_ga,
        position_y_ga = position_y_ga,
        position_z_ga = position_z_ga,
        position_x_aa = Matrix{Float64}(0.5 .* (position_x_aa .+ transpose(position_x_aa))),
        position_y_aa = Matrix{Float64}(0.5 .* (position_y_aa .+ transpose(position_y_aa))),
        position_z_aa = Matrix{Float64}(0.5 .* (position_z_aa .+ transpose(position_z_aa))),
        x2_x_ga = x2_x_ga,
        x2_y_ga = x2_y_ga,
        x2_z_ga = x2_z_ga,
        x2_x_aa = Matrix{Float64}(0.5 .* (x2_x_aa .+ transpose(x2_x_aa))),
        x2_y_aa = Matrix{Float64}(0.5 .* (x2_y_aa .+ transpose(x2_y_aa))),
        x2_z_aa = Matrix{Float64}(0.5 .* (x2_z_aa .+ transpose(x2_z_aa))),
    )
end

function _qwrg_contract_parent_ga_matrix(
    contraction::AbstractMatrix{<:Real},
    parent_ga::AbstractMatrix{<:Real},
)
    size(contraction, 1) == size(parent_ga, 1) || throw(
        DimensionMismatch("nested fixed-block contraction requires parent fixed rows to match the parent raw fixed-Gaussian block"),
    )
    return Matrix{Float64}(transpose(contraction) * parent_ga)
end

function _qwrg_contract_parent_symmetric_matrix(
    contraction::AbstractMatrix{<:Real},
    parent_matrix::AbstractMatrix{<:Real},
)
    contracted = Matrix{Float64}(transpose(contraction) * parent_matrix * contraction)
    return Matrix{Float64}(0.5 .* (contracted .+ transpose(contracted)))
end

function _qwrg_contract_parent_ga_terms(
    contraction::AbstractMatrix{<:Real},
    parent_terms::AbstractVector{<:AbstractMatrix{<:Real}},
)
    return [_qwrg_contract_parent_ga_matrix(contraction, term) for term in parent_terms]
end

function _qwrg_fixed_block_one_body_matrix(
    fixed_block::_NestedFixedBlock3D,
    expansion::CoulombGaussianExpansion;
    Z::Real,
)
    hamiltonian = Matrix{Float64}(fixed_block.kinetic)
    isnothing(fixed_block.gaussian_sum) && throw(
        ArgumentError("nested fixed-block QW one-body assembly requires gaussian_sum"),
    )
    hamiltonian .-= Float64(Z) .* fixed_block.gaussian_sum
    return Matrix{Float64}(0.5 .* (hamiltonian .+ transpose(hamiltonian)))
end

function _qwrg_fixed_block_interaction_matrix(
    fixed_block::_NestedFixedBlock3D,
    expansion::CoulombGaussianExpansion,
)
    isnothing(fixed_block.pair_sum) && throw(
        ArgumentError("nested fixed-block QW interaction assembly requires pair_sum"),
    )
    interaction = Matrix{Float64}(fixed_block.pair_sum)
    return Matrix{Float64}(0.5 .* (interaction .+ transpose(interaction)))
end
