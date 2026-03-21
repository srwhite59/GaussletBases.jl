"""
    QiuWhiteHybridOrbital3D

One orbital index in the paper-faithful Qiu-White residual-Gaussian reference
path.

The object records whether the orbital is a Cartesian gausslet product orbital
or a residual Gaussian together with the associated center and, for residual
Gaussians, the matched widths used in the MWG diagnostics.
"""
struct QiuWhiteHybridOrbital3D
    index::Int
    kind::Symbol
    label::String
    x::Float64
    y::Float64
    z::Float64
    wx::Float64
    wy::Float64
    wz::Float64
end

function Base.show(io::IO, orbital::QiuWhiteHybridOrbital3D)
    print(
        io,
        "QiuWhiteHybridOrbital3D(index=",
        orbital.index,
        ", kind=:",
        orbital.kind,
        ", label=\"",
        orbital.label,
        "\", center=(",
        orbital.x,
        ", ",
        orbital.y,
        ", ",
        orbital.z,
        ")",
    )
    if orbital.kind == :residual_gaussian
        print(
            io,
            ", widths=(",
            orbital.wx,
            ", ",
            orbital.wy,
            ", ",
            orbital.wz,
            ")",
        )
    end
    print(io, ")")
end

"""
    QiuWhiteResidualGaussianOperators

Paper-faithful Qiu-White residual-Gaussian ordinary Cartesian reference
Hamiltonian.

This object keeps the final basis as the full 3D gausslet product basis plus
orthonormalized 3D residual Gaussians. The one-body matrices are built exactly
in the raw gausslet-plus-GTO space and transformed into the final basis, while
the two-electron interaction stays in the same two-index integral-diagonal
approximation (IDA) representation used for the gausslet channel.
"""
struct QiuWhiteResidualGaussianOperators{B,D}
    basis::B
    gaussian_data::D
    gausslet_backend::Symbol
    interaction_treatment::Symbol
    expansion::CoulombGaussianExpansion
    overlap::Matrix{Float64}
    one_body_hamiltonian::Matrix{Float64}
    interaction_matrix::Matrix{Float64}
    orbital_data::Vector{QiuWhiteHybridOrbital3D}
    gausslet_count::Int
    residual_count::Int
    raw_to_final::Matrix{Float64}
    residual_centers::Matrix{Float64}
    residual_widths::Matrix{Float64}
end

function _qwrg_elapsed_seconds(start_ns::UInt64)
    return (time_ns() - start_ns) / 1.0e9
end

function _qwrg_maybe_print_timings(
    io::IO,
    timings::AbstractVector{<:Pair{<:AbstractString,<:Real}},
)
    println(io, "Qiu-White reference timings")
    total = 0.0
    for timing in timings
        value = Float64(timing.second)
        total += value
        println(io, "  ", timing.first, ": ", value, " s")
    end
    println(io, "  total: ", total, " s")
    return nothing
end

function _qwrg_record_timing!(
    io::IO,
    timings::Vector{Pair{String,Float64}},
    label::AbstractString,
    start_ns::UInt64,
)
    value = _qwrg_elapsed_seconds(start_ns)
    push!(timings, String(label) => value)
    println(io, "QW-RG timing  ", label, ": ", value, " s")
    flush(io)
    return value
end

function Base.show(io::IO, operators::QiuWhiteResidualGaussianOperators)
    print(
        io,
        "QiuWhiteResidualGaussianOperators(gausslet_backend=:",
        operators.gausslet_backend,
        ", interaction=:",
        operators.interaction_treatment,
        ", ngausslet=",
        operators.gausslet_count,
        ", nresidual=",
        operators.residual_count,
        ", reference=true)",
    )
end

orbitals(operators::QiuWhiteResidualGaussianOperators) = operators.orbital_data

function ordinary_cartesian_vee_expectation(
    operators::QiuWhiteResidualGaussianOperators,
    orbital::AbstractVector;
    overlap_tol::Real = 1.0e-8,
)
    length(orbital) == length(operators.orbital_data) ||
        throw(ArgumentError("orbital length must match the Qiu-White ordinary Cartesian orbital dimension"))
    overlap_error = norm(operators.overlap - I, Inf)
    overlap_error <= Float64(overlap_tol) || throw(
        ArgumentError(
            "ordinary_cartesian_vee_expectation currently requires an orthonormal final basis; got overlap error $(overlap_error)",
        ),
    )
    weights = Float64[abs2(coefficient) for coefficient in orbital]
    norm2 = sum(weights)
    norm2 > 0.0 || throw(ArgumentError("orbital must have nonzero norm"))
    weights ./= norm2
    return Float64(real(dot(weights, operators.interaction_matrix * weights)))
end

function ordinary_cartesian_1s2_check(
    operators::QiuWhiteResidualGaussianOperators;
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
        overlap_error = norm(operators.overlap - I, Inf),
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

function _qwrg_gausslet_1d_blocks(bundle::_MappedOrdinaryGausslet1DBundle)
    return (
        overlap_gg = bundle.overlap,
        kinetic_gg = bundle.kinetic,
        position_gg = bundle.position,
        x2_gg = bundle.x2,
        factor_gg = bundle.gaussian_factors,
        pair_gg = bundle.pair_factors,
        weight_gg = bundle.weights,
        center_gg = bundle.centers,
    )
end

function _qwrg_cross_1d_blocks(
    basis::MappedUniformBasis,
    gaussians::AbstractVector{<:Gaussian},
    expansion::CoulombGaussianExpansion;
    h = nothing,
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
        distance2 = (points .- transpose(points)) .^ 2
        pair_factors = Matrix{Float64}[]
        for exponent in exponent_values
            kernel = exp.(-exponent .* distance2)
            push!(pair_factors, Matrix{Float64}(transpose(weighted_basis) * (kernel * weighted_gaussians)))
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
            maxdiff = maximum(
                vcat(
                    norm(current.overlap_ga - previous.overlap_ga, Inf),
                    norm(current.kinetic_ga - previous.kinetic_ga, Inf),
                    norm(current.position_ga - previous.position_ga, Inf),
                    norm(current.x2_ga - previous.x2_ga, Inf),
                    [norm(current.factor_ga[index] - previous.factor_ga[index], Inf) for index in eachindex(current.factor_ga)]...,
                    [norm(current.pair_ga[index] - previous.pair_ga[index], Inf) for index in eachindex(current.pair_ga)]...,
                ),
            )
            maxdiff <= _PRIMITIVE_MATRIX_TOL && return current
        end
        previous = current
        h_try /= 2.0
    end
    return current
end

function _qwrg_split_block_matrices(
    gausslet_bundle::_MappedOrdinaryGausslet1DBundle,
    gaussians::AbstractVector{<:Gaussian},
    expansion::CoulombGaussianExpansion,
)
    # Alg QW-RG step 6: Build split 1D raw-space blocks in the legacy style:
    # gausslet-gausslet from the existing gausslet side, Gaussian-Gaussian
    # analytically, and gausslet-Gaussian by a dedicated cross-block route.
    # See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
    gausslet_blocks = _qwrg_gausslet_1d_blocks(gausslet_bundle)
    gaussian_blocks = _qwrg_gaussian_analytic_blocks(gaussians, expansion)
    cross_blocks = _qwrg_cross_1d_blocks(gausslet_bundle.basis, gaussians, expansion)
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

    for i in 1:ngaussian, j in i:ngaussian
        value = 3.0 * data.kinetic_aa[i, j] * data.overlap_aa[i, j]^2
        for term in eachindex(expansion.coefficients)
            value -= z_value * expansion.coefficients[term] * data.factor_aa[term][i, j]^3
        end
        one_body_aa[i, j] = value
        one_body_aa[j, i] = value
    end

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
    for term in eachindex(expansion.coefficients)
        factor = data.factor_gg[term]
        _qwrg_fill_product_matrix!(scratch, factor, factor, factor)
        hamiltonian .-= Float64(Z) * expansion.coefficients[term] .* scratch
    end
    return 0.5 .* (hamiltonian .+ transpose(hamiltonian))
end

function _qwrg_gausslet_interaction_matrix(data, expansion::CoulombGaussianExpansion)
    weight_outer = data.weight_gg * transpose(data.weight_gg)
    pair_factors = [factor ./ weight_outer for factor in data.pair_gg]
    ngausslet3d = length(_mapped_cartesian_orbitals(1:size(data.overlap_gg, 1)))
    interaction = zeros(Float64, ngausslet3d, ngausslet3d)
    scratch = zeros(Float64, ngausslet3d, ngausslet3d)
    for term in eachindex(expansion.coefficients)
        factor = pair_factors[term]
        _qwrg_fill_product_matrix!(scratch, factor, factor, factor)
        interaction .+= expansion.coefficients[term] .* scratch
    end
    return 0.5 .* (interaction .+ transpose(interaction))
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

# Alg QW-RG step 4: Define residual Gaussians by orthogonalizing 3D GTOs
# to the full 3D gausslet space.
# See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
function _qwrg_residual_space(
    gausslet_overlap::AbstractMatrix{<:Real},
    overlap_ga::AbstractMatrix{<:Real},
    overlap_aa::AbstractMatrix{<:Real},
)
    gausslet_overlap_value = Matrix{Float64}(gausslet_overlap)
    gausslet_overlap_error = norm(gausslet_overlap_value - I, Inf)
    gausslet_overlap_error <= 1.0e-8 || throw(
        ArgumentError(
            "Qiu-White reference path requires the fixed gausslet 3D block to be orthonormal without COMX; got overlap error $(gausslet_overlap_error)",
        ),
    )

    ngausslet = size(gausslet_overlap_value, 1)
    ngaussian = size(overlap_aa, 1)
    raw_overlap = [
        gausslet_overlap_value Matrix{Float64}(overlap_ga)
        Matrix{Float64}(transpose(overlap_ga)) Matrix{Float64}(overlap_aa)
    ]
    seed_projector = vcat(
        -(gausslet_overlap_value \ Matrix{Float64}(overlap_ga)),
        Matrix{Float64}(I, ngaussian, ngaussian),
    )
    residual_overlap = Matrix{Float64}(transpose(seed_projector) * raw_overlap * seed_projector)
    decomposition = eigen(Symmetric(residual_overlap))
    keep = findall(>(1.0e-10), decomposition.values)
    isempty(keep) && throw(
        ArgumentError("Qiu-White residual-Gaussian construction produced no nontrivial 3D residual directions"),
    )
    residual_coefficients =
        seed_projector *
        decomposition.vectors[:, keep] *
        Diagonal(1.0 ./ sqrt.(decomposition.values[keep]))
    gausslet_coefficients = vcat(
        Matrix{Float64}(I, ngausslet, ngausslet),
        zeros(Float64, ngaussian, ngausslet),
    )
    raw_to_final = hcat(gausslet_coefficients, residual_coefficients)
    final_overlap = Matrix{Float64}(transpose(raw_to_final) * raw_overlap * raw_to_final)
    return (
        raw_overlap = raw_overlap,
        raw_to_final = raw_to_final,
        residual_coefficients = residual_coefficients,
        final_overlap = final_overlap,
    )
end

# Alg QW-RG step 6: Construct exact raw-space one-body matrices and transform
# them into the final {G, R} basis.
# See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
function _qwrg_one_body_matrices(
    gausslet_one_body::AbstractMatrix{<:Real},
    one_body_ga::AbstractMatrix{<:Real},
    one_body_aa::AbstractMatrix{<:Real},
    raw_to_final::AbstractMatrix{<:Real},
)
    raw_one_body = [
        Matrix{Float64}(gausslet_one_body) Matrix{Float64}(one_body_ga)
        Matrix{Float64}(transpose(one_body_ga)) Matrix{Float64}(one_body_aa)
    ]
    final_one_body = Matrix{Float64}(transpose(raw_to_final) * raw_one_body * raw_to_final)
    return raw_one_body, 0.5 .* (final_one_body .+ transpose(final_one_body))
end

function _qwrg_residual_moment_data(
    raw_overlap::AbstractMatrix{<:Real},
    x_raw::AbstractMatrix{<:Real},
    x2_raw::AbstractMatrix{<:Real},
    y_raw::AbstractMatrix{<:Real},
    y2_raw::AbstractMatrix{<:Real},
    z_raw::AbstractMatrix{<:Real},
    z2_raw::AbstractMatrix{<:Real},
    raw_to_final::AbstractMatrix{<:Real},
    ngausslet::Int,
)
    residual_coefficients = Matrix{Float64}(raw_to_final[:, (ngausslet + 1):end])
    nresidual = size(residual_coefficients, 2)
    centers = zeros(Float64, nresidual, 3)
    widths = zeros(Float64, nresidual, 3)

    overlap_residual = Matrix{Float64}(transpose(residual_coefficients) * raw_overlap * residual_coefficients)
    for index in 1:nresidual
        vector = view(residual_coefficients, :, index)
        norm_value = Float64(dot(vector, raw_overlap * vector))
        norm_value > 1.0e-12 || throw(
            ArgumentError("residual moment extraction requires nonzero residual norm"),
        )

        x1 = Float64(dot(vector, x_raw * vector) / norm_value)
        y1 = Float64(dot(vector, y_raw * vector) / norm_value)
        z1 = Float64(dot(vector, z_raw * vector) / norm_value)
        x2 = Float64(dot(vector, x2_raw * vector) / norm_value)
        y2 = Float64(dot(vector, y2_raw * vector) / norm_value)
        z2 = Float64(dot(vector, z2_raw * vector) / norm_value)

        varx = x2 - x1^2
        vary = y2 - y1^2
        varz = z2 - z1^2
        min(varx, vary, varz) > 1.0e-12 || throw(
            ArgumentError("MWG residual moment extraction requires positive residual variances"),
        )

        centers[index, 1] = x1
        centers[index, 2] = y1
        centers[index, 3] = z1
        widths[index, 1] = sqrt(2.0 * varx)
        widths[index, 2] = sqrt(2.0 * vary)
        widths[index, 3] = sqrt(2.0 * varz)
    end

    overlap_error = norm(overlap_residual - I, Inf)
    return (
        centers = centers,
        widths = widths,
        overlap_error = overlap_error,
    )
end

function _qwrg_nearest_indices(
    gausslet_orbitals::AbstractVector{<:CartesianProductOrbital3D},
    residual_centers::AbstractMatrix{<:Real},
)
    return Int[
        argmin(
            [
                (orbital.x - residual_centers[index, 1])^2 +
                (orbital.y - residual_centers[index, 2])^2 +
                (orbital.z - residual_centers[index, 3])^2 for orbital in gausslet_orbitals
            ],
        ) for index in axes(residual_centers, 1)
    ]
end

function _qwrg_effective_gaussians(
    residual_centers::AbstractMatrix{<:Real},
    residual_widths::AbstractMatrix{<:Real},
)
    nresidual = size(residual_centers, 1)
    x_gaussians = Gaussian[]
    y_gaussians = Gaussian[]
    z_gaussians = Gaussian[]
    for index in 1:nresidual
        push!(
            x_gaussians,
            Gaussian(center = residual_centers[index, 1], width = residual_widths[index, 1]),
        )
        push!(
            y_gaussians,
            Gaussian(center = residual_centers[index, 2], width = residual_widths[index, 2]),
        )
        push!(
            z_gaussians,
            Gaussian(center = residual_centers[index, 3], width = residual_widths[index, 3]),
        )
    end
    return x_gaussians, y_gaussians, z_gaussians
end

function _qwrg_same_gaussians(
    left::AbstractVector{<:Gaussian},
    right::AbstractVector{<:Gaussian};
    tol::Real = 1.0e-12,
)
    length(left) == length(right) || return false
    for index in eachindex(left, right)
        abs(left[index].center_value - right[index].center_value) <= tol || return false
        abs(left[index].width - right[index].width) <= tol || return false
    end
    return true
end

function _qwrg_gaussian_analytic_blocks(
    gaussians::AbstractVector{<:Gaussian},
    expansion::CoulombGaussianExpansion,
)
    gaussian_set = PrimitiveSet1D(
        AbstractPrimitiveFunction1D[gaussian for gaussian in gaussians];
        name = :qiu_white_analytic_gaussians,
    )
    factors = Matrix{Float64}[
        Matrix{Float64}(gaussian_factor_matrix(
            gaussian_set;
            exponent = exponent,
            center = 0.0,
        )) for exponent in expansion.exponents
    ]
    pair_factors = _primitive_pair_gaussian_factor_matrices(
        gaussian_set,
        _select_primitive_matrix_backend(gaussian_set);
        exponents = expansion.exponents,
    )
    return (
        overlap = Matrix{Float64}(overlap_matrix(gaussian_set)),
        kinetic = Matrix{Float64}(kinetic_matrix(gaussian_set)),
        position = Matrix{Float64}(position_matrix(gaussian_set)),
        x2 = Matrix{Float64}(_x2_matrix(gaussian_set)),
        factors = factors,
        pair_factors = pair_factors,
    )
end

# Alg QW-RG step 8a: Build nearest-center / GGT residual-Gaussian interaction
# data in the same two-index IDA representation.
# See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
function _qwrg_interaction_matrix_nearest(
    gausslet_interaction::AbstractMatrix{<:Real},
    gausslet_orbitals::AbstractVector{<:CartesianProductOrbital3D},
    residual_centers::AbstractMatrix{<:Real},
)
    ngausslet = size(gausslet_interaction, 1)
    nresidual = size(residual_centers, 1)
    interaction = zeros(Float64, ngausslet + nresidual, ngausslet + nresidual)
    interaction[1:ngausslet, 1:ngausslet] .= Matrix{Float64}(gausslet_interaction)
    nearest = _qwrg_nearest_indices(gausslet_orbitals, residual_centers)
    for residual in 1:nresidual
        index = nearest[residual]
        interaction[1:ngausslet, ngausslet + residual] .= gausslet_interaction[:, index]
        interaction[ngausslet + residual, 1:ngausslet] .= transpose(gausslet_interaction[:, index])
    end
    for i in 1:nresidual, j in i:nresidual
        value = gausslet_interaction[nearest[i], nearest[j]]
        interaction[ngausslet + i, ngausslet + j] = value
        interaction[ngausslet + j, ngausslet + i] = value
    end
    return interaction
end

# Alg QW-RG step 8b, 8c, and 8d: Match exact residual moments by effective
# Gaussians and keep RG terms in the same two-index IDA representation.
# See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
function _qwrg_interaction_matrix_mwg(
    gausslet_bundle::_MappedOrdinaryGausslet1DBundle,
    gausslet_interaction::AbstractMatrix{<:Real},
    expansion::CoulombGaussianExpansion,
    residual_centers::AbstractMatrix{<:Real},
    residual_widths::AbstractMatrix{<:Real},
)
    ngausslet = size(gausslet_interaction, 1)
    nresidual = size(residual_centers, 1)
    interaction = zeros(Float64, ngausslet + nresidual, ngausslet + nresidual)
    interaction[1:ngausslet, 1:ngausslet] .= Matrix{Float64}(gausslet_interaction)

    x_gaussians, y_gaussians, z_gaussians = _qwrg_effective_gaussians(residual_centers, residual_widths)
    # Alg QW-RG step 8d: Evaluate RG-gausslet terms from exact contracted
    # 1D raw-space pair-factor blocks, then assemble the 3D IDA interaction.
    # See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
    pair_x = _qwrg_split_block_matrices(
        gausslet_bundle,
        x_gaussians,
        expansion,
    )
    pair_y = _qwrg_same_gaussians(x_gaussians, y_gaussians) ? pair_x :
        _qwrg_split_block_matrices(gausslet_bundle, y_gaussians, expansion)
    pair_z = if _qwrg_same_gaussians(x_gaussians, z_gaussians)
        pair_x
    elseif _qwrg_same_gaussians(y_gaussians, z_gaussians)
        pair_y
    else
        _qwrg_split_block_matrices(gausslet_bundle, z_gaussians, expansion)
    end
    analytic_x = _qwrg_gaussian_analytic_blocks(x_gaussians, expansion)
    analytic_y = _qwrg_same_gaussians(x_gaussians, y_gaussians) ? analytic_x : _qwrg_gaussian_analytic_blocks(y_gaussians, expansion)
    analytic_z = if _qwrg_same_gaussians(x_gaussians, z_gaussians)
        analytic_x
    elseif _qwrg_same_gaussians(y_gaussians, z_gaussians)
        analytic_y
    else
        _qwrg_gaussian_analytic_blocks(z_gaussians, expansion)
    end

    for residual in 1:nresidual
        column = zeros(Float64, ngausslet)
        scratch = zeros(Float64, ngausslet)
        for term in eachindex(expansion.coefficients)
            fx = view(pair_x.pair_ga[term], :, residual)
            fy = view(pair_y.pair_ga[term], :, residual)
            fz = view(pair_z.pair_ga[term], :, residual)
            _qwrg_fill_product_column!(scratch, fx, fy, fz)
            column .+= expansion.coefficients[term] .* scratch
        end
        interaction[1:ngausslet, ngausslet + residual] .= column
        interaction[ngausslet + residual, 1:ngausslet] .= column
    end

    for i in 1:nresidual, j in i:nresidual
        value = 0.0
        for term in eachindex(expansion.coefficients)
            value += expansion.coefficients[term] *
                analytic_x.pair_factors[term][i, j] *
                analytic_y.pair_factors[term][i, j] *
                analytic_z.pair_factors[term][i, j]
        end
        interaction[ngausslet + i, ngausslet + j] = value
        interaction[ngausslet + j, ngausslet + i] = value
    end

    return interaction
end

function _qwrg_orbital_data(
    gausslet_orbitals::AbstractVector{<:CartesianProductOrbital3D},
    residual_centers::AbstractMatrix{<:Real},
    residual_widths::AbstractMatrix{<:Real},
)
    orbitals_out = QiuWhiteHybridOrbital3D[]
    for orbital in gausslet_orbitals
        push!(
            orbitals_out,
            QiuWhiteHybridOrbital3D(
                orbital.index,
                :gausslet,
                "g($(orbital.ix),$(orbital.iy),$(orbital.iz))",
                orbital.x,
                orbital.y,
                orbital.z,
                NaN,
                NaN,
                NaN,
            ),
        )
    end
    base_index = length(gausslet_orbitals)
    for index in axes(residual_centers, 1)
        push!(
            orbitals_out,
            QiuWhiteHybridOrbital3D(
                base_index + index,
                :residual_gaussian,
                "rg$index",
                residual_centers[index, 1],
                residual_centers[index, 2],
                residual_centers[index, 3],
                residual_widths[index, 1],
                residual_widths[index, 2],
                residual_widths[index, 3],
            ),
        )
    end
    return orbitals_out
end

"""
    ordinary_cartesian_qiu_white_operators(
        basis::MappedUniformBasis,
        gaussian_data::LegacySGaussianData;
        expansion = coulomb_gaussian_expansion(doacc = false),
        Z = 2.0,
        interaction_treatment = :mwg,
        gausslet_backend = :numerical_reference,
        timing = false,
    )

Build the paper-faithful Qiu-White residual-Gaussian ordinary Cartesian
reference Hamiltonian.

This path is intentionally separate from the current COMX/localized hybrid
route. It:

- keeps the full 3D gausslet product basis fixed
- orthogonalizes the added 3D Gaussian supplement against that full 3D space
- builds the one-body matrices exactly in the raw gausslet-plus-GTO space
- keeps the two-electron interaction in the same two-index integral-diagonal
  approximation (IDA) representation used for the gausslet channel

Allowed `interaction_treatment` values are:

- `:ggt_nearest`
- `:mwg`

Set `timing = true` to print a coarse constructor-only phase breakdown for
debugging the reference implementation. This is intentionally narrow and does
not add timing noise to the broader library.

This is a reference path for validating the Qiu-White formulation, not yet a
claim that the ordinary branch is solver-ready.
"""
function ordinary_cartesian_qiu_white_operators(
    basis::MappedUniformBasis,
    gaussian_data::LegacySGaussianData;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
    interaction_treatment::Symbol = :mwg,
    gausslet_backend::Symbol = :numerical_reference,
    timing::Bool = false,
    )
    gausslet_backend == :numerical_reference || throw(
        ArgumentError("ordinary_cartesian_qiu_white_operators currently supports only gausslet_backend = :numerical_reference"),
    )
    timings = Pair{String,Float64}[]
    timing_io = stdout

    start_ns = time_ns()
    gausslet_bundle = _mapped_ordinary_gausslet_1d_bundle(
        basis,
        exponents = expansion.exponents,
        center = 0.0,
        backend = gausslet_backend,
    )
    timing && _qwrg_record_timing!(timing_io, timings, "shared gausslet-side 1D bundle", start_ns)

    start_ns = time_ns()
    blocks = _qwrg_split_block_matrices(
        gausslet_bundle,
        gaussian_data.gaussians,
        expansion,
    )
    timing && _qwrg_record_timing!(timing_io, timings, "split gausslet-Gaussian raw-block assembly", start_ns)

    gausslet_orbitals = _mapped_cartesian_orbitals(gausslet_bundle.centers)
    gausslet_count = length(gausslet_orbitals)
    start_ns = time_ns()
    gausslet_overlap_3d = zeros(Float64, length(gausslet_orbitals), length(gausslet_orbitals))
    _qwrg_fill_product_matrix!(
        gausslet_overlap_3d,
        blocks.overlap_gg,
        blocks.overlap_gg,
        blocks.overlap_gg,
    )
    timing && _qwrg_record_timing!(timing_io, timings, "3D gausslet overlap assembly", start_ns)

    start_ns = time_ns()
    overlap_ga, overlap_aa = _qwrg_raw_overlap_blocks((
        overlap_gg = blocks.overlap_gg,
        overlap_ga = blocks.overlap_ga,
        overlap_aa = blocks.overlap_aa,
    ))
    residual_data = _qwrg_residual_space(gausslet_overlap_3d, overlap_ga, overlap_aa)
    timing && _qwrg_record_timing!(timing_io, timings, "residual-space construction", start_ns)

    start_ns = time_ns()
    one_body_ga, one_body_aa = _qwrg_raw_one_body_blocks((
        overlap_ga = blocks.overlap_ga,
        overlap_aa = blocks.overlap_aa,
        kinetic_ga = blocks.kinetic_ga,
        kinetic_aa = blocks.kinetic_aa,
        factor_ga = blocks.factor_ga,
        factor_aa = blocks.factor_aa,
    ), expansion; Z = Z)
    gausslet_one_body = _qwrg_gausslet_one_body_matrix(blocks, expansion; Z = Z)
    timing && _qwrg_record_timing!(timing_io, timings, "3D gausslet one-body assembly", start_ns)

    start_ns = time_ns()
    _, final_one_body = _qwrg_one_body_matrices(
        gausslet_one_body,
        one_body_ga,
        one_body_aa,
        residual_data.raw_to_final,
    )
    timing && _qwrg_record_timing!(timing_io, timings, "raw one-body transform", start_ns)

    start_ns = time_ns()
    x_gg, x_ga, x_aa = _qwrg_raw_axis_blocks((
        overlap_gg = blocks.overlap_gg,
        overlap_ga = blocks.overlap_ga,
        overlap_aa = blocks.overlap_aa,
        position_gg = blocks.position_gg,
        position_ga = blocks.position_ga,
        position_aa = blocks.position_aa,
        x2_gg = blocks.x2_gg,
        x2_ga = blocks.x2_ga,
        x2_aa = blocks.x2_aa,
    ), :x)
    y_gg, y_ga, y_aa = _qwrg_raw_axis_blocks((
        overlap_gg = blocks.overlap_gg,
        overlap_ga = blocks.overlap_ga,
        overlap_aa = blocks.overlap_aa,
        position_gg = blocks.position_gg,
        position_ga = blocks.position_ga,
        position_aa = blocks.position_aa,
        x2_gg = blocks.x2_gg,
        x2_ga = blocks.x2_ga,
        x2_aa = blocks.x2_aa,
    ), :y)
    z_gg, z_ga, z_aa = _qwrg_raw_axis_blocks((
        overlap_gg = blocks.overlap_gg,
        overlap_ga = blocks.overlap_ga,
        overlap_aa = blocks.overlap_aa,
        position_gg = blocks.position_gg,
        position_ga = blocks.position_ga,
        position_aa = blocks.position_aa,
        x2_gg = blocks.x2_gg,
        x2_ga = blocks.x2_ga,
        x2_aa = blocks.x2_aa,
    ), :z)
    x2_gg, x2_ga, x2_aa = _qwrg_raw_axis_blocks((
        overlap_gg = blocks.overlap_gg,
        overlap_ga = blocks.overlap_ga,
        overlap_aa = blocks.overlap_aa,
        position_gg = blocks.position_gg,
        position_ga = blocks.position_ga,
        position_aa = blocks.position_aa,
        x2_gg = blocks.x2_gg,
        x2_ga = blocks.x2_ga,
        x2_aa = blocks.x2_aa,
    ), :x; squared = true)
    y2_gg, y2_ga, y2_aa = _qwrg_raw_axis_blocks((
        overlap_gg = blocks.overlap_gg,
        overlap_ga = blocks.overlap_ga,
        overlap_aa = blocks.overlap_aa,
        position_gg = blocks.position_gg,
        position_ga = blocks.position_ga,
        position_aa = blocks.position_aa,
        x2_gg = blocks.x2_gg,
        x2_ga = blocks.x2_ga,
        x2_aa = blocks.x2_aa,
    ), :y; squared = true)
    z2_gg, z2_ga, z2_aa = _qwrg_raw_axis_blocks((
        overlap_gg = blocks.overlap_gg,
        overlap_ga = blocks.overlap_ga,
        overlap_aa = blocks.overlap_aa,
        position_gg = blocks.position_gg,
        position_ga = blocks.position_ga,
        position_aa = blocks.position_aa,
        x2_gg = blocks.x2_gg,
        x2_ga = blocks.x2_ga,
        x2_aa = blocks.x2_aa,
    ), :z; squared = true)

    x_raw = [x_gg x_ga; transpose(x_ga) x_aa]
    y_raw = [y_gg y_ga; transpose(y_ga) y_aa]
    z_raw = [z_gg z_ga; transpose(z_ga) z_aa]
    x2_raw = [x2_gg x2_ga; transpose(x2_ga) x2_aa]
    y2_raw = [y2_gg y2_ga; transpose(y2_ga) y2_aa]
    z2_raw = [z2_gg z2_ga; transpose(z2_ga) z2_aa]

    moment_data = _qwrg_residual_moment_data(
        residual_data.raw_overlap,
        x_raw,
        x2_raw,
        y_raw,
        y2_raw,
        z_raw,
        z2_raw,
        residual_data.raw_to_final,
        gausslet_count,
    )
    timing && _qwrg_record_timing!(timing_io, timings, "raw moment-matrix assembly", start_ns)

    start_ns = time_ns()
    gausslet_interaction = _qwrg_gausslet_interaction_matrix(blocks, expansion)
    timing && _qwrg_record_timing!(timing_io, timings, "3D gausslet interaction assembly", start_ns)

    start_ns = time_ns()
    interaction_matrix = if interaction_treatment == :ggt_nearest
        _qwrg_interaction_matrix_nearest(
            gausslet_interaction,
            gausslet_orbitals,
            moment_data.centers,
        )
    elseif interaction_treatment == :mwg
        _qwrg_interaction_matrix_mwg(
            gausslet_bundle,
            gausslet_interaction,
            expansion,
            moment_data.centers,
            moment_data.widths,
        )
    else
        throw(ArgumentError("Qiu-White interaction_treatment must be :ggt_nearest or :mwg"))
    end
    timing && _qwrg_record_timing!(timing_io, timings, "RG interaction assembly", start_ns)
    timing && _qwrg_maybe_print_timings(timing_io, timings)

    return QiuWhiteResidualGaussianOperators(
        basis,
        gaussian_data,
        gausslet_backend,
        interaction_treatment,
        expansion,
        Matrix{Float64}(residual_data.final_overlap),
        final_one_body,
        Matrix{Float64}(0.5 .* (interaction_matrix .+ transpose(interaction_matrix))),
        _qwrg_orbital_data(
            gausslet_orbitals,
            moment_data.centers,
            moment_data.widths,
        ),
        gausslet_count,
        size(moment_data.centers, 1),
        Matrix{Float64}(residual_data.raw_to_final),
        Matrix{Float64}(moment_data.centers),
        Matrix{Float64}(moment_data.widths),
    )
end
