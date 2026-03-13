const _BASIS_GRID_SIGMA = 3.0
const _BASIS_GRID_S0 = 6.5
const _BASIS_GRID_MARGIN = 20.0
const _BASIS_EIG_TOL = 1.0e-10
const _BASIS_STENCIL_TOL = 1.0e-12

_as_family(family_value::GaussletFamily) = family_value
_as_family(family_value::Symbol) = GaussletFamily(family_value)

function _identity_mapping(mapping_value::AbstractCoordinateMapping)
    return mapping_value isa IdentityMapping
end

function _basis_grid_h(reference_spacing::Float64)
    return min(0.02, reference_spacing / 25.0)
end

function _make_erf_grid(; h::Real, rmax::Real, sigma::Real = _BASIS_GRID_SIGMA, s0::Real = _BASIS_GRID_S0)
    h > 0 || throw(ArgumentError("internal grid spacing must be positive"))
    rmax > 0 || throw(ArgumentError("internal grid extent must be positive"))

    points = Float64[]
    weights = Float64[]
    k = 0
    root_pi = sqrt(pi)
    while true
        s = k + 0.5
        z = (s / sigma) - s0
        gate = 0.5 * _erfc(-z)
        gate_prime = exp(-(z * z)) / (root_pi * sigma)
        point = h * s * gate
        weight = h * (gate + s * gate_prime)
        push!(points, point)
        push!(weights, weight)
        point >= rmax && break
        k += 1
    end
    return points, weights
end

_integrate(values::AbstractVector{<:Real}, weights::AbstractVector{<:Real}) = sum(weights .* values)
_integrate_x(values::AbstractVector{<:Real}, x::AbstractVector{<:Real}, weights::AbstractVector{<:Real}) =
    sum(weights .* x .* values)

function _seed_scalar_integrals(phi::Function, js::Vector{Int}, x::Vector{Float64}, weights::Vector{Float64})
    nseed = length(js)
    norms = zeros(Float64, nseed)
    totals = zeros(Float64, nseed)
    sampled = zeros(Float64, length(x), nseed)
    for column in 1:nseed
        values = [phi(js[column], xi) for xi in x]
        sampled[:, column] = values
        norms[column] = _integrate(values .* values, weights)
        totals[column] = _integrate(values, weights)
    end
    return norms, totals, sampled
end

function _assemble_halfline_seed(
    phi::Function,
    js::Vector{Int},
    x::Vector{Float64},
    weights::Vector{Float64};
    scale::Union{Nothing, Vector{Float64}} = nothing,
)
    nseed = length(js)
    overlap = zeros(Float64, nseed, nseed)
    position = zeros(Float64, nseed, nseed)
    for a in 1:nseed
        va = [phi(js[a], xi) for xi in x]
        if scale !== nothing
            va .*= scale[a]
        end
        for b in a:nseed
            vb = [phi(js[b], xi) for xi in x]
            if scale !== nothing
                vb .*= scale[b]
            end
            product = va .* vb
            s_ab = _integrate(product, weights)
            x_ab = _integrate_x(product, x, weights)
            overlap[a, b] = s_ab
            overlap[b, a] = s_ab
            position[a, b] = x_ab
            position[b, a] = x_ab
        end
    end
    return overlap, position
end

function _pivot_projection_transform(seed::Gausslet, js::Vector{Int}, scale::Vector{Float64}, spacing::Float64)
    nseed = length(js)
    values_at_origin = zeros(Float64, nseed)
    for index in 1:nseed
        values_at_origin[index] = scale[index] * seed(-js[index] * spacing)
    end
    pivot = argmax(abs.(values_at_origin))
    pivot_value = values_at_origin[pivot]
    pivot_value == 0.0 && throw(ArgumentError("failed to identify a delta-at-origin pivot direction"))

    transform = zeros(Float64, nseed, nseed - 1)
    column = 1
    for row in 1:nseed
        row == pivot && continue
        transform[row, column] = 1.0
        transform[pivot, column] -= values_at_origin[row] / pivot_value
        column += 1
    end
    return transform
end

function _s_invsqrt_reduced(overlap::AbstractMatrix{<:Real}; tol::Float64 = _BASIS_EIG_TOL)
    decomposition = eigen(Symmetric(Matrix{Float64}(overlap)))
    keep = findall(>(tol), decomposition.values)
    isempty(keep) && throw(ArgumentError("no overlap eigenvalues exceeded the internal cutoff"))
    vectors = decomposition.vectors[:, keep]
    invhalf = Diagonal(1.0 ./ sqrt.(decomposition.values[keep]))
    return vectors, invhalf
end

function _comx_reduced(vectors::AbstractMatrix{<:Real}, invhalf::Diagonal, position::AbstractMatrix{<:Real})
    projected = Matrix(Symmetric(vectors' * position * vectors))
    localized = Matrix(Symmetric(invhalf * projected * invhalf))
    decomposition = eigen(Symmetric(localized))
    return decomposition.vectors, decomposition.values
end

function _sign_fix_columns!(coefficients_matrix::Matrix{Float64}, sampled_values::Matrix{Float64}, weights::Vector{Float64})
    integrals = vec((weights' * sampled_values) * coefficients_matrix)
    for column in 1:size(coefficients_matrix, 2)
        if integrals[column] < 0.0
            coefficients_matrix[:, column] .*= -1.0
        end
    end
    return coefficients_matrix
end

function _finalize_localized_basis(
    base_coefficients::Matrix{Float64},
    sampled_values::Matrix{Float64},
    position_seed::Matrix{Float64},
    weights::Vector{Float64},
)
    weighted_values = sampled_values * base_coefficients
    weighted_values .*= sqrt.(weights)
    factorization = qr(weighted_values)
    ortho_coefficients = base_coefficients / Matrix(factorization.R)
    localized_position = Matrix(Symmetric(ortho_coefficients' * position_seed * ortho_coefficients))
    decomposition = eigen(Symmetric(localized_position))
    final_coefficients = ortho_coefficients * decomposition.vectors
    _sign_fix_columns!(final_coefficients, sampled_values, weights)
    return final_coefficients, decomposition.values
end

function _halfline_seed_primitive_layer(js::Vector{Int}, family_value::GaussletFamily, spacing::Float64)
    full_coefficients = _full_family_coefficients(family_value.positive_coefficients)
    radius = length(family_value.positive_coefficients) - 1
    primitive_indices = collect((3 * first(js) - radius):(3 * last(js) + radius))
    matrix = zeros(Float64, length(primitive_indices), length(js))
    scale = 1.0 / sqrt(2.0 * pi * spacing)
    stencil_coefficients = Float64[scale * coefficient for coefficient in full_coefficients]

    for column in eachindex(js)
        j = js[column]
        row = 1
        for offset in primitive_indices
            local_offset = offset - 3 * j
            if -radius <= local_offset <= radius
                matrix[row, column] = stencil_coefficients[local_offset + radius + 1]
            end
            row += 1
        end
    end

    primitive_spacing = spacing / 3.0
    primitive_list = AbstractPrimitiveFunction1D[
        HalfLineGaussian(center = index * primitive_spacing, width = primitive_spacing)
        for index in primitive_indices
    ]
    return matrix, primitive_list
end

function _with_mapping(primitives_ref::Vector{AbstractPrimitiveFunction1D}, mapping_value::AbstractCoordinateMapping)
    if _identity_mapping(mapping_value)
        return primitives_ref
    end
    return AbstractPrimitiveFunction1D[Distorted(primitive, mapping_value) for primitive in primitives_ref]
end

function _stencil_from_column(
    coefficients_column::AbstractVector{<:Real},
    primitive_list::Vector{AbstractPrimitiveFunction1D},
)
    keep = findall(coefficient -> abs(coefficient) > _BASIS_STENCIL_TOL, coefficients_column)
    if isempty(keep)
        keep = [argmax(abs.(coefficients_column))]
    end
    return FunctionStencil(coefficients_column[keep], primitive_list[keep])
end

function _integral_weights_from_stencils(stencils::Vector{FunctionStencil})
    weights = Float64[]
    for stencil_value in stencils
        total = 0.0
        for term in terms(stencil_value)
            total += term.coefficient * integral_weight(term.primitive)
        end
        push!(weights, total)
    end
    return weights
end

function _build_halfline_coefficients(spec)
    spacing = spec.reference_spacing
    reference_max = uofx(spec.mapping_value, spec.xmax)
    reference_max > 0.0 || throw(ArgumentError("HalfLineBasisSpec requires xmax to map to a positive reference range"))

    jneg = spec.tails
    jpos_save = max(1, ceil(Int, reference_max / spacing))
    jpos_build = jpos_save + spec.tails

    seed = Gausslet(spec.family_value; center = 0.0, spacing = spacing)
    phi = (j::Int, x::Float64) -> seed(x - j * spacing)
    grid_h = _basis_grid_h(spacing)
    grid_max = max(reference_max + _BASIS_GRID_MARGIN * spacing, (jpos_build + _BASIS_GRID_MARGIN) * spacing)
    xgrid, weights = _make_erf_grid(; h = grid_h, rmax = grid_max)
    js = collect(-jneg:jpos_build)

    norms, _, sampled_values = _seed_scalar_integrals(phi, js, xgrid, weights)
    scale = 1.0 ./ sqrt.(norms)
    projection = _pivot_projection_transform(seed, js, scale, spacing)
    overlap_seed, position_seed = _assemble_halfline_seed(phi, js, xgrid, weights; scale = scale)

    reduced_overlap = Matrix(Symmetric(projection' * overlap_seed * projection))
    reduced_position = Matrix(Symmetric(projection' * position_seed * projection))
    vectors, invhalf = _s_invsqrt_reduced(reduced_overlap)
    localizer, reference_centers_full = _comx_reduced(vectors, invhalf, reduced_position)
    coefficients_full = Diagonal(scale) * (projection * (vectors * (invhalf * localizer)))
    _sign_fix_columns!(coefficients_full, sampled_values, weights)

    nkeep = jneg + jpos_save
    return js, coefficients_full[:, 1:nkeep], reference_centers_full[1:nkeep]
end

function _xgaussian_sample_matrix(xgaussians::Vector{XGaussian}, xgrid::Vector{Float64})
    matrix = zeros(Float64, length(xgrid), length(xgaussians))
    for column in eachindex(xgaussians)
        matrix[:, column] = [value(xgaussians[column], x) for x in xgrid]
    end
    return matrix
end

function _build_radial_coefficients(spec)
    spacing = spec.reference_spacing
    ninj = length(spec.xgaussians)
    target_odd_count =
        if spec.count !== nothing
            spec.count > spec.odd_even_kmax + ninj ||
                throw(ArgumentError("RadialBasisSpec count must exceed odd_even_kmax + length(xgaussians)"))
            spec.count - spec.odd_even_kmax - ninj
        else
            reference_max = uofx(spec.mapping_value, spec.rmax)
            reference_max > 0.0 || throw(ArgumentError("RadialBasisSpec requires rmax to map to a positive reference range"))
            max(1, ceil(Int, reference_max / spacing))
        end

    lseed = target_odd_count + spec.tails
    lseed >= 1 || throw(ArgumentError("radial construction requires at least one odd mode"))
    spec.odd_even_kmax <= lseed ||
        throw(ArgumentError("odd_even_kmax must not exceed the internal seed-window size"))

    seed = Gausslet(spec.family_value; center = 0.0, spacing = spacing)
    phi = (j::Int, x::Float64) -> seed(x - j * spacing)
    grid_h = _basis_grid_h(spacing)
    reference_target =
        spec.count === nothing ? uofx(spec.mapping_value, spec.rmax) : target_odd_count * spacing
    grid_max = max(reference_target + _BASIS_GRID_MARGIN * spacing, (lseed + _BASIS_GRID_MARGIN) * spacing)
    xgrid, weights = _make_erf_grid(; h = grid_h, rmax = grid_max)
    js = collect(-lseed:lseed)

    _, seed_integrals, sampled_seed = _seed_scalar_integrals(phi, js, xgrid, weights)
    overlap_seed = sampled_seed' * (weights .* sampled_seed)
    position_seed = sampled_seed' * ((xgrid .* weights) .* sampled_seed)

    nseed = length(js)
    odd_block = zeros(Float64, nseed, lseed)
    for k in 1:lseed
        ip = (lseed + 1) + k
        im = (lseed + 1) - k
        odd_block[ip, k] = 1.0
        odd_block[im, k] = -1.0
    end

    even_block = zeros(Float64, nseed, spec.odd_even_kmax + 1)
    even_block[lseed + 1, 1] = 1.0
    for k in 1:spec.odd_even_kmax
        ip = (lseed + 1) + k
        im = (lseed + 1) - k
        even_block[ip, k + 1] = 1.0
        even_block[im, k + 1] = 1.0
    end

    function normalize_columns(block::Matrix{Float64})
        overlap = overlap_seed * block
        norms = vec(sum(block .* overlap; dims = 1))
        any(norms .<= 0.0) && throw(ArgumentError("encountered a nonpositive block norm during radial construction"))
        return block * Diagonal(1.0 ./ sqrt.(norms))
    end

    odd_normalized = normalize_columns(odd_block)
    odd_overlap = Matrix(Symmetric(odd_normalized' * overlap_seed * odd_normalized))
    odd_position = Matrix(Symmetric(odd_normalized' * position_seed * odd_normalized))
    odd_vectors, odd_invhalf = _s_invsqrt_reduced(odd_overlap)
    odd_localizer, _ = _comx_reduced(odd_vectors, odd_invhalf, odd_position)
    odd_coefficients = odd_normalized * (odd_vectors * (odd_invhalf * odd_localizer))

    even_clean = Matrix{Float64}(undef, nseed, 0)
    if spec.odd_even_kmax > 0
        even_normalized = normalize_columns(even_block)
        even_projection = odd_coefficients' * (overlap_seed * even_normalized)
        even_clean = even_normalized .- (odd_coefficients * even_projection)

        delta_values = [seed(-j * spacing) for j in js]
        delta_direction = vec(delta_values' * even_clean)
        denom = dot(delta_direction, delta_direction)
        if denom > 0.0
            unit = delta_direction ./ sqrt(denom)
            sign_value = unit[1] >= 0.0 ? 1.0 : -1.0
            reflector = copy(unit)
            reflector[1] += sign_value
            beta = dot(reflector, reflector)
            beta == 0.0 && throw(ArgumentError("failed to remove the radial delta direction"))
            householder = Matrix{Float64}(I, length(unit), length(unit)) .-
                          (2.0 / beta) .* (reflector * reflector')
            even_clean = even_clean * householder[:, 2:end]
        end
    end

    combined_base = hcat(odd_coefficients, even_clean)
    sampled_combined = sampled_seed
    position_combined = position_seed
    integral_combined = seed_integrals

    if ninj > 0
        sampled_x = _xgaussian_sample_matrix(spec.xgaussians, xgrid)
        position_sx = sampled_seed' * ((xgrid .* weights) .* sampled_x)
        position_xx = sampled_x' * ((xgrid .* weights) .* sampled_x)
        position_combined = [position_seed position_sx; position_sx' position_xx]
        sampled_combined = hcat(sampled_seed, sampled_x)
        integral_combined = vcat(seed_integrals, vec(sampled_x' * weights))
        combined_base = hcat(
            vcat(combined_base, zeros(Float64, ninj, size(combined_base, 2))),
            vcat(zeros(Float64, nseed, ninj), Matrix{Float64}(I, ninj, ninj)),
        )
    end

    final_coefficients, reference_centers_full = _finalize_localized_basis(
        combined_base,
        sampled_combined,
        position_combined,
        weights,
    )

    target_count =
        spec.count === nothing ?
        target_odd_count + spec.odd_even_kmax + ninj :
        spec.count
    return js, final_coefficients[:, 1:target_count], reference_centers_full[1:target_count]
end

function _boundary_stencils(
    js::Vector{Int},
    coefficients_seed::Matrix{Float64},
    family_value::GaussletFamily,
    spacing::Float64,
    mapping_value::AbstractCoordinateMapping,
)
    primitive_matrix, primitive_ref = _halfline_seed_primitive_layer(js, family_value, spacing)
    public_primitives = _with_mapping(primitive_ref, mapping_value)
    primitive_coefficients = primitive_matrix * coefficients_seed
    return [_stencil_from_column(primitive_coefficients[:, column], public_primitives)
            for column in 1:size(primitive_coefficients, 2)]
end

function _radial_stencils(
    js::Vector{Int},
    coefficients_full::Matrix{Float64},
    family_value::GaussletFamily,
    spacing::Float64,
    xgaussians::Vector{XGaussian},
    mapping_value::AbstractCoordinateMapping,
)
    nseed = length(js)
    ninj = length(xgaussians)
    primitive_matrix, primitive_ref_halfline = _halfline_seed_primitive_layer(js, family_value, spacing)
    seed_coefficients = coefficients_full[1:nseed, :]
    radial_coefficients = primitive_matrix * seed_coefficients
    primitive_ref = AbstractPrimitiveFunction1D[primitive_ref_halfline...]

    if ninj > 0
        inj_coefficients = coefficients_full[(nseed + 1):end, :]
        radial_coefficients = vcat(radial_coefficients, inj_coefficients)
        append!(primitive_ref, AbstractPrimitiveFunction1D[x for x in xgaussians])
    end

    public_primitives = _with_mapping(primitive_ref, mapping_value)
    return [_stencil_from_column(radial_coefficients[:, column], public_primitives)
            for column in 1:size(radial_coefficients, 2)]
end

function _physical_centers(reference_center_data::Vector{Float64}, mapping_value::AbstractCoordinateMapping)
    if _identity_mapping(mapping_value)
        return copy(reference_center_data)
    end
    return Float64[xofu(mapping_value, value) for value in reference_center_data]
end

"""
    UniformBasisSpec(family; xmin, xmax, spacing=1.0)

Recipe for building a uniform full-line gausslet basis over `[xmin, xmax]`.
"""
struct UniformBasisSpec <: AbstractBasisSpec
    family_value::GaussletFamily
    xmin::Float64
    xmax::Float64
    spacing::Float64

    function UniformBasisSpec(family_value::GaussletFamily, xmin::Real, xmax::Real, spacing::Real)
        xmin_value = Float64(xmin)
        xmax_value = Float64(xmax)
        spacing_value = Float64(spacing)
        xmax_value >= xmin_value || throw(ArgumentError("UniformBasisSpec requires xmax >= xmin"))
        spacing_value > 0.0 || throw(ArgumentError("UniformBasisSpec requires spacing > 0"))
        new(family_value, xmin_value, xmax_value, spacing_value)
    end
end

UniformBasisSpec(family_value::Union{GaussletFamily, Symbol}; xmin::Real, xmax::Real, spacing::Real = 1.0) =
    UniformBasisSpec(_as_family(family_value), xmin, xmax, spacing)

"""
    HalfLineBasisSpec(family; xmax, reference_spacing=1.0, tails=6,
                      mapping=IdentityMapping())

Recipe for building a boundary-correct half-line gausslet basis on `[0, xmax]`
or, when mapped, over approximately the corresponding physical range.
"""
struct HalfLineBasisSpec{M <: AbstractCoordinateMapping} <: AbstractBasisSpec
    family_value::GaussletFamily
    xmax::Float64
    reference_spacing::Float64
    tails::Int
    mapping_value::M

    function HalfLineBasisSpec(
        family_value::GaussletFamily,
        xmax::Real,
        reference_spacing::Real,
        tails::Int,
        mapping_value::M,
    ) where {M <: AbstractCoordinateMapping}
        xmax_value = Float64(xmax)
        spacing_value = Float64(reference_spacing)
        xmax_value > 0.0 || throw(ArgumentError("HalfLineBasisSpec requires xmax > 0"))
        spacing_value > 0.0 || throw(ArgumentError("HalfLineBasisSpec requires reference_spacing > 0"))
        tails >= 0 || throw(ArgumentError("HalfLineBasisSpec requires tails >= 0"))
        new{M}(family_value, xmax_value, spacing_value, tails, mapping_value)
    end
end

HalfLineBasisSpec(
    family_value::Union{GaussletFamily, Symbol};
    xmax::Real,
    reference_spacing::Real = 1.0,
    tails::Int = 6,
    mapping::AbstractCoordinateMapping = IdentityMapping(),
) = HalfLineBasisSpec(_as_family(family_value), xmax, reference_spacing, tails, mapping)

"""
    RadialBasisSpec(family; rmax, mapping, reference_spacing=1.0, tails=6,
                    odd_even_kmax=6, xgaussians=XGaussian[])
    RadialBasisSpec(family; count, mapping, reference_spacing=1.0, tails=6,
                    odd_even_kmax=6, xgaussians=XGaussian[])

Recipe for building a radial gausslet basis for the reduced radial function
`u(r) = r R(r)`.
"""
struct RadialBasisSpec{M <: AbstractCoordinateMapping} <: AbstractBasisSpec
    family_value::GaussletFamily
    rmax::Union{Nothing, Float64}
    count::Union{Nothing, Int}
    mapping_value::M
    reference_spacing::Float64
    tails::Int
    odd_even_kmax::Int
    xgaussians::Vector{XGaussian}

    function RadialBasisSpec(
        family_value::GaussletFamily,
        rmax::Union{Nothing, Real},
        count::Union{Nothing, Int},
        mapping_value::M,
        reference_spacing::Real,
        tails::Int,
        odd_even_kmax::Int,
        xgaussians::AbstractVector{XGaussian},
    ) where {M <: AbstractCoordinateMapping}
        ((rmax === nothing) == (count === nothing)) &&
            throw(ArgumentError("provide exactly one of rmax or count"))
        spacing_value = Float64(reference_spacing)
        spacing_value > 0.0 || throw(ArgumentError("RadialBasisSpec requires reference_spacing > 0"))
        tails >= 0 || throw(ArgumentError("RadialBasisSpec requires tails >= 0"))
        odd_even_kmax >= 0 || throw(ArgumentError("RadialBasisSpec requires odd_even_kmax >= 0"))
        rmax_value = rmax === nothing ? nothing : Float64(rmax)
        rmax_value === nothing || rmax_value > 0.0 || throw(ArgumentError("RadialBasisSpec requires rmax > 0"))
        count === nothing || count > 0 || throw(ArgumentError("RadialBasisSpec requires count > 0"))
        new{M}(family_value, rmax_value, count, mapping_value, spacing_value, tails, odd_even_kmax, collect(xgaussians))
    end
end

RadialBasisSpec(
    family_value::Union{GaussletFamily, Symbol};
    rmax::Union{Nothing, Real} = nothing,
    count::Union{Nothing, Int} = nothing,
    mapping::AbstractCoordinateMapping,
    reference_spacing::Real = 1.0,
    tails::Int = 6,
    odd_even_kmax::Int = 6,
    xgaussians::AbstractVector{XGaussian} = XGaussian[],
) = RadialBasisSpec(_as_family(family_value), rmax, count, mapping, reference_spacing, tails, odd_even_kmax, xgaussians)

"""
    UniformBasis

Concrete uniform full-line gausslet basis.
"""
struct UniformBasis
    spec::UniformBasisSpec
    center_data::Vector{Float64}
    integral_weight_data::Vector{Float64}
end

"""
    HalfLineBasis

Concrete boundary-correct half-line gausslet basis.
"""
struct HalfLineBasis{M <: AbstractCoordinateMapping}
    spec::HalfLineBasisSpec{M}
    stencil_data::Vector{FunctionStencil}
    reference_center_data::Vector{Float64}
    center_data::Vector{Float64}
    integral_weight_data::Vector{Float64}
end

"""
    RadialBasis

Concrete radial gausslet basis for the reduced radial function `u(r) = r R(r)`.
"""
struct RadialBasis{M <: AbstractCoordinateMapping}
    spec::RadialBasisSpec{M}
    stencil_data::Vector{FunctionStencil}
    reference_center_data::Vector{Float64}
    center_data::Vector{Float64}
    integral_weight_data::Vector{Float64}
end

"""
    BoundaryGausslet

Lightweight callable basis-function view returned by `HalfLineBasis[i]`.
"""
struct BoundaryGausslet{B <: HalfLineBasis} <: AbstractBasisFunction1D
    basis::B
    index::Int
end

"""
    RadialGausslet

Lightweight callable basis-function view returned by `RadialBasis[i]`.
"""
struct RadialGausslet{B <: RadialBasis} <: AbstractBasisFunction1D
    basis::B
    index::Int
end

Base.length(basis::UniformBasis) = length(basis.center_data)
Base.length(basis::HalfLineBasis) = length(basis.stencil_data)
Base.length(basis::RadialBasis) = length(basis.stencil_data)

function Base.getindex(basis::UniformBasis, index::Integer)
    return Gausslet(basis.spec.family_value; center = basis.center_data[index], spacing = basis.spec.spacing)
end

Base.getindex(basis::HalfLineBasis, index::Integer) = BoundaryGausslet(basis, index)
Base.getindex(basis::RadialBasis, index::Integer) = RadialGausslet(basis, index)

basis_spec(basis::UniformBasis) = basis.spec
basis_spec(basis::HalfLineBasis) = basis.spec
basis_spec(basis::RadialBasis) = basis.spec

family(basis::UniformBasis) = basis.spec.family_value
family(basis::HalfLineBasis) = basis.spec.family_value
family(basis::RadialBasis) = basis.spec.family_value

mapping(::UniformBasis) = IdentityMapping()
mapping(basis::HalfLineBasis) = basis.spec.mapping_value
mapping(basis::RadialBasis) = basis.spec.mapping_value

centers(basis::UniformBasis) = basis.center_data
centers(basis::HalfLineBasis) = basis.center_data
centers(basis::RadialBasis) = basis.center_data

reference_centers(basis::UniformBasis) = basis.center_data
reference_centers(basis::HalfLineBasis) = basis.reference_center_data
reference_centers(basis::RadialBasis) = basis.reference_center_data

integral_weights(basis::UniformBasis) = basis.integral_weight_data
integral_weights(basis::HalfLineBasis) = basis.integral_weight_data
integral_weights(basis::RadialBasis) = basis.integral_weight_data

stencil(f::BoundaryGausslet) = f.basis.stencil_data[f.index]
stencil(f::RadialGausslet) = f.basis.stencil_data[f.index]

center(f::BoundaryGausslet) = f.basis.center_data[f.index]
center(f::RadialGausslet) = f.basis.center_data[f.index]

reference_center(f::BoundaryGausslet) = f.basis.reference_center_data[f.index]
reference_center(f::RadialGausslet) = f.basis.reference_center_data[f.index]

integral_weight(f::BoundaryGausslet) = f.basis.integral_weight_data[f.index]
integral_weight(f::RadialGausslet) = f.basis.integral_weight_data[f.index]

function _uniform_centers(xmin::Float64, xmax::Float64, spacing::Float64)
    centers_out = Float64[]
    x = xmin
    while x <= xmax + 1.0e-12 * spacing
        push!(centers_out, x)
        x += spacing
    end
    return centers_out
end

"""
    build_basis(spec::UniformBasisSpec)
    build_basis(spec::HalfLineBasisSpec)
    build_basis(spec::RadialBasisSpec)

Build the concrete basis described by `spec`.
"""
function build_basis(spec::UniformBasisSpec)
    center_data = _uniform_centers(spec.xmin, spec.xmax, spec.spacing)
    weight = integral_weight(Gausslet(spec.family_value; center = 0.0, spacing = spec.spacing))
    integral_weight_data = fill(weight, length(center_data))
    return UniformBasis(spec, center_data, integral_weight_data)
end

function build_basis(spec::HalfLineBasisSpec)
    js, coefficients_seed, reference_center_data = _build_halfline_coefficients(spec)
    stencil_data = _boundary_stencils(js, coefficients_seed, spec.family_value, spec.reference_spacing, spec.mapping_value)
    center_data = _physical_centers(reference_center_data, spec.mapping_value)
    integral_weight_data = _integral_weights_from_stencils(stencil_data)
    return HalfLineBasis(spec, stencil_data, collect(reference_center_data), center_data, integral_weight_data)
end

function build_basis(spec::RadialBasisSpec)
    js, coefficients_full, reference_center_data = _build_radial_coefficients(spec)
    stencil_data = _radial_stencils(
        js,
        coefficients_full,
        spec.family_value,
        spec.reference_spacing,
        spec.xgaussians,
        spec.mapping_value,
    )
    center_data = _physical_centers(reference_center_data, spec.mapping_value)
    integral_weight_data = _integral_weights_from_stencils(stencil_data)
    return RadialBasis(spec, stencil_data, collect(reference_center_data), center_data, integral_weight_data)
end
