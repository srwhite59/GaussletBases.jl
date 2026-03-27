module _RadialPrototypeHighPrecFamilies
include(joinpath(@__DIR__, "internal", "families_high_prec.jl"))
end

const _RADIAL_BOUNDARY_PROTOTYPE_NAME = :paper_parity_g10_k6_x2
const _RADIAL_BOUNDARY_PROTOTYPE_H = 1.0e-3
const _RADIAL_BOUNDARY_PROTOTYPE_SIGMA = 3.0
const _RADIAL_BOUNDARY_PROTOTYPE_S0 = 6.5
const _RADIAL_BOUNDARY_PROTOTYPE_RMAX_INT = 80.0
const _RADIAL_BOUNDARY_PROTOTYPE_REFERENCE_SPACING = 1.0
const _RADIAL_BOUNDARY_PROTOTYPE_ODD_HALF_WIDTH = 24
const _RADIAL_BOUNDARY_PROTOTYPE_EVEN_K = 6
const _RADIAL_BOUNDARY_PROTOTYPE_EXPECTED_DIM = 32
const _RADIAL_BOUNDARY_PROTOTYPE_PRECISION_BITS = 256
const _RADIAL_BOUNDARY_PROTOTYPE_RELAXED_EIG_FLOOR = 1.0e-40

const _RADIAL_BOUNDARY_PROTOTYPE_CACHE =
    Ref{Union{Nothing, Dict{Symbol,Any}}}(nothing)

struct RadialBoundaryPrototype
    name::Symbol
    family_value::GaussletFamily
    reference_spacing::Float64
    odd_seed_half_width::Int
    even_tail_kmax::Int
    xgaussians::Vector{XGaussian}
    seed_indices::Vector{Int}
    canonical_coefficients_big::Matrix{BigFloat}
    canonical_coefficients::Matrix{Float64}
    reference_centers_big::Vector{BigFloat}
    reference_centers::Vector{Float64}
    runtime_primitives::Vector{AbstractPrimitiveFunction1D}
    runtime_coefficients::Matrix{Float64}
    stage_dimensions::NamedTuple
    diagnostics::NamedTuple
    checksums::NamedTuple
    provenance::NamedTuple
end

function Base.show(io::IO, prototype::RadialBoundaryPrototype)
    print(
        io,
        "RadialBoundaryPrototype(name=:",
        prototype.name,
        ", family=",
        repr(prototype.family_value.name),
        ", dim=",
        prototype.stage_dimensions.final_dimension,
        ", strict=true)",
    )
end

radial_boundary_prototype_names() = [_RADIAL_BOUNDARY_PROTOTYPE_NAME]

function _paper_parity_xgaussians()
    return XGaussian[
        XGaussian(alpha = 0.09358986806),
        XGaussian(alpha = 0.02357750369),
    ]
end

function _paper_parity_runtime_spec(mapping_value::AbstractCoordinateMapping)
    return RadialBasisSpec(
        :G10;
        count = _RADIAL_BOUNDARY_PROTOTYPE_EXPECTED_DIM,
        mapping = mapping_value,
        reference_spacing = _RADIAL_BOUNDARY_PROTOTYPE_REFERENCE_SPACING,
        tails = 0,
        odd_even_kmax = _RADIAL_BOUNDARY_PROTOTYPE_EVEN_K,
        xgaussians = _paper_parity_xgaussians(),
        rmax_count_policy = :legacy_strict_trim,
    )
end

function _paper_parity_full_runtime_spec(
    prototype::RadialBoundaryPrototype,
    rmax::Real,
    mapping_value::AbstractCoordinateMapping;
    rmax_count_policy::Symbol = :legacy_strict_trim,
)
    return RadialBasisSpec(
        prototype.family_value;
        rmax = rmax,
        mapping = mapping_value,
        reference_spacing = prototype.reference_spacing,
        tails = prototype.even_tail_kmax,
        odd_even_kmax = prototype.even_tail_kmax,
        xgaussians = prototype.xgaussians,
        rmax_count_policy = rmax_count_policy,
    )
end

function _radial_boundary_prototype_path(name::Symbol)
    name == _RADIAL_BOUNDARY_PROTOTYPE_NAME ||
        throw(
            ArgumentError(
                "no cached radial boundary prototype is available for $(repr(name)); available names are $(join(string.(radial_boundary_prototype_names()), ", "))",
            ),
        )
    return normpath(
        joinpath(@__DIR__, "..", "data", "radial", "paper_parity_g10_k6_x2.jld2"),
    )
end

function _prototype_digest_float_vector(values::AbstractVector{<:Real})
    text = join((bitstring(Float64(value)) for value in values), "\n")
    return bytes2hex(SHA.sha1(text))
end

function _prototype_digest_float_matrix(values::AbstractMatrix{<:Real})
    text =
        join((bitstring(Float64(value)) for value in vec(Matrix{Float64}(values))), "\n")
    return bytes2hex(SHA.sha1(text))
end

function _prototype_digest_sample_matrix(basis, points::AbstractVector{Float64})
    return _prototype_digest_float_matrix(_basis_values_matrix(basis, points))
end

function _paper_parity_high_prec_positive_coefficients()
    coeffs =
        _RadialPrototypeHighPrecFamilies._HIGH_PREC_POSITIVE_COEFFICIENTS[:G10]
    return BigFloat[BigFloat(value) for value in coeffs]
end

function _high_prec_full_family_coefficients(
    positive_coefficients::AbstractVector{BigFloat},
)
    tail = reverse(positive_coefficients[2:end])
    return vcat(tail, collect(positive_coefficients))
end

function _paper_parity_seed_primitive_layer_big()
    js = collect(
        -_RADIAL_BOUNDARY_PROTOTYPE_ODD_HALF_WIDTH:
        _RADIAL_BOUNDARY_PROTOTYPE_ODD_HALF_WIDTH,
    )
    positive = _paper_parity_high_prec_positive_coefficients()
    full = _high_prec_full_family_coefficients(positive)
    radius = length(positive) - 1
    primitive_indices = collect((3 * first(js) - radius):(3 * last(js) + radius))
    primitive_spacing = BigFloat(_RADIAL_BOUNDARY_PROTOTYPE_REFERENCE_SPACING) / 3
    scale = inv(sqrt(2 * big(pi) * BigFloat(_RADIAL_BOUNDARY_PROTOTYPE_REFERENCE_SPACING)))

    matrix = zeros(BigFloat, length(primitive_indices), length(js))
    for column in eachindex(js)
        j = js[column]
        for row in eachindex(primitive_indices)
            local_offset = primitive_indices[row] - 3 * j
            if -radius <= local_offset <= radius
                matrix[row, column] =
                    scale * full[local_offset + radius + 1]
            end
        end
    end

    primitive_centers =
        [BigFloat(index) * primitive_spacing for index in primitive_indices]
    return (
        js = js,
        primitive_indices = primitive_indices,
        primitive_centers = collect(primitive_centers),
        primitive_width = primitive_spacing,
        primitive_matrix = matrix,
    )
end

function _paper_parity_runtime_primitives()
    seed = _paper_parity_seed_primitive_layer_big()
    primitive_width = Float64(seed.primitive_width)
    primitives = AbstractPrimitiveFunction1D[
        HalfLineGaussian(center = Float64(center_value), width = primitive_width) for
        center_value in seed.primitive_centers
    ]
    append!(primitives, AbstractPrimitiveFunction1D[x for x in _paper_parity_xgaussians()])
    return primitives
end

function _bigfloat_sqrt_pi_over_two()
    return sqrt(big(pi) / 2)
end

function _bigfloat_sqrt_two()
    return sqrt(BigFloat(2))
end

function _truncated_gaussian_moments(mu::BigFloat, sigma::BigFloat)
    sigma > 0 || throw(ArgumentError("truncated Gaussian moments require sigma > 0"))
    expo = exp(-(mu^2) / (2 * sigma^2))
    moment0 =
        sigma * _bigfloat_sqrt_pi_over_two() *
        erfc(-mu / (_bigfloat_sqrt_two() * sigma))
    moment1 = mu * moment0 + sigma^2 * expo
    moment2 = (mu^2 + sigma^2) * moment0 + mu * sigma^2 * expo
    moment3 =
        (mu^3 + 3 * mu * sigma^2) * moment0 +
        sigma^2 * (mu^2 + 2 * sigma^2) * expo
    return moment0, moment1, moment2, moment3
end

function _gaussian_product_parameters(
    center_a::Real,
    width_a::Real,
    center_b::Real,
    width_b::Real,
)
    ca = BigFloat(center_a)
    wa = BigFloat(width_a)
    cb = BigFloat(center_b)
    wb = BigFloat(width_b)
    sigma2 = (wa^2 * wb^2) / (wa^2 + wb^2)
    mu = (ca * wb^2 + cb * wa^2) / (wa^2 + wb^2)
    prefactor = exp(-((ca - cb)^2) / (2 * (wa^2 + wb^2)))
    return mu, sqrt(sigma2), prefactor
end

function _analytic_overlap_position(
    a::HalfLineGaussian,
    b::HalfLineGaussian,
)
    mu, sigma, prefactor =
        _gaussian_product_parameters(a.center_value, a.width, b.center_value, b.width)
    moment0, moment1, _, _ = _truncated_gaussian_moments(mu, sigma)
    return prefactor * moment0, prefactor * moment1
end

function _analytic_overlap_position(a::HalfLineGaussian, b::XGaussian)
    mu, sigma, prefactor =
        _gaussian_product_parameters(a.center_value, a.width, 0.0, b.alpha)
    _, moment1, moment2, _ = _truncated_gaussian_moments(mu, sigma)
    return prefactor * moment1, prefactor * moment2
end

function _analytic_overlap_position(a::XGaussian, b::HalfLineGaussian)
    overlap_value, position_value = _analytic_overlap_position(b, a)
    return overlap_value, position_value
end

function _analytic_overlap_position(a::XGaussian, b::XGaussian)
    mu, sigma, prefactor = _gaussian_product_parameters(0.0, a.alpha, 0.0, b.alpha)
    _, _, moment2, moment3 = _truncated_gaussian_moments(mu, sigma)
    return prefactor * moment2, prefactor * moment3
end

function _paper_parity_analytic_primitive_matrices(
    primitive_data::AbstractVector{<:AbstractPrimitiveFunction1D},
)
    nprimitive = length(primitive_data)
    overlap = zeros(BigFloat, nprimitive, nprimitive)
    position = zeros(BigFloat, nprimitive, nprimitive)
    for a in 1:nprimitive
        pa = primitive_data[a]
        for b in a:nprimitive
            pb = primitive_data[b]
            overlap_ab, position_ab = _analytic_overlap_position(pa, pb)
            overlap[a, b] = overlap_ab
            overlap[b, a] = overlap_ab
            position[a, b] = position_ab
            position[b, a] = position_ab
        end
    end
    return overlap, position
end

function _paper_parity_primitive_integral_weights(
    primitive_data::AbstractVector{<:AbstractPrimitiveFunction1D},
)
    weights = Vector{BigFloat}(undef, length(primitive_data))
    sqrt_pi_over_two = _bigfloat_sqrt_pi_over_two()
    sqrt_two = _bigfloat_sqrt_two()
    for index in eachindex(primitive_data)
        primitive = primitive_data[index]
        if primitive isa HalfLineGaussian
            weights[index] =
                BigFloat(primitive.width) *
                sqrt_pi_over_two *
                erfc((-BigFloat(primitive.center_value)) /
                     (sqrt_two * BigFloat(primitive.width)))
        elseif primitive isa XGaussian
            weights[index] = BigFloat(primitive.alpha)^2
        else
            throw(
                ArgumentError(
                    "paper-parity prototype only supports HalfLineGaussian and XGaussian primitives",
                ),
            )
        end
    end
    return weights
end

function _paper_parity_primitive_values_at_zero(
    primitive_data::AbstractVector{<:AbstractPrimitiveFunction1D},
)
    values = Vector{BigFloat}(undef, length(primitive_data))
    for index in eachindex(primitive_data)
        primitive = primitive_data[index]
        if primitive isa HalfLineGaussian
            z = -BigFloat(primitive.center_value) / BigFloat(primitive.width)
            values[index] = exp(-(z^2) / 2)
        elseif primitive isa XGaussian
            values[index] = BigFloat(0)
        else
            throw(
                ArgumentError(
                    "paper-parity prototype only supports HalfLineGaussian and XGaussian primitives",
                ),
            )
        end
    end
    return values
end

function _prototype_eigen_floor()
    return setprecision(_RADIAL_BOUNDARY_PROTOTYPE_PRECISION_BITS) do
        BigFloat(_RADIAL_BOUNDARY_PROTOTYPE_RELAXED_EIG_FLOOR)
    end
end

function _jacobi_eigen_symmetric(
    matrix::AbstractMatrix{BigFloat};
    tol::Union{Nothing, BigFloat} = nothing,
    max_sweeps::Int = 256,
)
    size(matrix, 1) == size(matrix, 2) ||
        throw(ArgumentError("Jacobi eigensolver requires a square matrix"))
    n = size(matrix, 1)
    work = Matrix{BigFloat}(matrix)
    vectors = Matrix{BigFloat}(I, n, n)
    tolerance = tol === nothing ? sqrt(eps(BigFloat)) : tol

    for _ in 1:max_sweeps
        offdiag = zero(BigFloat)
        for p in 1:(n - 1), q in (p + 1):n
            offdiag = max(offdiag, abs(work[p, q]))
        end
        offdiag <= tolerance && break

        for p in 1:(n - 1), q in (p + 1):n
            apq = work[p, q]
            abs(apq) <= tolerance && continue

            app = work[p, p]
            aqq = work[q, q]
            tau = (aqq - app) / (2 * apq)
            t =
                tau == 0 ?
                one(BigFloat) :
                sign(tau) / (abs(tau) + sqrt(one(BigFloat) + tau^2))
            c = inv(sqrt(one(BigFloat) + t^2))
            s = t * c

            for k in 1:n
                if k != p && k != q
                    akp = work[k, p]
                    akq = work[k, q]
                    work[k, p] = c * akp - s * akq
                    work[p, k] = work[k, p]
                    work[k, q] = s * akp + c * akq
                    work[q, k] = work[k, q]
                end
            end

            work[p, p] = c^2 * app - 2 * s * c * apq + s^2 * aqq
            work[q, q] = s^2 * app + 2 * s * c * apq + c^2 * aqq
            work[p, q] = zero(BigFloat)
            work[q, p] = zero(BigFloat)

            for k in 1:n
                vkp = vectors[k, p]
                vkq = vectors[k, q]
                vectors[k, p] = c * vkp - s * vkq
                vectors[k, q] = s * vkp + c * vkq
            end
        end
    end

    values = BigFloat[work[index, index] for index in 1:n]
    ordering = sortperm(values)
    return Matrix{BigFloat}(vectors[:, ordering]), values[ordering]
end

function _strict_inverse_sqrt(overlap::AbstractMatrix{BigFloat}; label::AbstractString)
    vectors, values = _jacobi_eigen_symmetric(Matrix{BigFloat}(overlap))
    smallest = minimum(values)
    smallest > _prototype_eigen_floor() || throw(
        ArgumentError(
            string(
                label,
                " overlap lost a mode under the strict manuscript prototype contract; smallest eigenvalues are ",
                repr(Float64.(values[1:min(end, 6)])),
            ),
        ),
    )
    invhalf = Diagonal(inv.(sqrt.(values)))
    return vectors, invhalf, values
end

function _sign_fix_columns_big!(
    coefficient_matrix::Matrix{BigFloat},
    basis_integrals::Vector{BigFloat},
)
    integrals = vec(transpose(basis_integrals) * coefficient_matrix)
    for column in 1:size(coefficient_matrix, 2)
        if integrals[column] < 0
            coefficient_matrix[:, column] .*= -1
        end
    end
    return coefficient_matrix
end

function _paper_parity_stage_dimensions(runtime_primitive_count::Int)
    return (
        seed_count = 2 * _RADIAL_BOUNDARY_PROTOTYPE_ODD_HALF_WIDTH + 1,
        raw_odd_count = _RADIAL_BOUNDARY_PROTOTYPE_ODD_HALF_WIDTH,
        raw_even_count = _RADIAL_BOUNDARY_PROTOTYPE_EVEN_K + 1,
        cleaned_even_count = _RADIAL_BOUNDARY_PROTOTYPE_EVEN_K,
        xgaussian_count = length(_paper_parity_xgaussians()),
        canonical_base_count = 2 * _RADIAL_BOUNDARY_PROTOTYPE_ODD_HALF_WIDTH + 1 + length(_paper_parity_xgaussians()),
        final_dimension = _RADIAL_BOUNDARY_PROTOTYPE_EXPECTED_DIM,
        runtime_primitive_count = runtime_primitive_count,
    )
end

function _paper_parity_sample_points(; step::Float64 = 0.01, umax::Float64 = 10.0)
    return Float64[step * i for i in 0:round(Int, umax / step)]
end

function _paper_parity_prototype_runtime_basis(
    runtime_primitives::Vector{AbstractPrimitiveFunction1D},
    runtime_coefficients::AbstractMatrix{<:Real},
    reference_centers::AbstractVector{<:Real};
    mapping_value::AbstractCoordinateMapping = IdentityMapping(),
    spec::RadialBasisSpec = _paper_parity_runtime_spec(mapping_value),
    build_umax::Float64 = _RADIAL_BOUNDARY_PROTOTYPE_RMAX_INT,
)
    coefficient_matrix = Matrix{Float64}(runtime_coefficients)
    reference_center_data = Float64[Float64(value) for value in reference_centers]
    primitive_data = _with_mapping(runtime_primitives, mapping_value)
    center_data = _physical_centers(reference_center_data, mapping_value)
    integral_weight_data =
        _integral_weights_from_representation(primitive_data, coefficient_matrix)
    return RadialBasis(
        spec,
        primitive_data,
        coefficient_matrix,
        reference_center_data,
        center_data,
        integral_weight_data,
        build_umax,
    )
end

function _paper_parity_halfline_primitive_count(prototype::RadialBoundaryPrototype)
    return length(prototype.runtime_primitives) - length(prototype.xgaussians)
end

function _halfline_primitive_grid_index(
    primitive::HalfLineGaussian,
    primitive_spacing::Float64,
)
    grid_index = round(Int, primitive.center_value / primitive_spacing)
    abs(primitive.center_value - grid_index * primitive_spacing) <= 1.0e-12 ||
        throw(
            ArgumentError(
                "paper-parity radial tail extension expected half-line primitive centers on the primitive grid; got center=$(primitive.center_value), primitive_spacing=$(primitive_spacing)",
            ),
        )
    return grid_index
end

function _paper_parity_tail_extended_runtime_basis(
    prototype::RadialBoundaryPrototype;
    mapping_value::AbstractCoordinateMapping,
    rmax::Real,
    rmax_count_policy::Symbol = :legacy_strict_trim,
)
    full_spec = _paper_parity_full_runtime_spec(
        prototype,
        rmax,
        mapping_value;
        rmax_count_policy = rmax_count_policy,
    )
    target_odd_count =
        _radial_rmax_target_odd_count(full_spec, prototype.reference_spacing)
    ninj = length(prototype.xgaussians)
    target_count = target_odd_count + prototype.even_tail_kmax + ninj
    prototype_dimension = prototype.stage_dimensions.final_dimension
    target_count >= 1 ||
        throw(ArgumentError("paper-parity radial extension requires a positive target count"))

    if target_count <= prototype_dimension
        truncated_coefficients = prototype.runtime_coefficients[:, 1:target_count]
        truncated_centers = prototype.reference_centers[1:target_count]
        build_umax = max(
            _RADIAL_BOUNDARY_PROTOTYPE_RMAX_INT,
            _gausslet_reference_support_umax(
                truncated_centers,
                prototype.family_value,
                prototype.reference_spacing,
            ),
            _xgaussian_reference_support_umax(prototype.xgaussians),
        )
        return _paper_parity_prototype_runtime_basis(
            prototype.runtime_primitives,
            truncated_coefficients,
            truncated_centers;
            mapping_value = mapping_value,
            spec = full_spec,
            build_umax = build_umax,
        )
    end

    target_odd_count > prototype.odd_seed_half_width ||
        throw(ArgumentError("paper-parity radial extension expected the odd count to exceed the cached boundary prototype before tail appending"))

    js_tail = collect((prototype.odd_seed_half_width + 1):target_odd_count)
    tail_primitive_matrix, tail_primitives_ref =
        _halfline_seed_primitive_layer(
            js_tail,
            prototype.family_value,
            prototype.reference_spacing,
        )

    primitive_spacing = prototype.reference_spacing / 3.0
    nprototype_halfline = _paper_parity_halfline_primitive_count(prototype)
    prototype_halfline = prototype.runtime_primitives[1:nprototype_halfline]
    prototype_xgaussians = prototype.runtime_primitives[(nprototype_halfline + 1):end]

    merged_halfline = AbstractPrimitiveFunction1D[primitive for primitive in prototype_halfline]
    merged_halfline_index = Dict{Int, Int}()
    for primitive_index in eachindex(prototype_halfline)
        primitive = prototype_halfline[primitive_index]
        primitive isa HalfLineGaussian ||
            throw(ArgumentError("paper-parity cached runtime primitives must begin with half-line primitives"))
        merged_halfline_index[_halfline_primitive_grid_index(primitive, primitive_spacing)] =
            primitive_index
    end

    tail_row_to_merged = Vector{Int}(undef, length(tail_primitives_ref))
    for tail_row in eachindex(tail_primitives_ref)
        primitive = tail_primitives_ref[tail_row]
        primitive isa HalfLineGaussian ||
            throw(ArgumentError("paper-parity tail extension expected half-line primitive data"))
        grid_index = _halfline_primitive_grid_index(primitive, primitive_spacing)
        merged_index = get(merged_halfline_index, grid_index, 0)
        if merged_index == 0
            push!(merged_halfline, primitive)
            merged_index = length(merged_halfline)
            merged_halfline_index[grid_index] = merged_index
        end
        tail_row_to_merged[tail_row] = merged_index
    end

    merged_primitives_ref = vcat(merged_halfline, prototype_xgaussians)
    merged_coefficients =
        zeros(
            Float64,
            length(merged_primitives_ref),
            prototype_dimension + length(js_tail),
        )

    merged_coefficients[1:nprototype_halfline, 1:prototype_dimension] .=
        prototype.runtime_coefficients[1:nprototype_halfline, :]
    merged_x_offset = length(merged_halfline)
    merged_coefficients[(merged_x_offset + 1):end, 1:prototype_dimension] .=
        prototype.runtime_coefficients[(nprototype_halfline + 1):end, :]

    for tail_row in eachindex(tail_row_to_merged)
        merged_row = tail_row_to_merged[tail_row]
        @views merged_coefficients[merged_row, (prototype_dimension + 1):end] .+=
            tail_primitive_matrix[tail_row, :]
    end

    extended_reference_centers = vcat(
        prototype.reference_centers,
        Float64[prototype.reference_spacing * index for index in js_tail],
    )
    build_umax = max(
        _RADIAL_BOUNDARY_PROTOTYPE_RMAX_INT,
        _gausslet_reference_support_umax(
            extended_reference_centers,
            prototype.family_value,
            prototype.reference_spacing,
        ),
        _xgaussian_reference_support_umax(prototype.xgaussians),
    )
    return _paper_parity_prototype_runtime_basis(
        merged_primitives_ref,
        merged_coefficients,
        extended_reference_centers;
        mapping_value = mapping_value,
        spec = full_spec,
        build_umax = build_umax,
    )
end

function _paper_parity_prototype_grid(mapping_value::AbstractCoordinateMapping)
    points, weights = _make_erf_grid(
        ;
        h = _RADIAL_BOUNDARY_PROTOTYPE_H,
        rmax = _RADIAL_BOUNDARY_PROTOTYPE_RMAX_INT,
        sigma = _RADIAL_BOUNDARY_PROTOTYPE_SIGMA,
        s0 = _RADIAL_BOUNDARY_PROTOTYPE_S0,
    )
    return RadialQuadratureGrid(points, weights; mapping = mapping_value)
end

function _build_paper_parity_canonical_analytic()
    setprecision(_RADIAL_BOUNDARY_PROTOTYPE_PRECISION_BITS) do
        seed = _paper_parity_seed_primitive_layer_big()
        runtime_primitives = _paper_parity_runtime_primitives()
        primitive_overlap, primitive_position =
            _paper_parity_analytic_primitive_matrices(runtime_primitives)
        primitive_integrals =
            _paper_parity_primitive_integral_weights(runtime_primitives)
        primitive_zero_values =
            _paper_parity_primitive_values_at_zero(runtime_primitives)

        nseed = length(seed.js)
        ninj = length(_paper_parity_xgaussians())
        primitive_seed_matrix = seed.primitive_matrix
        seedfull_matrix = zeros(
            BigFloat,
            size(primitive_seed_matrix, 1) + ninj,
            nseed + ninj,
        )
        seedfull_matrix[1:size(primitive_seed_matrix, 1), 1:nseed] .=
            primitive_seed_matrix
        seedfull_matrix[(size(primitive_seed_matrix, 1) + 1):end, (nseed + 1):end] .=
            Matrix{BigFloat}(I, ninj, ninj)

        base_overlap =
            Matrix{BigFloat}(transpose(seedfull_matrix) * primitive_overlap * seedfull_matrix)
        base_position =
            Matrix{BigFloat}(transpose(seedfull_matrix) * primitive_position * seedfull_matrix)
        base_integrals = vec(transpose(seedfull_matrix) * primitive_integrals)
        seed_zero_values =
            vec(transpose(primitive_seed_matrix) * primitive_zero_values[1:size(primitive_seed_matrix, 1)])

        lseed = _RADIAL_BOUNDARY_PROTOTYPE_ODD_HALF_WIDTH
        odd_block = zeros(BigFloat, nseed, lseed)
        for k in 1:lseed
            ip = (lseed + 1) + k
            im = (lseed + 1) - k
            odd_block[ip, k] = 1
            odd_block[im, k] = -1
        end

        even_block = zeros(BigFloat, nseed, _RADIAL_BOUNDARY_PROTOTYPE_EVEN_K + 1)
        even_block[lseed + 1, 1] = 1
        for k in 1:_RADIAL_BOUNDARY_PROTOTYPE_EVEN_K
            ip = (lseed + 1) + k
            im = (lseed + 1) - k
            even_block[ip, k + 1] = 1
            even_block[im, k + 1] = 1
        end

        seed_overlap = @view base_overlap[1:nseed, 1:nseed]
        seed_position = @view base_position[1:nseed, 1:nseed]

        function normalize_block(block::Matrix{BigFloat})
            overlap = seed_overlap * block
            norms = vec(sum(block .* overlap; dims = 1))
            any(norms .<= 0) && throw(
                ArgumentError(
                    "paper-parity prototype encountered a nonpositive block norm",
                ),
            )
            return block * Diagonal(inv.(sqrt.(norms)))
        end

        odd_normalized = normalize_block(odd_block)
        odd_overlap = Matrix{BigFloat}(transpose(odd_normalized) * seed_overlap * odd_normalized)
        odd_position = Matrix{BigFloat}(transpose(odd_normalized) * seed_position * odd_normalized)
        odd_vectors, odd_invhalf, _ = _strict_inverse_sqrt(
            odd_overlap;
            label = "paper-parity odd block",
        )
        odd_localizer_matrix =
            Matrix{BigFloat}(odd_invhalf * (transpose(odd_vectors) * odd_position * odd_vectors) * odd_invhalf)
        odd_localizer_vectors, _ = _jacobi_eigen_symmetric(odd_localizer_matrix)
        odd_coefficients =
            odd_normalized *
            (odd_vectors * (odd_invhalf * odd_localizer_vectors))

        even_normalized = normalize_block(even_block)
        even_projection =
            transpose(odd_coefficients) * seed_overlap * even_normalized
        even_clean = even_normalized .- (odd_coefficients * even_projection)

        delta_direction = vec(transpose(seed_zero_values) * even_clean)
        delta_norm = dot(delta_direction, delta_direction)
        delta_norm > 0 || throw(
            ArgumentError(
                "paper-parity prototype failed to identify the even-tail delta direction",
            ),
        )
        unit = delta_direction ./ sqrt(delta_norm)
        sign_value = unit[1] >= 0 ? BigFloat(1) : BigFloat(-1)
        reflector = copy(unit)
        reflector[1] += sign_value
        beta = dot(reflector, reflector)
        beta > 0 || throw(
            ArgumentError(
                "paper-parity prototype failed while removing the delta direction",
            ),
        )
        householder =
            Matrix{BigFloat}(I, length(unit), length(unit)) .-
            ((BigFloat(2) / beta) .* (reflector * transpose(reflector)))
        even_clean = even_clean * householder[:, 2:end]

        combined_base = hcat(
            vcat(odd_coefficients, zeros(BigFloat, ninj, size(odd_coefficients, 2))),
            vcat(even_clean, zeros(BigFloat, ninj, size(even_clean, 2))),
            vcat(zeros(BigFloat, nseed, ninj), Matrix{BigFloat}(I, ninj, ninj)),
        )

        combined_overlap =
            Matrix{BigFloat}(transpose(combined_base) * base_overlap * combined_base)
        combined_position =
            Matrix{BigFloat}(transpose(combined_base) * base_position * combined_base)
        overlap_vectors, overlap_invhalf, overlap_values = _strict_inverse_sqrt(
            combined_overlap;
            label = "paper-parity combined pre-COMX block",
        )
        orthonormal_coefficients =
            combined_base * (overlap_vectors * overlap_invhalf)
        localized_position_matrix =
            Matrix{BigFloat}(transpose(orthonormal_coefficients) * base_position * orthonormal_coefficients)
        localized_position_vectors, localized_position_values =
            _jacobi_eigen_symmetric(localized_position_matrix)
        final_coefficients =
            orthonormal_coefficients * localized_position_vectors
        _sign_fix_columns_big!(final_coefficients, base_integrals)

        final_overlap =
            Matrix{BigFloat}(transpose(final_coefficients) * base_overlap * final_coefficients)
        overlap_identity_error =
            Float64(maximum(abs, final_overlap - Matrix{BigFloat}(I, size(final_overlap, 1), size(final_overlap, 2))))
        runtime_coefficients_big = seedfull_matrix * final_coefficients

        runtime_basis_identity = _paper_parity_prototype_runtime_basis(
            _paper_parity_runtime_primitives(),
            Float64.(runtime_coefficients_big),
            Float64.(localized_position_values);
            mapping_value = IdentityMapping(),
        )
        diagnostics_grid = _paper_parity_prototype_grid(IdentityMapping())
        runtime_diag = basis_diagnostics(runtime_basis_identity, diagnostics_grid)
        sample_points = _paper_parity_sample_points()

        stage_dimensions = _paper_parity_stage_dimensions(length(runtime_primitives))
        checksums = (
            reference_centers = _prototype_digest_float_vector(Float64.(localized_position_values)),
            canonical_coefficients = _prototype_digest_float_matrix(Float64.(final_coefficients)),
            runtime_coefficients = _prototype_digest_float_matrix(Float64.(runtime_coefficients_big)),
            sampled_basis_u0_10_du0_01 = _prototype_digest_sample_matrix(runtime_basis_identity, sample_points),
        )
        diagnostics = (
            expected_final_dimension = _RADIAL_BOUNDARY_PROTOTYPE_EXPECTED_DIM,
            retained_dimension = size(final_coefficients, 2),
            mode_drop_count = 0,
            smallest_overlap_eigenvalue = Float64(minimum(overlap_values)),
            overlap_condition = Float64(maximum(overlap_values) / minimum(overlap_values)),
            overlap_identity_error = overlap_identity_error,
            evaluation_overlap_identity_error = runtime_diag.overlap_error,
            D = runtime_diag.D,
            centers_monotone = issorted(Float64.(localized_position_values)),
        )
        provenance = (
            family = :G10,
            reference_spacing = _RADIAL_BOUNDARY_PROTOTYPE_REFERENCE_SPACING,
            odd_seed_half_width = _RADIAL_BOUNDARY_PROTOTYPE_ODD_HALF_WIDTH,
            even_tail_kmax = _RADIAL_BOUNDARY_PROTOTYPE_EVEN_K,
            xgaussian_alphas = Float64[x.alpha for x in _paper_parity_xgaussians()],
            h = _RADIAL_BOUNDARY_PROTOTYPE_H,
            sigma = _RADIAL_BOUNDARY_PROTOTYPE_SIGMA,
            s0 = _RADIAL_BOUNDARY_PROTOTYPE_S0,
            rmax_int = _RADIAL_BOUNDARY_PROTOTYPE_RMAX_INT,
            precision_bits = _RADIAL_BOUNDARY_PROTOTYPE_PRECISION_BITS,
            construction = "analytic_strict",
        )

        return (
            seed_indices = seed.js,
            canonical_coefficients_big = final_coefficients,
            canonical_coefficients = Float64.(final_coefficients),
            reference_centers_big = Vector{BigFloat}(localized_position_values),
            reference_centers = Float64.(localized_position_values),
            runtime_primitives = _paper_parity_runtime_primitives(),
            runtime_coefficients = Float64.(runtime_coefficients_big),
            stage_dimensions = stage_dimensions,
            diagnostics = diagnostics,
            checksums = checksums,
            provenance = provenance,
            primitive_centers = Float64[Float64(value) for value in seed.primitive_centers],
            primitive_width = Float64(seed.primitive_width),
        )
    end
end

function _write_paper_parity_radial_prototype_cache(path::AbstractString)
    data = _build_paper_parity_canonical_analytic()
    mkpath(dirname(path))
    jldopen(path, "w") do file
        file["meta/name"] = String(_RADIAL_BOUNDARY_PROTOTYPE_NAME)
        file["meta/family"] = "G10"
        file["meta/reference_spacing"] =
            _RADIAL_BOUNDARY_PROTOTYPE_REFERENCE_SPACING
        file["meta/odd_seed_half_width"] =
            _RADIAL_BOUNDARY_PROTOTYPE_ODD_HALF_WIDTH
        file["meta/even_tail_kmax"] = _RADIAL_BOUNDARY_PROTOTYPE_EVEN_K
        file["meta/xgaussian_alphas"] = Float64[x.alpha for x in _paper_parity_xgaussians()]
        file["meta/h"] = _RADIAL_BOUNDARY_PROTOTYPE_H
        file["meta/sigma"] = _RADIAL_BOUNDARY_PROTOTYPE_SIGMA
        file["meta/s0"] = _RADIAL_BOUNDARY_PROTOTYPE_S0
        file["meta/rmax_int"] = _RADIAL_BOUNDARY_PROTOTYPE_RMAX_INT
        file["meta/precision_bits"] = _RADIAL_BOUNDARY_PROTOTYPE_PRECISION_BITS

        file["canonical/js"] = data.seed_indices
        file["canonical/final_coefficients_big"] = data.canonical_coefficients_big
        file["canonical/final_coefficients_f64"] = data.canonical_coefficients
        file["canonical/reference_centers_big"] = data.reference_centers_big
        file["canonical/reference_centers_f64"] = data.reference_centers

        file["runtime/primitive_centers"] = data.primitive_centers
        file["runtime/primitive_width"] = data.primitive_width
        file["runtime/final_coefficients_f64"] = data.runtime_coefficients

        for (key, value) in pairs(data.stage_dimensions)
            file["stage_dimensions/$(key)"] = value
        end
        for (key, value) in pairs(data.diagnostics)
            file["diagnostics/$(key)"] = value
        end
        for (key, value) in pairs(data.checksums)
            file["checksums/$(key)"] = value
        end
        for (key, value) in pairs(data.provenance)
            file["provenance/$(key)"] = value
        end
    end
    return path
end

function _load_paper_parity_radial_prototype(path::AbstractString)
    return jldopen(path, "r") do file
        name = Symbol(String(file["meta/name"]))
        seed_indices = Vector{Int}(file["canonical/js"])
        canonical_coefficients_big =
            Matrix{BigFloat}(file["canonical/final_coefficients_big"])
        canonical_coefficients =
            Matrix{Float64}(file["canonical/final_coefficients_f64"])
        reference_centers_big =
            Vector{BigFloat}(file["canonical/reference_centers_big"])
        reference_centers =
            Vector{Float64}(file["canonical/reference_centers_f64"])
        primitive_centers = Vector{Float64}(file["runtime/primitive_centers"])
        primitive_width = Float64(file["runtime/primitive_width"])
        runtime_coefficients =
            Matrix{Float64}(file["runtime/final_coefficients_f64"])

        runtime_primitives = AbstractPrimitiveFunction1D[
            HalfLineGaussian(center = center_value, width = primitive_width) for
            center_value in primitive_centers
        ]
        append!(runtime_primitives, AbstractPrimitiveFunction1D[x for x in _paper_parity_xgaussians()])

        stage_dimensions = (
            seed_count = Int(file["stage_dimensions/seed_count"]),
            raw_odd_count = Int(file["stage_dimensions/raw_odd_count"]),
            raw_even_count = Int(file["stage_dimensions/raw_even_count"]),
            cleaned_even_count = Int(file["stage_dimensions/cleaned_even_count"]),
            xgaussian_count = Int(file["stage_dimensions/xgaussian_count"]),
            canonical_base_count = Int(file["stage_dimensions/canonical_base_count"]),
            final_dimension = Int(file["stage_dimensions/final_dimension"]),
            runtime_primitive_count = Int(file["stage_dimensions/runtime_primitive_count"]),
        )
        diagnostics = (
            expected_final_dimension = Int(file["diagnostics/expected_final_dimension"]),
            retained_dimension = Int(file["diagnostics/retained_dimension"]),
            mode_drop_count = Int(file["diagnostics/mode_drop_count"]),
            smallest_overlap_eigenvalue = Float64(file["diagnostics/smallest_overlap_eigenvalue"]),
            overlap_condition = Float64(file["diagnostics/overlap_condition"]),
            overlap_identity_error = Float64(file["diagnostics/overlap_identity_error"]),
            evaluation_overlap_identity_error = Float64(file["diagnostics/evaluation_overlap_identity_error"]),
            D = Float64(file["diagnostics/D"]),
            centers_monotone = Bool(file["diagnostics/centers_monotone"]),
        )
        diagnostics.retained_dimension == diagnostics.expected_final_dimension ||
            throw(
                ArgumentError(
                    "cached paper-parity radial prototype lost modes; retained $(diagnostics.retained_dimension) but expected $(diagnostics.expected_final_dimension)",
                ),
            )
        diagnostics.mode_drop_count == 0 ||
            throw(
                ArgumentError(
                    "cached paper-parity radial prototype reports nonzero mode_drop_count=$(diagnostics.mode_drop_count)",
                ),
            )

        checksums = (
            reference_centers = String(file["checksums/reference_centers"]),
            canonical_coefficients = String(file["checksums/canonical_coefficients"]),
            runtime_coefficients = String(file["checksums/runtime_coefficients"]),
            sampled_basis_u0_10_du0_01 = String(file["checksums/sampled_basis_u0_10_du0_01"]),
        )
        provenance = (
            family = Symbol(String(file["provenance/family"])),
            reference_spacing = Float64(file["provenance/reference_spacing"]),
            odd_seed_half_width = Int(file["provenance/odd_seed_half_width"]),
            even_tail_kmax = Int(file["provenance/even_tail_kmax"]),
            xgaussian_alphas = Vector{Float64}(file["provenance/xgaussian_alphas"]),
            h = Float64(file["provenance/h"]),
            sigma = Float64(file["provenance/sigma"]),
            s0 = Float64(file["provenance/s0"]),
            rmax_int = Float64(file["provenance/rmax_int"]),
            precision_bits = Int(file["provenance/precision_bits"]),
            construction = String(file["provenance/construction"]),
        )

        return RadialBoundaryPrototype(
            name,
            GaussletFamily(:G10),
            Float64(file["meta/reference_spacing"]),
            Int(file["meta/odd_seed_half_width"]),
            Int(file["meta/even_tail_kmax"]),
            _paper_parity_xgaussians(),
            seed_indices,
            canonical_coefficients_big,
            canonical_coefficients,
            reference_centers_big,
            reference_centers,
            runtime_primitives,
            runtime_coefficients,
            stage_dimensions,
            diagnostics,
            checksums,
            provenance,
        )
    end
end

function _radial_boundary_prototype_cache()
    cache = _RADIAL_BOUNDARY_PROTOTYPE_CACHE[]
    if cache === nothing
        prototype =
            _load_paper_parity_radial_prototype(
                _radial_boundary_prototype_path(_RADIAL_BOUNDARY_PROTOTYPE_NAME),
            )
        cache = Dict{Symbol, Any}(_RADIAL_BOUNDARY_PROTOTYPE_NAME => prototype)
        _RADIAL_BOUNDARY_PROTOTYPE_CACHE[] = cache
    end
    return cache
end

function radial_boundary_prototype(name::Symbol = _RADIAL_BOUNDARY_PROTOTYPE_NAME)
    cache = _radial_boundary_prototype_cache()
    prototype = get(cache, name, nothing)
    prototype === nothing &&
        throw(
            ArgumentError(
                "no cached radial boundary prototype is available for $(repr(name)); available names are $(join(string.(collect(keys(cache))), ", "))",
            ),
        )
    return prototype
end

function build_paper_parity_radial_basis(
    prototype::RadialBoundaryPrototype;
    rmax::Real,
    mapping::AbstractCoordinateMapping = IdentityMapping(),
    rmax_count_policy::Symbol = :legacy_strict_trim,
)
    return _paper_parity_tail_extended_runtime_basis(
        prototype;
        mapping_value = mapping,
        rmax = rmax,
        rmax_count_policy = rmax_count_policy,
    )
end

function build_paper_parity_radial_basis(
    ;
    rmax::Real,
    mapping::AbstractCoordinateMapping = IdentityMapping(),
    prototype_name::Symbol = _RADIAL_BOUNDARY_PROTOTYPE_NAME,
    rmax_count_policy::Symbol = :legacy_strict_trim,
)
    return build_paper_parity_radial_basis(
        radial_boundary_prototype(prototype_name);
        rmax = rmax,
        mapping = mapping,
        rmax_count_policy = rmax_count_policy,
    )
end

function build_basis(
    prototype::RadialBoundaryPrototype;
    mapping::AbstractCoordinateMapping = IdentityMapping(),
    rmax::Union{Nothing, Real} = nothing,
    rmax_count_policy::Symbol = :legacy_strict_trim,
)
    if rmax !== nothing
        return build_paper_parity_radial_basis(
            prototype;
            rmax = rmax,
            mapping = mapping,
            rmax_count_policy = rmax_count_policy,
        )
    end
    return _paper_parity_prototype_runtime_basis(
        prototype.runtime_primitives,
        prototype.runtime_coefficients,
        prototype.reference_centers;
        mapping_value = mapping,
    )
end

function _paper_parity_numerical_reference_basis(
    ;
    mapping::AbstractCoordinateMapping = IdentityMapping(),
)
    spacing = _RADIAL_BOUNDARY_PROTOTYPE_REFERENCE_SPACING
    xgaussians = _paper_parity_xgaussians()
    lseed = _RADIAL_BOUNDARY_PROTOTYPE_ODD_HALF_WIDTH
    js = collect(-lseed:lseed)

    xgrid, weights = _make_erf_grid(
        ;
        h = _RADIAL_BOUNDARY_PROTOTYPE_H,
        rmax = _RADIAL_BOUNDARY_PROTOTYPE_RMAX_INT,
        sigma = _RADIAL_BOUNDARY_PROTOTYPE_SIGMA,
        s0 = _RADIAL_BOUNDARY_PROTOTYPE_S0,
    )
    xweights = xgrid .* weights
    seed_sampler, sampled_seed_intervals =
        _sample_shifted_gausslets(GaussletFamily(:G10), js, xgrid, spacing)
    overlap_seed = _interval_gram_matrix(sampled_seed_intervals, weights)
    position_seed = _interval_gram_matrix(sampled_seed_intervals, xweights)

    nseed = length(js)
    odd_block = zeros(Float64, nseed, lseed)
    for k in 1:lseed
        ip = (lseed + 1) + k
        im = (lseed + 1) - k
        odd_block[ip, k] = 1.0
        odd_block[im, k] = -1.0
    end

    even_block = zeros(Float64, nseed, _RADIAL_BOUNDARY_PROTOTYPE_EVEN_K + 1)
    even_block[lseed + 1, 1] = 1.0
    for k in 1:_RADIAL_BOUNDARY_PROTOTYPE_EVEN_K
        ip = (lseed + 1) + k
        im = (lseed + 1) - k
        even_block[ip, k + 1] = 1.0
        even_block[im, k + 1] = 1.0
    end

    function normalize_block(block::Matrix{Float64})
        overlap = overlap_seed * block
        norms = vec(sum(block .* overlap; dims = 1))
        return block * Diagonal(1.0 ./ sqrt.(norms))
    end

    odd_normalized = normalize_block(odd_block)
    odd_overlap = Matrix(Symmetric(transpose(odd_normalized) * overlap_seed * odd_normalized))
    odd_position = Matrix(Symmetric(transpose(odd_normalized) * position_seed * odd_normalized))
    odd_vectors, odd_invhalf = _s_invsqrt_reduced(odd_overlap)
    odd_localizer, _ = _comx_reduced(odd_vectors, odd_invhalf, odd_position)
    odd_coefficients = odd_normalized * (odd_vectors * (odd_invhalf * odd_localizer))

    even_normalized = normalize_block(even_block)
    even_projection = transpose(odd_coefficients) * overlap_seed * even_normalized
    even_clean = even_normalized .- (odd_coefficients * even_projection)
    delta_values = [_shifted_gausslet_value(seed_sampler, j, 0.0) for j in js]
    delta_direction = vec(delta_values' * even_clean)
    unit = delta_direction ./ sqrt(dot(delta_direction, delta_direction))
    sign_value = unit[1] >= 0.0 ? 1.0 : -1.0
    reflector = copy(unit)
    reflector[1] += sign_value
    beta = dot(reflector, reflector)
    householder =
        Matrix{Float64}(I, length(unit), length(unit)) .-
        (2.0 / beta) .* (reflector * transpose(reflector))
    even_clean = even_clean * householder[:, 2:end]

    combined_base = hcat(odd_coefficients, even_clean)
    sampled_combined = _interval_sample_matrix(sampled_seed_intervals, length(xgrid))
    position_combined = position_seed

    sampled_x_intervals = _sample_xgaussian_intervals(xgaussians, xgrid)
    sampled_x = _interval_sample_matrix(sampled_x_intervals, length(xgrid))
    position_sx =
        _interval_cross_gram_matrix(sampled_seed_intervals, sampled_x_intervals, xweights)
    position_xx = _interval_gram_matrix(sampled_x_intervals, xweights)
    position_combined = [position_seed position_sx; position_sx' position_xx]
    sampled_combined = hcat(sampled_combined, sampled_x)
    combined_base = hcat(
        vcat(combined_base, zeros(Float64, length(xgaussians), size(combined_base, 2))),
        vcat(zeros(Float64, nseed, length(xgaussians)), Matrix{Float64}(I, length(xgaussians), length(xgaussians))),
    )

    final_coefficients, reference_centers = _finalize_localized_basis(
        combined_base,
        sampled_combined,
        position_combined,
        weights,
    )
    primitive_data, coefficient_matrix = _radial_primitive_layer(
        js,
        final_coefficients,
        GaussletFamily(:G10),
        spacing,
        xgaussians,
        mapping,
    )
    center_data = _physical_centers(reference_centers, mapping)
    integral_weight_data =
        _integral_weights_from_representation(primitive_data, coefficient_matrix)
    spec = _paper_parity_runtime_spec(mapping)
    return RadialBasis(
        spec,
        primitive_data,
        coefficient_matrix,
        reference_centers,
        center_data,
        integral_weight_data,
        _RADIAL_BOUNDARY_PROTOTYPE_RMAX_INT,
    )
end
