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

"""
    BondAlignedDiatomicQWBasis3D

Narrow mixed-axis basis container for the first bond-aligned diatomic QW
reference route.

The first supported geometry family is a bond-aligned homonuclear diatomic:

- one combined multi-center mapping on the distinguished bond axis
- one shared single-center mapping on the two transverse axes
- one rectangular 3D product basis built from those three one-dimensional
  mapped bases
"""
struct BondAlignedDiatomicQWBasis3D{B<:MappedUniformBasis}
    bond_axis::Symbol
    basis_x::B
    basis_y::B
    basis_z::B
    nuclei::Vector{NTuple{3,Float64}}
    target_core_spacing::Float64
end

function Base.show(io::IO, basis::BondAlignedDiatomicQWBasis3D)
    print(
        io,
        "BondAlignedDiatomicQWBasis3D(bond_axis=:",
        basis.bond_axis,
        ", nx=",
        length(basis.basis_x),
        ", ny=",
        length(basis.basis_y),
        ", nz=",
        length(basis.basis_z),
        ", nuclei=",
        basis.nuclei,
        ", target_core_spacing=",
        basis.target_core_spacing,
        ")",
    )
end

const _QWRG_RESIDUAL_KEEP_ABS_TOL = 1.0e-8
const _QWRG_RESIDUAL_KEEP_REL_TOL = 1.0e-1

function _qwrg_axis_coordinate(
    point::NTuple{3,Float64},
    axis::Symbol,
)
    axis == :x && return point[1]
    axis == :y && return point[2]
    axis == :z && return point[3]
    throw(ArgumentError("axis must be :x, :y, or :z"))
end

function _qwrg_mapped_odd_count_for_extent(
    mapping_value::AbstractCoordinateMapping,
    xmax::Real;
    reference_spacing::Real = 1.0,
)
    xmax_value = Float64(xmax)
    spacing_value = Float64(reference_spacing)
    xmax_value > 0.0 || throw(ArgumentError("mapped extent helper requires xmax > 0"))
    spacing_value > 0.0 || throw(ArgumentError("mapped extent helper requires reference_spacing > 0"))
    uedge = uofx(mapping_value, xmax_value)
    count = 2 * ceil(Int, uedge / spacing_value) + 1
    isodd(count) || throw(ArgumentError("mapped extent helper must produce an odd count"))
    return count
end

function _qwrg_bond_aligned_homonuclear_nuclei(
    bond_length::Real,
    bond_axis::Symbol,
)
    half = 0.5 * Float64(bond_length)
    bond_axis == :x && return [(-half, 0.0, 0.0), (half, 0.0, 0.0)]
    bond_axis == :y && return [(0.0, -half, 0.0), (0.0, half, 0.0)]
    bond_axis == :z && return [(0.0, 0.0, -half), (0.0, 0.0, half)]
    throw(ArgumentError("bond_axis must be :x, :y, or :z"))
end

"""
    bond_aligned_homonuclear_qw_basis(; ...)

Build the first bond-aligned homonuclear diatomic 3D product basis for the
ordinary QW reference line.

The bond axis uses a combined multi-center inverse-sqrt-density mapping, while
the two transverse axes share a single-center inverse-sqrt mapping at the
common transverse projection.
"""
function bond_aligned_homonuclear_qw_basis(;
    family = :G10,
    bond_length::Real,
    core_spacing::Real = 0.5,
    xmax_parallel::Real = 8.0,
    xmax_transverse::Real = 6.0,
    bond_axis::Symbol = :z,
    nuclear_charge::Real = 1.0,
    reference_spacing::Real = 1.0,
    tail_spacing::Real = 10.0,
)
    core_spacing_value = Float64(core_spacing)
    nuclear_charge_value = Float64(nuclear_charge)
    core_spacing_value > 0.0 || throw(ArgumentError("bond_aligned_homonuclear_qw_basis requires core_spacing > 0"))
    nuclear_charge_value > 0.0 || throw(ArgumentError("bond_aligned_homonuclear_qw_basis requires nuclear_charge > 0"))

    nuclei = _qwrg_bond_aligned_homonuclear_nuclei(bond_length, bond_axis)
    parallel_centers = Float64[_qwrg_axis_coordinate(nucleus, bond_axis) for nucleus in nuclei]
    core_range = sqrt(core_spacing_value / nuclear_charge_value)

    # Alg Nested-Diatomic-Map step 4 and 6: use a combined inverse-sqrt map on
    # the bond axis and one shared single-center map on the transverse axes.
    # See docs/src/algorithms/cartesian_nested_diatomic_coordinate_distortion.md.
    parallel_mapping = fit_combined_invsqrt_mapping(
        centers = parallel_centers,
        core_ranges = fill(core_range, length(parallel_centers)),
        target_spacings = fill(core_spacing_value, length(parallel_centers)),
        tail_spacing = tail_spacing,
    )
    transverse_mapping = fit_combined_invsqrt_mapping(
        centers = [0.0],
        core_ranges = [core_range],
        target_spacings = [core_spacing_value],
        tail_spacing = tail_spacing,
    )

    count_parallel = _qwrg_mapped_odd_count_for_extent(
        parallel_mapping,
        xmax_parallel;
        reference_spacing = reference_spacing,
    )
    count_transverse = _qwrg_mapped_odd_count_for_extent(
        transverse_mapping,
        xmax_transverse;
        reference_spacing = reference_spacing,
    )

    parallel_basis = build_basis(MappedUniformBasisSpec(
        family;
        count = count_parallel,
        mapping = parallel_mapping,
        reference_spacing = reference_spacing,
    ))
    transverse_basis = build_basis(MappedUniformBasisSpec(
        family;
        count = count_transverse,
        mapping = transverse_mapping,
        reference_spacing = reference_spacing,
    ))

    basis_x = bond_axis == :x ? parallel_basis : transverse_basis
    basis_y = bond_axis == :y ? parallel_basis : transverse_basis
    basis_z = bond_axis == :z ? parallel_basis : transverse_basis
    return BondAlignedDiatomicQWBasis3D(
        bond_axis,
        basis_x,
        basis_y,
        basis_z,
        nuclei,
        core_spacing_value,
    )
end

function _qwrg_basis_for_axis(
    basis::BondAlignedDiatomicQWBasis3D,
    axis::Symbol,
)
    axis == :x && return basis.basis_x
    axis == :y && return basis.basis_y
    axis == :z && return basis.basis_z
    throw(ArgumentError("axis must be :x, :y, or :z"))
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

function _qwrg_residual_keep_tol(values::AbstractVector{<:Real})
    isempty(values) && return _QWRG_RESIDUAL_KEEP_ABS_TOL
    return max(_QWRG_RESIDUAL_KEEP_ABS_TOL, _QWRG_RESIDUAL_KEEP_REL_TOL * maximum(values))
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

_qwrg_gaussian_exponent(gaussian::Gaussian) = 1.0 / (2.0 * gaussian.width^2)

function _qwrg_atomic_shell_prefactor(exponent::Float64, power::Int)
    power >= 0 || throw(ArgumentError("atomic shell prefactor requires power >= 0"))
    numerator = 2.0^(2 * power + 0.5) * exponent^(power + 0.5)
    denominator = sqrt(pi) * _qwrg_doublefactorial(2 * power - 1)
    return sqrt(numerator / denominator)
end

function _qwrg_doublefactorial(n::Int)
    n <= 0 && return 1.0
    value = 1.0
    for k in n:-2:1
        value *= k
    end
    return value
end

function _qwrg_shifted_gaussian_moment(gamma::Float64, power::Int)
    isodd(power) && return 0.0
    k = power ÷ 2
    return sqrt(pi) * _qwrg_doublefactorial(2 * k - 1) / (2.0^k * gamma^(k + 0.5))
end

function _qwrg_poly_shift_multiply(
    coefficients::Vector{Float64},
    shift::Float64,
    power::Int,
)
    power >= 0 || throw(ArgumentError("polynomial shift multiply requires power >= 0"))
    result = copy(coefficients)
    for _ in 1:power
        next = zeros(Float64, length(result) + 1)
        for degree in eachindex(result)
            value = result[degree]
            next[degree] += shift * value
            next[degree + 1] += value
        end
        result = next
    end
    return result
end

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
    gamma > 0.0 || throw(ArgumentError("atomic shell integral requires positive total Gaussian exponent"))
    weighted_center =
        (alpha_left * center_left + alpha_right * center_right + extra_exponent * extra_center) / gamma
    constant =
        alpha_left * center_left^2 + alpha_right * center_right^2 +
        extra_exponent * extra_center^2 - gamma * weighted_center^2

    polynomial = Float64[1.0]
    polynomial = _qwrg_poly_shift_multiply(polynomial, weighted_center - center_left, power_left)
    polynomial = _qwrg_poly_shift_multiply(polynomial, weighted_center - center_right, power_right)
    xpower > 0 && (polynomial = _qwrg_poly_shift_multiply(polynomial, weighted_center, xpower))

    value = 0.0
    for degree in eachindex(polynomial)
        value += polynomial[degree] * _qwrg_shifted_gaussian_moment(gamma, degree - 1)
    end
    return prefactor_left * prefactor_right * exp(-constant) * value
end

function _qwrg_atomic_derivative_terms(power::Int, exponent::Float64)
    power == 0 && return ((1, -2.0 * exponent),)
    power == 1 && return ((0, 1.0), (2, -2.0 * exponent))
    throw(ArgumentError("atomic shell derivative terms currently support only powers 0 and 1"))
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
    value = 0.0
    for (derived_left_power, left_scale) in _qwrg_atomic_derivative_terms(power_left, alpha_left)
        for (derived_right_power, right_scale) in _qwrg_atomic_derivative_terms(power_right, alpha_right)
            value += 0.5 * left_scale * right_scale * _qwrg_atomic_basic_integral(
                alpha_left,
                center_left,
                derived_left_power,
                prefactor_left,
                alpha_right,
                center_right,
                derived_right_power,
                prefactor_right,
            )
        end
    end
    return value
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

# Alg QW-RG step 4: Define residual Gaussians by orthogonalizing 3D GTOs
# to the full 3D fixed working space. For the active PGDG-mediated route,
# the carried 1D PGDG auxiliary line includes the COMX cleanup step inside the
# bundle, so this fixed block should already be orthonormal to numerical
# precision before the residual construction.
# See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
function _qwrg_residual_space(
    gausslet_overlap::AbstractMatrix{<:Real},
    overlap_ga::AbstractMatrix{<:Real},
    overlap_aa::AbstractMatrix{<:Real},
)
    gausslet_overlap_value = Matrix{Float64}(gausslet_overlap)
    overlap_error = norm(
        gausslet_overlap_value - Matrix{Float64}(I, size(gausslet_overlap_value, 1), size(gausslet_overlap_value, 2)),
        Inf,
    )
    overlap_error <= 1.0e-8 || throw(
        ArgumentError(
            "Qiu-White / QW-PGDG residual construction requires an orthonormal fixed 3D overlap block",
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
    # The retained residual directions should be physically meaningful compared
    # with the case's own residual scale, while true near-null modes are
    # dropped. A pure absolute cutoff was too brittle across count = 9/11/13.
    keep_tol = _qwrg_residual_keep_tol(decomposition.values)
    keep = findall(>(keep_tol), decomposition.values)
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

function _qwrg_residual_center_data(
    raw_overlap::AbstractMatrix{<:Real},
    x_raw::AbstractMatrix{<:Real},
    y_raw::AbstractMatrix{<:Real},
    z_raw::AbstractMatrix{<:Real},
    raw_to_final::AbstractMatrix{<:Real},
    ngausslet::Int,
)
    residual_coefficients = Matrix{Float64}(raw_to_final[:, (ngausslet + 1):end])
    nresidual = size(residual_coefficients, 2)
    centers = zeros(Float64, nresidual, 3)
    norms = zeros(Float64, nresidual)

    overlap_residual = Matrix{Float64}(transpose(residual_coefficients) * raw_overlap * residual_coefficients)
    for index in 1:nresidual
        vector = view(residual_coefficients, :, index)
        norm_value = Float64(dot(vector, raw_overlap * vector))
        norm_value > 1.0e-12 || throw(
            ArgumentError("residual center extraction requires nonzero residual norm"),
        )

        norms[index] = norm_value
        centers[index, 1] = Float64(dot(vector, x_raw * vector) / norm_value)
        centers[index, 2] = Float64(dot(vector, y_raw * vector) / norm_value)
        centers[index, 3] = Float64(dot(vector, z_raw * vector) / norm_value)
    end

    overlap_error = norm(overlap_residual - I, Inf)
    return (
        centers = centers,
        overlap_error = overlap_error,
        residual_coefficients = residual_coefficients,
        norms = norms,
    )
end

function _qwrg_residual_width_data(
    raw_overlap::AbstractMatrix{<:Real},
    x2_raw::AbstractMatrix{<:Real},
    y2_raw::AbstractMatrix{<:Real},
    z2_raw::AbstractMatrix{<:Real},
    center_data,
)
    residual_coefficients = center_data.residual_coefficients
    centers = center_data.centers
    norms = center_data.norms
    nresidual = size(residual_coefficients, 2)
    widths = zeros(Float64, nresidual, 3)

    for index in 1:nresidual
        vector = view(residual_coefficients, :, index)
        norm_value = norms[index]
        x1 = centers[index, 1]
        y1 = centers[index, 2]
        z1 = centers[index, 3]
        x2 = Float64(dot(vector, x2_raw * vector) / norm_value)
        y2 = Float64(dot(vector, y2_raw * vector) / norm_value)
        z2 = Float64(dot(vector, z2_raw * vector) / norm_value)

        varx = x2 - x1^2
        vary = y2 - y1^2
        varz = z2 - z1^2
        min(varx, vary, varz) > 1.0e-12 || throw(
            ArgumentError("MWG residual moment extraction requires positive residual variances"),
        )

        widths[index, 1] = sqrt(2.0 * varx)
        widths[index, 2] = sqrt(2.0 * vary)
        widths[index, 3] = sqrt(2.0 * varz)
    end

    overlap_residual = Matrix{Float64}(transpose(residual_coefficients) * raw_overlap * residual_coefficients)
    overlap_error = norm(overlap_residual - I, Inf)
    return (
        centers = centers,
        widths = widths,
        overlap_error = overlap_error,
    )
end

function _qwrg_residual_moment_data(
    raw_overlap::AbstractMatrix{<:Real},
    x_raw::AbstractMatrix{<:Real},
    x2_raw::AbstractMatrix{<:Real},
    y2_raw::AbstractMatrix{<:Real},
    z2_raw::AbstractMatrix{<:Real},
    center_data,
)
    width_data = _qwrg_residual_width_data(
        raw_overlap,
        x2_raw,
        y2_raw,
        z2_raw,
        center_data,
    )
    return (
        centers = center_data.centers,
        widths = width_data.widths,
        overlap_error = width_data.overlap_error,
    )
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
    center_data = _qwrg_residual_center_data(
        raw_overlap,
        x_raw,
        y_raw,
        z_raw,
        raw_to_final,
        ngausslet,
    )
    return _qwrg_residual_moment_data(
        raw_overlap,
        x_raw,
        x2_raw,
        y_raw,
        y2_raw,
        z_raw,
        z2_raw,
        center_data,
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

function _qwrg_nearest_indices(
    fixed_centers::AbstractMatrix{<:Real},
    residual_centers::AbstractMatrix{<:Real},
)
    size(fixed_centers, 2) == 3 || throw(
        ArgumentError("nearest-center selection for a nested fixed block requires centers as an n×3 matrix"),
    )
    return Int[
        argmin(
            [
                (fixed_centers[fixed, 1] - residual_centers[index, 1])^2 +
                (fixed_centers[fixed, 2] - residual_centers[index, 2])^2 +
                (fixed_centers[fixed, 3] - residual_centers[index, 3])^2 for fixed in axes(fixed_centers, 1)
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
        interaction[ngausslet + residual, 1:ngausslet] .= gausslet_interaction[:, index]
    end
    for i in 1:nresidual, j in i:nresidual
        value = gausslet_interaction[nearest[i], nearest[j]]
        interaction[ngausslet + i, ngausslet + j] = value
        interaction[ngausslet + j, ngausslet + i] = value
    end
    return interaction
end

function _qwrg_interaction_matrix_nearest(
    fixed_interaction::AbstractMatrix{<:Real},
    fixed_centers::AbstractMatrix{<:Real},
    residual_centers::AbstractMatrix{<:Real},
)
    nfixed = size(fixed_interaction, 1)
    nresidual = size(residual_centers, 1)
    interaction = zeros(Float64, nfixed + nresidual, nfixed + nresidual)
    interaction[1:nfixed, 1:nfixed] .= Matrix{Float64}(fixed_interaction)
    nearest = _qwrg_nearest_indices(fixed_centers, residual_centers)
    for residual in 1:nresidual
        index = nearest[residual]
        interaction[1:nfixed, nfixed + residual] .= fixed_interaction[:, index]
        interaction[nfixed + residual, 1:nfixed] .= fixed_interaction[:, index]
    end
    for i in 1:nresidual, j in i:nresidual
        value = fixed_interaction[nearest[i], nearest[j]]
        interaction[nfixed + i, nfixed + j] = value
        interaction[nfixed + j, nfixed + i] = value
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
                analytic_x.pair_aa[term][i, j] *
                analytic_y.pair_aa[term][i, j] *
                analytic_z.pair_aa[term][i, j]
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

function _qwrg_orbital_data(
    fixed_centers::AbstractMatrix{<:Real},
    residual_centers::AbstractMatrix{<:Real},
    residual_widths::AbstractMatrix{<:Real};
    fixed_kind::Symbol = :nested_fixed,
    fixed_label_prefix::AbstractString = "nf",
)
    size(fixed_centers, 2) == 3 || throw(
        ArgumentError("nested fixed-block orbital data requires an n×3 center matrix"),
    )
    orbitals_out = QiuWhiteHybridOrbital3D[]
    for index in axes(fixed_centers, 1)
        push!(
            orbitals_out,
            QiuWhiteHybridOrbital3D(
                index,
                fixed_kind,
                string(fixed_label_prefix, index),
                fixed_centers[index, 1],
                fixed_centers[index, 2],
                fixed_centers[index, 3],
                NaN,
                NaN,
                NaN,
            ),
        )
    end
    base_index = size(fixed_centers, 1)
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

function _qwrg_diatomic_one_body_matrix(
    basis::BondAlignedDiatomicQWBasis3D,
    bundle_x::_MappedOrdinaryGausslet1DBundle,
    bundle_y::_MappedOrdinaryGausslet1DBundle,
    bundle_z::_MappedOrdinaryGausslet1DBundle,
    expansion::CoulombGaussianExpansion,
    nuclear_charges::AbstractVector{<:Real},
)
    length(nuclear_charges) == length(basis.nuclei) || throw(
        ArgumentError("bond-aligned diatomic QW path requires one nuclear charge per nucleus"),
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

    for (nucleus, charge_value_raw) in zip(basis.nuclei, nuclear_charges)
        charge_value = Float64(charge_value_raw)
        charge_value > 0.0 || throw(ArgumentError("bond-aligned diatomic QW path requires positive nuclear charges"))
        matrix .+= _mapped_coulomb_expanded_symmetric_matrix(
            -charge_value .* expansion.coefficients,
            factor_x[nucleus[1]],
            factor_y[nucleus[2]],
            factor_z[nucleus[3]],
        )
    end
    return 0.5 .* (matrix .+ transpose(matrix))
end

"""
    ordinary_cartesian_qiu_white_operators(
        basis::BondAlignedDiatomicQWBasis3D;
        nuclear_charges = fill(1.0, length(basis.nuclei)),
        expansion = coulomb_gaussian_expansion(doacc = false),
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        timing = false,
    )

Build the first bond-aligned diatomic ordinary QW reference Hamiltonian on the
diatomic distortion path.

This first molecular pass is intentionally narrow:

- one bond-aligned homonuclear diatomic basis object
- no nested fixed block yet
- no molecular Gaussian supplement yet
- the final basis is the distorted 3D gausslet product basis itself, so the
  residual-Gaussian sector is empty
"""
function ordinary_cartesian_qiu_white_operators(
    basis::BondAlignedDiatomicQWBasis3D;
    nuclear_charges::AbstractVector{<:Real} = fill(1.0, length(basis.nuclei)),
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    timing::Bool = false,
)
    gausslet_backend == :numerical_reference || throw(
        ArgumentError("bond-aligned diatomic ordinary_cartesian_qiu_white_operators currently supports only gausslet_backend = :numerical_reference"),
    )
    interaction_treatment in (:ggt_nearest, :mwg) || throw(
        ArgumentError("bond-aligned diatomic ordinary_cartesian_qiu_white_operators requires interaction_treatment = :ggt_nearest or :mwg"),
    )
    timing && println("QW-RG timing  note: bond-aligned diatomic path currently has no residual-Gaussian sector")

    bundle_x = _mapped_ordinary_gausslet_1d_bundle(
        basis.basis_x;
        exponents = expansion.exponents,
        center = 0.0,
        backend = gausslet_backend,
    )
    bundle_y = _mapped_ordinary_gausslet_1d_bundle(
        basis.basis_y;
        exponents = expansion.exponents,
        center = 0.0,
        backend = gausslet_backend,
    )
    bundle_z = _mapped_ordinary_gausslet_1d_bundle(
        basis.basis_z;
        exponents = expansion.exponents,
        center = 0.0,
        backend = gausslet_backend,
    )

    overlap = _qwrg_diatomic_overlap_matrix(bundle_x, bundle_y, bundle_z)
    one_body_hamiltonian = _qwrg_diatomic_one_body_matrix(
        basis,
        bundle_x,
        bundle_y,
        bundle_z,
        expansion,
        nuclear_charges,
    )
    interaction_matrix = _qwrg_diatomic_interaction_matrix(
        bundle_x,
        bundle_y,
        bundle_z,
        expansion,
    )

    gausslet_orbitals = _mapped_cartesian_orbitals(
        centers(basis.basis_x),
        centers(basis.basis_y),
        centers(basis.basis_z),
    )
    gausslet_count = length(gausslet_orbitals)
    zero_residual_centers = zeros(Float64, 0, 3)
    zero_residual_widths = zeros(Float64, 0, 3)

    return QiuWhiteResidualGaussianOperators(
        basis,
        nothing,
        gausslet_backend,
        interaction_treatment,
        expansion,
        overlap,
        one_body_hamiltonian,
        interaction_matrix,
        _qwrg_orbital_data(
            gausslet_orbitals,
            zero_residual_centers,
            zero_residual_widths,
        ),
        gausslet_count,
        0,
        Matrix{Float64}(I, gausslet_count, gausslet_count),
        zero_residual_centers,
        zero_residual_widths,
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
    for term in eachindex(expansion.coefficients)
        hamiltonian .-= Float64(Z) * expansion.coefficients[term] .* @view(fixed_block.gaussian_terms[term, :, :])
    end
    return Matrix{Float64}(0.5 .* (hamiltonian .+ transpose(hamiltonian)))
end

function _qwrg_fixed_block_interaction_matrix(
    fixed_block::_NestedFixedBlock3D,
    expansion::CoulombGaussianExpansion,
)
    interaction = zeros(Float64, size(fixed_block.pair_terms, 2), size(fixed_block.pair_terms, 3))
    for term in eachindex(expansion.coefficients)
        interaction .+= expansion.coefficients[term] .* @view(fixed_block.pair_terms[term, :, :])
    end
    return Matrix{Float64}(0.5 .* (interaction .+ transpose(interaction)))
end

function _ordinary_cartesian_qiu_white_operators_atomic_shell_3d(
    basis::MappedUniformBasis,
    gaussian_data::LegacyAtomicGaussianSupplement;
    expansion::CoulombGaussianExpansion,
    Z::Real,
    interaction_treatment::Symbol,
    gausslet_backend::Symbol,
    timing::Bool,
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

    supplement3d = _atomic_cartesian_shell_supplement_3d(gaussian_data)
    gg_blocks = _qwrg_gausslet_1d_blocks(gausslet_bundle)

    start_ns = time_ns()
    blocks = _qwrg_atomic_cartesian_blocks_3d(
        gausslet_bundle,
        supplement3d,
        expansion,
    )
    timing && _qwrg_record_timing!(timing_io, timings, "explicit 3D atomic raw-block assembly", start_ns)

    gausslet_orbitals = _mapped_cartesian_orbitals(gausslet_bundle.pgdg_intermediate.centers)
    gausslet_count = length(gausslet_orbitals)

    start_ns = time_ns()
    gausslet_overlap_3d = zeros(Float64, gausslet_count, gausslet_count)
    _qwrg_fill_product_matrix!(
        gausslet_overlap_3d,
        gg_blocks.overlap_gg,
        gg_blocks.overlap_gg,
        gg_blocks.overlap_gg,
    )
    timing && _qwrg_record_timing!(timing_io, timings, "3D gausslet overlap assembly", start_ns)

    start_ns = time_ns()
    residual_data = _qwrg_residual_space(gausslet_overlap_3d, blocks.overlap_ga, blocks.overlap_aa)
    timing && _qwrg_record_timing!(timing_io, timings, "residual-space construction", start_ns)

    start_ns = time_ns()
    gausslet_one_body = _qwrg_gausslet_one_body_matrix(gg_blocks, expansion; Z = Z)
    one_body_ga = Matrix{Float64}(blocks.kinetic_ga)
    for term in eachindex(expansion.coefficients)
        one_body_ga .-= Float64(Z) * expansion.coefficients[term] .* blocks.factor_ga[term]
    end
    one_body_aa = _qwrg_atomic_cartesian_one_body_aa(blocks, expansion; Z = Z)
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
    x_gg = _qwrg_gausslet_axis_matrix(gg_blocks, :x)
    y_gg = _qwrg_gausslet_axis_matrix(gg_blocks, :y)
    z_gg = _qwrg_gausslet_axis_matrix(gg_blocks, :z)
    x_raw = [x_gg blocks.position_x_ga; transpose(blocks.position_x_ga) blocks.position_x_aa]
    y_raw = [y_gg blocks.position_y_ga; transpose(blocks.position_y_ga) blocks.position_y_aa]
    z_raw = [z_gg blocks.position_z_ga; transpose(blocks.position_z_ga) blocks.position_z_aa]
    center_data = _qwrg_residual_center_data(
        residual_data.raw_overlap,
        x_raw,
        y_raw,
        z_raw,
        residual_data.raw_to_final,
        gausslet_count,
    )
    residual_centers = center_data.centers
    residual_widths = fill(NaN, size(residual_centers, 1), 3)
    center_data.overlap_error <= 1.0e-8 || throw(
        ArgumentError("Qiu-White residual-center extraction requires an orthonormal residual block"),
    )

    if interaction_treatment == :mwg
        gg_x2 = (
            overlap_gg = gg_blocks.overlap_gg,
            position_gg = gg_blocks.position_gg,
            x2_gg = Matrix{Float64}(_x2_matrix(gausslet_bundle.basis)),
        )
        x2_x_gg = _qwrg_gausslet_axis_matrix(gg_x2, :x; squared = true)
        x2_y_gg = _qwrg_gausslet_axis_matrix(gg_x2, :y; squared = true)
        x2_z_gg = _qwrg_gausslet_axis_matrix(gg_x2, :z; squared = true)
        x2_raw = [x2_x_gg blocks.x2_x_ga; transpose(blocks.x2_x_ga) blocks.x2_x_aa]
        y2_raw = [x2_y_gg blocks.x2_y_ga; transpose(blocks.x2_y_ga) blocks.x2_y_aa]
        z2_raw = [x2_z_gg blocks.x2_z_ga; transpose(blocks.x2_z_ga) blocks.x2_z_aa]
        moment_data = _qwrg_residual_moment_data(
            residual_data.raw_overlap,
            x_raw,
            x2_raw,
            y_raw,
            y2_raw,
            z_raw,
            z2_raw,
            center_data,
        )
        residual_centers = moment_data.centers
        residual_widths = moment_data.widths
    end
    timing && _qwrg_record_timing!(timing_io, timings, interaction_treatment == :mwg ? "raw moment-matrix assembly" : "raw center-matrix assembly", start_ns)

    start_ns = time_ns()
    gausslet_interaction = _qwrg_gausslet_interaction_matrix(gg_blocks, expansion)
    timing && _qwrg_record_timing!(timing_io, timings, "3D gausslet interaction assembly", start_ns)

    start_ns = time_ns()
    interaction_matrix = if interaction_treatment == :ggt_nearest
        _qwrg_interaction_matrix_nearest(
            gausslet_interaction,
            gausslet_orbitals,
            residual_centers,
        )
    elseif interaction_treatment == :mwg
        _qwrg_interaction_matrix_mwg(
            gausslet_bundle,
            gausslet_interaction,
            expansion,
            residual_centers,
            residual_widths,
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
            residual_centers,
            residual_widths,
        ),
        gausslet_count,
        size(residual_centers, 1),
        Matrix{Float64}(residual_data.raw_to_final),
        Matrix{Float64}(residual_centers),
        Matrix{Float64}(residual_widths),
    )
end

function _ordinary_cartesian_qiu_white_operators_nested_atomic_shell_3d(
    fixed_block::_NestedFixedBlock3D,
    gaussian_data::LegacyAtomicGaussianSupplement;
    expansion::CoulombGaussianExpansion,
    Z::Real,
    gausslet_backend::Symbol,
    timing::Bool,
)
    timings = Pair{String,Float64}[]
    timing_io = stdout

    start_ns = time_ns()
    gausslet_bundle = _mapped_ordinary_gausslet_1d_bundle(
        fixed_block.parent_basis;
        exponents = expansion.exponents,
        center = 0.0,
        backend = gausslet_backend,
    )
    timing && _qwrg_record_timing!(timing_io, timings, "shared parent 1D bundle", start_ns)

    supplement3d = _atomic_cartesian_shell_supplement_3d(gaussian_data)

    start_ns = time_ns()
    blocks = _qwrg_atomic_cartesian_blocks_3d(
        gausslet_bundle,
        supplement3d,
        expansion,
    )
    timing && _qwrg_record_timing!(timing_io, timings, "explicit 3D atomic parent raw-block assembly", start_ns)

    contraction = fixed_block.coefficient_matrix
    fixed_count = size(fixed_block.overlap, 1)

    start_ns = time_ns()
    overlap_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.overlap_ga)
    residual_data = _qwrg_residual_space(fixed_block.overlap, overlap_fg, blocks.overlap_aa)
    timing && _qwrg_record_timing!(timing_io, timings, "nested residual-space construction", start_ns)

    start_ns = time_ns()
    kinetic_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.kinetic_ga)
    factor_fg = _qwrg_contract_parent_ga_terms(contraction, blocks.factor_ga)
    one_body_fg = Matrix{Float64}(kinetic_fg)
    for term in eachindex(expansion.coefficients)
        one_body_fg .-= Float64(Z) * expansion.coefficients[term] .* factor_fg[term]
    end
    one_body_aa = _qwrg_atomic_cartesian_one_body_aa(blocks, expansion; Z = Z)
    fixed_one_body = _qwrg_fixed_block_one_body_matrix(fixed_block, expansion; Z = Z)
    _, final_one_body = _qwrg_one_body_matrices(
        fixed_one_body,
        one_body_fg,
        one_body_aa,
        residual_data.raw_to_final,
    )
    timing && _qwrg_record_timing!(timing_io, timings, "nested raw one-body assembly and transform", start_ns)

    start_ns = time_ns()
    x_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.position_x_ga)
    y_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.position_y_ga)
    z_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.position_z_ga)
    x_raw = [Matrix{Float64}(fixed_block.position_x) x_fg; transpose(x_fg) blocks.position_x_aa]
    y_raw = [Matrix{Float64}(fixed_block.position_y) y_fg; transpose(y_fg) blocks.position_y_aa]
    z_raw = [Matrix{Float64}(fixed_block.position_z) z_fg; transpose(z_fg) blocks.position_z_aa]
    center_data = _qwrg_residual_center_data(
        residual_data.raw_overlap,
        x_raw,
        y_raw,
        z_raw,
        residual_data.raw_to_final,
        fixed_count,
    )
    residual_centers = center_data.centers
    residual_widths = fill(NaN, size(residual_centers, 1), 3)
    center_data.overlap_error <= 1.0e-8 || throw(
        ArgumentError("nested QW-PGDG residual-center extraction requires an orthonormal residual block"),
    )
    timing && _qwrg_record_timing!(timing_io, timings, "nested raw center-matrix assembly", start_ns)

    start_ns = time_ns()
    fixed_interaction = _qwrg_fixed_block_interaction_matrix(fixed_block, expansion)
    interaction_matrix = _qwrg_interaction_matrix_nearest(
        fixed_interaction,
        fixed_block.fixed_centers,
        residual_centers,
    )
    timing && _qwrg_record_timing!(timing_io, timings, "nested RG interaction assembly", start_ns)
    timing && _qwrg_maybe_print_timings(timing_io, timings)

    return QiuWhiteResidualGaussianOperators(
        fixed_block,
        gaussian_data,
        gausslet_backend,
        :ggt_nearest,
        expansion,
        Matrix{Float64}(residual_data.final_overlap),
        final_one_body,
        Matrix{Float64}(0.5 .* (interaction_matrix .+ transpose(interaction_matrix))),
        _qwrg_orbital_data(
            fixed_block.fixed_centers,
            residual_centers,
            residual_widths;
            fixed_kind = :nested_fixed,
            fixed_label_prefix = "nf",
        ),
        fixed_count,
        size(residual_centers, 1),
        Matrix{Float64}(residual_data.raw_to_final),
        Matrix{Float64}(residual_centers),
        Matrix{Float64}(residual_widths),
    )
end

"""
    ordinary_cartesian_qiu_white_operators(
        basis::MappedUniformBasis,
        gaussian_data::LegacyAtomicGaussianSupplement;
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

This path now has two active supplement modes:

- the earlier centered separable `s` route for `lmax = 0`
- an explicit atomic-centered 3D Cartesian shell route for `s` and `p`
  (`lmax <= 1`)
"""
function ordinary_cartesian_qiu_white_operators(
    basis::MappedUniformBasis,
    gaussian_data::LegacyAtomicGaussianSupplement;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
    interaction_treatment::Symbol = :mwg,
    gausslet_backend::Symbol = :numerical_reference,
    timing::Bool = false,
    )
    gausslet_backend == :numerical_reference || throw(
        ArgumentError("ordinary_cartesian_qiu_white_operators currently supports only gausslet_backend = :numerical_reference"),
    )
    if _legacy_atomic_has_nonseparable_shells(gaussian_data)
        return _ordinary_cartesian_qiu_white_operators_atomic_shell_3d(
            basis,
            gaussian_data;
            expansion = expansion,
            Z = Z,
            interaction_treatment = interaction_treatment,
            gausslet_backend = gausslet_backend,
            timing = timing,
        )
    end
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
        gaussian_data,
        expansion,
    )
    timing && _qwrg_record_timing!(timing_io, timings, "split gausslet-Gaussian raw-block assembly", start_ns)

    gausslet_orbitals = _mapped_cartesian_orbitals(gausslet_bundle.pgdg_intermediate.centers)
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

    x_raw = [x_gg x_ga; transpose(x_ga) x_aa]
    y_raw = [y_gg y_ga; transpose(y_ga) y_aa]
    z_raw = [z_gg z_ga; transpose(z_ga) z_aa]
    center_data = _qwrg_residual_center_data(
        residual_data.raw_overlap,
        x_raw,
        y_raw,
        z_raw,
        residual_data.raw_to_final,
        gausslet_count,
    )

    residual_centers = center_data.centers
    residual_widths = fill(NaN, size(residual_centers, 1), 3)
    center_data.overlap_error <= 1.0e-8 || throw(
        ArgumentError("Qiu-White residual-center extraction requires an orthonormal residual block"),
    )

    if interaction_treatment == :mwg
        x2_gg_exact = Matrix{Float64}(_x2_matrix(gausslet_bundle.basis))
        x2_gg, x2_ga, x2_aa = _qwrg_raw_axis_blocks((
            overlap_gg = blocks.overlap_gg,
            overlap_ga = blocks.overlap_ga,
            overlap_aa = blocks.overlap_aa,
            position_gg = blocks.position_gg,
            position_ga = blocks.position_ga,
            position_aa = blocks.position_aa,
            x2_gg = x2_gg_exact,
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
            x2_gg = x2_gg_exact,
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
            x2_gg = x2_gg_exact,
            x2_ga = blocks.x2_ga,
            x2_aa = blocks.x2_aa,
        ), :z; squared = true)

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
            center_data,
        )
        residual_centers = moment_data.centers
        residual_widths = moment_data.widths
    end
    timing && _qwrg_record_timing!(timing_io, timings, interaction_treatment == :mwg ? "raw moment-matrix assembly" : "raw center-matrix assembly", start_ns)

    start_ns = time_ns()
    gausslet_interaction = _qwrg_gausslet_interaction_matrix(blocks, expansion)
    timing && _qwrg_record_timing!(timing_io, timings, "3D gausslet interaction assembly", start_ns)

    start_ns = time_ns()
    interaction_matrix = if interaction_treatment == :ggt_nearest
        _qwrg_interaction_matrix_nearest(
            gausslet_interaction,
            gausslet_orbitals,
            residual_centers,
        )
    elseif interaction_treatment == :mwg
        _qwrg_interaction_matrix_mwg(
            gausslet_bundle,
            gausslet_interaction,
            expansion,
            residual_centers,
            residual_widths,
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
            residual_centers,
            residual_widths,
        ),
        gausslet_count,
        size(residual_centers, 1),
        Matrix{Float64}(residual_data.raw_to_final),
        Matrix{Float64}(residual_centers),
        Matrix{Float64}(residual_widths),
    )
end

"""
    ordinary_cartesian_qiu_white_operators(
        fixed_block::_NestedFixedBlock3D,
        gaussian_data::LegacyAtomicGaussianSupplement;
        expansion = coulomb_gaussian_expansion(doacc = false),
        Z = 2.0,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        timing = false,
    )

Build the first nested fixed-block QW-PGDG consumer path.

This overload reuses the stabilized parent PGDG raw fixed-to-Gaussian blocks,
contracts them through the supplied shell-level fixed map, and then runs the
same downstream residual-space / one-body / nearest-GGT algebra as the
unnested QW-PGDG route.

This first adapter is intentionally narrow:

- it consumes an already-assembled nonseparable 3D fixed packet
- it keeps the shell packet as the fixed-fixed block directly
- it supports only `interaction_treatment = :ggt_nearest`

It also shares the same supplement split as the ordinary path:

- the earlier centered separable `s` route for `lmax = 0`
- an explicit atomic-centered 3D Cartesian shell route for `s` and `p`
  (`lmax <= 1`)
"""
function ordinary_cartesian_qiu_white_operators(
    fixed_block::_NestedFixedBlock3D,
    gaussian_data::LegacyAtomicGaussianSupplement;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    timing::Bool = false,
)
    gausslet_backend == :numerical_reference || throw(
        ArgumentError("nested ordinary_cartesian_qiu_white_operators currently supports only gausslet_backend = :numerical_reference"),
    )
    interaction_treatment == :ggt_nearest || throw(
        ArgumentError("nested ordinary_cartesian_qiu_white_operators currently supports only interaction_treatment = :ggt_nearest"),
    )
    if _legacy_atomic_has_nonseparable_shells(gaussian_data)
        return _ordinary_cartesian_qiu_white_operators_nested_atomic_shell_3d(
            fixed_block,
            gaussian_data;
            expansion = expansion,
            Z = Z,
            gausslet_backend = gausslet_backend,
            timing = timing,
        )
    end

    timings = Pair{String,Float64}[]
    timing_io = stdout

    start_ns = time_ns()
    gausslet_bundle = _mapped_ordinary_gausslet_1d_bundle(
        fixed_block.parent_basis;
        exponents = expansion.exponents,
        center = 0.0,
        backend = gausslet_backend,
    )
    timing && _qwrg_record_timing!(timing_io, timings, "shared parent 1D bundle", start_ns)

    start_ns = time_ns()
    blocks = _qwrg_split_block_matrices(
        gausslet_bundle,
        gaussian_data,
        expansion,
    )
    timing && _qwrg_record_timing!(timing_io, timings, "parent split raw-block assembly", start_ns)

    contraction = fixed_block.coefficient_matrix
    fixed_count = size(fixed_block.overlap, 1)

    start_ns = time_ns()
    overlap_parent_ga, overlap_aa = _qwrg_raw_overlap_blocks((
        overlap_gg = blocks.overlap_gg,
        overlap_ga = blocks.overlap_ga,
        overlap_aa = blocks.overlap_aa,
    ))
    overlap_fg = _qwrg_contract_parent_ga_matrix(contraction, overlap_parent_ga)
    residual_data = _qwrg_residual_space(fixed_block.overlap, overlap_fg, overlap_aa)
    timing && _qwrg_record_timing!(timing_io, timings, "nested residual-space construction", start_ns)

    start_ns = time_ns()
    kinetic_parent_ga = _qwrg_raw_kinetic_cross_block((
        overlap_ga = blocks.overlap_ga,
        kinetic_ga = blocks.kinetic_ga,
    ))
    factor_parent_ga = _qwrg_raw_factor_cross_blocks((
        overlap_ga = blocks.overlap_ga,
        factor_ga = blocks.factor_ga,
    ), expansion)
    kinetic_fg = _qwrg_contract_parent_ga_matrix(contraction, kinetic_parent_ga)
    factor_fg = _qwrg_contract_parent_ga_terms(contraction, factor_parent_ga)
    one_body_fg = Matrix{Float64}(kinetic_fg)
    for term in eachindex(expansion.coefficients)
        one_body_fg .-= Float64(Z) * expansion.coefficients[term] .* factor_fg[term]
    end
    one_body_aa = _qwrg_raw_one_body_aa_block((
        overlap_aa = blocks.overlap_aa,
        kinetic_aa = blocks.kinetic_aa,
        factor_aa = blocks.factor_aa,
    ), expansion; Z = Z)
    fixed_one_body = _qwrg_fixed_block_one_body_matrix(fixed_block, expansion; Z = Z)
    _, final_one_body = _qwrg_one_body_matrices(
        fixed_one_body,
        one_body_fg,
        one_body_aa,
        residual_data.raw_to_final,
    )
    timing && _qwrg_record_timing!(timing_io, timings, "nested raw one-body assembly and transform", start_ns)

    start_ns = time_ns()
    x_parent_gg, x_parent_ga, x_aa = _qwrg_raw_axis_blocks((
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
    y_parent_gg, y_parent_ga, y_aa = _qwrg_raw_axis_blocks((
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
    z_parent_gg, z_parent_ga, z_aa = _qwrg_raw_axis_blocks((
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
    x_fg = _qwrg_contract_parent_ga_matrix(contraction, x_parent_ga)
    y_fg = _qwrg_contract_parent_ga_matrix(contraction, y_parent_ga)
    z_fg = _qwrg_contract_parent_ga_matrix(contraction, z_parent_ga)

    x2_parent_gg, x2_parent_ga, x2_aa = _qwrg_raw_axis_blocks((
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
    y2_parent_gg, y2_parent_ga, y2_aa = _qwrg_raw_axis_blocks((
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
    z2_parent_gg, z2_parent_ga, z2_aa = _qwrg_raw_axis_blocks((
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
    x2_fg = _qwrg_contract_parent_ga_matrix(contraction, x2_parent_ga)
    y2_fg = _qwrg_contract_parent_ga_matrix(contraction, y2_parent_ga)
    z2_fg = _qwrg_contract_parent_ga_matrix(contraction, z2_parent_ga)

    x_raw = [Matrix{Float64}(fixed_block.position_x) x_fg; transpose(x_fg) x_aa]
    y_raw = [Matrix{Float64}(fixed_block.position_y) y_fg; transpose(y_fg) y_aa]
    z_raw = [Matrix{Float64}(fixed_block.position_z) z_fg; transpose(z_fg) z_aa]
    center_data = _qwrg_residual_center_data(
        residual_data.raw_overlap,
        x_raw,
        y_raw,
        z_raw,
        residual_data.raw_to_final,
        fixed_count,
    )
    residual_centers = center_data.centers
    residual_widths = fill(NaN, size(residual_centers, 1), 3)
    center_data.overlap_error <= 1.0e-8 || throw(
        ArgumentError("nested QW-PGDG residual-center extraction requires an orthonormal residual block"),
    )
    timing && _qwrg_record_timing!(timing_io, timings, "nested raw center-matrix assembly", start_ns)

    start_ns = time_ns()
    fixed_interaction = _qwrg_fixed_block_interaction_matrix(fixed_block, expansion)
    interaction_matrix = _qwrg_interaction_matrix_nearest(
        fixed_interaction,
        fixed_block.fixed_centers,
        residual_centers,
    )
    timing && _qwrg_record_timing!(timing_io, timings, "nested RG interaction assembly", start_ns)
    timing && _qwrg_maybe_print_timings(timing_io, timings)

    return QiuWhiteResidualGaussianOperators(
        fixed_block,
        gaussian_data,
        gausslet_backend,
        interaction_treatment,
        expansion,
        Matrix{Float64}(residual_data.final_overlap),
        final_one_body,
        Matrix{Float64}(0.5 .* (interaction_matrix .+ transpose(interaction_matrix))),
        _qwrg_orbital_data(
            fixed_block.fixed_centers,
            residual_centers,
            residual_widths;
            fixed_kind = :nested_fixed,
            fixed_label_prefix = "nf",
        ),
        fixed_count,
        size(residual_centers, 1),
        Matrix{Float64}(residual_data.raw_to_final),
        Matrix{Float64}(residual_centers),
        Matrix{Float64}(residual_widths),
    )
end
