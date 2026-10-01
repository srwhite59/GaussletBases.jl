function _same_mapping(::IdentityMapping, ::IdentityMapping)
    return true
end

function _same_mapping(a::AsinhMapping, b::AsinhMapping)
    return a.a == b.a && a.s == b.s && a.tail_spacing == b.tail_spacing
end

function _same_mapping(a::AbstractCoordinateMapping, b::AbstractCoordinateMapping)
    return a === b
end

function _validate_radial_operator_grid(basis::RadialBasis, grid::RadialQuadratureGrid)
    points = quadrature_points(grid)
    weights = quadrature_weights(grid)
    length(points) == length(weights) || throw(ArgumentError("quadrature point and weight counts must match"))
    isempty(points) && throw(ArgumentError("quadrature grid must not be empty"))
    issorted(points) || throw(ArgumentError("quadrature points must be sorted in increasing order"))
    any(weight -> weight <= 0.0, weights) && throw(ArgumentError("quadrature weights must be positive"))

    if grid.mapping_value !== nothing && !_same_mapping(mapping(basis), grid.mapping_value)
        throw(ArgumentError("quadrature grid mapping is incompatible with the supplied radial basis"))
    end
    return points, weights
end

function _validate_radial_operator_grid(::PrimitiveSet1D, grid::RadialQuadratureGrid)
    points = quadrature_points(grid)
    weights = quadrature_weights(grid)
    length(points) == length(weights) || throw(ArgumentError("quadrature point and weight counts must match"))
    isempty(points) && throw(ArgumentError("quadrature grid must not be empty"))
    issorted(points) || throw(ArgumentError("quadrature points must be sorted in increasing order"))
    any(weight -> weight <= 0.0, weights) && throw(ArgumentError("quadrature weights must be positive"))
    return points, weights
end

function _symmetrize_matrix(matrix::AbstractMatrix{<:Real})
    return 0.5 .* (Matrix{Float64}(matrix) .+ Matrix{Float64}(transpose(matrix)))
end

function _primitive_derivative_sample_matrix(
    primitive_data::Vector{AbstractPrimitiveFunction1D},
    points::AbstractVector{Float64};
    order::Int = 1,
)
    samples = zeros(Float64, length(points), length(primitive_data))
    for mu in eachindex(primitive_data)
        samples[:, mu] = [derivative(primitive_data[mu], point; order = order) for point in points]
    end
    return samples
end

function _basis_derivative_matrix(
    basis::RadialBasis,
    points::AbstractVector{Float64};
    order::Int = 1,
)
    primitive_derivatives = _primitive_derivative_sample_matrix(primitives(basis), points; order = order)
    return primitive_derivatives * stencil_matrix(basis)
end

function _weighted_basis_gram(
    left_values::AbstractMatrix{<:Real},
    right_values::AbstractMatrix{<:Real},
    weights::AbstractVector{Float64},
)
    return Matrix(transpose(left_values) * (weights .* right_values))
end

function _radial_basis_integral_weights(
    values::AbstractMatrix{<:Real},
    weights::AbstractVector{Float64},
)
    return vec(transpose(values) * weights)
end

function _check_integral_weights(weight_data::AbstractVector{Float64})
    any(weight -> !isfinite(weight), weight_data) &&
        throw(ArgumentError("basis integral weights on the supplied grid must be finite"))
    any(weight -> abs(weight) <= sqrt(eps(Float64)), weight_data) &&
        throw(ArgumentError("basis integral weights on the supplied grid are too small for IntegralDiagonal"))
    return weight_data
end

"""
    IntegralDiagonal()

Two-index integral-diagonal approximation for the radial Coulomb multipole
matrices.

For `multipole_matrix(...; approximation = IntegralDiagonal())`, the returned
matrix has entries

    V_ab^(L) = [integral integral chi_a(r) K^(L)(r, r') chi_b(r') dr dr'] /
               [integral chi_a(r) dr * integral chi_b(r) dr]

with `K^(L)(r, r') = r_<^L / r_>^(L + 1)`.
"""
struct IntegralDiagonal <: AbstractDiagonalApproximation
end

"""
    overlap_matrix(basis::RadialBasis, grid::RadialQuadratureGrid)

Build the radial overlap matrix of `basis` on the supplied quadrature `grid`.

The quadrature grid is used directly. No hidden grid is constructed inside this
builder.
"""
function overlap_matrix(basis::RadialBasis, grid::RadialQuadratureGrid)
    points, weights = _validate_radial_operator_grid(basis, grid)
    values = _basis_values_matrix(basis, points)
    return _symmetrize_matrix(_weighted_basis_gram(values, values, weights))
end

"""
    overlap_matrix(set::PrimitiveSet1D, grid::RadialQuadratureGrid)

Build the primitive-space radial overlap matrix of `set` on the supplied
quadrature `grid`.

This is intended for the primitive layer behind a `RadialBasis`, for example
through `primitive_set(rb)`.
"""
function overlap_matrix(set::PrimitiveSet1D, grid::RadialQuadratureGrid)
    points, weights = _validate_radial_operator_grid(set, grid)
    values = _primitive_sample_matrix(set, points)
    return _symmetrize_matrix(_weighted_basis_gram(values, values, weights))
end

"""
    kinetic_matrix(basis::RadialBasis, grid::RadialQuadratureGrid)

Build the reduced-radial kinetic-energy matrix

    <chi_a | -0.5 d^2/dr^2 | chi_b>

on the supplied quadrature `grid`.
"""
function kinetic_matrix(basis::RadialBasis, grid::RadialQuadratureGrid)
    points, weights = _validate_radial_operator_grid(basis, grid)
    derivatives = _basis_derivative_matrix(basis, points)
    return _symmetrize_matrix(0.5 .* _weighted_basis_gram(derivatives, derivatives, weights))
end

"""
    kinetic_matrix(set::PrimitiveSet1D, grid::RadialQuadratureGrid)

Build the primitive-space reduced-radial kinetic-energy matrix

    <phi_mu | -0.5 d^2/dr^2 | phi_nu>

on the supplied quadrature `grid`.
"""
function kinetic_matrix(set::PrimitiveSet1D, grid::RadialQuadratureGrid)
    points, weights = _validate_radial_operator_grid(set, grid)
    derivatives = _primitive_sample_matrix(set, points; derivative_order = 1)
    return _symmetrize_matrix(0.5 .* _weighted_basis_gram(derivatives, derivatives, weights))
end

"""
    nuclear_matrix(basis::RadialBasis, grid::RadialQuadratureGrid; Z)

Build the reduced-radial nuclear attraction matrix

    <chi_a | -Z / r | chi_b>

on the supplied quadrature `grid`.
"""
function nuclear_matrix(basis::RadialBasis, grid::RadialQuadratureGrid; Z::Real)
    points, weights = _validate_radial_operator_grid(basis, grid)
    any(point -> point <= 0.0, points) && throw(ArgumentError("nuclear_matrix requires quadrature points strictly above zero"))
    values = _basis_values_matrix(basis, points)
    radial_factor = (-Float64(Z)) ./ points
    return _symmetrize_matrix(_weighted_basis_gram(values, values, weights .* radial_factor))
end

"""
    nuclear_matrix(set::PrimitiveSet1D, grid::RadialQuadratureGrid; Z)

Build the primitive-space reduced-radial nuclear attraction matrix

    <phi_mu | -Z / r | phi_nu>

on the supplied quadrature `grid`.
"""
function nuclear_matrix(set::PrimitiveSet1D, grid::RadialQuadratureGrid; Z::Real)
    points, weights = _validate_radial_operator_grid(set, grid)
    any(point -> point <= 0.0, points) && throw(ArgumentError("nuclear_matrix requires quadrature points strictly above zero"))
    values = _primitive_sample_matrix(set, points)
    radial_factor = (-Float64(Z)) ./ points
    return _symmetrize_matrix(_weighted_basis_gram(values, values, weights .* radial_factor))
end

"""
    centrifugal_matrix(basis::RadialBasis, grid::RadialQuadratureGrid; l)

Build the reduced-radial centrifugal matrix

    <chi_a | l(l + 1) / (2 r^2) | chi_b>

on the supplied quadrature `grid`.
"""
function centrifugal_matrix(basis::RadialBasis, grid::RadialQuadratureGrid; l::Int)
    l >= 0 || throw(ArgumentError("centrifugal_matrix requires l >= 0"))
    points, weights = _validate_radial_operator_grid(basis, grid)
    any(point -> point <= 0.0, points) && throw(ArgumentError("centrifugal_matrix requires quadrature points strictly above zero"))

    nbasis = length(basis)
    l == 0 && return zeros(Float64, nbasis, nbasis)

    values = _basis_values_matrix(basis, points)
    radial_factor = (0.5 * l * (l + 1.0)) ./ (points .^ 2)
    return _symmetrize_matrix(_weighted_basis_gram(values, values, weights .* radial_factor))
end

"""
    centrifugal_matrix(set::PrimitiveSet1D, grid::RadialQuadratureGrid; l)

Build the primitive-space reduced-radial centrifugal matrix

    <phi_mu | l(l + 1) / (2 r^2) | phi_nu>

on the supplied quadrature `grid`.
"""
function centrifugal_matrix(set::PrimitiveSet1D, grid::RadialQuadratureGrid; l::Int)
    l >= 0 || throw(ArgumentError("centrifugal_matrix requires l >= 0"))
    points, weights = _validate_radial_operator_grid(set, grid)
    any(point -> point <= 0.0, points) && throw(ArgumentError("centrifugal_matrix requires quadrature points strictly above zero"))

    nprimitive = length(set)
    l == 0 && return zeros(Float64, nprimitive, nprimitive)

    values = _primitive_sample_matrix(set, points)
    radial_factor = (0.5 * l * (l + 1.0)) ./ (points .^ 2)
    return _symmetrize_matrix(_weighted_basis_gram(values, values, weights .* radial_factor))
end

function _integral_diagonal_kernel_matrix_raw(
    values::AbstractMatrix{<:Real},
    points::AbstractVector{Float64},
    weights::AbstractVector{Float64},
    L::Int,
)
    any(point -> point <= 0.0, points) && throw(ArgumentError("multipole_matrix requires quadrature points strictly above zero"))

    rpow = L == 0 ? ones(Float64, length(points)) : points .^ L
    invrpow = 1.0 ./ (points .^ (L + 1))

    weighted_prefix = (weights .* rpow) .* values
    prefix = cumsum(weighted_prefix; dims = 1)

    weighted_suffix = (weights .* invrpow) .* values
    suffix_inclusive = reverse(cumsum(reverse(weighted_suffix; dims = 1); dims = 1); dims = 1)
    suffix_exclusive = suffix_inclusive .- weighted_suffix

    inner = (invrpow .* prefix) .+ (rpow .* suffix_exclusive)
    numerator = _weighted_basis_gram(values, inner, weights)

    integral_weights = _check_integral_weights(_radial_basis_integral_weights(values, weights))
    inv_integral_weights = 1.0 ./ integral_weights
    return _symmetrize_matrix(Diagonal(inv_integral_weights) * numerator * Diagonal(inv_integral_weights))
end

# Adjacent-point ratio powers rho[p] = (r[p-1] / r[p])^L in (0, 1] for the running-normalized
# multipole recursion. log1p of the relative step keeps the per-step relative error at a few ulps
# even for L in the hundreds; rho[1] is unused and set to 1.
function _radial_multipole_adjacent_ratio_powers(points::AbstractVector{Float64}, L::Int)
    L >= 0 || throw(ArgumentError("multipole_matrix requires L >= 0"))
    rho = ones(Float64, length(points))
    L == 0 && return rho
    @inbounds for p in 2:length(points)
        rho[p] = exp(L * log1p((points[p - 1] - points[p]) / points[p]))
    end
    return rho
end

# Integral-diagonal radial multipole kernel
#
#     R_L(a, b) = [sum_p sum_q chi_a(r_p) W_p W_q chi_b(r_q) r_<^L / r_>^(L + 1)] / (w_a w_b)
#
# evaluated with a running normalization instead of global r^L / r^-(L+1) scale factors. With
# x_p = W_p chi(r_p) and rho_p = (r_{p-1} / r_p)^L the prefix and suffix sums are
#
#     A_p = rho_p A_{p-1} + x_p                    = sum_{q <= p} x_q (r_q / r_p)^L
#     B_p = rho_{p+1} (B_{p+1} + x_{p+1} / r_{p+1}) = sum_{q > p} x_q (r_p / r_q)^L / r_q
#
# and inner_p = A_p / r_p + B_p. Local ratios avoid artificial global-scale loss, but do not
# remove genuine Float64 range or cancellation limitations. (The previous form scaled
# r^L by r_max^L and r^-(L+1) by
# r_min^-(L+1); their product underflowed for L * log(r_max / r_min) beyond ~700, i.e. L >= 17 on
# production atomic grids, and silently returned zeros.)
function _integral_diagonal_kernel_matrix(
    values::AbstractMatrix{<:Real},
    points::AbstractVector{Float64},
    weights::AbstractVector{Float64},
    L::Int,
)
    all(point -> isfinite(point) && point > 0.0, points) ||
        throw(ArgumentError("multipole_matrix requires finite quadrature points strictly above zero"))
    issorted(points) || throw(ArgumentError("multipole_matrix requires sorted quadrature points"))
    npoints, nbasis = size(values)
    length(points) == npoints == length(weights) ||
        throw(DimensionMismatch("multipole_matrix requires one basis-value row per quadrature point"))

    rho = _radial_multipole_adjacent_ratio_powers(points, L)
    weighted_values = Matrix{Float64}(weights .* values)
    inner = Matrix{Float64}(undef, npoints, nbasis)
    @inbounds for j in 1:nbasis
        prefix = 0.0
        for p in 1:npoints
            prefix = rho[p] * prefix + weighted_values[p, j]
            inner[p, j] = prefix / points[p]
        end
        suffix = 0.0
        for p in (npoints - 1):-1:1
            suffix = rho[p + 1] * (suffix + weighted_values[p + 1, j] / points[p + 1])
            inner[p, j] += suffix
        end
    end
    all(isfinite, inner) || throw(
        ArgumentError(
            "multipole_matrix produced non-finite inner data for L=$(L) over radial range [$(minimum(points)), $(maximum(points))]",
        ),
    )

    numerator = _weighted_basis_gram(values, inner, weights)
    integral_weights = _check_integral_weights(_radial_basis_integral_weights(values, weights))
    inv_integral_weights = 1.0 ./ integral_weights
    return _symmetrize_matrix(Diagonal(inv_integral_weights) * numerator * Diagonal(inv_integral_weights))
end

"""
    multipole_matrix(basis::RadialBasis, grid::RadialQuadratureGrid;
                     L::Int,
                     approximation::AbstractDiagonalApproximation = IntegralDiagonal())

Build the v0 supported two-index radial multipole matrix on the supplied
quadrature `grid`.

In this release, `multipole_matrix` means the IDA-style two-index radial object
used together with separate angular Gaunt/Ylm machinery. It does not represent
an exact four-index electron-electron tensor.
"""
function multipole_matrix(
    basis::RadialBasis,
    grid::RadialQuadratureGrid;
    L::Int,
    approximation::AbstractDiagonalApproximation = IntegralDiagonal(),
)
    L >= 0 || throw(ArgumentError("multipole_matrix requires L >= 0"))
    return _radial_multipole_from_samples(_radial_multipole_samples(basis, grid), L, approximation)
end

# Quadrature samples behind the radial multipole kernel. `atomic_operators` keeps them so that
# multipoles beyond the stored range can be evaluated later with exactly the same kernel.
struct _RadialMultipoleSamples
    points::Vector{Float64}
    weights::Vector{Float64}
    values::Matrix{Float64}
end

function _radial_multipole_samples(basis::RadialBasis, grid::RadialQuadratureGrid)
    points, weights = _validate_radial_operator_grid(basis, grid)
    values = _basis_values_matrix(basis, points)
    return _RadialMultipoleSamples(
        Vector{Float64}(points),
        Vector{Float64}(weights),
        Matrix{Float64}(values),
    )
end

function _radial_multipole_from_samples(
    samples::_RadialMultipoleSamples,
    L::Int,
    approximation::AbstractDiagonalApproximation,
)
    L >= 0 || throw(ArgumentError("multipole_matrix requires L >= 0"))
    if approximation isa IntegralDiagonal
        return _integral_diagonal_kernel_matrix(samples.values, samples.points, samples.weights, L)
    end
    throw(ArgumentError("unsupported diagonal approximation $(typeof(approximation))"))
end

struct _AtomicRadialSourceManifest{S <: RadialBasisSpec}
    basis_spec::S
    nuclear_charge::Float64
end

"""
    RadialAtomicOperators

Bundle of radial one-body matrices and precomputed radial multipole tables built
for a `RadialBasis` on an explicit `RadialQuadratureGrid`.

Use `ops.overlap`, `ops.kinetic`, and `ops.nuclear` for the fixed one-body
matrices, and `centrifugal(ops, l)` / `multipole(ops, L)` for the indexed
families.
"""
struct RadialAtomicOperators{A <: AbstractDiagonalApproximation,S <: _AtomicRadialSourceManifest}
    overlap::Matrix{Float64}
    kinetic::Matrix{Float64}
    nuclear::Matrix{Float64}
    centrifugal_data::Vector{Matrix{Float64}}
    multipole_data::Vector{Matrix{Float64}}
    shell_centers_r::Vector{Float64}
    source_manifest::S
    approximation::A
    multipole_samples::Union{Nothing,_RadialMultipoleSamples}
end

# Pre-2026-10 positional layout: no retained quadrature samples, so multipoles cannot be
# extended beyond `multipole_data`.
RadialAtomicOperators(args::Vararg{Any,8}) = RadialAtomicOperators(args..., nothing)
RadialAtomicOperators{A,S}(args::Vararg{Any,8}) where {A,S} =
    RadialAtomicOperators{A,S}(args..., nothing)

function Base.show(io::IO, ops::RadialAtomicOperators)
    print(
        io,
        "RadialAtomicOperators(size=",
        size(ops.overlap),
        ", lmax=",
        length(ops.centrifugal_data) - 1,
        ", Lmax=",
        length(ops.multipole_data) - 1,
        ", nradial=",
        length(ops.shell_centers_r),
        ", approximation=",
    )
    show(io, ops.approximation)
    print(io, ")")
end

"""
    atomic_operators(basis::RadialBasis, grid::RadialQuadratureGrid;
                     Z,
                     lmax::Int = 0,
                     multipole_lmax::Int = 2 * lmax,
                     approximation::AbstractDiagonalApproximation = IntegralDiagonal())

Build the high-level radial operator bundle for `basis` on the supplied
quadrature `grid`.

The bundle stores:
- `ops.overlap`
- `ops.kinetic`
- `ops.nuclear`
- `centrifugal(ops, l)` for `l = 0:lmax`
- `multipole(ops, L)` for `L = 0:multipole_lmax`

`lmax` is the one-electron (centrifugal) angular-momentum range. The default
`multipole_lmax = 2 * lmax` is the product range needed by `Y_lm` channels with
`l <= lmax` (`atomic_ida_operators`). Shell-local angular interactions usually
need far more multipoles than that (up to the shell interaction moment `lcap`,
for example 46 for `NΩ = 98`). The bundle therefore also keeps the quadrature
samples, and the angular interaction builders evaluate any missing multipoles
on demand with the same kernel. Pass a larger `multipole_lmax` to store them
up front.
"""
function atomic_operators(
    basis::RadialBasis,
    grid::RadialQuadratureGrid;
    Z::Real,
    lmax::Int = 0,
    multipole_lmax::Int = 2 * lmax,
    approximation::AbstractDiagonalApproximation = IntegralDiagonal(),
)
    lmax >= 0 || throw(ArgumentError("atomic_operators requires lmax >= 0"))
    multipole_lmax >= 0 || throw(ArgumentError("atomic_operators requires multipole_lmax >= 0"))

    overlap = overlap_matrix(basis, grid)
    kinetic = kinetic_matrix(basis, grid)
    nuclear = nuclear_matrix(basis, grid; Z = Z)
    centrifugal_data = Matrix{Float64}[centrifugal_matrix(basis, grid; l = l) for l in 0:lmax]
    samples = _radial_multipole_samples(basis, grid)
    multipole_data = Matrix{Float64}[
        _radial_multipole_from_samples(samples, L, approximation) for L in 0:multipole_lmax
    ]
    shell_centers_r = Float64[Float64(value) for value in centers(basis)]
    source_manifest = _AtomicRadialSourceManifest(basis.spec, Float64(Z))
    return RadialAtomicOperators(
        overlap,
        kinetic,
        nuclear,
        centrifugal_data,
        multipole_data,
        shell_centers_r,
        source_manifest,
        approximation,
        samples,
    )
end

"""
    centrifugal(ops::RadialAtomicOperators, l::Int)

Return the precomputed centrifugal matrix for angular momentum `l`.
"""
function centrifugal(ops::RadialAtomicOperators, l::Int)
    l >= 0 || throw(ArgumentError("centrifugal requires l >= 0"))
    l < length(ops.centrifugal_data) || throw(BoundsError(ops.centrifugal_data, l + 1))
    return ops.centrifugal_data[l + 1]
end

"""
    multipole(ops::RadialAtomicOperators, L::Int)

Return the precomputed radial two-index IDA multipole matrix for multipole
order `L`.
"""
function multipole(ops::RadialAtomicOperators, L::Int)
    L >= 0 || throw(ArgumentError("multipole requires L >= 0"))
    L < length(ops.multipole_data) || throw(BoundsError(ops.multipole_data, L + 1))
    return ops.multipole_data[L + 1]
end

_stored_multipole_lmax(ops::RadialAtomicOperators) = length(ops.multipole_data) - 1
_radial_multipoles_extendable(ops::RadialAtomicOperators) = ops.multipole_samples !== nothing

# Stored multipole when available; otherwise evaluated from the retained quadrature samples
# with the same kernel (identical to what `atomic_operators(...; multipole_lmax = L)` stores).
function _radial_multipole_on_demand(ops::RadialAtomicOperators, L::Int)
    L >= 0 || throw(ArgumentError("multipole requires L >= 0"))
    L <= _stored_multipole_lmax(ops) && return ops.multipole_data[L + 1]
    _radial_multipoles_extendable(ops) || throw(
        ArgumentError(
            "radial multipole L=$(L) is not stored (stored L <= $(_stored_multipole_lmax(ops))) and these " *
            "RadialAtomicOperators carry no quadrature samples; rebuild them with " *
            "atomic_operators(...; multipole_lmax = $(L))",
        ),
    )
    return _radial_multipole_from_samples(ops.multipole_samples, L, ops.approximation)
end
