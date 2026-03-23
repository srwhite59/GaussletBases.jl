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

function _integral_diagonal_kernel_matrix(
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
    points, weights = _validate_radial_operator_grid(basis, grid)
    values = _basis_values_matrix(basis, points)

    if approximation isa IntegralDiagonal
        return _integral_diagonal_kernel_matrix(values, points, weights, L)
    end

    throw(ArgumentError("unsupported diagonal approximation $(typeof(approximation))"))
end

"""
    RadialAtomicOperators

Bundle of radial one-body matrices and precomputed radial multipole tables built
for a `RadialBasis` on an explicit `RadialQuadratureGrid`.

Use `ops.overlap`, `ops.kinetic`, and `ops.nuclear` for the fixed one-body
matrices, and `centrifugal(ops, l)` / `multipole(ops, L)` for the indexed
families.
"""
struct RadialAtomicOperators{A <: AbstractDiagonalApproximation}
    overlap::Matrix{Float64}
    kinetic::Matrix{Float64}
    nuclear::Matrix{Float64}
    centrifugal_data::Vector{Matrix{Float64}}
    multipole_data::Vector{Matrix{Float64}}
    shell_centers_r::Vector{Float64}
    approximation::A
end

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
                     approximation::AbstractDiagonalApproximation = IntegralDiagonal())

Build the high-level radial operator bundle for `basis` on the supplied
quadrature `grid`.

The bundle stores:
- `ops.overlap`
- `ops.kinetic`
- `ops.nuclear`
- `centrifugal(ops, l)` for `l = 0:lmax`
- `multipole(ops, L)` for `L = 0:(2 * lmax)`
"""
function atomic_operators(
    basis::RadialBasis,
    grid::RadialQuadratureGrid;
    Z::Real,
    lmax::Int = 0,
    approximation::AbstractDiagonalApproximation = IntegralDiagonal(),
)
    lmax >= 0 || throw(ArgumentError("atomic_operators requires lmax >= 0"))

    overlap = overlap_matrix(basis, grid)
    kinetic = kinetic_matrix(basis, grid)
    nuclear = nuclear_matrix(basis, grid; Z = Z)
    centrifugal_data = Matrix{Float64}[centrifugal_matrix(basis, grid; l = l) for l in 0:lmax]
    multipole_data = Matrix{Float64}[multipole_matrix(basis, grid; L = L, approximation = approximation) for L in 0:(2 * lmax)]
    shell_centers_r = Float64[Float64(value) for value in centers(basis)]
    return RadialAtomicOperators(
        overlap,
        kinetic,
        nuclear,
        centrifugal_data,
        multipole_data,
        shell_centers_r,
        approximation,
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
