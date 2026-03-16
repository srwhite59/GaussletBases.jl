const _PRIMITIVE_MATRIX_TOL = 1.0e-10
const _PRIMITIVE_MATRIX_MAXITER = 5

abstract type _PrimitiveMatrixBackend end

struct _AnalyticPrimitiveMatrixBackend <: _PrimitiveMatrixBackend
end

struct _NumericalPrimitiveMatrixBackend <: _PrimitiveMatrixBackend
end

"""
    PrimitiveSet1D(primitives; name=nothing, labels=nothing)

Explicit ordered set of lowest-level primitives used as a common matrix-building
layer for ordinary, mapped, and later nested constructions.

The set may contain plain primitives such as `Gaussian` or explicit mapped
primitives such as `Distorted(Gaussian(...), mapping)`. Matrix builders such as
`overlap_matrix` and `kinetic_matrix` choose analytic or numerical evaluation
internally, depending on what the primitive content allows.
"""
struct PrimitiveSet1D
    primitive_data::Vector{AbstractPrimitiveFunction1D}
    name_value::Union{Nothing, Symbol}
    label_data::Vector{String}
end

function PrimitiveSet1D(
    primitive_data::AbstractVector{<:AbstractPrimitiveFunction1D};
    name::Union{Nothing, Symbol} = nothing,
    labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    primitive_list = AbstractPrimitiveFunction1D[primitive for primitive in primitive_data]
    label_list =
        if labels === nothing
            ["mu$(i)" for i in eachindex(primitive_list)]
        else
            length(labels) == length(primitive_list) ||
                throw(ArgumentError("primitive labels must match the primitive count"))
            String[String(label) for label in labels]
        end
    return PrimitiveSet1D(primitive_list, name, label_list)
end

Base.length(set::PrimitiveSet1D) = length(set.primitive_data)
Base.getindex(set::PrimitiveSet1D, index::Integer) = set.primitive_data[index]
primitives(set::PrimitiveSet1D) = set.primitive_data
primitive_set(set::PrimitiveSet1D) = set

function Base.show(io::IO, set::PrimitiveSet1D)
    print(io, "PrimitiveSet1D(length=", length(set))
    set.name_value === nothing || print(io, ", name=:", set.name_value)
    print(io, ")")
end

"""
    BasisMetadata1D

Minimal basis-level metadata bundle for downstream export or reconstruction of
matrix workflows from the shared primitive layer.
"""
struct BasisMetadata1D{M}
    basis_kind::Symbol
    family_name::Union{Nothing, Symbol}
    mapping_value::M
    center_data::Vector{Float64}
    reference_center_data::Vector{Float64}
    integral_weight_data::Vector{Float64}
    basis_labels::Vector{String}
    primitive_set::PrimitiveSet1D
    coefficient_matrix::Matrix{Float64}
end

function Base.show(io::IO, metadata::BasisMetadata1D)
    print(
        io,
        "BasisMetadata1D(kind=:",
        metadata.basis_kind,
        ", nbasis=",
        length(metadata.basis_labels),
        ", nprimitive=",
        length(metadata.primitive_set),
        ")",
    )
end

"""
    BasisRepresentation1D

Small in-memory downstream-consumer object bundling basis metadata, the shared
primitive layer, the contraction matrix, and selected primitive- and
basis-level matrices built through that layer.
"""
struct BasisRepresentation1D{M, PT <: NamedTuple, BT <: NamedTuple}
    metadata::BasisMetadata1D{M}
    primitive_set::PrimitiveSet1D
    coefficient_matrix::Matrix{Float64}
    primitive_matrices::PT
    basis_matrices::BT
end

function Base.show(io::IO, representation::BasisRepresentation1D)
    print(
        io,
        "BasisRepresentation1D(kind=:",
        representation.metadata.basis_kind,
        ", nbasis=",
        size(representation.coefficient_matrix, 2),
        ", nprimitive=",
        length(representation.primitive_set),
        ", operators=",
        collect(keys(representation.basis_matrices)),
        ")",
    )
end

function _basis_kind(::UniformBasis)
    return :uniform
end

function _basis_kind(::MappedUniformBasis)
    return :mapped_uniform
end

function _basis_kind(::HalfLineBasis)
    return :halfline
end

function _basis_kind(::RadialBasis)
    return :radial
end

function _basis_family_name(basis)
    return family(basis).name
end

function _basis_label_vector(length_value::Int)
    return ["phi$(i)" for i in 1:length_value]
end

"""
    primitive_set(basis)

Return the explicit shared primitive layer underlying `basis` as a
`PrimitiveSet1D`.

Use `stencil_matrix(basis)` as the contraction map from this primitive layer to
the final basis functions.
"""
function primitive_set(basis::Union{UniformBasis, MappedUniformBasis, HalfLineBasis, RadialBasis})
    return PrimitiveSet1D(
        primitives(basis);
        name = Symbol(_basis_kind(basis), "_primitives"),
    )
end

function basis_metadata(basis::Union{UniformBasis, MappedUniformBasis, HalfLineBasis, RadialBasis})
    primitive_layer = primitive_set(basis)
    return BasisMetadata1D(
        _basis_kind(basis),
        _basis_family_name(basis),
        mapping(basis),
        Float64[Float64(value) for value in centers(basis)],
        Float64[Float64(value) for value in reference_centers(basis)],
        Float64[Float64(value) for value in integral_weights(basis)],
        _basis_label_vector(length(basis)),
        primitive_layer,
        Matrix{Float64}(stencil_matrix(basis)),
    )
end

basis_metadata(representation::BasisRepresentation1D) = representation.metadata
primitive_set(representation::BasisRepresentation1D) = representation.primitive_set
stencil_matrix(representation::BasisRepresentation1D) = representation.coefficient_matrix

function _representation_operator_matrix(set::PrimitiveSet1D, operator_name::Symbol)
    operator_name === :overlap && return overlap_matrix(set)
    operator_name === :position && return position_matrix(set)
    operator_name === :kinetic && return kinetic_matrix(set)
    throw(ArgumentError("unsupported basis representation operator :$operator_name"))
end

"""
    basis_representation(basis; operators = (:overlap, :position, :kinetic))

Build a small in-memory representation object for `basis` that a downstream
consumer can use without depending on GaussletBases internals.

The returned `BasisRepresentation1D` stores:

- `metadata`
- `primitive_set`
- `coefficient_matrix`
- `primitive_matrices`
- `basis_matrices`

The operator matrices are built through the shared primitive layer and
contracted upward with the existing basis stencil matrix.
"""
function basis_representation(
    basis::Union{UniformBasis, MappedUniformBasis, HalfLineBasis, RadialBasis};
    operators = (:overlap, :position, :kinetic),
)
    operator_names = Tuple(Symbol(operator_name) for operator_name in operators)
    length(unique(operator_names)) == length(operator_names) ||
        throw(ArgumentError("basis representation operator names must be unique"))

    metadata = basis_metadata(basis)
    primitive_layer = metadata.primitive_set
    coefficient_matrix = metadata.coefficient_matrix

    primitive_matrix_values =
        Tuple(_representation_operator_matrix(primitive_layer, operator_name) for operator_name in operator_names)
    basis_matrix_values =
        Tuple(contract_primitive_matrix(basis, matrix) for matrix in primitive_matrix_values)

    return BasisRepresentation1D(
        metadata,
        primitive_layer,
        coefficient_matrix,
        NamedTuple{operator_names}(primitive_matrix_values),
        NamedTuple{operator_names}(basis_matrix_values),
    )
end

function _primitive_support_bounds(primitive::Gaussian)
    return _reference_bounds(primitive)
end

function _primitive_support_bounds(primitive::HalfLineGaussian)
    return _reference_bounds(primitive)
end

function _primitive_support_bounds(primitive::XGaussian)
    return _reference_bounds(primitive)
end

function _primitive_support_bounds(primitive::Distorted)
    ulo, uhi = _reference_bounds(primitive.primitive)
    return xofu(primitive.mapping, ulo), xofu(primitive.mapping, uhi)
end

function _primitive_length_scale(primitive::Gaussian)
    return primitive.width
end

function _primitive_length_scale(primitive::HalfLineGaussian)
    return primitive.width
end

function _primitive_length_scale(primitive::XGaussian)
    return primitive.alpha
end

function _primitive_length_scale(primitive::Distorted)
    local_scale = _primitive_length_scale(primitive.primitive)
    return local_scale / dudx(primitive.mapping, center(primitive))
end

function _primitive_set_bounds(set::PrimitiveSet1D)
    lowers = Float64[]
    uppers = Float64[]
    for primitive in primitives(set)
        xlo, xhi = _primitive_support_bounds(primitive)
        push!(lowers, xlo)
        push!(uppers, xhi)
    end
    return minimum(lowers), maximum(uppers)
end

function _primitive_matrix_start_h(set::PrimitiveSet1D)
    scales = Float64[_primitive_length_scale(primitive) for primitive in primitives(set)]
    xlo, xhi = _primitive_set_bounds(set)
    span = xhi - xlo
    span > 0.0 || throw(ArgumentError("primitive set support must have positive width"))
    return min(minimum(scales) / 10.0, span / 600.0)
end

function _primitive_sample_matrix(
    set::PrimitiveSet1D,
    points::AbstractVector{Float64};
    derivative_order::Int = 0,
)
    samples = zeros(Float64, length(points), length(set))
    for mu in eachindex(primitives(set))
        if derivative_order == 0
            samples[:, mu] = [value(primitives(set)[mu], point) for point in points]
        else
            samples[:, mu] = [
                derivative(primitives(set)[mu], point; order = derivative_order)
                for point in points
            ]
        end
    end
    return samples
end

function _symmetrize_primitive_matrix(matrix::AbstractMatrix{<:Real})
    return 0.5 .* (Matrix{Float64}(matrix) .+ Matrix{Float64}(transpose(matrix)))
end

function _primitive_overlap_matrix(
    set::PrimitiveSet1D,
    ::_NumericalPrimitiveMatrixBackend;
    h = nothing,
)
    xlo, xhi = _primitive_set_bounds(set)
    h_try = h === nothing ? _primitive_matrix_start_h(set) : Float64(h)
    h_try > 0.0 || throw(ArgumentError("numerical primitive overlap requires h > 0"))

    previous = nothing
    current = nothing
    for _ in 1:_PRIMITIVE_MATRIX_MAXITER
        points, weights = _make_midpoint_grid(xlo, xhi, h_try)
        values = _primitive_sample_matrix(set, points)
        current = _symmetrize_primitive_matrix(transpose(values) * (weights .* values))
        if previous !== nothing && norm(current - previous, Inf) <= _PRIMITIVE_MATRIX_TOL
            return current
        end
        previous = current
        h_try /= 2.0
    end
    return current
end

function _primitive_kinetic_matrix(
    set::PrimitiveSet1D,
    ::_NumericalPrimitiveMatrixBackend;
    h = nothing,
)
    xlo, xhi = _primitive_set_bounds(set)
    h_try = h === nothing ? _primitive_matrix_start_h(set) : Float64(h)
    h_try > 0.0 || throw(ArgumentError("numerical primitive kinetic matrix requires h > 0"))

    previous = nothing
    current = nothing
    for _ in 1:_PRIMITIVE_MATRIX_MAXITER
        points, weights = _make_midpoint_grid(xlo, xhi, h_try)
        derivatives = _primitive_sample_matrix(set, points; derivative_order = 1)
        current = _symmetrize_primitive_matrix(0.5 .* (transpose(derivatives) * (weights .* derivatives)))
        if previous !== nothing && norm(current - previous, Inf) <= _PRIMITIVE_MATRIX_TOL
            return current
        end
        previous = current
        h_try /= 2.0
    end
    return current
end

function _primitive_position_matrix(
    set::PrimitiveSet1D,
    ::_NumericalPrimitiveMatrixBackend;
    h = nothing,
)
    xlo, xhi = _primitive_set_bounds(set)
    h_try = h === nothing ? _primitive_matrix_start_h(set) : Float64(h)
    h_try > 0.0 || throw(ArgumentError("numerical primitive position matrix requires h > 0"))

    previous = nothing
    current = nothing
    for _ in 1:_PRIMITIVE_MATRIX_MAXITER
        points, weights = _make_midpoint_grid(xlo, xhi, h_try)
        values = _primitive_sample_matrix(set, points)
        current = _symmetrize_primitive_matrix(
            transpose(values) * (((weights .* points)) .* values),
        )
        if previous !== nothing && norm(current - previous, Inf) <= _PRIMITIVE_MATRIX_TOL
            return current
        end
        previous = current
        h_try /= 2.0
    end
    return current
end

function _gaussian_overlap(a::Gaussian, b::Gaussian)
    sigma2 = a.width^2 + b.width^2
    prefactor = sqrt(2.0 * pi) * a.width * b.width / sqrt(sigma2)
    return prefactor * exp(-0.5 * (a.center_value - b.center_value)^2 / sigma2)
end

function _gaussian_kinetic(a::Gaussian, b::Gaussian)
    sigma2 = a.width^2 + b.width^2
    overlap_value = _gaussian_overlap(a, b)
    delta = a.center_value - b.center_value
    return 0.5 * overlap_value * (1.0 / sigma2) * (1.0 - delta^2 / sigma2)
end

function _gaussian_position(a::Gaussian, b::Gaussian)
    sigma2 = a.width^2 + b.width^2
    weighted_center =
        (a.center_value * b.width^2 + b.center_value * a.width^2) / sigma2
    return weighted_center * _gaussian_overlap(a, b)
end

function _primitive_overlap_matrix(set::PrimitiveSet1D, ::_AnalyticPrimitiveMatrixBackend)
    matrix = zeros(Float64, length(set), length(set))
    for a in 1:length(set)
        pa = primitives(set)[a]
        for b in a:length(set)
            pb = primitives(set)[b]
            value_ab = _gaussian_overlap(pa, pb)
            matrix[a, b] = value_ab
            matrix[b, a] = value_ab
        end
    end
    return matrix
end

function _primitive_position_matrix(set::PrimitiveSet1D, ::_AnalyticPrimitiveMatrixBackend)
    matrix = zeros(Float64, length(set), length(set))
    for a in 1:length(set)
        pa = primitives(set)[a]
        for b in a:length(set)
            pb = primitives(set)[b]
            value_ab = _gaussian_position(pa, pb)
            matrix[a, b] = value_ab
            matrix[b, a] = value_ab
        end
    end
    return matrix
end

function _primitive_kinetic_matrix(set::PrimitiveSet1D, ::_AnalyticPrimitiveMatrixBackend)
    matrix = zeros(Float64, length(set), length(set))
    for a in 1:length(set)
        pa = primitives(set)[a]
        for b in a:length(set)
            pb = primitives(set)[b]
            value_ab = _gaussian_kinetic(pa, pb)
            matrix[a, b] = value_ab
            matrix[b, a] = value_ab
        end
    end
    return matrix
end

function _supports_analytic_gaussian_backend(set::PrimitiveSet1D)
    return all(primitive -> primitive isa Gaussian, primitives(set))
end

function _select_primitive_matrix_backend(set::PrimitiveSet1D)
    return _supports_analytic_gaussian_backend(set) ?
           _AnalyticPrimitiveMatrixBackend() :
           _NumericalPrimitiveMatrixBackend()
end

"""
    overlap_matrix(set::PrimitiveSet1D)

Build the primitive overlap matrix for `set`.

The public call chooses the available backend automatically. Plain full-line
Gaussian sets use analytic formulas. Distorted or otherwise unsupported
primitive content falls back to numerical quadrature over the explicit
primitive support.
"""
function overlap_matrix(set::PrimitiveSet1D)
    return _primitive_overlap_matrix(set, _select_primitive_matrix_backend(set))
end

"""
    position_matrix(set::PrimitiveSet1D)

Build the primitive one-body position matrix

    <phi_mu | x | phi_nu>

for `set`.

As with `overlap_matrix(set)`, the public call chooses the available backend
automatically.
"""
function position_matrix(set::PrimitiveSet1D)
    return _primitive_position_matrix(set, _select_primitive_matrix_backend(set))
end

"""
    kinetic_matrix(set::PrimitiveSet1D)

Build the primitive one-body kinetic matrix

    <phi_mu | -0.5 d^2/dx^2 | phi_nu>

for `set`.

As with `overlap_matrix(set)`, the public call chooses the available backend
automatically.
"""
function kinetic_matrix(set::PrimitiveSet1D)
    return _primitive_kinetic_matrix(set, _select_primitive_matrix_backend(set))
end
