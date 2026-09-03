"""
    StencilTerm

Coefficient and primitive pair from a `FunctionStencil`.
"""
struct StencilTerm
    coefficient::Float64
    primitive::AbstractPrimitiveFunction1D

    function StencilTerm(coefficient::Real, primitive::AbstractPrimitiveFunction1D)
        new(Float64(coefficient), primitive)
    end
end

"""
    FunctionStencil

Exact public stencil of a one-dimensional function in terms of the lowest-level
primitive functions exposed by the library.
"""
struct FunctionStencil{C <: AbstractVector{Float64}, P <: AbstractVector{AbstractPrimitiveFunction1D}}
    coefficient_data::C
    primitive_data::P
end

function FunctionStencil(
    coefficient_data::AbstractVector{<:Real},
    primitive_data::AbstractVector{<:AbstractPrimitiveFunction1D},
)
    length(coefficient_data) == length(primitive_data) ||
        throw(ArgumentError("coefficient and primitive counts must match"))
    coeffs =
        coefficient_data isa AbstractVector{Float64} ?
        coefficient_data :
        Float64[Float64(c) for c in coefficient_data]
    prims =
        primitive_data isa AbstractVector{AbstractPrimitiveFunction1D} ?
        primitive_data :
        AbstractPrimitiveFunction1D[p for p in primitive_data]
    return FunctionStencil{typeof(coeffs), typeof(prims)}(coeffs, prims)
end

Base.length(st::FunctionStencil) = length(st.coefficient_data)

function Base.getindex(st::FunctionStencil, index::Integer)
    return StencilTerm(st.coefficient_data[index], st.primitive_data[index])
end

coefficients(st::FunctionStencil) = st.coefficient_data
primitives(st::FunctionStencil) = st.primitive_data

function terms(st::FunctionStencil)
    return [StencilTerm(st.coefficient_data[i], st.primitive_data[i]) for i in eachindex(st.coefficient_data)]
end

function (st::FunctionStencil)(x::Real)
    total = 0.0
    for i in eachindex(st.coefficient_data)
        total += st.coefficient_data[i] * value(st.primitive_data[i], x)
    end
    return total
end
