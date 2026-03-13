module Gausslets

using LinearAlgebra

export AbstractFunction1D,
       AbstractPrimitiveFunction1D,
       AbstractBasisFunction1D,
       AbstractCoordinateMapping,
       AbstractBasisSpec,
       UniformBasisSpec,
       HalfLineBasisSpec,
       RadialBasisSpec,
       Gaussian,
       HalfLineGaussian,
       XGaussian,
       Distorted,
       GaussletFamily,
       Gausslet,
       UniformBasis,
       HalfLineBasis,
       RadialBasis,
       BoundaryGausslet,
       RadialGausslet,
       StencilTerm,
       FunctionStencil,
       value,
       direct_value,
       derivative,
       center,
       reference_center,
       integral_weight,
       stencil,
       build_basis,
       basis_spec,
       family,
       mapping,
       centers,
       reference_centers,
       integral_weights,
       coefficients,
       primitives,
       terms,
       IdentityMapping,
       AsinhMapping,
       uofx,
       xofu,
       dudx,
       du2dx2

"""
    AbstractFunction1D

Abstract supertype for callable one-dimensional function objects.
"""
abstract type AbstractFunction1D end

"""
    AbstractPrimitiveFunction1D <: AbstractFunction1D

Abstract supertype for lowest-level primitive function objects.
"""
abstract type AbstractPrimitiveFunction1D <: AbstractFunction1D end

"""
    AbstractBasisFunction1D <: AbstractFunction1D

Abstract supertype for higher-level callable basis functions.
"""
abstract type AbstractBasisFunction1D <: AbstractFunction1D end

"""
    AbstractCoordinateMapping

Abstract supertype for coordinate maps between physical `x` and reference `u`.
"""
abstract type AbstractCoordinateMapping end

"""
    AbstractBasisSpec

Placeholder abstract supertype for later basis-construction recipes.
"""
abstract type AbstractBasisSpec end

function value end
function direct_value end
function derivative end
function center end
function reference_center end
function integral_weight end
function stencil end
function build_basis end

function basis_spec end
function family end
function mapping end
function centers end
function reference_centers end
function integral_weights end

function coefficients end
function primitives end
function terms end

function uofx end
function xofu end
function dudx end
function du2dx2 end

(f::AbstractFunction1D)(x::Real) = value(f, x)
(mapping::AbstractCoordinateMapping)(x::Real) = uofx(mapping, x)

value(f::AbstractFunction1D, x::Real) = direct_value(f, x)
direct_value(f::AbstractFunction1D, x::Real) = stencil(f)(x)

function derivative(f::AbstractFunction1D, x::Real; order::Int = 1)
    order >= 0 || throw(ArgumentError("derivative order must be nonnegative"))
    order == 0 && return value(f, x)
    st = stencil(f)
    total = 0.0
    for i in eachindex(coefficients(st))
        total += coefficients(st)[i] * derivative(primitives(st)[i], x; order = order)
    end
    return total
end

reference_center(f::AbstractFunction1D) = center(f)

function integral_weight(f::AbstractFunction1D)
    st = stencil(f)
    total = 0.0
    for term in terms(st)
        total += term.coefficient * integral_weight(term.primitive)
    end
    return total
end

include("mappings.jl")
include("stencils.jl")
include("functions.jl")
include("families.jl")
include("bases.jl")

end
