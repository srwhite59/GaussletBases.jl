module Gausslets

export AbstractFunction1D,
       AbstractPrimitiveFunction1D,
       AbstractBasisFunction1D,
       AbstractCoordinateMapping,
       AbstractBasisSpec,
       Gaussian,
       HalfLineGaussian,
       XGaussian,
       Distorted,
       GaussletFamily,
       Gausslet,
       StencilTerm,
       FunctionStencil,
       value,
       direct_value,
       derivative,
       center,
       integral_weight,
       stencil,
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
function integral_weight end
function stencil end

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
    throw(ArgumentError("derivative order $(order) is not implemented for $(typeof(f))"))
end

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

end
