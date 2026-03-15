module GaussletBases

using LinearAlgebra

export AbstractFunction1D,
       AbstractPrimitiveFunction1D,
       AbstractBasisFunction1D,
       AbstractCoordinateMapping,
       AbstractBasisSpec,
       AbstractDiagonalApproximation,
       PrimitiveSet1D,
       BasisMetadata1D,
       BasisRepresentation1D,
       BasisBox1D,
       BasisPartition1D,
       HierarchicalBasisBox1D,
       HierarchicalBasisPartition1D,
       LeafLocalPGDG1D,
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
       moment_center,
       integral_weight,
       stencil,
       stencil_matrix,
       build_basis,
       basis_metadata,
       basis_representation,
       basis_partition,
       hierarchical_partition,
       build_leaf_pgdg,
       refine_partition,
        primitive_set,
        boxes,
        leaf_boxes,
        leaf_primitive_indices,
        box_indices,
        box_level,
        box_parent,
       box_children,
       box_block,
       box_coupling,
       basis_spec,
       family,
       mapping,
       centers,
       reference_centers,
       integral_weights,
       contract_primitive_vector,
       contract_primitive_diagonal,
       contract_primitive_matrix,
       RadialQuadratureGrid,
       radial_quadrature,
       quadrature_points,
       quadrature_weights,
       basis_diagnostics,
       IntegralDiagonal,
       overlap_matrix,
       position_matrix,
       kinetic_matrix,
       nuclear_matrix,
       centrifugal_matrix,
       multipole_matrix,
       RadialAtomicOperators,
       atomic_operators,
       centrifugal,
       multipole,
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

Abstract supertype for public basis-construction recipes.
"""
abstract type AbstractBasisSpec end

"""
    AbstractDiagonalApproximation

Abstract supertype for supported diagonal-approximation choices in the radial
electron-electron operator layer.
"""
abstract type AbstractDiagonalApproximation end

function value end
function direct_value end
function derivative end
function center end
function reference_center end
function moment_center end
function integral_weight end
function stencil end
function stencil_matrix end
function build_basis end
function basis_metadata end
function basis_representation end
function basis_partition end
function hierarchical_partition end
function build_leaf_pgdg end
function refine_partition end
function primitive_set end
function boxes end
function leaf_boxes end
function leaf_primitive_indices end
function box_indices end
function box_level end
function box_parent end
function box_children end
function box_block end
function box_coupling end

function basis_spec end
function family end
function mapping end
function centers end
function reference_centers end
function integral_weights end
function contract_primitive_vector end
function contract_primitive_diagonal end
function contract_primitive_matrix end
function radial_quadrature end
function quadrature_points end
function quadrature_weights end
function basis_diagnostics end
function overlap_matrix end
function position_matrix end
function kinetic_matrix end
function nuclear_matrix end
function centrifugal_matrix end
function multipole_matrix end
function atomic_operators end
function centrifugal end
function multipole end

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
include("internal/wavelet_filters.jl")
include("families.jl")
include("bases.jl")
include("quadrature.jl")
include("primitive_sets.jl")
include("partitions.jl")
include("hierarchical_partitions.jl")
include("leaf_pgdg.jl")
include("diagnostics.jl")
include("operators.jl")

end
