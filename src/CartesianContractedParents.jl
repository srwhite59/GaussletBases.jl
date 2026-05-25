module CartesianContractedParents

import ..GaussletBases: _CartesianCoefficientMap,
                         _NestedFixedBlock3D,
                         _cartesian_coefficient_map_storage
import ..GaussletBases.CartesianParentGaussletBases:
    CartesianParentGaussletBasis3D,
    cartesian_parent_gausslet_basis,
    parent_dimension

export CartesianContractionUnit3D,
       CartesianContractedParent3D,
       CartesianContractedParentStructuralAudit,
       cartesian_contracted_parent,
       contracted_parent_basis,
       contracted_parent_coefficients,
       contracted_parent_units,
       contracted_parent_metadata,
       contracted_parent_parent_dimension,
       contracted_parent_dimension,
       contraction_unit_role,
       contraction_unit_support_indices,
       contraction_unit_column_range,
       contraction_unit_metadata,
       contracted_parent_unit_column_ranges,
       contracted_parent_unit_support_indices,
       contracted_parent_support_indices,
       contracted_parent_structural_audit

"""
    CartesianContractionUnit3D

Internal provenance record for a group of contracted-parent columns.

The global coefficient matrix on `CartesianContractedParent3D` remains the
source of truth. Unit support and column ranges are structural metadata only:
support may overlap, may be incomplete, and does not imply orthonormality.
"""
struct CartesianContractionUnit3D{M}
    role::Symbol
    support_indices::Vector{Int}
    column_range::UnitRange{Int}
    metadata::M
end

function CartesianContractionUnit3D(
    role::Symbol,
    support_indices::AbstractVector{<:Integer},
    column_range::UnitRange{<:Integer};
    metadata = (;),
)
    isempty(column_range) && throw(
        ArgumentError("Cartesian contraction units require a non-empty column range"),
    )
    return CartesianContractionUnit3D(
        role,
        Int[Int(index) for index in support_indices],
        Int(first(column_range)):Int(last(column_range)),
        metadata,
    )
end

"""
    CartesianContractedParent3D

Internal identity object for columns formed as linear combinations of a full
Cartesian parent gausslet lattice.

This object deliberately stores no backend state, overlap/H/V matrices,
Gaussian supplements, residual columns, or operator packets. The coefficient
matrix is the source of truth; units are provenance only.
"""
struct CartesianContractedParent3D{P<:CartesianParentGaussletBasis3D,C<:_CartesianCoefficientMap,U<:AbstractVector,M}
    parent::P
    coefficient_matrix::C
    units::U
    metadata::M
end

function _validate_unit_column_ranges(
    units::AbstractVector{<:CartesianContractionUnit3D},
    contracted_dimension::Int,
)
    for unit in units
        first(unit.column_range) >= 1 || throw(
            ArgumentError("contraction unit column ranges must start inside contracted columns"),
        )
        last(unit.column_range) <= contracted_dimension || throw(
            ArgumentError("contraction unit column ranges must lie inside contracted columns"),
        )
    end
    return units
end

function CartesianContractedParent3D(
    parent::CartesianParentGaussletBasis3D,
    coefficients::AbstractMatrix{<:Real};
    units::AbstractVector{<:CartesianContractionUnit3D} = CartesianContractionUnit3D[],
    metadata = (;),
)
    coefficient_matrix = _cartesian_coefficient_map_storage(coefficients)
    size(coefficient_matrix, 1) == parent_dimension(parent) || throw(
        DimensionMismatch("contracted parent coefficient rows must match parent dimension"),
    )
    unit_values = CartesianContractionUnit3D[units...]
    _validate_unit_column_ranges(unit_values, size(coefficient_matrix, 2))
    return CartesianContractedParent3D(
        parent,
        coefficient_matrix,
        unit_values,
        metadata,
    )
end

function CartesianContractedParent3D(
    fixed_block::_NestedFixedBlock3D;
    metadata = (source = :nested_fixed_block,),
)
    parent = cartesian_parent_gausslet_basis(fixed_block)
    coefficients = fixed_block.coefficient_matrix
    unit = CartesianContractionUnit3D(
        :nested_fixed_block,
        fixed_block.support_indices,
        1:size(coefficients, 2);
        metadata = (source = :nested_fixed_block,),
    )
    return CartesianContractedParent3D(
        parent,
        coefficients;
        units = [unit],
        metadata,
    )
end

cartesian_contracted_parent(
    parent::CartesianParentGaussletBasis3D,
    coefficients::AbstractMatrix{<:Real};
    kwargs...,
) = CartesianContractedParent3D(parent, coefficients; kwargs...)

cartesian_contracted_parent(fixed_block::_NestedFixedBlock3D; kwargs...) =
    CartesianContractedParent3D(fixed_block; kwargs...)

contracted_parent_basis(parent::CartesianContractedParent3D) = parent.parent
contracted_parent_coefficients(parent::CartesianContractedParent3D) = parent.coefficient_matrix
contracted_parent_units(parent::CartesianContractedParent3D) = parent.units
contracted_parent_metadata(parent::CartesianContractedParent3D) = parent.metadata
contracted_parent_parent_dimension(parent::CartesianContractedParent3D) =
    parent_dimension(parent.parent)
contracted_parent_dimension(parent::CartesianContractedParent3D) =
    size(parent.coefficient_matrix, 2)

contraction_unit_role(unit::CartesianContractionUnit3D) = unit.role
contraction_unit_support_indices(unit::CartesianContractionUnit3D) = unit.support_indices
contraction_unit_column_range(unit::CartesianContractionUnit3D) = unit.column_range
contraction_unit_metadata(unit::CartesianContractionUnit3D) = unit.metadata

contracted_parent_unit_column_ranges(parent::CartesianContractedParent3D) =
    UnitRange{Int}[unit.column_range for unit in parent.units]

contracted_parent_unit_support_indices(parent::CartesianContractedParent3D) =
    Vector{Int}[copy(unit.support_indices) for unit in parent.units]

function contracted_parent_support_indices(parent::CartesianContractedParent3D)
    indices = Int[]
    for unit in parent.units
        append!(indices, unit.support_indices)
    end
    return indices
end

struct CartesianContractedParentStructuralAudit
    parent_dimension::Int
    contracted_dimension::Int
    unit_count::Int
    support_entry_count::Int
    unique_support_count::Int
    duplicate_support_count::Int
    missing_support_count::Int
    outside_support_count::Int
    column_entry_count::Int
    unique_column_count::Int
    duplicate_column_count::Int
    missing_column_count::Int
    outside_column_count::Int
    support_complete::Bool
    column_ranges_cover_contract::Bool
    structural_ok::Bool
end

function _entry_counts(values::AbstractVector{Int}, valid_range::UnitRange{Int})
    inside = Int[value for value in values if value in valid_range]
    outside_count = length(values) - length(inside)
    unique_inside = unique(inside)
    duplicate_count = length(inside) - length(unique_inside)
    missing_count = length(valid_range) - length(unique_inside)
    return (
        entry_count = length(values),
        unique_count = length(unique_inside),
        duplicate_count = duplicate_count,
        missing_count = missing_count,
        outside_count = outside_count,
    )
end

function contracted_parent_structural_audit(parent::CartesianContractedParent3D)
    parent_dim = contracted_parent_parent_dimension(parent)
    contracted_dim = contracted_parent_dimension(parent)
    support_counts = _entry_counts(
        contracted_parent_support_indices(parent),
        1:parent_dim,
    )
    columns = Int[]
    for unit in parent.units
        append!(columns, collect(unit.column_range))
    end
    column_counts = _entry_counts(columns, 1:contracted_dim)
    support_complete = support_counts.missing_count == 0 && support_counts.outside_count == 0
    column_ranges_cover_contract =
        column_counts.missing_count == 0 &&
        column_counts.outside_count == 0 &&
        column_counts.duplicate_count == 0
    structural_ok = support_counts.outside_count == 0 && column_ranges_cover_contract
    return CartesianContractedParentStructuralAudit(
        parent_dim,
        contracted_dim,
        length(parent.units),
        support_counts.entry_count,
        support_counts.unique_count,
        support_counts.duplicate_count,
        support_counts.missing_count,
        support_counts.outside_count,
        column_counts.entry_count,
        column_counts.unique_count,
        column_counts.duplicate_count,
        column_counts.missing_count,
        column_counts.outside_count,
        support_complete,
        column_ranges_cover_contract,
        structural_ok,
    )
end

end
