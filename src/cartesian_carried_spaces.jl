module CartesianCarriedSpaces

import ..GaussletBases:
    AxisAlignedHomonuclearSquareLatticeQWBasis3D,
    BondAlignedDiatomicQWBasis3D,
    BondAlignedHomonuclearChainQWBasis3D,
    CartesianBasisRepresentation3D,
    _NestedFixedBlock3D,
    basis_representation
import ..GaussletBases.CartesianContractedParents:
    CartesianContractedParent3D,
    cartesian_contracted_parent,
    contracted_parent_basis,
    contracted_parent_coefficients,
    contracted_parent_dimension,
    contracted_parent_metadata,
    contracted_parent_parent_dimension,
    contracted_parent_units
import ..GaussletBases.CartesianParentGaussletBases:
    CartesianParentGaussletBasis3D,
    cartesian_parent_gausslet_basis,
    parent_axes,
    parent_axis_counts,
    parent_box,
    parent_dimension

export CartesianCarriedSpace3D,
       cartesian_carried_space,
       carried_space_parent,
       carried_space_contracted_parent,
       carried_space_representation,
       carried_space_diagnostics,
       carried_space_provenance

"""
    CartesianCarriedSpace3D

Internal normalized view of an existing Cartesian carried space.

The adapter packages the route-neutral parent identity, an optional contracted
parent, the existing Cartesian representation, and lightweight provenance. It
does not build overlap/Hamiltonian packets, does not own backend state, and does
not change geometry or construction policy.
"""
struct CartesianCarriedSpace3D{P,C,R,D,V}
    parent::P
    contracted_parent::C
    representation::R
    diagnostics::D
    provenance::V
end

carried_space_parent(space::CartesianCarriedSpace3D) = space.parent
carried_space_contracted_parent(space::CartesianCarriedSpace3D) =
    space.contracted_parent
carried_space_representation(space::CartesianCarriedSpace3D) =
    space.representation
carried_space_diagnostics(space::CartesianCarriedSpace3D) =
    space.diagnostics
carried_space_provenance(space::CartesianCarriedSpace3D) =
    space.provenance

_maybe_property(object, name::Symbol, default = nothing) =
    hasproperty(object, name) ? getproperty(object, name) : default

function _coefficient_storage_name(coefficients)
    isnothing(coefficients) && return nothing
    return nameof(typeof(coefficients))
end

function _carried_space_diagnostics(
    parent::CartesianParentGaussletBasis3D,
    contracted_parent::Union{Nothing,CartesianContractedParent3D},
    representation::CartesianBasisRepresentation3D,
)
    metadata = representation.metadata
    contracted_metadata =
        isnothing(contracted_parent) ? (;) : contracted_parent_metadata(contracted_parent)
    contracted_dimension_value =
        isnothing(contracted_parent) ? nothing : contracted_parent_dimension(contracted_parent)
    contracted_parent_dimension_value =
        isnothing(contracted_parent) ? nothing : contracted_parent_parent_dimension(contracted_parent)

    return (
        parent_dimension = parent_dimension(parent),
        parent_axis_counts = parent_axis_counts(parent),
        representation_basis_kind = metadata.basis_kind,
        representation_final_dimension = metadata.final_dimension,
        representation_contraction_kind = representation.contraction_kind,
        representation_coefficient_storage =
            _coefficient_storage_name(representation.coefficient_matrix),
        support_count =
            isnothing(representation.support_indices) ? nothing : length(representation.support_indices),
        has_contracted_parent = !isnothing(contracted_parent),
        contracted_dimension = contracted_dimension_value,
        contracted_parent_dimension = contracted_parent_dimension_value,
        contracted_parent_dimension_matches_parent =
            isnothing(contracted_parent) ? true :
            contracted_parent_dimension_value == parent_dimension(parent),
        contracted_dimension_matches_representation =
            isnothing(contracted_parent) ? true :
            contracted_dimension_value == metadata.final_dimension,
        contracted_coefficient_storage =
            isnothing(contracted_parent) ? nothing :
            _coefficient_storage_name(contracted_parent_coefficients(contracted_parent)),
        contraction_unit_count =
            isnothing(contracted_parent) ? 0 : length(contracted_parent_units(contracted_parent)),
        has_staged_sidecar = hasproperty(contracted_metadata, :staged_by_center_sidecar),
        staged_by_center_path =
            _maybe_property(contracted_metadata, :staged_by_center_path),
        dense_parent_matrix_used = false,
        heavy_metric_packet_built = false,
    )
end

function _carried_space_provenance(
    input_kind::Symbol,
    source,
    representation::CartesianBasisRepresentation3D,
)
    return (
        source = :cartesian_carried_space_adapter,
        input_kind = input_kind,
        source_type = nameof(typeof(source)),
        representation_basis_kind = representation.metadata.basis_kind,
        route_metadata = representation.metadata.route_metadata,
    )
end

function _same_parent_identity(
    left::CartesianParentGaussletBasis3D,
    right::CartesianParentGaussletBasis3D,
)
    left === right && return true
    return parent_box(left) == parent_box(right) &&
           parent_axis_counts(left) == parent_axis_counts(right) &&
           left.axis_sharing == right.axis_sharing &&
           isequal(parent_axes(left), parent_axes(right)) &&
           isequal(left.metadata, right.metadata)
end

function _validate_carried_parts(
    parent::CartesianParentGaussletBasis3D,
    contracted_parent::Union{Nothing,CartesianContractedParent3D},
    representation::CartesianBasisRepresentation3D,
)
    representation.metadata.parent_dimension == parent_dimension(parent) || throw(
        DimensionMismatch(
            "Cartesian carried-space representation parent dimension does not match parent identity",
        ),
    )
    isnothing(contracted_parent) && return nothing
    _same_parent_identity(contracted_parent_basis(contracted_parent), parent) || throw(
        ArgumentError("contracted parent and parent identity describe different Cartesian parents"),
    )
    contracted_parent_parent_dimension(contracted_parent) == parent_dimension(parent) || throw(
        DimensionMismatch("contracted parent rows do not match parent identity dimension"),
    )
    contracted_parent_dimension(contracted_parent) == representation.metadata.final_dimension || throw(
        DimensionMismatch(
            "contracted parent columns do not match Cartesian representation final dimension",
        ),
    )
    return nothing
end

function _cartesian_carried_space(
    source,
    input_kind::Symbol,
    parent::CartesianParentGaussletBasis3D,
    contracted_parent::Union{Nothing,CartesianContractedParent3D},
    representation::CartesianBasisRepresentation3D,
)
    _validate_carried_parts(parent, contracted_parent, representation)
    diagnostics = _carried_space_diagnostics(parent, contracted_parent, representation)
    provenance = _carried_space_provenance(input_kind, source, representation)
    return CartesianCarriedSpace3D(
        parent,
        contracted_parent,
        representation,
        diagnostics,
        provenance,
    )
end

"""
    cartesian_carried_space(fixed_block::_NestedFixedBlock3D)

Return the normalized Cartesian carried-space view of a nested fixed block.

The fixed block's staged by-center sidecar, when present, is preserved through
the `CartesianContractedParent3D` metadata and units so existing staged metric
contraction remains available.
"""
function cartesian_carried_space(fixed_block::_NestedFixedBlock3D)
    parent = cartesian_parent_gausslet_basis(fixed_block)
    contracted_parent = cartesian_contracted_parent(fixed_block)
    representation = basis_representation(fixed_block)
    return _cartesian_carried_space(
        fixed_block,
        :nested_fixed_block,
        parent,
        contracted_parent,
        representation,
    )
end

function _direct_product_carried_space(source, input_kind::Symbol)
    parent = cartesian_parent_gausslet_basis(source)
    representation = basis_representation(source)
    return _cartesian_carried_space(
        source,
        input_kind,
        parent,
        nothing,
        representation,
    )
end

cartesian_carried_space(basis::BondAlignedDiatomicQWBasis3D) =
    _direct_product_carried_space(basis, :bond_aligned_diatomic_qw_basis)

cartesian_carried_space(basis::BondAlignedHomonuclearChainQWBasis3D) =
    _direct_product_carried_space(basis, :bond_aligned_homonuclear_chain_qw_basis)

cartesian_carried_space(basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D) =
    _direct_product_carried_space(basis, :axis_aligned_homonuclear_square_lattice_qw_basis)

"""
    cartesian_carried_space(representation; parent, contracted_parent=nothing)

Normalize an already-built Cartesian representation when the caller can supply
the unambiguous parent identity and optional contracted parent. This path is for
internal plumbing only; ambiguous representations should stay explicit rather
than trying to infer parent axes from labels or centers.
"""
function cartesian_carried_space(
    representation::CartesianBasisRepresentation3D;
    parent::CartesianParentGaussletBasis3D,
    contracted_parent::Union{Nothing,CartesianContractedParent3D} = nothing,
    input_kind::Symbol = :cartesian_basis_representation,
)
    return _cartesian_carried_space(
        representation,
        input_kind,
        parent,
        contracted_parent,
        representation,
    )
end

end
