module CartesianParentGaussletBases

import ..GaussletBases: MappedUniformBasis,
                         BondAlignedDiatomicQWBasis3D,
                         BondAlignedHomonuclearChainQWBasis3D,
                         AxisAlignedHomonuclearSquareLatticeQWBasis3D,
                         _NestedFixedBlock3D,
                         centers

export CartesianParentGaussletBasis3D,
       parent_axes,
       axis_basis,
       parent_box,
       parent_axis_counts,
       parent_dimension,
       parent_flat_index,
       parent_unflat_index,
       parent_center

"""
    CartesianParentGaussletBasis3D

Route-neutral identity object for a full Cartesian product gausslet parent.

The parent basis owns only the axis bases, the full parent lattice index
domain, Cartesian flattening/center lookup policy, axis sharing, and
scientific route metadata. Backend state, operator packets, residuals,
supplements, nested shell boxes, and working-box data belong to route-specific
objects, not this parent identity.
"""
struct CartesianParentGaussletBasis3D{A<:NamedTuple{(:x,:y,:z)},M}
    axes::A
    parent_box::NTuple{3,UnitRange{Int}}
    axis_sharing::Symbol
    metadata::M
end

function _parent_axis_signature(axis::MappedUniformBasis)
    return (
        family_name = axis.spec.family_value.name,
        reference_spacing = axis.spec.reference_spacing,
        reference_centers = axis.reference_center_data,
        centers = axis.center_data,
        integral_weights = axis.integral_weight_data,
        coefficient_matrix = axis.coefficient_matrix,
    )
end

function _parent_axis_sharing(axes::NamedTuple{(:x,:y,:z)})
    x_signature = _parent_axis_signature(axes.x)
    y_signature = _parent_axis_signature(axes.y)
    z_signature = _parent_axis_signature(axes.z)
    xy = isequal(x_signature, y_signature)
    xz = isequal(x_signature, z_signature)
    yz = isequal(y_signature, z_signature)
    xy && yz && return :shared_xyz
    xy && return :shared_xy
    xz && return :shared_xz
    yz && return :shared_yz
    return :separate_axes
end

function _parent_box_for_axes(axes::NamedTuple{(:x,:y,:z)})
    return (1:length(axes.x), 1:length(axes.y), 1:length(axes.z))
end

function _validate_parent_axes(axes::NamedTuple{(:x,:y,:z)})
    axes.x isa MappedUniformBasis || throw(ArgumentError("x parent axis must be a MappedUniformBasis"))
    axes.y isa MappedUniformBasis || throw(ArgumentError("y parent axis must be a MappedUniformBasis"))
    axes.z isa MappedUniformBasis || throw(ArgumentError("z parent axis must be a MappedUniformBasis"))
    return axes
end

function _validate_parent_box(
    axes::NamedTuple{(:x,:y,:z)},
    box::NTuple{3,UnitRange{Int}},
)
    expected = _parent_box_for_axes(axes)
    for axis_index in 1:3
        first(box[axis_index]) == 1 || throw(
            ArgumentError("parent_box ranges must start at 1 to preserve Cartesian parent indexing"),
        )
        length(box[axis_index]) == length(expected[axis_index]) || throw(
            ArgumentError("parent_box lengths must match parent axis lengths"),
        )
    end
    return box
end

function _axis_metadata(axis::MappedUniformBasis)
    return (
        basis_kind = :mapped_uniform,
        family_name = axis.spec.family_value.name,
        count = length(axis),
        reference_spacing = axis.spec.reference_spacing,
    )
end

function _default_parent_metadata(axes::NamedTuple{(:x,:y,:z)})
    return (
        basis_family = :mapped_uniform_cartesian_product,
        axis_metadata = (
            x = _axis_metadata(axes.x),
            y = _axis_metadata(axes.y),
            z = _axis_metadata(axes.z),
        ),
    )
end

function CartesianParentGaussletBasis3D(
    axes::NamedTuple{(:x,:y,:z)};
    parent_box::Union{Nothing,NTuple{3,UnitRange{Int}}} = nothing,
    metadata = _default_parent_metadata(axes),
)
    checked_axes = _validate_parent_axes(axes)
    checked_box = _validate_parent_box(
        checked_axes,
        parent_box === nothing ? _parent_box_for_axes(checked_axes) : parent_box,
    )
    return CartesianParentGaussletBasis3D(
        checked_axes,
        checked_box,
        _parent_axis_sharing(checked_axes),
        metadata,
    )
end

function CartesianParentGaussletBasis3D(
    basis_x::MappedUniformBasis,
    basis_y::MappedUniformBasis,
    basis_z::MappedUniformBasis;
    parent_box::Union{Nothing,NTuple{3,UnitRange{Int}}} = nothing,
    metadata = _default_parent_metadata((x = basis_x, y = basis_y, z = basis_z)),
)
    return CartesianParentGaussletBasis3D(
        (x = basis_x, y = basis_y, z = basis_z);
        parent_box,
        metadata,
    )
end

function CartesianParentGaussletBasis3D(basis::MappedUniformBasis)
    axes = (x = basis, y = basis, z = basis)
    return CartesianParentGaussletBasis3D(
        axes;
        metadata = (
            basis_family = :mapped_uniform_same_axis,
            axis_metadata = (
                x = _axis_metadata(basis),
                y = _axis_metadata(basis),
                z = _axis_metadata(basis),
            ),
        ),
    )
end

function CartesianParentGaussletBasis3D(basis::BondAlignedDiatomicQWBasis3D)
    return CartesianParentGaussletBasis3D(
        basis.basis_x,
        basis.basis_y,
        basis.basis_z;
        metadata = (
            basis_family = :bond_aligned_diatomic,
            bond_axis = basis.bond_axis,
            nuclei = copy(basis.nuclei),
            nuclear_charges = copy(basis.nuclear_charges),
            target_core_spacing = basis.target_core_spacing,
        ),
    )
end

function CartesianParentGaussletBasis3D(basis::BondAlignedHomonuclearChainQWBasis3D)
    return CartesianParentGaussletBasis3D(
        basis.basis_x,
        basis.basis_y,
        basis.basis_z;
        metadata = (
            basis_family = :bond_aligned_homonuclear_chain,
            chain_axis = basis.chain_axis,
            chain_coordinates = Float64[Float64(value) for value in basis.chain_coordinates],
            nuclei = copy(basis.nuclei),
            nuclear_charges = copy(basis.nuclear_charges),
            target_core_spacing = basis.target_core_spacing,
        ),
    )
end

function CartesianParentGaussletBasis3D(basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D)
    return CartesianParentGaussletBasis3D(
        basis.basis_x,
        basis.basis_y,
        basis.basis_z;
        metadata = (
            basis_family = :axis_aligned_homonuclear_square_lattice,
            lattice_size = basis.lattice_size,
            x_coordinates = Float64[Float64(value) for value in basis.x_coordinates],
            y_coordinates = Float64[Float64(value) for value in basis.y_coordinates],
            nuclei = copy(basis.nuclei),
            nuclear_charges = copy(basis.nuclear_charges),
            target_core_spacing = basis.target_core_spacing,
        ),
    )
end

CartesianParentGaussletBasis3D(fixed_block::_NestedFixedBlock3D) =
    CartesianParentGaussletBasis3D(fixed_block.parent_basis)

parent_axes(parent::CartesianParentGaussletBasis3D) = parent.axes
parent_box(parent::CartesianParentGaussletBasis3D) = parent.parent_box

function axis_basis(parent::CartesianParentGaussletBasis3D, axis::Symbol)
    axis == :x && return parent.axes.x
    axis == :y && return parent.axes.y
    axis == :z && return parent.axes.z
    throw(ArgumentError("axis must be :x, :y, or :z"))
end

function parent_axis_counts(parent::CartesianParentGaussletBasis3D)
    return (
        length(parent.parent_box[1]),
        length(parent.parent_box[2]),
        length(parent.parent_box[3]),
    )
end

parent_dimension(parent::CartesianParentGaussletBasis3D) = prod(parent_axis_counts(parent))

function parent_flat_index(
    parent::CartesianParentGaussletBasis3D,
    ix::Integer,
    iy::Integer,
    iz::Integer,
)
    x_index = Int(ix)
    y_index = Int(iy)
    z_index = Int(iz)
    box = parent.parent_box
    x_index in box[1] || throw(ArgumentError("x index must lie inside parent_box"))
    y_index in box[2] || throw(ArgumentError("y index must lie inside parent_box"))
    z_index in box[3] || throw(ArgumentError("z index must lie inside parent_box"))
    _nx, ny, nz = parent_axis_counts(parent)
    return (x_index - 1) * ny * nz + (y_index - 1) * nz + z_index
end

function parent_unflat_index(parent::CartesianParentGaussletBasis3D, index::Integer)
    flat_index = Int(index)
    nx, ny, nz = parent_axis_counts(parent)
    1 <= flat_index <= nx * ny * nz || throw(
        ArgumentError("flat parent index must lie inside parent_box"),
    )
    shifted = flat_index - 1
    plane = ny * nz
    ix = shifted ÷ plane + 1
    remainder = shifted % plane
    iy = remainder ÷ nz + 1
    iz = remainder % nz + 1
    return (ix, iy, iz)
end

function parent_center(
    parent::CartesianParentGaussletBasis3D,
    ix::Integer,
    iy::Integer,
    iz::Integer,
)
    state = (Int(ix), Int(iy), Int(iz))
    # Reuse parent_flat_index for bounds checking and to keep one ordering contract.
    parent_flat_index(parent, state...)
    return (
        Float64(centers(parent.axes.x)[state[1]]),
        Float64(centers(parent.axes.y)[state[2]]),
        Float64(centers(parent.axes.z)[state[3]]),
    )
end

parent_center(parent::CartesianParentGaussletBasis3D, state::NTuple{3,Int}) =
    parent_center(parent, state...)

end
