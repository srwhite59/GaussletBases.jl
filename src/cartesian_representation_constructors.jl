function _cartesian_basis_labels(prefix::AbstractString, count::Int)
    return [string(prefix, index) for index in 1:count]
end

function _cartesian_product_labels(
    orbitals::AbstractVector{<:CartesianProductOrbital3D},
)
    return [string("g(", orbital.ix, ",", orbital.iy, ",", orbital.iz, ")") for orbital in orbitals]
end

function _cartesian_labels(orbitals::AbstractVector{<:OrdinaryCartesianOrbital3D})
    return String[String(orbital.label) for orbital in orbitals]
end

function _cartesian_labels(orbitals::AbstractVector{<:_AtomicCartesianShellOrbital3D})
    return String[String(orbital.label) for orbital in orbitals]
end

function _cartesian_center_matrix(
    orbitals::AbstractVector{<:CartesianProductOrbital3D},
)
    centers_out = Matrix{Float64}(undef, length(orbitals), 3)
    for (row, orbital) in pairs(orbitals)
        centers_out[row, 1] = orbital.x
        centers_out[row, 2] = orbital.y
        centers_out[row, 3] = orbital.z
    end
    return centers_out
end

function _cartesian_center_matrix(
    orbitals::AbstractVector{<:OrdinaryCartesianOrbital3D},
)
    centers_out = Matrix{Float64}(undef, length(orbitals), 3)
    for (row, orbital) in pairs(orbitals)
        centers_out[row, 1] = orbital.x
        centers_out[row, 2] = orbital.y
        centers_out[row, 3] = orbital.z
    end
    return centers_out
end

function _cartesian_center_matrix(
    orbitals::AbstractVector{<:_AtomicCartesianShellOrbital3D},
)
    centers_out = Matrix{Float64}(undef, length(orbitals), 3)
    for (row, orbital) in pairs(orbitals)
        centers_out[row, 1] = orbital.center[1]
        centers_out[row, 2] = orbital.center[2]
        centers_out[row, 3] = orbital.center[3]
    end
    return centers_out
end

function _cartesian_axis_representation(basis::MappedUniformBasis)
    # The mapped ordinary Cartesian lane should expose the post-proxy,
    # post-cleanup working layer rather than the upstream distorted primitives.
    return _mapped_ordinary_working_basis_representation(basis)
end

function _cartesian_axis_representation(basis)
    return basis_representation(basis)
end

function _cartesian_axis_representations(basis::MappedUniformBasis)
    representation = _cartesian_axis_representation(basis)
    return (x = representation, y = representation, z = representation)
end

function _cartesian_axis_representations(basis::AbstractBondAlignedOrdinaryQWBasis3D)
    return (
        x = _cartesian_axis_representation(basis.basis_x),
        y = _cartesian_axis_representation(basis.basis_y),
        z = _cartesian_axis_representation(basis.basis_z),
    )
end

function _cartesian_axis_metadata(axis_representations::NamedTuple)
    return (
        x = axis_representations.x.metadata,
        y = axis_representations.y.metadata,
        z = axis_representations.z.metadata,
    )
end

function _cartesian_axis_signature(representation::BasisRepresentation1D)
    metadata = representation.metadata
    return (
        metadata.basis_kind,
        metadata.family_name,
        metadata.center_data,
        metadata.reference_center_data,
        metadata.integral_weight_data,
        metadata.basis_labels,
        metadata.coefficient_matrix,
    )
end

function _cartesian_axis_sharing(axis_representations::NamedTuple)
    x_signature = _cartesian_axis_signature(axis_representations.x)
    y_signature = _cartesian_axis_signature(axis_representations.y)
    z_signature = _cartesian_axis_signature(axis_representations.z)
    xy = isequal(x_signature, y_signature)
    xz = isequal(x_signature, z_signature)
    yz = isequal(y_signature, z_signature)
    xy && yz && return :shared_xyz
    xy && return :shared_xy
    xz && return :shared_xz
    yz && return :shared_yz
    return :separate_axes
end

function _cartesian_axis_counts(axis_representations::NamedTuple)
    return (
        length(axis_representations.x.metadata.center_data),
        length(axis_representations.y.metadata.center_data),
        length(axis_representations.z.metadata.center_data),
    )
end

function _cartesian_axis_centers(axis_representations::NamedTuple)
    return (
        axis_representations.x.metadata.center_data,
        axis_representations.y.metadata.center_data,
        axis_representations.z.metadata.center_data,
    )
end

function _cartesian_direct_product_representation(
    axis_representations::NamedTuple;
    route_metadata::NamedTuple = (;),
)
    axis_metadata = _cartesian_axis_metadata(axis_representations)
    axis_counts = _cartesian_axis_counts(axis_representations)
    x_centers, y_centers, z_centers = _cartesian_axis_centers(axis_representations)
    orbitals = _mapped_cartesian_orbitals(x_centers, y_centers, z_centers)
    labels = _cartesian_product_labels(orbitals)
    centers_matrix = _cartesian_center_matrix(orbitals)
    metadata = CartesianBasisMetadata3D(
        :direct_product,
        _cartesian_axis_sharing(axis_representations),
        axis_metadata,
        :cartesian_product_basis,
        axis_counts,
        length(orbitals),
        length(orbitals),
        nothing,
        labels,
        centers_matrix,
        route_metadata,
    )
    return CartesianBasisRepresentation3D(
        metadata,
        axis_representations,
        :identity,
        nothing,
        copy(labels),
        copy(centers_matrix),
        nothing,
        nothing,
        (;),
    )
end

function _cartesian_direct_product_representation(
    basis::MappedUniformBasis;
    route_metadata::NamedTuple = (;),
)
    return _cartesian_direct_product_representation(
        _cartesian_axis_representations(basis);
        route_metadata = route_metadata,
    )
end

function _cartesian_direct_product_representation(
    basis::AbstractBondAlignedOrdinaryQWBasis3D;
    route_metadata::NamedTuple = (;),
)
    return _cartesian_direct_product_representation(
        _cartesian_axis_representations(basis);
        route_metadata = route_metadata,
    )
end

function _cartesian_basis_route_metadata(basis::BondAlignedDiatomicQWBasis3D)
    return (
        basis_family = :bond_aligned_diatomic,
        bond_axis = basis.bond_axis,
        nuclei = copy(basis.nuclei),
        nuclear_charges = copy(basis.nuclear_charges),
        target_core_spacing = basis.target_core_spacing,
    )
end

function _cartesian_basis_route_metadata(basis::BondAlignedHomonuclearChainQWBasis3D)
    return (
        basis_family = :bond_aligned_homonuclear_chain,
        chain_axis = basis.chain_axis,
        chain_coordinates = Float64[Float64(value) for value in basis.chain_coordinates],
        nuclei = copy(basis.nuclei),
        nuclear_charges = copy(basis.nuclear_charges),
        target_core_spacing = basis.target_core_spacing,
    )
end

function _cartesian_basis_route_metadata(basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D)
    return (
        basis_family = :axis_aligned_homonuclear_square_lattice,
        lattice_size = basis.lattice_size,
        x_coordinates = Float64[Float64(value) for value in basis.x_coordinates],
        y_coordinates = Float64[Float64(value) for value in basis.y_coordinates],
        nuclei = copy(basis.nuclei),
        nuclear_charges = copy(basis.nuclear_charges),
        target_core_spacing = basis.target_core_spacing,
    )
end

function basis_representation(basis::BondAlignedDiatomicQWBasis3D)
    return _cartesian_direct_product_representation(
        basis;
        route_metadata = _cartesian_basis_route_metadata(basis),
    )
end

function basis_representation(basis::BondAlignedHomonuclearChainQWBasis3D)
    return _cartesian_direct_product_representation(
        basis;
        route_metadata = _cartesian_basis_route_metadata(basis),
    )
end

function basis_representation(basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D)
    return _cartesian_direct_product_representation(
        basis;
        route_metadata = _cartesian_basis_route_metadata(basis),
    )
end

basis_metadata(basis::AbstractBondAlignedOrdinaryQWBasis3D) = basis_representation(basis).metadata

function _cartesian_bond_aligned_build_metadata(
    basis::AbstractBondAlignedOrdinaryQWBasis3D,
)
    representation = basis_representation(basis)
    route_metadata = representation.metadata.route_metadata
    return (
        carried_representation = representation,
        parent_representation = representation,
        carried_metadata = representation.metadata,
        parent_metadata = representation.metadata,
        carried_route_metadata = route_metadata,
        parent_route_metadata = route_metadata,
        basis_family = route_metadata.basis_family,
    )
end

function _cartesian_parent_axis_representations(parent_basis)
    parent_basis isa MappedUniformBasis && return _cartesian_axis_representations(parent_basis)
    parent_basis isa AbstractBondAlignedOrdinaryQWBasis3D &&
        return _cartesian_axis_representations(parent_basis)
    throw(
        ArgumentError(
            "Cartesian basis representation does not yet support parent basis $(typeof(parent_basis))",
        ),
    )
end

function _cartesian_full_parent_box(axis_counts::NTuple{3,Int})
    return (1:axis_counts[1], 1:axis_counts[2], 1:axis_counts[3])
end

function _cartesian_shell_kind(::Any)
    return :unknown_shell
end

function _cartesian_shell_kind(::_CartesianNestedShell3D)
    return :shell
end

function _cartesian_shell_kind(::_CartesianNestedCompleteShell3D)
    return :complete_shell
end

function _cartesian_shell_kind(::_CartesianNestedShellPlusCore3D)
    return :shell_plus_core
end

function _cartesian_shell_kind(::_CartesianNestedShellSequence3D)
    return :shell_sequence
end

function _cartesian_shell_working_box(shell)
    return hasproperty(shell, :working_box) ? shell.working_box : nothing
end

function _cartesian_working_box_profile(
    working_box::Union{Nothing,NTuple{3,UnitRange{Int}}},
    axis_counts::NTuple{3,Int},
)
    working_box === nothing && return nothing
    return working_box == _cartesian_full_parent_box(axis_counts) ? :full_parent : :explicit_inner_box
end

function _cartesian_maybe_complete_shell_nside(shell)
    shell isa _CartesianNestedShellSequence3D || return nothing
    working_box = shell.working_box
    length(working_box[1]) == length(working_box[2]) == length(working_box[3]) || return nothing
    core_side = round(Int, cbrt(length(shell.core_indices)))
    core_side^3 == length(shell.core_indices) || return nothing
    expected_increment = _one_center_atomic_shell_increment(core_side)
    all(length(range) == expected_increment for range in shell.layer_column_ranges) || return nothing
    return core_side
end

function _cartesian_support_states(
    fixed_block::_NestedFixedBlock3D,
    axis_counts::NTuple{3,Int},
)
    if hasproperty(fixed_block.shell, :support_states)
        return NTuple{3,Int}[state for state in fixed_block.shell.support_states]
    end
    return NTuple{3,Int}[
        _cartesian_unflat_index(index, axis_counts) for index in fixed_block.support_indices
    ]
end

function _cartesian_fixed_block_route_metadata(
    fixed_block::_NestedFixedBlock3D,
    axis_counts::NTuple{3,Int},
)
    working_box = _cartesian_shell_working_box(fixed_block.shell)
    return (
        shell_kind = _cartesian_shell_kind(fixed_block.shell),
        working_box_profile = _cartesian_working_box_profile(working_box, axis_counts),
        nside = _cartesian_maybe_complete_shell_nside(fixed_block.shell),
        support_count = length(fixed_block.support_indices),
    )
end

function basis_representation(fixed_block::_NestedFixedBlock3D)
    axis_representations = _cartesian_parent_axis_representations(fixed_block.parent_basis)
    axis_metadata = _cartesian_axis_metadata(axis_representations)
    axis_counts = _cartesian_axis_counts(axis_representations)
    x_centers, y_centers, z_centers = _cartesian_axis_centers(axis_representations)
    parent_orbitals = _mapped_cartesian_orbitals(x_centers, y_centers, z_centers)
    parent_labels = _cartesian_product_labels(parent_orbitals)
    parent_centers = _cartesian_center_matrix(parent_orbitals)
    support_states = _cartesian_support_states(fixed_block, axis_counts)
    final_labels = _cartesian_basis_labels("nf", size(fixed_block.coefficient_matrix, 2))
    metadata = CartesianBasisMetadata3D(
        :nested_fixed_block,
        _cartesian_axis_sharing(axis_representations),
        axis_metadata,
        :cartesian_product_basis,
        axis_counts,
        length(parent_orbitals),
        size(fixed_block.coefficient_matrix, 2),
        _cartesian_shell_working_box(fixed_block.shell),
        final_labels,
        Matrix{Float64}(fixed_block.fixed_centers),
        _cartesian_fixed_block_route_metadata(fixed_block, axis_counts),
    )
    parent_data =
        isnothing(fixed_block.factorized_cartesian_parent_basis[]) ?
        (;) :
        (;
            factorized_cartesian_parent_basis =
                fixed_block.factorized_cartesian_parent_basis[],
        )
    return CartesianBasisRepresentation3D(
        metadata,
        axis_representations,
        :dense,
        _cartesian_coefficient_map_storage(fixed_block.coefficient_matrix),
        parent_labels,
        parent_centers,
        Vector{Int}(fixed_block.support_indices),
        support_states,
        parent_data,
    )
end

basis_metadata(fixed_block::_NestedFixedBlock3D) = basis_representation(fixed_block).metadata

function _cartesian_bond_aligned_build_metadata(
    fixed_block::_NestedFixedBlock3D{<:AbstractBondAlignedOrdinaryQWBasis3D},
)
    carried_representation = basis_representation(fixed_block)
    parent_representation = basis_representation(fixed_block.parent_basis)
    return (
        carried_representation = carried_representation,
        parent_representation = parent_representation,
        carried_metadata = carried_representation.metadata,
        parent_metadata = parent_representation.metadata,
        carried_route_metadata = carried_representation.metadata.route_metadata,
        parent_route_metadata = parent_representation.metadata.route_metadata,
        basis_family = parent_representation.metadata.route_metadata.basis_family,
    )
end
