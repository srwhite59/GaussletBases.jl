"""
    _CartesianNestedDoSide1D

First local one-dimensional nested-face contraction primitive on the finalized
QW-PGDG fixed line.

The object records the contracted side basis for one interval:

- the source interval in the finalized 1D fixed basis
- local overlap / position / weight / center data
- the retained local side basis in interval coordinates
- the same basis embedded back into the full 1D fixed line
- the localized side centers and signed local weights
"""
struct _CartesianNestedDoSide1D
    interval::UnitRange{Int}
    retained_count::Int
    local_overlap::Matrix{Float64}
    local_position::Matrix{Float64}
    local_weights::Vector{Float64}
    local_centers::Vector{Float64}
    local_coefficients::Matrix{Float64}
    coefficient_matrix::Matrix{Float64}
    localized_centers::Vector{Float64}
    localized_weights::Vector{Float64}
end

"""
    _CartesianNestedXYFace3D

First simple nested-face product object built from two tangential `doside`
spaces on one `x-y` face.

The face is attached to one fixed `z` index. The tangential directions are
contracted by local `doside` constructions and then combined as a product
space. The face support is restricted to the supplied tangential intervals and
the fixed `z` index, so opposite faces remain disjoint by construction.
"""
struct _CartesianNestedXYFace3D
    z_index::Int
    side_x::_CartesianNestedDoSide1D
    side_y::_CartesianNestedDoSide1D
    coefficient_matrix::Matrix{Float64}
    support_indices::Vector{Int}
end

abstract type _AbstractCartesianNestedShellLayer3D end

"""
    _CartesianNestedShellPacket3D

First shell-level transformed operator packet carried by the narrow nested
Cartesian face construction.

The packet stores the shell-basis transforms of the currently most important
fixed-block ingredients:

- overlap
- kinetic
- Cartesian position operators
- Cartesian second-moment operators
- contracted integral weights for the nested IDA transfer
- Gaussian-factor term packet
- pair-factor term packet
"""
struct _CartesianNestedShellPacket3D
    overlap::Matrix{Float64}
    kinetic::Matrix{Float64}
    position_x::Matrix{Float64}
    position_y::Matrix{Float64}
    position_z::Matrix{Float64}
    x2_x::Matrix{Float64}
    x2_y::Matrix{Float64}
    x2_z::Matrix{Float64}
    weights::Vector{Float64}
    gaussian_terms::Array{Float64,3}
    pair_terms::Array{Float64,3}
end

"""
    _CartesianNestedXYShell3D

First shell-level nested fixed-space object assembled from one opposite-face
pair of `x-y` faces.

This is the first shell object built on top of the local `doside` and
face-product primitives. It carries:

- the two opposite face objects
- the assembled shell contraction matrix
- the disjoint support rows in the parent Cartesian fixed block
- the first transformed shell-level operator packet
"""
struct _CartesianNestedXYShell3D
    faces::NTuple{2,_CartesianNestedXYFace3D}
    coefficient_matrix::Matrix{Float64}
    support_indices::Vector{Int}
    support_states::Vector{NTuple{3,Int}}
    packet::_CartesianNestedShellPacket3D
end

"""
    _CartesianNestedFace3D

Uniform shell-face object for the first generalized shell-packet interface.

Each face stores:

- the face kind `:xy`, `:xz`, or `:yz`
- the fixed axis and whether it is the low or high face
- the two tangential `doside` contractions
- the full parent-space coefficient matrix for that face piece
- the parent-space support rows for the face interior
"""
struct _CartesianNestedFace3D
    face_kind::Symbol
    fixed_axis::Symbol
    fixed_side::Symbol
    fixed_index::Int
    side_first::_CartesianNestedDoSide1D
    side_second::_CartesianNestedDoSide1D
    coefficient_matrix::Matrix{Float64}
    support_indices::Vector{Int}
end

"""
    _CartesianNestedShell3D

First generalized shell-packet interface built from a uniform collection of
shell faces.

The object is intended to look like a plausible future consumer input for the
existing Cartesian/QW-PGDG assembly style: one shell-level fixed basis plus one
transformed shell-level packet carrying the same operator ingredients as the
current fixed block.
"""
struct _CartesianNestedShell3D <: _AbstractCartesianNestedShellLayer3D
    faces::Vector{_CartesianNestedFace3D}
    face_column_ranges::Vector{UnitRange{Int}}
    coefficient_matrix::Matrix{Float64}
    support_indices::Vector{Int}
    support_states::Vector{NTuple{3,Int}}
    packet::_CartesianNestedShellPacket3D
end

"""
    _CartesianNestedEdge3D

First codimension-2 shell-edge primitive for the complete shell-layer
decomposition.

Each edge stores:

- the free axis along the edge
- the two fixed axes and whether they are on the low or high boundary
- one one-dimensional `doside` contraction on the free interval
- the parent-space coefficient matrix for that edge piece
- the parent-space support rows for the open edge
"""
struct _CartesianNestedEdge3D
    free_axis::Symbol
    fixed_axes::NTuple{2,Symbol}
    fixed_sides::NTuple{2,Symbol}
    fixed_indices::NTuple{2,Int}
    side::_CartesianNestedDoSide1D
    coefficient_matrix::Matrix{Float64}
    support_indices::Vector{Int}
end

"""
    _CartesianNestedCorner3D

First codimension-3 shell-corner primitive for the complete shell-layer
decomposition.

The first pass keeps corners as direct retained pieces.
"""
struct _CartesianNestedCorner3D
    fixed_sides::NTuple{3,Symbol}
    fixed_indices::NTuple{3,Int}
    coefficient_matrix::Matrix{Float64}
    support_indices::Vector{Int}
end

"""
    _CartesianNestedCompleteShell3D

First complete nonrecursive shell-layer object.

This augments the earlier face-only shell language with explicit edge and
corner pieces so the shell annulus is partitioned without leftovers.
"""
struct _CartesianNestedCompleteShell3D <: _AbstractCartesianNestedShellLayer3D
    faces::Vector{_CartesianNestedFace3D}
    face_column_ranges::Vector{UnitRange{Int}}
    edges::Vector{_CartesianNestedEdge3D}
    edge_column_ranges::Vector{UnitRange{Int}}
    corners::Vector{_CartesianNestedCorner3D}
    corner_column_ranges::Vector{UnitRange{Int}}
    coefficient_matrix::Matrix{Float64}
    support_indices::Vector{Int}
    support_states::Vector{NTuple{3,Int}}
    packet::_CartesianNestedShellPacket3D
end

"""
    _CartesianNestedShellPlusCore3D

First nonrecursive shell-plus-core fixed-space object.

This augments the existing shell-face packet with a direct interior block taken
from the parent fixed basis itself. The shell faces keep their current disjoint
face-interior role, while the core block fills the missing interior volume.
"""
struct _CartesianNestedShellPlusCore3D
    shell::_CartesianNestedShell3D
    core_indices::Vector{Int}
    core_states::Vector{NTuple{3,Int}}
    core_column_range::UnitRange{Int}
    shell_column_ranges::Vector{UnitRange{Int}}
    coefficient_matrix::Matrix{Float64}
    support_indices::Vector{Int}
    support_states::Vector{NTuple{3,Int}}
    packet::_CartesianNestedShellPacket3D
end

"""
    _CartesianNestedShellSequence3D

First nonrecursive multi-shell fixed-space source in the same fixed-block
language as the existing shell-plus-core adapter path.

The object keeps one retained interior core block and an ordered list of shell
layers, together with one combined coefficient matrix and one propagated packet
that downstream fixed-block consumers can read without any new consumer logic.
"""
struct _CartesianNestedShellSequence3D{S<:_AbstractCartesianNestedShellLayer3D}
    core_indices::Vector{Int}
    core_states::Vector{NTuple{3,Int}}
    core_column_range::UnitRange{Int}
    shell_layers::Vector{S}
    layer_column_ranges::Vector{UnitRange{Int}}
    working_box::NTuple{3,UnitRange{Int}}
    coefficient_matrix::Matrix{Float64}
    support_indices::Vector{Int}
    support_states::Vector{NTuple{3,Int}}
    packet::_CartesianNestedShellPacket3D
end

"""
    _NestedFixedBlock3D

Adapter object exposing one generalized nested shell packet in the same
fixed-block language that a downstream Cartesian/QW-PGDG consumer can read.

The object keeps:

- the parent mapped basis used to define the raw fixed-to-Gaussian blocks
- the shell-level contraction map
- the propagated fixed-fixed packet on the nested shell basis
- the contracted fixed-block integral weights used by the IDA interaction
  representation
- simple shell-center metadata for nearest/GGT diagnostics
"""
struct _NestedFixedBlock3D{B,S}
    parent_basis::B
    shell::S
    coefficient_matrix::Matrix{Float64}
    support_indices::Vector{Int}
    overlap::Matrix{Float64}
    kinetic::Matrix{Float64}
    position_x::Matrix{Float64}
    position_y::Matrix{Float64}
    position_z::Matrix{Float64}
    x2_x::Matrix{Float64}
    x2_y::Matrix{Float64}
    x2_z::Matrix{Float64}
    weights::Vector{Float64}
    gaussian_terms::Array{Float64,3}
    pair_terms::Array{Float64,3}
    fixed_centers::Matrix{Float64}
end

"""
    _CartesianNestedAxisBundles3D

Narrow mixed-axis parent bundle container for the first bond-aligned diatomic
nested fixed-block route.

The atomic shell language already assumes localized 1D PGDG data on each axis.
This object lifts that assumption to three explicit axes so the shell language
can be reused on a rectangular parent box with unequal axis lengths.
"""
struct _CartesianNestedAxisBundles3D{BX,BY,BZ}
    bundle_x::BX
    bundle_y::BY
    bundle_z::BZ
end

"""
    _BondAlignedDiatomicSplitGeometry3D

Geometry report for the first bond-aligned diatomic split/no-split decision.

The split is attached to the original parent grid:

- `working_box` is the shared box remaining after the outer shared shell stage
- `split_index` is the bond-axis parent-grid index nearest the bond midpoint
- `shared_midpoint_box` is the direct shared midpoint slab for the odd-length
  homonuclear case
- `child_boxes` are the two nonoverlapping child boxes if the split is allowed
- `child_physical_widths` records the mapped physical widths of those children
"""
struct _BondAlignedDiatomicSplitGeometry3D
    parent_box::NTuple{3,UnitRange{Int}}
    working_box::NTuple{3,UnitRange{Int}}
    bond_axis::Symbol
    midpoint::Float64
    split_index::Int
    count_eligible::Bool
    shape_eligible::Bool
    did_split::Bool
    shared_midpoint_box::Union{Nothing,NTuple{3,UnitRange{Int}}}
    child_boxes::Vector{NTuple{3,UnitRange{Int}}}
    child_physical_widths::Vector{NTuple{3,Float64}}
end

struct CartesianNestedCompleteShellRetentionContract
    nside::Int
    retain_xy::NTuple{2,Int}
    retain_xz::NTuple{2,Int}
    retain_yz::NTuple{2,Int}
    retain_x_edge::Int
    retain_y_edge::Int
    retain_z_edge::Int
    face_retained_count::Int
    edge_retained_count::Int
    corner_retained_count::Int
    shell_increment::Int
    matches_nside_default::Bool
end

struct CartesianNestedSequenceContractAudit
    parent_dims::NTuple{3,Int}
    working_box::NTuple{3,UnitRange{Int}}
    full_parent_working_box::Bool
    support_count::Int
    expected_support_count::Int
    missing_row_count::Int
    ownership_group_count_min::Int
    ownership_group_count_max::Int
    ownership_unowned_row_count::Int
    ownership_multi_owned_row_count::Int
end

struct NestedFixedBlockBuildTimingSummary
    records::Vector{Pair{String,Float64}}
end

struct TimedNestedFixedBlockBuild{F}
    fixed_block::F
    timings::NestedFixedBlockBuildTimingSummary
end

mutable struct _NestedFixedBlockTimingCollector
    records::Vector{Pair{String,Float64}}
end

struct _CartesianNestedSupportAxes3D
    x::Vector{Int}
    y::Vector{Int}
    z::Vector{Int}
end

struct _CartesianNestedFactorizedBasis3D
    dims::NTuple{3,Int}
    x_functions::Matrix{Float64}
    y_functions::Matrix{Float64}
    z_functions::Matrix{Float64}
    basis_triplets::Vector{NTuple{3,Int}}
    basis_amplitudes::Vector{Float64}
    reconstruction_max_error::Float64
end

struct _CartesianNestedFactorizedAxisBaseTables
    overlap::Matrix{Float64}
    kinetic::Matrix{Float64}
    position::Matrix{Float64}
    x2::Matrix{Float64}
end

"""
    _CartesianNestedBondAlignedDiatomicSource3D

First bond-aligned diatomic nested fixed-space source built on top of the
existing atomic shell language.

The source keeps:

- the mixed-axis parent bundle data
- the shared-box split/no-split geometry decision
- the outer shared shell layers
- the child atomic-style subtrees after the bond-axis split
- the merged shell-sequence object used to build the fixed block
"""
struct _CartesianNestedBondAlignedDiatomicSource3D{B}
    basis::B
    axis_bundles::_CartesianNestedAxisBundles3D
    nside::Int
    child_shell_retention_contract::CartesianNestedCompleteShellRetentionContract
    shared_shell_retention_contract::CartesianNestedCompleteShellRetentionContract
    geometry::_BondAlignedDiatomicSplitGeometry3D
    shared_shell_layers::Vector{_CartesianNestedCompleteShell3D}
    child_sequences::Vector{_CartesianNestedShellSequence3D}
    child_column_ranges::Vector{UnitRange{Int}}
    midpoint_slab_column_range::Union{Nothing,UnitRange{Int}}
    sequence::_CartesianNestedShellSequence3D
end

"""
    _BondAlignedHomonuclearChainSplitCandidate3D

Experimental split candidate record for one homonuclear chain subtree.

This keeps the candidate midpoint policy explicit rather than burying it in the
tree builder.
"""
struct _BondAlignedHomonuclearChainSplitCandidate3D
    split_kind::Symbol
    nucleus_ranges::Vector{UnitRange{Int}}
    midpoint_values::Vector{Float64}
    split_indices::Vector{Int}
    child_boxes::Vector{NTuple{3,UnitRange{Int}}}
    child_physical_widths::Vector{NTuple{3,Float64}}
    child_parallel_counts::Vector{Int}
    child_parallel_to_transverse_ratios::Vector{Float64}
    count_eligible::Bool
    shape_eligible::Bool
    did_split::Bool
end

"""
    _BondAlignedHomonuclearChainNodeGeometry3D

Experimental recursive split-tree geometry record for the first homonuclear
chain nesting pass.
"""
struct _BondAlignedHomonuclearChainNodeGeometry3D
    node_label::String
    parent_box::NTuple{3,UnitRange{Int}}
    working_box::NTuple{3,UnitRange{Int}}
    chain_axis::Symbol
    odd_chain_policy::Symbol
    nside::Int
    min_parallel_to_transverse_ratio::Float64
    nucleus_range::UnitRange{Int}
    chain_coordinates::Vector{Float64}
    shared_shell_count::Int
    shared_shell_dimensions::Vector{Int}
    candidate_splits::Vector{_BondAlignedHomonuclearChainSplitCandidate3D}
    accepted_candidate_index::Union{Nothing,Int}
    local_resolution_warning::Bool
    child_nodes::Vector{_BondAlignedHomonuclearChainNodeGeometry3D}
    subtree_fixed_dimension::Int
end

"""
    _CartesianNestedBondAlignedHomonuclearChainSource3D

Experimental fixed-source container for the first homonuclear chain nested
geometry line.
"""
struct _CartesianNestedBondAlignedHomonuclearChainSource3D{B}
    basis::B
    axis_bundles::_CartesianNestedAxisBundles3D
    shell_retention_contract::CartesianNestedCompleteShellRetentionContract
    root_geometry::_BondAlignedHomonuclearChainNodeGeometry3D
    leaf_sequences::Vector{_CartesianNestedShellSequence3D}
    sequence::_CartesianNestedShellSequence3D
end

"""
    _AxisAlignedHomonuclearSquareLatticeSplitCandidate3D

Experimental planar split candidate record for one square-lattice subtree.

The first pass keeps the candidate family explicit:

- strip splits along `x`
- strip splits along `y`
- binary for even site counts
- centered ternary for odd site counts
"""
struct _AxisAlignedHomonuclearSquareLatticeSplitCandidate3D
    split_family::Symbol
    split_axis::Symbol
    x_coordinate_ranges::Vector{UnitRange{Int}}
    y_coordinate_ranges::Vector{UnitRange{Int}}
    split_values::Vector{Float64}
    split_indices::Vector{Int}
    child_boxes::Vector{NTuple{3,UnitRange{Int}}}
    child_physical_widths::Vector{NTuple{3,Float64}}
    child_planar_counts::Vector{NTuple{2,Int}}
    child_in_plane_aspect_ratios::Vector{Float64}
    count_eligible::Bool
    shape_eligible::Bool
    symmetry_preserving::Bool
    did_split::Bool
end

"""
    _AxisAlignedHomonuclearSquareLatticeNodeGeometry3D

Experimental recursive split-tree geometry record for the first homonuclear
square-lattice planar nesting pass.
"""
struct _AxisAlignedHomonuclearSquareLatticeNodeGeometry3D
    node_label::String
    parent_box::NTuple{3,UnitRange{Int}}
    working_box::NTuple{3,UnitRange{Int}}
    lattice_size::Int
    nside::Int
    min_in_plane_aspect_ratio::Float64
    x_coordinate_range::UnitRange{Int}
    y_coordinate_range::UnitRange{Int}
    x_coordinates::Vector{Float64}
    y_coordinates::Vector{Float64}
    shared_shell_count::Int
    shared_shell_dimensions::Vector{Int}
    candidate_splits::Vector{_AxisAlignedHomonuclearSquareLatticeSplitCandidate3D}
    accepted_candidate_index::Union{Nothing,Int}
    local_resolution_warning::Bool
    child_nodes::Vector{_AxisAlignedHomonuclearSquareLatticeNodeGeometry3D}
    subtree_fixed_dimension::Int
end

"""
    _CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D

Experimental fixed-source container for the first homonuclear square-lattice
planar geometry line.
"""
struct _CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D{B}
    basis::B
    axis_bundles::_CartesianNestedAxisBundles3D
    shell_retention_contract::CartesianNestedCompleteShellRetentionContract
    root_geometry::_AxisAlignedHomonuclearSquareLatticeNodeGeometry3D
    leaf_sequences::Vector{_CartesianNestedShellSequence3D}
    sequence::_CartesianNestedShellSequence3D
end

"""
    _CartesianNestedDoSideTrace1D

Structured diagnostic record for one local 1D `doside` / COMX contraction
used in the current nested Cartesian source language.

This keeps enough information to make the local center-loss question explicit:

- where the contraction was used
- which parent interval it came from
- the physical centers before contraction
- the localized centers and signed weights after COMX cleanup
- whether the parent interval is symmetric about zero
- whether the localized centers retain a near-zero center
"""
struct _CartesianNestedDoSideTrace1D
    context_label::String
    group_kind::Symbol
    layer_index::Int
    piece_kind::Symbol
    axis::Symbol
    usage_label::String
    interval::UnitRange{Int}
    parent_centers::Vector{Float64}
    retained_count::Int
    localized_centers::Vector{Float64}
    localized_weights::Vector{Float64}
    symmetric_about_zero::Bool
    symmetry_error::Float64
    contains_near_zero_center::Bool
    even_retained_count::Bool
end

function Base.show(io::IO, side::_CartesianNestedDoSide1D)
    print(
        io,
        "_CartesianNestedDoSide1D(interval=",
        side.interval,
        ", retained=",
        side.retained_count,
        ")",
    )
end

function Base.show(io::IO, face::_CartesianNestedXYFace3D)
    print(
        io,
        "_CartesianNestedXYFace3D(z_index=",
        face.z_index,
        ", nx=",
        size(face.side_x.coefficient_matrix, 2),
        ", ny=",
        size(face.side_y.coefficient_matrix, 2),
        ")",
    )
end

function Base.show(io::IO, shell::_CartesianNestedXYShell3D)
    print(
        io,
        "_CartesianNestedXYShell3D(nfaces=2, nshell=",
        size(shell.coefficient_matrix, 2),
        ", nsupport=",
        length(shell.support_indices),
        ")",
    )
end

function Base.show(io::IO, face::_CartesianNestedFace3D)
    print(
        io,
        "_CartesianNestedFace3D(",
        face.face_kind,
        ", ",
        face.fixed_side,
        ", ncols=",
        size(face.coefficient_matrix, 2),
        ")",
    )
end

function Base.show(io::IO, shell::_CartesianNestedShell3D)
    print(
        io,
        "_CartesianNestedShell3D(nfaces=",
        length(shell.faces),
        ", nshell=",
        size(shell.coefficient_matrix, 2),
        ", nsupport=",
        length(shell.support_indices),
        ")",
    )
end

function Base.show(io::IO, edge::_CartesianNestedEdge3D)
    print(
        io,
        "_CartesianNestedEdge3D(",
        edge.free_axis,
        ", ",
        edge.fixed_sides,
        ", ncols=",
        size(edge.coefficient_matrix, 2),
        ")",
    )
end

function Base.show(io::IO, corner::_CartesianNestedCorner3D)
    print(
        io,
        "_CartesianNestedCorner3D(indices=",
        corner.fixed_indices,
        ")",
    )
end

function Base.show(io::IO, shell::_CartesianNestedCompleteShell3D)
    print(
        io,
        "_CartesianNestedCompleteShell3D(nfaces=",
        length(shell.faces),
        ", nedges=",
        length(shell.edges),
        ", ncorners=",
        length(shell.corners),
        ", nshell=",
        size(shell.coefficient_matrix, 2),
        ")",
    )
end

function Base.show(io::IO, shell::_CartesianNestedShellPlusCore3D)
    print(
        io,
        "_CartesianNestedShellPlusCore3D(ncore=",
        length(shell.core_indices),
        ", nfaces=",
        length(shell.shell.faces),
        ", nshell=",
        size(shell.coefficient_matrix, 2),
        ")",
    )
end

function Base.show(io::IO, shell::_CartesianNestedShellSequence3D)
    print(
        io,
        "_CartesianNestedShellSequence3D(ncore=",
        length(shell.core_indices),
        ", nlayers=",
        length(shell.shell_layers),
        ", nshell=",
        size(shell.coefficient_matrix, 2),
        ")",
    )
end

function Base.show(io::IO, fixed_block::_NestedFixedBlock3D)
    print(
        io,
        "_NestedFixedBlock3D(nfixed=",
        size(fixed_block.overlap, 1),
        ", nsupport=",
        length(fixed_block.support_indices),
        ")",
    )
end

function Base.show(io::IO, bundles::_CartesianNestedAxisBundles3D)
    dims = _nested_axis_lengths(bundles)
    print(
        io,
        "_CartesianNestedAxisBundles3D(nx=",
        dims[1],
        ", ny=",
        dims[2],
        ", nz=",
        dims[3],
        ")",
    )
end

function Base.show(io::IO, geometry::_BondAlignedDiatomicSplitGeometry3D)
    print(
        io,
        "_BondAlignedDiatomicSplitGeometry3D(axis=:",
        geometry.bond_axis,
        ", working_box=",
        geometry.working_box,
        ", split_index=",
        geometry.split_index,
        ", midpoint_slab=",
        isnothing(geometry.shared_midpoint_box) ? "nothing" : geometry.shared_midpoint_box,
        ", did_split=",
        geometry.did_split,
        ")",
    )
end

function Base.show(io::IO, source::_CartesianNestedBondAlignedDiatomicSource3D)
    print(
        io,
        "_CartesianNestedBondAlignedDiatomicSource3D(nshared=",
        length(source.shared_shell_layers),
        ", nchild=",
        length(source.child_sequences),
        ", nfixed=",
        size(source.sequence.coefficient_matrix, 2),
        ", nside=",
        source.nside,
        ", midpoint_slab=",
        !isnothing(source.midpoint_slab_column_range),
        ", did_split=",
        source.geometry.did_split,
        ")",
    )
end

function Base.show(io::IO, candidate::_BondAlignedHomonuclearChainSplitCandidate3D)
    print(
        io,
        "_BondAlignedHomonuclearChainSplitCandidate3D(kind=:",
        candidate.split_kind,
        ", nchild=",
        length(candidate.child_boxes),
        ", did_split=",
        candidate.did_split,
        ")",
    )
end

function Base.show(io::IO, node::_BondAlignedHomonuclearChainNodeGeometry3D)
    print(
        io,
        "_BondAlignedHomonuclearChainNodeGeometry3D(node=",
        node.node_label,
        ", nuclei=",
        node.nucleus_range,
        ", nchild=",
        length(node.child_nodes),
        ", did_split=",
        !isnothing(node.accepted_candidate_index),
        ", nfixed=",
        node.subtree_fixed_dimension,
        ")",
    )
end

function Base.show(io::IO, source::_CartesianNestedBondAlignedHomonuclearChainSource3D)
    print(
        io,
        "_CartesianNestedBondAlignedHomonuclearChainSource3D(nleaf=",
        length(source.leaf_sequences),
        ", nfixed=",
        size(source.sequence.coefficient_matrix, 2),
        ", nside=",
        source.shell_retention_contract.nside,
        ", did_split=",
        !isnothing(source.root_geometry.accepted_candidate_index),
        ")",
    )
end

function Base.show(io::IO, candidate::_AxisAlignedHomonuclearSquareLatticeSplitCandidate3D)
    print(
        io,
        "_AxisAlignedHomonuclearSquareLatticeSplitCandidate3D(family=:",
        candidate.split_family,
        ", nchild=",
        length(candidate.child_boxes),
        ", did_split=",
        candidate.did_split,
        ")",
    )
end

function Base.show(io::IO, node::_AxisAlignedHomonuclearSquareLatticeNodeGeometry3D)
    print(
        io,
        "_AxisAlignedHomonuclearSquareLatticeNodeGeometry3D(node=",
        node.node_label,
        ", x_range=",
        node.x_coordinate_range,
        ", y_range=",
        node.y_coordinate_range,
        ", nchild=",
        length(node.child_nodes),
        ", did_split=",
        !isnothing(node.accepted_candidate_index),
        ", nfixed=",
        node.subtree_fixed_dimension,
        ")",
    )
end

function Base.show(io::IO, source::_CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D)
    print(
        io,
        "_CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D(nleaf=",
        length(source.leaf_sequences),
        ", nfixed=",
        size(source.sequence.coefficient_matrix, 2),
        ", nside=",
        source.shell_retention_contract.nside,
        ", did_split=",
        !isnothing(source.root_geometry.accepted_candidate_index),
        ")",
    )
end

function Base.show(io::IO, timings::NestedFixedBlockBuildTimingSummary)
    total = nested_fixed_block_timing_seconds(timings, "fixed_block.total")
    print(
        io,
        "NestedFixedBlockBuildTimingSummary(ntimings=",
        length(timings.records),
        ", total=",
        total,
        "s)",
    )
end

function Base.show(io::IO, timed::TimedNestedFixedBlockBuild)
    print(
        io,
        "TimedNestedFixedBlockBuild(nfixed=",
        size(timed.fixed_block.overlap, 1),
        ", total=",
        nested_fixed_block_timing_seconds(timed.timings, "fixed_block.total"),
        "s)",
    )
end

function Base.show(io::IO, trace::_CartesianNestedDoSideTrace1D)
    print(
        io,
        "_CartesianNestedDoSideTrace1D(context=",
        trace.context_label,
        ", axis=:",
        trace.axis,
        ", retained=",
        trace.retained_count,
        ", near_zero=",
        trace.contains_near_zero_center,
        ")",
    )
end

function _nested_axis_bundle(
    bundles::_CartesianNestedAxisBundles3D,
    axis::Symbol,
)
    axis == :x && return bundles.bundle_x
    axis == :y && return bundles.bundle_y
    axis == :z && return bundles.bundle_z
    throw(ArgumentError("nested axis bundle lookup requires axis :x, :y, or :z"))
end

function _nested_axis_pgdg(
    bundles::_CartesianNestedAxisBundles3D,
    axis::Symbol,
)
    return _nested_axis_bundle(bundles, axis).pgdg_intermediate
end

function _nested_axis_lengths(bundles::_CartesianNestedAxisBundles3D)
    return (
        size(_nested_axis_pgdg(bundles, :x).overlap, 1),
        size(_nested_axis_pgdg(bundles, :y).overlap, 1),
        size(_nested_axis_pgdg(bundles, :z).overlap, 1),
    )
end

function _nested_elapsed_seconds(start_ns::UInt64)
    return (time_ns() - start_ns) / 1.0e9
end

function _nested_timing_enabled(timing)
    return timing === true || timing === :report
end

function _nested_new_timing_collector()
    return _NestedFixedBlockTimingCollector(Pair{String,Float64}[])
end

function _nested_record_timing!(
    collector::Union{Nothing,_NestedFixedBlockTimingCollector},
    label::AbstractString,
    start_ns::UInt64,
)
    isnothing(collector) && return nothing
    push!(collector.records, String(label) => _nested_elapsed_seconds(start_ns))
    return nothing
end

function _nested_timing_summary(
    collector::Union{Nothing,_NestedFixedBlockTimingCollector},
)
    return NestedFixedBlockBuildTimingSummary(
        isnothing(collector) ? Pair{String,Float64}[] : copy(collector.records),
    )
end

function nested_fixed_block_timing_seconds(
    timings::NestedFixedBlockBuildTimingSummary,
    label::AbstractString,
)
    return sum(record.second for record in timings.records if record.first == label)
end

function nested_fixed_block_timing_report(
    io::IO,
    timings::NestedFixedBlockBuildTimingSummary,
)
    println(io, "Nested fixed-block timings")
    labels = (
        "fixed_block.total",
        "fixed_block.parent_bundle",
        "fixed_block.sequence_build",
        "fixed_block.adapter",
        "shell_layer.nonpacket",
        "sequence_merge.nonpacket",
        "packet.setup",
        "packet.total",
        "packet.base.overlap",
        "packet.base.kinetic",
        "packet.base.position_x",
        "packet.base.position_y",
        "packet.base.position_z",
        "packet.base.x2_x",
        "packet.base.x2_y",
        "packet.base.x2_z",
        "packet.gaussian_terms",
        "packet.pair_terms",
    )
    for label in labels
        seconds = nested_fixed_block_timing_seconds(timings, label)
        seconds > 0.0 || continue
        println(io, "  ", rpad(label, 28), " ", string(round(seconds; digits = 6)), " s")
    end
    return timings
end

function nested_fixed_block_timing_report(
    timings::NestedFixedBlockBuildTimingSummary,
)
    return sprint(io -> nested_fixed_block_timing_report(io, timings))
end

function nested_fixed_block_timing_report(
    io::IO,
    timed::TimedNestedFixedBlockBuild,
)
    return nested_fixed_block_timing_report(io, timed.timings)
end

function nested_fixed_block_timing_report(
    timed::TimedNestedFixedBlockBuild,
)
    return nested_fixed_block_timing_report(timed.timings)
end

function _nested_metric_norm(
    vector::AbstractVector{<:Real},
    overlap::AbstractMatrix{<:Real},
)
    value = Float64(dot(vector, overlap * vector))
    return sqrt(abs(value))
end

function _nested_metric_normalize(
    vector::AbstractVector{<:Real},
    overlap::AbstractMatrix{<:Real};
    tol::Float64 = 1.0e-12,
)
    norm_value = _nested_metric_norm(vector, overlap)
    norm_value > tol || throw(ArgumentError("nested doside construction encountered a near-null local direction"))
    return Float64.(vector) ./ norm_value
end

function _nested_metric_orthogonalize(
    vector::AbstractVector{<:Real},
    basis::AbstractMatrix{<:Real},
    overlap::AbstractMatrix{<:Real},
)
    size(basis, 2) == 0 && return Float64.(vector)
    projected = Float64.(vector) .- Matrix{Float64}(basis) * vec(transpose(basis) * (overlap * vector))
    projected .-= Matrix{Float64}(basis) * vec(transpose(basis) * (overlap * projected))
    return projected
end

function _nested_standard_basis_direction(
    overlap::AbstractMatrix{<:Real},
    basis::AbstractMatrix{<:Real},
)
    nlocal = size(overlap, 1)
    for index in 1:nlocal
        candidate = zeros(Float64, nlocal)
        candidate[index] = 1.0
        candidate = _nested_metric_orthogonalize(candidate, basis, overlap)
        _nested_metric_norm(candidate, overlap) > 1.0e-12 && return candidate
    end
    throw(ArgumentError("nested doside construction could not find an additional independent local direction"))
end

function _nested_retained_span(
    local_weights::AbstractVector{<:Real},
    local_centers::AbstractVector{<:Real},
    local_position::AbstractMatrix{<:Real},
    local_overlap::AbstractMatrix{<:Real},
    retained_count::Int,
)
    nlocal = length(local_weights)
    retained_count >= 1 || throw(ArgumentError("nested doside construction requires retained_count >= 1"))
    retained_count <= nlocal || throw(ArgumentError("nested doside retained_count must not exceed the interval size"))

    start_vector = Float64.(local_weights)
    if _nested_metric_norm(start_vector, local_overlap) <= 1.0e-12
        start_vector = sqrt.(abs.(Float64.(local_weights)))
    end
    if _nested_metric_norm(start_vector, local_overlap) <= 1.0e-12
        start_vector .= 1.0
    end

    raw_basis = zeros(Float64, nlocal, retained_count)
    raw_basis[:, 1] .= _nested_metric_normalize(start_vector, local_overlap)

    for column in 2:retained_count
        previous = view(raw_basis, :, column - 1)
        candidate = Float64.(local_position * previous)
        candidate = _nested_metric_orthogonalize(candidate, view(raw_basis, :, 1:(column - 1)), local_overlap)
        if _nested_metric_norm(candidate, local_overlap) <= 1.0e-12
            power_direction = (Float64.(local_centers) .^ (column - 1)) .* start_vector
            candidate = _nested_metric_orthogonalize(power_direction, view(raw_basis, :, 1:(column - 1)), local_overlap)
        end
        if _nested_metric_norm(candidate, local_overlap) <= 1.0e-12
            candidate = _nested_standard_basis_direction(local_overlap, view(raw_basis, :, 1:(column - 1)))
        end
        raw_basis[:, column] .= _nested_metric_normalize(candidate, local_overlap)
    end
    return raw_basis
end

function _nested_interval_data(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    interval::UnitRange{Int},
)
    n1d = size(pgdg.overlap, 1)
    first(interval) >= 1 || throw(ArgumentError("nested doside interval must start inside the finalized fixed line"))
    last(interval) <= n1d || throw(ArgumentError("nested doside interval must end inside the finalized fixed line"))
    local_overlap = Matrix{Float64}(pgdg.overlap[interval, interval])
    local_position = Matrix{Float64}(pgdg.position[interval, interval])
    local_weights = Float64[pgdg.weights[index] for index in interval]
    local_centers = Float64[pgdg.centers[index] for index in interval]
    return (
        overlap = local_overlap,
        position = local_position,
        weights = local_weights,
        centers = local_centers,
        n1d = n1d,
    )
end

function _nested_zero_symmetry_error(values::AbstractVector{<:Real})
    isempty(values) && return 0.0
    errors = Float64[
        abs(Float64(values[index]) + Float64(values[end - index + 1])) for
        index in 1:length(values)
    ]
    return maximum(errors)
end

function _nested_is_symmetric_about_zero(
    values::AbstractVector{<:Real};
    tol::Float64 = 1.0e-8,
)
    error = _nested_zero_symmetry_error(values)
    return error <= tol, error
end

function _nested_contains_near_zero(
    values::AbstractVector{<:Real};
    tol::Float64 = 1.0e-8,
)
    return any(abs(Float64(value)) <= tol for value in values)
end

function _nested_doside_retained_count(
    local_centers::AbstractVector{<:Real},
    provisional_retained_count::Int,
)
    provisional_retained_count >= 1 || throw(
        ArgumentError("nested doside retained local count must be at least 1"),
    )
    symmetric_about_zero, _ = _nested_is_symmetric_about_zero(local_centers)
    if symmetric_about_zero && iseven(provisional_retained_count) && provisional_retained_count > 1
        return provisional_retained_count - 1
    end
    return provisional_retained_count
end

function _nested_doside_trace(
    side::_CartesianNestedDoSide1D;
    context_label::AbstractString,
    group_kind::Symbol,
    layer_index::Integer,
    piece_kind::Symbol,
    axis::Symbol,
    usage_label::AbstractString,
    symmetry_tol::Float64 = 1.0e-8,
    zero_tol::Float64 = 1.0e-8,
)
    symmetric_about_zero, symmetry_error = _nested_is_symmetric_about_zero(
        side.local_centers;
        tol = symmetry_tol,
    )
    contains_near_zero_center = _nested_contains_near_zero(
        side.localized_centers;
        tol = zero_tol,
    )
    return _CartesianNestedDoSideTrace1D(
        String(context_label),
        group_kind,
        Int(layer_index),
        piece_kind,
        axis,
        String(usage_label),
        side.interval,
        copy(side.local_centers),
        side.retained_count,
        copy(side.localized_centers),
        copy(side.localized_weights),
        symmetric_about_zero,
        symmetry_error,
        contains_near_zero_center,
        iseven(side.retained_count),
    )
end

function _nested_first_matching_face(
    shell::_CartesianNestedCompleteShell3D,
    face_kind::Symbol,
    fixed_side::Symbol,
)
    index = findfirst(face -> face.face_kind == face_kind && face.fixed_side == fixed_side, shell.faces)
    isnothing(index) && throw(ArgumentError("nested doside trace requires a $(face_kind) $(fixed_side) face"))
    return shell.faces[index]
end

function _nested_first_matching_edge(
    shell::_CartesianNestedCompleteShell3D,
    free_axis::Symbol,
)
    index = findfirst(edge -> edge.free_axis == free_axis, shell.edges)
    isnothing(index) && throw(ArgumentError("nested doside trace requires a $(free_axis)-edge representative"))
    return shell.edges[index]
end

function _nested_complete_shell_doside_traces(
    shell::_CartesianNestedCompleteShell3D,
    context_prefix::AbstractString,
    group_kind::Symbol,
    layer_index::Integer;
    symmetry_tol::Float64 = 1.0e-8,
    zero_tol::Float64 = 1.0e-8,
)
    traces = _CartesianNestedDoSideTrace1D[]

    face_xy = _nested_first_matching_face(shell, :xy, :low)
    push!(
        traces,
        _nested_doside_trace(
            face_xy.side_first;
            context_label = string(context_prefix, "/face_xy/tangential_x"),
            group_kind = group_kind,
            layer_index = layer_index,
            piece_kind = :face_tangent,
            axis = :x,
            usage_label = "face_kind=:xy shared_by=low/high tangential_axis=:x",
            symmetry_tol = symmetry_tol,
            zero_tol = zero_tol,
        ),
    )
    push!(
        traces,
        _nested_doside_trace(
            face_xy.side_second;
            context_label = string(context_prefix, "/face_xy/tangential_y"),
            group_kind = group_kind,
            layer_index = layer_index,
            piece_kind = :face_tangent,
            axis = :y,
            usage_label = "face_kind=:xy shared_by=low/high tangential_axis=:y",
            symmetry_tol = symmetry_tol,
            zero_tol = zero_tol,
        ),
    )

    face_xz = _nested_first_matching_face(shell, :xz, :low)
    push!(
        traces,
        _nested_doside_trace(
            face_xz.side_first;
            context_label = string(context_prefix, "/face_xz/tangential_x"),
            group_kind = group_kind,
            layer_index = layer_index,
            piece_kind = :face_tangent,
            axis = :x,
            usage_label = "face_kind=:xz shared_by=low/high tangential_axis=:x",
            symmetry_tol = symmetry_tol,
            zero_tol = zero_tol,
        ),
    )
    push!(
        traces,
        _nested_doside_trace(
            face_xz.side_second;
            context_label = string(context_prefix, "/face_xz/tangential_z"),
            group_kind = group_kind,
            layer_index = layer_index,
            piece_kind = :face_tangent,
            axis = :z,
            usage_label = "face_kind=:xz shared_by=low/high tangential_axis=:z",
            symmetry_tol = symmetry_tol,
            zero_tol = zero_tol,
        ),
    )

    face_yz = _nested_first_matching_face(shell, :yz, :low)
    push!(
        traces,
        _nested_doside_trace(
            face_yz.side_first;
            context_label = string(context_prefix, "/face_yz/tangential_y"),
            group_kind = group_kind,
            layer_index = layer_index,
            piece_kind = :face_tangent,
            axis = :y,
            usage_label = "face_kind=:yz shared_by=low/high tangential_axis=:y",
            symmetry_tol = symmetry_tol,
            zero_tol = zero_tol,
        ),
    )
    push!(
        traces,
        _nested_doside_trace(
            face_yz.side_second;
            context_label = string(context_prefix, "/face_yz/tangential_z"),
            group_kind = group_kind,
            layer_index = layer_index,
            piece_kind = :face_tangent,
            axis = :z,
            usage_label = "face_kind=:yz shared_by=low/high tangential_axis=:z",
            symmetry_tol = symmetry_tol,
            zero_tol = zero_tol,
        ),
    )

    edge_x = _nested_first_matching_edge(shell, :x)
    push!(
        traces,
        _nested_doside_trace(
            edge_x.side;
            context_label = string(context_prefix, "/edge_x/free_axis_x"),
            group_kind = group_kind,
            layer_index = layer_index,
            piece_kind = :edge_free,
            axis = :x,
            usage_label = "free_axis=:x shared_by=all_boundary_sign_pairs",
            symmetry_tol = symmetry_tol,
            zero_tol = zero_tol,
        ),
    )
    edge_y = _nested_first_matching_edge(shell, :y)
    push!(
        traces,
        _nested_doside_trace(
            edge_y.side;
            context_label = string(context_prefix, "/edge_y/free_axis_y"),
            group_kind = group_kind,
            layer_index = layer_index,
            piece_kind = :edge_free,
            axis = :y,
            usage_label = "free_axis=:y shared_by=all_boundary_sign_pairs",
            symmetry_tol = symmetry_tol,
            zero_tol = zero_tol,
        ),
    )
    edge_z = _nested_first_matching_edge(shell, :z)
    push!(
        traces,
        _nested_doside_trace(
            edge_z.side;
            context_label = string(context_prefix, "/edge_z/free_axis_z"),
            group_kind = group_kind,
            layer_index = layer_index,
            piece_kind = :edge_free,
            axis = :z,
            usage_label = "free_axis=:z shared_by=all_boundary_sign_pairs",
            symmetry_tol = symmetry_tol,
            zero_tol = zero_tol,
        ),
    )

    return traces
end

function _nested_sequence_doside_traces(
    sequence::_CartesianNestedShellSequence3D,
    region_label::AbstractString,
    group_kind::Symbol;
    symmetry_tol::Float64 = 1.0e-8,
    zero_tol::Float64 = 1.0e-8,
)
    traces = _CartesianNestedDoSideTrace1D[]
    for (layer_index, layer) in pairs(sequence.shell_layers)
        layer isa _CartesianNestedCompleteShell3D || continue
        append!(
            traces,
            _nested_complete_shell_doside_traces(
                layer,
                string(region_label, "/layer_", layer_index),
                group_kind,
                layer_index;
                symmetry_tol = symmetry_tol,
                zero_tol = zero_tol,
            ),
        )
    end
    return traces
end

function _bond_aligned_diatomic_doside_traces(
    source::_CartesianNestedBondAlignedDiatomicSource3D;
    symmetry_tol::Float64 = 1.0e-8,
    zero_tol::Float64 = 1.0e-8,
)
    traces = _CartesianNestedDoSideTrace1D[]
    for (layer_index, layer) in pairs(source.shared_shell_layers)
        append!(
            traces,
            _nested_complete_shell_doside_traces(
                layer,
                string("shared_shell/layer_", layer_index),
                :shared_shell,
                layer_index;
                symmetry_tol = symmetry_tol,
                zero_tol = zero_tol,
            ),
        )
    end
    if source.geometry.did_split
        append!(
            traces,
            _nested_sequence_doside_traces(
                source.child_sequences[1],
                "left_child",
                :left_child;
                symmetry_tol = symmetry_tol,
                zero_tol = zero_tol,
            ),
        )
        append!(
            traces,
            _nested_sequence_doside_traces(
                source.child_sequences[2],
                "right_child",
                :right_child;
                symmetry_tol = symmetry_tol,
                zero_tol = zero_tol,
            ),
        )
    else
        append!(
            traces,
            _nested_sequence_doside_traces(
                source.child_sequences[1],
                "shared_child",
                :shared_child;
                symmetry_tol = symmetry_tol,
                zero_tol = zero_tol,
            ),
        )
    end
    return traces
end

function _nested_trace_vector_string(values::AbstractVector{<:Real})
    return "[" * join((string(Float64(value)) for value in values), ", ") * "]"
end

function _bond_aligned_diatomic_doside_trace_notes(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    notes = String[]
    isempty(source.shared_shell_layers) && push!(
        notes,
        "# note shared_shell has no local side contractions",
    )
    if source.geometry.did_split
        isempty(source.child_sequences[1].shell_layers) && push!(
            notes,
            "# note left_child has no local side contractions; it remains a direct core block",
        )
        isempty(source.child_sequences[2].shell_layers) && push!(
            notes,
            "# note right_child has no local side contractions; it remains a direct core block",
        )
    else
        isempty(source.child_sequences[1].shell_layers) && push!(
            notes,
            "# note shared_child has no local side contractions; it remains a direct core block",
        )
    end
    return notes
end

function _write_bond_aligned_diatomic_doside_trace(
    path::AbstractString,
    source::_CartesianNestedBondAlignedDiatomicSource3D;
    symmetry_tol::Float64 = 1.0e-8,
    zero_tol::Float64 = 1.0e-8,
)
    traces = _bond_aligned_diatomic_doside_traces(
        source;
        symmetry_tol = symmetry_tol,
        zero_tol = zero_tol,
    )
    mkpath(dirname(String(path)))
    open(path, "w") do io
        write(io, "# GaussletBases bond-aligned diatomic doside/COMX trace\n")
        write(io, "# bond_axis = $(source.basis.bond_axis)\n")
        write(io, "# working_box = $(source.geometry.working_box)\n")
        write(io, "# did_split = $(source.geometry.did_split)\n")
        if !isnothing(source.geometry.shared_midpoint_box)
            write(io, "# shared_midpoint_box = $(source.geometry.shared_midpoint_box)\n")
        end
        write(io, "# symmetry_tol = $(symmetry_tol)\n")
        write(io, "# zero_tol = $(zero_tol)\n")
        write(io, "# trace_count = $(length(traces))\n")
        for note in _bond_aligned_diatomic_doside_trace_notes(source)
            write(io, note, "\n")
        end
        for (index, trace) in pairs(traces)
            write(io, "\n[trace $(index)]\n")
            write(io, "context_label = $(trace.context_label)\n")
            write(io, "group_kind = $(trace.group_kind)\n")
            write(io, "layer_index = $(trace.layer_index)\n")
            write(io, "piece_kind = $(trace.piece_kind)\n")
            write(io, "axis = $(trace.axis)\n")
            write(io, "usage = $(trace.usage_label)\n")
            write(io, "interval = $(first(trace.interval)):$(last(trace.interval))\n")
            write(io, "parent_centers = $(_nested_trace_vector_string(trace.parent_centers))\n")
            write(io, "retained_count = $(trace.retained_count)\n")
            write(io, "localized_centers = $(_nested_trace_vector_string(trace.localized_centers))\n")
            write(io, "localized_weights = $(_nested_trace_vector_string(trace.localized_weights))\n")
            write(io, "symmetric_about_zero = $(trace.symmetric_about_zero)\n")
            write(io, "symmetry_error = $(trace.symmetry_error)\n")
            write(io, "contains_near_zero_center = $(trace.contains_near_zero_center)\n")
            write(io, "even_retained_count = $(trace.even_retained_count)\n")
        end
    end
    return traces
end

function _embed_local_side_coefficients(
    local_coefficients::AbstractMatrix{<:Real},
    interval::UnitRange{Int},
    n1d::Int,
)
    full_coefficients = zeros(Float64, n1d, size(local_coefficients, 2))
    full_coefficients[interval, :] .= Matrix{Float64}(local_coefficients)
    return full_coefficients
end

# Alg Nested-Face step 3: Build a local 1D doside contraction on one interval,
# forcing odd retained counts on intervals symmetric about zero before COMX.
# See docs/src/algorithms/cartesian_nested_face_construction.md.
function _nested_doside_1d(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    interval::UnitRange{Int},
    retained_count::Int,
)
    interval_data = _nested_interval_data(pgdg, interval)
    retained_count = _nested_doside_retained_count(interval_data.centers, retained_count)
    raw_basis = _nested_retained_span(
        interval_data.weights,
        interval_data.centers,
        interval_data.position,
        interval_data.overlap,
        retained_count,
    )
    overlap_seed = Matrix{Float64}(transpose(raw_basis) * interval_data.overlap * raw_basis)
    position_seed = Matrix{Float64}(transpose(raw_basis) * interval_data.position * raw_basis)
    sign_vector = vec(transpose(raw_basis) * interval_data.weights)
    transform, localized_centers = _cleanup_comx_transform(overlap_seed, position_seed, sign_vector)
    local_coefficients = Matrix{Float64}(raw_basis * transform)
    localized_weights = vec(transpose(interval_data.weights) * local_coefficients)
    coefficient_matrix = _embed_local_side_coefficients(local_coefficients, interval, interval_data.n1d)
    return _CartesianNestedDoSide1D(
        interval,
        retained_count,
        interval_data.overlap,
        interval_data.position,
        interval_data.weights,
        interval_data.centers,
        local_coefficients,
        coefficient_matrix,
        localized_centers,
        localized_weights,
    )
end

function _nested_doside_1d(
    bundle::_MappedOrdinaryGausslet1DBundle,
    interval::UnitRange{Int},
    retained_count::Int,
)
    return _nested_doside_1d(bundle.pgdg_intermediate, interval, retained_count)
end

function _cartesian_flat_index(
    ix::Int,
    iy::Int,
    iz::Int,
    dims::NTuple{3,Int},
)
    nx, ny, nz = dims
    1 <= ix <= nx || throw(ArgumentError("x index must lie inside the parent Cartesian box"))
    1 <= iy <= ny || throw(ArgumentError("y index must lie inside the parent Cartesian box"))
    1 <= iz <= nz || throw(ArgumentError("z index must lie inside the parent Cartesian box"))
    return (ix - 1) * ny * nz + (iy - 1) * nz + iz
end

function _cartesian_flat_index(ix::Int, iy::Int, iz::Int, n1d::Int)
    return _cartesian_flat_index(ix, iy, iz, (n1d, n1d, n1d))
end

function _cartesian_unflat_index(index::Int, dims::NTuple{3,Int})
    nx, ny, nz = dims
    1 <= index <= nx * ny * nz || throw(
        ArgumentError("flat parent index must lie inside the parent Cartesian box"),
    )
    shifted = index - 1
    plane = ny * nz
    ix = shifted ÷ plane + 1
    remainder = shifted % plane
    iy = remainder ÷ nz + 1
    iz = remainder % nz + 1
    return (ix, iy, iz)
end

function _cartesian_unflat_index(index::Int, n1d::Int)
    return _cartesian_unflat_index(index, (n1d, n1d, n1d))
end

function _nested_xy_face_support_indices(
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_index::Int,
    dims::NTuple{3,Int},
)
    support = Int[]
    for ix in x_interval, iy in y_interval
        push!(support, _cartesian_flat_index(ix, iy, z_index, dims))
    end
    return support
end

function _nested_xy_face_support_indices(
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_index::Int,
    n1d::Int,
)
    return _nested_xy_face_support_indices(x_interval, y_interval, z_index, (n1d, n1d, n1d))
end

# Alg Nested-Face steps 5-7: Build one simple x-y face product from two local
# side spaces, keeping only the supplied face interior intervals so opposite
# faces remain disjoint.
# See docs/src/algorithms/cartesian_nested_face_construction.md.
function _nested_xy_face_product(
    side_x::_CartesianNestedDoSide1D,
    side_y::_CartesianNestedDoSide1D,
    z_index::Int,
    n1d::Int,
)
    1 <= z_index <= n1d || throw(ArgumentError("nested x-y face requires a fixed z index inside the finalized Cartesian line"))
    nx = size(side_x.coefficient_matrix, 2)
    ny = size(side_y.coefficient_matrix, 2)
    coefficients = zeros(Float64, n1d^3, nx * ny)
    column = 0
    for ix_side in 1:nx, iy_side in 1:ny
        column += 1
        for ix in side_x.interval
            xvalue = side_x.coefficient_matrix[ix, ix_side]
            iszero(xvalue) && continue
            for iy in side_y.interval
                yvalue = side_y.coefficient_matrix[iy, iy_side]
                iszero(yvalue) && continue
                coefficients[_cartesian_flat_index(ix, iy, z_index, n1d), column] = xvalue * yvalue
            end
        end
    end
    support_indices = _nested_xy_face_support_indices(side_x.interval, side_y.interval, z_index, n1d)
    return _CartesianNestedXYFace3D(
        z_index,
        side_x,
        side_y,
        coefficients,
        support_indices,
    )
end

function _nested_xy_face_product(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_index::Int;
    retain_x::Int,
    retain_y::Int,
)
    side_x = _nested_doside_1d(pgdg, x_interval, retain_x)
    side_y = _nested_doside_1d(pgdg, y_interval, retain_y)
    n1d = size(pgdg.overlap, 1)
    return _nested_xy_face_product(side_x, side_y, z_index, n1d)
end

function _nested_xy_face_overlap(
    face::_CartesianNestedXYFace3D,
    overlap_1d::AbstractMatrix{<:Real},
)
    x_overlap = Matrix{Float64}(transpose(face.side_x.coefficient_matrix) * overlap_1d * face.side_x.coefficient_matrix)
    y_overlap = Matrix{Float64}(transpose(face.side_y.coefficient_matrix) * overlap_1d * face.side_y.coefficient_matrix)
    z_norm = Float64(overlap_1d[face.z_index, face.z_index])
    return z_norm .* kron(x_overlap, y_overlap)
end

function _nested_xy_face_cross_overlap(
    face_left::_CartesianNestedXYFace3D,
    face_right::_CartesianNestedXYFace3D,
    overlap_1d::AbstractMatrix{<:Real},
)
    x_overlap = Matrix{Float64}(transpose(face_left.side_x.coefficient_matrix) * overlap_1d * face_right.side_x.coefficient_matrix)
    y_overlap = Matrix{Float64}(transpose(face_left.side_y.coefficient_matrix) * overlap_1d * face_right.side_y.coefficient_matrix)
    z_overlap = Float64(overlap_1d[face_left.z_index, face_right.z_index])
    return z_overlap .* kron(x_overlap, y_overlap)
end

function _nested_face_axes(face_kind::Symbol)
    if face_kind == :xy
        return ((:x, :y), :z)
    elseif face_kind == :xz
        return ((:x, :z), :y)
    elseif face_kind == :yz
        return ((:y, :z), :x)
    else
        throw(ArgumentError("nested face kind must be :xy, :xz, or :yz"))
    end
end

function _nested_face_support_indices(
    face_kind::Symbol,
    interval_first::UnitRange{Int},
    interval_second::UnitRange{Int},
    fixed_index::Int,
    dims::NTuple{3,Int},
)
    support = Int[]
    if face_kind == :xy
        for ix in interval_first, iy in interval_second
            push!(support, _cartesian_flat_index(ix, iy, fixed_index, dims))
        end
    elseif face_kind == :xz
        for ix in interval_first, iz in interval_second
            push!(support, _cartesian_flat_index(ix, fixed_index, iz, dims))
        end
    elseif face_kind == :yz
        for iy in interval_first, iz in interval_second
            push!(support, _cartesian_flat_index(fixed_index, iy, iz, dims))
        end
    else
        throw(ArgumentError("nested face support requires face kind :xy, :xz, or :yz"))
    end
    return support
end

function _nested_face_support_indices(
    face_kind::Symbol,
    interval_first::UnitRange{Int},
    interval_second::UnitRange{Int},
    fixed_index::Int,
    n1d::Int,
)
    return _nested_face_support_indices(face_kind, interval_first, interval_second, fixed_index, (n1d, n1d, n1d))
end

function _nested_edge_support_indices(
    free_axis::Symbol,
    free_interval::UnitRange{Int},
    fixed_indices::NTuple{2,Int},
    dims::NTuple{3,Int},
)
    support = Int[]
    if free_axis == :x
        iy, iz = fixed_indices
        for ix in free_interval
            push!(support, _cartesian_flat_index(ix, iy, iz, dims))
        end
    elseif free_axis == :y
        ix, iz = fixed_indices
        for iy in free_interval
            push!(support, _cartesian_flat_index(ix, iy, iz, dims))
        end
    elseif free_axis == :z
        ix, iy = fixed_indices
        for iz in free_interval
            push!(support, _cartesian_flat_index(ix, iy, iz, dims))
        end
    else
        throw(ArgumentError("nested edge support requires free axis :x, :y, or :z"))
    end
    return support
end

function _nested_edge_support_indices(
    free_axis::Symbol,
    free_interval::UnitRange{Int},
    fixed_indices::NTuple{2,Int},
    n1d::Int,
)
    return _nested_edge_support_indices(free_axis, free_interval, fixed_indices, (n1d, n1d, n1d))
end

function _nested_edge_product(
    free_axis::Symbol,
    fixed_sides::NTuple{2,Symbol},
    side::_CartesianNestedDoSide1D,
    fixed_indices::NTuple{2,Int},
    dims::NTuple{3,Int},
)
    limits =
        free_axis == :x ? (dims[2], dims[3]) :
        free_axis == :y ? (dims[1], dims[3]) :
        (dims[1], dims[2])
    1 <= fixed_indices[1] <= limits[1] || throw(
        ArgumentError("nested edge requires the first fixed index inside the finalized Cartesian line"),
    )
    1 <= fixed_indices[2] <= limits[2] || throw(
        ArgumentError("nested edge requires the second fixed index inside the finalized Cartesian line"),
    )
    coefficient_matrix = zeros(Float64, prod(dims), size(side.coefficient_matrix, 2))
    for col in 1:size(side.coefficient_matrix, 2)
        for free_index in side.interval
            value = side.coefficient_matrix[free_index, col]
            iszero(value) && continue
            flat =
                free_axis == :x ? _cartesian_flat_index(free_index, fixed_indices[1], fixed_indices[2], dims) :
                free_axis == :y ? _cartesian_flat_index(fixed_indices[1], free_index, fixed_indices[2], dims) :
                _cartesian_flat_index(fixed_indices[1], fixed_indices[2], free_index, dims)
            coefficient_matrix[flat, col] = value
        end
    end
    fixed_axes =
        free_axis == :x ? (:y, :z) :
        free_axis == :y ? (:x, :z) :
        (:x, :y)
    support_indices = _nested_edge_support_indices(free_axis, side.interval, fixed_indices, dims)
    return _CartesianNestedEdge3D(
        free_axis,
        fixed_axes,
        fixed_sides,
        fixed_indices,
        side,
        coefficient_matrix,
        support_indices,
    )
end

function _nested_edge_product(
    free_axis::Symbol,
    fixed_sides::NTuple{2,Symbol},
    side::_CartesianNestedDoSide1D,
    fixed_indices::NTuple{2,Int},
    n1d::Int,
)
    return _nested_edge_product(free_axis, fixed_sides, side, fixed_indices, (n1d, n1d, n1d))
end

function _nested_corner_piece(
    fixed_sides::NTuple{3,Symbol},
    fixed_indices::NTuple{3,Int},
    dims::NTuple{3,Int},
)
    flat = _cartesian_flat_index(fixed_indices[1], fixed_indices[2], fixed_indices[3], dims)
    coefficient_matrix = zeros(Float64, prod(dims), 1)
    coefficient_matrix[flat, 1] = 1.0
    return _CartesianNestedCorner3D(
        fixed_sides,
        fixed_indices,
        coefficient_matrix,
        [flat],
    )
end

function _nested_corner_piece(
    fixed_sides::NTuple{3,Symbol},
    fixed_indices::NTuple{3,Int},
    n1d::Int,
)
    return _nested_corner_piece(fixed_sides, fixed_indices, (n1d, n1d, n1d))
end

function _nested_box_support_indices(
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int},
    dims::NTuple{3,Int},
)
    support = Int[]
    for ix in x_interval, iy in y_interval, iz in z_interval
        push!(support, _cartesian_flat_index(ix, iy, iz, dims))
    end
    return support
end

function _nested_box_support_indices(
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int},
    n1d::Int,
)
    return _nested_box_support_indices(x_interval, y_interval, z_interval, (n1d, n1d, n1d))
end

function _nested_direct_core_coefficients(
    support_indices::AbstractVector{Int},
    nparent::Int,
)
    coefficients = zeros(Float64, nparent, length(support_indices))
    for (column, index) in enumerate(support_indices)
        coefficients[index, column] = 1.0
    end
    return coefficients
end

function _nested_product_coefficients(
    x_side::_CartesianNestedDoSide1D,
    y_side::_CartesianNestedDoSide1D,
    z_side::_CartesianNestedDoSide1D,
    dims::NTuple{3,Int},
)
    nparent = prod(dims)
    ncols = size(x_side.coefficient_matrix, 2) * size(y_side.coefficient_matrix, 2) * size(z_side.coefficient_matrix, 2)
    coefficients = zeros(Float64, nparent, ncols)
    column = 0
    for ixcol in 1:size(x_side.coefficient_matrix, 2),
        iycol in 1:size(y_side.coefficient_matrix, 2),
        izcol in 1:size(z_side.coefficient_matrix, 2)
        column += 1
        for ix in x_side.interval
            vx = x_side.coefficient_matrix[ix, ixcol]
            iszero(vx) && continue
            for iy in y_side.interval
                vy = y_side.coefficient_matrix[iy, iycol]
                iszero(vy) && continue
                for iz in z_side.interval
                    vz = z_side.coefficient_matrix[iz, izcol]
                    iszero(vz) && continue
                    flat = _cartesian_flat_index(ix, iy, iz, dims)
                    coefficients[flat, column] = vx * vy * vz
                end
            end
        end
    end
    return coefficients
end

function _nested_product_coefficients(
    x_side::_CartesianNestedDoSide1D,
    y_side::_CartesianNestedDoSide1D,
    z_side::_CartesianNestedDoSide1D,
    n1d::Int,
)
    return _nested_product_coefficients(x_side, y_side, z_side, (n1d, n1d, n1d))
end

function _nested_contracted_core_coefficients(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int};
    retain_x::Int,
    retain_y::Int,
    retain_z::Int,
)
    n1d = size(pgdg.overlap, 1)
    x_side = _nested_doside_1d(pgdg, x_interval, retain_x)
    y_side = _nested_doside_1d(pgdg, y_interval, retain_y)
    z_side = _nested_doside_1d(pgdg, z_interval, retain_z)
    support_indices = _nested_box_support_indices(x_interval, y_interval, z_interval, n1d)
    coefficients = _nested_product_coefficients(x_side, y_side, z_side, n1d)
    return (
        support_indices = support_indices,
        coefficient_matrix = coefficients,
    )
end

function _nested_sequence_support_indices(
    core_indices::AbstractVector{Int},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D},
)
    support = Int[]
    seen = Set{Int}()
    for index in core_indices
        if index in seen
            throw(ArgumentError("nested shell-sequence construction requires distinct core support rows"))
        end
        push!(support, index)
        push!(seen, index)
    end
    for shell in shell_layers
        for index in shell.support_indices
            if index in seen
                throw(ArgumentError("nested shell-sequence construction requires disjoint core and shell-layer supports"))
            end
            push!(support, index)
            push!(seen, index)
        end
    end
    sort!(support)
    return support
end

function _nested_sequence_working_box(
    core_indices::AbstractVector{Int},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D},
    dims::NTuple{3,Int},
)
    xmin = typemax(Int)
    ymin = typemax(Int)
    zmin = typemax(Int)
    xmax = typemin(Int)
    ymax = typemin(Int)
    zmax = typemin(Int)

    function update_bounds(state::NTuple{3,Int})
        ix, iy, iz = state
        xmin = min(xmin, ix)
        ymin = min(ymin, iy)
        zmin = min(zmin, iz)
        xmax = max(xmax, ix)
        ymax = max(ymax, iy)
        zmax = max(zmax, iz)
        return nothing
    end

    for index in core_indices
        update_bounds(_cartesian_unflat_index(index, dims))
    end
    for shell in shell_layers, state in shell.support_states
        update_bounds(state)
    end

    xmin <= xmax || throw(ArgumentError("nested shell-sequence construction requires at least one retained parent row"))
    return (xmin:xmax, ymin:ymax, zmin:zmax)
end

function _nested_sequence_working_box(
    core_indices::AbstractVector{Int},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D},
    n1d::Int,
)
    return _nested_sequence_working_box(core_indices, shell_layers, (n1d, n1d, n1d))
end

function _nested_assert_sequence_coverage(
    core_indices::AbstractVector{Int},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D},
    support_indices::AbstractVector{Int},
    dims::NTuple{3,Int},
)
    working_box = _nested_sequence_working_box(core_indices, shell_layers, dims)
    target_indices = _nested_box_support_indices(working_box..., dims)
    if target_indices != support_indices
        support_set = Set(support_indices)
        target_set = Set(target_indices)
        missing = length(setdiff(target_set, support_set))
        extra = length(setdiff(support_set, target_set))
        throw(ArgumentError("nested shell-sequence construction requires full coverage of the inferred working box $(working_box): missing $missing parent rows and extra $extra rows"))
    end
    return working_box
end

function _nested_assert_sequence_coverage(
    core_indices::AbstractVector{Int},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D},
    support_indices::AbstractVector{Int},
    n1d::Int,
)
    return _nested_assert_sequence_coverage(core_indices, shell_layers, support_indices, (n1d, n1d, n1d))
end

function _nested_shell_sequence_piece_ownership_audit(
    sequence::_CartesianNestedShellSequence3D,
)
    range_groups = UnitRange{Int}[sequence.core_column_range]
    append!(range_groups, sequence.layer_column_ranges)
    group_counts = Int[]
    coefficient_matrix = sequence.coefficient_matrix
    for row in axes(coefficient_matrix, 1)
        nzcols = findall(!iszero, @view coefficient_matrix[row, :])
        touched_groups = 0
        for range in range_groups
            any(col -> col in range, nzcols) && (touched_groups += 1)
        end
        push!(group_counts, touched_groups)
    end
    return (
        min_group_count = minimum(group_counts),
        max_group_count = maximum(group_counts),
        unowned_row_count = count(iszero, group_counts),
        multi_owned_row_count = count(>(1), group_counts),
    )
end

function _nested_shell_sequence_contract_audit(
    sequence::_CartesianNestedShellSequence3D,
    parent_dims::NTuple{3,Int},
)
    expected_box = (
        1:parent_dims[1],
        1:parent_dims[2],
        1:parent_dims[3],
    )
    ownership = _nested_shell_sequence_piece_ownership_audit(sequence)
    expected_support_count = prod(parent_dims)
    support_count = length(sequence.support_indices)
    return CartesianNestedSequenceContractAudit(
        parent_dims,
        sequence.working_box,
        sequence.working_box == expected_box,
        support_count,
        expected_support_count,
        expected_support_count - support_count,
        ownership.min_group_count,
        ownership.max_group_count,
        ownership.unowned_row_count,
        ownership.multi_owned_row_count,
    )
end

function _nested_complete_shell_retention_from_nside(nside::Int)
    nside >= 3 || throw(ArgumentError("nested complete-shell retention requires nside >= 3"))
    retained_side = nside - 2
    return CartesianNestedCompleteShellRetentionContract(
        nside,
        (retained_side, retained_side),
        (retained_side, retained_side),
        (retained_side, retained_side),
        retained_side,
        retained_side,
        retained_side,
        6 * retained_side^2,
        12 * retained_side,
        8,
        nside^3 - (nside - 2)^3,
        true,
    )
end

function _nested_resolve_complete_shell_retention(
    nside::Int;
    retain_xy::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_xz::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_yz::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_x_edge::Union{Nothing,Int} = nothing,
    retain_y_edge::Union{Nothing,Int} = nothing,
    retain_z_edge::Union{Nothing,Int} = nothing,
)
    default_contract = _nested_complete_shell_retention_from_nside(nside)
    actual_retain_xy = something(retain_xy, default_contract.retain_xy)
    actual_retain_xz = something(retain_xz, default_contract.retain_xz)
    actual_retain_yz = something(retain_yz, default_contract.retain_yz)
    actual_retain_x_edge = something(retain_x_edge, default_contract.retain_x_edge)
    actual_retain_y_edge = something(retain_y_edge, default_contract.retain_y_edge)
    actual_retain_z_edge = something(retain_z_edge, default_contract.retain_z_edge)
    face_retained_count =
        2 * (actual_retain_xy[1] * actual_retain_xy[2] +
             actual_retain_xz[1] * actual_retain_xz[2] +
             actual_retain_yz[1] * actual_retain_yz[2])
    edge_retained_count =
        4 * (actual_retain_x_edge + actual_retain_y_edge + actual_retain_z_edge)
    return CartesianNestedCompleteShellRetentionContract(
        nside,
        actual_retain_xy,
        actual_retain_xz,
        actual_retain_yz,
        actual_retain_x_edge,
        actual_retain_y_edge,
        actual_retain_z_edge,
        face_retained_count,
        edge_retained_count,
        8,
        face_retained_count + edge_retained_count + 8,
        actual_retain_xy == default_contract.retain_xy &&
        actual_retain_xz == default_contract.retain_xz &&
        actual_retain_yz == default_contract.retain_yz &&
        actual_retain_x_edge == default_contract.retain_x_edge &&
        actual_retain_y_edge == default_contract.retain_y_edge &&
        actual_retain_z_edge == default_contract.retain_z_edge,
    )
end

function _nested_shrunk_interval(
    interval::UnitRange{Int},
    nlayers::Integer;
    nside::Int,
)
    nside >= 1 || throw(ArgumentError("nested fixed-nside policy requires nside >= 1"))
    length(interval) >= nside || throw(
        ArgumentError("nested fixed-nside policy requires the starting direct core interval to have length at least nside"),
    )
    max_shrinks = max(0, (length(interval) - nside) ÷ 2)
    nshrinks = min(Int(nlayers), max_shrinks)
    return (first(interval) + nshrinks):(last(interval) - nshrinks)
end

function _nested_shrunk_box(
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int},
    nlayers::Integer;
    nside::Int,
)
    return (
        _nested_shrunk_interval(x_interval, nlayers; nside = nside),
        _nested_shrunk_interval(y_interval, nlayers; nside = nside),
        _nested_shrunk_interval(z_interval, nlayers; nside = nside),
    )
end

function _nested_face_product(
    face_kind::Symbol,
    fixed_side::Symbol,
    side_first::_CartesianNestedDoSide1D,
    side_second::_CartesianNestedDoSide1D,
    fixed_index::Int,
    dims::NTuple{3,Int},
)
    max_fixed =
        face_kind == :xy ? dims[3] :
        face_kind == :xz ? dims[2] :
        dims[1]
    1 <= fixed_index <= max_fixed || throw(ArgumentError("nested face requires a fixed index inside the finalized Cartesian line"))
    (fixed_side == :low || fixed_side == :high) || throw(ArgumentError("nested face fixed_side must be :low or :high"))
    _, fixed_axis = _nested_face_axes(face_kind)
    nfirst = size(side_first.coefficient_matrix, 2)
    nsecond = size(side_second.coefficient_matrix, 2)
    coefficients = zeros(Float64, prod(dims), nfirst * nsecond)
    column = 0
    for ifirst in 1:nfirst, isecond in 1:nsecond
        column += 1
        for index_first in side_first.interval
            value_first = side_first.coefficient_matrix[index_first, ifirst]
            iszero(value_first) && continue
            for index_second in side_second.interval
                value_second = side_second.coefficient_matrix[index_second, isecond]
                iszero(value_second) && continue
                flat_index =
                    face_kind == :xy ? _cartesian_flat_index(index_first, index_second, fixed_index, dims) :
                    face_kind == :xz ? _cartesian_flat_index(index_first, fixed_index, index_second, dims) :
                    _cartesian_flat_index(fixed_index, index_first, index_second, dims)
                coefficients[flat_index, column] = value_first * value_second
            end
        end
    end
    support_indices = _nested_face_support_indices(
        face_kind,
        side_first.interval,
        side_second.interval,
        fixed_index,
        dims,
    )
    return _CartesianNestedFace3D(
        face_kind,
        fixed_axis,
        fixed_side,
        fixed_index,
        side_first,
        side_second,
        coefficients,
        support_indices,
    )
end

function _nested_face_product(
    face_kind::Symbol,
    fixed_side::Symbol,
    side_first::_CartesianNestedDoSide1D,
    side_second::_CartesianNestedDoSide1D,
    fixed_index::Int,
    n1d::Int,
)
    return _nested_face_product(face_kind, fixed_side, side_first, side_second, fixed_index, (n1d, n1d, n1d))
end

function _nested_support_product_matrix(
    support_states::AbstractVector{<:NTuple{3,Int}},
    operator_x::AbstractMatrix{<:Real},
    operator_y::AbstractMatrix{<:Real},
    operator_z::AbstractMatrix{<:Real},
)
    # Retained as a convenience/reference wrapper. The hot packet-assembly path
    # now reuses one workspace through `_nested_fill_support_product_matrix!`
    # and `_nested_contract_support_product!` instead of allocating several
    # support-scale matrices at once.
    matrix = Matrix{Float64}(undef, length(support_states), length(support_states))
    return _nested_fill_support_product_matrix!(
        matrix,
        support_states,
        operator_x,
        operator_y,
        operator_z,
    )
end

function _nested_support_axes(
    support_states::AbstractVector{<:NTuple{3,Int}},
)
    x = Vector{Int}(undef, length(support_states))
    y = Vector{Int}(undef, length(support_states))
    z = Vector{Int}(undef, length(support_states))
    @inbounds for index in eachindex(support_states)
        ix, iy, iz = support_states[index]
        x[index] = ix
        y[index] = iy
        z[index] = iz
    end
    return _CartesianNestedSupportAxes3D(x, y, z)
end

function _nested_normalize_packet_kernel(packet_kernel::Symbol)
    packet_kernel in (:support_reference, :factorized_direct) || throw(
        ArgumentError("nested packet kernel must be :support_reference or :factorized_direct"),
    )
    return packet_kernel
end

function _nested_zero_small!(values::AbstractVector{Float64}; atol::Float64)
    @inbounds for index in eachindex(values)
        abs(values[index]) <= atol && (values[index] = 0.0)
    end
    return values
end

function _nested_find_or_push_axis_function!(
    functions::Vector{Vector{Float64}},
    candidate::Vector{Float64};
    atol::Float64,
)
    for (index, existing) in enumerate(functions)
        length(existing) == length(candidate) || continue
        maximum(abs.(existing .- candidate)) <= atol && return index
    end
    push!(functions, candidate)
    return length(functions)
end

function _nested_extract_factorized_basis(
    coefficient_matrix::AbstractMatrix{<:Real},
    dims::NTuple{3,Int};
    atol::Float64 = 1.0e-12,
)
    nparent = prod(dims)
    size(coefficient_matrix, 1) == nparent || throw(
        ArgumentError("nested factorized-basis extraction requires parent rows matching the Cartesian box volume"),
    )
    nfixed = size(coefficient_matrix, 2)
    x_functions = Vector{Vector{Float64}}()
    y_functions = Vector{Vector{Float64}}()
    z_functions = Vector{Vector{Float64}}()
    basis_triplets = Vector{NTuple{3,Int}}(undef, nfixed)
    basis_amplitudes = Vector{Float64}(undef, nfixed)
    reconstruction_max_error = 0.0

    for column in 1:nfixed
        coefficients = @view coefficient_matrix[:, column]
        anchor = findfirst(value -> abs(value) > atol, coefficients)
        isnothing(anchor) && throw(
            ArgumentError("nested factorized-basis extraction requires every retained fixed column to have at least one nonzero parent row"),
        )
        ix0, iy0, iz0 = _cartesian_unflat_index(anchor, dims)
        amplitude = Float64(coefficients[anchor])
        x_vector = Vector{Float64}(undef, dims[1])
        y_vector = Vector{Float64}(undef, dims[2])
        z_vector = Vector{Float64}(undef, dims[3])
        @inbounds for ix in 1:dims[1]
            x_vector[ix] = Float64(coefficients[_cartesian_flat_index(ix, iy0, iz0, dims)]) / amplitude
        end
        @inbounds for iy in 1:dims[2]
            y_vector[iy] = Float64(coefficients[_cartesian_flat_index(ix0, iy, iz0, dims)]) / amplitude
        end
        @inbounds for iz in 1:dims[3]
            z_vector[iz] = Float64(coefficients[_cartesian_flat_index(ix0, iy0, iz, dims)]) / amplitude
        end
        _nested_zero_small!(x_vector; atol = atol)
        _nested_zero_small!(y_vector; atol = atol)
        _nested_zero_small!(z_vector; atol = atol)

        column_error = 0.0
        @inbounds for iz in 1:dims[3], iy in 1:dims[2], ix in 1:dims[1]
            flat = _cartesian_flat_index(ix, iy, iz, dims)
            expected = amplitude * x_vector[ix] * y_vector[iy] * z_vector[iz]
            column_error = max(column_error, abs(expected - Float64(coefficients[flat])))
        end
        reconstruction_max_error = max(reconstruction_max_error, column_error)
        column_error <= 1.0e3 * atol || throw(
            ArgumentError("nested factorized-basis extraction failed to reconstruct fixed column $column to roundoff (max error = $column_error)"),
        )

        x_index = _nested_find_or_push_axis_function!(x_functions, x_vector; atol = atol)
        y_index = _nested_find_or_push_axis_function!(y_functions, y_vector; atol = atol)
        z_index = _nested_find_or_push_axis_function!(z_functions, z_vector; atol = atol)
        basis_triplets[column] = (x_index, y_index, z_index)
        basis_amplitudes[column] = amplitude
    end

    return _CartesianNestedFactorizedBasis3D(
        dims,
        hcat(x_functions...),
        hcat(y_functions...),
        hcat(z_functions...),
        basis_triplets,
        basis_amplitudes,
        reconstruction_max_error,
    )
end

function _nested_reconstruct_factorized_coefficients(
    factorized_basis::_CartesianNestedFactorizedBasis3D,
)
    dims = factorized_basis.dims
    nparent = prod(dims)
    nfixed = length(factorized_basis.basis_triplets)
    coefficients = zeros(Float64, nparent, nfixed)
    for column in 1:nfixed
        ix_function, iy_function, iz_function = factorized_basis.basis_triplets[column]
        amplitude = factorized_basis.basis_amplitudes[column]
        x_vector = @view factorized_basis.x_functions[:, ix_function]
        y_vector = @view factorized_basis.y_functions[:, iy_function]
        z_vector = @view factorized_basis.z_functions[:, iz_function]
        @inbounds for iz in 1:dims[3], iy in 1:dims[2], ix in 1:dims[1]
            coefficients[_cartesian_flat_index(ix, iy, iz, dims), column] =
                amplitude * x_vector[ix] * y_vector[iy] * z_vector[iz]
        end
    end
    return coefficients
end

function _nested_factorized_axis_weight_projections(
    axis_functions::AbstractMatrix{<:Real},
    weights::AbstractVector{<:Real},
)
    size(axis_functions, 1) == length(weights) || throw(
        ArgumentError("nested factorized axis-weight projection requires one weight per parent-axis site"),
    )
    return vec(transpose(axis_functions) * weights)
end

function _nested_factorized_axis_term_tables(
    operator_terms::Array{Float64,3},
    axis_functions::AbstractMatrix{<:Real},
)
    nterms = size(operator_terms, 1)
    nfunctions = size(axis_functions, 2)
    left_scratch = Matrix{Float64}(undef, nfunctions, size(axis_functions, 1))
    term_tables = Array{Float64,3}(undef, nterms, nfunctions, nfunctions)
    for term in 1:nterms
        mul!(left_scratch, transpose(axis_functions), @view(operator_terms[term, :, :]))
        mul!(@view(term_tables[term, :, :]), left_scratch, axis_functions)
    end
    return term_tables
end

function _nested_factorized_axis_matrix_table(
    operator::AbstractMatrix{<:Real},
    axis_functions::AbstractMatrix{<:Real},
    left_scratch::AbstractMatrix{<:Real},
)
    nfunctions = size(axis_functions, 2)
    size(left_scratch) == (nfunctions, size(axis_functions, 1)) || throw(
        ArgumentError("nested factorized axis matrix tables require scratch sized to the intermediate-function count and parent-axis length"),
    )
    table = Matrix{Float64}(undef, nfunctions, nfunctions)
    mul!(left_scratch, transpose(axis_functions), operator)
    mul!(table, left_scratch, axis_functions)
    return table
end

function _nested_factorized_axis_base_tables(
    axis_functions::AbstractMatrix{<:Real},
    overlap::AbstractMatrix{<:Real},
    kinetic::AbstractMatrix{<:Real},
    position::AbstractMatrix{<:Real},
    x2::AbstractMatrix{<:Real},
)
    nfunctions = size(axis_functions, 2)
    left_scratch = Matrix{Float64}(undef, nfunctions, size(axis_functions, 1))
    return _CartesianNestedFactorizedAxisBaseTables(
        _nested_factorized_axis_matrix_table(overlap, axis_functions, left_scratch),
        _nested_factorized_axis_matrix_table(kinetic, axis_functions, left_scratch),
        _nested_factorized_axis_matrix_table(position, axis_functions, left_scratch),
        _nested_factorized_axis_matrix_table(x2, axis_functions, left_scratch),
    )
end

function _nested_factorized_normalized_pair_term_tables(
    raw_term_tables::Array{Float64,3},
    axis_weight_projections::AbstractVector{<:Real},
)
    nterms, nfunctions_left, nfunctions_right = size(raw_term_tables)
    nfunctions_left == length(axis_weight_projections) == nfunctions_right || throw(
        ArgumentError("nested factorized pair normalization requires one axis weight per unique intermediate function"),
    )
    normalized = similar(raw_term_tables)
    @inbounds for j in 1:nfunctions_right, i in 1:nfunctions_left
        scale = Float64(axis_weight_projections[i]) * Float64(axis_weight_projections[j])
        abs(scale) > 1.0e-14 || throw(
            ArgumentError("nested factorized pair normalization requires nonzero axis-weight projection pairs"),
        )
        for term in 1:nterms
            normalized[term, i, j] = raw_term_tables[term, i, j] / scale
        end
    end
    return normalized
end

function _nested_fill_factorized_term_family!(
    destination_terms::Array{Float64,3},
    factorized_basis::_CartesianNestedFactorizedBasis3D,
    operator_terms_x::Array{Float64,3},
    operator_terms_y::Array{Float64,3},
    operator_terms_z::Array{Float64,3};
    include_basis_amplitudes::Bool,
)
    nterms = size(destination_terms, 1)
    nbasis = length(factorized_basis.basis_triplets)
    size(destination_terms, 2) == nbasis == size(destination_terms, 3) || throw(
        ArgumentError("nested factorized term-family fill requires square output sized to the retained fixed basis"),
    )
    nterms == size(operator_terms_x, 1) == size(operator_terms_y, 1) == size(operator_terms_z, 1) || throw(
        ArgumentError("nested factorized term-family fill requires matching term counts"),
    )
    amplitudes = factorized_basis.basis_amplitudes
    triplets = factorized_basis.basis_triplets
    @inbounds for column in 1:nbasis
        xj, yj, zj = triplets[column]
        amplitude_j = include_basis_amplitudes ? amplitudes[column] : 1.0
        for row in 1:column
            xi, yi, zi = triplets[row]
            scale = include_basis_amplitudes ? amplitudes[row] * amplitude_j : 1.0
            for term in 1:nterms
                value =
                    scale *
                    operator_terms_x[term, xi, xj] *
                    operator_terms_y[term, yi, yj] *
                    operator_terms_z[term, zi, zj]
                destination_terms[term, row, column] = value
                destination_terms[term, column, row] = value
            end
        end
    end
    return destination_terms
end

function _nested_fill_factorized_product_matrix!(
    destination::AbstractMatrix{<:Real},
    factorized_basis::_CartesianNestedFactorizedBasis3D,
    operator_x::AbstractMatrix{<:Real},
    operator_y::AbstractMatrix{<:Real},
    operator_z::AbstractMatrix{<:Real};
    include_basis_amplitudes::Bool = true,
)
    nbasis = length(factorized_basis.basis_triplets)
    size(destination) == (nbasis, nbasis) || throw(
        ArgumentError("nested factorized product fill requires square output sized to the retained fixed basis"),
    )
    amplitudes = factorized_basis.basis_amplitudes
    triplets = factorized_basis.basis_triplets
    @inbounds for column in 1:nbasis
        xj, yj, zj = triplets[column]
        amplitude_j = include_basis_amplitudes ? amplitudes[column] : 1.0
        for row in 1:column
            xi, yi, zi = triplets[row]
            scale = include_basis_amplitudes ? amplitudes[row] * amplitude_j : 1.0
            value =
                scale *
                operator_x[xi, xj] *
                operator_y[yi, yj] *
                operator_z[zi, zj]
            destination[row, column] = value
            destination[column, row] = value
        end
    end
    return destination
end

function _nested_fill_factorized_sum_of_products!(
    destination::AbstractMatrix{<:Real},
    factorized_basis::_CartesianNestedFactorizedBasis3D,
    terms;
    include_basis_amplitudes::Bool = true,
)
    nbasis = length(factorized_basis.basis_triplets)
    size(destination) == (nbasis, nbasis) || throw(
        ArgumentError("nested factorized sum-of-products fill requires square output sized to the retained fixed basis"),
    )
    amplitudes = factorized_basis.basis_amplitudes
    triplets = factorized_basis.basis_triplets
    @inbounds for column in 1:nbasis
        xj, yj, zj = triplets[column]
        amplitude_j = include_basis_amplitudes ? amplitudes[column] : 1.0
        for row in 1:column
            xi, yi, zi = triplets[row]
            scale = include_basis_amplitudes ? amplitudes[row] * amplitude_j : 1.0
            value = 0.0
            for term in terms
                value +=
                    term[1][xi, xj] *
                    term[2][yi, yj] *
                    term[3][zi, zj]
            end
            value *= scale
            destination[row, column] = value
            destination[column, row] = value
        end
    end
    return destination
end

function _nested_factorized_product_matrix(
    factorized_basis::_CartesianNestedFactorizedBasis3D,
    operator_x::AbstractMatrix{<:Real},
    operator_y::AbstractMatrix{<:Real},
    operator_z::AbstractMatrix{<:Real},
    timing_collector::Union{Nothing,_NestedFixedBlockTimingCollector} = nothing,
    timing_label::Union{Nothing,String} = nothing;
    include_basis_amplitudes::Bool = true,
)
    start_ns = time_ns()
    nbasis = length(factorized_basis.basis_triplets)
    matrix = Matrix{Float64}(undef, nbasis, nbasis)
    _nested_fill_factorized_product_matrix!(
        matrix,
        factorized_basis,
        operator_x,
        operator_y,
        operator_z;
        include_basis_amplitudes = include_basis_amplitudes,
    )
    !isnothing(timing_label) && _nested_record_timing!(timing_collector, timing_label, start_ns)
    return matrix
end

function _nested_factorized_sum_of_products(
    factorized_basis::_CartesianNestedFactorizedBasis3D,
    terms,
    timing_collector::Union{Nothing,_NestedFixedBlockTimingCollector} = nothing,
    timing_label::Union{Nothing,String} = nothing;
    include_basis_amplitudes::Bool = true,
)
    start_ns = time_ns()
    nbasis = length(factorized_basis.basis_triplets)
    matrix = Matrix{Float64}(undef, nbasis, nbasis)
    _nested_fill_factorized_sum_of_products!(
        matrix,
        factorized_basis,
        terms;
        include_basis_amplitudes = include_basis_amplitudes,
    )
    !isnothing(timing_label) && _nested_record_timing!(timing_collector, timing_label, start_ns)
    return matrix
end

function _nested_factorized_gaussian_terms(
    factorized_basis::_CartesianNestedFactorizedBasis3D,
    gaussian_terms_x::Array{Float64,3},
    gaussian_terms_y::Array{Float64,3},
    gaussian_terms_z::Array{Float64,3},
    timing_collector::Union{Nothing,_NestedFixedBlockTimingCollector} = nothing,
)
    start_ns = time_ns()
    axis_term_tables_x = _nested_factorized_axis_term_tables(
        gaussian_terms_x,
        factorized_basis.x_functions,
    )
    axis_term_tables_y = _nested_factorized_axis_term_tables(
        gaussian_terms_y,
        factorized_basis.y_functions,
    )
    axis_term_tables_z = _nested_factorized_axis_term_tables(
        gaussian_terms_z,
        factorized_basis.z_functions,
    )
    nbasis = length(factorized_basis.basis_triplets)
    gaussian_terms = zeros(Float64, size(gaussian_terms_x, 1), nbasis, nbasis)
    _nested_fill_factorized_term_family!(
        gaussian_terms,
        factorized_basis,
        axis_term_tables_x,
        axis_term_tables_y,
        axis_term_tables_z;
        include_basis_amplitudes = true,
    )
    _nested_record_timing!(timing_collector, "packet.gaussian_terms", start_ns)
    return gaussian_terms
end

function _nested_factorized_weight_aware_pair_terms(
    factorized_basis::_CartesianNestedFactorizedBasis3D,
    weights_x::AbstractVector{<:Real},
    weights_y::AbstractVector{<:Real},
    weights_z::AbstractVector{<:Real},
    pair_terms_x::Array{Float64,3},
    pair_terms_y::Array{Float64,3},
    pair_terms_z::Array{Float64,3},
    timing_collector::Union{Nothing,_NestedFixedBlockTimingCollector} = nothing,
)
    start_ns = time_ns()
    axis_weight_x = _nested_factorized_axis_weight_projections(factorized_basis.x_functions, weights_x)
    axis_weight_y = _nested_factorized_axis_weight_projections(factorized_basis.y_functions, weights_y)
    axis_weight_z = _nested_factorized_axis_weight_projections(factorized_basis.z_functions, weights_z)
    basis_weights = Vector{Float64}(undef, length(factorized_basis.basis_triplets))
    @inbounds for index in eachindex(factorized_basis.basis_triplets)
        ix, iy, iz = factorized_basis.basis_triplets[index]
        basis_weights[index] =
            factorized_basis.basis_amplitudes[index] *
            axis_weight_x[ix] *
            axis_weight_y[iy] *
            axis_weight_z[iz]
    end
    all(isfinite, basis_weights) || throw(
        ArgumentError("nested factorized pair contraction requires finite retained fixed integral weights"),
    )
    minimum(basis_weights) > 0.0 || throw(
        ArgumentError("nested factorized pair contraction requires positive retained fixed integral weights"),
    )
    axis_term_tables_x = _nested_factorized_normalized_pair_term_tables(
        _nested_factorized_axis_term_tables(pair_terms_x, factorized_basis.x_functions),
        axis_weight_x,
    )
    axis_term_tables_y = _nested_factorized_normalized_pair_term_tables(
        _nested_factorized_axis_term_tables(pair_terms_y, factorized_basis.y_functions),
        axis_weight_y,
    )
    axis_term_tables_z = _nested_factorized_normalized_pair_term_tables(
        _nested_factorized_axis_term_tables(pair_terms_z, factorized_basis.z_functions),
        axis_weight_z,
    )
    nbasis = length(factorized_basis.basis_triplets)
    pair_terms = zeros(Float64, size(pair_terms_x, 1), nbasis, nbasis)
    _nested_fill_factorized_term_family!(
        pair_terms,
        factorized_basis,
        axis_term_tables_x,
        axis_term_tables_y,
        axis_term_tables_z;
        include_basis_amplitudes = false,
    )
    _nested_record_timing!(timing_collector, "packet.pair_terms", start_ns)
    return (
        weights = basis_weights,
        pair_terms = pair_terms,
    )
end

function _nested_fill_support_product_matrix!(
    workspace::AbstractMatrix{<:Real},
    support_states::AbstractVector{<:NTuple{3,Int}},
    operator_x::AbstractMatrix{<:Real},
    operator_y::AbstractMatrix{<:Real},
    operator_z::AbstractMatrix{<:Real},
)
    nsupport = length(support_states)
    size(workspace) == (nsupport, nsupport) || throw(
        ArgumentError("nested support workspace must have size ($(nsupport), $(nsupport))"),
    )
    @inbounds for row in 1:nsupport
        ix, iy, iz = support_states[row]
        for col in 1:nsupport
            jx, jy, jz = support_states[col]
            workspace[row, col] =
                Float64(operator_x[ix, jx]) *
                Float64(operator_y[iy, jy]) *
                Float64(operator_z[iz, jz])
        end
    end
    return workspace
end

function _nested_fill_support_product_matrix!(
    workspace::AbstractMatrix{<:Real},
    support_axes::_CartesianNestedSupportAxes3D,
    operator_x::AbstractMatrix{<:Real},
    operator_y::AbstractMatrix{<:Real},
    operator_z::AbstractMatrix{<:Real},
)
    nsupport = length(support_axes.x)
    size(workspace) == (nsupport, nsupport) || throw(
        ArgumentError("nested support workspace must have size ($(nsupport), $(nsupport))"),
    )
    xstates = support_axes.x
    ystates = support_axes.y
    zstates = support_axes.z
    @inbounds for row in 1:nsupport
        ix = xstates[row]
        iy = ystates[row]
        iz = zstates[row]
        for col in 1:nsupport
            workspace[row, col] =
                Float64(operator_x[ix, xstates[col]]) *
                Float64(operator_y[iy, ystates[col]]) *
                Float64(operator_z[iz, zstates[col]])
        end
    end
    return workspace
end

function _nested_fill_symmetric_support_product_matrix!(
    workspace::AbstractMatrix{<:Real},
    support_axes::_CartesianNestedSupportAxes3D,
    operator_x::AbstractMatrix{<:Real},
    operator_y::AbstractMatrix{<:Real},
    operator_z::AbstractMatrix{<:Real},
)
    nsupport = length(support_axes.x)
    size(workspace) == (nsupport, nsupport) || throw(
        ArgumentError("nested support workspace must have size ($(nsupport), $(nsupport))"),
    )
    xstates = support_axes.x
    ystates = support_axes.y
    zstates = support_axes.z
    @inbounds for col in 1:nsupport
        jx = xstates[col]
        jy = ystates[col]
        jz = zstates[col]
        for row in 1:col
            value =
                Float64(operator_x[xstates[row], jx]) *
                Float64(operator_y[ystates[row], jy]) *
                Float64(operator_z[zstates[row], jz])
            workspace[row, col] = value
            workspace[col, row] = value
        end
    end
    return workspace
end

function _nested_contract_support_product!(
    destination::AbstractMatrix{<:Real},
    workspace::AbstractMatrix{<:Real},
    contraction_scratch::AbstractMatrix{<:Real},
    support_states::AbstractVector{<:NTuple{3,Int}},
    support_coefficients::AbstractMatrix{<:Real},
    operator_x::AbstractMatrix{<:Real},
    operator_y::AbstractMatrix{<:Real},
    operator_z::AbstractMatrix{<:Real};
    alpha::Float64 = 1.0,
    beta::Float64 = 0.0,
    assume_symmetric::Bool = false,
)
    nsupport = length(support_states)
    nfixed = size(support_coefficients, 2)
    size(destination) == (nfixed, nfixed) || throw(
        ArgumentError("nested support contraction destination must have size ($(nfixed), $(nfixed))"),
    )
    size(workspace) == (nsupport, nsupport) || throw(
        ArgumentError("nested support contraction workspace must have size ($(nsupport), $(nsupport))"),
    )
    size(contraction_scratch) == (nfixed, nsupport) || throw(
        ArgumentError("nested support contraction scratch must have size ($(nfixed), $(nsupport))"),
    )
    if assume_symmetric
        _nested_fill_symmetric_support_product_matrix!(
            workspace,
            _nested_support_axes(support_states),
            operator_x,
            operator_y,
            operator_z,
        )
    else
        _nested_fill_support_product_matrix!(
            workspace,
            support_states,
            operator_x,
            operator_y,
            operator_z,
        )
    end
    mul!(contraction_scratch, transpose(support_coefficients), workspace)
    mul!(destination, contraction_scratch, support_coefficients, alpha, beta)
    return destination
end

function _nested_contract_support_product!(
    destination::AbstractMatrix{<:Real},
    workspace::AbstractMatrix{<:Real},
    contraction_scratch::AbstractMatrix{<:Real},
    support_axes::_CartesianNestedSupportAxes3D,
    support_coefficients::AbstractMatrix{<:Real},
    operator_x::AbstractMatrix{<:Real},
    operator_y::AbstractMatrix{<:Real},
    operator_z::AbstractMatrix{<:Real};
    alpha::Float64 = 1.0,
    beta::Float64 = 0.0,
    assume_symmetric::Bool = false,
)
    nsupport = length(support_axes.x)
    nfixed = size(support_coefficients, 2)
    size(destination) == (nfixed, nfixed) || throw(
        ArgumentError("nested support contraction destination must have size ($(nfixed), $(nfixed))"),
    )
    size(workspace) == (nsupport, nsupport) || throw(
        ArgumentError("nested support contraction workspace must have size ($(nsupport), $(nsupport))"),
    )
    size(contraction_scratch) == (nfixed, nsupport) || throw(
        ArgumentError("nested support contraction scratch must have size ($(nfixed), $(nsupport))"),
    )
    if assume_symmetric
        _nested_fill_symmetric_support_product_matrix!(
            workspace,
            support_axes,
            operator_x,
            operator_y,
            operator_z,
        )
    else
        _nested_fill_support_product_matrix!(
            workspace,
            support_axes,
            operator_x,
            operator_y,
            operator_z,
        )
    end
    mul!(contraction_scratch, transpose(support_coefficients), workspace)
    mul!(destination, contraction_scratch, support_coefficients, alpha, beta)
    return destination
end

function _nested_contract_sum_of_support_products!(
    destination::AbstractMatrix{<:Real},
    workspace::AbstractMatrix{<:Real},
    contraction_scratch::AbstractMatrix{<:Real},
    support_states::AbstractVector{<:NTuple{3,Int}},
    support_coefficients::AbstractMatrix{<:Real},
    terms;
    beta::Float64 = 0.0,
    assume_symmetric::Bool = false,
)
    isempty(terms) && return (iszero(beta) ? fill!(destination, 0.0) : rmul!(destination, beta))
    first_term = true
    for term in terms
        _nested_contract_support_product!(
            destination,
            workspace,
            contraction_scratch,
            support_states,
            support_coefficients,
            term[1],
            term[2],
            term[3];
            alpha = 1.0,
            beta = first_term ? beta : 1.0,
            assume_symmetric = assume_symmetric,
        )
        first_term = false
    end
    return destination
end

function _nested_contract_sum_of_support_products!(
    destination::AbstractMatrix{<:Real},
    workspace::AbstractMatrix{<:Real},
    contraction_scratch::AbstractMatrix{<:Real},
    support_axes::_CartesianNestedSupportAxes3D,
    support_coefficients::AbstractMatrix{<:Real},
    terms;
    beta::Float64 = 0.0,
    assume_symmetric::Bool = false,
)
    isempty(terms) && return (iszero(beta) ? fill!(destination, 0.0) : rmul!(destination, beta))
    first_term = true
    for term in terms
        _nested_contract_support_product!(
            destination,
            workspace,
            contraction_scratch,
            support_axes,
            support_coefficients,
            term[1],
            term[2],
            term[3];
            alpha = 1.0,
            beta = first_term ? beta : 1.0,
            assume_symmetric = assume_symmetric,
        )
        first_term = false
    end
    return destination
end

function _nested_contract_support_term_family!(
    destination_terms::Array{Float64,3},
    workspace::AbstractMatrix{<:Real},
    contraction_scratch::AbstractMatrix{<:Real},
    support_axes::_CartesianNestedSupportAxes3D,
    support_coefficients::AbstractMatrix{<:Real},
    operator_terms_x::Array{Float64,3},
    operator_terms_y::Array{Float64,3},
    operator_terms_z::Array{Float64,3},
    assume_symmetric::Bool = false,
)
    nterms = size(destination_terms, 1)
    nterms == size(operator_terms_x, 1) == size(operator_terms_y, 1) == size(operator_terms_z, 1) || throw(
        ArgumentError("nested support term-family contraction requires matching term counts"),
    )
    for term in 1:nterms
        _nested_contract_support_product!(
            @view(destination_terms[term, :, :]),
            workspace,
            contraction_scratch,
            support_axes,
            support_coefficients,
            @view(operator_terms_x[term, :, :]),
            @view(operator_terms_y[term, :, :]),
            @view(operator_terms_z[term, :, :]);
            beta = 0.0,
            assume_symmetric = assume_symmetric,
        )
    end
    return destination_terms
end

function _nested_symmetrize_matrix!(matrix::AbstractMatrix{<:Real})
    size(matrix, 1) == size(matrix, 2) || throw(
        ArgumentError("nested symmetrization requires a square matrix"),
    )
    n = size(matrix, 1)
    @inbounds for col in 1:n
        for row in 1:col
            value = 0.5 * (Float64(matrix[row, col]) + Float64(matrix[col, row]))
            matrix[row, col] = value
            matrix[col, row] = value
        end
    end
    return matrix
end

function _nested_sum_of_support_products(
    support_states::AbstractVector{<:NTuple{3,Int}},
    terms,
)
    # Retained as an allocating reference helper. Kinetic assembly in
    # `_nested_shell_packet(...)` now sums directly at fixed-block scale via
    # `_nested_contract_sum_of_support_products!`.
    isempty(terms) && return zeros(Float64, length(support_states), length(support_states))
    accumulator = zeros(Float64, length(support_states), length(support_states))
    for term in terms
        accumulator .+= _nested_support_product_matrix(
            support_states,
            term[1],
            term[2],
            term[3],
        )
    end
    return accumulator
end

function _nested_support_weights(
    support_states::AbstractVector{<:NTuple{3,Int}},
    weights_1d::AbstractVector{<:Real},
)
    return _nested_support_weights(support_states, weights_1d, weights_1d, weights_1d)
end

function _nested_support_weights(
    support_states::AbstractVector{<:NTuple{3,Int}},
    weights_x::AbstractVector{<:Real},
    weights_y::AbstractVector{<:Real},
    weights_z::AbstractVector{<:Real},
)
    weights = zeros(Float64, length(support_states))
    for (index, state) in enumerate(support_states)
        ix, iy, iz = state
        weights[index] =
            Float64(weights_x[ix]) *
            Float64(weights_y[iy]) *
            Float64(weights_z[iz])
    end
    return weights
end

function _nested_weight_aware_pair_terms(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    support_states::AbstractVector{<:NTuple{3,Int}},
    support_axes::_CartesianNestedSupportAxes3D,
    support_coefficients::AbstractMatrix{<:Real},
    support_workspace::AbstractMatrix{<:Real},
    contraction_scratch::AbstractMatrix{<:Real},
    timing_collector::Union{Nothing,_NestedFixedBlockTimingCollector} = nothing,
)
    start_ns = time_ns()
    nterms = size(pgdg.pair_factor_terms, 1)
    nfixed = size(support_coefficients, 2)
    support_weights = _nested_support_weights(support_states, pgdg.weights)
    fixed_weights = vec(transpose(support_coefficients) * support_weights)
    all(isfinite, fixed_weights) || throw(
        ArgumentError("nested fixed-block IDA transfer requires finite contracted integral weights"),
    )
    minimum(fixed_weights) > 0.0 || throw(
        ArgumentError("nested fixed-block IDA transfer requires positive contracted integral weights"),
    )
    weighted_support_coefficients = support_coefficients .* reshape(1.0 ./ fixed_weights, 1, :)
    pair_terms = zeros(Float64, nterms, nfixed, nfixed)

    _nested_contract_support_term_family!(
        pair_terms,
        support_workspace,
        contraction_scratch,
        support_axes,
        weighted_support_coefficients,
        pgdg.pair_factor_terms_raw,
        pgdg.pair_factor_terms_raw,
        pgdg.pair_factor_terms_raw,
        true,
    )

    _nested_record_timing!(timing_collector, "packet.pair_terms", start_ns)

    return (
        weights = fixed_weights,
        pair_terms = pair_terms,
    )
end

function _nested_weight_aware_pair_terms(
    bundles::_CartesianNestedAxisBundles3D,
    support_states::AbstractVector{<:NTuple{3,Int}},
    support_axes::_CartesianNestedSupportAxes3D,
    support_coefficients::AbstractMatrix{<:Real},
    support_workspace::AbstractMatrix{<:Real},
    contraction_scratch::AbstractMatrix{<:Real},
    timing_collector::Union{Nothing,_NestedFixedBlockTimingCollector} = nothing,
)
    start_ns = time_ns()
    pgdg_x = _nested_axis_pgdg(bundles, :x)
    pgdg_y = _nested_axis_pgdg(bundles, :y)
    pgdg_z = _nested_axis_pgdg(bundles, :z)
    nterms = size(pgdg_x.pair_factor_terms, 1)
    nterms == size(pgdg_y.pair_factor_terms, 1) == size(pgdg_z.pair_factor_terms, 1) || throw(
        ArgumentError("mixed-axis nested IDA transfer requires the same Gaussian expansion term count on all axes"),
    )
    nfixed = size(support_coefficients, 2)
    support_weights = _nested_support_weights(
        support_states,
        pgdg_x.weights,
        pgdg_y.weights,
        pgdg_z.weights,
    )
    fixed_weights = vec(transpose(support_coefficients) * support_weights)
    all(isfinite, fixed_weights) || throw(
        ArgumentError("mixed-axis nested fixed-block IDA transfer requires finite contracted integral weights"),
    )
    minimum(fixed_weights) > 0.0 || throw(
        ArgumentError("mixed-axis nested fixed-block IDA transfer requires positive contracted integral weights"),
    )
    weighted_support_coefficients = support_coefficients .* reshape(1.0 ./ fixed_weights, 1, :)
    pair_terms = zeros(Float64, nterms, nfixed, nfixed)

    raw_pair_terms_x = pgdg_x.pair_factor_terms_raw
    raw_pair_terms_y = pgdg_y.pair_factor_terms_raw
    raw_pair_terms_z = pgdg_z.pair_factor_terms_raw
    _nested_contract_support_term_family!(
        pair_terms,
        support_workspace,
        contraction_scratch,
        support_axes,
        weighted_support_coefficients,
        raw_pair_terms_x,
        raw_pair_terms_y,
        raw_pair_terms_z,
        true,
    )

    _nested_record_timing!(timing_collector, "packet.pair_terms", start_ns)

    return (
        weights = fixed_weights,
        pair_terms = pair_terms,
    )
end

function _nested_weight_aware_pair_terms(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    support_states::AbstractVector{<:NTuple{3,Int}},
    support_coefficients::AbstractMatrix{<:Real},
)
    nsupport = length(support_states)
    nfixed = size(support_coefficients, 2)
    support_workspace = Matrix{Float64}(undef, nsupport, nsupport)
    contraction_scratch = Matrix{Float64}(undef, nfixed, nsupport)
    support_axes = _nested_support_axes(support_states)
    return _nested_weight_aware_pair_terms(
        pgdg,
        support_states,
        support_axes,
        support_coefficients,
        support_workspace,
        contraction_scratch,
        nothing,
    )
end

function _nested_weight_aware_pair_terms(
    bundles::_CartesianNestedAxisBundles3D,
    support_states::AbstractVector{<:NTuple{3,Int}},
    support_coefficients::AbstractMatrix{<:Real},
)
    nsupport = length(support_states)
    nfixed = size(support_coefficients, 2)
    support_workspace = Matrix{Float64}(undef, nsupport, nsupport)
    contraction_scratch = Matrix{Float64}(undef, nfixed, nsupport)
    support_axes = _nested_support_axes(support_states)
    return _nested_weight_aware_pair_terms(
        bundles,
        support_states,
        support_axes,
        support_coefficients,
        support_workspace,
        contraction_scratch,
        nothing,
    )
end

function _nested_shell_support_indices(
    faces::NTuple{2,_CartesianNestedXYFace3D},
)
    support = Int[]
    seen = Set{Int}()
    for face in faces
        for index in face.support_indices
            if index in seen
                throw(ArgumentError("nested shell assembly requires disjoint face interiors"))
            end
            push!(support, index)
            push!(seen, index)
        end
    end
    sort!(support)
    return support
end

function _nested_shell_packet(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    coefficient_matrix::AbstractMatrix{<:Real},
    support_indices::AbstractVector{Int},
    timing_collector::Union{Nothing,_NestedFixedBlockTimingCollector} = nothing;
    packet_kernel::Symbol = :support_reference,
)
    total_start_ns = time_ns()
    setup_start_ns = time_ns()
    packet_kernel = _nested_normalize_packet_kernel(packet_kernel)
    support_states = [_cartesian_unflat_index(index, size(pgdg.overlap, 1)) for index in support_indices]
    nshell = size(coefficient_matrix, 2)
    nsupport = length(support_states)
    nterms = size(pgdg.gaussian_factor_terms, 1)
    support_axes = packet_kernel == :support_reference ? _nested_support_axes(support_states) : nothing
    support_coefficients =
        packet_kernel == :support_reference ? Matrix{Float64}(coefficient_matrix[support_indices, :]) : nothing
    factorized_basis =
        packet_kernel == :factorized_direct ?
        _nested_extract_factorized_basis(
            coefficient_matrix,
            (size(pgdg.overlap, 1), size(pgdg.overlap, 1), size(pgdg.overlap, 1)),
        ) :
        nothing
    factorized_base_tables =
        packet_kernel == :factorized_direct ?
        _nested_factorized_axis_base_tables(
            factorized_basis.x_functions,
            pgdg.overlap,
            pgdg.kinetic,
            pgdg.position,
            pgdg.x2,
        ) :
        nothing
    support_workspace =
        packet_kernel == :support_reference ? Matrix{Float64}(undef, nsupport, nsupport) : Matrix{Float64}(undef, 0, 0)
    contraction_scratch =
        packet_kernel == :support_reference ? Matrix{Float64}(undef, nshell, nsupport) : Matrix{Float64}(undef, 0, 0)
    _nested_record_timing!(timing_collector, "packet.setup", setup_start_ns)

    overlap =
        packet_kernel == :factorized_direct ?
        _nested_factorized_product_matrix(
            factorized_basis,
            factorized_base_tables.overlap,
            factorized_base_tables.overlap,
            factorized_base_tables.overlap,
            timing_collector,
            "packet.base.overlap",
        ) :
        begin
            overlap_local = Matrix{Float64}(undef, nshell, nshell)
            start_ns = time_ns()
            _nested_contract_support_product!(
                overlap_local,
                support_workspace,
                contraction_scratch,
                support_axes,
                support_coefficients,
                pgdg.overlap,
                pgdg.overlap,
                pgdg.overlap;
                beta = 0.0,
                assume_symmetric = true,
            )
            _nested_record_timing!(timing_collector, "packet.base.overlap", start_ns)
            overlap_local
        end

    kinetic =
        packet_kernel == :factorized_direct ?
        _nested_factorized_sum_of_products(
            factorized_basis,
            (
                (factorized_base_tables.kinetic, factorized_base_tables.overlap, factorized_base_tables.overlap),
                (factorized_base_tables.overlap, factorized_base_tables.kinetic, factorized_base_tables.overlap),
                (factorized_base_tables.overlap, factorized_base_tables.overlap, factorized_base_tables.kinetic),
            ),
            timing_collector,
            "packet.base.kinetic",
        ) :
        begin
            kinetic_local = Matrix{Float64}(undef, nshell, nshell)
            start_ns = time_ns()
            _nested_contract_sum_of_support_products!(
                kinetic_local,
                support_workspace,
                contraction_scratch,
                support_axes,
                support_coefficients,
                (
                    (pgdg.kinetic, pgdg.overlap, pgdg.overlap),
                    (pgdg.overlap, pgdg.kinetic, pgdg.overlap),
                    (pgdg.overlap, pgdg.overlap, pgdg.kinetic),
                );
                beta = 0.0,
                assume_symmetric = true,
            )
            _nested_record_timing!(timing_collector, "packet.base.kinetic", start_ns)
            kinetic_local
        end

    position_x =
        packet_kernel == :factorized_direct ?
        _nested_factorized_product_matrix(
            factorized_basis,
            factorized_base_tables.position,
            factorized_base_tables.overlap,
            factorized_base_tables.overlap,
            timing_collector,
            "packet.base.position_x",
        ) :
        begin
            position_x_local = Matrix{Float64}(undef, nshell, nshell)
            start_ns = time_ns()
            _nested_contract_support_product!(
                position_x_local,
                support_workspace,
                contraction_scratch,
                support_axes,
                support_coefficients,
                pgdg.position,
                pgdg.overlap,
                pgdg.overlap;
                beta = 0.0,
                assume_symmetric = true,
            )
            _nested_record_timing!(timing_collector, "packet.base.position_x", start_ns)
            position_x_local
        end

    position_y =
        packet_kernel == :factorized_direct ?
        _nested_factorized_product_matrix(
            factorized_basis,
            factorized_base_tables.overlap,
            factorized_base_tables.position,
            factorized_base_tables.overlap,
            timing_collector,
            "packet.base.position_y",
        ) :
        begin
            position_y_local = Matrix{Float64}(undef, nshell, nshell)
            start_ns = time_ns()
            _nested_contract_support_product!(
                position_y_local,
                support_workspace,
                contraction_scratch,
                support_axes,
                support_coefficients,
                pgdg.overlap,
                pgdg.position,
                pgdg.overlap;
                beta = 0.0,
                assume_symmetric = true,
            )
            _nested_record_timing!(timing_collector, "packet.base.position_y", start_ns)
            position_y_local
        end

    position_z =
        packet_kernel == :factorized_direct ?
        _nested_factorized_product_matrix(
            factorized_basis,
            factorized_base_tables.overlap,
            factorized_base_tables.overlap,
            factorized_base_tables.position,
            timing_collector,
            "packet.base.position_z",
        ) :
        begin
            position_z_local = Matrix{Float64}(undef, nshell, nshell)
            start_ns = time_ns()
            _nested_contract_support_product!(
                position_z_local,
                support_workspace,
                contraction_scratch,
                support_axes,
                support_coefficients,
                pgdg.overlap,
                pgdg.overlap,
                pgdg.position;
                beta = 0.0,
                assume_symmetric = true,
            )
            _nested_record_timing!(timing_collector, "packet.base.position_z", start_ns)
            position_z_local
        end

    x2_x =
        packet_kernel == :factorized_direct ?
        _nested_factorized_product_matrix(
            factorized_basis,
            factorized_base_tables.x2,
            factorized_base_tables.overlap,
            factorized_base_tables.overlap,
            timing_collector,
            "packet.base.x2_x",
        ) :
        begin
            x2_x_local = Matrix{Float64}(undef, nshell, nshell)
            start_ns = time_ns()
            _nested_contract_support_product!(
                x2_x_local,
                support_workspace,
                contraction_scratch,
                support_axes,
                support_coefficients,
                pgdg.x2,
                pgdg.overlap,
                pgdg.overlap;
                beta = 0.0,
                assume_symmetric = true,
            )
            _nested_record_timing!(timing_collector, "packet.base.x2_x", start_ns)
            x2_x_local
        end

    x2_y =
        packet_kernel == :factorized_direct ?
        _nested_factorized_product_matrix(
            factorized_basis,
            factorized_base_tables.overlap,
            factorized_base_tables.x2,
            factorized_base_tables.overlap,
            timing_collector,
            "packet.base.x2_y",
        ) :
        begin
            x2_y_local = Matrix{Float64}(undef, nshell, nshell)
            start_ns = time_ns()
            _nested_contract_support_product!(
                x2_y_local,
                support_workspace,
                contraction_scratch,
                support_axes,
                support_coefficients,
                pgdg.overlap,
                pgdg.x2,
                pgdg.overlap;
                beta = 0.0,
                assume_symmetric = true,
            )
            _nested_record_timing!(timing_collector, "packet.base.x2_y", start_ns)
            x2_y_local
        end

    x2_z =
        packet_kernel == :factorized_direct ?
        _nested_factorized_product_matrix(
            factorized_basis,
            factorized_base_tables.overlap,
            factorized_base_tables.overlap,
            factorized_base_tables.x2,
            timing_collector,
            "packet.base.x2_z",
        ) :
        begin
            x2_z_local = Matrix{Float64}(undef, nshell, nshell)
            start_ns = time_ns()
            _nested_contract_support_product!(
                x2_z_local,
                support_workspace,
                contraction_scratch,
                support_axes,
                support_coefficients,
                pgdg.overlap,
                pgdg.overlap,
                pgdg.x2;
                beta = 0.0,
                assume_symmetric = true,
            )
            _nested_record_timing!(timing_collector, "packet.base.x2_z", start_ns)
            x2_z_local
        end

    pair_data =
        packet_kernel == :factorized_direct ?
        _nested_factorized_weight_aware_pair_terms(
            factorized_basis,
            pgdg.weights,
            pgdg.weights,
            pgdg.weights,
            pgdg.pair_factor_terms_raw,
            pgdg.pair_factor_terms_raw,
            pgdg.pair_factor_terms_raw,
            timing_collector,
        ) :
        _nested_weight_aware_pair_terms(
            pgdg,
            support_states,
            support_axes,
            support_coefficients,
            support_workspace,
            contraction_scratch,
            timing_collector,
        )
    gaussian_terms =
        packet_kernel == :factorized_direct ?
        _nested_factorized_gaussian_terms(
            factorized_basis,
            pgdg.gaussian_factor_terms,
            pgdg.gaussian_factor_terms,
            pgdg.gaussian_factor_terms,
            timing_collector,
        ) :
        begin
            gaussian_terms_local = zeros(Float64, nterms, nshell, nshell)
            start_ns = time_ns()
            _nested_contract_support_term_family!(
                gaussian_terms_local,
                support_workspace,
                contraction_scratch,
                support_axes,
                support_coefficients,
                pgdg.gaussian_factor_terms,
                pgdg.gaussian_factor_terms,
                pgdg.gaussian_factor_terms,
                true,
            )
            _nested_record_timing!(timing_collector, "packet.gaussian_terms", start_ns)
            gaussian_terms_local
        end
    _nested_record_timing!(timing_collector, "packet.total", total_start_ns)

    return (
        packet = _CartesianNestedShellPacket3D(
            overlap,
            kinetic,
            position_x,
            position_y,
            position_z,
            x2_x,
            x2_y,
            x2_z,
            pair_data.weights,
            gaussian_terms,
            pair_data.pair_terms,
        ),
        support_states = support_states,
    )
end

function _nested_shell_packet(
    bundles::_CartesianNestedAxisBundles3D,
    coefficient_matrix::AbstractMatrix{<:Real},
    support_indices::AbstractVector{Int},
    timing_collector::Union{Nothing,_NestedFixedBlockTimingCollector} = nothing;
    packet_kernel::Symbol = :support_reference,
)
    total_start_ns = time_ns()
    setup_start_ns = time_ns()
    packet_kernel = _nested_normalize_packet_kernel(packet_kernel)
    dims = _nested_axis_lengths(bundles)
    pgdg_x = _nested_axis_pgdg(bundles, :x)
    pgdg_y = _nested_axis_pgdg(bundles, :y)
    pgdg_z = _nested_axis_pgdg(bundles, :z)
    support_states = [_cartesian_unflat_index(index, dims) for index in support_indices]
    nshell = size(coefficient_matrix, 2)
    nsupport = length(support_states)
    nterms = size(pgdg_x.gaussian_factor_terms, 1)
    nterms == size(pgdg_y.gaussian_factor_terms, 1) == size(pgdg_z.gaussian_factor_terms, 1) || throw(
        ArgumentError("mixed-axis nested shell packets require the same Gaussian expansion term count on all axes"),
    )
    support_axes = packet_kernel == :support_reference ? _nested_support_axes(support_states) : nothing
    support_coefficients =
        packet_kernel == :support_reference ? Matrix{Float64}(coefficient_matrix[support_indices, :]) : nothing
    factorized_basis =
        packet_kernel == :factorized_direct ?
        _nested_extract_factorized_basis(coefficient_matrix, dims) :
        nothing
    factorized_base_tables_x =
        packet_kernel == :factorized_direct ?
        _nested_factorized_axis_base_tables(
            factorized_basis.x_functions,
            pgdg_x.overlap,
            pgdg_x.kinetic,
            pgdg_x.position,
            pgdg_x.x2,
        ) :
        nothing
    factorized_base_tables_y =
        packet_kernel == :factorized_direct ?
        _nested_factorized_axis_base_tables(
            factorized_basis.y_functions,
            pgdg_y.overlap,
            pgdg_y.kinetic,
            pgdg_y.position,
            pgdg_y.x2,
        ) :
        nothing
    factorized_base_tables_z =
        packet_kernel == :factorized_direct ?
        _nested_factorized_axis_base_tables(
            factorized_basis.z_functions,
            pgdg_z.overlap,
            pgdg_z.kinetic,
            pgdg_z.position,
            pgdg_z.x2,
        ) :
        nothing
    support_workspace =
        packet_kernel == :support_reference ? Matrix{Float64}(undef, nsupport, nsupport) : Matrix{Float64}(undef, 0, 0)
    contraction_scratch =
        packet_kernel == :support_reference ? Matrix{Float64}(undef, nshell, nsupport) : Matrix{Float64}(undef, 0, 0)
    _nested_record_timing!(timing_collector, "packet.setup", setup_start_ns)

    overlap =
        packet_kernel == :factorized_direct ?
        _nested_factorized_product_matrix(
            factorized_basis,
            factorized_base_tables_x.overlap,
            factorized_base_tables_y.overlap,
            factorized_base_tables_z.overlap,
            timing_collector,
            "packet.base.overlap",
        ) :
        begin
            overlap_local = Matrix{Float64}(undef, nshell, nshell)
            start_ns = time_ns()
            _nested_contract_support_product!(
                overlap_local,
                support_workspace,
                contraction_scratch,
                support_axes,
                support_coefficients,
                pgdg_x.overlap,
                pgdg_y.overlap,
                pgdg_z.overlap;
                beta = 0.0,
                assume_symmetric = true,
            )
            _nested_record_timing!(timing_collector, "packet.base.overlap", start_ns)
            overlap_local
        end

    kinetic =
        packet_kernel == :factorized_direct ?
        _nested_factorized_sum_of_products(
            factorized_basis,
            (
                (factorized_base_tables_x.kinetic, factorized_base_tables_y.overlap, factorized_base_tables_z.overlap),
                (factorized_base_tables_x.overlap, factorized_base_tables_y.kinetic, factorized_base_tables_z.overlap),
                (factorized_base_tables_x.overlap, factorized_base_tables_y.overlap, factorized_base_tables_z.kinetic),
            ),
            timing_collector,
            "packet.base.kinetic",
        ) :
        begin
            kinetic_local = Matrix{Float64}(undef, nshell, nshell)
            start_ns = time_ns()
            _nested_contract_sum_of_support_products!(
                kinetic_local,
                support_workspace,
                contraction_scratch,
                support_axes,
                support_coefficients,
                (
                    (pgdg_x.kinetic, pgdg_y.overlap, pgdg_z.overlap),
                    (pgdg_x.overlap, pgdg_y.kinetic, pgdg_z.overlap),
                    (pgdg_x.overlap, pgdg_y.overlap, pgdg_z.kinetic),
                );
                beta = 0.0,
                assume_symmetric = true,
            )
            _nested_record_timing!(timing_collector, "packet.base.kinetic", start_ns)
            kinetic_local
        end

    position_x =
        packet_kernel == :factorized_direct ?
        _nested_factorized_product_matrix(
            factorized_basis,
            factorized_base_tables_x.position,
            factorized_base_tables_y.overlap,
            factorized_base_tables_z.overlap,
            timing_collector,
            "packet.base.position_x",
        ) :
        begin
            position_x_local = Matrix{Float64}(undef, nshell, nshell)
            start_ns = time_ns()
            _nested_contract_support_product!(
                position_x_local,
                support_workspace,
                contraction_scratch,
                support_axes,
                support_coefficients,
                pgdg_x.position,
                pgdg_y.overlap,
                pgdg_z.overlap;
                beta = 0.0,
                assume_symmetric = true,
            )
            _nested_record_timing!(timing_collector, "packet.base.position_x", start_ns)
            position_x_local
        end

    position_y =
        packet_kernel == :factorized_direct ?
        _nested_factorized_product_matrix(
            factorized_basis,
            factorized_base_tables_x.overlap,
            factorized_base_tables_y.position,
            factorized_base_tables_z.overlap,
            timing_collector,
            "packet.base.position_y",
        ) :
        begin
            position_y_local = Matrix{Float64}(undef, nshell, nshell)
            start_ns = time_ns()
            _nested_contract_support_product!(
                position_y_local,
                support_workspace,
                contraction_scratch,
                support_axes,
                support_coefficients,
                pgdg_x.overlap,
                pgdg_y.position,
                pgdg_z.overlap;
                beta = 0.0,
                assume_symmetric = true,
            )
            _nested_record_timing!(timing_collector, "packet.base.position_y", start_ns)
            position_y_local
        end

    position_z =
        packet_kernel == :factorized_direct ?
        _nested_factorized_product_matrix(
            factorized_basis,
            factorized_base_tables_x.overlap,
            factorized_base_tables_y.overlap,
            factorized_base_tables_z.position,
            timing_collector,
            "packet.base.position_z",
        ) :
        begin
            position_z_local = Matrix{Float64}(undef, nshell, nshell)
            start_ns = time_ns()
            _nested_contract_support_product!(
                position_z_local,
                support_workspace,
                contraction_scratch,
                support_axes,
                support_coefficients,
                pgdg_x.overlap,
                pgdg_y.overlap,
                pgdg_z.position;
                beta = 0.0,
                assume_symmetric = true,
            )
            _nested_record_timing!(timing_collector, "packet.base.position_z", start_ns)
            position_z_local
        end

    x2_x =
        packet_kernel == :factorized_direct ?
        _nested_factorized_product_matrix(
            factorized_basis,
            factorized_base_tables_x.x2,
            factorized_base_tables_y.overlap,
            factorized_base_tables_z.overlap,
            timing_collector,
            "packet.base.x2_x",
        ) :
        begin
            x2_x_local = Matrix{Float64}(undef, nshell, nshell)
            start_ns = time_ns()
            _nested_contract_support_product!(
                x2_x_local,
                support_workspace,
                contraction_scratch,
                support_axes,
                support_coefficients,
                pgdg_x.x2,
                pgdg_y.overlap,
                pgdg_z.overlap;
                beta = 0.0,
                assume_symmetric = true,
            )
            _nested_record_timing!(timing_collector, "packet.base.x2_x", start_ns)
            x2_x_local
        end

    x2_y =
        packet_kernel == :factorized_direct ?
        _nested_factorized_product_matrix(
            factorized_basis,
            factorized_base_tables_x.overlap,
            factorized_base_tables_y.x2,
            factorized_base_tables_z.overlap,
            timing_collector,
            "packet.base.x2_y",
        ) :
        begin
            x2_y_local = Matrix{Float64}(undef, nshell, nshell)
            start_ns = time_ns()
            _nested_contract_support_product!(
                x2_y_local,
                support_workspace,
                contraction_scratch,
                support_axes,
                support_coefficients,
                pgdg_x.overlap,
                pgdg_y.x2,
                pgdg_z.overlap;
                beta = 0.0,
                assume_symmetric = true,
            )
            _nested_record_timing!(timing_collector, "packet.base.x2_y", start_ns)
            x2_y_local
        end

    x2_z =
        packet_kernel == :factorized_direct ?
        _nested_factorized_product_matrix(
            factorized_basis,
            factorized_base_tables_x.overlap,
            factorized_base_tables_y.overlap,
            factorized_base_tables_z.x2,
            timing_collector,
            "packet.base.x2_z",
        ) :
        begin
            x2_z_local = Matrix{Float64}(undef, nshell, nshell)
            start_ns = time_ns()
            _nested_contract_support_product!(
                x2_z_local,
                support_workspace,
                contraction_scratch,
                support_axes,
                support_coefficients,
                pgdg_x.overlap,
                pgdg_y.overlap,
                pgdg_z.x2;
                beta = 0.0,
                assume_symmetric = true,
            )
            _nested_record_timing!(timing_collector, "packet.base.x2_z", start_ns)
            x2_z_local
        end

    pair_data =
        packet_kernel == :factorized_direct ?
        _nested_factorized_weight_aware_pair_terms(
            factorized_basis,
            pgdg_x.weights,
            pgdg_y.weights,
            pgdg_z.weights,
            pgdg_x.pair_factor_terms_raw,
            pgdg_y.pair_factor_terms_raw,
            pgdg_z.pair_factor_terms_raw,
            timing_collector,
        ) :
        _nested_weight_aware_pair_terms(
            bundles,
            support_states,
            support_axes,
            support_coefficients,
            support_workspace,
            contraction_scratch,
            timing_collector,
        )
    gaussian_terms =
        packet_kernel == :factorized_direct ?
        _nested_factorized_gaussian_terms(
            factorized_basis,
            pgdg_x.gaussian_factor_terms,
            pgdg_y.gaussian_factor_terms,
            pgdg_z.gaussian_factor_terms,
            timing_collector,
        ) :
        begin
            gaussian_terms_local = zeros(Float64, nterms, nshell, nshell)
            start_ns = time_ns()
            _nested_contract_support_term_family!(
                gaussian_terms_local,
                support_workspace,
                contraction_scratch,
                support_axes,
                support_coefficients,
                pgdg_x.gaussian_factor_terms,
                pgdg_y.gaussian_factor_terms,
                pgdg_z.gaussian_factor_terms,
                true,
            )
            _nested_record_timing!(timing_collector, "packet.gaussian_terms", start_ns)
            gaussian_terms_local
        end
    _nested_record_timing!(timing_collector, "packet.total", total_start_ns)

    return (
        packet = _CartesianNestedShellPacket3D(
            overlap,
            kinetic,
            position_x,
            position_y,
            position_z,
            x2_x,
            x2_y,
            x2_z,
            pair_data.weights,
            gaussian_terms,
            pair_data.pair_terms,
        ),
        support_states = support_states,
    )
end

# Alg Nested-Face steps 8-9: Assemble one first shell-level fixed space from an
# opposite-face pair and propagate the carried operator packet through the same
# contractions.
# See docs/src/algorithms/cartesian_nested_face_construction.md.
function _nested_xy_shell_pair(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int};
    retain_x::Int,
    retain_y::Int,
)
    n1d = size(pgdg.overlap, 1)
    face_low = _nested_xy_face_product(
        pgdg,
        x_interval,
        y_interval,
        1;
        retain_x = retain_x,
        retain_y = retain_y,
    )
    face_high = _nested_xy_face_product(
        pgdg,
        x_interval,
        y_interval,
        n1d;
        retain_x = retain_x,
        retain_y = retain_y,
    )
    faces = (face_low, face_high)
    coefficient_matrix = hcat(face_low.coefficient_matrix, face_high.coefficient_matrix)
    support_indices = _nested_shell_support_indices(faces)
    shell_data = _nested_shell_packet(pgdg, coefficient_matrix, support_indices)
    return _CartesianNestedXYShell3D(
        faces,
        coefficient_matrix,
        support_indices,
        shell_data.support_states,
        shell_data.packet,
    )
end

function _nested_xy_shell_pair(
    bundle::_MappedOrdinaryGausslet1DBundle,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int};
    retain_x::Int,
    retain_y::Int,
)
    return _nested_xy_shell_pair(
        bundle.pgdg_intermediate,
        x_interval,
        y_interval;
        retain_x = retain_x,
        retain_y = retain_y,
    )
end

function _nested_shell_support_indices(
    faces::AbstractVector{<:_CartesianNestedFace3D},
)
    support = Int[]
    seen = Set{Int}()
    for face in faces
        for index in face.support_indices
            if index in seen
                throw(ArgumentError("nested shell assembly requires disjoint face interiors"))
            end
            push!(support, index)
            push!(seen, index)
        end
    end
    sort!(support)
    return support
end

function _nested_complete_shell_support_indices(
    faces::AbstractVector{<:_CartesianNestedFace3D},
    edges::AbstractVector{<:_CartesianNestedEdge3D},
    corners::AbstractVector{<:_CartesianNestedCorner3D},
)
    support = Int[]
    seen = Set{Int}()
    for piece in Iterators.flatten((faces, edges, corners))
        for index in piece.support_indices
            if index in seen
                throw(ArgumentError("nested complete shell assembly requires disjoint face, edge, and corner supports"))
            end
            push!(support, index)
            push!(seen, index)
        end
    end
    sort!(support)
    return support
end

# Alg Nested-Face steps 8-9: Assemble one first generalized shell-level fixed
# space from all six face interiors and propagate the carried operator packet
# through the same contractions.
# See docs/src/algorithms/cartesian_nested_face_construction.md.
function _nested_rectangular_shell(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int};
    retain_xy::Tuple{Int,Int} = (4, 3),
    retain_xz::Tuple{Int,Int} = (4, 3),
    retain_yz::Tuple{Int,Int} = (4, 3),
    x_fixed::Tuple{Int,Int} = (1, size(pgdg.overlap, 1)),
    y_fixed::Tuple{Int,Int} = (1, size(pgdg.overlap, 1)),
    z_fixed::Tuple{Int,Int} = (1, size(pgdg.overlap, 1)),
    packet_kernel::Symbol = :support_reference,
)
    n1d = size(pgdg.overlap, 1)
    side_x_xy = _nested_doside_1d(pgdg, x_interval, retain_xy[1])
    side_y_xy = _nested_doside_1d(pgdg, y_interval, retain_xy[2])
    side_x_xz = _nested_doside_1d(pgdg, x_interval, retain_xz[1])
    side_z_xz = _nested_doside_1d(pgdg, z_interval, retain_xz[2])
    side_y_yz = _nested_doside_1d(pgdg, y_interval, retain_yz[1])
    side_z_yz = _nested_doside_1d(pgdg, z_interval, retain_yz[2])

    faces = _CartesianNestedFace3D[
        _nested_face_product(:xy, :low, side_x_xy, side_y_xy, z_fixed[1], n1d),
        _nested_face_product(:xy, :high, side_x_xy, side_y_xy, z_fixed[2], n1d),
        _nested_face_product(:xz, :low, side_x_xz, side_z_xz, y_fixed[1], n1d),
        _nested_face_product(:xz, :high, side_x_xz, side_z_xz, y_fixed[2], n1d),
        _nested_face_product(:yz, :low, side_y_yz, side_z_yz, x_fixed[1], n1d),
        _nested_face_product(:yz, :high, side_y_yz, side_z_yz, x_fixed[2], n1d),
    ]

    coefficient_blocks = [face.coefficient_matrix for face in faces]
    coefficient_matrix = hcat(coefficient_blocks...)
    face_column_ranges = UnitRange{Int}[]
    column_start = 1
    for face in faces
        ncols = size(face.coefficient_matrix, 2)
        push!(face_column_ranges, column_start:(column_start + ncols - 1))
        column_start += ncols
    end

    support_indices = _nested_shell_support_indices(faces)
    shell_data = _nested_shell_packet(
        pgdg,
        coefficient_matrix,
        support_indices;
        packet_kernel = packet_kernel,
    )
    return _CartesianNestedShell3D(
        faces,
        face_column_ranges,
        coefficient_matrix,
        support_indices,
        shell_data.support_states,
        shell_data.packet,
    )
end

function _nested_rectangular_shell(
    bundle::_MappedOrdinaryGausslet1DBundle,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int};
    kwargs...,
)
    return _nested_rectangular_shell(bundle.pgdg_intermediate, x_interval, y_interval, z_interval; kwargs...)
end

# Alg Nested-Face completeness step: Build one complete nonrecursive shell
# layer with disjoint face, edge, and corner pieces.
# See docs/cartesian_nested_representation_completeness.md.
function _nested_complete_rectangular_shell(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int};
    retain_xy::Tuple{Int,Int} = (4, 3),
    retain_xz::Tuple{Int,Int} = (4, 3),
    retain_yz::Tuple{Int,Int} = (4, 3),
    retain_x_edge::Int = 3,
    retain_y_edge::Int = 3,
    retain_z_edge::Int = 3,
    x_fixed::Tuple{Int,Int} = (1, size(pgdg.overlap, 1)),
    y_fixed::Tuple{Int,Int} = (1, size(pgdg.overlap, 1)),
    z_fixed::Tuple{Int,Int} = (1, size(pgdg.overlap, 1)),
    timing_collector::Union{Nothing,_NestedFixedBlockTimingCollector} = nothing,
    packet_kernel::Symbol = :support_reference,
)
    prepacket_start_ns = time_ns()
    n1d = size(pgdg.overlap, 1)

    shell_faces = _nested_rectangular_shell(
        pgdg,
        x_interval,
        y_interval,
        z_interval;
        retain_xy = retain_xy,
        retain_xz = retain_xz,
        retain_yz = retain_yz,
        x_fixed = x_fixed,
        y_fixed = y_fixed,
        z_fixed = z_fixed,
        packet_kernel = packet_kernel,
    )

    side_x_edge = _nested_doside_1d(pgdg, x_interval, retain_x_edge)
    side_y_edge = _nested_doside_1d(pgdg, y_interval, retain_y_edge)
    side_z_edge = _nested_doside_1d(pgdg, z_interval, retain_z_edge)

    edges = _CartesianNestedEdge3D[
        _nested_edge_product(:x, (:low, :low), side_x_edge, (y_fixed[1], z_fixed[1]), n1d),
        _nested_edge_product(:x, (:low, :high), side_x_edge, (y_fixed[1], z_fixed[2]), n1d),
        _nested_edge_product(:x, (:high, :low), side_x_edge, (y_fixed[2], z_fixed[1]), n1d),
        _nested_edge_product(:x, (:high, :high), side_x_edge, (y_fixed[2], z_fixed[2]), n1d),
        _nested_edge_product(:y, (:low, :low), side_y_edge, (x_fixed[1], z_fixed[1]), n1d),
        _nested_edge_product(:y, (:low, :high), side_y_edge, (x_fixed[1], z_fixed[2]), n1d),
        _nested_edge_product(:y, (:high, :low), side_y_edge, (x_fixed[2], z_fixed[1]), n1d),
        _nested_edge_product(:y, (:high, :high), side_y_edge, (x_fixed[2], z_fixed[2]), n1d),
        _nested_edge_product(:z, (:low, :low), side_z_edge, (x_fixed[1], y_fixed[1]), n1d),
        _nested_edge_product(:z, (:low, :high), side_z_edge, (x_fixed[1], y_fixed[2]), n1d),
        _nested_edge_product(:z, (:high, :low), side_z_edge, (x_fixed[2], y_fixed[1]), n1d),
        _nested_edge_product(:z, (:high, :high), side_z_edge, (x_fixed[2], y_fixed[2]), n1d),
    ]

    corners = _CartesianNestedCorner3D[
        _nested_corner_piece((:low, :low, :low), (x_fixed[1], y_fixed[1], z_fixed[1]), n1d),
        _nested_corner_piece((:low, :low, :high), (x_fixed[1], y_fixed[1], z_fixed[2]), n1d),
        _nested_corner_piece((:low, :high, :low), (x_fixed[1], y_fixed[2], z_fixed[1]), n1d),
        _nested_corner_piece((:low, :high, :high), (x_fixed[1], y_fixed[2], z_fixed[2]), n1d),
        _nested_corner_piece((:high, :low, :low), (x_fixed[2], y_fixed[1], z_fixed[1]), n1d),
        _nested_corner_piece((:high, :low, :high), (x_fixed[2], y_fixed[1], z_fixed[2]), n1d),
        _nested_corner_piece((:high, :high, :low), (x_fixed[2], y_fixed[2], z_fixed[1]), n1d),
        _nested_corner_piece((:high, :high, :high), (x_fixed[2], y_fixed[2], z_fixed[2]), n1d),
    ]

    coefficient_blocks = Matrix{Float64}[face.coefficient_matrix for face in shell_faces.faces]
    append!(coefficient_blocks, [edge.coefficient_matrix for edge in edges])
    append!(coefficient_blocks, [corner.coefficient_matrix for corner in corners])
    coefficient_matrix = hcat(coefficient_blocks...)

    face_column_ranges = shell_faces.face_column_ranges
    edge_column_ranges = UnitRange{Int}[]
    column_start = size(shell_faces.coefficient_matrix, 2) + 1
    for edge in edges
        ncols = size(edge.coefficient_matrix, 2)
        push!(edge_column_ranges, column_start:(column_start + ncols - 1))
        column_start += ncols
    end
    corner_column_ranges = UnitRange{Int}[]
    for corner in corners
        ncols = size(corner.coefficient_matrix, 2)
        push!(corner_column_ranges, column_start:(column_start + ncols - 1))
        column_start += ncols
    end

    support_indices = _nested_complete_shell_support_indices(shell_faces.faces, edges, corners)
    _nested_record_timing!(timing_collector, "shell_layer.nonpacket", prepacket_start_ns)
    shell_data = _nested_shell_packet(
        pgdg,
        coefficient_matrix,
        support_indices,
        timing_collector;
        packet_kernel = packet_kernel,
    )

    return _CartesianNestedCompleteShell3D(
        shell_faces.faces,
        face_column_ranges,
        edges,
        edge_column_ranges,
        corners,
        corner_column_ranges,
        coefficient_matrix,
        support_indices,
        shell_data.support_states,
        shell_data.packet,
    )
end

function _nested_complete_rectangular_shell(
    bundle::_MappedOrdinaryGausslet1DBundle,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int};
    kwargs...,
)
    return _nested_complete_rectangular_shell(bundle.pgdg_intermediate, x_interval, y_interval, z_interval; kwargs...)
end

function _nested_fixed_block(
    shell::_CartesianNestedShell3D,
    parent_basis::MappedUniformBasis,
)
    packet = shell.packet
    fixed_centers = hcat(
        diag(packet.position_x),
        diag(packet.position_y),
        diag(packet.position_z),
    )
    return _NestedFixedBlock3D(
        parent_basis,
        shell,
        shell.coefficient_matrix,
        shell.support_indices,
        packet.overlap,
        packet.kinetic,
        packet.position_x,
        packet.position_y,
        packet.position_z,
        packet.x2_x,
        packet.x2_y,
        packet.x2_z,
        packet.weights,
        packet.gaussian_terms,
        packet.pair_terms,
        Matrix{Float64}(fixed_centers),
    )
end

function _nested_shell_plus_core(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    shell::_CartesianNestedShell3D,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int},
    ;
    packet_kernel::Symbol = :support_reference,
)
    n1d = size(pgdg.overlap, 1)
    core_indices = _nested_box_support_indices(x_interval, y_interval, z_interval, n1d)
    isempty(intersect(core_indices, shell.support_indices)) || throw(
        ArgumentError("nested shell-plus-core construction requires the direct core block to stay disjoint from the shell-face supports"),
    )
    core_coefficients = _nested_direct_core_coefficients(core_indices, n1d^3)
    coefficient_matrix = hcat(core_coefficients, shell.coefficient_matrix)
    support_indices = sort(vcat(core_indices, shell.support_indices))
    shell_data = _nested_shell_packet(
        pgdg,
        coefficient_matrix,
        support_indices;
        packet_kernel = packet_kernel,
    )
    ncore = length(core_indices)
    shell_column_ranges = UnitRange{Int}[
        (ncore + first(range)):(ncore + last(range)) for range in shell.face_column_ranges
    ]
    return _CartesianNestedShellPlusCore3D(
        shell,
        core_indices,
        [_cartesian_unflat_index(index, n1d) for index in core_indices],
        1:ncore,
        shell_column_ranges,
        coefficient_matrix,
        support_indices,
        shell_data.support_states,
        shell_data.packet,
    )
end

function _nested_shell_plus_core(
    bundle::_MappedOrdinaryGausslet1DBundle,
    shell::_CartesianNestedShell3D,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int},
)
    return _nested_shell_plus_core(bundle.pgdg_intermediate, shell, x_interval, y_interval, z_interval)
end

function _nested_shell_sequence(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D};
    enforce_coverage::Bool = true,
    timing_collector::Union{Nothing,_NestedFixedBlockTimingCollector} = nothing,
    packet_kernel::Symbol = :support_reference,
)
    n1d = size(pgdg.overlap, 1)
    core_indices = _nested_box_support_indices(x_interval, y_interval, z_interval, n1d)
    core_coefficients = _nested_direct_core_coefficients(core_indices, n1d^3)
    return _nested_shell_sequence_from_core_block(
        pgdg,
        core_indices,
        core_coefficients,
        shell_layers;
        enforce_coverage = enforce_coverage,
        timing_collector = timing_collector,
        packet_kernel = packet_kernel,
    )
end

function _nested_shell_sequence_from_core_block(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    core_indices::AbstractVector{Int},
    core_coefficients::AbstractMatrix{<:Real},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D};
    enforce_coverage::Bool = true,
    timing_collector::Union{Nothing,_NestedFixedBlockTimingCollector} = nothing,
    packet_kernel::Symbol = :support_reference,
)
    prepacket_start_ns = time_ns()
    n1d = size(pgdg.overlap, 1)
    support_indices = _nested_sequence_support_indices(core_indices, shell_layers)
    working_box =
        enforce_coverage ?
        _nested_assert_sequence_coverage(core_indices, shell_layers, support_indices, n1d) :
        _nested_sequence_working_box(core_indices, shell_layers, n1d)
    coefficient_blocks = Matrix{Float64}[Matrix{Float64}(core_coefficients)]
    append!(coefficient_blocks, [shell.coefficient_matrix for shell in shell_layers])
    coefficient_matrix = hcat(coefficient_blocks...)
    _nested_record_timing!(timing_collector, "sequence_merge.nonpacket", prepacket_start_ns)
    shell_data = _nested_shell_packet(
        pgdg,
        coefficient_matrix,
        support_indices,
        timing_collector;
        packet_kernel = packet_kernel,
    )

    ncore = size(core_coefficients, 2)
    layer_column_ranges = UnitRange{Int}[]
    column_start = ncore + 1
    for shell in shell_layers
        ncols = size(shell.coefficient_matrix, 2)
        push!(layer_column_ranges, column_start:(column_start + ncols - 1))
        column_start += ncols
    end

    return _CartesianNestedShellSequence3D(
        core_indices,
        [_cartesian_unflat_index(index, n1d) for index in core_indices],
        1:ncore,
        collect(shell_layers),
        layer_column_ranges,
        working_box,
        coefficient_matrix,
        support_indices,
        shell_data.support_states,
        shell_data.packet,
    )
end

function _nested_rectangular_shell(
    bundles::_CartesianNestedAxisBundles3D,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int};
    retain_xy::Tuple{Int,Int} = (4, 3),
    retain_xz::Tuple{Int,Int} = (4, 3),
    retain_yz::Tuple{Int,Int} = (4, 3),
    x_fixed::Tuple{Int,Int} = (1, _nested_axis_lengths(bundles)[1]),
    y_fixed::Tuple{Int,Int} = (1, _nested_axis_lengths(bundles)[2]),
    z_fixed::Tuple{Int,Int} = (1, _nested_axis_lengths(bundles)[3]),
    timing_collector::Union{Nothing,_NestedFixedBlockTimingCollector} = nothing,
    packet_kernel::Symbol = :support_reference,
)
    dims = _nested_axis_lengths(bundles)
    side_x_xy = _nested_doside_1d(_nested_axis_pgdg(bundles, :x), x_interval, retain_xy[1])
    side_y_xy = _nested_doside_1d(_nested_axis_pgdg(bundles, :y), y_interval, retain_xy[2])
    side_x_xz = _nested_doside_1d(_nested_axis_pgdg(bundles, :x), x_interval, retain_xz[1])
    side_z_xz = _nested_doside_1d(_nested_axis_pgdg(bundles, :z), z_interval, retain_xz[2])
    side_y_yz = _nested_doside_1d(_nested_axis_pgdg(bundles, :y), y_interval, retain_yz[1])
    side_z_yz = _nested_doside_1d(_nested_axis_pgdg(bundles, :z), z_interval, retain_yz[2])

    faces = _CartesianNestedFace3D[
        _nested_face_product(:xy, :low, side_x_xy, side_y_xy, z_fixed[1], dims),
        _nested_face_product(:xy, :high, side_x_xy, side_y_xy, z_fixed[2], dims),
        _nested_face_product(:xz, :low, side_x_xz, side_z_xz, y_fixed[1], dims),
        _nested_face_product(:xz, :high, side_x_xz, side_z_xz, y_fixed[2], dims),
        _nested_face_product(:yz, :low, side_y_yz, side_z_yz, x_fixed[1], dims),
        _nested_face_product(:yz, :high, side_y_yz, side_z_yz, x_fixed[2], dims),
    ]

    coefficient_blocks = [face.coefficient_matrix for face in faces]
    coefficient_matrix = hcat(coefficient_blocks...)
    face_column_ranges = UnitRange{Int}[]
    column_start = 1
    for face in faces
        ncols = size(face.coefficient_matrix, 2)
        push!(face_column_ranges, column_start:(column_start + ncols - 1))
        column_start += ncols
    end

    support_indices = _nested_shell_support_indices(faces)
    shell_data = _nested_shell_packet(
        bundles,
        coefficient_matrix,
        support_indices,
        timing_collector;
        packet_kernel = packet_kernel,
    )
    return _CartesianNestedShell3D(
        faces,
        face_column_ranges,
        coefficient_matrix,
        support_indices,
        shell_data.support_states,
        shell_data.packet,
    )
end

function _nested_complete_rectangular_shell(
    bundles::_CartesianNestedAxisBundles3D,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int};
    retain_xy::Tuple{Int,Int} = (4, 3),
    retain_xz::Tuple{Int,Int} = (4, 3),
    retain_yz::Tuple{Int,Int} = (4, 3),
    retain_x_edge::Int = 3,
    retain_y_edge::Int = 3,
    retain_z_edge::Int = 3,
    x_fixed::Tuple{Int,Int} = (1, _nested_axis_lengths(bundles)[1]),
    y_fixed::Tuple{Int,Int} = (1, _nested_axis_lengths(bundles)[2]),
    z_fixed::Tuple{Int,Int} = (1, _nested_axis_lengths(bundles)[3]),
    timing_collector::Union{Nothing,_NestedFixedBlockTimingCollector} = nothing,
    packet_kernel::Symbol = :support_reference,
)
    prepacket_start_ns = time_ns()
    dims = _nested_axis_lengths(bundles)
    shell_faces = _nested_rectangular_shell(
        bundles,
        x_interval,
        y_interval,
        z_interval;
        retain_xy = retain_xy,
        retain_xz = retain_xz,
        retain_yz = retain_yz,
        x_fixed = x_fixed,
        y_fixed = y_fixed,
        z_fixed = z_fixed,
        packet_kernel = packet_kernel,
    )

    side_x_edge = _nested_doside_1d(_nested_axis_pgdg(bundles, :x), x_interval, retain_x_edge)
    side_y_edge = _nested_doside_1d(_nested_axis_pgdg(bundles, :y), y_interval, retain_y_edge)
    side_z_edge = _nested_doside_1d(_nested_axis_pgdg(bundles, :z), z_interval, retain_z_edge)

    edges = _CartesianNestedEdge3D[
        _nested_edge_product(:x, (:low, :low), side_x_edge, (y_fixed[1], z_fixed[1]), dims),
        _nested_edge_product(:x, (:low, :high), side_x_edge, (y_fixed[1], z_fixed[2]), dims),
        _nested_edge_product(:x, (:high, :low), side_x_edge, (y_fixed[2], z_fixed[1]), dims),
        _nested_edge_product(:x, (:high, :high), side_x_edge, (y_fixed[2], z_fixed[2]), dims),
        _nested_edge_product(:y, (:low, :low), side_y_edge, (x_fixed[1], z_fixed[1]), dims),
        _nested_edge_product(:y, (:low, :high), side_y_edge, (x_fixed[1], z_fixed[2]), dims),
        _nested_edge_product(:y, (:high, :low), side_y_edge, (x_fixed[2], z_fixed[1]), dims),
        _nested_edge_product(:y, (:high, :high), side_y_edge, (x_fixed[2], z_fixed[2]), dims),
        _nested_edge_product(:z, (:low, :low), side_z_edge, (x_fixed[1], y_fixed[1]), dims),
        _nested_edge_product(:z, (:low, :high), side_z_edge, (x_fixed[1], y_fixed[2]), dims),
        _nested_edge_product(:z, (:high, :low), side_z_edge, (x_fixed[2], y_fixed[1]), dims),
        _nested_edge_product(:z, (:high, :high), side_z_edge, (x_fixed[2], y_fixed[2]), dims),
    ]

    corners = _CartesianNestedCorner3D[
        _nested_corner_piece((:low, :low, :low), (x_fixed[1], y_fixed[1], z_fixed[1]), dims),
        _nested_corner_piece((:low, :low, :high), (x_fixed[1], y_fixed[1], z_fixed[2]), dims),
        _nested_corner_piece((:low, :high, :low), (x_fixed[1], y_fixed[2], z_fixed[1]), dims),
        _nested_corner_piece((:low, :high, :high), (x_fixed[1], y_fixed[2], z_fixed[2]), dims),
        _nested_corner_piece((:high, :low, :low), (x_fixed[2], y_fixed[1], z_fixed[1]), dims),
        _nested_corner_piece((:high, :low, :high), (x_fixed[2], y_fixed[1], z_fixed[2]), dims),
        _nested_corner_piece((:high, :high, :low), (x_fixed[2], y_fixed[2], z_fixed[1]), dims),
        _nested_corner_piece((:high, :high, :high), (x_fixed[2], y_fixed[2], z_fixed[2]), dims),
    ]

    coefficient_blocks = Matrix{Float64}[face.coefficient_matrix for face in shell_faces.faces]
    append!(coefficient_blocks, [edge.coefficient_matrix for edge in edges])
    append!(coefficient_blocks, [corner.coefficient_matrix for corner in corners])
    coefficient_matrix = hcat(coefficient_blocks...)

    face_column_ranges = shell_faces.face_column_ranges
    edge_column_ranges = UnitRange{Int}[]
    column_start = size(shell_faces.coefficient_matrix, 2) + 1
    for edge in edges
        ncols = size(edge.coefficient_matrix, 2)
        push!(edge_column_ranges, column_start:(column_start + ncols - 1))
        column_start += ncols
    end
    corner_column_ranges = UnitRange{Int}[]
    for corner in corners
        ncols = size(corner.coefficient_matrix, 2)
        push!(corner_column_ranges, column_start:(column_start + ncols - 1))
        column_start += ncols
    end

    support_indices = _nested_complete_shell_support_indices(shell_faces.faces, edges, corners)
    _nested_record_timing!(timing_collector, "shell_layer.nonpacket", prepacket_start_ns)
    shell_data = _nested_shell_packet(
        bundles,
        coefficient_matrix,
        support_indices,
        timing_collector;
        packet_kernel = packet_kernel,
    )

    return _CartesianNestedCompleteShell3D(
        shell_faces.faces,
        face_column_ranges,
        edges,
        edge_column_ranges,
        corners,
        corner_column_ranges,
        coefficient_matrix,
        support_indices,
        shell_data.support_states,
        shell_data.packet,
    )
end

function _nested_shell_sequence(
    bundles::_CartesianNestedAxisBundles3D,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D};
    enforce_coverage::Bool = true,
    timing_collector::Union{Nothing,_NestedFixedBlockTimingCollector} = nothing,
    packet_kernel::Symbol = :support_reference,
)
    dims = _nested_axis_lengths(bundles)
    core_indices = _nested_box_support_indices(x_interval, y_interval, z_interval, dims)
    core_coefficients = _nested_direct_core_coefficients(core_indices, prod(dims))
    return _nested_shell_sequence_from_core_block(
        bundles,
        core_indices,
        core_coefficients,
        shell_layers;
        enforce_coverage = enforce_coverage,
        timing_collector = timing_collector,
        packet_kernel = packet_kernel,
    )
end

function _nested_shell_sequence_from_core_block(
    bundles::_CartesianNestedAxisBundles3D,
    core_indices::AbstractVector{Int},
    core_coefficients::AbstractMatrix{<:Real},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D};
    enforce_coverage::Bool = true,
    timing_collector::Union{Nothing,_NestedFixedBlockTimingCollector} = nothing,
    packet_kernel::Symbol = :support_reference,
)
    prepacket_start_ns = time_ns()
    dims = _nested_axis_lengths(bundles)
    support_indices = _nested_sequence_support_indices(core_indices, shell_layers)
    working_box =
        enforce_coverage ?
        _nested_assert_sequence_coverage(core_indices, shell_layers, support_indices, dims) :
        _nested_sequence_working_box(core_indices, shell_layers, dims)
    coefficient_blocks = Matrix{Float64}[Matrix{Float64}(core_coefficients)]
    append!(coefficient_blocks, [shell.coefficient_matrix for shell in shell_layers])
    coefficient_matrix = hcat(coefficient_blocks...)
    _nested_record_timing!(timing_collector, "sequence_merge.nonpacket", prepacket_start_ns)
    shell_data = _nested_shell_packet(
        bundles,
        coefficient_matrix,
        support_indices,
        timing_collector;
        packet_kernel = packet_kernel,
    )

    ncore = size(core_coefficients, 2)
    layer_column_ranges = UnitRange{Int}[]
    column_start = ncore + 1
    for shell in shell_layers
        ncols = size(shell.coefficient_matrix, 2)
        push!(layer_column_ranges, column_start:(column_start + ncols - 1))
        column_start += ncols
    end

    return _CartesianNestedShellSequence3D(
        collect(core_indices),
        [_cartesian_unflat_index(index, dims) for index in core_indices],
        1:ncore,
        collect(shell_layers),
        layer_column_ranges,
        working_box,
        coefficient_matrix,
        support_indices,
        shell_data.support_states,
        shell_data.packet,
    )
end

function _nested_can_shrink_box(box::NTuple{3,UnitRange{Int}})
    return all(length(interval) >= 3 for interval in box)
end

function _nested_inner_box(box::NTuple{3,UnitRange{Int}})
    _nested_can_shrink_box(box) || throw(
        ArgumentError("nested box shrink requires at least three raw sites along each axis"),
    )
    return (
        (first(box[1]) + 1):(last(box[1]) - 1),
        (first(box[2]) + 1):(last(box[2]) - 1),
        (first(box[3]) + 1):(last(box[3]) - 1),
    )
end

function _nested_complete_shell_sequence_for_box(
    bundles::_CartesianNestedAxisBundles3D,
    box::NTuple{3,UnitRange{Int}};
    nside::Int = 5,
    retain_xy::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_xz::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_yz::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_x_edge::Union{Nothing,Int} = nothing,
    retain_y_edge::Union{Nothing,Int} = nothing,
    retain_z_edge::Union{Nothing,Int} = nothing,
    timing_collector::Union{Nothing,_NestedFixedBlockTimingCollector} = nothing,
    packet_kernel::Symbol = :support_reference,
)
    retention = _nested_resolve_complete_shell_retention(
        nside;
        retain_xy = retain_xy,
        retain_xz = retain_xz,
        retain_yz = retain_yz,
        retain_x_edge = retain_x_edge,
        retain_y_edge = retain_y_edge,
        retain_z_edge = retain_z_edge,
    )
    current_box = box
    shell_layers = _CartesianNestedCompleteShell3D[]
    while minimum(length.(current_box)) > nside
        _nested_can_shrink_box(current_box) || break
        inner_box = _nested_inner_box(current_box)
        push!(
            shell_layers,
            _nested_complete_rectangular_shell(
                bundles,
                inner_box...;
                retain_xy = retention.retain_xy,
                retain_xz = retention.retain_xz,
                retain_yz = retention.retain_yz,
                retain_x_edge = retention.retain_x_edge,
                retain_y_edge = retention.retain_y_edge,
                retain_z_edge = retention.retain_z_edge,
                x_fixed = (first(current_box[1]), last(current_box[1])),
                y_fixed = (first(current_box[2]), last(current_box[2])),
                z_fixed = (first(current_box[3]), last(current_box[3])),
                timing_collector = timing_collector,
                packet_kernel = packet_kernel,
            ),
        )
        current_box = inner_box
    end
    return _nested_shell_sequence(
        bundles,
        current_box...,
        shell_layers,
        timing_collector = timing_collector,
        packet_kernel = packet_kernel,
    )
end

function _nested_source_contract_audit(source::_CartesianNestedBondAlignedDiatomicSource3D)
    return _nested_shell_sequence_contract_audit(source.sequence, _nested_axis_lengths(source.axis_bundles))
end

function _nested_source_contract_audit(source::_CartesianNestedBondAlignedHomonuclearChainSource3D)
    return _nested_shell_sequence_contract_audit(source.sequence, _nested_axis_lengths(source.axis_bundles))
end

function _nested_source_contract_audit(source::_CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D)
    return _nested_shell_sequence_contract_audit(source.sequence, _nested_axis_lengths(source.axis_bundles))
end

function _one_center_atomic_shell_increment(nside::Int)
    nside >= 3 || throw(ArgumentError("one-center atomic shell contract requires nside >= 3"))
    return nside^3 - (nside - 2)^3
end

function _one_center_atomic_complete_shell_retention(nside::Int)
    return _nested_complete_shell_retention_from_nside(nside)
end

function _one_center_atomic_shell_layer_count(working_box_side_count::Int, nside::Int)
    working_box_side_count >= nside || throw(
        ArgumentError("one-center atomic structure diagnostics require working_box_side_count >= nside"),
    )
    current_side = working_box_side_count
    nlayers = 0
    while current_side > nside
        current_side -= 2
        nlayers += 1
    end
    return nlayers, current_side
end

function _one_center_atomic_legacy_profile_working_box(
    parent_side_count::Int,
    working_box::UnitRange{Int},
)
    return _one_center_atomic_legacy_profile_working_box(
        parent_side_count,
        (working_box, working_box, working_box),
    )
end

function _one_center_atomic_legacy_profile_working_box(
    parent_side_count::Int,
    working_box::NTuple{3,UnitRange{Int}},
)
    expected_parent = 1:parent_side_count
    for interval in working_box
        interval == intersect(interval, expected_parent) || throw(
            ArgumentError("one-center atomic legacy profile working box must lie inside 1:$parent_side_count"),
        )
    end
    working_box_sides = Tuple(length.(working_box))
    (working_box_sides[1] == working_box_sides[2] && working_box_sides[2] == working_box_sides[3]) || throw(
        ArgumentError("one-center atomic legacy profile requires a cubic working box"),
    )
    return working_box
end

struct OneCenterAtomicNestedLayerStructure
    layer_index::Int
    face_retained_count::Int
    edge_retained_count::Int
    corner_retained_count::Int
    retained_dimension::Int
end

struct OneCenterAtomicNestedStructureDiagnostics
    parent_side_count::Int
    working_box_side_count::Int
    nside::Int
    core_side_count::Int
    shell_layer_count::Int
    expected_shell_increment::Int
    expected_face_retained_count::Int
    expected_edge_retained_count::Int
    expected_corner_retained_count::Int
    layer_structures::Vector{OneCenterAtomicNestedLayerStructure}
    total_face_retained_count::Int
    total_edge_retained_count::Int
    total_corner_retained_count::Int
    total_expected_gausslet_count::Int
    total_actual_gausslet_count::Int
    layers_match_expected::Bool
end

function _build_one_center_atomic_shell_sequence(
    bundle::_MappedOrdinaryGausslet1DBundle,
    working_box::NTuple{3,UnitRange{Int}};
    nside::Int,
    timing_collector::Union{Nothing,_NestedFixedBlockTimingCollector} = nothing,
    packet_kernel::Symbol = :factorized_direct,
)
    bundles = _CartesianNestedAxisBundles3D(bundle, bundle, bundle)
    retention = _one_center_atomic_complete_shell_retention(nside)
    return _nested_complete_shell_sequence_for_box(
        bundles,
        working_box;
        nside = nside,
        retain_xy = retention.retain_xy,
        retain_xz = retention.retain_xz,
        retain_yz = retention.retain_yz,
        retain_x_edge = retention.retain_x_edge,
        retain_y_edge = retention.retain_y_edge,
        retain_z_edge = retention.retain_z_edge,
        timing_collector = timing_collector,
        packet_kernel = packet_kernel,
    )
end

"""
    build_one_center_atomic_full_parent_shell_sequence(
        bundle::_MappedOrdinaryGausslet1DBundle;
        nside,
    )

Build the canonical one-center atomic nested shell sequence on the full parent
cube of the supplied mapped 1D gausslet bundle.

This helper is the supported atomic one-center backbone:

- it always uses full parent coverage
- it always uses `working_box = (1:n, 1:n, 1:n)`
- it peels complete shells until the direct inner cube reaches `nside`
- it uses the legacy/W&L complete-shell contract
  - shell increment `= nside^3 - (nside - 2)^3`
  - faces retain `(nside - 2) × (nside - 2)`
  - edges retain `nside - 2`
  - corners are carried directly

It should be used in place of any older central-box atomic diagnostic fixture.
"""
function build_one_center_atomic_full_parent_shell_sequence(
    bundle::_MappedOrdinaryGausslet1DBundle;
    nside::Int,
    timing_collector::Union{Nothing,_NestedFixedBlockTimingCollector} = nothing,
)
    n = length(bundle.basis)
    return _build_one_center_atomic_shell_sequence(
        bundle,
        (1:n, 1:n, 1:n);
        nside = nside,
        timing_collector = timing_collector,
    )
end

function build_one_center_atomic_full_parent_shell_sequence(
    basis::MappedUniformBasis;
    exponents::AbstractVector{<:Real} = Float64[],
    center::Real = 0.0,
    gausslet_backend::Symbol = :numerical_reference,
    refinement_levels::Integer = 0,
    kwargs...,
)
    bundle = _mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = exponents,
        center = center,
        backend = gausslet_backend,
        refinement_levels = refinement_levels,
    )
    return build_one_center_atomic_full_parent_shell_sequence(bundle; kwargs...)
end

"""
    build_one_center_atomic_legacy_profile_shell_sequence(
        bundle::_MappedOrdinaryGausslet1DBundle;
        working_box,
        nside,
    )

Build the explicit legacy-profile one-center atomic nested shell sequence on a
chosen inner working box of the supplied parent lattice.

This helper is intentionally separate from the modern canonical full-parent
path. It keeps the same exact complete-shell retention contract:

- shell increment `= nside^3 - (nside - 2)^3`
- faces retain `(nside - 2) × (nside - 2)`
- edges retain `nside - 2`
- corners are carried directly

but applies it on the explicit inner working box supplied by the caller, for
example `(2:28, 2:28, 2:28)` on a `29^3` parent lattice.
"""
function build_one_center_atomic_legacy_profile_shell_sequence(
    bundle::_MappedOrdinaryGausslet1DBundle;
    working_box::Union{UnitRange{Int},NTuple{3,UnitRange{Int}}},
    nside::Int,
    timing_collector::Union{Nothing,_NestedFixedBlockTimingCollector} = nothing,
)
    n = length(bundle.basis)
    normalized_working_box = _one_center_atomic_legacy_profile_working_box(n, working_box)
    return _build_one_center_atomic_shell_sequence(
        bundle,
        normalized_working_box;
        nside = nside,
        timing_collector = timing_collector,
    )
end

function build_one_center_atomic_legacy_profile_shell_sequence(
    basis::MappedUniformBasis;
    exponents::AbstractVector{<:Real} = Float64[],
    center::Real = 0.0,
    gausslet_backend::Symbol = :numerical_reference,
    refinement_levels::Integer = 0,
    kwargs...,
)
    bundle = _mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = exponents,
        center = center,
        backend = gausslet_backend,
        refinement_levels = refinement_levels,
    )
    return build_one_center_atomic_legacy_profile_shell_sequence(bundle; kwargs...)
end

"""
    one_center_atomic_full_parent_fixed_block(
        bundle::_MappedOrdinaryGausslet1DBundle;
        nside,
        kwargs...,
    )

Build the canonical one-center atomic nested fixed block on the full parent
cube of the supplied mapped 1D gausslet bundle.

This is a thin convenience wrapper around
[`build_one_center_atomic_full_parent_shell_sequence`](@ref) followed by
`_nested_fixed_block(...)`.
"""
function one_center_atomic_full_parent_fixed_block(
    bundle::_MappedOrdinaryGausslet1DBundle;
    timing::Union{Bool,Symbol} = false,
    timing_io::IO = stdout,
    kwargs...,
)
    (timing === false || _nested_timing_enabled(timing)) || throw(
        ArgumentError("one-center atomic fixed-block timing must be false, true, or :report"),
    )
    collector = _nested_timing_enabled(timing) ? _nested_new_timing_collector() : nothing
    total_start_ns = time_ns()
    sequence_start_ns = time_ns()
    sequence = build_one_center_atomic_full_parent_shell_sequence(
        bundle;
        timing_collector = collector,
        kwargs...,
    )
    _nested_record_timing!(collector, "fixed_block.sequence_build", sequence_start_ns)
    adapter_start_ns = time_ns()
    fixed_block = _nested_fixed_block(sequence, bundle)
    _nested_record_timing!(collector, "fixed_block.adapter", adapter_start_ns)
    _nested_record_timing!(collector, "fixed_block.total", total_start_ns)
    isnothing(collector) && return fixed_block
    timings = _nested_timing_summary(collector)
    timing === :report && nested_fixed_block_timing_report(timing_io, timings)
    return TimedNestedFixedBlockBuild(fixed_block, timings)
end

function one_center_atomic_full_parent_fixed_block(
    basis::MappedUniformBasis;
    exponents::AbstractVector{<:Real} = Float64[],
    center::Real = 0.0,
    gausslet_backend::Symbol = :numerical_reference,
    refinement_levels::Integer = 0,
    timing::Union{Bool,Symbol} = false,
    timing_io::IO = stdout,
    kwargs...,
)
    (timing === false || _nested_timing_enabled(timing)) || throw(
        ArgumentError("one-center atomic fixed-block timing must be false, true, or :report"),
    )
    collector = _nested_timing_enabled(timing) ? _nested_new_timing_collector() : nothing
    total_start_ns = time_ns()
    bundle_start_ns = time_ns()
    bundle = _mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = exponents,
        center = center,
        backend = gausslet_backend,
        refinement_levels = refinement_levels,
    )
    _nested_record_timing!(collector, "fixed_block.parent_bundle", bundle_start_ns)
    if isnothing(collector)
        return one_center_atomic_full_parent_fixed_block(bundle; timing = false, kwargs...)
    end
    sequence_start_ns = time_ns()
    sequence = build_one_center_atomic_full_parent_shell_sequence(
        bundle;
        timing_collector = collector,
        kwargs...,
    )
    _nested_record_timing!(collector, "fixed_block.sequence_build", sequence_start_ns)
    adapter_start_ns = time_ns()
    fixed_block = _nested_fixed_block(sequence, bundle)
    _nested_record_timing!(collector, "fixed_block.adapter", adapter_start_ns)
    _nested_record_timing!(collector, "fixed_block.total", total_start_ns)
    timings = _nested_timing_summary(collector)
    timing === :report && nested_fixed_block_timing_report(timing_io, timings)
    return TimedNestedFixedBlockBuild(fixed_block, timings)
end

"""
    one_center_atomic_legacy_profile_fixed_block(
        bundle::_MappedOrdinaryGausslet1DBundle;
        working_box,
        nside,
        kwargs...,
    )

Build the legacy-profile one-center atomic nested fixed block on an explicit
inner working box of the supplied parent lattice.

This is a thin convenience wrapper around
[`build_one_center_atomic_legacy_profile_shell_sequence`](@ref) followed by
`_nested_fixed_block(...)`.
"""
function one_center_atomic_legacy_profile_fixed_block(
    bundle::_MappedOrdinaryGausslet1DBundle;
    timing::Union{Bool,Symbol} = false,
    timing_io::IO = stdout,
    kwargs...,
)
    (timing === false || _nested_timing_enabled(timing)) || throw(
        ArgumentError("one-center atomic fixed-block timing must be false, true, or :report"),
    )
    collector = _nested_timing_enabled(timing) ? _nested_new_timing_collector() : nothing
    total_start_ns = time_ns()
    sequence_start_ns = time_ns()
    sequence = build_one_center_atomic_legacy_profile_shell_sequence(
        bundle;
        timing_collector = collector,
        kwargs...,
    )
    _nested_record_timing!(collector, "fixed_block.sequence_build", sequence_start_ns)
    adapter_start_ns = time_ns()
    fixed_block = _nested_fixed_block(sequence, bundle)
    _nested_record_timing!(collector, "fixed_block.adapter", adapter_start_ns)
    _nested_record_timing!(collector, "fixed_block.total", total_start_ns)
    isnothing(collector) && return fixed_block
    timings = _nested_timing_summary(collector)
    timing === :report && nested_fixed_block_timing_report(timing_io, timings)
    return TimedNestedFixedBlockBuild(fixed_block, timings)
end

function one_center_atomic_legacy_profile_fixed_block(
    basis::MappedUniformBasis;
    exponents::AbstractVector{<:Real} = Float64[],
    center::Real = 0.0,
    gausslet_backend::Symbol = :numerical_reference,
    refinement_levels::Integer = 0,
    timing::Union{Bool,Symbol} = false,
    timing_io::IO = stdout,
    kwargs...,
)
    (timing === false || _nested_timing_enabled(timing)) || throw(
        ArgumentError("one-center atomic fixed-block timing must be false, true, or :report"),
    )
    collector = _nested_timing_enabled(timing) ? _nested_new_timing_collector() : nothing
    total_start_ns = time_ns()
    bundle_start_ns = time_ns()
    bundle = _mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = exponents,
        center = center,
        backend = gausslet_backend,
        refinement_levels = refinement_levels,
    )
    _nested_record_timing!(collector, "fixed_block.parent_bundle", bundle_start_ns)
    if isnothing(collector)
        return one_center_atomic_legacy_profile_fixed_block(bundle; timing = false, kwargs...)
    end
    sequence_start_ns = time_ns()
    sequence = build_one_center_atomic_legacy_profile_shell_sequence(
        bundle;
        timing_collector = collector,
        kwargs...,
    )
    _nested_record_timing!(collector, "fixed_block.sequence_build", sequence_start_ns)
    adapter_start_ns = time_ns()
    fixed_block = _nested_fixed_block(sequence, bundle)
    _nested_record_timing!(collector, "fixed_block.adapter", adapter_start_ns)
    _nested_record_timing!(collector, "fixed_block.total", total_start_ns)
    timings = _nested_timing_summary(collector)
    timing === :report && nested_fixed_block_timing_report(timing_io, timings)
    return TimedNestedFixedBlockBuild(fixed_block, timings)
end

function _one_center_atomic_nested_layer_structure(
    shell::_CartesianNestedCompleteShell3D,
    layer_index::Int,
)
    face_retained_count = sum(length, shell.face_column_ranges)
    edge_retained_count = sum(length, shell.edge_column_ranges)
    corner_retained_count = sum(length, shell.corner_column_ranges)
    return OneCenterAtomicNestedLayerStructure(
        layer_index,
        face_retained_count,
        edge_retained_count,
        corner_retained_count,
        size(shell.coefficient_matrix, 2),
    )
end

function one_center_atomic_nested_structure_diagnostics(
    parent_side_count::Int;
    nside::Int,
    working_box_side_count::Int = parent_side_count,
)
    nlayers, core_side_count = _one_center_atomic_shell_layer_count(working_box_side_count, nside)
    retention = _one_center_atomic_complete_shell_retention(nside)
    layer_structures = OneCenterAtomicNestedLayerStructure[
        OneCenterAtomicNestedLayerStructure(
            layer,
            retention.face_retained_count,
            retention.edge_retained_count,
            retention.corner_retained_count,
            retention.shell_increment,
        ) for layer in 1:nlayers
    ]
    total_face_retained_count = nlayers * retention.face_retained_count
    total_edge_retained_count = nlayers * retention.edge_retained_count
    total_corner_retained_count = nlayers * retention.corner_retained_count
    total_expected = core_side_count^3 + nlayers * retention.shell_increment
    return OneCenterAtomicNestedStructureDiagnostics(
        parent_side_count,
        working_box_side_count,
        nside,
        core_side_count,
        nlayers,
        retention.shell_increment,
        retention.face_retained_count,
        retention.edge_retained_count,
        retention.corner_retained_count,
        layer_structures,
        total_face_retained_count,
        total_edge_retained_count,
        total_corner_retained_count,
        total_expected,
        total_expected,
        true,
    )
end

function one_center_atomic_nested_structure_diagnostics(
    sequence::_CartesianNestedShellSequence3D;
    parent_side_count::Int,
    nside::Int,
)
    working_box_sides = Tuple(length.(sequence.working_box))
    working_box_sides[1] == working_box_sides[2] == working_box_sides[3] || throw(
        ArgumentError("one-center atomic structure diagnostics require a cubic working box"),
    )
    core_side_count = round(Int, cbrt(length(sequence.core_indices)))
    layer_structures = OneCenterAtomicNestedLayerStructure[]
    for (layer_index, shell) in enumerate(sequence.shell_layers)
        shell isa _CartesianNestedCompleteShell3D || throw(
            ArgumentError("one-center atomic structure diagnostics require complete shell layers"),
        )
        push!(layer_structures, _one_center_atomic_nested_layer_structure(shell, layer_index))
    end
    total_face_retained_count = sum(layer.face_retained_count for layer in layer_structures)
    total_edge_retained_count = sum(layer.edge_retained_count for layer in layer_structures)
    total_corner_retained_count = sum(layer.corner_retained_count for layer in layer_structures)
    retention = _one_center_atomic_complete_shell_retention(nside)
    total_expected = core_side_count^3 + length(layer_structures) * retention.shell_increment
    return OneCenterAtomicNestedStructureDiagnostics(
        parent_side_count,
        working_box_sides[1],
        nside,
        core_side_count,
        length(layer_structures),
        retention.shell_increment,
        retention.face_retained_count,
        retention.edge_retained_count,
        retention.corner_retained_count,
        layer_structures,
        total_face_retained_count,
        total_edge_retained_count,
        total_corner_retained_count,
        total_expected,
        size(sequence.coefficient_matrix, 2),
        all(
            layer.face_retained_count == retention.face_retained_count &&
            layer.edge_retained_count == retention.edge_retained_count &&
            layer.corner_retained_count == retention.corner_retained_count &&
            layer.retained_dimension == retention.shell_increment
            for layer in layer_structures
        ),
    )
end

function one_center_atomic_nested_structure_diagnostics(
    fixed_block::_NestedFixedBlock3D;
    nside::Int,
)
    parent_side_count = length(fixed_block.parent_basis)
    return one_center_atomic_nested_structure_diagnostics(
        fixed_block.shell;
        parent_side_count = parent_side_count,
        nside = nside,
    )
end

function one_center_atomic_nested_structure_diagnostics(
    bundle::_MappedOrdinaryGausslet1DBundle;
    nside::Int,
)
    sequence = build_one_center_atomic_full_parent_shell_sequence(bundle; nside = nside)
    return one_center_atomic_nested_structure_diagnostics(
        sequence;
        parent_side_count = length(bundle.basis),
        nside = nside,
    )
end

function one_center_atomic_nested_structure_diagnostics(
    basis::MappedUniformBasis;
    exponents::AbstractVector{<:Real} = Float64[],
    center::Real = 0.0,
    gausslet_backend::Symbol = :numerical_reference,
    refinement_levels::Integer = 0,
    nside::Int,
)
    bundle = _mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = exponents,
        center = center,
        backend = gausslet_backend,
        refinement_levels = refinement_levels,
    )
    return one_center_atomic_nested_structure_diagnostics(bundle; nside = nside)
end

function one_center_atomic_nested_structure_report(
    diagnostics::OneCenterAtomicNestedStructureDiagnostics;
    supplement_orbital_count::Union{Nothing,Int} = nothing,
    total_expected_basis_count::Union{Nothing,Int} = nothing,
    total_actual_basis_count::Union{Nothing,Int} = nothing,
    low_one_body_eigenvalues::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    lines = String[
        "parent_side_count = $(diagnostics.parent_side_count)",
        "working_box_side_count = $(diagnostics.working_box_side_count)",
        "nside = $(diagnostics.nside)",
        "core_side_count = $(diagnostics.core_side_count)",
        "shell_layers = $(diagnostics.shell_layer_count)",
        "expected_shell_increment = $(diagnostics.expected_shell_increment)",
        "expected_face_retained_count = $(diagnostics.expected_face_retained_count)",
        "expected_edge_retained_count = $(diagnostics.expected_edge_retained_count)",
        "expected_corner_retained_count = $(diagnostics.expected_corner_retained_count)",
        "total_face_retained_count = $(diagnostics.total_face_retained_count)",
        "total_edge_retained_count = $(diagnostics.total_edge_retained_count)",
        "total_corner_retained_count = $(diagnostics.total_corner_retained_count)",
        "total_expected_gausslet_count = $(diagnostics.total_expected_gausslet_count)",
        "total_actual_gausslet_count = $(diagnostics.total_actual_gausslet_count)",
        "layers_match_expected = $(diagnostics.layers_match_expected)",
    ]
    if !isnothing(supplement_orbital_count)
        push!(lines, "supplement_orbital_count = $(supplement_orbital_count)")
    end
    if !isnothing(total_expected_basis_count)
        push!(lines, "total_expected_basis_count = $(total_expected_basis_count)")
    end
    if !isnothing(total_actual_basis_count)
        push!(lines, "total_actual_basis_count = $(total_actual_basis_count)")
    end
    if !isnothing(low_one_body_eigenvalues)
        push!(lines, "low_one_body_eigenvalues = $(repr(Float64[low_one_body_eigenvalues...]))")
    end
    for layer in diagnostics.layer_structures
        push!(lines, "layer_$(layer.layer_index)_faces = $(layer.face_retained_count)")
        push!(lines, "layer_$(layer.layer_index)_edges = $(layer.edge_retained_count)")
        push!(lines, "layer_$(layer.layer_index)_corners = $(layer.corner_retained_count)")
        push!(lines, "layer_$(layer.layer_index)_retained_dimension = $(layer.retained_dimension)")
    end
    return join(lines, "\n")
end

function _nested_interval_physical_width(
    centers_axis::AbstractVector{<:Real},
    interval::UnitRange{Int},
)
    length(interval) <= 1 && return 0.0
    return Float64(centers_axis[last(interval)] - centers_axis[first(interval)])
end

function _nested_box_physical_widths(
    bundles::_CartesianNestedAxisBundles3D,
    box::NTuple{3,UnitRange{Int}},
)
    return (
        _nested_interval_physical_width(_nested_axis_pgdg(bundles, :x).centers, box[1]),
        _nested_interval_physical_width(_nested_axis_pgdg(bundles, :y).centers, box[2]),
        _nested_interval_physical_width(_nested_axis_pgdg(bundles, :z).centers, box[3]),
    )
end

function _nested_diatomic_midpoint_row_index(
    centers_axis::AbstractVector{<:Real},
    interval::UnitRange{Int},
    midpoint::Real,
)
    length(interval) >= 2 || throw(
        ArgumentError("diatomic midpoint splitting requires at least two raw sites on the bond axis"),
    )
    candidates = collect(first(interval):(last(interval) - 1))
    _, local_index = findmin(abs.(Float64.(centers_axis[candidates]) .- Float64(midpoint)))
    return candidates[local_index]
end

function _nested_diatomic_split_plane_index(
    centers_axis::AbstractVector{<:Real},
    interval::UnitRange{Int},
    midpoint::Real;
    prefer_midpoint_tie_side::Symbol = :left,
    atol::Float64 = 1.0e-12,
    rtol::Float64 = 1.0e-10,
)
    length(interval) >= 2 || throw(
        ArgumentError("diatomic split-plane selection requires at least two raw sites on the bond axis"),
    )
    prefer_midpoint_tie_side in (:left, :right) || throw(
        ArgumentError("diatomic split-plane selection requires prefer_midpoint_tie_side = :left or :right"),
    )
    candidates = collect(first(interval):(last(interval) - 1))
    plane_positions = Float64[
        0.5 * (Float64(centers_axis[index]) + Float64(centers_axis[index + 1])) for index in candidates
    ]
    distances = abs.(plane_positions .- Float64(midpoint))
    minimum_distance = minimum(distances)
    tied = Int[
        candidates[index] for index in eachindex(candidates) if
        isapprox(distances[index], minimum_distance; atol = atol, rtol = rtol)
    ]
    length(tied) >= 1 || throw(ArgumentError("diatomic split-plane selection failed to find a nearest candidate"))
    if prefer_midpoint_tie_side == :left
        return maximum(tied)
    else
        return minimum(tied)
    end
end

function _nested_diatomic_child_boxes(
    box::NTuple{3,UnitRange{Int}},
    bond_axis::Symbol,
    split_index::Int,
)
    axis = bond_axis == :x ? 1 : bond_axis == :y ? 2 : bond_axis == :z ? 3 : 0
    axis != 0 || throw(ArgumentError("bond-axis child-box construction requires bond_axis = :x, :y, or :z"))
    interval = box[axis]
    first(interval) <= split_index < last(interval) || throw(
        ArgumentError("diatomic child-box construction requires the split index to lie strictly inside the working box"),
    )
    left_axis = first(interval):split_index
    right_axis = (split_index + 1):last(interval)
    left_box =
        axis == 1 ? (left_axis, box[2], box[3]) :
        axis == 2 ? (box[1], left_axis, box[3]) :
        (box[1], box[2], left_axis)
    right_box =
        axis == 1 ? (right_axis, box[2], box[3]) :
        axis == 2 ? (box[1], right_axis, box[3]) :
        (box[1], box[2], right_axis)
    return left_box, right_box
end

function _nested_diatomic_midpoint_slab_split(
    box::NTuple{3,UnitRange{Int}},
    bond_axis::Symbol,
    split_index::Int,
)
    axis = bond_axis == :x ? 1 : bond_axis == :y ? 2 : bond_axis == :z ? 3 : 0
    axis != 0 || throw(ArgumentError("bond-axis midpoint-slab construction requires bond_axis = :x, :y, or :z"))
    interval = box[axis]
    first(interval) < split_index < last(interval) || throw(
        ArgumentError("diatomic midpoint-slab construction requires the split index to lie strictly inside the working box"),
    )
    left_axis = first(interval):(split_index - 1)
    slab_axis = split_index:split_index
    right_axis = (split_index + 1):last(interval)
    left_box =
        axis == 1 ? (left_axis, box[2], box[3]) :
        axis == 2 ? (box[1], left_axis, box[3]) :
        (box[1], box[2], left_axis)
    slab_box =
        axis == 1 ? (slab_axis, box[2], box[3]) :
        axis == 2 ? (box[1], slab_axis, box[3]) :
        (box[1], box[2], slab_axis)
    right_box =
        axis == 1 ? (right_axis, box[2], box[3]) :
        axis == 2 ? (box[1], right_axis, box[3]) :
        (box[1], box[2], right_axis)
    return left_box, slab_box, right_box
end

function _nested_direct_box_coefficients(
    dims::NTuple{3,Int},
    support_indices::AbstractVector{Int},
)
    coefficients = zeros(Float64, prod(dims), length(support_indices))
    for (column, index) in pairs(support_indices)
        coefficients[index, column] = 1.0
    end
    return coefficients
end

function _nested_direct_box_coefficients(
    bundles::_CartesianNestedAxisBundles3D,
    box::NTuple{3,UnitRange{Int}},
)
    dims = _nested_axis_lengths(bundles)
    support_indices = _nested_box_support_indices(box..., dims)
    return (
        support_indices = support_indices,
        coefficient_matrix = _nested_direct_box_coefficients(dims, support_indices),
    )
end

function _nested_diatomic_children_are_roughly_cubic(
    bundles::_CartesianNestedAxisBundles3D,
    child_boxes::AbstractVector{<:NTuple{3,UnitRange{Int}}},
    bond_axis::Symbol;
    min_parallel_to_transverse_ratio::Float64 = 0.4,
)
    min_parallel_to_transverse_ratio > 0.0 || throw(
        ArgumentError("diatomic anti-sliver check requires min_parallel_to_transverse_ratio > 0"),
    )
    axis = bond_axis == :x ? 1 : bond_axis == :y ? 2 : bond_axis == :z ? 3 : 0
    axis != 0 || throw(ArgumentError("bond-axis anti-sliver check requires bond_axis = :x, :y, or :z"))
    for child_box in child_boxes
        widths = _nested_box_physical_widths(bundles, child_box)
        parallel = widths[axis]
        transverse = maximum(widths[index] for index in 1:3 if index != axis)
        parallel > 0.0 || return false
        transverse > 0.0 || return false
        parallel >= min_parallel_to_transverse_ratio * transverse || return false
    end
    return true
end

# Alg Nested-Diatomic step 5 and 6: Choose the bond-axis split plane at the
# parent-grid index nearest the midpoint, then reject it if the child boxes are
# too short or too thin in physical coordinates.
# See docs/src/algorithms/cartesian_nested_diatomic_box_policy.md.
function _nested_bond_aligned_diatomic_split_geometry(
    bundles::_CartesianNestedAxisBundles3D,
    parent_box::NTuple{3,UnitRange{Int}},
    working_box::NTuple{3,UnitRange{Int}};
    bond_axis::Symbol = :z,
    midpoint::Real = 0.0,
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    use_midpoint_slab::Bool = true,
    prefer_midpoint_tie_side::Symbol = :left,
)
    axis = bond_axis == :x ? 1 : bond_axis == :y ? 2 : bond_axis == :z ? 3 : 0
    axis != 0 || throw(ArgumentError("diatomic split geometry requires bond_axis = :x, :y, or :z"))
    parallel_interval = working_box[axis]
    parallel_centers = _nested_axis_pgdg(bundles, bond_axis).centers
    use_slab = use_midpoint_slab && isodd(length(parallel_interval))
    split_index = use_slab ?
        _nested_diatomic_midpoint_row_index(parallel_centers, parallel_interval, midpoint) :
        _nested_diatomic_split_plane_index(
            parallel_centers,
            parallel_interval,
            midpoint;
            prefer_midpoint_tie_side = prefer_midpoint_tie_side,
        )
    left_box, midpoint_slab_box, right_box = if use_slab
        _nested_diatomic_midpoint_slab_split(working_box, bond_axis, split_index)
    else
        left_box, right_box = _nested_diatomic_child_boxes(working_box, bond_axis, split_index)
        (left_box, nothing, right_box)
    end
    child_boxes = [left_box, right_box]
    count_eligible =
        length(parallel_interval) > 2 * nside &&
        minimum(length(box[axis]) for box in child_boxes) >= nside
    shape_eligible =
        count_eligible &&
        _nested_diatomic_children_are_roughly_cubic(
            bundles,
            child_boxes,
            bond_axis;
            min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        )
    did_split = count_eligible && shape_eligible
    return _BondAlignedDiatomicSplitGeometry3D(
        parent_box,
        working_box,
        bond_axis,
        Float64(midpoint),
        split_index,
        count_eligible,
        shape_eligible,
        did_split,
        did_split ? midpoint_slab_box : nothing,
        child_boxes,
        [_nested_box_physical_widths(bundles, box) for box in child_boxes],
    )
end

function _nested_bond_aligned_diatomic_source(
    basis,
    bundles::_CartesianNestedAxisBundles3D;
    bond_axis::Symbol = :z,
    midpoint::Real = 0.0,
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    use_midpoint_slab::Bool = true,
    prefer_midpoint_tie_side::Symbol = :left,
    shared_shell_retain_xy::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_xz::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_yz::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_xy::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_xz::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_yz::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_x_edge::Union{Nothing,Int} = nothing,
    retain_y_edge::Union{Nothing,Int} = nothing,
    retain_z_edge::Union{Nothing,Int} = nothing,
)
    child_retention = _nested_resolve_complete_shell_retention(
        nside;
        retain_xy = retain_xy,
        retain_xz = retain_xz,
        retain_yz = retain_yz,
        retain_x_edge = retain_x_edge,
        retain_y_edge = retain_y_edge,
        retain_z_edge = retain_z_edge,
    )
    shared_retention = _nested_resolve_complete_shell_retention(
        nside;
        retain_xy = shared_shell_retain_xy,
        retain_xz = shared_shell_retain_xz,
        retain_yz = shared_shell_retain_yz,
        retain_x_edge = child_retention.retain_x_edge,
        retain_y_edge = child_retention.retain_y_edge,
        retain_z_edge = child_retention.retain_z_edge,
    )
    dims = _nested_axis_lengths(bundles)
    parent_box = (1:dims[1], 1:dims[2], 1:dims[3])
    shared_shell_layers = _CartesianNestedCompleteShell3D[]
    current_box = parent_box
    geometry = _nested_bond_aligned_diatomic_split_geometry(
        bundles,
        parent_box,
        current_box;
        bond_axis = bond_axis,
        midpoint = midpoint,
        nside = nside,
        min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        use_midpoint_slab = use_midpoint_slab,
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
    )

    while true
        parallel_length = length(current_box[bond_axis == :x ? 1 : bond_axis == :y ? 2 : 3])
        if parallel_length <= 2 * nside || minimum(length.(current_box)) <= nside || !_nested_can_shrink_box(current_box)
            break
        end
        inner_box = _nested_inner_box(current_box)
        push!(
            shared_shell_layers,
            _nested_complete_rectangular_shell(
                bundles,
                inner_box...;
                retain_xy = shared_retention.retain_xy,
                retain_xz = shared_retention.retain_xz,
                retain_yz = shared_retention.retain_yz,
                retain_x_edge = shared_retention.retain_x_edge,
                retain_y_edge = shared_retention.retain_y_edge,
                retain_z_edge = shared_retention.retain_z_edge,
                x_fixed = (first(current_box[1]), last(current_box[1])),
                y_fixed = (first(current_box[2]), last(current_box[2])),
                z_fixed = (first(current_box[3]), last(current_box[3])),
            ),
        )
        current_box = inner_box
        geometry = _nested_bond_aligned_diatomic_split_geometry(
            bundles,
            parent_box,
            current_box;
            bond_axis = bond_axis,
            midpoint = midpoint,
            nside = nside,
            min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
            use_midpoint_slab = use_midpoint_slab,
            prefer_midpoint_tie_side = prefer_midpoint_tie_side,
        )
        geometry.did_split && break
    end

    child_sequences = _CartesianNestedShellSequence3D[]
    child_column_ranges = UnitRange{Int}[]
    midpoint_slab_column_range = nothing
    merged_sequence = nothing
    if geometry.did_split
        for child_box in geometry.child_boxes
            push!(
                child_sequences,
                _nested_complete_shell_sequence_for_box(
                    bundles,
                    child_box;
                    nside = nside,
                    retain_xy = child_retention.retain_xy,
                    retain_xz = child_retention.retain_xz,
                    retain_yz = child_retention.retain_yz,
                    retain_x_edge = child_retention.retain_x_edge,
                    retain_y_edge = child_retention.retain_y_edge,
                    retain_z_edge = child_retention.retain_z_edge,
                ),
            )
        end
        core_support_blocks = Vector{Vector{Int}}()
        core_coefficient_blocks = Matrix{Float64}[]
        push!(core_support_blocks, child_sequences[1].support_indices)
        push!(core_coefficient_blocks, child_sequences[1].coefficient_matrix)
        if !isnothing(geometry.shared_midpoint_box)
            slab_data = _nested_direct_box_coefficients(bundles, geometry.shared_midpoint_box)
            push!(core_support_blocks, slab_data.support_indices)
            push!(core_coefficient_blocks, slab_data.coefficient_matrix)
        end
        push!(core_support_blocks, child_sequences[2].support_indices)
        push!(core_coefficient_blocks, child_sequences[2].coefficient_matrix)
        child_support = vcat(core_support_blocks...)
        child_coefficients = hcat(core_coefficient_blocks...)
        merged_sequence = _nested_shell_sequence_from_core_block(
            bundles,
            child_support,
            child_coefficients,
            shared_shell_layers,
        )
        column_start = first(merged_sequence.core_column_range)
        left_columns = size(child_sequences[1].coefficient_matrix, 2)
        push!(child_column_ranges, column_start:(column_start + left_columns - 1))
        column_start = last(child_column_ranges[end]) + 1
        if !isnothing(geometry.shared_midpoint_box)
            slab_columns = prod(length.(geometry.shared_midpoint_box))
            midpoint_slab_column_range = column_start:(column_start + slab_columns - 1)
            column_start = last(midpoint_slab_column_range) + 1
        end
        right_columns = size(child_sequences[2].coefficient_matrix, 2)
        push!(child_column_ranges, column_start:(column_start + right_columns - 1))
    else
        shared_child = _nested_complete_shell_sequence_for_box(
            bundles,
            current_box;
            nside = nside,
            retain_xy = child_retention.retain_xy,
            retain_xz = child_retention.retain_xz,
            retain_yz = child_retention.retain_yz,
            retain_x_edge = child_retention.retain_x_edge,
            retain_y_edge = child_retention.retain_y_edge,
            retain_z_edge = child_retention.retain_z_edge,
        )
        push!(child_sequences, shared_child)
        merged_sequence =
            isempty(shared_shell_layers) ? shared_child :
            _nested_shell_sequence_from_core_block(
                bundles,
                shared_child.support_indices,
                shared_child.coefficient_matrix,
                shared_shell_layers,
            )
        push!(child_column_ranges, merged_sequence.core_column_range)
    end

    return _CartesianNestedBondAlignedDiatomicSource3D(
        basis,
        bundles,
        nside,
        child_retention,
        shared_retention,
        geometry,
        shared_shell_layers,
        child_sequences,
        child_column_ranges,
        midpoint_slab_column_range,
        merged_sequence,
    )
end

function _nested_axis_index(axis::Symbol)
    axis == :x && return 1
    axis == :y && return 2
    axis == :z && return 3
    throw(ArgumentError("nested axis lookup requires axis = :x, :y, or :z"))
end

function _nested_chain_three_child_boxes(
    box::NTuple{3,UnitRange{Int}},
    chain_axis::Symbol,
    left_split_index::Int,
    right_split_index::Int,
)
    axis = _nested_axis_index(chain_axis)
    interval = box[axis]
    first(interval) <= left_split_index < right_split_index < last(interval) || throw(
        ArgumentError("chain three-way child-box construction requires strictly ordered interior split indices"),
    )
    left_axis = first(interval):left_split_index
    middle_axis = (left_split_index + 1):right_split_index
    right_axis = (right_split_index + 1):last(interval)
    left_box =
        axis == 1 ? (left_axis, box[2], box[3]) :
        axis == 2 ? (box[1], left_axis, box[3]) :
        (box[1], box[2], left_axis)
    middle_box =
        axis == 1 ? (middle_axis, box[2], box[3]) :
        axis == 2 ? (box[1], middle_axis, box[3]) :
        (box[1], box[2], middle_axis)
    right_box =
        axis == 1 ? (right_axis, box[2], box[3]) :
        axis == 2 ? (box[1], right_axis, box[3]) :
        (box[1], box[2], right_axis)
    return [left_box, middle_box, right_box]
end

function _nested_odd_chain_policy_thresholds(
    odd_chain_policy::Symbol,
    nside::Int,
    min_parallel_to_transverse_ratio::Float64,
)
    nside >= 1 || throw(ArgumentError("odd-chain nested policy requires nside >= 1"))
    min_parallel_to_transverse_ratio > 0.0 || throw(
        ArgumentError("odd-chain nested policy requires min_parallel_to_transverse_ratio > 0"),
    )
    if odd_chain_policy == :strict_current
        return (
            label = odd_chain_policy,
            outer_parallel_count_min = nside,
            center_parallel_count_min = 1,
            total_parallel_count_min = 2 * nside + 2,
            outer_parallel_to_transverse_ratio_min = min_parallel_to_transverse_ratio,
            center_parallel_to_transverse_ratio_min = min_parallel_to_transverse_ratio,
        )
    elseif odd_chain_policy == :central_ternary_relaxed
        return (
            label = odd_chain_policy,
            outer_parallel_count_min = max(nside - 1, 3),
            center_parallel_count_min = 3,
            total_parallel_count_min = 2 * nside + 1,
            outer_parallel_to_transverse_ratio_min = min_parallel_to_transverse_ratio,
            center_parallel_to_transverse_ratio_min = min(
                min_parallel_to_transverse_ratio,
                0.35,
            ),
        )
    else
        throw(ArgumentError("unknown odd-chain nested policy $(repr(odd_chain_policy)); expected :strict_current or :central_ternary_relaxed"))
    end
end

function _nested_chain_parallel_diagnostics(
    bundles::_CartesianNestedAxisBundles3D,
    child_box::NTuple{3,UnitRange{Int}},
    chain_axis::Symbol,
)
    axis = _nested_axis_index(chain_axis)
    widths = _nested_box_physical_widths(bundles, child_box)
    parallel = widths[axis]
    transverse = maximum(widths[index] for index in 1:3 if index != axis)
    ratio =
        transverse > 0.0 ?
        parallel / transverse :
        (parallel > 0.0 ? Inf : 0.0)
    return (
        widths = widths,
        parallel_count = length(child_box[axis]),
        parallel_to_transverse_ratio = ratio,
    )
end

function _nested_chain_children_are_roughly_cubic(
    bundles::_CartesianNestedAxisBundles3D,
    child_boxes::AbstractVector{<:NTuple{3,UnitRange{Int}}},
    child_nucleus_ranges::AbstractVector{<:UnitRange{Int}},
    chain_axis::Symbol;
    min_parallel_to_transverse_ratio::Float64 = 0.4,
)
    min_parallel_to_transverse_ratio > 0.0 || throw(
        ArgumentError("chain anti-sliver check requires min_parallel_to_transverse_ratio > 0"),
    )
    axis = _nested_axis_index(chain_axis)
    for (child_box, nucleus_range) in zip(child_boxes, child_nucleus_ranges)
        if length(child_box[axis]) == 1 && length(nucleus_range) == 1
            continue
        end
        widths = _nested_box_physical_widths(bundles, child_box)
        parallel = widths[axis]
        transverse = maximum(widths[index] for index in 1:3 if index != axis)
        parallel > 0.0 || return false
        transverse > 0.0 || return false
        parallel >= min_parallel_to_transverse_ratio * transverse || return false
    end
    return true
end

function _nested_chain_binary_candidate(
    bundles::_CartesianNestedAxisBundles3D,
    working_box::NTuple{3,UnitRange{Int}},
    chain_axis::Symbol,
    nucleus_range::UnitRange{Int},
    chain_coordinates::AbstractVector{<:Real},
    split_after::Int;
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    prefer_midpoint_tie_side::Symbol = :left,
)
    midpoint = 0.5 * (Float64(chain_coordinates[split_after]) + Float64(chain_coordinates[split_after + 1]))
    axis = _nested_axis_index(chain_axis)
    parallel_centers = _nested_axis_pgdg(bundles, chain_axis).centers
    split_index = _nested_diatomic_split_plane_index(
        parallel_centers,
        working_box[axis],
        midpoint;
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
    )
    left_box, right_box = _nested_diatomic_child_boxes(working_box, chain_axis, split_index)
    child_boxes = [left_box, right_box]
    child_diagnostics = [
        _nested_chain_parallel_diagnostics(bundles, box, chain_axis) for box in child_boxes
    ]
    nucleus_ranges = [first(nucleus_range):(first(nucleus_range) + split_after - 1), (first(nucleus_range) + split_after):last(nucleus_range)]
    count_eligible =
        length(working_box[axis]) > 2 * nside &&
        minimum(length(box[axis]) for box in child_boxes) >= nside
    shape_eligible =
        count_eligible &&
        _nested_chain_children_are_roughly_cubic(
            bundles,
            child_boxes,
            nucleus_ranges,
            chain_axis;
            min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        )
    return _BondAlignedHomonuclearChainSplitCandidate3D(
        :binary,
        collect(nucleus_ranges),
        [midpoint],
        [split_index],
        child_boxes,
        [entry.widths for entry in child_diagnostics],
        [entry.parallel_count for entry in child_diagnostics],
        [entry.parallel_to_transverse_ratio for entry in child_diagnostics],
        count_eligible,
        shape_eligible,
        count_eligible && shape_eligible,
    )
end

function _nested_chain_ternary_candidate(
    bundles::_CartesianNestedAxisBundles3D,
    working_box::NTuple{3,UnitRange{Int}},
    chain_axis::Symbol,
    nucleus_range::UnitRange{Int},
    chain_coordinates::AbstractVector{<:Real},
    center_index::Int;
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    odd_chain_policy::Symbol = :strict_current,
    prefer_midpoint_tie_side::Symbol = :left,
)
    left_midpoint = 0.5 * (Float64(chain_coordinates[center_index - 1]) + Float64(chain_coordinates[center_index]))
    right_midpoint = 0.5 * (Float64(chain_coordinates[center_index]) + Float64(chain_coordinates[center_index + 1]))
    axis = _nested_axis_index(chain_axis)
    parallel_centers = _nested_axis_pgdg(bundles, chain_axis).centers
    left_split_index = _nested_diatomic_split_plane_index(
        parallel_centers,
        working_box[axis],
        left_midpoint;
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
    )
    right_split_index = _nested_diatomic_split_plane_index(
        parallel_centers,
        working_box[axis],
        right_midpoint;
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
    )
    child_boxes =
        left_split_index < right_split_index ?
        _nested_chain_three_child_boxes(working_box, chain_axis, left_split_index, right_split_index) :
        NTuple{3,UnitRange{Int}}[]
    child_diagnostics = [
        _nested_chain_parallel_diagnostics(bundles, box, chain_axis) for box in child_boxes
    ]
    nucleus_ranges = [
        first(nucleus_range):(first(nucleus_range) + center_index - 2),
        (first(nucleus_range) + center_index - 1):(first(nucleus_range) + center_index - 1),
        (first(nucleus_range) + center_index):last(nucleus_range),
    ]
    thresholds = _nested_odd_chain_policy_thresholds(
        odd_chain_policy,
        nside,
        min_parallel_to_transverse_ratio,
    )
    child_parallel_counts = [entry.parallel_count for entry in child_diagnostics]
    child_parallel_to_transverse_ratios = [
        entry.parallel_to_transverse_ratio for entry in child_diagnostics
    ]
    count_eligible =
        !isempty(child_boxes) &&
        length(working_box[axis]) >= thresholds.total_parallel_count_min &&
        min(child_parallel_counts[1], child_parallel_counts[end]) >= thresholds.outer_parallel_count_min &&
        child_parallel_counts[2] >= thresholds.center_parallel_count_min
    shape_eligible =
        count_eligible &&
        child_parallel_to_transverse_ratios[1] >= thresholds.outer_parallel_to_transverse_ratio_min &&
        child_parallel_to_transverse_ratios[end] >= thresholds.outer_parallel_to_transverse_ratio_min &&
        child_parallel_to_transverse_ratios[2] >= thresholds.center_parallel_to_transverse_ratio_min
    return _BondAlignedHomonuclearChainSplitCandidate3D(
        :ternary,
        collect(nucleus_ranges),
        [left_midpoint, right_midpoint],
        [left_split_index, right_split_index],
        collect(child_boxes),
        [entry.widths for entry in child_diagnostics],
        child_parallel_counts,
        child_parallel_to_transverse_ratios,
        count_eligible,
        shape_eligible,
        count_eligible && shape_eligible,
    )
end

function _nested_chain_split_candidates(
    bundles::_CartesianNestedAxisBundles3D,
    working_box::NTuple{3,UnitRange{Int}},
    chain_axis::Symbol,
    nucleus_range::UnitRange{Int},
    chain_coordinates::AbstractVector{<:Real};
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    odd_chain_policy::Symbol = :strict_current,
    prefer_midpoint_tie_side::Symbol = :left,
)
    natoms = length(chain_coordinates)
    natoms <= 1 && return _BondAlignedHomonuclearChainSplitCandidate3D[]
    candidates = _BondAlignedHomonuclearChainSplitCandidate3D[]
    if iseven(natoms)
        for split_after in 1:(natoms - 1)
            push!(
                candidates,
                _nested_chain_binary_candidate(
                    bundles,
                    working_box,
                    chain_axis,
                    nucleus_range,
                    chain_coordinates,
                    split_after;
                    nside = nside,
                    min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
                    prefer_midpoint_tie_side = prefer_midpoint_tie_side,
                ),
            )
        end
    else
        for center_index in 2:(natoms - 1)
            push!(
                candidates,
                _nested_chain_ternary_candidate(
                    bundles,
                    working_box,
                    chain_axis,
                    nucleus_range,
                    chain_coordinates,
                    center_index;
                    nside = nside,
                    min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
                    odd_chain_policy = odd_chain_policy,
                    prefer_midpoint_tie_side = prefer_midpoint_tie_side,
                ),
            )
        end
    end
    return candidates
end

function _nested_chain_candidate_center(candidate::_BondAlignedHomonuclearChainSplitCandidate3D)
    if candidate.split_kind == :binary
        return candidate.midpoint_values[1]
    else
        return sum(candidate.midpoint_values) / length(candidate.midpoint_values)
    end
end

function _nested_choose_chain_candidate(
    candidates::AbstractVector{_BondAlignedHomonuclearChainSplitCandidate3D},
    chain_coordinates::AbstractVector{<:Real};
    atol::Float64 = 1.0e-12,
    rtol::Float64 = 1.0e-10,
)
    accepted = Int[index for index in eachindex(candidates) if candidates[index].did_split]
    isempty(accepted) && return nothing
    target = 0.5 * (Float64(first(chain_coordinates)) + Float64(last(chain_coordinates)))
    distances = Float64[abs(_nested_chain_candidate_center(candidates[index]) - target) for index in accepted]
    minimum_distance = minimum(distances)
    tied = Int[
        accepted[index] for index in eachindex(accepted) if
        isapprox(distances[index], minimum_distance; atol = atol, rtol = rtol)
    ]
    return minimum(tied)
end

function _nested_chain_node_summary(
    node::_BondAlignedHomonuclearChainNodeGeometry3D,
)
    policy_thresholds = _nested_odd_chain_policy_thresholds(
        node.odd_chain_policy,
        node.nside,
        node.min_parallel_to_transverse_ratio,
    )
    candidate_summaries = [
        (
            split_kind = candidate.split_kind,
            nucleus_ranges = candidate.nucleus_ranges,
            midpoint_values = candidate.midpoint_values,
            split_indices = candidate.split_indices,
            child_boxes = candidate.child_boxes,
            child_physical_widths = candidate.child_physical_widths,
            child_parallel_counts = candidate.child_parallel_counts,
            child_parallel_to_transverse_ratios = candidate.child_parallel_to_transverse_ratios,
            count_eligible = candidate.count_eligible,
            shape_eligible = candidate.shape_eligible,
            did_split = candidate.did_split,
            accepted = index == node.accepted_candidate_index,
        ) for (index, candidate) in pairs(node.candidate_splits)
    ]
    return (
        node_label = node.node_label,
        parent_box = node.parent_box,
        working_box = node.working_box,
        chain_axis = node.chain_axis,
        odd_chain_policy = node.odd_chain_policy,
        odd_chain_policy_thresholds = policy_thresholds,
        nucleus_range = node.nucleus_range,
        chain_coordinates = node.chain_coordinates,
        shared_shell_count = node.shared_shell_count,
        shared_shell_dimensions = node.shared_shell_dimensions,
        accepted_candidate_index = node.accepted_candidate_index,
        did_split = !isnothing(node.accepted_candidate_index),
        local_resolution_warning = node.local_resolution_warning,
        child_count = length(node.child_nodes),
        subtree_fixed_dimension = node.subtree_fixed_dimension,
        candidate_summaries = candidate_summaries,
    )
end

function _nested_chain_collect_node_summaries(
    node::_BondAlignedHomonuclearChainNodeGeometry3D,
)
    summaries = Any[_nested_chain_node_summary(node)]
    for child in node.child_nodes
        append!(summaries, _nested_chain_collect_node_summaries(child))
    end
    return summaries
end

function _nested_bond_aligned_homonuclear_chain_node(
    basis,
    bundles::_CartesianNestedAxisBundles3D,
    parent_box::NTuple{3,UnitRange{Int}},
    nucleus_range::UnitRange{Int},
    node_label::AbstractString;
    chain_axis::Symbol = :z,
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    odd_chain_policy::Symbol = :strict_current,
    prefer_midpoint_tie_side::Symbol = :left,
    retain_xy::Tuple{Int,Int} = (4, 3),
    retain_xz::Tuple{Int,Int} = (4, 3),
    retain_yz::Tuple{Int,Int} = (4, 3),
    retain_x_edge::Int = 3,
    retain_y_edge::Int = 3,
    retain_z_edge::Int = 3,
)
    current_box = parent_box
    shared_shell_layers = _CartesianNestedCompleteShell3D[]
    shared_shell_dimensions = Int[]
    local_coordinates = basis.chain_coordinates[nucleus_range]

    while true
        candidates = _nested_chain_split_candidates(
            bundles,
            current_box,
            chain_axis,
            nucleus_range,
            local_coordinates;
            nside = nside,
            min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
            odd_chain_policy = odd_chain_policy,
            prefer_midpoint_tie_side = prefer_midpoint_tie_side,
        )
        accepted_candidate = _nested_choose_chain_candidate(candidates, local_coordinates)
        (accepted_candidate !== nothing || minimum(length.(current_box)) <= nside || !_nested_can_shrink_box(current_box)) && break

        inner_box = _nested_inner_box(current_box)
        shell = _nested_complete_rectangular_shell(
            bundles,
            inner_box...;
            retain_xy = retain_xy,
            retain_xz = retain_xz,
            retain_yz = retain_yz,
            retain_x_edge = retain_x_edge,
            retain_y_edge = retain_y_edge,
            retain_z_edge = retain_z_edge,
            x_fixed = (first(current_box[1]), last(current_box[1])),
            y_fixed = (first(current_box[2]), last(current_box[2])),
            z_fixed = (first(current_box[3]), last(current_box[3])),
        )
        push!(shared_shell_layers, shell)
        push!(shared_shell_dimensions, size(shell.coefficient_matrix, 2))
        current_box = inner_box
    end

    candidates = _nested_chain_split_candidates(
        bundles,
        current_box,
        chain_axis,
        nucleus_range,
        local_coordinates;
        nside = nside,
        min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        odd_chain_policy = odd_chain_policy,
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
    )
    accepted_candidate = _nested_choose_chain_candidate(candidates, local_coordinates)

    if isnothing(accepted_candidate)
        leaf_sequence = _nested_complete_shell_sequence_for_box(
            bundles,
            current_box;
            nside = nside,
            retain_xy = retain_xy,
            retain_xz = retain_xz,
            retain_yz = retain_yz,
            retain_x_edge = retain_x_edge,
            retain_y_edge = retain_y_edge,
            retain_z_edge = retain_z_edge,
        )
        merged_sequence =
            isempty(shared_shell_layers) ? leaf_sequence :
            _nested_shell_sequence_from_core_block(
                bundles,
                leaf_sequence.support_indices,
                leaf_sequence.coefficient_matrix,
                shared_shell_layers,
            )
        node = _BondAlignedHomonuclearChainNodeGeometry3D(
            String(node_label),
            parent_box,
            current_box,
            chain_axis,
            odd_chain_policy,
            nside,
            min_parallel_to_transverse_ratio,
            nucleus_range,
            collect(local_coordinates),
            length(shared_shell_layers),
            shared_shell_dimensions,
            candidates,
            nothing,
            !isempty(candidates),
            _BondAlignedHomonuclearChainNodeGeometry3D[],
            size(merged_sequence.coefficient_matrix, 2),
        )
        return (
            geometry = node,
            sequence = merged_sequence,
            leaf_sequences = _CartesianNestedShellSequence3D[leaf_sequence],
        )
    end

    candidate = candidates[accepted_candidate]
    child_results = [
        _nested_bond_aligned_homonuclear_chain_node(
            basis,
            bundles,
            candidate.child_boxes[index],
            candidate.nucleus_ranges[index],
            string(node_label, "_", index);
            chain_axis = chain_axis,
            nside = nside,
            min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
            odd_chain_policy = odd_chain_policy,
            prefer_midpoint_tie_side = prefer_midpoint_tie_side,
            retain_xy = retain_xy,
            retain_xz = retain_xz,
            retain_yz = retain_yz,
            retain_x_edge = retain_x_edge,
            retain_y_edge = retain_y_edge,
            retain_z_edge = retain_z_edge,
        ) for index in eachindex(candidate.child_boxes)
    ]
    child_nodes = [result.geometry for result in child_results]
    child_leaf_sequences = _CartesianNestedShellSequence3D[]
    for result in child_results
        append!(child_leaf_sequences, result.leaf_sequences)
    end
    core_support = reduce(vcat, [result.sequence.support_indices for result in child_results])
    core_coefficients = reduce(hcat, [result.sequence.coefficient_matrix for result in child_results])
    merged_sequence =
        isempty(shared_shell_layers) ? _nested_shell_sequence_from_core_block(
            bundles,
            core_support,
            core_coefficients,
            _AbstractCartesianNestedShellLayer3D[];
            enforce_coverage = false,
        ) :
        _nested_shell_sequence_from_core_block(
            bundles,
            core_support,
            core_coefficients,
            shared_shell_layers,
        )
    node = _BondAlignedHomonuclearChainNodeGeometry3D(
        String(node_label),
        parent_box,
        current_box,
        chain_axis,
        odd_chain_policy,
        nside,
        min_parallel_to_transverse_ratio,
        nucleus_range,
        collect(local_coordinates),
        length(shared_shell_layers),
        shared_shell_dimensions,
        candidates,
        accepted_candidate,
        false,
        child_nodes,
        size(merged_sequence.coefficient_matrix, 2),
    )
    return (
        geometry = node,
        sequence = merged_sequence,
        leaf_sequences = child_leaf_sequences,
    )
end

function _nested_bond_aligned_homonuclear_chain_source(
    basis,
    bundles::_CartesianNestedAxisBundles3D;
    chain_axis::Symbol = :z,
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    odd_chain_policy::Symbol = :strict_current,
    prefer_midpoint_tie_side::Symbol = :left,
    retain_xy::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_xz::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_yz::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_x_edge::Union{Nothing,Int} = nothing,
    retain_y_edge::Union{Nothing,Int} = nothing,
    retain_z_edge::Union{Nothing,Int} = nothing,
)
    retention = _nested_resolve_complete_shell_retention(
        nside;
        retain_xy = retain_xy,
        retain_xz = retain_xz,
        retain_yz = retain_yz,
        retain_x_edge = retain_x_edge,
        retain_y_edge = retain_y_edge,
        retain_z_edge = retain_z_edge,
    )
    dims = _nested_axis_lengths(bundles)
    parent_box = (1:dims[1], 1:dims[2], 1:dims[3])
    result = _nested_bond_aligned_homonuclear_chain_node(
        basis,
        bundles,
        parent_box,
        1:length(basis.chain_coordinates),
        "root";
        chain_axis = chain_axis,
        nside = nside,
        min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        odd_chain_policy = odd_chain_policy,
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
        retain_xy = retention.retain_xy,
        retain_xz = retention.retain_xz,
        retain_yz = retention.retain_yz,
        retain_x_edge = retention.retain_x_edge,
        retain_y_edge = retention.retain_y_edge,
        retain_z_edge = retention.retain_z_edge,
    )
    return _CartesianNestedBondAlignedHomonuclearChainSource3D(
        basis,
        bundles,
        retention,
        result.geometry,
        result.leaf_sequences,
        result.sequence,
    )
end

function _nested_square_lattice_planar_diagnostics(
    bundles::_CartesianNestedAxisBundles3D,
    child_box::NTuple{3,UnitRange{Int}},
)
    widths = _nested_box_physical_widths(bundles, child_box)
    planar_x = widths[1]
    planar_y = widths[2]
    longer = max(planar_x, planar_y)
    aspect =
        longer > 0.0 ?
        min(planar_x, planar_y) / longer :
        0.0
    return (
        widths = widths,
        planar_count = (length(child_box[1]), length(child_box[2])),
        in_plane_aspect_ratio = aspect,
    )
end

function _nested_square_lattice_symmetry_preserving(
    child_planar_counts::AbstractVector{<:NTuple{2,Int}},
)
    if length(child_planar_counts) == 2
        return (
            abs(child_planar_counts[1][1] - child_planar_counts[2][1]) <= 1 &&
            abs(child_planar_counts[1][2] - child_planar_counts[2][2]) <= 1
        )
    elseif length(child_planar_counts) == 3
        return (
            abs(child_planar_counts[1][1] - child_planar_counts[3][1]) <= 1 &&
            abs(child_planar_counts[1][2] - child_planar_counts[3][2]) <= 1
        )
    else
        return false
    end
end

function _nested_square_lattice_binary_candidate(
    bundles::_CartesianNestedAxisBundles3D,
    working_box::NTuple{3,UnitRange{Int}},
    split_axis::Symbol,
    x_coordinate_range::UnitRange{Int},
    y_coordinate_range::UnitRange{Int},
    x_coordinates::AbstractVector{<:Real},
    y_coordinates::AbstractVector{<:Real};
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
    prefer_midpoint_tie_side::Symbol = :left,
)
    axis = _nested_axis_index(split_axis)
    coordinate_range = split_axis == :x ? x_coordinate_range : y_coordinate_range
    coordinates = split_axis == :x ? x_coordinates : y_coordinates
    local_count = length(coordinate_range)
    local_count >= 2 || throw(ArgumentError("square-lattice binary split requires at least two coordinates on the chosen axis"))
    split_after = local_count ÷ 2
    midpoint = 0.5 * (
        Float64(coordinates[first(coordinate_range) + split_after - 1]) +
        Float64(coordinates[first(coordinate_range) + split_after])
    )
    split_index = _nested_diatomic_split_plane_index(
        _nested_axis_pgdg(bundles, split_axis).centers,
        working_box[axis],
        midpoint;
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
    )
    child_boxes = collect(_nested_diatomic_child_boxes(working_box, split_axis, split_index))
    split_axis_ranges = [
        first(coordinate_range):(first(coordinate_range) + split_after - 1),
        (first(coordinate_range) + split_after):last(coordinate_range),
    ]
    x_ranges = split_axis == :x ? split_axis_ranges : [x_coordinate_range, x_coordinate_range]
    y_ranges = split_axis == :y ? split_axis_ranges : [y_coordinate_range, y_coordinate_range]
    child_diagnostics = [
        _nested_square_lattice_planar_diagnostics(bundles, box) for box in child_boxes
    ]
    child_planar_counts = [entry.planar_count for entry in child_diagnostics]
    count_eligible =
        length(working_box[axis]) > 2 * nside &&
        minimum(length(box[axis]) for box in child_boxes) >= nside
    shape_eligible =
        count_eligible &&
        all(
            entry.in_plane_aspect_ratio >= min_in_plane_aspect_ratio for
            entry in child_diagnostics
        )
    split_family = split_axis == :x ? :split_x_binary : :split_y_binary
    return _AxisAlignedHomonuclearSquareLatticeSplitCandidate3D(
        split_family,
        split_axis,
        x_ranges,
        y_ranges,
        [midpoint],
        [split_index],
        child_boxes,
        [entry.widths for entry in child_diagnostics],
        child_planar_counts,
        [entry.in_plane_aspect_ratio for entry in child_diagnostics],
        count_eligible,
        shape_eligible,
        _nested_square_lattice_symmetry_preserving(child_planar_counts),
        count_eligible && shape_eligible,
    )
end

function _nested_square_lattice_ternary_candidate(
    bundles::_CartesianNestedAxisBundles3D,
    working_box::NTuple{3,UnitRange{Int}},
    split_axis::Symbol,
    x_coordinate_range::UnitRange{Int},
    y_coordinate_range::UnitRange{Int},
    x_coordinates::AbstractVector{<:Real},
    y_coordinates::AbstractVector{<:Real};
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
    prefer_midpoint_tie_side::Symbol = :left,
)
    axis = _nested_axis_index(split_axis)
    coordinate_range = split_axis == :x ? x_coordinate_range : y_coordinate_range
    coordinates = split_axis == :x ? x_coordinates : y_coordinates
    local_count = length(coordinate_range)
    isodd(local_count) && local_count >= 3 || throw(
        ArgumentError("square-lattice ternary split requires an odd coordinate count of at least three on the chosen axis"),
    )
    center_offset = (local_count + 1) ÷ 2
    center_index = first(coordinate_range) + center_offset - 1
    left_midpoint = 0.5 * (
        Float64(coordinates[center_index - 1]) + Float64(coordinates[center_index])
    )
    right_midpoint = 0.5 * (
        Float64(coordinates[center_index]) + Float64(coordinates[center_index + 1])
    )
    left_split_index = _nested_diatomic_split_plane_index(
        _nested_axis_pgdg(bundles, split_axis).centers,
        working_box[axis],
        left_midpoint;
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
    )
    right_split_index = _nested_diatomic_split_plane_index(
        _nested_axis_pgdg(bundles, split_axis).centers,
        working_box[axis],
        right_midpoint;
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
    )
    child_boxes =
        left_split_index < right_split_index ?
        _nested_chain_three_child_boxes(
            working_box,
            split_axis,
            left_split_index,
            right_split_index,
        ) :
        NTuple{3,UnitRange{Int}}[]
    split_axis_ranges = [
        first(coordinate_range):(center_index - 1),
        center_index:center_index,
        (center_index + 1):last(coordinate_range),
    ]
    x_ranges = split_axis == :x ? split_axis_ranges : [x_coordinate_range, x_coordinate_range, x_coordinate_range]
    y_ranges = split_axis == :y ? split_axis_ranges : [y_coordinate_range, y_coordinate_range, y_coordinate_range]
    child_diagnostics = [
        _nested_square_lattice_planar_diagnostics(bundles, box) for box in child_boxes
    ]
    child_planar_counts = [entry.planar_count for entry in child_diagnostics]
    outer_count_min = max(nside - 1, 3)
    center_count_min = 3
    count_eligible =
        !isempty(child_boxes) &&
        length(working_box[axis]) >= 2 * nside + 1 &&
        min(length(child_boxes[1][axis]), length(child_boxes[end][axis])) >= outer_count_min &&
        length(child_boxes[2][axis]) >= center_count_min
    shape_eligible =
        count_eligible &&
        all(
            entry.in_plane_aspect_ratio >= min_in_plane_aspect_ratio for
            entry in child_diagnostics
        )
    split_family = split_axis == :x ? :split_x_ternary : :split_y_ternary
    return _AxisAlignedHomonuclearSquareLatticeSplitCandidate3D(
        split_family,
        split_axis,
        x_ranges,
        y_ranges,
        [left_midpoint, right_midpoint],
        [left_split_index, right_split_index],
        collect(child_boxes),
        [entry.widths for entry in child_diagnostics],
        child_planar_counts,
        [entry.in_plane_aspect_ratio for entry in child_diagnostics],
        count_eligible,
        shape_eligible,
        _nested_square_lattice_symmetry_preserving(child_planar_counts),
        count_eligible && shape_eligible,
    )
end

function _nested_square_lattice_split_candidates(
    bundles::_CartesianNestedAxisBundles3D,
    working_box::NTuple{3,UnitRange{Int}},
    x_coordinate_range::UnitRange{Int},
    y_coordinate_range::UnitRange{Int},
    x_coordinates::AbstractVector{<:Real},
    y_coordinates::AbstractVector{<:Real};
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
    prefer_midpoint_tie_side::Symbol = :left,
)
    candidates = _AxisAlignedHomonuclearSquareLatticeSplitCandidate3D[]
    nx_sites = length(x_coordinate_range)
    ny_sites = length(y_coordinate_range)

    if nx_sites > 1
        push!(
            candidates,
            iseven(nx_sites) ?
            _nested_square_lattice_binary_candidate(
                bundles,
                working_box,
                :x,
                x_coordinate_range,
                y_coordinate_range,
                x_coordinates,
                y_coordinates;
                nside = nside,
                min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
                prefer_midpoint_tie_side = prefer_midpoint_tie_side,
            ) :
            _nested_square_lattice_ternary_candidate(
                bundles,
                working_box,
                :x,
                x_coordinate_range,
                y_coordinate_range,
                x_coordinates,
                y_coordinates;
                nside = nside,
                min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
                prefer_midpoint_tie_side = prefer_midpoint_tie_side,
            ),
        )
    end

    if ny_sites > 1
        push!(
            candidates,
            iseven(ny_sites) ?
            _nested_square_lattice_binary_candidate(
                bundles,
                working_box,
                :y,
                x_coordinate_range,
                y_coordinate_range,
                x_coordinates,
                y_coordinates;
                nside = nside,
                min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
                prefer_midpoint_tie_side = prefer_midpoint_tie_side,
            ) :
            _nested_square_lattice_ternary_candidate(
                bundles,
                working_box,
                :y,
                x_coordinate_range,
                y_coordinate_range,
                x_coordinates,
                y_coordinates;
                nside = nside,
                min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
                prefer_midpoint_tie_side = prefer_midpoint_tie_side,
            ),
        )
    end

    return candidates
end

function _nested_square_lattice_candidate_center(
    candidate::_AxisAlignedHomonuclearSquareLatticeSplitCandidate3D,
    x_coordinates::AbstractVector{<:Real},
    y_coordinates::AbstractVector{<:Real},
)
    x_mid = 0.5 * (Float64(first(x_coordinates)) + Float64(last(x_coordinates)))
    y_mid = 0.5 * (Float64(first(y_coordinates)) + Float64(last(y_coordinates)))
    split_center = sum(candidate.split_values) / length(candidate.split_values)
    return candidate.split_axis == :x ? (split_center, y_mid) : (x_mid, split_center)
end

function _nested_choose_square_lattice_candidate(
    candidates::AbstractVector{_AxisAlignedHomonuclearSquareLatticeSplitCandidate3D},
    x_coordinates::AbstractVector{<:Real},
    y_coordinates::AbstractVector{<:Real};
    atol::Float64 = 1.0e-12,
    rtol::Float64 = 1.0e-10,
)
    accepted = Int[index for index in eachindex(candidates) if candidates[index].did_split]
    isempty(accepted) && return nothing
    target = (
        0.5 * (Float64(first(x_coordinates)) + Float64(last(x_coordinates))),
        0.5 * (Float64(first(y_coordinates)) + Float64(last(y_coordinates))),
    )
    distances = Float64[
        hypot(
            _nested_square_lattice_candidate_center(candidates[index], x_coordinates, y_coordinates)[1] - target[1],
            _nested_square_lattice_candidate_center(candidates[index], x_coordinates, y_coordinates)[2] - target[2],
        ) for index in accepted
    ]
    minimum_distance = minimum(distances)
    tied = Int[
        accepted[index] for index in eachindex(accepted) if
        isapprox(distances[index], minimum_distance; atol = atol, rtol = rtol)
    ]
    return minimum(tied)
end

function _nested_square_lattice_node_summary(
    node::_AxisAlignedHomonuclearSquareLatticeNodeGeometry3D,
)
    candidate_summaries = [
        (
            split_family = candidate.split_family,
            split_axis = candidate.split_axis,
            x_coordinate_ranges = candidate.x_coordinate_ranges,
            y_coordinate_ranges = candidate.y_coordinate_ranges,
            split_values = candidate.split_values,
            split_indices = candidate.split_indices,
            child_boxes = candidate.child_boxes,
            child_physical_widths = candidate.child_physical_widths,
            child_planar_counts = candidate.child_planar_counts,
            child_in_plane_aspect_ratios = candidate.child_in_plane_aspect_ratios,
            count_eligible = candidate.count_eligible,
            shape_eligible = candidate.shape_eligible,
            symmetry_preserving = candidate.symmetry_preserving,
            did_split = candidate.did_split,
            accepted = index == node.accepted_candidate_index,
        ) for (index, candidate) in pairs(node.candidate_splits)
    ]
    return (
        node_label = node.node_label,
        parent_box = node.parent_box,
        working_box = node.working_box,
        lattice_size = node.lattice_size,
        min_in_plane_aspect_ratio = node.min_in_plane_aspect_ratio,
        x_coordinate_range = node.x_coordinate_range,
        y_coordinate_range = node.y_coordinate_range,
        x_coordinates = node.x_coordinates,
        y_coordinates = node.y_coordinates,
        shared_shell_count = node.shared_shell_count,
        shared_shell_dimensions = node.shared_shell_dimensions,
        accepted_candidate_index = node.accepted_candidate_index,
        did_split = !isnothing(node.accepted_candidate_index),
        local_resolution_warning = node.local_resolution_warning,
        child_count = length(node.child_nodes),
        subtree_fixed_dimension = node.subtree_fixed_dimension,
        candidate_summaries = candidate_summaries,
    )
end

function _nested_square_lattice_collect_node_summaries(
    node::_AxisAlignedHomonuclearSquareLatticeNodeGeometry3D,
)
    summaries = Any[_nested_square_lattice_node_summary(node)]
    for child in node.child_nodes
        append!(summaries, _nested_square_lattice_collect_node_summaries(child))
    end
    return summaries
end

function _nested_axis_aligned_homonuclear_square_lattice_node(
    basis,
    bundles::_CartesianNestedAxisBundles3D,
    parent_box::NTuple{3,UnitRange{Int}},
    x_coordinate_range::UnitRange{Int},
    y_coordinate_range::UnitRange{Int},
    node_label::AbstractString;
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
    prefer_midpoint_tie_side::Symbol = :left,
    retain_xy::Tuple{Int,Int} = (4, 3),
    retain_xz::Tuple{Int,Int} = (4, 3),
    retain_yz::Tuple{Int,Int} = (4, 3),
    retain_x_edge::Int = 3,
    retain_y_edge::Int = 3,
    retain_z_edge::Int = 3,
)
    current_box = parent_box
    shared_shell_layers = _CartesianNestedCompleteShell3D[]
    shared_shell_dimensions = Int[]
    local_x_coordinates = basis.x_coordinates[x_coordinate_range]
    local_y_coordinates = basis.y_coordinates[y_coordinate_range]

    while true
        candidates = _nested_square_lattice_split_candidates(
            bundles,
            current_box,
            x_coordinate_range,
            y_coordinate_range,
            basis.x_coordinates,
            basis.y_coordinates;
            nside = nside,
            min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
            prefer_midpoint_tie_side = prefer_midpoint_tie_side,
        )
        accepted_candidate = _nested_choose_square_lattice_candidate(
            candidates,
            local_x_coordinates,
            local_y_coordinates,
        )
        (accepted_candidate !== nothing || minimum(length.(current_box)) <= nside || !_nested_can_shrink_box(current_box)) && break

        inner_box = _nested_inner_box(current_box)
        shell = _nested_complete_rectangular_shell(
            bundles,
            inner_box...;
            retain_xy = retain_xy,
            retain_xz = retain_xz,
            retain_yz = retain_yz,
            retain_x_edge = retain_x_edge,
            retain_y_edge = retain_y_edge,
            retain_z_edge = retain_z_edge,
            x_fixed = (first(current_box[1]), last(current_box[1])),
            y_fixed = (first(current_box[2]), last(current_box[2])),
            z_fixed = (first(current_box[3]), last(current_box[3])),
        )
        push!(shared_shell_layers, shell)
        push!(shared_shell_dimensions, size(shell.coefficient_matrix, 2))
        current_box = inner_box
    end

    candidates = _nested_square_lattice_split_candidates(
        bundles,
        current_box,
        x_coordinate_range,
        y_coordinate_range,
        basis.x_coordinates,
        basis.y_coordinates;
        nside = nside,
        min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
    )
    accepted_candidate = _nested_choose_square_lattice_candidate(
        candidates,
        local_x_coordinates,
        local_y_coordinates,
    )

    if isnothing(accepted_candidate)
        leaf_sequence = _nested_complete_shell_sequence_for_box(
            bundles,
            current_box;
            nside = nside,
            retain_xy = retain_xy,
            retain_xz = retain_xz,
            retain_yz = retain_yz,
            retain_x_edge = retain_x_edge,
            retain_y_edge = retain_y_edge,
            retain_z_edge = retain_z_edge,
        )
        merged_sequence =
            isempty(shared_shell_layers) ? leaf_sequence :
            _nested_shell_sequence_from_core_block(
                bundles,
                leaf_sequence.support_indices,
                leaf_sequence.coefficient_matrix,
                shared_shell_layers,
            )
        node = _AxisAlignedHomonuclearSquareLatticeNodeGeometry3D(
            String(node_label),
            parent_box,
            current_box,
            basis.lattice_size,
            nside,
            min_in_plane_aspect_ratio,
            x_coordinate_range,
            y_coordinate_range,
            collect(local_x_coordinates),
            collect(local_y_coordinates),
            length(shared_shell_layers),
            shared_shell_dimensions,
            candidates,
            nothing,
            !isempty(candidates),
            _AxisAlignedHomonuclearSquareLatticeNodeGeometry3D[],
            size(merged_sequence.coefficient_matrix, 2),
        )
        return (
            geometry = node,
            sequence = merged_sequence,
            leaf_sequences = _CartesianNestedShellSequence3D[leaf_sequence],
        )
    end

    candidate = candidates[accepted_candidate]
    child_results = [
        _nested_axis_aligned_homonuclear_square_lattice_node(
            basis,
            bundles,
            candidate.child_boxes[index],
            candidate.x_coordinate_ranges[index],
            candidate.y_coordinate_ranges[index],
            string(node_label, "_", index);
            nside = nside,
            min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
            prefer_midpoint_tie_side = prefer_midpoint_tie_side,
            retain_xy = retain_xy,
            retain_xz = retain_xz,
            retain_yz = retain_yz,
            retain_x_edge = retain_x_edge,
            retain_y_edge = retain_y_edge,
            retain_z_edge = retain_z_edge,
        ) for index in eachindex(candidate.child_boxes)
    ]
    child_nodes = [result.geometry for result in child_results]
    child_leaf_sequences = _CartesianNestedShellSequence3D[]
    for result in child_results
        append!(child_leaf_sequences, result.leaf_sequences)
    end
    core_support = reduce(vcat, [result.sequence.support_indices for result in child_results])
    core_coefficients = reduce(hcat, [result.sequence.coefficient_matrix for result in child_results])
    merged_sequence =
        isempty(shared_shell_layers) ? _nested_shell_sequence_from_core_block(
            bundles,
            core_support,
            core_coefficients,
            _AbstractCartesianNestedShellLayer3D[];
            enforce_coverage = false,
        ) :
        _nested_shell_sequence_from_core_block(
            bundles,
            core_support,
            core_coefficients,
            shared_shell_layers,
        )
    node = _AxisAlignedHomonuclearSquareLatticeNodeGeometry3D(
        String(node_label),
        parent_box,
        current_box,
        basis.lattice_size,
        nside,
        min_in_plane_aspect_ratio,
        x_coordinate_range,
        y_coordinate_range,
        collect(local_x_coordinates),
        collect(local_y_coordinates),
        length(shared_shell_layers),
        shared_shell_dimensions,
        candidates,
        accepted_candidate,
        false,
        child_nodes,
        size(merged_sequence.coefficient_matrix, 2),
    )
    return (
        geometry = node,
        sequence = merged_sequence,
        leaf_sequences = child_leaf_sequences,
    )
end

function _nested_axis_aligned_homonuclear_square_lattice_source(
    basis,
    bundles::_CartesianNestedAxisBundles3D;
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
    prefer_midpoint_tie_side::Symbol = :left,
    retain_xy::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_xz::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_yz::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_x_edge::Union{Nothing,Int} = nothing,
    retain_y_edge::Union{Nothing,Int} = nothing,
    retain_z_edge::Union{Nothing,Int} = nothing,
)
    retention = _nested_resolve_complete_shell_retention(
        nside;
        retain_xy = retain_xy,
        retain_xz = retain_xz,
        retain_yz = retain_yz,
        retain_x_edge = retain_x_edge,
        retain_y_edge = retain_y_edge,
        retain_z_edge = retain_z_edge,
    )
    dims = _nested_axis_lengths(bundles)
    parent_box = (1:dims[1], 1:dims[2], 1:dims[3])
    result = _nested_axis_aligned_homonuclear_square_lattice_node(
        basis,
        bundles,
        parent_box,
        1:length(basis.x_coordinates),
        1:length(basis.y_coordinates),
        "root";
        nside = nside,
        min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
        retain_xy = retention.retain_xy,
        retain_xz = retention.retain_xz,
        retain_yz = retention.retain_yz,
        retain_x_edge = retention.retain_x_edge,
        retain_y_edge = retention.retain_y_edge,
        retain_z_edge = retention.retain_z_edge,
    )
    return _CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D(
        basis,
        bundles,
        retention,
        result.geometry,
        result.leaf_sequences,
        result.sequence,
    )
end

# Alg Nested-Face hierarchy step: Refine only the retained direct core block
# inside a trusted nonrecursive shell anchor, using the original parent-space
# functions assigned to that core region rather than re-coarsening already
# renormalized shell functions.
function _nested_hierarchical_core_refinement(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int};
    retain_xy::Tuple{Int,Int} = (2, 2),
    retain_xz::Tuple{Int,Int} = (2, 2),
    retain_yz::Tuple{Int,Int} = (2, 2),
    retain_x_edge::Int = 2,
    retain_y_edge::Int = 2,
    retain_z_edge::Int = 2,
)
    length(x_interval) >= 5 || throw(
        ArgumentError("hierarchical core refinement requires an x interval of length at least 5"),
    )
    length(y_interval) >= 5 || throw(
        ArgumentError("hierarchical core refinement requires a y interval of length at least 5"),
    )
    length(z_interval) >= 5 || throw(
        ArgumentError("hierarchical core refinement requires a z interval of length at least 5"),
    )

    inner_x = (first(x_interval) + 1):(last(x_interval) - 1)
    inner_y = (first(y_interval) + 1):(last(y_interval) - 1)
    inner_z = (first(z_interval) + 1):(last(z_interval) - 1)
    inner_shell = _nested_complete_rectangular_shell(
        pgdg,
        inner_x,
        inner_y,
        inner_z;
        retain_xy = retain_xy,
        retain_xz = retain_xz,
        retain_yz = retain_yz,
        retain_x_edge = retain_x_edge,
        retain_y_edge = retain_y_edge,
        retain_z_edge = retain_z_edge,
        x_fixed = (first(x_interval), last(x_interval)),
        y_fixed = (first(y_interval), last(y_interval)),
        z_fixed = (first(z_interval), last(z_interval)),
    )
    return _nested_shell_sequence(
        pgdg,
        inner_x,
        inner_y,
        inner_z,
        [inner_shell],
    )
end

function _nested_hierarchical_core_refinement(
    bundle::_MappedOrdinaryGausslet1DBundle,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int};
    kwargs...,
)
    return _nested_hierarchical_core_refinement(
        bundle.pgdg_intermediate,
        x_interval,
        y_interval,
        z_interval;
        kwargs...,
    )
end

function _nested_shell_sequence_with_hierarchical_core_refinement(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D};
    kwargs...,
)
    refined_core = _nested_hierarchical_core_refinement(
        pgdg,
        x_interval,
        y_interval,
        z_interval;
        kwargs...,
    )
    sequence = _nested_shell_sequence_from_core_block(
        pgdg,
        refined_core.support_indices,
        refined_core.coefficient_matrix,
        shell_layers,
    )
    return (
        refined_core = refined_core,
        sequence = sequence,
    )
end

function _nested_shell_sequence_with_hierarchical_core_refinement(
    bundle::_MappedOrdinaryGausslet1DBundle,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D};
    kwargs...,
)
    return _nested_shell_sequence_with_hierarchical_core_refinement(
        bundle.pgdg_intermediate,
        x_interval,
        y_interval,
        z_interval,
        shell_layers;
        kwargs...,
    )
end

function _nested_shell_sequence(
    bundle::_MappedOrdinaryGausslet1DBundle,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D};
    enforce_coverage::Bool = true,
)
    return _nested_shell_sequence(
        bundle.pgdg_intermediate,
        x_interval,
        y_interval,
        z_interval,
        shell_layers;
        enforce_coverage = enforce_coverage,
    )
end

# Alg Nested-Face policy step: Apply the first legacy-style fixed-nside
# replacement rule to the direct core block while leaving the shell faces and
# downstream fixed-block consumer unchanged.
# See docs/src/algorithms/cartesian_nested_face_construction.md.
function _nested_nside_shell_sequence(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D};
    nside::Int = 5,
)
    shrunk_x, shrunk_y, shrunk_z = _nested_shrunk_box(
        x_interval,
        y_interval,
        z_interval,
        length(shell_layers);
        nside = nside,
    )
    if length(shrunk_x) != nside || length(shrunk_y) != nside || length(shrunk_z) != nside
        throw(ArgumentError("nested fixed-nside policy requires enough shell layers to reduce the raw interior to an nside × nside × nside box before the final contracted core step"))
    end
    retained = max(1, nside - 2)
    core_data = _nested_contracted_core_coefficients(
        pgdg,
        shrunk_x,
        shrunk_y,
        shrunk_z;
        retain_x = retained,
        retain_y = retained,
        retain_z = retained,
    )
    return _nested_shell_sequence_from_core_block(
        pgdg,
        core_data.support_indices,
        core_data.coefficient_matrix,
        shell_layers;
        enforce_coverage = false,
    )
end

function _nested_nside_shell_sequence(
    bundle::_MappedOrdinaryGausslet1DBundle,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D};
    nside::Int = 5,
)
    return _nested_nside_shell_sequence(
        bundle.pgdg_intermediate,
        x_interval,
        y_interval,
        z_interval,
        shell_layers;
        nside = nside,
    )
end

function _nested_fixed_block(
    shell::_CartesianNestedShellPlusCore3D,
    parent_basis::MappedUniformBasis,
)
    packet = shell.packet
    fixed_centers = hcat(
        diag(packet.position_x),
        diag(packet.position_y),
        diag(packet.position_z),
    )
    return _NestedFixedBlock3D(
        parent_basis,
        shell,
        shell.coefficient_matrix,
        shell.support_indices,
        packet.overlap,
        packet.kinetic,
        packet.position_x,
        packet.position_y,
        packet.position_z,
        packet.x2_x,
        packet.x2_y,
        packet.x2_z,
        packet.weights,
        packet.gaussian_terms,
        packet.pair_terms,
        Matrix{Float64}(fixed_centers),
    )
end

function _nested_fixed_block(
    shell::_CartesianNestedShell3D,
    bundle::_MappedOrdinaryGausslet1DBundle,
)
    return _nested_fixed_block(shell, bundle.basis)
end

function _nested_fixed_block(
    shell::_CartesianNestedCompleteShell3D,
    parent_basis::MappedUniformBasis,
)
    packet = shell.packet
    fixed_centers = hcat(
        diag(packet.position_x),
        diag(packet.position_y),
        diag(packet.position_z),
    )
    return _NestedFixedBlock3D(
        parent_basis,
        shell,
        shell.coefficient_matrix,
        shell.support_indices,
        packet.overlap,
        packet.kinetic,
        packet.position_x,
        packet.position_y,
        packet.position_z,
        packet.x2_x,
        packet.x2_y,
        packet.x2_z,
        packet.weights,
        packet.gaussian_terms,
        packet.pair_terms,
        Matrix{Float64}(fixed_centers),
    )
end

function _nested_fixed_block(
    shell::_CartesianNestedCompleteShell3D,
    bundle::_MappedOrdinaryGausslet1DBundle,
)
    return _nested_fixed_block(shell, bundle.basis)
end

function _nested_fixed_block(
    shell::_CartesianNestedShellPlusCore3D,
    bundle::_MappedOrdinaryGausslet1DBundle,
)
    return _nested_fixed_block(shell, bundle.basis)
end

function _nested_fixed_block(
    shell::_CartesianNestedShellSequence3D,
    parent_basis::MappedUniformBasis,
)
    packet = shell.packet
    fixed_centers = hcat(
        diag(packet.position_x),
        diag(packet.position_y),
        diag(packet.position_z),
    )
    return _NestedFixedBlock3D(
        parent_basis,
        shell,
        shell.coefficient_matrix,
        shell.support_indices,
        packet.overlap,
        packet.kinetic,
        packet.position_x,
        packet.position_y,
        packet.position_z,
        packet.x2_x,
        packet.x2_y,
        packet.x2_z,
        packet.weights,
        packet.gaussian_terms,
        packet.pair_terms,
        Matrix{Float64}(fixed_centers),
    )
end

function _nested_fixed_block(
    shell::_CartesianNestedShellSequence3D,
    parent_basis,
)
    packet = shell.packet
    fixed_centers = hcat(
        diag(packet.position_x),
        diag(packet.position_y),
        diag(packet.position_z),
    )
    return _NestedFixedBlock3D(
        parent_basis,
        shell,
        shell.coefficient_matrix,
        shell.support_indices,
        packet.overlap,
        packet.kinetic,
        packet.position_x,
        packet.position_y,
        packet.position_z,
        packet.x2_x,
        packet.x2_y,
        packet.x2_z,
        packet.weights,
        packet.gaussian_terms,
        packet.pair_terms,
        Matrix{Float64}(fixed_centers),
    )
end

function _nested_fixed_block(
    shell::_CartesianNestedShellSequence3D,
    bundle::_MappedOrdinaryGausslet1DBundle,
)
    return _nested_fixed_block(shell, bundle.basis)
end

function _nested_fixed_block(source::_CartesianNestedBondAlignedDiatomicSource3D)
    return _nested_fixed_block(source.sequence, source.basis)
end

function _nested_fixed_block(source::_CartesianNestedBondAlignedHomonuclearChainSource3D)
    return _nested_fixed_block(source.sequence, source.basis)
end

function _nested_fixed_block(source::_CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D)
    return _nested_fixed_block(source.sequence, source.basis)
end
