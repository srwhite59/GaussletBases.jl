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
    geometry::_BondAlignedDiatomicSplitGeometry3D
    shared_shell_layers::Vector{_CartesianNestedCompleteShell3D}
    child_sequences::Vector{_CartesianNestedShellSequence3D}
    child_column_ranges::Vector{UnitRange{Int}}
    midpoint_slab_column_range::Union{Nothing,UnitRange{Int}}
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
        ", midpoint_slab=",
        !isnothing(source.midpoint_slab_column_range),
        ", did_split=",
        source.geometry.did_split,
        ")",
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

# Alg Nested-Face step 3: Build a local 1D doside contraction on one interval.
# See docs/src/algorithms/cartesian_nested_face_construction.md.
function _nested_doside_1d(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    interval::UnitRange{Int},
    retained_count::Int,
)
    interval_data = _nested_interval_data(pgdg, interval)
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
    nsupport = length(support_states)
    matrix = zeros(Float64, nsupport, nsupport)
    for row in 1:nsupport
        ix, iy, iz = support_states[row]
        for col in 1:nsupport
            jx, jy, jz = support_states[col]
            matrix[row, col] =
                Float64(operator_x[ix, jx]) *
                Float64(operator_y[iy, jy]) *
                Float64(operator_z[iz, jz])
        end
    end
    return matrix
end

function _nested_sum_of_support_products(
    support_states::AbstractVector{<:NTuple{3,Int}},
    terms,
)
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
    support_coefficients::AbstractMatrix{<:Real},
)
    nterms = size(pgdg.pair_factor_terms, 1)
    nfixed = size(support_coefficients, 2)
    parent_weight_outer = pgdg.weights * transpose(pgdg.weights)
    support_weights = _nested_support_weights(support_states, pgdg.weights)
    fixed_weights = vec(transpose(support_coefficients) * support_weights)
    all(isfinite, fixed_weights) || throw(
        ArgumentError("nested fixed-block IDA transfer requires finite contracted integral weights"),
    )
    minimum(fixed_weights) > 0.0 || throw(
        ArgumentError("nested fixed-block IDA transfer requires positive contracted integral weights"),
    )
    fixed_weight_outer = fixed_weights * transpose(fixed_weights)
    pair_terms = zeros(Float64, nterms, nfixed, nfixed)

    for term in 1:nterms
        raw_1d = @view(pgdg.pair_factor_terms[term, :, :]) .* parent_weight_outer
        raw_support = _nested_support_product_matrix(
            support_states,
            raw_1d,
            raw_1d,
            raw_1d,
        )
        raw_contracted = transpose(support_coefficients) * raw_support * support_coefficients
        matrix = Matrix{Float64}(raw_contracted ./ fixed_weight_outer)
        @views pair_terms[term, :, :] .= 0.5 .* (matrix .+ transpose(matrix))
    end

    return (
        weights = fixed_weights,
        pair_terms = pair_terms,
    )
end

function _nested_weight_aware_pair_terms(
    bundles::_CartesianNestedAxisBundles3D,
    support_states::AbstractVector{<:NTuple{3,Int}},
    support_coefficients::AbstractMatrix{<:Real},
)
    pgdg_x = _nested_axis_pgdg(bundles, :x)
    pgdg_y = _nested_axis_pgdg(bundles, :y)
    pgdg_z = _nested_axis_pgdg(bundles, :z)
    nterms = size(pgdg_x.pair_factor_terms, 1)
    nterms == size(pgdg_y.pair_factor_terms, 1) == size(pgdg_z.pair_factor_terms, 1) || throw(
        ArgumentError("mixed-axis nested IDA transfer requires the same Gaussian expansion term count on all axes"),
    )
    nfixed = size(support_coefficients, 2)
    parent_weight_outer_x = pgdg_x.weights * transpose(pgdg_x.weights)
    parent_weight_outer_y = pgdg_y.weights * transpose(pgdg_y.weights)
    parent_weight_outer_z = pgdg_z.weights * transpose(pgdg_z.weights)
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
    fixed_weight_outer = fixed_weights * transpose(fixed_weights)
    pair_terms = zeros(Float64, nterms, nfixed, nfixed)

    for term in 1:nterms
        raw_x = @view(pgdg_x.pair_factor_terms[term, :, :]) .* parent_weight_outer_x
        raw_y = @view(pgdg_y.pair_factor_terms[term, :, :]) .* parent_weight_outer_y
        raw_z = @view(pgdg_z.pair_factor_terms[term, :, :]) .* parent_weight_outer_z
        raw_support = _nested_support_product_matrix(
            support_states,
            raw_x,
            raw_y,
            raw_z,
        )
        raw_contracted = transpose(support_coefficients) * raw_support * support_coefficients
        matrix = Matrix{Float64}(raw_contracted ./ fixed_weight_outer)
        @views pair_terms[term, :, :] .= 0.5 .* (matrix .+ transpose(matrix))
    end

    return (
        weights = fixed_weights,
        pair_terms = pair_terms,
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
)
    support_states = [_cartesian_unflat_index(index, size(pgdg.overlap, 1)) for index in support_indices]
    support_coefficients = Matrix{Float64}(coefficient_matrix[support_indices, :])
    overlap_support = _nested_support_product_matrix(
        support_states,
        pgdg.overlap,
        pgdg.overlap,
        pgdg.overlap,
    )
    kinetic_support = _nested_sum_of_support_products(
        support_states,
        (
            (pgdg.kinetic, pgdg.overlap, pgdg.overlap),
            (pgdg.overlap, pgdg.kinetic, pgdg.overlap),
            (pgdg.overlap, pgdg.overlap, pgdg.kinetic),
        ),
    )
    position_x_support = _nested_support_product_matrix(support_states, pgdg.position, pgdg.overlap, pgdg.overlap)
    position_y_support = _nested_support_product_matrix(support_states, pgdg.overlap, pgdg.position, pgdg.overlap)
    position_z_support = _nested_support_product_matrix(support_states, pgdg.overlap, pgdg.overlap, pgdg.position)
    x2_x_support = _nested_support_product_matrix(support_states, pgdg.x2, pgdg.overlap, pgdg.overlap)
    x2_y_support = _nested_support_product_matrix(support_states, pgdg.overlap, pgdg.x2, pgdg.overlap)
    x2_z_support = _nested_support_product_matrix(support_states, pgdg.overlap, pgdg.overlap, pgdg.x2)

    nshell = size(coefficient_matrix, 2)
    nterms = size(pgdg.gaussian_factor_terms, 1)
    pair_data = _nested_weight_aware_pair_terms(pgdg, support_states, support_coefficients)
    gaussian_terms = zeros(Float64, nterms, nshell, nshell)
    for term in 1:nterms
        factor_support = _nested_support_product_matrix(
            support_states,
            @view(pgdg.gaussian_factor_terms[term, :, :]),
            @view(pgdg.gaussian_factor_terms[term, :, :]),
            @view(pgdg.gaussian_factor_terms[term, :, :]),
        )
        @views gaussian_terms[term, :, :] .= transpose(support_coefficients) * factor_support * support_coefficients
    end

    return (
        packet = _CartesianNestedShellPacket3D(
            transpose(support_coefficients) * overlap_support * support_coefficients,
            transpose(support_coefficients) * kinetic_support * support_coefficients,
            transpose(support_coefficients) * position_x_support * support_coefficients,
            transpose(support_coefficients) * position_y_support * support_coefficients,
            transpose(support_coefficients) * position_z_support * support_coefficients,
            transpose(support_coefficients) * x2_x_support * support_coefficients,
            transpose(support_coefficients) * x2_y_support * support_coefficients,
            transpose(support_coefficients) * x2_z_support * support_coefficients,
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
)
    dims = _nested_axis_lengths(bundles)
    pgdg_x = _nested_axis_pgdg(bundles, :x)
    pgdg_y = _nested_axis_pgdg(bundles, :y)
    pgdg_z = _nested_axis_pgdg(bundles, :z)
    support_states = [_cartesian_unflat_index(index, dims) for index in support_indices]
    support_coefficients = Matrix{Float64}(coefficient_matrix[support_indices, :])
    overlap_support = _nested_support_product_matrix(
        support_states,
        pgdg_x.overlap,
        pgdg_y.overlap,
        pgdg_z.overlap,
    )
    kinetic_support = _nested_sum_of_support_products(
        support_states,
        (
            (pgdg_x.kinetic, pgdg_y.overlap, pgdg_z.overlap),
            (pgdg_x.overlap, pgdg_y.kinetic, pgdg_z.overlap),
            (pgdg_x.overlap, pgdg_y.overlap, pgdg_z.kinetic),
        ),
    )
    position_x_support = _nested_support_product_matrix(
        support_states,
        pgdg_x.position,
        pgdg_y.overlap,
        pgdg_z.overlap,
    )
    position_y_support = _nested_support_product_matrix(
        support_states,
        pgdg_x.overlap,
        pgdg_y.position,
        pgdg_z.overlap,
    )
    position_z_support = _nested_support_product_matrix(
        support_states,
        pgdg_x.overlap,
        pgdg_y.overlap,
        pgdg_z.position,
    )
    x2_x_support = _nested_support_product_matrix(
        support_states,
        pgdg_x.x2,
        pgdg_y.overlap,
        pgdg_z.overlap,
    )
    x2_y_support = _nested_support_product_matrix(
        support_states,
        pgdg_x.overlap,
        pgdg_y.x2,
        pgdg_z.overlap,
    )
    x2_z_support = _nested_support_product_matrix(
        support_states,
        pgdg_x.overlap,
        pgdg_y.overlap,
        pgdg_z.x2,
    )

    nshell = size(coefficient_matrix, 2)
    nterms = size(pgdg_x.gaussian_factor_terms, 1)
    nterms == size(pgdg_y.gaussian_factor_terms, 1) == size(pgdg_z.gaussian_factor_terms, 1) || throw(
        ArgumentError("mixed-axis nested shell packets require the same Gaussian expansion term count on all axes"),
    )
    pair_data = _nested_weight_aware_pair_terms(bundles, support_states, support_coefficients)
    gaussian_terms = zeros(Float64, nterms, nshell, nshell)
    for term in 1:nterms
        factor_support = _nested_support_product_matrix(
            support_states,
            @view(pgdg_x.gaussian_factor_terms[term, :, :]),
            @view(pgdg_y.gaussian_factor_terms[term, :, :]),
            @view(pgdg_z.gaussian_factor_terms[term, :, :]),
        )
        @views gaussian_terms[term, :, :] .= transpose(support_coefficients) * factor_support * support_coefficients
    end

    return (
        packet = _CartesianNestedShellPacket3D(
            transpose(support_coefficients) * overlap_support * support_coefficients,
            transpose(support_coefficients) * kinetic_support * support_coefficients,
            transpose(support_coefficients) * position_x_support * support_coefficients,
            transpose(support_coefficients) * position_y_support * support_coefficients,
            transpose(support_coefficients) * position_z_support * support_coefficients,
            transpose(support_coefficients) * x2_x_support * support_coefficients,
            transpose(support_coefficients) * x2_y_support * support_coefficients,
            transpose(support_coefficients) * x2_z_support * support_coefficients,
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
    shell_data = _nested_shell_packet(pgdg, coefficient_matrix, support_indices)
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
)
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
    shell_data = _nested_shell_packet(pgdg, coefficient_matrix, support_indices)

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
)
    n1d = size(pgdg.overlap, 1)
    core_indices = _nested_box_support_indices(x_interval, y_interval, z_interval, n1d)
    isempty(intersect(core_indices, shell.support_indices)) || throw(
        ArgumentError("nested shell-plus-core construction requires the direct core block to stay disjoint from the shell-face supports"),
    )
    core_coefficients = _nested_direct_core_coefficients(core_indices, n1d^3)
    coefficient_matrix = hcat(core_coefficients, shell.coefficient_matrix)
    support_indices = sort(vcat(core_indices, shell.support_indices))
    shell_data = _nested_shell_packet(pgdg, coefficient_matrix, support_indices)
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
    )
end

function _nested_shell_sequence_from_core_block(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    core_indices::AbstractVector{Int},
    core_coefficients::AbstractMatrix{<:Real},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D};
    enforce_coverage::Bool = true,
)
    n1d = size(pgdg.overlap, 1)
    support_indices = _nested_sequence_support_indices(core_indices, shell_layers)
    working_box =
        enforce_coverage ?
        _nested_assert_sequence_coverage(core_indices, shell_layers, support_indices, n1d) :
        _nested_sequence_working_box(core_indices, shell_layers, n1d)
    coefficient_blocks = Matrix{Float64}[Matrix{Float64}(core_coefficients)]
    append!(coefficient_blocks, [shell.coefficient_matrix for shell in shell_layers])
    coefficient_matrix = hcat(coefficient_blocks...)
    shell_data = _nested_shell_packet(pgdg, coefficient_matrix, support_indices)

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
    shell_data = _nested_shell_packet(bundles, coefficient_matrix, support_indices)
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
)
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
    shell_data = _nested_shell_packet(bundles, coefficient_matrix, support_indices)

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
    )
end

function _nested_shell_sequence_from_core_block(
    bundles::_CartesianNestedAxisBundles3D,
    core_indices::AbstractVector{Int},
    core_coefficients::AbstractMatrix{<:Real},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D};
    enforce_coverage::Bool = true,
)
    dims = _nested_axis_lengths(bundles)
    support_indices = _nested_sequence_support_indices(core_indices, shell_layers)
    working_box =
        enforce_coverage ?
        _nested_assert_sequence_coverage(core_indices, shell_layers, support_indices, dims) :
        _nested_sequence_working_box(core_indices, shell_layers, dims)
    coefficient_blocks = Matrix{Float64}[Matrix{Float64}(core_coefficients)]
    append!(coefficient_blocks, [shell.coefficient_matrix for shell in shell_layers])
    coefficient_matrix = hcat(coefficient_blocks...)
    shell_data = _nested_shell_packet(bundles, coefficient_matrix, support_indices)

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
    retain_xy::Tuple{Int,Int} = (4, 3),
    retain_xz::Tuple{Int,Int} = (4, 3),
    retain_yz::Tuple{Int,Int} = (4, 3),
    retain_x_edge::Int = 3,
    retain_y_edge::Int = 3,
    retain_z_edge::Int = 3,
)
    nside >= 1 || throw(ArgumentError("diatomic child shell sequence requires nside >= 1"))
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
                retain_xy = retain_xy,
                retain_xz = retain_xz,
                retain_yz = retain_yz,
                retain_x_edge = retain_x_edge,
                retain_y_edge = retain_y_edge,
                retain_z_edge = retain_z_edge,
                x_fixed = (first(current_box[1]), last(current_box[1])),
                y_fixed = (first(current_box[2]), last(current_box[2])),
                z_fixed = (first(current_box[3]), last(current_box[3])),
            ),
        )
        current_box = inner_box
    end
    return _nested_shell_sequence(
        bundles,
        current_box...,
        shell_layers,
    )
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

function _nested_diatomic_split_index(
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

function _nested_symmetric_interval_retained_count(
    bundles::_CartesianNestedAxisBundles3D,
    axis::Symbol,
    interval::UnitRange{Int},
    provisional_retained_count::Int,
)
    provisional_retained_count >= 1 || throw(
        ArgumentError("shared-shell retained local count must be at least 1"),
    )
    centers = Float64[_nested_axis_pgdg(bundles, axis).centers[index] for index in interval]
    symmetric_about_zero, _ = _nested_is_symmetric_about_zero(centers)
    if symmetric_about_zero && iseven(provisional_retained_count) && provisional_retained_count > 1
        return provisional_retained_count - 1
    end
    return provisional_retained_count
end

# Alg Nested-Diatomic step 7: In the current homonuclear shared-shell path,
# keep tangential retained local counts odd on intervals symmetric about zero.
# This preserves the near-zero localized center without freezing a general
# adaptive retain-count rule.
# See docs/src/algorithms/cartesian_nested_diatomic_box_policy.md.
function _nested_homonuclear_shared_shell_face_retains(
    bundles::_CartesianNestedAxisBundles3D,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int};
    retain_xy::Tuple{Int,Int},
    retain_xz::Tuple{Int,Int},
    retain_yz::Tuple{Int,Int},
)
    return (
        retain_xy = (
            _nested_symmetric_interval_retained_count(bundles, :x, x_interval, retain_xy[1]),
            _nested_symmetric_interval_retained_count(bundles, :y, y_interval, retain_xy[2]),
        ),
        retain_xz = (
            _nested_symmetric_interval_retained_count(bundles, :x, x_interval, retain_xz[1]),
            _nested_symmetric_interval_retained_count(bundles, :z, z_interval, retain_xz[2]),
        ),
        retain_yz = (
            _nested_symmetric_interval_retained_count(bundles, :y, y_interval, retain_yz[1]),
            _nested_symmetric_interval_retained_count(bundles, :z, z_interval, retain_yz[2]),
        ),
    )
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
)
    axis = bond_axis == :x ? 1 : bond_axis == :y ? 2 : bond_axis == :z ? 3 : 0
    axis != 0 || throw(ArgumentError("diatomic split geometry requires bond_axis = :x, :y, or :z"))
    parallel_interval = working_box[axis]
    parallel_centers = _nested_axis_pgdg(bundles, bond_axis).centers
    split_index = _nested_diatomic_split_index(parallel_centers, parallel_interval, midpoint)
    use_slab = use_midpoint_slab && isodd(length(parallel_interval))
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
    shared_shell_retain_xy::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_xz::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_yz::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_xy::Tuple{Int,Int} = (4, 3),
    retain_xz::Tuple{Int,Int} = (4, 3),
    retain_yz::Tuple{Int,Int} = (4, 3),
    retain_x_edge::Int = 3,
    retain_y_edge::Int = 3,
    retain_z_edge::Int = 3,
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
    )

    while true
        parallel_length = length(current_box[bond_axis == :x ? 1 : bond_axis == :y ? 2 : 3])
        if parallel_length <= 2 * nside || minimum(length.(current_box)) <= nside || !_nested_can_shrink_box(current_box)
            break
        end
        inner_box = _nested_inner_box(current_box)
        shared_default_retains = _nested_homonuclear_shared_shell_face_retains(
            bundles,
            inner_box...;
            retain_xy = retain_xy,
            retain_xz = retain_xz,
            retain_yz = retain_yz,
        )
        shared_retain_xy =
            isnothing(shared_shell_retain_xy) ? shared_default_retains.retain_xy : shared_shell_retain_xy
        shared_retain_xz =
            isnothing(shared_shell_retain_xz) ? shared_default_retains.retain_xz : shared_shell_retain_xz
        shared_retain_yz =
            isnothing(shared_shell_retain_yz) ? shared_default_retains.retain_yz : shared_shell_retain_yz
        push!(
            shared_shell_layers,
            _nested_complete_rectangular_shell(
                bundles,
                inner_box...;
                retain_xy = shared_retain_xy,
                retain_xz = shared_retain_xz,
                retain_yz = shared_retain_yz,
                retain_x_edge = retain_x_edge,
                retain_y_edge = retain_y_edge,
                retain_z_edge = retain_z_edge,
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
                    retain_xy = retain_xy,
                    retain_xz = retain_xz,
                    retain_yz = retain_yz,
                    retain_x_edge = retain_x_edge,
                    retain_y_edge = retain_y_edge,
                    retain_z_edge = retain_z_edge,
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
            retain_xy = retain_xy,
            retain_xz = retain_xz,
            retain_yz = retain_yz,
            retain_x_edge = retain_x_edge,
            retain_y_edge = retain_y_edge,
            retain_z_edge = retain_z_edge,
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
        geometry,
        shared_shell_layers,
        child_sequences,
        child_column_ranges,
        midpoint_slab_column_range,
        merged_sequence,
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
