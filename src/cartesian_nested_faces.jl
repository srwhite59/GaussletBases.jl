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

function _cartesian_flat_index(ix::Int, iy::Int, iz::Int, n1d::Int)
    return (ix - 1) * n1d * n1d + (iy - 1) * n1d + iz
end

function _cartesian_unflat_index(index::Int, n1d::Int)
    shifted = index - 1
    plane = n1d * n1d
    ix = shifted ÷ plane + 1
    remainder = shifted % plane
    iy = remainder ÷ n1d + 1
    iz = remainder % n1d + 1
    return (ix, iy, iz)
end

function _nested_xy_face_support_indices(
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_index::Int,
    n1d::Int,
)
    support = Int[]
    for ix in x_interval, iy in y_interval
        push!(support, _cartesian_flat_index(ix, iy, z_index, n1d))
    end
    return support
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
    n1d::Int,
)
    support = Int[]
    if face_kind == :xy
        for ix in interval_first, iy in interval_second
            push!(support, _cartesian_flat_index(ix, iy, fixed_index, n1d))
        end
    elseif face_kind == :xz
        for ix in interval_first, iz in interval_second
            push!(support, _cartesian_flat_index(ix, fixed_index, iz, n1d))
        end
    elseif face_kind == :yz
        for iy in interval_first, iz in interval_second
            push!(support, _cartesian_flat_index(fixed_index, iy, iz, n1d))
        end
    else
        throw(ArgumentError("nested face support requires face kind :xy, :xz, or :yz"))
    end
    return support
end

function _nested_edge_support_indices(
    free_axis::Symbol,
    free_interval::UnitRange{Int},
    fixed_indices::NTuple{2,Int},
    n1d::Int,
)
    support = Int[]
    if free_axis == :x
        iy, iz = fixed_indices
        for ix in free_interval
            push!(support, _cartesian_flat_index(ix, iy, iz, n1d))
        end
    elseif free_axis == :y
        ix, iz = fixed_indices
        for iy in free_interval
            push!(support, _cartesian_flat_index(ix, iy, iz, n1d))
        end
    elseif free_axis == :z
        ix, iy = fixed_indices
        for iz in free_interval
            push!(support, _cartesian_flat_index(ix, iy, iz, n1d))
        end
    else
        throw(ArgumentError("nested edge support requires free axis :x, :y, or :z"))
    end
    return support
end

function _nested_edge_product(
    free_axis::Symbol,
    fixed_sides::NTuple{2,Symbol},
    side::_CartesianNestedDoSide1D,
    fixed_indices::NTuple{2,Int},
    n1d::Int,
)
    1 <= minimum(fixed_indices) <= maximum(fixed_indices) <= n1d || throw(
        ArgumentError("nested edge requires fixed indices inside the finalized Cartesian line"),
    )
    coefficient_matrix = zeros(Float64, n1d^3, size(side.coefficient_matrix, 2))
    for col in 1:size(side.coefficient_matrix, 2)
        for free_index in side.interval
            value = side.coefficient_matrix[free_index, col]
            iszero(value) && continue
            flat =
                free_axis == :x ? _cartesian_flat_index(free_index, fixed_indices[1], fixed_indices[2], n1d) :
                free_axis == :y ? _cartesian_flat_index(fixed_indices[1], free_index, fixed_indices[2], n1d) :
                _cartesian_flat_index(fixed_indices[1], fixed_indices[2], free_index, n1d)
            coefficient_matrix[flat, col] = value
        end
    end
    fixed_axes =
        free_axis == :x ? (:y, :z) :
        free_axis == :y ? (:x, :z) :
        (:x, :y)
    support_indices = _nested_edge_support_indices(free_axis, side.interval, fixed_indices, n1d)
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

function _nested_corner_piece(
    fixed_sides::NTuple{3,Symbol},
    fixed_indices::NTuple{3,Int},
    n1d::Int,
)
    flat = _cartesian_flat_index(fixed_indices[1], fixed_indices[2], fixed_indices[3], n1d)
    coefficient_matrix = zeros(Float64, n1d^3, 1)
    coefficient_matrix[flat, 1] = 1.0
    return _CartesianNestedCorner3D(
        fixed_sides,
        fixed_indices,
        coefficient_matrix,
        [flat],
    )
end

function _nested_box_support_indices(
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int},
    n1d::Int,
)
    support = Int[]
    for ix in x_interval, iy in y_interval, iz in z_interval
        push!(support, _cartesian_flat_index(ix, iy, iz, n1d))
    end
    return support
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
    n1d::Int,
)
    nparent = n1d^3
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
                    flat = _cartesian_flat_index(ix, iy, iz, n1d)
                    coefficients[flat, column] = vx * vy * vz
                end
            end
        end
    end
    return coefficients
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
    n1d::Int,
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
        update_bounds(_cartesian_unflat_index(index, n1d))
    end
    for shell in shell_layers, state in shell.support_states
        update_bounds(state)
    end

    xmin <= xmax || throw(ArgumentError("nested shell-sequence construction requires at least one retained parent row"))
    return (xmin:xmax, ymin:ymax, zmin:zmax)
end

function _nested_assert_sequence_coverage(
    core_indices::AbstractVector{Int},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D},
    support_indices::AbstractVector{Int},
    n1d::Int,
)
    working_box = _nested_sequence_working_box(core_indices, shell_layers, n1d)
    target_indices = _nested_box_support_indices(working_box..., n1d)
    if target_indices != support_indices
        support_set = Set(support_indices)
        target_set = Set(target_indices)
        missing = length(setdiff(target_set, support_set))
        extra = length(setdiff(support_set, target_set))
        throw(ArgumentError("nested shell-sequence construction requires full coverage of the inferred working box $(working_box): missing $missing parent rows and extra $extra rows"))
    end
    return working_box
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
    n1d::Int,
)
    1 <= fixed_index <= n1d || throw(ArgumentError("nested face requires a fixed index inside the finalized Cartesian line"))
    (fixed_side == :low || fixed_side == :high) || throw(ArgumentError("nested face fixed_side must be :low or :high"))
    _, fixed_axis = _nested_face_axes(face_kind)
    nfirst = size(side_first.coefficient_matrix, 2)
    nsecond = size(side_second.coefficient_matrix, 2)
    coefficients = zeros(Float64, n1d^3, nfirst * nsecond)
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
                    face_kind == :xy ? _cartesian_flat_index(index_first, index_second, fixed_index, n1d) :
                    face_kind == :xz ? _cartesian_flat_index(index_first, fixed_index, index_second, n1d) :
                    _cartesian_flat_index(fixed_index, index_first, index_second, n1d)
                coefficients[flat_index, column] = value_first * value_second
            end
        end
    end
    support_indices = _nested_face_support_indices(
        face_kind,
        side_first.interval,
        side_second.interval,
        fixed_index,
        n1d,
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
    weights = zeros(Float64, length(support_states))
    for (index, state) in enumerate(support_states)
        ix, iy, iz = state
        weights[index] =
            Float64(weights_1d[ix]) *
            Float64(weights_1d[iy]) *
            Float64(weights_1d[iz])
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
    bundle::_MappedOrdinaryGausslet1DBundle,
)
    return _nested_fixed_block(shell, bundle.basis)
end
