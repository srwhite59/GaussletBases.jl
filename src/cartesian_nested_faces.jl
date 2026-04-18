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
const _CartesianCoefficientMap =
    Union{Matrix{Float64},SparseArrays.SparseMatrixCSC{Float64,Int}}

function _cartesian_coefficient_map_storage(
    value::AbstractMatrix{<:Real},
)::_CartesianCoefficientMap
    if value isa SparseArrays.SparseMatrixCSC{Float64,Int}
        return copy(value)
    elseif SparseArrays.issparse(value)
        rows, cols, values = SparseArrays.findnz(value)
        return SparseArrays.sparse(
            Vector{Int}(rows),
            Vector{Int}(cols),
            Float64.(values),
            size(value, 1),
            size(value, 2),
        )
    end
    return Matrix{Float64}(value)
end

function _nested_sparse_coefficient_map(
    row_indices::AbstractVector{<:Integer},
    col_indices::AbstractVector{<:Integer},
    values::AbstractVector{<:Real},
    nrows::Integer,
    ncols::Integer,
)::_CartesianCoefficientMap
    isempty(values) && return SparseArrays.spzeros(Float64, Int(nrows), Int(ncols))
    return SparseArrays.sparse(
        Int.(row_indices),
        Int.(col_indices),
        Float64.(values),
        Int(nrows),
        Int(ncols),
    )
end

function _nested_hcat_sparse_csc_maps(
    blocks::AbstractVector,
)::_CartesianCoefficientMap
    isempty(blocks) && return SparseArrays.spzeros(Float64, 0, 0)
    nrows = size(first(blocks), 1)
    all(size(block, 1) == nrows for block in blocks) || throw(
        ArgumentError("nested sparse coefficient-map concatenation requires a common parent row count"),
    )
    total_columns = sum(size(block, 2) for block in blocks)
    total_nnz = sum(SparseArrays.nnz(block) for block in blocks)
    colptr = Vector{Int}(undef, total_columns + 1)
    rowvals = Vector{Int}(undef, total_nnz)
    nzvals = Vector{Float64}(undef, total_nnz)
    colptr[1] = 1
    column_offset = 0
    nz_offset = 0
    for block in blocks
        block = block::SparseArrays.SparseMatrixCSC{Float64,Int}
        block_nnz = SparseArrays.nnz(block)
        if block_nnz > 0
            copyto!(rowvals, nz_offset + 1, SparseArrays.rowvals(block), 1, block_nnz)
            copyto!(nzvals, nz_offset + 1, SparseArrays.nonzeros(block), 1, block_nnz)
        end
        block_colptr = SparseArrays.getcolptr(block)
        for column in 1:size(block, 2)
            colptr[column_offset + column + 1] = nz_offset + block_colptr[column + 1]
        end
        column_offset += size(block, 2)
        nz_offset += block_nnz
    end
    return SparseArrays.SparseMatrixCSC{Float64,Int}(nrows, total_columns, colptr, rowvals, nzvals)
end

function _nested_hcat_coefficient_maps(
    blocks::AbstractVector{<:AbstractMatrix{<:Real}},
)::_CartesianCoefficientMap
    isempty(blocks) && return zeros(Float64, 0, 0)
    nrows = size(first(blocks), 1)
    all(size(block, 1) == nrows for block in blocks) || throw(
        ArgumentError("nested coefficient-map concatenation requires a common parent row count"),
    )
    if all(block -> !SparseArrays.issparse(block), blocks)
        return Matrix{Float64}(hcat(blocks...))
    end
    if all(block -> block isa SparseArrays.SparseMatrixCSC{Float64,Int}, blocks)
        return _nested_hcat_sparse_csc_maps(blocks)
    end

    row_indices = Int[]
    col_indices = Int[]
    values = Float64[]
    total_columns = 0
    for block in blocks
        if SparseArrays.issparse(block)
            rows, cols, nzvals = SparseArrays.findnz(block)
            append!(row_indices, Int.(rows))
            append!(col_indices, Int.(cols) .+ total_columns)
            append!(values, Float64.(nzvals))
        else
            for col in axes(block, 2), row in axes(block, 1)
                value = Float64(block[row, col])
                iszero(value) && continue
                push!(row_indices, row)
                push!(col_indices, col + total_columns)
                push!(values, value)
            end
        end
        total_columns += size(block, 2)
    end
    return _nested_sparse_coefficient_map(row_indices, col_indices, values, nrows, total_columns)
end

function _nested_support_coefficient_slice(
    coefficient_matrix::_CartesianCoefficientMap,
    support_indices::AbstractVector{Int},
)
    return coefficient_matrix[support_indices, :]
end

function _nested_support_coefficient_slice(
    coefficient_matrix::AbstractMatrix{<:Real},
    support_indices::AbstractVector{Int},
)::_CartesianCoefficientMap
    return Matrix{Float64}(coefficient_matrix[support_indices, :])
end

const _NESTED_SPARSE_SUPPORT_CHUNK_ROWS = 96

function _nested_support_reference_chunk_row_count(
    support_coefficients::SparseArrays.SparseMatrixCSC{Float64,Int},
    nsupport::Int,
    _nfixed::Int,
)
    return min(nsupport, _NESTED_SPARSE_SUPPORT_CHUNK_ROWS)
end

function _nested_support_reference_chunk_row_count(
    _support_coefficients::AbstractMatrix{<:Real},
    nsupport::Int,
    _nfixed::Int,
)
    return nsupport
end

function _nested_support_reference_workspaces(
    support_coefficients::SparseArrays.SparseMatrixCSC{Float64,Int},
    nsupport::Int,
    nfixed::Int,
)
    _nested_support_reference_chunk_row_count(support_coefficients, nsupport, nfixed)
    return Matrix{Float64}(undef, 0, 0), Matrix{Float64}(undef, 0, 0)
end

function _nested_support_reference_workspaces(
    support_coefficients::AbstractMatrix{<:Real},
    nsupport::Int,
    nfixed::Int,
)
    _nested_support_reference_chunk_row_count(support_coefficients, nsupport, nfixed)
    return Matrix{Float64}(undef, nsupport, nsupport), Matrix{Float64}(undef, nfixed, nsupport)
end

struct _CartesianNestedDoSide1D
    interval::UnitRange{Int}
    retained_count::Int
    local_overlap::Matrix{Float64}
    local_position::Matrix{Float64}
    local_weights::Vector{Float64}
    local_centers::Vector{Float64}
    local_coefficients::Matrix{Float64}
    coefficient_matrix::_CartesianCoefficientMap
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
    coefficient_matrix::_CartesianCoefficientMap
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
    gaussian_sum::Union{Nothing,Matrix{Float64}}
    pair_sum::Union{Nothing,Matrix{Float64}}
    gaussian_terms::Union{Nothing,Array{Float64,3}}
    pair_terms::Union{Nothing,Array{Float64,3}}
    term_storage::Symbol
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
    coefficient_matrix::_CartesianCoefficientMap
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
    coefficient_matrix::_CartesianCoefficientMap
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
    coefficient_matrix::_CartesianCoefficientMap
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
    coefficient_matrix::_CartesianCoefficientMap
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
    coefficient_matrix::_CartesianCoefficientMap
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
    coefficient_matrix::_CartesianCoefficientMap
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
    coefficient_matrix::_CartesianCoefficientMap
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
    coefficient_matrix::_CartesianCoefficientMap
    support_indices::Vector{Int}
    support_states::Union{Nothing,Vector{NTuple{3,Int}}}
    packet::Union{Nothing,_CartesianNestedShellPacket3D}
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
    coefficient_matrix::_CartesianCoefficientMap
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
    gaussian_sum::Union{Nothing,Matrix{Float64}}
    pair_sum::Union{Nothing,Matrix{Float64}}
    gaussian_terms::Union{Nothing,Array{Float64,3}}
    pair_terms::Union{Nothing,Array{Float64,3}}
    term_storage::Symbol
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

struct TimedNestedFixedBlockBuild{F}
    fixed_block::F
    timings::TimeG.TimingReport
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

function _embed_local_side_coefficients(
    local_coefficients::AbstractMatrix{<:Real},
    interval::UnitRange{Int},
    n1d::Int,
)
    row_indices = Int[]
    col_indices = Int[]
    values = Float64[]
    for col in axes(local_coefficients, 2), (local_row, full_row) in enumerate(interval)
        value = Float64(local_coefficients[local_row, col])
        iszero(value) && continue
        push!(row_indices, full_row)
        push!(col_indices, col)
        push!(values, value)
    end
    return _nested_sparse_coefficient_map(
        row_indices,
        col_indices,
        values,
        n1d,
        size(local_coefficients, 2),
    )
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
    row_indices = Int[]
    col_indices = Int[]
    values = Float64[]
    column = 0
    for ix_side in 1:nx, iy_side in 1:ny
        column += 1
        for (ix_local, ix) in enumerate(side_x.interval)
            xvalue = Float64(side_x.local_coefficients[ix_local, ix_side])
            iszero(xvalue) && continue
            for (iy_local, iy) in enumerate(side_y.interval)
                yvalue = Float64(side_y.local_coefficients[iy_local, iy_side])
                iszero(yvalue) && continue
                push!(row_indices, _cartesian_flat_index(ix, iy, z_index, n1d))
                push!(col_indices, column)
                push!(values, xvalue * yvalue)
            end
        end
    end
    coefficients = _nested_sparse_coefficient_map(row_indices, col_indices, values, n1d^3, nx * ny)
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
    row_indices = Int[]
    col_indices = Int[]
    values = Float64[]
    for col in 1:size(side.coefficient_matrix, 2)
        for (local_index, free_index) in enumerate(side.interval)
            value = Float64(side.local_coefficients[local_index, col])
            iszero(value) && continue
            flat =
                free_axis == :x ? _cartesian_flat_index(free_index, fixed_indices[1], fixed_indices[2], dims) :
                free_axis == :y ? _cartesian_flat_index(fixed_indices[1], free_index, fixed_indices[2], dims) :
                _cartesian_flat_index(fixed_indices[1], fixed_indices[2], free_index, dims)
            push!(row_indices, flat)
            push!(col_indices, col)
            push!(values, value)
        end
    end
    coefficient_matrix = _nested_sparse_coefficient_map(
        row_indices,
        col_indices,
        values,
        prod(dims),
        size(side.coefficient_matrix, 2),
    )
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
    coefficient_matrix = _nested_sparse_coefficient_map([flat], [1], [1.0], prod(dims), 1)
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
    return _nested_sparse_coefficient_map(
        Int.(support_indices),
        collect(1:length(support_indices)),
        ones(Float64, length(support_indices)),
        nparent,
        length(support_indices),
    )
end

function _nested_product_coefficients(
    x_side::_CartesianNestedDoSide1D,
    y_side::_CartesianNestedDoSide1D,
    z_side::_CartesianNestedDoSide1D,
    dims::NTuple{3,Int},
)
    nparent = prod(dims)
    ncols = size(x_side.coefficient_matrix, 2) * size(y_side.coefficient_matrix, 2) * size(z_side.coefficient_matrix, 2)
    row_indices = Int[]
    col_indices = Int[]
    values = Float64[]
    column = 0
    for ixcol in 1:size(x_side.coefficient_matrix, 2),
        iycol in 1:size(y_side.coefficient_matrix, 2),
        izcol in 1:size(z_side.coefficient_matrix, 2)
        column += 1
        for (ix_local, ix) in enumerate(x_side.interval)
            vx = Float64(x_side.local_coefficients[ix_local, ixcol])
            iszero(vx) && continue
            for (iy_local, iy) in enumerate(y_side.interval)
                vy = Float64(y_side.local_coefficients[iy_local, iycol])
                iszero(vy) && continue
                for (iz_local, iz) in enumerate(z_side.interval)
                    vz = Float64(z_side.local_coefficients[iz_local, izcol])
                    iszero(vz) && continue
                    flat = _cartesian_flat_index(ix, iy, iz, dims)
                    push!(row_indices, flat)
                    push!(col_indices, column)
                    push!(values, vx * vy * vz)
                end
            end
        end
    end
    return _nested_sparse_coefficient_map(row_indices, col_indices, values, nparent, ncols)
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
    row_indices = Int[]
    col_indices = Int[]
    values = Float64[]
    column = 0
    for ifirst in 1:nfirst, isecond in 1:nsecond
        column += 1
        for (local_first, index_first) in enumerate(side_first.interval)
            value_first = Float64(side_first.local_coefficients[local_first, ifirst])
            iszero(value_first) && continue
            for (local_second, index_second) in enumerate(side_second.interval)
                value_second = Float64(side_second.local_coefficients[local_second, isecond])
                iszero(value_second) && continue
                flat_index =
                    face_kind == :xy ? _cartesian_flat_index(index_first, index_second, fixed_index, dims) :
                    face_kind == :xz ? _cartesian_flat_index(index_first, fixed_index, index_second, dims) :
                    _cartesian_flat_index(fixed_index, index_first, index_second, dims)
                push!(row_indices, flat_index)
                push!(col_indices, column)
                push!(values, value_first * value_second)
            end
        end
    end
    coefficients = _nested_sparse_coefficient_map(
        row_indices,
        col_indices,
        values,
        prod(dims),
        nfirst * nsecond,
    )
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

function _nested_normalize_term_storage(term_storage::Symbol)
    term_storage == :compact_production || throw(
        ArgumentError("nested term storage must be :compact_production"),
    )
    return term_storage
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

function _nested_fill_factorized_weighted_term_sum!(
    destination::AbstractMatrix{<:Real},
    factorized_basis::_CartesianNestedFactorizedBasis3D,
    term_coefficients::AbstractVector{<:Real},
    operator_terms_x::Array{Float64,3},
    operator_terms_y::Array{Float64,3},
    operator_terms_z::Array{Float64,3};
    include_basis_amplitudes::Bool,
)
    nterms = length(term_coefficients)
    nbasis = length(factorized_basis.basis_triplets)
    size(destination) == (nbasis, nbasis) || throw(
        ArgumentError("nested factorized weighted term-sum fill requires square output sized to the retained fixed basis"),
    )
    nterms == size(operator_terms_x, 1) == size(operator_terms_y, 1) == size(operator_terms_z, 1) || throw(
        ArgumentError("nested factorized weighted term-sum fill requires matching term counts"),
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
            for term in 1:nterms
                value +=
                    Float64(term_coefficients[term]) *
                    operator_terms_x[term, xi, xj] *
                    operator_terms_y[term, yi, yj] *
                    operator_terms_z[term, zi, zj]
            end
            value *= scale
            destination[row, column] = value
            destination[column, row] = value
        end
    end
    return destination
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
    include_basis_amplitudes::Bool = true,
)
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
    return matrix
end

function _nested_factorized_sum_of_products(
    factorized_basis::_CartesianNestedFactorizedBasis3D,
    terms,
    include_basis_amplitudes::Bool = true,
)
    nbasis = length(factorized_basis.basis_triplets)
    matrix = Matrix{Float64}(undef, nbasis, nbasis)
    _nested_fill_factorized_sum_of_products!(
        matrix,
        factorized_basis,
        terms;
        include_basis_amplitudes = include_basis_amplitudes,
    )
    return matrix
end

function _nested_factorized_gaussian_terms(
    factorized_basis::_CartesianNestedFactorizedBasis3D,
    gaussian_terms_x::Array{Float64,3},
    gaussian_terms_y::Array{Float64,3},
    gaussian_terms_z::Array{Float64,3},
    ;
    term_storage::Symbol = :compact_production,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    _nested_normalize_term_storage(term_storage)
    axis_term_tables_x, axis_term_tables_y, axis_term_tables_z, nbasis = @timeg "diatomic.packet.gaussian_terms.setup" begin
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
        (
            axis_term_tables_x,
            axis_term_tables_y,
            axis_term_tables_z,
            nbasis,
        )
    end
    isnothing(term_coefficients) && throw(
        ArgumentError("compact nested Gaussian-term storage requires explicit term coefficients"),
    )
    gaussian_sum = @timeg "diatomic.packet.gaussian_terms.contract" begin
        gaussian_sum_local = zeros(Float64, nbasis, nbasis)
        _nested_fill_factorized_weighted_term_sum!(
            gaussian_sum_local,
            factorized_basis,
            term_coefficients,
            axis_term_tables_x,
            axis_term_tables_y,
            axis_term_tables_z;
            include_basis_amplitudes = true,
        )
        gaussian_sum_local
    end
    return (
        gaussian_terms = nothing,
        gaussian_sum = gaussian_sum,
    )
end

function _nested_factorized_weight_aware_pair_terms(
    factorized_basis::_CartesianNestedFactorizedBasis3D,
    weights_x::AbstractVector{<:Real},
    weights_y::AbstractVector{<:Real},
    weights_z::AbstractVector{<:Real},
    pair_terms_x::Array{Float64,3},
    pair_terms_y::Array{Float64,3},
    pair_terms_z::Array{Float64,3},
    ;
    term_storage::Symbol = :compact_production,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    _nested_normalize_term_storage(term_storage)
    axis_weight_x, axis_weight_y, axis_weight_z, basis_weights = @timeg "diatomic.packet.pair_terms.weights" begin
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
        (
            axis_weight_x,
            axis_weight_y,
            axis_weight_z,
            basis_weights,
        )
    end
    axis_term_tables_x, axis_term_tables_y, axis_term_tables_z = @timeg "diatomic.packet.pair_terms.coefficients" begin
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
        (
            axis_term_tables_x,
            axis_term_tables_y,
            axis_term_tables_z,
        )
    end
    nbasis = length(factorized_basis.basis_triplets)
    isnothing(term_coefficients) && throw(
        ArgumentError("compact nested pair-term storage requires explicit term coefficients"),
    )
    pair_sum = @timeg "diatomic.packet.pair_terms.contract" begin
        pair_sum_local = zeros(Float64, nbasis, nbasis)
        _nested_fill_factorized_weighted_term_sum!(
            pair_sum_local,
            factorized_basis,
            term_coefficients,
            axis_term_tables_x,
            axis_term_tables_y,
            axis_term_tables_z;
            include_basis_amplitudes = false,
        )
        pair_sum_local
    end
    return (
        weights = basis_weights,
        pair_terms = nothing,
        pair_sum = pair_sum,
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

function _nested_fill_support_weighted_term_sum_matrix!(
    workspace::AbstractMatrix{<:Real},
    support_axes::_CartesianNestedSupportAxes3D,
    term_coefficients::AbstractVector{<:Real},
    operator_terms_x::Array{Float64,3},
    operator_terms_y::Array{Float64,3},
    operator_terms_z::Array{Float64,3};
    assume_symmetric::Bool = false,
)
    nsupport = length(support_axes.x)
    size(workspace) == (nsupport, nsupport) || throw(
        ArgumentError("nested support weighted term-sum workspace must have size ($(nsupport), $(nsupport))"),
    )
    nterms = size(operator_terms_x, 1)
    nterms == size(operator_terms_y, 1) == size(operator_terms_z, 1) || throw(
        ArgumentError("nested support weighted term sum requires matching term counts"),
    )
    length(term_coefficients) == nterms || throw(
        DimensionMismatch("nested support weighted term sum requires one coefficient per term slice"),
    )
    xstates = support_axes.x
    ystates = support_axes.y
    zstates = support_axes.z
    if assume_symmetric
        @inbounds for col in 1:nsupport
            jx = xstates[col]
            jy = ystates[col]
            jz = zstates[col]
            for row in 1:col
                ix = xstates[row]
                iy = ystates[row]
                iz = zstates[row]
                value = 0.0
                for term in 1:nterms
                    value +=
                        Float64(term_coefficients[term]) *
                        Float64(operator_terms_x[term, ix, jx]) *
                        Float64(operator_terms_y[term, iy, jy]) *
                        Float64(operator_terms_z[term, iz, jz])
                end
                workspace[row, col] = value
                workspace[col, row] = value
            end
        end
    else
        @inbounds for row in 1:nsupport
            ix = xstates[row]
            iy = ystates[row]
            iz = zstates[row]
            for col in 1:nsupport
                value = 0.0
                for term in 1:nterms
                    value +=
                        Float64(term_coefficients[term]) *
                        Float64(operator_terms_x[term, ix, xstates[col]]) *
                        Float64(operator_terms_y[term, iy, ystates[col]]) *
                        Float64(operator_terms_z[term, iz, zstates[col]])
                end
                workspace[row, col] = value
            end
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

function _nested_fill_support_product_matrix_chunk!(
    workspace::AbstractMatrix{<:Real},
    support_axes::_CartesianNestedSupportAxes3D,
    row_range::UnitRange{Int},
    operator_x::AbstractMatrix{<:Real},
    operator_y::AbstractMatrix{<:Real},
    operator_z::AbstractMatrix{<:Real},
)
    nsupport = length(support_axes.x)
    size(workspace) == (length(row_range), nsupport) || throw(
        ArgumentError("nested support chunk workspace must have size ($(length(row_range)), $(nsupport))"),
    )
    xstates = support_axes.x
    ystates = support_axes.y
    zstates = support_axes.z
    @inbounds for (local_row, row) in enumerate(row_range)
        ix = xstates[row]
        iy = ystates[row]
        iz = zstates[row]
        for col in 1:nsupport
            workspace[local_row, col] =
                Float64(operator_x[ix, xstates[col]]) *
                Float64(operator_y[iy, ystates[col]]) *
                Float64(operator_z[iz, zstates[col]])
        end
    end
    return workspace
end

function _nested_fill_support_weighted_term_sum_matrix_chunk!(
    workspace::AbstractMatrix{<:Real},
    support_axes::_CartesianNestedSupportAxes3D,
    row_range::UnitRange{Int},
    term_coefficients::AbstractVector{<:Real},
    operator_terms_x::Array{Float64,3},
    operator_terms_y::Array{Float64,3},
    operator_terms_z::Array{Float64,3},
)
    nsupport = length(support_axes.x)
    size(workspace) == (length(row_range), nsupport) || throw(
        ArgumentError("nested support weighted term-sum chunk workspace must have size ($(length(row_range)), $(nsupport))"),
    )
    nterms = size(operator_terms_x, 1)
    nterms == size(operator_terms_y, 1) == size(operator_terms_z, 1) || throw(
        ArgumentError("nested support weighted term sum requires matching term counts"),
    )
    length(term_coefficients) == nterms || throw(
        DimensionMismatch("nested support weighted term sum requires one coefficient per term slice"),
    )
    xstates = support_axes.x
    ystates = support_axes.y
    zstates = support_axes.z
    @inbounds for (local_row, row) in enumerate(row_range)
        ix = xstates[row]
        iy = ystates[row]
        iz = zstates[row]
        for col in 1:nsupport
            value = 0.0
            for term in 1:nterms
                value +=
                    Float64(term_coefficients[term]) *
                    Float64(operator_terms_x[term, ix, xstates[col]]) *
                    Float64(operator_terms_y[term, iy, ystates[col]]) *
                    Float64(operator_terms_z[term, iz, zstates[col]])
            end
            workspace[local_row, col] = value
        end
    end
    return workspace
end

function _nested_support_reference_weighted_term_sum(
    support_axes::_CartesianNestedSupportAxes3D,
    support_coefficients::AbstractMatrix{<:Real},
    support_workspace::AbstractMatrix{<:Real},
    contraction_scratch::AbstractMatrix{<:Real},
    term_coefficients::AbstractVector{<:Real},
    operator_terms_x::Array{Float64,3},
    operator_terms_y::Array{Float64,3},
    operator_terms_z::Array{Float64,3},
    sum_fill_label::AbstractString,
    project_label::AbstractString;
    assume_symmetric::Bool = false,
)
    nfixed = size(support_coefficients, 2)
    destination = Matrix{Float64}(undef, nfixed, nfixed)
    @timeg sum_fill_label begin
        _nested_fill_support_weighted_term_sum_matrix!(
            support_workspace,
            support_axes,
            term_coefficients,
            operator_terms_x,
            operator_terms_y,
            operator_terms_z;
            assume_symmetric = assume_symmetric,
        )
    end
    @timeg project_label begin
        mul!(contraction_scratch, transpose(support_coefficients), support_workspace)
        mul!(destination, contraction_scratch, support_coefficients, 1.0, 0.0)
    end
    return destination
end

function _nested_support_reference_weighted_term_sum(
    support_axes::_CartesianNestedSupportAxes3D,
    support_coefficients::SparseArrays.SparseMatrixCSC{Float64,Int},
    _support_workspace::AbstractMatrix{<:Real},
    _contraction_scratch::AbstractMatrix{<:Real},
    term_coefficients::AbstractVector{<:Real},
    operator_terms_x::Array{Float64,3},
    operator_terms_y::Array{Float64,3},
    operator_terms_z::Array{Float64,3},
    sum_fill_label::AbstractString,
    project_label::AbstractString;
    assume_symmetric::Bool = false,
)
    assume_symmetric
    nsupport = length(support_axes.x)
    nfixed = size(support_coefficients, 2)
    destination = Matrix{Float64}(undef, nfixed, nfixed)
    chunk_rows = _nested_support_reference_chunk_row_count(support_coefficients, nsupport, nfixed)
    workspace_chunk = Matrix{Float64}(undef, chunk_rows, nsupport)
    contraction_chunk = Matrix{Float64}(undef, chunk_rows, nfixed)
    first_chunk = true
    row_start = 1
    while row_start <= nsupport
        row_stop = min(nsupport, row_start + chunk_rows - 1)
        row_range = row_start:row_stop
        workspace_view = @view(workspace_chunk[1:length(row_range), :])
        contraction_view = @view(contraction_chunk[1:length(row_range), :])
        @timeg sum_fill_label begin
            _nested_fill_support_weighted_term_sum_matrix_chunk!(
                workspace_view,
                support_axes,
                row_range,
                term_coefficients,
                operator_terms_x,
                operator_terms_y,
                operator_terms_z,
            )
        end
        @timeg project_label begin
            mul!(contraction_view, workspace_view, support_coefficients)
            support_coefficients_chunk = support_coefficients[row_range, :]
            mul!(
                destination,
                transpose(support_coefficients_chunk),
                contraction_view,
                1.0,
                first_chunk ? 0.0 : 1.0,
            )
        end
        first_chunk = false
        row_start = row_stop + 1
    end
    return destination
end

function _nested_support_reference_gaussian_sum(
    support_axes::_CartesianNestedSupportAxes3D,
    support_coefficients::AbstractMatrix{<:Real},
    support_workspace::AbstractMatrix{<:Real},
    contraction_scratch::AbstractMatrix{<:Real},
    term_coefficients::AbstractVector{<:Real},
    operator_terms_x::Array{Float64,3},
    operator_terms_y::Array{Float64,3},
    operator_terms_z::Array{Float64,3},
)
    return _nested_support_reference_weighted_term_sum(
        support_axes,
        support_coefficients,
        support_workspace,
        contraction_scratch,
        term_coefficients,
        operator_terms_x,
        operator_terms_y,
        operator_terms_z,
        "diatomic.packet.gaussian_terms.sum_fill",
        "diatomic.packet.gaussian_terms.project";
        assume_symmetric = true,
    )
end

function _nested_support_reference_pair_sum(
    support_axes::_CartesianNestedSupportAxes3D,
    support_coefficients::AbstractMatrix{<:Real},
    support_workspace::AbstractMatrix{<:Real},
    contraction_scratch::AbstractMatrix{<:Real},
    term_coefficients::AbstractVector{<:Real},
    operator_terms_x::Array{Float64,3},
    operator_terms_y::Array{Float64,3},
    operator_terms_z::Array{Float64,3},
)
    return _nested_support_reference_weighted_term_sum(
        support_axes,
        support_coefficients,
        support_workspace,
        contraction_scratch,
        term_coefficients,
        operator_terms_x,
        operator_terms_y,
        operator_terms_z,
        "diatomic.packet.pair_terms.sum_fill",
        "diatomic.packet.pair_terms.project";
        assume_symmetric = true,
    )
end

function _nested_contract_support_product_sparse!(
    destination::AbstractMatrix{<:Real},
    support_axes::_CartesianNestedSupportAxes3D,
    support_coefficients::SparseArrays.SparseMatrixCSC{Float64,Int},
    operator_x::AbstractMatrix{<:Real},
    operator_y::AbstractMatrix{<:Real},
    operator_z::AbstractMatrix{<:Real};
    alpha::Float64 = 1.0,
    beta::Float64 = 0.0,
)
    nsupport = length(support_axes.x)
    nfixed = size(support_coefficients, 2)
    size(destination) == (nfixed, nfixed) || throw(
        ArgumentError("nested support contraction destination must have size ($(nfixed), $(nfixed))"),
    )
    chunk_rows = _nested_support_reference_chunk_row_count(support_coefficients, nsupport, nfixed)
    workspace_chunk = Matrix{Float64}(undef, chunk_rows, nsupport)
    contraction_chunk = Matrix{Float64}(undef, chunk_rows, nfixed)
    first_chunk = true
    row_start = 1
    while row_start <= nsupport
        row_stop = min(nsupport, row_start + chunk_rows - 1)
        row_range = row_start:row_stop
        workspace_view = @view(workspace_chunk[1:length(row_range), :])
        contraction_view = @view(contraction_chunk[1:length(row_range), :])
        _nested_fill_support_product_matrix_chunk!(
            workspace_view,
            support_axes,
            row_range,
            operator_x,
            operator_y,
            operator_z,
        )
        mul!(contraction_view, workspace_view, support_coefficients)
        support_coefficients_chunk = support_coefficients[row_range, :]
        mul!(
            destination,
            transpose(support_coefficients_chunk),
            contraction_view,
            alpha,
            first_chunk ? beta : 1.0,
        )
        first_chunk = false
        row_start = row_stop + 1
    end
    return destination
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
    _workspace::AbstractMatrix{<:Real},
    _contraction_scratch::AbstractMatrix{<:Real},
    support_states::AbstractVector{<:NTuple{3,Int}},
    support_coefficients::SparseArrays.SparseMatrixCSC{Float64,Int},
    operator_x::AbstractMatrix{<:Real},
    operator_y::AbstractMatrix{<:Real},
    operator_z::AbstractMatrix{<:Real};
    alpha::Float64 = 1.0,
    beta::Float64 = 0.0,
    assume_symmetric::Bool = false,
)
    assume_symmetric
    return _nested_contract_support_product_sparse!(
        destination,
        _nested_support_axes(support_states),
        support_coefficients,
        operator_x,
        operator_y,
        operator_z;
        alpha = alpha,
        beta = beta,
    )
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

function _nested_contract_support_product!(
    destination::AbstractMatrix{<:Real},
    _workspace::AbstractMatrix{<:Real},
    _contraction_scratch::AbstractMatrix{<:Real},
    support_axes::_CartesianNestedSupportAxes3D,
    support_coefficients::SparseArrays.SparseMatrixCSC{Float64,Int},
    operator_x::AbstractMatrix{<:Real},
    operator_y::AbstractMatrix{<:Real},
    operator_z::AbstractMatrix{<:Real};
    alpha::Float64 = 1.0,
    beta::Float64 = 0.0,
    assume_symmetric::Bool = false,
)
    assume_symmetric
    return _nested_contract_support_product_sparse!(
        destination,
        support_axes,
        support_coefficients,
        operator_x,
        operator_y,
        operator_z;
        alpha = alpha,
        beta = beta,
    )
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
    ;
    term_storage::Symbol = :compact_production,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    _nested_normalize_term_storage(term_storage)
    nterms = size(pgdg.pair_factor_terms, 1)
    nfixed = size(support_coefficients, 2)
    fixed_weights = @timeg "diatomic.packet.pair_terms.weights" begin
        support_weights = _nested_support_weights(support_states, pgdg.weights)
        fixed_weights_local = vec(transpose(support_coefficients) * support_weights)
        all(isfinite, fixed_weights_local) || throw(
            ArgumentError("nested fixed-block IDA transfer requires finite contracted integral weights"),
        )
        minimum(fixed_weights_local) > 0.0 || throw(
            ArgumentError("nested fixed-block IDA transfer requires positive contracted integral weights"),
        )
        fixed_weights_local
    end
    weighted_support_coefficients = @timeg "diatomic.packet.pair_terms.coefficients" begin
        support_coefficients .* reshape(1.0 ./ fixed_weights, 1, :)
    end
    isnothing(term_coefficients) && throw(
        ArgumentError("compact nested pair-term storage requires explicit term coefficients"),
    )
    pair_sum = _nested_support_reference_pair_sum(
        support_axes,
        weighted_support_coefficients,
        support_workspace,
        contraction_scratch,
        term_coefficients,
        pgdg.pair_factor_terms_raw,
        pgdg.pair_factor_terms_raw,
        pgdg.pair_factor_terms_raw,
    )

    return (
        weights = fixed_weights,
        pair_terms = nothing,
        pair_sum = pair_sum,
    )
end

function _nested_weight_aware_pair_terms(
    bundles::_CartesianNestedAxisBundles3D,
    support_states::AbstractVector{<:NTuple{3,Int}},
    support_axes::_CartesianNestedSupportAxes3D,
    support_coefficients::AbstractMatrix{<:Real},
    support_workspace::AbstractMatrix{<:Real},
    contraction_scratch::AbstractMatrix{<:Real},
    ;
    term_storage::Symbol = :compact_production,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    _nested_normalize_term_storage(term_storage)
    pgdg_x = _nested_axis_pgdg(bundles, :x)
    pgdg_y = _nested_axis_pgdg(bundles, :y)
    pgdg_z = _nested_axis_pgdg(bundles, :z)
    nterms = size(pgdg_x.pair_factor_terms, 1)
    nterms == size(pgdg_y.pair_factor_terms, 1) == size(pgdg_z.pair_factor_terms, 1) || throw(
        ArgumentError("mixed-axis nested IDA transfer requires the same Gaussian expansion term count on all axes"),
    )
    nfixed = size(support_coefficients, 2)
    fixed_weights = @timeg "diatomic.packet.pair_terms.weights" begin
        support_weights = _nested_support_weights(
            support_states,
            pgdg_x.weights,
            pgdg_y.weights,
            pgdg_z.weights,
        )
        fixed_weights_local = vec(transpose(support_coefficients) * support_weights)
        all(isfinite, fixed_weights_local) || throw(
            ArgumentError("mixed-axis nested fixed-block IDA transfer requires finite contracted integral weights"),
        )
        minimum(fixed_weights_local) > 0.0 || throw(
            ArgumentError("mixed-axis nested fixed-block IDA transfer requires positive contracted integral weights"),
        )
        fixed_weights_local
    end
    weighted_support_coefficients = @timeg "diatomic.packet.pair_terms.coefficients" begin
        support_coefficients .* reshape(1.0 ./ fixed_weights, 1, :)
    end

    raw_pair_terms_x = pgdg_x.pair_factor_terms_raw
    raw_pair_terms_y = pgdg_y.pair_factor_terms_raw
    raw_pair_terms_z = pgdg_z.pair_factor_terms_raw
    isnothing(term_coefficients) && throw(
        ArgumentError("compact nested pair-term storage requires explicit term coefficients"),
    )
    pair_sum = _nested_support_reference_pair_sum(
        support_axes,
        weighted_support_coefficients,
        support_workspace,
        contraction_scratch,
        term_coefficients,
        raw_pair_terms_x,
        raw_pair_terms_y,
        raw_pair_terms_z,
    )

    return (
        weights = fixed_weights,
        pair_terms = nothing,
        pair_sum = pair_sum,
    )
end

function _nested_weight_aware_pair_terms(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    support_states::AbstractVector{<:NTuple{3,Int}},
    support_coefficients::AbstractMatrix{<:Real},
    ;
    term_storage::Symbol = :compact_production,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    nsupport = length(support_states)
    nfixed = size(support_coefficients, 2)
    support_workspace, contraction_scratch =
        _nested_support_reference_workspaces(support_coefficients, nsupport, nfixed)
    support_axes = _nested_support_axes(support_states)
    return _nested_weight_aware_pair_terms(
        pgdg,
        support_states,
        support_axes,
        support_coefficients,
        support_workspace,
        contraction_scratch,
        term_storage = term_storage,
        term_coefficients = term_coefficients,
    )
end

function _nested_weight_aware_pair_terms(
    bundles::_CartesianNestedAxisBundles3D,
    support_states::AbstractVector{<:NTuple{3,Int}},
    support_coefficients::AbstractMatrix{<:Real},
    ;
    term_storage::Symbol = :compact_production,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    nsupport = length(support_states)
    nfixed = size(support_coefficients, 2)
    support_workspace, contraction_scratch =
        _nested_support_reference_workspaces(support_coefficients, nsupport, nfixed)
    support_axes = _nested_support_axes(support_states)
    return _nested_weight_aware_pair_terms(
        bundles,
        support_states,
        support_axes,
        support_coefficients,
        support_workspace,
        contraction_scratch,
        term_storage = term_storage,
        term_coefficients = term_coefficients,
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
    ;
    packet_kernel::Symbol = :support_reference,
    term_storage::Symbol = :compact_production,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    return @timeg "diatomic.packet.total" begin
        packet_kernel, term_storage_value, support_states, nshell, nsupport, nterms, support_axes, support_coefficients, factorized_basis, factorized_base_tables, support_workspace, contraction_scratch = @timeg "diatomic.packet.setup" begin
            packet_kernel = _nested_normalize_packet_kernel(packet_kernel)
            term_storage_value = _nested_normalize_term_storage(term_storage)
            support_states = [_cartesian_unflat_index(index, size(pgdg.overlap, 1)) for index in support_indices]
            nshell = size(coefficient_matrix, 2)
            nsupport = length(support_states)
            nterms = size(pgdg.gaussian_factor_terms, 1)
            support_axes = packet_kernel == :support_reference ? _nested_support_axes(support_states) : nothing
            support_coefficients =
                packet_kernel == :support_reference ?
                _nested_support_coefficient_slice(coefficient_matrix, support_indices) :
                nothing
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
            support_workspace, contraction_scratch =
                packet_kernel == :support_reference ?
                _nested_support_reference_workspaces(support_coefficients, nsupport, nshell) :
                (Matrix{Float64}(undef, 0, 0), Matrix{Float64}(undef, 0, 0))
            (
                packet_kernel,
                term_storage_value,
                support_states,
                nshell,
                nsupport,
                nterms,
                support_axes,
                support_coefficients,
                factorized_basis,
                factorized_base_tables,
                support_workspace,
                contraction_scratch,
            )
        end

        overlap = @timeg "diatomic.packet.base.overlap" begin
            packet_kernel == :factorized_direct ?
            _nested_factorized_product_matrix(
                factorized_basis,
                factorized_base_tables.overlap,
                factorized_base_tables.overlap,
                factorized_base_tables.overlap,
            ) :
            begin
                overlap_local = Matrix{Float64}(undef, nshell, nshell)
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
                overlap_local
            end
        end

        kinetic = @timeg "diatomic.packet.base.kinetic" begin
            packet_kernel == :factorized_direct ?
            _nested_factorized_sum_of_products(
                factorized_basis,
                (
                    (factorized_base_tables.kinetic, factorized_base_tables.overlap, factorized_base_tables.overlap),
                    (factorized_base_tables.overlap, factorized_base_tables.kinetic, factorized_base_tables.overlap),
                    (factorized_base_tables.overlap, factorized_base_tables.overlap, factorized_base_tables.kinetic),
                ),
            ) :
            begin
                kinetic_local = Matrix{Float64}(undef, nshell, nshell)
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
                kinetic_local
            end
        end

        position_x = @timeg "diatomic.packet.base.position_x" begin
            packet_kernel == :factorized_direct ?
            _nested_factorized_product_matrix(
                factorized_basis,
                factorized_base_tables.position,
                factorized_base_tables.overlap,
                factorized_base_tables.overlap,
            ) :
            begin
                position_x_local = Matrix{Float64}(undef, nshell, nshell)
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
                position_x_local
            end
        end

        position_y = @timeg "diatomic.packet.base.position_y" begin
            packet_kernel == :factorized_direct ?
            _nested_factorized_product_matrix(
                factorized_basis,
                factorized_base_tables.overlap,
                factorized_base_tables.position,
                factorized_base_tables.overlap,
            ) :
            begin
                position_y_local = Matrix{Float64}(undef, nshell, nshell)
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
                position_y_local
            end
        end

        position_z = @timeg "diatomic.packet.base.position_z" begin
            packet_kernel == :factorized_direct ?
            _nested_factorized_product_matrix(
                factorized_basis,
                factorized_base_tables.overlap,
                factorized_base_tables.overlap,
                factorized_base_tables.position,
            ) :
            begin
                position_z_local = Matrix{Float64}(undef, nshell, nshell)
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
                position_z_local
            end
        end

        x2_x = @timeg "diatomic.packet.base.x2_x" begin
            packet_kernel == :factorized_direct ?
            _nested_factorized_product_matrix(
                factorized_basis,
                factorized_base_tables.x2,
                factorized_base_tables.overlap,
                factorized_base_tables.overlap,
            ) :
            begin
                x2_x_local = Matrix{Float64}(undef, nshell, nshell)
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
                x2_x_local
            end
        end

        x2_y = @timeg "diatomic.packet.base.x2_y" begin
            packet_kernel == :factorized_direct ?
            _nested_factorized_product_matrix(
                factorized_basis,
                factorized_base_tables.overlap,
                factorized_base_tables.x2,
                factorized_base_tables.overlap,
            ) :
            begin
                x2_y_local = Matrix{Float64}(undef, nshell, nshell)
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
                x2_y_local
            end
        end

        x2_z = @timeg "diatomic.packet.base.x2_z" begin
            packet_kernel == :factorized_direct ?
            _nested_factorized_product_matrix(
                factorized_basis,
                factorized_base_tables.overlap,
                factorized_base_tables.overlap,
                factorized_base_tables.x2,
            ) :
            begin
                x2_z_local = Matrix{Float64}(undef, nshell, nshell)
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
                x2_z_local
            end
        end

        pair_data = @timeg "diatomic.packet.pair_terms" begin
            packet_kernel == :factorized_direct ?
            _nested_factorized_weight_aware_pair_terms(
                factorized_basis,
                pgdg.weights,
                pgdg.weights,
                pgdg.weights,
                pgdg.pair_factor_terms_raw,
                pgdg.pair_factor_terms_raw,
                pgdg.pair_factor_terms_raw,
                term_storage = term_storage_value,
                term_coefficients = term_coefficients,
            ) :
            _nested_weight_aware_pair_terms(
                pgdg,
                support_states,
                support_axes,
                support_coefficients,
                support_workspace,
                contraction_scratch,
                term_storage = term_storage_value,
                term_coefficients = term_coefficients,
            )
        end
        gaussian_terms = @timeg "diatomic.packet.gaussian_terms" begin
            packet_kernel == :factorized_direct ?
            _nested_factorized_gaussian_terms(
                factorized_basis,
                pgdg.gaussian_factor_terms,
                pgdg.gaussian_factor_terms,
                pgdg.gaussian_factor_terms,
                term_storage = term_storage_value,
                term_coefficients = term_coefficients,
            ) :
            begin
                isnothing(term_coefficients) && throw(
                    ArgumentError("compact nested Gaussian-term storage requires explicit term coefficients"),
                )
                gaussian_sum = _nested_support_reference_gaussian_sum(
                    support_axes,
                    support_coefficients,
                    support_workspace,
                    contraction_scratch,
                    term_coefficients,
                    pgdg.gaussian_factor_terms,
                    pgdg.gaussian_factor_terms,
                    pgdg.gaussian_factor_terms,
                )
                (
                    gaussian_terms = nothing,
                    gaussian_sum = gaussian_sum,
                )
            end
        end

        (
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
                gaussian_terms.gaussian_sum,
                pair_data.pair_sum,
                gaussian_terms.gaussian_terms,
                pair_data.pair_terms,
                term_storage_value,
            ),
            support_states = support_states,
        )
    end
end

function _nested_shell_packet(
    bundles::_CartesianNestedAxisBundles3D,
    coefficient_matrix::AbstractMatrix{<:Real},
    support_indices::AbstractVector{Int},
    ;
    packet_kernel::Symbol = :support_reference,
    term_storage::Symbol = :compact_production,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    return @timeg "diatomic.packet.total" begin
        packet_kernel, term_storage_value, pgdg_x, pgdg_y, pgdg_z, support_states, nshell, nsupport, nterms, support_axes, support_coefficients, factorized_basis, factorized_base_tables_x, factorized_base_tables_y, factorized_base_tables_z, support_workspace, contraction_scratch = @timeg "diatomic.packet.setup" begin
            packet_kernel = _nested_normalize_packet_kernel(packet_kernel)
            term_storage_value = _nested_normalize_term_storage(term_storage)
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
                packet_kernel == :support_reference ?
                _nested_support_coefficient_slice(coefficient_matrix, support_indices) :
                nothing
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
            support_workspace, contraction_scratch =
                packet_kernel == :support_reference ?
                _nested_support_reference_workspaces(support_coefficients, nsupport, nshell) :
                (Matrix{Float64}(undef, 0, 0), Matrix{Float64}(undef, 0, 0))
            (
                packet_kernel,
                term_storage_value,
                pgdg_x,
                pgdg_y,
                pgdg_z,
                support_states,
                nshell,
                nsupport,
                nterms,
                support_axes,
                support_coefficients,
                factorized_basis,
                factorized_base_tables_x,
                factorized_base_tables_y,
                factorized_base_tables_z,
                support_workspace,
                contraction_scratch,
            )
        end

        overlap = @timeg "diatomic.packet.base.overlap" begin
            packet_kernel == :factorized_direct ?
            _nested_factorized_product_matrix(
                factorized_basis,
                factorized_base_tables_x.overlap,
                factorized_base_tables_y.overlap,
                factorized_base_tables_z.overlap,
            ) :
            begin
                overlap_local = Matrix{Float64}(undef, nshell, nshell)
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
                overlap_local
            end
        end

        kinetic = @timeg "diatomic.packet.base.kinetic" begin
            packet_kernel == :factorized_direct ?
            _nested_factorized_sum_of_products(
                factorized_basis,
                (
                    (factorized_base_tables_x.kinetic, factorized_base_tables_y.overlap, factorized_base_tables_z.overlap),
                    (factorized_base_tables_x.overlap, factorized_base_tables_y.kinetic, factorized_base_tables_z.overlap),
                    (factorized_base_tables_x.overlap, factorized_base_tables_y.overlap, factorized_base_tables_z.kinetic),
                ),
            ) :
            begin
                kinetic_local = Matrix{Float64}(undef, nshell, nshell)
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
                kinetic_local
            end
        end

        position_x = @timeg "diatomic.packet.base.position_x" begin
            packet_kernel == :factorized_direct ?
            _nested_factorized_product_matrix(
                factorized_basis,
                factorized_base_tables_x.position,
                factorized_base_tables_y.overlap,
                factorized_base_tables_z.overlap,
            ) :
            begin
                position_x_local = Matrix{Float64}(undef, nshell, nshell)
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
                position_x_local
            end
        end

        position_y = @timeg "diatomic.packet.base.position_y" begin
            packet_kernel == :factorized_direct ?
            _nested_factorized_product_matrix(
                factorized_basis,
                factorized_base_tables_x.overlap,
                factorized_base_tables_y.position,
                factorized_base_tables_z.overlap,
            ) :
            begin
                position_y_local = Matrix{Float64}(undef, nshell, nshell)
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
                position_y_local
            end
        end

        position_z = @timeg "diatomic.packet.base.position_z" begin
            packet_kernel == :factorized_direct ?
            _nested_factorized_product_matrix(
                factorized_basis,
                factorized_base_tables_x.overlap,
                factorized_base_tables_y.overlap,
                factorized_base_tables_z.position,
            ) :
            begin
                position_z_local = Matrix{Float64}(undef, nshell, nshell)
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
                position_z_local
            end
        end

        x2_x = @timeg "diatomic.packet.base.x2_x" begin
            packet_kernel == :factorized_direct ?
            _nested_factorized_product_matrix(
                factorized_basis,
                factorized_base_tables_x.x2,
                factorized_base_tables_y.overlap,
                factorized_base_tables_z.overlap,
            ) :
            begin
                x2_x_local = Matrix{Float64}(undef, nshell, nshell)
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
                x2_x_local
            end
        end

        x2_y = @timeg "diatomic.packet.base.x2_y" begin
            packet_kernel == :factorized_direct ?
            _nested_factorized_product_matrix(
                factorized_basis,
                factorized_base_tables_x.overlap,
                factorized_base_tables_y.x2,
                factorized_base_tables_z.overlap,
            ) :
            begin
                x2_y_local = Matrix{Float64}(undef, nshell, nshell)
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
                x2_y_local
            end
        end

        x2_z = @timeg "diatomic.packet.base.x2_z" begin
            packet_kernel == :factorized_direct ?
            _nested_factorized_product_matrix(
                factorized_basis,
                factorized_base_tables_x.overlap,
                factorized_base_tables_y.overlap,
                factorized_base_tables_z.x2,
            ) :
            begin
                x2_z_local = Matrix{Float64}(undef, nshell, nshell)
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
                x2_z_local
            end
        end

        pair_data = @timeg "diatomic.packet.pair_terms" begin
            packet_kernel == :factorized_direct ?
            _nested_factorized_weight_aware_pair_terms(
                factorized_basis,
                pgdg_x.weights,
                pgdg_y.weights,
                pgdg_z.weights,
                pgdg_x.pair_factor_terms_raw,
                pgdg_y.pair_factor_terms_raw,
                pgdg_z.pair_factor_terms_raw,
                term_storage = term_storage_value,
                term_coefficients = term_coefficients,
            ) :
            _nested_weight_aware_pair_terms(
                bundles,
                support_states,
                support_axes,
                support_coefficients,
                support_workspace,
                contraction_scratch,
                term_storage = term_storage_value,
                term_coefficients = term_coefficients,
            )
        end
        gaussian_terms = @timeg "diatomic.packet.gaussian_terms" begin
            packet_kernel == :factorized_direct ?
            _nested_factorized_gaussian_terms(
                factorized_basis,
                pgdg_x.gaussian_factor_terms,
                pgdg_y.gaussian_factor_terms,
                pgdg_z.gaussian_factor_terms,
                term_storage = term_storage_value,
                term_coefficients = term_coefficients,
            ) :
            begin
                isnothing(term_coefficients) && throw(
                    ArgumentError("compact nested Gaussian-term storage requires explicit term coefficients"),
                )
                gaussian_sum = _nested_support_reference_gaussian_sum(
                    support_axes,
                    support_coefficients,
                    support_workspace,
                    contraction_scratch,
                    term_coefficients,
                    pgdg_x.gaussian_factor_terms,
                    pgdg_y.gaussian_factor_terms,
                    pgdg_z.gaussian_factor_terms,
                )
                (
                    gaussian_terms = nothing,
                    gaussian_sum = gaussian_sum,
                )
            end
        end

        (
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
                gaussian_terms.gaussian_sum,
                pair_data.pair_sum,
                gaussian_terms.gaussian_terms,
                pair_data.pair_terms,
                term_storage_value,
            ),
            support_states = support_states,
        )
    end
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
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
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
    coefficient_matrix = _nested_hcat_coefficient_maps([face_low.coefficient_matrix, face_high.coefficient_matrix])
    support_indices = _nested_shell_support_indices(faces)
    shell_data = _nested_shell_packet(
        pgdg,
        coefficient_matrix,
        support_indices;
        term_coefficients = term_coefficients,
    )
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
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    return _nested_xy_shell_pair(
        bundle.pgdg_intermediate,
        x_interval,
        y_interval;
        retain_x = retain_x,
        retain_y = retain_y,
        term_coefficients = term_coefficients,
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
    term_storage::Symbol = :compact_production,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
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

    coefficient_blocks = AbstractMatrix{Float64}[face.coefficient_matrix for face in faces]
    coefficient_matrix = _nested_hcat_coefficient_maps(coefficient_blocks)
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
        term_storage = term_storage,
        term_coefficients = term_coefficients,
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
    packet_kernel::Symbol = :support_reference,
    term_storage::Symbol = :compact_production,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    n1d = size(pgdg.overlap, 1)
    shell_faces, edges, corners, coefficient_matrix, face_column_ranges, edge_column_ranges, corner_column_ranges, support_indices = @timeg "shell_layer.nonpacket" begin
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
            term_storage = term_storage,
            term_coefficients = term_coefficients,
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

        coefficient_blocks = AbstractMatrix{Float64}[face.coefficient_matrix for face in shell_faces.faces]
        append!(coefficient_blocks, [edge.coefficient_matrix for edge in edges])
        append!(coefficient_blocks, [corner.coefficient_matrix for corner in corners])
        coefficient_matrix = _nested_hcat_coefficient_maps(coefficient_blocks)

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
        (
            shell_faces,
            edges,
            corners,
            coefficient_matrix,
            face_column_ranges,
            edge_column_ranges,
            corner_column_ranges,
            support_indices,
        )
    end
    shell_data = _nested_shell_packet(
        pgdg,
        coefficient_matrix,
        support_indices,
        packet_kernel = packet_kernel,
        term_storage = term_storage,
        term_coefficients = term_coefficients,
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
        packet.gaussian_sum,
        packet.pair_sum,
        packet.gaussian_terms,
        packet.pair_terms,
        packet.term_storage,
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
    term_storage::Symbol = :compact_production,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    n1d = size(pgdg.overlap, 1)
    core_indices = _nested_box_support_indices(x_interval, y_interval, z_interval, n1d)
    isempty(intersect(core_indices, shell.support_indices)) || throw(
        ArgumentError("nested shell-plus-core construction requires the direct core block to stay disjoint from the shell-face supports"),
    )
    core_coefficients = _nested_direct_core_coefficients(core_indices, n1d^3)
    coefficient_matrix = _nested_hcat_coefficient_maps([core_coefficients, shell.coefficient_matrix])
    support_indices = sort(vcat(core_indices, shell.support_indices))
    shell_data = _nested_shell_packet(
        pgdg,
        coefficient_matrix,
        support_indices;
        packet_kernel = packet_kernel,
        term_storage = term_storage,
        term_coefficients = term_coefficients,
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
    ;
    kwargs...,
)
    return _nested_shell_plus_core(bundle.pgdg_intermediate, shell, x_interval, y_interval, z_interval; kwargs...)
end

function _nested_shell_sequence(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D};
    enforce_coverage::Bool = true,
    packet_kernel::Symbol = :support_reference,
    term_storage::Symbol = :compact_production,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
    build_packet::Bool = true,
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
        packet_kernel = packet_kernel,
        term_storage = term_storage,
        term_coefficients = term_coefficients,
        build_packet = build_packet,
    )
end

function _nested_shell_sequence_from_core_block(
    pgdg::_MappedOrdinaryPGDGIntermediate1D,
    core_indices::AbstractVector{Int},
    core_coefficients::AbstractMatrix{<:Real},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D};
    enforce_coverage::Bool = true,
    packet_kernel::Symbol = :support_reference,
    term_storage::Symbol = :compact_production,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
    build_packet::Bool = true,
)
    n1d = size(pgdg.overlap, 1)
    support_indices, working_box, coefficient_matrix = @timeg "sequence_merge.nonpacket" begin
        support_indices = _nested_sequence_support_indices(core_indices, shell_layers)
        working_box =
            enforce_coverage ?
            _nested_assert_sequence_coverage(core_indices, shell_layers, support_indices, n1d) :
            _nested_sequence_working_box(core_indices, shell_layers, n1d)
        coefficient_blocks = _CartesianCoefficientMap[core_coefficients]
        append!(coefficient_blocks, [shell.coefficient_matrix for shell in shell_layers])
        coefficient_matrix = @timeg "diatomic.sequence.coefficient_concat" begin
            _nested_hcat_coefficient_maps(coefficient_blocks)
        end
        (support_indices, working_box, coefficient_matrix)
    end
    shell_data = if build_packet
        @timeg "diatomic.sequence.packet" begin
            _nested_shell_packet(
                pgdg,
                coefficient_matrix,
                support_indices,
                packet_kernel = packet_kernel,
                term_storage = term_storage,
                term_coefficients = term_coefficients,
            )
        end
    else
        nothing
    end

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
        isnothing(shell_data) ? nothing : shell_data.support_states,
        isnothing(shell_data) ? nothing : shell_data.packet,
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
    packet_kernel::Symbol = :support_reference,
    term_storage::Symbol = :compact_production,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
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

    coefficient_blocks = AbstractMatrix{Float64}[face.coefficient_matrix for face in faces]
    coefficient_matrix = _nested_hcat_coefficient_maps(coefficient_blocks)
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
        packet_kernel = packet_kernel,
        term_storage = term_storage,
        term_coefficients = term_coefficients,
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
    packet_kernel::Symbol = :support_reference,
    term_storage::Symbol = :compact_production,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    dims = _nested_axis_lengths(bundles)
    shell_faces, edges, corners, coefficient_matrix, face_column_ranges, edge_column_ranges, corner_column_ranges, support_indices = @timeg "shell_layer.nonpacket" begin
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
            term_storage = term_storage,
            term_coefficients = term_coefficients,
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

        coefficient_blocks = AbstractMatrix{Float64}[face.coefficient_matrix for face in shell_faces.faces]
        append!(coefficient_blocks, [edge.coefficient_matrix for edge in edges])
        append!(coefficient_blocks, [corner.coefficient_matrix for corner in corners])
        coefficient_matrix = _nested_hcat_coefficient_maps(coefficient_blocks)

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
        (
            shell_faces,
            edges,
            corners,
            coefficient_matrix,
            face_column_ranges,
            edge_column_ranges,
            corner_column_ranges,
            support_indices,
        )
    end
    shell_data = _nested_shell_packet(
        bundles,
        coefficient_matrix,
        support_indices,
        packet_kernel = packet_kernel,
        term_storage = term_storage,
        term_coefficients = term_coefficients,
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
    packet_kernel::Symbol = :support_reference,
    term_storage::Symbol = :compact_production,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
    build_packet::Bool = true,
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
        packet_kernel = packet_kernel,
        term_storage = term_storage,
        term_coefficients = term_coefficients,
        build_packet = build_packet,
    )
end

function _nested_shell_sequence_from_core_block(
    bundles::_CartesianNestedAxisBundles3D,
    core_indices::AbstractVector{Int},
    core_coefficients::AbstractMatrix{<:Real},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D};
    enforce_coverage::Bool = true,
    packet_kernel::Symbol = :support_reference,
    term_storage::Symbol = :compact_production,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
    build_packet::Bool = true,
)
    dims = _nested_axis_lengths(bundles)
    support_indices, working_box, coefficient_matrix = @timeg "sequence_merge.nonpacket" begin
        support_indices = _nested_sequence_support_indices(core_indices, shell_layers)
        working_box =
            enforce_coverage ?
            _nested_assert_sequence_coverage(core_indices, shell_layers, support_indices, dims) :
            _nested_sequence_working_box(core_indices, shell_layers, dims)
        coefficient_blocks = _CartesianCoefficientMap[core_coefficients]
        append!(coefficient_blocks, [shell.coefficient_matrix for shell in shell_layers])
        coefficient_matrix = @timeg "diatomic.sequence.coefficient_concat" begin
            _nested_hcat_coefficient_maps(coefficient_blocks)
        end
        (support_indices, working_box, coefficient_matrix)
    end
    shell_data = if build_packet
        @timeg "diatomic.sequence.packet" begin
            _nested_shell_packet(
                bundles,
                coefficient_matrix,
                support_indices,
                packet_kernel = packet_kernel,
                term_storage = term_storage,
                term_coefficients = term_coefficients,
            )
        end
    else
        nothing
    end

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
        isnothing(shell_data) ? nothing : shell_data.support_states,
        isnothing(shell_data) ? nothing : shell_data.packet,
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
    packet_kernel::Symbol = :support_reference,
    term_storage::Symbol = :compact_production,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
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
                packet_kernel = packet_kernel,
                term_storage = term_storage,
                term_coefficients = term_coefficients,
            ),
        )
        current_box = inner_box
    end
    return _nested_shell_sequence(
        bundles,
        current_box...,
        shell_layers,
        packet_kernel = packet_kernel,
        term_storage = term_storage,
        term_coefficients = term_coefficients,
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

function _nested_axis_index(axis::Symbol)
    axis == :x && return 1
    axis == :y && return 2
    axis == :z && return 3
    throw(ArgumentError("nested axis lookup requires axis = :x, :y, or :z"))
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
    kwargs...,
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
        kwargs...,
    )
    return _nested_shell_sequence(
        pgdg,
        inner_x,
        inner_y,
        inner_z,
        [inner_shell],
        kwargs...,
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
        kwargs...,
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
    kwargs...,
)
    return _nested_shell_sequence(
        bundle.pgdg_intermediate,
        x_interval,
        y_interval,
        z_interval,
        shell_layers;
        enforce_coverage = enforce_coverage,
        kwargs...,
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
    kwargs...,
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
        kwargs...,
    )
end

function _nested_nside_shell_sequence(
    bundle::_MappedOrdinaryGausslet1DBundle,
    x_interval::UnitRange{Int},
    y_interval::UnitRange{Int},
    z_interval::UnitRange{Int},
    shell_layers::AbstractVector{<:_AbstractCartesianNestedShellLayer3D};
    nside::Int = 5,
    kwargs...,
)
    return _nested_nside_shell_sequence(
        bundle.pgdg_intermediate,
        x_interval,
        y_interval,
        z_interval,
        shell_layers;
        nside = nside,
        kwargs...,
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
        packet.gaussian_sum,
        packet.pair_sum,
        packet.gaussian_terms,
        packet.pair_terms,
        packet.term_storage,
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
        packet.gaussian_sum,
        packet.pair_sum,
        packet.gaussian_terms,
        packet.pair_terms,
        packet.term_storage,
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
    isnothing(shell.packet) && throw(
        ArgumentError("nested fixed-block construction requires a shell sequence with an assembled packet"),
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
        packet.gaussian_sum,
        packet.pair_sum,
        packet.gaussian_terms,
        packet.pair_terms,
        packet.term_storage,
        Matrix{Float64}(fixed_centers),
    )
end

function _nested_fixed_block(
    shell::_CartesianNestedShellSequence3D,
    parent_basis,
)
    isnothing(shell.packet) && throw(
        ArgumentError("nested fixed-block construction requires a shell sequence with an assembled packet"),
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
        packet.gaussian_sum,
        packet.pair_sum,
        packet.gaussian_terms,
        packet.pair_terms,
        packet.term_storage,
        Matrix{Float64}(fixed_centers),
    )
end

function _nested_fixed_block(
    shell::_CartesianNestedShellSequence3D,
    bundle::_MappedOrdinaryGausslet1DBundle,
)
    return _nested_fixed_block(shell, bundle.basis)
end
