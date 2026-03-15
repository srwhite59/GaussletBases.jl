"""
    HierarchicalBasisBox1D

Single node in a 1D hierarchical basis partition tree.
"""
struct HierarchicalBasisBox1D
    index::Int
    parent_index::Union{Nothing, Int}
    child_indices::Vector{Int}
    level::Int
    xmin::Float64
    xmax::Float64
    basis_indices::Vector{Int}
end

function Base.show(io::IO, box::HierarchicalBasisBox1D)
    print(
        io,
        "HierarchicalBasisBox1D(index=",
        box.index,
        ", level=",
        box.level,
        ", range=(",
        box.xmin,
        ", ",
        box.xmax,
        "), n=",
        length(box.basis_indices),
        ", children=",
        box.child_indices,
        ")",
    )
end

"""
    HierarchicalBasisPartition1D

Small parent-child tree built on top of `BasisPartition1D`.
"""
struct HierarchicalBasisPartition1D
    center_data::Vector{Float64}
    box_data::Vector{HierarchicalBasisBox1D}
    root_indices::Vector{Int}
end

function Base.show(io::IO, hierarchy::HierarchicalBasisPartition1D)
    print(
        io,
        "HierarchicalBasisPartition1D(nbasis=",
        length(hierarchy.center_data),
        ", nboxes=",
        length(hierarchy.box_data),
        ", nleaves=",
        length(leaf_boxes(hierarchy)),
        ")",
    )
end

boxes(hierarchy::HierarchicalBasisPartition1D) = hierarchy.box_data
leaf_boxes(hierarchy::HierarchicalBasisPartition1D) =
    [box for box in hierarchy.box_data if isempty(box.child_indices)]
box_indices(hierarchy::HierarchicalBasisPartition1D, box_index::Integer) = hierarchy.box_data[box_index].basis_indices
box_level(hierarchy::HierarchicalBasisPartition1D, box_index::Integer) = hierarchy.box_data[box_index].level
box_parent(hierarchy::HierarchicalBasisPartition1D, box_index::Integer) = hierarchy.box_data[box_index].parent_index
box_children(hierarchy::HierarchicalBasisPartition1D, box_index::Integer) = hierarchy.box_data[box_index].child_indices

function _hierarchical_box_assignment(center_value::Float64, edges::Vector{Float64})
    return _box_assignment(center_value, edges)
end

function _build_hierarchical_box(
    centers_data::Vector{Float64},
    basis_indices_data::Vector{Int},
    edges_data::Vector{Float64},
    level_value::Int,
    parent_index_value::Union{Nothing, Int},
    start_index::Int,
)
    assignments = [Int[] for _ in 1:(length(edges_data) - 1)]
    for basis_index in basis_indices_data
        box_index = _hierarchical_box_assignment(centers_data[basis_index], edges_data)
        box_index === nothing &&
            throw(ArgumentError("basis center $(centers_data[basis_index]) is outside the refinement edges"))
        push!(assignments[box_index], basis_index)
    end

    return HierarchicalBasisBox1D[
        HierarchicalBasisBox1D(
            start_index + local_index - 1,
            parent_index_value,
            Int[],
            level_value,
            edges_data[local_index],
            edges_data[local_index + 1],
            assignments[local_index],
        ) for local_index in 1:length(assignments)
    ]
end

"""
    hierarchical_partition(partition)
    hierarchical_partition(basis, edges)
    hierarchical_partition(representation, edges)

Build a small geometric parent-child tree on top of an existing 1D interval
partition.
"""
function hierarchical_partition(partition::BasisPartition1D)
    box_data = HierarchicalBasisBox1D[
        HierarchicalBasisBox1D(box.index, nothing, Int[], 0, box.xmin, box.xmax, copy(box.basis_indices))
        for box in boxes(partition)
    ]
    return HierarchicalBasisPartition1D(copy(partition.center_data), box_data, collect(1:length(box_data)))
end

function hierarchical_partition(
    basis::Union{UniformBasis, HalfLineBasis, RadialBasis},
    edges::AbstractVector{<:Real},
)
    return hierarchical_partition(basis_partition(basis, edges))
end

function hierarchical_partition(
    representation::BasisRepresentation1D,
    edges::AbstractVector{<:Real},
)
    return hierarchical_partition(basis_partition(representation, edges))
end

function _refinement_edges(
    hierarchy::HierarchicalBasisPartition1D,
    box_index::Integer;
    child_edges = nothing,
)
    parent_box = hierarchy.box_data[box_index]
    if child_edges === nothing
        midpoint = 0.5 * (parent_box.xmin + parent_box.xmax)
        return [parent_box.xmin, midpoint, parent_box.xmax]
    end
    edges_data = Float64[Float64(edge) for edge in child_edges]
    length(edges_data) >= 2 || throw(ArgumentError("refinement needs at least two child edges"))
    issorted(edges_data) || throw(ArgumentError("refinement child edges must be sorted"))
    allunique(edges_data) || throw(ArgumentError("refinement child edges must be distinct"))
    edges_data[1] == parent_box.xmin ||
        throw(ArgumentError("refinement child edges must start at the parent xmin"))
    edges_data[end] == parent_box.xmax ||
        throw(ArgumentError("refinement child edges must end at the parent xmax"))
    return edges_data
end

"""
    refine_partition(hierarchy, box_index; child_edges=nothing)

Refine one leaf box of a hierarchical 1D partition.

If `child_edges` is omitted, the box is split at its midpoint. Otherwise the
provided child edges are used.
"""
function refine_partition(
    hierarchy::HierarchicalBasisPartition1D,
    box_index::Integer;
    child_edges = nothing,
)
    1 <= box_index <= length(hierarchy.box_data) ||
        throw(BoundsError(hierarchy.box_data, box_index))
    parent_box = hierarchy.box_data[box_index]
    isempty(parent_box.child_indices) ||
        throw(ArgumentError("only leaf boxes can be refined"))

    edges_data = _refinement_edges(hierarchy, box_index; child_edges = child_edges)
    new_boxes = copy(hierarchy.box_data)
    child_box_data = _build_hierarchical_box(
        hierarchy.center_data,
        parent_box.basis_indices,
        edges_data,
        parent_box.level + 1,
        box_index,
        length(new_boxes) + 1,
    )
    child_ids = [box.index for box in child_box_data]
    new_boxes[box_index] = HierarchicalBasisBox1D(
        parent_box.index,
        parent_box.parent_index,
        child_ids,
        parent_box.level,
        parent_box.xmin,
        parent_box.xmax,
        copy(parent_box.basis_indices),
    )
    append!(new_boxes, child_box_data)
    return HierarchicalBasisPartition1D(copy(hierarchy.center_data), new_boxes, copy(hierarchy.root_indices))
end

function _check_hierarchical_matrix(matrix::AbstractMatrix{<:Real}, hierarchy::HierarchicalBasisPartition1D)
    size(matrix, 1) == size(matrix, 2) ||
        throw(ArgumentError("hierarchical box extraction requires a square basis matrix"))
    size(matrix, 1) == length(hierarchy.center_data) ||
        throw(ArgumentError("matrix size does not match the hierarchical basis size"))
end

function box_block(matrix::AbstractMatrix{<:Real}, hierarchy::HierarchicalBasisPartition1D, box_index::Integer)
    _check_hierarchical_matrix(matrix, hierarchy)
    indices = box_indices(hierarchy, box_index)
    return Matrix{Float64}(matrix[indices, indices])
end

function box_coupling(
    matrix::AbstractMatrix{<:Real},
    hierarchy::HierarchicalBasisPartition1D,
    box_i::Integer,
    box_j::Integer,
)
    _check_hierarchical_matrix(matrix, hierarchy)
    indices_i = box_indices(hierarchy, box_i)
    indices_j = box_indices(hierarchy, box_j)
    return Matrix{Float64}(matrix[indices_i, indices_j])
end

function box_block(
    representation::BasisRepresentation1D,
    hierarchy::HierarchicalBasisPartition1D,
    operator_name::Symbol,
    box_index::Integer,
)
    return box_block(_representation_basis_matrix(representation, operator_name), hierarchy, box_index)
end

function box_coupling(
    representation::BasisRepresentation1D,
    hierarchy::HierarchicalBasisPartition1D,
    operator_name::Symbol,
    box_i::Integer,
    box_j::Integer,
)
    return box_coupling(
        _representation_basis_matrix(representation, operator_name),
        hierarchy,
        box_i,
        box_j,
    )
end
