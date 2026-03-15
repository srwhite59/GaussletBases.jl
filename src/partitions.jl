"""
    BasisBox1D

Single interval box in a one-dimensional basis partition.
"""
struct BasisBox1D
    index::Int
    xmin::Float64
    xmax::Float64
    basis_indices::Vector{Int}
end

function Base.show(io::IO, box::BasisBox1D)
    print(
        io,
        "BasisBox1D(index=",
        box.index,
        ", range=(",
        box.xmin,
        ", ",
        box.xmax,
        "), n=",
        length(box.basis_indices),
        ")",
    )
end

"""
    BasisPartition1D

Small 1D grouping layer that partitions basis functions into interval boxes
using their physical-space centers.
"""
struct BasisPartition1D
    center_data::Vector{Float64}
    box_data::Vector{BasisBox1D}
end

function Base.show(io::IO, partition::BasisPartition1D)
    print(
        io,
        "BasisPartition1D(nbasis=",
        length(partition.center_data),
        ", nboxes=",
        length(partition.box_data),
        ")",
    )
end

boxes(partition::BasisPartition1D) = partition.box_data
box_indices(partition::BasisPartition1D, box_index::Integer) = partition.box_data[box_index].basis_indices

function _partition_centers_and_count(basis::Union{UniformBasis, HalfLineBasis, RadialBasis})
    return Float64[Float64(value) for value in centers(basis)], length(basis)
end

function _partition_centers_and_count(representation::BasisRepresentation1D)
    return representation.metadata.center_data, size(representation.coefficient_matrix, 2)
end

function _box_assignment(center_value::Float64, edges::Vector{Float64})
    for box_index in 1:(length(edges) - 2)
        edges[box_index] <= center_value < edges[box_index + 1] && return box_index
    end
    edges[end - 1] <= center_value <= edges[end] && return length(edges) - 1
    return nothing
end

function _build_basis_partition(centers_data::Vector{Float64}, basis_count::Int, edges_data::Vector{Float64})
    length(edges_data) >= 2 || throw(ArgumentError("basis partitions require at least two interval edges"))
    issorted(edges_data) || throw(ArgumentError("basis partition edges must be sorted"))
    allunique(edges_data) || throw(ArgumentError("basis partition edges must be distinct"))

    assignments = [Int[] for _ in 1:(length(edges_data) - 1)]
    for basis_index in 1:basis_count
        box_index = _box_assignment(centers_data[basis_index], edges_data)
        box_index === nothing &&
            throw(ArgumentError("basis center $(centers_data[basis_index]) is outside the partition edges"))
        push!(assignments[box_index], basis_index)
    end

    box_data = BasisBox1D[
        BasisBox1D(box_index, edges_data[box_index], edges_data[box_index + 1], assignments[box_index])
        for box_index in 1:length(assignments)
    ]
    return BasisPartition1D(centers_data, box_data)
end

"""
    basis_partition(basis, edges)
    basis_partition(representation, edges)

Partition basis functions into 1D interval boxes using their physical-space
centers.

The intervals are half-open on the right except for the final box, which
includes its upper edge.
"""
function basis_partition(
    basis::Union{UniformBasis, HalfLineBasis, RadialBasis},
    edges::AbstractVector{<:Real},
)
    centers_data, basis_count = _partition_centers_and_count(basis)
    return _build_basis_partition(centers_data, basis_count, Float64[Float64(edge) for edge in edges])
end

function basis_partition(
    representation::BasisRepresentation1D,
    edges::AbstractVector{<:Real},
)
    centers_data, basis_count = _partition_centers_and_count(representation)
    return _build_basis_partition(centers_data, basis_count, Float64[Float64(edge) for edge in edges])
end

function _check_partition_matrix(matrix::AbstractMatrix{<:Real}, partition::BasisPartition1D)
    size(matrix, 1) == size(matrix, 2) ||
        throw(ArgumentError("box extraction requires a square basis matrix"))
    size(matrix, 1) == length(partition.center_data) ||
        throw(ArgumentError("matrix size does not match the partition basis size"))
end

function box_block(matrix::AbstractMatrix{<:Real}, partition::BasisPartition1D, box_index::Integer)
    _check_partition_matrix(matrix, partition)
    indices = box_indices(partition, box_index)
    return Matrix{Float64}(matrix[indices, indices])
end

function box_coupling(
    matrix::AbstractMatrix{<:Real},
    partition::BasisPartition1D,
    box_i::Integer,
    box_j::Integer,
)
    _check_partition_matrix(matrix, partition)
    indices_i = box_indices(partition, box_i)
    indices_j = box_indices(partition, box_j)
    return Matrix{Float64}(matrix[indices_i, indices_j])
end

function _representation_basis_matrix(representation::BasisRepresentation1D, operator_name::Symbol)
    hasproperty(representation.basis_matrices, operator_name) ||
        throw(ArgumentError("basis representation does not contain operator :$operator_name"))
    return getproperty(representation.basis_matrices, operator_name)
end

function box_block(
    representation::BasisRepresentation1D,
    partition::BasisPartition1D,
    operator_name::Symbol,
    box_index::Integer,
)
    return box_block(_representation_basis_matrix(representation, operator_name), partition, box_index)
end

function box_coupling(
    representation::BasisRepresentation1D,
    partition::BasisPartition1D,
    operator_name::Symbol,
    box_i::Integer,
    box_j::Integer,
)
    return box_coupling(
        _representation_basis_matrix(representation, operator_name),
        partition,
        box_i,
        box_j,
    )
end
