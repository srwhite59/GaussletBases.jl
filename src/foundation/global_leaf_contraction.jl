"""
    GlobalMappedPrimitiveLayer1D

Common 1D globally mapped primitive layer over a region.
"""
struct GlobalMappedPrimitiveLayer1D{M <: AbstractCoordinateMapping}
    xmin::Float64
    xmax::Float64
    reference_spacing::Float64
    width_scale::Float64
    mapping::M
    primitive_layer::PrimitiveSet1D
    coefficient_matrix::Matrix{Float64}
    metadata::BasisMetadata1D{M}
end

function Base.show(io::IO, layer::GlobalMappedPrimitiveLayer1D)
    print(
        io,
        "GlobalMappedPrimitiveLayer1D(nbasis=",
        size(layer.coefficient_matrix, 2),
        ", xmin=",
        layer.xmin,
        ", xmax=",
        layer.xmax,
        ", reference_spacing=",
        layer.reference_spacing,
        ")",
    )
end

primitive_set(layer::GlobalMappedPrimitiveLayer1D) = layer.primitive_layer
stencil_matrix(layer::GlobalMappedPrimitiveLayer1D) = layer.coefficient_matrix
basis_metadata(layer::GlobalMappedPrimitiveLayer1D) = layer.metadata

function basis_representation(
    layer::GlobalMappedPrimitiveLayer1D;
    operators = (:overlap, :position, :kinetic),
)
    operator_names = Tuple(Symbol(operator_name) for operator_name in operators)
    primitive_matrix_values =
        Tuple(_representation_operator_matrix(layer.primitive_layer, operator_name) for operator_name in operator_names)
    basis_matrix_values =
        Tuple(contract_primitive_matrix(layer, matrix) for matrix in primitive_matrix_values)

    return BasisRepresentation1D(
        layer.metadata,
        layer.primitive_layer,
        layer.coefficient_matrix,
        NamedTuple{operator_names}(primitive_matrix_values),
        NamedTuple{operator_names}(basis_matrix_values),
    )
end

function basis_partition(
    layer::GlobalMappedPrimitiveLayer1D,
    edges::AbstractVector{<:Real},
)
    return _build_basis_partition(layer.metadata.center_data, size(layer.coefficient_matrix, 2), Float64[Float64(edge) for edge in edges])
end

function hierarchical_partition(
    layer::GlobalMappedPrimitiveLayer1D,
    edges::AbstractVector{<:Real},
)
    return hierarchical_partition(basis_partition(layer, edges))
end

"""
    build_global_mapped_primitive_layer(; xmin, xmax, mapping, reference_spacing=1.0, width_scale=1.0)

Build one common 1D primitive layer over `[xmin, xmax]` using a single global
coordinate map.
"""
function build_global_mapped_primitive_layer(;
    xmin::Real,
    xmax::Real,
    mapping::AbstractCoordinateMapping,
    reference_spacing::Real = 1.0,
    width_scale::Real = 1.0,
)
    xmin_value = Float64(xmin)
    xmax_value = Float64(xmax)
    xmin_value < xmax_value || throw(ArgumentError("global mapped primitive layer requires xmin < xmax"))
    spacing_value = Float64(reference_spacing)
    spacing_value > 0.0 || throw(ArgumentError("global mapped primitive layer requires reference_spacing > 0"))
    width_scale_value = Float64(width_scale)
    width_scale_value > 0.0 || throw(ArgumentError("global mapped primitive layer requires width_scale > 0"))

    umin = uofx(mapping, xmin_value)
    umax = uofx(mapping, xmax_value)
    reference_centers = _uniform_centers(umin, umax, spacing_value)
    primitive_width = width_scale_value * spacing_value
    primitive_ref = AbstractPrimitiveFunction1D[
        Distorted(Gaussian(center = ucenter, width = primitive_width), mapping)
        for ucenter in reference_centers
    ]
    primitive_layer = PrimitiveSet1D(primitive_ref; name = :global_mapped_primitive_1d)
    coefficient_matrix = Matrix{Float64}(I, length(primitive_layer), length(primitive_layer))
    center_data = Float64[xofu(mapping, value) for value in reference_centers]
    integral_weight_data = Float64[integral_weight(primitive) for primitive in primitive_ref]
    labels = ["chi$(i)" for i in 1:length(primitive_layer)]
    metadata = BasisMetadata1D(
        :global_mapped_primitive_1d,
        nothing,
        mapping,
        center_data,
        copy(reference_centers),
        integral_weight_data,
        labels,
        primitive_layer,
        coefficient_matrix,
    )
    return GlobalMappedPrimitiveLayer1D(
        xmin_value,
        xmax_value,
        spacing_value,
        width_scale_value,
        mapping,
        primitive_layer,
        coefficient_matrix,
        metadata,
    )
end

"""
    LeafBoxContraction1D

Single leaf-local contraction built from a globally defined primitive layer.
"""
struct LeafBoxContraction1D
    leaf_box_index::Int
    primitive_indices::Vector{Int}
    local_coefficient_matrix::Matrix{Float64}
    retained_centers::Vector{Float64}
end

function Base.show(io::IO, contraction::LeafBoxContraction1D)
    print(
        io,
        "LeafBoxContraction1D(leaf=",
        contraction.leaf_box_index,
        ", nprimitive=",
        length(contraction.primitive_indices),
        ", nretained=",
        size(contraction.local_coefficient_matrix, 2),
        ")",
    )
end

"""
    LeafBoxContractionLayer1D

Optional leaf-local contraction layer built from one common globally mapped
primitive layer.
"""
struct LeafBoxContractionLayer1D{M <: AbstractCoordinateMapping}
    global_layer::GlobalMappedPrimitiveLayer1D{M}
    hierarchy::HierarchicalBasisPartition1D
    contraction_data::Vector{LeafBoxContraction1D}
    coefficient_matrix::Matrix{Float64}
    metadata::BasisMetadata1D{M}
    retained_per_leaf::Int
end

function Base.show(io::IO, layer::LeafBoxContractionLayer1D)
    print(
        io,
        "LeafBoxContractionLayer1D(nleaves=",
        length(layer.contraction_data),
        ", nbasis=",
        size(layer.coefficient_matrix, 2),
        ", retained_per_leaf=",
        layer.retained_per_leaf,
        ")",
    )
end

primitive_set(layer::LeafBoxContractionLayer1D) = primitive_set(layer.global_layer)
stencil_matrix(layer::LeafBoxContractionLayer1D) = layer.coefficient_matrix
basis_metadata(layer::LeafBoxContractionLayer1D) = layer.metadata
leaf_contractions(layer::LeafBoxContractionLayer1D) = layer.contraction_data

function basis_representation(
    layer::LeafBoxContractionLayer1D;
    operators = (:overlap, :position, :kinetic),
)
    operator_names = Tuple(Symbol(operator_name) for operator_name in operators)
    primitive_matrix_values =
        Tuple(_representation_operator_matrix(primitive_set(layer), operator_name) for operator_name in operator_names)
    basis_matrix_values =
        Tuple(contract_primitive_matrix(layer, matrix) for matrix in primitive_matrix_values)

    return BasisRepresentation1D(
        layer.metadata,
        primitive_set(layer),
        layer.coefficient_matrix,
        NamedTuple{operator_names}(primitive_matrix_values),
        NamedTuple{operator_names}(basis_matrix_values),
    )
end

function basis_partition(
    layer::LeafBoxContractionLayer1D,
    edges::AbstractVector{<:Real},
)
    return _build_basis_partition(layer.metadata.center_data, size(layer.coefficient_matrix, 2), Float64[Float64(edge) for edge in edges])
end

function hierarchical_partition(
    layer::LeafBoxContractionLayer1D,
    edges::AbstractVector{<:Real},
)
    return hierarchical_partition(basis_partition(layer, edges))
end

function _local_orthonormal_coefficients(overlap::Matrix{Float64})
    eigenvalues, eigenvectors = eigen(Symmetric(overlap))
    maximum(eigenvalues) > 0.0 || throw(ArgumentError("local contraction requires a positive local overlap spectrum"))
    keep = findall(value -> value > 1.0e-10 * maximum(eigenvalues), eigenvalues)
    isempty(keep) && throw(ArgumentError("local contraction found no linearly independent local primitive directions"))
    return eigenvectors[:, keep] * Diagonal(1.0 ./ sqrt.(eigenvalues[keep]))
end

function _selected_local_modes(position_values::Vector{Float64}, midpoint::Float64, retained_per_leaf::Int)
    retained_per_leaf >= length(position_values) && return collect(1:length(position_values))
    ordering = sortperm(abs.(position_values .- midpoint))
    selected = ordering[1:retained_per_leaf]
    return sort(selected; by = index -> position_values[index])
end

"""
    contract_leaf_boxes(global_layer, hierarchy; retained_per_leaf=1)

Build a leaf-local contraction layer from one common globally mapped primitive
layer and an existing hierarchy.
"""
function contract_leaf_boxes(
    global_layer::GlobalMappedPrimitiveLayer1D,
    hierarchy::HierarchicalBasisPartition1D;
    retained_per_leaf::Int = 1,
)
    retained_per_leaf >= 1 || throw(ArgumentError("contract_leaf_boxes requires retained_per_leaf >= 1"))
    length(global_layer.metadata.center_data) == length(hierarchy.center_data) ||
        throw(ArgumentError("leaf contraction requires a hierarchy defined on the global layer centers"))

    global_representation = basis_representation(global_layer)
    overlap_global = global_representation.basis_matrices.overlap
    position_global = global_representation.basis_matrices.position
    primitive_weights = global_layer.metadata.integral_weight_data

    local_contractions = LeafBoxContraction1D[]
    retained_centers = Float64[]
    retained_reference_centers = Float64[]
    basis_labels = String[]
    total_columns = 0

    for leaf in sort(leaf_boxes(hierarchy); by = box -> box.xmin)
        primitive_indices = copy(leaf.basis_indices)
        overlap_local = overlap_global[primitive_indices, primitive_indices]
        position_local = position_global[primitive_indices, primitive_indices]
        orthonormal_coefficients = _local_orthonormal_coefficients(overlap_local)
        position_orthonormal = Symmetric(transpose(orthonormal_coefficients) * position_local * orthonormal_coefficients)
        local_positions, local_vectors = eigen(position_orthonormal)
        local_coefficients = orthonormal_coefficients * local_vectors
        selected_modes = _selected_local_modes(local_positions, 0.5 * (leaf.xmin + leaf.xmax), retained_per_leaf)
        retained_local_coefficients = local_coefficients[:, selected_modes]
        retained_local_centers = collect(local_positions[selected_modes])

        push!(
            local_contractions,
            LeafBoxContraction1D(
                leaf.index,
                primitive_indices,
                retained_local_coefficients,
                retained_local_centers,
            ),
        )

        append!(retained_centers, retained_local_centers)
        append!(retained_reference_centers, [uofx(global_layer.mapping, value) for value in retained_local_centers])
        append!(basis_labels, ["leaf$(leaf.index)_phi$(i)" for i in 1:length(retained_local_centers)])
        total_columns += length(retained_local_centers)
    end

    coefficient_matrix = zeros(Float64, length(primitive_set(global_layer)), total_columns)
    column_offset = 0
    for contraction in local_contractions
        nlocal = size(contraction.local_coefficient_matrix, 2)
        coefficient_matrix[contraction.primitive_indices, (column_offset + 1):(column_offset + nlocal)] .=
            contraction.local_coefficient_matrix
        column_offset += nlocal
    end

    integral_weight_data = vec(transpose(coefficient_matrix) * primitive_weights)
    metadata = BasisMetadata1D(
        :leaf_box_contraction_1d,
        nothing,
        global_layer.mapping,
        retained_centers,
        retained_reference_centers,
        integral_weight_data,
        basis_labels,
        primitive_set(global_layer),
        coefficient_matrix,
    )
    return LeafBoxContractionLayer1D(
        global_layer,
        hierarchy,
        local_contractions,
        coefficient_matrix,
        metadata,
        retained_per_leaf,
    )
end
