"""
    LeafLocalPGDG1D

Small 1D hierarchy-driven local Gaussian construction.
"""
struct LeafLocalPGDG1D
    hierarchy::HierarchicalBasisPartition1D
    leaf_box_ids::Vector{Int}
    leaf_primitive_map::Dict{Int, UnitRange{Int}}
    primitive_layer::PrimitiveSet1D
    coefficient_matrix::Matrix{Float64}
    metadata::BasisMetadata1D{Nothing}
    primitives_per_leaf::Int
    width_scale::Float64
end

function Base.show(io::IO, generator::LeafLocalPGDG1D)
    print(
        io,
        "LeafLocalPGDG1D(nleaves=",
        length(generator.leaf_box_ids),
        ", nbasis=",
        size(generator.coefficient_matrix, 2),
        ", primitives_per_leaf=",
        generator.primitives_per_leaf,
        ")",
    )
end

primitive_set(generator::LeafLocalPGDG1D) = generator.primitive_layer
stencil_matrix(generator::LeafLocalPGDG1D) = generator.coefficient_matrix
basis_metadata(generator::LeafLocalPGDG1D) = generator.metadata
leaf_primitive_indices(generator::LeafLocalPGDG1D, leaf_box_index::Integer) = generator.leaf_primitive_map[leaf_box_index]

function _leaf_local_centers(xmin::Float64, xmax::Float64, count::Int)
    span = xmax - xmin
    span > 0.0 || throw(ArgumentError("leaf PGDG generation requires a positive box width"))
    count >= 1 || throw(ArgumentError("leaf PGDG generation requires at least one primitive per leaf"))
    if count == 1
        return [0.5 * (xmin + xmax)]
    end
    spacing = span / (count + 1)
    return collect((xmin + spacing):spacing:(xmax - spacing))
end

function _leaf_local_width(xmin::Float64, xmax::Float64, count::Int, width_scale::Float64)
    span = xmax - xmin
    spacing = count == 1 ? span : span / (count + 1)
    width = width_scale * spacing
    width > 0.0 || throw(ArgumentError("leaf PGDG generation requires a positive local width"))
    return width
end

function _leaf_pgdg_representation(
    metadata::BasisMetadata1D{Nothing},
    generator::LeafLocalPGDG1D,
    operators,
)
    operator_names = Tuple(Symbol(operator_name) for operator_name in operators)
    length(unique(operator_names)) == length(operator_names) ||
        throw(ArgumentError("basis representation operator names must be unique"))

    primitive_matrix_values =
        Tuple(_representation_operator_matrix(generator.primitive_layer, operator_name) for operator_name in operator_names)
    basis_matrix_values =
        Tuple(contract_primitive_matrix(generator, matrix) for matrix in primitive_matrix_values)

    return BasisRepresentation1D(
        metadata,
        generator.primitive_layer,
        generator.coefficient_matrix,
        NamedTuple{operator_names}(primitive_matrix_values),
        NamedTuple{operator_names}(basis_matrix_values),
    )
end

function basis_representation(
    generator::LeafLocalPGDG1D;
    operators = (:overlap, :position, :kinetic),
)
    return _leaf_pgdg_representation(generator.metadata, generator, operators)
end

"""
    build_leaf_pgdg(hierarchy; primitives_per_leaf=2, width_scale=0.75)

Build a small 1D hierarchy-driven local Gaussian construction from the leaf
boxes of `hierarchy`.

Each leaf contributes a fixed local family of full-line Gaussians whose centers
are evenly spaced inside the leaf interval.
"""
function build_leaf_pgdg(
    hierarchy::HierarchicalBasisPartition1D;
    primitives_per_leaf::Int = 2,
    width_scale::Real = 0.75,
)
    primitives_per_leaf >= 1 ||
        throw(ArgumentError("build_leaf_pgdg requires primitives_per_leaf >= 1"))
    width_scale_value = Float64(width_scale)
    width_scale_value > 0.0 || throw(ArgumentError("build_leaf_pgdg requires width_scale > 0"))

    leaf_data = sort(leaf_boxes(hierarchy); by = box -> box.xmin)
    primitive_data = AbstractPrimitiveFunction1D[]
    primitive_labels = String[]
    center_data = Float64[]
    integral_weight_data = Float64[]
    leaf_primitive_map = Dict{Int, UnitRange{Int}}()
    next_index = 1

    for leaf_box in leaf_data
        local_centers = _leaf_local_centers(leaf_box.xmin, leaf_box.xmax, primitives_per_leaf)
        local_width = _leaf_local_width(leaf_box.xmin, leaf_box.xmax, primitives_per_leaf, width_scale_value)
        start_index = next_index
        for (local_index, local_center) in enumerate(local_centers)
            primitive = Gaussian(center = local_center, width = local_width)
            push!(primitive_data, primitive)
            push!(primitive_labels, "leaf$(leaf_box.index)_mu$(local_index)")
            push!(center_data, center(primitive))
            push!(integral_weight_data, integral_weight(primitive))
            next_index += 1
        end
        leaf_primitive_map[leaf_box.index] = start_index:(next_index - 1)
    end

    primitive_layer = PrimitiveSet1D(primitive_data; name = :leaf_pgdg_1d, labels = primitive_labels)
    coefficient_matrix = Matrix{Float64}(I, length(primitive_layer), length(primitive_layer))
    metadata = BasisMetadata1D(
        :leaf_pgdg_1d,
        nothing,
        nothing,
        center_data,
        copy(center_data),
        integral_weight_data,
        copy(primitive_labels),
        primitive_layer,
        coefficient_matrix,
    )
    return LeafLocalPGDG1D(
        hierarchy,
        [box.index for box in leaf_data],
        leaf_primitive_map,
        primitive_layer,
        coefficient_matrix,
        metadata,
        primitives_per_leaf,
        width_scale_value,
    )
end
