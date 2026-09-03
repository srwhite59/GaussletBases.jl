"""
    LeafGaussianSpec1D(; relative_position, width_scale)

Small user-supplied local Gaussian augmentation spec for one leaf interval.

`relative_position` is measured inside the leaf interval, from `0.0` at the
left edge to `1.0` at the right edge.
"""
struct LeafGaussianSpec1D
    relative_position::Float64
    width_scale::Float64

    function LeafGaussianSpec1D(; relative_position::Real, width_scale::Real)
        position_value = Float64(relative_position)
        width_value = Float64(width_scale)
        0.0 <= position_value <= 1.0 ||
            throw(ArgumentError("LeafGaussianSpec1D requires 0 <= relative_position <= 1"))
        width_value > 0.0 ||
            throw(ArgumentError("LeafGaussianSpec1D requires width_scale > 0"))
        new(position_value, width_value)
    end
end

function Base.show(io::IO, spec::LeafGaussianSpec1D)
    print(
        io,
        "LeafGaussianSpec1D(relative_position=",
        spec.relative_position,
        ", width_scale=",
        spec.width_scale,
        ")",
    )
end

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
    primitive_origin_data::Vector{Symbol}
    primitive_leaf_box_data::Vector{Int}
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
        ", naugmented=",
        count(==(:augmented), generator.primitive_origin_data),
        ")",
    )
end

primitive_set(generator::LeafLocalPGDG1D) = generator.primitive_layer
stencil_matrix(generator::LeafLocalPGDG1D) = generator.coefficient_matrix
basis_metadata(generator::LeafLocalPGDG1D) = generator.metadata
leaf_primitive_indices(generator::LeafLocalPGDG1D, leaf_box_index::Integer) = generator.leaf_primitive_map[leaf_box_index]
primitive_origins(generator::LeafLocalPGDG1D) = generator.primitive_origin_data
primitive_leaf_boxes(generator::LeafLocalPGDG1D) = generator.primitive_leaf_box_data

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

function _leaf_local_augmented_gaussian(
    leaf_box::HierarchicalBasisBox1D,
    spec::LeafGaussianSpec1D,
)
    span = leaf_box.xmax - leaf_box.xmin
    center_value = leaf_box.xmin + spec.relative_position * span
    width_value = spec.width_scale * span
    return Gaussian(center = center_value, width = width_value)
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

function _build_leaf_pgdg(
    hierarchy::HierarchicalBasisPartition1D,
    primitives_per_leaf::Int,
    width_scale_value::Float64,
    every_leaf_specs::Vector{LeafGaussianSpec1D},
    by_leaf_specs::Dict{Int, Vector{LeafGaussianSpec1D}},
)
    leaf_data = sort(leaf_boxes(hierarchy); by = box -> box.xmin)
    primitive_data = AbstractPrimitiveFunction1D[]
    primitive_labels = String[]
    center_data = Float64[]
    integral_weight_data = Float64[]
    primitive_origin_data = Symbol[]
    primitive_leaf_box_data = Int[]
    leaf_primitive_map = Dict{Int, UnitRange{Int}}()
    next_index = 1

    for leaf_box in leaf_data
        start_index = next_index

        local_centers = _leaf_local_centers(leaf_box.xmin, leaf_box.xmax, primitives_per_leaf)
        local_width = _leaf_local_width(leaf_box.xmin, leaf_box.xmax, primitives_per_leaf, width_scale_value)
        for (local_index, local_center) in enumerate(local_centers)
            primitive = Gaussian(center = local_center, width = local_width)
            push!(primitive_data, primitive)
            push!(primitive_labels, "leaf$(leaf_box.index)_mu$(local_index)")
            push!(center_data, center(primitive))
            push!(integral_weight_data, integral_weight(primitive))
            push!(primitive_origin_data, :generated)
            push!(primitive_leaf_box_data, leaf_box.index)
            next_index += 1
        end

        local_specs = vcat(every_leaf_specs, get(by_leaf_specs, leaf_box.index, LeafGaussianSpec1D[]))
        for (aug_index, spec) in enumerate(local_specs)
            primitive = _leaf_local_augmented_gaussian(leaf_box, spec)
            push!(primitive_data, primitive)
            push!(primitive_labels, "leaf$(leaf_box.index)_aug$(aug_index)")
            push!(center_data, center(primitive))
            push!(integral_weight_data, integral_weight(primitive))
            push!(primitive_origin_data, :augmented)
            push!(primitive_leaf_box_data, leaf_box.index)
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
        primitive_origin_data,
        primitive_leaf_box_data,
        primitives_per_leaf,
        width_scale_value,
    )
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
    return _build_leaf_pgdg(
        hierarchy,
        primitives_per_leaf,
        width_scale_value,
        LeafGaussianSpec1D[],
        Dict{Int, Vector{LeafGaussianSpec1D}}(),
    )
end

"""
    augment_leaf_pgdg(generator; by_leaf=Dict{Int,Vector{LeafGaussianSpec1D}}(), every_leaf=LeafGaussianSpec1D[])

Add optional user-supplied Gaussian augmentation to an existing
`LeafLocalPGDG1D`.

`every_leaf` applies the same local augmentation specs to every leaf. `by_leaf`
adds leaf-specific specs keyed by hierarchical leaf-box index.
"""
function augment_leaf_pgdg(
    generator::LeafLocalPGDG1D;
    by_leaf = Dict{Int, Vector{LeafGaussianSpec1D}}(),
    every_leaf = LeafGaussianSpec1D[],
)
    by_leaf_specs = Dict{Int, Vector{LeafGaussianSpec1D}}(
        Int(box_index) => LeafGaussianSpec1D[spec for spec in specs]
        for (box_index, specs) in pairs(by_leaf)
    )
    every_leaf_specs = LeafGaussianSpec1D[spec for spec in every_leaf]

    for box_index in keys(by_leaf_specs)
        box_index in generator.leaf_box_ids ||
            throw(ArgumentError("augmentation leaf $box_index is not a leaf of the generator hierarchy"))
    end

    return _build_leaf_pgdg(
        generator.hierarchy,
        generator.primitives_per_leaf,
        generator.width_scale,
        every_leaf_specs,
        by_leaf_specs,
    )
end
