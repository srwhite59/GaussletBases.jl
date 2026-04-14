function _cartesian_bundle_supported_basis(
    representation::CartesianBasisRepresentation3D,
)
    representation.metadata.parent_kind == :cartesian_product_basis || throw(
        ArgumentError(
            "cartesian basis bundle export does not yet support hybrid/QW residual representations because this first bundle pass only covers pure Cartesian raw-space contracts",
        ),
    )
    return representation
end

function _cartesian_bundle_representation(
    representation::CartesianBasisRepresentation3D,
)
    return _cartesian_bundle_supported_basis(representation)
end

function _cartesian_bundle_representation(
    basis::Union{
        BondAlignedDiatomicQWBasis3D,
        BondAlignedHomonuclearChainQWBasis3D,
        AxisAlignedHomonuclearSquareLatticeQWBasis3D,
        _NestedFixedBlock3D,
    },
)
    return _cartesian_bundle_supported_basis(basis_representation(basis))
end

function _cartesian_bundle_representation(operators::OrdinaryCartesianIDAOperators)
    return _cartesian_bundle_supported_basis(basis_representation(operators.basis))
end

function _cartesian_qw_bundle_is_pure(
    operators::QiuWhiteResidualGaussianOperators,
)
    return operators.gaussian_data === nothing && operators.residual_count == 0
end

function _cartesian_bundle_representation(
    operators::QiuWhiteResidualGaussianOperators,
)
    _cartesian_qw_bundle_is_pure(operators) || throw(
        ArgumentError(
            "cartesian basis bundle export does not yet support hybrid/QW residual operators because the current first-pass bundle contract only covers pure Cartesian basis/raw-space families",
        ),
    )
    return _cartesian_bundle_representation(operators.basis)
end

function _cartesian_store_value!(
    dest::Dict{String,Any},
    prefix::AbstractString,
    value,
)
    if value === nothing
        dest[string(prefix, "/is_nothing")] = true
    elseif value isa NamedTuple
        for (key, child) in pairs(value)
            _cartesian_store_value!(dest, string(prefix, "/", key), child)
        end
    elseif value isa AbstractDict
        for (key, child) in pairs(value)
            _cartesian_store_value!(dest, string(prefix, "/", key), child)
        end
    elseif value isa Symbol
        dest[prefix] = String(value)
    elseif value isa UnitRange{<:Integer}
        dest[prefix] = Int[first(value), last(value)]
    elseif value isa Tuple && all(item -> item isa Symbol, value)
        dest[prefix] = String[String(item) for item in value]
    elseif value isa Tuple && all(item -> item isa Integer, value)
        dest[prefix] = Int[Int(item) for item in value]
    elseif value isa AbstractArray{<:Symbol}
        dest[prefix] = String[String(item) for item in value]
    elseif value isa AbstractArray
        dest[prefix] = copy(value)
    elseif value isa Number || value isa AbstractString || value isa Bool
        dest[prefix] = value
    else
        dest[prefix] = value
    end
    return dest
end

function _cartesian_store_axis_representation!(
    dest::Dict{String,Any},
    prefix::AbstractString,
    representation::BasisRepresentation1D,
)
    metadata = representation.metadata
    primitive_layer = representation.primitive_set

    dest[string(prefix, "/format")] = "basis_representation_1d_v1"
    dest[string(prefix, "/version")] = 1
    dest[string(prefix, "/metadata/basis_kind")] = String(metadata.basis_kind)
    dest[string(prefix, "/metadata/has_family_name")] = metadata.family_name !== nothing
    metadata.family_name !== nothing &&
        (dest[string(prefix, "/metadata/family_name")] = String(metadata.family_name))
    dest[string(prefix, "/metadata/mapping_type")] = string(nameof(typeof(metadata.mapping_value)))
    dest[string(prefix, "/metadata/mapping_object")] = metadata.mapping_value
    dest[string(prefix, "/metadata/centers")] = copy(metadata.center_data)
    dest[string(prefix, "/metadata/reference_centers")] = copy(metadata.reference_center_data)
    dest[string(prefix, "/metadata/integral_weights")] = copy(metadata.integral_weight_data)
    dest[string(prefix, "/metadata/basis_labels")] = copy(metadata.basis_labels)
    dest[string(prefix, "/primitive_set/has_name")] = primitive_layer.name_value !== nothing
    primitive_layer.name_value !== nothing &&
        (dest[string(prefix, "/primitive_set/name")] = String(primitive_layer.name_value))
    dest[string(prefix, "/primitive_set/labels")] = copy(primitive_layer.label_data)
    dest[string(prefix, "/primitive_set/primitives")] =
        AbstractPrimitiveFunction1D[primitive for primitive in primitives(primitive_layer)]
    dest[string(prefix, "/coefficient_matrix")] = Matrix{Float64}(representation.coefficient_matrix)
    return dest
end

function _cartesian_parent_integral_weights(
    representation::CartesianBasisRepresentation3D,
)
    weights_x = representation.axis_representations.x.metadata.integral_weight_data
    weights_y = representation.axis_representations.y.metadata.integral_weight_data
    weights_z = representation.axis_representations.z.metadata.integral_weight_data
    return _mapped_cartesian_weights(weights_x, weights_y, weights_z)
end

function _cartesian_representation_integral_weights(
    representation::CartesianBasisRepresentation3D,
)
    _cartesian_bundle_supported_basis(representation)
    parent_weights = _cartesian_parent_integral_weights(representation)

    if representation.contraction_kind == :identity
        return Vector{Float64}(parent_weights)
    elseif representation.contraction_kind == :dense
        support_weights =
            representation.support_indices === nothing ?
            parent_weights :
            parent_weights[representation.support_indices]
        return Vector{Float64}(transpose(representation.coefficient_matrix) * support_weights)
    end

    throw(
        ArgumentError(
            "cartesian basis bundle export does not yet support contraction kind :$(representation.contraction_kind)",
        ),
    )
end

function _cartesian_bundle_integral_weights(
    representation::CartesianBasisRepresentation3D,
)
    return _cartesian_representation_integral_weights(representation)
end

function _cartesian_bundle_integral_weights(
    fixed_block::_NestedFixedBlock3D,
)
    return Vector{Float64}(fixed_block.weights)
end

function _cartesian_bundle_integral_weights(
    basis::Union{
        BondAlignedDiatomicQWBasis3D,
        BondAlignedHomonuclearChainQWBasis3D,
        AxisAlignedHomonuclearSquareLatticeQWBasis3D,
    },
)
    return _cartesian_representation_integral_weights(basis_representation(basis))
end

function _cartesian_bundle_integral_weights(
    operators::OrdinaryCartesianIDAOperators,
)
    return Vector{Float64}(operators.weight_3d)
end

function _cartesian_bundle_integral_weights(
    operators::QiuWhiteResidualGaussianOperators,
)
    _cartesian_qw_bundle_is_pure(operators) || throw(
        ArgumentError(
            "cartesian basis bundle export does not yet support hybrid/QW residual operators because the current first-pass bundle contract only covers pure Cartesian basis/raw-space families",
        ),
    )
    return _cartesian_bundle_integral_weights(operators.basis)
end

function _cartesian_support_state_matrix(
    states::Union{Nothing,AbstractVector{<:NTuple{3,Int}}},
)
    states === nothing && return zeros(Int, 0, 3)
    matrix = Matrix{Int}(undef, length(states), 3)
    for (row, state) in pairs(states)
        matrix[row, 1] = state[1]
        matrix[row, 2] = state[2]
        matrix[row, 3] = state[3]
    end
    return matrix
end

function _cartesian_working_box_matrix(
    working_box::Union{Nothing,NTuple{3,UnitRange{Int}}},
)
    working_box === nothing && return zeros(Int, 0, 2)
    matrix = Matrix{Int}(undef, 3, 2)
    for axis in 1:3
        matrix[axis, 1] = first(working_box[axis])
        matrix[axis, 2] = last(working_box[axis])
    end
    return matrix
end

function _cartesian_basis_values(
    representation::CartesianBasisRepresentation3D;
    final_integral_weights::AbstractVector{<:Real},
)
    basis_values = Dict{String,Any}(
        "format" => "cartesian_basis_bundle_v1",
        "version" => 1,
        "basis_kind" => String(representation.metadata.basis_kind),
        "parent_kind" => String(representation.metadata.parent_kind),
        "axis_sharing" => String(representation.metadata.axis_sharing),
        "parent_axis_counts" => Int[representation.metadata.parent_axis_counts...],
        "parent_dimension" => representation.metadata.parent_dimension,
        "final_dimension" => representation.metadata.final_dimension,
        "working_box_present" => representation.metadata.working_box !== nothing,
        "working_box_bounds" => _cartesian_working_box_matrix(representation.metadata.working_box),
        "contraction_kind" => String(representation.contraction_kind),
        "basis_labels" => copy(representation.metadata.basis_labels),
        "basis_centers" => Matrix{Float64}(representation.metadata.basis_centers),
        "final_integral_weights" => Float64[Float64(value) for value in final_integral_weights],
        "parent_labels" => copy(representation.parent_labels),
        "parent_centers" => Matrix{Float64}(representation.parent_centers),
        "support_indices_present" => representation.support_indices !== nothing,
        "support_indices" =>
            representation.support_indices === nothing ? zeros(Int, 0) : Vector{Int}(representation.support_indices),
        "support_states_present" => representation.support_states !== nothing,
        "support_states" => _cartesian_support_state_matrix(representation.support_states),
    )

    representation.coefficient_matrix !== nothing &&
        (basis_values["coefficient_matrix"] = Matrix{Float64}(representation.coefficient_matrix))

    _cartesian_store_value!(basis_values, "metadata/route", representation.metadata.route_metadata)
    _cartesian_store_axis_representation!(basis_values, "axes/x", representation.axis_representations.x)
    _cartesian_store_axis_representation!(basis_values, "axes/y", representation.axis_representations.y)
    _cartesian_store_axis_representation!(basis_values, "axes/z", representation.axis_representations.z)

    return basis_values
end

function _cartesian_ham_values(
    ::CartesianBasisRepresentation3D,
    representation::CartesianBasisRepresentation3D,
)
    _cartesian_bundle_supported_basis(representation)
    return nothing
end

function _cartesian_ham_values(
    basis::Union{
        BondAlignedDiatomicQWBasis3D,
        BondAlignedHomonuclearChainQWBasis3D,
        AxisAlignedHomonuclearSquareLatticeQWBasis3D,
        _NestedFixedBlock3D,
    },
    representation::CartesianBasisRepresentation3D,
)
    _cartesian_bundle_supported_basis(representation)
    return nothing
end

function _cartesian_ham_values(
    operators::OrdinaryCartesianIDAOperators,
    representation::CartesianBasisRepresentation3D,
)
    return Dict{String,Any}(
        "format" => "cartesian_hamiltonian_bundle_v1",
        "version" => 1,
        "object_type" => string(nameof(typeof(operators))),
        "model_kind" => "ordinary_cartesian_ida",
        "interaction_model" => "density_density",
        "interaction_treatment" => String(operators.interaction_treatment),
        "backend" => String(operators.backend),
        "overlap_key" => "overlap",
        "onebody_key" => "one_body_hamiltonian",
        "interaction_key" => "interaction_matrix",
        "basis_integral_weights_key" => "basis/final_integral_weights",
        "overlap" => Matrix{Float64}(operators.overlap_3d),
        "one_body_hamiltonian" => Matrix{Float64}(operators.one_body_hamiltonian),
        "interaction_matrix" => Matrix{Float64}(operators.interaction_matrix),
        "orbital_labels" => String[String(orbital.label) for orbital in orbitals(operators)],
        "basis_centers" => Matrix{Float64}(representation.metadata.basis_centers),
        "basis_integral_weights" => Vector{Float64}(operators.weight_3d),
        "expansion/exponents" => copy(operators.expansion.exponents),
        "expansion/coefficients" => copy(operators.expansion.coefficients),
    )
end

function _cartesian_ham_values(
    operators::QiuWhiteResidualGaussianOperators,
    representation::CartesianBasisRepresentation3D,
)
    _cartesian_qw_bundle_is_pure(operators) || throw(
        ArgumentError(
            "cartesian basis bundle export does not yet support hybrid/QW residual operators because the current first-pass bundle contract only covers pure Cartesian basis/raw-space families",
        ),
    )
    return Dict{String,Any}(
        "format" => "cartesian_hamiltonian_bundle_v1",
        "version" => 1,
        "object_type" => string(nameof(typeof(operators))),
        "model_kind" => "ordinary_cartesian_qiu_white",
        "interaction_model" => "density_density",
        "interaction_treatment" => String(operators.interaction_treatment),
        "gausslet_backend" => String(operators.gausslet_backend),
        "overlap_key" => "overlap",
        "onebody_key" => "one_body_hamiltonian",
        "interaction_key" => "interaction_matrix",
        "basis_integral_weights_key" => "basis/final_integral_weights",
        "overlap" => Matrix{Float64}(operators.overlap),
        "one_body_hamiltonian" => Matrix{Float64}(operators.one_body_hamiltonian),
        "interaction_matrix" => Matrix{Float64}(operators.interaction_matrix),
        "orbital_labels" => String[String(orbital.label) for orbital in orbitals(operators)],
        "basis_centers" => Matrix{Float64}(representation.metadata.basis_centers),
        "basis_integral_weights" => _cartesian_bundle_integral_weights(operators),
        "expansion/exponents" => copy(operators.expansion.exponents),
        "expansion/coefficients" => copy(operators.expansion.coefficients),
        "residual_count" => operators.residual_count,
    )
end

function _cartesian_bundle_meta_values(
    object,
    representation::CartesianBasisRepresentation3D;
    include_ham::Bool,
    producer_entrypoint::AbstractString,
    meta = nothing,
)
    meta_values = _normalize_meta_dict(meta)
    merge!(
        meta_values,
        Dict{String,Any}(
            "producer" => producer_entrypoint,
            "producer_type" => string(nameof(typeof(object))),
            "basis_kind" => String(representation.metadata.basis_kind),
            "parent_kind" => String(representation.metadata.parent_kind),
            "has_ham" => include_ham,
            "final_dimension" => representation.metadata.final_dimension,
            "parent_dimension" => representation.metadata.parent_dimension,
            "manifest/producer/package" => "GaussletBases",
            "manifest/producer/version" => string(Base.pkgversion(@__MODULE__)),
            "manifest/producer/entrypoint" => producer_entrypoint,
            "manifest/contract/format" => "cartesian_basis_bundle_v1",
            "manifest/contract/version" => 1,
            "manifest/contract/status" => "public_first_pass",
            "manifest/contract/scope" => "cartesian_basis_and_optional_hamiltonian_bundle",
            "manifest/basis/kind" => String(representation.metadata.basis_kind),
            "manifest/basis/parent_kind" => String(representation.metadata.parent_kind),
            "manifest/basis/final_dimension" => representation.metadata.final_dimension,
            "manifest/basis/working_box_present" => representation.metadata.working_box !== nothing,
            "manifest/ham/present" => include_ham,
        ),
    )
    return meta_values
end

"""
    cartesian_basis_bundle_payload(object; include_ham = true, meta = nothing)

Build one public first-pass Cartesian basis/Hamiltonian bundle payload in
memory without writing a JLD2 file.

The returned named tuple contains:

- `basis`
- `ham`
- `meta`

The writer [`write_cartesian_basis_bundle_jld2`](@ref) stores these under the
grouped JLD2 contract:

- `basis/...`
- `ham/...` when present
- `meta/...`

This first pass supports the pure Cartesian representation families already
covered by the public exact-overlap layer. Hybrid/QW residual final bases are
rejected explicitly for now.
"""
function cartesian_basis_bundle_payload(
    object;
    include_ham::Bool = true,
    meta = nothing,
)
    representation = _cartesian_bundle_representation(object)
    basis_values = _cartesian_basis_values(
        representation;
        final_integral_weights = _cartesian_bundle_integral_weights(object, representation),
    )
    ham_values = include_ham ? _cartesian_ham_values(object, representation) : nothing
    meta_values = _cartesian_bundle_meta_values(
        object,
        representation;
        include_ham = ham_values !== nothing,
        producer_entrypoint = "GaussletBases.write_cartesian_basis_bundle_jld2",
        meta = meta,
    )
    return (
        basis = basis_values,
        ham = ham_values,
        meta = meta_values,
    )
end

function _cartesian_bundle_integral_weights(
    object,
    representation::CartesianBasisRepresentation3D,
)
    object isa CartesianBasisRepresentation3D && return _cartesian_bundle_integral_weights(representation)
    return _cartesian_bundle_integral_weights(object)
end

"""
    write_cartesian_basis_bundle_jld2(path, object; include_ham = true, meta = nothing)

Write one public first-pass Cartesian basis/Hamiltonian bundle to JLD2.

The written file uses a grouped contract:

- `basis/...`
- `ham/...` when a native Hamiltonian payload is available and requested
- `meta/...`
"""
function write_cartesian_basis_bundle_jld2(
    path::AbstractString,
    object;
    include_ham::Bool = true,
    meta = nothing,
)
    data = cartesian_basis_bundle_payload(object; include_ham = include_ham, meta = meta)
    jldopen(path, "w") do file
        _write_prefixed_values!(file, "basis", data.basis)
        data.ham !== nothing && _write_prefixed_values!(file, "ham", data.ham)
        _write_prefixed_values!(file, "meta", data.meta)
    end
    return path
end
