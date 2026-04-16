function _cartesian_bundle_supported_basis(
    representation::CartesianBasisRepresentation3D,
)
    if representation.metadata.parent_kind == :cartesian_product_basis
        return representation
    elseif _cartesian_supports_exact_hybrid_overlap(representation)
        return representation
    elseif representation.metadata.parent_kind == :cartesian_plus_supplement_raw
        throw(
            ArgumentError(
                "cartesian basis bundle export requires explicit mixed raw-space identity for hybrid residual-Gaussian representations; this representation does not carry the exact Cartesian parent and supplement representation data needed for export",
            ),
        )
    end
    throw(
        ArgumentError(
            "cartesian basis bundle export does not yet support parent kind :$(representation.metadata.parent_kind)",
        ),
    )
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
    operators::OrdinaryCartesianOperators3D,
)
    return operators.gaussian_data === nothing && operators.residual_count == 0
end

function _cartesian_bundle_representation(
    operators::OrdinaryCartesianOperators3D,
)
    if _cartesian_qw_bundle_is_pure(operators)
        return _cartesian_bundle_representation(operators.basis)
    end
    return _cartesian_bundle_supported_basis(basis_representation(operators))
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

function _cartesian_write_value!(
    file,
    prefix::AbstractString,
    value,
)
    if value === nothing
        file[string(prefix, "/is_nothing")] = true
    elseif value isa NamedTuple
        for (key, child) in pairs(value)
            _cartesian_write_value!(file, string(prefix, "/", key), child)
        end
    elseif value isa AbstractDict
        for (key, child) in pairs(value)
            _cartesian_write_value!(file, string(prefix, "/", key), child)
        end
    elseif value isa Symbol
        file[prefix] = String(value)
    elseif value isa UnitRange{<:Integer}
        file[prefix] = Int[first(value), last(value)]
    elseif value isa Tuple && all(item -> item isa Symbol, value)
        file[prefix] = String[String(item) for item in value]
    elseif value isa Tuple && all(item -> item isa Integer, value)
        file[prefix] = Int[Int(item) for item in value]
    elseif value isa AbstractArray{<:Symbol}
        file[prefix] = String[String(item) for item in value]
    else
        file[prefix] = value
    end
    return nothing
end

function _cartesian_write_sparse_safe_matrix!(
    file,
    prefix::AbstractString,
    value::AbstractMatrix{<:Real},
)
    if value isa Matrix{Float64} || value isa SparseArrays.SparseMatrixCSC{Float64,Int}
        file[prefix] = value
    else
        file[prefix] = _cartesian_coefficient_map_storage(value)
    end
    return nothing
end

function _cartesian_write_dense_matrix!(
    file,
    prefix::AbstractString,
    value::Matrix{Float64},
)
    file[prefix] = value
    return nothing
end

function _cartesian_write_dense_matrix!(
    file,
    prefix::AbstractString,
    value::AbstractMatrix{<:Real},
)
    file[prefix] = Matrix{Float64}(value)
    return nothing
end

function _cartesian_write_float_vector!(
    file,
    prefix::AbstractString,
    value::Vector{Float64},
)
    file[prefix] = value
    return nothing
end

function _cartesian_write_float_vector!(
    file,
    prefix::AbstractString,
    value::AbstractVector{<:Real},
)
    file[prefix] = Float64[Float64(entry) for entry in value]
    return nothing
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
    haskey(representation.primitive_matrices, :overlap) &&
        (dest[string(prefix, "/primitive_matrices/overlap")] =
            Matrix{Float64}(representation.primitive_matrices.overlap))
    haskey(representation.basis_matrices, :overlap) &&
        (dest[string(prefix, "/basis_matrices/overlap")] =
            Matrix{Float64}(representation.basis_matrices.overlap))
    return dest
end

function _write_cartesian_axis_representation!(
    file,
    prefix::AbstractString,
    representation::BasisRepresentation1D,
)
    metadata = representation.metadata
    primitive_layer = representation.primitive_set

    file[string(prefix, "/format")] = "basis_representation_1d_v1"
    file[string(prefix, "/version")] = 1
    file[string(prefix, "/metadata/basis_kind")] = String(metadata.basis_kind)
    file[string(prefix, "/metadata/has_family_name")] = metadata.family_name !== nothing
    metadata.family_name !== nothing &&
        (file[string(prefix, "/metadata/family_name")] = String(metadata.family_name))
    file[string(prefix, "/metadata/mapping_type")] = string(nameof(typeof(metadata.mapping_value)))
    file[string(prefix, "/metadata/mapping_object")] = metadata.mapping_value
    file[string(prefix, "/metadata/centers")] = metadata.center_data
    file[string(prefix, "/metadata/reference_centers")] = metadata.reference_center_data
    file[string(prefix, "/metadata/integral_weights")] = metadata.integral_weight_data
    file[string(prefix, "/metadata/basis_labels")] = metadata.basis_labels
    file[string(prefix, "/primitive_set/has_name")] = primitive_layer.name_value !== nothing
    primitive_layer.name_value !== nothing &&
        (file[string(prefix, "/primitive_set/name")] = String(primitive_layer.name_value))
    file[string(prefix, "/primitive_set/labels")] = primitive_layer.label_data
    file[string(prefix, "/primitive_set/primitives")] =
        AbstractPrimitiveFunction1D[primitive for primitive in primitives(primitive_layer)]
    file[string(prefix, "/coefficient_matrix")] = Matrix{Float64}(representation.coefficient_matrix)
    haskey(representation.primitive_matrices, :overlap) &&
        (file[string(prefix, "/primitive_matrices/overlap")] =
            Matrix{Float64}(representation.primitive_matrices.overlap))
    haskey(representation.basis_matrices, :overlap) &&
        (file[string(prefix, "/basis_matrices/overlap")] =
            Matrix{Float64}(representation.basis_matrices.overlap))
    return nothing
end

function _cartesian_parent_integral_weights(
    representation::CartesianBasisRepresentation3D,
)
    weights_x = representation.axis_representations.x.metadata.integral_weight_data
    weights_y = representation.axis_representations.y.metadata.integral_weight_data
    weights_z = representation.axis_representations.z.metadata.integral_weight_data
    return _mapped_cartesian_weights(weights_x, weights_y, weights_z)
end

function _cartesian_supplement_axis_integral(
    orbital::CartesianGaussianShellOrbitalRepresentation3D,
    axis::Symbol,
)
    axis_index = axis == :x ? 1 : axis == :y ? 2 : 3
    power = orbital.angular_powers[axis_index]
    value = 0.0
    for (exponent, coefficient) in zip(orbital.exponents, orbital.coefficients)
        prefactor = _qwrg_atomic_shell_prefactor(Float64(exponent), power)
        value +=
            Float64(coefficient) *
            prefactor *
            _qwrg_shifted_gaussian_moment(Float64(exponent), power)
    end
    return value
end

function _cartesian_supplement_orbital_integral(
    orbital::CartesianGaussianShellOrbitalRepresentation3D,
)
    return _cartesian_supplement_axis_integral(orbital, :x) *
           _cartesian_supplement_axis_integral(orbital, :y) *
           _cartesian_supplement_axis_integral(orbital, :z)
end

function _cartesian_supplement_integral_weights(
    supplement::CartesianGaussianShellSupplementRepresentation3D,
)
    return Float64[
        _cartesian_supplement_orbital_integral(orbital) for orbital in supplement.orbitals
    ]
end

function _cartesian_raw_integral_weights(
    representation::CartesianBasisRepresentation3D,
)
    if representation.metadata.parent_kind == :cartesian_product_basis
        return _cartesian_parent_integral_weights(representation)
    elseif representation.metadata.parent_kind == :cartesian_plus_supplement_raw
        _cartesian_bundle_supported_basis(representation)
        cartesian_weights = _cartesian_representation_integral_weights(
            representation.parent_data.cartesian_parent_representation,
        )
        supplement_weights = _cartesian_supplement_integral_weights(
            representation.parent_data.supplement_representation,
        )
        return vcat(cartesian_weights, supplement_weights)
    end
    throw(
        ArgumentError(
            "cartesian basis bundle export does not yet support parent kind :$(representation.metadata.parent_kind)",
        ),
    )
end

function _cartesian_representation_integral_weights(
    representation::CartesianBasisRepresentation3D,
)
    _cartesian_bundle_supported_basis(representation)
    parent_weights = _cartesian_raw_integral_weights(representation)

    if representation.contraction_kind == :identity
        return Vector{Float64}(parent_weights)
    elseif representation.contraction_kind == :dense
        contraction =
            representation.support_indices === nothing ?
            representation.coefficient_matrix :
            representation.coefficient_matrix[representation.support_indices, :]
        support_weights =
            representation.support_indices === nothing ?
            parent_weights :
            parent_weights[representation.support_indices]
        return Vector{Float64}(transpose(contraction) * support_weights)
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
    operators::OrdinaryCartesianOperators3D,
)
    if _cartesian_qw_bundle_is_pure(operators)
        return _cartesian_bundle_integral_weights(operators.basis)
    end
    return _cartesian_representation_integral_weights(basis_representation(operators))
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
        (basis_values["coefficient_matrix"] =
            _cartesian_coefficient_map_storage(representation.coefficient_matrix))

    _cartesian_store_value!(basis_values, "metadata/route", representation.metadata.route_metadata)
    _cartesian_store_axis_representation!(basis_values, "axes/x", representation.axis_representations.x)
    _cartesian_store_axis_representation!(basis_values, "axes/y", representation.axis_representations.y)
    _cartesian_store_axis_representation!(basis_values, "axes/z", representation.axis_representations.z)

    if representation.metadata.parent_kind == :cartesian_plus_supplement_raw
        _cartesian_bundle_supported_basis(representation)
        basis_values["parent/format"] = "cartesian_plus_supplement_raw_v1"
        basis_values["parent/version"] = 1
        hasproperty(representation.parent_data, :hybrid_overlap_kind) &&
            (basis_values["parent/hybrid_overlap_kind"] =
                String(representation.parent_data.hybrid_overlap_kind))
        parent_cartesian_values = _cartesian_basis_values(
            representation.parent_data.cartesian_parent_representation;
            final_integral_weights = _cartesian_representation_integral_weights(
                representation.parent_data.cartesian_parent_representation,
            ),
        )
        for (key, value) in pairs(parent_cartesian_values)
            basis_values[string("parent/cartesian/", key)] = value
        end
        _cartesian_store_supplement_representation!(
            basis_values,
            "parent/supplement",
            representation.parent_data.supplement_representation,
        )
        hasproperty(representation.parent_data, :cartesian_supplement_axis_tables) &&
            _cartesian_store_value!(
                basis_values,
                "parent/cartesian_supplement_axis_tables",
                representation.parent_data.cartesian_supplement_axis_tables,
            )
        hasproperty(representation.parent_data, :exact_cartesian_supplement_overlap) &&
            (basis_values["parent/exact_cartesian_supplement_overlap"] =
                Matrix{Float64}(representation.parent_data.exact_cartesian_supplement_overlap))
        hasproperty(representation.parent_data, :exact_supplement_overlap) &&
            (basis_values["parent/exact_supplement_overlap"] =
                Matrix{Float64}(representation.parent_data.exact_supplement_overlap))
        hasproperty(representation.parent_data, :cartesian_supplement_overlap) &&
            (basis_values["parent/cartesian_supplement_overlap"] =
                Matrix{Float64}(representation.parent_data.cartesian_supplement_overlap))
        hasproperty(representation.parent_data, :supplement_overlap) &&
            (basis_values["parent/supplement_overlap"] =
                Matrix{Float64}(representation.parent_data.supplement_overlap))
    end

    return basis_values
end

function _write_cartesian_basis_group!(
    file,
    prefix::AbstractString,
    representation::CartesianBasisRepresentation3D;
    final_integral_weights::AbstractVector{<:Real},
)
    file[string(prefix, "/format")] = "cartesian_basis_bundle_v1"
    file[string(prefix, "/version")] = 1
    file[string(prefix, "/basis_kind")] = String(representation.metadata.basis_kind)
    file[string(prefix, "/parent_kind")] = String(representation.metadata.parent_kind)
    file[string(prefix, "/axis_sharing")] = String(representation.metadata.axis_sharing)
    file[string(prefix, "/parent_axis_counts")] = Int[representation.metadata.parent_axis_counts...]
    file[string(prefix, "/parent_dimension")] = representation.metadata.parent_dimension
    file[string(prefix, "/final_dimension")] = representation.metadata.final_dimension
    file[string(prefix, "/working_box_present")] = representation.metadata.working_box !== nothing
    file[string(prefix, "/working_box_bounds")] =
        _cartesian_working_box_matrix(representation.metadata.working_box)
    file[string(prefix, "/contraction_kind")] = String(representation.contraction_kind)
    file[string(prefix, "/basis_labels")] = representation.metadata.basis_labels
    file[string(prefix, "/basis_centers")] = representation.metadata.basis_centers
    file[string(prefix, "/final_integral_weights")] =
        Float64[Float64(value) for value in final_integral_weights]
    file[string(prefix, "/parent_labels")] = representation.parent_labels
    file[string(prefix, "/parent_centers")] = representation.parent_centers
    file[string(prefix, "/support_indices_present")] = representation.support_indices !== nothing
    file[string(prefix, "/support_indices")] =
        representation.support_indices === nothing ? zeros(Int, 0) : Vector{Int}(representation.support_indices)
    file[string(prefix, "/support_states_present")] = representation.support_states !== nothing
    file[string(prefix, "/support_states")] = _cartesian_support_state_matrix(representation.support_states)

    if representation.coefficient_matrix !== nothing
        _cartesian_write_sparse_safe_matrix!(
            file,
            string(prefix, "/coefficient_matrix"),
            representation.coefficient_matrix,
        )
    end

    _cartesian_write_value!(file, string(prefix, "/metadata/route"), representation.metadata.route_metadata)
    _write_cartesian_axis_representation!(file, string(prefix, "/axes/x"), representation.axis_representations.x)
    _write_cartesian_axis_representation!(file, string(prefix, "/axes/y"), representation.axis_representations.y)
    _write_cartesian_axis_representation!(file, string(prefix, "/axes/z"), representation.axis_representations.z)

    if representation.metadata.parent_kind == :cartesian_plus_supplement_raw
        _cartesian_bundle_supported_basis(representation)
        file[string(prefix, "/parent/format")] = "cartesian_plus_supplement_raw_v1"
        file[string(prefix, "/parent/version")] = 1
        hasproperty(representation.parent_data, :hybrid_overlap_kind) &&
            (file[string(prefix, "/parent/hybrid_overlap_kind")] =
                String(representation.parent_data.hybrid_overlap_kind))
        _write_cartesian_basis_group!(
            file,
            string(prefix, "/parent/cartesian"),
            representation.parent_data.cartesian_parent_representation;
            final_integral_weights = _cartesian_representation_integral_weights(
                representation.parent_data.cartesian_parent_representation,
            ),
        )
        _write_cartesian_supplement_representation!(
            file,
            string(prefix, "/parent/supplement"),
            representation.parent_data.supplement_representation,
        )
        hasproperty(representation.parent_data, :cartesian_supplement_axis_tables) &&
            _cartesian_write_value!(
                file,
                string(prefix, "/parent/cartesian_supplement_axis_tables"),
                representation.parent_data.cartesian_supplement_axis_tables,
            )
        hasproperty(representation.parent_data, :exact_cartesian_supplement_overlap) &&
            (file[string(prefix, "/parent/exact_cartesian_supplement_overlap")] =
                representation.parent_data.exact_cartesian_supplement_overlap)
        hasproperty(representation.parent_data, :exact_supplement_overlap) &&
            (file[string(prefix, "/parent/exact_supplement_overlap")] =
                representation.parent_data.exact_supplement_overlap)
        hasproperty(representation.parent_data, :cartesian_supplement_overlap) &&
            (file[string(prefix, "/parent/cartesian_supplement_overlap")] =
                representation.parent_data.cartesian_supplement_overlap)
        hasproperty(representation.parent_data, :supplement_overlap) &&
            (file[string(prefix, "/parent/supplement_overlap")] =
                representation.parent_data.supplement_overlap)
    end

    return nothing
end

function _cartesian_store_supplement_representation!(
    dest::Dict{String,Any},
    prefix::AbstractString,
    supplement::CartesianGaussianShellSupplementRepresentation3D,
)
    dest[string(prefix, "/format")] = "cartesian_gaussian_shell_supplement_v1"
    dest[string(prefix, "/version")] = 1
    dest[string(prefix, "/supplement_kind")] = String(supplement.supplement_kind)
    _cartesian_store_value!(dest, string(prefix, "/metadata"), supplement.metadata)
    dest[string(prefix, "/orbital_count")] = length(supplement.orbitals)
    for (index, orbital) in pairs(supplement.orbitals)
        orbital_prefix = string(prefix, "/orbitals/", index)
        dest[string(orbital_prefix, "/label")] = orbital.label
        dest[string(orbital_prefix, "/angular_powers")] = Int[orbital.angular_powers...]
        dest[string(orbital_prefix, "/center")] = Float64[orbital.center...]
        dest[string(orbital_prefix, "/exponents")] = Float64[orbital.exponents...]
        dest[string(orbital_prefix, "/coefficients")] = Float64[orbital.coefficients...]
        dest[string(orbital_prefix, "/primitive_normalization")] =
            String(orbital.primitive_normalization)
    end
    return dest
end

function _write_cartesian_supplement_representation!(
    file,
    prefix::AbstractString,
    supplement::CartesianGaussianShellSupplementRepresentation3D,
)
    file[string(prefix, "/format")] = "cartesian_gaussian_shell_supplement_v1"
    file[string(prefix, "/version")] = 1
    file[string(prefix, "/supplement_kind")] = String(supplement.supplement_kind)
    _cartesian_write_value!(file, string(prefix, "/metadata"), supplement.metadata)
    file[string(prefix, "/orbital_count")] = length(supplement.orbitals)
    for (index, orbital) in pairs(supplement.orbitals)
        orbital_prefix = string(prefix, "/orbitals/", index)
        file[string(orbital_prefix, "/label")] = orbital.label
        file[string(orbital_prefix, "/angular_powers")] = Int[orbital.angular_powers...]
        file[string(orbital_prefix, "/center")] = Float64[orbital.center...]
        file[string(orbital_prefix, "/exponents")] = Float64[orbital.exponents...]
        file[string(orbital_prefix, "/coefficients")] = Float64[orbital.coefficients...]
        file[string(orbital_prefix, "/primitive_normalization")] =
            String(orbital.primitive_normalization)
    end
    return nothing
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
    operators::OrdinaryCartesianOperators3D,
    representation::CartesianBasisRepresentation3D,
)
    values = Dict{String,Any}(
        "format" => "cartesian_hamiltonian_bundle_v1",
        "version" => 1,
        "object_type" => string(nameof(typeof(operators))),
        "model_kind" => "ordinary_cartesian_operators",
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
        "basis_integral_weights" => _cartesian_bundle_integral_weights(operators, representation),
        "expansion/exponents" => copy(operators.expansion.exponents),
        "expansion/coefficients" => copy(operators.expansion.coefficients),
        "residual_count" => operators.residual_count,
    )
    if !isnothing(operators.nuclear_charges)
        values["default_nuclear_charges"] = copy(operators.nuclear_charges)
        values["nuclear_term_storage"] = String(operators.nuclear_term_storage)
    end
    if !isnothing(operators.kinetic_one_body)
        values["kinetic_one_body"] = Matrix{Float64}(operators.kinetic_one_body)
    end
    if !isnothing(operators.nuclear_one_body_by_center)
        values["nuclear_one_body_by_center/count"] = length(operators.nuclear_one_body_by_center)
        for (index, matrix) in pairs(operators.nuclear_one_body_by_center)
            values["nuclear_one_body_by_center/$(index)"] = Matrix{Float64}(matrix)
        end
    end
    return values
end

function _write_cartesian_ham_group!(
    file,
    ::CartesianBasisRepresentation3D,
    representation::CartesianBasisRepresentation3D,
)
    _cartesian_bundle_supported_basis(representation)
    return false
end

function _write_cartesian_ham_group!(
    file,
    basis::Union{
        BondAlignedDiatomicQWBasis3D,
        BondAlignedHomonuclearChainQWBasis3D,
        AxisAlignedHomonuclearSquareLatticeQWBasis3D,
        _NestedFixedBlock3D,
    },
    representation::CartesianBasisRepresentation3D,
)
    _cartesian_bundle_supported_basis(representation)
    return false
end

function _write_cartesian_ham_group!(
    file,
    operators::OrdinaryCartesianIDAOperators,
    representation::CartesianBasisRepresentation3D,
)
    file["ham/format"] = "cartesian_hamiltonian_bundle_v1"
    file["ham/version"] = 1
    file["ham/object_type"] = string(nameof(typeof(operators)))
    file["ham/model_kind"] = "ordinary_cartesian_ida"
    file["ham/interaction_model"] = "density_density"
    file["ham/interaction_treatment"] = String(operators.interaction_treatment)
    file["ham/backend"] = String(operators.backend)
    file["ham/overlap_key"] = "overlap"
    file["ham/onebody_key"] = "one_body_hamiltonian"
    file["ham/interaction_key"] = "interaction_matrix"
    file["ham/basis_integral_weights_key"] = "basis/final_integral_weights"
    _cartesian_write_dense_matrix!(file, "ham/overlap", operators.overlap_3d)
    _cartesian_write_dense_matrix!(file, "ham/one_body_hamiltonian", operators.one_body_hamiltonian)
    _cartesian_write_dense_matrix!(file, "ham/interaction_matrix", operators.interaction_matrix)
    file["ham/orbital_labels"] = String[String(orbital.label) for orbital in orbitals(operators)]
    _cartesian_write_dense_matrix!(file, "ham/basis_centers", representation.metadata.basis_centers)
    _cartesian_write_float_vector!(file, "ham/basis_integral_weights", operators.weight_3d)
    _cartesian_write_float_vector!(file, "ham/expansion/exponents", operators.expansion.exponents)
    _cartesian_write_float_vector!(file, "ham/expansion/coefficients", operators.expansion.coefficients)
    return true
end

function _write_cartesian_ham_group!(
    file,
    operators::OrdinaryCartesianOperators3D,
    representation::CartesianBasisRepresentation3D,
)
    file["ham/format"] = "cartesian_hamiltonian_bundle_v1"
    file["ham/version"] = 1
    file["ham/object_type"] = string(nameof(typeof(operators)))
    file["ham/model_kind"] = "ordinary_cartesian_operators"
    file["ham/interaction_model"] = "density_density"
    file["ham/interaction_treatment"] = String(operators.interaction_treatment)
    file["ham/gausslet_backend"] = String(operators.gausslet_backend)
    file["ham/overlap_key"] = "overlap"
    file["ham/onebody_key"] = "one_body_hamiltonian"
    file["ham/interaction_key"] = "interaction_matrix"
    file["ham/basis_integral_weights_key"] = "basis/final_integral_weights"
    _cartesian_write_dense_matrix!(file, "ham/overlap", operators.overlap)
    _cartesian_write_dense_matrix!(file, "ham/one_body_hamiltonian", operators.one_body_hamiltonian)
    _cartesian_write_dense_matrix!(file, "ham/interaction_matrix", operators.interaction_matrix)
    file["ham/orbital_labels"] = String[String(orbital.label) for orbital in orbitals(operators)]
    _cartesian_write_dense_matrix!(file, "ham/basis_centers", representation.metadata.basis_centers)
    _cartesian_write_float_vector!(
        file,
        "ham/basis_integral_weights",
        _cartesian_bundle_integral_weights(operators, representation),
    )
    _cartesian_write_float_vector!(file, "ham/expansion/exponents", operators.expansion.exponents)
    _cartesian_write_float_vector!(file, "ham/expansion/coefficients", operators.expansion.coefficients)
    file["ham/residual_count"] = operators.residual_count
    if !isnothing(operators.nuclear_charges)
        _cartesian_write_float_vector!(file, "ham/default_nuclear_charges", operators.nuclear_charges)
        file["ham/nuclear_term_storage"] = String(operators.nuclear_term_storage)
    end
    if !isnothing(operators.kinetic_one_body)
        _cartesian_write_dense_matrix!(file, "ham/kinetic_one_body", operators.kinetic_one_body)
    end
    if !isnothing(operators.nuclear_one_body_by_center)
        file["ham/nuclear_one_body_by_center/count"] = length(operators.nuclear_one_body_by_center)
        for (index, matrix) in pairs(operators.nuclear_one_body_by_center)
            _cartesian_write_dense_matrix!(
                file,
                "ham/nuclear_one_body_by_center/$(index)",
                matrix,
            )
        end
    end
    return true
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

function _write_cartesian_meta_group!(
    file,
    object,
    representation::CartesianBasisRepresentation3D;
    include_ham::Bool,
    producer_entrypoint::AbstractString,
    meta = nothing,
)
    for (key, value) in pairs(_normalize_meta_dict(meta))
        _cartesian_write_value!(file, string("meta/", key), value)
    end
    _cartesian_write_value!(file, "meta/producer", producer_entrypoint)
    _cartesian_write_value!(file, "meta/producer_type", string(nameof(typeof(object))))
    _cartesian_write_value!(file, "meta/basis_kind", String(representation.metadata.basis_kind))
    _cartesian_write_value!(file, "meta/parent_kind", String(representation.metadata.parent_kind))
    _cartesian_write_value!(file, "meta/has_ham", include_ham)
    _cartesian_write_value!(file, "meta/final_dimension", representation.metadata.final_dimension)
    _cartesian_write_value!(file, "meta/parent_dimension", representation.metadata.parent_dimension)
    _cartesian_write_value!(file, "meta/manifest/producer/package", "GaussletBases")
    _cartesian_write_value!(file, "meta/manifest/producer/version", string(Base.pkgversion(@__MODULE__)))
    _cartesian_write_value!(file, "meta/manifest/producer/entrypoint", producer_entrypoint)
    _cartesian_write_value!(file, "meta/manifest/contract/format", "cartesian_basis_bundle_v1")
    _cartesian_write_value!(file, "meta/manifest/contract/version", 1)
    _cartesian_write_value!(file, "meta/manifest/contract/status", "public_first_pass")
    _cartesian_write_value!(
        file,
        "meta/manifest/contract/scope",
        "cartesian_basis_and_optional_hamiltonian_bundle",
    )
    _cartesian_write_value!(file, "meta/manifest/basis/kind", String(representation.metadata.basis_kind))
    _cartesian_write_value!(file, "meta/manifest/basis/parent_kind", String(representation.metadata.parent_kind))
    _cartesian_write_value!(file, "meta/manifest/basis/final_dimension", representation.metadata.final_dimension)
    _cartesian_write_value!(
        file,
        "meta/manifest/basis/working_box_present",
        representation.metadata.working_box !== nothing,
    )
    _cartesian_write_value!(file, "meta/manifest/ham/present", include_ham)
    return nothing
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

This first pass supports the Cartesian representation families already covered
by the public exact-overlap layer, including hybrid residual-Gaussian final
bases whose public representation carries explicit mixed raw-space identity.
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
    representation = _cartesian_bundle_representation(object)
    final_integral_weights = _cartesian_bundle_integral_weights(object, representation)
    jldopen(path, "w") do file
        _write_cartesian_basis_group!(
            file,
            "basis",
            representation;
            final_integral_weights = final_integral_weights,
        )
        ham_written = include_ham ? _write_cartesian_ham_group!(file, object, representation) : false
        _write_cartesian_meta_group!(
            file,
            object,
            representation;
            include_ham = ham_written,
            producer_entrypoint = "GaussletBases.write_cartesian_basis_bundle_jld2",
            meta = meta,
        )
    end
    return path
end
