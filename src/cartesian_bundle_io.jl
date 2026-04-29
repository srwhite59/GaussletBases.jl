"""
    CartesianBasisBundle3D

Read-only in-memory view of one first-pass Cartesian JLD2 basis bundle.

The object stores:

- the source bundle path
- the reconstructed public `CartesianBasisRepresentation3D`
- optional `ham/...` payload values when present
- `meta/...` values
- light bundle diagnostics useful for disk-level overlap/projector debugging

This read-side surface intentionally stays thin. Exact overlap and projector
logic continue to live on the public in-memory representation API.
"""
struct CartesianBasisBundle3D
    path::String
    basis::CartesianBasisRepresentation3D
    ham::Union{Nothing,Dict{String,Any}}
    meta::Dict{String,Any}
    diagnostics::NamedTuple
end

function Base.show(io::IO, bundle::CartesianBasisBundle3D)
    print(
        io,
        "CartesianBasisBundle3D(path=",
        repr(bundle.path),
        ", kind=:",
        bundle.basis.metadata.basis_kind,
        ", nfinal=",
        bundle.basis.metadata.final_dimension,
        ", has_ham=",
        bundle.ham !== nothing,
        ")",
    )
end

basis_representation(bundle::CartesianBasisBundle3D) = bundle.basis
basis_metadata(bundle::CartesianBasisBundle3D) = bundle.basis.metadata

function _cartesian_bundle_default_nuclear_charges(
    ham::Union{Nothing,AbstractDict{String,Any}},
)
    isnothing(ham) && return nothing
    haskey(ham, "default_nuclear_charges") || return nothing
    return Float64[Float64(value) for value in ham["default_nuclear_charges"]]
end

function _cartesian_bundle_nuclear_one_body_by_center(
    ham::AbstractDict{String,Any},
)
    haskey(ham, "nuclear_one_body_by_center/count") || return nothing
    count = Int(ham["nuclear_one_body_by_center/count"])
    return [Matrix{Float64}(ham["nuclear_one_body_by_center/$(index)"]) for index in 1:count]
end

"""
    assembled_one_body_hamiltonian(bundle::CartesianBasisBundle3D;
                                   nuclear_charges = default bundle charges)

Reassemble a stored one-body Hamiltonian from the bundle's final-basis kinetic
and per-center nuclear terms when those terms were retained.

If the bundle does not carry per-center nuclear terms, the stored total
`ham["one_body_hamiltonian"]` can only be returned for the original default
nuclear charges.
"""
function assembled_one_body_hamiltonian(
    bundle::CartesianBasisBundle3D;
    nuclear_charges = _cartesian_bundle_default_nuclear_charges(bundle.ham),
)
    bundle.ham === nothing && throw(
        ArgumentError("this Cartesian basis bundle does not carry a Hamiltonian payload"),
    )
    default_charges = _cartesian_bundle_default_nuclear_charges(bundle.ham)
    kinetic_one_body =
        haskey(bundle.ham, "kinetic_one_body") ?
        Matrix{Float64}(bundle.ham["kinetic_one_body"]) :
        nothing
    nuclear_one_body_by_center = _cartesian_bundle_nuclear_one_body_by_center(bundle.ham)
    if isnothing(kinetic_one_body) || isnothing(nuclear_one_body_by_center)
        if nuclear_charges === default_charges ||
           _same_nuclear_charge_configuration(nuclear_charges, default_charges)
            return Matrix{Float64}(bundle.ham["one_body_hamiltonian"])
        end
        throw(
            ArgumentError(
                "this bundle does not retain per-center nuclear one-body terms; rebuild or export with nuclear_term_storage = :by_center",
            ),
        )
    end
    charges = isnothing(nuclear_charges) ? default_charges : nuclear_charges
    isnothing(charges) && return Matrix{Float64}(bundle.ham["one_body_hamiltonian"])
    return _assemble_one_body_hamiltonian(
        kinetic_one_body,
        nuclear_one_body_by_center,
        charges,
    )
end

function _cartesian_jld_key(raw_key)
    return raw_key isa AbstractVector ? join(string.(raw_key), "/") : string(raw_key)
end

function _cartesian_jld_has_group(container, name::AbstractString)
    return any(_cartesian_jld_key(key) == name for key in keys(container))
end

function _cartesian_jld_group_values(group::JLD2.Group)
    values = Dict{String,Any}()
    for raw_key in keys(group)
        key = _cartesian_jld_key(raw_key)
        value = group[key]
        if value isa JLD2.Group
            for (child_key, child_value) in pairs(_cartesian_jld_group_values(value))
                values[string(key, "/", child_key)] = child_value
            end
        else
            values[key] = value
        end
    end
    return values
end

function _cartesian_bundle_required_value(
    values::AbstractDict{String},
    key::AbstractString,
)
    haskey(values, key) || throw(
        ArgumentError("Cartesian basis bundle is missing required key $(repr(key))"),
    )
    return values[key]
end

function _cartesian_bundle_prefix_values(
    values::AbstractDict{String},
    prefix::AbstractString,
)
    prefix_with_sep = string(prefix, "/")
    out = Dict{String,Any}()
    for (key, value) in pairs(values)
        startswith(key, prefix_with_sep) || continue
        out[key[length(prefix_with_sep)+1:end]] = value
    end
    return out
end

function _cartesian_bundle_restore_value(value; strings_to_symbols::Bool = false)
    if strings_to_symbols && value isa AbstractString
        return Symbol(String(value))
    elseif strings_to_symbols && value isa AbstractVector{<:AbstractString}
        return Symbol[Symbol(String(item)) for item in value]
    elseif value isa AbstractVector
        return copy(value)
    elseif value isa AbstractMatrix
        return copy(value)
    end
    return value
end

function _cartesian_bundle_tree_to_namedtuple(tree::Dict{Symbol,Any})
    if length(tree) == 1 && get(tree, :is_nothing, false) === true
        return nothing
    end
    keys_sorted = sort!(collect(keys(tree)); by = string)
    values = map(keys_sorted) do key
        value = tree[key]
        value isa Dict{Symbol,Any} ? _cartesian_bundle_tree_to_namedtuple(value) : value
    end
    return NamedTuple{Tuple(keys_sorted)}(Tuple(values))
end

function _cartesian_bundle_namedtuple(
    values::AbstractDict{String};
    strings_to_symbols::Bool = false,
)
    isempty(values) && return NamedTuple()
    tree = Dict{Symbol,Any}()
    for key in sort!(collect(keys(values)))
        parts = Symbol.(split(key, '/'))
        node = tree
        for part in parts[1:end-1]
            child = get!(node, part) do
                Dict{Symbol,Any}()
            end
            child isa Dict{Symbol,Any} || throw(
                ArgumentError(
                    "Cartesian basis bundle key collision while reconstructing $(repr(key))",
                ),
            )
            node = child
        end
        node[parts[end]] = _cartesian_bundle_restore_value(
            values[key];
            strings_to_symbols = strings_to_symbols,
        )
    end
    return _cartesian_bundle_tree_to_namedtuple(tree)
end

function _cartesian_bundle_axis_representation(
    values::AbstractDict{String},
    axis_label::AbstractString,
)
    axis_values = _cartesian_bundle_prefix_values(values, string("axes/", axis_label))
    String(_cartesian_bundle_required_value(axis_values, "format")) == "basis_representation_1d_v1" || throw(
        ArgumentError(
            "Cartesian basis bundle axis $(axis_label) format is unsupported; expected basis_representation_1d_v1",
        ),
    )
    Int(_cartesian_bundle_required_value(axis_values, "version")) == 1 || throw(
        ArgumentError(
            "Cartesian basis bundle axis $(axis_label) version is unsupported; expected version 1",
        ),
    )

    basis_kind = Symbol(String(_cartesian_bundle_required_value(axis_values, "metadata/basis_kind")))
    family_name =
        Bool(_cartesian_bundle_required_value(axis_values, "metadata/has_family_name")) ?
        Symbol(String(_cartesian_bundle_required_value(axis_values, "metadata/family_name"))) :
        nothing
    mapping_value = _cartesian_bundle_required_value(axis_values, "metadata/mapping_object")
    centers =
        Float64[Float64(value) for value in _cartesian_bundle_required_value(axis_values, "metadata/centers")]
    reference_centers = Float64[
        Float64(value) for value in _cartesian_bundle_required_value(axis_values, "metadata/reference_centers")
    ]
    integral_weights = Float64[
        Float64(value) for value in _cartesian_bundle_required_value(axis_values, "metadata/integral_weights")
    ]
    basis_labels =
        String[String(label) for label in _cartesian_bundle_required_value(axis_values, "metadata/basis_labels")]
    primitive_name =
        Bool(_cartesian_bundle_required_value(axis_values, "primitive_set/has_name")) ?
        Symbol(String(_cartesian_bundle_required_value(axis_values, "primitive_set/name"))) :
        nothing
    primitive_labels = String[
        String(label) for label in _cartesian_bundle_required_value(axis_values, "primitive_set/labels")
    ]
    primitive_data = AbstractPrimitiveFunction1D[
        primitive for primitive in _cartesian_bundle_required_value(axis_values, "primitive_set/primitives")
    ]
    primitive_layer = PrimitiveSet1D(
        primitive_data;
        name = primitive_name,
        labels = primitive_labels,
    )
    coefficient_matrix =
        Matrix{Float64}(_cartesian_bundle_required_value(axis_values, "coefficient_matrix"))
    metadata = BasisMetadata1D(
        basis_kind,
        family_name,
        mapping_value,
        centers,
        reference_centers,
        integral_weights,
        basis_labels,
        primitive_layer,
        coefficient_matrix,
    )
    primitive_matrices =
        haskey(axis_values, "primitive_matrices/overlap") ?
        (
            overlap = Matrix{Float64}(
                _cartesian_bundle_required_value(axis_values, "primitive_matrices/overlap"),
            ),
        ) :
        NamedTuple()
    basis_matrices =
        haskey(axis_values, "basis_matrices/overlap") ?
        (
            overlap = Matrix{Float64}(
                _cartesian_bundle_required_value(axis_values, "basis_matrices/overlap"),
            ),
        ) :
        NamedTuple()
    return BasisRepresentation1D(
        metadata,
        primitive_layer,
        coefficient_matrix,
        primitive_matrices,
        basis_matrices,
    )
end

function _cartesian_bundle_working_box(values::AbstractDict{String})
    Bool(_cartesian_bundle_required_value(values, "working_box_present")) || return nothing
    bounds = Matrix{Int}(_cartesian_bundle_required_value(values, "working_box_bounds"))
    size(bounds) == (3, 2) || throw(
        ArgumentError("Cartesian basis bundle working_box_bounds must have size (3, 2)"),
    )
    return (
        Int(bounds[1, 1]):Int(bounds[1, 2]),
        Int(bounds[2, 1]):Int(bounds[2, 2]),
        Int(bounds[3, 1]):Int(bounds[3, 2]),
    )
end

function _cartesian_bundle_support_states(values::AbstractDict{String})
    Bool(_cartesian_bundle_required_value(values, "support_states_present")) || return nothing
    matrix = Matrix{Int}(_cartesian_bundle_required_value(values, "support_states"))
    size(matrix, 2) == 3 || throw(
        ArgumentError("Cartesian basis bundle support_states must have three columns"),
    )
    return NTuple{3,Int}[
        (Int(matrix[row, 1]), Int(matrix[row, 2]), Int(matrix[row, 3])) for row in axes(matrix, 1)
    ]
end

function _cartesian_bundle_basis_representation(
    values::AbstractDict{String},
)
    String(_cartesian_bundle_required_value(values, "format")) == "cartesian_basis_bundle_v1" || throw(
        ArgumentError(
            "Cartesian basis bundle format is unsupported; expected cartesian_basis_bundle_v1",
        ),
    )
    Int(_cartesian_bundle_required_value(values, "version")) == 1 || throw(
        ArgumentError("Cartesian basis bundle version is unsupported; expected version 1"),
    )

    axis_representations = (
        x = _cartesian_bundle_axis_representation(values, "x"),
        y = _cartesian_bundle_axis_representation(values, "y"),
        z = _cartesian_bundle_axis_representation(values, "z"),
    )
    axis_metadata = (
        x = axis_representations.x.metadata,
        y = axis_representations.y.metadata,
        z = axis_representations.z.metadata,
    )
    route_metadata = _cartesian_bundle_namedtuple(
        _cartesian_bundle_prefix_values(values, "metadata/route");
        strings_to_symbols = true,
    )
    metadata = CartesianBasisMetadata3D(
        Symbol(String(_cartesian_bundle_required_value(values, "basis_kind"))),
        Symbol(String(_cartesian_bundle_required_value(values, "axis_sharing"))),
        axis_metadata,
        Symbol(String(_cartesian_bundle_required_value(values, "parent_kind"))),
        Tuple(Int.(Vector{Int}(_cartesian_bundle_required_value(values, "parent_axis_counts")))),
        Int(_cartesian_bundle_required_value(values, "parent_dimension")),
        Int(_cartesian_bundle_required_value(values, "final_dimension")),
        _cartesian_bundle_working_box(values),
        String[String(label) for label in _cartesian_bundle_required_value(values, "basis_labels")],
        Matrix{Float64}(_cartesian_bundle_required_value(values, "basis_centers")),
        route_metadata,
    )
    support_indices =
        Bool(_cartesian_bundle_required_value(values, "support_indices_present")) ?
        Vector{Int}(_cartesian_bundle_required_value(values, "support_indices")) :
        nothing
    coefficient_matrix =
        haskey(values, "coefficient_matrix") ?
        _cartesian_coefficient_map_storage(
            _cartesian_bundle_required_value(values, "coefficient_matrix"),
        ) :
        nothing
    parent_data =
        metadata.parent_kind == :cartesian_plus_supplement_raw ?
        let
            cartesian_parent_representation = _cartesian_bundle_basis_representation(
                _cartesian_bundle_prefix_values(values, "parent/cartesian"),
            )
            supplement_representation = _cartesian_bundle_supplement_representation(
                _cartesian_bundle_prefix_values(values, "parent/supplement"),
            )
            exact_cartesian_supplement_overlap =
                haskey(values, "parent/exact_cartesian_supplement_overlap") ?
                Matrix{Float64}(
                    _cartesian_bundle_required_value(values, "parent/exact_cartesian_supplement_overlap"),
                ) :
                haskey(values, "parent/cartesian_supplement_overlap") ?
                Matrix{Float64}(
                    _cartesian_bundle_required_value(values, "parent/cartesian_supplement_overlap"),
                ) :
                nothing
            exact_supplement_overlap =
                haskey(values, "parent/exact_supplement_overlap") ?
                Matrix{Float64}(
                    _cartesian_bundle_required_value(values, "parent/exact_supplement_overlap"),
                ) :
                haskey(values, "parent/supplement_overlap") ?
                Matrix{Float64}(
                    _cartesian_bundle_required_value(values, "parent/supplement_overlap"),
                ) :
                nothing
            cartesian_supplement_axis_tables =
                haskey(values, "parent/cartesian_supplement_axis_tables/x") ?
                (
                    x = Matrix{Float64}(
                        _cartesian_bundle_required_value(
                            values,
                            "parent/cartesian_supplement_axis_tables/x",
                        ),
                    ),
                    y = Matrix{Float64}(
                        _cartesian_bundle_required_value(
                            values,
                            "parent/cartesian_supplement_axis_tables/y",
                        ),
                    ),
                    z = Matrix{Float64}(
                        _cartesian_bundle_required_value(
                            values,
                            "parent/cartesian_supplement_axis_tables/z",
                        ),
                    ),
                ) :
                nothing
            (
                ;
                cartesian_parent_representation = cartesian_parent_representation,
                supplement_representation = supplement_representation,
                hybrid_overlap_kind =
                    haskey(values, "parent/hybrid_overlap_kind") ?
                    Symbol(String(_cartesian_bundle_required_value(values, "parent/hybrid_overlap_kind"))) :
                    :unspecified,
                factorized_cartesian_parent_basis =
                    _cartesian_factorized_parent_basis(cartesian_parent_representation),
                cartesian_supplement_axis_tables = cartesian_supplement_axis_tables,
                exact_cartesian_supplement_overlap = exact_cartesian_supplement_overlap,
                exact_supplement_overlap = exact_supplement_overlap,
                cartesian_supplement_overlap = exact_cartesian_supplement_overlap,
                supplement_overlap = exact_supplement_overlap,
            )
        end :
        (;)
    return CartesianBasisRepresentation3D(
        metadata,
        axis_representations,
        Symbol(String(_cartesian_bundle_required_value(values, "contraction_kind"))),
        coefficient_matrix,
        String[String(label) for label in _cartesian_bundle_required_value(values, "parent_labels")],
        Matrix{Float64}(_cartesian_bundle_required_value(values, "parent_centers")),
        support_indices,
        _cartesian_bundle_support_states(values),
        parent_data,
    )
end

function _cartesian_bundle_supplement_representation(
    values::AbstractDict{String},
)
    String(_cartesian_bundle_required_value(values, "format")) ==
    "cartesian_gaussian_shell_supplement_v1" || throw(
        ArgumentError(
            "Cartesian basis bundle supplement format is unsupported; expected cartesian_gaussian_shell_supplement_v1",
        ),
    )
    Int(_cartesian_bundle_required_value(values, "version")) == 1 || throw(
        ArgumentError(
            "Cartesian basis bundle supplement version is unsupported; expected version 1",
        ),
    )
    metadata = _cartesian_bundle_namedtuple(
        _cartesian_bundle_prefix_values(values, "metadata");
        strings_to_symbols = true,
    )
    orbital_count = Int(_cartesian_bundle_required_value(values, "orbital_count"))
    orbitals = CartesianGaussianShellOrbitalRepresentation3D[]
    sizehint!(orbitals, orbital_count)
    for index in 1:orbital_count
        orbital_values = _cartesian_bundle_prefix_values(values, string("orbitals/", index))
        push!(
            orbitals,
            CartesianGaussianShellOrbitalRepresentation3D(
                String(_cartesian_bundle_required_value(orbital_values, "label")),
                Tuple(
                    Int.(
                        Vector{Int}(
                            _cartesian_bundle_required_value(
                                orbital_values,
                                "angular_powers",
                            ),
                        ),
                    ),
                ),
                Tuple(
                    Float64.(
                        Vector{Float64}(
                            _cartesian_bundle_required_value(orbital_values, "center"),
                        ),
                    ),
                ),
                Float64[
                    Float64(value) for
                    value in _cartesian_bundle_required_value(orbital_values, "exponents")
                ],
                Float64[
                    Float64(value) for
                    value in _cartesian_bundle_required_value(orbital_values, "coefficients")
                ],
                Symbol(
                    String(
                        _cartesian_bundle_required_value(
                            orbital_values,
                            "primitive_normalization",
                        ),
                    ),
                ),
            ),
        )
    end
    return CartesianGaussianShellSupplementRepresentation3D(
        Symbol(String(_cartesian_bundle_required_value(values, "supplement_kind"))),
        orbitals,
        metadata,
    )
end

function _cartesian_bundle_diagnostics(
    path::AbstractString,
    basis::CartesianBasisRepresentation3D,
    ham::Union{Nothing,Dict{String,Any}},
)
    return (
        path = path,
        basis_kind = basis.metadata.basis_kind,
        parent_kind = basis.metadata.parent_kind,
        final_dimension = basis.metadata.final_dimension,
        parent_dimension = basis.metadata.parent_dimension,
        has_ham = ham !== nothing,
        ham_model_kind = ham === nothing ? nothing : get(ham, "model_kind", nothing),
    )
end

"""
    read_cartesian_basis_bundle(path)

Read one first-pass Cartesian JLD2 bundle from `path` and reconstruct its
public `CartesianBasisRepresentation3D`.

The returned [`CartesianBasisBundle3D`](@ref) keeps the loaded basis
representation plus any `ham/...` and `meta/...` payloads present in the file.
Exact overlap/projector work continues to use the public in-memory basis API on
the reconstructed representation.
"""
function read_cartesian_basis_bundle(path::AbstractString)
    bundle_path = abspath(path)
    return jldopen(path, "r") do file
        _cartesian_jld_has_group(file, "basis") || throw(
            ArgumentError("Cartesian basis bundle file $(repr(bundle_path)) is missing the basis/ group"),
        )
        basis_values = _cartesian_jld_group_values(file["basis"])
        ham_values =
            _cartesian_jld_has_group(file, "ham") ?
            _cartesian_jld_group_values(file["ham"]) :
            nothing
        meta_values =
            _cartesian_jld_has_group(file, "meta") ?
            _cartesian_jld_group_values(file["meta"]) :
            Dict{String,Any}()
        basis = _cartesian_bundle_basis_representation(basis_values)
        diagnostics = _cartesian_bundle_diagnostics(bundle_path, basis, ham_values)
        return CartesianBasisBundle3D(bundle_path, basis, ham_values, meta_values, diagnostics)
    end
end

"""
    load_cartesian_basis_representation(path)

Read `path` as a first-pass Cartesian JLD2 bundle and return the reconstructed
public `CartesianBasisRepresentation3D`.
"""
function load_cartesian_basis_representation(path::AbstractString)
    return read_cartesian_basis_bundle(path).basis
end

function _cartesian_require_exact_disk_bundle_support(
    bundle::CartesianBasisBundle3D,
    operation::AbstractString,
)
    if _cartesian_supports_exact_hybrid_overlap(bundle.basis)
        return bundle
    end
    throw(
        ArgumentError(
            "$(operation) does not yet support Cartesian bundle $(repr(bundle.path)) because its basis subtree does not carry a supported exact raw-space identity for parent kind :$(bundle.basis.metadata.parent_kind)",
        ),
    )
end

function cross_overlap(
    left::CartesianBasisBundle3D,
    right::CartesianBasisBundle3D,
)
    _cartesian_require_exact_disk_bundle_support(left, "disk-level Cartesian cross overlap")
    _cartesian_require_exact_disk_bundle_support(right, "disk-level Cartesian cross overlap")
    return cross_overlap(left.basis, right.basis)
end

function cross_overlap(
    left_path::AbstractString,
    right_path::AbstractString,
)
    return cross_overlap(
        read_cartesian_basis_bundle(left_path),
        read_cartesian_basis_bundle(right_path),
    )
end

function basis_projector(
    source::CartesianBasisBundle3D,
    target::CartesianBasisBundle3D,
)
    _cartesian_require_exact_disk_bundle_support(source, "disk-level Cartesian basis projector")
    _cartesian_require_exact_disk_bundle_support(target, "disk-level Cartesian basis projector")
    return GaussletBases._cartesian_basis_projector_with_stage_timings(
        source.basis,
        target.basis,
    ).projector
end

function basis_projector(
    source_path::AbstractString,
    target_path::AbstractString,
)
    return basis_projector(
        read_cartesian_basis_bundle(source_path),
        read_cartesian_basis_bundle(target_path),
    )
end

"""
    transfer_orbitals(source_coefficients, source_bundle, target_bundle; materialize_projector = true)
    transfer_orbitals(source_coefficients, source_path, target_path; materialize_projector = true)

Transfer coefficients between Cartesian basis bundles using the same final-basis
contract as representation-level transfer: `C_B = S_BA * C_A`.

The default `materialize_projector = true` preserves the explicit-projector
route and returns a result carrying the materialized projector. When
`materialize_projector = false` is supported by the loaded source/target basis
representations, the transfer is applied directly to the coefficient block and
the returned result may have `projector === nothing`.
"""
function transfer_orbitals(
    source_coefficients::AbstractVecOrMat{<:Real},
    source::CartesianBasisBundle3D,
    target::CartesianBasisBundle3D;
    materialize_projector::Bool = true,
)
    _cartesian_require_exact_disk_bundle_support(source, "disk-level Cartesian orbital transfer")
    _cartesian_require_exact_disk_bundle_support(target, "disk-level Cartesian orbital transfer")
    if !materialize_projector
        return transfer_orbitals(
            source_coefficients,
            source.basis,
            target.basis;
            materialize_projector = false,
        )
    end
    projector = basis_projector(source, target)
    return transfer_orbitals(source_coefficients, projector)
end

function transfer_orbitals(
    source_coefficients::AbstractVecOrMat{<:Real},
    source_path::AbstractString,
    target_path::AbstractString;
    materialize_projector::Bool = true,
)
    return transfer_orbitals(
        source_coefficients,
        read_cartesian_basis_bundle(source_path),
        read_cartesian_basis_bundle(target_path);
        materialize_projector = materialize_projector,
    )
end
