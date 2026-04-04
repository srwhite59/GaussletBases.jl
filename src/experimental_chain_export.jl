function _experimental_chain_nuclear_coordinate_matrix(
    nuclei::AbstractVector{<:NTuple{3,<:Real}},
)
    matrix = zeros(Float64, length(nuclei), 3)
    for (i, nucleus) in pairs(nuclei)
        matrix[i, 1] = Float64(nucleus[1])
        matrix[i, 2] = Float64(nucleus[2])
        matrix[i, 3] = Float64(nucleus[3])
    end
    return matrix
end

function _experimental_square_lattice_uniform_spacing(
    values::AbstractVector{<:Real};
    atol::Float64 = 1.0e-12,
    rtol::Float64 = 1.0e-10,
)
    length(values) <= 1 && return nothing
    diffs = diff(Float64[Float64(value) for value in values])
    reference = diffs[1]
    all(isapprox(delta, reference; atol = atol, rtol = rtol) for delta in diffs) || return nothing
    return reference
end

function _experimental_square_lattice_spacing(
    basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D;
    atol::Float64 = 1.0e-12,
    rtol::Float64 = 1.0e-10,
)
    x_spacing = _experimental_square_lattice_uniform_spacing(basis.x_coordinates; atol = atol, rtol = rtol)
    y_spacing = _experimental_square_lattice_uniform_spacing(basis.y_coordinates; atol = atol, rtol = rtol)
    x_spacing === nothing && return nothing
    y_spacing === nothing && return nothing
    isapprox(x_spacing, y_spacing; atol = atol, rtol = rtol) || return nothing
    return x_spacing
end

function _experimental_square_lattice_geometry_report_text(
    path::ExperimentalAxisAlignedHomonuclearSquareLatticeNestedQWPath,
)
    return join(_square_lattice_nested_geometry_report_lines(path.source), "\n") * "\n"
end

function _experimental_square_lattice_geometry_checksum(
    path::ExperimentalAxisAlignedHomonuclearSquareLatticeNestedQWPath,
)
    return bytes2hex(SHA.sha1(_experimental_square_lattice_geometry_report_text(path)))
end

function _experimental_square_lattice_accepted_root_candidate(root)
    accepted_index = root.accepted_candidate_index
    accepted_index === nothing && return nothing
    return root.candidate_summaries[accepted_index]
end

function _experimental_square_lattice_child_planar_count_matrix(
    candidate,
)
    candidate === nothing && return zeros(Int, 0, 2)
    counts = candidate.child_planar_counts
    matrix = zeros(Int, length(counts), 2)
    for (i, count) in pairs(counts)
        matrix[i, 1] = Int(count[1])
        matrix[i, 2] = Int(count[2])
    end
    return matrix
end

function _experimental_square_lattice_child_width_matrix(
    candidate,
)
    candidate === nothing && return zeros(Float64, 0, 3)
    widths = candidate.child_physical_widths
    matrix = zeros(Float64, length(widths), 3)
    for (i, width) in pairs(widths)
        matrix[i, 1] = Float64(width[1])
        matrix[i, 2] = Float64(width[2])
        matrix[i, 3] = Float64(width[3])
    end
    return matrix
end

function _experimental_square_lattice_export_meta_values(
    path::ExperimentalAxisAlignedHomonuclearSquareLatticeNestedQWPath;
    meta = nothing,
)
    diagnostics = path.diagnostics
    root = diagnostics.root_node
    spacing = _experimental_square_lattice_spacing(path.basis)
    meta_values = _normalize_meta_dict(meta)
    merge!(
        meta_values,
        Dict{String,Any}(
            "producer" => "GaussletBases.write_experimental_homonuclear_square_lattice_nested_dense_jld2",
            "manifest/producer/package" => "GaussletBases",
            "manifest/producer/version" => string(Base.pkgversion(@__MODULE__)),
            "manifest/producer/entrypoint" =>
                "GaussletBases.write_experimental_homonuclear_square_lattice_nested_dense_jld2",
            "manifest/producer/object_type" =>
                "ExperimentalAxisAlignedHomonuclearSquareLatticeNestedQWPath",
            "manifest/contract/status" => "experimental",
            "manifest/contract/scope" => "producer_only",
            "manifest/contract/line" => "homonuclear_square_lattice",
            "manifest/source/branch" => "experimental_homonuclear_square_lattice_nested_qw",
            "manifest/source/homonuclear" => true,
            "manifest/source/lattice_size" => path.basis.lattice_size,
            "manifest/source/natoms" => length(path.basis.nuclei),
            "manifest/source/coordinate_provenance" =>
                spacing === nothing ? "explicit_axis_coordinates" : "uniform_square_spacing",
            "manifest/source/min_in_plane_aspect_ratio" => path.min_in_plane_aspect_ratio,
            "manifest/source/root_did_split" => root.did_split,
            "manifest/source/root_child_count" => root.child_count,
            "manifest/source/leaf_count" => diagnostics.leaf_count,
            "manifest/source/fixed_dimension" => diagnostics.fixed_dimension,
            "manifest/source/residual_sector_empty" => path.operators.residual_count == 0,
            "manifest/source/nested_geometry_checksum" =>
                _experimental_square_lattice_geometry_checksum(path),
        ),
    )
    spacing === nothing || (meta_values["manifest/source/lattice_spacing"] = spacing)
    return meta_values
end

function _experimental_chain_geometry_report_text(
    path::ExperimentalBondAlignedHomonuclearChainNestedQWPath,
)
    return join(_chain_nested_geometry_report_lines(path.source), "\n") * "\n"
end

function _experimental_chain_geometry_checksum(
    path::ExperimentalBondAlignedHomonuclearChainNestedQWPath,
)
    return bytes2hex(SHA.sha1(_experimental_chain_geometry_report_text(path)))
end

function _experimental_chain_export_meta_values(
    path::ExperimentalBondAlignedHomonuclearChainNestedQWPath;
    meta = nothing,
)
    diagnostics = path.diagnostics
    root = diagnostics.root_node
    meta_values = _normalize_meta_dict(meta)
    merge!(
        meta_values,
        Dict{String,Any}(
            "producer" => "GaussletBases.write_experimental_homonuclear_chain_nested_dense_jld2",
            "manifest/producer/package" => "GaussletBases",
            "manifest/producer/version" => string(Base.pkgversion(@__MODULE__)),
            "manifest/producer/entrypoint" =>
                "GaussletBases.write_experimental_homonuclear_chain_nested_dense_jld2",
            "manifest/producer/object_type" =>
                "ExperimentalBondAlignedHomonuclearChainNestedQWPath",
            "manifest/contract/status" => "experimental",
            "manifest/contract/scope" => "producer_only",
            "manifest/contract/line" => "homonuclear_linear_chain",
            "manifest/source/branch" => "experimental_homonuclear_chain_nested_qw",
            "manifest/source/homonuclear" => true,
            "manifest/source/chain_axis" => string(path.basis.chain_axis),
            "manifest/source/natoms" => length(path.basis.nuclei),
            "manifest/source/odd_chain_policy" => string(path.odd_chain_policy),
            "manifest/source/root_did_split" => root.did_split,
            "manifest/source/root_child_count" => root.child_count,
            "manifest/source/leaf_count" => diagnostics.leaf_count,
            "manifest/source/fixed_dimension" => diagnostics.fixed_dimension,
            "manifest/source/residual_sector_empty" => path.operators.residual_count == 0,
            "manifest/source/nested_geometry_checksum" =>
                _experimental_chain_geometry_checksum(path),
        ),
    )
    return meta_values
end

"""
    experimental_homonuclear_chain_nested_dense_payload(path; meta=nothing)

Build the narrow experimental dense export payload for one
`ExperimentalBondAlignedHomonuclearChainNestedQWPath` without writing a file.

This payload is intentionally policy-explicit and producer-facing. It records
the chain geometry, the odd-chain policy used to build the nested path, the
split-tree checksum, and the dense `S`, `H1`, and `Vee` objects needed by
downstream experimental HF/DMRG consumers.
"""
function experimental_homonuclear_chain_nested_dense_payload(
    path::ExperimentalBondAlignedHomonuclearChainNestedQWPath;
    meta = nothing,
)
    diagnostics = path.diagnostics
    root = diagnostics.root_node
    operators = path.operators
    fixed_block = path.fixed_block
    orbital_labels = String[orbital.label for orbital in orbitals(operators)]
    geometry_report_text = _experimental_chain_geometry_report_text(path)
    bridge_meta = Dict{String,Any}(
        "format" => "experimental_homonuclear_chain_nested_dense_v1",
        "version" => 1,
        "experimental" => true,
        "producer_contract_scope" => "experimental_producer_only",
        "site_type" => "Electron",
        "interaction_model" => "density_density",
        "model_detail" => "nested_fixed_block_qw",
        "source_branch" => "experimental_homonuclear_chain_nested_qw",
        "onebody_key" => "H1",
        "interaction_key" => "Vee",
        "overlap_key" => "S",
        "norb" => size(operators.one_body_hamiltonian, 1),
        "natoms" => length(path.basis.nuclei),
        "chain_axis" => string(path.basis.chain_axis),
        "homonuclear" => true,
        "odd_chain_policy" => string(path.odd_chain_policy),
        "odd_chain_policy_is_default_reference" => path.odd_chain_policy == :strict_current,
        "residual_sector_empty" => operators.residual_count == 0,
        "gausslet_count" => operators.gausslet_count,
        "residual_count" => operators.residual_count,
        "fixed_dimension" => size(fixed_block.overlap, 1),
        "leaf_count" => diagnostics.leaf_count,
        "root_did_split" => root.did_split,
        "root_has_accepted_candidate" => !isnothing(root.accepted_candidate_index),
        "root_accepted_candidate_index" => something(root.accepted_candidate_index, 0),
        "root_child_count" => root.child_count,
        "root_local_resolution_warning" => root.local_resolution_warning,
        "nested_geometry_checksum" => _experimental_chain_geometry_checksum(path),
        "orbital_label_count" => length(orbital_labels),
    )
    payload = Dict{String,Any}(
        "S" => Matrix{Float64}(operators.overlap),
        "H1" => Matrix{Float64}(operators.one_body_hamiltonian),
        "Vee" => Matrix{Float64}(operators.interaction_matrix),
        "basis_centers" => Matrix{Float64}(fixed_block.fixed_centers),
        "nuclear_coordinates_xyz" => _experimental_chain_nuclear_coordinate_matrix(path.basis.nuclei),
        "nuclear_charges" => copy(path.nuclear_charges),
        "orbital_labels" => orbital_labels,
        "raw_to_final" => Matrix{Float64}(operators.raw_to_final),
        "residual_centers" => Matrix{Float64}(operators.residual_centers),
        "residual_widths" => Matrix{Float64}(operators.residual_widths),
        "geometry_report_text" => geometry_report_text,
    )
    return (
        payload = payload,
        bridge_meta = bridge_meta,
        meta_values = _experimental_chain_export_meta_values(path; meta = meta),
    )
end

"""
    write_experimental_homonuclear_chain_nested_dense_jld2(path, chain_path; meta=nothing)

Write one experimental dense producer-side artifact for the homonuclear nested
chain line.

This writer stores the dense `S`, `H1`, and `Vee` operators together with the
explicit odd-chain policy, the chain geometry, and a stable split-tree report
checksum so downstream experimental consumers can audit exactly which branch
they received.
"""
function write_experimental_homonuclear_chain_nested_dense_jld2(
    path::AbstractString,
    chain_path::ExperimentalBondAlignedHomonuclearChainNestedQWPath;
    meta = nothing,
)
    data = experimental_homonuclear_chain_nested_dense_payload(chain_path; meta = meta)
    jldopen(path, "w") do file
        for (key, value) in data.payload
            file[key] = value
        end
        _write_prefixed_values!(file, "bridge", data.bridge_meta)
        _write_prefixed_values!(file, "meta", data.meta_values)
    end
    return path
end

"""
    experimental_homonuclear_square_lattice_nested_dense_payload(path; meta=nothing)

Build the narrow experimental dense export payload for one
`ExperimentalAxisAlignedHomonuclearSquareLatticeNestedQWPath` without writing a
file.

This payload is intentionally planar-policy-explicit and producer-facing. It
records the lattice geometry, the accepted planar root split metadata, the
in-plane aspect threshold, the split-tree checksum, and the dense `S`, `H1`,
and `Vee` objects needed by downstream experimental HF/DMRG consumers.
"""
function experimental_homonuclear_square_lattice_nested_dense_payload(
    path::ExperimentalAxisAlignedHomonuclearSquareLatticeNestedQWPath;
    meta = nothing,
)
    diagnostics = path.diagnostics
    root = diagnostics.root_node
    accepted = _experimental_square_lattice_accepted_root_candidate(root)
    operators = path.operators
    fixed_block = path.fixed_block
    orbital_labels = String[orbital.label for orbital in orbitals(operators)]
    geometry_report_text = _experimental_square_lattice_geometry_report_text(path)
    spacing = _experimental_square_lattice_spacing(path.basis)
    bridge_meta = Dict{String,Any}(
        "format" => "experimental_homonuclear_square_lattice_nested_dense_v1",
        "version" => 1,
        "experimental" => true,
        "producer_contract_scope" => "experimental_producer_only",
        "site_type" => "Electron",
        "interaction_model" => "density_density",
        "model_detail" => "nested_fixed_block_qw",
        "source_branch" => "experimental_homonuclear_square_lattice_nested_qw",
        "onebody_key" => "H1",
        "interaction_key" => "Vee",
        "overlap_key" => "S",
        "norb" => size(operators.one_body_hamiltonian, 1),
        "natoms" => length(path.basis.nuclei),
        "lattice_size" => path.basis.lattice_size,
        "homonuclear" => true,
        "coordinate_provenance" =>
            spacing === nothing ? "explicit_axis_coordinates" : "uniform_square_spacing",
        "residual_sector_empty" => operators.residual_count == 0,
        "gausslet_count" => operators.gausslet_count,
        "residual_count" => operators.residual_count,
        "fixed_dimension" => size(fixed_block.overlap, 1),
        "leaf_count" => diagnostics.leaf_count,
        "root_did_split" => root.did_split,
        "root_has_accepted_candidate" => !isnothing(root.accepted_candidate_index),
        "root_accepted_candidate_index" => something(root.accepted_candidate_index, 0),
        "root_accepted_split_family" => accepted === nothing ? "none" : string(accepted.split_family),
        "root_accepted_split_axis" => accepted === nothing ? "none" : string(accepted.split_axis),
        "root_child_count" => root.child_count,
        "root_local_resolution_warning" => root.local_resolution_warning,
        "min_in_plane_aspect_ratio" => path.min_in_plane_aspect_ratio,
        "nested_geometry_checksum" => _experimental_square_lattice_geometry_checksum(path),
        "orbital_label_count" => length(orbital_labels),
    )
    spacing === nothing || (bridge_meta["lattice_spacing"] = spacing)
    payload = Dict{String,Any}(
        "S" => Matrix{Float64}(operators.overlap),
        "H1" => Matrix{Float64}(operators.one_body_hamiltonian),
        "Vee" => Matrix{Float64}(operators.interaction_matrix),
        "basis_centers" => Matrix{Float64}(fixed_block.fixed_centers),
        "nuclear_coordinates_xyz" => _experimental_chain_nuclear_coordinate_matrix(path.basis.nuclei),
        "nuclear_charges" => copy(path.nuclear_charges),
        "lattice_x_coordinates" => Float64[Float64(value) for value in path.basis.x_coordinates],
        "lattice_y_coordinates" => Float64[Float64(value) for value in path.basis.y_coordinates],
        "orbital_labels" => orbital_labels,
        "raw_to_final" => Matrix{Float64}(operators.raw_to_final),
        "residual_centers" => Matrix{Float64}(operators.residual_centers),
        "residual_widths" => Matrix{Float64}(operators.residual_widths),
        "root_accepted_child_planar_counts" =>
            _experimental_square_lattice_child_planar_count_matrix(accepted),
        "root_accepted_child_physical_widths" =>
            _experimental_square_lattice_child_width_matrix(accepted),
        "root_accepted_child_in_plane_aspect_ratios" =>
            accepted === nothing ? Float64[] : Float64[Float64(value) for value in accepted.child_in_plane_aspect_ratios],
        "root_accepted_split_values" =>
            accepted === nothing ? Float64[] : Float64[Float64(value) for value in accepted.split_values],
        "root_accepted_split_indices" =>
            accepted === nothing ? Int[] : Int[Int(value) for value in accepted.split_indices],
        "geometry_report_text" => geometry_report_text,
    )
    return (
        payload = payload,
        bridge_meta = bridge_meta,
        meta_values = _experimental_square_lattice_export_meta_values(path; meta = meta),
    )
end

"""
    write_experimental_homonuclear_square_lattice_nested_dense_jld2(path, lattice_path; meta=nothing)

Write one experimental dense producer-side artifact for the homonuclear nested
square-lattice line.

This writer stores the dense `S`, `H1`, and `Vee` operators together with the
explicit planar-policy provenance, the lattice geometry, and a stable split-tree
report checksum so downstream experimental consumers can audit exactly which
branch they received.
"""
function write_experimental_homonuclear_square_lattice_nested_dense_jld2(
    path::AbstractString,
    lattice_path::ExperimentalAxisAlignedHomonuclearSquareLatticeNestedQWPath;
    meta = nothing,
)
    data = experimental_homonuclear_square_lattice_nested_dense_payload(lattice_path; meta = meta)
    jldopen(path, "w") do file
        for (key, value) in data.payload
            file[key] = value
        end
        _write_prefixed_values!(file, "bridge", data.bridge_meta)
        _write_prefixed_values!(file, "meta", data.meta_values)
    end
    return path
end
