"""
    ExperimentalBondAlignedHomonuclearChainNestedQWPath

First experimental nested fixed-block ordinary-QW chain consumer payload.

This keeps the operator object together with the chain nested geometry source
and the explicit odd-chain policy diagnostics used to build it.
"""
struct ExperimentalBondAlignedHomonuclearChainNestedQWPath{B,S,F,O,D}
    basis::B
    source::S
    fixed_block::F
    operators::O
    diagnostics::D
    nuclear_charges::Vector{Float64}
    odd_chain_policy::Symbol
end

"""
    ExperimentalAxisAlignedHomonuclearSquareLatticeNestedQWPath

First experimental nested fixed-block ordinary-QW square-lattice consumer
payload.

This keeps the operator object together with the exploratory planar split-tree
geometry source and the explicit in-plane aspect threshold used to build it.
"""
struct ExperimentalAxisAlignedHomonuclearSquareLatticeNestedQWPath{B,S,F,O,D}
    basis::B
    source::S
    fixed_block::F
    operators::O
    diagnostics::D
    nuclear_charges::Vector{Float64}
    min_in_plane_aspect_ratio::Float64
end

function _normalized_nested_source_frontend_context(
    basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    return _QWRGNestedSourceFrontendContext(
        basis,
        expansion,
        gausslet_backend,
        (
            nside = nside,
            min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
            term_coefficients = term_coefficients,
        ),
        (
            route_label = "axis-aligned homonuclear square-lattice nested fixed source",
            total_timing_label = nothing,
            axis_bundles_timing_label = nothing,
            source_assembly_timing_label = nothing,
        ),
    )
end

function _normalized_nested_source_frontend_context(
    basis::BondAlignedHomonuclearChainQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    odd_chain_policy::Symbol = :strict_current,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    return _QWRGNestedSourceFrontendContext(
        basis,
        expansion,
        gausslet_backend,
        (
            nside = nside,
            min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
            odd_chain_policy = odd_chain_policy,
            term_coefficients = term_coefficients,
        ),
        (
            route_label = "bond-aligned homonuclear chain nested fixed source",
            total_timing_label = nothing,
            axis_bundles_timing_label = nothing,
            source_assembly_timing_label = nothing,
        ),
    )
end

function _nested_source_from_frontend_context(
    context::_QWRGNestedSourceFrontendContext{<:AxisAlignedHomonuclearSquareLatticeQWBasis3D},
    bundles::_CartesianNestedAxisBundles3D,
    term_coefficients::AbstractVector{Float64},
)
    options = context.build_options
    return _nested_axis_aligned_homonuclear_square_lattice_source(
        context.basis,
        bundles;
        nside = options.nside,
        min_in_plane_aspect_ratio = options.min_in_plane_aspect_ratio,
        term_coefficients = term_coefficients,
    )
end

function _nested_source_from_frontend_context(
    context::_QWRGNestedSourceFrontendContext{<:BondAlignedHomonuclearChainQWBasis3D},
    bundles::_CartesianNestedAxisBundles3D,
    term_coefficients::AbstractVector{Float64},
)
    options = context.build_options
    return _nested_bond_aligned_homonuclear_chain_source(
        context.basis,
        bundles;
        chain_axis = context.basis.chain_axis,
        nside = options.nside,
        min_parallel_to_transverse_ratio = options.min_parallel_to_transverse_ratio,
        odd_chain_policy = options.odd_chain_policy,
        term_coefficients = term_coefficients,
    )
end

function _axis_aligned_homonuclear_square_lattice_nested_fixed_source(
    basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
)
    context = _normalized_nested_source_frontend_context(
        basis;
        expansion = expansion,
        gausslet_backend = gausslet_backend,
        nside = nside,
        min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
    )
    return _nested_source_frontend_source(context)
end

function _axis_aligned_homonuclear_square_lattice_nested_fixed_block(
    basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
)
    context = _normalized_nested_source_frontend_context(
        basis;
        expansion = expansion,
        gausslet_backend = gausslet_backend,
        nside = nside,
        min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
    )
    return _nested_source_frontend_fixed_block(context)
end

function _square_lattice_nested_geometry_report_lines(
    source::_CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D,
)
    retention = source.shell_retention_contract
    audit = _nested_source_contract_audit(source)
    lines = String[
        "# GaussletBases axis-aligned homonuclear square-lattice nested geometry report",
        "# lattice_size = $(source.basis.lattice_size)",
        "# natoms = $(length(source.basis.nuclei))",
        "# nside = $(retention.nside)",
        "# retain_xy = $(retention.retain_xy)",
        "# retain_xz = $(retention.retain_xz)",
        "# retain_yz = $(retention.retain_yz)",
        "# retain_x_edge = $(retention.retain_x_edge)",
        "# retain_y_edge = $(retention.retain_y_edge)",
        "# retain_z_edge = $(retention.retain_z_edge)",
        "# shell_increment = $(retention.shell_increment)",
        "# matches_nside_default = $(retention.matches_nside_default)",
        "# full_parent_working_box = $(audit.full_parent_working_box)",
        "# support_count = $(audit.support_count)",
        "# expected_support_count = $(audit.expected_support_count)",
        "# missing_row_count = $(audit.missing_row_count)",
        "# ownership_group_count_min = $(audit.ownership_group_count_min)",
        "# ownership_group_count_max = $(audit.ownership_group_count_max)",
        "# ownership_unowned_row_count = $(audit.ownership_unowned_row_count)",
        "# ownership_multi_owned_row_count = $(audit.ownership_multi_owned_row_count)",
        "# nleaf = $(length(source.leaf_sequences))",
        "# nfixed = $(size(source.sequence.coefficient_matrix, 2))",
    ]
    for node in _nested_square_lattice_collect_node_summaries(source.root_geometry)
        push!(lines, "")
        push!(lines, "[node $(node.node_label)]")
        push!(lines, "x_coordinate_range = $(node.x_coordinate_range)")
        push!(lines, "y_coordinate_range = $(node.y_coordinate_range)")
        push!(lines, "working_box = $(node.working_box)")
        push!(lines, "min_in_plane_aspect_ratio = $(node.min_in_plane_aspect_ratio)")
        push!(lines, "shared_shell_count = $(node.shared_shell_count)")
        push!(lines, "shared_shell_dimensions = $(node.shared_shell_dimensions)")
        for (index, provenance) in pairs(node.shared_shell_provenance)
            push!(lines, "shared_shell[$index].source_box = $(provenance.source_box)")
            push!(lines, "shared_shell[$index].next_inner_box = $(provenance.next_inner_box)")
            push!(lines, "shared_shell[$index].source_point_count = $(provenance.source_point_count)")
            push!(lines, "shared_shell[$index].retained_fixed_count = $(provenance.retained_fixed_count)")
        end
        push!(lines, "accepted_candidate_index = $(node.accepted_candidate_index)")
        push!(lines, "local_resolution_warning = $(node.local_resolution_warning)")
        push!(lines, "child_count = $(node.child_count)")
        push!(lines, "subtree_fixed_dimension = $(node.subtree_fixed_dimension)")
        for (index, candidate) in pairs(node.candidate_summaries)
            push!(lines, "candidate[$index].split_family = $(candidate.split_family)")
            push!(lines, "candidate[$index].split_axis = $(candidate.split_axis)")
            push!(lines, "candidate[$index].x_coordinate_ranges = $(candidate.x_coordinate_ranges)")
            push!(lines, "candidate[$index].y_coordinate_ranges = $(candidate.y_coordinate_ranges)")
            push!(lines, "candidate[$index].split_values = $(candidate.split_values)")
            push!(lines, "candidate[$index].split_indices = $(candidate.split_indices)")
            push!(lines, "candidate[$index].child_boxes = $(candidate.child_boxes)")
            push!(lines, "candidate[$index].child_planar_counts = $(candidate.child_planar_counts)")
            push!(lines, "candidate[$index].child_physical_widths = $(candidate.child_physical_widths)")
            push!(lines, "candidate[$index].child_in_plane_aspect_ratios = $(candidate.child_in_plane_aspect_ratios)")
            push!(lines, "candidate[$index].count_eligible = $(candidate.count_eligible)")
            push!(lines, "candidate[$index].shape_eligible = $(candidate.shape_eligible)")
            push!(lines, "candidate[$index].symmetry_preserving = $(candidate.symmetry_preserving)")
            push!(lines, "candidate[$index].did_split = $(candidate.did_split)")
            push!(lines, "candidate[$index].accepted = $(candidate.accepted)")
        end
    end
    return lines
end

function _axis_aligned_homonuclear_square_lattice_nested_geometry_diagnostics(
    source::_CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D,
)
    node_summaries = _nested_square_lattice_collect_node_summaries(source.root_geometry)
    common_contract = _nested_source_common_contract(source)
    return (
        source = source,
        root_node = _nested_square_lattice_node_summary(source.root_geometry),
        node_summaries = node_summaries,
        nside = source.shell_retention_contract.nside,
        retention_contract = source.shell_retention_contract,
        shared_shell_dimensions = common_contract.shared_shell_dimensions,
        shared_shell_provenance = common_contract.shared_shell_provenance,
        shared_shells_match_contract =
            all(
                shell -> shell.retained_fixed_count == source.shell_retention_contract.shell_increment,
                common_contract.shared_shell_provenance,
            ),
        contract_audit = common_contract.contract_audit,
        leaf_count = common_contract.leaf_count,
        fixed_dimension = common_contract.fixed_dimension,
    )
end

function _nested_source_geometry_diagnostics(
    source::_CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D,
)
    return _axis_aligned_homonuclear_square_lattice_nested_geometry_diagnostics(source)
end

function axis_aligned_homonuclear_square_lattice_nested_geometry_diagnostics(
    basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
)
    context = _normalized_nested_source_frontend_context(
        basis;
        expansion = expansion,
        gausslet_backend = gausslet_backend,
        nside = nside,
        min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
    )
    return _nested_source_frontend_geometry_diagnostics(context)
end

function write_axis_aligned_homonuclear_square_lattice_nested_geometry_report(
    path::AbstractString,
    basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D;
    kwargs...,
)
    diagnostics = axis_aligned_homonuclear_square_lattice_nested_geometry_diagnostics(
        basis;
        kwargs...,
    )
    mkpath(dirname(String(path)))
    open(path, "w") do io
        for line in _square_lattice_nested_geometry_report_lines(diagnostics.source)
            write(io, line, "\n")
        end
    end
    return diagnostics
end

function _nested_geometry_shared_shell_dimensions(
    node_summaries,
)
    dimensions = Int[]
    for node in node_summaries
        append!(dimensions, Int.(node.shared_shell_dimensions))
    end
    return dimensions
end

function _nested_geometry_shared_shell_provenance(
    node_summaries,
)
    provenance = _CartesianNestedShellLayerProvenance3D[]
    for node in node_summaries
        append!(provenance, node.shared_shell_provenance)
    end
    return provenance
end

function _nested_source_leaf_count(
    source::_CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D,
)
    return length(source.leaf_sequences)
end

function _nested_source_shared_shell_dimensions(
    source::_CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D,
)
    return _nested_geometry_shared_shell_dimensions(
        _nested_square_lattice_collect_node_summaries(source.root_geometry),
    )
end

function _nested_source_shared_shell_provenance(
    source::_CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D,
)
    return _nested_geometry_shared_shell_provenance(
        _nested_square_lattice_collect_node_summaries(source.root_geometry),
    )
end

function _nested_source_leaf_count(
    source::_CartesianNestedBondAlignedHomonuclearChainSource3D,
)
    return length(source.leaf_sequences)
end

function _nested_source_shared_shell_dimensions(
    source::_CartesianNestedBondAlignedHomonuclearChainSource3D,
)
    return _nested_geometry_shared_shell_dimensions(
        _nested_chain_collect_node_summaries(source.root_geometry),
    )
end

function _nested_source_shared_shell_provenance(
    source::_CartesianNestedBondAlignedHomonuclearChainSource3D,
)
    return _nested_geometry_shared_shell_provenance(
        _nested_chain_collect_node_summaries(source.root_geometry),
    )
end

function _bond_aligned_homonuclear_chain_nested_fixed_source(
    basis::BondAlignedHomonuclearChainQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    odd_chain_policy::Symbol = :strict_current,
)
    context = _normalized_nested_source_frontend_context(
        basis;
        expansion = expansion,
        gausslet_backend = gausslet_backend,
        nside = nside,
        min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        odd_chain_policy = odd_chain_policy,
    )
    return _nested_source_frontend_source(context)
end

function _bond_aligned_homonuclear_chain_nested_fixed_block(
    basis::BondAlignedHomonuclearChainQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    odd_chain_policy::Symbol = :strict_current,
)
    context = _normalized_nested_source_frontend_context(
        basis;
        expansion = expansion,
        gausslet_backend = gausslet_backend,
        nside = nside,
        min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        odd_chain_policy = odd_chain_policy,
    )
    return _nested_source_frontend_fixed_block(context)
end

function _chain_nested_geometry_report_lines(
    source::_CartesianNestedBondAlignedHomonuclearChainSource3D,
)
    retention = source.shell_retention_contract
    audit = _nested_source_contract_audit(source)
    lines = String[
        "# GaussletBases bond-aligned homonuclear chain nested geometry report",
        "# chain_axis = $(source.basis.chain_axis)",
        "# natoms = $(length(source.basis.nuclei))",
        "# nside = $(retention.nside)",
        "# retain_xy = $(retention.retain_xy)",
        "# retain_xz = $(retention.retain_xz)",
        "# retain_yz = $(retention.retain_yz)",
        "# retain_x_edge = $(retention.retain_x_edge)",
        "# retain_y_edge = $(retention.retain_y_edge)",
        "# retain_z_edge = $(retention.retain_z_edge)",
        "# shell_increment = $(retention.shell_increment)",
        "# matches_nside_default = $(retention.matches_nside_default)",
        "# full_parent_working_box = $(audit.full_parent_working_box)",
        "# support_count = $(audit.support_count)",
        "# expected_support_count = $(audit.expected_support_count)",
        "# missing_row_count = $(audit.missing_row_count)",
        "# ownership_group_count_min = $(audit.ownership_group_count_min)",
        "# ownership_group_count_max = $(audit.ownership_group_count_max)",
        "# ownership_unowned_row_count = $(audit.ownership_unowned_row_count)",
        "# ownership_multi_owned_row_count = $(audit.ownership_multi_owned_row_count)",
        "# nleaf = $(length(source.leaf_sequences))",
        "# nfixed = $(size(source.sequence.coefficient_matrix, 2))",
    ]
    for node in _nested_chain_collect_node_summaries(source.root_geometry)
        push!(lines, "")
        push!(lines, "[node $(node.node_label)]")
        push!(lines, "nucleus_range = $(node.nucleus_range)")
        push!(lines, "working_box = $(node.working_box)")
        push!(lines, "odd_chain_policy = $(node.odd_chain_policy)")
        push!(lines, "odd_policy.outer_parallel_count_min = $(node.odd_chain_policy_thresholds.outer_parallel_count_min)")
        push!(lines, "odd_policy.center_parallel_count_min = $(node.odd_chain_policy_thresholds.center_parallel_count_min)")
        push!(lines, "odd_policy.total_parallel_count_min = $(node.odd_chain_policy_thresholds.total_parallel_count_min)")
        push!(lines, "odd_policy.outer_parallel_to_transverse_ratio_min = $(node.odd_chain_policy_thresholds.outer_parallel_to_transverse_ratio_min)")
        push!(lines, "odd_policy.center_parallel_to_transverse_ratio_min = $(node.odd_chain_policy_thresholds.center_parallel_to_transverse_ratio_min)")
        push!(lines, "shared_shell_count = $(node.shared_shell_count)")
        push!(lines, "shared_shell_dimensions = $(node.shared_shell_dimensions)")
        for (index, provenance) in pairs(node.shared_shell_provenance)
            push!(lines, "shared_shell[$index].source_box = $(provenance.source_box)")
            push!(lines, "shared_shell[$index].next_inner_box = $(provenance.next_inner_box)")
            push!(lines, "shared_shell[$index].source_point_count = $(provenance.source_point_count)")
            push!(lines, "shared_shell[$index].retained_fixed_count = $(provenance.retained_fixed_count)")
        end
        push!(lines, "accepted_candidate_index = $(node.accepted_candidate_index)")
        push!(lines, "local_resolution_warning = $(node.local_resolution_warning)")
        push!(lines, "child_count = $(node.child_count)")
        push!(lines, "subtree_fixed_dimension = $(node.subtree_fixed_dimension)")
        for (index, candidate) in pairs(node.candidate_summaries)
            push!(lines, "candidate[$index].split_kind = $(candidate.split_kind)")
            push!(lines, "candidate[$index].nucleus_ranges = $(candidate.nucleus_ranges)")
            push!(lines, "candidate[$index].midpoint_values = $(candidate.midpoint_values)")
            push!(lines, "candidate[$index].split_indices = $(candidate.split_indices)")
            push!(lines, "candidate[$index].child_boxes = $(candidate.child_boxes)")
            push!(lines, "candidate[$index].child_parallel_counts = $(candidate.child_parallel_counts)")
            push!(lines, "candidate[$index].child_parallel_to_transverse_ratios = $(candidate.child_parallel_to_transverse_ratios)")
            push!(lines, "candidate[$index].count_eligible = $(candidate.count_eligible)")
            push!(lines, "candidate[$index].shape_eligible = $(candidate.shape_eligible)")
            push!(lines, "candidate[$index].did_split = $(candidate.did_split)")
            push!(lines, "candidate[$index].accepted = $(candidate.accepted)")
        end
    end
    return lines
end

function _bond_aligned_homonuclear_chain_nested_geometry_diagnostics(
    source::_CartesianNestedBondAlignedHomonuclearChainSource3D,
)
    node_summaries = _nested_chain_collect_node_summaries(source.root_geometry)
    common_contract = _nested_source_common_contract(source)
    return (
        source = source,
        root_node = _nested_chain_node_summary(source.root_geometry),
        node_summaries = node_summaries,
        nside = source.shell_retention_contract.nside,
        retention_contract = source.shell_retention_contract,
        shared_shell_dimensions = common_contract.shared_shell_dimensions,
        shared_shell_provenance = common_contract.shared_shell_provenance,
        shared_shells_match_contract =
            all(
                shell -> shell.retained_fixed_count == source.shell_retention_contract.shell_increment,
                common_contract.shared_shell_provenance,
            ),
        contract_audit = common_contract.contract_audit,
        leaf_count = common_contract.leaf_count,
        fixed_dimension = common_contract.fixed_dimension,
    )
end

function _nested_source_geometry_diagnostics(
    source::_CartesianNestedBondAlignedHomonuclearChainSource3D,
)
    return _bond_aligned_homonuclear_chain_nested_geometry_diagnostics(source)
end

function bond_aligned_homonuclear_chain_nested_geometry_diagnostics(
    basis::BondAlignedHomonuclearChainQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    odd_chain_policy::Symbol = :strict_current,
)
    context = _normalized_nested_source_frontend_context(
        basis;
        expansion = expansion,
        gausslet_backend = gausslet_backend,
        nside = nside,
        min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        odd_chain_policy = odd_chain_policy,
    )
    return _nested_source_frontend_geometry_diagnostics(context)
end

function write_bond_aligned_homonuclear_chain_nested_geometry_report(
    path::AbstractString,
    basis::BondAlignedHomonuclearChainQWBasis3D;
    kwargs...,
)
    diagnostics = bond_aligned_homonuclear_chain_nested_geometry_diagnostics(
        basis;
        kwargs...,
    )
    mkpath(dirname(String(path)))
    open(path, "w") do io
        for line in _chain_nested_geometry_report_lines(diagnostics.source)
            write(io, line, "\n")
        end
    end
    return diagnostics
end

function Base.show(io::IO, path::ExperimentalBondAlignedHomonuclearChainNestedQWPath)
    print(
        io,
        "ExperimentalBondAlignedHomonuclearChainNestedQWPath(odd_chain_policy=:",
        path.odd_chain_policy,
        ", nfixed=",
        size(path.fixed_block.overlap, 1),
        ", nleaf=",
        path.diagnostics.leaf_count,
        ", did_split=",
        path.diagnostics.root_node.did_split,
        ")",
    )
end

function Base.show(io::IO, path::ExperimentalAxisAlignedHomonuclearSquareLatticeNestedQWPath)
    print(
        io,
        "ExperimentalAxisAlignedHomonuclearSquareLatticeNestedQWPath(min_in_plane_aspect_ratio=",
        path.min_in_plane_aspect_ratio,
        ", nfixed=",
        size(path.fixed_block.overlap, 1),
        ", nleaf=",
        path.diagnostics.leaf_count,
        ", did_split=",
        path.diagnostics.root_node.did_split,
        ")",
    )
end

function _experimental_nested_source_backed_path(
    context::_QWRGNestedSourceFrontendContext;
    nuclear_charges::AbstractVector{<:Real},
    nuclear_term_storage::Symbol = :auto,
    interaction_treatment::Symbol = :ggt_nearest,
    timing::Bool = false,
)
    source_fixed = _nested_source_frontend_fixed_block(context)
    source = source_fixed.source
    fixed_block = source_fixed.fixed_block
    common_contract = _nested_source_common_contract(source)
    size(fixed_block.overlap, 1) == common_contract.fixed_dimension || throw(
        ArgumentError(
            "experimental nested source-backed wrapper received fixed block inconsistent with the common nested-source contract",
        ),
    )
    diagnostics = _nested_source_geometry_diagnostics(source)
    diagnostics.fixed_dimension == common_contract.fixed_dimension || throw(
        ArgumentError(
            "experimental nested source-backed wrapper diagnostics drifted from the common nested-source contract on fixed_dimension",
        ),
    )
    diagnostics.shared_shell_dimensions == common_contract.shared_shell_dimensions || throw(
        ArgumentError(
            "experimental nested source-backed wrapper diagnostics drifted from the common nested-source contract on shared_shell_dimensions",
        ),
    )
    diagnostics.shared_shell_provenance == common_contract.shared_shell_provenance || throw(
        ArgumentError(
            "experimental nested source-backed wrapper diagnostics drifted from the common nested-source contract on shared_shell_provenance",
        ),
    )
    diagnostics.leaf_count == common_contract.leaf_count || throw(
        ArgumentError(
            "experimental nested source-backed wrapper diagnostics drifted from the common nested-source contract on leaf_count",
        ),
    )
    receipt = CartesianQWOperatorCarriedSpaces.cartesian_qw_operator_construction_receipt(
        fixed_block;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = context.expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = context.gausslet_backend,
        timing = timing,
    )
    receipt_diagnostics =
        CartesianQWOperatorCarriedSpaces.qw_operator_construction_receipt_diagnostics(receipt)
    receipt_diagnostics.delegated_to_existing_builder || throw(
        ArgumentError(
            "experimental nested source-backed wrapper receipt did not delegate to the existing QW builder",
        ),
    )
    receipt_diagnostics.source_sidecar_agree || throw(
        ArgumentError(
            "experimental nested source-backed wrapper receipt source/sidecar mismatch for fields " *
            join(string.(receipt_diagnostics.mismatch_fields), ", "),
        ),
    )
    receipt_diagnostics.new_hamiltonian_kernel_used == false || throw(
        ArgumentError(
            "experimental nested source-backed wrapper receipt unexpectedly used a new Hamiltonian kernel",
        ),
    )
    receipt_diagnostics.numerical_outputs_changed == false || throw(
        ArgumentError(
            "experimental nested source-backed wrapper receipt reported changed numerical outputs",
        ),
    )
    operators =
        CartesianQWOperatorCarriedSpaces.qw_operator_construction_receipt_operators(receipt)
    return _experimental_nested_source_backed_path(
        context,
        source,
        fixed_block,
        operators,
        diagnostics,
        nuclear_charges,
    )
end

function _experimental_nested_source_backed_path(
    context::_QWRGNestedSourceFrontendContext{<:BondAlignedHomonuclearChainQWBasis3D},
    source::_CartesianNestedBondAlignedHomonuclearChainSource3D,
    fixed_block::_NestedFixedBlock3D,
    operators,
    diagnostics,
    nuclear_charges::AbstractVector{<:Real},
)
    return ExperimentalBondAlignedHomonuclearChainNestedQWPath(
        context.basis,
        source,
        fixed_block,
        operators,
        diagnostics,
        Float64[Float64(value) for value in nuclear_charges],
        context.build_options.odd_chain_policy,
    )
end

function _experimental_nested_source_backed_path(
    context::_QWRGNestedSourceFrontendContext{<:AxisAlignedHomonuclearSquareLatticeQWBasis3D},
    source::_CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D,
    fixed_block::_NestedFixedBlock3D,
    operators,
    diagnostics,
    nuclear_charges::AbstractVector{<:Real},
)
    return ExperimentalAxisAlignedHomonuclearSquareLatticeNestedQWPath(
        context.basis,
        source,
        fixed_block,
        operators,
        diagnostics,
        Float64[Float64(value) for value in nuclear_charges],
        context.build_options.min_in_plane_aspect_ratio,
    )
end

function _experimental_high_order_cr_sp_supplement_fixture()
    return mktemp() do path, io
        write(
            io,
            "#BASIS SET: Cr repo-cr-sp\n" *
            "Cr    S\n" *
            "      4.0000000              1.0000000\n" *
            "Cr    P\n" *
            "      2.5000000              1.0000000\n" *
            "END\n",
        )
        close(io)
        legacy_atomic_gaussian_supplement("Cr", "repo-cr-sp"; lmax = 1, basisfile = path)
    end
end

function _experimental_high_order_axis_x2_matrix(axis_data::_ExperimentalHighOrderAxisData1D)
    isnothing(axis_data.pgdg_intermediate) && throw(
        ArgumentError("high-order stack QW diagnostic requires PGDG axis x2 metadata"),
    )
    return Matrix{Float64}(axis_data.pgdg_intermediate.x2)
end

function _experimental_high_order_parent_kinetic_3d(axis_data::_ExperimentalHighOrderAxisData1D)
    overlap = axis_data.overlap
    kinetic = axis_data.kinetic
    return _symmetrize_ida_matrix(
        kron(kinetic, kron(overlap, overlap)) +
        kron(overlap, kron(kinetic, overlap)) +
        kron(overlap, kron(overlap, kinetic)),
    )
end

function _experimental_high_order_parent_position_matrices_3d(axis_data::_ExperimentalHighOrderAxisData1D)
    overlap = axis_data.overlap
    position = axis_data.position
    return (
        x = _symmetrize_ida_matrix(kron(position, kron(overlap, overlap))),
        y = _symmetrize_ida_matrix(kron(overlap, kron(position, overlap))),
        z = _symmetrize_ida_matrix(kron(overlap, kron(overlap, position))),
    )
end

function _experimental_high_order_parent_x2_matrices_3d(axis_data::_ExperimentalHighOrderAxisData1D)
    overlap = axis_data.overlap
    x2 = _experimental_high_order_axis_x2_matrix(axis_data)
    return (
        x = _symmetrize_ida_matrix(kron(x2, kron(overlap, overlap))),
        y = _symmetrize_ida_matrix(kron(overlap, kron(x2, overlap))),
        z = _symmetrize_ida_matrix(kron(overlap, kron(overlap, x2))),
    )
end

function _experimental_high_order_parent_gaussian_sum_3d(
    axis_data::_ExperimentalHighOrderAxisData1D,
    expansion::CoulombGaussianExpansion,
)
    length(axis_data.gaussian_factors) == length(expansion) || throw(
        ArgumentError("high-order stack QW diagnostic requires one PGDG Gaussian factor per Coulomb-expansion term"),
    )
    terms = _term_tensor(axis_data.gaussian_factors)
    return _mapped_coulomb_expanded_symmetric_matrix(expansion.coefficients, terms, terms, terms)
end

function _experimental_high_order_parent_pair_sum_3d(
    axis_data::_ExperimentalHighOrderAxisData1D,
    expansion::CoulombGaussianExpansion,
)
    length(axis_data.pair_factors_1d) == length(expansion) || throw(
        ArgumentError("high-order stack QW diagnostic requires one PGDG pair factor per Coulomb-expansion term"),
    )
    terms = _term_tensor(axis_data.pair_factors_1d)
    return _mapped_coulomb_expanded_symmetric_matrix(expansion.coefficients, terms, terms, terms)
end

function _experimental_high_order_project_parent_matrix(
    coefficients::AbstractMatrix{<:Real},
    parent_matrix::AbstractMatrix{<:Real},
)
    coefficient_value = Matrix{Float64}(coefficients)
    return _symmetrize_ida_matrix(transpose(coefficient_value) * parent_matrix * coefficient_value)
end

function _experimental_high_order_active_support_indices(
    coefficients::AbstractMatrix{<:Real};
    tol::Real = 1.0e-14,
)
    coefficient_value = Matrix{Float64}(coefficients)
    tol_value = Float64(tol)
    return Int[
        row for row in axes(coefficient_value, 1) if maximum(abs, view(coefficient_value, row, :)) > tol_value
    ]
end

function _experimental_high_order_stack_fixed_block_for_atomic_qw(
    stack::ExperimentalHighOrderDosideStack3D,
    axis_data::_ExperimentalHighOrderAxisData1D,
    expansion::CoulombGaussianExpansion,
)
    stack.backend == :pgdg_localized_experimental || throw(
        ArgumentError("high-order stack QW diagnostic requires a PGDG stack; numerical-reference stacks are rejected"),
    )
    axis_data.backend == stack.backend || throw(
        ArgumentError("high-order stack QW diagnostic axis backend must match stack backend"),
    )
    axis_data.basis === stack.parent_basis || throw(
        ArgumentError("high-order stack QW diagnostic axis data must come from the stack parent basis"),
    )

    coefficients = Matrix{Float64}(stack.coefficient_matrix)
    parent_overlap = _experimental_high_order_parent_overlap_3d(axis_data)
    parent_kinetic = _experimental_high_order_parent_kinetic_3d(axis_data)
    positions = _experimental_high_order_parent_position_matrices_3d(axis_data)
    x2 = _experimental_high_order_parent_x2_matrices_3d(axis_data)
    parent_gaussian_sum = _experimental_high_order_parent_gaussian_sum_3d(axis_data, expansion)
    parent_pair_sum = _experimental_high_order_parent_pair_sum_3d(axis_data, expansion)
    fixed_position_x = _experimental_high_order_project_parent_matrix(coefficients, positions.x)
    fixed_position_y = _experimental_high_order_project_parent_matrix(coefficients, positions.y)
    fixed_position_z = _experimental_high_order_project_parent_matrix(coefficients, positions.z)

    return _NestedFixedBlock3D(
        stack.parent_basis,
        stack,
        stack.backend,
        coefficients,
        _experimental_high_order_active_support_indices(coefficients),
        _experimental_high_order_project_parent_matrix(coefficients, parent_overlap),
        _experimental_high_order_project_parent_matrix(coefficients, parent_kinetic),
        fixed_position_x,
        fixed_position_y,
        fixed_position_z,
        _experimental_high_order_project_parent_matrix(coefficients, x2.x),
        _experimental_high_order_project_parent_matrix(coefficients, x2.y),
        _experimental_high_order_project_parent_matrix(coefficients, x2.z),
        copy(stack.contracted_weights),
        _experimental_high_order_project_parent_matrix(coefficients, parent_gaussian_sum),
        _experimental_high_order_project_parent_matrix(coefficients, parent_pair_sum),
        Matrix{Float64}(hcat(diag(fixed_position_x), diag(fixed_position_y), diag(fixed_position_z))),
        Ref{Any}(nothing),
        Ref{Any}(nothing),
    )
end

function _experimental_high_order_stack_to_atomic_qw_operator_diagnostic(
    smoke = nothing;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    supplement::Union{Nothing,LegacyAtomicGaussianSupplement} = nothing,
    Z::Real = 24.0,
    interaction_treatment::Symbol = :mwg,
    gausslet_backend::Symbol = :pgdg_localized_experimental,
)
    gausslet_backend == :pgdg_localized_experimental || throw(
        ArgumentError("high-order stack QW diagnostic requires gausslet_backend = :pgdg_localized_experimental; no numerical-reference fallback is allowed"),
    )
    interaction_treatment == :mwg || throw(
        ArgumentError("high-order stack QW diagnostic currently validates only interaction_treatment = :mwg"),
    )

    smoke_timed = if isnothing(smoke)
        @timed _experimental_high_order_cr_map_count7_pgdg_smoke_diagnostic()
    else
        (value = smoke, time = 0.0, bytes = 0, gctime = 0.0, gcstats = nothing)
    end
    smoke_value = smoke_timed.value
    smoke_value.stack.backend == gausslet_backend || throw(
        ArgumentError("high-order stack QW diagnostic smoke stack backend must be :pgdg_localized_experimental"),
    )

    axis_timed = @timed _experimental_high_order_axis_data_1d(
        smoke_value.basis;
        backend = gausslet_backend,
        one_body_exponents = expansion.exponents,
        one_body_center = 0.0,
        include_pair_factors = true,
    )
    axis_data = axis_timed.value
    stack_timed = @timed _experimental_high_order_doside_stack_3d(
        smoke_value.basis;
        axis_data,
        backend = gausslet_backend,
        doside = smoke_value.stack.doside,
        sides = smoke_value.stack.sides,
    )
    stack = stack_timed.value
    fixed_timed = @timed _experimental_high_order_stack_fixed_block_for_atomic_qw(
        stack,
        axis_data,
        expansion,
    )
    fixed_block = fixed_timed.value
    supplement_value = isnothing(supplement) ? _experimental_high_order_cr_sp_supplement_fixture() : supplement
    supplement_dimension = length(_atomic_cartesian_shell_supplement_3d(supplement_value).orbitals)
    operators_timed = @timed ordinary_cartesian_qiu_white_operators(
        fixed_block,
        supplement_value;
        expansion,
        Z,
        interaction_treatment,
        gausslet_backend,
    )
    operators = operators_timed.value

    residual_widths_positive =
        operators.residual_count > 0 &&
        all(isfinite, operators.residual_widths) &&
        all(>(0.0), vec(operators.residual_widths))
    fixed_overlap_error = norm(fixed_block.overlap - I, Inf)
    operator_overlap_error = norm(operators.overlap - I, Inf)
    h_symmetry_error = norm(operators.one_body_hamiltonian - transpose(operators.one_body_hamiltonian), Inf)
    v_symmetry_error = norm(operators.interaction_matrix - transpose(operators.interaction_matrix), Inf)

    return (
        smoke = smoke_value,
        axis_data = axis_data,
        stack = stack,
        fixed_block = fixed_block,
        supplement = supplement_value,
        operators = operators,
        diagnostics = (
            route = :cr_map_count7_high_order_stack_to_atomic_qw_operator_smoke_diagnostic,
            classification = :diagnostic_only,
            operator_construction_scope = :smoke_plumbing_only_not_same_density_route_validation,
            interaction_treatment = interaction_treatment,
            parent_side = length(stack.parent_basis),
            parent_dimension = length(stack.parent_basis)^3,
            route_comparability = :smoke_only_not_ordinary_ns7,
            ordinary_ns7_comparable = false,
            ordinary_ns7_reference_parent_side = 27,
            ordinary_ns7_reference_parent_dimension = 27^3,
            high_order_retained_dimension = size(stack.coefficient_matrix, 2),
            support_count = length(fixed_block.support_indices),
            supplement_dimension = supplement_dimension,
            final_operator_dimension = size(operators.overlap, 1),
            residual_count = operators.residual_count,
            stack_backend = stack.backend,
            fixed_backend = fixed_block.gausslet_backend,
            adapter_backend = gausslet_backend,
            operator_backend = operators.gausslet_backend,
            contracted_weight_zeroish_count = count(abs.(stack.contracted_weights) .<= 1.0e-14),
            contracted_weight_negative_count = count(stack.contracted_weights .< -1.0e-14),
            contracted_weight_minimum = minimum(stack.contracted_weights),
            contracted_weight_maximum = maximum(stack.contracted_weights),
            fixed_overlap_error = fixed_overlap_error,
            operator_overlap_error = operator_overlap_error,
            h_symmetry_error = h_symmetry_error,
            v_symmetry_error = v_symmetry_error,
            residual_widths_finite_positive = residual_widths_positive,
            same_density_two_electron_evaluation =
                :blocked_smoke_only_not_route_validation,
            same_density_route_comparison = :blocked_smoke_only_not_ordinary_ns7,
            occupied_capture_gate = smoke_value.diagnostics.occupied_capture_gate,
            occupied_capture_status = smoke_value.diagnostics.occupied_capture_status,
            smallest_missing_interface = :ordinary_ns7_extent_high_order_route,
            timing_seconds = (
                smoke = smoke_timed.time,
                axis_data = axis_timed.time,
                stack = stack_timed.time,
                fixed_block = fixed_timed.time,
                operators = operators_timed.time,
                total = smoke_timed.time + axis_timed.time + stack_timed.time + fixed_timed.time + operators_timed.time,
            ),
            allocation_bytes = (
                smoke = smoke_timed.bytes,
                axis_data = axis_timed.bytes,
                stack = stack_timed.bytes,
                fixed_block = fixed_timed.bytes,
                operators = operators_timed.bytes,
                total = smoke_timed.bytes + axis_timed.bytes + stack_timed.bytes + fixed_timed.bytes + operators_timed.bytes,
            ),
        ),
    )
end

"""
    experimental_bond_aligned_homonuclear_chain_nested_qw_operators(
        basis::BondAlignedHomonuclearChainQWBasis3D;
        nuclear_charges = basis.nuclear_charges,
        expansion = coulomb_gaussian_expansion(doacc = false),
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        nside = 5,
        min_parallel_to_transverse_ratio = 0.4,
        odd_chain_policy = :central_ternary_relaxed,
        timing = false,
    )

Build the first experimental nested-chain ordinary-QW path on top of the
chain nested fixed-block geometry.

For odd chains, this milestone uses `odd_chain_policy = :central_ternary_relaxed`
explicitly by default. The stricter `:strict_current` branch remains available
through the lower-level geometry/fixed-block diagnostics and is not replaced as
the conservative reference policy.
"""
function experimental_bond_aligned_homonuclear_chain_nested_qw_operators(
    basis::BondAlignedHomonuclearChainQWBasis3D;
    nuclear_charges::AbstractVector{<:Real} = basis.nuclear_charges,
    nuclear_term_storage::Symbol = :auto,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    odd_chain_policy::Symbol = :central_ternary_relaxed,
    timing::Bool = false,
)
    context = _normalized_nested_source_frontend_context(
        basis;
        expansion = expansion,
        gausslet_backend = gausslet_backend,
        nside = nside,
        min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        odd_chain_policy = odd_chain_policy,
    )
    return _experimental_nested_source_backed_path(
        context;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        interaction_treatment = interaction_treatment,
        timing = timing,
    )
end

"""
    experimental_axis_aligned_homonuclear_square_lattice_nested_qw_operators(
        basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D;
        nuclear_charges = basis.nuclear_charges,
        expansion = coulomb_gaussian_expansion(doacc = false),
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        nside = 5,
        min_in_plane_aspect_ratio = 0.15,
        timing = false,
    )

Build the first experimental nested square-lattice ordinary-QW path on top of
the planar nested fixed-block geometry.

This remains explicitly exploratory. In particular, the center-strip aspect
threshold is carried as an exposed policy knob rather than a hidden settled
contract.
"""
function experimental_axis_aligned_homonuclear_square_lattice_nested_qw_operators(
    basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D;
    nuclear_charges::AbstractVector{<:Real} = basis.nuclear_charges,
    nuclear_term_storage::Symbol = :auto,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
    timing::Bool = false,
)
    context = _normalized_nested_source_frontend_context(
        basis;
        expansion = expansion,
        gausslet_backend = gausslet_backend,
        nside = nside,
        min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
    )
    return _experimental_nested_source_backed_path(
        context;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        interaction_treatment = interaction_treatment,
        timing = timing,
    )
end
