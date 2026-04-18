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

function _axis_aligned_homonuclear_square_lattice_nested_fixed_source(
    basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
)
    gausslet_backend == :numerical_reference || throw(
        ArgumentError("axis-aligned homonuclear square-lattice nested fixed source currently supports only gausslet_backend = :numerical_reference"),
    )
    bundles = _qwrg_bond_aligned_axis_bundles(
        basis,
        expansion;
        gausslet_backend = gausslet_backend,
    )
    term_coefficients = _resolved_nested_term_coefficients(expansion, nothing)
    return _nested_axis_aligned_homonuclear_square_lattice_source(
        basis,
        bundles;
        nside = nside,
        min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
        term_coefficients = term_coefficients,
    )
end

function _axis_aligned_homonuclear_square_lattice_nested_fixed_block(
    basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
)
    source = _axis_aligned_homonuclear_square_lattice_nested_fixed_source(
        basis;
        expansion = expansion,
        gausslet_backend = gausslet_backend,
        nside = nside,
        min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
    )
    return (
        source = source,
        fixed_block = _nested_fixed_block(source),
    )
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
    shared_shell_dimensions = _nested_geometry_shared_shell_dimensions(node_summaries)
    return (
        source = source,
        root_node = _nested_square_lattice_node_summary(source.root_geometry),
        node_summaries = node_summaries,
        nside = source.shell_retention_contract.nside,
        retention_contract = source.shell_retention_contract,
        shared_shell_dimensions = shared_shell_dimensions,
        shared_shells_match_contract =
            all(==(source.shell_retention_contract.shell_increment), shared_shell_dimensions),
        contract_audit = _nested_source_contract_audit(source),
        leaf_count = length(source.leaf_sequences),
        fixed_dimension = size(source.sequence.coefficient_matrix, 2),
    )
end

function axis_aligned_homonuclear_square_lattice_nested_geometry_diagnostics(
    basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
)
    source = _axis_aligned_homonuclear_square_lattice_nested_fixed_source(
        basis;
        expansion = expansion,
        gausslet_backend = gausslet_backend,
        nside = nside,
        min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
    )
    return _axis_aligned_homonuclear_square_lattice_nested_geometry_diagnostics(source)
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

function _bond_aligned_homonuclear_chain_nested_fixed_source(
    basis::BondAlignedHomonuclearChainQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    odd_chain_policy::Symbol = :strict_current,
)
    gausslet_backend == :numerical_reference || throw(
        ArgumentError("bond-aligned homonuclear chain nested fixed source currently supports only gausslet_backend = :numerical_reference"),
    )
    bundles = _qwrg_bond_aligned_axis_bundles(
        basis,
        expansion;
        gausslet_backend = gausslet_backend,
    )
    term_coefficients = _resolved_nested_term_coefficients(expansion, nothing)
    return _nested_bond_aligned_homonuclear_chain_source(
        basis,
        bundles;
        chain_axis = basis.chain_axis,
        nside = nside,
        min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        odd_chain_policy = odd_chain_policy,
        term_coefficients = term_coefficients,
    )
end

function _bond_aligned_homonuclear_chain_nested_fixed_block(
    basis::BondAlignedHomonuclearChainQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    odd_chain_policy::Symbol = :strict_current,
)
    source = _bond_aligned_homonuclear_chain_nested_fixed_source(
        basis;
        expansion = expansion,
        gausslet_backend = gausslet_backend,
        nside = nside,
        min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        odd_chain_policy = odd_chain_policy,
    )
    return (
        source = source,
        fixed_block = _nested_fixed_block(source),
    )
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
    shared_shell_dimensions = _nested_geometry_shared_shell_dimensions(node_summaries)
    return (
        source = source,
        root_node = _nested_chain_node_summary(source.root_geometry),
        node_summaries = node_summaries,
        nside = source.shell_retention_contract.nside,
        retention_contract = source.shell_retention_contract,
        shared_shell_dimensions = shared_shell_dimensions,
        shared_shells_match_contract =
            all(==(source.shell_retention_contract.shell_increment), shared_shell_dimensions),
        contract_audit = _nested_source_contract_audit(source),
        leaf_count = length(source.leaf_sequences),
        fixed_dimension = size(source.sequence.coefficient_matrix, 2),
    )
end

function bond_aligned_homonuclear_chain_nested_geometry_diagnostics(
    basis::BondAlignedHomonuclearChainQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    odd_chain_policy::Symbol = :strict_current,
)
    source = _bond_aligned_homonuclear_chain_nested_fixed_source(
        basis;
        expansion = expansion,
        gausslet_backend = gausslet_backend,
        nside = nside,
        min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        odd_chain_policy = odd_chain_policy,
    )
    return _bond_aligned_homonuclear_chain_nested_geometry_diagnostics(source)
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

"""
    experimental_bond_aligned_homonuclear_chain_nested_qw_operators(
        basis::BondAlignedHomonuclearChainQWBasis3D;
        nuclear_charges = fill(1.0, length(basis.nuclei)),
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
    nuclear_charges::AbstractVector{<:Real} = fill(1.0, length(basis.nuclei)),
    nuclear_term_storage::Symbol = :auto,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    odd_chain_policy::Symbol = :central_ternary_relaxed,
    timing::Bool = false,
)
    source = _bond_aligned_homonuclear_chain_nested_fixed_source(
        basis;
        expansion = expansion,
        gausslet_backend = gausslet_backend,
        nside = nside,
        min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        odd_chain_policy = odd_chain_policy,
    )
    diagnostics = _bond_aligned_homonuclear_chain_nested_geometry_diagnostics(source)
    fixed_block = _nested_fixed_block(source)
    operators = ordinary_cartesian_qiu_white_operators(
        fixed_block;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        timing = timing,
    )
    return ExperimentalBondAlignedHomonuclearChainNestedQWPath(
        basis,
        source,
        fixed_block,
        operators,
        diagnostics,
        Float64[Float64(value) for value in nuclear_charges],
        odd_chain_policy,
    )
end

"""
    experimental_axis_aligned_homonuclear_square_lattice_nested_qw_operators(
        basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D;
        nuclear_charges = fill(1.0, length(basis.nuclei)),
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
    nuclear_charges::AbstractVector{<:Real} = fill(1.0, length(basis.nuclei)),
    nuclear_term_storage::Symbol = :auto,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
    timing::Bool = false,
)
    source = _axis_aligned_homonuclear_square_lattice_nested_fixed_source(
        basis;
        expansion = expansion,
        gausslet_backend = gausslet_backend,
        nside = nside,
        min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
    )
    diagnostics = _axis_aligned_homonuclear_square_lattice_nested_geometry_diagnostics(source)
    fixed_block = _nested_fixed_block(source)
    operators = ordinary_cartesian_qiu_white_operators(
        fixed_block;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        timing = timing,
    )
    return ExperimentalAxisAlignedHomonuclearSquareLatticeNestedQWPath(
        basis,
        source,
        fixed_block,
        operators,
        diagnostics,
        Float64[Float64(value) for value in nuclear_charges],
        min_in_plane_aspect_ratio,
    )
end
