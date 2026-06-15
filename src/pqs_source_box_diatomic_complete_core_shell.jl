struct _PQSDiatomicCompleteCoreShellHamReadinessPayload
    status::Symbol
    blocker
    route_family::Symbol
    system_classification
    bond_axis
    center_summary
    parent_axis_bundle_object_available::Bool
    route_skeleton_summary
    source_box_summary
    retained_unit_summary
    pair_inventory_summary
    available_objects::Tuple
    missing_objects::Tuple
    summary
    metadata
end

function _pqs_source_box_route_driver_diatomic_center_summary(parent)
    center_table =
        hasproperty(parent, :center_table) ? Tuple(parent.center_table) : ()
    return (;
        center_count = length(center_table),
        center_keys = Tuple(
            _pqs_source_box_route_driver_descriptor_property(
                center,
                :center_key,
                index,
            ) for (index, center) in pairs(center_table)
        ),
        nuclear_charges = Tuple(
            _pqs_source_box_route_driver_descriptor_property(
                center,
                :nuclear_charge,
            ) for center in center_table
        ),
        locations = Tuple(
            _pqs_source_box_route_driver_descriptor_property(center, :location)
            for center in center_table
        ),
    )
end

function _pqs_source_box_route_driver_diatomic_source_box_summary(route_skeleton)
    source_boxes =
        hasproperty(route_skeleton, :source_boxes) ?
        route_skeleton.source_boxes :
        (;)
    source_dimensions =
        hasproperty(route_skeleton, :source_dimensions) ?
        route_skeleton.source_dimensions :
        (;)
    source_box_keys = Tuple(keys(source_boxes))
    return (;
        source_box_count = length(source_box_keys),
        source_box_keys,
        source_dimensions,
        source_boxes_available = !isempty(source_box_keys),
    )
end

function _pqs_source_box_route_driver_diatomic_retained_unit_summary(route_skeleton)
    retained_units =
        hasproperty(route_skeleton, :retained_units) ?
        Tuple(route_skeleton.retained_units) :
        ()
    retained_ranges =
        Tuple(
            _pqs_source_box_route_driver_descriptor_property(
                unit,
                :retained_range,
            ) for unit in retained_units
        )
    return (;
        retained_unit_count = length(retained_units),
        unit_keys = Tuple(
            _pqs_source_box_route_driver_descriptor_property(unit, :unit_key)
            for unit in retained_units
        ),
        retained_unit_kinds = Tuple(
            _pqs_source_box_route_driver_descriptor_property(
                unit,
                :retained_unit_kind,
            ) for unit in retained_units
        ),
        retained_counts =
            hasproperty(route_skeleton, :retained_counts) ?
            route_skeleton.retained_counts :
            nothing,
        retained_ranges_available =
            !isempty(retained_ranges) && all(!isnothing, retained_ranges),
        retained_dimension =
            hasproperty(route_skeleton, :retained_dimension) ?
            route_skeleton.retained_dimension :
            nothing,
        weight_semantics = Tuple(
            _pqs_source_box_route_driver_descriptor_property(
                unit,
                :weight_semantics,
            ) for unit in retained_units
        ),
    )
end

function _pqs_source_box_route_driver_diatomic_pair_inventory_summary(route_skeleton)
    pair_entries =
        hasproperty(route_skeleton, :pair_entries) ?
        Tuple(route_skeleton.pair_entries) :
        ()
    source_box_algorithmic_flags = Tuple(
        _pqs_source_box_route_driver_descriptor_property(
            entry,
            :source_box_algorithmic_path,
            false,
        ) for entry in pair_entries
    )
    return (;
        pair_count = length(pair_entries),
        pair_family_counts =
            hasproperty(route_skeleton, :pair_family_counts) ?
            route_skeleton.pair_family_counts :
            nothing,
        pair_families = Tuple(
            unique(
                _pqs_source_box_route_driver_descriptor_property(
                    entry,
                    :pair_family,
                    :unknown,
                ) for entry in pair_entries
            ),
        ),
        helper_by_pair_family =
            hasproperty(route_skeleton, :helper_by_pair_family) ?
            route_skeleton.helper_by_pair_family :
            nothing,
        source_box_algorithmic_path_for_all_pairs =
            !isempty(source_box_algorithmic_flags) &&
            all(source_box_algorithmic_flags),
    )
end

function _pqs_source_box_route_driver_diatomic_route_skeleton_summary(
    route_skeleton,
    retained_unit_summary,
)
    return (;
        object_kind =
            _pqs_source_box_route_driver_descriptor_property(
                route_skeleton,
                :object_kind,
            ),
        status =
            _pqs_source_box_route_driver_descriptor_property(
                route_skeleton,
                :status,
            ),
        route_shape =
            _pqs_source_box_route_driver_descriptor_property(
                route_skeleton,
                :route_shape,
            ),
        retained_dimension = retained_unit_summary.retained_dimension,
    )
end

struct _PQSDiatomicCompleteCoreShellSupportWindowPayload
    status::Symbol
    blocker
    route_family::Symbol
    system_classification
    bond_axis
    parent_dims
    parent_axis_bundle_object_available::Bool
    source_box_windows
    source_mode_dims
    retained_order::Tuple
    candidate_core_then_shell_support_order::Tuple
    retained_to_support_order_permutation_required::Bool
    support_counts
    missing_objects::Tuple
    summary
    metadata
end

struct _PQSDiatomicRawBoxRoutePayload
    status::Symbol
    blocker
    producer
    producer_status::Symbol
    descriptor_summary
    raw_product_box_plan_summary
    raw_pqs_plan_summary
    product_unit_summary
    pair_inventory_summary
    support_window_payload_status::Symbol
    available_objects::Tuple
    missing_objects::Tuple
    summary
    metadata
end

struct _PQSDiatomicCompleteCoreShellSourceRealizationPayload
    status::Symbol
    blocker
    route_family::Symbol
    system_classification
    bond_axis
    parent_axis_bundle_object_available::Bool
    support_window_payload_status::Symbol
    raw_box_route_payload_status::Symbol
    core_unit_key::Symbol
    shell_unit_keys::Tuple
    retained_order::Tuple
    support_order::Tuple
    retained_to_support_order_permutation_required::Bool
    route_retained_ranges
    source_plan_precleanup_ranges
    core_support_count
    shell_support_counts
    shell_support_count
    shell_retained_counts
    shell_retained_count
    precleanup_retained_dimension
    shell_final_coefficients_shape
    shell_coefficient_block_structure::Symbol
    bundles_role::Symbol
    object_kind_claim::Symbol
    available_objects::Tuple
    missing_objects::Tuple
    summary
    metadata
end

struct _PQSDiatomicCompleteCoreShellSourcePlan
    object_kind::Symbol
    status::Symbol
    blocker
    bundles
    metrics
    core_unit_key::Symbol
    shell_unit_keys::Tuple
    core_support_indices::Vector{Int}
    core_support_states::Vector{NTuple{3,Int}}
    shell_support_indices::Vector{Int}
    shell_support_states::Vector{NTuple{3,Int}}
    shell_final_coefficients::Matrix{Float64}
    support_order::Tuple
    route_retained_order::Tuple
    retained_pre_final_map
    source_unit_summaries
    convention_labels
    summary
    metadata
end

struct _PQSDiatomicCompleteCoreShellFinalBasisPayload
    status::Symbol
    blocker
    route_family::Symbol
    source_plan
    source_plan_status::Symbol
    final_basis
    final_basis_status::Symbol
    available_objects::Tuple
    missing_objects::Tuple
    summary
    metadata
end

struct _PQSDiatomicCompleteCoreShellH1Payload
    status::Symbol
    blocker
    route_family::Symbol
    source_plan
    source_plan_status::Symbol
    final_basis
    final_basis_status::Symbol
    support_kinetic
    support_kinetic_status::Symbol
    support_electron_nuclear_by_center
    support_electron_nuclear_status::Symbol
    final_kinetic
    final_kinetic_status::Symbol
    final_electron_nuclear_by_center
    final_electron_nuclear_status::Symbol
    final_hamiltonian
    final_hamiltonian_status::Symbol
    h1
    h1_status::Symbol
    available_objects::Tuple
    missing_objects::Tuple
    summary
    metadata
end

struct _PQSDiatomicCompleteCoreShellHamInputPayload
    status::Symbol
    blocker
    route_family::Symbol
    source_plan
    source_plan_status::Symbol
    final_basis
    final_basis_status::Symbol
    h1_payload
    h1_payload_status::Symbol
    density_provenance
    density_provenance_status::Symbol
    support_weights
    support_weights_status::Symbol
    raw_pair_factor_terms
    raw_pair_factor_status::Symbol
    support_pair_raw_numerator
    support_pair_raw_numerator_status::Symbol
    density_interaction
    density_interaction_status::Symbol
    available_objects::Tuple
    missing_objects::Tuple
    summary
    metadata
end

struct _PQSDiatomicCompleteCoreShellHamiltonianHandoffPayload
    status::Symbol
    blocker
    route_family::Symbol
    source_plan
    source_plan_status::Symbol
    final_basis
    final_basis_status::Symbol
    h1_payload
    h1_payload_status::Symbol
    ham_input_payload
    ham_input_payload_status::Symbol
    one_body_hamiltonian
    one_body_hamiltonian_status::Symbol
    density_interaction
    density_interaction_status::Symbol
    pre_final_pair_matrix
    final_to_pre_final_coefficients
    pre_final_weights
    support_weights
    support_pair_raw_numerator
    raw_pair_factor_terms
    center_records
    center_metadata
    nuclear_repulsion_status::Symbol
    nuclear_repulsion
    electron_count_status::Symbol
    electron_count
    spin_sector_status::Symbol
    spin_sector
    ordering
    conventions
    available_objects::Tuple
    missing_objects::Tuple
    summary
    metadata
end

struct _PQSDiatomicCompleteCoreShellHamiltonianConsumerContractPayload
    status::Symbol
    blocker
    route_family::Symbol
    source_handoff
    source_handoff_status::Symbol
    readiness
    available_objects::Tuple
    missing_objects::Tuple
    summary
    metadata
end

struct _PQSDiatomicPhysicalGaussletCoreShellTargetPayload
    status::Symbol
    blocker
    route_family::Symbol
    route_kind::Symbol
    parent_axis_counts
    support_units::Tuple
    retained_units::Tuple
    support_counts
    retained_counts
    retained_order::Tuple
    expected_final_dimension
    retained_atom_core_interiors::Bool
    source_plan_role::Symbol
    supplement_policy::Symbol
    available_objects::Tuple
    missing_objects::Tuple
    summary
    metadata
end

struct _PQSIndependentH2PQSSupplementSupportPartitionPayload
    status::Symbol
    blocker
    route_family::Symbol
    route_kind::Symbol
    support_order::Tuple
    retained_order::Tuple
    support_counts
    retained_counts
    retained_ranges
    unit_partitions::Tuple
    total_tile_count::Int
    total_support_count::Int
    coverage_complete::Bool
    duplicate_parent_row_count::Int
    missing_parent_row_count::Int
    outside_parent_row_count::Int
    source_backed_fixed_source_oracle_used::Bool
    fake_pqs_enabled::Bool
    retained_transform_authority::Symbol
    summary
    metadata
end

function _pqs_source_box_route_driver_support_tile(
    unit_key::Symbol,
    source_region_role::Symbol,
    tile_index::Int,
    box::NTuple{3,UnitRange{Int}},
    parent_dims::NTuple{3,Int},
)
    tile_key = Symbol(unit_key, "_", source_region_role, "_tile_", tile_index)
    support_indices = _nested_box_support_indices(box[1], box[2], box[3], parent_dims)
    support_states = NTuple{3,Int}[
        _cartesian_unflat_index(index, parent_dims) for index in support_indices
    ]
    cpb = CartesianCPB.cpb(
        box...;
        role = tile_key,
        metadata = (;
            unit_key,
            source_region_role,
            source = :independent_h2_pqs_supplement_support_partition,
        ),
    )
    return (;
        tile_key,
        unit_key,
        source_region_role,
        support_tile_role = :rectangular_parent_row_support_tile,
        cpb,
        cpb_role = tile_key,
        intervals = box,
        support_indices,
        support_states,
        support_count = length(support_indices),
        parent_row_min = minimum(support_indices),
        parent_row_max = maximum(support_indices),
        provider_tile_ready = true,
    )
end

function _pqs_source_box_route_driver_nonempty_box_push!(
    boxes,
    x::UnitRange{Int},
    y::UnitRange{Int},
    z::UnitRange{Int},
)
    (isempty(x) || isempty(y) || isempty(z)) && return boxes
    push!(boxes, (x, y, z))
    return boxes
end

function _pqs_source_box_route_driver_rectangular_difference_boxes(
    outer_box::NTuple{3,UnitRange{Int}},
    inner_box::NTuple{3,UnitRange{Int}},
)
    boxes = NTuple{3,UnitRange{Int}}[]
    _pqs_source_box_route_driver_nonempty_box_push!(
        boxes,
        first(outer_box[1]):(first(inner_box[1]) - 1),
        outer_box[2],
        outer_box[3],
    )
    _pqs_source_box_route_driver_nonempty_box_push!(
        boxes,
        (last(inner_box[1]) + 1):last(outer_box[1]),
        outer_box[2],
        outer_box[3],
    )
    _pqs_source_box_route_driver_nonempty_box_push!(
        boxes,
        inner_box[1],
        first(outer_box[2]):(first(inner_box[2]) - 1),
        outer_box[3],
    )
    _pqs_source_box_route_driver_nonempty_box_push!(
        boxes,
        inner_box[1],
        (last(inner_box[2]) + 1):last(outer_box[2]),
        outer_box[3],
    )
    _pqs_source_box_route_driver_nonempty_box_push!(
        boxes,
        inner_box[1],
        inner_box[2],
        first(outer_box[3]):(first(inner_box[3]) - 1),
    )
    _pqs_source_box_route_driver_nonempty_box_push!(
        boxes,
        inner_box[1],
        inner_box[2],
        (last(inner_box[3]) + 1):last(outer_box[3]),
    )
    return Tuple(boxes)
end

function _pqs_source_box_route_driver_shared_shell_support_tiles(
    role::Symbol,
    outer_box::NTuple{3,UnitRange{Int}},
    inner_box::NTuple{3,UnitRange{Int}},
    parent_dims::NTuple{3,Int},
)
    outer_cpb = CartesianCPB.filled_cpb(
        outer_box...;
        role = Symbol(role, "_outer_support_box"),
    )
    inner_cpb = CartesianCPB.filled_cpb(
        inner_box...;
        role = Symbol(role, "_inner_exclusion_box"),
    )
    boxes = try
        strata = CartesianCPB.complete_shell_boundary_strata(outer_cpb, inner_cpb)
        Tuple(CartesianCPB.intervals(cpb) for cpb in strata.all_strata)
    catch
        _pqs_source_box_route_driver_rectangular_difference_boxes(outer_box, inner_box)
    end
    return Tuple(
        _pqs_source_box_route_driver_support_tile(
            role,
            :shared_molecular_shell,
            index,
            box,
            parent_dims,
        ) for (index, box) in pairs(boxes)
    )
end

function _pqs_source_box_route_driver_independent_h2_support_region_plan(
    parent,
    route_skeleton,
    recipe,
)
    route_kind =
        hasproperty(route_skeleton, :route_kind) ?
        route_skeleton.route_kind :
        recipe.route_kind
    route_kind === :bond_aligned_diatomic_independent_pqs_source_box_core_shell ||
        return nothing
    blocked(blocker) = (;
        status = :blocked_independent_pqs_support_region_plan,
        blocker,
        authority = :cartesian_shellification_route_geometry,
        support_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2),
        support_counts = (;),
        counts_generated = false,
        counts_source = :support_region_materializer_blocked,
        coverage_complete = false,
        duplicate_count = nothing,
        missing_count = nothing,
        outside_count = nothing,
    )
    parent_axis_counts =
        hasproperty(route_skeleton, :parent_axis_counts) ?
        route_skeleton.parent_axis_counts :
        hasproperty(parent, :axis_counts) ? parent.axis_counts : nothing
    axis_count_tuple =
        _pqs_source_box_route_driver_axis_counts_tuple(parent_axis_counts)
    axis_count_tuple == (9, 9, 15) ||
        return blocked(:unexpected_independent_h2_parent_axis_counts)
    bundles =
        hasproperty(parent, :parent_axis_bundle_object) ?
        parent.parent_axis_bundle_object :
        nothing
    isnothing(bundles) && return blocked(:missing_parent_axis_bundle_object)
    axes = Tuple(collect(_nested_axis_pgdg(bundles, axis).centers)
                 for axis in (:x, :y, :z))
    locations =
        hasproperty(parent, :atom_locations) ? Tuple(parent.atom_locations) : ()
    length(locations) == 2 || return blocked(:missing_diatomic_atom_locations)
    bond_axis = hasproperty(parent, :bond_axis) ? parent.bond_axis : :z
    raw = CartesianShellification.raw_terminal_geometry(
        axes,
        Tuple(Tuple(Float64(value) for value in location) for location in locations);
        core_side = 5,
        q = 5,
        bond_axis,
    )
    regions = raw.regions
    atom_regions = Tuple(region for region in regions if region.role == :atom_local_core)
    midpoint_regions = Tuple(region for region in regions if region.role == :midpoint_slab)
    shared_regions = sort(
        collect(region for region in regions if region.role == :shared_molecular_shell);
        by = region -> prod(length.(region.outer_box)),
        rev = true,
    )
    length(atom_regions) == 2 && length(midpoint_regions) == 1 &&
        length(shared_regions) == 2 ||
        return blocked(:unexpected_independent_h2_terminal_region_shape)
    support_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
    support_counts = (;
        atom_contact_core =
            sum(region.support_count for region in (atom_regions..., midpoint_regions...)),
        shared_shell_1 = shared_regions[1].support_count,
        shared_shell_2 = shared_regions[2].support_count,
    )
    atom_contact_regions = (atom_regions..., midpoint_regions...)
    atom_contact_support(region) =
        _nested_box_support_indices(
            region.box[1],
            region.box[2],
            region.box[3],
            axis_count_tuple,
        )
    atom_contact_core_support_indices =
        reduce(vcat, (atom_contact_support(region) for region in atom_contact_regions))
    atom_contact_core_support_states = NTuple{3,Int}[
        _cartesian_unflat_index(index, axis_count_tuple) for
        index in atom_contact_core_support_indices
    ]
    atom_contact_support_tiles = Tuple(
        _pqs_source_box_route_driver_support_tile(
            :atom_contact_core,
            region.role,
            index,
            region.box,
            axis_count_tuple,
        ) for (index, region) in pairs(atom_contact_regions)
    )
    shared_shell_descriptors =
        _pqs_source_box_route_driver_independent_h2_shared_shell_region_descriptors(
            shared_regions,
            axis_count_tuple,
        )
    coverage = raw.coverage
    counts_match =
        Tuple(values(support_counts)) == (275, 578, 362)
    status =
        coverage.coverage_complete && counts_match ?
        :available_independent_pqs_support_region_plan :
        :blocked_independent_pqs_support_region_plan
    blocker =
        status === :available_independent_pqs_support_region_plan ?
        nothing :
        !counts_match ? :independent_h2_support_count_mismatch :
        :independent_h2_support_coverage_mismatch
    return (;
        status,
        blocker,
        authority = :cartesian_shellification_route_geometry,
        support_order,
        support_counts,
        counts_generated = status === :available_independent_pqs_support_region_plan,
        counts_source = :cartesian_shellification_route_geometry,
        coverage_complete = coverage.coverage_complete,
        duplicate_count = coverage.duplicate_site_count,
        missing_count = coverage.missing_site_count,
        outside_count = 0,
        nuclear_indices = raw.nuclear_indices,
        bond_axis = raw.bond_axis,
        atom_contact_core_descriptor = (;
            role = :atom_contact_core,
            source_region_roles = Tuple(region.role for region in atom_contact_regions),
            support_tiles = atom_contact_support_tiles,
            tile_count = length(atom_contact_support_tiles),
            tile_support_counts =
                Tuple(tile.support_count for tile in atom_contact_support_tiles),
            support_indices = atom_contact_core_support_indices,
            support_states = atom_contact_core_support_states,
            support_count = length(atom_contact_core_support_indices),
            coefficient_representation = :sparse_parent_row_direct_selector,
        ),
        shared_shell_descriptors,
    )
end

function _pqs_source_box_route_driver_independent_h2_shared_shell_region_descriptors(
    shared_regions,
    parent_dims::NTuple{3,Int},
)
    length(shared_regions) == 2 ||
        throw(ArgumentError("independent H2 PQS requires exactly two shared shell regions"))
    descriptor(role, region) = begin
        raw = hasproperty(region, :raw_region) ? region.raw_region : region
        region_key = hasproperty(region, :key) ? region.key : role
        isnothing(raw.inner_exclusion_box) &&
            throw(ArgumentError("independent H2 PQS shared shell requires inner exclusion box"))
        outer_indices =
            _nested_box_support_indices(raw.outer_box[1], raw.outer_box[2], raw.outer_box[3], parent_dims)
        inner = Set(
            _nested_box_support_indices(
                raw.inner_exclusion_box[1],
                raw.inner_exclusion_box[2],
                raw.inner_exclusion_box[3],
                parent_dims,
            ),
        )
        support_indices = Int[index for index in outer_indices if !(index in inner)]
        support_states = NTuple{3,Int}[
            _cartesian_unflat_index(index, parent_dims) for index in support_indices
        ]
        support_tiles =
            _pqs_source_box_route_driver_shared_shell_support_tiles(
                role,
                raw.outer_box,
                raw.inner_exclusion_box,
                parent_dims,
            )
        source_cpb = CartesianCPB.filled_cpb(
            raw.outer_box...;
            role = Symbol(role, "_pqs_filled_source_cpb"),
            metadata = (;
                terminal_region_key = region_key,
                route_unit_role = role,
                source = :independent_h2_pqs_shared_shell_region_descriptor,
            ),
        )
        (;
            role,
            terminal_region_key = region_key,
            current_box = raw.outer_box,
            inner_box = raw.inner_exclusion_box,
            source_cpb,
            support_tiles,
            tile_count = length(support_tiles),
            tile_support_counts = Tuple(tile.support_count for tile in support_tiles),
            support_indices,
            support_states,
            support_count = length(support_indices),
        )
    end
    return (;
        shared_shell_1 = descriptor(:shared_shell_1, shared_regions[1]),
        shared_shell_2 = descriptor(:shared_shell_2, shared_regions[2]),
    )
end

function _pqs_source_box_route_driver_independent_h2_retained_rule_plan(
    support_plan,
)
    support_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
    blocked(blocker) = (;
        object_kind = :pqs_independent_h2_core_shell_retained_rule_plan,
        status = :blocked_independent_pqs_retained_rule_readiness_plan,
        blocker,
        authority = :independent_pqs_route_owned_retained_rule_readiness,
        support_order,
        support_counts = (;),
        retained_order = (),
        retained_counts = (;),
        expected_final_dimension = nothing,
        retained_counts_generated = false,
        retained_counts_source = :retained_rule_plan_blocked,
        per_unit_authority = (;),
        missing_objects = (blocker,),
    )
    isnothing(support_plan) &&
        return blocked(:missing_independent_pqs_support_region_plan)
    support_plan.status === :available_independent_pqs_support_region_plan ||
        return blocked(support_plan.blocker)
    support_plan.support_order == support_order ||
        return blocked(:independent_pqs_support_order_mismatch)
    support_counts = support_plan.support_counts
    Tuple(values(support_counts)) == (275, 578, 362) ||
        return blocked(:independent_pqs_support_count_mismatch)

    shared_shell_rule =
        CartesianRawProductSources.pqs_boundary_product_mode_retained_rule(
            (5, 5, 5);
            source_key = :independent_h2_shared_shell_q5_source,
            metadata = (;
                route_unit_role = :shared_molecular_shell,
                source = :independent_h2_retained_rule_readiness,
            ),
        )
    retained_counts = (;
        atom_contact_core = support_counts.atom_contact_core,
        shared_shell_1 = shared_shell_rule.retained_count,
        shared_shell_2 = shared_shell_rule.retained_count,
    )
    counts_match = Tuple(values(retained_counts)) == (275, 98, 98)
    status =
        counts_match ?
        :available_independent_pqs_retained_rule_readiness_plan :
        :blocked_independent_pqs_retained_rule_readiness_plan
    blocker =
        counts_match ? nothing : :independent_pqs_retained_rule_count_mismatch
    return (;
        object_kind = :pqs_independent_h2_core_shell_retained_rule_plan,
        status,
        blocker,
        authority = :independent_pqs_route_owned_retained_rule_readiness,
        support_order,
        support_counts,
        retained_order = support_order,
        retained_counts,
        expected_final_dimension = sum(values(retained_counts)),
        retained_counts_generated = counts_match,
        retained_counts_source = :route_owned_retained_rule_plan,
        per_unit_authority = (;
            atom_contact_core = :direct_source_modes,
            shared_shell_1 = :pqs_boundary_comx_product_modes,
            shared_shell_2 = :pqs_boundary_comx_product_modes,
        ),
        per_unit_transform_kind = (;
            atom_contact_core = :identity_source_modes,
            shared_shell_1 = shared_shell_rule.transform_kind,
            shared_shell_2 = shared_shell_rule.transform_kind,
        ),
        shared_shell_source_mode_dims = shared_shell_rule.source_mode_dims,
        shared_shell_retained_rule_kind = shared_shell_rule.retained_rule_kind,
        missing_objects =
            status === :available_independent_pqs_retained_rule_readiness_plan ?
            (:independent_pqs_physical_source_plan_materializer,) :
            (blocker,),
    )
end

struct _PQSDiatomicPhysicalGaussletSupplementRequestPayload
    status::Symbol
    blocker
    route_family::Symbol
    route_kind::Symbol
    fixture_label::Symbol
    supplement_policy::Symbol
    atom_symbols::Tuple
    nuclear_charges::Tuple
    atom_locations::Tuple
    bond_axis
    bond_length
    basis_name
    lmax
    uncontracted
    residual_keep_policy::Symbol
    residual_drop_tolerance
    required_provider_blocks::Tuple
    available_fact_labels::Tuple
    missing_fact_labels::Tuple
    summary
    metadata
end

struct _PQSDiatomicPhysicalGaussletSupplementRepresentationPayload
    status::Symbol
    blocker
    route_family::Symbol
    route_kind::Symbol
    fixture_label::Symbol
    supplement_policy::Symbol
    object_kind::Symbol
    supplement
    representation
    basis_name
    lmax
    uncontracted
    atom_symbols::Tuple
    center_count
    orbital_count
    summary
    metadata
end

struct _PQSDiatomicPhysicalGaussletSupplementPreflightPayload
    status::Symbol
    blocker
    route_family::Symbol
    route_kind::Symbol
    fixture_label::Symbol
    support_counts
    retained_counts
    retained_order::Tuple
    retained_transform_kind::Symbol
    gausslet_final_dimension
    supplement_policy::Symbol
    required_fact_labels::Tuple
    available_fact_labels::Tuple
    missing_fact_labels::Tuple
    summary
    metadata
end

struct _PQSH2WLGaussletOnlyReferenceCandidatePayload
    status::Symbol
    blocker
    route_family::Symbol
    route_kind::Symbol
    label::String
    retained_transform_kind::Symbol
    supplement_policy::Symbol
    final_dimension
    conditions
    summary
    metadata
end

struct _PQSDiatomicPhysicalGaussletCoreShellSourcePlanPayload
    status::Symbol
    blocker
    route_family::Symbol
    source_plan
    available_objects::Tuple
    missing_objects::Tuple
    summary
    metadata
end

struct _PQSDiatomicPhysicalGaussletFinalBasisPayload
    status::Symbol
    blocker
    route_family::Symbol
    source_plan
    source_plan_status::Symbol
    final_basis
    final_basis_status::Symbol
    available_objects::Tuple
    missing_objects::Tuple
    summary
    metadata
end

struct _PQSDiatomicPhysicalGaussletH1Payload
    status::Symbol
    blocker
    route_family::Symbol
    source_plan
    source_plan_status::Symbol
    final_basis
    final_basis_status::Symbol
    support_kinetic
    support_kinetic_status::Symbol
    support_electron_nuclear_by_center
    support_electron_nuclear_status::Symbol
    final_kinetic
    final_kinetic_status::Symbol
    final_electron_nuclear_by_center
    final_electron_nuclear_status::Symbol
    final_hamiltonian
    final_hamiltonian_status::Symbol
    h1
    h1_status::Symbol
    available_objects::Tuple
    missing_objects::Tuple
    summary
    metadata
end

struct _PQSDiatomicPhysicalGaussletH1JPayload
    status::Symbol
    blocker
    route_family::Symbol
    source_plan
    source_plan_status::Symbol
    final_basis
    final_basis_status::Symbol
    h1_payload
    h1_payload_status::Symbol
    density_provenance
    density_provenance_status::Symbol
    support_weights
    support_weights_status::Symbol
    raw_pair_factor_terms
    raw_pair_factor_status::Symbol
    support_pair_raw_numerator
    support_pair_raw_numerator_status::Symbol
    density_interaction
    density_interaction_status::Symbol
    h1_j_diagnostic
    h1_j_status::Symbol
    available_objects::Tuple
    missing_objects::Tuple
    summary
    metadata
end

struct _PQSDiatomicPhysicalGaussletRHFInputContractPayload
    status::Symbol
    blocker
    route_family::Symbol
    source_plan
    source_plan_status::Symbol
    final_basis
    final_basis_status::Symbol
    h1_payload
    h1_payload_status::Symbol
    h1_j_payload
    h1_j_payload_status::Symbol
    electron_count
    occupation
    available_objects::Tuple
    missing_objects::Tuple
    summary
    metadata
end

struct _PQSDiatomicPhysicalGaussletRHFExecutionPayload
    status::Symbol
    blocker
    route_family::Symbol
    input_contract
    input_contract_status::Symbol
    h1_payload
    h1_payload_status::Symbol
    h1_j_payload
    h1_j_payload_status::Symbol
    summary
    metadata
end

struct _PQSDiatomicPhysicalGaussletSourcePlanCandidatePayload
    status::Symbol
    blocker
    candidate_source::Symbol
    candidate
    counts_match::Bool
    authority_status::Symbol
    available_objects::Tuple
    missing_objects::Tuple
    summary
    metadata
end

struct _PQSDiatomicPhysicalGaussletCoreShellSourcePlan
    object_kind::Symbol
    status::Symbol
    parent_basis
    axis_bundles
    atom_contact_core_support_indices::Vector{Int}
    atom_contact_core_support_states
    shared_shell_support_indices::Tuple
    shared_shell_support_states::Tuple
    core_coefficient_matrix
    shared_shell_coefficient_matrices::Tuple
    support_order::Tuple
    retained_order::Tuple
    retained_ranges
    final_dimension::Int
    convention_labels
    summary
    metadata
end

function _pqs_source_box_route_driver_axis_counts_tuple(axis_counts)
    isnothing(axis_counts) && return nothing
    if axis_counts isa NamedTuple
        return (Int(axis_counts.x), Int(axis_counts.y), Int(axis_counts.z))
    end
    axis_counts_tuple = Tuple(axis_counts)
    length(axis_counts_tuple) == 3 || return nothing
    return ntuple(axis -> Int(axis_counts_tuple[axis]), 3)
end

function _pqs_source_box_route_driver_ordered_count_tuple(counts, order)
    isnothing(counts) && return nothing
    isempty(order) && return ()
    counts isa NamedTuple || return nothing
    all(hasproperty(counts, key) for key in order) || return nothing
    return Tuple(Int(getproperty(counts, key)) for key in order)
end

function _pqs_source_box_route_driver_location_tuple(location)
    isnothing(location) && return nothing
    location_tuple = Tuple(location)
    length(location_tuple) == 3 || return nothing
    return ntuple(axis -> Float64(location_tuple[axis]), 3)
end

function _pqs_source_box_route_driver_h2_r4_center_match(parent)
    summary = _pqs_source_box_route_driver_diatomic_center_summary(parent)
    summary.center_count == 2 || return false
    all(charge -> isapprox(Float64(charge), 1.0; atol = 1.0e-12),
        summary.nuclear_charges) || return false
    locations =
        _pqs_source_box_route_driver_location_tuple.(summary.locations)
    all(!isnothing, locations) || return false
    ordered = sort(collect(locations); by = location -> location[3])
    expected = ((0.0, 0.0, -2.0), (0.0, 0.0, 2.0))
    return all(
        isapprox(ordered[index][axis], expected[index][axis]; atol = 1.0e-12)
        for index in 1:2 for axis in 1:3
    )
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_target_payload(
    parent,
    route_skeleton,
    recipe,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    route_kind =
        hasproperty(route_skeleton, :route_kind) ?
        route_skeleton.route_kind :
        recipe.route_kind
    system_classification =
        hasproperty(parent, :system_classification) ?
        parent.system_classification :
        nothing
    parent_axis_counts =
        hasproperty(route_skeleton, :parent_axis_counts) ?
        route_skeleton.parent_axis_counts :
        hasproperty(parent, :axis_counts) ?
        parent.axis_counts :
        nothing
    parent_axis_count_tuple =
        _pqs_source_box_route_driver_axis_counts_tuple(parent_axis_counts)
    inventory =
        hasproperty(route_skeleton, :physical_target_inventory) ?
        route_skeleton.physical_target_inventory :
        nothing
    independent_target = route_kind ===
                         :bond_aligned_diatomic_independent_pqs_source_box_core_shell
    expected_route_kinds = (
        :bond_aligned_diatomic_physical_gausslet_core_shell_pqs,
        :bond_aligned_diatomic_independent_pqs_source_box_core_shell,
    )
    expected_parent_axis_counts = (9, 9, 15)

    if route_family !== :pqs_source_box
        status = :not_applicable_physical_gausslet_target_non_pqs_route
        blocker = nothing
    elseif !(route_kind in expected_route_kinds)
        status = :not_applicable_physical_gausslet_target_route_kind
        blocker = nothing
    elseif system_classification !== :bond_aligned_diatomic
        status = :blocked_physical_gausslet_target_inventory
        blocker = :not_bond_aligned_diatomic
    elseif isnothing(inventory)
        status = :blocked_physical_gausslet_target_inventory
        blocker = :missing_physical_gausslet_target_inventory
    elseif parent_axis_count_tuple != expected_parent_axis_counts
        status = :blocked_physical_gausslet_target_inventory
        blocker = :unexpected_physical_gausslet_parent_axis_counts
    elseif independent_target
        status = inventory.status
        blocker = inventory.blocker
    else
        status = :available_physical_gausslet_core_shell_target_inventory
        blocker = nothing
    end

    support_units =
        isnothing(inventory) ? () : Tuple(inventory.support_units)
    retained_units =
        isnothing(inventory) ? () : Tuple(inventory.retained_units)
    generated_support_plan =
        independent_target ?
        _pqs_source_box_route_driver_independent_h2_support_region_plan(
            parent,
            route_skeleton,
            recipe,
        ) :
        nothing
    generated_support_available =
        !isnothing(generated_support_plan) &&
        generated_support_plan.status === :available_independent_pqs_support_region_plan
    retained_rule_plan =
        independent_target ?
        _pqs_source_box_route_driver_independent_h2_retained_rule_plan(
            generated_support_plan,
        ) :
        nothing
    retained_rule_plan_available =
        !isnothing(retained_rule_plan) &&
        retained_rule_plan.status ===
        :available_independent_pqs_retained_rule_readiness_plan
    if independent_target && retained_rule_plan_available
        status = :available_physical_gausslet_core_shell_target_inventory
        blocker = nothing
    end
    support_counts =
        generated_support_available ? generated_support_plan.support_counts :
        isnothing(inventory) ? (;) : inventory.support_counts
    retained_units =
        retained_rule_plan_available ? retained_rule_plan.retained_order :
        retained_units
    retained_counts =
        retained_rule_plan_available ? retained_rule_plan.retained_counts :
        isnothing(inventory) ? (;) : inventory.retained_counts
    retained_order =
        retained_rule_plan_available ? retained_rule_plan.retained_order :
        isnothing(inventory) ? () : Tuple(inventory.retained_order)
    expected_final_dimension =
        retained_rule_plan_available ? retained_rule_plan.expected_final_dimension :
        isnothing(inventory) ? nothing : inventory.expected_final_dimension
    retained_atom_core_interiors =
        retained_rule_plan_available ||
        (!isnothing(inventory) && inventory.retained_atom_core_interiors)
    source_plan_role =
        isnothing(inventory) ?
        :not_available :
        inventory.source_plan_role
    source_backed_fixed_source_oracle_used =
        !isnothing(inventory) && get(inventory, :source_backed_fixed_source_oracle_used, false)
    source_plan_blocker =
        isnothing(inventory) ? nothing : get(inventory, :source_plan_blocker, nothing)
    retained_transform_authority =
        retained_rule_plan_available ?
        retained_rule_plan.authority :
        isnothing(inventory) ?
        :not_available :
        get(inventory, :retained_transform_authority, :not_available)
    primary_blocker =
        retained_rule_plan_available ?
        source_plan_blocker :
        isnothing(inventory) ? nothing : get(inventory, :primary_blocker, blocker)
    secondary_blocker =
        retained_rule_plan_available ?
        nothing :
        isnothing(inventory) ? nothing : get(inventory, :secondary_blocker, nothing)
    support_plan =
        !isnothing(generated_support_plan) ? generated_support_plan :
        isnothing(inventory) ? nothing : get(inventory, :support_plan, nothing)
    supplement_policy =
        isnothing(inventory) ?
        :not_available :
        something(get(recipe, :supplement_policy, nothing), inventory.supplement_policy)
    available_objects =
        status === :available_physical_gausslet_core_shell_target_inventory ?
        (:physical_gausslet_core_shell_target_inventory,) :
        ()
    missing_objects =
        independent_target && retained_rule_plan_available ?
        retained_rule_plan.missing_objects :
        status === :available_physical_gausslet_core_shell_target_inventory ?
        (
            :physical_gausslet_source_plan_producer,
            :physical_gausslet_final_basis_builder,
            :physical_gausslet_h1_builder,
        ) :
        isnothing(blocker) ? () : (blocker,)
    summary = (;
        status,
        blocker,
        route_family,
        route_kind,
        parent_axis_counts,
        support_units,
        support_counts,
        retained_units,
        retained_counts,
        retained_order,
        expected_final_dimension,
        retained_atom_core_interiors,
        source_plan_role,
        supplement_policy,
        target_inventory_available =
            status === :available_physical_gausslet_core_shell_target_inventory,
        source_plan_materialized = false,
        source_backed_fixed_source_oracle_used,
        retained_transform_authority,
        primary_blocker,
        secondary_blocker,
        independent_source_plan_blocker = source_plan_blocker,
        support_plan,
        retained_rule_plan,
        retained_rule_plan_status =
            isnothing(retained_rule_plan) ? :not_available : retained_rule_plan.status,
        retained_rule_plan_blocker =
            isnothing(retained_rule_plan) ? nothing : retained_rule_plan.blocker,
        retained_counts_generated =
            !isnothing(retained_rule_plan) && retained_rule_plan.retained_counts_generated,
        retained_counts_source =
            isnothing(retained_rule_plan) ?
            :not_available :
            retained_rule_plan.retained_counts_source,
        final_basis_materialized = false,
        h1_materialized = false,
        h1_j_materialized = false,
        rhf_materialized = false,
        missing_objects,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_physical_gausslet_target_payload,
        route_kind,
        provenance =
            isnothing(inventory) ? :missing_inventory : inventory.provenance,
        target_inventory_hard_coded_from_reviewed_contract = !independent_target,
        reviewed_contract_pass = independent_target ? nothing : 200,
        old_wl_qw_fixed_block_size = independent_target ? nothing : (1215, 463),
        source_plan_role,
    )

    return _PQSDiatomicPhysicalGaussletCoreShellTargetPayload(
        status,
        blocker,
        route_family,
        route_kind,
        parent_axis_counts,
        support_units,
        retained_units,
        support_counts,
        retained_counts,
        retained_order,
        expected_final_dimension,
        retained_atom_core_interiors,
        source_plan_role,
        supplement_policy,
        available_objects,
        missing_objects,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_diatomic_bond_length(atom_locations)
    length(atom_locations) == 2 || return nothing
    a = atom_locations[1]
    b = atom_locations[2]
    length(a) == 3 && length(b) == 3 || return nothing
    return sqrt(sum((Float64(a[axis]) - Float64(b[axis]))^2 for axis in 1:3))
end

function _pqs_source_box_route_driver_supplement_constructor_basis_name(basis_name)
    isnothing(basis_name) && return basis_name
    parts = split(String(basis_name), '/')
    return String(parts[end])
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_request_payload(
    parent,
    target_payload,
)
    target_available =
        !isnothing(target_payload) &&
        target_payload.status === :available_physical_gausslet_core_shell_target_inventory
    center_table =
        hasproperty(parent, :center_table) ? Tuple(parent.center_table) : ()
    atom_symbols =
        Tuple(
            _pqs_source_box_route_driver_descriptor_property(
                center,
                :atom_symbol,
                :unknown,
            ) for center in center_table
        )
    nuclear_charges =
        Tuple(
            _pqs_source_box_route_driver_descriptor_property(
                center,
                :nuclear_charge,
                nothing,
            ) for center in center_table
        )
    atom_locations =
        Tuple(
            _pqs_source_box_route_driver_location_tuple(
                _pqs_source_box_route_driver_descriptor_property(
                    center,
                    :location,
                    nothing,
                ),
            ) for center in center_table
        )
    bond_axis = hasproperty(parent, :bond_axis) ? parent.bond_axis : nothing
    bond_length = _pqs_source_box_route_driver_diatomic_bond_length(atom_locations)
    fixture_label = :h2_r4_physical_gausslet_q5
    supplement_policy =
        target_available ? target_payload.supplement_policy : :not_available
    route_family = target_available ? target_payload.route_family : :not_available
    route_kind = target_available ? target_payload.route_kind : :not_available
    basis_name =
        supplement_policy === :mwg_residual_gto ? "H/cc-pVTZ" : nothing
    lmax = supplement_policy === :mwg_residual_gto ? 1 : nothing
    uncontracted = supplement_policy === :mwg_residual_gto ? false : nothing
    residual_keep_policy =
        supplement_policy === :mwg_residual_gto ?
        :route_private_mwg_residual_gto_preflight_only :
        :not_applicable
    residual_drop_tolerance =
        supplement_policy === :mwg_residual_gto ? 1.0e-10 : nothing
    required_provider_blocks =
        supplement_policy === :mwg_residual_gto ?
        (
            :mixed_gausslet_gto_blocks,
            :gto_gto_blocks,
            :combined_raw_moment_matrices,
            :residual_mwg_representation,
            :combined_density_density_readiness,
        ) :
        ()

    if !target_available
        status = :blocked_pqs_physical_gausslet_supplement_request
        blocker = :missing_physical_gausslet_target_inventory
        available_fact_labels = ()
        missing_fact_labels = (blocker,)
    elseif supplement_policy === :none
        status = :not_requested
        blocker = nothing
        available_fact_labels = ()
        missing_fact_labels = ()
    elseif supplement_policy === :mwg_residual_gto
        status = :available_pqs_physical_gausslet_supplement_request
        blocker = nothing
        available_fact_labels = (
            :physical_gausslet_core_shell_target_inventory,
            :route_owned_gto_supplement_request,
            :h2_r4_geometry,
            :h_cc_pvtz_lmax1_supplement_identity,
            :gto_supplement_representation,
        )
        missing_fact_labels = (
            :missing_provider_gto_supplement_blocks,
        )
    else
        status = :blocked_pqs_physical_gausslet_supplement_request
        blocker = :unsupported_physical_gausslet_supplement_policy
        available_fact_labels = ()
        missing_fact_labels = (blocker,)
    end

    summary = (;
        status,
        blocker,
        route_family,
        route_kind,
        fixture_label,
        supplement_policy,
        atom_symbols,
        nuclear_charges,
        atom_locations,
        bond_axis,
        bond_length,
        basis_name,
        lmax,
        uncontracted,
        residual_keep_policy,
        residual_drop_tolerance,
        required_provider_blocks,
        available_fact_labels,
        missing_fact_labels,
        matrices_materialized = false,
        supplemented_values_materialized = false,
        public_api = false,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_request_payload,
        boundary = :request_metadata_only,
        route_private = true,
        provider_blocks_materialized = false,
        gto_mwg_materialization = false,
    )
    return _PQSDiatomicPhysicalGaussletSupplementRequestPayload(
        status,
        blocker,
        route_family,
        route_kind,
        fixture_label,
        supplement_policy,
        atom_symbols,
        nuclear_charges,
        atom_locations,
        bond_axis,
        bond_length,
        basis_name,
        lmax,
        uncontracted,
        residual_keep_policy,
        residual_drop_tolerance,
        required_provider_blocks,
        available_fact_labels,
        missing_fact_labels,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_representation_payload(
    request_payload,
)
    if isnothing(request_payload)
        status = :not_available
        blocker = :missing_physical_gausslet_supplement_request_payload
        route_family = :not_available
        route_kind = :not_available
        fixture_label = :h2_r4_physical_gausslet_q5
        supplement_policy = :not_available
        object_kind = :not_available
        supplement = nothing
        representation = nothing
        basis_name = nothing
        lmax = nothing
        uncontracted = nothing
        atom_symbols = ()
    else
        route_family = request_payload.route_family
        route_kind = request_payload.route_kind
        fixture_label = request_payload.fixture_label
        supplement_policy = request_payload.supplement_policy
        basis_name = request_payload.basis_name
        lmax = request_payload.lmax
        uncontracted = request_payload.uncontracted
        atom_symbols = request_payload.atom_symbols
        if supplement_policy === :none
            status = :not_requested
            blocker = nothing
            object_kind = :not_requested
            supplement = nothing
            representation = nothing
        elseif request_payload.status !==
               :available_pqs_physical_gausslet_supplement_request
            status = :blocked_pqs_physical_gausslet_gto_supplement_representation
            blocker = request_payload.blocker
            object_kind = :not_available
            supplement = nothing
            representation = nothing
        else
            locations = request_payload.atom_locations
            locations_valid =
                length(locations) == 2 &&
                all(
                    location -> !isnothing(location) && length(location) == 3,
                    locations,
                )
            if !locations_valid
                status = :blocked_pqs_physical_gausslet_gto_supplement_representation
                blocker = :missing_h2_supplement_nuclei
                object_kind = :not_available
                supplement = nothing
                representation = nothing
            else
                nuclei = NTuple{3,Float64}[
                    ntuple(axis -> Float64(location[axis]), 3)
                    for location in locations
                ]
                supplement = legacy_bond_aligned_diatomic_gaussian_supplement(
                    "H",
                    _pqs_source_box_route_driver_supplement_constructor_basis_name(
                        basis_name,
                    ),
                    nuclei;
                    lmax,
                    uncontracted,
                )
                representation = basis_representation(supplement)
                status =
                    :available_pqs_physical_gausslet_gto_supplement_representation
                blocker = nothing
                object_kind =
                    :cartesian_gaussian_shell_supplement_representation
            end
        end
    end

    center_count = isnothing(representation) ? 0 : length(representation.metadata.nuclei)
    orbital_count = isnothing(representation) ? 0 : length(representation.orbitals)
    summary = (;
        status,
        blocker,
        route_family,
        route_kind,
        fixture_label,
        supplement_policy,
        object_kind,
        basis_name,
        lmax,
        uncontracted,
        atom_symbols,
        center_count,
        orbital_count,
        matrices_materialized = false,
        provider_blocks_materialized = false,
        public_api = false,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_representation_payload,
        boundary = :route_private_matrix_free_representation,
        representation_materialized =
            status ===
            :available_pqs_physical_gausslet_gto_supplement_representation,
        provider_blocks_materialized = false,
        matrices_materialized = false,
    )
    return _PQSDiatomicPhysicalGaussletSupplementRepresentationPayload(
        status,
        blocker,
        route_family,
        route_kind,
        fixture_label,
        supplement_policy,
        object_kind,
        supplement,
        representation,
        basis_name,
        lmax,
        uncontracted,
        atom_symbols,
        center_count,
        orbital_count,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_preflight_payload(
    target_payload,
    request_payload = nothing,
    representation_payload = nothing,
)
    target_available =
        !isnothing(target_payload) &&
        target_payload.status === :available_physical_gausslet_core_shell_target_inventory
    fixture_label = :h2_r4_physical_gausslet_q5
    retained_transform_kind = :pqs
    required_fact_labels =
        target_available && target_payload.supplement_policy === :mwg_residual_gto ?
        (
            :gto_supplement_representation,
            :provider_gto_supplement_blocks,
            :mixed_gausslet_gto_blocks,
            :gto_gto_blocks,
            :combined_raw_moment_matrices,
            :residual_mwg_representation,
            :combined_density_density_readiness,
        ) :
        ()
    available_fact_labels =
        target_available ?
        (
            :physical_gausslet_core_shell_target_inventory,
            :pqs_retained_transform_kind,
            :gausslet_only_final_dimension,
            (
                !isnothing(representation_payload) &&
                representation_payload.status ===
                :available_pqs_physical_gausslet_gto_supplement_representation
            ) ? :gto_supplement_representation : nothing,
        ) :
        ()
    available_fact_labels =
        Tuple(label for label in available_fact_labels if !isnothing(label))
    if !target_available
        status = :blocked_pqs_physical_gausslet_supplement_preflight
        blocker = :missing_physical_gausslet_target_inventory
        missing_fact_labels = (blocker,)
    elseif target_payload.supplement_policy === :none
        status = :not_requested
        blocker = nothing
        missing_fact_labels = ()
    elseif target_payload.supplement_policy === :mwg_residual_gto
        status = :blocked_pqs_physical_gausslet_mwg_residual_gto_preflight
        request_blocker =
            !isnothing(request_payload) &&
            request_payload.status !== :available_pqs_physical_gausslet_supplement_request ?
            request_payload.blocker :
            nothing
        representation_blocker =
            isnothing(request_blocker) &&
            (
                isnothing(representation_payload) ||
                representation_payload.status !==
                :available_pqs_physical_gausslet_gto_supplement_representation
            ) ?
            :missing_gto_supplement_representation :
            nothing
        blocker =
            !isnothing(request_blocker) ?
            request_blocker :
            !isnothing(representation_blocker) ?
            representation_blocker :
            :missing_provider_gto_supplement_blocks
        missing_fact_labels =
            !isnothing(representation_blocker) ?
            (
                :missing_gto_supplement_representation,
                :missing_provider_gto_supplement_blocks,
                :missing_mixed_gausslet_gto_blocks,
                :missing_gto_gto_blocks,
                :missing_combined_raw_moment_matrices,
                :missing_residual_mwg_representation,
                :missing_combined_density_density_readiness,
            ) :
            (
                :missing_provider_gto_supplement_blocks,
                :missing_mixed_gausslet_gto_blocks,
                :missing_gto_gto_blocks,
                :missing_combined_raw_moment_matrices,
                :missing_residual_mwg_representation,
                :missing_combined_density_density_readiness,
            )
    else
        status = :blocked_pqs_physical_gausslet_supplement_preflight
        blocker = :unsupported_physical_gausslet_supplement_policy
        missing_fact_labels = (blocker,)
    end
    support_counts =
        target_available ?
        _pqs_source_box_route_driver_ordered_count_tuple(
            target_payload.support_counts,
            target_payload.support_units,
        ) :
        nothing
    retained_counts =
        target_available ?
        _pqs_source_box_route_driver_ordered_count_tuple(
            target_payload.retained_counts,
            target_payload.retained_order,
        ) :
        nothing
    retained_order = target_available ? target_payload.retained_order : ()
    gausslet_final_dimension =
        target_available ? target_payload.expected_final_dimension : nothing
    supplement_policy =
        target_available ? target_payload.supplement_policy : :not_available
    summary = (;
        status,
        blocker,
        route_family =
            target_available ? target_payload.route_family : :not_available,
        route_kind =
            target_available ? target_payload.route_kind : :not_available,
        fixture_label,
        support_counts,
        retained_counts,
        retained_order,
        retained_transform_kind,
        gausslet_final_dimension,
        supplement_policy,
        required_fact_labels,
        available_fact_labels,
        missing_fact_labels,
        matrices_materialized = false,
        supplemented_values_materialized = false,
        public_api = false,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_preflight_payload,
        boundary = :metadata_readiness_only,
        gto_mwg_materialization = false,
    )
    return _PQSDiatomicPhysicalGaussletSupplementPreflightPayload(
        status,
        blocker,
        summary.route_family,
        summary.route_kind,
        fixture_label,
        support_counts,
        retained_counts,
        retained_order,
        retained_transform_kind,
        gausslet_final_dimension,
        supplement_policy,
        required_fact_labels,
        available_fact_labels,
        missing_fact_labels,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_h2_wl_gausslet_only_reference_candidate(
    parent,
    route_skeleton,
    recipe,
    target_payload,
    final_basis_payload = nothing,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    route_kind =
        hasproperty(route_skeleton, :route_kind) ?
        route_skeleton.route_kind :
        recipe.route_kind
    final_basis_summary =
        isnothing(final_basis_payload) || !hasproperty(final_basis_payload, :summary) ?
        (;) :
        final_basis_payload.summary
    actual_final_dimension = get(final_basis_summary, :final_dimension, nothing)
    reference_label = "WL/QW H2 R=4 gausslet-only 463"
    retained_transform_kind = :white_lindsey_old_qw_gausslet_retained_transform
    expected_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
    support_counts =
        _pqs_source_box_route_driver_ordered_count_tuple(
            target_payload.support_counts,
            expected_order,
        )
    retained_counts =
        _pqs_source_box_route_driver_ordered_count_tuple(
            target_payload.retained_counts,
            expected_order,
        )
    final_dimension =
        isnothing(actual_final_dimension) ?
        target_payload.expected_final_dimension :
        actual_final_dimension
    conditions = (;
        route_kind =
            route_family === :pqs_source_box &&
            route_kind === :bond_aligned_diatomic_physical_gausslet_core_shell_pqs,
        geometry = _pqs_source_box_route_driver_h2_r4_center_match(parent),
        parent_axis_counts =
            _pqs_source_box_route_driver_axis_counts_tuple(
                target_payload.parent_axis_counts,
            ) == (9, 9, 15),
        support_counts = support_counts == (275, 578, 362),
        retained_counts = retained_counts == (251, 98, 114),
        retained_order = target_payload.retained_order == expected_order,
        supplement_policy = target_payload.supplement_policy === :none,
        final_dimension =
            target_payload.expected_final_dimension == 463 &&
            (isnothing(actual_final_dimension) || actual_final_dimension == 463),
        retained_transform_kind =
            retained_transform_kind ===
            :white_lindsey_old_qw_gausslet_retained_transform,
        reference_label = reference_label == "WL/QW H2 R=4 gausslet-only 463",
    )
    mismatches = Tuple(key for key in keys(conditions) if !getproperty(conditions, key))
    status =
        isempty(mismatches) ?
        :available_wl_h2_gausslet_only_reference_candidate :
        :blocked_wl_h2_gausslet_only_reference_candidate
    blocker =
        isempty(mismatches) ?
        nothing :
        :wl_h2_gausslet_only_reference_candidate_mismatch
    summary = (;
        status,
        blocker,
        route_family,
        route_kind,
        geometry_label = :h2_r4,
        parent_axis_counts =
            _pqs_source_box_route_driver_axis_counts_tuple(
                target_payload.parent_axis_counts,
            ),
        support_counts,
        retained_counts,
        retained_order = target_payload.retained_order,
        supplement_policy = target_payload.supplement_policy,
        final_dimension,
        retained_transform_kind,
        label = reference_label,
        old_supplemented_scalar_references_blocked = true,
        conditions,
        mismatches,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_h2_wl_gausslet_only_reference_candidate,
        old_wl_qw_fixed_block_size = (1215, 463),
        reviewed_contract_pass = 218,
        matrix_materialized = false,
        scalar_reference_materialized = false,
    )
    return _PQSH2WLGaussletOnlyReferenceCandidatePayload(
        status,
        blocker,
        route_family,
        route_kind,
        reference_label,
        retained_transform_kind,
        target_payload.supplement_policy,
        final_dimension,
        conditions,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_candidate_payload(
    parent,
    route_skeleton,
    _recipe,
    target_payload,
)
    expected_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
    expected_support_counts = (275, 578, 362)
    expected_retained_counts = (251, 98, 114)
    expected_final_dimension = 463
    target_available =
        !isnothing(target_payload) &&
        target_payload.status === :available_physical_gausslet_core_shell_target_inventory
    parent_basis =
        hasproperty(parent, :parent_qw_basis_object) ?
        parent.parent_qw_basis_object :
        nothing
    if !target_available || isnothing(parent_basis)
        blocker =
            target_available ?
            :missing_physical_gausslet_source_plan_materializer :
            :missing_physical_gausslet_target_inventory
        summary = (;
            candidate_status = :not_available,
            candidate_source = :not_available,
            candidate_counts_match = false,
            source_plan_authority_status = :not_available,
            blocker,
        )
        return _PQSDiatomicPhysicalGaussletSourcePlanCandidatePayload(
            :not_available_physical_gausslet_source_plan_candidate,
            blocker,
            :not_available,
            nothing,
            false,
            :not_available,
            (),
            (blocker,),
            summary,
            (; source = :pqs_diatomic_physical_gausslet_source_plan_candidate),
        )
    end

    source = bond_aligned_diatomic_nested_fixed_source(parent_basis; nside = 5)
    atom_support_count = sum(
        child -> length(child.support_indices),
        source.child_sequences;
        init = 0,
    )
    if !isnothing(source.geometry.shared_midpoint_box)
        dims = Tuple(length.(source.geometry.parent_box))
        atom_support_count += length(
            _nested_box_support_indices(source.geometry.shared_midpoint_box..., dims),
        )
    end
    support_counts = (
        atom_support_count,
        (length(layer.support_indices) for layer in source.shared_shell_layers)...,
    )
    retained_counts = (
        length(source.sequence.core_column_range),
        (length(range) for range in source.sequence.layer_column_ranges)...,
    )
    final_dimension = size(source.sequence.coefficient_matrix, 2)
    counts_match =
        support_counts == expected_support_counts &&
        retained_counts == expected_retained_counts &&
        final_dimension == expected_final_dimension
    retained_order = expected_order
    order_match =
        target_payload.support_units == expected_order &&
        target_payload.retained_order == retained_order
    no_supplement = target_payload.supplement_policy === :none
    no_diagnostic_object_kind =
        !isa(source, _PQSDiatomicCompleteCoreShellSourcePlan)
    status =
        counts_match && order_match && no_supplement && no_diagnostic_object_kind ?
        :available_physical_gausslet_source_plan_candidate :
        :blocked_physical_gausslet_source_plan_candidate
    blocker =
        status === :available_physical_gausslet_source_plan_candidate ?
        :source_plan_candidate_not_route_authority :
        :physical_gausslet_source_plan_count_mismatch
    authority_status = :candidate_not_route_authority
    summary = (;
        candidate_status = status,
        candidate_source = :source_backed_fixed_source_oracle,
        candidate_counts_match = counts_match,
        source_plan_authority_status = authority_status,
        blocker,
        support_order = expected_order,
        retained_order,
        support_counts,
        retained_counts,
        final_dimension,
        no_supplement,
        no_h2_221_diagnostic_source_plan = no_diagnostic_object_kind,
    )
    return _PQSDiatomicPhysicalGaussletSourcePlanCandidatePayload(
        status,
        blocker,
        :source_backed_fixed_source_oracle,
        source,
        counts_match,
        authority_status,
        (:source_backed_fixed_source_oracle,),
        status === :available_physical_gausslet_source_plan_candidate ?
        (:physical_gausslet_source_plan_route_authority,) :
        (:physical_gausslet_source_plan_count_contract,),
        summary,
        (;
            source =
                :pqs_diatomic_physical_gausslet_source_plan_candidate,
            route_owned_authority = false,
            diagnostic_221_source_plan_reused = false,
        ),
    )
end

function _pqs_source_box_route_driver_physical_gausslet_source_plan_from_candidate(
    candidate_payload,
)
    source = candidate_payload.candidate
    shared_shell_support_indices =
        Tuple(Vector{Int}(layer.support_indices) for layer in source.shared_shell_layers)
    shared_shell_support_states =
        Tuple(layer.support_states for layer in source.shared_shell_layers)
    retained_ranges = (
        atom_contact_core = source.sequence.core_column_range,
        shared_shell_1 = source.sequence.layer_column_ranges[1],
        shared_shell_2 = source.sequence.layer_column_ranges[2],
    )
    core_coefficient_matrix =
        source.sequence.coefficient_matrix[:, retained_ranges.atom_contact_core]
    shared_shell_coefficient_matrices = (
        source.sequence.coefficient_matrix[:, retained_ranges.shared_shell_1],
        source.sequence.coefficient_matrix[:, retained_ranges.shared_shell_2],
    )
    support_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
    retained_order = support_order
    support_counts = (
        length(source.sequence.core_indices),
        length(shared_shell_support_indices[1]),
        length(shared_shell_support_indices[2]),
    )
    retained_counts = (
        length(retained_ranges.atom_contact_core),
        length(retained_ranges.shared_shell_1),
        length(retained_ranges.shared_shell_2),
    )
    final_dimension = size(source.sequence.coefficient_matrix, 2)
    convention_labels = (;
        source_plan_family = :fake_pqs_source_backed_fixed_source_adapter,
        source_backed_adapter = true,
        source_backed_candidate_source = candidate_payload.candidate_source,
        retained_transform_kind = :wl_qw_source_backed_retained_transform,
        independent_pqs_transform = false,
        fake_pqs_source = :source_backed_fixed_source_oracle,
        fake_pqs_warning =
            :retained_transform_imported_from_wl_qw_fixed_source_oracle,
        route_owned_authority = true,
        supplement_policy = :none,
        h2_221_diagnostic_source_plan_reused = false,
        final_basis_materialized = false,
        h1_materialized = false,
        h1_j_materialized = false,
        rhf_materialized = false,
        exports_materialized = false,
        public_api = false,
    )
    summary = (;
        object_kind = :pqs_diatomic_physical_gausslet_core_shell_source_plan,
        status = :available_pqs_diatomic_physical_gausslet_core_shell_source_plan,
        support_order,
        retained_order,
        support_counts,
        retained_counts,
        final_dimension,
        retained_ranges,
        source_plan_authority_status = :fake_pqs_private_source_backed_adapter_authority,
        source_backed_candidate_source = candidate_payload.candidate_source,
        source_backed_adapter = true,
        retained_transform_kind = :wl_qw_source_backed_retained_transform,
        independent_pqs_transform = false,
        route_owned_authority = true,
        supplement_policy = :none,
    )
    return _PQSDiatomicPhysicalGaussletCoreShellSourcePlan(
        :pqs_diatomic_physical_gausslet_core_shell_source_plan,
        :available_pqs_diatomic_physical_gausslet_core_shell_source_plan,
        source.basis,
        source.axis_bundles,
        Vector{Int}(source.sequence.core_indices),
        source.sequence.core_states,
        shared_shell_support_indices,
        shared_shell_support_states,
        core_coefficient_matrix,
        shared_shell_coefficient_matrices,
        support_order,
        retained_order,
        retained_ranges,
        final_dimension,
        convention_labels,
        summary,
        (;
            source =
                :pqs_source_box_route_driver_physical_gausslet_source_plan_from_candidate,
            source_backed_adapter = true,
        ),
    )
end

function _pqs_source_box_route_driver_independent_h2_physical_source_plan_descriptor(
    target_payload,
)
    descriptor_blocker = :missing_independent_pqs_source_plan_numerical_materialization
    blocked(blocker) = (;
        object_kind = :pqs_independent_h2_physical_source_plan_descriptor,
        status = :blocked_independent_pqs_physical_source_plan_descriptor,
        blocker,
        source_plan_family = :independent_pqs_physical_source_box_descriptor,
        source_plan_authority_status =
            :blocked_independent_pqs_route_owned_source_plan_descriptor,
        source_coefficients_materialized = false,
        final_basis_materialized = false,
        source_backed_fixed_source_oracle_used = false,
        fake_pqs_enabled = false,
        missing_objects = (blocker,),
    )
    isnothing(target_payload) &&
        return blocked(:missing_independent_pqs_target_payload)
    target_payload.source_plan_role === :independent_pqs_source_box_construction ||
        return blocked(:not_independent_pqs_source_box_construction)
    target_payload.status === :available_physical_gausslet_core_shell_target_inventory ||
        return blocked(target_payload.blocker)
    support_plan = get(target_payload.summary, :support_plan, nothing)
    retained_rule_plan = get(target_payload.summary, :retained_rule_plan, nothing)
    isnothing(support_plan) &&
        return blocked(:missing_independent_pqs_support_region_plan)
    isnothing(retained_rule_plan) &&
        return blocked(:missing_independent_pqs_retained_rule_plan)
    support_counts = target_payload.support_counts
    retained_counts = target_payload.retained_counts
    Tuple(values(support_counts)) == (275, 578, 362) ||
        return blocked(:independent_pqs_source_plan_support_count_mismatch)
    Tuple(values(retained_counts)) == (275, 98, 98) ||
        return blocked(:independent_pqs_source_plan_retained_count_mismatch)

    shared_source(role) = (;
        source_family = :pqs_filled_source_cpb,
        source_mode_dims = retained_rule_plan.shared_shell_source_mode_dims,
        source_mode_count = prod(retained_rule_plan.shared_shell_source_mode_dims),
        support_count = getproperty(support_counts, role),
        transform_kind = getproperty(retained_rule_plan.per_unit_transform_kind, role),
        coefficient_matrix_materialized = false,
    )
    retained_rule(role, rule, materialized) = (;
        retained_rule = rule,
        retained_count = getproperty(retained_counts, role),
        retained_column_selector_materialized = materialized,
    )
    descriptor_summary = (;
        source_plan_family = :independent_pqs_physical_source_box_descriptor,
        source_plan_authority_status =
            :independent_pqs_route_owned_source_plan_descriptor,
        support_counts,
        retained_counts,
        final_dimension = target_payload.expected_final_dimension,
        source_coefficients_materialized = false,
        final_basis_materialized = false,
    )
    return (;
        object_kind = :pqs_independent_h2_physical_source_plan_descriptor,
        status = :available_independent_pqs_physical_source_plan_descriptor,
        blocker = nothing,
        source_plan_family = :independent_pqs_physical_source_box_descriptor,
        source_plan_authority_status =
            :independent_pqs_route_owned_source_plan_descriptor,
        support_order = target_payload.support_units,
        retained_order = target_payload.retained_order,
        support_counts,
        retained_counts,
        final_dimension = target_payload.expected_final_dimension,
        expected_final_dimension = target_payload.expected_final_dimension,
        unit_source_descriptors = (;
            atom_contact_core = (;
                source_family = :direct_terminal_source_modes,
                source_region_roles =
                    (:atom_local_core, :atom_local_core, :midpoint_slab),
                support_count = support_counts.atom_contact_core,
                source_mode_count = support_counts.atom_contact_core,
                transform_kind = :identity_source_modes,
                coefficient_matrix_materialized = false,
            ),
            shared_shell_1 = shared_source(:shared_shell_1),
            shared_shell_2 = shared_source(:shared_shell_2),
        ),
        unit_retained_rule_descriptors = (;
            atom_contact_core =
                retained_rule(:atom_contact_core, :direct_source_modes, false),
            shared_shell_1 = retained_rule(
                :shared_shell_1,
                retained_rule_plan.shared_shell_retained_rule_kind,
                true,
            ),
            shared_shell_2 = retained_rule(
                :shared_shell_2,
                retained_rule_plan.shared_shell_retained_rule_kind,
                true,
            ),
        ),
        source_coefficients_materialized = false,
        final_basis_materialized = false,
        source_backed_fixed_source_oracle_used = false,
        fake_pqs_enabled = false,
        next_blocker = descriptor_blocker,
        support_plan,
        retained_rule_plan,
        missing_objects = (descriptor_blocker,),
        summary = descriptor_summary,
        metadata = (;
            source =
                :pqs_source_box_route_driver_independent_h2_physical_source_plan_descriptor,
            descriptor_only = true,
        ),
    )
end

function _pqs_source_box_route_driver_independent_h2_shared_shell_realization_payload(
    parent,
    target_payload,
    descriptor,
)
    blocked(blocker) = (;
        object_kind = :pqs_independent_h2_shared_shell_realization_payload,
        status = :blocked_independent_pqs_shared_shell_realization_payload,
        blocker,
        shared_shell_realization_counts = (),
        shared_shell_realization_identity_errors = (),
        shared_shell_realization_materialized = false,
        source_backed_fixed_source_oracle_used = false,
        fake_pqs_enabled = false,
        shared_shells = (;),
        summary = (;
            status = :blocked_independent_pqs_shared_shell_realization_payload,
            blocker,
            shared_shell_realization_counts = (),
            shared_shell_realization_identity_errors = (),
            shared_shell_realization_materialized = false,
        ),
    )
    isnothing(descriptor) &&
        return blocked(:missing_independent_pqs_source_plan_descriptor)
    descriptor.status === :available_independent_pqs_physical_source_plan_descriptor ||
        return blocked(descriptor.blocker)
    support_plan = get(target_payload.summary, :support_plan, nothing)
    retained_rule_plan = get(target_payload.summary, :retained_rule_plan, nothing)
    isnothing(support_plan) &&
        return blocked(:missing_independent_pqs_support_region_plan)
    isnothing(retained_rule_plan) &&
        return blocked(:missing_independent_pqs_retained_rule_plan)
    hasproperty(support_plan, :shared_shell_descriptors) ||
        return blocked(:missing_independent_pqs_shared_shell_support_descriptors)
    bundles =
        hasproperty(parent, :parent_axis_bundle_object) ?
        parent.parent_axis_bundle_object :
        nothing
    isnothing(bundles) && return blocked(:missing_parent_axis_bundle_object)
    metrics = _pqs_multilayer_axis_metrics(bundles)
    source_mode_dims = retained_rule_plan.shared_shell_source_mode_dims
    shells = map(
        role -> _pqs_source_box_route_driver_independent_h2_shared_shell_realization(
            role,
            getproperty(support_plan.shared_shell_descriptors, role),
            bundles,
            metrics,
            source_mode_dims,
            support_plan.bond_axis,
        ),
        (:shared_shell_1, :shared_shell_2),
    )
    blocker = findfirst(shell -> shell.status !==
                                 :available_independent_pqs_shared_shell_realization,
                        shells)
    isnothing(blocker) || return blocked(shells[blocker].blocker)
    counts = Tuple(shell.retained_count for shell in shells)
    identity_errors = Tuple(shell.realized_overlap_identity_error for shell in shells)
    summary = (;
        status = :available_independent_pqs_shared_shell_realization_payload,
        blocker = nothing,
        shared_shell_realization_counts = counts,
        shared_shell_realization_identity_errors = identity_errors,
        coefficient_shapes = Tuple(shell.coefficient_shape for shell in shells),
        shell_projection_materialized = true,
        lowdin_cleanup_materialized = true,
        shared_shell_realization_materialized = true,
        source_backed_fixed_source_oracle_used = false,
        fake_pqs_enabled = false,
    )
    return (;
        object_kind = :pqs_independent_h2_shared_shell_realization_payload,
        status = summary.status,
        blocker = nothing,
        shared_shell_realization_counts = counts,
        shared_shell_realization_identity_errors = identity_errors,
        shared_shell_realization_materialized = true,
        source_backed_fixed_source_oracle_used = false,
        fake_pqs_enabled = false,
        shared_shells = (; shared_shell_1 = shells[1], shared_shell_2 = shells[2]),
        summary,
    )
end

function _pqs_source_box_route_driver_independent_h2_shared_shell_realization(
    role::Symbol,
    shell_descriptor,
    bundles,
    metrics,
    source_mode_dims::NTuple{3,Int},
    bond_axis::Symbol,
)
    blocked(blocker) = (;
        role,
        status = :blocked_independent_pqs_shared_shell_realization,
        blocker,
        retained_count = nothing,
        coefficient_shape = nothing,
        realized_overlap_identity_error = nothing,
    )
    try
        raw_plan = CartesianRawProductSources.raw_product_box_plan(
            shell_descriptor.source_cpb;
            source_mode_dims,
            source_key = Symbol("independent_h2_", role, "_q5_source"),
            metadata = (; route_unit_role = role),
        )
        retained_rule =
            CartesianRawProductSources.pqs_boundary_product_mode_retained_rule(
                raw_plan;
                metadata = (; route_unit_role = role),
            )
        retained_rule.retained_count == 98 ||
            return blocked(:independent_pqs_shared_shell_retained_count_mismatch)
        layer = _nested_projected_q_shell_layer(
            bundles,
            shell_descriptor.current_box,
            shell_descriptor.inner_box;
            bond_axis,
            q = source_mode_dims[1],
            L = source_mode_dims[1],
            raw_source_dims = source_mode_dims,
            selected_q = source_mode_dims[1],
        )
        projected_descriptor = _nested_projected_q_shell_staged_unit_descriptor(layer)
        projected_descriptor.support_count == shell_descriptor.support_count ||
            return blocked(:independent_pqs_shared_shell_support_count_mismatch)
        shell_plan =
            CartesianContractedParentMetrics._pqs_shell_realization_plan(
                projected_descriptor,
                metrics,
            )
        coefficients = shell_plan.shell_projection_matrix * shell_plan.lowdin_cleanup
        size(coefficients) == (shell_descriptor.support_count, 98) ||
            return blocked(:independent_pqs_shared_shell_coefficient_shape_mismatch)
        return (;
            role,
            status = :available_independent_pqs_shared_shell_realization,
            blocker = nothing,
            raw_source_plan = raw_plan,
            retained_rule,
            support_indices = shell_descriptor.support_indices,
            support_states = shell_descriptor.support_states,
            shell_projection = shell_plan.shell_projection_matrix,
            lowdin_cleanup = shell_plan.lowdin_cleanup,
            shell_final_coefficients = coefficients,
            retained_count = size(coefficients, 2),
            coefficient_shape = size(coefficients),
            realized_overlap_identity_error = shell_plan.isometry_error,
            shell_projection_materialized = true,
            lowdin_cleanup_materialized = true,
            source_backed_fixed_source_oracle_used = false,
            fake_pqs_enabled = false,
        )
    catch err
        return blocked(Symbol(:independent_pqs_shared_shell_realization_error))
    end
end

function _pqs_source_box_route_driver_parent_row_sparse_coefficients(
    support_indices::AbstractVector{Int},
    local_coefficients::AbstractMatrix{<:Real},
    parent_count::Int,
)
    row_indices = Int[]
    col_indices = Int[]
    values = Float64[]
    for column in axes(local_coefficients, 2)
        for (local_row, parent_row) in pairs(support_indices)
            value = local_coefficients[local_row, column]
            iszero(value) && continue
            push!(row_indices, Int(parent_row))
            push!(col_indices, Int(column))
            push!(values, Float64(value))
        end
    end
    return _nested_sparse_coefficient_map(
        row_indices,
        col_indices,
        values,
        parent_count,
        size(local_coefficients, 2),
    )
end

function _pqs_source_box_route_driver_independent_h2_complete_core_shell_source_plan(
    parent,
    target_payload,
    descriptor,
    shared_shell_realization,
)
    blocked(blocker) = (;
        status = :blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan,
        blocker,
        source_plan = nothing,
    )
    isnothing(target_payload) && return blocked(:missing_independent_pqs_target_payload)
    descriptor.status === :available_independent_pqs_physical_source_plan_descriptor ||
        return blocked(descriptor.blocker)
    shared_shell_realization.status ===
        :available_independent_pqs_shared_shell_realization_payload ||
        return blocked(shared_shell_realization.blocker)
    support_plan = get(target_payload.summary, :support_plan, nothing)
    isnothing(support_plan) &&
        return blocked(:missing_independent_pqs_support_region_plan)
    hasproperty(support_plan, :atom_contact_core_descriptor) ||
        return blocked(:missing_independent_pqs_atom_contact_core_descriptor)
    bundles =
        hasproperty(parent, :parent_axis_bundle_object) ?
        parent.parent_axis_bundle_object :
        nothing
    isnothing(bundles) && return blocked(:missing_parent_axis_bundle_object)

    support_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
    retained_order = support_order
    support_counts = (275, 578, 362)
    retained_counts = (275, 98, 98)
    target_payload.support_units == support_order ||
        return blocked(:independent_pqs_source_plan_support_order_mismatch)
    target_payload.retained_order == retained_order ||
        return blocked(:independent_pqs_source_plan_retained_order_mismatch)
    Tuple(values(target_payload.support_counts)) == support_counts ||
        return blocked(:independent_pqs_source_plan_support_count_mismatch)
    Tuple(values(target_payload.retained_counts)) == retained_counts ||
        return blocked(:independent_pqs_source_plan_retained_count_mismatch)

    core = support_plan.atom_contact_core_descriptor
    core.support_count == retained_counts[1] ||
        return blocked(:independent_pqs_atom_contact_core_count_mismatch)
    shell_1 = shared_shell_realization.shared_shells.shared_shell_1
    shell_2 = shared_shell_realization.shared_shells.shared_shell_2
    parent_dims = _pqs_source_box_route_driver_axis_counts_tuple(
        target_payload.parent_axis_counts,
    )
    isnothing(parent_dims) &&
        return blocked(:missing_independent_pqs_parent_axis_counts)
    parent_count = prod(parent_dims)
    core_indices = Vector{Int}(core.support_indices)
    core_coefficient_matrix = _nested_sparse_coefficient_map(
        core_indices,
        collect(1:length(core_indices)),
        ones(Float64, length(core_indices)),
        parent_count,
        length(core_indices),
    )
    shared_shell_coefficient_matrices = (
        _pqs_source_box_route_driver_parent_row_sparse_coefficients(
            shell_1.support_indices,
            shell_1.shell_final_coefficients,
            parent_count,
        ),
        _pqs_source_box_route_driver_parent_row_sparse_coefficients(
            shell_2.support_indices,
            shell_2.shell_final_coefficients,
            parent_count,
        ),
    )
    retained_ranges = (;
        atom_contact_core = 1:retained_counts[1],
        shared_shell_1 = (retained_counts[1] + 1):sum(retained_counts[1:2]),
        shared_shell_2 = (sum(retained_counts[1:2]) + 1):sum(retained_counts),
    )
    convention_labels = (;
        source_plan_family = :independent_pqs_physical_source_box_core_shell,
        source_backed_adapter = false,
        source_backed_fixed_source_oracle_used = false,
        retained_transform_kind = :pqs_source_box_retained_transform,
        independent_pqs_transform = true,
        fake_pqs_enabled = false,
        retained_transform_authority = :pqs_source_box_construction,
        source_plan_authority_status = :independent_pqs_route_owned_source_plan,
        route_owned_authority = true,
        core_coefficient_representation = :sparse_parent_row_direct_selector,
        shared_shell_coefficient_representation = :sparse_parent_row_realization,
        final_basis_materialized = false,
        h1_materialized = false,
        h1_j_materialized = false,
        rhf_materialized = false,
        exports_materialized = false,
        public_api = false,
    )
    summary = (;
        object_kind = :pqs_diatomic_physical_gausslet_core_shell_source_plan,
        status = :available_pqs_diatomic_physical_gausslet_core_shell_source_plan,
        blocker = nothing,
        support_order,
        retained_order,
        support_counts,
        retained_counts,
        final_dimension = sum(retained_counts),
        retained_ranges,
        source_plan_family = :independent_pqs_physical_source_box_core_shell,
        source_plan_authority_status = :independent_pqs_route_owned_source_plan,
        source_backed_fixed_source_oracle_used = false,
        fake_pqs_enabled = false,
        retained_transform_authority = :pqs_source_box_construction,
        source_coefficients_materialized = true,
        final_basis_materialized = false,
        shared_shell_realization_counts =
            shared_shell_realization.shared_shell_realization_counts,
        shared_shell_realization_identity_errors =
            shared_shell_realization.shared_shell_realization_identity_errors,
    )
    parent_basis =
        hasproperty(parent, :parent_qw_basis_object) ? parent.parent_qw_basis_object :
        hasproperty(parent, :parent_basis_object) ? parent.parent_basis_object :
        nothing
    source_plan = _PQSDiatomicPhysicalGaussletCoreShellSourcePlan(
        :pqs_diatomic_physical_gausslet_core_shell_source_plan,
        :available_pqs_diatomic_physical_gausslet_core_shell_source_plan,
        parent_basis,
        bundles,
        core_indices,
        core.support_states,
        (Vector{Int}(shell_1.support_indices), Vector{Int}(shell_2.support_indices)),
        (shell_1.support_states, shell_2.support_states),
        core_coefficient_matrix,
        shared_shell_coefficient_matrices,
        support_order,
        retained_order,
        retained_ranges,
        sum(retained_counts),
        convention_labels,
        summary,
        (;
            source =
                :pqs_source_box_route_driver_independent_h2_complete_core_shell_source_plan,
            route_owned = true,
            source_backed_adapter = false,
        ),
    )
    return (;
        status = source_plan.status,
        blocker = nothing,
        source_plan,
    )
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_payload(
    parent,
    target_payload,
    candidate_payload = nothing,
)
    target_available =
        !isnothing(target_payload) &&
        target_payload.status === :available_physical_gausslet_core_shell_target_inventory
    independent_target =
        !isnothing(target_payload) &&
        target_payload.source_plan_role === :independent_pqs_source_box_construction
    source_plan =
        !isnothing(candidate_payload) &&
        candidate_payload.status === :available_physical_gausslet_source_plan_candidate ?
        _pqs_source_box_route_driver_physical_gausslet_source_plan_from_candidate(
            candidate_payload,
        ) :
        nothing
    descriptor =
        independent_target ?
        _pqs_source_box_route_driver_independent_h2_physical_source_plan_descriptor(
            target_payload,
        ) :
        nothing
    descriptor_available =
        !isnothing(descriptor) &&
        descriptor.status === :available_independent_pqs_physical_source_plan_descriptor
    shared_shell_realization =
        descriptor_available ?
        _pqs_source_box_route_driver_independent_h2_shared_shell_realization_payload(
            parent,
            target_payload,
            descriptor,
        ) :
        nothing
    shared_shell_realization_available =
        !isnothing(shared_shell_realization) &&
        shared_shell_realization.status ===
        :available_independent_pqs_shared_shell_realization_payload
    independent_source_plan =
        shared_shell_realization_available ?
        _pqs_source_box_route_driver_independent_h2_complete_core_shell_source_plan(
            parent,
            target_payload,
            descriptor,
            shared_shell_realization,
        ) :
        nothing
    independent_source_plan_available =
        !isnothing(independent_source_plan) &&
        independent_source_plan.status ===
        :available_pqs_diatomic_physical_gausslet_core_shell_source_plan
    source_plan =
        independent_source_plan_available ?
        independent_source_plan.source_plan :
        descriptor_available ? nothing : source_plan
    independent_source_plan_blocker =
        !isnothing(target_payload) ?
        get(
            target_payload.summary,
            :independent_source_plan_blocker,
            :missing_independent_pqs_physical_source_plan_materializer,
        ) :
        :missing_independent_pqs_physical_source_plan_materializer
    status =
        independent_source_plan_available ?
        source_plan.status :
        descriptor_available ?
        :blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan :
        isnothing(source_plan) ?
        :blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan :
        source_plan.status
    blocker =
        independent_source_plan_available ?
        nothing :
        descriptor_available ?
        shared_shell_realization_available ?
        independent_source_plan.blocker :
        shared_shell_realization.blocker :
        !isnothing(source_plan) ?
        nothing :
        independent_target ?
        independent_source_plan_blocker :
        target_available ?
        :missing_atom_contact_core_support_rows :
        :missing_physical_gausslet_target_inventory
    missing_objects =
        independent_source_plan_available ?
        () :
        descriptor_available ?
        (blocker,) :
        !isnothing(source_plan) ?
        () :
        independent_target ?
        (independent_source_plan_blocker,) :
        target_available ?
        (
            :atom_contact_core_support_rows,
            :shared_shell_1_coefficients,
            :shared_shell_2_coefficients,
        ) :
        (blocker,)
    summary = (;
        object_kind = :pqs_diatomic_physical_gausslet_core_shell_source_plan,
        status,
        blocker,
        target_status = isnothing(target_payload) ? :not_available : target_payload.status,
        parent_axis_counts = isnothing(target_payload) ? nothing : target_payload.parent_axis_counts,
        support_order = isnothing(target_payload) ? () : target_payload.support_units,
        retained_order = isnothing(target_payload) ? () : target_payload.retained_order,
        support_counts = isnothing(target_payload) ? (;) : target_payload.support_counts,
        retained_counts = isnothing(target_payload) ? (;) : target_payload.retained_counts,
        expected_final_dimension =
            isnothing(target_payload) ? nothing : target_payload.expected_final_dimension,
        retained_atom_core_interiors =
            !isnothing(target_payload) && target_payload.retained_atom_core_interiors,
        source_plan_role =
            isnothing(target_payload) ? :not_available : target_payload.source_plan_role,
        supplement_policy =
            isnothing(target_payload) ? :not_available : target_payload.supplement_policy,
        source_plan_materialized = !isnothing(source_plan),
        source_plan_descriptor_status =
            isnothing(descriptor) ? :not_available : descriptor.status,
        source_plan_descriptor_blocker =
            isnothing(descriptor) ? nothing : descriptor.blocker,
        shared_shell_realization_status =
            isnothing(shared_shell_realization) ?
            :not_available :
            shared_shell_realization.status,
        shared_shell_realization_blocker =
            isnothing(shared_shell_realization) ?
            nothing :
            shared_shell_realization.blocker,
        shared_shell_realization_counts =
            isnothing(shared_shell_realization) ?
            () :
            shared_shell_realization.shared_shell_realization_counts,
        shared_shell_realization_identity_errors =
            isnothing(shared_shell_realization) ?
            () :
            shared_shell_realization.shared_shell_realization_identity_errors,
        shared_shell_realization_materialized = shared_shell_realization_available,
        source_plan_descriptor_available = descriptor_available,
        source_plan_family =
            independent_source_plan_available ?
            source_plan.summary.source_plan_family :
            isnothing(source_plan) ?
            :not_available :
            hasproperty(source_plan, :source_plan_family) ?
            source_plan.source_plan_family :
            :not_available,
        source_plan_authority_status =
            independent_source_plan_available ?
            source_plan.summary.source_plan_authority_status :
            descriptor_available ?
            descriptor.source_plan_authority_status :
            isnothing(source_plan) ?
            independent_target ?
            :blocked_pqs_source_box_construction_authority :
            isnothing(candidate_payload) ? :not_available : candidate_payload.authority_status :
            :fake_pqs_private_source_backed_adapter_authority,
        source_backed_fixed_source_oracle_used =
            !isnothing(target_payload) &&
            get(target_payload.summary, :source_backed_fixed_source_oracle_used, false),
        retained_transform_authority =
            independent_source_plan_available ?
            source_plan.summary.retained_transform_authority :
            isnothing(target_payload) ?
            :not_available :
            get(target_payload.summary, :retained_transform_authority, :not_available),
        source_coefficients_materialized =
            independent_source_plan_available,
        secondary_blocker =
            isnothing(target_payload) ? nothing : get(target_payload.summary, :secondary_blocker, nothing),
        independent_source_plan_blocker =
            isnothing(target_payload) ?
            nothing :
            get(target_payload.summary, :independent_source_plan_blocker, nothing),
        missing_objects,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_payload,
        route_owned = true,
        diagnostic_221_source_plan_reused = false,
        placeholders_synthesized = false,
    )

    return _PQSDiatomicPhysicalGaussletCoreShellSourcePlanPayload(
        status,
        blocker,
        isnothing(target_payload) ? :not_available : target_payload.route_family,
        source_plan,
        (),
        missing_objects,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_support_partition_rows(
    source_plan,
    unit_key::Symbol,
)
    unit_key === :atom_contact_core &&
        return Vector{Int}(source_plan.atom_contact_core_support_indices)
    unit_key === :shared_shell_1 &&
        return Vector{Int}(source_plan.shared_shell_support_indices[1])
    unit_key === :shared_shell_2 &&
        return Vector{Int}(source_plan.shared_shell_support_indices[2])
    return Int[]
end

function _pqs_source_box_route_driver_support_partition_local_rows(
    support_indices::AbstractVector{Int},
    tile_indices::AbstractVector{Int},
)
    local_by_parent = Dict{Int,Int}(
        Int(parent_row) => Int(local_row) for (local_row, parent_row) in pairs(support_indices)
    )
    return Tuple(local_by_parent[Int(parent_row)] for parent_row in tile_indices)
end

function _pqs_source_box_route_driver_contiguous_unit_range(rows::Tuple)
    isempty(rows) && return nothing
    first_row = first(rows)
    expected = first_row:(first_row + length(rows) - 1)
    return Tuple(expected) == rows ? expected : nothing
end

function _pqs_source_box_route_driver_support_partition_unit(
    unit_key::Symbol,
    support_indices::AbstractVector{Int},
    retained_range,
    support_tiles::Tuple,
)
    observed = Int[]
    tile_fingerprints = map(support_tiles) do tile
        tile_rows = Vector{Int}(tile.support_indices)
        append!(observed, tile_rows)
        unit_local_rows =
            _pqs_source_box_route_driver_support_partition_local_rows(
                support_indices,
                tile_rows,
            )
        unit_local_row_range =
            _pqs_source_box_route_driver_contiguous_unit_range(unit_local_rows)
        (;
            tile_key = tile.tile_key,
            unit_key,
            source_region_role = tile.source_region_role,
            support_tile_role = tile.support_tile_role,
            cpb_role = tile.cpb_role,
            intervals = tile.intervals,
            support_count = tile.support_count,
            parent_row_count = length(tile_rows),
            parent_row_min = tile.parent_row_min,
            parent_row_max = tile.parent_row_max,
            parent_row_ids = Tuple(tile_rows),
            unit_local_row_range,
            unit_local_rows,
            retained_range,
            provider_tile_ready = tile.provider_tile_ready,
        )
    end
    expected_set = Set(Int.(support_indices))
    observed_set = Set(observed)
    duplicate_count = length(observed) - length(observed_set)
    missing_count = length(setdiff(expected_set, observed_set))
    outside_count = length(setdiff(observed_set, expected_set))
    coverage_complete =
        duplicate_count == 0 && missing_count == 0 && outside_count == 0 &&
        length(observed_set) == length(expected_set)
    return (;
        unit_key,
        support_count = length(support_indices),
        retained_range,
        tile_count = length(support_tiles),
        tile_support_counts = Tuple(tile.support_count for tile in support_tiles),
        parent_row_count = length(observed_set),
        parent_row_min = isempty(support_indices) ? nothing : minimum(support_indices),
        parent_row_max = isempty(support_indices) ? nothing : maximum(support_indices),
        duplicate_parent_row_count = duplicate_count,
        missing_parent_row_count = missing_count,
        outside_parent_row_count = outside_count,
        coverage_complete,
        tiles = support_tiles,
        tile_fingerprints = Tuple(tile_fingerprints),
    )
end

function _pqs_source_box_route_driver_partition_summary_unit(partition)
    return (;
        unit_key = partition.unit_key,
        support_count = partition.support_count,
        retained_range = partition.retained_range,
        tile_count = partition.tile_count,
        tile_support_counts = partition.tile_support_counts,
        parent_row_count = partition.parent_row_count,
        parent_row_min = partition.parent_row_min,
        parent_row_max = partition.parent_row_max,
        duplicate_parent_row_count = partition.duplicate_parent_row_count,
        missing_parent_row_count = partition.missing_parent_row_count,
        outside_parent_row_count = partition.outside_parent_row_count,
        coverage_complete = partition.coverage_complete,
        tile_fingerprints = partition.tile_fingerprints,
    )
end

function _pqs_source_box_route_driver_independent_h2_pqs_supplement_support_partition_payload(
    target_payload,
    source_plan_payload,
)
    support_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
    blocked(blocker) = _PQSIndependentH2PQSSupplementSupportPartitionPayload(
        :blocked_independent_h2_pqs_supplement_support_partition,
        blocker,
        isnothing(target_payload) ? :not_available : target_payload.route_family,
        isnothing(target_payload) ? :not_available : target_payload.route_kind,
        support_order,
        (),
        nothing,
        nothing,
        nothing,
        (),
        0,
        0,
        false,
        0,
        0,
        0,
        false,
        false,
        :not_available,
        (; status = :blocked_independent_h2_pqs_supplement_support_partition, blocker),
        (;
            source =
                :pqs_source_box_route_driver_independent_h2_pqs_supplement_support_partition_payload,
            provider_blocks_materialized = false,
            route_global_matrix_materialized = false,
        ),
    )
    isnothing(target_payload) && return blocked(:missing_independent_h2_target_payload)
    target_payload.status === :available_physical_gausslet_core_shell_target_inventory ||
        return blocked(target_payload.blocker)
    support_plan = get(target_payload.summary, :support_plan, nothing)
    isnothing(support_plan) && return blocked(:missing_independent_pqs_support_region_plan)
    source_plan =
        !isnothing(source_plan_payload) && hasproperty(source_plan_payload, :source_plan) ?
        source_plan_payload.source_plan :
        nothing
    isnothing(source_plan) && return blocked(:missing_independent_pqs_source_plan)
    source_plan.status ===
        :available_pqs_diatomic_physical_gausslet_core_shell_source_plan ||
        return blocked(source_plan.status)
    support_plan.support_order == support_order ||
        return blocked(:independent_pqs_support_partition_order_mismatch)
    source_plan.support_order == support_order ||
        return blocked(:independent_pqs_source_plan_support_order_mismatch)

    atom_descriptor = support_plan.atom_contact_core_descriptor
    shared_descriptors = support_plan.shared_shell_descriptors
    retained_ranges = source_plan.retained_ranges
    tile_sources = (;
        atom_contact_core = atom_descriptor.support_tiles,
        shared_shell_1 = shared_descriptors.shared_shell_1.support_tiles,
        shared_shell_2 = shared_descriptors.shared_shell_2.support_tiles,
    )
    unit_partitions = Tuple(
        _pqs_source_box_route_driver_support_partition_unit(
            unit_key,
            _pqs_source_box_route_driver_support_partition_rows(source_plan, unit_key),
            getproperty(retained_ranges, unit_key),
            getproperty(tile_sources, unit_key),
        ) for unit_key in support_order
    )
    duplicate_count =
        sum(partition.duplicate_parent_row_count for partition in unit_partitions)
    missing_count =
        sum(partition.missing_parent_row_count for partition in unit_partitions)
    outside_count =
        sum(partition.outside_parent_row_count for partition in unit_partitions)
    coverage_complete =
        all(partition.coverage_complete for partition in unit_partitions) &&
        duplicate_count == 0 && missing_count == 0 && outside_count == 0
    total_tile_count = sum(partition.tile_count for partition in unit_partitions)
    total_support_count = sum(partition.support_count for partition in unit_partitions)
    status =
        coverage_complete ?
        :available_independent_h2_pqs_supplement_support_partition :
        :blocked_independent_h2_pqs_supplement_support_partition
    blocker =
        coverage_complete ? nothing : :independent_h2_pqs_support_partition_coverage_mismatch
    support_counts =
        _pqs_source_box_route_driver_ordered_count_tuple(
            target_payload.support_counts,
            support_order,
        )
    retained_counts =
        _pqs_source_box_route_driver_ordered_count_tuple(
            target_payload.retained_counts,
            source_plan.retained_order,
        )
    summary = (;
        status,
        blocker,
        route_family = target_payload.route_family,
        route_kind = target_payload.route_kind,
        support_order,
        retained_order = source_plan.retained_order,
        support_counts,
        retained_counts,
        retained_ranges,
        unit_partitions =
            Tuple(_pqs_source_box_route_driver_partition_summary_unit(partition)
                  for partition in unit_partitions),
        total_tile_count,
        total_support_count,
        coverage_complete,
        duplicate_parent_row_count = duplicate_count,
        missing_parent_row_count = missing_count,
        outside_parent_row_count = outside_count,
        source_backed_fixed_source_oracle_used = false,
        fake_pqs_enabled = false,
        retained_transform_authority = :pqs_source_box_construction,
        provider_blocks_materialized = false,
        route_global_matrix_materialized = false,
    )
    return _PQSIndependentH2PQSSupplementSupportPartitionPayload(
        status,
        blocker,
        target_payload.route_family,
        target_payload.route_kind,
        support_order,
        source_plan.retained_order,
        support_counts,
        retained_counts,
        retained_ranges,
        unit_partitions,
        total_tile_count,
        total_support_count,
        coverage_complete,
        duplicate_count,
        missing_count,
        outside_count,
        false,
        false,
        :pqs_source_box_construction,
        summary,
        (;
            source =
                :pqs_source_box_route_driver_independent_h2_pqs_supplement_support_partition_payload,
            provider_blocks_materialized = false,
            route_global_matrix_materialized = false,
        ),
    )
end

function _pqs_source_box_route_driver_physical_gausslet_axis_metrics(axis_bundles)
    isnothing(axis_bundles) && return nothing
    pgdg_x = _nested_axis_pgdg(axis_bundles, :x)
    pgdg_y = _nested_axis_pgdg(axis_bundles, :y)
    pgdg_z = _nested_axis_pgdg(axis_bundles, :z)
    return (;
        x = (; overlap = pgdg_x.overlap, kinetic = pgdg_x.kinetic,
            source = :nested_pgdg_axis),
        y = (; overlap = pgdg_y.overlap, kinetic = pgdg_y.kinetic,
            source = :nested_pgdg_axis),
        z = (; overlap = pgdg_z.overlap, kinetic = pgdg_z.kinetic,
            source = :nested_pgdg_axis),
    )
end

function _pqs_source_box_route_driver_physical_gausslet_local_coefficients(
    coefficient_matrix,
    support_indices,
    retained_count::Int,
    blocker::Symbol,
)
    coefficients = Matrix{Float64}(coefficient_matrix)
    support = Vector{Int}(support_indices)
    shape = (length(support), retained_count)
    size(coefficients, 2) == retained_count || return (; blocker, coefficients = nothing)
    isempty(support) && return (; blocker, coefficients = nothing)
    maximum(support) <= size(coefficients, 1) || return (; blocker, coefficients = nothing)
    local_coefficients = coefficients[support, :]
    size(local_coefficients) == shape || return (; blocker, coefficients = nothing)
    all(isfinite, local_coefficients) || return (;
        blocker = :physical_gausslet_coefficient_nonfinite,
        coefficients = nothing,
    )
    return (; blocker = nothing, coefficients = local_coefficients)
end

function _pqs_source_box_route_driver_physical_gausslet_overlap_diagnostics(
    overlap::AbstractMatrix{<:Real};
    identity_atol::Real,
    rank_atol::Real,
)
    matrix = Matrix{Float64}(overlap)
    size(matrix, 1) == size(matrix, 2) || return (;
        blocker = :physical_gausslet_final_overlap_not_square,
        summary = (;),
        cleanup = nothing,
        final_overlap = nothing,
        final_identity_error = nothing,
    )
    symmetric_matrix = Symmetric((matrix + transpose(matrix)) ./ 2)
    symmetry_error = norm(matrix - transpose(matrix), Inf)
    decomposition = eigen(symmetric_matrix)
    eigenvalues = decomposition.values
    threshold = max(Float64(rank_atol), eps(Float64) * max(size(matrix, 1), 1))
    rank = count(value -> value > threshold, eigenvalues)
    summary = (;
        shape = size(matrix),
        symmetry_error,
        eigenvalue_min = minimum(eigenvalues),
        eigenvalue_max = maximum(eigenvalues),
        rank,
        expected_dimension = size(matrix, 1),
        rank_atol = Float64(rank_atol),
        identity_atol = Float64(identity_atol),
        full_rank = rank == size(matrix, 1),
    )
    rank == size(matrix, 1) || return (;
        blocker = :physical_gausslet_final_overlap_rank_deficient,
        summary,
        cleanup = nothing,
        final_overlap = nothing,
        final_identity_error = nothing,
    )
    cleanup =
        decomposition.vectors * Diagonal(1.0 ./ sqrt.(decomposition.values))
    final_overlap = transpose(cleanup) * matrix * cleanup
    final_identity_error = norm(
        final_overlap -
        Matrix{Float64}(I, size(final_overlap, 1), size(final_overlap, 2)),
        Inf,
    )
    blocker =
        final_identity_error <= Float64(identity_atol) ?
        nothing :
        :physical_gausslet_final_overlap_not_identity
    return (; blocker, summary, cleanup, final_overlap, final_identity_error)
end

function _pqs_source_box_route_driver_physical_gausslet_final_basis(
    source_plan;
    identity_atol::Real = 1.0e-8,
    rank_atol::Real = 1.0e-10,
)
    support_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
    source_plan.support_order == support_order ||
        return (; status = :blocked_pqs_physical_gausslet_final_basis,
            blocker = :physical_gausslet_support_order_mismatch)
    source_plan.retained_order == support_order ||
        return (; status = :blocked_pqs_physical_gausslet_final_basis,
            blocker = :physical_gausslet_retained_order_mismatch)

    support_indices = (
        source_plan.atom_contact_core_support_indices,
        source_plan.shared_shell_support_indices[1],
        source_plan.shared_shell_support_indices[2],
    )
    support_states = (
        source_plan.atom_contact_core_support_states,
        source_plan.shared_shell_support_states[1],
        source_plan.shared_shell_support_states[2],
    )
    coefficient_matrices = (
        source_plan.core_coefficient_matrix,
        source_plan.shared_shell_coefficient_matrices[1],
        source_plan.shared_shell_coefficient_matrices[2],
    )
    retained_counts = (
        length(source_plan.retained_ranges.atom_contact_core),
        length(source_plan.retained_ranges.shared_shell_1),
        length(source_plan.retained_ranges.shared_shell_2),
    )
    support_counts = map(length, support_indices)
    support_counts == (275, 578, 362) ||
        return (; status = :blocked_pqs_physical_gausslet_final_basis,
            blocker = :physical_gausslet_support_count_mismatch)
    expected_retained_counts =
        get(source_plan.summary, :retained_counts, retained_counts)
    expected_retained_counts =
        expected_retained_counts isa Tuple ?
        Tuple(Int(count) for count in expected_retained_counts) :
        expected_retained_counts
    retained_counts == expected_retained_counts ||
        return (; status = :blocked_pqs_physical_gausslet_final_basis,
            blocker = :physical_gausslet_retained_count_mismatch)

    local_blocks = Matrix{Float64}[]
    for (index, matrix) in pairs(coefficient_matrices)
        blocker =
            index == 1 ?
            :physical_gausslet_core_coefficient_shape_mismatch :
            :physical_gausslet_shared_shell_coefficient_shape_mismatch
        local_coefficients =
            _pqs_source_box_route_driver_physical_gausslet_local_coefficients(
                matrix,
                support_indices[index],
                retained_counts[index],
                blocker,
            )
        isnothing(local_coefficients.blocker) || return (;
            status = :blocked_pqs_physical_gausslet_final_basis,
            blocker = local_coefficients.blocker,
        )
        push!(local_blocks, local_coefficients.coefficients)
    end

    metrics =
        _pqs_source_box_route_driver_physical_gausslet_axis_metrics(
            source_plan.axis_bundles,
        )
    isnothing(metrics) && return (;
        status = :blocked_pqs_physical_gausslet_final_basis,
        blocker = :missing_physical_gausslet_support_overlap,
    )

    combined_support_states = vcat(
        Vector{NTuple{3,Int}}(support_states[1]),
        Vector{NTuple{3,Int}}(support_states[2]),
        Vector{NTuple{3,Int}}(support_states[3]),
    )
    support_overlap = _pqs_multilayer_support_product_matrix(
        combined_support_states,
        combined_support_states,
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    total_support_count = sum(support_counts)
    final_dimension = sum(retained_counts)
    pre_final_coefficients = zeros(Float64, total_support_count, final_dimension)
    support_row_ranges = (
        atom_contact_core = 1:support_counts[1],
        shared_shell_1 = (support_counts[1] + 1):(support_counts[1] + support_counts[2]),
        shared_shell_2 =
            (support_counts[1] + support_counts[2] + 1):total_support_count,
    )
    retained_ranges = (
        atom_contact_core = 1:retained_counts[1],
        shared_shell_1 = (retained_counts[1] + 1):(retained_counts[1] + retained_counts[2]),
        shared_shell_2 =
            (retained_counts[1] + retained_counts[2] + 1):final_dimension,
    )
    pre_final_coefficients[support_row_ranges.atom_contact_core, retained_ranges.atom_contact_core] .=
        local_blocks[1]
    pre_final_coefficients[support_row_ranges.shared_shell_1, retained_ranges.shared_shell_1] .=
        local_blocks[2]
    pre_final_coefficients[support_row_ranges.shared_shell_2, retained_ranges.shared_shell_2] .=
        local_blocks[3]
    pre_final_overlap =
        transpose(pre_final_coefficients) * support_overlap * pre_final_coefficients
    pre_final_identity_error = norm(
        pre_final_overlap -
        Matrix{Float64}(I, final_dimension, final_dimension),
        Inf,
    )
    diagnostics =
        _pqs_source_box_route_driver_physical_gausslet_overlap_diagnostics(
            pre_final_overlap;
            identity_atol,
            rank_atol,
        )
    isnothing(diagnostics.blocker) || return (;
        object_kind = :pqs_physical_gausslet_final_basis,
        status = :blocked_pqs_physical_gausslet_final_basis,
        blocker = diagnostics.blocker,
        final_retained_count = final_dimension,
        final_overlap_identity_error = diagnostics.final_identity_error,
        pre_final_overlap_identity_error = pre_final_identity_error,
        support_counts,
        retained_counts,
        overlap_diagnostics = diagnostics.summary,
    )
    final_coefficients = pre_final_coefficients * diagnostics.cleanup
    return (;
        object_kind = :pqs_physical_gausslet_final_basis,
        status = :available_pqs_physical_gausslet_final_basis,
        blocker = nothing,
        source_plan_object_kind = source_plan.object_kind,
        support_order,
        retained_order = support_order,
        support_counts,
        retained_counts,
        support_row_ranges,
        retained_ranges,
        final_retained_count = final_dimension,
        final_dimension,
        pre_final_coefficients,
        combined_lowdin_cleanup = diagnostics.cleanup,
        final_coefficients,
        pre_final_overlap_identity_error = pre_final_identity_error,
        final_overlap_identity_error = diagnostics.final_identity_error,
        final_overlap_is_identity = true,
        overlap_diagnostics = diagnostics.summary,
        transform_source_plan_provenance = source_plan.convention_labels,
        source_backed_fixed_source_oracle_used =
            get(source_plan.convention_labels, :source_backed_fixed_source_oracle_used, false),
        fake_pqs_enabled =
            get(source_plan.convention_labels, :fake_pqs_enabled, false),
        retained_transform_authority =
            get(
                source_plan.convention_labels,
                :retained_transform_authority,
                :not_available,
            ),
        final_basis_materialized = true,
        h1_materialized = false,
        h1_j_materialized = false,
        density_interaction_materialized = false,
        rhf_materialized = false,
        exports_materialized = false,
        public_api = false,
        metadata = (;
            source = :pqs_source_box_route_driver_physical_gausslet_final_basis,
            support_decomposition = :shared_physical_gausslet_core_shell,
            source_backed_adapter =
                get(source_plan.convention_labels, :source_backed_adapter, false),
            route_owned_authority = true,
            diagnostic_materialization = true,
            self_overlap_stored_for_downstream = false,
        ),
    )
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_final_basis_payload(
    route_skeleton,
    recipe,
    source_plan_payload = nothing,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    source_plan =
        !isnothing(source_plan_payload) && hasproperty(source_plan_payload, :source_plan) ?
        source_plan_payload.source_plan :
        nothing
    source_plan_status =
        isnothing(source_plan) ? :not_available : source_plan.status
    final_basis_requested =
        get(recipe, :run_final_basis, false) ||
        get(recipe, :run_h1, false) ||
        get(recipe, :run_h1_j, false) ||
        get(get(recipe, :private_rhf_inputs, (;)), :run_private_rhf, false)
    final_basis = nothing
    final_basis_status = :not_materialized_pqs_physical_gausslet_final_basis
    available = Symbol[]
    missing = Symbol[]

    if route_family !== :pqs_source_box
        status = :not_applicable_pqs_physical_gausslet_final_basis_non_pqs_route
        blocker = nothing
    elseif source_plan_status !==
           :available_pqs_diatomic_physical_gausslet_core_shell_source_plan
        status = :blocked_pqs_physical_gausslet_final_basis_payload
        blocker = :missing_pqs_diatomic_physical_gausslet_source_plan
        push!(missing, :pqs_diatomic_physical_gausslet_source_plan)
    elseif !final_basis_requested
        status = :not_requested_pqs_physical_gausslet_final_basis_payload
        blocker = :physical_gausslet_final_basis_request_not_enabled
        push!(available, :pqs_diatomic_physical_gausslet_source_plan)
        push!(missing, :physical_gausslet_final_basis_request)
    else
        final_basis =
            _pqs_source_box_route_driver_physical_gausslet_final_basis(
                source_plan,
            )
        final_basis_status = final_basis.status
        push!(available, :pqs_diatomic_physical_gausslet_source_plan)
        if final_basis_status === :available_pqs_physical_gausslet_final_basis
            status = :available_pqs_physical_gausslet_final_basis_payload
            blocker = nothing
            push!(available, :pqs_physical_gausslet_final_basis)
            push!(missing, :physical_gausslet_h1_builder)
        else
            status = :blocked_pqs_physical_gausslet_final_basis_payload
            blocker =
                isnothing(final_basis.blocker) ?
                :physical_gausslet_final_basis_blocked :
                final_basis.blocker
            push!(missing, :pqs_physical_gausslet_final_basis)
        end
    end

    available_objects = Tuple(unique(available))
    missing_objects = Tuple(unique(missing))
    support_counts =
        isnothing(source_plan) ? nothing : source_plan.summary.support_counts
    retained_counts =
        isnothing(source_plan) ? nothing : source_plan.summary.retained_counts
    final_dimension =
        !isnothing(final_basis) && hasproperty(final_basis, :final_dimension) ?
        final_basis.final_dimension :
        isnothing(source_plan) ? nothing : source_plan.final_dimension
    final_overlap_identity_error =
        !isnothing(final_basis) &&
        hasproperty(final_basis, :final_overlap_identity_error) ?
        final_basis.final_overlap_identity_error :
        nothing
    overlap_diagnostics =
        !isnothing(final_basis) &&
        hasproperty(final_basis, :overlap_diagnostics) ?
        final_basis.overlap_diagnostics :
        (;)
    summary = (;
        status,
        blocker,
        route_family,
        source_plan_status,
        final_basis_status,
        support_order =
            isnothing(source_plan) ? () : source_plan.support_order,
        retained_order =
            isnothing(source_plan) ? () : source_plan.retained_order,
        support_counts,
        retained_counts,
        final_dimension,
        final_overlap_identity_error,
        pre_final_overlap_identity_error =
            !isnothing(final_basis) &&
            hasproperty(final_basis, :pre_final_overlap_identity_error) ?
            final_basis.pre_final_overlap_identity_error :
            nothing,
        final_overlap_rank = get(overlap_diagnostics, :rank, nothing),
        final_overlap_full_rank = get(overlap_diagnostics, :full_rank, nothing),
        final_overlap_eigenvalue_min =
            get(overlap_diagnostics, :eigenvalue_min, nothing),
        final_overlap_eigenvalue_max =
            get(overlap_diagnostics, :eigenvalue_max, nothing),
        source_backed_fixed_source_oracle_used =
            !isnothing(final_basis) &&
            hasproperty(final_basis, :source_backed_fixed_source_oracle_used) ?
            final_basis.source_backed_fixed_source_oracle_used :
            false,
        fake_pqs_enabled =
            !isnothing(final_basis) &&
            hasproperty(final_basis, :fake_pqs_enabled) ?
            final_basis.fake_pqs_enabled :
            false,
        retained_transform_authority =
            !isnothing(final_basis) &&
            hasproperty(final_basis, :retained_transform_authority) ?
            final_basis.retained_transform_authority :
            isnothing(source_plan) ?
            :not_available :
            get(source_plan.summary, :retained_transform_authority, :not_available),
        final_basis_materialized =
            final_basis_status === :available_pqs_physical_gausslet_final_basis,
        endpoint_blocker =
            final_basis_status === :available_pqs_physical_gausslet_final_basis ?
            :missing_physical_gausslet_h1_builder :
            blocker,
        h1_materialized = false,
        h1_j_materialized = false,
        density_interaction_materialized = false,
        rhf_materialized = false,
        exports_materialized = false,
        public_api = false,
        available_objects,
        missing_objects,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_physical_gausslet_final_basis_payload,
        route_kind = recipe.route_kind,
        route_family,
        source_plan_status,
        final_basis_status,
        support_decomposition = :shared_physical_gausslet_core_shell,
        h2_221_diagnostic_source_plan_reused = false,
    )
    return _PQSDiatomicPhysicalGaussletFinalBasisPayload(
        status,
        blocker,
        route_family,
        source_plan,
        source_plan_status,
        final_basis,
        final_basis_status,
        available_objects,
        missing_objects,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_physical_gausslet_support_states(source_plan)
    return vcat(
        Vector{NTuple{3,Int}}(source_plan.atom_contact_core_support_states),
        Vector{NTuple{3,Int}}(source_plan.shared_shell_support_states[1]),
        Vector{NTuple{3,Int}}(source_plan.shared_shell_support_states[2]),
    )
end

function _pqs_source_box_route_driver_physical_gausslet_support_kinetic(source_plan)
    states = _pqs_source_box_route_driver_physical_gausslet_support_states(source_plan)
    metrics =
        _pqs_source_box_route_driver_physical_gausslet_axis_metrics(
            source_plan.axis_bundles,
        )
    support_operator =
        _pqs_multilayer_support_product_matrix(
            states,
            states,
            metrics.x.kinetic,
            metrics.y.overlap,
            metrics.z.overlap,
        ) +
        _pqs_multilayer_support_product_matrix(
            states,
            states,
            metrics.x.overlap,
            metrics.y.kinetic,
            metrics.z.overlap,
        ) +
        _pqs_multilayer_support_product_matrix(
            states,
            states,
            metrics.x.overlap,
            metrics.y.overlap,
            metrics.z.kinetic,
        )
    return (;
        object_kind = :pqs_physical_gausslet_support_kinetic_matrix,
        status = :materialized_pqs_physical_gausslet_support_kinetic_matrix,
        blocker = nothing,
        term = :kinetic,
        support_operator,
        support_operator_shape = size(support_operator),
        support_operator_finite = all(isfinite, support_operator),
        support_state_count = length(states),
        support_order = source_plan.support_order,
    )
end

function _pqs_source_box_route_driver_physical_gausslet_support_electron_nuclear(
    source_plan;
    coulomb_expansion,
    center_records,
    axis_layers,
)
    centers = _pqs_multilayer_center_records(center_records)
    coefficients = Float64.(coulomb_expansion.coefficients)
    states = _pqs_source_box_route_driver_physical_gausslet_support_states(source_plan)
    records = map(enumerate(centers)) do (_center_index, center_record)
        center = _pqs_multilayer_center_summary(center_record)
        factor_x, factor_y, factor_z, factor_source =
            _pqs_multilayer_centered_factor_terms(
                axis_layers,
                coulomb_expansion,
                center,
            )
        axis_terms =
            _pqs_multilayer_validate_factor_terms(
                (factor_x, factor_y, factor_z),
                length(coefficients),
            )
        support_operator =
            _pqs_multilayer_support_electron_nuclear_matrix(
                states,
                axis_terms,
                coefficients,
            )
        (;
            object_kind =
                :pqs_physical_gausslet_support_electron_nuclear_by_center_matrix,
            status =
                :materialized_pqs_physical_gausslet_support_electron_nuclear_by_center_matrix,
            blocker = nothing,
            term = :electron_nuclear_by_center,
            support_operator,
            support_operator_shape = size(support_operator),
            support_operator_finite = all(isfinite, support_operator),
            support_state_count = length(states),
            support_order = source_plan.support_order,
            center_key = center.center_key,
            center_index = center.center_index,
            center_location = center.location,
            nuclear_charge = center.nuclear_charge,
            nuclear_charge_recorded = true,
            nuclear_charge_applied = false,
            centers_summed = false,
            uncharged_by_center_convention = true,
            gaussian_factor_terms_source = factor_source,
            metadata = (;
                source =
                    :pqs_source_box_route_driver_physical_gausslet_support_electron_nuclear,
                physical_operator = :electron_nuclear_attraction,
                nuclear_charge = center.nuclear_charge,
                nuclear_charge_recorded = true,
                nuclear_charge_applied = false,
                centers_summed = false,
                center_key = center.center_key,
                center_index = center.center_index,
                center_location = center.location,
            ),
        )
    end
    return (;
        object_kind =
            :pqs_physical_gausslet_support_electron_nuclear_by_center_matrix_set,
        status =
            :materialized_pqs_physical_gausslet_support_electron_nuclear_by_center_matrix_set,
        blocker = nothing,
        term = :electron_nuclear_by_center,
        records,
        center_count = length(records),
        support_state_count = length(states),
        support_order = source_plan.support_order,
    )
end

function _pqs_source_box_route_driver_physical_gausslet_final_one_body(
    final_basis,
    support_operator;
    term::Symbol,
    center_record = nothing,
    symmetry_atol::Real = 1.0e-8,
)
    operator = Matrix{Float64}(support_operator)
    coefficients = Matrix{Float64}(final_basis.final_coefficients)
    size(operator) == (size(coefficients, 1), size(coefficients, 1)) ||
        throw(DimensionMismatch("physical gausslet support one-body shape mismatch"))
    all(isfinite, operator) ||
        throw(ArgumentError("physical gausslet support one-body matrix is nonfinite"))
    support_symmetry_error = norm(operator - transpose(operator), Inf)
    support_symmetry_error <= Float64(symmetry_atol) ||
        throw(ArgumentError("physical gausslet support one-body matrix is not symmetric"))
    final_operator = transpose(coefficients) * operator * coefficients
    final_symmetry_error = norm(final_operator - transpose(final_operator), Inf)
    final_symmetry_error <= Float64(symmetry_atol) ||
        throw(ArgumentError("physical gausslet final one-body matrix is not symmetric"))
    center_metadata =
        isnothing(center_record) ?
        (;) :
        (;
            center_key = center_record.center_key,
            center_index = center_record.center_index,
            center_location = center_record.center_location,
            nuclear_charge = center_record.nuclear_charge,
            nuclear_charge_recorded = true,
            nuclear_charge_applied = false,
            centers_summed = false,
            uncharged_by_center_convention = true,
        )
    return (;
        object_kind = :pqs_physical_gausslet_final_one_body_matrix,
        status = :materialized_pqs_physical_gausslet_final_one_body_matrix,
        blocker = nothing,
        term,
        final_basis_object_kind = final_basis.object_kind,
        final_basis_status = final_basis.status,
        support_order = final_basis.support_order,
        support_operator_shape = size(operator),
        support_operator_finite = all(isfinite, operator),
        support_operator_symmetry_error = support_symmetry_error,
        final_operator,
        final_operator_shape = size(final_operator),
        final_operator_finite = all(isfinite, final_operator),
        final_operator_symmetry_error = final_symmetry_error,
        final_retained_count = final_basis.final_retained_count,
        metadata = merge(
            (;
                source =
                    :pqs_source_box_route_driver_physical_gausslet_final_one_body,
                physical_final_basis = true,
            ),
            center_metadata,
        ),
    )
end

function _pqs_source_box_route_driver_physical_gausslet_h1_hamiltonian(
    final_kinetic,
    final_nuclear_by_center;
    symmetry_atol::Real = 1.0e-8,
)
    kinetic = Matrix{Float64}(final_kinetic.final_operator)
    dimension = size(kinetic, 1)
    charged_nuclear_sum = zeros(Float64, dimension, dimension)
    for result in Tuple(final_nuclear_by_center)
        nuclear = Matrix{Float64}(result.final_operator)
        size(nuclear) == (dimension, dimension) ||
            throw(DimensionMismatch("physical gausslet nuclear one-body shape mismatch"))
        metadata = result.metadata
        charge = Float64(metadata.nuclear_charge)
        charged_nuclear_sum .+= charge .* nuclear
    end
    hamiltonian_matrix = kinetic + charged_nuclear_sum
    symmetry_error = norm(hamiltonian_matrix - transpose(hamiltonian_matrix), Inf)
    symmetry_error <= Float64(symmetry_atol) ||
        throw(ArgumentError("physical gausslet H1 Hamiltonian is not symmetric"))
    return (;
        object_kind = :pqs_physical_gausslet_final_one_electron_hamiltonian,
        status = :materialized_pqs_physical_gausslet_final_one_electron_hamiltonian,
        blocker = nothing,
        final_dimension = dimension,
        kinetic_matrix = kinetic,
        charged_nuclear_matrix = charged_nuclear_sum,
        hamiltonian_matrix,
        hamiltonian_matrix_shape = size(hamiltonian_matrix),
        hamiltonian_matrix_finite = all(isfinite, hamiltonian_matrix),
        hamiltonian_matrix_symmetry_error = symmetry_error,
        center_count = length(Tuple(final_nuclear_by_center)),
    )
end

function _pqs_source_box_route_driver_physical_gausslet_h1_solve(final_hamiltonian)
    matrix = Matrix{Float64}(final_hamiltonian.hamiltonian_matrix)
    all(isfinite, matrix) ||
        throw(ArgumentError("physical gausslet H1 matrix contains non-finite entries"))
    solution = eigen(Symmetric((matrix + transpose(matrix)) ./ 2))
    eigenvalues = solution.values
    energy = first(eigenvalues)
    return (;
        object_kind = :pqs_physical_gausslet_h1_solve,
        status = :materialized_pqs_physical_gausslet_h1_solve,
        blocker = nothing,
        solve_kind = :ordinary_symmetric,
        final_dimension = size(matrix, 1),
        lowest_energy = energy,
        lowest_orbital_coefficients = Vector{Float64}(solution.vectors[:, 1]),
        eigenvalue_min = minimum(eigenvalues),
        eigenvalue_max = maximum(eigenvalues),
        h1_solve_materialized = true,
    )
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_h1_payload(
    parent,
    route_skeleton,
    recipe,
    source_plan_payload = nothing,
    final_basis_payload = nothing,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    source_plan =
        !isnothing(source_plan_payload) && hasproperty(source_plan_payload, :source_plan) ?
        source_plan_payload.source_plan :
        nothing
    final_basis =
        !isnothing(final_basis_payload) && hasproperty(final_basis_payload, :final_basis) ?
        final_basis_payload.final_basis :
        nothing
    source_plan_status = isnothing(source_plan) ? :not_available : source_plan.status
    final_basis_status = isnothing(final_basis) ? :not_available : final_basis.status
    h1_requested =
        get(recipe, :run_h1, false) ||
        get(recipe, :run_h1_j, false) ||
        get(get(recipe, :private_rhf_inputs, (;)), :run_private_rhf, false)
    support_kinetic = support_nuclear = final_kinetic = final_nuclear = nothing
    final_hamiltonian = h1 = nothing
    support_kinetic_status = :not_materialized_pqs_physical_gausslet_support_kinetic
    support_nuclear_status = :not_materialized_pqs_physical_gausslet_support_electron_nuclear
    final_kinetic_status = :not_materialized_pqs_physical_gausslet_final_kinetic
    final_nuclear_status = :not_materialized_pqs_physical_gausslet_final_electron_nuclear
    final_hamiltonian_status = :not_materialized_pqs_physical_gausslet_h1_hamiltonian
    h1_status = :not_materialized_pqs_physical_gausslet_h1_solve
    available = Symbol[]
    missing = Symbol[]

    if route_family !== :pqs_source_box
        status = :not_applicable_pqs_physical_gausslet_h1_non_pqs_route
        blocker = nothing
    elseif source_plan_status !==
           :available_pqs_diatomic_physical_gausslet_core_shell_source_plan
        status = :blocked_pqs_physical_gausslet_h1_payload
        blocker = :missing_pqs_diatomic_physical_gausslet_source_plan
        push!(missing, :pqs_diatomic_physical_gausslet_source_plan)
    elseif final_basis_status !== :available_pqs_physical_gausslet_final_basis
        status = :blocked_pqs_physical_gausslet_h1_payload
        blocker = :missing_pqs_physical_gausslet_final_basis
        push!(available, :pqs_diatomic_physical_gausslet_source_plan)
        push!(missing, :pqs_physical_gausslet_final_basis)
    elseif !h1_requested
        status = :not_requested_pqs_physical_gausslet_h1_payload
        blocker = :physical_gausslet_h1_request_not_enabled
        push!(available, :pqs_diatomic_physical_gausslet_source_plan)
        push!(available, :pqs_physical_gausslet_final_basis)
        push!(missing, :physical_gausslet_h1_request)
    else
        push!(available, :pqs_diatomic_physical_gausslet_source_plan)
        push!(available, :pqs_physical_gausslet_final_basis)
        center_records, missing_centers =
            _pqs_source_box_route_driver_complete_core_shell_center_records(parent)
        axis_layers, missing_axis_layers =
            _pqs_source_box_route_driver_complete_core_shell_axis_layers(
                source_plan.axis_bundles,
            )
        if !isempty((missing_centers..., missing_axis_layers...))
            status = :blocked_pqs_physical_gausslet_h1_payload
            blocker = :missing_physical_gausslet_h1_inputs
            append!(missing, (missing_centers..., missing_axis_layers...))
        else
            try
                coulomb_expansion = coulomb_gaussian_expansion(doacc = false)
                support_kinetic =
                    _pqs_source_box_route_driver_physical_gausslet_support_kinetic(
                        source_plan,
                    )
                support_kinetic_status = support_kinetic.status
                support_nuclear =
                    _pqs_source_box_route_driver_physical_gausslet_support_electron_nuclear(
                        source_plan;
                        coulomb_expansion,
                        center_records,
                        axis_layers,
                    )
                support_nuclear_status = support_nuclear.status
                final_kinetic =
                    _pqs_source_box_route_driver_physical_gausslet_final_one_body(
                        final_basis,
                        support_kinetic.support_operator;
                        term = :kinetic,
                    )
                final_kinetic_status = final_kinetic.status
                final_nuclear = map(support_nuclear.records) do record
                    _pqs_source_box_route_driver_physical_gausslet_final_one_body(
                        final_basis,
                        record.support_operator;
                        term = :electron_nuclear_by_center,
                        center_record = record,
                    )
                end
                final_nuclear_status =
                    :materialized_pqs_physical_gausslet_final_electron_nuclear_by_center
                final_hamiltonian =
                    _pqs_source_box_route_driver_physical_gausslet_h1_hamiltonian(
                        final_kinetic,
                        final_nuclear,
                    )
                final_hamiltonian_status = final_hamiltonian.status
                h1 = _pqs_source_box_route_driver_physical_gausslet_h1_solve(
                    final_hamiltonian,
                )
                h1_status = h1.status
                if !(isfinite(h1.lowest_energy) && h1.lowest_energy < 0)
                    status = :blocked_pqs_physical_gausslet_h1_payload
                    blocker = :physical_gausslet_h1_lowest_energy_nonnegative
                    push!(missing, :physical_gausslet_h1_sanity)
                else
                    status = :available_pqs_physical_gausslet_h1_payload
                    blocker = nothing
                    push!(available, :pqs_physical_gausslet_h1_payload)
                    push!(missing, :physical_gausslet_h1_j_builder)
                end
            catch error
                error isa ArgumentError || error isa DimensionMismatch || rethrow()
                status = :blocked_pqs_physical_gausslet_h1_payload
                blocker = :physical_gausslet_h1_builder_error
                h1_status = :blocked_pqs_physical_gausslet_h1_solve
                push!(missing, :pqs_physical_gausslet_h1_payload)
            end
        end
    end

    available_objects = Tuple(unique(available))
    missing_objects = Tuple(unique(missing))
    h1_materialized = status === :available_pqs_physical_gausslet_h1_payload
    summary = (;
        status,
        blocker,
        route_family,
        source_plan_status,
        final_basis_status,
        support_kinetic_status,
        support_electron_nuclear_status = support_nuclear_status,
        final_kinetic_status,
        final_electron_nuclear_status = final_nuclear_status,
        final_hamiltonian_status,
        h1_status,
        final_dimension =
            !isnothing(h1) ? h1.final_dimension :
            !isnothing(final_basis) ? final_basis.final_dimension : nothing,
        lowest_energy = isnothing(h1) ? nothing : h1.lowest_energy,
        h1_hamiltonian_matrix_finite =
            isnothing(final_hamiltonian) ? nothing :
            final_hamiltonian.hamiltonian_matrix_finite,
        h1_hamiltonian_symmetry_error =
            isnothing(final_hamiltonian) ? nothing :
            final_hamiltonian.hamiltonian_matrix_symmetry_error,
        center_count = isnothing(support_nuclear) ? nothing : support_nuclear.center_count,
        h1_materialized,
        endpoint_blocker =
            h1_materialized ? :missing_physical_gausslet_h1_j_builder : blocker,
        h1_j_materialized = false,
        density_interaction_materialized = false,
        rhf_materialized = false,
        exports_materialized = false,
        public_api = false,
        available_objects,
        missing_objects,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_physical_gausslet_h1_payload,
        route_kind = recipe.route_kind,
        physical_final_basis = true,
        h2_221_diagnostic_source_plan_reused = false,
    )
    return _PQSDiatomicPhysicalGaussletH1Payload(
        status,
        blocker,
        route_family,
        source_plan,
        source_plan_status,
        final_basis,
        final_basis_status,
        support_kinetic,
        support_kinetic_status,
        support_nuclear,
        support_nuclear_status,
        final_kinetic,
        final_kinetic_status,
        final_nuclear,
        final_nuclear_status,
        final_hamiltonian,
        final_hamiltonian_status,
        h1,
        h1_status,
        available_objects,
        missing_objects,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_physical_gausslet_density_provenance(
    source_plan,
    coulomb_expansion,
)
    expected_term_count = length(coulomb_expansion.coefficients)
    provenance =
        CartesianContractedParentMetrics._pqs_source_box_ida_factor_provenance(
            source_plan.axis_bundles;
            expected_term_count,
        )
    hasproperty(provenance, :axis_weights) ||
        throw(ArgumentError("physical gausslet density provenance missing axis weights"))
    hasproperty(provenance, :raw_axis_pair_factor_terms) ||
        throw(ArgumentError("physical gausslet density provenance missing raw pair factor terms"))
    return (;
        object_kind = :pqs_physical_gausslet_density_provenance,
        status = :available_pqs_physical_gausslet_density_provenance,
        blocker = nothing,
        provenance,
        axis_weights = provenance.axis_weights,
        raw_pair_factor_terms = provenance.raw_axis_pair_factor_terms,
        term_count = provenance.term_count,
        factor_dimensions = provenance.factor_dimensions,
        metadata = (;
            source =
                :pqs_source_box_route_driver_physical_gausslet_density_provenance,
            provenance_source = :pqs_source_box_ida_factor_provenance,
            expected_term_count,
            retained_diagnostic_weights_are_ida_weights = false,
        ),
    )
end

function _pqs_source_box_route_driver_physical_gausslet_support_weights(
    source_plan;
    axis_weights,
)
    states = _pqs_source_box_route_driver_physical_gausslet_support_states(source_plan)
    weights_x, weights_y, weights_z =
        _pqs_multilayer_common_or_axis_tuple(axis_weights, "axis_weights")
    weights = Vector{Float64}(undef, length(states))
    @inbounds for (index, (ix, iy, iz)) in pairs(states)
        weights[index] =
            Float64(weights_x[ix]) * Float64(weights_y[iy]) * Float64(weights_z[iz])
    end
    return (;
        object_kind = :pqs_physical_gausslet_support_weights,
        status = :materialized_pqs_physical_gausslet_support_weights,
        blocker = nothing,
        support_weights = weights,
        support_weight_count = length(weights),
        support_weights_all_finite = all(isfinite, weights),
        support_weights_all_positive = all(>(0.0), weights),
        support_weight_min = minimum(weights),
        support_weight_max = maximum(weights),
        support_weight_sum = sum(weights),
        support_order = source_plan.support_order,
        metadata = (;
            source =
                :pqs_source_box_route_driver_physical_gausslet_support_weights,
            support_order = source_plan.support_order,
            density_gauge = :pre_final_localized_positive_weight,
        ),
    )
end

function _pqs_source_box_route_driver_physical_gausslet_support_pair_raw_numerator_matrix(
    source_plan;
    raw_pair_factor_terms,
    coulomb_expansion,
)
    states = _pqs_source_box_route_driver_physical_gausslet_support_states(source_plan)
    coefficients = Float64.(coulomb_expansion.coefficients)
    all(>(0.0), coefficients) ||
        throw(ArgumentError("physical gausslet raw pair numerator requires positive Coulomb expansion coefficients"))
    axis_terms =
        _pqs_multilayer_validate_raw_pair_factor_terms(
            _pqs_multilayer_raw_pair_factor_terms(raw_pair_factor_terms),
            length(coefficients),
        )
    support_pair_raw_numerator = zeros(Float64, length(states), length(states))
    for term_index in eachindex(coefficients)
        support_pair_raw_numerator .+=
            coefficients[term_index] *
            _pqs_multilayer_support_product_matrix(
                states,
                states,
                @view(axis_terms[1][term_index, :, :]),
                @view(axis_terms[2][term_index, :, :]),
                @view(axis_terms[3][term_index, :, :]),
            )
    end
    symmetry_error =
        norm(support_pair_raw_numerator - transpose(support_pair_raw_numerator), Inf)
    return (;
        object_kind = :pqs_physical_gausslet_support_raw_pair_numerator_matrix,
        status =
            :materialized_pqs_physical_gausslet_support_raw_pair_numerator_matrix,
        blocker = nothing,
        support_pair_raw_numerator,
        support_pair_raw_numerator_shape = size(support_pair_raw_numerator),
        support_pair_raw_numerator_finite = all(isfinite, support_pair_raw_numerator),
        support_pair_raw_numerator_symmetry_error = symmetry_error,
        support_order = source_plan.support_order,
        raw_pair_factor_convention = :raw_numerator,
        metadata = (;
            source =
                :pqs_source_box_route_driver_physical_gausslet_support_pair_raw_numerator_matrix,
            raw_pair_factor_convention = :raw_numerator,
        ),
    )
end

function _pqs_source_box_route_driver_physical_gausslet_pre_final_density_interaction(
    final_basis,
    support_pair_raw_numerator,
    support_weights;
    near_zero_atol::Real = 1.0e-12,
    symmetry_atol::Real = 1.0e-8,
    metadata = (;),
)
    get(final_basis, :object_kind, nothing) === :pqs_physical_gausslet_final_basis ||
        throw(ArgumentError("physical gausslet density interaction requires physical final basis"))
    get(final_basis, :status, nothing) === :available_pqs_physical_gausslet_final_basis ||
        throw(ArgumentError("physical gausslet density interaction requires available final basis"))
    get(final_basis, :final_basis_materialized, false) ||
        throw(ArgumentError("physical gausslet final basis is not materialized"))
    support_count = sum(final_basis.support_counts)
    raw_pair = Matrix{Float64}(support_pair_raw_numerator)
    size(raw_pair) == (support_count, support_count) ||
        throw(DimensionMismatch("physical gausslet raw pair numerator shape mismatch"))
    all(isfinite, raw_pair) ||
        throw(ArgumentError("physical gausslet raw pair numerator contains non-finite entries"))
    raw_pair_symmetry_error = norm(raw_pair - transpose(raw_pair), Inf)
    raw_pair_symmetry_error <= Float64(symmetry_atol) ||
        throw(ArgumentError("physical gausslet raw pair numerator must be symmetric"))

    weights = Float64[Float64(weight) for weight in support_weights]
    length(weights) == support_count ||
        throw(DimensionMismatch("physical gausslet support weight length mismatch"))
    all(isfinite, weights) ||
        throw(ArgumentError("physical gausslet support weights contain non-finite entries"))

    pre_final_coefficients = Matrix{Float64}(final_basis.pre_final_coefficients)
    final_coefficients = Matrix{Float64}(final_basis.final_coefficients)
    cleanup = Matrix{Float64}(final_basis.combined_lowdin_cleanup)
    size(pre_final_coefficients, 1) == support_count ||
        throw(DimensionMismatch("physical gausslet pre-final coefficient row mismatch"))
    size(pre_final_coefficients, 2) == final_basis.final_retained_count ||
        throw(DimensionMismatch("physical gausslet pre-final coefficient column mismatch"))

    pre_final_weights = vec(transpose(pre_final_coefficients) * weights)
    near_zero_threshold = Float64(near_zero_atol)
    near_zero_count = count(weight -> abs(weight) <= near_zero_threshold, pre_final_weights)
    negative_count = count(<(0.0), pre_final_weights)
    positive_count = count(>(0.0), pre_final_weights)
    all_finite = all(isfinite, pre_final_weights)
    positive_weight_gauge =
        all_finite && near_zero_count == 0 && negative_count == 0 &&
        positive_count == length(pre_final_weights)
    reconstruction_error =
        norm(pre_final_coefficients * cleanup - final_coefficients, Inf)

    if !positive_weight_gauge
        return (;
            object_kind = :pqs_physical_gausslet_pre_final_density_interaction,
            status =
                :blocked_pqs_physical_gausslet_pre_final_density_interaction,
            blocker = :physical_gausslet_pre_final_density_weights_not_positive,
            final_basis_object_kind = final_basis.object_kind,
            final_basis_status = final_basis.status,
            support_order = final_basis.support_order,
            density_gauge = :pre_final_localized_positive_weight,
            support_weight_count = length(weights),
            pre_final_weight_count = length(pre_final_weights),
            pre_final_weights,
            pre_final_weight_min = minimum(pre_final_weights),
            pre_final_weight_max = maximum(pre_final_weights),
            pre_final_weight_sum = sum(pre_final_weights),
            pre_final_weight_positive_count = positive_count,
            pre_final_weight_negative_count = negative_count,
            pre_final_weight_near_zero_count = near_zero_count,
            pre_final_weights_all_finite = all_finite,
            pre_final_weights_all_positive = false,
            final_to_pre_final_coefficients = cleanup,
            final_to_pre_final_reconstruction_error = reconstruction_error,
            raw_pair_numerator_shape = size(raw_pair),
            raw_pair_numerator_symmetry_error = raw_pair_symmetry_error,
            pre_final_pair_matrix = nothing,
            pre_final_pair_matrix_shape = nothing,
            pre_final_pair_matrix_finite = false,
            pre_final_pair_matrix_symmetry_error = nothing,
            metadata = merge(NamedTuple(metadata), (;
                source =
                    :pqs_source_box_route_driver_physical_gausslet_pre_final_density_interaction,
                raw_pair_factor_convention = :raw_numerator,
            )),
        )
    end

    weighted_coefficients =
        pre_final_coefficients .* reshape(1.0 ./ pre_final_weights, 1, :)
    pre_final_pair_matrix =
        transpose(weighted_coefficients) * raw_pair * weighted_coefficients
    pre_final_symmetry_error =
        norm(pre_final_pair_matrix - transpose(pre_final_pair_matrix), Inf)
    return (;
        object_kind = :pqs_physical_gausslet_pre_final_density_interaction,
        status = :materialized_pqs_physical_gausslet_pre_final_density_interaction,
        blocker = nothing,
        final_basis_object_kind = final_basis.object_kind,
        final_basis_status = final_basis.status,
        support_order = final_basis.support_order,
        density_gauge = :pre_final_localized_positive_weight,
        support_weight_count = length(weights),
        pre_final_weight_count = length(pre_final_weights),
        pre_final_weights,
        pre_final_weight_min = minimum(pre_final_weights),
        pre_final_weight_max = maximum(pre_final_weights),
        pre_final_weight_sum = sum(pre_final_weights),
        pre_final_weight_positive_count = positive_count,
        pre_final_weight_negative_count = negative_count,
        pre_final_weight_near_zero_count = near_zero_count,
        pre_final_weights_all_finite = all_finite,
        pre_final_weights_all_positive = true,
        final_to_pre_final_coefficients = cleanup,
        final_to_pre_final_reconstruction_error = reconstruction_error,
        raw_pair_numerator_shape = size(raw_pair),
        raw_pair_numerator_symmetry_error = raw_pair_symmetry_error,
        pre_final_pair_matrix,
        pre_final_pair_matrix_shape = size(pre_final_pair_matrix),
        pre_final_pair_matrix_finite = all(isfinite, pre_final_pair_matrix),
        pre_final_pair_matrix_symmetry_error = pre_final_symmetry_error,
        pre_final_weight_division_applied = true,
        signed_final_weight_division_used = false,
        raw_no_division_used = false,
        fixed_block_pair_data_authority_used = false,
        density_density_materialized = true,
        metadata = merge(NamedTuple(metadata), (;
            source =
                :pqs_source_box_route_driver_physical_gausslet_pre_final_density_interaction,
            raw_pair_factor_convention = :raw_numerator,
            weight_application_stage = :pre_final_density_interaction_boundary,
            final_orbital_consumption_rule =
                :combined_lowdin_cleanup_times_final_coefficients,
        )),
    )
end

function _pqs_source_box_route_driver_physical_gausslet_h1_j_diagnostic(
    density_interaction,
    final_orbital_coefficients,
)
    get(density_interaction, :status, nothing) ===
        :materialized_pqs_physical_gausslet_pre_final_density_interaction ||
        throw(ArgumentError("physical gausslet H1-J diagnostic requires materialized density interaction"))
    coefficients = Float64[Float64(value) for value in final_orbital_coefficients]
    length(coefficients) == size(density_interaction.final_to_pre_final_coefficients, 2) ||
        throw(DimensionMismatch("physical gausslet H1 orbital coefficient length mismatch"))
    final_to_pre_final =
        Matrix{Float64}(density_interaction.final_to_pre_final_coefficients)
    orbital = final_to_pre_final * coefficients
    pair_matrix = Matrix{Float64}(density_interaction.pre_final_pair_matrix)
    density = orbital * transpose(orbital)
    rho = 0.5 .* (density .+ transpose(density))
    v = 0.5 .* (pair_matrix .+ transpose(pair_matrix))
    occupations = vec(diag(rho))
    direct = 2.0 * dot(occupations, v * occupations)
    exchange = dot(vec(rho), vec(v .* rho))
    self_coulomb = direct - exchange
    return (;
        object_kind = :pqs_physical_gausslet_h1_j_diagnostic,
        status = :materialized_pqs_physical_gausslet_h1_j_diagnostic,
        blocker = nothing,
        density_gauge = density_interaction.density_gauge,
        raw_pair_factor_convention = :raw_numerator,
        self_coulomb,
        h1_j_self_coulomb = self_coulomb,
        final_orbital_coefficient_count = length(coefficients),
        pre_final_orbital_coefficient_count = length(orbital),
        h1_self_coulomb_contraction =
            :restricted_one_orbital_direct_minus_exchange,
    )
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_h1_j_payload(
    route_skeleton,
    recipe,
    source_plan_payload = nothing,
    final_basis_payload = nothing,
    h1_payload = nothing,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    source_plan =
        !isnothing(source_plan_payload) && hasproperty(source_plan_payload, :source_plan) ?
        source_plan_payload.source_plan :
        nothing
    final_basis =
        !isnothing(final_basis_payload) && hasproperty(final_basis_payload, :final_basis) ?
        final_basis_payload.final_basis :
        nothing
    source_plan_status = isnothing(source_plan) ? :not_available : source_plan.status
    final_basis_status = isnothing(final_basis) ? :not_available : final_basis.status
    h1_payload_status = isnothing(h1_payload) ? :not_available : h1_payload.status
    h1_j_requested =
        get(recipe, :run_h1_j, false) ||
        get(get(recipe, :private_rhf_inputs, (;)), :run_private_rhf, false)
    density_provenance = support_weights = raw_pair_factor_terms = nothing
    support_pair_raw_numerator = density_interaction = h1_j_diagnostic = nothing
    density_provenance_status =
        :not_materialized_pqs_physical_gausslet_density_provenance
    support_weights_status =
        :not_materialized_pqs_physical_gausslet_support_weights
    raw_pair_factor_status =
        :not_materialized_pqs_physical_gausslet_raw_pair_factor_terms
    support_pair_raw_numerator_status =
        :not_materialized_pqs_physical_gausslet_support_raw_pair_numerator_matrix
    density_interaction_status =
        :not_materialized_pqs_physical_gausslet_pre_final_density_interaction
    h1_j_status = :not_materialized_pqs_physical_gausslet_h1_j_payload
    available = Symbol[]
    missing = Symbol[]

    if route_family !== :pqs_source_box
        status = :not_applicable_pqs_physical_gausslet_h1_j_non_pqs_route
        blocker = nothing
    elseif source_plan_status !==
           :available_pqs_diatomic_physical_gausslet_core_shell_source_plan
        status = :blocked_pqs_physical_gausslet_h1_j_payload
        blocker = :missing_pqs_diatomic_physical_gausslet_source_plan
        push!(missing, :pqs_diatomic_physical_gausslet_source_plan)
    elseif final_basis_status !== :available_pqs_physical_gausslet_final_basis
        status = :blocked_pqs_physical_gausslet_h1_j_payload
        blocker = :missing_pqs_physical_gausslet_final_basis
        push!(missing, :pqs_physical_gausslet_final_basis)
    elseif h1_payload_status !== :available_pqs_physical_gausslet_h1_payload
        status = :blocked_pqs_physical_gausslet_h1_j_payload
        blocker = :missing_pqs_physical_gausslet_h1_payload
        push!(missing, :pqs_physical_gausslet_h1_payload)
    elseif !h1_j_requested
        status = :not_requested_pqs_physical_gausslet_h1_j_payload
        blocker = :physical_gausslet_h1_j_request_not_enabled
        append!(
            available,
            (
                :pqs_diatomic_physical_gausslet_source_plan,
                :pqs_physical_gausslet_final_basis,
                :pqs_physical_gausslet_h1_payload,
            ),
        )
        push!(missing, :physical_gausslet_h1_j_request)
    else
        append!(
            available,
            (
                :pqs_diatomic_physical_gausslet_source_plan,
                :pqs_physical_gausslet_final_basis,
                :pqs_physical_gausslet_h1_payload,
            ),
        )
        coulomb_expansion = coulomb_gaussian_expansion(doacc = false)
        try
            density_provenance =
                _pqs_source_box_route_driver_physical_gausslet_density_provenance(
                    source_plan,
                    coulomb_expansion,
                )
            density_provenance_status = density_provenance.status
            raw_pair_factor_terms = density_provenance.raw_pair_factor_terms
            raw_pair_factor_status =
                :available_pqs_physical_gausslet_raw_pair_factor_terms
            support_weights =
                _pqs_source_box_route_driver_physical_gausslet_support_weights(
                    source_plan;
                    axis_weights = density_provenance.axis_weights,
                )
            support_weights_status = support_weights.status
            support_weights.support_weights_all_positive ||
                throw(ArgumentError("physical gausslet support weights must be positive"))
            support_pair_raw_numerator =
                _pqs_source_box_route_driver_physical_gausslet_support_pair_raw_numerator_matrix(
                    source_plan;
                    raw_pair_factor_terms,
                    coulomb_expansion,
                )
            support_pair_raw_numerator_status =
                support_pair_raw_numerator.status
            density_interaction =
                _pqs_source_box_route_driver_physical_gausslet_pre_final_density_interaction(
                    final_basis,
                    support_pair_raw_numerator.support_pair_raw_numerator,
                    support_weights.support_weights;
                    metadata = (;
                        source =
                            :pqs_source_box_route_driver_diatomic_physical_gausslet_h1_j_payload,
                        raw_pair_factor_convention = :raw_numerator,
                    ),
                )
            density_interaction_status = density_interaction.status
            if density_interaction_status !==
               :materialized_pqs_physical_gausslet_pre_final_density_interaction
                status = :blocked_pqs_physical_gausslet_h1_j_payload
                blocker = :physical_gausslet_pre_final_density_interaction_blocked
                push!(missing, :pqs_physical_gausslet_pre_final_density_interaction)
            else
                h1_j_diagnostic =
                    _pqs_source_box_route_driver_physical_gausslet_h1_j_diagnostic(
                        density_interaction,
                        h1_payload.h1.lowest_orbital_coefficients,
                    )
                h1_j_status = :materialized_pqs_physical_gausslet_h1_j_payload
                status = h1_j_status
                blocker = nothing
                push!(available, :pqs_physical_gausslet_h1_j_payload)
                push!(missing, :physical_gausslet_rhf_or_solver_contract)
            end
        catch error
            error isa ArgumentError || error isa DimensionMismatch || rethrow()
            status = :blocked_pqs_physical_gausslet_h1_j_payload
            message = sprint(showerror, error)
            if occursin("axis weights", message)
                blocker = :missing_physical_gausslet_axis_weights
            elseif occursin("raw pair", message)
                blocker = :missing_physical_gausslet_raw_pair_factor_terms
            elseif occursin("support weights", message)
                blocker = :physical_gausslet_support_weight_nonpositive
            else
                blocker = :physical_gausslet_pre_final_density_interaction_blocked
            end
            push!(missing, blocker)
        end
    end

    available_objects = Tuple(unique(available))
    missing_objects = Tuple(unique(missing))
    h1_j_materialized =
        status === :materialized_pqs_physical_gausslet_h1_j_payload
    density_interaction_materialized =
        density_interaction_status ===
        :materialized_pqs_physical_gausslet_pre_final_density_interaction
    summary = (;
        status,
        blocker,
        route_family,
        source_plan_status,
        final_basis_status,
        h1_payload_status,
        density_provenance_status,
        support_weights_status,
        raw_pair_factor_status,
        support_pair_raw_numerator_status,
        density_interaction_status,
        h1_j_status,
        density_gauge =
            isnothing(density_interaction) ? nothing : density_interaction.density_gauge,
        raw_pair_factor_convention =
            isnothing(support_pair_raw_numerator) ? nothing :
            support_pair_raw_numerator.raw_pair_factor_convention,
        support_weight_count =
            isnothing(support_weights) ? nothing : support_weights.support_weight_count,
        support_weights_all_positive =
            isnothing(support_weights) ? nothing :
            support_weights.support_weights_all_positive,
        support_raw_pair_shape =
            isnothing(support_pair_raw_numerator) ? nothing :
            support_pair_raw_numerator.support_pair_raw_numerator_shape,
        support_raw_pair_finite =
            isnothing(support_pair_raw_numerator) ? nothing :
            support_pair_raw_numerator.support_pair_raw_numerator_finite,
        pre_final_pair_matrix_shape =
            isnothing(density_interaction) ? nothing :
            density_interaction.pre_final_pair_matrix_shape,
        pre_final_pair_matrix_finite =
            isnothing(density_interaction) ? nothing :
            density_interaction.pre_final_pair_matrix_finite,
        pre_final_pair_matrix_symmetry_error =
            isnothing(density_interaction) ? nothing :
            density_interaction.pre_final_pair_matrix_symmetry_error,
        self_coulomb =
            isnothing(h1_j_diagnostic) ? nothing : h1_j_diagnostic.self_coulomb,
        h1_j_materialized,
        density_interaction_materialized,
        endpoint_blocker =
            h1_j_materialized ?
            :missing_physical_gausslet_rhf_or_solver_contract :
            blocker,
        rhf_materialized = false,
        exports_materialized = false,
        public_api = false,
        available_objects,
        missing_objects,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_physical_gausslet_h1_j_payload,
        route_kind = recipe.route_kind,
        physical_final_basis = true,
        density_gauge = summary.density_gauge,
        raw_pair_factor_convention = summary.raw_pair_factor_convention,
        h2_221_diagnostic_source_plan_reused = false,
    )
    return _PQSDiatomicPhysicalGaussletH1JPayload(
        status,
        blocker,
        route_family,
        source_plan,
        source_plan_status,
        final_basis,
        final_basis_status,
        h1_payload,
        h1_payload_status,
        density_provenance,
        density_provenance_status,
        support_weights,
        support_weights_status,
        raw_pair_factor_terms,
        raw_pair_factor_status,
        support_pair_raw_numerator,
        support_pair_raw_numerator_status,
        density_interaction,
        density_interaction_status,
        h1_j_diagnostic,
        h1_j_status,
        available_objects,
        missing_objects,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_rhf_input_contract(
    parent,
    route_skeleton,
    recipe,
    source_plan_payload = nothing,
    final_basis_payload = nothing,
    h1_payload = nothing,
    h1_j_payload = nothing,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    route_kind =
        hasproperty(route_skeleton, :route_kind) ?
        route_skeleton.route_kind :
        recipe.route_kind
    fixture_role = get(recipe, :artifact_role, nothing)
    source_plan =
        !isnothing(source_plan_payload) && hasproperty(source_plan_payload, :source_plan) ?
        source_plan_payload.source_plan :
        nothing
    final_basis =
        !isnothing(final_basis_payload) && hasproperty(final_basis_payload, :final_basis) ?
        final_basis_payload.final_basis :
        nothing
    source_plan_status = isnothing(source_plan) ? :not_available : source_plan.status
    final_basis_status = isnothing(final_basis) ? :not_available : final_basis.status
    h1_payload_status = isnothing(h1_payload) ? :not_available : h1_payload.status
    h1_j_payload_status =
        isnothing(h1_j_payload) ? :not_available : h1_j_payload.status
    private_rhf_inputs = get(recipe, :private_rhf_inputs, (;))
    electron_count =
        get(private_rhf_inputs, :private_rhf_electron_count, nothing)
    occupation =
        electron_count === 2 ?
        (; policy = :closed_shell_rhf, nocc = 1, occupancy = 2.0) :
        nothing
    h1_matrix =
        !isnothing(h1_payload) &&
        hasproperty(h1_payload, :final_hamiltonian) &&
        !isnothing(h1_payload.final_hamiltonian) ?
        h1_payload.final_hamiltonian.hamiltonian_matrix :
        nothing
    density_interaction =
        !isnothing(h1_j_payload) && hasproperty(h1_j_payload, :density_interaction) ?
        h1_j_payload.density_interaction :
        nothing
    final_to_pre_final =
        !isnothing(density_interaction) &&
        hasproperty(density_interaction, :final_to_pre_final_coefficients) ?
        density_interaction.final_to_pre_final_coefficients :
        nothing
    pre_final_pair_matrix =
        !isnothing(density_interaction) &&
        hasproperty(density_interaction, :pre_final_pair_matrix) ?
        density_interaction.pre_final_pair_matrix :
        nothing
    available = Symbol[]
    missing = Symbol[]

    if route_family !== :pqs_source_box
        status = :not_applicable_pqs_physical_gausslet_rhf_input_contract
        blocker = nothing
    elseif !(
        route_kind === :bond_aligned_diatomic_physical_gausslet_core_shell_pqs &&
        fixture_role in (
            :physical_gausslet_endpoint_target,
            :fake_pqs_source_backed_wl_reproduction,
        ) ||
        route_kind === :bond_aligned_diatomic_independent_pqs_source_box_core_shell &&
        _pqs_source_box_route_driver_independent_h2_pqs_artifact_role(
            fixture_role,
        )
    )
        status = :blocked_pqs_physical_gausslet_rhf_input_contract
        blocker = :unsupported_physical_gausslet_fixture_role
        push!(missing, :physical_gausslet_fixture_role)
    elseif source_plan_status !==
           :available_pqs_diatomic_physical_gausslet_core_shell_source_plan
        status = :blocked_pqs_physical_gausslet_rhf_input_contract
        blocker = :missing_pqs_diatomic_physical_gausslet_source_plan
        push!(missing, :pqs_diatomic_physical_gausslet_source_plan)
    elseif final_basis_status !== :available_pqs_physical_gausslet_final_basis
        status = :blocked_pqs_physical_gausslet_rhf_input_contract
        blocker = :missing_pqs_physical_gausslet_final_basis
        push!(missing, :pqs_physical_gausslet_final_basis)
    elseif h1_payload_status !== :available_pqs_physical_gausslet_h1_payload ||
           isnothing(h1_matrix)
        status = :blocked_pqs_physical_gausslet_rhf_input_contract
        blocker = :missing_physical_gausslet_h1_payload
        push!(missing, :pqs_physical_gausslet_h1_payload)
    elseif h1_j_payload_status !== :materialized_pqs_physical_gausslet_h1_j_payload ||
           isnothing(density_interaction)
        status = :blocked_pqs_physical_gausslet_rhf_input_contract
        blocker = :missing_physical_gausslet_density_interaction
        push!(missing, :pqs_physical_gausslet_density_interaction)
    elseif isnothing(final_to_pre_final)
        status = :blocked_pqs_physical_gausslet_rhf_input_contract
        blocker = :missing_physical_gausslet_final_to_pre_final_transform
        push!(missing, :physical_gausslet_final_to_pre_final_transform)
    elseif isnothing(pre_final_pair_matrix)
        status = :blocked_pqs_physical_gausslet_rhf_input_contract
        blocker = :missing_physical_gausslet_pre_final_pair_matrix
        push!(missing, :physical_gausslet_pre_final_pair_matrix)
    elseif isnothing(electron_count)
        status = :blocked_pqs_physical_gausslet_rhf_input_contract
        blocker = :missing_physical_gausslet_electron_count
        push!(missing, :physical_gausslet_electron_count)
    elseif electron_count !== 2 || isnothing(occupation)
        status = :blocked_pqs_physical_gausslet_rhf_input_contract
        blocker = :unsupported_physical_gausslet_electron_count
        push!(missing, :physical_gausslet_electron_count)
    else
        h1 = Matrix{Float64}(h1_matrix)
        pair = Matrix{Float64}(pre_final_pair_matrix)
        transform = Matrix{Float64}(final_to_pre_final)
        h1_symmetry_error = norm(h1 - transpose(h1), Inf)
        pair_symmetry_error = norm(pair - transpose(pair), Inf)
        if size(h1) != (final_basis.final_dimension, final_basis.final_dimension) ||
           !all(isfinite, h1) ||
           h1_symmetry_error > 1.0e-8
            status = :blocked_pqs_physical_gausslet_rhf_input_contract
            blocker = :missing_physical_gausslet_h1_payload
            push!(missing, :physical_gausslet_h1_matrix)
        elseif get(density_interaction, :density_gauge, nothing) !==
               :pre_final_localized_positive_weight
            status = :blocked_pqs_physical_gausslet_rhf_input_contract
            blocker = :physical_gausslet_rhf_input_contract_unreviewed
            push!(missing, :physical_gausslet_density_gauge)
        elseif get(density_interaction.metadata, :raw_pair_factor_convention, nothing) !==
               :raw_numerator
            status = :blocked_pqs_physical_gausslet_rhf_input_contract
            blocker = :physical_gausslet_rhf_input_contract_unreviewed
            push!(missing, :physical_gausslet_raw_pair_factor_convention)
        elseif size(pair, 1) != size(pair, 2) ||
               size(transform, 1) != size(pair, 1) ||
               size(transform, 2) != size(h1, 1) ||
               !all(isfinite, pair) ||
               pair_symmetry_error > 1.0e-8
            status = :blocked_pqs_physical_gausslet_rhf_input_contract
            blocker = :missing_physical_gausslet_pre_final_pair_matrix
            push!(missing, :physical_gausslet_pre_final_pair_matrix)
        else
            status = :available_pqs_physical_gausslet_rhf_input_contract
            blocker = nothing
            append!(
                available,
                (
                    :pqs_diatomic_physical_gausslet_source_plan,
                    :pqs_physical_gausslet_final_basis,
                    :pqs_physical_gausslet_h1_payload,
                    :pqs_physical_gausslet_h1_j_payload,
                    :pqs_physical_gausslet_rhf_input_contract,
                ),
            )
        end
    end

    available_objects = Tuple(unique(available))
    missing_objects = Tuple(unique(missing))
    contract_available =
        status === :available_pqs_physical_gausslet_rhf_input_contract
    h1_matrix_finite =
        isnothing(h1_matrix) ? nothing : all(isfinite, Matrix{Float64}(h1_matrix))
    h1_matrix_symmetry_error =
        isnothing(h1_matrix) ? nothing :
        norm(Matrix{Float64}(h1_matrix) - transpose(Matrix{Float64}(h1_matrix)), Inf)
    pre_final_pair_matrix_finite =
        isnothing(pre_final_pair_matrix) ? nothing :
        all(isfinite, Matrix{Float64}(pre_final_pair_matrix))
    pre_final_pair_matrix_symmetry_error =
        isnothing(pre_final_pair_matrix) ? nothing :
        norm(
            Matrix{Float64}(pre_final_pair_matrix) -
            transpose(Matrix{Float64}(pre_final_pair_matrix)),
            Inf,
        )
    summary = (;
        status,
        blocker,
        route_family,
        route_kind,
        fixture_role,
        source_plan_status,
        final_basis_status,
        h1_payload_status,
        h1_j_payload_status,
        input_contract_available = contract_available,
        electron_count,
        occupation_policy = isnothing(occupation) ? nothing : occupation.policy,
        occupation_nocc = isnothing(occupation) ? nothing : occupation.nocc,
        final_dimension =
            isnothing(final_basis) ? nothing : final_basis.final_dimension,
        h1_matrix_available = !isnothing(h1_matrix),
        h1_matrix_finite,
        h1_matrix_symmetry_error,
        density_interaction_available = !isnothing(density_interaction),
        density_gauge =
            isnothing(density_interaction) ? nothing : density_interaction.density_gauge,
        raw_pair_factor_convention =
            isnothing(density_interaction) ? nothing :
            get(density_interaction.metadata, :raw_pair_factor_convention, nothing),
        final_to_pre_final_transform_available = !isnothing(final_to_pre_final),
        pre_final_pair_matrix_available = !isnothing(pre_final_pair_matrix),
        pre_final_pair_matrix_finite,
        pre_final_pair_matrix_symmetry_error,
        private_diagnostic_only = true,
        private_rhf_materialized = false,
        endpoint_blocker =
            contract_available ?
            :missing_physical_gausslet_private_rhf_execution :
            blocker,
        available_objects,
        missing_objects,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_physical_gausslet_rhf_input_contract,
        existing_complete_core_shell_rhf_helpers_adapted = false,
        scf_executed = false,
        private_diagnostic_only = true,
    )
    return _PQSDiatomicPhysicalGaussletRHFInputContractPayload(
        status,
        blocker,
        route_family,
        source_plan,
        source_plan_status,
        final_basis,
        final_basis_status,
        h1_payload,
        h1_payload_status,
        h1_j_payload,
        h1_j_payload_status,
        electron_count,
        occupation,
        available_objects,
        missing_objects,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_rhf_execution_payload(
    input_contract = nothing,
    h1_payload = nothing,
    h1_j_payload = nothing,
)
    input_contract_status =
        isnothing(input_contract) ? :not_available : input_contract.status
    h1_payload_status = isnothing(h1_payload) ? :not_available : h1_payload.status
    h1_j_payload_status =
        isnothing(h1_j_payload) ? :not_available : h1_j_payload.status
    contract_summary =
        isnothing(input_contract) || !hasproperty(input_contract, :summary) ?
        (;) :
        input_contract.summary
    route_family =
        isnothing(input_contract) ? :unknown : input_contract.route_family
    electron_count = get(contract_summary, :electron_count, nothing)
    occupation_nocc = get(contract_summary, :occupation_nocc, nothing)
    scf_payload = nothing
    if input_contract_status !==
       :available_pqs_physical_gausslet_rhf_input_contract
        status = :blocked_pqs_physical_gausslet_private_rhf_execution
        blocker =
            isnothing(input_contract) ?
            :missing_physical_gausslet_rhf_input_contract :
            input_contract.blocker
        executed = false
        materialized = false
        converged = false
        total_energy = nothing
        one_body_energy = nothing
        two_body_energy = nothing
        iteration_count = 0
        density_trace = nothing
        idempotency_residual = nothing
        commutator_residual = nothing
        energy_delta = nothing
        consistency_status =
            :not_evaluated_missing_physical_gausslet_rhf_input_contract
        endpoint_blocker = blocker
    else
        scf_payload = _pqs_multilayer_complete_core_shell_rhf_scf_payload(
            ;
            input_contract,
            h1_payload,
            h1_j_payload,
            mixing_kind = :fock_diis,
            max_iterations = 25,
            density_atol = 1.0e-8,
            energy_atol = 1.0e-10,
            residual_atol = 1.0e-8,
            trace_atol = 1.0e-8,
            idempotency_atol = 1.0e-8,
            metadata = (;
                source =
                    :pqs_source_box_route_driver_diatomic_physical_gausslet_rhf_execution_payload,
                private_diagnostic_only = true,
            ),
        )
        scf_summary = scf_payload.summary
        residual_diagnostics = get(scf_summary, :residual_diagnostics, (;))
        final_one_step = scf_payload.final_one_step_payload
        converged =
            scf_payload.status ===
            :materialized_pqs_multilayer_complete_core_shell_rhf_scf_payload
        materialized = converged
        executed = true
        status =
            converged ?
            :materialized_pqs_physical_gausslet_private_rhf_execution :
            :blocked_pqs_physical_gausslet_private_rhf_execution
        blocker =
            converged ? nothing :
            scf_payload.blocker === :scf_not_converged ?
            :physical_gausslet_private_rhf_not_converged :
            scf_payload.blocker
        total_energy =
            isnothing(final_one_step) ? nothing : final_one_step.total_energy
        one_body_energy =
            isnothing(final_one_step) ? nothing : final_one_step.one_body_energy
        two_body_energy =
            isnothing(final_one_step) ? nothing : final_one_step.two_body_energy
        iteration_count = get(scf_summary, :iteration_count, 0)
        density_trace = get(residual_diagnostics, :density_trace, nothing)
        idempotency_residual =
            get(residual_diagnostics, :closed_shell_idempotency_error, nothing)
        commutator_residual =
            get(residual_diagnostics, :commutator_residual, nothing)
        energy_delta = get(scf_summary, :final_energy_change, nothing)
        consistency_status =
            get(scf_summary, :final_one_step_recomputed, false) ?
            :reviewed_recomputed :
            :not_evaluated_missing_recomputed_final_one_step
        endpoint_blocker =
            converged ? :missing_h2_gausslet_only_reference_comparison : blocker
    end
    summary = (;
        status,
        blocker,
        input_contract_status,
        h1_payload_status,
        h1_j_payload_status,
        input_contract_available =
            input_contract_status ===
            :available_pqs_physical_gausslet_rhf_input_contract,
        executed,
        materialized,
        converged,
        electron_count,
        electron_count_source = :explicit_private_rhf_electron_count,
        occupation_policy = get(contract_summary, :occupation_policy, nothing),
        occupation_nocc,
        total_energy,
        one_body_energy,
        two_body_energy,
        iteration_count,
        density_trace,
        idempotency_residual,
        commutator_residual,
        energy_delta,
        final_density_one_step_consistency_status = consistency_status,
        endpoint_blocker,
        private_diagnostic_only = true,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_physical_gausslet_rhf_execution_payload,
        existing_complete_core_shell_rhf_helpers_adapted = true,
        execution_attempted = summary.executed,
        scf_payload_status = isnothing(scf_payload) ? nothing : scf_payload.status,
        scf_payload_blocker =
            isnothing(scf_payload) ? nothing : scf_payload.blocker,
        private_diagnostic_only = true,
    )
    return _PQSDiatomicPhysicalGaussletRHFExecutionPayload(
        status,
        blocker,
        route_family,
        input_contract,
        input_contract_status,
        h1_payload,
        h1_payload_status,
        h1_j_payload,
        h1_j_payload_status,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_source_box_window(source_box)
    isnothing(source_box) && return nothing
    if source_box isa NamedTuple
        return (source_box.x, source_box.y, source_box.z)
    end
    source_box_tuple = Tuple(source_box)
    length(source_box_tuple) == 3 || return nothing
    return ntuple(axis -> source_box_tuple[axis], 3)
end

function _pqs_source_box_route_driver_window_support_count(window)
    isnothing(window) && return nothing
    return prod(length(window[axis]) for axis in 1:3)
end

function _pqs_source_box_route_driver_diatomic_axis_metrics(
    parent_axis_bundle_object,
)
    axis_pgdg = ntuple(
        axis -> _nested_axis_pgdg(parent_axis_bundle_object, (:x, :y, :z)[axis]),
        3,
    )
    return (;
        x = (;
            overlap = axis_pgdg[1].overlap,
            position = axis_pgdg[1].position,
            x2 = axis_pgdg[1].x2,
            weights = axis_pgdg[1].weights,
            centers = axis_pgdg[1].centers,
            kinetic = axis_pgdg[1].kinetic,
            source = :nested_pgdg_axis,
        ),
        y = (;
            overlap = axis_pgdg[2].overlap,
            position = axis_pgdg[2].position,
            x2 = axis_pgdg[2].x2,
            weights = axis_pgdg[2].weights,
            centers = axis_pgdg[2].centers,
            kinetic = axis_pgdg[2].kinetic,
            source = :nested_pgdg_axis,
        ),
        z = (;
            overlap = axis_pgdg[3].overlap,
            position = axis_pgdg[3].position,
            x2 = axis_pgdg[3].x2,
            weights = axis_pgdg[3].weights,
            centers = axis_pgdg[3].centers,
            kinetic = axis_pgdg[3].kinetic,
            source = :nested_pgdg_axis,
        ),
    )
end

function _pqs_source_box_route_driver_raw_product_box_plan_summary(raw_plan)
    return (;
        object_kind =
            _pqs_source_box_route_driver_descriptor_property(
                raw_plan,
                :object_kind,
            ),
        source_mode_dims =
            _pqs_source_box_route_driver_descriptor_property(
                raw_plan,
                :source_mode_dims,
            ),
        source_mode_count =
            _pqs_source_box_route_driver_descriptor_property(
                raw_plan,
                :source_mode_count,
            ),
        source_mode_ordering =
            _pqs_source_box_route_driver_descriptor_property(
                raw_plan,
                :source_mode_ordering,
            ),
        axis_local_coefficient_shapes = ntuple(
            axis -> size(raw_plan.axis_local_coefficients[axis]),
            3,
        ),
    )
end

function _pqs_source_box_route_driver_raw_pqs_plan_summary(raw_pqs_plan)
    raw_plan =
        CartesianContractedParentMetrics._pqs_raw_product_box_plan_view(
            raw_pqs_plan,
        )
    return (;
        representation =
            _pqs_source_box_route_driver_descriptor_property(
                raw_plan,
                :representation,
            ),
        source_mode_dims =
            _pqs_source_box_route_driver_descriptor_property(
                raw_plan,
                :source_mode_dims,
            ),
        boundary_selected_count =
            raw_plan.boundary_selector.selected_count,
        axis_local_coefficient_shapes = ntuple(
            axis -> size(raw_plan.axis_local_coefficients[axis]),
            3,
        ),
    )
end

function _pqs_source_box_route_driver_product_unit_summary(product_unit)
    return (;
        kind =
            _pqs_source_box_route_driver_descriptor_property(
                product_unit,
                :kind,
            ),
        support_count =
            hasproperty(product_unit, :support_indices) ?
            length(product_unit.support_indices) :
            nothing,
        support_state_count =
            hasproperty(product_unit, :support_states) ?
            length(product_unit.support_states) :
            nothing,
        coefficient_matrix_shape =
            hasproperty(product_unit, :coefficient_matrix) ?
            size(product_unit.coefficient_matrix) :
            nothing,
    )
end

function _pqs_source_box_route_driver_raw_box_pair_inventory_summary(
    pair_inventory,
)
    diagnostics =
        hasproperty(pair_inventory, :diagnostics) ?
        pair_inventory.diagnostics :
        (;)
    return (;
        pair_count =
            hasproperty(pair_inventory, :pair_entries) ?
            length(pair_inventory.pair_entries) :
            _pqs_source_box_route_driver_descriptor_property(
                diagnostics,
                :source_box_algorithmic_pair_count,
            ),
        pair_family_counts =
            hasproperty(pair_inventory, :pair_family_counts) ?
            pair_inventory.pair_family_counts :
            nothing,
        every_pair_uses_source_box_algorithmic_policy =
            _pqs_source_box_route_driver_descriptor_property(
                diagnostics,
                :every_pair_uses_source_box_algorithmic_policy,
                false,
            ),
        source_box_algorithmic_pair_count =
            _pqs_source_box_route_driver_descriptor_property(
                diagnostics,
                :source_box_algorithmic_pair_count,
            ),
    )
end

function _pqs_source_box_route_driver_raw_box_descriptor_summary(descriptor)
    return (;
        object_kind =
            _pqs_source_box_route_driver_descriptor_property(
                descriptor,
                :object_kind,
            ),
        roles =
            _pqs_source_box_route_driver_descriptor_property(
                descriptor,
                :roles,
            ),
        retained_dimension =
            _pqs_source_box_route_driver_descriptor_property(
                descriptor,
                :retained_dimension,
            ),
        expected_pair_count =
            _pqs_source_box_route_driver_descriptor_property(
                descriptor,
                :expected_pair_count,
            ),
    )
end

function _pqs_source_box_route_driver_diatomic_complete_core_shell_support_window_payload(
    parent,
    route_skeleton,
    recipe,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    system_classification =
        hasproperty(parent, :system_classification) ?
        parent.system_classification :
        nothing
    bond_axis =
        hasproperty(parent, :bond_axis) ? parent.bond_axis : nothing
    parent_axis_bundle_object_available =
        hasproperty(parent, :parent_axis_bundle_object_available) ?
        parent.parent_axis_bundle_object_available :
        (
            hasproperty(parent, :parent_axis_bundle_object) &&
            !isnothing(parent.parent_axis_bundle_object)
        )
    parent_dims =
        hasproperty(route_skeleton, :parent_axis_counts) ?
        _pqs_source_box_route_driver_axis_counts_tuple(
            route_skeleton.parent_axis_counts,
        ) :
        hasproperty(parent, :axis_counts) ?
        _pqs_source_box_route_driver_axis_counts_tuple(parent.axis_counts) :
        nothing
    source_boxes =
        hasproperty(route_skeleton, :source_boxes) ?
        route_skeleton.source_boxes :
        (;)
    source_dimensions =
        hasproperty(route_skeleton, :source_dimensions) ?
        route_skeleton.source_dimensions :
        (;)
    retained_units =
        hasproperty(route_skeleton, :retained_units) ?
        Tuple(route_skeleton.retained_units) :
        ()
    retained_order = Tuple(
        _pqs_source_box_route_driver_descriptor_property(unit, :unit_key)
        for unit in retained_units
    )
    expected_window_keys = (:pqs_left, :product, :pqs_right)
    source_box_windows = (;
        pqs_left =
            hasproperty(source_boxes, :pqs_left) ?
            _pqs_source_box_route_driver_source_box_window(source_boxes.pqs_left) :
            nothing,
        product =
            hasproperty(source_boxes, :product) ?
            _pqs_source_box_route_driver_source_box_window(source_boxes.product) :
            nothing,
        pqs_right =
            hasproperty(source_boxes, :pqs_right) ?
            _pqs_source_box_route_driver_source_box_window(source_boxes.pqs_right) :
            nothing,
    )
    source_mode_dims = (;
        pqs_left =
            hasproperty(source_dimensions, :pqs_left) ?
            Tuple(source_dimensions.pqs_left) :
            nothing,
        product =
            hasproperty(source_dimensions, :product) ?
            Tuple(source_dimensions.product) :
            nothing,
        pqs_right =
            hasproperty(source_dimensions, :pqs_right) ?
            Tuple(source_dimensions.pqs_right) :
            nothing,
    )
    support_counts = (;
        pqs_left =
            _pqs_source_box_route_driver_window_support_count(
                source_box_windows.pqs_left,
            ),
        product =
            _pqs_source_box_route_driver_window_support_count(
                source_box_windows.product,
            ),
        pqs_right =
            _pqs_source_box_route_driver_window_support_count(
                source_box_windows.pqs_right,
            ),
    )
    candidate_core_then_shell_support_order = (:product, :pqs_left, :pqs_right)
    retained_to_support_order_permutation_required =
        retained_order != candidate_core_then_shell_support_order

    missing_support_inputs = Symbol[]
    isnothing(parent_dims) && push!(missing_support_inputs, :parent_dims)
    for key in expected_window_keys
        isnothing(getproperty(source_box_windows, key)) &&
            push!(missing_support_inputs, Symbol(key, :_source_box_window))
        isnothing(getproperty(source_mode_dims, key)) &&
            push!(missing_support_inputs, Symbol(key, :_source_mode_dims))
        isnothing(getproperty(support_counts, key)) &&
            push!(missing_support_inputs, Symbol(key, :_support_count))
    end
    isempty(retained_order) && push!(missing_support_inputs, :retained_order)

    materializer_missing_objects = (
        :raw_product_box_plan_objects,
        :pqs_axis_local_coefficients,
        :diatomic_complete_core_shell_source_plan_materializer,
    )
    if route_family !== :pqs_source_box
        status =
            :not_applicable_diatomic_complete_core_shell_support_windows_non_pqs_route
        blocker = nothing
        missing_objects = ()
    elseif system_classification !== :bond_aligned_diatomic
        status =
            :not_applicable_diatomic_complete_core_shell_support_windows_non_diatomic
        blocker = nothing
        missing_objects = ()
    elseif isempty(missing_support_inputs)
        status = :available_diatomic_complete_core_shell_support_windows
        blocker = nothing
        missing_objects = materializer_missing_objects
    else
        status = :blocked_diatomic_complete_core_shell_support_windows
        blocker = first(missing_support_inputs)
        missing_objects = Tuple((missing_support_inputs..., materializer_missing_objects...))
    end

    summary = (;
        status,
        blocker,
        route_family,
        system_classification,
        bond_axis,
        parent_dims,
        parent_axis_bundle_object_available,
        source_box_keys = expected_window_keys,
        retained_order,
        candidate_core_then_shell_support_order,
        retained_to_support_order_permutation_required,
        support_counts,
        missing_objects,
        support_states_materialized = false,
        raw_product_box_plans_materialized = false,
        source_coefficients_materialized = false,
        source_plan_materialized = false,
        final_basis_materialized = false,
        h1_materialized = false,
        h1_j_materialized = false,
        ham_payload_materialized = false,
        route_driver_public_surface = false,
        exports_materialized = false,
        artifacts_materialized = false,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_complete_core_shell_support_window_payload,
        route_kind = recipe.route_kind,
        route_family,
        source_box_first = true,
        metadata_only_support_windows = true,
        support_states_materialized = false,
        raw_product_box_plans_materialized = false,
        source_coefficients_materialized = false,
        raw_product_box_probe_authority = false,
        shell_support_row_contraction_authority = false,
        retained_diagnostic_weights_are_ida_weights = false,
    )
    return _PQSDiatomicCompleteCoreShellSupportWindowPayload(
        status,
        blocker,
        route_family,
        system_classification,
        bond_axis,
        parent_dims,
        parent_axis_bundle_object_available,
        source_box_windows,
        source_mode_dims,
        retained_order,
        candidate_core_then_shell_support_order,
        retained_to_support_order_permutation_required,
        support_counts,
        missing_objects,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_diatomic_raw_box_route_payload(
    parent,
    route_skeleton,
    recipe,
    support_window_payload = nothing,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    system_classification =
        hasproperty(parent, :system_classification) ?
        parent.system_classification :
        nothing
    parent_axis_bundle_object =
        hasproperty(parent, :parent_axis_bundle_object) ?
        parent.parent_axis_bundle_object :
        nothing
    parent_axis_bundle_object_available =
        hasproperty(parent, :parent_axis_bundle_object_available) ?
        parent.parent_axis_bundle_object_available :
        !isnothing(parent_axis_bundle_object)
    support_window_payload_value =
        isnothing(support_window_payload) ?
        _pqs_source_box_route_driver_diatomic_complete_core_shell_support_window_payload(
            parent,
            route_skeleton,
            recipe,
        ) :
        support_window_payload
    support_window_payload_status = support_window_payload_value.status

    producer = nothing
    producer_status = :not_available
    descriptor_summary = nothing
    raw_product_box_plan_summary = nothing
    raw_pqs_plan_summary = nothing
    product_unit_summary = nothing
    pair_inventory_summary = nothing

    available = Symbol[]
    missing = Symbol[]
    if route_family !== :pqs_source_box
        status = :not_applicable_diatomic_raw_box_route_payload_non_pqs_route
        blocker = nothing
    elseif system_classification !== :bond_aligned_diatomic
        status = :not_applicable_diatomic_raw_box_route_payload_non_diatomic
        blocker = nothing
    elseif !parent_axis_bundle_object_available
        status = :blocked_diatomic_raw_box_route_payload
        blocker = :missing_parent_axis_bundle_object
        support_window_payload_status ===
            :available_diatomic_complete_core_shell_support_windows &&
            push!(available, :diatomic_complete_core_shell_support_windows)
        append!(
            missing,
            (
                :parent_axis_bundle_object,
                :raw_product_box_plan_objects,
                :pqs_axis_local_coefficients,
                :product_doside_unit,
                :pair_inventory,
            ),
        )
    elseif support_window_payload_status !==
           :available_diatomic_complete_core_shell_support_windows
        status = :blocked_diatomic_raw_box_route_payload
        blocker = :missing_diatomic_complete_core_shell_support_windows
        push!(available, :parent_axis_bundle_object)
        push!(missing, :diatomic_complete_core_shell_support_windows)
    else
        metrics =
            _pqs_source_box_route_driver_diatomic_axis_metrics(
                parent_axis_bundle_object,
            )
        producer =
            CartesianContractedParentMetrics._pqs_pqs_product_raw_box_route_producer(
                parent_axis_bundle_object,
                support_window_payload_value.source_box_windows.pqs_left,
                support_window_payload_value.source_box_windows.pqs_right,
                support_window_payload_value.source_box_windows.product,
                metrics;
                source_mode_dims =
                    support_window_payload_value.source_mode_dims.pqs_left,
                route_name = :diatomic_raw_box_route_payload,
                parent_dims = support_window_payload_value.parent_dims,
                bond_axis =
                    hasproperty(parent, :bond_axis) ? parent.bond_axis : nothing,
                metadata = (;
                    source = :diatomic_raw_box_route_payload,
                    route_kind = recipe.route_kind,
                ),
                provenance = (source = :diatomic_raw_box_route_payload,),
            )
        producer_status =
            _pqs_source_box_route_driver_descriptor_property(
                producer,
                :status,
                :unknown,
            )
        descriptor_summary =
            _pqs_source_box_route_driver_raw_box_descriptor_summary(
                producer.descriptor,
            )
        raw_product_box_plan_summary = (;
            pqs_left =
                _pqs_source_box_route_driver_raw_product_box_plan_summary(
                    producer.raw_product_box_plans.pqs_left,
                ),
            pqs_right =
                _pqs_source_box_route_driver_raw_product_box_plan_summary(
                    producer.raw_product_box_plans.pqs_right,
                ),
        )
        raw_pqs_plan_summary = (;
            pqs_left =
                _pqs_source_box_route_driver_raw_pqs_plan_summary(
                    producer.raw_pqs_plans.pqs_left,
                ),
            pqs_right =
                _pqs_source_box_route_driver_raw_pqs_plan_summary(
                    producer.raw_pqs_plans.pqs_right,
                ),
        )
        product_unit_summary =
            _pqs_source_box_route_driver_product_unit_summary(
                producer.product_unit,
            )
        pair_inventory_summary =
            _pqs_source_box_route_driver_raw_box_pair_inventory_summary(
                producer.all_pairs_inventory,
            )
        status = :available_diatomic_raw_box_route_payload
        blocker = nothing
        append!(
            available,
            (
                :parent_axis_bundle_object,
                :diatomic_complete_core_shell_support_windows,
                :raw_box_route_producer,
                :raw_product_box_plan_objects,
                :pqs_axis_local_coefficients,
                :product_doside_unit,
                :pair_inventory,
            ),
        )
        push!(missing, :diatomic_complete_core_shell_source_plan_materializer)
    end

    available_objects = Tuple(available)
    missing_objects = Tuple(missing)
    summary = (;
        status,
        blocker,
        producer_status,
        support_window_payload_status,
        descriptor_summary,
        raw_product_box_plan_summary,
        raw_pqs_plan_summary,
        product_unit_summary,
        pair_inventory_summary,
        available_objects,
        missing_objects,
        private_candidate_only = true,
        raw_product_box_probe_authority = false,
        source_plan_materialized = false,
        final_basis_materialized = false,
        h1_materialized = false,
        h1_j_materialized = false,
        ham_payload_materialized = false,
        route_driver_public_surface = false,
        exports_materialized = false,
        artifacts_materialized = false,
    )
    metadata = (;
        source = :pqs_source_box_route_driver_diatomic_raw_box_route_payload,
        route_kind = recipe.route_kind,
        route_family,
        support_window_payload_status,
        private_candidate_only = true,
        raw_product_box_probe_authority = false,
        source_box_first = true,
        shell_support_row_contraction_authority = false,
        retained_diagnostic_weights_are_ida_weights = false,
    )
    return _PQSDiatomicRawBoxRoutePayload(
        status,
        blocker,
        producer,
        producer_status,
        descriptor_summary,
        raw_product_box_plan_summary,
        raw_pqs_plan_summary,
        product_unit_summary,
        pair_inventory_summary,
        support_window_payload_status,
        available_objects,
        missing_objects,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_diatomic_complete_core_shell_source_realization_payload(
    parent,
    route_skeleton,
    recipe,
    support_window_payload = nothing,
    raw_box_route_payload = nothing,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    system_classification =
        hasproperty(parent, :system_classification) ?
        parent.system_classification :
        nothing
    bond_axis =
        hasproperty(parent, :bond_axis) ? parent.bond_axis : nothing
    parent_axis_bundle_object_available =
        hasproperty(parent, :parent_axis_bundle_object_available) ?
        parent.parent_axis_bundle_object_available :
        (
            hasproperty(parent, :parent_axis_bundle_object) &&
            !isnothing(parent.parent_axis_bundle_object)
        )
    support_window_payload_value =
        isnothing(support_window_payload) ?
        _pqs_source_box_route_driver_diatomic_complete_core_shell_support_window_payload(
            parent,
            route_skeleton,
            recipe,
        ) :
        support_window_payload
    support_window_payload_status = support_window_payload_value.status
    raw_box_route_payload_status =
        isnothing(raw_box_route_payload) ?
        :not_available :
        raw_box_route_payload.status

    core_unit_key = :product
    shell_unit_keys = (:pqs_left, :pqs_right)
    retained_order = support_window_payload_value.retained_order
    support_order = support_window_payload_value.candidate_core_then_shell_support_order
    retained_to_support_order_permutation_required =
        support_window_payload_value.retained_to_support_order_permutation_required
    support_counts = support_window_payload_value.support_counts
    core_support_count = support_counts.product
    shell_support_counts = (;
        pqs_left = support_counts.pqs_left,
        pqs_right = support_counts.pqs_right,
    )
    shell_support_count =
        isnothing(shell_support_counts.pqs_left) ||
        isnothing(shell_support_counts.pqs_right) ?
        nothing :
        shell_support_counts.pqs_left + shell_support_counts.pqs_right
    shell_retained_counts = (pqs_left = nothing, pqs_right = nothing)
    shell_retained_count = nothing
    route_retained_ranges = nothing
    source_plan_precleanup_ranges = nothing
    precleanup_retained_dimension = nothing
    shell_final_coefficients_shape = nothing

    available = Symbol[]
    missing = Symbol[]
    if route_family !== :pqs_source_box
        status =
            :not_applicable_diatomic_complete_core_shell_source_realization_non_pqs_route
        blocker = nothing
    elseif system_classification !== :bond_aligned_diatomic
        status =
            :not_applicable_diatomic_complete_core_shell_source_realization_non_diatomic
        blocker = nothing
    elseif support_window_payload_status !==
           :available_diatomic_complete_core_shell_support_windows
        status = :blocked_diatomic_complete_core_shell_source_realization
        blocker = :missing_diatomic_complete_core_shell_support_windows
        push!(missing, :diatomic_complete_core_shell_support_windows)
    elseif raw_box_route_payload_status !== :available_diatomic_raw_box_route_payload
        status = :blocked_diatomic_complete_core_shell_source_realization
        blocker =
            !isnothing(raw_box_route_payload) &&
            raw_box_route_payload.blocker === :missing_parent_axis_bundle_object ?
            :missing_parent_axis_bundle_object :
            :missing_diatomic_raw_box_route_payload
        push!(available, :diatomic_complete_core_shell_support_windows)
        if !isnothing(raw_box_route_payload)
            append!(missing, raw_box_route_payload.missing_objects)
        else
            push!(missing, :diatomic_raw_box_route_payload)
        end
    else
        status = :available_diatomic_complete_core_shell_source_realization
        blocker = nothing
        producer = raw_box_route_payload.producer
        route_retained_ranges = producer.descriptor.expected_ranges
        left_retained_count =
            raw_box_route_payload.raw_pqs_plan_summary.pqs_left.boundary_selected_count
        right_retained_count =
            raw_box_route_payload.raw_pqs_plan_summary.pqs_right.boundary_selected_count
        shell_retained_counts = (;
            pqs_left = left_retained_count,
            pqs_right = right_retained_count,
        )
        shell_retained_count = left_retained_count + right_retained_count
        precleanup_retained_dimension = core_support_count + shell_retained_count
        source_plan_precleanup_ranges = (;
            product = 1:core_support_count,
            pqs_left = (core_support_count + 1):(core_support_count + left_retained_count),
            pqs_right =
                (core_support_count + left_retained_count + 1):
                precleanup_retained_dimension,
        )
        shell_final_coefficients_shape =
            (shell_support_count, shell_retained_count)
        append!(
            available,
            (
                :parent_axis_bundle_object,
                :diatomic_complete_core_shell_support_windows,
                :diatomic_raw_box_route_payload,
                :raw_box_route_producer,
                :raw_product_box_plan_objects,
                :pqs_axis_local_coefficients,
                :product_doside_unit,
                :pair_inventory,
                :diatomic_complete_core_shell_source_realization,
            ),
        )
        push!(missing, :pqs_multilayer_shell_source_plan_adapter_contract)
    end

    if status !== :available_diatomic_complete_core_shell_source_realization
        parent_axis_bundle_object_available ||
            push!(missing, :parent_axis_bundle_object)
        push!(missing, :diatomic_raw_box_route_payload)
    end
    available_objects = Tuple(unique(available))
    missing_objects = Tuple(unique(missing))
    shell_coefficient_block_structure = :block_diagonal_left_right_pqs
    bundles_role = :parent_axis_bundle
    object_kind_claim = :not_pqs_multilayer_shell_source_plan
    summary = (;
        status,
        blocker,
        route_family,
        system_classification,
        bond_axis,
        parent_axis_bundle_object_available,
        support_window_payload_status,
        raw_box_route_payload_status,
        core_unit_key,
        shell_unit_keys,
        retained_order,
        support_order,
        retained_to_support_order_permutation_required,
        route_retained_ranges,
        source_plan_precleanup_ranges,
        core_support_count,
        shell_support_counts,
        shell_support_count,
        shell_retained_counts,
        shell_retained_count,
        precleanup_retained_dimension,
        shell_final_coefficients_shape,
        shell_coefficient_block_structure,
        bundles_role,
        object_kind_claim,
        available_objects,
        missing_objects,
        source_plan_materialized = false,
        returns_pqs_multilayer_shell_source_plan = false,
        final_basis_materialized = false,
        h1_materialized = false,
        h1_j_materialized = false,
        ham_payload_materialized = false,
        route_driver_public_surface = false,
        exports_materialized = false,
        artifacts_materialized = false,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_complete_core_shell_source_realization_payload,
        route_kind = recipe.route_kind,
        route_family,
        support_window_payload_status,
        raw_box_route_payload_status,
        source_box_first = true,
        private_internal_only = true,
        bundles_role,
        object_kind_claim,
        source_plan_materialized = false,
        returns_pqs_multilayer_shell_source_plan = false,
        support_matrices_materialized = false,
        shell_final_coefficients_materialized = false,
        final_basis_materialized = false,
        h1_materialized = false,
        h1_j_materialized = false,
        ham_payload_materialized = false,
        shell_support_row_contraction_authority = false,
        retained_diagnostic_weights_are_ida_weights = false,
    )
    return _PQSDiatomicCompleteCoreShellSourceRealizationPayload(
        status,
        blocker,
        route_family,
        system_classification,
        bond_axis,
        parent_axis_bundle_object_available,
        support_window_payload_status,
        raw_box_route_payload_status,
        core_unit_key,
        shell_unit_keys,
        retained_order,
        support_order,
        retained_to_support_order_permutation_required,
        route_retained_ranges,
        source_plan_precleanup_ranges,
        core_support_count,
        shell_support_counts,
        shell_support_count,
        shell_retained_counts,
        shell_retained_count,
        precleanup_retained_dimension,
        shell_final_coefficients_shape,
        shell_coefficient_block_structure,
        bundles_role,
        object_kind_claim,
        available_objects,
        missing_objects,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_raw_plan_support_rows(
    raw_pqs_plan,
    parent_dims::NTuple{3,Int},
)
    plan =
        CartesianContractedParentMetrics._pqs_raw_product_box_plan_view(
            raw_pqs_plan,
        )
    intervals = plan.axis_intervals
    support_states = NTuple{3,Int}[]
    for ix in intervals[1], iy in intervals[2], iz in intervals[3]
        push!(support_states, (ix, iy, iz))
    end
    support_indices = [
        _cartesian_flat_index(state[1], state[2], state[3], parent_dims) for
        state in support_states
    ]
    return (; support_indices, support_states)
end

function _pqs_source_box_route_driver_raw_plan_support_coefficients(
    raw_pqs_plan,
    support_indices::AbstractVector{Int},
    parent_dims::NTuple{3,Int},
)
    parent_coefficients =
        CartesianContractedParentMetrics._pqs_parent_coefficient_matrix_from_raw_plan(
            raw_pqs_plan,
            parent_dims,
        )
    return Matrix{Float64}(parent_coefficients[support_indices, :])
end

function _pqs_source_box_route_driver_diatomic_complete_core_shell_source_plan(
    parent,
    recipe,
    support_window_payload,
    raw_box_route_payload,
    source_realization_payload,
)
    parent_axis_bundle_object =
        hasproperty(parent, :parent_axis_bundle_object) ?
        parent.parent_axis_bundle_object :
        nothing
    isnothing(parent_axis_bundle_object) &&
        throw(ArgumentError("diatomic complete core/shell source plan requires parent axis bundle object"))
    raw_box_route_payload.status === :available_diatomic_raw_box_route_payload ||
        throw(ArgumentError("diatomic complete core/shell source plan requires available raw-box route payload"))
    source_realization_payload.status ===
        :available_diatomic_complete_core_shell_source_realization ||
        throw(ArgumentError("diatomic complete core/shell source plan requires available source realization payload"))
    parent_dims = support_window_payload.parent_dims
    parent_dims isa NTuple{3,Int} ||
        throw(ArgumentError("diatomic complete core/shell source plan requires parent dimensions"))

    producer = raw_box_route_payload.producer
    core_support_indices = Int[index for index in producer.product_unit.support_indices]
    core_support_states =
        NTuple{3,Int}[state for state in producer.product_unit.support_states]
    left_support =
        _pqs_source_box_route_driver_raw_plan_support_rows(
            producer.raw_pqs_plans.pqs_left,
            parent_dims,
        )
    right_support =
        _pqs_source_box_route_driver_raw_plan_support_rows(
            producer.raw_pqs_plans.pqs_right,
            parent_dims,
        )
    shell_support_indices = vcat(
        Int[index for index in left_support.support_indices],
        Int[index for index in right_support.support_indices],
    )
    shell_support_states = vcat(
        NTuple{3,Int}[state for state in left_support.support_states],
        NTuple{3,Int}[state for state in right_support.support_states],
    )
    left_coefficients =
        _pqs_source_box_route_driver_raw_plan_support_coefficients(
            producer.raw_pqs_plans.pqs_left,
            left_support.support_indices,
            parent_dims,
        )
    right_coefficients =
        _pqs_source_box_route_driver_raw_plan_support_coefficients(
            producer.raw_pqs_plans.pqs_right,
            right_support.support_indices,
            parent_dims,
        )
    shell_support_count = length(shell_support_indices)
    left_retained_count = size(left_coefficients, 2)
    right_retained_count = size(right_coefficients, 2)
    shell_retained_count = left_retained_count + right_retained_count
    shell_final_coefficients = zeros(
        Float64,
        shell_support_count,
        shell_retained_count,
    )
    left_rows = 1:length(left_support.support_indices)
    right_rows = (
        last(left_rows) + 1
    ):(length(left_support.support_indices) + length(right_support.support_indices))
    left_columns = 1:left_retained_count
    right_columns = (left_retained_count + 1):shell_retained_count
    shell_final_coefficients[left_rows, left_columns] .= left_coefficients
    shell_final_coefficients[right_rows, right_columns] .= right_coefficients

    core_shell_duplicate_count =
        length(intersect(core_support_indices, shell_support_indices))
    shell_duplicate_count =
        length(shell_support_indices) - length(unique(shell_support_indices))
    blocker =
        shell_duplicate_count != 0 ? :duplicate_diatomic_shell_support :
        core_shell_duplicate_count != 0 ? :diatomic_core_shell_support_overlap :
        nothing
    status =
        isnothing(blocker) ?
        :available_pqs_diatomic_complete_core_shell_source_plan :
        :blocked_pqs_diatomic_complete_core_shell_source_plan
    metrics =
        _pqs_source_box_route_driver_diatomic_axis_metrics(
            parent_axis_bundle_object,
        )
    retained_pre_final_map = (;
        support_order = source_realization_payload.support_order,
        route_retained_order = source_realization_payload.retained_order,
        route_retained_ranges = source_realization_payload.route_retained_ranges,
        precleanup_ranges = source_realization_payload.source_plan_precleanup_ranges,
        retained_to_support_order_permutation_required =
            source_realization_payload.retained_to_support_order_permutation_required,
    )
    source_unit_summaries = (;
        pqs_left = raw_box_route_payload.raw_pqs_plan_summary.pqs_left,
        pqs_right = raw_box_route_payload.raw_pqs_plan_summary.pqs_right,
        product = raw_box_route_payload.product_unit_summary,
    )
    convention_labels = (;
        source_plan_family = :diatomic_complete_core_shell_source_plan,
        old_source_plan_object_kind = false,
        support_row_order = :core_product_then_shell_left_right_pqs,
        shell_coefficient_block_structure = :block_diagonal_left_right_pqs,
        source_box_first = true,
        shell_support_row_contraction_authority = false,
        retained_diagnostic_weights_are_ida_weights = false,
        final_basis_materialized = false,
        h1_materialized = false,
        h1_j_materialized = false,
        ham_payload_materialized = false,
        route_driver_public_surface = false,
        exports_materialized = false,
        artifacts_materialized = false,
    )
    summary = (;
        status,
        blocker,
        object_kind = :pqs_diatomic_complete_core_shell_source_plan,
        old_source_plan_object_kind = false,
        core_support_count = length(core_support_indices),
        shell_support_count,
        shell_retained_count,
        precleanup_retained_dimension =
            length(core_support_indices) + shell_retained_count,
        shell_final_coefficients_shape = size(shell_final_coefficients),
        shell_coefficient_block_structure = :block_diagonal_left_right_pqs,
        support_order = source_realization_payload.support_order,
        route_retained_order = source_realization_payload.retained_order,
        core_shell_duplicate_count,
        shell_duplicate_count,
        source_plan_materialized = true,
        final_basis_materialized = false,
        h1_materialized = false,
        h1_j_materialized = false,
        ham_payload_materialized = false,
        route_driver_public_surface = false,
        exports_materialized = false,
        artifacts_materialized = false,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_complete_core_shell_source_plan,
        route_kind = recipe.route_kind,
        route_family = source_realization_payload.route_family,
        source_realization_status = source_realization_payload.status,
        raw_box_route_payload_status = raw_box_route_payload.status,
        support_window_payload_status = support_window_payload.status,
        source_box_first = true,
        old_source_plan_object_kind = false,
        returns_pqs_multilayer_shell_source_plan = false,
        final_basis_materialized = false,
        h1_materialized = false,
        h1_j_materialized = false,
        ham_payload_materialized = false,
        shell_support_row_contraction_authority = false,
        retained_diagnostic_weights_are_ida_weights = false,
    )
    return _PQSDiatomicCompleteCoreShellSourcePlan(
        :pqs_diatomic_complete_core_shell_source_plan,
        status,
        blocker,
        parent_axis_bundle_object,
        metrics,
        source_realization_payload.core_unit_key,
        source_realization_payload.shell_unit_keys,
        core_support_indices,
        core_support_states,
        shell_support_indices,
        shell_support_states,
        shell_final_coefficients,
        source_realization_payload.support_order,
        source_realization_payload.retained_order,
        retained_pre_final_map,
        source_unit_summaries,
        convention_labels,
        summary,
        metadata,
    )
end

struct _PQSDiatomicCompleteCoreShellSourcePlanPayload
    status::Symbol
    blocker
    route_family::Symbol
    system_classification
    bond_axis
    parent_axis_bundle_object_available::Bool
    route_skeleton_summary
    source_box_summary
    retained_unit_summary
    pair_inventory_summary
    center_summary
    coulomb_expansion_summary
    support_window_payload
    source_realization_payload
    source_plan
    source_plan_status::Symbol
    available_objects::Tuple
    missing_objects::Tuple
    summary
    metadata
end

function _pqs_source_box_route_driver_diatomic_coulomb_expansion_summary(
    route_family,
)
    route_family === :pqs_source_box || return (;
        status = :not_applicable,
        coefficient_count = 0,
        coefficients_positive = false,
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    return (;
        status = :available_coulomb_gaussian_expansion,
        coefficient_count = length(expansion.coefficients),
        coefficients_positive = all(>(0.0), expansion.coefficients),
    )
end

function _pqs_source_box_route_driver_diatomic_complete_core_shell_source_plan_payload(
    parent,
    route_skeleton,
    recipe,
    support_window_payload = nothing,
    raw_box_route_payload = nothing,
    source_realization_payload = nothing,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    system_classification =
        hasproperty(parent, :system_classification) ?
        parent.system_classification :
        nothing
    bond_axis =
        hasproperty(parent, :bond_axis) ? parent.bond_axis : nothing
    parent_axis_bundle_object_available =
        hasproperty(parent, :parent_axis_bundle_object_available) ?
        parent.parent_axis_bundle_object_available :
        (
            hasproperty(parent, :parent_axis_bundle_object) &&
            !isnothing(parent.parent_axis_bundle_object)
        )
    center_summary =
        _pqs_source_box_route_driver_diatomic_center_summary(parent)
    source_box_summary =
        _pqs_source_box_route_driver_diatomic_source_box_summary(route_skeleton)
    retained_unit_summary =
        _pqs_source_box_route_driver_diatomic_retained_unit_summary(route_skeleton)
    pair_inventory_summary =
        _pqs_source_box_route_driver_diatomic_pair_inventory_summary(route_skeleton)
    route_skeleton_summary =
        _pqs_source_box_route_driver_diatomic_route_skeleton_summary(
            route_skeleton,
            retained_unit_summary,
        )
    coulomb_expansion_summary =
        _pqs_source_box_route_driver_diatomic_coulomb_expansion_summary(
            route_family,
        )
    support_window_payload_value =
        isnothing(support_window_payload) ?
        _pqs_source_box_route_driver_diatomic_complete_core_shell_support_window_payload(
            parent,
            route_skeleton,
            recipe,
        ) :
        support_window_payload
    support_window_payload_status = support_window_payload_value.status
    raw_box_route_payload_status =
        isnothing(raw_box_route_payload) ?
        :not_available :
        raw_box_route_payload.status
    source_realization_payload_value =
        isnothing(source_realization_payload) ?
        _pqs_source_box_route_driver_diatomic_complete_core_shell_source_realization_payload(
            parent,
            route_skeleton,
            recipe,
            support_window_payload_value,
            raw_box_route_payload,
        ) :
        source_realization_payload
    source_realization_payload_status = source_realization_payload_value.status

    source_plan =
        source_realization_payload_status ===
        :available_diatomic_complete_core_shell_source_realization &&
        raw_box_route_payload_status === :available_diatomic_raw_box_route_payload ?
        _pqs_source_box_route_driver_diatomic_complete_core_shell_source_plan(
            parent,
            recipe,
            support_window_payload_value,
            raw_box_route_payload,
            source_realization_payload_value,
        ) :
        nothing
    source_plan_status =
        isnothing(source_plan) ?
        :not_materialized_diatomic_complete_core_shell_source_plan :
        source_plan.status
    if route_family !== :pqs_source_box
        status =
            :not_applicable_diatomic_complete_core_shell_source_plan_non_pqs_route
        blocker = nothing
        available_objects = ()
        missing_objects = ()
        source_plan_status = :not_applicable
    elseif system_classification !== :bond_aligned_diatomic
        status =
            :not_applicable_diatomic_complete_core_shell_source_plan_non_diatomic
        blocker = nothing
        available_objects = ()
        missing_objects = ()
        source_plan_status = :not_applicable
    else
        status = :blocked_diatomic_complete_core_shell_source_plan
        blocker =
            source_plan_status ===
            :available_pqs_diatomic_complete_core_shell_source_plan ?
            :missing_diatomic_complete_core_shell_final_basis_consumer :
            source_realization_payload_status ===
            :available_diatomic_complete_core_shell_source_realization ?
            :missing_pqs_diatomic_complete_core_shell_source_plan :
            parent_axis_bundle_object_available ?
            :missing_diatomic_complete_core_shell_source_realization_payload :
            :missing_parent_axis_bundle_object
        available = Symbol[
            :route_skeleton,
            :diatomic_center_metadata,
            :source_boxes,
            :retained_units,
            :pair_inventory,
            :coulomb_expansion,
        ]
        support_window_payload_value.status ===
            :available_diatomic_complete_core_shell_support_windows &&
            push!(available, :diatomic_complete_core_shell_support_windows)
        raw_box_route_payload_status === :available_diatomic_raw_box_route_payload &&
            push!(available, :diatomic_raw_box_route_payload)
        source_realization_payload_status ===
            :available_diatomic_complete_core_shell_source_realization &&
            push!(available, :diatomic_complete_core_shell_source_realization)
        source_plan_status ===
            :available_pqs_diatomic_complete_core_shell_source_plan &&
            push!(available, :pqs_diatomic_complete_core_shell_source_plan)
        parent_axis_bundle_object_available &&
            push!(available, :parent_axis_bundle_object)
        missing =
            source_plan_status ===
            :available_pqs_diatomic_complete_core_shell_source_plan ?
            Symbol[:diatomic_complete_core_shell_final_basis_consumer] :
            source_realization_payload_status ===
            :available_diatomic_complete_core_shell_source_realization ?
            Symbol[:pqs_diatomic_complete_core_shell_source_plan] :
            Symbol[:diatomic_complete_core_shell_source_realization_payload]
        parent_axis_bundle_object_available ||
            push!(missing, :parent_axis_bundle_object)
        available_objects = Tuple(available)
        missing_objects = Tuple(missing)
    end

    summary = (;
        status,
        blocker,
        route_family,
        system_classification,
        bond_axis,
        parent_axis_bundle_object_available,
        center_count = center_summary.center_count,
        source_box_count = source_box_summary.source_box_count,
        retained_unit_count = retained_unit_summary.retained_unit_count,
        pair_count = pair_inventory_summary.pair_count,
        coulomb_expansion_status = coulomb_expansion_summary.status,
        support_window_payload_status,
        raw_box_route_payload_status,
        source_realization_payload_status,
        source_plan_status,
        missing_objects,
        source_plan_materialized = !isnothing(source_plan),
        final_basis_materialized = false,
        h1_materialized = false,
        h1_j_materialized = false,
        ham_payload_materialized = false,
        route_driver_public_surface = false,
        exports_materialized = false,
        artifacts_materialized = false,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_complete_core_shell_source_plan_payload,
        route_kind = recipe.route_kind,
        route_family,
        support_window_payload_status,
        raw_box_route_payload_status,
        source_realization_payload_status,
        source_box_first = true,
        returns_pqs_multilayer_shell_source_plan = false,
        source_plan_materialized = !isnothing(source_plan),
        shell_support_row_contraction_authority = false,
        retained_diagnostic_weights_are_ida_weights = false,
    )
    return _PQSDiatomicCompleteCoreShellSourcePlanPayload(
        status,
        blocker,
        route_family,
        system_classification,
        bond_axis,
        parent_axis_bundle_object_available,
        route_skeleton_summary,
        source_box_summary,
        retained_unit_summary,
        pair_inventory_summary,
        center_summary,
        coulomb_expansion_summary,
        support_window_payload_value,
        source_realization_payload_value,
        source_plan,
        source_plan_status,
        available_objects,
        missing_objects,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_diatomic_complete_core_shell_final_basis_payload(
    route_skeleton,
    recipe,
    support_window_payload,
    raw_box_route_payload,
    source_realization_payload,
    source_plan_payload,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    source_plan =
        !isnothing(source_plan_payload) && hasproperty(source_plan_payload, :source_plan) ?
        source_plan_payload.source_plan :
        nothing
    source_plan_status =
        isnothing(source_plan) ?
        :not_available :
        source_plan.status
    final_basis = nothing
    final_basis_status = :not_materialized_diatomic_complete_core_shell_final_basis
    final_basis_requested =
        get(recipe, :run_final_basis, false) ||
        get(recipe, :run_h1, false) ||
        get(recipe, :run_h1_j, false) ||
        get(get(recipe, :private_rhf_inputs, (;)), :run_private_rhf, false)
    available = Symbol[]
    missing = Symbol[]

    if route_family !== :pqs_source_box
        status =
            :not_applicable_diatomic_complete_core_shell_final_basis_non_pqs_route
        blocker = nothing
    elseif source_plan_status !==
           :available_pqs_diatomic_complete_core_shell_source_plan
        status = :blocked_diatomic_complete_core_shell_final_basis_payload
        blocker =
            source_plan_status === :not_available ?
            :missing_pqs_diatomic_complete_core_shell_source_plan :
            :pqs_diatomic_complete_core_shell_source_plan_not_available
        push!(missing, :pqs_diatomic_complete_core_shell_source_plan)
        !isnothing(source_plan_payload) &&
            append!(missing, source_plan_payload.missing_objects)
    elseif !final_basis_requested
        status = :not_requested_diatomic_complete_core_shell_final_basis_payload
        blocker = nothing
        push!(available, :pqs_diatomic_complete_core_shell_source_plan)
        push!(missing, :diatomic_complete_core_shell_final_basis_request)
    else
        metrics = source_plan.metrics
        core_overlap = _pqs_multilayer_support_product_matrix(
            source_plan.core_support_states,
            source_plan.core_support_states,
            metrics.x.overlap,
            metrics.y.overlap,
            metrics.z.overlap,
        )
        core_shell_overlap = _pqs_multilayer_support_product_matrix(
            source_plan.core_support_states,
            source_plan.shell_support_states,
            metrics.x.overlap,
            metrics.y.overlap,
            metrics.z.overlap,
        )
        shell_overlap = _pqs_multilayer_support_product_matrix(
            source_plan.shell_support_states,
            source_plan.shell_support_states,
            metrics.x.overlap,
            metrics.y.overlap,
            metrics.z.overlap,
        )
        final_basis =
            CartesianFinalBasisRealization.pqs_complete_core_shell_final_basis(
                core_support_indices = source_plan.core_support_indices,
                shell_support_indices = source_plan.shell_support_indices,
                core_overlap = core_overlap,
                core_shell_overlap = core_shell_overlap,
                shell_overlap = shell_overlap,
                shell_final_coefficients =
                    source_plan.shell_final_coefficients,
                metadata = merge(
                    NamedTuple(source_plan.metadata),
                    (;
                        source =
                            :pqs_source_box_route_driver_diatomic_complete_core_shell_final_basis_payload,
                        input_source_plan =
                            :pqs_diatomic_complete_core_shell_source_plan,
                        retained_pre_final_map =
                            source_plan.retained_pre_final_map,
                        old_source_plan_object_kind = false,
                        final_basis_materialized = true,
                        h1_materialized = false,
                        h1_j_materialized = false,
                        ham_payload_materialized = false,
                        route_driver_public_surface = false,
                        exports_materialized = false,
                        artifacts_materialized = false,
                    ),
                ),
            )
        final_basis_status = final_basis.status
        if final_basis_status === :available_pqs_complete_core_shell_final_basis
            status =
                :available_diatomic_complete_core_shell_final_basis_payload
            blocker = nothing
            append!(
                available,
                (
                    :pqs_diatomic_complete_core_shell_source_plan,
                    :diatomic_complete_core_shell_final_basis,
                ),
            )
            push!(missing, :diatomic_complete_core_shell_h1_consumer)
        else
            status = :blocked_diatomic_complete_core_shell_final_basis_payload
            blocker =
                isnothing(final_basis.blocker) ?
                :diatomic_complete_core_shell_final_basis_blocked :
                final_basis.blocker
            push!(available, :pqs_diatomic_complete_core_shell_source_plan)
            push!(missing, :diatomic_complete_core_shell_final_basis)
        end
    end

    available_objects = Tuple(unique(available))
    missing_objects = Tuple(unique(missing))
    final_dimension =
        !isnothing(final_basis) &&
        hasproperty(final_basis, :final_retained_count) ?
        final_basis.final_retained_count :
        nothing
    core_support_count =
        !isnothing(source_plan) ? length(source_plan.core_support_indices) : nothing
    shell_support_count =
        !isnothing(source_plan) ? length(source_plan.shell_support_indices) : nothing
    shell_retained_count =
        !isnothing(source_plan) ? size(source_plan.shell_final_coefficients, 2) : nothing
    precleanup_retained_dimension =
        !isnothing(source_plan) &&
        hasproperty(source_plan, :summary) ?
        source_plan.summary.precleanup_retained_dimension :
        nothing
    support_row_order =
        !isnothing(final_basis) &&
        hasproperty(final_basis, :support_row_order) ?
        final_basis.support_row_order :
        nothing
    summary = (;
        status,
        blocker,
        route_family,
        source_plan_status,
        final_basis_status,
        final_dimension,
        precleanup_retained_dimension,
        core_support_count,
        shell_support_count,
        shell_retained_count,
        support_row_order,
        old_source_plan_object_kind = false,
        shell_coefficient_block_structure =
            !isnothing(source_plan) ?
            source_plan.summary.shell_coefficient_block_structure :
            nothing,
        final_basis_materialized =
            final_basis_status === :available_pqs_complete_core_shell_final_basis,
        h1_materialized = false,
        h1_j_materialized = false,
        ham_payload_materialized = false,
        route_driver_public_surface = false,
        exports_materialized = false,
        artifacts_materialized = false,
        available_objects,
        missing_objects,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_complete_core_shell_final_basis_payload,
        route_kind = recipe.route_kind,
        route_family,
        source_plan_status,
        final_basis_status,
        support_window_payload_status =
            isnothing(support_window_payload) ? :not_available : support_window_payload.status,
        raw_box_route_payload_status =
            isnothing(raw_box_route_payload) ? :not_available : raw_box_route_payload.status,
        source_realization_payload_status =
            isnothing(source_realization_payload) ? :not_available : source_realization_payload.status,
        input_source_plan = :pqs_diatomic_complete_core_shell_source_plan,
        old_source_plan_object_kind = false,
        final_basis_materialized = summary.final_basis_materialized,
        h1_materialized = false,
        h1_j_materialized = false,
        ham_payload_materialized = false,
        route_driver_public_surface = false,
        exports_materialized = false,
        artifacts_materialized = false,
    )
    return _PQSDiatomicCompleteCoreShellFinalBasisPayload(
        status,
        blocker,
        route_family,
        source_plan,
        source_plan_status,
        final_basis,
        final_basis_status,
        available_objects,
        missing_objects,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_diatomic_support_kinetic_matrix(source_plan)
    states = vcat(source_plan.core_support_states, source_plan.shell_support_states)
    metrics = source_plan.metrics
    support_operator =
        _pqs_multilayer_support_product_matrix(
            states,
            states,
            metrics.x.kinetic,
            metrics.y.overlap,
            metrics.z.overlap,
        ) +
        _pqs_multilayer_support_product_matrix(
            states,
            states,
            metrics.x.overlap,
            metrics.y.kinetic,
            metrics.z.overlap,
        ) +
        _pqs_multilayer_support_product_matrix(
            states,
            states,
            metrics.x.overlap,
            metrics.y.overlap,
            metrics.z.kinetic,
        )
    return (;
        object_kind = :pqs_diatomic_complete_core_shell_support_kinetic_matrix,
        status = :materialized_diatomic_complete_core_shell_support_kinetic_matrix,
        blocker = nothing,
        term = :kinetic,
        support_operator,
        support_operator_shape = size(support_operator),
        support_operator_finite = all(isfinite, support_operator),
        support_state_count = length(states),
        support_ordering = :core_product_then_shell_left_right_pqs,
        final_basis_transfer_materialized = false,
        h1_materialized = false,
        h1_j_materialized = false,
        density_density_materialized = false,
        rhf_materialized = false,
        route_driver_public_surface = false,
        exports_materialized = false,
        artifacts_materialized = false,
        metadata = (;
            source = :pqs_source_box_route_driver_diatomic_support_kinetic_matrix,
            input_source_plan = :pqs_diatomic_complete_core_shell_source_plan,
            support_row_order = :core_product_then_shell_left_right_pqs,
            old_source_plan_object_kind = false,
            shell_support_row_contraction_authority = false,
        ),
    )
end

function _pqs_source_box_route_driver_diatomic_support_electron_nuclear_by_center_matrices(
    source_plan;
    coulomb_expansion,
    center_records,
    axis_layers,
)
    centers = _pqs_multilayer_center_records(center_records)
    isempty(centers) &&
        throw(ArgumentError("diatomic complete core/shell H1 requires nuclear centers"))
    coefficients = Float64.(coulomb_expansion.coefficients)
    states = vcat(source_plan.core_support_states, source_plan.shell_support_states)
    records = map(enumerate(centers)) do (center_index, center_record)
        center = _pqs_multilayer_center_summary(center_record)
        factor_x, factor_y, factor_z, factor_source =
            _pqs_multilayer_centered_factor_terms(
                axis_layers,
                coulomb_expansion,
                center,
            )
        axis_terms =
            _pqs_multilayer_validate_factor_terms(
                (factor_x, factor_y, factor_z),
                length(coefficients),
            )
        support_operator =
            _pqs_multilayer_support_electron_nuclear_matrix(
                states,
                axis_terms,
                coefficients,
            )
        (;
            object_kind =
                :pqs_diatomic_complete_core_shell_support_electron_nuclear_by_center_matrix,
            status =
                :materialized_diatomic_complete_core_shell_support_electron_nuclear_by_center_matrix,
            blocker = nothing,
            term = :electron_nuclear_by_center,
            support_operator,
            support_operator_shape = size(support_operator),
            support_operator_finite = all(isfinite, support_operator),
            support_state_count = length(states),
            support_ordering = :core_product_then_shell_left_right_pqs,
            by_center = true,
            center_key = center.center_key,
            center_index = center.center_index,
            center_location = center.location,
            location = center.location,
            charge = center.nuclear_charge,
            nuclear_charge = center.nuclear_charge,
            nuclear_charge_recorded = true,
            nuclear_charge_applied = false,
            centers_summed = false,
            center_summation = false,
            uncharged_by_center_convention = true,
            gaussian_factor_terms_source = factor_source,
            gaussian_term_count = length(coefficients),
            final_basis_transfer_materialized = false,
            h1_materialized = false,
            h1_j_materialized = false,
            density_density_materialized = false,
            rhf_materialized = false,
            route_driver_public_surface = false,
            exports_materialized = false,
            artifacts_materialized = false,
            metadata = (;
                source =
                    :pqs_source_box_route_driver_diatomic_support_electron_nuclear_by_center_matrices,
                input_source_plan =
                    :pqs_diatomic_complete_core_shell_source_plan,
                physical_operator = :electron_nuclear_attraction,
                negative_unit_charge_attraction = true,
                nuclear_charge_recorded = true,
                nuclear_charge_applied = false,
                centers_summed = false,
                uncharged_by_center_convention = true,
                charge_application_stage = :hamiltonian_assembly,
                gaussian_factor_terms_source = factor_source,
                old_source_plan_object_kind = false,
                shell_support_row_contraction_authority = false,
                wl_matrix_authority_used = false,
                center_key = center.center_key,
                center_index = center.center_index,
                center_location = center.location,
                nuclear_charge = center.nuclear_charge,
            ),
        )
    end
    return (;
        object_kind =
            :pqs_diatomic_complete_core_shell_support_electron_nuclear_by_center_matrix_set,
        status =
            :materialized_diatomic_complete_core_shell_support_electron_nuclear_by_center_matrix_set,
        blocker = nothing,
        term = :electron_nuclear_by_center,
        records,
        center_count = length(records),
        support_state_count = length(states),
        support_ordering = :core_product_then_shell_left_right_pqs,
        by_center = true,
        nuclear_charge_recorded =
            all(record -> record.nuclear_charge_recorded, records),
        nuclear_charge_applied = false,
        centers_summed = false,
        uncharged_by_center_convention = true,
        final_basis_transfer_materialized = false,
        h1_materialized = false,
        h1_j_materialized = false,
        density_density_materialized = false,
        rhf_materialized = false,
        route_driver_public_surface = false,
        exports_materialized = false,
        artifacts_materialized = false,
    )
end

function _pqs_source_box_route_driver_diatomic_complete_core_shell_h1_payload(
    parent,
    route_skeleton,
    recipe,
    source_plan_payload = nothing,
    final_basis_payload = nothing,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    source_plan =
        !isnothing(source_plan_payload) && hasproperty(source_plan_payload, :source_plan) ?
        source_plan_payload.source_plan :
        nothing
    source_plan_status =
        isnothing(source_plan) ?
        :not_available :
        source_plan.status
    final_basis =
        !isnothing(final_basis_payload) && hasproperty(final_basis_payload, :final_basis) ?
        final_basis_payload.final_basis :
        nothing
    final_basis_status =
        isnothing(final_basis) ?
        :not_available :
        final_basis.status

    support_kinetic = nothing
    support_electron_nuclear_by_center = nothing
    final_kinetic = nothing
    final_electron_nuclear_by_center = nothing
    final_hamiltonian = nothing
    h1 = nothing
    support_kinetic_status = :not_materialized_diatomic_support_kinetic_matrix
    support_electron_nuclear_status =
        :not_materialized_diatomic_support_electron_nuclear_by_center
    final_kinetic_status = :not_materialized_diatomic_final_kinetic
    final_electron_nuclear_status =
        :not_materialized_diatomic_final_electron_nuclear_by_center
    final_hamiltonian_status = :not_materialized_diatomic_final_h1_hamiltonian
    h1_status = :not_materialized_diatomic_complete_core_shell_h1
    h1_requested =
        get(recipe, :run_h1, false) ||
        get(recipe, :run_h1_j, false) ||
        get(get(recipe, :private_rhf_inputs, (;)), :run_private_rhf, false)
    available = Symbol[]
    missing = Symbol[]

    if route_family !== :pqs_source_box
        status = :not_applicable_diatomic_complete_core_shell_h1_non_pqs_route
        blocker = nothing
    elseif isnothing(source_plan) ||
           source_plan_status !==
           :available_pqs_diatomic_complete_core_shell_source_plan
        status = :blocked_diatomic_complete_core_shell_h1_payload
        blocker = :missing_pqs_diatomic_complete_core_shell_source_plan
        push!(missing, :pqs_diatomic_complete_core_shell_source_plan)
    elseif source_plan.object_kind !== :pqs_diatomic_complete_core_shell_source_plan
        status = :blocked_diatomic_complete_core_shell_h1_payload
        blocker = :unexpected_diatomic_complete_core_shell_source_plan_object_kind
        push!(missing, :pqs_diatomic_complete_core_shell_source_plan)
    elseif source_plan.support_order !== (:product, :pqs_left, :pqs_right)
        status = :blocked_diatomic_complete_core_shell_h1_payload
        blocker = :unexpected_diatomic_complete_core_shell_support_order
        push!(missing, :diatomic_complete_core_shell_support_order)
    elseif get(source_plan.convention_labels, :support_row_order, nothing) !==
           :core_product_then_shell_left_right_pqs
        status = :blocked_diatomic_complete_core_shell_h1_payload
        blocker = :unexpected_diatomic_complete_core_shell_support_row_order
        push!(missing, :diatomic_complete_core_shell_support_row_order)
    elseif isnothing(final_basis) ||
           final_basis_status !== :available_pqs_complete_core_shell_final_basis
        status = :blocked_diatomic_complete_core_shell_h1_payload
        blocker = :missing_diatomic_complete_core_shell_final_basis
        push!(available, :pqs_diatomic_complete_core_shell_source_plan)
        push!(missing, :diatomic_complete_core_shell_final_basis)
    elseif !h1_requested
        status = :not_requested_diatomic_complete_core_shell_h1_payload
        blocker = nothing
        push!(available, :pqs_diatomic_complete_core_shell_source_plan)
        push!(available, :diatomic_complete_core_shell_final_basis)
        push!(missing, :diatomic_complete_core_shell_h1_request)
    elseif get(final_basis, :support_row_order, nothing) !== :core_then_shell
        status = :blocked_diatomic_complete_core_shell_h1_payload
        blocker = :diatomic_complete_core_shell_final_basis_support_order_mismatch
        push!(available, :pqs_diatomic_complete_core_shell_source_plan)
        push!(missing, :diatomic_complete_core_shell_final_basis_support_order)
    else
        push!(
            available,
            :pqs_diatomic_complete_core_shell_source_plan,
            :diatomic_complete_core_shell_final_basis,
        )
        center_records, missing_centers =
            _pqs_source_box_route_driver_complete_core_shell_center_records(parent)
        axis_layers, missing_axis_layers =
            _pqs_source_box_route_driver_complete_core_shell_axis_layers(
                source_plan.bundles,
            )
        missing_inputs = (missing_centers..., missing_axis_layers...)
        if !isempty(missing_inputs)
            status = :blocked_diatomic_complete_core_shell_h1_payload
            blocker = :missing_diatomic_complete_core_shell_h1_inputs
            append!(missing, missing_inputs)
        else
            coulomb_expansion = coulomb_gaussian_expansion(doacc = false)
            metadata = (;
                source =
                    :pqs_source_box_route_driver_diatomic_complete_core_shell_h1_payload,
                route_kind = recipe.route_kind,
                route_family,
                input_source_plan =
                    :pqs_diatomic_complete_core_shell_source_plan,
                old_source_plan_object_kind = false,
                source_box_first = true,
                shell_support_row_contraction_authority = false,
                retained_diagnostic_weights_are_ida_weights = false,
            )
            try
                support_kinetic =
                    _pqs_source_box_route_driver_diatomic_support_kinetic_matrix(
                        source_plan,
                    )
                support_kinetic_status = support_kinetic.status
                support_electron_nuclear_by_center =
                    _pqs_source_box_route_driver_diatomic_support_electron_nuclear_by_center_matrices(
                        source_plan;
                        coulomb_expansion,
                        center_records,
                        axis_layers,
                    )
                support_electron_nuclear_status =
                    support_electron_nuclear_by_center.status
                final_kinetic =
                    CartesianFinalBasisRealization.pqs_complete_core_shell_final_one_body_matrix(
                        final_basis,
                        support_kinetic.support_operator;
                        term = :kinetic,
                        metadata,
                    )
                final_kinetic_status = final_kinetic.status
                final_electron_nuclear_by_center =
                    map(support_electron_nuclear_by_center.records) do record
                        CartesianFinalBasisRealization.pqs_complete_core_shell_final_one_body_matrix(
                            final_basis,
                            record.support_operator;
                            term = :electron_nuclear_by_center,
                            center_record = record,
                            metadata = merge(
                                record.metadata,
                                (;
                                    nuclear_factor_source =
                                        :centered_axis_layer_gaussian_factor_terms,
                                    support_gaussian_factor_terms_source =
                                        record.gaussian_factor_terms_source,
                                    input_source_plan =
                                        :pqs_diatomic_complete_core_shell_source_plan,
                                    old_source_plan_object_kind = false,
                                ),
                            ),
                        )
                    end
                final_electron_nuclear_status =
                    all(
                        result ->
                            result.status ===
                            :materialized_pqs_complete_core_shell_final_one_body_matrix,
                        final_electron_nuclear_by_center,
                    ) ?
                    :materialized_diatomic_final_electron_nuclear_by_center :
                    :blocked_diatomic_final_electron_nuclear_by_center
                final_hamiltonian =
                    CartesianFinalBasisRealization.pqs_complete_core_shell_final_one_electron_hamiltonian(
                        final_kinetic,
                        final_electron_nuclear_by_center,
                    )
                final_hamiltonian_status = final_hamiltonian.status
                h1 =
                    CartesianFinalBasisRealization.pqs_complete_core_shell_final_h1_solve(
                        final_hamiltonian,
                    )
                h1_status = h1.status
                status = :available_diatomic_complete_core_shell_h1_payload
                blocker = nothing
                push!(available, :diatomic_complete_core_shell_h1_payload)
            catch error
                error isa ArgumentError || error isa DimensionMismatch || rethrow()
                status = :blocked_diatomic_complete_core_shell_h1_payload
                blocker = :diatomic_complete_core_shell_h1_payload_error
                push!(missing, :diatomic_complete_core_shell_h1_payload)
                h1_status = :blocked_argument_error
                metadata = merge(
                    metadata,
                    (error_message = sprint(showerror, error),),
                )
            end
        end
    end

    available_objects = Tuple(unique(available))
    missing_objects = Tuple(unique(missing))
    final_dimension =
        !isnothing(h1) && hasproperty(h1, :final_dimension) ?
        h1.final_dimension :
        !isnothing(final_basis) && hasproperty(final_basis, :final_retained_count) ?
        final_basis.final_retained_count :
        nothing
    lowest_energy =
        !isnothing(h1) && hasproperty(h1, :lowest_energy) ? h1.lowest_energy : nothing
    summary = (;
        status,
        blocker,
        route_family,
        source_plan_status,
        final_basis_status,
        support_kinetic_status,
        support_electron_nuclear_status,
        final_kinetic_status,
        final_electron_nuclear_status,
        final_hamiltonian_status,
        h1_status,
        final_dimension,
        lowest_energy,
        center_count =
            !isnothing(support_electron_nuclear_by_center) ?
            support_electron_nuclear_by_center.center_count :
            nothing,
        support_row_order =
            !isnothing(final_basis) && hasproperty(final_basis, :support_row_order) ?
            final_basis.support_row_order :
            nothing,
        source_plan_support_row_order =
            !isnothing(source_plan) &&
            hasproperty(source_plan, :convention_labels) ?
            get(source_plan.convention_labels, :support_row_order, nothing) :
            nothing,
        old_source_plan_object_kind = false,
        h1_materialized =
            status === :available_diatomic_complete_core_shell_h1_payload,
        h1_j_materialized = false,
        density_density_materialized = false,
        ham_payload_materialized = false,
        rhf_materialized = false,
        route_driver_public_surface = false,
        public_api = false,
        exports_materialized = false,
        artifacts_materialized = false,
        available_objects,
        missing_objects,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_complete_core_shell_h1_payload,
        route_kind = recipe.route_kind,
        route_family,
        source_plan_status,
        final_basis_status,
        support_kinetic_status,
        support_electron_nuclear_status,
        final_kinetic_status,
        final_electron_nuclear_status,
        final_hamiltonian_status,
        h1_status,
        source_box_first = true,
        old_source_plan_object_kind = false,
        h1_materialized = summary.h1_materialized,
        h1_j_materialized = false,
        density_density_materialized = false,
        ham_payload_materialized = false,
        rhf_materialized = false,
        route_driver_public_surface = false,
        public_api = false,
        exports_materialized = false,
        artifacts_materialized = false,
        shell_support_row_contraction_authority = false,
        retained_diagnostic_weights_are_ida_weights = false,
    )
    return _PQSDiatomicCompleteCoreShellH1Payload(
        status,
        blocker,
        route_family,
        source_plan,
        source_plan_status,
        final_basis,
        final_basis_status,
        support_kinetic,
        support_kinetic_status,
        support_electron_nuclear_by_center,
        support_electron_nuclear_status,
        final_kinetic,
        final_kinetic_status,
        final_electron_nuclear_by_center,
        final_electron_nuclear_status,
        final_hamiltonian,
        final_hamiltonian_status,
        h1,
        h1_status,
        available_objects,
        missing_objects,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_diatomic_support_states(source_plan)
    return vcat(source_plan.core_support_states, source_plan.shell_support_states)
end

function _pqs_source_box_route_driver_diatomic_support_weights(
    source_plan;
    axis_weights,
)
    states = _pqs_source_box_route_driver_diatomic_support_states(source_plan)
    weights_x, weights_y, weights_z =
        _pqs_multilayer_common_or_axis_tuple(axis_weights, "axis_weights")
    weights = Vector{Float64}(undef, length(states))
    @inbounds for (index, (ix, iy, iz)) in pairs(states)
        weights[index] =
            Float64(weights_x[ix]) * Float64(weights_y[iy]) * Float64(weights_z[iz])
    end
    return (;
        object_kind = :pqs_diatomic_complete_core_shell_support_weights,
        status = :materialized_diatomic_complete_core_shell_support_weights,
        blocker = nothing,
        support_weights = weights,
        support_weight_count = length(weights),
        support_weights_all_finite = all(isfinite, weights),
        support_weights_all_positive = all(>(0.0), weights),
        support_weight_min = minimum(weights),
        support_weight_max = maximum(weights),
        support_weight_sum = sum(weights),
        support_state_count = length(states),
        support_ordering = :core_product_then_shell_left_right_pqs,
        final_basis_transfer_materialized = false,
        h1_j_materialized = false,
        density_density_materialized = false,
        rhf_materialized = false,
        route_driver_public_surface = false,
        exports_materialized = false,
        artifacts_materialized = false,
        metadata = (;
            source = :pqs_source_box_route_driver_diatomic_support_weights,
            input_source_plan = :pqs_diatomic_complete_core_shell_source_plan,
            support_row_order = :core_product_then_shell_left_right_pqs,
            old_source_plan_object_kind = false,
            retained_diagnostic_weights_are_ida_weights = false,
        ),
    )
end

function _pqs_source_box_route_driver_diatomic_support_pair_raw_numerator_matrix(
    source_plan;
    raw_pair_factor_terms,
    coulomb_expansion,
)
    states = _pqs_source_box_route_driver_diatomic_support_states(source_plan)
    coefficients = Float64.(coulomb_expansion.coefficients)
    all(>(0.0), coefficients) ||
        throw(ArgumentError("diatomic support raw pair numerator requires positive Coulomb expansion coefficients"))
    axis_terms =
        _pqs_multilayer_validate_raw_pair_factor_terms(
            _pqs_multilayer_raw_pair_factor_terms(raw_pair_factor_terms),
            length(coefficients),
        )
    support_pair_raw_numerator = zeros(Float64, length(states), length(states))
    for term_index in eachindex(coefficients)
        support_pair_raw_numerator .+=
            coefficients[term_index] *
            _pqs_multilayer_support_product_matrix(
                states,
                states,
                @view(axis_terms[1][term_index, :, :]),
                @view(axis_terms[2][term_index, :, :]),
                @view(axis_terms[3][term_index, :, :]),
            )
    end
    symmetry_error =
        norm(support_pair_raw_numerator - transpose(support_pair_raw_numerator), Inf)
    return (;
        object_kind =
            :pqs_diatomic_complete_core_shell_support_raw_pair_numerator_matrix,
        status =
            :materialized_diatomic_complete_core_shell_support_raw_pair_numerator_matrix,
        blocker = nothing,
        support_pair_raw_numerator,
        support_pair_raw_numerator_shape = size(support_pair_raw_numerator),
        support_pair_raw_numerator_finite = all(isfinite, support_pair_raw_numerator),
        support_pair_raw_numerator_symmetry_error = symmetry_error,
        support_state_count = length(states),
        support_ordering = :core_product_then_shell_left_right_pqs,
        raw_pair_factor_convention = :raw_numerator,
        density_normalized_pair_terms_used_as_authority = false,
        h1_j_materialized = false,
        density_density_materialized = false,
        rhf_materialized = false,
        route_driver_public_surface = false,
        exports_materialized = false,
        artifacts_materialized = false,
        metadata = (;
            source =
                :pqs_source_box_route_driver_diatomic_support_pair_raw_numerator_matrix,
            input_source_plan = :pqs_diatomic_complete_core_shell_source_plan,
            raw_pair_factor_convention = :raw_numerator,
            support_row_order = :core_product_then_shell_left_right_pqs,
            old_source_plan_object_kind = false,
            density_normalized_pair_terms_used_as_authority = false,
        ),
    )
end

function _pqs_source_box_route_driver_diatomic_density_provenance(
    source_plan,
    coulomb_expansion,
)
    expected_term_count = length(coulomb_expansion.coefficients)
    provenance =
        CartesianContractedParentMetrics._pqs_source_box_ida_factor_provenance(
            source_plan.bundles;
            expected_term_count,
        )
    hasproperty(provenance, :axis_weights) ||
        throw(ArgumentError("diatomic support density provenance missing axis weights"))
    hasproperty(provenance, :raw_axis_pair_factor_terms) ||
        throw(ArgumentError("diatomic support density provenance missing raw pair factor terms"))
    return (;
        object_kind = :pqs_diatomic_complete_core_shell_density_provenance,
        status = :available_diatomic_complete_core_shell_density_provenance,
        blocker = nothing,
        provenance,
        axis_weights = provenance.axis_weights,
        raw_pair_factor_terms = provenance.raw_axis_pair_factor_terms,
        term_count = provenance.term_count,
        factor_dimensions = provenance.factor_dimensions,
        axis_weights_available = true,
        raw_pair_factor_terms_available = true,
        private_diagnostic_only = true,
        retained_pqs_weights_used = false,
        density_normalized_pair_terms_used_as_authority = false,
        metadata = (;
            source = :pqs_source_box_route_driver_diatomic_density_provenance,
            provenance_source = :pqs_source_box_ida_factor_provenance,
            input_source_plan = :pqs_diatomic_complete_core_shell_source_plan,
            expected_term_count,
            retained_diagnostic_weights_are_ida_weights = false,
        ),
    )
end

function _pqs_source_box_route_driver_diatomic_complete_core_shell_ham_input_payload(
    route_skeleton,
    recipe,
    source_plan_payload = nothing,
    final_basis_payload = nothing,
    h1_payload = nothing,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    source_plan =
        !isnothing(source_plan_payload) && hasproperty(source_plan_payload, :source_plan) ?
        source_plan_payload.source_plan :
        nothing
    source_plan_status =
        isnothing(source_plan) ? :not_available : source_plan.status
    final_basis =
        !isnothing(final_basis_payload) && hasproperty(final_basis_payload, :final_basis) ?
        final_basis_payload.final_basis :
        nothing
    final_basis_status =
        isnothing(final_basis) ? :not_available : final_basis.status
    h1_payload_status =
        isnothing(h1_payload) ? :not_available : h1_payload.status

    density_provenance = nothing
    support_weights = nothing
    raw_pair_factor_terms = nothing
    support_pair_raw_numerator = nothing
    density_interaction = nothing
    density_provenance_status =
        :not_materialized_diatomic_complete_core_shell_density_provenance
    support_weights_status =
        :not_materialized_diatomic_complete_core_shell_support_weights
    raw_pair_factor_status = :not_materialized_diatomic_raw_pair_factor_terms
    support_pair_raw_numerator_status =
        :not_materialized_diatomic_support_raw_pair_numerator_matrix
    density_interaction_status =
        :not_materialized_diatomic_pre_final_density_interaction
    available = Symbol[]
    missing = Symbol[]

    if route_family !== :pqs_source_box
        status =
            :not_applicable_diatomic_complete_core_shell_ham_input_non_pqs_route
        blocker = nothing
    elseif isnothing(source_plan) ||
           source_plan_status !==
           :available_pqs_diatomic_complete_core_shell_source_plan
        status = :blocked_diatomic_complete_core_shell_ham_input_payload
        blocker = :missing_diatomic_complete_core_shell_source_plan
        push!(missing, :pqs_diatomic_complete_core_shell_source_plan)
    elseif isnothing(final_basis) ||
           final_basis_status !== :available_pqs_complete_core_shell_final_basis
        status = :blocked_diatomic_complete_core_shell_ham_input_payload
        blocker = :missing_diatomic_complete_core_shell_final_basis
        push!(available, :pqs_diatomic_complete_core_shell_source_plan)
        push!(missing, :diatomic_complete_core_shell_final_basis)
    elseif isnothing(h1_payload) ||
           h1_payload_status !== :available_diatomic_complete_core_shell_h1_payload
        status = :blocked_diatomic_complete_core_shell_ham_input_payload
        blocker = :missing_diatomic_complete_core_shell_h1_payload
        append!(
            available,
            (
                :pqs_diatomic_complete_core_shell_source_plan,
                :diatomic_complete_core_shell_final_basis,
            ),
        )
        push!(missing, :diatomic_complete_core_shell_h1_payload)
    elseif get(final_basis, :support_row_order, nothing) !== :core_then_shell
        status = :blocked_diatomic_complete_core_shell_ham_input_payload
        blocker = :diatomic_complete_core_shell_final_basis_support_order_mismatch
        push!(missing, :diatomic_complete_core_shell_final_basis_support_order)
    elseif get(source_plan.convention_labels, :support_row_order, nothing) !==
           :core_product_then_shell_left_right_pqs
        status = :blocked_diatomic_complete_core_shell_ham_input_payload
        blocker = :unexpected_diatomic_complete_core_shell_support_row_order
        push!(missing, :diatomic_complete_core_shell_support_row_order)
    else
        append!(
            available,
            (
                :pqs_diatomic_complete_core_shell_source_plan,
                :diatomic_complete_core_shell_final_basis,
                :diatomic_complete_core_shell_h1_payload,
            ),
        )
        coulomb_expansion = coulomb_gaussian_expansion(doacc = false)
        try
            density_provenance =
                _pqs_source_box_route_driver_diatomic_density_provenance(
                    source_plan,
                    coulomb_expansion,
                )
            density_provenance_status = density_provenance.status
            push!(available, :diatomic_complete_core_shell_density_provenance)
            raw_pair_factor_terms = density_provenance.raw_pair_factor_terms
            raw_pair_factor_status = :available_diatomic_raw_pair_factor_terms
            support_weights =
                _pqs_source_box_route_driver_diatomic_support_weights(
                    source_plan;
                    axis_weights = density_provenance.axis_weights,
                )
            support_weights_status = support_weights.status
            push!(available, :diatomic_complete_core_shell_support_weights)
            support_pair_raw_numerator =
                _pqs_source_box_route_driver_diatomic_support_pair_raw_numerator_matrix(
                    source_plan;
                    raw_pair_factor_terms,
                    coulomb_expansion,
                )
            support_pair_raw_numerator_status =
                support_pair_raw_numerator.status
            push!(
                available,
                :diatomic_complete_core_shell_support_raw_pair_numerator,
            )
            density_interaction =
                CartesianFinalBasisRealization.pqs_complete_core_shell_pre_final_density_interaction(
                    final_basis,
                    support_pair_raw_numerator.support_pair_raw_numerator,
                    support_weights.support_weights;
                    metadata = (;
                        source =
                            :pqs_source_box_route_driver_diatomic_complete_core_shell_ham_input_payload,
                        input_source_plan =
                            :pqs_diatomic_complete_core_shell_source_plan,
                        support_density_input_source =
                            :pqs_source_box_route_driver_diatomic_density_provenance,
                        raw_pair_factor_convention = :raw_numerator,
                        old_source_plan_object_kind = false,
                        source_box_first = true,
                        retained_diagnostic_weights_are_ida_weights = false,
                    ),
                )
            density_interaction_status = density_interaction.status
            if density_interaction_status ===
               :materialized_pqs_complete_core_shell_pre_final_density_interaction
                status =
                    :available_diatomic_complete_core_shell_ham_input_payload
                blocker = nothing
                push!(
                    available,
                    :diatomic_complete_core_shell_ham_input_payload,
                    :diatomic_complete_core_shell_pre_final_density_interaction,
                )
            else
                status = :blocked_diatomic_complete_core_shell_ham_input_payload
                blocker = :blocked_diatomic_pre_final_density_interaction
                push!(missing, :diatomic_complete_core_shell_pre_final_density_interaction)
            end
        catch error
            error isa ArgumentError || error isa DimensionMismatch || rethrow()
            status = :blocked_diatomic_complete_core_shell_ham_input_payload
            message = sprint(showerror, error)
            if occursin("axis weights", message)
                blocker = :missing_diatomic_support_weights
                push!(missing, :diatomic_complete_core_shell_support_weights)
            elseif occursin("raw pair factor", message)
                blocker = :missing_diatomic_raw_pair_factor_terms
                push!(missing, :diatomic_raw_pair_factor_terms)
            else
                blocker = :missing_diatomic_support_density_provenance
                push!(missing, :diatomic_complete_core_shell_density_provenance)
            end
        end
    end

    available_objects = Tuple(unique(available))
    missing_objects = Tuple(unique(missing))
    final_dimension =
        !isnothing(final_basis) && hasproperty(final_basis, :final_retained_count) ?
        final_basis.final_retained_count :
        nothing
    pre_final_dimension =
        !isnothing(density_interaction) &&
        hasproperty(density_interaction, :pre_final_weight_count) ?
        density_interaction.pre_final_weight_count :
        nothing
    density_gauge =
        !isnothing(density_interaction) &&
        hasproperty(density_interaction, :density_gauge) ?
        density_interaction.density_gauge :
        nothing
    raw_pair_factor_convention =
        !isnothing(support_pair_raw_numerator) ?
        support_pair_raw_numerator.raw_pair_factor_convention :
        nothing
    summary = (;
        status,
        blocker,
        route_family,
        source_plan_status,
        final_basis_status,
        h1_payload_status,
        density_provenance_status,
        support_weights_status,
        raw_pair_factor_status,
        support_pair_raw_numerator_status,
        density_interaction_status,
        density_gauge,
        raw_pair_factor_convention,
        final_dimension,
        pre_final_dimension,
        support_row_order =
            !isnothing(final_basis) && hasproperty(final_basis, :support_row_order) ?
            final_basis.support_row_order :
            nothing,
        source_plan_support_row_order =
            !isnothing(source_plan) &&
            hasproperty(source_plan, :convention_labels) ?
            get(source_plan.convention_labels, :support_row_order, nothing) :
            nothing,
        pre_final_pair_matrix_shape =
            !isnothing(density_interaction) &&
            hasproperty(density_interaction, :pre_final_pair_matrix_shape) ?
            density_interaction.pre_final_pair_matrix_shape :
            nothing,
        support_weight_count =
            !isnothing(support_weights) ? support_weights.support_weight_count : nothing,
        support_raw_pair_shape =
            !isnothing(support_pair_raw_numerator) ?
            support_pair_raw_numerator.support_pair_raw_numerator_shape :
            nothing,
        ham_input_materialized =
            status === :available_diatomic_complete_core_shell_ham_input_payload,
        available_objects,
        missing_objects,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_complete_core_shell_ham_input_payload,
        route_kind = recipe.route_kind,
        route_family,
        source_plan_status,
        final_basis_status,
        h1_payload_status,
        density_provenance_status,
        support_weights_status,
        raw_pair_factor_status,
        support_pair_raw_numerator_status,
        density_interaction_status,
        electron_electron_representation = :pre_final_density_interaction,
        density_gauge,
        raw_pair_factor_convention,
        source_box_first = true,
        old_source_plan_object_kind = false,
        ham_input_materialized = summary.ham_input_materialized,
    )
    return _PQSDiatomicCompleteCoreShellHamInputPayload(
        status,
        blocker,
        route_family,
        source_plan,
        source_plan_status,
        final_basis,
        final_basis_status,
        h1_payload,
        h1_payload_status,
        density_provenance,
        density_provenance_status,
        support_weights,
        support_weights_status,
        raw_pair_factor_terms,
        raw_pair_factor_status,
        support_pair_raw_numerator,
        support_pair_raw_numerator_status,
        density_interaction,
        density_interaction_status,
        available_objects,
        missing_objects,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_diatomic_handoff_center_metadata(parent)
    center_records, missing_centers =
        _pqs_source_box_route_driver_complete_core_shell_center_records(parent)
    if !isempty(missing_centers) || isnothing(center_records)
        return (;
            status = :blocked_diatomic_handoff_center_metadata,
            blocker = :missing_diatomic_center_records,
            center_records = nothing,
            center_count = 0,
            nuclear_charges = (),
            nuclear_coordinates = (),
            nuclear_repulsion_status = :missing_diatomic_nuclear_repulsion,
            nuclear_repulsion = nothing,
            nuclear_repulsion_source = nothing,
            electron_count_status = :missing_diatomic_electron_count_convention,
            electron_count = nothing,
            electron_count_source = nothing,
            spin_sector_status = :missing_diatomic_spin_sector_convention,
            spin_sector = nothing,
            spin_sector_source = nothing,
            missing_objects = (
                :diatomic_center_records,
                :diatomic_nuclear_repulsion,
                :diatomic_electron_count_convention,
                :diatomic_spin_sector_convention,
            ),
        )
    end

    records = Tuple(center_records)
    nuclear_charges = Tuple(Float64(record.nuclear_charge) for record in records)
    nuclear_coordinates = Tuple(
        ntuple(axis -> Float64(record.location[axis]), 3) for record in records
    )
    missing = Symbol[]
    nuclear_repulsion_status = :missing_diatomic_nuclear_repulsion
    nuclear_repulsion = nothing
    nuclear_repulsion_source = nothing
    if length(records) == 2
        distance = sqrt(
            sum(
                (
                    nuclear_coordinates[1][axis] -
                    nuclear_coordinates[2][axis]
                )^2 for axis in 1:3
            ),
        )
        if isfinite(distance) && distance > 0.0
            nuclear_repulsion =
                nuclear_charges[1] * nuclear_charges[2] / distance
            nuclear_repulsion_status =
                :available_diatomic_nuclear_repulsion
            nuclear_repulsion_source = :center_record_charge_distance
        else
            push!(missing, :diatomic_nuclear_repulsion)
        end
    else
        push!(missing, :diatomic_nuclear_repulsion)
    end

    electron_count_status = :missing_diatomic_electron_count_convention
    electron_count = nothing
    electron_count_source = nothing
    spin_sector_status = :missing_diatomic_spin_sector_convention
    spin_sector = nothing
    spin_sector_source = nothing
    if length(records) == 2 && nuclear_charges == (4.0, 4.0)
        electron_count = 8
        electron_count_status = :available_diatomic_electron_count_convention
        electron_count_source = :neutral_sum_nuclear_charges_private_route_smoke
        spin_sector = :closed_shell_singlet
        spin_sector_status = :available_diatomic_spin_sector_convention
        spin_sector_source = :private_be2_route_smoke_default
    else
        push!(missing, :diatomic_electron_count_convention)
        push!(missing, :diatomic_spin_sector_convention)
    end

    status =
        isempty(missing) &&
        nuclear_repulsion_status === :available_diatomic_nuclear_repulsion &&
        electron_count_status ===
        :available_diatomic_electron_count_convention &&
        spin_sector_status === :available_diatomic_spin_sector_convention ?
        :available_diatomic_handoff_center_metadata :
        :blocked_diatomic_handoff_center_metadata
    blocker =
        status === :available_diatomic_handoff_center_metadata ?
        nothing :
        first(missing)
    return (;
        status,
        blocker,
        center_records = records,
        center_count = length(records),
        nuclear_charges,
        nuclear_coordinates,
        nuclear_repulsion_status,
        nuclear_repulsion,
        nuclear_repulsion_source,
        electron_count_status,
        electron_count,
        electron_count_source,
        spin_sector_status,
        spin_sector,
        spin_sector_source,
        missing_objects = Tuple(unique(missing)),
    )
end

function _pqs_source_box_route_driver_diatomic_complete_core_shell_hamiltonian_handoff_payload(
    parent,
    route_skeleton,
    recipe,
    source_plan_payload = nothing,
    final_basis_payload = nothing,
    h1_payload = nothing,
    ham_input_payload = nothing,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    source_plan =
        !isnothing(source_plan_payload) && hasproperty(source_plan_payload, :source_plan) ?
        source_plan_payload.source_plan :
        nothing
    source_plan_status =
        isnothing(source_plan) ? :not_available : source_plan.status
    final_basis =
        !isnothing(final_basis_payload) && hasproperty(final_basis_payload, :final_basis) ?
        final_basis_payload.final_basis :
        nothing
    final_basis_status =
        isnothing(final_basis) ? :not_available : final_basis.status
    h1_payload_status =
        isnothing(h1_payload) ? :not_available : h1_payload.status
    ham_input_payload_status =
        isnothing(ham_input_payload) ? :not_available : ham_input_payload.status

    one_body_hamiltonian = nothing
    one_body_hamiltonian_status = :not_available_diatomic_one_body_hamiltonian
    density_interaction = nothing
    density_interaction_status = :not_available_diatomic_density_interaction
    pre_final_pair_matrix = nothing
    final_to_pre_final_coefficients = nothing
    pre_final_weights = nothing
    support_weights = nothing
    support_pair_raw_numerator = nothing
    raw_pair_factor_terms = nothing
    center_records = nothing
    center_metadata = nothing
    nuclear_repulsion_status = :missing_diatomic_nuclear_repulsion
    nuclear_repulsion = nothing
    electron_count_status = :missing_diatomic_electron_count_convention
    electron_count = nothing
    spin_sector_status = :missing_diatomic_spin_sector_convention
    spin_sector = nothing
    available = Symbol[]
    missing = Symbol[]

    if route_family !== :pqs_source_box
        status =
            :not_applicable_diatomic_complete_core_shell_hamiltonian_handoff_non_pqs_route
        blocker = nothing
    elseif isnothing(source_plan) ||
           source_plan_status !==
           :available_pqs_diatomic_complete_core_shell_source_plan
        status =
            :blocked_diatomic_complete_core_shell_hamiltonian_handoff_payload
        blocker = :missing_diatomic_complete_core_shell_source_plan
        push!(missing, :pqs_diatomic_complete_core_shell_source_plan)
    elseif isnothing(final_basis) ||
           final_basis_status !== :available_pqs_complete_core_shell_final_basis
        status =
            :blocked_diatomic_complete_core_shell_hamiltonian_handoff_payload
        blocker = :missing_diatomic_complete_core_shell_final_basis
        push!(available, :pqs_diatomic_complete_core_shell_source_plan)
        push!(missing, :diatomic_complete_core_shell_final_basis)
    elseif isnothing(h1_payload) ||
           h1_payload_status !== :available_diatomic_complete_core_shell_h1_payload
        status =
            :blocked_diatomic_complete_core_shell_hamiltonian_handoff_payload
        blocker = :missing_diatomic_complete_core_shell_h1_payload
        append!(
            available,
            (
                :pqs_diatomic_complete_core_shell_source_plan,
                :diatomic_complete_core_shell_final_basis,
            ),
        )
        push!(missing, :diatomic_complete_core_shell_h1_payload)
    elseif isnothing(ham_input_payload) ||
           ham_input_payload_status !==
           :available_diatomic_complete_core_shell_ham_input_payload
        status =
            :blocked_diatomic_complete_core_shell_hamiltonian_handoff_payload
        blocker = :missing_diatomic_complete_core_shell_ham_input_payload
        append!(
            available,
            (
                :pqs_diatomic_complete_core_shell_source_plan,
                :diatomic_complete_core_shell_final_basis,
                :diatomic_complete_core_shell_h1_payload,
            ),
        )
        push!(missing, :diatomic_complete_core_shell_ham_input_payload)
    else
        append!(
            available,
            (
                :pqs_diatomic_complete_core_shell_source_plan,
                :diatomic_complete_core_shell_final_basis,
                :diatomic_complete_core_shell_h1_payload,
                :diatomic_complete_core_shell_ham_input_payload,
            ),
        )
        one_body_hamiltonian =
            !isnothing(h1_payload.final_hamiltonian) &&
            hasproperty(h1_payload.final_hamiltonian, :hamiltonian_matrix) ?
            h1_payload.final_hamiltonian.hamiltonian_matrix :
            nothing
        if isnothing(one_body_hamiltonian)
            one_body_hamiltonian_status =
                :missing_diatomic_one_body_hamiltonian_reference
            push!(missing, :diatomic_one_body_hamiltonian_reference)
        else
            one_body_hamiltonian_status =
                :available_diatomic_one_body_hamiltonian_reference
            push!(available, :diatomic_one_body_hamiltonian_reference)
        end

        density_interaction = ham_input_payload.density_interaction
        density_interaction_status = ham_input_payload.density_interaction_status
        if density_interaction_status !==
           :materialized_pqs_complete_core_shell_pre_final_density_interaction
            push!(missing, :diatomic_pre_final_density_interaction)
        else
            push!(available, :diatomic_pre_final_density_interaction)
            pre_final_pair_matrix =
                hasproperty(density_interaction, :pre_final_pair_matrix) ?
                density_interaction.pre_final_pair_matrix :
                nothing
            final_to_pre_final_coefficients =
                hasproperty(density_interaction, :final_to_pre_final_coefficients) ?
                density_interaction.final_to_pre_final_coefficients :
                nothing
            pre_final_weights =
                hasproperty(density_interaction, :pre_final_weights) ?
                density_interaction.pre_final_weights :
                nothing
        end
        support_weights =
            !isnothing(ham_input_payload.support_weights) &&
            hasproperty(ham_input_payload.support_weights, :support_weights) ?
            ham_input_payload.support_weights.support_weights :
            nothing
        support_pair_raw_numerator =
            !isnothing(ham_input_payload.support_pair_raw_numerator) &&
            hasproperty(
                ham_input_payload.support_pair_raw_numerator,
                :support_pair_raw_numerator,
            ) ?
            ham_input_payload.support_pair_raw_numerator.support_pair_raw_numerator :
            nothing
        raw_pair_factor_terms = ham_input_payload.raw_pair_factor_terms

        center_metadata =
            _pqs_source_box_route_driver_diatomic_handoff_center_metadata(parent)
        center_records = center_metadata.center_records
        nuclear_repulsion_status = center_metadata.nuclear_repulsion_status
        nuclear_repulsion = center_metadata.nuclear_repulsion
        electron_count_status = center_metadata.electron_count_status
        electron_count = center_metadata.electron_count
        spin_sector_status = center_metadata.spin_sector_status
        spin_sector = center_metadata.spin_sector
        if center_metadata.status === :available_diatomic_handoff_center_metadata
            push!(available, :diatomic_handoff_center_metadata)
        else
            append!(missing, center_metadata.missing_objects)
        end

        if isempty(missing)
            status =
                :available_diatomic_complete_core_shell_hamiltonian_handoff_payload
            blocker = nothing
            push!(
                available,
                :diatomic_complete_core_shell_hamiltonian_handoff_payload,
            )
        else
            status =
                :blocked_diatomic_complete_core_shell_hamiltonian_handoff_payload
            blocker = first(missing)
        end
    end

    available_objects = Tuple(unique(available))
    missing_objects = Tuple(unique(missing))
    final_dimension =
        !isnothing(final_basis) && hasproperty(final_basis, :final_retained_count) ?
        final_basis.final_retained_count :
        nothing
    pre_final_dimension =
        !isnothing(density_interaction) &&
        hasproperty(density_interaction, :pre_final_weight_count) ?
        density_interaction.pre_final_weight_count :
        nothing
    density_gauge =
        !isnothing(density_interaction) &&
        hasproperty(density_interaction, :density_gauge) ?
        density_interaction.density_gauge :
        nothing
    raw_pair_factor_convention =
        !isnothing(ham_input_payload) &&
        hasproperty(ham_input_payload.summary, :raw_pair_factor_convention) ?
        ham_input_payload.summary.raw_pair_factor_convention :
        nothing
    ordering = (;
        support_row_order =
            !isnothing(final_basis) && hasproperty(final_basis, :support_row_order) ?
            final_basis.support_row_order :
            nothing,
        source_plan_support_row_order =
            !isnothing(source_plan) &&
            hasproperty(source_plan, :convention_labels) ?
            get(source_plan.convention_labels, :support_row_order, nothing) :
            nothing,
        source_plan_support_order =
            !isnothing(source_plan) && hasproperty(source_plan, :support_order) ?
            source_plan.support_order :
            nothing,
        route_retained_order =
            !isnothing(source_plan) && hasproperty(source_plan, :route_retained_order) ?
            source_plan.route_retained_order :
            nothing,
        pre_final_order = :pqs_complete_core_shell_pre_final_density_order,
        final_order = :pqs_complete_core_shell_final_basis_order,
    )
    conventions = (;
        density_gauge,
        raw_pair_factor_convention,
        electron_count_source =
            isnothing(center_metadata) ? nothing : center_metadata.electron_count_source,
        spin_sector_source =
            isnothing(center_metadata) ? nothing : center_metadata.spin_sector_source,
        nuclear_repulsion_source =
            isnothing(center_metadata) ? nothing : center_metadata.nuclear_repulsion_source,
        private_inspect_only = true,
    )
    summary = (;
        status,
        blocker,
        route_family,
        source_plan_status,
        final_basis_status,
        h1_payload_status,
        ham_input_payload_status,
        one_body_hamiltonian_status,
        density_interaction_status,
        final_dimension,
        pre_final_dimension,
        support_weight_count =
            !isnothing(ham_input_payload) &&
            hasproperty(ham_input_payload.summary, :support_weight_count) ?
            ham_input_payload.summary.support_weight_count :
            nothing,
        pre_final_pair_matrix_shape =
            !isnothing(density_interaction) &&
            hasproperty(density_interaction, :pre_final_pair_matrix_shape) ?
            density_interaction.pre_final_pair_matrix_shape :
            nothing,
        final_to_pre_final_coefficient_shape =
            isnothing(final_to_pre_final_coefficients) ?
            nothing :
            size(final_to_pre_final_coefficients),
        density_gauge,
        raw_pair_factor_convention,
        support_row_order = ordering.support_row_order,
        source_plan_support_row_order = ordering.source_plan_support_row_order,
        pre_final_order = ordering.pre_final_order,
        final_order = ordering.final_order,
        center_count =
            isnothing(center_metadata) ? nothing : center_metadata.center_count,
        nuclear_charges =
            isnothing(center_metadata) ? () : center_metadata.nuclear_charges,
        nuclear_coordinates =
            isnothing(center_metadata) ? () : center_metadata.nuclear_coordinates,
        nuclear_repulsion_status,
        nuclear_repulsion,
        electron_count_status,
        electron_count,
        spin_sector_status,
        spin_sector,
        available_objects,
        missing_objects,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_complete_core_shell_hamiltonian_handoff_payload,
        route_kind = recipe.route_kind,
        route_family,
        source_plan_status,
        final_basis_status,
        h1_payload_status,
        ham_input_payload_status,
        source_box_first = true,
        shell_support_row_contraction_authority = false,
        retained_diagnostic_weights_are_ida_weights = false,
    )
    return _PQSDiatomicCompleteCoreShellHamiltonianHandoffPayload(
        status,
        blocker,
        route_family,
        source_plan,
        source_plan_status,
        final_basis,
        final_basis_status,
        h1_payload,
        h1_payload_status,
        ham_input_payload,
        ham_input_payload_status,
        one_body_hamiltonian,
        one_body_hamiltonian_status,
        density_interaction,
        density_interaction_status,
        pre_final_pair_matrix,
        final_to_pre_final_coefficients,
        pre_final_weights,
        support_weights,
        support_pair_raw_numerator,
        raw_pair_factor_terms,
        center_records,
        center_metadata,
        nuclear_repulsion_status,
        nuclear_repulsion,
        electron_count_status,
        electron_count,
        spin_sector_status,
        spin_sector,
        ordering,
        conventions,
        available_objects,
        missing_objects,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_diatomic_downstream_ham_nonclaims()
    return (;
        hfdmrg_density_density_ready = false,
        hfdmrg_sliced_ready = false,
        hamv6_export_ready = false,
        cr2_ready = false,
        public_api = false,
        exports_materialized = false,
        artifacts_materialized = false,
        downstream_missing_objects = (
            :hfdmrg_density_density_contract,
            :hfdmrg_sliced_integrals,
            :hamv6_export_contract,
            :cr2_handoff_format,
        ),
    )
end

function _pqs_source_box_route_driver_diatomic_complete_core_shell_hamiltonian_consumer_contract_payload(
    route_skeleton,
    recipe,
    hamiltonian_handoff_payload = nothing,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    source_handoff_status =
        isnothing(hamiltonian_handoff_payload) ?
        :not_available :
        hamiltonian_handoff_payload.status
    source_handoff_available =
        source_handoff_status ===
        :available_diatomic_complete_core_shell_hamiltonian_handoff_payload
    nonclaims = _pqs_source_box_route_driver_diatomic_downstream_ham_nonclaims()
    available = Symbol[]

    if route_family !== :pqs_source_box
        status =
            :not_applicable_diatomic_complete_core_shell_hamiltonian_consumer_contract_non_pqs_route
        blocker = nothing
        missing_objects = ()
    elseif source_handoff_available
        status =
            :available_diatomic_complete_core_shell_hamiltonian_consumer_contract_payload
        blocker = nothing
        push!(
            available,
            :diatomic_complete_core_shell_hamiltonian_handoff_payload,
            :diatomic_hamiltonian_consumer_contract,
        )
        missing_objects = nonclaims.downstream_missing_objects
    else
        status =
            :blocked_diatomic_complete_core_shell_hamiltonian_consumer_contract_payload
        blocker = :missing_diatomic_complete_core_shell_hamiltonian_handoff_payload
        missing_objects = (:diatomic_complete_core_shell_hamiltonian_handoff_payload,)
    end

    available_objects = Tuple(unique(available))
    private_inspector_ready =
        status ===
        :available_diatomic_complete_core_shell_hamiltonian_consumer_contract_payload
    nonclaim_flags = (;
        hfdmrg_density_density_ready = nonclaims.hfdmrg_density_density_ready,
        hfdmrg_sliced_ready = nonclaims.hfdmrg_sliced_ready,
        hamv6_export_ready = nonclaims.hamv6_export_ready,
        cr2_ready = nonclaims.cr2_ready,
        public_api = nonclaims.public_api,
        exports_materialized = nonclaims.exports_materialized,
        artifacts_materialized = nonclaims.artifacts_materialized,
    )
    handoff_summary =
        source_handoff_available ? hamiltonian_handoff_payload.summary : nothing
    cr2_inspection = (;
        cr2_read_only_inspector_ready = private_inspector_ready,
        cr2_solver_ready = false,
        cr2_export_ready = false,
        cr2_handoff_blocker =
            private_inspector_ready ? :missing_cr2_solver_handoff_format : blocker,
        two_body_representation_kind =
            private_inspector_ready ? :pre_final_density_interaction : nothing,
        density_gauge = isnothing(handoff_summary) ? nothing : handoff_summary.density_gauge,
        raw_pair_factor_convention =
            isnothing(handoff_summary) ? nothing : handoff_summary.raw_pair_factor_convention,
    )
    readiness = (;
        private_inspector_ready,
        nonclaim_flags...,
        cr2_inspection...,
        downstream_blocker =
            private_inspector_ready ?
            :missing_hfdmrg_density_density_contract :
            blocker,
        missing_objects,
    )
    summary = (;
        status,
        blocker,
        route_family,
        source_handoff_status,
        private_inspector_ready,
        nonclaim_flags...,
        available_objects,
        missing_objects,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_complete_core_shell_hamiltonian_consumer_contract_payload,
        route_kind = recipe.route_kind,
        route_family,
        source_handoff_status,
        source_box_first = true,
        private_inspect_only = true,
        expiration_condition =
            :replace_when_real_downstream_hamiltonian_consumer_contract_is_chosen,
        shell_support_row_contraction_authority = false,
        retained_diagnostic_weights_are_ida_weights = false,
        nonclaim_flags...,
    )
    return _PQSDiatomicCompleteCoreShellHamiltonianConsumerContractPayload(
        status,
        blocker,
        route_family,
        hamiltonian_handoff_payload,
        source_handoff_status,
        readiness,
        available_objects,
        missing_objects,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_tsv_cell(value)
    value === nothing && return ""
    return replace(string(value), '\t' => ' ', '\n' => ' ')
end

function _pqs_source_box_route_driver_be2_cr2_inspection_bundle_payload(assembly)
    handoff = assembly.diatomic_complete_core_shell_hamiltonian_handoff_payload
    readiness =
        assembly.diatomic_complete_core_shell_hamiltonian_consumer_contract_payload.readiness
    h1_payload = assembly.diatomic_complete_core_shell_h1_payload
    h1_matrix = Matrix{Float64}(handoff.one_body_hamiltonian)
    pair_matrix = Matrix{Float64}(handoff.pre_final_pair_matrix)
    coefficients = Matrix{Float64}(handoff.final_to_pre_final_coefficients)
    final_interaction_matrix =
        transpose(coefficients) * pair_matrix * coefficients
    support_raw = Matrix{Float64}(handoff.support_pair_raw_numerator)
    hf_convention_blocker =
        :missing_reviewed_density_density_hf_fock_energy_convention
    h1_symmetry_defect = norm(h1_matrix - h1_matrix')
    two_body_symmetry_defect = norm(pair_matrix - pair_matrix')
    final_two_body_symmetry_defect =
        norm(final_interaction_matrix - final_interaction_matrix')
    pqs_fingerprint = (route_label = :pqs_source_box, status = handoff.status,
        blocker = handoff.blocker, final_dimension = size(h1_matrix, 1),
        pre_final_dimension = size(pair_matrix, 1),
        support_weight_count = length(handoff.support_weights),
        one_body_shape = size(h1_matrix), two_body_shape = size(pair_matrix),
        h1_lowest = h1_payload.summary.lowest_energy, h1_symmetry_defect,
        one_body_finite = all(isfinite, h1_matrix),
        two_body_finite = all(isfinite, pair_matrix),
        density_gauge = handoff.summary.density_gauge,
        raw_pair_factor_convention = handoff.summary.raw_pair_factor_convention,
        cr2_read_only_inspector_ready = readiness.cr2_read_only_inspector_ready,
        cr2_solver_ready = readiness.cr2_solver_ready,
        cr2_export_ready = readiness.cr2_export_ready,
        cr2_handoff_blocker = readiness.cr2_handoff_blocker,
        nuclear_repulsion = handoff.nuclear_repulsion,
        electron_count = handoff.electron_count, spin_sector = handoff.spin_sector)
    wl_fingerprint = (route_label = :white_lindsey, status = :unavailable,
        blocker = :white_lindsey_inspection_not_materialized,
        final_dimension = nothing, pre_final_dimension = nothing,
        support_weight_count = nothing, one_body_shape = nothing,
        two_body_shape = nothing, h1_lowest = nothing,
        h1_symmetry_defect = nothing, one_body_finite = false,
        two_body_finite = false, density_gauge = :not_applicable,
        raw_pair_factor_convention = :not_applicable,
        cr2_read_only_inspector_ready = false, cr2_solver_ready = false,
        cr2_export_ready = false,
        cr2_handoff_blocker = :white_lindsey_inspection_not_materialized,
        nuclear_repulsion = handoff.nuclear_repulsion,
        electron_count = handoff.electron_count, spin_sector = handoff.spin_sector)
    coordinates = reduce(
        vcat,
        (permutedims(collect(coord)) for coord in handoff.summary.nuclear_coordinates),
    )
    unavailable_wl = (status = :unavailable, blocker = :white_lindsey_inspection_not_materialized)
    return (;
        jld2_values = Dict{String,Any}(
            "schema/name" => "be2_wl_pqs_handoff_inspection_bundle",
            "schema/version" => 1,
            "bundle/purpose" => "cr2_read_only_hamiltonian_inspection",
            "producer" => (;
                package = "GaussletBases",
                repo_commit = :unavailable,
                dirty = :unavailable,
                generated_at = :unavailable,
            ),
            "routes/pqs_source_box/route" => (;
                label = :pqs_source_box,
                family = handoff.route_family,
                kind = get(assembly, :low_order_assembly_route_kind, :unavailable),
                status = handoff.status,
                blocker = handoff.blocker,
                density_density_hf_convention_status = hf_convention_blocker,
                density_density_hf_convention_blocker = hf_convention_blocker,
            ),
            "routes/pqs_source_box/readiness" => readiness,
            "routes/pqs_source_box/system" => (;
                nuclear_charges = collect(handoff.summary.nuclear_charges),
                nuclear_coordinates = coordinates,
                nuclear_repulsion = handoff.nuclear_repulsion,
                electron_count = handoff.electron_count,
                spin_sector = handoff.spin_sector,
            ),
            "routes/pqs_source_box/final_basis" => (;
                final_dimension = size(h1_matrix, 1),
                order_label = handoff.ordering.final_order,
                overlap_convention = :orthonormal_identity_by_contract,
                overlap_matrix_stored = false,
                overlap_identity_defect = 0.0,
            ),
            "routes/pqs_source_box/one_body/hamiltonian" => h1_matrix,
            "routes/pqs_source_box/h1/lowest_energy" => h1_payload.summary.lowest_energy,
            "routes/pqs_source_box/two_body" => (;
                representation_kind = :pre_final_density_interaction,
                interaction_matrix = final_interaction_matrix,
                interaction_matrix_representation_kind =
                    :final_basis_density_density_matrix,
                interaction_matrix_derivation =
                    :final_to_pre_final_density_congruence,
                interaction_matrix_formula =
                    :transpose_final_to_pre_final_times_pre_final_pair_times_final_to_pre_final,
                interaction_matrix_shape = size(final_interaction_matrix),
                interaction_matrix_symmetry_defect =
                    final_two_body_symmetry_defect,
                interaction_matrix_finite = all(isfinite, final_interaction_matrix),
                pre_final_pair_matrix = pair_matrix,
                final_to_pre_final_coefficients = coefficients,
                pre_final_weights = Vector{Float64}(handoff.pre_final_weights),
                support_weights = Vector{Float64}(handoff.support_weights),
                support_raw_pair_numerator = support_raw,
                density_gauge = handoff.summary.density_gauge,
                raw_pair_factor_convention =
                    handoff.summary.raw_pair_factor_convention,
            ),
            "routes/pqs_source_box/hf_convention" => (;
                density_density_hf_convention_status = hf_convention_blocker,
                density_density_hf_convention_blocker = hf_convention_blocker,
            ),
            "routes/pqs_source_box/validation" => (;
                h1_symmetry_defect,
                two_body_symmetry_defect,
                final_two_body_symmetry_defect,
                one_body_finite = pqs_fingerprint.one_body_finite,
                two_body_finite = pqs_fingerprint.two_body_finite,
                final_two_body_finite = all(isfinite, final_interaction_matrix),
            ),
            "routes/white_lindsey/route" => (;
                label = :white_lindsey,
                family = :white_lindsey_low_order,
                kind = :not_materialized,
                status = :unavailable,
                blocker = :white_lindsey_inspection_not_materialized,
                density_density_hf_convention_status = hf_convention_blocker,
                density_density_hf_convention_blocker = hf_convention_blocker,
            ),
            "routes/white_lindsey/readiness" => (;
                cr2_read_only_inspector_ready = false,
                cr2_solver_ready = false,
                cr2_export_ready = false,
                cr2_handoff_blocker = :white_lindsey_inspection_not_materialized,
            ),
            "routes/white_lindsey/system" => merge(
                unavailable_wl,
                (;
                    nuclear_charges = collect(handoff.summary.nuclear_charges),
                    nuclear_coordinates = coordinates,
                    nuclear_repulsion = handoff.nuclear_repulsion,
                    electron_count = handoff.electron_count,
                    spin_sector = handoff.spin_sector,
                ),
            ),
            "routes/white_lindsey/final_basis" => merge(
                unavailable_wl,
                (final_dimension = nothing, order_label = :not_applicable),
            ),
            "routes/white_lindsey/one_body" => merge(
                unavailable_wl,
                (hamiltonian = Float64[], representation_kind = :not_materialized),
            ),
            "routes/white_lindsey/two_body" => merge(
                unavailable_wl,
                (;
                    representation_kind = :not_applicable,
                    interaction_matrix_representation_kind =
                        :final_basis_density_density_matrix,
                    pre_final_pair_matrix = Float64[],
                    final_to_pre_final_coefficients = Float64[],
                    pre_final_weights = Float64[],
                    support_weights = Float64[],
                    support_raw_pair_numerator = Float64[],
                    density_gauge = :not_applicable,
                    raw_pair_factor_convention = :not_applicable,
                ),
            ),
            "routes/white_lindsey/hf_convention" => (;
                density_density_hf_convention_status = hf_convention_blocker,
                density_density_hf_convention_blocker = hf_convention_blocker,
            ),
            "routes/white_lindsey/validation" => merge(
                unavailable_wl,
                (;
                    h1_symmetry_defect = nothing,
                    two_body_symmetry_defect = nothing,
                    one_body_finite = false,
                    two_body_finite = false,
                ),
            ),
        ),
        rows = (pqs_fingerprint, wl_fingerprint),
        fingerprint_columns = propertynames(pqs_fingerprint),
    )
end

function _pqs_source_box_route_driver_write_be2_cr2_inspection_bundle(
    jld2_path, tsv_path, payload)
    jldopen(jld2_path, "w") do file
        for (key, value) in payload.jld2_values
            _cartesian_write_value!(file, key, value)
        end
    end
    open(tsv_path, "w") do io
        columns = payload.fingerprint_columns
        println(io, join(String.(columns), '\t'))
        for row in payload.rows
            println(io, join((
                _pqs_source_box_route_driver_tsv_cell(getproperty(row, col)) for
                col in columns), '\t'))
        end
    end
    return (jld2_path = jld2_path, tsv_path = tsv_path)
end

function _pqs_source_box_route_driver_diatomic_complete_core_shell_ham_readiness_payload(
    parent,
    route_skeleton,
    recipe,
    source_plan_payload = nothing,
    final_basis_payload = nothing,
    h1_payload = nothing,
    ham_input_payload = nothing,
    hamiltonian_handoff_payload = nothing,
    hamiltonian_consumer_contract_payload = nothing,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    system_classification =
        hasproperty(parent, :system_classification) ?
        parent.system_classification :
        nothing
    bond_axis =
        hasproperty(parent, :bond_axis) ? parent.bond_axis : nothing
    parent_axis_bundle_object_available =
        hasproperty(parent, :parent_axis_bundle_object_available) ?
        parent.parent_axis_bundle_object_available :
        (
            hasproperty(parent, :parent_axis_bundle_object) &&
            !isnothing(parent.parent_axis_bundle_object)
        )

    center_summary =
        isnothing(source_plan_payload) ?
        _pqs_source_box_route_driver_diatomic_center_summary(parent) :
        source_plan_payload.center_summary
    source_box_summary =
        isnothing(source_plan_payload) ?
        _pqs_source_box_route_driver_diatomic_source_box_summary(route_skeleton) :
        source_plan_payload.source_box_summary
    retained_unit_summary =
        isnothing(source_plan_payload) ?
        _pqs_source_box_route_driver_diatomic_retained_unit_summary(route_skeleton) :
        source_plan_payload.retained_unit_summary
    pair_inventory_summary =
        isnothing(source_plan_payload) ?
        _pqs_source_box_route_driver_diatomic_pair_inventory_summary(route_skeleton) :
        source_plan_payload.pair_inventory_summary
    route_skeleton_summary =
        isnothing(source_plan_payload) ?
        _pqs_source_box_route_driver_diatomic_route_skeleton_summary(
            route_skeleton,
            retained_unit_summary,
        ) :
        source_plan_payload.route_skeleton_summary
    source_plan_payload_status =
        isnothing(source_plan_payload) ? :not_available : source_plan_payload.status
    source_plan_status =
        !isnothing(source_plan_payload) &&
        hasproperty(source_plan_payload, :source_plan_status) ?
        source_plan_payload.source_plan_status :
        :not_available
    final_basis_payload_status =
        isnothing(final_basis_payload) ? :not_available : final_basis_payload.status
    final_basis_status =
        !isnothing(final_basis_payload) &&
        hasproperty(final_basis_payload, :final_basis_status) ?
        final_basis_payload.final_basis_status :
        :not_available
    h1_payload_status =
        isnothing(h1_payload) ? :not_available : h1_payload.status
    h1_status =
        !isnothing(h1_payload) && hasproperty(h1_payload, :h1_status) ?
        h1_payload.h1_status :
        :not_available
    h1_available =
        h1_payload_status === :available_diatomic_complete_core_shell_h1_payload
    ham_input_payload_status =
        isnothing(ham_input_payload) ? :not_available : ham_input_payload.status
    ham_input_available =
        ham_input_payload_status ===
        :available_diatomic_complete_core_shell_ham_input_payload
    hamiltonian_handoff_payload_status =
        isnothing(hamiltonian_handoff_payload) ?
        :not_available :
        hamiltonian_handoff_payload.status
    hamiltonian_handoff_available =
        hamiltonian_handoff_payload_status ===
        :available_diatomic_complete_core_shell_hamiltonian_handoff_payload
    hamiltonian_consumer_contract_payload_status =
        isnothing(hamiltonian_consumer_contract_payload) ?
        :not_available :
        hamiltonian_consumer_contract_payload.status
    hamiltonian_consumer_contract_available =
        hamiltonian_consumer_contract_payload_status ===
        :available_diatomic_complete_core_shell_hamiltonian_consumer_contract_payload
    nonclaims = _pqs_source_box_route_driver_diatomic_downstream_ham_nonclaims()
    nonclaim_flags = (;
        hfdmrg_density_density_ready = nonclaims.hfdmrg_density_density_ready,
        hfdmrg_sliced_ready = nonclaims.hfdmrg_sliced_ready,
        hamv6_export_ready = nonclaims.hamv6_export_ready,
        cr2_ready = nonclaims.cr2_ready,
        public_api = nonclaims.public_api,
        exports_materialized = nonclaims.exports_materialized,
        artifacts_materialized = nonclaims.artifacts_materialized,
    )

    if route_family !== :pqs_source_box
        status =
            :not_applicable_diatomic_complete_core_shell_ham_readiness_non_pqs_route
        blocker = nothing
        available_objects = ()
        missing_objects = ()
    elseif system_classification !== :bond_aligned_diatomic
        status =
            :not_applicable_diatomic_complete_core_shell_ham_readiness_non_diatomic
        blocker = nothing
        available_objects = ()
        missing_objects = ()
    else
        status = :blocked_diatomic_complete_core_shell_ham_readiness
        blocker =
            hamiltonian_consumer_contract_available ?
            :missing_hfdmrg_density_density_contract :
            hamiltonian_handoff_available ?
            :missing_diatomic_hamiltonian_consumer_contract :
            ham_input_available ?
            :missing_diatomic_complete_core_shell_hamiltonian_handoff_payload :
            h1_available ?
            :missing_diatomic_complete_core_shell_ham_input_payload :
            final_basis_status ===
            :available_pqs_complete_core_shell_final_basis ?
            :missing_diatomic_complete_core_shell_h1_consumer :
            source_plan_status ===
            :available_pqs_diatomic_complete_core_shell_source_plan ?
            :missing_diatomic_complete_core_shell_final_basis_consumer :
            :missing_diatomic_complete_core_shell_source_plan_producer
        available = Symbol[
            :route_skeleton,
            :diatomic_center_metadata,
            :source_boxes,
            :retained_units,
            :pair_inventory,
        ]
        source_plan_status ===
            :available_pqs_diatomic_complete_core_shell_source_plan &&
            push!(available, :pqs_diatomic_complete_core_shell_source_plan)
        final_basis_status === :available_pqs_complete_core_shell_final_basis &&
            push!(available, :diatomic_complete_core_shell_final_basis)
        h1_available &&
            push!(available, :diatomic_complete_core_shell_h1_payload)
        ham_input_available &&
            push!(available, :diatomic_complete_core_shell_ham_input_payload)
        hamiltonian_handoff_available &&
            push!(
                available,
                :diatomic_complete_core_shell_hamiltonian_handoff_payload,
            )
        hamiltonian_consumer_contract_available &&
            push!(available, :diatomic_hamiltonian_consumer_contract)
        parent_axis_bundle_object_available &&
            push!(available, :parent_axis_bundle_object)
        missing =
            hamiltonian_consumer_contract_available ?
            Symbol[nonclaims.downstream_missing_objects...] :
            hamiltonian_handoff_available ?
            Symbol[:diatomic_hamiltonian_consumer_contract] :
            ham_input_available ?
            Symbol[:diatomic_complete_core_shell_hamiltonian_handoff_payload] :
            h1_available ?
            Symbol[:diatomic_complete_core_shell_ham_input_payload] :
            final_basis_status ===
            :available_pqs_complete_core_shell_final_basis ?
            Symbol[
                :diatomic_complete_core_shell_h1_consumer,
                :complete_core_shell_density_inputs,
                :complete_core_shell_h1_j_diagnostic_payload,
            ] :
            source_plan_status ===
            :available_pqs_diatomic_complete_core_shell_source_plan ?
            Symbol[:diatomic_complete_core_shell_final_basis_consumer] :
            Symbol[
                :diatomic_complete_core_shell_source_plan_producer,
                :pqs_multilayer_shell_region_plan,
                :pqs_multilayer_shell_source_plan,
                :pqs_multilayer_complete_core_shell_final_basis,
                :pqs_multilayer_complete_core_shell_h1_payload,
                :complete_core_shell_density_inputs,
                :complete_core_shell_h1_j_diagnostic_payload,
            ]
        parent_axis_bundle_object_available ||
            push!(missing, :parent_axis_bundle_object)
        available_objects = Tuple(available)
        missing_objects = Tuple(missing)
    end

    summary = (;
        status,
        blocker,
        route_family,
        system_classification,
        bond_axis,
        center_count = center_summary.center_count,
        source_box_count = source_box_summary.source_box_count,
        retained_unit_count = retained_unit_summary.retained_unit_count,
        pair_count = pair_inventory_summary.pair_count,
        parent_axis_bundle_object_available,
        source_plan_payload_status,
        source_plan_status,
        final_basis_payload_status,
        final_basis_status,
        h1_payload_status,
        h1_status,
        ham_input_payload_status,
        hamiltonian_handoff_payload_status,
        hamiltonian_consumer_contract_payload_status,
        missing_objects,
        private_internal_only = true,
        final_basis_materialized =
            final_basis_status === :available_pqs_complete_core_shell_final_basis,
        h1_materialized = h1_available,
        ham_input_materialized = ham_input_available,
        hamiltonian_handoff_materialized = hamiltonian_handoff_available,
        hamiltonian_consumer_contract_materialized =
            hamiltonian_consumer_contract_available,
        nonclaim_flags...,
        h1_j_materialized = false,
        ham_payload_materialized = false,
        rhf_materialized = false,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_complete_core_shell_ham_readiness_payload,
        route_kind = recipe.route_kind,
        route_family,
        source_plan_payload_status,
        final_basis_payload_status,
        h1_payload_status,
        ham_input_payload_status,
        hamiltonian_handoff_payload_status,
        hamiltonian_consumer_contract_payload_status,
        source_box_first = true,
        shell_support_row_contraction_authority = false,
        retained_diagnostic_weights_are_ida_weights = false,
    )
    return _PQSDiatomicCompleteCoreShellHamReadinessPayload(
        status,
        blocker,
        route_family,
        system_classification,
        bond_axis,
        center_summary,
        parent_axis_bundle_object_available,
        route_skeleton_summary,
        source_box_summary,
        retained_unit_summary,
        pair_inventory_summary,
        available_objects,
        missing_objects,
        summary,
        metadata,
    )
end
