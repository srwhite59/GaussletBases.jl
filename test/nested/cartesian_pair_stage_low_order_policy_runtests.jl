using Test
using GaussletBases

function _cartesian_pair_stage_low_order_policy_fixture(;
    probe_parent_axis_construction = :auto,
    atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
    parent_axis_counts = (x = 9, y = 7, z = 9),
)
    system_inputs = (;
        atom_symbols = ("Be", "Be"),
        nuclear_charges = (4, 4),
        atom_locations,
        radius = 15.0,
        parent_axis_counts,
        map_backend = :pgdg_localized_experimental,
    )
    spacing_inputs = (;
        q = 5,
        n_s = 5,
        reference_spacing = 1.0,
        tail_spacing = 10.0,
        q_to_core_spacing_rule = :standard_pqs_ns_equals_q,
        core_spacing = 0.15,
    )
    parent_inputs = (;
        probe_parent_axis_construction,
        parent_axis_probe_backend = :pgdg_localized_experimental,
        parent_axis_probe_family = :G10,
    )
    route_probe_inputs = (;
        probe_raw_product_box_plans = :auto,
        raw_product_box_probe_backend = :pgdg_localized_experimental,
    )
    route_inputs = (;
        route_family = :white_lindsey_low_order,
        route_kind = :be2_cartesian_pair_stage_low_order_policy,
        route_shape = (:pqs_left, :product, :pqs_right),
        product_body_rule = :centered_single_z_slab,
        pqs_retained_rule = :boundary_comx_product_mode_selection,
        product_retained_rule = :product_doside_retained_unit,
        terms = (
            :overlap,
            :position_x,
            :position_y,
            :position_z,
            :x2_x,
            :x2_y,
            :x2_z,
            :kinetic,
        ),
        pair_factor_normalization = :density_normalized,
        support_dense_direct_allowed = false,
        reference_only_authorities =
            (:support_row_oracle, :dense_parent_projection),
        white_lindsey_route_shape =
            (:standard_cartesian_units, :low_order_comx_coarsening),
        white_lindsey_mapping_rule = :standard_unit_backbone_mapping_family,
        white_lindsey_nesting_rule = :unit_box_low_order_comx_coarsening,
        white_lindsey_retained_rule = :low_order_unit_comx_retained_basis,
        white_lindsey_operator_rule = :low_order_unit_operator_blocks,
        white_lindsey_benchmark_role =
            :published_cartesian_baseline_for_pqs_comparison,
    )

    system = GaussletBases.cartesian_system(system_inputs)
    recipe = GaussletBases.cartesian_recipe(route_inputs)
    parent = GaussletBases.cartesian_parent(
        system,
        spacing_inputs,
        parent_inputs,
        recipe,
    )
    return (; parent, spacing_inputs, recipe, route_probe_inputs)
end

function _cartesian_pair_stage_low_order_policy_pairs(fixture; policy = nothing)
    shells = isnothing(policy) ?
        GaussletBases.cartesian_shells(
            fixture.parent,
            fixture.spacing_inputs,
            fixture.recipe,
        ) :
        GaussletBases.cartesian_shells(
            fixture.parent,
            fixture.spacing_inputs,
            fixture.recipe;
            low_order_shellization_policy = policy,
        )
    units = GaussletBases.cartesian_units(
        fixture.parent,
        shells,
        fixture.route_probe_inputs,
        fixture.recipe,
    )
    transforms = GaussletBases.cartesian_transforms(units, fixture.recipe)
    pairs = GaussletBases.cartesian_pair_terms(
        units,
        transforms,
        fixture.recipe,
    )
    return (; shells, units, transforms, pairs)
end

function _cartesian_pair_stage_count_by_field(entries, field, value)
    for entry in entries
        getproperty(entry, field) == value && return entry.pair_count
    end
    return 0
end

function _cartesian_pair_stage_property(record, field, default = nothing)
    return hasproperty(record, field) ? getproperty(record, field) : default
end

function _cartesian_pair_stage_fingerprint(record, spec)
    names = Tuple(item[1] for item in spec)
    values = Tuple(
        _cartesian_pair_stage_property(record, item[2], item[3])
        for item in spec
    )
    return NamedTuple{names}(values)
end

function _cartesian_pair_stage_contract_counts_fingerprint(entries)
    return Tuple(
        (;
            unit_key = _cartesian_pair_stage_property(entry, :unit_key),
            lowering_contract_count = _cartesian_pair_stage_property(
                entry,
                :lowering_contract_count,
            ),
            selected_contract_count = _cartesian_pair_stage_property(
                entry,
                :selected_contract_count,
            ),
        )
        for entry in entries
    )
end

const _CARTESIAN_PAIR_STAGE_TERMINAL_UNIT_FIELDS = (
    (:available, :terminal_shellification_unit_inventory_available, false),
    (:region_count, :terminal_shellification_region_count, 0),
    (:unit_count, :terminal_shellification_unit_count, 0),
    (:unit_keys, :terminal_shellification_unit_keys, ()),
    (:unit_roles, :terminal_shellification_unit_roles, ()),
    (:unit_kinds, :terminal_shellification_unit_kinds, ()),
    (:support_counts, :terminal_shellification_unit_support_counts, ()),
    (:central_gap_region_count,
        :terminal_shellification_central_gap_region_count, 0),
    (:central_midpoint_slab_count,
        :terminal_shellification_central_midpoint_slab_count, 0),
    (:central_distorted_product_box_count,
        :terminal_shellification_central_distorted_product_box_count, 0),
)

function _cartesian_pair_stage_terminal_unit_inventory_fingerprint(stage)
    return _cartesian_pair_stage_fingerprint(
        stage,
        _CARTESIAN_PAIR_STAGE_TERMINAL_UNIT_FIELDS,
    )
end

const _CARTESIAN_PAIR_STAGE_LOWERING_CONTRACT_FIELDS = (
    (:available,
        :terminal_shellification_lowering_contract_inventory_available, false),
    (:status,
        :terminal_shellification_lowering_contract_inventory_status, nothing),
    (:count, :terminal_shellification_lowering_contract_count, 0),
    (:kinds, :terminal_shellification_lowering_contract_kinds, ()),
    (:kind_counts,
        :terminal_shellification_lowering_contract_kind_counts, ()),
    (:lw_complete_shell_cpb_count,
        :terminal_shellification_lw_complete_shell_cpb_count, 0),
    (:lw_complete_shell_cpb_family_counts,
        :terminal_shellification_lw_complete_shell_cpb_family_counts, ()),
)

function _cartesian_pair_stage_lowering_contract_fingerprint(stage)
    base = _cartesian_pair_stage_fingerprint(
        stage,
        _CARTESIAN_PAIR_STAGE_LOWERING_CONTRACT_FIELDS,
    )
    return merge(
        base,
        (;
            counts_by_unit = _cartesian_pair_stage_contract_counts_fingerprint(
                _cartesian_pair_stage_property(
                    stage,
                    :terminal_shellification_contract_counts_by_unit,
                    (),
                ),
            ),
        ),
    )
end

const _CARTESIAN_PAIR_STAGE_SELECTED_CONTRACT_FIELDS = (
    (:available,
        :terminal_shellification_selected_lowering_contract_inventory_available,
        false),
    (:status,
        :terminal_shellification_selected_lowering_contract_inventory_status,
        nothing),
    (:family, :terminal_shellification_selected_lowering_family, nothing),
    (:count, :terminal_shellification_selected_contract_count, 0),
    (:kinds, :terminal_shellification_selected_contract_kinds, ()),
    (:kind_counts, :terminal_shellification_selected_contract_kind_counts, ()),
    (:all_units_have_exactly_one_selected_contract,
        :terminal_shellification_all_units_have_exactly_one_selected_contract,
        false),
    (:unselected_count, :terminal_shellification_unselected_contract_count, 0),
    (:unselected_kinds, :terminal_shellification_unselected_contract_kinds, ()),
)

function _cartesian_pair_stage_selected_contract_fingerprint(stage)
    base = _cartesian_pair_stage_fingerprint(
        stage,
        _CARTESIAN_PAIR_STAGE_SELECTED_CONTRACT_FIELDS,
    )
    return merge(
        base,
        (;
            counts_by_unit = _cartesian_pair_stage_contract_counts_fingerprint(
                _cartesian_pair_stage_property(
                    stage,
                    :terminal_shellification_selected_contract_counts_by_unit,
                    (),
                ),
            ),
        ),
    )
end

const _CARTESIAN_PAIR_STAGE_CRC_SIDECAR_FIELDS = (
    (:object_kind, :object_kind, nothing),
    (:status, :status, nothing),
    (:selected_contract_count, :selected_contract_count, 0),
    (:sidecar_available_count, :sidecar_available_count, 0),
    (:sidecar_missing_count, :sidecar_missing_count, 0),
    (:sidecar_inventory_complete, :sidecar_inventory_complete, false),
    (:missing_sidecar_reasons, :missing_sidecar_reasons, ()),
    (:missing_sidecar_kinds, :missing_sidecar_kinds, ()),
    (:final_retained_unit_inventory_available,
        :final_retained_unit_inventory_available, false),
    (:pair_inventory_available, :pair_inventory_available, false),
    (:pair_inventory_status, :pair_inventory_status, nothing),
    (:operator_blocks_materialized, :operator_blocks_materialized, false),
    (:pair_operator_blocks_materialized,
        :pair_operator_blocks_materialized, false),
    (:hamiltonian_data_materialized, :hamiltonian_data_materialized, false),
    (:artifacts_materialized, :artifacts_materialized, false),
)

function _cartesian_pair_stage_crc_sidecar_summary_fingerprint(summary)
    return _cartesian_pair_stage_fingerprint(
        summary,
        _CARTESIAN_PAIR_STAGE_CRC_SIDECAR_FIELDS,
    )
end

const _CARTESIAN_PAIR_STAGE_PAIR_METADATA_FIELDS = (
    (:pair_inventory_known, :pair_inventory_known, false),
    (:pair_inventory_source, :pair_inventory_source, nothing),
    (:pair_family_counts, :pair_family_counts, ()),
    (:helper_by_pair_family, :helper_by_pair_family, ()),
    (:pair_operator_helper_by_family, :pair_operator_helper_by_family, ()),
    (:pair_helper_status_by_family, :pair_helper_status_by_family, ()),
    (:pair_operator_blocks_materialized,
        :pair_operator_blocks_materialized, false),
    (:operator_pairs_materialized, :operator_pairs_materialized, false),
    (:route_core_pair_inventory_available,
        :route_core_pair_inventory_available, false),
    (:route_core_pair_inventory_status,
        :route_core_pair_inventory_status, nothing),
    (:route_core_pair_count, :route_core_pair_count, 0),
    (:route_core_pair_operator_ready, :route_core_pair_operator_ready, false),
    (:route_core_pair_operator_readiness_status,
        :route_core_pair_operator_readiness_status, nothing),
    (:route_core_pair_operator_blocker,
        :route_core_pair_operator_blocker, nothing),
)

function _cartesian_pair_stage_pair_metadata_fingerprint(stage)
    pair_entries = _cartesian_pair_stage_property(stage, :pair_entries, ())
    route_core_pair_keys =
        _cartesian_pair_stage_property(stage, :route_core_pair_keys, ())
    base = _cartesian_pair_stage_fingerprint(
        stage,
        _CARTESIAN_PAIR_STAGE_PAIR_METADATA_FIELDS,
    )
    return merge(
        base,
        (;
            pair_entry_count = length(pair_entries),
            pair_count = _cartesian_pair_stage_property(
                stage,
                :pair_count,
                length(pair_entries),
            ),
            route_core_pair_key_count = length(route_core_pair_keys),
            route_core_pair_keys,
        ),
    )
end

function _assert_terminal_lowering_contract_fields_match_pair_stage(
    pair_stage,
    transforms,
)
    @test _cartesian_pair_stage_lowering_contract_fingerprint(pair_stage) ==
          _cartesian_pair_stage_lowering_contract_fingerprint(transforms)
    @test _cartesian_pair_stage_selected_contract_fingerprint(pair_stage) ==
          _cartesian_pair_stage_selected_contract_fingerprint(transforms)
    transform_selected_crc_sidecars =
        transforms.terminal_shellification_selected_crc_sidecar_summary
    pair_stage_selected_crc_sidecars =
        pair_stage.terminal_shellification_selected_crc_sidecar_summary
    @test _cartesian_pair_stage_crc_sidecar_summary_fingerprint(
        pair_stage_selected_crc_sidecars,
    ) == _cartesian_pair_stage_crc_sidecar_summary_fingerprint(
        transform_selected_crc_sidecars,
    )
    selected_inventory =
        pair_stage.terminal_shellification_selected_lowering_contract_inventory
    @test selected_inventory.selected_contract_count ==
          pair_stage.terminal_shellification_unit_count
    @test selected_inventory.route_lowering_family == :white_lindsey_low_order
    @test selected_inventory.all_units_have_exactly_one_selected_contract
    @test !selected_inventory.final_retained_unit_inventory_available
    @test !selected_inventory.pair_inventory_available
    @test !selected_inventory.coefficient_maps_materialized
    @test !selected_inventory.transform_contracts_materialized
    @test !selected_inventory.retained_spaces_materialized
    @test !selected_inventory.operator_blocks_materialized
    @test !selected_inventory.pair_operator_blocks_materialized
    @test !selected_inventory.hamiltonian_data_materialized
    @test !selected_inventory.artifacts_materialized
end

@testset "cartesian pair-stage compact fingerprint helpers" begin
    synthetic = (;
        terminal_shellification_unit_inventory_available = true,
        terminal_shellification_region_count = 2,
        terminal_shellification_unit_count = 2,
        terminal_shellification_unit_keys = (:left, :right),
        terminal_shellification_unit_roles = (:atom_local, :shared_shell),
        terminal_shellification_unit_kinds = (:direct_core, :complete_shell),
        terminal_shellification_unit_support_counts = (8, 26),
        terminal_shellification_central_gap_region_count = 1,
        terminal_shellification_central_midpoint_slab_count = 1,
        terminal_shellification_central_distorted_product_box_count = 0,
        terminal_shellification_lowering_contract_inventory_available = true,
        terminal_shellification_lowering_contract_inventory_status =
            :available,
        terminal_shellification_lowering_contract_count = 3,
        terminal_shellification_lowering_contract_kinds =
            (:direct_core_identity_cpb, :white_lindsey_boundary_strata),
        terminal_shellification_lowering_contract_kind_counts = (;
            direct_core_identity_cpb_count = 1,
            white_lindsey_boundary_strata_count = 2,
        ),
        terminal_shellification_contract_counts_by_unit = (
            (; unit_key = :left, lowering_contract_count = 1),
        ),
        terminal_shellification_lw_complete_shell_cpb_count = 26,
        terminal_shellification_lw_complete_shell_cpb_family_counts = (;
            facet_cpb_count = 6,
            edge_cpb_count = 12,
            corner_cpb_count = 8,
        ),
        terminal_shellification_selected_lowering_contract_inventory_available =
            true,
        terminal_shellification_selected_lowering_contract_inventory_status =
            :available,
        terminal_shellification_selected_lowering_family =
            :white_lindsey_low_order,
        terminal_shellification_selected_contract_count = 2,
        terminal_shellification_selected_contract_kinds =
            (:direct_core_identity_cpb, :white_lindsey_boundary_strata),
        terminal_shellification_selected_contract_kind_counts = (;
            direct_core_identity_cpb_count = 1,
            white_lindsey_boundary_strata_count = 1,
        ),
        terminal_shellification_selected_contract_counts_by_unit =
            ((; unit_key = :right, selected_contract_count = 1),),
        terminal_shellification_all_units_have_exactly_one_selected_contract =
            true,
        terminal_shellification_unselected_contract_count = 1,
        terminal_shellification_unselected_contract_kinds =
            (:pqs_filled_source_cpb,),
        object_kind =
            :cartesian_unit_stage_selected_terminal_lowering_crc_sidecar_summary,
        status = :available,
        selected_contract_count = 2,
        sidecar_available_count = 2,
        sidecar_missing_count = 0,
        sidecar_inventory_complete = true,
        missing_sidecar_reasons = (),
        missing_sidecar_kinds = (),
        final_retained_unit_inventory_available = false,
        pair_inventory_available = false,
        pair_inventory_status = :deferred,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        pair_inventory_known = true,
        pair_inventory_source = :final_retained_units,
        pair_entries = ((; pair_key = (:left, :left)),),
        pair_count = 1,
        pair_family_counts = (; direct_pair_count = 1),
        helper_by_pair_family = (; direct_pair = :metadata_only),
        pair_operator_helper_by_family = (; direct_pair = :metadata_only),
        pair_helper_status_by_family = (; direct_pair = :metadata_only),
        pair_operator_blocks_materialized = false,
        operator_pairs_materialized = false,
        route_core_pair_inventory_available = true,
        route_core_pair_inventory_status = :available,
        route_core_pair_count = 1,
        route_core_pair_keys = ((:left, :left),),
        route_core_pair_operator_ready = true,
        route_core_pair_operator_readiness_status = :ready,
        route_core_pair_operator_blocker = nothing,
    )

    unit_fingerprint =
        _cartesian_pair_stage_terminal_unit_inventory_fingerprint(synthetic)
    @test unit_fingerprint.unit_count == 2
    @test unit_fingerprint.unit_keys == (:left, :right)

    lowering_fingerprint =
        _cartesian_pair_stage_lowering_contract_fingerprint(synthetic)
    @test lowering_fingerprint.count == 3
    @test lowering_fingerprint.counts_by_unit == (
        (;
            unit_key = :left,
            lowering_contract_count = 1,
            selected_contract_count = nothing,
        ),
    )

    selected_fingerprint =
        _cartesian_pair_stage_selected_contract_fingerprint(synthetic)
    @test selected_fingerprint.family == :white_lindsey_low_order
    @test selected_fingerprint.unselected_kinds == (:pqs_filled_source_cpb,)
    @test selected_fingerprint.all_units_have_exactly_one_selected_contract

    crc_fingerprint =
        _cartesian_pair_stage_crc_sidecar_summary_fingerprint(synthetic)
    @test crc_fingerprint.sidecar_inventory_complete
    @test !crc_fingerprint.pair_inventory_available

    pair_fingerprint =
        _cartesian_pair_stage_pair_metadata_fingerprint(synthetic)
    @test pair_fingerprint.pair_entry_count == 1
    @test pair_fingerprint.route_core_pair_key_count == 1
    @test pair_fingerprint.route_core_pair_operator_ready
    @test !pair_fingerprint.pair_operator_blocks_materialized
end

@testset "cartesian pair stage carries selected low-order policy (slow integration)" begin
    fixture = _cartesian_pair_stage_low_order_policy_fixture()

    default_stages = _cartesian_pair_stage_low_order_policy_pairs(fixture)
    default_pairs = default_stages.pairs
    default_summary = default_pairs.low_order_pairs
    @test default_pairs.object_kind == :cartesian_pair_terms
    @test default_summary.object_kind == :cartesian_pair_stage_low_order_summary
    @test default_summary.low_order_shellization_policy_resolved ==
          :legacy_diatomic_source
    @test default_summary.shellization_source ==
          :route_configured_bond_aligned_diatomic_source
    @test default_summary.transform_route_kind ==
          :legacy_diatomic_source_low_order_transforms
    @test default_summary.pair_route_kind ==
          :legacy_diatomic_source_low_order_pairs
    @test default_summary.legacy_source_pairs_selected
    @test default_summary.active_source_authority
    @test !default_summary.atom_growth_pairs_selected
    @test !default_summary.terminal_shellification_pairs_selected
    @test !default_summary.pair_operator_blocks_materialized
    @test default_summary.pair_inventory_known
    @test default_summary.pair_inventory_source ==
          :route_skeleton_pair_entries_only
    @test default_summary.route_skeleton_pair_entry_count ==
          length(default_pairs.pair_entries)
    @test default_summary.route_skeleton_pair_family_counts ==
          default_pairs.pair_family_counts
    @test !default_summary.independent_atom_growth_pair_inventory_available
    @test !default_summary.route_core_pair_inventory_available
    @test default_summary.route_core_pair_inventory_status ==
          :not_selected_legacy_source_pairs
    @test default_summary.route_core_pair_count == 0
    @test isempty(default_summary.route_core_pair_keys)
    @test !default_summary.route_core_pair_operator_ready
    @test default_summary.route_core_pair_operator_readiness_status ==
          :not_selected_legacy_source_pairs
    @test default_summary.route_core_pair_operator_blocker ==
          :not_selected_legacy_source_pairs
    @test default_summary.route_core_pair_operator_preflight_available
    @test default_summary.route_core_pair_operator_preflight_status ==
          :blocked_route_core_pair_operator_preflight
    @test default_summary.route_core_pair_operator_preflight_blocker ==
          :not_selected_legacy_source_pairs
    @test default_summary.route_core_pair_operator_preflight.operator_blocks_materialized ==
          false
    @test default_summary.route_core_pair_operator_plan_available
    @test default_summary.route_core_pair_operator_plan_status ==
          :blocked_route_core_pair_operator_plan
    @test default_summary.route_core_pair_operator_plan_blocker ==
          :not_selected_legacy_source_pairs
    @test !default_summary.route_core_pair_operator_plan.operator_blocks_materialized
    @test default_summary.pair_stage_fields_preserved
    @test length(default_pairs.pair_entries) ==
          length(default_stages.units.route_skeleton.pair_entries)
    @test default_pairs.pair_family_counts ==
          default_stages.units.route_skeleton.pair_family_counts
    @test default_pairs.helper_by_pair_family ==
          default_stages.units.route_skeleton.helper_by_pair_family
    @test default_pairs.pair_operator_helper_by_family ==
          default_stages.units.route_skeleton.helper_by_pair_family
    @test isnothing(default_pairs.pair_helper_status_by_family)
    @test !default_pairs.route_core_pair_inventory_available
    @test default_pairs.route_core_pair_inventory_status ==
          :not_selected_legacy_source_pairs
    @test !default_pairs.route_core_pair_operator_ready
    @test default_pairs.pair_stage == :pair_operator_terms_described

    atom_growth_stages = _cartesian_pair_stage_low_order_policy_pairs(
        fixture;
        policy = :atom_growth_complete_rectangular,
    )
    atom_growth_pairs = atom_growth_stages.pairs
    atom_growth_summary = atom_growth_pairs.low_order_pairs
    @test atom_growth_pairs.low_order_pair_route_kind ==
          :atom_growth_complete_rectangular_low_order_pairs
    @test atom_growth_pairs.atom_growth_pairs_selected
    @test !atom_growth_pairs.pair_operator_blocks_materialized
    @test !atom_growth_pairs.operator_pairs_materialized
    @test atom_growth_pairs.pair_inventory_source ==
          :atom_growth_unit_inventory
    @test atom_growth_pairs.independent_atom_growth_pair_inventory_available
    @test !atom_growth_pairs.active_source_authority
    @test atom_growth_summary.low_order_shellization_policy_resolved ==
          :atom_growth_complete_rectangular
    @test atom_growth_summary.shellization_source ==
          :bond_aligned_diatomic_atom_growth_construction_plan
    @test atom_growth_summary.shellization_kind ==
          :atom_growth_complete_rectangular
    @test atom_growth_summary.unit_route_kind ==
          :atom_growth_complete_rectangular_low_order_units
    @test atom_growth_summary.transform_route_kind ==
          :atom_growth_complete_rectangular_low_order_transforms
    @test atom_growth_summary.pair_route_kind ==
          :atom_growth_complete_rectangular_low_order_pairs
    @test atom_growth_summary.atom_growth_pairs_selected
    @test !atom_growth_summary.terminal_shellification_pairs_selected
    @test !atom_growth_summary.legacy_source_pairs_selected
    @test atom_growth_summary.plan_authority
    @test !atom_growth_summary.active_source_authority
    @test !atom_growth_summary.legacy_source_authority
    @test !atom_growth_summary.pair_operator_blocks_materialized
    @test !atom_growth_summary.operator_pairs_materialized
    @test atom_growth_summary.pair_inventory_known
    @test atom_growth_summary.pair_inventory_source ==
          :atom_growth_unit_inventory
    @test atom_growth_summary.independent_atom_growth_pair_inventory_available
    @test atom_growth_summary.route_core_pair_inventory_available
    @test atom_growth_summary.route_core_pair_inventory_status ==
          :available_route_core_unit_pair_inventory
    @test atom_growth_summary.route_core_pair_inventory isa
          GaussletBases.CartesianRouteCore.UnitPairInventory
    @test atom_growth_summary.route_core_pair_count == 36
    @test atom_growth_summary.route_core_pair_order_comparison_source ==
          :atom_growth_pair_inventory
    @test atom_growth_summary.route_core_pair_family_count_source ==
          :crc_final_unit_lowering_recipes
    @test atom_growth_summary.route_core_pair_operator_ready
    @test atom_growth_summary.route_core_pair_operator_readiness_status ==
          :ready_route_core_pair_operator_metadata
    @test isnothing(atom_growth_summary.route_core_pair_operator_blocker)
    @test atom_growth_summary.route_core_pair_operator_readiness_requirements ==
          (
              :complete_crc_final_unit_inventory,
              :available_crc_pair_inventory,
              :positive_crc_pair_count,
              :crc_pair_order_matches_staged,
              :crc_pair_family_metadata_available,
          )
    @test atom_growth_summary.route_core_pair_operator_preflight_available
    @test atom_growth_summary.route_core_pair_operator_preflight_status ==
          :ready_route_core_pair_operator_preflight
    @test isnothing(atom_growth_summary.route_core_pair_operator_preflight_blocker)
    atom_growth_preflight =
        atom_growth_summary.route_core_pair_operator_preflight
    @test atom_growth_preflight.route_core_final_unit_count == 8
    @test atom_growth_preflight.route_core_pair_count == 36
    @test atom_growth_preflight.route_core_pair_order_matches_staged
    @test !atom_growth_preflight.operator_blocks_materialized
    @test !atom_growth_preflight.pair_operator_blocks_materialized
    @test !atom_growth_preflight.source_operator_blocks_materialized
    @test atom_growth_summary.route_core_pair_operator_plan_available
    @test atom_growth_summary.route_core_pair_operator_plan_status ==
          :ready_route_core_pair_operator_plan
    @test isnothing(atom_growth_summary.route_core_pair_operator_plan_blocker)
    atom_growth_plan = atom_growth_summary.route_core_pair_operator_plan
    @test atom_growth_plan.route_core_final_unit_count == 8
    @test atom_growth_plan.route_core_pair_count == 36
    @test atom_growth_plan.route_core_pair_order_matches_staged
    @test !atom_growth_plan.source_operator_blocks_materialized
    @test !atom_growth_plan.pair_operator_blocks_materialized
    @test !atom_growth_plan.operator_blocks_materialized
    @test !atom_growth_plan.hamiltonian_matrices_materialized
    @test !isempty(atom_growth_plan.operator_block_family_plan)
    @test atom_growth_summary.route_core_typed_pair_operator_plan_inventory_available
    @test atom_growth_summary.route_core_typed_pair_operator_plan_inventory_status ==
          :blocked_route_core_pair_operator_plan_inventory
    @test atom_growth_summary.route_core_typed_pair_operator_plan_blocker ==
          :aggregate_subtree_operator_plan_required
    @test atom_growth_summary.route_core_typed_pair_operator_plan_count == 36
    @test atom_growth_summary.route_core_typed_pair_operator_plan_blocked_count == 15
    @test !atom_growth_summary.route_core_typed_pair_operator_plan_materialized
    @test _cartesian_pair_stage_count_by_field(
        atom_growth_summary.route_core_typed_pair_operator_source_path_counts,
        :source_operator_path,
        :aggregate_subtree_adapter_required,
    ) == 15
    @test _cartesian_pair_stage_count_by_field(
        atom_growth_summary.route_core_typed_pair_operator_materialization_status_counts,
        :materialization_status,
        :metadata_only_not_materialized,
    ) == 21
    @test _cartesian_pair_stage_count_by_field(
        atom_growth_summary.route_core_typed_pair_operator_materialization_status_counts,
        :materialization_status,
        :blocked_metadata_only_not_materialized,
    ) == 15
    @test _cartesian_pair_stage_count_by_field(
        atom_growth_summary.route_core_typed_pair_operator_blocker_counts,
        :blocker,
        nothing,
    ) == 21
    @test _cartesian_pair_stage_count_by_field(
        atom_growth_summary.route_core_typed_pair_operator_blocker_counts,
        :blocker,
        :aggregate_subtree_operator_plan_required,
    ) == 15
    @test sum(
        entry.pair_count for entry in
        atom_growth_summary.route_core_typed_pair_operator_plan_family_counts
    ) == 36
    @test all(
        entry -> !entry.materialized,
        atom_growth_summary.route_core_typed_pair_operator_plan_family_counts,
    )
    @test !atom_growth_summary.route_core_typed_pair_operator_materialization_ready
    @test atom_growth_summary.route_core_typed_pair_operator_materialization_readiness_status ==
          :blocked_pair_operator_materialization
    @test atom_growth_summary.route_core_typed_pair_operator_materialization_readiness_blocker ==
          :aggregate_subtree_operator_plan_required
    @test atom_growth_summary.route_core_typed_pair_operator_materialization_readiness_plan_count ==
          36
    @test atom_growth_summary.route_core_typed_pair_operator_materialization_readiness_blocked_count ==
          15
    @test atom_growth_summary.route_core_typed_pair_operator_materialization_readiness_materialized_count ==
          0
    @test atom_growth_pairs.route_core_pair_inventory_available
    @test atom_growth_pairs.route_core_pair_count ==
          atom_growth_summary.route_core_pair_count
    @test atom_growth_pairs.route_core_pair_order_matches_staged ==
          atom_growth_summary.route_core_pair_order_matches_staged
    @test atom_growth_pairs.route_core_pair_operator_ready
    @test atom_growth_pairs.route_core_pair_operator_preflight_status ==
          atom_growth_preflight.status
    @test atom_growth_pairs.route_core_pair_operator_plan_status ==
          atom_growth_plan.status
    @test atom_growth_pairs.route_core_typed_pair_operator_plan_count ==
          atom_growth_summary.route_core_typed_pair_operator_plan_count
    @test atom_growth_pairs.route_core_typed_pair_operator_plan_blocked_count ==
          atom_growth_summary.route_core_typed_pair_operator_plan_blocked_count
    @test atom_growth_pairs.route_core_typed_pair_operator_plan_materialized ==
          atom_growth_summary.route_core_typed_pair_operator_plan_materialized
    @test atom_growth_pairs.route_core_typed_pair_operator_materialization_ready ==
          atom_growth_summary.route_core_typed_pair_operator_materialization_ready
    @test atom_growth_pairs.route_core_typed_pair_operator_materialization_readiness_status ==
          atom_growth_summary.route_core_typed_pair_operator_materialization_readiness_status
    @test atom_growth_pairs.route_core_typed_pair_operator_materialization_readiness_blocker ==
          atom_growth_summary.route_core_typed_pair_operator_materialization_readiness_blocker
    pair_inventory = atom_growth_summary.pair_inventory
    unit_inventory = atom_growth_stages.units.plan_unit_inventory
    @test pair_inventory.object_kind == :cartesian_atom_growth_plan_pair_inventory
    @test pair_inventory.status == :available_atom_growth_pair_inventory
    @test pair_inventory.pair_inventory_source == :atom_growth_unit_inventory
    @test pair_inventory.unit_count == unit_inventory.unit_count
    @test pair_inventory.unit_count == 8
    @test pair_inventory.pair_count ==
          pair_inventory.unit_count * (pair_inventory.unit_count + 1) ÷ 2
    @test pair_inventory.pair_count == 36
    @test length(atom_growth_pairs.pair_entries) == pair_inventory.pair_count
    @test atom_growth_summary.pair_count == length(atom_growth_pairs.pair_entries)
    @test atom_growth_summary.pair_count == pair_inventory.pair_count
    @test atom_growth_summary.pair_family_counts ==
          pair_inventory.pair_family_counts
    @test pair_inventory.pair_family_counts.white_lindsey_low_order_atom_growth_unit_pair ==
          pair_inventory.pair_count
    @test all(
        pair -> pair.pair_family == :white_lindsey_low_order_atom_growth_unit_pair,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> pair.pair_contract == :planned_low_order_unit_pair_operator_block,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> pair.pair_inventory_source == :atom_growth_unit_inventory,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> pair.pair_planning_source == :final_retained_units,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> pair.left_final_retained_unit.object_kind ==
                :cartesian_final_retained_unit_contract,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> pair.right_final_retained_unit.object_kind ==
                :cartesian_final_retained_unit_contract,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> pair.left_final_retained_unit.downstream_of_lowering,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> pair.right_final_retained_unit.downstream_of_lowering,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> !pair.left_final_retained_unit.direct_shellification_region_alias,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> !pair.right_final_retained_unit.direct_shellification_region_alias,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> pair.left_owned_support.object_kind ==
                :cartesian_owned_support_region3d,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> pair.right_owned_support.object_kind ==
                :cartesian_owned_support_region3d,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> pair.left_unit_index <= pair.right_unit_index,
        pair_inventory.pair_entries,
    )
    expected_pair_keys = Tuple(
        (
            unit_inventory.plan_units[left_index].unit_key,
            unit_inventory.plan_units[right_index].unit_key,
        )
        for left_index in 1:unit_inventory.unit_count
        for right_index in left_index:unit_inventory.unit_count
    )
    @test Tuple(pair.pair_key for pair in pair_inventory.pair_entries) ==
          expected_pair_keys
    @test atom_growth_summary.route_core_pair_keys == expected_pair_keys
    @test atom_growth_summary.route_core_pair_order_matches_staged
    @test atom_growth_pairs.route_core_pair_keys == expected_pair_keys
    route_core_family_counts =
        atom_growth_summary.route_core_pair_family_counts
    @test !isempty(route_core_family_counts)
    @test all(
        record -> record.pair_family_source == :crc_final_unit_lowering_recipes,
        route_core_family_counts,
    )
    @test all(record -> length(record.pair_family) == 2, route_core_family_counts)
    @test sum(record.pair_count for record in route_core_family_counts) ==
          atom_growth_summary.route_core_pair_count
    @test all(
        pair -> !pair.operator_pair_block_materialized,
        pair_inventory.pair_entries,
    )
    @test all(pair -> !pair.operator_block_materialized, pair_inventory.pair_entries)
    @test !pair_inventory.operator_pairs_materialized
    @test !pair_inventory.pair_operator_blocks_materialized
    @test atom_growth_summary.route_skeleton_pair_entry_count ==
          length(atom_growth_stages.units.route_skeleton.pair_entries)
    @test atom_growth_summary.route_skeleton_pair_family_counts ==
          atom_growth_stages.units.route_skeleton.pair_family_counts
    @test atom_growth_summary.route_skeleton_pair_inventory_source ==
          :route_skeleton_compatibility_fields
    @test !atom_growth_summary.summary_only
    @test atom_growth_summary.pair_stage_fields_preserved
    @test length(atom_growth_pairs.route_skeleton_pair_entries) ==
          length(atom_growth_stages.units.route_skeleton.pair_entries)
    @test atom_growth_pairs.route_skeleton_pair_family_counts ==
          atom_growth_stages.units.route_skeleton.pair_family_counts
    expected_atom_growth_helper_status = (
        white_lindsey_low_order_atom_growth_unit_pair =
            :deferred_no_pair_operator_block_helper,
    )
    @test atom_growth_pairs.helper_by_pair_family ==
          expected_atom_growth_helper_status
    @test atom_growth_pairs.pair_operator_helper_by_family ==
          expected_atom_growth_helper_status
    @test atom_growth_pairs.pair_helper_status_by_family ==
          expected_atom_growth_helper_status
    @test !(:white_lindsey_low_order in keys(atom_growth_pairs.helper_by_pair_family))
    @test atom_growth_pairs.route_skeleton_helper_by_pair_family ==
          atom_growth_stages.units.route_skeleton.helper_by_pair_family
    @test atom_growth_pairs.pair_stage == :pair_operator_terms_described

    terminal_fixture =
        _cartesian_pair_stage_low_order_policy_fixture(
            probe_parent_axis_construction = false,
            atom_locations = ((-4.0, 0.0, 0.0), (4.0, 0.0, 0.0)),
            parent_axis_counts = (x = 13, y = 7, z = 7),
        )
    terminal_stages = _cartesian_pair_stage_low_order_policy_pairs(
        terminal_fixture;
        policy = :terminal_cartesian_shellification_geometry,
    )
    terminal_pairs = terminal_stages.pairs
    terminal_summary = terminal_pairs.low_order_pairs
    @test terminal_pairs.low_order_pair_route_kind ==
          :terminal_shellification_low_order_pairs
    @test terminal_pairs.terminal_shellification_pairs_selected
    @test terminal_pairs.terminal_shellification_pair_summary_available
    @test terminal_pairs.terminal_shellification_scaffold_available
    @test _cartesian_pair_stage_terminal_unit_inventory_fingerprint(
        terminal_pairs,
    ) == _cartesian_pair_stage_terminal_unit_inventory_fingerprint(
        terminal_stages.transforms,
    )
    _assert_terminal_lowering_contract_fields_match_pair_stage(
        terminal_pairs,
        terminal_stages.transforms,
    )
    @test !terminal_pairs.terminal_shellification_final_retained_unit_inventory_available
    @test !terminal_pairs.terminal_shellification_transform_contracts_available
    @test !terminal_pairs.terminal_shellification_pair_inventory_available
    @test terminal_pairs.terminal_shellification_pair_inventory_status ==
          :deferred_terminal_shellification_pair_inventory
    @test terminal_pairs.terminal_shellification_pair_materialization_status ==
          :deferred_terminal_shellification_pair_materialization
    @test !terminal_pairs.independent_atom_growth_pair_inventory_available
    @test !terminal_pairs.pair_operator_blocks_materialized
    @test !terminal_pairs.operator_pairs_materialized
    @test terminal_pairs.pair_inventory === nothing
    @test terminal_pairs.pair_inventory_source ==
          :terminal_shellification_scaffold
    @test terminal_pairs.pair_entries == ()
    @test terminal_pairs.pair_family_counts == ()
    @test terminal_pairs.helper_by_pair_family == ()
    @test terminal_pairs.pair_operator_helper_by_family == ()
    @test terminal_pairs.pair_helper_status_by_family == ()
    @test !terminal_pairs.route_core_pair_inventory_available
    @test terminal_pairs.route_core_pair_inventory_status ==
          :deferred_terminal_shellification_pair_inventory
    @test !terminal_pairs.route_core_pair_operator_ready
    @test terminal_pairs.route_core_pair_operator_readiness_status ==
          :deferred_terminal_shellification_pair_inventory
    @test terminal_pairs.route_core_pair_operator_blocker ==
          :deferred_terminal_shellification_pair_inventory
    @test terminal_pairs.route_core_typed_pair_operator_plan_inventory_status ==
          :deferred_terminal_shellification_typed_pair_operator_plan_inventory
    @test terminal_pairs.route_core_typed_pair_operator_plan_blocker ==
          :deferred_terminal_shellification_pair_inventory
    @test !terminal_pairs.route_core_typed_pair_operator_plan_inventory_available
    @test terminal_pairs.route_core_typed_pair_operator_plan_count == 0
    @test terminal_pairs.route_core_typed_pair_operator_plan_blocked_count == 0
    @test !terminal_pairs.route_core_typed_pair_operator_plan_materialized
    @test !terminal_pairs.route_core_typed_pair_operator_materialization_ready
    @test terminal_pairs.terminal_shellification_central_gap_region_count ==
          terminal_stages.transforms.terminal_shellification_central_gap_region_count
    @test terminal_pairs.terminal_shellification_central_midpoint_slab_count ==
          terminal_stages.transforms.terminal_shellification_central_midpoint_slab_count
    @test terminal_pairs.terminal_shellification_central_distorted_product_box_count ==
          terminal_stages.transforms.terminal_shellification_central_distorted_product_box_count
    @test terminal_pairs.terminal_shellification_central_gap_region_count == 3
    @test terminal_pairs.terminal_shellification_central_midpoint_slab_count == 3
    @test terminal_pairs.terminal_shellification_central_distorted_product_box_count ==
          0
    @test terminal_summary.object_kind == :cartesian_pair_stage_low_order_summary
    @test terminal_summary.status ==
          :deferred_terminal_shellification_pair_inventory
    @test terminal_summary.low_order_shellization_policy_resolved ==
          :terminal_cartesian_shellification_geometry
    @test terminal_summary.shellization_source ==
          :terminal_cartesian_shellification_geometry
    @test terminal_summary.shellization_kind ==
          :terminal_cartesian_shellification_geometry
    @test terminal_summary.unit_route_kind ==
          :terminal_shellification_low_order_units
    @test terminal_summary.transform_route_kind ==
          :terminal_shellification_low_order_transforms
    @test terminal_summary.pair_route_kind ==
          :terminal_shellification_low_order_pairs
    @test terminal_summary.terminal_shellification_pairs_selected
    @test !terminal_summary.atom_growth_pairs_selected
    @test !terminal_summary.legacy_source_pairs_selected
    @test terminal_summary.terminal_shellification_pair_summary_available
    @test terminal_summary.terminal_shellification_scaffold_available
    @test _cartesian_pair_stage_terminal_unit_inventory_fingerprint(
        terminal_summary,
    ) == _cartesian_pair_stage_terminal_unit_inventory_fingerprint(
        terminal_stages.transforms,
    )
    _assert_terminal_lowering_contract_fields_match_pair_stage(
        terminal_summary,
        terminal_stages.transforms,
    )
    @test !terminal_summary.terminal_shellification_final_retained_unit_inventory_available
    @test !terminal_summary.terminal_shellification_transform_contracts_available
    @test !terminal_summary.terminal_shellification_pair_inventory_available
    @test terminal_summary.terminal_shellification_pair_inventory_status ==
          :deferred_terminal_shellification_pair_inventory
    @test terminal_summary.terminal_shellification_pair_materialization_status ==
          :deferred_terminal_shellification_pair_materialization
    @test !terminal_summary.independent_atom_growth_pair_inventory_available
    @test !terminal_summary.pair_inventory_known
    @test terminal_summary.pair_inventory_source ==
          :terminal_shellification_scaffold
    @test terminal_summary.pair_inventory === nothing
    @test terminal_summary.pair_entries == ()
    @test terminal_summary.pair_count == 0
    @test terminal_summary.pair_family_counts == ()
    @test terminal_summary.helper_by_pair_family == ()
    @test terminal_summary.pair_operator_helper_by_family == ()
    @test terminal_summary.pair_helper_status_by_family == ()
    @test !terminal_summary.pair_operator_blocks_materialized
    @test !terminal_summary.operator_pairs_materialized
    @test terminal_summary.plan_authority
    @test !terminal_summary.active_source_authority
    @test !terminal_summary.legacy_source_authority
    @test terminal_summary.route_core_pair_inventory_status ==
          :deferred_terminal_shellification_pair_inventory
    @test terminal_summary.route_core_pair_inventory === nothing
    @test terminal_summary.route_core_pair_count == 0
    @test terminal_summary.route_core_pair_keys == ()
    @test terminal_summary.route_core_pair_family_counts == ()
    @test terminal_summary.route_core_summary_status ==
          :deferred_terminal_shellification_pair_inventory
    @test !terminal_summary.route_core_pair_operator_ready
    @test terminal_summary.route_core_pair_operator_preflight_available
    @test terminal_summary.route_core_pair_operator_preflight_status ==
          :blocked_route_core_pair_operator_preflight
    @test terminal_summary.route_core_pair_operator_preflight_blocker ==
          :deferred_terminal_shellification_pair_inventory
    @test !terminal_summary.route_core_pair_operator_preflight.operator_blocks_materialized
    @test terminal_summary.route_core_pair_operator_plan_available
    @test terminal_summary.route_core_pair_operator_plan_status ==
          :blocked_route_core_pair_operator_plan
    @test terminal_summary.route_core_pair_operator_plan_blocker ==
          :deferred_terminal_shellification_pair_inventory
    @test !terminal_summary.route_core_pair_operator_plan.operator_blocks_materialized
    @test !terminal_summary.route_core_typed_pair_operator_plan_inventory_available
    @test terminal_summary.route_core_typed_pair_operator_plan_inventory_status ==
          :deferred_terminal_shellification_typed_pair_operator_plan_inventory
    @test terminal_summary.route_core_typed_pair_operator_plan_blocker ==
          :deferred_terminal_shellification_pair_inventory
    @test length(terminal_summary.route_skeleton_pair_entries) ==
          length(terminal_stages.units.route_skeleton.pair_entries)
    @test terminal_summary.route_skeleton_pair_family_counts ==
          terminal_stages.units.route_skeleton.pair_family_counts
    @test terminal_summary.summary_only
    @test all(
        record -> !record.owned_support_is_cpb,
        terminal_summary.terminal_shellification_unit_inventory.terminal_region_units,
    )
    @test all(
        record -> !record.shellification_region_is_cpb,
        terminal_summary.terminal_shellification_unit_inventory.terminal_region_units,
    )
    @test terminal_summary.pair_stage_fields_preserved
    @test _cartesian_pair_stage_pair_metadata_fingerprint(terminal_summary) ==
          _cartesian_pair_stage_pair_metadata_fingerprint(terminal_pairs)
    @test length(terminal_pairs.route_skeleton_pair_entries) ==
          length(terminal_stages.units.route_skeleton.pair_entries)
    @test terminal_pairs.route_skeleton_pair_family_counts ==
          terminal_stages.units.route_skeleton.pair_family_counts
    @test terminal_pairs.route_skeleton_helper_by_pair_family ==
          terminal_stages.units.route_skeleton.helper_by_pair_family
    @test terminal_pairs.pair_stage == :pair_operator_terms_described

    blocked_units = merge(
        atom_growth_stages.units,
        (; plan_unit_inventory = nothing),
    )
    blocked_summary =
        GaussletBases._pqs_source_box_route_driver_pair_stage_low_order_summary(
            blocked_units,
            atom_growth_stages.transforms,
            atom_growth_stages.units.route_skeleton,
        )
    @test blocked_summary.atom_growth_pairs_selected
    @test !blocked_summary.independent_atom_growth_pair_inventory_available
    @test !blocked_summary.pair_inventory_known
    @test blocked_summary.pair_inventory_source ==
          :blocked_missing_plan_unit_inventory
    @test blocked_summary.pair_entries == ()
    @test blocked_summary.pair_count == 0
    @test blocked_summary.pair_family_counts.white_lindsey_low_order_atom_growth_unit_pair ==
          0
    @test blocked_summary.route_core_pair_inventory_available
    @test blocked_summary.route_core_pair_inventory_status ==
          :available_route_core_unit_pair_inventory
    @test blocked_summary.route_core_pair_count == 36
    @test blocked_summary.route_core_pair_order_matches_staged
    @test blocked_summary.route_core_pair_order_comparison_source ==
          :route_core_sidecar_staged_pair_keys
    @test blocked_summary.route_core_pair_operator_ready
    @test blocked_summary.helper_by_pair_family ==
          expected_atom_growth_helper_status
    @test length(blocked_summary.route_skeleton_pair_entries) ==
          length(atom_growth_stages.units.route_skeleton.pair_entries)
    @test blocked_summary.route_skeleton_pair_family_counts ==
          atom_growth_stages.units.route_skeleton.pair_family_counts
    @test blocked_summary.route_skeleton_helper_by_pair_family ==
          atom_growth_stages.units.route_skeleton.helper_by_pair_family

    missing_crc_units = merge(
        atom_growth_stages.units,
        (; route_core_sidecar_inventory = nothing),
    )
    missing_crc_summary =
        GaussletBases._pqs_source_box_route_driver_pair_stage_low_order_summary(
            missing_crc_units,
            atom_growth_stages.transforms,
            atom_growth_stages.units.route_skeleton,
        )
    @test missing_crc_summary.independent_atom_growth_pair_inventory_available
    @test !missing_crc_summary.route_core_pair_inventory_available
    @test !missing_crc_summary.route_core_pair_operator_ready
    @test missing_crc_summary.route_core_pair_operator_readiness_status ==
          :blocked_missing_route_core_sidecar_inventory
    @test missing_crc_summary.route_core_pair_operator_blocker ==
          :missing_route_core_sidecar_inventory
    @test missing_crc_summary.route_core_pair_operator_preflight_available
    @test missing_crc_summary.route_core_pair_operator_preflight_status ==
          :blocked_route_core_pair_operator_preflight
    @test missing_crc_summary.route_core_pair_operator_preflight_blocker ==
          :missing_route_core_sidecar_inventory
    @test !missing_crc_summary.route_core_pair_operator_preflight.operator_blocks_materialized
    @test missing_crc_summary.route_core_pair_operator_plan_available
    @test missing_crc_summary.route_core_pair_operator_plan_status ==
          :blocked_route_core_pair_operator_plan
    @test missing_crc_summary.route_core_pair_operator_plan_blocker ==
          :missing_route_core_sidecar_inventory
    @test !missing_crc_summary.route_core_pair_operator_plan.operator_blocks_materialized
end
