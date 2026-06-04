using Test
using GaussletBases

const _CARTESIAN_PARENT_CONTRACT_FIELDS = (
    :object_kind,
    :status,
    :system,
    :spacing_inputs,
    :parent_inputs,
    :standard_setup,
    :parent_axis,
    :parent_axis_readiness,
    :parent_axis_probe,
    :route_axis_counts,
    :atom_count,
    :atom_symbols,
    :nuclear_charges,
    :atom_locations,
    :center_table,
    :center_count,
    :center_axis_metadata,
    :system_classification,
    :system_classification_status,
    :bond_axis,
    :chain_axis,
    :axis_counts,
    :axis_counts_source,
    :axis_counts_status,
    :physical_box,
    :physical_box_rule,
    :parent_object_carry,
    :parent_basis_object,
    :parent_qw_basis_object,
    :parent_axis_bundle_object,
    :parent_basis_object_available,
    :parent_qw_basis_object_available,
    :parent_axis_bundle_object_available,
    :parent_basis_object_type_label,
    :parent_qw_basis_object_type_label,
    :parent_axis_bundle_object_type_label,
    :parent_materialization_plan,
    :parent_materialization_plan_status,
    :parent_materialization_planning_family,
    :parent_materialization_blocker,
    :parent_basis_materialization,
    :parent_basis_materialization_status,
    :parent_basis_materialized,
    :parent_axis_metadata_constructed,
    :axis_bundle_materialized,
)

function _cartesian_parent_contract_recipe()
    return GaussletBases.cartesian_recipe(
        (;
            route_family = :pqs_source_box,
            route_kind = :cartesian_parent_contract_probe,
            route_shape = (:pqs_left, :product, :pqs_right),
            product_body_rule = :centered_single_z_slab,
            pqs_retained_rule = :boundary_comx_product_mode_selection,
            product_retained_rule = :product_doside_retained_unit,
            terms = (:overlap, :kinetic),
            pair_factor_normalization = :density_normalized,
            support_dense_direct_allowed = false,
            reference_only_authorities = (:support_row_oracle,),
            white_lindsey_route_shape =
                (:standard_cartesian_units, :low_order_comx_coarsening),
            white_lindsey_mapping_rule = :standard_unit_backbone_mapping_family,
            white_lindsey_nesting_rule = :unit_box_low_order_comx_coarsening,
            white_lindsey_retained_rule = :low_order_unit_comx_retained_basis,
            white_lindsey_operator_rule = :low_order_unit_operator_blocks,
            white_lindsey_benchmark_role =
                :published_cartesian_baseline_for_pqs_comparison,
        ),
    )
end

function _cartesian_parent_contract_parent(;
    atom_symbols,
    nuclear_charges,
    atom_locations,
    parent_axis_counts,
    probe_parent_axis_construction = false,
)
    system = GaussletBases.cartesian_system(
        (;
            atom_symbols,
            nuclear_charges,
            atom_locations,
            radius = 15.0,
            parent_axis_counts,
            map_backend = :pgdg_localized_experimental,
        ),
    )
    spacing_inputs = (;
        q = 5,
        n_s = 5,
        reference_spacing = 1.0,
        tail_spacing = 10.0,
        q_to_core_spacing_rule = :standard_pqs_ns_equals_q,
        core_spacing = nothing,
    )
    parent_inputs = (;
        probe_parent_axis_construction,
        parent_axis_probe_backend = :pgdg_localized_experimental,
        parent_axis_probe_family = :G10,
    )
    return GaussletBases.cartesian_parent(
        system,
        spacing_inputs,
        parent_inputs,
        _cartesian_parent_contract_recipe(),
    )
end

@testset "Cartesian parent shared contract" begin
    one_center = _cartesian_parent_contract_parent(
        atom_symbols = ("Be",),
        nuclear_charges = (4,),
        atom_locations = ((0.0, 0.0, 0.0),),
        parent_axis_counts = (x = 7, y = 7, z = 7),
    )
    scalar_one_center = _cartesian_parent_contract_parent(
        atom_symbols = "Be",
        nuclear_charges = 4,
        atom_locations = ((0.0, 0.0, 0.0),),
        parent_axis_counts = (x = 7, y = 7, z = 7),
    )
    one_center_rectangular = _cartesian_parent_contract_parent(
        atom_symbols = ("Be",),
        nuclear_charges = (4,),
        atom_locations = ((0.0, 0.0, 0.0),),
        parent_axis_counts = (x = 7, y = 5, z = 9),
    )
    one_center_shifted = _cartesian_parent_contract_parent(
        atom_symbols = ("Be",),
        nuclear_charges = (4,),
        atom_locations = ((1.0, 0.0, 0.0),),
        parent_axis_counts = (x = 7, y = 7, z = 7),
    )
    be2 = _cartesian_parent_contract_parent(
        atom_symbols = ("Be", "Be"),
        nuclear_charges = (4, 4),
        atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        parent_axis_counts = (x = 9, y = 7, z = 9),
    )
    probed_be2 = _cartesian_parent_contract_parent(
        atom_symbols = ("Be", "Be"),
        nuclear_charges = (4, 4),
        atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        parent_axis_counts = nothing,
        probe_parent_axis_construction = true,
    )
    chain = _cartesian_parent_contract_parent(
        atom_symbols = ("Be", "Be", "Be"),
        nuclear_charges = (4, 4, 4),
        atom_locations = ((-4.0, 0.0, 0.0), (0.0, 0.0, 0.0), (4.0, 0.0, 0.0)),
        parent_axis_counts = (x = 11, y = 7, z = 7),
    )
    heteronuclear_chain = _cartesian_parent_contract_parent(
        atom_symbols = ("Li", "Be", "B"),
        nuclear_charges = (3, 4, 5),
        atom_locations = ((-4.0, 0.0, 0.0), (0.0, 0.0, 0.0), (4.0, 0.0, 0.0)),
        parent_axis_counts = (x = 11, y = 7, z = 7),
    )
    non_axis_aligned = _cartesian_parent_contract_parent(
        atom_symbols = ("Be", "Be"),
        nuclear_charges = (4, 4),
        atom_locations = ((-2.0, 0.0, 0.0), (2.0, 1.0, 0.0)),
        parent_axis_counts = (x = 9, y = 7, z = 9),
    )

    @test Tuple(propertynames(one_center)) == _CARTESIAN_PARENT_CONTRACT_FIELDS
    @test Tuple(propertynames(scalar_one_center)) == _CARTESIAN_PARENT_CONTRACT_FIELDS
    @test Tuple(propertynames(one_center_rectangular)) ==
          _CARTESIAN_PARENT_CONTRACT_FIELDS
    @test Tuple(propertynames(one_center_shifted)) ==
          _CARTESIAN_PARENT_CONTRACT_FIELDS
    @test Tuple(propertynames(be2)) == _CARTESIAN_PARENT_CONTRACT_FIELDS
    @test Tuple(propertynames(probed_be2)) == _CARTESIAN_PARENT_CONTRACT_FIELDS
    @test Tuple(propertynames(chain)) == _CARTESIAN_PARENT_CONTRACT_FIELDS
    @test Tuple(propertynames(heteronuclear_chain)) == _CARTESIAN_PARENT_CONTRACT_FIELDS
    @test Tuple(propertynames(non_axis_aligned)) == _CARTESIAN_PARENT_CONTRACT_FIELDS

    @test one_center.object_kind == :cartesian_route_parent
    @test scalar_one_center.object_kind == :cartesian_route_parent
    @test one_center_rectangular.object_kind == :cartesian_route_parent
    @test one_center_shifted.object_kind == :cartesian_route_parent
    @test be2.object_kind == :cartesian_route_parent
    @test probed_be2.object_kind == :cartesian_route_parent
    @test chain.object_kind == :cartesian_route_parent
    @test heteronuclear_chain.object_kind == :cartesian_route_parent
    @test non_axis_aligned.object_kind == :cartesian_route_parent

    @test one_center.axis_counts == (x = 7, y = 7, z = 7)
    @test one_center_rectangular.axis_counts == (x = 7, y = 5, z = 9)
    @test one_center_shifted.axis_counts == (x = 7, y = 7, z = 7)
    @test be2.axis_counts == (x = 9, y = 7, z = 9)
    @test probed_be2.axis_counts == (x = 31, y = 17, z = 17)
    @test chain.axis_counts == (x = 11, y = 7, z = 7)
    @test heteronuclear_chain.axis_counts == (x = 11, y = 7, z = 7)
    @test non_axis_aligned.axis_counts == (x = 9, y = 7, z = 9)
    @test one_center.axis_counts_source == :manual_fixture
    @test one_center_rectangular.axis_counts_source == :manual_fixture
    @test one_center_shifted.axis_counts_source == :manual_fixture
    @test be2.axis_counts_source == :manual_fixture
    @test probed_be2.axis_counts_source == :constructed_parent_axis_probe
    @test chain.axis_counts_source == :manual_fixture
    @test heteronuclear_chain.axis_counts_source == :manual_fixture
    @test non_axis_aligned.axis_counts_source == :manual_fixture
    @test one_center.axis_counts_status == :available
    @test one_center_rectangular.axis_counts_status == :available
    @test one_center_shifted.axis_counts_status == :available
    @test be2.axis_counts_status == :available
    @test probed_be2.axis_counts_status == :available
    @test chain.axis_counts_status == :available
    @test heteronuclear_chain.axis_counts_status == :available
    @test non_axis_aligned.axis_counts_status == :available

    @test one_center.atom_count == 1
    @test scalar_one_center.atom_count == 1
    @test one_center_rectangular.atom_count == 1
    @test one_center_shifted.atom_count == 1
    @test scalar_one_center.atom_symbols == ("Be",)
    @test scalar_one_center.nuclear_charges == (4.0,)
    @test be2.atom_count == 2
    @test probed_be2.atom_count == 2
    @test chain.atom_count == 3
    @test heteronuclear_chain.atom_count == 3
    @test non_axis_aligned.atom_count == 2
    @test one_center.center_count == one_center.atom_count
    @test one_center_rectangular.center_count == one_center_rectangular.atom_count
    @test one_center_shifted.center_count == one_center_shifted.atom_count
    @test be2.center_count == be2.atom_count
    @test probed_be2.center_count == probed_be2.atom_count
    @test chain.center_count == chain.atom_count
    @test heteronuclear_chain.center_count == heteronuclear_chain.atom_count
    @test non_axis_aligned.center_count == non_axis_aligned.atom_count
    @test Tuple(center.center_index for center in be2.center_table) == (1, 2)
    @test Tuple(center.atom_symbol for center in be2.center_table) == ("Be", "Be")
    @test Tuple(center.nuclear_charge for center in be2.center_table) == (4, 4)

    @test one_center.system_classification == :one_center
    @test one_center_rectangular.system_classification == :one_center
    @test one_center_shifted.system_classification == :one_center
    @test one_center.system_classification_status == :explicit_atom_count_one
    @test one_center.bond_axis === nothing
    @test one_center.chain_axis === nothing
    @test isempty(one_center.center_axis_metadata.active_axes)

    @test be2.system_classification == :bond_aligned_diatomic
    @test be2.system_classification_status ==
          :explicit_two_atom_single_axis_separation
    @test be2.bond_axis == :x
    @test be2.chain_axis == :x
    @test be2.center_axis_metadata.axis_aligned
    @test be2.center_axis_metadata.center_axis_coordinates == (-2.0, 2.0)
    @test probed_be2.system_classification == be2.system_classification
    @test probed_be2.bond_axis == :x
    @test probed_be2.chain_axis == :x

    @test chain.system_classification == :axis_aligned_chain_metadata_only
    @test chain.system_classification_status ==
          :multi_center_single_axis_chain_not_materialized
    @test chain.bond_axis === nothing
    @test chain.chain_axis == :x
    @test chain.center_axis_metadata.axis_aligned
    @test chain.center_axis_metadata.center_axis_coordinates == (-4.0, 0.0, 4.0)

    @test heteronuclear_chain.system_classification == :axis_aligned_chain_metadata_only
    @test heteronuclear_chain.chain_axis == :x
    @test heteronuclear_chain.center_axis_metadata.axis_aligned

    @test non_axis_aligned.system_classification == :pending_system_classification
    @test non_axis_aligned.system_classification_status ==
          :diatomic_not_axis_aligned_by_metadata
    @test non_axis_aligned.bond_axis === nothing
    @test non_axis_aligned.chain_axis === nothing
    @test !non_axis_aligned.center_axis_metadata.axis_aligned
    @test non_axis_aligned.center_axis_metadata.active_axes == (:x, :y)

    @test one_center.parent_axis_readiness.parent_axis_counts_status == :manual_fixture
    @test one_center_rectangular.parent_axis_readiness.parent_axis_counts_status ==
          :manual_fixture
    @test one_center_shifted.parent_axis_readiness.parent_axis_counts_status ==
          :manual_fixture
    @test be2.parent_axis_readiness.parent_axis_counts_status == :manual_fixture
    @test probed_be2.parent_axis_readiness.parent_axis_counts_status ==
          :pending_helper_or_documented_rule
    @test chain.parent_axis_readiness.parent_axis_counts_status == :manual_fixture
    @test heteronuclear_chain.parent_axis_readiness.parent_axis_counts_status ==
          :manual_fixture
    @test non_axis_aligned.parent_axis_readiness.parent_axis_counts_status ==
          :manual_fixture

    @test one_center.parent_materialization_plan.object_kind ==
          :cartesian_parent_materialization_plan
    @test one_center.parent_materialization_plan_status ==
          :materialized_parent_objects_available
    @test one_center.parent_materialization_planning_family ==
          :one_center_parent_lattice
    @test one_center.parent_materialization_blocker === nothing
    @test one_center.parent_materialization_plan.one_center_compatible
    @test one_center.parent_materialization_plan.constructs_basis_now
    @test one_center.parent_materialization_plan.constructs_axis_bundle_now
    @test one_center.parent_basis_object_available
    @test !one_center.parent_qw_basis_object_available
    @test one_center.parent_axis_bundle_object_available
    @test one_center.parent_basis_object !== nothing
    @test one_center.parent_qw_basis_object === nothing
    @test one_center.parent_axis_bundle_object !== nothing
    @test one_center.parent_object_carry.parent_axis_bundle_object_source ==
          :one_center_mapped_ordinary_axis_bundles
    @test one_center.parent_basis_object_type_label ==
          "CartesianParentGaussletBasis3D"
    @test one_center.parent_qw_basis_object_type_label == "unavailable"
    @test one_center.parent_axis_bundle_object_type_label ==
          "_CartesianNestedAxisBundles3D"
    @test GaussletBases.CartesianParentGaussletBases.parent_axis_counts(
        one_center.parent_basis_object,
    ) == (7, 7, 7)
    @test GaussletBases.CartesianParentGaussletBases.parent_axis_counts(
        one_center_rectangular.parent_basis_object,
    ) == (7, 5, 9)
    @test one_center_rectangular.parent_basis_object.axis_sharing == :separate_axes
    @test one_center_rectangular.parent_axis_bundle_object_available
    @test GaussletBases._nested_axis_lengths(
        one_center_rectangular.parent_axis_bundle_object,
    ) == (7, 5, 9)

    @test one_center_shifted.parent_materialization_plan_status ==
          :metadata_only_pending_one_center_parent_axis_builder
    @test one_center_shifted.parent_materialization_blocker ==
          :pending_one_center_parent_axis_builder
    @test !one_center_shifted.parent_materialization_plan.constructs_basis_now
    @test !one_center_shifted.parent_materialization_plan.constructs_axis_bundle_now
    @test !one_center_shifted.parent_basis_object_available
    @test !one_center_shifted.parent_qw_basis_object_available
    @test !one_center_shifted.parent_axis_bundle_object_available
    @test one_center_shifted.parent_basis_object === nothing
    @test one_center_shifted.parent_object_carry.parent_basis_object_source ==
          :pending_non_origin_one_center_parent_mapping

    @test be2.parent_materialization_plan_status ==
          :metadata_only_diatomic_parent_api_candidate
    @test be2.parent_materialization_planning_family ==
          :bond_aligned_diatomic_parent_lattice
    @test be2.parent_materialization_plan.intended_parent_constructor ==
          :bond_aligned_homonuclear_qw_basis
    @test be2.parent_materialization_plan.intended_axis_bundle_helper ==
          :_qwrg_bond_aligned_axis_bundles
    @test be2.parent_materialization_plan.bond_aligned_diatomic_compatible
    @test !be2.parent_materialization_plan.constructs_basis_now
    @test !be2.parent_materialization_plan.constructs_axis_bundle_now
    @test !be2.parent_basis_object_available
    @test !be2.parent_axis_bundle_object_available
    @test be2.parent_basis_object === nothing
    @test be2.parent_axis_bundle_object === nothing

    @test probed_be2.parent_axis_probe.carry_objects_requested
    @test probed_be2.parent_axis_probe.basis_object_available
    @test probed_be2.parent_axis_probe.axis_bundle_object_available
    @test probed_be2.parent_basis_object_available
    @test probed_be2.parent_qw_basis_object_available
    @test probed_be2.parent_axis_bundle_object_available
    @test probed_be2.parent_basis_object !== nothing
    @test probed_be2.parent_qw_basis_object !== nothing
    @test probed_be2.parent_axis_bundle_object !== nothing
    @test probed_be2.parent_basis_object_type_label ==
          "CartesianParentGaussletBasis3D"
    @test probed_be2.parent_qw_basis_object_type_label ==
          "BondAlignedDiatomicQWBasis3D"
    @test probed_be2.parent_axis_bundle_object_type_label ==
          "_CartesianNestedAxisBundles3D"
    @test probed_be2.parent_materialization_plan_status ==
          :materialized_parent_objects_available
    @test probed_be2.parent_materialization_blocker === nothing
    @test probed_be2.parent_materialization_plan.constructs_basis_now
    @test probed_be2.parent_materialization_plan.constructs_axis_bundle_now

    @test chain.parent_materialization_plan_status ==
          :metadata_only_chain_parent_constructor_candidate
    @test chain.parent_materialization_planning_family ==
          :axis_aligned_chain_parent_lattice
    @test chain.parent_materialization_plan.intended_parent_constructor ==
          :bond_aligned_homonuclear_chain_qw_basis
    @test chain.parent_materialization_plan.axis_aligned_chain_compatible
    @test chain.parent_materialization_blocker ==
          :chain_parent_materializer_not_connected

    @test heteronuclear_chain.parent_materialization_plan_status ==
          :blocked_unsupported_heteronuclear_chain_parent_materializer
    @test heteronuclear_chain.parent_materialization_planning_family ==
          :axis_aligned_chain_parent_lattice
    @test heteronuclear_chain.parent_materialization_plan.intended_parent_constructor ===
          nothing
    @test heteronuclear_chain.parent_materialization_plan.axis_aligned_chain_compatible
    @test heteronuclear_chain.parent_materialization_blocker ==
          :heteronuclear_chain_parent_materializer_not_available

    @test non_axis_aligned.parent_materialization_plan_status ==
          :blocked_pending_system_classification
    @test non_axis_aligned.parent_materialization_planning_family ==
          :pending_system_classification_parent_lattice
    @test non_axis_aligned.parent_materialization_blocker ==
          :pending_system_classification
    @test non_axis_aligned.parent_materialization_plan.blocked

    @test one_center.parent_basis_materialization_status ==
          :materialized_parent_objects_available
    @test one_center_shifted.parent_basis_materialization_status ==
          :metadata_only_not_materialized
    @test be2.parent_basis_materialization_status ==
          :metadata_only_not_materialized
    @test probed_be2.parent_basis_materialization_status ==
          :materialized_parent_objects_available
    @test chain.parent_basis_materialization_status ==
          :metadata_only_not_materialized
    @test heteronuclear_chain.parent_basis_materialization_status ==
          :metadata_only_not_materialized
    @test non_axis_aligned.parent_basis_materialization_status ==
          :metadata_only_not_materialized
    @test one_center.parent_basis_materialized
    @test !one_center_shifted.parent_basis_materialized
    @test !be2.parent_basis_materialized
    @test probed_be2.parent_basis_materialized
    @test !chain.parent_basis_materialized
    @test !heteronuclear_chain.parent_basis_materialized
    @test !non_axis_aligned.parent_basis_materialized
    @test one_center.axis_bundle_materialized
    @test !one_center_shifted.axis_bundle_materialized
    @test !be2.axis_bundle_materialized
    @test probed_be2.axis_bundle_materialized
    @test !chain.axis_bundle_materialized
    @test !heteronuclear_chain.axis_bundle_materialized
    @test !non_axis_aligned.axis_bundle_materialized
end
