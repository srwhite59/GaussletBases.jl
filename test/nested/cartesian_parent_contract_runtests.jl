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
        probe_parent_axis_construction = false,
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
    be2 = _cartesian_parent_contract_parent(
        atom_symbols = ("Be", "Be"),
        nuclear_charges = (4, 4),
        atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        parent_axis_counts = (x = 9, y = 7, z = 9),
    )
    chain = _cartesian_parent_contract_parent(
        atom_symbols = ("Be", "Be", "Be"),
        nuclear_charges = (4, 4, 4),
        atom_locations = ((-4.0, 0.0, 0.0), (0.0, 0.0, 0.0), (4.0, 0.0, 0.0)),
        parent_axis_counts = (x = 11, y = 7, z = 7),
    )

    @test Tuple(propertynames(one_center)) == _CARTESIAN_PARENT_CONTRACT_FIELDS
    @test Tuple(propertynames(scalar_one_center)) == _CARTESIAN_PARENT_CONTRACT_FIELDS
    @test Tuple(propertynames(be2)) == _CARTESIAN_PARENT_CONTRACT_FIELDS
    @test Tuple(propertynames(chain)) == _CARTESIAN_PARENT_CONTRACT_FIELDS

    @test one_center.object_kind == :cartesian_route_parent
    @test scalar_one_center.object_kind == :cartesian_route_parent
    @test be2.object_kind == :cartesian_route_parent
    @test chain.object_kind == :cartesian_route_parent

    @test one_center.axis_counts == (x = 7, y = 7, z = 7)
    @test be2.axis_counts == (x = 9, y = 7, z = 9)
    @test chain.axis_counts == (x = 11, y = 7, z = 7)
    @test one_center.axis_counts_source == :manual_fixture
    @test be2.axis_counts_source == :manual_fixture
    @test chain.axis_counts_source == :manual_fixture
    @test one_center.axis_counts_status == :available
    @test be2.axis_counts_status == :available
    @test chain.axis_counts_status == :available

    @test one_center.atom_count == 1
    @test scalar_one_center.atom_count == 1
    @test scalar_one_center.atom_symbols == ("Be",)
    @test scalar_one_center.nuclear_charges == (4.0,)
    @test be2.atom_count == 2
    @test chain.atom_count == 3
    @test one_center.center_count == one_center.atom_count
    @test be2.center_count == be2.atom_count
    @test chain.center_count == chain.atom_count
    @test Tuple(center.center_index for center in be2.center_table) == (1, 2)
    @test Tuple(center.atom_symbol for center in be2.center_table) == ("Be", "Be")
    @test Tuple(center.nuclear_charge for center in be2.center_table) == (4, 4)

    @test one_center.system_classification == :one_center
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

    @test chain.system_classification == :axis_aligned_chain_metadata_only
    @test chain.system_classification_status ==
          :multi_center_single_axis_chain_not_materialized
    @test chain.bond_axis === nothing
    @test chain.chain_axis == :x
    @test chain.center_axis_metadata.axis_aligned
    @test chain.center_axis_metadata.center_axis_coordinates == (-4.0, 0.0, 4.0)

    @test one_center.parent_axis_readiness.parent_axis_counts_status == :manual_fixture
    @test be2.parent_axis_readiness.parent_axis_counts_status == :manual_fixture
    @test chain.parent_axis_readiness.parent_axis_counts_status == :manual_fixture
    @test one_center.parent_basis_materialization_status ==
          :metadata_only_not_materialized
    @test be2.parent_basis_materialization_status ==
          :metadata_only_not_materialized
    @test chain.parent_basis_materialization_status ==
          :metadata_only_not_materialized
    @test !one_center.parent_basis_materialized
    @test !be2.parent_basis_materialized
    @test !chain.parent_basis_materialized
    @test !one_center.axis_bundle_materialized
    @test !be2.axis_bundle_materialized
    @test !chain.axis_bundle_materialized
end
