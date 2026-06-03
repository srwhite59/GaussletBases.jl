using Test
using GaussletBases

@testset "PQS source-box route skeleton helper" begin
    metrics_module = GaussletBases.CartesianContractedParentMetrics

    skeleton =
        metrics_module._pqs_pqs_product_source_box_route_skeleton(
            q = 5,
            parent_axis_counts = (x = 9, y = 7, z = 9),
            route_shape = (:pqs_left, :product, :pqs_right),
            product_body_rule = :centered_single_z_slab,
            pair_factor_normalization = :density_normalized,
        )

    @test skeleton.object_kind == :pqs_pqs_product_source_box_route_skeleton
    @test skeleton.status == :private_development_skeleton
    @test skeleton.route_shape == (:pqs_left, :product, :pqs_right)
    @test skeleton.retained_unit_order == (:pqs_left, :pqs_right, :product)
    @test skeleton.source_boxes.pqs_left == (x = 1:5, y = 1:5, z = 1:5)
    @test skeleton.source_boxes.pqs_right == (x = 1:5, y = 1:5, z = 5:9)
    @test skeleton.source_boxes.product == (x = 1:5, y = 1:5, z = 5:5)
    @test skeleton.source_dimensions ==
          (pqs_left = (5, 5, 5), pqs_right = (5, 5, 5), product = (5, 5, 1))
    @test skeleton.retained_counts == (pqs_left = 98, pqs_right = 98, product = 25)
    @test skeleton.ranges == (pqs_left = 1:98, pqs_right = 99:196, product = 197:221)
    @test skeleton.retained_dimension == 221
    @test length(skeleton.pair_entries) == 6
    @test skeleton.pair_family_counts ==
          (pqs_pqs = 3, pqs_product = 2, product_pqs = 0, product_product = 1)
    @test skeleton.helper_by_pair_family.pqs_pqs ==
          :_pqs_pqs_source_box_density_density_interaction_block
    @test skeleton.helper_by_pair_family.pqs_product ==
          :_pqs_product_source_box_density_density_interaction_block
    @test skeleton.helper_by_pair_family.product_product ==
          :_product_doside_source_box_density_density_interaction_block
    @test skeleton.product_body.length == 1
    @test skeleton.product_body.derivation ==
          :derived_from_centered_single_z_slab_rule
    @test skeleton.diagnostics.pqs_retained_count_derivation ==
          :boundary_comx_product_mode_selection
    @test skeleton.diagnostics.product_retained_count_derivation ==
          :product_source_dimension_under_product_doside_rule
    @test skeleton.diagnostics.derived_retained_counts
    @test skeleton.diagnostics.derived_pair_inventory
    @test skeleton.diagnostics.source_box_first
    @test skeleton.diagnostics.source_box_algorithmic_path_true_for_every_pair
    @test !skeleton.diagnostics.packet_adoption
    @test !skeleton.diagnostics.public_default_consumes
    @test !skeleton.diagnostics.retained_weight_division_allowed
    @test !skeleton.diagnostics.repo_side_ray_id

    raw_skeleton =
        metrics_module._pqs_pqs_product_source_box_route_skeleton(
            q = 5,
            parent_axis_counts = (9, 7, 9),
            pair_factor_normalization = :raw_weighted,
        )
    @test raw_skeleton.helper_by_pair_family.pqs_pqs ==
          :_pqs_pqs_source_box_raw_weighted_density_density_interaction_block
    @test raw_skeleton.helper_by_pair_family.pqs_product ==
          :_pqs_product_source_box_raw_weighted_density_density_interaction_block
    @test raw_skeleton.helper_by_pair_family.product_product ==
          :_product_doside_source_box_raw_weighted_density_density_interaction_block
    @test raw_skeleton.diagnostics.raw_weight_division_owner ==
          :explicit_source_quadrature_weight_outer_products
    @test raw_skeleton.retained_counts == skeleton.retained_counts
    @test raw_skeleton.ranges == skeleton.ranges

    manual_length_skeleton =
        metrics_module._pqs_pqs_product_source_box_route_skeleton(
            q = 5,
            parent_axis_counts = (x = 9, y = 7, z = 11),
            product_body_rule = (kind = :centered_z_slab, length = 3),
        )
    @test manual_length_skeleton.source_boxes.product == (x = 1:5, y = 1:5, z = 5:7)
    @test manual_length_skeleton.retained_counts.product == 75
    @test manual_length_skeleton.product_body.derivation ==
          :manual_fixture_named_rule_length

    @test_throws ArgumentError metrics_module._pqs_pqs_product_source_box_route_skeleton(
        q = 1,
        parent_axis_counts = (x = 9, y = 7, z = 9),
    )
    @test_throws ArgumentError metrics_module._pqs_pqs_product_source_box_route_skeleton(
        q = 5,
        parent_axis_counts = (x = 4, y = 7, z = 9),
    )
end
