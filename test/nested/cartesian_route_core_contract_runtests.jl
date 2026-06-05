using Test
using GaussletBases

const CRC = GaussletBases.CartesianRouteCore

@testset "CartesianRouteCore CPB contract" begin
    filled = CRC.filled_cpb(1:5, 1:5, 1:5; role = :filled_source_cpb)
    @test CRC.shape(filled) == (5, 5, 5)
    @test CRC.codimension(filled) == 0
    @test CRC.support_count(filled) == 125

    facet = CRC.cpb(1:1, 2:4, 2:4; role = :x_low_facet)
    edge = CRC.cpb(1:1, 1:1, 2:4; role = :x_low_y_low_edge)
    corner = CRC.cpb(1:1, 1:1, 1:1; role = :x_low_y_low_z_low_corner)
    @test CRC.codimension(facet) == 1
    @test CRC.codimension(edge) == 2
    @test CRC.codimension(corner) == 3

    outer = CRC.filled_cpb(1:5, 1:5, 1:5; role = :outer_box)
    inner = CRC.filled_cpb(2:4, 2:4, 2:4; role = :inner_exclusion_box)
    shell_support = CRC.complete_shell_support(outer, inner)
    @test shell_support isa CRC.OwnedSupport
    @test !(shell_support isa CRC.CoordinateProductBox)
    @test CRC.support_count(shell_support) == 98

    strata = CRC.complete_shell_boundary_strata(outer, inner)
    @test length(strata.facets) == 6
    @test length(strata.edges) == 12
    @test length(strata.corners) == 8
    @test length(strata.all_strata) == 26
    @test strata.stratum_support_count == CRC.support_count(shell_support)
    @test all(cpb -> CRC.codimension(cpb) == 1, strata.facets)
    @test all(cpb -> CRC.codimension(cpb) == 2, strata.edges)
    @test all(cpb -> CRC.codimension(cpb) == 3, strata.corners)

    region = CRC.shellification_region(
        :regular_shared_molecular_shell,
        shell_support;
        metadata = (; policy = :atom_growth),
    )

    lw_lowering = CRC.white_lindsey_boundary_strata_lowering(region, strata.all_strata)
    @test CRC.lowering_recipe(lw_lowering) == :white_lindsey_boundary_strata
    @test length(CRC.source_cpbs(lw_lowering)) == 26

    lw_intermediate = CRC.intermediate_retained_space(
        lw_lowering;
        retained_rule = :white_lindsey_boundary_stratum_product,
        dimension = 98,
        materialized = false,
    )
    lw_realization = CRC.trivial_shell_realization(lw_intermediate, region)
    lw_unit = CRC.final_retained_unit(
        :lw_shell,
        :regular_shared_molecular_shell,
        lw_lowering,
        lw_intermediate,
        lw_realization;
        dimension = 98,
    )
    @test lw_unit.dimension == 98

    pqs_lowering = CRC.pqs_filled_source_lowering(region, filled)
    @test CRC.lowering_recipe(pqs_lowering) == :pqs_filled_source_cpb
    @test only(CRC.source_cpbs(pqs_lowering)) === filled

    pqs_intermediate = CRC.intermediate_retained_space(
        pqs_lowering;
        retained_rule = :pqs_boundary_comx_product_modes,
        source_mode_dims = (5, 5, 5),
    )
    @test pqs_intermediate.dimension == 98
    pqs_realization = CRC.pqs_shell_realization(pqs_intermediate, region)
    @test pqs_realization.realization_kind == :shell_projection_lowdin
    pqs_unit = CRC.final_retained_unit(
        :pqs_shell,
        :regular_shared_molecular_shell,
        pqs_lowering,
        pqs_intermediate,
        pqs_realization;
        dimension = 98,
    )

    direct_support = CRC.owned_cpb(CRC.slab_cpb(1:5, 1:5, 3:3; role = :midpoint_slab))
    direct_region = CRC.shellification_region(:midpoint_slab, direct_support)
    direct_lowering = CRC.lowering_source(
        :direct_identity_cpb,
        direct_region,
        direct_support.cpbs,
    )
    direct_intermediate = CRC.intermediate_retained_space(
        direct_lowering;
        retained_rule = :identity_source_modes,
        source_mode_dims = (5, 5, 1),
    )
    direct_realization = CRC.trivial_shell_realization(direct_intermediate, direct_region)
    direct_unit = CRC.final_retained_unit(
        :midpoint,
        :midpoint_slab,
        direct_lowering,
        direct_intermediate,
        direct_realization;
        dimension = 25,
    )

    inventory = CRC.unit_pair_inventory((lw_unit, pqs_unit, direct_unit))
    @test CRC.unit_keys(inventory) == (:lw_shell, :pqs_shell, :midpoint)
    @test length(CRC.pair_entries(inventory)) == 6
    @test (:lw_shell, :pqs_shell) in CRC.pair_keys(inventory)

    @test_throws ArgumentError CRC.unit_pair_inventory((filled,))
    @test_throws ArgumentError CRC.pqs_filled_source_lowering(region, facet)
end
