using Test
using GaussletBases
using GaussletBases.CartesianRouteCore

@testset "CartesianRouteCore documented examples" begin
    outer = filled_cpb(1:5, 1:5, 1:5; role = :outer_box)
    inner = filled_cpb(2:4, 2:4, 2:4; role = :inner_box)
    support = complete_shell_support(outer, inner)
    region = shellification_region(:atom_local_shell, support)
    strata = complete_shell_boundary_strata(outer, inner)
    lowering = white_lindsey_boundary_strata_lowering(region, strata.all_strata)

    @test lowering isa LoweringSource
    @test lowering_recipe(lowering) == :white_lindsey_boundary_strata
    @test support isa OwnedSupport
    @test support_count(support) == 98
    @test length(source_cpbs(lowering)) == 26

    source = filled_cpb(1:5, 1:5, 1:5; role = :pqs_source_box)
    lowering = pqs_filled_source_lowering(region, source)
    space = intermediate_retained_space(
        lowering;
        retained_rule = :pqs_boundary_comx_product_modes,
        source_mode_dims = (5, 5, 5),
    )
    realization = pqs_shell_realization(space, region)

    @test lowering_recipe(lowering) == :pqs_filled_source_cpb
    @test only(source_cpbs(lowering)) === source
    @test space isa IntermediateRetainedSpace
    @test space.dimension == 98
    @test realization isa ShellRealization
    @test realization.realization_kind == :shell_projection_lowdin
end
