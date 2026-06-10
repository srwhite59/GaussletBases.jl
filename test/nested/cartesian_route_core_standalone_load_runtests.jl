using Test

const MODULES_DIR = normpath(joinpath(@__DIR__, "..", "..", "modules"))
if !(MODULES_DIR in LOAD_PATH)
    pushfirst!(LOAD_PATH, MODULES_DIR)
end

using CartesianCPB
using CartesianRouteCore

@testset "CartesianRouteCore standalone load-path module" begin
    @test Base.find_package("CartesianCPB") == joinpath(MODULES_DIR, "CartesianCPB.jl")
    @test Base.find_package("CartesianRouteCore") == joinpath(MODULES_DIR, "CartesianRouteCore.jl")

    filled = CartesianRouteCore.filled_cpb(1:3, 1:3, 1:3; role = :standalone_filled)
    facet = CartesianRouteCore.cpb(1:1, 1:3, 1:3; role = :standalone_facet)

    @test CartesianRouteCore.shape(filled) == (3, 3, 3)
    @test CartesianRouteCore.support_count(filled) == 27
    @test CartesianRouteCore.codimension(facet) == 1
    @test CartesianRouteCore.role(facet) == :standalone_facet

    outer = CartesianRouteCore.filled_cpb(1:5, 1:5, 1:5; role = :outer)
    inner = CartesianRouteCore.filled_cpb(2:4, 2:4, 2:4; role = :inner)
    support = CartesianRouteCore.complete_shell_support(outer, inner)
    region = CartesianRouteCore.shellification_region(:standalone_shell, support)
    strata = CartesianRouteCore.complete_shell_boundary_strata(outer, inner)
    lowering = CartesianRouteCore.white_lindsey_boundary_strata_lowering(
        region,
        strata.all_strata,
    )

    @test CartesianRouteCore.support_count(support) == 98
    @test CartesianRouteCore.role(region) == :standalone_shell
    @test CartesianRouteCore.lowering_recipe(lowering) == :white_lindsey_boundary_strata
    @test length(CartesianRouteCore.source_cpbs(lowering)) == 26

    direct_support = CartesianRouteCore.owned_cpb(facet)
    direct_region = CartesianRouteCore.shellification_region(:direct_facet, direct_support)
    direct_lowering = CartesianRouteCore.lowering_source(
        :direct_identity_cpb,
        direct_region,
        direct_support.cpbs,
    )
    retained = CartesianRouteCore.intermediate_retained_space(
        direct_lowering;
        retained_rule = :identity_source_modes,
        source_mode_dims = (1, 3, 3),
    )
    realization = CartesianRouteCore.trivial_shell_realization(retained, direct_region)
    unit = CartesianRouteCore.final_retained_unit(
        :direct_facet_unit,
        :direct_facet,
        direct_lowering,
        retained,
        realization;
        dimension = 9,
    )

    @test retained.dimension == 9
    @test unit.dimension == 9
    @test CartesianRouteCore.role(unit) == :direct_facet
    @test CartesianRouteCore.unit_keys(CartesianRouteCore.unit_pair_inventory((unit,))) ==
          (:direct_facet_unit,)
end
