using Test
using GaussletBases

const CPB = GaussletBases.CartesianCPB

@testset "CartesianCPB contract" begin
    filled = CPB.filled_cpb(1:5, 1:5, 1:5; role = :filled_box)
    slab = CPB.slab_cpb(1:5, 1:5, 3:3; role = :midpoint_slab)
    corner = CPB.cpb(1:1, 1:1, 1:1; role = :corner)

    @test CPB.codimension(filled) == 0
    @test CPB.codimension(slab) == 1
    @test CPB.codimension(corner) == 3

    outer = CPB.filled_cpb(1:5, 1:5, 1:5; role = :outer_box)
    inner = CPB.filled_cpb(2:4, 2:4, 2:4; role = :inner_box)
    strata = CPB.complete_shell_boundary_strata(outer, inner)

    @test length(strata.facets) == 6
    @test length(strata.edges) == 12
    @test length(strata.corners) == 8
    @test length(strata.all_strata) == 26
    @test sum(CPB.support_count, strata.all_strata; init = 0) ==
          CPB.support_count(outer) - CPB.support_count(inner)
    @test strata.stratum_support_count == strata.shell_support_count
end
