using Test

const MODULES_DIR = normpath(joinpath(@__DIR__, "..", "..", "modules"))
if !(MODULES_DIR in LOAD_PATH)
    pushfirst!(LOAD_PATH, MODULES_DIR)
end

using CartesianCPB

@testset "CartesianCPB standalone load-path module" begin
    @test Base.find_package("CartesianCPB") == joinpath(MODULES_DIR, "CartesianCPB.jl")

    box = filled_cpb(1:2, 3:4, 5:6; role = :standalone_probe)
    @test intervals(box) == (1:2, 3:4, 5:6)
    @test shape(box) == (2, 2, 2)
    @test codimension(box) == 0
    @test role(box) == :standalone_probe
    @test support_count(box) == 8
end
