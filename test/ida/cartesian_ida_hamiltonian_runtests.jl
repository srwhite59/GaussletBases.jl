using JLD2
using Test

using GaussletBases

@testset "Public Cartesian IDA Hamiltonian contract" begin
    kinetic = [1.0 0.1; 0.1 1.4]
    unit_nuclear = Matrix{Float64}[
        [-0.8 -0.05; -0.05 -0.3],
        [-0.2 -0.02; -0.02 -0.6],
    ]
    vee = [0.7 0.2; 0.2 0.9]
    charges = [2.0, 1.0]
    positions = [0.0 0.0 0.0; 0.0 0.0 2.0]

    ham = CartesianIDAHamiltonian(
        kinetic,
        unit_nuclear,
        vee,
        1,
        1;
        nuclear_charges = charges,
        nuclear_positions = positions,
    )

    @test ham.nuclear_repulsion ≈ 1.0 atol = 0.0 rtol = 0.0
    @test one_body_hamiltonian(ham) ≈
          kinetic + charges[1] .* unit_nuclear[1] + charges[2] .* unit_nuclear[2]
    @test one_body_hamiltonian(ham; center_weights = [1.0, 0.0]) ≈
          kinetic + charges[1] .* unit_nuclear[1]
    @test nuclear_repulsion(ham; center_weights = [1.0, 0.0]) == 0.0

    mktempdir() do dir
        path = joinpath(dir, "cartesian_ida_hamiltonian.jld2")
        @test write_cartesian_ida_hamiltonian(path, ham) == path
        jldopen(path, "r") do file
            @test file["artifact_kind"] === :cartesian_ida_hamiltonian
            @test Int(file["format_version"]) == 1
            @test size(file["nuclear_attraction_unit_by_center"]) == (2, 2, 2)
            @test !haskey(file, "nuclear_repulsion")
        end
        reloaded = read_cartesian_ida_hamiltonian(path)
        @test reloaded.kinetic ≈ ham.kinetic
        @test length(reloaded.nuclear_attraction_unit_by_center) ==
              length(ham.nuclear_attraction_unit_by_center)
        for center in eachindex(ham.nuclear_attraction_unit_by_center)
            @test reloaded.nuclear_attraction_unit_by_center[center] ≈
                  ham.nuclear_attraction_unit_by_center[center]
        end
        @test reloaded.electron_electron_ida ≈ ham.electron_electron_ida
        @test reloaded.nuclear_charges ≈ ham.nuclear_charges
        @test reloaded.nuclear_positions ≈ ham.nuclear_positions
        @test reloaded.nup == ham.nup
        @test reloaded.ndn == ham.ndn
        @test reloaded.nuclear_repulsion ≈ ham.nuclear_repulsion
    end
end
