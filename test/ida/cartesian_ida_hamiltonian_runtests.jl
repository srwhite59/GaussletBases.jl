using JLD2
using LinearAlgebra
using Test

using GaussletBases

@testset "Matrix EGOI numerical product rank" begin
    Q = Matrix(qr([1. 3.; 2. -2.; 3. 1.; 4. -4.; 5. .5; 6. 2.; 7. -1.; 8. 6.]).Q)
    P = egoi_target_product_matrix(Q)
    A = [1 / (1 + abs(i - j)) for i in 1:8, j in 1:8]
    T = P' * A * P
    retained = Matrix(qr(P[:, [1, 2, 4]]).Q)
    R = retained * retained'
    zero_reg = egoi_density_density_correction(zeros(8, 8), Q, T; regularization = 0.0)
    D = zero_reg.interaction_delta
    @test zero_reg.diagnostics.product_rank == 3
    @test opnorm(D - R * A * R, Inf) <= 1e-10
    @test opnorm(P' * zero_reg.interaction_matrix * P - T, Inf) <= 1e-10
    @test opnorm((I - R) * D, Inf) <= 1e-10
    @test opnorm(D * (I - R), Inf) <= 1e-10
    F = svd(P)
    products = F.S * F.S'
    smooth = F.U * ((F.V' * T * F.V) .* products ./ (1e-18 .+ products.^2)) * F.U'
    @test opnorm(egoi_density_density_correction(zeros(8, 8), Q, T).interaction_delta - smooth, Inf) <= 1e-10
end

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
    charges[1] = 99.0
    positions[2, 3] = 99.0
    @test ham.nuclear_charges ≈ [2.0, 1.0]
    @test ham.nuclear_positions ≈ [0.0 0.0 0.0; 0.0 0.0 2.0]
    @test ham.nuclear_repulsion ≈ 1.0 atol = 0.0 rtol = 0.0
    @test nuclear_repulsion(ham) ≈ ham.nuclear_repulsion atol = 0.0 rtol = 0.0

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
