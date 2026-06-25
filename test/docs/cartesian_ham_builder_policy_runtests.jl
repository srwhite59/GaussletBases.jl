@testset "Canonical Cartesian Ham Builder Policy" begin
    driver_path = joinpath(_PROJECT_ROOT, "bin", "cartesian_ham_builder.jl")
    driver = read(driver_path, String)

    @test occursin("Canonical human-facing Cartesian producer", driver)
    @test !occursin("private_rhf", driver)
    @test !occursin("stop_after_stage", driver)
    @test !occursin("GaussletBases._", driver)
    @test !occursin("residual_gto_provider_blocks", driver)
    @test !occursin(":one_body_and_density_provider", driver)

    retired_calls = [
        join(("GaussletBases.", "cartesian_", "materialization", "(")),
        join(("GaussletBases.", "cartesian_", "print_summary", "(")),
        join(("GaussletBases.", "cartesian_", "print_details", "(")),
        join(("GaussletBases.", "cartesian_", "save", "(")),
    ]
    for call in retired_calls
        @test !occursin(call, driver)
    end
end
