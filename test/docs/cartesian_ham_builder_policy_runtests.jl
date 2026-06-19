@testset "Canonical Cartesian Ham Builder Policy" begin
    driver_path = joinpath(_PROJECT_ROOT, "bin", "cartesian_ham_builder.jl")
    driver = read(driver_path, String)

    @test occursin("Canonical human-facing Cartesian producer template", driver)
    @test !occursin("private_rhf", driver)
    @test !occursin("stop_after_stage", driver)
    @test !occursin("GaussletBases._", driver)
    @test !occursin("residual_gto_provider_blocks", driver)
    @test !occursin(":one_body_and_density_provider", driver)
    @test !occursin("Meta.parse", driver)

    stage_calls = [
        "GaussletBases.cartesian_system(",
        "GaussletBases.cartesian_recipe(",
        "GaussletBases.cartesian_parent(",
        "GaussletBases.cartesian_shells(",
        "GaussletBases.cartesian_units(",
        "GaussletBases.cartesian_transforms(",
        "GaussletBases.cartesian_pair_terms(",
        "GaussletBases.cartesian_assembly(",
        "GaussletBases.cartesian_report(",
        "GaussletBases.cartesian_materialization(",
        "GaussletBases.cartesian_print_summary(",
        "GaussletBases.cartesian_save(",
    ]
    positions = Int[]
    for call in stage_calls
        location = findfirst(call, driver)
        @test !isnothing(location)
        isnothing(location) || push!(positions, first(location))
    end
    @test issorted(positions; lt = <)
end
