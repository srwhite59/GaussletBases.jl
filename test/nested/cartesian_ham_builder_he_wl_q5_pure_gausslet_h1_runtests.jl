using Test
using JLD2

repo_root = abspath(joinpath(@__DIR__, "..", ".."))
hamfile = tempname() * ".jld2"

redirect_stdout(devnull) do
    empty!(ARGS)
    append!(
        ARGS,
        [
            joinpath("test", "driver_inputs", "he_wl_q5_pure_gausslet_h1.jl"),
            "save_ham_artifact=true",
            "hamfile=$(repr(hamfile))",
        ],
    )
    include(joinpath(repo_root, "bin", "cartesian_ham_builder.jl"))
end

@test materialization.status == :materialized_white_lindsey_atomic_pure_gausslet
@test materialization.retained_dimension == 419
@test materialization.h1_finite
@test materialization.h1_symmetry_error <= 1.0e-12
@test materialization.overlap_identity_error <= 1.0e-12
@test materialization.h1_lowest ≈ -1.991344469963435 atol = 1.0e-12
@test materialization.jld2_ham_artifact_status == :written
@test isfile(hamfile)

jldopen(hamfile, "r") do file
    @test file["artifact_kind"] == :white_lindsey_atomic_pure_gausslet_hamiltonian
    h1 = file["one_body_hamiltonian"]
    @test size(h1) == (419, 419)
    @test all(isfinite, h1)
    @test maximum(abs.(h1 .- transpose(h1))) <= 1.0e-12
end
