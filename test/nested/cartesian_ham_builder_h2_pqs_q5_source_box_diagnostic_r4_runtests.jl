using Test, JLD2

const _H2_PQS_DRIVER = normpath(joinpath(@__DIR__, "..", "..", "bin", "cartesian_ham_builder.jl"))
const _H2_PQS_INPUT = normpath(joinpath(@__DIR__, "..", "driver_inputs", "h2_pqs_q5_source_box_diagnostic_r4.jl"))

@testset "cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4 artifact" begin
    mktempdir() do dir
        outfile = joinpath(dir, "h2_pqs_q5_source_box_diagnostic_r4.jld2")
        tsvfile = joinpath(dir, "h2_pqs_q5_source_box_diagnostic_r4.tsv")
        saved_args = copy(ARGS)
        empty!(ARGS)
        append!(ARGS, [_H2_PQS_INPUT, "outfile=$(repr(outfile))", "tsvfile=$(repr(tsvfile))"])
        try include(_H2_PQS_DRIVER) finally empty!(ARGS); append!(ARGS, saved_args) end

        @test isfile(outfile)
        jldopen(outfile, "r") do file
            @test file["route/artifact_role"] === :source_box_diagnostic
            @test file["physics/endpoint_ready"] == false
            @test file["physics/endpoint_blocker"] === :retained_atom_core_interiors_missing
            @test file["comparison/ready"] == false
            @test file["basis/final_dimension"] == 221
            @test all(!haskey(file, key) for key in ("comparison/wl_h1_lowest", "comparison/delta_h1", "comparison/wl_rhf_total", "comparison/delta_rhf"))
            @test file["private_rhf/requested"] == false
            @test file["route/h1_materialized"] == true
            @test isfinite(file["physics/h1_lowest"])
            @test file["physics/h1_hamiltonian_matrix_finite"] == true
            @test file["physics/h1_hamiltonian_symmetry_error"] < 1e-10
        end
    end
end
