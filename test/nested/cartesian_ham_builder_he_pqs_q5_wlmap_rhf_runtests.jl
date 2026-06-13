using Test
using JLD2

const _HE_PQS_RHF_DRIVER =
    normpath(joinpath(@__DIR__, "..", "..", "bin", "cartesian_ham_builder.jl"))
const _HE_PQS_RHF_INPUT =
    normpath(joinpath(@__DIR__, "..", "driver_inputs", "he_pqs_q5_wlmap.jl"))

@testset "cartesian_ham_builder_he_pqs_q5_wlmap_rhf artifact" begin
    mktempdir() do dir
        outfile = joinpath(dir, "he_pqs_q5_wlmap_rhf.jld2")
        tsvfile = joinpath(dir, "he_pqs_q5_wlmap_rhf.tsv")
        saved_args = copy(ARGS)
        empty!(ARGS)
        append!(
            ARGS,
            [
                _HE_PQS_RHF_INPUT,
                "run_private_rhf=true",
                "wl_rhf_total=-2.85080350301779",
                "outfile=$(repr(outfile))",
                "tsvfile=$(repr(tsvfile))",
            ],
        )
        try
            include(_HE_PQS_RHF_DRIVER)
        finally
            empty!(ARGS)
            append!(ARGS, saved_args)
        end

        @test isfile(outfile)
        jldopen(outfile, "r") do file
            @test file["private_rhf/converged"]
            @test file["private_rhf/total_energy"] ≈ -2.850817886618113 atol = 1.0e-10 rtol = 0.0
            @test file["comparison/wl_rhf_total"] ≈ -2.85080350301779 atol = 1.0e-12 rtol = 0.0
            @test file["comparison/delta_rhf"] ≈ -1.4383600322798173e-5 atol = 1.0e-10 rtol = 0.0
        end
    end
end
