include("h2_pqs_q5_independent_source_box_r4_supplement_preflight.jl")

materialize_route = true
save_basis_artifact = true
save_ham_artifact = true
basisfile = joinpath(
    tempdir(),
    "h2_pqs_q5_independent_source_box_r4_gto_basis.jld2",
)
hamfile = joinpath(
    tempdir(),
    "h2_pqs_q5_independent_source_box_r4_gto_ham.jld2",
)
