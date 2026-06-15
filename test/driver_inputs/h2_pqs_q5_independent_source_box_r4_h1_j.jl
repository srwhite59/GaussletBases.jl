include("h2_pqs_q5_independent_source_box_r4.jl")

artifact_role = :independent_h2_pqs_h1_j_density_diagnostic
physics_endpoint_blocker = :missing_physical_gausslet_rhf_or_solver_contract
run_final_basis = true
run_h1 = true
run_h1_j = true
run_private_rhf = false
