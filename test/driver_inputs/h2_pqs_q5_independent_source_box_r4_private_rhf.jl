include("h2_pqs_q5_independent_source_box_r4_h1_j.jl")

physics_endpoint_ready = false
physics_endpoint_blocker = :private_rhf_diagnostic_not_public_solver_contract
run_final_basis = true
run_h1 = true
run_h1_j = true
run_private_rhf = true
private_rhf_electron_count = 2
