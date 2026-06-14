route_family = :pqs_source_box
route_kind = :bond_aligned_diatomic_fixed_q_complete_core_shell

atom_symbols = ("H", "H")
nuclear_charges = (1, 1)
atom_locations = ((0.0, 0.0, -2.0), (0.0, 0.0, 2.0))
bond_axis = :z
bond_length = 4.0
radius = 4.0

parent_axis_family = :G10
core_spacing = 0.5
xmax_parallel = 6.0
xmax_transverse = 4.0
parent_axis_counts = nothing
probe_parent_axis_construction = true
probe_raw_product_box_plans = false

q = 5
n_s = 5
fixed_source_mode_shape = true

supplement_policy = :none
comparison_ready = false
comparison_blocker = :supplemented_reference_not_comparable_to_gausslet_only
comparison_reference_label = "WL/QW H2 R=4 supplemented reference not used"

run_final_basis = true
run_h1 = true
run_h1_j = false
run_private_rhf = false
save_artifact = true
save_tsv = true
