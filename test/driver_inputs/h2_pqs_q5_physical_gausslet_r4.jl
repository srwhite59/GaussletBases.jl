route_family = :pqs_source_box
route_kind = :bond_aligned_diatomic_physical_gausslet_core_shell_pqs
route_shape = (:atom_contact_core, :shared_shell_1, :shared_shell_2)

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
comparison_blocker = :gausslet_only_reference_not_selected
comparison_reference_label = "WL/QW H2 R=4 gausslet-only target inventory only"
artifact_role = :physical_gausslet_endpoint_target
physics_endpoint_ready = false
physics_endpoint_blocker = :missing_atom_contact_core_support_rows
retained_atom_core_interiors = true
source_plan_role = :atom_contact_core_plus_pqs_shared_shells

run_final_basis = false
run_h1 = false
run_h1_j = false
run_private_rhf = false
save_artifact = true
save_tsv = true
