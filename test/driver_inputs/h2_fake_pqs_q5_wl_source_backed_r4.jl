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
comparison_ready = true
comparison_blocker = nothing
comparison_reference_label = "WL/QW H2 R=4 gausslet-only 463 reproduction"
artifact_role = :fake_pqs_source_backed_wl_reproduction
physics_endpoint_ready = false
physics_endpoint_blocker = :fake_pqs_source_backed_wl_reproduction_not_independent_pqs
retained_atom_core_interiors = true
source_plan_role = :fake_pqs_source_backed_fixed_source_adapter

wl_h1_lowest = -0.7946609179724673
wl_h1_self_coulomb = 0.45696639804337047
wl_rhf_one_electron_energy = -1.5611571934181985
wl_rhf_electron_electron_energy = 0.40220533775308426
wl_rhf_electronic_energy = -1.1589518556651142
wl_rhf_nuclear_repulsion = 0.25
wl_rhf_total_with_nuclear_repulsion = -0.9089518556651142

run_final_basis = true
run_h1 = true
run_h1_j = true
run_private_rhf = true
private_rhf_electron_count = 2
save_artifact = true
save_tsv = true
