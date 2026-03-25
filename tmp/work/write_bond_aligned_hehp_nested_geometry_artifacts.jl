using GaussletBases

const DEBUG_TOL = 1.0e-5
const BOND_LENGTH = 1.45

basis = bond_aligned_heteronuclear_qw_basis(
    atoms = ("He", "H"),
    bond_length = BOND_LENGTH,
    core_spacings = (0.25, 0.5),
    nuclear_charges = (2.0, 1.0),
    xmax_parallel = 6.0,
    xmax_transverse = 4.0,
    bond_axis = :z,
)

expansion = coulomb_gaussian_expansion(doacc = false)
source = GaussletBases._bond_aligned_diatomic_nested_fixed_source(
    basis;
    expansion = expansion,
)
fixed_block = GaussletBases._nested_fixed_block(source)

supplement = legacy_bond_aligned_heteronuclear_gaussian_supplement(
    "He",
    "cc-pVTZ",
    "H",
    "cc-pVTZ",
    basis.nuclei;
    lmax = 1,
)

ops = ordinary_cartesian_qiu_white_operators(
    fixed_block,
    supplement;
    nuclear_charges = [2.0, 1.0],
    expansion = expansion,
    interaction_treatment = :ggt_nearest,
)

fixed_payload = bond_aligned_diatomic_geometry_payload(ops, source)
source_payload = bond_aligned_diatomic_source_geometry_payload(source)

fixed_projection_path = joinpath(
    @__DIR__,
    "..",
    "plots",
    "hehp_bond_aligned_R1p45_nested_hybrid_xz_tol1e-5.dat",
)
fixed_centers_path = joinpath(
    @__DIR__,
    "..",
    "plots",
    "hehp_bond_aligned_R1p45_nested_hybrid_centers3d.dat",
)
source_projection_path = joinpath(
    @__DIR__,
    "..",
    "plots",
    "hehp_bond_aligned_R1p45_nested_source_xz_tol1e-5.dat",
)
source_regions_path = joinpath(
    @__DIR__,
    "..",
    "plots",
    "hehp_bond_aligned_R1p45_nested_source_regions3d.dat",
)

fixed_slice = write_bond_aligned_diatomic_plane_projection(
    fixed_projection_path,
    fixed_payload;
    plane_axis = :y,
    plane_value = 0.0,
    plane_tol = DEBUG_TOL,
)
source_slice = write_bond_aligned_diatomic_plane_projection(
    source_projection_path,
    source_payload;
    plane_axis = :y,
    plane_value = 0.0,
    plane_tol = DEBUG_TOL,
    include_box_outlines = true,
)
write_bond_aligned_diatomic_points3d(fixed_centers_path, fixed_payload)
write_bond_aligned_diatomic_points3d(source_regions_path, source_payload)

check = GaussletBases.ordinary_cartesian_1s2_check(ops)
println("fixed_projection_path=", fixed_projection_path)
println("fixed_centers_path=", fixed_centers_path)
println("source_projection_path=", source_projection_path)
println("source_regions_path=", source_regions_path)
println("fixed_selected=", fixed_slice.selected_count, "/", fixed_slice.total_count)
println("source_selected=", source_slice.selected_count, "/", source_slice.total_count)
println("plane_axis=:y plane_value=0.0 plane_tol=", DEBUG_TOL)
println("split_index=", source.geometry.split_index)
println("shared_midpoint_box=", source.geometry.shared_midpoint_box)
println("child_boxes=", source.geometry.child_boxes)
println("child_physical_widths=", source.geometry.child_physical_widths)
println("fixed_count=", size(fixed_block.overlap, 1))
println("residual_count=", ops.residual_count)
println("overlap_error=", check.overlap_error)
println("orbital_energy=", check.orbital_energy)
println("vee=", check.vee_expectation)
