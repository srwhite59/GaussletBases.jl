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

supplement = legacy_bond_aligned_heteronuclear_gaussian_supplement(
    "He",
    "cc-pVTZ",
    "H",
    "cc-pVTZ",
    basis.nuclei;
    lmax = 1,
)

ops = ordinary_cartesian_qiu_white_operators(
    basis,
    supplement;
    nuclear_charges = [2.0, 1.0],
    interaction_treatment = :ggt_nearest,
)

payload = bond_aligned_diatomic_geometry_payload(ops)

projection_path = joinpath(
    @__DIR__,
    "..",
    "plots",
    "hehp_bond_aligned_R1p45_ordinary_hybrid_xz_tol1e-5.dat",
)
centers_path = joinpath(
    @__DIR__,
    "..",
    "plots",
    "hehp_bond_aligned_R1p45_ordinary_hybrid_centers3d.dat",
)

slice = write_bond_aligned_diatomic_plane_projection(
    projection_path,
    payload;
    plane_axis = :y,
    plane_value = 0.0,
    plane_tol = DEBUG_TOL,
)

write_bond_aligned_diatomic_points3d(centers_path, payload)

check = GaussletBases.ordinary_cartesian_1s2_check(ops)
println("projection_path=", projection_path)
println("centers_path=", centers_path)
println("selected_count=", slice.selected_count, "/", slice.total_count)
println("plane_axis=:y plane_value=0.0 plane_tol=", DEBUG_TOL)
println("parallel_spacing_left=", 1.0 / dudx(mapping(basis.basis_z), -0.5 * BOND_LENGTH))
println("parallel_spacing_right=", 1.0 / dudx(mapping(basis.basis_z), 0.5 * BOND_LENGTH))
println("transverse_spacing=", 1.0 / dudx(mapping(basis.basis_x), 0.0))
println("residual_count=", ops.residual_count)
println("overlap_error=", check.overlap_error)
println("orbital_energy=", check.orbital_energy)
println("vee=", check.vee_expectation)
