using GaussletBases

function build_ordinary_h2_payload(; bond_length::Float64 = 1.4)
    basis = bond_aligned_homonuclear_qw_basis(
        bond_length = bond_length,
        core_spacing = 0.5,
        xmax_parallel = 8,
        xmax_transverse = 5,
        bond_axis = :z,
    )
    supplement = legacy_bond_aligned_diatomic_gaussian_supplement(
        "H",
        "cc-pVTZ",
        basis.nuclei;
        lmax = 1,
    )
    ops = ordinary_cartesian_qiu_white_operators(
        basis,
        supplement;
        nuclear_charges = [1.0, 1.0],
        interaction_treatment = :ggt_nearest,
    )
    return bond_aligned_diatomic_geometry_payload(ops)
end

function build_nested_h2_payloads(; bond_length::Float64 = 1.4)
    basis = bond_aligned_homonuclear_qw_basis(
        bond_length = bond_length,
        core_spacing = 0.5,
        xmax_parallel = 8,
        xmax_transverse = 5,
        bond_axis = :z,
    )
    source = GaussletBases._bond_aligned_diatomic_nested_fixed_source(basis)
    fixed_block = GaussletBases._nested_fixed_block(source)
    supplement = legacy_bond_aligned_diatomic_gaussian_supplement(
        "H",
        "cc-pVTZ",
        basis.nuclei;
        lmax = 1,
    )
    ops = ordinary_cartesian_qiu_white_operators(
        fixed_block,
        supplement;
        nuclear_charges = [1.0, 1.0],
        interaction_treatment = :ggt_nearest,
    )
    return (
        fixed_payload = bond_aligned_diatomic_geometry_payload(ops, source),
        source_payload = bond_aligned_diatomic_source_geometry_payload(source),
    )
end

const DEBUG_TOL = 1.0e-5

ordinary_payload = build_ordinary_h2_payload()
nested_payloads = build_nested_h2_payloads()

ordinary_projection_path = joinpath(
    @__DIR__,
    "..",
    "plots",
    "h2_bond_aligned_R1p4_ordinary_hybrid_xz_tol1e-5.dat",
)
nested_projection_path = joinpath(
    @__DIR__,
    "..",
    "plots",
    "h2_bond_aligned_R1p4_nested_hybrid_xz_tol1e-5.dat",
)
nested_centers_path = joinpath(
    @__DIR__,
    "..",
    "plots",
    "h2_bond_aligned_R1p4_nested_hybrid_centers3d.dat",
)
nested_source_projection_path = joinpath(
    @__DIR__,
    "..",
    "plots",
    "h2_bond_aligned_R1p4_nested_source_xz_tol1e-5.dat",
)
nested_source_path = joinpath(
    @__DIR__,
    "..",
    "plots",
    "h2_bond_aligned_R1p4_nested_source_regions3d.dat",
)

ordinary_slice = write_bond_aligned_diatomic_plane_projection(
    ordinary_projection_path,
    ordinary_payload;
    plane_axis = :y,
    plane_value = 0.0,
    plane_tol = DEBUG_TOL,
)
nested_slice = write_bond_aligned_diatomic_plane_projection(
    nested_projection_path,
    nested_payloads.fixed_payload;
    plane_axis = :y,
    plane_value = 0.0,
    plane_tol = DEBUG_TOL,
)
nested_source_slice = write_bond_aligned_diatomic_plane_projection(
    nested_source_projection_path,
    nested_payloads.source_payload;
    plane_axis = :y,
    plane_value = 0.0,
    plane_tol = DEBUG_TOL,
    include_box_outlines = true,
)
write_bond_aligned_diatomic_points3d(nested_centers_path, nested_payloads.fixed_payload)
write_bond_aligned_diatomic_points3d(nested_source_path, nested_payloads.source_payload)

strict_nested_slice = bond_aligned_diatomic_plane_slice(
    nested_payloads.fixed_payload;
    plane_axis = :y,
    plane_value = 0.0,
    plane_tol = 1.0e-12,
)
source_debug_slice = bond_aligned_diatomic_plane_slice(
    nested_payloads.source_payload;
    plane_axis = :y,
    plane_value = 0.0,
    plane_tol = DEBUG_TOL,
)

println("ordinary_projection_path=", ordinary_projection_path)
println("ordinary_selected=", ordinary_slice.selected_count, "/", ordinary_slice.total_count)
println("nested_projection_path=", nested_projection_path)
println("nested_selected_debug=", nested_slice.selected_count, "/", nested_slice.total_count)
println("nested_selected_strict=", strict_nested_slice.selected_count, "/", strict_nested_slice.total_count)
println("nested_centers_path=", nested_centers_path)
println("nested_source_projection_path=", nested_source_projection_path)
println("nested_source_projection_selected_debug=", nested_source_slice.selected_count, "/", nested_source_slice.total_count)
println("nested_source_path=", nested_source_path)
println("nested_source_selected_debug=", source_debug_slice.selected_count, "/", source_debug_slice.total_count)
println("plane_axis=:y plane_value=0.0 plane_tol=", DEBUG_TOL)
