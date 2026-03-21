using GaussletBases

basis_name = "cc-pVTZ"
count = 9
s = 0.8
xmax = 6.0
z_value = 2.0
reference_value = 1.25

source_basis = build_basis(MappedUniformBasisSpec(:G10;
    count = count,
    mapping = fit_asinh_mapping_for_strength(s = s, npoints = count, xmax = xmax),
    reference_spacing = 1.0,
))
legacy = legacy_s_gaussian_data("He", basis_name)
expansion = coulomb_gaussian_expansion(doacc = false)

hybrid_basis = hybrid_mapped_ordinary_basis(
    source_basis;
    core_gaussians = legacy,
    backend = :pgdg_localized_experimental,
)
surrogate_mwg = ordinary_cartesian_ida_operators(
    hybrid_basis;
    expansion = expansion,
    Z = z_value,
    interaction_treatment = :residual_gaussian_mwg,
)
qiu_white_nearest = ordinary_cartesian_qiu_white_operators(
    source_basis,
    legacy;
    expansion = expansion,
    Z = z_value,
    interaction_treatment = :ggt_nearest,
)
qiu_white_mwg = ordinary_cartesian_qiu_white_operators(
    source_basis,
    legacy;
    expansion = expansion,
    Z = z_value,
    interaction_treatment = :mwg,
)

surrogate_check = GaussletBases.ordinary_cartesian_1s2_check(surrogate_mwg)
nearest_check = GaussletBases.ordinary_cartesian_1s2_check(qiu_white_nearest)
mwg_check = GaussletBases.ordinary_cartesian_1s2_check(qiu_white_mwg)

sample_count = min(3, qiu_white_mwg.residual_count)

println("Qiu-White residual-Gaussian reference check")
println("  basis: ", source_basis)
println("  He basis supplement: ", legacy.basis_name)
println("  primitive exponents: ", legacy.primitive_exponents)
println("  primitive widths: ", legacy.primitive_widths)
println("  contraction matrix size: ", size(legacy.contraction_matrix))
println("  contracted widths: ", legacy.widths)
println("  count = ", count, ", s = ", s, ", xmax = ", xmax)
println("  note: this is a slow full-expansion light reference run, not a quick smoke check.")
println("  hydrogenic 1s^2 target: ", reference_value)
println()

println("Current surrogate MWG path")
println("  basis: ", hybrid_basis)
println("  E1: ", surrogate_check.orbital_energy)
println("  <Vee>: ", surrogate_check.vee_expectation)
println("  difference from 1.25: ", surrogate_check.vee_expectation - reference_value)
println()

println("Paper-faithful Qiu-White nearest / GGT path")
println("  E1: ", nearest_check.orbital_energy)
println("  <Vee>: ", nearest_check.vee_expectation)
println("  difference from 1.25: ", nearest_check.vee_expectation - reference_value)
println()

println("Paper-faithful Qiu-White MWG path")
println("  E1: ", mwg_check.orbital_energy)
println("  <Vee>: ", mwg_check.vee_expectation)
println("  difference from 1.25: ", mwg_check.vee_expectation - reference_value)
println("  representative residual centers:")
println(qiu_white_mwg.residual_centers[1:sample_count, :])
println("  representative residual widths:")
println(qiu_white_mwg.residual_widths[1:sample_count, :])
