# Legacy/internal experimental example.
# The old 1D COMX-cleaned hybrid route is not part of the supported public
# GaussletBases workflow. This file is kept only for surrogate regression and
# historical comparison.

using GaussletBases

Z = 2.0
source_basis = build_basis(MappedUniformBasisSpec(:G10;
    count = 7,
    mapping = fit_asinh_mapping_for_strength(s = 0.2, npoints = 7, xmax = 6.0),
))
core_gaussians = [
    Gaussian(center = 0.0, width = 0.2),
    Gaussian(center = 0.0, width = 0.6),
]

pure_ops = ordinary_cartesian_ida_operators(
    source_basis;
    expansion = coulomb_gaussian_expansion(doacc = false),
    Z = Z,
    backend = :pgdg_localized_experimental,
)
hybrid_basis = GaussletBases.hybrid_mapped_ordinary_basis(
    source_basis;
    core_gaussians = core_gaussians,
    backend = :pgdg_localized_experimental,
)
hybrid_ops = ordinary_cartesian_ida_operators(
    hybrid_basis;
    expansion = coulomb_gaussian_expansion(doacc = false),
    Z = Z,
)

pure_check = GaussletBases.ordinary_cartesian_1s2_check(pure_ops)
hybrid_check = GaussletBases.ordinary_cartesian_1s2_check(hybrid_ops)
reference_value = (5.0 / 8.0) * Z

println("Hybrid ordinary Cartesian IDA 1s^2 Vee check")
println("  source basis: ", source_basis)
println("  hybrid basis: ", hybrid_basis)
println("  core Gaussian widths: ", [gaussian.width for gaussian in core_gaussians])
println("  pure mapped dimension: ", length(orbitals(pure_ops)))
println("  hybrid dimension: ", length(orbitals(hybrid_ops)))
println("  pure mapped overlap error: ", pure_check.overlap_error)
println("  hybrid overlap error: ", hybrid_check.overlap_error)
println("  pure mapped lowest one-body orbital energy: ", pure_check.orbital_energy)
println("  hybrid lowest one-body orbital energy: ", hybrid_check.orbital_energy)
println("  pure mapped 1s^2 IDA Vee expectation: ", pure_check.vee_expectation)
println("  hybrid 1s^2 IDA Vee expectation: ", hybrid_check.vee_expectation)
println("  hydrogenic reference (5/8)Z: ", reference_value)
println("  pure mapped difference: ", pure_check.vee_expectation - reference_value)
println("  hybrid difference: ", hybrid_check.vee_expectation - reference_value)
