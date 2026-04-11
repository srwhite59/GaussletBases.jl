# Legacy/internal experimental example.
# The old 1D COMX-cleaned hybrid route is not part of the supported public
# GaussletBases workflow. This file is kept only for surrogate regression and
# historical comparison.

using LinearAlgebra
using GaussletBases

function truncate_expansion(expansion::CoulombGaussianExpansion, nterms::Int)
    return CoulombGaussianExpansion(
        expansion.coefficients[1:nterms],
        expansion.exponents[1:nterms];
        del = expansion.del,
        s = expansion.s,
        c = expansion.c,
        maxu = expansion.maxu,
    )
end

full_expansion = coulomb_gaussian_expansion(doacc = false)
hydrogen_expansion = truncate_expansion(full_expansion, 5)

basis = build_basis(MappedUniformBasisSpec(:G10;
    count = 5,
    mapping = fit_asinh_mapping_for_strength(s = 0.2, npoints = 5, xmax = 6.0),
))
core_gaussians = [
    Gaussian(center = 0.0, width = 0.2),
    Gaussian(center = 0.0, width = 0.6),
]

hybrid_reference = GaussletBases.hybrid_mapped_ordinary_basis(
    basis;
    core_gaussians = core_gaussians,
    backend = :numerical_reference,
)
hybrid_analytic = GaussletBases.hybrid_mapped_ordinary_basis(
    basis;
    core_gaussians = core_gaussians,
    backend = :pgdg_localized_experimental,
)

reference_1d = mapped_ordinary_one_body_operators(
    hybrid_reference;
    exponents = hydrogen_expansion.exponents,
)
analytic_1d = mapped_ordinary_one_body_operators(
    hybrid_analytic;
    exponents = hydrogen_expansion.exponents,
)

factor_error = maximum(
    norm(analytic_1d.gaussian_factors[index] - reference_1d.gaussian_factors[index], Inf)
    for index in eachindex(analytic_1d.gaussian_factors)
)

energy_reference = mapped_cartesian_hydrogen_energy(
    hybrid_reference;
    expansion = hydrogen_expansion,
    Z = 1.0,
)
energy_analytic = mapped_cartesian_hydrogen_energy(
    hybrid_analytic;
    expansion = hydrogen_expansion,
    Z = 1.0,
)

println("Hybrid mapped Cartesian hydrogen")
println("  source basis: ", basis)
println("  hybrid reference: ", hybrid_reference)
println("  hybrid analytic: ", hybrid_analytic)
println("  core Gaussian widths: ", [gaussian.width for gaussian in core_gaussians])
println("  1D overlap difference: ", norm(analytic_1d.overlap - reference_1d.overlap, Inf))
println("  1D kinetic difference: ", norm(analytic_1d.kinetic - reference_1d.kinetic, Inf))
println("  max Gaussian-factor difference: ", factor_error)
println("  hydrogen reference energy: ", energy_reference)
println("  hydrogen analytic energy: ", energy_analytic)
println("  hydrogen energy difference: ", abs(energy_analytic - energy_reference))
