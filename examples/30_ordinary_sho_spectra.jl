# Legacy/internal experimental example.
# The old 1D COMX-cleaned hybrid route is not part of the supported public
# GaussletBases workflow. This file is kept only for surrogate regression and
# historical comparison.

using GaussletBases

function print_case(label, reference, analytic)
    diff = maximum(abs.(reference.eigenvalues .- analytic.eigenvalues))
    ref_exact = maximum(abs.(reference.eigenvalues .- reference.exact))
    ana_exact = maximum(abs.(analytic.eigenvalues .- analytic.exact))

    println(label)
    println("  exact:      ", reference.exact)
    println("  numerical:  ", reference.eigenvalues)
    println("  analytic:   ", analytic.eigenvalues)
    println("  max |analytic - numerical|: ", diff)
    println("  max |numerical - exact|: ", ref_exact)
    println("  max |analytic - exact|: ", ana_exact)
    println("  <T>_0 numerical / analytic: ", reference.kinetic_expectation, " / ", analytic.kinetic_expectation)
    println(
        "  <(x-x0)^2>_0 numerical / analytic: ",
        reference.displacement2_expectation,
        " / ",
        analytic.displacement2_expectation,
    )
end

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

centered_reference = ordinary_sho_spectrum(hybrid_reference; omega = 0.25, center = 0.0, nev = 3)
centered_analytic = ordinary_sho_spectrum(hybrid_analytic; omega = 0.25, center = 0.0, nev = 3)

shifted_reference = ordinary_sho_spectrum(hybrid_reference; omega = 0.25, center = 1.0, nev = 3)
shifted_analytic = ordinary_sho_spectrum(hybrid_analytic; omega = 0.25, center = 1.0, nev = 3)

stiff_reference = ordinary_sho_spectrum(hybrid_reference; omega = 0.75, center = 0.0, nev = 3)
stiff_analytic = ordinary_sho_spectrum(hybrid_analytic; omega = 0.75, center = 0.0, nev = 3)

stress_basis = build_basis(MappedUniformBasisSpec(:G10;
    count = 5,
    mapping = fit_asinh_mapping_for_strength(s = 2.0, npoints = 5, xmax = 6.0),
))
stress_reference = ordinary_sho_spectrum(
    stress_basis;
    omega = 0.5,
    center = 0.0,
    nev = 3,
    backend = :numerical_reference,
)
stress_analytic = ordinary_sho_spectrum(
    stress_basis;
    omega = 0.5,
    center = 0.0,
    nev = 3,
    backend = :pgdg_localized_experimental,
)

println("Ordinary mapped SHO spectral comparison")
println("  source basis: ", basis)
println("  hybrid reference: ", hybrid_reference)
println("  hybrid analytic: ", hybrid_analytic)
println("  core Gaussian widths: ", [gaussian.width for gaussian in core_gaussians])
print_case("  centered mild hybrid (omega = 0.25, center = 0.0)", centered_reference, centered_analytic)
print_case("  shifted mild hybrid (omega = 0.25, center = 1.0)", shifted_reference, shifted_analytic)
print_case("  stiffer mild hybrid (omega = 0.75, center = 0.0)", stiff_reference, stiff_analytic)
print_case("  pure mapped stress test (omega = 0.5, center = 0.0)", stress_reference, stress_analytic)
