using GaussletBases

Z = 1.0
npoints = 5
rmax = 6.0

mapping = fit_asinh_mapping_for_strength(s = 0.5, npoints = npoints, xmax = rmax)
basis = build_basis(MappedUniformBasisSpec(:G10;
    count = npoints,
    mapping = mapping,
    reference_spacing = 1.0,
))

expansion = coulomb_gaussian_expansion(doacc = false)
analytic = mapped_ordinary_one_body_operators(
    basis;
    exponents = expansion.exponents,
    backend = :pgdg_experimental,
)
reference = mapped_ordinary_one_body_operators(
    basis;
    exponents = expansion.exponents,
    backend = :numerical_reference,
)

energy_analytic = mapped_cartesian_hydrogen_energy(analytic, expansion; Z = Z)
energy_reference = mapped_cartesian_hydrogen_energy(reference, expansion; Z = Z)

println("Mapped Cartesian hydrogen with backend comparison")
println("  basis: ", basis)
println("  mapping: ", mapping)
println("  analytic backend: ", analytic)
println("  reference backend: ", reference)
println("  analytic energy: ", energy_analytic)
println("  numerical reference energy: ", energy_reference)
println("  difference: ", energy_analytic - energy_reference)
