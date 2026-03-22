using GaussletBases

function run_case(label::String, core_gaussians)
    source_basis = build_basis(MappedUniformBasisSpec(:G10;
        count = 11,
        mapping = fit_asinh_mapping_for_strength(s = 0.6, npoints = 11, xmax = 6.0),
    ))
    operators = if label == "pure ordinary"
        ordinary_cartesian_ida_operators(
            source_basis;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
            backend = :pgdg_localized_experimental,
        )
    else
        hybrid_basis = hybrid_mapped_ordinary_basis(
            source_basis;
            core_gaussians = core_gaussians,
            backend = :pgdg_localized_experimental,
        )
        ordinary_cartesian_ida_operators(
            hybrid_basis;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
        )
    end
    check = GaussletBases.ordinary_cartesian_1s2_check(operators)
    return (
        label = label,
        E1 = check.orbital_energy,
        Vee = check.vee_expectation,
        delta = check.vee_expectation - 1.25,
    )
end

function print_legacy(label::String, data::LegacySGaussianData)
    println(label)
    println("  primitive exponents: ", data.primitive_exponents)
    println("  primitive widths: ", data.primitive_widths)
    println("  contraction matrix size: ", size(data.contraction_matrix))
    println("  contracted widths: ", data.widths)
end

toy_gaussians = [
    Gaussian(center = 0.0, width = 0.2),
    Gaussian(center = 0.0, width = 0.6),
]
vtz = legacy_s_gaussian_data("He", "cc-pVTZ")
vqz = legacy_s_gaussian_data("He", "cc-pVQZ")

println("Hybrid ordinary Cartesian He s-basis supplement check")
println("  count = 11")
println("  s = 0.6")
println("  xmax = 6.0")
println("  backend = :pgdg_localized_experimental")
println("  using contracted He s shells, with analytic primitive integrals underneath")
println()

print_legacy("He cc-pVTZ s primitives", vtz)
println()
print_legacy("He cc-pVQZ s primitives", vqz)
println()

for result in (
    run_case("pure ordinary", Gaussian[]),
    run_case("toy widths [0.2, 0.6]", toy_gaussians),
    run_case("He cc-pVTZ contracted s supplement", vtz),
    run_case("He cc-pVQZ contracted s supplement", vqz),
)
    println(result.label)
    println("  E1: ", result.E1)
    println("  <Vee>: ", result.Vee)
    println("  difference from 1.25: ", result.delta)
    println()
end
