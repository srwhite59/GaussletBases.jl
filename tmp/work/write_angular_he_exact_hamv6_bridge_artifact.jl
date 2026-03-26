using GaussletBases

const OUTPUT_PATH = normpath(joinpath(
    @__DIR__,
    "..",
    "work",
    "angular_bridge",
    "he_angular_exact_common_hamv6.jld2",
))

function build_he_benchmark()
    Z = 2.0
    s = 0.15
    basis = build_basis(
        RadialBasisSpec(
            :G10;
            count = 6,
            mapping = AsinhMapping(c = s, s = s),
            reference_spacing = 1.0,
            tails = 3,
            odd_even_kmax = 2,
            xgaussians = [XGaussian(alpha = 0.2)],
        ),
    )
    grid = radial_quadrature(basis; accuracy = :medium, quadrature_rmax = 12.0)
    radial_ops = atomic_operators(basis, grid; Z = Z, lmax = 2)
    return build_atomic_injected_angular_small_ed_benchmark(radial_ops)
end

function main()
    benchmark = build_he_benchmark()
    mkpath(dirname(OUTPUT_PATH))
    write_angular_benchmark_exact_hamv6_jld2(
        OUTPUT_PATH,
        benchmark;
        nelec = 2,
        meta = (
            example = "angular_bridge_he",
            benchmark_line = "atomic_injected_angular_small_ed",
        ),
    )

    diagnostics = atomic_injected_angular_small_ed_diagnostics(benchmark)
    println("artifact_path=$(OUTPUT_PATH)")
    println("full_energy=$(diagnostics.full_energy)")
    println("exact_reference_energy=$(diagnostics.exact_reference_energy)")
    println("exact_common_lmax=$(diagnostics.exact_common_lmax)")
    println("shell_orders=$(join(diagnostics.shell_orders, ","))")
end

main()
