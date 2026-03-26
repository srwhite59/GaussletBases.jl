using GaussletBases

const HFDMRG_PATH = "/Users/srw/Dropbox/codexhome/work/hfdmrg"
if !(HFDMRG_PATH in LOAD_PATH)
    push!(LOAD_PATH, HFDMRG_PATH)
end
using HFDMRG

function build_open_shell_benchmark()
    Z = 3.0
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
    return build_atomic_injected_angular_hf_style_benchmark(radial_ops)
end

function main()
    benchmark = build_open_shell_benchmark()
    default_seeds = build_atomic_injected_angular_hfdmrg_hf_seeds(benchmark; nup = 2, ndn = 1)
    explicit_psiup0 = default_seeds.psiup0[:, [2, 1]]
    explicit_psidn0 = default_seeds.psidn0

    default_adapter = build_atomic_injected_angular_hfdmrg_hf_adapter(
        benchmark;
        nup = 2,
        ndn = 1,
    )
    explicit_adapter = build_atomic_injected_angular_hfdmrg_hf_adapter(
        benchmark;
        nup = 2,
        ndn = 1,
        psiup0 = explicit_psiup0,
        psidn0 = explicit_psidn0,
    )

    default_result = run_atomic_injected_angular_hfdmrg_hf(
        default_adapter;
        hfmod = HFDMRG,
        maxiter = 40,
        blocksize = 16,
        cutoff = 1.0e-10,
        scf_cutoff = 1.0e-11,
        verbose = false,
    )
    explicit_result = run_atomic_injected_angular_hfdmrg_hf(
        explicit_adapter;
        hfmod = HFDMRG,
        maxiter = 40,
        blocksize = 16,
        cutoff = 1.0e-10,
        scf_cutoff = 1.0e-11,
        verbose = false,
    )

    println("default_route=$(default_adapter.route)")
    println("default_nup=$(default_adapter.nup)")
    println("default_ndn=$(default_adapter.ndn)")
    println("default_psiup0_source=$(default_adapter.psiup0_source)")
    println("default_psidn0_source=$(default_adapter.psidn0_source)")
    println("default_energy=$(default_result.energy)")
    println("explicit_psiup0_source=$(explicit_adapter.psiup0_source)")
    println("explicit_psidn0_source=$(explicit_adapter.psidn0_source)")
    println("explicit_energy=$(explicit_result.energy)")
    println("energy_difference=$(explicit_result.energy - default_result.energy)")
end

main()
