using GaussletBases
using LinearAlgebra

const HFDMRG_PATH = "/Users/srw/Dropbox/codexhome/work/hfdmrg"
if !(HFDMRG_PATH in LOAD_PATH)
    push!(LOAD_PATH, HFDMRG_PATH)
end
using HFDMRG

function build_he_hf_benchmark()
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
    return build_atomic_injected_angular_hf_style_benchmark(radial_ops)
end

function main()
    benchmark = build_he_hf_benchmark()
    adapter = build_atomic_injected_angular_hfdmrg_hf_adapter(benchmark)
    adapter_diag = atomic_injected_angular_hfdmrg_hf_adapter_diagnostics(adapter)
    result = run_atomic_injected_angular_hfdmrg_hf(
        adapter;
        hfmod = HFDMRG,
        maxiter = 40,
        blocksize = 16,
        cutoff = 1.0e-10,
        scf_cutoff = 1.0e-11,
        verbose = false,
    )

    println("route=$(adapter.route)")
    println("basis_dim=$(adapter_diag.basis_dim)")
    println("nup=$(adapter.nup)")
    println("ndn=$(adapter.ndn)")
    println("overlap_identity_error=$(adapter_diag.overlap_identity_error)")
    println("benchmark_full_energy=$(adapter_diag.benchmark_full_energy)")
    println("benchmark_exact_energy=$(adapter_diag.benchmark_exact_energy)")
    println("hfdmrg_energy=$(result.energy)")
end

main()
