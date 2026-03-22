using LinearAlgebra
using GaussletBases

Z = 2.0
s = 0.5
lmax = 0

rb = build_basis(RadialBasisSpec(:G10;
    rmax = 8.0,
    mapping = AsinhMapping(c = s / (2Z), s = s),
    reference_spacing = 1.0,
    tails = 6,
    odd_even_kmax = 6,
))

grid = radial_quadrature(rb; accuracy = :high)
radial_ops = atomic_operators(rb, grid; Z = Z, lmax = lmax)
ida = atomic_ida_operators(radial_ops; lmax = lmax)
scf = uhf_scf(ida; nalpha = 1, nbeta = 1, maxiter = 80, damping = 0.25, tol = 1.0e-10)
exact_problem = atomic_ida_two_electron_problem(ida)
exact_energy = ground_state_energy(exact_problem)

println("He-like UHF in the current atomic IDA model")
println("  spatial orbitals: ", length(orbitals(ida)))
println("  alpha electrons: 1")
println("  beta electrons: 1")
println("  converged: ", scf.converged)
println("  iterations: ", scf.iterations)
println("  final UHF energy: ", scf.energy)
println("  exact tiny-model energy: ", exact_energy)
println("  UHF - exact: ", scf.energy - exact_energy)
println("  final density residual: ", scf.residuals[end])
println("  Falpha - Fbeta inf-norm: ", norm(scf.fock_alpha - scf.fock_beta, Inf))
println()
println("Iteration history:")
for iteration in eachindex(scf.energies)
    println(
        "  iter ",
        lpad(iteration, 2),
        ": energy = ",
        scf.energies[iteration],
        ", residual = ",
        scf.residuals[iteration],
    )
end
