using GaussletBases
using LinearAlgebra

Z = 2.0
s = 2.0
lmax = 1

rb = build_basis(RadialBasisSpec(:G10;
    rmax = 5.0,
    mapping = AsinhMapping(c = s / (2Z), s = s),
    reference_spacing = 1.0,
    tails = 6,
    odd_even_kmax = 6,
))

grid = radial_quadrature(rb; accuracy = :medium)
radial_ops = atomic_operators(rb, grid; Z = Z, lmax = lmax)
ida = atomic_ida_operators(radial_ops; lmax = lmax)
problem = atomic_ida_two_electron_problem(ida)

orbital_overlap_error = norm(problem.orbital_overlap - Matrix{Float64}(I, size(problem.orbital_overlap)...), Inf)
two_electron_overlap_error = norm(problem.overlap - Matrix{Float64}(I, size(problem.overlap)...), Inf)
result = lanczos_ground_state(problem; krylovdim = 200, maxiter = 200, tol = 1.0e-7)

println("channel set: ", ida.one_body.channels)
println("radial basis functions per channel: ", length(rb))
println("spatial orbitals: ", length(orbitals(problem)))
println("two-electron product states: ", length(two_electron_states(problem)))
println("Hamiltonian size: ", size(problem.hamiltonian))
println("orbital overlap error: ", orbital_overlap_error)
println("two-electron overlap error: ", two_electron_overlap_error)
println("Lanczos converged: ", result.converged)
println("Lanczos iterations: ", result.iterations)
println("Lanczos residual: ", result.residual)
println("ground-state energy: ", result.value)
