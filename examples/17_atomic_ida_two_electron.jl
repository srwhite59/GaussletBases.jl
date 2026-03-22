using GaussletBases
using LinearAlgebra

Z = 2.0
s = 0.5
lmax = 0

rb = build_basis(RadialBasisSpec(:G10;
    rmax = 8.0,
    mapping = AsinhMapping(c = s / (2Z), s = s),
))

grid = radial_quadrature(rb; accuracy = :high)
radial_ops = atomic_operators(rb, grid; Z = Z, lmax = lmax)
ida = atomic_ida_operators(radial_ops; lmax = lmax)
problem = atomic_ida_two_electron_problem(ida)
orbital_overlap_error = norm(problem.orbital_overlap - Matrix{Float64}(I, size(problem.orbital_overlap)...), Inf)
two_electron_overlap_error = norm(problem.overlap - Matrix{Float64}(I, size(problem.overlap)...), Inf)
two_body_offdiag_error = norm(problem.two_body - Diagonal(diag(problem.two_body)), Inf)

println("channel set: ", ida.one_body.channels)
println("radial basis functions per channel: ", length(rb))
println("spatial orbitals: ", length(orbitals(ida)))
println("two-electron product states: ", length(two_electron_states(problem)))
println("Hamiltonian size: ", size(problem.hamiltonian))
println("orbital overlap error: ", orbital_overlap_error)
println("two-electron overlap error: ", two_electron_overlap_error)
println("two-body off-diagonal norm: ", two_body_offdiag_error)
println("ground-state energy: ", ground_state_energy(problem))
println("first three product states:")
for state in two_electron_states(problem)[1:3]
    println("  ", state)
end
