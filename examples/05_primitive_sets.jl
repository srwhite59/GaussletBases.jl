using LinearAlgebra
using GaussletBases

g = Gausslet(:G10; center = 0.0, spacing = 1.0)
plain_set = PrimitiveSet1D(primitives(stencil(g)); name = :plain_gausslet_stencil)

S_plain = overlap_matrix(plain_set)
T_plain = kinetic_matrix(plain_set)

map = AsinhMapping(c = 0.15, s = 0.15)
distorted_set = PrimitiveSet1D(
    [Distorted(primitive, map) for primitive in primitives(plain_set)];
    name = :distorted_gausslet_stencil,
)
S_distorted = overlap_matrix(distorted_set)

ub = build_basis(UniformBasisSpec(:G10; xmin = -1.0, xmax = 1.0, spacing = 1.0))
meta = basis_metadata(ub)

println(plain_set)
println("plain overlap size: ", size(S_plain))
println("plain kinetic size: ", size(T_plain))
println("distorted overlap size: ", size(S_distorted))
println("basis metadata kind: ", meta.basis_kind)
println("basis metadata primitive count: ", length(meta.primitive_set))
println("basis metadata coefficient matrix size: ", size(meta.coefficient_matrix))
println("plain overlap is symmetric: ", isapprox(S_plain, transpose(S_plain), atol = 1.0e-12, rtol = 1.0e-12))
