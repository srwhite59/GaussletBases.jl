using GaussletBases

map = AsinhMapping(c = 0.15, s = 0.15)
spec = RadialBasisSpec(:G10;
    count = 6,
    mapping = map,
    reference_spacing = 1.0,
    tails = 3,
    odd_even_kmax = 2,
    xgaussians = [XGaussian(alpha = 0.2)],
)

rb = build_basis(spec)
f = rb[2]
grid = radial_quadrature(rb; refine = 24, rmax = 12.0)
diag = basis_diagnostics(rb, grid)

println(rb)
println("basis length: ", length(rb))
println("reference center: ", reference_center(f))
println("physical center: ", center(f))
println("value at 0.2: ", f(0.2))
println("moment center: ", moment_center(f, grid))
println("overlap error: ", diag.overlap_error)
