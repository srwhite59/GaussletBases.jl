using GaussletBases

map = AsinhMapping(c = 0.15, s = 0.15)
spec = RadialBasisSpec(:G10;
    count = 9,
    mapping = map,
)

rb = build_basis(spec)
f = rb[2]
grid = radial_quadrature(rb)
diag = basis_diagnostics(rb)

println(rb)
println("basis length: ", length(rb))
println("reference center: ", reference_center(f))
println("physical center: ", center(f))
println("value at 0.2: ", f(0.2))
println("moment center: ", moment_center(f, grid))
println("overlap error: ", diag.overlap_error)
println("aggregate center mismatch D: ", diag.D)
println("largest center mismatch: ", maximum(abs, diag.center_mismatches))
