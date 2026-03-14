using GaussletBases

g = Gausslet(:G10; center = 0.0, spacing = 1.0)
x = 0.2
st = stencil(g)

println("gausslet value at ", x, ": ", g(x))
println("named value matches: ", value(g, x))
println("center: ", center(g))
println("integral weight: ", integral_weight(g))
println("stencil length: ", length(st))
println("first primitive type: ", typeof(primitives(st)[1]))
