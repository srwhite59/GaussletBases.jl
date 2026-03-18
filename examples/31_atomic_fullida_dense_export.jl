using GaussletBases
using JLD2

Z = 2.0
s = 0.2
lmax = 1

rb = build_basis(RadialBasisSpec(:G10;
    rmax = 20.0,
    mapping = AsinhMapping(c = s / (2Z), s = s),
    reference_spacing = 1.0,
    tails = 6,
    odd_even_kmax = 6,
    xgaussians = XGaussian[],
))

grid = radial_quadrature(rb)
radial_ops = atomic_operators(rb, grid; Z = Z, lmax = lmax)
ida = atomic_ida_operators(radial_ops; lmax = lmax)

mktempdir() do dir
    path = joinpath(dir, "atomic_fullida_dense.jld2")
    write_fullida_dense_jld2(
        path,
        ida;
        nelec = 2,
        meta = (Z = Z, example = "31_atomic_fullida_dense_export.jl"),
    )

    jldopen(path, "r") do file
        key_list = sort!(String[key for key in keys(file)])
        println("wrote: ", path)
        println("top-level keys:")
        for key in key_list
            println("  ", key)
        end
        println("H1 size: ", size(file["H1"]))
        println("Vee size: ", size(file["Vee"]))
        println("dims_per_shell: ", Int.(file["dims_per_shell"]))
        println("orders: ", Int.(file["orders"]))
        println("bridge/format: ", String(file["bridge/format"]))
        println("bridge/interaction_model: ", String(file["bridge/interaction_model"]))
    end
end
