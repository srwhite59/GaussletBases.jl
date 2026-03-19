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
    path = joinpath(dir, "atomic_sliced_ham.jld2")
    write_sliced_ham_jld2(
        path,
        ida;
        nelec = 2,
        meta = (Z = Z, example = "32_atomic_sliced_export.jl"),
    )

    jldopen(path, "r") do file
        top_keys = sort!(String[
            key isa AbstractString ? String(key) : join(string.(key), "/") for key in keys(file)
        ])
        H1blocks = file["onebody/H1blocks"]
        Vblocks = file["twobody/Vblocks"]
        println("wrote: ", path)
        println("top-level groups:")
        for key in top_keys
            println("  ", key)
        end
        println("layout/nslices: ", Int(file["layout/nslices"]))
        println("layout/dims: ", Int.(file["layout/dims"]))
        println("ordering/within_slice: ", String(file["ordering/within_slice"]))
        println("onebody/stored: ", String(file["onebody/stored"]))
        println("twobody/convention: ", String(file["twobody/convention"]))
        println("H1blocks outer sizes: ", (length(H1blocks), length(H1blocks[1])))
        println("Vblocks outer sizes: ", (length(Vblocks), length(Vblocks[1])))
        println("nnz(H1blocks[1][1]): ", length(H1blocks[1][1].vals))
        println("nnz(Vblocks[1][2]): ", length(Vblocks[1][2].vals))
        println("meta/interaction_model: ", String(file["meta/interaction_model"]))
    end
end
