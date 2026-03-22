using GaussletBases
using JLD2
using LinearAlgebra

Z = 2.0
s = 0.2
lmax = 1

rb = build_basis(RadialBasisSpec(:G10;
    rmax = 20.0,
    mapping = AsinhMapping(c = s / (2Z), s = s),
    # Keep the export demo focused on file structure rather than xgaussian tuning.
    xgaussian_count = 0,
))

grid = radial_quadrature(rb; accuracy = :medium)
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
        H1 = file["H1"]
        Vee = file["Vee"]
        dims_per_shell = Int.(file["dims_per_shell"])
        orders = Int.(file["orders"])
        hermiticity_error = norm(H1 - H1', Inf)
        max_vee = maximum(abs, Vee)
        sum(dims_per_shell) == size(H1, 1) ||
            error("exported dims_per_shell does not sum to the H1 dimension")
        size(H1, 1) == size(H1, 2) ||
            error("exported H1 is not square")
        all(size(Vee) .== size(H1, 1)) ||
            error("exported Vee tensor does not match the H1 dimension")

        println("wrote: ", path)
        println("top-level keys:")
        for key in key_list
            println("  ", key)
        end
        println("H1 size: ", size(H1))
        println("Vee size: ", size(Vee))
        println("dims_per_shell: ", dims_per_shell)
        println("orders: ", orders)
        println("bridge/format: ", String(file["bridge/format"]))
        println("bridge/interaction_model: ", String(file["bridge/interaction_model"]))
        println("sanity checks:")
        println("  H1 Hermiticity error: ", hermiticity_error)
        println("  max |Vee|: ", max_vee)
        println("  shell dimension sum matches H1: ", true)
    end
end
