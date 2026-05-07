#!/usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__)

if isempty(ARGS) || ARGS[1] != "--describe"
    using GLMakie
end

include(joinpath(@__DIR__, "show_bond_aligned_diatomic_points3d.jl"))

main_show_bond_aligned_diatomic_points3d(ARGS)
