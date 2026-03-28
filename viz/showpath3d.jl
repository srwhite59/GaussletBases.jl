#!/usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__)

include(joinpath(@__DIR__, "show_bond_aligned_diatomic_points3d.jl"))

main_show_ordering_path3d(ARGS)
