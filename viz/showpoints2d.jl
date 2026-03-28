#!/usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__)

include(joinpath(@__DIR__, "render_bond_aligned_diatomic_projection.jl"))

main_render_bond_aligned_diatomic_projection(ARGS)
