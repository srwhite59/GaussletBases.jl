"""
    CartesianCPB

Coordinate Product Box (CPB) geometry module.

This module owns only CPB primitives and CPB-only geometry helpers:
axis-aligned integer product boxes, codimension classification, support counts,
and complete-shell boundary-stratum decomposition.

It deliberately does not depend on route-core, shellification, retained spaces,
pair plans, operators, Hamiltonians, or report plumbing.
"""
module CartesianCPB

export CoordinateProductBox,
       cpb,
       filled_cpb,
       slab_cpb,
       complete_shell_boundary_strata,
       intervals,
       shape,
       codimension,
       role,
       support_count

include(joinpath(@__DIR__, "..", "src", "cartesian_cpb", "coordinate_product_boxes.jl"))

end # module CartesianCPB
