"""
    CartesianRawProductSources

Internal metadata-only raw product source module.

This module owns raw product source-box facts:

source CPB + source-mode dimensions + source-mode ordering
    -> raw product-box source facts
    -> optional PQS source-mode retained boundary rule

It depends only on `CartesianCPB`. It owns the source-mode boundary selector
metadata needed for the first PQS source-box route. It does not own shell
projection or Lowdin cleanup, final retained units, pair blocks, IDA weights,
support-row oracle contraction, Hamiltonians, Ham bundles, exports, or
artifacts.

There is deliberately no policy object yet. The module currently supports one
source-mode ordering, `:x_major_y_major_z_fast`; a policy type should be added
only when there is a real second ordering or source-mode convention to select.
"""
module CartesianRawProductSources

using ..CartesianCPB

const CPB = CartesianCPB

export RawProductBoxPlan,
       AxisSourceTransformFact,
       PQSBoundaryProductModeRetainedRule,
       raw_product_box_plan,
       pqs_boundary_product_mode_retained_rule,
       source_mode_indices,
       source_mode_count,
       source_mode_dims,
       source_cpb,
       retained_mode_indices,
       retained_column_indices,
       axis_transform_facts,
       summary,
       unavailable_summary

# File organization:
#
# records.jl
#     Metadata-only raw product source records.
#
# source_mode_indices.jl
#     Deterministic source-mode ordering helpers.
#
# axis_transform_facts.jl
#     Metadata-only axis transform fact defaults.
#
# summaries.jl
#     Compact summaries for tests and reports.
include("records.jl")
include("source_mode_indices.jl")
include("axis_transform_facts.jl")
include("summaries.jl")

end # module CartesianRawProductSources
