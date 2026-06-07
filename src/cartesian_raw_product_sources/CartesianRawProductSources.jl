"""
    CartesianRawProductSources

Internal metadata-only raw product source module.

This module owns raw product source-box facts:

    source CPB + source-mode dimensions + source-mode ordering
    -> raw product-box source facts

It depends only on `CartesianCPB`. It does not own retained rules, shell
projection or Lowdin cleanup, final retained units, pair blocks, IDA weights,
support-row oracle contraction, Hamiltonians, Ham bundles, exports, or
artifacts.
"""
module CartesianRawProductSources

using ..CartesianCPB

const CPB = CartesianCPB

export RawProductBoxPlan,
       AxisSourceTransformFact,
       raw_product_box_plan,
       source_mode_indices,
       source_mode_count,
       source_mode_dims,
       source_cpb,
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
