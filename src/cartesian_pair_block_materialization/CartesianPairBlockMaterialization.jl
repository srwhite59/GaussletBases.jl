"""
    CartesianPairBlockMaterialization

Internal owner of PQS source-axis transform facts used by mapped-COMX terminal
realization. Pair planning and pair-block materialization are retired.
"""
module CartesianPairBlockMaterialization

using ..CartesianRawProductSources
using LinearAlgebra

const CRPS = CartesianRawProductSources
const ParentGaussletBases = Base.parentmodule(@__MODULE__)

export pqs_source_axis_transform_facts_from_pgdg_axes

include("pqs_source_axis_transforms.jl")

end # module CartesianPairBlockMaterialization
