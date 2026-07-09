module CartesianReferenceDensity

using JLD2
using LinearAlgebra
using SHA
using SpecialFunctions

const _GB_PARENT = parentmodule(@__MODULE__)

include("atomic_hf_reference_packets.jl")
include("screened_hartree_correction.jl")

end # module CartesianReferenceDensity
