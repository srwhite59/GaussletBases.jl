"""
    CartesianFinalBasisRealization

Internal terminal/final-basis realization and operator-transfer seam.

This module owns the common terminal-basis object consumed by the staged
Cartesian Hamiltonian producer. It realizes PQS terminal shells, the narrow
White-Lindsey low-order terminal-basis seam, terminal one-body products, IDA
interaction assembly, and Residual-Gaussian augmented-operator compatibility

It does not build raw source operator blocks, RHF results, driver public inputs,
exports, or artifact schemas.
"""
module CartesianFinalBasisRealization

using ..CartesianRawProductSources
using ..CartesianResidualGaussians
using ..CartesianCPB
using LinearAlgebra
using SparseArrays
import ..GaussletBases: _cartesian_flat_index, _cartesian_unflat_index, _nested_axis_lengths,
       _nested_axis_pgdg, _nested_box_support_indices, _nested_doside_1d,
       _nested_face_product, _nested_product_coefficients,
       _nested_projected_q_shell_boundary_comx_product_modes,
       _nested_projected_q_shell_full_sides, gto_overlap_matrix,
       CoulombGaussianExpansion, _ParentGaussianDirectResource,
       _coulomb_expansion_fingerprint, _parent_gaussian_direct_resource,
       _parent_gaussian_direct_value

const CRPS = CartesianRawProductSources
const CRG = CartesianResidualGaussians

include("terminal_face_product_blocks.jl")
include("pqs_terminal_basis_realization.jl")
include("white_lindsey_terminal_basis_realization.jl")
include("pqs_terminal_residual_gto.jl")
include("pqs_terminal_one_body.jl")
include("pqs_terminal_ida.jl")

end # module CartesianFinalBasisRealization
