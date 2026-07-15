"""
    CartesianFinalBasisRealization

Internal terminal/final-basis realization and operator-transfer seam.

This module owns the common terminal-basis object consumed by the staged
Cartesian Hamiltonian producer. It realizes PQS terminal shells, the narrow
White-Lindsey low-order terminal-basis seam, terminal one-body products, IDA
interaction assembly, and Residual-Gaussian augmented-operator compatibility
helpers. It also keeps older PQS shell-support projection and complete
core/shell final-basis utilities used by live reference and migration paths.

It does not build raw source operator blocks, pair-block materialization
records, by-center nuclear CPBM result blocks, RHF results, driver public
inputs, exports, or artifact schemas.
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

export pqs_source_shell_realization_final_basis,
       pqs_source_shell_projected_one_body_matrix,
       pqs_source_shell_final_one_body_from_boundary_matrix

# pqs_source_shell_final_basis.jl
#     PQS shell/support realization, Lowdin final basis diagnostics,
#     shell-support oracle projection, and retained-boundary one-body transfer.
include("pqs_source_shell_final_basis.jl")

# pqs_complete_core_shell_final_basis.jl
#     Direct-core plus surrounding-shell final-basis realization, one-body
#     transfer, H1 solve, and localized IDA density-interaction seam.
include("pqs_complete_core_shell_final_basis.jl")

include("terminal_face_product_blocks.jl")
include("pqs_terminal_basis_realization.jl")
include("white_lindsey_terminal_basis_realization.jl")
include("pqs_terminal_residual_gto.jl")
include("pqs_terminal_one_body.jl")
include("pqs_terminal_ida.jl")

end # module CartesianFinalBasisRealization
