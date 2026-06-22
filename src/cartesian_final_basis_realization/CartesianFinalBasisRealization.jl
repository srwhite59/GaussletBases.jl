"""
    CartesianFinalBasisRealization

Internal final-basis realization and operator-transfer seam for PQS work.

This module consumes raw product-source facts, PQS source-mode boundary rules,
caller-supplied shell-realization data, and retained-boundary one-body
operators. It produces shell-realized final-basis data, overlap/isometry
diagnostics, oracle shell-support projections, and retained-boundary
one-body transfers into the final basis. For the complete core/shell route it
also owns the narrow localized IDA density-interaction seam used to audit
final-basis orbital consumption before RHF.

It does not build raw source operator blocks, pair-block materialization
records, by-center nuclear CPBM result blocks, RHF results, driver wiring,
exports, or artifacts.
"""
module CartesianFinalBasisRealization

using ..CartesianRawProductSources
using LinearAlgebra
import ..GaussletBases: _cartesian_unflat_index, _nested_axis_lengths,
       _nested_axis_pgdg, _nested_box_support_indices, _nested_product_coefficients,
       _nested_projected_q_shell_boundary_comx_product_modes,
       _nested_projected_q_shell_full_sides, gto_overlap_matrix

const CRPS = CartesianRawProductSources

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

include("pqs_terminal_basis_realization.jl")
include("pqs_terminal_residual_gto.jl")
include("pqs_terminal_one_body.jl")
include("pqs_terminal_ida.jl")

end # module CartesianFinalBasisRealization
