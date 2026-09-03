"""
    CartesianTerminalLowering

Internal terminal-region lowering module for Cartesian routes.

This module owns one step:

    terminal shellification regions -> source CPBs + lowering contracts

Inputs are `CartesianShellification.ShellificationPlan` objects or their
terminal regions. Outputs are metadata-only lowering contracts and compact
summaries. This module does not build coefficient maps, COMX matrices, Lowdin
matrices, retained spaces, final retained units, pair inventories, operator
blocks, Hamiltonians, artifacts, or reports.
"""
module CartesianTerminalLowering

using ..CartesianCPB
using ..CartesianShellification

const CPB = CartesianCPB
const CSH = CartesianShellification

export TerminalLoweringPolicy,
       WhiteLindseyLowering,
       PQSLowering,
       TerminalLoweringContract,
       TerminalLoweringPlan,
       lower_terminal_regions,
       available_contracts,
       selected_contracts,
       contracts,
       summary,
       source_cpbs,
       lowering_kind

# File organization:
#
# policies.jl
#     Terminal lowering policy objects.
#
# contracts.jl
#     TerminalLoweringContract and TerminalLoweringPlan records.
#
# region_contracts.jl
#     Terminal-region to available lowering-contract rules.
#
# selection.jl
#     Route-policy selection of one contract per terminal region.
#
# summaries.jl
#     Compact, stable summaries for reports and tests.
include("policies.jl")
include("contracts.jl")
include("region_contracts.jl")
include("selection.jl")
include("summaries.jl")

end # module CartesianTerminalLowering
