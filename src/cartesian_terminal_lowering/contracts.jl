# Terminal lowering contract records.

"""
    TerminalLoweringContract

Metadata-only contract describing how one terminal shellification region will
be lowered into source CPBs and a construction recipe. It is not a transform,
retained space, final retained unit, pair block, or Hamiltonian object.
"""
struct TerminalLoweringContract
    contract_key::Symbol
    terminal_region_key::Symbol
    terminal_region_role::Symbol
    terminal_region_kind::Symbol
    lowering_kind::Symbol
    owned_support::Any
    source_cpbs::Tuple{Vararg{CPB.CoordinateProductBox}}
    retained_rule::Symbol
    realization_rule::Union{Symbol,Nothing}
    final_unit_granularity::Symbol
    materialized::Bool
    metadata::NamedTuple
end

"""
    TerminalLoweringPlan

Selected terminal lowering plan for one shellification plan and lowering
policy. `available_contracts` may include alternative contracts; `contracts`
contains the selected route contracts downstream construction should consume.
"""
struct TerminalLoweringPlan
    policy::TerminalLoweringPolicy
    available_contracts::Tuple{Vararg{TerminalLoweringContract}}
    contracts::Tuple{Vararg{TerminalLoweringContract}}
    summary::NamedTuple
    metadata::NamedTuple
end

available_contracts(plan::TerminalLoweringPlan) = plan.available_contracts
selected_contracts(plan::TerminalLoweringPlan) = plan.contracts
contracts(plan::TerminalLoweringPlan) = plan.contracts
summary(plan::TerminalLoweringPlan) = plan.summary

source_cpbs(contract::TerminalLoweringContract) = contract.source_cpbs
lowering_kind(contract::TerminalLoweringContract) = contract.lowering_kind

function _terminal_lowering_contract(;
    contract_key::Symbol,
    terminal_region,
    lowering_kind::Symbol,
    source_cpbs,
    retained_rule::Symbol,
    realization_rule,
    final_unit_granularity::Symbol,
    metadata = (;),
)
    cpb_tuple = Tuple(source_cpbs)
    isempty(cpb_tuple) &&
        throw(ArgumentError("terminal lowering contract requires at least one source CPB"))
    all(cpb -> cpb isa CPB.CoordinateProductBox, cpb_tuple) ||
        throw(ArgumentError("terminal lowering source_cpbs must be CoordinateProductBox objects"))

    return TerminalLoweringContract(
        contract_key,
        terminal_region.key,
        terminal_region.role,
        terminal_region.region_kind,
        lowering_kind,
        terminal_region.owned_support,
        cpb_tuple,
        retained_rule,
        realization_rule,
        final_unit_granularity,
        false,
        NamedTuple(metadata),
    )
end
