# Route-specific lowering selection.

function _selected_lowering_kind(region, ::WhiteLindseyLowering)
    region.region_kind == :complete_shell && return :white_lindsey_boundary_strata
    region.region_kind == :central_distorted_product_box &&
        return :distorted_product_box_comx
    return _direct_lowering_kind(region)
end

function _selected_lowering_kind(region, ::PQSLowering)
    region.region_kind == :complete_shell && return :pqs_filled_source_cpb
    region.region_kind == :central_distorted_product_box &&
        return :distorted_product_box_comx
    return _direct_lowering_kind(region)
end

function _contract_with_kind(contracts, kind::Symbol)
    matches = Tuple(contract for contract in contracts if contract.lowering_kind == kind)
    length(matches) == 1 ||
        throw(ArgumentError("expected exactly one terminal lowering contract with kind $kind"))
    return only(matches)
end

function _available_contracts(region, ::WhiteLindseyLowering)
    return available_contracts(region)
end

function _available_contracts(region, policy::PQSLowering)
    return available_contracts(region; pqs_q = policy.q)
end

function _selected_contract(region, policy::TerminalLoweringPolicy)
    candidates = _available_contracts(region, policy)
    return _contract_with_kind(candidates, _selected_lowering_kind(region, policy))
end

"""
    lower_terminal_regions(shellification_plan, policy)

Select metadata-only lowering contracts for the terminal regions in a
`CartesianShellification.ShellificationPlan`.
"""
function lower_terminal_regions(
    shellification_plan::CSH.ShellificationPlan,
    policy::TerminalLoweringPolicy = WhiteLindseyLowering();
    metadata = (;),
)
    regions = CSH.terminal_regions(shellification_plan)
    available = Tuple(
        contract
        for region in regions
        for contract in _available_contracts(region, policy)
    )
    selected = Tuple(_selected_contract(region, policy) for region in regions)
    return TerminalLoweringPlan(
        policy,
        available,
        selected,
        _terminal_lowering_summary(policy, regions, available, selected),
        NamedTuple(metadata),
    )
end
