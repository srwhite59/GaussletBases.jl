# Route-specific lowering selection.

function _direct_selected_contract(region, ::WhiteLindseyLowering)
    region.region_kind == :complete_shell &&
        return _white_lindsey_complete_shell_contract(region)
    region.region_kind == :central_distorted_product_box &&
        return _distorted_product_contract(region)
    region.region_kind in (:direct_midpoint_slab, :outer_mismatch_slab, :angular_z_extension_slab) &&
        return _thin_slab_contract(region)
    return _direct_terminal_contract(region)
end

function _direct_selected_contract(region, policy::PQSLowering)
    region.region_kind == :complete_shell &&
        return _pqs_complete_shell_contract(region, policy)
    region.region_kind == :central_distorted_product_box &&
        return _distorted_product_contract(region)
    region.region_kind in (:direct_midpoint_slab, :outer_mismatch_slab, :angular_z_extension_slab) &&
        return _thin_slab_contract(region)
    return _direct_terminal_contract(region)
end

function _available_contracts(region, ::WhiteLindseyLowering)
    return available_contracts(region)
end

function _available_contracts(region, policy::PQSLowering)
    return available_contracts(region; pqs_q = policy.q)
end

function _selected_contract(region, policy::TerminalLoweringPolicy)
    return _direct_selected_contract(region, policy)
end

"""
    lower_terminal_regions(shellification_plan, policy)

Select metadata-only lowering contracts for the terminal regions in a
`CartesianShellification.ShellificationPlan`.
"""
function lower_terminal_regions(
    shellification_plan::CartesianShellification.ShellificationPlan,
    policy::TerminalLoweringPolicy = WhiteLindseyLowering();
    enumerate_available_contracts::Bool = false,
    metadata = (;),
)
    regions = CartesianShellification.terminal_regions(shellification_plan)
    selected = TerminalLoweringContract[
        _selected_contract(region, policy) for region in regions]
    available =
        enumerate_available_contracts ?
        TerminalLoweringContract[
            contract
            for region in regions
            for contract in _available_contracts(region, policy)
        ] :
        copy(selected)
    return TerminalLoweringPlan(
        policy,
        available,
        selected,
        _terminal_lowering_summary(policy, regions, available, selected),
        NamedTuple(metadata),
    )
end
