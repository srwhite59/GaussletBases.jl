function _wl_terminal_cpb_support!(indices, states, cpb, bundles)
    intervals = CartesianCPB.intervals(cpb)
    dims = _nested_axis_lengths(bundles)
    cpb_indices = _nested_box_support_indices(intervals[1], intervals[2], intervals[3], dims)
    append!(indices, cpb_indices)
    append!(states, NTuple{3,Int}[_cartesian_unflat_index(index, dims) for index in cpb_indices])
    return indices, states
end

function _wl_terminal_source_support(source_cpbs, bundles)
    indices = Int[]
    states = NTuple{3,Int}[]
    for cpb in source_cpbs
        _wl_terminal_cpb_support!(indices, states, cpb, bundles)
    end
    return indices, states
end

function _wl_terminal_support_record(unit, indices, states)
    return (;
        support_indices = indices,
        support_states = states,
    )
end

function _append_white_lindsey_identity_unit!(
    blocks,
    unit,
    contract,
    bundles,
    overlaps,
    seen_support,
    nextcol,
    identity_atol,
)
    is_direct =
        unit.unit_kind === :direct_cpb_retained_unit &&
        contract.transform_path === :direct_identity_transform_contract
    is_stratum =
        unit.unit_kind === :white_lindsey_boundary_stratum_retained_unit &&
        contract.transform_path === :white_lindsey_boundary_stratum_product_contract
    (is_direct || is_stratum) ||
        throw(ArgumentError("White-Lindsey terminal basis requires direct or boundary-stratum units"))
    unit.source_cpbs == contract.source_cpbs ||
        throw(ArgumentError("White-Lindsey retained unit and transform contract source CPB mismatch"))
    indices, states = _wl_terminal_source_support(unit.source_cpbs, bundles)
    nextcol = _append_direct!(
        blocks,
        (; support_record = (; unit_key = unit.unit_key, support_indices = indices, support_states = states)),
        bundles,
        overlaps,
        nextcol,
        identity_atol,
    )
    _validate_block_support!(
        last(blocks),
        _wl_terminal_support_record(unit, indices, states),
        seen_support,
    )
    return nextcol
end

function white_lindsey_terminal_basis_realization(
    retained_units,
    transform_contracts,
    bundles;
    identity_atol::Real = 1.0e-8,
)
    contracts = Dict(contract.unit_key => contract for contract in transform_contracts)
    overlaps = (_nested_axis_pgdg(bundles, :x).overlap,
        _nested_axis_pgdg(bundles, :y).overlap,
        _nested_axis_pgdg(bundles, :z).overlap)
    blocks = CartesianTerminalBasisBlock[]
    seen_support = Set{Int}()
    nextcol = 1
    for unit in retained_units
        contract = get(contracts, unit.unit_key, nothing)
        isnothing(contract) &&
            throw(ArgumentError("missing White-Lindsey terminal transform contract"))
        nextcol = _append_white_lindsey_identity_unit!(
            blocks, unit, contract, bundles, overlaps, seen_support, nextcol,
            Float64(identity_atol))
    end
    return CartesianTerminalBasisRealization(blocks, nextcol - 1, 0.0)
end
