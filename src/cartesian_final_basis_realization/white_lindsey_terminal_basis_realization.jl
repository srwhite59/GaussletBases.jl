const ParentGaussletBases = Base.parentmodule(@__MODULE__)

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

_wl_axis_index(axis::Symbol) = _terminal_face_axis_index(axis)

_wl_edge_fixed_axes(axis::Symbol) =
    axis === :x ? (:y, :z) :
    axis === :y ? (:x, :z) :
    axis === :z ? (:x, :y) :
    throw(ArgumentError("unsupported White-Lindsey edge free axis"))

function _wl_metadata_value(metadata::NamedTuple, key::Symbol)
    haskey(metadata, key) ||
        throw(ArgumentError("White-Lindsey boundary metadata missing $key"))
    return getfield(metadata, key)
end

function _wl_block_from_product(product, bundles)
    dims = _nested_axis_lengths(bundles)
    indices = Int.(product.support_indices)
    states = NTuple{3,Int}[_cartesian_unflat_index(index, dims) for index in indices]
    coefficients = Matrix{Float64}(product.coefficient_matrix[indices, :])
    return indices, states, coefficients
end

function _wl_boundary_stratum_block(unit, bundles)
    length(unit.source_cpbs) == 1 ||
        throw(ArgumentError("White-Lindsey boundary unit requires one source CPB"))
    cpb = only(unit.source_cpbs)
    intervals = CartesianCPB.intervals(cpb)
    dims = _nested_axis_lengths(bundles)
    q = _wl_metadata_value(unit.metadata, :white_lindsey_retained_count_1d)
    q isa Integer && q > 0 ||
        throw(ArgumentError("White-Lindsey boundary retained count must be positive"))
    side(axis) = ParentGaussletBases._nested_doside_1d(
        _nested_axis_pgdg(bundles, axis), intervals[_wl_axis_index(axis)], Int(q);
        enforce_symmetric_odd = false)
    stratum_kind = _wl_metadata_value(cpb.metadata, :stratum_kind)
    if stratum_kind === :facet_cpb
        fixed_axis = _wl_metadata_value(cpb.metadata, :axis)
        return _terminal_face_product_block(
            cpb,
            bundles;
            normal_axis = fixed_axis,
            fixed_indices = (first(intervals[_wl_axis_index(fixed_axis)]),),
            retained_count = Int(q),
            fixed_side = _wl_metadata_value(cpb.metadata, :side),
        )
    elseif stratum_kind === :edge_cpb
        free_axis = _wl_metadata_value(cpb.metadata, :free_axis)
        fixed_axes = _wl_edge_fixed_axes(free_axis)
        meta_axes = _wl_metadata_value(cpb.metadata, :fixed_axes)
        meta_sides = _wl_metadata_value(cpb.metadata, :sides)
        fixed_sides = ntuple(
            i -> meta_sides[findfirst(==(fixed_axes[i]), meta_axes)],
            2,
        )
        fixed_indices = ntuple(i -> first(intervals[_wl_axis_index(fixed_axes[i])]), 2)
        product = ParentGaussletBases._nested_edge_product(
            free_axis, fixed_sides, side(free_axis), fixed_indices, dims)
        return _wl_block_from_product(product, bundles)
    elseif stratum_kind === :corner_cpb
        indices, states = _wl_terminal_source_support((cpb,), bundles)
        length(indices) == 1 ||
            throw(ArgumentError("White-Lindsey corner CPB must have one support row"))
        return indices, states, reshape([1.0], 1, 1)
    end
    throw(ArgumentError("unsupported White-Lindsey boundary stratum $stratum_kind"))
end

function _wl_validate_coefficients(states, coefficients, overlaps, identity_atol)
    gram = transpose(coefficients) *
           _support_action(states, states, coefficients, overlaps)
    error = _matrix_identity_error(Symmetric((gram + transpose(gram)) ./ 2))
    error <= identity_atol ||
        throw(ArgumentError("White-Lindsey boundary block overlap is not identity"))
    return nothing
end

function _append_white_lindsey_unit!(
    blocks,
    unit,
    contract,
    bundles,
    overlaps,
    seen_support,
    nextcol,
    identity_atol,
)
    unit.source_cpbs == contract.source_cpbs ||
        throw(ArgumentError("White-Lindsey retained unit and transform contract source CPB mismatch"))
    if unit.unit_kind === :direct_cpb_retained_unit
        contract.transform_path === :direct_identity_transform_contract ||
            throw(ArgumentError("White-Lindsey direct unit requires identity transform"))
        indices, states = _wl_terminal_source_support(unit.source_cpbs, bundles)
        nextcol = _append_direct!(
            blocks,
            (; support_record = (; unit_key = unit.unit_key, support_indices = indices, support_states = states)),
            bundles,
            overlaps,
            nextcol,
            identity_atol,
        )
    elseif unit.unit_kind === :white_lindsey_boundary_stratum_retained_unit
        contract.transform_path === :white_lindsey_boundary_stratum_product_contract ||
            throw(ArgumentError("White-Lindsey boundary unit requires product transform"))
        indices, states, coefficients = _wl_boundary_stratum_block(unit, bundles)
        _wl_validate_coefficients(states, coefficients, overlaps, identity_atol)
        nextcol = _push_block!(blocks, unit.unit_key, indices, states, coefficients, nextcol)
    elseif unit.unit_kind === :compact_thin_slab_retained_unit
        contract.transform_path === :compact_thin_slab_face_product_contract ||
            throw(ArgumentError("compact thin slab unit requires face-product transform"))
        indices, states, coefficients =
            _terminal_compact_thin_slab_block(contract.source_cpbs, contract.metadata, bundles)
        _wl_validate_coefficients(states, coefficients, overlaps, identity_atol)
        nextcol = _push_block!(blocks, unit.unit_key, indices, states, coefficients, nextcol)
    else
        throw(ArgumentError("unsupported White-Lindsey terminal retained unit"))
    end
    _validate_block_support!(
        last(blocks),
        (; support_indices = indices, support_states = states),
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
        nextcol = _append_white_lindsey_unit!(
            blocks, unit, contract, bundles, overlaps, seen_support, nextcol,
            Float64(identity_atol))
    end
    return CartesianTerminalBasisRealization(blocks, nextcol - 1, 0.0)
end
