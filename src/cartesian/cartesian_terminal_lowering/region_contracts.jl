# Region-to-contract lowering rules.

function _terminal_source_role(region, suffix::Symbol)
    return Symbol(String(region.key), "_", String(suffix))
end

function _thin_slab_contract(region)
    raw = region.raw_region
    axis = get(raw.metadata, :slab_normal_axis, get(raw.metadata, :bond_axis, nothing))
    axis_index = findfirst(==(axis), (:x, :y, :z))
    thickness = length(raw.outer_box[axis_index])
    q = get(raw.metadata, :thin_slab_retained_count_1d, nothing)
    q isa Integer && q > 0 ||
        throw(ArgumentError("thin-slab lowering requires public ns retained count"))
    thickness <= q ||
        throw(ArgumentError("thin-slab thickness exceeds public ns retained count"))
    source_cpb = CartesianCPB.cpb(
        raw.outer_box;
        role = _terminal_source_role(region, :thin_slab_source_cpb),
        metadata = (; terminal_region_key = region.key, slab_normal_axis = axis),
    )
    return _terminal_lowering_contract(
        contract_key = Symbol(String(region.key), "_compact_thin_slab_product"),
        terminal_region = region,
        lowering_kind = :compact_thin_slab_product_cpb,
        source_cpbs = (source_cpb,),
        retained_rule = :compact_thin_slab_face_product,
        realization_rule = :compact_thin_slab_face_stack,
        final_unit_granularity = :one_terminal_region,
        metadata = merge(raw.metadata, (; slab_normal_axis = axis, slab_thickness = thickness)),
    )
end

function _direct_terminal_contract(region)
    region.region_kind in (:direct_core, :direct_atom_contact_core) ||
        throw(ArgumentError("not a direct terminal lowering region: $(region.region_kind)"))
    source_cpb = CartesianCPB.cpb(
        region.raw_region.outer_box;
        role = _terminal_source_role(region, :direct_source_cpb),
        metadata = (;
            terminal_region_key = region.key,
            terminal_region_kind = region.region_kind,
        ),
    )
    return _terminal_lowering_contract(
        contract_key = Symbol(String(region.key), "_direct_core_identity_cpb"),
        terminal_region = region,
        lowering_kind = :direct_core_identity_cpb,
        source_cpbs = (source_cpb,),
        retained_rule = :direct_source_modes,
        realization_rule = :direct_or_trivial_embedding,
        final_unit_granularity = :one_terminal_region,
        metadata = (;
            identity_like = true,
            source_cpb_equals_owned_support = true,
        ),
    )
end

function _white_lindsey_complete_shell_contract(region)
    raw = region.raw_region
    isnothing(raw.inner_exclusion_box) &&
        throw(ArgumentError("complete-shell lowering requires an inner exclusion box"))

    outer = CartesianCPB.filled_cpb(
        raw.outer_box...;
        role = _terminal_source_role(region, :outer_box),
        metadata = (; terminal_region_key = region.key),
    )
    inner = CartesianCPB.filled_cpb(
        raw.inner_exclusion_box...;
        role = _terminal_source_role(region, :inner_box),
        metadata = (; terminal_region_key = region.key),
    )
    strata = CartesianCPB.complete_shell_boundary_strata(outer, inner)
    return _terminal_lowering_contract(
        contract_key = Symbol(String(region.key), "_white_lindsey_boundary_strata"),
        terminal_region = region,
        lowering_kind = :white_lindsey_boundary_strata,
        source_cpbs = strata.all_strata,
        retained_rule = :white_lindsey_boundary_stratum_product,
        realization_rule = :white_lindsey_boundary_stratum_product,
        final_unit_granularity = :boundary_stratum_children,
        metadata = (;
            facet_count = length(strata.facets),
            edge_count = length(strata.edges),
            corner_count = length(strata.corners),
            total_source_cpb_count = length(strata.all_strata),
            shell_support_count = strata.shell_support_count,
            stratum_support_count = strata.stratum_support_count,
        ),
    )
end

function _pqs_complete_shell_contract(region, policy::PQSLowering)
    source = CartesianCPB.filled_cpb(
        region.raw_region.outer_box...;
        role = _terminal_source_role(region, :pqs_filled_source_cpb),
        metadata = (; terminal_region_key = region.key),
    )
    return _terminal_lowering_contract(
        contract_key = Symbol(String(region.key), "_pqs_filled_source_cpb"),
        terminal_region = region,
        lowering_kind = :pqs_filled_source_cpb,
        source_cpbs = (source,),
        retained_rule = :pqs_boundary_comx_product_modes,
        realization_rule = :shell_projection_lowdin,
        final_unit_granularity = :one_terminal_region,
        metadata = (;
            q = policy.q,
            source_mode_shape = ntuple(_ -> policy.q, 3),
            source_box_shape = CartesianCPB.shape(source),
            face_edge_corner_decomposition_required = false,
        ),
    )
end

function _pqs_complete_shell_contract(region; q = nothing)
    isnothing(q) && return _pqs_unparameterized_complete_shell_contract(region)
    return _pqs_complete_shell_contract(region, PQSLowering(q = q))
end

function _pqs_unparameterized_complete_shell_contract(region)
    source = CartesianCPB.filled_cpb(
        region.raw_region.outer_box...;
        role = _terminal_source_role(region, :pqs_filled_source_cpb),
        metadata = (; terminal_region_key = region.key),
    )
    return _terminal_lowering_contract(
        contract_key = Symbol(String(region.key), "_pqs_filled_source_cpb"),
        terminal_region = region,
        lowering_kind = :pqs_filled_source_cpb,
        source_cpbs = (source,),
        retained_rule = :pqs_boundary_comx_product_modes,
        realization_rule = :shell_projection_lowdin,
        final_unit_granularity = :one_terminal_region,
        metadata = (;
            q = nothing,
            parameter_status = :available_but_unparameterized,
            source_mode_shape = nothing,
            source_box_shape = CartesianCPB.shape(source),
            face_edge_corner_decomposition_required = false,
        ),
    )
end

function _distorted_product_contract(region)
    raw = region.raw_region
    source = CartesianCPB.cpb(
        raw.outer_box;
        role = _terminal_source_role(region, :distorted_product_source_cpb),
        metadata = (; terminal_region_key = region.key),
    )
    return _terminal_lowering_contract(
        contract_key = Symbol(String(region.key), "_distorted_product_box_comx"),
        terminal_region = region,
        lowering_kind = :distorted_product_box_comx,
        source_cpbs = (source,),
        retained_rule = :distorted_product_comx_all_axes,
        realization_rule = :distorted_product_realization_planned,
        final_unit_granularity = :one_terminal_region,
        metadata = (;
            q = raw.metadata.q,
            L = raw.metadata.L,
            source_mode_shape = raw.metadata.source_mode_shape,
            aspect_ratio = raw.metadata.aspect_ratio,
            identity_like = false,
        ),
    )
end

function available_contracts(region::CartesianShellification.TerminalRegion; pqs_q = nothing)
    if region.region_kind in (
        :direct_core,
        :direct_atom_contact_core,
    )
        return (_direct_terminal_contract(region),)
    end
    if region.region_kind in (:direct_midpoint_slab, :outer_mismatch_slab, :angular_z_extension_slab)
        return (_thin_slab_contract(region),)
    end
    if region.region_kind == :complete_shell
        return (
            _white_lindsey_complete_shell_contract(region),
            _pqs_complete_shell_contract(region; q = pqs_q),
        )
    end
    if region.region_kind == :central_distorted_product_box
        return (_distorted_product_contract(region),)
    end
    throw(ArgumentError("unsupported terminal region kind $(region.region_kind)"))
end
