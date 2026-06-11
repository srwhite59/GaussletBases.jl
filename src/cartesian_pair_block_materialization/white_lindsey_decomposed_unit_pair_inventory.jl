# Route-owned decomposed White-Lindsey unit-pair inventory.
#
# This validates retained-unit pair metadata that later acceptance assembly can
# consume. For shellification-derived retained units it materializes retained
# dimensions and column ranges from the White-Lindsey retained-count policy. It
# does not build local blocks, assemble matrices, or call the legacy direct
# Cartesian operator path.

const _WHITE_LINDSEY_ACCEPTANCE_ONE_BODY_TERMS = (
    :overlap,
    :kinetic,
    :electron_nuclear_by_center,
)
const _WHITE_LINDSEY_COMPACT_PAIR_SUMMARY_LIMIT = 1024

function white_lindsey_shellification_decomposed_unit_pair_inventory(
    parent_axes::NTuple{3,<:AbstractVector},
    nuclear_positions;
    shellification_policy = CSH.OneCenterShellification(core_side = 5, q = 5),
    lowering_policy = CTL.WhiteLindseyLowering(),
    supported_terms = _WHITE_LINDSEY_ACCEPTANCE_ONE_BODY_TERMS,
    metadata = (;),
    parent_axis_counts = ntuple(index -> length(parent_axes[index]), 3),
    parent_axis_bundle_object = nothing,
)
    metadata_tuple = NamedTuple(metadata)
    try
        shellification_plan = CSH.shellify(
            parent_axes,
            nuclear_positions;
            policy = shellification_policy,
        )
        lowering_plan =
            CTL.lower_terminal_regions(shellification_plan, lowering_policy)
        retained_unit_plan = CRU.retained_unit_plan(lowering_plan)
        inventory = _white_lindsey_decomposed_unit_pair_inventory_from_shellification_retained_plan(
            retained_unit_plan;
            supported_terms,
            metadata = merge(
                metadata_tuple,
                (;
                    parent_axis_counts,
                    parent_axis_bundle_object,
                    shellification_policy =
                        Symbol(nameof(typeof(shellification_policy))),
                    lowering_policy = Symbol(nameof(typeof(lowering_policy))),
                    driver_path_source =
                        :shellification_lowering_retained_units_lightweight_unit_pairs,
                ),
            ),
        )
        return (;
            object_kind =
                :white_lindsey_shellification_decomposed_unit_pair_inventory,
            status = inventory.status,
            blocker = inventory.blocker,
            shellification_plan,
            lowering_plan,
            retained_unit_plan,
            unit_pair_plan = nothing,
            unit_pair_source_kind = :upper_triangular_unit_index_table,
            inventory,
            shellification_backed_decomposed_wl_inventory =
                inventory.source_kind ===
                :cartesian_shellification_retained_unit_pair_plan,
            low_order_materialized_seed_inventory_used = false,
            full_parent_window_cpb_used = false,
            direct_cartesian_product_assembly_used = false,
            ordinary_cartesian_ida_operators_used = false,
        )
    catch err
        return (;
            object_kind =
                :white_lindsey_shellification_decomposed_unit_pair_inventory,
            status =
                :blocked_white_lindsey_shellification_decomposed_unit_pair_inventory,
            blocker = :shellification_decomposed_inventory_construction_failed,
            error = sprint(showerror, err),
            shellification_plan = nothing,
            lowering_plan = nothing,
            retained_unit_plan = nothing,
            unit_pair_plan = nothing,
            unit_pair_source_kind = :unavailable,
            inventory = _white_lindsey_decomposed_unit_pair_inventory_result(
                :blocked_white_lindsey_decomposed_unit_pair_inventory,
                :shellification_decomposed_inventory_construction_failed,
                ();
                source_kind = :cartesian_shellification_retained_unit_pair_plan,
                supported_terms = Tuple(supported_terms),
                metadata = metadata_tuple,
            ),
            shellification_backed_decomposed_wl_inventory = false,
            low_order_materialized_seed_inventory_used = false,
            full_parent_window_cpb_used = false,
            direct_cartesian_product_assembly_used = false,
            ordinary_cartesian_ida_operators_used = false,
        )
    end
end

function white_lindsey_decomposed_unit_pair_inventory(
    plan::CUP.UnitPairPlan;
    supported_terms = _WHITE_LINDSEY_ACCEPTANCE_ONE_BODY_TERMS,
    metadata = (;),
)
    materialized_pairs =
        _white_lindsey_shellification_unit_pair_records_with_ranges(
            plan;
            metadata = NamedTuple(metadata),
        )
    if materialized_pairs.status ===
       :available_shellification_decomposed_wl_unit_pair_records
        return white_lindsey_decomposed_unit_pair_inventory(
            materialized_pairs.unit_pairs;
            source_kind = :cartesian_shellification_retained_unit_pair_plan,
            supported_terms,
            metadata = merge(
                NamedTuple(metadata),
                (;
                    shellification_unit_pair_plan_source =
                        :cartesian_unit_pair_plan,
                    retained_range_source =
                        :shellification_source_cpb_support_order,
                ),
            ),
        )
    end
    return white_lindsey_decomposed_unit_pair_inventory(
        CUP.unit_pairs(plan);
        source_kind = :cartesian_unit_pair_plan,
        supported_terms,
        metadata,
    )
end

function _white_lindsey_decomposed_unit_pair_inventory_from_shellification_retained_plan(
    retained_unit_plan::CRU.RetainedUnitPlan;
    supported_terms = _WHITE_LINDSEY_ACCEPTANCE_ONE_BODY_TERMS,
    metadata = (;),
)
    metadata_tuple = NamedTuple(metadata)
    materialized_units =
        _white_lindsey_shellification_retained_units_with_ranges(
            retained_unit_plan;
            metadata = metadata_tuple,
        )
    if materialized_units.status !==
       :available_shellification_decomposed_wl_retained_units
        return _white_lindsey_decomposed_unit_pair_inventory_result(
            :blocked_white_lindsey_decomposed_unit_pair_inventory,
            materialized_units.blocker,
            ();
            source_kind = :cartesian_shellification_retained_unit_pair_plan,
            supported_terms = Tuple(supported_terms),
            metadata = metadata_tuple,
        )
    end
    pairs = WhiteLindseyUnitPairIndexTable(
        materialized_units.retained_units,
        (; retained_range_source = :shellification_source_cpb_support_order),
    )
    return _white_lindsey_decomposed_unit_pair_inventory_result(
        :available_white_lindsey_decomposed_unit_pair_inventory,
        nothing,
        pairs;
        source_kind = :cartesian_shellification_retained_unit_pair_plan,
        supported_terms = Tuple(supported_terms),
        metadata = merge(
            metadata_tuple,
            (;
                shellification_unit_pair_plan_source =
                    :upper_triangular_unit_index_table,
                retained_range_source =
                    :shellification_source_cpb_support_order,
                rich_unit_pair_records_stored = false,
                route_core_pair_sidecars_duplicated = false,
            ),
        ),
        units = materialized_units.retained_units,
    )
end

function white_lindsey_decomposed_unit_pair_inventory(
    source::NamedTuple;
    supported_terms = _WHITE_LINDSEY_ACCEPTANCE_ONE_BODY_TERMS,
    metadata = (;),
)
    if _white_lindsey_seed_source_available(source)
        metadata_tuple = NamedTuple(metadata)
        pairs = _white_lindsey_seed_decomposed_unit_pairs(
            source;
            metadata = metadata_tuple,
        )
        return white_lindsey_decomposed_unit_pair_inventory(
            pairs;
            source_kind = :white_lindsey_low_order_materialized_seed_ranges,
            supported_terms,
            metadata = merge(
                metadata_tuple,
                (;
                    seed_inventory_source =
                        :white_lindsey_low_order_materialized_seed_inventory,
                    fixed_block_operator_matrices_used = false,
                ),
            ),
        )
    end
    return _white_lindsey_decomposed_unit_pair_inventory_result(
        :blocked_white_lindsey_decomposed_unit_pair_inventory,
        :unsupported_decomposed_wl_unit_pair_inventory_source,
        ();
        source_kind = :unsupported_named_tuple_source,
        supported_terms = Tuple(supported_terms),
        metadata,
    )
end

function white_lindsey_decomposed_unit_pair_inventory(
    pairs::AbstractVector{<:CUP.UnitPairRecord};
    kwargs...,
)
    return _white_lindsey_decomposed_unit_pair_inventory_from_pairs(
        pairs;
        kwargs...,
    )
end

function white_lindsey_decomposed_unit_pair_inventory(
    pairs::Tuple{Vararg{CUP.UnitPairRecord}};
    kwargs...,
)
    return _white_lindsey_decomposed_unit_pair_inventory_from_pairs(
        pairs;
        kwargs...,
    )
end

function _white_lindsey_decomposed_unit_pair_inventory_from_pairs(
    pairs;
    source_kind::Symbol = :unit_pair_records,
    supported_terms = _WHITE_LINDSEY_ACCEPTANCE_ONE_BODY_TERMS,
    metadata = (;),
)
    term_tuple = Tuple(supported_terms)
    if isempty(pairs)
        return _white_lindsey_decomposed_unit_pair_inventory_result(
            :blocked_white_lindsey_decomposed_unit_pair_inventory,
            :empty_decomposed_wl_unit_pair_inventory,
            ();
            source_kind,
            supported_terms = term_tuple,
            metadata,
        )
    end

    units = _white_lindsey_decomposed_inventory_units(pairs)
    blocker = _white_lindsey_decomposed_inventory_blocker(units)
    status =
        isnothing(blocker) ?
        :available_white_lindsey_decomposed_unit_pair_inventory :
        :blocked_white_lindsey_decomposed_unit_pair_inventory

    return _white_lindsey_decomposed_unit_pair_inventory_result(
        status,
        blocker,
        pairs;
        source_kind,
        supported_terms = term_tuple,
        metadata,
        units,
    )
end

function white_lindsey_decomposed_unit_pair_inventory(
    source;
    source_kind::Symbol = :missing_unit_pair_inventory,
    supported_terms = _WHITE_LINDSEY_ACCEPTANCE_ONE_BODY_TERMS,
    metadata = (;),
)
    isnothing(source) ||
        return _white_lindsey_decomposed_unit_pair_inventory_result(
            :blocked_white_lindsey_decomposed_unit_pair_inventory,
            :unsupported_decomposed_wl_unit_pair_inventory_source,
            ();
            source_kind,
            supported_terms = Tuple(supported_terms),
            metadata,
        )
    return _white_lindsey_decomposed_unit_pair_inventory_result(
        :blocked_white_lindsey_decomposed_unit_pair_inventory,
        :missing_decomposed_wl_unit_pair_inventory_source,
        ();
        source_kind,
        supported_terms = Tuple(supported_terms),
        metadata,
    )
end

function _white_lindsey_seed_source_available(source::NamedTuple)
    hasproperty(source, :object_kind) &&
        getproperty(source, :object_kind) ===
        :white_lindsey_low_order_materialized_seed_report &&
        hasproperty(source, :fixture) &&
        hasproperty(source.fixture, :sequence) &&
        hasproperty(source, :inventory)
end

function _white_lindsey_seed_decomposed_unit_pairs(
    source::NamedTuple;
    metadata = (;),
)
    units = _white_lindsey_seed_decomposed_units(source; metadata)
    pair_records = CUP.UnitPairRecord[]
    pair_index = 0
    for left_index in eachindex(units)
        for right_index in left_index:length(units)
            pair_index += 1
            left = units[left_index]
            right = units[right_index]
            push!(
                pair_records,
                CUP.UnitPairRecord(
                    (left.unit_key, right.unit_key),
                    pair_index,
                    Symbol(String(left.unit_kind), "__", String(right.unit_kind)),
                    left,
                    right,
                    left_index,
                    right_index,
                    left.unit_key,
                    right.unit_key,
                    left.unit_kind,
                    right.unit_kind,
                    nothing,
                    false,
                    (;
                        source =
                            :white_lindsey_low_order_materialized_seed_ranges,
                        fixed_block_operator_matrices_used = false,
                    ),
                ),
            )
        end
    end
    return Tuple(pair_records)
end

function _white_lindsey_shellification_unit_pair_records_with_ranges(
    plan::CUP.UnitPairPlan;
    metadata::NamedTuple,
)
    materialized =
        _white_lindsey_shellification_retained_units_with_ranges(
            plan.retained_unit_plan;
            metadata,
        )
    materialized.status === :available_shellification_decomposed_wl_retained_units ||
        return (;
            status = :blocked_shellification_decomposed_wl_unit_pair_records,
            blocker = materialized.blocker,
            blocked_unit_key = materialized.blocked_unit_key,
            retained_units = (),
            unit_pairs = (),
        )
    materialized_units = materialized.retained_units
    by_key = Dict(unit.unit_key => unit for unit in materialized_units)
    materialized_pairs = CUP.UnitPairRecord[]
    for pair in CUP.unit_pairs(plan)
        push!(
            materialized_pairs,
            CUP.UnitPairRecord(
                pair.pair_key,
                pair.pair_index,
                pair.pair_family,
                by_key[pair.left_unit_key],
                by_key[pair.right_unit_key],
                pair.left_index,
                pair.right_index,
                pair.left_unit_key,
                pair.right_unit_key,
                pair.left_unit_kind,
                pair.right_unit_kind,
                pair.route_core_pair_sidecar,
                false,
                merge(
                    pair.metadata,
                    (;
                        retained_range_source =
                            :shellification_source_cpb_support_order,
                    ),
                ),
            ),
        )
    end
    return (;
        status = :available_shellification_decomposed_wl_unit_pair_records,
        blocker = nothing,
        blocked_unit_key = nothing,
        retained_units = Tuple(materialized_units),
        unit_pairs = materialized_pairs,
    )
end

function _white_lindsey_shellification_retained_units_with_ranges(
    retained_unit_plan::CRU.RetainedUnitPlan;
    metadata::NamedTuple,
)
    materialized_units = CRU.RetainedUnitRecord[]
    next_column = 1
    for unit in CRU.units(retained_unit_plan)
        materialized =
            _white_lindsey_shellification_retained_unit_with_range(
                unit,
                next_column;
                metadata,
            )
        materialized.status ===
        :available_shellification_decomposed_wl_retained_unit ||
            return (;
                status = :blocked_shellification_decomposed_wl_retained_units,
                blocker = materialized.blocker,
                blocked_unit_key = unit.unit_key,
                retained_units = (),
            )
        push!(materialized_units, materialized.unit)
        next_column = last(materialized.unit.column_range) + 1
    end
    return (;
        status = :available_shellification_decomposed_wl_retained_units,
        blocker = nothing,
        blocked_unit_key = nothing,
        retained_units = Tuple(materialized_units),
    )
end

function _white_lindsey_shellification_retained_unit_with_range(
    unit::CRU.RetainedUnitRecord,
    next_column::Int;
    metadata::NamedTuple,
)
    source_cpb_count = length(unit.source_cpbs)
    source_cpb_count == 1 || return (;
        status = :blocked_shellification_decomposed_wl_retained_unit,
        blocker = :white_lindsey_unit_source_cpb_count_not_one,
        unit = nothing,
    )
    source_cpb = only(unit.source_cpbs)
    stratum_kind = _white_lindsey_shellification_unit_stratum_kind(unit)
    _white_lindsey_shellification_unit_kind_supported(unit, stratum_kind) ||
        return (;
            status = :blocked_shellification_decomposed_wl_retained_unit,
            blocker = :unsupported_decomposed_wl_unit_kind,
            unit = nothing,
        )
    shape = CPB.shape(source_cpb)
    retained_counts =
        _white_lindsey_shellification_unit_retained_counts(
            stratum_kind,
            shape,
            metadata,
        )
    isnothing(retained_counts) && return (;
        status = :blocked_shellification_decomposed_wl_retained_unit,
        blocker = :missing_shellification_wl_retained_count_policy,
        unit = nothing,
    )
    dimension = retained_counts.x * retained_counts.y * retained_counts.z
    column_range = next_column:(next_column + dimension - 1)
    unit_context = _white_lindsey_seed_unit_context_metadata(metadata)
    materialized_metadata = merge(
        unit.metadata,
        unit_context,
        (;
            stratum_kind,
            retained_counts,
            retained_range_source = :shellification_source_cpb_support_order,
            shellification_backed_decomposed_wl_inventory = true,
            low_order_materialized_seed_inventory_used = false,
        ),
    )
    return (;
        status = :available_shellification_decomposed_wl_retained_unit,
        blocker = nothing,
        unit = CRU.RetainedUnitRecord(
            unit.unit_key,
            unit.unit_index,
            unit.unit_kind,
            unit.source_contract_key,
            unit.terminal_region_key,
            unit.terminal_region_role,
            unit.terminal_region_kind,
            unit.lowering_kind,
            unit.retained_rule,
            unit.realization_rule,
            unit.owned_support,
            unit.source_cpbs,
            unit.source_cpb_index,
            :available,
            dimension,
            :available,
            column_range,
            unit.route_core_final_unit,
            false,
            materialized_metadata,
        ),
    )
end

function _white_lindsey_shellification_unit_retained_counts(
    stratum_kind,
    shape,
    metadata::NamedTuple,
)
    shape_counts = (; x = Int(shape[1]), y = Int(shape[2]), z = Int(shape[3]))
    stratum_kind === :direct_core && return shape_counts
    stratum_kind === :corner_cpb && return (; x = 1, y = 1, z = 1)
    nside = _white_lindsey_shellification_nside(metadata)
    isnothing(nside) && return nothing
    nside >= 3 || return nothing
    retained_side = nside - 2
    axes = (:x, :y, :z)
    retained = Dict(axis => 1 for axis in axes)
    for (axis, count) in zip(axes, shape)
        if count > 1
            Int(count) >= nside || return nothing
            retained[axis] = retained_side
        end
    end
    return (; x = retained[:x], y = retained[:y], z = retained[:z])
end

function _white_lindsey_shellification_nside(metadata::NamedTuple)
    haskey(metadata, :nside) && return Int(metadata.nside)
    haskey(metadata, :q) && return Int(metadata.q)
    haskey(metadata, :n_s) && return Int(metadata.n_s)
    haskey(metadata, :ns) && return Int(metadata.ns)
    return nothing
end

function _white_lindsey_shellification_unit_stratum_kind(
    unit::CRU.RetainedUnitRecord,
)
    metadata_kind = _white_lindsey_unit_metadata_value(unit, :stratum_kind)
    !isnothing(metadata_kind) && return metadata_kind
    unit.unit_kind === :direct_cpb_retained_unit && return :direct_core
    return nothing
end

function _white_lindsey_shellification_unit_kind_supported(unit, stratum_kind)
    unit.unit_kind === :white_lindsey_boundary_stratum_retained_unit &&
        stratum_kind in (:facet_cpb, :face_cpb, :edge_cpb, :corner_cpb) &&
        return true
    unit.unit_kind === :direct_cpb_retained_unit && stratum_kind === :direct_core &&
        return true
    return false
end

function _white_lindsey_seed_decomposed_units(source::NamedTuple; metadata = (;))
    shell = only(source.fixture.sequence.shell_layers)
    inventory = source.inventory
    unit_context = _white_lindsey_seed_unit_context_metadata(metadata)
    core_unit = _white_lindsey_seed_core_unit(inventory, unit_context)
    face_units = _white_lindsey_seed_units_for_strata(
        :face,
        shell.faces,
        _white_lindsey_seed_source_cpbs(:face, shell.faces),
        inventory.retained_ranges.faces,
        0,
        unit_context,
    )
    edge_units = _white_lindsey_seed_units_for_strata(
        :edge,
        shell.edges,
        _white_lindsey_seed_source_cpbs(:edge, shell.edges),
        inventory.retained_ranges.edges,
        length(face_units),
        unit_context,
    )
    corner_units = _white_lindsey_seed_units_for_strata(
        :corner,
        shell.corners,
        _white_lindsey_seed_source_cpbs(:corner, shell.corners),
        inventory.retained_ranges.corners,
        length(face_units) + length(edge_units),
        unit_context,
    )
    return Tuple(
        vcat(
            [core_unit],
            collect(face_units),
            collect(edge_units),
            collect(corner_units),
        ),
    )
end

function _white_lindsey_seed_core_unit(inventory, unit_context)
    source_side_count = Int(inventory.source_side_count)
    intervals = ntuple(_ -> 2:(source_side_count - 1), 3)
    retained_counts = ntuple(axis -> length(intervals[axis]), 3)
    source_cpb = CPB.cpb(
        intervals[1],
        intervals[2],
        intervals[3];
        role = :white_lindsey_seed_direct_core_source_cpb,
        metadata = (;
            stratum_kind = :direct_core,
            source_cpb_index = 1,
            retained_counts,
        ),
    )
    column_range = inventory.retained_ranges.core
    dimension = length(column_range)
    metadata = merge(
        (;
            stratum_kind = :direct_core,
            source_cpb_index = 1,
            shell_piece_signature = :direct_core,
            retained_counts = CPB.shape(source_cpb),
            direct_core = true,
        ),
        unit_context,
    )
    return CRU.RetainedUnitRecord(
        :white_lindsey_seed_direct_core,
        1,
        :white_lindsey_boundary_stratum_retained_unit,
        :white_lindsey_seed_direct_core_contract,
        :white_lindsey_seed_direct_core_terminal_region,
        :white_lindsey_seed_complete_shell,
        :synthetic_terminal_region,
        :white_lindsey_boundary_strata,
        :white_lindsey_boundary_stratum_product,
        :direct_or_trivial_embedding,
        CRC.owned_cpb(source_cpb),
        (source_cpb,),
        1,
        :available,
        dimension,
        :available,
        column_range,
        nothing,
        false,
        metadata,
    )
end

function _white_lindsey_seed_source_cpbs(group::Symbol, shell_pieces)
    return Tuple(
        _white_lindsey_seed_source_cpb(group, index, shell_pieces[index])
        for index in eachindex(shell_pieces)
    )
end

function _white_lindsey_seed_source_cpb(
    group::Symbol,
    index::Int,
    shell_piece,
)
    intervals = _white_lindsey_seed_source_cpb_intervals(shell_piece)
    return CPB.cpb(
        intervals.x,
        intervals.y,
        intervals.z;
        role = Symbol(:white_lindsey_seed_, group, :_source_cpb_, index),
        metadata = _white_lindsey_seed_source_cpb_metadata(
            group,
            index,
            shell_piece,
        ),
    )
end

function _white_lindsey_seed_source_cpb_metadata(
    group::Symbol,
    index::Int,
    shell_piece,
)
    base = (;
        stratum_kind = _white_lindsey_seed_stratum_kind(group),
        source_cpb_index = index,
    )
    if hasproperty(shell_piece, :fixed_axis) &&
       hasproperty(shell_piece, :fixed_side)
        return merge(
            base,
            (;
                axis = shell_piece.fixed_axis,
                side = shell_piece.fixed_side,
            ),
        )
    elseif hasproperty(shell_piece, :fixed_axes) &&
           hasproperty(shell_piece, :fixed_sides)
        return merge(
            base,
            (;
                fixed_axes = shell_piece.fixed_axes,
                sides = shell_piece.fixed_sides,
            ),
        )
    end
    return base
end

function _white_lindsey_seed_source_cpb_intervals(shell_piece)
    axes = (:x, :y, :z)
    intervals = Dict(axis => 1:1 for axis in axes)
    if hasproperty(shell_piece, :fixed_axis) &&
       hasproperty(shell_piece, :side_first) &&
       hasproperty(shell_piece, :side_second)
        active_axes = Tuple(axis for axis in axes if axis !== shell_piece.fixed_axis)
        intervals[shell_piece.fixed_axis] =
            shell_piece.fixed_index:shell_piece.fixed_index
        intervals[active_axes[1]] = shell_piece.side_first.interval
        intervals[active_axes[2]] = shell_piece.side_second.interval
    elseif hasproperty(shell_piece, :free_axis) &&
           hasproperty(shell_piece, :fixed_axes) &&
           hasproperty(shell_piece, :fixed_indices)
        intervals[shell_piece.free_axis] = shell_piece.side.interval
        for (axis, coordinate) in
            zip(shell_piece.fixed_axes, shell_piece.fixed_indices)
            intervals[axis] = coordinate:coordinate
        end
    elseif hasproperty(shell_piece, :fixed_indices)
        for (axis, coordinate) in zip(axes, shell_piece.fixed_indices)
            intervals[axis] = coordinate:coordinate
        end
    end
    return (; x = intervals[:x], y = intervals[:y], z = intervals[:z])
end

function _white_lindsey_seed_units_for_strata(
    group::Symbol,
    shell_pieces,
    source_cpbs,
    column_ranges,
    offset::Int,
    unit_context::NamedTuple,
)
    length(source_cpbs) == length(column_ranges) ||
        throw(
            ArgumentError(
                "White-Lindsey seed source CPB count does not match retained ranges",
            ),
        )
    length(shell_pieces) == length(column_ranges) ||
        throw(
            ArgumentError(
                "White-Lindsey seed shell piece count does not match retained ranges",
            ),
        )
    return Tuple(
        _white_lindsey_seed_retained_unit(
            group,
            index + offset,
            source_cpbs[index],
            column_ranges[index],
            shell_pieces[index],
            unit_context,
        ) for index in eachindex(source_cpbs)
    )
end

function _white_lindsey_seed_retained_unit(
    group::Symbol,
    unit_index::Int,
    source_cpb,
    column_range,
    shell_piece,
    unit_context::NamedTuple,
)
    unit_key = Symbol(:white_lindsey_seed_, group, :_, unit_index)
    dimension = length(column_range)
    unit_metadata = merge(
        (;
            stratum_kind = _white_lindsey_seed_stratum_kind(group),
            source_cpb_index = unit_index,
            shell_piece_signature =
                _white_lindsey_seed_shell_piece_signature(shell_piece),
            retained_counts =
                _white_lindsey_seed_unit_retained_counts(shell_piece, source_cpb),
            retained_range_source =
                :white_lindsey_low_order_materialized_seed_inventory,
            fixed_block_operator_matrices_used = false,
        ),
        unit_context,
    )
    return CRU.RetainedUnitRecord(
        unit_key,
        unit_index,
        :white_lindsey_boundary_stratum_retained_unit,
        Symbol(unit_key, :_contract),
        Symbol(unit_key, :_terminal_region),
        :white_lindsey_seed_complete_shell,
        :complete_shell,
        :white_lindsey_boundary_strata,
        :white_lindsey_boundary_stratum_product,
        :direct_or_trivial_embedding,
        CRC.owned_cpb(source_cpb),
        (source_cpb,),
        unit_index,
        :available,
        dimension,
        :available,
        column_range,
        nothing,
        false,
        unit_metadata,
    )
end

function _white_lindsey_seed_unit_retained_counts(shell_piece, source_cpb)
    axes = (:x, :y, :z)
    counts = Dict(axis => 1 for axis in axes)
    if hasproperty(shell_piece, :fixed_axis) &&
       hasproperty(shell_piece, :side_first) &&
       hasproperty(shell_piece, :side_second)
        active_axes = Tuple(axis for axis in axes if axis !== shell_piece.fixed_axis)
        counts[active_axes[1]] = Int(shell_piece.side_first.retained_count)
        counts[active_axes[2]] = Int(shell_piece.side_second.retained_count)
    elseif hasproperty(shell_piece, :free_axis) && hasproperty(shell_piece, :side)
        counts[shell_piece.free_axis] = Int(shell_piece.side.retained_count)
    else
        source_shape = CPB.shape(source_cpb)
        counts[:x] = Int(source_shape[1])
        counts[:y] = Int(source_shape[2])
        counts[:z] = Int(source_shape[3])
    end
    return (; x = counts[:x], y = counts[:y], z = counts[:z])
end

function _white_lindsey_seed_unit_context_metadata(metadata::NamedTuple)
    context = NamedTuple()
    parent_dims =
        haskey(metadata, :parent_dims) ?
        metadata.parent_dims :
        (
            haskey(metadata, :parent_axis_counts) ?
            metadata.parent_axis_counts :
            nothing
        )
    doside_source =
        haskey(metadata, :doside_source_1d) ?
        metadata.doside_source_1d :
        (
            haskey(metadata, :parent_axis_bundle_object) ?
            metadata.parent_axis_bundle_object :
            nothing
        )
    isnothing(parent_dims) ||
        (context = merge(context, (; parent_dims = parent_dims)))
    isnothing(doside_source) ||
        (context = merge(context, (; doside_source_1d = doside_source)))
    return context
end

function _white_lindsey_seed_stratum_kind(group::Symbol)
    group === :face && return :facet_cpb
    group === :edge && return :edge_cpb
    group === :corner && return :corner_cpb
    return :unknown_cpb
end

function _white_lindsey_seed_shell_piece_signature(piece)
    fields = propertynames(piece)
    pairs = NamedTuple()
    for key in (
        :face_kind,
        :fixed_axis,
        :fixed_side,
        :fixed_index,
        :free_axis,
        :fixed_axes,
        :fixed_sides,
        :fixed_indices,
    )
        if key in fields
            pairs = merge(pairs, NamedTuple{(key,)}((getproperty(piece, key),)))
        end
    end
    return pairs
end

function _white_lindsey_decomposed_inventory_units(pairs)
    by_key = Dict{Symbol,Any}()
    ordered = Any[]
    for pair in pairs
        for unit in (pair.left_unit, pair.right_unit)
            if !haskey(by_key, unit.unit_key)
                by_key[unit.unit_key] = unit
                push!(ordered, unit)
            end
        end
    end
    return Tuple(ordered)
end

function _white_lindsey_decomposed_inventory_blocker(units)
    all(
        unit -> unit.unit_kind in (
            :direct_cpb_retained_unit,
            :white_lindsey_boundary_stratum_retained_unit,
        ),
        units,
    ) || return :unsupported_decomposed_wl_unit_kind
    all(unit -> unit.dimension_status === :available && !isnothing(unit.dimension), units) ||
        return :missing_retained_unit_dimension
    all(unit -> unit.column_range_status === :available && !isnothing(unit.column_range), units) ||
        return :missing_retained_unit_column_ranges
    all(
        unit -> length(unit.column_range) == unit.dimension,
        units,
    ) || return :retained_unit_range_dimension_mismatch
    return nothing
end

function _white_lindsey_decomposed_inventory_global_dimension(units, blocker)
    !isnothing(blocker) && return nothing
    isempty(units) && return nothing
    return maximum(last(unit.column_range) for unit in units)
end

function _white_lindsey_decomposed_inventory_unit_summaries(units)
    return Tuple(
        (;
            unit_key = unit.unit_key,
            unit_index = unit.unit_index,
            unit_kind = unit.unit_kind,
            stratum_kind =
                _white_lindsey_unit_metadata_value(unit, :stratum_kind),
            source_cpb_index = unit.source_cpb_index,
            source_cpb_count = length(unit.source_cpbs),
            dimension_status = unit.dimension_status,
            dimension = unit.dimension,
            column_range_status = unit.column_range_status,
            column_range = unit.column_range,
        ) for unit in units
    )
end

function _white_lindsey_decomposed_inventory_pair_summaries(pairs, blocker)
    length(pairs) > _WHITE_LINDSEY_COMPACT_PAIR_SUMMARY_LIMIT &&
        return :omitted_large_decomposed_wl_pair_inventory
    return Tuple(
        (;
            pair_key = pair.pair_key,
            pair_index = pair.pair_index,
            pair_family = pair.pair_family,
            left_unit_key = pair.left_unit_key,
            right_unit_key = pair.right_unit_key,
            left_column_range =
                isnothing(blocker) ? pair.left_unit.column_range : nothing,
            right_column_range =
                isnothing(blocker) ? pair.right_unit.column_range : nothing,
            left_dimension =
                isnothing(blocker) ? pair.left_unit.dimension : nothing,
            right_dimension =
                isnothing(blocker) ? pair.right_unit.dimension : nothing,
        ) for pair in pairs
    )
end

function _white_lindsey_decomposed_inventory_pair_keys(pairs)
    length(pairs) > _WHITE_LINDSEY_COMPACT_PAIR_SUMMARY_LIMIT &&
        return :omitted_large_decomposed_wl_pair_inventory
    return Tuple(pair.pair_key for pair in pairs)
end

function _white_lindsey_decomposed_unit_pair_inventory_result(
    status::Symbol,
    blocker,
    pairs;
    source_kind::Symbol,
    supported_terms,
    metadata,
    units = (),
)
    available = status === :available_white_lindsey_decomposed_unit_pair_inventory
    unit_summaries = _white_lindsey_decomposed_inventory_unit_summaries(units)
    pair_summaries =
        _white_lindsey_decomposed_inventory_pair_summaries(pairs, blocker)
    retained_dimension =
        _white_lindsey_decomposed_inventory_global_dimension(units, blocker)

    return (;
        object_kind = :white_lindsey_decomposed_unit_pair_inventory,
        status,
        blocker,
        source_kind,
        supported_terms,
        term_compatibility = (;
            overlap = :overlap in supported_terms,
            kinetic = :kinetic in supported_terms,
            electron_nuclear_by_center =
                :electron_nuclear_by_center in supported_terms,
        ),
        unit_count = length(units),
        pair_count = length(pairs),
        retained_units = units,
        unit_pairs = pairs,
        unit_keys = Tuple(unit.unit_key for unit in units),
        pair_keys = _white_lindsey_decomposed_inventory_pair_keys(pairs),
        unit_summaries,
        pair_summaries,
        retained_dimension,
        retained_dimension_status =
            available ?
            :available_from_decomposed_wl_unit_column_ranges :
            :not_available,
        retained_unit_column_ranges_materialized = available,
        decomposed_unit_pair_column_ranges_available = available,
        decomposed_unit_pair_records_available = available,
        decomposed_wl_unit_pair_inventory_available = available,
        route_owned_inventory = true,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        local_operator_blocks_materialized = false,
        global_matrices_materialized = false,
        hamiltonian_data_materialized = false,
        route_driver_wiring = false,
        metadata = NamedTuple(metadata),
    )
end
