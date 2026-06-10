# Route-owned decomposed White-Lindsey unit-pair inventory.
#
# This validates retained-unit pair metadata that later acceptance assembly can
# consume. It does not build local blocks, assign ranges, assemble matrices, or
# call the legacy direct Cartesian operator path.

const _WHITE_LINDSEY_ACCEPTANCE_ONE_BODY_TERMS = (
    :overlap,
    :kinetic,
    :electron_nuclear_by_center,
)

function white_lindsey_decomposed_unit_pair_inventory(
    plan::CUP.UnitPairPlan;
    supported_terms = _WHITE_LINDSEY_ACCEPTANCE_ONE_BODY_TERMS,
    metadata = (;),
)
    return white_lindsey_decomposed_unit_pair_inventory(
        CUP.unit_pairs(plan);
        source_kind = :cartesian_unit_pair_plan,
        supported_terms,
        metadata,
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
    return white_lindsey_decomposed_unit_pair_inventory(Tuple(pairs); kwargs...)
end

function white_lindsey_decomposed_unit_pair_inventory(
    pairs::Tuple{Vararg{CUP.UnitPairRecord}};
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
        unit -> unit.unit_kind === :white_lindsey_boundary_stratum_retained_unit,
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
        pair_keys = Tuple(pair.pair_key for pair in pairs),
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
