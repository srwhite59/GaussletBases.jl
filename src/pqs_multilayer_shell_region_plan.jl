# Shellification/lowering-backed region planning for multi-layer PQS shells.

function _pqs_multilayer_axis_metrics(bundles::_CartesianNestedAxisBundles3D)
    pgdg_x = _nested_axis_pgdg(bundles, :x)
    pgdg_y = _nested_axis_pgdg(bundles, :y)
    pgdg_z = _nested_axis_pgdg(bundles, :z)
    return (;
        x = (overlap = pgdg_x.overlap, kinetic = pgdg_x.kinetic),
        y = (overlap = pgdg_y.overlap, kinetic = pgdg_y.kinetic),
        z = (overlap = pgdg_z.overlap, kinetic = pgdg_z.kinetic),
    )
end

function _pqs_multilayer_box_depth(
    outer_box::NTuple{3,UnitRange{Int}},
    core_box::NTuple{3,UnitRange{Int}},
)
    lower_depths = ntuple(axis -> first(core_box[axis]) - first(outer_box[axis]), 3)
    upper_depths = ntuple(axis -> last(outer_box[axis]) - last(core_box[axis]), 3)
    all(depth -> depth >= 1, lower_depths) ||
        throw(ArgumentError("multi-layer PQS shell plan requires the core box to be strictly inside the outer box"))
    lower_depths == upper_depths ||
        throw(ArgumentError("multi-layer PQS shell plan currently requires symmetric shell depth on every axis"))
    all(depth -> depth == lower_depths[1], lower_depths) ||
        throw(ArgumentError("multi-layer PQS shell plan currently requires the same shell depth on every axis"))
    return lower_depths[1]
end

function _pqs_multilayer_core_box_at_depth(
    core_box::NTuple{3,UnitRange{Int}},
    depth::Int,
)
    return ntuple(axis -> (first(core_box[axis]) - depth):(last(core_box[axis]) + depth), 3)
end

function _pqs_multilayer_block_concatenate_shell_coefficients(shell_records)
    total_support = sum(record -> record.support_count, shell_records)
    total_retained = sum(record -> record.retained_count, shell_records)
    coefficients = zeros(Float64, total_support, total_retained)
    support_offset = 0
    retained_offset = 0
    for record in shell_records
        support_range = (support_offset + 1):(support_offset + record.support_count)
        retained_range = (retained_offset + 1):(retained_offset + record.retained_count)
        coefficients[support_range, retained_range] .= record.shell_final_coefficients
        support_offset += record.support_count
        retained_offset += record.retained_count
    end
    return coefficients
end

function _pqs_multilayer_duplicate_count(values)
    return length(values) - length(unique(values))
end

struct PQSMultilayerShellLayerRegion
    layer_index::Int
    terminal_region::Any
    lowering_contract::Any
    current_box::NTuple{3,UnitRange{Int}}
    inner_box::NTuple{3,UnitRange{Int}}
    source_cpbs::Vector{CartesianCPB.CoordinateProductBox}
    support_count::Int
    source_mode_shape::NTuple{3,Int}
    metadata::NamedTuple
end

struct PQSMultilayerShellRegionPlan
    status::Symbol
    blocker::Union{Nothing,Symbol}
    shellification_plan::Any
    lowering_plan::Any
    core_region::Any
    core_contract::Any
    core_box::NTuple{3,UnitRange{Int}}
    outer_box::NTuple{3,UnitRange{Int}}
    shell_layers::Vector{PQSMultilayerShellLayerRegion}
    coverage::NamedTuple
    summary::NamedTuple
    metadata::NamedTuple
end

function _pqs_multilayer_property(plan, key::Symbol, default = nothing)
    return hasproperty(plan, key) ? getproperty(plan, key) : default
end

function _pqs_multilayer_region_support_indices(region, dims)
    raw = region.raw_region
    outer = _nested_box_support_indices(raw.outer_box[1], raw.outer_box[2], raw.outer_box[3], dims)
    isnothing(raw.inner_exclusion_box) && return outer
    inner = Set(
        _nested_box_support_indices(
            raw.inner_exclusion_box[1],
            raw.inner_exclusion_box[2],
            raw.inner_exclusion_box[3],
            dims,
        ),
    )
    return [index for index in outer if !(index in inner)]
end

function _pqs_multilayer_contract_source_mode_shape(contract, raw_region)
    shape = get(contract.metadata, :source_mode_shape, nothing)
    isnothing(shape) && return Tuple(length.(raw_region.outer_box))
    return Tuple(Int(value) for value in shape)
end

"""
    pqs_multilayer_shell_region_plan(shellification_plan, lowering_plan; ...)

Build the geometry/lowering-owned region plan consumed by multi-layer PQS
source realization. The plan records direct-core and ordered complete-shell
regions, selected PQS lowering contracts, source CPBs, and compact coverage
fingerprints. It does not build PQS descriptors, Lowdin matrices, operators,
final-basis data, H1, IDA, density-density, RHF, driver data, exports, or
artifacts.
"""
function pqs_multilayer_shell_region_plan(
    shellification_plan::CartesianShellification.ShellificationPlan,
    lowering_plan::CartesianTerminalLowering.TerminalLoweringPlan;
    metadata = (;),
)
    regions = collect(CartesianShellification.terminal_regions(shellification_plan))
    contracts = collect(CartesianTerminalLowering.selected_contracts(lowering_plan))
    length(regions) == length(contracts) ||
        throw(ArgumentError("PQS shell region plan requires one selected lowering contract per terminal region"))

    contract_by_region = Dict(contract.terminal_region_key => contract for contract in contracts)
    core_regions = [region for region in regions if region.region_kind === :direct_core]
    shell_regions = sort!(
        [region for region in regions if region.region_kind === :complete_shell];
        by = region -> region.order_index,
    )
    length(core_regions) == 1 ||
        throw(ArgumentError("PQS multi-layer shell region plan currently requires exactly one direct core"))
    !isempty(shell_regions) ||
        throw(ArgumentError("PQS multi-layer shell region plan requires at least one complete shell layer"))

    core_region = only(core_regions)
    core_contract = contract_by_region[core_region.key]
    core_box = core_region.raw_region.outer_box
    dims = CartesianShellification.raw_plan(shellification_plan).parent_dims
    core_support_indices = _pqs_multilayer_region_support_indices(core_region, dims)

    shell_layers = PQSMultilayerShellLayerRegion[]
    shell_support_indices = Int[]
    for (layer_index, region) in enumerate(shell_regions)
        contract = contract_by_region[region.key]
        contract.lowering_kind === :pqs_filled_source_cpb ||
            throw(ArgumentError("complete shell region $(region.key) is not selected for PQS lowering"))
        raw = region.raw_region
        isnothing(raw.inner_exclusion_box) &&
            throw(ArgumentError("PQS complete shell region $(region.key) requires an inner exclusion box"))
        source_cpbs = collect(CartesianTerminalLowering.source_cpbs(contract))
        support_indices = _pqs_multilayer_region_support_indices(region, dims)
        append!(shell_support_indices, support_indices)
        push!(
            shell_layers,
            PQSMultilayerShellLayerRegion(
                layer_index,
                region,
                contract,
                raw.outer_box,
                raw.inner_exclusion_box,
                source_cpbs,
                length(support_indices),
                _pqs_multilayer_contract_source_mode_shape(contract, raw),
                (;
                    terminal_region_key = region.key,
                    terminal_region_role = region.role,
                    terminal_region_kind = region.region_kind,
                    lowering_kind = contract.lowering_kind,
                    source_cpb_count = length(source_cpbs),
                ),
            ),
        )
    end

    outer_box = last(shell_layers).current_box
    combined_support_indices = vcat(core_support_indices, shell_support_indices)
    intended_support_indices =
        _nested_box_support_indices(outer_box[1], outer_box[2], outer_box[3], dims)
    combined_unique = sort!(unique(combined_support_indices))
    intended_sorted = sort!(collect(intended_support_indices))
    shell_duplicate_count = _pqs_multilayer_duplicate_count(shell_support_indices)
    core_shell_duplicate_count = length(intersect(core_support_indices, shell_support_indices))
    covers_outer_box = combined_unique == intended_sorted
    blocker =
        shell_duplicate_count != 0 ? :duplicate_shell_support :
        core_shell_duplicate_count != 0 ? :core_shell_support_overlap :
        !covers_outer_box ? :combined_support_does_not_cover_outer_box :
        nothing
    status = isnothing(blocker) ?
        :available_pqs_multilayer_shell_region_plan :
        :blocked_pqs_multilayer_shell_region_plan
    coverage = (;
        shellification_coverage = CartesianShellification.coverage(shellification_plan),
        core_support_count = length(core_support_indices),
        shell_support_count = length(shell_support_indices),
        combined_support_count = length(combined_support_indices),
        intended_support_count = length(intended_sorted),
        shell_duplicate_count,
        core_shell_duplicate_count,
        covers_outer_box,
    )
    summary = (;
        status,
        blocker,
        source = :cartesian_shellification_and_terminal_lowering,
        core_region_key = core_region.key,
        core_support_count = length(core_support_indices),
        shell_layer_count = length(shell_layers),
        shell_support_count = length(shell_support_indices),
        outer_support_coverage = covers_outer_box,
        shell_duplicate_count,
        core_shell_duplicate_count,
        pqs_descriptors_materialized = false,
        lowdin_materialized = false,
        support_operator_blocks_materialized = false,
        h1_materialized = false,
        ida_materialized = false,
        rhf_materialized = false,
    )

    return PQSMultilayerShellRegionPlan(
        status,
        blocker,
        shellification_plan,
        lowering_plan,
        core_region,
        core_contract,
        core_box,
        outer_box,
        shell_layers,
        coverage,
        summary,
        NamedTuple(metadata),
    )
end
