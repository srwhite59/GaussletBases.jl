"""
    _BondAlignedDiatomicSplitGeometry3D

Geometry report for the first bond-aligned diatomic split/no-split decision.

The split is attached to the original parent grid:

- `working_box` is the shared box remaining after the outer shared shell stage
- `split_index` is the bond-axis parent-grid index nearest the bond midpoint
- `shared_midpoint_box` is the direct shared midpoint slab for the odd-length
  homonuclear case
- `child_boxes` are the two nonoverlapping child boxes if the split is allowed
- `child_physical_widths` records the mapped physical widths of those children
"""
struct _BondAlignedDiatomicSplitGeometry3D
    parent_box::NTuple{3,UnitRange{Int}}
    working_box::NTuple{3,UnitRange{Int}}
    bond_axis::Symbol
    midpoint::Float64
    split_index::Int
    count_eligible::Bool
    unsplit_aspect_eligible::Bool
    shape_eligible::Bool
    did_split::Bool
    shared_midpoint_box::Union{Nothing,NTuple{3,UnitRange{Int}}}
    child_boxes::Vector{NTuple{3,UnitRange{Int}}}
    child_physical_widths::Vector{NTuple{3,Float64}}
end

"""
    _CartesianNestedBondAlignedDiatomicSource3D

First bond-aligned diatomic nested fixed-space source built on top of the
existing atomic shell language.

The source keeps:

- the mixed-axis parent bundle data
- the shared-box split/no-split geometry decision
- the outer shared shell layers
- the child atomic-style subtrees after the bond-axis split
- the merged shell-sequence object used to build the fixed block
"""
struct _CartesianNestedBondAlignedDiatomicSource3D{B}
    basis::B
    axis_bundles::_CartesianNestedAxisBundles3D
    nside::Int
    child_shell_retention_contract::CartesianNestedCompleteShellRetentionContract
    shared_shell_retention_contract::CartesianNestedCompleteShellRetentionContract
    geometry::_BondAlignedDiatomicSplitGeometry3D
    shared_shell_layers::Vector{_CartesianNestedCompleteShell3D}
    child_sequences::Vector{_CartesianNestedShellSequence3D}
    child_column_ranges::Vector{UnitRange{Int}}
    midpoint_slab_column_range::Union{Nothing,UnitRange{Int}}
    sequence::_CartesianNestedShellSequence3D
end

function Base.show(io::IO, geometry::_BondAlignedDiatomicSplitGeometry3D)
    print(
        io,
        "_BondAlignedDiatomicSplitGeometry3D(axis=:",
        geometry.bond_axis,
        ", working_box=",
        geometry.working_box,
        ", split_index=",
        geometry.split_index,
        ", midpoint_slab=",
        isnothing(geometry.shared_midpoint_box) ? "nothing" : geometry.shared_midpoint_box,
        ", did_split=",
        geometry.did_split,
        ")",
    )
end

function Base.show(io::IO, source::_CartesianNestedBondAlignedDiatomicSource3D)
    print(
        io,
        "_CartesianNestedBondAlignedDiatomicSource3D(nshared=",
        length(source.shared_shell_layers),
        ", nchild=",
        length(source.child_sequences),
        ", nfixed=",
        size(source.sequence.coefficient_matrix, 2),
        ", nside=",
        source.nside,
        ", midpoint_slab=",
        !isnothing(source.midpoint_slab_column_range),
        ", did_split=",
        source.geometry.did_split,
        ")",
    )
end

function _nested_source_contract_audit(source::_CartesianNestedBondAlignedDiatomicSource3D)
    return _nested_shell_sequence_contract_audit(source.sequence, _nested_axis_lengths(source.axis_bundles))
end

function _nested_diatomic_midpoint_row_index(
    centers_axis::AbstractVector{<:Real},
    interval::UnitRange{Int},
    midpoint::Real,
)
    length(interval) >= 2 || throw(
        ArgumentError("diatomic midpoint splitting requires at least two raw sites on the bond axis"),
    )
    candidates = collect(first(interval):(last(interval) - 1))
    _, local_index = findmin(abs.(Float64.(centers_axis[candidates]) .- Float64(midpoint)))
    return candidates[local_index]
end

function _nested_diatomic_midpoint_slab_split(
    box::NTuple{3,UnitRange{Int}},
    bond_axis::Symbol,
    split_index::Int,
)
    axis = bond_axis == :x ? 1 : bond_axis == :y ? 2 : bond_axis == :z ? 3 : 0
    axis != 0 || throw(ArgumentError("bond-axis midpoint-slab construction requires bond_axis = :x, :y, or :z"))
    interval = box[axis]
    first(interval) < split_index < last(interval) || throw(
        ArgumentError("diatomic midpoint-slab construction requires the split index to lie strictly inside the working box"),
    )
    left_axis = first(interval):(split_index - 1)
    slab_axis = split_index:split_index
    right_axis = (split_index + 1):last(interval)
    left_box =
        axis == 1 ? (left_axis, box[2], box[3]) :
        axis == 2 ? (box[1], left_axis, box[3]) :
        (box[1], box[2], left_axis)
    slab_box =
        axis == 1 ? (slab_axis, box[2], box[3]) :
        axis == 2 ? (box[1], slab_axis, box[3]) :
        (box[1], box[2], slab_axis)
    right_box =
        axis == 1 ? (right_axis, box[2], box[3]) :
        axis == 2 ? (box[1], right_axis, box[3]) :
        (box[1], box[2], right_axis)
    return left_box, slab_box, right_box
end

function _nested_direct_box_coefficients(
    dims::NTuple{3,Int},
    support_indices::AbstractVector{Int},
)
    return _nested_sparse_coefficient_map(
        Int.(support_indices),
        collect(1:length(support_indices)),
        ones(Float64, length(support_indices)),
        prod(dims),
        length(support_indices),
    )
end

function _nested_direct_box_coefficients(
    bundles::_CartesianNestedAxisBundles3D,
    box::NTuple{3,UnitRange{Int}},
)
    dims = _nested_axis_lengths(bundles)
    support_indices = _nested_box_support_indices(box..., dims)
    return (
        support_indices = support_indices,
        coefficient_matrix = _nested_direct_box_coefficients(dims, support_indices),
    )
end

function _nested_diatomic_children_are_roughly_cubic(
    bundles::_CartesianNestedAxisBundles3D,
    child_boxes::AbstractVector{<:NTuple{3,UnitRange{Int}}},
    bond_axis::Symbol;
    min_parallel_to_transverse_ratio::Float64 = 0.4,
)
    min_parallel_to_transverse_ratio > 0.0 || throw(
        ArgumentError("diatomic anti-sliver check requires min_parallel_to_transverse_ratio > 0"),
    )
    axis = bond_axis == :x ? 1 : bond_axis == :y ? 2 : bond_axis == :z ? 3 : 0
    axis != 0 || throw(ArgumentError("bond-axis anti-sliver check requires bond_axis = :x, :y, or :z"))
    for child_box in child_boxes
        widths = _nested_box_physical_widths(bundles, child_box)
        parallel = widths[axis]
        transverse = maximum(widths[index] for index in 1:3 if index != axis)
        parallel > 0.0 || return false
        transverse > 0.0 || return false
        parallel >= min_parallel_to_transverse_ratio * transverse || return false
    end
    return true
end

function _nested_diatomic_resolve_core_near_nucleus_protect_rows(
    value::Union{Symbol,Integer},
    nside::Int,
)
    value === :auto && return max(0, div(nside, 2) - 1)
    value isa Integer && return max(0, Int(value))
    throw(
        ArgumentError(
            "core_near_nucleus_protect_rows must be :auto or an integer row count",
        ),
    )
end

function _nested_diatomic_axis_minimum_retain(
    retention::CartesianNestedCompleteShellRetentionContract,
    axis::Symbol,
)
    if axis == :x
        return maximum((retention.retain_xy[1], retention.retain_xz[1], retention.retain_x_edge))
    elseif axis == :y
        return maximum((retention.retain_xy[2], retention.retain_yz[1], retention.retain_y_edge))
    elseif axis == :z
        return maximum((retention.retain_xz[2], retention.retain_yz[2], retention.retain_z_edge))
    end
    throw(ArgumentError("diatomic adaptive retain lookup requires axis = :x, :y, or :z"))
end

function _nested_diatomic_candidate_counts(
    length_value::Int,
    minimum_retained::Int,
)
    length_value >= 1 || throw(
        ArgumentError("diatomic adaptive retain-count generation requires a positive line length"),
    )
    lower = min(length_value, max(1, minimum_retained))
    return collect(lower:length_value)
end

function _nested_diatomic_line_points_for_interval(
    bundles::_CartesianNestedAxisBundles3D,
    axis::Symbol,
    interval::UnitRange{Int};
    retained_count::Int,
    enforce_symmetric_odd::Bool = false,
)
    centers_axis = _nested_axis_pgdg(bundles, axis).centers
    if retained_count >= length(interval)
        return Float64[Float64(centers_axis[index]) for index in interval]
    end
    side = _nested_doside_1d(
        _nested_axis_pgdg(bundles, axis),
        interval,
        retained_count;
        enforce_symmetric_odd = enforce_symmetric_odd,
    )
    return Float64[Float64(value) for value in side.localized_centers]
end

function _nested_diatomic_uniform_line_points(
    low::Real,
    high::Real,
    retained_count::Int,
)
    retained_count >= 1 || throw(
        ArgumentError("diatomic ideal reference line construction requires retained_count >= 1"),
    )
    if retained_count == 1
        return Float64[0.5 * (Float64(low) + Float64(high))]
    end
    return collect(Float64, range(Float64(low), Float64(high); length = retained_count))
end

function _nested_diatomic_line_theta_stats(
    line_points::AbstractVector{<:Real},
    line_axis::Symbol,
    fixed_coords::NTuple{2,<:Real},
    basis,
)
    length(line_points) <= 1 && return (
        theta_min = NaN,
        theta_max = NaN,
    )
    theta_min = Inf
    theta_max = 0.0
    for index in 1:(length(line_points) - 1)
        left = Float64(line_points[index])
        right = Float64(line_points[index + 1])
        midpoint = 0.5 * (left + right)
        x, y, z =
            line_axis == :x ? (midpoint, Float64(fixed_coords[1]), Float64(fixed_coords[2])) :
            line_axis == :y ? (Float64(fixed_coords[1]), midpoint, Float64(fixed_coords[2])) :
            line_axis == :z ? (Float64(fixed_coords[1]), Float64(fixed_coords[2]), midpoint) :
            throw(
                ArgumentError(
                    "diatomic line angular diagnostics require line_axis = :x, :y, or :z",
                ),
            )
        gap = right - left
        rmin = minimum(
            sqrt(
                (x - Float64(nucleus[1]))^2 +
                (y - Float64(nucleus[2]))^2 +
                (z - Float64(nucleus[3]))^2,
            ) for nucleus in basis.nuclei
        )
        theta = rmin > 0.0 ? gap / rmin : Inf
        theta_min = min(theta_min, theta)
        theta_max = max(theta_max, theta)
    end
    return (
        theta_min = theta_min,
        theta_max = theta_max,
    )
end

function _nested_diatomic_aggregate_line_theta_stats(
    line_points::AbstractVector{<:Real},
    line_axis::Symbol,
    fixed_coord_pairs::AbstractVector{<:NTuple{2,<:Real}},
    basis,
)
    isempty(fixed_coord_pairs) && throw(
        ArgumentError(
            "diatomic aggregated line-theta statistics require at least one fixed-coordinate pair",
        ),
    )
    theta_min = Inf
    theta_max = 0.0
    for fixed_coords in fixed_coord_pairs
        stats = _nested_diatomic_line_theta_stats(line_points, line_axis, fixed_coords, basis)
        theta_min = min(theta_min, stats.theta_min)
        theta_max = max(theta_max, stats.theta_max)
    end
    return (
        theta_min = theta_min,
        theta_max = theta_max,
    )
end

function _nested_diatomic_unique_fixed_coord_pairs(
    pairs::AbstractVector{<:NTuple{2,<:Real}};
    atol::Float64 = 1.0e-12,
)
    unique_pairs = NTuple{2,Float64}[]
    for pair in pairs
        candidate = (Float64(pair[1]), Float64(pair[2]))
        any(
            existing ->
                isapprox(candidate[1], existing[1]; atol = atol, rtol = 0.0) &&
                isapprox(candidate[2], existing[2]; atol = atol, rtol = 0.0),
            unique_pairs,
        ) && continue
        push!(unique_pairs, candidate)
    end
    return unique_pairs
end

function _nested_diatomic_boundary_fixed_coord_pairs(
    bundles::_CartesianNestedAxisBundles3D,
    box::NTuple{3,UnitRange{Int}},
    line_axis::Symbol,
)
    centers_x = _nested_axis_pgdg(bundles, :x).centers
    centers_y = _nested_axis_pgdg(bundles, :y).centers
    centers_z = _nested_axis_pgdg(bundles, :z).centers
    pairs =
        line_axis == :x ? [(Float64(centers_y[iy]), Float64(centers_z[iz])) for iy in (first(box[2]), last(box[2])), iz in (first(box[3]), last(box[3]))] :
        line_axis == :y ? [(Float64(centers_x[ix]), Float64(centers_z[iz])) for ix in (first(box[1]), last(box[1])), iz in (first(box[3]), last(box[3]))] :
        line_axis == :z ? [(Float64(centers_x[ix]), Float64(centers_y[iy])) for ix in (first(box[1]), last(box[1])), iy in (first(box[2]), last(box[2]))] :
        throw(
            ArgumentError(
                "diatomic boundary fixed-coordinate selection requires line_axis = :x, :y, or :z",
            ),
        )
    return _nested_diatomic_unique_fixed_coord_pairs(vec(pairs))
end

function _nested_diatomic_reference_bounds(
    bundles::_CartesianNestedAxisBundles3D,
    box::NTuple{3,UnitRange{Int}},
)
    widths = _nested_box_physical_widths(bundles, box)
    side_length = minimum(widths)
    centers_x = _nested_axis_pgdg(bundles, :x).centers
    centers_y = _nested_axis_pgdg(bundles, :y).centers
    centers_z = _nested_axis_pgdg(bundles, :z).centers
    lows = (
        Float64(centers_x[first(box[1])]),
        Float64(centers_y[first(box[2])]),
        Float64(centers_z[first(box[3])]),
    )
    highs = (
        Float64(centers_x[last(box[1])]),
        Float64(centers_y[last(box[2])]),
        Float64(centers_z[last(box[3])]),
    )
    mids = ntuple(index -> 0.5 * (lows[index] + highs[index]), 3)
    return ntuple(index -> (mids[index] - 0.5 * side_length, mids[index] + 0.5 * side_length), 3)
end

function _nested_diatomic_reference_fixed_coord_pairs(
    reference_bounds::NTuple{3,NTuple{2,Float64}},
    line_axis::Symbol,
)
    pairs =
        line_axis == :x ? [(reference_bounds[2][iy], reference_bounds[3][iz]) for iy in 1:2, iz in 1:2] :
        line_axis == :y ? [(reference_bounds[1][ix], reference_bounds[3][iz]) for ix in 1:2, iz in 1:2] :
        line_axis == :z ? [(reference_bounds[1][ix], reference_bounds[2][iy]) for ix in 1:2, iy in 1:2] :
        throw(
            ArgumentError(
                "diatomic ideal-reference fixed-coordinate selection requires line_axis = :x, :y, or :z",
            ),
        )
    return _nested_diatomic_unique_fixed_coord_pairs(vec(pairs))
end

function _nested_diatomic_reference_band(
    basis,
    bundles::_CartesianNestedAxisBundles3D,
    box::NTuple{3,UnitRange{Int}};
    nside::Int,
    reference_fudge_factor::Float64,
)
    return @timeg "diatomic.source.reference_band" begin
        reference_fudge_factor > 0.0 || throw(
            ArgumentError("diatomic ideal-reference tolerance expansion requires reference_fudge_factor > 0"),
        )
        reference_retain = min(nside, minimum(length.(box)))
        reference_bounds = _nested_diatomic_reference_bounds(bundles, box)
        axis_stats = Pair{Symbol,NamedTuple}[]
        for axis in (:x, :y, :z)
            bounds = axis == :x ? reference_bounds[1] : axis == :y ? reference_bounds[2] : reference_bounds[3]
            points = _nested_diatomic_uniform_line_points(bounds[1], bounds[2], reference_retain)
            stats = _nested_diatomic_aggregate_line_theta_stats(
                points,
                axis,
                _nested_diatomic_reference_fixed_coord_pairs(reference_bounds, axis),
                basis,
            )
            push!(axis_stats, axis => stats)
        end
        ideal_theta_min = minimum(last(pair).theta_min for pair in axis_stats)
        ideal_theta_max = maximum(last(pair).theta_max for pair in axis_stats)
        (
            reference_bounds = reference_bounds,
            reference_retain = reference_retain,
            ideal_theta_min = ideal_theta_min,
            ideal_theta_max = ideal_theta_max,
            theta_min = ideal_theta_min / reference_fudge_factor,
            theta_max = ideal_theta_max * reference_fudge_factor,
            reference_fudge_factor = reference_fudge_factor,
        )
    end
end

function _nested_diatomic_shared_shell_reference_band(
    basis,
    bundles::_CartesianNestedAxisBundles3D,
    box::NTuple{3,UnitRange{Int}};
    nside::Int,
    angular_resolution_scale::Float64,
)
    # Shared-shell calibration now uses one symmetric angular-resolution scale
    # on the ideal reference limits instead of the old lower/upper band
    # expansion semantics used on the child/core path.
    return @timeg "diatomic.source.reference_band" begin
        angular_resolution_scale > 0.0 || throw(
            ArgumentError(
                "diatomic shared-shell angular-resolution scaling requires angular_resolution_scale > 0",
            ),
        )
        reference = _nested_diatomic_reference_band(
            basis,
            bundles,
            box;
            nside = nside,
            reference_fudge_factor = 1.0,
        )
        (
            reference_bounds = reference.reference_bounds,
            reference_retain = reference.reference_retain,
            ideal_theta_min = reference.ideal_theta_min,
            ideal_theta_max = reference.ideal_theta_max,
            theta_min = angular_resolution_scale * reference.ideal_theta_min,
            theta_max = angular_resolution_scale * reference.ideal_theta_max,
            angular_resolution_scale = angular_resolution_scale,
        )
    end
end

function _nested_diatomic_choose_candidate_from_stats(
    length_value::Int,
    candidate_stats::AbstractVector,
    reference_band,
)
    feasible = findall(
        stat -> stat.theta_min >= reference_band.theta_min && stat.theta_max <= reference_band.theta_max,
        candidate_stats,
    )
    if !isempty(feasible)
        chosen = candidate_stats[first(feasible)]
        return (
            retain = chosen.retained_count,
            mode = chosen.retained_count == length_value ? :direct : :doside,
            theta_min = chosen.theta_min,
            theta_max = chosen.theta_max,
            parent_limited = false,
            lower_band_limited = false,
            candidates = candidate_stats,
        )
    end
    upper_safe = findall(stat -> stat.theta_max <= reference_band.theta_max, candidate_stats)
    if !isempty(upper_safe)
        chosen = candidate_stats[first(upper_safe)]
        return (
            retain = chosen.retained_count,
            mode = chosen.retained_count == length_value ? :direct : :doside,
            theta_min = chosen.theta_min,
            theta_max = chosen.theta_max,
            parent_limited = false,
            lower_band_limited = true,
            candidates = candidate_stats,
        )
    end
    chosen = candidate_stats[end]
    return (
        retain = chosen.retained_count,
        mode = :direct,
        theta_min = chosen.theta_min,
        theta_max = chosen.theta_max,
        parent_limited = chosen.theta_max > reference_band.theta_max,
        lower_band_limited = false,
        candidates = candidate_stats,
    )
end

function _nested_diatomic_force_choice_direct(choice, length_value::Int)
    direct_candidates = [
        candidate for candidate in choice.candidates if candidate.retained_count == length_value
    ]
    direct_candidate = isempty(direct_candidates) ? choice.candidates[end] : direct_candidates[end]
    return (
        retain = direct_candidate.retained_count,
        mode = :direct,
        theta_min = direct_candidate.theta_min,
        theta_max = direct_candidate.theta_max,
        parent_limited = choice.parent_limited,
        lower_band_limited = choice.lower_band_limited,
        candidates = choice.candidates,
    )
end

function _nested_diatomic_choose_shell_axis_retain_count(
    basis,
    bundles::_CartesianNestedAxisBundles3D,
    line_axis::Symbol,
    line_interval::UnitRange{Int},
    boundary_box::NTuple{3,UnitRange{Int}},
    reference_band;
    minimum_retained::Int,
)
    candidate_stats = NamedTuple[]
    for count in _nested_diatomic_candidate_counts(length(line_interval), minimum_retained)
        points = _nested_diatomic_line_points_for_interval(
            bundles,
            line_axis,
            line_interval;
            retained_count = count,
        )
        stats = _nested_diatomic_aggregate_line_theta_stats(
            points,
            line_axis,
            _nested_diatomic_boundary_fixed_coord_pairs(bundles, boundary_box, line_axis),
            basis,
        )
        push!(
            candidate_stats,
            (
                retained_count = count,
                theta_min = stats.theta_min,
                theta_max = stats.theta_max,
            ),
        )
    end
    return _nested_diatomic_choose_candidate_from_stats(
        length(line_interval),
        candidate_stats,
        reference_band,
    )
end

function _nested_diatomic_adaptive_shell_retention(
    basis,
    bundles::_CartesianNestedAxisBundles3D,
    boundary_box::NTuple{3,UnitRange{Int}},
    inner_box::NTuple{3,UnitRange{Int}},
    retention::CartesianNestedCompleteShellRetentionContract;
    nside::Int,
    reference_fudge_factor::Union{Nothing,Float64} = nothing,
    shared_shell_angular_resolution_scale::Union{Nothing,Float64} = nothing,
)
    return @timeg "diatomic.source.adaptive_shell_retention" begin
        if isnothing(reference_fudge_factor) == isnothing(shared_shell_angular_resolution_scale)
            throw(
                ArgumentError(
                    "diatomic adaptive shell retention requires exactly one of reference_fudge_factor or shared_shell_angular_resolution_scale",
                ),
            )
        end
        reference =
            isnothing(shared_shell_angular_resolution_scale) ?
            _nested_diatomic_reference_band(
                basis,
                bundles,
                boundary_box;
                nside = nside,
                reference_fudge_factor = something(reference_fudge_factor),
            ) :
            _nested_diatomic_shared_shell_reference_band(
                basis,
                bundles,
                boundary_box;
                nside = nside,
                angular_resolution_scale = something(shared_shell_angular_resolution_scale),
            )
        chosen_x = @timeg "diatomic.source.axis_choice.x" begin
            _nested_diatomic_choose_shell_axis_retain_count(
                basis,
                bundles,
                :x,
                inner_box[1],
                boundary_box,
                reference;
                minimum_retained = _nested_diatomic_axis_minimum_retain(retention, :x),
            )
        end
        chosen_y = @timeg "diatomic.source.axis_choice.y" begin
            _nested_diatomic_choose_shell_axis_retain_count(
                basis,
                bundles,
                :y,
                inner_box[2],
                boundary_box,
                reference;
                minimum_retained = _nested_diatomic_axis_minimum_retain(retention, :y),
            )
        end
        chosen_z = @timeg "diatomic.source.axis_choice.z" begin
            _nested_diatomic_choose_shell_axis_retain_count(
                basis,
                bundles,
                :z,
                inner_box[3],
                boundary_box,
                reference;
                minimum_retained = _nested_diatomic_axis_minimum_retain(retention, :z),
            )
        end
        (
            reference = reference,
            chosen_x = chosen_x,
            chosen_y = chosen_y,
            chosen_z = chosen_z,
            retain_xy = (chosen_x.retain, chosen_y.retain),
            retain_xz = (chosen_x.retain, chosen_z.retain),
            retain_yz = (chosen_y.retain, chosen_z.retain),
            retain_x_edge = chosen_x.retain,
            retain_y_edge = chosen_y.retain,
            retain_z_edge = chosen_z.retain,
        )
    end
end

function _nested_diatomic_core_protected_radii(
    basis,
    bundles::_CartesianNestedAxisBundles3D,
    box::NTuple{3,UnitRange{Int}},
    bond_axis::Symbol,
    protect_rows::Int;
    atol::Float64 = 1.0e-10,
)
    protect_rows > 0 || return Float64[]
    centers_x = _nested_axis_pgdg(bundles, :x).centers
    centers_y = _nested_axis_pgdg(bundles, :y).centers
    centers_z = _nested_axis_pgdg(bundles, :z).centers
    axis_center_x = sum(Float64(nucleus[1]) for nucleus in basis.nuclei) / length(basis.nuclei)
    axis_center_y = sum(Float64(nucleus[2]) for nucleus in basis.nuclei) / length(basis.nuclei)
    axis_center_z = sum(Float64(nucleus[3]) for nucleus in basis.nuclei) / length(basis.nuclei)
    raw_distances =
        bond_axis == :x ? [hypot(Float64(centers_y[iy]) - axis_center_y, Float64(centers_z[iz]) - axis_center_z) for iy in box[2], iz in box[3]] :
        bond_axis == :y ? [hypot(Float64(centers_x[ix]) - axis_center_x, Float64(centers_z[iz]) - axis_center_z) for ix in box[1], iz in box[3]] :
        bond_axis == :z ? [hypot(Float64(centers_x[ix]) - axis_center_x, Float64(centers_y[iy]) - axis_center_y) for ix in box[1], iy in box[2]] :
        throw(
            ArgumentError(
                "diatomic core protection requires bond_axis = :x, :y, or :z",
            ),
        )
    distances = Float64[]
    for distance in sort(vec(Float64.(raw_distances)))
        if isempty(distances) || !isapprox(distance, last(distances); atol = atol, rtol = 0.0)
            push!(distances, distance)
        end
    end
    isempty(distances) && return Float64[]
    protected_count = isapprox(distances[1], 0.0; atol = atol, rtol = 0.0) ? protect_rows + 1 : protect_rows
    return distances[1:min(length(distances), protected_count)]
end

function _nested_diatomic_core_fixed_pair_radius(
    basis,
    bundles::_CartesianNestedAxisBundles3D,
    bond_axis::Symbol,
    fixed_indices::NTuple{2,Int},
)
    centers_x = _nested_axis_pgdg(bundles, :x).centers
    centers_y = _nested_axis_pgdg(bundles, :y).centers
    centers_z = _nested_axis_pgdg(bundles, :z).centers
    axis_center_x = sum(Float64(nucleus[1]) for nucleus in basis.nuclei) / length(basis.nuclei)
    axis_center_y = sum(Float64(nucleus[2]) for nucleus in basis.nuclei) / length(basis.nuclei)
    axis_center_z = sum(Float64(nucleus[3]) for nucleus in basis.nuclei) / length(basis.nuclei)
    return if bond_axis == :x
        hypot(
            Float64(centers_y[fixed_indices[1]]) - axis_center_y,
            Float64(centers_z[fixed_indices[2]]) - axis_center_z,
        )
    elseif bond_axis == :y
        hypot(
            Float64(centers_x[fixed_indices[1]]) - axis_center_x,
            Float64(centers_z[fixed_indices[2]]) - axis_center_z,
        )
    elseif bond_axis == :z
        hypot(
            Float64(centers_x[fixed_indices[1]]) - axis_center_x,
            Float64(centers_y[fixed_indices[2]]) - axis_center_y,
        )
    else
        throw(
            ArgumentError(
                "diatomic core fixed-pair radius requires bond_axis = :x, :y, or :z",
            ),
        )
    end
end

function _nested_diatomic_choose_fixed_parallel_line_retain_count(
    basis,
    bundles::_CartesianNestedAxisBundles3D,
    bond_axis::Symbol,
    parallel_interval::UnitRange{Int},
    fixed_indices::NTuple{2,Int},
    reference_band;
    minimum_retained::Int,
)
    line_axis = bond_axis
    fixed_coord_pairs =
        bond_axis == :x ?
        [(Float64(_nested_axis_pgdg(bundles, :y).centers[fixed_indices[1]]), Float64(_nested_axis_pgdg(bundles, :z).centers[fixed_indices[2]]))] :
        bond_axis == :y ?
        [(Float64(_nested_axis_pgdg(bundles, :x).centers[fixed_indices[1]]), Float64(_nested_axis_pgdg(bundles, :z).centers[fixed_indices[2]]))] :
        [(Float64(_nested_axis_pgdg(bundles, :x).centers[fixed_indices[1]]), Float64(_nested_axis_pgdg(bundles, :y).centers[fixed_indices[2]]))]
    candidate_stats = NamedTuple[]
    for count in _nested_diatomic_candidate_counts(length(parallel_interval), minimum_retained)
        points = _nested_diatomic_line_points_for_interval(
            bundles,
            line_axis,
            parallel_interval;
            retained_count = count,
        )
        stats = _nested_diatomic_aggregate_line_theta_stats(
            points,
            line_axis,
            fixed_coord_pairs,
            basis,
        )
        push!(
            candidate_stats,
            (
                retained_count = count,
                theta_min = stats.theta_min,
                theta_max = stats.theta_max,
            ),
        )
    end
    return _nested_diatomic_choose_candidate_from_stats(
        length(parallel_interval),
        candidate_stats,
        reference_band,
    )
end

function _nested_bond_aligned_diatomic_nonuniform_core_block(
    basis,
    bundles::_CartesianNestedAxisBundles3D,
    box::NTuple{3,UnitRange{Int}};
    bond_axis::Symbol,
    nside::Int,
    minimum_parallel_retain::Int,
    reference_fudge_factor::Float64,
    core_near_nucleus_protect_rows::Int,
)
    return @timeg "diatomic.source.nonuniform_core_block" begin
        reference = _nested_diatomic_reference_band(
            basis,
            bundles,
            box;
            nside = nside,
            reference_fudge_factor = reference_fudge_factor,
        )
        protected_radii = @timeg "diatomic.source.nonuniform_core.protected_rows" begin
            _nested_diatomic_core_protected_radii(
                basis,
                bundles,
                box,
                bond_axis,
                core_near_nucleus_protect_rows,
            )
        end
        dims = _nested_axis_lengths(bundles)
        coefficient_blocks = _CartesianCoefficientMap[]
        if bond_axis == :x
            @timeg "diatomic.source.nonuniform_core.line_emit.x" begin
                for iy in box[2], iz in box[3]
                    fixed_indices = (iy, iz)
                    choice = _nested_diatomic_choose_fixed_parallel_line_retain_count(
                        basis,
                        bundles,
                        bond_axis,
                        box[1],
                        fixed_indices,
                        reference;
                        minimum_retained = max(nside, minimum_parallel_retain),
                    )
                    radius = _nested_diatomic_core_fixed_pair_radius(basis, bundles, bond_axis, fixed_indices)
                    if any(level -> isapprox(radius, level; atol = 1.0e-10, rtol = 0.0), protected_radii)
                        choice = _nested_diatomic_force_choice_direct(choice, length(box[1]))
                    end
                    side = _nested_doside_1d(
                        _nested_axis_pgdg(bundles, :x),
                        box[1],
                        choice.retain;
                        enforce_symmetric_odd = false,
                    )
                    push!(
                        coefficient_blocks,
                        _nested_edge_product(:x, (:interior, :interior), side, fixed_indices, dims).coefficient_matrix,
                    )
                end
            end
        elseif bond_axis == :y
            @timeg "diatomic.source.nonuniform_core.line_emit.y" begin
                for ix in box[1], iz in box[3]
                    fixed_indices = (ix, iz)
                    choice = _nested_diatomic_choose_fixed_parallel_line_retain_count(
                        basis,
                        bundles,
                        bond_axis,
                        box[2],
                        fixed_indices,
                        reference;
                        minimum_retained = max(nside, minimum_parallel_retain),
                    )
                    radius = _nested_diatomic_core_fixed_pair_radius(basis, bundles, bond_axis, fixed_indices)
                    if any(level -> isapprox(radius, level; atol = 1.0e-10, rtol = 0.0), protected_radii)
                        choice = _nested_diatomic_force_choice_direct(choice, length(box[2]))
                    end
                    side = _nested_doside_1d(
                        _nested_axis_pgdg(bundles, :y),
                        box[2],
                        choice.retain;
                        enforce_symmetric_odd = false,
                    )
                    push!(
                        coefficient_blocks,
                        _nested_edge_product(:y, (:interior, :interior), side, fixed_indices, dims).coefficient_matrix,
                    )
                end
            end
        elseif bond_axis == :z
            @timeg "diatomic.source.nonuniform_core.line_emit.z" begin
                for ix in box[1], iy in box[2]
                    fixed_indices = (ix, iy)
                    choice = _nested_diatomic_choose_fixed_parallel_line_retain_count(
                        basis,
                        bundles,
                        bond_axis,
                        box[3],
                        fixed_indices,
                        reference;
                        minimum_retained = max(nside, minimum_parallel_retain),
                    )
                    radius = _nested_diatomic_core_fixed_pair_radius(basis, bundles, bond_axis, fixed_indices)
                    if any(level -> isapprox(radius, level; atol = 1.0e-10, rtol = 0.0), protected_radii)
                        choice = _nested_diatomic_force_choice_direct(choice, length(box[3]))
                    end
                    side = _nested_doside_1d(
                        _nested_axis_pgdg(bundles, :z),
                        box[3],
                        choice.retain;
                        enforce_symmetric_odd = false,
                    )
                    push!(
                        coefficient_blocks,
                        _nested_edge_product(:z, (:interior, :interior), side, fixed_indices, dims).coefficient_matrix,
                    )
                end
            end
        else
            throw(
                ArgumentError(
                    "diatomic nonuniform core construction requires bond_axis = :x, :y, or :z",
                ),
            )
        end
        (
            support_indices = _nested_box_support_indices(box..., dims),
            coefficient_matrix = (
                @timeg "diatomic.source.nonuniform_core.coefficient_merge" begin
                    _nested_hcat_coefficient_maps(coefficient_blocks)
                end
            ),
        )
    end
end

function _nested_bond_aligned_diatomic_sequence_for_box(
    basis,
    bundles::_CartesianNestedAxisBundles3D,
    box::NTuple{3,UnitRange{Int}},
    retention::CartesianNestedCompleteShellRetentionContract;
    bond_axis::Symbol,
    nside::Int,
    reference_fudge_factor::Float64,
    core_near_nucleus_protect_rows::Int,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
    packet_kernel::Symbol = :factorized_direct,
    build_packet::Bool = true,
)
    return @timeg "diatomic.source.child_sequence" begin
        isnothing(term_coefficients) && throw(
            ArgumentError("diatomic child-sequence assembly requires explicit term coefficients"),
        )
        current_box = box
        shell_layers = _CartesianNestedCompleteShell3D[]
        @timeg "diatomic.source.child_sequence.shell_construction" begin
            while minimum(length.(current_box)) > nside
                _nested_can_shrink_box(current_box) || break
                inner_box = _nested_inner_box(current_box)
                adaptive_retention = _nested_diatomic_adaptive_shell_retention(
                    basis,
                    bundles,
                    current_box,
                    inner_box,
                    retention;
                    nside = nside,
                    reference_fudge_factor = reference_fudge_factor,
                )
                push!(
                    shell_layers,
                    _nested_complete_rectangular_shell(
                        bundles,
                        inner_box...;
                        retain_xy = adaptive_retention.retain_xy,
                        retain_xz = adaptive_retention.retain_xz,
                        retain_yz = adaptive_retention.retain_yz,
                        retain_x_edge = adaptive_retention.retain_x_edge,
                        retain_y_edge = adaptive_retention.retain_y_edge,
                        retain_z_edge = adaptive_retention.retain_z_edge,
                        x_fixed = (first(current_box[1]), last(current_box[1])),
                        y_fixed = (first(current_box[2]), last(current_box[2])),
                        z_fixed = (first(current_box[3]), last(current_box[3])),
                        enforce_symmetric_odd = false,
                        term_storage = :compact_production,
                        term_coefficients = term_coefficients,
                        packet_kernel = packet_kernel,
                        verify_factorized_reconstruction = false,
                    ),
                )
                current_box = inner_box
            end
        end
        core_block = @timeg "diatomic.source.child_sequence.core_block" begin
            _nested_bond_aligned_diatomic_nonuniform_core_block(
                basis,
                bundles,
                current_box;
                bond_axis = bond_axis,
                nside = nside,
                minimum_parallel_retain = _nested_diatomic_axis_minimum_retain(retention, bond_axis),
                reference_fudge_factor = reference_fudge_factor,
                core_near_nucleus_protect_rows = core_near_nucleus_protect_rows,
            )
        end
        @timeg "diatomic.source.child_sequence.sequence_merge" begin
            _nested_shell_sequence_from_core_block(
                bundles,
                core_block.support_indices,
                core_block.coefficient_matrix,
                shell_layers,
                term_storage = :compact_production,
                term_coefficients = term_coefficients,
                packet_kernel = packet_kernel,
                build_packet = build_packet,
                verify_factorized_reconstruction = false,
            )
        end
    end
end

# Alg Nested-Diatomic step 5 and 6: Choose the bond-axis split plane at the
# parent-grid index nearest the midpoint, then reject it if the child boxes are
# too short or too thin in physical coordinates.
# See docs/src/algorithms/cartesian_nested_diatomic_box_policy.md.
function _nested_bond_aligned_diatomic_split_geometry(
    bundles::_CartesianNestedAxisBundles3D,
    parent_box::NTuple{3,UnitRange{Int}},
    working_box::NTuple{3,UnitRange{Int}};
    bond_axis::Symbol = :z,
    midpoint::Real = 0.0,
    nside::Int = 5,
    min_unsplit_parallel_to_transverse_ratio_for_split::Float64 = 3.0,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    use_midpoint_slab::Bool = true,
    prefer_midpoint_tie_side::Symbol = :left,
)
    return @timeg "diatomic.source.split_geometry" begin
        axis = bond_axis == :x ? 1 : bond_axis == :y ? 2 : bond_axis == :z ? 3 : 0
        axis != 0 || throw(ArgumentError("diatomic split geometry requires bond_axis = :x, :y, or :z"))
        min_unsplit_parallel_to_transverse_ratio_for_split > 0.0 || throw(
            ArgumentError(
                "diatomic split-eligibility guard requires min_unsplit_parallel_to_transverse_ratio_for_split > 0",
            ),
        )
        parallel_interval = working_box[axis]
        parallel_centers = _nested_axis_pgdg(bundles, bond_axis).centers
        working_widths = _nested_box_physical_widths(bundles, working_box)
        parallel_width = working_widths[axis]
        short_side_width = minimum(
            working_widths[index] for index in 1:3 if index != axis
        )
        use_slab = use_midpoint_slab && isodd(length(parallel_interval))
        split_index = use_slab ?
            _nested_diatomic_midpoint_row_index(parallel_centers, parallel_interval, midpoint) :
            _nested_diatomic_split_plane_index(
                parallel_centers,
                parallel_interval,
                midpoint;
                prefer_midpoint_tie_side = prefer_midpoint_tie_side,
            )
        left_box, midpoint_slab_box, right_box = if use_slab
            _nested_diatomic_midpoint_slab_split(working_box, bond_axis, split_index)
        else
            left_box, right_box = _nested_diatomic_child_boxes(working_box, bond_axis, split_index)
            (left_box, nothing, right_box)
        end
        child_boxes = [left_box, right_box]
        count_eligible =
            length(parallel_interval) > 2 * nside &&
            minimum(length(box[axis]) for box in child_boxes) >= nside
        unsplit_aspect_eligible =
            parallel_width > 0.0 &&
            short_side_width > 0.0 &&
            parallel_width > min_unsplit_parallel_to_transverse_ratio_for_split * short_side_width
        shape_eligible =
            count_eligible &&
            _nested_diatomic_children_are_roughly_cubic(
                bundles,
                child_boxes,
                bond_axis;
                min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
            )
        did_split = count_eligible && unsplit_aspect_eligible && shape_eligible
        _BondAlignedDiatomicSplitGeometry3D(
            parent_box,
            working_box,
            bond_axis,
            Float64(midpoint),
            split_index,
            count_eligible,
            unsplit_aspect_eligible,
            shape_eligible,
            did_split,
            did_split ? midpoint_slab_box : nothing,
            child_boxes,
            [_nested_box_physical_widths(bundles, box) for box in child_boxes],
        )
    end
end

function _nested_bond_aligned_diatomic_source(
    basis,
    bundles::_CartesianNestedAxisBundles3D;
    bond_axis::Symbol = :z,
    midpoint::Real = 0.0,
    nside::Int = 5,
    min_unsplit_parallel_to_transverse_ratio_for_split::Float64 = 3.0,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    reference_fudge_factor::Float64 = 1.2,
    shared_shell_angular_resolution_scale::Float64 = 1.4,
    core_near_nucleus_protect_rows::Union{Symbol,Integer} = :auto,
    use_midpoint_slab::Bool = true,
    prefer_midpoint_tie_side::Symbol = :left,
    shared_shell_retain_xy::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_xz::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_yz::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_xy::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_xz::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_yz::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_x_edge::Union{Nothing,Int} = nothing,
    retain_y_edge::Union{Nothing,Int} = nothing,
    retain_z_edge::Union{Nothing,Int} = nothing,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
    packet_kernel::Symbol = :factorized_direct,
)
    return @timeg "diatomic.source.total" begin
        isnothing(term_coefficients) && throw(
            ArgumentError("diatomic source assembly requires explicit term coefficients"),
        )
        child_retention = _nested_resolve_complete_shell_retention(
            nside;
            retain_xy = retain_xy,
            retain_xz = retain_xz,
            retain_yz = retain_yz,
            retain_x_edge = retain_x_edge,
            retain_y_edge = retain_y_edge,
            retain_z_edge = retain_z_edge,
        )
        shared_retention = _nested_resolve_complete_shell_retention(
            nside;
            retain_xy = shared_shell_retain_xy,
            retain_xz = shared_shell_retain_xz,
            retain_yz = shared_shell_retain_yz,
            retain_x_edge = child_retention.retain_x_edge,
            retain_y_edge = child_retention.retain_y_edge,
            retain_z_edge = child_retention.retain_z_edge,
        )
        protect_rows = _nested_diatomic_resolve_core_near_nucleus_protect_rows(
            core_near_nucleus_protect_rows,
            nside,
        )
        dims = _nested_axis_lengths(bundles)
        parent_box = (1:dims[1], 1:dims[2], 1:dims[3])
        shared_shell_layers = _CartesianNestedCompleteShell3D[]
        current_box = parent_box
        geometry = @timeg "diatomic.source.split_geometry.initial" begin
            _nested_bond_aligned_diatomic_split_geometry(
                bundles,
                parent_box,
                current_box;
                bond_axis = bond_axis,
                midpoint = midpoint,
                nside = nside,
                min_unsplit_parallel_to_transverse_ratio_for_split =
                    min_unsplit_parallel_to_transverse_ratio_for_split,
                min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
                use_midpoint_slab = use_midpoint_slab,
                prefer_midpoint_tie_side = prefer_midpoint_tie_side,
            )
        end

        @timeg "diatomic.source.shared_shell_construction" begin
            while true
                parallel_length = length(current_box[bond_axis == :x ? 1 : bond_axis == :y ? 2 : 3])
                if parallel_length <= 2 * nside || minimum(length.(current_box)) <= nside || !_nested_can_shrink_box(current_box)
                    break
                end
                inner_box = _nested_inner_box(current_box)
                adaptive_retention = _nested_diatomic_adaptive_shell_retention(
                    basis,
                    bundles,
                    current_box,
                    inner_box,
                    shared_retention;
                    nside = nside,
                    shared_shell_angular_resolution_scale = shared_shell_angular_resolution_scale,
                )
                push!(
                    shared_shell_layers,
                    _nested_complete_rectangular_shell(
                        bundles,
                        inner_box...;
                        retain_xy = adaptive_retention.retain_xy,
                        retain_xz = adaptive_retention.retain_xz,
                        retain_yz = adaptive_retention.retain_yz,
                        retain_x_edge = adaptive_retention.retain_x_edge,
                        retain_y_edge = adaptive_retention.retain_y_edge,
                        retain_z_edge = adaptive_retention.retain_z_edge,
                        x_fixed = (first(current_box[1]), last(current_box[1])),
                        y_fixed = (first(current_box[2]), last(current_box[2])),
                        z_fixed = (first(current_box[3]), last(current_box[3])),
                        enforce_symmetric_odd = false,
                        term_storage = :compact_production,
                        term_coefficients = term_coefficients,
                        packet_kernel = packet_kernel,
                        verify_factorized_reconstruction = false,
                    ),
                )
                current_box = inner_box
                geometry = @timeg "diatomic.source.split_geometry.rescan" begin
                    _nested_bond_aligned_diatomic_split_geometry(
                        bundles,
                        parent_box,
                        current_box;
                        bond_axis = bond_axis,
                        midpoint = midpoint,
                        nside = nside,
                        min_unsplit_parallel_to_transverse_ratio_for_split =
                            min_unsplit_parallel_to_transverse_ratio_for_split,
                        min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
                        use_midpoint_slab = use_midpoint_slab,
                        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
                    )
                end
                geometry.did_split && break
            end
        end

        child_sequences = _CartesianNestedShellSequence3D[]
        child_column_ranges = UnitRange{Int}[]
        midpoint_slab_column_range = nothing
        merged_sequence = nothing
        if geometry.did_split
            @timeg "diatomic.source.child_sequence_builds" begin
                for child_box in geometry.child_boxes
                    push!(
                        child_sequences,
                        _nested_bond_aligned_diatomic_sequence_for_box(
                            basis,
                            bundles,
                            child_box,
                            child_retention;
                            bond_axis = bond_axis,
                            nside = nside,
                            reference_fudge_factor = reference_fudge_factor,
                            core_near_nucleus_protect_rows = protect_rows,
                            term_coefficients = term_coefficients,
                            packet_kernel = packet_kernel,
                            build_packet = false,
                        ),
                    )
                end
            end
            merged_sequence = @timeg "diatomic.source.final_sequence_merge" begin
                core_support_blocks = Vector{Vector{Int}}()
                core_coefficient_blocks = _CartesianCoefficientMap[]
                push!(core_support_blocks, child_sequences[1].support_indices)
                push!(core_coefficient_blocks, child_sequences[1].coefficient_matrix)
                if !isnothing(geometry.shared_midpoint_box)
                    slab_data = _nested_direct_box_coefficients(bundles, geometry.shared_midpoint_box)
                    push!(core_support_blocks, slab_data.support_indices)
                    push!(core_coefficient_blocks, slab_data.coefficient_matrix)
                end
                push!(core_support_blocks, child_sequences[2].support_indices)
                push!(core_coefficient_blocks, child_sequences[2].coefficient_matrix)
                child_support = vcat(core_support_blocks...)
                child_coefficients = _nested_hcat_coefficient_maps(core_coefficient_blocks)
                _nested_shell_sequence_from_core_block(
                    bundles,
                    child_support,
                    child_coefficients,
                    shared_shell_layers,
                    term_storage = :compact_production,
                    term_coefficients = term_coefficients,
                    packet_kernel = packet_kernel,
                    verify_factorized_reconstruction = false,
                )
            end
            column_start = first(merged_sequence.core_column_range)
            left_columns = size(child_sequences[1].coefficient_matrix, 2)
            push!(child_column_ranges, column_start:(column_start + left_columns - 1))
            column_start = last(child_column_ranges[end]) + 1
            if !isnothing(geometry.shared_midpoint_box)
                slab_columns = prod(length.(geometry.shared_midpoint_box))
                midpoint_slab_column_range = column_start:(column_start + slab_columns - 1)
                column_start = last(midpoint_slab_column_range) + 1
            end
            right_columns = size(child_sequences[2].coefficient_matrix, 2)
            push!(child_column_ranges, column_start:(column_start + right_columns - 1))
        else
            shared_child = @timeg "diatomic.source.child_sequence_builds" begin
                _nested_bond_aligned_diatomic_sequence_for_box(
                    basis,
                    bundles,
                    current_box,
                    child_retention;
                    bond_axis = bond_axis,
                    nside = nside,
                    reference_fudge_factor = reference_fudge_factor,
                    core_near_nucleus_protect_rows = protect_rows,
                    term_coefficients = term_coefficients,
                    packet_kernel = packet_kernel,
                    build_packet = false,
                )
            end
            push!(child_sequences, shared_child)
            merged_sequence = @timeg "diatomic.source.final_sequence_merge" begin
                isempty(shared_shell_layers) ? shared_child :
                _nested_shell_sequence_from_core_block(
                    bundles,
                    shared_child.support_indices,
                    shared_child.coefficient_matrix,
                    shared_shell_layers,
                    term_storage = :compact_production,
                    term_coefficients = term_coefficients,
                    packet_kernel = packet_kernel,
                    verify_factorized_reconstruction = false,
                )
            end
            push!(child_column_ranges, merged_sequence.core_column_range)
        end

        _CartesianNestedBondAlignedDiatomicSource3D(
            basis,
            bundles,
            nside,
            child_retention,
            shared_retention,
            geometry,
            shared_shell_layers,
            child_sequences,
            child_column_ranges,
            midpoint_slab_column_range,
            merged_sequence,
        )
    end
end

function _nested_fixed_block(source::_CartesianNestedBondAlignedDiatomicSource3D)
    return @timeg "diatomic.fixed_block.contraction" begin
        _nested_fixed_block(source.sequence, source.basis)
    end
end
