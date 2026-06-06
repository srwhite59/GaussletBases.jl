"""
    _cartesian_terminal_shellification_geometry(parent_axes, nuclear_positions; kwargs...)

Private geometry-only terminal shellification for one atom or a bond-aligned
diatomic.

This helper owns only terminal shellification geometry and owned support. It
returns atom-local shells/cores, optional midpoint/contact slabs, shared
molecular shells, and outer mismatch slabs. It does not build coordinate
product box lowering objects, retained spaces, coefficient maps, operator
blocks, Hamiltonians, or artifacts.

`parent_axes` is an `NTuple{3}` of coordinate vectors, one per Cartesian axis.
Use `(collect(1:nx), collect(1:ny), collect(1:nz))` for index-only geometry.
`nuclear_positions` is either one 3-tuple position or an iterable of one or two
3-tuple positions in the coordinate system of `parent_axes`; each nucleus is
snapped to the nearest parent index along each axis.
"""
function _cartesian_terminal_shellification_geometry(
    parent_axes::NTuple{3,<:AbstractVector},
    nuclear_positions;
    core_side::Int = 5,
    q::Int = core_side,
    bond_axis::Symbol = :auto,
    audit_coverage::Bool = true,
)
    core_side > 0 && isodd(core_side) ||
        throw(ArgumentError("core_side must be a positive odd integer"))
    q > 0 || throw(ArgumentError("q must be a positive integer"))
    all(!isempty, parent_axes) ||
        throw(ArgumentError("parent axes must be nonempty"))

    parent_dims = Tuple(length(axis) for axis in parent_axes)
    parent_box = ntuple(axis -> 1:parent_dims[axis], 3)

    as_position_tuple(position) = begin
        length(position) == 3 ||
            throw(ArgumentError("each nuclear position must have length 3"))
        return (Float64(position[1]), Float64(position[2]), Float64(position[3]))
    end

    nuclei =
        nuclear_positions isa Tuple &&
        length(nuclear_positions) == 3 &&
        all(x -> x isa Real, nuclear_positions) ?
        (as_position_tuple(nuclear_positions),) :
        Tuple(as_position_tuple(position) for position in nuclear_positions)

    length(nuclei) in (1, 2) ||
        throw(
            ArgumentError(
                "_cartesian_terminal_shellification_geometry supports one atom or a diatomic",
            ),
        )

    nearest_index(axis_values, x) = begin
        distances = abs.(Float64.(axis_values) .- Float64(x))
        return findmin(distances)[2]
    end

    nuclear_indices = Tuple(
        ntuple(axis -> nearest_index(parent_axes[axis], nucleus[axis]), 3)
        for nucleus in nuclei
    )

    box_shape(box) = Tuple(length(box[axis]) for axis in 1:3)
    box_count(box) = prod(box_shape(box))
    same_box(a, b) = all(a[axis] == b[axis] for axis in 1:3)

    function box_inside_parent(box)
        return all(
            first(parent_box[axis]) <= first(box[axis]) &&
            last(box[axis]) <= last(parent_box[axis])
            for axis in 1:3
        )
    end

    function centered_box(center::NTuple{3,Int}, side::Int)
        radius = div(side - 1, 2)
        box = ntuple(axis -> (center[axis] - radius):(center[axis] + radius), 3)
        box_inside_parent(box) ||
            throw(
                ArgumentError(
                    "core box $box around nucleus index $center is outside parent box $parent_box",
                ),
            )
        return box
    end

    function can_expand(box)
        return all(
            first(parent_box[axis]) < first(box[axis]) &&
            last(box[axis]) < last(parent_box[axis])
            for axis in 1:3
        )
    end

    expand_box(box) =
        ntuple(axis -> (first(box[axis]) - 1):(last(box[axis]) + 1), 3)

    function hull_box(boxes...)
        return ntuple(
            axis -> begin
                lo = minimum(first(box[axis]) for box in boxes)
                hi = maximum(last(box[axis]) for box in boxes)
                lo:hi
            end,
            3,
        )
    end

    axis_symbol(axis) = (:x, :y, :z)[axis]
    shell_count(outer_box, inner_box) = box_count(outer_box) - box_count(inner_box)

    regions = NamedTuple[]

    function push_region!(;
        role::Symbol,
        region_kind::Symbol,
        outer_box,
        inner_box = nothing,
        atom_index = nothing,
        atom_side = nothing,
        shell_index = nothing,
        source = :_cartesian_terminal_shellification_geometry,
        metadata = (;),
    )
        support_count =
            isnothing(inner_box) ? box_count(outer_box) : shell_count(outer_box, inner_box)
        support_count > 0 ||
            throw(ArgumentError("terminal shellification region $role has nonpositive support"))
        push!(
            regions,
            (;
                object_kind = :cartesian_terminal_shellification_geometry_region,
                order_index = length(regions) + 1,
                terminal = true,
                role,
                region_kind,
                atom_index,
                atom_side,
                shell_index,
                outer_box,
                box = outer_box,
                inner_exclusion_box = inner_box,
                support_count,
                source,
                metadata = NamedTuple(metadata),
            ),
        )
    end

    function ordered_diatomic_indices(indices::Tuple)
        length(indices) == 2 ||
            throw(ArgumentError("ordered_diatomic_indices requires two nuclei"))
        axis = if bond_axis == :auto
            deltas = ntuple(a -> abs(indices[2][a] - indices[1][a]), 3)
            findmax(deltas)[2]
        elseif bond_axis == :x
            1
        elseif bond_axis == :y
            2
        elseif bond_axis == :z
            3
        else
            throw(ArgumentError("bond_axis must be :auto, :x, :y, or :z"))
        end

        for a in 1:3
            a == axis && continue
            indices[1][a] == indices[2][a] ||
                throw(
                    ArgumentError(
                        "diatomic nuclei must be bond-aligned on the parent grid; transverse index mismatch on axis $a",
                    ),
                )
        end

        if indices[1][axis] <= indices[2][axis]
            return axis, 1, 2
        else
            return axis, 2, 1
        end
    end

    gap_between(left_box, right_box, axis) =
        first(right_box[axis]) - last(left_box[axis]) - 1

    function central_gap_box_between(left_box, right_box, axis)
        gap = gap_between(left_box, right_box, axis)
        gap > 0 || return nothing
        ranges = Vector{UnitRange{Int}}(undef, 3)
        for a in 1:3
            if a == axis
                ranges[a] = (last(left_box[a]) + 1):(first(right_box[a]) - 1)
            else
                lo = min(first(left_box[a]), first(right_box[a]))
                hi = max(last(left_box[a]), last(right_box[a]))
                ranges[a] = lo:hi
            end
        end
        return Tuple(ranges)
    end

    function local_axis_spacing(axis_values, index::Int)
        n = length(axis_values)
        n == 1 && return 1.0
        if index == 1
            return abs(Float64(axis_values[2]) - Float64(axis_values[1]))
        elseif index == n
            return abs(Float64(axis_values[n]) - Float64(axis_values[n - 1]))
        else
            left = abs(Float64(axis_values[index]) - Float64(axis_values[index - 1]))
            right = abs(Float64(axis_values[index + 1]) - Float64(axis_values[index]))
            return (left + right) / 2
        end
    end

    function range_physical_size(axis_values, range::UnitRange{Int})
        length(range) == 0 && return 0.0
        if length(range) == 1
            return local_axis_spacing(axis_values, first(range))
        end

        lo = first(range)
        hi = last(range)
        span = abs(Float64(axis_values[hi]) - Float64(axis_values[lo]))
        left_spacing = local_axis_spacing(axis_values, lo)
        right_spacing = local_axis_spacing(axis_values, hi)
        return span + (left_spacing + right_spacing) / 2
    end

    function central_distorted_product_metadata(box, axis, gap_width)
        physical_sizes =
            ntuple(a -> range_physical_size(parent_axes[a], box[a]), 3)
        transverse_axes = Tuple(a for a in 1:3 if a != axis)
        transverse_size =
            (physical_sizes[transverse_axes[1]] + physical_sizes[transverse_axes[2]]) / 2
        transverse_size > 0 ||
            throw(
                ArgumentError(
                    "central distorted product box has nonpositive transverse size",
                ),
            )

        aspect_ratio = physical_sizes[axis] / transverse_size
        L = max(q, round(Int, q * aspect_ratio))
        source_mode_shape = ntuple(a -> a == axis ? L : q, 3)
        return (;
            bond_axis = axis_symbol(axis),
            central_gap_width = gap_width,
            q,
            L,
            source_mode_shape,
            physical_axis_sizes = physical_sizes,
            transverse_physical_size = transverse_size,
            aspect_ratio,
            lowering_hint = :distorted_comx_all_axes,
        )
    end

    function push_central_gap_regions!(left_box, right_box, axis)
        gap_width = gap_between(left_box, right_box, axis)
        gap_width == 0 && return nothing
        gap_width > 0 ||
            throw(ArgumentError("atom-local boxes overlap before central gap handling"))

        gap_box = central_gap_box_between(left_box, right_box, axis)
        if gap_width < q
            for (slice_index, slice) in enumerate(gap_box[axis])
                slab_box = ntuple(a -> a == axis ? (slice:slice) : gap_box[a], 3)
                push_region!(
                    role = :midpoint_slab,
                    region_kind = :direct_midpoint_slab,
                    outer_box = slab_box,
                    shell_index = slice_index,
                    metadata = (;
                        bond_axis = axis_symbol(axis),
                        central_gap_width = gap_width,
                        central_gap_slice_index = slice_index,
                    ),
                )
            end
        else
            push_region!(
                role = :central_distorted_product_box,
                region_kind = :central_distorted_product_box,
                outer_box = gap_box,
                metadata = central_distorted_product_metadata(
                    gap_box,
                    axis,
                    gap_width,
                ),
            )
        end
        return gap_box
    end

    function outer_mismatch_pieces(inner_box)
        same_box(inner_box, parent_box) && return ()
        pieces = NamedTuple[]

        # Axis-ordered disjoint decomposition of parent_box minus inner_box.
        # Earlier axes absorb edges/corners, so later slabs use the inner
        # intervals of earlier axes.
        for axis in 1:3
            low_range = first(parent_box[axis]):(first(inner_box[axis]) - 1)
            if !isempty(low_range)
                box = ntuple(
                    a -> begin
                        if a < axis
                            inner_box[a]
                        elseif a == axis
                            low_range
                        else
                            parent_box[a]
                        end
                    end,
                    3,
                )
                push!(
                    pieces,
                    (;
                        role =
                            Symbol((:x, :y, :z)[axis], "_low_outer_mismatch_slab"),
                        box,
                    ),
                )
            end

            high_range = (last(inner_box[axis]) + 1):last(parent_box[axis])
            if !isempty(high_range)
                box = ntuple(
                    a -> begin
                        if a < axis
                            inner_box[a]
                        elseif a == axis
                            high_range
                        else
                            parent_box[a]
                        end
                    end,
                    3,
                )
                push!(
                    pieces,
                    (;
                        role =
                            Symbol((:x, :y, :z)[axis], "_high_outer_mismatch_slab"),
                        box,
                    ),
                )
            end
        end
        return Tuple(pieces)
    end

    function grow_single_atom!(atom_index::Int, center::NTuple{3,Int})
        core = centered_box(center, core_side)
        push_region!(
            role = :atom_local_core,
            region_kind = :direct_core,
            outer_box = core,
            atom_index = atom_index,
            atom_side = :single,
            shell_index = 0,
        )

        current = core
        shell_index = 0
        while can_expand(current)
            next_box = expand_box(current)
            shell_index += 1
            push_region!(
                role = :atom_local_shell,
                region_kind = :complete_shell,
                outer_box = next_box,
                inner_box = current,
                atom_index = atom_index,
                atom_side = :single,
                shell_index = shell_index,
            )
            current = next_box
        end

        for piece in outer_mismatch_pieces(current)
            push_region!(
                role = piece.role,
                region_kind = :outer_mismatch_slab,
                outer_box = piece.box,
                atom_index = atom_index,
                atom_side = :single,
            )
        end
        return current
    end

    function grow_diatomic!()
        axis, left_atom, right_atom = ordered_diatomic_indices(nuclear_indices)

        left_core = centered_box(nuclear_indices[left_atom], core_side)
        right_core = centered_box(nuclear_indices[right_atom], core_side)
        initial_gap = gap_between(left_core, right_core, axis)
        initial_gap >= 0 ||
            throw(ArgumentError("initial atom core boxes overlap along the bond axis"))

        push_region!(
            role = :atom_local_core,
            region_kind = :direct_core,
            outer_box = left_core,
            atom_index = left_atom,
            atom_side = :left,
            shell_index = 0,
        )
        push_region!(
            role = :atom_local_core,
            region_kind = :direct_core,
            outer_box = right_core,
            atom_index = right_atom,
            atom_side = :right,
            shell_index = 0,
        )

        left_box = left_core
        right_box = right_core
        shell_index = 0

        while gap_between(left_box, right_box, axis) > 1 &&
              can_expand(left_box) &&
              can_expand(right_box)
            next_left = expand_box(left_box)
            next_right = expand_box(right_box)
            gap_between(next_left, next_right, axis) >= 0 || break

            shell_index += 1
            push_region!(
                role = :atom_local_shell,
                region_kind = :complete_shell,
                outer_box = next_left,
                inner_box = left_box,
                atom_index = left_atom,
                atom_side = :left,
                shell_index = shell_index,
            )
            push_region!(
                role = :atom_local_shell,
                region_kind = :complete_shell,
                outer_box = next_right,
                inner_box = right_box,
                atom_index = right_atom,
                atom_side = :right,
                shell_index = shell_index,
            )
            left_box = next_left
            right_box = next_right
        end

        central_gap = push_central_gap_regions!(left_box, right_box, axis)
        molecular_inner =
            isnothing(central_gap) ?
            hull_box(left_box, right_box) :
            hull_box(left_box, central_gap, right_box)

        current = molecular_inner
        shared_shell_index = 0
        while can_expand(current)
            next_box = expand_box(current)
            shared_shell_index += 1
            push_region!(
                role = :shared_molecular_shell,
                region_kind = :complete_shell,
                outer_box = next_box,
                inner_box = current,
                shell_index = shared_shell_index,
                metadata = (; bond_axis = axis_symbol(axis)),
            )
            current = next_box
        end

        for piece in outer_mismatch_pieces(current)
            push_region!(
                role = piece.role,
                region_kind = :outer_mismatch_slab,
                outer_box = piece.box,
                metadata = (; bond_axis = axis_symbol(axis)),
            )
        end

        return (; bond_axis = axis_symbol(axis), left_atom, right_atom)
    end

    function linear_indices_for_box(box)
        inds = Int[]
        sizehint!(inds, box_count(box))
        nx, ny, _ = parent_dims
        for k in box[3], j in box[2], i in box[1]
            push!(inds, i + (j - 1) * nx + (k - 1) * nx * ny)
        end
        return inds
    end

    function linear_indices_for_region(region)
        outer = Set(linear_indices_for_box(region.outer_box))
        if !isnothing(region.inner_exclusion_box)
            inner = Set(linear_indices_for_box(region.inner_exclusion_box))
            setdiff!(outer, inner)
        end
        return outer
    end

    function coverage_audit()
        expected = prod(parent_dims)
        seen = Set{Int}()
        duplicate_count = 0
        for region in regions
            inds = linear_indices_for_region(region)
            length(inds) == region.support_count ||
                throw(
                    ArgumentError(
                        "region $(region.order_index) support_count does not match explicit index count",
                    ),
                )
            for idx in inds
                if idx in seen
                    duplicate_count += 1
                else
                    push!(seen, idx)
                end
            end
        end
        return (;
            expected_parent_site_count = expected,
            covered_site_count = length(seen),
            region_support_count =
                sum(region.support_count for region in regions; init = 0),
            duplicate_site_count = duplicate_count,
            missing_site_count = expected - length(seen),
            coverage_complete = length(seen) == expected && duplicate_count == 0,
        )
    end

    system_kind =
        length(nuclear_indices) == 1 ? :one_center : :bond_aligned_diatomic

    diatomic_metadata = nothing
    resolved_bond_axis =
        if system_kind == :one_center
            grow_single_atom!(1, nuclear_indices[1])
            nothing
        else
            diatomic_metadata = grow_diatomic!()
            diatomic_metadata.bond_axis
        end

    regions_tuple = Tuple(regions)
    coverage = audit_coverage ? coverage_audit() : nothing
    if audit_coverage && !coverage.coverage_complete
        throw(
            ArgumentError(
                "terminal shellification coverage failed: " *
                "covered=$(coverage.covered_site_count), " *
                "expected=$(coverage.expected_parent_site_count), " *
                "duplicates=$(coverage.duplicate_site_count), " *
                "missing=$(coverage.missing_site_count)",
            ),
        )
    end

    return (;
        object_kind = :cartesian_terminal_shellification_geometry_plan,
        status = :terminal_regions_available,
        system_kind,
        parent_axes,
        parent_dims,
        parent_box,
        nuclear_positions = nuclei,
        nuclear_indices,
        bond_axis = resolved_bond_axis,
        core_side,
        q,
        regions = regions_tuple,
        region_count = length(regions_tuple),
        region_roles = Tuple(region.role for region in regions_tuple),
        region_kinds = Tuple(region.region_kind for region in regions_tuple),
        terminal_region_count = length(regions_tuple),
        aggregate_atom_boxes_emitted = false,
        diatomic_metadata,
        coverage,
        diagnostics = (;
            private_development_only = true,
            shellification_geometry_only = true,
            shellification_owns_geometry_and_owned_support_only = true,
            coordinate_product_box_lowering_materialized = false,
            retained_spaces_materialized = false,
            coefficient_maps_materialized = false,
            operator_blocks_materialized = false,
            pair_operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
            public_default_behavior_changed = false,
        ),
    )
end

function _cartesian_terminal_shellification_geometry_region_dependency(region)
    region.region_kind == :direct_core && return :plan_lowerable_direct_core
    if region.region_kind == :complete_shell
        region.role == :shared_molecular_shell &&
            return :plan_lowerable_shared_complete_shell
        return :plan_lowerable_complete_shell
    end
    region.region_kind == :direct_midpoint_slab &&
        return :plan_lowerable_direct_slab
    region.region_kind == :central_distorted_product_box &&
        return :planned_distorted_product_box_lowering
    region.region_kind == :outer_mismatch_slab &&
        return :plan_lowerable_direct_boundary_slab
    return :unsupported_terminal_shellification_region
end

function _cartesian_terminal_shellification_geometry_dependency_counts(dependencies)
    return (;
        plan_lowerable_direct_core_count =
            count(==(:plan_lowerable_direct_core), dependencies),
        plan_lowerable_complete_shell_count =
            count(==(:plan_lowerable_complete_shell), dependencies),
        plan_lowerable_shared_complete_shell_count =
            count(==(:plan_lowerable_shared_complete_shell), dependencies),
        plan_lowerable_direct_slab_count =
            count(==(:plan_lowerable_direct_slab), dependencies),
        planned_distorted_product_box_lowering_count =
            count(==(:planned_distorted_product_box_lowering), dependencies),
        plan_lowerable_direct_boundary_slab_count =
            count(==(:plan_lowerable_direct_boundary_slab), dependencies),
        unsupported_terminal_shellification_region_count =
            count(==(:unsupported_terminal_shellification_region), dependencies),
    )
end

function _cartesian_terminal_shellification_geometry_region_summary(region)
    dependency =
        _cartesian_terminal_shellification_geometry_region_dependency(region)
    return (;
        object_kind = :cartesian_terminal_shellification_geometry_region_summary,
        order_index = region.order_index,
        role = region.role,
        region_kind = region.region_kind,
        materialization_dependency = dependency,
        lowering_status = :planned_not_lowered,
        terminal = region.terminal,
        shellification_region_is_cpb = false,
        shellification_region_is_lowering_source = false,
        coordinate_product_box_lowering_materialized = false,
        support_count = region.support_count,
        outer_box = region.outer_box,
        inner_exclusion_box = region.inner_exclusion_box,
        metadata = region.metadata,
    )
end

function _cartesian_terminal_shellification_geometry_summary_coverage(plan)
    coverage = plan.coverage
    isnothing(coverage) && return (;
        coverage_status = :not_audited,
        coverage_complete = nothing,
        expected_parent_site_count = prod(plan.parent_dims),
        covered_site_count = nothing,
        region_support_count =
            sum(region.support_count for region in plan.regions; init = 0),
        duplicate_site_count = nothing,
        missing_site_count = nothing,
    )
    return (;
        coverage_status =
            coverage.coverage_complete ? :coverage_complete : :coverage_incomplete,
        coverage_complete = coverage.coverage_complete,
        expected_parent_site_count = coverage.expected_parent_site_count,
        covered_site_count = coverage.covered_site_count,
        region_support_count = coverage.region_support_count,
        duplicate_site_count = coverage.duplicate_site_count,
        missing_site_count = coverage.missing_site_count,
    )
end

function _cartesian_terminal_shellification_geometry_private_summary(plan)
    plan.object_kind == :cartesian_terminal_shellification_geometry_plan ||
        throw(
            ArgumentError(
                "terminal shellification geometry summary requires a cartesian_terminal_shellification_geometry_plan",
            ),
        )

    region_summaries = Tuple(
        _cartesian_terminal_shellification_geometry_region_summary(region)
        for region in plan.regions
    )
    dependencies =
        Tuple(region.materialization_dependency for region in region_summaries)
    central_midpoint_regions =
        Tuple(region for region in plan.regions if region.role == :midpoint_slab)
    central_distorted_regions =
        Tuple(
            region for region in plan.regions
            if region.role == :central_distorted_product_box
        )
    coverage_summary =
        _cartesian_terminal_shellification_geometry_summary_coverage(plan)

    return (;
        object_kind = :cartesian_terminal_shellification_geometry_private_summary,
        status = :metadata_only_summary_available,
        source_kind = :terminal_cartesian_shellification_geometry,
        source_object_kind = plan.object_kind,
        private_development_only = true,
        system_kind = plan.system_kind,
        parent_dims = plan.parent_dims,
        parent_box = plan.parent_box,
        nuclear_indices = plan.nuclear_indices,
        bond_axis = plan.bond_axis,
        core_side = plan.core_side,
        q = plan.q,
        ordered_terminal_region_roles =
            Tuple(region.role for region in plan.regions),
        ordered_terminal_region_kinds =
            Tuple(region.region_kind for region in plan.regions),
        ordered_materialization_dependencies = dependencies,
        region_count = plan.region_count,
        terminal_region_count = plan.terminal_region_count,
        region_support_counts =
            Tuple(region.support_count for region in plan.regions),
        total_support_count =
            sum(region.support_count for region in plan.regions; init = 0),
        materialization_dependency_counts =
            _cartesian_terminal_shellification_geometry_dependency_counts(
                dependencies,
            ),
        coverage = coverage_summary,
        coverage_status = coverage_summary.coverage_status,
        coverage_complete = coverage_summary.coverage_complete,
        central_gap_region_count =
            length(central_midpoint_regions) + length(central_distorted_regions),
        central_midpoint_slab_count = length(central_midpoint_regions),
        central_distorted_product_box_count = length(central_distorted_regions),
        central_distorted_product_box_metadata =
            Tuple(region.metadata for region in central_distorted_regions),
        regions = region_summaries,
        shellification_regions_are_cpbs = false,
        shellification_regions_are_lowering_sources = false,
        lowering_applied_by_summary = false,
        retained_spaces_materialized = false,
        coefficient_maps_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        public_default_behavior_changed = false,
    )
end

function _cartesian_terminal_shellification_geometry_region_lowering_family(
    region_summary,
)
    region_summary.region_kind == :direct_core && return :direct_product_core
    if region_summary.region_kind == :complete_shell
        region_summary.role == :shared_molecular_shell &&
            return :white_lindsey_adaptive_complete_shell
        return :white_lindsey_complete_shell
    end
    region_summary.region_kind == :direct_midpoint_slab &&
        return :direct_midpoint_slab
    region_summary.region_kind == :central_distorted_product_box &&
        return :distorted_comx_product_box_deferred
    region_summary.region_kind == :outer_mismatch_slab &&
        return :direct_boundary_slab
    return :unsupported_terminal_shellification_region
end

function _cartesian_terminal_shellification_geometry_independently_lowerable(
    dependency::Symbol,
)
    return dependency in (
        :plan_lowerable_direct_core,
        :plan_lowerable_complete_shell,
        :plan_lowerable_shared_complete_shell,
        :plan_lowerable_direct_slab,
        :plan_lowerable_direct_boundary_slab,
    )
end

function _cartesian_terminal_shellification_geometry_missing_lowering_reason(
    dependency::Symbol,
)
    dependency == :planned_distorted_product_box_lowering &&
        return :distorted_product_box_lowering_not_implemented
    dependency == :unsupported_terminal_shellification_region &&
        return :unsupported_terminal_shellification_region
    return nothing
end

function _cartesian_terminal_shellification_geometry_region_retirement_target(
    dependency::Symbol,
)
    _cartesian_terminal_shellification_geometry_independently_lowerable(
        dependency,
    ) && return :already_plan_lowered_region
    dependency == :planned_distorted_product_box_lowering &&
        return :pending_distorted_product_box_lowering_support
    return :requires_manager_review
end

function _cartesian_terminal_shellification_geometry_scaffold_region(
    region_summary,
)
    dependency = region_summary.materialization_dependency
    independently_lowerable =
        _cartesian_terminal_shellification_geometry_independently_lowerable(
            dependency,
        )
    return (;
        object_kind = :cartesian_terminal_shellification_scaffold_region,
        order_index = region_summary.order_index,
        role = region_summary.role,
        region_kind = region_summary.region_kind,
        box = region_summary.outer_box,
        outer_box = region_summary.outer_box,
        inner_exclusion_box = region_summary.inner_exclusion_box,
        support_count = region_summary.support_count,
        source_point_count = region_summary.support_count,
        materialization_dependency = dependency,
        lowering_family =
            _cartesian_terminal_shellification_geometry_region_lowering_family(
                region_summary,
            ),
        lowering_hint =
            hasproperty(region_summary.metadata, :lowering_hint) ?
            region_summary.metadata.lowering_hint :
            nothing,
        lowering_status = :planned_not_lowered,
        terminal = region_summary.terminal,
        shellification_region_is_cpb = false,
        shellification_region_is_lowering_source = false,
        coordinate_product_box_lowering_materialized = false,
        source_backed = false,
        independently_lowerable,
        missing_independent_lowering_reason =
            _cartesian_terminal_shellification_geometry_missing_lowering_reason(
                dependency,
            ),
        retirement_target =
            _cartesian_terminal_shellification_geometry_region_retirement_target(
                dependency,
            ),
        metadata = region_summary.metadata,
        provenance = (;
            source = :_cartesian_terminal_shellification_geometry_private_summary,
            summary_region_order_index = region_summary.order_index,
            terminal_region_role = region_summary.role,
            terminal_region_kind = region_summary.region_kind,
        ),
    )
end

function _cartesian_terminal_shellification_geometry_scaffold_materialization_status(
    regions,
)
    any(
        region ->
            region.materialization_dependency ==
            :unsupported_terminal_shellification_region,
        regions,
    ) && return :blocked_unsupported_regions
    any(!getproperty(region, :independently_lowerable) for region in regions) &&
        return :deferred_pending_distorted_product_box_lowering
    return :ready_supported_terminal_subset
end

function _cartesian_terminal_shellification_geometry_scaffold_coverage(summary)
    coverage = summary.coverage
    return (;
        object_kind = :cartesian_terminal_shellification_scaffold_coverage3d,
        expected_support_count = coverage.expected_parent_site_count,
        region_support_count = coverage.region_support_count,
        covered_support_count = coverage.covered_site_count,
        duplicate_count = coverage.duplicate_site_count,
        missing_count = coverage.missing_site_count,
        outside_count = 0,
        coverage_complete = coverage.coverage_complete,
    )
end

function _cartesian_terminal_shellification_geometry_spatial_policy_order(
    system_kind::Symbol,
)
    system_kind == :one_center && return :single_center_outward
    system_kind == :bond_aligned_diatomic && return :atom_outward
    return :terminal_geometry_order
end

function _cartesian_terminal_shellification_geometry_scaffold(
    plan;
    route_family::Symbol = :white_lindsey_low_order,
)
    summary =
        _cartesian_terminal_shellification_geometry_private_summary(plan)
    regions = Tuple(
        _cartesian_terminal_shellification_geometry_scaffold_region(region)
        for region in summary.regions
    )
    materialization_status =
        _cartesian_terminal_shellification_geometry_scaffold_materialization_status(
            regions,
        )
    dependency_counts = summary.materialization_dependency_counts
    deferred_regions =
        Tuple(region for region in regions if !region.independently_lowerable)
    coverage =
        _cartesian_terminal_shellification_geometry_scaffold_coverage(summary)

    return (;
        object_kind = :cartesian_terminal_shellification_scaffold3d,
        status = :planned_metadata_only,
        materialization_status,
        private_development_only = true,
        source_kind = :terminal_cartesian_shellification_geometry,
        route_family,
        system_classification = summary.system_kind,
        shellification_role = :terminal_cartesian_shellification_geometry,
        shellification_stage = :route_neutral_spatial_planning,
        lowering_stage = :not_lowered_by_shellification_plan,
        spatial_policy_order =
            _cartesian_terminal_shellification_geometry_spatial_policy_order(
                summary.system_kind,
            ),
        parent_box = summary.parent_box,
        working_box = summary.parent_box,
        full_parent_working_box = coverage.coverage_complete,
        parent_dims = summary.parent_dims,
        nuclear_indices = summary.nuclear_indices,
        bond_axis = summary.bond_axis,
        core_side = summary.core_side,
        q = summary.q,
        regions,
        region_count = length(regions),
        ordered_region_roles = Tuple(region.role for region in regions),
        ordered_region_kinds = Tuple(region.region_kind for region in regions),
        ordered_region_boxes = Tuple(region.box for region in regions),
        ordered_materialization_dependencies =
            Tuple(region.materialization_dependency for region in regions),
        materialization_dependency_counts = dependency_counts,
        independently_lowerable_region_count =
            count(region -> region.independently_lowerable, regions),
        deferred_region_count = length(deferred_regions),
        deferred_regions,
        unsupported_region_count =
            dependency_counts.unsupported_terminal_shellification_region_count,
        central_gap_region_count = summary.central_gap_region_count,
        central_midpoint_slab_count = summary.central_midpoint_slab_count,
        central_distorted_product_box_count =
            summary.central_distorted_product_box_count,
        central_distorted_product_box_metadata =
            summary.central_distorted_product_box_metadata,
        coverage,
        terminal_geometry_summary = summary,
        diagnostics = (;
            source = :_cartesian_terminal_shellification_geometry_scaffold,
            private_development_only = true,
            terminal_geometry_authority = true,
            route_neutral_spatial_planning = true,
            lowering_applied_by_plan = false,
            materialization_behavior_changed = false,
            public_default_behavior_changed = false,
            shellification_regions_are_cpbs = false,
            shellification_regions_are_lowering_sources = false,
            retained_spaces_materialized = false,
            coefficient_maps_materialized = false,
            operator_blocks_materialized = false,
            pair_operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
            cartesian_shells_behavior_changed = false,
        ),
    )
end

function _cartesian_terminal_region_unit_key(region)
    return Symbol(
        "terminal_region_",
        string(region.order_index),
        "_",
        String(region.role),
    )
end

function _cartesian_terminal_region_unit_role(region)
    return Symbol(String(region.role), "_unit")
end

function _cartesian_terminal_region_unit_required_metadata(region, field::Symbol)
    hasproperty(region.metadata, field) ||
        throw(
            ArgumentError(
                "terminal region $(region.role) is missing required metadata field $field",
            ),
        )
    return getproperty(region.metadata, field)
end

function _cartesian_terminal_region_complete_shell_lw_metadata(region)
    isnothing(region.inner_exclusion_box) &&
        throw(ArgumentError("complete shell terminal region is missing inner exclusion box"))
    return (;
        status = :planned_not_materialized,
        lowering_family_planned = :white_lindsey_boundary_stratum_cpbs,
        enumeration_policy = :complete_shell_boundary_strata_metadata,
        facet_count = 6,
        edge_count = 12,
        corner_count = 8,
        total_source_cpb_count = 26,
        source_cpbs_materialized = false,
        coefficient_maps_materialized = false,
        transform_contracts_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
    )
end

function _cartesian_terminal_region_complete_shell_pqs_metadata(region)
    isnothing(region.inner_exclusion_box) &&
        throw(ArgumentError("complete shell terminal region is missing inner exclusion box"))
    return (;
        status = :planned_not_materialized,
        lowering_family_planned = :pqs_filled_source_cpb,
        source_cpb_plan_status = :planned_not_materialized,
        source_cpb_plan_box = region.outer_box,
        source_cpb_plan_kind = :filled_source_cpb,
        source_cpb_count = 1,
        retained_rule = :boundary_comx_product_mode_selection,
        intermediate_retained_space_status = :not_materialized,
        shell_realization_status = :projection_lowdin_planned_not_materialized,
        face_edge_corner_decomposition_required = false,
        coefficient_maps_materialized = false,
        transform_contracts_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
    )
end

function _cartesian_terminal_region_direct_source_metadata(
    region,
    lowering_family_planned::Symbol,
)
    return (;
        status = :planned_not_materialized,
        lowering_family_planned,
        source_cpb_plan_status = :planned_not_materialized,
        source_cpb_plan_box = region.outer_box,
        source_cpb_plan_kind = :direct_coordinate_product_source,
        source_cpb_plan_equals_owned_support = true,
        source_cpb_count = 1,
        retained_rule = :direct_source_modes,
        coefficient_maps_materialized = false,
        transform_contracts_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
    )
end

function _cartesian_terminal_region_unit_mapping(region)
    if region.region_kind == :direct_core
        return (;
            unit_kind = :direct_core_unit,
            lowering_family_planned = :direct_core_identity_cpb,
            identity_lowering_planned = true,
            owned_support_is_cpb = false,
            owned_support_status = :owned_support_equals_region_box,
            source_cpb_plan =
                _cartesian_terminal_region_direct_source_metadata(
                    region,
                    :direct_core_identity_cpb,
                ),
            lw_complete_shell_lowering = nothing,
            pqs_complete_shell_lowering = nothing,
        )
    end
    if region.region_kind == :complete_shell
        unit_kind =
            region.role == :shared_molecular_shell ?
            :shared_molecular_complete_shell_unit :
            :atom_local_complete_shell_unit
        return (;
            unit_kind,
            lowering_family_planned = :complete_shell_lowering_recipe_choice_pending,
            identity_lowering_planned = false,
            owned_support_is_cpb = false,
            owned_support_status = :owned_shell_support_outer_minus_inner_exclusion,
            source_cpb_plan = nothing,
            lw_complete_shell_lowering =
                _cartesian_terminal_region_complete_shell_lw_metadata(region),
            pqs_complete_shell_lowering =
                _cartesian_terminal_region_complete_shell_pqs_metadata(region),
        )
    end
    if region.region_kind == :direct_midpoint_slab
        return (;
            unit_kind = :direct_midpoint_slab_unit,
            lowering_family_planned = :direct_slab_identity_cpb,
            identity_lowering_planned = true,
            owned_support_is_cpb = false,
            owned_support_status = :owned_support_equals_region_box,
            source_cpb_plan =
                _cartesian_terminal_region_direct_source_metadata(
                    region,
                    :direct_slab_identity_cpb,
                ),
            lw_complete_shell_lowering = nothing,
            pqs_complete_shell_lowering = nothing,
        )
    end
    if region.region_kind == :central_distorted_product_box
        source_mode_shape =
            _cartesian_terminal_region_unit_required_metadata(
                region,
                :source_mode_shape,
            )
        q = _cartesian_terminal_region_unit_required_metadata(region, :q)
        L = _cartesian_terminal_region_unit_required_metadata(region, :L)
        aspect_ratio =
            _cartesian_terminal_region_unit_required_metadata(region, :aspect_ratio)
        return (;
            unit_kind = :central_distorted_product_box_unit,
            lowering_family_planned = :distorted_product_box_comx,
            identity_lowering_planned = false,
            owned_support_is_cpb = false,
            owned_support_status = :owned_support_equals_region_box,
            source_cpb_plan = merge(
                _cartesian_terminal_region_direct_source_metadata(
                    region,
                    :distorted_product_box_comx,
                ),
                (;
                    source_cpb_plan_kind = :central_distorted_product_box_source,
                    retained_rule = :distorted_product_comx_all_axes,
                    source_mode_shape,
                    q,
                    L,
                    aspect_ratio,
                ),
            ),
            lw_complete_shell_lowering = nothing,
            pqs_complete_shell_lowering = nothing,
        )
    end
    if region.region_kind == :outer_mismatch_slab
        return (;
            unit_kind = :outer_mismatch_boundary_slab_unit,
            lowering_family_planned = :direct_boundary_slab_identity_cpb,
            identity_lowering_planned = true,
            owned_support_is_cpb = false,
            owned_support_status = :owned_support_equals_region_box,
            source_cpb_plan =
                _cartesian_terminal_region_direct_source_metadata(
                    region,
                    :direct_boundary_slab_identity_cpb,
                ),
            lw_complete_shell_lowering = nothing,
            pqs_complete_shell_lowering = nothing,
        )
    end
    throw(
        ArgumentError(
            "unsupported terminal shellification region kind $(region.region_kind)",
        ),
    )
end

function _cartesian_terminal_region_unit_record(region)
    mapping = _cartesian_terminal_region_unit_mapping(region)
    source_mode_shape =
        region.region_kind == :central_distorted_product_box ?
        mapping.source_cpb_plan.source_mode_shape :
        nothing
    q =
        region.region_kind == :central_distorted_product_box ?
        mapping.source_cpb_plan.q :
        nothing
    L =
        region.region_kind == :central_distorted_product_box ?
        mapping.source_cpb_plan.L :
        nothing
    aspect_ratio =
        region.region_kind == :central_distorted_product_box ?
        mapping.source_cpb_plan.aspect_ratio :
        nothing

    return (;
        object_kind = :cartesian_terminal_region_unit_record,
        unit_index = region.order_index,
        unit_key = _cartesian_terminal_region_unit_key(region),
        unit_role = _cartesian_terminal_region_unit_role(region),
        unit_kind = mapping.unit_kind,
        terminal_region_order_index = region.order_index,
        terminal_region_role = region.role,
        terminal_region_kind = region.region_kind,
        outer_box = region.outer_box,
        box = region.outer_box,
        inner_exclusion_box = region.inner_exclusion_box,
        support_count = region.support_count,
        owned_support_status = mapping.owned_support_status,
        owned_support_is_cpb = mapping.owned_support_is_cpb,
        shellification_region_is_cpb = false,
        shellification_region_is_lowering_source = false,
        lowering_family_planned = mapping.lowering_family_planned,
        lowering_recipe_status = :planned_not_materialized,
        source_cpb_plan_status = :planned_not_materialized,
        source_cpb_plan = mapping.source_cpb_plan,
        lw_complete_shell_lowering = mapping.lw_complete_shell_lowering,
        pqs_complete_shell_lowering = mapping.pqs_complete_shell_lowering,
        pqs_filled_source_cpb_plan =
            isnothing(mapping.pqs_complete_shell_lowering) ?
            nothing :
            mapping.pqs_complete_shell_lowering,
        source_mode_shape,
        q,
        L,
        aspect_ratio,
        identity_lowering_planned = mapping.identity_lowering_planned,
        retained_space_status = :not_materialized,
        coefficient_maps_materialized = false,
        transform_contracts_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        final_retained_unit_status =
            :not_available_terminal_region_metadata_only,
        metadata = region.metadata,
        provenance = (;
            source = :_cartesian_terminal_shellification_region_unit_inventory,
            scaffold_region_order_index = region.order_index,
            terminal_region_role = region.role,
            terminal_region_kind = region.region_kind,
        ),
    )
end

function _cartesian_terminal_shellification_region_unit_inventory(scaffold)
    scaffold.object_kind == :cartesian_terminal_shellification_scaffold3d ||
        throw(
            ArgumentError(
                "terminal-region unit inventory requires a cartesian_terminal_shellification_scaffold3d",
            ),
        )
    records = Tuple(
        _cartesian_terminal_region_unit_record(region)
        for region in scaffold.regions
    )
    complete_shell_records =
        Tuple(record for record in records if record.terminal_region_kind == :complete_shell)
    lw_complete_shell_cpb_count =
        sum(
            record -> record.lw_complete_shell_lowering.total_source_cpb_count,
            complete_shell_records;
            init = 0,
        )

    return (;
        object_kind = :cartesian_terminal_region_unit_inventory,
        status = :available_terminal_region_unit_inventory,
        inventory_source = :terminal_shellification_scaffold,
        private_development_only = true,
        terminal_region_count = scaffold.region_count,
        unit_count = length(records),
        terminal_region_units = records,
        unit_records = records,
        unit_keys = Tuple(record.unit_key for record in records),
        unit_roles = Tuple(record.unit_role for record in records),
        unit_kinds = Tuple(record.unit_kind for record in records),
        terminal_region_roles =
            Tuple(record.terminal_region_role for record in records),
        terminal_region_kinds =
            Tuple(record.terminal_region_kind for record in records),
        support_counts = Tuple(record.support_count for record in records),
        all_units_from_terminal_regions = length(records) == scaffold.region_count,
        complete_shell_unit_count = length(complete_shell_records),
        lw_complete_shell_cpb_enumeration_available =
            !isempty(complete_shell_records),
        lw_complete_shell_region_count = length(complete_shell_records),
        lw_complete_shell_cpb_count,
        lw_complete_shell_cpb_family_counts =
            (
                facet_cpb = 6 * length(complete_shell_records),
                edge_cpb = 12 * length(complete_shell_records),
                corner_cpb = 8 * length(complete_shell_records),
            ),
        final_retained_unit_inventory_available = false,
        pair_inventory_available = false,
        pair_inventory_status = :not_available_terminal_region_metadata_only,
        coefficient_maps_materialized = false,
        transform_contracts_materialized = false,
        retained_spaces_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        diagnostics = (;
            source = :_cartesian_terminal_shellification_region_unit_inventory,
            private_development_only = true,
            terminal_region_metadata_only = true,
            all_units_from_terminal_regions =
                length(records) == scaffold.region_count,
            no_aggregate_atom_box_units =
                all(
                    !(record.unit_key in (:left_atom_box, :right_atom_box))
                    for record in records
                ),
            no_atom_growth_route_units =
                all(
                    record.lowering_family_planned !=
                    :white_lindsey_atom_local_child_shellification
                    for record in records
                ),
            shellification_regions_are_cpbs = false,
            final_retained_unit_inventory_available = false,
            pair_inventory_available = false,
            materialization_behavior_changed = false,
            public_default_behavior_changed = false,
        ),
    )
end

function _cartesian_terminal_region_lowering_contract_base(
    unit,
    contract_index::Int,
    lowering_contract_kind::Symbol,
)
    return (;
        object_kind = :cartesian_terminal_region_lowering_contract,
        status = :planned_not_materialized,
        contract_index,
        contract_key =
            Symbol(String(unit.unit_key), "_", String(lowering_contract_kind)),
        unit_index = unit.unit_index,
        unit_key = unit.unit_key,
        unit_role = unit.unit_role,
        unit_kind = unit.unit_kind,
        terminal_region_order_index = unit.terminal_region_order_index,
        terminal_region_role = unit.terminal_region_role,
        terminal_region_kind = unit.terminal_region_kind,
        lowering_contract_kind,
        outer_box = unit.outer_box,
        box = unit.outer_box,
        inner_exclusion_box = unit.inner_exclusion_box,
        support_count = unit.support_count,
        owned_support_status = unit.owned_support_status,
        owned_support_is_cpb = unit.owned_support_is_cpb,
        shellification_region_is_cpb = unit.shellification_region_is_cpb,
        shellification_region_is_lowering_source =
            unit.shellification_region_is_lowering_source,
        coefficient_maps_materialized = false,
        transform_contracts_materialized = false,
        retained_spaces_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        final_retained_unit_status = :not_materialized,
        final_retained_unit_records_materialized = false,
        metadata = unit.metadata,
        provenance = (;
            source = :_cartesian_terminal_region_lowering_contract_inventory,
            unit_inventory_source =
                :_cartesian_terminal_shellification_region_unit_inventory,
            unit_key = unit.unit_key,
            terminal_region_role = unit.terminal_region_role,
            terminal_region_kind = unit.terminal_region_kind,
        ),
    )
end

function _cartesian_terminal_region_direct_lowering_contract(
    unit,
    contract_index::Int,
)
    plan = unit.source_cpb_plan
    return merge(
        _cartesian_terminal_region_lowering_contract_base(
            unit,
            contract_index,
            unit.lowering_family_planned,
        ),
        (;
            source_cpb_plan_status = plan.source_cpb_plan_status,
            source_cpb_plan_box = plan.source_cpb_plan_box,
            source_cpb_plan_kind = plan.source_cpb_plan_kind,
            source_cpb_plan_equals_owned_support =
                plan.source_cpb_plan_equals_owned_support,
            source_cpb_count = plan.source_cpb_count,
            source_cpb_family_counts =
                (facet_cpb = 0, edge_cpb = 0, corner_cpb = 0),
            identity_like_source_contract = true,
            retained_rule = plan.retained_rule,
            intermediate_retained_space_status = :not_materialized,
            shell_realization_status = :not_materialized,
            face_edge_corner_decomposition_required = false,
            final_unit_count_planned = 1,
            final_unit_granularity = :single_direct_source_cpb,
        ),
    )
end

function _cartesian_terminal_region_lw_complete_shell_lowering_contract(
    unit,
    contract_index::Int,
)
    lowering = unit.lw_complete_shell_lowering
    return merge(
        _cartesian_terminal_region_lowering_contract_base(
            unit,
            contract_index,
            :white_lindsey_boundary_strata,
        ),
        (;
            source_cpb_plan_status = :planned_not_materialized,
            source_cpb_plan_box = nothing,
            source_cpb_plan_kind = :complete_shell_boundary_strata,
            source_cpb_plan_equals_owned_support = false,
            source_cpb_count = lowering.total_source_cpb_count,
            source_cpb_family_counts =
                (
                    facet_cpb = lowering.facet_count,
                    edge_cpb = lowering.edge_count,
                    corner_cpb = lowering.corner_count,
                ),
            source_cpbs_materialized = false,
            identity_like_source_contract = false,
            retained_rule = :white_lindsey_boundary_strata_children,
            intermediate_retained_space_status = :not_materialized,
            shell_realization_status = :not_materialized,
            face_edge_corner_decomposition_required = true,
            final_unit_count_planned = lowering.total_source_cpb_count,
            final_unit_granularity =
                :white_lindsey_boundary_strata_children_planned,
        ),
    )
end

function _cartesian_terminal_region_pqs_complete_shell_lowering_contract(
    unit,
    contract_index::Int,
)
    lowering = unit.pqs_complete_shell_lowering
    return merge(
        _cartesian_terminal_region_lowering_contract_base(
            unit,
            contract_index,
            :pqs_filled_source_cpb,
        ),
        (;
            source_cpb_plan_status = lowering.source_cpb_plan_status,
            source_cpb_plan_box = lowering.source_cpb_plan_box,
            source_cpb_plan_kind = lowering.source_cpb_plan_kind,
            source_cpb_plan_equals_owned_support = false,
            source_cpb_count = lowering.source_cpb_count,
            source_cpb_family_counts =
                (facet_cpb = 0, edge_cpb = 0, corner_cpb = 0),
            identity_like_source_contract = false,
            retained_rule = lowering.retained_rule,
            intermediate_retained_space_status = :planned_not_materialized,
            shell_realization_status = lowering.shell_realization_status,
            face_edge_corner_decomposition_required =
                lowering.face_edge_corner_decomposition_required,
            final_unit_count_planned = 1,
            final_unit_granularity = :single_pqs_shell_realized_unit_planned,
        ),
    )
end

function _cartesian_terminal_region_distorted_product_box_lowering_contract(
    unit,
    contract_index::Int,
)
    plan = unit.source_cpb_plan
    return merge(
        _cartesian_terminal_region_lowering_contract_base(
            unit,
            contract_index,
            :distorted_product_box_comx,
        ),
        (;
            source_cpb_plan_status = plan.source_cpb_plan_status,
            source_cpb_plan_box = plan.source_cpb_plan_box,
            source_cpb_plan_kind = plan.source_cpb_plan_kind,
            source_cpb_plan_equals_owned_support =
                plan.source_cpb_plan_equals_owned_support,
            source_cpb_count = plan.source_cpb_count,
            source_cpb_family_counts =
                (facet_cpb = 0, edge_cpb = 0, corner_cpb = 0),
            identity_like_source_contract = false,
            retained_rule = plan.retained_rule,
            source_mode_shape = plan.source_mode_shape,
            q = plan.q,
            L = plan.L,
            aspect_ratio = plan.aspect_ratio,
            intermediate_retained_space_status = :planned_not_materialized,
            shell_realization_status = :not_materialized,
            face_edge_corner_decomposition_required = false,
            final_unit_count_planned = 1,
            final_unit_granularity =
                :single_distorted_product_box_unit_planned,
        ),
    )
end

function _cartesian_terminal_region_lowering_contracts_for_unit(
    unit,
    first_contract_index::Int,
)
    records = NamedTuple[]

    if unit.terminal_region_kind in (
        :direct_core,
        :direct_midpoint_slab,
        :outer_mismatch_slab,
    )
        push!(
            records,
            _cartesian_terminal_region_direct_lowering_contract(
                unit,
                first_contract_index,
            ),
        )
    elseif unit.terminal_region_kind == :complete_shell
        push!(
            records,
            _cartesian_terminal_region_lw_complete_shell_lowering_contract(
                unit,
                first_contract_index,
            ),
        )
        push!(
            records,
            _cartesian_terminal_region_pqs_complete_shell_lowering_contract(
                unit,
                first_contract_index + 1,
            ),
        )
    elseif unit.terminal_region_kind == :central_distorted_product_box
        push!(
            records,
            _cartesian_terminal_region_distorted_product_box_lowering_contract(
                unit,
                first_contract_index,
            ),
        )
    else
        throw(
            ArgumentError(
                "unsupported terminal-region unit kind $(unit.terminal_region_kind)",
            ),
        )
    end

    return Tuple(records)
end

function _cartesian_terminal_region_lowering_contract_kind_counts(contracts)
    return (;
        direct_core_identity_cpb_count =
            count(
                contract ->
                    contract.lowering_contract_kind == :direct_core_identity_cpb,
                contracts,
            ),
        direct_slab_identity_cpb_count =
            count(
                contract ->
                    contract.lowering_contract_kind == :direct_slab_identity_cpb,
                contracts,
            ),
        direct_boundary_slab_identity_cpb_count =
            count(
                contract ->
                    contract.lowering_contract_kind ==
                    :direct_boundary_slab_identity_cpb,
                contracts,
            ),
        white_lindsey_boundary_strata_count =
            count(
                contract ->
                    contract.lowering_contract_kind ==
                    :white_lindsey_boundary_strata,
                contracts,
            ),
        pqs_filled_source_cpb_count =
            count(
                contract ->
                    contract.lowering_contract_kind == :pqs_filled_source_cpb,
                contracts,
            ),
        distorted_product_box_comx_count =
            count(
                contract ->
                    contract.lowering_contract_kind == :distorted_product_box_comx,
                contracts,
            ),
    )
end

function _cartesian_terminal_region_lowering_contract_inventory(unit_inventory)
    unit_inventory.object_kind == :cartesian_terminal_region_unit_inventory ||
        throw(
            ArgumentError(
                "terminal-region lowering contracts require a cartesian_terminal_region_unit_inventory",
            ),
        )

    units = unit_inventory.terminal_region_units
    contracts = NamedTuple[]
    for unit in units
        append!(
            contracts,
            _cartesian_terminal_region_lowering_contracts_for_unit(
                unit,
                length(contracts) + 1,
            ),
        )
    end
    lowering_contracts = Tuple(contracts)
    contract_counts_by_unit = Tuple(
        (
            unit_key = unit.unit_key,
            lowering_contract_count =
                count(contract -> contract.unit_key == unit.unit_key, lowering_contracts),
        ) for unit in units
    )
    all(entry -> entry.lowering_contract_count >= 1, contract_counts_by_unit) ||
        throw(ArgumentError("every terminal-region unit must have a lowering contract"))

    complete_shell_contracts = Tuple(
        contract for contract in lowering_contracts
        if contract.terminal_region_kind == :complete_shell
    )
    lw_contracts = Tuple(
        contract for contract in lowering_contracts
        if contract.lowering_contract_kind == :white_lindsey_boundary_strata
    )

    return (;
        object_kind = :cartesian_terminal_region_lowering_contract_inventory,
        status = :available_terminal_region_lowering_contract_inventory,
        inventory_source = :terminal_region_unit_inventory,
        source_object_kind = unit_inventory.object_kind,
        private_development_only = true,
        terminal_region_unit_count = unit_inventory.unit_count,
        lowering_contract_count = length(lowering_contracts),
        lowering_contracts,
        contract_records = lowering_contracts,
        unit_keys = unit_inventory.unit_keys,
        unit_roles = unit_inventory.unit_roles,
        unit_kinds = unit_inventory.unit_kinds,
        terminal_region_roles = unit_inventory.terminal_region_roles,
        terminal_region_kinds = unit_inventory.terminal_region_kinds,
        support_counts = unit_inventory.support_counts,
        contract_counts_by_unit,
        lowering_contract_kinds =
            Tuple(contract.lowering_contract_kind for contract in lowering_contracts),
        lowering_contract_kind_counts =
            _cartesian_terminal_region_lowering_contract_kind_counts(
                lowering_contracts,
            ),
        complete_shell_unit_count = unit_inventory.complete_shell_unit_count,
        complete_shell_lowering_contract_count = length(complete_shell_contracts),
        lw_complete_shell_lowering_contract_count = length(lw_contracts),
        lw_complete_shell_cpb_count =
            sum(contract -> contract.source_cpb_count, lw_contracts; init = 0),
        lw_complete_shell_cpb_family_counts =
            (
                facet_cpb =
                    sum(
                        contract -> contract.source_cpb_family_counts.facet_cpb,
                        lw_contracts;
                        init = 0,
                    ),
                edge_cpb =
                    sum(
                        contract -> contract.source_cpb_family_counts.edge_cpb,
                        lw_contracts;
                        init = 0,
                    ),
                corner_cpb =
                    sum(
                        contract -> contract.source_cpb_family_counts.corner_cpb,
                        lw_contracts;
                        init = 0,
                    ),
            ),
        final_retained_unit_inventory_available = false,
        pair_inventory_available = false,
        pair_inventory_status = :not_available_lowering_contract_metadata_only,
        coefficient_maps_materialized = false,
        transform_contracts_materialized = false,
        retained_spaces_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        diagnostics = (;
            source = :_cartesian_terminal_region_lowering_contract_inventory,
            private_development_only = true,
            terminal_region_metadata_only = true,
            lowering_contracts_metadata_only = true,
            all_units_have_lowering_contracts = true,
            final_retained_unit_inventory_available = false,
            pair_inventory_available = false,
            coefficient_maps_materialized = false,
            transform_contracts_materialized = false,
            retained_spaces_materialized = false,
            operator_blocks_materialized = false,
            pair_operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
            artifacts_materialized = false,
            materialization_behavior_changed = false,
            public_default_behavior_changed = false,
        ),
    )
end
