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
