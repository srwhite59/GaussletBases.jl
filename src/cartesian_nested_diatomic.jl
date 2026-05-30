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
struct _CartesianNestedBondAlignedDiatomicSource3D{B,S<:_AbstractCartesianNestedShellLayer3D}
    basis::B
    axis_bundles::_CartesianNestedAxisBundles3D
    nside::Int
    child_shell_retention_contract::CartesianNestedCompleteShellRetentionContract
    shared_shell_retention_contract::CartesianNestedCompleteShellRetentionContract
    geometry::_BondAlignedDiatomicSplitGeometry3D
    shared_shell_layers::Vector{S}
    child_sequences::Vector{_CartesianNestedShellSequence3D}
    child_column_ranges::Vector{UnitRange{Int}}
    midpoint_slab_column_range::Union{Nothing,UnitRange{Int}}
    sequence::_CartesianNestedShellSequence3D
end

"""
    _BondAlignedDiatomicAtomGrowthRecipe3D

Internal high-order/general-q anatomy recipe for bond-aligned diatomic atom
growth. This is metadata only: it records the official box policy without
building source functions, QW operators, Hamiltonian kernels, or quadrature
objects.
"""
struct _BondAlignedDiatomicAtomGrowthRecipe3D
    parent_box::NTuple{3,UnitRange{Int}}
    bond_axis::Symbol
    atom_axis_indices::NTuple{2,Int}
    protected_atom_side_count::Int
    cover_parent::Bool
    contact_cap_policy::Symbol
    mismatch_absorption_policy::Symbol
end

"""
    _BondAlignedDiatomicAtomGrowthCoverage3D

Support-coverage audit for the atom-growth anatomy layer. The shared molecular
support is the parent-box exterior outside the atom/contact bounding box.
"""
struct _BondAlignedDiatomicAtomGrowthCoverage3D
    expected_support_count::Int
    atom_contact_support_count::Int
    shared_molecular_support_count::Int
    covered_support_count::Int
    duplicate_count::Int
    missing_count::Int
    outside_count::Int
    status::Symbol
    coverage_ok::Bool
end

"""
    _BondAlignedDiatomicAtomGrowthAnatomy3D

Internal report for the protected-atom growth policy. Atom-local boxes grow in
side-count steps of two until they touch or leave exactly one shared contact
layer. Any parity/spacing/q mismatch between the atom/contact interior and the
full parent is assigned to the outermost shared molecular shell.
"""
struct _BondAlignedDiatomicAtomGrowthAnatomy3D
    recipe::_BondAlignedDiatomicAtomGrowthRecipe3D
    atom_side_count_ladder::Vector{Int}
    final_atom_side_count::Int
    contact_gap_count::Int
    contact_policy::Symbol
    left_atom_box::NTuple{3,UnitRange{Int}}
    right_atom_box::NTuple{3,UnitRange{Int}}
    contact_box::Union{Nothing,NTuple{3,UnitRange{Int}}}
    inner_atom_contact_box::NTuple{3,UnitRange{Int}}
    outer_regular_start_box::NTuple{3,UnitRange{Int}}
    regular_shared_shell_count::Int
    outer_mismatch_low_counts::NTuple{3,Int}
    outer_mismatch_high_counts::NTuple{3,Int}
    support_coverage::_BondAlignedDiatomicAtomGrowthCoverage3D
end

"""
    _BondAlignedDiatomicAtomGrowthConstructionRegion3D

One ordered, disjoint support region in the recipe-backed atom-growth
construction plan. Regions are geometry/provenance only: they do not carry
coefficient maps, operators, backend state, or Hamiltonian packets.
"""
struct _BondAlignedDiatomicAtomGrowthConstructionRegion3D
    role::Symbol
    order_index::Int
    box::NTuple{3,UnitRange{Int}}
    inner_exclusion_box::Union{Nothing,NTuple{3,UnitRange{Int}}}
    support_indices::Vector{Int}
    metadata::NamedTuple
end

"""
    _BondAlignedDiatomicAtomGrowthConstructionCoverage3D

Support audit for atom-growth construction-plan regions. Duplicate counts
represent extra ownership occurrences beyond the first region claiming a
parent support index.
"""
struct _BondAlignedDiatomicAtomGrowthConstructionCoverage3D
    expected_support_count::Int
    region_support_count::Int
    covered_support_count::Int
    duplicate_count::Int
    missing_count::Int
    outside_count::Int
    coverage_ok::Bool
end

"""
    _BondAlignedDiatomicAtomGrowthConstructionPlan3D

Internal construction-region plan derived from
`_BondAlignedDiatomicAtomGrowthAnatomy3D`. The ordered regions are:

- optional outermost mismatch shell
- regular shared molecular shells, ordered outside-in
- left and right protected atom boxes
- optional one-layer contact cap

This is still a diagnostic/planning object and is not used by the active
diatomic source builder.
"""
struct _BondAlignedDiatomicAtomGrowthConstructionPlan3D
    anatomy::_BondAlignedDiatomicAtomGrowthAnatomy3D
    regions::Vector{_BondAlignedDiatomicAtomGrowthConstructionRegion3D}
    region_order::Vector{Symbol}
    support_coverage::_BondAlignedDiatomicAtomGrowthConstructionCoverage3D
    outer_mismatch_is_outermost::Bool
    middle_contact_clean::Bool
end

"""
    _BondAlignedDiatomicHighOrderRecipeRegionChoice3D

Internal per-region high-order recipe selection metadata. This records what a
future construction consumer selected by q/order and whether that piece is
currently buildable or only planned/experimental. It does not build any
coefficient maps or operators.
"""
struct _BondAlignedDiatomicHighOrderRecipeRegionChoice3D
    region_role::Symbol
    region_order_index::Int
    region_category::Symbol
    recipe_family::Symbol
    q::Int
    order::Int
    q_min::Int
    support_count::Int
    implementation_status::Symbol
    buildability_status::Symbol
    selected::Bool
    metadata::NamedTuple
end

"""
    _BondAlignedDiatomicHighOrderRecipePolicy3D

Internal q-organized high-order recipe policy layered on an atom-growth
construction plan. This is a diagnostic/planning object only; active diatomic
source construction and Hamiltonian paths do not consume it yet.
"""
struct _BondAlignedDiatomicHighOrderRecipePolicy3D
    construction_plan::_BondAlignedDiatomicAtomGrowthConstructionPlan3D
    recipe_label::Symbol
    q_min::Int
    region_choices::Vector{_BondAlignedDiatomicHighOrderRecipeRegionChoice3D}
    q_region_counts::Dict{Int,Int}
    buildable_region_count::Int
    planned_region_count::Int
    experimental_region_count::Int
    metadata::NamedTuple
end

"""
    _BondAlignedDiatomicHighOrderConstructionPiece3D

One disjoint, support-owning construction piece inside a selected high-order
recipe region. Pieces are realization metadata only: they describe how a future
opt-in builder can feed owned-unit/shell-like primitives without constructing
coefficients, operators, backend state, or Hamiltonian packets.
"""
struct _BondAlignedDiatomicHighOrderConstructionPiece3D
    role::Symbol
    piece_index::Int
    primitive_family::Symbol
    box::NTuple{3,UnitRange{Int}}
    inner_exclusion_box::Union{Nothing,NTuple{3,UnitRange{Int}}}
    support_indices::Vector{Int}
    metadata::NamedTuple
end

"""
    _BondAlignedDiatomicHighOrderConstructionPieceCoverage3D

Support audit for the pieces associated with one recipe region.
"""
struct _BondAlignedDiatomicHighOrderConstructionPieceCoverage3D
    expected_support_count::Int
    piece_support_count::Int
    covered_support_count::Int
    duplicate_count::Int
    missing_count::Int
    outside_count::Int
    coverage_ok::Bool
end

"""
    _BondAlignedDiatomicHighOrderRecipeRegionRealizationDescriptor3D

Construction-piece descriptor for one selected high-order recipe region. This
closes region-level planning gaps without making the active source builder
consume the policy.
"""
struct _BondAlignedDiatomicHighOrderRecipeRegionRealizationDescriptor3D
    region_role::Symbol
    region_order_index::Int
    recipe_family::Symbol
    q::Int
    order::Int
    parent_support_count::Int
    owned_unit_count::Int
    primitive_family::Symbol
    pieces::Vector{_BondAlignedDiatomicHighOrderConstructionPiece3D}
    support_coverage::_BondAlignedDiatomicHighOrderConstructionPieceCoverage3D
    exact_full_coverage::Bool
    metadata::NamedTuple
end

"""
    _BondAlignedDiatomicHighOrderRecipeRegionRealization3D

Internal build-readiness row for one selected high-order recipe region. The
row names the existing primitive a future opt-in builder would use, or the
missing implementation needed before the selected region can be consumed.
"""
struct _BondAlignedDiatomicHighOrderRecipeRegionRealization3D
    region_role::Symbol
    region_order_index::Int
    region_category::Symbol
    recipe_family::Symbol
    q::Int
    order::Int
    support_count::Int
    buildability_status::Symbol
    mapped_primitive::Union{Nothing,Symbol}
    mapped_primitive_status::Symbol
    missing_implementation::Union{Nothing,Symbol}
    active_builder_consumes::Bool
    existing_opt_in_route::Union{Nothing,Symbol}
    realization_descriptor::Union{
        Nothing,
        _BondAlignedDiatomicHighOrderRecipeRegionRealizationDescriptor3D,
    }
    metadata::NamedTuple
end

"""
    _BondAlignedDiatomicHighOrderRecipeRealizationAudit3D

Internal audit for mapping a q-organized recipe policy to currently available
construction primitives. This is metadata only; it does not call those
primitives or alter source construction.
"""
struct _BondAlignedDiatomicHighOrderRecipeRealizationAudit3D
    policy::_BondAlignedDiatomicHighOrderRecipePolicy3D
    region_realizations::Vector{_BondAlignedDiatomicHighOrderRecipeRegionRealization3D}
    mapped_region_count::Int
    missing_region_count::Int
    buildable_without_mapped_primitive_count::Int
    active_builder_consumed_region_count::Int
    ready_for_opt_in_builder::Bool
    metadata::NamedTuple
end

"""
    _BondAlignedDiatomicHighOrderRecipeRegionSourceBuild3D

One region consumed by the experimental recipe-backed source-piece builder.
This records which existing primitive was called and what support/columns it
produced. It is opt-in construction provenance only; the active diatomic source
builder does not use these rows.
"""
struct _BondAlignedDiatomicHighOrderRecipeRegionSourceBuild3D
    region_role::Symbol
    region_order_index::Int
    region_category::Symbol
    recipe_family::Symbol
    q::Int
    order::Int
    mapped_primitive::Symbol
    primitive_family::Symbol
    region_support_count::Int
    built_support_count::Int
    retained_count::Int
    column_range::UnitRange{Int}
    support_coverage::_BondAlignedDiatomicHighOrderConstructionPieceCoverage3D
    built::Bool
    unsupported_reason::Union{Nothing,Symbol}
    active_builder_consumes::Bool
    built_object::Any
    metadata::NamedTuple
end

"""
    _BondAlignedDiatomicHighOrderRecipeSourceConstruction3D

Experimental opt-in source-piece construction for a ready q4 atom-growth
recipe. It delegates to existing atom-box, endcap/panel, and direct-slab
primitives, then merges the pieces into a nested shell sequence. It does not
change public/default source construction or any Hamiltonian route.
"""
struct _BondAlignedDiatomicHighOrderRecipeSourceConstruction3D
    basis::Any
    axis_bundles::_CartesianNestedAxisBundles3D
    policy::_BondAlignedDiatomicHighOrderRecipePolicy3D
    realization_audit::_BondAlignedDiatomicHighOrderRecipeRealizationAudit3D
    region_builds::Vector{_BondAlignedDiatomicHighOrderRecipeRegionSourceBuild3D}
    sequence::_CartesianNestedShellSequence3D
    shared_shell_layers::Vector{_AbstractCartesianNestedShellLayer3D}
    core_support_indices::Vector{Int}
    core_coefficient_matrix::_CartesianCoefficientMap
    support_coverage::_BondAlignedDiatomicHighOrderConstructionPieceCoverage3D
    consumed_region_count::Int
    unsupported_region_count::Int
    active_builder_consumes::Bool
    metadata::NamedTuple
end

"""
    _BondAlignedDiatomicHighOrderRecipeSourceReadiness3D

Readiness audit for handing an opt-in recipe source construction to existing
nested source/fixed-block consumers. Fixed-block handoff is allowed only when
the construction already carries an assembled shell-sequence packet. Legacy
diatomic source-object handoff remains false until the recipe path owns the old
split-geometry/child-sequence fields explicitly.
"""
struct _BondAlignedDiatomicHighOrderRecipeSourceReadiness3D
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D
    parent_dimension::Int
    fixed_dimension::Int
    support_coverage::_BondAlignedDiatomicHighOrderConstructionPieceCoverage3D
    sequence_packet_available::Bool
    overlap_available::Bool
    weights_available::Bool
    overlap_error::Union{Nothing,Float64}
    can_produce_nested_source::Bool
    nested_source_missing_fields::Vector{Symbol}
    can_produce_fixed_block::Bool
    fixed_block_missing_fields::Vector{Symbol}
    fixed_block::Union{Nothing,_NestedFixedBlock3D}
    active_builder_consumes::Bool
    default_builders_unchanged::Bool
    metadata::NamedTuple
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

function Base.show(io::IO, recipe::_BondAlignedDiatomicAtomGrowthRecipe3D)
    print(
        io,
        "_BondAlignedDiatomicAtomGrowthRecipe3D(axis=:",
        recipe.bond_axis,
        ", atoms=",
        recipe.atom_axis_indices,
        ", core_side=",
        recipe.protected_atom_side_count,
        ", cover_parent=",
        recipe.cover_parent,
        ")",
    )
end

function Base.show(io::IO, anatomy::_BondAlignedDiatomicAtomGrowthAnatomy3D)
    print(
        io,
        "_BondAlignedDiatomicAtomGrowthAnatomy3D(axis=:",
        anatomy.recipe.bond_axis,
        ", final_atom_side=",
        anatomy.final_atom_side_count,
        ", contact_policy=:",
        anatomy.contact_policy,
        ", regular_shared_shells=",
        anatomy.regular_shared_shell_count,
        ", coverage=:",
        anatomy.support_coverage.status,
        ")",
    )
end

function Base.show(io::IO, plan::_BondAlignedDiatomicAtomGrowthConstructionPlan3D)
    print(
        io,
        "_BondAlignedDiatomicAtomGrowthConstructionPlan3D(nregions=",
        length(plan.regions),
        ", coverage_ok=",
        plan.support_coverage.coverage_ok,
        ", outer_mismatch_is_outermost=",
        plan.outer_mismatch_is_outermost,
        ")",
    )
end

function Base.show(io::IO, policy::_BondAlignedDiatomicHighOrderRecipePolicy3D)
    print(
        io,
        "_BondAlignedDiatomicHighOrderRecipePolicy3D(label=:",
        policy.recipe_label,
        ", q_min=",
        policy.q_min,
        ", nregions=",
        length(policy.region_choices),
        ", buildable=",
        policy.buildable_region_count,
        ")",
    )
end

function Base.show(io::IO, audit::_BondAlignedDiatomicHighOrderRecipeRealizationAudit3D)
    print(
        io,
        "_BondAlignedDiatomicHighOrderRecipeRealizationAudit3D(label=:",
        audit.policy.recipe_label,
        ", mapped=",
        audit.mapped_region_count,
        ", missing=",
        audit.missing_region_count,
        ", ready=",
        audit.ready_for_opt_in_builder,
        ")",
    )
end

function Base.show(
    io::IO,
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
)
    print(
        io,
        "_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D(label=:",
        construction.policy.recipe_label,
        ", consumed=",
        construction.consumed_region_count,
        ", unsupported=",
        construction.unsupported_region_count,
        ", nfixed=",
        size(construction.sequence.coefficient_matrix, 2),
        ")",
    )
end

function Base.show(
    io::IO,
    readiness::_BondAlignedDiatomicHighOrderRecipeSourceReadiness3D,
)
    print(
        io,
        "_BondAlignedDiatomicHighOrderRecipeSourceReadiness3D(fixed_ready=",
        readiness.can_produce_fixed_block,
        ", nested_source_ready=",
        readiness.can_produce_nested_source,
        ", nfixed=",
        readiness.fixed_dimension,
        ")",
    )
end

function _nested_source_contract_audit(source::_CartesianNestedBondAlignedDiatomicSource3D)
    return _nested_shell_sequence_contract_audit(source.sequence, _nested_axis_lengths(source.axis_bundles))
end

function _nested_validate_diatomic_atom_growth_parent_box(parent_box::NTuple{3,UnitRange{Int}})
    all(length(interval) >= 1 for interval in parent_box) || throw(
        ArgumentError("diatomic atom-growth recipe requires nonempty parent-box intervals"),
    )
    return parent_box
end

function _nested_bond_aligned_diatomic_atom_growth_recipe(
    parent_box::NTuple{3,UnitRange{Int}};
    bond_axis::Symbol = :z,
    atom_axis_indices::NTuple{2,<:Integer},
    protected_atom_side_count::Integer,
    cover_parent::Bool = true,
)
    _nested_validate_diatomic_atom_growth_parent_box(parent_box)
    axis = _nested_axis_index(bond_axis)
    side_count = Int(protected_atom_side_count)
    side_count >= 1 || throw(
        ArgumentError("diatomic atom-growth protected atom side count must be positive"),
    )
    atoms = (Int(atom_axis_indices[1]), Int(atom_axis_indices[2]))
    atoms[1] < atoms[2] || throw(
        ArgumentError("diatomic atom-growth atom_axis_indices must be ordered left to right"),
    )
    interval = parent_box[axis]
    first(interval) <= atoms[1] <= last(interval) || throw(
        ArgumentError("diatomic atom-growth left atom index must lie inside the parent bond-axis interval"),
    )
    first(interval) <= atoms[2] <= last(interval) || throw(
        ArgumentError("diatomic atom-growth right atom index must lie inside the parent bond-axis interval"),
    )
    return _BondAlignedDiatomicAtomGrowthRecipe3D(
        parent_box,
        bond_axis,
        atoms,
        side_count,
        cover_parent,
        :single_shared_contact_cap,
        :outermost_shared_molecular_shell,
    )
end

function _nested_diatomic_atom_growth_nearest_axis_index(
    centers_axis::AbstractVector{<:Real},
    coordinate::Real,
)
    _, index = findmin(abs.(Float64.(centers_axis) .- Float64(coordinate)))
    return Int(index)
end

function _nested_diatomic_atom_growth_bond_axis(basis, bond_axis::Union{Nothing,Symbol})
    if isnothing(bond_axis)
        hasproperty(basis, :bond_axis) || throw(
            ArgumentError("diatomic atom-growth basis recipe requires bond_axis or a basis carrying bond_axis"),
        )
        return getproperty(basis, :bond_axis)
    end
    return bond_axis
end

function _nested_diatomic_atom_growth_atom_axis_indices(
    basis,
    bundles::_CartesianNestedAxisBundles3D,
    bond_axis::Symbol,
)
    length(basis.nuclei) == 2 || throw(
        ArgumentError("diatomic atom-growth basis recipe requires exactly two nuclei"),
    )
    axis = _nested_axis_index(bond_axis)
    centers_axis = _nested_axis_pgdg(bundles, bond_axis).centers
    indices = sort(Int[
        _nested_diatomic_atom_growth_nearest_axis_index(centers_axis, nucleus[axis])
        for nucleus in basis.nuclei
    ])
    indices[1] < indices[2] || throw(
        ArgumentError("diatomic atom-growth basis recipe requires nuclei to map to distinct bond-axis parent sites"),
    )
    return (indices[1], indices[2])
end

function _nested_bond_aligned_diatomic_atom_growth_recipe(
    basis,
    bundles::_CartesianNestedAxisBundles3D;
    bond_axis::Union{Nothing,Symbol} = nothing,
    protected_atom_side_count::Integer,
    cover_parent::Bool = true,
)
    resolved_bond_axis = _nested_diatomic_atom_growth_bond_axis(basis, bond_axis)
    dims = _nested_axis_lengths(bundles)
    parent_box = (1:dims[1], 1:dims[2], 1:dims[3])
    return _nested_bond_aligned_diatomic_atom_growth_recipe(
        parent_box;
        bond_axis = resolved_bond_axis,
        atom_axis_indices = _nested_diatomic_atom_growth_atom_axis_indices(
            basis,
            bundles,
            resolved_bond_axis,
        ),
        protected_atom_side_count,
        cover_parent,
    )
end

function _nested_diatomic_atom_growth_interval(
    center::Int,
    side_count::Int,
    even_bias::Symbol,
)
    if isodd(side_count)
        half = div(side_count, 2)
        return (center - half):(center + half)
    end
    half = div(side_count, 2)
    even_bias == :lower && return (center - half):(center + half - 1)
    even_bias == :upper && return (center - half + 1):(center + half)
    throw(ArgumentError("diatomic atom-growth even interval bias must be :lower or :upper"))
end

function _nested_diatomic_atom_growth_transverse_center(interval::UnitRange{Int})
    return first(interval) + div(length(interval) - 1, 2)
end

function _nested_diatomic_atom_growth_box(
    parent_box::NTuple{3,UnitRange{Int}},
    bond_axis::Symbol,
    atom_axis_index::Int,
    side_count::Int,
    atom_side::Symbol,
)
    axis = _nested_axis_index(bond_axis)
    atom_side in (:left, :right) || throw(
        ArgumentError("diatomic atom-growth box side must be :left or :right"),
    )
    return ntuple(3) do index
        if index == axis
            bias = atom_side == :left ? :upper : :lower
            _nested_diatomic_atom_growth_interval(atom_axis_index, side_count, bias)
        else
            center = _nested_diatomic_atom_growth_transverse_center(parent_box[index])
            _nested_diatomic_atom_growth_interval(center, side_count, :upper)
        end
    end
end

function _nested_interval_inside_parent(
    interval::UnitRange{Int},
    parent_interval::UnitRange{Int},
)
    return first(parent_interval) <= first(interval) && last(interval) <= last(parent_interval)
end

function _nested_box_inside_parent(
    box::NTuple{3,UnitRange{Int}},
    parent_box::NTuple{3,UnitRange{Int}},
)
    return all(
        _nested_interval_inside_parent(box[index], parent_box[index]) for index in 1:3
    )
end

function _nested_diatomic_atom_growth_contact_box(
    left_atom_box::NTuple{3,UnitRange{Int}},
    right_atom_box::NTuple{3,UnitRange{Int}},
    bond_axis::Symbol,
)
    axis = _nested_axis_index(bond_axis)
    contact_axis = (last(left_atom_box[axis]) + 1):(first(right_atom_box[axis]) - 1)
    length(contact_axis) == 1 || throw(
        ArgumentError("diatomic atom-growth contact cap requires exactly one missing bond-axis layer"),
    )
    return ntuple(index -> index == axis ? contact_axis : left_atom_box[index], 3)
end

function _nested_diatomic_bounding_box(
    boxes::AbstractVector{<:NTuple{3,UnitRange{Int}}},
)
    isempty(boxes) && throw(ArgumentError("diatomic atom-growth bounding box requires at least one box"))
    return ntuple(3) do index
        minimum(first(box[index]) for box in boxes):maximum(last(box[index]) for box in boxes)
    end
end

function _nested_diatomic_box_support_count(box::NTuple{3,UnitRange{Int}})
    return prod(length.(box))
end

function _nested_diatomic_atom_growth_component_count(
    component_boxes::AbstractVector{<:NTuple{3,UnitRange{Int}}},
    parent_box::NTuple{3,UnitRange{Int}},
)
    seen = Set{NTuple{3,Int}}()
    duplicate_count = 0
    outside_count = 0
    for box in component_boxes, ix in box[1], iy in box[2], iz in box[3]
        state = (ix, iy, iz)
        if !(
            first(parent_box[1]) <= ix <= last(parent_box[1]) &&
            first(parent_box[2]) <= iy <= last(parent_box[2]) &&
            first(parent_box[3]) <= iz <= last(parent_box[3])
        )
            outside_count += 1
        end
        if state in seen
            duplicate_count += 1
        else
            push!(seen, state)
        end
    end
    return (
        unique_support_count = length(seen),
        duplicate_count = duplicate_count,
        outside_count = outside_count,
    )
end

function _nested_diatomic_atom_growth_outer_regular_start(
    parent_box::NTuple{3,UnitRange{Int}},
    inner_box::NTuple{3,UnitRange{Int}},
)
    low_margins = ntuple(index -> first(inner_box[index]) - first(parent_box[index]), 3)
    high_margins = ntuple(index -> last(parent_box[index]) - last(inner_box[index]), 3)
    all(value >= 0 for value in low_margins) && all(value >= 0 for value in high_margins) || throw(
        ArgumentError("diatomic atom-growth inner atom/contact box must lie inside the parent box"),
    )
    regular_shell_count = minimum((low_margins..., high_margins...))
    outer_regular_start_box = ntuple(
        index -> (first(inner_box[index]) - regular_shell_count):(last(inner_box[index]) + regular_shell_count),
        3,
    )
    outer_mismatch_low_counts = ntuple(index -> first(outer_regular_start_box[index]) - first(parent_box[index]), 3)
    outer_mismatch_high_counts = ntuple(index -> last(parent_box[index]) - last(outer_regular_start_box[index]), 3)
    return (
        outer_regular_start_box = outer_regular_start_box,
        regular_shared_shell_count = regular_shell_count,
        outer_mismatch_low_counts = outer_mismatch_low_counts,
        outer_mismatch_high_counts = outer_mismatch_high_counts,
    )
end

function _nested_diatomic_atom_growth_support_coverage(
    recipe::_BondAlignedDiatomicAtomGrowthRecipe3D,
    component_boxes::AbstractVector{<:NTuple{3,UnitRange{Int}}},
    inner_box::NTuple{3,UnitRange{Int}},
)
    parent_count = _nested_diatomic_box_support_count(recipe.parent_box)
    component_count = _nested_diatomic_atom_growth_component_count(component_boxes, recipe.parent_box)
    atom_contact_count = component_count.unique_support_count
    shared_count = recipe.cover_parent ? parent_count - _nested_diatomic_box_support_count(inner_box) : 0
    covered_count = atom_contact_count + shared_count
    missing_count = recipe.cover_parent ? 0 : max(parent_count - covered_count, 0)
    status =
        recipe.cover_parent && component_count.duplicate_count == 0 && component_count.outside_count == 0 && missing_count == 0 ?
        :full_parent_covered :
        recipe.cover_parent ? :invalid_full_parent_coverage : :cropped_parent
    return _BondAlignedDiatomicAtomGrowthCoverage3D(
        parent_count,
        atom_contact_count,
        shared_count,
        covered_count,
        component_count.duplicate_count,
        missing_count,
        component_count.outside_count,
        status,
        status == :full_parent_covered,
    )
end

function _nested_bond_aligned_diatomic_atom_growth_anatomy(
    recipe::_BondAlignedDiatomicAtomGrowthRecipe3D,
)
    axis = _nested_axis_index(recipe.bond_axis)
    side_count = recipe.protected_atom_side_count
    side_ladder = Int[]
    while true
        push!(side_ladder, side_count)
        left_atom_box = _nested_diatomic_atom_growth_box(
            recipe.parent_box,
            recipe.bond_axis,
            recipe.atom_axis_indices[1],
            side_count,
            :left,
        )
        right_atom_box = _nested_diatomic_atom_growth_box(
            recipe.parent_box,
            recipe.bond_axis,
            recipe.atom_axis_indices[2],
            side_count,
            :right,
        )
        if !(
            _nested_box_inside_parent(left_atom_box, recipe.parent_box) &&
            _nested_box_inside_parent(right_atom_box, recipe.parent_box)
        )
            throw(
                ArgumentError(
                    "diatomic atom-growth boxes cannot reach touch or one-layer contact before hitting the parent boundary",
                ),
            )
        end
        all(left_atom_box[index] == right_atom_box[index] for index in 1:3 if index != axis) || throw(
            ArgumentError("diatomic atom-growth atom boxes must share identical transverse ranges"),
        )
        contact_gap_count = first(right_atom_box[axis]) - last(left_atom_box[axis]) - 1
        contact_gap_count >= 0 || throw(
            ArgumentError("diatomic atom-growth protected atom boxes overlap before clean contact"),
        )
        if contact_gap_count <= 1
            contact_policy =
                contact_gap_count == 0 ? :touching_atom_boxes : recipe.contact_cap_policy
            contact_box =
                contact_gap_count == 1 ?
                _nested_diatomic_atom_growth_contact_box(left_atom_box, right_atom_box, recipe.bond_axis) :
                nothing
            component_boxes = isnothing(contact_box) ?
                [left_atom_box, right_atom_box] :
                [left_atom_box, right_atom_box, contact_box]
            inner_atom_contact_box = _nested_diatomic_bounding_box(component_boxes)
            outer = _nested_diatomic_atom_growth_outer_regular_start(
                recipe.parent_box,
                inner_atom_contact_box,
            )
            support_coverage = _nested_diatomic_atom_growth_support_coverage(
                recipe,
                component_boxes,
                inner_atom_contact_box,
            )
            return _BondAlignedDiatomicAtomGrowthAnatomy3D(
                recipe,
                side_ladder,
                side_count,
                contact_gap_count,
                contact_policy,
                left_atom_box,
                right_atom_box,
                contact_box,
                inner_atom_contact_box,
                outer.outer_regular_start_box,
                outer.regular_shared_shell_count,
                outer.outer_mismatch_low_counts,
                outer.outer_mismatch_high_counts,
                support_coverage,
            )
        end
        side_count += 2
    end
end

function _nested_bond_aligned_diatomic_atom_growth_anatomy(
    parent_box::NTuple{3,UnitRange{Int}};
    kwargs...,
)
    recipe = _nested_bond_aligned_diatomic_atom_growth_recipe(parent_box; kwargs...)
    return _nested_bond_aligned_diatomic_atom_growth_anatomy(recipe)
end

function _nested_bond_aligned_diatomic_atom_growth_anatomy(
    basis,
    bundles::_CartesianNestedAxisBundles3D;
    kwargs...,
)
    recipe = _nested_bond_aligned_diatomic_atom_growth_recipe(basis, bundles; kwargs...)
    return _nested_bond_aligned_diatomic_atom_growth_anatomy(recipe)
end

function _nested_diatomic_atom_growth_plan_dims(
    parent_box::NTuple{3,UnitRange{Int}},
)
    all(first(parent_box[index]) == 1 for index in 1:3) || throw(
        ArgumentError("diatomic atom-growth construction plan requires a 1-based parent index domain"),
    )
    return ntuple(index -> length(parent_box[index]), 3)
end

function _nested_diatomic_state_inside_box(
    state::NTuple{3,Int},
    box::NTuple{3,UnitRange{Int}},
)
    return state[1] in box[1] && state[2] in box[2] && state[3] in box[3]
end

function _nested_diatomic_shell_support_indices(
    outer_box::NTuple{3,UnitRange{Int}},
    inner_box::NTuple{3,UnitRange{Int}},
    dims::NTuple{3,Int},
)
    _nested_box_inside_parent(inner_box, outer_box) || throw(
        ArgumentError("diatomic atom-growth shell region requires inner_box inside outer_box"),
    )
    support = Int[]
    for ix in outer_box[1], iy in outer_box[2], iz in outer_box[3]
        state = (ix, iy, iz)
        _nested_diatomic_state_inside_box(state, inner_box) && continue
        push!(support, _cartesian_flat_index(ix, iy, iz, dims))
    end
    return support
end

function _nested_diatomic_expand_box(
    box::NTuple{3,UnitRange{Int}},
    amount::Int,
)
    amount >= 0 || throw(ArgumentError("diatomic atom-growth box expansion must be nonnegative"))
    return ntuple(
        index -> (first(box[index]) - amount):(last(box[index]) + amount),
        3,
    )
end

function _nested_diatomic_atom_growth_plan_region(
    role::Symbol,
    order_index::Int,
    box::NTuple{3,UnitRange{Int}},
    inner_exclusion_box::Union{Nothing,NTuple{3,UnitRange{Int}}},
    dims::NTuple{3,Int};
    metadata::NamedTuple = (;),
)
    support_indices = isnothing(inner_exclusion_box) ?
        _nested_box_support_indices(box..., dims) :
        _nested_diatomic_shell_support_indices(box, inner_exclusion_box, dims)
    return _BondAlignedDiatomicAtomGrowthConstructionRegion3D(
        role,
        order_index,
        box,
        inner_exclusion_box,
        support_indices,
        metadata,
    )
end

function _nested_diatomic_atom_growth_construction_coverage(
    parent_box::NTuple{3,UnitRange{Int}},
    regions::AbstractVector{_BondAlignedDiatomicAtomGrowthConstructionRegion3D},
)
    dims = _nested_diatomic_atom_growth_plan_dims(parent_box)
    expected = Set(_nested_box_support_indices(parent_box..., dims))
    owned_counts = Dict{Int,Int}()
    region_support_count = 0
    for region in regions
        region_support_count += length(region.support_indices)
        for index in region.support_indices
            owned_counts[index] = get(owned_counts, index, 0) + 1
        end
    end
    owned = Set(keys(owned_counts))
    duplicate_count = sum(max(count - 1, 0) for count in values(owned_counts))
    missing_count = length(setdiff(expected, owned))
    outside_count = length(setdiff(owned, expected))
    coverage_ok = duplicate_count == 0 && missing_count == 0 && outside_count == 0
    return _BondAlignedDiatomicAtomGrowthConstructionCoverage3D(
        length(expected),
        region_support_count,
        length(owned),
        duplicate_count,
        missing_count,
        outside_count,
        coverage_ok,
    )
end

function _nested_diatomic_atom_growth_outer_mismatch_is_outermost(
    regions::AbstractVector{_BondAlignedDiatomicAtomGrowthConstructionRegion3D},
)
    mismatch_indices = findall(
        region -> region.role == :outer_mismatch_shared_molecular_shell,
        regions,
    )
    isempty(mismatch_indices) && return true
    return mismatch_indices == [1]
end

function _nested_diatomic_atom_growth_middle_contact_clean(
    anatomy::_BondAlignedDiatomicAtomGrowthAnatomy3D,
    regions::AbstractVector{_BondAlignedDiatomicAtomGrowthConstructionRegion3D},
)
    atom_contact_roles = Set([:left_atom_box, :right_atom_box, :contact_cap])
    atom_contact_count = sum(
        length(region.support_indices) for region in regions if region.role in atom_contact_roles
    )
    contact_policy_ok =
        anatomy.contact_policy == :touching_atom_boxes ||
        any(region -> region.role == :contact_cap, regions)
    return contact_policy_ok && atom_contact_count == anatomy.support_coverage.atom_contact_support_count
end

function _nested_bond_aligned_diatomic_atom_growth_construction_plan(
    anatomy::_BondAlignedDiatomicAtomGrowthAnatomy3D,
)
    anatomy.support_coverage.coverage_ok || throw(
        ArgumentError("diatomic atom-growth construction plan requires clean anatomy support coverage"),
    )
    dims = _nested_diatomic_atom_growth_plan_dims(anatomy.recipe.parent_box)
    regions = _BondAlignedDiatomicAtomGrowthConstructionRegion3D[]
    order_index = 1

    outer_mismatch = _nested_diatomic_atom_growth_plan_region(
        :outer_mismatch_shared_molecular_shell,
        order_index,
        anatomy.recipe.parent_box,
        anatomy.outer_regular_start_box,
        dims;
        metadata = (
            mismatch_absorption_policy = anatomy.recipe.mismatch_absorption_policy,
            low_counts = anatomy.outer_mismatch_low_counts,
            high_counts = anatomy.outer_mismatch_high_counts,
            is_outermost = true,
        ),
    )
    if !isempty(outer_mismatch.support_indices)
        push!(regions, outer_mismatch)
        order_index += 1
    end

    for shell_offset in anatomy.regular_shared_shell_count:-1:1
        current_box = _nested_diatomic_expand_box(anatomy.inner_atom_contact_box, shell_offset)
        next_inner_box = _nested_diatomic_expand_box(anatomy.inner_atom_contact_box, shell_offset - 1)
        push!(
            regions,
            _nested_diatomic_atom_growth_plan_region(
                :regular_shared_molecular_shell,
                order_index,
                current_box,
                next_inner_box,
                dims;
                metadata = (
                    shell_offset = shell_offset,
                    outside_in_index = anatomy.regular_shared_shell_count - shell_offset + 1,
                    total_regular_shared_shells = anatomy.regular_shared_shell_count,
                ),
            ),
        )
        order_index += 1
    end

    push!(
        regions,
        _nested_diatomic_atom_growth_plan_region(
            :left_atom_box,
            order_index,
            anatomy.left_atom_box,
            nothing,
            dims;
            metadata = (
                atom_side = :left,
                atom_axis_index = anatomy.recipe.atom_axis_indices[1],
                final_atom_side_count = anatomy.final_atom_side_count,
            ),
        ),
    )
    order_index += 1
    push!(
        regions,
        _nested_diatomic_atom_growth_plan_region(
            :right_atom_box,
            order_index,
            anatomy.right_atom_box,
            nothing,
            dims;
            metadata = (
                atom_side = :right,
                atom_axis_index = anatomy.recipe.atom_axis_indices[2],
                final_atom_side_count = anatomy.final_atom_side_count,
            ),
        ),
    )
    order_index += 1
    if !isnothing(anatomy.contact_box)
        push!(
            regions,
            _nested_diatomic_atom_growth_plan_region(
                :contact_cap,
                order_index,
                anatomy.contact_box,
                nothing,
                dims;
                metadata = (
                    contact_policy = anatomy.contact_policy,
                    contact_gap_count = anatomy.contact_gap_count,
                ),
            ),
        )
    end

    support_coverage = _nested_diatomic_atom_growth_construction_coverage(
        anatomy.recipe.parent_box,
        regions,
    )
    return _BondAlignedDiatomicAtomGrowthConstructionPlan3D(
        anatomy,
        regions,
        [region.role for region in regions],
        support_coverage,
        _nested_diatomic_atom_growth_outer_mismatch_is_outermost(regions),
        _nested_diatomic_atom_growth_middle_contact_clean(anatomy, regions),
    )
end

function _nested_bond_aligned_diatomic_atom_growth_construction_plan(
    recipe::_BondAlignedDiatomicAtomGrowthRecipe3D,
)
    return _nested_bond_aligned_diatomic_atom_growth_construction_plan(
        _nested_bond_aligned_diatomic_atom_growth_anatomy(recipe),
    )
end

function _nested_bond_aligned_diatomic_atom_growth_construction_plan(
    parent_box::NTuple{3,UnitRange{Int}};
    kwargs...,
)
    return _nested_bond_aligned_diatomic_atom_growth_construction_plan(
        _nested_bond_aligned_diatomic_atom_growth_anatomy(parent_box; kwargs...),
    )
end

function _nested_diatomic_high_order_recipe_region_category(role::Symbol)
    role in (:left_atom_box, :right_atom_box) && return :atom_local
    role == :contact_cap && return :contact_cap
    role == :regular_shared_molecular_shell && return :shared_exterior
    role == :outer_mismatch_shared_molecular_shell && return :outer_mismatch
    throw(ArgumentError("unsupported diatomic high-order recipe region role $role"))
end

function _nested_diatomic_high_order_default_recipe_label(shared_exterior_family::Symbol)
    shared_exterior_family == :shared_endcap_panel_exterior &&
        return :mixed_atom_cubic_shared_endcap_panel
    shared_exterior_family == :transverse_annulus_exterior &&
        return :mixed_atom_cubic_shared_transverse_annulus
    throw(ArgumentError("unsupported shared-exterior recipe family $shared_exterior_family"))
end

function _nested_diatomic_high_order_recipe_status(recipe_family::Symbol)
    if recipe_family == :protected_atom_cubic_shell
        return (
            implementation_status = :validated_improvement,
            buildability_status = :buildable_now,
        )
    elseif recipe_family == :shared_endcap_panel_exterior
        return (
            implementation_status = :validated_improvement,
            buildability_status = :buildable_now,
        )
    elseif recipe_family == :transverse_annulus_exterior
        return (
            implementation_status = :promising_experimental,
            buildability_status = :planned_experimental,
        )
    elseif recipe_family == :shared_contact_cap
        return (
            implementation_status = :planned_policy_metadata,
            buildability_status = :planned_only,
        )
    elseif recipe_family == :outermost_mismatch_shared_molecular_shell
        return (
            implementation_status = :planned_policy_metadata,
            buildability_status = :planned_only,
        )
    end
    throw(ArgumentError("unsupported diatomic high-order recipe family $recipe_family"))
end

function _nested_diatomic_high_order_recipe_family(
    region_category::Symbol;
    atom_recipe_family::Symbol,
    shared_exterior_family::Symbol,
    contact_recipe_family::Symbol,
    outer_mismatch_recipe_family::Symbol,
)
    region_category == :atom_local && return atom_recipe_family
    region_category == :shared_exterior && return shared_exterior_family
    region_category == :contact_cap && return contact_recipe_family
    region_category == :outer_mismatch && return outer_mismatch_recipe_family
    throw(ArgumentError("unsupported diatomic high-order recipe region category $region_category"))
end

function _nested_diatomic_high_order_recipe_q_order(
    region_category::Symbol;
    q_min::Int,
    atom_q::Int,
    atom_order::Int,
    shared_q::Int,
    shared_order::Int,
    contact_q::Int,
    contact_order::Int,
    outer_mismatch_q::Int,
    outer_mismatch_order::Int,
)
    q, order =
        region_category == :atom_local ? (atom_q, atom_order) :
        region_category == :shared_exterior ? (shared_q, shared_order) :
        region_category == :contact_cap ? (contact_q, contact_order) :
        region_category == :outer_mismatch ? (outer_mismatch_q, outer_mismatch_order) :
        throw(ArgumentError("unsupported diatomic high-order recipe region category $region_category"))
    q >= q_min || throw(
        ArgumentError("diatomic high-order recipe region q=$q is below q_min=$q_min"),
    )
    order >= 1 || throw(
        ArgumentError("diatomic high-order recipe region order must be positive"),
    )
    return q, order
end

function _nested_diatomic_high_order_recipe_region_choice(
    region::_BondAlignedDiatomicAtomGrowthConstructionRegion3D;
    q_min::Int,
    atom_recipe_family::Symbol,
    shared_exterior_family::Symbol,
    contact_recipe_family::Symbol,
    outer_mismatch_recipe_family::Symbol,
    atom_q::Int,
    atom_order::Int,
    shared_q::Int,
    shared_order::Int,
    contact_q::Int,
    contact_order::Int,
    outer_mismatch_q::Int,
    outer_mismatch_order::Int,
)
    category = _nested_diatomic_high_order_recipe_region_category(region.role)
    family = _nested_diatomic_high_order_recipe_family(
        category;
        atom_recipe_family,
        shared_exterior_family,
        contact_recipe_family,
        outer_mismatch_recipe_family,
    )
    q, order = _nested_diatomic_high_order_recipe_q_order(
        category;
        q_min,
        atom_q,
        atom_order,
        shared_q,
        shared_order,
        contact_q,
        contact_order,
        outer_mismatch_q,
        outer_mismatch_order,
    )
    status = _nested_diatomic_high_order_recipe_status(family)
    return _BondAlignedDiatomicHighOrderRecipeRegionChoice3D(
        region.role,
        region.order_index,
        category,
        family,
        q,
        order,
        q_min,
        length(region.support_indices),
        status.implementation_status,
        status.buildability_status,
        true,
        (
            region_box = region.box,
            inner_exclusion_box = region.inner_exclusion_box,
            region_metadata = region.metadata,
        ),
    )
end

function _nested_diatomic_high_order_q_region_counts(
    choices::AbstractVector{_BondAlignedDiatomicHighOrderRecipeRegionChoice3D},
)
    counts = Dict{Int,Int}()
    for choice in choices
        counts[choice.q] = get(counts, choice.q, 0) + 1
    end
    return counts
end

function _nested_bond_aligned_diatomic_high_order_recipe_policy(
    plan::_BondAlignedDiatomicAtomGrowthConstructionPlan3D;
    recipe_label::Union{Nothing,Symbol} = nothing,
    q_min::Integer = 4,
    atom_q::Integer = q_min,
    atom_order::Integer = atom_q,
    shared_q::Integer = q_min,
    shared_order::Integer = shared_q,
    contact_q::Integer = shared_q,
    contact_order::Integer = contact_q,
    outer_mismatch_q::Integer = shared_q,
    outer_mismatch_order::Integer = outer_mismatch_q,
    atom_recipe_family::Symbol = :protected_atom_cubic_shell,
    shared_exterior_family::Symbol = :shared_endcap_panel_exterior,
    contact_recipe_family::Symbol = :shared_contact_cap,
    outer_mismatch_recipe_family::Symbol = :outermost_mismatch_shared_molecular_shell,
)
    q_min_int = Int(q_min)
    q_min_int >= 1 || throw(ArgumentError("diatomic high-order recipe q_min must be positive"))
    atom_q_int = Int(atom_q)
    atom_order_int = Int(atom_order)
    shared_q_int = Int(shared_q)
    shared_order_int = Int(shared_order)
    contact_q_int = Int(contact_q)
    contact_order_int = Int(contact_order)
    outer_mismatch_q_int = Int(outer_mismatch_q)
    outer_mismatch_order_int = Int(outer_mismatch_order)
    label = isnothing(recipe_label) ?
        _nested_diatomic_high_order_default_recipe_label(shared_exterior_family) :
        recipe_label

    choices = [
        _nested_diatomic_high_order_recipe_region_choice(
            region;
            q_min = q_min_int,
            atom_recipe_family,
            shared_exterior_family,
            contact_recipe_family,
            outer_mismatch_recipe_family,
            atom_q = atom_q_int,
            atom_order = atom_order_int,
            shared_q = shared_q_int,
            shared_order = shared_order_int,
            contact_q = contact_q_int,
            contact_order = contact_order_int,
            outer_mismatch_q = outer_mismatch_q_int,
            outer_mismatch_order = outer_mismatch_order_int,
        )
        for region in plan.regions
    ]
    buildable_count = count(
        choice -> choice.buildability_status == :buildable_now,
        choices,
    )
    planned_count = count(
        choice -> choice.buildability_status == :planned_only,
        choices,
    )
    experimental_count = count(
        choice -> choice.buildability_status == :planned_experimental,
        choices,
    )
    return _BondAlignedDiatomicHighOrderRecipePolicy3D(
        plan,
        label,
        q_min_int,
        choices,
        _nested_diatomic_high_order_q_region_counts(choices),
        buildable_count,
        planned_count,
        experimental_count,
        (
            active_builder_uses_policy = false,
            source_builder_changed = false,
            q_organized = true,
            support_coverage = plan.support_coverage,
            outer_mismatch_is_outermost = plan.outer_mismatch_is_outermost,
            middle_contact_clean = plan.middle_contact_clean,
        ),
    )
end

function _nested_bond_aligned_diatomic_high_order_recipe_policy(
    anatomy::_BondAlignedDiatomicAtomGrowthAnatomy3D;
    kwargs...,
)
    return _nested_bond_aligned_diatomic_high_order_recipe_policy(
        _nested_bond_aligned_diatomic_atom_growth_construction_plan(anatomy);
        kwargs...,
    )
end

function _nested_bond_aligned_diatomic_high_order_recipe_policy(
    parent_box::NTuple{3,UnitRange{Int}};
    bond_axis::Symbol = :z,
    atom_axis_indices::NTuple{2,<:Integer},
    protected_atom_side_count::Integer,
    cover_parent::Bool = true,
    kwargs...,
)
    plan = _nested_bond_aligned_diatomic_atom_growth_construction_plan(
        parent_box;
        bond_axis,
        atom_axis_indices,
        protected_atom_side_count,
        cover_parent,
    )
    return _nested_bond_aligned_diatomic_high_order_recipe_policy(
        plan;
        kwargs...,
    )
end

function _nested_bond_aligned_diatomic_high_order_recipe_policy(
    basis,
    bundles::_CartesianNestedAxisBundles3D;
    bond_axis::Union{Nothing,Symbol} = nothing,
    protected_atom_side_count::Integer,
    cover_parent::Bool = true,
    kwargs...,
)
    anatomy = _nested_bond_aligned_diatomic_atom_growth_anatomy(
        basis,
        bundles;
        bond_axis,
        protected_atom_side_count,
        cover_parent,
    )
    return _nested_bond_aligned_diatomic_high_order_recipe_policy(
        anatomy;
        kwargs...,
    )
end

function _nested_bond_aligned_diatomic_high_order_recipe_policy_diagnostics(
    policy::_BondAlignedDiatomicHighOrderRecipePolicy3D,
)
    return (
        recipe_label = policy.recipe_label,
        q_min = policy.q_min,
        q_region_counts = policy.q_region_counts,
        region_count = length(policy.region_choices),
        buildable_region_count = policy.buildable_region_count,
        planned_region_count = policy.planned_region_count,
        experimental_region_count = policy.experimental_region_count,
        region_choices = [
            (
                role = choice.region_role,
                order_index = choice.region_order_index,
                region_category = choice.region_category,
                recipe_family = choice.recipe_family,
                q = choice.q,
                order = choice.order,
                implementation_status = choice.implementation_status,
                buildability_status = choice.buildability_status,
                support_count = choice.support_count,
                selected = choice.selected,
            )
            for choice in policy.region_choices
        ],
        support_coverage = policy.construction_plan.support_coverage,
        outer_mismatch_is_outermost = policy.construction_plan.outer_mismatch_is_outermost,
        middle_contact_clean = policy.construction_plan.middle_contact_clean,
        active_builder_uses_policy = false,
    )
end

function _nested_diatomic_high_order_construction_piece_coverage(
    expected_support_indices::AbstractVector{Int},
    pieces::AbstractVector{_BondAlignedDiatomicHighOrderConstructionPiece3D},
)
    expected = Set(expected_support_indices)
    length(expected) == length(expected_support_indices) || throw(
        ArgumentError("diatomic high-order construction-piece coverage requires unique expected support"),
    )
    owned_counts = Dict{Int,Int}()
    piece_support_count = 0
    for piece in pieces
        piece_support_count += length(piece.support_indices)
        for index in piece.support_indices
            owned_counts[index] = get(owned_counts, index, 0) + 1
        end
    end
    owned = Set(keys(owned_counts))
    duplicate_count = sum(max(count - 1, 0) for count in values(owned_counts))
    missing_count = length(setdiff(expected, owned))
    outside_count = length(setdiff(owned, expected))
    coverage_ok = duplicate_count == 0 && missing_count == 0 && outside_count == 0
    return _BondAlignedDiatomicHighOrderConstructionPieceCoverage3D(
        length(expected),
        piece_support_count,
        length(owned),
        duplicate_count,
        missing_count,
        outside_count,
        coverage_ok,
    )
end

function _nested_diatomic_high_order_construction_piece(
    role::Symbol,
    piece_index::Int,
    primitive_family::Symbol,
    box::NTuple{3,UnitRange{Int}},
    inner_exclusion_box::Union{Nothing,NTuple{3,UnitRange{Int}}},
    dims::NTuple{3,Int};
    metadata::NamedTuple = (;),
)
    support_indices = isnothing(inner_exclusion_box) ?
        _nested_box_support_indices(box..., dims) :
        _nested_diatomic_shell_support_indices(box, inner_exclusion_box, dims)
    return _BondAlignedDiatomicHighOrderConstructionPiece3D(
        role,
        piece_index,
        primitive_family,
        box,
        inner_exclusion_box,
        support_indices,
        metadata,
    )
end

function _nested_diatomic_high_order_piece_descriptor(
    choice::_BondAlignedDiatomicHighOrderRecipeRegionChoice3D,
    region::_BondAlignedDiatomicAtomGrowthConstructionRegion3D,
    primitive_family::Symbol,
    pieces::Vector{_BondAlignedDiatomicHighOrderConstructionPiece3D};
    metadata::NamedTuple = (;),
)
    coverage = _nested_diatomic_high_order_construction_piece_coverage(
        region.support_indices,
        pieces,
    )
    return _BondAlignedDiatomicHighOrderRecipeRegionRealizationDescriptor3D(
        choice.region_role,
        choice.region_order_index,
        choice.recipe_family,
        choice.q,
        choice.order,
        length(region.support_indices),
        length(pieces),
        primitive_family,
        pieces,
        coverage,
        coverage.coverage_ok,
        metadata,
    )
end

function _nested_diatomic_high_order_contact_cap_descriptor(
    choice::_BondAlignedDiatomicHighOrderRecipeRegionChoice3D,
    region::_BondAlignedDiatomicAtomGrowthConstructionRegion3D,
    parent_box::NTuple{3,UnitRange{Int}},
)
    region.role == :contact_cap || throw(
        ArgumentError("diatomic high-order contact-cap descriptor requires a contact-cap region"),
    )
    dims = _nested_diatomic_atom_growth_plan_dims(parent_box)
    piece = _nested_diatomic_high_order_construction_piece(
        :contact_cap_slab,
        1,
        :contact_cap_owned_slab,
        region.box,
        nothing,
        dims;
        metadata = (
            contact_policy = get(region.metadata, :contact_policy, :unknown),
            contact_gap_count = get(region.metadata, :contact_gap_count, 0),
            exact_region_support = true,
        ),
    )
    return _nested_diatomic_high_order_piece_descriptor(
        choice,
        region,
        :contact_cap_owned_slab,
        [piece];
        metadata = (
            descriptor_scope = :middle_contact_cap,
            exact_full_coverage = true,
            feeds_future_primitive_family = :contact_cap_owned_slab,
        ),
    )
end

function _nested_diatomic_outer_mismatch_piece_role(axis_index::Int, side::Symbol)
    axis_symbol = (:x, :y, :z)[axis_index]
    return Symbol("outer_mismatch_", String(axis_symbol), "_", String(side), "_slab")
end

function _nested_diatomic_outer_mismatch_piece_box(
    outer_box::NTuple{3,UnitRange{Int}},
    inner_box::NTuple{3,UnitRange{Int}},
    axis_index::Int,
    side::Symbol,
)
    axis_range =
        side == :low ? (first(outer_box[axis_index]):(first(inner_box[axis_index]) - 1)) :
        side == :high ? ((last(inner_box[axis_index]) + 1):last(outer_box[axis_index])) :
        throw(ArgumentError("diatomic high-order outer mismatch side must be :low or :high"))
    isempty(axis_range) && return nothing
    return ntuple(3) do index
        index < axis_index ? inner_box[index] :
        index == axis_index ? axis_range :
        outer_box[index]
    end
end

function _nested_diatomic_high_order_outer_mismatch_descriptor(
    choice::_BondAlignedDiatomicHighOrderRecipeRegionChoice3D,
    region::_BondAlignedDiatomicAtomGrowthConstructionRegion3D,
    parent_box::NTuple{3,UnitRange{Int}},
)
    region.role == :outer_mismatch_shared_molecular_shell || throw(
        ArgumentError("diatomic high-order outer-mismatch descriptor requires an outer-mismatch region"),
    )
    inner_box = region.inner_exclusion_box
    isnothing(inner_box) && throw(
        ArgumentError("diatomic high-order outer-mismatch descriptor requires an inner exclusion box"),
    )
    dims = _nested_diatomic_atom_growth_plan_dims(parent_box)
    outer_box = region.box
    _nested_box_inside_parent(inner_box, outer_box) || throw(
        ArgumentError("diatomic high-order outer-mismatch inner box must be inside the outer box"),
    )
    pieces = _BondAlignedDiatomicHighOrderConstructionPiece3D[]
    piece_index = 1
    for axis_index in 1:3, side in (:low, :high)
        piece_box = _nested_diatomic_outer_mismatch_piece_box(
            outer_box,
            inner_box,
            axis_index,
            side,
        )
        isnothing(piece_box) && continue
        push!(
            pieces,
            _nested_diatomic_high_order_construction_piece(
                _nested_diatomic_outer_mismatch_piece_role(axis_index, side),
                piece_index,
                :outer_mismatch_boundary_slab,
                piece_box,
                nothing,
                dims;
                metadata = (
                    axis = (:x, :y, :z)[axis_index],
                    axis_index = axis_index,
                    side = side,
                    edge_assignment_policy = :earlier_axes_absorb_outer_edges,
                ),
            ),
        )
        piece_index += 1
    end
    return _nested_diatomic_high_order_piece_descriptor(
        choice,
        region,
        :outer_mismatch_boundary_slab_set,
        pieces;
        metadata = (
            descriptor_scope = :outermost_mismatch_shell,
            exact_full_coverage = true,
            feeds_future_primitive_family = :outer_mismatch_boundary_slab,
            slab_order_policy = :axis_ordered_disjoint_boundary_slabs,
        ),
    )
end

function _nested_diatomic_high_order_recipe_region_realization_descriptor(
    choice::_BondAlignedDiatomicHighOrderRecipeRegionChoice3D,
    region::_BondAlignedDiatomicAtomGrowthConstructionRegion3D,
    parent_box::NTuple{3,UnitRange{Int}},
)
    choice.recipe_family == :shared_contact_cap && return (
        _nested_diatomic_high_order_contact_cap_descriptor(choice, region, parent_box)
    )
    choice.recipe_family == :outermost_mismatch_shared_molecular_shell && return (
        _nested_diatomic_high_order_outer_mismatch_descriptor(choice, region, parent_box)
    )
    return nothing
end

function _nested_diatomic_high_order_projected_q_shell_candidate(
    choice::_BondAlignedDiatomicHighOrderRecipeRegionChoice3D,
)
    choice.region_role == :regular_shared_molecular_shell || return nothing
    choice.recipe_family == :shared_endcap_panel_exterior || return nothing
    return (
        realization_primitive = :projected_q_shell_boundary_comx_product_modes,
        mapped_primitive = :_nested_projected_q_shell_layer,
        mapped_primitive_status = :private_internal_metadata_candidate,
        primitive_family = :projected_q_shell,
        support_contract = :projected_q_shell_raw_boundary,
        coefficient_contract = :full_block_boundary_comx_product_mode_projection,
        seed_contract = :boundary_comx_product_modes_from_full_local_block_transform,
        mode_selection_rule = :any_axis_mode_index_first_or_last,
        cleanup_contract = :full_rank_symmetric_lowdin,
        sidecar_status = :not_yet_optimized_product_staged_for_pqs,
        active_builder_consumes = false,
        source_builder_consumes = false,
        fixed_block_consumes = false,
        qw_consumes = false,
        hamiltonian_consumes = false,
        current_transitional_implementation = :_nested_endcap_panel_shell_layer,
        current_transitional_recipe_family = :shared_endcap_panel_exterior,
        endcap_panel_status = :transitional_current_active_implementation,
        q = choice.q,
        order = choice.order,
    )
end

function _nested_diatomic_high_order_recipe_realization_mapping(
    choice::_BondAlignedDiatomicHighOrderRecipeRegionChoice3D,
    descriptor::Union{
        Nothing,
        _BondAlignedDiatomicHighOrderRecipeRegionRealizationDescriptor3D,
    } = nothing,
)
    if choice.recipe_family == :protected_atom_cubic_shell
        return (
            mapped_primitive = :_nested_bond_aligned_diatomic_sequence_for_box,
            mapped_primitive_status = :existing_internal_primitive,
            missing_implementation = nothing,
            existing_opt_in_route = nothing,
            notes = (
                primitive_scope = :atom_local_cubic_box,
                uses_existing_nested_complete_shell_sequence = true,
            ),
        )
    elseif choice.recipe_family == :shared_endcap_panel_exterior
        pqs_candidate = _nested_diatomic_high_order_projected_q_shell_candidate(choice)
        return (
            mapped_primitive = :_nested_endcap_panel_shell_layer,
            mapped_primitive_status = :existing_internal_primitive,
            missing_implementation = nothing,
            existing_opt_in_route = :shared_shell_layer_policy_endcap_panel_owned,
            notes = (
                primitive_scope = :one_cell_shared_exterior_shell,
                requires_product_doside = true,
                current_transitional_implementation = :_nested_endcap_panel_shell_layer,
                projected_q_shell_candidate = pqs_candidate,
            ),
        )
    elseif choice.recipe_family == :shared_contact_cap
        if !isnothing(descriptor)
            return (
                mapped_primitive = :_nested_diatomic_high_order_contact_cap_descriptor,
                mapped_primitive_status = :construction_piece_descriptor,
                missing_implementation = nothing,
                existing_opt_in_route = nothing,
                notes = (
                    primitive_scope = :middle_contact_cap,
                    descriptor_family = descriptor.primitive_family,
                    owned_unit_count = descriptor.owned_unit_count,
                ),
            )
        end
        return (
            mapped_primitive = nothing,
            mapped_primitive_status = :missing_region_primitive,
            missing_implementation = :contact_cap_region_constructor,
            existing_opt_in_route = nothing,
            notes = (
                primitive_scope = :middle_contact_cap,
                reason = :contact_cap_is_not_a_complete_shell_or_endcap_panel_exterior,
            ),
        )
    elseif choice.recipe_family == :outermost_mismatch_shared_molecular_shell
        if !isnothing(descriptor)
            return (
                mapped_primitive = :_nested_diatomic_high_order_outer_mismatch_descriptor,
                mapped_primitive_status = :construction_piece_descriptor,
                missing_implementation = nothing,
                existing_opt_in_route = nothing,
                notes = (
                    primitive_scope = :outermost_mismatch_shell,
                    descriptor_family = descriptor.primitive_family,
                    owned_unit_count = descriptor.owned_unit_count,
                ),
            )
        end
        return (
            mapped_primitive = nothing,
            mapped_primitive_status = :missing_region_primitive,
            missing_implementation = :outer_mismatch_shell_decomposition_or_owned_units,
            existing_opt_in_route = nothing,
            notes = (
                primitive_scope = :outermost_mismatch_shell,
                reason = :mismatch_region_can_be_multi_layer_or_directionally_unbalanced,
            ),
        )
    elseif choice.recipe_family == :transverse_annulus_exterior
        return (
            mapped_primitive = nothing,
            mapped_primitive_status = :missing_experimental_primitive,
            missing_implementation = :transverse_annulus_owned_unit_producer,
            existing_opt_in_route = nothing,
            notes = (
                primitive_scope = :shared_exterior_transverse_annulus,
                reason = :q5_annulus_is_optional_experimental_policy_metadata,
            ),
        )
    end
    throw(ArgumentError("unsupported diatomic high-order realization recipe family $(choice.recipe_family)"))
end

function _nested_diatomic_high_order_recipe_region_realization(
    choice::_BondAlignedDiatomicHighOrderRecipeRegionChoice3D,
    region::_BondAlignedDiatomicAtomGrowthConstructionRegion3D,
    parent_box::NTuple{3,UnitRange{Int}},
)
    choice.region_role == region.role && choice.region_order_index == region.order_index || throw(
        ArgumentError("diatomic high-order realization choice does not match its construction region"),
    )
    descriptor = _nested_diatomic_high_order_recipe_region_realization_descriptor(
        choice,
        region,
        parent_box,
    )
    mapping = _nested_diatomic_high_order_recipe_realization_mapping(choice, descriptor)
    return _BondAlignedDiatomicHighOrderRecipeRegionRealization3D(
        choice.region_role,
        choice.region_order_index,
        choice.region_category,
        choice.recipe_family,
        choice.q,
        choice.order,
        choice.support_count,
        choice.buildability_status,
        mapping.mapped_primitive,
        mapping.mapped_primitive_status,
        mapping.missing_implementation,
        false,
        mapping.existing_opt_in_route,
        descriptor,
        (
            implementation_status = choice.implementation_status,
            selected = choice.selected,
            realization_notes = mapping.notes,
            active_builder_uses_policy = false,
        ),
    )
end

function _nested_bond_aligned_diatomic_high_order_recipe_realization_audit(
    policy::_BondAlignedDiatomicHighOrderRecipePolicy3D,
)
    length(policy.region_choices) == length(policy.construction_plan.regions) || throw(
        ArgumentError("diatomic high-order realization audit requires one choice per construction region"),
    )
    realizations = [
        _nested_diatomic_high_order_recipe_region_realization(
            choice,
            region,
            policy.construction_plan.anatomy.recipe.parent_box,
        )
        for (choice, region) in zip(policy.region_choices, policy.construction_plan.regions)
    ]
    mapped_count = count(realization -> !isnothing(realization.mapped_primitive), realizations)
    missing_count = count(realization -> !isnothing(realization.missing_implementation), realizations)
    buildable_without_mapped = count(
        realization ->
            realization.buildability_status == :buildable_now &&
            isnothing(realization.mapped_primitive),
        realizations,
    )
    active_consumed = count(realization -> realization.active_builder_consumes, realizations)
    return _BondAlignedDiatomicHighOrderRecipeRealizationAudit3D(
        policy,
        realizations,
        mapped_count,
        missing_count,
        buildable_without_mapped,
        active_consumed,
        missing_count == 0,
        (
            active_builder_uses_policy = false,
            audit_only = true,
            no_source_builder_change = true,
            no_operator_change = true,
            support_coverage = policy.construction_plan.support_coverage,
        ),
    )
end

function _nested_diatomic_high_order_piece_coverage_diagnostics(
    coverage::_BondAlignedDiatomicHighOrderConstructionPieceCoverage3D,
)
    return (
        expected_support_count = coverage.expected_support_count,
        piece_support_count = coverage.piece_support_count,
        covered_support_count = coverage.covered_support_count,
        duplicate_count = coverage.duplicate_count,
        missing_count = coverage.missing_count,
        outside_count = coverage.outside_count,
        coverage_ok = coverage.coverage_ok,
    )
end

function _nested_diatomic_high_order_realization_descriptor_diagnostics(
    descriptor::Nothing,
)
    return nothing
end

function _nested_diatomic_high_order_realization_descriptor_diagnostics(
    descriptor::_BondAlignedDiatomicHighOrderRecipeRegionRealizationDescriptor3D,
)
    return (
        region_role = descriptor.region_role,
        region_order_index = descriptor.region_order_index,
        recipe_family = descriptor.recipe_family,
        q = descriptor.q,
        order = descriptor.order,
        parent_support_count = descriptor.parent_support_count,
        owned_unit_count = descriptor.owned_unit_count,
        primitive_family = descriptor.primitive_family,
        support_coverage =
            _nested_diatomic_high_order_piece_coverage_diagnostics(descriptor.support_coverage),
        exact_full_coverage = descriptor.exact_full_coverage,
        pieces = [
            (
                role = piece.role,
                piece_index = piece.piece_index,
                primitive_family = piece.primitive_family,
                support_count = length(piece.support_indices),
                box = piece.box,
                inner_exclusion_box = piece.inner_exclusion_box,
                metadata = piece.metadata,
            )
            for piece in descriptor.pieces
        ],
        metadata = descriptor.metadata,
    )
end

function _nested_bond_aligned_diatomic_high_order_recipe_realization_audit(
    plan::_BondAlignedDiatomicAtomGrowthConstructionPlan3D;
    kwargs...,
)
    return _nested_bond_aligned_diatomic_high_order_recipe_realization_audit(
        _nested_bond_aligned_diatomic_high_order_recipe_policy(plan; kwargs...),
    )
end

function _nested_bond_aligned_diatomic_high_order_recipe_realization_diagnostics(
    audit::_BondAlignedDiatomicHighOrderRecipeRealizationAudit3D,
)
    return (
        recipe_label = audit.policy.recipe_label,
        q_min = audit.policy.q_min,
        region_count = length(audit.region_realizations),
        mapped_region_count = audit.mapped_region_count,
        missing_region_count = audit.missing_region_count,
        buildable_without_mapped_primitive_count =
            audit.buildable_without_mapped_primitive_count,
        active_builder_consumed_region_count = audit.active_builder_consumed_region_count,
        descriptor_region_count =
            count(realization -> !isnothing(realization.realization_descriptor), audit.region_realizations),
        exact_descriptor_region_count =
            count(
                realization ->
                    !isnothing(realization.realization_descriptor) &&
                    realization.realization_descriptor.exact_full_coverage,
                audit.region_realizations,
            ),
        ready_for_opt_in_builder = audit.ready_for_opt_in_builder,
        region_realizations = [
            (
                role = realization.region_role,
                order_index = realization.region_order_index,
                region_category = realization.region_category,
                recipe_family = realization.recipe_family,
                q = realization.q,
                order = realization.order,
                support_count = realization.support_count,
                buildability_status = realization.buildability_status,
                mapped_primitive = realization.mapped_primitive,
                mapped_primitive_status = realization.mapped_primitive_status,
                missing_implementation = realization.missing_implementation,
                active_builder_consumes = realization.active_builder_consumes,
                existing_opt_in_route = realization.existing_opt_in_route,
                realization_descriptor =
                    _nested_diatomic_high_order_realization_descriptor_diagnostics(
                        realization.realization_descriptor,
                    ),
                metadata = realization.metadata,
            )
            for realization in audit.region_realizations
        ],
        support_coverage = audit.policy.construction_plan.support_coverage,
        active_builder_uses_policy = false,
    )
end

function _nested_diatomic_high_order_indices_coverage(
    expected_support_indices::AbstractVector{Int},
    owned_support_indices::AbstractVector{Int},
)
    expected = Set(expected_support_indices)
    length(expected) == length(expected_support_indices) || throw(
        ArgumentError("diatomic high-order support coverage requires unique expected support"),
    )
    owned_counts = Dict{Int,Int}()
    for index in owned_support_indices
        owned_counts[index] = get(owned_counts, index, 0) + 1
    end
    owned = Set(keys(owned_counts))
    duplicate_count = sum(max(count - 1, 0) for count in values(owned_counts))
    missing_count = length(setdiff(expected, owned))
    outside_count = length(setdiff(owned, expected))
    return _BondAlignedDiatomicHighOrderConstructionPieceCoverage3D(
        length(expected),
        length(owned_support_indices),
        length(owned),
        duplicate_count,
        missing_count,
        outside_count,
        duplicate_count == 0 && missing_count == 0 && outside_count == 0,
    )
end

function _nested_diatomic_assert_default_q4_opt_in_policy(
    policy::_BondAlignedDiatomicHighOrderRecipePolicy3D,
    audit::_BondAlignedDiatomicHighOrderRecipeRealizationAudit3D,
)
    policy.recipe_label == :mixed_atom_cubic_shared_endcap_panel || throw(
        ArgumentError(
            "experimental diatomic high-order source construction only consumes the atom-growth/endcap-panel policy",
        ),
    )
    policy.q_min == 4 || throw(
        ArgumentError("experimental diatomic high-order source construction requires q_min = 4"),
    )
    for choice in policy.region_choices
        if choice.region_category == :shared_exterior
            choice.recipe_family == :shared_endcap_panel_exterior || throw(
                ArgumentError(
                    "experimental diatomic high-order source construction only consumes shared endcap/panel exterior regions",
                ),
            )
            choice.q >= 4 || throw(
                ArgumentError("experimental diatomic high-order shared endcap/panel q must be >= 4"),
            )
            choice.order >= 1 || throw(
                ArgumentError("experimental diatomic high-order shared endcap/panel order must be positive"),
            )
        else
            choice.q == 4 && choice.order == 4 || throw(
                ArgumentError(
                    "experimental diatomic high-order source construction keeps non-shared regions at q=4, order=4",
                ),
            )
        end
    end
    audit.ready_for_opt_in_builder || throw(
        ArgumentError(
            "experimental diatomic high-order source construction requires a realization-ready policy",
        ),
    )
    for realization in audit.region_realizations
        isnothing(realization.missing_implementation) || throw(
            ArgumentError(
                "experimental diatomic high-order source construction cannot consume $(realization.region_role): missing $(realization.missing_implementation)",
            ),
        )
    end
    return nothing
end

function _nested_diatomic_validate_shared_shell_realization(realization::Symbol)
    realization in (:endcap_panel_owned, :projected_q_shell) || throw(
        ArgumentError(
            "experimental diatomic high-order shared-shell realization must be :endcap_panel_owned or :projected_q_shell",
        ),
    )
    return realization
end

function _nested_diatomic_projected_q_shell_retained_count(
    raw_source_dims::NTuple{3,Int},
)
    all(dim -> dim >= 3, raw_source_dims) || throw(
        ArgumentError("projected q-shell raw source dimensions must all be at least 3"),
    )
    inner_dims = ntuple(axis -> raw_source_dims[axis] - 2, 3)
    return prod(raw_source_dims) - prod(inner_dims)
end

function _nested_diatomic_projected_q_shell_adaptive_source_dimensions(
    basis,
    bundles::_CartesianNestedAxisBundles3D,
    region::_BondAlignedDiatomicAtomGrowthConstructionRegion3D,
    retention::CartesianNestedCompleteShellRetentionContract;
    bond_axis::Symbol,
    nside::Int,
    selected_q::Int,
    shared_shell_angular_resolution_scale::Float64,
)
    isnothing(region.inner_exclusion_box) && throw(
        ArgumentError("adaptive projected q-shell source dimensions require an inner exclusion box"),
    )
    adaptive_retention = _nested_diatomic_adaptive_shell_retention(
        basis,
        bundles,
        region.box,
        region.inner_exclusion_box,
        retention;
        nside = nside,
        shared_shell_angular_resolution_scale = shared_shell_angular_resolution_scale,
    )
    axis_selector_retained_counts = (
        adaptive_retention.chosen_x.retain,
        adaptive_retention.chosen_y.retain,
        adaptive_retention.chosen_z.retain,
    )
    raw_source_dims = ntuple(axis -> axis_selector_retained_counts[axis] + 2, 3)
    bond_axis_index = _nested_axis_index(bond_axis)
    transverse_raw_dims =
        Tuple(raw_source_dims[axis] for axis in 1:3 if axis != bond_axis_index)
    transverse_raw_dims[1] == transverse_raw_dims[2] || throw(
        ArgumentError("adaptive projected q-shell source dimensions require equal transverse raw dimensions"),
    )
    raw_q = transverse_raw_dims[1]
    raw_L = raw_source_dims[bond_axis_index]
    raw_q_matches_selected_q = raw_q == selected_q
    return (
        adaptive_retention = adaptive_retention,
        raw_source_dims = raw_source_dims,
        axis_selector_retained_counts = axis_selector_retained_counts,
        raw_q = raw_q,
        raw_L = raw_L,
        selected_q = selected_q,
        raw_q_matches_selected_q = raw_q_matches_selected_q,
        physical_box_lengths = length.(region.box),
        support_count = length(region.support_indices),
        pqs_retained_count =
            _nested_diatomic_projected_q_shell_retained_count(raw_source_dims),
        decomposition_status = raw_q_matches_selected_q ?
            :adaptive_broad_support_q_local_modes :
            :adaptive_raw_q_mismatch,
        broad_parent_boundary_reference = false,
        excluded_from_mvp_gate = !raw_q_matches_selected_q,
    )
end

function _nested_diatomic_high_order_descriptor_direct_coefficients(
    bundles::_CartesianNestedAxisBundles3D,
    descriptor::_BondAlignedDiatomicHighOrderRecipeRegionRealizationDescriptor3D,
)
    descriptor.exact_full_coverage || throw(
        ArgumentError("diatomic high-order descriptor construction requires exact region coverage"),
    )
    descriptor.support_coverage.coverage_ok || throw(
        ArgumentError("diatomic high-order descriptor construction requires clean piece support coverage"),
    )
    dims = _nested_axis_lengths(bundles)
    support_indices = reduce(
        vcat,
        (piece.support_indices for piece in descriptor.pieces);
        init = Int[],
    )
    return (
        support_indices = support_indices,
        coefficient_matrix = _nested_direct_box_coefficients(dims, support_indices),
    )
end

function _nested_diatomic_region_build_from_data(
    choice::_BondAlignedDiatomicHighOrderRecipeRegionChoice3D,
    realization::_BondAlignedDiatomicHighOrderRecipeRegionRealization3D,
    region::_BondAlignedDiatomicAtomGrowthConstructionRegion3D;
    mapped_primitive::Symbol,
    primitive_family::Symbol,
    built_support_indices::AbstractVector{Int},
    retained_count::Int,
    column_range::UnitRange{Int},
    built_object,
    metadata::NamedTuple,
)
    coverage = _nested_diatomic_high_order_indices_coverage(
        region.support_indices,
        built_support_indices,
    )
    coverage.coverage_ok || throw(
        ArgumentError(
            "experimental diatomic high-order source construction produced incomplete support for $(region.role)",
        ),
    )
    return _BondAlignedDiatomicHighOrderRecipeRegionSourceBuild3D(
        choice.region_role,
        choice.region_order_index,
        choice.region_category,
        choice.recipe_family,
        choice.q,
        choice.order,
        mapped_primitive,
        primitive_family,
        length(region.support_indices),
        length(built_support_indices),
        retained_count,
        column_range,
        coverage,
        true,
        nothing,
        true,
        built_object,
        merge(
            (
                realization_mapped_primitive_status =
                    realization.mapped_primitive_status,
                active_builder_uses_policy = true,
            ),
            metadata,
        ),
    )
end

function _nested_bond_aligned_diatomic_high_order_recipe_source_construction(
    basis,
    bundles::_CartesianNestedAxisBundles3D,
    policy::_BondAlignedDiatomicHighOrderRecipePolicy3D;
    nside::Int = 5,
    reference_fudge_factor::Float64 = 1.2,
    shared_shell_angular_resolution_scale::Float64 = 1.4,
    core_near_nucleus_protect_rows::Union{Symbol,Integer} = :auto,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
    packet_kernel::Symbol = :factorized_direct,
    build_sequence_packet::Bool = true,
    shared_shell_realization::Symbol = :endcap_panel_owned,
)
    return @timeg "diatomic.high_order_recipe_source.total" begin
        shared_shell_realization =
            _nested_diatomic_validate_shared_shell_realization(shared_shell_realization)
        normalized_packet_kernel = _nested_normalize_packet_kernel(packet_kernel)
        if shared_shell_realization == :projected_q_shell
            normalized_packet_kernel == :support_reference || throw(
                ArgumentError(
                    "projected q-shell opt-in source construction requires packet_kernel = :support_reference",
                ),
            )
            build_sequence_packet || throw(
                ArgumentError(
                    "projected q-shell opt-in source construction requires build_sequence_packet = true",
                ),
            )
        end
        isnothing(term_coefficients) && throw(
            ArgumentError("experimental diatomic high-order source construction requires explicit term coefficients"),
        )
        dims = _nested_axis_lengths(bundles)
        parent_box = (1:dims[1], 1:dims[2], 1:dims[3])
        policy.construction_plan.anatomy.recipe.parent_box == parent_box || throw(
            ArgumentError("experimental diatomic high-order source construction requires the policy parent_box to match the bundles"),
        )
        audit = _nested_bond_aligned_diatomic_high_order_recipe_realization_audit(policy)
        _nested_diatomic_assert_default_q4_opt_in_policy(policy, audit)

        bond_axis = policy.construction_plan.anatomy.recipe.bond_axis
        retention = _nested_resolve_complete_shell_retention(nside)
        protect_rows = _nested_diatomic_resolve_core_near_nucleus_protect_rows(
            core_near_nucleus_protect_rows,
            nside,
        )
        core_support_blocks = Vector{Vector{Int}}()
        core_coefficient_blocks = _CartesianCoefficientMap[]
        shared_shell_layers = _AbstractCartesianNestedShellLayer3D[]
        region_data = NamedTuple[]

        for (choice, region, realization) in zip(
            policy.region_choices,
            policy.construction_plan.regions,
            audit.region_realizations,
        )
            if choice.recipe_family == :protected_atom_cubic_shell
                sequence = _nested_bond_aligned_diatomic_sequence_for_box(
                    basis,
                    bundles,
                    region.box,
                    retention;
                    bond_axis = bond_axis,
                    nside = nside,
                    reference_fudge_factor = reference_fudge_factor,
                    core_near_nucleus_protect_rows = protect_rows,
                    term_coefficients = term_coefficients,
                    packet_kernel = packet_kernel,
                    build_packet = false,
                )
                push!(core_support_blocks, sequence.support_indices)
                push!(core_coefficient_blocks, sequence.coefficient_matrix)
                push!(
                    region_data,
                    (
                        kind = :core,
                        choice = choice,
                        region = region,
                        realization = realization,
                        mapped_primitive = :_nested_bond_aligned_diatomic_sequence_for_box,
                        primitive_family = :atom_local_complete_shell_sequence,
                        built_support_indices = sequence.support_indices,
                        retained_count = size(sequence.coefficient_matrix, 2),
                        built_object = sequence,
                        metadata = (
                            atom_side = get(region.metadata, :atom_side, :unknown),
                            nside = nside,
                            primitive_note = :existing_nested_sequence_for_box,
                        ),
                    ),
                )
            elseif choice.recipe_family == :shared_endcap_panel_exterior
                isnothing(region.inner_exclusion_box) && throw(
                    ArgumentError("shared exterior recipe region requires an inner exclusion box"),
                )
                if shared_shell_realization == :projected_q_shell
                    pqs_dimensions =
                        _nested_diatomic_projected_q_shell_adaptive_source_dimensions(
                            basis,
                            bundles,
                            region,
                            retention;
                            bond_axis = bond_axis,
                            nside = nside,
                            selected_q = choice.q,
                            shared_shell_angular_resolution_scale =
                                shared_shell_angular_resolution_scale,
                        )
                    layer = _nested_projected_q_shell_layer(
                        bundles,
                        region.box,
                        region.inner_exclusion_box;
                        bond_axis = bond_axis,
                        q = pqs_dimensions.raw_q,
                        L = pqs_dimensions.raw_L,
                        raw_source_dims = pqs_dimensions.raw_source_dims,
                        selected_q = pqs_dimensions.selected_q,
                        axis_selector_retained_counts =
                            pqs_dimensions.axis_selector_retained_counts,
                        decomposition_status = pqs_dimensions.decomposition_status,
                        broad_parent_boundary_reference =
                            pqs_dimensions.broad_parent_boundary_reference,
                        excluded_from_mvp_gate =
                            pqs_dimensions.excluded_from_mvp_gate,
                        packet_kernel = :support_reference,
                        term_coefficients = term_coefficients,
                        verify_factorized_reconstruction = false,
                    )
                    primitive_family = :projected_q_shell
                    mapped_primitive = :_nested_projected_q_shell_layer
                    layer_metadata = (
                        support_contract = layer.provenance.support_contract,
                        coefficient_contract = layer.provenance.coefficient_contract,
                        seed_contract =
                            layer.provenance.construction_contract,
                        cleanup_contract = :full_rank_symmetric_lowdin,
                        cleanup_method = layer.diagnostics.cleanup_method,
                        packet_kernel = layer.provenance.packet_kernel,
                        selected_q = pqs_dimensions.selected_q,
                        raw_source_dims = pqs_dimensions.raw_source_dims,
                        source_mode_dims = pqs_dimensions.raw_source_dims,
                        axis_selector_retained_counts =
                            pqs_dimensions.axis_selector_retained_counts,
                        raw_q = pqs_dimensions.raw_q,
                        raw_L = pqs_dimensions.raw_L,
                        raw_q_matches_selected_q =
                            pqs_dimensions.raw_q_matches_selected_q,
                        physical_box_lengths = pqs_dimensions.physical_box_lengths,
                        support_count = pqs_dimensions.support_count,
                        pqs_retained_count = pqs_dimensions.pqs_retained_count,
                        decomposition_status = pqs_dimensions.decomposition_status,
                        broad_parent_boundary_reference =
                            pqs_dimensions.broad_parent_boundary_reference,
                        excluded_from_mvp_gate =
                            pqs_dimensions.excluded_from_mvp_gate,
                        policy_q = choice.q,
                        policy_order = choice.order,
                        pqs_staged_unit_descriptor =
                            _nested_projected_q_shell_staged_unit_descriptor(layer),
                        pqs_product_staged_sidecar_available = false,
                        factorized_direct_allowed = false,
                        active_default_builder_changed = false,
                    )
                else
                    layer = _nested_endcap_panel_shell_layer(
                        bundles,
                        region.box,
                        region.inner_exclusion_box;
                        bond_axis = bond_axis,
                        q = choice.q,
                        L = choice.order,
                        packet_kernel = packet_kernel,
                        term_coefficients = term_coefficients,
                        verify_factorized_reconstruction = false,
                    )
                    primitive_family = :shared_endcap_panel_shell_layer
                    mapped_primitive = :_nested_endcap_panel_shell_layer
                    layer_metadata = (
                        support_contract = layer.provenance.support_contract,
                        coefficient_contract = layer.provenance.coefficient_contract,
                        packet_kernel = layer.provenance.packet_kernel,
                    )
                end
                push!(shared_shell_layers, layer)
                push!(
                    region_data,
                    (
                        kind = :shell,
                        choice = choice,
                        region = region,
                        realization = realization,
                        mapped_primitive = mapped_primitive,
                        primitive_family = primitive_family,
                        built_support_indices = layer.support_indices,
                        retained_count = size(layer.coefficient_matrix, 2),
                        built_object = layer,
                        metadata = layer_metadata,
                    ),
                )
            elseif choice.recipe_family in (
                :shared_contact_cap,
                :outermost_mismatch_shared_molecular_shell,
            )
                descriptor = realization.realization_descriptor
                isnothing(descriptor) && throw(
                    ArgumentError("descriptor-backed recipe region $(region.role) has no realization descriptor"),
                )
                direct = _nested_diatomic_high_order_descriptor_direct_coefficients(
                    bundles,
                    descriptor,
                )
                push!(core_support_blocks, direct.support_indices)
                push!(core_coefficient_blocks, direct.coefficient_matrix)
                push!(
                    region_data,
                    (
                        kind = :core,
                        choice = choice,
                        region = region,
                        realization = realization,
                        mapped_primitive = realization.mapped_primitive,
                        primitive_family = descriptor.primitive_family,
                        built_support_indices = direct.support_indices,
                        retained_count = size(direct.coefficient_matrix, 2),
                        built_object = direct,
                        metadata = (
                            descriptor_scope = get(descriptor.metadata, :descriptor_scope, :unknown),
                            owned_unit_count = descriptor.owned_unit_count,
                            exact_full_coverage = descriptor.exact_full_coverage,
                        ),
                    ),
                )
            else
                throw(
                    ArgumentError(
                        "experimental diatomic high-order source construction does not consume recipe family $(choice.recipe_family)",
                    ),
                )
            end
        end

        core_support_indices = reduce(vcat, core_support_blocks; init = Int[])
        core_coefficients = _nested_hcat_coefficient_maps(core_coefficient_blocks)
        sequence = _nested_shell_sequence_from_core_block(
            bundles,
            core_support_indices,
            core_coefficients,
            shared_shell_layers;
            term_coefficients = term_coefficients,
            packet_kernel = packet_kernel,
            build_packet = build_sequence_packet,
            verify_factorized_reconstruction = false,
        )

        core_column_start = 1
        core_ranges = UnitRange{Int}[]
        for block in core_coefficient_blocks
            column_count = size(block, 2)
            push!(core_ranges, core_column_start:(core_column_start + column_count - 1))
            core_column_start += column_count
        end
        shell_ranges = sequence.layer_column_ranges
        core_index = 0
        shell_index = 0
        region_builds = _BondAlignedDiatomicHighOrderRecipeRegionSourceBuild3D[]
        for data in region_data
            column_range = if data.kind == :core
                core_index += 1
                core_ranges[core_index]
            else
                shell_index += 1
                shell_ranges[shell_index]
            end
            push!(
                region_builds,
                _nested_diatomic_region_build_from_data(
                    data.choice,
                    data.realization,
                    data.region;
                    mapped_primitive = data.mapped_primitive,
                    primitive_family = data.primitive_family,
                    built_support_indices = data.built_support_indices,
                    retained_count = data.retained_count,
                    column_range = column_range,
                    built_object = data.built_object,
                    metadata = data.metadata,
                ),
            )
        end

        expected_support = _nested_box_support_indices(parent_box..., dims)
        support_coverage = _nested_diatomic_high_order_indices_coverage(
            expected_support,
            sequence.support_indices,
        )
        support_coverage.coverage_ok || throw(
            ArgumentError("experimental diatomic high-order source construction did not cover the parent support"),
        )
        shared_choices = [
            choice for choice in policy.region_choices if choice.region_category == :shared_exterior
        ]
        _BondAlignedDiatomicHighOrderRecipeSourceConstruction3D(
            basis,
            bundles,
            policy,
            audit,
            region_builds,
            sequence,
            shared_shell_layers,
            core_support_indices,
            core_coefficients,
            support_coverage,
            count(build -> build.built, region_builds),
            count(build -> !isnothing(build.unsupported_reason), region_builds),
            true,
            (
                active_builder_uses_policy = true,
                default_source_builder_changed = false,
                source_builder_entry = :_nested_bond_aligned_diatomic_high_order_recipe_source_construction,
                q_policy = :atom_growth_endcap_panel_shared_q_variable,
                q4_acceptance_fixture =
                    all(choice -> choice.q == 4 && choice.order == 4, policy.region_choices),
                non_shared_q_policy = :fixed_q4_order4,
                shared_q_values = Tuple(sort(unique(choice.q for choice in shared_choices))),
                shared_order_values =
                    Tuple(sort(unique(choice.order for choice in shared_choices))),
                packet_kernel = normalized_packet_kernel,
                build_sequence_packet = build_sequence_packet,
                shared_shell_realization = shared_shell_realization,
                projected_q_shell_opt_in = shared_shell_realization == :projected_q_shell,
            ),
        )
    end
end

function _nested_bond_aligned_diatomic_high_order_recipe_source_construction_diagnostics(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
)
    return (
        recipe_label = construction.policy.recipe_label,
        q_min = construction.policy.q_min,
        active_builder_consumes = construction.active_builder_consumes,
        active_builder_uses_policy = construction.active_builder_consumes,
        default_source_builder_changed = false,
        consumed_region_count = construction.consumed_region_count,
        unsupported_region_count = construction.unsupported_region_count,
        region_count = length(construction.region_builds),
        fixed_dimension = size(construction.sequence.coefficient_matrix, 2),
        parent_dimension = size(construction.sequence.coefficient_matrix, 1),
        support_coverage = _nested_diatomic_high_order_piece_coverage_diagnostics(
            construction.support_coverage,
        ),
        region_builds = [
            (
                role = build.region_role,
                order_index = build.region_order_index,
                region_category = build.region_category,
                recipe_family = build.recipe_family,
                q = build.q,
                order = build.order,
                mapped_primitive = build.mapped_primitive,
                primitive_family = build.primitive_family,
                region_support_count = build.region_support_count,
                built_support_count = build.built_support_count,
                retained_count = build.retained_count,
                column_range = build.column_range,
                support_coverage =
                    _nested_diatomic_high_order_piece_coverage_diagnostics(
                        build.support_coverage,
                    ),
                built = build.built,
                unsupported_reason = build.unsupported_reason,
                active_builder_consumes = build.active_builder_consumes,
                metadata = build.metadata,
            )
            for build in construction.region_builds
        ],
        metadata = construction.metadata,
    )
end

function _nested_diatomic_high_order_sequence_overlap_error(
    sequence::_CartesianNestedShellSequence3D,
)
    packet = sequence.packet
    isnothing(packet) && return nothing
    n = size(packet.overlap, 1)
    return norm(packet.overlap - Matrix{Float64}(I, n, n), Inf)
end

function _nested_diatomic_high_order_fixed_block_missing_fields(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
)
    missing = Symbol[]
    construction.support_coverage.coverage_ok || push!(missing, :support_coverage)
    sequence = construction.sequence
    size(sequence.coefficient_matrix, 1) == construction.support_coverage.expected_support_count ||
        push!(missing, :parent_dimension)
    size(sequence.coefficient_matrix, 2) > 0 || push!(missing, :fixed_dimension)
    isnothing(sequence.packet) && push!(missing, :sequence_packet)
    return missing
end

function _nested_diatomic_high_order_legacy_source_missing_fields(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
)
    missing = Symbol[]
    construction.active_builder_consumes || push!(missing, :active_builder_consumes)
    construction.support_coverage.coverage_ok || push!(missing, :support_coverage)
    append!(
        missing,
        [
            :split_geometry,
            :child_shell_retention_contract,
            :shared_shell_retention_contract,
            :child_sequences,
            :child_column_ranges,
            :midpoint_slab_column_range,
        ],
    )
    return unique(missing)
end

function _nested_bond_aligned_diatomic_high_order_recipe_source_fixed_block(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
)
    missing = _nested_diatomic_high_order_fixed_block_missing_fields(construction)
    isempty(missing) || throw(
        ArgumentError(
            "experimental diatomic high-order recipe source construction is not fixed-block ready; missing $(missing)",
        ),
    )
    return @timeg "diatomic.high_order_recipe_source.fixed_block" begin
        _nested_fixed_block(construction.sequence, construction.basis)
    end
end

function _nested_bond_aligned_diatomic_high_order_recipe_source_readiness(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D;
    build_fixed_block::Bool = false,
)
    sequence = construction.sequence
    packet = sequence.packet
    sequence_packet_available = !isnothing(packet)
    overlap_available =
        sequence_packet_available &&
        size(packet.overlap, 1) == size(sequence.coefficient_matrix, 2) &&
        size(packet.overlap, 2) == size(sequence.coefficient_matrix, 2)
    weights_available =
        sequence_packet_available &&
        length(packet.weights) == size(sequence.coefficient_matrix, 2) &&
        all(isfinite, packet.weights)
    fixed_missing = _nested_diatomic_high_order_fixed_block_missing_fields(construction)
    can_produce_fixed_block = isempty(fixed_missing)
    fixed_block =
        build_fixed_block && can_produce_fixed_block ?
        _nested_bond_aligned_diatomic_high_order_recipe_source_fixed_block(construction) :
        nothing
    nested_source_missing = _nested_diatomic_high_order_legacy_source_missing_fields(construction)
    return _BondAlignedDiatomicHighOrderRecipeSourceReadiness3D(
        construction,
        size(sequence.coefficient_matrix, 1),
        size(sequence.coefficient_matrix, 2),
        construction.support_coverage,
        sequence_packet_available,
        overlap_available,
        weights_available,
        _nested_diatomic_high_order_sequence_overlap_error(sequence),
        false,
        nested_source_missing,
        can_produce_fixed_block,
        fixed_missing,
        fixed_block,
        construction.active_builder_consumes,
        true,
        (
            source_object_contract = :legacy_diatomic_source_requires_split_geometry,
            fixed_block_contract = :existing_nested_shell_sequence_fixed_block,
            default_source_builder_changed = false,
            no_qw_or_hamiltonian_change = true,
        ),
    )
end

function _nested_bond_aligned_diatomic_high_order_recipe_source_readiness_diagnostics(
    readiness::_BondAlignedDiatomicHighOrderRecipeSourceReadiness3D,
)
    fixed_block = readiness.fixed_block
    return (
        parent_dimension = readiness.parent_dimension,
        fixed_dimension = readiness.fixed_dimension,
        support_coverage = _nested_diatomic_high_order_piece_coverage_diagnostics(
            readiness.support_coverage,
        ),
        sequence_packet_available = readiness.sequence_packet_available,
        overlap_available = readiness.overlap_available,
        weights_available = readiness.weights_available,
        overlap_error = readiness.overlap_error,
        can_produce_nested_source = readiness.can_produce_nested_source,
        nested_source_missing_fields = readiness.nested_source_missing_fields,
        can_produce_fixed_block = readiness.can_produce_fixed_block,
        fixed_block_missing_fields = readiness.fixed_block_missing_fields,
        fixed_block_built = !isnothing(fixed_block),
        fixed_block_backend = isnothing(fixed_block) ? nothing : fixed_block.gausslet_backend,
        fixed_block_dimension =
            isnothing(fixed_block) ? nothing : size(fixed_block.coefficient_matrix, 2),
        fixed_block_support_count =
            isnothing(fixed_block) ? nothing : length(fixed_block.support_indices),
        staged_by_center_sidecar_available =
            isnothing(fixed_block) ? nothing : !isnothing(fixed_block.staged_by_center_sidecar[]),
        active_builder_consumes = readiness.active_builder_consumes,
        default_builders_unchanged = readiness.default_builders_unchanged,
        metadata = readiness.metadata,
    )
end

function _nested_normalize_shared_shell_layer_policy(policy::Symbol)
    policy in (:complete_rectangular, :endcap_panel_owned) || throw(
        ArgumentError(
            "diatomic nested shared-shell layer policy must be :complete_rectangular or :endcap_panel_owned",
        ),
    )
    return policy
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
    shared_shell_layer_policy::Symbol = :complete_rectangular,
    shared_shell_endcap_panel_q::Int = 4,
    shared_shell_endcap_panel_L::Int = 4,
)
    return @timeg "diatomic.source.total" begin
        isnothing(term_coefficients) && throw(
            ArgumentError("diatomic source assembly requires explicit term coefficients"),
        )
        shared_layer_policy = _nested_normalize_shared_shell_layer_policy(shared_shell_layer_policy)
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
        shared_shell_layers =
            shared_layer_policy == :complete_rectangular ?
            _CartesianNestedCompleteShell3D[] :
            _AbstractCartesianNestedShellLayer3D[]
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
                if shared_layer_policy == :complete_rectangular
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
                            term_coefficients = term_coefficients,
                            packet_kernel = packet_kernel,
                            verify_factorized_reconstruction = false,
                        ),
                    )
                else
                    push!(
                        shared_shell_layers,
                        _nested_endcap_panel_shell_layer(
                            bundles,
                            current_box,
                            inner_box;
                            bond_axis = bond_axis,
                            q = shared_shell_endcap_panel_q,
                            L = shared_shell_endcap_panel_L,
                            term_coefficients = term_coefficients,
                            packet_kernel = packet_kernel,
                            verify_factorized_reconstruction = false,
                        ),
                    )
                end
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
