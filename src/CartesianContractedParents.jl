module CartesianContractedParents

import ..GaussletBases: _CartesianCoefficientMap,
                         _CartesianNestedProductStagedByCenterUnit3D,
                         _CartesianNestedProjectedQShellStagedUnitDescriptor3D,
                         _NestedFixedBlock3D,
                         _cartesian_coefficient_map_storage,
                         _nested_staged_by_center_sidecar
import ..GaussletBases.CartesianParentGaussletBases:
    CartesianParentGaussletBasis3D,
    cartesian_parent_gausslet_basis,
    parent_dimension

export CartesianContractionUnit3D,
       CartesianContractionRule3D,
       CartesianContractionRuleInventory3D,
       CartesianContractedParent3D,
       CartesianContractedParentStructuralAudit,
       cartesian_contraction_rule,
       cartesian_contraction_rule_inventory,
       cartesian_contraction_unit_from_rule,
       cartesian_contracted_parent,
       contracted_parent_contraction_rules,
       contracted_parent_rule_inventory,
       contracted_parent_basis,
       contracted_parent_coefficients,
       contracted_parent_units,
       contracted_parent_metadata,
       contracted_parent_parent_dimension,
       contracted_parent_dimension,
       contraction_unit_role,
       contraction_unit_support_indices,
       contraction_unit_column_range,
       contraction_unit_metadata,
       contraction_unit_rule,
       contraction_rule_family,
       contraction_rule_kind,
       contraction_rule_support_summary,
       contraction_rule_column_range,
       contraction_rule_source_dimension,
       contraction_rule_retained_dimension,
       contraction_rule_transform_rule,
       contraction_rule_cleanup_rule,
       contraction_rule_metric_capability,
       contracted_parent_unit_column_ranges,
       contracted_parent_unit_support_indices,
       contracted_parent_support_indices,
       contracted_parent_structural_audit

"""
    CartesianContractionRule3D

Internal metadata record describing how a contracted-parent unit is meant to
be constructed. This is deliberately descriptive only: coefficient matrices,
fixed-block builders, metric packets, QW operators, and Hamiltonians remain
owned by their existing route-specific implementations.
"""
struct CartesianContractionRule3D{S,L,D,P}
    rule_family::Symbol
    kind::Symbol
    role::Union{Nothing,Symbol}
    support_indices::Vector{Int}
    support_summary::S
    local_geometry::L
    column_range::Union{Nothing,UnitRange{Int}}
    source_dimension::Int
    retained_dimension::Int
    transform_rule::Symbol
    cleanup_rule::Symbol
    metric_capability::Symbol
    diagnostics::D
    provenance::P
end

"""
    CartesianContractionRuleInventory3D

Metadata-only parent/rule inventory for contracted-parent construction rules.
It summarizes rule families, support coverage, retained dimensions, and metric
capabilities without building coefficient maps, metric packets, or operators.
"""
struct CartesianContractionRuleInventory3D{R,F,C,S,D,P}
    parent_dimension::Int
    contracted_dimension::Union{Nothing,Int}
    unit_count::Int
    rule_count::Int
    rules::R
    rule_family_counts::F
    metric_capabilities::C
    total_source_dimension::Int
    total_retained_dimension::Int
    support_summary::S
    rule_support_summaries::Vector{Any}
    every_unit_has_rule_metadata::Bool
    every_unit_rule_derivable::Bool
    metadata_only_rule_count::Int
    prototype_rule_count::Int
    any_metadata_only_rule::Bool
    any_prototype_rule::Bool
    diagnostics::D
    provenance::P
end

function _contraction_rule_support_summary(
    support_indices::AbstractVector{<:Integer};
    parent_dimension::Union{Nothing,Int} = nothing,
)
    values = Int[Int(index) for index in support_indices]
    unique_values = unique(values)
    duplicate_count = length(values) - length(unique_values)
    if isnothing(parent_dimension)
        return (
            parent_dimension = nothing,
            entry_count = length(values),
            unique_count = length(unique_values),
            duplicate_count,
            outside_count = nothing,
            missing_count = nothing,
            support_complete = nothing,
            coverage_checked = false,
        )
    end
    valid_range = 1:Int(parent_dimension)
    inside = Int[value for value in values if value in valid_range]
    unique_inside = unique(inside)
    return (
        parent_dimension = Int(parent_dimension),
        entry_count = length(values),
        unique_count = length(unique_inside),
        duplicate_count = length(inside) - length(unique_inside),
        outside_count = length(values) - length(inside),
        missing_count = Int(parent_dimension) - length(unique_inside),
        support_complete = length(unique_inside) == Int(parent_dimension),
        coverage_checked = true,
    )
end

function _staged_axis_rule_summary(axis)
    return (
        kind = axis.kind,
        fixed_index = axis.fixed_index,
        interval = axis.interval,
        coefficient_shape = size(axis.coefficient_matrix),
    )
end

function _product_staged_rule_family(kind::Symbol)
    kind == :product_doside && return :product_owned_unit
    kind == :support_dense && return :support_dense_fallback
    return :staged_unit
end

function _product_staged_transform_rule(kind::Symbol)
    kind == :product_doside && return :two_active_axis_product_doside
    kind == :support_dense && return :explicit_support_dense_coefficients
    return :staged_unit_coefficients
end

function _product_staged_cleanup_rule(kind::Symbol)
    kind == :product_doside && return :locally_orthonormal_product_doside
    kind == :support_dense && return :external_or_already_cleaned
    return :unspecified
end

function _product_staged_metric_capability(kind::Symbol)
    kind == :product_doside && return :product_staged_metric_contraction
    kind == :support_dense && return :support_local_product
    return :support_local_product
end

function _symbol_count_pairs(values)
    counts = Dict{Symbol,Int}()
    for value in values
        counts[value] = get(counts, value, 0) + 1
    end
    return sort!(collect(pairs(counts)); by = pair -> string(first(pair)))
end

_sorted_unique_symbols(values) = sort!(collect(Set(values)); by = string)

function _diagnostic_bool(diagnostics, name::Symbol)
    hasproperty(diagnostics, name) || return false
    return getproperty(diagnostics, name) === true
end

function _rule_metadata_only(rule::CartesianContractionRule3D)
    _diagnostic_bool(rule.diagnostics, :metadata_only) && return true
    hasproperty(rule.diagnostics, :original_diagnostics) || return false
    return _diagnostic_bool(rule.diagnostics.original_diagnostics, :metadata_only)
end

function _rule_prototype_only(rule::CartesianContractionRule3D)
    _diagnostic_bool(rule.diagnostics, :prototype_only) && return true
    rule.metric_capability == :pqs_low_order_product_metric_prototype && return true
    hasproperty(rule.diagnostics, :original_diagnostics) || return false
    return _diagnostic_bool(rule.diagnostics.original_diagnostics, :prototype_only)
end

function cartesian_contraction_rule(
    unit::_CartesianNestedProductStagedByCenterUnit3D;
    parent_dimension::Union{Nothing,Int} = nothing,
)
    diagnostics = merge(
        (
            source = :nested_product_staged_by_center_unit,
            coefficient_shape = size(unit.coefficient_matrix),
            axis_function_count = length(unit.axis_function_indices),
            metadata_only = false,
            prototype_only = false,
            coefficient_contract = hasproperty(unit.provenance, :coefficient_contract) ?
                                   unit.provenance.coefficient_contract :
                                   nothing,
        ),
        unit.diagnostics,
    )
    local_geometry = (
        axes = map(_staged_axis_rule_summary, unit.axes),
        axis_function_index_count = length(unit.axis_function_indices),
    )
    return CartesianContractionRule3D(
        _product_staged_rule_family(unit.kind),
        unit.kind,
        unit.role,
        copy(unit.support_indices),
        _contraction_rule_support_summary(
            unit.support_indices;
            parent_dimension,
        ),
        local_geometry,
        unit.column_range,
        size(unit.coefficient_matrix, 1),
        length(unit.column_range),
        _product_staged_transform_rule(unit.kind),
        _product_staged_cleanup_rule(unit.kind),
        _product_staged_metric_capability(unit.kind),
        diagnostics,
        (
            source = :nested_product_staged_by_center_sidecar,
            staged_unit = unit,
            original_provenance = unit.provenance,
        ),
    )
end

function cartesian_contraction_rule(
    descriptor::_CartesianNestedProjectedQShellStagedUnitDescriptor3D;
    parent_dimension::Union{Nothing,Int} = nothing,
)
    full_block_dimension = prod(length.(descriptor.current_box))
    diagnostics = (
        source = :projected_q_shell_staged_unit_descriptor,
        support_count = descriptor.support_count,
        boundary_mode_count = descriptor.mode_count,
        retained_count = descriptor.retained_count,
        boundary_column_count = length(descriptor.boundary_column_indices),
        cleanup_rank_count = descriptor.cleanup_rank_count,
        cleanup_rank_drop_count = descriptor.cleanup_rank_drop_count,
        cleanup_cutoff = descriptor.cleanup_cutoff,
        metadata_only = true,
        prototype_only = true,
        contracted_parent_unit_installed = false,
        non_contracts = descriptor.non_contracts,
        active_consumption = descriptor.active_consumption,
        original_diagnostics = descriptor.diagnostics,
    )
    local_geometry = (
        current_box = descriptor.current_box,
        inner_box = descriptor.inner_box,
        bond_axis = descriptor.bond_axis,
        q = descriptor.q,
        L = descriptor.L,
        axis_intervals = descriptor.axis_intervals,
        axis_local_coefficient_shapes = map(size, descriptor.axis_local_coefficients),
        cleanup_matrix_size = descriptor.cleanup_matrix_size,
    )
    return CartesianContractionRule3D(
        :projected_q_shell_boundary_modes,
        descriptor.kind,
        descriptor.role,
        copy(descriptor.support_indices),
        _contraction_rule_support_summary(
            descriptor.support_indices;
            parent_dimension,
        ),
        local_geometry,
        nothing,
        full_block_dimension,
        descriptor.retained_count,
        :boundary_comx_product_modes_raw_boundary_projection,
        :full_rank_symmetric_lowdin,
        :pqs_low_order_product_metric_prototype,
        diagnostics,
        (
            source = :projected_q_shell_staged_unit_descriptor,
            descriptor,
        ),
    )
end

function cartesian_contraction_rule_inventory(
    rules::AbstractVector{<:CartesianContractionRule3D};
    parent_dimension::Integer,
    contracted_dimension::Union{Nothing,Integer} = nothing,
    unit_count::Integer = length(rules),
    every_unit_has_rule_metadata::Bool = false,
    every_unit_rule_derivable::Bool = false,
    provenance = (; source = :cartesian_contraction_rule_collection),
)
    rule_values = collect(rules)
    parent_dim = Int(parent_dimension)
    contracted_dim = isnothing(contracted_dimension) ? nothing : Int(contracted_dimension)
    support_indices = Int[]
    for rule in rule_values
        append!(support_indices, rule.support_indices)
    end
    family_counts = _symbol_count_pairs(rule.rule_family for rule in rule_values)
    metric_capabilities = _sorted_unique_symbols(rule.metric_capability for rule in rule_values)
    metadata_only_count = count(_rule_metadata_only, rule_values)
    prototype_count = count(_rule_prototype_only, rule_values)
    all_rules_have_column_ranges = all(rule -> !isnothing(rule.column_range), rule_values)
    diagnostics = (
        source = :cartesian_contraction_rule_inventory,
        parent_level_unit_inventory =
            every_unit_rule_derivable &&
            Int(unit_count) == length(rule_values) &&
            all_rules_have_column_ranges,
        all_rules_have_column_ranges,
        rule_family_counts = family_counts,
        metric_capabilities,
        metadata_only_rule_count = metadata_only_count,
        prototype_rule_count = prototype_count,
        any_metadata_only_rule = metadata_only_count > 0,
        any_prototype_rule = prototype_count > 0,
        q_shell_rule_present = any(
            rule -> rule.rule_family == :projected_q_shell_boundary_modes,
            rule_values,
        ),
        q_shell_installed_as_contracted_parent_unit = any(
            rule ->
                rule.rule_family == :projected_q_shell_boundary_modes &&
                !isnothing(rule.column_range),
            rule_values,
        ),
    )
    return CartesianContractionRuleInventory3D(
        parent_dim,
        contracted_dim,
        Int(unit_count),
        length(rule_values),
        rule_values,
        family_counts,
        metric_capabilities,
        sum(rule.source_dimension for rule in rule_values),
        sum(rule.retained_dimension for rule in rule_values),
        _contraction_rule_support_summary(support_indices; parent_dimension = parent_dim),
        Any[rule.support_summary for rule in rule_values],
        every_unit_has_rule_metadata,
        every_unit_rule_derivable,
        metadata_only_count,
        prototype_count,
        metadata_only_count > 0,
        prototype_count > 0,
        diagnostics,
        provenance,
    )
end

"""
    CartesianContractionUnit3D

Internal provenance record for a group of contracted-parent columns.

The global coefficient matrix on `CartesianContractedParent3D` remains the
source of truth. Unit support and column ranges are structural metadata only:
support may overlap, may be incomplete, and does not imply orthonormality.
"""
struct CartesianContractionUnit3D{M}
    role::Symbol
    support_indices::Vector{Int}
    column_range::UnitRange{Int}
    metadata::M
end

function CartesianContractionUnit3D(
    role::Symbol,
    support_indices::AbstractVector{<:Integer},
    column_range::UnitRange{<:Integer};
    metadata = (;),
)
    isempty(column_range) && throw(
        ArgumentError("Cartesian contraction units require a non-empty column range"),
    )
    return CartesianContractionUnit3D(
        role,
        Int[Int(index) for index in support_indices],
        Int(first(column_range)):Int(last(column_range)),
        metadata,
    )
end

function cartesian_contraction_unit_from_rule(
    rule::CartesianContractionRule3D,
    payload::_CartesianNestedProductStagedByCenterUnit3D;
    metadata = (;),
)
    rule.kind == payload.kind || throw(
        ArgumentError("contraction rule kind $(rule.kind) does not match staged payload kind $(payload.kind)"),
    )
    rule.role == payload.role || throw(
        ArgumentError("contraction rule role $(rule.role) does not match staged payload role $(payload.role)"),
    )
    rule.support_indices == payload.support_indices || throw(
        ArgumentError("contraction rule support does not match staged payload support"),
    )
    rule.column_range == payload.column_range || throw(
        ArgumentError("contraction rule column range does not match staged payload column range"),
    )
    rule.retained_dimension == length(payload.column_range) || throw(
        ArgumentError("contraction rule retained dimension does not match staged payload column range"),
    )
    return CartesianContractionUnit3D(
        payload.role,
        payload.support_indices,
        payload.column_range;
        metadata = merge(
            (
                source = :nested_product_staged_by_center_sidecar,
                rule_driven_unit_creation = true,
                rule_family = rule.rule_family,
                rule_kind = rule.kind,
                rule_metric_capability = rule.metric_capability,
                staged_by_center_kind = payload.kind,
                staged_by_center_unit = payload,
                contraction_rule = rule,
            ),
            metadata,
        ),
    )
end

function cartesian_contraction_rule(
    unit::CartesianContractionUnit3D;
    parent_dimension::Union{Nothing,Int} = nothing,
)
    if hasproperty(unit.metadata, :contraction_rule)
        return unit.metadata.contraction_rule
    elseif hasproperty(unit.metadata, :staged_by_center_unit)
        return cartesian_contraction_rule(
            unit.metadata.staged_by_center_unit;
            parent_dimension,
        )
    end
    diagnostics = (
        source = :cartesian_contraction_unit,
        metadata = unit.metadata,
        metadata_only = false,
        prototype_only = false,
    )
    return CartesianContractionRule3D(
        :support_dense_fallback,
        :support_dense,
        unit.role,
        copy(unit.support_indices),
        _contraction_rule_support_summary(
            unit.support_indices;
            parent_dimension,
        ),
        (;),
        unit.column_range,
        length(unit.support_indices),
        length(unit.column_range),
        :explicit_support_dense_coefficients,
        :external_or_already_cleaned,
        :support_local_product,
        diagnostics,
        (source = :cartesian_contraction_unit,),
    )
end

"""
    CartesianContractedParent3D

Internal identity object for columns formed as linear combinations of a full
Cartesian parent gausslet lattice.

This object deliberately stores no backend state, overlap/H/V matrices,
Gaussian supplements, residual columns, or operator packets. The coefficient
matrix is the source of truth; units are provenance only.
"""
struct CartesianContractedParent3D{P<:CartesianParentGaussletBasis3D,C<:_CartesianCoefficientMap,U<:AbstractVector,M}
    parent::P
    coefficient_matrix::C
    units::U
    metadata::M
end

function _validate_unit_column_ranges(
    units::AbstractVector{<:CartesianContractionUnit3D},
    contracted_dimension::Int,
)
    for unit in units
        first(unit.column_range) >= 1 || throw(
            ArgumentError("contraction unit column ranges must start inside contracted columns"),
        )
        last(unit.column_range) <= contracted_dimension || throw(
            ArgumentError("contraction unit column ranges must lie inside contracted columns"),
        )
    end
    return units
end

function CartesianContractedParent3D(
    parent::CartesianParentGaussletBasis3D,
    coefficients::AbstractMatrix{<:Real};
    units::AbstractVector{<:CartesianContractionUnit3D} = CartesianContractionUnit3D[],
    metadata = (;),
)
    coefficient_matrix = _cartesian_coefficient_map_storage(coefficients)
    size(coefficient_matrix, 1) == parent_dimension(parent) || throw(
        DimensionMismatch("contracted parent coefficient rows must match parent dimension"),
    )
    unit_values = CartesianContractionUnit3D[units...]
    _validate_unit_column_ranges(unit_values, size(coefficient_matrix, 2))
    return CartesianContractedParent3D(
        parent,
        coefficient_matrix,
        unit_values,
        metadata,
    )
end

function _contracted_parent_units_from_staged_sidecar(sidecar)
    if hasproperty(sidecar, :units)
        parent_dim = hasproperty(sidecar, :dims) ? prod(sidecar.dims) : nothing
        return CartesianContractionUnit3D[
            cartesian_contraction_unit_from_rule(
                cartesian_contraction_rule(unit; parent_dimension = parent_dim),
                unit,
            ) for unit in sidecar.units
        ]
    elseif hasproperty(sidecar, :block_column_ranges)
        return CartesianContractionUnit3D[
            CartesianContractionUnit3D(
                Symbol(:staged_block_, block_index),
                sidecar.block_support_indices[block_index],
                sidecar.block_column_ranges[block_index];
                metadata = (
                    source = :nested_staged_by_center_sidecar,
                    staged_by_center_kind = :support_dense,
                    staged_by_center_block_index = block_index,
                    staged_by_center_coefficients = sidecar.block_coefficients[block_index],
                    staged_by_center_support_states = sidecar.block_support_states[block_index],
                ),
            ) for block_index in eachindex(sidecar.block_column_ranges)
        ]
    end
    throw(ArgumentError("unsupported staged by-center sidecar for Cartesian contracted parent adapter"))
end

function CartesianContractedParent3D(
    fixed_block::_NestedFixedBlock3D;
    metadata = (source = :nested_fixed_block,),
)
    parent = cartesian_parent_gausslet_basis(fixed_block)
    coefficients = fixed_block.coefficient_matrix
    staged_sidecar = _nested_staged_by_center_sidecar(fixed_block)
    if !isnothing(staged_sidecar)
        path = hasproperty(staged_sidecar, :units) ? :product_staged_factorized : :staged_factorized
        return CartesianContractedParent3D(
            parent,
            coefficients;
            units = _contracted_parent_units_from_staged_sidecar(staged_sidecar),
            metadata = merge(
                metadata,
                (
                    staged_by_center_path = path,
                    staged_by_center_sidecar = staged_sidecar,
                ),
            ),
        )
    end
    unit = CartesianContractionUnit3D(
        :nested_fixed_block,
        fixed_block.support_indices,
        1:size(coefficients, 2);
        metadata = (source = :nested_fixed_block,),
    )
    return CartesianContractedParent3D(
        parent,
        coefficients;
        units = [unit],
        metadata,
    )
end

cartesian_contracted_parent(
    parent::CartesianParentGaussletBasis3D,
    coefficients::AbstractMatrix{<:Real};
    kwargs...,
) = CartesianContractedParent3D(parent, coefficients; kwargs...)

cartesian_contracted_parent(fixed_block::_NestedFixedBlock3D; kwargs...) =
    CartesianContractedParent3D(fixed_block; kwargs...)

contracted_parent_basis(parent::CartesianContractedParent3D) = parent.parent
contracted_parent_coefficients(parent::CartesianContractedParent3D) = parent.coefficient_matrix
contracted_parent_units(parent::CartesianContractedParent3D) = parent.units
contracted_parent_metadata(parent::CartesianContractedParent3D) = parent.metadata
contracted_parent_parent_dimension(parent::CartesianContractedParent3D) =
    parent_dimension(parent.parent)
contracted_parent_dimension(parent::CartesianContractedParent3D) =
    size(parent.coefficient_matrix, 2)

function contracted_parent_contraction_rules(parent::CartesianContractedParent3D)
    parent_dim = contracted_parent_parent_dimension(parent)
    return [
        contraction_unit_rule(unit; parent_dimension = parent_dim) for unit in parent.units
    ]
end

function contracted_parent_rule_inventory(parent::CartesianContractedParent3D)
    units = contracted_parent_units(parent)
    return cartesian_contraction_rule_inventory(
        contracted_parent_contraction_rules(parent);
        parent_dimension = contracted_parent_parent_dimension(parent),
        contracted_dimension = contracted_parent_dimension(parent),
        unit_count = length(units),
        every_unit_has_rule_metadata = all(
            unit -> hasproperty(unit.metadata, :contraction_rule),
            units,
        ),
        every_unit_rule_derivable = true,
        provenance = (
            source = :cartesian_contracted_parent,
            contracted_parent_metadata = parent.metadata,
        ),
    )
end

contraction_unit_role(unit::CartesianContractionUnit3D) = unit.role
contraction_unit_support_indices(unit::CartesianContractionUnit3D) = unit.support_indices
contraction_unit_column_range(unit::CartesianContractionUnit3D) = unit.column_range
contraction_unit_metadata(unit::CartesianContractionUnit3D) = unit.metadata
contraction_unit_rule(
    unit::CartesianContractionUnit3D;
    parent_dimension::Union{Nothing,Int} = nothing,
) = cartesian_contraction_rule(unit; parent_dimension)

contraction_rule_family(rule::CartesianContractionRule3D) = rule.rule_family
contraction_rule_kind(rule::CartesianContractionRule3D) = rule.kind
contraction_rule_support_summary(rule::CartesianContractionRule3D) = rule.support_summary
contraction_rule_column_range(rule::CartesianContractionRule3D) = rule.column_range
contraction_rule_source_dimension(rule::CartesianContractionRule3D) = rule.source_dimension
contraction_rule_retained_dimension(rule::CartesianContractionRule3D) = rule.retained_dimension
contraction_rule_transform_rule(rule::CartesianContractionRule3D) = rule.transform_rule
contraction_rule_cleanup_rule(rule::CartesianContractionRule3D) = rule.cleanup_rule
contraction_rule_metric_capability(rule::CartesianContractionRule3D) = rule.metric_capability

contracted_parent_unit_column_ranges(parent::CartesianContractedParent3D) =
    UnitRange{Int}[unit.column_range for unit in parent.units]

contracted_parent_unit_support_indices(parent::CartesianContractedParent3D) =
    Vector{Int}[copy(unit.support_indices) for unit in parent.units]

function contracted_parent_support_indices(parent::CartesianContractedParent3D)
    indices = Int[]
    for unit in parent.units
        append!(indices, unit.support_indices)
    end
    return indices
end

struct CartesianContractedParentStructuralAudit
    parent_dimension::Int
    contracted_dimension::Int
    unit_count::Int
    support_entry_count::Int
    unique_support_count::Int
    duplicate_support_count::Int
    missing_support_count::Int
    outside_support_count::Int
    column_entry_count::Int
    unique_column_count::Int
    duplicate_column_count::Int
    missing_column_count::Int
    outside_column_count::Int
    support_complete::Bool
    column_ranges_cover_contract::Bool
    structural_ok::Bool
end

function _entry_counts(values::AbstractVector{Int}, valid_range::UnitRange{Int})
    inside = Int[value for value in values if value in valid_range]
    outside_count = length(values) - length(inside)
    unique_inside = unique(inside)
    duplicate_count = length(inside) - length(unique_inside)
    missing_count = length(valid_range) - length(unique_inside)
    return (
        entry_count = length(values),
        unique_count = length(unique_inside),
        duplicate_count = duplicate_count,
        missing_count = missing_count,
        outside_count = outside_count,
    )
end

function contracted_parent_structural_audit(parent::CartesianContractedParent3D)
    parent_dim = contracted_parent_parent_dimension(parent)
    contracted_dim = contracted_parent_dimension(parent)
    support_counts = _entry_counts(
        contracted_parent_support_indices(parent),
        1:parent_dim,
    )
    columns = Int[]
    for unit in parent.units
        append!(columns, collect(unit.column_range))
    end
    column_counts = _entry_counts(columns, 1:contracted_dim)
    support_complete = support_counts.missing_count == 0 && support_counts.outside_count == 0
    column_ranges_cover_contract =
        column_counts.missing_count == 0 &&
        column_counts.outside_count == 0 &&
        column_counts.duplicate_count == 0
    structural_ok = support_counts.outside_count == 0 && column_ranges_cover_contract
    return CartesianContractedParentStructuralAudit(
        parent_dim,
        contracted_dim,
        length(parent.units),
        support_counts.entry_count,
        support_counts.unique_count,
        support_counts.duplicate_count,
        support_counts.missing_count,
        support_counts.outside_count,
        column_counts.entry_count,
        column_counts.unique_count,
        column_counts.duplicate_count,
        column_counts.missing_count,
        column_counts.outside_count,
        support_complete,
        column_ranges_cover_contract,
        structural_ok,
    )
end

end
