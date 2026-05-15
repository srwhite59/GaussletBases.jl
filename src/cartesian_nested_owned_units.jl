"""
    _CartesianNestedOwnedUnit3D

Internal experimental scaffolding for future nested owned-unit/endcap-panel
construction. An owned unit records one declared support region and the local
contraction map on that support. It is not a public OPCU API and is not wired
into nested source construction.
"""
struct _CartesianNestedOwnedUnit3D{M}
    role::Symbol
    support_indices::Vector{Int}
    coefficient_matrix::_CartesianCoefficientMap
    metadata::M
end

function _CartesianNestedOwnedUnit3D(
    role::Symbol,
    support_indices::AbstractVector{<:Integer},
    coefficient_matrix::AbstractMatrix{<:Real};
    metadata = (;),
)
    support = Int.(support_indices)
    coefficients = _cartesian_coefficient_map_storage(coefficient_matrix)
    size(coefficients, 1) == length(support) || throw(
        DimensionMismatch("nested owned-unit coefficient rows must match support-index count"),
    )
    return _CartesianNestedOwnedUnit3D{typeof(metadata)}(
        role,
        support,
        coefficients,
        metadata,
    )
end

"""
    _CartesianNestedOwnedUnitCoverageAudit

Compact report for owned-unit support coverage. `owned_support_count` counts
unique owned support indices. `duplicate_count` counts extra support-owner
occurrences beyond the first owner of each support index.
"""
struct _CartesianNestedOwnedUnitCoverageAudit
    expected_support_count::Int
    owned_support_count::Int
    duplicate_count::Int
    missing_count::Int
    outside_count::Int
    retained_count::Int
    coverage_ok::Bool
end

function _nested_owned_unit_coefficient_values(unit::_CartesianNestedOwnedUnit3D)
    matrix = unit.coefficient_matrix
    return matrix isa SparseArrays.SparseMatrixCSC ? SparseArrays.nonzeros(matrix) : vec(matrix)
end

function _nested_validate_owned_unit_coefficients(unit::_CartesianNestedOwnedUnit3D)
    retained_count = size(unit.coefficient_matrix, 2)
    retained_count >= 1 || throw(
        ArgumentError("nested owned-unit $(unit.role) must retain at least one column"),
    )
    all(isfinite, _nested_owned_unit_coefficient_values(unit)) || throw(
        ArgumentError("nested owned-unit $(unit.role) coefficient map must contain only finite values"),
    )
    return retained_count
end

"""
    _nested_owned_unit_coverage_audit(units, expected_support_indices)

Audit exact support ownership for future endcap/panel shell units. The helper
checks only declared support coverage and local contraction-map sanity; it does
not assemble packets or mutate source builders.
"""
function _nested_owned_unit_coverage_audit(
    units::AbstractVector{<:_CartesianNestedOwnedUnit3D},
    expected_support_indices::AbstractVector{<:Integer},
)
    expected_counts = Dict{Int,Int}()
    for index in expected_support_indices
        expected_index = Int(index)
        expected_counts[expected_index] = get(expected_counts, expected_index, 0) + 1
    end
    all(count == 1 for count in values(expected_counts)) || throw(
        ArgumentError("nested owned-unit coverage audit requires unique expected support indices"),
    )
    expected = Set(keys(expected_counts))

    owned_counts = Dict{Int,Int}()
    retained_count = 0
    for unit in units
        retained_count += _nested_validate_owned_unit_coefficients(unit)
        for index in unit.support_indices
            owned_counts[index] = get(owned_counts, index, 0) + 1
        end
    end
    owned = Set(keys(owned_counts))
    duplicate_count = sum(max(count - 1, 0) for count in values(owned_counts))
    missing_count = length(setdiff(expected, owned))
    outside_count = length(setdiff(owned, expected))
    coverage_ok = duplicate_count == 0 && missing_count == 0 && outside_count == 0
    return _CartesianNestedOwnedUnitCoverageAudit(
        length(expected),
        length(owned),
        duplicate_count,
        missing_count,
        outside_count,
        retained_count,
        coverage_ok,
    )
end
