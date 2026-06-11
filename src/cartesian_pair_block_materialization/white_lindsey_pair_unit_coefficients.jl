# Pair-level White--Lindsey unit coefficient preparation.

"""
    white_lindsey_boundary_stratum_pair_unit_coefficients(unit_pair)
    white_lindsey_boundary_stratum_pair_unit_coefficients(left_unit, right_unit; ...)

Gather left/right White--Lindsey boundary-stratum unit coefficient maps for one
retained-unit pair. This prepares later pair-operator passes without building
overlap blocks, one-body blocks, Hamiltonian data, exports, artifacts, IDA/MWG
data, or Coulomb.
"""
function white_lindsey_boundary_stratum_pair_unit_coefficients(
    unit_pair::CUP.UnitPairRecord,
)
    cache = Dict{Symbol,Any}()
    left_coefficients =
        _white_lindsey_unit_coefficients_from_local_cache(cache, unit_pair.left_unit)
    right_coefficients =
        _white_lindsey_unit_coefficients_from_local_cache(cache, unit_pair.right_unit)
    return _white_lindsey_pair_unit_coefficients_result(
        unit_pair,
        left_coefficients,
        right_coefficients,
        length(cache),
    )
end

function white_lindsey_boundary_stratum_pair_unit_coefficients(
    left_unit_or_descriptor,
    right_unit_or_descriptor;
    pair_key = _white_lindsey_explicit_pair_key(
        left_unit_or_descriptor,
        right_unit_or_descriptor,
    ),
    pair_index = nothing,
    pair_family = _white_lindsey_explicit_pair_family(
        left_unit_or_descriptor,
        right_unit_or_descriptor,
    ),
)
    cache = Dict{Symbol,Any}()
    left_coefficients =
        _white_lindsey_unit_coefficients_from_local_cache(
            cache,
            left_unit_or_descriptor,
        )
    right_coefficients =
        _white_lindsey_unit_coefficients_from_local_cache(
            cache,
            right_unit_or_descriptor,
        )
    pair_metadata = (;
        pair_key,
        pair_index,
        pair_family,
        left_unit_key = pair_key[1],
        right_unit_key = pair_key[2],
        left_unit_kind =
            _white_lindsey_descriptor_property(left_unit_or_descriptor, :unit_kind),
        right_unit_kind =
            _white_lindsey_descriptor_property(right_unit_or_descriptor, :unit_kind),
    )
    return _white_lindsey_pair_unit_coefficients_result(
        pair_metadata,
        left_coefficients,
        right_coefficients,
        length(cache),
    )
end

function white_lindsey_boundary_stratum_pair_unit_coefficients(
    unit_pairs::Tuple{Vararg{CUP.UnitPairRecord}},
)
    return _white_lindsey_boundary_stratum_pair_unit_coefficients(unit_pairs)
end

function white_lindsey_boundary_stratum_pair_unit_coefficients(
    unit_pairs::AbstractVector{<:CUP.UnitPairRecord},
)
    return _white_lindsey_boundary_stratum_pair_unit_coefficients(unit_pairs)
end

function _white_lindsey_boundary_stratum_pair_unit_coefficients(unit_pairs)
    cache = Dict{Symbol,Any}()
    results = map(
        unit_pair ->
            _white_lindsey_pair_unit_coefficients_result(
                unit_pair,
                _white_lindsey_unit_coefficients_from_local_cache(
                    cache,
                    unit_pair.left_unit,
                ),
                _white_lindsey_unit_coefficients_from_local_cache(
                    cache,
                    unit_pair.right_unit,
                ),
                length(cache),
            ),
        unit_pairs,
    )
    if unit_pairs isa Tuple
        results = Tuple(results)
    end
    return _white_lindsey_pair_unit_coefficients_batch_result(
        results,
        length(cache),
    )
end

function _white_lindsey_unit_coefficients_from_local_cache(cache, unit_or_descriptor)
    unit_key =
        _white_lindsey_descriptor_property(unit_or_descriptor, :unit_key)
    if !isnothing(unit_key) && haskey(cache, unit_key)
        return cache[unit_key]
    end
    coefficients =
        white_lindsey_boundary_stratum_unit_coefficients(unit_or_descriptor)
    isnothing(unit_key) || (cache[unit_key] = coefficients)
    return coefficients
end

function _white_lindsey_pair_unit_coefficients_result(
    pair_source,
    left_coefficients,
    right_coefficients,
    cache_entry_count::Int,
)
    left_materialized =
        _white_lindsey_unit_coefficients_materialized(left_coefficients)
    right_materialized =
        _white_lindsey_unit_coefficients_materialized(right_coefficients)
    status, blocker =
        _white_lindsey_pair_unit_coefficients_status(
            left_materialized,
            right_materialized,
        )
    return (;
        object_kind =
            :white_lindsey_boundary_stratum_pair_unit_coefficients,
        status,
        blocker,
        pair_key = _white_lindsey_descriptor_property(pair_source, :pair_key),
        pair_index = _white_lindsey_descriptor_property(pair_source, :pair_index),
        pair_family =
            _white_lindsey_descriptor_property(pair_source, :pair_family),
        left_unit_key =
            _white_lindsey_descriptor_property(pair_source, :left_unit_key),
        right_unit_key =
            _white_lindsey_descriptor_property(pair_source, :right_unit_key),
        left_unit_kind =
            _white_lindsey_descriptor_property(pair_source, :left_unit_kind),
        right_unit_kind =
            _white_lindsey_descriptor_property(pair_source, :right_unit_kind),
        left_stratum_kind = left_coefficients.stratum_kind,
        right_stratum_kind = right_coefficients.stratum_kind,
        left_coefficient_status = left_coefficients.status,
        right_coefficient_status = right_coefficients.status,
        left_blocker = left_coefficients.blocker,
        right_blocker = right_coefficients.blocker,
        left_coefficient_space = left_coefficients.coefficient_space,
        right_coefficient_space = right_coefficients.coefficient_space,
        left_coefficient_matrix =
            left_materialized ? left_coefficients.coefficient_matrix : nothing,
        right_coefficient_matrix =
            right_materialized ? right_coefficients.coefficient_matrix : nothing,
        left_support_indices =
            left_materialized ? left_coefficients.support_indices : nothing,
        right_support_indices =
            right_materialized ? right_coefficients.support_indices : nothing,
        left_retained_column_count =
            left_materialized ? left_coefficients.retained_column_count : nothing,
        right_retained_column_count =
            right_materialized ? right_coefficients.retained_column_count : nothing,
        left_source_support_row_count =
            left_materialized ? left_coefficients.source_support_row_count : nothing,
        right_source_support_row_count =
            right_materialized ? right_coefficients.source_support_row_count : nothing,
        left_nonzero_count =
            left_materialized ? left_coefficients.nonzero_count : 0,
        right_nonzero_count =
            right_materialized ? right_coefficients.nonzero_count : 0,
        unit_coefficient_cache_scope = :local_pair_or_batch_call,
        unit_coefficient_cache_entry_count = cache_entry_count,
        coefficient_maps_materialized =
            left_materialized && right_materialized,
        pair_unit_coefficient_maps_materialized =
            left_materialized && right_materialized,
        pair_blocks_materialized = false,
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

function _white_lindsey_pair_unit_coefficients_batch_result(
    results,
    cache_entry_count::Int,
)
    materialized_count =
        count(result -> result.pair_unit_coefficient_maps_materialized, results)
    blocked_count = length(results) - materialized_count
    status, blocker =
        blocked_count == 0 ?
        (:materialized_white_lindsey_pair_unit_coefficients_batch, nothing) :
        (
            :blocked_white_lindsey_pair_unit_coefficients_batch,
            first(
                result.blocker for result in results
                if !result.pair_unit_coefficient_maps_materialized
            ),
        )
    return (;
        object_kind =
            :white_lindsey_boundary_stratum_pair_unit_coefficients_batch,
        status,
        blocker,
        record_count = length(results),
        materialized_count,
        blocked_count,
        results,
        unit_coefficient_cache_scope = :local_pair_batch_call,
        unit_coefficient_cache_entry_count = cache_entry_count,
        coefficient_maps_materialized = blocked_count == 0,
        pair_blocks_materialized = false,
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

function _white_lindsey_pair_unit_coefficients_status(
    left_materialized::Bool,
    right_materialized::Bool,
)
    !left_materialized && return (
        :blocked_white_lindsey_pair_unit_coefficients,
        :left_white_lindsey_unit_coefficients_not_materialized,
    )
    !right_materialized && return (
        :blocked_white_lindsey_pair_unit_coefficients,
        :right_white_lindsey_unit_coefficients_not_materialized,
    )
    return :materialized_white_lindsey_pair_unit_coefficients, nothing
end

function _white_lindsey_unit_coefficients_materialized(coefficients)
    return _white_lindsey_descriptor_property(
        coefficients,
        :coefficient_maps_materialized,
        false,
    ) === true
end

function _white_lindsey_explicit_pair_key(left_unit, right_unit)
    return (
        _white_lindsey_descriptor_property(left_unit, :unit_key, :left_unit),
        _white_lindsey_descriptor_property(right_unit, :unit_key, :right_unit),
    )
end

function _white_lindsey_explicit_pair_family(left_unit, right_unit)
    left_kind = _white_lindsey_descriptor_property(left_unit, :unit_kind)
    right_kind = _white_lindsey_descriptor_property(right_unit, :unit_kind)
    (isnothing(left_kind) || isnothing(right_kind)) && return nothing
    return Symbol(String(left_kind), "__", String(right_kind))
end
