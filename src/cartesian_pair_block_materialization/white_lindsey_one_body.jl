# Selector over White--Lindsey one-body pair-block pilots.

const _WHITE_LINDSEY_ONE_BODY_TERMS = (
    :overlap,
    :position_x,
    :position_y,
    :position_z,
    :x2_x,
    :x2_y,
    :x2_z,
    :kinetic,
)

"""
    white_lindsey_boundary_stratum_one_body_block(pair_unit_coefficients, term; ...)
    white_lindsey_boundary_stratum_one_body_block(unit_pair, term; ...)

Dispatch one supported White--Lindsey boundary-stratum one-body safe-term
pilot. Supported terms are overlap, position_x/y/z, x2_x/y/z, and kinetic.
This selector does not build plan-level operators, Hamiltonian data, exports,
artifacts, IDA/MWG data, or Coulomb.
"""
function white_lindsey_boundary_stratum_one_body_block(
    pair_unit_coefficients,
    term;
    parent_axis_counts,
    overlap_1d = nothing,
    position_1d = nothing,
    x2_1d = nothing,
    kinetic_1d = nothing,
)
    term = _white_lindsey_one_body_term(term)
    if term === :overlap
        _white_lindsey_require_one_body_factor(term, :overlap_1d, overlap_1d)
        return white_lindsey_boundary_stratum_overlap_block(
            pair_unit_coefficients;
            parent_axis_counts,
            overlap_1d,
        )
    elseif term in (:position_x, :position_y, :position_z)
        _white_lindsey_require_one_body_factor(term, :overlap_1d, overlap_1d)
        _white_lindsey_require_one_body_factor(term, :position_1d, position_1d)
        return white_lindsey_boundary_stratum_position_block(
            pair_unit_coefficients;
            axis = _white_lindsey_term_axis(term, "position_"),
            parent_axis_counts,
            overlap_1d,
            position_1d,
        )
    elseif term in (:x2_x, :x2_y, :x2_z)
        _white_lindsey_require_one_body_factor(term, :overlap_1d, overlap_1d)
        _white_lindsey_require_one_body_factor(term, :x2_1d, x2_1d)
        return white_lindsey_boundary_stratum_x2_block(
            pair_unit_coefficients;
            axis = _white_lindsey_term_axis(term, "x2_"),
            parent_axis_counts,
            overlap_1d,
            x2_1d,
        )
    elseif term === :kinetic
        _white_lindsey_require_one_body_factor(term, :overlap_1d, overlap_1d)
        _white_lindsey_require_one_body_factor(term, :kinetic_1d, kinetic_1d)
        return white_lindsey_boundary_stratum_kinetic_block(
            pair_unit_coefficients;
            parent_axis_counts,
            overlap_1d,
            kinetic_1d,
        )
    end
    throw(ArgumentError("unsupported White--Lindsey one-body term $(term)"))
end

function white_lindsey_boundary_stratum_one_body_block(
    unit_pair::CUP.UnitPairRecord,
    term;
    parent_axis_counts,
    overlap_1d = nothing,
    position_1d = nothing,
    x2_1d = nothing,
    kinetic_1d = nothing,
)
    return white_lindsey_boundary_stratum_one_body_block(
        white_lindsey_boundary_stratum_pair_unit_coefficients(unit_pair),
        term;
        parent_axis_counts,
        overlap_1d,
        position_1d,
        x2_1d,
        kinetic_1d,
    )
end

"""
    white_lindsey_boundary_stratum_one_body_blocks(records_or_pairs, term; ...)

Materialize one supported White--Lindsey boundary-stratum one-body term over a
local collection of prepared pair-unit coefficient records or `UnitPairRecord`s.
Unsupported or blocked pair-unit coefficient inputs are returned as compact
skipped summaries. Missing factor groups for the requested term throw
`ArgumentError`.
"""
function white_lindsey_boundary_stratum_one_body_blocks(
    records_or_pairs,
    term;
    parent_axis_counts,
    overlap_1d = nothing,
    position_1d = nothing,
    x2_1d = nothing,
    kinetic_1d = nothing,
)
    term = _white_lindsey_one_body_term(term)
    _white_lindsey_require_one_body_factors(
        term;
        overlap_1d,
        position_1d,
        x2_1d,
        kinetic_1d,
    )
    pair_unit_coefficients, input_kind, cache_entry_count =
        _white_lindsey_one_body_batch_pair_inputs(records_or_pairs)

    results = PairBlockMaterializationResult[]
    skipped = NamedTuple[]
    for pair_unit_coefficient in pair_unit_coefficients
        if _is_ready_white_lindsey_pair_unit_coefficients(pair_unit_coefficient)
            push!(
                results,
                white_lindsey_boundary_stratum_one_body_block(
                    pair_unit_coefficient,
                    term;
                    parent_axis_counts,
                    overlap_1d,
                    position_1d,
                    x2_1d,
                    kinetic_1d,
                ),
            )
        else
            push!(
                skipped,
                _skipped_white_lindsey_pair_unit_coefficients_summary(
                    pair_unit_coefficient,
                ),
            )
        end
    end

    result_tuple = Tuple(results)
    skipped_tuple = Tuple(skipped)
    any_materialized = !isempty(result_tuple)
    return PairBlockMaterializationBatchResult(
        term,
        result_tuple,
        skipped_tuple,
        length(result_tuple),
        length(skipped_tuple),
        any_materialized,
        any_materialized,
        any_materialized,
        false,
        false,
        false,
        (;
            materialization_path =
                :white_lindsey_boundary_stratum_one_body_batch_selector,
            selector_helper = :white_lindsey_boundary_stratum_one_body_block,
            pair_input_kind = input_kind,
            pair_unit_coefficient_record_count = length(pair_unit_coefficients),
            unit_coefficient_cache_entry_count = cache_entry_count,
            local_pair_block_materialized = any_materialized,
            source_operator_blocks_materialized = any_materialized,
            final_pair_blocks_materialized = any_materialized,
            operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
            artifacts_materialized = false,
            ida_mwg_data_materialized = false,
            dense_parent_fallback_used = false,
            old_white_lindsey_seed_route_authority = false,
        ),
    )
end

function _white_lindsey_one_body_term(term)
    term isa Symbol || throw(
        ArgumentError("White--Lindsey one-body term must be a Symbol"),
    )
    term in _WHITE_LINDSEY_ONE_BODY_TERMS || throw(
        ArgumentError("unsupported White--Lindsey one-body term $(term)"),
    )
    return term
end

function _white_lindsey_require_one_body_factors(
    term::Symbol;
    overlap_1d,
    position_1d,
    x2_1d,
    kinetic_1d,
)
    if term === :overlap
        _white_lindsey_require_one_body_factor(term, :overlap_1d, overlap_1d)
    elseif term in (:position_x, :position_y, :position_z)
        _white_lindsey_require_one_body_factor(term, :overlap_1d, overlap_1d)
        _white_lindsey_require_one_body_factor(term, :position_1d, position_1d)
    elseif term in (:x2_x, :x2_y, :x2_z)
        _white_lindsey_require_one_body_factor(term, :overlap_1d, overlap_1d)
        _white_lindsey_require_one_body_factor(term, :x2_1d, x2_1d)
    elseif term === :kinetic
        _white_lindsey_require_one_body_factor(term, :overlap_1d, overlap_1d)
        _white_lindsey_require_one_body_factor(term, :kinetic_1d, kinetic_1d)
    else
        throw(ArgumentError("unsupported White--Lindsey one-body term $(term)"))
    end
    return nothing
end

function _white_lindsey_require_one_body_factor(term::Symbol, name::Symbol, value)
    isnothing(value) && throw(
        ArgumentError("White--Lindsey $(term) requires $(name)"),
    )
    return nothing
end

function _white_lindsey_one_body_batch_pair_inputs(records_or_pairs)
    if _white_lindsey_descriptor_property(records_or_pairs, :object_kind) ===
       :white_lindsey_boundary_stratum_pair_unit_coefficients_batch
        return (
            records_or_pairs.results,
            :white_lindsey_pair_unit_coefficients_batch,
            _white_lindsey_descriptor_property(
                records_or_pairs,
                :unit_coefficient_cache_entry_count,
            ),
        )
    end

    (records_or_pairs isa Tuple || records_or_pairs isa AbstractVector) ||
        throw(ArgumentError("White--Lindsey one-body batch expects a tuple or vector"))
    input_tuple = Tuple(records_or_pairs)
    if all(input -> input isa CUP.UnitPairRecord, input_tuple)
        batch =
            white_lindsey_boundary_stratum_pair_unit_coefficients(input_tuple)
        return (
            batch.results,
            :unit_pair_records,
            batch.unit_coefficient_cache_entry_count,
        )
    end
    return (input_tuple, :prepared_pair_unit_coefficients, nothing)
end

function _is_ready_white_lindsey_pair_unit_coefficients(pair_unit_coefficients)
    return _white_lindsey_descriptor_property(pair_unit_coefficients, :object_kind) ===
           :white_lindsey_boundary_stratum_pair_unit_coefficients &&
           _white_lindsey_descriptor_property(pair_unit_coefficients, :status) ===
           :materialized_white_lindsey_pair_unit_coefficients &&
           _white_lindsey_descriptor_property(
               pair_unit_coefficients,
               :pair_unit_coefficient_maps_materialized,
               false,
           ) === true
end

function _skipped_white_lindsey_pair_unit_coefficients_summary(
    pair_unit_coefficients,
)
    blocker = _white_lindsey_descriptor_property(pair_unit_coefficients, :blocker)
    if isnothing(blocker)
        blocker =
            _white_lindsey_descriptor_property(
                pair_unit_coefficients,
                :object_kind,
            ) === :white_lindsey_boundary_stratum_pair_unit_coefficients ?
            :white_lindsey_pair_unit_coefficients_not_materialized :
            :unsupported_white_lindsey_pair_unit_coefficients_input
    end
    return (;
        pair_key = _white_lindsey_descriptor_property(
            pair_unit_coefficients,
            :pair_key,
        ),
        pair_index = _white_lindsey_descriptor_property(
            pair_unit_coefficients,
            :pair_index,
        ),
        pair_family = _white_lindsey_descriptor_property(
            pair_unit_coefficients,
            :pair_family,
        ),
        object_kind = _white_lindsey_descriptor_property(
            pair_unit_coefficients,
            :object_kind,
        ),
        status = _white_lindsey_descriptor_property(
            pair_unit_coefficients,
            :status,
            :unsupported_white_lindsey_pair_unit_coefficients_input,
        ),
        blocker,
        left_coefficient_status = _white_lindsey_descriptor_property(
            pair_unit_coefficients,
            :left_coefficient_status,
        ),
        right_coefficient_status = _white_lindsey_descriptor_property(
            pair_unit_coefficients,
            :right_coefficient_status,
        ),
        coefficient_maps_materialized = _white_lindsey_descriptor_property(
            pair_unit_coefficients,
            :coefficient_maps_materialized,
            false,
        ),
        pair_unit_coefficient_maps_materialized =
            _white_lindsey_descriptor_property(
                pair_unit_coefficients,
                :pair_unit_coefficient_maps_materialized,
                false,
            ),
    )
end

function _white_lindsey_term_axis(term::Symbol, prefix::AbstractString)
    text = String(term)
    startswith(text, prefix) || throw(
        ArgumentError("White--Lindsey term $(term) does not have prefix $(prefix)"),
    )
    axis = Symbol(text[(lastindex(prefix) + 1):end])
    axis in (:x, :y, :z) || throw(
        ArgumentError("White--Lindsey term $(term) does not name axis x, y, or z"),
    )
    return axis
end
