# PQS/PQS raw source-space safe-term helpers.

"""
    pqs_source_pair_overlap_block(record; overlap_1d)

Materialize one PQS/PQS raw source-space overlap block from explicit 1D overlap
matrices. This consumes only a ready PQS source-pair preflight record and does
not build shell projection, Lowdin realization, final retained pair blocks,
Hamiltonian data, IDA data, exports, or artifacts.
"""
function pqs_source_pair_overlap_block(
    record::PairBlockMaterializationRecord;
    overlap_1d,
)
    descriptor = _supported_pqs_source_safe_term_descriptor(:overlap)
    overlap_x, overlap_y, overlap_z = _overlap_1d_tuple(overlap_1d)
    overlap_axes = (overlap_x, overlap_y, overlap_z)
    left_dims, right_dims = _pqs_source_pair_dims(record)
    _assert_pqs_source_axis_sizes(overlap_axes, left_dims, right_dims, "overlap_1d")
    return _pqs_source_pair_product_result(
        record,
        descriptor.source_term,
        overlap_axes,
        (;),
    )
end

"""
    pqs_source_pair_overlap_blocks(plan; overlap_1d)

Materialize raw source-space overlap blocks only for ready PQS/PQS source-pair
records in a pair-block materialization plan. Unsupported or blocked records
are returned as compact skipped summaries.
"""
function pqs_source_pair_overlap_blocks(
    plan::PairBlockMaterializationPlan;
    overlap_1d,
)
    descriptor = _supported_pqs_source_safe_term_descriptor(:overlap)
    return _pqs_source_pair_batch_results(
        record -> pqs_source_pair_overlap_block(record; overlap_1d),
        plan,
        descriptor.source_term,
        descriptor.batch_materialization_path,
        descriptor.unsupported_record_blocker,
        (;),
    )
end

"""
    pqs_source_pair_position_block(record; axis, overlap_1d, position_1d)

Materialize one PQS/PQS raw source-space position block for
`axis in (:x, :y, :z)` from explicit 1D source-mode factors. This remains
before shell projection, Lowdin realization, and final retained pair blocks.
"""
function pqs_source_pair_position_block(
    record::PairBlockMaterializationRecord;
    axis,
    overlap_1d,
    position_1d,
)
    descriptor = _pqs_source_axis_safe_term_descriptor(:position, axis)
    left_dims, right_dims = _pqs_source_pair_dims(record)
    overlap_x, overlap_y, overlap_z = _overlap_1d_tuple(overlap_1d)
    position_x, position_y, position_z =
        _operator_1d_tuple(position_1d, "position_1d")
    overlap_axes = (overlap_x, overlap_y, overlap_z)
    position_axes = (position_x, position_y, position_z)
    _assert_pqs_source_axis_sizes(overlap_axes, left_dims, right_dims, "overlap_1d")
    _assert_pqs_source_axis_sizes(position_axes, left_dims, right_dims, "position_1d")

    operator_axes =
        axis === :x ? (position_x, overlap_y, overlap_z) :
        axis === :y ? (overlap_x, position_y, overlap_z) :
        (overlap_x, overlap_y, position_z)
    return _pqs_source_pair_product_result(
        record,
        descriptor.source_term,
        operator_axes,
        _pqs_source_axis_metadata(descriptor),
    )
end

"""
    pqs_source_pair_position_blocks(plan; axis, overlap_1d, position_1d)

Materialize raw source-space position blocks only for ready PQS/PQS source-pair
records in a pair-block materialization plan.
"""
function pqs_source_pair_position_blocks(
    plan::PairBlockMaterializationPlan;
    axis,
    overlap_1d,
    position_1d,
)
    descriptor = _pqs_source_axis_safe_term_descriptor(:position, axis)
    return _pqs_source_pair_batch_results(
        record -> pqs_source_pair_position_block(
            record;
            axis,
            overlap_1d,
            position_1d,
        ),
        plan,
        descriptor.source_term,
        descriptor.batch_materialization_path,
        descriptor.unsupported_record_blocker,
        _pqs_source_axis_metadata(descriptor),
    )
end

"""
    pqs_source_pair_x2_block(record; axis, overlap_1d, x2_1d)

Materialize one PQS/PQS raw source-space x2 block for
`axis in (:x, :y, :z)` from explicit 1D source-mode factors. This remains
before shell projection, Lowdin realization, and final retained pair blocks.
"""
function pqs_source_pair_x2_block(
    record::PairBlockMaterializationRecord;
    axis,
    overlap_1d,
    x2_1d,
)
    descriptor = _pqs_source_axis_safe_term_descriptor(:x2, axis)
    left_dims, right_dims = _pqs_source_pair_dims(record)
    overlap_x, overlap_y, overlap_z = _overlap_1d_tuple(overlap_1d)
    x2_x, x2_y, x2_z = _operator_1d_tuple(x2_1d, "x2_1d")
    overlap_axes = (overlap_x, overlap_y, overlap_z)
    x2_axes = (x2_x, x2_y, x2_z)
    _assert_pqs_source_axis_sizes(overlap_axes, left_dims, right_dims, "overlap_1d")
    _assert_pqs_source_axis_sizes(x2_axes, left_dims, right_dims, "x2_1d")

    operator_axes =
        axis === :x ? (x2_x, overlap_y, overlap_z) :
        axis === :y ? (overlap_x, x2_y, overlap_z) :
        (overlap_x, overlap_y, x2_z)
    return _pqs_source_pair_product_result(
        record,
        descriptor.source_term,
        operator_axes,
        _pqs_source_axis_metadata(descriptor),
    )
end

"""
    pqs_source_pair_x2_blocks(plan; axis, overlap_1d, x2_1d)

Materialize raw source-space x2 blocks only for ready PQS/PQS source-pair
records in a pair-block materialization plan.
"""
function pqs_source_pair_x2_blocks(
    plan::PairBlockMaterializationPlan;
    axis,
    overlap_1d,
    x2_1d,
)
    descriptor = _pqs_source_axis_safe_term_descriptor(:x2, axis)
    return _pqs_source_pair_batch_results(
        record -> pqs_source_pair_x2_block(
            record;
            axis,
            overlap_1d,
            x2_1d,
        ),
        plan,
        descriptor.source_term,
        descriptor.batch_materialization_path,
        descriptor.unsupported_record_blocker,
        _pqs_source_axis_metadata(descriptor),
    )
end

"""
    pqs_source_pair_kinetic_block(record; overlap_1d, kinetic_1d)

Materialize one PQS/PQS raw source-space kinetic block from explicit 1D
source-mode factors. The caller-supplied 1D kinetic factors own signs and
prefactors. This remains before shell projection, Lowdin realization, and final
retained pair blocks.
"""
function pqs_source_pair_kinetic_block(
    record::PairBlockMaterializationRecord;
    overlap_1d,
    kinetic_1d,
)
    descriptor = _supported_pqs_source_safe_term_descriptor(:kinetic)
    left_dims, right_dims = _pqs_source_pair_dims(record)
    overlap_x, overlap_y, overlap_z = _overlap_1d_tuple(overlap_1d)
    kinetic_x, kinetic_y, kinetic_z =
        _operator_1d_tuple(kinetic_1d, "kinetic_1d")
    overlap_axes = (overlap_x, overlap_y, overlap_z)
    kinetic_axes = (kinetic_x, kinetic_y, kinetic_z)
    _assert_pqs_source_axis_sizes(overlap_axes, left_dims, right_dims, "overlap_1d")
    _assert_pqs_source_axis_sizes(kinetic_axes, left_dims, right_dims, "kinetic_1d")

    kinetic_factor_form = _pqs_source_kinetic_factor_form()
    kinetic_x_result = _pqs_source_pair_product_result(
        record,
        descriptor.source_term,
        (kinetic_x, overlap_y, overlap_z),
        (; kinetic_factor_form),
    )
    kinetic_y_result = _pqs_source_pair_product_result(
        record,
        descriptor.source_term,
        (overlap_x, kinetic_y, overlap_z),
        (; kinetic_factor_form),
    )
    kinetic_z_result = _pqs_source_pair_product_result(
        record,
        descriptor.source_term,
        (overlap_x, overlap_y, kinetic_z),
        (; kinetic_factor_form),
    )

    return PairBlockMaterializationResult(
        descriptor.source_term,
        record.pair_key,
        kinetic_x_result.block + kinetic_y_result.block + kinetic_z_result.block,
        true,
        true,
        false,
        false,
        false,
        false,
        kinetic_x_result.metadata,
    )
end

"""
    pqs_source_pair_kinetic_blocks(plan; overlap_1d, kinetic_1d)

Materialize raw source-space kinetic blocks only for ready PQS/PQS source-pair
records in a pair-block materialization plan.
"""
function pqs_source_pair_kinetic_blocks(
    plan::PairBlockMaterializationPlan;
    overlap_1d,
    kinetic_1d,
)
    descriptor = _supported_pqs_source_safe_term_descriptor(:kinetic)
    return _pqs_source_pair_batch_results(
        record -> pqs_source_pair_kinetic_block(record; overlap_1d, kinetic_1d),
        plan,
        descriptor.source_term,
        descriptor.batch_materialization_path,
        descriptor.unsupported_record_blocker,
        (; kinetic_factor_form = _pqs_source_kinetic_factor_form()),
    )
end

function _source_position_term(axis)
    return _pqs_source_axis_safe_term_descriptor(:position, axis).source_term
end

function _source_x2_term(axis)
    return _pqs_source_axis_safe_term_descriptor(:x2, axis).source_term
end

function _pqs_source_kinetic_factor_form()
    return (
        (:kinetic, :overlap, :overlap),
        (:overlap, :kinetic, :overlap),
        (:overlap, :overlap, :kinetic),
    )
end

function _pqs_source_safe_term_descriptor(term::Symbol)
    if term === :overlap
        return _pqs_source_safe_term_descriptor(
            term,
            :source_overlap,
            :overlap,
            nothing,
            "overlap_1d",
            :ready_pqs_source_overlap_blocks_only,
            :unsupported_pqs_source_overlap_materialization_record,
            nothing,
        )
    end

    position_axis = _position_axis_for_term(term)
    if !isnothing(position_axis)
        return _pqs_source_safe_term_descriptor(
            term,
            Symbol(:source_, term),
            :position,
            position_axis,
            "position_1d",
            :ready_pqs_source_position_blocks_only,
            :unsupported_pqs_source_position_materialization_record,
            :position_axis,
        )
    end

    x2_axis = _x2_axis_for_term(term)
    if !isnothing(x2_axis)
        return _pqs_source_safe_term_descriptor(
            term,
            Symbol(:source_, term),
            :x2,
            x2_axis,
            "x2_1d",
            :ready_pqs_source_x2_blocks_only,
            :unsupported_pqs_source_x2_materialization_record,
            :x2_axis,
        )
    end

    if term === :kinetic
        return _pqs_source_safe_term_descriptor(
            term,
            :source_kinetic,
            :kinetic,
            nothing,
            "kinetic_1d",
            :ready_pqs_source_kinetic_blocks_only,
            :unsupported_pqs_source_kinetic_materialization_record,
            nothing,
        )
    end

    return (;
        status = :blocked_unsupported_pqs_source_safe_term,
        blocker = :unsupported_pqs_source_one_body_term,
        requested_term = term,
        source_term = nothing,
        family = :unsupported,
        axis = nothing,
        required_factor_name = nothing,
        batch_materialization_path = :blocked_pqs_source_safe_term_blocks,
        unsupported_record_blocker = :unsupported_pqs_source_materialization_record,
        axis_metadata_key = nothing,
    )
end

function _pqs_source_safe_term_descriptor(
    requested_term::Symbol,
    source_term::Symbol,
    family::Symbol,
    axis,
    required_factor_name::AbstractString,
    batch_materialization_path::Symbol,
    unsupported_record_blocker::Symbol,
    axis_metadata_key,
)
    return (;
        status = :available_pqs_source_safe_term,
        blocker = nothing,
        requested_term,
        source_term,
        family,
        axis,
        required_factor_name = String(required_factor_name),
        batch_materialization_path,
        unsupported_record_blocker,
        axis_metadata_key,
    )
end

function _supported_pqs_source_safe_term_descriptor(term::Symbol)
    descriptor = _pqs_source_safe_term_descriptor(term)
    isnothing(descriptor.blocker) ||
        throw(ArgumentError("unsupported PQS source one-body term: $(term)"))
    return descriptor
end

function _pqs_source_axis_safe_term_descriptor(family::Symbol, axis)
    requested_term =
        family === :position ? _position_term(axis) :
        family === :x2 ? _x2_term(axis) :
        throw(ArgumentError("PQS source axis family must be :position or :x2"))
    return _supported_pqs_source_safe_term_descriptor(requested_term)
end

function _pqs_source_axis_metadata(descriptor)
    isnothing(descriptor.axis_metadata_key) && return (;)
    return (; descriptor.axis_metadata_key => descriptor.axis)
end

function _is_ready_pqs_source_pair_record(record::PairBlockMaterializationRecord)
    return record.materialization_path === :pqs_source_pair_preflight &&
           record.readiness_status === :ready_metadata_only_not_materialized &&
           isnothing(record.blocker)
end

function _skipped_pqs_source_record_summary(
    record::PairBlockMaterializationRecord,
    unsupported_blocker::Symbol,
)
    return (;
        pair_key = record.pair_key,
        pair_index = record.pair_index,
        pair_family = record.pair_family,
        materialization_path = record.materialization_path,
        readiness_status = record.readiness_status,
        blocker = isnothing(record.blocker) ? unsupported_blocker : record.blocker,
    )
end

function _pqs_source_pair_batch_results(
    materialize_record,
    plan::PairBlockMaterializationPlan,
    term::Symbol,
    materialization_path::Symbol,
    unsupported_blocker::Symbol,
    metadata,
)
    results = PairBlockMaterializationResult[]
    skipped = NamedTuple[]

    for record in pair_block_materialization_records(plan)
        if _is_ready_pqs_source_pair_record(record)
            push!(results, materialize_record(record))
        else
            push!(
                skipped,
                _skipped_pqs_source_record_summary(record, unsupported_blocker),
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
        false,
        false,
        false,
        false,
        merge(
            (; materialization_path),
            NamedTuple(metadata),
            (; pair_block_record_count = length(pair_block_materialization_records(plan))),
        ),
    )
end

function _pqs_source_pair_product_result(
    record::PairBlockMaterializationRecord,
    term::Symbol,
    operator_axes,
    metadata,
)
    _assert_pqs_source_pair_record(record)

    left_dims = _pqs_source_mode_dims(record, :left)
    right_dims = _pqs_source_mode_dims(record, :right)
    left_count = _pqs_source_mode_count(record, :left, left_dims)
    right_count = _pqs_source_mode_count(record, :right, right_dims)
    ordering = _pqs_source_mode_ordering(record)

    left_modes = CRPS.source_mode_indices(left_dims; source_mode_ordering = ordering)
    right_modes = CRPS.source_mode_indices(right_dims; source_mode_ordering = ordering)
    length(left_modes) == left_count ||
        throw(ArgumentError("left source-mode count does not match left source-mode dims"))
    length(right_modes) == right_count ||
        throw(ArgumentError("right source-mode count does not match right source-mode dims"))

    block = Matrix{Float64}(undef, left_count, right_count)
    _fill_source_mode_product_block!(
        block,
        left_modes,
        right_modes,
        operator_axes[1],
        operator_axes[2],
        operator_axes[3],
    )

    return PairBlockMaterializationResult(
        term,
        record.pair_key,
        block,
        true,
        true,
        false,
        false,
        false,
        false,
        merge(
            (;
                materialization_path = record.materialization_path,
                readiness_status_before_materialization = record.readiness_status,
                block_space = :raw_product_source_modes,
                source_mode_ordering = ordering,
                left_source_mode_dims = left_dims,
                right_source_mode_dims = right_dims,
                left_source_mode_count = left_count,
                right_source_mode_count = right_count,
                transform_contract_keys =
                    _pqs_source_pair_metadata_value(
                        record,
                        :transform_contract_keys,
                        (; left = nothing, right = nothing),
                    ),
                source_contract_keys =
                    _pqs_source_pair_metadata_value(
                        record,
                        :source_contract_keys,
                        (; left = nothing, right = nothing),
                    ),
                transform_paths = record.transform_path,
                realization_paths = record.realization_path,
                source_operator_blocks_materialized = true,
                final_pair_blocks_materialized = false,
                shell_realization_materialized = false,
                operator_blocks_materialized = false,
                hamiltonian_data_materialized = false,
                artifacts_materialized = false,
            ),
            NamedTuple(metadata),
        ),
    )
end

function _pqs_source_pair_metadata_value(
    record::PairBlockMaterializationRecord,
    key::Symbol,
    default = nothing,
)
    return haskey(record.metadata, key) ? getfield(record.metadata, key) : default
end

function _assert_pqs_source_pair_overlap_record(record::PairBlockMaterializationRecord)
    return _assert_pqs_source_pair_record(record)
end

function _assert_pqs_source_pair_record(record::PairBlockMaterializationRecord)
    record.materialization_path === :pqs_source_pair_preflight ||
        throw(ArgumentError("PQS source block requires a PQS source-pair preflight record"))
    record.readiness_status === :ready_metadata_only_not_materialized ||
        throw(
            ArgumentError(
                "PQS source block requires a ready metadata-only PQS source-pair preflight record",
            ),
        )
    isnothing(record.blocker) ||
        throw(ArgumentError("PQS source block record is blocked"))
    _pqs_source_mode_ordering(record)
    _pqs_source_mode_dims(record, :left)
    _pqs_source_mode_dims(record, :right)
    return nothing
end

function _pqs_source_pair_dims(record::PairBlockMaterializationRecord)
    _assert_pqs_source_pair_record(record)
    return (
        _pqs_source_mode_dims(record, :left),
        _pqs_source_mode_dims(record, :right),
    )
end

function _pqs_source_mode_dims(record::PairBlockMaterializationRecord, side::Symbol)
    key =
        side === :left ? :left_source_mode_dims :
        side === :right ? :right_source_mode_dims :
        throw(ArgumentError("PQS source side must be :left or :right"))
    haskey(record.metadata, key) ||
        throw(ArgumentError("PQS source block record is missing $(key)"))
    dims = getfield(record.metadata, key)
    isnothing(dims) &&
        throw(ArgumentError("PQS source block record has unavailable $(key)"))
    return _pqs_source_mode_dims_tuple(dims, String(key))
end

function _pqs_source_mode_dims_tuple(dims, name::AbstractString)
    length(dims) == 3 ||
        throw(ArgumentError("$(name) must contain exactly three dimensions"))
    normalized = map(dims) do dim
        dim isa Integer ||
            throw(ArgumentError("$(name) entries must be integers"))
        dim > 0 ||
            throw(ArgumentError("$(name) entries must be positive"))
        Int(dim)
    end
    return Tuple(normalized)::NTuple{3,Int}
end

function _pqs_source_mode_count(
    record::PairBlockMaterializationRecord,
    side::Symbol,
    dims::NTuple{3,Int},
)
    key =
        side === :left ? :left_source_mode_count :
        side === :right ? :right_source_mode_count :
        throw(ArgumentError("PQS source side must be :left or :right"))
    haskey(record.metadata, key) ||
        throw(ArgumentError("PQS source block record is missing $(key)"))
    count = getfield(record.metadata, key)
    count isa Integer ||
        throw(ArgumentError("$(key) must be an integer"))
    normalized_count = Int(count)
    normalized_count == prod(dims) ||
        throw(ArgumentError("$(key) does not match $(side)_source_mode_dims"))
    return normalized_count
end

function _pqs_source_mode_ordering(record::PairBlockMaterializationRecord)
    haskey(record.metadata, :source_mode_ordering) ||
        throw(ArgumentError("PQS source block record is missing source_mode_ordering"))
    ordering = record.metadata.source_mode_ordering
    ordering isa Symbol ||
        throw(ArgumentError("PQS source block record has unavailable source_mode_ordering"))
    left_ordering =
        haskey(record.metadata, :left_source_mode_ordering) ?
        record.metadata.left_source_mode_ordering :
        ordering
    right_ordering =
        haskey(record.metadata, :right_source_mode_ordering) ?
        record.metadata.right_source_mode_ordering :
        ordering
    left_ordering === ordering && right_ordering === ordering ||
        throw(ArgumentError("PQS source block requires compatible source-mode ordering"))
    return ordering
end

function _assert_pqs_source_axis_sizes(
    overlap_axes,
    left_dims::NTuple{3,Int},
    right_dims::NTuple{3,Int},
    name::AbstractString,
)
    axis_names = (:x, :y, :z)
    for axis_index in 1:3
        overlap = overlap_axes[axis_index]
        overlap isa AbstractMatrix{<:Real} ||
            throw(ArgumentError("$(name) entries must be real matrices"))
        size(overlap) == (left_dims[axis_index], right_dims[axis_index]) ||
            throw(
                ArgumentError(
                    "$(name).$(axis_names[axis_index]) has incompatible source-mode size",
                ),
            )
    end
    return nothing
end

function _assert_pqs_source_overlap_axis_sizes(
    overlap_axes,
    left_dims::NTuple{3,Int},
    right_dims::NTuple{3,Int},
)
    return _assert_pqs_source_axis_sizes(overlap_axes, left_dims, right_dims, "overlap_1d")
end
