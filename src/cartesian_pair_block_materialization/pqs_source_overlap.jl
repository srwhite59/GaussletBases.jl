# First tiny PQS/PQS raw source-space overlap pilot.

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
    _assert_pqs_source_pair_overlap_record(record)

    left_dims = _pqs_source_mode_dims(record, :left)
    right_dims = _pqs_source_mode_dims(record, :right)
    left_count = _pqs_source_mode_count(record, :left, left_dims)
    right_count = _pqs_source_mode_count(record, :right, right_dims)
    ordering = _pqs_source_mode_ordering(record)
    overlap_x, overlap_y, overlap_z = _overlap_1d_tuple(overlap_1d)
    overlap_axes = (overlap_x, overlap_y, overlap_z)
    _assert_pqs_source_overlap_axis_sizes(overlap_axes, left_dims, right_dims)

    left_modes = CRPS.source_mode_indices(left_dims; source_mode_ordering = ordering)
    right_modes = CRPS.source_mode_indices(right_dims; source_mode_ordering = ordering)
    length(left_modes) == left_count ||
        throw(ArgumentError("left source-mode count does not match left source-mode dims"))
    length(right_modes) == right_count ||
        throw(ArgumentError("right source-mode count does not match right source-mode dims"))

    block = Matrix{Float64}(undef, left_count, right_count)
    _fill_direct_direct_product_block!(
        block,
        left_modes,
        right_modes,
        overlap_axes[1],
        overlap_axes[2],
        overlap_axes[3],
    )

    return PairBlockMaterializationResult(
        :source_overlap,
        record.pair_key,
        block,
        true,
        true,
        false,
        false,
        false,
        false,
        (;
            materialization_path = record.materialization_path,
            readiness_status_before_materialization = record.readiness_status,
            block_space = :raw_product_source_modes,
            source_mode_ordering = ordering,
            left_source_mode_dims = left_dims,
            right_source_mode_dims = right_dims,
            left_source_mode_count = left_count,
            right_source_mode_count = right_count,
            source_operator_blocks_materialized = true,
            final_pair_blocks_materialized = false,
            operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
            artifacts_materialized = false,
        ),
    )
end

function _assert_pqs_source_pair_overlap_record(record::PairBlockMaterializationRecord)
    record.materialization_path === :pqs_source_pair_preflight ||
        throw(ArgumentError("PQS source overlap requires a PQS source-pair preflight record"))
    record.readiness_status === :ready_metadata_only_not_materialized ||
        throw(
            ArgumentError(
                "PQS source overlap requires a ready metadata-only PQS source-pair preflight record",
            ),
        )
    isnothing(record.blocker) ||
        throw(ArgumentError("PQS source overlap record is blocked"))
    _pqs_source_mode_ordering(record)
    _pqs_source_mode_dims(record, :left)
    _pqs_source_mode_dims(record, :right)
    return nothing
end

function _pqs_source_mode_dims(record::PairBlockMaterializationRecord, side::Symbol)
    key =
        side === :left ? :left_source_mode_dims :
        side === :right ? :right_source_mode_dims :
        throw(ArgumentError("PQS source side must be :left or :right"))
    haskey(record.metadata, key) ||
        throw(ArgumentError("PQS source overlap record is missing $(key)"))
    dims = getfield(record.metadata, key)
    isnothing(dims) &&
        throw(ArgumentError("PQS source overlap record has unavailable $(key)"))
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
        throw(ArgumentError("PQS source overlap record is missing $(key)"))
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
        throw(ArgumentError("PQS source overlap record is missing source_mode_ordering"))
    ordering = record.metadata.source_mode_ordering
    ordering isa Symbol ||
        throw(ArgumentError("PQS source overlap record has unavailable source_mode_ordering"))
    left_ordering =
        haskey(record.metadata, :left_source_mode_ordering) ?
        record.metadata.left_source_mode_ordering :
        ordering
    right_ordering =
        haskey(record.metadata, :right_source_mode_ordering) ?
        record.metadata.right_source_mode_ordering :
        ordering
    left_ordering === ordering && right_ordering === ordering ||
        throw(ArgumentError("PQS source overlap requires compatible source-mode ordering"))
    return ordering
end

function _assert_pqs_source_overlap_axis_sizes(
    overlap_axes,
    left_dims::NTuple{3,Int},
    right_dims::NTuple{3,Int},
)
    axis_names = (:x, :y, :z)
    for axis_index in 1:3
        overlap = overlap_axes[axis_index]
        overlap isa AbstractMatrix{<:Real} ||
            throw(ArgumentError("overlap_1d entries must be real matrices"))
        size(overlap) == (left_dims[axis_index], right_dims[axis_index]) ||
            throw(
                ArgumentError(
                    "overlap_1d.$(axis_names[axis_index]) has incompatible source-mode size",
                ),
            )
    end
    return nothing
end
