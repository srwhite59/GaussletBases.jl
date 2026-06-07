# Metadata-only axis transform fact defaults.
#
# These facts document what a raw product source would need on each axis. They
# do not construct source-mode coefficient matrices. Numerical axis transforms
# remain owned by legacy adapters or future materialization modules.

function _default_axis_transform_fact(
    axis::Int,
    source_interval::UnitRange{Int},
    source_mode_dim::Int,
)
    return AxisSourceTransformFact(
        axis,
        source_interval,
        source_mode_dim,
        :not_materialized,
        nothing,
        (; object_kind = :axis_source_transform_fact, materialized = false),
    )
end

function _default_axis_transform_facts(
    source_intervals::NTuple{3,UnitRange{Int}},
    source_mode_dims::NTuple{3,Int},
)
    return ntuple(
        axis -> _default_axis_transform_fact(
            axis,
            source_intervals[axis],
            source_mode_dims[axis],
        ),
        3,
    )
end
