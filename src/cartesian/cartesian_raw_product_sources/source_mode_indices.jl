# Deterministic source-mode ordering helpers.

const _SUPPORTED_SOURCE_MODE_ORDERING = :x_major_y_major_z_fast

function _assert_supported_source_mode_ordering(ordering::Symbol)
    ordering === _SUPPORTED_SOURCE_MODE_ORDERING ||
        throw(ArgumentError("unsupported source_mode_ordering: $(ordering)"))
    return nothing
end

"""
    source_mode_indices(source_mode_dims; source_mode_ordering = :x_major_y_major_z_fast)

Return deterministic source-mode tuples in the current x-major, z-fast order:
`(ix, iy, iz) for ix in 1:nx, iy in 1:ny, iz in 1:nz`.
"""
function source_mode_indices(
    dims;
    source_mode_ordering::Symbol = _SUPPORTED_SOURCE_MODE_ORDERING,
)
    normalized_dims = _normalize_source_mode_dims(dims)
    _assert_supported_source_mode_ordering(source_mode_ordering)
    nx, ny, nz = normalized_dims
    return NTuple{3,Int}[
        (ix, iy, iz)
        for ix in 1:nx
        for iy in 1:ny
        for iz in 1:nz
    ]
end

function source_mode_count(dims)
    return prod(_normalize_source_mode_dims(dims))
end
