"""
    SpherePointSetProvenance

Minimal provenance record for one curated sphere-point-set source.

This is an experimental research-track data container, not a frozen public
format contract.
"""
struct SpherePointSetProvenance
    source_tag::String
    source_project::String
    source_artifact::String
    source_note::String
end

function Base.show(io::IO, provenance::SpherePointSetProvenance)
    print(
        io,
        "SpherePointSetProvenance(source_tag=",
        repr(provenance.source_tag),
        ", source_artifact=",
        repr(provenance.source_artifact),
        ")",
    )
end

"""
    CuratedSpherePointSet

Read-only curated sphere-point-set entry for the experimental angular research
track.

The current repo-side scaffold intentionally stops at point-set access and
provenance. It does not yet build shell-local injected angular bases.
"""
struct CuratedSpherePointSet
    order::Int
    cardinality::Int
    coordinates::Matrix{Float64}
    nn_ratio::Float64
    provenance::SpherePointSetProvenance

    function CuratedSpherePointSet(
        order::Int,
        cardinality::Int,
        coordinates::AbstractMatrix{<:Real},
        nn_ratio::Real,
        provenance::SpherePointSetProvenance,
    )
        order >= 1 || throw(ArgumentError("CuratedSpherePointSet requires order >= 1"))
        cardinality >= 1 || throw(ArgumentError("CuratedSpherePointSet requires cardinality >= 1"))
        size(coordinates, 1) == cardinality ||
            throw(ArgumentError("coordinate row count must match the recorded cardinality"))
        size(coordinates, 2) == 3 ||
            throw(ArgumentError("sphere-point coordinates must be stored as an N x 3 matrix"))
        return new(
            order,
            cardinality,
            Matrix{Float64}(coordinates),
            Float64(nn_ratio),
            provenance,
        )
    end
end

function Base.show(io::IO, set::CuratedSpherePointSet)
    print(
        io,
        "CuratedSpherePointSet(order=",
        set.order,
        ", cardinality=",
        set.cardinality,
        ", source_tag=",
        repr(set.provenance.source_tag),
        ")",
    )
end

_curated_sphere_point_sets_path() =
    normpath(joinpath(@__DIR__, "..", "data", "angular", "curated_sphere_points.toml"))

const _CURATED_SPHERE_POINT_SET_CACHE =
    Ref{Union{Nothing,NamedTuple{(:orders, :sets),Tuple{Vector{Int},Dict{Int,CuratedSpherePointSet}}}}}(nothing)

function _row_coordinates_matrix(rows)
    nrows = length(rows)
    coordinates = Matrix{Float64}(undef, nrows, 3)
    for i in 1:nrows
        row = rows[i]
        length(row) == 3 || throw(ArgumentError("each curated sphere-point row must have length 3"))
        for j in 1:3
            coordinates[i, j] = Float64(row[j])
        end
    end
    return coordinates
end

function _load_curated_sphere_point_sets()
    raw = TOML.parsefile(_curated_sphere_point_sets_path())
    meta = raw["meta"]
    provenance = SpherePointSetProvenance(
        String(meta["source_tag"]),
        String(meta["source_project"]),
        String(meta["source_artifact"]),
        String(meta["source_note"]),
    )

    sets = Dict{Int,CuratedSpherePointSet}()
    for (key, entry) in raw["sets"]
        order = parse(Int, key)
        set = CuratedSpherePointSet(
            Int(entry["order"]),
            Int(entry["cardinality"]),
            _row_coordinates_matrix(entry["coordinates"]),
            Float64(entry["nn_ratio"]),
            provenance,
        )
        set.order == order || throw(ArgumentError("point-set key/order mismatch for entry $(repr(key))"))
        sets[order] = set
    end

    orders = sort!(Int[Int(value) for value in meta["orders"]])
    Set(orders) == Set(keys(sets)) ||
        throw(ArgumentError("curated sphere-point-set metadata orders do not match the stored entries"))

    return (orders = orders, sets = sets)
end

function _curated_sphere_point_set_cache()
    cache = _CURATED_SPHERE_POINT_SET_CACHE[]
    if cache === nothing
        cache = _load_curated_sphere_point_sets()
        _CURATED_SPHERE_POINT_SET_CACHE[] = cache
    end
    return cache
end

"""
    curated_sphere_point_set_orders()

Return the available curated sphere-point-set orders currently vendored in the
experimental angular research-track scaffold.
"""
function curated_sphere_point_set_orders()
    return copy(_curated_sphere_point_set_cache().orders)
end

"""
    curated_sphere_point_set(order::Int)

Load one read-only curated sphere-point set by order/cardinality.

This helper is intentionally narrow: it exposes only the curated subset
vendored with `GaussletBases` for the active experimental angular line.
"""
function curated_sphere_point_set(order::Int)
    cache = _curated_sphere_point_set_cache()
    set = get(cache.sets, order, nothing)
    set === nothing && throw(
        ArgumentError(
            "no curated sphere-point set is vendored for order $order; available orders are $(join(cache.orders, ", "))",
        ),
    )
    return set
end
