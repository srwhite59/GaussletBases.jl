"""
    SpherePointSetProvenance

Minimal provenance record for one vendored sphere-point-set source.

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
    SpherePointSet

Read-only vendored sphere-point-set entry for the experimental angular
research track.

The normal/default angular backing store is the full vendored JLD2 point-set
collection under `data/angular/SpherePoints.jld2`. The smaller curated subset
remains available through explicit fixture helpers.
"""
struct SpherePointSet
    order::Int
    cardinality::Int
    coordinates::Matrix{Float64}
    nn_ratio::Float64
    provenance::SpherePointSetProvenance

    function SpherePointSet(
        order::Int,
        cardinality::Int,
        coordinates::AbstractMatrix{<:Real},
        nn_ratio::Real,
        provenance::SpherePointSetProvenance,
    )
        order >= 1 || throw(ArgumentError("SpherePointSet requires order >= 1"))
        cardinality >= 1 || throw(ArgumentError("SpherePointSet requires cardinality >= 1"))
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

"""
    CuratedSpherePointSet

Compatibility alias for the explicit curated fixture subset.
"""
const CuratedSpherePointSet = SpherePointSet

function Base.show(io::IO, set::SpherePointSet)
    print(
        io,
        "SpherePointSet(order=",
        set.order,
        ", cardinality=",
        set.cardinality,
        ", source_tag=",
        repr(set.provenance.source_tag),
        ")",
    )
end

_sphere_point_sets_path() =
    normpath(joinpath(@__DIR__, "..", "data", "angular", "SpherePoints.jld2"))
_sphere_point_sets_manifest_path() =
    normpath(joinpath(@__DIR__, "..", "data", "angular", "SpherePoints_manifest.toml"))
_curated_sphere_point_sets_path() =
    normpath(joinpath(@__DIR__, "..", "data", "angular", "curated_sphere_points.toml"))

const _SPHERE_POINT_SET_CACHE =
    Ref{Union{Nothing,NamedTuple{(:orders, :sets),Tuple{Vector{Int},Dict{Int,SpherePointSet}}}}}(nothing)
const _CURATED_SPHERE_POINT_SET_CACHE =
    Ref{Union{Nothing,NamedTuple{(:orders, :sets),Tuple{Vector{Int},Dict{Int,SpherePointSet}}}}}(nothing)

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

function _point_set_provenance(meta)
    return SpherePointSetProvenance(
        String(meta["source_tag"]),
        String(meta["source_project"]),
        String(meta["source_artifact"]),
        String(meta["source_note"]),
    )
end

function _load_sphere_point_sets()
    manifest = TOML.parsefile(_sphere_point_sets_manifest_path())
    meta = manifest["meta"]
    provenance = _point_set_provenance(meta)
    raw_sets = JLD2.load(_sphere_point_sets_path(), "xyzsets")

    sets = Dict{Int,SpherePointSet}()
    for order in eachindex(raw_sets)
        entry = raw_sets[order]
        entry === nothing && continue
        coordinates = entry[1]
        nn_ratio = entry[2]
        sets[order] = SpherePointSet(
            order,
            size(coordinates, 1),
            coordinates,
            nn_ratio,
            provenance,
        )
    end

    orders = sort!(collect(keys(sets)))
    manifest_orders = sort!(Int[Int(value) for value in meta["orders"]])
    Set(orders) == Set(manifest_orders) ||
        throw(ArgumentError("vendored sphere-point-set manifest orders do not match the JLD2 contents"))
    length(orders) == Int(meta["order_count"]) ||
        throw(ArgumentError("vendored sphere-point-set manifest order_count does not match the JLD2 contents"))
    first(orders) == Int(meta["min_order"]) ||
        throw(ArgumentError("vendored sphere-point-set manifest min_order does not match the JLD2 contents"))
    last(orders) == Int(meta["max_order"]) ||
        throw(ArgumentError("vendored sphere-point-set manifest max_order does not match the JLD2 contents"))
    return (orders = orders, sets = sets)
end

function _load_curated_sphere_point_sets()
    raw = TOML.parsefile(_curated_sphere_point_sets_path())
    meta = raw["meta"]
    provenance = _point_set_provenance(meta)

    sets = Dict{Int,SpherePointSet}()
    for (key, entry) in raw["sets"]
        order = parse(Int, key)
        set = SpherePointSet(
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

function _sphere_point_set_cache()
    cache = _SPHERE_POINT_SET_CACHE[]
    if cache === nothing
        cache = _load_sphere_point_sets()
        _SPHERE_POINT_SET_CACHE[] = cache
    end
    return cache
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
    sphere_point_set_orders()

Return the populated orders available in the full vendored sphere-point-set
collection.

This is the normal/default angular order pool.
"""
function sphere_point_set_orders()
    return copy(_sphere_point_set_cache().orders)
end

"""
    sphere_point_set(order::Int)

Load one read-only sphere point set from the full vendored JLD2 collection.

This is the primary point-set access path for normal angular use.
"""
function sphere_point_set(order::Int)
    cache = _sphere_point_set_cache()
    set = get(cache.sets, order, nothing)
    set === nothing && throw(
        ArgumentError(
            "no vendored sphere-point set is available for order $order; the populated order pool runs from $(first(cache.orders)) to $(last(cache.orders)) with $(length(cache.orders)) entries",
        ),
    )
    return set
end

"""
    curated_sphere_point_set_orders()

Return the explicit curated fixture subset currently vendored in
`data/angular/curated_sphere_points.toml`.
"""
function curated_sphere_point_set_orders()
    return copy(_curated_sphere_point_set_cache().orders)
end

"""
    curated_sphere_point_set(order::Int)

Load one read-only curated sphere-point-set fixture by order/cardinality.

This helper is intentionally narrow: it exposes only the pinned curated subset
used for tiny tests, examples, and paper-stable fixtures.
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
