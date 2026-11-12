"""
    CubatureSpherePointSet

Read-only positive-weight spherical cubature entry for coiflet angular shells.
"""
struct CubatureSpherePointSet
    name::String
    family::Symbol
    point_set::SpherePointSet
    cubature_degree::Int
    default_l_inject::Int
    default_l_flat::Int
    default_beta::Float64
    cubature_weights::Vector{Float64}
end

const _CUBATURE_SPHERE_CACHE =
    Ref{Union{Nothing,Tuple{Vector{String},Dict{String,CubatureSpherePointSet}}}}(nothing)

function _cubature_sphere_points_toml_path()
    return _package_data_path("angular", "cubature_sphere_points.toml")
end

function _load_cubature_sphere_point_sets(path::AbstractString = _cubature_sphere_points_toml_path())
    raw = TOML.parsefile(path)
    meta, sets_raw = raw["meta"], raw["sets"]
    names = String[String(x) for x in meta["names"]]
    entries = Dict{String,CubatureSpherePointSet}()
    for name in names
        item = sets_raw[name]
        order = Int(item["order"])
        rows = item["coordinates"]
        coords = [Float64(rows[i][d]) for i in 1:order, d in 1:3]
        weights = Float64[Float64(w) for w in item["cubature_weights"]]
        minimum(weights) > 0.0 || throw(ArgumentError("cubature weights must be strictly positive for $name"))
        prov = SpherePointSetProvenance(
            String(meta["source_tag"]), String(meta["source_project"]),
            String(meta["source_artifact"]), String(item["note"]),
        )
        pset = SpherePointSet(order, order, coords, Float64(item["nn_ratio"]), prov)
        entries[name] = CubatureSpherePointSet(
            name, Symbol(String(item["family"])), pset, Int(item["cubature_degree"]),
            Int(item["default_l_inject"]), Int(item["default_l_flat"]), Float64(item["default_beta"]), weights,
        )
    end
    return (names, entries)
end

function _cubature_sphere_point_set_cache()
    _CUBATURE_SPHERE_CACHE[] === nothing && (_CUBATURE_SPHERE_CACHE[] = _load_cubature_sphere_point_sets())
    return _CUBATURE_SPHERE_CACHE[]
end

"""
    cubature_sphere_point_set_names()

Return a copied list of available positive-cubature sphere point-set names.
"""
cubature_sphere_point_set_names() = copy(_cubature_sphere_point_set_cache()[1])

"""
    cubature_sphere_point_set(name)

Return one cached positive-cubature entry. Treat its arrays as read-only.
"""
function cubature_sphere_point_set(name::Union{Symbol,AbstractString})
    names, entries = _cubature_sphere_point_set_cache()
    entry = get(entries, String(name), nothing)
    entry === nothing && throw(ArgumentError("unknown cubature sphere point set $(repr(String(name))); available: $(join(names, ", "))"))
    return entry
end
