"""
    SpherePointSetProvenance

Minimal provenance record for one sphere-point-set source or derivation path.

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

Read-only sphere-point-set entry for the experimental angular research track.

The normal/default angular backing store is the full vendored JLD2 point-set
collection under `data/angular/SpherePoints.jld2`. The smaller curated subset
remains available through explicit fixture helpers, while deterministic
Fibonacci and explicit optimization paths remain opt-in.
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

function _sphere_point_set_nn_ratio(coordinates::AbstractMatrix{<:Real})
    npoints = size(coordinates, 1)
    npoints >= 2 || return 1.0
    min_nn = Inf
    max_nn = 0.0
    for i in 1:npoints
        best = Inf
        xi = view(coordinates, i, :)
        for j in 1:npoints
            i == j && continue
            best = min(best, norm(xi .- view(coordinates, j, :)))
        end
        min_nn = min(min_nn, best)
        max_nn = max(max_nn, best)
    end
    return max_nn / min_nn
end

function _fibonacci_sphere_coordinates(order::Int)
    order >= 1 || throw(ArgumentError("fibonacci_sphere_point_set requires order >= 1"))
    phi_golden = (1.0 + sqrt(5.0)) / 2.0
    coordinates = Matrix{Float64}(undef, order, 3)
    for i in 0:(order - 1)
        z = 1.0 - 2.0 * (i + 0.5) / order
        r = sqrt(max(0.0, 1.0 - z * z))
        phi = 2.0 * pi * i / phi_golden
        coordinates[i + 1, 1] = r * cos(phi)
        coordinates[i + 1, 2] = r * sin(phi)
        coordinates[i + 1, 3] = z
    end
    return coordinates
end

function _pack_sphere_point_angles(coordinates::AbstractMatrix{<:Real})
    npoints = size(coordinates, 1)
    angles = Vector{Float64}(undef, 2 * npoints)
    for i in 1:npoints
        x, y, z = coordinates[i, 1], coordinates[i, 2], coordinates[i, 3]
        z_clamped = clamp(z, -1 + 1.0e-15, 1 - 1.0e-15)
        angles[2 * i - 1] = atanh(z_clamped)
        angles[2 * i] = atan(y, x)
    end
    return angles
end

function _wrap_azimuth(phi::Real)
    return mod(Float64(phi) + pi, 2pi) - pi
end

function _unpack_sphere_point_angles(angles::AbstractVector{<:Real})
    iseven(length(angles)) ||
        throw(ArgumentError("packed sphere-point angles must have even length"))
    npoints = length(angles) ÷ 2
    coordinates = Matrix{Float64}(undef, npoints, 3)
    for i in 1:npoints
        xi = Float64(angles[2 * i - 1])
        phi = _wrap_azimuth(angles[2 * i])
        z = tanh(xi)
        r = sqrt(max(0.0, 1.0 - z * z))
        coordinates[i, 1] = r * cos(phi)
        coordinates[i, 2] = r * sin(phi)
        coordinates[i, 3] = z
    end
    return coordinates
end

@inline function _sphere_point_sinhc_and_deriv(x::Float64)
    ax = abs(x)
    if ax < 1.0e-8
        x2 = x * x
        s = 1.0 + x2 / 6.0 + x2 * x2 / 120.0
        ds = x / 3.0 + x2 * x / 30.0
        return (s, ds)
    end
    return (sinh(x) / x, (x * cosh(x) - sinh(x)) / (x * x))
end

@inline function _sphere_point_sinhc_scalar(x::Float64)
    s, _ = _sphere_point_sinhc_and_deriv(x)
    return s
end

@inline function _sphere_point_kappa_from_beta(order::Int, beta::Real)
    theta_nn = sqrt(4.0 * pi / order)
    return (Float64(beta) / theta_nn)^2
end

function _sphere_point_gaussian_gram(
    coordinates::AbstractMatrix{<:Real},
    kappa::Real;
    return_r::Bool = false,
)
    npoints = size(coordinates, 1)
    dot_products = Matrix{Float64}(coordinates * transpose(coordinates))
    for j in 1:npoints
        for i in 1:npoints
            dot_products[i, j] = clamp(dot_products[i, j], -1.0 + 1.0e-12, 1.0)
        end
    end

    radii = Matrix{Float64}(undef, npoints, npoints)
    kappaf = Float64(kappa)
    for j in 1:npoints
        for i in 1:npoints
            t = max(0.0, 2.0 * (1.0 + dot_products[i, j]))
            radii[i, j] = kappaf * sqrt(t)
        end
    end

    gram = Matrix{Float64}(undef, npoints, npoints)
    prefactor = 4.0 * pi * exp(-2.0 * kappaf)
    for j in 1:npoints
        for i in 1:npoints
            gram[i, j] = prefactor * _sphere_point_sinhc_scalar(radii[i, j])
        end
    end
    for i in 1:npoints
        gram[i, i] += 1.0e-15
    end

    if return_r
        return gram, radii
    end
    return gram
end

function _sphere_point_logdet(gram::AbstractMatrix{<:Real})
    eig = eigvals(Symmetric(Matrix{Float64}(gram)))
    shift = minimum(eig) < 0.0 ? max(1.0e-15, -2.0 * minimum(eig)) : 0.0
    stabilized = shift == 0.0 ? Matrix{Float64}(gram) : Matrix{Float64}(gram) + Diagonal(fill(shift, size(gram, 1)))
    factor = cholesky(Symmetric(stabilized); check = true)
    return 2.0 * sum(log, diag(factor.U))
end

function _sphere_point_objective_and_gradient(
    angles::AbstractVector{<:Real},
    kappa::Real;
    compute_gradient::Bool = true,
)
    npoints = length(angles) ÷ 2
    z = Vector{Float64}(undef, npoints)
    r = Vector{Float64}(undef, npoints)
    cphi = Vector{Float64}(undef, npoints)
    sphi = Vector{Float64}(undef, npoints)
    coordinates = Matrix{Float64}(undef, npoints, 3)

    for i in 1:npoints
        xi = Float64(angles[2 * i - 1])
        phi = _wrap_azimuth(angles[2 * i])
        zi = tanh(xi)
        ri = sqrt(max(0.0, 1.0 - zi * zi))
        cph = cos(phi)
        sph = sin(phi)
        z[i] = zi
        r[i] = ri
        cphi[i] = cph
        sphi[i] = sph
        coordinates[i, 1] = ri * cph
        coordinates[i, 2] = ri * sph
        coordinates[i, 3] = zi
    end

    gram, radii = _sphere_point_gaussian_gram(coordinates, kappa; return_r = true)
    objective = -_sphere_point_logdet(gram)
    compute_gradient || return (objective = objective, gradient = nothing, coordinates = coordinates)

    eig = eigvals(Symmetric(gram))
    shift = minimum(eig) < 0.0 ? max(1.0e-15, -2.0 * minimum(eig)) : 0.0
    stabilized = shift == 0.0 ? gram : gram + Diagonal(fill(shift, size(gram, 1)))
    factor = cholesky(Symmetric(stabilized); check = true)
    inv_u = inv(factor.U)
    gram_inv = inv_u * transpose(inv_u)

    gradient = Vector{Float64}(undef, length(angles))
    prefactor = 4.0 * pi * exp(-2.0 * Float64(kappa))
    kappa2 = Float64(kappa)^2

    for i in 1:npoints
        gi1 = 0.0
        gi2 = 0.0
        gi3 = 0.0
        for j in 1:npoints
            rij = radii[i, j]
            _, ds = _sphere_point_sinhc_and_deriv(rij)
            aij = prefactor * ds * (kappa2 / rij)
            coeff = 2.0 * gram_inv[i, j] * aij
            gi1 += coeff * coordinates[j, 1]
            gi2 += coeff * coordinates[j, 2]
            gi3 += coeff * coordinates[j, 3]
        end

        gi1 = -gi1
        gi2 = -gi2
        gi3 = -gi3

        zi = z[i]
        ri = r[i]
        cph = cphi[i]
        sph = sphi[i]
        ds_dxi = -zi * ri
        dz_dxi = ri * ri

        gradient[2 * i - 1] = gi1 * (ds_dxi * cph) + gi2 * (ds_dxi * sph) + gi3 * dz_dxi
        gradient[2 * i] = gi1 * (-ri * sph) + gi2 * (ri * cph)
    end

    return (objective = objective, gradient = gradient, coordinates = coordinates)
end

function _optimized_sphere_point_coordinates(
    coordinates::AbstractMatrix{<:Real};
    beta::Real = 2.0,
    iters::Int = 200,
    gtol::Real = 1.0e-8,
    trace::Bool = false,
)
    npoints = size(coordinates, 1)
    iters >= 0 || throw(ArgumentError("optimize_sphere_point_set requires iters >= 0"))
    gtol > 0 || throw(ArgumentError("optimize_sphere_point_set requires gtol > 0"))
    beta > 0 || throw(ArgumentError("optimize_sphere_point_set requires beta > 0"))

    angles = _pack_sphere_point_angles(coordinates)
    kappa = _sphere_point_kappa_from_beta(npoints, beta)
    current = _sphere_point_objective_and_gradient(angles, kappa)
    trace &&
        @info(
            "sphere-point optimization start",
            order = npoints,
            beta = Float64(beta),
            objective = current.objective,
            gradient_supnorm = norm(current.gradient, Inf),
        )

    for iteration in 1:iters
        gradient_supnorm = norm(current.gradient, Inf)
        gradient_supnorm <= gtol && break

        direction = -current.gradient
        slope = dot(current.gradient, direction)
        step = 1.0
        accepted = false
        best = current
        best_angles = angles

        while step >= 1.0e-12
            trial_angles = copy(angles)
            trial_angles .+= step .* direction
            for i in 1:npoints
                trial_angles[2 * i] = _wrap_azimuth(trial_angles[2 * i])
            end

            trial = _sphere_point_objective_and_gradient(trial_angles, kappa)
            if trial.objective <= current.objective + 1.0e-4 * step * slope
                accepted = true
                best = trial
                best_angles = trial_angles
                break
            end
            step *= 0.5
        end

        accepted || break
        angles = best_angles
        current = best
        trace &&
            (iteration == 1 || iteration % 10 == 0 || iteration == iters) &&
            @info(
                "sphere-point optimization step",
                iteration = iteration,
                objective = current.objective,
                gradient_supnorm = norm(current.gradient, Inf),
            )
    end

    return current.coordinates
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
    fibonacci_sphere_point_set(order::Int)

Build one deterministic Fibonacci sphere point set using the canonical
`z = 1 - 2*(i + 0.5)/N`, `phi = 2π*i/φ_golden` convention.

This is an explicit debug/experimental path. It does not affect the vendored
default `sphere_point_set(order)` behavior.
"""
function fibonacci_sphere_point_set(order::Int)
    coordinates = _fibonacci_sphere_coordinates(order)
    provenance = SpherePointSetProvenance(
        "deterministic_fibonacci_seed",
        "GaussletBases",
        "in_memory_generation",
        "canonical Fibonacci sphere seed with i=0:N-1, z=1-2*(i+0.5)/N, phi=2π*i/φ_golden; no randomization or rotation",
    )
    return SpherePointSet(
        order,
        order,
        coordinates,
        _sphere_point_set_nn_ratio(coordinates),
        provenance,
    )
end

"""
    optimize_sphere_point_set(
        set::SpherePointSet;
        beta=2.0,
        iters=200,
        gtol=1e-8,
        trace=false,
    )

Optimize an existing sphere point set explicitly on request.

This is an experimental deterministic optimization path layered on top of an
existing point set. It does not run implicitly inside `sphere_point_set(order)`.
"""
function optimize_sphere_point_set(
    set::SpherePointSet;
    beta::Real = 2.0,
    iters::Int = 200,
    gtol::Real = 1.0e-8,
    trace::Bool = false,
)
    coordinates = _optimized_sphere_point_coordinates(
        set.coordinates;
        beta = beta,
        iters = iters,
        gtol = gtol,
        trace = trace,
    )
    provenance = SpherePointSetProvenance(
        "optimized_from_input_point_set",
        "GaussletBases",
        "in_memory_optimization",
        "explicit Gaussian-Gram logdet optimization from $(repr(set.provenance.source_tag)) with beta=$(Float64(beta)), iters=$(iters), gtol=$(Float64(gtol))",
    )
    return SpherePointSet(
        set.order,
        set.cardinality,
        coordinates,
        _sphere_point_set_nn_ratio(coordinates),
        provenance,
    )
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
