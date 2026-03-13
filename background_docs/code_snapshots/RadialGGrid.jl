module RadialGGrid

# ==============================================================================
# RadialGGrid.jl — radial boundary-gausslet machinery (r ≥ 0)
# ==============================================================================
#
# What this module provides
# - A localized, (approximately) orthonormal radial basis on the half-line r ≥ 0
#   built from full-line Gausslets.jl gausslets by enforcing the r=0 boundary
#   condition and COMX-localizing the span.
# - A quadrature grid on r ≥ 0 built by mapping a uniform auxiliary coordinate s
#   through a smooth erf "switch-on" (to avoid a hard endpoint at r=0) and then
#   through a user radial CoordinateMapping.
# - One-electron radial matrices on that basis:
#     S, T, V = -Z/r, Vcentr = 1/r^2 (centrifugal prefactor).
# - Coulomb radial multipole tables `Veel[L]` for use by Ylm/Gaunt machinery.
#   Two builders exist:
#     :trap  → legacy trapezoid/prefix-integral route (fast, lower order).
#     :cheb  → higher-accuracy Clenshaw–Curtis on piecewise-u segments (default).
#
# Boundary-gausslet cache
# - A JLD2 cache (default `BoundaryGausslets.jld2`) stores the half-line basis as
#   coefficients `C` in a finite window of raw full-line gausslets indexed by
#   integers `js`.  Each column of `C` is one boundary basis function.
# - For the default `boundary_method=:odd_even` the cache also stores a piecewise
#   Chebyshev representation of χ_a(u) and dχ/du for fast/accurate evaluation.
#
# Parameter guidance (high level)
# - Accuracy knobs typically come in two layers:
#   (1) physics/basis knobs: `corespacing`, `s`, `Rmax`, boundary cache size.
#   (2) integration knobs: `hgrid` (runtime quadrature), and the cache-build grid.
# - `build_radial_and_one_electron` can refine `hgrid` automatically by halving
#   it until `||S - I||_Inf ≤ hgrid_tol` (default on) so silent integration
#   regressions are less likely.

export RadialData, radial_moments, xdiag_finalize, make_erf_mapped_grid, build_radial_and_one_electron,
       build_Veel_L, build_Veel_L_cheb, build_Veel_tables, radial_states

using LinearAlgebra
using CoordinateMapping
using Gausslets
using JLD2
using SpecialFunctions
using Logging
using Util

using PiecewiseChebyshev

"""
    RadialData

Container returned by `build_radial_and_one_electron`.

The intent is that higher-level drivers (HF, Ylm Coulomb, post-HF, etc.) can
depend only on this payload without needing to know how the radial basis was
constructed.
"""
struct RadialData
    r::Vector{Float64}        # physical radii for the erf-mapped grid (not the gausslet centers)
    w::Vector{Float64}        # quadrature weights associated with r
    chiv::Matrix{Float64}     # boundary-gausslet values chi_a(r_i) on the grid
    chivp::Matrix{Float64}    # radial derivatives chi'_a(r_i) evaluated on the grid
    invr::Vector{Float64}     # pointwise 1/r_i values
    S::Matrix{Float64}        # overlap matrix <chi_a|chi_b>
    T::Matrix{Float64}        # kinetic matrix <chi_a| -0.5 d^2/dr^2 |chi_b>
    V::Matrix{Float64}        # nuclear attraction <chi_a| -Z/r |chi_b>
    Vcentr::Matrix{Float64}   # centrifugal <chi_a| ell(ell+1)/(2 r^2) |chi_b> (ell-independent factor)
    wts::Vector{Float64}      # integrated weights of the gausslets (needed for IDA scaling)
    N::Int                    # number of boundary-gausslet functions kept
    Ngr::Int                  # number of radial grid points
    centers::Vector{Float64}  # physical centers of each boundary-gausslet
    evalcache::Any            # optional runtime eval payload for higher-order quadratures (e.g. Vee)
end

"""
    radial_moments(data::RadialData) -> (R1, R2)

Compute radial moment matrices <r> and <r^2> from the quadrature grid.
"""
function radial_moments(data::RadialData)
    r = data.r
    w = data.w
    chiv = data.chiv
    rw = r .* w
    r2w = r .* rw
    R1 = transpose(chiv) * (rw .* chiv)
    R2 = transpose(chiv) * (r2w .* chiv)
    R1 .= 0.5 .* (R1 .+ transpose(R1))
    R2 .= 0.5 .* (R2 .+ transpose(R2))
    return R1, R2
end

"""
    xdiag_finalize(data::RadialData; veel=nothing, rotate_evalcache=false, verbose=false) -> RadialData

Diagonalize X = <r> in the mapped basis (no generalized metric), rotate the one-electron
matrices and basis values accordingly, and optionally rotate IDA Vee tables using
the correct wts-scaled similarity transform.

Notes:
- If `rotate_evalcache=false` (default), the returned `RadialData` has
  `evalcache = nothing` to prevent accidental reuse.
- If `rotate_evalcache=true`, the Chebyshev/gaussian evalcache is rotated so
  `build_Veel_tables` can be called on the rotated basis.
"""
function xdiag_finalize(data::RadialData;
        veel::Union{Nothing,Vector{Matrix{Float64}}}=nothing,
        rotate_evalcache::Bool=false,
        verbose::Bool=false)
    R1, _ = radial_moments(data)
    X = Matrix(Symmetric(R1))
    evals, U = eigen(X)
    wts_raw = U' * data.wts
    signs = sign.(wts_raw)
    @inbounds for i in eachindex(signs)
        signs[i] == 0.0 && (signs[i] = 1.0)
    end
    U = U * Diagonal(signs)

    if verbose
        off = norm(X - Diagonal(diag(X)))
        @info "RadialGGrid: final X diag (offdiag norm=$(off))"
    end

    wts_new = abs.(wts_raw)
    chiv  = data.chiv * U
    chivp = data.chivp * U
    S     = U' * data.S * U
    T     = U' * data.T * U
    V     = U' * data.V * U
    Vcentr = U' * data.Vcentr * U
    wts   = wts_new
    centers = evals

    for M in (S, T, V, Vcentr)
        M .= 0.5 .* (M .+ transpose(M))
    end

    if veel !== nothing
        wts_old = data.wts
        any(wts_old .<= 0.0) && error("RadialGGrid: xdiag Vee rotation requires positive wts; min(wts)=$(minimum(wts_old))")
        Dold = Diagonal(wts_old)
        Dnew_inv = Diagonal(1.0 ./ wts_new)
        for L in eachindex(veel)
            B = Dold * veel[L] * Dold
            Brot = U' * B * U
            VL = Dnew_inv * Brot * Dnew_inv
            VL .= 0.5 .* (VL .+ transpose(VL))
            veel[L] = VL
        end
    end

    evalcache = nothing
    if rotate_evalcache && data.evalcache !== nothing
        ec = data.evalcache
        cheb = ec.cheb
        gaussian = ec.gaussian

        cheb_rot = nothing
        if cheb !== nothing
            coeffs = cheb.coeffs
            s1, s2, n = size(coeffs)
            if n == size(U, 1)
                coeffs2 = reshape(coeffs, s1 * s2, n) * U
                coeffs_new = reshape(coeffs2, s1, s2, n)
                cheb_rot = (; cheb..., coeffs=coeffs_new)
            else
                verbose && @warn "RadialGGrid: X-diag dropped cheb evalcache (ncheb=$(n) != N=$(size(U,1)))"
                cheb_rot = nothing
            end
        end

        gaussian_rot = nothing
        if gaussian !== nothing
            G = gaussian.G
            size(G, 2) == size(U, 1) || error("RadialGGrid: gaussian evalcache width=$(size(G,2)) != N=$(size(U,1))")
            gaussian_rot = (; gaussian..., G=G * U)
        end

        inj_rot = nothing
        if hasproperty(ec, :inj) && ec.inj !== nothing
            inj = ec.inj
            if size(inj.Cinj, 2) == size(U, 1)
                inj_rot = (; inj..., Cinj=inj.Cinj * U)
            else
                verbose && @warn "RadialGGrid: X-diag dropped inj evalcache (ninj=$(size(inj.Cinj,2)) != N=$(size(U,1)))"
            end
        end

        evalcache = (; ec..., cheb=cheb_rot, gaussian=gaussian_rot, inj=inj_rot, N=ec.N)
    end

    return RadialData(
        data.r,
        data.w,
        chiv,
        chivp,
        data.invr,
        S,
        T,
        V,
        Vcentr,
        wts,
        data.N,
        data.Ngr,
        centers,
        evalcache
    )
end

# ------------------------------------------------------------------------------
# Defaults and tunables
# ------------------------------------------------------------------------------
#
# Boundary cache defaults (used when building/reading BoundaryGausslets.jld2)
# - `_BOUNDARY_*`: control how the *cached* half-line basis is constructed.
#   These are not performance-critical at runtime (unless you rebuild the cache).
# - `jpos_save` (aka `L` for odd/even): increases cache size and far-tail coverage.
# - `neven_add`: adds even-Taylor content near the origin; too large can increase
#   cancellation without materially improving physics targets.
# - `grid_h/grid_sigma/grid_s0/rgridmax`: internal quadrature for cache building.
#
# Boundary Chebyshev defaults
# - `_BOUNDARY_CHEB_*`: control the piecewise Chebyshev fit of χ(u) and dχ/du
#   stored in the cache (used at runtime for fast/accurate evaluation).
#
# Vee defaults
# - `_VEE_*`: control the high-accuracy Coulomb table builder (default `:cheb`).
#   Increase `deg` or decrease segment widths if you see Vee-integration drift.
const DEFAULT_BOUNDARY_CACHE = "BoundaryGausslets.jld2"

const _BOUNDARY_H            = 1.0                  # seed gausslet spacing in u-space (Gausslets.jl default is 1)
const _BOUNDARY_FAMILY       = :cf1092              # Gausslets.jl stencil family for seed gausslets
const _BOUNDARY_JPOS_BUILD   = 40                   # legacy build window max +j (for :add_tails only)
const _BOUNDARY_JPOS_SAVE    = 24                   # cached window size: store functions with centers <≈ jpos_save (u-space)
const _BOUNDARY_JNEG_LIST    = 2:8                  # legacy left-tail sizes to prebuild (for :add_tails groups jneg_k)
const _BOUNDARY_DEFAULT_METHOD = :odd_even          # default boundary cache method (:odd_even or :add_tails)
const _BOUNDARY_ODD_EVEN_NEVEN_ADD = 8              # default number of even-complement modes added in :odd_even
const _BOUNDARY_ODD_EVEN_TAIL_K = 6                 # default max even index K (keep E_k for k=0:K) in :odd_even_tail
const _BOUNDARY_ODD_EVEN_TAIL_VERSION = 2           # bump when :odd_even_tail construction changes (forces cache rebuild)
const _BOUNDARY_GRID_H       = 0.01                 # cache-build quadrature: base step size in u(s)=h*s*g(s)
const _BOUNDARY_GRID_SIGMA   = 3.0                  # cache-build quadrature: erf switch width (sigma in g(s))
const _BOUNDARY_GRID_S0      = 6.5                  # cache-build quadrature: erf switch offset (s0 in g(s))
const _BOUNDARY_RGRIDMAX     = 80.0                 # cache-build quadrature: integrate out to this physical radius
const _BOUNDARY_S_EIG_TOL    = 1e-12                # overlap eigenvalue cutoff when orthonormalizing cache spans

# Boundary injection defaults (extra gpoly primitives added during cache build)
const _BOUNDARY_INJECT_VERSION = 2
const _BOUNDARY_INJECT_COUNT = 2                    # number of injected gaussians (0,1,2)
const _BOUNDARY_INJECT_MODE = :gpoly                # currently only :gpoly is supported
const _BOUNDARY_INJECT_SEARCH = true                # search alphas during cache build
const _BOUNDARY_INJECT_ALPHA1 = 0.07
const _BOUNDARY_INJECT_ALPHA2 = 0.02
const _BOUNDARY_INJECT_ALPHA1_RANGE = (0.03, 0.12)
const _BOUNDARY_INJECT_ALPHA2_RANGE = (0.005, 0.05)
const _BOUNDARY_INJECT_ALPHA1_GRID_N = 8
const _BOUNDARY_INJECT_ALPHA2_GRID_N = 8
const _BOUNDARY_INJECT_ALPHA_REFINE_STEP = 0.01
const _BOUNDARY_INJECT_ALPHA_REFINE_TOL = 1e-4
const _BOUNDARY_INJECT_ALPHA_REFINE_MAXITER = 40
const _BOUNDARY_INJECT_MULTI_START = false
const _BOUNDARY_INJECT_MULTI_START_N = 6
const _BOUNDARY_INJECT_NM_MAXITER = 200
const _BOUNDARY_INJECT_VERBOSE = true

const _BOUNDARY_CHEB_DEG     = 64                   # piecewise Chebyshev degree for cached χ(u)
const _BOUNDARY_CHEB_UPAD    = 20.0                 # extra u-range beyond max(js) to cover tail evaluation
const _BOUNDARY_CHEB_TOL_REL = 1e-14                # relative tolerance for χ(u) fit vs high-precision truth
const _BOUNDARY_CHEB_TOL_REL_DPHI = 1e-12           # relative tolerance for dχ/du fit vs truth
const _BOUNDARY_CHEB_BIG_PREC = 256                 # BigFloat precision (bits) used while fitting χ(u)
const _BOUNDARY_CHEB_MAXSEG  = 1024                 # max number of Chebyshev segments allowed in cache
const _BOUNDARY_CHEB_MINWIDTH = 0.05                # minimum segment width in u (prevents over-splitting)
const _BOUNDARY_CHEB_CHECK_DEG = 128                # higher degree used for fit verification on each segment
const _BOUNDARY_CHEB_VERSION = 4                    # cache format version for cheb payload (bump if on-disk schema changes)

const _VEE_DEFAULT_METHOD        = :cheb             # default build_Veel_tables method (:cheb or :trap)
const _VEE_CHEB_DEG              = 64                # CC degree per u-segment for Vee tables
const _VEE_CHEB_DEG_MAX          = 256               # max CC degree if degree-refinement is enabled
const _VEE_CHEB_TAIL_DU          = 0.2               # tail segment width in u beyond cache-aligned region
const _VEE_CHEB_GEOM_RATIO       = 1.5               # geometric spacing ratio near u≈0 (stability for large L)
const _VEE_CHEB_GEOM_UMAX        = 1.0               # end of near-0 geometric region in u
const _VEE_CHEB_REFINE_MINWIDTH  = 0.10              # split any segment wider than this during edge refinement
const _VEE_CHEB_REFINE_MAXITER   = 6                 # edge-refine rounds (split long segments each round)
const _VEE_CHEB_TOL_REL          = 1e-12             # relative tolerance for degree-refinement consistency check
const _VEE_CHEB_DEG_MAXITER      = 0                 # number of degree-doubling rounds; 0 disables self-check
const _VEE_CHEB_MAX_EDGES        = 8192              # cap on number of u-edges (segments = edges-1)

# boundary_cache_group: map cache settings → JLD2 group name.
function boundary_cache_group(method::Symbol;
        jneg::Int,
        jpos_save::Int,
        neven_add::Int,
        even_tail_K::Int)
    if method == :add_tails
        return "jneg_$(jneg)"
    elseif method == :odd_even
        return "odd_even_L$(jpos_save)_e$(neven_add)"
    elseif method == :odd_even_tail
        return "odd_even_tail_L$(jpos_save)_K$(even_tail_K)"
    else
        error("Unknown boundary cache method: $(method). Expected :odd_even, :odd_even_tail, or :add_tails.")
    end
end

# cache metadata helpers for injected gpoly primitives
_float_match(a::Real, b::Real) = isapprox(a, b; rtol=1e-12, atol=0.0)

function _inject_meta_ok(f, grp::AbstractString;
        inject_count::Int,
        inject_mode::Symbol,
        inject_search::Bool,
        inject_alpha1::Real,
        inject_alpha2::Real,
        inject_alpha1_range::Tuple{Real,Real},
        inject_alpha2_range::Tuple{Real,Real},
        inject_alpha1_grid_n::Int,
        inject_alpha2_grid_n::Int,
        inject_alpha_refine_step::Real,
        inject_alpha_refine_tol::Real,
        inject_alpha_refine_maxiter::Int,
        inject_multi_start::Bool,
        inject_multi_start_n::Int,
        inject_nm_maxiter::Int)

    key = "$grp/inject/version"
    haskey(f, key) || return false
    Int(f[key]) == _BOUNDARY_INJECT_VERSION || return false
    Int(f["$grp/inject/count"]) == inject_count || return false
    String(f["$grp/inject/mode"]) == String(inject_mode) || return false
    Bool(f["$grp/inject/search"]) == inject_search || return false
    inject_count == 0 && return true

    if inject_search
        ar1 = Vector{Float64}(f["$grp/inject/alpha1_range"])
        ar2 = Vector{Float64}(f["$grp/inject/alpha2_range"])
        length(ar1) == 2 || return false
        length(ar2) == 2 || return false
        _float_match(ar1[1], inject_alpha1_range[1]) || return false
        _float_match(ar1[2], inject_alpha1_range[2]) || return false
        _float_match(ar2[1], inject_alpha2_range[1]) || return false
        _float_match(ar2[2], inject_alpha2_range[2]) || return false
        Int(f["$grp/inject/alpha1_grid_n"]) == inject_alpha1_grid_n || return false
        Int(f["$grp/inject/alpha2_grid_n"]) == inject_alpha2_grid_n || return false
        _float_match(f["$grp/inject/alpha_refine_step"], inject_alpha_refine_step) || return false
        _float_match(f["$grp/inject/alpha_refine_tol"], inject_alpha_refine_tol) || return false
        Int(f["$grp/inject/alpha_refine_maxiter"]) == inject_alpha_refine_maxiter || return false
        Bool(f["$grp/inject/multi_start"]) == inject_multi_start || return false
        Int(f["$grp/inject/multi_start_n"]) == inject_multi_start_n || return false
        Int(f["$grp/inject/nm_maxiter"]) == inject_nm_maxiter || return false
    else
        _float_match(f["$grp/inject/alpha1"], inject_alpha1) || return false
        if inject_count >= 2
            _float_match(f["$grp/inject/alpha2"], inject_alpha2) || return false
        end
    end
    return true
end

function _write_inject_meta!(f, grp::AbstractString;
        inject_count::Int,
        inject_mode::Symbol,
        inject_search::Bool,
        alpha1::Real,
        alpha2::Real,
        inject_alpha1_range::Tuple{Real,Real},
        inject_alpha2_range::Tuple{Real,Real},
        inject_alpha1_grid_n::Int,
        inject_alpha2_grid_n::Int,
        inject_alpha_refine_step::Real,
        inject_alpha_refine_tol::Real,
        inject_alpha_refine_maxiter::Int,
        inject_multi_start::Bool,
        inject_multi_start_n::Int,
        inject_nm_maxiter::Int)

    f["$grp/inject/version"] = _BOUNDARY_INJECT_VERSION
    f["$grp/inject/count"] = inject_count
    f["$grp/inject/mode"] = String(inject_mode)
    f["$grp/inject/search"] = inject_search
    f["$grp/inject/alpha1"] = Float64(alpha1)
    f["$grp/inject/alpha2"] = Float64(alpha2)
    f["$grp/inject/alpha1_range"] = [Float64(inject_alpha1_range[1]), Float64(inject_alpha1_range[2])]
    f["$grp/inject/alpha2_range"] = [Float64(inject_alpha2_range[1]), Float64(inject_alpha2_range[2])]
    f["$grp/inject/alpha1_grid_n"] = inject_alpha1_grid_n
    f["$grp/inject/alpha2_grid_n"] = inject_alpha2_grid_n
    f["$grp/inject/alpha_refine_step"] = Float64(inject_alpha_refine_step)
    f["$grp/inject/alpha_refine_tol"] = Float64(inject_alpha_refine_tol)
    f["$grp/inject/alpha_refine_maxiter"] = inject_alpha_refine_maxiter
    f["$grp/inject/multi_start"] = inject_multi_start
    f["$grp/inject/multi_start_n"] = inject_multi_start_n
    f["$grp/inject/nm_maxiter"] = inject_nm_maxiter
end

# ensure_boundary_cache!: verify the cache exists (and has the needed group), else rebuild it.
function ensure_boundary_cache!(filename::AbstractString;
        method::Symbol=_BOUNDARY_DEFAULT_METHOD,
        jneg::Int=last(_BOUNDARY_JNEG_LIST),
        jpos_save::Int=_BOUNDARY_JPOS_SAVE,
        neven_add::Int=_BOUNDARY_ODD_EVEN_NEVEN_ADD,
        even_tail_K::Int=_BOUNDARY_ODD_EVEN_TAIL_K,
        inject_count::Int=_BOUNDARY_INJECT_COUNT,
        inject_mode::Symbol=_BOUNDARY_INJECT_MODE,
        inject_search::Bool=_BOUNDARY_INJECT_SEARCH,
        inject_verbose::Bool=_BOUNDARY_INJECT_VERBOSE,
        inject_alpha1::Real=_BOUNDARY_INJECT_ALPHA1,
        inject_alpha2::Real=_BOUNDARY_INJECT_ALPHA2,
        inject_alpha1_range::Tuple{Real,Real}=_BOUNDARY_INJECT_ALPHA1_RANGE,
        inject_alpha2_range::Tuple{Real,Real}=_BOUNDARY_INJECT_ALPHA2_RANGE,
        inject_alpha1_grid_n::Int=_BOUNDARY_INJECT_ALPHA1_GRID_N,
        inject_alpha2_grid_n::Int=_BOUNDARY_INJECT_ALPHA2_GRID_N,
        inject_alpha_refine_step::Real=_BOUNDARY_INJECT_ALPHA_REFINE_STEP,
        inject_alpha_refine_tol::Real=_BOUNDARY_INJECT_ALPHA_REFINE_TOL,
        inject_alpha_refine_maxiter::Int=_BOUNDARY_INJECT_ALPHA_REFINE_MAXITER,
        inject_multi_start::Bool=_BOUNDARY_INJECT_MULTI_START,
        inject_multi_start_n::Int=_BOUNDARY_INJECT_MULTI_START_N,
        inject_nm_maxiter::Int=_BOUNDARY_INJECT_NM_MAXITER)
    build_jneg_list = (method == :add_tails) ? unique!(sort!(vcat(collect(_BOUNDARY_JNEG_LIST), [jneg]))) : _BOUNDARY_JNEG_LIST
    needs_cheb = (method == :odd_even) || (method == :odd_even_tail)
    if isfile(filename)
        grp = boundary_cache_group(method; jneg, jpos_save, neven_add, even_tail_K)
        (has_group, has_cheb, cheb_kind_ok, cheb_version_ok, tail_version_ok, inject_ok) = jldopen(filename, "r") do f
            has_grp = haskey(f, grp)
            has_edges = haskey(f, "$grp/cheb/edges")
            kind_ok = true
            ver_ok = true
            tail_ok = true
            inj_ok = false
            if needs_cheb && has_edges
                if haskey(f, "$grp/cheb/dcoeffs_kind")
                    kind_ok = String(f["$grp/cheb/dcoeffs_kind"]) == "dpdu"
                else
                    kind_ok = false
                end
                if haskey(f, "$grp/cheb/version")
                    ver_ok = Int(f["$grp/cheb/version"]) == _BOUNDARY_CHEB_VERSION
                else
                    ver_ok = false
                end
            end
            if method == :odd_even_tail
                if haskey(f, "meta/odd_even_tail/version")
                    tail_ok = Int(f["meta/odd_even_tail/version"]) == _BOUNDARY_ODD_EVEN_TAIL_VERSION
                else
                    tail_ok = false
                end
            end
            if has_grp
                if method == :add_tails
                    inj_ok = true
                else
                    inj_ok = _inject_meta_ok(f, grp;
                        inject_count, inject_mode, inject_search,
                        inject_alpha1, inject_alpha2,
                        inject_alpha1_range, inject_alpha2_range,
                        inject_alpha1_grid_n, inject_alpha2_grid_n,
                        inject_alpha_refine_step, inject_alpha_refine_tol, inject_alpha_refine_maxiter,
                        inject_multi_start, inject_multi_start_n, inject_nm_maxiter)
                end
            end
            has_grp, has_edges, kind_ok, ver_ok, tail_ok, inj_ok
        end
        if has_group && (!needs_cheb || (has_cheb && cheb_kind_ok && cheb_version_ok)) && tail_version_ok && inject_ok
            return
        end
        if !has_group
            @info "RadialGGrid: boundary cache missing group=$(grp) at $(filename); regenerating."
        elseif needs_cheb && !has_cheb
            @info "RadialGGrid: boundary cache group=$(grp) missing Chebyshev data at $(filename); regenerating."
        elseif needs_cheb && !cheb_kind_ok
            @info "RadialGGrid: boundary cache group=$(grp) has legacy Chebyshev format; regenerating."
        elseif needs_cheb && !cheb_version_ok
            @info "RadialGGrid: boundary cache group=$(grp) has stale Chebyshev version; regenerating."
        elseif method == :odd_even_tail && !tail_version_ok
            @info "RadialGGrid: boundary cache group=$(grp) has stale :odd_even_tail basis version; regenerating."
        elseif !inject_ok
            @info "RadialGGrid: boundary cache group=$(grp) has stale injection settings; regenerating."
        else
            @info "RadialGGrid: boundary cache group=$(grp) invalid; regenerating."
        end
        build_boundary_cache!(filename; jneg_list=build_jneg_list, jpos_save, neven_add, even_tail_K,
            inject_count, inject_mode, inject_search, inject_verbose,
            inject_alpha1, inject_alpha2,
            inject_alpha1_range, inject_alpha2_range,
            inject_alpha1_grid_n, inject_alpha2_grid_n,
            inject_alpha_refine_step, inject_alpha_refine_tol, inject_alpha_refine_maxiter,
            inject_multi_start, inject_multi_start_n, inject_nm_maxiter)
        return
    end
    @info "RadialGGrid: boundary cache missing at $(filename); generating."
    build_boundary_cache!(filename; jneg_list=build_jneg_list, jpos_save, neven_add, even_tail_K,
        inject_count, inject_mode, inject_search, inject_verbose,
        inject_alpha1, inject_alpha2,
        inject_alpha1_range, inject_alpha2_range,
        inject_alpha1_grid_n, inject_alpha2_grid_n,
        inject_alpha_refine_step, inject_alpha_refine_tol, inject_alpha_refine_maxiter,
        inject_multi_start, inject_multi_start_n, inject_nm_maxiter)
    isfile(filename) || error("Boundary cache generation failed; $(filename) not found.")
end

# ------------------------------------------------------------------------------
# Standalone erf-grid helpers (moved from Erfgrid.jl)
# ------------------------------------------------------------------------------

"""
    make_erf_grid(; h, rmax, sigma=1.0, s0=6.5) -> (r, w)

Build a smooth half-infinite quadrature grid on `[0, rmax]` suitable for
midpoint-style rules (no hard endpoint at `r=0`).

Implementation sketch
- Uses midpoint samples `s = k + 0.5` for `k=0,1,2,...`.
- Defines a smooth “switch-on” factor

    `g(s) = 0.5*erfc(s0 - s/sigma)`

  so `g(0.5)` is tiny when `s0≈6.5` and the effective first radius is extremely
  close to zero but not exactly at an endpoint.
- Uses `r(s) = h*s*g(s)` and returns:
  - `r[k] = r(s_k)`
  - `w[k] = dr/ds(s_k)` (the quadrature weight for unit steps in `s`)

Parameter meanings / tips
- `h`: step size in the auxiliary coordinate; smaller → more points (higher cost,
  higher accuracy).  For cache-building grids, values like `0.01` are typical.
- `rmax`: stop once `r ≥ rmax`.
- `sigma`: controls how quickly the switch turns on; larger → broader switch.
- `s0`: where the switch turns on. `s0≈6.5` is chosen so `erfc(s0)` is below
  double precision for O(1) numbers, effectively pushing the hard edge “to
  infinity” in `s`.
"""
function make_erf_grid(; h::Real, rmax::Real, sigma::Real=1.0, s0::Real=6.5)
    (h > 0)    || error("h must be > 0")
    (rmax > 0) || error("rmax must be > 0")
    r = Float64[]
    w = Float64[]
    k = 0
    rtpi = sqrt(pi)
    while true
        s  = k + 0.5
        z  = (s / sigma) - s0
        g  = 0.5 * SpecialFunctions.erfc(-z)
        gp = exp(-(z*z)) / (rtpi * sigma)
        rk  = h * s * g
        rpk = h * (g + s*gp)
        push!(r, rk)
        push!(w, rpk)
        rk >= rmax && break
        k += 1
    end
    return r, w
end

# integrate: convenience wrapper for ∫ v(r) dr given grid points `r` and weights `w`.
integrate(v, r, w) = sum(w .* v)

# integrate_r_p: helper for ∫ r^p v(r) dr on a nonuniform quadrature grid.
function integrate_r_p(v, r, w, p::Real; r_eps::Real=1e-14, clamp::Bool=true)
    if p < 0 && !clamp && any(==(0.0), r)
        error("integrate_r_p: r contains zero and p<0 with clamp=false")
    end
    if p == 0
        return sum(w .* v)
    elseif p == 1
        return sum(w .* r .* v)
    else
        rr = (p < 0 && clamp) ? max.(r, r_eps) : r
        return sum(w .* (rr .^ p) .* v)
    end
end

# ------------------------------------------------------------------------------
# Boundary Gausslet cache construction (inlined from make_boundary_cache.jl)
# ------------------------------------------------------------------------------

# phi_eval_factory: build a closure φ(j,x)=g(x-jH) for seed gausslet evaluation.
phi_eval_factory(g, H) = (j::Int, x::Float64) -> g(x - j * H)

# assemble_halfline_seed: build overlap/position matrices for a seed set restricted to r≥0.
function assemble_halfline_seed(phi, js::Vector{Int}, r, w; D::Union{Nothing,Vector{Float64}}=nothing)
    n = length(js)
    S = zeros(n, n)
    X = zeros(n, n)
    for a in 1:n
        ja = js[a]
        va = phi.(Ref(ja), r)
        if D !== nothing
            va .= va .* D[a]
        end
        for b in a:n
            jb = js[b]
            vb = phi.(Ref(jb), r)
            if D !== nothing
                vb .= vb .* D[b]
            end
            vv = va .* vb
            S_ab = integrate(vv, r, w)
            X_ab = integrate_r_p(vv, r, w, 1)
            S[a,b]=S_ab; S[b,a]=S_ab
            X[a,b]=X_ab; X[b,a]=X_ab
        end
    end
    return S, X
end

# seed_scalar_integrals: compute ||seed||^2 and ∫seed for each seed column, plus sampled values.
function seed_scalar_integrals(phi, js::Vector{Int}, r, w)
    n = length(js)
    sdiag = zeros(n)
    wseed = zeros(n)
    seedv = zeros(length(r), n)
    for a in 1:n
        ja = js[a]
        va = phi.(Ref(ja), r)
        seedv[:, a] = va
        sdiag[a] = integrate(va .* va, r, w)
        wseed[a] = integrate(va, r, w)
    end
    return sdiag, wseed, seedv
end

# pivot_projection_transform: remove the (approximate) delta-at-origin direction via a pivoted projection.
function pivot_projection_transform(g, js::Vector{Int}, D::Vector{Float64}, H)
    n = length(js)
    v = zeros(n)
    for (k,j) in enumerate(js)
        v[k] = D[k] * g(-j * H)
    end
    p = argmax(abs.(v))
    vp = v[p]
    vp == 0 && error("pivot_projection_transform: pivot v[p]==0")
    T = zeros(n, n - 1)
    col = 1
    for jidx in 1:n
        if jidx == p
            continue
        end
        T[jidx, col] = 1.0
        T[p, col]   -= v[jidx] / vp
        col += 1
    end
    return p, v, T
end

# s_invsqrt_reduced: compute Q*D^{-1/2} for an overlap matrix, dropping tiny modes.
function s_invsqrt_reduced(Sred::AbstractMatrix; tol=_BOUNDARY_S_EIG_TOL)
    λ, Q = eig(Matrix{Float64}(Sred))
    keep = findall(>(tol), λ)
    isempty(keep) && error("s_invsqrt_reduced: no eigenvalues > tol")
    λk = λ[keep]; Qk = Q[:, keep]
    Dk_invhalf = Diagonal(1.0 ./ sqrt.(λk))
    return Qk, Dk_invhalf
end

# comx_reduced: COMX-localize by diagonalizing X in an orthonormalized subspace.
function comx_reduced(Qk::AbstractMatrix, Dk_invhalf::Diagonal, Xred::AbstractMatrix)
    Xk = Matrix(Symmetric(Qk' * Xred * Qk))
    Xc = Matrix(Symmetric(Dk_invhalf * Xk * Dk_invhalf))
    xc, Uk = eig(Xc)
    return Uk, xc
end

function _gpoly_injection_coeff(alpha::Real, H::Real, r::AbstractVector{Float64},
        w::AbstractVector{Float64}, seedv::AbstractMatrix{Float64},
        S_seed::AbstractMatrix{Float64}; tol::Real=_BOUNDARY_S_EIG_TOL)
    (alpha > 0) || error("gpoly injection: alpha must be > 0, got $(alpha)")
    x = r ./ H
    f = x .* exp.(-0.5 .* (x ./ alpha) .^ 2)
    b = seedv' * (w .* f)
    λ, Q = eig(Matrix{Float64}(S_seed))
    keep = findall(>(tol), λ)
    isempty(keep) && error("gpoly injection: no overlap eigenvalues > tol=$(tol)")
    Qk = Q[:, keep]
    λk = λ[keep]
    c = Qk * (Diagonal(1.0 ./ λk) * (Qk' * b))
    nrm2 = dot(c, S_seed * c)
    (nrm2 > 0) || error("gpoly injection: nonpositive norm encountered")
    return c ./ sqrt(nrm2)
end

function _gpoly_injection_values(alpha::Real, H::Real, r::AbstractVector{Float64})
    (alpha > 0) || error("gpoly injection: alpha must be > 0, got $(alpha)")
    x = r ./ H
    return x .* exp.(-0.5 .* (x ./ alpha) .^ 2)
end

function _gpoly_injection_matrix(alphas::AbstractVector{Float64}, H::Real, r::AbstractVector{Float64})
    n = length(alphas)
    vals = Matrix{Float64}(undef, length(r), n)
    @inbounds for k in 1:n
        vals[:, k] = _gpoly_injection_values(alphas[k], H, r)
    end
    return vals
end

function _gpoly_eval!(vals::Vector{Float64}, dvals::Vector{Float64},
        u::Float64, H::Real, alphas::AbstractVector{Float64})
    x = u / H
    invH = 1.0 / H
    @inbounds for k in eachindex(alphas)
        a = alphas[k]
        t = x / a
        ex = exp(-0.5 * t * t)
        vals[k] = x * ex
        dvals[k] = invH * ex * (1.0 - t * t)
    end
    return nothing
end

function _collect_inject_alphas(alpha1::Real, alpha2::Real)
    alphas = Float64[]
    if isfinite(alpha1) && !isnan(alpha1)
        push!(alphas, Float64(alpha1))
        if isfinite(alpha2) && !isnan(alpha2)
            push!(alphas, Float64(alpha2))
        end
    end
    return alphas
end

function _basis_merit(C::AbstractMatrix{Float64}, xc::AbstractVector{Float64},
        seedv::AbstractMatrix{Float64}, r::AbstractVector{Float64},
        w::AbstractVector{Float64})
    phi = seedv * C
    merit = 0.0
    for k in 1:size(phi, 2)
        χ = phi[:, k]
        w0 = sum(w .* χ)
        w0 == 0.0 && continue
        xbar = sum(w .* (r .* χ)) / w0
        d = xc[k] - xbar
        merit += d * d
    end
    return merit
end

function _basis_merit_stats(C::AbstractMatrix{Float64}, xc::AbstractVector{Float64},
        seedv::AbstractMatrix{Float64}, r::AbstractVector{Float64},
        w::AbstractVector{Float64})
    phi = seedv * C
    merit = 0.0
    max_abs = 0.0
    sum_abs = 0.0
    min_abs = Inf
    min_w0 = Inf
    max_w0 = -Inf
    min_abs_w0 = Inf
    for k in 1:size(phi, 2)
        χ = phi[:, k]
        w0 = sum(w .* χ)
        if w0 == 0.0
            continue
        end
        xbar = sum(w .* (r .* χ)) / w0
        d = xc[k] - xbar
        ad = abs(d)
        merit += d * d
        max_abs = max(max_abs, ad)
        min_abs = min(min_abs, ad)
        sum_abs += ad
        min_w0 = min(min_w0, w0)
        max_w0 = max(max_w0, w0)
        min_abs_w0 = min(min_abs_w0, abs(w0))
    end
    mean_abs = sum_abs / max(1, size(phi, 2))
    return (; merit, max_abs, min_abs, mean_abs, min_w0, max_w0, min_abs_w0)
end

"""
    build_odd_even_for(; jpos_save, neven_add, grid_h, grid_sigma, grid_s0, rgridmax, ...) -> NamedTuple

Construct the default (new) boundary basis (`boundary_method=:odd_even`):

- Start from a symmetric full-line window of raw gausslets `G_j` with `j=-L:L`
  where `L=jpos_save`.
- Form odd and even linear combinations, restrict to the half-line, then:
  - COMX-localize the odd subspace (already a good radial basis for `u(r)`).
  - Project the odd subspace out of the even span (since `Θ(r)` breaks parity).
  - Also project out the “delta at the origin” direction from the even complement.
  - Take the largest-eigenvalue modes of the even-complement overlap as
    completeness fixes (`neven_add` of them).
- Orthonormalize the combined span and COMX-localize again.

Returns `(; js_trim, Ctrim, Cinj, xctrim, params)` where
- `js_trim` is the seed index window (`-L:L`),
- `Ctrim` maps seed gausslets → boundary basis functions (columns),
- `Cinj` maps injected primitives → boundary basis functions (columns; empty if none),
- `xctrim` are COMX centers in the seed coordinate `u` (used as cache metadata).

Parameter tips
- `jpos_save` (`L`): increases the number of stored boundary functions and pushes
  the far-tail completeness outwards.  This is the “size” of the boundary cache.
- `neven_add`: 0 gives odd-only; small numbers (2–8) add even-Taylor content near
  the origin but can increase cancellation if pushed too high.
- `grid_h/grid_sigma/grid_s0/rgridmax`: internal quadrature used only while
  constructing/orthonormalizing the cached boundary basis.
"""
function build_odd_even_for(;
        H::Real=_BOUNDARY_H,
        stencil = getproperty(Gausslets, _BOUNDARY_FAMILY),
        jpos_save::Int=_BOUNDARY_JPOS_SAVE,
        neven_add::Int=_BOUNDARY_ODD_EVEN_NEVEN_ADD,
        inject_count::Int=_BOUNDARY_INJECT_COUNT,
        inject_mode::Symbol=_BOUNDARY_INJECT_MODE,
        inject_search::Bool=_BOUNDARY_INJECT_SEARCH,
        inject_verbose::Bool=_BOUNDARY_INJECT_VERBOSE,
        inject_alpha1::Real=_BOUNDARY_INJECT_ALPHA1,
        inject_alpha2::Real=_BOUNDARY_INJECT_ALPHA2,
        inject_alpha1_range::Tuple{Real,Real}=_BOUNDARY_INJECT_ALPHA1_RANGE,
        inject_alpha2_range::Tuple{Real,Real}=_BOUNDARY_INJECT_ALPHA2_RANGE,
        inject_alpha1_grid_n::Int=_BOUNDARY_INJECT_ALPHA1_GRID_N,
        inject_alpha2_grid_n::Int=_BOUNDARY_INJECT_ALPHA2_GRID_N,
        inject_alpha_refine_step::Real=_BOUNDARY_INJECT_ALPHA_REFINE_STEP,
        inject_alpha_refine_tol::Real=_BOUNDARY_INJECT_ALPHA_REFINE_TOL,
        inject_alpha_refine_maxiter::Int=_BOUNDARY_INJECT_ALPHA_REFINE_MAXITER,
        inject_multi_start::Bool=_BOUNDARY_INJECT_MULTI_START,
        inject_multi_start_n::Int=_BOUNDARY_INJECT_MULTI_START_N,
        inject_nm_maxiter::Int=_BOUNDARY_INJECT_NM_MAXITER,
        grid_h::Real=_BOUNDARY_GRID_H,
        grid_sigma::Real=_BOUNDARY_GRID_SIGMA,
        grid_s0::Real=_BOUNDARY_GRID_S0,
        rgridmax::Real=_BOUNDARY_RGRIDMAX)

    L = jpos_save
    (L >= 1) || error("build_odd_even_for: jpos_save (L) must be >= 1, got $(L)")
    (neven_add >= 0) || error("build_odd_even_for: neven_add must be >= 0, got $(neven_add)")

    r, w = make_erf_grid(; h=grid_h, rmax=rgridmax, sigma=grid_sigma, s0=grid_s0)
    js = collect(-L:L)
    nseed = length(js)

    g = Gausslets.gausslet(stencil, H)
    phi = phi_eval_factory(g, H)

    _, wseed, seedv = seed_scalar_integrals(phi, js, r, w)
    S_seed = seedv' * (w .* seedv)
    X_seed = seedv' * ((r .* w) .* seedv)

    # Build raw odd/even coefficient blocks in the seed basis.
    nodd  = L
    neven = L + 1
    Bodd  = zeros(Float64, nseed, nodd)
    Beven = zeros(Float64, nseed, neven)

    # E_0 = G_0
    Beven[L + 1, 1] = 1.0
    for k in 1:L
        ip = (L + 1) + k
        im = (L + 1) - k
        Bodd[ip, k]  = 1.0
        Bodd[im, k]  = -1.0
        Beven[ip, k + 1] = 1.0
        Beven[im, k + 1] = 1.0
    end

    # Normalize columns on the half-line: ||f||^2 = c' S_seed c.
	    function normalize_columns(B::Matrix{Float64})
	        SB = S_seed * B
	        norms2 = vec(sum(B .* SB; dims=1))
	        any(<=(0.0), norms2) && error("build_odd_even_for: nonpositive column norm^2 encountered")
	        return B * Diagonal(1.0 ./ sqrt.(norms2))
	    end

    BoddN  = normalize_columns(Bodd)
    BevenN = normalize_columns(Beven)

    # Orthonormalize odd span and COMX-localize.
    Sodd = Matrix(Symmetric(BoddN' * S_seed * BoddN))
    Xodd = Matrix(Symmetric(BoddN' * X_seed * BoddN))

    Qk, Dk_invhalf = s_invsqrt_reduced(Sodd)
    Uk, xc_odd = comx_reduced(Qk, Dk_invhalf, Xodd)
    Todd = Qk * (Dk_invhalf * Uk)
    Codd = BoddN * Todd

    # If neven_add=0 and no injection, return the odd-only COMX basis.
    if neven_add == 0 && inject_count == 0
        Cfinal = Codd
        xc = xc_odd
        wchi = vec(wseed' * Cfinal)
        for m in 1:size(Cfinal, 2)
            if wchi[m] < 0.0
                Cfinal[:, m] .*= -1.0
                wchi[m] = -wchi[m]
            end
        end
        Cgauss = Cfinal
        Cinj = zeros(Float64, 0, size(Cgauss, 2))
        params = (; jpos_save, neven_add, grid_h, grid_sigma, grid_s0, rgridmax,
                  inject_count, inject_mode, inject_search,
                  alpha1=NaN, alpha2=NaN,
                  inject_alpha1_range, inject_alpha2_range,
                  inject_alpha1_grid_n, inject_alpha2_grid_n,
                  inject_alpha_refine_step, inject_alpha_refine_tol, inject_alpha_refine_maxiter,
                  inject_multi_start, inject_multi_start_n, inject_nm_maxiter)
        return (; js_trim=js, Ctrim=Cgauss, Cinj, xctrim=xc, params)
    end

    Ceven_add = Matrix{Float64}(undef, nseed, 0)
    if neven_add > 0
        # Project odd subspace out of the even span.
        Cproj = Codd' * (S_seed * BevenN)
        Ceven_ortho = BevenN .- (Codd * Cproj)

        # Project out the delta-at-origin direction in the even complement.
        seed0 = [g(-j * H) for j in js]
        v0 = vec(seed0' * Ceven_ortho)
        denom = dot(v0, v0)
        if denom > 0.0
            Pdelta = Matrix{Float64}(I, neven, neven) .- (v0 * v0') ./ denom
            Ceven_ortho = Ceven_ortho * Pdelta
        end

        Seven = Matrix(Symmetric(Ceven_ortho' * S_seed * Ceven_ortho))
        lamE, vecE = eig(Seven)
        permE = sortperm(lamE; rev=true)

        nadd = min(neven_add, length(lamE))
        for p in 1:nadd
            idx = permE[p]
            lam = lamE[idx]
            if !(lam > 0.0)
                continue
            end
            v = vecE[:, idx]
            Ce = (Ceven_ortho * v) ./ sqrt(lam)
            Ceven_add = hcat(Ceven_add, Ce)
        end
    end

    C0 = hcat(Codd, Ceven_add)
    nseed = length(js)
    function finalize_basis(Cbase::Matrix{Float64},
            S_seed::AbstractMatrix{Float64},
            X_seed::AbstractMatrix{Float64},
            wseed::AbstractVector{Float64})
        S0 = Matrix(Symmetric(Cbase' * S_seed * Cbase))
        X0 = Matrix(Symmetric(Cbase' * X_seed * Cbase))
        Q0, D0_invhalf = s_invsqrt_reduced(S0)
        U0, xc = comx_reduced(Q0, D0_invhalf, X0)
        T0 = Q0 * (D0_invhalf * U0)
        Cfinal = Cbase * T0
        wchi = vec(wseed' * Cfinal)
        for m in 1:size(Cfinal, 2)
            if wchi[m] < 0.0
                Cfinal[:, m] .*= -1.0
                wchi[m] = -wchi[m]
            end
        end
        return Cfinal, xc
    end

    inject_count < 0 && error("build_odd_even_for: inject_count must be >= 0")
    (inject_count <= 2) || error("build_odd_even_for: inject_count must be 0, 1, or 2")
    (inject_mode == :gpoly || inject_count == 0) || error("build_odd_even_for: only inject_mode=:gpoly is supported")
    alpha1_used = NaN
    alpha2_used = NaN

    function build_Cbase(alpha1::Float64, alpha2::Float64)
        if inject_count == 0
            return (Cbase=C0,
                    S_seed=S_seed,
                    X_seed=X_seed,
                    wseed=wseed,
                    seedv=seedv,
                    alphas=Float64[])
        end
        alphas = inject_count == 1 ? [alpha1] : [alpha1, alpha2]
        any(isnan, alphas) && error("build_odd_even_for: inject alphas must be finite")
        ninj = length(alphas)
        injv = _gpoly_injection_matrix(alphas, H, r)
        Sgi = seedv' * (w .* injv)
        Sii = injv' * (w .* injv)
        Xgi = seedv' * ((r .* w) .* injv)
        Xii = injv' * ((r .* w) .* injv)
        S_ext = [S_seed Sgi; Sgi' Sii]
        X_ext = [X_seed Xgi; Xgi' Xii]
        w_inj = vec(injv' * w)
        seedv_ext = hcat(seedv, injv)
        wseed_ext = vcat(wseed, w_inj)
        C0_ext = vcat(C0, zeros(Float64, ninj, size(C0, 2)))
        inj_basis = Matrix{Float64}(I, ninj, ninj)
        Cbase = hcat(C0_ext, vcat(zeros(Float64, nseed, ninj), inj_basis))
        return (Cbase=Cbase,
                S_seed=S_ext,
                X_seed=X_ext,
                wseed=wseed_ext,
                seedv=seedv_ext,
                alphas=alphas)
    end

    function merit_for(alpha1::Float64, alpha2::Float64)
        data = build_Cbase(alpha1, alpha2)
        Cfinal, xc = finalize_basis(data.Cbase, data.S_seed, data.X_seed, data.wseed)
        return _basis_merit(Cfinal, xc, data.seedv, r, w)
    end

    if inject_count == 0
        Cfinal, xc = finalize_basis(C0, S_seed, X_seed, wseed)
        Cgauss = Cfinal
        Cinj = zeros(Float64, 0, size(Cgauss, 2))
        params = (; jpos_save, neven_add, grid_h, grid_sigma, grid_s0, rgridmax,
                  inject_count, inject_mode, inject_search,
                  alpha1=NaN, alpha2=NaN,
                  inject_alpha1_range, inject_alpha2_range,
                  inject_alpha1_grid_n, inject_alpha2_grid_n,
                  inject_alpha_refine_step, inject_alpha_refine_tol, inject_alpha_refine_maxiter,
                  inject_multi_start, inject_multi_start_n, inject_nm_maxiter)
        return (; js_trim=js, Ctrim=Cgauss, Cinj, xctrim=xc, params)
    elseif inject_count == 1
        if inject_search
            a1_min, a1_max = inject_alpha1_range
            grid = collect(range(a1_min, a1_max, length=inject_alpha1_grid_n))
            merits = [merit_for(a, NaN) for a in grid]
            best_idx = argmin(merits)
            best_a = grid[best_idx]
            best_m = merits[best_idx]
            if inject_verbose
                @info "RadialGGrid: inject1 grid best alpha=$(best_a) merit=$(best_m)"
            end
            step = inject_alpha_refine_step
            iter = 0
            while step > inject_alpha_refine_tol && iter < inject_alpha_refine_maxiter
                iter += 1
                a_left = max(a1_min, best_a - step)
                a_right = min(a1_max, best_a + step)
                m_left = (a_left == best_a) ? Inf : merit_for(a_left, NaN)
                m_right = (a_right == best_a) ? Inf : merit_for(a_right, NaN)
                if m_left < best_m || m_right < best_m
                    if m_left <= m_right
                        best_a = a_left
                        best_m = m_left
                    else
                        best_a = a_right
                        best_m = m_right
                    end
                else
                    step *= 0.5
                end
            end
            alpha1_used = best_a
            if inject_verbose
                @info "RadialGGrid: inject1 refined alpha=$(alpha1_used) merit=$(best_m)"
            end
        else
            alpha1_used = Float64(inject_alpha1)
        end
        data = build_Cbase(alpha1_used, NaN)
        Cfinal, xc = finalize_basis(data.Cbase, data.S_seed, data.X_seed, data.wseed)
    else
        if !inject_search
            alpha1_used = Float64(inject_alpha1)
            alpha2_used = Float64(inject_alpha2)
            alpha2_used < alpha1_used || error("build_odd_even_for: require alpha2 < alpha1")
        else
            try
                @eval import Optim
            catch
                error("build_odd_even_for: Optim.jl is required for inject_count=2 search; install Optim or set inject_search=false.")
            end
            a1_min, a1_max = inject_alpha1_range
            a2_min, a2_max = inject_alpha2_range
            function merit_pair(x::Vector{Float64})
                a1 = clamp(x[1], a1_min, a1_max)
                a2 = clamp(x[2], a2_min, a2_max)
                if a2 >= a1
                    return 1.0e6 + 1.0e4 * (a2 - a1 + 1.0e-6)^2
                end
                return merit_for(a1, a2)
            end
            method = Optim.NelderMead()
            opts = Optim.Options(iterations=inject_nm_maxiter, show_trace=false)
            seeds = Vector{Vector{Float64}}()
            if inject_multi_start
                g1 = collect(range(a1_min, a1_max, length=inject_alpha1_grid_n))
                g2 = collect(range(a2_min, a2_max, length=inject_alpha2_grid_n))
                grid = Vector{Tuple{Float64,Float64,Float64}}()
                for a1 in g1, a2 in g2
                    a2 < a1 || continue
                    m = merit_for(a1, a2)
                    push!(grid, (m, a1, a2))
                end
                sort!(grid, by = t -> t[1])
                nsel = min(inject_multi_start_n, length(grid))
                for i in 1:nsel
                    push!(seeds, [grid[i][2], grid[i][3]])
                end
            end
            if isempty(seeds)
                push!(seeds, [Float64(inject_alpha1), Float64(inject_alpha2)])
            end
            if inject_verbose
                @info "RadialGGrid: inject2 search seeds=$(length(seeds)) a1_range=($(a1_min), $(a1_max)) a2_range=($(a2_min), $(a2_max))"
            end
            best_m = Inf
            for seed in seeds
                res = Optim.optimize(merit_pair, seed, method, opts)
                xmin = Optim.minimizer(res)
                a1 = clamp(xmin[1], a1_min, a1_max)
                a2 = clamp(xmin[2], a2_min, a2_max)
                a2 < a1 || continue
                m = Optim.minimum(res)
                if inject_verbose
                    @info "RadialGGrid: inject2 seed=$(seed) -> alpha1=$(a1) alpha2=$(a2) merit=$(m) iters=$(Optim.iterations(res)) converged=$(Optim.converged(res))"
                end
                if m < best_m
                    best_m = m
                    alpha1_used = a1
                    alpha2_used = a2
                end
            end
            isnan(alpha1_used) && error("build_odd_even_for: inject_count=2 search failed")
        end
        data = build_Cbase(alpha1_used, alpha2_used)
        Cfinal, xc = finalize_basis(data.Cbase, data.S_seed, data.X_seed, data.wseed)
    end

    nalpha = inject_count == 0 ? 0 : (inject_count == 1 ? 1 : 2)
    Cgauss = Cfinal[1:nseed, :]
    Cinj = nalpha == 0 ? zeros(Float64, 0, size(Cgauss, 2)) : Cfinal[(nseed + 1):end, :]
    params = (; jpos_save, neven_add, grid_h, grid_sigma, grid_s0, rgridmax,
              inject_count, inject_mode, inject_search,
              alpha1=alpha1_used, alpha2=alpha2_used,
              inject_alpha1_range, inject_alpha2_range,
              inject_alpha1_grid_n, inject_alpha2_grid_n,
              inject_alpha_refine_step, inject_alpha_refine_tol, inject_alpha_refine_maxiter,
              inject_multi_start, inject_multi_start_n, inject_nm_maxiter)
    if inject_verbose && inject_count > 0
        data = build_Cbase(alpha1_used, alpha2_used)
        stats = _basis_merit_stats(Cfinal, xc, data.seedv, r, w)
        @info "RadialGGrid: merit stats (base_h=$(grid_h)) merit=$(stats.merit) max|xc-xbar|=$(stats.max_abs) mean|xc-xbar|=$(stats.mean_abs) min|w0|=$(stats.min_abs_w0) max|w0|=$(stats.max_w0)"
        rchk, wchk = make_erf_grid(; h=grid_h / 2, rmax=rgridmax, sigma=grid_sigma, s0=grid_s0)
        _, _, seedv_chk = seed_scalar_integrals(phi, js, rchk, wchk)
        if inject_count > 0
            alphas = data.alphas
            injv_chk = _gpoly_injection_matrix(alphas, H, rchk)
            seedv_chk = hcat(seedv_chk, injv_chk)
        end
        stats_chk = _basis_merit_stats(Cfinal, xc, seedv_chk, rchk, wchk)
        @info "RadialGGrid: merit stats (check_h=$(grid_h / 2)) merit=$(stats_chk.merit) max|xc-xbar|=$(stats_chk.max_abs) mean|xc-xbar|=$(stats_chk.mean_abs) min|w0|=$(stats_chk.min_abs_w0) max|w0|=$(stats_chk.max_w0)"
    end
    return (; js_trim=js, Ctrim=Cgauss, Cinj, xctrim=xc, params)
end

"""
    build_odd_even_tail_for(; jpos_save, even_tail_K, ...) -> NamedTuple

Hybrid boundary basis (`boundary_method=:odd_even_tail`):

- Build and COMX-localize the full odd span (as in `build_odd_even_for`) from a
  symmetric seed window `j=-L:L` with `L=jpos_save`.
- For the even side, keep only the “near-origin tail window” `E_k` for `k=0:K`,
  where `K=even_tail_K`:
    - `E_0 = G_0`
    - `E_k = G_k + G_-k`
  then project these evens against the odd basis on the half-line and *remove*
  the delta-at-origin direction (dropping one column so the count becomes `K`).
- Combine `[odd ; even_tail_clean]`, orthonormalize, and COMX-localize again.

The motivation is to add a controlled set of even “tail-like” directions without
solving an eigenproblem in the full even-complement space.
"""
function build_odd_even_tail_for(;
        H::Real=_BOUNDARY_H,
        stencil = getproperty(Gausslets, _BOUNDARY_FAMILY),
        jpos_save::Int=_BOUNDARY_JPOS_SAVE,
        even_tail_K::Int=_BOUNDARY_ODD_EVEN_TAIL_K,
        inject_count::Int=_BOUNDARY_INJECT_COUNT,
        inject_mode::Symbol=_BOUNDARY_INJECT_MODE,
        inject_search::Bool=_BOUNDARY_INJECT_SEARCH,
        inject_verbose::Bool=_BOUNDARY_INJECT_VERBOSE,
        inject_alpha1::Real=_BOUNDARY_INJECT_ALPHA1,
        inject_alpha2::Real=_BOUNDARY_INJECT_ALPHA2,
        inject_alpha1_range::Tuple{Real,Real}=_BOUNDARY_INJECT_ALPHA1_RANGE,
        inject_alpha2_range::Tuple{Real,Real}=_BOUNDARY_INJECT_ALPHA2_RANGE,
        inject_alpha1_grid_n::Int=_BOUNDARY_INJECT_ALPHA1_GRID_N,
        inject_alpha2_grid_n::Int=_BOUNDARY_INJECT_ALPHA2_GRID_N,
        inject_alpha_refine_step::Real=_BOUNDARY_INJECT_ALPHA_REFINE_STEP,
        inject_alpha_refine_tol::Real=_BOUNDARY_INJECT_ALPHA_REFINE_TOL,
        inject_alpha_refine_maxiter::Int=_BOUNDARY_INJECT_ALPHA_REFINE_MAXITER,
        inject_multi_start::Bool=_BOUNDARY_INJECT_MULTI_START,
        inject_multi_start_n::Int=_BOUNDARY_INJECT_MULTI_START_N,
        inject_nm_maxiter::Int=_BOUNDARY_INJECT_NM_MAXITER,
        grid_h::Real=_BOUNDARY_GRID_H,
        grid_sigma::Real=_BOUNDARY_GRID_SIGMA,
        grid_s0::Real=_BOUNDARY_GRID_S0,
        rgridmax::Real=_BOUNDARY_RGRIDMAX)

    L = jpos_save
    (L >= 1) || error("build_odd_even_tail_for: jpos_save (L) must be >= 1, got $(L)")
    (0 <= even_tail_K <= L) || error("build_odd_even_tail_for: need 0 <= even_tail_K <= L, got K=$(even_tail_K) L=$(L)")

    r, w = make_erf_grid(; h=grid_h, rmax=rgridmax, sigma=grid_sigma, s0=grid_s0)
    js = collect(-L:L)
    nseed = length(js)

    g = Gausslets.gausslet(stencil, H)
    phi = phi_eval_factory(g, H)

    _, wseed, seedv = seed_scalar_integrals(phi, js, r, w)
    S_seed = seedv' * (w .* seedv)
    X_seed = seedv' * ((r .* w) .* seedv)

    # Odd block: O_k = G_k - G_-k, k=1..L
    nodd = L
    Bodd = zeros(Float64, nseed, nodd)
    for k in 1:L
        ip = (L + 1) + k
        im = (L + 1) - k
        Bodd[ip, k] = 1.0
        Bodd[im, k] = -1.0
    end

    # Even tail block: E_0, E_1..E_K
    neven0 = even_tail_K + 1
    Beven = zeros(Float64, nseed, neven0)
    Beven[L + 1, 1] = 1.0  # E_0 = G_0
    for k in 1:even_tail_K
        ip = (L + 1) + k
        im = (L + 1) - k
        Beven[ip, k + 1] = 1.0
        Beven[im, k + 1] = 1.0
    end

    function normalize_columns(B::Matrix{Float64})
        SB = S_seed * B
        norms2 = vec(sum(B .* SB; dims=1))
        any(<=(0.0), norms2) && error("build_odd_even_tail_for: nonpositive column norm^2 encountered")
        return B * Diagonal(1.0 ./ sqrt.(norms2))
    end

    BoddN = normalize_columns(Bodd)
    BevenN = normalize_columns(Beven)

    # Odd: orthonormalize and COMX-localize.
    Sodd = Matrix(Symmetric(BoddN' * S_seed * BoddN))
    Xodd = Matrix(Symmetric(BoddN' * X_seed * BoddN))

    Qk, Dk_invhalf = s_invsqrt_reduced(Sodd)
    Uk, xc_odd = comx_reduced(Qk, Dk_invhalf, Xodd)
    Todd = Qk * (Dk_invhalf * Uk)
    Codd = BoddN * Todd

    # K=0 and no injection => odd-only.
    if even_tail_K == 0 && inject_count == 0
        Cfinal = Codd
        xc = xc_odd
        wchi = vec(wseed' * Cfinal)
        for m in 1:size(Cfinal, 2)
            if wchi[m] < 0.0
                Cfinal[:, m] .*= -1.0
                wchi[m] = -wchi[m]
            end
        end
        Cgauss = Cfinal
        Cinj = zeros(Float64, 0, size(Cgauss, 2))
        params = (; jpos_save, even_tail_K, grid_h, grid_sigma, grid_s0, rgridmax,
                  inject_count, inject_mode, inject_search,
                  alpha1=NaN, alpha2=NaN,
                  inject_alpha1_range, inject_alpha2_range,
                  inject_alpha1_grid_n, inject_alpha2_grid_n,
                  inject_alpha_refine_step, inject_alpha_refine_tol, inject_alpha_refine_maxiter,
                  inject_multi_start, inject_multi_start_n, inject_nm_maxiter)
        return (; js_trim=js, Ctrim=Cgauss, Cinj, xctrim=xc, params)
    end

    Ceven_ortho = Matrix{Float64}(undef, nseed, 0)
    if even_tail_K > 0
        # Even tail: project out odd subspace (half-line metric), then remove delta-at-origin by dropping 1 column.
        Cproj = Codd' * (S_seed * BevenN)               # (nodd x neven0)
        Ceven_ortho = BevenN .- (Codd * Cproj)          # (nseed x neven0)

        # delta direction in coefficient space (values at u=0).
        seed0 = [g(-j * H) for j in js]
        v0 = vec(seed0' * Ceven_ortho)                  # length neven0
        denom = dot(v0, v0)
        if denom > 0.0
            u = v0 ./ sqrt(denom)                       # unit vector
            sgn = (u[1] >= 0.0) ? 1.0 : -1.0
            v = copy(u)
            v[1] += sgn
            beta = dot(v, v)
            beta == 0.0 && error("build_odd_even_tail_for: Householder beta==0 (unexpected)")
            Hh = Matrix{Float64}(I, neven0, neven0) .- (2.0 / beta) .* (v * v')
            Qperp = Hh[:, 2:end]                        # neven0×(neven0-1)
            Ceven_ortho = Ceven_ortho * Qperp           # nseed×K
        end
    end

    C0 = hcat(Codd, Ceven_ortho)

    function finalize_basis(Cbase::Matrix{Float64},
            seedv::AbstractMatrix{Float64},
            X_seed::AbstractMatrix{Float64},
            wseed::AbstractVector{Float64})
        # Orthonormalize the combined span in the half-line metric using a weighted
        # value-space QR. This avoids forming S0 for tail windows.
        B = seedv * Cbase
        B .*= sqrt.(w)                         # incorporate weights into Euclidean inner product
        F = qr(B)
        R = F.R                                # (ncol×ncol) upper triangular
        Cortho = Cbase / R                     # ensures Cortho' * S_seed * Cortho ≈ I
        Xortho = Matrix(Symmetric(Cortho' * X_seed * Cortho))
        xc, Ux = eig(Xortho)
        Cfinal = Cortho * Ux
        wchi = vec(wseed' * Cfinal)
        for m in 1:size(Cfinal, 2)
            if wchi[m] < 0.0
                Cfinal[:, m] .*= -1.0
                wchi[m] = -wchi[m]
            end
        end
        return Cfinal, xc
    end

    inject_count < 0 && error("build_odd_even_tail_for: inject_count must be >= 0")
    (inject_count <= 2) || error("build_odd_even_tail_for: inject_count must be 0, 1, or 2")
    (inject_mode == :gpoly || inject_count == 0) || error("build_odd_even_tail_for: only inject_mode=:gpoly is supported")
    alpha1_used = NaN
    alpha2_used = NaN

    function build_Cbase(alpha1::Float64, alpha2::Float64)
        if inject_count == 0
            return (Cbase=C0,
                    seedv=seedv,
                    X_seed=X_seed,
                    wseed=wseed,
                    alphas=Float64[])
        end
        alphas = inject_count == 1 ? [alpha1] : [alpha1, alpha2]
        any(isnan, alphas) && error("build_odd_even_tail_for: inject alphas must be finite")
        ninj = length(alphas)
        injv = _gpoly_injection_matrix(alphas, H, r)
        Xgi = seedv' * ((r .* w) .* injv)
        Xii = injv' * ((r .* w) .* injv)
        X_ext = [X_seed Xgi; Xgi' Xii]
        w_inj = vec(injv' * w)
        seedv_ext = hcat(seedv, injv)
        wseed_ext = vcat(wseed, w_inj)
        C0_ext = vcat(C0, zeros(Float64, ninj, size(C0, 2)))
        inj_basis = Matrix{Float64}(I, ninj, ninj)
        Cbase = hcat(C0_ext, vcat(zeros(Float64, nseed, ninj), inj_basis))
        return (Cbase=Cbase,
                seedv=seedv_ext,
                X_seed=X_ext,
                wseed=wseed_ext,
                alphas=alphas)
    end

    function merit_for(alpha1::Float64, alpha2::Float64)
        data = build_Cbase(alpha1, alpha2)
        Cfinal, xc = finalize_basis(data.Cbase, data.seedv, data.X_seed, data.wseed)
        return _basis_merit(Cfinal, xc, data.seedv, r, w)
    end

    if inject_count == 0
        Cfinal, xc = finalize_basis(C0, seedv, X_seed, wseed)
        Cgauss = Cfinal
        Cinj = zeros(Float64, 0, size(Cgauss, 2))
        params = (; jpos_save, even_tail_K, grid_h, grid_sigma, grid_s0, rgridmax,
                  inject_count, inject_mode, inject_search,
                  alpha1=NaN, alpha2=NaN,
                  inject_alpha1_range, inject_alpha2_range,
                  inject_alpha1_grid_n, inject_alpha2_grid_n,
                  inject_alpha_refine_step, inject_alpha_refine_tol, inject_alpha_refine_maxiter,
                  inject_multi_start, inject_multi_start_n, inject_nm_maxiter)
        return (; js_trim=js, Ctrim=Cgauss, Cinj, xctrim=xc, params)
    elseif inject_count == 1
        if inject_search
            a1_min, a1_max = inject_alpha1_range
            grid = collect(range(a1_min, a1_max, length=inject_alpha1_grid_n))
            merits = [merit_for(a, NaN) for a in grid]
            best_idx = argmin(merits)
            best_a = grid[best_idx]
            best_m = merits[best_idx]
            if inject_verbose
                @info "RadialGGrid: inject1(grid) best alpha=$(best_a) merit=$(best_m)"
            end
            step = inject_alpha_refine_step
            iter = 0
            while step > inject_alpha_refine_tol && iter < inject_alpha_refine_maxiter
                iter += 1
                a_left = max(a1_min, best_a - step)
                a_right = min(a1_max, best_a + step)
                m_left = (a_left == best_a) ? Inf : merit_for(a_left, NaN)
                m_right = (a_right == best_a) ? Inf : merit_for(a_right, NaN)
                if m_left < best_m || m_right < best_m
                    if m_left <= m_right
                        best_a = a_left
                        best_m = m_left
                    else
                        best_a = a_right
                        best_m = m_right
                    end
                else
                    step *= 0.5
                end
            end
            alpha1_used = best_a
            if inject_verbose
                @info "RadialGGrid: inject1 refined alpha=$(alpha1_used) merit=$(best_m)"
            end
        else
            alpha1_used = Float64(inject_alpha1)
        end
        data = build_Cbase(alpha1_used, NaN)
        Cfinal, xc = finalize_basis(data.Cbase, data.seedv, data.X_seed, data.wseed)
    else
        if !inject_search
            alpha1_used = Float64(inject_alpha1)
            alpha2_used = Float64(inject_alpha2)
            alpha2_used < alpha1_used || error("build_odd_even_tail_for: require alpha2 < alpha1")
        else
            try
                @eval import Optim
            catch
                error("build_odd_even_tail_for: Optim.jl is required for inject_count=2 search; install Optim or set inject_search=false.")
            end
            a1_min, a1_max = inject_alpha1_range
            a2_min, a2_max = inject_alpha2_range
            function merit_pair(x::Vector{Float64})
                a1 = clamp(x[1], a1_min, a1_max)
                a2 = clamp(x[2], a2_min, a2_max)
                if a2 >= a1
                    return 1.0e6 + 1.0e4 * (a2 - a1 + 1.0e-6)^2
                end
                return merit_for(a1, a2)
            end
            method = Optim.NelderMead()
            opts = Optim.Options(iterations=inject_nm_maxiter, show_trace=false)
            seeds = Vector{Vector{Float64}}()
            if inject_multi_start
                g1 = collect(range(a1_min, a1_max, length=inject_alpha1_grid_n))
                g2 = collect(range(a2_min, a2_max, length=inject_alpha2_grid_n))
                grid = Vector{Tuple{Float64,Float64,Float64}}()
                for a1 in g1, a2 in g2
                    a2 < a1 || continue
                    m = merit_for(a1, a2)
                    push!(grid, (m, a1, a2))
                end
                sort!(grid, by = t -> t[1])
                nsel = min(inject_multi_start_n, length(grid))
                for i in 1:nsel
                    push!(seeds, [grid[i][2], grid[i][3]])
                end
            end
            if isempty(seeds)
                push!(seeds, [Float64(inject_alpha1), Float64(inject_alpha2)])
            end
            if inject_verbose
                @info "RadialGGrid: inject2 search seeds=$(length(seeds)) a1_range=($(a1_min), $(a1_max)) a2_range=($(a2_min), $(a2_max))"
            end
            best_m = Inf
            for seed in seeds
                res = Optim.optimize(merit_pair, seed, method, opts)
                xmin = Optim.minimizer(res)
                a1 = clamp(xmin[1], a1_min, a1_max)
                a2 = clamp(xmin[2], a2_min, a2_max)
                a2 < a1 || continue
                m = Optim.minimum(res)
                if inject_verbose
                    @info "RadialGGrid: inject2 seed=$(seed) -> alpha1=$(a1) alpha2=$(a2) merit=$(m) iters=$(Optim.iterations(res)) converged=$(Optim.converged(res))"
                end
                if m < best_m
                    best_m = m
                    alpha1_used = a1
                    alpha2_used = a2
                end
            end
            isnan(alpha1_used) && error("build_odd_even_tail_for: inject_count=2 search failed")
        end
        data = build_Cbase(alpha1_used, alpha2_used)
        Cfinal, xc = finalize_basis(data.Cbase, data.seedv, data.X_seed, data.wseed)
    end

    nalpha = inject_count == 0 ? 0 : (inject_count == 1 ? 1 : 2)
    Cgauss = Cfinal[1:nseed, :]
    Cinj = nalpha == 0 ? zeros(Float64, 0, size(Cgauss, 2)) : Cfinal[(nseed + 1):end, :]
    params = (; jpos_save, even_tail_K, grid_h, grid_sigma, grid_s0, rgridmax,
              inject_count, inject_mode, inject_search,
              alpha1=alpha1_used, alpha2=alpha2_used,
              inject_alpha1_range, inject_alpha2_range,
              inject_alpha1_grid_n, inject_alpha2_grid_n,
              inject_alpha_refine_step, inject_alpha_refine_tol, inject_alpha_refine_maxiter,
              inject_multi_start, inject_multi_start_n, inject_nm_maxiter)
    if inject_verbose && inject_count > 0
        data = build_Cbase(alpha1_used, alpha2_used)
        stats = _basis_merit_stats(Cfinal, xc, data.seedv, r, w)
        @info "RadialGGrid: merit stats (base_h=$(grid_h)) merit=$(stats.merit) max|xc-xbar|=$(stats.max_abs) mean|xc-xbar|=$(stats.mean_abs) min|w0|=$(stats.min_abs_w0) max|w0|=$(stats.max_w0)"
        rchk, wchk = make_erf_grid(; h=grid_h / 2, rmax=rgridmax, sigma=grid_sigma, s0=grid_s0)
        _, _, seedv_chk = seed_scalar_integrals(phi, js, rchk, wchk)
        if inject_count > 0
            alphas = data.alphas
            injv_chk = _gpoly_injection_matrix(alphas, H, rchk)
            seedv_chk = hcat(seedv_chk, injv_chk)
        end
        stats_chk = _basis_merit_stats(Cfinal, xc, seedv_chk, rchk, wchk)
        @info "RadialGGrid: merit stats (check_h=$(grid_h / 2)) merit=$(stats_chk.merit) max|xc-xbar|=$(stats_chk.max_abs) mean|xc-xbar|=$(stats_chk.mean_abs) min|w0|=$(stats_chk.min_abs_w0) max|w0|=$(stats_chk.max_w0)"
    end
    return (; js_trim=js, Ctrim=Cgauss, Cinj, xctrim=xc, params)
end

"""
    build_and_trim_for(jneg; jpos_build, jpos_save, ...) -> NamedTuple

Legacy boundary cache builder (`boundary_method=:add_tails`):

- Build seed gausslets on the half-line with indices `j=-jneg:jpos_build`.
- Project out the pivoted “delta at the origin” direction.
- Orthonormalize and COMX-localize.
- Save a trimmed set keeping indices `j=-jneg:jpos_save`.

This exists mainly for backward comparisons; new work should prefer
`build_odd_even_for` and `boundary_method=:odd_even`.
"""
function build_and_trim_for(jneg::Int;
        H::Real=_BOUNDARY_H,
        stencil = getproperty(Gausslets, _BOUNDARY_FAMILY),
        jpos_build::Int=_BOUNDARY_JPOS_BUILD,
        jpos_save::Int=_BOUNDARY_JPOS_SAVE,
        grid_h::Real=_BOUNDARY_GRID_H,
        grid_sigma::Real=_BOUNDARY_GRID_SIGMA,
        grid_s0::Real=_BOUNDARY_GRID_S0,
        rgridmax::Real=_BOUNDARY_RGRIDMAX)

    r, w = make_erf_grid(; h=grid_h, rmax=rgridmax, sigma=grid_sigma, s0=grid_s0)
    js_full = collect(-jneg:jpos_build)
    g   = Gausslets.gausslet(stencil, H)
    phi = phi_eval_factory(g, H)

    sdiag_raw, _, seedv = seed_scalar_integrals(phi, js_full, r, w)
    D = 1.0 ./ sqrt.(sdiag_raw)

    _, _, Tproj = pivot_projection_transform(g, js_full, D, H)

    S_seed, X_seed = assemble_halfline_seed(phi, js_full, r, w; D=D)

    Sred = Matrix(Symmetric(Tproj' * S_seed * Tproj))
    Xred = Matrix(Symmetric(Tproj' * X_seed * Tproj))

    Qk, Dk_invhalf = s_invsqrt_reduced(Sred)
    Uk, xc = comx_reduced(Qk, Dk_invhalf, Xred)

    Cfull = D .* (Tproj * (Qk * (Dk_invhalf * Uk)))

    Cseed = seedv * Cfull
    for m in 1:size(Cfull, 2)
        if sum(Cseed[:, m] .* w) < 0.0
            Cfull[:, m] .*= -1.0
        end
    end

    js_trim = collect(-jneg:jpos_save)
    rowmask = [j in js_trim for j in js_full]
    Crows   = Cfull[rowmask, :]

    nm_trim = length(js_trim) - 1
    Ctrim   = Crows[:, 1:nm_trim]
    xctrim  = xc[1:nm_trim]

    return (; js_trim, Ctrim, xctrim, params=(; jneg, jpos_save, grid_h, grid_sigma, grid_s0, rgridmax))
end

# _eval_gaussian_combo_big!: high-precision evaluation of a Gaussian expansion and its du-derivative.
function _eval_gaussian_combo_big!(
        vals::Vector{Float64},
        dvals::Vector{Float64},
        u::Float64,
        G::Matrix{Float64},
        kmin::Int,
        sp::Float64,
        wg::Float64,
        ev_buf::Vector{BigFloat},
        dev_buf::Vector{BigFloat};
        win_half::Int=12)

    N = size(G, 2)
    length(vals) == N || resize!(vals, N)
    length(dvals) == N || resize!(dvals, N)

    fill!(vals, 0.0)
    fill!(dvals, 0.0)

    kmax = kmin + size(G, 1) - 1

    ub = BigFloat(u)
    spb = BigFloat(sp)
    wgb = BigFloat(wg)
    inv_w2sp = inv(wgb * wgb * spb)

    t = ub / spb
    ii = Int(round(t))
    klo = max(ii - win_half, kmin)
    khi = min(ii + win_half, kmax)
    nk = khi - klo + 1
    row_lo = klo - kmin + 1

    k = klo
    @inbounds for idx in 1:nk
        dk = t - BigFloat(k)
        ev = exp(-BigFloat(0.5) * (dk / wgb)^2)
        ev_buf[idx] = ev
        dev_buf[idx] = ev * (-(dk) * inv_w2sp)
        k += 1
    end

    @inbounds for a in 1:N
        s0 = BigFloat(0.0)
        s1 = BigFloat(0.0)
        row = row_lo
        for idx in 1:nk
            ga = BigFloat(G[row, a])
            s0 += ga * ev_buf[idx]
            s1 += ga * dev_buf[idx]
            row += 1
        end
        vals[a] = Float64(s0)
        dvals[a] = Float64(s1)
    end
    return nothing
end

# _build_boundary_cheb_representation: fit a piecewise Chebyshev representation χ(u) and dχ/du for cached boundary functions.
function _build_boundary_cheb_representation(C::Matrix{Float64}, js::Vector{Int}, g;
        deg::Int=_BOUNDARY_CHEB_DEG,
        check_deg::Int=_BOUNDARY_CHEB_CHECK_DEG,
        u_max::Real=(maximum(js) + _BOUNDARY_CHEB_UPAD),
        tol_rel::Real=_BOUNDARY_CHEB_TOL_REL,
        tol_rel_dphi::Real=_BOUNDARY_CHEB_TOL_REL_DPHI,
        big_prec::Int=_BOUNDARY_CHEB_BIG_PREC,
        maxseg::Int=_BOUNDARY_CHEB_MAXSEG,
        minwidth::Real=_BOUNDARY_CHEB_MINWIDTH)

    (deg >= 2) || error("Chebyshev degree must be >= 2, got $(deg)")
    (check_deg >= deg) || error("check_deg must be >= deg")
    (u_max > 0) || error("u_max must be > 0")

    G, kmin = gausslet_to_gaussian_coeffs(C, js, g)
    sp = g.sp
    wg = g.w
    Nfunc = size(C, 2)

    x_fit = cheb_lobatto_nodes(deg, Float64)
    x_chk = cheb_lobatto_nodes(check_deg, Float64)

    maxabs_phi = zeros(Float64, Nfunc)
    maxabs_dphi = zeros(Float64, Nfunc)

    ev_buf = Vector{BigFloat}(undef, 2 * 12 + 1)
    dev_buf = Vector{BigFloat}(undef, 2 * 12 + 1)
    phi_tmp = zeros(Float64, Nfunc)
    dphi_tmp = zeros(Float64, Nfunc)

    return setprecision(BigFloat, big_prec) do
        function eval_truth!(phi_out::Vector{Float64}, dphi_out::Vector{Float64}, u::Float64)
            _eval_gaussian_combo_big!(phi_out, dphi_out, u, G, kmin, sp, wg, ev_buf, dev_buf)
            @inbounds for a in 1:Nfunc
                v = abs(phi_out[a])
                if v > maxabs_phi[a]
                    maxabs_phi[a] = v
                end
                dv = abs(dphi_out[a])
                if dv > maxabs_dphi[a]
                    maxabs_dphi[a] = dv
                end
            end
            return nothing
        end

        # Prepass: estimate global max amplitudes so early-tail segments do not
        # over-split due to tiny "max so far" values for far-out basis functions.
        du_pre = 0.1
        u_last = Float64(u_max)
        for u in 0.0:du_pre:u_last
            eval_truth!(phi_tmp, dphi_tmp, u)
        end

        coeffs_list = Vector{Matrix{Float64}}()
        dcoeffs_list = Vector{Matrix{Float64}}()
        edges = Float64[0.0]
        naccepted = Ref(0)
        nforced = Ref(0)
        forced_maxerr_phi = Ref(0.0)
        forced_maxerr_dphi = Ref(0.0)

        function fit_segment(a::Float64, b::Float64)
            width = b - a
            mid = 0.5 * (a + b)
            half = 0.5 * width

            vals = Matrix{Float64}(undef, deg + 1, Nfunc)
            dvals = Matrix{Float64}(undef, deg + 1, Nfunc)
            @inbounds for j in 1:(deg + 1)
                u = mid + half * x_fit[j]
                eval_truth!(phi_tmp, dphi_tmp, u)
                for aidx in 1:Nfunc
                    vals[j, aidx] = phi_tmp[aidx]
                    dvals[j, aidx] = dphi_tmp[aidx]
                end
            end

            cseg = Matrix{Float64}(undef, deg + 1, Nfunc)
            dcseg = Matrix{Float64}(undef, deg + 1, Nfunc)
            @inbounds for aidx in 1:Nfunc
                c = cheb_coeffs_lobatto(view(vals, :, aidx))
                cseg[:, aidx] = c
                dcseg[:, aidx] = cheb_coeffs_lobatto(view(dvals, :, aidx))
            end

            max_err_phi = 0.0
            max_err_dphi = 0.0
            ok = true

            @inbounds for j in 1:(check_deg + 1)
                u = mid + half * x_chk[j]
                eval_truth!(phi_tmp, dphi_tmp, u)
                x = x_chk[j]
                for aidx in 1:Nfunc
                    p = cheb_eval(view(cseg, :, aidx), x)
                    dpdu = cheb_eval(view(dcseg, :, aidx), x)
                    err = abs(p - phi_tmp[aidx])
                    derr = abs(dpdu - dphi_tmp[aidx])
                    if err > max_err_phi
                        max_err_phi = err
                    end
                    if derr > max_err_dphi
                        max_err_dphi = derr
                    end
                    scale = maxabs_phi[aidx]
                    dscale = maxabs_dphi[aidx]
                    atol = tol_rel * (scale == 0.0 ? 1.0 : scale)
                    datol = tol_rel_dphi * (dscale == 0.0 ? 1.0 : dscale)
                    if err > atol || derr > datol
                        ok = false
                        break
                    end
                end
                ok || break
            end

            if ok || width <= minwidth || (length(edges) - 1) >= maxseg
                push!(coeffs_list, cseg)
                push!(dcoeffs_list, dcseg)
                push!(edges, b)
                naccepted[] += 1
                if !ok
                    nforced[] += 1
                    forced_maxerr_phi[] = max(forced_maxerr_phi[], max_err_phi)
                    forced_maxerr_dphi[] = max(forced_maxerr_dphi[], max_err_dphi)
                end
                if naccepted[] == 1 || naccepted[] % 50 == 0
                    @info "RadialGGrid: cheb fit accepted segments=$(naccepted[]) last=[ $(round(a; digits=4)), $(round(b; digits=4)) ] maxerr=(φ=$(max_err_phi), dφ=$(max_err_dphi))"
                end
                return
            end

            m = 0.5 * (a + b)
            fit_segment(a, m)
            fit_segment(m, b)
            return
        end

        fit_segment(0.0, Float64(u_max))

        nseg = length(edges) - 1
        coeffs = Array{Float64,3}(undef, deg + 1, nseg, Nfunc)
        dcoeffs = Array{Float64,3}(undef, deg + 1, nseg, Nfunc)
        @inbounds for seg in 1:nseg
            cseg = coeffs_list[seg]
            dcseg = dcoeffs_list[seg]
            for aidx in 1:Nfunc
                coeffs[:, seg, aidx] = cseg[:, aidx]
                dcoeffs[:, seg, aidx] = dcseg[:, aidx]
            end
        end

        if nforced[] > 0
            @warn "RadialGGrid: cheb fit forced segments=$(nforced[]) (minwidth=$(minwidth), maxseg=$(maxseg)); worst maxerr=(φ=$(forced_maxerr_phi[]), dφ=$(forced_maxerr_dphi[]))"
        end

        return (; edges, deg, coeffs, dcoeffs, maxabs_phi, maxabs_dphi, u_max=Float64(u_max))
    end
end

# _boundary_overlap_deviation: diagnostic for cache-build quadrature accuracy (||S-I||∞ on a check grid).
function _boundary_overlap_deviation(C::Matrix{Float64}, js::Vector{Int}, g, H;
        grid_h::Real,
        grid_sigma::Real,
        grid_s0::Real,
        rgridmax::Real,
        Cinj::Union{Nothing,Matrix{Float64}}=nothing,
        inj_alphas::AbstractVector{Float64}=Float64[])

    r, w = make_erf_grid(; h=grid_h, rmax=rgridmax, sigma=grid_sigma, s0=grid_s0)
    phi = phi_eval_factory(g, H)
    _, _, seedv = seed_scalar_integrals(phi, js, r, w)
    B = seedv * C
    if Cinj !== nothing && !isempty(inj_alphas)
        injv = _gpoly_injection_matrix(collect(inj_alphas), H, r)
        B .+= injv * Cinj
    end
    S = B' * (w .* B)
    devS = norm(S - I, Inf)
    return devS
end

# _refine_boundary_grid_h: choose an internal cache-build grid_h by checking ||S-I||∞ on a halved grid.
function _refine_boundary_grid_h(;
        H::Real,
        stencil,
        jpos_save::Int,
        neven_add::Int,
        even_tail_K::Int=_BOUNDARY_ODD_EVEN_TAIL_K,
        inject_count::Int=_BOUNDARY_INJECT_COUNT,
        inject_mode::Symbol=_BOUNDARY_INJECT_MODE,
        inject_search::Bool=_BOUNDARY_INJECT_SEARCH,
        inject_alpha1::Real=_BOUNDARY_INJECT_ALPHA1,
        inject_alpha2::Real=_BOUNDARY_INJECT_ALPHA2,
        inject_alpha1_range::Tuple{Real,Real}=_BOUNDARY_INJECT_ALPHA1_RANGE,
        inject_alpha2_range::Tuple{Real,Real}=_BOUNDARY_INJECT_ALPHA2_RANGE,
        inject_alpha1_grid_n::Int=_BOUNDARY_INJECT_ALPHA1_GRID_N,
        inject_alpha2_grid_n::Int=_BOUNDARY_INJECT_ALPHA2_GRID_N,
        inject_alpha_refine_step::Real=_BOUNDARY_INJECT_ALPHA_REFINE_STEP,
        inject_alpha_refine_tol::Real=_BOUNDARY_INJECT_ALPHA_REFINE_TOL,
        inject_alpha_refine_maxiter::Int=_BOUNDARY_INJECT_ALPHA_REFINE_MAXITER,
        inject_multi_start::Bool=_BOUNDARY_INJECT_MULTI_START,
        inject_multi_start_n::Int=_BOUNDARY_INJECT_MULTI_START_N,
        inject_nm_maxiter::Int=_BOUNDARY_INJECT_NM_MAXITER,
        grid_h::Real,
        grid_sigma::Real,
        grid_s0::Real,
        rgridmax::Real,
        tol::Real=1e-10,
        maxiter::Int=6)

    g = Gausslets.gausslet(stencil, H)

    best_h = Float64(grid_h)
    best_dev = Inf
    h_try = Float64(grid_h)

    @info "RadialGGrid: refining boundary cache grid_h (start=$(grid_h), tol=$(tol))"
    for iter in 1:maxiter
        odd = build_odd_even_for(; H,
            stencil,
            jpos_save,
            neven_add,
            inject_count,
            inject_mode,
            inject_search=false,
            inject_alpha1,
            inject_alpha2,
            inject_alpha1_range,
            inject_alpha2_range,
            inject_alpha1_grid_n,
            inject_alpha2_grid_n,
            inject_alpha_refine_step,
            inject_alpha_refine_tol,
            inject_alpha_refine_maxiter,
            inject_multi_start,
            inject_multi_start_n,
            inject_nm_maxiter,
            grid_h=h_try,
            grid_sigma,
            grid_s0,
            rgridmax)

        tail = build_odd_even_tail_for(; H,
            stencil,
            jpos_save,
            even_tail_K,
            inject_count,
            inject_mode,
            inject_search=false,
            inject_alpha1,
            inject_alpha2,
            inject_alpha1_range,
            inject_alpha2_range,
            inject_alpha1_grid_n,
            inject_alpha2_grid_n,
            inject_alpha_refine_step,
            inject_alpha_refine_tol,
            inject_alpha_refine_maxiter,
            inject_multi_start,
            inject_multi_start_n,
            inject_nm_maxiter,
            grid_h=h_try,
            grid_sigma,
            grid_s0,
            rgridmax)

        h_check = h_try / 2
        dev_odd = _boundary_overlap_deviation(odd.Ctrim, odd.js_trim, g, H;
            grid_h=h_check,
            grid_sigma,
            grid_s0,
            rgridmax,
            Cinj=odd.Cinj,
            inj_alphas=_collect_inject_alphas(odd.params.alpha1, odd.params.alpha2))
        dev_tail = _boundary_overlap_deviation(tail.Ctrim, tail.js_trim, g, H;
            grid_h=h_check,
            grid_sigma,
            grid_s0,
            rgridmax,
            Cinj=tail.Cinj,
            inj_alphas=_collect_inject_alphas(tail.params.alpha1, tail.params.alpha2))
        dev = max(dev_odd, dev_tail)

        if dev < best_dev
            best_dev = dev
            best_h = h_try
        end
        @info "RadialGGrid: grid_h=$(h_try) check_h=$(h_check) devS=max(odd=$(dev_odd), tail=$(dev_tail))=$(dev)"

        dev <= tol && return h_try
        h_try /= 2
    end
    @warn "RadialGGrid: grid_h refinement hit maxiter=$(maxiter); best devS=$(best_dev) at grid_h=$(best_h)"
    return best_h
end

"""
    build_boundary_cache!(filename; kwargs...)

Build / overwrite the boundary-gausslet cache JLD2 file.

This is a “one-time (or occasional)” preprocessing step. The runtime builders
(`build_radial_and_one_electron`, `build_Veel_tables`) will call
`ensure_boundary_cache!` to create/refresh the needed cache groups.

Cache contents
- For `boundary_method=:odd_even`:
  - Group name: `odd_even_L{jpos_save}_e{neven_add}`
  - Stored datasets: `js`, `C`, `C_inj`, `centers`, and `cheb/*` (piecewise Chebyshev
    representation of χ(u) and dχ/du used for fast/accurate evaluation).
- For `boundary_method=:odd_even_tail`:
  - Group name: `odd_even_tail_L{jpos_save}_K{even_tail_K}`
  - Stored datasets: `js`, `C`, `C_inj`, `centers`, and `cheb/*` (same format as `:odd_even`).
- For `boundary_method=:add_tails` (legacy):
  - Group name: `jneg_{jneg}` for each `jneg` in `jneg_list`.

Key parameters / tips
- `jpos_save`: size of the symmetric seed window for the odd/even build; larger
  pushes completeness further out and increases cache size/time.
- `neven_add`: how many even-complement modes to add back; small values (2–8) are
  usually the sweet spot before cancellation grows.
- `even_tail_K`: (only for `:odd_even_tail`) maximum even index K to keep in the
  “tail” window (keeps E_k for k=0:K, then removes the delta direction leaving K).
- `grid_h/grid_sigma/grid_s0/rgridmax`: internal half-line quadrature used for
  orthonormalization/COMX while *building* the cache (independent of runtime
  `hgrid` used later in physics calculations).
- `inject_count`: number of injected gpoly primitives to add during cache build
  (0, 1, or 2). When 2, a Nelder-Mead search is run to choose `(alpha1, alpha2)`.
- `refine_grid_h=true`: automatically halves `grid_h` until the cache’s own
  overlap check meets `grid_h_tol` (helps avoid silent regressions).
- `cheb_*`: controls the piecewise Chebyshev fit. `cheb_deg=64` is a good
  default; tolerances around `1e-14` target ~14–15 digits pointwise relative to
  each function’s max.
"""
function build_boundary_cache!(filename::AbstractString;
        jneg_list=_BOUNDARY_JNEG_LIST,
        family::Symbol=_BOUNDARY_FAMILY,
        H::Real=_BOUNDARY_H,
        jpos_build::Int=_BOUNDARY_JPOS_BUILD,
        jpos_save::Int=_BOUNDARY_JPOS_SAVE,
        neven_add::Int=_BOUNDARY_ODD_EVEN_NEVEN_ADD,
        even_tail_K::Int=_BOUNDARY_ODD_EVEN_TAIL_K,
        inject_count::Int=_BOUNDARY_INJECT_COUNT,
        inject_mode::Symbol=_BOUNDARY_INJECT_MODE,
        inject_search::Bool=_BOUNDARY_INJECT_SEARCH,
        inject_verbose::Bool=_BOUNDARY_INJECT_VERBOSE,
        inject_alpha1::Real=_BOUNDARY_INJECT_ALPHA1,
        inject_alpha2::Real=_BOUNDARY_INJECT_ALPHA2,
        inject_alpha1_range::Tuple{Real,Real}=_BOUNDARY_INJECT_ALPHA1_RANGE,
        inject_alpha2_range::Tuple{Real,Real}=_BOUNDARY_INJECT_ALPHA2_RANGE,
        inject_alpha1_grid_n::Int=_BOUNDARY_INJECT_ALPHA1_GRID_N,
        inject_alpha2_grid_n::Int=_BOUNDARY_INJECT_ALPHA2_GRID_N,
        inject_alpha_refine_step::Real=_BOUNDARY_INJECT_ALPHA_REFINE_STEP,
        inject_alpha_refine_tol::Real=_BOUNDARY_INJECT_ALPHA_REFINE_TOL,
        inject_alpha_refine_maxiter::Int=_BOUNDARY_INJECT_ALPHA_REFINE_MAXITER,
        inject_multi_start::Bool=_BOUNDARY_INJECT_MULTI_START,
        inject_multi_start_n::Int=_BOUNDARY_INJECT_MULTI_START_N,
        inject_nm_maxiter::Int=_BOUNDARY_INJECT_NM_MAXITER,
        grid_h::Real=_BOUNDARY_GRID_H,
        grid_sigma::Real=_BOUNDARY_GRID_SIGMA,
        grid_s0::Real=_BOUNDARY_GRID_S0,
        rgridmax::Real=_BOUNDARY_RGRIDMAX,
        refine_grid_h::Bool=true,
        grid_h_tol::Real=1e-10,
        grid_h_maxiter::Int=6,
        cheb_deg::Int=_BOUNDARY_CHEB_DEG,
        cheb_u_pad::Real=_BOUNDARY_CHEB_UPAD,
        cheb_tol_rel::Real=_BOUNDARY_CHEB_TOL_REL,
        cheb_tol_rel_dphi::Real=_BOUNDARY_CHEB_TOL_REL_DPHI,
        cheb_big_prec::Int=_BOUNDARY_CHEB_BIG_PREC,
        cheb_maxseg::Int=_BOUNDARY_CHEB_MAXSEG,
        cheb_minwidth::Real=_BOUNDARY_CHEB_MINWIDTH,
        cheb_check_deg::Int=_BOUNDARY_CHEB_CHECK_DEG)

    if inject_count > 0 && inject_search
        @info "RadialGGrid: injection search enabled (count=$(inject_count), mode=$(inject_mode), alpha1=$(inject_alpha1), alpha2=$(inject_alpha2))"
    end

    stencil = getproperty(Gausslets, family)

    grid_h_use = Float64(grid_h)
    if refine_grid_h
        grid_h_use = _refine_boundary_grid_h(; H, stencil, jpos_save, neven_add,
            even_tail_K,
            inject_count,
            inject_mode,
            inject_search,
            inject_alpha1,
            inject_alpha2,
            inject_alpha1_range,
            inject_alpha2_range,
            inject_alpha1_grid_n,
            inject_alpha2_grid_n,
            inject_alpha_refine_step,
            inject_alpha_refine_tol,
            inject_alpha_refine_maxiter,
            inject_multi_start,
            inject_multi_start_n,
            inject_nm_maxiter,
            grid_h=grid_h_use, grid_sigma, grid_s0, rgridmax,
            tol=grid_h_tol, maxiter=grid_h_maxiter)
    end

    @info "RadialGGrid: building boundary cache → $(filename)"
	    jldopen(filename, "w") do f
        f["meta/H"]        = H
        f["meta/family"]   = String(family)
	        f["meta/desc"]     = "Boundary Gausslet C maps on [0,∞) for $(family); grid: erf-switch g(s)=0.5*erfc(s0 - s/sigma) with (sigma=$(grid_sigma), s0=$(grid_s0)), h=$(grid_h_use). Columns are COMX modes ordered by increasing <r>."
        f["meta/JPOS_SAVE"] = jpos_save
        f["meta/default_method"] = String(_BOUNDARY_DEFAULT_METHOD)
        f["meta/odd_even/neven_add"] = neven_add
        f["meta/odd_even/jpos_save"] = jpos_save
        f["meta/odd_even_tail/K"] = even_tail_K
        f["meta/odd_even_tail/jpos_save"] = jpos_save
        f["meta/odd_even_tail/version"] = _BOUNDARY_ODD_EVEN_TAIL_VERSION
        f["meta/grid_h"] = grid_h_use
        f["meta/grid_h_refined"] = refine_grid_h
        f["meta/cheb/deg"] = cheb_deg
        f["meta/cheb/u_pad"] = cheb_u_pad
        f["meta/cheb/tol_rel"] = cheb_tol_rel
        f["meta/cheb/tol_rel_dphi"] = cheb_tol_rel_dphi
        f["meta/cheb/big_prec"] = cheb_big_prec
        f["meta/cheb/maxseg"] = cheb_maxseg
        f["meta/cheb/minwidth"] = cheb_minwidth
        f["meta/cheb/check_deg"] = cheb_check_deg
        f["meta/cheb/version"] = _BOUNDARY_CHEB_VERSION

        odd = build_odd_even_for(; H,
                                 stencil,
                                 jpos_save,
                                 neven_add,
                                 inject_count,
                                 inject_mode,
                                 inject_search,
                                 inject_verbose,
                                 inject_alpha1,
                                 inject_alpha2,
                                 inject_alpha1_range,
                                 inject_alpha2_range,
                                 inject_alpha1_grid_n,
                                 inject_alpha2_grid_n,
                                 inject_alpha_refine_step,
                                 inject_alpha_refine_tol,
                                 inject_alpha_refine_maxiter,
                                 inject_multi_start,
                                 inject_multi_start_n,
                                 inject_nm_maxiter,
                                 grid_h=grid_h_use,
                                 grid_sigma,
                                 grid_s0,
                                 rgridmax)
        if inject_count > 0 && inject_search
            @info "RadialGGrid: inject result (odd_even) alpha1=$(odd.params.alpha1) alpha2=$(odd.params.alpha2)"
        end
        grp = boundary_cache_group(:odd_even; jneg=last(_BOUNDARY_JNEG_LIST), jpos_save, neven_add, even_tail_K)
        f["$grp/js"]      = odd.js_trim
        f["$grp/C"]       = odd.Ctrim
        f["$grp/C_inj"]   = odd.Cinj
        f["$grp/inject/alphas"] = _collect_inject_alphas(odd.params.alpha1, odd.params.alpha2)
        f["$grp/centers"] = odd.xctrim
        f["$grp/params"]  = string(odd.params)
        _write_inject_meta!(f, grp;
            inject_count,
            inject_mode,
            inject_search,
            alpha1=odd.params.alpha1,
            alpha2=odd.params.alpha2,
            inject_alpha1_range,
            inject_alpha2_range,
            inject_alpha1_grid_n,
            inject_alpha2_grid_n,
            inject_alpha_refine_step,
            inject_alpha_refine_tol,
            inject_alpha_refine_maxiter,
            inject_multi_start,
            inject_multi_start_n,
            inject_nm_maxiter)

        g = Gausslets.gausslet(stencil, H)
        cheb_u_max = jpos_save + cheb_u_pad
        @info "RadialGGrid: building Chebyshev cache for group=$(grp) (deg=$(cheb_deg), u_max=$(cheb_u_max), tol_rel=$(cheb_tol_rel), tol_rel_dphi=$(cheb_tol_rel_dphi))"
        t0 = time_ns()
        cheb = _build_boundary_cheb_representation(odd.Ctrim, odd.js_trim, g;
            deg=cheb_deg,
            check_deg=cheb_check_deg,
            u_max=cheb_u_max,
            tol_rel=cheb_tol_rel,
            tol_rel_dphi=cheb_tol_rel_dphi,
            big_prec=cheb_big_prec,
            maxseg=cheb_maxseg,
            minwidth=cheb_minwidth)
        dt = (time_ns() - t0) * 1e-9
        @info "RadialGGrid: Chebyshev cache built (segments=$(length(cheb.edges)-1)) in $(round(dt; digits=2)) s"
        f["$grp/cheb/edges"] = cheb.edges
        f["$grp/cheb/deg"] = cheb.deg
        f["$grp/cheb/u_max"] = cheb.u_max
        f["$grp/cheb/coeffs"] = cheb.coeffs
        f["$grp/cheb/dcoeffs"] = cheb.dcoeffs
        f["$grp/cheb/dcoeffs_kind"] = "dpdu"
        f["$grp/cheb/version"] = _BOUNDARY_CHEB_VERSION
        f["$grp/cheb/maxabs_phi"] = cheb.maxabs_phi
        f["$grp/cheb/maxabs_dphi"] = cheb.maxabs_dphi

        tail = build_odd_even_tail_for(; H,
                                       stencil,
                                       jpos_save,
                                       even_tail_K,
                                       inject_count,
                                       inject_mode,
                                       inject_search,
                                       inject_verbose,
                                       inject_alpha1,
                                       inject_alpha2,
                                       inject_alpha1_range,
                                       inject_alpha2_range,
                                       inject_alpha1_grid_n,
                                       inject_alpha2_grid_n,
                                       inject_alpha_refine_step,
                                       inject_alpha_refine_tol,
                                       inject_alpha_refine_maxiter,
                                       inject_multi_start,
                                       inject_multi_start_n,
                                       inject_nm_maxiter,
                                       grid_h=grid_h_use,
                                       grid_sigma,
                                       grid_s0,
                                       rgridmax)
        if inject_count > 0 && inject_search
            @info "RadialGGrid: inject result (odd_even_tail) alpha1=$(tail.params.alpha1) alpha2=$(tail.params.alpha2)"
        end
        grp = boundary_cache_group(:odd_even_tail; jneg=last(_BOUNDARY_JNEG_LIST), jpos_save, neven_add, even_tail_K)
        f["$grp/js"]      = tail.js_trim
        f["$grp/C"]       = tail.Ctrim
        f["$grp/C_inj"]   = tail.Cinj
        f["$grp/inject/alphas"] = _collect_inject_alphas(tail.params.alpha1, tail.params.alpha2)
        f["$grp/centers"] = tail.xctrim
        f["$grp/params"]  = string(tail.params)
        _write_inject_meta!(f, grp;
            inject_count,
            inject_mode,
            inject_search,
            alpha1=tail.params.alpha1,
            alpha2=tail.params.alpha2,
            inject_alpha1_range,
            inject_alpha2_range,
            inject_alpha1_grid_n,
            inject_alpha2_grid_n,
            inject_alpha_refine_step,
            inject_alpha_refine_tol,
            inject_alpha_refine_maxiter,
            inject_multi_start,
            inject_multi_start_n,
            inject_nm_maxiter)

        @info "RadialGGrid: building Chebyshev cache for group=$(grp) (deg=$(cheb_deg), u_max=$(cheb_u_max), tol_rel=$(cheb_tol_rel), tol_rel_dphi=$(cheb_tol_rel_dphi))"
        t0 = time_ns()
        cheb = _build_boundary_cheb_representation(tail.Ctrim, tail.js_trim, g;
            deg=cheb_deg,
            check_deg=cheb_check_deg,
            u_max=cheb_u_max,
            tol_rel=cheb_tol_rel,
            tol_rel_dphi=cheb_tol_rel_dphi,
            big_prec=cheb_big_prec,
            maxseg=cheb_maxseg,
            minwidth=cheb_minwidth)
        dt = (time_ns() - t0) * 1e-9
        @info "RadialGGrid: Chebyshev cache built (segments=$(length(cheb.edges)-1)) in $(round(dt; digits=2)) s"
        f["$grp/cheb/edges"] = cheb.edges
        f["$grp/cheb/deg"] = cheb.deg
        f["$grp/cheb/u_max"] = cheb.u_max
        f["$grp/cheb/coeffs"] = cheb.coeffs
        f["$grp/cheb/dcoeffs"] = cheb.dcoeffs
        f["$grp/cheb/dcoeffs_kind"] = "dpdu"
        f["$grp/cheb/version"] = _BOUNDARY_CHEB_VERSION
        f["$grp/cheb/maxabs_phi"] = cheb.maxabs_phi
        f["$grp/cheb/maxabs_dphi"] = cheb.maxabs_dphi

        for jneg in jneg_list
            data = build_and_trim_for(jneg;
                                      H=H,
                                      stencil=stencil,
                                      jpos_build=jpos_build,
                                      jpos_save=jpos_save,
                                      grid_h=grid_h_use,
                                      grid_sigma=grid_sigma,
                                      grid_s0=grid_s0,
                                      rgridmax=rgridmax)
            grp = "jneg_$(jneg)"
            f["$grp/js"]      = data.js_trim
            f["$grp/C"]       = data.Ctrim
            f["$grp/centers"] = data.xctrim
            f["$grp/params"]  = string(data.params)
        end
    end
    @info "RadialGGrid: boundary cache ready."
end

"""
    make_erf_mapped_grid(; h, xmax, mapping, sigma=3.0, s0=6.5) -> (x, w)

Build a semi-infinite quadrature grid on `r∈[0, xmax]` by:
1) sampling midpoints `s = k + 0.5` with unit step in `s`,
2) applying the erf "switch-on" to build an auxiliary coordinate `u(s)`,
3) mapping to the physical radius `x = xofu(u, mapping)` using a
   `CoordinateMapping` object.

The returned weights `w` incorporate both Jacobians (`du/ds` and `dx/du`).

Parameter meanings / tips
- `h`: base step size in `u` once the switch is fully on (smaller → more points).
  In most drivers, this is the `hgrid` knob and is refined automatically in
  `build_radial_and_one_electron`.
- `xmax`: physical cutoff for the quadrature.
- `mapping`: a `CoordinateMapping` used to concentrate points near the core and
  become coarser at large `r`.
- `sigma`, `s0`: shape/offset of the switch-on `g(s)=0.5*erfc(s0 - s/sigma)`.
  `s0≈6.5` pushes the hard endpoint far enough that midpoint rules behave well.
"""
function make_erf_mapped_grid(; h::Real,
                               xmax::Real,
                               mapping,
                               sigma::Real=3.0,
                               s0::Real=6.5)

    (h > 0)    || error("h must be > 0")
    (xmax > 0) || error("xmax must be > 0")

    x = Float64[]
    w = Float64[]

    k = 0
    rtpi = sqrt(pi)

    while true
        s  = k + 0.5
        z  = (s / sigma) - s0

        g  = 0.5 * SpecialFunctions.erfc(-z)
        gp = exp(-z*z) / (rtpi * sigma)

        u  = h * s * g
        xu = xofu(u, mapping)
        rho_x = dudx(xu, mapping)
        wx = (h * (g + s*gp)) / rho_x

        push!(x, xu)
        push!(w, wx)

        xu >= xmax && break
        k += 1
    end

    return x, w
end

# Expand a seed-gausslet coefficient map (rows indexed by seed indices `js`,
# columns are basis functions) into coefficients in the underlying Gaussian
# basis used by Gausslets.jl.  Returned `kmin` is the minimum Gaussian index
# (in units of `g.sp`) for row 1 of `G`.
function gausslet_to_gaussian_coeffs(C::AbstractMatrix{Float64}, js::Vector{Int}, g)
    cfg = g.cfg
    shift = round(Int, 1 / g.sp)  # gausslet spacing (1) in Gaussian-grid units
    imin = Gausslets.mini(cfg)
    imax = Gausslets.maxi(cfg)

    jmin = minimum(js)
    jmax = maximum(js)
    kmin = shift * jmin + imin
    kmax = shift * jmax + imax

    N = size(C, 2)
    G = zeros(Float64, kmax - kmin + 1, N)
    cfg_range = Gausslets.arange(cfg)

    @inbounds for seed_idx in 1:length(js)
        j = js[seed_idx]
        base = shift * j - kmin
        for i in cfg_range
            rowG = base + i + 1
            ci = cfg[i]
            for a in 1:N
                G[rowG, a] += ci * C[seed_idx, a]
            end
        end
    end

    return G, kmin
end

"""
    build_radial_and_one_electron(; Z, corespacing, wi, s, Rmax, hgrid, sigma, kwargs...) -> RadialData

Main entry point: build the radial boundary-gausslet basis on `r≥0`, evaluate it
on an erf-mapped quadrature grid, and assemble one-electron matrices.

Returns a `RadialData` containing:
- quadrature grid `r`, weights `w`
- basis values `chiv[a](r_i)` and derivatives `chivp[a](r_i)`
- one-electron matrices `S, T, V, Vcentr` on the boundary basis
- IDA weights `wts[a] = ∫ χ_a(r) dr`
- basis centers and an `evalcache` used by the Chebyshev Vee builder

## Coordinate mapping / grid parameters
- `corespacing`: target physical spacing near the core. Smaller improves core
  resolution (and diagonal-approx behavior near the origin) but increases basis
  size and can worsen conditioning.
- `s`: shape parameter in the coordinate mapping (used together with
  `corespacing` as `c = corespacing/s` in `CoordinateMapping.Invsqrt`).
  Typical values in this repo are ~`0.1–0.2`; tune with hydrogenic checks.
- `wi`: outer-length knob for the mapping (used as `iw = 1/wi` in a constant
  mapping term).  Typical value is `wi≈10`.
- `Rmax`: physical radius used to *trim the basis* (keeps functions with centers
  `< Rmax`).  Increase to capture diffuse tails; decrease for speed.
- `grid_pad`: extra padding beyond `Rmax` for the *quadrature grid* so tails are
  integrated accurately. If you reduce it, check `||S-I||` and energies.

## Quadrature accuracy knobs
- `hgrid`: base step size used by `make_erf_mapped_grid` (smaller → more points).
- `sigma`, `S0`: parameters for the erf switch-on in the grid construction.
  `S0≈6.5` is chosen so the effective left endpoint is pushed far enough out
  that midpoint-style integration behaves well.
- `refine_hgrid=true`: by default, the function recursively halves `hgrid` until
  `||S-I||_Inf ≤ hgrid_tol` (or `hgrid_maxiter` is reached). This makes grid
  integration failures much harder to miss.
  - `hgrid_tol`: typical `1e-6`–`1e-8` for Float64 work; tighter costs more.
  - `hgrid_maxiter`: number of halvings allowed.

## Boundary cache / basis parameters
- `BGfilename`: JLD2 cache file to read/write boundary basis data.
- `boundary_method`:
  - `:odd_even` (default): symmetric `j=-L:L` seed window with odd/even build +
    optional even-complement modes (`neven_add`).
  - `:odd_even_tail`: symmetric `j=-L:L` seed window, but adds only a small even
    “tail window” `E_k` for `k=0:K` (controlled by `even_tail_K`), then removes
    the delta-at-origin direction and COMX-localizes again.
  - `:add_tails` (legacy): older asymmetric `j=-jneg:jpos_build` construction.
- `jpos_save`: `L` for the symmetric cache window size (used for `:odd_even` and `:odd_even_tail`).
- `neven_add`: number of even-complement eigenmodes to add (only used for `:odd_even`).
- `even_tail_K`: maximum even index `K` to keep (only used for `:odd_even_tail`).
- `jneg`: left-tail size for the legacy `:add_tails` cache (ignored for `:odd_even*`).
- `inject_count`: number of injected gpoly primitives added during cache build (0, 1, or 2).
- `inject_search`: search alpha(s) during cache build; if false uses `inject_alpha1/2`.
- `use_cheb_eval`: use cached piecewise Chebyshev χ(u) evaluation when available
  (recommended true).
- `use_gaussian_eval`: also build/store an underlying Gaussian expansion of χ(u)
  as a fallback for u beyond the Chebyshev range (recommended true).

Notes
- The basis is intended to be orthonormal in the continuum; any non-unitness of
  `S` is a proxy for numerical (quadrature/evaluation) error.
"""
function build_radial_and_one_electron(; Z::Real, corespacing::Real, wi::Real, s::Real, Rmax::Real,
    hgrid::Real, sigma::Real, jneg::Int=last(_BOUNDARY_JNEG_LIST),
    boundary_method::Symbol=_BOUNDARY_DEFAULT_METHOD,
    jpos_save::Int=_BOUNDARY_JPOS_SAVE,
    neven_add::Int=_BOUNDARY_ODD_EVEN_NEVEN_ADD,
    even_tail_K::Int=_BOUNDARY_ODD_EVEN_TAIL_K,
    inject_count::Int=_BOUNDARY_INJECT_COUNT,
    inject_mode::Symbol=_BOUNDARY_INJECT_MODE,
    inject_search::Bool=_BOUNDARY_INJECT_SEARCH,
    inject_alpha1::Real=_BOUNDARY_INJECT_ALPHA1,
    inject_alpha2::Real=_BOUNDARY_INJECT_ALPHA2,
    inject_alpha1_range::Tuple{Real,Real}=_BOUNDARY_INJECT_ALPHA1_RANGE,
    inject_alpha2_range::Tuple{Real,Real}=_BOUNDARY_INJECT_ALPHA2_RANGE,
    inject_alpha1_grid_n::Int=_BOUNDARY_INJECT_ALPHA1_GRID_N,
    inject_alpha2_grid_n::Int=_BOUNDARY_INJECT_ALPHA2_GRID_N,
    inject_alpha_refine_step::Real=_BOUNDARY_INJECT_ALPHA_REFINE_STEP,
    inject_alpha_refine_tol::Real=_BOUNDARY_INJECT_ALPHA_REFINE_TOL,
    inject_alpha_refine_maxiter::Int=_BOUNDARY_INJECT_ALPHA_REFINE_MAXITER,
    inject_multi_start::Bool=_BOUNDARY_INJECT_MULTI_START,
    inject_multi_start_n::Int=_BOUNDARY_INJECT_MULTI_START_N,
    inject_nm_maxiter::Int=_BOUNDARY_INJECT_NM_MAXITER,
    S0::Real=6.5, grid_pad::Real=80.0,
    BGfilename::AbstractString="BoundaryGausslets.jld2",
    use_cheb_eval::Bool=true,
    use_gaussian_eval::Bool=true,
    refine_hgrid::Bool=true,
    hgrid_tol::Real=1e-6,
    hgrid_maxiter::Int=4)
    if refine_hgrid
        @info "RadialGGrid: building radial data (boundary_method=$(boundary_method), jneg=$(jneg), neven_add=$(neven_add), even_tail_K=$(even_tail_K))"
        best_data = nothing
        best_dev = Inf
        h_try = Float64(hgrid)
        for iter in 1:hgrid_maxiter
            data = build_radial_and_one_electron(; Z,
                corespacing, wi, s, Rmax,
                hgrid=h_try, sigma, jneg,
                boundary_method, jpos_save, neven_add, even_tail_K,
                inject_count, inject_mode, inject_search,
                inject_alpha1, inject_alpha2,
                inject_alpha1_range, inject_alpha2_range,
                inject_alpha1_grid_n, inject_alpha2_grid_n,
                inject_alpha_refine_step, inject_alpha_refine_tol, inject_alpha_refine_maxiter,
                inject_multi_start, inject_multi_start_n, inject_nm_maxiter,
                S0, grid_pad,
                BGfilename,
                use_cheb_eval,
                use_gaussian_eval,
                refine_hgrid=false,
                hgrid_tol,
                hgrid_maxiter)
            dev = norm(data.S - I, Inf)
            @info "RadialGGrid: hgrid refine iter=$(iter) hgrid=$(h_try) ||S-I||_Inf=$(dev)"
            if dev < best_dev
                best_dev = dev
                best_data = data
            end
            dev <= hgrid_tol && return data
            h_try /= 2
        end
        @warn "RadialGGrid: hgrid refinement hit maxiter=$(hgrid_maxiter); best devS=$(best_dev)"
        return best_data::RadialData
    end

    iw = 1 / wi
    c  = corespacing / s
    mapping = CoordinateMapping.CombinedMappingFn(
        CoordinateMapping.Constfun(iw, 0.0),
        CoordinateMapping.Invsqrt[CoordinateMapping.Invsqrt(0.0, c, s)]
    )

    rmax = Rmax + grid_pad
    r, w = make_erf_mapped_grid(; h=hgrid, sigma, s0=S0, xmax=rmax, mapping)

    g = Gausslets.gausslet(Gausslets.cf1092)
    u_of_r(x) = CoordinateMapping.uofx(x, mapping)
    rho(r)         = CoordinateMapping.dudx(r, mapping)
    rho_prime(r)   = CoordinateMapping.du2dx2(r, mapping)
    jacobian(r)    = 1.0 / rho(r)

		    ensure_boundary_cache!(BGfilename; method=boundary_method, jneg, jpos_save, neven_add, even_tail_K,
                inject_count, inject_mode, inject_search,
                inject_alpha1, inject_alpha2,
                inject_alpha1_range, inject_alpha2_range,
                inject_alpha1_grid_n, inject_alpha2_grid_n,
                inject_alpha_refine_step, inject_alpha_refine_tol, inject_alpha_refine_maxiter,
                inject_multi_start, inject_multi_start_n, inject_nm_maxiter)

	        C, Cinj, inj_alphas, js, centersU, cheb = jldopen(BGfilename, "r") do f
	            grp_name = boundary_cache_group(boundary_method; jneg, jpos_save, neven_add, even_tail_K)
	            grp = f[grp_name]
	            C0 = Matrix(grp["C"])
	            Cinj0 = haskey(grp, "C_inj") ? Matrix(grp["C_inj"]) : zeros(Float64, 0, size(C0, 2))
	            inj_alphas0 = haskey(grp, "inject/alphas") ? Vector{Float64}(grp["inject/alphas"]) : Float64[]
	            js0 = Vector{Int}(grp["js"])
	            centers0 = Vector{Float64}(grp["centers"])
	            cheb0 = nothing
	            if use_cheb_eval && haskey(f, "$grp_name/cheb/edges")
	                edges = Vector{Float64}(f["$grp_name/cheb/edges"])
	                coeffs = Array{Float64,3}(f["$grp_name/cheb/coeffs"])
	                dcoeffs = Array{Float64,3}(f["$grp_name/cheb/dcoeffs"])
	                u_max = Float64(f["$grp_name/cheb/u_max"])
                kind = haskey(f, "$grp_name/cheb/dcoeffs_kind") ?
                    String(f["$grp_name/cheb/dcoeffs_kind"]) : "dpdx"
                cheb0 = (; edges, coeffs, dcoeffs, u_max, dcoeffs_kind=kind)
            end
	            (C0, Cinj0, inj_alphas0, js0, centers0, cheb0)
	        end
        nboundary_cache = size(C, 2)
	    realmax = round(Int, u_of_r(Rmax))
	    while js[end] < realmax
	        push!(js, js[end] + 1)
	        push!(centersU, float(js[end]))
	    end
	    centers = [CoordinateMapping.xofu(cu, mapping) for cu in centersU]

	    nseed = maximum(js) - minimum(js) + 1
	    if size(C,1) < nseed
	        old_rows = size(C, 1)
	        old_cols = size(C, 2)
	        nm = nseed - old_rows
	        Cbig = zeros(Float64, nseed, old_cols + nm)
	        Cbig[1:old_rows, 1:old_cols] .= C
	        Cbig[end-nm+1:end, old_cols+1:end] .= Matrix{Float64}(I, nm, nm)
	        C = Cbig
            if size(Cinj, 1) > 0
                Cinj = hcat(Cinj, zeros(Float64, size(Cinj, 1), nm))
            end
	    else
	        keeps   = findall(centers .< Rmax)
	        C       = C[:, keeps]
	        centers = centers[keeps]
            if size(Cinj, 1) > 0
                Cinj = Cinj[:, keeps]
            end
            if cheb !== nothing
                cheb = (; cheb..., coeffs=cheb.coeffs[:, :, keeps], dcoeffs=cheb.dcoeffs[:, :, keeps])
                nboundary_cache = size(C, 2)
            end
	    end

    chiv  = Matrix{Float64}(undef, length(r), size(C, 2))
    chivp = similar(chiv)

    ninj = size(Cinj, 1)
    has_inj = (ninj > 0) && (length(inj_alphas) == ninj)
    if ninj > 0 && length(inj_alphas) != ninj
        @warn "RadialGGrid: inject alphas length mismatch (ninj=$(ninj), nalpha=$(length(inj_alphas))); ignoring injection."
    end
    inj_vals = has_inj ? zeros(Float64, ninj) : Float64[]
    inj_dvals = has_inj ? zeros(Float64, ninj) : Float64[]
    N = size(C, 2)

    function add_injection!(ir::Int, u::Float64, sqrt_rho::Float64, rho_32::Float64, add_fac::Float64)
        has_inj || return nothing
        _gpoly_eval!(inj_vals, inj_dvals, u, _BOUNDARY_H, inj_alphas)
        @inbounds for a in 1:N
            s0 = 0.0
            s1 = 0.0
            for k in 1:ninj
                ck = Cinj[k, a]
                s0 += ck * inj_vals[k]
                s1 += ck * inj_dvals[k]
            end
            chiv[ir, a] += sqrt_rho * s0
            chivp[ir, a] += rho_32 * s1 + add_fac * s0
        end
        return nothing
    end

    if use_gaussian_eval
        ncheb = (cheb === nothing) ? 0 : size(cheb.coeffs, 3)
        u_max_grid = u_of_r(r[end])
        need_gaussian = (cheb === nothing) || (ncheb < size(C, 2)) || (u_max_grid > cheb.u_max)

        mids = Float64[]
        inv_halfwidths = Float64[]
        edges = Float64[]
        coeffs = Array{Float64,3}(undef, 0, 0, 0)
        dcoeffs = Array{Float64,3}(undef, 0, 0, 0)
        if cheb !== nothing
            edges = cheb.edges
            nseg = length(edges) - 1
            mids = Vector{Float64}(undef, nseg)
            inv_halfwidths = Vector{Float64}(undef, nseg)
            @inbounds for seg in 1:nseg
                a = edges[seg]
                b = edges[seg + 1]
                mids[seg] = 0.5 * (a + b)
                inv_halfwidths[seg] = 2.0 / (b - a)
            end
            coeffs = cheb.coeffs
            dcoeffs = cheb.dcoeffs
        end
        dcoeffs_dpdu = (cheb !== nothing) && (cheb.dcoeffs_kind == "dpdu")

        if !need_gaussian && cheb !== nothing && ncheb == size(C, 2)
            @inbounds for ir in 1:length(r)
                rr = r[ir]
                u = u_of_r(rr)
                seg = (u == edges[end]) ? (length(edges) - 1) : searchsortedlast(edges, u)
                seg < 1 && (seg = 1)
                seg > (length(edges) - 1) && (seg = length(edges) - 1)
                x = (u - mids[seg]) * inv_halfwidths[seg]

                rho_val = rho(rr)
                rho_p_val = rho_prime(rr)
                sqrt_rho = sqrt(rho_val)
                rho_32 = rho_val * sqrt_rho
                add_fac = 0.5 * (rho_p_val / sqrt_rho)

                for a in 1:ncheb
                    s0 = cheb_eval(view(coeffs, :, seg, a), x)
                    dp = cheb_eval(view(dcoeffs, :, seg, a), x)
                    s1 = dcoeffs_dpdu ? dp : (inv_halfwidths[seg] * dp)
                    chiv[ir, a]  = sqrt_rho * s0
                    chivp[ir, a] = rho_32 * s1 + add_fac * s0
                end
                add_injection!(ir, u, sqrt_rho, rho_32, add_fac)
            end
        else
        G, kmin = gausslet_to_gaussian_coeffs(C, js, g)
        kmax = kmin + size(G, 1) - 1
        N = size(C, 2)
        sp = g.sp
        wg = g.w
        inv_w2sp = 1.0 / (wg * wg * sp)

        # Reusable buffers for Gaussian evaluations (window size matches Gausslets.fg/dfg).
        win = 25
        ev_buf = Vector{Float64}(undef, win)
        dev_buf = Vector{Float64}(undef, win)

        @inbounds for ir in 1:length(r)
            rr = r[ir]
            u = u_of_r(rr)
            t = u / sp
            ii = round(Int, t)
            klo = max(ii - 12, kmin)
            khi = min(ii + 12, kmax)
            nk = khi - klo + 1

            k = klo
            for idx in 1:nk
                dk = t - k
                ev = exp(-0.5 * (dk / wg)^2)
                ev_buf[idx] = ev
                dev_buf[idx] = ev * (-(dk) * inv_w2sp)
                k += 1
            end
            row_lo = klo - kmin + 1

            rho_val = rho(rr)
            rho_p_val = rho_prime(rr)
            sqrt_rho = sqrt(rho_val)
            rho_32 = rho_val * sqrt_rho
            add_fac = 0.5 * (rho_p_val / sqrt_rho)

            use_cheb_here = (cheb !== nothing) && (u <= cheb.u_max)
            if use_cheb_here
                seg = (u == edges[end]) ? (length(edges) - 1) : searchsortedlast(edges, u)
                seg < 1 && (seg = 1)
                seg > (length(edges) - 1) && (seg = length(edges) - 1)
                x = (u - mids[seg]) * inv_halfwidths[seg]
                for a in 1:ncheb
                    s0 = cheb_eval(view(coeffs, :, seg, a), x)
                    dp = cheb_eval(view(dcoeffs, :, seg, a), x)
                    s1 = dcoeffs_dpdu ? dp : (inv_halfwidths[seg] * dp)
                    chiv[ir, a]  = sqrt_rho * s0
                    chivp[ir, a] = rho_32 * s1 + add_fac * s0
                end
            end

            a_start = use_cheb_here ? (ncheb + 1) : 1
            for a in a_start:N
                s0 = 0.0
                s1 = 0.0
                row = row_lo
                for idx in 1:nk
                    ga = G[row, a]
                    s0 += ga * ev_buf[idx]
                    s1 += ga * dev_buf[idx]
                    row += 1
                end
                chiv[ir, a]  = sqrt_rho * s0
                chivp[ir, a] = rho_32 * s1 + add_fac * s0
            end
            add_injection!(ir, u, sqrt_rho, rho_32, add_fac)
        end
        end
    else
        phi(j, x)  = g(x - j)
        dphi(j, x) = Gausslets.dfg(g, x - j)

        phiv = [phi(js[k], u_of_r(rr)) / sqrt(jacobian(rr)) for rr in r, k=1:size(C, 1)]
        chiv .= phiv * C

        function dchi_col(k, rr)
            u = u_of_r(rr)
            rho_val = rho(rr)
            rho_p_val = rho_prime(rr)
            phi_val = phi(js[k], u)
            phi_deriv = dphi(js[k], u)
            (rho_val^1.5) * phi_deriv + 0.5 * (rho_p_val / sqrt(rho_val)) * phi_val
        end

        phivp = [dchi_col(k, rr) for rr in r, k=1:size(C, 1)]
        chivp .= phivp * C

        if has_inj
            @inbounds for ir in 1:length(r)
                rr = r[ir]
                u = u_of_r(rr)
                rho_val = rho(rr)
                rho_p_val = rho_prime(rr)
                sqrt_rho = sqrt(rho_val)
                rho_32 = rho_val * sqrt_rho
                add_fac = 0.5 * (rho_p_val / sqrt_rho)
                add_injection!(ir, u, sqrt_rho, rho_32, add_fac)
            end
        end
    end

    S      = chiv' * (w .* chiv)
    T      = 0.5 * (chivp' * (w .* chivp))
    invr   = 1.0 ./ r
    V      = -Z * (chiv' * ((invr     .* w) .* chiv))
    Vcentr =        (chiv' * ((invr.^2 .* w) .* chiv))
    wts    = vec(chiv' * w)

    # Store enough information to re-evaluate χ(u) on alternate quadrature grids
    # (used by the Chebyshev Vee builder).  Kept as an `Any` to avoid making
    # RadialData parametric on mapping types.
    cheb_eval_data = nothing
    if cheb !== nothing
        edgesC = cheb.edges
        nsegC = length(edgesC) - 1
        midsC = Vector{Float64}(undef, nsegC)
        inv_halfwidthsC = Vector{Float64}(undef, nsegC)
        @inbounds for seg in 1:nsegC
            a = edgesC[seg]
            b = edgesC[seg + 1]
            midsC[seg] = 0.5 * (a + b)
            inv_halfwidthsC[seg] = 2.0 / (b - a)
        end
        cheb_eval_data = (; edges=edgesC,
                            mids=midsC,
                            inv_halfwidths=inv_halfwidthsC,
                            coeffs=cheb.coeffs,
                            u_max=cheb.u_max)
    end

    G_eval, kmin_eval = gausslet_to_gaussian_coeffs(C, js, g)
    gaussian_eval_data = (; G=G_eval,
                            kmin=kmin_eval,
                            kmax=kmin_eval + size(G_eval, 1) - 1,
                            sp=g.sp,
                            wg=g.w)

    inj_eval_data = has_inj ? (; alphas=inj_alphas, Cinj, H=_BOUNDARY_H) : nothing
    evalcache = (; mapping,
                  cheb=cheb_eval_data,
                  gaussian=gaussian_eval_data,
                  inj=inj_eval_data,
                  N=size(C, 2))

    return RadialData(
        collect(r),
        collect(w),
        Matrix{Float64}(chiv),
        Matrix{Float64}(chivp),
        collect(invr),
        Matrix{Float64}(S),
        Matrix{Float64}(T),
        Matrix{Float64}(V),
        Matrix{Float64}(Vcentr),
        collect(wts),
        size(C, 2),
        length(r),
        collect(centers),
        evalcache
    )
end

# ------------------------------------------------------------------------------
# Chebyshev/Clenshaw–Curtis Vee builder (high-accuracy alternative to trapezoid)
# ------------------------------------------------------------------------------

struct VeeChebSegment
    half::Float64             # (b-a)/2 in u-space
    drdu::Vector{Float64}     # dr/du at quadrature nodes
    wdr::Vector{Float64}      # dr weights: (half*w_cc).*drdu
    log_r::Vector{Float64}    # log(max(r,eps))
    chi::Matrix{Float64}      # χ(u) values (nodes × N)
end

# _xofu_newton: invert u(x) to x(u) robustly (bracket + Newton).
function _xofu_newton(u::Float64, mapping; rho0::Float64, tol_u::Float64=1e-14, maxiter::Int=30)
    u >= 0.0 || error("RadialGGrid: expected u>=0 for radial mapping, got u=$(u)")
    u == 0.0 && return 0.0

    lo = 0.0
    hi = max(u / rho0, 1e-12)
    uhi = CoordinateMapping.uofx(hi, mapping)
    grow = 0
    while uhi < u && grow < 80
        hi *= 2.0
        uhi = CoordinateMapping.uofx(hi, mapping)
        grow += 1
    end
    uhi < u && error("RadialGGrid: failed to bracket x(u): u=$(u) uofx(hi)=$(uhi) after grow=$(grow)")

    x = min(max(u / rho0, lo), hi)
    tol_abs = tol_u * max(1.0, u)
    for _ in 1:maxiter
        fx = CoordinateMapping.uofx(x, mapping) - u
        abs(fx) <= tol_abs && return x
        if fx > 0
            hi = x
        else
            lo = x
        end
        dfx = CoordinateMapping.dudx(x, mapping)
        x_new = x - fx / dfx
        if !(lo < x_new < hi) || !isfinite(x_new)
            x_new = 0.5 * (lo + hi)
        end
        x = x_new
    end
    return x
end

# _eval_phi_gaussian!: evaluate χ(u) in the underlying Gaussian basis (fast local window).
function _eval_phi_gaussian!(phi_out::Vector{Float64}, u::Float64, gaussian, ev_buf::Vector{Float64})
    G = gaussian.G
    kmin = gaussian.kmin
    kmax = gaussian.kmax
    sp = gaussian.sp
    wg = gaussian.wg
    N = size(G, 2)
    length(phi_out) == N || resize!(phi_out, N)

    t = u / sp
    ii = round(Int, t)
    klo = max(ii - 12, kmin)
    khi = min(ii + 12, kmax)
    nk = khi - klo + 1
    row_lo = klo - kmin + 1

    k = klo
    @inbounds for idx in 1:nk
        dk = t - k
        ev_buf[idx] = exp(-0.5 * (dk / wg)^2)
        k += 1
    end

    @inbounds for a in 1:N
        s0 = 0.0
        row = row_lo
        for idx in 1:nk
            s0 += G[row, a] * ev_buf[idx]
            row += 1
        end
        phi_out[a] = s0
    end
    return nothing
end

# _build_vee_u_edges: choose u-segment edges for the Chebyshev/CC Vee builder (geometric near 0 + cache-aligned edges + uniform tail).
function _build_vee_u_edges(evalcache, u_min::Float64, u_end::Float64;
        tail_du::Float64,
        geom_ratio::Float64=_VEE_CHEB_GEOM_RATIO,
        geom_u_max::Float64=_VEE_CHEB_GEOM_UMAX)
    (u_end > u_min) || error("_build_vee_u_edges: need u_end>u_min, got $(u_end) <= $(u_min)")
    tail_du > 0 || error("_build_vee_u_edges: tail_du must be > 0, got $(tail_du)")
    geom_ratio > 1 || error("_build_vee_u_edges: geom_ratio must be > 1, got $(geom_ratio)")
    geom_u_max > 0 || error("_build_vee_u_edges: geom_u_max must be > 0, got $(geom_u_max)")

    edges = Float64[u_min]

    # Near u=0, r(u) is small and r^L can vary by huge dynamic range for large L.
    # A geometric spacing up to `geom_u_max` keeps prefix integrals stable.
    u_geom_end = min(u_end, geom_u_max)
    if u_min < u_geom_end
        u = u_min
        while u < u_geom_end
            u2 = min(u * geom_ratio, u_geom_end)
            u2 <= u && break
            push!(edges, u2)
            u = u2
        end
    end

    cheb = evalcache.cheb
    if cheb !== nothing
        for e in cheb.edges
            ee = Float64(e)
            (ee <= u_geom_end) && continue
            (ee >= u_end) && break
            push!(edges, ee)
        end
    end

    last = maximum(edges)
    if last < u_end
        while last + tail_du < u_end
            last += tail_du
            push!(edges, last)
        end
        push!(edges, u_end)
    elseif last > u_end
        edges[end] = u_end
    end

    # Ensure strictly increasing edges (guard against duplicates from mixed sources).
    sort!(edges)
    edges = unique!(edges)
    if edges[1] != u_min
        edges[1] = u_min
    end
    if edges[end] != u_end
        push!(edges, u_end)
    end
    all(diff(edges) .> 0) || error("_build_vee_u_edges: non-increasing edges encountered")
    return edges
end

# _refine_edges: one refinement pass splitting segments wider than minwidth.
function _refine_edges(edges::Vector{Float64}; minwidth::Float64, max_edges::Int=8192)
    minwidth > 0 || error("_refine_edges: minwidth must be > 0, got $(minwidth)")
    new_edges = Float64[edges[1]]
    for i in 1:(length(edges) - 1)
        a = edges[i]
        b = edges[i + 1]
        if (b - a) > minwidth
            push!(new_edges, 0.5 * (a + b))
        end
        push!(new_edges, b)
        if length(new_edges) > max_edges
            @warn "RadialGGrid: edge refinement hit max_edges=$(max_edges); truncating refinement."
            break
        end
    end
    new_edges = unique!(sort!(new_edges))
    all(diff(new_edges) .> 0) || error("_refine_edges: produced non-increasing edges")
    return new_edges
end

# _refine_edges_to_minwidth: repeat edge refinement up to maxiter.
function _refine_edges_to_minwidth(edges::Vector{Float64}; minwidth::Float64, maxiter::Int, max_edges::Int)
    (minwidth > 0) || error("_refine_edges_to_minwidth: minwidth must be > 0, got $(minwidth)")
    (maxiter >= 0) || error("_refine_edges_to_minwidth: maxiter must be >= 0, got $(maxiter)")
    eds = edges
    for _ in 1:maxiter
        eds2 = _refine_edges(eds; minwidth=minwidth, max_edges=max_edges)
        if length(eds2) == length(eds)
            return eds
        end
        eds = eds2
    end
    return eds
end

# _build_vee_segments: precompute per-segment χ(u) values and CC weights used by the Vee builder.
function _build_vee_segments(evalcache, edges::Vector{Float64}; deg::Int)
    mapping = evalcache.mapping
    cheb = evalcache.cheb
    gaussian = evalcache.gaussian
    N = Int(evalcache.N)
    inj = hasproperty(evalcache, :inj) ? evalcache.inj : nothing
    has_inj = (inj !== nothing) && !isempty(inj.alphas)
    ninj = has_inj ? length(inj.alphas) : 0
    if has_inj && size(inj.Cinj, 1) != ninj
        @warn "RadialGGrid: inject alphas length mismatch in evalcache (ninj=$(size(inj.Cinj,1)) nalpha=$(ninj)); ignoring injection."
        has_inj = false
        ninj = 0
    end
    inj_vals = has_inj ? zeros(Float64, ninj) : Float64[]
    inj_dvals = has_inj ? zeros(Float64, ninj) : Float64[]
    Cinj = has_inj ? inj.Cinj : zeros(Float64, 0, N)
    H_inj = has_inj ? inj.H : _BOUNDARY_H
    rho0 = CoordinateMapping.dudx(0.0, mapping)

    xnodes = reverse(cheb_lobatto_nodes(deg, Float64))  # ascending: -1 → +1
    wcc = clenshaw_curtis_weights(deg, Float64; order=:ascending)
    P = cheb_prefix_integral_matrix(deg, Float64; order=:ascending)

    nseg = length(edges) - 1
    segments = Vector{VeeChebSegment}(undef, nseg)

    wts = zeros(Float64, N)
    tmpN = zeros(Float64, N)
    phi_tmp = zeros(Float64, N)
    ev_buf = Vector{Float64}(undef, 25)
    log_r_min = Inf
    log_r_max = -Inf

    # Cached Cheb lookup helpers for boundary functions.
    cheb_edges = cheb === nothing ? Float64[] : cheb.edges
    cheb_mids = cheb === nothing ? Float64[] : cheb.mids
    cheb_invhalf = cheb === nothing ? Float64[] : cheb.inv_halfwidths
    cheb_coeffs = cheb === nothing ? Array{Float64,3}(undef, 0, 0, 0) : cheb.coeffs
    cheb_u_max = cheb === nothing ? -Inf : Float64(cheb.u_max)
    ncheb = cheb === nothing ? 0 : size(cheb_coeffs, 3)

    @inbounds for s in 1:nseg
        a = edges[s]
        b = edges[s + 1]
        mid = 0.5 * (a + b)
        half = 0.5 * (b - a)

        u_nodes = mid .+ half .* xnodes
        nn = length(u_nodes)
        drdu = Vector{Float64}(undef, nn)
        wdr = Vector{Float64}(undef, nn)
        log_r = Vector{Float64}(undef, nn)
        rho = Vector{Float64}(undef, nn)

        for j in 1:nn
            u = u_nodes[j]
            r = _xofu_newton(u, mapping; rho0)
            rho_val = CoordinateMapping.dudx(r, mapping)
            rho[j] = rho_val
            drdu[j] = 1.0 / rho_val
            wdr[j] = half * wcc[j] * drdu[j]
            lr = log(max(r, eps(Float64)))
            log_r[j] = lr
            if lr < log_r_min
                log_r_min = lr
            end
            if lr > log_r_max
                log_r_max = lr
            end
        end

        chi = Matrix{Float64}(undef, nn, N)
        for j in 1:nn
            u = u_nodes[j]
            sqrt_rho = sqrt(rho[j])

            _eval_phi_gaussian!(phi_tmp, u, gaussian, ev_buf)

            if ncheb > 0 && u <= cheb_u_max
                seg = (u == cheb_edges[end]) ? (length(cheb_edges) - 1) : searchsortedlast(cheb_edges, u)
                seg < 1 && (seg = 1)
                seg > (length(cheb_edges) - 1) && (seg = length(cheb_edges) - 1)
                x = (u - cheb_mids[seg]) * cheb_invhalf[seg]
                for aidx in 1:ncheb
                    phi_tmp[aidx] = cheb_eval(view(cheb_coeffs, :, seg, aidx), x)
                end
            end

            if has_inj
                _gpoly_eval!(inj_vals, inj_dvals, u, H_inj, inj.alphas)
                @inbounds for aidx in 1:N
                    s0 = 0.0
                    for k in 1:ninj
                        s0 += Cinj[k, aidx] * inj_vals[k]
                    end
                    phi_tmp[aidx] += s0
                end
            end

            for aidx in 1:N
                chi[j, aidx] = sqrt_rho * phi_tmp[aidx]
            end
        end

        # Accumulate wts = ∫ χ(r) dr on the same quadrature.
        mul!(tmpN, chi', wdr)
        wts .+= tmpN

        segments[s] = VeeChebSegment(half, drdu, wdr, log_r, chi)
    end

    return (; segments, P, wts, log_r_min, log_r_max)
end

# _recover_value: undo log-scaling/shifts safely when reconstructing very large/small Vee contributions.
function _recover_value(sum_scaled::Float64, shift_total::Float64, log_extra::Float64,
        log_floatmin::Float64, log_floatmax::Float64)
    sum_scaled == 0.0 && return 0.0
    logv = shift_total + log_extra + log(abs(sum_scaled))
    if logv < log_floatmin
        return 0.0
    elseif logv > log_floatmax
        return copysign(Inf, sum_scaled)
    else
        return sign(sum_scaled) * exp(logv)
    end
end

# _build_Veel_cheb_on_edges: core CC/prefix-integral kernel building Veel[0:LmaxC] on a fixed u-edge partition.
function _build_Veel_cheb_on_edges(evalcache, edges::Vector{Float64}, LmaxC::Int; deg::Int)
    seginfo = _build_vee_segments(evalcache, edges; deg)
    segments = seginfo.segments
    P = seginfo.P
    wts = seginfo.wts
    log_r_min = seginfo.log_r_min
    log_r_max = seginfo.log_r_max

    N = Int(evalcache.N)
    any(wts .<= 0.0) && error("RadialGGrid: cheb Vee build got nonpositive wts; min(wts)=$(minimum(wts))")
    inv_wts = 1.0 ./ wts
    log_inv_wts = log.(inv_wts)

    log_floatmax = log(floatmax(Float64))
    log_floatmin = log(floatmin(Float64))

    Veel = [zeros(Float64, N, N) for _ in 0:LmaxC]

    nn = deg + 1
    H = Matrix{Float64}(undef, nn, N)
    Iinc = similar(H)
    Iseg = similar(H)
    WI = similar(H)
    tmpA = Matrix{Float64}(undef, N, N)
    A = zeros(Float64, N, N)
    Iprev = zeros(Float64, N)
    rpow = Vector{Float64}(undef, nn)
    invr = Vector{Float64}(undef, nn)
    w_outer = Vector{Float64}(undef, nn)

    @inbounds for L in 0:LmaxC
        fill!(A, 0.0)
        fill!(Iprev, 0.0)

        shift_r = (L == 0) ? 0.0 : (L * log_r_max)
        shift_invr = (-(L + 1)) * log_r_min
        shift_sum = shift_r + shift_invr

        for seg in segments
            log_r = seg.log_r
            chi = seg.chi
            drdu = seg.drdu
            wdr = seg.wdr

            for j in 1:nn
                lr = log_r[j]
                rpow[j] = (L == 0) ? 1.0 : exp(L * lr - shift_r)
                invr[j] = exp(-(L + 1) * lr - shift_invr)
                w_outer[j] = wdr[j] * invr[j]
            end

            for j in 1:nn
                scale = drdu[j] * rpow[j]
                for aidx in 1:N
                    H[j, aidx] = chi[j, aidx] * scale
                end
            end

            mul!(Iinc, P, H)
            Iinc .*= seg.half
            @views fill!(Iinc[1, :], 0.0)

            for j in 1:nn
                for aidx in 1:N
                    Iseg[j, aidx] = Iinc[j, aidx] + Iprev[aidx]
                end
            end

            for j in 1:nn
                wj = w_outer[j]
                for aidx in 1:N
                    WI[j, aidx] = Iseg[j, aidx] * wj
                end
            end

            mul!(tmpA, chi', WI)
            A .+= tmpA

            @views Iprev .= Iseg[end, :]
        end

        Sscaled = A .+ A'
        VL = Veel[L + 1]
        for i in 1:N
            for j in i:N
                log_extra = log_inv_wts[i] + log_inv_wts[j]
                v = _recover_value(Sscaled[i, j], shift_sum, log_extra, log_floatmin, log_floatmax)
                VL[i, j] = v
                VL[j, i] = v
            end
        end
    end

    return Veel
end

"""
    build_Veel_L_cheb(data, LmaxC; kwargs...) -> Vector{Matrix{Float64}}

High-accuracy builder for the radial multipole tables `Veel[L]` used by Ylm/Gaunt
angular machinery.

This uses piecewise Clenshaw–Curtis quadrature on segments in the auxiliary
coordinate `u` (not the physical radius `r`).  Segment edges are chosen to:
- resolve the near-origin dynamic range for large multipoles (`L`) via a short
  geometric progression in `u`,
- reuse/align with cached boundary-basis Chebyshev edges when available, and
- use a uniform tail spacing `tail_du` out to the end of the runtime grid.

Key parameters / tips
- `deg` (default 64): CC polynomial degree per segment. Larger → more accuracy,
  higher cost.  64 is usually a good Float64 target.
- `tail_du` (default 0.2): tail segment width in `u`. Smaller increases segments.
- `geom_ratio`, `geom_u_max`: near-0 geometric refinement; only matters for large
  `LmaxC` and tight core resolution.
- `refine_minwidth` + `refine_maxiter`: optional further splitting of long
  segments; helps avoid under-resolving oscillatory/core features.
- `deg_maxiter`: optional “self-check” that doubles the degree and compares the
  result; default 0 (off) for speed.
- `tol_rel`: relative tolerance for the degree-refinement check when enabled.
"""
function build_Veel_L_cheb(data::RadialData, LmaxC::Int;
        deg::Int=_VEE_CHEB_DEG,
        deg_max::Int=_VEE_CHEB_DEG_MAX,
        tail_du::Real=_VEE_CHEB_TAIL_DU,
        geom_ratio::Real=_VEE_CHEB_GEOM_RATIO,
        geom_u_max::Real=_VEE_CHEB_GEOM_UMAX,
        tol_rel::Real=_VEE_CHEB_TOL_REL,
        refine_minwidth::Real=_VEE_CHEB_REFINE_MINWIDTH,
        refine_maxiter::Int=_VEE_CHEB_REFINE_MAXITER,
        deg_maxiter::Int=_VEE_CHEB_DEG_MAXITER,
        max_edges::Int=_VEE_CHEB_MAX_EDGES)

    evalcache = data.evalcache
    evalcache === nothing && error("RadialGGrid: build_Veel_L_cheb requires data.evalcache (got nothing)")
    mapping = evalcache.mapping
    u_min = CoordinateMapping.uofx(data.r[1], mapping)
    u_end = CoordinateMapping.uofx(data.r[end], mapping)

    edges0 = _build_vee_u_edges(evalcache, Float64(u_min), Float64(u_end);
        tail_du=Float64(tail_du),
        geom_ratio=Float64(geom_ratio),
        geom_u_max=Float64(geom_u_max))
    edges = refine_minwidth > 0 ?
        _refine_edges_to_minwidth(edges0; minwidth=Float64(refine_minwidth), maxiter=refine_maxiter, max_edges=max_edges) :
        edges0
    @info "RadialGGrid: Vee(cheb) build: LmaxC=$(LmaxC) deg=$(deg) segments=$(length(edges)-1) tail_du=$(tail_du) geom_ratio=$(geom_ratio) geom_u_max=$(geom_u_max) minwidth=$(refine_minwidth)"

    deg_try = Int(deg)
    deg_try >= 8 || error("RadialGGrid: Vee(cheb) deg must be >= 8, got $(deg_try)")
    deg_max = max(Int(deg_try), Int(deg_max))

    t0 = time_ns()
    Veel_lo = _build_Veel_cheb_on_edges(evalcache, edges, LmaxC; deg=deg_try)
    dt = (time_ns() - t0) * 1e-9
    @info "RadialGGrid: Vee(cheb) deg=$(deg_try) done in $(round(dt; digits=2)) s"

    deg_maxiter <= 0 && return Veel_lo

    best = Veel_lo
    for iter in 1:deg_maxiter
        deg_next = min(2 * deg_try, deg_max)
        deg_next == deg_try && break
        t1 = time_ns()
        Veel_hi = _build_Veel_cheb_on_edges(evalcache, edges, LmaxC; deg=deg_next)
        dt_iter = (time_ns() - t1) * 1e-9

        max_rel = 0.0
        worst_L = 0
        for L in 0:LmaxC
            Vh = Veel_hi[L + 1]
            Vl = Veel_lo[L + 1]
            denom = max(norm(Vh, Inf), eps(Float64))
            rel = norm(Vh - Vl, Inf) / denom
            if rel > max_rel
                max_rel = rel
                worst_L = L
            end
        end
        @info "RadialGGrid: Vee(cheb) deg refine $(deg_try)→$(deg_next) in $(round(dt_iter; digits=2)) s max_rel=$(max_rel) (worst L=$(worst_L))"
        best = Veel_hi
        (max_rel <= tol_rel) && return best
        deg_try = deg_next
        Veel_lo = Veel_hi
    end

    @warn "RadialGGrid: Vee(cheb) did not reach tol_rel=$(tol_rel) within deg_maxiter=$(deg_maxiter); returning best (deg=$(deg_try))."
    return best
end

# build_Veel_L: legacy trapezoid/prefix-integral builder for the Coulomb IDA tensor blocks.
function build_Veel_L(r, w, chiv, invr, wts, N::Int, Ngr::Int, LmaxC::Int)
    inv_wts = 1.0 ./ wts
    log_inv_wts = log.(inv_wts)
    log_r = log.(max.(r, eps(Float64)))
    log_floatmax = log(floatmax(Float64))
    log_floatmin = log(floatmin(Float64))

    scaled_pow(logvals::Vector{Float64}, exponent::Int) =
        exponent == 0 ? (ones(Float64, length(logvals)), 0.0) :
        let logs = exponent .* logvals,
            shift = begin
                m = maximum(logs)
                isfinite(m) ? m : 0.0
            end
            vals = exp.(logs .- shift)
            (vals, shift)
        end

    function recover_value(sum_scaled::Float64, shift_total::Float64,
                           log_extra::Float64)
        if sum_scaled == 0.0
            return 0.0
        end
        logv = shift_total + log_extra + log(abs(sum_scaled))
        if logv < log_floatmin
            return 0.0
        elseif logv > log_floatmax
            return copysign(Inf, sum_scaled)
        else
            return sign(sum_scaled) * exp(logv)
        end
    end

    Veel = [zeros(Float64, N, N) for _ in 0:LmaxC]

    chivint = zeros(Float64, Ngr, N)
    @inbounds for i in 1:N
        chivint[1, i] = 0.5 * chiv[1, i] * w[1]
        for ir in 2:Ngr
            chivint[ir, i] = chivint[ir - 1, i] +
                0.5 * (chiv[ir - 1, i] * w[ir - 1] + chiv[ir, i] * w[ir])
        end
    end
    V0 = Veel[1]
    @inbounds for i in 1:N, j in i:N
        s = 0.0
        for ir in 1:Ngr
            s += w[ir] * (chiv[ir, i] * chivint[ir, j] + chivint[ir, i] * chiv[ir, j]) * invr[ir]
        end
        v = s * inv_wts[i] * inv_wts[j]
        V0[i, j] = V0[j, i] = v
    end

    for L in 1:LmaxC
        rpow_scaled, shift_r = scaled_pow(log_r, L)
        chivintL = zeros(Float64, Ngr, N)
        @inbounds for i in 1:N
            chivintL[1, i] = 0.5 * chiv[1, i] * w[1] * rpow_scaled[1]
            for ir in 2:Ngr
                chivintL[ir, i] = chivintL[ir - 1, i] +
                    0.5 * (chiv[ir - 1, i] * w[ir - 1] * rpow_scaled[ir - 1] +
                           chiv[ir, i] * w[ir] * rpow_scaled[ir])
            end
        end
        invr_scaled, shift_invr = scaled_pow(-log_r, L + 1)
        shift_sum = shift_r + shift_invr
        VL = Veel[L + 1]
        @inbounds for i in 1:N, j in i:N
            s = 0.0
            for ir in 1:Ngr
                s += w[ir] * invr_scaled[ir] *
                     (chiv[ir, i] * chivintL[ir, j] + chivintL[ir, i] * chiv[ir, j])
            end
            log_extra = log_inv_wts[i] + log_inv_wts[j]
            v = recover_value(s, shift_sum, log_extra)
            VL[i, j] = VL[j, i] = v
        end
    end
    return Veel
end

"""
    build_Veel_tables(data, Lc; method=:cheb, kwargs...) -> Vector{Matrix{Float64}}

Dispatch wrapper for the Coulomb radial multipole tables.
- `method=:cheb` (default): `build_Veel_L_cheb` (higher accuracy, robust for large `L`).
- `method=:trap`: `build_Veel_L` (legacy trapezoid rule; faster but lower order).
"""
function build_Veel_tables(data::RadialData, Lc::Int; method::Symbol=_VEE_DEFAULT_METHOD, kwargs...)
    if method == :trap
        return build_Veel_L(data.r, data.w, data.chiv, data.invr, data.wts, data.N, data.Ngr, Lc)
    elseif method == :cheb
        return build_Veel_L_cheb(data, Lc; kwargs...)
    else
        error("build_Veel_tables: unknown method=$(method); use :trap or :cheb")
    end
end

# radial_states: quick diagnostic eigenpairs for ℓ=0 and ℓ=1 radial one-electron blocks.
function radial_states(data::RadialData)
    H0 = data.T .+ data.V
    H1 = H0 .+ data.Vcentr
    eig0 = eigen(Symmetric(H0))
    eig1 = eigen(Symmetric(H1))
    return (; evals0 = eig0.values, vecs0 = eig0.vectors,
              evals1 = eig1.values, vecs1 = eig1.vectors)
end

end # module
