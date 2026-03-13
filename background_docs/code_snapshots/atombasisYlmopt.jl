#!/usr/bin/env julia
# atombasisYlm.jl — write a HamV6 (V6) JLD2 file using COMPLEX spherical harmonics Y_{ℓm}
#
# Purpose
# - Build the radial half-line basis (RadialGGrid) and its one-electron matrices.
# - Build the radial Coulomb multipole tables Veel[L] (IDA, L=0..2*lmax).
# - Build the angular Coulomb kernels QL[L] in the COMPLEX Y_{ℓm} basis (DiagYlm).
# - Write a JLD2 file matching `HamV6.addendum.txt` (COO sparse blocks, explicit (ℓ,m) ordering).
#
# Basis layout
# - Slices are radial basis functions (nslices = Nr = data.N).
# - Within each slice: orbitals are ordered by
#     (i)  m in zigzag order: 0, -1, +1, -2, +2, ...
#     (ii) for fixed m, increasing ℓ: ℓ = |m|, |m|+1, ..., lmax
#
# CLI overrides
#   julia atombasisYlm.jl Z=10 s=0.2 corespacing=0.02 lmax=3
#
# Output file name (default)
#   HamV6.<Elem>.s=<s>_l=<lmax>.jld2
#
# Notes
# - The radial basis is intended orthonormal; treat ||S-I|| as numerical error (no generalized eigensolve).
# - Two-body blocks are stored only for m-conserving pair sectors; strict m symmetry is expected in complex Ylm.

push!(LOAD_PATH, @__DIR__)
push!(LOAD_PATH, joinpath(@__DIR__, "sandbox_run", "modules"))

using LinearAlgebra
using Printf
using Dates
using JLD2

include("RadialGGrid.jl")
using .RadialGGrid

include("GauntTables.jl")
using .GauntTables

include("DiagYlm.jl")
using .DiagYlm

# ------------------------------- Defaults -------------------------------------

Z           = 10.0
corespacing = NaN             # auto -> 0.2/Z
wi          = 10.0
s           = 0.2
Rmax        = 12.0
hgrid       = 0.002
sigma       = 2.0
jneg        = 8
grid_pad    = 80.0
S0          = 6.5
cache_path  = RadialGGrid.DEFAULT_BOUNDARY_CACHE

boundary_method = :odd_even_tail   # good default for recent RadialGGrid
even_tail_K     = 6
neven_add       = 8
inject_count    = 2

lmax        = 1
drop_tol    = 1e-14
threads     = false

outfile     = ""   # if empty, auto-name

# ------------------------------- CLI parsing ----------------------------------

function _parse_arg_value(val::AbstractString)
    v = strip(String(val))
    if isempty(v)
        return ""
    end
    # quoted strings
    if length(v) >= 2 && ((v[1] == '"' && v[end] == '"') || (v[1] == '\'' && v[end] == '\''))
        return v[2:end-1]
    end
    # symbols like :odd_even_tail
    if startswith(v, ":")
        return Symbol(v[2:end])
    end
    lower = lowercase(v)
    if lower == "true";  return true;  end
    if lower == "false"; return false; end
    if occursin(r"^[+-]?\d+$", v); return parse(Int, v); end
    if occursin(r"^[+-]?(\d+\.\d*|\.\d+|\d+)([eE][+-]?\d+)?$", v); return parse(Float64, v); end
    return v
end

assigned = Dict{Symbol,Any}()
for arg in ARGS
    if occursin('=', arg)
        name, value = split(arg, '=', limit=2)
        assigned[Symbol(lowercase(strip(name)))] = _parse_arg_value(strip(value))
    else
        assigned[Symbol(lowercase(strip(arg)))] = true
    end
end

if haskey(assigned, :z);            Z           = Float64(assigned[:z]); end
if haskey(assigned, :s);            s           = Float64(assigned[:s]); end
if isnan(corespacing)
    corespacing = s * 0.5 / Z
end
if haskey(assigned, :corespacing);  corespacing = Float64(assigned[:corespacing]); end
if haskey(assigned, :wi);           wi          = Float64(assigned[:wi]); end
if haskey(assigned, :rmax);         Rmax        = Float64(assigned[:rmax]); end
if haskey(assigned, :hgrid);        hgrid       = Float64(assigned[:hgrid]); end
if haskey(assigned, :sigma);        sigma       = Float64(assigned[:sigma]); end
if haskey(assigned, :jneg);         jneg        = Int(assigned[:jneg]); end
if haskey(assigned, :grid_pad);     grid_pad    = Float64(assigned[:grid_pad]); end
if haskey(assigned, :s0);           S0          = Float64(assigned[:s0]); end
if haskey(assigned, :cache_path);   cache_path  = String(assigned[:cache_path]); end

if haskey(assigned, :boundary_method); boundary_method = Symbol(assigned[:boundary_method]); end
if haskey(assigned, :even_tail_k);     even_tail_K     = Int(assigned[:even_tail_k]); end
if haskey(assigned, :neven_add);       neven_add       = Int(assigned[:neven_add]); end
if haskey(assigned, :inject_count);    inject_count    = Int(assigned[:inject_count]); end

if haskey(assigned, :lmax);        lmax     = Int(assigned[:lmax]); end
if haskey(assigned, :drop_tol);    drop_tol = Float64(assigned[:drop_tol]); end
if haskey(assigned, :threads);     threads  = Bool(assigned[:threads]); end
if haskey(assigned, :outfile);     outfile  = String(assigned[:outfile]); end

# ------------------------------- Helpers --------------------------------------

const _ELEM = [
    "", "H","He","Li","Be","B","C","N","O","F","Ne",
    "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
    "In","Sn","Sb","Te","I","Xe",
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
    "Tl","Pb","Bi","Po","At","Rn",
    "Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
    "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn",
    "Nh","Fl","Mc","Lv","Ts","Og"
]

element_symbol(Z::Real) = begin
    Zi = Int(round(Z))
    if Zi >= 1 && Zi <= length(_ELEM) - 1 && isapprox(Z, Zi; atol=0, rtol=0)
        return _ELEM[Zi+1]
    end
    return "Z$(Zi)"
end

fmt_g(x::Real) = replace(@sprintf("%.6g", x), "+" => "")  # filename-friendly, modest precision

function mzigzag_list(lmax::Int)
    ms = Int[]
    push!(ms, 0)
    for k in 1:lmax
        push!(ms, -k)
        push!(ms, +k)
    end
    return ms
end

function lm_order_mzigzag_then_l(lmax::Int)
    l_list = Int[]
    m_list = Int[]
    for m in mzigzag_list(lmax)
        for l in abs(m):lmax
            push!(l_list, l)
            push!(m_list, m)
        end
    end
    @assert length(l_list) == (lmax + 1)^2
    return l_list, m_list
end

function pair_groups_by_msum(m_list::Vector{Int})
    nang = length(m_list)
    groups = Dict{Int, Vector{NTuple{3,Int}}}()
    for a in 1:nang, b in 1:nang
        msum = m_list[a] + m_list[b]
        row  = (a - 1) * nang + b
        push!(get!(groups, msum, NTuple{3,Int}[]), (a, b, row))
    end
    msums = sort(collect(keys(groups)))
    return msums, groups
end

function coo_block_from_diag(diag_vals::Vector{Float64}; drop_tol::Float64)
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    for i in eachindex(diag_vals)
        v = diag_vals[i]
        abs(v) > drop_tol || continue
        push!(rows, i)
        push!(cols, i)
        push!(vals, v)
    end
    return (rows=rows, cols=cols, vals=vals)
end

@inline function swap_pair_index(idx::Int, dim::Int)
    p = (idx - 1) ÷ dim + 1
    q = (idx - 1) % dim + 1
    return (q - 1) * dim + p
end

# ------------------------------- Main -----------------------------------------

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
let
    outname = outfile
    elem = element_symbol(Z)
    if isempty(outname)
        outname = "HamV6.$(elem).s=$(fmt_g(s))_l=$(lmax).jld2"
    end

    @printf("atombasisYlm: Z=%.6g (%s)  corespacing=%.6g  s=%.6g  wi=%.6g  Rmax=%.6g  hgrid=%.6g  sigma=%.6g  jneg=%d\n",
            Z, elem, corespacing, s, wi, Rmax, hgrid, sigma, jneg)
    @printf("             boundary_method=%s  neven_add=%d  even_tail_K=%d  inject_count=%d\n",
            string(boundary_method), neven_add, even_tail_K, inject_count)
    @printf("             lmax=%d  LmaxC=%d  drop_tol=%.1e  threads=%s\n",
            lmax, 2*lmax, drop_tol, string(threads))
    @printf("             outfile=%s\n", outname)

    # (ℓ,m) ordering required by HamV6.addendum.txt
    l_list, m_list = lm_order_mzigzag_then_l(lmax)
    nang = length(l_list)
    @printf("Angular orbitals per slice: Nang=%d (= (lmax+1)^2)\n", nang)

    # Map our (ℓ,m) ordering to DiagYlm's packed ordering (ℓ panels, m runs fastest)
    lmtoj, _, Nj, _ = DiagYlm.lm_indexers(lmax)
    @assert Nj == nang
    perm = [lmtoj(l_list[i], m_list[i]) for i in 1:nang]  # our_index -> diagYlm_index

    # Precompute pair sectors by m_sum (to enforce / exploit strict m conservation)
    msums, pairs_by_msum = pair_groups_by_msum(m_list)
    @printf("m-sum sectors for pair indices: [%s]\n", join(string.(msums), ", "))

    # Radial payload
    t_rad = time()
    data = RadialGGrid.build_radial_and_one_electron(; Z=Z,
        corespacing=corespacing, wi=wi, s=s, Rmax=Rmax,
        hgrid=hgrid, sigma=sigma, jneg=jneg,
        boundary_method=boundary_method, even_tail_K=even_tail_K,
        neven_add=neven_add, inject_count=inject_count,
        S0=S0, grid_pad=grid_pad, BGfilename=cache_path)
    @printf("RadialGGrid: Nr=%d  grid_pts=%d  ||S-I||_Inf=%.3e  build_time=%.2fs\n",
            data.N, data.Ngr, norm(data.S - I, Inf), time() - t_rad)

    nslices = data.N
    dims = fill(nang, nslices)
    offs = Vector{Int}(undef, nslices + 1)
    offs[1] = 1
    for n in 1:nslices
        offs[n+1] = offs[n] + dims[n]
    end
    Norb = offs[end] - 1
    @printf("Total orbitals (radial slices × angular): Norb=%d\n", Norb)

    # Radial multipoles and angular kernels
    LmaxC = 2*lmax
    t_veel = time()
    Veel = RadialGGrid.build_Veel_tables(data, LmaxC)
    @printf("RadialGGrid: built Veel[0:%d] in %.2fs\n", LmaxC, time() - t_veel)

    t_q = time()
    QL, _, _, _, _, _, _ = DiagYlm.build_Qkernels(lmax; Lmax=LmaxC, basis=:complex, threads=threads)
    @printf("DiagYlm: built QL[0:%d] (basis=:complex) in %.2fs\n", LmaxC, time() - t_q)

    # Precompute Q blocks on each m-sum sector in our pair ordering.
    # This avoids repeated 4D indexing inside the expensive (radial-pair × pair-pair) loops.
    t_qsec = time()
    nsectors = length(msums)
    sector_rows = Vector{Vector{Int}}(undef, nsectors)
    sector_ad   = Vector{Vector{Int}}(undef, nsectors)
    sector_bd   = Vector{Vector{Int}}(undef, nsectors)
    for (si, msum) in enumerate(msums)
        plist = pairs_by_msum[msum]  # Vector{(α,β,row)} in our within-slice ordering
        n = length(plist)
        rowids = Vector{Int}(undef, n)
        ad = Vector{Int}(undef, n)
        bd = Vector{Int}(undef, n)
        @inbounds for t in 1:n
            α, β, row = plist[t]
            rowids[t] = row
            ad[t] = perm[α]   # to DiagYlm packed index
            bd[t] = perm[β]
        end
        sector_rows[si] = rowids
        sector_ad[si] = ad
        sector_bd[si] = bd
    end

    Qsec = [Vector{Matrix{Float64}}(undef, nsectors) for _ in 1:(LmaxC + 1)]
    for Lp in 1:(LmaxC + 1)
        Q = QL[Lp]
        for si in 1:nsectors
            ad = sector_ad[si]
            bd = sector_bd[si]
            n = length(ad)
            M = Matrix{Float64}(undef, n, n)
            @inbounds for j in 1:n
                adj = ad[j]
                bdj = bd[j]
                for i in 1:n
                    M[i, j] = Q[ad[i], adj, bd[i], bdj]
                end
            end
            Qsec[Lp][si] = M
        end
    end
    nnz_full_per_block = sum(length(sector_rows[si])^2 for si in 1:nsectors)
    @printf("Prepared sector Q blocks in %.2fs  (m-sectors=%d, nnz_full_per_block=%d)\n",
            time() - t_qsec, nsectors, nnz_full_per_block)

    # One-body COO blocks
    t_1e = time()
    H1blocks = [Vector{NamedTuple{(:rows,:cols,:vals),Tuple{Vector{Int},Vector{Int},Vector{Float64}}}}(undef, nslices) for _ in 1:nslices]
    Tblocks  = [Vector{NamedTuple{(:rows,:cols,:vals),Tuple{Vector{Int},Vector{Int},Vector{Float64}}}}(undef, nslices) for _ in 1:nslices]
    Vnucblocks = [Vector{NamedTuple{(:rows,:cols,:vals),Tuple{Vector{Int},Vector{Int},Vector{Float64}}}}(undef, nslices) for _ in 1:nslices]

    for a in 1:nslices, b in 1:nslices
        Tnm = data.T[a,b]
        Vnm = data.V[a,b]
        Vcnm = data.Vcentr[a,b]
        hdiag = Vector{Float64}(undef, nang)
        tdiag = fill(Tnm, nang)
        vdiag = fill(Vnm, nang)
        @inbounds for i in 1:nang
            l = l_list[i]
            hdiag[i] = Tnm + Vnm + 0.5 * l * (l + 1) * Vcnm
        end
        H1blocks[a][b] = coo_block_from_diag(hdiag; drop_tol=drop_tol)
        Tblocks[a][b]  = coo_block_from_diag(tdiag; drop_tol=drop_tol)
        Vnucblocks[a][b] = coo_block_from_diag(vdiag; drop_tol=drop_tol)
    end
    @printf("Built one-body COO blocks in %.2fs\n", time() - t_1e)

    # Two-body COO blocks
    t_2e = time()
    Vblocks = [Vector{NamedTuple{(:rows,:cols,:vals),Tuple{Vector{Int},Vector{Int},Vector{Float64}}}}(undef, nslices) for _ in 1:nslices]
    prefL = Vector{Float64}(undef, LmaxC + 1)
    for a in 1:nslices
        for b in a:nslices
            rows = Int[]
            cols = Int[]
            vals = Float64[]
            sizehint!(rows, nnz_full_per_block)
            sizehint!(cols, nnz_full_per_block)
            sizehint!(vals, nnz_full_per_block)

            # prefactors per L for this (a,b)
            maxpref = 0.0
            @inbounds for Lp in 1:(LmaxC + 1)
                v = Veel[Lp][a, b]
                prefL[Lp] = v
                av = abs(v)
                av > maxpref && (maxpref = av)
            end
            # short-circuit if everything is numerically zero
            if maxpref <= drop_tol
                Vblocks[a][b] = (rows=rows, cols=cols, vals=vals)
                if b != a
                    Vblocks[b][a] = (rows=Int[], cols=Int[], vals=Float64[])
                end
                continue
            end

            # Exploit (pq|rs)=(rs|pq): each (a,b) pair-block is symmetric, so compute only upper-tri entries
            # within each m-sum sector and mirror entries into COO.
            @inbounds for si in 1:nsectors
                rowids = sector_rows[si]
                n = length(rowids)
                for j in 1:n
                    colj = rowids[j]
                    for i in 1:j
                        v = 0.0
                        for Lp in 1:(LmaxC + 1)
                            p = prefL[Lp]
                            p == 0.0 && continue
                            v += p * Qsec[Lp][si][i, j]
                        end
                        abs(v) > drop_tol || continue
                        rowi = rowids[i]
                        push!(rows, rowi)
                        push!(cols, colj)
                        push!(vals, v)
                        if rowi != colj
                            push!(rows, colj)
                            push!(cols, rowi)
                            push!(vals, v)
                        end
                    end
                end
            end

            Vblocks[a][b] = (rows=rows, cols=cols, vals=vals)

            # Fill the opposite slice-pair block (b,a) without recomputation.
            # Using (pq|rs) = (qp|sr) under swap of electron labels 1↔2:
            #   row=(p,q) -> (q,p),  col=(r,s) -> (s,r)
            if b != a
                rows2 = similar(rows)
                cols2 = similar(cols)
                @inbounds for k in eachindex(rows)
                    rows2[k] = swap_pair_index(rows[k], nang)
                    cols2[k] = swap_pair_index(cols[k], nang)
                end
                # vals are identical but MUST NOT be aliased across blocks (downstream code may mutate COO buffers).
                Vblocks[b][a] = (rows=rows2, cols=cols2, vals=copy(vals))
            end
        end
    end

    nstored_total = sum(length(Vblocks[a][b].vals) for a in 1:nslices for b in 1:nslices)
    @printf("Built two-body COO blocks in %.2fs  (stored nnz total=%d)\n", time() - t_2e, nstored_total)

    # Basis metadata (per-slice vectors; keep explicit and redundant for simplicity downstream)
    m_by_slice = [copy(m_list) for _ in 1:nslices]
    l_by_slice = [copy(l_list) for _ in 1:nslices]
    labels_by_slice = [ [@sprintf("l=%d m=%d", l_list[i], m_list[i]) for i in 1:nang] for _ in 1:nslices ]

    created = Dates.format(Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS")
    @printf("Writing %s ...\n", outname)
    t_write = time()
    jldopen(outname, "w") do f
        # meta
        f["meta/version"] = "HamV6.v0.2.coo"
        f["meta/created"] = created
        f["meta/units"] = "Hartree"
        f["meta/eltype"] = "Float64"
        f["meta/twobody_convention"] = "physicists: (pq|rs), V = 1/2 sum (pq|rs) c†_p c†_q c_s c_r"
        f["meta/notes"] = "complex Ylm; slices are radial basis functions; within slice mzigzag_then_l"
        f["meta/Z"] = Z
        f["meta/s"] = s
        f["meta/lmax"] = lmax
        f["meta/drop_tol"] = drop_tol

        # layout
        f["layout/nslices"] = nslices
        f["layout/dims"] = dims
        f["layout/offs"] = offs
        f["layout/slice_coord"] = data.centers

        # ordering
        f["ordering/within_slice"] = "mzigzag_then_l"
        f["ordering/description"] = "slice-major; within slice m=0,-1,+1,-2,+2,...; for each m, ℓ increasing"

        # basis
        f["basis/m_by_slice"] = m_by_slice
        f["basis/l_by_slice"] = l_by_slice
        f["basis/m_flat"] = vcat(m_by_slice...)
        f["basis/l_flat"] = vcat(l_by_slice...)
        f["basis/labels_by_slice"] = labels_by_slice

        # radial params (for provenance)
        f["radial/corespacing"] = corespacing
        f["radial/wi"] = wi
        f["radial/s"] = s
        f["radial/Rmax"] = Rmax
        f["radial/hgrid"] = hgrid
        f["radial/sigma"] = sigma
        f["radial/jneg"] = jneg
        f["radial/boundary_method"] = string(boundary_method)
        f["radial/even_tail_K"] = even_tail_K
        f["radial/neven_add"] = neven_add
        f["radial/inject_count"] = inject_count
        f["radial/||S-I||_Inf"] = norm(data.S - I, Inf)

        # one-body
        f["onebody/stored"] = "coo"
        f["onebody/is_hermitian"] = true
        f["onebody/H1blocks"] = H1blocks
        f["onebody/Tblocks"] = Tblocks
        f["onebody/Vnucblocks"] = Vnucblocks

        # two-body
        f["twobody/stored"] = "coo_all"
        f["twobody/convention"] = "COO over pair indices: row=(p,q) with p in slice n, q in slice m; col=(r,s) with r in slice n, s in slice m; vals=(pq|rs)."
        f["twobody/symmetry"] = "real_coulomb in complex Ylm: exact m conservation; (pq|rs)=(qp|sr)=(rs|pq)"
        f["twobody/Vblocks"] = Vblocks
    end
    @printf("Wrote %s in %.2fs\n", outname, time() - t_write)
end
end
