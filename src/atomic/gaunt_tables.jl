"""
    GauntEntry{T}

One nonzero Gaunt coefficient entry in a sparse `(L, l1, l2)` block.

The stored indices are:
- `m1` for the left one-electron channel
- `M` for the multipole channel
- `m2` for the right one-electron channel
"""
struct GauntEntry{T <: Real}
    m1::Int
    M::Int
    m2::Int
    val::T
end

"""
    GauntTable{T}

Sparse block Gaunt table up to fixed `lmax` and `Lmax`.

For each `L`, the table stores an `(lmax + 1) × (lmax + 1)` matrix whose
`(l1 + 1, l2 + 1)` entry is a `Vector{GauntEntry{T}}`.
"""
struct GauntTable{T <: Real}
    lmax::Int
    Lmax::Int
    blocks::Vector{Matrix{Vector{GauntEntry{T}}}}
end

gaunt_lmax(table::GauntTable) = table.lmax
gaunt_Lmax(table::GauntTable) = table.Lmax

@inline function gaunt_legal_triple(l1::Int, l2::Int, L::Int)
    return (abs(l1 - l2) <= L) && (l1 + l2 >= L) && iseven(l1 + l2 + L)
end

@inline function gaunt_legal_ms(l1::Int, m1::Int, l2::Int, m2::Int, L::Int, M::Int)
    return abs(m1) <= l1 && abs(m2) <= l2 && abs(M) <= L
end

@inline function gaunt_block(table::GauntTable, L::Int, l1::Int, l2::Int)
    0 <= L <= table.Lmax || throw(BoundsError(table.blocks, L + 1))
    0 <= l1 <= table.lmax || throw(BoundsError(table.blocks[L + 1], l1 + 1))
    0 <= l2 <= table.lmax || throw(BoundsError(table.blocks[L + 1], l2 + 1))
    return table.blocks[L + 1][l1 + 1, l2 + 1]
end

gaunt_hasblock(table::GauntTable, L::Int, l1::Int, l2::Int) = !isempty(gaunt_block(table, L, l1, l2))
gaunt_nnz(table::GauntTable) = sum(gaunt_nnz(table, L) for L in 0:table.Lmax)
gaunt_nnz(table::GauntTable, L::Int) = sum(length, table.blocks[L + 1])
gaunt_nnz(table::GauntTable, L::Int, l1::Int, l2::Int) = length(gaunt_block(table, L, l1, l2))

function gaunt_each_block(table::GauntTable, L::Int)
    blocks = table.blocks[L + 1]
    return ((l1 - 1, l2 - 1, blocks[l1, l2])
        for l1 in axes(blocks, 1), l2 in axes(blocks, 2)
        if !isempty(blocks[l1, l2]))
end

function gaunt_each_nonzero(table::GauntTable, L::Int, l1::Int, l2::Int)
    return ((entry.m1, entry.M, entry.m2, entry.val) for entry in gaunt_block(table, L, l1, l2))
end

function gaunt_value(
    table::GauntTable{T},
    L::Int,
    l1::Int,
    m1::Int,
    l2::Int,
    m2::Int,
    M::Int,
) where {T}
    for entry in gaunt_block(table, L, l1, l2)
        if entry.m1 == m1 && entry.M == M && entry.m2 == m2
            return entry.val
        end
    end
    return zero(T)
end

function gaunt_allowed_triples(lmax::Int; Lmax::Int = 2 * lmax)
    triples = NTuple{3,Int}[]
    for L in 0:Lmax, l1 in 0:lmax, l2 in 0:lmax
        gaunt_legal_triple(l1, l2, L) && push!(triples, (l1, l2, L))
    end
    return triples
end

function _factorial_table(nmax::Int)
    return [BigFloat(factorial(big(n))) for n in 0:nmax]
end

@inline _fact(table::AbstractVector{BigFloat}, n::Int) = table[n + 1]

function _wigner3j(
    fact_table::AbstractVector{BigFloat},
    j1::Int,
    j2::Int,
    j3::Int,
    m1::Int,
    m2::Int,
    m3::Int,
)
    m1 + m2 + m3 == 0 || return 0.0
    abs(m1) <= j1 || return 0.0
    abs(m2) <= j2 || return 0.0
    abs(m3) <= j3 || return 0.0
    gaunt_legal_triple(j1, j2, j3) || return 0.0

    delta = sqrt(
        _fact(fact_table, j1 + j2 - j3) *
        _fact(fact_table, j1 - j2 + j3) *
        _fact(fact_table, -j1 + j2 + j3) /
        _fact(fact_table, j1 + j2 + j3 + 1),
    )

    prefactor = (isodd(j1 - j2 - m3) ? -one(BigFloat) : one(BigFloat)) * delta
    prefactor *= sqrt(
        _fact(fact_table, j1 + m1) *
        _fact(fact_table, j1 - m1) *
        _fact(fact_table, j2 + m2) *
        _fact(fact_table, j2 - m2) *
        _fact(fact_table, j3 + m3) *
        _fact(fact_table, j3 - m3),
    )

    tmin = max(0, j2 - j3 - m1, j1 - j3 + m2)
    tmax = min(j1 + j2 - j3, j1 - m1, j2 + m2)
    tmin <= tmax || return 0.0

    total = zero(BigFloat)
    for t in tmin:tmax
        denominator =
            _fact(fact_table, t) *
            _fact(fact_table, j1 + j2 - j3 - t) *
            _fact(fact_table, j1 - m1 - t) *
            _fact(fact_table, j2 + m2 - t) *
            _fact(fact_table, j3 - j2 + m1 + t) *
            _fact(fact_table, j3 - j1 - m2 + t)
        total += (isodd(t) ? -one(BigFloat) : one(BigFloat)) / denominator
    end

    return Float64(prefactor * total)
end

@inline function _gaunt_complex(
    fact_table::AbstractVector{BigFloat},
    l1::Int,
    m1::Int,
    l2::Int,
    m2::Int,
    L::Int,
    M::Int,
)
    gaunt_legal_triple(l1, l2, L) || return 0.0
    -m1 + M + m2 == 0 || return 0.0

    prefactor = (isodd(m1) ? -1.0 : 1.0) *
                sqrt(((2 * l1 + 1) * (2 * l2 + 1) * (2 * L + 1)) / (4 * pi))
    return prefactor *
           _wigner3j(fact_table, l1, L, l2, 0, 0, 0) *
           _wigner3j(fact_table, l1, L, l2, -m1, M, m2)
end

@inline function _to_real_pairs(l::Int, μ::Int)
    if μ == 0
        return ((0, 1.0, false),)
    elseif μ > 0
        m = μ
        coefficient = inv(sqrt(2.0))
        return (
            (m, coefficient, false),
            (-m, (iseven(m) ? 1.0 : -1.0) * coefficient, false),
        )
    else
        m = -μ
        coefficient = inv(sqrt(2.0))
        return (
            (m, coefficient, true),
            (-m, -(iseven(m) ? 1.0 : -1.0) * coefficient, true),
        )
    end
end

@inline function _gaunt_real(
    fact_table::AbstractVector{BigFloat},
    l1::Int,
    μ1::Int,
    l2::Int,
    μ2::Int,
    L::Int,
    Μ::Int,
)
    real_part = 0.0
    imag_part = 0.0

    for (m1, c1raw, imag1) in _to_real_pairs(l1, μ1)
        c1 = imag1 ? -c1raw : c1raw
        for (m2, c2, imag2) in _to_real_pairs(l2, μ2)
            for (mM, cM, imagM) in _to_real_pairs(L, Μ)
                value = _gaunt_complex(fact_table, l1, m1, l2, m2, L, mM)
                coefficient = c1 * c2 * cM
                imag_count = (imag1 ? 1 : 0) + (imag2 ? 1 : 0) + (imagM ? 1 : 0)

                if iseven(imag_count)
                    real_part += (imag_count % 4 == 0 ? coefficient : -coefficient) * value
                else
                    imag_part += (imag_count % 4 == 1 ? coefficient : -coefficient) * value
                end
            end
        end
    end

    abs(imag_part) <= 1.0e-12 * (abs(real_part) + 1.0) ||
        throw(ArgumentError("real Gaunt construction produced a non-negligible imaginary residue"))
    return real_part
end

function _empty_gaunt_blocks(::Type{T}, lmax::Int) where {T <: Real}
    return Matrix{Vector{GauntEntry{T}}}([GauntEntry{T}[] for _ in 1:(lmax + 1), __ in 1:(lmax + 1)])
end

function _build_gaunt_blocks_for_L(::Type{T}, lmax::Int, L::Int, atol::Real, basis::Symbol) where {T <: Real}
    blocks = _empty_gaunt_blocks(T, lmax)
    fact_table = _factorial_table(max(1, 4 * max(lmax, L) + 1))

    for l1 in 0:lmax, l2 in 0:lmax
        gaunt_legal_triple(l1, l2, L) || continue
        entries = GauntEntry{T}[]

        if basis === :complex
            for m1 in -l1:l1, m2 in -l2:l2
                M = m1 - m2
                abs(M) <= L || continue
                value = _gaunt_complex(fact_table, l1, m1, l2, m2, L, M)
                abs(value) > atol || continue
                push!(entries, GauntEntry{T}(m1, M, m2, T(value)))
            end
        elseif basis === :real
            for M in -L:L, m1 in -l1:l1, m2 in -l2:l2
                gaunt_legal_ms(l1, m1, l2, m2, L, M) || continue
                value = _gaunt_real(fact_table, l1, m1, l2, m2, L, M)
                abs(value) > atol || continue
                push!(entries, GauntEntry{T}(m1, M, m2, T(value)))
            end
        else
            throw(ArgumentError("build_gaunt_table basis must be :complex or :real"))
        end

        if !isempty(entries)
            sort!(entries; by = entry -> (entry.M, entry.m1, entry.m2))
            blocks[l1 + 1, l2 + 1] = entries
        end
    end

    return blocks
end

"""
    build_gaunt_table(lmax; Lmax = 2 * lmax, atol = 0.0, basis = :complex, T = Float64)

Build a sparse/block Gaunt table using the package's internal Wigner-3j
evaluation rather than an external angular-special-function dependency.
"""
function build_gaunt_table(
    lmax::Int;
    Lmax::Int = 2 * lmax,
    atol::Real = 0.0,
    basis::Symbol = :complex,
    T::Type{<:Real} = Float64,
)
    lmax >= 0 || throw(ArgumentError("build_gaunt_table requires lmax >= 0"))
    Lmax >= 0 || throw(ArgumentError("build_gaunt_table requires Lmax >= 0"))
    basis in (:complex, :real) || throw(ArgumentError("build_gaunt_table basis must be :complex or :real"))

    blocks = Vector{Matrix{Vector{GauntEntry{T}}}}(undef, Lmax + 1)
    for L in 0:Lmax
        blocks[L + 1] = _build_gaunt_blocks_for_L(T, lmax, L, atol, basis)
    end
    return GauntTable{T}(lmax, Lmax, blocks)
end

