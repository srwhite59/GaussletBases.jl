module DiagYlm

using LinearAlgebra
# Prefer an already-loaded `Main.GauntTables` (e.g. a local dev copy included by a driver)
# to avoid accidentally picking up a different GauntTables on LOAD_PATH.
if isdefined(Main, :GauntTables)
    import Main: GauntTables
else
    using GauntTables
end

export build_Qkernels, getJ, getK, multV, lm_indexers

"""
    lm_indexers(lmax; maxm=-1)

Return `(lmtoj, jtolm, Nj, ang_offset)` for Y_{ℓm} with 0≤ℓ≤lmax and
m=-mmax(ℓ):mmax(ℓ), where `mmax(ℓ) = (maxm < 0) ? ℓ : min(ℓ, maxm)`.
- `lmtoj(ℓ,m)` → packed j (m runs fastest per ℓ)
- `jtolm(j)`   → (ℓ,m)
- `Nj`         → total number of (ℓ,m)
- `ang_offset` → base index per ℓ panel
"""
function lm_indexers(lmax; maxm::Int=-1)
    ang_offset = Int[]
    mmax_list = Int[]
    base = 1
    for l in 0:lmax
        mmax_l = (maxm < 0) ? l : min(l, maxm)
        push!(ang_offset, base)
        push!(mmax_list, mmax_l)
        base += 2*mmax_l + 1
    end
    Nj = base - 1
    lmtoj = (l,m)->begin
        @assert 0 ≤ l ≤ lmax "l out of range"
        mmax_l = mmax_list[l+1]
        @assert -mmax_l ≤ m ≤ mmax_l "m out of range"
        ang_offset[l+1] + (m + mmax_l)
    end
    jtolm = j->begin
        # find largest l with base ≤ j
        l = searchsortedlast(ang_offset, j) - 1
        mmax_l = mmax_list[l+1]
        m = j - ang_offset[l+1] - mmax_l
        (l, m)
    end
    return lmtoj, jtolm, Nj, ang_offset
end

"""
    build_Qkernels(lmax; Lmax=2*lmax, maxm=-1, atol=0.0, threads=false, basis=:real)

Build M‑summed **angular kernels** for each multipole L=0..Lmax:

- `QL[L+1][α,α′,β,β′] = (4π/(2L+1)) Σ_M 𝒢^{(L)}_{α,α′;M} 𝒢̃^{(L)}_{β,β′;M}`
  This is the “direct” kernel with **bra/ket grouped per electron**.

`basis` controls the spherical-harmonic convention:
- `basis=:real` (default): real spherical harmonics; 𝒢̃≡𝒢.
- `basis=:complex`: complex spherical harmonics; the Coulomb addition theorem uses
  `Y_{LM}(Ω₁) Y*_{LM}(Ω₂)`, so 𝒢̃ uses an `M→-M` map with a phase `(-1)^M`, enforcing
  total-m conservation in Q.

- `QVL[L+1][α,β,α′,β′] = permutedims(QL[L+1], (1,3,2,4))`
  This is the 4‑index ordering convenient for applying V to a two‑electron amplitude
  `ψ[a,α,b,β]` via `multV(ψ, QVL, Veel)`.

- `QxL[L+1]` is returned in the same ordering as `QVL` for use by `getK`.
  `getK` uses a *crossed* packing of the density tensor, so the explicit “exchange”
  permutation lives in the `ρ` reshaping rather than in `QxL`.

Also returns `lmtoj`/`jtolm`/`Nj` for 0..lmax (optionally restricted by `maxm`),
and `NJ` for 0..Lmax.
"""
function build_Qkernels(lmax; Lmax=2*lmax, maxm::Int=-1, atol=0.0, threads=false, basis::Symbol=:real)
    basis in (:real, :complex) || throw(ArgumentError("basis must be :real or :complex"))
    Gtab = GauntTables.build(lmax; Lmax=Lmax, atol=atol, threads=threads, basis=basis)

    # indexers for one-electron Y_{ℓm} and multipole Y_{LM}
    lmtoj, _,  Nj,  _   = lm_indexers(lmax; maxm=maxm)
    LMtoj, _, NJ, _LM   = lm_indexers(Lmax)

    # Gpair[L+1][i, j, J] = 𝒢^{(L)}_{i j; M} with J≡(L,M), i≡(ℓ1,m1), j≡(ℓ2,m2)
    Gpair = [zeros(Nj, Nj, NJ) for _ in 0:Lmax]
    for L in 0:Lmax, l1 in 0:lmax, l2 in 0:lmax
        mmax1 = (maxm < 0) ? l1 : min(l1, maxm)
        mmax2 = (maxm < 0) ? l2 : min(l2, maxm)
        for e in GauntTables.block(Gtab, L, l1, l2)
            (abs(e.m1) <= mmax1 && abs(e.m2) <= mmax2) || continue
            J  = LMtoj(L, e.M)
            i  = lmtoj(l1, e.m1)
            j  = lmtoj(l2, e.m2)
            Gpair[L+1][i, j, J] = e.val
        end
    end

    # Q kernels (M-summed). One GEMM per L.
    QL  = [zeros(Nj, Nj, Nj, Nj) for _ in 0:Lmax]   # (α,α′,β,β′)
    QxL = [zeros(Nj, Nj, Nj, Nj) for _ in 0:Lmax]   # (α,β,α′,β′) for getK
    QVL = [zeros(Nj, Nj, Nj, Nj) for _ in 0:Lmax]   # (α,β,α′,β′) for multV

    for L in 0:Lmax
        pref = 4*pi/(2L+1)
        # rows = (i,j), cols = J  →  (Nj^2 × NJ)
        GLij = reshape(Gpair[L+1], Nj*Nj, NJ)
        if basis === :complex
            # For Coulomb: Y_{LM}(Ω₁) Y*_{LM}(Ω₂) = (-1)^M Y_{LM}(Ω₁) Y_{L,-M}(Ω₂)
            GLtil = zeros(eltype(GLij), size(GLij))
            for M in -L:L
                J    = LMtoj(L, M)
                Jbar = LMtoj(L, -M)
                phase = isodd(M) ? -1.0 : 1.0   # (-1)^M for integer M
                @views GLtil[:, J] .= phase .* GLij[:, Jbar]
            end
            Qdir = pref .* (GLij * GLtil')                   # (Nj^2 × Nj^2)
        else
            # Qdir[(i,j), (k,l)] = pref * Σ_J 𝒢(i,j;J) 𝒢(k,l;J)
            Qdir = pref .* (GLij * GLij')                    # (Nj^2 × Nj^2)
        end
        Q4   = reshape(Qdir, Nj, Nj, Nj, Nj)                # (α,α′,β,β′)

        QL[L+1]  = Q4                                       # (α,α′,β,β′)
        QxL[L+1] = permutedims(Q4, (1,3,2,4))               # (α,β,α′,β′)
        QVL[L+1] = permutedims(Q4, (1,3,2,4))               # (α,β,α′,β′)
    end

    # L = 0 must satisfy Q^{(0)}_{ij;kl} = δ_{ij} δ_{kl}
    Q0 = QL[1]
    for i in 1:Nj, j in 1:Nj, k in 1:Nj, l in 1:Nj
        expected = ((i == j) && (k == l)) ? 1.0 : 0.0
        @assert isapprox(Q0[i,j,k,l], expected; atol=1e-12)
    end

    return QL, QxL, QVL, lmtoj, Nj, NJ, Gtab
end

#=
function build_Qkernels(lmax; Lmax=2*lmax, atol=0.0, threads=false)
    G = GauntTables.build(lmax; Lmax=Lmax, atol=atol, threads=threads)

    # two independent (ℓ,m) packings: one for basis (Nj) and one for (L,M) (NJ)
    lmtoj, _,  Nj,  _  = lm_indexers(lmax)
    LMtoj, _, NJ, _LM  = lm_indexers(Lmax)

    # Gfull[L+1][i,k,J] stores G^{(L)}_{i,k;M} with J≡(L,M)
    Gfull = [zeros(Nj, Nj, NJ) for _ in 0:Lmax]
    for L in 0:Lmax, l1 in 0:lmax, l2 in 0:lmax
        for e in GauntTables.block(G, L, l1, l2)
            J  = LMtoj(L, e.M)
            i  = lmtoj(l1, e.m1)
            k  = lmtoj(l2, e.m2)
            Gfull[L+1][i, k, J] = e.val
        end
    end

    # Q^{(L)} and Qx^{(L)} (M-summed). One dense GEMM per L.
    QL  = [zeros(Nj, Nj, Nj, Nj) for _ in 0:Lmax]
    QxL = [zeros(Nj, Nj, Nj, Nj) for _ in 0:Lmax]
    QVL = [zeros(Nj, Nj, Nj, Nj) for _ in 0:Lmax]
    for L in 0:Lmax
        pref = 4*pi/(2L+1)
        GL   = reshape(Gfull[L+1], Nj*Nj, NJ)        # rows ≡ (i,k), cols ≡ J
        Qraw = GL * GL'                               # [(i,k) × (j,l)] with Σ_J
        Qraw .*= pref
        Q4   = reshape(Qraw, Nj, Nj, Nj, Nj)          # (i,k,j,l)
        # permute to (i,j,k,l) for the direct kernel
        QL[L+1]  = permutedims(Q4, (1,3,2,4))         # (i,j,k,l)
        QxL[L+1] = permutedims(QL[L+1], (1,4,3,2))    # (i,l,k,j)
        QVL[L+1]  = permutedims(Q4, (1,3,2,4))         # (i,j,k,l)
    end

    return QL, QxL, QVL, lmtoj, Nj, NJ, G
end
=#

# ---- Contractions (explicit shapes; DA semantics) --------------------------

"""
    getJ(rho, QL, VL)

Direct Coulomb under the **local DA**.

Inputs:
- `rho`  :: Array (Nr, Nj, Nr, Nj)   with entries ρ[a,i,b,j] = γ_{(a,i),(b,j)}
- `QL`   :: Vector{Array{…}} length Lmax+1; each QL[L+1] is (Nj,Nj,Nj,Nj)
- `VL`   :: Vector{Matrix} length Lmax+1; each VL[L+1] is (Nr,Nr)

Returns:
- `J`    :: Array (Nr, Nj, Nj) giving J[a, i, j] (nonzero only on left radial block a=a)
"""
function getJ(rho, QL, VL)
    Nr, Nj = size(rho,1), size(rho,2)
    # collect γ_{(c,ℓ),(c,k)} for each c into columns (Nj^2 × Nr)
    rhod = zeros(Nj*Nj, Nr)
    @inbounds for c in 1:Nr
        rhod[:, c] = vec(@view rho[c, 1:Nj, c, 1:Nj])
    end
    acc = zeros(Nr, Nj*Nj)
    @inbounds for Lp in eachindex(QL)                     # Lp = L+1
        Qmat = reshape(QL[Lp], Nj*Nj, Nj*Nj)              # (ij, kl)
        acc .+= VL[Lp] * (Qmat * rhod)'                   # Nr×Nr * (Nr×Nj^2)ᵗ → Nr×Nj^2
    end
    reshape(acc, Nr, Nj, Nj)                              # J[a, i, j]
end

"""
    getK(rho, QxL, VL)

Exchange (two-site on the left).

Inputs/returns as in `getJ`, but output is full 4D:
- `K` :: Array (Nr, Nj, Nr, Nj) giving K[a, i, b, j]
"""
function getK(rho, QxL, VL)
    Nr, Nj = size(rho,1), size(rho,2)
    # ρ̃[(k,ℓ), a, b] = γ_{(b,ℓ),(a,k)} reshaped as a matrix (Nj^2 × Nr^2)
    rhotil = reshape(permutedims(rho, (2,4,1,3)), Nj*Nj, Nr*Nr)
    acc = zeros(Nj*Nj, Nr, Nr)
    @inbounds for Lp in eachindex(QxL)
        Qmat = reshape(QxL[Lp], Nj*Nj, Nj*Nj)             # (i j, k ℓ)
        tmp  = reshape(Qmat * rhotil, Nj*Nj, Nr, Nr)      # (i j, a, b)
        # scale by VL[L][a,b] per (a,b) (diagonal approx in radial block)
        for a in 1:Nr, b in 1:Nr
            @views tmp[:, a, b] .*= VL[Lp][a, b]
        end
        acc .+= tmp
    end
    # reshape back to (a,i,b,j)
    permutedims(reshape(acc, Nj, Nj, Nr, Nr), (3,1,4,2))
end

"""
    multV(psi, QVL, VL)

Apply V to a two-electron amplitude ψ with annihilators on the right:
- `psi[a,i,b,j]` corresponds to |(a,i);(b,j)⟩ right indices (i,j).
Returns an array with the same shape as `psi`.
"""
function multV(psi, QVL, VL)
    Nr, Nj = size(psi,1), size(psi,2)
    # ψ̃[(k,ℓ), a, b] = ψ[a, k, b, ℓ]
    psitil = reshape(permutedims(psi, (2,4,1,3)), Nj*Nj, Nr * Nr)
    acc = zeros(Nj*Nj, Nr, Nr)
    @inbounds for Lp in eachindex(QVL)
        Qmat = reshape(QVL[Lp], Nj*Nj, Nj*Nj)              # (i j, k ℓ)
        tmp  = reshape(Qmat * psitil,Nj*Nj,Nr,Nr)
        for a in 1:Nr, b in 1:Nr
            @views tmp[:, a, b] .*= VL[Lp][a, b]
        end
        acc .+= tmp
    end
    permutedims(reshape(acc, Nj, Nj, Nr, Nr), (3,1,4,2))
end

"""
    build_V6_tensor(Veel, QL, Nr, lmax; maxm=-1)

Materialize the six-index Coulomb tensor
V(a,α,α′, b,β,β′) with α=(ℓ,m) packed over 0≤ℓ≤lmax, m=-mmax(ℓ):mmax(ℓ).

Inputs
- `Veel` :: Vector of radial multipole matrices (length Lmax+1)
- `QL`   :: Vector of DiagYlm kernels (length Lmax+1)
- `Nr`   :: Number of radial functions
- `lmax` :: Maximum orbital angular momentum included
- `maxm` :: Max |m| retained per ℓ (negative -> full m)

Returns `(V6, ang_offset)` where `V6` has shape `(Nr, Nang, Nang, Nr, Nang, Nang)`.
The tensor is symmetric by construction (no post symmetrization).
"""
function build_V6_tensor(Veel, QL, Nr::Int, lmax::Int; maxm::Int=-1)
    _, _, nang, ang_offset = lm_indexers(lmax; maxm=maxm)
    ang_offset = vcat(ang_offset, nang + 1)

    V6 = zeros(Float64, Nr, nang, nang, Nr, nang, nang)
    for Lp in eachindex(QL)
        VL = Veel[Lp]
        Q = QL[Lp]
        for a in 1:Nr, b in 1:Nr
            pref = VL[a,b]
            pref == 0.0 && continue
            for α in 1:nang, β in 1:nang, αp in 1:nang, βp in 1:nang
                V6[a, α, αp, b, β, βp] += pref * Q[α, αp, β, βp]
            end
        end
    end

    return V6, ang_offset
end

end # module
