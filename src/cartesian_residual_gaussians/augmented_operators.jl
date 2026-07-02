symmetrize_operator(matrix) = Matrix{Float64}(0.5 .* (matrix .+ transpose(matrix)))
_augmented_operator_block(O_GG, O_GA, O_AA, L_G, L_A, R_G, R_A) = transpose(L_G) * O_GG * R_G + transpose(L_G) * O_GA * R_A + transpose(L_A) * transpose(O_GA) * R_G + transpose(L_A) * O_AA * R_A

function transform_augmented_operator(O_GG, O_GA, O_AA, residual)
    T_G, T_A = residual.T_G, residual.T_A
    if isnothing(residual.injected_G)
        O_FF = O_GG; O_FR = O_GG * T_G + O_GA * T_A
    else
        Y, Qp = residual.injected_A, injection_complement(residual)
        O_YY, O_YQ, O_QQ = transpose(Y) * O_AA * Y, transpose(Y) * transpose(O_GA) * Qp, transpose(Qp) * O_GG * Qp
        O_FF = [O_YY O_YQ; transpose(O_YQ) O_QQ]
        O_FR = vcat(transpose(Y) * (transpose(O_GA) * T_G + O_AA * T_A),
            transpose(Qp) * (O_GG * T_G + O_GA * T_A))
    end
    O_RR = _augmented_operator_block(O_GG, O_GA, O_AA, T_G, T_A, T_G, T_A)
    nG, nR = residual.base_dimension, residual.residual_dimension
    out = zeros(Float64, nG + nR, nG + nR)
    residual_range = (nG + 1):(nG + nR)
    out[1:nG, 1:nG] .= O_FF
    out[1:nG, residual_range] .= O_FR
    out[residual_range, 1:nG] .= transpose(O_FR)
    out[residual_range, residual_range] .= O_RR
    return symmetrize_operator(out)
end

function protected_original_fixed_sector_components(geometry)
    T_G, T_A, Z, B = geometry.T_G, geometry.T_A, geometry.Z, geometry.B
    nG, nR, nZ, nM = size(T_G, 1), size(T_A, 2), size(Z, 2), size(B, 1)
    size(T_G, 2) == nR || throw(DimensionMismatch("protected fixed-sector T_G/T_A size mismatch"))
    size(B, 2) == nZ || throw(DimensionMismatch("protected fixed-sector B/Z size mismatch"))
    nM == nG + nR || throw(DimensionMismatch("protected fixed-sector M dimension mismatch"))
    qfull = Matrix(qr(B).Q * Matrix{Float64}(I, nM, nM))
    Qp = Matrix{Float64}(view(qfull, :, (nZ + 1):nM))
    Qp_G = view(Qp, 1:nG, :)
    Qp_R = view(Qp, (nG + 1):nM, :)
    G_perp = Matrix{Float64}(Qp_G) + T_G * Qp_R
    A_perp = T_A * Qp_R
    return (; Z, Qp, G_perp, A_perp,
        protected_count = hasproperty(geometry, :Z_protected) ?
            size(geometry.Z_protected, 2) : getproperty(geometry, :protected_original_count),
        z_dimension = nZ, f_dimension = nM)
end

function transform_protected_original_fixed_sector_operator(O_GG, O_GA, O_AA, components)
    Z, Gp, Ap = components.Z, components.G_perp, components.A_perp
    O_ZZ = transpose(Z) * O_AA * Z
    O_ZQ = transpose(Z) * (transpose(O_GA) * Gp + O_AA * Ap)
    O_QQ = _augmented_operator_block(O_GG, O_GA, O_AA, Gp, Ap, Gp, Ap)
    nZ, nQ = size(Z, 2), size(Gp, 2)
    out = zeros(Float64, nZ + nQ, nZ + nQ)
    qrange = (nZ + 1):(nZ + nQ)
    out[1:nZ, 1:nZ] .= O_ZZ
    out[1:nZ, qrange] .= O_ZQ
    out[qrange, 1:nZ] .= transpose(O_ZQ)
    out[qrange, qrange] .= O_QQ
    return symmetrize_operator(out)
end

function transform_protected_original_fixed_sector_one_body(kinetic, unit_nuclear_by_center,
    nuclear_charges, geometry)
    components = protected_original_fixed_sector_components(geometry)
    K = transform_protected_original_fixed_sector_operator(
        kinetic.GG, kinetic.GA, kinetic.AA, components)
    U = Matrix{Float64}[transform_protected_original_fixed_sector_operator(
        unit.GG, unit.GA, unit.AA, components) for unit in unit_nuclear_by_center]
    length(U) == length(nuclear_charges) ||
        throw(DimensionMismatch("protected fixed-sector nuclear charge count mismatch"))
    H1 = copy(K)
    for (charge, unit) in zip(nuclear_charges, U)
        H1 .+= Float64(charge) .* unit
    end
    return (; kinetic = K, nuclear_attraction_unit_by_center = U,
        one_body_hamiltonian = symmetrize_operator(H1),
        protected_count = components.protected_count,
        z_dimension = components.z_dimension,
        f_dimension = components.f_dimension)
end
