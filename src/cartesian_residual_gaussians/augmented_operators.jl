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
