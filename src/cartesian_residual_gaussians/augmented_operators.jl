symmetrize_operator(matrix) = Matrix{Float64}(0.5 .* (matrix .+ transpose(matrix)))

function transform_augmented_operator(O_GG, O_GA, O_AA, residual)
    T_G, T_A = residual.T_G, residual.T_A
    O_GR = O_GG * T_G + O_GA * T_A
    O_RR = transpose(T_G) * O_GG * T_G +
           transpose(T_G) * O_GA * T_A +
           transpose(T_A) * transpose(O_GA) * T_G +
           transpose(T_A) * O_AA * T_A
    nG, nR = residual.base_dimension, residual.residual_dimension
    out = zeros(Float64, nG + nR, nG + nR)
    residual_range = (nG + 1):(nG + nR)
    out[1:nG, 1:nG] .= O_GG
    out[1:nG, residual_range] .= O_GR
    out[residual_range, 1:nG] .= transpose(O_GR)
    out[residual_range, residual_range] .= O_RR
    return symmetrize_operator(out)
end
