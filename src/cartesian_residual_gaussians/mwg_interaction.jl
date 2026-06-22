const _RG_PARENT = parentmodule(@__MODULE__)
function moment_matched_gaussians(operators, residual)
    nG, nR = residual.base_dimension, residual.residual_dimension
    centers, widths = zeros(Float64, nR, 3), zeros(Float64, nR, 3)
    for (axis_index, axis) in pairs((:x, :y, :z))
        position, second = operators.position[axis], operators.x2[axis]
        for residual_index in 1:nR
            row = nG + residual_index
            center = Float64(position[row, row])
            variance = Float64(second[row, row]) - center^2
            isfinite(center) && isfinite(variance) && variance > 0.0 || throw(
                ArgumentError("residual MWG moment variance must be finite positive"))
            width = sqrt(2.0 * variance)
            isfinite(width) && width > 0.0 || throw(ArgumentError("residual MWG width must be finite positive"))
            centers[residual_index, axis_index] = center; widths[residual_index, axis_index] = width
        end
    end
    return centers, widths
end
function _mwg_axis_pairs(bundles, expansion, residual_centers, residual_widths)
    bundle(axis) = axis == 1 ? bundles.bundle_x : axis == 2 ? bundles.bundle_y : bundles.bundle_z
    effective = getfield(_RG_PARENT, :_qwrg_effective_gaussians)(residual_centers, residual_widths)
    split = ntuple(axis -> getfield(_RG_PARENT, :_qwrg_split_block_matrices)(
        bundle(axis), effective[axis], expansion), 3)
    analytic = ntuple(axis -> getfield(_RG_PARENT, :_qwrg_gaussian_analytic_blocks)(effective[axis], expansion), 3)
    weights = ntuple(axis -> getfield(_RG_PARENT, :_qwrg_supplement_integral_weights)(effective[axis]), 3)
    normalize = getfield(_RG_PARENT, :_qwrg_density_normalized_pair_matrices)
    return (; ga = ntuple(axis -> normalize(split[axis].pair_ga, split[axis].weight_gg,
        weights[axis]; label = "R3-B MWG $((:x, :y, :z)[axis]) G-M"), 3),
        aa = ntuple(axis -> normalize(analytic[axis].pair_aa, weights[axis],
        weights[axis]; label = "R3-B MWG $((:x, :y, :z)[axis]) M-M"), 3))
end
function _mwg_support_weights(states, bundles)
    axis_pgdg = getfield(_RG_PARENT, :_nested_axis_pgdg)
    wx, wy, wz = axis_pgdg(bundles, :x).weights, axis_pgdg(bundles, :y).weights, axis_pgdg(bundles, :z).weights
    return [wx[s[1]] * wy[s[2]] * wz[s[3]] for s in states]
end
function _terminal_mwg_fixed_residual(basis, bundles, pair_terms, coefficients)
    nR = size(first(pair_terms.ga[1]), 2)
    V_GM = zeros(Float64, basis.final_dimension, nR)
    for block in basis.blocks
        local_values = zeros(Float64, length(block.support_states), nR)
        for residual in 1:nR
            column = view(local_values, :, residual)
            @inbounds for (row, state) in pairs(block.support_states)
                ix, iy, iz = state
                value = 0.0
                for term in eachindex(coefficients)
                    value += Float64(coefficients[term]) * pair_terms.ga[1][term][ix, residual] *
                        pair_terms.ga[2][term][iy, residual] * pair_terms.ga[3][term][iz, residual]
                end
                column[row] = value
            end
        end
        if block.coefficients === nothing
            V_GM[block.column_range, :] .= local_values
        else
            support_weights = _mwg_support_weights(block.support_states, bundles)
            final_weights = vec(transpose(block.coefficients) * support_weights)
            all(weight -> isfinite(weight) && weight > 1.0e-12, final_weights) ||
                throw(ArgumentError("residual MWG final density weights must be finite positive"))
            density_coefficients = block.coefficients .* reshape(support_weights, :, 1) ./ reshape(final_weights, 1, :)
            mul!(view(V_GM, block.column_range, :), transpose(density_coefficients), local_values)
        end
    end
    return V_GM
end
function _mwg_residual_residual(pair_terms, coefficients)
    nR = size(first(pair_terms.aa[1]), 1)
    V_MM = zeros(Float64, nR, nR)
    for i in 1:nR, j in i:nR
        value = sum(Float64(coefficients[term]) * pair_terms.aa[1][term][i, j] *
            pair_terms.aa[2][term][i, j] * pair_terms.aa[3][term][i, j] for term in eachindex(coefficients))
        V_MM[i, j] = value; V_MM[j, i] = value
    end
    return V_MM
end
_mwg_default_expansion(expansion) = isnothing(expansion) ? getfield(_RG_PARENT, :coulomb_gaussian_expansion)(doacc = false) : throw(ArgumentError("residual MWG interaction rejects custom expansion because base V_GG has no expansion provenance"))
function assemble_residual_ida_interaction(base_V_GG, basis, bundles, residual, augmented_operators; expansion = nothing)
    nG, nR = residual.base_dimension, residual.residual_dimension
    size(base_V_GG) == (nG, nG) || throw(DimensionMismatch("residual MWG base V_GG dimension mismatch"))
    centers, widths = moment_matched_gaussians(augmented_operators, residual)
    expansion_value = _mwg_default_expansion(expansion)
    pair_terms = _mwg_axis_pairs(bundles, expansion_value, centers, widths)
    V_GM = _terminal_mwg_fixed_residual(basis, bundles, pair_terms, expansion_value.coefficients)
    V_MM = _mwg_residual_residual(pair_terms, expansion_value.coefficients)
    V = zeros(Float64, nG + nR, nG + nR); residual_range = (nG + 1):(nG + nR)
    V[1:nG, 1:nG] .= base_V_GG; V[1:nG, residual_range] .= V_GM
    V[residual_range, 1:nG] .= transpose(V_GM); V[residual_range, residual_range] .= V_MM
    return V
end
