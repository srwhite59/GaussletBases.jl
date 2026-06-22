const _GB_PARENT = parentmodule(@__MODULE__)

_float_center(center) = (Float64(center[1]), Float64(center[2]), Float64(center[3]))
_symmetrize_raw_block(matrix) = Matrix{Float64}(0.5 .* (matrix .+ transpose(matrix)))

function gaussian_nuclear_raw_blocks_by_center(proxy, supplement, expansion, atom_locations)
    centers = [_float_center(center) for center in atom_locations]
    ncart, norbital = proxy.ncart, length(supplement.orbitals)
    ga = [zeros(Float64, ncart, norbital) for _ in centers]
    aa = [zeros(Float64, norbital, norbital) for _ in centers]

    parent = _GB_PARENT
    fill_product! = getfield(parent, :_qwrg_fill_product_column!)
    axis_cross = getfield(parent, :_qwrg_atomic_axis_factor_cross_data)
    cross_cache = Dict{Tuple{Int,Symbol,Float64},Any}()
    scratch = zeros(Float64, ncart)
    for (orbital_index, orbital) in pairs(supplement.orbitals)
        for (center_index, center) in pairs(centers)
            factors = (
                get!(cross_cache, (orbital_index, :x, center[1])) do
                    axis_cross(proxy.x, orbital, :x, expansion, center[1])
                end,
                get!(cross_cache, (orbital_index, :y, center[2])) do
                    axis_cross(proxy.y, orbital, :y, expansion, center[2])
                end,
                get!(cross_cache, (orbital_index, :z, center[3])) do
                    axis_cross(proxy.z, orbital, :z, expansion, center[3])
                end,
            )
            column = view(ga[center_index], :, orbital_index)
            for term in eachindex(expansion.coefficients),
                    primitive in eachindex(orbital.coefficients)
                fill_product!(scratch,
                    view(factors[1][term], :, primitive),
                    view(factors[2][term], :, primitive),
                    view(factors[3][term], :, primitive))
                scale = expansion.coefficients[term] *
                    Float64(orbital.coefficients[primitive])
                @inbounds for row in eachindex(scratch)
                    column[row] -= scale * scratch[row]
                end
            end
        end
    end

    axis_aa = getfield(parent, :_qwrg_atomic_axis_factor_aa_data)
    weighted = getfield(parent, :_qwrg_atomic_weighted_hadamard)
    for left_index in eachindex(supplement.orbitals)
        left = supplement.orbitals[left_index]
        for right_index in left_index:length(supplement.orbitals)
            right = supplement.orbitals[right_index]
            local_cache = Dict{Tuple{Symbol,Float64},Any}()
            for (center_index, center) in pairs(centers)
                fx = get!(local_cache, (:x, center[1])) do
                    axis_aa(left, right, :x, expansion, center[1])
                end
                fy = get!(local_cache, (:y, center[2])) do
                    axis_aa(left, right, :y, expansion, center[2])
                end
                fz = get!(local_cache, (:z, center[3])) do
                    axis_aa(left, right, :z, expansion, center[3])
                end
                value = 0.0
                for term in eachindex(expansion.coefficients)
                    value -= expansion.coefficients[term] *
                        weighted(left.coefficients, fx[term], fy[term], fz[term],
                            right.coefficients)
                end
                aa[center_index][left_index, right_index] = value
                aa[center_index][right_index, left_index] = value
            end
        end
    end
    return (; ga, aa = [_symmetrize_raw_block(matrix) for matrix in aa])
end
