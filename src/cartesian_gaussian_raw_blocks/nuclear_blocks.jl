const _GB_PARENT = parentmodule(@__MODULE__)

_float_center(center) = (Float64(center[1]), Float64(center[2]), Float64(center[3]))
_symmetrize_raw_block(matrix) = Matrix{Float64}(0.5 .* (matrix .+ transpose(matrix)))
_axis_coordinate(orbital, axis) =
    axis === :x ? orbital.center[1] :
    axis === :y ? orbital.center[2] :
    axis === :z ? orbital.center[3] :
    throw(ArgumentError("axis must be :x, :y, or :z"))
_axis_power(orbital, axis) =
    axis === :x ? orbital.lx :
    axis === :y ? orbital.ly :
    axis === :z ? orbital.lz :
    throw(ArgumentError("axis must be :x, :y, or :z"))

function _axis_factor_inputs(orbital, axis)
    parent = _GB_PARENT
    exponents = Float64.(orbital.exponents)
    powers = fill(_axis_power(orbital, axis), length(exponents))
    return (;
        exponents,
        centers = fill(Float64(_axis_coordinate(orbital, axis)), length(exponents)),
        powers,
        prefactors = getfield(parent, :_cartesian_gaussian_axis_prefactors)(
            exponents, powers),
    )
end

function _unique_axis_centers(centers, axis)
    values = Float64[]
    lookup = Int[]
    for center in centers
        value = center[axis]
        index = findfirst(==(value), values)
        if isnothing(index)
            push!(values, value)
            push!(lookup, length(values))
        else
            push!(lookup, index)
        end
    end
    return values, lookup
end

function _fill_axis_factor_table!(
    destination,
    left,
    right,
    factor_exponent::Float64,
    factor_center::Float64,
)
    integral = getfield(_GB_PARENT, :_cartesian_gaussian_axis_integral)
    @inbounds for column in axes(destination, 2), row in axes(destination, 1)
        destination[row, column] = integral(
            left.exponents[row],
            left.centers[row],
            left.powers[row],
            left.prefactors[row],
            right.exponents[column],
            right.centers[column],
            right.powers[column],
            right.prefactors[column],
            :factor;
            factor_exponent,
            factor_center,
        )
    end
    return destination
end

function _weighted_hadamard3(
    left_coefficients,
    x,
    y,
    z,
    right_coefficients,
)
    value = 0.0
    @inbounds for column in eachindex(right_coefficients), row in eachindex(left_coefficients)
        value += Float64(left_coefficients[row]) * Float64(right_coefficients[column]) *
            x[row, column] * y[row, column] * z[row, column]
    end
    return value
end

function _accumulate_aa_nuclear_pair!(
    aa,
    left_index::Int,
    right_index::Int,
    left,
    right,
    expansion,
    centers,
)
    left_axes = (_axis_factor_inputs(left, :x), _axis_factor_inputs(left, :y),
        _axis_factor_inputs(left, :z))
    right_axes = (_axis_factor_inputs(right, :x), _axis_factor_inputs(right, :y),
        _axis_factor_inputs(right, :z))
    unique_centers = ntuple(axis -> _unique_axis_centers(centers, axis), 3)
    tables = ntuple(axis -> [zeros(Float64, length(left.coefficients),
        length(right.coefficients)) for _ in unique_centers[axis][1]], 3)

    values = zeros(Float64, length(centers))
    for (term_index, coefficient) in pairs(expansion.coefficients)
        exponent = Float64(expansion.exponents[term_index])
        for axis in 1:3, (center_index, center) in pairs(unique_centers[axis][1])
            _fill_axis_factor_table!(
                tables[axis][center_index], left_axes[axis], right_axes[axis],
                exponent, center)
        end
        for center_index in eachindex(centers)
            ix = unique_centers[1][2][center_index]
            iy = unique_centers[2][2][center_index]
            iz = unique_centers[3][2][center_index]
            values[center_index] -= coefficient * _weighted_hadamard3(
                left.coefficients, tables[1][ix], tables[2][iy], tables[3][iz],
                right.coefficients)
        end
    end
    for (center_index, value) in pairs(values)
        aa[center_index][left_index, right_index] = value
        aa[center_index][right_index, left_index] = value
    end
    return nothing
end

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

    for left_index in eachindex(supplement.orbitals)
        left = supplement.orbitals[left_index]
        for right_index in left_index:length(supplement.orbitals)
            right = supplement.orbitals[right_index]
            _accumulate_aa_nuclear_pair!(
                aa, left_index, right_index, left, right, expansion, centers)
        end
    end
    return (; ga, aa = [_symmetrize_raw_block(matrix) for matrix in aa])
end
