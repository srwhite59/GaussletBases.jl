const _GB_PARENT = parentmodule(@__MODULE__)

_float_center(center) = (Float64(center[1]), Float64(center[2]), Float64(center[3]))
_symmetrize_raw_block(matrix) = Matrix{Float64}(0.5 .* (matrix .+ transpose(matrix)))
_axis_layer(proxy, axis::Int) = axis == 1 ? proxy.x : axis == 2 ? proxy.y : proxy.z
_axis_coordinate(orbital, axis) =
    axis == 1 || axis === :x ? orbital.center[1] :
    axis == 2 || axis === :y ? orbital.center[2] :
    axis == 3 || axis === :z ? orbital.center[3] :
    throw(ArgumentError("axis must be :x, :y, or :z"))
_axis_power(orbital, axis) =
    axis == 1 || axis === :x ? orbital.lx :
    axis == 2 || axis === :y ? orbital.ly :
    axis == 3 || axis === :z ? orbital.lz :
    throw(ArgumentError("axis must be :x, :y, or :z"))

function _proxy_axis_factor_inputs(proxy_layer)
    parent = _GB_PARENT
    proxy_gaussians = getfield(parent, :_qwrg_proxy_gaussian_primitives)(proxy_layer)
    gaussian_exponent = getfield(parent.GaussianAnalyticIntegrals, :gaussian_exponent)
    return (; exponents = Float64[gaussian_exponent(g) for g in proxy_gaussians],
        centers = Float64[gaussian.center_value for gaussian in proxy_gaussians],
        powers = zeros(Int, length(proxy_gaussians)),
        prefactors = ones(Float64, length(proxy_gaussians)))
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

function _factor_axis_integral(left_exponent, left_center, left_power, left_prefactor,
    right_exponent, right_center, right_power, right_prefactor, factor_exponent,
    factor_center)
    gamma = left_exponent + right_exponent + factor_exponent
    weighted_center = (left_exponent * left_center + right_exponent * right_center +
        factor_exponent * factor_center) / gamma
    constant = (
        left_exponent * right_exponent * (left_center - right_center)^2 +
        left_exponent * factor_exponent * (left_center - factor_center)^2 +
        right_exponent * factor_exponent * (right_center - factor_center)^2
    ) / gamma
    left_shift = weighted_center - left_center
    right_shift = weighted_center - right_center
    value = 0.0
    for i in 0:left_power, j in 0:right_power
        moment_power = i + j
        isodd(moment_power) && continue
        moment = sqrt(pi / gamma)
        for k in 1:div(moment_power, 2)
            moment *= (2 * k - 1) / (2 * gamma)
        end
        value += binomial(left_power, i) * binomial(right_power, j) *
            left_shift^(left_power - i) * right_shift^(right_power - j) * moment
    end
    return left_prefactor * right_prefactor * exp(-constant) * value
end

function _fill_axis_factor_table!(destination, left, right, factor_exponent::Float64,
    factor_center::Float64)
    @inbounds for column in axes(destination, 2), row in axes(destination, 1)
        destination[row, column] = _factor_axis_integral(
            left.exponents[row],
            left.centers[row],
            left.powers[row],
            left.prefactors[row],
            right.exponents[column],
            right.centers[column],
            right.powers[column],
            right.prefactors[column],
            factor_exponent,
            factor_center,
        )
    end
    return destination
end

function _weighted_hadamard3(left_coefficients, x, x_reversed, y, y_reversed, z,
    z_reversed, right_coefficients)
    value = 0.0
    @inbounds for column in eachindex(right_coefficients), row in eachindex(left_coefficients)
        value += Float64(left_coefficients[row]) * Float64(right_coefficients[column]) *
            (x_reversed ? x[column, row] : x[row, column]) *
            (y_reversed ? y[column, row] : y[row, column]) *
            (z_reversed ? z[column, row] : z[row, column])
    end
    return value
end

function _axis_family_inventory(supplement)
    maps = [zeros(Int, length(supplement.orbitals)) for _ in 1:3]
    families = [Any[] for _ in 1:3]
    for axis in 1:3
        lookup = Dict{Any,Int}()
        for (orbital_index, orbital) in pairs(supplement.orbitals)
            exponents = Float64.(orbital.exponents)
            powers = fill(_axis_power(orbital, axis), length(exponents))
            prefactors = getfield(_GB_PARENT, :_cartesian_gaussian_axis_prefactors)(
                exponents, powers)
            key = (Tuple(exponents), Float64(_axis_coordinate(orbital, axis)),
                powers[1], Tuple(prefactors))
            maps[axis][orbital_index] = get!(lookup, key) do
                push!(families[axis], (; exponents,
                    centers = fill(Float64(_axis_coordinate(orbital, axis)),
                        length(exponents)),
                    powers, prefactors))
                length(families[axis])
            end
        end
    end
    return (; maps, families)
end

function _make_ga_tables(proxy, inventory, axis_centers)
    stencils = [Matrix{Float64}(_GB_PARENT.stencil_matrix(_axis_layer(proxy, axis)))
        for axis in 1:3]
    proxy_inputs = [_proxy_axis_factor_inputs(_axis_layer(proxy, axis)) for axis in 1:3]
    tables = [Matrix{Float64}[zeros(Float64, size(stencils[axis], 2),
        length(family.exponents)) for family in inventory.families[axis],
        _ in axis_centers[axis]] for axis in 1:3]
    primitive = [Matrix{Float64}[zeros(Float64, length(proxy_inputs[axis].exponents),
        length(family.exponents)) for family in inventory.families[axis],
        _ in axis_centers[axis]] for axis in 1:3]
    return (; stencils, proxy_inputs, tables, primitive)
end

function _make_aa_tables(inventory, axis_centers)
    tables = [Matrix{Float64}[] for _ in 1:3]
    specs = [Tuple{Int,Int,Int}[] for _ in 1:3]
    ids = [zeros(Int, length(inventory.families[axis]),
        length(inventory.families[axis]), length(axis_centers[axis])) for axis in 1:3]
    for axis in 1:3, left in eachindex(inventory.families[axis]),
            right in left:length(inventory.families[axis]), center in eachindex(axis_centers[axis])
        push!(tables[axis], zeros(Float64, length(inventory.families[axis][left].exponents),
            length(inventory.families[axis][right].exponents)))
        push!(specs[axis], (left, right, center))
        id = length(tables[axis])
        ids[axis][left, right, center] = id
        ids[axis][right, left, center] = id
    end
    return (; tables, specs, ids)
end

function gaussian_nuclear_raw_blocks_by_center(proxy, supplement, expansion, atom_locations)
    centers = [_float_center(center) for center in atom_locations]
    ncart, norbital = proxy.ncart, length(supplement.orbitals)
    ga = [zeros(Float64, ncart, norbital) for _ in centers]
    aa = [zeros(Float64, norbital, norbital) for _ in centers]

    parent = _GB_PARENT
    fill_product! = getfield(parent, :_qwrg_fill_product_column!)
    left_contract! = getfield(parent, :mul!)
    inventory = _axis_family_inventory(supplement)
    axis_centers = [first(_unique_axis_centers(centers, axis)) for axis in 1:3]
    center_lookup = [last(_unique_axis_centers(centers, axis)) for axis in 1:3]
    ga_tables = _make_ga_tables(proxy, inventory, axis_centers)
    aa_tables = _make_aa_tables(inventory, axis_centers)
    ga_plans = [(orbital, center, (inventory.maps[1][orbital],
        inventory.maps[2][orbital], inventory.maps[3][orbital]),
        (center_lookup[1][center], center_lookup[2][center],
            center_lookup[3][center])) for orbital in 1:norbital for center in eachindex(centers)]
    aa_plans = [(left, right, center,
        ntuple(axis -> aa_tables.ids[axis][inventory.maps[axis][left],
            inventory.maps[axis][right], center_lookup[axis][center]], 3),
        ntuple(axis -> inventory.maps[axis][left] > inventory.maps[axis][right], 3))
        for left in 1:norbital for right in left:norbital for center in eachindex(centers)]
    scratch = zeros(Float64, ncart)

    for (term_index, coefficient) in pairs(expansion.coefficients)
        exponent = Float64(expansion.exponents[term_index])
        for axis in 1:3
            for family in eachindex(inventory.families[axis]), center in eachindex(axis_centers[axis])
                _fill_axis_factor_table!(ga_tables.primitive[axis][family, center],
                    ga_tables.proxy_inputs[axis], inventory.families[axis][family],
                    exponent, axis_centers[axis][center])
                left_contract!(ga_tables.tables[axis][family, center],
                    transpose(ga_tables.stencils[axis]),
                    ga_tables.primitive[axis][family, center])
            end
            for (id, spec) in pairs(aa_tables.specs[axis])
                _fill_axis_factor_table!(aa_tables.tables[axis][id],
                    inventory.families[axis][spec[1]],
                    inventory.families[axis][spec[2]], exponent,
                    axis_centers[axis][spec[3]])
            end
        end

        for (orbital_index, center_index, families, center_ids) in ga_plans
            orbital = supplement.orbitals[orbital_index]
            column = view(ga[center_index], :, orbital_index)
            factors = ntuple(axis -> ga_tables.tables[axis][families[axis],
                center_ids[axis]], 3)
            for primitive in eachindex(orbital.coefficients)
                fill_product!(scratch,
                    view(factors[1], :, primitive),
                    view(factors[2], :, primitive),
                    view(factors[3], :, primitive))
                scale = coefficient * Float64(orbital.coefficients[primitive])
                @inbounds for row in eachindex(scratch)
                    column[row] -= scale * scratch[row]
                end
            end
        end

        for (left_index, right_index, center_index, ids, reversed) in aa_plans
            left = supplement.orbitals[left_index]
            right = supplement.orbitals[right_index]
            value = _weighted_hadamard3(left.coefficients,
                aa_tables.tables[1][ids[1]], reversed[1],
                aa_tables.tables[2][ids[2]], reversed[2],
                aa_tables.tables[3][ids[3]], reversed[3],
                right.coefficients)
            aa[center_index][left_index, right_index] -= coefficient * value
        end
    end
    for matrix in aa, row in 1:norbital, column in (row + 1):norbital
        matrix[column, row] = matrix[row, column]
    end
    return (; ga, aa = [_symmetrize_raw_block(matrix) for matrix in aa])
end
