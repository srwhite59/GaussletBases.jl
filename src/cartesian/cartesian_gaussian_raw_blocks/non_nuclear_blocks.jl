function _non_nuclear_overlap_axis_data(proxy_layer, orbital, axis::Symbol)
    parent = _GB_PARENT
    proxy_gaussians = parent._qwrg_proxy_gaussian_primitives(proxy_layer)
    gaussian_exponent = getfield(parent.GaussianAnalyticIntegrals, :gaussian_exponent)
    center_value = axis == :x ? orbital.center[1] :
        axis == :y ? orbital.center[2] : orbital.center[3]
    power = axis == :x ? orbital.lx : axis == :y ? orbital.ly : orbital.lz
    nproxy = length(proxy_gaussians)
    nprimitive = length(orbital.exponents)
    left_exponents = Float64[gaussian_exponent(gaussian) for gaussian in proxy_gaussians]
    right_exponents = Float64.(orbital.exponents)
    table = parent._cartesian_gaussian_axis_integral_table(
        left_exponents,
        Float64[gaussian.center_value for gaussian in proxy_gaussians],
        zeros(Int, nproxy),
        ones(Float64, nproxy),
        right_exponents,
        fill(Float64(center_value), nprimitive),
        fill(power, nprimitive),
        parent._cartesian_gaussian_axis_prefactors(
            right_exponents, fill(power, nprimitive)),
        :overlap,
    )
    return parent._qwrg_left_contract_cross_matrix(proxy_layer, table)
end

function gaussian_non_nuclear_overlap_blocks(proxy, supplement)
    ncart, norbital = proxy.ncart, length(supplement.orbitals)
    overlap_ga = zeros(Float64, ncart, norbital)
    scratch = zeros(Float64, ncart)
    fill_product! = getfield(_GB_PARENT, :_qwrg_fill_product_column!)
    for (orbital_index, orbital) in pairs(supplement.orbitals)
        x_data = _non_nuclear_overlap_axis_data(proxy.x, orbital, :x)
        y_data = _non_nuclear_overlap_axis_data(proxy.y, orbital, :y)
        z_data = _non_nuclear_overlap_axis_data(proxy.z, orbital, :z)
        column = view(overlap_ga, :, orbital_index)
        for primitive in eachindex(orbital.coefficients)
            fill_product!(scratch, view(x_data, :, primitive),
                view(y_data, :, primitive), view(z_data, :, primitive))
            column .+= Float64(orbital.coefficients[primitive]) .* scratch
        end
    end
    return (; ga = (; overlap = overlap_ga,))
end

function _non_nuclear_aa_axis_tables(inventory, axis::Int)
    parent = _GB_PARENT
    families = inventory.families[axis]
    overlap = Matrix{Float64}[]
    kinetic = Matrix{Float64}[]
    position = Matrix{Float64}[]
    x2 = Matrix{Float64}[]
    ids = zeros(Int, length(families), length(families))
    reversed = falses(length(families), length(families))
    for left in eachindex(families), right in left:length(families)
        left_family = families[left]
        right_family = families[right]
        push!(overlap, parent._cartesian_gaussian_axis_integral_table(
            left_family.exponents, left_family.centers, left_family.powers,
            left_family.prefactors, right_family.exponents, right_family.centers,
            right_family.powers, right_family.prefactors, :overlap))
        push!(kinetic, parent._cartesian_gaussian_axis_integral_table(
            left_family.exponents, left_family.centers, left_family.powers,
            left_family.prefactors, right_family.exponents, right_family.centers,
            right_family.powers, right_family.prefactors, :kinetic))
        push!(position, parent._cartesian_gaussian_axis_integral_table(
            left_family.exponents, left_family.centers, left_family.powers,
            left_family.prefactors, right_family.exponents, right_family.centers,
            right_family.powers, right_family.prefactors, :position))
        push!(x2, parent._cartesian_gaussian_axis_integral_table(
            left_family.exponents, left_family.centers, left_family.powers,
            left_family.prefactors, right_family.exponents, right_family.centers,
            right_family.powers, right_family.prefactors, :x2))
        id = length(overlap)
        ids[left, right] = id
        ids[right, left] = id
        reversed[right, left] = left != right
    end
    return (; overlap, kinetic, position, x2, ids, reversed)
end

function _non_nuclear_weighted_hadamard3(
    left_coefficients, x, x_reversed, y, y_reversed, z, z_reversed, right_coefficients)
    return _GB_PARENT._cartesian_weighted_hadamard3(left_coefficients,
        right_coefficients, x_reversed ? transpose(x) : x,
        y_reversed ? transpose(y) : y, z_reversed ? transpose(z) : z)
end

function _non_nuclear_aa_values(left, right, tables, ids, reversed)
    x, y, z = tables[1], tables[2], tables[3]
    overlap = _non_nuclear_weighted_hadamard3(left.coefficients,
        x.overlap[ids[1]], reversed[1], y.overlap[ids[2]], reversed[2],
        z.overlap[ids[3]], reversed[3], right.coefficients)
    kinetic =
        _non_nuclear_weighted_hadamard3(left.coefficients, x.kinetic[ids[1]],
            reversed[1], y.overlap[ids[2]], reversed[2], z.overlap[ids[3]],
            reversed[3], right.coefficients) +
        _non_nuclear_weighted_hadamard3(left.coefficients, x.overlap[ids[1]],
            reversed[1], y.kinetic[ids[2]], reversed[2], z.overlap[ids[3]],
            reversed[3], right.coefficients) +
        _non_nuclear_weighted_hadamard3(left.coefficients, x.overlap[ids[1]],
            reversed[1], y.overlap[ids[2]], reversed[2], z.kinetic[ids[3]],
            reversed[3], right.coefficients)
    px = _non_nuclear_weighted_hadamard3(left.coefficients, x.position[ids[1]],
        reversed[1], y.overlap[ids[2]], reversed[2], z.overlap[ids[3]],
        reversed[3], right.coefficients)
    py = _non_nuclear_weighted_hadamard3(left.coefficients, x.overlap[ids[1]],
        reversed[1], y.position[ids[2]], reversed[2], z.overlap[ids[3]],
        reversed[3], right.coefficients)
    pz = _non_nuclear_weighted_hadamard3(left.coefficients, x.overlap[ids[1]],
        reversed[1], y.overlap[ids[2]], reversed[2], z.position[ids[3]],
        reversed[3], right.coefficients)
    x2x = _non_nuclear_weighted_hadamard3(left.coefficients, x.x2[ids[1]],
        reversed[1], y.overlap[ids[2]], reversed[2], z.overlap[ids[3]],
        reversed[3], right.coefficients)
    x2y = _non_nuclear_weighted_hadamard3(left.coefficients, x.overlap[ids[1]],
        reversed[1], y.x2[ids[2]], reversed[2], z.overlap[ids[3]],
        reversed[3], right.coefficients)
    x2z = _non_nuclear_weighted_hadamard3(left.coefficients, x.overlap[ids[1]],
        reversed[1], y.overlap[ids[2]], reversed[2], z.x2[ids[3]],
        reversed[3], right.coefficients)
    return (; overlap, kinetic, px, py, pz, x2x, x2y, x2z)
end

function _non_nuclear_aa_block_matrices(supplement,
    inventory = _axis_family_inventory(supplement))
    norbital = length(supplement.orbitals)
    overlap = zeros(Float64, norbital, norbital)
    kinetic = zeros(Float64, norbital, norbital)
    px = zeros(Float64, norbital, norbital)
    py = zeros(Float64, norbital, norbital)
    pz = zeros(Float64, norbital, norbital)
    x2x = zeros(Float64, norbital, norbital)
    x2y = zeros(Float64, norbital, norbital)
    x2z = zeros(Float64, norbital, norbital)
    tables = [_non_nuclear_aa_axis_tables(inventory, axis) for axis in 1:3]

    for left_index in 1:norbital, right_index in left_index:norbital
        left = supplement.orbitals[left_index]
        right = supplement.orbitals[right_index]
        ids = ntuple(axis -> tables[axis].ids[inventory.maps[axis][left_index],
            inventory.maps[axis][right_index]], 3)
        reversed = ntuple(axis -> tables[axis].reversed[
            inventory.maps[axis][left_index], inventory.maps[axis][right_index]], 3)
        values = _non_nuclear_aa_values(left, right, tables, ids, reversed)
        if left_index != right_index
            reverse_ids = ntuple(axis -> tables[axis].ids[
                inventory.maps[axis][right_index],
                inventory.maps[axis][left_index]], 3)
            reverse_reversed = ntuple(axis -> tables[axis].reversed[
                inventory.maps[axis][right_index],
                inventory.maps[axis][left_index]], 3)
            reverse_values = _non_nuclear_aa_values(
                right, left, tables, reverse_ids, reverse_reversed)
            values = map((a, b) -> 0.5 * (a + b), values, reverse_values)
        end
        overlap[left_index, right_index] = values.overlap
        kinetic[left_index, right_index] = values.kinetic
        px[left_index, right_index] = values.px
        py[left_index, right_index] = values.py
        pz[left_index, right_index] = values.pz
        x2x[left_index, right_index] = values.x2x
        x2y[left_index, right_index] = values.x2y
        x2z[left_index, right_index] = values.x2z
        overlap[right_index, left_index] = values.overlap
        kinetic[right_index, left_index] = values.kinetic
        px[right_index, left_index] = values.px
        py[right_index, left_index] = values.py
        pz[right_index, left_index] = values.pz
        x2x[right_index, left_index] = values.x2x
        x2y[right_index, left_index] = values.x2y
        x2z[right_index, left_index] = values.x2z
    end

    return (; overlap, kinetic, position = (x = px, y = py, z = pz),
        x2 = (x = x2x, y = x2y, z = x2z))
end

function _non_nuclear_ga_axis_tables(proxy, inventory, axis::Int)
    parent = _GB_PARENT
    layer = _axis_layer(proxy, axis)
    inputs = _proxy_axis_factor_inputs(layer)
    overlap = Matrix{Float64}[]
    kinetic = Matrix{Float64}[]
    position = Matrix{Float64}[]
    x2 = Matrix{Float64}[]
    for family in inventory.families[axis]
        push!(overlap, parent._qwrg_left_contract_cross_matrix(layer,
            parent._cartesian_gaussian_axis_integral_table(inputs.exponents,
                inputs.centers, inputs.powers, inputs.prefactors,
                family.exponents, family.centers, family.powers,
                family.prefactors, :overlap)))
        push!(kinetic, parent._qwrg_left_contract_cross_matrix(layer,
            parent._cartesian_gaussian_axis_integral_table(inputs.exponents,
                inputs.centers, inputs.powers, inputs.prefactors,
                family.exponents, family.centers, family.powers,
                family.prefactors, :kinetic)))
        push!(position, parent._qwrg_left_contract_cross_matrix(layer,
            parent._cartesian_gaussian_axis_integral_table(inputs.exponents,
                inputs.centers, inputs.powers, inputs.prefactors,
                family.exponents, family.centers, family.powers,
                family.prefactors, :position)))
        push!(x2, parent._qwrg_left_contract_cross_matrix(layer,
            parent._cartesian_gaussian_axis_integral_table(inputs.exponents,
                inputs.centers, inputs.powers, inputs.prefactors,
                family.exponents, family.centers, family.powers,
                family.prefactors, :x2)))
    end
    return (; overlap, kinetic, position, x2)
end

function _non_nuclear_ga_block_matrices(proxy, supplement, inventory)
    ncart, norbital = proxy.ncart, length(supplement.orbitals)
    overlap = zeros(Float64, ncart, norbital)
    kinetic = zeros(Float64, ncart, norbital)
    px = zeros(Float64, ncart, norbital)
    py = zeros(Float64, ncart, norbital)
    pz = zeros(Float64, ncart, norbital)
    x2x = zeros(Float64, ncart, norbital)
    x2y = zeros(Float64, ncart, norbital)
    x2z = zeros(Float64, ncart, norbital)
    tables = [_non_nuclear_ga_axis_tables(proxy, inventory, axis) for axis in 1:3]
    fill_product! = getfield(_GB_PARENT, :_qwrg_fill_product_column!)
    scratch = zeros(Float64, ncart)

    for (orbital_index, orbital) in pairs(supplement.orbitals)
        ids = ntuple(axis -> inventory.maps[axis][orbital_index], 3)
        x, y, z = tables[1], tables[2], tables[3]
        for primitive in eachindex(orbital.coefficients)
            coefficient = Float64(orbital.coefficients[primitive])
            fill_product!(scratch, view(x.overlap[ids[1]], :, primitive),
                view(y.overlap[ids[2]], :, primitive),
                view(z.overlap[ids[3]], :, primitive))
            view(overlap, :, orbital_index) .+= coefficient .* scratch
            fill_product!(scratch, view(x.kinetic[ids[1]], :, primitive),
                view(y.overlap[ids[2]], :, primitive),
                view(z.overlap[ids[3]], :, primitive))
            view(kinetic, :, orbital_index) .+= coefficient .* scratch
            fill_product!(scratch, view(x.overlap[ids[1]], :, primitive),
                view(y.kinetic[ids[2]], :, primitive),
                view(z.overlap[ids[3]], :, primitive))
            view(kinetic, :, orbital_index) .+= coefficient .* scratch
            fill_product!(scratch, view(x.overlap[ids[1]], :, primitive),
                view(y.overlap[ids[2]], :, primitive),
                view(z.kinetic[ids[3]], :, primitive))
            view(kinetic, :, orbital_index) .+= coefficient .* scratch
            fill_product!(scratch, view(x.position[ids[1]], :, primitive),
                view(y.overlap[ids[2]], :, primitive),
                view(z.overlap[ids[3]], :, primitive))
            view(px, :, orbital_index) .+= coefficient .* scratch
            fill_product!(scratch, view(x.overlap[ids[1]], :, primitive),
                view(y.position[ids[2]], :, primitive),
                view(z.overlap[ids[3]], :, primitive))
            view(py, :, orbital_index) .+= coefficient .* scratch
            fill_product!(scratch, view(x.overlap[ids[1]], :, primitive),
                view(y.overlap[ids[2]], :, primitive),
                view(z.position[ids[3]], :, primitive))
            view(pz, :, orbital_index) .+= coefficient .* scratch
            fill_product!(scratch, view(x.x2[ids[1]], :, primitive),
                view(y.overlap[ids[2]], :, primitive),
                view(z.overlap[ids[3]], :, primitive))
            view(x2x, :, orbital_index) .+= coefficient .* scratch
            fill_product!(scratch, view(x.overlap[ids[1]], :, primitive),
                view(y.x2[ids[2]], :, primitive),
                view(z.overlap[ids[3]], :, primitive))
            view(x2y, :, orbital_index) .+= coefficient .* scratch
            fill_product!(scratch, view(x.overlap[ids[1]], :, primitive),
                view(y.overlap[ids[2]], :, primitive),
                view(z.x2[ids[3]], :, primitive))
            view(x2z, :, orbital_index) .+= coefficient .* scratch
        end
    end
    return (; overlap, kinetic, position = (x = px, y = py, z = pz),
        x2 = (x = x2x, y = x2y, z = x2z))
end

function gaussian_non_nuclear_raw_blocks(proxy, supplement, expansion)
    inventory = _axis_family_inventory(supplement)
    ga = _non_nuclear_ga_block_matrices(proxy, supplement, inventory)
    aa = _non_nuclear_aa_block_matrices(supplement, inventory)
    return (;
        ga = (;
            overlap = ga.overlap,
            kinetic = ga.kinetic,
            position = ga.position,
            x2 = ga.x2,
        ),
        aa = (;
            overlap = aa.overlap,
            kinetic = aa.kinetic,
            position = aa.position,
            x2 = aa.x2,
        ),
    )
end
