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

function gaussian_non_nuclear_raw_blocks(proxy, supplement, expansion)
    cross = getfield(_GB_PARENT, :_qwrg_cartesian_shell_cross_moment_blocks_3d)(
        (x = proxy.x, y = proxy.y, z = proxy.z),
        supplement,
        expansion,
        proxy.ncart;
        include_factor_terms = false,
    )
    self = getfield(_GB_PARENT, :_qwrg_cartesian_shell_self_moment_blocks_3d)(
        supplement,
        expansion;
        include_factor_terms = false,
    )
    return (;
        ga = (;
            overlap = cross.overlap_ga,
            kinetic = cross.kinetic_ga,
            position = (x = cross.position_x_ga, y = cross.position_y_ga,
                z = cross.position_z_ga),
            x2 = (x = cross.x2_x_ga, y = cross.x2_y_ga, z = cross.x2_z_ga),
        ),
        aa = (;
            overlap = _symmetrize_raw_block(self.overlap_aa),
            kinetic = _symmetrize_raw_block(self.kinetic_aa),
            position = (x = _symmetrize_raw_block(self.position_x_aa),
                y = _symmetrize_raw_block(self.position_y_aa),
                z = _symmetrize_raw_block(self.position_z_aa)),
            x2 = (x = _symmetrize_raw_block(self.x2_x_aa),
                y = _symmetrize_raw_block(self.x2_y_aa),
                z = _symmetrize_raw_block(self.x2_z_aa)),
        ),
    )
end
