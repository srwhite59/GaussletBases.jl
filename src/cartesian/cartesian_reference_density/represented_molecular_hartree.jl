struct RepresentedMixedDensity{B,U,S,R,D}
    basis::B
    bundles::U
    supplement::S
    residual::R
    final_alpha::Matrix{Float64}
    final_beta::Matrix{Float64}
    occupations_alpha::Vector{Float64}
    occupations_beta::Vector{Float64}
    C_G_alpha::Matrix{Float64}
    C_A_alpha::Matrix{Float64}
    C_G_beta::Matrix{Float64}
    C_A_beta::Matrix{Float64}
    X::Matrix{Float64}
    S_AA::Matrix{Float64}
    fingerprint::String
    diagnostics::D
end
struct ContractedMixedDensityPairAction{D,P,S,F,Q}
    density::D
    proxy::P
    target_supplement::S
    axis_functions::F
    component_axis_indices::Matrix{Int}
    coefficients_alpha::Matrix{Float64}
    coefficients_beta::Matrix{Float64}
    proxy_fingerprint::String
    fingerprint::String
    diagnostics::Q
end
struct DirectRepresentedHartreePotential{A,E}
    action::A
    expansion::E
    expansion_fingerprint::String
end
struct RepresentedHartreeField{P,D}
    GG::Matrix{Float64}
    GA::Matrix{Float64}
    AA::Matrix{Float64}
    native::Matrix{Float64}
    potential::P
    no_half_E0::Float64
    hartree_energy::Float64
    fingerprint::String
    diagnostics::D
end
_represented_array_fingerprint(array) = bytes2hex(sha256(
    reinterpret(UInt8, vec(Array(array)))))
_represented_record_fingerprint(value) = bytes2hex(sha256(codeunits(repr(value))))
_represented_basis_fingerprint(basis) = _represented_record_fingerprint((
    basis.final_dimension, [(block.unit_key, block.support_indices, block.support_states,
        isnothing(block.coefficients) ? nothing : _represented_array_fingerprint(block.coefficients),
        block.column_range) for block in basis.blocks]))
_represented_supplement_fingerprint(supplement) = _represented_record_fingerprint([
    (orbital.label, orbital.angular_powers, orbital.center,
        orbital.exponents, orbital.coefficients, orbital.primitive_normalization)
    for orbital in supplement.orbitals])
_represented_density_fingerprint(basis, supplement, residual, arrays...) =
    _represented_record_fingerprint((_represented_basis_fingerprint(basis),
        _represented_supplement_fingerprint(supplement), residual.residual_labels,
        map(_represented_array_fingerprint, (residual.T_G, residual.T_A, arrays...))))

function _represented_proxy_fingerprint(proxy)
    raw = _GB_PARENT.CartesianGaussianRawBlocks
    return _represented_record_fingerprint(ntuple(3) do axis
        layer = getfield(raw, :_axis_layer)(proxy, axis)
        inputs = getfield(raw, :_proxy_axis_factor_inputs)(layer)
        map(_represented_array_fingerprint, (_GB_PARENT.stencil_matrix(layer),
            inputs.exponents, inputs.centers, inputs.powers, inputs.prefactors))
    end)
end

function _represented_spin_reconstruction(final, occupations, nG, residual, X, S_AA,
    label, tolerance)
    coefficients = Matrix{Float64}(final)
    occupation = Vector{Float64}(occupations)
    nR = residual.residual_dimension
    size(coefficients, 1) == nG + nR && size(coefficients, 2) == length(occupation) ||
        throw(DimensionMismatch("$label final coefficient/occupation dimensions are inconsistent"))
    all(isfinite, coefficients) && all(isfinite, occupation) && all(>=(0.0), occupation) ||
        throw(ArgumentError("$label coefficients and occupations must be finite and nonnegative"))
    C_R = view(coefficients, (nG + 1):(nG + nR), :)
    C_G = Matrix{Float64}(view(coefficients, 1:nG, :) + residual.T_G * C_R)
    C_A = Matrix{Float64}(residual.T_A * C_R)
    final_gram = transpose(coefficients) * coefficients
    mixed_gram = transpose(C_G) * C_G + transpose(C_G) * X * C_A +
        transpose(C_A) * transpose(X) * C_G + transpose(C_A) * S_AA * C_A
    identity_error = isempty(occupation) ? 0.0 : norm(final_gram - I, Inf)
    recovery_error = isempty(occupation) ? 0.0 : norm(mixed_gram - final_gram, Inf)
    max(identity_error, recovery_error) <= tolerance ||
        throw(ArgumentError("$label represented mixed-basis recovery failed"))
    final_charge = dot(occupation, diag(final_gram))
    mixed_charge = dot(occupation, diag(mixed_gram))
    abs(final_charge - mixed_charge) <= tolerance * max(1.0, abs(final_charge)) ||
        throw(ArgumentError("$label represented charge recovery failed"))
    return (; coefficients, occupation, C_G, C_A, final_charge, mixed_charge,
        identity_error, recovery_error)
end
function represented_mixed_density(basis, bundles, supplement, residual,
    final_alpha, occupations_alpha, final_beta, occupations_beta;
    proxy = nothing, tolerance::Real = 1.0e-10)
    tol = Float64(tolerance)
    isfinite(tol) && tol > 0.0 || throw(ArgumentError("tolerance must be finite and positive"))
    nG, nA = basis.final_dimension, length(supplement.orbitals)
    residual.base_dimension == nG && residual.candidate_count == nA ||
        throw(DimensionMismatch("represented density basis/residual dimensions disagree"))
    size(residual.T_G) == (nG, residual.residual_dimension) &&
        size(residual.T_A) == (nA, residual.residual_dimension) ||
        throw(DimensionMismatch("represented density residual transforms have wrong dimensions"))
    isnothing(residual.injected_G) || throw(ArgumentError(
        "represented density requires native [G,R] ordering without residual-sector injection"))
    realizer = _GB_PARENT.CartesianFinalBasisRealization
    X = if isnothing(proxy)
        Matrix{Float64}(getfield(realizer, :_terminal_residual_mixed_overlap)(
            basis, bundles, supplement))
    else
        donor = getfield(realizer, :_r3a_qw_supplement)(supplement)
        parent_overlap = getfield(_GB_PARENT.CartesianGaussianRawBlocks,
            :gaussian_non_nuclear_overlap_blocks)(proxy, donor).ga.overlap
        Matrix{Float64}(getfield(realizer, :_r3a_project_parent_ga)(basis, parent_overlap))
    end
    S_AA = Matrix{Float64}(getfield(_GB_PARENT, :_cartesian_supplement_cross_overlap)(
        supplement, supplement))
    size(X) == (nG, nA) && size(S_AA) == (nA, nA) ||
        throw(DimensionMismatch("represented density mixed metric dimensions disagree"))
    norm(S_AA - transpose(S_AA), Inf) <= tol && all(isfinite, X) && all(isfinite, S_AA) ||
        throw(ArgumentError("represented density mixed metric is invalid"))
    all(isfinite, residual.T_G) && all(isfinite, residual.T_A) ||
        throw(ArgumentError("represented density residual transform is not finite"))
    G_R = residual.T_G + X * residual.T_A
    R_R = transpose(residual.T_G) * G_R + transpose(residual.T_A) *
        (transpose(X) * residual.T_G + S_AA * residual.T_A)
    residual_cross_error = isempty(G_R) ? 0.0 : norm(G_R, Inf)
    residual_metric_error = isempty(R_R) ? 0.0 : norm(R_R - Matrix{Float64}(
        getfield(_GB_PARENT, :I), residual.residual_dimension, residual.residual_dimension), Inf)
    max(residual_cross_error, residual_metric_error) <= tol ||
        throw(ArgumentError("represented density residual transform fails the complete mixed metric"))
    alpha = _represented_spin_reconstruction(final_alpha, occupations_alpha,
        nG, residual, X, S_AA, "alpha", tol)
    beta = _represented_spin_reconstruction(final_beta, occupations_beta,
        nG, residual, X, S_AA, "beta", tol)
    fingerprints = (; basis = _represented_basis_fingerprint(basis),
        supplement = _represented_supplement_fingerprint(supplement),
        residual = _represented_record_fingerprint((map(_represented_array_fingerprint,
            (residual.T_G, residual.T_A)), residual.residual_labels)),
        final_order = _represented_record_fingerprint(map(_represented_array_fingerprint,
            (alpha.coefficients, beta.coefficients))),
        metric = _represented_record_fingerprint(map(_represented_array_fingerprint, (X, S_AA))))
    fingerprint = _represented_density_fingerprint(basis, supplement, residual,
        alpha.coefficients, beta.coefficients, alpha.occupation, beta.occupation,
        alpha.C_G, alpha.C_A, beta.C_G, beta.C_A, X, S_AA)
    diagnostics = (; fingerprints,
        alpha_charge = alpha.mixed_charge, beta_charge = beta.mixed_charge,
        total_charge = alpha.mixed_charge + beta.mixed_charge,
        residual_cross_error, residual_metric_error,
        alpha_gram_error = alpha.recovery_error, beta_gram_error = beta.recovery_error,
        alpha_identity_error = alpha.identity_error, beta_identity_error = beta.identity_error)
    return RepresentedMixedDensity(basis, bundles, supplement, residual,
        alpha.coefficients, beta.coefficients, alpha.occupation, beta.occupation,
        alpha.C_G, alpha.C_A, beta.C_G, beta.C_A, X, S_AA, fingerprint, diagnostics)
end
function _represented_parent_coefficients(basis, coefficients, parent_count)
    result = zeros(Float64, parent_count, size(coefficients, 2))
    for block in basis.blocks
        target = view(result, block.support_indices, :)
        source = view(coefficients, block.column_range, :)
        isnothing(block.coefficients) ? (target .+= source) : mul!(
            target, block.coefficients, source, 1.0, 1.0)
    end
    return result
end
_represented_axis_function(inputs, coefficients) = (; exponents = inputs.exponents,
    centers = inputs.centers, powers = inputs.powers, prefactors = inputs.prefactors,
    coefficients = Vector{Float64}(coefficients))

function contracted_mixed_density_pair_action(density::RepresentedMixedDensity; proxy = nothing)
    current_density_fingerprint = _represented_density_fingerprint(density.basis,
        density.supplement, density.residual, density.final_alpha, density.final_beta,
        density.occupations_alpha, density.occupations_beta, density.C_G_alpha,
        density.C_A_alpha, density.C_G_beta, density.C_A_beta, density.X, density.S_AA)
    current_density_fingerprint == density.fingerprint ||
        throw(ArgumentError("represented density changed after validation"))
    realizer = _GB_PARENT.CartesianFinalBasisRealization
    raw = _GB_PARENT.CartesianGaussianRawBlocks
    actual_proxy = isnothing(proxy) ? getfield(realizer, :_r3a_qw_proxy_layers)(density.bundles) : proxy
    donor = getfield(realizer, :_r3a_qw_supplement)(density.supplement)
    inventory = getfield(raw, :_axis_family_inventory)(donor)
    dims = getfield(_GB_PARENT, :_nested_axis_lengths)(density.bundles)
    actual_proxy.ncart == prod(dims) || throw(DimensionMismatch(
        "represented Hartree proxy dimension does not match parent axes"))
    axis_functions = ntuple(3) do axis
        layer = getfield(raw, :_axis_layer)(actual_proxy, axis)
        inputs = getfield(raw, :_proxy_axis_factor_inputs)(layer)
        stencil = Matrix{Float64}(_GB_PARENT.stencil_matrix(layer))
        [_represented_axis_function(inputs, view(stencil, :, column))
            for column in axes(stencil, 2)]
    end
    parent_alpha = _represented_parent_coefficients(
        density.basis, density.C_G_alpha, prod(dims))
    parent_beta = _represented_parent_coefficients(
        density.basis, density.C_G_beta, prod(dims))
    component_indices = NTuple{3,Int}[]
    alpha_rows = Vector{Float64}[]
    beta_rows = Vector{Float64}[]
    for parent_index in 1:prod(dims)
        alpha_row, beta_row = parent_alpha[parent_index, :], parent_beta[parent_index, :]
        any(!iszero, alpha_row) || any(!iszero, beta_row) || continue
        push!(component_indices, getfield(_GB_PARENT, :_cartesian_unflat_index)(parent_index, dims))
        push!(alpha_rows, Vector{Float64}(alpha_row)); push!(beta_rows, Vector{Float64}(beta_row))
    end
    shell_prefactor = getfield(_GB_PARENT.GaussianAnalyticIntegrals,
        :polynomial_gaussian_shell_prefactor)
    for (orbital_index, orbital) in pairs(density.supplement.orbitals),
            primitive in eachindex(orbital.coefficients)
        alpha_row = Float64(orbital.coefficients[primitive]) .* density.C_A_alpha[orbital_index, :]
        beta_row = Float64(orbital.coefficients[primitive]) .* density.C_A_beta[orbital_index, :]
        any(!iszero, alpha_row) || any(!iszero, beta_row) || continue
        ids = ntuple(3) do axis
            exponent = Float64(orbital.exponents[primitive])
            input = (; exponents = [exponent], centers = [Float64(orbital.center[axis])],
                powers = [orbital.angular_powers[axis]],
                prefactors = [shell_prefactor(exponent, orbital.angular_powers[axis])])
            push!(axis_functions[axis], _represented_axis_function(input, [1.0]))
            length(axis_functions[axis])
        end
        push!(component_indices, ids); push!(alpha_rows, alpha_row); push!(beta_rows, beta_row)
    end
    component_count = length(component_indices)
    component_count > 0 || throw(ArgumentError("represented density has no nonzero source components"))
    coefficients_alpha = reduce(vcat, transpose.(alpha_rows))
    coefficients_beta = reduce(vcat, transpose.(beta_rows))
    indices = reduce(vcat, (reshape(collect(index), 1, 3) for index in component_indices))
    nG, nA = density.basis.final_dimension, length(density.supplement.orbitals)
    parent_ga_shape = (prod(dims), nA)
    axis_table_entries = ntuple(3) do axis
        layer = getfield(raw, :_axis_layer)(actual_proxy, axis)
        nparent = size(_GB_PARENT.stencil_matrix(layer), 2)
        family_sizes = [length(family.exponents) for family in inventory.families[axis]]
        nfamily = sum(family_sizes)
        nparent^2 + nparent * nfamily + nfamily^2
    end
    primitive_scratch_entries = ntuple(3) do axis
        layer = getfield(raw, :_axis_layer)(actual_proxy, axis)
        nproxy = size(_GB_PARENT.stencil_matrix(layer), 1)
        nparent = size(_GB_PARENT.stencil_matrix(layer), 2)
        nfamily = maximum((length(family.exponents)
            for family in inventory.families[axis]); init = 0)
        max(nproxy^2 + nparent * nproxy, nfamily^2 + (nproxy + nparent) * nfamily,
            4 * nfamily^2)
    end
    output_bytes = 8 * (nG^2 + nG * nA + nA^2 + (nG + density.residual.residual_dimension)^2)
    terminal_buffer_entries = 3 * maximum(length(block.support_states) for block in density.basis.blocks)^2
    diagnostics = (; occupied_rank = size(coefficients_alpha, 2) + size(coefficients_beta, 2),
        source_component_count = component_count,
        source_block_pair_count = div(component_count * (component_count + 1), 2),
        retained_one_dimensional_cache_entries = 0,
        parent_axis_dimensions = dims, target_dimensions = (nG, nA), parent_ga_shape,
        axis_table_entries, primitive_scratch_entries, terminal_buffer_entries,
        output_bytes,
        peak_additional_workspace_bytes = 8 * (prod(parent_ga_shape) + prod(dims) +
            sum(axis_table_entries) + maximum(primitive_scratch_entries) + terminal_buffer_entries))
    proxy_fingerprint = _represented_proxy_fingerprint(actual_proxy)
    fingerprint = _represented_record_fingerprint((current_density_fingerprint,
        proxy_fingerprint, map(_represented_array_fingerprint,
            (indices, coefficients_alpha, coefficients_beta))))
    return ContractedMixedDensityPairAction(density, actual_proxy, donor, axis_functions,
        indices, coefficients_alpha, coefficients_beta, proxy_fingerprint, fingerprint, diagnostics)
end
function direct_represented_hartree_potential(density::RepresentedMixedDensity, expansion;
    proxy = nothing)
    expansion isa getfield(_GB_PARENT, :CoulombGaussianExpansion) ||
        throw(ArgumentError("represented Hartree requires a producer-owned Coulomb expansion"))
    all(isfinite, expansion.coefficients) && all(isfinite, expansion.exponents) &&
        all(>(0.0), expansion.exponents) ||
        throw(ArgumentError("represented Hartree expansion must be finite with positive exponents"))
    action = contracted_mixed_density_pair_action(density; proxy)
    expansion_fingerprint = getfield(_GB_PARENT.CartesianFinalBasisRealization,
        :_coulomb_expansion_fingerprint)(expansion)
    return DirectRepresentedHartreePotential(action, expansion, expansion_fingerprint)
end
function _represented_density_pair_weight(action, left, right)
    value = 0.0
    @inbounds for column in eachindex(action.density.occupations_alpha)
        value += action.density.occupations_alpha[column] *
            action.coefficients_alpha[left, column] * action.coefficients_alpha[right, column]
    end
    @inbounds for column in eachindex(action.density.occupations_beta)
        value += action.density.occupations_beta[column] *
            action.coefficients_beta[left, column] * action.coefficients_beta[right, column]
    end
    return left == right ? value : 2.0 * value
end
_represented_expectation(C, occupations, matrix) = sum((occupations[column] *
    dot(view(C, :, column), matrix, view(C, :, column)) for column in eachindex(occupations)); init = 0.0)
_represented_raw_expectation(C_G, C_A, occupations, GG, GA, AA) = sum((
    occupations[column] * (dot(view(C_G, :, column), GG, view(C_G, :, column)) +
        2.0 * dot(view(C_G, :, column), GA, view(C_A, :, column)) +
        dot(view(C_A, :, column), AA, view(C_A, :, column)))
    for column in eachindex(occupations)); init = 0.0)
function represented_hartree_field(potential::DirectRepresentedHartreePotential;
    tolerance::Real = 1.0e-10)
    action, density, expansion = potential.action, potential.action.density, potential.expansion
    density_fingerprint = _represented_density_fingerprint(density.basis, density.supplement,
        density.residual, density.final_alpha, density.final_beta, density.occupations_alpha,
        density.occupations_beta, density.C_G_alpha, density.C_A_alpha, density.C_G_beta,
        density.C_A_beta, density.X, density.S_AA)
    expansion_fingerprint = getfield(_GB_PARENT.CartesianFinalBasisRealization,
        :_coulomb_expansion_fingerprint)(expansion)
    proxy_fingerprint = _represented_proxy_fingerprint(action.proxy)
    density_fingerprint == density.fingerprint &&
        expansion_fingerprint == potential.expansion_fingerprint &&
        proxy_fingerprint == action.proxy_fingerprint ||
        throw(ArgumentError("represented Hartree inputs changed after potential construction"))
    tol = Float64(tolerance)
    isfinite(tol) && tol > 0.0 || throw(ArgumentError("tolerance must be finite and positive"))
    raw = _GB_PARENT.CartesianGaussianRawBlocks
    realizer = _GB_PARENT.CartesianFinalBasisRealization
    inventory = getfield(raw, :_axis_family_inventory)(action.target_supplement)
    nG, nA = density.basis.final_dimension, length(density.supplement.orbitals)
    GG = zeros(Float64, nG, nG)
    parent_GA = zeros(Float64, action.proxy.ncart, nA)
    AA = zeros(Float64, nA, nA)
    scratch = zeros(Float64, action.proxy.ncart)
    action_buffer = Ref(Matrix{Float64}(undef, 0, 0))
    tile_buffer = Ref(Matrix{Float64}(undef, 0, 0))
    block_buffer = Ref(Matrix{Float64}(undef, 0, 0))
    ncomponent = size(action.component_axis_indices, 1)
    for left in 1:ncomponent, right in left:ncomponent
        density_weight = _represented_density_pair_weight(action, left, right)
        density_weight == 0.0 && continue
        for expansion_index in eachindex(expansion.coefficients)
            scale = density_weight * Float64(expansion.coefficients[expansion_index])
            gamma = Float64(expansion.exponents[expansion_index])
            axis_tables = ntuple(3) do axis
                layer = getfield(raw, :_axis_layer)(action.proxy, axis)
                source_left = action.axis_functions[axis][action.component_axis_indices[left, axis]]
                source_right = action.axis_functions[axis][action.component_axis_indices[right, axis]]
                getfield(raw, :_mixed_hartree_pair_target_axis_tables)(layer,
                    inventory.families[axis], source_left, source_right, gamma)
            end
            getfield(realizer, :_assemble_terminal_product_operator!)(GG, density.basis,
                axis_tables[1].gg, axis_tables[2].gg, axis_tables[3].gg,
                action_buffer, tile_buffer, block_buffer; scale)
            getfield(raw, :_accumulate_mixed_hartree_pair_targets!)(parent_GA, AA,
                action.target_supplement, inventory, axis_tables, scale, scratch)
        end
    end
    GG = Matrix{Float64}(0.5 .* (GG .+ transpose(GG)))
    GA = Matrix{Float64}(getfield(realizer, :_r3a_project_parent_ga)(density.basis, parent_GA))
    AA = Matrix{Float64}(0.5 .* (AA .+ transpose(AA)))
    all(isfinite, GG) && all(isfinite, GA) && all(isfinite, AA) ||
        throw(ArgumentError("represented Hartree raw field is not finite"))
    native = getfield(_GB_PARENT.CartesianResidualGaussians,
        :transform_augmented_operator)(GG, GA, AA, density.residual)
    E0 = _represented_expectation(density.final_alpha, density.occupations_alpha, native) +
        _represented_expectation(density.final_beta, density.occupations_beta, native)
    raw_E0 = _represented_raw_expectation(density.C_G_alpha, density.C_A_alpha,
        density.occupations_alpha, GG, GA, AA) + _represented_raw_expectation(
        density.C_G_beta, density.C_A_beta, density.occupations_beta, GG, GA, AA)
    closure_error = abs(E0 - raw_E0)
    closure_error <= tol * max(1.0, abs(E0)) ||
        throw(ArgumentError("represented Hartree native/raw energy closure failed"))
    symmetry_error = norm(native - transpose(native), Inf)
    symmetry_error <= tol * max(1.0, norm(native, Inf)) ||
        throw(ArgumentError("represented Hartree native field is not symmetric"))
    fingerprint = _represented_record_fingerprint((action.fingerprint,
        potential.expansion_fingerprint,
        _represented_array_fingerprint(GG), _represented_array_fingerprint(GA),
        _represented_array_fingerprint(AA), _represented_array_fingerprint(native)))
    diagnostics = (; raw_native_energy_closure = closure_error,
        raw_symmetry_error = max(norm(GG - transpose(GG), Inf),
            norm(AA - transpose(AA), Inf)), native_symmetry_error = symmetry_error,
        resource = action.diagnostics)
    return RepresentedHartreeField(GG, GA, AA, native, potential, E0, 0.5 * E0,
        fingerprint, diagnostics)
end
