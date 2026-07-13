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

function numerical_complete_reference_blocks_in_augmented_basis(
    residual, X, S_AA, occupied_blocks, occupations;
    recovery_atol::Real = 1.0e-10)
    length(occupied_blocks) == length(occupations) ||
        throw(DimensionMismatch("occupied block/occupation count mismatch"))
    tolerance = Float64(recovery_atol)
    isfinite(tolerance) && tolerance >= 0.0 ||
        throw(ArgumentError("recovery_atol must be finite and nonnegative"))
    residual.occupation_cutoff == 1.0e-10 &&
        residual.residual_injection_cutoff == 0.0 &&
        isnothing(residual.injected_G) &&
        isnothing(residual.compact_source_candidate_indices) ||
        throw(ArgumentError("reference representation requires the numerical-complete residual policy"))
    nG, nA = residual.base_dimension, residual.candidate_count
    size(X) == (nG, nA) ||
        throw(DimensionMismatch("numerical-complete mixed overlap dimension mismatch"))
    size(S_AA) == (nA, nA) ||
        throw(DimensionMismatch("numerical-complete supplement overlap dimension mismatch"))
    coefficient_blocks = Matrix{Float64}[]
    diagnostics = NamedTuple[]
    for (block_raw, occupations_raw) in zip(occupied_blocks, occupations)
        block = Matrix{Float64}(block_raw)
        occ = Vector{Float64}(occupations_raw)
        size(block, 1) == nA && size(block, 2) == length(occ) ||
            throw(DimensionMismatch("occupied packet block dimension mismatch"))
        all(isfinite, block) && all(isfinite, occ) && all(>=(0.0), occ) ||
            throw(ArgumentError("occupied packet block and occupations must be finite and nonnegative"))
        C_G = X * block
        C_R = transpose(residual.T_G) * C_G +
            transpose(residual.T_A) * S_AA * block
        C_M = vcat(C_G, C_R)
        singulars = Float64[svdvals(C_M)...]
        recovery_loss = maximum(abs.(1.0 .- singulars))
        gram_error = norm(transpose(C_M) * C_M - I, Inf)
        electron_count = dot(occ, vec(sum(abs2, C_M; dims = 1)))
        expected_electron_count = sum(occ)
        electron_trace_error = abs(electron_count - expected_electron_count)
        maximum((recovery_loss, gram_error, electron_trace_error)) <= tolerance ||
            throw(ArgumentError("numerical-complete occupied packet recovery failed"))
        push!(coefficient_blocks, C_M)
        push!(diagnostics, (; recovery_singular_values = singulars, recovery_loss,
            gram_error, electron_count, expected_electron_count,
            electron_trace_error))
    end
    return (; coefficient_blocks, diagnostics)
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

function _transform_protected_original_fixed_sector_operator_result(
    O_GG, O_GA, O_AA, components)
    Z, Gp, Ap = components.Z, components.G_perp, components.A_perp
    nG, nA = size(Gp, 1), size(Z, 1)
    size(O_GG) == (nG, nG) ||
        throw(DimensionMismatch("protected fixed-sector operator GG dimension mismatch"))
    size(O_GA) == (nG, nA) ||
        throw(DimensionMismatch("protected fixed-sector operator GA dimension mismatch"))
    size(O_AA) == (nA, nA) ||
        throw(DimensionMismatch("protected fixed-sector operator AA dimension mismatch"))
    all(isfinite, O_GG) && all(isfinite, O_GA) && all(isfinite, O_AA) ||
        throw(ArgumentError("protected fixed-sector operator raw blocks must be finite"))
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
    size(out) == (components.f_dimension, components.f_dimension) ||
        throw(DimensionMismatch("protected fixed-sector operator dimension mismatch"))
    all(isfinite, out) ||
        throw(ArgumentError("protected fixed-sector transformed operator must be finite"))
    operator = symmetrize_operator(out)
    all(isfinite, operator) ||
        throw(ArgumentError("protected fixed-sector symmetrized operator must be finite"))
    return (; operator, diagnostics = (;
        raw_gg_symmetry_error = norm(O_GG - transpose(O_GG), Inf),
        raw_aa_symmetry_error = norm(O_AA - transpose(O_AA), Inf),
        transformed_pre_cleanup_symmetry_error = norm(out - transpose(out), Inf),
        trace = tr(out)))
end

function transform_protected_original_fixed_sector_operator(O_GG, O_GA, O_AA, components)
    return _transform_protected_original_fixed_sector_operator_result(
        O_GG, O_GA, O_AA, components).operator
end

function transform_protected_original_fixed_sector_one_body(kinetic, unit_nuclear_by_center,
    nuclear_charges, geometry)
    components = protected_original_fixed_sector_components(geometry)
    length(unit_nuclear_by_center) == length(nuclear_charges) ||
        throw(DimensionMismatch("protected fixed-sector nuclear charge count mismatch"))
    charges = Float64.(nuclear_charges)
    all(isfinite, charges) ||
        throw(ArgumentError("protected fixed-sector nuclear charges must be finite"))
    kinetic_result = _transform_protected_original_fixed_sector_operator_result(
        kinetic.GG, kinetic.GA, kinetic.AA, components)
    unit_results = [_transform_protected_original_fixed_sector_operator_result(
        unit.GG, unit.GA, unit.AA, components) for unit in unit_nuclear_by_center]
    K = kinetic_result.operator
    U = Matrix{Float64}[result.operator for result in unit_results]
    expected_size = (components.f_dimension, components.f_dimension)
    size(K) == expected_size && all(unit -> size(unit) == expected_size, U) ||
        throw(DimensionMismatch("protected fixed-sector one-body matrix dimension mismatch"))
    all(isfinite, K) && all(unit -> all(isfinite, unit), U) ||
        throw(ArgumentError("protected fixed-sector one-body matrices must be finite"))
    H1 = copy(K)
    for (charge, unit) in zip(charges, U)
        H1 .+= charge .* unit
    end
    size(H1) == expected_size ||
        throw(DimensionMismatch("protected fixed-sector assembled H1 dimension mismatch"))
    all(isfinite, H1) ||
        throw(ArgumentError("protected fixed-sector assembled H1 must be finite"))
    H1_symmetry_error = norm(H1 - transpose(H1), Inf)
    return (; kinetic = K, nuclear_attraction_unit_by_center = U,
        one_body_hamiltonian = symmetrize_operator(H1),
        protected_count = components.protected_count,
        z_dimension = components.z_dimension,
        f_dimension = components.f_dimension,
        diagnostics = (;
            dimensions = (; base = size(components.G_perp, 1),
                supplement = size(components.Z, 1),
                fixed_sector = components.f_dimension,
                center_count = length(U)),
            kinetic = kinetic_result.diagnostics,
            unit_nuclear_by_center = [result.diagnostics for result in unit_results],
            one_body_hamiltonian = (;
                transformed_pre_cleanup_symmetry_error = H1_symmetry_error,
                trace = tr(H1)),
            geometry = (;
                f_s_f_identity_block_max = geometry.f_s_f_identity_block_max,
                z_m_qperp_max = geometry.z_m_qperp_max,
                qperp_identity_sample_max = geometry.qperp_identity_sample_max,
                protected_span_min_sv = geometry.protected_span_min_sv)))
end

function protected_localized_inherited_site_transform(geometry)
    components = protected_original_fixed_sector_components(geometry)
    C = Matrix{Float64}(geometry.B)
    Qp = Matrix{Float64}(components.Qp)
    A = vcat(transpose(C), transpose(Qp))
    Sbar = symmetrize_operator(transpose(A) * A)
    values, vectors = eigen(Symmetric(Sbar))
    minimum(values) > 0.0 ||
        throw(ArgumentError("protected-localized metric is not positive definite"))
    invsqrt = vectors * Diagonal(1.0 ./ sqrt.(values)) * transpose(vectors)
    W = A * invsqrt
    ML = hcat(C, Qp) * W
    off = copy(ML)
    for i in axes(off, 1)
        off[i, i] = 0.0
    end
    return (; W, C, Qp, Sbar_values = values,
        L_identity_error = norm(transpose(W) * W - I, Inf),
        M_L_diag_min = minimum(diag(ML)),
        M_L_diag_max = maximum(diag(ML)),
        M_L_diag_delta_max = maximum(abs.(diag(ML) .- 1.0)),
        M_L_offdiag_max = maximum(abs, off),
        M_L_fro_delta = norm(ML - I) / sqrt(length(ML)))
end

function protected_localized_inherited_site_hamiltonian(
    protected_fixed_one_body,
    inherited_site_vee,
    geometry,
)
    loc = protected_localized_inherited_site_transform(geometry)
    H1 = symmetrize_operator(transpose(loc.W) * protected_fixed_one_body * loc.W)
    Vee = symmetrize_operator(inherited_site_vee)
    size(H1) == size(Vee) ||
        throw(DimensionMismatch("protected-localized H1/Vee dimension mismatch"))
    all(isfinite, H1) && all(isfinite, Vee) ||
        throw(ArgumentError("protected-localized matrices must be finite"))
    return (; H1_L = H1, Vee_L = Vee, transform = loc,
        diagnostics = (;
            L_identity_error = loc.L_identity_error,
            M_L_diag_delta_max = loc.M_L_diag_delta_max,
            M_L_offdiag_max = loc.M_L_offdiag_max,
            M_L_fro_delta = loc.M_L_fro_delta,
            H1_L_symmetry_error = norm(H1 - transpose(H1), Inf),
            Vee_L_symmetry_error = norm(Vee - transpose(Vee), Inf)))
end
function protected_original_reference_blocks_in_local_basis(
    G_L, A_L, X, S_AA, loc, occupied_blocks; recovery_atol::Real = 1.0e-10)
    coefficient_blocks = Matrix{Float64}[]
    diagnostics = NamedTuple[]
    for block_raw in occupied_blocks
        block = Matrix{Float64}(block_raw)
        C_L = transpose(G_L) * X * block + transpose(A_L) * S_AA * block
        C_F = loc.W * C_L
        f_singulars, l_singulars = Float64[svdvals(C_F)...], Float64[svdvals(C_L)...]
        f_loss, l_loss = maximum(abs.(1.0 .- f_singulars)),
            maximum(abs.(1.0 .- l_singulars))
        orthogonality_error = norm(transpose(C_L) * C_L - I, Inf)
        maximum((f_loss, l_loss, orthogonality_error)) <= recovery_atol ||
            throw(ArgumentError("protected occupied block is not recovered in F/L"))
        push!(coefficient_blocks, C_L)
        push!(diagnostics, (; f_recovery_singular_values = f_singulars,
            f_recovery_loss = f_loss, l_recovery_singular_values = l_singulars,
            l_recovery_loss = l_loss, l_orthogonality_error = orthogonality_error))
    end
    return (; coefficient_blocks, diagnostics)
end

function _protected_localized_axis_expectations(operator, ML, name::AbstractString)
    dense = symmetrize_operator(operator)
    size(dense) == (size(ML, 1), size(ML, 1)) ||
        throw(DimensionMismatch("$(name) dimension mismatch"))
    all(isfinite, dense) || throw(ArgumentError("$(name) must be finite"))
    return vec(sum(ML .* (dense * ML), dims = 1))
end

function _protected_localized_native_sector_rows(sector_counts, dimension::Int)
    base = Int(sector_counts.base)
    compact_R = Int(sector_counts.compact_R)
    base + compact_R == dimension ||
        throw(DimensionMismatch("protected-localized sector count mismatch"))
    labels = Vector{String}(undef, dimension)
    indices = Vector{Int}(undef, dimension)
    for row in 1:dimension
        if row <= base
            labels[row] = "base_G"
            indices[row] = row
        else
            labels[row] = "compact_R"
            indices[row] = row - base
        end
    end
    return labels, indices
end

function protected_localized_row_locality(
    loc,
    position;
    sector_counts,
    x2 = nothing,
)
    ML = hcat(loc.C, loc.Qp) * loc.W
    center_x = _protected_localized_axis_expectations(position.x, ML, "position.x")
    center_y = _protected_localized_axis_expectations(position.y, ML, "position.y")
    center_z = _protected_localized_axis_expectations(position.z, ML, "position.z")
    z_order_to_native = sort(collect(eachindex(center_z)), by = row -> (center_z[row], row))
    native_to_z_order = zeros(Int, length(center_z))
    for (z_index, native_index) in pairs(z_order_to_native)
        native_to_z_order[native_index] = z_index
    end
    sector_label, native_sector_index =
        _protected_localized_native_sector_rows(sector_counts, length(center_z))
    row_locality = (; center_x, center_y, center_z,
        native_to_z_order, z_order_to_native, sector_label, native_sector_index)
    isnothing(x2) && return row_locality
    second_x = _protected_localized_axis_expectations(x2.x, ML, "x2.x")
    second_y = _protected_localized_axis_expectations(x2.y, ML, "x2.y")
    second_z = _protected_localized_axis_expectations(x2.z, ML, "x2.z")
    spreads = (; spread_x = sqrt.(max.(0.0, second_x .- center_x .^ 2)),
        spread_y = sqrt.(max.(0.0, second_y .- center_y .^ 2)),
        spread_z = sqrt.(max.(0.0, second_z .- center_z .^ 2)))
    return merge(row_locality, spreads)
end

function _exact_hartree_raw_block_matrices(raw_blocks)
    GG = Matrix{Float64}(raw_blocks.GG)
    GA = Matrix{Float64}(raw_blocks.GA)
    AA = Matrix{Float64}(raw_blocks.AA)
    return (; GG, GA, AA)
end

function _validate_exact_hartree_raw_blocks(raw, components)
    nG, nA = size(components.G_perp, 1), size(components.Z, 1)
    size(raw.GG) == (nG, nG) ||
        throw(DimensionMismatch("exact Hartree GG dimension mismatch"))
    size(raw.GA) == (nG, nA) ||
        throw(DimensionMismatch("exact Hartree GA dimension mismatch"))
    size(raw.AA) == (nA, nA) ||
        throw(DimensionMismatch("exact Hartree AA dimension mismatch"))
    all(isfinite, raw.GG) && all(isfinite, raw.GA) && all(isfinite, raw.AA) ||
        throw(ArgumentError("exact Hartree raw blocks must be finite"))
    return nothing
end

function _protected_original_hartree_transform_diagnostics(raw, hartree, geometry, components)
    b_singulars = hasproperty(geometry, :B) ? svdvals(Matrix{Float64}(geometry.B)) : Float64[]
    b_sorted = sort(b_singulars)
    b_median = isempty(b_sorted) ? NaN : b_sorted[cld(length(b_sorted), 2)]
    return (;
        base_dimension = size(raw.GG, 1),
        supplement_dimension = size(raw.AA, 1),
        z_dimension = components.z_dimension,
        qperp_dimension = size(components.G_perp, 2),
        f_dimension = size(hartree, 1),
        protected_count = components.protected_count,
        raw_gg_symmetry_error = norm(raw.GG - transpose(raw.GG), Inf),
        raw_aa_symmetry_error = norm(raw.AA - transpose(raw.AA), Inf),
        raw_blocks_finite = all(isfinite, raw.GG) && all(isfinite, raw.GA) &&
            all(isfinite, raw.AA),
        hartree_trace = tr(hartree),
        hartree_symmetry_error = norm(hartree - transpose(hartree), Inf),
        hartree_finite = all(isfinite, hartree),
        qperp_identity_error = norm(transpose(components.Qp) * components.Qp -
            Matrix{Float64}(I, size(components.Qp, 2), size(components.Qp, 2)), Inf),
        b_min = isempty(b_singulars) ? NaN : minimum(b_singulars),
        b_median,
        b_max = isempty(b_singulars) ? NaN : maximum(b_singulars),
        b_below_0p999 = count(<(0.999), b_singulars),
        b_below_0p99 = count(<(0.99), b_singulars),
        b_below_0p95 = count(<(0.95), b_singulars),
        b_below_0p9 = count(<(0.9), b_singulars))
end

function transform_protected_original_fixed_sector_exact_hartree(raw_blocks, geometry)
    components = protected_original_fixed_sector_components(geometry)
    raw = _exact_hartree_raw_block_matrices(raw_blocks)
    _validate_exact_hartree_raw_blocks(raw, components)
    hartree = transform_protected_original_fixed_sector_operator(
        raw.GG, raw.GA, raw.AA, components)
    diagnostics = _protected_original_hartree_transform_diagnostics(
        raw, hartree, geometry, components)
    return (; hartree, exact_hartree = hartree, diagnostics,
        protected_count = components.protected_count,
        z_dimension = components.z_dimension,
        f_dimension = components.f_dimension)
end

function transform_protected_original_localized_exact_hartree(raw_blocks, geometry, loc)
    fixed = transform_protected_original_fixed_sector_exact_hartree(raw_blocks, geometry)
    localized = symmetrize_operator(transpose(loc.W) * fixed.hartree * loc.W)
    return (; J0_L = localized,
        diagnostics = merge(fixed.diagnostics, (; localized_symmetry_error =
            norm(localized - transpose(localized), Inf),
            localized_finite = all(isfinite, localized))))
end

function atomic_reference_protected_original_fixed_sector_exact_hartree(
    basis, bundles, proxy, supplement, reference_supplement, reference_density, geometry;
    expansion = getfield(parentmodule(@__MODULE__), :coulomb_gaussian_expansion)(doacc = false))
    raw_owner = getfield(parentmodule(@__MODULE__), :CartesianGaussianRawBlocks)
    gg = raw_owner.atomic_reference_hartree_gg_block(
        basis, bundles, reference_supplement, reference_density; expansion)
    ga_aa = raw_owner.atomic_reference_hartree_ga_aa_blocks(
        proxy, supplement, reference_supplement, reference_density; expansion)
    final_ga = getfield(getfield(parentmodule(@__MODULE__), :CartesianFinalBasisRealization),
        :_r3a_project_parent_ga)(basis, ga_aa.GA)
    result = transform_protected_original_fixed_sector_exact_hartree(
        (; GG = gg.GG, GA = final_ga, AA = ga_aa.AA), geometry)
    return merge(result, (; raw_parent_ga_dimension = size(ga_aa.GA),
        raw_diagnostics = (; gg = gg.diagnostics, ga_aa = ga_aa.diagnostics)))
end

function _hartree_anchor_density(density, dimension::Int, name::AbstractString)
    matrix = Matrix{Float64}(density)
    size(matrix) == (dimension, dimension) ||
        throw(DimensionMismatch("$(name) must have size $((dimension, dimension))"))
    all(isfinite, matrix) ||
        throw(ArgumentError("$(name) contains non-finite entries"))
    norm(matrix - transpose(matrix), Inf) <= 1.0e-8 ||
        throw(ArgumentError("$(name) must be symmetric"))
    return 0.5 .* (matrix .+ transpose(matrix))
end

function _hartree_anchor_matrix(matrix, dimension::Int, name::AbstractString)
    dense = _hartree_anchor_density(matrix, dimension, name)
    norm(dense - transpose(dense), Inf) <= 1.0e-8 ||
        throw(ArgumentError("$(name) must be symmetric"))
    return dense
end

function _hartree_anchor_trace(density, matrix)
    return Float64(sum(density .* matrix))
end

function _hartree_anchor_diagnostics(exact, fapp_direct, delta, density_alpha,
    density_beta, exact_energy, app_direct_energy, c0)
    alpha_trace = tr(density_alpha)
    beta_trace = tr(density_beta)
    total_density = density_alpha + density_beta
    new_energy = app_direct_energy + _hartree_anchor_trace(total_density, delta) + c0
    direct_anchor = fapp_direct + delta - exact
    eigs = eigvals(Symmetric(delta))
    return (;
        dimension = size(exact, 1),
        trace_alpha = alpha_trace,
        trace_beta = beta_trace,
        trace_total = alpha_trace + beta_trace,
        exact_hartree_energy = exact_energy,
        app_direct_energy = app_direct_energy,
        anchored_direct_energy = new_energy,
        energy_anchor_error = new_energy - exact_energy,
        f_anchor_direct_error = norm(direct_anchor, Inf),
        exact_hartree_symmetry_error = norm(exact - transpose(exact), Inf),
        f_app_direct_symmetry_error = norm(fapp_direct - transpose(fapp_direct), Inf),
        delta_symmetry_error = norm(delta - transpose(delta), Inf),
        exact_hartree_finite = all(isfinite, exact),
        f_app_finite = all(isfinite, fapp_direct),
        delta_finite = all(isfinite, delta),
        delta_eig_min = minimum(eigs),
        delta_eig_max = maximum(eigs),
        delta_diag_min = minimum(diag(delta)),
        delta_diag_max = maximum(diag(delta)),
        delta_reference_expectation = _hartree_anchor_trace(total_density, delta),
        delta_trace_normalized_expectation =
            alpha_trace + beta_trace == 0.0 ? NaN :
            _hartree_anchor_trace(total_density, delta) / (alpha_trace + beta_trace))
end

function hartree_reference_correction_anchor(
    ham,
    exact_hartree,
    density_alpha,
    density_beta;
    exact_hartree_energy = nothing,
)
    dimension = size(ham.kinetic, 1)
    size(ham.electron_electron_ida) == (dimension, dimension) ||
        throw(DimensionMismatch("Hamiltonian interaction dimension mismatch"))
    exact = _hartree_anchor_matrix(exact_hartree, dimension, "exact_hartree")
    alpha = _hartree_anchor_density(density_alpha, dimension, "density_alpha")
    beta = _hartree_anchor_density(density_beta, dimension, "density_beta")
    parent = parentmodule(@__MODULE__)
    app_direct_energy = getfield(parent, :_cartesian_ida_approximate_direct_interaction_energy)(
        ham, alpha, beta)
    fapp_direct = getfield(parent, :_cartesian_ida_approximate_direct_interaction_fock)(
        ham, alpha, beta)
    delta = symmetrize_operator(exact - fapp_direct)
    exact_energy = isnothing(exact_hartree_energy) ?
        0.5 * _hartree_anchor_trace(alpha + beta, exact) :
        Float64(exact_hartree_energy)
    isfinite(exact_energy) ||
        throw(ArgumentError("exact_hartree_energy must be finite"))
    c0 = exact_energy - app_direct_energy - _hartree_anchor_trace(alpha + beta, delta)
    diagnostics = _hartree_anchor_diagnostics(exact, fapp_direct, delta,
        alpha, beta, exact_energy, app_direct_energy, c0)
    return (; exact_hartree = exact,
        f_app_direct = fapp_direct,
        delta_J0 = delta,
        C0_J = c0,
        exact_hartree_energy = exact_energy,
        app_direct_energy = app_direct_energy,
        diagnostics)
end
