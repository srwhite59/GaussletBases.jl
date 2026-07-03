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

function transform_protected_original_fixed_sector_operator(O_GG, O_GA, O_AA, components)
    Z, Gp, Ap = components.Z, components.G_perp, components.A_perp
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
    return symmetrize_operator(out)
end

function transform_protected_original_fixed_sector_one_body(kinetic, unit_nuclear_by_center,
    nuclear_charges, geometry)
    components = protected_original_fixed_sector_components(geometry)
    K = transform_protected_original_fixed_sector_operator(
        kinetic.GG, kinetic.GA, kinetic.AA, components)
    U = Matrix{Float64}[transform_protected_original_fixed_sector_operator(
        unit.GG, unit.GA, unit.AA, components) for unit in unit_nuclear_by_center]
    length(U) == length(nuclear_charges) ||
        throw(DimensionMismatch("protected fixed-sector nuclear charge count mismatch"))
    H1 = copy(K)
    for (charge, unit) in zip(nuclear_charges, U)
        H1 .+= Float64(charge) .* unit
    end
    return (; kinetic = K, nuclear_attraction_unit_by_center = U,
        one_body_hamiltonian = symmetrize_operator(H1),
        protected_count = components.protected_count,
        z_dimension = components.z_dimension,
        f_dimension = components.f_dimension)
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

function _hartree_anchor_diagnostics(exact, fapp_alpha, fapp_beta, delta_alpha,
    delta_beta, density_alpha, density_beta, exact_energy, app_energy, c0)
    alpha_trace = tr(density_alpha)
    beta_trace = tr(density_beta)
    new_energy = app_energy + _hartree_anchor_trace(density_alpha, delta_alpha) +
        _hartree_anchor_trace(density_beta, delta_beta) + c0
    alpha_anchor = fapp_alpha + delta_alpha - exact
    beta_anchor = fapp_beta + delta_beta - exact
    alpha_eigs = eigvals(Symmetric(delta_alpha))
    beta_eigs = eigvals(Symmetric(delta_beta))
    return (;
        dimension = size(exact, 1),
        trace_alpha = alpha_trace,
        trace_beta = beta_trace,
        exact_hartree_energy = exact_energy,
        app_interaction_energy = app_energy,
        anchored_interaction_energy = new_energy,
        energy_anchor_error = new_energy - exact_energy,
        f_anchor_alpha_error = norm(alpha_anchor, Inf),
        f_anchor_beta_error = norm(beta_anchor, Inf),
        exact_hartree_symmetry_error = norm(exact - transpose(exact), Inf),
        f_app_alpha_symmetry_error = norm(fapp_alpha - transpose(fapp_alpha), Inf),
        f_app_beta_symmetry_error = norm(fapp_beta - transpose(fapp_beta), Inf),
        delta_alpha_symmetry_error = norm(delta_alpha - transpose(delta_alpha), Inf),
        delta_beta_symmetry_error = norm(delta_beta - transpose(delta_beta), Inf),
        exact_hartree_finite = all(isfinite, exact),
        f_app_finite = all(isfinite, fapp_alpha) && all(isfinite, fapp_beta),
        delta_finite = all(isfinite, delta_alpha) && all(isfinite, delta_beta),
        delta_alpha_eig_min = minimum(alpha_eigs),
        delta_alpha_eig_max = maximum(alpha_eigs),
        delta_beta_eig_min = minimum(beta_eigs),
        delta_beta_eig_max = maximum(beta_eigs),
        delta_alpha_diag_min = minimum(diag(delta_alpha)),
        delta_alpha_diag_max = maximum(diag(delta_alpha)),
        delta_beta_diag_min = minimum(diag(delta_beta)),
        delta_beta_diag_max = maximum(diag(delta_beta)),
        delta_alpha_reference_expectation = _hartree_anchor_trace(
            density_alpha, delta_alpha),
        delta_beta_reference_expectation = _hartree_anchor_trace(
            density_beta, delta_beta),
        delta_alpha_trace_normalized_expectation =
            alpha_trace == 0.0 ? NaN :
            _hartree_anchor_trace(density_alpha, delta_alpha) / alpha_trace,
        delta_beta_trace_normalized_expectation =
            beta_trace == 0.0 ? NaN :
            _hartree_anchor_trace(density_beta, delta_beta) / beta_trace)
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
    app_energy = getfield(parent, :_cartesian_ida_approximate_interaction_energy)(
        ham, alpha, beta)
    fapp_alpha = getfield(parent, :_cartesian_ida_approximate_interaction_fock_alpha)(
        ham, alpha, beta)
    fapp_beta = getfield(parent, :_cartesian_ida_approximate_interaction_fock_beta)(
        ham, alpha, beta)
    delta_alpha = symmetrize_operator(exact - fapp_alpha)
    delta_beta = symmetrize_operator(exact - fapp_beta)
    exact_energy = isnothing(exact_hartree_energy) ?
        0.5 * _hartree_anchor_trace(alpha + beta, exact) :
        Float64(exact_hartree_energy)
    isfinite(exact_energy) ||
        throw(ArgumentError("exact_hartree_energy must be finite"))
    c0 = exact_energy - app_energy -
        _hartree_anchor_trace(alpha, delta_alpha) -
        _hartree_anchor_trace(beta, delta_beta)
    diagnostics = _hartree_anchor_diagnostics(exact, fapp_alpha, fapp_beta,
        delta_alpha, delta_beta, alpha, beta, exact_energy, app_energy, c0)
    return (; exact_hartree = exact,
        f_app_interaction_alpha = fapp_alpha,
        f_app_interaction_beta = fapp_beta,
        delta_F0_alpha = delta_alpha,
        delta_F0_beta = delta_beta,
        C0 = c0,
        exact_hartree_energy = exact_energy,
        app_interaction_energy = app_energy,
        diagnostics)
end
