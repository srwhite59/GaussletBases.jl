using GaussletBases
using JLD2
using LinearAlgebra
using Test

const NUCLEI = NTuple{3,Float64}[(0.0, 0.0, -2.0), (0.0, 0.0, 2.0)]
const CRG = GaussletBases.CartesianResidualGaussians
const EXPECTED_LABELS = [
    "a_s1", "a_s2", "a_s3", "a_px1", "a_py1", "a_pz1", "a_px2", "a_py2", "a_pz2",
    "b_s1", "b_s2", "b_s3", "b_px1", "b_py1", "b_pz1", "b_px2", "b_py2", "b_pz2",
]
const R3B_OWNER_LOCAL_SELF_COULOMB = 0.4574265214362075

function terminal_h2()
    route_inputs = (;
        route_family = :pqs_source_box,
        route_kind = :bond_aligned_diatomic_independent_pqs_source_box_core_shell,
        route_shape = (:atom_contact_core, :shared_shell_1, :shared_shell_2),
        product_body_rule = :centered_single_z_slab,
        pqs_retained_rule = :boundary_comx_product_mode_selection,
        product_retained_rule = :product_doside_retained_unit,
        terms = (:overlap,),
        pair_factor_normalization = :density_normalized,
        white_lindsey_route_shape = (:standard_cartesian_units, :low_order_comx_coarsening),
        white_lindsey_mapping_rule = :standard_unit_backbone_mapping_family,
        white_lindsey_nesting_rule = :unit_box_low_order_comx_coarsening,
        white_lindsey_retained_rule = :low_order_unit_comx_retained_basis,
        white_lindsey_operator_rule = :low_order_unit_operator_blocks,
        supplement_policy = nothing,
        run_final_basis = true,
        run_h1 = false,
        run_h1_j = false,
    )
    system_inputs = (;
        atom_symbols = ("H", "H"),
        nuclear_charges = (1, 1),
        atom_locations = Tuple(NUCLEI),
        nup = 1,
        ndn = 1,
        bond_axis = :z,
        bond_length = 4.0,
        radius = 4.0,
        parent_axis_counts = nothing,
        map_backend = :pgdg_localized_experimental,
    )
    spacing_inputs = (;
        q = 5,
        n_s = 5,
        reference_spacing = 1.0,
        tail_spacing = 10.0,
        q_to_core_spacing_rule = :standard_pqs_ns_equals_q,
        core_spacing = 0.5,
        xmax_parallel = 6.0,
        xmax_transverse = 4.0,
    )
    parent_inputs = (;
        parent_axis_bundle_backend = :pgdg_localized_experimental,
        parent_axis_family = :G10,
        parent_mapping_rule = :identity_mapping,
        parent_mapping_Z = nothing,
        parent_mapping_d = nothing,
        parent_mapping_tail_spacing = 10.0,
    )
    system = GaussletBases.cartesian_system(system_inputs)
    recipe = GaussletBases.cartesian_recipe(route_inputs)
    parent = GaussletBases.cartesian_parent(system, spacing_inputs, parent_inputs, recipe)
    shells = GaussletBases.cartesian_shells(parent, spacing_inputs, recipe)
    units = GaussletBases.cartesian_units(parent, shells, recipe)
    transforms = GaussletBases.cartesian_transforms(units, recipe)
    return parent, transforms.terminal_basis_realization
end

function product_matrix(C, basis, ax, ay, az)
    matrix = zeros(Float64, basis.final_dimension, basis.final_dimension)
    C.assemble_terminal_product_operator!(matrix, basis, ax, ay, az)
    return matrix
end

function centered_factor_terms(axis, expansion, center)
    center == axis.center && return axis.gaussian_factor_terms
    return GaussletBases.mapped_ordinary_one_body_operators(
        axis.basis; exponents = expansion.exponents, center, backend = axis.backend).gaussian_factors
end

function base_blocks(C, parent, basis)
    pgdg = Tuple(GaussletBases._nested_axis_pgdg(parent.parent_axis_bundle_object, axis)
        for axis in (:x, :y, :z))
    S = Tuple(axis.overlap for axis in pgdg)
    K =
        product_matrix(C, basis, pgdg[1].kinetic, S[2], S[3]) +
        product_matrix(C, basis, S[1], pgdg[2].kinetic, S[3]) +
        product_matrix(C, basis, S[1], S[2], pgdg[3].kinetic)
    expansion = GaussletBases.coulomb_gaussian_expansion(doacc = false)
    U = Matrix{Float64}[]
    for location in NUCLEI
        matrix = zeros(Float64, basis.final_dimension, basis.final_dimension)
        factors = ntuple(axis -> centered_factor_terms(pgdg[axis], expansion, location[axis]), 3)
        C._accumulate_terminal_gaussian_sum!(
            matrix, basis, expansion.coefficients, factors[1], factors[2], factors[3])
        push!(U, matrix)
    end
    position = (;
        x = product_matrix(C, basis, pgdg[1].position, S[2], S[3]),
        y = product_matrix(C, basis, S[1], pgdg[2].position, S[3]),
        z = product_matrix(C, basis, S[1], S[2], pgdg[3].position),
    )
    x2 = (;
        x = product_matrix(C, basis, pgdg[1].x2, S[2], S[3]),
        y = product_matrix(C, basis, S[1], pgdg[2].x2, S[3]),
        z = product_matrix(C, basis, S[1], S[2], pgdg[3].x2),
    )
    return (; K, U, position, x2)
end

function symmetry_error(matrix)
    @test all(isfinite, matrix)
    return norm(matrix - transpose(matrix), Inf)
end

function self_coulomb(V, orbital)
    density = orbital * transpose(orbital)
    rho = 0.5 .* (density .+ transpose(density))
    occupations = vec(diag(rho))
    sym = 0.5 .* (V .+ transpose(V))
    return 2.0 * dot(occupations, sym * occupations) -
           dot(vec(rho), vec(sym .* rho))
end

function r3b_support_values(block, pair_terms, coefficients)
    nR = size(first(pair_terms.ga[1]), 2)
    values = zeros(Float64, length(block.support_states), nR)
    @inbounds for residual in 1:nR, (row, state) in pairs(block.support_states)
        ix, iy, iz = state
        value = 0.0
        for term in eachindex(coefficients)
            value += Float64(coefficients[term]) *
                pair_terms.ga[1][term][ix, residual] *
                pair_terms.ga[2][term][iy, residual] *
                pair_terms.ga[3][term][iz, residual]
        end
        values[row, residual] = value
    end
    return values
end

function independent_weight_aware_vgm(C, basis, bundles, pair_terms, coefficients)
    nR = size(first(pair_terms.ga[1]), 2)
    V = zeros(Float64, basis.final_dimension, nR)
    for block in basis.blocks
        values = r3b_support_values(block, pair_terms, coefficients)
        if isnothing(block.coefficients)
            V[block.column_range, :] .= values
        else
            support_weights = C._support_weights(block.support_states, bundles)
            final_weights = vec(transpose(block.coefficients) * support_weights)
            @test all(weight -> isfinite(weight) && weight > 1.0e-12, final_weights)
            C_density = block.coefficients .* reshape(support_weights, :, 1) ./
                        reshape(final_weights, 1, :)
            V[block.column_range, :] .= transpose(C_density) * values
        end
    end
    return V
end

function block_vgm_errors(basis, actual, expected)
    direct = 0.0
    pqs = 0.0
    for block in basis.blocks
        error = norm(actual[block.column_range, :] - expected[block.column_range, :], Inf)
        isnothing(block.coefficients) ? (direct = max(direct, error)) :
            (pqs = max(pqs, error))
    end
    return direct, pqs
end

const FACADE_SYSTEM = (;
    atom_symbols = ["H", "H"],
    nuclear_charges = [1.0, 1.0],
    atom_locations = NUCLEI,
    nup = 1,
    ndn = 1,
)
const FACADE_BASIS = (;
    q = 5,
    core_spacing = 0.5,
    xmax_parallel = 6.0,
    xmax_transverse = 4.0,
)
const FACADE_SUPPLEMENT = (;
    basis_by_center = ["cc-pVTZ", "cc-pVTZ"],
    lmax = 1,
    uncontracted = false,
    width_filtering = nothing,
)

elapsed = @elapsed @testset "R3-A H2 augmented one-body and moments" begin
    raw_supplement = legacy_bond_aligned_diatomic_gaussian_supplement(
        "H", "cc-pVTZ", NUCLEI; lmax = 1, uncontracted = false, max_width = nothing)
    supplement = basis_representation(raw_supplement)
    parent, basis = terminal_h2()
    C = GaussletBases.CartesianFinalBasisRealization
    residual = C.pqs_terminal_residual_gto_augmentation(
        basis, parent.parent_axis_bundle_object, supplement, NUCLEI)
    operators = C.pqs_terminal_residual_gto_augmented_operators(
        basis, parent.parent_axis_bundle_object, parent.parent_basis_object,
        supplement, residual, NUCLEI, [1.0, 1.0])

    X = C._terminal_residual_mixed_overlap(basis, parent.parent_axis_bundle_object, supplement)
    S_AA = GaussletBases._cartesian_supplement_cross_overlap(supplement, supplement)
    RSR = transpose(residual.T_G) * residual.T_G +
          transpose(residual.T_G) * X * residual.T_A +
          transpose(residual.T_A) * transpose(X) * residual.T_G +
          transpose(residual.T_A) * S_AA * residual.T_A

    @test residual.candidate_labels == EXPECTED_LABELS
    @test Dict(owner => count(==(owner), residual.candidate_owner_indices)
        for owner in unique(residual.candidate_owner_indices)) == Dict(1 => 9, 2 => 9)
    @test residual.residual_source_owner_indices == vcat(fill(1, 9), fill(2, 9))
    @test residual.owner_retained_counts == [9, 9]
    @test residual.occupation_cutoff == 1.0e-6
    @test residual.base_dimension == 471
    @test residual.residual_dimension == 18
    @test size(operators.kinetic) == (489, 489)
    @test minimum(residual.residual_occupations) ≈ 5.809976169475061e-4 atol = 1.0e-14
    @test maximum(residual.residual_occupations) ≈ 1.3455087681603393e-2 atol = 1.0e-14
    @test norm(residual.T_G + X * residual.T_A, Inf) <= 1.0e-10
    @test norm(RSR - I, Inf) <= 1.0e-10

    for matrix in (
            operators.kinetic,
            operators.nuclear_attraction_unit_by_center...,
            operators.position.x, operators.position.y, operators.position.z,
            operators.x2.x, operators.x2.y, operators.x2.z)
        @test symmetry_error(matrix) <= 1.0e-9
    end

    base = base_blocks(C, parent, basis)
    nG = residual.base_dimension
    @test norm(operators.kinetic[1:nG, 1:nG] - base.K, Inf) <= 1.0e-10
    @test norm(operators.nuclear_attraction_unit_by_center[1][1:nG, 1:nG] - base.U[1], Inf) <= 1.0e-10
    @test norm(operators.nuclear_attraction_unit_by_center[2][1:nG, 1:nG] - base.U[2], Inf) <= 1.0e-10
    for axis in (:x, :y, :z)
        @test norm(operators.position[axis][1:nG, 1:nG] - base.position[axis], Inf) <= 1.0e-10
        @test norm(operators.x2[axis][1:nG, 1:nG] - base.x2[axis], Inf) <= 1.0e-10
    end

    E_base = minimum(eigvals(Symmetric(base.K + base.U[1] + base.U[2])))
    H_aug = operators.kinetic + operators.nuclear_attraction_unit_by_center[1] +
            operators.nuclear_attraction_unit_by_center[2]
    E_aug = minimum(eigvals(Symmetric(H_aug)))
    @test E_base ≈ -0.7946037173365925 atol = 1.0e-10
    @test E_aug ≈ -0.7959028345077851 atol = 1.0e-10
    @test E_aug <= E_base + 1.0e-10

    base_ham = GaussletBases._cartesian_base_ida_hamiltonian(
        basis, parent.parent_axis_bundle_object, NUCLEI, [1.0, 1.0], 1, 1)
    ham = C.pqs_terminal_residual_gto_augmented_hamiltonian(
        base_ham, basis, parent.parent_axis_bundle_object, residual, operators)
    @test ham isa CartesianIDAHamiltonian{Float64}
    @test size(ham.electron_electron_ida) == (489, 489)
    @test symmetry_error(ham.electron_electron_ida) <= 1.0e-10
    @test norm(ham.electron_electron_ida[1:nG, 1:nG] -
               base_ham.electron_electron_ida, Inf) <= 1.0e-12

    expansion = GaussletBases.coulomb_gaussian_expansion(doacc = false)
    centers, widths = CRG.moment_matched_gaussians(operators, residual)
    pair_terms = CRG._mwg_axis_pairs(parent.parent_axis_bundle_object, expansion,
        centers, widths)
    V_GM_expected = independent_weight_aware_vgm(
        C, basis, parent.parent_axis_bundle_object, pair_terms, expansion.coefficients)
    V_GM_actual = ham.electron_electron_ida[1:nG, (nG + 1):end]
    direct_vgm_error, pqs_vgm_error = block_vgm_errors(basis, V_GM_actual, V_GM_expected)
    @test direct_vgm_error <= 1.0e-12
    @test pqs_vgm_error <= 1.0e-12

    H_r3b = one_body_hamiltonian(ham)
    eig = eigen(Symmetric(H_r3b))
    orbital = eig.vectors[:, argmin(eig.values)]
    r3b_self_coulomb = self_coulomb(ham.electron_electron_ida, orbital)
    @test r3b_self_coulomb ≈ R3B_OWNER_LOCAL_SELF_COULOMB atol = 1.0e-10

    println("r3a_h2_augmented_dimensions base=471 residual=18 augmented=489")
    println("r3a_h2_augmented_energies E_base=", E_base, " E_aug=", E_aug)
    println("r3b_h2_self_coulomb=", r3b_self_coulomb,
        " delta=", r3b_self_coulomb - R3B_OWNER_LOCAL_SELF_COULOMB)
    println("r3b_h2_vgm_errors direct=", direct_vgm_error, " pqs=", pqs_vgm_error)
end

println("cartesian_r3a_h2_augmented_one_body_elapsed_s=", elapsed)

facade_elapsed = @elapsed @testset "R3 H2 supplemented Hamiltonian facade" begin
    artifact = joinpath(mktempdir(), "r3_h2_supplemented.jld2")
    ham = GaussletBases.cartesian_residual_gto_mwg_hamiltonian(
        FACADE_SYSTEM; basis = FACADE_BASIS, supplement = FACADE_SUPPLEMENT,
        hamfile = artifact)
    @test ham isa CartesianIDAHamiltonian{Float64}
    @test size(ham.electron_electron_ida) == (489, 489)
    H = one_body_hamiltonian(ham)
    eig = eigen(Symmetric(H))
    orbital = eig.vectors[:, argmin(eig.values)]
    facade_self_coulomb = self_coulomb(ham.electron_electron_ida, orbital)
    @test facade_self_coulomb ≈ R3B_OWNER_LOCAL_SELF_COULOMB atol = 1.0e-10

    readback = read_cartesian_ida_hamiltonian(artifact)
    kinetic_delta = norm(readback.kinetic - ham.kinetic, Inf)
    V_delta = norm(readback.electron_electron_ida - ham.electron_electron_ida, Inf)
    one_body_delta = norm(one_body_hamiltonian(readback) - H, Inf)
    @test kinetic_delta <= 1.0e-12
    @test V_delta <= 1.0e-12
    @test one_body_delta <= 1.0e-12
    unit_delta = 0.0
    for (left, right) in zip(readback.nuclear_attraction_unit_by_center,
                             ham.nuclear_attraction_unit_by_center)
        unit_delta = max(unit_delta, norm(left - right, Inf))
    end
    @test unit_delta <= 1.0e-12

    required = (
        :provenance_version, :producer, :supplement_policy, :basis_by_center,
        :lmax, :uncontracted, :width_filtering, :candidate_count, :owner_counts,
        :base_dimension, :residual_dimension, :augmented_dimension,
        :augmented_basis_order, :residual_basis_convention, :rank_rule,
        :occupation_cutoff, :tau_neg_abs, :tau_neg_rel, :tau_merge_abs,
        :tau_merge_rel, :mwg_convention_version,
        :mwg_convention, :one_body_source, :interaction_source,
        :validation_check_labels, :h2_self_coulomb_reference,
    )
    JLD2.jldopen(artifact, "r") do file
        @test haskey(file, "supplement_provenance")
        @test haskey(file, "recipe_provenance")
        @test haskey(file, "hamiltonian_manifest")
        @test haskey(file, "hamiltonian_manifest/manifest_version")
        @test file["hamiltonian_manifest/manifest_version"] == 1
        @test file["recipe_provenance/route"] === :z_axis_diatomic_pqs_residual_gto_mwg
        @test file["recipe_provenance/nesting"] === :pqs
        @test file["recipe_provenance/producer"] === :cartesian_residual_gto_mwg_hamiltonian

        values = Dict{Symbol,Any}()
        for key in required
            path = "supplement_provenance/$(key)"
            @test haskey(file, path)
            values[key] = file[path]
        end
        @test values[:producer] === :cartesian_residual_gto_mwg_augmentation
        @test values[:supplement_policy] === :mwg_residual_gto
        @test values[:basis_by_center] == FACADE_SUPPLEMENT.basis_by_center
        @test values[:lmax] == 1
        @test values[:uncontracted] == false
        @test values[:width_filtering] === nothing
        @test values[:candidate_count] == 18
        @test values[:owner_counts] == [9, 9]
        @test values[:base_dimension] == 471
        @test values[:residual_dimension] == 18
        @test values[:augmented_dimension] == 489
        @test values[:augmented_basis_order] === :base_then_residual
        @test values[:residual_basis_convention] === :owner_local_residual_occupation_final_merge_lowdin
        @test values[:rank_rule] === :owner_local_residual_occupation
        @test values[:occupation_cutoff] == 1.0e-6
        @test values[:tau_merge_abs] == 1.0e-12
        @test values[:tau_merge_rel] == 1.0e-12
        @test values[:interaction_source] === :weight_aware_residual_mwg_ida_blocks
        @test values[:h2_self_coulomb_reference] == R3B_OWNER_LOCAL_SELF_COULOMB
    end

    @test_throws ArgumentError GaussletBases.cartesian_residual_gto_mwg_hamiltonian(
        merge(FACADE_SYSTEM, (; extra = true));
        basis = FACADE_BASIS, supplement = FACADE_SUPPLEMENT)
    @test_throws ArgumentError GaussletBases.cartesian_residual_gto_mwg_hamiltonian(
        FACADE_SYSTEM; basis = merge(FACADE_BASIS, (; radius = 4.0)),
        supplement = FACADE_SUPPLEMENT)
    @test_throws ArgumentError GaussletBases.cartesian_residual_gto_mwg_hamiltonian(
        FACADE_SYSTEM; basis = FACADE_BASIS,
        supplement = merge(FACADE_SUPPLEMENT, (; extra = true)))
    @test_throws ArgumentError GaussletBases.cartesian_residual_gto_mwg_hamiltonian(
        FACADE_SYSTEM; basis = FACADE_BASIS,
        supplement = merge(FACADE_SUPPLEMENT, (; width_filtering = (; min_width = 1.0))))
    @test_throws ArgumentError GaussletBases.cartesian_residual_gto_mwg_hamiltonian(
        merge(FACADE_SYSTEM, (; atom_locations = [(-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)]));
        basis = FACADE_BASIS, supplement = FACADE_SUPPLEMENT)
    @test_throws ArgumentError GaussletBases.cartesian_residual_gto_mwg_hamiltonian(
        merge(FACADE_SYSTEM, (; atom_locations = [(1.0, 0.0, -2.0), (1.0, 0.0, 2.0)]));
        basis = FACADE_BASIS, supplement = FACADE_SUPPLEMENT)
    @test_throws ArgumentError GaussletBases.cartesian_residual_gto_mwg_hamiltonian(
        (; atom_symbols = ["Cr", "Cr"], nuclear_charges = [24.0, 24.0],
         atom_locations = NUCLEI, nup = 24, ndn = 24);
        basis = FACADE_BASIS, supplement = FACADE_SUPPLEMENT)

    println("r3u_h2_facade_self_coulomb=", facade_self_coulomb,
        " delta=", facade_self_coulomb - R3B_OWNER_LOCAL_SELF_COULOMB)
    println("r3u_h2_facade_readback_deltas kinetic=", kinetic_delta,
        " unit_U=", unit_delta, " one_body=", one_body_delta, " V=", V_delta)
    println("r3u_h2_facade_artifact=", artifact)
end

println("cartesian_r3u_h2_facade_elapsed_s=", facade_elapsed)
