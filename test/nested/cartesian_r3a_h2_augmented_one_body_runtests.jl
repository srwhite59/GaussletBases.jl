using GaussletBases
using LinearAlgebra
using Test

const NUCLEI = NTuple{3,Float64}[(0.0, 0.0, -2.0), (0.0, 0.0, 2.0)]
const EXPECTED_LABELS = [
    "a_s1", "a_s2", "a_s3", "a_px1", "a_py1", "a_pz1", "a_px2", "a_py2", "a_pz2",
    "b_s1", "b_s2", "b_s3", "b_px1", "b_py1", "b_pz1", "b_px2", "b_py2", "b_pz2",
]

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
    @test residual.base_dimension == 471
    @test residual.residual_dimension == 18
    @test size(operators.kinetic) == (489, 489)
    @test minimum(residual.residual_metric_eigenvalues) ≈ 3.0488355008683734e-4 atol = 1.0e-14
    @test maximum(residual.residual_metric_eigenvalues) ≈ 1.3512432621413795e-2 atol = 1.0e-14
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

    println("r3a_h2_augmented_dimensions base=471 residual=18 augmented=489")
    println("r3a_h2_augmented_energies E_base=", E_base, " E_aug=", E_aug)
end

println("cartesian_r3a_h2_augmented_one_body_elapsed_s=", elapsed)
