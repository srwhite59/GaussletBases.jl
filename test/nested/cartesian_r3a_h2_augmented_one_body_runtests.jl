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
const R3B_OWNER_LOCAL_SELF_COULOMB = 0.4574161883692301

function product_matrix(C, basis, ax, ay, az)
    matrix = zeros(Float64, basis.final_dimension, basis.final_dimension)
    C.assemble_terminal_product_operator!(matrix, basis, ax, ay, az)
    return matrix
end

function base_blocks(C, base)
    parent = base.parent
    basis = base.terminal_basis
    pgdg = Tuple(GaussletBases._nested_axis_pgdg(parent.parent_axis_bundle_object, axis)
        for axis in (:x, :y, :z))
    S = Tuple(axis.overlap for axis in pgdg)
    K = GaussletBases.cartesian_base_products(base).kinetic
    U = GaussletBases.cartesian_base_unit_nuclear(base)
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

function due_row_geometry(row)
    return (;
        row.terminal_order, row.role, row.region_kind, row.shell_index,
        row.index_ranges, row.physical_ranges, row.outer_box, row.inner_box,
        row.support_rows, row.slab_axis, row.slab_side, row.slab_thickness,
        row.slab_stack_index, row.slab_stack_count)
end

function check_override_geometry_parity(ordinary, refined, role, shell_index)
    due0 = ordinary.terminal_due_diligence
    due1 = refined.terminal_due_diligence
    @test due0.geometry.parent_axis_counts == due1.geometry.parent_axis_counts
    @test due0.geometry.parent_physical_bounds == due1.geometry.parent_physical_bounds
    @test due0.geometry.q == due1.geometry.q
    @test length(due0.terminal_rows) == length(due1.terminal_rows)
    matched = Int[]
    for (index, (row0, row1)) in enumerate(zip(due0.terminal_rows, due1.terminal_rows))
        @test due_row_geometry(row0) == due_row_geometry(row1)
        ismatch = row0.role === role && row0.shell_index == shell_index
        ismatch && push!(matched, index)
        if !ismatch
            @test row0.source_mode_shape == row1.source_mode_shape
            @test row0.retained_count == row1.retained_count
        end
    end
    return matched
end

function due_summary(due)
    return (;
        bounds = due.geometry.parent_physical_bounds,
        axes = due.geometry.parent_axis_counts,
        q = due.geometry.q,
        final_dimension = due.dimensions.base_final_dimension,
        retained = [row.retained_count for row in due.terminal_rows],
        topology = [row.region_kind for row in due.terminal_rows],
        source_shapes = [row.source_mode_shape for row in due.terminal_rows],
        row_warnings = [row.warning_summary for row in due.terminal_rows],
        warnings = due.warnings)
end

boundary_count(shape) = prod(shape) - prod(max(value - 2, 0) for value in shape)

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
    base_working = GaussletBases.cartesian_base_working_basis(
        FACADE_SYSTEM; basis = FACADE_BASIS, supplemented = true)
    parent, basis = base_working.parent, base_working.terminal_basis
    expansion = base_working.coulomb_expansion
    C = GaussletBases.CartesianFinalBasisRealization
    residual = C.pqs_terminal_residual_gto_augmentation(
        basis, parent.parent_axis_bundle_object, supplement, NUCLEI)
    operators = C.pqs_terminal_residual_gto_augmented_operators(
        basis, parent.parent_axis_bundle_object, parent.parent_basis_object,
        supplement, residual, NUCLEI, [1.0, 1.0]; expansion)

    X = C._terminal_residual_mixed_overlap(basis, parent.parent_axis_bundle_object, supplement)
    S_AA = GaussletBases._cartesian_supplement_cross_overlap(supplement, supplement)
    RSR = transpose(residual.T_G) * residual.T_G +
          transpose(residual.T_G) * X * residual.T_A +
          transpose(residual.T_A) * transpose(X) * residual.T_G +
          transpose(residual.T_A) * S_AA * residual.T_A
    numerical_residual = CRG.build_residual_gaussian_basis(
        basis.final_dimension, X, S_AA, residual.candidate_labels,
        residual.candidate_centers, residual.candidate_owner_indices;
        residual_occupation_cutoff = 1.0e-10,
        residual_injection_cutoff = 0.0,
        residual_compactness = nothing)
    numerical_operators = C.pqs_terminal_residual_gto_augmented_operators(
        basis, parent.parent_axis_bundle_object, parent.parent_basis_object,
        supplement, numerical_residual, NUCLEI, [1.0, 1.0]; expansion)

    @test residual.candidate_labels == EXPECTED_LABELS
    @test Dict(owner => count(==(owner), residual.candidate_owner_indices)
        for owner in unique(residual.candidate_owner_indices)) == Dict(1 => 9, 2 => 9)
    @test residual.residual_source_owner_indices == vcat(fill(1, 9), fill(2, 9))
    @test residual.owner_retained_counts == [9, 9]
    @test residual.occupation_cutoff == 1.0e-6
    @test residual.base_dimension == 487
    @test residual.residual_dimension == 18
    @test size(operators.kinetic) == (505, 505)
    @test minimum(residual.residual_occupations) ≈ 5.102500905664382e-4 atol = 1.0e-14
    @test maximum(residual.residual_occupations) ≈ 1.2243126230584132e-2 atol = 1.0e-14
    @test norm(residual.T_G + X * residual.T_A, Inf) <= 1.0e-10
    @test norm(RSR - I, Inf) <= 1.0e-10
    @test numerical_residual.occupation_cutoff == 1.0e-10
    @test numerical_residual.residual_injection_cutoff == 0.0
    @test isnothing(numerical_residual.injected_G)
    @test numerical_residual.residual_dimension == residual.residual_dimension
    @test numerical_residual.owner_retained_counts == residual.owner_retained_counts
    @test numerical_residual.T_G == residual.T_G
    @test numerical_residual.T_A == residual.T_A
    @test norm(numerical_operators.kinetic - operators.kinetic, Inf) <= 1.0e-10
    @test maximum(norm(a - b, Inf) for (a, b) in zip(
        numerical_operators.nuclear_attraction_unit_by_center,
        operators.nuclear_attraction_unit_by_center)) <= 1.0e-10
    Y_probe = zeros(Float64, size(S_AA, 1), 1)
    Y_probe[1] = inv(sqrt(S_AA[1, 1]))
    represented_probe = CRG.numerical_complete_reference_blocks_in_augmented_basis(
        numerical_residual, X, S_AA, [Y_probe], [[2.0]])
    @test represented_probe.diagnostics[1].recovery_loss <= 1.0e-10
    @test represented_probe.diagnostics[1].gram_error <= 1.0e-10
    @test represented_probe.diagnostics[1].electron_trace_error <= 1.0e-10
    @test_throws ArgumentError CRG.numerical_complete_reference_blocks_in_augmented_basis(
        residual, X, S_AA, [Y_probe], [[2.0]])

    for matrix in (
            operators.kinetic,
            operators.nuclear_attraction_unit_by_center...,
            operators.position.x, operators.position.y, operators.position.z,
            operators.x2.x, operators.x2.y, operators.x2.z)
        @test symmetry_error(matrix) <= 1.0e-9
    end

    base = base_blocks(C, base_working)
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

    base_ham = GaussletBases.cartesian_base_hamiltonian_assembly(base_working)
    ham = C.pqs_terminal_residual_gto_augmented_hamiltonian(
        base_ham, basis, parent.parent_axis_bundle_object, residual, operators;
        expansion)
    numerical_ham = C.pqs_terminal_residual_gto_augmented_hamiltonian(
        base_ham, basis, parent.parent_axis_bundle_object, numerical_residual,
        numerical_operators; expansion)
    @test ham isa CartesianIDAHamiltonian{Float64}
    @test size(ham.electron_electron_ida) == (505, 505)
    @test symmetry_error(ham.electron_electron_ida) <= 1.0e-10
    @test norm(ham.electron_electron_ida[1:nG, 1:nG] -
               base_ham.electron_electron_ida, Inf) <= 1.0e-12
    @test norm(one_body_hamiltonian(numerical_ham) - one_body_hamiltonian(ham), Inf) <= 1.0e-10
    @test norm(numerical_ham.electron_electron_ida -
        ham.electron_electron_ida, Inf) <= 1.0e-10

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

    println("r3a_h2_augmented_dimensions base=487 residual=18 augmented=505")
    println("r3a_h2_augmented_energies E_base=", E_base, " E_aug=", E_aug)
    println("r3b_h2_self_coulomb=", r3b_self_coulomb,
        " delta=", r3b_self_coulomb - R3B_OWNER_LOCAL_SELF_COULOMB)
    println("r3b_h2_vgm_errors direct=", direct_vgm_error, " pqs=", pqs_vgm_error)
    due = base_working.terminal_due_diligence
    println("r3a_h2_due_diligence=", (;
        bounds = due.geometry.parent_physical_bounds,
        axes = due.geometry.parent_axis_counts,
        padding = due.geometry.xmax_transverse,
        final_dimension = due.dimensions.base_final_dimension,
        retained = [row.retained_count for row in due.terminal_rows],
        topology = [row.region_kind for row in due.terminal_rows],
        source_shapes = [row.source_mode_shape for row in due.terminal_rows],
        row_warnings = [row.warning_summary for row in due.terminal_rows],
        warnings = due.warnings))

    empty_working = GaussletBases.cartesian_base_working_basis(
        FACADE_SYSTEM; basis = FACADE_BASIS, supplemented = true,
        source_mode_overrides = NamedTuple[])
    empty_ham = GaussletBases.cartesian_base_hamiltonian_assembly(empty_working)
    @test empty_working.terminal_basis.final_dimension == basis.final_dimension
    @test due_summary(empty_working.terminal_due_diligence) ==
        due_summary(base_working.terminal_due_diligence)
    @test one_body_hamiltonian(empty_ham) == one_body_hamiltonian(base_ham)
    @test empty_ham.electron_electron_ida == base_ham.electron_electron_ida

    shared_override = [(; role = :shared_molecular_shell,
        shell_index = 1, owner = :all, source_q = 6)]
    shared_working = GaussletBases.cartesian_base_working_basis(
        FACADE_SYSTEM; basis = FACADE_BASIS, supplemented = true,
        source_mode_overrides = shared_override)
    shared_rows = check_override_geometry_parity(
        base_working, shared_working, :shared_molecular_shell, 1)
    @test length(shared_rows) == 1
    shared_row = shared_working.terminal_due_diligence.terminal_rows[only(shared_rows)]
    @test shared_row.source_mode_shape == (6, 6, 7)
    @test shared_row.retained_count == boundary_count(shared_row.source_mode_shape)
    shared_ham = GaussletBases.cartesian_base_hamiltonian_assembly(shared_working)
    @test symmetry_error(one_body_hamiltonian(shared_ham)) <= 1.0e-10
    @test symmetry_error(shared_ham.electron_electron_ida) <= 1.0e-10

    separated_system = merge(FACADE_SYSTEM, (;
        atom_locations = [(0.0, 0.0, -6.0), (0.0, 0.0, 6.0)]))
    separated_basis = merge(FACADE_BASIS, (; xmax_parallel = 14.0))
    separated = GaussletBases.cartesian_base_working_basis(
        separated_system; basis = separated_basis, supplemented = true)
    atom_override = [(; role = :atom_local_shell,
        shell_index = 1, owner = :all, source_q = 6)]
    separated_refined = GaussletBases.cartesian_base_working_basis(
        separated_system; basis = separated_basis, supplemented = true,
        source_mode_overrides = atom_override)
    atom_rows = check_override_geometry_parity(
        separated, separated_refined, :atom_local_shell, 1)
    @test length(atom_rows) == 2
    @test all(index -> separated_refined.terminal_due_diligence.terminal_rows[index].source_mode_shape ==
        (6, 6, 6), atom_rows)
    @test all(index -> begin
        row = separated_refined.terminal_due_diligence.terminal_rows[index]
        row.retained_count == boundary_count(row.source_mode_shape)
    end, atom_rows)
    separated_ham = GaussletBases.cartesian_base_hamiltonian_assembly(separated_refined)
    @test symmetry_error(one_body_hamiltonian(separated_ham)) <= 1.0e-10
    @test symmetry_error(separated_ham.electron_electron_ida) <= 1.0e-10

    normalize = GaussletBases._pqs_source_box_route_driver_source_mode_overrides
    ordered = normalize(vcat(atom_override, shared_override), 5)
    @test ordered == normalize(vcat(shared_override, atom_override), 5)
    invalid_overrides = Any[
        [(; role = :atom_local_shell, shell_index = 0, owner = :all, source_q = 6)],
        [(; role = :atom_local_shell, shell_index = true, owner = :all, source_q = 6)],
        [(; role = :atom_local_shell, shell_index = 1.5, owner = :all, source_q = 6)],
        [(; role = :atom_local_shell, shell_index = -1, owner = :all, source_q = 6)],
        [(; role = :atom_local_shell, shell_index = 1, owner = :all, source_q = 5)],
        [(; role = :atom_local_shell, shell_index = 1, owner = :all, source_q = -1)],
        [(; role = :atom_local_shell, shell_index = 1, owner = :all, source_q = 6.5)],
        [(; role = :atom_local_shell, shell_index = 1, owner = :all, source_q = true)],
        [(; role = :atom_local_shell, shell_index = 1, owner = :left, source_q = 6)],
        [(; role = :midpoint_slab, shell_index = 1, owner = :all, source_q = 6)],
        [(; role = :atom_local_shell, shell_index = 1, owner = :all, source_q = 6, L = 7)],
        [(; role = :atom_local_shell, shell_index = 1, owner = :all,
            source_q = 6, source_mode_shape = (6, 6, 6))],
        [(; role = :atom_local_shell, shell_index = 1, owner = :all,
            source_q = 6, terminal_region_key = :terminal_region_3)],
        [(; role = :atom_local_shell, shell_index = 1, owner = :all,
            source_q = 6, order_index = 3)],
        vcat(atom_override, atom_override),
    ]
    for overrides in invalid_overrides
        @test_throws ArgumentError normalize(overrides, 5)
    end
    @test_throws ArgumentError GaussletBases.cartesian_base_working_basis(
        FACADE_SYSTEM; basis = FACADE_BASIS, supplemented = true,
        source_mode_overrides = [(; role = :shared_molecular_shell,
            shell_index = 99, owner = :all, source_q = 6)])
    @test_throws ArgumentError GaussletBases.cartesian_base_working_basis(
        FACADE_SYSTEM; basis = merge(FACADE_BASIS, (; nesting = :wl)),
        supplemented = true, source_mode_overrides = shared_override)
    @test_throws ArgumentError GaussletBases.cartesian_base_working_basis(
        FACADE_SYSTEM; basis = merge(FACADE_BASIS, (; source_span = :mapped_comx)),
        supplemented = true, source_mode_overrides = shared_override)
    atom_system = (; atom_symbols = ["H"], nuclear_charges = [1.0],
        atom_locations = [(0.0, 0.0, 0.0)], nup = 1, ndn = 0)
    @test_throws ArgumentError GaussletBases.cartesian_base_working_basis(
        atom_system; basis = (; q = 5, core_spacing = 0.5, radius = 4.0),
        supplemented = true, source_mode_overrides = atom_override)
    @test_throws ArgumentError GaussletBases._plb_build_member(
        (; source_mode_overrides = shared_override), NamedTuple[])
    @test_throws ArgumentError GaussletBases._plb_build_numerical_complete_additive_reference_member(
        (; source_mode_overrides = shared_override), NamedTuple[], Any[])
    @test_throws ArgumentError GaussletBases._plb_build_numerical_complete_additive_reference_member(
        (; source_mode_overrides = shared_override, source_span = :mapped_comx),
        NamedTuple[], Any[nothing])

    println("source_q_override_h2_due_diligence=", (;
        shared = due_summary(shared_working.terminal_due_diligence),
        atom_local = due_summary(separated_refined.terminal_due_diligence)))
end

println("cartesian_r3a_h2_augmented_one_body_elapsed_s=", elapsed)

facade_elapsed = @elapsed @testset "R3 H2 supplemented Hamiltonian facade" begin
    artifact = joinpath(mktempdir(), "r3_h2_supplemented.jld2")
    ham = GaussletBases.cartesian_residual_gto_mwg_hamiltonian(
        FACADE_SYSTEM; basis = FACADE_BASIS, supplement = FACADE_SUPPLEMENT,
        hamfile = artifact)
    @test ham isa CartesianIDAHamiltonian{Float64}
    @test size(ham.electron_electron_ida) == (505, 505)
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
        @test file["coulomb_expansion/policy"] === :compact
        @test file["coulomb_expansion/term_count"] == 45

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
        @test values[:base_dimension] == 487
        @test values[:residual_dimension] == 18
        @test values[:augmented_dimension] == 505
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
