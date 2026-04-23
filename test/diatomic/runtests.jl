@testset "Bond-aligned diatomic QW reference path" begin
    basis14, operators14, check14 = _bond_aligned_diatomic_qw_fixture(; bond_length = 1.4)
    basis20, operators20, check20 = _bond_aligned_diatomic_qw_fixture(; bond_length = 2.0)

    @test basis14 isa BondAlignedDiatomicQWBasis3D
    @test operators14 isa QiuWhiteResidualGaussianOperators
    @test operators14.gaussian_data === nothing
    @test operators14.residual_count == 0
    @test operators14.gausslet_count == length(basis14.basis_x) * length(basis14.basis_y) * length(basis14.basis_z)
    @test norm(operators14.overlap - I, Inf) < 1.0e-8
    @test norm(operators20.overlap - I, Inf) < 1.0e-8
    @test operators14.one_body_hamiltonian ≈ transpose(operators14.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test operators14.interaction_matrix ≈ transpose(operators14.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test isfinite(check14.orbital_energy)
    @test isfinite(check14.vee_expectation)
    @test isfinite(check20.orbital_energy)
    @test isfinite(check20.vee_expectation)
    @test check14.orbital_energy < -1.0
    @test check20.orbital_energy < -1.0
    @test 0.5 < check14.vee_expectation < 1.0
    @test 0.5 < check20.vee_expectation < 1.0
    @test check14.orbital_energy < check20.orbital_energy
    @test length(basis14.basis_x) == length(basis14.basis_y)
    @test length(basis14.basis_z) > length(basis14.basis_x)
    @test length(basis20.basis_z) >= length(basis14.basis_z)
end

@testset "Per-center nuclear one-body reassembly on diatomic routes" begin
    basis, _operators, _check = _bond_aligned_diatomic_qw_fixture(; bond_length = 1.4)
    full_ops = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :by_center,
        interaction_treatment = :ggt_nearest,
    )
    atom_a_ops = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = [1.0, 0.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
    )
    atom_b_ops = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = [0.0, 1.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
    )

    @test full_ops.nuclear_term_storage == :by_center
    @test full_ops.nuclear_charges == [1.0, 1.0]
    @test !isnothing(full_ops.kinetic_one_body)
    @test !isnothing(full_ops.nuclear_one_body_by_center)
    @test length(full_ops.nuclear_one_body_by_center) == 2
    @test assembled_one_body_hamiltonian(full_ops) ≈
          full_ops.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test assembled_one_body_hamiltonian(full_ops; nuclear_charges = [1.0, 0.0]) ≈
          atom_a_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
    @test assembled_one_body_hamiltonian(full_ops; nuclear_charges = [0.0, 1.0]) ≈
          atom_b_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10

    full_ops_localized = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :by_center,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
    )
    atom_a_ops_localized = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = [1.0, 0.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
    )
    atom_b_ops_localized = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = [0.0, 1.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
    )

    @test full_ops_localized.gausslet_backend == :pgdg_localized_experimental
    @test full_ops_localized.nuclear_term_storage == :by_center
    @test !isnothing(full_ops_localized.kinetic_one_body)
    @test !isnothing(full_ops_localized.nuclear_one_body_by_center)
    @test length(full_ops_localized.nuclear_one_body_by_center) == 2
    @test assembled_one_body_hamiltonian(full_ops_localized) ≈
          full_ops_localized.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test assembled_one_body_hamiltonian(full_ops_localized; nuclear_charges = [1.0, 0.0]) ≈
          atom_a_ops_localized.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
    @test assembled_one_body_hamiltonian(full_ops_localized; nuclear_charges = [0.0, 1.0]) ≈
          atom_b_ops_localized.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10

    (
        _nested_basis,
        _nested_parent_ops,
        _nested_parent_check,
        expansion,
        _nested_source,
        fixed_block,
        _parent_modes,
        _parent_ground,
        _projected,
        _projected_vee,
        _capture,
        _projected_energy,
    ) = _bond_aligned_diatomic_nested_fixed_block_fixture(; bond_length = 1.4)

    full_nested_ops_localized = ordinary_cartesian_qiu_white_operators(
        fixed_block;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :by_center,
        expansion = expansion,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
    )
    atom_a_nested_ops_localized = ordinary_cartesian_qiu_white_operators(
        fixed_block;
        nuclear_charges = [1.0, 0.0],
        nuclear_term_storage = :total_only,
        expansion = expansion,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
    )
    atom_b_nested_ops_localized = ordinary_cartesian_qiu_white_operators(
        fixed_block;
        nuclear_charges = [0.0, 1.0],
        nuclear_term_storage = :total_only,
        expansion = expansion,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
    )

    @test full_nested_ops_localized.gausslet_backend == :pgdg_localized_experimental
    @test full_nested_ops_localized.nuclear_term_storage == :by_center
    @test !isnothing(full_nested_ops_localized.kinetic_one_body)
    @test !isnothing(full_nested_ops_localized.nuclear_one_body_by_center)
    @test length(full_nested_ops_localized.nuclear_one_body_by_center) == 2
    localized_nested_bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(
        fixed_block.parent_basis,
        expansion;
        gausslet_backend = :pgdg_localized_experimental,
    )
    localized_nested_factorized_basis =
        GaussletBases._nested_factorized_parent_basis(fixed_block)
    direct_contracted_nested_nuclear =
        GaussletBases._qwrg_bond_aligned_direct_contracted_nuclear_one_body_by_center(
            fixed_block.parent_basis,
            localized_nested_factorized_basis,
            localized_nested_bundles.bundle_x,
            localized_nested_bundles.bundle_y,
            localized_nested_bundles.bundle_z,
            expansion,
        )
    reference_nested_parent_nuclear = GaussletBases._qwrg_diatomic_nuclear_one_body_by_center(
        fixed_block.parent_basis,
        localized_nested_bundles.bundle_x,
        localized_nested_bundles.bundle_y,
        localized_nested_bundles.bundle_z,
        expansion,
    )
    reference_nested_contracted_nuclear = [
        GaussletBases._qwrg_contract_parent_symmetric_matrix(
            fixed_block.coefficient_matrix,
            matrix,
        ) for matrix in reference_nested_parent_nuclear
    ]
    @test length(direct_contracted_nested_nuclear) == length(reference_nested_contracted_nuclear)
    for nucleus_index in eachindex(direct_contracted_nested_nuclear, reference_nested_contracted_nuclear)
        @test direct_contracted_nested_nuclear[nucleus_index] ≈
              reference_nested_contracted_nuclear[nucleus_index] atol = 1.0e-10 rtol = 1.0e-10
    end
    @test assembled_one_body_hamiltonian(full_nested_ops_localized) ≈
          full_nested_ops_localized.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test assembled_one_body_hamiltonian(full_nested_ops_localized; nuclear_charges = [1.0, 0.0]) ≈
          atom_a_nested_ops_localized.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
    @test assembled_one_body_hamiltonian(full_nested_ops_localized; nuclear_charges = [0.0, 1.0]) ≈
          atom_b_nested_ops_localized.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
end

@testset "Nuclear term storage auto stays lightweight on longer-center routes" begin
    basis, _operators, _diagnostics = _bond_aligned_homonuclear_chain_qw_fixture(; natoms = 3)
    auto_ops = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = fill(1.0, length(basis.nuclei)),
        nuclear_term_storage = :auto,
        interaction_treatment = :ggt_nearest,
    )

    @test auto_ops.nuclear_term_storage == :total_only
    @test isnothing(auto_ops.kinetic_one_body)
    @test isnothing(auto_ops.nuclear_one_body_by_center)
    @test_throws ArgumentError assembled_one_body_hamiltonian(
        auto_ops;
        nuclear_charges = [1.0, 0.0, 1.0],
    )
end

@testset "Ordinary Cartesian naming surface distinguishes geometry from supplement" begin
    @test OrdinaryCartesianOrbital3D === QiuWhiteHybridOrbital3D
    @test OrdinaryCartesianOperators3D === QiuWhiteResidualGaussianOperators

    (
        direct_basis,
        _direct_parent_ops,
        _direct_parent_check,
        direct_supplement,
        direct_hybrid_ops,
        _direct_hybrid_check,
    ) = _bond_aligned_diatomic_hybrid_qw_fixture(; bond_length = 1.4)

    direct_via_clear_name = ordinary_cartesian_product_operators(
        direct_basis,
        direct_supplement;
        nuclear_charges = [1.0, 1.0],
        interaction_treatment = :ggt_nearest,
    )

    @test direct_via_clear_name isa OrdinaryCartesianOperators3D
    @test direct_via_clear_name.overlap ≈ direct_hybrid_ops.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test direct_via_clear_name.one_body_hamiltonian ≈ direct_hybrid_ops.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test direct_via_clear_name.interaction_matrix ≈ direct_hybrid_ops.interaction_matrix atol = 1.0e-12 rtol = 1.0e-12
    @test direct_via_clear_name.residual_count == direct_hybrid_ops.residual_count

    (
        _nested_basis,
        _nested_parent_ops,
        _nested_parent_check,
        nested_expansion,
        _nested_source,
        nested_fixed_block,
        nested_pure_ops,
        _nested_pure_check,
    ) = _bond_aligned_diatomic_nested_qw_fixture(; bond_length = 1.4)

    nested_via_clear_name = nested_cartesian_operators(
        nested_fixed_block;
        nuclear_charges = [1.0, 1.0],
        expansion = nested_expansion,
        interaction_treatment = :ggt_nearest,
    )

    @test nested_via_clear_name isa OrdinaryCartesianOperators3D
    @test nested_via_clear_name.gaussian_data === nothing
    @test nested_via_clear_name.overlap ≈ nested_pure_ops.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test nested_via_clear_name.one_body_hamiltonian ≈ nested_pure_ops.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test nested_via_clear_name.interaction_matrix ≈ nested_pure_ops.interaction_matrix atol = 1.0e-12 rtol = 1.0e-12
    @test nested_via_clear_name.residual_count == 0
end

@testset "Bond-aligned diatomic split geometry" begin
    basis, _operators, _check = _bond_aligned_diatomic_qw_fixture(; bond_length = 1.4)
    expansion = coulomb_gaussian_expansion(doacc = false)
    bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(basis, expansion)
    parent_box = (
        1:length(basis.basis_x),
        1:length(basis.basis_y),
        1:length(basis.basis_z),
    )
    working_box = (2:8, 2:8, 2:12)
    midpoint = 0.0

    geometry = GaussletBases._nested_bond_aligned_diatomic_split_geometry(
        bundles,
        parent_box,
        working_box;
        bond_axis = :z,
        midpoint = midpoint,
        nside = 5,
        min_parallel_to_transverse_ratio = 0.4,
    )
    sliver_geometry = GaussletBases._nested_bond_aligned_diatomic_split_geometry(
        bundles,
        parent_box,
        working_box;
        bond_axis = :z,
        midpoint = midpoint,
        nside = 5,
        min_parallel_to_transverse_ratio = 0.75,
    )
    short_geometry = GaussletBases._nested_bond_aligned_diatomic_split_geometry(
        bundles,
        parent_box,
        (3:7, 3:7, 4:12);
        bond_axis = :z,
        midpoint = midpoint,
        nside = 5,
        min_parallel_to_transverse_ratio = 0.4,
    )
    split_basis = bond_aligned_homonuclear_qw_basis(
        family = :G10,
        bond_length = 1.4,
        core_spacing = 0.3,
        xmax_parallel = 8.0,
        xmax_transverse = 6.0,
        bond_axis = :z,
        nuclear_charge = 2.0,
    )
    split_diagnostics = bond_aligned_diatomic_nested_geometry_diagnostics(
        split_basis;
        nside = 5,
    )
    split_bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(split_basis, expansion)
    split_parent_box = (
        1:length(split_basis.basis_x),
        1:length(split_basis.basis_y),
        1:length(split_basis.basis_z),
    )
    split_geometry = GaussletBases._nested_bond_aligned_diatomic_split_geometry(
        split_bundles,
        split_parent_box,
        split_diagnostics.geometry.working_box;
        bond_axis = :z,
        midpoint = midpoint,
        nside = 5,
        min_parallel_to_transverse_ratio = 0.4,
    )

    @test geometry.count_eligible
    @test !geometry.unsplit_aspect_eligible
    @test geometry.shape_eligible
    @test !geometry.did_split
    @test geometry.split_index == 7
    @test geometry.working_box == working_box
    @test isnothing(geometry.shared_midpoint_box)
    @test geometry.child_boxes == [(2:8, 2:8, 2:6), (2:8, 2:8, 8:12)]
    @test all(widths[3] >= 0.4 * max(widths[1], widths[2]) for widths in geometry.child_physical_widths)

    @test split_geometry.count_eligible
    @test split_geometry.unsplit_aspect_eligible
    @test split_geometry.shape_eligible
    @test split_geometry.did_split
    @test !isnothing(split_geometry.shared_midpoint_box)
    @test length(split_geometry.child_boxes) == 2

    @test sliver_geometry.count_eligible
    @test !sliver_geometry.shape_eligible
    @test !sliver_geometry.did_split

    @test !short_geometry.count_eligible
    @test !short_geometry.did_split
end

@testset "Bond-aligned diatomic adaptive retained counts are parity-neutral" begin
    @test GaussletBases._nested_diatomic_candidate_counts(11, 4) == collect(4:11)
    @test GaussletBases._nested_diatomic_candidate_counts(11, 5) == collect(5:11)

    symmetric_centers = [-2.0, -1.0, 0.0, 1.0, 2.0]
    @test GaussletBases._nested_doside_retained_count(symmetric_centers, 4) == 3
    @test GaussletBases._nested_doside_retained_count(
        symmetric_centers,
        4;
        enforce_symmetric_odd = false,
    ) == 4

    basis = bond_aligned_homonuclear_qw_basis(
        bond_length = 2.0,
        core_spacing = 0.5,
        xmax_parallel = 4.0,
        xmax_transverse = 3.0,
        bond_axis = :z,
    )
    source_nside5 = bond_aligned_diatomic_nested_fixed_source(basis; nside = 5)
    source_nside6 = bond_aligned_diatomic_nested_fixed_source(basis; nside = 6)
    counts_nside5 = sort!(
        unique(trace.retained_count for trace in GaussletBases._bond_aligned_diatomic_doside_traces(source_nside5)),
    )
    counts_nside6 = sort!(
        unique(trace.retained_count for trace in GaussletBases._bond_aligned_diatomic_doside_traces(source_nside6)),
    )

    @test any(iseven, counts_nside5)
    @test any(isodd, counts_nside5)
    @test any(iseven, counts_nside6)
    @test any(isodd, counts_nside6)
    @test 4 in counts_nside5
    @test 4 in counts_nside6
end

@testset "Bond-aligned heteronuclear split geometry" begin
    centers_axis = collect(-5.0:1.0:5.0)
    @test GaussletBases._nested_diatomic_split_plane_index(
        centers_axis,
        1:11,
        0.0;
        prefer_midpoint_tie_side = :left,
    ) == 6
    @test GaussletBases._nested_diatomic_split_plane_index(
        centers_axis,
        1:11,
        0.0;
        prefer_midpoint_tie_side = :right,
    ) == 5

    basis = bond_aligned_heteronuclear_qw_basis(
        atoms = ("He", "H"),
        bond_length = 1.45,
        core_spacings = (0.25, 0.5),
        nuclear_charges = (2.0, 1.0),
        xmax_parallel = 6.0,
        xmax_transverse = 4.0,
        bond_axis = :z,
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(basis, expansion)
    parent_box = (
        1:length(basis.basis_x),
        1:length(basis.basis_y),
        1:length(basis.basis_z),
    )
    working_box = (2:(length(basis.basis_x) - 1), 2:(length(basis.basis_y) - 1), 2:(length(basis.basis_z) - 1))
    midpoint = sum(GaussletBases._qwrg_axis_coordinate(nucleus, :z) for nucleus in basis.nuclei) / 2
    geometry = GaussletBases._nested_bond_aligned_diatomic_split_geometry(
        bundles,
        parent_box,
        working_box;
        bond_axis = :z,
        midpoint = midpoint,
        nside = 5,
        min_parallel_to_transverse_ratio = 0.4,
        use_midpoint_slab = false,
        prefer_midpoint_tie_side = :left,
    )

    @test geometry.count_eligible
    @test !geometry.unsplit_aspect_eligible
    @test geometry.shape_eligible
    @test !geometry.did_split
    @test geometry.shared_midpoint_box === nothing
    @test length(geometry.child_boxes) == 2
    @test sum(length(box[3]) for box in geometry.child_boxes) == length(working_box[3])
    @test all(widths[3] >= 0.4 * max(widths[1], widths[2]) for widths in geometry.child_physical_widths)
end

@testset "Bond-aligned diatomic nested fixed block" begin
    (
        basis,
        operators,
        check,
        expansion,
        source,
        fixed_block,
        parent_modes,
        _parent_ground,
        _projected,
        projected_vee,
        capture,
        projected_energy,
    ) = _bond_aligned_diatomic_nested_fixed_block_fixture(; bond_length = 1.4)

    @test source isa GaussletBases._CartesianNestedBondAlignedDiatomicSource3D
    @test fixed_block isa GaussletBases._NestedFixedBlock3D
    @test !source.geometry.did_split
    @test length(source.shared_shell_layers) >= 1
    @test length(source.child_sequences) == 1
    @test isnothing(source.geometry.shared_midpoint_box)
    @test isnothing(source.midpoint_slab_column_range)
    @test length(source.child_column_ranges) == 1
    @test length(source.child_column_ranges[1]) == size(source.child_sequences[1].coefficient_matrix, 2)
    @test size(fixed_block.overlap, 1) == size(source.sequence.coefficient_matrix, 2)
    @test source.sequence.working_box == (
        1:length(basis.basis_x),
        1:length(basis.basis_y),
        1:length(basis.basis_z),
    )
    @test length(source.sequence.support_indices) ==
        length(basis.basis_x) * length(basis.basis_y) * length(basis.basis_z)
    @test size(fixed_block.coefficient_matrix, 1) ==
        length(basis.basis_x) * length(basis.basis_y) * length(basis.basis_z)
    @test size(fixed_block.coefficient_matrix, 2) < size(fixed_block.coefficient_matrix, 1)
    @test source.sequence.packet.term_storage == :compact_production
    @test fixed_block.term_storage == :compact_production
    @test isnothing(fixed_block.gaussian_terms)
    @test isnothing(fixed_block.pair_terms)
    @test !isnothing(fixed_block.gaussian_sum)
    @test !isnothing(fixed_block.pair_sum)
    @test norm(fixed_block.overlap - I, Inf) < 1.0e-10
    @test all(isfinite, fixed_block.weights)
    @test minimum(fixed_block.weights) > 0.0
    @test all(isfinite, fixed_block.fixed_centers)
    @test capture > 0.998
    @test projected_energy < -1.2
    @test abs(projected_energy - parent_modes.values[1]) < 0.03
    @test 0.7 < projected_vee < 0.8
    @test abs(projected_vee - check.vee_expectation) < 5.0e-4
end

@testset "Bond-aligned diatomic compact nested fixed-block contract" begin
    (
        basis,
        _operators,
        _check,
        expansion,
        default_source,
        default_fixed_block,
        _parent_modes,
        _parent_ground,
        _projected,
        _projected_vee,
        _capture,
        _projected_energy,
    ) = _bond_aligned_diatomic_nested_fixed_block_fixture(; bond_length = 1.4)

    term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
    explicit_nested = bond_aligned_diatomic_nested_fixed_block(
        basis;
        expansion = expansion,
        term_coefficients = term_coefficients,
    )
    explicit_fixed_block = explicit_nested.fixed_block

    @test default_source.sequence.packet.term_storage == :compact_production
    @test default_fixed_block.term_storage == :compact_production
    @test isnothing(default_fixed_block.gaussian_terms)
    @test isnothing(default_fixed_block.pair_terms)
    @test !isnothing(default_fixed_block.gaussian_sum)
    @test !isnothing(default_fixed_block.pair_sum)

    @test explicit_fixed_block.term_storage == :compact_production
    @test isnothing(explicit_fixed_block.gaussian_terms)
    @test isnothing(explicit_fixed_block.pair_terms)
    @test !isnothing(explicit_fixed_block.gaussian_sum)
    @test !isnothing(explicit_fixed_block.pair_sum)

    @test explicit_nested.source.sequence.packet.term_storage == :compact_production

    @test norm(default_fixed_block.overlap - explicit_fixed_block.overlap, Inf) < 1.0e-12
    @test norm(default_fixed_block.coefficient_matrix - explicit_fixed_block.coefficient_matrix, Inf) < 1.0e-12
    @test norm(default_fixed_block.gaussian_sum - explicit_fixed_block.gaussian_sum, Inf) < 1.0e-12
    @test norm(default_fixed_block.pair_sum - explicit_fixed_block.pair_sum, Inf) < 1.0e-12

    @test GaussletBases._qwrg_fixed_block_one_body_matrix(default_fixed_block, expansion; Z = 1.0) ≈
        GaussletBases._qwrg_fixed_block_one_body_matrix(explicit_fixed_block, expansion; Z = 1.0) atol =
        1.0e-10 rtol = 1.0e-10
    @test GaussletBases._qwrg_fixed_block_interaction_matrix(default_fixed_block, expansion) ≈
        GaussletBases._qwrg_fixed_block_interaction_matrix(explicit_fixed_block, expansion) atol =
        1.0e-10 rtol = 1.0e-10

    @test_throws MethodError bond_aligned_diatomic_nested_fixed_source(
        basis;
        expansion = expansion,
        term_storage = :full_debug,
    )
    @test_throws MethodError bond_aligned_diatomic_nested_fixed_block(
        basis;
        expansion = expansion,
        term_storage = :full_debug,
    )
end

@testset "Bond-aligned diatomic nested QW consumer path" begin
    (
        _basis,
        parent_ops,
        parent_check,
        _expansion,
        _source,
        fixed_block,
        nested_ops,
        nested_check,
    ) = _bond_aligned_diatomic_nested_qw_fixture(; bond_length = 1.4)

    @test nested_ops isa GaussletBases.QiuWhiteResidualGaussianOperators
    @test nested_ops.basis === fixed_block
    @test nested_ops.gaussian_data === nothing
    @test nested_ops.interaction_treatment == :ggt_nearest
    @test nested_ops.residual_count == 0
    @test nested_ops.gausslet_count == size(fixed_block.overlap, 1)
    @test size(nested_ops.residual_centers) == (0, 3)
    @test size(nested_ops.residual_widths) == (0, 3)
    @test norm(nested_ops.overlap - I, Inf) < 1.0e-10
    @test norm(nested_ops.one_body_hamiltonian - transpose(nested_ops.one_body_hamiltonian), Inf) < 1.0e-12
    @test norm(nested_ops.interaction_matrix - transpose(nested_ops.interaction_matrix), Inf) < 1.0e-12
    @test all(orbital.kind == :nested_fixed for orbital in orbitals(nested_ops))
    @test all(isfinite, fixed_block.weights)
    @test minimum(fixed_block.weights) > 0.0
    @test nested_check.overlap_error < 1.0e-10
    @test nested_check.orbital_energy < -1.2
    @test 0.7 < nested_check.vee_expectation < 0.8
    @test abs(nested_check.orbital_energy - parent_check.orbital_energy) < 0.03
    @test abs(nested_check.vee_expectation - parent_check.vee_expectation) < 0.03
    @test parent_ops.residual_count == 0
end

@testset "Bond-aligned diatomic molecular supplement ordinary QW path" begin
    fixture = _bond_aligned_diatomic_hybrid_qw_fixture(; bond_length = 1.4)
    if fixture === nothing
        @test true
    else
        (
            _basis,
            parent_ops,
            parent_check,
            supplement,
            ordinary_ops,
            ordinary_check,
        ) = fixture

        @test supplement isa LegacyBondAlignedDiatomicGaussianSupplement
        @test supplement.atomic_source.lmax == 1
        @test ordinary_ops.gaussian_data === supplement
        @test ordinary_ops.interaction_treatment == :ggt_nearest
        @test ordinary_ops.residual_count > 0
        @test size(ordinary_ops.residual_centers, 1) == ordinary_ops.residual_count
        @test size(ordinary_ops.residual_widths, 1) == ordinary_ops.residual_count
        @test norm(ordinary_ops.overlap - I, Inf) < 1.0e-10
        @test ordinary_check.overlap_error < 1.0e-10
        @test ordinary_check.orbital_energy < -1.25
        @test 0.75 < ordinary_check.vee_expectation < 0.81
        @test abs(ordinary_check.orbital_energy - parent_check.orbital_energy) < 0.01
        @test abs(ordinary_check.vee_expectation - parent_check.vee_expectation) < 0.01
        @test parent_ops.residual_count == 0
    end
end

@testset "Bond-aligned heteronuclear molecular supplement ordinary QW path" begin
    (
        basis,
        parent_ops,
        parent_check,
        supplement,
        ordinary_ops,
        ordinary_check,
    ) = _bond_aligned_heteronuclear_hybrid_qw_fixture(; bond_length = 1.45)

    @test supplement isa LegacyBondAlignedHeteronuclearGaussianSupplement
    @test supplement.atomic_sources[1].atom == "He"
    @test supplement.atomic_sources[2].atom == "H"
    @test supplement.atomic_sources[1].basis_name == "cc-pVTZ"
    @test supplement.atomic_sources[2].basis_name == "cc-pVTZ"
    @test supplement.atomic_sources[1].lmax == 1
    @test supplement.atomic_sources[2].lmax == 1
    @test ordinary_ops.gaussian_data === supplement
    @test ordinary_ops.interaction_treatment == :ggt_nearest
    @test ordinary_ops.residual_count > 0
    @test size(ordinary_ops.residual_centers, 1) == ordinary_ops.residual_count
    @test size(ordinary_ops.residual_widths, 1) == ordinary_ops.residual_count
    @test norm(ordinary_ops.overlap - I, Inf) < 1.0e-10
    @test ordinary_check.overlap_error < 1.0e-10
    @test ordinary_check.orbital_energy < -2.0
    @test 0.8 < ordinary_check.vee_expectation < 1.5
    @test ordinary_check.orbital_energy < parent_check.orbital_energy
    @test ordinary_check.vee_expectation > parent_check.vee_expectation

    payload = bond_aligned_diatomic_geometry_payload(ordinary_ops)
    slice = bond_aligned_diatomic_plane_slice(
        payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-12,
    )
    @test payload.bond_axis == :z
    @test length(payload.nuclei) == 2
    @test count(point -> point.group_kind == :residual_gaussian, payload.points) == ordinary_ops.residual_count
    @test slice.selected_count > 0
    @test all(abs(point.y) <= slice.plane_tol for point in slice.points)
    @test parent_ops.residual_count == 0
end

@testset "Bond-aligned heteronuclear nested fixed block" begin
    (
        basis,
        parent_ops,
        parent_check,
        _supplement,
        ordinary_ops,
        ordinary_check,
        _expansion,
        source,
        fixed_block,
        _parent_modes,
        _parent_ground,
        _projected,
        projected_vee,
        capture,
        projected_energy,
    ) = _bond_aligned_heteronuclear_nested_fixed_block_fixture(; bond_length = 1.45)

    @test source isa GaussletBases._CartesianNestedBondAlignedDiatomicSource3D
    @test fixed_block isa GaussletBases._NestedFixedBlock3D
    @test !source.geometry.did_split
    @test source.geometry.shared_midpoint_box === nothing
    @test isnothing(source.midpoint_slab_column_range)
    @test length(source.shared_shell_layers) >= 1
    @test length(source.child_sequences) == 1
    @test length(source.child_column_ranges) == 1
    @test source.sequence.working_box == (
        1:length(basis.basis_x),
        1:length(basis.basis_y),
        1:length(basis.basis_z),
    )
    @test norm(fixed_block.overlap - I, Inf) < 1.0e-10
    @test all(isfinite, fixed_block.weights)
    @test minimum(fixed_block.weights) > 0.0
    @test all(isfinite, fixed_block.fixed_centers)
    @test parent_check.overlap_error < 1.0e-10
    @test ordinary_check.overlap_error < 1.0e-10
    @test capture > 0.99
    @test projected_energy < -2.0
    @test projected_vee > 0.8
    @test size(fixed_block.overlap, 1) < ordinary_ops.gausslet_count
    @test size(fixed_block.overlap, 1) < parent_ops.gausslet_count
    @test source.geometry.child_boxes[1][3] != source.geometry.child_boxes[2][3]
end

@testset "Bond-aligned heteronuclear nested QW consumer path" begin
    (
        _basis,
        parent_ops,
        parent_check,
        source,
        fixed_block,
        supplement,
        nested_ops,
        nested_check,
    ) = _bond_aligned_heteronuclear_nested_hybrid_qw_fixture(; bond_length = 1.45)

    @test nested_ops isa GaussletBases.QiuWhiteResidualGaussianOperators
    @test nested_ops.basis === fixed_block
    @test nested_ops.gaussian_data === supplement
    @test nested_ops.interaction_treatment == :ggt_nearest
    @test nested_ops.residual_count > 0
    @test nested_check.overlap_error < 1.0e-8
    @test nested_check.orbital_energy < -2.0
    @test nested_check.vee_expectation > 0.8
    @test abs(nested_check.orbital_energy - parent_check.orbital_energy) < 0.01
    @test nested_check.vee_expectation > parent_check.vee_expectation
    payload = bond_aligned_diatomic_geometry_payload(nested_ops, source)
    slice = bond_aligned_diatomic_plane_slice(
        payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-5,
    )
    @test count(point -> point.group_kind == :shared_midpoint_slab, payload.points) == 0
    @test count(point -> point.group_kind == :residual_gaussian, payload.points) == nested_ops.residual_count
    @test slice.selected_count > 0
    @test parent_ops.residual_count == 0
end

@testset "Bond-aligned diatomic molecular supplement nested QW path" begin
    fixture = _bond_aligned_diatomic_nested_hybrid_qw_fixture(; bond_length = 1.4)
    if fixture === nothing
        @test true
    else
        (
            _basis,
            parent_ops,
            parent_check,
            _source,
            fixed_block,
            supplement,
            nested_ops,
            nested_check,
        ) = fixture

        @test supplement isa LegacyBondAlignedDiatomicGaussianSupplement
        @test nested_ops.basis === fixed_block
        @test nested_ops.gaussian_data === supplement
        @test nested_ops.interaction_treatment == :ggt_nearest
        @test nested_ops.residual_count > 0
        @test size(nested_ops.residual_centers, 1) == nested_ops.residual_count
        @test size(nested_ops.residual_widths, 1) == nested_ops.residual_count
        @test norm(nested_ops.overlap - I, Inf) < 1.0e-10
        @test nested_check.overlap_error < 1.0e-10
        @test all(isfinite, fixed_block.weights)
        @test minimum(fixed_block.weights) > 0.0
        @test nested_check.orbital_energy < -1.25
        @test 0.75 < nested_check.vee_expectation < 0.81
        @test abs(nested_check.orbital_energy - parent_check.orbital_energy) < 0.03
        @test abs(nested_check.vee_expectation - parent_check.vee_expectation) < 0.01
        @test parent_ops.residual_count == 0
    end
end

@testset "Bond-aligned diatomic geometry payloads and plane slices" begin
    basis, ordinary_ops, _check = _bond_aligned_diatomic_qw_fixture(; bond_length = 1.4)
    basis_payload = bond_aligned_diatomic_geometry_payload(basis)
    ordinary_payload = bond_aligned_diatomic_geometry_payload(ordinary_ops)
    basis_slice = bond_aligned_diatomic_plane_slice(
        basis_payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-12,
    )

    @test basis_payload isa GaussletBases.BondAlignedDiatomicGeometryPayload3D
    @test length(basis_payload.nuclei) == 2
    @test basis_payload.bond_axis == :z
    @test length(basis_payload.points) ==
        length(basis.basis_x) * length(basis.basis_y) * length(basis.basis_z)
    @test Set(point.group_kind for point in basis_payload.points) == Set([:gausslet_product])
    @test length(basis_payload.box_outlines) == 1
    @test basis_payload.box_outlines[1].group_kind == :basis_box
    @test ordinary_payload.bond_axis == :z
    @test length(ordinary_payload.points) == ordinary_ops.gausslet_count
    @test Set(point.group_kind for point in ordinary_payload.points) == Set([:gausslet_product])
    @test basis_slice.plane_axis == :y
    @test basis_slice.plane_value == 0.0
    @test basis_slice.plane_tol == 1.0e-12
    @test basis_slice.total_count == length(basis_payload.points)
    expected_y_count = count(y -> abs(y) <= 1.0e-12, centers(basis.basis_y))
    @test basis_slice.selected_count == length(basis.basis_x) * expected_y_count * length(basis.basis_z)
    @test all(abs(point.y) <= basis_slice.plane_tol for point in basis_slice.points)
    @test length(basis_slice.nuclei) == 2

    (
        _nested_basis,
        _parent_ops,
        _parent_check,
        _expansion,
        source,
        fixed_block,
        _nested_ops,
        _nested_check,
    ) = _bond_aligned_diatomic_nested_qw_fixture(; bond_length = 1.4)
    nested_payload = bond_aligned_diatomic_geometry_payload(fixed_block, source)
    expected_nested_groups =
        source.geometry.did_split ?
        Set([:shared_shell_layer, :left_child, :shared_midpoint_slab, :right_child]) :
        Set([:shared_shell_layer, :shared_child])
    expected_nested_box_count =
        2 + length(source.geometry.child_boxes) + (isnothing(source.geometry.shared_midpoint_box) ? 0 : 1)

    @test nested_payload isa GaussletBases.BondAlignedDiatomicGeometryPayload3D
    @test length(nested_payload.points) == size(fixed_block.fixed_centers, 1)
    @test Set(point.group_kind for point in nested_payload.points) == expected_nested_groups
    @test length(nested_payload.box_outlines) == expected_nested_box_count
    @test length(nested_payload.shell_provenance) == length(source.shared_shell_layers)
    @test nested_payload.shell_provenance[1].source_box == source.shared_shell_layers[1].provenance.source_box
    @test nested_payload.shell_provenance[1].next_inner_box ==
        source.shared_shell_layers[1].provenance.next_inner_box
    @test nested_payload.box_outlines[1].group_kind == :parent_box
    @test nested_payload.box_outlines[2].group_kind == :working_box
    @test count(box -> box.group_kind == :child_box, nested_payload.box_outlines) == length(source.geometry.child_boxes)
    @test any(box -> box.group_kind == :shared_midpoint_slab_box, nested_payload.box_outlines) ==
        !isnothing(source.geometry.shared_midpoint_box)

    hybrid_fixture = _bond_aligned_diatomic_hybrid_qw_fixture(; bond_length = 1.4)
    @test hybrid_fixture !== nothing
    if hybrid_fixture !== nothing
        (
            _hybrid_basis,
            _hybrid_parent_ops,
            _hybrid_parent_check,
            _supplement,
            hybrid_ops,
            _hybrid_check,
        ) = hybrid_fixture
        hybrid_payload = bond_aligned_diatomic_geometry_payload(hybrid_ops)
        @test count(point -> point.group_kind == :residual_gaussian, hybrid_payload.points) == hybrid_ops.residual_count
    end

    nested_hybrid_fixture = _bond_aligned_diatomic_nested_hybrid_qw_fixture(; bond_length = 1.4)
    @test nested_hybrid_fixture !== nothing
    if nested_hybrid_fixture !== nothing
        (
            _basis2,
            _parent_ops2,
            _parent_check2,
            hybrid_source,
            _hybrid_fixed_block,
            _hybrid_supplement,
            hybrid_nested_ops,
            _hybrid_nested_check,
        ) = nested_hybrid_fixture
        hybrid_nested_payload = bond_aligned_diatomic_geometry_payload(hybrid_nested_ops, hybrid_source)
        expected_hybrid_nested_groups =
            hybrid_source.geometry.did_split ?
            Set([:shared_shell_layer, :left_child, :shared_midpoint_slab, :right_child, :residual_gaussian]) :
            Set([:shared_shell_layer, :shared_child, :residual_gaussian])
        @test count(point -> point.group_kind == :residual_gaussian, hybrid_nested_payload.points) ==
            hybrid_nested_ops.residual_count
        @test Set(point.group_kind for point in hybrid_nested_payload.points) ==
            expected_hybrid_nested_groups
    end
end

@testset "Bond-aligned diatomic raw source geometry and 3d export" begin
    nested_hybrid_fixture = _bond_aligned_diatomic_nested_hybrid_qw_fixture(; bond_length = 1.4)
    @test nested_hybrid_fixture !== nothing
    if nested_hybrid_fixture !== nothing
        (
            _basis,
            _parent_ops,
            _parent_check,
            source,
            _fixed_block,
            _supplement,
            hybrid_nested_ops,
            _hybrid_nested_check,
        ) = nested_hybrid_fixture

        fixed_payload = bond_aligned_diatomic_geometry_payload(hybrid_nested_ops, source)
        source_payload = bond_aligned_diatomic_source_geometry_payload(source)
        expected_source_groups =
            source.geometry.did_split ?
            Set([:shared_shell_region, :left_child_region, :shared_midpoint_slab_region, :right_child_region]) :
            Set([:shared_shell_region, :shared_child_region])
        expected_child_box_count = length(source.geometry.child_boxes)
        expected_shared_midpoint_box = !isnothing(source.geometry.shared_midpoint_box)
        expected_shared_shell_points = sum(length(layer.support_indices) for layer in source.shared_shell_layers)
        expected_left_child_points = source.geometry.did_split ? length(source.child_sequences[1].support_indices) : 0
        expected_right_child_points = source.geometry.did_split ? length(source.child_sequences[2].support_indices) : 0
        expected_shared_child_points = source.geometry.did_split ? 0 : length(source.child_sequences[1].support_indices)
        expected_midpoint_points =
            isnothing(source.geometry.shared_midpoint_box) ?
            0 :
            prod(length.(source.geometry.shared_midpoint_box))
        first_shell = source.shared_shell_layers[1].provenance
        first_payload_shell = source_payload.shell_provenance[1]
        first_shell_source_dims = string(
            length(first_shell.source_box[1]), "x",
            length(first_shell.source_box[2]), "x",
            length(first_shell.source_box[3]),
        )
        first_shell_next_inner_dims = string(
            length(first_shell.next_inner_box[1]), "x",
            length(first_shell.next_inner_box[2]), "x",
            length(first_shell.next_inner_box[3]),
        )

        @test all(isnothing(sequence.support_states) for sequence in source.child_sequences)
        @test Set(point.group_kind for point in source_payload.points) == expected_source_groups
        @test length(source_payload.points) == prod(length.(source.geometry.parent_box))
        @test length(source_payload.shell_provenance) == length(source.shared_shell_layers)
        @test first_shell.source_box == source.geometry.parent_box
        @test first_shell.source_point_count ==
            prod(length.(first_shell.source_box)) - prod(length.(first_shell.next_inner_box))
        @test first_shell.retained_fixed_count == length(source.sequence.layer_column_ranges[1])
        @test first_payload_shell.source_box == first_shell.source_box
        @test first_payload_shell.next_inner_box == first_shell.next_inner_box
        @test first_payload_shell.source_point_count == first_shell.source_point_count
        @test first_payload_shell.retained_fixed_count == first_shell.retained_fixed_count
        @test count(point -> point.group_kind == :shared_shell_region, source_payload.points) ==
            expected_shared_shell_points
        @test count(point -> point.group_kind == :left_child_region, source_payload.points) ==
            expected_left_child_points
        @test count(point -> point.group_kind == :right_child_region, source_payload.points) ==
            expected_right_child_points
        @test count(point -> point.group_kind == :shared_child_region, source_payload.points) ==
            expected_shared_child_points
        @test count(point -> point.group_kind == :shared_midpoint_slab_region, source_payload.points) ==
            expected_midpoint_points
        @test any(box.group_kind == :parent_box for box in source_payload.box_outlines)
        @test any(box.group_kind == :working_box for box in source_payload.box_outlines)
        @test count(box -> box.group_kind == :child_box, source_payload.box_outlines) == expected_child_box_count
        @test any(box -> box.group_kind == :shared_midpoint_slab_box, source_payload.box_outlines) ==
            expected_shared_midpoint_box

        fixed_slice = bond_aligned_diatomic_plane_slice(
            fixed_payload;
            plane_axis = :y,
            plane_value = 0.0,
            plane_tol = 1.0e-12,
        )
        debug_slice = bond_aligned_diatomic_plane_slice(
            fixed_payload;
            plane_axis = :y,
            plane_value = 0.0,
            plane_tol = 1.0e-5,
        )
        source_slice = bond_aligned_diatomic_plane_slice(
            source_payload;
            plane_axis = :y,
            plane_value = 0.0,
            plane_tol = 1.0e-5,
        )
        @test debug_slice.selected_count >= fixed_slice.selected_count
        @test source_slice.selected_count > debug_slice.selected_count
        @test source_slice.shell_provenance == source_payload.shell_provenance

        mktemp() do path, io
            close(io)
            payload = write_bond_aligned_diatomic_points3d(path, source_payload)
            text = read(path, String)
            @test payload === source_payload
            @test occursin("# bond_axis = z", text)
            @test occursin("# point_count = $(length(source_payload.points))", text)
            @test occursin("# nucleus_count = 2", text)
            @test occursin("# columns = x y z role kind group_kind group_id label", text)
            @test occursin("# shell label=shared_shell_1", text)
            @test occursin("source_box=$(first_shell_source_dims)", text)
            @test occursin("source_points=$(first_shell.source_point_count)", text)
            @test occursin("retained_fixed_count=$(first_shell.retained_fixed_count)", text)
            @test occursin("next_inner_box=$(first_shell_next_inner_dims)", text)
            @test occursin("# box label=parent_box", text)
            @test occursin("# box label=working_box", text)
            if source.geometry.did_split
                @test occursin("# box label=left_child_box", text)
                @test occursin("# box label=right_child_box", text)
                @test occursin("# box label=shared_midpoint_slab_box", text)
            else
                @test occursin("# box label=shared_child_box", text)
                @test !occursin("# box label=shared_midpoint_slab_box", text)
            end
            @test occursin("\tpoint\tsource_region\tshared_shell_region\t1\t", text)
            if source.geometry.did_split
                @test occursin("\tpoint\tsource_region\tleft_child_region\t1\t", text)
                @test occursin("\tpoint\tsource_region\tshared_midpoint_slab_region\t1\t", text)
                @test occursin("\tpoint\tsource_region\tright_child_region\t2\t", text)
            else
                @test occursin("\tpoint\tsource_region\tshared_child_region\t1\t", text)
            end
            @test occursin("\tnucleus\tnucleus\tnucleus\t1\tA", text)
            @test occursin("\tnucleus\tnucleus\tnucleus\t2\tB", text)
        end

        mktemp() do path, io
            close(io)
            slice = write_bond_aligned_diatomic_plane_projection(
                path,
                source_payload;
                plane_axis = :y,
                plane_value = 0.0,
                plane_tol = 1.0e-5,
            )
            text = read(path, String)
            @test slice.selected_count == source_slice.selected_count
            @test occursin("# shell label=shared_shell_1", text)
            @test occursin("source_box=$(first_shell_source_dims)", text)
            @test occursin("next_inner_box=$(first_shell_next_inner_dims)", text)
        end
    end
end

@testset "Bond-aligned diatomic doside / COMX trace diagnostics" begin
    (
        _basis,
        _parent_ops,
        _parent_check,
        _expansion,
        source,
        _fixed_block,
        _nested_ops,
        _nested_check,
    ) = _bond_aligned_diatomic_nested_qw_fixture(; bond_length = 1.4)

    traces = GaussletBases._bond_aligned_diatomic_doside_traces(
        source;
        symmetry_tol = 1.0e-8,
        zero_tol = 1.0e-8,
    )
    lost_center = filter(
        trace -> trace.symmetric_about_zero && !trace.contains_near_zero_center,
        traces,
    )
    expected_context_suffixes = Set([
        "face_xy/tangential_x",
        "face_xy/tangential_y",
        "face_xz/tangential_x",
        "face_xz/tangential_z",
        "face_yz/tangential_y",
        "face_yz/tangential_z",
        "edge_x/free_axis_x",
        "edge_y/free_axis_y",
        "edge_z/free_axis_z",
    ])

    @test length(traces) == 9 * length(source.shared_shell_layers)
    @test all(trace.group_kind == :shared_shell for trace in traces)
    @test Set(trace.layer_index for trace in traces) ==
        Set(1:length(source.shared_shell_layers))
    @test all(trace.symmetric_about_zero for trace in traces)
    @test isempty(lost_center)
    @test all(trace.contains_near_zero_center for trace in traces)
    for layer_index in 1:length(source.shared_shell_layers)
        layer_traces = filter(trace -> trace.layer_index == layer_index, traces)
        @test length(layer_traces) == 9
        @test Set(
            replace(
                trace.context_label,
                "shared_shell/layer_$(layer_index)/" => "",
            ) for trace in layer_traces
        ) == expected_context_suffixes
    end

    mktemp() do path, io
        close(io)
        written = GaussletBases._write_bond_aligned_diatomic_doside_trace(
            path,
            source;
            symmetry_tol = 1.0e-8,
            zero_tol = 1.0e-8,
        )
        text = read(path, String)
        @test length(written) == length(traces)
        @test occursin("# trace_count = $(length(traces))", text)
        @test occursin("# shared_shell layer=1 source_box=", text)
        @test occursin(" source_points=", text)
        @test occursin(" retained_fixed_count=", text)
        @test occursin(" next_inner_box=", text)
        @test occursin(
            "# note shared_child has no local side contractions; it remains a direct core block",
            text,
        )
        @test !occursin(
            "# note left_child has no local side contractions; it remains a direct core block",
            text,
        )
        @test !occursin(
            "# note right_child has no local side contractions; it remains a direct core block",
            text,
        )
        @test occursin("context_label = shared_shell/layer_1/face_xy/tangential_x", text)
        @test occursin(
            "context_label = shared_shell/layer_$(length(source.shared_shell_layers))/edge_z/free_axis_z",
            text,
        )
        @test occursin("parent_centers = [", text)
        @test occursin("localized_centers = [", text)
        @test occursin("contains_near_zero_center = true", text)
        @test occursin("even_retained_count =", text)
    end
end

@testset "Bond-aligned diatomic shared-shell parity experiment at R=1.4" begin
    (
        _basis,
        _parent_ops,
        _parent_check,
        _expansion,
        baseline_source,
        baseline_fixed_block,
        _baseline_parent_modes,
        _baseline_parent_ground,
        _baseline_projected,
        baseline_projected_vee,
        baseline_capture,
        baseline_projected_energy,
        _baseline_supplement,
        baseline_ops,
        baseline_check,
    ) = _bond_aligned_diatomic_nested_hybrid_qw_shared_shell_experiment_fixture(
        ;
        bond_length = 1.4,
        shared_shell_retain_xy = (4, 3),
        shared_shell_retain_xz = (4, 3),
        shared_shell_retain_yz = (4, 3),
    )
    (
        _basis2,
        _parent_ops2,
        _parent_check2,
        _expansion2,
        experiment_source,
        experiment_fixed_block,
        _experiment_parent_modes,
        _experiment_parent_ground,
        _experiment_projected,
        experiment_projected_vee,
        experiment_capture,
        experiment_projected_energy,
        _experiment_supplement,
        experiment_ops,
        experiment_check,
    ) = _bond_aligned_diatomic_nested_hybrid_qw_shared_shell_experiment_fixture(
        ;
        bond_length = 1.4,
    )

    baseline_traces = GaussletBases._bond_aligned_diatomic_doside_traces(baseline_source)
    experiment_traces = GaussletBases._bond_aligned_diatomic_doside_traces(experiment_source)
    baseline_trace_map = Dict(trace.context_label => trace for trace in baseline_traces)
    trace_map = Dict(trace.context_label => trace for trace in experiment_traces)
    targeted_contexts = (
        "shared_shell/layer_1/face_xy/tangential_x",
        "shared_shell/layer_1/face_xz/tangential_x",
        "shared_shell/layer_1/face_yz/tangential_y",
    )
    fixed_payload = bond_aligned_diatomic_geometry_payload(baseline_ops, baseline_source)
    experiment_payload = bond_aligned_diatomic_geometry_payload(experiment_ops, experiment_source)
    fixed_slice = bond_aligned_diatomic_plane_slice(
        fixed_payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-5,
    )
    experiment_slice = bond_aligned_diatomic_plane_slice(
        experiment_payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-5,
    )

    @test norm(experiment_fixed_block.overlap - I, Inf) < 1.0e-10
    @test all(isfinite, experiment_fixed_block.weights)
    @test minimum(experiment_fixed_block.weights) > 0.0
    @test size(experiment_fixed_block.overlap, 1) < size(baseline_fixed_block.overlap, 1)
    @test Set(baseline_trace_map[context].contains_near_zero_center for context in targeted_contexts) ==
        Set([false])
    @test Set(trace_map[context].contains_near_zero_center for context in targeted_contexts) ==
        Set([true])
    @test Set(baseline_trace_map[context].retained_count for context in targeted_contexts) ==
        Set([4])
    @test Set(trace_map[context].retained_count for context in targeted_contexts) ==
        Set([3])
    @test experiment_slice.selected_count > fixed_slice.selected_count > 0
    @test count(point -> point.group_kind == :shared_shell_layer, experiment_slice.points) >
        count(point -> point.group_kind == :shared_shell_layer, fixed_slice.points)
    @test count(point -> point.group_kind == :shared_shell_layer, experiment_slice.points) > 0
    @test abs(experiment_projected_vee - baseline_projected_vee) < 2.0e-5
    @test abs(experiment_capture - baseline_capture) < 2.0e-5
    @test abs(experiment_projected_energy - baseline_projected_energy) < 2.0e-5
    @test abs(experiment_check.orbital_energy - baseline_check.orbital_energy) < 2.0e-5
    @test abs(experiment_check.vee_expectation - baseline_check.vee_expectation) < 2.0e-5
end

@testset "Bond-aligned diatomic shared-shell parity experiment at R=2.0" begin
    (
        _basis,
        _parent_ops,
        _parent_check,
        _expansion,
        baseline_source,
        baseline_fixed_block,
        _baseline_parent_modes,
        _baseline_parent_ground,
        _baseline_projected,
        baseline_projected_vee,
        baseline_capture,
        baseline_projected_energy,
        _baseline_supplement,
        baseline_ops,
        baseline_check,
    ) = _bond_aligned_diatomic_nested_hybrid_qw_shared_shell_experiment_fixture(
        ;
        bond_length = 2.0,
        shared_shell_retain_xy = (4, 3),
        shared_shell_retain_xz = (4, 3),
        shared_shell_retain_yz = (4, 3),
    )
    (
        _basis2,
        _parent_ops2,
        _parent_check2,
        _expansion2,
        experiment_source,
        experiment_fixed_block,
        _experiment_parent_modes,
        _experiment_parent_ground,
        _experiment_projected,
        experiment_projected_vee,
        experiment_capture,
        experiment_projected_energy,
        _experiment_supplement,
        experiment_ops,
        experiment_check,
    ) = _bond_aligned_diatomic_nested_hybrid_qw_shared_shell_experiment_fixture(
        ;
        bond_length = 2.0,
    )

    baseline_traces = GaussletBases._bond_aligned_diatomic_doside_traces(baseline_source)
    experiment_traces = GaussletBases._bond_aligned_diatomic_doside_traces(experiment_source)
    baseline_trace_map = Dict(trace.context_label => trace for trace in baseline_traces)
    trace_map = Dict(trace.context_label => trace for trace in experiment_traces)
    targeted_contexts = (
        "shared_shell/layer_1/face_xy/tangential_x",
        "shared_shell/layer_1/face_xz/tangential_x",
        "shared_shell/layer_1/face_yz/tangential_y",
    )
    baseline_payload = bond_aligned_diatomic_geometry_payload(baseline_ops, baseline_source)
    experiment_payload = bond_aligned_diatomic_geometry_payload(experiment_ops, experiment_source)
    baseline_slice = bond_aligned_diatomic_plane_slice(
        baseline_payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-5,
    )
    experiment_slice = bond_aligned_diatomic_plane_slice(
        experiment_payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-5,
    )

    @test norm(experiment_fixed_block.overlap - I, Inf) < 1.0e-10
    @test all(isfinite, experiment_fixed_block.weights)
    @test minimum(experiment_fixed_block.weights) > 0.0
    @test size(experiment_fixed_block.overlap, 1) < size(baseline_fixed_block.overlap, 1)
    @test Set(baseline_trace_map[context].contains_near_zero_center for context in targeted_contexts) ==
        Set([false])
    @test Set(trace_map[context].contains_near_zero_center for context in targeted_contexts) ==
        Set([true])
    @test Set(baseline_trace_map[context].retained_count for context in targeted_contexts) ==
        Set([4])
    @test Set(trace_map[context].retained_count for context in targeted_contexts) ==
        Set([3])
    @test experiment_slice.selected_count > baseline_slice.selected_count > 0
    @test count(point -> point.group_kind == :shared_shell_layer, experiment_slice.points) >
        count(point -> point.group_kind == :shared_shell_layer, baseline_slice.points)
    @test count(point -> point.group_kind == :shared_shell_layer, experiment_slice.points) > 0
    @test abs(experiment_projected_vee - baseline_projected_vee) < 2.0e-5
    @test abs(experiment_capture - baseline_capture) < 2.0e-5
    @test abs(experiment_projected_energy - baseline_projected_energy) < 2.0e-5
    @test abs(experiment_check.orbital_energy - baseline_check.orbital_energy) < 2.0e-5
    @test abs(experiment_check.vee_expectation - baseline_check.vee_expectation) < 2.0e-5
end

@testset "Bond-aligned diatomic doside boundary correction on larger debug box" begin
    basis = bond_aligned_homonuclear_qw_basis(
        bond_length = 1.4,
        core_spacing = 0.5,
        xmax_parallel = 8,
        xmax_transverse = 5,
        bond_axis = :z,
    )
    source = bond_aligned_diatomic_nested_fixed_source(basis)
    traces = GaussletBases._bond_aligned_diatomic_doside_traces(source)

    @test length(traces) == 27
    symmetric_traces = filter(trace -> trace.symmetric_about_zero, traces)
    @test !isempty(symmetric_traces)
    @test any(trace -> iseven(trace.retained_count), symmetric_traces)
    @test any(trace -> isodd(trace.retained_count), symmetric_traces)
    @test any(trace -> trace.contains_near_zero_center, symmetric_traces)
    @test any(trace -> !trace.contains_near_zero_center, symmetric_traces)
    layer_1_shared_shell_traces = filter(
        trace -> startswith(trace.context_label, "shared_shell/layer_1/"),
        symmetric_traces,
    )
    @test Set(trace.context_label for trace in layer_1_shared_shell_traces if trace.contains_near_zero_center) ==
        Set([
            "shared_shell/layer_1/edge_x/free_axis_x",
            "shared_shell/layer_1/edge_y/free_axis_y",
            "shared_shell/layer_1/face_xy/tangential_x",
            "shared_shell/layer_1/face_xy/tangential_y",
            "shared_shell/layer_1/face_xz/tangential_x",
            "shared_shell/layer_1/face_yz/tangential_y",
        ])
    @test Set(trace.context_label for trace in layer_1_shared_shell_traces if !trace.contains_near_zero_center) ==
        Set([
            "shared_shell/layer_1/edge_z/free_axis_z",
            "shared_shell/layer_1/face_xz/tangential_z",
            "shared_shell/layer_1/face_yz/tangential_z",
        ])

    fixed_block = GaussletBases._nested_fixed_block(source)
    supplement = legacy_bond_aligned_diatomic_gaussian_supplement(
        "H",
        "cc-pVTZ",
        basis.nuclei;
        lmax = 1,
    )
    ops = ordinary_cartesian_qiu_white_operators(
        fixed_block,
        supplement;
        nuclear_charges = [1.0, 1.0],
        interaction_treatment = :ggt_nearest,
    )
    payload = bond_aligned_diatomic_geometry_payload(ops, source)
    slice = bond_aligned_diatomic_plane_slice(
        payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-5,
    )
    group_counts = Dict{Symbol,Int}()
    for point in slice.points
        group_counts[point.group_kind] = get(group_counts, point.group_kind, 0) + 1
    end

    @test slice.selected_count == 119
    @test get(group_counts, :shared_child, 0) == 43
    @test get(group_counts, :shared_shell_layer, 0) == 58
    @test get(group_counts, :residual_gaussian, 0) == 18
    @test get(group_counts, :left_child, 0) == 0
    @test get(group_counts, :shared_midpoint_slab, 0) == 0
    @test get(group_counts, :right_child, 0) == 0
end

@testset "Bond-aligned diatomic plane projection export" begin
    hybrid_fixture = _bond_aligned_diatomic_hybrid_qw_fixture(; bond_length = 1.4)
    @test hybrid_fixture !== nothing
    if hybrid_fixture !== nothing
        (
            _basis,
            _parent_ops,
            _parent_check,
            _supplement,
            hybrid_ops,
            _hybrid_check,
        ) = hybrid_fixture
        payload = bond_aligned_diatomic_geometry_payload(hybrid_ops)
        mktemp() do path, io
            close(io)
            slice = write_bond_aligned_diatomic_plane_projection(
                path,
                payload;
                plane_axis = :y,
                plane_value = 0.0,
                plane_tol = 1.0e-12,
            )
            text = read(path, String)
            @test slice.plane_axis == :y
            @test slice.plane_value == 0.0
            @test slice.plane_tol == 1.0e-12
            @test occursin("# plane_axis = y", text)
            @test occursin("# plane_value = 0.0", text)
            @test occursin("# plane_tol = 1.0e-12", text)
            @test occursin("# bond_axis = z", text)
            @test occursin("# selected_count = $(slice.selected_count)", text)
            @test occursin("# total_count = $(slice.total_count)", text)
            @test occursin("# projection_axes = x z", text)
            @test occursin("role=point group_kind=gausslet_product group_id=1", text)
            selected_residual_count = count(
                point -> point.group_kind == :residual_gaussian,
                slice.points,
            )
            @test count(==('@'), text) == 1 + selected_residual_count + 2
            @test sum(occursin("role=point group_kind=residual_gaussian", line) for line in split(text, '\n')) ==
                selected_residual_count
            gausslet_pos = findfirst("role=point group_kind=gausslet_product", text)
            residual_pos = findfirst("role=point group_kind=residual_gaussian", text)
            nucleus_pos = findfirst("role=nucleus group_kind=nucleus", text)
            @test gausslet_pos !== nothing
            @test residual_pos !== nothing
            @test nucleus_pos !== nothing
            @test gausslet_pos < residual_pos < nucleus_pos
        end
    end

    nested_hybrid_fixture = _bond_aligned_diatomic_nested_hybrid_qw_fixture(; bond_length = 1.4)
    @test nested_hybrid_fixture !== nothing
    if nested_hybrid_fixture !== nothing
        (
            _basis2,
            _parent_ops2,
            _parent_check2,
            source,
            _fixed_block,
            _supplement2,
            hybrid_nested_ops,
            _nested_check2,
        ) = nested_hybrid_fixture
        payload = bond_aligned_diatomic_geometry_payload(hybrid_nested_ops, source)
        mktemp() do path, io
            close(io)
            slice = write_bond_aligned_diatomic_plane_projection(
                path,
                payload;
                plane_axis = :y,
                plane_value = 0.0,
                plane_tol = 1.0e-12,
            )
            text = read(path, String)
            @test slice.selected_count <= slice.total_count
            @test occursin("role=point group_kind=shared_child group_id=1", text)
            @test occursin("role=point group_kind=shared_shell_layer group_id=1", text)
            selected_residual_count = count(
                point -> point.group_kind == :residual_gaussian,
                slice.points,
            )
            @test sum(occursin("role=point group_kind=residual_gaussian", line) for line in split(text, '\n')) ==
                selected_residual_count
            shared_child_pos = findfirst("role=point group_kind=shared_child", text)
            shared_pos = findfirst("role=point group_kind=shared_shell_layer", text)
            residual_pos = findfirst("role=point group_kind=residual_gaussian", text)
            nucleus_pos = findfirst("role=nucleus group_kind=nucleus", text)
            @test shared_child_pos !== nothing
            @test shared_pos !== nothing
            @test residual_pos !== nothing
            @test nucleus_pos !== nothing
            @test shared_child_pos < shared_pos < residual_pos < nucleus_pos
            @test count(==('@'), text) == 3 + selected_residual_count + 2
        end
    end
end
