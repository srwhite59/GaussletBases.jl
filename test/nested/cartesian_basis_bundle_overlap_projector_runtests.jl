# Integration/slow test. Do not include in default nested runner.

@testset "Cartesian basis bundle overlap and projector" begin
    square_basis, _square_source, square_fixed_block, _square_diagnostics =
        _axis_aligned_homonuclear_square_lattice_nested_fixture()
    square_basis_rep = basis_representation(square_basis)
    square_fixed_rep = basis_representation(square_fixed_block)

    diatomic_basis14, diatomic_ops14, _diatomic_check14 =
        _bond_aligned_diatomic_qw_fixture(; bond_length = 1.4)
    diatomic_basis20, _diatomic_ops20, _diatomic_check20 =
        _bond_aligned_diatomic_qw_fixture(; bond_length = 2.0)
    diatomic_rep14 = basis_representation(diatomic_basis14)
    diatomic_rep20 = basis_representation(diatomic_basis20)
    bond_aligned_hybrid_fixture = _bond_aligned_diatomic_nested_hybrid_bundle_fixture()

    dir = mktempdir()
    try
        square_path = joinpath(dir, "square_basis.jld2")
        square_fixed_path = joinpath(dir, "square_fixed.jld2")
        diatomic14_path = joinpath(dir, "diatomic14.jld2")
        diatomic20_path = joinpath(dir, "diatomic20.jld2")
        diatomic_ops_path = joinpath(dir, "diatomic_ops.jld2")
        atomic_fixed_full_path = joinpath(dir, "atomic_fixed_full.jld2")
        atomic_hybrid_full_path = joinpath(dir, "atomic_hybrid_full.jld2")
        atomic_hybrid_legacy_path = joinpath(dir, "atomic_hybrid_legacy.jld2")
        bond_aligned_hybrid_fixed_path = joinpath(dir, "bond_aligned_hybrid_fixed.jld2")
        bond_aligned_hybrid_path = joinpath(dir, "bond_aligned_hybrid_ops.jld2")

        write_cartesian_basis_bundle_jld2(square_path, square_basis)
        write_cartesian_basis_bundle_jld2(square_fixed_path, square_fixed_block)
        write_cartesian_basis_bundle_jld2(diatomic14_path, diatomic_basis14)
        write_cartesian_basis_bundle_jld2(diatomic20_path, diatomic_basis20)
        write_cartesian_basis_bundle_jld2(diatomic_ops_path, diatomic_ops14)
        hybrid_fixture = _atomic_hybrid_cartesian_representation_fixture()
        write_cartesian_basis_bundle_jld2(atomic_fixed_full_path, hybrid_fixture.fixed_full)
        write_cartesian_basis_bundle_jld2(atomic_hybrid_full_path, hybrid_fixture.full_ops)
        write_cartesian_basis_bundle_jld2(atomic_hybrid_legacy_path, hybrid_fixture.legacy_ops)
        write_cartesian_basis_bundle_jld2(
            bond_aligned_hybrid_fixed_path,
            bond_aligned_hybrid_fixture.fixed_block,
        )
        write_cartesian_basis_bundle_jld2(
            bond_aligned_hybrid_path,
            bond_aligned_hybrid_fixture.hybrid_ops,
        )

        square_bundle = read_cartesian_basis_bundle(square_path)
        square_fixed_bundle = read_cartesian_basis_bundle(square_fixed_path)
        diatomic14_bundle = read_cartesian_basis_bundle(diatomic14_path)
        diatomic20_bundle = read_cartesian_basis_bundle(diatomic20_path)
        diatomic_ops_bundle = read_cartesian_basis_bundle(diatomic_ops_path)
        atomic_hybrid_full_bundle = read_cartesian_basis_bundle(atomic_hybrid_full_path)
        atomic_hybrid_legacy_bundle = read_cartesian_basis_bundle(atomic_hybrid_legacy_path)
        bond_aligned_hybrid_bundle = read_cartesian_basis_bundle(bond_aligned_hybrid_path)

        @test square_bundle.path == abspath(square_path)
        @test square_bundle.diagnostics.basis_kind == :direct_product
        @test square_bundle.diagnostics.final_dimension == square_basis_rep.metadata.final_dimension
        @test square_bundle.ham === nothing
        @test diatomic_ops_bundle.ham !== nothing
        @test diatomic_ops_bundle.ham["model_kind"] == "ordinary_cartesian_operators"
        @test diatomic_ops_bundle.diagnostics.has_ham

        diatomic_atom_a_ops = ordinary_cartesian_qiu_white_operators(
            diatomic_basis14;
            nuclear_charges = [1.0, 0.0],
            nuclear_term_storage = :total_only,
            interaction_treatment = :ggt_nearest,
        )
        diatomic_atom_b_ops = ordinary_cartesian_qiu_white_operators(
            diatomic_basis14;
            nuclear_charges = [0.0, 1.0],
            nuclear_term_storage = :total_only,
            interaction_treatment = :ggt_nearest,
        )
        @test assembled_one_body_hamiltonian(diatomic_ops14) ≈
              diatomic_ops14.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
        @test assembled_one_body_hamiltonian(diatomic_ops_bundle) ≈
              diatomic_ops14.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
        @test assembled_one_body_hamiltonian(diatomic_ops14; nuclear_charges = [1.0, 0.0]) ≈
              diatomic_atom_a_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(diatomic_ops_bundle; nuclear_charges = [1.0, 0.0]) ≈
              diatomic_atom_a_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(diatomic_ops14; nuclear_charges = [0.0, 1.0]) ≈
              diatomic_atom_b_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(diatomic_ops_bundle; nuclear_charges = [0.0, 1.0]) ≈
              diatomic_atom_b_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10

        hybrid_atom_a_ops = ordinary_cartesian_qiu_white_operators(
            bond_aligned_hybrid_fixture.fixed_block,
            bond_aligned_hybrid_fixture.supplement;
            nuclear_charges = [1.0, 0.0],
            nuclear_term_storage = :total_only,
            interaction_treatment = :ggt_nearest,
        )
        hybrid_atom_b_ops = ordinary_cartesian_qiu_white_operators(
            bond_aligned_hybrid_fixture.fixed_block,
            bond_aligned_hybrid_fixture.supplement;
            nuclear_charges = [0.0, 1.0],
            nuclear_term_storage = :total_only,
            interaction_treatment = :ggt_nearest,
        )
        @test assembled_one_body_hamiltonian(bond_aligned_hybrid_fixture.hybrid_ops) ≈
              bond_aligned_hybrid_fixture.hybrid_ops.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
        @test assembled_one_body_hamiltonian(bond_aligned_hybrid_bundle) ≈
              bond_aligned_hybrid_fixture.hybrid_ops.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
        @test assembled_one_body_hamiltonian(
                  bond_aligned_hybrid_fixture.hybrid_ops;
                  nuclear_charges = [1.0, 0.0],
              ) ≈ hybrid_atom_a_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(
                  bond_aligned_hybrid_bundle;
                  nuclear_charges = [1.0, 0.0],
              ) ≈ hybrid_atom_a_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(
                  bond_aligned_hybrid_fixture.hybrid_ops;
                  nuclear_charges = [0.0, 1.0],
              ) ≈ hybrid_atom_b_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(
                  bond_aligned_hybrid_bundle;
                  nuclear_charges = [0.0, 1.0],
              ) ≈ hybrid_atom_b_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10

        loaded_square_rep = load_cartesian_basis_representation(square_path)
        @test loaded_square_rep.metadata.basis_kind == square_basis_rep.metadata.basis_kind
        @test loaded_square_rep.metadata.final_dimension == square_basis_rep.metadata.final_dimension
        @test loaded_square_rep.metadata.parent_kind == square_basis_rep.metadata.parent_kind

        direct_self_disk = cross_overlap(square_bundle, square_bundle)
        direct_cross_disk = cross_overlap(diatomic14_bundle, diatomic20_bundle)
        nested_cross_disk = cross_overlap(square_path, square_fixed_path)

        @test direct_self_disk ≈ cross_overlap(square_basis_rep, square_basis_rep) atol = 1.0e-10 rtol = 1.0e-10
        @test direct_cross_disk ≈ cross_overlap(diatomic_rep14, diatomic_rep20) atol = 1.0e-10 rtol = 1.0e-10
        @test nested_cross_disk ≈ cross_overlap(square_basis_rep, square_fixed_rep) atol = 1.0e-10 rtol = 1.0e-10

        disk_projector = basis_projector(square_fixed_path, square_path)
        memory_projector = basis_projector(square_fixed_rep, square_basis_rep)
        @test disk_projector.matrix ≈ memory_projector.matrix atol = 1.0e-10 rtol = 1.0e-10
        @test disk_projector.diagnostics.transfer_path == memory_projector.diagnostics.transfer_path

        fixed_coefficients = cos.(Float64.(1:square_fixed_rep.metadata.final_dimension))
        disk_transfer = transfer_orbitals(fixed_coefficients, square_fixed_path, square_path)
        memory_transfer = transfer_orbitals(fixed_coefficients, square_fixed_rep, square_basis_rep)
        @test disk_transfer.coefficients ≈ memory_transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test disk_transfer.diagnostics.transferred_residual_inf < 1.0e-10

        @test atomic_hybrid_full_bundle.diagnostics.parent_kind == :cartesian_plus_supplement_raw
        @test atomic_hybrid_full_bundle.ham !== nothing
        @test hasproperty(
            atomic_hybrid_full_bundle.basis.parent_data,
            :cartesian_supplement_axis_tables,
        )
        @test bond_aligned_hybrid_bundle.diagnostics.parent_kind == :cartesian_plus_supplement_raw
        @test bond_aligned_hybrid_bundle.ham !== nothing
        @test hasproperty(
            bond_aligned_hybrid_bundle.basis.parent_data,
            :cartesian_supplement_axis_tables,
        )

        disk_hybrid_self = cross_overlap(atomic_hybrid_full_path, atomic_hybrid_full_path)
        memory_hybrid_self = cross_overlap(hybrid_fixture.full_rep, hybrid_fixture.full_rep)
        @test disk_hybrid_self ≈ memory_hybrid_self atol = 1.0e-10 rtol = 1.0e-10

        disk_hybrid_cross = cross_overlap(atomic_hybrid_full_path, atomic_hybrid_legacy_path)
        memory_hybrid_cross = cross_overlap(hybrid_fixture.full_rep, hybrid_fixture.legacy_rep)
        @test disk_hybrid_cross ≈ memory_hybrid_cross atol = 1.0e-10 rtol = 1.0e-10

        disk_hybrid_parent = cross_overlap(atomic_fixed_full_path, atomic_hybrid_full_path)
        memory_hybrid_parent = cross_overlap(hybrid_fixture.fixed_full_rep, hybrid_fixture.full_rep)
        @test disk_hybrid_parent ≈ memory_hybrid_parent atol = 1.0e-10 rtol = 1.0e-10

        disk_hybrid_projector =
            basis_projector(atomic_hybrid_full_path, atomic_hybrid_legacy_path)
        memory_hybrid_projector =
            basis_projector(hybrid_fixture.full_rep, hybrid_fixture.legacy_rep)
        @test disk_hybrid_projector.matrix ≈ memory_hybrid_projector.matrix atol = 1.0e-10 rtol = 1.0e-10
        @test disk_hybrid_projector.diagnostics.transfer_path ==
            memory_hybrid_projector.diagnostics.transfer_path

        hybrid_coefficients = cos.(Float64.(1:hybrid_fixture.full_rep.metadata.final_dimension))
        disk_hybrid_transfer = transfer_orbitals(
            hybrid_coefficients,
            atomic_hybrid_full_path,
            atomic_hybrid_legacy_path,
        )
        memory_hybrid_transfer = transfer_orbitals(
            hybrid_coefficients,
            hybrid_fixture.full_rep,
            hybrid_fixture.legacy_rep,
        )
        disk_hybrid_fast_transfer = transfer_orbitals(
            hybrid_coefficients,
            atomic_hybrid_full_path,
            atomic_hybrid_legacy_path;
            materialize_projector = false,
        )
        bundle_hybrid_fast_transfer = transfer_orbitals(
            hybrid_coefficients,
            atomic_hybrid_full_bundle,
            atomic_hybrid_legacy_bundle;
            materialize_projector = false,
        )
        memory_hybrid_fast_transfer = transfer_orbitals(
            hybrid_coefficients,
            hybrid_fixture.full_rep,
            hybrid_fixture.legacy_rep;
            materialize_projector = false,
        )
        @test disk_hybrid_transfer.coefficients ≈
            memory_hybrid_transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test disk_hybrid_transfer.projector !== nothing
        @test memory_hybrid_transfer.projector !== nothing
        @test disk_hybrid_transfer.projector.matrix ≈
            memory_hybrid_transfer.projector.matrix atol = 1.0e-10 rtol = 1.0e-10
        @test disk_hybrid_transfer.diagnostics.transferred_residual_inf < 1.0e-10
        @test memory_hybrid_fast_transfer.coefficients ≈
            memory_hybrid_transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test disk_hybrid_fast_transfer.coefficients ≈
            memory_hybrid_fast_transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test bundle_hybrid_fast_transfer.coefficients ≈
            memory_hybrid_fast_transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test disk_hybrid_fast_transfer.projector === nothing
        @test bundle_hybrid_fast_transfer.projector === nothing
        @test memory_hybrid_fast_transfer.projector === nothing
        @test disk_hybrid_fast_transfer.diagnostics.transfer_path == :hybrid_mixed_raw_cross_overlap_transfer
        @test bundle_hybrid_fast_transfer.diagnostics.transfer_path == :hybrid_mixed_raw_cross_overlap_transfer
        @test isnan(disk_hybrid_fast_transfer.diagnostics.projector_residual_inf)
        @test isnan(bundle_hybrid_fast_transfer.diagnostics.projector_residual_inf)
        @test isnan(disk_hybrid_fast_transfer.diagnostics.transferred_residual_inf)
        @test isnan(bundle_hybrid_fast_transfer.diagnostics.transferred_residual_inf)

        disk_bond_aligned_hybrid_self =
            cross_overlap(bond_aligned_hybrid_path, bond_aligned_hybrid_path)
        memory_bond_aligned_hybrid_self =
            cross_overlap(
                bond_aligned_hybrid_fixture.hybrid_rep,
                bond_aligned_hybrid_fixture.hybrid_rep,
            )
        @test disk_bond_aligned_hybrid_self ≈
            memory_bond_aligned_hybrid_self atol = 1.0e-10 rtol = 1.0e-10

        disk_bond_aligned_hybrid_cross =
            cross_overlap(bond_aligned_hybrid_fixed_path, bond_aligned_hybrid_path)
        memory_bond_aligned_hybrid_cross =
            cross_overlap(
                bond_aligned_hybrid_fixture.fixed_rep,
                bond_aligned_hybrid_fixture.hybrid_rep,
            )
        @test disk_bond_aligned_hybrid_cross ≈
            memory_bond_aligned_hybrid_cross atol = 1.0e-10 rtol = 1.0e-10

        disk_bond_aligned_hybrid_projector =
            basis_projector(bond_aligned_hybrid_fixed_path, bond_aligned_hybrid_path)
        memory_bond_aligned_hybrid_projector =
            basis_projector(
                bond_aligned_hybrid_fixture.fixed_rep,
                bond_aligned_hybrid_fixture.hybrid_rep,
            )
        @test disk_bond_aligned_hybrid_projector.matrix ≈
            memory_bond_aligned_hybrid_projector.matrix atol = 1.0e-10 rtol = 1.0e-10
        @test disk_bond_aligned_hybrid_projector.diagnostics.transfer_path ==
            memory_bond_aligned_hybrid_projector.diagnostics.transfer_path

        bond_aligned_coefficients =
            cos.(Float64.(1:bond_aligned_hybrid_fixture.fixed_rep.metadata.final_dimension))
        disk_bond_aligned_hybrid_transfer = transfer_orbitals(
            bond_aligned_coefficients,
            bond_aligned_hybrid_fixed_path,
            bond_aligned_hybrid_path,
        )
        memory_bond_aligned_hybrid_transfer = transfer_orbitals(
            bond_aligned_coefficients,
            bond_aligned_hybrid_fixture.fixed_rep,
            bond_aligned_hybrid_fixture.hybrid_rep,
        )
        @test disk_bond_aligned_hybrid_transfer.coefficients ≈
            memory_bond_aligned_hybrid_transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test disk_bond_aligned_hybrid_transfer.projector.matrix ≈
            memory_bond_aligned_hybrid_transfer.projector.matrix atol = 1.0e-10 rtol = 1.0e-10
        @test disk_bond_aligned_hybrid_transfer.diagnostics.transferred_residual_inf < 1.0e-10
    finally
        rm(dir; recursive = true, force = true)
    end
end
