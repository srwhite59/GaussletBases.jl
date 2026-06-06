@testset "Cartesian basis representation for atomic QW residual bases" begin
    fixture = _atomic_hybrid_cartesian_representation_fixture()
    operators = fixture.full_ops
    representation = fixture.full_rep
    metadata = basis_metadata(representation)
    supplement_representation = representation.parent_data.supplement_representation

    @test representation isa CartesianBasisRepresentation3D
    @test metadata.basis_kind == :hybrid_residual
    @test metadata.parent_kind == :cartesian_plus_supplement_raw
    @test metadata.final_dimension == length(operators.orbital_data)
    @test metadata.final_dimension == size(operators.raw_to_final, 2)
    @test metadata.parent_dimension == size(operators.raw_to_final, 1)
    @test metadata.route_metadata.gausslet_count == operators.gausslet_count
    @test metadata.route_metadata.residual_count == operators.residual_count
    @test metadata.route_metadata.supplement_kind == :atomic_cartesian_shell
    @test metadata.route_metadata.supplement_lmax == fixture.supplement.lmax
    @test size(representation.coefficient_matrix) == size(operators.raw_to_final)
    @test length(representation.parent_labels) == size(operators.raw_to_final, 1)
    @test size(representation.parent_centers, 1) == size(operators.raw_to_final, 1)
    @test hasproperty(representation.parent_data, :cartesian_parent_representation)
    @test representation.parent_data.cartesian_parent_representation.metadata.basis_kind ==
        :nested_fixed_block
    @test representation.parent_data.cartesian_parent_representation.metadata.final_dimension ==
        operators.gausslet_count
    @test hasproperty(representation.parent_data, :supplement_representation)
    @test hasproperty(representation.parent_data, :factorized_cartesian_parent_basis)
    @test hasproperty(representation.parent_data, :cartesian_supplement_axis_tables)
    @test supplement_representation isa CartesianGaussianShellSupplementRepresentation3D
    @test supplement_representation.supplement_kind == :atomic_cartesian_shell
    @test length(supplement_representation.orbitals) ==
        size(operators.raw_to_final, 1) - operators.gausslet_count
    @test size(representation.parent_data.cartesian_supplement_axis_tables.x, 2) ==
        length(supplement_representation.orbitals)
    @test size(representation.parent_data.cartesian_supplement_axis_tables.y, 2) ==
        length(supplement_representation.orbitals)
    @test size(representation.parent_data.cartesian_supplement_axis_tables.z, 2) ==
        length(supplement_representation.orbitals)
    @test any(
        orbital -> sum(orbital.angular_powers) > 0,
        supplement_representation.orbitals,
    )
end
