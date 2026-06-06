@testset "Cartesian basis representation for nested fixed blocks" begin
    CP = GaussletBases.CartesianParentGaussletBases
    CCP = GaussletBases.CartesianContractedParents
    CCS = GaussletBases.CartesianCarriedSpaces
    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 13,
            mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    fixed_block = one_center_atomic_full_parent_fixed_block(
        basis;
        expansion,
        nside = 5,
    )
    fixed_parent = CP.cartesian_parent_gausslet_basis(fixed_block)
    fixed_contracted_parent = CCP.cartesian_contracted_parent(fixed_block)
    fixed_contracted_audit = CCP.contracted_parent_structural_audit(fixed_contracted_parent)
    representation = basis_representation(fixed_block)
    metadata = basis_metadata(representation)

    @test representation isa CartesianBasisRepresentation3D
    @test metadata.basis_kind == :nested_fixed_block
    @test metadata.parent_kind == :cartesian_product_basis
    @test metadata.parent_axis_counts == (13, 13, 13)
    @test metadata.parent_axis_counts == CP.parent_axis_counts(fixed_parent)
    @test metadata.parent_dimension == 13^3
    @test metadata.parent_dimension == CP.parent_dimension(fixed_parent)
    @test metadata.final_dimension == size(fixed_block.coefficient_matrix, 2)
    @test metadata.working_box == (1:13, 1:13, 1:13)
    @test metadata.route_metadata.shell_kind == :shell_sequence
    @test metadata.route_metadata.working_box_profile == :full_parent
    @test metadata.route_metadata.nside == 5
    @test metadata.route_metadata.support_count == length(fixed_block.support_indices)
    @test size(representation.coefficient_matrix) == size(fixed_block.coefficient_matrix)
    @test representation.support_indices == fixed_block.support_indices
    @test length(representation.support_states) == length(fixed_block.support_indices)
    @test size(metadata.basis_centers) == size(fixed_block.fixed_centers)
    @test CP.parent_center(fixed_parent, (1, 1, 1)) == (
        centers(basis)[1],
        centers(basis)[1],
        centers(basis)[1],
    )
    @test !hasproperty(fixed_parent, :gausslet_backend)
    @test !hasproperty(fixed_parent, :backend)
    @test CCP.contracted_parent_basis(fixed_contracted_parent).parent_box ==
        fixed_parent.parent_box
    @test CCP.contracted_parent_parent_dimension(fixed_contracted_parent) ==
        metadata.parent_dimension
    @test CCP.contracted_parent_dimension(fixed_contracted_parent) ==
        size(fixed_block.coefficient_matrix, 2)
    @test CCP.contracted_parent_coefficients(fixed_contracted_parent) isa
        SparseMatrixCSC{Float64,Int}
    @test CCP.contracted_parent_coefficients(fixed_contracted_parent) ==
        Matrix{Float64}(fixed_block.coefficient_matrix)
    @test only(CCP.contracted_parent_units(fixed_contracted_parent)).role ==
        :nested_fixed_block
    @test only(CCP.contracted_parent_unit_column_ranges(fixed_contracted_parent)) ==
        1:size(fixed_block.coefficient_matrix, 2)
    @test fixed_contracted_audit.outside_support_count == 0
    @test fixed_contracted_audit.column_ranges_cover_contract
    @test fixed_contracted_audit.structural_ok
    @test !hasproperty(fixed_contracted_parent, :gausslet_backend)
    @test !hasproperty(fixed_contracted_parent, :backend)
    @test !hasproperty(fixed_contracted_parent, :interaction_matrix)
    @test !isnothing(fixed_block.factorized_cartesian_parent_basis[])
    @test hasproperty(representation.parent_data, :factorized_cartesian_parent_basis)
    @test representation.parent_data.factorized_cartesian_parent_basis ===
          fixed_block.factorized_cartesian_parent_basis[]
    fixed_carried = CCS.cartesian_carried_space(fixed_block)
    @test CCS.carried_space_parent(fixed_carried) isa CP.CartesianParentGaussletBasis3D
    @test CCS.carried_space_contracted_parent(fixed_carried) isa CCP.CartesianContractedParent3D
    @test CCS.carried_space_representation(fixed_carried) isa CartesianBasisRepresentation3D
    @test CCS.carried_space_diagnostics(fixed_carried).parent_dimension ==
        metadata.parent_dimension
    @test CCS.carried_space_diagnostics(fixed_carried).contracted_dimension ==
        metadata.final_dimension
    @test CCS.carried_space_diagnostics(fixed_carried).contracted_dimension_matches_representation
    @test CCS.carried_space_diagnostics(fixed_carried).contracted_parent_dimension_matches_parent
    @test CCS.carried_space_provenance(fixed_carried).input_kind == :nested_fixed_block

    square_basis, _source, square_fixed_block, _diagnostics =
        _axis_aligned_homonuclear_square_lattice_nested_fixture()
    square_representation = basis_representation(square_fixed_block)
    square_metadata = basis_metadata(square_representation)
    @test square_metadata.basis_kind == :nested_fixed_block
    @test square_metadata.parent_axis_counts == (
        length(square_basis.basis_x),
        length(square_basis.basis_y),
        length(square_basis.basis_z),
    )
    @test square_metadata.parent_dimension == prod(square_metadata.parent_axis_counts)
    @test square_metadata.final_dimension == size(square_fixed_block.coefficient_matrix, 2)
    @test size(square_representation.coefficient_matrix) == size(square_fixed_block.coefficient_matrix)
    @test square_metadata.working_box == square_fixed_block.shell.working_box
    @test square_metadata.route_metadata.support_count == length(square_fixed_block.support_indices)
    @test !isnothing(square_fixed_block.factorized_cartesian_parent_basis[])
    @test hasproperty(square_representation.parent_data, :factorized_cartesian_parent_basis)
    @test square_representation.parent_data.factorized_cartesian_parent_basis ===
          square_fixed_block.factorized_cartesian_parent_basis[]
    square_carried = CCS.cartesian_carried_space(square_fixed_block)
    @test CCS.carried_space_parent(square_carried) isa CP.CartesianParentGaussletBasis3D
    @test CCS.carried_space_contracted_parent(square_carried) isa CCP.CartesianContractedParent3D
    @test CCS.carried_space_diagnostics(square_carried).parent_axis_counts ==
        square_metadata.parent_axis_counts
    @test CCS.carried_space_diagnostics(square_carried).contracted_dimension ==
        square_metadata.final_dimension
    @test CCS.carried_space_diagnostics(square_carried).contracted_dimension_matches_representation
    @test CCS.carried_space_provenance(square_carried).input_kind == :nested_fixed_block
end
