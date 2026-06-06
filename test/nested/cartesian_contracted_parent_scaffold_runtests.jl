@testset "Cartesian contracted parent scaffold" begin
    CP = GaussletBases.CartesianParentGaussletBases
    CCP = GaussletBases.CartesianContractedParents
    axis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 3,
            mapping = IdentityMapping(),
            reference_spacing = 1.0,
        ),
    )
    parent = CP.cartesian_parent_gausslet_basis(axis)
    parent_dim = CP.parent_dimension(parent)

    coefficients = zeros(Float64, parent_dim, 4)
    coefficients[1, 1] = 1.0
    coefficients[2, 2] = 1.0
    coefficients[2, 3] = 2.0
    coefficients[5, 4] = 1.0
    coefficients[6, 4] = -1.0
    unit_a = CCP.CartesianContractionUnit3D(
        :cube_a,
        [1, 2, 3],
        1:2;
        metadata = (shape = :cube,),
    )
    unit_b = CCP.CartesianContractionUnit3D(
        :cube_b,
        [2, 5, 6],
        3:4;
        metadata = (shape = :cube,),
    )
    contracted = CCP.CartesianContractedParent3D(
        parent,
        coefficients;
        units = [unit_a, unit_b],
        metadata = (source = :synthetic,),
    )
    audit = CCP.contracted_parent_structural_audit(contracted)

    @test CCP.contracted_parent_basis(contracted) === parent
    @test CCP.contracted_parent_coefficients(contracted) == coefficients
    @test CCP.contracted_parent_parent_dimension(contracted) == parent_dim
    @test CCP.contracted_parent_dimension(contracted) == 4
    @test CCP.contracted_parent_metadata(contracted).source == :synthetic
    @test CCP.contracted_parent_units(contracted) == [unit_a, unit_b]
    @test CCP.contracted_parent_unit_column_ranges(contracted) == [1:2, 3:4]
    @test CCP.contracted_parent_unit_support_indices(contracted) == [[1, 2, 3], [2, 5, 6]]
    @test CCP.contracted_parent_support_indices(contracted) == [1, 2, 3, 2, 5, 6]
    @test CCP.contraction_unit_role(unit_a) == :cube_a
    @test CCP.contraction_unit_support_indices(unit_a) == [1, 2, 3]
    @test CCP.contraction_unit_column_range(unit_a) == 1:2
    @test CCP.contraction_unit_metadata(unit_a).shape == :cube

    @test coefficients[:, 3] == 2.0 .* coefficients[:, 2]
    @test audit.parent_dimension == parent_dim
    @test audit.contracted_dimension == 4
    @test audit.unit_count == 2
    @test audit.support_entry_count == 6
    @test audit.unique_support_count == 5
    @test audit.duplicate_support_count == 1
    @test audit.missing_support_count == parent_dim - 5
    @test audit.outside_support_count == 0
    @test !audit.support_complete
    @test audit.column_entry_count == 4
    @test audit.unique_column_count == 4
    @test audit.duplicate_column_count == 0
    @test audit.missing_column_count == 0
    @test audit.outside_column_count == 0
    @test audit.column_ranges_cover_contract
    @test audit.structural_ok
    @test !hasproperty(contracted, :gausslet_backend)
    @test !hasproperty(contracted, :backend)
    @test !hasproperty(contracted, :overlap)
    @test !hasproperty(contracted, :interaction_matrix)

    sparse_coefficients = sparse(coefficients)
    sparse_contracted = CCP.CartesianContractedParent3D(
        parent,
        sparse_coefficients;
        units = [unit_a, unit_b],
    )
    @test CCP.contracted_parent_coefficients(sparse_contracted) isa SparseMatrixCSC{Float64,Int}
    @test CCP.contracted_parent_coefficients(sparse_contracted) == sparse_coefficients

    outside_unit = CCP.CartesianContractionUnit3D(:outside, [1, parent_dim + 1], 1:1)
    outside = CCP.CartesianContractedParent3D(
        parent,
        coefficients[:, 1:1];
        units = [outside_unit],
    )
    outside_audit = CCP.contracted_parent_structural_audit(outside)
    @test outside_audit.outside_support_count == 1
    @test !outside_audit.support_complete
    @test !outside_audit.structural_ok

    overlapping_columns = CCP.CartesianContractedParent3D(
        parent,
        coefficients;
        units = [
            CCP.CartesianContractionUnit3D(:left, [1], 1:2),
            CCP.CartesianContractionUnit3D(:right, [2], 2:4),
        ],
    )
    overlapping_column_audit = CCP.contracted_parent_structural_audit(overlapping_columns)
    @test overlapping_column_audit.duplicate_column_count == 1
    @test !overlapping_column_audit.column_ranges_cover_contract
    @test !overlapping_column_audit.structural_ok

    missing_columns = CCP.CartesianContractedParent3D(
        parent,
        coefficients;
        units = [
            CCP.CartesianContractionUnit3D(:first, [1], 1:1),
            CCP.CartesianContractionUnit3D(:last, [2], 3:4),
        ],
    )
    missing_column_audit = CCP.contracted_parent_structural_audit(missing_columns)
    @test missing_column_audit.missing_column_count == 1
    @test !missing_column_audit.column_ranges_cover_contract
    @test !missing_column_audit.structural_ok

    @test_throws ArgumentError CCP.CartesianContractionUnit3D(:empty, [1], 1:0)
    @test_throws ArgumentError CCP.CartesianContractedParent3D(
        parent,
        coefficients;
        units = [CCP.CartesianContractionUnit3D(:bad_columns, [1], 4:5)],
    )
    @test_throws DimensionMismatch CCP.CartesianContractedParent3D(
        parent,
        coefficients[1:(end - 1), :],
    )
end
