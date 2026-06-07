# Integration/slow test. Do not include in default nested runner.

@testset "Cartesian contracted parent metric packet" begin
    CP = GaussletBases.CartesianParentGaussletBases
    CCP = GaussletBases.CartesianContractedParents
    CCPM = GaussletBases.CartesianContractedParentMetrics
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
    coefficients = zeros(Float64, parent_dim, 3)
    coefficients[1, 1] = 1.0
    coefficients[2, 1] = 0.25
    coefficients[5, 2] = 1.0
    coefficients[14, 2] = -0.5
    coefficients[parent_dim, 3] = 1.0
    units = [
        CCP.CartesianContractionUnit3D(:left, [1, 2], 1:1),
        CCP.CartesianContractionUnit3D(:middle, [5, 14], 2:2),
        CCP.CartesianContractionUnit3D(:right, [parent_dim], 3:3),
    ]
    contracted = CCP.CartesianContractedParent3D(parent, coefficients; units)
    packet = CCPM.cartesian_contracted_parent_metric_packet(contracted)
    reference = CCPM.cartesian_contracted_parent_metric_packet_dense_reference(contracted)

    @test packet isa CCPM.CartesianContractedParentMetricPacket3D
    @test CCPM.contracted_parent_metric_packet_parent(packet) === contracted
    @test packet.diagnostics.construction_path == :support_local_product
    @test packet.diagnostics.dense_parent_matrix_used == false
    @test reference.diagnostics.construction_path == :dense_reference_oracle
    @test reference.diagnostics.dense_parent_matrix_used == true
    @test size(packet.overlap) == (3, 3)
    @test size(packet.centers) == (3, 3)
    @test length(packet.weights) == 3
    @test packet.overlap ≈ reference.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test packet.position_x ≈ reference.position_x atol = 1.0e-12 rtol = 1.0e-12
    @test packet.position_y ≈ reference.position_y atol = 1.0e-12 rtol = 1.0e-12
    @test packet.position_z ≈ reference.position_z atol = 1.0e-12 rtol = 1.0e-12
    @test packet.weights ≈ reference.weights atol = 1.0e-12 rtol = 1.0e-12
    @test packet.first_moments ≈ reference.first_moments atol = 1.0e-12 rtol = 1.0e-12
    @test packet.centers ≈ reference.centers atol = 1.0e-12 rtol = 1.0e-12
    @test isfinite(packet.diagnostics.overlap_symmetry_error)
    @test isfinite(packet.diagnostics.overlap_identity_error)

    sparse_coefficients = sparse([1, 7, 13, 27], [1, 1, 2, 2], [0.5, -0.25, 1.0, 0.125], parent_dim, 2)
    sparse_contracted = CCP.CartesianContractedParent3D(parent, sparse_coefficients)
    sparse_packet = CCPM.cartesian_contracted_parent_metric_packet(sparse_contracted)
    sparse_reference = CCPM.cartesian_contracted_parent_metric_packet_dense_reference(sparse_contracted)
    @test CCP.contracted_parent_coefficients(sparse_contracted) isa SparseMatrixCSC{Float64,Int}
    @test sparse_packet.diagnostics.dense_parent_matrix_used == false
    @test sparse_packet.diagnostics.coefficient_storage == :SparseMatrixCSC
    @test sparse_packet.overlap ≈ sparse_reference.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test sparse_packet.weights ≈ sparse_reference.weights atol = 1.0e-12 rtol = 1.0e-12

    larger_axis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 5,
            mapping = IdentityMapping(),
            reference_spacing = 1.0,
        ),
    )
    larger_parent = CP.cartesian_parent_gausslet_basis(larger_axis)
    larger_coefficients = sparse([1, 125], [1, 2], [1.0, 1.0], 125, 2)
    larger_contracted = CCP.CartesianContractedParent3D(larger_parent, larger_coefficients)
    larger_packet = CCPM.cartesian_contracted_parent_metric_packet(larger_contracted)
    @test larger_packet.diagnostics.dense_parent_matrix_used == false
    @test_throws ArgumentError CCPM.cartesian_contracted_parent_metric_packet_dense_reference(
        larger_contracted;
        max_parent_dimension = 64,
    )

    distorted_axis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 3,
            mapping = AsinhMapping(a = 0.25, s = 0.5, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    distorted_parent = CP.cartesian_parent_gausslet_basis(distorted_axis)
    distorted_coefficients = sparse([1, 27], [1, 2], [1.0, 1.0], 27, 2)
    distorted_contracted = CCP.CartesianContractedParent3D(
        distorted_parent,
        distorted_coefficients,
    )
    explicit_axis_metric = (
        overlap = Matrix{Float64}(I, 3, 3),
        position = Matrix(Diagonal(centers(distorted_axis))),
        weights = ones(Float64, 3),
        centers = Float64.(centers(distorted_axis)),
        source = :test_explicit_no_quadrature,
    )
    explicit_packet = CCPM.cartesian_contracted_parent_metric_packet(
        distorted_contracted;
        axis_metrics = (
            x = explicit_axis_metric,
            y = explicit_axis_metric,
            z = explicit_axis_metric,
        ),
    )
    @test_throws ArgumentError CCPM.cartesian_contracted_parent_metric_packet(distorted_contracted)
    @test explicit_packet.diagnostics.dense_parent_matrix_used == false
    @test explicit_packet.diagnostics.axis_metric_sources.x == :test_explicit_no_quadrature
    @test_throws ArgumentError CCPM.cartesian_contracted_parent_metric_packet(
        sparse_contracted;
        construction_path = :product_staged_metric_contraction,
    )
end
