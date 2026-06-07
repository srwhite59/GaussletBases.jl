# Integration/slow test. Do not include in default nested runner.

@testset "Cartesian basis representation for direct-product QW bases" begin
    CP = GaussletBases.CartesianParentGaussletBases
    CCS = GaussletBases.CartesianCarriedSpaces
    QWCS = GaussletBases.CartesianQWOperatorCarriedSpaces
    basis, operators, _check = _bond_aligned_diatomic_qw_fixture()
    representation = basis_representation(basis)
    metadata = basis_metadata(representation)
    chain_basis, _chain_ops, _chain_diagnostics = _bond_aligned_homonuclear_chain_qw_fixture()
    square_basis, _square_ops, _square_diagnostics, _square_check =
        _axis_aligned_homonuclear_square_lattice_qw_fixture()
    chain_representation = basis_representation(chain_basis)
    square_representation = basis_representation(square_basis)

    @test representation isa CartesianBasisRepresentation3D
    @test metadata.basis_kind == :direct_product
    @test metadata.parent_kind == :cartesian_product_basis
    @test metadata.parent_axis_counts == (length(basis.basis_x), length(basis.basis_y), length(basis.basis_z))
    @test metadata.parent_dimension == prod(metadata.parent_axis_counts)
    @test metadata.final_dimension == prod(metadata.parent_axis_counts)
    @test metadata.axis_sharing == :shared_xy
    @test metadata.route_metadata.basis_family == :bond_aligned_diatomic
    @test metadata.route_metadata.bond_axis == basis.bond_axis
    @test metadata.route_metadata.nuclei == basis.nuclei
    @test representation.contraction_kind == :identity
    @test isnothing(representation.coefficient_matrix)
    @test isnothing(representation.support_indices)
    @test isnothing(representation.support_states)
    @test metadata.basis_labels == representation.parent_labels
    @test metadata.basis_centers == representation.parent_centers
    @test size(metadata.basis_centers, 1) == metadata.final_dimension
    @test size(metadata.basis_centers, 2) == 3
    @test chain_representation.metadata.basis_kind == :direct_product
    @test chain_representation.metadata.route_metadata.basis_family ==
        :bond_aligned_homonuclear_chain
    @test square_representation.metadata.basis_kind == :direct_product
    @test square_representation.metadata.route_metadata.basis_family ==
        :axis_aligned_homonuclear_square_lattice

    carried = CCS.cartesian_carried_space(basis)
    chain_carried = CCS.cartesian_carried_space(chain_basis)
    square_carried = CCS.cartesian_carried_space(square_basis)
    @test CCS.carried_space_parent(carried) isa CP.CartesianParentGaussletBasis3D
    @test isnothing(CCS.carried_space_contracted_parent(carried))
    @test CCS.carried_space_representation(carried) isa CartesianBasisRepresentation3D
    @test CCS.carried_space_diagnostics(carried).parent_axis_counts ==
        metadata.parent_axis_counts
    @test CCS.carried_space_diagnostics(carried).parent_dimension ==
        metadata.parent_dimension
    @test CCS.carried_space_diagnostics(carried).representation_final_dimension ==
        metadata.final_dimension
    @test CCS.carried_space_diagnostics(carried).has_contracted_parent == false
    @test CCS.carried_space_diagnostics(carried).dense_parent_matrix_used == false
    @test CCS.carried_space_diagnostics(carried).heavy_metric_packet_built == false
    @test CCS.carried_space_provenance(carried).input_kind ==
        :bond_aligned_diatomic_qw_basis
    @test CCS.carried_space_diagnostics(chain_carried).parent_axis_counts ==
        chain_representation.metadata.parent_axis_counts
    @test CCS.carried_space_provenance(chain_carried).input_kind ==
        :bond_aligned_homonuclear_chain_qw_basis
    @test CCS.carried_space_diagnostics(square_carried).parent_axis_counts ==
        square_representation.metadata.parent_axis_counts
    @test CCS.carried_space_provenance(square_carried).input_kind ==
        :axis_aligned_homonuclear_square_lattice_qw_basis

    overlap_before = copy(operators.overlap)
    one_body_before = copy(operators.one_body_hamiltonian)
    interaction_before = copy(operators.interaction_matrix)
    operator_sidecar = QWCS.cartesian_qw_operator_carried_space_sidecar(operators)
    operator_carried = QWCS.qw_operator_carried_space(operator_sidecar)
    operator_representation = QWCS.qw_operator_basis_representation(operator_sidecar)
    operator_diagnostics = QWCS.qw_operator_carried_space_diagnostics(operator_sidecar)
    @test QWCS.qw_operator_carried_space_provenance(operator_sidecar).input_kind ==
        :bond_aligned_direct_product_operator
    @test operator_carried isa CCS.CartesianCarriedSpace3D
    @test operator_representation isa CartesianBasisRepresentation3D
    @test operator_diagnostics.operator_dimension == size(operators.overlap, 1)
    @test operator_diagnostics.operator_gausslet_count == operators.gausslet_count
    @test operator_diagnostics.operator_residual_count == operators.residual_count
    @test operator_diagnostics.carried_dimension == metadata.final_dimension
    @test operator_diagnostics.carried_dimension_matches_operator_gausslet_count
    @test operator_diagnostics.operator_representation_matches_operator_dimension
    @test operator_diagnostics.carried_has_contracted_parent == false
    @test operator_diagnostics.carried_has_staged_sidecar == false
    @test operator_diagnostics.dense_parent_matrix_used == false
    @test operator_diagnostics.heavy_metric_packet_built == false
    @test operators.overlap == overlap_before
    @test operators.one_body_hamiltonian == one_body_before
    @test operators.interaction_matrix == interaction_before

    build_source = QWCS.cartesian_qw_operator_build_source(
        basis;
        nuclear_charges = operators.nuclear_charges,
        nuclear_term_storage = :auto,
        interaction_treatment = operators.interaction_treatment,
        gausslet_backend = operators.gausslet_backend,
    )
    build_diagnostics = QWCS.operator_build_source_diagnostics(build_source)
    @test QWCS.operator_build_source_carried_space(build_source) isa
        CCS.CartesianCarriedSpace3D
    @test QWCS.operator_build_source_provenance(build_source).input_kind ==
        :bond_aligned_direct_product_input
    @test build_source.basis_family == :bond_aligned_diatomic
    @test build_source.carried_space_kind == :direct_product
    @test build_source.nuclear_charges == operators.nuclear_charges
    @test build_source.gausslet_backend == operators.gausslet_backend
    @test build_source.interaction_treatment == operators.interaction_treatment
    @test build_source.nuclear_term_storage == operators.nuclear_term_storage
    @test build_diagnostics.carried_dimension == operator_diagnostics.carried_dimension
    @test build_diagnostics.carried_has_contracted_parent ==
        operator_diagnostics.carried_has_contracted_parent
    @test build_diagnostics.carried_has_staged_sidecar ==
        operator_diagnostics.carried_has_staged_sidecar
    @test build_diagnostics.dense_parent_matrix_used == false
    @test build_diagnostics.heavy_metric_packet_built == false
    @test build_diagnostics.operator_built == false

    construction_record =
        QWCS.cartesian_qw_operator_construction_record(build_source, operators)
    record_diagnostics =
        QWCS.qw_operator_construction_record_diagnostics(construction_record)
    @test record_diagnostics.source_sidecar_agree
    @test isempty(record_diagnostics.mismatch_fields)
    @test isempty(record_diagnostics.ambiguous_mismatch_fields)
    @test :operator_input_kind in record_diagnostics.compared_fields
    @test :gausslet_backend in record_diagnostics.compared_fields
    @test :interaction_treatment in record_diagnostics.compared_fields
    @test :nuclear_charges in record_diagnostics.compared_fields
    @test :carried_parent_axis_counts in record_diagnostics.compared_fields
    @test :carried_parent_dimension in record_diagnostics.compared_fields
    @test :carried_representation_basis_kind in record_diagnostics.compared_fields
    @test :carried_representation_parent_kind in record_diagnostics.compared_fields
    @test :carried_representation_final_dimension in record_diagnostics.compared_fields
    @test :carried_axis_sharing in record_diagnostics.compared_fields
    @test :carried_provenance_input_kind in record_diagnostics.compared_fields
    @test :carried_provenance_route_metadata in record_diagnostics.compared_fields
    @test :carried_has_staged_sidecar in record_diagnostics.compared_fields
    @test record_diagnostics.source_basis_family == :bond_aligned_diatomic
    @test record_diagnostics.source_carried_space_kind == :direct_product
    @test record_diagnostics.sidecar_input_kind == :bond_aligned_direct_product_operator
    @test record_diagnostics.source_parent_axis_counts ==
        record_diagnostics.sidecar_parent_axis_counts
    @test record_diagnostics.source_parent_dimension ==
        record_diagnostics.sidecar_parent_dimension
    @test :coefficient_matrix_values in
        record_diagnostics.intentionally_not_compared
    @test :interaction_matrix_values in
        record_diagnostics.intentionally_not_compared
    @test record_diagnostics.numerical_outputs_changed == false
    @test record_diagnostics.dense_parent_matrix_used == false
    @test record_diagnostics.heavy_metric_packet_built == false
    @test record_diagnostics.operator_built == false
    @test QWCS.qw_operator_construction_record_sidecar(construction_record) isa
        QWCS.CartesianQWOperatorCarriedSpaceSidecar
    @test QWCS.qw_operator_construction_record_provenance(construction_record).source ==
        :cartesian_qw_operator_construction_record

    construction_receipt = QWCS.cartesian_qw_operator_construction_receipt(
        basis;
        nuclear_charges = operators.nuclear_charges,
        nuclear_term_storage = operators.nuclear_term_storage,
        interaction_treatment = operators.interaction_treatment,
        gausslet_backend = operators.gausslet_backend,
    )
    receipt_operators =
        QWCS.qw_operator_construction_receipt_operators(construction_receipt)
    receipt_record =
        QWCS.qw_operator_construction_receipt_record(construction_receipt)
    receipt_diagnostics =
        QWCS.qw_operator_construction_receipt_diagnostics(construction_receipt)
    @test QWCS.qw_operator_construction_receipt_source(construction_receipt) isa
        QWCS.CartesianOperatorBuildSource3D
    @test receipt_record isa QWCS.CartesianQWOperatorConstructionRecord3D
    @test receipt_diagnostics.delegated_to_existing_builder
    @test receipt_diagnostics.builder == :ordinary_cartesian_qiu_white_operators
    @test receipt_diagnostics.source_sidecar_agree
    @test isempty(receipt_diagnostics.mismatch_fields)
    @test receipt_diagnostics.operator_built
    @test receipt_diagnostics.new_hamiltonian_kernel_used == false
    @test receipt_diagnostics.dense_parent_matrix_used == false
    @test receipt_diagnostics.heavy_metric_packet_built == false
    @test receipt_diagnostics.numerical_outputs_changed == false
    @test QWCS.qw_operator_construction_receipt_provenance(construction_receipt).source ==
        :cartesian_qw_operator_construction_receipt
    @test receipt_operators.overlap == operators.overlap
    @test receipt_operators.one_body_hamiltonian == operators.one_body_hamiltonian
    @test receipt_operators.interaction_matrix == operators.interaction_matrix
    @test receipt_operators.gausslet_count == operators.gausslet_count
    @test receipt_operators.residual_count == operators.residual_count
    @test receipt_operators.gausslet_backend == operators.gausslet_backend
    @test receipt_operators.interaction_treatment == operators.interaction_treatment
    @test receipt_operators.nuclear_term_storage == operators.nuclear_term_storage

    mismatched_source = QWCS.cartesian_qw_operator_build_source(
        basis;
        nuclear_charges = [
            operators.nuclear_charges[1] + 0.5,
            operators.nuclear_charges[2],
        ],
        nuclear_term_storage = :auto,
        interaction_treatment = operators.interaction_treatment,
        gausslet_backend = operators.gausslet_backend,
    )
    mismatch_record = QWCS.cartesian_qw_operator_construction_record(
        mismatched_source,
        operators;
        throw_on_mismatch = false,
    )
    mismatch_diagnostics =
        QWCS.qw_operator_construction_record_diagnostics(mismatch_record)
    @test !mismatch_diagnostics.source_sidecar_agree
    @test :nuclear_charges in mismatch_diagnostics.mismatch_fields
    @test_throws ArgumentError QWCS.cartesian_qw_operator_construction_record(
        mismatched_source,
        operators,
    )
    mismatched_carried_source = QWCS.CartesianOperatorBuildSource3D(
        chain_carried,
        build_source.basis_family,
        build_source.carried_space_kind,
        build_source.nuclei,
        build_source.nuclear_charges,
        build_source.gausslet_backend,
        build_source.requested_gausslet_backend,
        build_source.interaction_treatment,
        build_source.nuclear_term_storage,
        build_source.requested_nuclear_term_storage,
        build_source.capabilities,
        build_source.diagnostics,
        build_source.provenance,
    )
    carried_mismatch_record = QWCS.cartesian_qw_operator_construction_record(
        mismatched_carried_source,
        operators;
        throw_on_mismatch = false,
    )
    carried_mismatch_diagnostics =
        QWCS.qw_operator_construction_record_diagnostics(carried_mismatch_record)
    @test !carried_mismatch_diagnostics.source_sidecar_agree
    @test !isempty(
        intersect(
            carried_mismatch_diagnostics.mismatch_fields,
            [
                :carried_parent_axis_counts,
                :carried_parent_dimension,
                :carried_representation_final_dimension,
                :carried_provenance_input_kind,
            ],
        ),
    )
    @test_throws ArgumentError QWCS.cartesian_qw_operator_construction_record(
        mismatched_carried_source,
        operators,
    )
end
