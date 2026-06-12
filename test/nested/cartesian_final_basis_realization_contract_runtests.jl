using Test
using GaussletBases
using LinearAlgebra: I

const CFBR = GaussletBases.CartesianFinalBasisRealization
const CPBForFinalBasis = GaussletBases.CartesianCPB
const CRPSForFinalBasis = GaussletBases.CartesianRawProductSources

@testset "CartesianFinalBasisRealization PQS final-basis contract" begin
    source_box = CPBForFinalBasis.filled_cpb(
        1:5,
        1:5,
        1:5;
        role = :pqs_final_basis_contract_source,
    )
    raw_plan = CRPSForFinalBasis.raw_product_box_plan(
        source_box;
        source_mode_dims = (5, 5, 5),
        source_key = :pqs_final_basis_contract_source,
    )
    retained_rule =
        CRPSForFinalBasis.pqs_boundary_product_mode_retained_rule(raw_plan)
    boundary_count = retained_rule.retained_count
    identity = Matrix{Float64}(I, boundary_count, boundary_count)

    @test raw_plan.source_mode_dims == (5, 5, 5)
    @test raw_plan.source_mode_count == 125
    @test boundary_count == 98

    final_basis = CFBR.pqs_source_shell_realization_final_basis(
        raw_plan,
        retained_rule;
        shell_support_indices = collect(1:boundary_count),
        shell_overlap = identity,
        shell_projection = identity,
        lowdin_cleanup = identity,
    )

    @test final_basis.object_kind == :pqs_source_shell_realization_final_basis
    @test final_basis.status == :available_pqs_shell_realization_final_basis
    @test final_basis.blocker === nothing
    @test final_basis.boundary_source_mode_count == 98
    @test final_basis.final_retained_count == 98
    @test final_basis.final_overlap == identity
    @test final_basis.final_overlap_identity_error == 0.0
    @test final_basis.final_overlap_is_identity
    @test !final_basis.ida_data_materialized
    @test !final_basis.rhf_materialized
    @test !final_basis.driver_route_materialized
    @test !final_basis.artifacts_materialized

    shell_operator = copy(identity)
    shell_operator[1, 1] = 2.0
    shell_operator[2, 2] = 3.0
    shell_operator[1, 2] = 0.25
    shell_operator[2, 1] = 0.25
    projected_shell_operator =
        CFBR.pqs_source_shell_projected_one_body_matrix(
            final_basis,
            shell_operator;
            term = :overlap,
        )

    @test projected_shell_operator.status ==
          :materialized_pqs_shell_projected_one_body_matrix
    @test projected_shell_operator.term == :overlap
    @test projected_shell_operator.boundary_operator == shell_operator
    @test projected_shell_operator.final_operator == shell_operator
    @test projected_shell_operator.lowdin_boundary_crosscheck_error == 0.0
    @test !projected_shell_operator.metadata.shell_support_operator_generated

    boundary_overlap_result = (;
        object_kind = :pqs_retained_source_one_body_matrix,
        status = :materialized_pqs_retained_source_one_body_matrix,
        blocker = nothing,
        term = :retained_source_overlap,
        matrix = identity,
        matrix_space = :retained_pqs_source_modes,
        retained_dimension = boundary_count,
        matrix_materialized = true,
    )
    final_overlap_from_boundary =
        CFBR.pqs_source_shell_final_one_body_from_boundary_matrix(
            final_basis,
            boundary_overlap_result,
        )

    @test final_overlap_from_boundary.status ==
          :materialized_pqs_shell_final_one_body_from_boundary_matrix
    @test final_overlap_from_boundary.term == :overlap
    @test final_overlap_from_boundary.final_operator == identity
    @test final_overlap_from_boundary.retained_boundary_operator_input_used
    @test !final_overlap_from_boundary.raw_source_operator_input_used

    boundary_kinetic_result = merge(
        boundary_overlap_result,
        (term = :retained_source_kinetic, matrix = shell_operator),
    )
    final_kinetic_from_boundary =
        CFBR.pqs_source_shell_final_one_body_from_boundary_matrix(
            final_basis,
            boundary_kinetic_result;
            term = :kinetic,
        )

    @test final_kinetic_from_boundary.term == :kinetic
    @test final_kinetic_from_boundary.boundary_operator == shell_operator
    @test final_kinetic_from_boundary.final_operator == shell_operator
    @test final_kinetic_from_boundary.boundary_operator_symmetry_error == 0.0
    @test final_kinetic_from_boundary.final_operator_symmetry_error == 0.0
end
