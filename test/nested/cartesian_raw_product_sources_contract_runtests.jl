using Test
using GaussletBases

const CRPS = GaussletBases.CartesianRawProductSources
const CPBForRawSources = GaussletBases.CartesianCPB

@testset "CartesianRawProductSources raw product box plan contract" begin
    source_box = CPBForRawSources.filled_cpb(
        10:11,
        20:22,
        30:33;
        role = :raw_product_source_fixture,
    )

    plan = CRPS.raw_product_box_plan(
        source_box;
        source_key = :fixture_source,
        source_mode_dims = (3, 4, 5),
        metadata = (; fixture = :raw_product_sources_contract),
    )

    @test plan.source_key == :fixture_source
    @test CRPS.source_cpb(plan) === source_box
    @test plan.source_intervals == CPBForRawSources.intervals(source_box)
    @test plan.source_shape == (2, 3, 4)
    @test CRPS.source_mode_dims(plan) == (3, 4, 5)
    @test plan.source_shape != CRPS.source_mode_dims(plan)
    @test CRPS.source_mode_count(plan) == 60

    indices = CRPS.source_mode_indices(plan)
    @test length(indices) == 60
    @test first(indices) == (1, 1, 1)
    @test indices[2] == (1, 1, 2)
    @test indices[5] == (1, 1, 5)
    @test indices[6] == (1, 2, 1)
    @test last(indices) == (3, 4, 5)
    @test plan.source_mode_column_indices == Tuple(1:60)
    @test CRPS.source_mode_indices((3, 4, 5)) == indices

    facts = CRPS.axis_transform_facts(plan)
    @test length(facts) == 3
    @test Tuple(fact.axis for fact in facts) == (1, 2, 3)
    @test Tuple(fact.source_mode_dim for fact in facts) == (3, 4, 5)
    @test Tuple(fact.source_interval for fact in facts) ==
          CPBForRawSources.intervals(source_box)
    @test all(fact -> fact.coefficient_status === :not_materialized, facts)
    @test all(fact -> isnothing(fact.coefficient_matrix), facts)

    compact = CRPS.summary(plan)
    @test compact.object_kind == :raw_product_box_plan_summary
    @test compact.status == :available_raw_product_box_plan
    @test compact.source_key == :fixture_source
    @test compact.source_shape == (2, 3, 4)
    @test compact.source_mode_dims == (3, 4, 5)
    @test compact.source_mode_count == 60
    @test compact.source_mode_ordering == :x_major_y_major_z_fast
    @test compact.axis_transform_statuses ==
          (:not_materialized, :not_materialized, :not_materialized)
    @test compact.materialized == false
    @test compact.retained_rule_materialized == false
    @test compact.shell_realization_materialized == false
    @test compact.pair_blocks_materialized == false
    @test compact.hamiltonian_data_materialized == false
    @test compact.artifacts_materialized == false
end

@testset "CartesianRawProductSources materialized source-axis transforms" begin
    source_box = CPBForRawSources.filled_cpb(
        10:11,
        20:22,
        30:33;
        role = :raw_product_source_axis_transform_fixture,
    )
    transform_x = [1.0 0.1 0.2; 0.3 1.1 0.4]
    transform_y = [
        1.2 0.5 0.6 0.7
        0.8 1.3 0.9 1.0
        1.1 1.2 1.4 1.5
    ]
    transform_z = [
        1.6 0.2 0.3 0.4 0.5
        0.6 1.7 0.7 0.8 0.9
        1.0 1.1 1.8 1.2 1.3
        1.4 1.5 1.6 1.9 1.7
    ]

    plan = CRPS.raw_product_box_plan(
        source_box;
        source_key = :materialized_axis_transform_fixture,
        source_mode_dims = (3, 4, 5),
        axis_transform_matrices = (;
            x = transform_x,
            y = transform_y,
            z = transform_z,
        ),
    )
    facts = CRPS.axis_transform_facts(plan)

    @test Tuple(fact.coefficient_status for fact in facts) ==
          (:materialized, :materialized, :materialized)
    @test facts[1].coefficient_matrix == transform_x
    @test facts[2].coefficient_matrix == transform_y
    @test facts[3].coefficient_matrix == transform_z
    @test Tuple(size(fact.coefficient_matrix) for fact in facts) ==
          ((2, 3), (3, 4), (4, 5))
    @test all(fact -> fact.metadata.materialized, facts)
    @test CRPS.summary(plan).axis_transform_statuses ==
          (:materialized, :materialized, :materialized)
    @test CRPS.summary(plan).materialized

    fact_plan = CRPS.raw_product_box_plan(
        source_box;
        source_key = :materialized_axis_transform_fact_fixture,
        source_mode_dims = (3, 4, 5),
        axis_transform_facts = facts,
    )
    @test CRPS.axis_transform_facts(fact_plan) === facts

    @test_throws ArgumentError CRPS.raw_product_box_plan(
        source_box;
        source_mode_dims = (3, 4, 5),
        axis_transform_matrices = (transform_x, transform_y, zeros(3, 5)),
    )
end

@testset "CartesianRawProductSources PQS boundary retained source-mode rule" begin
    rule = CRPS.pqs_boundary_product_mode_retained_rule(
        (5, 5, 5);
        source_key = :pqs_boundary_fixture,
    )

    @test rule isa CRPS.PQSBoundaryProductModeRetainedRule
    @test rule.source_key == :pqs_boundary_fixture
    @test rule.source_mode_dims == (5, 5, 5)
    @test rule.source_mode_ordering == :x_major_y_major_z_fast
    @test rule.retained_rule_kind == :boundary_comx_product_mode_selection
    @test rule.transform_kind == :source_mode_column_selector
    @test rule.retained_count == 98
    @test !rule.shell_realization_materialized
    @test !rule.lowdin_cleanup_used

    retained_modes = CRPS.retained_mode_indices(rule)
    retained_columns = CRPS.retained_column_indices(rule)
    source_modes = CRPS.source_mode_indices((5, 5, 5))
    expected_columns = [
        column
        for (column, mode) in pairs(source_modes)
        if any(axis -> mode[axis] == 1 || mode[axis] == 5, 1:3)
    ]

    @test retained_columns == expected_columns
    @test retained_modes == collect(source_modes[expected_columns])
    @test first(retained_modes) == (1, 1, 1)
    @test retained_modes[5] == (1, 1, 5)
    @test retained_modes[6] == (1, 2, 1)
    @test last(retained_modes) == (5, 5, 5)
    @test all(mode -> any(axis -> mode[axis] in (1, 5), 1:3), retained_modes)
    @test !any(mode -> all(axis -> 1 < mode[axis] < 5, 1:3), retained_modes)

    compact = CRPS.summary(rule)
    @test compact.object_kind == :pqs_boundary_product_mode_retained_rule_summary
    @test compact.status == :available_pqs_boundary_product_mode_retained_rule
    @test compact.source_key == :pqs_boundary_fixture
    @test compact.source_mode_dims == (5, 5, 5)
    @test compact.retained_rule_kind == :boundary_comx_product_mode_selection
    @test compact.retained_count == 98
    @test compact.transform_kind == :source_mode_column_selector
    @test compact.shell_realization_materialized == false
    @test compact.lowdin_cleanup_used == false
    @test compact.transforms_materialized == false
    @test compact.coefficient_maps_materialized == false
    @test compact.pair_blocks_materialized == false
    @test compact.hamiltonian_data_materialized == false
    @test compact.artifacts_materialized == false
end

@testset "CartesianRawProductSources validation" begin
    source_box = CPBForRawSources.filled_cpb(1:2, 1:2, 1:2)

    @test_throws ArgumentError CRPS.raw_product_box_plan(
        source_box;
        source_mode_dims = (3, 4),
    )
    @test_throws ArgumentError CRPS.raw_product_box_plan(
        source_box;
        source_mode_dims = (3, 0, 5),
    )
    @test_throws ArgumentError CRPS.raw_product_box_plan(
        source_box;
        source_mode_dims = (3, 4.0, 5),
    )
    @test_throws ArgumentError CRPS.raw_product_box_plan(
        source_box;
        source_mode_dims = (3, 4, 5),
        source_mode_ordering = :unsupported_ordering,
    )
    @test_throws ArgumentError CRPS.raw_product_box_plan(
        :not_a_cpb;
        source_mode_dims = (3, 4, 5),
    )

    unavailable = CRPS.unavailable_summary(:not_selected, :not_selected_route)
    @test unavailable.object_kind == :raw_product_box_plan_summary
    @test unavailable.status == :not_selected
    @test unavailable.blocker == :not_selected_route
    @test unavailable.source_mode_count == 0
    @test unavailable.axis_transform_statuses == ()
    @test unavailable.materialized == false
    @test unavailable.retained_rule_materialized == false
    @test unavailable.shell_realization_materialized == false
    @test unavailable.pair_blocks_materialized == false
    @test unavailable.hamiltonian_data_materialized == false
    @test unavailable.artifacts_materialized == false
end

@testset "old raw product source-mode indices adapter" begin
    old_indices = GaussletBases._cartesian_raw_product_box_source_mode_indices((3, 4, 5))
    new_indices = CRPS.source_mode_indices((3, 4, 5))

    @test old_indices isa Vector{NTuple{3,Int}}
    @test old_indices == collect(new_indices)
    @test first(old_indices) == first(new_indices)
    @test old_indices[2] == (1, 1, 2)
    @test old_indices[6] == (1, 2, 1)
    @test last(old_indices) == last(new_indices)

    @test_throws ArgumentError GaussletBases._cartesian_raw_product_box_source_mode_indices(
        (3, 0, 5),
    )
end
