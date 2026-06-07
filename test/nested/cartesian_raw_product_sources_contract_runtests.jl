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
