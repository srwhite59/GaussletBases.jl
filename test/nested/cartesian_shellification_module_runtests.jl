using Test
using GaussletBases

const CSH = GaussletBases.CartesianShellification
const CRCForShellification = GaussletBases.CartesianRouteCore

function _shellification_module_compact_fingerprint(plan)
    raw = CSH.raw_plan(plan)
    return (;
        plan_type = typeof(plan),
        raw_region_count = raw.region_count,
        typed_region_count = length(CSH.terminal_regions(plan)),
        coverage_complete = CSH.coverage(plan).coverage_complete,
        summary_region_count = CSH.summary(plan).region_count,
        region_kinds = Tuple(region.region_kind for region in CSH.terminal_regions(plan)),
        support_counts =
            Tuple(CRCForShellification.support_count(region.owned_support)
                  for region in CSH.terminal_regions(plan)),
        route_region_roles = Tuple(
            CRCForShellification.role(region.route_region)
            for region in CSH.terminal_regions(plan)
        ),
        raw_region_roles = raw.region_roles,
    )
end

@testset "CartesianShellification module facade" begin
    parent_axes = (collect(1:9), collect(1:9), collect(1:13))
    policy = CSH.AtomOutwardShellification(;
        core_side = 3,
        q = 3,
        bond_axis = :z,
    )
    plan = CSH.shellify(parent_axes, ((5, 5, 4), (5, 5, 10)), policy)
    raw = CSH.raw_plan(plan)
    fingerprint = _shellification_module_compact_fingerprint(plan)

    @test plan isa CSH.ShellificationPlan
    @test plan.policy === policy
    @test fingerprint.raw_region_count == raw.region_count
    @test fingerprint.typed_region_count == raw.region_count
    @test fingerprint.summary_region_count == raw.region_count
    @test fingerprint.coverage_complete
    @test fingerprint.raw_region_roles == fingerprint.route_region_roles
    @test fingerprint.support_counts ==
          Tuple(region.support_count for region in raw.regions)
    @test all(region -> region isa CSH.TerminalRegion, CSH.terminal_regions(plan))
    @test all(
        region -> region.owned_support isa CRCForShellification.OwnedSupport,
        CSH.terminal_regions(plan),
    )
    @test all(
        region -> region.route_region isa CRCForShellification.ShellificationRegion,
        CSH.terminal_regions(plan),
    )
    @test !raw.aggregate_atom_boxes_emitted
    @test raw.diagnostics.shellification_geometry_only
    @test !raw.diagnostics.coordinate_product_box_lowering_materialized

    one_center = CSH.shellify(
        parent_axes,
        (5, 5, 7);
        policy = CSH.OneCenterShellification(; core_side = 3, q = 3),
    )
    @test one_center.policy isa CSH.OneCenterShellification
    @test CSH.raw_plan(one_center).system_kind == :one_center
    @test CSH.coverage(one_center).coverage_complete
    @test length(CSH.terminal_regions(one_center)) ==
          CSH.raw_plan(one_center).region_count
end
