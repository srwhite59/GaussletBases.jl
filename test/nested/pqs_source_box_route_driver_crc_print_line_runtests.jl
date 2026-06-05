using Test
using GaussletBases

@testset "CRC typed pair-operator print line handles unavailable inventory" begin
    report = (;
        low_order_route_core_typed_pair_operator_plan_inventory_available =
            false,
        low_order_route_core_typed_pair_operator_plan_inventory_status =
            :not_available,
        low_order_route_core_typed_pair_operator_plan_count = 0,
        low_order_route_core_typed_pair_operator_plan_blocked_count = 0,
        low_order_route_core_typed_pair_operator_plan_materialized = false,
        low_order_route_core_typed_pair_operator_plan_blocker =
            :not_available,
    )

    @test GaussletBases._pqs_source_box_route_driver_crc_typed_pair_operator_plan_print_line(
        report,
    ) ==
          "CRC typed pair-operator inventory: unavailable (:not_available), typed plans 0, blocked 0, materialized no, blocker :not_available"
end
