# Integration/slow test. Do not include in default nested runner.

@testset "Atomic direct-product He extent change is not an outer-only identity" begin
    fixture = _atomic_direct_product_he_extent_change_contract_fixture()

    @test fixture.source_count == 3
    @test fixture.target_count == 5
    @test fixture.shared_slice == 2:4
    @test fixture.centers_subset
    @test !fixture.weights_subset
    @test !fixture.coefficient_core_match
end
