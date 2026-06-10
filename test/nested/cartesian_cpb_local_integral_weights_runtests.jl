# Runtime role: CPB-local product integral weight maker test.
#
# This validates provider-level local product weights built from parent-axis
# integral weights. It does not wire these weights into Coulomb pair factors,
# route/global placement, Hamiltonian assembly, IDA/MWG/PQS semantics, exports,
# or artifacts.

using Test
using GaussletBases

const CPBWeightTest = GaussletBases.CartesianCPB
const CPBProviderWeightTest = GaussletBases.CartesianCPBBlockProviders

function _weight_test_axis_bundle(; x = [2.0, 3.0, 5.0], y = [7.0, 11.0, 13.0], z = [17.0, 19.0, 23.0])
    return (;
        x = (; pgdg_intermediate = (; weights = x,)),
        y = (; pgdg_intermediate = (; weights = y,)),
        z = (; pgdg_intermediate = (; weights = z,)),
    )
end

function _expected_local_weights(axis_bundle, cpb)
    intervals = CPBWeightTest.intervals(cpb)
    weights = axis_bundle
    values = Float64[]
    for ix in intervals[1], iy in intervals[2], iz in intervals[3]
        push!(
            values,
            weights.x.pgdg_intermediate.weights[ix] *
            weights.y.pgdg_intermediate.weights[iy] *
            weights.z.pgdg_intermediate.weights[iz],
        )
    end
    return values
end

function _test_weight_nonclaims(summary)
    @test summary.retained_pqs_weights === false
    @test summary.ida_mwg_semantics === false
    @test summary.route_driver_wiring === false
    @test summary.route_global_matrix_materialized === false
    @test summary.global_matrix_materialized === false
    @test summary.hamiltonian_assembly === false
    @test summary.hamiltonian_data_materialized === false
    @test summary.exports_or_artifacts === false
end

@testset "CPB local integral weights" begin
    axis_bundle = _weight_test_axis_bundle()

    filled_cpb = CPBWeightTest.filled_cpb(
        1:2,
        2:3,
        1:2;
        role = :filled_weight_window,
    )
    filled_weights = CPBProviderWeightTest.cpb_local_integral_weights(
        axis_bundle,
        filled_cpb,
    )
    filled_summary = CPBProviderWeightTest.summary(filled_weights)

    @test filled_summary.status === :available_cpb_local_integral_weights
    @test isnothing(filled_summary.blocker)
    @test filled_summary.weight_source === :axis_pgdg_intermediate_weights
    @test filled_summary.weight_kind === :basis_function_integral_weights
    @test filled_summary.squared_self_integral === false
    @test filled_summary.local_shape == (x = 2, y = 2, z = 2)
    @test filled_summary.local_ordering ===
          :parent_compatible_x_slowest_z_fastest
    @test filled_summary.local_weight_count == 8
    @test filled_summary.axis_weight_counts == (x = 3, y = 3, z = 3)
    @test filled_summary.local_weights_available
    @test filled_weights.local_weights == _expected_local_weights(axis_bundle, filled_cpb)
    @test !hasproperty(filled_summary, :local_weights)
    _test_weight_nonclaims(filled_summary)

    face_cpb = CPBWeightTest.slab_cpb(
        1:3,
        2:2,
        1:2;
        role = :face_weight_window,
    )
    face_weights = CPBProviderWeightTest.cpb_local_integral_weights(
        axis_bundle;
        cpb = face_cpb,
    )
    face_summary = CPBProviderWeightTest.summary(face_weights)

    @test face_summary.status === :available_cpb_local_integral_weights
    @test face_summary.local_shape == (x = 3, y = 1, z = 2)
    @test face_summary.local_weight_count == 6
    @test face_weights.local_weights == _expected_local_weights(axis_bundle, face_cpb)
    _test_weight_nonclaims(face_summary)

    edge_cpb = CPBWeightTest.cpb(
        2:2,
        1:1,
        1:3;
        role = :edge_weight_window,
    )
    edge_weights = CPBProviderWeightTest.cpb_local_integral_weights(
        axis_bundle,
        edge_cpb,
    )
    edge_summary = CPBProviderWeightTest.summary(edge_weights)

    @test edge_summary.status === :available_cpb_local_integral_weights
    @test edge_summary.local_shape == (x = 1, y = 1, z = 3)
    @test edge_summary.local_weight_count == 3
    @test edge_weights.local_weights == _expected_local_weights(axis_bundle, edge_cpb)
    _test_weight_nonclaims(edge_summary)

    point_cpb = CPBWeightTest.cpb(
        3:3,
        2:2,
        1:1;
        role = :point_weight_window,
    )
    point_weights = CPBProviderWeightTest.cpb_local_integral_weights(
        axis_bundle,
        point_cpb,
    )
    point_summary = CPBProviderWeightTest.summary(point_weights)

    @test point_summary.status === :available_cpb_local_integral_weights
    @test point_summary.local_shape == (x = 1, y = 1, z = 1)
    @test point_summary.local_weight_count == 1
    @test point_weights.local_weights == _expected_local_weights(axis_bundle, point_cpb)
    _test_weight_nonclaims(point_summary)

    missing_axis_bundle = merge(axis_bundle, (y = (; pgdg_intermediate = (;)),))
    missing_weights = CPBProviderWeightTest.cpb_local_integral_weights(
        missing_axis_bundle,
        filled_cpb,
    )
    missing_summary = CPBProviderWeightTest.summary(missing_weights)

    @test missing_summary.status === :blocked_cpb_local_integral_weights
    @test missing_summary.blocker === :missing_y_axis_integral_weights
    @test isnothing(missing_weights.local_weights)
    @test missing_summary.local_weight_count == 0
    @test missing_summary.axis_weight_counts == (x = 3, y = :unavailable, z = 3)
    _test_weight_nonclaims(missing_summary)

    short_axis_bundle = _weight_test_axis_bundle(; z = [17.0])
    mismatch_weights = CPBProviderWeightTest.cpb_local_integral_weights(
        short_axis_bundle,
        filled_cpb,
    )
    mismatch_summary = CPBProviderWeightTest.summary(mismatch_weights)

    @test mismatch_summary.status === :blocked_cpb_local_integral_weights
    @test mismatch_summary.blocker === :z_axis_integral_weight_length_mismatch
    @test isnothing(mismatch_weights.local_weights)
    @test mismatch_summary.axis_weight_counts == (x = 3, y = 3, z = 1)
    _test_weight_nonclaims(mismatch_summary)
end
