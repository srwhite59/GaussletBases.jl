include("cartesian_white_lindsey_adapter_fixture_helpers.jl")

# Manual timing probe for the White--Lindsey boundary-stratum adapter path.
#
# This file is intentionally not included by default test runners. It measures
# phase costs for the local LW adapter fixture so later oracle-validation work
# can avoid duplicating expensive setup.

function _lw_timed(label::AbstractString, thunk)
    result_ref = Ref{Any}()
    elapsed = @elapsed result_ref[] = thunk()
    println("lw_timing.", label, "_s=", elapsed)
    return result_ref[]
end

facet_unit, edge_unit, corner_unit =
    _lw_timed("descriptor_unit_setup", () ->
        _lw_adapter_descriptor_units(; prefix = "lw_timing"))

facet_descriptor, edge_descriptor, corner_descriptor =
    _lw_timed("unit_descriptor_metadata", () -> (
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            facet_unit,
        ),
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            edge_unit,
        ),
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            corner_unit,
        ),
    ))

_ = _lw_timed("blocked_context_metadata", () -> (
    CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficient_context(
        facet_descriptor,
    ),
    CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficient_context(
        edge_descriptor,
    ),
))

corner_coefficients =
    _lw_timed("corner_unit_coefficients", () ->
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficients(
            corner_descriptor,
        ))

doside_source_1d = _lw_timed("doside_source_setup", _lw_adapter_doside_source_1d)

real_units = _lw_timed("real_unit_descriptor_setup", () ->
    _lw_adapter_real_units(doside_source_1d; prefix = "lw_timing"))

real_facet_coefficients =
    _lw_timed("facet_unit_coefficients", () ->
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficients(
            real_units.real_facet_descriptor,
        ))

real_edge_coefficients =
    _lw_timed("edge_unit_coefficients", () ->
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficients(
            real_units.real_edge_descriptor,
        ))

real_pair_coefficients =
    _lw_timed("pair_coefficients_facet_edge", () ->
        CPBMForLWAdapter.white_lindsey_boundary_stratum_pair_unit_coefficients(
            _lw_adapter_unit_pair(
                real_units.real_facet_unit,
                real_units.real_edge_unit,
                1,
            ),
        ))

edge_corner_coefficients =
    _lw_timed("pair_coefficients_edge_corner", () ->
        CPBMForLWAdapter.white_lindsey_boundary_stratum_pair_unit_coefficients(
            _lw_adapter_unit_pair(real_units.real_edge_unit, corner_unit, 2),
        ))

factors = _lw_timed("one_body_factor_setup", _lw_adapter_one_body_factors)

overlap_result = _lw_timed("overlap_block", () ->
    CPBMForLWAdapter.white_lindsey_boundary_stratum_overlap_block(
        real_pair_coefficients;
        parent_axis_counts = (7, 7, 7),
        overlap_1d = factors.overlap_1d,
    ))

position_results = _lw_timed("position_blocks_xyz", () ->
    Tuple(
        CPBMForLWAdapter.white_lindsey_boundary_stratum_position_block(
            real_pair_coefficients;
            axis,
            parent_axis_counts = (7, 7, 7),
            overlap_1d = factors.overlap_1d,
            position_1d = factors.position_1d,
        )
        for axis in (:x, :y, :z)
    ))

x2_results = _lw_timed("x2_blocks_xyz", () ->
    Tuple(
        CPBMForLWAdapter.white_lindsey_boundary_stratum_x2_block(
            real_pair_coefficients;
            axis,
            parent_axis_counts = (7, 7, 7),
            overlap_1d = factors.overlap_1d,
            x2_1d = factors.x2_1d,
        )
        for axis in (:x, :y, :z)
    ))

kinetic_result = _lw_timed("kinetic_block", () ->
    CPBMForLWAdapter.white_lindsey_boundary_stratum_kinetic_block(
        real_pair_coefficients;
        parent_axis_counts = (7, 7, 7),
        overlap_1d = factors.overlap_1d,
        kinetic_1d = factors.kinetic_1d,
    ))

selector_results = _lw_timed("record_selector_all_terms", () ->
    Tuple(
        CPBMForLWAdapter.white_lindsey_boundary_stratum_one_body_block(
            real_pair_coefficients,
            term;
            parent_axis_counts = (7, 7, 7),
            overlap_1d = factors.overlap_1d,
            position_1d = factors.position_1d,
            x2_1d = factors.x2_1d,
            kinetic_1d = factors.kinetic_1d,
        )
        for term in (
            :overlap,
            :position_x,
            :position_y,
            :position_z,
            :x2_x,
            :x2_y,
            :x2_z,
            :kinetic,
        )
    ))

batch_overlap = _lw_timed("batch_selector_overlap", () ->
    CPBMForLWAdapter.white_lindsey_boundary_stratum_one_body_blocks(
        (real_pair_coefficients,),
        :overlap;
        parent_axis_counts = (7, 7, 7),
        overlap_1d = factors.overlap_1d,
    ))

seed_oracle_summary = _lw_timed("seed_oracle_summary", () ->
    CPBMForLWAdapter.white_lindsey_materialized_seed_oracle_summary())

println("lw_timing.summary_facet_coeff_shape=", size(real_facet_coefficients.coefficient_matrix))
println("lw_timing.summary_edge_coeff_shape=", size(real_edge_coefficients.coefficient_matrix))
println("lw_timing.summary_corner_coeff_shape=", size(corner_coefficients.coefficient_matrix))
println("lw_timing.summary_facet_edge_pair_shape=", size(overlap_result.block))
println("lw_timing.summary_edge_corner_ready=", edge_corner_coefficients.pair_unit_coefficient_maps_materialized)
println("lw_timing.summary_position_terms=", Tuple(result.term for result in position_results))
println("lw_timing.summary_x2_terms=", Tuple(result.term for result in x2_results))
println("lw_timing.summary_kinetic_term=", kinetic_result.term)
println("lw_timing.summary_selector_result_count=", length(selector_results))
println("lw_timing.summary_batch_counts=", (batch_overlap.materialized_count, batch_overlap.skipped_count))
println("lw_timing.summary_seed_retained_dimension=", seed_oracle_summary.retained_dimension)
