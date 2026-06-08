using GaussletBases

# Manual timing probe for the White--Lindsey boundary-stratum adapter path.
#
# This file is intentionally not included by default test runners. It measures
# phase costs for the local LW adapter fixture so later oracle-validation work
# can avoid duplicating expensive setup.

const CPBMForLWTiming = GaussletBases.CartesianPairBlockMaterialization
const CRUForLWTiming = GaussletBases.CartesianRetainedUnits
const CRCForLWTiming = GaussletBases.CartesianRouteCore
const CUPForLWTiming = GaussletBases.CartesianUnitPairs
const CPBForLWTiming = GaussletBases.CartesianCPB

function _lw_timing_retained_unit(
    unit_key::Symbol,
    unit_index::Int,
    source_cpb,
    stratum_kind::Symbol,
    source_cpb_index::Int;
    dimension_status::Symbol = :not_materialized,
    dimension = nothing,
    extra_metadata = (;),
)
    return CRUForLWTiming.RetainedUnitRecord(
        unit_key,
        unit_index,
        :white_lindsey_boundary_stratum_retained_unit,
        Symbol(unit_key, "_contract"),
        Symbol(unit_key, "_terminal_region"),
        :synthetic_terminal_region,
        :synthetic_terminal_region,
        :white_lindsey_boundary_strata,
        :white_lindsey_boundary_stratum_product,
        :direct_or_trivial_embedding,
        CRCForLWTiming.owned_cpb(source_cpb),
        (source_cpb,),
        source_cpb_index,
        dimension_status,
        dimension,
        :not_materialized,
        nothing,
        nothing,
        false,
        merge((; stratum_kind, source_cpb_index), extra_metadata),
    )
end

function _lw_timing_unit_pair(left_unit, right_unit, pair_index::Int)
    return CUPForLWTiming.UnitPairRecord(
        (left_unit.unit_key, right_unit.unit_key),
        pair_index,
        Symbol(String(left_unit.unit_kind), "__", String(right_unit.unit_kind)),
        left_unit,
        right_unit,
        left_unit.unit_index,
        right_unit.unit_index,
        left_unit.unit_key,
        right_unit.unit_key,
        left_unit.unit_kind,
        right_unit.unit_kind,
        nothing,
        false,
        (;),
    )
end

function _lw_timing_descriptor_units()
    facet_source = CPBForLWTiming.slab_cpb(
        1:1,
        1:3,
        1:3;
        role = :lw_timing_facet_source_cpb,
        metadata = (; stratum_kind = :facet_cpb, source_cpb_index = 1),
    )
    edge_source = CPBForLWTiming.cpb(
        4:4,
        2:2,
        1:3;
        role = :lw_timing_edge_source_cpb,
        metadata = (; stratum_kind = :edge_cpb, source_cpb_index = 2),
    )
    corner_source = CPBForLWTiming.cpb(
        4:4,
        3:3,
        3:3;
        role = :lw_timing_corner_source_cpb,
        metadata = (; stratum_kind = :corner_cpb, source_cpb_index = 3),
    )
    return (
        _lw_timing_retained_unit(
            :lw_timing_facet_unit,
            1,
            facet_source,
            :facet_cpb,
            1,
        ),
        _lw_timing_retained_unit(
            :lw_timing_edge_unit,
            2,
            edge_source,
            :edge_cpb,
            2,
        ),
        _lw_timing_retained_unit(
            :lw_timing_corner_unit,
            3,
            corner_source,
            :corner_cpb,
            3,
        ),
    )
end

function _lw_timing_doside_source_1d()
    count = 7
    endpoint = (count - 1) / 2
    a = 0.25
    xmax = 10.0
    tail_spacing = 10.0
    s = asinh(xmax / a) / (endpoint - xmax / tail_spacing)
    basis = build_basis(MappedUniformBasisSpec(:G10;
        count,
        mapping = AsinhMapping(; a, s, tail_spacing),
        reference_spacing = 1.0,
    ))
    expansion = coulomb_gaussian_expansion(doacc = false)
    return GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
end

function _lw_timing_real_units(doside_source_1d)
    real_facet_source = CPBForLWTiming.slab_cpb(
        1:1,
        2:6,
        2:6;
        role = :lw_timing_real_facet_source_cpb,
        metadata = (;
            stratum_kind = :facet_cpb,
            source_cpb_index = 5,
            fixed_axes = (:x,),
            sides = (:low,),
        ),
    )
    real_facet_unit = _lw_timing_retained_unit(
        :lw_timing_real_facet_unit,
        5,
        real_facet_source,
        :facet_cpb,
        5;
        dimension_status = :available,
        dimension = 3,
        extra_metadata = (;
            parent_dims = (7, 7, 7),
            doside_source_1d,
        ),
    )
    real_facet_descriptor = merge(
        CPBMForLWTiming.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            real_facet_unit,
        ),
        (;
            retained_count = 3,
            parent_dims = (7, 7, 7),
            doside_source_1d,
        ),
    )

    real_edge_source = CPBForLWTiming.cpb(
        7:7,
        1:1,
        2:6;
        role = :lw_timing_real_edge_source_cpb,
        metadata = (;
            stratum_kind = :edge_cpb,
            source_cpb_index = 4,
            fixed_axes = (:x, :y),
            sides = (:high, :low),
        ),
    )
    real_edge_unit = _lw_timing_retained_unit(
        :lw_timing_real_edge_unit,
        4,
        real_edge_source,
        :edge_cpb,
        4;
        dimension_status = :available,
        dimension = 3,
        extra_metadata = (;
            parent_dims = (7, 7, 7),
            doside_source_1d,
        ),
    )
    real_edge_descriptor = merge(
        CPBMForLWTiming.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            real_edge_unit,
        ),
        (;
            retained_count = 3,
            parent_dims = (7, 7, 7),
            doside_source_1d,
        ),
    )

    return (;
        real_facet_unit,
        real_edge_unit,
        real_facet_descriptor,
        real_edge_descriptor,
    )
end

function _lw_timing_one_body_factors()
    return (;
        overlap_1d = (;
            x = ones(Float64, 7, 7),
            y = ones(Float64, 7, 7),
            z = ones(Float64, 7, 7),
        ),
        position_1d = (;
            x = [Float64(i + j) for i in 1:7, j in 1:7],
            y = [Float64(2i + j) for i in 1:7, j in 1:7],
            z = [Float64(i + 2j) for i in 1:7, j in 1:7],
        ),
        x2_1d = (;
            x = [Float64((i + j)^2) for i in 1:7, j in 1:7],
            y = [Float64((2i + j)^2) for i in 1:7, j in 1:7],
            z = [Float64((i + 2j)^2) for i in 1:7, j in 1:7],
        ),
        kinetic_1d = (;
            x = [Float64(i + 3j) for i in 1:7, j in 1:7],
            y = [Float64(3i + j) for i in 1:7, j in 1:7],
            z = [Float64(2i + 2j) for i in 1:7, j in 1:7],
        ),
    )
end

function _lw_timed(label::AbstractString, thunk)
    result_ref = Ref{Any}()
    elapsed = @elapsed result_ref[] = thunk()
    println("lw_timing.", label, "_s=", elapsed)
    return result_ref[]
end

facet_unit, edge_unit, corner_unit =
    _lw_timed("descriptor_unit_setup", _lw_timing_descriptor_units)

facet_descriptor, edge_descriptor, corner_descriptor =
    _lw_timed("unit_descriptor_metadata", () -> (
        CPBMForLWTiming.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            facet_unit,
        ),
        CPBMForLWTiming.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            edge_unit,
        ),
        CPBMForLWTiming.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            corner_unit,
        ),
    ))

_ = _lw_timed("blocked_context_metadata", () -> (
    CPBMForLWTiming.white_lindsey_boundary_stratum_unit_coefficient_context(
        facet_descriptor,
    ),
    CPBMForLWTiming.white_lindsey_boundary_stratum_unit_coefficient_context(
        edge_descriptor,
    ),
))

corner_coefficients =
    _lw_timed("corner_unit_coefficients", () ->
        CPBMForLWTiming.white_lindsey_boundary_stratum_unit_coefficients(
            corner_descriptor,
        ))

doside_source_1d = _lw_timed("doside_source_setup", _lw_timing_doside_source_1d)

real_units = _lw_timed("real_unit_descriptor_setup", () ->
    _lw_timing_real_units(doside_source_1d))

real_facet_coefficients =
    _lw_timed("facet_unit_coefficients", () ->
        CPBMForLWTiming.white_lindsey_boundary_stratum_unit_coefficients(
            real_units.real_facet_descriptor,
        ))

real_edge_coefficients =
    _lw_timed("edge_unit_coefficients", () ->
        CPBMForLWTiming.white_lindsey_boundary_stratum_unit_coefficients(
            real_units.real_edge_descriptor,
        ))

real_pair_coefficients =
    _lw_timed("pair_coefficients_facet_edge", () ->
        CPBMForLWTiming.white_lindsey_boundary_stratum_pair_unit_coefficients(
            _lw_timing_unit_pair(
                real_units.real_facet_unit,
                real_units.real_edge_unit,
                1,
            ),
        ))

edge_corner_coefficients =
    _lw_timed("pair_coefficients_edge_corner", () ->
        CPBMForLWTiming.white_lindsey_boundary_stratum_pair_unit_coefficients(
            _lw_timing_unit_pair(real_units.real_edge_unit, corner_unit, 2),
        ))

factors = _lw_timed("one_body_factor_setup", _lw_timing_one_body_factors)

overlap_result = _lw_timed("overlap_block", () ->
    CPBMForLWTiming.white_lindsey_boundary_stratum_overlap_block(
        real_pair_coefficients;
        parent_axis_counts = (7, 7, 7),
        overlap_1d = factors.overlap_1d,
    ))

position_results = _lw_timed("position_blocks_xyz", () ->
    Tuple(
        CPBMForLWTiming.white_lindsey_boundary_stratum_position_block(
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
        CPBMForLWTiming.white_lindsey_boundary_stratum_x2_block(
            real_pair_coefficients;
            axis,
            parent_axis_counts = (7, 7, 7),
            overlap_1d = factors.overlap_1d,
            x2_1d = factors.x2_1d,
        )
        for axis in (:x, :y, :z)
    ))

kinetic_result = _lw_timed("kinetic_block", () ->
    CPBMForLWTiming.white_lindsey_boundary_stratum_kinetic_block(
        real_pair_coefficients;
        parent_axis_counts = (7, 7, 7),
        overlap_1d = factors.overlap_1d,
        kinetic_1d = factors.kinetic_1d,
    ))

selector_results = _lw_timed("record_selector_all_terms", () ->
    Tuple(
        CPBMForLWTiming.white_lindsey_boundary_stratum_one_body_block(
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
    CPBMForLWTiming.white_lindsey_boundary_stratum_one_body_blocks(
        (real_pair_coefficients,),
        :overlap;
        parent_axis_counts = (7, 7, 7),
        overlap_1d = factors.overlap_1d,
    ))

seed_oracle_summary = _lw_timed("seed_oracle_summary", () ->
    CPBMForLWTiming.white_lindsey_materialized_seed_oracle_summary())

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
