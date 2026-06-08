using GaussletBases

const CPBMForLWAdapter = GaussletBases.CartesianPairBlockMaterialization
const CRUForLWAdapter = GaussletBases.CartesianRetainedUnits
const CRCForLWAdapter = GaussletBases.CartesianRouteCore
const CUPForLWAdapter = GaussletBases.CartesianUnitPairs
const CPBForLWAdapter = GaussletBases.CartesianCPB

function _lw_adapter_retained_unit(
    unit_key::Symbol,
    unit_index::Int,
    source_cpb,
    stratum_kind::Symbol,
    source_cpb_index::Int;
    dimension_status::Symbol = :not_materialized,
    dimension = nothing,
    extra_metadata = (;),
)
    return CRUForLWAdapter.RetainedUnitRecord(
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
        CRCForLWAdapter.owned_cpb(source_cpb),
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

function _lw_adapter_unit_pair(left_unit, right_unit, pair_index::Int)
    return CUPForLWAdapter.UnitPairRecord(
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

function _lw_adapter_descriptor_units(; prefix::AbstractString = "lw_adapter_test")
    facet_source = CPBForLWAdapter.slab_cpb(
        1:1,
        1:3,
        1:3;
        role = Symbol(prefix, "_facet_source_cpb"),
        metadata = (; stratum_kind = :facet_cpb, source_cpb_index = 1),
    )
    edge_source = CPBForLWAdapter.cpb(
        4:4,
        2:2,
        1:3;
        role = Symbol(prefix, "_edge_source_cpb"),
        metadata = (; stratum_kind = :edge_cpb, source_cpb_index = 2),
    )
    corner_source = CPBForLWAdapter.cpb(
        4:4,
        3:3,
        3:3;
        role = Symbol(prefix, "_corner_source_cpb"),
        metadata = (; stratum_kind = :corner_cpb, source_cpb_index = 3),
    )
    return (
        _lw_adapter_retained_unit(
            Symbol(prefix, "_facet_unit"),
            1,
            facet_source,
            :facet_cpb,
            1,
        ),
        _lw_adapter_retained_unit(
            Symbol(prefix, "_edge_unit"),
            2,
            edge_source,
            :edge_cpb,
            2,
        ),
        _lw_adapter_retained_unit(
            Symbol(prefix, "_corner_unit"),
            3,
            corner_source,
            :corner_cpb,
            3,
        ),
    )
end

function _lw_adapter_doside_source_1d()
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

function _lw_adapter_real_units(
    doside_source_1d;
    prefix::AbstractString = "lw_adapter_test",
)
    real_facet_source = CPBForLWAdapter.slab_cpb(
        1:1,
        2:6,
        2:6;
        role = Symbol(prefix, "_real_facet_source_cpb"),
        metadata = (;
            stratum_kind = :facet_cpb,
            source_cpb_index = 5,
            fixed_axes = (:x,),
            sides = (:low,),
        ),
    )
    real_facet_unit = _lw_adapter_retained_unit(
        Symbol(prefix, "_real_facet_unit"),
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
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            real_facet_unit,
        ),
        (;
            retained_count = 3,
            parent_dims = (7, 7, 7),
            doside_source_1d,
        ),
    )

    real_edge_source = CPBForLWAdapter.cpb(
        7:7,
        1:1,
        2:6;
        role = Symbol(prefix, "_real_edge_source_cpb"),
        metadata = (;
            stratum_kind = :edge_cpb,
            source_cpb_index = 4,
            fixed_axes = (:x, :y),
            sides = (:high, :low),
        ),
    )
    real_edge_unit = _lw_adapter_retained_unit(
        Symbol(prefix, "_real_edge_unit"),
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
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            real_edge_unit,
        ),
        (;
            retained_count = 3,
            parent_dims = (7, 7, 7),
            doside_source_1d,
        ),
    )

    return (;
        real_facet_source,
        real_facet_unit,
        real_facet_descriptor,
        real_edge_source,
        real_edge_unit,
        real_edge_descriptor,
    )
end

function _lw_adapter_one_body_factors()
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

function _lw_adapter_prepared_facet_edge_fixture(;
    prefix::AbstractString = "lw_adapter_test",
)
    facet_unit, edge_unit, corner_unit =
        _lw_adapter_descriptor_units(; prefix)
    facet_descriptor =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            facet_unit,
        )
    edge_descriptor =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            edge_unit,
        )
    corner_descriptor =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_adapter_descriptor(
            corner_unit,
        )
    doside_source_1d = _lw_adapter_doside_source_1d()
    real_units = _lw_adapter_real_units(doside_source_1d; prefix)
    real_facet_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficients(
            real_units.real_facet_descriptor,
        )
    real_edge_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_unit_coefficients(
            real_units.real_edge_descriptor,
        )
    real_pair_coefficients =
        CPBMForLWAdapter.white_lindsey_boundary_stratum_pair_unit_coefficients(
            _lw_adapter_unit_pair(
                real_units.real_facet_unit,
                real_units.real_edge_unit,
                1,
            ),
        )
    return (;
        facet_unit,
        edge_unit,
        corner_unit,
        facet_descriptor,
        edge_descriptor,
        corner_descriptor,
        doside_source_1d,
        real_units,
        real_facet_coefficients,
        real_edge_coefficients,
        real_pair_coefficients,
        factors = _lw_adapter_one_body_factors(),
    )
end
