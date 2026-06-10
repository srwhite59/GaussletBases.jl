# Recipe-specific lowering-source records.
#
# Lowering is the first recipe-specific step after shellification. It chooses
# which CPBs are used to construct retained functions from a region's owned
# support.

"""
    LoweringSource

Recipe-specific source-domain record that connects a shellification region to
the CPBs used to construct retained functions.
"""
struct LoweringSource
    recipe::Symbol
    owned_region::ShellificationRegion
    source_cpbs::Tuple{Vararg{CoordinateProductBox}}
    metadata::NamedTuple
end

"""
    lowering_source(recipe, owned_region, source_cpbs; metadata = (;))

Construct the first recipe-specific object after shellification.

A `LoweringSource` answers: which CPBs are used to construct retained functions
from this owned support? It is not yet a numerical contraction.
"""
function lowering_source(
    recipe::Symbol,
    owned_region::ShellificationRegion,
    source_cpbs;
    metadata = (;),
)
    cpb_tuple = Tuple(source_cpbs)
    isempty(cpb_tuple) && throw(ArgumentError("a LoweringSource requires at least one source CPB"))
    all(cpb -> cpb isa CoordinateProductBox, cpb_tuple) ||
        throw(ArgumentError("LoweringSource source_cpbs must be CoordinateProductBox objects"))
    return LoweringSource(recipe, owned_region, cpb_tuple, NamedTuple(metadata))
end

"""
    source_cpbs(source_or_unit)

Return the CPBs used as lowering sources for a lowering source or final
retained unit.
"""
source_cpbs(source::LoweringSource) = source.source_cpbs

"""
    lowering_recipe(source_or_unit)

Return the symbolic lowering recipe, such as `:white_lindsey_boundary_strata`
or `:pqs_filled_source_cpb`.
"""
lowering_recipe(source::LoweringSource) = source.recipe
owned_support(source::LoweringSource) = owned_support(source.owned_region)
support_count(source::LoweringSource) = sum(support_count, source.source_cpbs; init = 0)

"""
    white_lindsey_boundary_strata_lowering(region, strata; metadata = (;))

Construct a White--Lindsey lowering source from boundary-stratum CPBs.

The `strata` are usually the 6 facet CPBs, 12 edge CPBs, and 8 corner CPBs of a
complete shell. This is LW-specific; PQS should not use this as its primary
shell construction.
"""
function white_lindsey_boundary_strata_lowering(
    region::ShellificationRegion,
    strata;
    metadata = (;),
)
    strata_tuple = Tuple(strata)
    isempty(strata_tuple) &&
        throw(ArgumentError("White-Lindsey boundary-stratum lowering requires CPB strata"))
    all(cpb -> cpb isa CoordinateProductBox, strata_tuple) ||
        throw(ArgumentError("White-Lindsey boundary-stratum lowering requires CPBs"))
    all(cpb -> CartesianCPB.codimension(cpb) >= 1, strata_tuple) ||
        throw(ArgumentError("White-Lindsey boundary-stratum CPBs should be codimension >= 1"))
    return lowering_source(
        :white_lindsey_boundary_strata,
        region,
        strata_tuple;
        metadata = merge((; lowering_family = :boundary_stratum_cpbs), NamedTuple(metadata)),
    )
end

"""
    pqs_filled_source_lowering(region, source_cpb; metadata = (;))

Construct a PQS lowering source from one filled source CPB.

The source CPB is the filled product box used to build boundary COMX-product
modes. It is not the shell itself. Shell projection and Lowdin cleanup happen
later in `pqs_shell_realization`.
"""
function pqs_filled_source_lowering(
    region::ShellificationRegion,
    source_cpb::CoordinateProductBox;
    metadata = (;),
)
    CartesianCPB.codimension(source_cpb) == 0 ||
        throw(ArgumentError("PQS filled-source lowering requires a codimension-0 source CPB"))
    return lowering_source(
        :pqs_filled_source_cpb,
        region,
        (source_cpb,);
        metadata = merge((; lowering_family = :pqs_filled_source_cpb), NamedTuple(metadata)),
    )
end
