# Shellification-owned support records.
#
# ShellificationRegion owns disjoint parent support. It does not own the CPBs
# chosen later by LW/PQS lowering, and it does not own COMX/product/doside
# transforms.

struct OwnedSupport
    support_kind::Symbol
    cpbs::Tuple{Vararg{CoordinateProductBox}}
    outer_box::Union{CoordinateProductBox,Nothing}
    inner_exclusion_box::Union{CoordinateProductBox,Nothing}
    metadata::NamedTuple
end

function owned_cpb(cpb::CoordinateProductBox; support_kind::Symbol = :single_cpb_support, metadata = (;))
    return OwnedSupport(support_kind, (cpb,), nothing, nothing, NamedTuple(metadata))
end

function owned_cpb_union(cpbs; support_kind::Symbol = :cpb_union_support, metadata = (;))
    cpb_tuple = Tuple(cpbs)
    all(cpb -> cpb isa CoordinateProductBox, cpb_tuple) ||
        throw(ArgumentError("owned_cpb_union requires CoordinateProductBox entries"))
    isempty(cpb_tuple) && throw(ArgumentError("owned_cpb_union requires at least one CPB"))
    return OwnedSupport(support_kind, cpb_tuple, nothing, nothing, NamedTuple(metadata))
end

function complete_shell_support(
    outer_box::CoordinateProductBox,
    inner_exclusion_box::CoordinateProductBox;
    support_kind::Symbol = :complete_shell_support,
    metadata = (;),
)
    _assert_complete_shell_boxes(outer_box, inner_exclusion_box)
    return OwnedSupport(
        support_kind,
        (),
        outer_box,
        inner_exclusion_box,
        NamedTuple(metadata),
    )
end

function support_count(support::OwnedSupport)
    if !isnothing(support.outer_box)
        inner_count = isnothing(support.inner_exclusion_box) ? 0 : support_count(support.inner_exclusion_box)
        return support_count(support.outer_box) - inner_count
    end
    return sum(support_count, support.cpbs; init = 0)
end

struct ShellificationRegion
    role::Symbol
    owned_support::OwnedSupport
    metadata::NamedTuple
end

function shellification_region(role::Symbol, owned_support::OwnedSupport; metadata = (;))
    return ShellificationRegion(role, owned_support, NamedTuple(metadata))
end

owned_support(region::ShellificationRegion) = region.owned_support
role(region::ShellificationRegion) = region.role
support_count(region::ShellificationRegion) = support_count(region.owned_support)
