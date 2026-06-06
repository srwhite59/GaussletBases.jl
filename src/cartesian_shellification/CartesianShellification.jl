"""
    CartesianShellification

Internal shellification module for Cartesian route geometry.

Shellification has one job: given parent geometry and nuclei, decide which
parent rows/sites are owned by which terminal route regions. This module is
route-neutral: it does not choose White--Lindsey boundary strata, PQS source
boxes, retained modes, transforms, operator blocks, Hamiltonians, or report
plumbing.

This first module pass is a typed facade over the existing terminal
shellification implementation. It gives later extraction work a stable target
without changing numerical or staged-route behavior.
"""
module CartesianShellification

using ..CartesianRouteCore

const CRC = CartesianRouteCore

export ShellificationPolicy,
       AtomOutwardShellification,
       OneCenterShellification,
       TerminalRegion,
       ShellificationPlan,
       shellify,
       raw_terminal_geometry,
       terminal_regions,
       coverage,
       summary,
       private_summary,
       scaffold,
       raw_plan

abstract type ShellificationPolicy end

"""
    AtomOutwardShellification(; core_side = 5, q = core_side,
                              bond_axis = :auto, audit_coverage = true)

Policy object for atom-outward terminal shellification. It supports the current
one-center and bond-aligned diatomic geometry path.
"""
struct AtomOutwardShellification <: ShellificationPolicy
    core_side::Int
    q::Int
    bond_axis::Symbol
    audit_coverage::Bool

    function AtomOutwardShellification(
        core_side::Integer,
        q::Integer,
        bond_axis::Symbol,
        audit_coverage::Bool,
    )
        core_side > 0 && isodd(core_side) ||
            throw(ArgumentError("core_side must be a positive odd integer"))
        q > 0 || throw(ArgumentError("q must be a positive integer"))
        bond_axis in (:auto, :x, :y, :z) ||
            throw(ArgumentError("bond_axis must be :auto, :x, :y, or :z"))
        return new(Int(core_side), Int(q), bond_axis, audit_coverage)
    end
end

AtomOutwardShellification(;
    core_side::Integer = 5,
    q::Integer = core_side,
    bond_axis::Symbol = :auto,
    audit_coverage::Bool = true,
) = AtomOutwardShellification(core_side, q, bond_axis, audit_coverage)

"""
    OneCenterShellification(; core_side = 5, q = core_side,
                            audit_coverage = true)

Policy object for the one-center specialization of the atom-outward terminal
shellification path.
"""
struct OneCenterShellification <: ShellificationPolicy
    core_side::Int
    q::Int
    audit_coverage::Bool

    function OneCenterShellification(
        core_side::Integer,
        q::Integer,
        audit_coverage::Bool,
    )
        core_side > 0 && isodd(core_side) ||
            throw(ArgumentError("core_side must be a positive odd integer"))
        q > 0 || throw(ArgumentError("q must be a positive integer"))
        return new(Int(core_side), Int(q), audit_coverage)
    end
end

OneCenterShellification(;
    core_side::Integer = 5,
    q::Integer = core_side,
    audit_coverage::Bool = true,
) = OneCenterShellification(core_side, q, audit_coverage)

"""
    TerminalRegion

Typed view of one terminal shellification region. `route_region` is the
corresponding `CartesianRouteCore.ShellificationRegion`; `raw_region` preserves
the existing compatibility metadata during migration.
"""
struct TerminalRegion
    key::Symbol
    order_index::Int
    role::Symbol
    region_kind::Symbol
    owned_support::CRC.OwnedSupport
    route_region::CRC.ShellificationRegion
    raw_region::Any
    metadata::NamedTuple
end

"""
    ShellificationPlan

Typed shellification plan facade. It preserves the existing raw plan for
compatibility while exposing typed terminal regions and compact accessors.
"""
struct ShellificationPlan
    policy::ShellificationPolicy
    raw_plan::Any
    terminal_regions::Tuple{Vararg{TerminalRegion}}
    coverage::Any
    summary::Any
    diagnostics::Any
    metadata::NamedTuple
end

include("terminal_geometry.jl")

function _cpb_for_box(box; role::Symbol, metadata = (;))
    return CRC.cpb(box; role, metadata)
end

function _owned_support_for_region(region)
    metadata = (;
        terminal_region_order_index = region.order_index,
        terminal_region_role = region.role,
        terminal_region_kind = region.region_kind,
    )
    outer = _cpb_for_box(
        region.outer_box;
        role = Symbol(region.role, "_outer_box"),
        metadata,
    )
    if isnothing(region.inner_exclusion_box)
        return CRC.owned_cpb(
            outer;
            support_kind = Symbol(region.region_kind, "_owned_support"),
            metadata,
        )
    end

    inner = _cpb_for_box(
        region.inner_exclusion_box;
        role = Symbol(region.role, "_inner_exclusion_box"),
        metadata,
    )
    return CRC.complete_shell_support(
        outer,
        inner;
        support_kind = :complete_shell_support,
        metadata,
    )
end

function _terminal_region(raw_region)
    owned_support = _owned_support_for_region(raw_region)
    route_region = CRC.shellification_region(
        raw_region.role,
        owned_support;
        metadata = raw_region.metadata,
    )
    return TerminalRegion(
        Symbol("terminal_region_", raw_region.order_index),
        raw_region.order_index,
        raw_region.role,
        raw_region.region_kind,
        owned_support,
        route_region,
        raw_region,
        raw_region.metadata,
    )
end

function _shellification_plan(raw, policy::ShellificationPolicy; metadata = (;))
    terminal_region_tuple = Tuple(_terminal_region(region) for region in raw.regions)
    return ShellificationPlan(
        policy,
        raw,
        terminal_region_tuple,
        raw.coverage,
        private_summary(raw),
        raw.diagnostics,
        NamedTuple(metadata),
    )
end

"""
    shellify(parent_axes, nuclear_positions; policy = AtomOutwardShellification())
    shellify(parent_axes, nuclear_positions, policy)

Build a typed terminal shellification plan. This is geometry only: output
regions own parent support, but no lowering or operator construction has been
applied.
"""
function shellify(
    parent_axes::NTuple{3,<:AbstractVector},
    nuclear_positions;
    policy::ShellificationPolicy = AtomOutwardShellification(),
    metadata = (;),
)
    return shellify(parent_axes, nuclear_positions, policy; metadata)
end

function shellify(
    parent_axes::NTuple{3,<:AbstractVector},
    nuclear_positions,
    policy::AtomOutwardShellification;
    metadata = (;),
)
    raw = raw_terminal_geometry(
        parent_axes,
        nuclear_positions;
        core_side = policy.core_side,
        q = policy.q,
        bond_axis = policy.bond_axis,
        audit_coverage = policy.audit_coverage,
    )
    return _shellification_plan(raw, policy; metadata)
end

function shellify(
    parent_axes::NTuple{3,<:AbstractVector},
    nuclear_positions,
    policy::OneCenterShellification;
    metadata = (;),
)
    raw = raw_terminal_geometry(
        parent_axes,
        nuclear_positions;
        core_side = policy.core_side,
        q = policy.q,
        bond_axis = :auto,
        audit_coverage = policy.audit_coverage,
    )
    return _shellification_plan(raw, policy; metadata)
end

terminal_regions(plan::ShellificationPlan) = plan.terminal_regions
coverage(plan::ShellificationPlan) = plan.coverage
summary(plan::ShellificationPlan) = plan.summary
raw_plan(plan::ShellificationPlan) = plan.raw_plan

private_summary(plan::ShellificationPlan) = summary(plan)

function scaffold(
    plan::ShellificationPlan;
    route_family::Symbol = :white_lindsey_low_order,
)
    return scaffold(raw_plan(plan); route_family)
end

end # module CartesianShellification
