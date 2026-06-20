#!/usr/bin/env julia

# Developer smoke for the active H2 PQS terminal-topology blocked route.
# It reuses the existing harness once, then asserts compact manager-review facts.

fixture = joinpath(
    dirname(@__DIR__),
    "test/driver_inputs/h2_pqs_q5_independent_source_box_r4_supplement_preflight.jl",
)

empty!(ARGS)
append!(ARGS, [
    fixture,
    "save_artifact=false",
    "save_tsv=false",
    "save_basis_artifact=false",
    "save_ham_artifact=false",
    "materialize_route=false",
])

harness_elapsed = @elapsed include(joinpath(@__DIR__, "cartesian_driver_harness.jl"))

function _get(obj, key::Symbol, default = nothing)
    isnothing(obj) && return default
    hasproperty(obj, key) && return getproperty(obj, key)
    obj isa NamedTuple && haskey(obj, key) && return getfield(obj, key)
    return default
end

function _check(label, observed, expected)
    observed == expected || error("$label expected $expected, got $observed")
    println("  ", label, " = ", observed)
    return nothing
end

function _topology_facts(shells)
    scaffold = shells.shellification_scaffold
    coverage = scaffold.coverage
    return (;
        roles = scaffold.ordered_region_roles,
        support_counts = Tuple(region.support_count for region in scaffold.regions),
        duplicate_count = coverage.duplicate_count,
        missing_count = coverage.missing_count,
        outside_count = coverage.outside_count,
        q = scaffold.q,
        core_side = scaffold.core_side,
    )
end

println("h2_pqs_terminal_stage_smoke_checks")
facts = _topology_facts(shells)
_check(
    "terminal_region_roles",
    facts.roles,
    (:atom_contact_core, :shared_molecular_shell, :shared_molecular_shell),
)
_check("support_counts", facts.support_counts, (275, 362, 578))
_check("coverage_duplicate_count", facts.duplicate_count, 0)
_check("coverage_missing_count", facts.missing_count, 0)
_check("coverage_outside_count", facts.outside_count, 0)
_check("q", facts.q, 5)
_check("core_side", facts.core_side, 5)

target = report.physical_gausslet_target_summary
source_plan = report.physical_gausslet_source_plan_summary
preflight = _get(source_plan, :terminal_source_realization_preflight_summary)

_check("target_status", _get(target, :status),
    :blocked_physical_gausslet_target_inventory)
_check("target_retained_transform_authority",
    _get(target, :retained_transform_authority),
    :terminal_retained_rule_preflight)
_check("source_plan_status", _get(source_plan, :status),
    :blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan)
_check("source_plan_blocker", _get(source_plan, :blocker),
    :missing_terminal_shell_projection)
_check("source_plan_materialized", _get(source_plan, :source_plan_materialized),
    false)
_check("source_coefficients_materialized",
    _get(source_plan, :source_coefficients_materialized), false)
_check("source_plan_descriptor_available",
    _get(source_plan, :source_plan_descriptor_available), false)
_check("shared_shell_realization_materialized",
    _get(source_plan, :shared_shell_realization_materialized), false)
_check("terminal_realization_status", _get(preflight, :status),
    :blocked_terminal_source_realization_preflight)
_check("terminal_realization_blocker", _get(preflight, :blocker),
    :missing_terminal_shell_projection)
_check("terminal_realization_total_retained_dimension",
    _get(_get(preflight, :dimension_summary), :total_retained_dimension),
    471,
)
println("h2_pqs_terminal_stage_smoke_elapsed_s=", harness_elapsed)
