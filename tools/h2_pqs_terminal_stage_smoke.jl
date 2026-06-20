#!/usr/bin/env julia

# Developer smoke for the active H2 PQS terminal-topology/materialization path.
# It reuses the existing harness once, then asserts compact manager-review facts.

fixture = joinpath(
    dirname(@__DIR__),
    "test/driver_inputs/h2_pqs_q5_independent_source_box_r4_supplement_materialized.jl",
)

empty!(ARGS)
append!(ARGS, [
    fixture,
    "save_artifact=false",
    "save_tsv=false",
    "save_basis_artifact=false",
    "save_ham_artifact=false",
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

function _check_approx(label, observed, expected; atol)
    observed isa Number || error("$label expected numeric value, got $observed")
    isapprox(observed, expected; atol, rtol = 0.0) ||
        error("$label expected $expected +/- $atol, got $observed")
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

final_basis = report.physical_gausslet_final_basis_summary
h1 = report.physical_gausslet_h1_summary
h1_j = report.physical_gausslet_h1_j_summary
_check("final_dimension", _get(final_basis, :final_dimension), 471)
_check_approx(
    "final_overlap_identity_error",
    _get(final_basis, :final_overlap_identity_error),
    5.29668900282789e-14;
    atol = 1.0e-12,
)
_check_approx("h1_lowest", _get(h1, :lowest_energy), -0.7946037173365863;
    atol = 1.0e-12)
_check_approx("h1_j_self_coulomb", _get(h1_j, :self_coulomb), 0.4569117646737212;
    atol = 1.0e-12)
_check("residual_rank", _get(materialization, :residual_rank), 18)
_check("ida_orbital_dimension", _get(materialization, :ida_orbital_dimension), 489)
_check("augmented_dimension", _get(materialization, :augmented_dimension), 489)
println("h2_pqs_terminal_stage_smoke_elapsed_s=", harness_elapsed)
