#!/usr/bin/env julia

# Developer smoke for the active H2 PQS terminal-basis route.
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

function _check_le(label, observed, bound)
    observed <= bound || error("$label expected <= $bound, got $observed")
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

basis = transforms.terminal_basis_realization
_check("terminal_basis_block_count", length(basis.blocks), 3)
_check("terminal_basis_support_counts",
    Tuple(length(block.support_indices) for block in basis.blocks),
    (275, 637, 1215))
_check("terminal_basis_column_ranges",
    Tuple(block.column_range for block in basis.blocks),
    (1:275, 276:373, 374:471))
_check("terminal_basis_final_dimension", basis.final_dimension, 471)
_check_le("terminal_basis_max_cross_overlap", basis.max_cross_overlap, 1.0e-8)
println("h2_pqs_terminal_stage_smoke_elapsed_s=", harness_elapsed)
