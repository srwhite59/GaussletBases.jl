# Julia Style Notes

These are local Julia style preferences. They are meant to guide small cleanup passes and new code edits without changing algorithms or method contracts.

## Prefer Module-Owned Concepts Over Flat Field Clouds

When new code introduces a stable concept, prefer a small internal module with
explicit records and query functions over large flat `NamedTuple` field clouds.

A good Julia module should own one coherent idea, for example:

```text
geometry primitive
shellification policy
lowering contract
retained-unit plan
unit-pair inventory
pair-block materialization
local block collection
```

The module should expose a narrow interface:

```julia
plan = build_plan(input; options...)
summary(plan)
records(plan)
status(plan)
```

or similarly named constructors and accessors. Other code should use that
interface instead of reaching through many nested fields or manually copying
dozens of scalar fields.

Use structs for concepts that cross module boundaries or are consumed by later
stages. Use `NamedTuple`s mainly for compact summaries, local metadata, and
temporary compatibility adapters.

Avoid this pattern:

```julia
return (;
    terminal_shellification_unit_count = ...,
    terminal_shellification_unit_keys = ...,
    terminal_shellification_lowering_contract_count = ...,
    terminal_shellification_lowering_contract_kinds = ...,
    terminal_shellification_pair_inventory_available = ...,
    terminal_shellification_pair_operator_status = ...,
    ...
)
```

Prefer this pattern:

```julia
state = TerminalRouteState(...)
return (;
    terminal_route_state = state,
    terminal_route_summary = summary(state),
)
```

If old report fields must be preserved, derive them in one compatibility helper
at the report boundary. Do not pass compatibility aliases through every
intermediate stage.

A useful rule of thumb: if three or more related fields are being carried
together through multiple functions, consider whether they are really one object
owned by a module. Stable concepts crossing stage or module boundaries deserve a
module or typed object; local temporary values can stay ordinary local variables
or small `NamedTuple`s.

## Make Internal Modules Human-Facing

Even private/internal modules should be readable by a human who opens the file
cold.

Each module should have a top-level docstring that says:

1. what concept the module owns;
2. what inputs it consumes;
3. what outputs it produces;
4. what it explicitly does not do.

Prefer a short file map in the module front door:

```julia
# records.jl
#     Public/private record types.
#
# planning.jl
#     Constructors and plan-building logic.
#
# summaries.jl
#     Compact summaries for tests and reports.
```

Export only the intended interface. Keep helper functions private unless another
module really needs them. If another module needs a private-looking helper,
promote it to a deliberately named internal API rather than importing an
underscore function casually.

## Avoid Repeated Status-Flag Bundles

Do not repeatedly hand-write long bundles of status flags such as:

```julia
operator_blocks_materialized = false,
hamiltonian_data_materialized = false,
artifacts_materialized = false,
pqs_lowdin_materialized = false,
full_white_lindsey_route_assembled = false,
```

If the same flag bundle appears in several summaries, define a helper:

```julia
_no_global_assembly_flags()
_no_pqs_realization_flags()
_no_hamiltonian_artifact_flags()
```

or a compact status object.

Repeated flags are sometimes useful in final reports, but they should be
generated from one helper so that the meaning stays consistent and the code does
not become a field cloud.

## Keep Glue Thin And Directional

Glue code may connect modules, but it should not become the owner of the
algorithm.

Good glue:

```text
call module A
call module B
combine their summaries
return a structured state object
```

Bad glue:

```text
copy dozens of fields from A
copy dozens of fields from B
invent new scalar aliases
make downstream decisions from those aliases
```

Dependency direction should usually flow downward:

```text
driver / route state
-> shellification
-> lowering
-> retained units
-> transforms
-> pairs
-> pair-block materialization
-> assembly
```

Lower modules should not call back upward into driver helpers or report code. If
they need information, pass a small explicit object.

## Keep Boundary Comments Current

Comments that describe module boundaries are part of the architecture. When a
file changes from metadata-only to numerical, or from planning to materializing,
update the opening comments and docstrings in the same pass.

A stale boundary comment is worse than no comment, because later work may follow
the old contract.

## Test Through Public Or Module-Level Interfaces

Prefer tests that exercise the module's intended interface:

```julia
plan = retained_unit_plan(lowering_plan)
summary(plan)
records(plan)
```

Avoid tests that depend on deep internal field chains unless the test is
specifically about that internal representation.

Do not compare large staged objects with `==` or `===`. Compare compact
summaries, counts, keys, statuses, dimensions, and selected numerical values.

## Prefer Comprehensions For Simple Filled Arrays

When an array is allocated only to be filled by a simple elementwise formula, use a comprehension.

```julia
# Prefer
pp = [phi[k, i] * phi[k, j] for k = 1:lb, i = 1:m, j = 1:m]

# Instead of
pp = Array{T,3}(undef, lb, m, m)
for k = 1:lb, i = 1:m, j = 1:m
    pp[k, i, j] = phi[k, i] * phi[k, j]
end
```

Do not replace scratch buffers that are reused as `mul!` destinations or accumulator workspaces.

## Use Array Operations For Standard Prefix Data

Use Julia's array syntax for standard constructions such as prefix offsets.

```julia
# Prefer
offs = [0; cumsum(dims)]

# Instead of
offs = Vector{Int}(undef, length(dims) + 1)
offs[1] = 0
for i = 1:length(dims)
    offs[i + 1] = offs[i] + dims[i]
end
```

## Use Views For Reused Slices

When the same slice is used more than once and a copy is not needed, bind a view.

```julia
# Prefer
psiup_occ = @view psiup[:, 1:Nup]
rhoup = psiup_occ * psiup_occ'

# Instead of
rhoup = psiup[:, 1:Nup] * psiup[:, 1:Nup]'
```

Keep explicit copies when the code needs an owned array, for example assigning selected eigenvectors to the next orbital matrix.

## Combine Simple Nested Loop Headers

For simple Cartesian iteration with no setup between loop levels, put the iterators on one `for` line.

```julia
# Prefer
for np = 1:ns, nq = 1:ns, a = 1:nj, c = 1:nj
    p = (np - 1) * nj + a
    q = (nq - 1) * nj + c
    # ...
end

# Instead of
for np = 1:ns, nq = 1:ns
    for a = 1:nj, c = 1:nj
        p = (np - 1) * nj + a
        q = (nq - 1) * nj + c
        # ...
    end
end
```

Leave loop levels separate when the outer loop has meaningful setup, skip logic, timing sections, or temporary allocations.

## Use Keyword Shorthand

When forwarding a keyword with the same local variable name, use Julia's shorthand.

```julia
# Prefer
getblocksizes(N, m, blocksize; verbose)
_h1cache(range, phi, m, H, raH1; dofull)
solve_hfdmrg(H, V, psiup0; maxiter, blocksize, cutoff, verbose = false)

# Instead of
getblocksizes(N, m, blocksize; verbose = verbose)
_h1cache(range, phi, m, H, raH1; dofull = dofull)
solve_hfdmrg(H, V, psiup0;
    maxiter = maxiter, blocksize = blocksize, cutoff = cutoff, verbose = false)
```

Keep explicit `keyword = value` when the names differ or when a literal is being passed.

## Do Not Over-Compact Numerical Kernels

Compact syntax is not automatically better. Keep explicit loops when they make tensor contractions, Fock accumulation, cache construction, or sweep bookkeeping easier to audit. Style cleanups should be local and should not alter allocation strategy in hot paths unless that is the explicit goal.
