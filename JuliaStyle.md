# Julia Style Notes

These are local Julia style preferences. They are meant to guide small cleanup passes and new code edits without changing algorithms or method contracts.

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
