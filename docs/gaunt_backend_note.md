# Angular Gaunt backend note

This note explains the next small angular-structure correction in the atomic line.

The goal is not to port the old Ylm/Coulomb stack wholesale. The goal is narrower:

- keep the present radial and atomic public objects
- replace the current hand-built dense Gaunt construction in `src/atomic_ida.jl`
- use a cleaner sparse/block Gaunt backend underneath

## 1. What should be ported from `GauntTables`

The highest-value parts of the legacy `GauntTables.jl` code are:

- sparse storage of Gaunt coefficients by `(L, l_1, l_2)` block
- exact triangle and parity selection-rule pruning
- explicit iteration over nonzero blocks and entries
- simple lookup helpers for tests and diagnostics

Those pieces are structurally clean and reusable. They improve the present
atomic line without forcing the repository to adopt the larger legacy Ylm
application layer.

## 2. What should be adapted for this repository

The port here should be light rather than literal.

In particular:

- keep the sparse/block table structure
- keep the selection-rule organization
- keep the distinction between `basis = :complex` and `basis = :real`

but:

- do **not** pull in the old larger Ylm/Coulomb modules
- do **not** require a new external angular-special-function dependency
- reuse the repository’s current internal Wigner-3j evaluation instead

So the result is best understood as a lightly adapted internal equivalent of
`GauntTables`, not a wholesale transplant.

## 3. What should be replaced in `src/atomic_ida.jl`

The present file builds Gaunt tensors by:

1. looping directly over all channel pairs
2. evaluating coefficients on the fly
3. then building dense `Q_L` kernels from those dense tensors

The part that should be replaced is the first step:

- the direct dense construction of Gaunt tensors from channel-pair loops

Instead, `src/atomic_ida.jl` should:

1. build a sparse/block Gaunt table once
2. form the current dense `gaunt_tensor` output from that table
3. build the current dense `angular_kernel` output from the same sparse data

That keeps the present public layer while giving the internal angular backend a
much cleaner structure.

## 4. What should stay unchanged publicly

The public atomic story should remain the same.

In particular, it is desirable to keep:

- `AtomicIDAOperators`
- `gaunt_tensor`
- `gaunt_coefficient`
- `angular_kernel`

with the same scientific meaning as before.

This pass is about strengthening the internal angular backend, not about making
users learn a new public layer.

## 5. Why this is the right middle path

This repository is now beyond the point where the current dense/custom Gaunt
construction should remain the permanent solution.

But it is still too early to port the whole legacy Ylm/Coulomb stack.

So the right middle path is:

- port the Gaunt table backend
- keep the present atomic objects
- defer larger `DiagYlm` / `VeeYlm` application-layer ideas

That gives better selection-rule structure and a better future scaling path
without dragging in the full older infrastructure.

## 6. What may become worthwhile after this

Once this backend is in place, two small later helper ports become more
attractive:

- sector organization by conserved `m` sums
- packed orbital/block indexing helpers

Those ideas appear in the older Ylm code and could later help the atomic IDA
layer scale better.

But they should still come **after** the Gaunt table backend itself, not before.
