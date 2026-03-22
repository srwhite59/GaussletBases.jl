## Cartesian Nested Fixed-Block Adapter

The generalized shell packet is already a valid nonseparable three-dimensional
fixed packet candidate:

- it carries one shell-level coefficient matrix
- it carries the propagated fixed-fixed overlap / one-body packet
- and it is already localized on disjoint shell-face interiors

So the next pre-recursion step is not another basis primitive. It is a narrow
consumer adapter.

The first target is the QW-PGDG nearest/GGT route. That path already has the
right downstream algebra:

- residual-space construction from a fixed block plus an added Gaussian channel
- exact raw one-body assembly
- and nearest/GGT interaction assembly in the same two-index IDA form

What changes in this pass is only the fixed block.

Instead of rebuilding a separable 3D fixed block from one finalized 1D line,
the consumer reads a preassembled shell-level fixed block:

- fixed-fixed data come directly from the nested shell packet
- fixed-Gaussian raw blocks are built by contracting the parent raw
  fixed-to-Gaussian blocks through the shell coefficient matrix
- Gaussian-Gaussian blocks stay on the existing analytic route

So the first nested consumer should look like:

- same downstream QW-PGDG residual / one-body / nearest-GGT algebra
- different fixed block source
- no recursion yet
- no MWG yet

If this adapter works cleanly, the next remaining pre-recursion step is to let
the existing Cartesian/QW-PGDG fixed-block consumers read this shell-packet
interface more uniformly, rather than through one special nested adapter.
