# Nested Cartesian Fixed-Block Memory Refactor Plan

Date: 2026-04-12

## Scope

This note records a narrow repo-side plan to reduce memory use in the nested
Cartesian fixed-block build. The target is the support-space packet assembly in
`src/cartesian_nested_faces.jl`, not the later QW operator merge or HF layer.

The immediate goal is a simpler refactor first:

- keep the existing mathematics and support-space construction
- do not attempt a fully streamed / matrix-free redesign yet
- materially lower peak memory by creating at most one large
  `nsupport × nsupport` support matrix at a time and contracting it
  immediately

## Current Peak-Memory Source

The large allocations come from:

- `_nested_support_product_matrix(...)`
- `_nested_sum_of_support_products(...)`
- `_nested_weight_aware_pair_terms(...)`
- `_nested_shell_packet(...)`

The key issue is not the final fixed-block size. It is the temporary
support-space materialization before contraction to the `nfixed × nfixed`
packet.

For a large full-parent atomic case such as corrected one-center `ns9`:

- `nsupport ≈ 35937`
- one dense `Float64` support matrix is about
  `35937^2 * 8 ≈ 10.4 GiB`

The current code materializes several such objects in one shell-packet build:

- overlap support
- kinetic support accumulator
- `x`, `y`, `z`
- `x^2`, `y^2`, `z^2`
- one Gaussian-factor support matrix per expansion term
- one raw pair support matrix per expansion term

That is why the fixed-block stage dominates memory and runtime.

## Smallest Effective Refactor Boundary

The smallest effective refactor is:

- keep `_nested_support_product_matrix(...)` semantics
- but stop holding several support matrices live at once
- instead:
  - build one support matrix
  - contract it immediately with the nesting isometry
  - accumulate the resulting `nfixed × nfixed` matrix
  - reuse the support workspace for the next operator

This should be done before any deeper streaming redesign.

## Recommended Phase 1

### 1. Add a support-workspace fill helper

Refactor `_nested_support_product_matrix(...)` into a fill-style form, e.g.:

- `_nested_fill_support_product_matrix!(workspace, support_states, opx, opy, opz)`

with the current allocating wrapper preserved if needed for convenience.

This lets the packet builder reuse one large support workspace.

### 2. Add one immediate-contraction helper

Add a helper along the lines of:

- `_nested_contract_support_product!(dest, workspace, tmp, support_states, C, opx, opy, opz; alpha=1.0, beta=0.0)`

where:

- `workspace` is `nsupport × nsupport`
- `tmp` is `nfixed × nsupport`
- `C` is the support-restricted contraction matrix
- `dest` is `nfixed × nfixed`

Operationally:

1. fill `workspace`
2. compute `tmp = C' * workspace`
3. compute `dest = tmp * C`

This preserves the current algebra exactly while capping peak memory to one
large support matrix plus one moderate intermediate.

### 3. Rewrite `_nested_shell_packet(...)` to contract operator-by-operator

Do not hold these simultaneously at support scale:

- overlap
- kinetic
- `x`, `y`, `z`
- `x^2`, `y^2`, `z^2`

Instead:

- allocate the final `nfixed × nfixed` packet matrices up front
- reuse one support workspace and one contraction scratch
- contract each operator triple directly into its final fixed-block matrix

In particular, remove the current pattern:

- build `overlap_support`
- build `kinetic_support`
- build `position_x_support`
- ...
- then contract all of them later

That pattern is what creates the large live set.

### 4. Rewrite `_nested_sum_of_support_products(...)` away from support-scale accumulation

The current kinetic path is especially wasteful because it accumulates a dense
support-scale matrix.

Instead:

- contract each kinetic contribution directly into the final `nfixed × nfixed`
  destination
- sum at fixed-block scale, not support scale

So kinetic should become:

- final `K = contract(T,S,S) + contract(S,T,S) + contract(S,S,T)`

with each `contract(...)` using the shared support workspace.

### 5. Apply the same one-big-matrix-at-a-time policy to Gaussian terms

The Gaussian-term loop is already one-term-at-a-time, but it still allocates a
fresh `factor_support` every iteration.

Change it to:

- reuse the same support workspace
- contract directly into `gaussian_terms[term, :, :]`

### 6. Apply the same one-big-matrix-at-a-time policy to pair terms

The pair-term path should similarly:

- reuse the same support workspace
- avoid keeping a support-scale accumulator
- contract each term directly into the final `pair_terms[term, :, :]`

The `parent_weight_outer_*` pieces are small compared with the support matrix
and do not need to be the first optimization target.

## Expected Memory Effect

For large full-parent atomic cases, this should reduce peak memory from
"many large support matrices live" to roughly:

- one `nsupport × nsupport` workspace
- one `nfixed × nsupport` intermediate
- one support-restricted contraction matrix
- final `nfixed × nfixed` packet outputs

For the observed corrected `ns9` atomic scale, that means roughly:

- support workspace: about `10.4 GiB`
- `nfixed × nsupport` intermediate: about `0.5 GiB`
- support-restricted contraction matrix: about `0.5 GiB`
- final matrices: tens to low hundreds of MB

So Phase 1 should plausibly cut peak memory from the currently observed
`36–54 GiB` regime down into something much closer to the low-teens GiB range.

That is a significant improvement without changing the contract.

## Recommended Phase 2, only if needed later

If Phase 1 is still too large for the intended Cr atomic / Cr2 work, the next
step would be a deeper refactor:

- blocked or streamed contraction that never forms the full support matrix

But that is not the first move. The first move should be the simpler workspace
reuse / immediate-contraction rewrite above.

## Validation Plan

After the Phase 1 refactor, validate on one corrected full-parent atomic anchor
such as:

- one-center `ns9_family`, `Z = 10`, `rmax = 10`
or a nearby large case already used by the current diagnostics

Need to record:

- fixed-block construction wall time
- max RSS
- exact overlap / one-body agreement with the pre-refactor implementation

The correctness requirement is strict:

- the new assembly must reproduce the same fixed-block packet to numerical
  roundoff

## Recommendation

The next repo-doer pass should implement Phase 1 only:

- one support workspace at a time
- immediate contraction to `nfixed × nfixed`
- no full streaming redesign yet

That is the smallest refactor that should significantly improve nested
Cartesian memory footprint while keeping the current interface and math intact.
