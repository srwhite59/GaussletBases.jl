# Post-Decomposition Efficiency Resume Plan

## Purpose

This note records the efficiency-cleanup status after the structural
decomposition work and the first post-decomposition cleanup passes.

It is a resume memo:

- what has already been completed
- what cleanup remains
- what order to do it in
- what is still suspicious but not yet justified for removal

This is intentionally narrower than the earlier
[`cartesian_nested_decomposition_plan.md`](cartesian_nested_decomposition_plan.md).
That plan was about file ownership and extraction. This note is about the next
round of code-surface reduction now that the split is finished.

## Current Status

### Completed structural decomposition

The seven planned decomposition chunks are complete:

1. nested timing/reporting extracted
2. one-center atomic nested core isolated
3. bond-aligned diatomic geometry/source isolated
4. experimental chain/square nested producers isolated
5. QW basis/front-door layer isolated
6. QW residual diagnostics/stabilization isolated
7. QW raw block builders split from final operator assembly

The previous mixed-responsibility files are no longer the main shape problem:

- `src/cartesian_nested_faces.jl`
- `src/ordinary_qiu_white_rg.jl`

In particular, `src/ordinary_qiu_white_rg.jl` has now been retired and deleted.

### Completed post-decomposition cleanup

Two follow-up cleanup passes are already done:

1. first dead-surface cleanup
2. `TimeG`-only timing cleanup

Completed removals include:

- public/shared full-debug nested packet mode
- public `retain_term_tensors` seams
- sum-to-term fallback loops in operator assembly
- dead underscored diatomic fixed-source/fixed-block wrappers
- stale compatibility shell `src/ordinary_qiu_white_rg.jl`
- legacy nested timing collector/report machinery
- old ordinary/QW coarse phase timer helpers

Timing now has one authority:

- `TimeG`

The only remaining user-facing timing report surface is a thin `TimeG`-backed
wrapper.

## Best Next Cleanup Order

The next passes should be:

1. broader metadata/options sweep
2. targeted redundancy cleanup via semantically owned helpers
3. optional internal-module pass only if the value still justifies it

That ordering matters. The repo should delete stale behaviorless surface before
trying to tidy namespaces.

## Next Pass: Metadata And Single-Valued Option Sweep

This is the next highest-value cleanup.

### Primary targets

- `term_storage` fields and route metadata that may no longer carry real
  behavior
- old timing/report metadata that survived the collector removal
- single-valued mode selectors that still exist in signatures or internal route
  plumbing
- struct fields that are now always `nothing`, always one fixed mode, or only
  retained for stale compatibility

### Files most likely in scope

- `src/cartesian_nested_faces.jl`
- `src/cartesian_nested_reporting.jl`
- `src/cartesian_nested_atomic.jl`
- `src/cartesian_nested_diatomic.jl`
- `src/cartesian_basis_representation.jl`
- `src/ordinary_qw_nested_frontends.jl`
- `src/ordinary_qw_raw_blocks.jl`
- `src/ordinary_qw_operator_assembly.jl`
- `test/runtests.jl`

### Specific suspicious surfaces

These are suspicious, but were not yet removed because they still had live
representation or test footprint:

- `term_storage` metadata on nested packet/fixed-block structs
- route metadata that still threads through `basis_representation(...)`
- legacy `packet.gaussian_terms` / `packet.pair_terms` compatibility in tests
  where those hits now refer to data fields rather than timing labels

### Rule for removal

Remove a field, option, or metadata surface only if all of the following are
true:

- no production path uses it
- no experimental path meaningfully changes behavior because of it
- tests that mention it are only pinning compatibility residue rather than real
  capability

If the dependency picture is still ambiguous, leave it in place and record the
reason.

## Pass After That: Targeted Redundancy Cleanup

This should not be a generic "make helper functions" pass.

The goal is to remove repeated setup/policy logic while keeping hot numeric
paths visible and specialization-friendly.

### Good targets

- repeated source/fixed-block/diagnostic wrapper patterns in
  `src/ordinary_qw_nested_frontends.jl`
- repeated supplement/raw-block plumbing patterns in
  `src/ordinary_qw_raw_blocks.jl`
- repeated operator-preparation and reassembly patterns in
  `src/ordinary_qw_operator_assembly.jl`
- repeated geometry/report formatting scaffolding in
  `src/cartesian_nested_reporting.jl` and
  `src/cartesian_nested_experimental_geometries.jl`

### Bad targets

Do not extract common helpers just because two hot kernels look similar.

In particular:

- do not merge hot packet kernels for stylistic reasons
- do not hide distinct geometry policies behind overly generic utility layers
- do not replace clear ownership with grab-bag `utils` files

Helper extraction should follow semantic ownership, not line-count reduction by
itself.

## Internal Modules: Defer Unless Still Needed

Internal modules may still be useful later, but they should not be the next
step.

The file split already delivered most of the ownership benefit with less
complexity than adding new module/import boundaries.

Only reconsider internal modules if, after the metadata/options cleanup and the
redundancy pass, there is still clear value in enforcing namespace boundaries.

If that becomes worthwhile, the likely candidates are:

- one internal nested Cartesian block around the `cartesian_nested_*` files
- one internal ordinary/QW block around the `ordinary_qw_*` files

That should be treated as optional, not assumed.

## Performance Notes

The most likely shape-related performance wins from here are not from helper
extraction. They are from:

- deleting stale branches
- deleting behaviorless mode selectors
- deleting dead metadata threading
- keeping compact production paths explicit

Potential remaining pressure points to watch during cleanup:

- runtime `Symbol` branches that now have only one real production mode
- dense/sparse or compact/noncompact storage unions that remain on hot structs
- wrapper layers that obscure source/bundle reuse and encourage rebuilding

The main rule is:

- deduplicate setup/policy layers first
- leave hot math visible unless there is a clear measured benefit

## Validation Surfaces For Resume

When resuming cleanup, the smallest meaningful validation set should still hit:

### Nested / timing / atomic

- `One-center atomic full-parent nested contract`
- `One-center atomic legacy-profile nested contract`
- `One-center atomic fixed-block timing surface`
- `Global timing macro surface`

### Diatomic

- `Bond-aligned diatomic nested source reuse path`
- `Bond-aligned diatomic nested fixed block`
- `Bond-aligned diatomic compact nested fixed-block contract`
- `Bond-aligned diatomic nested QW consumer path`
- `Bond-aligned diatomic doside / COMX trace diagnostics`

### Experimental

- bond-aligned homonuclear chain ordinary/QW reference path
- axis-aligned homonuclear square-lattice ordinary/QW reference path

### Operator assembly

- `assembled_one_body_hamiltonian`
- `ordinary_cartesian_vee_expectation`
- `ordinary_cartesian_1s2_check`

Direct local harnesses are acceptable if broad group runs are too slow, but the
report should state exactly what was covered.

## Resume Point

The decomposition project itself is complete.

The next concrete resume point is:

1. run the metadata/options sweep
2. remove any unambiguously dead fields/keywords/route metadata
3. rewrite tests that only pin stale compatibility surface
4. only then evaluate whether a targeted redundancy pass is still worth doing

That is the lowest-risk path that will continue reducing code bloat without
hiding real behavior behind a new abstraction layer.
