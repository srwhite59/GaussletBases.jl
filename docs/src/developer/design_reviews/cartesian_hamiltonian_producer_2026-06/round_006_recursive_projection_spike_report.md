# Round 006 Recursive Projection Spike Report

Reporter: repo-doer@macmini

Date: 2026-06-20

Status: uncommitted numerical spike; only ignored `tmp/work/*.jl` files touched.
No tracked files changed and no commit was made by the spike.

## Commands Reported

- `julia --project=. tmp/work/terminal_projection_recursive_spike.jl`
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git status --short --branch`

Doer reported validation passed and a clean tracked worktree on
`design/cartesian-hamiltonian-producer`.

## Entry Point

All three cases used the same terminal realization inputs once reached:

```text
terminal_retained_rule_plan
retained_unit_transform_contract_plan
existing shell-local projection kernel
```

Caveat: one-center atomic still needs the old skeleton route shape
`(:pqs_left, :product, :pqs_right)` to get through current public staging. The
terminal records are usable after that, but this route-shape mismatch should not
be frozen as the intended public contract.

## Projection Settings

- `projection_atol = 1e-12`
- true recursive accepted blocks were used
- direct blocks stayed implicit identity on support
- PQS blocks carried recursively accepted effective support and coefficients
- no global parent overlap or global final overlap was formed

With this tolerance, every residual was skipped; no subtraction was numerically
justified. Raw residuals were all near roundoff.

## One-Center Atomic

- time: `21.643 s`
- records: `4`
- roles: `(:atom_local_core, :atom_local_shell, :atom_local_shell, :atom_local_shell)`
- shell supports: `218`, `386`, `602`
- retained ranks: all `98`
- raw residuals: `1.00e-16 .. 7.89e-16`, all skipped
- final cross overlaps: max `8.073e-16`
- Gram min: `0.9960 .. 0.9994`
- effective supports: unchanged
- PQS coefficient memory: `0.902 MiB`
- largest workspace: `2.765 MiB`

## H2 Contact-Core

- time: `18.766 s`
- records: `3`
- roles: `(:atom_contact_core, :shared_molecular_shell, :shared_molecular_shell)`
- shell supports: `362`, `578`
- retained ranks: both `98`
- raw residuals: `6.38e-16 .. 6.88e-15`, all skipped
- final cross overlaps: max `3.095e-15`
- Gram min: `0.9500`, `0.9521`
- effective supports: unchanged
- PQS coefficient memory: `0.703 MiB`
- largest workspace: `2.549 MiB`

## Cr2 Separated

- time: `83.568 s`
- records: `19`
- roles: `2` atom-local cores, `8` atom-local shells, midpoint slab,
  `6` shared molecular shells, `2` outer mismatch slabs
- total retained dimension: `4291`
- shell ranks: all `98`
- raw residuals: about `1e-17 .. 4.61e-14`, all skipped
- final cross overlaps: max `5.463e-14`
- Gram min range: `0.9025 .. 0.9996`
- effective supports: unchanged
- PQS coefficient memory: `17.911 MiB`
- largest workspace: `175.928 MiB`, from the largest candidate self-overlap
- later direct cross-overlaps:
  - midpoint slab vs prior blocks: `5.267e-15`
  - low outer slab vs prior blocks: `3.010e-15`
  - high outer slab vs prior blocks: `2.432e-15`

## Design Findings

- Slice A is numerically stable for these three cases with
  `projection_atol = 1e-12`.
- A tighter threshold could subtract roundoff residuals and artificially grow
  effective support.
- Cr2's largest local workspace is already about `176 MiB`; production should
  tile or stream local overlap actions when they exceed the approved cap.
- One-center and diatomic terminal records can share the realization entry
  point, but one-center public staging still has route-shape/skeleton cleanup
  pending before the contract is clean.
