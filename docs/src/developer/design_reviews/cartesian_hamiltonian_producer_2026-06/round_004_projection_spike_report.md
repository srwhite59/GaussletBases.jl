# Round 004 Terminal Projection Spike Report

Reporter: repo-doer@macmini

Date: 2026-06-20

Status: uncommitted numerical spike; no source, test, tool, bin, docs, or policy
edits by the spike.

## Commands Reported

- `julia --project=. tmp/work/terminal_projection_introspection.jl`
- `julia --project=. tmp/work/terminal_projection_spike.jl`
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git status --short --branch`

Doer reported validation passed and no tracked worktree changes.

## One-Center Atomic PQS

- elapsed: `15.591 s`
- terminal records: `4`
- roles: `(:atom_local_core, :atom_local_shell, :atom_local_shell, :atom_local_shell)`
- retained dimension: `419`
- direct core support: `125`
- direct identity error: `3.331e-15`
- direct IDA weight range: `1.607e-1 .. 5.446e-1`
- PQS shell retained count: `98` each
- raw/projected overlap norms:
  - `7.892e-16 / 1.972e-31`
  - `7.265e-16 / 3.944e-31`
  - `2.646e-16 / 9.861e-32`
- Gram min range: `0.9960 .. 0.9994`
- ranks: `98`
- coefficient memory per shell: `0.163 .. 0.450 MiB`

## H2 Contact-Core Terminal Topology

- elapsed: `10.701 s`
- terminal records: `3`
- roles: `(:atom_contact_core, :shared_molecular_shell, :shared_molecular_shell)`
- retained dimension: `471`
- direct core support: `275`
- direct identity error: `1.787e-14`
- direct IDA weight range: `3.592e-1 .. 9.693e-1`
- PQS shell retained count: `98` each
- raw/projected overlap norms:
  - `6.881e-15 / 2.367e-30`
  - `1.524e-15 / 5.916e-31`
- Gram min range: `0.9500 .. 0.9521`
- ranks: `98`
- coefficient memory per shell: `0.271 .. 0.432 MiB`

## Cr2 Separated Terminal Topology

- elapsed: `75.295 s`
- terminal records: `19`
- roles: `2` atom-local cores, `8` atom-local shells, `1` midpoint slab,
  `6` shared molecular shells, `2` outer mismatch slabs
- retained dimension: `4291`
- direct supports: atom cores `125 + 125`, midpoint slab `169`, outer slabs
  `1250 + 1250`
- direct identity errors: `1.215e-13 .. 1.893e-13`
- direct IDA weight ranges:
  - atom cores: `4.055e-3 .. 1.615e-2`
  - midpoint slab: `1.935e-2 .. 8.113e-1`
  - outer slabs: `5.279e-2 .. 1.419e1`
- atom-local PQS shells: raw overlaps `8.373e-15 .. 4.612e-14`,
  projected overlaps `3.155e-30 .. 9.466e-30`, Gram min `0.9951 .. 0.9996`,
  ranks `98`, coefficient memory `0.163 .. 0.647 MiB`
- shared molecular PQS shells: raw overlaps `6.422e-15 .. 3.761e-14`,
  projected overlaps `1.578e-30 .. 1.578e-29`, Gram min `0.9025 .. 0.9455`,
  ranks `98`, coefficient memory `1.467 .. 3.590 MiB`

## Design Findings

- The same terminal retained/support/transform record shape is usable for
  one-center, H2, and Cr2.
- One-center still required the old skeleton route-shape input to reach
  terminal records; this must be clarified before freezing Slice A.
- Previous-block projection appeared numerically stable in this spike.
- Raw cross-overlaps were already small and projected overlaps dropped to
  roundoff.
- All PQS shell ranks stayed at the expected `98`.
- Effective support growth was modest in this approximation: each PQS shell
  stayed at its own support size.
- The spike used existing shell-local coefficients for previous PQS blocks, not
  fully recursively expanded projected coefficients. Production Slice A still
  needs the true accumulated-coefficient/support audit.
- Dense scratch projection was already a `75 s` Cr2 spike, so production Slice A
  should use factorized or incremental overlap construction.
