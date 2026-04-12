# Fixed-Radial Angular Sequence Export Milestone

This note records the repository milestone introduced by commit `4d04591`
(`Add fixed-radial angular sequence export`).

This is a real separate milestone from the earlier angular HF / consumer-boundary
work. It is also distinct from later DMRG-side continuation logic. The repo now
owns a concrete producer contract for an increasing-`N_sph` fixed-radial
sequence line.

## New producer contract

The repo now supports one fixed-radial increasing-`N_sph` sequence line with:

- one dense native level artifact per `N_sph`
- one adjacent shell-local overlap sidecar per neighboring
  `N_sph[k] -> N_sph[k+1]` pair
- stable radial shell ids across the sequence
- stable physical radial shell centers across the sequence
- deterministic per-level shell dimensions
- deterministic within-shell labels from the cached shell-local angular
  profiles

The sequence line is producer-side only. It is intended to hand clean per-level
Hamiltonians and compact adjacent-overlap data to external continuation
consumers.

## Key scientific/software commitment

Under the current fixed-profile shell-local construction, the shell-local
angular basis is radius-independent for fixed:

- point-set source / `N_sph`
- `beta`
- `l_inject`
- `tau`
- `whiten`
- gauge version

That means the adjacent shell-local profile overlap for
`N_sph[k] -> N_sph[k+1]` is common across shells. The repo therefore exports
that overlap once as a shell-independent sidecar rather than as a shell-by-shell
family.

This is a real contract commitment, not just an implementation accident.

## Plumbing introduced by this milestone

The main implementation surface is:

- `src/angular_sequence_export.jl`

The sequence line depends on the shell-local profile/gauge reuse added to:

- `src/angular_shell_basis.jl`

and the supporting assembly reuse in:

- `src/angular_shell_assembly.jl`

The public exports are wired through:

- `src/GaussletBases.jl`

## Public interface

The main public entry points are:

- `build_atomic_fixed_radial_angular_sequence`
- `atomic_fixed_radial_angular_level_dense_payload`
- `write_atomic_fixed_radial_angular_level_jld2`
- `atomic_fixed_radial_angular_overlap_sidecar_payload`
- `write_atomic_fixed_radial_angular_overlap_sidecar_jld2`

These entry points define the current repo-owned producer surface for the
fixed-radial increasing-`N_sph` line.

## Current trust boundary

This is a producer-side continuation contract.

It does **not** yet include:

- common-target embedding export
- smaller-to-max lift matrices
- DMRG-side Givens generation
- restart ladders
- many-body orchestration

It also does **not** by itself make `GaussletBases` the demonstrated
producer-of-record for the full capstone DMRG path. It only settles the
producer-side fixed-radial sequence export contract.

## Validation

Validation for this milestone is currently narrow and explicit.

The trust story includes:

- He-focused targeted regression coverage in `test/runtests.jl`
- targeted checks for:
  - shared radial basis id across sequence levels
  - stable shell ids / shell centers across levels
  - cached profile reuse
  - dense level payload shape/metadata
  - adjacent-overlap sidecar shape/metadata
  - JLD2 roundtrips for one level artifact and one adjacent sidecar

Docs build status:

- passed with:
  - `env JULIA_DEPOT_PATH=/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/tmp/julia_depot julia --project=docs docs/make.jl`

Full test status:

- no fresh full grouped angular test rerun is claimed here beyond the targeted
  fixed-radial sequence coverage already reported for this line

So the present trust statement is:

- producer contract and He-targeted regression coverage: yes
- broad full-angular grouped validation for the entire branch: not claimed by
  this milestone note

## Relation to earlier and later work

- earlier angular HF / consumer-boundary work established narrow benchmark and
  payload/export boundaries for the experimental angular line
- this milestone is specifically the producer-side fixed-radial increasing-`N_sph`
  sequence contract
- later DMRG-side lifting, restart ladders, and capstone continuation logic are
  downstream of this milestone and are not part of it

## Practical consequence

The repo is now ready to support He-first adjacent-`N_sph` continuation studies
from the producer side.

After this milestone, the main next risk shifts away from the producer export
contract and toward the external DMRG-side lifting / restart logic that consumes
the exported levels and adjacent overlaps.
