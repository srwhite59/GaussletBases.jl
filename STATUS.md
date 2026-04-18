# Current Status

This file is the short current **capability / trust matrix** for the repo.

It is not a changelog and it is not a roadmap. The question here is simply:

- what is mature enough to trust
- what is real but still experimental
- what is legacy or quarantined

## Mature

### Radial / atomic line

Trust this as the primary onboarding path.

Current mature surface includes:

- ordinary 1D, half-line, and radial gausslet bases
- explicit radial quadrature and diagnostics
- radial one-body operators
- explicit atomic `(l,m)` channels
- static atomic IDA ingredients
- direct / exchange / Fock helpers
- minimal UHF kernel
- dense and sliced atomic export for the current density-density model

Primary docs:

- `README.md`
- `docs/src/explanations/current_atomic_branch.md`
- `docs/src/manual/index.md`

### Exact Cartesian overlap / projector / transfer primitives

Trust this as a real current workflow primitive, not as a side note.

Current mature surface includes:

- `cross_overlap`
- `basis_projector`
- `transfer_orbitals`

These exact Cartesian transfer contracts are now part of the real current repo
story and support both validation and workflow handoff across current Cartesian
representation families.

## Real but experimental

### Ordinary / Cartesian mapped and Qiu-White routes

This is a real current workflow surface, but newer and less settled than the
radial line.

Current supported surface includes:

- mapped ordinary Cartesian basis construction
- exact overlap / projector / transfer primitives
- current Qiu-White residual-Gaussian route
- current supplement and compact packet contracts

Primary docs:

- `docs/src/explanations/current_ordinary_branch.md`
- `docs/src/algorithms/qiu_white_residual_gaussian_route.md`

### One-center nested Cartesian path

This is real code, not a sketch.

Current supported surface includes:

- one-center shell-sequence construction
- one-center fixed-block construction
- nested Cartesian operator consumers
- compact-only packet construction on the production path

This line is still evolving scientifically, especially around frozen-core and
exactification questions, but it is no longer accurate to describe it as absent
or purely aspirational.

### Bond-aligned diatomic workflow

This is also real code, not a placeholder.

Current supported surface includes:

- bond-aligned diatomic nested fixed-source construction
- supported one-build source reuse
- source-level geometry diagnostics
- source geometry payloads and plane slices
- bond-aligned diatomic Qiu-White workflow support
- current diatomic geometry policy and diagnostics

This line is still experimental in the scientific sense, especially around
policy and workflow maturation, but it is a real current repo surface.

### Experimental chain / square-lattice nested producers

These routes exist and can be used for producer-side experiments, but they are
not yet settled as broad public workflows.

### Angular research track

The angular line is real and has:

- shell-local injected basis construction
- shell-to-atom assembly
- one-electron, HF-style, and small-ED benchmark paths
- direct in-memory HFDMRG payload handshake

But it should still be read as experimental/manuscript-facing rather than as a
stable general user workflow.

## Legacy / quarantined

### Old 1D COMX-cleaned hybrid ordinary route

This should be treated as legacy/internal code for surrogate comparisons and
historical regression checks.

It is **not** the supported current ordinary / Cartesian workflow.

### Flat supporting-note history in `docs/`

The broader flat `docs/*.md` note tree remains useful as supporting or
historical material, but it is not the primary current authority for repo
status. Use the root docs and rendered `docs/src/*` pages first.

## Not yet broad or stable

The following are still outside the current mature/public workflow contract:

- broad exact four-index electron-electron workflows
- a broad stabilized molecular HF / post-HF workflow
- a broad stabilized molecule-scale nested workflow beyond the current
  one-center and bond-aligned diatomic surfaces
- named chemistry basis-set libraries
- large solver layers such as DMRG as a native repo workflow
- Python / Fortran interoperability layers

Those absences are now narrower than the old root docs implied, but they are
still real.
