# R3 Usability Supplemented Workflow

Status: approved internal supported facade for the residual-GTO/MWG
supplemented PQS Hamiltonian path. This is not an exported public API.

This amendment closes the immediate usability gap after R3-A/B/C: callers
should not have to manually reconstruct the base stages, load a supplement,
call the R3 same-construction function, and then call the R3-C writer.
Implementation of this facade should wait until the owner-local residual
selection correction in `r3_residual_gto_mwg_augmentation.md` has measured and
recorded the corrected H2 MWG scalar.

It approves only a narrow module-qualified internal facade for z-axis H2 and
internal/performance-supported z-axis Be2. Cr2 remains deferred.

## Approved IDs

- `HP-R3U-FILE-01` - exact source and test files for the usability facade.
- `HP-R3U-FN-01` - non-exported supplemented Hamiltonian facade.
- `HP-R3U-WIRE-01` - same-construction base/R3/R3-C wiring.
- `HP-R3U-TEST-01` - standalone H2 usability endpoint validation.

## Decision

The first usability surface is an internal supported function, not an exported
public API and not a driver/tool workflow. The function name is:

```julia
cartesian_residual_gto_mwg_hamiltonian
```

It may be called as a module-qualified internal function, but it must not be
added to the public export list in `src/GaussletBases.jl`. A later public API
may choose a broader name such as `cartesian_hamiltonian`; this amendment does
not freeze a public supplemented-Hamiltonian name.

Approved call shape:

```julia
cartesian_residual_gto_mwg_hamiltonian(
    system::NamedTuple;
    basis::NamedTuple,
    supplement::NamedTuple,
    hamfile::Union{Nothing,AbstractString} = nothing,
)::CartesianIDAHamiltonian{Float64}
```

The return is the existing `CartesianIDAHamiltonian{Float64}` directly. The
function must not return a wrapper, status object, report object, payload, or
`(value, metadata)` pair.

## Approved Files

Primary owner file:

```text
src/cartesian_base_hamiltonian.jl
```

This file may add the non-exported facade, input validation, supplement-spec
normalization, and the one-call wiring from base spec to R3 construction.

Existing R3 owner file:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

This file may be touched only if the implementation needs a small local seam
to reuse the same-construction path and R3-C writer without recomputing the
residual object. It must not add a new artifact schema, public API, status
object, or broad provider/cache object.

Approved validation file:

```text
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
```

`HP-R3U-TEST-01` may extend this existing standalone endpoint gate. It does
not approve a new committed test file or inclusion in `test/runtests.jl`.

No edit is approved for `src/GaussletBases.jl`; no export is approved. If the
implementation cannot avoid a root include/export change, stop and return for
a docs-only amendment.

## System Input

`system` uses the same plain `NamedTuple` style as R1. Required keys:

- `atom_symbols::AbstractVector`
- `nuclear_charges::AbstractVector{<:Real}`
- `atom_locations::AbstractVector{<:NTuple{3}}`
- `nup`
- `ndn`

Center-sized collections must be vectors or other `AbstractVector` values, not
variable-size tuples. Unknown keys must throw `ArgumentError`.

Supported first systems:

- z-axis H2:
  - `atom_symbols = ["H", "H"]`;
  - `nuclear_charges = [1.0, 1.0]`;
  - `nup = 1`, `ndn = 1`;
  - both centers have `x = 0` and `y = 0`;
  - distinct finite `z` coordinates.
- z-axis Be2, internal/performance-supported only:
  - `atom_symbols = ["Be", "Be"]`;
  - `nuclear_charges = [4.0, 4.0]`;
  - `nup = 4`, `ndn = 4`;
  - both centers have `x = 0` and `y = 0`;
  - distinct finite `z` coordinates.

Be2 support is for internal usability and performance sanity, not a public API
guarantee and not a committed validation gate in the first implementation.

Unsupported systems must throw clear `ArgumentError`s before expensive
construction where practical. Cr, Cr2, ECP systems, heteronuclear molecules,
non-z-axis diatomics, translated/rotated general molecules, RHF/solver
handoff, and non-base Hamiltonian variants are not approved.

## Base Basis Input

`basis` is a plain `NamedTuple`. Required keys for the first diatomic scope:

- `q`
- `core_spacing`
- `xmax_parallel`
- `xmax_transverse`

Optional keys and defaults:

- `parent_axis_family = :G10`
- `reference_spacing = 1.0`
- `tail_spacing = 10.0`

Validation:

- `q` must be a positive integer;
- spacing and extents must be finite and positive;
- `parent_axis_family` must remain the validated R1 value `:G10`;
- unknown keys throw `ArgumentError`.

The usability facade has no public `method`, `route`, `n_s`, `bond_axis`,
`bond_length`, `radius`, `d`, mapping backend, parent-axis-count, or output
group selector. Those are either derived internally or unsupported.

## Supplement Input

`supplement` is a plain `NamedTuple`. Required keys:

- `basis_by_center::AbstractVector{<:AbstractString}`
- `lmax`

Optional keys and defaults:

- `uncontracted = false`
- `width_filtering = nothing`

Validation and normalization:

- `basis_by_center` length must equal the number of centers.
- First usability scope is homonuclear: all atom symbols must match and all
  `basis_by_center` values must match. This keeps the implementation on the
  existing `legacy_bond_aligned_diatomic_gaussian_supplement` route. A
  heteronuclear basis-by-center route requires a later amendment.
- `lmax` must be an integer with `0 <= lmax <= 6`, matching the existing
  legacy Cartesian shell cap.
- `uncontracted` must be `Bool`.
- `width_filtering` must be either `nothing` or a `NamedTuple` with exactly
  `max_width`, where `max_width` is finite and positive. This maps to the
  existing legacy `max_width` filter.
- Unknown supplement keys throw `ArgumentError`.

First H2 validation fixture:

```julia
supplement = (;
    basis_by_center = ["cc-pVTZ", "cc-pVTZ"],
    lmax = 1,
    uncontracted = false,
    width_filtering = nothing,
)
```

Basis-file lookup uses the existing legacy named-basis resolution order. This
amendment does not approve a public basis-file path selector.

## Wiring Contract

The facade must construct all base and supplement objects inside one call:

```text
validated system/basis/supplement spec
-> R1-style/base producer normalization and base stages
-> base CartesianIDAHamiltonian plus same-construction terminal basis/bundles
-> legacy named-basis supplement loading
-> basis_representation(supplement)
-> R3 same-construction augmented Hamiltonian path
-> optional R3-C artifact writer
-> CartesianIDAHamiltonian{Float64}
```

The facade must reuse the R1/base producer construction path where applicable.
For H2, it should reuse the existing R1 validation/stage helpers rather than
creating a second public-shaped base route. For Be2, it may add a generalized
internal z-axis homonuclear diatomic normalization in the same owner file. It
must not call the public `cartesian_base_hamiltonian` and then reconstruct
terminal basis state separately. The base Hamiltonian, terminal basis
realization, and parent axis bundle used for R3 must come from the same base
construction call.

The facade must reuse the R3 same-construction augmented Hamiltonian path and
R3-C provenance writer. It may refactor local R3 internals enough to avoid
recomputing residual objects for artifact writing, but it must not introduce a
persistent raw-block bundle, cache object, status object, payload, report
field, or public stage object.

The facade must not expose or return terminal basis realizations, bundles,
residual objects, augmented-operator objects, MWG descriptors, pair factors,
or provenance payloads.

## Artifact Contract

If `hamfile === nothing`, no artifact is written.

If `hamfile !== nothing`, the facade writes the supplemented Hamiltonian using
the existing Cartesian IDA Hamiltonian artifact shape plus the approved
`HP-R3-ART-01` `supplement_provenance/` group. It still returns the in-memory
Hamiltonian.

An empty `hamfile` must throw `ArgumentError`. File-system errors from the
existing writer should propagate normally.

Production must not require readback. Readback is validation-only.

This amendment does not approve:

- a new Hamiltonian wrapper;
- a new artifact format;
- new artifact keys beyond `supplement_provenance/`;
- a separate manifest;
- broad residual-basis serialization;
- a public provenance reader.

## Validation

Approved committed validation remains a standalone endpoint gate:

```text
julia --project=. test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
```

`HP-R3U-TEST-01` may extend that file with one usability-facade section that:

- calls `cartesian_residual_gto_mwg_hamiltonian` on the H2 fixture above;
- writes to `mktempdir()`;
- validates returned type `CartesianIDAHamiltonian{Float64}`;
- validates augmented dimension `489`;
- validates lowest augmented one-body orbital IDA self-Coulomb against the
  remeasured owner-local residual-selection scalar within the approved
  tolerance;
- validates readback Hamiltonian matrices against the returned Hamiltonian
  within tight deltas;
- validates `supplement_provenance/` keys match the normalized supplement
  spec;
- checks unknown keys, malformed supplement width filtering, unsupported
  orientation, and unsupported Cr2 rejection without private stage assertions.

An ignored Be2 timing/proxy script under `tmp/work` is allowed but not
committed. It may report candidate count, residual rank, dimension, elapsed
time, allocations, and artifact write/readback sanity for the z-axis Be2
internal performance proxy. It must not become a gate in this amendment.

## Forbidden

This usability amendment does not approve:

- public export;
- Cr2 full run or Cr2 artifact;
- ECP, EGOI, RHF, solver, or HamV6 export;
- driver/bin/tool workflow;
- report/status/payload object;
- new artifact shape or artifact keys beyond R3-C provenance;
- exposing internal stage objects;
- pair/assembly public workflow;
- parent-stage fields or persistent caches;
- new source file;
- new committed test file.
