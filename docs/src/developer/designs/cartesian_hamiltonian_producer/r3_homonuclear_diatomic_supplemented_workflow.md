# R3 Homonuclear Z-Axis Diatomic Supplemented Workflow

Status: approved narrow molecule-scope relaxation for the residual-GTO/MWG
supplemented facade and canonical driver workflow. This is not element-specific
Cr2 authority and not general molecule support.

## Decision

The R3 usability facade originally supported only z-axis H2 and
internal/performance-supported z-axis Be2. That was correct while the
supplemented path was still being validated, but it now blocks the canonical
driver from acting as a useful molecule driver for explicit homonuclear
diatomic inputs.

Approve a narrow relaxation from hardcoded H2/Be2 guards to explicit
homonuclear two-center z-axis diatomic validation. Cr2 may be used as an
ignored/user-run stress or usability case through this generic path after H2
and Be2 pass; it must not become a Cr2-specific branch, default, fixture, or
committed gate.

## Approved IDs

- `HP-R3U-ZDI-FN-01` - homonuclear z-axis diatomic supplemented facade scope.
- `HP-R3U-ZDI-WIRE-01` - canonical driver supplemented-mode wiring to the
  supported facade.
- `HP-R3U-ZDI-TEST-01` - validation gates for the molecule-scope relaxation.

## Approved Files

Approved source files:

```text
src/cartesian_base_hamiltonian.jl
bin/cartesian_ham_builder.jl
```

`src/cartesian_base_hamiltonian.jl` owns the existing non-exported
`cartesian_residual_gto_mwg_hamiltonian(...)` facade and may replace H/Be
hardcoded guarded helpers with explicit homonuclear z-axis validation.

`bin/cartesian_ham_builder.jl` may wire supplemented mode to that supported
facade through the `HP-DRV-*` compact driver workflow.

No other `src`, `test`, `tools`, `bin`, or committed input-fixture file is
approved by this molecule-scope lane.

## Supported Systems

Approved system scope:

- exactly two centers;
- homonuclear: equal atom symbols and equal finite positive nuclear charges;
- z-axis aligned: both centers have `x == 0`, `y == 0`, and distinct finite
  `z` coordinates;
- neutral all-electron count: `nup` and `ndn` are nonnegative integers and
  `nup + ndn == round(Int, sum(nuclear_charges))`;
- no element-specific defaults, no element-specific branch, and no Cr2-specific
  path.

Required explicit `system` keys remain:

- `atom_symbols`;
- `nuclear_charges`;
- `atom_locations`;
- `nup`;
- `ndn`.

Unsupported systems must throw clear `ArgumentError`s before expensive
construction where practical.

Unsupported without a later amendment:

- heteronuclear systems;
- non-z-axis, translated, rotated, or generally oriented molecules;
- one-center supplemented atoms through this lane;
- ECP inputs;
- charged systems;
- solver/RHF workflow;
- public API/export redesign.

## Basis And Supplement Inputs

The base `basis` input remains explicit. For the diatomic scope it must provide
the existing base-lattice controls such as `q`, `core_spacing`,
`xmax_parallel`, and `xmax_transverse`, with only previously approved defaults
for nonphysical construction controls.

The supplement spec remains explicit and must provide:

- `basis_by_center`;
- `lmax`;
- `uncontracted` when not using the default;
- `width_filtering` when not using the default;
- optional `basisfile`.

Additional supplement validation:

- `basis_by_center` length must equal `2`;
- homonuclear scope requires equal basis labels for both centers;
- optional `basisfile` is either `nothing` or an `AbstractString` naming a
  trusted local/project basis source used by the existing basis-loading path;
- no element-specific basis default is approved;
- no ECP or pseudopotential basis behavior is approved.

## Wiring Contract

`HP-R3U-ZDI-FN-01` may replace hardcoded H/Be checks in the current
`_cartesian_r3_diatomic_inputs(...)` family with explicit homonuclear z-axis
validation. The facade must still construct the base Hamiltonian, terminal
basis, bundles, supplement, residual Gaussian basis, residual MWG interaction,
and optional artifact within the approved same-construction path.

`HP-R3U-ZDI-WIRE-01` may let the canonical driver `:supplemented` mode call the
supported facade instead of composing package-internal helper paths. The driver
must stay within `HP-DRV-*`: visible defaults, optional trusted input file,
command-line overrides, compact timing/summary, artifact write, and optional
readback.

The implementation must not add route objects, status/report/payload fields,
metadata carriers, new artifact schema, public exports, committed fixtures, or
special Cr2 workflow.

## Validation

`HP-R3U-ZDI-TEST-01` approves validation only for this molecule-scope
relaxation:

- H2 supplemented facade/driver artifact path remains unchanged;
- Be2 supplemented facade/driver artifact path remains unchanged and acts as
  the non-H correctness/performance gate;
- optional Cr2 ignored/user-run stress or usability run after H2/Be2 pass;
- no committed Cr2 fixture, committed Cr2 test, or Cr2-specific branch.

Recommended implementation validation:

```text
git diff --check
julia --project=. -e 'using GaussletBases; println("load ok")'
H2 supplemented facade/driver artifact write/readback
Be2 supplemented facade/driver artifact write/readback
optional ignored Cr2 supplemented driver stress/usability run
```

## Line Budget And Failure Rule

Line budget:

- at most `100` added `src`/`bin` lines total;
- net simplification is expected where H/Be-specific guards are replaced by
  explicit homonuclear validation;
- no new committed test, tool, or input-fixture file.

Failure rule: if this requires new route objects, metadata/status/report
fields, artifact schema changes, public exports, committed Cr2 fixtures,
non-z-axis support, heteronuclear support, ECP behavior, solver workflow, or
Cr2-specific branching, stop and request a separate docs-only amendment.
