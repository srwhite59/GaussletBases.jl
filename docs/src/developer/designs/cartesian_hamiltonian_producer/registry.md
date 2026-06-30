# Registry

Only entries marked approved/implemented authorize work on the exact surface
they describe. Measurement-only entries do not authorize production source
edits. Candidate or rejected entries do not authorize implementation.

## Approved And Implemented

### HP-OBJ-01 — `CartesianTerminalBasisBlock`

Owner: `CartesianFinalBasisRealization`.

Exact fields:

```julia
struct CartesianTerminalBasisBlock
    unit_key::Symbol
    support_indices::Vector{Int}
    support_states::Vector{NTuple{3,Int}}
    coefficients::Union{Nothing,Matrix{Float64}}
    column_range::UnitRange{Int}
end
```

`coefficients === nothing` means direct identity on the listed support rows.
Otherwise `coefficients` is support-local rows by retained columns.
`support_indices` and `support_states` are the authoritative owned terminal
rows for the block. They must not mean a post-projection enlarged support, and
PQS shell realization must not grow a block onto previous terminal regions.
Terminal blocks are block-sparse by representation; direct blocks remain
implicit identity on their owned support.

### HP-OBJ-02 — `CartesianTerminalBasisRealization`

Owner: `CartesianFinalBasisRealization`.

Exact fields:

```julia
struct CartesianTerminalBasisRealization
    blocks::Vector{CartesianTerminalBasisBlock}
    final_dimension::Int
    max_cross_overlap::Float64
end
```

The object does not store parent bundles, stage objects, summaries, metadata,
global coefficients, or global self-overlap matrices. Its blocks are
represented on disjoint owned terminal regions. This block-sparse basis
representation does not imply block-diagonal kinetic, nuclear, or IDA
operators. Overlap between distinct owned-support blocks is structurally zero.

The implemented `max_cross_overlap` field is legacy implementation debt after
the structural-support correction. It must not be treated as a physical
cross-block residual or construction repair signal. Source cleanup should
replace production cross-overlap audit plumbing with structural support checks
under a separate implementation handoff.

Under `HP-HAM-MANIFEST-SRC-FN-01`, a source pass may add one optional compact
terminal source-mode provenance carrier if the implementation chooses to attach
the manifest seam to this realization object. That carrier is artifact
provenance only and must not become a basis, operator, route-report, or
algorithm input.

### HP-FILE-01 — terminal realization file

Approved file:

```text
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
```

### HP-FN-00 — block-local terminal shell realization

Approved as a file-local/internal helper under HP-FILE-01. It realizes one PQS
terminal shell by:

```text
full source-box product modes
-> boundary product-mode column selection
-> restrict rows to support.support_indices / support.support_states
-> shell-local Gram on that owned support
-> symmetric shell-local Lowdin
-> final sign canonicalization
-> append block with unchanged owned support
```

Previous-block projection, recursive projection, projection basis, and
effective-support growth are forbidden. Full source-box modes are used only to
generate boundary product-mode columns before row restriction to the
shell-owned support.

### HP-FN-01 — terminal basis realizer

Approved internal surface:

```julia
pqs_terminal_basis_realization(
    support_records,
    retained_records,
    transform_contracts,
    bundles;
    identity_atol = 1.0e-8,
    cross_atol = 1.0e-8,
    weight_atol = 1.0e-14,
)
```

Returns `CartesianTerminalBasisRealization` on success. It owns direct-sector
checks, PQS shell realization, positive final-integral sign canonicalization,
and construction of `CartesianTerminalBasisBlock`.

### HP-FN-02 — structural terminal support checks

Corrected design authority: parent gausslet rows are orthonormal to machine
precision and terminal regions own disjoint parent rows. Therefore block-local
terminal basis supports are structurally orthogonal across blocks, and
cross-block overlap is zero by construction.

Production validation should check:

- every block support equals its authoritative terminal support;
- terminal support sets are pairwise disjoint;
- each shell-local overlap is identity within tolerance.

A nonzero structural overlap means duplicated support rows, incorrect row
restriction, wrong support ownership, or an indexing error. It is not a
physical residual to compute or repair, and production code must not mix
coefficients into previous supports to reduce it. This correction does not
approve source cleanup by itself; removing `max_cross_overlap` and replacing
the audit plumbing needs a separate implementation blurb or approved cleanup
surface.

### HP-WIRE-01 — terminal-basis stage integration

Approved owner: `cartesian_transforms`.

It connects supported PQS terminal plans to
`pqs_terminal_basis_realization(...)` from typed terminal support, retained, and
transform records. It must serve one-center atomic, contact-core diatomic, and
separated diatomic terminal plans through the same entry point.

## Approved For White-Lindsey Terminal Basis Implementation

This section approves only the narrow terminal-basis seam recorded in
`white_lindsey_terminal_basis_realization.md`.

### HP-WLTERM-FILE-01 — optional WL terminal realization file

Approved source files:

```text
src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
```

This ID is optional. It approves creating a small WL-specific terminal
realization sibling and adding its include to `CartesianFinalBasisRealization`
only if extending `pqs_terminal_basis_realization.jl` directly would obscure
the distinct WL boundary-stratum construction. No public export, root include,
new module, new basis object, artifact, report, or status/payload object is
approved.

### HP-WLTERM-FN-01 — WL low-order terminal basis realization

Approved source files:

```text
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl
```

Approved target: realize terminal final-basis blocks for the existing
`:white_lindsey_low_order` route family and return the existing
`CartesianTerminalBasisRealization`.

Allowed behavior:

- support direct identity terminal blocks on their authoritative owned rows;
- realize WL boundary-stratum/product terminal blocks from existing native
  terminal support, retained-rule, and transform records;
- use only support-local coefficients on `support.support_indices` /
  `support.support_states`;
- preserve deterministic terminal support, lowering, retained-record, and
  transform-contract order;
- validate disjoint owned supports and block-local identity overlaps under the
  same structural terminal-basis policy as PQS.

This ID does not approve old WL H1/H1+J materialization adaptation, new
Hamiltonian objects, new route-stage objects, route reports, status/result
payloads, diagnostics, artifact/schema changes, public API/export changes,
raw-block changes, Residual Gaussian/MWG/IDA changes, terminal
shellification-policy changes, retained-selection-policy changes, route
skeleton construction changes, or source files outside the approved surfaces.

Failure rule: if WL boundary-stratum final basis cannot be materialized from
existing terminal lowering, retained-unit, and transform records without
broader route redesign, make no source commit and report the exact missing
native fact.

### HP-WLTERM-WIRE-01 — WL route helper terminal-basis wiring

Approved source file:

```text
src/pqs_source_box_route_driver_helpers.jl
```

Approved behavior:

- remove or narrow the PQS-only terminal-basis guard so
  `route_family = :white_lindsey_low_order` can call the approved WL terminal
  realizer when native terminal records are available;
- keep `route_family = :pqs_source_box` behavior unchanged;
- keep route skeleton construction semantics, route recipe semantics,
  shellification behavior, terminal lowering order, retained-rule order,
  public driver contract, and artifact schema unchanged;
- reject or return a clear unsupported route error when the WL route lacks the
  native facts required for terminal-basis realization.

This ID does not approve adapting old WL materialization, broad route redesign,
supplemented WL behavior, driver input changes, diagnostics/status/report
expansion, raw-block switches, stop-after controls, or source files outside the
approved route helper file.

### HP-WLTERM-TEST-01 — WL terminal-basis validation

Approved validation:

- `git diff --check`;
- package load;
- current default `nesting = :pqs` atom or H2 base artifact/readback remains
  unchanged;
- `nesting = :wl` base atom artifact/readback;
- `nesting = :wl` base H2 artifact/readback if the existing WL diatomic route
  has sufficient native terminal records;
- H2 residual-GTO/MWG PQS endpoint remains unchanged if terminal realization
  code is touched;
- clear unsupported-input/blocker report if WL H2 cannot be realized from
  current native records.

No Cr2 run, supplemented WL run, committed fixture, committed test file,
solver/RHF/ECP/EGOI workflow, artifact schema validation, or broad WL workflow
validation is approved.

## Approved First Composition Lane: WL Z-Axis Diatomic Base

This section promotes the first
`nesting_supplement_composition_plan.md` placeholder. It approves only the
`geometry = z-axis diatomic`, `nesting = :wl`, `supplement = off` base path.

### HP-COMP-WLDIAT-FN-01 — WL z-axis diatomic base terminal records

Approved source files:

```text
src/pqs_source_box_diatomic_complete_core_shell.jl
src/cartesian_terminal_shellification_geometry.jl
src/cartesian_terminal_lowering/selection.jl
src/cartesian_terminal_lowering/region_contracts.jl
src/pqs_source_box_route_driver_helpers.jl
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
src/cartesian_base_hamiltonian.jl
```

Approved goal:

```text
Natom = 2
nesting = :wl
basisname = nothing
```

must produce an existing `CartesianIDAHamiltonian{Float64}` artifact/readback
through the same `CartesianTerminalBasisRealization` and staged base
Hamiltonian path used by the PQS producer.

Allowed behavior:

- produce native White-Lindsey z-axis diatomic terminal support,
  shellification, lowering, retained-rule, and transform records needed by the
  existing WL terminal realizer;
- preserve the existing `:white_lindsey_low_order` construction-family route
  and deterministic support/lowering/retained/transform order;
- realize WL boundary-stratum/product terminal blocks as owned-support
  terminal blocks in the existing `CartesianTerminalBasisRealization`;
- reuse the existing staged product, unit-nuclear, IDA, Hamiltonian
  construction, writer, and reader path;
- in `src/cartesian_base_hamiltonian.jl`, add only narrow staged/facade wiring
  required by WL z-axis diatomic base construction and the truthful route
  provenance value `:z_axis_diatomic_wl_base`.

The `:z_axis_diatomic_wl_base` route value is approved as a value under the
existing `producer_provenance/route` and `recipe_provenance/route` keys. It is
not an artifact schema change.

Forbidden:

- driver public input changes or driver special cases;
- old WL H1/H1+J materialization revival or adaptation;
- artifact schema changes, matrix-key changes, reader behavior changes,
  manifest shape changes, public API/export changes, or new Hamiltonian
  wrapper/result objects;
- Residual Gaussian, MWG/IDA, supplement, ECP, solver/RHF, or Cr2 workflow
  work;
- route diagnostics, stop-after controls, report/status/payload fields,
  raw-block switches, route-stage labels, or broad route skeleton redesign;
- committed tests, committed fixtures, committed driver input files, or source
  files outside the approved list.

Failure rule: if WL z-axis diatomic base construction requires adapting the old
WL H1/H1+J materialization path, changing artifact schema/reader behavior,
adding driver special cases, or creating a parallel Hamiltonian builder, make
no source commit and report the blocker.

Line budget: at most `250` added `src` lines, with deletion or simplification
of obsolete blocker-only WL diatomic guards expected where practical. Stop for
a new amendment if the pass needs broader route skeleton redesign, source
files outside the approved list, or persistent payload/cache objects.

### HP-COMP-WLDIAT-TEST-01 — WL z-axis diatomic base validation

Approved validation:

- `git diff --check`;
- package load;
- current default `nesting = :pqs` H2 base artifact/readback remains
  unchanged;
- existing `nesting = :wl` one-center atom artifact/readback remains
  unchanged;
- `nesting = :wl` z-axis H2 base artifact/readback succeeds through the staged
  base Hamiltonian path;
- direct provenance inspection confirms `nesting = :wl` and
  `route = :z_axis_diatomic_wl_base` for the WL H2 artifact;
- H2 residual-GTO/MWG PQS endpoint remains unchanged if terminal realization
  code is touched;
- no Cr2 run.

No supplemented WL run, committed test file, committed fixture, driver
contract test, solver/RHF/ECP/EGOI validation, route-diagnostic validation, or
Cr2 fixture is approved.

## Approved Correction Lane: WL Z-Axis Diatomic Compact Retained Basis

This section records the follow-up design decision after the WL diatomic
terminal-record endpoint exposed the remaining placeholder-like retained-basis
shape. The current `nesting = :wl`, `Natom = 2` path can be mechanically
realized, but it still follows:

```text
elongated shared complete shell
-> boundary CPB strata
-> one retained identity unit per stratum
-> identity terminal blocks
```

That is not the intended compact White-Lindsey retained basis and should not be
used as the production PQS/WL comparison story. The observed audit evidence was
that bounded WL diatomic `ns = 4/5` examples built an elongated shell with
support-size scale `(5,5,9) - (3,3,7) = 162`, rather than the cubic shell-size
scales `4^3 - 2^3 = 56` and `5^3 - 3^3 = 98`; WL retained-unit lowering then
split that shell into 26 boundary-stratum units that the terminal realizer
appended as full-support identity blocks. For `ns = 6`, contact-core geometry
can consume the bounded support and collapse the terminal basis to one direct
identity block.

### HP-WLDIAT-COMPACT-FN-01 — WL diatomic compact retained basis

Approved source files:

```text
src/cartesian_shellification/terminal_geometry.jl
src/cartesian_terminal_lowering/region_contracts.jl
src/cartesian_retained_units/lower_contract_units.jl
src/cartesian_retained_unit_transform_contracts/unit_contracts.jl
src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl
src/pqs_source_box_route_driver_helpers.jl
```

`src/pqs_source_box_route_driver_helpers.jl` is approved only for narrow route
wiring if the compact WL retained-unit facts must be passed to the existing WL
terminal-basis seam.

Approved behavior:

- preserve the WL unit-based implementation model: faces, edges, corners, and
  small boundary units after shellification;
- do not force a persistent shell object after retained-unit splitting;
- make each WL unit carry or realize the intended compact retained basis from
  products of one-dimensional contractions on the authoritative owned unit
  support;
- treat identity realization as valid only for true direct/core identity
  units, not for WL boundary-stratum retained units;
- use deleted WL coefficient helpers only as historical donor/reference
  material for the compact CPB-local product-of-1D coefficient primitive;
- preserve deterministic geometry, lowering, retained-unit, transform-contract,
  and terminal-block ordering;
- keep the same public `ns` as the fair starting input for PQS/WL comparison,
  while allowing WL-specific geometry/contact cases and not promising equal
  final dimensions.

Forbidden:

- driver changes;
- artifact schema, provenance, matrix-key, reader, or manifest changes;
- PQS behavior changes;
- Hamiltonian assembly changes;
- raw-block, Residual Gaussian, MWG/IDA, Qiu-White, supplement, solver/ECP,
  or Cr2 workflow changes;
- old WL route-global stack, reports, adapters, or H1/H1+J materialization
  revival or adaptation;
- broad route diagnostics, report/status/payload fields, raw-block switches,
  retained-rule dumps, or route-stage labels;
- fake compactness by dropping support rows, relabeling full-support identity
  units, or changing the driver comparison;
- committed tests, committed fixtures, or committed driver input files.

Failure rule: if compact WL retained units require construction-native facts
that are not currently available, make no source commit and report the exact
missing fact. Do not fake compactness by deleting rows, changing public input
semantics, or rerouting through old WL materialization. If an essential
primitive exists only in deleted WL files, restore or re-express only that
primitive behind the current terminal-basis boundary; do not restore the old
route-global framework around it.

Line budget: at most `250` added `src` lines unless a later source blurb
narrows or revises the budget after auditing the exact live callers.

### HP-WLDIAT-COMPACT-TEST-01 — WL compact-basis validation

Approved validation:

- `git diff --check`;
- package load;
- small H2 or Be2 WL base artifact/readback;
- small H2 or Be2 WL supplemented artifact/readback if the compact base path
  works through the existing supplemented boundary;
- PQS base and supplemented smokes remain unchanged;
- WL retained dimension is compared against expected shell-size scale for
  bounded `ns = 4/5` examples;
- finite/symmetric `K` and `V` checks;
- no Cr2 run.

No committed test file, committed fixture, driver contract test,
solver/RHF/ECP/EGOI validation, route-diagnostic validation, or Cr2 fixture is
approved.

## Approved Correction Lane: WL Boundary-Stratum Retained-Count Parity

This section records the follow-up policy correction after the compact WL
diatomic terminal-basis source pass. The full-support identity bug is fixed,
but `nesting = :wl`, `ns = 4` still follows the inherited symmetric-odd donor
rule and produces 26 boundary columns rather than the nominal shell count
`4^3 - 2^3 = 56`.

That remaining behavior is not an acceptable WL policy. The odd-side rule is a
direct core-block centering requirement: a nucleus-centered direct core should
have odd side length so the nucleus is centered. Boundary shells and boundary
strata outside that core do not require odd side counts and must retain the
requested shell contraction count.

### HP-WLDIAT-PARITY-FN-01 — WL boundary retained-count parity

Approved source file:

```text
src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl
```

Approved behavior:

- preserve odd-side enforcement for true direct nucleus-centered core blocks;
- for WL boundary shell strata, use the requested boundary retained count
  without symmetric-odd coercion;
- for public `nesting = :wl`, `ns = 4`, route-local `q = 2` must retain the
  shell count `4^3 - 2^3 = 56`;
- for public `nesting = :wl`, `ns = 5`, retain `5^3 - 3^3 = 98`;
- keep the compact WL product-of-1D coefficient construction and existing
  terminal-basis boundary.

Forbidden:

- driver changes;
- public `ns` normalization or route-local `q` rule changes;
- route skeleton, shellification, terminal lowering, retained-unit metadata
  shape, or contract-plan changes;
- direct/core identity behavior changes;
- artifact schema/provenance, matrix-key, reader, or manifest changes;
- PQS behavior changes;
- Hamiltonian assembly changes;
- raw-block, Residual Gaussian, MWG/IDA, Qiu-White, supplement, solver/ECP,
  or Cr2 workflow changes;
- old WL route-global stack, reports, adapters, or H1/H1+J materialization
  revival or adaptation;
- broad route diagnostics, report/status/payload fields, raw-block switches,
  retained-rule dumps, or route-stage labels;
- committed tests, committed fixtures, or committed driver input files.

Failure rule: if fixing boundary parity requires changing source files outside
`src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`,
public `ns` semantics, route/shellification/terminal-lowering contracts,
artifact schema, or old WL materialization, make no source commit and report
the exact blocker.

Line budget: target under `30` added `src` lines, with no new persistent shape.

### HP-WLDIAT-PARITY-TEST-01 — WL parity validation

Approved validation:

- `git diff --check`;
- package load;
- WL H2 or Be2 z-axis diatomic `ns = 4` retained boundary count / dimension
  demonstrates 56 boundary columns rather than 26;
- WL H2 or Be2 z-axis diatomic `ns = 5` retained boundary count demonstrates
  98 boundary columns;
- small WL base artifact/readback smoke;
- small WL supplemented artifact/readback smoke if bounded by the existing
  supplemented boundary;
- PQS H2 residual-GTO/MWG endpoint remains unchanged;
- finite/symmetric `K` and `V` checks for the WL smoke;
- no Cr2 run.

No committed test file, committed fixture, driver contract test,
solver/RHF/ECP/EGOI validation, route-diagnostic validation, or Cr2 fixture is
approved.

## Approved Composition Cleanup: Public ns Direct-Core Side Parity

This lane corrects an even-`ns` PQS/WL comparison bug. Public `ns` is the
requested cube/source/nesting size, but shared route setup still derives the
nucleus-centered direct/core side from route-local `q`. Since WL derives
`q = ns - 2`, even-`ns` same-size PQS/WL comparisons can get different direct
core boxes even though the direct core centering rule is not a route-family
physics difference.

The odd-side parity rule is necessary only for direct nucleus-centered core
identity blocks. It must not apply to boundary shells, WL boundary-stratum
retained products, or non-direct support regions.

### HP-COMP-NSCORE-FN-01 — public ns direct core side

Approved source files:

```text
src/pqs_source_box_route_driver_helpers.jl
src/cartesian_base_hamiltonian.jl
```

`src/cartesian_base_hamiltonian.jl` is approved only if needed to keep
one-center parent minimum sizing consistent with the same direct-core rule.

Approved behavior:

- keep public `ns` as the requested cube/source/nesting size;
- keep route-local `q` derivation unchanged:
  - PQS: `q = ns`;
  - WL: `q = ns - 2`;
- derive only direct nucleus-centered core side from public `ns`:
  `direct_core_side = isodd(ns) ? ns : ns + 1`;
- do not apply this oddization rule to boundary shells, WL boundary-stratum
  retained products, or non-direct support regions;
- preserve odd-`ns` dimensions where this cleanup should be a no-op;
- preserve WL boundary retained count policy, including `ns = 4 -> 56`,
  `ns = 5 -> 98`, and `ns = 6 -> 152`;
- update internal provenance or summary labels such as `:odd_q_core_side` to
  an `ns`-truthful rule if they are still written.

Forbidden:

- driver changes;
- public input changes;
- route skeleton redesign;
- shellification geometry rewrites beyond direct-core side authority;
- terminal lowering, retained-unit, or terminal-realizer changes;
- artifact schema changes;
- manifest expansion;
- PQS/WL special cases in the driver;
- old WL materialization revival;
- committed tests or fixtures;
- Cr2-specific workflow.

Failure rule: if fixing the parity requires changing route skeleton semantics,
terminal lowering, retained-unit records, WL boundary coefficient construction,
artifact schema, or driver inputs, make no source commit and report the exact
blocker.

### HP-COMP-NSCORE-TEST-01 — public ns direct-core parity validation

Approved validation:

- `git diff --check`;
- package load;
- small one-center atom base artifact/readback for `ns = 5`, `6`, and `7` with
  `nesting = :pqs` and `nesting = :wl`;
- same-`ns` PQS/WL atom dimensions match for `ns = 5`, `6`, and `7` under a
  bounded common fixture;
- `ns = 6` no longer has the reported `66`-row skew;
- small H2 or Be2 smoke to confirm the diatomic path still constructs;
- existing H2 Residual Gaussian endpoint smoke if touched code crosses the
  supplemented path;
- no Cr2 run.

No committed test file, committed fixture, driver contract test,
solver/RHF/ECP/EGOI validation, route-diagnostic validation, artifact schema
validation, or Cr2 fixture is approved.

## Approved Composition Cleanup: Common Terminal Shell Decomposition

This lane freezes the two-geometry distinction needed to keep PQS/WL
composition from drifting.

Common terminal shell decomposition is the first geometry:

```text
parent lattice + nuclear centers + direct core side
-> direct core regions
-> terminal shell regions
-> owned support rows
```

PQS and White-Lindsey must share this first step. The family split begins only
after common shell records exist:

```text
PQS: common shell support + full source CPB -> retained source-box modes
WL:  common shell support -> faces/edges/corners/strata -> 1D contractions
```

### HP-COMP-SHELLGEOM-FN-01 — common shell decomposition audit/cleanup

Approved source files:

```text
src/cartesian_shellification/terminal_geometry.jl
src/pqs_source_box_route_driver_helpers.jl
```

Approved behavior:

- audit whether first-step shell/core region construction is already identical
  for `nesting = :pqs` and `nesting = :wl`;
- keep direct core regions, shell regions, owned support rows, ordering, and
  coverage route-family-free;
- remove route-family branching from common terminal shell decomposition if it
  exists;
- rename or locally clarify shellifier parameters and summary labels so common
  shell geometry is not presented as governed by PQS `q` or WL inner side;
- keep direct-core side tied to public `ns` through `HP-COMP-NSCORE-*`;
- leave PQS and White-Lindsey lowering/realization separate after common shell
  records are produced.

Approved source surface details:

- `src/cartesian_shellification/terminal_geometry.jl` owns common
  route-family-free shell/core region decomposition;
- `src/pqs_source_box_route_driver_helpers.jl` may change only narrow caller
  plumbing and summary/provenance wording needed to pass common shell inputs
  before selecting PQS or White-Lindsey lowering.

For one-center atoms, same public system, parent extent, and `ns` must produce
the same direct core and shell-owned support regions before family-specific
lowering.

For z-axis diatomics, this lane may audit whether central-gap and shared-shell
planning use the same common shell decomposition. It must not change central
gap/contact policy unless the fix is the same route-family-free shell input
cleanup and does not touch lowering, retained units, or WL/PQS realization.

### HP-COMP-SHELLGEOM-DIAT-FN-01 — diatomic common shellifier entry

Approved source files:

```text
src/cartesian_shellification/terminal_geometry.jl
src/pqs_source_box_route_driver_helpers.jl
```

Approved behavior:

- for fixed public z-axis diatomic system, parent axes, public `ns`, direct
  core side, nuclear centers, and bond axis, `nesting = :pqs` and
  `nesting = :wl` must call the same common shellifier with the same
  first-step arguments;
- central-gap, contact-core, shared-shell, and outer-mismatch region ownership
  are common shell decomposition facts;
- PQS `q` and WL inner side are retained-construction inputs after common
  shell records exist and must not change common diatomic shellifier entry;
- if `raw_terminal_geometry(...)` currently uses a parameter named `q` for
  common central-gap or shared-shell decisions, source work may rename or
  reinterpret that shellifier boundary parameter as common `ns`;
- keep the central-gap/contact algorithm unchanged except for replacing
  route-family-dependent inputs with the shared common inputs.

This does not approve changing terminal lowering, retained units, PQS
source-box realization, WL face/edge/corner coefficient construction, route
skeletons, artifacts, driver inputs, or public API.

Failure rule: if same-function/same-argument z-axis diatomic entry requires
changing the central-gap/contact algorithm rather than only its
route-family-independent inputs, make no source commit and request a separate
docs-only amendment.

### HP-COMP-SHELLGEOM-DIAT-TEST-01 — diatomic shellifier-entry validation

Approved validation:

- `git diff --check`;
- package load;
- focused audit showing z-axis diatomic PQS/WL calls enter the common
  shellifier with the same parent axes, nuclear centers, direct core side,
  public `ns`, and bond axis;
- same-`ns` PQS/WL z-axis diatomic direct core, central-gap/contact,
  shared-shell, and outer-mismatch region counts match before
  family-specific lowering for a bounded H2 or Be2 fixture;
- base artifact/readback smoke for the same bounded fixture under
  `nesting = :pqs` and `nesting = :wl` if the WL retained-basis path remains
  available;
- existing H2 Residual Gaussian endpoint smoke only if touched code crosses
  supplemented path;
- no Cr2 run.

No committed test file, committed fixture, driver contract test,
solver/RHF/ECP/EGOI validation, route-diagnostic validation, artifact schema
validation, or Cr2 fixture is approved.

### HP-COMP-OUTERMM-FN-01 / HP-COMP-OUTERMM-TEST-01 — superseded subset

Status: superseded by `HP-COMP-THINSLAB-FN-01` and
`HP-COMP-THINSLAB-TEST-01`.

The outer-mismatch-only lane correctly identified that identity lowering was
wrong, but the same direct-identity mistake also applies to central midpoint
slabs. Do not implement `HP-COMP-OUTERMM-*` as a separate source lane.

### HP-COMP-THINSLAB-FN-01 — common thin-slab stack compact lowering

Status: approved.

Approved source files:

```text
src/cartesian_terminal_lowering/selection.jl
src/cartesian_terminal_lowering/region_contracts.jl
src/cartesian_retained_units/lower_contract_units.jl
src/cartesian_retained_unit_transform_contracts/unit_contracts.jl
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl
```

Optional only if directly required for native slab-axis/thickness metadata or
existing summaries/records:

```text
src/cartesian_shellification/terminal_geometry.jl
src/pqs_source_box_route_driver_helpers.jl
src/pqs_source_box_diatomic_complete_core_shell.jl
```

Problem:

- `src/cartesian_shellification/terminal_geometry.jl` creates
  `region_kind = :direct_midpoint_slab` for central-gap one-slice slabs and
  `region_kind = :outer_mismatch_slab` for boundary mismatch slabs;
- `src/cartesian_terminal_lowering/region_contracts.jl` currently maps
  midpoint slabs to `:direct_slab_identity_cpb`;
- `src/cartesian_terminal_lowering/region_contracts.jl` currently maps that
  outer-mismatch region to `:direct_boundary_slab_identity_cpb`;
- `src/cartesian_retained_units/lower_contract_units.jl` then treats the
  regions as direct retained units;
- terminal realization therefore writes a full identity block.

For z-axis diatomics, thin slabs are not production identity sectors for either
PQS or White-Lindsey. They are face-like compact slab objects.
Outer-mismatch evidence from CR2 found `z_low_outer_mismatch_slab` and
`z_high_outer_mismatch_slab` contributing `7605` final rows each, or `15210`
rows total. That is a producer basis-size bug, not a CR2/HFDMRG consumer
issue. The central midpoint slab has the same lowering category: a one-slice
slab between atom-local boxes, not a direct core.

Approved behavior:

- for `PQSLowering` and `WhiteLindseyLowering`, `:direct_midpoint_slab` and
  `:outer_mismatch_slab` must not lower to `:direct_slab_identity_cpb` or
  `:direct_boundary_slab_identity_cpb`;
- use the same compact slab lowering function for both construction families,
  with the same terminal region, public `ns`, slab normal axis, slab
  thickness, and native source/support facts, when those facts are sufficient;
- keep real shell-region lowering route-specific after common shellification:
  PQS full source-box shell projection and WL face/edge/corner product
  contractions remain different constructions;
- keep direct nucleus-centered core and atom-contact core sectors identity
  sectors;
- the compact unit slice should have retained scale `ns x ns x 1` after
  standard one-dimensional COMX/product compression, not full identity support
  rows;
- an outer-mismatch region of thickness `t <= ns` should be decomposed or
  realized as an oriented stack of compact slices with scale about
  `t * ns * ns`, not as one identity block;
- if slab thickness exceeds `ns`, source work must stop and report whether a
  whole-block `ns x ns x ns` compression or a setup-error policy needs
  separate approval;
- preserve shellification coverage, owned support disjointness, deterministic
  ordering, common PQS/WL first-step geometry, Hamiltonian assembly, artifacts,
  driver inputs, Residual Gaussian, MWG/IDA, and reader behavior;
- keep PQS and White-Lindsey thin-slab lowering identical at this boundary.
  Route-family-specific thin-slab lowering requires a later amendment.

Forbidden:

- driver changes;
- artifact, manifest, provenance, schema, or reader changes;
- residual Gaussian, MWG, IDA, Hamiltonian assembly, raw-block, or solver
  changes;
- route-family-specific thin-slab behavior;
- route skeleton redesign;
- broad terminal realization redesign beyond the shared thin-slab consumer;
- retained-unit record changes beyond the shared thin-slab retained object;
- retained-unit transform changes beyond the shared thin-slab transform
  contract;
- old high-order workflow revival;
- committed Cr2 tests or fixtures;
- direct slab deletion unless a separate approval proves the slabs are
  nonphysical padding regions.

Failure rule: if the slab cannot be compacted through existing compact slab
machinery without changing shellification semantics, route skeletons, artifact
schema, driver inputs, or real-shell policy, make no source commit and report
the exact missing native fact. Do not parse region role names to infer slab
normal or thickness; add native shellification metadata under the optional
source surface if needed. If slab thickness exceeds `ns`, do not install a
silent fallback under this ID.

### HP-COMP-THINSLAB-TEST-01 — common thin-slab stack validation

Status: approved.

Approved validation:

- `git diff --check`;
- package load;
- bounded H2 or Be2 artifact/readback under `nesting = :pqs` and
  `nesting = :wl` where midpoint slabs and outer-mismatch slabs are present or
  explicitly probed;
- confirm no `direct_slab_identity_cpb` rows for `:direct_midpoint_slab`
  regions under either lowering family;
- confirm no `direct_boundary_slab_identity_cpb` rows for
  `:outer_mismatch_slab` regions under either lowering family;
- focused audit showing PQS and White-Lindsey call the same compact
  thin-slab lowering function with the same region/public-`ns` inputs for
  matched slab regions;
- confirm retained thin-slab count is compact relative to support count, with
  normal oriented unit-slice target `ns x ns x 1` and thickness-`t <= ns`
  outer-mismatch scale about `t * ns * ns`;
- existing H2 Residual Gaussian endpoint smoke if touched code crosses
  supplemented construction;
- optional CR2 user-run inventory only, not a committed gate.

No committed test file, committed fixture, driver contract test,
solver/RHF/ECP/EGOI validation, route-diagnostic validation, artifact schema
validation, or Cr2 fixture is approved.

Forbidden:

- driver changes;
- public input changes;
- route skeleton redesign;
- terminal lowering redesign beyond the narrow common thin-slab repair in
  `HP-COMP-THINSLAB-FN-01`;
- retained-unit record changes beyond the shared thin-slab retained object;
- retained-unit transform changes beyond the shared thin-slab transform
  contract;
- PQS source-box retained-mode realization changes;
- WL face/edge/corner coefficient or retained-basis changes;
- direct-core parity changes beyond `HP-COMP-NSCORE-*`;
- central-gap/contact policy redesign;
- artifact schema, manifest, reader, or provenance expansion;
- Hamiltonian, one-body, IDA, MWG, Residual Gaussian, raw-block, solver, ECP,
  or Cr2 workflow changes;
- old WL materialization revival;
- committed tests or fixtures.

Failure rule: if common shell decomposition cannot be made route-family-free
without changing terminal lowering, retained-unit records, PQS retained-mode
realization, WL boundary coefficient construction, route skeleton semantics,
artifact schema, or driver inputs, make no source commit and report the exact
blocker.

### HP-COMP-THINSLAB-META-FN-01 — thin-slab metadata inventory cleanup

Status: approved.

Approved source file:

```text
src/cartesian_terminal_shellification_geometry.jl
```

Only in support of the already approved thin-slab lowering pass:

```text
src/pqs_source_box_route_driver_helpers.jl
src/pqs_source_box_diatomic_complete_core_shell.jl
```

Problem:

`src/cartesian_terminal_shellification_geometry.jl` is not the new
shellifier owner, but it is still a live metadata/scaffold inventory path.
The route helper calls
`_cartesian_terminal_shellification_region_unit_inventory(...)`, and
`_cartesian_terminal_region_unit_mapping(region)` still encodes the old direct
identity slab contract:

```text
:direct_midpoint_slab -> :direct_slab_identity_cpb
:outer_mismatch_slab -> :direct_boundary_slab_identity_cpb
```

It also has no planned `:angular_z_extension_slab` case. This makes route
inventory disagree with the approved compact thin-slab lowering contract.

Approved behavior:

- update `_cartesian_terminal_region_unit_mapping(...)` so midpoint slabs,
  outer-mismatch fallback slabs, and angular z-extension slabs map to the
  compact thin-slab lowering category, not direct identity categories;
- add only the minimal compact thin-slab inventory/count vocabulary needed for
  existing route summaries to agree with terminal lowering;
- recognize `:angular_z_extension_slab` as a planned compact thin-slab region;
- preserve the metadata-only nature of the file: it describes planned
  lowering consistently, but does not materialize coefficients or construct
  Hamiltonian data;
- keep direct core and atom-contact core identity mappings unchanged.

Forbidden:

- driver changes;
- artifact, manifest, provenance, schema, or reader changes;
- shellification algorithm changes;
- route skeleton redesign;
- Residual Gaussian, MWG, IDA, Hamiltonian, raw-block, or solver changes;
- committed tests or fixtures;
- Cr2 workflow;
- new reporting framework;
- reintroduction of direct identity slab lowering under a new name;
- broad cleanup or deletion of the terminal-shellification metadata file in
  this pass.

Failure rule: if updating this metadata inventory requires materializing
retained units, adding artifact/report payloads, changing shellification
geometry, or redesigning route skeletons, make no source commit and report the
blocker. This ID is only for keeping the live metadata/scaffold inventory
consistent with compact thin-slab lowering.

### HP-COMP-THINSLAB-META-TEST-01 — thin-slab metadata inventory validation

Status: approved.

Approved validation:

- `git diff --check`;
- package load;
- existing angular geometry audit still passes;
- thin-slab inventory/probe no longer blocks on unsupported terminal
  shellification region kind `angular_z_extension_slab`;
- focused scan confirms no planned direct identity lowering for midpoint,
  outer-mismatch, or angular z-extension slabs in this metadata inventory;
- bounded H2/Be2 artifact/readback validation remains under the main
  `HP-COMP-THINSLAB-*` implementation pass;
- no Cr2 run.

No committed test file, committed fixture, public driver test, artifact schema
test, or Cr2 fixture is approved.

### HP-COMP-FACEPROD-FN-01 — neutral terminal face-product helper

Status: approved.

Approved source files:

```text
src/cartesian_final_basis_realization/terminal_face_product_blocks.jl
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
```

Problem:

`HP-COMP-THINSLAB-*` needs a compact product-block seam for slabs, but the
first source attempt placed a reusable helper in
`white_lindsey_terminal_basis_realization.jl` and then called it from PQS.
That makes shared numerical work look White-Lindsey-owned. The correct owner is
a neutral internal helper under `CartesianFinalBasisRealization`.

Approved behavior:

- add a private/module-internal neutral helper for compact face-product
  terminal blocks;
- reuse `_nested_doside_1d(...)` and `_nested_face_product(...)`;
- support normal axes `:x`, `:y`, and `:z`;
- support one fixed normal-axis index for a face-like block;
- support an ordered stack of fixed normal-axis indices for a thickness-`t`
  slab;
- take the retained count from the caller, normally public `ns`;
- refactor White-Lindsey facet terminal realization to use the neutral helper;
- allow later `HP-COMP-THINSLAB-*` implementation to use the same helper for
  midpoint, outer-mismatch fallback, and angular z-extension thin slabs;
- preserve existing support validation, overlap identity validation,
  deterministic ordering, and owned-support disjointness.

This helper seam does not replace `HP-COMP-THINSLAB-*`. It only approves the
shared face-product coefficient assembly needed so thin-slab lowering can be
implemented through reuse instead of duplicated PQS/WL terminal code.

Forbidden:

- driver changes;
- public API/export;
- artifact, manifest, provenance, schema, or reader changes;
- shellification algorithm changes;
- terminal lowering policy changes by itself beyond enabling the later
  `HP-COMP-THINSLAB-*` pass;
- route skeleton changes;
- Residual Gaussian, MWG, IDA, Hamiltonian, raw-block, or solver changes;
- old high-order workflow revival;
- committed tests or fixtures;
- Cr2 workflow;
- duplicate implementation of face-product coefficient assembly;
- PQS-specific thin-slab projection path;
- treating thin slabs as White-Lindsey boundary strata for naming convenience.

Line budget: target at most `80` added source lines for the neutral helper and
WL facet refactor. If the helper needs substantially more, stop and report the
missing abstraction before continuing.

Failure rule: if a neutral helper cannot serve both current White-Lindsey
facets and future thin slabs without changing numerical semantics, make no
source commit and report whether the blocker is helper signature,
support-record shape, retained-unit metadata, or terminal-realization
ownership.

### HP-COMP-FACEPROD-TEST-01 — neutral terminal face-product validation

Status: approved.

Approved validation:

- `git diff --check`;
- package load;
- focused ignored helper probe if needed;
- White-Lindsey facet parity before/after the refactor, showing identical
  dimensions and coefficients or roundoff agreement;
- H2 or Be2 base artifact/readback for `nesting = :wl` to confirm WL facet
  behavior is unchanged;
- H2 or Be2 base artifact/readback for `nesting = :pqs` if PQS terminal
  realization imports or consumes the helper;
- no Cr2 run.

No committed test file, committed fixture, public driver test, artifact schema
test, or Cr2 fixture is approved.

### HP-COMP-ANGBOX-FN-01 — angular-balanced z-axis diatomic shellification

Status: approved.

Approved source file:

```text
src/cartesian_shellification/terminal_geometry.jl
```

Optional only if directly required for existing summary/caller plumbing:

```text
src/pqs_source_box_route_driver_helpers.jl
src/pqs_source_box_diatomic_complete_core_shell.jl
```

Approved behavior:

- each shared z-axis diatomic molecular shellification step should compute an
  angular-balanced target box from the outer nuclei in physical parent-axis
  coordinates;
- compare the physical bond-axis longitudinal margin from each outer nucleus
  to the selected transverse physical scale; if `x` and `y` transverse scales
  differ, use the smaller scale unless a later amendment approves another
  convention;
- treat this as the operational `:outer_nucleus_45_degree` shellification
  rule;
- when the angular-balanced target requires bond-axis-only extension beyond
  the ordinary index-layer shell body, emit that difference as planned axial
  thin-slab stacks with native metadata. The ordinary body plus planned
  z-extension slabs, not the ordinary body alone, realizes the target
  coverage;
- specifically, when ordinary shared-shell expansion stops with transverse
  axes saturated and bond-axis parent support remaining, emit the bond-axis
  leftovers as planned `:angular_z_extension_slab` stack regions;
- do not treat planned z-extension slabs as route-family-specific shell
  regions or direct identity sectors;
- apply the same thin-slab category to central midpoint slabs, planned
  non-boundary angular z-extension slabs, planned boundary angular z-extension
  slabs, and
  unexpected outer-mismatch fallback slabs;
- planned z-extension stacks with total thickness greater than `ns` should be
  split into multiple ordered compact slab units with thickness `<= ns`;
- unexpected fallback slabs with thickness greater than `ns` remain a
  setup/shellification failure unless a later policy approves a whole-block
  compression.

Planned angular z-extension metadata should include:

```text
slab_kind = :angular_z_extension_slab
slab_normal_axis
slab_side
slab_thickness
slab_stack_index
slab_stack_count
bond_axis
reference_nucleus_index
angular_balance_rule = :outer_nucleus_45_degree
longitudinal_margin_physical
transverse_scale_physical
angular_extension_physical
```

This lane must not parse region labels to infer slab geometry, must not change
real shell retained policy, and must not turn axial z-extension slabs into
identity sectors. Thin-slab lowering remains under `HP-COMP-THINSLAB-*`.

Forbidden:

- driver changes;
- public input changes;
- artifact, manifest, provenance, schema, or reader changes;
- terminal lowering, retained-unit, transform-contract, or terminal-realization
  changes;
- Hamiltonian, one-body, IDA, MWG, Residual Gaussian, raw-block, solver, ECP,
  or Cr2 workflow changes;
- route skeleton redesign;
- route-family-specific PQS/WL geometry;
- PQS/WL real-shell retained-policy changes;
- direct slab deletion;
- Cr2-specific branches;
- committed tests or fixtures.

Failure rule: if angular-balanced shellification
requires changing terminal lowering, retained-unit records, route skeletons,
artifact schema, driver inputs, or real-shell retained policy beyond emitting
native thin-slab stack regions, do not commit source work and report the exact
missing native fact or policy.

### HP-COMP-ANGBOX-TEST-01 — angular-balanced shellification validation

Status: approved.

Approved validation:

- `git diff --check`;
- package load;
- ignored geometry audit output for bounded H2/Be2 and user-run Cr2-style
  fixtures showing parent physical endpoints/counts, snapped nuclear indices,
  core boxes, molecular inner box, each proposed shared-shell expansion,
  transverse scale, low/high longitudinal margins, angular-balance ratios,
  planned z-extension slab stacks, and residual outer mismatch if any;
- ignored angular geometry audit showing planned z-extension support, zero
  residual z mismatch after classification, and PQS/WL geometry parity;
- prove planned z-extension stack slices use the same thin-slab lowering
  category for PQS and WL;
- no artifact/readback is required while lowering is intentionally deferred;
- no committed Cr2 fixture or Cr2 run.

Separate existing `HP-COMP-SHELLGEOM-DIAT-*` note:
central-gap/contact algorithm redesign remains deferred. That lane approves
only shared shellifier entry and route-family-independent inputs to the current
algorithm. `HP-COMP-ANGBOX-*` does not change that central-gap/contact policy.

### HP-COMP-SHELLGEOM-TEST-01 — common shell decomposition validation

Approved validation:

- `git diff --check`;
- package load;
- focused audit showing the one-center atom common shell decomposition is
  route-family-free for the same public `ns`, parent extent, and center;
- same-`ns` PQS/WL one-center atom direct core and shell-owned support counts
  match before family-specific lowering;
- same-`ns` PQS/WL one-center atom base artifact/readback still works for a
  bounded fixture;
- H2 or Be2 smoke to confirm the diatomic path still constructs;
- existing H2 Residual Gaussian endpoint smoke only if touched code crosses
  supplemented path;
- no Cr2 run.

No committed test file, committed fixture, driver contract test,
solver/RHF/ECP/EGOI validation, route-diagnostic validation, artifact schema
validation, or Cr2 fixture is approved.

## Approved Source-Span Facility: Mapped-COMX

This section approves a narrow mainline source-span option based on
high-order-manager scratch evidence. The high-order branch is an experimental
proving ground and benchmark consumer. It does not own the mainline
implementation shape, and its scripts, route wrappers, reports, and benchmark
scaffolding must not be copied into production.

Source-review correction: mapped-COMX is a new source-span option inside the
existing nested doside / COMX path, not a new numerical facility under
`CartesianRawProductSources`.

Approved implementation path:

```text
pqs_source_axis_transform_facts_from_pgdg_axes(...)
-> _nested_doside_1d(...)
-> _nested_retained_span(...)
-> _cleanup_comx_transform(...)
```

The option changes only the raw span passed into the existing physical COMX
cleanup.

Approved first numerical rule:

```text
protected physical P2
+ mapped Chebyshev enrichment T_k(s_lambda(u))
+ lambda = 0.5
+ no sqrtJ
+ physical-u COMX localization
```

with:

```text
s_lambda(u) = (1 + lambda) * u / (1 + lambda * u^2)
```

Here `u` is the dimensionless local coordinate on the source interval,
`u = (x - x_mid) / x_half`, not raw physical center `x`. The existing COMX
cleanup still uses the physical position matrix. Applying `s_lambda` directly
to physical centers is forbidden.

The ordinary polynomial source-span path remains available and unchanged. This
lane installs an additional source-span facility; it does not change defaults
or make high-order benchmark results production acceptance tests.

Post-installation evidence from high-order-manager He/PQS testing on
2026-06-26 found that the installed `n_s = 5` mapped-COMX recipe is not robust
enough for all-electron scalar capture. It remains an opt-in construction
choice and must not be promoted to a default route from the first angular proxy
or driver-wiring success. The next mapped-COMX promotion evidence should be
bounded He `n_s = 6` and `n_s = 7` H1/IDA testing with shell-restricted scalar
capture diagnostics.

### HP-MCOMX-FILE-01 — mapped-COMX source-span files

Approved source files:

```text
src/cartesian_nested_faces.jl
src/cartesian_pair_block_materialization/pqs_source_axis_transforms.jl
src/cartesian_raw_product_sources/axis_transform_facts.jl
src/cartesian_raw_product_sources/records.jl
```

Primary owner: the existing nested doside / COMX source-span seam in
`src/cartesian_nested_faces.jl`.

`src/cartesian_pair_block_materialization/pqs_source_axis_transforms.jl` is
approved only for narrow keyword/spec plumbing into the existing doside seam
and for reporting the returned source-span facts.

`CartesianRawProductSources` files are approved only for compact provenance or
accessors on existing `AxisSourceTransformFact` records if needed. They must
not own the numerical mapped-COMX span builder.

No new source file is approved. In particular,
`src/cartesian_raw_product_sources/mapped_comx_source_span.jl` is not approved.

### HP-MCOMX-OBJ-01 — mapped-COMX source-span specification

Approved object: a compact `MappedCOMXSourceSpec` or equivalent typed
source-span specification.

Approved fields/semantics:

- `protected_degree`, initially `2`;
- `lambda`, initially `0.5`;
- `mapped_family`, initially `:chebyshev_s`;
- `include_sqrt_jacobian = false`;
- `localization_coordinate = :physical_u`;
- requested and resolved one-dimensional source mode count;
- deterministic mapped-order list or equivalent compact source-span metadata;
- compact rank and overlap/orthogonality diagnostics.

The object is source-span construction data at the doside seam, not a route
report, status payload, public input object, artifact schema,
`CartesianRawProductSources` numerical builder, or high-order benchmark record.

### HP-MCOMX-FN-01 — mapped-COMX source-span construction

Approved behavior:

- extend `_nested_doside_1d(...)` / `_nested_retained_span(...)` with an
  internal keyword or spec that changes only the raw source-span columns;
- keep the current ordinary span as the default behavior;
- build protected physical polynomial columns `1, u, u^2` using normalized
  local `u`;
- add mapped Chebyshev columns `T_k(s_lambda(u))` until the requested source
  mode count is reached;
- project mapped columns against the protected physical block in the local
  parent metric;
- orthonormalize the combined one-dimensional source span;
- continue through the existing `_cleanup_comx_transform(...)` using the
  physical position matrix;
- return materialized `AxisSourceTransformFact`s compatible with existing
  `RawProductBoxPlan` and PQS boundary product-mode retained rules;
- preserve deterministic ordering and ordinary polynomial source-span behavior.

The first implementation is restricted to `protected_degree = 2`. General
protected degrees require a later parity-balanced mapped-order fill rule; they
must not be implied by blindly adding `T_1`, `T_2`, and so on.

### HP-MCOMX-WIRE-01 — raw-source and PQS axis-transform wiring

Approved behavior:

- make the mapped-COMX option reachable through internal construction controls
  needed by the approved validation gates;
- carry descriptive source-span provenance such as
  `source_span_family = :mapped_comx`, protected degree, lambda, mapped family,
  mapped orders, `include_sqrt_jacobian = false`, and
  `localization_coordinate = :physical_u`;
- keep Hamiltonian, operator, and artifact layers consuming the usual
  carried-space / raw product source facts;
- avoid branching downstream Hamiltonian construction on ordinary polynomial
  versus mapped-COMX source spans except for descriptive provenance already
  carried by source facts.

Internal consumers must not read provenance metadata as a data bus. If mapped
orders or source-span diagnostics are needed after construction, return them as
real result fields or accessors; metadata remains reporting/provenance.

Forbidden for all `HP-MCOMX-*` IDs:

- changing default source spans;
- public API or export changes;
- canonical driver input changes;
- artifact schema, manifest, or reader changes;
- Hamiltonian, one-body, IDA, MWG, Residual Gaussian, raw Gaussian block, or
  solver changes;
- ECP, EGOI, RHF, ED, DMRG, or Cr2 workflow;
- explicit `Y_lm` / angular injection;
- `sqrtJ` weighting;
- mapped-`s` COMX as the production localization gauge;
- applying `s_lambda` to raw physical centers;
- `protected_degree != 2` in the first implementation;
- a parallel mapped-COMX axis-transform route outside the existing doside /
  COMX path;
- numerical source-span builders under `CartesianRawProductSources`;
- `src/cartesian_raw_product_sources/mapped_comx_source_span.jl`;
- importing high-order branch scaffolding, scripts, route wrappers, status
  objects, diagnostics, or reports;
- a duplicate high-order-maintained implementation of the same mainline
  option;
- committed Cr fixtures or broad high-order benchmark fixtures.

Failure rule: if the source option cannot be installed as a small branch in
the existing doside source-span seam before `_cleanup_comx_transform(...)`,
make no source commit and report the exact missing mainline seam. If the pass
needs a new source file, a second COMX wrapper, a
`CartesianRawProductSources` numerical builder, Hamiltonian assembly changes,
artifact schemas, public driver inputs, or high-order-specific workflow,
request a separate amendment.

### HP-MCOMX-TEST-01 — mapped-COMX validation

Approved validation:

- `git diff --check`;
- package load;
- local source-span validation for `n_s = 5`, `6`, and `7`:
  - full retained rank;
  - protected `P2` span preserved;
  - mapped columns use normalized local `u in [-1, 1]`;
  - source columns/centers match the high-order scratch convention within a
    reviewed tolerance;
  - source-axis overlap approximately identity after construction;
  - physical-`u` COMX off-diagonal residual reported;
  - metadata records the approved source-span rule;
- bounded cubic H and He+ one-electron gate comparing ordinary polynomial and
  mapped-COMX source spans with fixed support and retained count;
- bounded He `1s^2` fixed-orbital IDA gate if the already-supported analytic
  path can run without new solver or artifact workflow;
- high-order-manager consumer benchmarks on the installed mainline option for
  Cr occupied capture, reported back as evidence rather than committed mainline
  fixtures;
- no Cr2 run.

No committed test file is approved by default. A later implementation blurb may
name a small standalone script or ignored probe for the approved gates.
Committed fixtures, public driver tests, solver tests, and Cr/Cr2 benchmark
fixtures require a separate amendment.

### HP-MCOMX-TERM-FN-01 — mapped-COMX terminal-basis wiring

Approved source files:

```text
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
```

`CartesianFinalBasisRealization.jl` is approved only for import/include cleanup
if directly required.

Approved behavior:

- in `_shell_seed(...)`, prefer
  `contract.metadata.raw_product_source_axis_transform_facts` when present;
- validate exactly three axis facts;
- validate each fact is an `AxisSourceTransformFact`;
- validate `coefficient_status === :materialized`;
- validate source intervals match `support.outer_box`;
- validate source mode dimensions match `source_shape`;
- validate coefficient matrix row/column sizes match interval length and
  source mode dimension;
- build `full_coefficients` from the materialized axis coefficient matrices;
- keep existing boundary mode selection, support restriction,
  shell-local Lowdin, canonicalization, and support validation;
- preserve ordinary fallback through
  `_nested_projected_q_shell_full_sides(...)` when materialized facts are
  absent.

This ID makes carried mapped-COMX axis facts basis-defining at terminal shell
realization. It does not approve changing source-span construction, shell
ownership, retained-rule semantics, Lowdin realization, artifact schemas,
driver inputs, Hamiltonian assembly, IDA, MWG, RG, raw Gaussian blocks, solver,
EGOI, Cr2, high-order workflow, or defaults.

Failure rule: if terminal realization cannot consume carried axis facts without
changing shell ownership, retained-rule semantics, Lowdin realization, artifact
schema, or driver inputs, make no source commit and report the exact blocker.

### HP-MCOMX-TERM-TEST-01 — mapped-COMX terminal seam validation

Approved validation:

- `git diff --check`;
- package load;
- ordinary PQS H2 endpoint/regression unchanged;
- mapped source-span probe still passes;
- focused He or H terminal seam check showing mapped terminal shell
  coefficients differ from ordinary and match the carried materialized axis
  facts;
- H2 supplemented RG endpoint if the touched path crosses it;
- no Cr2 run.

No driver input test, artifact schema test, committed test file, committed
fixture, high-order benchmark fixture, or Cr2 fixture is approved.

### HP-MCOMX-DRV-FN-01 — mapped-COMX canonical driver selection

Approved source files:

```text
bin/cartesian_ham_builder.jl
src/cartesian_base_hamiltonian.jl
src/pqs_source_box_route_driver_helpers.jl
```

`src/pqs_source_box_route_driver_helpers.jl` is approved only for narrow
propagation of a normalized source-span selector to the already-approved
source-axis transform fact path. This ID must not add route records,
route-stage diagnostics, terminal-lowering contracts, artifact fields, or a
new COMX implementation.

Approved behavior:

- add a visible canonical-driver input `source_span`;
- support trusted input-file and command-line override handling through the
  existing compact driver input mechanism;
- print `source_span` in the compact public contract when `print_contract` is
  enabled;
- accept only `:ordinary` and `:mapped_comx` as driver-facing values;
- default to `:ordinary`;
- normalize and validate the selector in `src/cartesian_base_hamiltonian.jl`;
- pass `:mapped_comx` through the existing PQS source-box path so terminal
  realization receives materialized mapped-COMX axis facts;
- reject `source_span = :mapped_comx` clearly for `nesting = :wl` until a
  future WL-specific amendment approves otherwise;
- keep ordinary driver artifact/readback behavior unchanged.

This is a public construction choice, not a diagnostic route switch.

Forbidden:

- new driver hooks beyond `source_span`;
- route skeleton, route record, terminal-lowering, retained-rule, shell
  ownership, artifact schema, manifest, reader, Hamiltonian, IDA, MWG, RG, raw
  Gaussian block, solver, EGOI, Cr2, high-order workflow, or source-default
  changes;
- another COMX path, route-stage diagnostics, stop-after controls, raw-block
  switches, allocation probes, committed tests, or committed fixtures.

Failure rule: if making `source_span` driver-selectable requires new route
records, terminal-lowering changes, artifact schema changes, or another COMX
path, make no source commit and report the exact blocker.

### HP-MCOMX-DRV-TEST-01 — mapped-COMX driver validation

Approved validation:

- `git diff --check`;
- package load;
- default ordinary driver artifact/readback still passes;
- mapped-COMX He or H PQS driver smoke proves carried facts are
  basis-defining;
- ordinary versus mapped He supplemented/MWG/IDA comparison through the real
  driver if bounded;
- H2 RG endpoint still passes;
- no Cr2 run.

No committed test file, committed fixture, high-order benchmark fixture,
artifact schema test, route-diagnostic test, or Cr2 fixture is approved.

## Approved Composition Lane: Base Homonuclear Z-Axis Diatomics

This section promotes the base z-axis diatomic validation relaxation. It
approves only the `geometry = z-axis diatomic`, `supplement = off` public input
contract in the base facade. It does not approve route/shellification changes.

### HP-COMP-BASEDIAT-FN-01 — base homonuclear z-axis diatomic validation

Approved source file:

```text
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- relax the `_cartesian_base_inputs(...)` two-center branch from H2-only to
  explicit homonuclear z-axis all-electron diatomics;
- require exactly two atom symbols and exactly two nuclear charges;
- require equal atom symbols;
- require equal finite positive nuclear charges that are integer-valued after
  validation;
- require two finite centers on the Cartesian z axis (`x = y = 0`) with
  distinct `z`;
- require explicit `nup` and `ndn` nonnegative integers with
  `nup + ndn == sum(nuclear_charges)` after charge validation;
- keep atom symbol as provenance label only and nuclear charge as authority;
- keep the basis contract unchanged except for `HP-COMP-NS-*` public size
  naming: `ns`, `core_spacing`, `xmax_parallel`,
  `xmax_transverse`, optional `parent_axis_family`, `reference_spacing`,
  `tail_spacing`, and `nesting`;
- preserve both `nesting = :pqs` and `nesting = :wl` as public construction
  family choices;
- preserve current H2 behavior and existing route provenance labels
  (`:z_axis_diatomic_pqs_base`, and `:z_axis_diatomic_wl_base` only when the
  WL diatomic terminal-record lane succeeds).

The `nesting = :wl` path still depends on `HP-COMP-WLDIAT-FN-01` for native WL
diatomic terminal records. This ID may allow validated non-H diatomic inputs to
reach the same WL route path, but it must not implement terminal records,
shellification, or route lowering outside `src/cartesian_base_hamiltonian.jl`.

Forbidden:

- driver changes;
- source files outside `src/cartesian_base_hamiltonian.jl`;
- supplement, Residual Gaussian, MWG/IDA, or artifact-manifest changes;
- route skeleton, shellification, terminal lowering, raw-block, writer,
  reader, public API/export, solver/ECP, diagnostics, status/report payload,
  or Cr2-specific workflow changes;
- element lookup/default tables or element-inferred electron counts;
- heteronuclear, translated, non-z-axis, or general-geometry support;
- committed tests or committed input fixtures.

Failure rule: if non-H base diatomic construction requires route/shellification
changes outside `src/cartesian_base_hamiltonian.jl`, make no source commit and
report the exact blocker.

Line budget: target under `60` added `src` lines.

### HP-COMP-BASEDIAT-TEST-01 — base homonuclear z-axis diatomic validation

Approved validation:

- `git diff --check`;
- package load;
- H2 base artifact/readback unchanged for `nesting = :pqs`;
- H2 base artifact/readback for `nesting = :wl` if the WL diatomic terminal
  lane is implemented in the working tree; otherwise report the existing
  `HP-COMP-WLDIAT-*` dependency;
- one small explicit non-H homonuclear z-axis base artifact/readback,
  preferably Be2 or N2 with bounded `q`, for `nesting = :pqs`;
- the same small non-H base artifact/readback for `nesting = :wl` if runtime is
  bounded and the WL diatomic terminal lane is implemented in the working tree;
- clear rejection for heteronuclear symbols, unequal charges, non-neutral
  electron count, and non-z-axis centers;
- no Cr2 run.

No committed test file, committed fixture, driver contract test, supplemented
run, solver/RHF/ECP/EGOI validation, route-diagnostic validation, or Cr2
fixture is approved.

## Approved Composition Lane: Supplemented White-Lindsey Z-Axis Diatomics

This section promotes the supplemented WL z-axis diatomic composition lane.
The goal is to make `Natom = 2`, `basisname !== nothing`, and
`nesting = :wl` use the same supplemented homonuclear z-axis diatomic staged
facade as `nesting = :pqs`, after the WL base route has produced a valid
`CartesianTerminalBasisRealization`.

### HP-COMP-SUPPWL-FN-01 — supplemented WL z-axis diatomic composition

Approved source files:

```text
src/cartesian_base_hamiltonian.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

`src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` is
optional in this lane and may be edited only if a direct genericity blocker
appears in the existing RG/MWG compatibility entry point. The default expected
source change is in `src/cartesian_base_hamiltonian.jl`.

Approved behavior:

- allow `nesting = :wl` through the existing supplemented homonuclear z-axis
  diatomic facade/staged path;
- remove the early supplemented-`nesting = :wl` blockers in
  `cartesian_base_working_basis(...; supplemented = true)` and
  `cartesian_residual_gto_supplement_basis(...)` only if the existing
  Residual Gaussian/MWG path works with the WL `CartesianTerminalBasisRealization`;
- preserve the existing supplement contract: `basis_by_center`, `lmax`,
  optional `uncontracted`, `width_filtering`, and `basisfile`;
- preserve residual selection, exact augmented operators, residual MWG/IDA
  interaction, base K/U reuse, artifact keys, manifest/provenance, driver
  inputs, and stage labels;
- keep `nesting` as a construction-family choice, not a diagnostic route
  switch.

Forbidden:

- driver changes;
- supplemented atoms;
- route skeleton, shellification, terminal lowering, raw-block,
  residual-selection, MWG/IDA convention, artifact schema, writer, reader,
  public API/export, solver/ECP, diagnostics, status/report payload, or
  Cr2-specific workflow changes;
- old White-Lindsey H1/H1+J materialization revival or adaptation;
- committed tests or committed input fixtures.

Failure rule: if supplemented WL requires new terminal records, route lowering
semantics, residual-selection changes, artifact/schema changes, or source files
outside the approved surface, make no source commit and report the exact
blocker.

Line budget: target under `80` added `src` lines, with deletion or
simplification of the early supplemented-WL blockers expected where practical.

### HP-COMP-SUPPWL-TEST-01 — supplemented WL validation

Approved validation:

- `git diff --check`;
- package load;
- H2 supplemented artifact/readback with `nesting = :pqs`;
- small H2 supplemented artifact/readback with `nesting = :wl`;
- small non-H homonuclear z-axis supplemented artifact/readback with
  `nesting = :wl` if bounded;
- finite/symmetric `K` and `V` checks plus readback deltas for WL supplemented
  output;
- clear rejection remains for supplemented atoms and invalid diatomic inputs;
- no Cr2 run.

No committed test file, committed fixture, driver contract test,
solver/RHF/ECP/EGOI validation, route-diagnostic validation, or Cr2 fixture is
approved.

## Approved Composition Lane: Supplemented One-Center Atoms

This section promotes the supplemented one-center atom composition lane. The
goal is to make `Natom = 1`, `basisname !== nothing`, and either
`nesting = :pqs` or `nesting = :wl` use the same staged supplemented
Hamiltonian machinery as supported diatomics, with one physical owner center.

### HP-COMP-SUPPATOM-FN-01 — supplemented atom composition

Approved source files:

```text
src/cartesian_base_hamiltonian.jl
bin/cartesian_ham_builder.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

`src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` is
optional in this lane and may be edited only if a direct one-owner RG/MWG
genericity blocker appears. The default expected source changes are in
`src/cartesian_base_hamiltonian.jl` and `bin/cartesian_ham_builder.jl`.

Approved behavior:

- allow origin-centered one-center all-electron atoms with supplement enabled;
- preserve both `nesting = :pqs` and `nesting = :wl`;
- use the existing base atom validation and terminal basis construction;
- use the existing residual Gaussian augmentation, exact augmented operators,
  residual MWG/IDA interaction, base K/U reuse, assembly, writer, readback,
  manifest, and provenance;
- use `legacy_atomic_gaussian_supplement(...)` when `length(input.symbols) == 1`;
- keep `legacy_bond_aligned_diatomic_gaussian_supplement(...)` when
  `length(input.symbols) == 2`;
- relax only the canonical driver's supplemented `Natom == 2` guard so
  `basisname !== nothing` works for both `Natom = 1` and `Natom = 2`;
- keep driver public inputs, ordering, comments, hooks, spacing/layout, stage
  labels, and artifact contract otherwise unchanged.

Forbidden:

- separate atom-only Hamiltonian builder or materialization path;
- new driver inputs, route switches, diagnostics, stop-after controls, or
  stage labels;
- route skeleton, shellification, terminal lowering, raw-block,
  residual-selection, MWG/IDA convention, artifact schema, writer, reader,
  public API/export, solver/ECP, status/report payload, heteronuclear/general
  geometry, translated atom, or Cr2-specific workflow changes;
- committed tests or committed input fixtures.

Failure rule: if supplemented atoms require new residual-selection semantics,
atom-only MWG/IDA conventions, terminal-lowering changes, route changes,
artifact schema changes, or a separate atom Hamiltonian path, make no source
commit and report the exact blocker.

Line budget: target under `80` added `src`/`bin` lines, with deletion or
simplification expected where guards become obsolete.

### HP-COMP-SUPPATOM-TEST-01 — supplemented atom validation

Approved validation:

- `git diff --check`;
- package load;
- base atom artifact/readback still passes;
- H atom supplemented artifact/readback with `nesting = :pqs`;
- H atom supplemented artifact/readback with `nesting = :wl`;
- small Be atom supplemented artifact/readback for `nesting = :pqs` or
  `nesting = :wl` if bounded;
- H2 supplemented PQS and WL smoke still pass;
- translated atom still rejects clearly;
- supplement basis count mismatch still rejects clearly;
- no Cr2 run.

No committed test file, committed fixture, driver contract test,
solver/RHF/ECP/EGOI validation, route-diagnostic validation, or Cr2 fixture is
approved.

## Approved Contract Lane: One-Center Atom Parent Sizing

This section promotes a narrow atom box-sizing correction. The current
one-center atom producer still contains an old minimal-test artifact in which
`parent_axis_counts = (2*q + 1, 2*q + 1, 2*q + 1)`. That makes public
`padding` / `basis.radius` non-authoritative for atom box size. The approved
contract is that one-center atom sizing uses the same physical-size idea as
z-axis diatomics: the public physical extent and spacing policy determine the
parent box/counts. Public size-parameter naming is separately normalized under
`HP-COMP-NS-*`.

### HP-COMP-ATOMBOX-FN-01 — one-center atom parent sizing

Approved source file:

```text
src/cartesian_base_hamiltonian.jl
```

No optional helper file is approved in this amendment. If a later source audit
finds an existing parent/system sizing helper that should own the shared rule,
that exact file must be named in a separate docs-only amendment.

Approved behavior:

- remove the fixed `2*q + 1` atom parent-axis count artifact;
- use `input.radius` / public `basis.radius` as the one-center physical box
  extent authority;
- make atom parent axis counts depend on radius plus `core_spacing` / existing
  spacing policy, analogous to the z-axis diatomic physical-extent sizing;
- preserve the public size/nesting control as resolution metadata, not direct
  box side count;
- preserve origin-centered atom validation, explicit charge/electron-count
  validation, `nesting = :pqs` and `nesting = :wl`, supplemented atoms,
  artifact keys, manifest/provenance, and canonical driver inputs;
- ensure different atom padding/radius values can produce different atom
  dimensions when the physical extent changes enough.

Forbidden:

- driver changes;
- route-family switches, raw-block changes, residual-selection changes,
  MWG/IDA convention changes, artifact schema changes, writer/reader changes,
  public API/export changes, solver/ECP work, diagnostics/status/report
  payloads, committed tests, or Cr2-specific workflow changes;
- translated atoms, non-origin atom support, or element lookup/default tables;
- broad parent-construction rewrites or diatomic sizing changes.

Failure rule: if fixing atom parent sizing requires broad parent-construction
redesign, route semantics changes, driver contract changes, artifact schema
changes, translated-atom support, or source files outside the approved surface,
make no source commit and report the exact blocker.

Line budget: target under `80` added `src` lines.

### HP-COMP-ATOMBOX-TEST-01 — atom parent-sizing validation

Approved validation:

- `git diff --check`;
- package load;
- H atom base artifact/readback still passes;
- H atom supplemented PQS and WL artifact/readback still pass;
- small Be atom base and supplemented smoke if bounded;
- atom padding sensitivity check: same atom inputs with two padding/radius
  values should show changed parent counts/dimension or a clear explanation if
  both values fall in the same count bin;
- H2/Be2 diatomic smoke to confirm no diatomic sizing regression;
- no Cr2 run.

No committed test file, committed fixture, driver contract test,
solver/RHF/ECP/EGOI validation, route-diagnostic validation, or Cr2 fixture is
approved.

## Approved Contract Lane: Public `ns` Input And Derived `q`

This section approves a narrow public input cleanup. The prior public contract
used `q` for the requested cube/source/nesting size, but White-Lindsey uses a
route-local `q` that is naturally two smaller than the requested `ns` source
size. That makes `q` an ambiguous public knob across the composed
`nesting = :pqs | :wl` workflow.

The durable public field is now `ns`: the requested cube/source/nesting size.
Route-local `q` is derived after the construction family is selected.

### HP-COMP-NS-FN-01 — public `ns` normalization

Approved source files:

```text
src/cartesian_base_hamiltonian.jl
bin/cartesian_ham_builder.jl
```

Approved behavior:

- prefer `basis.ns` / driver `ns` in public examples, templates, and driver
  contract construction;
- derive route-local `q` from the normalized public `ns` and `nesting`:
  - `nesting = :pqs`: `q = ns`;
  - `nesting = :wl`: `q = ns - 2`;
- reject `nesting = :wl` with `ns < 3`;
- keep `ns` separate from physical box controls: atom `radius` / driver
  `padding`, diatomic extents, `core_spacing`, `reference_spacing`, and
  `tail_spacing`;
- temporarily accept legacy public `q` only as compatibility:
  - if `ns` is absent, derive `ns = q` for `nesting = :pqs`;
  - if `ns` is absent, derive `ns = q + 2` for `nesting = :wl`;
  - if both `ns` and `q` are present, require consistency with the selected
    nesting or throw `ArgumentError`;
- record compact provenance in existing groups:
  - `ns`: normalized public requested source size;
  - `q`: derived route-local value passed to route construction;
  - `q_rule`: `:pqs_ns_equals_q` or `:wl_ns_minus_2`;
  - `ns_source`: `:public_ns` or `:legacy_q_compatibility`.

Forbidden:

- route skeleton, shellification, terminal lowering, raw-block,
  residual-selection, MWG/IDA, numerical-kernel, solver/ECP, or Cr2-specific
  workflow changes;
- public API/export redesign beyond accepting/prefering `ns` in the existing
  public groups;
- artifact matrix-key changes, reader behavior changes, new artifact format,
  route diagnostics, status/report payloads, raw-block switches, or new driver
  hooks/stage labels;
- committed tests or fixtures.

Failure rule: if this cleanup requires route-stage redesign, terminal
construction changes, artifact reader changes, or a broad compatibility layer
outside the approved files, make no source commit and report the blocker.

Line budget: target under `80` added source/bin lines, with net simplification
where legacy `q` plumbing becomes local normalization.

### HP-COMP-NS-TEST-01 — public `ns` validation

Approved validation:

- `git diff --check`;
- package load;
- atom and z-axis diatomic base artifact/readback using public `ns` under
  `nesting = :pqs`;
- small atom and z-axis diatomic base artifact/readback using public `ns` under
  `nesting = :wl` where the corresponding WL construction cell is already
  supported;
- supplemented atom/diatomic smoke using public `ns` only where the
  corresponding supplemented composition cell is already supported and bounded;
- legacy `q` compatibility smoke if the implementation keeps it;
- clear rejection for inconsistent `ns`/`q`, `nesting = :wl` with `ns < 3`,
  and unsupported geometry/supplement combinations;
- direct provenance inspection for `ns`, derived `q`, `q_rule`, `ns_source`,
  and `nesting`;
- no Cr2 run.

No committed test file, committed fixture, driver contract test,
solver/RHF/ECP/EGOI validation, route-diagnostic validation, or Cr2 fixture is
approved.

## Approved Contract Lane: WL Diatomic `ns` Early Rejection

This section approves a narrow cleanup of the public `ns` contract for
White-Lindsey z-axis diatomics. Current `nesting = :wl`, `Natom = 2`,
`ns = 3` normalizes to route-local `q = 1` and `core_cube_side = 1`, then
fails later in terminal shellification because a complete-shell inner box
cannot be formed. That input should be rejected early as unsupported WL
diatomic input, not discovered after route construction.

The same audit showed that working WL diatomic `ns` ranges may change
route-local `q`, shellification block decomposition, and final row order
without changing the retained support set over a fixed physical parent extent.
That is retained-support saturation, not ignored public input.

### HP-COMP-WLNS-FN-01 — WL diatomic `ns` contract cleanup

Approved source file:

```text
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- for homonuclear z-axis diatomics with `nesting = :wl`, reject normalized
  `ns < 4` early with a clear `ArgumentError`;
- preserve current public size normalization:
  - `nesting = :pqs`: route-local `q = ns`;
  - `nesting = :wl`: route-local `q = ns - 2`;
- preserve all working cells in the current composition matrix: atom and
  z-axis diatomic geometry, base and supplemented mode, and `nesting = :pqs`
  / `:wl`;
- document and preserve the fact that WL diatomic retained support may
  saturate over `ns` ranges when the physical parent extent dominates; the
  same public `ns` value is not a fair PQS/WL basis-size comparison.

Forbidden:

- driver public-input changes, layout/comment/stage/timing changes, or
  compatibility aliases;
- route skeleton, shellification, terminal lowering, retained-unit,
  terminal-basis, raw-block, Residual Gaussian, MWG/IDA, artifact schema,
  writer, reader, public API/export, solver/ECP, diagnostics, report/status
  payload, or Cr2-specific workflow changes;
- attempts to "fix" WL `ns = 4:7` retained-support saturation in this lane;
- committed tests or committed input fixtures.

Failure rule: if early WL diatomic `ns` rejection cannot be implemented in
`src/cartesian_base_hamiltonian.jl` without touching driver, shellification,
terminal lowering, route semantics, or artifact schema, make no source commit
and report the exact blocker.

Line budget: target under `60` added `src` lines.

### HP-COMP-WLNS-TEST-01 — WL diatomic `ns` validation

Approved validation:

- `git diff --check`;
- package load;
- WL Be2 or H2 z-axis diatomic with `ns = 3` rejects early;
- WL z-axis diatomic `ns = 4` base artifact/readback still passes;
- WL z-axis diatomic supplemented smoke still passes where that composition
  cell is already implemented and bounded;
- PQS atom and z-axis diatomic smoke still pass;
- no Cr2 run.

No committed test file, committed fixture, driver contract test,
solver/RHF/ECP/EGOI validation, route-diagnostic validation, or Cr2 fixture is
approved.

### HP-FN-03 — blockwise one-body assembly

Approved file:

```text
src/cartesian_final_basis_realization/pqs_terminal_one_body.jl
```

Approved public internal helper:

```julia
assemble_terminal_product_operator!(
    destination,
    basis::CartesianTerminalBasisRealization,
    axis_x,
    axis_y,
    axis_z;
    scale = 1.0,
)
```

Approved file-local Gaussian-sum helper:

```julia
_accumulate_terminal_gaussian_sum!(
    destination,
    basis::CartesianTerminalBasisRealization,
    coefficients,
    factors_x,
    factors_y,
    factors_z;
    scale = -1.0,
)
```

No `K`/`U_A` payload, stage field, report object, persistent cache, or
orchestration API is approved by HP-FN-03.

### HP-FN-04 — localized IDA assembly

Approved file:

```text
src/cartesian_final_basis_realization/pqs_terminal_ida.jl
```

Approved function:

```julia
assemble_terminal_ida_interaction!(
    destination,
    basis::CartesianTerminalBasisRealization,
    coefficients,
    raw_pair_terms_x,
    raw_pair_terms_y,
    raw_pair_terms_z,
    weights_x,
    weights_y,
    weights_z;
    weight_atol = 1.0e-12,
    symmetry_atol = 1.0e-10,
)
```

This is Slice C1 only: it produces final-basis `electron_electron_ida`. It does
not authorize Hamiltonian construction, route wiring, artifacts, or a pair
payload/cache.

### HP-FN-05 — final Hamiltonian construction

Approved as the narrow Slice C2 construction boundary for the existing
`CartesianIDAHamiltonian`.

Conceptual boundary:

```julia
build_cartesian_ida_hamiltonian(
    kinetic,
    nuclear_attraction_unit_by_center,
    electron_electron_ida,
    nup,
    ndn;
    nuclear_charges,
    nuclear_positions,
)::CartesianIDAHamiltonian{Float64}
```

Implementation may call the existing `CartesianIDAHamiltonian(...)`
constructor directly if no helper is needed.

### HP-WIRE-02 — historical direct materialization Hamiltonian handoff

Historically approved and implemented Slice D wrapper boundary:

```julia
cartesian_materialization(
    report,
    terminal_basis_realization,
    materialization_inputs,
)::Union{Nothing,CartesianIDAHamiltonian{Float64}}
```

This old route-driver wrapper workflow is now approved for retirement under
`HP-RETIRE-DRV-MAT-*`. Current canonical producer work should use the staged
driver-facing producer functions and `CartesianIDAHamiltonian` artifact path,
not add new callers to `cartesian_materialization`.

The call site passes `transforms.terminal_basis_realization` directly. The
terminal basis must not be embedded in `cartesian_report`, reconstructed from
summaries, or recovered by passing the full `transforms` stage.

Return contract:

- no request returns `nothing`;
- requested base PQS materialization returns `CartesianIDAHamiltonian{Float64}`;
- `save_ham_artifact = true` writes with `write_cartesian_ida_hamiltonian` and
  still returns the same Hamiltonian;
- no materialization wrapper, `result_kind`, `materialized`, status mirror, or
  `ida_hamiltonian` field is approved.

## Approved For R1 Implementation

These entries authorize only the R1 public base producer scope recorded in
`r1_public_base_producer.md`. They do not approve broad driver polish,
additional routes, new artifact shapes beyond the approved `HP-R1-ART-01`
`producer_provenance/` keys in the final Hamiltonian file, solver work,
supplements, corrections, or status/report/payload expansion.

### HP-R1-FILE-01 — public base producer source file

Approved file:

```text
src/cartesian_base_hamiltonian.jl
```

Approved owner: top-level `GaussletBases` public API.

The only new export approved by R1 is `cartesian_base_hamiltonian` from
`src/GaussletBases.jl`.

### HP-R1-FN-01 — public base Hamiltonian producer facade

Approved public call shape:

```julia
cartesian_base_hamiltonian(
    system::NamedTuple;
    basis::NamedTuple,
    hamfile::Union{Nothing,AbstractString} = nothing,
)::CartesianIDAHamiltonian{Float64}
```

Scope:

- origin-centered base H and Cartesian z-axis aligned base H2 first;
- plain `NamedTuple` input groups only;
- no public `method`, `route`, or output group in R1;
- `n_s`, bond length, and private H2 radius are derived internally;
- `core_spacing` is the authoritative physical near-nucleus spacing for each
  center;
- public `d` is deprecated and is not part of the durable producer contract;
- one-center H maps `reference_spacing` and `core_spacing` separately:
  `spacing_inputs.reference_spacing = basis.reference_spacing`,
  `parent_inputs.parent_mapping_rule = :white_lindsey_atomic_mapping`, and
  `parent_inputs.parent_mapping_d = basis.core_spacing`;
- the reviewed H baseline uses explicit public `core_spacing = 0.3` and
  `reference_spacing = 1.0`;
- if a temporary compatibility path accepts public `d`, it must require
  `d == resolved core_spacing`; z-axis H2 and durable public examples must not
  use `d`;
- x/y-aligned diatomics, shifted-parallel diatomics, generally oriented
  molecules, translation, and rotation are deferred;
- center-sized public collections must be vectors or other `AbstractVector`
  values, not variable-size tuples;
- unknown public input keys throw `ArgumentError`;
- scalar inputs must be positive and finite where applicable;
- symbols, charges, coordinates, and electron counts must match the approved
  H/H2 scope;
- return the existing `CartesianIDAHamiltonian{Float64}` directly;
- no wrapper, payload, status object, report mirror, or new artifact shape
  except the approved `HP-R1-ART-01` `producer_provenance/` keys in the final
  Hamiltonian file;
- non-`nothing` `hamfile` writes with existing
  `write_cartesian_ida_hamiltonian`; production does not automatically
  read back the artifact.

### HP-R1-CORE-FN-01 — unified core-spacing producer contract

Approved source file:

```text
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- `core_spacing` is the single public physical near-nucleus spacing for each
  center after explicit input or preset resolution;
- public `d` is deprecated as a producer field;
- a temporary compatibility acceptance of `d` may exist only when
  `d == resolved core_spacing`; mismatches must throw `ArgumentError`;
- one-center White-Lindsey atom wiring sets
  `parent_inputs.parent_mapping_d = resolved core_spacing`;
- the White-Lindsey `Z` dependence is an internal mapping-shape rule:
  `core_range = sqrt(core_spacing / Z)` and
  `mapping_strength = sqrt(core_spacing * Z)`;
- multi-center mappings use the same per-center resolved-core-spacing model
  before applying the combined/neighbor mapping effects;
- future automatic presets may derive `core_spacing = core_scale(q or n_s) / Z`
  or an equivalent fixed `core_spacing * Z` family, but once resolved,
  `core_spacing` remains the authoritative scale.
- canonical driver or project-input defaults may choose visible editable values
  such as `core_spacing = 0.3`; these are explicit resolved inputs, not hidden
  universal producer defaults, and normal overrides such as quick-test
  `core_spacing = 0.5` remain allowed.
- routine correctness tests may override driver physics defaults, but any
  asserted scalar must be tied to the exact test input and not described as a
  physics-default result.

This ID does not approve a public `d`, public `parent_mapping_d`, public
mapping-strength/range knobs, element-table defaults, ECP behavior, solver
workflow, artifact schema changes, or source files outside
`src/cartesian_base_hamiltonian.jl`. `reference_spacing`, `tail_spacing`, and
physical box padding remain separate concepts.

### HP-R1-WIRE-01 — report-free base producer wiring

Approved wiring for the R1 public facade:

```text
system / specification
-> parent and route geometry
-> terminal basis realization
-> Hamiltonian production
-> optional existing Hamiltonian artifact
```

The recommended base-public path must not require `cartesian_pair_terms` or
`cartesian_assembly`. Existing stages may remain temporarily for legacy
script/report compatibility, but this R1 authority does not approve adding a
new base-route consumer to either stage.

Approved private shared constructor seam:

```julia
_cartesian_base_ida_hamiltonian(
    terminal_basis_realization,
    parent_axis_bundle_object,
    atom_locations::Vector{NTuple{3,Float64}},
    nuclear_charges::Vector{Float64},
    nup::Int,
    ndn::Int,
)::CartesianIDAHamiltonian{Float64}
```

Approved owner file:

```text
src/pqs_source_box_low_order_materialization.jl
```

This PQS-named owner is acceptable for the explicitly PQS-only R1 migration.
It is not permanent method-neutral ownership for R2.

Allowed caller files:

- `src/pqs_source_box_low_order_materialization.jl`
- `src/cartesian_base_hamiltonian.jl`

The existing report-bound materialization path and the new public facade should
call this same private constructor. The implementation must remove or bypass
the current report dependency where `cartesian_assembly` exists chiefly to
supply `route_skeleton` and a low-order shellification summary to
`cartesian_report`, and delete the duplicated report-bound numerical
orchestration it replaces. It must not fabricate a private report object,
duplicate the Hamiltonian builder in the new public file, leave two parallel
base Hamiltonian construction paths, or replace the dependency with a new
report field cloud, status payload, or metadata-carried numerical data.

### HP-ROUTE-RECIPE-FN-01 — family-selective route recipe cleanup

Approved source files:

```text
src/pqs_source_box_route_driver_helpers.jl
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- `cartesian_recipe(route_inputs)` may construct only the subrecipe selected by
  `route_inputs.route_family`;
- for `route_family = :pqs_source_box`, route inputs must not require inactive
  `white_lindsey_*` fields; the produced recipe may set the inactive
  `white_lindsey` subrecipe to `nothing` while retaining the existing field
  name for caller compatibility;
- for `route_family = :white_lindsey_low_order`, explicit White-Lindsey route
  support must be preserved and the selected `white_lindsey` subrecipe must
  continue to be built from the existing WL route fields; the inactive
  `source_box` subrecipe may be `nothing` if no live WL caller requires it;
- `_cartesian_base_route(kind)` in `src/cartesian_base_hamiltonian.jl` may
  remove unused `white_lindsey_*` fields because the live base producer route
  uses `route_family = :pqs_source_box`;
- existing precomposed recipes that already provide `source_box` and
  `white_lindsey` fields may remain accepted if that compatibility path is
  still live, but it must not force new PQS-only route inputs to carry inactive
  WL vocabulary.

This ID preserves real WL/PQS algorithm differences while removing inactive WL
route-family fields from the current PQS base producer contract. It does not
approve canonical-driver changes, numerical kernel changes, terminal lowering
policy changes, shellification behavior changes, materialization or artifact
schema changes, route-stage diagnostics, status/report expansion, deletion of
WL materialization, or source files outside the two approved files.

Line budget: at most `80` added `src` lines, with net simplification expected.

Failure rule: if `cartesian_recipe(...)` cannot be made family-selective
without broader route-driver, report, materialization, or stage-object changes,
make no source commit and report the blocker.

### HP-ROUTE-RECIPE-TEST-01 — route recipe cleanup validation

Approved validation:

- `git diff --check`;
- package load;
- H atom/base artifact readback;
- H2 base artifact readback;
- compact H2 supplemented facade or driver path;
- focused route recipe smoke for explicit `:white_lindsey_low_order` if still
  practical, or a report of the exact live test/tool callers that block further
  WL route-input cleanup.

Existing committed tests may be adjusted only where they directly construct
route inputs that now no longer need inactive family fields. Known direct
route-recipe tests such as
`test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` may drop inactive
WL fields if required by the source cleanup. No new committed test file, Cr2
run, driver workflow change, route diagnostic, or physics-reference scalar is
approved by this ID.

### HP-ROUTE-INV-FN-01 — retained-unit route inventory type cleanup

Approved source file:

```text
src/pqs_source_box_route_driver_helpers.jl
```

Approved cleanup targets:

- `_pqs_source_box_route_driver_named_tuple_from_units(...)`;
- runtime-keyed retained-unit inventory fields derived from unit labels,
  including `source_boxes`, `source_dimensions`, `retained_counts`, and
  `ranges`;
- runtime-keyed `pair_family_counts = NamedTuple{families}(...)`;
- same-file internal consumers that currently expect those runtime-keyed
  `NamedTuple` shapes.

Approved replacements:

- vector-backed records or tables with stable field names;
- stable dictionaries keyed by unit or pair-family labels only where lookup by
  label is genuinely needed;
- helper accessors that hide the storage shape from same-file callers;
- compact summaries that expose counts/order without encoding route size in the
  concrete type.

The retained-unit vector remains the ordered inventory authority. Unit labels
and pair-family labels may remain data values, but they must not become type
parameters.

This ID does not approve edits to `RawProductBoxPlan.source_mode_indices`,
`source_mode_column_indices`, `source_mode_indices(...)`,
`TerminalLoweringPlan.available_contracts`, `TerminalLoweringPlan.contracts`,
or `RetainedUnitTransformContractPlan.contracts`. It also does not approve
public input `NamedTuple` changes, fixed `NTuple{3,Int}` coordinate/dimension
changes, artifact sidecar table changes, numerical kernels, route recipe
behavior changes, shellification, terminal lowering, terminal basis, Residual
Gaussian, raw-block changes, canonical driver changes, Hamiltonian object
changes, matrix-key changes, reader changes, public API/export changes,
report/status/payload expansion, compatibility adapters, new committed tests,
Cr2 runs, or Cr2-specific workflow.

Line budget: at most `120` added `src` lines, with net simplification expected.
Failure rule: if the cleanup requires source files outside the approved file,
broader route/stage rewiring, public API changes, artifact changes, or an
adapter that preserves the old runtime-keyed type shape, make no source commit
and report the blocker.

### HP-ROUTE-INV-TEST-01 — retained-unit route inventory cleanup validation

Approved validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader or
  canonical driver path;
- focused search confirming no `NamedTuple{unit_keys}` or
  `NamedTuple{families}` route inventory remains in
  `src/pqs_source_box_route_driver_helpers.jl`;
- no Cr2 run.

Existing committed tests may be adjusted only if they directly assert the old
runtime-keyed inventory shape. No new committed test file, Cr2 fixture,
driver-input fixture, benchmark, or route-diagnostic test is approved by this
ID.

### HP-RAW-SRCMODE-FN-01 — raw product source-mode inventory cleanup

Approved source files:

```text
src/cartesian_raw_product_sources/records.jl
src/cartesian_raw_product_sources/source_mode_indices.jl
src/cartesian_raw_product_sources/summaries.jl
```

Approved narrow consumer files, only as required by the storage change:

```text
src/cartesian_retained_unit_transform_contracts/unit_contracts.jl
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
src/cartesian_base_hamiltonian.jl
```

Approved cleanup targets:

- `RawProductBoxPlan.source_mode_indices::Tuple{Vararg{NTuple{3,Int}}}`;
- `RawProductBoxPlan.source_mode_column_indices::Tuple{Vararg{Int}}`;
- `source_mode_indices(...)` / source-mode summary accessors only to the extent
  required to hide vector-backed storage from approved callers;
- same-file and listed narrow consumers that currently depend on the
  tuple-backed storage shape.

Approved replacement:

- vector-backed source-mode coordinate storage;
- vector-backed source-mode column storage, or no stored column vector when the
  column sequence is exactly `1:count` and accessors provide the same ordered
  column numbers;
- stable accessors preserving deterministic source-mode order, mode values,
  length, indexing/iteration where currently used, and retained-rule parity.

The fixed `NTuple{3,Int}` coordinate and dimension values remain valid. The
variable-length source-mode inventory must not be encoded in `RawProductBoxPlan`
field types or accessor return types. Accessor compatibility means same facts
and order, not the old concrete `Tuple{Vararg{...}}` return shape.

This ID does not approve terminal-lowering `contracts` /
`available_contracts` tuple cleanup, retained-unit transform-contract tuple
cleanup outside the listed narrow consumer wiring, broad pair-block/source-box
rewrites, numerical kernel changes, route semantic changes, public API/export
changes, canonical driver changes, Hamiltonian object changes, matrix-key
changes, reader changes, artifact schema changes, route-stage objects,
report/status/payload expansion, persistent caches, compatibility layers that
preserve the old tuple-backed shape, new committed tests, Cr2 runs, or
Cr2-specific workflow.

Line budget: at most `150` added `src` lines, with net simplification expected.
Failure rule: if vectorizing the raw product plan requires source files outside
the approved surfaces, broad pair-block/source-box rewrites, public API or
artifact changes, numerical changes, or compatibility layers preserving the old
tuple-backed shape, make no source commit and report the exact caller/blocker.

### HP-RAW-SRCMODE-TEST-01 — raw product source-mode inventory validation

Approved validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader;
- H2 R3 endpoint;
- focused raw-product source order and retained-rule parity;
- manifest source-mode and final-basis source-relation inspection;
- focused search confirming `RawProductBoxPlan` no longer stores source-mode
  inventories as `Tuple{Vararg{...}}`;
- no Cr2 run.

Existing committed tests may be adjusted only if they directly assert the old
tuple-backed source-mode inventory shape. No new committed test file, Cr2
fixture, driver-input fixture, benchmark, or route-diagnostic test is approved
by this ID.

### HP-CONTRACT-VEC-FN-01 — contract-plan vector cleanup

Approved source files:

```text
src/cartesian_terminal_lowering/contracts.jl
src/cartesian_terminal_lowering/selection.jl
src/cartesian_terminal_lowering/summaries.jl
src/cartesian_retained_unit_transform_contracts/records.jl
src/cartesian_retained_unit_transform_contracts/unit_contracts.jl
src/cartesian_retained_unit_transform_contracts/summaries.jl
```

Approved narrow consumer files, only as required by the storage change:

```text
src/pqs_source_box_route_driver_helpers.jl
src/cartesian_base_hamiltonian.jl
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
```

Approved cleanup targets:

- `TerminalLoweringPlan.available_contracts::Tuple{Vararg{TerminalLoweringContract}}`;
- `TerminalLoweringPlan.contracts::Tuple{Vararg{TerminalLoweringContract}}`;
- `RetainedUnitTransformContractPlan.contracts::Tuple{Vararg{RetainedUnitTransformContract}}`;
- same-file and listed narrow consumers that currently depend on those
  tuple-backed plan field shapes.

Approved replacement:

- vector-backed terminal-lowering available-contract storage;
- vector-backed terminal-lowering selected-contract storage;
- vector-backed retained-unit transform-contract storage;
- stable accessors preserving ordered contract facts and current behavior:
  `available_contracts(plan)`, `selected_contracts(plan)`, `contracts(plan)`,
  and `transform_contracts(plan)`.

Accessor compatibility means same ordered facts, iteration order, selected
contract behavior, transform-contract behavior, and summaries. It does not
mean preserving variable-length `Tuple` concrete field types or accessor return
types.

This ID does not approve changing
`source_cpbs::Tuple{Vararg{CoordinateProductBox}}`, raw product source-mode
storage, retained-unit route inventories, public input `NamedTuple`s, fixed
coordinate/product-box value objects, numerical kernels, route semantic
changes, shellification behavior changes, public API/export changes, canonical
driver changes, Hamiltonian object changes, matrix-key changes, reader changes,
artifact/manifest schema changes, route-stage objects,
report/status/payload expansion, persistent caches, compatibility layers that
preserve the old tuple-backed plan field types, new committed tests, Cr2 runs,
or Cr2-specific workflow.

Line budget: at most `150` added `src` lines, with net simplification expected.
Failure rule: if vectorizing the plan inventories requires source files outside
the approved surfaces, broad route/stage rewrites, public API or artifact
changes, numerical changes, or compatibility layers preserving the old
tuple-backed plan field types, make no source commit and report the exact
caller/blocker.

### HP-CONTRACT-VEC-TEST-01 — contract-plan vector cleanup validation

Approved validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader;
- H2 R3 endpoint;
- focused terminal-lowering contract order parity;
- focused retained-unit transform-contract order parity;
- focused search confirming targeted plan inventories no longer store
  contracts as `Tuple{Vararg{...}}`;
- no Cr2 run.

Existing committed tests may be adjusted only if they directly assert the old
tuple-backed contract-plan field shape. No new committed test file, Cr2
fixture, driver-input fixture, benchmark, or route-diagnostic test is approved
by this ID.

### HP-ROUTE-STAGE-TYPE-FN-01 — route/stage type-surface cleanup

Approved source files:

```text
src/pqs_source_box_route_driver_helpers.jl
src/cartesian_terminal_shellification_geometry.jl
```

Approved cleanup targets:

- `_pqs_source_box_route_driver_terminal_lowering_contract_inventory_from_plan`;
- `cartesian_units` route/stage return surfaces that carry oversized
  compatibility inventories;
- `_pqs_source_box_route_driver_transform_stage_low_order_summary`;
- `cartesian_transforms` route/stage return surfaces that carry oversized
  compatibility inventories;
- `_cartesian_terminal_shellification_region_unit_inventory`;
- related terminal-region lowering inventory summary surfaces in
  `src/cartesian_terminal_shellification_geometry.jl` only where the same
  runtime-sized type-surface pattern appears.

Approved replacement/deletion shapes:

- delete stale route/stage compatibility inventories with no active approved
  caller;
- replace remaining runtime-sized `NamedTuple` / `Tuple` carriers with
  vector-backed compact internal objects, stable dictionaries, accessors, or
  smaller summaries;
- shrink wide internal stage return signatures only where all live approved
  callers can be updated within the approved source files;
- preserve deterministic terminal shellification/lowering order and existing
  behavior.

Required preservation:

- H2 base artifact/readback behavior;
- H2 supplemented artifact/readback behavior;
- deterministic terminal shellification/lowering order;
- existing public driver contract;
- existing artifact schema and manifest behavior;
- existing numerical matrices.

This ID does not approve source files outside the approved boundary, driver
changes, artifact schema or manifest changes, public API/export changes,
numerical kernel changes, matrix value changes, raw-block changes, Residual
Gaussian/MWG/IDA semantic changes, route semantic changes, shellification
behavior changes, route diagnostic/status/report expansion, broad route-stage
redesign, new public contracts, PackageCompiler/PrecompileTools/sysimage or
precompile workload work, new committed tests, Cr2 runs, or Cr2-specific
workflow. No compatibility adapter may preserve the old runtime-sized type
surface merely under a new name.

Line budget: at most `200` added `src` lines, with net simplification expected.
Failure rule: if cleanup requires source files outside the approved boundary,
broad route-stage redesign, new public contracts, artifact changes, numerical
changes, or a precompile/sysimage mechanism, make no source commit and report
the exact blocker.

### HP-ROUTE-STAGE-TYPE-TEST-01 — route/stage type-surface cleanup validation

Approved validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader;
- H2 R3 endpoint if the pass touches terminal realization behavior;
- focused terminal shellification/lowering order parity;
- focused scan for newly introduced `NamedTuple{...}`, variable-size
  `Tuple(...)`, `Tuple{Vararg{...}}`, and runtime-keyed inventories in the
  approved files;
- optional Be2 q5 compile/timing comparison after correctness passes;
- no Cr2 run.

Existing committed tests may be adjusted only if they directly assert the old
stale compatibility inventory shape. No new committed test file, Cr2 fixture,
driver-input fixture, benchmark, route-diagnostic test, or precompile workload
is approved by this ID.

### HP-ROUTE-STAGE-CARRIER-FN-01 — route/stage carrier cleanup

Approved source files:

```text
src/pqs_source_box_route_driver_helpers.jl
src/pqs_source_box_diatomic_complete_core_shell.jl
```

Optional source file, only if directly required to slim the approved path:

```text
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
```

Approved cleanup targets:

- `cartesian_shells` stage carrier and return signature;
- `cartesian_units` stage carrier and return signature;
- `cartesian_transforms` stage carrier and return signature;
- terminal topology support-region planning in
  `src/pqs_source_box_diatomic_complete_core_shell.jl`;
- terminal retained-rule planning in
  `src/pqs_source_box_diatomic_complete_core_shell.jl`;
- terminal realization plan carriers in
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl` only
  where directly required to avoid reintroducing a large stage-carried plan
  shape through the approved route/stage path.

Approved replacement/deletion shapes:

- stop carrying giant shellification, route-skeleton, support-plan,
  retained-rule-plan, and terminal-plan `NamedTuple` / tuple shapes across the
  approved stage function signatures;
- replace necessary carriers with compact typed/vector-backed records, stable
  dictionaries, accessors, or smaller summaries;
- recompute small derived summaries from canonical objects inside the approved
  path where that is simpler than carrying wide stage payloads;
- delete stale compatibility carriers with no active approved caller;
- preserve deterministic terminal support, shellification, and lowering order.

Route skeleton construction semantics are not changed by this ID. Edits to
`src/pqs_source_box_route_driver_skeletons.jl` are not approved.

Required preservation:

- H2 base artifact/readback behavior;
- H2 supplemented artifact/readback behavior;
- H2 R3 endpoint if terminal realization is touched;
- deterministic terminal support/shellification/lowering order;
- existing public driver contract;
- existing artifact schema and manifest behavior;
- existing route semantics and numerical matrices.

This ID does not approve source files outside the approved boundary, edits to
`src/pqs_source_box_route_driver_skeletons.jl`, driver changes, artifact schema
or manifest changes, public API/export changes, numerical kernel changes,
matrix value changes, raw-block changes, Residual Gaussian/MWG/IDA semantic
changes, route semantic changes, shellification behavior changes, route
diagnostic/status/report expansion, broad route-stage redesign, new public
contracts, PackageCompiler/PrecompileTools/sysimage or precompile workload
work, new committed tests, Cr2 runs, or Cr2-specific workflow. No compatibility
adapter may preserve the old runtime-sized carrier merely under a new name.

Line budget: at most `250` added `src` lines, with net simplification expected.
Failure rule: if cleanup requires source files outside the approved boundary,
broad route-stage redesign, public API changes, artifact changes, numerical
changes, or precompile/sysimage machinery, make no source commit and report the
exact blocker.

### HP-ROUTE-STAGE-CARRIER-TEST-01 — route/stage carrier cleanup validation

Approved validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader;
- H2 R3 endpoint if terminal realization is touched;
- focused terminal support/shellification/lowering order parity;
- focused scan for newly introduced runtime-sized `NamedTuple{...}`,
  `Tuple(...)`, `Tuple{Vararg{...}}`, and runtime-keyed inventories in the
  approved files;
- optional Be2 q5 post-cleanup compile/timing comparison after correctness
  passes;
- no Cr2 run.

Existing committed tests may be adjusted only if they directly assert the old
stale carrier shape. No new committed test file, Cr2 fixture,
driver-input fixture, benchmark, route-diagnostic test, or precompile workload
is approved by this ID.

### HP-R1-ART-01 — public base producer artifact provenance

Approved artifact extension for R1 public facade writes only. When
`hamfile !== nothing`, the final Hamiltonian file may add the following JLD2
keys under `producer_provenance/`:

```text
provenance_version
producer
nesting
route
ns
q
q_rule
ns_source
core_spacing
reference_spacing
tail_spacing
parent_axis_family
parent_axis_counts
mapping_kind
mapping_d
radius
xmax_parallel
xmax_transverse
atom_symbols
nuclear_charges
atom_locations
nup
ndn
final_dimension
```

Exact values and route-specific `nothing` fields are defined in
`r1_public_base_producer.md`. For one-center H, the provenance must record
`core_spacing = 0.3`, `reference_spacing = 1.0`,
`mapping_kind = :white_lindsey_atomic_mapping`, and
`mapping_d = 0.3` for the reviewed endpoint. The `mapping_d` provenance value
is the resolved internal White-Lindsey mapping parameter and equals
`core_spacing` for one-center atoms; it is not a separate public input.
`nesting` must record the public construction-family input, and `route` must be
the truthful base route label derived from `(input.kind, input.nesting)`, not a
PQS-oriented default string. `ns` records the normalized public requested
source size. `q` records the derived route-local value consumed by route
construction, not a common public cube-size knob. `q_rule` records the
derivation rule, and `ns_source` records whether `ns` came from the new public
field or a temporary legacy-`q` compatibility path.
Existing
`read_cartesian_ida_hamiltonian` must continue reading the Hamiltonian matrices
while ignoring these extra keys. This ID does not approve a separate manifest,
provenance file, wrapper object, public provenance reader, status/report
payload, or use of provenance as staged algorithm data after initial
lattice/parent construction.

### HP-R1-TEST-01 — public base producer endpoint test/example

Approved committed validation surface for R1.

Approved standalone integration gate:

```text
test/driver_public/cartesian_base_hamiltonian_runtests.jl
```

Invocation:

```text
julia --project=. test/driver_public/cartesian_base_hamiltonian_runtests.jl
```

The test/example should exercise the public facade for one-center H and
z-axis H2, verify `CartesianIDAHamiltonian{Float64}` output, validate the
reviewed H baseline using explicit public `core_spacing = 0.3` with
`reference_spacing = 1.0` and H2 endpoint facts, validate unknown-key and
malformed input errors including mismatched temporary H `d` and rejected H2
`d`, validate
x/y-aligned, shifted-parallel, and generally oriented H2 rejection before
expensive construction, and validate existing Hamiltonian artifact
write/readback plus `HP-R1-ART-01` provenance using `mktempdir()`. It is a
standalone integration/endpoint gate, not ordinary tiny unit coverage, and this
ID does not approve adding it to `test/runtests.jl`. It must not assert private
route-stage fields, report mirrors, status/blocker symbols, terminal role
vocabulary, or pair inventories.

## Approved For R1 One-Center Base Atoms

This section approves only the explicit origin-centered one-center
all-electron atom relaxation recorded in `r1_one_center_base_atoms.md`. It
extends the existing base facade scope without adding a new public function,
new export, new artifact schema, new route vocabulary, or supplemented atom
authority.

### HP-R1-ATOM-FN-01 — explicit one-center all-electron base atom facade

Approved source file:

```text
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- accept exactly one origin-centered atom in the existing
  `cartesian_base_hamiltonian(system; basis, hamfile)` call shape;
- require explicit vector-valued `atom_symbols`, `nuclear_charges`,
  `atom_locations`, and explicit integer `nup`, `ndn`;
- require finite positive integer-valued nuclear charge supplied by the caller;
- require neutral all-electron count
  `nup + ndn == round(Int, only(nuclear_charges))`;
- treat the atom symbol as provenance/user labeling only, not as a source of
  charge, spin, basis, or ECP defaults;
- keep required one-center basis fields `ns`, `core_spacing`, and `radius`
  after `HP-COMP-NS-*` normalization;
- treat public `d`, if temporarily accepted, as a deprecated compatibility
  alias that must equal resolved `core_spacing`.

This ID does not approve translated atoms, element lookup/default tables,
inferred charge or spin, ECP, solver workflow, supplemented atom construction
under the base facade, public API redesign, or new artifact fields. Supported
supplemented one-center atoms are governed separately by `HP-COMP-SUPPATOM-*`.

### HP-R1-ATOM-WIRE-01 — one-center atom shared workflow wiring

Approved source file:

```text
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- map public `only(system.nuclear_charges)` to the existing private
  White-Lindsey atomic mapping `Z`;
- map the resolved public `core_spacing` to the private White-Lindsey
  `parent_mapping_d`;
- keep `reference_spacing`, `tail_spacing`, and box/domain controls separate
  from `core_spacing`;
- feed atom geometry/shellification normalization into the same terminal-basis,
  one-body, IDA, `CartesianIDAHamiltonian`, artifact-writing, and provenance
  machinery used by the base producer;
- preserve existing `HP-R1-ART-01` `producer_provenance/` keys with
  `route = :one_center_pqs_base`.

Atoms and diatomics must share the same producer workflow after the narrow
geometry/shellification differences. This ID does not approve an atom-only
Hamiltonian builder, parallel atom materialization path, atom route-stage
object, atom report/status/payload object, or metadata/provenance carrier used
as algorithmic data.

Line budget for `HP-R1-ATOM-FN-01` plus `HP-R1-ATOM-WIRE-01`: at most `80`
added `src` lines. If implementation needs source edits outside
`src/cartesian_base_hamiltonian.jl`, changes to private materialization
owners, atom-only materialization, new artifact keys, translated atoms,
ECP behavior, solver workflow, element lookup/default tables, committed
fixtures/tests, route/report/status/payload expansion, or supplemented atom
work outside `HP-COMP-SUPPATOM-*`, stop and request a new docs-only amendment.

### HP-R1-ATOM-TEST-01 — one-center base atom validation

Approved validation:

- existing origin-centered H public facade endpoint remains unchanged,
  now expressed as `core_spacing = 0.3`, `reference_spacing = 1.0`, and
  internal `parent_mapping_d = core_spacing`;
- optional ignored/user-run Be or Cr one-center base atom artifact
  write/readback using explicit charge, spin sectors, origin geometry, and
  basis controls;
- finite/symmetric `K`, unit `U_A`, and IDA `V` for ignored/user-run non-H
  atom checks;
- clear `ArgumentError` for translated atom input, mismatched temporary `d`,
  noninteger or nonpositive charge, nonneutral electron count, or
  element-table/default requests where practical.

No new committed test file, committed non-H atom fixture, public non-H
reference scalar, solver run, supplemented atom endpoint under this base-atom
lane, ECP gate, translated-atom gate, or driver change is approved by this ID.
Supported supplemented atom validation is governed by
`HP-COMP-SUPPATOM-TEST-01`.

## Approved For R3/RG Implementation

The R3 labels remain approved compatibility and endpoint-history IDs. Current
Residual Gaussian algorithm authority lives in
`residual_gaussian_domain_module.md`; do not copy that algorithm into this
registry. This section records approved IDs, source owners, function surfaces,
artifact keys, and validation gates.

### R3 Compatibility Boundary

Approved R3 owner/path for compatibility entry points, artifact writing, and the
non-exported usability facade:

```text
Owner module: CartesianFinalBasisRealization
Source file: src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

R3 compatibility surfaces must delegate current residual-basis, exact-operator,
and residual-interaction physics to the Residual Gaussian module where those
module functions exist. They must not preserve retired global raw-candidate
selection, old density gauges, or duplicate RG algorithms.

### HP-R3-OBJ-01 — residual-GTO augmentation object

Approved historical/object authority for the first H2 residual-GTO endpoint.
Current domain object fields and semantics are recorded under `HP-RG-OBJ-01` and
`residual_gaussian_domain_module.md`. Compatibility names may remain only where
live callers still use them.

### HP-R3-FN-01 — residual-basis construction

Approved historical R3-A basis-construction surface. Current production logic is
`build_residual_gaussian_basis(...)` under `HP-RG-FN-01`.

Binding guardrails: residual basis directions are selected separately on each
physical owner atom; residual occupation is not numerical rank; owner-local
sectors are merged once; global raw-candidate Lowdin and global raw-column
pivoted-Cholesky selection are not approved current algorithms.

### HP-R3-FN-02 — exact augmented one-body and moment assembly

Approved historical R3-A exact-operator surface. Current production logic is
`transform_augmented_operator(...)` under `HP-RG-FN-02`.

The exact transformed operators are kinetic, every uncharged by-center nuclear
attraction, `x`, `y`, `z`, `x^2`, `y^2`, and `z^2`. This exact transformation is
not the MWG approximation and must not be replaced by moment-matched interaction
logic.

### HP-R3-FN-03 — residual MWG/IDA and in-memory Hamiltonian

Approved R3-B compatibility entry point:

```text
pqs_terminal_residual_gto_augmented_hamiltonian(...)
```

Output is the existing `CartesianIDAHamiltonian{Float64}`. Current residual MWG
descriptor and interaction math is owned by `moment_matched_gaussians(...)` and
`assemble_residual_ida_interaction(...)` under `HP-RG-FN-03` and `HP-RG-FN-04`.

The accepted H2 owner-local endpoint has augmented dimension `489` and
lowest-orbital IDA self-Coulomb `0.4574265214362075` within `1.0e-10`. Older
R3-B scalars from global-selection or retired density-gauge diagnostics are
historical evidence only and are not active targets.

### HP-R3-ART-01 — compact supplemented artifact provenance

Approved R3-C internal artifact helper remains outside the RG module. It may
write an existing `CartesianIDAHamiltonian` file and add compact
`supplement_provenance/` provenance. RG does not own artifact writing,
artifact schema, JLD2 workflow, or provenance readers.

Approved source owner/path:

```text
Owner module: CartesianFinalBasisRealization
Source file: src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

Approved `supplement_provenance/` keys:

- `provenance_version`;
- `producer = :cartesian_residual_gto_mwg_augmentation`;
- `supplement_policy = :mwg_residual_gto`;
- `basis_by_center`;
- `lmax`;
- `uncontracted`;
- `width_filtering`;
- `candidate_count`;
- `owner_counts`;
- `base_dimension`;
- `residual_dimension`;
- `augmented_dimension`;
- `augmented_basis_order = :base_then_residual`;
- `residual_basis_convention = :owner_local_residual_occupation_final_merge_lowdin`;
- `rank_rule` / owner-local selection rule;
- `occupation_cutoff = 1.0e-6`;
- `tau_neg_abs`, `tau_neg_rel`;
- `tau_merge_abs = 1.0e-12`, `tau_merge_rel = 1.0e-12`;
- `mwg_convention_version`;
- `mwg_convention = :separable_moment_matched_density_normalized`;
- `one_body_source`;
- `interaction_source = :weight_aware_residual_mwg_ida_blocks`;
- compact validation labels and H2 reference value when supplied.

Do not serialize full residual bases, dense moments, `T_G`, `T_A`, MWG centers,
MWG widths, or broad construction inputs in this compact group.

## Approved For Compact Hamiltonian Artifact Manifest

This section approves only artifact sidecar groups for existing
`CartesianIDAHamiltonian{Float64}` JLD2 files. It does not approve a new
Hamiltonian object, new matrix keys, public reader API, driver public input
change, solver workflow, CR2-consumer-specific field, Cr2-specific field,
report/status payload, or artifact schema dump in the driver.

### HP-HAM-MANIFEST-FN-01 — compact Hamiltonian artifact manifest

Approved source files:

```text
src/cartesian_base_hamiltonian.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/cartesian_ida_hamiltonian.jl
```

`src/cartesian_ida_hamiltonian.jl` is approved only for a small unexported
sidecar writer/helper if needed; it must not change
`write_cartesian_ida_hamiltonian` matrix keys or
`read_cartesian_ida_hamiltonian` behavior.

Approved sidecar groups:

```text
hamiltonian_manifest/
hamiltonian_manifest/final_basis_labels/
hamiltonian_manifest/final_basis_source_relations/
hamiltonian_manifest/source_shells/
hamiltonian_manifest/source_modes/
recipe_provenance/
```

`hamiltonian_manifest/` reuses the source/fixed-column provenance model from
`docs/src/developer/projected_q_shell_policy.md`. Basis identity is a
status-bearing construction label, not `center_xyz`. Centers are
representative metadata only.

Approved `hamiltonian_manifest/` root key:

- `manifest_version = 1`;

Approved required `hamiltonian_manifest/final_basis_labels/` fields:

- `final_basis_col`;
- `sector`;
- `unit_label`;
- `unit_kind`;
- `source_region_label`;
- `source_region_label_status`;
- `source_box_label`;
- `source_box_label_status`;
- `owner_nucleus_index`;
- `owner_label_status`;
- `shell_label_status`;
- `shell_index`;
- `ray_label_status`;
- `ray_id`;
- `ray_family_label`;
- `radial_order_status`;
- `radial_order`;
- `center_x`;
- `center_y`;
- `center_z`;
- `center_definition`;
- `center_status`;
- `lowdin_correction_applied`;
- `supplement_label`;
- `angular_power_x`;
- `angular_power_y`;
- `angular_power_z`;
- `inferred_from_centers`;
- `inferred_from_nearest_grid`;
- `inferred_from_support_order`;
- `inferred_from_support_indices`;
- `inferred_from_raw_to_final_support`.

The final-basis label rows must follow the exact matrix row/column order.
`owner_nucleus_index` uses one-based physical nucleus indices and `0` when no
owner is meaningful; unavailable integer shell/ray/radial labels use `0`;
unavailable angular powers use `-1`; and all `inferred_from_*` flags must be
`false` for production manifest rows. Approved sectors are `:base`,
`:residual`, and `:supplement_derived`.

Approved optional `hamiltonian_manifest/final_basis_source_relations/` fields:

- `final_basis_col`;
- `relation_index`;
- `relation_kind`;
- `source_shell_id`;
- `source_mode_label`;
- `local_axis_x`;
- `local_axis_y`;
- `local_axis_z`;
- `relation_status`;
- `shell_label_status`;
- `ray_label_status`;
- `radial_order_status`;
- `coefficient_status`;
- `weight_status`;
- `span_status`;
- `inferred_from_centers`;
- `inferred_from_nearest_grid`;
- `inferred_from_support_order`;
- `inferred_from_support_indices`;
- `inferred_from_raw_to_final_support`.

Relations, `source_shells/`, and `source_modes/` may be populated only when
the construction natively defines those facts. Source-mode identity is
`(source_shell_id, local_axis_x, local_axis_y, local_axis_z)` in shell-local
coordinates. It is a label, not a coefficient map, support row, parent row, or
operator reconstruction claim. Missing shell, ray, radial, source-box, or
relation facts must be status-bearing `:unavailable` or `:mixed`, not inferred
from centers, nearest-grid snapping, support order, support indices, or
raw-to-final support.

Approved `recipe_provenance/` keys:

- `provenance_version = 1`;
- `producer`;
- `nesting`;
- `route`;
- `ns`;
- `q`;
- `q_rule`;
- `ns_source`;
- `core_spacing`;
- `padding`;
- `radius`;
- `xmax_parallel`;
- `xmax_transverse`;
- `extent_source`;
- `parent_axis_counts`;
- `atom_symbols`;
- `nuclear_charges`;
- `atom_locations`;
- `nup`;
- `ndn`;
- `basisname`;
- `basisfile`;
- `lmax`;
- `uncontracted`;
- `width_filtering`;
- `base_dimension`;
- `residual_dimension`;
- `augmented_dimension`.

The recipe group may repeat facts from `producer_provenance/` and
`supplement_provenance/` so consumers have one uniform location. Values must
come from the validated public construction contract and produced dimensions,
not route reports, element tables, solver assumptions, or private diagnostics.
`nesting` records the public construction family (`:pqs` or `:wl`), and
`route` records the truthful base route label derived from `(input.kind,
input.nesting)`. `ns` records the normalized public requested source size.
`q` records the derived route-local value, not a common public cube-size knob.
`q_rule` records the derivation rule, and `ns_source` records whether `ns` came
from the new public field or a temporary legacy-`q` compatibility path.

Center conventions and construction labels must be derived from existing
terminal basis blocks, parent axes, residual metadata, and augmented moment/MWG
descriptors. If the implementation cannot derive a required center convention
or source label from those existing objects without adding algorithmic
metadata, it must stop and report the missing seam.

This ID does not approve `T_G`, `T_A`, dense residual transforms, coefficients,
dense moment matrices, raw inventories, allocation probes, report/status
payloads, public reader APIs, public exports, driver public input changes,
artifact schema dumps in the driver, solver-specific fields,
CR2-consumer-specific fields, Cr2-specific fields, committed Cr2 fixtures, or
Cr2-specific branches.

One-center atom padding is provenance-only for this manifest lane. Source work
under `HP-HAM-MANIFEST-*` must not change atom size policy or parent-axis
counts. Production atom parent sizing is separately governed by
`HP-COMP-ATOMBOX-*`; diatomic padding-derived extents continue to use the
existing facade contract.

Line budget: at most `150` added `src` lines. Stop for a separate amendment if
the implementation needs source files outside the approved surfaces, new
algorithmic metadata, a public reader, driver contract changes, or changes to
the Hamiltonian matrix writer/reader contract.

### HP-HAM-MANIFEST-TEST-01 — artifact manifest validation

Approved validation:

- `git diff --check`;
- package load;
- H atom or H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader;
- direct JLD2 checks that `hamiltonian_manifest/final_basis_labels/` rows
  match matrix dimension and approved status-bearing fields;
- direct JLD2 checks that unavailable or mixed shell/ray/radial/source labels
  are explicit and that no production row is marked inferred from centers,
  nearest grid, support order, support indices, or raw-to-final support;
- direct JLD2 checks that `recipe_provenance/` records validated public
  system, basis, supplement, route, parent-axis counts, and dimensions;
- optional practical Be2 supplemented artifact manifest inspection;
- no Cr2 run.

No new committed test file, public reader API, artifact schema dump, driver
public input change, Cr2 fixture, or solver/CR2 workflow validation is approved
by this ID.

### HP-HAM-MANIFEST-SRC-FN-01 — source-mode provenance seam

Approved purpose: carry compact construction-native source-mode provenance
from terminal lowering / retained-unit / raw-product source planning to the
base working basis manifest context so artifact writing can populate optional
source provenance groups without route reports or center inference.

Approved source files:

```text
src/cartesian_terminal_lowering/contracts.jl
src/cartesian_terminal_lowering/region_contracts.jl
src/cartesian_raw_product_sources/CartesianRawProductSources.jl
src/cartesian_raw_product_sources/records.jl
src/cartesian_raw_product_sources/source_mode_indices.jl
src/cartesian_retained_units/CartesianRetainedUnits.jl
src/cartesian_retained_units/records.jl
src/cartesian_retained_units/lower_contract_units.jl
src/cartesian_retained_unit_transform_contracts/CartesianRetainedUnitTransformContracts.jl
src/cartesian_retained_unit_transform_contracts/records.jl
src/cartesian_retained_unit_transform_contracts/unit_contracts.jl
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
src/cartesian_base_hamiltonian.jl
```

Module wrapper files are approved only for internal exports/includes required
by the compact provenance seam. They are not public API authority.

Approved carrier:

- preferred: one internal `source_mode_provenance` field on the
  `cartesian_base_working_basis(...)` return value;
- allowed if cleaner: one optional compact source-mode provenance field on
  `CartesianTerminalBasisRealization`;
- no other stage object, report, status payload, route summary, artifact
  wrapper, or persistent cache is approved.

Approved source facts:

- source shell IDs, unit links, source-box/source-region labels, construction
  kind, source intervals, source-mode dimensions, source-mode ordering, center
  definition/status, Lowdin-correction status, and shell/ray/radial label
  statuses;
- source mode identities
  `(source_shell_id, local_axis_x, local_axis_y, local_axis_z)`, local
  source-mode ordering, native parent-lattice coordinates only when already
  available, representative center metadata, and status flags;
- final-basis source-relation rows only where the relation is construction
  native: direct identity, boundary source-mode selection, product-axis tuple,
  or explicit `:mixed` / `:unavailable`;
- final-basis label improvements only where the final basis column is directly
  and natively tied to the row's unit/source mode.

Ray, cone, shell, and radial labels may be written only when already natively
defined by the construction. This ID does not approve a repo-chosen ray/cone
grouping policy or inferred labels from centers, nearest-grid snapping, support
order, support indices, or raw-to-final support.

This ID does not approve coefficients, dense transforms, `T_G`, `T_A`, raw
candidate inventories, raw pair inventories, allocation probes, route reports,
diagnostic payloads, metadata/status field clouds beyond the compact approved
provenance object, Hamiltonian object changes, matrix-key changes,
`read_cartesian_ida_hamiltonian` changes, public API/export changes, driver
changes, artifact schema dumps, solver fields, CR2-consumer-specific fields,
Cr2-specific fields, committed Cr2 fixtures, or Cr2-specific branches.

Line budget: at most `180` added `src` lines. Stop for a separate amendment if
the implementation needs new source files, public exports, dense payloads,
driver changes, reader changes, route report plumbing, or a source-mode/ray
policy not already present in construction metadata.

### HP-HAM-MANIFEST-SRC-TEST-01 — source-mode provenance seam validation

Approved validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader;
- direct JLD2 checks that optional `source_shells/`, `source_modes/`, and
  `final_basis_source_relations/` are present only for construction-native
  provenance rows;
- checks that absent shell/ray/radial/source facts are explicitly
  `:unavailable` or `:mixed` and not inferred;
- optional practical Be2 supplemented manifest inspection;
- no Cr2 run.

No new committed test file, public reader API, driver public input change,
artifact schema dump, Cr2 fixture, solver/CR2 workflow validation, or broad
route/report validation is approved by this ID.

### HP-NEST-ART-FN-01 — nesting artifact-truth cleanup

Approved source files:

```text
src/cartesian_base_hamiltonian.jl
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
```

`src/cartesian_base_hamiltonian.jl` is approved only for artifact provenance
truth for the public `nesting` construction-family input. It may:

- record `nesting` in `producer_provenance/` and `recipe_provenance/`;
- choose truthful base route labels from `(input.kind, input.nesting)`, with
  approved labels `:one_center_pqs_base`, `:one_center_wl_base`, and
  `:z_axis_diatomic_pqs_base`;
- historically enforced a supplemented-WL early stop under this narrow
  provenance lane; that stop is superseded for the supported z-axis diatomic
  supplemented WL composition cell by `HP-COMP-SUPPWL-*`;
- leave unsupported WL H2 without a provenance label until that path succeeds
  under separate authority.

`src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl` is
approved only for a docstring correction so the module description no longer
says it is exclusively PQS-specific now that the WL terminal-basis seam uses the
same final-basis boundary.

This ID does not approve driver public input changes, route skeleton changes,
shellification changes, terminal lowering changes, raw-block changes,
Residual Gaussian/MWG/IDA changes, artifact matrix-key changes,
`read_cartesian_ida_hamiltonian` behavior changes, public API/export changes,
Cr2 workflow, committed tests, diagnostic/report changes, or WL H2 support.

Failure rule: if truthful nesting provenance requires changing reader behavior,
artifact matrix keys, or the broader manifest structure, make no source commit
and report the blocker.

### HP-NEST-ART-TEST-01 — nesting artifact-truth validation

Approved validation:

- `git diff --check`;
- package load;
- small `nesting = :pqs` base artifact write/readback plus direct provenance
  inspection;
- small `nesting = :wl` one-center atom artifact write/readback plus direct
  provenance inspection;
- the historical supplemented-WL rejection from this lane is superseded by
  `HP-COMP-SUPPWL-*` for the supported z-axis diatomic supplemented WL cell;
- no Cr2 run.

No new committed test file, driver-input fixture, public reader API, artifact
schema dump, WL H2 validation, supplemented WL validation, or Cr2 fixture was
approved by this ID. Later supplemented WL validation is owned by
`HP-COMP-SUPPWL-TEST-01`.

### HP-R3U-FILE-01 — supplemented workflow source and validation files

Approved non-exported usability owner:

```text
src/cartesian_base_hamiltonian.jl
```

Allowed companion surfaces:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
```

No public export, new source file, new committed test file, driver/bin/tool
workflow, report/status/payload object, or artifact shape beyond
`supplement_provenance/` is approved.

### HP-R3U-FN-01 — non-exported supplemented Hamiltonian facade

Approved internal call shape:

```julia
cartesian_residual_gto_mwg_hamiltonian(
    system::NamedTuple;
    basis::NamedTuple,
    supplement::NamedTuple,
    hamfile::Union{Nothing,AbstractString} = nothing,
)::CartesianIDAHamiltonian{Float64}
```

Original first systems were z-axis H2 and z-axis Be2. `HP-R3U-ZDI-FN-01`
relaxes that guard to explicit homonuclear two-center z-axis diatomics. Be2
remains an internal/performance-supported proxy. Cr2 is permitted only as an
explicit generic homonuclear z-axis ignored/user-run stress or usability case
after H2/Be2 validation. Heteronuclear systems, non-z-axis or arbitrary
orientations, charged systems, ECP inputs, solver/RHF workflow, public export,
and Cr2-specific branches remain unapproved.

### HP-R3U-WIRE-01 — base-to-RG same-construction workflow

Approved wiring:

```text
validated system/basis/supplement spec
-> R1-style/base normalization and base stages
-> base CartesianIDAHamiltonian plus same-construction terminal basis/bundles
-> legacy named-basis supplement loading
-> basis_representation(supplement)
-> RG/R3 same-construction augmented Hamiltonian path
-> optional R3-C artifact writer
-> CartesianIDAHamiltonian{Float64}
```

The base Hamiltonian, terminal basis realization, and parent axis bundle must
come from the same base construction call. The facade must not expose terminal
basis realizations, bundles, residual objects, augmented-operator objects, MWG
descriptors, pair factors, or provenance payloads.

### HP-R3-TEST-01 / HP-R3U-TEST-01 — standalone H2 endpoint gate

Approved standalone validation file:

```text
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
```

Invocation:

```text
julia --project=. test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
```

The gate covers the H2 residual-GTO/MWG endpoint family, including residual
basis checks, exact augmented one-body/moment checks, independent weight-aware
`V_GM` comparison, in-memory Hamiltonian checks, and the non-exported usability
facade/artifact readback section. It is not approved for `test/runtests.jl`,
Be2 committed validation, Cr2 validation, or private route/status assertions.

## Approved For Residual Gaussian Domain Migration

Current Residual Gaussian domain algorithm authority is
`residual_gaussian_domain_module.md`. This registry records the approved module
files and function surfaces only.

The Residual Gaussian module does not own compact supplemented artifact writing
or `supplement_provenance/`. Artifact/facade hooks remain outside RG unless a
later amendment names a real duplication or consumer reason.

### HP-RG-FILE-01 — Residual Gaussian module files

Approved internal module and files:

```text
src/cartesian_residual_gaussians/CartesianResidualGaussians.jl
src/cartesian_residual_gaussians/residual_basis.jl
src/cartesian_residual_gaussians/augmented_operators.jl
src/cartesian_residual_gaussians/mwg_interaction.jl
```

No public export is approved.

### HP-RG-OBJ-01 — residual Gaussian basis object

Approved domain object: a numerical residual Gaussian basis object carrying base
dimension, candidate count, residual dimension, candidate owner indices,
residual source owner indices, owner retained counts, residual occupations,
cutoff/tolerance policy, selection/orientation/sign rules, and final
`T_G::Matrix{Float64}` / `T_A::Matrix{Float64}` transforms.

It must not be a status/result payload and must not carry route metadata,
report fields, status flags, artifact data, MWG descriptors, or public API
state.

### HP-RG-FN-01 — residual Gaussian basis construction

Approved production name:

```julia
build_residual_gaussian_basis(...)
```

Candidate owner indices are required. The current algorithm is the owner-local
residual occupation construction defined in
`residual_gaussian_domain_module.md`.

### HP-RG-FN-02 — exact augmented operator transformation

Approved production name:

```julia
transform_augmented_operator(...)
```

This owns exact `[G,A] -> [G,R]` transformation for kinetic, uncharged
by-center nuclear attraction, and first/second Cartesian moment matrices.

### HP-RG-FN-03 — moment-matched Gaussian descriptors

Approved production name:

```julia
moment_matched_gaussians(...)
```

Descriptors are computed from final merged residual functions and exact moment
matrices. They are residual-containing interaction descriptors only, not exact
residual-GTO Coulomb integrals.

### HP-RG-FN-04 — residual IDA interaction assembly

Approved production name:

```julia
assemble_residual_ida_interaction(...)
```

This owns residual-containing MWG/IDA blocks `V_GM` and `V_MM`, combined with
unchanged base `V_GG`. `V_GM` uses weight-aware final-basis density
normalization for PQS shell blocks.

### HP-RG-WIRE-01 — migration from terminal residual file

`src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` may keep
small compatibility wrappers for live callers and artifact/facade hooks. Moved
physics helpers should be deleted from that file once callers use RG-domain
helpers directly.

### HP-RG-TEST-01 — migration validation

Approved validation is the existing standalone H2 residual-GTO/MWG endpoint
with augmented dimension `489`, self-Coulomb `0.4574265214362075`, exact
one-body/moment checks, independent weight-aware `V_GM` check, and optional
ignored Be2 usability/performance measurement when a source pass changes the
interaction path or facade wiring.

No new committed test file, Be2 committed gate, or Cr2 full
Hamiltonian/artifact/facade validation is approved.

### HP-RG-ORTHO-FN-01 — residual final-orthogonality robustness

Approved source files:

```text
src/cartesian_residual_gaussians/residual_basis.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

`residual_basis.jl` is the primary owner. The terminal residual file is
approved only for narrow existing internal keyword plumbing if an approved
tolerance/check option must be passed through a compatibility entry point.

Approved target: make final residual orthogonalization and final
`R' S R - I` identity validation robust for small floating-point overshoots
with healthy owner-local selection and final merge spectra. The motivating
strict N2 q5 p10 case at `core_spacing = 0.042857` has `max |G' S R| =
1.776e-14`, `max |R' S R - I| = 1.673e-10`, retained counts `9,9`, and final
merge eigenvalues `7.232e-2 .. 1.928`.

Allowed changes:

- stable symmetric final merge normalization/check;
- explicitly symmetrized final residual-overlap validation;
- combined absolute/relative final identity check
  `err_RR <= 1.0e-10 + 1.0e-10 * max(1, scale_RR)`;
- no public API unless the option is already routed through internal keywords.

This ID does not approve blind broad tolerance relaxation, residual-selection
semantic changes, global residual selection, occupation-cutoff changes,
negative-eigenvalue tolerance changes, final merge eigenvalue flooring, width
filtering as a conditioning repair, MWG/IDA/nuclear/raw-block changes, artifact
schema changes, driver changes, status/report fields, public API/export
changes, new committed tests, Cr2 workflow, or source files outside the two
approved files.

Failure rule: if the strict N2 case requires changing residual selection,
supplement construction, or final-basis construction, make no source commit and
report the blocker.

### HP-RG-ORTHO-TEST-01 — residual final-orthogonality validation

Approved validation:

- existing H2 residual-GTO/MWG endpoint;
- H2 base/supplemented readback if touched through the facade or compatibility
  file;
- ignored strict N2 q5 p10 residual audit or artifact smoke at
  `core_spacing = 0.042857`;
- one passing N2 comparison at `core_spacing = 0.05` or `0.075`;
- report `max |G' S R|`, `max |R' S R - I|`, retained owner counts, and final
  merge eigenvalue range/condition.

No committed fixture/test, Cr2 full Hamiltonian, Cr2 artifact, Cr2 facade
support, driver workflow, artifact schema change, solver/RHF, ECP, or EGOI
work is approved.

### HP-RG-IDTOL-FN-01 — residual final-identity tolerance default

Approved source files:

```text
src/cartesian_residual_gaussians/residual_basis.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

`residual_basis.jl` is the primary owner. The terminal residual file is
approved only for narrow compatibility keyword default plumbing if needed.

Approved target: set the default final residual `R' S R` identity validation
tolerance to `1.0e-8` when owner-local selection, final merge metric checks,
and `G' S R` orthogonality remain healthy. This is a final
validation/cleanup tolerance only. It is not a residual direction-selection
criterion.

Evidence: Be atom cc-pV5Z `lmax = 1`, `ns = 5`,
`core_spacing = 0.075`, `radius = 20`, `Z = 4`, `nup = 2`, `ndn = 2`
failed only with `max |R' S R - I| = 2.183e-10` against an allowed error of
about `2.000e-10`. The same run had retained residual count `21`, minimum
retained occupation `6.151e-6`, final merge condition `1.0`, and
`max |G' S R| = 1.776e-14`.

Required policy:

- keep the then-current `residual_occupation_cutoff = 1.0e-8` for the Be
  tolerance pass only; this older production default is superseded by
  `HP-RG-CUTOFF-FN-01`;
- keep width/zeta filtering explicit and user-controlled;
- keep owner-local metric checks, final merge metric checks, and `G' S R`
  orthogonality checks active;
- do not drop retained directions to pass final identity validation.

This ID does not approve driver changes, artifact schema/provenance/reader/
manifest changes, residual-selection algorithm changes, width/zeta filtering
default changes, owner grouping changes, merge metric failure-rule changes,
MWG/IDA convention changes, Gaussian raw-block changes, terminal-basis changes,
WL/PQS route changes, shellification changes, Hamiltonian assembly changes,
committed tests/fixtures, Cr2 workflow, or source files outside the two
approved files.

Failure rule: if Be cc-pV5Z cannot pass by changing only the final identity
tolerance default, make no source commit and report the exact blocker.

### HP-RG-IDTOL-TEST-01 — residual final-identity tolerance validation

Approved validation:

- Be atom cc-pV5Z `lmax = 1` residual audit/artifact path passes with the same
  `21` retained residual directions;
- Be atom cc-pVDZ `lmax = 1` still passes;
- H2 residual-GTO/MWG endpoint remains unchanged;
- report `max |R' S R - I|`, allowed tolerance, retained count, minimum
  retained occupation, final merge condition, and `max |G' S R|`;
- no Cr2 run.

No committed fixture/test, driver workflow, artifact schema change,
solver/RHF, ECP, EGOI, Cr2 full Hamiltonian, Cr2 artifact, or Cr2 facade
support is approved.

### HP-RG-CUTOFF-FN-01 — residual occupation cutoff and identity tolerance defaults

Approved source files:

```text
src/cartesian_residual_gaussians/residual_basis.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

`residual_basis.jl` is the primary owner. The terminal residual file is
approved only for narrow compatibility keyword default plumbing if needed.

Approved target: set the default Residual Gaussian owner-local residual
selection cutoff and final identity validation tolerance to:

```text
residual_occupation_cutoff = 5.0e-8
identity_atol = 5.0e-8
```

This was the production default after the Cr atom marginal-direction pass. The
residual occupation cutoff is superseded by `HP-RG-CUTOFF-FN-02`; the final
identity validation tolerance remains `5.0e-8`.

Evidence: Cr atom PQS `basis_ns = 9`, `map_ns = 11`, `lmax = 1` retained a
marginal residual direction at occupation `3.637e-8`. If the production policy
is to discard such marginal residual directions, the default cutoff must say
so explicitly rather than preserving the direction through the older
`1.0e-8` cutoff.

Approved behavior for this now-superseded default:

- discard owner-local residual directions below `5.0e-8` by default;
- keep the final `R' S R` identity validation default aligned at
  `identity_atol = 5.0e-8`;
- preserve explicit caller overrides where already supported;
- keep width/zeta filtering explicit and user-controlled;
- keep owner-local metric checks, negative-eigenvalue tolerances, final merge
  metric checks, and `G' S R` orthogonality checks active;
- keep owner grouping, final merge failure rules, MWG/IDA conventions,
  artifact schema/provenance/reader/manifest, driver workflow, and public API
  unchanged.

This ID supersedes the older `1.0e-8` production defaults recorded under
`HP-RG-IDTOL-FN-01`. It does not change the residual-selection algorithm; it
changes the default retained-occupation policy.

Forbidden:

- residual selection algorithm changes;
- global raw-candidate residual selection;
- global raw-column pivoted-Cholesky residual selection;
- owner grouping changes;
- negative-eigenvalue tolerance changes;
- final merge metric failure-rule changes or merge eigenvalue flooring;
- width/zeta filtering default changes or width filtering as conditioning
  repair;
- MWG/IDA, nuclear, raw-block, exact-operator, terminal-basis, WL/PQS route,
  shellification, Hamiltonian assembly, artifact schema, reader, manifest,
  driver, public API/export, solver/RHF, ECP, EGOI, or Cr2 workflow changes;
- committed fixtures/tests except the exact existing H2 endpoint assertion
  update named under `HP-RG-CUTOFF-TEST-01`.

Failure rule: if the Cr atom case cannot pass or cleanly drop the marginal
direction by changing only the two approved defaults, make no source commit in
the later implementation pass and report the exact blocker.

### HP-RG-CUTOFF-TEST-01 — residual cutoff/tolerance validation

Approved validation:

- Cr atom PQS `basis_ns = 9`, `map_ns = 11`, `lmax = 1` residual construction
  passes or cleanly drops the marginal `s4` direction at occupation
  `3.637e-8` as intended;
- Be atom cc-pV5Z still passes;
- H2 residual-GTO/MWG endpoint remains unchanged;
- exactly update the existing committed H2 endpoint test
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` so both cutoff
  assertions expect `5.0e-8` instead of `1.0e-8`: the in-memory
  `residual.occupation_cutoff` assertion and the artifact/provenance
  `values[:occupation_cutoff]` assertion;
- report retained counts, minimum retained occupation, `max |G' S R|`,
  `max |R' S R - I|`, allowed tolerance, and final merge condition;
- no Cr2 run.

No other committed fixture/test, driver workflow, artifact schema change,
solver/RHF, ECP, EGOI, Cr2 full Hamiltonian, Cr2 artifact, or Cr2 facade
support is approved.

### HP-RG-CUTOFF-FN-02 — production residual cutoff tightening

Status: approved.

Approved source files:

```text
src/cartesian_residual_gaussians/residual_basis.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

`residual_basis.jl` is the primary owner. The terminal residual file is
approved only for narrow compatibility keyword default plumbing if needed.

Approved target:

```text
residual_occupation_cutoff = 1.0e-6
identity_atol = 5.0e-8
```

Evidence: Cr2 residual spectra show the worst low-H1 modes are built from
marginal owner-local residual directions with occupations around
`1.27e-7` to `8.98e-7`. The current `5.0e-8` cutoff retains those directions;
`1.0e-6` drops `6` directions per owner in the cited audit. Broad residual
widths still matter, but broad width alone is not the first production rule
because one-center atoms can have broad candidates without the same bad
`H1_RR` sector.

Approved behavior:

- set the default owner-local residual occupation cutoff to `1.0e-6`;
- keep the final `R' S R` identity validation default at
  `identity_atol = 5.0e-8`;
- preserve explicit caller overrides where already supported;
- keep width/zeta filtering explicit and user-controlled;
- keep owner-local metric checks, negative-eigenvalue tolerances, final merge
  metric checks, and `G' S R` orthogonality checks active;
- keep owner grouping, final merge failure rules, MWG/IDA conventions,
  artifact schema/provenance/reader/manifest, driver workflow, and public API
  unchanged.

This ID supersedes only the `HP-RG-CUTOFF-FN-01` residual occupation default.
It does not change the residual-selection algorithm, identity tolerance,
negative-eigenvalue tolerances, or merge policy.

Forbidden:

- residual selection algorithm changes;
- kinetic or `H1_RR` spectral guards;
- global raw-candidate residual selection;
- global raw-column pivoted-Cholesky residual selection;
- owner grouping changes;
- negative-eigenvalue tolerance changes;
- final merge metric failure-rule changes or merge eigenvalue flooring;
- width/zeta filtering default changes or width filtering as conditioning
  repair;
- MWG/IDA, nuclear, raw-block, exact-operator, terminal-basis, WL/PQS route,
  shellification, Hamiltonian assembly, artifact schema, reader, manifest,
  driver, public API/export, solver/RHF, ECP, EGOI, Cr2 workflow, Cr2 artifact,
  or full HF changes;
- committed fixtures/tests except the exact existing H2 endpoint assertion
  update named under `HP-RG-CUTOFF-TEST-02`.

Failure rule: if the Cr2 residual-only audit does not cleanly remove the
identified low-occupation directions or if low-`H1_RR` ghost modes remain
after the cutoff change, make no further source changes in this lane. Report
the residual spectra and request separate kinetic/`H1_RR` spectral-guard
authority if needed.

### HP-RG-CUTOFF-TEST-02 — production residual cutoff validation

Status: approved.

Approved validation:

- Cr2 residual-only audit, not full HF or a Cr2 artifact/workflow run;
- Cr2 owner retained counts should drop from `68 + 68` to `62 + 62`;
- recompute and report residual spectra:
  - `min eig(K_RR)`;
  - `min eig(H1_RR)`;
  - low-mode candidate composition;
- if low-H1 ghost modes remain, stop and request separate kinetic/`H1_RR`
  spectral-guard authority;
- Be high-zeta residual construction still passes;
- H2 residual-GTO/MWG endpoint remains unchanged except for the approved
  default cutoff/provenance value;
- exactly update the existing committed H2 endpoint test
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` so both cutoff
  assertions expect `1.0e-6` instead of `5.0e-8`: the in-memory
  `residual.occupation_cutoff` assertion and the artifact/provenance
  `values[:occupation_cutoff]` assertion.

No other committed fixture/test, driver workflow, artifact schema change,
solver/RHF, ECP, EGOI, full HF, Cr2 full Hamiltonian, Cr2 artifact, or Cr2
facade support is approved.

### HP-RG-SPECTRAL-AUDIT-01 — residual-sector spectral audit

Status: approved measurement-only authority. This is not production source
authority.

Evidence after `HP-RG-CUTOFF-FN-02`: the tightened
`residual_occupation_cutoff = 1.0e-6` performs the intended first cleanup,
dropping Cr2 retained residuals from `68 + 68` to `62 + 62`, but residual-only
spectra still show a low two-owner residual mode:

```text
min eig(K_RR)  =  0.3700413519
min eig(H1_RR) = -7.1647854052
owner weights  = approximately 0.5 / 0.5
```

Approved behavior:

- measurement-only residual-sector audit;
- compute retained residual count by owner;
- compute low eigenvalues of `K_RR`;
- compute low eigenvalues of
  `H1_RR = K_RR + sum_A Z_A U_A_RR`;
- report owner weights for low or flagged eigenvectors;
- report residual occupation composition of low or flagged eigenvectors;
- compare Cr2 residual spectra against one-center atom residual baselines when
  available;
- classify whether low modes are dominated by the smallest retained
  occupations or by otherwise healthy retained modes.

Approved surfaces:

- ignored `tmp/work/*.jl` probes only;
- durable text/TSV output under `/Users/srw/dmrgtmp/...` or CR2 run
  directories.

Forbidden:

- production source changes;
- committed tests or fixtures;
- artifact schema/provenance/reader/manifest changes;
- driver changes;
- MWG/IDA changes;
- dense Vee, full HF, or solver workflow;
- automatic residual pruning;
- kinetic or `H1_RR` spectral-guard implementation;
- cutoff, tolerance, owner grouping, residual-selection, or merge-policy
  changes.

Validation for later audit:

- `git diff --check`;
- package load;
- residual-only audit for one-center Cr atom baseline if available;
- residual-only audit for the current Cr2 fixture;
- no full HF and no new Hamiltonian artifact.

Failure rule: if the audit cannot reconstruct `K_RR`/`H1_RR` cheaply enough
from existing construction seams, stop and report the exact missing reusable
seam. Do not add source instrumentation as part of this lane.

## Approved For Cartesian Gaussian Raw-Block Nuclear Owner

This section approves only the neutral uncharged by-center nuclear raw-block
slice recorded in `cartesian_gaussian_raw_blocks_nuclear.md`. It does not
approve a broad Gaussian raw-block framework.

### HP-CGRB-FILE-01 — neutral Cartesian Gaussian raw-block module files

Approved internal module and files:

```text
src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl
src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl
```

`src/GaussletBases.jl` may add only the internal include needed to load the
module, with include-order changes limited to immediate dependency/caller
needs. No public export is approved.

### HP-CGRB-FN-01 — exact uncharged Gaussian nuclear raw blocks

Approved internal kernel family:

```text
cartesian Gaussian nuclear raw blocks by center
```

The exact Julia names may follow local module style, but they must be neutral:
no `r3`, `residual`, `qw`, route-stage, report, or status vocabulary.

Approved numerical outputs:

- parent-supplement `G-A` matrices by nuclear center;
- supplement-supplement `A-A` matrices by nuclear center;
- uncharged unit attraction convention, `U_A = -1/r_A`.

Approved construction details:

- analytic one-dimensional nuclear factor construction;
- reuse across unique nuclear coordinates;
- upper-triangular `A-A` assembly and mirroring;
- function-local scratch/workspace reuse;
- term-first contraction over the Gaussian expansion.

The kernel must not apply physical nuclear charges, perform terminal
projection, transform into residual bases, create overlap/kinetic/moment
blocks, assemble Hamiltonians, or create persistent caches/bundles.

### HP-CGRB-FN-02 — nuclear one-dimensional axis-family reuse

Approved source owner:

```text
src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl
```

Approved optimization target: reorganize exact uncharged nuclear raw-block
construction around unique supplement one-dimensional axis families rather
than flattened 3D orbital labels.

Approved concepts:

- unique supplement axis-family inventory independent of flattened 3D orbital
  labels;
- integer map `orbital_axis_family[orbital, axis] -> family_id`;
- unique `G-A` table keys
  `(axis, supplement_axis_family, nuclear_axis_coordinate)`;
- unique `A-A` table keys
  `(axis, canonical(left_family, right_family), nuclear_axis_coordinate)`;
- transpose/orientation flags when canonical `A-A` family order is reversed;
- term-first filling of each required one-dimensional table at most once per
  Coulomb Gaussian term;
- reuse of those tables across all 3D orbitals and orbital pairs that reference
  the same axis families;
- coupled primitive-pair contraction
  `sum_pq c_p c_q Ix[p,q] Iy[p,q] Iz[p,q]`;
- function-local workspaces and integer lookup plans only.

Independent contraction of x/y/z axis tables into separate scalar contractions
is forbidden. The kernel must not introduce persistent caches, metadata,
status/report fields, route objects, payload structs, public API/export,
artifact changes, Residual Gaussian algorithm changes, Qiu-White route
semantic changes, overlap/kinetic/moment migration, Cr2 facade support, or Cr2
artifact workflow.

### HP-CGAI-FN-01 — optional Cartesian Gaussian axis helper

Approved source owner:

```text
src/cartesian_gaussian_axis_integrals.jl
```

Status: superseded as a performance endpoint. It remains an optional helper
surface only if needed by `HP-CGRB-FN-02`.

Optional internal helper concept:

```julia
_cartesian_gaussian_axis_integral_table!(
    destination,
    left_exponents,
    left_centers,
    left_powers,
    left_prefactors,
    right_exponents,
    right_centers,
    right_powers,
    right_prefactors,
    term;
    factor_exponent = 0.0,
    factor_center = 0.0,
)
```

The helper fills an already allocated destination matrix and must return the
same values as `_cartesian_gaussian_axis_integral_table(...)` without
allocating the result matrix. The existing scalar
`_cartesian_gaussian_axis_integral(...)` behavior is unchanged. The allocating
helper may delegate to the in-place helper if that preserves behavior cleanly.
The same owner may also add a specialized nonallocating nuclear-factor scalar
integral for the `:factor` term if needed by the `HP-CGRB-FN-02` family-reuse
kernel.

Allowed consumer surface:

```text
src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl
```

That file may consume the optional helper only inside the exact uncharged
by-center Gaussian nuclear raw-block construction under `HP-CGRB-FN-02`. Do
not treat result-matrix allocation removal as the accepted Cr2 optimization
target. No public API, export, persistent cache, raw-block payload,
metadata/status/report field, artifact change, route object, Residual Gaussian
algorithm change, Qiu-White route semantic change, overlap/kinetic/moment
migration, Cr2 facade, or Cr2 artifact workflow is approved by this ID.

### HP-CGRB-WIRE-01 — Residual Gaussian and Qiu-White rewiring

Approved caller rewiring surfaces:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/ordinary_qw_raw_blocks.jl
src/ordinary_qw_operator_assembly.jl
```

The implementation sequence is binding:

1. Extract current nuclear `G-A`/`A-A` behavior preserving conventions.
2. Rewire Residual Gaussian and Qiu-White callers to the neutral kernel.
3. Delete duplicate route-local nuclear loops once parity is established.
4. Optimize allocation inside the neutral owner only after extraction parity.

No Qiu-White route objects, Residual Gaussian selection logic, augmented
operator transforms, terminal projection, parent construction, artifact
workflow, report/status/payload fields, public API, or Cr2 facade/artifact
workflow may be added under this ID.

### HP-CGRB-TEST-01 — nuclear extraction validation

Approved validation:

- existing H2 Residual Gaussian endpoint unchanged;
- Be2 Residual Gaussian endpoint unchanged, as ignored validation if needed;
- Cr2 q4 exact nuclear `G-A`/`A-A` blocks match the current implementation at
  roundoff, as ignored measurement only;
- one small Qiu-White nuclear parity fixture.

Approved committed standalone parity file, if needed:

```text
test/nested/cartesian_gaussian_raw_blocks_nuclear_runtests.jl
```

This file must not be added to `test/runtests.jl` without a later amendment.
It must not validate Cr2 workflow, artifacts, public API, route statuses,
report mirrors, payload fields, or Residual Gaussian internals.

## Approved For Cartesian Gaussian Raw-Block Non-Nuclear Owner

This section approves only the non-nuclear raw-block slice recorded in
`cartesian_gaussian_raw_blocks_non_nuclear.md`. It extends the existing neutral
Cartesian Gaussian raw-block owner under new `HP-CGRB-NN-*` IDs. It must not be
implemented under `HP-CGRB-FN-02`.

### HP-CGRB-NN-FILE-01 — non-nuclear raw-block file

Approved owner file:

```text
src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl
```

Allowed module plumbing:

```text
src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl
```

Only the include needed to load `non_nuclear_blocks.jl` is approved there.
Root include changes in `src/GaussletBases.jl` are not expected and are not
approved unless a later amendment identifies a real include-order blocker.

### HP-CGRB-NN-FN-01 — exact non-nuclear Gaussian raw blocks

Approved internal kernel family:

```text
cartesian Gaussian non-nuclear raw blocks
```

The exact Julia names may follow local module style, but they must be neutral:
no `r3`, `residual`, `qw`, route-stage, report, or status vocabulary.

Approved numerical outputs:

- parent-supplement `G-A` overlap, kinetic, coordinate moments, and second
  moments;
- supplement-supplement `A-A` overlap, kinetic, coordinate moments, and second
  moments.

Approved construction details:

- analytic one-dimensional table construction;
- unique supplement axis-family reuse;
- canonical `A-A` family-pair table keys and orientation handling;
- upper-triangular `A-A` assembly and mirroring;
- function-local scratch/workspace reuse;
- coupled product-axis contraction preserving existing Qiu-White values;
- reuse of once-built overlap `G-A` for residual setup mixed overlap
  `X = G' S A` and exact augmented-operator assembly when both are built in
  the same local construction call.

The kernel may return a compact fixed-field internal result containing only the
approved raw matrices. It must not be a status object, route stage, report
payload, metadata carrier, persistent cache, broad provider bundle, or artifact
data.

### HP-CGRB-NN-WIRE-01 — Residual Gaussian and Qiu-White non-nuclear rewiring

Approved caller rewiring surfaces:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/ordinary_qw_raw_blocks.jl
src/ordinary_qw_operator_assembly.jl
```

The implementation sequence is binding:

1. Extract current Qiu-White non-nuclear `G-A`/`A-A` behavior preserving
   conventions.
2. Rewire Residual Gaussian exact-operator construction and residual mixed
   overlap setup to consume the neutral output.
3. Rewire Qiu-White consumers to the neutral output.
4. Delete duplicate route-local non-nuclear loops once parity is established.
5. Optimize allocation inside the neutral owner only after extraction parity.

No nuclear raw-block changes, final-basis `G-G` product-matrix optimization,
terminal projection, Residual Gaussian algorithm changes, augmented-operator
transform changes, Qiu-White semantic changes, Qiu-White route objects, parent
construction, persistent cache, report/status/payload fields, public API,
artifact workflow, or Cr2 facade/artifact workflow may be added under this ID.

### HP-CGRB-NN-TEST-01 — non-nuclear extraction validation

Approved validation:

- existing H2 Residual Gaussian endpoint unchanged;
- Be2 Residual Gaussian endpoint unchanged, as ignored validation if needed;
- Cr2 q4 non-nuclear `G-A`/`A-A` overlap, kinetic, coordinate moment, and
  second-moment blocks match the current implementation at roundoff, as
  ignored measurement only;
- residual setup mixed overlap `X` matches the current construction at
  roundoff;
- one small Qiu-White non-nuclear parity fixture.

Approved committed standalone parity file, if needed:

```text
test/nested/cartesian_gaussian_raw_blocks_non_nuclear_runtests.jl
```

This file must not be added to `test/runtests.jl` without a later amendment.
It must not validate Cr2 workflow, artifacts, public API, route statuses,
report mirrors, payload fields, or Residual Gaussian internals.

## Approved For R3 Terminal G-G Product Matrices

This section approves only the terminal final-basis `G-G` product-matrix
optimization recorded in `r3_terminal_gg_product_matrices.md`. It is owned by
`CartesianFinalBasisRealization`, not by `CartesianGaussianRawBlocks`.

### HP-R3GG-FN-01 — R3/RG terminal G-G product-matrix optimization

Approved owner:

```text
Owner module: CartesianFinalBasisRealization
```

Approved source files:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/cartesian_final_basis_realization/pqs_terminal_one_body.jl
```

The first implementation should prefer editing only
`pqs_terminal_residual_gto.jl`. Edits to `pqs_terminal_one_body.jl` are
approved only for a small internal terminal-product workspace or multi-product
helper needed to reuse function-local buffers across consecutive product
assemblies.

Approved product matrices:

- kinetic `K_GG`;
- coordinate moments `x_GG`, `y_GG`, `z_GG`;
- second moments `x2_GG`, `y2_GG`, `z2_GG`.

Approved implementation shapes:

- accumulate the three kinetic-axis product contributions into one destination;
- reuse an already constructed base Hamiltonian kinetic `G-G` block in the
  same-construction path when available and validated equal;
- build coordinate and second-moment `G-G` products one axis at a time and
  transform immediately;
- share function-local scratch/workspace across consecutive terminal product
  assemblies;
- delete or simplify `_r3a_product_matrix(...)` when replaced and no live
  caller remains.

This ID does not approve `G-A`/`A-A` raw-block changes, nuclear raw-block
changes, unit-nuclear `U_A` Gaussian-sum changes, terminal basis realization
changes, residual Gaussian algorithm changes, Qiu-White semantic changes,
IDA/MWG changes, parent construction, persistent caches, metadata,
report/status/payload fields, public API/export, artifact changes, Cr2 facade
support, or Cr2 artifact workflow.

Line budget: at most `100` added `src` lines total. If implementation needs a
broad product-operator framework, persistent workspace/cache object, files
outside the approved source files, or a public/internal payload, stop and
request a new docs-only amendment.

### HP-R3GG-TEST-01 — terminal G-G product validation

Approved validation:

- existing H2 Residual Gaussian endpoint unchanged;
- Be2 Residual Gaussian usability/performance measurement unchanged except for
  allowed timing/allocation improvement;
- Cr2 q4 `K_GG`, coordinate moment `G-G`, and second-moment `G-G` products
  match the current construction at roundoff as ignored validation;
- augmented exact operators remain finite and symmetric;
- base `G-G` block equality checks in the existing H2 endpoint still pass;
- Cr2 q4 exact-operator allocation is remeasured after parity.

No new committed test file is approved by this ID.

## Approved For R3 Unit-Nuclear U_GG Gaussian Sum

This section approves only the terminal final-basis unit-nuclear `U_GG`
Gaussian-sum optimization recorded in
`r3_unit_nuclear_ugg_gaussian_sum.md`. It is owned by
`CartesianFinalBasisRealization`, not by `CartesianGaussianRawBlocks`.

### HP-R3UN-FN-01 — R3/RG unit-nuclear U_GG Gaussian-sum optimization

Approved owner:

```text
Owner module: CartesianFinalBasisRealization
```

Approved source files:

```text
src/cartesian_final_basis_realization/pqs_terminal_one_body.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

The first implementation should prefer `pqs_terminal_one_body.jl`. Edits to
`pqs_terminal_residual_gto.jl` are approved only for narrow R3 exact-operator
caller wiring needed to use function-local scratch or the optimized helper.

Approved target functions:

```text
_accumulate_terminal_gaussian_sum!
_terminal_gaussian_sum_action
```

Approved implementation shapes:

- reuse function-local scratch/workspace across Gaussian-sum terms and center
  calls;
- accumulate terminal Gaussian-sum contributions in-place into the caller's
  destination;
- reduce avoidable allocation in factor lookup and terminal Gaussian-sum action
  construction;
- add small internal scratch arguments or file-local helpers only if they remain
  inside `CartesianFinalBasisRealization` and create no persistent state;
- simplify or delete allocation-heavy helper code inside the targeted
  Gaussian-sum path after parity.

This ID does not approve neutral raw-block changes, terminal kinetic/moment
`G-G` product changes, residual Gaussian selection/orientation/transform
changes, MWG/IDA changes, Qiu-White semantic changes, route/stage setup
cleanup, raw-block setup cleanup, parent construction, terminal basis
realization changes, persistent caches/workspaces, broad Gaussian-sum
frameworks, metadata/report/status/payload fields, artifacts, public
API/export, Cr2 facade support, or Cr2 artifact workflow.

Line budget: at most `100` added `src` lines total. If implementation needs a
broad Gaussian-sum framework, persistent cache/workspace object, files outside
the approved source files, or source edits outside the terminal unit-nuclear
`U_GG` path, stop and request a new docs-only amendment.

### HP-R3UN-TEST-01 — unit-nuclear U_GG validation

Approved validation:

- existing H2 Residual Gaussian endpoint unchanged;
- Be2 Residual Gaussian facade/readback unchanged except for allowed
  timing/allocation improvement;
- Cr2 q4 exact-operator audit reports before/after unit-nuclear `U_GG`
  allocation and total wrapper allocation;
- Cr2 q4 unit-nuclear `U_GG` block replay parity and final exact augmented
  operator parity at roundoff;
- exact operators remain finite and symmetric.

No new committed test file is approved by this ID.

## Approved For R3 Same-Construction Base K/U Reuse

This section approves only narrow reuse of already-built same-construction
base final-basis kinetic and unit-nuclear blocks in supplemented residual-GTO
/ MWG exact augmented operators. It is an orchestration reuse lane, not a
terminal product, Gaussian-sum, raw-block, residual-basis, or interaction
algorithm lane.

Evidence recorded before approval: a replay that reused same-construction base
`K_GG` and unit `U_GG[A]` blocks had exact operator delta `0.0` and reduced the
exact augmented-operator replay to `0.8620s / 1237.136 MiB`.

### HP-R3BASE-FN-01 — same-construction base K/U reuse

Approved owner:

```text
Owner module: CartesianFinalBasisRealization plus narrow caller wiring
```

Approved source files:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- `pqs_terminal_residual_gto_augmented_products(...)`, or its approved caller
  wrapper, may accept a trusted same-construction base kinetic matrix and use
  it as the `G-G` kinetic block for `transform_augmented_operator`;
- `pqs_terminal_residual_gto_augmented_unit_nuclear(...)`, or its approved
  caller wrapper, may accept trusted same-construction unit nuclear
  `U_GG[A]` matrices and use them as the `G-G` unit blocks for
  `transform_augmented_operator`;
- `cartesian_residual_gto_mwg_hamiltonian(...)` and staged helpers in
  `src/cartesian_base_hamiltonian.jl` may pass `base_ham.kinetic` and
  `base_ham.nuclear_attraction_unit_by_center` into the augmented operator
  construction;
- current behavior must be preserved when trusted base blocks are not supplied.

Trust condition:

- the base Hamiltonian, terminal basis realization, parent bundles, residual
  basis, and supplement must come from the same
  `cartesian_base_working_basis(...)` construction path;
- implementation must validate matrix dimensions and center count before
  reuse;
- no provenance payload, metadata proof, report field, status object, or
  persistent cache is required or approved for this trust check.

This ID does not approve public API/export changes, canonical-driver changes,
raw-block changes, residual selection/orientation/transform changes, MWG/IDA
convention changes, terminal product or Gaussian-sum kernel rewrites,
persistent cache/workspace objects, metadata/status/report/artifact schema
fields, route/stage setup cleanup, committed tests, Cr2 workflow, or source
files outside the two approved files.

Line budget: target under `100` added `src` lines. If trusted
same-construction provenance cannot be guaranteed by local call shape plus
dimension/center validation, or if implementation needs public payloads,
metadata, or stage objects, make no source commit and report the blocker.

### HP-R3BASE-TEST-01 — base K/U reuse validation

Approved validation:

- `git diff --check`;
- package load;
- existing H2 R3 endpoint unchanged;
- Be2 supplemented facade/readback unchanged except allowed
  timing/allocation improvement;
- Cr2 exact-operator attribution audit or focused ignored replay showing base
  `K_GG` / unit `U_GG[A]` reuse parity and allocation effect;
- final exact operators finite and symmetric.

No new committed test file, Cr2 artifact, Cr2 workflow, public API/export,
driver workflow change, metadata/status/report field, or artifact schema
change is approved by this ID.

### HP-R3BASE-DRV-WIRE-01 — canonical driver K/U reuse call-site wiring

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

Approved behavior:

- in supplemented mode only, pass `base_ham.kinetic` into
  `cartesian_residual_gto_augmented_products(...)` as `base_kinetic`;
- in supplemented mode only, pass
  `base_ham.nuclear_attraction_unit_by_center` into
  `cartesian_residual_gto_augmented_unit_nuclear(...)` as
  `base_unit_nuclear`;
- keep the current public inputs, hooks, timing labels, visible stage sequence,
  artifact schema, and driver contract unchanged.

This ID is only call-site wiring so the canonical driver uses the already
approved same-construction base K/U reuse path. It does not approve source or
kernel changes, diagnostics, new hooks, new timing labels, public API/export
changes, artifact changes, tests/fixtures, Cr2 workflow, or edits outside
`bin/cartesian_ham_builder.jl`.

Failure rule: if the driver call-site update needs any visible driver contract
change, make no source commit and report the blocker.

### HP-R3BASE-DRV-TEST-01 — driver K/U reuse validation

Approved validation:

- `git diff --check`;
- package load;
- H2 supplemented driver artifact/readback;
- Be2 supplemented driver artifact/readback if practical;
- no Cr2 run.

No new committed test file, fixture, diagnostic, hook, timing label, public
input, artifact schema, or Cr2 workflow is approved by this ID.

## Approved For Canonical Cartesian Driver Usability

This section approves only the compact artifact-producing canonical driver
workflow recorded in `cartesian_driver_usability_workflow.md`. It is workflow
authority over approved producer surfaces, not algorithm, kernel, solver,
artifact-schema, or diagnostic authority.

### HP-DRV-FILE-01 — canonical driver file

Approved file:

```text
bin/cartesian_ham_builder.jl
```

No other `bin`, `tools`, `src`, `test`, or committed driver-input fixture file
is approved by this ID.

### HP-DRV-FN-01 — compact functional driver workflow

Approved invocation shape:

```text
julia --project=. bin/cartesian_ham_builder.jl [input.jl] [key=value ...]
```

Approved behavior:

- visible editable defaults near the top of the driver;
- optional trusted local Julia input file for project-specific defaults;
- later command-line `key=value` overrides;
- visible public `system`, `basis`, and optional `supplement` contract
  construction before calling an approved facade;
- compact normalized run summary;
- coarse user-facing phase timing;
- visible physics-level construction stages through the staged producer surface;
- base or supported supplemented Hamiltonian construction through approved
  producer surfaces;
- artifact write;
- optional readback check.

Approved configuration concepts are `basisname`, `system`, base `basis`,
optional `supplement`, `nesting`, `hamfile`, `padding`, `check_file`,
`print_contract`, `print_timing`, and `expected_dimension`.

Compact summary printing and artifact readback checks remain allowed workflow
behavior, but they are not open-ended hooks and must not introduce route,
diagnostic, artifact-schema, or solver controls.

`basisname = nothing` selects base mode. `basisname !== nothing` selects a
supported supplemented mode and is the visible supplement basis label. The
original driver workflow lane covered supplemented diatomics only;
`HP-COMP-SUPPATOM-*` separately approves relaxing the old `Natom == 1`
rejection.

`padding` is a public physical box-padding control. For one-center atoms it
maps to the base facade `radius`. For z-axis diatomics it maps to the existing
facade extents as padding around the two nuclei; under the current origin-based
z-axis contract this means `xmax_parallel = max(abs(z_i)) + padding` and
`xmax_transverse = padding`.

Approved hooks are only `check_file`, `print_contract`, `print_timing`, and
`expected_dimension`. They may support human expert review and
Codex-controlled artifact checks. They must not expose route internals,
stop-after stages, raw-block switches, allocation probes, artifact schema
dumps, solver controls, Cr2-specific workflow, or private helper calls.
`check_file` may contain compact public contract facts, artifact path, final
dimension, expected-dimension result, readback deltas, and coarse timing only.

This ID does not approve private route-stage controls, stop-after internals,
ladder probes, stage markers, fixture hacks, diagnostic knobs, underscored
package helper calls from the driver, raw-block provider switches,
report/status/payload dumps, metadata field clouds, allocation probes,
benchmark harness behavior, solver/RHF/ECP/EGOI/HamV6 workflow, public
API/export changes, artifact schema changes or dumps, committed test files,
committed driver-input fixtures, unsupported atom/supplement combinations, or
Cr2-specific workflow support. Generic explicit homonuclear z-axis Cr2 stress through
`HP-R3U-ZDI-WIRE-01` is separate ignored/user-run validation authority, not
driver-owned Cr2 support.

Line budget: at most `150` added `bin` lines. If implementation needs a parser
framework, source files outside the approved driver and staged producer
surfaces, committed input fixtures, route-stage diagnostics,
status/report/payload expansion, artifact schema changes, or Cr2-specific
workflow support, stop and request a new docs-only amendment.

### HP-DRV-NEST-FN-01 — construction-family driver input

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

Approved behavior:

- add a visible public driver input `nesting`;
- accepted values are `:pqs` and `:wl`;
- default is `nesting = :pqs`;
- `nesting = :pqs` means the PQS source-box construction family;
- `nesting = :wl` means the White-Lindsey low-order construction family;
- include `nesting` in public contract construction, optional
  `print_contract`, and optional `check_file` output as a public contract fact.

This is a first-class construction-family choice, not a diagnostic route
switch. The driver must not expose internal route-family names, route
skeletons, retained-rule plans, raw-block switches, stop-after controls,
diagnostic knobs, old route-stage labels, private helper calls, allocation
probes, or route reports.

This ID does not approve public API/export changes, artifact schema changes,
stage-label changes, solver/ECP workflow, Cr2-specific behavior, broad driver
diagnostics, committed fixtures/tests, or source files outside the canonical
driver.

### HP-DRV-NEST-WIRE-01 — construction-family route mapping

Approved source files:

```text
bin/cartesian_ham_builder.jl
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- map public `nesting = :pqs` to the existing internal `:pqs_source_box`
  route family;
- map public `nesting = :wl` to the existing internal
  `:white_lindsey_low_order` route family;
- keep route skeletons, retained rules, raw-block switches, stop-after
  controls, diagnostics, and internal route-stage vocabulary hidden;
- preserve the existing public stage labels, Hamiltonian object, matrix keys,
  artifact schema, driver hooks, and solver-free workflow;
- reject unsupported combinations with clear `ArgumentError`s.

Supplemented `nesting = :wl` is governed by `HP-COMP-SUPPWL-*` for the
supported homonuclear z-axis diatomic composition cell. Unsupported geometry
or supplement combinations must still reject clearly rather than adding broad
White-Lindsey route behavior.

This ID does not approve new route algorithms, route-skeleton construction
changes, White-Lindsey materialization deletion, terminal lowering policy
changes, shellification behavior changes, numerical kernel changes, raw-block
changes, Residual Gaussian/MWG/IDA changes, artifact/provenance schema changes,
public API/export changes, committed tests, Cr2-specific workflow, or source
files outside the two approved files.

Line budget: at most `80` added source/bin lines, with net simplification
preferred where old hidden assumptions can be removed.

Failure rule: if `nesting = :wl` cannot produce a small base artifact/readback
through the existing White-Lindsey low-order route without broader route or
materialization work, make no source commit and report the exact blocker.

### HP-DRV-NEST-TEST-01 — construction-family validation

Approved validation:

- `git diff --check`;
- package load;
- current default `nesting = :pqs` base driver or facade artifact/readback;
- current default `nesting = :pqs` supplemented H2 driver/facade path if
  supplemented-mode input plumbing is touched;
- one small `nesting = :wl` base artifact/readback using a currently supported
  base geometry;
- explicit negative check or ignored smoke showing unsupported supplemented
  `nesting = :wl` combinations fail clearly outside the supported
  `HP-COMP-SUPPWL-*` cell;
- no Cr2 run.

No new committed test file, committed input fixture, artifact schema
validation, solver run, Cr2-specific driver run, or broad White-Lindsey
workflow validation is approved.

### HP-DRV-STAGE-FN-01 — visible physics-stage producer surface

Approved source file:

```text
src/cartesian_base_hamiltonian.jl
src/pqs_source_box_low_order_materialization.jl
src/cartesian_final_basis_realization/pqs_terminal_one_body.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

Approved purpose: expose a small non-exported, non-underscored,
driver-facing staged producer surface so the canonical driver can execute and
time visible physics-level construction stages without calling private
underscored helpers.

`src/cartesian_base_hamiltonian.jl` remains the primary driver-facing owner.
The lower-level files listed above are approved only for behavior-preserving
operator-class stage factoring in their existing domains; they are not
approved for new algorithms, raw-block changes, or diagnostics.

Approved visible stages are:

- construct public `system`, `basis`, and optional `supplement`;
- build the base working basis / terminal realization;
- build base product/moment operators;
- build base unit-nuclear attraction operators;
- build base electron-electron / IDA interaction;
- assemble the base Hamiltonian;
- load or build the Gaussian supplement basis when `basisname !== nothing`;
- build residual Gaussian augmentation;
- build augmented product/moment operators;
- build augmented unit-nuclear attraction operators;
- build augmented electron-electron / residual-MWG interaction;
- assemble the supplemented Hamiltonian;
- write and check the artifact.

Approved physical operator classes are:

- product/moment operators: kinetic `K`, Cartesian coordinate moments
  `x`/`y`/`z`, and second moments `x^2`/`y^2`/`z^2` where present;
- unit-nuclear attraction: uncharged by-center `U_A` / `Vnuc` matrices before
  applying physical nuclear charges;
- electron-electron interaction: base localized IDA `Vee` and supplemented
  residual-MWG/IDA `Vee`.

The first and last stages remain driver/writer responsibilities. This ID
approves source factoring needed for the base working-basis/terminal
realization, base product/moment, base unit-nuclear, base IDA, Gaussian
supplement, residual augmentation, augmented product/moment, augmented
unit-nuclear, residual-MWG/IDA, and Hamiltonian assembly stages.

The staged surface may factor the existing `cartesian_base_hamiltonian(...)`
and `cartesian_residual_gto_mwg_hamiltonian(...)` bodies so that those facades
can remain wrappers over the same implementation. It may return existing
domain objects and small fixed-key ephemeral stage products required by the
next approved stage.

The staged surface must be a set of separate named construction-stage
functions. It must not be a single opaque replacement wrapper that hides the
same construction sequence under a new name. The canonical driver must be able
to bind visible local variables for the base realization, base products,
base unit-nuclear operators, base `Vee`, base Hamiltonian, supplement basis,
residual augmentation, augmented products, augmented unit-nuclear operators,
augmented `Vee`, and final Hamiltonian assembly.

This ID does not approve public exports, public API redesign, route-stage
objects, reports, status/result payloads, metadata field clouds, runtime-keyed
field groups, persistent caches, raw-block switches, allocation probes,
per-kernel timing frameworks, solver/ECP workflow, artifact schema changes, or
source files outside the four approved paths listed above.

Line budget: at most `200` added `src` lines across the approved source files.
If the staged surface requires a new module, source files outside the four
approved paths, a broad payload object, committed tests, new artifact keys,
raw-block changes, or kernel rewrites, stop and request a new docs-only
amendment.

### HP-DRV-STAGE-WIRE-01 — canonical driver staged wiring

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

The canonical driver may call the `HP-DRV-STAGE-FN-01` staged producer surface
as separate top-level stage calls and assign local variables using the approved
physics-stage names. It may print or record coarse user-facing timings for
those stages. Replacing the current facade call with one all-in-one staged
wrapper call is not approved for the canonical driver.

Driver timing should expose the three physical operator classes separately:
product/moment, unit-nuclear, and electron-electron. These timings are
user-facing stage timings only. They must not become allocation probes,
raw-block timing controls, per-kernel instrumentation, or diagnostic stop
points.

This ID does not approve calls from the driver to underscored package helpers,
old route stages such as `cartesian_parent`, `cartesian_shells`,
`cartesian_units`, `cartesian_pair_terms`, or `cartesian_assembly`, raw-block
provider switches, stop-after controls, route diagnostics, allocation probes,
artifact schema dumps, solver controls, Cr2-specific workflow, or new
committed fixtures/tests.

### HP-DRV-STAGE-TEST-01 — staged driver validation

Approved validation:

- package load;
- H atom or H2 base driver run with visible base-stage timing/summary;
- H2 supplemented driver run with visible supplement/residual/operator/
  Hamiltonian stage timing/summary;
- artifact write/readback still passes for those runs;
- `expected_dimension`, `print_contract`, and `check_file` behavior still
  uses only public contract and coarse stage facts.

No committed test file, committed input fixture, Cr2-specific driver run,
solver run, or diagnostic harness is approved by this ID.

### HP-DRV-INV-FN-01 — canonical driver terminal-region inventory

Status: approved.

Approved source files:

```text
bin/cartesian_ham_builder.jl
src/cartesian_base_hamiltonian.jl
```

Optional only if a compact accessor is directly required:

```text
src/pqs_source_box_route_driver_helpers.jl
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
src/cartesian_final_basis_realization/terminal_face_product_blocks.jl
```

Problem:

Large accidental identity sectors, such as the Cr2 z-end slab blowup, are
hard to notice from the current canonical driver output. Normal driver users
should see a compact basis-region inventory without running ignored probes and
without receiving a route-debug dump.

Approved behavior:

- print a compact terminal-region / shellification inventory as part of
  canonical driver output;
- include it for base construction;
- for supplemented construction, report at least the base terminal inventory
  plus final supplemented dimension;
- keep the output bounded and human-facing;
- preserve the existing driver stage sequence and public inputs;
- preserve artifact schema, matrix keys, reader behavior, and readback checks.

Minimum useful columns:

- region key/index or compact label;
- region kind;
- lowering kind or final realization kind;
- shell index or explicit unavailable status;
- support row count;
- retained/final column count;
- compression ratio;
- identity versus compact/product realization;
- index ranges for each axis, `x = i:j`, `y = k:l`, `z = m:n`;
- physical coordinate ranges for each axis, `x`, `y`, and `z`;
- slab normal axis, side, thickness, and stack index/count when applicable.

The summary should also print total base final dimension, supplemented final
dimension when applicable, and a clear count or visible rows showing any
direct identity slab sectors if they exist.

The geometry columns are part of the canonical inventory contract. They are
needed to catch shellification errors where every region is compact but the
z-axis slab stack is emitted only after the final shared shell instead of being
interleaved with the angular-balanced shell steps. Physical `x`/`y` ranges are
required, not only `z`, because the angular-balance rule compares the
transverse physical scale against the bond-axis margin.

This ID does not approve route skeleton exposure, source-mode inventories,
pair inventories, raw-block details, all-row listings, full metadata dumps,
recursive route-stage dumps, new driver inputs, flags, stop-after controls,
route switches, diagnostic switches, solver settings, broad status/report
payloads, artifact schema changes, reader changes, public API/export changes,
numerical construction changes, shellification changes, terminal lowering,
retained-unit changes, transform-contract changes, terminal-realization
changes, Residual Gaussian, MWG, IDA, Hamiltonian assembly, raw-block changes,
Cr2-specific workflow, committed Cr2 fixtures, or committed tests.

Line budget: target at most `80` added `src`/`bin` lines. This should be
formatting plus compact accessor work, not a reporting subsystem.

Failure rule: if the driver cannot print this from existing compact stage or
final-basis summaries without adding a broad payload, artifact fields, or a
route-report framework, make no source commit and report the missing summary
seam.

### HP-DRV-INV-TEST-01 — terminal-region inventory validation

Status: approved.

Approved validation:

- `git diff --check`;
- package load;
- bounded H2 or Be2 driver run showing the summary for `nesting = :pqs`;
- bounded H2 or Be2 driver run showing the summary for `nesting = :wl`;
- supplemented smoke if the printed summary touches supplemented-stage
  objects;
- artifact/readback deltas unchanged;
- output remains bounded and excludes source modes, pair inventories,
  raw-block details, all-row listings, and full metadata;
- output includes shell index or explicit unavailable status, index ranges for
  all axes, physical coordinate ranges for all axes, and slab stack facts when
  applicable;
- no Cr2 run required; optional user-side Cr2 run only.

No committed test file, committed driver-input fixture, Cr2-specific driver
run, artifact schema validation, solver run, or diagnostic harness is approved
by this ID.

### HP-DRV-TEST-01 — driver workflow validation

Approved validation:

- package load;
- public contract construction and optional `print_contract`/`check_file`
  output for at least one base run when driver construction code changes;
- H2 base driver run writes a `CartesianIDAHamiltonian` artifact and optional
  readback passes;
- H2 supplemented driver run writes a supplemented `CartesianIDAHamiltonian`
  artifact with approved compact `supplement_provenance/` and optional readback
  passes;
- optional ignored Be2 usability run if the implementation touches
  supplemented mode.

Validation input files, if needed, must be ignored `tmp/work` files. No
committed test file, committed driver-input fixture, Cr2-specific driver run,
or solver run is approved by this ID.

## Approved For Canonical Driver Atom Workflow

This section approves only the base atom workflow recorded in
`cartesian_driver_atom_workflow.md`. It is driver authority over the existing
base facade, not new atom physics, Residual Gaussian, artifact-schema, or
solver authority.

### HP-DRV-ATOM-FN-01 — explicit base atom driver workflow

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

Approved behavior:

- accept explicit one-center atom input in `mode = :base`;
- require `atom_symbols`, `nuclear_charges`, `atom_locations`, `nup`, and
  `ndn`;
- require exactly one center at `(0.0, 0.0, 0.0)`;
- require finite positive explicit nuclear charge and neutral all-electron
  count `nup + ndn == round(Int, only(nuclear_charges))`;
- pass explicit one-center base `basis` fields, including `core_spacing`, to
  the existing base facade;
- allow visible, easily edited driver/project defaults such as
  `core_spacing = 0.3` and template `padding`, while treating them as explicit
  resolved input values that may be overridden for quick tests;
- use clear `ArgumentError`s for unsupported atom workflow inputs where
  practical.

Current driver validation remains origin-centered H. This driver ID does not
approve changing `src/cartesian_base_hamiltonian.jl`; producer-side
one-center atom support is governed separately by `HP-R1-ATOM-*`.

### HP-DRV-ATOM-WIRE-01 — driver atom-to-base-facade wiring

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

Approved behavior:

- the canonical driver may call
  `cartesian_base_hamiltonian(system; basis, hamfile)` for one-center base
  atom construction;
- artifact write/readback uses the existing base facade and existing
  `producer_provenance/` schema;
- no package-internal route-stage helper, terminal basis object, raw-block
  provider, report/status/payload object, or new artifact field is approved.

This base-atom driver wiring ID did not itself approve supplemented atom
Hamiltonians. Supported origin-centered one-center supplemented atom wiring is
now governed by `HP-COMP-SUPPATOM-*`. If a requested atom is outside the
existing base or supported supplemented facade scope, the implementation must
stop at a clear unsupported-input error rather than adding broader atom
construction.

Line budget for `HP-DRV-ATOM-FN-01` plus `HP-DRV-ATOM-WIRE-01`: at most `80`
added `bin` lines, with no committed test, tool, or input-fixture file.

### HP-DRV-ATOM-TEST-01 — base atom driver validation

Approved validation:

- package load;
- origin-centered H base driver artifact write/readback with explicit system
  and one-center basis fields;
- optional ignored negative checks for non-origin atom input, nonneutral
  electron count, mismatched temporary `d` if accepted, or unsupported atom
  input.

No supplemented atom endpoint was approved by this base-atom driver test ID;
supported supplemented atom validation is owned by
`HP-COMP-SUPPATOM-TEST-01`. No translated-atom gate, committed atom fixture,
new committed test file, solver run, artifact-schema validation, or broader
base-atom validation is approved by this ID.

### HP-DRV-ATOM-CLEAN-01 — remove hidden atom `d` driver residue

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

Approved behavior:

- remove the hidden one-center atom basis field `d = vars[:core_spacing]`;
- keep the visible driver atom basis in terms of `ns`, `core_spacing`,
  `radius`, and existing optional public fields only;
- keep public inputs, defaults, overrides, hooks, timing labels, visible stage
  sequence, artifact schema, and driver contract unchanged.

This ID exists because the producer no longer requires public `d` for
one-center atoms. It does not approve source/kernel changes, diagnostics, new
hooks, new timing labels, public input changes, committed tests/fixtures, Cr2
workflow, old `:white_lindsey_low_order` retirement, test/tool route-input
cleanup, or edits outside `bin/cartesian_ham_builder.jl`.

Failure rule: if removing the hidden `d` field requires any visible driver
contract change or producer/source change, make no source commit and report the
blocker.

## Approved For Complete-Core-Shell RHF Retirement

This section approves only the deletion lane recorded in
`complete_core_shell_rhf_retirement.md`. The old complete-core-shell RHF stack
is stale route-era workflow machinery, not a live Cartesian Hamiltonian
producer path.

### HP-RETIRE-CCS-RHF-FN-01 — remove stale RHF payload stack

Approved source files:

```text
src/GaussletBases.jl
src/pqs_multilayer_complete_core_shell_rhf.jl
```

Approved behavior:

- remove the `pqs_multilayer_complete_core_shell_rhf.jl` include from
  `src/GaussletBases.jl`;
- delete `src/pqs_multilayer_complete_core_shell_rhf.jl`;
- remove only docs/index references that describe this RHF stack as active
  current code, if encountered during the deletion pass;
- add no replacements, adapters, compatibility wrappers, status objects,
  payload objects, reports, or tests.

Evidence: focused search found no live `src`, `bin`, `test`, or `tool` caller
outside the file itself and the root include. The file carries old payload and
blocked-status vocabulary such as
`pqs_multilayer_complete_core_shell_rhf_input_contract`,
`pqs_multilayer_complete_core_shell_rhf_scf_payload`, and
`pqs_multilayer_complete_core_shell_rhf_one_step_payload`. The current
CR2-facing producer path is `bin/cartesian_ham_builder.jl`, staged producer
functions, and `CartesianIDAHamiltonian` artifacts.

This ID does not approve canonical driver changes, source changes outside the
approved file/include except minimal stale active-reference cleanup, changes to
`pqs_multilayer_complete_core_shell_h1.jl`,
`pqs_complete_core_shell_final_basis.jl`, or
`pqs_source_box_low_order_materialization.jl`, ordinary/Qiu-White donor-kernel
changes, artifact schema/provenance/reader changes, route/shellification/
terminal-lowering/raw-block/RG/MWG/IDA changes, Hamiltonian assembly changes,
committed tests/fixtures, or Cr2 workflow.

Failure rule: if any live `src`, `bin`, `test`, or `tool` caller depends on the
RHF stack, make no source commit and report the exact caller. Do not preserve
the path through an adapter.

### HP-RETIRE-CCS-RHF-TEST-01 — retirement validation

Approved validation:

- `git diff --check`;
- package load;
- focused `rg` showing no remaining live references to
  `pqs_multilayer_complete_core_shell_rhf`;
- canonical small base artifact/readback smoke;
- canonical small supplemented artifact/readback smoke;
- H2 Residual Gaussian endpoint remains unchanged;
- no Cr2 run.

No committed test, fixture, replacement path, adapter, status/report/payload
object, artifact-schema validation, or Cr2 workflow is approved.

## Approved For Route-Driver Materialization Workflow Retirement

This section approves only the retirement/quarantine lane recorded in
`route_driver_materialization_retirement.md`. The old route-driver
materialization/report/save wrapper workflow is not the canonical Cartesian
producer path; the current path is the staged human-facing driver plus
`CartesianIDAHamiltonian` artifacts.

### HP-RETIRE-DRV-MAT-FN-01 — remove old materialization/report/save wrappers

Approved source files:

```text
src/pqs_source_box_route_driver_helpers.jl
src/pqs_source_box_low_order_materialization.jl
src/pqs_source_box_route_driver_reporting.jl
src/GaussletBases.jl
```

`src/GaussletBases.jl` is approved only if the reporting include becomes unused
after the wrapper/report/save path is removed.

Approved behavior:

- remove `cartesian_materialization`, `cartesian_print_summary`,
  `cartesian_print_details`, and `cartesian_save`;
- remove `_pqs_source_box_route_driver_materialization`,
  `_pqs_source_box_route_driver_print_materialization`, and
  `_pqs_source_box_route_driver_save` if they are used only by those wrappers;
- remove a related old White-Lindsey atomic pure-gausslet materialization
  helper only if it becomes uncalled;
- do not add replacement wrappers, adapters, status fields, payload objects, or
  tests.

Evidence: focused search found no hits for the audited names in
`bin/cartesian_ham_builder.jl`. Current CR2-facing artifact workflow uses the
canonical staged producer and `CartesianIDAHamiltonian` artifacts. Live hits are
old wrapper definitions, old tools/harnesses, stale docs-policy assertions, and
stale compact-doc references.

This ID does not approve canonical driver changes, current staged producer
function changes, artifact schema/provenance/reader/manifest changes, route,
shellification, terminal-lowering, raw-block, Residual Gaussian, MWG, IDA, or
Hamiltonian assembly changes, changes to
`pqs_multilayer_complete_core_shell_h1.jl`,
`pqs_complete_core_shell_final_basis.jl`, broad ordinary/Qiu-White donor-kernel
retirement, replacement wrappers, adapters, status fields, payloads, new tests,
or Cr2 workflow.

Failure rule: if any current canonical producer path or public artifact
workflow depends on these wrappers, make no source commit and report the exact
dependency. Do not preserve the wrapper workflow through compatibility
adapters.

### HP-RETIRE-DRV-MAT-TOOL-01 — old wrapper-tool quarantine

Approved tool files:

```text
tools/cartesian_driver_harness.jl
tools/cr2_cartesian_ida_stage_probe.jl
tools/cartesian_driver_ladder_lib.jl
tools/h2_pqs_base_hamiltonian_smoke.jl
```

Approved behavior:

- delete or quarantine only old tools that exist to drive the retired wrapper
  workflow;
- do not move route diagnostics, ladder probing, stage stops, or wrapper
  behavior into `bin/cartesian_ham_builder.jl`;
- do not create replacement tool frameworks.

### HP-RETIRE-DRV-MAT-DOC-01 — active docs cleanup

Approved docs files:

```text
docs/src/developer/algorithm_implementation_index.md
docs/src/developer/designs/cartesian_hamiltonian_producer/implementation_slices.md
docs/src/developer/designs/cartesian_hamiltonian_producer/registry.md
docs/src/developer/designs/cartesian_hamiltonian_producer/r1_public_base_producer.md
docs/src/developer/pqs_manager_running_log.md
```

Approved behavior:

- stop describing the old wrapper workflow as canonical or active current
  authority;
- keep historical references historical;
- keep the canonical staged driver and current producer/artifact path
  unchanged.

### HP-RETIRE-DRV-MAT-TEST-01 — retirement validation

Approved test file:

```text
test/docs/cartesian_ham_builder_policy_runtests.jl
```

Approved behavior:

- remove or update only old route-stage wrapper assertions that assume the
  canonical driver should avoid calling these now-retired wrappers.

Approved validation:

- `git diff --check`;
- package load;
- focused `rg` over `src`, `bin`, `test`, and `tools` showing no remaining
  live references to the retired wrapper names;
- canonical small base artifact/readback smoke;
- canonical small supplemented artifact/readback smoke;
- H2 Residual Gaussian endpoint remains unchanged;
- docs-policy test either passes after update or the stale assertion is removed
  under this authority;
- no Cr2 run.

No new committed test or fixture is approved.

### HP-RETIRE-LADDER-RUNNERS-FN-01 — delete dangling ladder runners

Approved tool files:

```text
tools/run_cartesian_driver_ladder.jl
tools/run_cartesian_line_ladder.jl
```

Approved behavior:

- delete both runner scripts;
- do not add replacements;
- do not modify `bin/cartesian_ham_builder.jl`;
- do not modify `tools/cartesian_driver_ladder_lib.jl` unless a later
  docs-only amendment explicitly approves deleting that quarantined library.

These scripts are only entrypoints into the retired route-driver ladder
workflow after `HP-RETIRE-DRV-MAT-*`. This ID does not approve canonical driver
changes, source changes, test changes except validation scans, artifact/
provenance/reader changes, route/shellification/terminal-lowering/raw-block/
RG/MWG/IDA/Hamiltonian assembly changes, new wrappers, adapters, status fields,
payloads, reports, tools, tests, or Cr2 workflow.

Failure rule: if any live source, canonical workflow, or approved tool still
depends on these runner scripts, make no commit and report the exact
dependency. Do not preserve them through an adapter.

### HP-RETIRE-LADDER-RUNNERS-TEST-01 — ladder runner deletion validation

Approved validation:

- `git diff --check`;
- package load;
- focused `rg` over `src`, `bin`, `test`, and `tools` for
  `run_cartesian_driver_ladder`, `run_cartesian_line_ladder`, and
  `cartesian_driver_ladder_lib`;
- canonical small base artifact/readback smoke;
- no Cr2 run.

No committed test, replacement tool, adapter, or Cr2 workflow is approved.
After this deletion pass, pause the cleanup lane unless a later amendment names
another stale surface.

## Approved For Homonuclear Z-Axis Diatomic Supplemented Workflow

This section approves only the molecule-scope relaxation recorded in
`r3_homonuclear_diatomic_supplemented_workflow.md`. It is generic
homonuclear z-axis diatomic authority, not element-specific Cr2 authority.

### HP-R3U-ZDI-FN-01 — homonuclear z-axis diatomic supplemented facade

Approved source file:

```text
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- replace hardcoded H/Be supplemented guards with explicit homonuclear z-axis
  diatomic validation;
- require explicit atom symbols, nuclear charges, `nup`, `ndn`, geometry, base
  basis parameters, supplement basis labels, and optional supplement
  `basisfile`;
- support exactly two equal-symbol/equal-charge centers on the Cartesian
  z-axis with distinct finite `z` coordinates;
- require neutral all-electron count
  `nup + ndn == round(Int, sum(nuclear_charges))`;
- throw clear `ArgumentError`s for unsupported systems before expensive
  construction where practical.

This ID does not approve heteronuclear systems, non-z-axis/general orientation,
charged systems, ECP, solver/RHF workflow, public API/export redesign,
artifact schema changes, route diagnostics, metadata/status/report fields, or
Cr2-specific branches/defaults/fixtures.

### HP-R3U-ZDI-WIRE-01 — canonical driver supplemented-mode wiring

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

Approved behavior:

- canonical driver `:supplemented` mode may call the supported
  `cartesian_residual_gto_mwg_hamiltonian(...)` facade;
- driver inputs may carry explicit homonuclear z-axis diatomic system, base
  basis, supplement labels, optional `basisfile`, and `hamfile`;
- Cr2 may be an ignored/user-run stress or usability case through the generic
  path only after H2/Be2 validation.

This ID does not approve package-internal helper composition from the driver,
Cr2-specific workflow, committed Cr2 fixtures, route diagnostics, artifact
schema changes, public exports, solver workflow, or broad driver feature
growth.

Line budget for `HP-R3U-ZDI-FN-01` plus `HP-R3U-ZDI-WIRE-01`: at most `100`
added `src`/`bin` lines total, with net simplification expected where
H/Be-specific checks are removed.

### HP-R3U-ZDI-TEST-01 — homonuclear diatomic validation

Approved validation:

- H2 supplemented facade/driver artifact path remains unchanged;
- Be2 supplemented facade/driver artifact path remains unchanged and acts as
  the non-H correctness/performance gate;
- optional ignored/user-run Cr2 stress or usability run after H2/Be2 pass.

No committed Cr2 fixture, committed Cr2 test, new committed test file,
heteronuclear gate, non-z-axis gate, solver run, or artifact schema validation
is approved by this ID.

## Approved Measurement-Only Authority

These entries authorize ignored measurement/probe work only. They do not
authorize production source edits, committed tests, source files, persistent
objects, metadata/report/status/payload fields, artifacts, public API, or Cr2
workflow support.

### HP-COMP-ANGBOX-AUDIT-01 — angular-balanced shellification geometry audit

Approved scope:

- use ignored `tmp/work` probes only to measure z-axis diatomic
  shellification against the angular-balanced molecular box rule;
- report parent axis physical endpoints and counts, snapped nuclear indices,
  core boxes, molecular inner box, each proposed shared-shell expansion,
  transverse physical scale, low/high longitudinal margins from outer nuclei,
  angular-balance ratios, planned non-boundary and boundary z-extension slab
  stacks, and residual outer mismatch if any;
- classify whether the CR2-style thickness-5 axial slabs are planned angular
  z-extension stacks or unexplained fallback outer mismatch;
- recommend a later source lane only if exact files, functions, forbidden
  surfaces, validation, and failure rules are clear.

This ID does not approve production source edits, shellification repair,
thin-slab lowering changes, driver changes, artifact/schema/reader changes,
route skeleton changes, RG/MWG/IDA/Hamiltonian/raw-block changes, public API,
committed tests/fixtures, Cr2-specific workflow, or Cr2 Hamiltonian runs.

### HP-R3REM-AUDIT-01 — remaining exact-operator allocation audit

Approved scope:

- measure the Cr2 q4 R3/RG exact augmented-operator allocation remaining after
  `954c86cd` and the terminal `G-G` product-workspace optimization;
- separate total wrapper allocation from neutral raw-block construction,
  terminal `G-G` product buffers, unit-nuclear `U_GG` Gaussian-sum
  construction, exact augmented nuclear transforms, route/stage setup, and
  audit/replay overhead;
- use ignored `tmp/work` probes only, with H2/Be2 sanity if needed.

Required outcome:

- classify the dominant remaining allocation bucket;
- recommend a future source lane only if the owner, files, functions,
  forbidden surfaces, validation gates, line budget, deletion/simplification
  expectation, and failure rule are specific enough for a separate docs-only
  amendment.

This ID does not extend `HP-R3GG-FN-01` and does not approve unit-nuclear
`U_GG` Gaussian-sum optimization, route/raw-block setup cleanup, final-basis
`G-G` changes, `G-A`/`A-A` raw-block changes, residual Gaussian algorithm
changes, IDA/MWG changes, Qiu-White semantic changes, parent or terminal-basis
changes, persistent caches/workspaces, artifacts, public API/export, Cr2 facade
support, or Cr2 artifact workflow.

## Rejected Or Deferred

### HP-RES-01 — terminal basis build result — rejected

Do not introduce a persistent terminal-basis result wrapper. The realizer
returns `CartesianTerminalBasisRealization` on success.

### HP-CHANGE-01 — return shell overlap from existing shell plan — rejected/deferred

This can be a helper detail under HP-FN-00, but it is not standalone authority.

### HP-OBJ-03 — generic build-result wrapper — rejected

Do not introduce `CartesianHamiltonianBuildResult`, another payload, or a broad
status wrapper around `CartesianIDAHamiltonian`.

### HP-TEST-01 — new committed terminal smoke — rejected

No new committed terminal smoke/probe is approved. Use existing smokes or
ignored `tmp/work` validation unless a later design explicitly approves a test.

### R3 Be2/Cr2 readiness guardrail — deferred

R3-A/B/C are implemented for the narrow H2 residual-GTO/MWG endpoint and
compact supplemented artifact. The following items are closed for that narrow
path and should not be listed as future blockers by default:

- same-construction internal path for the accepted H2 R3 construction;
- deterministic rank-deficient handling in the legacy global-selection
  implementation, now superseded by the approved owner-local selection source
  correction;
- one-shot parent-by-supplement analytic exact-block organization for R3-A
  mixed/self blocks, avoiding repeated CPB-per-terminal-block construction on
  the Be2 proxy;
- independent weight-aware final-basis `V_GM` validation for R3-B;
- compact `supplement_provenance/` artifact group for R3-C.

Do not present Cr2 or broader residual-GTO/MWG supplement support as approved
until a later docs-only amendment chooses and closes the next lane. The
non-exported H2/Be2 R3 usability facade is approved separately by
`HP-R3U-FILE-01`, `HP-R3U-FN-01`, `HP-R3U-WIRE-01`, and `HP-R3U-TEST-01`.
Remaining deferred lanes are:

- implementation and validation of the approved owner-local source correction,
  including the updated H2 scalar and no full Cr2 Hamiltonian;
- public export, public examples, or driver workflow beyond the non-exported
  R3 usability facade;
- a Cr2-readiness lane for measurement-only candidate/rank/memory forecasting,
  with no full Cr2 Hamiltonian yet;
- a basis/supplement-realism lane for validated supplement choices, basis
  labels, and filtering policy beyond the first H2 fixture;
- bounded or streamed residual MWG term storage if higher residual rank makes
  the current dense residual term storage costly;
- allocation-free or bounded-allocation validation reductions for large dense
  matrices.

### Nesting/supplement composition target

The target producer contract is the three-choice composition documented in
`nesting_supplement_composition_plan.md`:

```text
geometry:   atom | z-axis diatomic
nesting:    :pqs | :wl
supplement: off | on
```

This registry section is planning only except for the promoted
`HP-COMP-WLDIAT-*`, `HP-COMP-BASEDIAT-*`, `HP-COMP-SUPPWL-*`,
`HP-COMP-SUPPATOM-*`, `HP-COMP-ATOMBOX-*`, `HP-COMP-WLNS-*`,
`HP-WLDIAT-COMPACT-*`, and `HP-WLDIAT-PARITY-*` pairs above.
Current support remains partial:

- atom / no supplement / `:pqs`: implemented for explicit origin-centered base
  atoms;
- atom / no supplement / `:wl`: implemented for one-center base atoms;
- atom / supplement / either nesting: approved under `HP-COMP-SUPPATOM-*`;
- z-axis diatomic / no supplement / `:pqs`: implemented for explicit
  homonuclear z-axis all-electron inputs;
- z-axis diatomic / no supplement / `:wl`: implemented for explicit
  homonuclear z-axis all-electron inputs through native WL terminal records,
  with compact retained-basis correction approved separately under
  `HP-WLDIAT-COMPACT-*`;
- z-axis diatomic / supplement / `:pqs`: supported for explicit homonuclear
  z-axis diatomics through the residual-GTO/MWG path;
- z-axis diatomic / supplement / `:wl`: supported for explicit homonuclear
  z-axis diatomics through the same residual-GTO/MWG boundary after WL terminal
  realization.

For `nesting = :wl` z-axis diatomics, `HP-COMP-WLNS-*` additionally approves
early rejection of normalized `ns < 4` and records that final retained support
may saturate across working `ns` ranges.
`HP-WLDIAT-COMPACT-*` records that the current boundary-stratum identity
realization is not the production compact retained-basis contract.
`HP-WLDIAT-PARITY-*` records that boundary strata retain the requested shell
count without nucleus-centered symmetric-odd coercion.

Composition IDs:

- `HP-COMP-WLDIAT-FN-01` / `HP-COMP-WLDIAT-TEST-01`: approved WL z-axis
  diatomic base terminal records and artifact path;
- `HP-COMP-BASEDIAT-FN-01` / `HP-COMP-BASEDIAT-TEST-01`: approved base
  homonuclear z-axis diatomic input validation relaxation in
  `src/cartesian_base_hamiltonian.jl`;
- `HP-COMP-SUPPWL-FN-01` / `HP-COMP-SUPPWL-TEST-01`: approved supplemented
  White-Lindsey z-axis diatomic composition through the existing RG/MWG
  boundary;
- `HP-COMP-SUPPATOM-FN-01` / `HP-COMP-SUPPATOM-TEST-01`: supplemented
  one-center atom path through common Residual Gaussian augmentation;
- `HP-COMP-WLNS-FN-01` / `HP-COMP-WLNS-TEST-01`: WL z-axis diatomic `ns`
  early rejection and retained-support saturation wording.
- `HP-WLDIAT-COMPACT-FN-01` / `HP-WLDIAT-COMPACT-TEST-01`: WL z-axis
  diatomic compact retained-basis correction.
- `HP-WLDIAT-PARITY-FN-01` / `HP-WLDIAT-PARITY-TEST-01`: WL boundary-stratum
  retained-count parity correction.

The initial explicit `atom | z-axis diatomic`, `:pqs | :wl`,
`supplement = off | on` composition lanes now all have approved implementation
authority. Remaining geometry, solver, ECP, public export, and Cr2-specific
work still need later docs-only amendments before implementation may begin.
