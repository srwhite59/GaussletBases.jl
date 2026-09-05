# R1 Public Base Producer

Status: implemented exported public facade for unsupplemented, uncorrected,
all-electron Cartesian IDA Hamiltonians. This page is the canonical facade,
input, composition, and base-producer provenance contract. Explicit charged
electron sectors are implemented under `HP-R1-ESECTOR-*`.

## Owned IDs

- `HP-R1-FILE-01` - implemented public producer source owner;
- `HP-R1-FN-01` - implemented exported facade;
- `HP-R1-CORE-FN-01` - implemented public near-nucleus spacing convention;
- `HP-R1-WIRE-01` - implemented report-free staged composition;
- `HP-R1-ART-01` - implemented base-producer provenance sidecar;
- `HP-R1-TEST-01` - implemented standalone public endpoint gate.
- `HP-R1-ESECTOR-FN-01` - implemented explicit charged-sector relaxation;
- `HP-R1-ESECTOR-TEST-01` - completed sector-independence validation.
- `HP-PQS-READER-DOC-01` - implemented reader-facing PQS documentation entrance;
- `HP-PQS-READER-TEST-01` - completed public example and quick-smoke validation.

Exact lifecycle, permission, source, test, and dependency metadata is recorded
in [the registry](registry.md).

## Meaning Of Base

`base` means the Hamiltonian before Gaussian supplements, screened-Hartree or
other corrections, fragment/counterpoise changes, and solver processing. The
output uses the terminal Cartesian basis, exact assembled one-body matrices,
and the producer's localized IDA electron-electron matrix.

Supplements, corrected Hamiltonians, ECPs, arbitrary molecular geometry, and
solver workflows are separate contracts. They must not be added by redefining
this facade.

## Public Interface

The function is exported from `GaussletBases`:

```julia
cartesian_base_hamiltonian(
    system::NamedTuple;
    basis::NamedTuple,
    hamfile::Union{Nothing,AbstractString} = nothing,
)::CartesianIDAHamiltonian{Float64}
```

Implementation owner:

```text
src/cartesian/cartesian_base_hamiltonian.jl
```

Export/include owner:

```text
src/GaussletBases.jl
```

The return is the existing `CartesianIDAHamiltonian{Float64}` directly. The
facade never returns a wrapper, status object, report, payload, or partial
result.

`HP-R1-FILE-01` also maintains the existing root export
`cartesian_base_working_basis`. It is an expert/unstable staged constructor
used by current drivers, source consumers, and public/nested validation. Its
export is maintained, but its returned stage fields are not a stable public
schema and it is not a second base-Hamiltonian facade. PRF-specific wrappers
do not inherit public status from this working-basis export.

## System Input

`system` must contain exactly:

- `atom_symbols::AbstractVector`;
- `nuclear_charges::AbstractVector`;
- `atom_locations::AbstractVector`;
- `nup`;
- `ndn`.

Center collections must have equal lengths. Each location is a fixed
three-element tuple with finite coordinates. Nuclear charges are converted to
`Float64` and must be finite and positive. `nup` and `ndn` must be nonnegative
integers and may not be `Bool`. Unknown or missing keys fail.

### One-Center Atom

The approved atom scope is:

- exactly one center at `(0.0, 0.0, 0.0)`;
- finite positive nuclear charge;
- explicit nonnegative `nup` and `ndn` valid for the realized orbital
  dimension, with positive total electron count under the current
  compatibility container;
- atom symbol used as an explicit provenance label, not as charge, electron,
  basis, or ECP authority.

The broader atom rationale and exclusions are canonical in
[R1 one-center base atoms](r1_one_center_base_atoms.md).

### Homonuclear Z-Axis Diatomic

The approved molecular scope is:

- exactly two centers;
- equal atom-symbol labels and equal finite positive charges;
- both centers on the Cartesian z axis with distinct finite z coordinates;
- explicit nonnegative `nup` and `ndn` valid for the realized orbital
  dimension, with positive total electron count under the current
  compatibility container.

There is no translation or rotation step. Heteronuclear, x/y-aligned,
shifted-parallel, generally oriented, or ECP systems fail validation.

## Electron-Sector Independence

For fixed nuclear charges, nuclear positions, basis input, and numerical
policy, the bare Cartesian construction is independent of `nup` and `ndn`.
The terminal basis, kinetic matrix, unit-nuclear matrices, assembled physical
one-body matrix, IDA interaction, and nuclear repulsion must therefore be
identical when only the electron sector changes. The net molecular charge is
derived rather than independently supplied:

```text
Q_net = sum(nuclear_charges) - (nup + ndn).
```

The compatibility facade continues to require explicit `nup` and `ndn` and
continues to return `CartesianIDAHamiltonian`, which stores those counts as
sector metadata. This pass does not add an optional neutral default,
`net_charge` input, derived-charge field, public accessor, artifact key, or
cache-key convention. Existing neutral calls remain unchanged.

The same invariant composes through the already-implemented supplemented
atom/diatomic path: changing only the sector must not change supplement
loading, residual selection, augmented one-body matrices, or MWG interaction.
This authorizes no supplement algorithm change.

A future architecture may separate a particle-number-independent Cartesian
operator from an electron-sector problem object. That is accepted design
direction only: `cartesian_operator`, `with_electron_sector`, replacement
types, artifact migration, and generalized correction/cache machinery require
separate authority. Sector- or reference-dependent corrections must remain
explicit objects and must not alter the meaning of the bare producer.

## Basis Input

`basis` is a plain `NamedTuple`. It must contain `core_spacing` and at least one
of the supported public size inputs `q` or `ns`. Size values must be positive
integers and may not be `Bool`.

Geometry-specific required fields:

| System | Required extent fields |
| --- | --- |
| one-center atom | `radius` |
| z-axis diatomic | `xmax_parallel`, `xmax_transverse` |

All spacings and extents must be finite and positive.

Implemented optional fields and defaults:

| Field | Default / rule |
| --- | --- |
| `parent_axis_family` | `:G10`; no other family is accepted |
| `reference_spacing` | `1.0` |
| `tail_spacing` | `10.0` |
| `nesting` | `:pqs`; accepted values are `:pqs` and `:wl` |
| `source_span` | `:ordinary`; `:mapped_comx` is PQS-only |
| `s_factor` | `1.0`, finite and positive |
| `coulomb_accuracy` | `:compact`; current source also accepts `:high` |

`nesting`, `source_span`, and `coulomb_accuracy` may be supplied as symbols or
strings and are normalized to symbols before construction.

Public `q` and `ns` semantics belong to
[nesting and supplement composition](nesting_supplement_composition_plan.md).
`q` is the scientific local-order parameter and is the natural lead input for
a PQS-only construction. `ns` is the matched-family parameter and is the
natural lead input when comparing PQS with White-Lindsey. Normalization uses
`q = ns` for PQS and `q = ns - 2` for White-Lindsey; conversely, explicit `q`
implies `ns = q` or `ns = q + 2`, respectively. When both are present, they
must agree with `nesting`.

The existing provenance spelling `ns_source = :legacy_q_compatibility` for an
explicit-`q` call is retained unchanged. Renaming or interpreting that metadata
is a separate compatibility decision and is not part of the public input
classification.

`core_spacing` is the one public near-nucleus physical scale. It is not
`reference_spacing`, box extent, or a hidden element default. One-center
mapping uses the resolved `core_spacing` as its internal mapping `d`.

Public `d` is deprecated compatibility input for one-center atoms only. If
present, it must be finite, positive, and exactly equal to `core_spacing`.
Diatomics reject `d`; public `parent_mapping_d` is unsupported. The expert
mapping-strength exception is separately canonical in
[PQS mapping `s_factor`](pqs_mapping_s_factor.md).

The complete three-tier `:compact | :standard | :high` design belongs to
[Coulomb accuracy policy](coulomb_accuracy_policy.md). At the Pass 378 source
baseline, this facade still validates only `:compact | :high`; `:standard` is
approved there but has not yet landed in `src/cartesian/cartesian_base_hamiltonian.jl`.
This R1 reconciliation does not change that implementation lifecycle.

Unknown basis keys, invalid types, incompatible geometry-specific fields, and
unsupported policy combinations fail rather than being ignored.

## Report-Free Composition

The facade uses the current staged producer, without route reports or old
pair/assembly materialization wrappers:

```text
cartesian_base_working_basis
-> cartesian_base_products
-> cartesian_base_unit_nuclear
-> cartesian_base_vee
-> cartesian_base_hamiltonian_assembly
-> optional base artifact write
```

One resolved input and Coulomb expansion are carried through the parent,
terminal basis, one-body, unit-nuclear, and IDA interaction construction. The
assembly returns the existing Hamiltonian with explicit charges, positions,
and electron counts.

Numerical terminal realization, blockwise one-body construction, localized
IDA, and Hamiltonian assembly are canonical in
[Terminal basis and base assembly](terminal_basis_and_base_assembly.md). This
facade does not duplicate those algorithms or expose their internal objects.

## Artifact Behavior

If `hamfile === nothing`, no file is written. A nonempty path writes the
existing version-1 `:cartesian_ida_hamiltonian` artifact and the facade still
returns the same in-memory Hamiltonian. Empty paths fail; file-system errors
propagate normally. Production does not require readback.

The minimal matrix payload and `read_cartesian_ida_hamiltonian(...)` are owned
by `src/cartesian/cartesian_ida_hamiltonian.jl`. The reader reconstructs the Hamiltonian
from its standard matrix and physical metadata keys and ignores additional
provenance groups.

### `producer_provenance/`

`HP-R1-ART-01` owns one fixed base-producer sidecar in the same file. Current
keys are:

```text
provenance_version        producer
nesting                   route
ns                        q
q_rule                    ns_source
core_spacing              s_factor
reference_spacing         tail_spacing
parent_axis_family        parent_axis_counts
mapping_kind              mapping_d
mapping_s_factor          mapping_s_standard
mapping_s_effective       radius
xmax_parallel             xmax_transverse
atom_symbols              nuclear_charges
atom_locations            nup
ndn                       final_dimension
```

`producer` is `:cartesian_base_hamiltonian`. `route` is derived truthfully
from system kind and nesting. `ns`, `q`, `q_rule`, and `ns_source` preserve the
public-size normalization. Atom artifacts record
`:white_lindsey_atomic_mapping`, resolved `mapping_d`, and `radius`; diatomic
artifacts record `mapping_d = nothing` and their explicit extents.

Mapping-strength provenance records the requested factor, standard value, and
effective value. The sidecar is consumer provenance only; construction stages
must not read it back as numerical input.

The broader `hamiltonian_manifest/`, `recipe_provenance/`, and
`coulomb_expansion/` groups are governed by
[Cartesian Hamiltonian artifact manifest](cartesian_hamiltonian_artifact_manifest.md)
and the Coulomb policy. They are not duplicated R1 schemas.

## Library Facade Versus Driver

The public library facade is the function documented above. The human-facing
canonical driver is the separate script:

```text
bin/cartesian_ham_builder.jl
```

The driver owns editable defaults, trusted input-file/override handling,
terminal due-diligence presentation, coarse stage timing, artifact checks, and
optional supplemented workflow. It currently calls the same staged producer
functions directly so timings remain visible. Driver variables are not facade
keywords, and driver defaults are not hidden library defaults.

Driver authority is canonical in
[Cartesian driver usability workflow](cartesian_driver_usability_workflow.md).

## Reader-Facing PQS Entrance

`HP-PQS-READER-DOC-01` owns the implemented public entrance for the PQS and
White-Lindsey choices. It changes no producer source, export, input, output,
default, or numerical policy. The maintained reader surfaces are:

```text
README.md
docs/make.jl
docs/src/index.md
docs/src/algorithms/cartesian_ida_overview.md
docs/src/algorithms/ida_hamiltonian_and_counterpoise.md
docs/src/algorithms/pqs_shell_construction.md
docs/src/manual/index.md
docs/src/manual/projected_q_shells.md
docs/src/examples/index.md
docs/src/howto/example_guide.md
```

`Projected q-shells (PQS)` preserves the radial workflow as the general
beginner starting point and links rather than repeats the detailed algorithm
pages. Both terminal constructions begin from one common parent and common
physical shell ownership. White-Lindsey uses face/edge/corner products of
one-dimensional contractions. PQS uses a filled source box, boundary product
modes, restriction to shell-owned rows, and only the shell-local symmetric
Lowdin step.

The page states the bounded public geometry: origin-centered one-center atoms
and homonuclear diatomics whose distinct centers lie on the Cartesian z axis.
Its H2+ example at `z = +/-1` bohr is construction smoke evidence, not basis
convergence or publication evidence. SCF and correlated solvers, supplements,
PRFs, screening, and paper-scale campaigns remain consumer or external work.

Rendered navigation must show the new manual page and the existing
`cartesian_ida_overview.md`, `pqs_shell_construction.md`, and
`ida_hamiltonian_and_counterpoise.md` algorithm pages. Their implementation
and retirement status must remain current without changing their mathematics.
README, documentation home, manual, examples landing, and the full example
guide must continue to reach this path. No new docs framework, duplicated
algorithm account, manuscript result, or completeness claim belongs here.

`HP-PQS-READER-TEST-01` owns exactly one public example:

```text
examples/39_pqs_h2plus.jl
```

The example imports only `GaussletBases` and `LinearAlgebra`. It uses
only the exported `cartesian_base_hamiltonian`, `CartesianIDAHamiltonian`,
`one_body_hamiltonian`, and `nuclear_repulsion` interface plus the public
Hamiltonian fields. It constructs two complete unsupplemented Hamiltonians
from one common input, changing only `nesting = :pqs` versus `:wl`:

```text
H/H at z = -1,+1; nup = 1; ndn = 0
ns = 4; core_spacing = 0.6
xmax_parallel = 3.0; xmax_transverse = 2.0
tail_spacing = 2.8; compact Coulomb default
```

It requires dimension `293` for both methods, stored charges/positions and
electron sector, finite one-body and IDA interaction matrices,
`E_nn = 0.5`, symmetry error at most `1e-10`, and lowest-eigenpair residual at
most `1e-10`. It prints each one-body energy but must not assert an energy
ordering or present either value as converged or publication evidence.

The existing quick-example test invokes this file, the focused docs owner
protects its links and public-name boundary, and the Julia-1.10 CI selection
includes the existing `examples` group. No new test file, workflow framework,
other CI change, source/API surface, private call, artifact, solver,
supplement, PRF, screening, or numerical-order gate is part of maintenance.

## Public Cartesian Front-Door Documentation

`HP-PUBLIC-CARTESIAN-FRONTDOOR-DOC-FN-01` owns one bounded correction of the
existing public facade documentation. It may change only the
`cartesian_base_hamiltonian` docstring, the Export reference, and the PQS manual
page. The docstring must state the exact `system` key set, the atom and
diatomic `basis` schemas, all implemented optional defaults and restrictions,
the public `q`/`ns` normalization above, and the atom-only `d` rule. Its
geometry statement must describe origin-centered one-center atoms and
equal-label/equal-charge homonuclear z-axis diatomics; H and H2 are examples,
not the element boundary.

The Export reference keeps its intentional public-`q` examples. The Hydrogen
example must omit `d` or set it equal to `core_spacing`, and no prose may claim
that those values are independent. The PQS manual gains one compact,
copyable PQS-only H2+ construction using `q = 4`; the matched Example 39 keeps
`ns = 4` because `ns` is the comparison parameter. The manual must identify
algorithm-page `Vee` with `ham.electron_electron_ida` and must not present the
construction as a solver or convergence result.

`HP-PUBLIC-CARTESIAN-FRONTDOOR-DOC-TEST-01` adds focused checks only in the
existing docs owner. They must syntax-parse the canonical fenced Julia block,
check its exact system/basis schema and public call shape, and tie its physical
inputs and normalized PQS size to the already-executed Example 39. They must
also protect the corrected `d`, geometry, `q`/`ns`, and `Vee` identities.
Example 39 remains the sole numerical execution for this documentation path;
the docs tests and Documenter must not construct another Cartesian
Hamiltonian.

Preferred/hard implementation budgets are `45/65` added source-docstring
lines, `45/70` added reader-documentation lines, and `20/30` added docs-test
lines, with no new file. No example, definition, method, API, default,
provenance spelling, numerical behavior, workflow, artifact format, version,
tag, or release change is authorized. If accurate copyable documentation or
its structural protection requires another public helper, a numerical doctest,
or broader test machinery, implementation must stop and report the exact gap.

A main-branch implementation updates `/dev/`. `/stable/` remains the immutable
`v0.2.0` documentation until a separately authorized patch release or
documentation-policy decision; this repair grants neither.

## Current Status And Onboarding Documentation

`HP-PUBLIC-STATUS-ONBOARDING-DOC-01` owns one editorial reconciliation of
README, STATUS, DESIGN, ROADMAP, the docs home/manual overview, and the current
ordinary branch page with the released v0.2 interfaces. Link to existing
manuals instead of repeating schemas or tutorials. Inspect the developer
architecture summary and edit only contradictory package/Cartesian claims;
do not redesign its internal architecture account.

Present PQS as the standard current Cartesian construction and default base
route. Preserve supported `q` and `ns` and their accepted relation;
White-Lindsey remains a matched comparison/alternative. Radial construction
remains the recommended beginner path. Separate supported interface contracts
from scientific maturity: publication or release does not establish arbitrary
molecular geometry, convergence, or complete solver support.

Consistently distinguish the released surfaces:

- Base Hamiltonian: unsupplemented, uncorrected all-electron localized IDA for
  origin-centered atoms and equal-label/equal-charge z-axis diatomics, within
  the facade contract, with no general molecular solver claim.
- Reference-density Hartree screening: supplied-field correction assembly in
  the same orthonormal basis/order plus a separate scalar. No automatic field
  generation, fitting, exchange, or SCF claim. PQS and screening remain
  distinct method surfaces.
- External Cartesian GTO transfer: version-1 validated packets, optional
  checkpoint-only PySCF/NumPy exporter, raw projection with capture loss, and
  separate caller-thresholded determinant cleanup. Exporter scope does not
  broaden destination-basis geometry or provide a Hamiltonian/restart.

Remove the completed-repair warning from STATUS, and the obsolete pending
repair and governance ID from the reader branch page. Shorten ROADMAP to at
most 60 lines of current maintenance/research boundaries with no new work
commitments, dates, paper status, or compatibility promises. DESIGN directs
usage readers to manuals/reference before optional roadmap context. Preserve
accurate expert/research distinctions and existing useful links.

Fix README's First radial workflow URL to the existing stable tutorial path.
In `docs/make.jl`, move exactly the First radial workflow and Recommended
atomic setup entries from hidden entries to visible Manual entries after
Overview. Preserve labels/paths, every other entry, and all build settings.
No tutorial-body changes or duplicate navigation entries are authorized.

Across the eight Markdown implementation pages (1,070 baseline lines),
preferred/hard added lines are 180/260, with net non-increasing total lines;
prefer at least 80 lines removed net. Optional architecture edits are capped
at 35 added lines within that total. Navigation is exactly +2/-2 lines.
No new file or test edit. Validate by editorial review, existing docs_fast/full
docs, package load, link/navigation inspection, authority/self-test, generated
view parity, Documenter, and diff checks. Root-document implementation changes
trigger the unchanged full three-gate CI plus Docs; authority/closeout use
docs-only checks. Add no prose/title-equality tests or numerical doctests.

Main corrects `/dev/`; immutable v0.2.0 `/stable/` remains unchanged. Preserve
README's release pin and stable-link policy except the mislabeled tutorial
target. Exclude source/API/numerics, reference expansion, example scheduling,
algorithm pointers, orphan-page deletion, workflows, deployment/version
policy, releases, and other review items. If truthful wording requires an
interface, test, policy, or broader path change, stop and report before commit.

## Failure Behavior

Malformed public requests throw, normally with `ArgumentError`, before
expensive construction where practical. This includes:

- missing or unknown keys;
- non-vector center collections or non-tuple locations;
- mismatched center counts;
- nonfinite/nonpositive charges, spacings, or extents;
- invalid `ns`, `q`, electron counts, geometry, nesting, source span, mapping,
  or Coulomb policy;
- translated atoms, unsupported diatomics, mismatched atom `d`, any diatomic
  `d`, or empty `hamfile`.

Numerical construction and file-system failures propagate. The facade never
converts failure into `nothing`, readiness flags, blocker symbols, or partial
artifacts.

## Validation

The committed standalone gate is:

```text
test/driver_public/cartesian_base_hamiltonian_runtests.jl
```

It validates the exported facade for one-center H and z-axis H2, direct
Hamiltonian return, finite/symmetric matrices, endpoint values, compact/high
Coulomb behavior, atom and molecular geometry rejection, key/type/spacing
validation, deprecated `d` behavior, artifact write/readback, and base
producer provenance. The current H2 fixture has final dimension `487`.

`HP-R1-ESECTOR-TEST-01` additionally approves exact neutral/charged He and
He2 operator parity in this existing file, supplemented He2 parity in
`test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`, and charged
artifact/readback sector preservation. It adds no new test or fixture file.

This is a focused public integration gate, not default tiny-unit coverage. It
must not assert private route-stage inventories or teach the human-facing
driver as the facade implementation.

## Non-Goals

This contract does not authorize:

- supplement-algorithm changes, numerical-complete/protected bases, or
  corrections;
- solver, RHF/UHF, fragment, counterpoise, ECP, or Cr2 workflow;
- translated atoms, arbitrary orientation, or heteronuclear molecules;
- new public symbols, input objects, result wrappers, reports, or status
  payloads;
- new artifact kinds, matrix keys, public provenance readers, or schema
  expansion;
- mapping, nesting, Coulomb, terminal-basis, or driver policy changes.
