# PQS v0.2 Public Surface

Status: implemented and validated under:

- `HP-PQS-PUBLIC-DOC-01`;
- `HP-PQS-PUBLIC-MATCHED-FN-01` and
  `HP-PQS-PUBLIC-MATCHED-TEST-01`;
- `HP-PQS-PUBLIC-SCREEN-FN-01` and
  `HP-PQS-PUBLIC-SCREEN-TEST-01`.

Implementation commit: `058ee54f45c759949f70b54a699ccc318476f8ac`.
Clean replay candidate: `72c46f9ea0dd6b2da7a6a302d34ea1c501d18647`.

This page owns the narrow supported interface intended for GaussletBases
v0.2.0. Exact candidate and release authority are isolated below; the final
conditional tag-and-publication transaction is approved but not yet executed.
The existing numerical
contracts remain in [R1 public base producer](r1_public_base_producer.md),
[matched PQS/White-Lindsey aspect modes](pqs_complete_shell_aspect_source_modes.md),
and [screened-Hartree correction assembly](screened_hartree_correction_assembly.md).

## Purpose And Boundary

The public surface has two independent parts:

1. one bounded H2+ comparison that reconstructs the common full parent and
   matched PQS/White-Lindsey terminal spaces used for Table I; and
2. correction assembly from a consumer-supplied same-basis reference Hartree
   field and represented fixed density.

Neither part exposes the staged producer, atomic packet construction, general
mixed-basis field evaluation, a solver, or paper-local state exactification.
The existing `cartesian_base_hamiltonian` remains the general supported
PQS/White-Lindsey constructor. `cartesian_base_working_basis` remains an
expert/unstable construction surface whose result fields are not public
schema.

## Matched H2+ Public Comparison

The approved exported names are:

```julia
PQSH2PlusRow
PQSH2PlusComparison
pqs_h2plus_comparison
```

The single constructor is:

```julia
pqs_h2plus_comparison(;
    independent_reference_total_energy_hartree::Real,
)::PQSH2PlusComparison
```

The independent reference is required and is used only to form the reported
total errors. The Table I example supplies
`-0.6026342144949465` hartree. The producer does not relabel that external
reference as a repository-computed value.

`PQSH2PlusRow` is a small immutable result row with stable fields:

```text
basis
dimension
capture
electronic_energy_hartree
nuclear_repulsion_energy_hartree
total_energy_hartree
contraction_error_hartree
total_error_hartree
```

The basis values and row order are fixed as `:parent`, `:pqs`, and
`:white_lindsey`. Parent contraction error is zero by definition.

`PQSH2PlusComparison` contains the fixed three rows, the declared independent
reference, and compact common-construction facts:

```text
rows
independent_reference_total_energy_hartree
parent_axis_counts
parent_dimension
direct_core_columns
complete_shell_count
complete_shell_columns
slab_count
slab_columns
```

No internal parent, shell, block, transform, coefficient, due-diligence row,
or fingerprint object is returned. The constructor must fail before returning
unless its live objects establish one common physical construction:

- H2+ at `R = 2` bohr;
- public `ns = 5`, `core_spacing = 0.30` bohr, `s_factor = 1`;
- padding `10` bohr, tail spacing `2.8` bohr, ordinary source span, and
  high135 Coulomb representation;
- route-local PQS `q = 5` and White-Lindsey `q = 3` on the same physical
  parent;
- parent axes `21 x 21 x 29`, dimension `12789`;
- identical parent centers, weights, mapped-axis operators, physical bounds,
  and same-run numerical fingerprints;
- identical physical region keys, shell ownership/supports, direct core, and
  slabs;
- eight complete shells, two slabs, and equal aggregate column accounting
  `275 + 960 + 50 = 1285` for both terminal methods; and
- the current axis-specific White-Lindsey inner construction derived from the
  shared `(ns,ns,L)` shape, never the superseded scalar shell policy.

The full parent solve remains matrix-free. The implementation may call the
accepted staged and apply-based numerical owners internally, but those names
do not become public. It constructs no `Vee`, artifact, or solver state.

### Private Parent-Apply Workspace

Pass 530 authorized and Pass 531 accepted one allocation-focused replacement
inside the private matrix-free parent solve. The Julia `1.12.6` audit measured
`_pqs_h2plus_product_apply` at `34,106` calls, `10.652` seconds, and
`21.914 GiB` of cumulative allocation. The `248` `raw_nuclear` calls accounted
for `10.589` seconds and `21.538 GiB` inclusively; their `33,480` nested product
applications accounted for `21.512 GiB`. These figures overlap and must not be
added.

Commit `fbef95e60ca3aadafe2082b4a63f8522a724e7be` replaced the allocating
three-axis action with a file-private mutating action using explicit
caller-supplied output and scratch storage. One invocation of
`_pqs_h2plus_parent_solution` allocates and exclusively owns its complete
workspace. Every output and intermediate is fully overwritten before use,
active input, output, accumulator, and scratch arrays must not alias, and no
view backed by scratch may escape. `raw_nuclear` may reuse one product output
while accumulating all `135` terms into a distinct caller-owned result. The
orthogonalization, kinetic, final callback, and capture callers may be adjusted
only as needed to consume the same replacement without retaining a second
allocating hot path.

Workspace ownership is lexical and per comparison call. No mutable global or
module constant, task-local store, pool, public workspace, retained result
field, metadata entry, or hidden persistent state is permitted. Concurrent
public comparison calls must own disjoint storage. The current Lanczos owner
invokes its nonescaping callback synchronously; the implementation must not
change that owner or rely on cross-call sharing. If safe reuse requires public
API expansion, another source owner, callback/lifecycle changes, locking, a
cache, or a reusable storage framework, implementation stops without a source
commit and returns the exact obstacle.

This is a storage-only optimization. Preserve the existing axis-product order,
matrix order, Gaussian-term order and signs, Float64 behavior, and complete
operator action. Deterministic old/new product-action probes must agree exactly.
The parent/PQS/White-Lindsey dimensions, parent and terminal fingerprints,
`275 + 960 + 50` accounting, shell/slab topology, bounds, captures, energies,
eigen-residuals, symmetry checks, due-diligence warnings, and all release
tolerances must remain unchanged.

Performance acceptance used Julia `1.12.6`, one Julia thread, and eight BLAS
threads, with package loading excluded. The warmed isolated action allocated
`1088` bytes per call and no parent-sized storage. Paired clean baseline and
candidate comparisons reduced cumulative allocation by `21.606 GiB` fresh and
`21.963 GiB` warm; warmed elapsed time improved from `33.074` to `32.795`
seconds. A transient two-task probe used independent workspaces and matched the
serial product actions exactly. No benchmark, baseline helper, or performance
fixture is committed. Maintenance must preserve the absence of parent-sized
per-action allocation and must not reintroduce cross-call storage.

Maintenance is confined to `src/pqs_matched_h2plus.jl`. The accepted
implementation added/deleted `55/32` source lines, stayed within the `60`-line
hard and `30`-line net limits, and removed the allocating private helpers rather
than retaining parallel implementations. No test edit, new file, type, export,
public binding, dependency, metadata, workflow, example, fixture, or release
change is authorized. Package load, the unchanged `18/18` matched-H2+ owner,
normal public CI, terminal due-diligence inspection, and repository
authority/docs gates remain the maintenance owners.

The supported example writes a TSV with exactly the stable row fields above.
Additional timing, allocation, Git, path, shell-ledger, or report fields from
the private paper driver are not part of this API.

## Screened-Hartree Public Assembly

The approved exported names are:

```julia
ExactRepresentedHartreeField
FittedReferenceHartreeField
ScreenedHartreeCorrection
screened_hartree_correction
screened_hartree_delta_one_body
screened_hartree_energy_constant
screened_hartree_consistency_error
screened_hartree_field_kind
```

The field constructors are:

```julia
ExactRepresentedHartreeField(
    reference_hartree_field,
    reference_coulomb_self_integral;
    provenance,
)

FittedReferenceHartreeField(
    fitted_reference_hartree_field,
    density_coulomb_self_integral;
    provenance,
)
```

Both matrices are copied finite symmetric matrices. Provenance is a required,
nonempty descriptive string. The two concrete types are intentionally
distinct: a fitted potential field cannot silently enter the exact represented
route.

The public assembly entry point is:

```julia
screened_hartree_correction(
    V_IDA,
    field,
    reference_coefficients,
    occupations;
    representation_atol = 1.0e-8,
    density_nonnegativity_atol = 1.0e-12,
    symmetry_atol = 1.0e-10,
    closure_atol = 1.0e-8,
)::ScreenedHartreeCorrection
```

`V_IDA`, the reference-field matrix, and the rows of
`reference_coefficients` all use the same orthonormal one-particle basis and
ordering. The columns of `reference_coefficients` are represented reference
orbitals. `occupations` are finite, nonnegative, spin-summed occupations;
fractional values are allowed.

Define the spin-summed reference density and localized charges by:

```text
P0 = reference_coefficients * Diagonal(occupations) *
     reference_coefficients'
q0 = diag(P0)
```

`representation_atol` owns orbital orthonormality, trace, and occupation-sum
agreement. `density_nonnegativity_atol` is only the numerical allowance for a
small negative entry of `q0`; it is not a neutrality or charge-state test.
The reference particle number is `sum(occupations)`. No relationship to the
sum of nuclear charges is required.

Public documentation names the exact constructor's scalar input
`reference_coulomb_self_integral` and defines it as:

```text
(rho0|rho0) = Tr(P0 * J0)
```

It is twice the reference Hartree energy. The fitted constructor's
`density_coulomb_self_integral` is the corresponding density-representation
quantity; the fitted potential does not redefine it.

The returned correction is:

```text
Delta_J0 = J0 - Diagonal(V_IDA * q0)
C        = 0.5 * q0' * V_IDA * q0 - 0.5 * (rho0|rho0)
```

`Delta_J0` and `C` are added to direct electron-electron/Hartree energy
accounting. They do not replace or mutate the physical kinetic-plus-nuclear
one-body operator. The scalar remains a separate total-energy constant. A
consumer may represent it as `(C/N)I` only inside a fixed nonzero
particle-number sector `N`; that optional rewrite is not provided by v0.2 and
is not valid across particle sectors.

The stable accessors return the correction matrix, scalar, signed consistency
error, and field kind. The signed convention is always:

```text
screened_hartree_consistency_error =
    Tr(P0 * J0) - (rho0|rho0)
```

For `ExactRepresentedHartreeField`, its magnitude must not exceed
`closure_atol`. Exact means exact for the declared represented density and
Coulomb representation, not continuum exactness. For
`FittedReferenceHartreeField`, the same signed value is a reported potential-
fit approximation and is not by itself a rejection condition. Finiteness,
symmetry, dimensions, represented-density validity, and the derivative/action
identity remain hard failures for both routes.

The result type is supported through the four accessors above. Its internal
field layout, packet summary, and diagnostic `NamedTuple` are not public
schema and receive no serialization promise.

This v0.2 surface assembles a correction from a supplied same-basis field. It
does not publicly construct atomic density fits, fit a potential, evaluate a
fitted potential on a terminal-plus-supplement basis, place atomic packets, or
form additive molecular self/cross energies. Those capabilities remain
internal under their existing contracts.

No exchange field or exchange correction is returned. The words "exact" and
"fitted" classify the supplied direct Hartree field only.

## User Documentation And Layering

The user documentation must distinguish:

1. stable PQS/White-Lindsey construction through
   `cartesian_base_hamiltonian` and the bounded matched H2+ comparison;
2. stable correction assembly from a supplied same-basis reference Hartree
   field; and
3. paper-local Ximg/XHF finite-reference exactification.

Ximg/XHF remains provenance-only. It is neither reconstruction of a physical
mixed-basis four-index interaction nor transferable exchange/correlation
theory. No Ximg/XHF name, constructor, result, or accessor is exported.

The implemented user surfaces are:

```text
docs/src/manual/projected_q_shells.md
docs/src/manual/reference_density_hartree_screening.md
docs/src/manual/index.md
docs/src/examples/index.md
docs/src/howto/example_guide.md
examples/40_screened_hartree_fixed_density.jl
examples/41_pqs_h2plus_table1.jl
```

The screening manual page and both examples are source-backed and may be
rendered as executable public documentation. Their clean candidate replay is
accepted under the protocol below. This validation does not itself create a
version bump or release tag.

The screening example uses two bounded pure-GTO orbitals and an explicitly
declared one-center, two-electron fixed density. It obtains the direct IDA
matrix from a small `CartesianIDAHamiltonian`, keeps
`one_body_hamiltonian(ham)` separate and unchanged, and obtains the accurate
represented reference field from the existing bounded pure-GTO Coulomb
oracle. It demonstrates energy and occupied-action closure with a nontrivial
correction. It is an executable assembly example, not reproduction of the
historical He calculation or a new He endpoint.

## Reader-Facing Documentation Parity

`HP-PQS-PUBLIC-DOC-PARITY-FN-01` maintains the reader-document repair
implemented by `c78defcc9b299ee5f32cf42910812f5581657d93`. The curated export
reference resolves all eleven accepted v0.2 names:

```text
PQSH2PlusRow
PQSH2PlusComparison
pqs_h2plus_comparison
ExactRepresentedHartreeField
FittedReferenceHartreeField
ScreenedHartreeCorrection
screened_hartree_correction
screened_hartree_delta_one_body
screened_hartree_energy_constant
screened_hartree_consistency_error
screened_hartree_field_kind
```

Documenter's `@docs` resolution requires attached Julia docstrings. Concise
public docstrings now exist for the ten names that previously lacked them:
both result types, both field types, `ScreenedHartreeCorrection`, the
`screened_hartree_correction` generic covering both accepted overloads, and
the four accessors. The pre-existing `pqs_h2plus_comparison` docstring remained
unchanged.

Those docstrings must describe only the accepted interface on this page. The
row and comparison docstrings may name their stable public fields. The field
docstrings must preserve exact-versus-fitted route identity, copied matrix
ownership, and Coulomb-self-integral terminology. The correction constructor
must state the common orthonormal-basis convention, spin-summed density with
fractional occupations, and supplied-field boundary. Its signed consistency
error remains `Tr(P0 * J0) - (rho0|rho0)`. The correction-result docstring must
direct callers to the four accessors rather than make its internal fields
public schema. No definition, signature, field, export, dispatch, executable
expression, numerical behavior, or compatibility promise may change.

The two rendered nested-diatomic algorithm pages must preserve their private
historical provenance only as path-free prose. They must expose no personal
absolute path, Dropbox path, invented URL, imported private file, or new
bibliographic claim. `examples/28_ordinary_one_body_fidelity.jl` remains
unchanged and must be labeled explicitly as an internal, non-contractual
diagnostic in both reader guides that present it. This does not promote or
retire the example.

`HP-PQS-PUBLIC-DOC-PARITY-TEST-01` owns only a focused extension of the
existing docs consistency test: use one fixed eleven-name list to require
agreement between the existing root exports and curated reference entries,
reject personal paths in the two rendered pages, and require the example-28
internal classification. Documenter remains the executable `@docs` resolution
gate. No new test file, parser, generalized export inventory, example
execution, or broad documentation scan is authorized.

The accepted implementation added `80` source lines, all docstrings in the two
existing owners, `23` reader-documentation lines, and `23` lines in the
existing docs test. It added no new file, test file, metadata, executable
behavior, or public name. Maintenance remains confined to the same owned paths
and fixed eleven-name inventory; broader public-documentation or API work
requires a separate amendment.

## Validation Contract

Ordinary pull-request coverage uses small dimensions and requires:

- exported-name and constructor smoke coverage;
- exact structural dimensions and spin-summed reference charge;
- matrix symmetry at `1e-10`;
- energy and occupied-action closure at `1e-11`;
- a nonzero `Delta_J0` in the two-function fixture;
- the signed exact consistency error within `1e-11`;
- a fitted-field control that reports rather than rejects a finite nonzero
  signed consistency error;
- unchanged physical one-body data; and
- malformed dimensions, asymmetry, nonfinite data, negative occupations,
  representation loss, and exact-route closure failure rejected with
  `ArgumentError` or `DimensionMismatch` as appropriate.

The existing small `examples/39_pqs_h2plus.jl` remains the ordinary PQS/WL
construction smoke. The new fixed-density screening example may join the
existing quick-example set if its measured runtime remains bounded.

A separate release/nightly test owns the full Table I comparison. It requires:

- exact parent/PQS/White-Lindsey dimensions and structural counts;
- exact same-run parent and physical-shell identity between PQS and
  White-Lindsey;
- capture absolute tolerance `2e-8`, relative tolerance zero;
- electronic, total, contraction-error, and total-error absolute tolerance
  `1e-7` hartree, relative tolerance zero;
- same-run parent residual and capture closure at `1e-9`;
- same-run symmetry at `1e-10`; and
- no cross-platform raw-byte or matrix-fingerprint equality requirement.

Commit `b0dbd9ea37317590334a24883ef0667bdb0195a5` preserves this complete
`18`-assertion contract while removing its duplicate Example 41 execution. The
standalone example returns the same comparison it prints and writes, and the
release owner applies these checks to that returned value. The successful
release path therefore constructs the complete comparison exactly once;
nonfinite-input rejection remains an early failure check rather than a second
successful construction.

The slower test is not part of every local edit gate. It is run by a focused
read-only release/nightly workflow and manually before the candidate release
is frozen.

## Clean Candidate Replay Environment

GaussletBases is a Julia package, and its root `Manifest.toml` is intentionally
ignored. A committed root manifest is not required for the candidate replay
and must not be added solely for this campaign.

The clean replay starts from one exact detached candidate commit with no
private manifest and an empty machine-local depot. It may run
`Pkg.instantiate()` against the committed root `Project.toml`, allowing Julia
to perform one fresh compatible resolution. Before any example or release test
runs, preserve the generated `Manifest.toml` outside the repository together
with:

- the exact candidate commit and tracked-clean status;
- Julia version, platform, and `versioninfo()` output;
- the complete committed `Project.toml` SHA-256;
- the complete generated `Manifest.toml` SHA-256, `project_hash`,
  `julia_version`, and full bytes;
- registry name/UUID and available revision identity;
- the instantiate command and log; and
- hashes of the example outputs and focused release-test log.

The generated ignored manifest may remain in the detached worktree while the
replay runs, but it is environment evidence rather than a source change. The
tracked tree must remain unchanged and pass `git diff --check`. Do not copy a
private ignored manifest from another checkout, resolve inside a Dropbox depot,
commit the generated manifest, or reinterpret the fresh resolution as proof
that every future dependency resolution is numerically identical. The archived
manifest pins this candidate replay; lower-bound Julia compatibility remains a
separate CI fact.

Candidate `72c46f9ea0dd6b2da7a6a302d34ea1c501d18647` passed this protocol on
Julia `1.12.6` for arm64 macOS. The archived fresh resolution has manifest
SHA-256 `28d952b22ab0685bcc56fb11adb150fdbb0b7faf46e79a7a9f76ecf554f2c342`
and project hash `4036d1ee56daf08aad6eb3a169bfb75510483a9e`; the General registry tree
was `d274819867e3891d366fad62967ab1b7d1ca283a`. Both public examples and
the focused `18/18` release gate passed with a tracked-clean candidate tree.
This accepts the bounded replay only; release review, versioning, and tagging
remain separate decisions.

## Dependency Compatibility Lifecycle

`HP-PQS-PUBLIC-COMPAT-FN-01` maintains the implemented package-metadata
declaration. `HP-PQS-PUBLIC-COMPAT-TEST-01` records its completed acceptance
evidence but grants no test or workflow edit. The root `Project.toml` must
retain `julia = "1.10"` and exactly these six direct-dependency ranges:

```toml
JLD2 = "0.6.4"
LinearAlgebra = "1.10"
SHA = "0.7"
SparseArrays = "1.10"
SpecialFunctions = "2.8"
TOML = "1"
```

These are Julia package compatibility ranges, not exact dependency pins.
Julia `1.10` remains the minimum supported line; Julia `1.12.6` remains the
canonical current replay environment. Package-version ownership is separate
from dependency compatibility. The candidate-preparation decision below may
change the package version without changing any compatibility range or the
Julia floor.

Commit `873b6e27a0f3e049fda70b4ff16dc0682354efce` added exactly those six
lines. A fresh Julia `1.12.6` resolution selected JLD2 `0.6.5`,
LinearAlgebra `1.12.0`, SHA `0.7.0`, SparseArrays `1.12.0`, SpecialFunctions
`2.9.0`, and TOML `1.0.3`; every version satisfies its declared range. The
public examples and numerical gates remained unchanged, Julia `1.10` CI
passed, and no root manifest became tracked.

Maintenance is limited to preserving these six ranges and the Julia `1.10`
floor. It must not change source code, exports, tests, examples, documentation
content, workflows, the docs project, or manifest policy. No root manifest may
be added or tracked. Any widened range, lowered floor, added compatibility key,
tag, or release requires a separate docs-only decision; package-version changes
require their own candidate/release record.

Compatibility changes require a fresh Julia `1.12.6` resolution and
instantiation from the edited committed project in an empty machine-local
depot. Archive the new generated manifest, project and manifest hashes,
`project_hash`, runtime and
registry identity, direct dependency versions, commands, logs, and clean-tree
state. Verify that every resolved direct dependency satisfies its declared
range, then run package load, both bounded public examples, the focused release
gate, docs resolution/build, authority check/self-test, docs tests, and the
bounded numerical groups. Existing numerical tolerances govern output parity;
raw-byte identity is not required. Julia `1.10` CI remains the lower-line gate.

Stop without widening or lowering a bound if either fresh resolution, package
load, public numerical output, docs, or Julia `1.10` CI fails. General
registration, a root manifest, candidate publication, tagging, and release
publication are not authorized by the compatibility record; their exact
lifecycle owners are separated below.

## v0.2.0-rc1 Candidate Preparation

`HP-PQS-PUBLIC-RC1-FN-01` owns maintenance of the prepared candidate identity
at commit `41fa897ae919510b179a425027a8ce2d4a2167b3`. The accepted
implementation made exactly two repository changes:

1. change the root `Project.toml` version from `0.1.0` to `0.2.0-rc1`; and
2. add one concise root `CHANGELOG.md`.

The changelog is reader-facing and grouped under `Added`, `Changed`, `Fixed`,
`Public-surface reduction`, and `Scope`. It must summarize only:

- the bounded matched-PQS/White-Lindsey H2+ types and constructor, supplied-
  field screened-Hartree types/constructor/accessors, examples 39-41, and
  bounded CI/release validation and versioned documentation support;
- dependency compatibility, current Cartesian/nested documentation and API
  reference, and least-privilege `/dev/` plus version-tag deployment;
- legal unsplit-H2 packet construction, rectangular/endcap provenance, and
  forwarding of existing endcap policy, `q`, and `L` diagnostics;
- removal of six PRF-specific root exports while preserving their private
  diagnostic/provenance implementation; and
- the package boundary: basis/operator construction and bounded screening
  assembly are included, while external SCF, correlated solvers, and paper-
  scale campaigns are not.

These claims trace respectively to the accepted reader/public records and
commits `01f0fe002`, `058ee54f4`, `873b6e27a`, `c78defcc9`, `abee269ee`,
`6a3656991`, and `820c2fc8b`. The implementation must also compare the fixed
eleven-name curated reference with live root exports. The changelog must not
contain manager-pass chronology, internal retirement inventories, detailed
fingerprints or energies, private paper-validation history, manuscript
headline claims, a publication date, DOI, journal status, final paper citation,
or content imported from a paper workspace. `CITATION.cff` remains blocked on
separate approved wording.

The resulting changelog is `39` lines, within its preferred `45`-line and hard
`70`-line bounds. The preparation changed no source, test, example, dependency,
workflow, API, scientific document, badge, stable link, or homepage, and no
root manifest became tracked.

`HP-PQS-PUBLIC-RC1-TEST-01` is completed evidence and grants no file edit. The
accepted preparation passed:

- exact TOML parsing with package version `0.2.0-rc1` and unchanged dependency
  compatibility;
- a fresh Julia `1.12.6` GitHub installation from the exact pushed preparation
  commit, followed by package load;
- public examples `01_first_gausslet.jl`, `39_pqs_h2plus.jl`,
  `40_screened_hartree_fixed_density.jl`, and `41_pqs_h2plus_table1.jl`, plus
  the focused H2+ release test passing `18/18`;
- successful Julia `1.10` CI run `32279994531`, authority check/self-test,
  generated-view parity, docs tests, Documenter, manager-log bound, and
  `git diff --check`;
- simulated `v0.2.0-rc1` documentation classification and canonical path, with
  standard prerelease behavior not selecting `stable`; and
- a clean `git archive` with SHA-256
  `afc7fcfbe9871272a70b06c06340a99d60346aa44ed2d5621cdf88195e2c0a65`,
  containing `Project.toml` and `CHANGELOG.md` but no root manifest or untracked
  handoff.

The focused H2+ release test passed `18/18` with dimensions
`12789/1285/1285`; examples 01, 39, 40, and 41 also passed. The implementation
commit's Docs run failed only because the newly tracked changelog still had
planned-path authority; this lifecycle closeout resolves that expected state
transition without changing the candidate files. The completed annotated tag
and approved GitHub prerelease each have separate lifecycle owners below.
General registration, citation metadata, and final-release authorization still
require separate explicit decisions.

## v0.2.0-rc1 Annotated Tag Lifecycle

`HP-PQS-PUBLIC-RC1-TAG-FN-01` records the completed tag-only publication.
The immutable target is the accepted candidate closeout commit
`1546c18d3058cce2b5051b50788cda3c12585e51`, with tree
`1b53a9eb51d11cfc31b8b0356349c62f0de8915f` and full tracked `git archive`
SHA-256 `2a0b6938d3b341900d73668e7f0644c34b8a851e1b823356c53c2866fd19522a`.
The later authority commit is deliberately not the tag target.

The frozen target contains package version `0.2.0-rc1`, the accepted 39-line
changelog and six dependency bounds, the accepted eleven-name public surface,
reader parity, examples and tests, and tag-aware documentation deployment.
Relative to implementation commit `41fa897ae919510b179a425027a8ce2d4a2167b3`,
only the six accepted Pass 493 lifecycle documents changed; `Project.toml`,
`CHANGELOG.md`, source and exports, tests and examples, dependency declarations,
the Docs workflow, `docs/make.jl`, and reader-facing documentation are
byte-identical.

The tag operation is complete. Annotated tag object
`a4284f0bf448fb9d717de26ccbe1e9fc16db5ed2` exists locally and remotely and
peels to the frozen target. Its tree and archive match the frozen values.
Tag-triggered Docs run `32295705338` passed, `/v0.2.0-rc1/` is live with the
exact canonical URL, `/stable/` remains absent, and `/dev/` is byte-identical
to its pre-tag state. No further tag mutation is authorized.

`HP-PQS-PUBLIC-RC1-TAG-TEST-01` is completed. Commit
`31caa87d3b83599de7f7295678ee599209113552` added the exact RC1 self-mapping
under `HP-PQS-DOCS-TAGDEPLOY-FN-01/TEST-01`. Main-deployment Docs run
`32302304167` passed; live `versions.js` lists RC1 and `dev`, `/stable/`
remains absent, and both canonical folders remain intact. Documenter's
internal `DOCUMENTER_STABLE` fallback names RC1 because no final release
exists, but no `stable` selector entry, symlink, path, or release status was
created. The repair did not move, replace, delete, or recreate the tag.

The tag lifecycle is closed. The exact RC1 GitHub prerelease operation is
separately bounded below. General registration, `CITATION.cff`, final
`v0.2.0`, homepage or stable-link edits, a custom release framework, and
tracked source/API/dependency/workflow/numerical/manuscript mutations remain
unauthorized.

## v0.2.0-rc1 GitHub Prerelease Lifecycle

`HP-PQS-PUBLIC-RC1-RELEASE-FN-01` records the completed GitHub prerelease for
the existing immutable `v0.2.0-rc1` tag. Published release `373460389` has:

- tag `v0.2.0-rc1`, verified rather than created by the release operation;
- title `GaussletBases v0.2.0-rc1`;
- prerelease status `true`, draft status `false`, and latest status `false`;
- no uploaded assets, so only GitHub's automatic source archives appear; and
- the exact release body below.

~~~~markdown
GaussletBases v0.2.0-rc1 is a public-testing prerelease of the GaussletBases
v0.2 package surface.

## Projected q-shells (PQS)

- Bounded public PQS/White-Lindsey H2+ comparison types and constructor.
- Public examples 39 and 41 cover the supported PQS construction and bounded
  comparison.

## Reference-density Hartree screening

- Supplied-field correction assembly with stable public result accessors.
- Public example 40 demonstrates the bounded screening algebra.

PQS and reference-density Hartree screening are distinct method surfaces.
Screening is not a layer of the public PQS workflow, and this prerelease does
not present them as one combined method or one combined-paper validation.

## Compatibility and scope

- Julia 1.10 is the supported minimum; the canonical clean replay used Julia
  1.12.6.
- Six PRF-specific wrappers remain available internally for diagnostics and
  provenance but are not exported.
- GaussletBases supplies basis/operator construction and bounded screening
  algebra; it is not a complete SCF or correlated molecular solver.

## Install

Install this prerelease directly from its immutable tag:

```julia
using Pkg
Pkg.add(url = "https://github.com/srwhite59/GaussletBases.jl", rev = "v0.2.0-rc1")
```

- [Versioned v0.2.0-rc1 documentation](https://srwhite59.github.io/GaussletBases.jl/v0.2.0-rc1/)
- [Changelog at the immutable tag](https://github.com/srwhite59/GaussletBases.jl/blob/v0.2.0-rc1/CHANGELOG.md)
~~~~

Because the body contains an inner Julia fence, the exact text between the
outer fence markers was copied into a temporary notes file outside the
repository. The completed noninteractive operation was:

```bash
gh release create v0.2.0-rc1 \
  --repo srwhite59/GaussletBases.jl \
  --verify-tag \
  --title "GaussletBases v0.2.0-rc1" \
  --prerelease \
  --latest=false \
  --notes-file /private/tmp/gaussletbases_v0.2.0-rc1_release_notes.md
```

No filename or pattern followed the tag; no generated notes, discussion,
draft, or `--target` was used. `CHANGELOG.md` remained a read-only tracked body
source. Preflight verified synchronized `main`, only the two established
untracked handoffs, release absence, exact local/remote tag identity and peel,
package-version/tag agreement, and immutable links before publication.

`HP-PQS-PUBLIC-RC1-RELEASE-TEST-01` is completed and grants no edit or
publication authority. GitHub API and rendered-page checks confirmed the
exact tag, title, `1455`-byte body, `draft=false`, `prerelease=true`, no
uploaded assets, and automatic source archives tied to the accepted tag. The
latest-release endpoint remains absent. Both archives validated, and a fresh
isolated Julia `1.12.6` installation loaded GaussletBases `v0.2.0-rc1` with
manifest tree `1b53a9eb51d11cfc31b8b0356349c62f0de8915f`. RC1 and `/dev/`
remain live with exact canonical URLs, `versions.js` lists RC1 and `dev`, and
`/stable/` remains absent. Evidence is preserved under
`/Users/srw/dmrgtmp/release_replays/gaussletbases_v0p2_rc1_release_2026-08-19/`.

The release lifecycle is closed. Do not edit, delete, recreate, retarget, add
assets to, or silently retry the accepted release. Publication changed no
tracked file.

No General registration, `CITATION.cff`, final `v0.2.0`, paper title, DOI,
journal status, manuscript headline number, publication-timing promise,
combined-paper framing, custom asset, binary, manifest, paper workspace,
validation log, data bundle, repository setting, metadata edit, tag mutation,
stable link, source/API/test/example/dependency/workflow/numerical change, or
broader release framework is authorized. Future citation metadata must permit
separate PQS-paper, screening-paper, and package-software citations; exact
titles and preferred citations remain deferred.

## Package Export Integrity

`HP-PUBLIC-EXPORT-INTEGRITY-FN-01/TEST-01` own the implemented bounded cleanup
of package export declarations discovered after the RC1 release. At baseline
`41ab3e13121f5af1e145775500e91f9ac61c9760`, runtime inspection of the root
module and every direct package-owned child module finds exactly four defects:

- root export `CartesianBasisBundle3D` is undefined;
- root export `nested_fixed_block_timing_report` names a bare generic with no
  methods;
- internal `CartesianRouteCore` exports `final_units` and `unit_keys`, neither
  of which is defined by that module.

No committed source, test, driver, tool, or example called these names. Commit
`b72500f7e619db5875918e3290ed2b306be51f43` removed the two root export
entries, the bare `nested_fixed_block_timing_report` generic declaration, and
the two internal route-core export entries without replacement. The exact
source delta was `+1/-6`, net `-5`; the sole added line retained valid export
`lowering_recipe` without the trailing comma that became syntactically invalid.
Maintenance keeps these invalid names and the bare generic absent while
preserving valid exports and all numerical behavior. It adds no replacement
definition, alias, shim, deprecation, helper, dependency, metadata, or public
name.

One regression in existing owner `test/core/runtests.jl` audits
`GaussletBases` plus each defined direct child module discovered with
`names(GaussletBases; all=true, imported=false)` where
`parentmodule(child) === GaussletBases`. For every exported name
returned by `names(module; all=false, imported=false)`, it requires a defined
binding; if that binding is a `Function`, it requires at least one method. This
checks export integrity without freezing an export inventory or traversing
external/imported modules. The accepted `22`-line test passed `1386/1386`
dynamic checks. It remains the compact maintenance gate; no export snapshot,
new test file, fixture, compatibility assertion, or numerical assertion is
authorized.

Package load, core, IDA, Cartesian, examples, docs `87/87` and `10/10`, the
authority self-test, Documenter, and diff checks passed. GitHub CI run
`32444750148` and Docs run `32444750177` passed on the implementation commit.

Commit `6e7bcbb7dae4e865dbdc0362b8f39ffd23f0a468` implements the bounded
post-v0.2 continuation under the same function record. It deletes the unused
four-line `TimedNestedFixedBlockBuild` definition and removes its root export.
It also removes `OneCenterAtomicNestedLayerStructure` and
`QiuWhiteResidualGaussianOperators` from the root export while retaining the
unchanged internal diagnostic-storage type and one-line alias to
`OrdinaryCartesianOperators3D`. The exact delta is `+0/-7` across
`src/GaussletBases.jl` and `src/cartesian_nested_faces.jl`. The root-surface
reduction is the useful cleanup; retaining the alias avoids compatibility risk
without meaningful carrying cost. Ignored historical probes and conflicted
copies remain untouched.

The existing dynamic export-integrity regression remained unchanged. It proved
`TimedNestedFixedBlockBuild` absent, both retained names defined but
unexported, and the undocumented exported-binding backlog reduced from `74` to
`71`. The complete core owner, including existing nested/QW regressions,
passed. Docs passed `115/115` plus `10/10`; CI run `33143556046` passed all
three gates and Docs run `33143556067` passed. Package load, authority
check/self-test, Documenter, and diff checks also passed.

Maintenance adds no shim, warning, replacement API, alias deletion, docstring,
test, helper, metadata, file, or numerical behavior. `ShellLocalAngularProfileKey`,
`@timeg`, `CuratedSpherePointSet`, `LegacySGaussianData`, and
`QiuWhiteHybridOrbital3D` remain unchanged. Package versions, immutable v0.2.0
tags/releases, workflows, and release state also remain unchanged.

## Package-Root Data And Include Indirection

`HP-PACKAGE-ROOT-PATH-FN-01` maintains the source-layout Step 0 implemented in
commit `cadf02f7c6785fdf8fa792358837c5e83010376b`: package data and one nested
include are independent of the directory containing their current source file.
The 2026-09-02 external follow-up review and metrics ledger remain audit
evidence, not authority. The fresh baseline audit at
`5b844e619b317633d6747974af54d8893df653f7` corrected their count to five data
paths plus one nested include.

The implementation owns exactly one private package-root constant in
`GaussletBases` and one small private data-path helper:

```julia
const _PACKAGE_ROOT = normpath(joinpath(@__DIR__, ".."))
_package_data_path(parts...) = normpath(joinpath(_PACKAGE_ROOT, "data", parts...))
```

The five existing radial, angular, and legacy path functions must use that
helper without changing their signatures or callers. The nested
`_RadialPrototypeHighPrecFamilies` module does not inherit private parent
bindings. Its `families_high_prec.jl` include must therefore refer to the same
root through an explicit parent-qualified expression such as
`parentmodule(@__MODULE__)._PACKAGE_ROOT`, not through an assumed unqualified
binding. Include order remains unchanged.

All six targets and their baseline SHA-256 values are fixed:

| Target | SHA-256 |
| --- | --- |
| `data/radial/paper_parity_g10_k6_x2.jld2` | `8821715de3b4490049295afa25aaf4cd340b37e95a90d71da00e48fe99abe462` |
| `data/angular/SpherePoints.jld2` | `c73392f970fdb42b0e250ce5f5beb9281e1cc0352434de24905fb17137dad7fe` |
| `data/angular/SpherePoints_manifest.toml` | `afc39d4a34b4171bbdaecf52cc4f5dfe0c4f583d1ff4eb8e72252054fb2fedb2` |
| `data/angular/curated_sphere_points.toml` | `4666e89fa2f68397bad8d031817afe460094e7fde94fa1f7cfeb110f48f30f71` |
| `data/legacy/BasisSets` | `6796d34f2c813e5b627b03d003ddd29e0701cfd18e1b4c2ee7d263ce671dbede` |
| `src/foundation/internal/families_high_prec.jl` | `df1585417300676aed7d585cf74071ae4a5568a1d3b076948693ba1b5f72cb87` |

The accepted source change was confined to `src/GaussletBases.jl`,
`src/radial/radial_boundary_prototypes.jl`, `src/angular/angular_point_sets.jl`, and
`src/ordinary/legacy_basis_adapter.jl`: exactly two private declarations and six direct
substitutions produced a `+8/-8` diff. No test, file, data/fixture, include
order, export, API, dependency, metadata framework, cache, numerical policy,
workflow, release state, or later source-layout step changed.

The completed validation resolved all six package-relative targets to the
paths above with exact hash parity, loaded all `122` full angular orders and
curated orders `[15,32,51]`, and loaded the radial prototype, legacy data, and
nested high-precision family source. Unchanged `core`, `radial` (`322/322`),
`angular_public` (`83/83`), `misc` (`59/59`), and docs (`157/157` plus
`10/10`) passed, as did package load, authority checks, Documenter, CI run
`33670968148`, and Docs run `33670968155`. Any future target/byte change, test
edit, move, path framework, elaborate nested-module workaround, or later
source-layout step requires separate authority.

## Radial And Atomic Source Relocation

`HP-SOURCE-LAYOUT-RADIAL-ATOMIC-FN-01` and
`HP-SOURCE-LAYOUT-RADIAL-ATOMIC-TEST-01` maintain the path-only relocation
accepted in commit `812a44d4476c598448bfa3ac733054f009c55b89` from baseline
`5643e3376bd9f8a1a17ded24f05c4d2c363f89cc`. The external follow-up review is
audit evidence rather than authority. Independent inspection fixed its Step 1
accounting at `14` files and `5,469` lines and found two maintained authority
path owners for the radial prototype, not the review's one:
`HP-PACKAGE-ROOT-PATH-FN-01` and
`HP-PUBLIC-RADIAL-PARITY-DOC-FN-01`.

The accepted implementation used `git mv` for exactly this map:

| Pre-move path | Current target | Pre-move SHA-256 |
| --- | --- | --- |
| `src/radial_boundary_prototypes.jl` | `src/radial/radial_boundary_prototypes.jl` | `f843d9c95328f890bfbd276439022b425eb62dfd5b896fc354d9c420d8cc5be7` |
| `src/operators.jl` | `src/radial/operators.jl` | `9e2e954bcb96b4d33434f535b8c12d6a027750e72e412b7c27d0fc383e41df24` |
| `src/radial_ylm_gto_bridge.jl` | `src/radial/radial_ylm_gto_bridge.jl` | `6aa1a3f1533ba0b733a920bd82136cac14a759e5fc6d26fe0442dfa92c5e23b4` |
| `src/atomic_ylm.jl` | `src/atomic/atomic_ylm.jl` | `1a8072648007c08f56987f063fb83e33cf0909ee850a1746dc3a7def0c73ce21` |
| `src/atomic_angular_sectors.jl` | `src/atomic/atomic_angular_sectors.jl` | `df7dc874c4a570a8546d68c54dbc032b6dfa6458fabbae0e202589bbf3b9a703` |
| `src/gaunt_tables.jl` | `src/atomic/gaunt_tables.jl` | `b87f1417175aee60e3615a3ef1972a79dffc1ca64a766f662456256ca6805a79` |
| `src/atomic_ida.jl` | `src/atomic/atomic_ida.jl` | `544159c9f86c375131bbfd61f35dac48b989c54c96ce5d40889b1a3d4b314591` |
| `src/atomic_ida_direct.jl` | `src/atomic/atomic_ida_direct.jl` | `756baf892a32dd8a263b1101880d9f4ebb9ef5dc9db4a0e20441ebea0c293a37` |
| `src/atomic_ida_exchange.jl` | `src/atomic/atomic_ida_exchange.jl` | `fc1883ca05dab6c9a289291be47e7d1a4790ec0ae9abaecfa093de2c1b052b86` |
| `src/atomic_ida_fock.jl` | `src/atomic/atomic_ida_fock.jl` | `d534b50f39630e1842d8d59ab443ce7100fa28a5559378a5ce0fb0529d7df1f4` |
| `src/atomic_ida_uhf.jl` | `src/atomic/atomic_ida_uhf.jl` | `7e428c59fb7415136e9662dd5587aad7e2b64ae0a2fd6421556e9e504986c1c5` |
| `src/atomic_ida_two_electron.jl` | `src/atomic/atomic_ida_two_electron.jl` | `c8ef8ff1e489b96c46c5bb13fe8cb0962b1a38fb0b585fd57d243a94cb9040c8` |
| `src/fullida_dense_export.jl` | `src/atomic/fullida_dense_export.jl` | `46b1b96252ee2c3e07ad1e5f54a467364d6b44ab1be5009846cdc6147657a8d0` |
| `src/sliced_ham_export.jl` | `src/atomic/sliced_ham_export.jl` | `408602299e66e84152b866e76c55beb0312565c4226f591694cfa6a57e59bb12` |

Every moved file remains byte-identical. In `src/GaussletBases.jl`, only the
`14` quoted include paths changed, in their existing positions. The complete
include sequence, including the unmoved files interleaved with these lines,
remains exact. No grouping comment or include reordering was introduced.

Pass 578 removed the only containing-file-location dependencies relevant to
this move. The radial nested module already resolves
`src/foundation/internal/families_high_prec.jl` through the parent-qualified package root;
that include line must move byte-for-byte with its containing file. Independent
inspection found no other `@__DIR__`, nested include, or repository data path
inside the `14` files and no test that directly includes an old path.

The source implementation commit also completed the following exact path
reconciliation so authority and documentation remain valid on the moved tree:

- replace the radial prototype path in both maintained owners named above,
  including the path text in the radial-parity scope;
- replace this relocation record's `14` pre-move `existing` entries with its
  `14` canonical target paths in `existing` state, without changing its
  lifecycle, grant, scope, or evidence. Target paths are not recorded as
  `planned` before the move because the authority schema deliberately rejects
  a missing path whose new parent directory does not yet exist;
- replace the `12` current code pointers across
  `docs/gaunt_backend_note.md`,
  `docs/radial_multipole_stabilization_milestone.md`,
  `docs/radial_paper_parity_prototype_milestone.md`,
  `docs/src/algorithms/atomic_ida_exchange_angular_sectors.md`,
  `docs/src/developer/cartesian_parent_factors_and_cpb_kernels.md`, this
  document's Package-Root and Radial Paper-Parity sections, and
  `sliced_hydrogen_chain.md`;
- regenerate the checked `registry.md` and `execution_whitelist.md` views.
  Historical commit evidence and archived documents may retain their historical
  paths.

The moved-file content delta is zero and the root include substitution is
exactly `+14/-14`, for net-zero production source. No wrapper, module, symbol,
API, dispatch, dependency, cache, fixture, test, workflow, numerical policy,
release change, empty-directory cleanup, or later layout step was added.
Angular, foundation, ordinary, Cartesian, and existing submodule directories
remain where they are.

All `14` target SHA-256 values match the table, and
`git diff --find-renames=100%` reports exactly `14` complete renames. Package
load, full and curated angular-resource loading, and unchanged `core`,
`radial`, `ida`, and `misc` groups passed. Docs tests passed `157/157` plus
`10/10`; Documenter, authority check/self-test, generated-view parity,
manager-log bound, and `git diff --check` passed. CI run `33679026378` and Docs
run `33679026359` passed against the exact implementation commit. Future file
moves, body edits, include reordering, or broader layout work require separate
authority.

## Angular Source Relocation

`HP-SOURCE-LAYOUT-ANGULAR-FN-01` and
`HP-SOURCE-LAYOUT-ANGULAR-TEST-01` maintain the path-only Step 2 relocation
accepted in commit `f8421494e19b3a4553a7f89e20049520c38e3b2c` from
baseline `82b8d97951fa08bd2b295f8d4a9cf6316d6d3b6c`. The external follow-up
review is audit evidence, not authority. The accepted map is:

| Pre-move path | Current target | SHA-256 |
| --- | --- | --- |
| `src/angular_point_sets.jl` | `src/angular/angular_point_sets.jl` | `8b7db8ebcde142212fc72dda1999d013628428a4a8ce7135eb1ed387faac40d9` |
| `src/angular_shell_basis.jl` | `src/angular/angular_shell_basis.jl` | `6799c7ab7ac3c2e9ac6823719348d0143b24b3929e499a326df956c741ea3039` |
| `src/angular_shell_assembly.jl` | `src/angular/angular_shell_assembly.jl` | `3e60d20b2e79a074b7af4db6954428b917d64cb1d705934a5c4448bfa91b548b` |
| `src/angular_atomic_benchmark.jl` | `src/angular/angular_atomic_benchmark.jl` | `d678baf70c4e705112490d63bcdf9aa8865080f5dcc8bf76170bcc0a483917d3` |
| `src/angular_sequence_export.jl` | `src/angular/angular_sequence_export.jl` | `5377c4b047e8d12615504dbd58689c20ff86c1507ee5fdff4e985d112ca6170a` |

All five files remain byte-identical and appear as `100%` renames. Only their
five quoted include paths in `src/GaussletBases.jl` changed, in place, and the
complete include order remains exact. Current authority and documentation
pointers were reconciled mechanically. The accepted tree is
`61aa2190735cbb31267195722a1197ff93ddaebf`; the complete change is
`+31/-31`, with no source-body, test, workflow, API, or numerical edit.

`_lanczos_ground_state_apply` remains defined in
`src/angular/angular_atomic_benchmark.jl` with its exact name, dispatch,
visibility, and load behavior. Its matched-H2+ and angular call sites remain
unchanged; relocation or generalization of that helper requires separate
authority.

Package/resource loading, `angular_public` (`83/83`), the one required complete
angular run (`61,812/61,812`), the isolated HFDMRG adapter (`66/66`), and the
matched H2+ owner (`18/18`) passed. Docs passed `157/157` plus `10/10`;
authority checks, Documenter, and diff checks passed. CI run `33711417479`
passed all three numerical gates and Docs run `33711417480` passed. The full
angular owner is not repeated during lifecycle closeout. Future file-body,
test, workflow, API, numerical, helper-relocation, release, or later-layout
work remains outside this maintenance boundary.

## Foundation Source Relocation

`HP-SOURCE-LAYOUT-FOUNDATION-FN-01` and
`HP-SOURCE-LAYOUT-FOUNDATION-TEST-01` authorize only Step 3 of the external
follow-up review from baseline
`f8421494e19b3a4553a7f89e20049520c38e3b2c`. The review is audit evidence,
not authority. Independent inspection fixes the destination at
`src/foundation/`, matching the architecture's “Broad foundation” layer, and
fixes the inventory at `16` files and `6,001` lines: the fourteen flat
foundation owners plus both tracked files under `src/internal/`.

Thirteen files must move with unchanged bytes and appear as `100%` renames:

| Pre-move path | Authorized target | Required SHA-256 |
| --- | --- | --- |
| `src/mappings.jl` | `src/foundation/mappings.jl` | `605f2dc5b61fd06753e1952756c6fac4881e0fb8cd1c54da2f8fb3999fea9914` |
| `src/stencils.jl` | `src/foundation/stencils.jl` | `109c12436f0edb0a1bd68332bca86bb2b91b07c1bc725b6cb2488bc87e084891` |
| `src/functions.jl` | `src/foundation/functions.jl` | `4c04e03046066c2a6eed32c451380a21998510a5991e2d444446066c4d667995` |
| `src/bases.jl` | `src/foundation/bases.jl` | `442fcb2e25d72692e2c3db99dcca3525e745cac44441468b408512844b58b0ec` |
| `src/quadrature.jl` | `src/foundation/quadrature.jl` | `b8cb02661521abe60f9026320a65e08809d12e5436deae20ed7b7b4fab61ce80` |
| `src/primitive_sets.jl` | `src/foundation/primitive_sets.jl` | `edfd611468516c1707e26a6ecc954e59cd3860b69c160a0eeb9ae09c599a2a99` |
| `src/partitions.jl` | `src/foundation/partitions.jl` | `8fef22f9ee424b9b12cd0e7f613a6b7ff2479b334244cef4d0749203d69861c4` |
| `src/hierarchical_partitions.jl` | `src/foundation/hierarchical_partitions.jl` | `e18b6f1fc3594e9ea47a65db375c96181cdb347e6cfaa9265f7c4f012e7bb183` |
| `src/leaf_pgdg.jl` | `src/foundation/leaf_pgdg.jl` | `bb48f1492ebf32b228b5456e7792e938ef610f013a28f5a72c085189605337cc` |
| `src/global_leaf_contraction.jl` | `src/foundation/global_leaf_contraction.jl` | `56510159720aaa08b6d8d9891574d45fa573cd8de920d4594c3ddff196eb1004` |
| `src/diagnostics.jl` | `src/foundation/diagnostics.jl` | `69fd53ad714765663a5009c7d9fbe714b36f596b36e501a0a77949ba683f75d7` |
| `src/timing.jl` | `src/foundation/timing.jl` | `6edc25510ad3d3819912583b575d23d3e7b028db583582b3db2650cb1db684c4` |
| `src/GaussianAnalyticIntegrals.jl` | `src/foundation/GaussianAnalyticIntegrals.jl` | `b593aff9394ea49d380cc071b1c7b5b8eb3e8165ba07f22e67e5d8b0acc3d304` |

Three moved files contain four current path comments. Leaving those comments
stale would violate the path-only purpose, so only these literal substitutions
are authorized; no executable statement may change:

| Pre-move path | Authorized target | Baseline SHA-256 | Required target SHA-256 |
| --- | --- | --- | --- |
| `src/families.jl` | `src/foundation/families.jl` | `29dc311da6cfc7bb280dc7fc429e3cd5de2af1072c249509a3c522978a58fd20` | `c1550b78b597c734a13a398a790f56009349559019e18b89d39652e709e2738d` |
| `src/internal/families_high_prec.jl` | `src/foundation/internal/families_high_prec.jl` | `4c7ee5dee90cfb3296b392a41ba3cbd30423bd30aaea84935bb60ebd23608877` | `df1585417300676aed7d585cf74071ae4a5568a1d3b076948693ba1b5f72cb87` |
| `src/internal/wavelet_filters.jl` | `src/foundation/internal/wavelet_filters.jl` | `6162f7328ae1badb742f292fdec8b0bb9d2f7233642103251e47152d5be11cfe` | `34f006a0f20f0f1ade3c7efb26cf86e88c1751ac19ce63de8af5be10b7939043` |

The exact four comment changes are the two `src/internal/...` pointers in
`families.jl` and the one `src/families.jl` pointer in each internal file.
`GaussianAnalyticIntegrals` and `TimeG` remain nested modules loaded at their
existing positions. The active call from `primitive_sets.jl` to the later
ordinary mapped-backend helpers remains unchanged; this transaction neither
relocates nor repairs that cross-track dependency.

In `src/GaussletBases.jl`, only the fifteen quoted include paths for these
runtime-loaded files may change, in place. The complete include sequence and
the adjacent `using .TimeG` block must remain exact. The radial prototype's
parent-qualified include may change only from `src/internal/families_high_prec.jl`
to `src/foundation/internal/families_high_prec.jl`; its required post-change
SHA-256 is `fd00b5ee833a15b03f394db320b76e36e9cf2b73a4359fedd59be8852e57397d`.
The same source-path substitution is authorized in `test/runtests.jl` and
`test/radial/runtests.jl`, with no assertion or runner change.

Current-pointer reconciliation is mechanical: replace the `26` affected path
tokens across the nine hand-edited current documents named by the authority
record; replace the three active authority paths for mappings and the shared
Gaussian analytic owner; and regenerate the corresponding registry entries and
both whole-file views, `registry.md` and `execution_whitelist.md`. Also reconcile
the Package-Root target path/hash and the Radial/Atomic maintenance wording for
the one changed parent-qualified include. Historical evidence, pre-move mapping
columns, and completed validation facts retain their historical paths and
hashes.

Implementation target card:

```text
Target: relocate exactly 16 foundation owners into src/foundation/
Physics endpoint: none changed; package load plus core/radial are the local oracles
Allowed files: 16 git mv operations, four path comments, 15 root includes, one radial include, two test paths, and exact current-pointer reconciliation
Forbidden files/surfaces: all other source bodies, tests, workflows, APIs, numerics, releases, and later layout work
Must delete or simplify: remove the flat foundation locations and tracked src/internal location; delete no implementation
Forbidden additions: wrappers, modules, symbols, aliases, helpers, caches, fixtures, metadata, dependencies, tests, or files beyond destination directories
Success condition: 13 byte-identical R100 moves, three exact comment-only moves, exact path substitutions, and unchanged behavior
Validation: hashes/diffs, include order, package-root targets, package load, core/radial, docs/authority/Documenter/diff, normal CI/Docs
Line budget: source path substitutions +20/-20; tests +2/-2; current documentation and authority path substitutions only
Failure rule: make no implementation commit if any hash, executable line, include position, package target, or unchanged test differs
```

The complete angular owner is not part of Step 3 validation and must not be
repeated. No module wrapper, symbol move, include reorder, API, dispatch,
dependency, fixture, cache, numerical policy, workflow, release work, Step 4,
or later source-layout authority is granted.

## Ordinary Source Relocation

`HP-SOURCE-LAYOUT-ORDINARY-FN-01` maintains the Ordinary half of source-layout
Step 4 accepted in commit `96802b00cc380707306106b529ff22a09ccb415e`
from baseline `24512b0beda75174572fc2b728010b432331c05c`;
`HP-SOURCE-LAYOUT-ORDINARY-TEST-01` is completed with no implementation grant.
The external follow-up review is audit evidence, not authority. Independent
inspection split Step 4 because the Ordinary and Cartesian path closures are
substantially separable. This transaction owns exactly `15` flat files and
`14,430` lines.

The implementation must use `git mv` for exactly this frozen map:

| Pre-move path | Authorized target | Git blob | SHA-256 |
| --- | --- | --- | --- |
| `src/legacy_basis_adapter.jl` | `src/ordinary/legacy_basis_adapter.jl` | `b630ce9544297b76879a0013a984e9b64ba7ea32` | `5469b107c1a270cdfb169cd621beb656c0d8485b432ea7a1a68928d5bcd890f6` |
| `src/ordinary_cartesian_ida.jl` | `src/ordinary/ordinary_cartesian_ida.jl` | `4f85f37c92108626ff99c5971a196202bae22c14` | `18fa23fb09c3034d1791005e1f3e514721ed639c700718ef0f301265a56902fe` |
| `src/ordinary_coulomb.jl` | `src/ordinary/ordinary_coulomb.jl` | `94fd3f0d1c5f671316365c7264ba179a703ccf88` | `469efa553ff3a44913a4ea00170f0a122b0f920c412c3061290dd9a3baf44183` |
| `src/ordinary_hybrid.jl` | `src/ordinary/ordinary_hybrid.jl` | `32bfc7230711a8a8bf7dabaf9e7173a12b50e4a3` | `154ac16119779fcbf2ebb0f22b1873cbad8268773bafaeef765019e25c22ac77` |
| `src/ordinary_mapped_backends.jl` | `src/ordinary/ordinary_mapped_backends.jl` | `b17a84858a60405c1491b01d18a7188bbb993b57` | `de450f2142e37e10cf5951e2dff5a41812c43935274c04c55c91a5b02736df9d` |
| `src/ordinary_pgdg.jl` | `src/ordinary/ordinary_pgdg.jl` | `97032896bf29b3fafc09c1ba8a5742c2f4520a99` | `a25d23e4077e9ca9656e397308315b6692f913450a42121793fc56b55f55f149` |
| `src/ordinary_pgdg_refinement_masks.jl` | `src/ordinary/ordinary_pgdg_refinement_masks.jl` | `845d9b74f0541a7e9b7d4997e9004e78f94bdea3` | `326348f60b31f827339bec36c4ba9ca12cad3e85366dea22e9184ccfb91cef05` |
| `src/ordinary_qw_corrections.jl` | `src/ordinary/ordinary_qw_corrections.jl` | `046196521cdd70e7cb9621919f067dd92c9e9372` | `d07cb98640861dbbeb7b0d3f3e5852e249bed747c8749bbfe2a1994d8fdbd272` |
| `src/ordinary_qw_nested_frontends.jl` | `src/ordinary/ordinary_qw_nested_frontends.jl` | `7ca32b580aadffc0079ee83e0207b777bbd948ac` | `0bb862f36d8ef9bf6e080fee8e8cf9ce95e99bdb5b49e7016610dceb5f8e3fd0` |
| `src/ordinary_qw_operator_assembly.jl` | `src/ordinary/ordinary_qw_operator_assembly.jl` | `d321d07d187b452ff6b2ddda42f782b297f2a8af` | `57eb1411f07e99383b66687172a5c12d835121c0cb249d01651539edbede968f` |
| `src/ordinary_qw_raw_blocks.jl` | `src/ordinary/ordinary_qw_raw_blocks.jl` | `bb2aacde6adbd26993618b989f24ea59fcf32093` | `3c724e76a44866406d6d20a8a4ce60388fd9415c2155ba7480273df6c36b79da` |
| `src/ordinary_qw_residuals.jl` | `src/ordinary/ordinary_qw_residuals.jl` | `19c7f081429370d5b9cf51b592649eccfd82ff5e` | `1ac7e69771ddb101e7f55f30c3e153890c4ffad749534a90ef362a77a440061a` |
| `src/ordinary_qw_types_and_bases.jl` | `src/ordinary/ordinary_qw_types_and_bases.jl` | `7dbca0e51e85c3a884b8862a4ac19521e8a2dd03` | `333183ed6957b16c1fcf2a14dde832b916e69bf96ccf673c63c3b20a7bac5a97` |
| `src/ordinary_sho.jl` | `src/ordinary/ordinary_sho.jl` | `9a5e3fd844e82b5a914cbca4354205ee37819e59` | `0218364dc8af4a62ad16f84515f71de2fbf24d4e7f9b134c402b2e8415d47ef5` |
| `src/sliced_hydrogen_chain.jl` | `src/ordinary/sliced_hydrogen_chain.jl` | `462cbdfaabd53ccff9558095cb376b6cbcb8b2f5` | `186044e0bdd951a9db7380e0c364fb96a1d7b2981357c601df7962429649301a` |

Every moved file must remain byte-identical and appear as a `100%` rename. In
`src/GaussletBases.jl`, only the corresponding `15` quoted include paths may
change, in their existing positions. The complete include order remains exact;
no grouping comments are added.

The frozen owners contain no nested include and no source-location-dependent
data path. Vendored legacy data already uses `_package_data_path`; its existing
explicit caller path and home-directory fallback are behavior, not relocation
dependencies. Outside the root `src/GaussletBases.jl` loader, no source, test,
example, script, tool, or workflow directly includes any of these files. The
known call from
`src/foundation/primitive_sets.jl` into later ordinary mapped-backend helpers
remains unchanged and valid under the existing root include sequence.

Current-pointer reconciliation is mechanical and part of the same source
commit:

- replace `14` existing path entries across the ten pre-existing authority
  records and five path literals in three of their scope strings;
- replace this relocation record's `15` source paths with their targets while
  preserving lifecycle, grant, scope, evidence, and all nonpath text;
- replace `86` current code-pointer occurrences across the `17` hand-edited
  documents named by this record; and
- update the canonical document digest and regenerate `registry.md` and
  `execution_whitelist.md` as whole files.

The explicit pre-move-path column above, archived reviews, completed historical
evidence, and the `65` classified historical occurrences retain their old path
text. The implementation must not perform a repository-wide blind replacement.

Implementation target card:

```text
Target: relocate exactly 15 Ordinary owners into src/ordinary/
Physics endpoint: none changed; the existing public and ordinary owners remain the oracles
Allowed files: 15 git mv operations, 15 root include substitutions, and exact current authority/documentation path reconciliation
Forbidden files/surfaces: every file body, test, workflow, API, numerical contract, release, and Cartesian/PQS move
Must delete or simplify: remove the 15 flat source locations; delete no implementation
Forbidden additions: modules, namespaces, wrappers, aliases, shims, helpers, caches, fixtures, dependencies, tests, or files beyond the destination directory
Success condition: 15 byte-identical R100 moves, exact path substitutions, unchanged include order and behavior
Validation: hashes/renames, include order and direct-path closure, package load, Supported floor, docs/authority/Documenter/diff, normal CI/Docs
Line budget: moved bodies +0/-0; src/GaussletBases.jl +15/-15; current authority/documentation path substitutions only; tests +0/-0
Failure rule: make no implementation commit if any hash, rename score, include position, body, test, path count, package load, or unchanged owner differs
```

The unchanged Supported-floor selection must pass, including `core`, `ida`,
`cartesian`, `examples`, `radial`, `misc`, `angular_public`, and `docs_fast`.
The source-bearing implementation commit must classify as full and pass all
three existing CI jobs plus Docs. No complete angular rerun is required.
Cartesian/PQS Step 4B is governed separately below. Decomposition, Lanczos
relocation, foundation dependency inversion, orphan-test work, conditioning
diagnostics, prose-test cleanup, empty-directory cleanup, API/export changes,
and release work remain outside this maintenance authority.

The independently governed
[Cartesian/PQS source relocation](cartesian_pqs_source_relocation.md) follows
the accepted Ordinary move and changes no public or numerical contract.

## Foundational Basis And Mapping Documentation

`HP-PUBLIC-FOUNDATION-DOC-FN-01` and
`HP-PUBLIC-FOUNDATION-DOC-TEST-01` maintain the implemented documentation for
ten existing exported generic functions. Commit
`03eb75c66e5cc8b24714d3769097ea81f0d74b15` added no binding or behavior. The
maintained meanings are:

| Binding | Public meaning |
| --- | --- |
| `uofx` | Map physical coordinate `x` to reference coordinate `u`. |
| `xofu` | Apply the inverse map from reference `u` to physical `x`. |
| `dudx` | Evaluate the first derivative of `u(x)` with respect to `x`. |
| `du2dx2` | Evaluate the second derivative of `u(x)` with respect to `x`. |
| `basis_spec` | Return the construction specification associated with a basis-like object. |
| `family` | Identify the object's gausslet family. |
| `mapping` | Return its physical/reference coordinate map. |
| `centers` | Return physical-coordinate basis centers. |
| `reference_centers` | Return reference-coordinate basis centers. |
| `integral_weights` | Return basis-function integrals, not quadrature weights. |

The implementation added `58` source-docstring lines, `19` reference lines,
and `8` docs-test lines, all within the preferred limits. The existing
evaluation/introspection section contains `@docs` entries for exactly these
ten bindings. Family-scoped checks require every binding to remain exported,
documented, and present in that page. The undocumented exported-binding
backlog fell exactly from `71` to `61`; no global allowlist or `checkdocs`
change was added.

Focused radial tests passed `322/322`, docs tests passed `118/118` plus
`10/10`, and package load, authority check/self-test, Documenter, and diff
checks passed. Remote CI run `33181875262` and Docs run `33181875247` passed.
Both records are in maintenance. Preserve every definition, method, signature,
dispatch, export, return value, default, numerical result, compatibility floor,
workflow, tag, and release artifact. Do not broaden the page into evaluation
algorithms, stencils, partitions, atomic operators, angular work, further API
reduction, dependencies, examples, or release work.

## Function Evaluation And Stencil Documentation

`HP-PUBLIC-FUNCTION-STENCIL-DOC-FN-01` and
`HP-PUBLIC-FUNCTION-STENCIL-DOC-TEST-01` maintain the implemented documentation
for nine existing exported generic functions. Commit
`624718972627b5ebf74beaa57ad73b8725aac587` added no binding or behavior. The
maintained meanings are:

| Binding | Public meaning |
| --- | --- |
| `value` | Normal function evaluation; calling `f(x)` delegates to this route. |
| `direct_value` | Evaluate the object's stencil representation exactly as represented, not as an independent analytic oracle. |
| `derivative` | Evaluate the requested nonnegative derivative order with respect to physical `x`. |
| `center` | Return the function's physical center. |
| `reference_center` | Return the corresponding pre-mapping/reference-coordinate center. |
| `integral_weight` | Return one function's integral, distinct from basis-level `integral_weights` and quadrature weights. |
| `stencil` | Expose the ordered primitive expansion. |
| `coefficients` | Return stored coefficient data; for `FunctionStencil`, order matches `primitives(stencil)`. |
| `terms` | Materialize ordered coefficient/primitive pairs as `StencilTerm` values. |

The implementation added `57` source-docstring lines, `18` reference lines,
and `8` docs-test lines, all within the preferred limits. The existing
function-evaluation/stencil section contains `@docs` entries for exactly these
nine bindings. Family-scoped checks require every binding to remain exported,
documented, and present in that page. The undocumented exported-binding
backlog fell exactly from `61` to `52`; no global allowlist or `checkdocs`
change was added.

Focused core/radial tests and radial `322/322` passed; docs tests passed
`121/121` plus `10/10`. Package load, authority check/self-test, Documenter,
and diff checks passed. Remote CI run `33213490510` and Docs run `33213490543`
passed. Both records are in maintenance. Preserve every definition, method,
signature, dispatch, export, numerical result, default, allocation policy,
compatibility floor, workflow, tag, and release artifact. Do not broaden the
page into basis metadata, partitions, operators, angular work, further API
reduction, dependencies, examples, or release work.

## Partition Hierarchy And Leaf-Local Accessor Documentation

`HP-PUBLIC-PARTITION-LEAF-DOC-FN-01` and
`HP-PUBLIC-PARTITION-LEAF-DOC-TEST-01` maintain the implemented documentation
for twelve existing exported partition-hierarchy and leaf-local accessors.
Commit `b24c0dbeeb0a6b5526fde9004cbc7319c1395ec4` added no binding or behavior.
The maintained meanings are:

| Binding | Public meaning |
| --- | --- |
| `boxes` | Return stored ordered box records for a flat or hierarchical partition. |
| `leaf_boxes` | Select hierarchy nodes with no children while retaining stored order. |
| `box_indices` | Return basis indices owned by the selected box. |
| `box_level` | Return the selected hierarchy node's level. |
| `box_parent` | Return its parent identifier, with `nothing` for a root. |
| `box_children` | Return its child identifiers. |
| `box_block` | Materialize a `Float64` copy of one diagonal box block. |
| `box_coupling` | Materialize a `Float64` copy of the rectangular block between two boxes. |
| `leaf_primitive_indices` | Map a leaf-box identifier to its global primitive-index range. |
| `primitive_origins` | Return primitive-aligned provenance labels. |
| `primitive_leaf_boxes` | Return primitive-aligned leaf-box provenance. |
| `leaf_contractions` | Return stored ordered `LeafBoxContraction1D` records. |

The block accessors accept either a square matrix directly or a named operator
from `BasisRepresentation1D`. Accessors exposing stored vectors or records must
tell callers to treat those objects as read-only; this grant does not change
them into copies or add mutation protection.

The compact implemented `Partitions and leaf-local layers` section in
`docs/src/reference/bases_and_mappings.md` includes those accessors and this
already-documented public context without duplicating prose:

- `BasisBox1D`, `BasisPartition1D`, and `basis_partition`;
- `HierarchicalBasisBox1D`, `HierarchicalBasisPartition1D`,
  `hierarchical_partition`, and `refine_partition`;
- `LeafGaussianSpec1D`, `LeafLocalPGDG1D`, `build_leaf_pgdg`, and
  `augment_leaf_pgdg`;
- `GlobalMappedPrimitiveLayer1D` and `build_global_mapped_primitive_layer`;
- `LeafBoxContraction1D`, `LeafBoxContractionLayer1D`, and
  `contract_leaf_boxes`.

The implementation added `80` source-docstring lines, `35` reference lines,
and `15` docs-test lines, all within the authorized hard limits and with no new
file. Family-scoped checks require all twelve accessors to remain exported,
documented, and present with the listed context. The undocumented exported-
binding backlog, excluding the module self-binding, fell exactly from `52` to
`40`; no global allowlist or `checkdocs` change was added.

The focused core suite and docs tests `123/123` plus `10/10` passed. Package
load, authority check/self-test, Documenter, and diff checks passed; remote CI
run `33219963234` and Docs run `33219963250` passed. Both records are in
maintenance. Preserve every definition, method, signature, dispatch, export,
array ownership, ordering, numerical result, default, compatibility floor,
workflow, tag, and release artifact. Do not change returned objects into
copies, add mutation protection, alter indexing, or broaden into basis
metadata, operators, angular work, API reduction, dependencies, examples, or
release work.

## Atomic IDA Inspection And Tiny Two-Electron Reference Documentation

`HP-PUBLIC-ATOMIC-IDA-DOC-FN-01` and
`HP-PUBLIC-ATOMIC-IDA-DOC-TEST-01` maintain documentation for ten existing
exported inspection and tiny-reference functions, implemented by commit
`8f84afb491ff7cc0f818a6ed982eca638a509b65`:

| Binding | Required public meaning |
| --- | --- |
| `orbitals` | Return an owning object's stored ordered orbital records; callers treat returned storage as read-only, and ordering is only the owning construction's documented ordering. |
| `two_electron_states` | Return the tiny problem's stored ordered state records under the same read-only and owner-specific-order rule. |
| `radial_multipole` | Access the requested radial multipole matrix from `AtomicIDAOperators`. |
| `gaunt_coefficient` | Return one complex spherical-harmonic Gaunt coefficient. |
| `gaunt_tensor` | Reconstruct the requested dense Gaunt tensor on demand from the internal sparse representation. |
| `angular_kernel` | Reconstruct the requested dense angular kernel on demand from the internal sectorized representation. |
| `apply_overlap` | Validate the coefficient dimension and apply the tiny problem's stored dense overlap matrix. |
| `apply_hamiltonian` | Validate the coefficient dimension and apply its stored dense Hamiltonian matrix. |
| `ground_state_energy` | Return the lowest energy from the stored dense Hermitian problem. |
| `lanczos_ground_state` | Run the existing small fully reorthogonalized reference Lanczos routine with its current controls and return `value`, `vector`, `residual`, `iterations`, and `converged`. |

The compact section in `docs/src/reference/atomic_and_ordinary.md` places these
functions with the already-documented `AtomicOrbital`,
`AtomicIDAOperators`, `atomic_ida_operators`, `AtomicIDATwoElectronState`,
`AtomicIDATwoElectronProblem`, and `atomic_ida_two_electron_problem` context.
It does not duplicate explanatory prose. `gaunt_tensor` and `angular_kernel`
promise neither retained dense storage nor caching or scalable large-angular-
momentum use. `lanczos_ground_state` is a small reference routine, not a
general many-electron or production eigensolver, and
`AtomicIDATwoElectronProblem` remains only the tiny one-up/one-down consumer.

The implementation added `70` source-docstring lines, `26` reference lines,
and `15` focused docs-test lines. These additions are within all hard limits,
with no new file or executable change. Family-scoped checks require the exact
ten names to remain exported, documented, and present in the intended section,
and preserve the read-only/order, on-demand reconstruction, dense-reference,
and non-production limits. Documenter remains the executable `@docs`
resolution gate. The undocumented exported-binding backlog, excluding the
module self-binding, fell exactly from `40` to `30`.

The required classification stop was refined by the bounded documentation and
read-only caller audits that followed it. The six supported-public candidates
are now documented under the next section. At snapshot
`0832bdce0f3ff9bdd07d8580bedcda0c47dadbe8`, the remaining audit retains
`bond_aligned_diatomic_geometry_payload` as expert/experimental and reserves
four original de-export candidates in two next-minor transactions described
below. The already-documented `QWRGResidualSpaceDiagnostics` type is coupled
to the residual-diagnostics transaction. No name is deletion-ready, and this
classification grants no export change.

Docs tests passed `126/126` plus `10/10`; the focused core group, package load,
authority check/self-test, Documenter, manager-log bound, and diff checks also
passed. Remote CI run `33230195333` and Docs run `33230195291` passed. Both
records are in maintenance. The remaining classification is not an automatic
documentation-debt queue.

Preserve every implementation, method, signature, dispatch, export, numerical
behavior, solver policy, storage and cache choice, compatibility floor,
workflow, and release artifact. Add no generalized solver, dense-tensor cache,
global `checkdocs` policy or allowlist, ordinary-test change, angular-research
work, source behavior, dependency, example, API reduction, or release work.

## Experimental Angular Producer Documentation

`HP-PUBLIC-ANGULAR-PRODUCER-DOC-FN-01` and
`HP-PUBLIC-ANGULAR-PRODUCER-DOC-TEST-01` maintain the accepted documentation
for exactly these seven durable experimental producer bindings, implemented by
commit `9b410ea14356cee6085ad2de51e73849695a0d97`:

- `ShellLocalAngularProfile`;
- `ShellLocalAngularProfileOverlap`;
- `shell_local_angular_profile`;
- `adjacent_shell_local_angular_profile_overlap`;
- `AtomicFixedRadialAngularSequenceLevel`;
- `AtomicFixedRadialAngularSequenceOverlapSidecar`; and
- `AtomicFixedRadialAngularSequence`.

This remains an experimental producer-side surface, not a completed angular
application workflow. A shell-local profile records an ordered exact injected-
Ylm block followed by the mixed complement, together with deterministic labels
and gauge metadata. `profile_id`, level IDs, and sequence IDs are provenance
and integrity identities; they are not tolerance-, convention-, or numerical-
equivalence claims. The cache's `ShellLocalAngularProfileKey` is not a consumer
identity and must not be promoted by this packet.

A `ShellLocalAngularProfileOverlap` describes the source-to-target shell-local
basis overlap. It is shell-independent. Fixed-radial sequence sidecars retain
the distinction between adjacent level pairs and the complete non-adjacent
direct upper triangle. Every sequence level shares one radial basis, shell
inventory, and shell centers while changing `N_sph` and its angular profile.
The final working bases are intended to be orthonormal: these overlap sidecars
are continuation and transfer data, not authority for generalized-overlap
reinterpretation. Reader documentation must direct consumers to the documented
IDs, labels, gauge metadata, diagnostics, and existing payload access rather
than to the profile cache key.

Maintenance is limited to the accepted concise docstrings in
`src/angular/angular_shell_basis.jl` and `src/angular/angular_sequence_export.jl`, one compact
reader-facing addition in `docs/src/explanations/angular_research_track.md`,
one curated seven-binding section in `docs/src/reference/export.md`, and
focused checks in `test/docs/runtests.jl`. The accepted implementation added
`58` source-docstring lines, changed reader/reference documentation by
`+42/-10`, and changed the docs tests by `+18/-1`, with no new file.

All seven names remain exported, documented, and in the curated reference.
`ShellLocalAngularProfileKey` remains exported temporarily but has neither a
docstring nor a reference entry. The undocumented exported-binding count,
excluding the module self-binding, is exactly `16`, reduced from `23`. The key
remains the sole angular next-minor de-export candidate; no de-export is
authorized here. The focused `angular_public` owner passed `83/83`; the roughly
`930`-second full angular suite was not repeated and is not a routine
maintenance gate for this documentation.

Restart orchestration, common-target lifts, Givens transformations, downstream
campaigns, cache API, definition or field changes, signatures, dispatch,
exports, cache behavior, numerical policy or results, workflows, dependencies,
versions, releases, and key de-export remain outside maintenance authority.

## Radial Paper-Parity Reproducibility Documentation

`HP-PUBLIC-RADIAL-PARITY-DOC-FN-01` and
`HP-PUBLIC-RADIAL-PARITY-DOC-TEST-01` maintain the accepted documentation for
exactly these four exported bindings:

- `RadialBoundaryPrototype`;
- `radial_boundary_prototype`;
- `radial_boundary_prototype_names`; and
- `build_paper_parity_radial_basis`.

This is an expert reproducibility and provenance surface, not the recommended
beginner construction path. `RadialBoundaryPrototype` is a validated,
artifact-backed, cached description of a frozen radial boundary construction.
The only currently supported name is `:paper_parity_g10_k6_x2`. Prototype
records and their contained arrays are shared cached data and must be treated
as read-only inspection data. Detailed high-precision coefficients, primitive
storage, diagnostics, checksums, and provenance fields do not form a general
mutable construction API.

The frozen prototype uses exact widths `0.09358986806` and `0.02357750369`;
these are distinct from the rounded normal-workflow presets `0.0936` and
`0.0236`. `radial_boundary_prototype_names()` reports the available frozen
prototypes but is not a general registry or plugin mechanism.
`radial_boundary_prototype(name)` returns the validated cached prototype and
rejects unknown names.

`build_paper_parity_radial_basis` holds the frozen boundary construction fixed
while extending it to the caller's requested `rmax` and mapping. Its default
count policy is `rmax_count_policy=:legacy_strict_trim`. Reader documentation
must also explain the corresponding `build_basis(prototype; ...)` route,
including the frozen mapped prototype when `rmax` is omitted and the same
paper-parity extension when `rmax` is supplied. These routes support controlled
paper-parity comparisons; they neither promise exact reproduction of an entire
paper calculation nor replace the normal `RadialBasisSpec` front door.

The accepted implementation in commit `84ebdc5792c954b4c09f08b53ba86ef9332bf6b8`
is limited to concise docstrings in
`src/radial/radial_boundary_prototypes.jl`, one compact expert/paper-parity section in
`docs/src/reference/bases_and_mappings.md`, a short distinction from the normal
recommended route in `docs/src/howto/recommended_atomic_setup.md`, and focused
checks in `test/docs/runtests.jl`. It added `33` source-docstring lines, `32`
combined reader/reference lines, and `14/-1` docs-test lines, with no new file
or executable change.

All four names remain exported, have source docstrings, and appear together in
the curated expert/paper-parity reference section. The recommended setup page
retains `RadialBasisSpec` as the ordinary front door and explains the
exact-versus-rounded width distinction. The undocumented exported-binding
count, excluding the module self-binding, is exactly `12` after the accepted
`16 -> 12` reduction. Existing radial tests and the unchanged Supported-floor
gate remain the numerical maintenance owners.

Coefficients, artifacts, definitions, fields, methods, signatures, dispatch,
defaults, cache behavior, numerical results, exports, dependencies, workflows,
versions, and releases remain outside this packet. The remaining expert
families and every next-minor de-export candidate are excluded from maintenance.

## Experimental QW Geometry Diagnostics Documentation

`HP-PUBLIC-QW-GEOMETRY-DOC-FN-01` and
`HP-PUBLIC-QW-GEOMETRY-DOC-TEST-01` maintain the accepted documentation from
commit `c553c40c2a464f639c4385b6155d37963facdc3d` for exactly these three
exported route-specific expert diagnostics:

- `bond_aligned_diatomic_nested_geometry_diagnostics`;
- `bond_aligned_homonuclear_chain_geometry_diagnostics`; and
- `axis_aligned_homonuclear_square_lattice_geometry_diagnostics`.

These functions inspect the current QW and nested geometry conventions. They
are not a general molecular-geometry API or a stable serialization schema.
Returned vectors and route-specific named records are inspection data and must
be treated as read-only.

The nested-diatomic source overload inspects and preserves the exact supplied
nested source; it must not rebuild or normalize another source. Its basis
overload constructs the normalized nested source under the documented frontend
controls and is therefore not a cheap geometry-only query. Documentation must
identify its principal results: geometry, child and shared-shell retention
contracts, shell counts and dimensions, actual shared-shell provenance, fixed
dimension, contract audit, and atom-growth anatomy. The existing forwarding of
shared-shell policy, endcap/panel `q` and `L`, packet kernel, retention
overrides, split guards, and protection controls remains part of the signature
contract. Reader prose need not duplicate every default already visible there.

The chain diagnostic only inspects an existing experimental homonuclear chain
basis. It reports chain coordinates, mapped-axis and transverse centers, local
spacings at nuclei and midpoints, mapping kinds, monotonicity, and axis-center
symmetry error. The square-lattice diagnostic only inspects an existing
axis-aligned homonuclear square basis. It reports x/y coordinates, Cartesian
basis centers, nuclear-coordinate, representative-midpoint, and plane-center
spacings, mapping kinds, monotonicity and symmetry errors, and x-y center
agreement. Neither route constructs a nested source, Hamiltonian, solver input,
arbitrary-orientation model, or heteronuclear model.

The accepted implementation is limited to concise docstrings in
`src/ordinary/ordinary_qw_nested_frontends.jl` and
`src/ordinary/ordinary_qw_types_and_bases.jl`, one compact section titled
"Experimental QW geometry diagnostics" in
`docs/src/reference/operators_and_diagnostics.md`, and focused checks in
`test/docs/runtests.jl`. It added `31` source-docstring and `24` reader-
documentation lines and changed the docs tests by `+13/-1`, with no new file
or executable change.

All three names remain exported, carry Julia docstrings, and appear together
in the curated reference. The checks distinguish exact source reuse from
basis-overload source construction and lightweight chain/square inspection
from nested-diatomic construction. The undocumented exported-binding count,
excluding the module self-binding, is exactly `9` after the accepted
`12 -> 9` reduction. Relevant core geometry testsets passed `6/6`, `14/14`,
`14/14`, and `35/35`; docs passed `151/151` plus `10/10`. Remote CI run
`33545468306` passed all three numerical gates, and Docs run `33545468193`
passed and deployed the verified live reference.

No returned field, definition, signature, default, method, dispatch,
construction behavior, export, ordinary numerical test, workflow, dependency,
numerical policy, version, tag, or release may change. This packet excludes
`cartesian_base_working_basis`, sliced-chain operators, every next-minor
de-export candidate, global documentation policy, and broader QW work. If
maintenance requires any excluded change, it needs a separate authority
amendment.

## Experimental Sliced Hydrogen-Chain Operator Documentation

`HP-PUBLIC-SLICED-HCHAIN-DOC-FN-01` and
`HP-PUBLIC-SLICED-HCHAIN-DOC-TEST-01` maintain the accepted documentation from
commit `ec3b8fa12a68118c922ae87ec3cacae821633d5b` for exactly these three
exported operator functions:

- `sliced_h1`;
- `sliced_vee`; and
- `sliced_row!`.

The compact expert-reference section must present them together with the
already-documented companion exports `SlicedHydrogenChain`,
`sliced_hydrogen_chain`, and `sliced_h1_bandwidth`. This is one experimental
minimal hydrogen-chain producer, not a general molecular Hamiltonian, solver,
MPO, electron-sector framework, or sliced-basis framework.

`sliced_h1(chain)` returns a read-only, lazily evaluated
`AbstractMatrix{Float64}` view of the electronic one-body operator. Nuclear
repulsion remains the separate `chain.nuclear_repulsion` scalar. The
represented structural half-bandwidth is returned by
`sliced_h1_bandwidth(chain)`. Entries beyond that band are exactly zero in the
accepted representation; this is not a claim that omitted continuum matrix
elements vanish physically.

`sliced_vee(chain)` returns a read-only, lazily evaluated two-index
density-density interaction view, not a four-index ERI tensor. Both operator
views retain their chain owner and expose matrix size, scalar indexing, and
`mul!` without owning an `N x N` dense matrix. Their concrete view types remain
private and non-contractual; consumers obtain them only through `sliced_h1`
and `sliced_vee`. Explicit `Matrix(view)` materializes a dense matrix and is
suitable only for bounded diagnostics.

`sliced_row!(destination, view, row)` validates destination length and row
bounds, completely overwrites the correctly sized caller-owned destination,
and returns that destination. Preserve the existing zero-allocation path for a
compatible caller buffer. H1 row extraction may use its represented stored
band; Vee row extraction covers every column.

Existing stable chain fields and diagnostics remain read-only consumer data.
Internal H1 bands, Coulomb blocks, longitudinal coefficients, transverse work
arrays, and operator-view concrete types remain non-contractual. Documentation
must not turn those internals into a storage, ownership, or mutation promise.

The accepted implementation is limited to concise docstrings in
`src/ordinary/sliced_hydrogen_chain.jl`, one compact section in
`docs/src/reference/atomic_and_ordinary.md`, and focused checks in
`test/docs/runtests.jl`. It added `36` source-docstring and `32`
reader-reference lines and changed docs tests by `+12/-1`, with no new file or
executable change.

All six exports now appear as one coherent expert interface, the three
previously undocumented functions have docstrings, and the undocumented
exported-binding count, excluding the module self-binding, is exactly `6`
after the accepted `9 -> 6` reduction. Existing sliced-chain numerical,
allocation, bandwidth, finite-chain, periodic-template, singleton, and
long-range behavior remains unchanged. The minimal sliced hydrogen-chain owner
passed `72/72`, docs passed `154/154` plus `10/10`, and package load,
authority check/self-test,
Documenter, and diff checks passed. Remote CI run `33565530332` passed all three
gates; Docs run `33565530379` passed and deployed the verified live reference.
Existing IDA validation and Supported-floor CI remain the numerical owners; no
new numerical assertion is authorized.

No type, field, definition, method, signature, dispatch, coefficient,
tolerance, storage, allocation policy, export, workflow, dependency, version,
tag, or release may change. This packet excludes `cartesian_base_working_basis`,
every next-minor de-export candidate, API reduction, global documentation
policy, new examples or fixtures, and broader sliced-chain work. If truthful
documentation requires an excluded change or exceeds a hard budget, make no
implementation commit and report the exact mismatch.

## Expert Staged Cartesian Construction Documentation

`HP-PUBLIC-BASE-WORKING-DOC-FN-01` and
`HP-PUBLIC-BASE-WORKING-DOC-TEST-01` maintain the accepted documentation from
commit `5c4f5f363cb5c955a7046e9f287c0dbb67f19abb` for the existing exported
`cartesian_base_working_basis` function. The export remains
expert/unstable. It is a staged inspection and composition boundary, not a
second supported Hamiltonian facade; ordinary users should call
`cartesian_base_hamiltonian`.

The function accepts the same validated `system::NamedTuple` and
`basis::NamedTuple` descriptions used by the current Cartesian base producer,
plus its existing expert `supplemented` and `source_mode_overrides` controls.
It resolves those inputs and the Coulomb expansion, constructs the parent, and
realizes the terminal basis in memory. It does not itself assemble the complete
one-body operator, nuclear attraction, `Vee`, or
`CartesianIDAHamiltonian{Float64}`.

The current staged return carries these fields:

- `input`;
- `parent`;
- `terminal_basis`;
- `coulomb_expansion`;
- `source_mode_provenance`;
- `terminal_inventory`; and
- `terminal_due_diligence`.

This list describes present inspection/composition state, not a stable public
schema. Field names, nested private types, source-mode records, terminal
inventory rows, and due-diligence report layouts are expert/unstable and are
not compatibility promises. `terminal_inventory` and
`terminal_due_diligence` may be `nothing` when the selected expert route has
no terminal realization. The `source_mode_overrides` hook remains
expert/research-only; documenting the constructor does not promote its current
record layout into public policy.

The constructor writes no Hamiltonian artifact and provides no solver,
correction, PRF wrapper, automatic paper workflow, or complete operator
assembly. PRF-specific definitions remain unexported and receive no public
status through this staged constructor.

The accepted implementation added a `21`-line docstring in
`src/cartesian_base_hamiltonian.jl`, an `18`-line paragraph and curated entry
beside `cartesian_base_hamiltonian` in `docs/src/reference/export.md`, and
changed focused checks in `test/docs/runtests.jl` by `+13/-1`. All additions
remain within the authorized hard limits, with no new file.

The documentation check replaces the former count-only assertion with
exact equality between `Docs.undocumented_names(GaussletBases)`, excluding the
module self-binding, and this temporary five-name set:

- `OneCenterAtomicNestedStructureDiagnostics`;
- `one_center_atomic_nested_structure_diagnostics`;
- `one_center_atomic_nested_structure_report`;
- `diagnose_qwrg_residual_space`; and
- `ShellLocalAngularProfileKey`.

Exact set equality must reject substitution by any newly undocumented export.
This is a temporary regression boundary pending explicitly authorized
next-minor namespace transactions, not a general allowlist or de-export grant.
All five implementations remain available for qualified internal and research
use, and none may be de-exported by this packet. This namespace state does not
begin v0.3 development.

`cartesian_base_working_basis` remains exported, has a Julia docstring, and
appears in the curated reference with repeated, unambiguous expert/unstable and
non-schema wording. The undocumented-export backlog fell exactly `6 -> 5`.
Package load, docs `157/157` and `10/10`, public Cartesian `232/232` and
`80/80`, authority check/self-test, local Documenter, and diff checks passed.
Remote CI run `33579964076` passed Screening, PQS, and Supported floor; Docs run
`33579964078` passed and deployed the verified live reference. No definition,
signature, field, route, default, method, dispatch, staged result, PRF status,
numerical result, test owner, workflow, dependency, version, tag, or release
may change under maintenance. Any such work, any change to the exact five-name
set, or any de-export or v0.3 action requires separate docs-only authority.

## Supported Public Surface Documentation

`HP-PUBLIC-SUPPORTED-SURFACE-DOC-FN-01` and
`HP-PUBLIC-SUPPORTED-SURFACE-DOC-TEST-01` maintain documentation for the six
supported-public bindings separated from the remaining export backlog,
implemented by commit `8161f131aa962fef979f8ef09c14d23231eb14e4`:

| Binding | Required public meaning |
| --- | --- |
| `BondAlignedDiatomicQWBasis3D` | Describe the narrow bond-aligned mapped-product container and recommend its public homonuclear and heteronuclear builders without implying arbitrary molecular geometry. Move the currently misattached docstring from `AbstractBondAlignedOrdinaryQWBasis3D` to this concrete type rather than duplicating it. |
| `CoulombGaussianExpansion` | Represent the finite sum `sum_i coefficients[i] * exp(-exponents[i] * r^2)` with matching counts, positive exponents, copied `Float64` inputs, callable evaluation, and recorded `del`, `s`, `c`, and `maxu` generator parameters. It is neither an exact Coulomb operator nor a universal-interval accuracy claim. |
| `basis_metadata` | Return metadata associated with a supported basis, layer, representation, or supplement. Concrete type, fields, and ownership/copy behavior remain owner-specific; no universal schema is promised. |
| `cartesian_base_hamiltonian` | Return `CartesianIDAHamiltonian{Float64}` and optionally write the existing artifact for only the stated origin-centered H and z-aligned H2 base routes. General molecules, supplements, corrections, route controls, and solver behavior remain outside this facade. |
| `external_gto_ordering_fingerprint` | Provide a strict packet-integrity hash over ordered AO labels and stored angular powers, centers, primitive data, and normalization convention. It is not permutation-, tolerance-, or convention-invariant. |
| `external_gto_overlap_fingerprint` | Provide a strict packet-integrity hash over column-major `Matrix{Float64}` values. Dimensions and scientific consistency are validated elsewhere; the hash is not a numerical-equivalence, tolerance-invariant, permutation-invariant, or convention-invariant test. |

Maintenance is limited to the accepted docstrings in
`src/ordinary/ordinary_qw_types_and_bases.jl`, `src/ordinary/ordinary_coulomb.jl`,
`src/GaussletBases.jl`, `src/cartesian_base_hamiltonian.jl`, and
`src/cartesian_external_gto_import.jl`, and their curated entries or short
limiting prose in `docs/src/reference/bases_and_mappings.md`,
`docs/src/reference/operators_and_diagnostics.md`, and
`docs/src/reference/export.md`. The implementation added `52` and removed `10`
source-docstring lines, added `20` reader-documentation lines, and added no
file or executable behavior. The misplaced abstract-type docstring was removed
rather than duplicated.

The accepted `18`-line family-scoped check addition in
`test/docs/runtests.jl` requires all six names to remain exported, documented,
and present on their intended reference pages with the geometry, finite-
expansion, owner-specific metadata, narrow facade, and strict fingerprint
limits above. Documenter remains the executable `@docs` resolution gate. The
undocumented-export count, excluding the module self-binding, fell exactly
from `30` to `24`.

Docs tests passed `134/134` plus `10/10`; package load, authority check/self-
test, Documenter, and diff checks also passed. Remote CI run `33235720583`
passed all three numerical gates and Docs run `33235720615` passed. Both
records are in maintenance.

Preserve every implementation, constructor, field, method, signature,
dispatch, export, numerical behavior, artifact, hashing algorithm,
compatibility floor, workflow, tag, and release artifact. Add no global
`checkdocs` policy or allowlist, source behavior, test outside the docs family,
dependency, example, schema, API reduction, or release work. The remaining
classified bindings are outside this maintenance record; their refined
disposition follows and grants no automatic documentation or removal.

## Expert Geometry Documentation And Namespace Reservations

The read-only de-export audit at snapshot
`0832bdce0f3ff9bdd07d8580bedcda0c47dadbe8` refines the remaining public-
surface classification without granting any export change.

`bond_aligned_diatomic_geometry_payload` remains exported as an
expert/experimental interface. It is the central constructor in a coherent
bond-aligned geometry family with active high-order consumers, and owns
overloads for the supported bond-aligned QW basis, ordinary operators, nested
source and fixed-block combinations, and residual-augmented operators.
Removing its export would split that family rather than remove an orphan.

Two future next-minor namespace transactions are reserved but not authorized:

1. `OneCenterAtomicNestedStructureDiagnostics`,
   `one_center_atomic_nested_structure_diagnostics`, and
   `one_center_atomic_nested_structure_report` must move together; and
2. `diagnose_qwrg_residual_space` and its already-documented result type
   `QWRGResidualSpaceDiagnostics` must move together.

Every implementation remains available for qualified internal and research
use. These are namespace reductions, not source-code deletion, numerical,
performance, migration, or release work. They authorize no de-export, version
change, changelog entry, v0.2.x compatibility change, or v0.3 development.

Implementation commit `312f69a4b757e8d559387e0ef14450f031691039`
added the accepted concise docstring on the existing
`bond_aligned_diatomic_geometry_payload` generic declaration in
`src/GaussletBases.jl` and the compact section titled "Expert/experimental
bond-aligned geometry inspection" in
`docs/src/reference/bases_and_mappings.md`. The section curates the central
function with its already-documented exported geometry records and companion
operations without duplicating their docstrings:

- `BondAlignedDiatomicGeometryPoint3D`;
- `BondAlignedDiatomicGeometryNucleus3D`;
- `BondAlignedDiatomicGeometryBox3D`;
- `BondAlignedDiatomicGeometryPayload3D`;
- `BondAlignedDiatomicGeometryPlaneSlice3D`;
- `bond_aligned_diatomic_source_geometry_payload`; and
- `bond_aligned_diatomic_plane_slice`.

The documentation must state that the function constructs inspection and
visualization data without constructing or mutating a basis or Hamiltonian.
Nested payloads show compressed fixed centers, while
`bond_aligned_diatomic_source_geometry_payload` separately shows raw source-
region points. Multi-object overloads require the same basis geometry where
the current methods enforce it. The interface remains tied to current bond-
aligned QW/nested geometry and provenance conventions; it is not a general
molecular-geometry API, plotting framework, or stable serialization format.
Returned records and contained vectors are inspection data, not mutable
construction state.

The implementation added `18` source-docstring, `21` reader-documentation,
and `13` focused docs-test lines, all within preferred limits and with no new
file. The focused checks cover documentation presence, curated placement,
limiting language, and expert/experimental classification. Documenter remains
the executable `@docs` resolution gate. The undocumented-export count,
excluding the module self-binding, fell exactly from `24` to `23` and then
stopped. The accepted remaining state is nineteen undocumented expert/
experimental bindings plus the four original de-export candidates above,
with the already-documented `QWRGResidualSpaceDiagnostics` joined to the
future residual transaction.

Docs tests passed `138/138` plus `10/10`; radial passed `322/322`; package
load, authority check/self-test, Documenter, and diff checks passed. Remote CI
run `33326473189` passed all three gates, Docs run `33326473214` passed, and
the `/dev/` page was verified at the implementation commit. Both records are
now in maintenance.

Preserve every method, dispatch, field, payload value, export, numerical
behavior, plotting dependency, serialization behavior, ordinary test,
example, compatibility floor, workflow, tag, and release artifact. Add no
global documentation policy, implementation helper, metadata, file, export
change, migration note, version change, changelog entry, or release work. If
truthful documentation requires behavior or interface work, stop without an
implementation commit and report the exact gap.

## Paper-Aligned CI Boundary

`HP-PUBLIC-PAPER-CI-FN-01` and `HP-PUBLIC-PAPER-CI-TEST-01` own one
implemented, bounded test/workflow contract. PQS and reference-density Hartree
screening are separate publication surfaces and must appear as independently
named CI gates. This boundary adds no scientific or release authority.

Implementation commit `e65764377bd4640916e80342071da754d80aca32` landed the
three-row matrix, both release groups, and the mechanical screening-test
relocation. Commit `a4e85e820fd4056e985a18e20da87180f370ef66` changed only
the inherited workflow timeout from `15` to `30` minutes. Remote CI run
`32606409367` then passed all three unchanged gates: `Supported floor` on
Julia `1.10.12` passed `14,133/14,133` in `5m21s`; `Screening paper` on Julia
`1.12.7` passed `22/22` plus its example smoke in `1m03s`; and `PQS paper` on
Julia `1.12.7` passed `18/18` plus its example smoke in `18m53s` total. Docs
run `32606409363` also passed. The earlier PQS cancellation was therefore an
operational timeout, not a numerical failure. Both records are in maintenance.

The CI matrix has exactly three rows:

| Gate name | Julia | Test selection |
| --- | --- | --- |
| `Supported floor` | `1.10` | `core,ida,cartesian,examples,radial,misc,angular_public` |
| `PQS paper` | `1.12` | `pqs_release` group |
| `Screening paper` | `1.12` | `screening_release` group |

Julia `1.10` remains the declared compatibility-floor evidence. Julia `1.12`
is the canonical paper-gate environment; this does not add Julia `1.11`,
nightly, a new compatibility declaration, or a manifest policy.

Pass 542 extended only the Supported-floor selection shown above. Commit
`15676153aec1569f5224ffa6ff5ed67b054c837f` replaced the existing group-list
value and added one exact workflow-policy assertion, for a total `+2/-1` delta.
Clean Julia `1.10` evidence passed `radial` at `322/322`, `misc` at `59/59`,
and the combined selection at `381/381`, with about `70` seconds of testset
time beyond normal fresh-environment preparation. Remote CI run `33141930944`
passed the expanded Supported-floor gate and the unchanged PQS and Screening
gates; Docs run `33141930932` passed. All other workflow behavior, including
the three job names and rows, Julia versions, timeout, commands, classifier,
documentation-only markers, triggers, permissions, PQS and Screening groups,
and tag lane remains unchanged. Both records are back in maintenance.

The complete `angular` research group remains explicitly excluded. A Julia
`1.12` audit reached its first package-owned angular one-body fixture and then
spent `13m49s` without reaching an assertion, so it is neither a fast group nor
a per-push CI selection. Commit
`db156966a4a3a7bf2f685fa0f89312afca7b4280` implemented the separately
authorized `[10, 15, 32]` fixed-radial extraction: an `80`-line,
`83`-assertion `angular_public` owner replaced the `163`-line research-suite
block. The complete test delta was `+88/-165`, net `-77` lines. The owner is
explicitly selectable and included by `all`, while both `angular_public` and
bare `angular` remain outside `fast`; only `angular_public` was appended to the
Supported-floor row above. It passed `83/83` on Julia `1.10` in `11.17`
seconds. The single accepted Julia `1.12` combined `angular,angular_public` run
passed in `930.34` seconds and is not repeated during lifecycle closeout.
Remote CI run `33509253422` passed all three named gates and visibly ran the
focused owner `83/83` in Supported floor. The weekly Cartesian suite,
occupied-first private owner, represented-Hartree blocked owner, broader
HFDMRG portability, documentation triage, and export triage remain outside
this amendment. Add no CI row, job, dependency, source or numerical-policy
change, or release action.

The `pqs_release` group originally ran the complete comparison twice: once in
`test/pqs_h2plus_table1_release_runtests.jl` and once in an Example 41
subprocess. The two executions together measured `110.56 s / 19.342 GiB`.
Commit `b0dbd9ea37317590334a24883ef0667bdb0195a5` implements the bounded
replacement. Example 41 remains independently executable, writes the same
three-row/eight-column TSV and concise summary, and returns its already-computed
comparison as its final top-level value. The release owner includes that
example exactly once and applies all existing `18` assertions to the returned
comparison. The runner deletes only the duplicate Example 41 subprocess
assertion. Result types, tolerances, topology, accounting, public API checks,
nonfinite-reference rejection, required-check identity, and documentation links
remain unchanged. The accepted gate passed `18/18` at `54.91 s / 9.716 GiB`,
down from `110.56 s / 19.342 GiB`; the three-file delta is `+4/-5`.
The emitted TSV has SHA-256
`9a540cffff3b9ed87076063ed75626c55c2fda90457f4ac5a01e303ab3289d83`.
All three remote gates passed in CI run `33079914157`, and Docs run
`33079914164` passed.

The `screening_release` group must move the complete existing
`Public supplied-field screened Hartree` testset mechanically from
`test/ida/runtests.jl` to
`test/driver_public/screened_hartree_runtests.jl`. Every assertion, numerical
tolerance, oracle, malformed-input case, and physical-H1 check must remain
unchanged, and the old copy must be deleted in the same commit. Add exactly
one runner assertion that executes
`examples/40_screened_hartree_fixed_density.jl` through the existing helper.
Do not duplicate the public testset or substitute the `481`-line nested
screening owner.

Pass 525 reopened these records only for a path-aware routing repair. Commit
`3cce96d40f9e4f06f23a190f782834271fca884b` implemented that repair with
`74` workflow lines and `22` focused documentation-test lines. It changed no
source, numerical test owner, dependency, helper, or lifecycle vocabulary. The
implementation commit changed `.github/workflows/ci.yml`, correctly classified
itself as full, and passed the three unchanged jobs in CI run `32895536531`:
`Supported floor` in `5m58s`, `PQS paper` in `17m52s`, and `Screening paper`
in `45s`. Independent Docs run `32895536562` passed.

The three matrix rows, job names, Julia versions, test groups, `30`-minute
timeout, job-level `contents: read`, package instantiate/load boundary,
disabled slow tests, and pull-request behavior remain exact. Pull requests
always run all three rows. The accepted single-execution change affects only
the bounded test/example composition above; no workflow or numerical contract
changes with it.

For a push to `main`, the classifier must begin with `scope=full`. It may emit
`scope=docs_only` only when all of these conditions hold:

1. checkout used `fetch-depth: 0`, or an equally narrow explicit fetch made
   the exact `github.event.before` object available;
2. `github.event.before` is nonzero and both it and `github.sha` resolve as
   commit objects;
3. `git diff --name-only --no-renames <before> <sha> --` succeeds and returns
   at least one path; and
4. every changed path is in this exact allowlist:

```text
AGENTS.md
docs/**
test/docs/**
.github/workflows/docs.yml
```

Missing, zero, unfetchable, malformed, empty, renamed-to-hidden, or unknown
state remains full-matrix. In particular, `Project.toml`, root `README.md`,
root `CHANGELOG.md`, `.github/workflows/ci.yml`, source, ordinary tests,
examples, and every unlisted path are full-matrix. The implementation must not
use workflow-level path filters or add a second workflow.

Classifier errors must not suppress numerical validation through ordinary
dependency skipping. The matrix job must retain its three existing names and
use an `always()`-aware job condition that expands all three rows for every
non-tag event. A shared fail-closed step predicate skips numerical work only
for the single explicit successful `docs_only` value; every other result runs
the unchanged numerical steps. The lightweight job runs only after an explicit
successful `docs_only` result and performs package load plus the existing docs
test group. An absent or invalid classifier output is not a third lane.

`v*` pushes use a separate tag-CI lane and do not rerun the numerical matrix.
That lane must explicitly fetch the remote tag ref and prove that
`refs/tags/<tag>^{tag}` exists as an annotated tag object. It may only then
check the `^{commit}` peel, tree, canonical semantic-version and
`Project.toml` agreement, fresh installation from the remote tag, and package
load. Checkout `HEAD` and the event SHA are not substitutes for preserving or
verifying the tag object. Malformed or lightweight tags fail.

Tag CI does not own deployed-documentation acceptance. After tag CI succeeds,
repo-manager must independently wait for the already-triggered Docs/Pages
workflow and inspect the versioned canonical URL, `versions.js`, `/dev/`, and
the exact prerelease/final `stable` policy before any separately authorized
publication. Add no workflow chaining, `workflow_run` trigger, polling action,
or checked release helper. This CI amendment authorizes no tag, release,
registration, citation, stable alias, or final-v0.2 operation.

The identified `481`-line screening owner and five related research/internal
suites under `test/nested` remain unwired. Their classification as maintained
scheduled coverage, quarantine, or deletion is a separate bloat-review
decision. No claim about their value or lifecycle is made here.

Before the bounded G4 extension below, maintenance ownership for the
path-aware routing was limited to:

```text
.github/workflows/ci.yml
test/docs/runtests.jl
```

Production source and numerical-test changes are zero. The G4 extension below
supersedes only the former prohibition on one mechanical test group and one
shared test owner. Any other scope expansion or change to the three scientific
gates requires new design authority.

Transition acceptance was intentionally staged. The first closeout probe,
commit `72446880603c7e554f6ae71b2de2dc6edf28b31b`, correctly reached the
documentation-only classification, and independent Docs run `32899962286`
passed. CI run `32899962285` nevertheless rejected closure for two mechanical
reasons:

1. GitHub evaluated the matrix job's false job-level condition before matrix
   expansion, yielding one skipped job named `matrix.gate` rather than the
   three accepted gate names; and
2. the lightweight job instantiated only the root environment, so the existing
   isolated Documenter version-expansion fixture failed at `105/106` because
   the docs environment had not been instantiated.

Pass 528 reopened only the same two records for this exact repair. Commit
`9ddc689c1bc806c7ec899cac7a39d77cb7fad3bf` implemented it with `13` added
and one deleted workflow line plus `9` focused policy-test lines. The matrix
now expands on every non-tag event. On an explicit successful `docs_only`
classification, each of its three named rows runs one visible no-op marker and
skips checkout, Julia setup/cache, root instantiation/load, and numerical test
execution; every other non-tag state runs the unchanged full steps. The
separate lightweight job instantiates the existing docs environment with the
same local-package develop/instantiate operation used by
`.github/workflows/docs.yml` before invoking the existing docs test group. No
fourth gate, workflow, group, or numerical command was added.

Because the repair commit changed `.github/workflows/ci.yml`, it correctly took
the full path. CI run `32901325992` passed `Screening paper` in `50s`,
`Supported floor` in `6m`, and `PQS paper` in `17m01s`; Docs run `32901326008`
passed. Local policy checks passed `110/110` and `10/10` together with YAML,
package load, authority check/self-test, Documenter, manager-log, and diff
checks. Both records return to maintenance after the corrected docs-only
transition below succeeds.

The corrected transition remains staged:

1. the docs-only authority commit ran the legacy full matrix;
2. workflow implementation commit `3cce96d40` changed
   `.github/workflows/ci.yml`, classified itself as full, and passed all three
   unchanged named jobs;
3. repair implementation commit `9ddc689c1` changed
   `.github/workflows/ci.yml`, classified itself full, and passed all three
   unchanged numerical jobs plus Docs;
4. this docs-only lifecycle closeout must pass its lightweight
   package/docs checks, show `Supported floor`, `PQS paper`, and `Screening
   paper` as successful jobs whose numerical steps are visibly skipped, and
   pass the independent Docs workflow; and
5. the tag-lane commands were rehearsed read-only against immutable
   `v0.2.0-rc2`: annotated object `7c8a21b998a838d245e0b5a7f4915910e2a091bc`
   peeled to commit `2b3c23970144aa030ae52b875a5cf01b32886b6e` and tree
   `7a4b51aec25f62436620f4ff938262d0f6b2fd62`, and a fresh remote-tag
   installation loaded the package without moving, deleting, recreating, or
   pushing a tag.

Acceptance also requires YAML parsing, exact allowlist and synthetic
classification checks, authority check/self-test, generated-view parity, docs
tests, package load, Documenter, manager-log bound, remote CI and Docs at the
stages above, and `git diff --check`. No new tag is authorized. The accepted
release-efficiency process reduces future candidate publication to three
design-manager lifecycle passes, but each exact tag/publication transaction
still requires separate version-specific authority with frozen candidate and
release identities.

### Lightweight Source-Push Public-Surface Integrity

Pass 585 reopens `HP-PUBLIC-PAPER-CI-FN-01` and
`HP-PUBLIC-PAPER-CI-TEST-01` only for the independently reviewed G4 boundary.
The implementation adds one `docs_fast` group to the existing Julia `1.10`
`Supported floor` selection. It changes no matrix row, job or required-check
name, Julia version, timeout, classifier, marker, command, trigger, permission,
PQS/Screening gate, tag rule, or `fast` alias.

One new shared owner, `test/docs/public_surface_runtests.jl`, must contain only
the existing mechanical public-surface contract: the `95` named documented
exports remain exported, each has documentation, and the complete
undocumented-export set excluding the module self-binding is exactly:

```text
OneCenterAtomicNestedStructureDiagnostics
one_center_atomic_nested_structure_diagnostics
one_center_atomic_nested_structure_report
diagnose_qwrg_residual_space
ShellLocalAngularProfileKey
```

Move the corresponding inventories and mechanical assertions out of
`test/docs/runtests.jl`; do not duplicate them. Reference-page placement,
reader wording, version-menu behavior, Documenter fixtures, authority checks,
deployment policy, and every other prose assertion remain in the full `docs`
owner. The shared owner must require only the root test environment and package
load, never the docs environment or Documenter.

Julia `1.10.12` exposes neither `Base.Docs.hasdoc` nor
`Base.Docs.undocumented_names`. A small test-local semantic fallback may inspect
Base Docs binding metadata and resolve aliases; no production compatibility
shim or dependency is permitted. On Julia `1.12.6`, it must agree with
`Base.Docs.hasdoc` for every current public name, reproduce
`Base.Docs.undocumented_names` after excluding the module self-binding, and
pass an explicit documented-alias check. On Julia `1.10.12`, it must yield the
exact five-name set above.

`test/runtests.jl` must expose `docs_fast` for explicit selection and `all`,
and include the shared owner exactly once whenever either `docs_fast` or
`docs` is selected. Thus `docs_fast`, `docs`, `docs_fast,docs`, and `all` each
execute the mechanical owner once. The existing `fast` alias remains byte-for-
byte unchanged. In `.github/workflows/ci.yml`, append only `docs_fast` to the
existing Supported-floor group string; bare `docs` and every numerical group
remain unchanged.

Budget: one new test file at preferred `92-98`, hard `110` lines; move about
`58` inventory lines and remove or simplify about `24` mechanical assertion
lines from the full docs owner; runner additions preferred `4`, hard `8` lines;
workflow delta exactly `+1/-1`; total tracked-test growth preferred `15-25`,
hard `35` lines. Broad phrase-test pruning is forbidden. If cross-version alias
or exact-set equivalence fails, execution duplicates, the Julia `1.10` package
test requires a docs environment, or any existing numerical selection changes,
make no implementation commit and report the mismatch.

Acceptance requires focused `docs_fast`, full `docs`, combined
`docs_fast,docs`, and `all` inclusion/control-flow evidence without running all
numerical groups merely to count the include; direct Julia `1.10.12` and
`1.12.6` semantic checks; unchanged normal docs/authority/Documenter checks;
and a full implementation-push CI because `.github/workflows/ci.yml` changes.
The Supported-floor log must visibly execute `docs_fast`; PQS and Screening
remain unchanged. No G3 whitelist relocation, prose cleanup, Step 4, source,
API, export, dependency, numerical-policy, release, or new workflow authority
is granted.

## Compatibility

Historical Table I records remain immutable provenance. The public replay must
agree semantically and within the declared tolerances, not reproduce the
composite TSV byte-for-byte. No old staged object, paper-driver report, or
internal packet object becomes a public serialized format.

The accepted historical He screening result remains data-only. The public
assembly example validates the stable algebra but does not claim its basis,
field, or energy reproduces that calculation. Retired moment-polish packet
provenance remains rejected. Existing internal packet readers retain only
their current internal compatibility contract.

HFDMRG and its solver-dependent paper rows remain separately versioned and do
not block this GaussletBases surface.

## Implementation Surface And Budgets

Approved source ownership:

```text
src/GaussletBases.jl
src/pqs_matched_h2plus.jl
src/cartesian_reference_density/CartesianReferenceDensity.jl
src/cartesian_reference_density/screened_hartree_correction.jl
```

Approved test ownership:

```text
test/ida/runtests.jl
test/runtests.jl
test/pqs_h2plus_table1_release_runtests.jl
```

The focused release/nightly workflow and two example files are allowed only
for the contracts above. No general example runner or release framework is
approved.

Stop-and-report source budgets, including root wiring, are:

- matched H2+ wrapper: preferred at most `300`, hard at most `360` added
  source lines;
- public screened-Hartree layer: preferred at most `90`, hard at most `130`
  added source lines;
- examples together: hard at most `180` lines;
- tests together: hard at most `220` lines.

These are review and stop conditions, not instructions to compress equations,
validation, or readable implementation. If a coherent implementation exceeds
a hard limit, make no source commit and return the smallest justified revision
to this design manager.

The private paper driver may remain because it owns broader historical report,
supplemented-H2, and fixed-state operations. Do not copy its input parser,
report schema, timing fields, or unrelated interaction path into the public
wrapper. Reuse accepted numerical owners directly and account for any
remaining duplicate H2+ mechanics in the handoff.

## Explicit Non-Goals And Failure Rule

This authority adds no public staged parent/shell/transform object, packet
builder or reader, mixed-basis Hartree producer, density fitter, potential
fitter/evaluator, solver, artifact, checkpoint, Ximg/XHF surface, HFDMRG
surface, interaction/exchange correction, general four-index engine, or
release machinery. The single annotated tag operation above is separately
owned and adds no public interface or tracked implementation.

If either example requires a new scientific approximation, public stage
schema, general field-construction engine, or paper-local data dependency,
stop without a source commit and report the exact missing operation to this
same design manager.

## v0.2.0-rc2 Reader Front Door And Candidate Preparation

`HP-PQS-PUBLIC-RC2-FN-01` owns maintenance of the exact candidate prepared by
commit `2b3c23970144aa030ae52b875a5cf01b32886b6e`.
`HP-PQS-PUBLIC-RC2-TEST-01` owns completed acceptance evidence and grants no
further execution. RC2 is required because the public residual-GTO/MWG
working-system facade and external Cartesian-GTO/PySCF transfer were added
after immutable RC1. This candidate record does not authorize a tag; the exact
annotated operation is isolated below. GitHub release, registration, citation
record, stable documentation, and final `v0.2.0` remain unauthorized.

### Reader Front Door

The root `README.md` names three distinct capabilities:

- Projected q-shells;
- reference-density Hartree screening; and
- external Cartesian GTO transfer.

The radial workflow remains the recommended beginner route. PQS and screening
must remain separate method surfaces. The external exporter must be described
only as a checkpoint reader; PySCF and NumPy remain optional external-command
dependencies. The README may add no tutorial duplication, benchmark or
manuscript claim, paper status, DOI, or publication promise. Its stop-and-
report budget is `14` preferred and `20` hard added lines.

### Candidate Files

The accepted candidate contains only these changes:

1. change the root `Project.toml` version from `0.2.0-rc1` to
   `0.2.0-rc2`;
2. add one concise `v0.2.0-rc2` section above the byte-unchanged RC1 section
   in `CHANGELOG.md`;
3. add the exact `v0.2.0-rc2 => v0.2.0-rc2` selector beside the retained RC1
   selector and `dev` in `docs/make.jl`;
4. make the existing deployment contract describe that exact selector; and
5. focused existing docs tests for the RC2/version-selector contract.

The implementation changed five existing files by `+48/-2`: `8` README
lines, an `18`-line RC2 changelog section, one version replacement, one
selector line, and `20/1` focused docs-test lines. The RC1 section is
byte-identical to its accepted predecessor; its section SHA-256 remains
`0c92db23146e520ff0365e7724eb9b9fa1fb4e8d453c1544be4f7cdaa6a3a42d`.
No source, API, numerical, dependency, example, workflow, fixture-format, or
manifest-policy behavior changed.

The RC2 changelog section may summarize only the public residual-GTO/MWG
working-system constructor, strict version-1 external Cartesian-GTO reader,
checkpoint-only PySCF exporter, caller-thresholded closest-determinant
operation, separate Supported-floor/PQS/Screening CI gates, and removal of
invalid package exports. Keep it grouped and reader-facing, within `28`
preferred and `40` hard added lines. Do not alter the RC1 section.

The exact RC2 selector is one added `docs/make.jl` entry; the hard bound is
`3` added lines. Deployment-contract additions are bounded at `18` preferred
and `28` hard lines. Existing docs-test additions are bounded at `20`
preferred and `35` hard lines. No new file is allowed. Readability and exact
policy checks take precedence over compressing to a limit.

No production source, API, export, numerical rule, dependency, example,
workflow, fixture format, manifest policy, or test owner may change. The two
established untracked handoffs remain outside the candidate and its archive.
The accepted C2/aug-cc-pV6Z replay remains canonical evidence and need not be
repeated because this pass changes no implementation.

### Candidate Acceptance

Acceptance preserved the exact candidate commit and established:

- clean installation and package load on Julia `1.10` and the verified Julia
  `1.12` line;
- all three public CI gates and Docs passing, including the dynamic export-
  integrity regression;
- the public residual-GTO/external-transfer owner and frozen H2/cc-pVTZ
  fixture checks;
- public examples 01, 39, 40, and 41 plus the focused H2+ `18/18` release
  gate;
- an archive/install replay containing no root `Manifest.toml` or untracked
  handoff; and
- RC2, RC1, and `dev` selector simulation with no `stable` entry, alias, or
  path.

Local docs passed `99/99` and `10/10`; the residual-GTO/external-transfer owner
passed `80/80`, H2+ passed `18/18`, and screening passed `22/22`. Remote Docs
run `32780718822` passed. CI run `32780718923` passed Supported floor in
`5m57s`, PQS paper in `17m09s`, and Screening paper in `1m43s`. The clean
archive contains `676` entries in `9,994,240` bytes, has SHA-256
`6728b80c1397f13b367c2d898fbdda3176c6cb87c39597817c590a0f41c1e2ac`,
and contains neither a root Manifest nor either handoff.

The RC2 selector coexists with RC1 and `dev`. The completed tag lifecycle below
publishes the canonical RC2 folder while `/dev/` remains live and `/stable/`
remains absent. Candidate acceptance grants no further publication authority.

## v0.2.0-rc2 Annotated Tag Lifecycle

`HP-PQS-PUBLIC-RC2-TAG-FN-01` records the completed tag-only publication. The
immutable target is candidate commit
`2b3c23970144aa030ae52b875a5cf01b32886b6e`, with tree
`7a4b51aec25f62436620f4ff938262d0f6b2fd62`. Its accepted tracked archive
contains `676` entries in `9,994,240` bytes and has SHA-256
`6728b80c1397f13b367c2d898fbdda3176c6cb87c39597817c590a0f41c1e2ac`.
Annotated tag object `7c8a21b998a838d245e0b5a7f4915910e2a091bc` has
message `GaussletBases v0.2.0-rc2` and peels exactly to that target and tree.
Only `refs/tags/v0.2.0-rc2` was pushed; no branch commit, file, or asset changed.

`HP-PQS-PUBLIC-RC2-TAG-TEST-01` is completed. Tag-triggered Docs run
`32798625043` passed, and CI run `32798625045` passed the PQS, Supported-floor,
and Screening gates. Pages deployment `32798719038` passed.
`/v0.2.0-rc2/` is live with its exact canonical URL, `versions.js` lists RC2,
RC1, and `dev`, `/stable/` remains absent, and `/dev/` remains intact. `main`
and `origin/main` remained at `707e5bcbb35f375a09c663cfd949f40377d37e19`,
with both handoffs untouched.

The tag lifecycle is closed. Never move, replace, delete, recreate, force, or
silently retry the immutable tag. The corresponding GitHub prerelease
lifecycle is also closed below.

## v0.2.0-rc2 GitHub Prerelease Lifecycle

`HP-PQS-PUBLIC-RC2-RELEASE-FN-01` records the completed GitHub prerelease for
the existing immutable `v0.2.0-rc2` tag. Published release
[`376503169`](https://github.com/srwhite59/GaussletBases.jl/releases/tag/v0.2.0-rc2)
has:

- tag `v0.2.0-rc2`, verified rather than created by the release operation;
- title `GaussletBases v0.2.0-rc2`;
- prerelease status `true`, draft status `false`, and latest status `false`;
- no uploaded assets, so only GitHub's automatic source archives appear; and
- the exact `2,200`-byte ASCII release body below, including its final newline.

~~~~markdown
GaussletBases v0.2.0-rc2 is a public-testing prerelease of the GaussletBases
v0.2 package surface. It retains the RC1 PQS and screening surfaces and adds a
supported external Cartesian-GTO orbital-transfer path.

## New in RC2

- A public `cartesian_residual_gto_mwg_system` constructor for the supported
  residual-GTO/MWG working object.
- Strict version-1 external Cartesian-GTO packet reading and validation.
- A checkpoint-only PySCF exporter using explicit Cartesian AOs.
- Optional caller-thresholded closest-determinant preparation, while preserving
  the unchanged raw projection and its capture diagnostics.
- Separate Supported-floor, PQS, and Screening CI gates.
- Removal of invalid package exports that did not name usable bindings.

The exporter runs no SCF calculation. PySCF and NumPy are optional dependencies
of that external command, not Julia package dependencies. External transfer
produces orbitals and diagnostics, not another solver's restart or a transformed
Hamiltonian.

## Projected q-shells (PQS)

- Bounded public PQS/White-Lindsey H2+ comparison types and constructor.
- Public examples 39 and 41 cover the supported PQS construction and comparison.

## Reference-density Hartree screening

- Supplied-field correction assembly with stable public result accessors.
- Public example 40 demonstrates the bounded screening algebra.

PQS and reference-density Hartree screening remain distinct method surfaces.
This prerelease does not present them as one combined method or one
combined-paper validation.

## Compatibility and scope

- Julia 1.10 is the supported minimum; RC2 was replayed on Julia 1.10.12 and
  1.12.6.
- GaussletBases supplies basis/operator construction, orbital transfer, and
  bounded screening algebra; it is not a complete SCF or correlated molecular
  solver.

## Install

Install this prerelease directly from its immutable tag:

```julia
using Pkg
Pkg.add(url = "https://github.com/srwhite59/GaussletBases.jl", rev = "v0.2.0-rc2")
```

- [Versioned v0.2.0-rc2 documentation](https://srwhite59.github.io/GaussletBases.jl/v0.2.0-rc2/)
- [Changelog at the immutable tag](https://github.com/srwhite59/GaussletBases.jl/blob/v0.2.0-rc2/CHANGELOG.md)
~~~~

Because the body contains an inner Julia fence, the text between the outer
fence markers was copied into a temporary notes file outside the repository.
The completed noninteractive operation was:

```bash
gh release create v0.2.0-rc2 \
  --repo srwhite59/GaussletBases.jl \
  --verify-tag \
  --title "GaussletBases v0.2.0-rc2" \
  --prerelease \
  --latest=false \
  --notes-file /private/tmp/gaussletbases_v0.2.0-rc2_release_notes.md
```

No filename or pattern followed the tag; no generated notes, discussion,
draft, `--target`, or asset was supplied. Preflight verified synchronized
`main`, only the two established untracked handoffs, release absence, and the
exact local and remote tag object, peeled commit, and tree. The published body
is byte-identical to the canonical `2,200` bytes and has SHA-256
`a2cbaaa2a349857e897d6d58fb728c6ccd9d731c7371bb39997bfc0360f3653a`.

`HP-PQS-PUBLIC-RC2-RELEASE-TEST-01` is completed and grants no edit or
publication authority. GitHub API and rendered-page checks confirmed the
exact tag and title, `draft=false`, `prerelease=true`, no latest final release,
and zero uploaded assets. The automatic tarball is `2,564,018` bytes with
SHA-256 `5e92245b92350865415facf1350ba186b58052f0befb6380a5397d1eadaef445`;
the zipball is `2,957,570` bytes with SHA-256
`84a5619ab679650ceb6e7c0ab0e728d078cb0785d1b26d4f83f6e81471b24d30`.
Both reconstruct tree `7a4b51aec25f62436620f4ff938262d0f6b2fd62`. A fresh
isolated Julia `1.12.6` installation from the immutable tag loaded
GaussletBases `0.2.0-rc2` and recorded that same tree in its manifest.

RC2, RC1, and `/dev/` remain live, `versions.js` lists RC2, RC1, and `dev`,
and `/stable/` remains absent. `main` and `origin/main` remained at
`b341be8814e654eb6039e18d1033c3a71936019b`, with both handoffs untouched.

The release lifecycle is closed. Preserve the accepted release and immutable
tag without editing, deleting, recreating, retargeting, adding assets, changing
the narrative, or silently retrying either operation. Publication changed no
tracked file.

No RC tag mutation, RC1 release edit, custom asset, registration, citation,
stable alias, final `v0.2.0`, repository-metadata change, or tracked source,
API, test, example, dependency, workflow, numerical, or manuscript mutation
is authorized.

## v0.2.0 Final Candidate And Conditional Publication Process

`HP-PQS-PUBLIC-V020-FN-01` and `HP-PQS-PUBLIC-V020-TEST-01` authorize one
direct final-candidate preparation from the accepted post-RC2 `main`. No RC3
is required. This is candidate authority only: it does not create, push, or
publish `v0.2.0`.

The final package claim is deliberately modest:

> v0.2.0 is the supported public package version closest to the software used
> during the PQS and reference-density Hartree-screening work. It is not an
> exact archive of either paper's complete computational history.

PQS and reference-density Hartree screening remain distinct method surfaces.
The final release must not present them as one combined method, one combined
paper, or an exact reproduction archive. External Cartesian-GTO transfer is a
separate public interoperability surface. GaussletBases provides basis and
operator construction, orbital transfer, and bounded screening algebra; it is
not a complete SCF, correlated-solver, or paper-campaign package.

### Candidate Edits

The candidate implementation may make exactly these tracked changes:

1. change the root `Project.toml` version from `0.2.0-rc2` to `0.2.0`;
2. prepend one concise `v0.2.0` section to `CHANGELOG.md`, leaving the accepted
   RC2 and RC1 sections byte-identical;
3. change the README installation command to the immutable future tag
   `rev = "v0.2.0"` and point its rendered-documentation links at `/stable/`
   while retaining the radial workflow as the recommended beginner route; and
4. update only the existing focused documentation tests for the final version,
   changelog heading, README links, exact final-tag canonical path, and standard
   Documenter final-release selector behavior.

The changelog section must group, rather than narrate commit by commit, the
post-RC2 changes:

- support-local PQS shell-seed construction and safe call-local workspace and
  terminal-buffer reuse;
- the exact-order four-element terminal Gaussian-sum acceleration;
- fail-closed path-aware CI that preserves full candidate/code gates while
  using lightweight package/docs checks for proven docs-only main pushes; and
- one matched-H2+ comparison shared by Example 41 output and all unchanged
  release assertions instead of duplicate execution.

The changelog must state that these changes preserve the accepted public
numerics. It must not contain manager-pass history, detailed allocation or
timing tables, manuscript headline numbers, private paper-validation history,
paper titles or status, DOI invention, registration claims, or promises about
paper publication.

`docs/make.jl` already contains the standard `"stable" => "v^"` and `"v#.#"`
policy plus exact RC2 and RC1 prerelease retention. The candidate must not edit
that policy or the Docs workflow. Focused tests must retain the prerelease-only
fixture proving that RC tags do not create `stable`, then add or update a
bounded final-folder fixture proving that `v0.2.0` selects the real `stable`
alias while RC2, RC1, and `dev` remain available.

Stop-and-report budgets are:

- `Project.toml`: exactly one version-line replacement;
- `CHANGELOG.md`: preferred at most `18`, hard at most `24` added lines;
- `README.md`: preferred at most `4`, hard at most `8` added lines, with URL
  changes otherwise performed as replacements;
- `test/docs/runtests.jl`: preferred at most `20`, hard at most `30` added
  lines; and
- no new file.

These are review bounds, not reasons to obscure the release statement or test
logic. If a truthful candidate requires source, API, dependency, numerical,
example, workflow, fixture-format, manifest-policy, or other test-owner work,
make no candidate commit and return the exact blocker.

### Candidate Acceptance And Freeze

The candidate commit changes non-documentation-only paths and therefore must
take the full CI route exactly once. Acceptance requires:

- synchronized `main` and `origin/main`, with only the two established
  untracked handoffs and no tracked root `Manifest.toml`;
- exact version `0.2.0`, concise changelog traceability, and unchanged public
  exports, dependencies, examples, implementation, and numerical policy;
- package installation and load on Julia `1.10` and Julia `1.12.6`;
- all three named public CI gates, including export integrity, the public
  residual-GTO/external-transfer owner and frozen H2/cc-pVTZ fixture, screening
  coverage, and the single-execution H2+ `18/18` release owner;
- public examples 01, 39, 40, and 41 without a second Example 41 comparison;
- docs tests, local Documenter, authority checks, manager-log bound, and
  `git diff --check`;
- one clean candidate archive and fresh archive installation, excluding a root
  Manifest, ignored state, and both handoffs; and
- final/prerelease selector simulation proving `v0.2.0` may become `stable`,
  RC tags alone cannot, and RC2, RC1, and `dev` remain represented.

Candidate acceptance freezes the exact commit SHA, tree, archive entry count,
archive byte count, and archive SHA-256. Repo-manager must also return a
proposed exact ASCII GitHub release body, including its final newline, byte
count, and SHA-256. That body must be package-centered, use the modest claim
above, keep PQS and screening separate, summarize external transfer and the
post-RC2 maintenance improvements, provide the exact-tag installation command,
and link to immutable `v0.2.0` documentation and changelog locations. It must
not be published during candidate preparation.

### Accepted Final Candidate

The final candidate is accepted at commit
`adfcaba32d4db06d9d796d947276433717bd2d89` with tree
`f64ba21e06ff57e2b5e78d91214398115afbe8de`. It changes only
`Project.toml`, `CHANGELOG.md`, `README.md`, and `test/docs/runtests.jl` from
the Pass 538 authority baseline. The RC2 and RC1 changelog content remains
byte-identical with SHA-256
`fa6b3c4a64921a78fe52469925a48cb2a8fecc3890bcfddf2e4588d5ebfa2380`.
The documentation-test diff is `+38/-9`: nine additions are one-for-one
version/link replacements and the remaining `29` are new validation lines,
within the hard `30`-line substantive-addition bound.

The clean candidate archive contains `677` entries, has `10,137,600` bytes,
and has SHA-256
`df09cc6fd7dc144daa168c9feb4a41be9b974ef450e1e81bf586787318ad1566`.
It contains no root `Manifest.toml` and neither established handoff. Direct
GitHub installation and package load passed on Julia `1.10.12` and `1.12.6`;
fresh archive installation passed on Julia `1.12.6`; examples 01, 39, 40,
and 41, the H2+ release owner `18/18`, screening, export integrity,
residual-GTO/external-transfer fixture, docs `114/114 + 10/10`, Documenter,
authority, and diff checks passed. CI run `33126022579` passed all three
numerical gates and Docs run `33126022531` passed.

The exact release body is ASCII, includes its final newline, contains `2,278`
bytes, and has SHA-256
`e9ae9bcdad74b33bb66fb3e7e6a149d26285cb9bcc2f4c9555ac713be8bc90d2`:

````text
# GaussletBases v0.2.0

GaussletBases v0.2.0 is the supported public package version closest to the
software used during the separate PQS and reference-density Hartree-screening
work. It is not an exact archive of either paper's complete computational
history.

## Projected q-shells (PQS)

- Provides the bounded public PQS and matched White-Lindsey H2+ construction.
- Example 39 demonstrates a small public basis and one-body calculation.
- Example 41 produces the matched three-row H2+ comparison and TSV used by the
  package release gate.

## Reference-density Hartree screening

- Provides typed supplied-field screening assembly for exact represented and
  explicitly fitted reference fields.
- Example 40 demonstrates fixed-density energy and occupied-action closure.
- Atomic fitting, packet placement, SCF, and correlated solvers remain outside
  this public screening surface.

## External Cartesian GTO transfer

- Provides a strict versioned Cartesian-GTO packet reader, raw orbital import,
  overlap diagnostics, and an explicit caller-thresholded closest-determinant
  operation.
- The optional PySCF checkpoint exporter records the Cartesian AO convention,
  source overlap, occupations, and producer attestation needed by the reader.

## Post-RC2 maintenance

- Reduced PQS construction cost with support-local shell seeds and safe
  call-local workspace and terminal-buffer reuse.
- Accelerated terminal Gaussian sums with an exact-order four-element path.
- Made path-aware CI fail closed while retaining lighter checks for proven
  documentation-only main pushes.
- Removed the duplicate matched-H2+ calculation from the Example 41 release
  gate without removing its output or any of the 18 assertions.

## Install

```julia
using Pkg
Pkg.add(url = "https://github.com/srwhite59/GaussletBases.jl", rev = "v0.2.0")
```

Documentation: https://srwhite59.github.io/GaussletBases.jl/v0.2.0/

Changelog: https://github.com/srwhite59/GaussletBases.jl/blob/v0.2.0/CHANGELOG.md

GaussletBases provides basis and operator construction, orbital transfer, and
bounded screening algebra. It is not a complete molecular SCF, correlated-
solver, or paper-campaign package. The PQS and reference-density Hartree-
screening work remain distinct method and citation lines.
````

### Conditional Final Transaction

`HP-PQS-PUBLIC-V020-RELEASE-FN-01` authorizes one version-specific,
exact-hash conditional tag-plus-publication transaction. It freezes:

- candidate commit, tree, and archive identity;
- annotated tag `v0.2.0` with message `GaussletBases v0.2.0`;
- exact release title `GaussletBases v0.2.0` and the `2,278`-byte body above;
- `draft=false`, `prerelease=false`, latest-final status `true`, and zero
  custom uploaded assets; and
- the stop order between tag verification, documentation deployment, and
  GitHub publication.

The transaction must first prove local and remote tag absence, GitHub-release
absence, synchronized `main` and `origin/main`, and the exact candidate,
tree, archive, version, and release-body identities above. It then creates one
annotated tag explicitly at the frozen candidate and pushes only that tag. It
must stop if the remote tag object is not annotated, cannot be dereferenced as
`refs/tags/v0.2.0^{tag}`, or has a different peel, tree, version, or clean
remote installation. Repo-manager then waits for the independently triggered
Docs/Pages workflow and verifies `/v0.2.0/`, the real `/stable/` alias,
`versions.js`, `/dev/`, RC2, and RC1 before publishing the exact GitHub release.
No workflow chaining or polling framework is added.

The exact final tag does not rerun the three numerical gates. Candidate SHA and
tree identity make that redundant. The tag lane owns annotated-object
identity, peel/tree/version agreement, remote-tag installation, and package
load; independent Docs/Pages owns deployed documentation. Post-publication
checks own exact title/tag/body bytes and hash, `draft=false`,
`prerelease=false`, latest-final status, zero assets, automatic archive
reconstruction to the frozen tree, and a fresh Julia `1.12.6` remote-tag
installation. If tag verification or documentation deployment fails after the
push, preserve the immutable tag, publish no release, and report. If release
publication partially succeeds or any post-publication check differs, preserve
the release and tag and report without moving, deleting, recreating, silently
editing, adding assets, or retrying.

#### Immutable-Tag Verification Recovery

The final annotated tag was created correctly but the first tag CI run
`33130411193` failed before verification. `actions/checkout` left a local
lightweight `v0.2.0` ref at the peeled commit; the lane then attempted to fetch
the remote annotated object into that same local tag name without force, and
Git rejected the replacement as `would clobber existing tag`. This is an
operational defect in the workflow frozen inside the immutable candidate, not
a tag-identity mismatch.

The immutable remote tag is preserved as object
`722e8e8752a9d23f45e95d2f88e1749f9f3002e4`, peeling to candidate
`adfcaba32d4db06d9d796d947276433717bd2d89` and tree
`f64ba21e06ff57e2b5e78d91214398115afbe8de`, with message
`GaussletBases v0.2.0`. A fresh Julia `1.12.6` installation from
`rev = "v0.2.0"` loaded package version `0.2.0` with the frozen tree. Docs run
`33130411176` and Pages run `33130489319` passed; `/v0.2.0/` and `/stable/`
are identical with the exact-version canonical URL, `versions.js` retains
stable, `v0.2`, RC2, RC1, and `dev`, and all prior folders remain intact. At
recovery authorization time, no GitHub final release existed.

Because the tag cannot be changed or made to rerun corrected frozen workflow
code, the transaction may recover through one bounded manual verification
lane. It must use a fresh machine-local scratch repository, fetch the remote
`refs/tags/v0.2.0` into a namespaced temporary ref such as
`refs/verification/v0.2.0`, and prove all of:

- the fetched ref resolves through `^{tag}` and has Git object type `tag`;
- its object ID, exact message, `^{commit}` peel, `^{tree}`, and
  `Project.toml` version match the frozen identities above;
- `git ls-remote` and the GitHub API report the same remote object and peel;
- a fresh isolated Julia `1.12.6` installation from the remote tag loads
  GaussletBases `0.2.0` with the frozen tree; and
- the verification scratch and any generated manifest remain outside the
  repository, while shared `main`, the immutable tag, and both handoffs remain
  unchanged.

This manual lane replaces only the failed mechanical tag-CI verification. It
does not waive any identity check, rerun the numerical matrix, create a new
repository event, dispatch or alter a workflow, or authorize force-fetching,
moving, deleting, or recreating the tag. Immediately before publication,
repo-manager must recheck the accepted Docs/Pages state and release absence.
If the manual lane passes, the exact GitHub release authorized above may be
published. If it fails, preserve the tag and stop. Repairing the main-branch
tag lane for later releases is a separate workflow-maintenance decision and
does not block this exact release.

#### Final Publication Closeout

The namespaced fresh-scratch lane passed. It independently resolved annotated
tag object `722e8e8752a9d23f45e95d2f88e1749f9f3002e4`, exact message
`GaussletBases v0.2.0`, candidate `adfcaba32d4db06d9d796d947276433717bd2d89`,
tree `f64ba21e06ff57e2b5e78d91214398115afbe8de`, and package version `0.2.0`.
`git ls-remote`, the GitHub API, and a fresh Julia `1.12.6` installation from
`rev = "v0.2.0"` agreed. The install loaded GaussletBases `0.2.0` with the
frozen tree. No numerical gate was rerun.

GitHub release `378216554` is published at the immutable `v0.2.0` tag with
title `GaussletBases v0.2.0`, `draft=false`, `prerelease=false`, latest-final
status, zero uploaded assets, and the exact `2,278`-byte ASCII body above with
final newline and SHA-256
`e9ae9bcdad74b33bb66fb3e7e6a149d26285cb9bcc2f4c9555ac713be8bc90d2`.
Both automatic source archives reconstruct the frozen tree. `/v0.2.0/` and
the real `/stable/` alias are identical and canonicalize to `/v0.2.0/`;
`versions.js` retains stable, `v0.2`, RC2, RC1, and `dev`, and all prior
documentation folders remain intact.

The final tag and release are now immutable maintenance evidence. Do not edit,
delete, recreate, retarget, move, add assets to, or silently retry either one.
General registration, `CITATION.cff`, package/paper citations, paper titles or
journal status, large reproduction bundles, future tag-lane repair, and any
later patch release remain separate decisions.
