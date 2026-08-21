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
v0.2.0. Exact candidate-tag and RC1 GitHub-prerelease authority are isolated
below; final-release publication remains separate. The existing numerical
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

## Post-RC1 Export Integrity Repair

`HP-PUBLIC-EXPORT-INTEGRITY-FN-01/TEST-01` authorize one bounded cleanup of
package export declarations discovered after the RC1 release. At baseline
`41ab3e13121f5af1e145775500e91f9ac61c9760`, runtime inspection of the root
module and every direct package-owned child module finds exactly four defects:

- root export `CartesianBasisBundle3D` is undefined;
- root export `nested_fixed_block_timing_report` names a bare generic with no
  methods;
- internal `CartesianRouteCore` exports `final_units` and `unit_keys`, neither
  of which is defined by that module.

No committed source, test, driver, tool, or example calls these names. The
implementation must remove the two root export entries, delete the bare
`nested_fixed_block_timing_report` generic declaration, and remove the two
internal route-core export entries. It must add no replacement definition,
alias, shim, deprecation, helper, file, dependency, metadata, or numerical
behavior. Julia does not permit the remaining final route-core export to keep
a trailing comma, so the exact coherent source delta is `+1/-6`, net `-5`:
the sole added line is existing export `lowering_recipe` without that comma.
The hard added-source limit is one line. If a live caller or implementation is
found, stop without committing rather than preserve an invalid surface through
compatibility glue.

One regression belongs in existing owner `test/core/runtests.jl`. It must audit
`GaussletBases` plus each defined direct child module discovered with
`names(GaussletBases; all=true, imported=false)` where
`parentmodule(child) === GaussletBases`. For every exported name
returned by `names(module; all=false, imported=false)`, require a defined
binding; if that binding is a `Function`, require at least one method. This
checks export integrity without freezing an export inventory or traversing
external/imported modules. The preferred added-test budget is `20` lines and
the hard limit is `30`; no new test file or fixture is authorized.

This repair deliberately leaves `@timeg`, `CuratedSpherePointSet`,
`LegacySGaussianData`, `QiuWhiteHybridOrbital3D`, and
`QiuWhiteResidualGaussianOperators` unchanged. Their compatibility or use
policy is not decided here. Package version, tags, GitHub releases, workflows,
and the immutable RC1 tag/release state also remain unchanged. After
implementation and lifecycle closeout, return control; paper-aligned
PQS-versus-screening CI/release wiring is a separate design boundary, not a
continuation of this cleanup.

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
