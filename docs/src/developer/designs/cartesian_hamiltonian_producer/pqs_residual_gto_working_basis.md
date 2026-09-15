# PQS Residual-GTO Same-Construction Working Basis

## Status And Authority

`HP-REP-PQS-RG-WORKING-FN-01` and
`HP-REP-PQS-RG-WORKING-TEST-01` own the implemented narrow expert construction
and its completed bounded validation. Commit `346589a6d` is the accepted source
implementation; both records now grant maintenance only.

This contract closes one representation boundary. It does not change residual
selection, one-body or interaction physics, external-packet semantics, or the
opaque `CartesianIDAHamiltonian` contract.

Pass 631 separately authorizes the [finite-collinear connection](#Finite-Collinear-Supplementation)
under HP-COLLINEAR-PQS-RG-FN-01/TEST-01. The original atom/diatomic contract
below remains unchanged; only that new overload returns accumulated matrices.

## Consumer Need

Before this boundary was implemented,
`cartesian_residual_gto_mwg_hamiltonian(...)` constructed a PQS terminal basis
`G`, a Gaussian supplement `A`, retained residuals `R`, exact augmented
one-body matrices, and the MWG interaction in one pass, then returned only
`CartesianIDAHamiltonian`. The implemented system constructor retains the
compact same-construction data needed to overlap that final basis with later
GTO probes.

The external importer correctly requires

```text
C_B = S_BX * C_X,
S_BX = <B|X>,
```

for an orthonormal final basis `B` and explicit external GTO basis `X`.
Reconstructing a representation from an opaque Hamiltonian is impossible and
is not authorized.

## Public Boundary

The package exports exactly one expert root constructor:

```julia
cartesian_residual_gto_mwg_system(
    system::NamedTuple;
    basis::NamedTuple,
    supplement::NamedTuple,
)
```

The call returns one concrete in-memory result from the same construction as
the Hamiltonian. The concrete result type remains unexported. Its supported
consumer surface is:

```julia
result.hamiltonian
gto_overlap_matrix(result, probes; block_indices = nothing)
import_external_gto_orbitals(result, packet)
```

For this NamedTuple-input overload, `result.hamiltonian` is the existing
`CartesianIDAHamiltonian{Float64}`. The
result retains only the terminal realization, factorized parent-axis donor,
explicit supplement representation, and existing residual object needed by
the overlap action. It must not retain the full base stage, operator stages,
or a dense parent-by-terminal coefficient matrix.

The existing non-exported
`cartesian_residual_gto_mwg_hamiltonian(...; hamfile)` call keeps its exact
signature, direct Hamiltonian return, artifact behavior, and numerical output.
Both public and compatibility calls must share one private construction rather
than rebuild independently.

No `basis_representation(result)` method is required or approved. This object
supports exact GTO cross overlap and external-orbital import; it is not a new
general Cartesian representation family.

## Cross-Overlap Contract

For terminal basis `G`, explicit supplement `A`, and retained residuals

```text
R = G*T_G + A*T_A,
B = [G,R],
```

reuse the existing factorized parent/GTO and contracted-GTO overlap kernels to
form

```text
S_GX = <G|X>,
S_AX = <A|X>,
S_RX = T_G' * S_GX + T_A' * S_AX,
S_BX = [S_GX; S_RX].
```

Rows are returned in native `[G,R]` Hamiltonian order. Probe columns retain
the packet/probe AO order. Existing `block_indices` ordering and bounds
semantics remain unchanged.

The terminal projection must use the support-blockwise terminal realization
and factorized parent data. The unavoidable final-by-probe output and a
parent-by-probe analytic temporary are allowed. A dense parent-by-terminal
map, dense final self-overlap, generalized final metric, or assembled
`[I T_G; 0 T_A]` raw-to-final matrix is forbidden.

Before returning or applying an overlap, validate:

- Hamiltonian dimension equals `nG + nR`;
- the residual base dimension equals `nG`;
- supplement orbital count equals the residual candidate count;
- `T_G` and `T_A` have exact `(nG,nR)` and `(nA,nR)` shapes;
- the public constructor accepts no separately supplied terminal, supplement,
  transform, representation, or Hamiltonian component; all retained pieces are
  captured inside the one producer call;
- requested probes and all computed cross-overlap entries are finite.

The result stores no external packet, imported orbitals, density, screening
field, solver state, or consumer policy.

## Implementation Boundary

Approved source owners are limited to:

- `src/GaussletBases.jl` for the one export;
- `src/cartesian/cartesian_base_hamiltonian.jl` for the compact result and shared
  same-construction assembly;
- `src/cartesian/cartesian_gto_probes.jl` for the exact cross-overlap dispatch;
- `docs/src/reference/export.md` for the exported function's doc entry.

Reuse the existing terminal projection, supplement-overlap, residual, and
external-import owners. Add no source file, module, exported type, accessor,
adapter, provider hierarchy, metadata key, status symbol, or artifact field.

Source budget, including the export and docstring:

- preferred: at most 110 added source lines;
- hard stop: at most 150 added source lines.

These are stop-and-report bounds, not reasons to obscure validation or merge
unrelated responsibilities. If coherent implementation exceeds the hard
bound, make no source commit and return the exact missing reusable operation.

## Accepted Implementation Validation

The original implementation validation used:

- `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` for exact
  same-construction, formula, import, artifact, and provenance evidence;
- `test/docs/runtests.jl` for the single export/reference parity entry.

That supplemented H2 gate caught the information-loss and transpose/order
defects that the older Hamiltonian tests could not detect:

1. Construct the result through the new public function and preserve exact
   Hamiltonian matrices, dimensions, due diligence, and direct-facade behavior.
2. Compare `gto_overlap_matrix(result, probes)` with a direct assembly of
   `S_GX`, `S_AX`, and `T_G'*S_GX + T_A'*S_AX`.
3. Verify native row order, `block_indices`, finiteness, and dimensions.
4. Build a bounded valid `ExternalGTOOrbitalPacket`; verify the unchanged
   importer returns exactly `S_BX*C_X` and its existing capture diagnostics.
5. Reject malformed probe identity, dimensions, nonfinite data, and any
   dimension/identity inconsistency at the private assembly boundary before
   import.
6. Inspect and report terminal due diligence for the endpoint.

Accepted implementation evidence from commit `346589a6d` is:

- the H2 basis has `487` terminal plus `18` residual functions, for native
  dimension `505`;
- the system and direct facades produce exactly equal Hamiltonian matrices,
  and artifact readback deltas remain zero;
- direct assembly and the implemented `S_GX`/`S_AX`/`S_RX` formula agree
  exactly, as do imported coefficients and `S_BX*C_X`;
- the imported occupied capture is within `1e-8` of unity;
- the complete existing H2 owner passed `464/464`, and the supplemented
  facade/import owner passed `92/92`;
- terminal due diligence remains `9x9x15` parent axes, approximately
  `+/-4.871`, `+/-4.871`, and `+/-6.110` bounds, `4`-bohr padding, retained
  rows `[275,114,98]`, and only the pre-existing shape, axis-display, and
  large-identity warnings.

## Public Cartesian CI Extraction

`HP-REP-PQS-RG-WORKING-CI-FN-01` owns the implemented public-CI wiring.
`HP-REP-PQS-RG-WORKING-CI-TEST-01` owns the completed public validation under
maintenance. Commit `50b5ee9a1` moved the original `21`-check contract out of
the unwired, mixed-purpose R3A suite; commit `c633d5db0` expanded that owner to
`47` checks while deleting the redundant general-import nested testset. Neither
commit changed production behavior.

Maintain exactly one narrow public owner:

```text
test/driver_public/cartesian_residual_gto_mwg_system_runtests.jl
```

Keep that owner included once in the existing `:cartesian` group in
`test/runtests.jl`. The existing Julia `1.10` Supported-floor row owns the
remote regression. Do not add a CI row, alter the PQS or screening paper
groups, or wire the complete R3A file or any other `test/nested` suite.

The new owner exercises only documented public behavior:

1. Construct the bounded H2 system with the accepted `q=5`, `0.5` core
   spacing, `6.0/4.0` extents, and contracted H/cc-pVTZ `lmax=1` supplement.
2. Treat the concrete result as opaque except for `result.hamiltonian`.
   Validate native dimension `487 + 18 = 505`, finite symmetric matrices,
   the `1/1` particle sector, nuclei, and the existing
   `0.4574161883692301` self-Coulomb fingerprint.
3. Build a normalized one-orbital external GTO probe entirely through public
   representation and packet constructors. Exercise full and indexed
   `gto_overlap_matrix`, either restricted or alpha/beta import, source
   orthogonality, capture, and stale packet fingerprint rejection. Do not call
   `_cartesian_supplement_cross_overlap`, inspect private result fields, or
   reconstruct the residual transform formula in this owner.
4. Reject a small representative set of malformed public `system`, `basis`,
   and `supplement` inputs. Do not reproduce every historical blocked-path
   assertion.

The old R3A owner must lose the public-facade assertions duplicated by this
gate. Preserve there only unique private exact-formula, assembly-boundary,
artifact/readback, provenance, and due-diligence evidence until that suite has
its separate maintenance-versus-quarantine review. The existing
`HP-REP-PQS-RG-WORKING-TEST-01` continues to own those private oracle and docs
checks; the new CI records own only the public extraction and runner wiring.

The accepted owner is now `136` lines. The original extraction added `70` and
deleted `75` tracked test lines. The external-import extraction then added and
deleted `75/7` lines in the public owner and deleted `111` nested lines. Its
total test delta is `+75/-118`, net `-43`. Maintenance must not expand the
fixture, restore duplicate nested assertions, or add a standalone fixture or
helper file.

Accepted evidence from commit `50b5ee9a1` is:

- the focused public owner passed `21/21` and the complete `:cartesian` group
  passed `253/253`;
- the retained nested owners passed `464/464` plus `64/64` private checks;
- remote CI run `32613923516` passed Supported floor, PQS paper, and Screening
  paper without changing either paper gate;
- terminal due diligence remains `9x9x15` parent axes, `4`-bohr padding, `487`
  terminal functions, `18` residuals, and final dimension `505`.

### Accepted External-GTO Public Coverage Extraction

Commit `c633d5db0` moved the unique general-import assertions from the first
testset in `test/nested/cartesian_external_gto_import_runtests.jl` into the
existing public owner. It reuses the already-constructed residual-GTO system;
no second Cartesian construction was added.

The accepted fixture selects two normalized, mutually orthogonal matched
`px`/`py` Cartesian Gaussian probes through root-public representation APIs.
Their common center, exponent, contraction, normalization, and angular parity
make `S_GG = I` a property of the fixture. The public owner calls no private
helper and inspects no internal representation.

The extracted assertions are limited to:

- one basic successful import;
- equal-occupation occupied-space rotation invariance;
- explicit alpha/beta import and invalid spin combinations;
- valid source-metric and ordering identity;
- stale ordering or fingerprint rejection;
- refingerprinted-but-wrong and nonsymmetric source metrics;
- rejection of nonorthonormal source coefficients.

The complete redundant `External GTO orbital import` testset is deleted from
the nested file. Its `Protected external GTO representation sidecar` testset,
artifact/tamper cases, direct-run status, and required shared helpers remain
unchanged. No duplicate restricted-import cluster remains.

Acceptance passed `47/47` public checks, `49/49` protected-sidecar checks, and
`279/279` in the wired `:cartesian` group. Terminal due diligence remained
`9x9x15` parent axes, approximately `+/-4.871`, `+/-4.871`, and `+/-6.110`
bounds, `4`-bohr padding, retained rows `[275,114,98]`, and native dimension
`487 + 18 = 505`, with unchanged advisory warnings. No test file, runner,
workflow, CI row, source/API surface, dependency, or numerical policy changed.

With repository acceptance complete, the separate REQ-101 consumer may rerun
only its frozen `R=2.35`, `ns=5`, PQS early case. It must use the unchanged
external
determinant and supplement, the public packet/importer, source metric
`C_X' S_XX C_X`, Euclidean final capture, the existing `1e-8` capture bounds,
and worst occupied loss at most `1e-4`. This is external acceptance evidence,
not a committed C2 fixture or a repo-owned paper campaign.

## Exclusions And Failure Rule

Do not change `CartesianIDAHamiltonian`, residual selection or cutoff, exact
one-body assembly, MWG/IDA interaction, external packet/importer behavior,
screened-Hartree algebra, public basis inputs, artifacts, drivers, solvers, or
defaults. Add no PySCF dependency or parser, C2/REQ-101 branch, molecular
screening automation, generalized overlap, final-final transfer, Hamiltonian
transform, interaction transform, persistent sidecar, or dense representation.

If the exact cross overlap cannot be produced from the same construction with
the existing factorized and supplement kernels, or if direct-facade parity,
the H2 gate, metric identities, or resource boundary fails, make no source
commit. Report the exact missing seam; do not reconstruct from the Hamiltonian,
fall back to a private consumer map, or add partial scaffolding.

## Finite Collinear Supplementation

Pass 631 authorizes HP-COLLINEAR-PQS-RG-FN-01/TEST-01 at
4cf1c7b125f6a8f11f3559f835eab24e70798f89. Named consumer hchain-doer remains
paused until implementation is independently accepted. All scientific H-chain
calculations must be supplemented; bare matrices here are intermediate data,
not permission for a bare scientific campaign.

### Target and interface

Target: connect the existing finite-collinear basis to residual-GTO selection,
complete supplemented H1/MWG operators and raw transfer. Physics endpoint:
small H3 and H-He-H supplemented matrices, not H10 accuracy or scaling.
Add exactly this overload, with both keywords required:

```julia
cartesian_residual_gto_mwg_system(working::_CartesianCollinearWorkingBasis,
    z, Z; supplement::CartesianGaussianShellSupplementRepresentation3D,
    expansion::CoulombGaussianExpansion)
```

Keep the spelling private to the existing opaque working type; add no export.
Reuse the existing five-field _CartesianResidualGTOMWGSystem, adding a concrete
Hamiltonian type parameter, not an Any field or new result type. Its new
.hamiltonian is exactly (; one_body, electron_electron_ida, nuclear_repulsion):
two Matrix{Float64} values and a Float64 scalar. Validate only this exact
three-field representation and CartesianIDAHamiltonian{Float64}, not arbitrary
duck-typed objects or generic NamedTuples. Preserve old dimension/error behavior;
check both new matrices' final square dimensions and finite values.
Reuse _cartesian_residual_gto_mwg_system_result and all same-construction checks.

The existing NamedTuple-input overload, two-center guard, construction and
CartesianIDAHamiltonian return, electron-sector and artifact/reweighting
semantics remain unchanged. The new overload provides none of those latter
semantics. No determinant cleanup is added or implied. Cross overlap and raw
import bodies remain byte-identical; the shared validator serves both carriers.

Use existing z/Z validation. Supplied ordered finite z and positive Z define
the potential AND candidate-center ownership, not the original basis geometry.
Do not compare them with the basis-construction nuclei or rebuild the basis.
Each candidate must exactly match one supplied center. Require nonempty,
consistent primitive arrays, finite coefficients/centers, positive finite
exponents, nonnegative Cartesian powers and supported axiswise normalization.
Support contracted candidates, not only one-primitive orbitals. Preserve existing
rank/metric failures, including no surviving residual direction; no bare fallback.

### Numerical and storage invariants

Use existing residual selection with occupation cutoff 1e-6, injection disabled,
negative/merge absolute and relative thresholds 1e-12, cross orthogonality 1e-10,
identity tolerance 5e-8. Expose no new threshold or approximation control.
Expansion/parent exponents must match under the existing validator.

Reuse the non-nuclear raw blocks and augmented kinetic/moment assembly.
The unused parent_basis_object argument may receive nothing; create no old
parent-stage wrapper. Process one nucleus at a time through the existing raw
nuclear kernel. Accumulate charged terminal/candidate and candidate/candidate
blocks, then apply the existing residual transform. Base H1 already includes
the potential; do not double-count it. Build MWG through the existing assembler.
Do not retain per-nucleus final matrices or a global parent-to-final map.
One temporary parent-by-candidate block, candidate self block and existing
moment workspace are allowed. Final data retain only existing five carrier
fields, with two complete operator matrices and scalar nuclear repulsion.

Preserve MWG w=sqrt(2*(<x^2>-<x>^2)), Gaussian exp(-(x-c)^2/(2*w^2)),
division by its integral, and terminal integral-weighted projection. Independent
R-R convolution uses 1+2*t*(w_i^2+w_j^2), not squared Gaussians. Preserve the
base G-G IDA block exactly. MWG is approximate two-index interaction, not
four-index ERIs or exact represented-density Hartree.

### Evidence and frozen acceptance

Reviewed report tmp/reviews/collinear-supplement-connection-2026-09-14/REPORT.md
SHA-256 7694dcd9a6cd0c44d6789c057072b05fdacef1b700d44d425dd9b72585b6cc08;
probe.jl SHA-256 c3f4dcd482ef487b02b5a2403cab073404e6c4639d978c71ebfe259cb39fd4f4.
The corrected independent oracle, not its initial squared-Gaussian mistake,
is evidence. No production change follows from that oracle correction.

Reuse compact45, spacing .6, padding 3/3, core_side=3, angular reference 3,
outer count 3, tail_spacing=2.8 and angular scale 1.4. H3 positions [-1.2,0,1.2],
charges [1,1,1]; H-He-H [-2.4,0,2.4], charges [1,2,1]. Normalized s/px/py/pz
candidates at each center with exponent .8 retain 12 residuals, four per owner.
Base/final dimensions 231/243 and 223/235; repulsion 2.0833333333333335/1.875 Ha.

Freeze max-entry absolute tolerances, rtol=0: actual complete supplemented metric
against identity <=5e-8; H1/MWG symmetry and complete H1 block-congruence <=1e-10;
independent AA/all-pair and selected GA overlap <=1e-12, H1 <=1e-10;
independent all R-R and selected all-row G-R MWG <=1e-10; repulsion <=1e-12 Ha.
Preserve all-row GA tests for representative s/px/py/pz and all 144 AA/R-R pairs.
Reuse compact GH/local-support and Gaussian-convolution oracles, not copies of
production contractions. Raw-block observed errors (~1e-15 to 1e-14) are distinct
from complete transformed H1 (~2.1e-12/5.8e-12) and metric (~3.2e-12/7.1e-12).
Do not advertise raw-block precision for complete transformed matrices.

Raw import equals cross overlap times input coefficients exactly; cross overlap
agrees with the independent augmented transform within 1e-10. The normalized
primitive px/py capture deviation <=5e-8 is only this fixture, not certification
of diffuse virtual-space convergence. Preserve separate px and py results.
Add one compact nontrivial contracted-candidate case using the shared
primitive-sum oracle and the same tolerances. Existing contracted cc-pVTZ
atom/diatomic evidence is reusable but does not alone test the new connection.
Add one case using a fixed working basis with different valid potential z/Z
and matching candidate centers; compare against the same independent blocks.
Do not require every nucleus to supply candidates if the existing selection
allows a subset. Check invalid nuclear data, off-center/invalid candidates,
mismatched expansion, failed residual selection and unsupported carrier shapes.

Reuse existing atom/diatomic residual-GTO owner and transfer fixture, finite-chain
public owner, relevant atomic/diatomic regression endpoints, package load,
authority/self-test/generated parity, docs and Documenter. Full source-bearing
CI must execute all three unchanged numerical jobs plus Docs. Do not repeat
the complete angular suite or duplicate Example 41's shared release calculation.
No test may depend on ignored scratch files. Inspect actual terminal bounds,
counts, ownership, residual counts, metric and warnings, not invented diagnostics.

Record one fresh/warm bounded fixture cost separately from reference checks.
Evidence: H3 supplemented assembly 6.354/.111 s, 1.502 GB/114 MB allocated;
retained 2.10 MB, H-He-H warm .071 s/98.6 MB, retained 1.93 MB. Keep final
dimension below 600 and retained storage below 512 MiB for these fixtures.
Investigate material regressions; no end-to-end or long-chain speed claim.
Parent-by-candidate storage and dense final quadratic scaling remain limitations.

### Exact scope, budgets and stop rule

Source: only src/cartesian/cartesian_base_hamiltonian.jl, preferred/hard
100/130 added lines including docstrings and relocated code.
Tests: only test/driver_public/cartesian_residual_gto_mwg_system_runtests.jl,
190/220 added lines. Shared compact oracles should cover primitive and contracted
cases without duplicated test bodies; preserve existing tests.
Reader guidance: docs/src/manual/projected_q_shells.md (20/30) and
docs/src/reference/export.md (5/10), total 25/40 added lines.
Normal authority/current/generated/log updates are documentation-only.
No new tracked file, public binding/type, field, helper framework, cache,
metadata, kernel, parser, dependency, workflow, approximation, threshold,
solver, screening, H10/H20 campaign, release, stable promotion or retirement.

Must simplify: connect through existing numerical owners and stream nuclear
blocks rather than adapt the old staged frontend. No deletion of live old
facades, ordinary/sliced capabilities or recursive code is authorized.
Private file-local composition/validation is allowed within the single source
budget; no shared numerical helper or new carrier is needed.

Failure rule: stop without an implementation commit if budgets, frozen
tolerances, rank rejection, old behavior, storage or two-representation validation
cannot be preserved, or another file/kernel/semantic change is needed.
Do not relax tolerances, clamp, silently fall back, add compatibility machinery
or release hchain-doer. Repo-manager begins only after this authority commit
and required checks pass; independent implementation acceptance is a separate
gate before any consumer resumption. No scientific campaign is authorized here.
