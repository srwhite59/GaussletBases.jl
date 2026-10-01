# Numerical Contracts

This page records small internal engineering policies that are easy for
generated code to miss. It is developer-facing rather than part of the public
manual.

## Radial Construction And Quadrature Convergence

Pass 622 accepts `a7ec77247b9bc45a349771e58cb76cb8019e08f8`, implementing
the Pass 621 repair under `HP-RADIAL-CONVERGENCE-FN-01/TEST-01`. The records
are implemented/completed, maintenance-only; no implementation grant remains.
The implementation boundary below records the completed grant. Construction
accuracy and quadrature convergence are separate tests, neither an energy-error
estimate. This supersedes the quiet exhausted-`:high` quadrature fallback,
not the existing refinement schedules or tolerances.

### Construction

In `src/foundation/bases.jl`, increase `_CONSTRUCTION_GRID_MAXITER` from four
to five. The shared `_select_construction_data` serves both radial and
half-line construction: preserve candidate order, half-grid spacing, the
`1e-6` overlap-deviation target, early exit at `deviation <= target`, and
best-candidate selection. Exhaustion returns that best candidate with one
warning carrying its achieved deviation and target. It must not report only
the final candidate when an earlier candidate was better. Explicit
`refine_grid_h=false` still builds once without evaluating refinement quality
or emitting this exhaustion warning. Preserve all other validation.

This is a numerical change: fixtures that previously exhausted can now have
different coefficients. Do not require old coefficient fingerprints for those
fixtures or weaken scientific tolerances to accommodate a failure. A selected
fifth candidate must agree exactly with direct construction at the same
`grid_h` and `refine_grid_h=false`. Already-converged inputs remain unchanged.

### Quadrature Notification

In `src/foundation/quadrature.jl`, exhausted `:high` refinement must return the
same last best-effort grid with one informative warning even when overlap
alone passes the looser public-quality fallback. Keep the `:medium` quiet
fallback exactly as before; `:veryhigh` and explicit-cutoff failures still
warn. Preserve the shortened-tail explanation for a cutoff below the
conservative tail bound. Do not add sampling or refine beyond the current
schedule to silence a warning.

An explicit positive `refine` is a starting hint, not a fixed grid or a
no-refinement switch. Preserve its doubling schedule and convergence checks;
exhausted explicit `:high` hints receive the same truthful warning. Preserve
`quadrature_rmax`, its legacy `rmax` spelling, mutual exclusion, and cutoff
validation. Successful refinement remains quiet. One call emits at most one
exhaustion warning; two independently requested grids may each warn.

The automatic high schedule stays `(24, 32, 48, 64)`. Its overlap, overlap
change, inverse-radius change, and moment-center change targets stay
`1e-6`, `1e-10`, `2e-10`, and `1e-12`, respectively. All medium/veryhigh
settings and the public `1.5e-5` overlap threshold remain unchanged.
Keep the actual penultimate metrics before replacing the latest metrics;
compute successive changes from these distinct grids, not latest minus itself.
Report the requested accuracy, starting/final refine values, unmet enabled
criteria, achieved measures and corresponding targets, cutoff, and tail bound.
These are transient logging values, not new stored metadata or a result type.
Replacing misleading `best_*` warning keys with truthful last-grid measures is
authorized; the returned grid and diagnostics APIs do not change.

Matrix changes use the existing maximum absolute entry norm. Inverse-radius
changes have units `bohr^-1`; moment-center changes have units bohr. They are
successive-grid diagnostics, not rigorous matrix-error bounds or energy-error
estimates. Existing `radial_quadrature(basis; accuracy=:high, refine=128)`
requests the schedule `128,256,512,1024`; it costs more and is not a universal
convergence guarantee. The construction opt-out does not disable quadrature
refinement. Update the two construction docstrings, quadrature docstring,
recommended atomic setup note, and radial algorithm's bounded-exit wording.

### Evidence And Limitation

Independent Julia 1.12.6 review reproduced both September 9 reports:
`radial-exits-2026-09-09/REPORT.md` SHA-256
`4c241350940bb12f57725daab7c3a5a4992dee9eef851815cce89d5b82051bac`
and `radial-fifth-attempt-2026-09-09/REPORT.md` SHA-256
`aa07189336b8e679fc1b5d0a3e9cc474b79426cad93dd9f16c84eb77b2b4add1`
under ignored `tmp/reviews/`. Larger sweeps remain evidence, not CI additions.
The supplemented Example 02 and recommended H fixtures improve construction
deviation from `2.313e-6` at `h=0.0025` to `3.542e-9` at `h=0.00125`, keeping
dimensions 9 and 35. Independent warmed construction costs were 0.406-0.466
to 0.804-0.900 seconds and 38.48 to 78.44 MB for Example 02; 1.970-2.030
to 4.049-4.063 seconds and 276.5-276.6 to 559.5-560.1 MB for H. These are
cumulative allocations, not peak memory. Accept this roughly doubled cost for
inputs needing attempt five, not a claim of free improvement.

On the improved fixed bases, default high quadrature still exhausts. Refining
64 to 512 changes the maximum inverse-radius entry by `0.0024226465 bohr^-1`
(Example 02) and `0.0036160909 bohr^-1` (H), both at `(1,1)`. For H that
entry changes from 486.6192516879 to 486.6228677788. The limitation is the
near-origin mixed boundary function, not a pure added Gaussian. The unchanged
default last-step changes are larger: 0.0048652696 and 0.0072618077
`bohr^-1`, respectively. These must not be falsely logged as zero. The
reported H construction energy shift is only about `-4.014e-13 Ha`; small
energy sensitivity does not establish convergence of every matrix entry.
The positive inverse-radius matrix is distinct from nuclear attraction `-Z V`.
Standard60, finite-expansion limitations, and other numerical findings remain
separate and held where previously held.

### Implementation Boundary And Acceptance

Only the two source owners above, `test/core/runtests.jl`,
`test/radial/runtests.jl`, `docs/src/howto/recommended_atomic_setup.md`, and
`docs/src/algorithms/radial_interval_sampled_build_and_extents.md` may change,
apart from required lifecycle records. The concrete scratch source proposal
is +38/-31 including docstrings (+21/-23 executable lines); preferred/hard
added-source budgets are 45/60. The 46-line scratch regression body passed
8 construction and 36 quadrature checks in about five seconds. Integrate it
compactly, replacing the two existing recommended-H quiet assertions with
expected single-warning assertions; estimated tests +48/-2, preferred/hard
55/70 added lines. Reader edits have preferred/hard additions 20/28. No new
file, helper framework, result type, stored metadata, preset, or adaptive policy.

Regressions must cover early success, nonmonotone best-on-exhaustion, explicit
construction opt-out, real fifth-candidate equivalence, quiet quadrature
success, exhausted automatic and explicit hints, actual penultimate measures,
enabled unmet criteria/targets, returned last grid, and shortened-tail warning.
Reuse existing core/radial owners for half-line, radial, energy, diagnostics,
and input coverage. Run these unchanged except the authorized additions and
two warning expectations, plus package load, docs, authority/self-test,
deterministic views, Documenter, log/diff checks, and normal full three-job CI
and Docs. Existing quick examples provide onboarding acceptance; do not add
duplicate expensive examples, the full angular suite, or audit sweeps.

Record fresh/warm construction cost and allocations, achieved deviations,
dimensions, and fixed-basis inverse-radius changes separately from energies.
If the repair needs other files, a changed schedule/tolerance, weakened
numerical assertions, broader machinery, or causes an unexplained material
regression beyond the measured fifth-attempt cost, make no implementation
commit and report. No export/API/default accuracy selection, dependency,
workflow, release, tag, or stable-documentation change is authorized.

## Orthonormal Blocks

When a block is constructed to be orthonormal, the intended contract is:

1. build it to be orthonormal in the relevant metric;
2. check that the resulting overlap is `I +` small unavoidable Float64 noise;
3. then treat the block as orthonormal.

Examples include finalized PGDG/COMX-cleaned blocks, nested fixed blocks, and
other internal blocks explicitly constructed as orthonormal bases.

Do not propagate a near-identity overlap matrix as meaningful working data. In
particular:

- do not store `S = I + eps` by default when `eps` is only Float64 residue;
- do not build downstream logic that repeatedly consults such an `S`;
- do not interpret `1e-12` to `1e-14` nonorthogonality as a structural feature.

Use self-overlaps only for construction, validation, assertions, or one final
cleanup when the deviation is too large. For transfer between final orthonormal
working bases, use the cross overlap only. Do not turn that path into a
generalized-overlap formulation.

For decomposed White--Lindsey plus a GTO supplement, the raw combined Galerkin
generalized solve is diagnostic. The working path projects the supplement into
the residual space, orthonormalizes retained residual directions, and solves an
ordinary Hermitian problem in the resulting final basis.

## Nested Fixed-Block Kinetic

`_NestedFixedBlock3D.kinetic` follows the nested packet contract:

- it is the kinetic matrix carried by the assembled nested packet;
- it is the kinetic payload downstream nested operators use;
- it is not automatically interchangeable with a later contraction of a
  separately assembled ordinary parent kinetic.

For current nested diatomic routes, these can differ measurably even on the same
final basis dimension.

## One-Body Reassembly

`assembled_one_body_hamiltonian(...)` reassembles from the payload actually
stored:

- `kinetic_one_body`;
- `nuclear_one_body_by_center`;
- requested `nuclear_charges`.

Compare reassembled one-body matrices only with operators built under the same
stored kinetic convention. Do not use this helper to claim equality between
routes that disagree about their kinetic payload.

## Current Cartesian Route Boundary

The active PQS/White--Lindsey producer realizes terminal blocks on disjoint
owned support and assembles exact one-body and IDA operators through current
terminal owners. Retained-unit and transform-contract records may feed that
construction, but the former unit-pair index, pair-operator-plan, source-safe
term, bridge/readiness, and pair-block-materialization ladder was
[retired](designs/cartesian_hamiltonian_producer/cartesian_pair_planning_materialization_retirement.md).

Old decomposed-WL pair-count ladders, timing tables, local block placement
results, and He/H2 acceptance numbers were development evidence for that
abandoned architecture. They are available in Git history, not current
numerical contracts.

The continuing rules are:

- shellification owns disjoint support;
- terminal realization owns shell-local rank, Lowdin, and sign conventions;
- cross-block overlap is structurally zero, but cross-block operators need not
  be zero;
- exact one-body and IDA matrices use the realized terminal basis directly;
- source/support summaries are diagnostics, never numerical authority;
- retained density normalization uses final function integrals at the explicit
  IDA boundary;
- uncharged nuclear matrices remain separated by center until physical charges
  are applied;
- no retired pair status, index table, placement plan, or materialization flag
  should be restored without a new source-backed consumer and docs-only
  authority.

The mapped-COMX axis-transform helper remains live at its existing path. Its
presence inside the reduced `CartesianPairBlockMaterialization` module is a
narrow ownership detail, not a surviving pair-materialization contract.

## Radial Multipole and Angular Interaction Repair

Pass 653 / `GB-ANGULAR-MULTIPOLE-20261001` accepted implementation
`35684b6ddeb0433f6f6e1d9ee324c9a6173d7cb9`. Both task-specific FN/TEST
execution grants are consumed (`completed`/`none`); the following records the
accepted contract and historical implementation boundary, not new authority.
The earlier policy describes the intended moment span, not proof that the
pre-repair assembler used it: that assembler stopped at the stored radial cap.
Reviewed candidates are packet patches 01, 02 and the improvements-on-02 diff;
00 duplicates 01+02. Neither patch 03 nor its high-l harmonic changes is granted.

**Target and implementation.** Repair artificial high-L radial underflow and
assemble the experimental mixed angular interaction over its full moment span.
Source edits are confined to `src/radial/operators.jl`,
`src/angular/angular_atomic_benchmark.jl`, and
`src/angular/angular_sequence_export.jl`; no new source file or export.
Replace global power scaling/recovery with adjacent-ratio prefix/suffix sums,
retaining the same quadrature formula, integral-weight normalization and
symmetrization. Require finite positive, nondecreasing radial points before the
recurrence; reject unsorted input rather than sorting it. Remove superseded
scaling helpers. Keep the useful benign raw-power test oracle. This removes
artificial global-scale underflow, not all possible floating-point limitations.

Keep `atomic_operators`' centrifugal range and default stored multipoles
`0:2*lmax`. Permit the reviewed nonnegative integer `multipole_lmax` keyword,
one private three-field sample owner (points, weights, basis values), and one
optional sample field on `RadialAtomicOperators`. Beyond-storage angular
multipoles use those exact retained samples without a new persistent cache.
Public `multipole(ops,L)` still bounds access to the stored range. Preserve
both inferred and parameterized eight-argument constructors with absent samples;
requesting unavailable on-demand data must fail clearly, not synthesize it.

Angular builders and the fixed-radial sequence use `interaction_lmax=:auto`
by default: include all L present in the assembly's interaction moment tables.
`:stored` selects the legacy stored cap; a nonnegative integer requests an
explicit cap, bounded by the available moment span. Reject other symbols,
negative integers and Boolean caps. `:stored` reproduces the cap convention,
not the old arithmetic bitwise. Exact radial/Ylm reference assembly stays on its
existing stored-range contract. Preserve the one-body, SCF and Lanczos policies.

Permit only the reviewed five-field compact interaction plan (used, required,
stored caps, truncation, mode), carried on the HF-style benchmark and HF adapter,
and surfaced in their existing diagnostics and fixed-radial level metadata.
Preserve inferred/parameterized seven-argument HF-style and thirteen-argument
adapter constructors with explicitly unknown provenance. Unknown cap or
truncation is `missing` in memory and an explicit `"unknown"` in dense metadata,
never zero, false or an inferred enforced cap. Known facts retain their current
integer/Boolean representation. No broader carrier or schema redesign.

For the one-body adapter overload accepting an external interaction, a `nothing`
request sentinel may distinguish omission from an explicitly requested cap.
Omission resolves to auto only when assembling an interaction. Reject an
explicit cap alongside an external matrix: this path cannot enforce that cap.
Reuse a carried plan only for the exact HF-style one-body/interaction objects
forwarded by the HF-style adapter; an unrelated `hf_style` must not attest a
caller's matrix. Otherwise record external/unknown provenance. Preserve existing
matrix-size, occupations, seeds, route and solver-mode checks.

Include the normalized cap request in fixed-radial sequence and level identities;
preserve shell/profile/gauge identities and adjacent versus direct sidecars.
Retain `manifest/source/multipole_Lmax` as stored range, adding only the reviewed
interaction used/required/mode/truncation metadata. Arithmetic can change the
radial multipole checksum and hence radial/sequence/level IDs; receipts require
owner review, not golden repinning. Existing dense/JLD2 payload arrays remain
the interchange boundary. No binary-carrier migration framework is granted;
stop if a current consumer needs old whole-carrier deserialization support.

**Tests and budgets.** Edit only `test/radial/runtests.jl`,
`test/angular/runtests.jl`, `test/runtests.jl` (small-ED fixture keyword only),
`test/driver_public/angular_fixed_radial_sequence_runtests.jl`, and
`test/docs/runtests.jl` (changelog placement only). The docs check must preserve
the exact Changelog header, permit an optional leading Unreleased section, and
require the first versioned release section to be v0.2.1; all other release and
version checks remain unchanged. No loose presence check or new helper.
Source additions including docstrings/compatibility: preferred 300, hard 330;
per-file ceilings radial 160, benchmark 160, sequence 30, subject to total 330.
Test additions: preferred 220, hard 250; radial 100, angular 114, runner 5,
angular_public 25, docs 6 (at most 3 deleted docs-test lines). Reader additions:
preferred 75, hard 85, only CHANGELOG,
the angular interaction note, radial stabilization milestone and
`docs/src/explanations/angular_research_track.md`. Necessary authority digest
and deterministic-view reconciliation is permitted, not lifecycle/scope changes.

Independent direct-pairwise or extended-precision quadrature oracles must cover
benign and high-L input. Add compact production-like erf-grid coverage at
L=17,18,24 that fails baseline; check nonzero accurate values, not only symmetry
or finiteness. Use normwise relative error at most 1e-12 with meaningful entry
checks, and sorted-point rejection. Reassemble angular interactions independently
from every moment block and independent radial kernels: auto agrees within
1e-12 relative/scaled absolute error. Compare legacy cap within roundoff, not
bitwise equality with old source. Cover explicit caps, unavailable samples,
legacy constructors, external unknown status and external-plus-cap rejection.
Correct both small-ED fixture and its adapter comparison to explicit `:stored`;
preserve the energy anchor and solver settings, not an auto golden repin.
Check cap-dependent level/sequence IDs and one metadata round trip.

Run existing radial, angular, ida, core, angular_public, docs_fast and misc
owners unchanged apart from these bounded tests. Complete angular runs once on
the exact final candidate; reuse that evidence during closeout. Record timings
and allocations for representative kernel/operator assembly without an end-to-end
speedup claim. Require package/resource load, full docs, authority/self-test,
two matching external generated views, Documenter, manager-log/diff gates, and
exact-head full source-bearing CI (all three numerical jobs) plus separate Docs.
Keep fast, CI classifications/rows, optional HFDMRG behavior and all gates intact.

Changelog distinguishes arithmetic repair from the intentionally changed angular
default and documents ID/receipt consequences. Do not claim general atomic or
molecular accuracy: quadrature bias and high-l spherical-harmonic limitations
remain. Ne is closed-shell with occupied p orbitals; its measured proxy shift
does not quantify Cr2. After acceptance, make a bounded read-only consumer
inventory using actual stored cap/multipole manifests where available; classify
unknowns and owner-reviewed rebuild/receipt needs, preserving all old artifacts.
Do not rebuild or alter consumers or other repositories.

**Failure and endpoint.** Stop without an implementation commit for scope/budget
expansion, unresolved compatibility/science, failed scientific gates or required
golden/tolerance relaxation. At most two in-scope correction rounds. Exclude
quadrature/spacing/kink changes, one-body inverse-r2/lmax redesign, adaptive-L
optimization, new frameworks, releases, H10/Hchain and any successor. Endpoint
is independently accepted A+B repair, green checks and lifecycle closeout
consuming implementation grants; adviser notification only, consumers stay paused.

Acceptance reused the exact once-only angular evidence and green source-bearing
CI/Docs; independent direct/256-bit radial checks pass at L17/18/24. A bounded
read-only consumer inventory confirms one old Be15 level with stored cap4 and
absent actual-cap provenance. Other inspected manifests/configurations require
owner review; do not infer complete caps or quantified energy changes. Preserve
old artifacts and review rebuilding/identity receipts as coherent families.
High-l harmonics and quadrature bias remain separate; no successor is granted.
