# Numerical Contracts

This page records small internal engineering policies that are easy for
generated code to miss. It is developer-facing rather than part of the public
manual.

## Radial Construction And Quadrature Convergence

Pass 621 authorizes the following bounded repair under
`HP-RADIAL-CONVERGENCE-FN-01/TEST-01`; implementation is pending. Construction
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
