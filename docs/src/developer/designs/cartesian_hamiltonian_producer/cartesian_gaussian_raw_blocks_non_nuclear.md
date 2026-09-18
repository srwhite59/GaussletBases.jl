# Cartesian Gaussian Raw Blocks - Non-Nuclear Slice

This page owns the implemented neutral `G-A` and `A-A` overlap, kinetic,
coordinate-moment, and second-moment raw blocks. It does not own residual
selection, final-basis transforms, terminal `G-G` products, or Qiu-White
provider semantics.

## Lifecycle

| ID | Lifecycle | Current boundary |
| --- | --- | --- |
| `HP-CGRB-NN-FILE-01` | Implemented | Non-nuclear owner file and module include |
| `HP-CGRB-NN-FN-01` | Implemented | Full and overlap-only neutral raw blocks |
| `HP-CGRB-NN-WIRE-01` | Implemented | Main diatomic Residual-Gaussian/Qiu-White callers use the neutral owner |
| `HP-CGRB-NN-TEST-01` | Validation completed; tracked coverage is indirect | H2 endpoint plus accepted QW/Be2/Cr2 parity evidence |

Implementation commits are `00d052e29` for extraction, `9fa0cc16d` for the
overlap-only path, `806b37e32` for `A-A` family reuse, and `71a89433c` for
`G-A` family reuse. Manager-log Passes 084B-086B record acceptance.

## Ownership And Callers

The owner is:

```text
src/cartesian/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl
```

It is loaded by
`src/cartesian/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl`. Shared
private axis families come from the same internal module, and analytic
overlap/kinetic/position/x2 tables come from
`src/cartesian/cartesian_gaussian_axis_integrals.jl`.

Implemented entry points are:

```julia
gaussian_non_nuclear_raw_blocks(proxy, supplement, expansion)
gaussian_non_nuclear_overlap_blocks(proxy, supplement)
```

Direct live callers are in:

- `src/cartesian/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`;
- `src/ordinary/ordinary_qw_raw_blocks.jl`.

The full helper supplies exact augmented operators and the main diatomic
Qiu-White wrapper. The overlap-only helper supplies the residual setup mixed
overlap without constructing unused kinetic, moment, or `A-A` blocks.

## Returned Blocks

Let `nG = proxy.ncart` and `nA = length(supplement.orbitals)`. The full helper
returns exactly:

```text
ga.overlap       :: nG x nA
ga.kinetic       :: nG x nA
ga.position.x/y/z:: nG x nA
ga.x2.x/y/z      :: nG x nA

aa.overlap       :: nA x nA
aa.kinetic       :: nA x nA
aa.position.x/y/z:: nA x nA
aa.x2.x/y/z      :: nA x nA
```

The overlap-only helper returns `(; ga = (; overlap))` with one `nG x nA`
matrix. Columns and `A-A` rows follow `supplement.orbitals` order exactly.
There is no sorting or owner regrouping. All `A-A` outputs are explicitly
symmetric; `G-A` outputs are rectangular.

The full helper retains an `expansion` argument for caller-signature parity,
but non-nuclear values do not depend on Coulomb expansion terms.

## Numerical Convention

One-dimensional tables use the shared private axis helper terms:

```text
:overlap
:kinetic
:position
:x2
```

The existing Gaussian polynomial prefactors and analytic-integral sign and
normalization conventions are unchanged. Three-dimensional primitive products
are coupled across axes. Kinetic is the sum of the three terms with one
kinetic axis and two overlap axes; each coordinate or second moment replaces
the corresponding axis overlap table only.

The full helper builds one supplement axis-family inventory keyed by exponents,
axis center, Cartesian power, and axis prefactors. Its `G-A` tables are reused
by family, while `A-A` tables use canonical family-pair IDs and reversal flags.
The overlap-only helper skips every unused operator and all `A-A` work. For an
off-diagonal `A-A` orbital pair, both requested orientations are evaluated and
averaged before the value is mirrored; this preserves the prior Qiu-White
roundoff/symmetry convention.

All tables, maps, and scratch matrices are function-local. The owner returns
raw matrices only, with no status, route, provider, cache, or artifact object.

Accepted Cr2 q4 evidence reduced full non-nuclear raw-block allocation from
about `10.9 GiB` to about `0.86 GiB`, with block parity at roundoff. The
overlap-only residual setup fell from about `10.99 GiB` to about `199 MiB`.
These are compact implementation evidence, not Cr2-specific behavior.

## Consumer Boundaries

`CartesianResidualGaussians` owns residual selection, basis orientation,
augmented transformations, moment-matched descriptors, and MWG interaction.
It consumes these raw matrices but does not own their analytic construction.

The main diatomic Qiu-White path consumes the neutral owner. The following
Qiu-White helpers remain live and must not be deleted under `HP-CGRB-NN-*`:

```text
_qwrg_cartesian_shell_cross_moment_blocks_3d
_qwrg_cartesian_shell_self_moment_blocks_3d
_qwrg_atomic_cartesian_blocks_3d
```

They still serve atomic Qiu-White reference/operator assembly, hybrid
representation sidecars, and dense-parent GTO probes in
`ordinary_qw_operator_assembly.jl`, `cartesian_qw_hybrid_representation.jl`,
and `cartesian_gto_probes.jl`. The orphaned
`CartesianCPBBlockProviders.jl` caller was removed in `181ed6968` and is not
part of this ownership contract. Their
`factor_ga`/`factor_aa` outputs are outside this overlap/kinetic/moment
contract. Rewiring or deleting them requires a separate caller/ownership audit.

Terminal final-basis `G-G` product matrices are separately owned by the R3
terminal optimization contract. Mixed-Hartree `GG/GA/AA` reference blocks are
separately owned by the reference-Hartree contract. Neither belongs here.

## Validation And Failure Boundary

Tracked endpoint coverage is:

```text
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
```

It exercises overlap setup and full raw blocks through residual construction
and exact augmented kinetic/position/x2 operators. Accepted ignored
Qiu-White, Be2, and Cr2 parity evidence remains in manager-log Passes 084B-086B;
no dedicated committed
`cartesian_gaussian_raw_blocks_non_nuclear_runtests.jl` file exists.

The internal helpers assume validated proxy/supplement data. Unsupported axis
terms, inconsistent primitive arrays, and linear-algebra dimension errors
throw or propagate; there is no fallback result.

This owner does not authorize nuclear or mixed-Hartree changes, final-basis
`G-G` optimization, residual selection/transforms, Qiu-White semantic changes,
factor blocks, persistent providers/caches, artifacts, drivers, solvers, or
Cr2 workflow.

## Basic Gaussian Integral Arithmetic Repair

Pass 649 accepts implementation `1dbc1856d7b699311dcf65ae6a68fb193f9c1009`
for task GB-BASIC-INTEGRAL-20260918. HP-GAUSSIAN-BASIC-ARITH-FN-01/TEST-01
are implemented/completed maintenance; Pass 648's implementation grants are
consumed. The accepted repair changes only one-particle integral arithmetic,
not residual algorithms or physical-basis qualification. Independent review
passed the committed 81 checks plus 33 high-precision scratch checks; full
source CI and Docs passed. No successor task is authorized.

The accepted change affects only `polynomial_gaussian_basic_integral` in
`src/foundation/GaussianAnalyticIntegrals.jl`: +10/-10 within the 30-added-line
budget. Relative-coordinate damping/shifts replace the absolute-center
damping subtraction and cancellation-prone polynomial shifts.
No other function, caller, signature, return type, kernel, helper or file changes.
The existing kinetic caller and two-particle repairs remain unchanged.

For left/right/extra exponents `a,b,h`, centers `c,d,e`, and `g=a+b+h`, use
the algebraically equivalent relative-coordinate construction:

```text
D = (a/g)*b*(c-d)^2 + (a/g)*h*(c-e)^2 + (b/g)*h*(d-e)^2
s_left  = (b/g)*(d-c) + (h/g)*(e-c)
s_right = (a/g)*(c-d) + (h/g)*(e-d)
mu = c + s_left
```

Use the two independently formed relative shifts for the primitive factors
and `mu` for the absolute-coordinate `xpower` factor. Preserve the existing
polynomial multiplication, central moments, summation, prefactors and errors.
Preserve the complete accepted argument domain: positive total exponent is
the existing guard, not positivity of each exponent. Signed/zero exponents
with positive total remain accepted. Negative primitive powers retain their
existing rejection; nonpositive `xpower` retains its current ignored-factor
behavior. Do not add finite-input restrictions, damping clamps or fallback
branches. Do not require translation invariance for positive `xpower`:
its moments transform by the binomial rule under common translation.

### Focused Acceptance

Accepted +35 lines in existing `test/core/runtests.jl`, within the 40-line
budget; no new test owner or assertions elsewhere. Maintain table-driven use of the existing
high-precision axis oracle, with independent binomial expansion for absolute
moments. These checks catch translated overlap and polynomial-moment errors
that existing endpoint tests neither isolate nor independently reference.

- Freeze the actual primitive case
  `(200.37578609661574,-8.074717920578998,0,1.0,1776.776,-8.1,0,5.7993350529830705)`:
  reference `0.2060379396832722746153333949198068445882`, absolute error at most
  `5e-16`; the old implementation must fail this regression.
- Include centered and large-common-translation cases, diffuse/ordinary/tight
  exponents, s/p and higher polynomial powers, optional extra Gaussian,
  positive absolute moments and their binomial translation law. For bounded
  nondegenerate fixtures use `rtol=5e-13, atol=1e-24`; kinetic/covariance
  cancellation checks may use `atol=1e-20` with the same relative tolerance.
- Preserve accepted signed/zero inputs, prefactor signs/zero, nonpositive
  `xpower` and invalid-total/negative-primitive-power errors. Validate kinetic
  against independent second-derivative integration, not itself.
- Run existing core, Cartesian public/residual-GTO, IDA (including sliced
  chain), atomic-packet and public screening owners; retain their assertions
  and tolerances. Reuse overlapping full-CI results rather than duplicate
  paper examples. Require all three normal numerical CI jobs and Docs.
- Package load, docs tests, authority check/self-test, deterministic generated
  views, Documenter, manager-log bound and diff checks must pass. Record small
  warmed scalar/moment/kinetic time and allocation comparisons, not end-to-end
  speed claims. Scratch stays machine-local; no q9 basis or old audit replay.

### Evidence And Stop Boundary

Baseline `68a05b6c1` and the completed 140-second audit are preserved in
`tmp/reviews/h10-q9-residual-audit-2026-09-18.md` (SHA-256
`56cc10241f63bf31b872b3d20e24e2740a2757bceb4c54e27a05a5ea18f356c7`).
Independent design review reproduced the primitive error `-1.0703271e-11`;
relative-coordinate scratch reduced it to `2.85e-17`. Its final moment-test
summary was not recovered and is not passing evidence. New focused acceptance
is mandatory. Existing rounded-input merge success does not certify physical
residuals; real base-metric defects and input errors remain distinct.

If accepted argument behavior, budget or gates cannot be preserved, make no
implementation commit and report the obstruction. Do not loosen tests or
adapt callers. Exclude residual merge/validation changes, cutoff/rank/tolerance
changes, q9 loading/rebuild/requalification, H10 operators/fields/HF, source
monkey-patching, new dependencies/caches/APIs, workflows and releases.
No q8 energy claim follows. Hchain-doer stays paused. Accurate-input residual
requalification, including actual base metrics and independent physical checks,
requires a separate approved task; this grant ends at kernel acceptance and
lifecycle closeout, with no automatic successor.
