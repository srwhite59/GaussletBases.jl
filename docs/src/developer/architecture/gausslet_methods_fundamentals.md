# Software Packet: Gausslet Methods Fundamentals

Date: 2026-06-03
Owner: archive-manager@faraday

This packet is a shared conceptual gate for agents that may interpret or change
Gausslet semantics. It is intentionally paper-centered: the goal is not to
replace the papers or repo docs, but to prevent category mistakes when managers
and doers work across ordinary gausslets, Qiu-White hybrids, White-Lindsey
nesting, PGDG, radial/angular variants, and the newer PQS source-box line.

Use this packet before semantic/API decisions in `GaussletBases`, high-order
nesting, CR2-facing basis work, or manuscript background work.

## Required Source Material

Read the relevant subset for the task, but a new manager touching Gausslet
semantics should at least skim all of these:

- Qiu-White / hybrid grid-basis source:
  `references/papers/canonical/QiuWhite.pdf`
- Qiu-White text companion:
  `references/papers/canonical/QiuWhite.txt`
- Qiu-White TeX companion:
  `references/papers/canonical/QiuWhite_source.tex`
- White-Lindsey nested Gausslet paper:
  `references/papers/canonical/White and Lindsey - 2023 - Nested Gausslet Basis Sets.pdf`
- White-Lindsey text companion:
  `references/papers/canonical/White and Lindsey - 2023 - Nested Gausslet Basis Sets.txt`
- Nested / pure-Gaussian-gausslet software lineage:
  `reports/software_reviews/nested_gausslet_software_lineage_2026-03-15.md`
- Radial Gausslets accepted/published manuscript source:
  `references/manuscripts/radial_gausslets/main.tex`
- Radial Gausslets accepted manuscript text:
  `references/manuscripts/radial_gausslets/radial_gausslets_AIP_accepted_manuscript_JCP26-AR-01490_2026-06-01.txt`
- Angular Gausslet injection manuscript material:
  `references/manuscripts/angular_spherical_injected/Injection.tex`
- Angular Gausslet injection PDF:
  `references/manuscripts/angular_spherical_injected/Injection.pdf`
- COMX angular legacy report, useful for a compact explanation of coordinate
  multiplication localization:
  `references/reports/angular_zm_legacy/ReportZloc.md`
- Current GaussletBases software packet:
  `handoff/software_packets/gaussletbases.md`

When this packet and a live repo doc disagree, do not improvise. Report the
conflict and ask which is authoritative for the current lane.

## Conceptual Lineage

### Ordinary Gausslets and IDA

The baseline gausslet idea is to build localized orthonormal functions with an
accurate integral diagonal approximation (IDA) for two-electron interactions.

IDA is not generic point quadrature. It is a structured approximation in a
particular localized basis. A quantity named `weight`, `center`, or `moment` is
not automatically a quadrature object that can be divided into arbitrary
operators. Whether a weight has quadrature meaning depends on the basis contract.

Core distinction:

- Valid IDA semantics come from the basis construction and operator contract.
- A retained vector produced by projection, cleanup, Lowdin rotation, or shell
  selection may have diagnostics called weights or moments without being a
  positive quadrature point.

This distinction is a common failure mode. If an agent treats retained PQS or
injected-space weights as ordinary quadrature weights, it is probably making a
category error.

### Qiu-White: Distorted Parent Plus GTO Supplement

The Qiu-White line is best remembered as the distorted-parent / hybrid
grid-basis line. Its important ingredients are:

- a parent grid or gausslet-like basis that can be distorted or adapted;
- Gaussian-type orbital (GTO) supplement functions;
- hybridization between grid/gausslet and GTO components;
- the ability to repair cusp/core/completeness weaknesses of a pure grid basis.

Qiu-White is not the place to attribute the whole later COMX/nesting story. It
is the right paper to consult when the question is how GTOs are added to a
gausslet/grid basis, how overlap/projection/hybridization is handled, or how
analytic Gaussian information enters a grid method.

GTOs in this lineage are external analytic/completeness/cusp ingredients. They
do not automatically inherit gausslet IDA semantics just because they are
included in a hybrid basis.

### White-Lindsey: COMX Contraction and Nesting

White-Lindsey introduces the paper-level COMX/nesting story. For current agents,
this is the central paper for understanding why local contractions and nesting
can reduce basis size while preserving much of the gausslet structure.

COMX is a foundational primitive. At the practical level, it means localizing a
finite span by diagonalizing a projected coordinate-multiplication operator in
the proper metric. It supplies centers, localized transformed functions, and
moment/polynomial structure. It is not just an implementation detail of one
high-order experiment.

Nesting then uses local COMX/product blocks and shell-like retained content to
avoid keeping every parent function everywhere. This is the conceptual root of
later source-box and shell constructions.

Do not blur these two papers:

- Qiu-White: distorted parent plus GTO supplement/hybridization.
- White-Lindsey: COMX contraction, nesting, and PGDG.

### PGDG: Pure Gaussian Distorted Gausslets

PGDG means pure Gaussian distorted gausslet. It is a White-Lindsey-side idea,
not just a generic GTO supplement.

The key idea is to avoid doing Hamiltonian integrals with literally distorted
Gaussian primitives. Instead, the distorted-gausslet functions are redefined as
exact linear combinations of undistorted Gaussian primitives with transformed
centers/widths. The resulting pure-Gaussian analogs can use analytic Gaussian
integrals, while the basis is still cleaned/localized through orthogonalization
and COMX.

This gives PGDG a different role from Qiu-White GTO supplementation:

- PGDG changes the representation of the gausslet-like parent functions so the
  distorted basis has analytic Gaussian integral machinery.
- GTO supplementation adds external quantum-chemistry Gaussian functions for
  cusp/core/completeness repair or reference comparison.

Important cautions:

- PGDG preserves an analytic integral route; it does not by itself prove
  chemical completeness or IDA quality at a chosen spacing.
- COMX cleanup/localization is part of making the pure-Gaussian analog retain
  gausslet-like structure.
- A hybrid PGDG/GTO basis has both roles at once; do not treat the supplement
  GTOs as if they were automatically ordinary gausslet IDA points.
- Production PGDG routes should not silently fall back to numerical quadrature
  or dense parent 3D matrix materialization.

### Radial Gausslets

The radial paper is the cleanest modern reference for radial gausslet
construction, radial IDA behavior, and atomic radial/Ylm calculations. It is
especially important when judging one-center atomic references, radial
resolution, and the meaning of radial function counts.

A radial/Ylm calculation can be a strong referee for atomic one-body or
many-electron energy checks, but that does not mean a Cartesian basis with the
same nominal scale is equivalent. Cartesian completeness, core resolution, and
nested/local contraction issues are separate.

### Angular Gausslets and Injection

The angular Gausslet paper is the main reference for the protected-span
injection idea.

Injection is a same-dimension replacement and relocalization construction, not
an append operation. Let `G` be an orthonormal old localized basis and let `Y`
be the smaller orthonormal subspace to protect exactly. The complete angular
construction has three stages:

1. Form the part of the old span exactly orthogonal to `Y`. With
   `C = G' * Y`, choose orthonormal columns `B` spanning `null(C')`, and set
   `Q = G * B`. The injected span is

   ```text
   F = span(Y) + span(Q) = Y + (span(G) intersect Y_perp).
   ```

   This requires `C` to have full target rank, so `B` has exactly
   `dim(G) - dim(Y)` columns and `F` has the original dimension. Rank failure
   is a stop condition, not permission to drop a protected direction.

   This replaces the approximate `Y` content of `G`; it does not append `Y` to
   all of `G`. In a nonorthogonal coefficient representation, use the stated
   physical overlap metric in every overlap and projector.

2. Recover localized representatives by projecting every old localized seed
   into the new span:

   ```text
   Gbar = P_F * G = [Y, Q] * [C'; B'].
   ```

   The columns `[Y, Q]` are a convenient orthonormal basis for the span, but
   they are not the final localized basis. This old-seed projection is an
   essential part of angular-style injection.

3. Symmetrically Lowdin-orthogonalize the projected seeds:

   ```text
   G_inj = Gbar * (Gbar' * Gbar)^(-1/2).
   ```

The final localized vectors need not be literal members of `Y`, but their span
must recover `Y` to the stated tolerance. Do not abbreviate the construction as
only "protect `Y`" or only `F = [Y, Q]`: those statements specify the span but
omit localization recovery. Also do not casually replace the old-span
nullspace construction by `(I - P_Y)G`; the latter need not lie in the old span
unless an equivalence has been established for the active geometry and metric.

The angular setting is special because the sphere has more uniform scale and
some exact angular structure, including special constant-function/Y00 behavior.
Do not copy those geometric assumptions into direct 3D Cartesian work without
checking them.

Angular injection is useful conceptual background for later ideas such as
one-body repair, protected sectors, and projected localized representatives. It
does not by itself make injected functions compatible with ordinary IDA.

### Hooke / Harmonic-Trap Fitting

The Hooke line is not a core Gausslet paper in the same sense, but it is
important method background. Its lesson is that fitting corrections to a
carefully chosen family of harmonic-trap tests can be very powerful when the
trap orbitals are exactly represented in the basis and edge traps are avoided.

For current agents, the safe takeaway is:

- exact representation of the probe/trap family matters;
- correction fitting can be highly effective;
- do not assume a correction fitted on one family repairs arbitrary missing
  locality, completeness, or IDA structure.

### MWG

MWG belongs in the approximate-interaction/repair category, not in the basic IDA
category. The exact expansion of the acronym should be taken from the active
lane docs if a code path depends on it, but the conceptual warning is stable:

- MWG may approximate interaction effects for Gaussian-like residual/protected
  sectors.
- MWG is not a proof that arbitrary injected or residual functions can be used
  as ordinary gausslet quadrature points.
- If MWG is the only thing making a construction viable, it needs direct
  numerical validation against exact or trusted reference interactions.

This distinction mattered in the generalized-gausslet S+P line: broad or
protected Gaussian sectors improved completeness, but could not automatically be
made IDA-clean.

## Current Promising Development: PQS Source Boxes

PQS means projected Q-shell. In the current repo line, it should be understood
as a source-box-first retained-rule construction, not as a support-row oracle.

Essential objects:

- Source box: the compact product-type parent block where local operators
  should be built from low-dimensional/product pieces.
- Retained rule: the rule selecting retained modes/functions from the source
  box, such as boundary/shell content.
- Source-box pair operator: the intended fast operator object assembled from
  the source box and retained transforms.
- Shell/support-row contraction: oracle/debug/compatibility evidence only,
  unless explicitly promoted by a reviewed framework decision.

Important negative rule:

- Do not define the algorithm by reverse-engineering whatever representation the
  current shell-realized fixture happens to expose.

This was a recent manager drift point. The corrected framing is:

1. Define the source-box objects and retained rules first.
2. Build pair operators from those objects.
3. Use shell/support-row oracles only to validate or debug.
4. Export facts to CR2 carefully, without prematurely defining ray/cone grouping
   semantics in the repo.

For near-term agents, PQS is the high-order topic to know. Do not try to learn
or restate every abandoned high-order branch before acting. The essential
background is COMX, nesting, source boxes, retained rules, and the oracle-vs-
algorithm boundary.

## GTO Roles Across Lines

GTOs appear in several roles, and confusing these roles causes bad designs:

- Cusp/core repair in Qiu-White-style hybrids.
- Atomic or molecular reference basis for comparison.
- Supplement functions projected into or against a gausslet basis.
- Candidate protected/injected subspace in generalized Gausslet ideas.
- Dense Gaussian reference calculations for validation.
- PGDG primitives, where undistorted Gaussians represent gausslet-like distorted
  parent functions analytically.

These roles are not equivalent. A GTO can be exactly representable or useful for
completeness without being a valid IDA point. If a method needs exact Gaussian
integrals, density fitting, MWG, projector corrections, or explicit local blocks,
state that instead of pretending the GTO sector is ordinary gausslet IDA.

## Core Failure Modes To Catch

These are the mistakes this packet is meant to prevent:

- Treating stored weights or self-integrals as quadrature weights without a
  contract.
- Treating COMX as an optional high-order detail rather than a basic W&L
  primitive.
- Attributing W&L COMX/nesting ideas to the earlier Qiu-White hybrid paper.
- Treating PGDG as just another GTO supplement, instead of as a pure-Gaussian
  representation of distorted gausslet-like functions cleaned with COMX.
- Treating GTO supplement functions as if they automatically preserve IDA.
- Treating angular injection as literal survival of injected basis vectors.
- Calling `[Y, Q]` the final localized injected basis, or omitting projection
  of the old localized seeds into the injected span before final Lowdin.
- Treating a shell/support-row oracle as the PQS source-box algorithm.
- Treating FSB/FBU agreement on reduced transformed blocks as proof of true
  distorted-cube completeness.
- Treating tiny correctness tests as enough for public readiness.
- Allowing CR2 consumer needs to define repo API before repo-manager review.

## Validation Philosophy

Correctness is necessary but not sufficient for nontrivial numerical routes.

Every nontrivial route should state:

- category: production, reference-only but usable, or diagnostic/prototype;
- what mathematical object it claims to compute;
- what fixture validates it;
- what the validation does not prove;
- expected scaling and representative timing/allocation when practical.

For basis/semantic work, also state whether the result changes:

- a trusted baseline;
- a public repo contract;
- a private implementation detail;
- a downstream CR2-facing workflow.

## Takeover Quiz

Use this quiz when onboarding repo managers, repo doers, high-order agents, or
CR2-facing agents that will touch Gausslet semantics. Answers should cite files
or papers when possible.

1. Which paper should an agent read for distorted-parent/GTO supplement
   semantics, and which paper should it read for COMX contraction and nesting?

2. Define COMX at the practical level. Why is it a basic W&L primitive rather
   than a minor high-order implementation detail?

3. What is PGDG? Explain how it differs from adding GTO supplement functions.

4. What is IDA, and why is it wrong to treat every stored `weight` as a
   quadrature weight?

5. A retained PQS vector has a stored weight-like diagnostic. What additional
   contract would be needed before using that value for IDA division?

6. Starting from an old localized orthonormal basis `G` and a protected
   orthonormal subspace `Y`, give all three stages of angular-style injection.
   How is the old-space complement formed, how are localized representatives
   recovered, and what is finally Lowdin-orthogonalized? What must be exactly
   recoverable, and what need not survive literally?

7. List three distinct roles GTOs can play in gausslet work. Why do those roles
   not imply the same operator semantics?

8. What is MWG allowed to be in a design argument, and what should it not be
   used to prove without validation?

9. What is the PQS source-box algorithmic direction, in one paragraph? Include
   the role of shell/support-row contraction.

10. Why does FSB/FBU agreement on a reduced transformed-block target not prove
   true distorted-cube completeness?

11. For a new public or CR2-facing numerical route, what evidence is needed
    beyond a tiny correctness test?

12. Create one additional quiz question that should be added to this
    fundamentals quiz. The question must target an important concept or failure
    mode from the packet that is not already well covered by questions 1-11.
    Include the expected answer and state what misunderstanding your question is
    meant to catch.

## Evaluator Key

1. Expected: Qiu-White for distorted parent plus GTO supplement/hybridization;
   White-Lindsey for COMX contraction and nesting. Catches paper-lineage drift.

2. Expected: COMX localizes a finite span by diagonalizing a coordinate-
   multiplication operator in the relevant metric, producing localized modes
   with centers/moment content. It underlies W&L contraction/nesting. Catches
   agents treating COMX as incidental code.

3. Expected: PGDG redefines distorted gausslet-like functions as exact linear
   combinations of undistorted Gaussian primitives with transformed centers and
   widths, then uses orthogonalization/COMX cleanup to recover gausslet-like
   localized structure while keeping analytic Gaussian integrals. GTO
   supplementation is instead an added external Gaussian subspace for
   cusp/core/completeness/reference roles. Catches PGDG/GTO conflation.

4. Expected: IDA is a structured integral diagonal approximation tied to a
   specific localized basis contract. A weight may be a self-integral,
   normalization, or diagnostic; it is not automatically a positive quadrature
   rule. Catches the retained-weight/quadrature category error.

5. Expected: need an explicit basis/operator contract saying retained functions
   behave as IDA points with appropriate positivity/normalization/moment
   behavior; otherwise weights are only metadata/diagnostics. Catches PQS
   retained-weight misuse.

6. Expected: require full target rank in `C = G'Y`, form an orthonormal
   `B = null(C')`, and
   `Q = G B`; define the same-dimension injected span
   `F = span(Y) + span(Q)`; project every old localized seed into that span,
   `Gbar = P_F G = [Y,Q][C';B']`; then symmetrically Lowdin-orthogonalize
   `Gbar`. The final span must recover `Y` exactly to tolerance, but individual
   final vectors need not equal the injected functions. Catches append-only,
   span-only, and over-literal injection readings.

7. Expected examples: cusp/core repair, reference GTO basis, supplement,
   protected/injected subspace, dense Gaussian validation. These do not share
   IDA semantics because representation/completeness is not operator diagonal
   structure. Catches GTO role conflation.

8. Expected: MWG can be an approximate interaction or repair candidate that
   needs validation. It cannot be used as a proof that arbitrary residual or
   injected sectors are IDA-compatible. Catches generalized-gausslet overreach.

9. Expected: define source-box objects and retained rules first; build pair
   operators from compact source-box/product structure; use shell/support-row
   contraction as oracle/debug unless explicitly promoted by reviewed design.
   Catches the recent PQS drift.

10. Expected: reduced transformed blocks are already compressed targets; true
   distorted cubes are larger raw parent spans. Equality on the reduced target
   does not prove completeness against the raw distorted cube or physical
   Gaussian/chemical targets. Catches high-order completeness overclaim.

11. Expected: category, mathematical object, baseline, realistic fixture,
    validation meaning, performance/scaling/timing/allocation, and explicit
    readiness decision. Catches minimal-test-is-enough behavior.

12. Expected: a real gap or under-tested concept, not a rephrasing of an
    existing question. Strong candidates include paper lineage, GTO role
    conflation, MWG overclaiming, Hooke trap exact-representation logic,
    radial/Ylm referee limits, performance-as-correctness, or role-boundary/API
    drift. The expected answer should be specific enough to evaluate, and the
    misunderstanding caught should be explicit.

## Recommended Use By Role

### Repo-manager

Read this before approving public API, docs, tests, or semantic changes
involving COMX, IDA, PGDG, GTO supplements, PQS, retained weights, or
source-box operators. Use the quiz when rotating managers.

### Repo-doer

Read this before implementing tasks where names such as `weight`, `center`,
`PQS`, `source box`, `GTO`, `PGDG`, `COMX`, `IDA`, or `MWG` affect the
algorithm. If a blurb uses those terms without enough context, ask for the
governing contract.

### High-order manager/doer

Read the COMX/W&L/PQS sections. Do not broaden into every high-order historical
branch unless the task requires it. Keep FSB/FBU reduced-target facts separate
from true distorted-cube or chemical completeness.

### CR2 agents

Use this to interpret repo capabilities safely. CR2 can propose requirements and
test downstream usefulness, but should not define repo API or retained-weight
operator semantics directly.
