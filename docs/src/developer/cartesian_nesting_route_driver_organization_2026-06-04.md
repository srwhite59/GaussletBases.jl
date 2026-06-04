# Cartesian Nesting Route Driver Organization Memo

## Purpose

Record the architectural conclusions from the June 2026 PQS route-driver
cleanup and White-Lindsey benchmark discussion.

This note is meant to be durable guidance for future repo-manager and doer
work. It is not a pass/fail status report for the current private driver.

Terminology:

- PQS means Projected q-Shell: a retained boundary object built by projecting
  raw product-box modes, not by assigning positive retained quadrature weights.
- Source-box route means the current PQS/source-box-first algorithmic
  direction, where retained blocks are meant to come from source-box pair
  transforms.
- White-Lindsey low-order route means the old published Cartesian nested
  gausslet baseline, used here as the prior-work comparison route. The durable
  spelling should be `Lindsey`.
- Standard units means the modern top-level partition into Cartesian
  source/product units. A one-atom problem may produce one unsplit unit; a
  diatomic or larger system may produce several units. The route code should
  treat these as outcomes of the same organization, not as separate world
  views.

## Main Judgment

The top-level route driver should be treated as a design tool, not only as a
script.

Now that the private driver is small enough to read, it exposes the desired
architecture better than the lower-level helper files do. Future cleanup should
let that top-level organization push downward:

- visible route-family selection;
- system, spacing, probe, and route-recipe bundles;
- standard setup;
- unit/box selection;
- route skeleton;
- route facts;
- route-specific contract metadata;
- diagnostics, reporting, and optional artifact save.

Bottom-up pressure is still valid. If a lower-level implementation makes a
top-level distinction awkward or false, the driver should adjust. But the
default pressure should come from the organized driver down into the helper and
source-file structure. If a helper file becomes much messier than the driver,
that is evidence that the helper layer is not yet organized around the right
concepts.

## Atom Count Should Not Be A Route Boundary

The distinction between one center, two centers, and more centers should
mostly disappear from the route layer.

The same standard construction idea should apply to a single atom. It simply
does not need to split into multiple product units. A diatomic may split, and a
larger system may split more, but those are decomposition outcomes under a
common unit/box organization.

The route layer should therefore avoid logic shaped like:

- atomic route;
- diatomic route;
- chain route;
- square route.

Those distinctions can still exist where they are real geometry policies, but
they should not define different high-level operator or benchmark routes unless
the mathematics genuinely differs. The more useful top-level abstraction is:

- construct or describe standard units;
- apply a route family to those units;
- report unit-level facts and route-level pending facts.

This applies to PQS as well. A single-atom PQS construction should be a simple
case of the standard PQS route, not a separate conceptual branch.

## White-Lindsey As A Low-Order Benchmark Route

The White-Lindsey work should be preserved as a published low-order Cartesian
benchmark, not as a commitment to clone historical code details.

The important comparison for a PQS paper is:

- old low-order Cartesian nesting baseline;
- modern PQS/source-box route;
- same or comparable physical problem and unit organization where possible.

The relatively unimportant part is the old code's exact rule for when to split
the system into boxes. That rule was not the scientific core of the method, and
it should not become a route contract in the current repo.

The durable reading is:

- White-Lindsey owns a low-order COMX-style retained-basis treatment inside
  chosen shells or units.
- The modern repo owns the outward unit/box organization.
- The benchmark route should fit the modern unit/box description, rather than
  forcing the modern code to reproduce old atom-count split heuristics.

This means the route skeleton can honestly say:

- standard unit/box organization is intended;
- atom-count special casing is not required;
- historical split heuristics are not preserved as a contract;
- materialized low-order unit transforms and operator blocks remain pending
  implementation work.

## Two Useful Extremes

PQS and the White-Lindsey low-order method should both remain visible because
they exercise different extremes of the architecture.

White-Lindsey low-order nesting has disjoint units. Once the units are chosen,
they are comparatively simple: they do not require a later projection,
orthogonalization, or Lowdin-style cleanup stage to resolve overlap between
retained boundary objects. The difficulty is pushed earlier into the breakup of
space into pieces, especially the edge, face, and corner bookkeeping needed by
a low-order Cartesian shell decomposition.

PQS points in the opposite direction. The unit description can be simpler and
more regular at the route level, but the retained objects are projected
boundary objects. That means the later projection, cleanup, and contract
checks are part of the route's real complexity.

These two extremes are useful design anchors:

- White-Lindsey: complicated spatial breakup, simple disjoint retained units.
- PQS: simpler route units, more involved projected retained units and cleanup.

A code organization that can honestly represent both should be able to cover
most intermediate Cartesian nesting ideas. This is another reason not to erase
the low-order route as merely a historical artifact, and also not to preserve
its exact old split heuristic as the modern contract.

## Driver As A Template

The current private route driver should be viewed as a template for
well-organized top-level drivers in this repo.

A good top-level driver should:

- show editable defaults at the top;
- group related inputs into small named bundles;
- keep the route order visible in the file;
- call helpers for field packing, report formatting, and artifact writing;
- avoid long one-argument-per-line walls when a compact grouped call is easier
  to scan;
- keep route-family differences visible at the recipe and contract level;
- avoid burying all meaningful work in one opaque `dry_run` helper.

The driver does not need to expose every implementation detail. But a reader
should be able to answer these questions without studying a large helper file:

- What route family is being run?
- What physical system and spacing policy are being used?
- What units or route skeleton will be reported?
- What is a real algorithmic path and what is only reference/debug metadata?
- What facts remain pending before this is more than a metadata dry run?

## Consequences For Helper And Source Organization

Helper code should mirror the driver stages. If the driver has clear stages but
the helper file is a flat pile of unrelated packing logic, that mismatch is a
refactor signal.

Useful follow-through directions:

- split helper code by route stage when it grows;
- promote recurring input bundles or route descriptions into named structs only
  when they reduce cognitive load;
- consider modules for stable route-family concepts when they become larger
  than one helper section;
- use selective typing in private helpers where explicit types enforce a real
  contract, but avoid visual type noise where it only obscures the algorithm;
- keep source-box/PQS and low-order benchmark contracts separate, even when
  they share standard units.

Typed structs and modules should not be added just because they are more
formal. They are worthwhile when they make the route organization easier to
read, prevent conceptual drift, or remove repetitive named-tuple plumbing.

## Contract Boundaries To Preserve

Several existing lane boundaries remain important:

- Shell/support-row contraction is an oracle/debug path, not the PQS
  source-box algorithm.
- Retained PQS columns are not positive quadrature weights.
- IDA/source-box comparisons and MWG/final-residual supplement logic remain
  separate concepts.
- A metadata skeleton is not a Hamiltonian implementation.
- White-Lindsey low-order metadata is a benchmark contract until materialized
  unit transforms, operator blocks, and validation comparisons exist.

These boundaries should be visible in diagnostics and reports. They should not
be hidden in helper names or assumed from context.

## Practical Next Steps

Near-term work should keep the driver reviewable while pushing organization
downward:

- keep the PQS/source-box route as the default route family;
- keep the White-Lindsey route as a low-order benchmark skeleton until it has
  materialized unit transforms;
- make lower-level helper names match the driver stages;
- avoid preserving historical atom-count split behavior unless a specific
  validation need demands it;
- add benchmark comparisons only after the low-order route reports enough
  common unit facts to compare fairly with PQS.

The desired end state is one clear Cartesian nesting driver shape that can
describe multiple route families without making single atoms, diatomics, and
larger systems feel like unrelated architectures.
