# Gausslet Algorithm Refresher

This is a short memory-refresh note for agents working on Cartesian, PQS,
White-Lindsey, Qiu-White, Residual Gaussian, or IDA/MWG code. It is not design
authority and does not approve implementation surfaces. If this note conflicts
with the compact Cartesian Hamiltonian producer authority, the compact
authority wins.

Use this page when you need the concrete algorithmic distinctions without
rereading the full design registry or the paper-centered fundamentals packet.

## One-Screen Map

Qiu-White, White-Lindsey, and PQS are related, but they are not the same
algorithm with different route labels.

- Qiu-White is the distorted-parent plus Gaussian supplement line. It combines
  a grid/gausslet-like parent with GTO supplement functions and exact Gaussian
  mixed/self blocks. In this repo it is often a donor for analytic Gaussian
  raw-block organization, not a route object to copy wholesale.
- White-Lindsey is the nested/COMX contraction line. The important step is
  local contraction: source geometry is turned into compact retained functions
  by products of one-dimensional contractions. Boundary-stratum CPBs are not
  already retained functions.
- PQS is the projected Q-shell/source-box line. It uses filled source boxes,
  boundary product-mode selection, retained-space projection/orthogonalization
  against already accepted local content where that retained rule requires it,
  and localized Lowdin cleanup to make compact terminal blocks.

Route names do not prove numerical semantics. A value named `:white_lindsey` or
`:pqs` only says which construction family was requested. The live code must
still carry the numerical objects that make that family true.

## White-Lindsey

White-Lindsey should be remembered as:

```text
nested local geometry
-> boundary faces / edges / corners / strata
-> local product-of-1D contraction
-> compact retained functions
-> terminal basis blocks
```

A boundary-stratum CPB is a geometry/source block. It says where source
content lives and how the boundary unit is shaped. It is not itself the final
retained basis.

The retained functions come from a coefficient map, conceptually:

```text
WL boundary-stratum CPB
-> local product-of-1D coefficient map
-> terminal block coefficients
```

The deleted route-global WL stack contained versions of this coefficient
primitive, but that stack also contained too much scaffolding: reports, status
vocabulary, route-global adapters, old materialization paths, and tests for
transitional surfaces. Do not revive that stack. If a source pass needs the
old numerical primitive, mine or re-express only the compact local coefficient
construction behind the current terminal-basis boundary.

Identity rows are valid for true direct/core identity units: a basis function
is exactly a parent row on authoritative owned support. Identity rows are not
valid for a metadata-only `:white_lindsey_boundary_stratum_retained_unit`.
Treating boundary-stratum support rows as retained functions is the recent WL
diatomic mistake.

## PQS Source Boxes

PQS should be remembered as:

```text
filled source CPB
-> boundary COMX-product mode selection
-> retained-space projection / orthogonalization
-> shell-owned row restriction
-> shell-local symmetric Lowdin
-> retained localized terminal block
```

The source box is filled because product-mode selection needs the local tensor
structure. The terminal block is localized because the final coefficients live
only on the authoritative owned shell support.

There is one important current-repo guardrail: do not implement previous
terminal-block projection that grows a shell into earlier terminal regions.
Older or more abstract PQS descriptions may talk about projection against
already accepted content. In the live Cartesian terminal-basis authority, final
blocks have disjoint owned supports and must not repair themselves by mixing
into earlier supports. The correct terminal realization is:

```text
full source-box product modes
-> boundary product-mode columns
-> restrict rows to shell-owned support
-> shell-local Gram
-> shell-local symmetric Lowdin
-> sign canonicalization
```

Cross-block overlap is zero by support structure. Cross-block kinetic,
nuclear-attraction, and IDA interactions can still be nonzero and are assembled
as operator matrix elements.

## Direct/Core/Slab Identity Sectors

Identity support rows are valid when the sector is intentionally direct:

- a complete/core block that keeps parent rows as final basis rows;
- a slab/direct unit whose final function really is the corresponding parent
  support row;
- a true identity terminal block on disjoint owned support.

Identity support rows are not valid merely because a record has support
indices. Support rows are where a function may live; they are not the function.
A retained unit, shell, stratum, route label, or provenance row is numerical
authority only when it carries or can construct the coefficients for the
retained functions.

Practical rule: if a unit name says "boundary", "stratum", "source", or
"retained", check whether coefficients have been constructed. Do not assume
the CPB support is the retained basis.

## COMX And PGDG

COMX is a coordinate-multiplication/localization primitive. Operationally, it
takes a local finite span and constructs localized modes using coordinate
operator structure. It is central to White-Lindsey-style contraction and to how
local product blocks become useful retained functions. It is not just a
diagnostic or a high-order decoration.

PGDG means pure Gaussian distorted gausslets. The point is to represent
distorted gausslet-like functions using Gaussian primitives so analytic
integral machinery can be used while retaining gausslet-like locality and COMX
cleanup. PGDG is therefore different from generic GTO supplementation:

- PGDG changes how the gausslet-like parent functions are represented.
- GTO supplementation adds extra candidate functions to repair or augment the
  basis.
- A hybrid route may use both ideas, but they answer different questions.

If a path uses Gaussian primitives, do not immediately call it a supplement.
Ask whether those Gaussians represent the parent gausslet-like basis itself or
extra functions added around that basis.

## IDA, MWG, And GTO Supplements

IDA is a two-electron approximation tied to a localized gausslet-style basis.
It is not generic point quadrature. A stored center, moment, or weight-like
diagnostic does not automatically authorize density normalization or Coulomb
assembly.

Residual GTO/MWG supplementation follows a different split:

- exact one-body operators stay exact through transformed overlap, kinetic,
  nuclear attraction, coordinate, and second-moment matrices;
- residual functions are selected by owner-local residual content and then
  merged into an orthonormal residual basis;
- matched-width Gaussians are an interaction approximation built from moments
  of the final residual functions;
- residual rotations matter because MWG descriptors are not invariant under
  arbitrary rotations.

Base gausslet IDA semantics do not automatically transfer to residual GTO
functions. Residual integral weights may be small, negative, or near zero.
MWG uses moment-matched positive Gaussian surrogates for residual-containing
interaction blocks and needs its own validation.

## Core Scale And Spacing Names

Current public convention:

- `core_spacing` is the physical near-core spacing after explicit input,
  driver default, or preset resolution;
- legacy White-Lindsey one-center `d` is the old name for the same physical
  scale in that mapping context;
- `reference_spacing` is the reference-grid spacing and stays separate;
- `tail_spacing` is a separate outer/tail spacing concept;
- public `ns` is the requested cube/source/nesting size; route-local `q` is
  derived from `ns` and `nesting`.

Do not reintroduce a second public `d` knob. If compatibility code accepts it
temporarily, it must agree with the resolved `core_spacing`.

## Common Mistakes

- Support rows are not retained basis functions. They are rows where a function
  may have coefficients.
- Metadata and provenance are not numerical authority. They help consumers
  interpret artifacts, but they do not replace coefficients, transforms, or
  operator matrices.
- Route labels do not imply numerical semantics. A route named WL is not a WL
  retained basis unless the WL contraction has happened.
- PQS diagnostic weights are not automatically quadrature weights. Check the
  basis and operator contract before dividing or normalizing by them.
- Shell/support-row oracle paths are not production algorithms. They may prove
  an endpoint or compare dimensions, but they do not define the live route.
- Old route-global scaffolding should not be revived. Restore only the small
  numerical primitive that is still needed, behind the current compact owner.
- Identity blocks are valid only for direct identity sectors. Boundary
  source/stratum units need their contraction or Lowdin construction first.

## When In Doubt

Ask two questions before editing:

1. What is the numerical object being constructed: support geometry,
   coefficients, transformed operators, IDA/MWG interaction, or provenance?
2. Which current owner is responsible for that object?

If the answer is "a route label" or "a metadata record", stop. Find the
coefficient map, transform, raw block, or operator contract that makes the
algorithm real.
