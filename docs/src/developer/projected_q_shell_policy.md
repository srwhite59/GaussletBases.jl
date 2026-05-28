# Projected q-Shell Policy

## Purpose

Projected q-Shell (PQS) is the intended high-order shell abstraction for
regular atom-local and molecular/exterior shell construction. It replaces the
earlier preferred mental model of stitching separate endcap, panel, annulus,
or leftover slab pieces in the regular middle of the construction.

This page records the policy direction only. It does not make PQS the default
builder route, does not change any Hamiltonian or QW construction path, and
does not claim production science validation.

## Definition

PQS is a raw-boundary projection of the full local block transform.

For a `q x q x L` local block:

1. Build or use the full local product/block transform for the whole block.
2. Select the boundary COMX-product modes: product columns whose local mode
   index is first or last on at least one axis.
3. Project those selected product modes onto the raw outer-shell coordinates.
4. Apply full-rank symmetric Lowdin cleanup to the projected boundary span.

The shell must not depend on how the interior block was contracted. It must
not be defined by subtracting a previously locked contracted inner cube, and
it must not depend on locked prior spans. Interior contraction is a separate
choice; PQS is defined by boundary COMX-product modes from the full local
block transform projected onto raw boundary rows.

The subtraction formula is a counting picture, not the construction contract:

```text
dim PQS(q, L) = q^2 L - (q - 2)^2 (L - 2)
```

For `q = 5`, this gives `16L + 18`. The cubic case is `L = q`, so
`PQS(5, 5)` has dimension `5^3 - 3^3 = 98`.

## Policy

- Cubic atom shells are `PQS(q, q)`.
- Rectangular molecular/exterior shells are `PQS(q, L)`.
- Regular shells should be represented as PQS boundary traces, not as
  separate endcap/panel stitching.
- PQS must be independent of the contraction chosen for the interior block.
- Fit-to-box adjustments happen only at the outermost parent boundary.
- Generous boxes may omit awkward outer padding, but diagnostics must say
  exactly what was omitted.
- Regular mid-region construction should not leave leftover slabs.
- Regular construction should not create far split-box shells.
- PGDG analytic integrals remain mandatory unless an explicit
  `:numerical_reference` or debug path is requested.

The existing q-row and endcap/panel work in mainline is transitional
infrastructure and validation scaffolding. It remains useful for diagnostics,
receipt wrapping, sidecar plumbing, and CR2 handoffs, but it is not the final
preferred shell abstraction.

## External First-Gate Evidence

High-order ran a q=5 unified block-shell C2 `R = 4.7` preflight. The unified
shell row matched the q5 annulus dimension, beat the q4 and q5 annulus rows on
RHF/EGOI for that point, had clean support accounting, avoided dense fallback,
used the full analytic PGDG-style path, had exact H/V symmetry, and had clean
by-center/total H agreement.

The reported rows were:

| route | dimension | RHF total | EGOI relative Frobenius |
|---|---:|---:|---:|
| q4 four-panel/direct5 | `1482` | `-74.881526646276` | `0.00135160` |
| q5 true annulus | `1554` | `-74.888408985519` | `0.00142188` |
| q5 unified block shell/PQS candidate | `1554` | `-74.891723662832` | `0.000866239` |

This is positive first-gate evidence for PQS as the intended abstraction. It
deserves an accepted-R ladder follow-up before any default-policy decision.
It is not enough by itself to make PQS the default route.

## Mainline Implementation Checkpoint

As of 2026-05-28, mainline has a private opt-in PQS construction smoke. The
implementation trail is:

- `06483f8`: added the private projected q-shell local layer helper
- `bd66e51`: exposed PQS as a realization candidate in metadata
- `e4fbcca`: added opt-in PQS source and fixed-block construction
- `42103e3`: added the pure nested PGDG QW smoke for the opt-in PQS fixed
  block

The helper follows the intended boundary-mode contract: build the full local
COMX transforms, select boundary COMX-product modes, project those modes onto
raw boundary rows, and apply full-rank symmetric Lowdin cleanup. It does not
subtract a contracted inner cube and does not project against locked prior
spans.

The path remains private and explicit opt-in. Existing defaults continue to
use the endcap/panel route where they did before. The first small fixture has
full-parent dimension `735`; that is construction and QW-smoke evidence only,
not a compression-quality, energy, CR2, or production science claim.

The pure nested QW smoke uses the existing receipt path with
`gausslet_backend = :pgdg_localized_experimental`,
`interaction_treatment = :ggt_nearest`, and `nuclear_term_storage =
:total_only`. It reports clean source/sidecar agreement, finite symmetric
overlap/one-body/interaction matrices, zero residuals, and no warning-level
numerical-quadrature logs for that path.

The missing next design problem is product-staged PQS sidecar and performance
work. PQS should not be adopted for by-center, supplement, or
performance-sensitive routes until that sidecar/performance contract exists.

## Consequences For Mainline

The next source-construction design should treat PQS as the regular shell
target:

- atom-local q shells should move toward `PQS(q, q)`;
- molecular/exterior q shells should move toward `PQS(q, L)`;
- existing endcap/panel helpers should be treated as transitional support and
  audit scaffolding unless explicitly retained for a boundary-only special
  case;
- diagnostics should distinguish construction by raw-boundary projection from
  any counting-only annulus formula;
- any route that trims or omits outer support for fit-to-box reasons must
  report that as an outer-boundary adjustment, not as a regular shell policy.

Before adoption, repo-side work still needs a dedicated implementation and
validation plan. That plan should compare PQS against the transitional
endcap/panel and annulus rows without changing public defaults, backend
policy, Hamiltonian kernels, or quadrature behavior.
