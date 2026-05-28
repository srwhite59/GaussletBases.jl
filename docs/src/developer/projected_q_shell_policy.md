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
2. Project that transform onto the raw outer-shell coordinates.
3. Equivalently for support accounting, zero or discard the raw inner
   `(q - 2) x (q - 2) x (L - 2)` coordinates.
4. Keep and clean/orthonormalize the resulting boundary span.

The shell must not depend on how the interior block was contracted. It must
not be defined by subtracting a previously locked contracted inner cube, and
it must not depend on locked prior spans. Interior contraction is a separate
choice; PQS is defined by the raw boundary trace of the full local block
transform.

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
