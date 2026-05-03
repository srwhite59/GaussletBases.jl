# Performance Review Contracts

Performance is part of correctness for public routes, advertised reference
utilities, and any implementation likely to be used by downstream scientific
drivers. A feature is not review-complete until its CPU time and allocation
behavior have been checked against a reasonable scale model.

This policy exists because a dense Gaussian Coulomb reference helper was
correct on tiny tests but was unusable on a still-small one-center `25`-orbital
`625 x 625` pair matrix. The output size was modest, so minutes of runtime and
billions of allocations indicated a review failure, not merely an optimization
opportunity.

## Route Categories

Every new or materially changed route should identify one of these categories:

- `production route`: expected to be used directly for real calculations.
- `reference-only but usable`: dense or diagnostic reference code that should
  work at the advertised small fixture size.
- `diagnostic/prototype`: exploratory code that may be slow, but must say so
  clearly and should have guardrails.

The category controls the performance evidence required for review. Public API
entry points should not silently be `diagnostic/prototype`.

## Required Scale Contract

For a new public route, a new algorithmic path, or a major performance-sensitive
change, the handoff or review note must state:

- expected asymptotic work and memory scaling.
- representative fixture size used for performance validation.
- observed wall time or CPU time.
- observed allocations or memory footprint.
- durable artifact path for timing/allocation output when practical.
- whether observed behavior matches the scale model.

The representative fixture should be large enough to exercise the intended
use. A tiny correctness fixture is not a performance validation.

## Back-Of-The-Envelope Review

Reviewers should compare the measured result to a simple model. If the model
says the output is small but the run allocates heavily or takes minutes, the
route is not ready, even if numerical tests pass.

Examples:

- A dense `N^4` Gaussian pair matrix with `N = 25` has a `625 x 625` output.
  That is a small dense matrix. Very large allocation counts indicate
  avoidable inner-loop allocation or an accidental algorithmic factor.
- A PGDG or mapped ordinary path that is intended to be analytic should not
  silently call numerical quadrature.
- A branch/counterpoise correction helper should validate the branch sizes and
  centers used by realistic molecular fixtures, not only a one-center atom.

## Handoff Checklist

Doer handoffs for relevant changes should include:

- `Performance category`: production route, reference-only but usable, or
  diagnostic/prototype.
- `Scale model`: short expected CPU and memory scaling.
- `Representative fixture`: problem size and why it is representative.
- `Measured performance`: time and allocations or memory.
- `Artifact`: path to timing/allocation output if available.
- `Decision`: ready, needs optimization, or deliberately prototype-only.

If performance is not measured, the report must say why and should not call the
feature fully validated.

## Guardrails

Guardrails are part of the public contract:

- Dense reference utilities should expose size limits or clear warnings.
- Slow diagnostic paths should have explicit names or docstrings marking them
  as diagnostic/prototype.
- Production routes should fail loudly rather than silently taking a much slower
  fallback path.
- Performance artifacts should be durable enough for the next agent or reviewer
  to inspect, not only terminal scrollback.

## Documentation Requirement

Public APIs with nontrivial scaling should document the scale and intended use
in the user-facing reference page. Developer-only caveats are not enough for a
discoverable public utility.
