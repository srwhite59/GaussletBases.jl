# Performance Review Contracts

Performance is part of correctness for public routes, advertised reference
utilities, and any implementation likely to be used by downstream scientific
drivers. A feature is not review-complete until its CPU time and allocation
behavior have been checked against a reasonable scale model. For nontrivial
algorithmic work, performance and code organization also belong in the design
stage, not only in end-stage certification.

This policy exists because a dense Gaussian Coulomb reference helper was
correct on tiny tests but was unusable on a still-small one-center `25`-orbital
`625 x 625` pair matrix. The output size was modest, so minutes of runtime and
billions of allocations indicated a review failure, not merely an optimization
opportunity.

## Expected Workflow

For nontrivial coding work, the expected workflow is:

1. define the computations, operators, or data transforms that are actually
   needed.
2. do an early performance/scaling design pass:
   - identify likely cost centers.
   - decide what should be tabulated, cached, reused, or contracted early.
   - choose representations that avoid avoidable dense work, repeated setup, or
     bad scaling.
3. do an early code-organization/reuse pass:
   - reuse existing kernels and contracts when they fit.
   - avoid duplicated formulas or parallel private implementations.
   - choose an implementation seam that can be reviewed and maintained.
4. implement the route.
5. do the end-stage performance/readiness review.

This does not mean every change needs an architecture exercise. Clearly trivial
edits should remain small. The point is that nontrivial numerical routes should
not be designed as a minimal correctness patch with performance and structure
left as afterthoughts.

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

- `Computation definition`: what operators/data transforms are being built.
- `Performance strategy`: expected cost centers and chosen representation.
- `Organization strategy`: reuse points, new seams, and duplication risks.
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
