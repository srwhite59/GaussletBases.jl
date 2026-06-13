# PQS Source-Box Fixture Policy

This note classifies PQS source-box fixtures before RHF/SCF work. It is a
private development policy, not a public API contract.

## Fixture Roles

- Route smoke: proves a route-owned object path materializes and preserves
  non-promotion flags. It should assert compact status, dimensions, blockers,
  and materialization flags, not physical accuracy.
- Convention diagnostic: checks a numerical convention such as the pre-final
  positive-weight density gauge or a self-Coulomb scalar. It does not promote
  that scalar to endpoint validation.
- Oracle/debug: compares against shell/support-row, fixed-block, or
  explicit-box paths. These paths help detect regressions but are not PQS route
  authority.
- Physics endpoint: validates a named physical target with explicit system,
  resolution, reference value, tolerance, and timing policy.

## Current Compact Fixtures

- Tracked Z=1 H1 seam:
  - file: `test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - parent box: `7 x 7 x 7`
  - source/core box: `5 x 5 x 5`
  - final dimension: `223`
  - H1 energy: `-0.48047934800387226`
  - role: route seam plus fixed-block H1 oracle/debug check.
- Direct structured H1/J convention probe:
  - same compact source/core family and final dimension `223`
  - H1 energy: `-0.48047934800387126`
  - self-Coulomb: `0.6397851751855723`
  - density gauge: `pre_final_localized_positive_weight`
  - role: convention diagnostic, not endpoint validation.
- One-center source-box driver H1/J dry-run:
  - metadata uses one center with `Z = 4`
  - `q = 5`, `n_s = 5`
  - parent box: `7 x 7 x 7`
  - source/core family: `5 x 5 x 5`
  - final dimension: `223`
  - H1 energy: `-5.6629907690725245`
  - self-Coulomb: `1.8691288063594704`
  - density gauge: `pre_final_localized_positive_weight`
  - role: private driver route smoke, not a neutral-Be physics endpoint.

## Nonclaims

- Compact H1/J route materialization is not RHF readiness.
- `final_dimension == 223` is not physics acceptance.
- A self-Coulomb scalar alone is not endpoint validation.
- Shell/support-row, fixed-block, and explicit-box paths are oracle/debug, not
  route authority.
- Retained diagnostic weights are not final IDA or quadrature weights.

## Coupled Fixture Parameters

Physics comparisons must move these as reviewed fixture families, not as
independent one-off knobs:

- `Z`, electron count, charge, and spin/closed-shell convention;
- spacing/core spacing, reference spacing, tail spacing, distortion, and
  parent radius;
- `q`, `n_s`, core/source side, shell depth, and retained rule;
- parent axis counts and source-box placement;
- Coulomb expansion and factor inputs;
- density gauge;
- fixture family, such as compact `7^3 / 5^3 / dim 223` versus side13 or a
  q-ladder row.

## Before RHF

Before implementing RHF/SCF/Fock:

- Decide whether the first RHF fixture is route-smoke-only or a physics
  endpoint.
- If it is a physics endpoint, choose the target system, `Z`, electron count,
  closed-shell rule, fixture size, reference value, error threshold, and timing
  threshold before code.
- If it is route smoke, assert only compact status and non-promotion facts.
- Keep H1/J diagnostic materialization separate from RHF readiness and route
  adoption.
