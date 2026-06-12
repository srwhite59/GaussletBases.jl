Purpose:
  Audit how the current PQS source-box/multi-layer/final-basis work maps onto
  the canonical Cartesian driver spine. This is an audit/planning pass only,
  intended to prevent the successful PQS probes from becoming a parallel private
  route island.

Context:
  Recent passes established:

  - `pqs_multilayer_shell_source_plan(...)` for repeated one-cell shell source
    layers on a fixed parent;
  - `CartesianFinalBasisRealization` for complete core/shell final-basis
    construction, final one-body transfer, H1 solve, final IDA weights, and the
    pre-final positive-weight density interaction seam;
  - direct retained PQS source one-body blocks for overlap, kinetic, position,
    x2, and by-center nuclear;
  - generic retained selector coverage for overlap, kinetic, position, and x2;
  - side13 PQS H1/J/RHF probes that are coherent but not accepted fixtures.

  The driver spine remains the route contract:

  ```text
  cartesian_system
  cartesian_recipe
  cartesian_parent
  cartesian_shells
  cartesian_units
  cartesian_transforms
  cartesian_pairs
  cartesian_assembly
  cartesian_report
  cartesian_materialization
  ```

Task:
  Read the current driver/helper code and write a concise audit. Do not change
  production source or tests unless you find a trivial stale doc typo that is
  directly in scope.

  The audit should answer:

  1. Which current PQS objects/helpers belong to each driver stage?
     For example: raw source facts, multi-layer shell source plan, retained
     source one-body blocks, final-basis realization, final one-body H,
     density interaction, and RHF probe logic.

  2. What is the first missing integration seam that prevents the current
     side13 PQS H1/J/RHF probe from being expressed through the normal driver
     lifecycle?

  3. Which existing private probes/oracles should remain private, and which
     module-owned objects are ready to be consumed by a future driver stage?

  4. What old or duplicate surfaces would become less necessary if that first
     integration seam were implemented?

  5. What should not be integrated yet because the physical fixture rule
     (`Z`, `d`, `s`, radius, core size, shell depth) is still unsettled?

Deliverable:
  Add a concise note to
  `docs/src/developer/pqs_near_term_final_basis_realization_plan.md` or a small
  adjacent developer note if the current file would become too crowded.

Do not:
  - implement driver wiring;
  - add tests;
  - run new physics probes;
  - add acceptance gates;
  - change exports;
  - study or codify a new `Z,d,s,ns` fixture rule;
  - promote the side13 RHF result.

Validation:
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - `git diff --check`

Deletion/shrinkage report required:
  - what old code, test, metadata, or compatibility path became unnecessary;
  - what was deleted or simplified;
  - if nothing was deleted, why no existing surface was made obsolete;
  - whether any new test replaces/shrinks older coverage or is genuinely new
    live-contract coverage;
  - any remaining stale or duplicate surfaces to retire next.

-- repo-manager@macmini
