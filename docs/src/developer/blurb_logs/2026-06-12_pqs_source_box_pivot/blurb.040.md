Purpose:
  Run a physically interpretable one-center PQS H1 diagnostic through the
  direct-retained final-basis path before starting IDA or RHF work.

Context:
  The current durable H1 gate proves the workflow seam:

  raw product source box -> boundary retained source modes -> direct retained
  overlap/kinetic/nuclear -> final basis -> one-electron Hamiltonian -> ordinary
  eigensolve.

  Its eigenvalue is a workflow/oracle diagnostic, not yet a meaningful
  physical H/He+ acceptance value. Before adding final IDA weights or
  density-density, we need a physical one-center H1 probe using the same direct
  retained path.

Exact task:
  1. Create or update an ignored `tmp/work` probe for a one-center physical PQS
     H1 diagnostic. Use the current direct-retained final-basis path and no
     `_pqs_current_route_safe_term_matrices(...)`.

  2. Choose a small reviewed fixture first. Prefer the current `q=5/L=5`
     source-box-style setup if it can be mapped to a clear one-center physical
     box. Use either:

     - H / Z=1, compare to exact H1 `-0.5`, or
     - He+ / Z=2, compare to exact H1 `-2.0`.

     If neither is straightforward from the current PQS source-box machinery,
     report the exact missing mapping/shell-realization input instead of
     inventing a broad route.

  3. Preserve the direct retained one-body checks:

     - overlap/kinetic/by-center nuclear direct-boundary flags are true;
     - raw source one-body blocks are not materialized on the active path;
     - final overlap identity error is small;
     - final Hamiltonian is finite/symmetric;
     - ordinary symmetric eigensolve, not generalized overlap solve.

  4. Compare against a shell-support oracle if already available from the same
     fixture. The oracle is for route validation; the physical exact value is
     for basis-quality interpretation.

  5. Write a concise ignored artifact, for example:

     `tmp/work/pqs_one_center_physical_h1_probe_summary.txt`

     Include fixture parameters, physical extent if known, basis/final
     dimension, H1 eigenvalue, exact-reference error, oracle error if available,
     direct-path flags, and coarse timings.

Do not:
  - add a permanent test in this pass;
  - add IDA, density-density, RHF, driver wiring, exports/artifacts, or GTO
    changes;
  - build a broad new route framework;
  - change source code unless a tiny fix is required to run the probe and is
    clearly justified.

Validation:
  - run the probe;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - probe artifact path;
  - fixture parameters and whether physical extent is known;
  - H1 eigenvalue and exact-reference error;
  - shell-support oracle comparison if available;
  - direct-retained path flags;
  - coarse timings;
  - validation run and result;
  - deletion/shrinkage report:
      - whether any old coverage/code was removed;
      - if no deletion, why this was probe-only;
      - exact next blocker or target for final IDA weights.

Continue the baton loop after writing `response.040.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini
