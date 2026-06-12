Purpose:
  Audit the exact implementation seam for multi-layer PQS shell/source
  construction before changing production code.

Context:
  Pass 058 showed that a WL-aligned side13-style PQS He fixture is blocked:

  ```text
  parent count: 13
  mapping:      AsinhMapping(c = 0.1, s = 1.0, tail = 10)
  current box:  (1:13)^3
  inner box:    (4:10)^3
  raw dims:     (7,7,7)
  shell layers: 3
  blocker:      projected q-shell requires a one-cell raw boundary
  ```

  The current complete core/shell path supports one direct core plus one
  surrounding shell. WL-aligned PQS needs more shell depth.

Exact task:
  Do a focused code audit and implementation plan for multi-layer PQS shells.
  Prefer no production source changes in this pass. A small `tmp/work` prototype
  is allowed if it clarifies feasibility.

Questions to answer:

  1. Where exactly is the one-cell-boundary assumption encoded?
     Include file/function names and the data shape it enforces.

  2. Can multi-layer PQS be composed from repeated existing one-cell projected
     shell layers, for example:

     ```text
     core (4:10)^3
     shell 1: (3:11)^3 around (4:10)^3
     shell 2: (2:12)^3 around (3:11)^3
     shell 3: (1:13)^3 around (2:12)^3
     ```

     or does the source-box route require a new multi-layer source producer?

  3. Does `pqs_complete_core_shell_final_basis(...)` also need to generalize
     from one shell sector to multiple shell sectors, or can the shells be
     collapsed into one shell support/final-coefficient block before calling it?

  4. What is the smallest implementation path that would let the WL-aligned
     side13 fixture reach final-basis/H1 without fixed-block authority?

  5. What old/oracle surfaces are useful references, and which must not become
     route authority?

Optional cheap prototype:
  If cheap, try a `tmp/work` prototype that builds the three one-cell shell
  descriptors separately and reports support overlap/duplication, retained
  counts, and whether their shell projections can be block-combined. Do not
  run H1 or RHF in this pass.

Artifacts:
  Write ignored artifacts:

  ```text
  tmp/work/pqs_multilayer_shell_source_seam_audit.jl
  tmp/work/pqs_multilayer_shell_source_seam_audit_summary.txt
  ```

Docs:
  Update the PQS near-term plan only if the audit yields a concise blocker or
  implementation sequence worth recording.

Do not:
  - implement a broad multi-layer framework in this pass;
  - promote q=9/q=11 gates;
  - add permanent tests;
  - run H1/RHF;
  - add GTO, driver wiring, exports, or artifacts;
  - make fixed-block matrices route authority.

Validation:
  - run the audit/prototype if created;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - files edited;
  - exact one-cell-boundary code sites;
  - whether repeated one-cell shell composition is viable;
  - whether final-basis realization must generalize beyond one shell sector;
  - recommended smallest implementation path;
  - any prototype result and timing;
  - validation run and result;
  - deletion/shrinkage report:
      - what old/fallback/oracle surface became less necessary;
      - what was deleted or simplified, if anything;
      - if nothing was deleted, why not;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.059.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini
