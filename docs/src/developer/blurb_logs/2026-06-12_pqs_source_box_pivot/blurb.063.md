Purpose:
  Run a probe-only side-13 PQS fixture ladder that varies direct core size and
  shell depth together on the same physical parent. This should answer whether
  the current side-13, core `(4:10)^3`, three-shell PQS RHF point is special, and
  should keep the WL-style fixture controls visible: box radius, central spacing,
  mapping/distortion, direct core size, and shell depth.

Context:
  Pass 062 gave a coherent side-13 multi-layer PQS He RHF result:

  - parent count 13;
  - `AsinhMapping(c = 0.1, s = 1.0, tail_spacing = 10.0)`;
  - endpoints about `+-8.565228460168399`;
  - direct core `(4:10)^3`, outer box `(1:13)^3`, three shell layers;
  - final dimension 1549;
  - Z = 2 H1 about `-1.9755618232013417`;
  - self-Coulomb J about `1.2169264388860319`;
  - RHF total about `-2.8372556463894707`;
  - about `+0.02442434922276826` Hartree above the He HF reference;
  - about `-0.0007576466885570454` Hartree relative to the WL side-13 RHF probe.

  This is a useful reference probe, not an accepted fixture. The relation among
  `Z`, core spacing `d`, mapping parameter `s`, parent radius, direct core size,
  and shell depth may need a more detailed study later, but do not try to settle
  that rule in this pass.

Exact task:
  Add and run an ignored `tmp/work` probe for the side-13 parent only:

  - parent count 13;
  - same axis mapping as pass 062: `AsinhMapping(c = 0.1, s = 1.0, tail_spacing = 10.0)`;
  - outer box `(1:13)^3`;
  - test direct core boxes:
      - `(5:9)^3`   # core side 5, four surrounding one-cell shell layers
      - `(4:10)^3`  # core side 7, three surrounding one-cell shell layers, repeat/reference
      - `(3:11)^3`  # core side 9, two surrounding one-cell shell layers
      - optionally `(2:12)^3` only if the first three are fast and stable

  For each point, report:

  - core side length and shell layer count;
  - final dimension;
  - final overlap identity error;
  - Z = 2 H1 and H1 error vs `-2`;
  - self-Coulomb J and error vs `1.25`;
  - He RHF convergence, iterations, one-electron energy, electron-electron energy,
    and total energy;
  - error vs He HF reference `-2.8616799956122388788`;
  - timing split, at least plan/final-basis, one-body/H1, density interaction,
    RHF, and total;
  - route flags: no full-parent CPB fallback, no old fixed-block pair authority,
    no ordinary Cartesian IDA fallback, no GTO/PQS driver/export/artifact.

Decision rule:
  This is a probe-only ladder. Do not add permanent tests and do not promote a
  fixture as an acceptance gate. If a point becomes very slow or fails in a way
  that suggests a route bug, stop after recording the completed points and the
  blocker. If results are nonmonotonic, report the likely fixture-control
  interpretation instead of continuing blindly.

Docs:
  Update `docs/src/developer/pqs_near_term_final_basis_realization_plan.md` with
  a concise result table and interpretation only. Keep it explicit that no new
  `Z,d,s,ns` fixture rule was adopted.

Validation:
  - Run the probe script.
  - Run `julia --project=. -e 'using GaussletBases; println("load ok")'`.
  - Run `git diff --check`.
  - No Julia test suite is required unless production source is changed, which
    should not happen in this pass.

Deletion/shrinkage report required:
  - what old code, test, metadata, or compatibility path became unnecessary;
  - what was deleted or simplified;
  - if nothing was deleted, why no existing surface was made obsolete;
  - whether any new test replaces/shrinks older coverage or is genuinely new
    live-contract coverage;
  - any remaining stale or duplicate surfaces to retire next.

-- repo-manager@macmini
