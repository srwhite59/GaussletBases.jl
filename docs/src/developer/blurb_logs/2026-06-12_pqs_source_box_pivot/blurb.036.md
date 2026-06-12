Purpose:
  Prevent PQS retained overlap/kinetic construction from adopting the full
  raw-source dense-block plus selector path as the scalable algorithm.

Context:
  The current retained source-mode helpers for overlap and kinetic first build
  a full raw source-space block, then select retained boundary rows/columns:

  `O_retained = O_source[left_columns, right_columns]`

  This is acceptable as an oracle for the current `5 x 5 x 5` fixture, but it
  is the wrong scalable habit. For a boundary source-mode retained rule, the
  retained block can be built directly from retained source-mode tuples and 1D
  axis factors.

Exact task:
  1. Add a direct retained-boundary product kernel for PQS overlap and kinetic
     in the CPBM PQS source safe-term area.

  2. The direct overlap kernel should iterate retained mode tuples from
     `CRPS.retained_mode_indices(rule)` and fill:

     `S[a,b] = Sx[ix_a, ix_b] * Sy[iy_a, iy_b] * Sz[iz_a, iz_b]`

  3. The direct kinetic kernel should fill:

     `K[a,b] =
         Kx[ix_a, ix_b] * Sy[iy_a, iy_b] * Sz[iz_a, iz_b]
       + Sx[ix_a, ix_b] * Ky[iy_a, iy_b] * Sz[iz_a, iz_b]
       + Sx[ix_a, ix_b] * Sy[iy_a, iy_b] * Kz[iz_a, iz_b]`

  4. Route the active retained overlap/kinetic helpers through this direct
     retained-boundary construction when the retained rules are
     `:source_mode_column_selector` boundary product rules.

  5. Keep the existing raw-source-block -> retained selector path as oracle and
     fallback. Do not delete `pqs_source_pair_retained_one_body_block(...)`.

  6. Update metadata/comments so the active overlap/kinetic retained helpers no
     longer claim they necessarily materialize a raw source-space input block.
     Preserve clear nonclaims: no shell realization, Lowdin, final basis, IDA,
     Hamiltonian, RHF, driver, export, or artifact.

Test policy:
  Add or update one compact CPBM module-contract check only. It should compare:

  - direct retained overlap vs existing raw-source overlap -> selector oracle;
  - direct retained kinetic vs existing raw-source kinetic -> selector oracle;
  - retained shape/count for the `5 x 5 x 5` boundary rule.

  Avoid broad metadata assertions. Check only the key contract fields needed to
  prove the algorithm shape and nonclaim boundary.

Do not:
  - add direct retained nuclear yet;
  - add shell realization, final basis, H1, IDA, density-density, RHF, driver
    wiring, exports/artifacts, or GTO changes;
  - add a new broad test file;
  - remove the raw-source selector oracle path.

Validation:
  - `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - `git diff --check`

Required response:
  - files edited;
  - names of direct retained-boundary helper(s) added or changed;
  - whether active retained overlap/kinetic now avoid full raw source blocks;
  - oracle comparison result and max differences;
  - validation run and result;
  - deletion/shrinkage report:
      - whether any raw-source selector path became oracle/fallback only;
      - what comments/metadata were simplified;
      - whether tests replaced/shrank old assertions or added only compact live
        kernel coverage;
      - remaining next optimization target.

Continue the baton loop after writing `response.036.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini
