Pass 183 - recover/compute the matching WL gausslet-only q=5/n_s=5 He baseline for comparison hygiene.

Current accepted PQS diagnostic result:

```text
PQS fixed-q one-center He
final dimension = 419
H1 energy = -1.9866819751748936
H1-J self-Coulomb = 1.2261626003119184
private RHF diagnostic total = -2.8493871224303167
```

The next question is comparison hygiene. Do not compare this 419-dimensional
gausslet-only PQS value directly to the old supplemented WL Fig.8-style final
447 value unless you state the mismatch clearly.

Goal:

```text
Find or compute the matching White--Lindsey gausslet-only q=5/n_s=5 He result:
  old Fig.8-style one-center WL construction
  d = 0.3
  fixed local 5^3 core
  three fixed-q shell sectors
  gausslet-only retained dimension = 419
  no AHGBS-9 residual/supplement functions
```

This is a no-tracked-source/test pass unless you hit a genuine blocker. Use
ignored `tmp/work` scripts/reports if computation is needed.

Do not use these as the comparison baseline:

- current compact decomposed WL 223 readiness route;
- current atom-growth/rectangular WL 2287 route;
- old supplemented AHGBS-9 WL final 447 value as if it were gausslet-only;
- Be2/Cr2 artifact routes.

You may inspect:

- `tmp/work/wl_old_nested_vs_decomposed_he_timing_probe.jl`
- `tmp/work/wl_old_nested_fixed_block_timeg_small_he.txt`
- `tmp/work/wl_old_public_qw_timeg_small_he.txt`
- old blurb log responses around WL reproduction, especially pass 178's cited `response.007.md` and `response.008.md` if present
- current WL-related source if needed to identify the old nested/QW constructor

If a matching result is already present in logs/artifacts, report it with the
source path and do not recompute. If not, run the smallest local ignored probe
needed to compute:

```text
WL gausslet-only final dimension
WL H1 energy
WL two-body/J or RHF electron-electron component
WL RHF total energy
geometry/mapping parameters
shell inventory
whether AHGBS residuals were excluded
elapsed time
```

Then compare only at the diagnostic level:

```text
PQS 419 RHF total - WL 419 gausslet-only RHF total
PQS H1 - WL H1, if both available
PQS two-body/J - WL two-body/J, if conventionally comparable
```

Interpretation constraints:

- Do not call either value a final physics endpoint.
- Do not tune PQS or WL in this pass.
- Do not run HFDMRG/DMRG/CR2.
- Do not edit tracked source/test/generator/artifacts.
- Do not revive Be2/Cr2 artifact work.

Line-budget/diff rule:

- `git diff --numstat -- src test` should be empty.
- If tracked `src`/`test` changes become necessary, stop with `ATTENTION.md`
  rather than editing.

Validation:

```text
<only the local ignored probe command if recomputation is needed>
git status --short --branch
git diff --numstat -- src test
```

Report:

- whether the WL 419 gausslet-only result was found or newly computed;
- exact source/probe path;
- WL shell inventory and final dimension;
- WL H1, two-body/J, RHF total;
- PQS-WL diagnostic deltas;
- explicit statement that supplemented WL 447 is separate;
- no tracked source/test changes.

Write the result to `.agent_handoffs/response.183.md` and copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.183.md
```

-- repo-manager@macmini
