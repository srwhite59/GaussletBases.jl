Purpose:

Shrink residual test bloat now that the decomposed WL He and MWG/GTO contracts
are clearer. This is a cleanup/deletion pass, not a feature pass.

Why now:

Recent physics tests exposed real convention bugs and route holes, while broad
helper/status assertions mostly increased carrying cost. Some test bloat has
leaked back in, especially around report-field policing and helper-vocabulary
coverage. We need to keep the valuable scientific/module contracts and remove
or demote assertions that no longer protect a live route contract.

Current state:

- The default nested runner is better organized than before, but it may still
  include compact internal audit coverage that preserves old helper vocabulary.
- The He acceptance test is high-value because it exercises decomposed WL He
  physics and corrected IDA convention.
- The He acceptance test still appears to assert too many report/status/timing
  fields and repeated nonclaim flags.
- Current long-term policy is scientific/workflow acceptance tests plus compact
  module-contract tests, not broad helper-vocabulary tests.

Exact task:

Audit and shrink active nested test bloat with a narrow deletion-oriented pass.

Start with:

- `test/nested/cartesian_wl_gausslet_he_atom_acceptance_runtests.jl`
- `test/nested/runtests.jl`
- `test/nested/integration_runtests.jl`, only if needed to understand runner
  boundaries

For the He acceptance test, keep the live scientific contract:

- H1 energy and error;
- 1s self-Coulomb / J and error;
- RHF one-electron, electron-electron, and total energy;
- density trace and electron count;
- overlap rank/conditioning;
- matrix symmetry/finite sanity where it localizes real failures;
- key IDA weight-boundary convention;
- anti-fallback boundary:
  - no full-parent CPB;
  - no direct Cartesian product assembly;
  - no ordinary Cartesian IDA path;
- basic dimensions:
  - retained dimension;
  - unit count;
  - pair count.

Remove or demote where no longer needed:

- exhaustive route/status fields;
- detailed pair-factor metadata assertions;
- timing nonnegativity assertions unless a timing value is an active
  performance contract;
- repeated nonclaim flags;
- helper-name and helper-vocabulary checks;
- broad report-shape checks;
- assertions that preserve old development scaffolding rather than the live
  physics contract.

For `test/nested/runtests.jl`, audit any inline "owned-unit coverage" or helper
coverage blocks. If they still protect a live module contract, leave them and
say why. If they mainly preserve old helper vocabulary, shrink or move them out
of the default runner.

Trust boundary:

This pass should delete or simplify tests. Do not add new tests unless one
replaces/shrinks older coverage and protects a live contract.

Do not change production source code unless a test cleanup exposes a trivial
stale comment or name that must be adjusted. If source code appears wrong,
stop and report instead of fixing it in this cleanup pass.

Do not change physics baselines.

Do not touch GTO/MWG implementation paths in this pass.

Do not run broad slow integration gates as routine validation.

Validation:

- run the specific test file(s) touched;
- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`;
- do not run `test/nested/integration_runtests.jl` unless you move runner
  boundaries there and explain why.

Decision rules:

- If a removed assertion would be the only guard for a non-obvious live
  numerical convention, keep a smaller assertion instead of deleting it.
- If a test block is unclear, prefer a small comment/classification over adding
  more assertions.
- If shrinking the He acceptance test would require changing production code,
  stop and report.
- If the cleanup would become broad, stop after the He test and report the next
  target separately.

Deletion/shrinkage report required:

- number of test lines deleted/simplified;
- what old helper/status/metadata checks were removed;
- what live physics/module-contract checks remain;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any test was added, and if so what older coverage it replaced or
  shrank;
- any remaining stale or duplicate test surfaces to retire next.

Report back:

- files changed;
- net diff stat;
- validation run;
- what was deleted or simplified;
- what remains intentionally protected;
- any remaining cleanup target.
