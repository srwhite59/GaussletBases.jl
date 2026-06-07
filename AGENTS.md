# Agent Notes

These repo-local instructions are meant to reduce friction and keep future
agent work aligned with the current engineering contract.

## Agent identity and signoff

When a `GaussletBases` agent is asked to reread startup docs, it should recover
its role identity and handback signature from these rules plus the current user
or manager assignment.

Default repo role signoffs are:

- manager handbacks: `-- repo-manager@<host>`
- implementation doer handbacks: `-- repo-doer@<host>`

Use a more specific role name only when the startup prompt, role identity file,
or user explicitly gives one. For example, an informal "repo-doer2" terminal is
still normally the implementation doer unless it has been assigned a distinct
signature.

Sign any message that hands control back to the user: final handoffs,
close-outs, ready-for-review reports, gate completions, audit results, and
safe-to-resume judgments. Do not sign progress updates while continuing work
without waiting on the user. The signoff must be the final line.

After compaction, session resume, or any sign that live context has been
replaced by a summary, treat the next turn as reentry before continuing
substantive work. Reread this file, the role prompt or current assignment, and
any governing framework doc cited by the task; then briefly restate role,
current task boundary, and signoff line.

## Julia execution

Prefer one of these two launch styles for routine work:

- plain `julia --project=. ...`
- `/Users/srw/Dropbox/codexhome/cjulia ...`

`cjulia` is preferred when you want one stable launcher path. It already:

- finds the repo root
- runs with `--project=<repo-root>`
- prefers the normal user Julia setup
- does **not** force a Dropbox-local depot

Do **not** use routine commands of the form:

- `env JULIA_DEPOT_PATH=... julia ...`

unless there is a specific reason that has already been justified.

In particular:

- do **not** use Dropbox-local Julia depots
- prefer the normal user depot in `~/.julia` / `~/.juliaup`
- prefer checked-in or `tmp/work/*.jl` scripts over long inline `julia -e '...'`
  commands

Reason:

- `env ... julia ...` causes unnecessary permission friction in this
  environment
- Dropbox-local depots have previously caused operational problems

For routine Julia test/probe timing, prefer Julia-level timing rather than
wrapping commands with `/usr/bin/time`:

- `julia --project=. -e '@time include("tmp/work/script.jl")'`
- `julia --project=. -e 't = @elapsed include("tmp/work/script.jl"); println("elapsed_s=", t)'`

Do not use `/usr/bin/time` routinely. If OS-level memory data such as maximum
RSS is genuinely needed, use a stable wrapper script such as `tools/time_julia`
and get that wrapper prefix approved once instead of requesting broad
`/usr/bin/time` approval.

## Runtime environment policy

Do not create generated runtime environments inside Dropbox-synced trees.

Keep in Dropbox:

- source code
- docs and handoffs
- configs
- `Project.toml`, `Manifest.toml`, `requirements.txt`, `pyproject.toml`,
  lockfiles, and small setup scripts

Keep out of Dropbox:

- Python `venv/`, `.venv/`, and `work/venvs/`
- Julia depots, `.julia_depot/`, `tmp/julia_depot/`, package artifacts, and
  compiled caches
- binary/runtime installs and package-manager caches

Use machine-local locations such as `~/.venvs`, `~/.julia`, or
`~/dmrgtmp/{envs,julia_depots,artifacts,compiled}`. If a special environment is
needed, create it outside Dropbox and document the activation command.

## Julia style

Use `JuliaStyle.md` as the local guide for small Julia cleanup passes and new
code edits. These preferences are low priority relative to correctness and
contract clarity; do not churn numerical kernels solely for style.

## Structured state, staged metadata, and test runtime rules

Short commandment:

- do not create flat field clouds for new route concepts
- if a concept has internal structure, carry a compact module-owned object and
  expose only a small summary or fingerprint
- before adding more than three related fields to a staged object, stop and
  define a compact object or summary
- before copying the same field group across two stages, stop and carry the
  object instead
- do not compare large staged metadata objects with `==` or `===`
- compare compact fingerprints or summaries: statuses, counts, keys, kinds,
  booleans, materialization flags, and short missing-reason tuples
- choose the smallest test that validates the edit
- before running any test expected to take more than 60 seconds, explain why it
  is necessary
- time routine Julia tests/probes with Julia-level timing such as `@elapsed` or
  `@time`, not broad `/usr/bin/time` wrappers
- do not use `test/nested/cartesian_pair_stage_low_order_policy_runtests.jl` as
  a per-pass gate; it is an integration gate

Detailed policy lives in
`docs/src/developer/test_suite_reorganization_plan.md`.

## Final-basis overlap policy

For final working bases that are intended to be orthonormal:

- self-overlaps `S_AA` and `S_BB` are diagnostic/debug objects only
- do **not** use them as normal downstream working data
- do **not** build generalized-overlap transfer logic on the final working-basis
  path

Normal final-basis transfer should use only the cross overlap:

- `S_BA = <B|A>`
- transfer with `C_B = S_BA * C_A`

Use self-overlaps only for:

- construction checks
- assertions
- difficult bug diagnostics
- or one final cleanup step if orthonormality is not yet good enough

See also:

- `docs/src/developer/numerical_contracts.md`
- `DESIGN.md`

## Code cleanup policy

When a contract changes:

- do not only patch the main implementation
- also remove or simplify stale helpers, diagnostics, benchmarks, and docstrings
  that still encode the old idea

Two different costs matter:

- code volume
- conceptual drift

Conceptual drift is worse.

Small code that preserves a wrong contract is more dangerous than large code
that is conceptually clean.

Default rule:

- if a path is no longer part of the intended contract, delete it or quarantine
  it clearly as debug/reference-only

See also:

- `docs/code_bloat_and_wrong_contract_cleanup_note.md`

## Basis bundle policy

It is acceptable for a basis bundle to carry:

- basis-defining core data
- optional basis-attached helpful sidecars

Examples of acceptable sidecars:

- factorized parent structure
- shell/block grouping
- cached axis tables
- ordering/permutation metadata

But those sidecars must remain:

- clearly auxiliary
- reproducible from the basis construction
- not confused with the basis definition itself

## Calculation and analysis preferences

For nontrivial calculations, benchmarks, transfer studies, and solver runs:

- keep coarse `TimeG` timing enabled by default
- do **not** turn coarse timing off for expensive runs unless there is a
  specific reason
- save at least a coarse timing record as a durable artifact when practical

The minimum coarse timing summary should prefer user-facing phase labels such
as:

- build nested basis
- assemble operators / Hamiltonian
- transfer orbitals or states
- Hartree-Fock / SCF solve
- total

Do not rely only on terminal scrollback for timing. If a run matters, preserve
the coarse timing summary in a file, note, or other durable artifact.

## Performance review policy

Performance is part of correctness for public APIs, algorithmic routes, and
advertised reference utilities. A change is not fully validated until CPU time
and allocation behavior have been checked against a reasonable scale model.

For nontrivial coding work, do not jump directly from "what result is needed"
to "smallest patch that passes a test." First do a short design pass that:

- defines the computation, operator, or data transform being implemented
- identifies likely cost centers and expected scaling
- decides what should be tabulated, cached, reused, or contracted early
- checks whether existing kernels or contracts should be reused
- chooses an implementation seam that is reviewable and not needlessly
  duplicated

For performance-sensitive changes, reports and handoffs should include:

- computation definition and performance/organization strategy
- performance category: production route, reference-only but usable, or
  diagnostic/prototype
- expected CPU and memory scaling
- representative fixture size
- measured time and allocations or memory footprint
- durable timing/allocation artifact path when practical
- explicit decision: ready, needs optimization, or prototype-only

Do not call a public feature ready based only on tiny correctness tests. If
performance was not measured, say that explicitly and treat it as remaining
validation work.

See also:

- `docs/src/developer/performance_review_contracts.md`

## Result-reporting language

When reporting calculation or analysis results back to the user:

- prefer user-facing language over internal implementation language
- check labels and summaries for repo-internal shorthand before sending them
- if a new term is necessary, define it briefly the first time it appears
- if an internal code label is not directly meaningful to the user, translate it
  into a plain-language description

Default rule:

- summaries should be understandable without requiring the user to know repo
  internals, Julia helper names, or private shorthand
