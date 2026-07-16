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
- bloat cleanup handbacks: `-- bloat-fixer@<host>`

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

## Long-running job monitoring

For long Julia/Python runners, avoid broad process scans and piped shell probes
that trigger approval stalls.

Preferred pattern:

- poll the active tool/session directly when possible
- make long runners print `getpid()` at startup and write a small PID file next
  to the runtime log or output artifact
- if process inspection is needed, use an exact PID command such as
  `ps -p <pid> -o pid,etime,rss,vsz,pcpu,pmem,command`
- do not use broad scans such as `ps ... | rg ...` when a PID file or active
  session exists
- for progress checks, read known logs directly with `tail`, `sed`, `wc`, or
  similar simple file commands
- if a reusable runner lacks PID/log reporting, prefer patching that runner to
  report PID and progress rather than repeatedly using ad hoc process searches

Long-run handbacks should report the log path, PID path if any, and the exact
polling command used or recommended.

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

Every line of source, test, documentation, compatibility glue, metadata, and
adapter code has carrying cost. New code should earn that cost by protecting a
live contract, improving clarity, reducing duplication, improving performance,
or enabling a current workflow. Do not optimize for the fewest lines when
clarity or numerical safety would suffer, but prefer net simplification when the
benefit is otherwise equal.

## Blurb style

Use `BlurbStyle.md` when drafting repo-manager to repo-doer blurbs. Blurbs
should include exact known code surfaces, explicit exclusions, decision rules,
and reporting requirements rather than expecting the doer to reconstruct the
manager's searches or reasoning.

Reports back to a human should use two layers:

1. Plain-language result first: what problem was solved, what physical object
   or calculation changed, what result changed, what remains uncertain, and
   what should happen next.
2. Technical evidence second: files, functions, dimensions, timings,
   validation commands, artifacts, and commit IDs.

Do not make the user infer the physical meaning from repo-internal vocabulary.
For example, say "the residual basis is selected separately on each atom"
before naming an owner-local residual-selection helper.

## Manager running log

For Cartesian/PQS manager-led work, read:

- `docs/src/developer/pqs_manager_running_log.md`

before drafting a blurb, accepting a pass, or resuming after compaction.

The live running log is a bounded manager-level decision ledger. It records
strategic interpretation, long-term goals, current medium-term goals,
guardrails, and remaining blockers. It does not replace per-pass doer
responses, baton reviews, or `state.md`. Older volumes are linked from its
archive index and are task-gated archaeology, not startup reading.

After each accepted substantive manager-reviewed pass, append a compact but real
entry to the running log. The default is not a two-sentence tick; use roughly
100-250 words or a short bullet list that preserves strategic interpretation.
The entry should include:

- pass number/title;
- accepted commit(s);
- short summary;
- validation actually used by doer and manager;
- goal advancement using LT/MT labels;
- medium-goal update, if any;
- risk or guardrail;
- remaining blocker or next step;
- line-count/complexity note when relevant.

For purely mechanical passes with no strategic change, a 1-3 sentence tick is
allowed, but it should explicitly say "no strategic change" and name the
current goal or guardrail it leaves unchanged.

Every 5 accepted Cartesian/PQS passes, add a medium-term goal checkpoint. The
checkpoint should classify current MT goals as active, completed, blocked,
stale, or needing refinement, and update MT wording when evidence warrants it.

Every 10-20 accepted passes, or after a major correction, add a strategic
compression entry summarizing durable decisions, stale stories, false starts,
and next lane direction. Follow the live log's rotation policy: archive old
entries verbatim, retain the recent decision window, and refresh the current
strategic preamble. Rotation must not rewrite or silently summarize accepted
history.

Do not duplicate the doer response. Use the running log to preserve strategic
interpretation and prevent drift.

## Cartesian Hamiltonian producer startup

For Cartesian Hamiltonian producer work, normal startup reading is:

- `docs/src/developer/designs/cartesian_hamiltonian_producer/README.md`
- `docs/src/developer/designs/cartesian_hamiltonian_producer/current.md`
- `docs/src/developer/designs/cartesian_hamiltonian_producer/invariants.md`

Then read only the assigned ID entry in
`docs/src/developer/designs/cartesian_hamiltonian_producer/registry.md` and its
linked subsystem contract. Before numerical implementation, also consult
`docs/src/developer/algorithm_implementation_index.md`. Do not traverse the
full registry or every subsystem page as startup reading.

The full historical design and review rounds remain available under:

- `docs/src/developer/designs/cartesian_hamiltonian_producer/history/`
- `docs/src/developer/designs/cartesian_hamiltonian_producer/reviews/`

Those history/review files are not normal startup reading unless a task
explicitly asks for design-history review or archaeology.

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
- do not store basis-size, shell-size, unit-size, pair-size, center-size, or
  route-inventory collections as variable-size `Tuple(...)` or runtime-keyed
  `NamedTuple` objects; use vectors, indexed/lazy views, or compact summaries
  unless the tuple shape is mathematically fixed and tiny
- `NamedTuple` is allowed for small stable records, explicit input
  normalization, and final serialization bundles, but not as a substitute for
  construction-stage objects or route plans with identity
- do not let a `NamedTuple` value crossing a hot construction-stage boundary
  have a concrete type that encodes `q`, molecule size, retained-unit count,
  source-mode count, shell count, candidate count, or other basis-size
  inventory dimensions
- if a record crosses a construction-stage boundary and carries nontrivial
  algorithmic data, prefer a small concrete struct with vector-backed fields
  over nested `NamedTuple` carriers
- choose the smallest test that validates the edit
- before running any test expected to take more than 60 seconds, explain why it
  is necessary
- time routine Julia tests/probes with Julia-level timing such as `@elapsed` or
  `@time`, not broad `/usr/bin/time` wrappers
- do not use `test/nested/cartesian_pair_stage_low_order_policy_runtests.jl` as
  a per-pass gate; it is an integration gate

Detailed policy lives in
`docs/src/developer/test_suite_reorganization_plan.md`.

## Test scope and deletion bias

Tests are code. They carry runtime cost, maintenance cost, and conceptual cost.
Do not treat more tests as automatically better.

Prefer this test hierarchy:

- Scientific endpoint tests: highest value. Prefer checks that exercise real
  user workflows and physically meaningful results, such as atom/molecule
  energies, matrix symmetries, convergence, and known reference values.
  Endpoint tests should still be representative and bounded; they are high value
  as acceptance gates, not default per-edit gates.
- Compact module-contract tests: useful for active kernels, public/module-level
  constructors, and boundaries likely to be edited.
- Oracle/reference tests: keep only when they validate a live replacement path
  or a difficult numerical convention.
- Development scaffolding tests: delete or quarantine once the transition they
  supported is complete.
- Exhaustive metadata tests: avoid unless the metadata itself is the active
  contract.

Before adding a new test, state what non-obvious bug it would catch and why an
existing endpoint, smoke, or module-contract test would not catch or localize
that bug well enough.

During cleanup or retirement work:

- Prefer deleting stale helpers and their tests over preserving them through
  adapters.
- If a path is obsolete and only tests call it, delete both the path and the
  test pressure.
- If deletion is blocked by live source callers, report the exact callers and
  stop rather than adding a bridge by default.
- Net line count should usually go down for retirement tasks. New scaffolding is
  justified only when it replaces more stale code than it adds.
- Do not add tests that mainly preserve old helper names, route-shadow
  vocabulary, all-pairs inventory details, or transitional metadata flags.
- Do not require slow old tests to pass merely because they used to protect the
  stale path being deleted. For deletion work, prefer validating active callers,
  replacement endpoints, package load/import, and the smallest surviving
  contract smoke over reviving obsolete test shape.

Stable code that is not expected to change soon usually needs only a small smoke
or endpoint check, not all development-era blocked-path and internal-vocabulary
tests.

When uncertain during cleanup, prefer a smaller active-contract smoke test plus
a real downstream endpoint check over preserving exhaustive development-era
tests.

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

## Bloat-fixer role

New bloat-fixer agents should be started only by the user or under
archive-manager/user coordination. The repo-local `bloat-fixer` role is
defined in
`docs/src/developer/bloat_fixer_role_memo_2026-06-30.md` and should read that
memo at startup when archive-manager or the user starts it.

`bloat-fixer` is not `repo-manager` and not the implementation `repo-doer`.
Repo-manager classifies targets as stable, active/new, or scientifically
sensitive, assigns exact cleanup surfaces, reviews every diff, and decides
whether to commit or push.

`bloat-fixer` may remove or simplify carrying-cost code inside manager-approved
stable contracts: redundant internal assertions, stale transition checks, dead
helpers/tests/probes, duplicated local boilerplate, and repeated metadata
plumbing. It must not change scientific contracts, numerical thresholds,
defaults, public APIs, artifact schemas, route semantics, or active-lane
behavior unless repo-manager gives an explicit close-supervision blurb.
Bloat-fixer must not create new abstractions to remove old ones unless
repo-manager explicitly approves that replacement and the net result is
smaller.

Durable cleanup rule:

- keep checks that prevent silent wrong science
- keep public/input and artifact boundary validation
- keep rank, near-singular metric, orthogonality, identity, and symmetry checks
  where computation could otherwise continue with useless results
- delete mature internal checks that only turn an inevitable Julia or linear
  algebra crash into a nicer crash
- delete tests or diagnostics that preserve obsolete helper names, stale
  transition vocabulary, or inactive metadata shape

## Physics-target work discipline

For Cartesian and other active algorithmic work, physics targets are the
default permission source for new code. Before adding a helper, adapter,
metadata field, summary object, or test, state which live physics target or
active module contract it advances.

Before running Cartesian/PQS endpoint probes or drafting blurbs for atom,
diatomic, screened-Hartree, EGOI, residual-GTO, or Cr/Cr2 work, read the short
operational memory sheet:

- `docs/src/developer/cartesian_pqs_operational_facts.md`

It records easy-to-forget numerical facts such as the standard `core_spacing`
ladder, driver-style diatomic padding, capture due diligence, and
screened-Hartree energy accounting. It is not design authority; use it to avoid
stale probe helpers and bad fixtures.

For Cartesian/PQS implementation-doer work, the terminal due-diligence report
is mandatory review material for every endpoint, energy, residual, injection,
screened-Hartree, EGOI, Be2/Cr2-style, or curve test/probe. Doer reports must
explicitly state that the due-diligence report was inspected and summarize the
relevant parent bounds, axis counts, padding/radius, final dimension, retained
counts, shell/slab topology, and warning flags. If the report is unavailable,
the doer must say so and explain what equivalent construction facts were
checked before interpreting numerical results.

Agents touching Gausslet semantics should use
`docs/src/developer/architecture/gausslet_algorithm_refresher.md` as the
short operational memory refresh for Qiu-White, White-Lindsey, PQS, COMX,
PGDG, angular-style injection, IDA/MWG, identity sectors, and common category
mistakes. Use
`docs/src/developer/architecture/gausslet_methods_fundamentals.md` as a
fresh-start conceptual reference when deeper paper-centered context is needed.
It is a paper-centered onboarding packet for
ordinary gausslets, IDA, COMX, Qiu-White hybrids, White-Lindsey nesting, PGDG,
radial/angular variants, MWG, and PQS source boxes. If that packet and a live
repo doc disagree, report the conflict instead of improvising.

The following answers are not enough without explicit manager approval:

- makes the architecture more complete
- adds coverage
- preserves old behavior
- might be useful later

Prefer tasks phrased as:

- make old surface X unnecessary for physics target Y

over tasks phrased as:

- add helper Z

For each coding handoff, include a compact deletion/shrinkage result:

- deleted:
- simplified:
- quarantined:
- not deleted because:
- exact remaining caller/blocker:

This does not require every pass to delete code, but it does require every pass
to account for carrying cost and stale-contract pressure.

Short physics target cards are encouraged at the start of a pass. They should
state:

- target
- physics endpoint
- allowed implementation surface
- forbidden surfaces
- success condition

Keep target cards short and task-local. Do not create a new design-note layer
unless the manager explicitly asks for one.

## Cartesian Hamiltonian producer design authority

For Cartesian Hamiltonian producer source work, use the compact current
authority path:

- `docs/src/developer/designs/cartesian_hamiltonian_producer/README.md`
- `docs/src/developer/designs/cartesian_hamiltonian_producer/current.md`
- `docs/src/developer/designs/cartesian_hamiltonian_producer/registry.md`
- `docs/src/developer/designs/cartesian_hamiltonian_producer/invariants.md`

Before Cartesian/PQS numerical implementation, also check
`docs/src/developer/algorithm_implementation_index.md` for existing kernels,
optimization notes, and oracle/reference paths. This index is navigation, not
new algorithm authority; use it to avoid reimplementing known optimized
patterns.

The full June 2026 design and reviews are historical material under
`docs/src/developer/designs/cartesian_hamiltonian_producer/history/` and
`docs/src/developer/designs/cartesian_hamiltonian_producer/reviews/`; do not
use them as normal startup reading.

The machine-readable
`docs/src/developer/designs/cartesian_hamiltonian_producer/authority.toml` is
the deny-by-default record-level authority for producer work. `registry.md`
and the marked whitelist below are deterministic checked views; they grant no
independent authority. The linked subsystem document owns the numerical and
behavioral contract. `current.md`, `invariants.md`, this surrounding policy,
and canonical contracts may restrict machine-authorized work but may never add
or broaden an ID, grant, surface, path, or scope.

An ID must be authorized by `authority.toml` and appear in the derived marked
block before it permits source-bearing work. `current.md` need not enumerate
every active ID; silence there is neutral. If the machine authority, generated
views, an explicit current-status statement, or canonical contract disagree,
make no producer edit and request a docs-only reconciliation. Historical,
superseded, rejected, measurement-only, and completed-retirement IDs may remain
in the generated registry without appearing in the execution whitelist.

A validation/test ID belongs in this whitelist only when its registry entry
names an exact committed or exact planned test path, or an explicit shared
test owner plus path. One-time smokes and completed evidence grant no work.
Ignored probes, reports, and docs/evidence maintenance remain outside this
execution whitelist even when their historical ID stays in the registry.

<!-- BEGIN CARTESIAN HAMILTONIAN PRODUCER EXECUTION WHITELIST -->
> **Generated authority view. Do not edit this block.**
> Source: `docs/src/developer/designs/cartesian_hamiltonian_producer/authority.toml`.
> Authority SHA-256: `3206229afb2d505441b5118f28866672faee0fda0db39a13b1347fab6e53b7a4`.

Cartesian Hamiltonian producer source work is currently authorized only for
these approved design IDs:

- `HP-CGRB-FILE-01`
- `HP-CGRB-FN-01`
- `HP-CGRB-FN-02`
- `HP-CGRB-NN-FILE-01`
- `HP-CGRB-NN-FN-01`
- `HP-CGRB-NN-TEST-01`
- `HP-CGRB-NN-WIRE-01`
- `HP-CGRB-TEST-01`
- `HP-CGRB-WIRE-01`
- `HP-COMP-ANGBOX-FN-01`
- `HP-COMP-ATOMBOX-FN-01`
- `HP-COMP-ATOMBOX-TEST-01`
- `HP-COMP-BASEDIAT-FN-01`
- `HP-COMP-BASEDIAT-TEST-01`
- `HP-COMP-FACEPROD-FN-01`
- `HP-COMP-NS-FN-01`
- `HP-COMP-NS-TEST-01`
- `HP-COMP-NSCORE-FN-01`
- `HP-COMP-NSCORE-TEST-01`
- `HP-COMP-SHELLGEOM-DIAT-FN-01`
- `HP-COMP-SHELLGEOM-FN-01`
- `HP-COMP-SUPPATOM-FN-01`
- `HP-COMP-SUPPATOM-TEST-01`
- `HP-COMP-SUPPWL-FN-01`
- `HP-COMP-SUPPWL-TEST-01`
- `HP-COMP-THINSLAB-FN-01`
- `HP-COMP-THINSLAB-META-FN-01`
- `HP-COMP-WLDIAT-FN-01`
- `HP-COMP-WLNS-FN-01`
- `HP-COMP-WLNS-TEST-01`
- `HP-CONTRACT-VEC-FN-01`
- `HP-DRV-ATOM-CLEAN-01`
- `HP-DRV-ATOM-FN-01`
- `HP-DRV-ATOM-WIRE-01`
- `HP-DRV-FILE-01`
- `HP-DRV-FN-01`
- `HP-DRV-INV-FN-01`
- `HP-DRV-NEST-FN-01`
- `HP-DRV-NEST-WIRE-01`
- `HP-DRV-SHELLDD-FN-01`
- `HP-DRV-STAGE-FN-01`
- `HP-DRV-STAGE-WIRE-01`
- `HP-FILE-01`
- `HP-FN-00`
- `HP-FN-01`
- `HP-FN-02`
- `HP-FN-03`
- `HP-FN-04`
- `HP-FN-05`
- `HP-HAM-MANIFEST-FN-01`
- `HP-HAM-MANIFEST-SRC-FN-01`
- `HP-HAM-MANIFEST-SRC-TEST-01`
- `HP-HAM-MANIFEST-TEST-01`
- `HP-MCOMX-DRV-FN-01`
- `HP-MCOMX-FILE-01`
- `HP-MCOMX-FN-01`
- `HP-MCOMX-OBJ-01`
- `HP-MCOMX-TERM-FN-01`
- `HP-MCOMX-WIRE-01`
- `HP-NEST-ART-FN-01`
- `HP-NEST-ART-TEST-01`
- `HP-OBJ-01`
- `HP-OBJ-02`
- `HP-PQS-ASPECTSHELL-FN-01`
- `HP-PQS-ATOMREF-PACKET-FN-01`
- `HP-PQS-ATOMREF-PACKET-TEST-01`
- `HP-PQS-COULOMB-ACCURACY-FN-01`
- `HP-PQS-COULOMB-ACCURACY-TEST-01`
- `HP-PQS-MAP-SFACTOR-FN-01`
- `HP-PQS-MAP-SFACTOR-TEST-01`
- `HP-PQS-PARENT-GDIRECT-FN-01`
- `HP-PQS-PARENT-GDIRECT-TEST-01`
- `HP-PQS-PRF-FN-01`
- `HP-PQS-PRF-INJECT-COMP-FN-01`
- `HP-PQS-PRF-INJECT-COMP-TEST-01`
- `HP-PQS-PRF-INJECT-INTERACT-FN-01`
- `HP-PQS-PRF-INJECT-INTERACT-TEST-01`
- `HP-PQS-PRF-TEST-01`
- `HP-PQS-SCREEN-HARTREE-CORR-FN-01`
- `HP-PQS-SCREEN-HARTREE-CORR-TEST-01`
- `HP-PQS-SHELLQ-OVERRIDE-FN-01`
- `HP-PQS-SHELLQ-OVERRIDE-TEST-01`
- `HP-R1-ART-01`
- `HP-R1-ATOM-FN-01`
- `HP-R1-ATOM-TEST-01`
- `HP-R1-ATOM-WIRE-01`
- `HP-R1-CORE-FN-01`
- `HP-R1-FILE-01`
- `HP-R1-FN-01`
- `HP-R1-TEST-01`
- `HP-R1-WIRE-01`
- `HP-R3-ART-01`
- `HP-R3-FN-01`
- `HP-R3-FN-02`
- `HP-R3-FN-03`
- `HP-R3-OBJ-01`
- `HP-R3-TEST-01`
- `HP-R3BASE-DRV-WIRE-01`
- `HP-R3BASE-FN-01`
- `HP-R3BASE-TEST-01`
- `HP-R3GG-FN-01`
- `HP-R3GG-TEST-01`
- `HP-R3U-FILE-01`
- `HP-R3U-FN-01`
- `HP-R3U-TEST-01`
- `HP-R3U-WIRE-01`
- `HP-R3U-ZDI-FN-01`
- `HP-R3U-ZDI-TEST-01`
- `HP-R3U-ZDI-WIRE-01`
- `HP-R3UN-FN-01`
- `HP-R3UN-TEST-01`
- `HP-RAW-SRCMODE-FN-01`
- `HP-REP-XGTO-IMPORT-FN-01`
- `HP-REP-XGTO-IMPORT-TEST-01`
- `HP-REP-XGTO-PROTECT-SIDECAR-FN-01`
- `HP-REP-XGTO-PROTECT-SIDECAR-TEST-01`
- `HP-RG-CUTOFF-FN-02`
- `HP-RG-CUTOFF-TEST-02`
- `HP-RG-FILE-01`
- `HP-RG-FN-01`
- `HP-RG-FN-02`
- `HP-RG-FN-03`
- `HP-RG-FN-04`
- `HP-RG-INJECT-FN-01`
- `HP-RG-NUMCOMP-FN-01`
- `HP-RG-NUMCOMP-TEST-01`
- `HP-RG-OBJ-01`
- `HP-RG-OCC-FIRST-INJECT-FN-01`
- `HP-RG-OCC-FIRST-INJECT-TEST-01`
- `HP-RG-ORTHO-FN-01`
- `HP-RG-ORTHO-TEST-01`
- `HP-RG-PROTECT-ADDREF-FN-01`
- `HP-RG-PROTECT-ADDREF-TEST-01`
- `HP-RG-PROTECT-ART-FN-01`
- `HP-RG-PROTECT-ARTLOC-FN-01`
- `HP-RG-PROTECT-EGOI-FN-01`
- `HP-RG-PROTECT-EGOI-TEST-01`
- `HP-RG-PROTECT-INJECT-FN-01`
- `HP-RG-PROTECT-LADDER-BUNDLE-FN-01`
- `HP-RG-PROTECT-ONEBODY-FN-01`
- `HP-RG-TEST-01`
- `HP-RG-WIRE-01`
- `HP-RHO0-MIXH-FEXACT-FN-01`
- `HP-RHO0-MIXH-GAAA-FN-01`
- `HP-RHO0-MIXH-GAAA-TEST-01`
- `HP-RHO0-MIXH-GG-FN-01`
- `HP-RHO0-MIXH-GG-TEST-01`
- `HP-ROUTE-INV-FN-01`
- `HP-ROUTE-RECIPE-FN-01`
- `HP-ROUTE-STAGE-CARRIER-FN-01`
- `HP-ROUTE-STAGE-TYPE-FN-01`
- `HP-WIRE-01`
- `HP-WLDIAT-COMPACT-FN-01`
- `HP-WLDIAT-PARITY-FN-01`
- `HP-WLTERM-FILE-01`
- `HP-WLTERM-FN-01`
- `HP-WLTERM-TEST-01`
- `HP-WLTERM-WIRE-01`

<!-- END CARTESIAN HAMILTONIAN PRODUCER EXECUTION WHITELIST -->

No other production surface may be added in this lane without a prior
documentation-only design amendment. This includes new structs, persistent
result shapes, modules/files, stage-return fields, metadata keys, status or
blocker symbols, report/artifact fields, committed probes/tests, and cross-file
or module-owned helpers.

The marked IDs are lexically sorted. Their order carries no lifecycle,
dependency, priority, or implementation-sequence meaning.

The hand-written guardrails below may restrict an approved ID further, but
they cannot add or broaden an ID, grant, surface, path, or scope.

Private file-local helpers are allowed only when they implement the assigned
canonical contract, create no persistent shape or vocabulary, stay within the
approved source surface, and are reported in the implementation handoff.

## Hard Cartesian/PQS anti-bloat gate

Until explicitly relaxed by the user, repo-manager must reject Cartesian/PQS
commits that fail this gate. This section applies to active work touching:

- `src/cartesian*`
- `src/pqs*`
- `src/*source_box*`
- `src/*route_driver*`
- `bin/cartesian_ham_builder.jl`
- `tools/cartesian*`
- tests that exercise these paths

No blocker-only commits:

- Do not accept a source commit whose main achievement is only exposing a later
  blocker, adding a preflight/status vocabulary, carrying a planned object
  farther through the pipeline, adding report/summary/metadata fields, or
  proving a future implementation is needed.
- A source commit must produce a real numerical object used by a physics
  endpoint, fix a demonstrated wrong result, delete obsolete code, or reduce
  measured cost/complexity.
- If a blocker is discovered but not crossed within the task constraints, make
  no commit and report the blocker.

Before coding, require this compact target card:

```text
Target:
Physics endpoint:
Allowed files:
Forbidden files/surfaces:
Must delete or simplify:
Forbidden additions:
Success condition:
Validation:
Line budget:
Failure rule:
```

The failure rule must say that if the goal cannot be met within the constraints,
the agent makes no commit and reports the obstacle.

Default forbidden additions unless explicitly approved before implementation:

- new `NamedTuple{...}` with runtime-generated keys
- new `NamedTuple` construction-stage objects, route plans, or provenance
  carriers whose concrete type encodes `q`, molecule size, retained-unit count,
  source-mode count, shell count, candidate count, or other basis-size
  inventory dimensions
- new variable-size `Tuple(...)` / `Tuple{Vararg{...}}` route inventories for
  basis-size, shell-size, unit-size, pair-size, center-size, or all-pairs data
- new `NamedTuple` wrappers around variable-size `Tuple{Vararg{...}}`
  inventories, even when the wrapper field names are stable
- new `.metadata` reads for algorithmic data
- new metadata fields that carry transforms, rules, matrices, source plans,
  dimensions, or coefficients
- new status/materialized flag clouds
- new payload structs with many untyped fields
- new broad `catch`
- new compatibility adapters
- new committed tests
- new report fields duplicating stage internals
- recursive stage returns that embed prior full stages

Metadata may contain provenance only. Summaries are disposable views, not
algorithmic APIs.

Line budgets use added source lines, not net lines:

- bug fix: at most 60 added `src` lines
- cleanup/refactor: at most 80 added `src` lines and net source decrease
  required
- new numerical capability: at most 150 added `src` lines, with replaced
  implementation deleted in the same commit
- status/preflight/report-only work: at most 20 added `src` lines

Reject over-budget patches unless the user explicitly approved the overage
before coding.

Replacement rule:

- Do not accept "new path plus old path still present" unless there is a listed
  live external caller.
- A replacement commit must add the replacement, switch live callers, delete the
  old implementation, and delete obsolete tests/probes/docstrings preserving the
  old contract.
- Git history is the archive.

Test rule:

- Do not accept new committed smoke tests in this lane unless explicitly
  approved.
- Tests should be real physics endpoints or very small numerical-kernel tests.
- Reject tests that mainly assert internal status symbols, blocker names,
  metadata flags, helper names, field presence, terminal role vocabulary, or
  all-pairs inventory details.
- Temporary probes belong in ignored `tmp/work` files and must not be
  committed.

Mechanical review comes before scientific review. For accepted Cartesian/PQS
commits, run a diff gate and paste the result or a compact summary into the
manager-log entry:

```bash
git diff --check

echo "SRC line changes:"
git diff --numstat HEAD~1..HEAD -- src bin tools test docs | sort -k1,1nr | head -40

echo "Suspicious added lines:"
git diff -U0 HEAD~1..HEAD -- src bin tools test |
  grep '^+' |
  grep -Ev '^\+\+\+' |
  grep -nE 'NamedTuple\{|Tuple\(|Tuple\{Vararg|tuple\(.*record|tuple\(.*unit|tuple\(.*pair|\.metadata|get\(.*metadata|haskey\(.*metadata|_materialized|status.*=|blocker.*=|catch$|catch err|::Any|Payload|summary.*=' || true

echo "New tests/files:"
git diff --name-status HEAD~1..HEAD |
  grep -E '^(A|R|C).*test|^(A|R|C).*tools' || true
```

If suspicious lines appear, the doer must justify each one. Default action is
reject.

Every accepted Cartesian/PQS manager-log entry must include:

- deleted:
- simplified:
- quarantined:
- not deleted because:
- exact remaining caller/blocker:
- added src lines:
- deleted src lines:
- new tests:
- new metadata/status fields:
- validation:

If `deleted` and `simplified` are empty, be suspicious.

Cleanup mode is stricter:

- no new committed tests
- no new public helper files
- no new payload structs
- no new status symbols
- no new metadata keys
- no new compatibility adapters
- no recursive stage embedding
- no report expansion

Each cleanup pass should either delete a stale surface or replace duplicate code
with a canonical existing implementation. Do not start a new framework to clean
up the old framework.

Current preferred deletion sequence:

1. Stop assembly/report objects from recursively embedding prior stages and
   duplicating route summaries.
2. Delete old route-skeleton retained/pair mirrors once the typed terminal plan
   is live.
3. Route H2 H1/J through the common H1/J helper and delete H2-local
   support-weight/raw-pair/IDA/self-Coulomb duplicates.
4. Quarantine or delete inactive pair-materialization scaffolding not used by a
   live physics endpoint.
5. Replace metadata-carried numerical data with typed fields on one canonical
   terminal/unit object, not more summaries.

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
- expand an acronym the first time it appears unless it is already common in
  the immediate conversation
- avoid unexplained internal status symbols
- write complete sentences for conclusions and recommendations
- explain the physical object before its data structure
- avoid words such as "lane", "seam", "surface", and "payload" in user-facing
  explanations unless the word is genuinely necessary; use plainer words such
  as "path", "boundary", "interface", or "object" when that is what is meant
- separate uncertainty from evidence: say what is still unknown and what
  measurement or review would resolve it

Default rule:

- summaries should be understandable without requiring the user to know repo
  internals, Julia helper names, or private shorthand

For subtle numerical code, comments should explain why the equation or
convention is being used, not narrate the loop. Prefer names whose meaning is
visible at the call site, such as `residual_occupation_cutoff`,
`owner_candidate_indices`, or `interowner_overlap`, over compact abbreviations
whose meaning must be reconstructed from surrounding context.
