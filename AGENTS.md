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

The running log is a manager-level decision ledger. It records strategic
interpretation, long-term goals, current medium-term goals, guardrails, and
remaining blockers. It does not replace per-pass doer responses, baton reviews,
or `state.md`.

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
and next lane direction. Do not prune prior running-log entries by default;
prefer append-only compression unless the user asks for archival
reorganization.

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
PGDG, IDA/MWG, identity sectors, and common category mistakes. Use
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

The list below is the deny-by-default execution whitelist for source-bearing
producer work. `registry.md` owns each ID's permission and lifecycle, while
the linked subsystem document owns its numerical/behavioral contract. An ID
must be present below and remain active in `registry.md` before it authorizes
source work. `current.md` need not enumerate every active ID; silence there is
neutral. If the whitelist, registry lifecycle, an explicit current-status
statement, and canonical contract disagree, make no source edit and request a
docs-only reconciliation. Historical, superseded, rejected, measurement-only,
and completed-retirement IDs may remain documented in the registry without
being source-authorized here.

Cartesian Hamiltonian producer source work is currently authorized only for
these approved design IDs:

- `HP-OBJ-01`
- `HP-OBJ-02`
- `HP-FILE-01`
- `HP-FN-00`
- `HP-FN-01`
- `HP-FN-02`
- `HP-WIRE-01`
- `HP-FN-03`
- `HP-FN-04`
- `HP-FN-05`
- `HP-R1-FILE-01`
- `HP-R1-FN-01`
- `HP-R1-CORE-FN-01`
- `HP-PQS-MAP-SFACTOR-FN-01`
- `HP-PQS-MAP-SFACTOR-TEST-01`
- `HP-PQS-COULOMB-ACCURACY-FN-01`
- `HP-PQS-COULOMB-ACCURACY-TEST-01`
- `HP-PQS-ATOMREF-PACKET-FN-01`
- `HP-PQS-ATOMREF-PACKET-TEST-01`
- `HP-REP-XGTO-IMPORT-FN-01`
- `HP-REP-XGTO-IMPORT-TEST-01`
- `HP-REP-XGTO-PROTECT-SIDECAR-FN-01`
- `HP-REP-XGTO-PROTECT-SIDECAR-TEST-01`
- `HP-PQS-SCREEN-HARTREE-CORR-FN-01`
- `HP-PQS-SCREEN-HARTREE-CORR-TEST-01`
- `HP-R1-WIRE-01`
- `HP-R1-ART-01`
- `HP-R1-TEST-01`
- `HP-R1-ATOM-FN-01`
- `HP-R1-ATOM-WIRE-01`
- `HP-R1-ATOM-TEST-01`
- `HP-ROUTE-RECIPE-FN-01`
- `HP-ROUTE-RECIPE-TEST-01`
- `HP-ROUTE-INV-FN-01`
- `HP-ROUTE-INV-TEST-01`
- `HP-RAW-SRCMODE-FN-01`
- `HP-RAW-SRCMODE-TEST-01`
- `HP-CONTRACT-VEC-FN-01`
- `HP-CONTRACT-VEC-TEST-01`
- `HP-ROUTE-STAGE-TYPE-FN-01`
- `HP-ROUTE-STAGE-TYPE-TEST-01`
- `HP-ROUTE-STAGE-CARRIER-FN-01`
- `HP-ROUTE-STAGE-CARRIER-TEST-01`
- `HP-WLTERM-FILE-01`
- `HP-WLTERM-FN-01`
- `HP-WLTERM-WIRE-01`
- `HP-WLTERM-TEST-01`
- `HP-R3-OBJ-01`
- `HP-R3-FN-01`
- `HP-R3-FN-02`
- `HP-R3-FN-03`
- `HP-R3-ART-01`
- `HP-R3-TEST-01`
- `HP-R3U-FILE-01`
- `HP-R3U-FN-01`
- `HP-R3U-WIRE-01`
- `HP-R3U-TEST-01`
- `HP-RG-FILE-01`
- `HP-RG-OBJ-01`
- `HP-RG-FN-01`
- `HP-RG-FN-02`
- `HP-RG-FN-03`
- `HP-RG-FN-04`
- `HP-RG-WIRE-01`
- `HP-RG-TEST-01`
- `HP-RG-ORTHO-FN-01`
- `HP-RG-ORTHO-TEST-01`
- `HP-RG-IDTOL-FN-01`
- `HP-RG-IDTOL-TEST-01`
- `HP-RG-CUTOFF-FN-01`
- `HP-RG-CUTOFF-TEST-01`
- `HP-RG-CUTOFF-FN-02`
- `HP-RG-CUTOFF-TEST-02`
- `HP-RG-INJECT-FN-01`
- `HP-RG-OCC-FIRST-INJECT-FN-01`
- `HP-RG-OCC-FIRST-INJECT-TEST-01`
- `HP-RG-PROTECT-ADDREF-FN-01`
- `HP-RG-PROTECT-ADDREF-TEST-01`
- `HP-RG-PROTECT-INJECT-FN-01`
- `HP-RG-PROTECT-INJECT-TEST-01`
- `HP-RG-PROTECT-ONEBODY-FN-01`
- `HP-RG-PROTECT-ONEBODY-TEST-01`
- `HP-RG-PROTECT-ART-FN-01`
- `HP-RG-PROTECT-ART-TEST-01`
- `HP-RG-PROTECT-ARTLOC-FN-01`
- `HP-RG-PROTECT-ARTLOC-TEST-01`
- `HP-RG-PROTECT-EGOI-FN-01`
- `HP-RG-PROTECT-EGOI-TEST-01`
- `HP-RG-PROTECT-LADDER-BUNDLE-FN-01`
- `HP-RG-PROTECT-LADDER-BUNDLE-TEST-01`
- `HP-RHO0-MIXH-GG-FN-01`
- `HP-RHO0-MIXH-GG-TEST-01`
- `HP-RHO0-MIXH-GAAA-FN-01`
- `HP-RHO0-MIXH-GAAA-TEST-01`
- `HP-RHO0-MIXH-FEXACT-FN-01`
- `HP-RHO0-MIXH-FEXACT-TEST-01`
- `HP-CGRB-FILE-01`
- `HP-CGRB-FN-01`
- `HP-CGRB-FN-02`
- `HP-CGAI-FN-01`
- `HP-CGRB-WIRE-01`
- `HP-CGRB-TEST-01`
- `HP-CGRB-NN-FILE-01`
- `HP-CGRB-NN-FN-01`
- `HP-CGRB-NN-WIRE-01`
- `HP-CGRB-NN-TEST-01`
- `HP-R3GG-FN-01`
- `HP-R3GG-TEST-01`
- `HP-R3UN-FN-01`
- `HP-R3UN-TEST-01`
- `HP-R3BASE-FN-01`
- `HP-R3BASE-TEST-01`
- `HP-R3BASE-DRV-WIRE-01`
- `HP-R3BASE-DRV-TEST-01`
- `HP-DRV-FILE-01`
- `HP-DRV-FN-01`
- `HP-DRV-TEST-01`
- `HP-DRV-NEST-FN-01`
- `HP-DRV-NEST-WIRE-01`
- `HP-DRV-NEST-TEST-01`
- `HP-DRV-ATOM-FN-01`
- `HP-DRV-ATOM-WIRE-01`
- `HP-DRV-ATOM-TEST-01`
- `HP-DRV-ATOM-CLEAN-01`
- `HP-R3U-ZDI-FN-01`
- `HP-R3U-ZDI-WIRE-01`
- `HP-R3U-ZDI-TEST-01`
- `HP-DRV-STAGE-FN-01`
- `HP-DRV-STAGE-WIRE-01`
- `HP-DRV-STAGE-TEST-01`
- `HP-DRV-INV-FN-01`
- `HP-DRV-INV-TEST-01`
- `HP-DRV-SHELLDD-FN-01`
- `HP-DRV-SHELLDD-TEST-01`
- `HP-PQS-ASPECTSHELL-FN-01`
- `HP-PQS-ASPECTSHELL-TEST-01`
- `HP-HAM-MANIFEST-FN-01`
- `HP-HAM-MANIFEST-TEST-01`
- `HP-HAM-MANIFEST-SRC-FN-01`
- `HP-HAM-MANIFEST-SRC-TEST-01`
- `HP-NEST-ART-FN-01`
- `HP-NEST-ART-TEST-01`
- `HP-COMP-WLDIAT-FN-01`
- `HP-COMP-WLDIAT-TEST-01`
- `HP-COMP-BASEDIAT-FN-01`
- `HP-COMP-BASEDIAT-TEST-01`
- `HP-COMP-SUPPWL-FN-01`
- `HP-COMP-SUPPWL-TEST-01`
- `HP-COMP-SUPPATOM-FN-01`
- `HP-COMP-SUPPATOM-TEST-01`
- `HP-COMP-ATOMBOX-FN-01`
- `HP-COMP-ATOMBOX-TEST-01`
- `HP-COMP-NS-FN-01`
- `HP-COMP-NS-TEST-01`
- `HP-COMP-NSCORE-FN-01`
- `HP-COMP-NSCORE-TEST-01`
- `HP-COMP-SHELLGEOM-FN-01`
- `HP-COMP-SHELLGEOM-TEST-01`
- `HP-COMP-SHELLGEOM-DIAT-FN-01`
- `HP-COMP-SHELLGEOM-DIAT-TEST-01`
- `HP-COMP-THINSLAB-FN-01`
- `HP-COMP-THINSLAB-TEST-01`
- `HP-COMP-THINSLAB-META-FN-01`
- `HP-COMP-THINSLAB-META-TEST-01`
- `HP-COMP-FACEPROD-FN-01`
- `HP-COMP-FACEPROD-TEST-01`
- `HP-COMP-ANGBOX-FN-01`
- `HP-COMP-ANGBOX-TEST-01`
- `HP-MCOMX-FILE-01`
- `HP-MCOMX-OBJ-01`
- `HP-MCOMX-FN-01`
- `HP-MCOMX-WIRE-01`
- `HP-MCOMX-TEST-01`
- `HP-MCOMX-TERM-FN-01`
- `HP-MCOMX-TERM-TEST-01`
- `HP-MCOMX-DRV-FN-01`
- `HP-MCOMX-DRV-TEST-01`
- `HP-COMP-WLNS-FN-01`
- `HP-COMP-WLNS-TEST-01`
- `HP-WLDIAT-COMPACT-FN-01`
- `HP-WLDIAT-COMPACT-TEST-01`
- `HP-WLDIAT-PARITY-FN-01`
- `HP-WLDIAT-PARITY-TEST-01`

No other production surface may be added in this lane without a prior
documentation-only design amendment. This includes new structs, persistent
result shapes, modules/files, stage-return fields, metadata keys, status or
blocker symbols, report/artifact fields, committed probes/tests, and cross-file
or module-owned helpers.

Private file-local helpers are allowed only when they implement the approved
Slice A/B/C/D pseudocode, create no persistent shape or vocabulary, stay within
the approved file and line budget, and are reported in the implementation
handoff.

`HP-FN-04` approves only the internal Slice C1 localized IDA matrix assembly
surface. `HP-FN-05` approves only the narrow Slice C2 construction boundary
for the existing `CartesianIDAHamiltonian`. `HP-WIRE-02` is historical: its
route-driver materialization wrapper was removed by `e2e164e9b` and is no
longer source authority. Do not restore that wrapper, its report/save
choreography, or compatibility adapters without a new docs-only amendment.

`HP-R1-FILE-01` approves only `src/cartesian_base_hamiltonian.jl`.
`HP-R1-FN-01` approves only the public `cartesian_base_hamiltonian` facade with
the approved signature. The reviewed one-center H endpoint uses explicit
public `core_spacing = 0.3` with `reference_spacing = 1.0`; resolved
`core_spacing` maps to the private one-center `parent_mapping_d`, while
`reference_spacing` remains the separate reference-grid spacing. Public `d` is
deprecated; if temporarily accepted, it must equal resolved `core_spacing`, and
z-axis H2 rejects `d`. Public `parent_mapping_d` remains unsupported.
`HP-R1-CORE-FN-01` freezes `core_spacing` as the single public physical
near-nucleus spacing after explicit input, driver default, or preset
resolution. White-Lindsey `Z` behavior is an internal mapping-shape/default
rule, not a second public knob. Driver/project defaults such as
`core_spacing = 0.3` are allowed only when visible and overrideable; once
resolved, they are ordinary explicit inputs. Routine correctness tests may
override driver physics defaults, but any asserted scalar must be tied to the
exact test input and not described as a physics-default result.
`HP-PQS-MAP-SFACTOR-FN-01` and `HP-PQS-MAP-SFACTOR-TEST-01` approve the only
current public mapping-strength exception: optional positive expert
`s_factor`, default `1.0`, with one-center
`effective_s = s_factor * sqrt(Z * core_spacing)` and mapping provenance
recording `mapping_s_factor`, `mapping_s_standard`, and
`mapping_s_effective`. This does not revive public `d`, public
`parent_mapping_d`, public `parent_mapping_Z`, route-specific mapping controls,
element defaults, automatic tuning, solver workflow, EGOI, rho0/P0, or
protected-localized convention changes. If multicenter combined-invsqrt
mapping cannot unambiguously support the same semantics, implement only the
one-center path and report the exact blocker.
`HP-PQS-COULOMB-ACCURACY-FN-01` and
`HP-PQS-COULOMB-ACCURACY-TEST-01` approve only the producer-wide
`coulomb_accuracy = :compact | :standard | :high` policy in
`docs/src/developer/designs/cartesian_hamiltonian_producer/coulomb_accuracy_policy.md`.
The default remains `:compact`; the presets are the legacy 45-term compact,
fixed analytic K60 standard, and 135-term high expansions. `:standard` is the
recommended opt-in accuracy/cost tier; `:high` remains reference-grade.
Policy, exact parameters, term count, and coefficient/exponent fingerprint
define the preset; `doacc` is only a legacy compatibility field. One resolved
`CoulombGaussianExpansion` must be carried from parent/PGDG construction
through base unit-nuclear/IDA, residual-GTO exact Coulomb-expanded blocks, and
MWG. New artifacts may add one Hamiltonian-wide compact expansion summary;
protected/ladder readback must not infer missing legacy provenance as
`:standard` or `:high`. Atomic packet RHF remains high accuracy while its
current density/self-energy and fitted-potential scaffold evaluations remain
role-qualified compact approximations. This authority narrowly allows
`src/ordinary_coulomb.jl` to own the exact fixed K60 construction without
changing compact/high bit patterns. It does not approve custom expansion
inputs or a standard/high default. It narrowly allows
`bin/cartesian_ham_builder.jl` to expose the same `coulomb_accuracy` symbol,
default `:compact`, validate `:compact | :standard | :high`, forward it in
`common_basis`, and print it; the driver must not resolve the expansion. No
other canonical driver/CLI changes, ordinary QW/legacy cleanup, solver
workflow, EGOI or screened-Hartree formula changes, or Cr2-specific behavior
are approved.
The same IDs narrowly add `src/GaussianAnalyticIntegrals.jl` and
`src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl` for algebraically
stable determinant and pairwise weighted-distance forms in
`gaussian_factor`, `gaussian_pair_factor`, and
`_factor_axis_integral`. They also allow one focused BigFloat-oracle test
in `test/core/runtests.jl`. This is not authority for clamping, exponent
truncation, scaled/log PGDG carriers, new stage objects, terminal contraction
changes, or a broad analytic-integral rewrite.
`HP-PQS-ATOMREF-PACKET-FN-01` and
`HP-PQS-ATOMREF-PACKET-TEST-01` authorize only the implemented one-center
packet facility recorded in the registry and
`docs/src/developer/designs/cartesian_hamiltonian_producer/atomic_hf_reference_packets.md`.
Unconverged packets are rejected. Fitted density/potential terms are not
protected orbitals. Density-fit `J0_G` must receive its role-qualified compact
Coulomb expansion explicitly. The stored overlap fingerprint is an exact
packet self-integrity check. Translated owner-local embedding keeps exact
atom/basis/count/owner/center/order checks but compares overlap numerically at
the unchanged `norm(..., Inf) <= 1e-10`; mapped raw-byte hash mismatch alone
is not a failure. Packet construction uses the ordinary density fit followed
by the ordinary radial potential fit. The density fit owns `E0`; the fitted
potential is an approximate `J0` evaluator whose radial, tail, matrix, and
`Tr(P0*J0_fit)-E0_fit` consistency errors are reported. Packets carrying
retired `potential_fit/moment_polish/*` provenance are rejected and require
regeneration. The retired `HP-PQS-ATOMREF-POTMOM-*` IDs are historical only
and do not authorize source work. No production correction, solver, public
default, exchange, EGOI, row-gauge rho0/P0, or Cr/Cr2 claim is approved.
`HP-REP-XGTO-IMPORT-FN-01` and `HP-REP-XGTO-IMPORT-TEST-01` own the implemented
external-GTO import and approve its narrow protected-member composition
extension. Import uses `S_FG = <F|G_external>` and `C_F = S_FG*C_G`; a
protected member uses native `S_LG = <L|G_external>` from the existing exact
handoff. External `S_GG` is validation and metric-aware capture data, not a
generalized final-basis metric. Direct imported coefficients are not
solver-orthonormalized. `HP-REP-XGTO-PROTECT-SIDECAR-FN-01` and
`HP-REP-XGTO-PROTECT-SIDECAR-TEST-01` additionally approve one standalone
native-order sidecar containing `S_LG`, direct spin imports, packet/member
fingerprints, provenance, and metric-aware capture diagnostics. The sidecar is
not a protected Hamiltonian field or ladder transfer/restart. No PySCF
dependency, raw-coefficient angular capture, raw `G_L/A_L` persistence,
Hamiltonian/`Vee` transform, `C' V C`, generalized final overlap, public
workflow, solver/HF loop, screened-Hartree/EGOI change, or Cr2 production claim
is approved.
`HP-PQS-SCREEN-HARTREE-CORR-FN-01` and
`HP-PQS-SCREEN-HARTREE-CORR-TEST-01` implement only the internal in-memory
screened direct-Hartree correction specified by
`docs/src/developer/designs/cartesian_hamiltonian_producer/screened_hartree_correction_assembly.md`
and owned by
`src/cartesian_reference_density/screened_hartree_correction.jl`, its module
wiring, and small correctness tests. The helper returns
`Delta_J0 = J0_G - Diagonal(V_IDA * q0)` and
`C = 0.5*q0'V_IDA*q0 - 0.5*E0_G` from represented, converged determinants and
same-basis fields. `Vnuc` remains Galerkin; the density fit owns `E0_G`; the
potential fit only evaluates `J0_G`; and `Delta_J0 + C` is direct
electron-electron accounting. Exact/density-fit energy identities and all
representation, finiteness, symmetry, convergence, and derivative/algebra
checks remain strict. An ordinary fitted-potential consistency error is
reported, not rejected solely because it exceeds `1e-8 Ha`.
No public workflow, corrected artifact, solver integration, exchange, EGOI,
source/interaction transform, `C' V C`, or Cr2 production claim is approved.
`HP-R1-WIRE-01` approves
only the report-free shared base constructor seam and the approved callers.
`HP-R1-ART-01` approves only the fixed `producer_provenance/` schema in the
final Hamiltonian file. `HP-R1-TEST-01` approves only
`test/driver_public/cartesian_base_hamiltonian_runtests.jl` as a standalone
integration/endpoint gate. Base R1 scope is origin-centered H, Cartesian
z-axis H2, and explicit origin-centered all-electron one-center atoms under
`HP-R1-ATOM-*`. No driver/bin/tool/report/payload/status/pair/assembly public
workflow expansion is approved, and no artifact expansion is approved except
the `HP-R1-ART-01` provenance keys.

`HP-R1-ATOM-FN-01` and `HP-R1-ATOM-WIRE-01` relax the one-center base facade in
`src/cartesian_base_hamiltonian.jl` from H-only validation to explicit
origin-centered all-electron atoms. The caller must provide vector-valued
`atom_symbols`, `nuclear_charges`, `atom_locations`, explicit `nup`/`ndn`, and
one-center basis controls including `core_spacing`. The atom symbol is
provenance only; charge, electron count, spin, basis, and ECP behavior must not
be inferred from element tables. Public charge maps to the private
White-Lindsey atomic mapping `Z`; resolved `core_spacing` maps to the private
`parent_mapping_d`; `reference_spacing`, tail spacing, and box/domain controls
remain separate.
Atoms and diatomics must share the same producer workflow after
geometry/shellification normalization. Do not add atom-only Hamiltonian
builders, parallel atom materialization, atom route/report/status payloads, or
atom-specific artifact shapes. `HP-R1-ATOM-TEST-01` approves only H regression
plus optional ignored/user-run Be or Cr one-center base atom artifact checks; it
does not approve committed non-H fixtures/tests, translated atoms, supplemented
atoms, ECP, solver workflow, new artifact keys, driver changes, or source files
outside `src/cartesian_base_hamiltonian.jl`.

`HP-ROUTE-RECIPE-FN-01` approves only family-selective route recipe cleanup in
`src/pqs_source_box_route_driver_helpers.jl` and
`src/cartesian_base_hamiltonian.jl`. `cartesian_recipe(...)` may build only the
selected `route_family` subrecipe and set the inactive subrecipe to `nothing`;
`:pqs_source_box` route inputs must not require inactive `white_lindsey_*`
fields, and `_cartesian_base_route(kind)` may delete those fields. Explicit
`:white_lindsey_low_order` route support must remain. This ID does not approve
canonical driver changes, numerical kernel changes, terminal lowering or
shellification changes, materialization/artifact schema changes, route-stage
diagnostics, status/report expansion, WL materialization deletion, Cr2 runs, or
new committed tests. `HP-ROUTE-RECIPE-TEST-01` permits only existing direct
route-input tests to be adjusted if necessary.

`HP-ROUTE-INV-FN-01` approves only retained-unit route inventory type cleanup
in `src/pqs_source_box_route_driver_helpers.jl`: remove runtime-keyed
`NamedTuple{unit_keys}` retained-unit inventory shapes and runtime-keyed
`pair_family_counts = NamedTuple{families}(...)`, replacing them with
vector-backed records/tables, stable dictionaries, or helper accessors with
stable concrete types. `HP-ROUTE-INV-TEST-01` approves only `git diff --check`,
package load, H2 base artifact readback, H2 supplemented artifact/readback or
canonical driver path, focused search for absence of those runtime-keyed route
inventories, and no Cr2 run. This lane does not approve raw product source-mode
tuple cleanup, terminal-lowering contract tuple cleanup, retained-unit
transform-contract tuple cleanup, public input `NamedTuple` changes, artifact
sidecar table changes, source files outside the approved file, numerical
kernels, driver changes, report/status/payload expansion, compatibility
adapters, new committed tests, or Cr2 workflow.

`HP-RAW-SRCMODE-FN-01` approves only raw product source-mode inventory cleanup
in `src/cartesian_raw_product_sources/records.jl`,
`src/cartesian_raw_product_sources/source_mode_indices.jl`, and
`src/cartesian_raw_product_sources/summaries.jl`, with narrow consumer wiring
only as needed in
`src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`,
`src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`, and
`src/cartesian_base_hamiltonian.jl`. The target is replacing
`RawProductBoxPlan.source_mode_indices` and `source_mode_column_indices`
variable-length tuple storage with vector-backed storage, or removing the
column storage when accessors can supply the same `1:count` ordering. Fixed
`NTuple{3,Int}` source-mode coordinates remain valid. Accessor compatibility
means preserving deterministic ordered facts, column associations, retained-rule
parity, and manifest source-mode/relation output, not preserving the old
variable-length `Tuple` concrete type. `HP-RAW-SRCMODE-TEST-01` approves only
`git diff --check`, package load, H2 base and supplemented artifact/readback,
H2 R3 endpoint, focused raw-product source order and retained-rule parity,
manifest source-mode/relation inspection, focused search for absence of
tuple-backed `RawProductBoxPlan` source-mode inventories, and no Cr2 run. This
lane does not approve terminal-lowering contract tuple cleanup, broader
retained-unit transform-contract tuple cleanup, broad pair-block/source-box
rewrites, public input `NamedTuple` changes, numerical kernel or route semantic
changes, driver changes, artifact schema changes, report/status/payload
expansion, compatibility adapters preserving the old tuple-backed shape, new
committed tests, or Cr2 workflow.

`HP-CONTRACT-VEC-FN-01` approves only contract-plan vector cleanup in
`src/cartesian_terminal_lowering/contracts.jl`,
`src/cartesian_terminal_lowering/selection.jl`,
`src/cartesian_terminal_lowering/summaries.jl`,
`src/cartesian_retained_unit_transform_contracts/records.jl`,
`src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`, and
`src/cartesian_retained_unit_transform_contracts/summaries.jl`, with narrow
consumer wiring only as needed in `src/pqs_source_box_route_driver_helpers.jl`,
`src/cartesian_base_hamiltonian.jl`, and
`src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`. The
target is replacing `TerminalLoweringPlan.available_contracts`,
`TerminalLoweringPlan.contracts`, and
`RetainedUnitTransformContractPlan.contracts` variable-length tuple storage
with vector-backed storage while preserving `available_contracts(plan)`,
`selected_contracts(plan)`, `contracts(plan)`, `transform_contracts(plan)`,
iteration order, summaries, and existing behavior.
`source_cpbs::Tuple{Vararg{CoordinateProductBox}}` is explicitly out of
scope. `HP-CONTRACT-VEC-TEST-01` approves only
`git diff --check`, package load, H2 base and supplemented artifact/readback,
H2 R3 endpoint, focused terminal-lowering and retained-unit transform-contract
order parity, focused search for absence of targeted tuple-backed plan
inventories, and no Cr2 run. This lane does not approve raw product
source-mode changes, numerical kernel or route semantic changes, driver
changes, artifact/manifest schema changes, report/status/payload expansion,
compatibility adapters preserving old tuple-backed plan field types, new
committed tests, or Cr2 workflow.

`HP-ROUTE-STAGE-TYPE-FN-01` approves only Be2 q5 compile-attributed
route/stage type-surface cleanup in `src/pqs_source_box_route_driver_helpers.jl`
and `src/cartesian_terminal_shellification_geometry.jl`. Approved targets are
`_pqs_source_box_route_driver_terminal_lowering_contract_inventory_from_plan`,
`cartesian_units`,
`_pqs_source_box_route_driver_transform_stage_low_order_summary`,
`cartesian_transforms`,
`_cartesian_terminal_shellification_region_unit_inventory`, and related
terminal-region lowering inventory summary surfaces in
`src/cartesian_terminal_shellification_geometry.jl` only where the same
runtime-sized type-surface pattern appears. Stale compatibility inventories may
be deleted, and remaining runtime-sized `NamedTuple` / `Tuple` carriers may be
replaced with vector-backed compact internal objects, stable dictionaries,
accessors, or smaller summaries, provided H2 base/supplemented artifact
behavior, terminal shellification/lowering order, public driver contract,
artifact/manifest schema, route semantics, and numerical matrices stay
unchanged. `HP-ROUTE-STAGE-TYPE-TEST-01` approves only `git diff --check`,
package load, H2 base and supplemented artifact/readback, H2 R3 endpoint if
terminal realization behavior is touched, focused terminal shellification/
lowering order parity, focused scan for newly introduced runtime-sized
`NamedTuple`/`Tuple` inventories in the approved files, optional Be2 q5
compile/timing comparison after correctness passes, and no Cr2 run. This lane
does not approve driver changes, artifact/manifest changes, public API/export
changes, numerical/raw-block/RG/MWG/IDA changes, route diagnostic/status/report
expansion, committed tests, PackageCompiler/PrecompileTools/sysimage work, or
Cr2 workflow.

`HP-ROUTE-STAGE-CARRIER-FN-01` approves only post-cleanup route/stage carrier
cleanup in `src/pqs_source_box_route_driver_helpers.jl` and
`src/pqs_source_box_diatomic_complete_core_shell.jl`, with
`src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
optional only where directly required to slim terminal realization plan
carriers in the approved path. Approved targets are `cartesian_shells`,
`cartesian_units`, `cartesian_transforms`, terminal topology support-region
planning, terminal retained-rule planning, and directly required terminal
realization plan carriers. Giant shellification, route-skeleton, support-plan,
retained-rule-plan, and terminal-plan `NamedTuple` / tuple shapes may stop
crossing approved stage signatures; necessary facts should move to compact
typed/vector-backed carriers, smaller summaries, accessors, or local
recomputation from canonical objects. Route skeleton construction semantics and
`src/pqs_source_box_route_driver_skeletons.jl` remain out of scope.
`HP-ROUTE-STAGE-CARRIER-TEST-01` approves only `git diff --check`, package
load, H2 base and supplemented artifact/readback, H2 R3 endpoint if terminal
realization is touched, focused terminal support/shellification/lowering order
parity, focused scan for newly introduced runtime-sized `NamedTuple`/`Tuple`
inventories in the approved files, optional Be2 q5 post-cleanup compile/timing
comparison after correctness passes, and no Cr2 run. This lane does not
approve driver changes, artifact/manifest changes, public API/export changes,
numerical/raw-block/RG/MWG/IDA changes, route diagnostic/status/report
expansion, broad route-stage redesign, committed tests,
PackageCompiler/PrecompileTools/sysimage work, or Cr2 workflow.

R3/RG current source authority is compact by design. Read
`docs/src/developer/designs/cartesian_hamiltonian_producer/residual_gaussian_domain_module.md`
for the Residual Gaussian algorithm contract; do not duplicate the full
algorithm in `AGENTS.md`.

Approved R3 compatibility and endpoint surfaces:

- `HP-R3-OBJ-01`, `HP-R3-FN-01`, `HP-R3-FN-02`, and `HP-R3-TEST-01` remain
  approved for the first H2 residual-GTO exact one-body/moment endpoint and its
  standalone gate.
- `HP-R3-FN-03` remains approved for the in-memory residual-MWG/IDA Hamiltonian
  compatibility entry point
  `pqs_terminal_residual_gto_augmented_hamiltonian(...)`.
- `HP-R3-ART-01` remains approved only for the compact supplemented artifact
  provenance writer that adds `supplement_provenance/` to the existing
  Hamiltonian file.
- `HP-HAM-MANIFEST-FN-01` approves only compact JLD2 sidecar groups
  `hamiltonian_manifest/` and `recipe_provenance/` for existing
  `CartesianIDAHamiltonian{Float64}` artifacts. Approved source files are
  `src/cartesian_base_hamiltonian.jl`,
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`, and
  `src/cartesian_ida_hamiltonian.jl` only for a small unexported sidecar
  writer/helper if needed. Existing Hamiltonian matrix keys and
  `read_cartesian_ida_hamiltonian` behavior must not change. The manifest must
  reuse the prior PQS fixed-column/source-mode provenance model:
  `hamiltonian_manifest/final_basis_labels/` has one status-bearing row per
  matrix-order final basis column, and optional
  `final_basis_source_relations/`, `source_shells/`, and `source_modes/`
  groups may be written only for native construction facts. Basis identity is a
  construction label, not a representative center; centers are metadata with
  explicit definition/status. Do not infer shell/ray/radial/source labels from
  centers, nearest-grid snapping, support order, support indices, or
  raw-to-final support.
- `HP-HAM-MANIFEST-TEST-01` approves only existing-reader artifact readback plus
  direct JLD2 sidecar checks for H atom or H2 base artifacts, H2 supplemented
  artifacts, optional practical Be2 supplemented artifacts, explicit
  unavailable/mixed status checks, no inferred-label checks, and no Cr2 run.
  This lane must not add `T_G`, `T_A`, dense transforms, coefficients, raw
  inventories, allocation probes, route reports, status/payload fields, public
  reader APIs, driver public input changes, artifact schema dumps in the
  driver, solver-specific, CR2-consumer-specific, Cr2-specific fields,
  committed Cr2 fixtures, or
  Cr2-specific branches. One-center atom padding is provenance-only in this
  lane; do not change atom parent counts or atom size policy under these IDs.
- `HP-NEST-ART-FN-01` approves only nesting artifact-truth cleanup in
  `src/cartesian_base_hamiltonian.jl`, plus a docstring-only correction in
  `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`.
  Base and recipe provenance must record public `nesting`, and route labels
  must be truthful values derived from `(input.kind, input.nesting)`, including
  `:one_center_pqs_base`, `:one_center_wl_base`, and
  `:z_axis_diatomic_pqs_base`. The original supplemented-WL early-rejection
  boundary in this lane is superseded for the supported z-axis diatomic
  supplemented WL composition cell by `HP-COMP-SUPPWL-*`. This does not
  approve driver public input changes, route skeleton/shellification/terminal-
  lowering changes, raw-block changes, RG/MWG/IDA changes, artifact matrix or
  reader changes, public API/export changes, diagnostics/reports, committed
  tests, or Cr2 workflow.
- `HP-NEST-ART-TEST-01` approves only `git diff --check`, package load, small
  `nesting = :pqs` base artifact/readback with provenance inspection, small
  `nesting = :wl` one-center atom artifact/readback with provenance
  inspection, historical supplemented-WL early rejection before
  `HP-COMP-SUPPWL-*`, and no Cr2 run.
- `HP-HAM-MANIFEST-SRC-FN-01` approves only a compact construction-native
  source-mode provenance seam for optional manifest groups
  `hamiltonian_manifest/source_shells/`,
  `hamiltonian_manifest/source_modes/`, and native
  `final_basis_source_relations/` / `final_basis_labels/` improvements.
  Approved source files are `src/cartesian_terminal_lowering/contracts.jl`,
  `src/cartesian_terminal_lowering/region_contracts.jl`,
  `src/cartesian_raw_product_sources/CartesianRawProductSources.jl`,
  `src/cartesian_raw_product_sources/records.jl`,
  `src/cartesian_raw_product_sources/source_mode_indices.jl`,
  `src/cartesian_retained_units/CartesianRetainedUnits.jl`,
  `src/cartesian_retained_units/records.jl`,
  `src/cartesian_retained_units/lower_contract_units.jl`,
  `src/cartesian_retained_unit_transform_contracts/CartesianRetainedUnitTransformContracts.jl`,
  `src/cartesian_retained_unit_transform_contracts/records.jl`,
  `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`,
  `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`,
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`,
  and `src/cartesian_base_hamiltonian.jl`. Preferred carrier is one internal
  `source_mode_provenance` field on the `cartesian_base_working_basis(...)`
  result; one optional `CartesianTerminalBasisRealization` field is allowed only
  if needed to avoid duplicated or lost terminal construction ordering. The
  seam must not add coefficients, dense transforms, `T_G`, `T_A`, raw
  inventories, route reports, allocation probes, diagnostic payloads, driver
  changes, reader changes, matrix-key changes, public API/export changes, or
  non-native ray/radial labels. Line budget is at most 180 added `src` lines.
- `HP-HAM-MANIFEST-SRC-TEST-01` approves only `git diff --check`, package load,
  H2 base and H2 supplemented artifact write/readback through the existing
  reader, direct JLD2 checks for optional source groups when native rows exist,
  unavailable/mixed status checks for missing labels, optional practical Be2
  manifest inspection, and no Cr2 run. No committed test file is approved.
- `HP-R3U-FILE-01`, `HP-R3U-FN-01`, `HP-R3U-WIRE-01`, and `HP-R3U-TEST-01`
  remain approved only for the non-exported supplemented usability facade and
  its existing standalone H2 validation section.
- The R3 compatibility/artifact owner file remains
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`; it may
  keep small wrappers and artifact/facade hooks, but moved physics helpers
  should delegate to or be deleted in favor of `CartesianResidualGaussians`.

Approved Residual Gaussian module surfaces:

- `HP-RG-FILE-01` approves only
  `src/cartesian_residual_gaussians/CartesianResidualGaussians.jl`,
  `src/cartesian_residual_gaussians/residual_basis.jl`,
  `src/cartesian_residual_gaussians/augmented_operators.jl`, and
  `src/cartesian_residual_gaussians/mwg_interaction.jl`.
- `HP-RG-OBJ-01` approves the residual Gaussian basis object.
- `HP-RG-FN-01` approves `build_residual_gaussian_basis(...)`.
- `HP-RG-FN-02` approves `transform_augmented_operator(...)`.
- `HP-RG-FN-03` approves `moment_matched_gaussians(...)`.
- `HP-RG-FN-04` approves `assemble_residual_ida_interaction(...)`.
- `HP-RG-WIRE-01` approves migration/delegation from the old terminal residual
  file.
- `HP-RG-TEST-01` approves only validation through the existing standalone H2
  residual-GTO/MWG endpoint and optional ignored Be2 measurement.
- `HP-RG-ORTHO-FN-01` approves only robust final residual
  orthogonalization/identity validation in
  `src/cartesian_residual_gaussians/residual_basis.jl`, with
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` allowed
  only for narrow internal keyword plumbing if needed.
- `HP-RG-ORTHO-TEST-01` approves only the existing H2 residual-GTO/MWG endpoint,
  H2 readback if the facade path is touched, ignored strict N2 q5 p10 residual
  audit/artifact smoke, and one passing N2 comparison.
- `HP-RG-IDTOL-FN-01` approves only the default final residual
  `R' S R` identity validation tolerance update to `1.0e-8` in
  `src/cartesian_residual_gaussians/residual_basis.jl`, with
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` allowed
  only for narrow compatibility keyword default plumbing if needed. This
  older default policy is superseded for production by `HP-RG-CUTOFF-FN-01`
  and then `HP-RG-CUTOFF-FN-02`.
- `HP-RG-IDTOL-TEST-01` approves only Be atom cc-pV5Z `lmax = 1`
  residual audit/artifact validation with the same `21` retained residual
  directions, Be atom cc-pVDZ `lmax = 1` comparison, the unchanged H2
  residual-GTO/MWG endpoint, and reporting of `max |G' S R|`,
  `max |R' S R - I|`, allowed tolerance, retained count, minimum retained
  occupation, and final merge condition.
- `HP-RG-CUTOFF-FN-01` supersedes the RG default cutoff/tolerance policy in
  `src/cartesian_residual_gaussians/residual_basis.jl`, with
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` allowed
  only for narrow compatibility keyword default plumbing if needed. The default
  `residual_occupation_cutoff` was `5.0e-8`, and the default final residual
  `R' S R` identity validation `identity_atol` was `5.0e-8`. This older
  production cutoff is superseded by `HP-RG-CUTOFF-FN-02`; the identity
  tolerance remains `5.0e-8`.
- `HP-RG-CUTOFF-TEST-01` approves only Cr atom
  `basis_ns = 9`, `map_ns = 11`, `lmax = 1` residual validation showing the
  marginal `3.637e-8` direction is dropped or the construction passes under the
  new policy, Be atom cc-pV5Z still passing, the unchanged H2 residual-GTO/MWG
  endpoint, and reporting of retained counts, minimum retained occupation,
  `max |G' S R|`, `max |R' S R - I|`, allowed tolerance, and final merge
  condition. It also approves exactly updating the existing H2 endpoint test
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` so both cutoff
  assertions expect `5.0e-8` instead of `1.0e-8`: the in-memory
  `residual.occupation_cutoff` assertion and the artifact/provenance
  `values[:occupation_cutoff]` assertion. No other committed test or fixture
  change is approved.
- `HP-RG-CUTOFF-FN-02` supersedes the production residual occupation cutoff in
  `src/cartesian_residual_gaussians/residual_basis.jl`, with
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` allowed
  only for narrow compatibility keyword default plumbing if needed. The default
  `residual_occupation_cutoff` is `1.0e-6`; the default final residual
  identity validation `identity_atol` remains `5.0e-8`.
- `HP-RG-CUTOFF-TEST-02` approves only residual-only validation after that
  cutoff change: Cr2 owner retained counts should drop from `68 + 68` to
  `62 + 62`; recompute and report residual spectra including `min eig(K_RR)`,
  `min eig(H1_RR)`, and low-mode candidate composition; Be high-zeta and H2
  residual-GTO/MWG endpoints must still pass; existing H2 cutoff/provenance
  assertions may be updated from `5.0e-8` to `1.0e-6`. It does not approve
  full HF, Cr2 artifact/workflow, kinetic/H1 spectral guards, width-filtering
  defaults, or new committed fixtures/tests.
- `HP-RG-SPECTRAL-AUDIT-01` is measurement-only authority after the
  `1.0e-6` cutoff cleanup. Ignored probes may report retained residual counts
  by owner, low `K_RR`, low `H1_RR = K_RR + sum_A Z_A U_A_RR`, low-mode owner
  weights, residual-occupation composition, and one-center atom baselines when
  available. It approves only ignored `tmp/work/*.jl` probes and durable
  text/TSV output under `/Users/srw/dmrgtmp/...` or CR2 run directories. It
  does not approve production source changes, committed tests, artifacts,
  driver changes, MWG/IDA changes, full HF, dense Vee/solver work, automatic
  residual pruning, kinetic/`H1_RR` guards, cutoff/tolerance changes, or
  source instrumentation.
- `HP-RG-INJECT-AUDIT-01` is measurement-only authority for the optional
  injection-plus-RG scheme recorded in
  `docs/src/developer/designs/cartesian_hamiltonian_producer/residual_gaussian_injection_hybrid.md`.
  Ignored probes may classify owner-local principal modes
  `y_i = A_tilde v_i`, sweep trial `residual_injection_cutoff` values,
  globally merge injected modes, report rank/condition of
  `B = G' S Y_inj`, residualize true RG candidates against
  `F = [Y, G Q_perp]`, and report true RG counts, `K_RR`/`H1_RR` spectra, and
  injected-sector one-body projection errors. It approves only ignored
  `tmp/work/*.jl` probes and durable text/TSV output under
  `/Users/srw/dmrgtmp/...` or CR2 run directories. It does not approve
  production source changes, source instrumentation, committed tests,
  artifacts, driver changes, public API/export changes, RG default changes,
  automatic pruning, spectral-guard implementation, MWG/IDA convention
  changes, full HF, dense Vee/solver work, Cr2 full Hamiltonian, Cr2 artifact,
  or Cr2-specific workflow.
- `HP-RG-INJECT-FN-01` approves only the historical default-off direct
  `G`-injection implementation of the injection-plus-RG hybrid in
  `src/cartesian_residual_gaussians/`:
  `residual_basis.jl`, `augmented_operators.jl`, and `mwg_interaction.jl`;
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` is
  allowed only for narrow internal keyword plumbing, same-construction
  validation, and compatibility wiring. The source path may classify
  owner-local principal modes, globally merge injected modes, construct a
  replacement base sector `F`, transform exact one-body operators into
  `[F, R]`, and inherit gausslet-sector IDA for injected directions. It must
  preserve current behavior when `residual_injection_cutoff <= 0`. It does not
  approve default-on injection, driver input, public API/export changes,
  artifact schema/provenance/reader/manifest changes, injection-enabled
  artifact writing, MWG channels for injected directions, global residual
  selection, spectral pruning, full HF, dense Vee/solver work, Cr2 artifact or
  workflow, route/shellification/raw-block changes, or committed tests.
  This remains a preservation-only compatibility contract for live
  default-off code; it is not authority to extend that algorithm or substitute
  it for the occupied-first/protected-main direction.
- `HP-RG-OCC-FIRST-INJECT-FN-01` and
  `HP-RG-OCC-FIRST-INJECT-TEST-01` authorize only the implemented helper
  governed by
  `docs/src/developer/designs/cartesian_hamiltonian_producer/occupied_first_injection.md`.
  Identified `Y_occ` is mandatory; pre-inclusion capture and post-inclusion
  recovery remain distinct; malformed complement/capture geometry fails; weak
  optional directions are rejected and never become MWG residual channels.
  The helper is not wired into the protected builder and is not a direct
  substitute for staged geometry over `M = [G, R_compact]`.
- `HP-RG-PROTECT-ADDREF-FN-01` and
  `HP-RG-PROTECT-ADDREF-TEST-01` authorize the implemented first internal
  protected-localized occupied-reference consumer. Build compact `R` once, make the full-rank
  union of all placed converged packet occupied spaces mandatory in
  `M = [G,R_compact]`, and apply current optional staged selection only after
  that mandatory block. Keep each original packet block separate for
  `P0 = sum_a P_a`; do not globally orthogonalize packet blocks when forming
  the reference density. Packet self-fingerprints and structural owner-local
  mapping remain exact hard failures; the mapped overlap block uses the
  unchanged numerical `1e-10` infinity-norm gate and reports any raw hash
  mismatch only in one nested internal summary. Build placed fitted-potential
  `GG/GA/AA` through the neutral raw-block owner, transform `J0` through the
  existing protected fixed-sector and localized `W` helpers, include
  `E0 = sum_a E_aa + 2*sum_{a<b}E_ab`, and call the existing in-memory
  `ScreenedHartreeCorrection` with native `Vee_L`. Approved files are
  `residual_basis.jl`, `augmented_operators.jl`,
  `cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl`, the two
  `cartesian_reference_density` implementation files, and narrow internal
  composition in `cartesian_protected_ladder_bundle.jl`. One native
  vector-backed residual source-index field is allowed to eliminate duplicate
  compact selection; label parsing is forbidden. Validation is small committed
  algebra coverage plus an ignored physically padded Be2 two-packet gate with
  terminal due diligence. Ordinary fitted-potential consumption reports total
  and self/cross `Tr(P0*J0_fit)-E0_fit` contributions; the decomposition is
  strict, while its magnitude is not a `1e-8 Ha` rejection gate. Retired
  polished packets must be regenerated. No public input, corrected artifact,
  protected atom, counterpoise, compact/high transfer, `Vee` rotation, solver,
  EGOI, exchange, or Cr2 production claim is approved.
- `HP-RG-PROTECT-INJECT-DESIGN-01`,
  `HP-RG-PROTECT-INJECT-FN-01` / `TEST-01`, and
  `HP-RG-PROTECT-ONEBODY-FN-01` / `TEST-01` are governed by
  `docs/src/developer/designs/cartesian_hamiltonian_producer/protected_localized_basis.md`.
  The source-backed internal/default-off convention builds
  `M = [G, R_compact]`, replaces a represented subspace with protected
  originals, transforms exact one-body operators through localized `L`, and
  inherits pre-injection site-order `Vee_M` as `Vee_L`. Gaussian Gram cleanup
  and compact-main representability are separate gates; unsupported broad
  directions are rejected and never become MWG residual channels. The
  completed `HP-RG-PROTECT-VEE-AUDIT-01` invalidated direct `C' V C`; it is
  not an alternative or diagnostic option. This contract does not approve
  public/default workflow, artifact/schema changes, EGOI, ladder behavior,
  screened-reference work, solver/HF, or Cr2 production claims.
- `HP-RG-PROTECT-ART-FN-01` and `HP-RG-PROTECT-ART-TEST-01` govern only the
  implemented opt-in protected-localized Hamiltonian artifact contract in
  `docs/src/developer/designs/cartesian_hamiltonian_producer/protected_localized_artifact.md`.
  Source ownership remains
  `src/cartesian_residual_gaussians/augmented_operators.jl` and
  `src/cartesian_ida_hamiltonian.jl`. The canonical page owns artifact
  identity, native sector/order law, readback compatibility, and rejection
  behavior. This lane does not approve ordinary artifact changes,
  public/default workflow, driver/API/solver behavior, selection changes,
  EGOI, rho0 or screened-Hartree, alternative interaction transforms,
  committed tests by default, or Cr2-specific production behavior.
- `HP-RG-PROTECT-ARTLOC-FN-01` and
  `HP-RG-PROTECT-ARTLOC-TEST-01` govern the implemented native-order locality
  metadata in the same canonical artifact contract. Centers come from actual
  `L`-basis position expectations; deterministic inverse z permutations are
  metadata only and never reorder canonical matrices or native ranges. This
  lane does not approve label-derived centers, new second-moment construction,
  matrix reordering, public/solver workflow, EGOI, additive-reference,
  screening, rho0, committed tests by default, or Cr2 production claims.
- `HP-RG-PROTECT-EGOI-AUDIT-01` is completed historical measurement
  evidence, not active implementation authority.
- `HP-RG-PROTECT-EGOI-FN-01` and `HP-RG-PROTECT-EGOI-TEST-01` are approved
  but pending; the canonical contract is
  `docs/src/developer/designs/cartesian_hamiltonian_producer/retained_gto_egoi.md`.
  No protected retained-GTO helper or focused test is implemented in committed
  source, and uncommitted `src/hamiltonian_corrections.jl` additions are not
  authoritative. The only approved first target is owner-balanced retained
  original `s1+s2`, local symmetric products with `AA-BB` Coulomb acceptance,
  the `M2` mask, and exactly zero disallowed/long-range `DeltaV`. This does not
  approve AB overlap products, `s3`/`p`/`d`, artifacts, public/solver workflow,
  rho0/screened-Hartree changes, or Cr2 production claims.
- `HP-RG-PROTECT-LADDER-XFER-AUDIT-01` is completed historical
  measurement evidence.
- `HP-RG-PROTECT-LADDER-BUNDLE-FN-01` and
  `HP-RG-PROTECT-LADDER-BUNDLE-TEST-01` govern the implemented internal
  opt-in facility in
  `docs/src/developer/designs/cartesian_hamiltonian_producer/protected_localized_ladder.md`.
  Transfers require shared parent lattice, identical supplement and Coulomb
  expansion, and use only `S_BA = <L_B|L_A>` with `C_B = S_BA*C_A`.
  Final self-overlaps are diagnostic only; evaluate transferred densities
  with target `H1_L` and target inherited-site `Vee_L`. Never use generalized
  overlap, source-Hamiltonian/`Vee` transforms, or interaction rotation.
  Restarts are native order. The lane does not approve new representation
  sidecars, solver/UHF continuation, EGOI, screened-Hartree/rho0 changes,
  public defaults, or Cr2 production claims.
- `HP-RG-RHO0-GAL-AUDIT-01`, `HP-RHO0-REFDENS-AUDIT-01`,
  `HP-RHO0-REFDENS-MIXH-AUDIT-01`, `HP-RHO0-FAPP-AUDIT-01`, and
  `HP-RHO0-CORR-AUDIT-01` are completed or superseded historical
  measurement evidence. Their durable interpretation is recorded in
  `docs/src/developer/designs/cartesian_hamiltonian_producer/rho0_reference_density_matrix.md`.
  They authorize no current source work.
- `HP-RHO0-MIXH-GG-FN-01` / `HP-RHO0-MIXH-GG-TEST-01`,
  `HP-RHO0-MIXH-GAAA-FN-01` / `HP-RHO0-MIXH-GAAA-TEST-01`, and
  `HP-RHO0-MIXH-FEXACT-FN-01` / `HP-RHO0-MIXH-FEXACT-TEST-01` remain
  implemented neutral exact-Hartree authority. Their source ownership,
  numerical contract, validation gates, and exclusions are canonical in
  `docs/src/developer/designs/cartesian_hamiltonian_producer/reference_hartree_numerics.md`.
  Correction-policy retirement does not retire these kernels or transforms.
- `HP-RHO0-REFDENS-FN-01` and `HP-RHO0-REFDENS-ERI-01` remain unapproved
  candidates and must not enter the approved source-ID list.
  `HP-RHO0-ANCHOR-FN-01` / `HP-RHO0-ANCHOR-TEST-01` are superseded and
  carry no authority.
- `HP-RHO0-FAPP-FN-01` / `HP-RHO0-FAPP-TEST-01` are implemented but
  caller-free dormant retirement candidates. `HP-RHO0-JANCHOR-FN-01` /
  `HP-RHO0-JANCHOR-TEST-01` are source-backed but superseded in use by the
  screened-Hartree implementation. Neither pair belongs in the approved
  source-ID list, authorizes new callers, or authorizes source work.
- Live `Delta_J0/C` physics is governed by
  `screened_hartree_residual_density.md`; its API by
  `screened_hartree_correction_assembly.md`; and additive molecular
  `P0/J0/E0` by `protected_additive_reference_correction.md`.
  `HP-RHO0-XPAIR-AUDIT-01` remains an approved but deferred H/Be/Be2
  measurement question only. It is not a current blocker or source lane.

Non-negotiable RG guardrails:

- residual directions are selected separately on each physical owner atom and
  merged once;
- residual occupation is not numerical rank, not an integral weight, and not a
  conditioning repair knob;
- exact augmented one-body/moment transformation is not the MWG approximation;
- MWG descriptors are not invariant under arbitrary residual rotations and must
  be computed from the final merged residual basis;
- `V_GM` uses weight-aware final-basis density normalization for PQS shell
  blocks;
- final residual identity validation may use the approved
  `HP-RG-ORTHO-FN-01` combined absolute/relative check only for small
  floating-point overshoots after owner-local selection and a healthy final
  merge; it must not change occupation cutoff, selection semantics, or merge
  failure rules;
- `HP-RG-IDTOL-FN-01` sets the production default final residual identity
  tolerance to `1.0e-8` in the older Be tolerance lane. This is now superseded
  by `HP-RG-CUTOFF-FN-01` and `HP-RG-CUTOFF-FN-02` for production defaults;
- `HP-RG-CUTOFF-FN-02` supersedes the production residual occupation cutoff:
  `residual_occupation_cutoff = 1.0e-6`, while
  `identity_atol = 5.0e-8` remains unchanged. This is an explicit owner-local
  residual selection policy for marginal Cr2 directions; it does not change
  owner grouping, merge checks,
  `G' S R` validation, width/zeta filtering, MWG/IDA, artifacts, driver
  workflow, public API, or source files outside the approved RG owner/plumbing
  surface;
- `HP-RG-SPECTRAL-AUDIT-01` is a follow-up residual-sector measurement lane,
  not a guard implementation. It may classify the remaining low Cr2
  residual-only `H1_RR` mode after the `1.0e-6` cutoff, but it must not change
  selection policy, add pruning, or add production instrumentation;
- RG does not own basis loading, parent lattice construction, terminal topology,
  raw analytic formula ownership, facade parsing, artifact writing,
  `supplement_provenance/`, report/status/payload objects, or public exports.

The active H2 owner-local residual-MWG endpoint has augmented dimension `489`
and lowest-orbital IDA self-Coulomb `0.4574265214362075` within `1.0e-10`.
Older R3-B scalars and global-selection construction paths are historical only.
Do not add width scaling, tolerance relaxation, global raw-candidate Lowdin,
global raw-column pivoted-Cholesky residual selection, public export,
Cr2-specific facade support, full Cr2 Hamiltonian/artifact, new committed
tests, driver/bin/tool workflow, artifact schema expansion, or further
unapproved tolerance changes without a prior docs-only amendment.

Approved neutral Cartesian Gaussian raw-block nuclear owner:

- `HP-CGRB-FILE-01` approves only
  `src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl`,
  `src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl`, and the internal
  `src/GaussletBases.jl` include needed to load the module. No public export is
  approved.
- `HP-CGRB-FN-01` approves only exact uncharged by-center Cartesian Gaussian
  nuclear parent-supplement `G-A` and supplement-supplement `A-A` raw-block
  construction, including analytic 1D nuclear factors, unique coordinate reuse,
  upper-triangular `A-A` assembly/mirroring, function-local scratch reuse, and
  term-first contraction.
- `HP-CGRB-FN-02` approves only
  `src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl` for reorganizing the
  neutral nuclear kernel around unique one-dimensional supplement axis-family
  reuse. It must not be used for non-nuclear overlap/kinetic/moment work.
- `HP-CGAI-FN-01` is optional helper authority only for
  `src/cartesian_gaussian_axis_integrals.jl` support needed by
  `HP-CGRB-FN-02`; it is not a broad raw-block or non-nuclear authority.
- `HP-CGRB-WIRE-01` approves only behavior-preserving rewiring of the Residual
  Gaussian and Qiu-White nuclear callers in
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`,
  `src/ordinary_qw_raw_blocks.jl`, and `src/ordinary_qw_operator_assembly.jl`,
  with duplicate route-local nuclear loops deleted after parity.
- `HP-CGRB-TEST-01` approves the existing H2 Residual Gaussian endpoint,
  ignored Be2 Residual Gaussian parity/performance if needed, ignored Cr2 q4
  exact nuclear block parity, and one small standalone Qiu-White nuclear parity
  fixture at `test/nested/cartesian_gaussian_raw_blocks_nuclear_runtests.jl` if
  no existing test can host it cleanly. Do not add it to `test/runtests.jl`
  without a later amendment.

Approved neutral Cartesian Gaussian raw-block non-nuclear owner:

- `HP-CGRB-NN-FILE-01` approves only
  `src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl` and the include in
  `src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl` needed to
  load it. Root include changes are not approved unless a later amendment names
  a real include-order blocker.
- `HP-CGRB-NN-FN-01` approves only exact non-nuclear Cartesian Gaussian
  parent-supplement `G-A` and supplement-supplement `A-A` raw-block
  construction for overlap, kinetic, coordinate moments `x`/`y`/`z`, and
  second moments `x^2`/`y^2`/`z^2`. It may use analytic 1D tables, unique
  supplement axis-family reuse, canonical `A-A` family-pair keys, orientation
  handling, upper-triangular `A-A` assembly/mirroring, function-local scratch
  reuse, and coupled product-axis contraction.
- `HP-CGRB-NN-WIRE-01` approves only behavior-preserving rewiring of Residual
  Gaussian exact-operator/mixed-overlap setup and Qiu-White non-nuclear callers
  in `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`,
  `src/ordinary_qw_raw_blocks.jl`, and
  `src/ordinary_qw_operator_assembly.jl`, with duplicate route-local
  non-nuclear loops deleted after parity.
- The main diatomic Qiu-White non-nuclear path has crossed this lane. Remaining
  QW-local non-nuclear cross/self helpers used by atomic QW reference,
  factor-term, hybrid sidecar, dense-parent probe, or CPB/provider surfaces are
  not dead duplicates under `HP-CGRB-NN-WIRE-01`. They remain retained
  reference/sidecar/provider surfaces until a later amendment either adds
  neutral factor-block ownership or explicitly approves rewiring those
  callers.
- `HP-CGRB-NN-TEST-01` approves the existing H2 Residual Gaussian endpoint,
  ignored Be2 Residual Gaussian parity/performance if needed, ignored Cr2 q4
  non-nuclear raw-block parity, residual mixed-overlap parity, and one small
  standalone Qiu-White non-nuclear parity fixture at
  `test/nested/cartesian_gaussian_raw_blocks_non_nuclear_runtests.jl` if no
  existing test can host it cleanly. Do not add it to `test/runtests.jl`
  without a later amendment.

The neutral raw-block owner must not construct pair factors, MWG interaction,
terminal projection, Residual Gaussian selection/transforms, Qiu-White route
objects, parent construction, final-basis `G-G` product-matrix optimization,
persistent caches, metadata/report/status/payload fields, artifacts, public
API, Cr2 facade support, or Cr2 artifact workflow. Physical nuclear charges
are applied by consumers, not by the neutral uncharged nuclear kernel.

Approved R3/RG terminal `G-G` product-matrix optimization:

- `HP-R3GG-FN-01` approves only
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` and, if
  needed for small internal terminal-product workspace/helper reuse,
  `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`.
- Approved product matrices are only the final-basis `G-G` kinetic,
  coordinate-moment, and second-moment matrices used by
  `pqs_terminal_residual_gto_augmented_operators(...)`.
- Allowed shapes are accumulation of the three kinetic-axis contributions into
  one destination, reuse of an already constructed same-construction base
  Hamiltonian kinetic block when available and validated equal, axis-by-axis
  product/transform for coordinate and second moments, function-local scratch
  reuse across consecutive product assemblies, and deletion/simplification of
  `_r3a_product_matrix(...)` when no live caller remains.
- `HP-R3GG-TEST-01` approves only the existing H2 Residual Gaussian endpoint,
  ignored Be2 measurement, ignored Cr2 q4 `K_GG`/moment `G-G` parity, exact
  operator finiteness/symmetry, existing base `G-G` block equality checks, and
  Cr2 q4 exact-operator allocation remeasurement. No new committed test file
  is approved.

This G-G lane must not change `G-A`/`A-A` raw blocks, nuclear raw blocks,
unit-nuclear Gaussian-sum construction, IDA/MWG, residual selection,
orientation, or transforms, terminal basis realization, Qiu-White semantics,
route setup, parent construction, persistent caches, metadata/report/status/
payload fields, artifacts, public API/export, Cr2 facade support, or Cr2
artifact workflow. Line budget is at most 100 added `src` lines.

Measurement-only current R3 allocation decision:

- `HP-R3REM-AUDIT-01` approves only ignored `tmp/work` measurement/probe work
  to classify remaining Cr2 q4 exact augmented-operator allocation after
  `954c86cd`.
- It is intentionally not listed as production source authority.

Approved R3 unit-nuclear `U_GG` Gaussian-sum optimization:

- `HP-R3UN-FN-01` approves only terminal final-basis unit-nuclear `U_GG`
  Gaussian-sum allocation reduction in
  `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`, with
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` allowed
  only for narrow caller wiring if needed.
- Target functions are `_accumulate_terminal_gaussian_sum!` and
  `_terminal_gaussian_sum_action`.
- Allowed changes are function-local scratch/workspace reuse, in-place
  accumulation into caller destinations, allocation reduction in factor lookup
  and terminal Gaussian-sum action, and deletion/simplification of obsolete
  allocation-heavy code in that path.
- `HP-R3UN-TEST-01` approves only H2 endpoint validation, Be2
  facade/readback measurement, Cr2 exact-operator audit with before/after
  `U_GG` allocation, Cr2 `U_GG` replay parity, and finite/symmetric exact
  operators. No committed test file is approved.
- This lane must not change neutral raw blocks, terminal kinetic/moment `G-G`
  products, residual Gaussian algorithms or transforms, IDA/MWG, Qiu-White
  semantics, route/stage setup, raw-block setup, parent construction, terminal
  basis realization, persistent caches/workspaces, broad Gaussian-sum
  frameworks, metadata/report/status/payload fields, artifacts, public
  API/export, Cr2 facade support, or Cr2 artifact workflow. Line budget is at
  most 100 added `src` lines.

Approved R3 same-construction base K/U reuse:

- `HP-R3BASE-FN-01` approves only reuse of already-built same-construction base
  final-basis kinetic `K_GG` and unit nuclear `U_GG[A]` blocks in supplemented
  exact augmented operators.
- Approved files are
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` and
  `src/cartesian_base_hamiltonian.jl`.
- Allowed wiring is passing `base_ham.kinetic` and
  `base_ham.nuclear_attraction_unit_by_center` from the same
  `cartesian_base_working_basis(...)` construction path into the augmented
  product/unit-nuclear construction, with matrix-dimension and center-count
  validation. Existing recomputation behavior must remain when trusted base
  blocks are not supplied.
- `HP-R3BASE-TEST-01` approves only H2 R3 endpoint validation, Be2
  facade/readback validation, ignored Cr2 exact-operator attribution or replay
  showing base K/U reuse parity and allocation effect, and finite/symmetric
  exact-operator checks. No committed test file is approved.
- This lane must not change public API/export, canonical driver, raw blocks,
  residual selection/orientation/transforms, MWG/IDA conventions, terminal
  product or Gaussian-sum kernels, route/stage setup, metadata/status/report/
  artifact schema, persistent cache/workspace objects, committed tests, Cr2
  workflow, or files outside the approved surfaces. Line budget target is under
  100 added `src` lines.

Approved canonical-driver call-site wiring for R3 same-construction K/U reuse:

- `HP-R3BASE-DRV-WIRE-01` approves only `bin/cartesian_ham_builder.jl`.
- In supplemented mode only, the driver may pass `base_ham.kinetic` as
  `base_kinetic` to `cartesian_residual_gto_augmented_products(...)` and
  `base_ham.nuclear_attraction_unit_by_center` as `base_unit_nuclear` to
  `cartesian_residual_gto_augmented_unit_nuclear(...)`.
- Public inputs, hooks, timing labels, visible stage sequence, artifact schema,
  and the driver contract must remain unchanged.
- `HP-R3BASE-DRV-TEST-01` approves only `git diff --check`, package load, H2
  supplemented driver artifact/readback, optional practical Be2 supplemented
  driver/readback, and no Cr2 run.
- This lane must not change source/kernels, diagnostics, hooks, timing labels,
  public inputs, artifacts, committed tests/fixtures, Cr2 workflow, or any file
  outside `bin/cartesian_ham_builder.jl`.

`HP-FN-03` specifically approves
`src/cartesian_final_basis_realization/pqs_terminal_one_body.jl` as the Slice B
source file. It does not approve a new K/U payload, stage-return field, report
object, persistent one-body orchestration API, or status vocabulary.

PQS terminal basis blocks must remain supported on owned terminal regions.
Previous-block projection, recursive projection, and effective-support growth
onto earlier terminal supports are not approved. Cross-block overlap is zero by
construction because parent rows are orthonormal and terminal support regions
are disjoint. A nonzero structural overlap means duplicated support rows,
incorrect row restriction, wrong support ownership, or an indexing error; it is
not a physical residual or repair path. Cross-block kinetic, nuclear, and IDA
operator terms may still be nonzero.

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

## Canonical Cartesian driver

`bin/cartesian_ham_builder.jl` is now the canonical compact, human-facing
Cartesian Hamiltonian producer driver. Its job is to prove the approved
producer paths work together by producing a Hamiltonian artifact directly.

Approved under `HP-DRV-FILE-01`, `HP-DRV-FN-01`, `HP-DRV-STAGE-FN-01`,
`HP-DRV-STAGE-WIRE-01`, `HP-DRV-STAGE-TEST-01`, `HP-DRV-INV-FN-01`,
`HP-DRV-INV-TEST-01`, and `HP-DRV-TEST-01`:

- visible editable defaults near the top of the file;
- optional trusted local Julia input file for project-specific defaults;
- command-line `key=value` overrides;
- visible public `system`, `basis`, and optional `supplement` contract
  construction before calling an approved facade;
- compact normalized run summary;
- visible physics-level construction stages and coarse user-facing phase
  timing;
- base or supported supplemented Hamiltonian construction through approved
  producer surfaces;
- artifact write and optional readback check.

The intended shape is compact and copyable. Consumers may copy the standard
driver for project-specific customization; copied local drivers are not
canonical repo surfaces.

The only approved compact run-level hooks are `check_file`, `print_contract`,
`print_timing`, and `expected_dimension`. They are for human expert review and
Codex-controlled artifact checks only. The driver may use `basisname = nothing`
as the base-mode selector; `basisname !== nothing` selects a supported
supplemented mode and is further governed by the composition IDs. The original
driver-stage lane covered supplemented diatomics only; `HP-COMP-SUPPATOM-*`
separately approves relaxing the old `Natom == 1` rejection. `padding` is
physical box padding around the atom or two nuclei and maps internally to the
existing facade fields. The driver must not expose private route-stage
choreography as a substitute for constructing public `system`, `basis`, and
optional `supplement` objects.

`HP-DRV-INV-FN-01` and `HP-DRV-INV-TEST-01` govern the implemented compact
terminal-region / shellification inventory summary in the canonical driver
output. Approved files are `bin/cartesian_ham_builder.jl` and
`src/cartesian_base_hamiltonian.jl`, with optional compact accessors in
`src/pqs_source_box_route_driver_helpers.jl`,
`src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`, and
`src/cartesian_final_basis_realization/terminal_face_product_blocks.jl` only if
directly required. The summary may print bounded rows with region label/index,
region kind, lowering or realization kind, support rows, final columns,
compression ratio, shell index or explicit unavailable status, index ranges
for `x`/`y`/`z`, physical coordinate ranges for `x`/`y`/`z`,
identity-vs-compact/product realization, and slab axis/side/thickness/stack
facts when applicable. Physical `x`/`y` ranges are required, not only `z`,
because angular-balance review compares transverse physical scale against the
bond-axis margin. It should include total base final dimension, supplemented
final dimension when applicable, and visible direct identity slab sectors if
any remain. This is human-facing driver output, not an artifact schema, route
diagnostic, public input, stop-after control, solver hook, broad status/report
payload, source-mode dump, pair inventory, raw-block dump, all-row listing, or
full metadata dump. It must not change numerical construction, shellification,
terminal lowering, retained units, transform contracts, terminal realization,
RG/MWG/IDA, Hamiltonian assembly, artifacts/readers, public exports, Cr2
workflow, stage sequence, or driver inputs. The implementation line-budget
target was `80` added `src`/`bin` lines.

`HP-DRV-SHELLDD-FN-01` and `HP-DRV-SHELLDD-TEST-01` govern the implemented
standard terminal due-diligence report for Cartesian/PQS terminal bases,
as recorded in
`docs/src/developer/designs/cartesian_hamiltonian_producer/terminal_shellification_due_diligence.md`.
The implementation extends
`src/cartesian_base_hamiltonian.jl`'s `_cartesian_terminal_inventory_rows(...)`
and joining existing terminal inventory rows with terminal retained-rule
plan/support records. `bin/cartesian_ham_builder.jl` prints the bounded
report through the canonical driver summary path. `src/pqs_source_box_route_driver_helpers.jl`
is optional only if a compact accessor is directly required. The report must
include normalized system/geometry facts, validated atom locations, bond axis/
length, snapped nuclear indices, parent physical bounds, parent axis counts,
1D center summaries or bounded tables, spacing summaries, gausslet/IDA weight
statistics, dimension/compression accounting, shell-by-shell order/key, role,
region kind, shell index, owner/contact/shared classification, index and
physical boxes, physical side lengths/aspect ratios, actual and expected
aspect-balanced source-mode shape, source-mode count, retained count, final
column range, lowering/retained/realization rule, slab metadata, and advisory
warning flags. Consumers are expected to inspect this report before
interpreting energies, residual
behavior, or injection behavior. Warning flags are advisory by default. This
lane does not approve artifact schema/provenance/reader changes, public input
or semantic changes, shellification policy changes, source-mode selection
changes, aspect-balanced source-mode implementation, numerical construction
changes, dense coefficient/transform/support dumps, source-mode/pair/
raw-block/all-row/full-metadata dumps, broad report payloads, Cr2 workflow, or
committed fixtures/tests by default. Gausslet/IDA weight summaries are
diagnostic only and are not residual integral weights, MWG weights, or proof of
quadrature quality. Later implementation line budget is target `180` added
`src`/`bin` lines.

`HP-PQS-ASPECTSHELL-FN-01` and `HP-PQS-ASPECTSHELL-TEST-01` approve the
separate future source-policy lane for z-axis diatomic PQS complete-shell
source modes, as recorded in
`docs/src/developer/designs/cartesian_hamiltonian_producer/pqs_complete_shell_aspect_source_modes.md`.
This lane is not part of due-diligence reporting: `HP-DRV-SHELLDD-*` may report
that a physically rectangular complete shell is represented by cubic
`(q,q,q)` source modes, but only `HP-PQS-ASPECTSHELL-*` may change the actual
basis construction to aspect-aware `(q,q,L)` source modes. The old
angular-resolution code to recover explicitly is in
`src/cartesian_nested_diatomic.jl` and
`src/cartesian_nested_faces.jl`, especially the
`_nested_diatomic_*reference_band`, adaptive retained-count, source-dimension
plan, and `_nested_projected_q_shell_layer(...)` helpers. The approved seam is
in `src/pqs_source_box_route_driver_helpers.jl`, after shellification has
produced complete-shell regions and parent/bundle facts but before
lowering-contract inventory, retained-unit plans, retained-unit transform
contracts, and terminal retained-rule plans are frozen. Additional approved
files are `src/cartesian_terminal_lowering/region_contracts.jl`,
`src/pqs_multilayer_shell_source_plan.jl`, and
`src/pqs_multilayer_shell_region_plan.jl`; the old nested helper files are
optional only if directly needed, and
`src/pqs_source_box_diatomic_complete_core_shell.jl` is optional only for
support-record consistency. `region_contracts.jl` is too early to choose `L`
by itself, while `pqs_multilayer_shell_source_plan.jl` is too late to be the
only fix. Later implementation may change retained counts, final dimensions,
Hamiltonian matrices, and energies, so old scalar targets tied to cubic
complete-shell source modes must be remeasured. This lane does not approve
artifact schema/provenance/reader changes, public input or driver semantic
changes, WL source-mode policy changes, thin-slab/angular
z-extension/direct-core/RG/MWG/IDA/global-injection changes, old route-global
materialization revival, broad source-mode/report/payload frameworks, Cr2
production claims, committed fixtures/tests by default, or more than target
`160` added `src` lines without a new amendment.

`HP-DRV-NEST-FN-01` and `HP-DRV-NEST-WIRE-01` approve one visible construction
family input, `nesting = :pqs` or `nesting = :wl`, in
`bin/cartesian_ham_builder.jl` plus narrow input plumbing in
`src/cartesian_base_hamiltonian.jl`. `:pqs` maps to the existing
`:pqs_source_box` route family and remains the default. `:wl` maps only to the
existing `:white_lindsey_low_order` route family. This is not a diagnostic
route switch: do not expose route skeletons, retained rules, raw-block
switches, stop-after controls, diagnostics, route reports, or route-stage
labels. Supplemented `nesting = :wl` is governed by `HP-COMP-SUPPWL-*` for the
supported homonuclear z-axis diatomic composition cell; unsupported geometry or
supplement combinations must still reject clearly.
`HP-DRV-NEST-TEST-01` approves default `:pqs` validation plus one small base
artifact/readback path with `nesting = :wl`, and no Cr2 run.

The target producer shape is the 2 x 2 x 2 composition of geometry
(`atom` or z-axis diatomic), `nesting` (`:pqs` or `:wl`), and supplement state
(`off` or `on`) recorded in
`docs/src/developer/designs/cartesian_hamiltonian_producer/nesting_supplement_composition_plan.md`.
This is planning authority except where a composition cell is explicitly
promoted below. Supplemented atoms and supplemented White-Lindsey are approved
only where `HP-COMP-SUPPATOM-*` and `HP-COMP-SUPPWL-*` name the exact cells;
do not implement any remaining missing cells as driver-level special cases or
parallel Hamiltonian builders.

`HP-COMP-WLDIAT-FN-01` and `HP-COMP-WLDIAT-TEST-01` promote the first
composition cell: `Natom = 2`, `nesting = :wl`, `basisname = nothing`.
Approved source files are
`src/pqs_source_box_diatomic_complete_core_shell.jl`,
`src/cartesian_terminal_shellification_geometry.jl`,
`src/cartesian_terminal_lowering/selection.jl`,
`src/cartesian_terminal_lowering/region_contracts.jl`,
`src/pqs_source_box_route_driver_helpers.jl`,
`src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`,
`src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`,
`src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`, and
`src/cartesian_base_hamiltonian.jl`. The pass may produce native WL z-axis
diatomic terminal records, route them through the existing
`CartesianTerminalBasisRealization` and staged base Hamiltonian path, and use
`:z_axis_diatomic_wl_base` as a truthful existing-schema route provenance
value. It must not add driver special cases, revive/adapt old WL H1/H1+J
materialization, change artifact schema/matrix keys/reader behavior, touch
RG/MWG/supplement work, add route diagnostics/status/report payloads, create a
parallel Hamiltonian builder, add committed tests, or run Cr2. Line budget is
at most 250 added `src` lines, with blocker-only guard cleanup expected where
practical.

`HP-WLDIAT-COMPACT-FN-01` and `HP-WLDIAT-COMPACT-TEST-01` approve the
narrow White-Lindsey z-axis diatomic compact retained-basis correction. The
existing WL diatomic artifact path may be mechanically reachable, but the
current elongated shared-shell boundary-stratum identity realization is not the
intended compact WL retained basis and is not a fair production comparison
against PQS. Approved files are
`src/cartesian_shellification/terminal_geometry.jl`,
`src/cartesian_terminal_lowering/region_contracts.jl`,
`src/cartesian_retained_units/lower_contract_units.jl`,
`src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`,
`src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`,
and narrow route wiring in `src/pqs_source_box_route_driver_helpers.jl` only if
needed. The pass must preserve the WL unit-based model of faces, edges,
corners, and small boundary units after shellification, but each WL unit must
carry or realize the compact retained basis generated by products of
one-dimensional contractions rather than retaining full-support identity rows.
Identity realization remains valid only for true direct/core identity units,
not for WL boundary-stratum retained units. Historical deleted WL coefficient
helpers may be used only as donor/reference material for the compact
CPB-local product-of-1D coefficient primitive; do not revive the old
route-global WL stack, reports, adapters, or H1/H1+J materialization. The pass
must not force a persistent shell object after splitting, fake compactness by
dropping rows or relabeling full-support units, change the driver, artifact
schema/provenance, PQS behavior, Hamiltonian assembly, raw blocks, RG/MWG/IDA,
diagnostics/status/report payloads, committed tests/fixtures, or Cr2 workflow.
The same public `ns` is the fair starting input for PQS/WL comparison after
this correction, but it is not a promise of identical dimensions.

`HP-WLDIAT-PARITY-FN-01` and `HP-WLDIAT-PARITY-TEST-01` approve the narrow
White-Lindsey boundary-stratum retained-count parity cleanup in
`src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`.
Direct nucleus-centered core blocks keep the odd-side requirement, but
boundary shells/strata and non-direct support regions must not inherit that
symmetric-odd rule. For public `nesting = :wl`, `ns = 4`, route-local `q = 2`
must retain the boundary shell count `4^3 - 2^3 = 56`, not the current
inherited-donor result of `26`; `ns = 5` should retain `5^3 - 3^3 = 98`. This
lane may only set the boundary-stratum product contraction to use the requested
retained count without symmetric-odd coercion. It must not change core/direct
identity rules, public `ns` normalization, route skeletons, shellification,
retained-unit metadata shape, driver behavior, artifacts/provenance, PQS
behavior, Hamiltonian assembly, raw blocks, RG/MWG/IDA, old WL materialization,
diagnostics/status/report payloads, committed tests/fixtures, or Cr2 workflow.

`HP-COMP-BASEDIAT-FN-01` and `HP-COMP-BASEDIAT-TEST-01` approve only
`src/cartesian_base_hamiltonian.jl` for relaxing the base two-center branch
from H2-only to explicit homonuclear z-axis all-electron diatomics. The input
contract must require equal symbols, equal finite positive integer-valued
nuclear charges, two finite distinct z-axis centers, and neutral
`nup + ndn == sum(charges)`. Symbols are labels only; charges and explicit
electron counts are authority. The basis contract stays unchanged and both
`nesting = :pqs` and `nesting = :wl` remain visible construction choices. This
lane does not approve driver changes, route skeleton/shellification/terminal
lowering changes, raw-block changes, supplement/RG/MWG changes, artifact
schema or reader changes, public API/export changes, element lookup/default
tables, heteronuclear or translated/general geometry support, committed tests,
or Cr2 workflow. Target line budget is under 60 added `src` lines.

`HP-COMP-SUPPWL-FN-01` and `HP-COMP-SUPPWL-TEST-01` approve only the
supplemented White-Lindsey z-axis diatomic composition lane. Approved source
surface is `src/cartesian_base_hamiltonian.jl`, with
`src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` allowed
only if a direct genericity blocker appears in the existing RG/MWG compatibility
entry point. The lane may remove the early supplemented-`nesting = :wl`
blockers only if the existing Residual Gaussian/MWG path works with the WL
`CartesianTerminalBasisRealization`. It must preserve the supplement contract,
residual selection, exact augmented operators, residual MWG/IDA interaction,
base K/U reuse, artifact keys, manifest/provenance, driver inputs, and stage
labels. It does not approve driver changes, supplemented atoms, route skeleton
or shellification changes, terminal lowering changes, raw-block changes,
residual-selection changes, MWG/IDA convention changes, artifact schema or
reader changes, public API/export changes, old WL H1/H1+J materialization,
solver/ECP work, diagnostics/status/report payloads, committed tests, or Cr2
workflow. Target line budget is under 80 added `src` lines.

`HP-COMP-SUPPATOM-FN-01` and `HP-COMP-SUPPATOM-TEST-01` approve only the
supplemented one-center atom composition lane. Approved implementation surfaces
are `src/cartesian_base_hamiltonian.jl` and `bin/cartesian_ham_builder.jl`,
with `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
allowed only for a direct one-owner RG/MWG genericity blocker. The lane may
allow origin-centered one-center all-electron atoms with `basisname !==
nothing` for both `nesting = :pqs` and `nesting = :wl`, using the existing base
atom validation, terminal basis construction, residual Gaussian augmentation,
exact augmented operators, MWG/IDA interaction, base K/U reuse, assembly,
writer, readback, manifest, and provenance. It may select
`legacy_atomic_gaussian_supplement(...)` for one-center inputs and keep the
existing diatomic supplement loader for two-center inputs. It does not approve
a separate atom-only Hamiltonian builder, new driver inputs, route switches,
diagnostics, stop-after controls, new stage labels, route/shellification/
terminal-lowering changes, raw-block changes, residual-selection changes,
MWG/IDA convention changes, artifact schema or reader changes, public
API/export changes, solver/ECP work, status/report payloads, heteronuclear or
general geometry, translated atoms, committed tests, or Cr2 workflow. Target
line budget is under 80 added `src`/`bin` lines.

`HP-COMP-ATOMBOX-FN-01` and `HP-COMP-ATOMBOX-TEST-01` approve only the
one-center atom parent-sizing correction in `src/cartesian_base_hamiltonian.jl`.
The fixed `2*q + 1` atom parent-axis count artifact must be removed; public
`basis.radius` is the one-center physical box extent authority, and parent
axis counts must depend on radius plus `core_spacing`/existing spacing policy
analogously to the z-axis diatomic physical-extent sizing. This
lane preserves origin-centered atom validation, explicit charge/electron-count
validation, `nesting = :pqs` and `nesting = :wl`, supplemented atoms, artifact
keys, manifest/provenance, and canonical driver inputs. It does not approve
driver changes, route-family switches, raw-block changes, residual-selection
changes, MWG/IDA convention changes, artifact schema or reader changes, public
API/export changes, solver/ECP work, diagnostics/status/report payloads,
committed tests, Cr2-specific workflow, translated atoms, non-origin atom
support, element lookup/default tables, broad parent-construction rewrites, or
diatomic sizing changes. Target line budget is under 80 added `src` lines.

`HP-COMP-NS-FN-01` and `HP-COMP-NS-TEST-01` approve only public size-parameter
normalization in `bin/cartesian_ham_builder.jl` and
`src/cartesian_base_hamiltonian.jl`. The durable public field is `ns`, the
requested cube/source/nesting size. Route-local `q` is derived after selecting
`nesting`: `q = ns` for `nesting = :pqs`, and `q = ns - 2` for
`nesting = :wl`. `nesting = :wl` must reject `ns < 3`. Legacy public `q` may
remain temporarily only as compatibility: if `ns` is absent, derive `ns` from
`q` and `nesting`; if both are present, require consistency with the selected
nesting or throw `ArgumentError`. Driver examples and new docs should use
`ns`, not `q`. Provenance may add compact `ns`, `q_rule`, and `ns_source`
entries next to the existing derived `q` in `producer_provenance/` and
`recipe_provenance/`. This lane must not change matrix keys, reader behavior,
artifact format, route skeletons, shellification, terminal lowering, raw
blocks, RG/MWG/IDA, solver/ECP workflow, Cr2 workflow, route diagnostics,
status/report payloads, driver hooks/stage labels, or committed tests.

`HP-COMP-NSCORE-FN-01` and `HP-COMP-NSCORE-TEST-01` approve only direct
nucleus-centered core side parity cleanup in
`src/pqs_source_box_route_driver_helpers.jl`, with
`src/cartesian_base_hamiltonian.jl` allowed only if needed for one-center parent
minimum sizing consistency. Route-local `q` derivation remains `q = ns` for
PQS and `q = ns - 2` for WL, but direct core side must come from public `ns`:
`direct_core_side = isodd(ns) ? ns : ns + 1`. This oddization rule is only for
direct nucleus-centered core identity blocks. It must not apply to boundary
shells, WL boundary-stratum retained products, or non-direct support regions.
WL boundary retained counts remain `ns = 4 -> 56`, `ns = 5 -> 98`, and
`ns = 6 -> 152`. This lane does not approve driver changes, public input
changes, route skeleton redesign, terminal lowering, retained-unit or terminal
realizer changes, artifact schema changes, manifest expansion, old WL
materialization revival, committed tests/fixtures, or Cr2 workflow.

`HP-COMP-SHELLGEOM-FN-01` and `HP-COMP-SHELLGEOM-TEST-01` approve only common
terminal shell decomposition audit/cleanup in
`src/cartesian_shellification/terminal_geometry.jl` and narrow caller plumbing
in `src/pqs_source_box_route_driver_helpers.jl`. Direct core regions, terminal
shell regions, owned support rows, ordering, and coverage are common geometry
and must be computed before the PQS/WL route-family split. PQS full source-box
geometry and WL face/edge/corner/stratum geometry are retained-construction
geometry after common shell records exist. This lane does not approve driver
changes, public input changes, route skeleton redesign, terminal lowering
redesign, retained-unit or transform changes, PQS retained-mode realization
changes, WL boundary coefficient changes, artifact/manifest/reader changes,
Hamiltonian/IDA/MWG/RG/raw-block changes, old WL materialization, committed
tests/fixtures, or Cr2 workflow.

`HP-COMP-SHELLGEOM-DIAT-FN-01` and
`HP-COMP-SHELLGEOM-DIAT-TEST-01` extend common shell decomposition authority to
z-axis diatomic shellifier entry. For a fixed public z-axis diatomic system,
parent axes, public `ns`, direct core side, nuclear centers, and bond axis,
PQS and WL must call the same common shellifier with the same first-step
arguments. Central-gap/contact, shared-shell, and outer-mismatch ownership are
common shell geometry. PQS `q` and WL inner side are retained-construction
inputs after common shell records exist. This lane allows only caller plumbing
and shellifier-boundary naming/input cleanup in
`src/cartesian_shellification/terminal_geometry.jl` and
`src/pqs_source_box_route_driver_helpers.jl`; it does not approve changing the
central-gap/contact algorithm, terminal lowering, retained units, PQS retained
realization, WL boundary coefficients, route skeletons, artifacts, driver
inputs, committed tests/fixtures, or Cr2 workflow.

`HP-COMP-THINSLAB-FN-01` and `HP-COMP-THINSLAB-TEST-01` supersede the
outer-mismatch-only `HP-COMP-OUTERMM-*` lane. They approve one unified
thin-slab stack compact-lowering repair for z-axis diatomics. Approved source
files are `src/cartesian_terminal_lowering/selection.jl`,
`src/cartesian_terminal_lowering/region_contracts.jl`,
`src/cartesian_retained_units/lower_contract_units.jl`,
`src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`,
`src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`, and
`src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`,
with `src/cartesian_shellification/terminal_geometry.jl` allowed only for
native slab-axis/thickness metadata if directly required, and route-driver
summary plumbing allowed only if directly required in
`src/pqs_source_box_route_driver_helpers.jl` or
`src/pqs_source_box_diatomic_complete_core_shell.jl`.
For both `PQSLowering` and `WhiteLindseyLowering`, `:direct_midpoint_slab` and
`:outer_mismatch_slab` must not lower to `:direct_slab_identity_cpb` or
`:direct_boundary_slab_identity_cpb`. They must lower through the same compact
retained-block function and inputs. The unit slice scale is `ns x ns x 1`
after standard one-dimensional COMX/product compression; an outer-mismatch
region of thickness `t <= ns` should be decomposed or realized as a stack with
scale about `t * ns * ns`, not as one identity block. This identical-lowering
rule is only for thin slabs; real shells remain route-specific after common
shellification, and direct/core identity sectors remain identity. If a slab
stack thickness exceeds `ns`, source work must stop and report whether an
`ns x ns x ns` whole-block compression or a setup-error policy needs separate
approval. This lane does not approve driver changes, artifact/schema/reader
changes, RG/MWG/IDA changes, route skeleton redesign, broad terminal
realization redesign, direct slab deletion, committed Cr2 fixtures/tests, or
Cr2 workflow.

`HP-COMP-THINSLAB-META-FN-01` and `HP-COMP-THINSLAB-META-TEST-01` approve
only the live terminal-shellification metadata/scaffold inventory update in
`src/cartesian_terminal_shellification_geometry.jl`, with
`src/pqs_source_box_route_driver_helpers.jl` and
`src/pqs_source_box_diatomic_complete_core_shell.jl` allowed only in support
of the already approved thin-slab lowering pass. The stale owner is
`_cartesian_terminal_region_unit_mapping(region)`: midpoint slabs,
outer-mismatch slabs, and `:angular_z_extension_slab` must map to the compact
thin-slab lowering category, not direct identity categories. This file remains
metadata-only: it may add minimal compact thin-slab inventory/count vocabulary
needed for route summaries to agree with terminal lowering, but it must not
materialize coefficients, construct Hamiltonian data, add artifact/report
payloads, change shellification geometry, redesign route skeletons, add a new
reporting framework, or reintroduce direct identity slab lowering under a new
name. Direct core and atom-contact core identity mappings remain unchanged.

`HP-COMP-FACEPROD-FN-01` and `HP-COMP-FACEPROD-TEST-01` approve only a
neutral `CartesianFinalBasisRealization` face-product terminal helper seam for
compact thin-slab lowering and WL facet reuse. Approved source files are
`src/cartesian_final_basis_realization/terminal_face_product_blocks.jl`,
`src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`,
`src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`,
and `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`.
The helper is private/module-internal, reuses `_nested_doside_1d(...)` and
`_nested_face_product(...)`, supports normal axes `:x`, `:y`, and `:z`, one
fixed normal-axis index or an ordered stack of fixed normal indices, and a
caller-supplied retained count normally equal to public `ns`. WL facet terminal
realization should use this neutral helper, and `HP-COMP-THINSLAB-*` may use it
for midpoint, outer-mismatch fallback, and angular z-extension slabs. This
lane does not approve driver changes, public API/export, artifact/schema/
reader/provenance changes, shellification changes, terminal lowering policy
changes by itself, route skeleton changes, RG/MWG/IDA/Hamiltonian/raw-block
changes, old high-order workflow revival, committed tests/fixtures, Cr2
workflow, duplicate face-product assembly, PQS-specific thin-slab projection,
or treating thin slabs as WL boundary strata for naming convenience.

`HP-COMP-ANGBOX-FN-01` and `HP-COMP-ANGBOX-TEST-01` approve only
z-axis diatomic angular-balanced shellification in
`src/cartesian_shellification/terminal_geometry.jl`, with narrow route-driver
or complete-core-shell summary plumbing only if directly required. The
ordinary index-layer shared-shell body remains valid, but when shared-shell
growth stops with transverse axes saturated and bond-axis parent support
remaining, shellification must emit the bond-axis leftovers as planned
`:angular_z_extension_slab` stack regions with native axis/side/thickness,
stack, bond-axis, angular-rule, margin, transverse-scale, and extension-size
metadata. The same thin-slab concept applies to midpoint slabs, planned
non-boundary angular z-extension slabs, planned boundary angular z-extension
slabs, and unexpected outer-mismatch fallback slabs. Planned z-extension
stacks lower through `HP-COMP-THINSLAB-*`, which remains a separate blocker;
real shells remain route-specific after common shellification. This lane does
not approve driver changes, artifact/schema/reader/provenance changes,
terminal lowering, retained units, transform contracts, terminal realization,
RG/MWG/IDA/Hamiltonian/raw-block changes, route-family-specific PQS/WL
geometry, committed tests/fixtures, Cr2-specific branches, or Cr2 workflow.

`HP-MCOMX-FILE-01`, `HP-MCOMX-OBJ-01`, `HP-MCOMX-FN-01`,
`HP-MCOMX-WIRE-01`, and `HP-MCOMX-TEST-01` approve only the mainline
mapped-COMX source-span facility described in
`docs/src/developer/designs/cartesian_hamiltonian_producer/mapped_comx_source_span.md`.
The approved owner is the existing nested doside / COMX source-span seam, with
allowed files `src/cartesian_nested_faces.jl`,
`src/cartesian_pair_block_materialization/pqs_source_axis_transforms.jl`, and
`src/cartesian_raw_product_sources/axis_transform_facts.jl` /
`src/cartesian_raw_product_sources/records.jl` only for compact
`AxisSourceTransformFact` provenance or accessors if needed.
The first source rule is protected physical `P2` plus mapped Chebyshev
enrichment `T_k(s_lambda(u))` with `lambda = 0.5`, normalized local
`u in [-1, 1]`, no `sqrtJ`, and existing physical-position COMX cleanup.
High-order is a consumer/benchmark lane for the installed option, not the owner
of a duplicate implementation. This lane does not approve
`src/cartesian_raw_product_sources/mapped_comx_source_span.jl`,
`CartesianRawProductSources` numerical builders, a second COMX wrapper,
source-default changes, public API/export changes, canonical driver input
changes, artifact/manifest/reader changes, Hamiltonian/one-body/IDA/MWG/RG/
raw-block/solver changes, `protected_degree != 2`, injection/Ylm mechanisms,
mapped-`s` localization as production gauge, high-order scaffolding imports,
committed Cr/Cr2 fixtures, or Cr2 workflow.

`HP-MCOMX-TERM-FN-01` approves only terminal-basis wiring for mapped-COMX
source-axis transform facts in
`src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`,
with
`src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`
allowed only for import/include cleanup if directly required. In `_shell_seed`,
production source may prefer carried
`contract.metadata.raw_product_source_axis_transform_facts` when present,
validate exactly three materialized `AxisSourceTransformFact`s whose intervals,
mode dimensions, and coefficient matrix shapes match the shell source box, and
build full shell seed coefficients from those axis coefficient matrices before
continuing through the existing boundary-mode selection, owned-support row
restriction, shell-local Lowdin, sign canonicalization, and support validation.
The ordinary fallback through `_nested_projected_q_shell_full_sides(...)` must
remain when materialized facts are absent. `HP-MCOMX-TERM-TEST-01` approves
only ordinary PQS H2 regression, mapped source-span probe validation, a focused
terminal seam check showing mapped shell coefficients are basis-defining, and
the H2 supplemented RG endpoint if the touched path crosses it. This lane does
not approve driver inputs, source defaults, artifacts/manifests/readers,
Hamiltonian/IDA/MWG/RG/raw-Gaussian/solver/EGOI/Cr2/high-order workflow
changes, a second COMX wrapper, or committed tests/fixtures.

`HP-MCOMX-DRV-FN-01` approves only a compact `source_span` construction choice
in the canonical driver and staged base/facade plumbing. Approved files are
`bin/cartesian_ham_builder.jl`, `src/cartesian_base_hamiltonian.jl`, and
`src/pqs_source_box_route_driver_helpers.jl` only for narrow propagation to the
already-approved mapped-COMX source-axis transform fact path. Public driver
values are `:ordinary` and `:mapped_comx`; the default is `:ordinary`.
`:mapped_comx` is currently a PQS source-box option and must reject clearly for
`nesting = :wl` unless a later WL-specific amendment approves otherwise. The
driver may include `source_span` in its editable defaults, trusted input-file
keys, command-line overrides, and compact contract printout. This lane must not
add route records, terminal-lowering changes, artifact/schema/manifest/reader
changes, Hamiltonian/IDA/MWG/RG/raw-Gaussian/solver/EGOI/Cr2/high-order
workflow changes, source-span default changes, another COMX path, or committed
tests/fixtures. `HP-MCOMX-DRV-TEST-01` approves only default ordinary driver
artifact/readback, mapped-COMX H or He PQS driver smoke proving carried facts
are basis-defining, bounded ordinary-vs-mapped He supplemented/MWG/IDA driver
comparison if practical, H2 RG endpoint, and no Cr2 run.
Post-installation He/PQS evidence from 2026-06-26 found `n_s = 5`
mapped-COMX too weak for robust all-electron scalar capture. Do not promote
mapped-COMX or `n_s = 5` mapped-COMX to default behavior without a later
docs-only amendment based on bounded `n_s = 6`/`7` He H1/IDA evidence.

`HP-WLTERM-FILE-01`, `HP-WLTERM-FN-01`, and `HP-WLTERM-WIRE-01` approve only
the narrow White-Lindsey terminal-basis seam needed by `nesting = :wl`.
Approved files are `src/pqs_source_box_route_driver_helpers.jl`,
`src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`,
optional
`src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`,
and the matching include in
`src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`.
The goal is for the existing `:white_lindsey_low_order` route to produce the
same `CartesianTerminalBasisRealization` consumed by the staged Hamiltonian
path, using existing terminal support, retained-rule, and transform records.
This does not approve old WL H1/H1+J materialization adaptation, route
skeleton changes, shellification or retained-selection changes, raw-block
changes, Residual Gaussian/MWG/IDA changes, supplemented WL behavior,
diagnostics/status/report expansion, artifact schema changes, public API/export
changes, or Cr2 workflow. `HP-WLTERM-TEST-01` approves default `:pqs`
artifact/readback, `nesting = :wl` base atom artifact/readback, optional
`nesting = :wl` base H2 artifact/readback only if current native records
support it, H2 R3 PQS endpoint if terminal realization is touched, and no Cr2
run.

`HP-DRV-STAGE-FN-01` approves only a narrow non-exported, non-underscored
driver-facing staged producer surface. Approved source files are
`src/cartesian_base_hamiltonian.jl`,
`src/pqs_source_box_low_order_materialization.jl`,
`src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`, and
`src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`.
`src/cartesian_base_hamiltonian.jl` remains the primary driver-facing owner;
the lower-level files are approved only for behavior-preserving
operator-class stage factoring in their existing domains. That surface may
factor the existing base and residual-GTO/MWG facade bodies so the canonical
driver can execute visible physics-level stages: base working basis / terminal
realization, product/moment operators, unit-nuclear attraction operators,
electron-electron / IDA or residual-MWG interactions, base Hamiltonian
assembly, Gaussian supplement, residual Gaussian augmentation, supplemented
Hamiltonian assembly, followed by the existing artifact writer/readback stage.
The staged surface must be separate named construction-stage functions, not one
all-in-one replacement wrapper. It may return existing domain objects and small
fixed-key ephemeral stage products needed by the next stage, but it must not
create a public API/export, route-stage object, report, status/result payload,
metadata field cloud, runtime-keyed field group, persistent cache, artifact
schema, raw-block switch, allocation probe, per-kernel timing framework, or
source files outside the approved paths.

`HP-DRV-STAGE-WIRE-01` allows `bin/cartesian_ham_builder.jl` to call that staged
producer surface as separate visible top-level stage calls and time/print the
visible physics stages. Driver timing must distinguish product/moment,
unit-nuclear, and electron-electron stages. It does not approve one opaque
staged wrapper call, driver calls to underscored package helpers, or old route
stages such as `cartesian_parent`, `cartesian_shells`, `cartesian_units`,
`cartesian_pair_terms`, or `cartesian_assembly`.

`HP-R3U-ZDI-FN-01` relaxes the supplemented facade from hardcoded H2/Be2 checks
to explicit homonuclear two-center z-axis validation in
`src/cartesian_base_hamiltonian.jl`. Required inputs remain explicit: atom
symbols, nuclear charges, `nup`, `ndn`, geometry, base basis controls,
supplement labels, and optional trusted supplement `basisfile`. This does not
approve heteronuclear systems, non-z-axis orientation, ECP, charged systems,
solver workflow, public export/API redesign, artifact schema changes,
metadata/status/report fields, or Cr2-specific branches/defaults/fixtures.

`HP-R3U-ZDI-WIRE-01` allows `bin/cartesian_ham_builder.jl` supplemented mode
to call the supported supplemented facade through the compact `HP-DRV-*`
workflow. `HP-R3U-ZDI-TEST-01` approves H2 and Be2 facade/driver artifact
write/readback plus optional ignored/user-run Cr2 stress after H2/Be2 pass. It
does not approve committed Cr2 fixtures, committed Cr2 tests, Cr2-specific
workflow, or new committed tests.

`HP-DRV-ATOM-FN-01` and `HP-DRV-ATOM-WIRE-01` approve only explicit
origin-centered one-center base atom workflow in `bin/cartesian_ham_builder.jl`.
The driver may normalize explicit `atom_symbols`, `nuclear_charges`,
`atom_locations`, `nup`, `ndn`, visible one-center basis fields, and call the
existing `cartesian_base_hamiltonian(system; basis, hamfile)` facade where that
facade already supports the atom. Current validation remains origin-centered H.
These driver IDs do not approve edits to `src/cartesian_base_hamiltonian.jl`;
producer-side atom support is governed separately by `HP-R1-ATOM-*`, and
supported supplemented one-center atoms are governed by
`HP-COMP-SUPPATOM-*`. They also do not approve translated atoms, element
lookup/default tables, ECP, solver workflow, artifact schema changes, public
API/export changes, route diagnostics, metadata/status/report fields,
committed atom fixtures, or committed tests. `HP-DRV-ATOM-TEST-01` approves
only H atom driver artifact write/readback and optional ignored negative
checks.

`HP-DRV-ATOM-CLEAN-01` approves only removing hidden
`d = vars[:core_spacing]` from one-center atom basis construction in
`bin/cartesian_ham_builder.jl`. Public inputs, defaults, overrides, hooks,
timing labels, visible stage sequence, artifact schema, and the driver contract
must remain unchanged. This cleanup does not approve source/kernel changes,
diagnostics, new hooks, committed tests/fixtures, Cr2 workflow,
`:white_lindsey_low_order` retirement, test/tool route-input cleanup, or edits
outside `bin/cartesian_ham_builder.jl`.

Do not add private route-stage controls, stop-after internals, ladder probes,
stage markers, fixture hacks, diagnostic knobs, underscored package helper
calls, raw-block provider switches, report/status/payload dumps, metadata field
clouds, allocation probes, benchmark harness behavior, solver/RHF/ECP/EGOI/
HamV6 workflow, public API/export changes, artifact schema changes, committed
driver-input fixtures, committed tests, artifact schema dumps, supplemented
atoms, Cr2-specific driver branches, or Cr2-specific workflow support. Generic
explicit Cr2 may be used only as an ignored/user-run homonuclear z-axis stress
through `HP-R3U-ZDI-*` after H2/Be2 validation. Diagnostics and ladder probing
belong in `tools/` or ignored `tmp/work` probes, not in the canonical driver.

The `HP-RETIRE-CCS-RHF-*`, `HP-RETIRE-DRV-MAT-*`, and
`HP-RETIRE-LADDER-RUNNERS-*` lanes are completed and closed. Commits
`28e9b2c84`, `e2e164e9b`, and `77fa2700b` removed the RHF payload stack, old
materialization/report/save wrapper workflow, and dangling ladder runners.
Their registry/design entries are historical deletion records, not active
source authority. Do not restore those files, wrappers, tools, reports,
payloads, or compatibility adapters without a new docs-only amendment.

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
