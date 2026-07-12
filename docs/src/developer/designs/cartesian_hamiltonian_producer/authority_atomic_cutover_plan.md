# Cartesian Authority Atomic Cutover Plan

Status: independently reviewed planning material. This page does **not**
authorize the cutover, flag changes, live-view replacement, or source work.

## Purpose

Replace the temporary dual authority system with one machine-readable
record-level authority and two deterministic human views in one Git commit.
The cutover must reduce documentation carrying cost without changing any
scientific contract, execution membership, owned surface, or current blocker.

The immutable semantic-review input is:

- repository commit
  `9b283e16cf75b3264b00bc3f34d15e21bbb92a96`;
- candidate SHA-256
  `30cf4ed840b00c09da39ba4e15b3cb6c3d2c1376263877a2e7778f5a15bef716`;
- schema version `3`;
- `232` authority records;
- `44` hashed referenced documents;
- `150` derived execution IDs;
- marked `AGENTS.md` block SHA-256
  `35c944b5bda9225e2e33a1fc0c712e236a5147bc0679fa1fd459d7bfb8914ec5`.

Commit `9b283e16...` is the semantic baseline, not the future cutover parent.
The execution-authorization handoff must pin the exact accepted Pass 395 commit
as the cutover parent. That parent must descend from the semantic baseline and
must retain byte-identical candidate, 44 hashed documents, prose registry,
`AGENTS.md` authority section, and legacy checkers/snapshots. Any authority
input difference requires a new explicit parity review. A regenerated digest
is not review.

## Post-Cutover Authority Model

`authority.toml` will be the sole record-level authority for:

- ID and title;
- lifecycle and grant;
- surfaces and exact owned paths;
- dependency IDs;
- canonical, historical, and evidence document references;
- evidence references;
- scope;
- derived execution membership.

`registry.md` and the marked `AGENTS.md` whitelist will be deterministic,
checked views. They grant no independent authority and cannot be reverse-parsed
to reconstruct or repair `authority.toml`.

Manually maintained `current.md`, `invariants.md`, repo operating policy in
`AGENTS.md`, and linked canonical contracts may restrict machine-authorized
work. They may never add or broaden an ID, grant, surface, path, or scope.
A missing fact, conflict, stale generated view, or checker failure means no
producer work and requires a docs-only repair.

Lifecycle alone never grants execution. An ID is execution-listed only when:

1. its grant is `implementation`, `maintenance`, `preservation`, or
   `retirement`;
2. it owns at least one `source`, `tests`, `tools`, `artifacts`, or `driver`
   surface; and
3. its lifecycle/grant combination passes the schema contract.

Dependencies remain non-authorizing and non-topological.

## Exact Promotion

Rename:

```text
authority_candidate.toml -> authority.toml
```

Preserve every `documents` and `records` value and their order exactly. Keep
schema version `3`. Change only:

```toml
artifact_kind = "cartesian_authority"
authoritative = true
authorization_complete = true
```

Canonical serialization of that exact promotion must have SHA-256:

```text
6057ef50fd48ff329c9e7226665ea29f59ae4e35db3f43cb71bf29f5c4c2175e
```

No schema bump, provenance field, stored execution boolean, section field,
owner alias, supersession field, or compatibility key is approved. Such a
change would create unreviewed semantics.

## Generated Registry

Replace the hand-written `registry.md` with a checked projection generated
from one captured `authority.toml` snapshot. It must contain:

- a generated/do-not-edit banner and authority SHA-256;
- one lexically ordered record per ID;
- every schema field, including empty dependencies/evidence/paths;
- canonical document file links plus the recorded heading text;
- exact execution-membership display derived from the grant/surface rule.

Use the proven compact renderer: one `### ID - title` heading per record,
file-level document links, and the recorded heading as visible text. The
checker must continue validating that each recorded heading exists exactly
once, and Documenter must resolve every file link. Do not add an index, stored
anchor, fragment field, or new slug algorithm.

Design decision: deep links are intentionally waived. Exact IDs remain direct
record headings and search keys; adding a second slug/fragment contract would
increase checker complexity without adding authority safety.

The current registry's non-record `R3 Be2/Cr2 readiness guardrail` is stale
chronology and stays only in Git history. The `Composition family index` is
already represented by `current.md` and the canonical composition contract.
Neither becomes a new machine field.

Git history, the Pass 392/394 reports, and the manager log preserve the old
registry. Do not copy the full 4,851-line prose registry into another tracked
archive.

## Generated AGENTS View

Generate only the marked execution-whitelist block from `authority.toml`. It
must contain the authority source path and digest plus exactly 150 unique,
lexically sorted IDs. All derived wording belongs inside the marked region so
the renderer owns every generated byte.

Keep repo operating policy in `AGENTS.md`, including:

- startup and identity rules;
- fail-closed authority behavior;
- exact test-path and private-helper limits;
- runtime, test, deletion, and physics discipline;
- anti-bloat, driver, artifact, performance, and reporting policy that remains
  architecture-wide rather than ID-specific.

The Pass 395 audit found no rule requiring relocation. At semantic baseline
`9b283e16...`, delete these inclusive ranges wholesale:

- lines `691-1120` (`430` lines), SHA-256
  `f61a971952df45daaf95d8f1530ab937d57627588de955867661823c6a4f3755`;
- lines `1288-1723` (`436` lines), SHA-256
  `4b90527f8af92e853d97912fbafc4f3d74e8f297c22c8219ed0551043ed2bb7f`.

The first range contains terminal/R1, route/R3, RG/protected/reference, raw
block, exact-operator, and endpoint per-ID prose. Schema-v3 records,
`current.md`, `invariants.md`, and linked contracts already own it. The second
contains driver/due-diligence, source-mode, composition/WL/atom, common-shell,
mapped-COMX, staged-producer, and retirement prose with the same complete
owners. The allowed relocation set is empty.

Two stale statements are deliberately deleted: copied local drivers are
already excluded by exact tracked-path ownership, and the blanket prohibition
on supplemented atoms contradicts implemented registered composition cells.
Do not relocate either statement.

All 44 hashed documents are frozen during cutover. A range-hash mismatch or an
architecture-wide rule found only in deleted prose stops cutover for separate
reconciliation. It may move only to unhashed `current.md` or `invariants.md`
after review; any canonical-contract or authority-record edit requires renewed
parity review. Do not preserve the 866-line duplicate as another active or
historical copy; Git history and the completed cutover record are the archive.

The final `AGENTS.md` wording must state that its generated block is a checked
view and that surrounding policy can restrict but never grant machine-absent
authority.

## Permanent Checker

Rename and refactor:

```text
docs/check_cartesian_authority_candidate.jl
    -> docs/check_cartesian_authority.jl
```

Use module `CartesianAuthority`; add no compatibility wrapper. The checker must
retain the reviewed schema/path/hash/heading/grant/lifecycle validation and:

- require schema `3`, artifact kind `cartesian_authority`, and both true flags;
- load one immutable authority snapshot per operation;
- validate canonical TOML serialization;
- validate document hashes and unique recorded headings;
- validate tracked path states, containment, and symlink rejection;
- validate dependency closure and grant/lifecycle/surface compatibility;
- derive the execution set;
- render the complete registry and exact marked whitelist block in memory;
- byte-compare both committed views in `--check` mode;
- re-read captured inputs before success;
- reject any legacy live transition/shadow artifact;
- fail rather than use generated views, Git parents, cached bytes, or prose as
  fallback authority.

Permanent operation modes are:

```text
--check          read-only authority and committed-view validation
--self-test      permanent negative fail-closed regression coverage
--render <dir>   deterministic external view rendering only
```

`--render` must reject repository-contained and symlinked destinations. It
writes a complete registry and marked-block payload only to a new external
staging directory. It never edits live repository files and is not activation.
Install reviewed bytes only in the isolated cutover worktree; Git commit/ref
activation remains the atomic transaction. Remove rehearsal, manifest,
transition-binding, live-writer, and compatibility modes.

Permanent self-tests must cover missing/corrupt authority, false flags, unknown
schema, malformed grants/paths/dependencies, path escape and symlink attempts,
document/hash/heading drift, marker corruption, stale/partial generated views,
legacy live files, and concurrent mutation. CI must run them.

## Filesystem And Git Transaction

Perform the cutover in a separate machine-local clean worktree pinned to the
exact accepted Pass 395 commit named by the later execution authorization. Do
not operate inside the Dropbox-synced shared worktree with unrelated WIP.
Before generation require:

- exact authorized cutover-parent commit;
- proof that its authority inputs are byte-identical to semantic baseline
  `9b283e16...` and candidate digest `30cf4ed8...`;
- empty tracked, staged, and untracked status;
- passing six legacy authority gates;
- exact `232/44/150` counts and marked-block hash.

Generate and validate all proposed outputs before creating one cutover commit.
Filesystem writes are preparation only. The single Git commit and subsequent
branch-ref update define activation.

The cutover commit must contain the complete transition. A partial subset is
forbidden.

### Promote Or Transform

- `.github/workflows/docs.yml`;
- `AGENTS.md`;
- `docs/src/developer/designs/cartesian_hamiltonian_producer/authority_candidate.toml`
  -> `authority.toml`;
- `docs/check_cartesian_authority_candidate.jl` ->
  `docs/check_cartesian_authority.jl`;
- `docs/make.jl`;
- `docs/src/developer/designs/cartesian_hamiltonian_producer/registry.md`;
- `docs/src/developer/designs/cartesian_hamiltonian_producer/README.md`;
- `docs/src/developer/designs/cartesian_hamiltonian_producer/current.md`;
- `docs/src/developer/designs/cartesian_hamiltonian_producer/invariants.md`;
- `docs/src/developer/designs/cartesian_hamiltonian_producer/reviews/README.md`;
- this `authority_atomic_cutover_plan.md`, compressed in the same commit to a
  short completed cutover record with old/new hashes, counts, validation, the
  exact parent OID, the phrase `this commit` for the cutover itself, and
  rollback instructions;
- `docs/src/developer/test_suite_reorganization_plan.md`;
- `docs/src/developer/pqs_manager_running_log.md`;
- `test/docs/runtests.jl`;

No separate cutover report is added. This plan's compressed completed record,
the manager log, Pass 392/394 reports, and Git history are the evidence.

### Delete

- `docs/src/developer/designs/cartesian_hamiltonian_producer/registry_whitelist_shadow.toml`;
- `docs/src/developer/designs/cartesian_hamiltonian_producer/authority_transition_snapshot.toml`;
- `docs/check_cartesian_authority_shadow.jl`;
- `docs/check_cartesian_authority_transition.jl`.

Do not retain aliases, shims, old include paths, transition CI, rehearsal
output, or a reverse parser. Historical review reports and prior manager-log
entries remain unchanged.

No source, scientific test, artifact, solver, driver, or any of the 44 hashed
document files belongs in the cutover commit. An uncaptured architecture-wide
guardrail stops cutover for separate reconciliation. The unrelated
`src/hamiltonian_corrections.jl` WIP and successor handoffs are excluded.

## Permanent CI And Tests

Replace the transition job with a required `cartesian-authority` job using
read-only contents permission and no package instantiation:

```sh
julia --startup-file=no --project=docs \
  docs/check_cartesian_authority.jl --check
julia --startup-file=no --project=docs \
  docs/check_cartesian_authority.jl --self-test
```

Make the existing `docs` job depend on `cartesian-authority`, so docs cannot
pass or deploy while authority fails. Remove all authority-checker includes
from `docs/make.jl`; the required dependency owns that CI gate. One permanent
checker is sufficient, and no rehearsal is generated in CI.

Reduce `test/docs/runtests.jl` to the permanent authority/workflow contract and
delete transition filename, warning, snapshot, marker-count, and stale
checker-source assertions. Keep compact assertions that both permanent commands
are wired. Delete the three existing stale phrase-inventory assertions rather
than rewriting them to another prose snapshot. All focused docs tests must
pass; do not hardcode a permanent test count.

## Acceptance Gates

Before commit:

1. prove parsed equality of all 44 document entries and 232 records against the
   reviewed candidate;
2. prove the authority promotion changes only artifact kind and two flags;
3. verify the exact promoted authority digest;
4. verify exact 150-ID execution membership and no set delta;
5. run `--check`, `--self-test`, and focused negative mutations of authority,
   registry, and `AGENTS.md` generated bytes;
6. render to two fresh external directories and require byte-identical output;
7. install those reviewed bytes in the isolated worktree, create the proposed
   commit, then require a fresh render to equal the committed views and produce
   zero Git diff;
8. validate unique generated ID headings, recorded heading uniqueness, and all
   generated file-level document links;
9. run focused docs tests and local Documenter;
10. require `git diff --cached --check`, exact name-status review, no unstaged or
   untracked files, and no active legacy references outside historical evidence;
11. require roughly 5,200 or more net lines removed; less shrinkage requires a
    carrying-cost review for retained transition/prose machinery.

The reviewed digests and `232/44/150` counts are one-time cutover assertions.
Do not hardcode them in permanent checker behavior.

Then create one commit. Validate that commit from a second fresh detached
checkout before advancing the live branch ref. The clean checkout must pass the
canonical authority job, docs tests, Documenter, deterministic rendering, and
legacy-absence scans.

Before activation, in another temporary checkout at the proposed cutover
commit, run `git revert --no-commit <cutover-commit>` and require the resulting
index tree (`git write-tree`) to equal the exact cutover-parent tree.

## Activation And Rollback

No agent may treat generated files, a dirty worktree, Dropbox synchronization,
or a local/staging/review commit as activation. The only live ref is
`refs/heads/main` on `origin`. The execution handoff must pin its exact expected
parent OID. Activation is an ordinary fast-forward push of the one complete
cutover commit to that ref; it must fail if the remote ref is no longer the
pinned parent. Review, PR, and temporary refs are non-activating.

Freeze producer authority edits and source work during the acceptance window.
If post-activation CI fails, perform no producer work and revert the complete
cutover commit. Rollback must restore the exact parent tree and all six legacy
gates. Never repair by flipping one flag, restoring one prose file, retaining a
mixed checker set, or falling back to last-known-good generated views.

Before activation, rollback means discarding the proposed commit. After later
authority changes exist, rollback requires a separately reviewed whole-tree
reconciliation rather than blindly restoring stale prose.

## Stop Rules

Stop without cutover if:

- the execution parent differs from the explicitly authorized Pass 395 commit;
- any authority input differs from semantic baseline `9b283e16...` or candidate
  digest `30cf4ed8...`;
- either audited `AGENTS.md` deletion-range hash differs;
- any document/record value changes beyond the three promotion fields;
- execution membership differs from the reviewed 150 IDs;
- generated views cannot be deterministic and byte-checked;
- any active permission remains only in prose;
- retained hand-written `AGENTS.md` text grants or broadens an ID-level surface;
- any old checker/shadow/transition path remains live;
- generated file links do not resolve or recorded headings are nonunique;
- CI cannot fail closed on partial or stale state;
- the commit overlaps scientific source, artifacts, drivers, or unrelated WIP;
- whole-commit revert cannot restore the parent tree.

## Explicit Non-Goals

- no scientific, numerical, source, test-surface, artifact, API, or workflow
  authority change;
- no new IDs, grants, paths, dependencies, or execution members;
- no schema redesign or metadata expansion;
- no public producer or solver behavior change;
- no authority inference from source callers or Git history;
- no archival copy of the full old registry or duplicated AGENTS authority;
- no cutover execution under this planning page.

## Decision Boundary

This plan passed independent transaction, generated-view, and
tooling/minimality review. A later explicit user/design-manager decision may
authorize one implementation pass that follows this transaction exactly and
pins the accepted Pass 395 commit as execution parent. Approval of this plan
alone must not be interpreted as approval to execute the cutover.
