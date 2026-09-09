# Documentation Deployment

## Released Documentation Refresh

Pass 615 accepts the Pass 614 documentation-only implementation and publication
for released v0.2.0 under `HP-PQS-DOCS-TAGDEPLOY-FN-01/TEST-01`.
The records are implemented/completed, maintenance-only; the single publication
grant is exhausted. No further snapshot publication or stable promotion is granted.
This section supersedes older stable-alias restrictions below and earlier
maintenance records' stable-preservation wording only for this transaction.
No package version, tag, release, artifact, API, source, dependency, example,
or numerical change is authorized.

Accepted implementation: `6fd436cda4953c3886280761556575a5fc4d8660`.
CI `34137792695`, bootstrap-hold Docs `34137792698`, publication
`34138301136`, and Pages `34138461292` passed. Published gh-pages commit:
`40912e64c2d66cf8928429b0687e89b787e72034`; frozen `release-0.2.0` tree:
`c76c64547dfbc2c00ec0a0e5efebee7ccb5cf33a`. Stable points to that snapshot;
the root redirects to stable. Original v0.2.0/RC1/RC2 trees remain respectively
`8990d998955e92b2761d8654d24863d578802a4b`,
`45568ee3d836b68ee8536b21c92f982eab7c35b4`, and
`255824282cab1670c3e3daed1a7e369a0d6febc0`.
The 25-URL and 14-README-destination checks passed; tags/releases are unchanged.
Closeout must verify these trees and the stable pin after its normal-main Docs
deployment. That deployment may update dev only; it is not another publication.

### Compatibility Evidence

At `a77e483d57f9b018566223f9352d14f95ceb9fac`, all 346 curated reference
entries resolve against the immutable v0.2.0 source archive, with 279 already
documented there. The remaining 67 acquire documentation on main, not new API.
Project.toml and every example are byte-identical to release. Source review
found documentation, ownership relocations, exact Lanczos consolidation, and
private retirements; none adds a reader-facing operation absent from release.
The three post-release removed root names are not promoted by the reference.

Julia 1.12.6 loaded the released source from isolated machine-local scratch,
using the existing resolved dependency manifest (not a fresh dependency
resolution). Fresh processes passed examples 01/02/03/04/15/39/40 in
5.26/7.50/8.06/9.51/10.74/28.07/7.16 seconds, and the public residual-GTO
interchange owner passed 80/80 in 39.21 seconds including startup.
No larger matched comparison or full angular suite was rerun. Evidence is in
`/Users/srw/dmrgtmp/stable_docs_20260907/`; the counts and conclusions here
are canonical, not dependent on that machine-local location.

### Exact Publication Boundary

Use Documenter's existing deployment mechanism, not a custom site copier.
Publish a separately built, approved documentation snapshot to the new
`release-0.2.0` gh-pages folder. Map `stable => release-0.2.0` explicitly
as the first version-selector entry. Keep the existing v0.2 minor selector,
RC1, RC2, and dev entries. The snapshot is documentation, not a package tag or
new package version. Its HTML canonical base is the public `/stable/` URL.
Preserve the entire original `/v0.2.0/`, both RC folders, and their canonical
URLs byte-for-byte. Never deploy directly through the existing stable symlink:
Documenter's copying operation can follow that link into the original folder.

Normal main builds still publish only dev with /dev/ canonical URLs. Future
tag builds using the amended policy publish only their own new exact version
folder, retaining the explicit stable mapping; a new release does not silently
advance stable. Existing version folders must not be overwritten by this
transaction. A later stable promotion requires its own compatibility decision.
Old frozen tag-workflow reruns are forbidden: their old policy is immutable and
cannot be hardened by a main-only edit. This grant neither reruns them nor
claims repository-setting enforcement against arbitrary historical workflows.

Add a narrowly validated workflow_dispatch context to the existing Docs
workflow and deployment job, with a required exact expected commit SHA.
Only dispatch on main at that SHA may request this v0.2.0 snapshot publication.
Reject wrong branch, SHA, version, mode, or malformed input before deploydocs.
Use existing deployment-step credentials/permissions; PRs remain build-only.
Serialize deployment writes with one job-level concurrency group and no
cancellation of an in-progress deployment. Add no workflow, job, secret,
credential, polling/chaining infrastructure, or configurable publication target.

The ordinary implementation push may build but must not publish while the
snapshot is absent: explicitly report this bootstrap hold, preserving the old
site/selector. The one manual dispatch then creates the snapshot and standard
Documenter selector/root redirect together. After that, ordinary main/tag
deployment must require and preserve the snapshot, not silently fall back to
the original release or dev. Validate a subsequent normal main deployment at
closeout. A missing snapshot after publication is an operational blocker, not
permission to regenerate it from arbitrary main. No second manual publication,
overwrite, force push, deletion, or repair is implicit if the first partially
succeeds; preserve state and report.

Build from the reviewed implementation revision, keeping the improved current
docstrings. Add concise reader labeling in docs/src/index.md and
docs/src/developer/index.md and, if needed, a build-context footer in make.jl:
this is revised documentation validated for package v0.2.0; Developer Notes,
private internals, source-layout pointers, and development history describe the
documented revision, not additional released API. Preserve experimental/expert
labels and distinguish interface support from scientific maturity. Source/edit
links must identify the reviewed documentation revision, not falsely pretend
new source paths or docstrings exist at the release tag. Do not turn the whole
site into an automatically moving stable alias to dev.

### Budgets And Acceptance

Allowed implementation paths: .github/workflows/docs.yml, docs/make.jl,
test/docs/runtests.jl, docs/src/index.md, docs/src/developer/index.md.
Preferred/hard additions: workflow 25/40, make.jl 50/80, focused tests 30/45,
reader labeling 12/20. No new tracked file, source edit, docs environment
dependency, general framework, numerical test, prose-lock suite, or public CI
matrix change. Replace affected old selector tests, do not retain conflicting
policies. Stop before implementation commit if reliable publication needs
broader machinery, source changes, or incompatible reader documentation.

Use installed Documenter behavior in isolated fixtures to check main, manual,
new tag, prerelease, missing snapshot, wrong SHA, canonical URLs, selector,
root redirect, and untouched original-folder hashes. Preflight already
simulated refresh/dev/v0.2.1/dev postprocessing with Documenter 1.17:
stable stayed on the snapshot and the original folder hash did not change.
Do not patch Documenter or mutate source tags for rehearsal.

Before publication freeze the implementation SHA, snapshot build identity,
original gh-pages subtrees, and remote annotated v0.2.0 identity:
object 722e8e8752a9d23f45e95d2f88e1749f9f3002e4, commit
adfcaba32d4db06d9d796d947276433717bd2d89, tree
f64ba21e06ff57e2b5e78d91214398115afbe8de. Require local package/docs,
authority/self-test/views, Documenter, YAML, log bound, and diff checks;
normal CI classification applies to the workflow edit, with no extra numerical
reruns beyond existing required CI. Preflight release examples need not repeat
when their source and reader contract remain unchanged.

After the single dispatch, require Docs/Pages success and verify via rendered
HTTP pages: stable corrected manual/reference content, exact stable canonicals,
root redirect to stable, stable/v0.2/RC2/RC1/dev selector, all README
destinations, and unchanged original v0.2.0/RC folders and tag/release identity.
Report publication run IDs and gh-pages commit, hashes, and any partial state.
Close this publication separately; repo-manager waits for this authority
commit and its checks before implementation.

## Historical Deployment Policy And Evidence

The following records describe previous deployments, not current permission
to reset stable. Their release identities and historical outcomes are retained.

The rendered Documenter site is the primary documentation surface:

- <https://srwhite59.github.io/GaussletBases.jl/dev/>

Local builds remain available with `julia --project=docs docs/make.jl`.

The `Docs` workflow enforces a least-privilege split:

- the Cartesian-authority job uses `contents: read`;
- pull requests build documentation with `contents: read` and never deploy;
- pushes to `main` use a separate deployment job with only `contents: write`;
- `GITHUB_TOKEN` and `GAUSSLETBASES_DOCS_DEPLOY=true` are exposed only to the
  deployment step.

## Tag-Aware Version Deployment

Tag-aware deployment is implemented and maintained under
`HP-PQS-DOCS-TAGDEPLOY-FN-01` and
`HP-PQS-DOCS-TAGDEPLOY-TEST-01`. This contract does not authorize a
package-version change, a tag, or a release.

The deployment classifier preserves three distinct contexts:

- a pull request is build-only with `contents: read`, receives no deployment
  credential, and cannot call `deploydocs`;
- a push to `main` builds and deploys the `dev` folder with canonical URL
  `https://srwhite59.github.io/GaussletBases.jl/dev/`;
- a pushed tag of the exact form `vMAJOR.MINOR.PATCH` or
  `vMAJOR.MINOR.PATCH-PRERELEASE`, accepted as a Julia/Documenter semantic
  version and containing no build-metadata suffix, builds and deploys the exact
  tag folder. For example, `v0.2.0-rc1` uses canonical URL
  `https://srwhite59.github.io/GaussletBases.jl/v0.2.0-rc1/`.

The workflow uses a `v*` trigger so that GitHub starts the check for a
candidate version tag, while `docs/make.jl` independently validates the full
tag before deployment. A non-version tag, a malformed `v` tag, an unsupported
event, or a mismatch between the resolved build context and deployment request
fails before `deploydocs`; it never falls back to `dev` or another
release folder.

Deployment uses Documenter's standard final-release selectors together with
explicit prerelease selector entries:

```julia
"v0.2.0-rc2" => "v0.2.0-rc2"
"v0.2.0-rc1" => "v0.2.0-rc1"
```

The RC1 self-mapping keeps its existing folder in `versions.js` across both
tag and later `main` deployments. Commit
`2b3c23970144aa030ae52b875a5cf01b32886b6e` implemented the RC2 self-mapping.
The completed tag lifecycle now retains the existing RC2 and RC1 folders on
later deployments. Neither entry creates an alias. Documenter still
excludes prerelease folders from the release set used to create or advance
`stable`. A later final `v0.2.0` tag may become `stable` under the standard
policy. No other prerelease entry, custom `stable` alias, or dynamic release
index is authorized. The frozen RC2 annotated-tag lifecycle is completed;
GitHub release publication remains unauthorized. These semantics follow
[Documenter's deployment criteria](https://documenter.juliadocs.org/stable/lib/public/#Documenter.deploy_folder)
and [versioned deployment policy](https://documenter.juliadocs.org/stable/man/hosting/#Documentation-Versions).

Maintenance is limited to `.github/workflows/docs.yml`, `docs/make.jl`, and
focused assertions in `test/docs/runtests.jl`. It may make this contract and
compact current status truthful, but it may not add a
helper file, release framework, custom credential, alternate host, manifest,
artifact, source/API change, version bump, citation/changelog change, or
scientific-document rewrite. The stop-and-report limits are:

- at most `30` added workflow lines;
- at most `40` added `docs/make.jl` lines;
- at most `50` added lines in the existing docs test;
- no new file.

Readability and fail-closed classification take precedence over approaching a
limit. If the behavior cannot fit without duplicated parsers, a new framework,
broader credentials, or hidden release policy, implementation must stop and
return to repo-design-manager.

Commit `31caa87d3b83599de7f7295678ee599209113552` implemented the Pass 495
repair with `7` added `docs/make.jl` lines and `22/1` existing-test
lines. `_DOCS_VERSIONS` is now the single explicit selector policy passed to
`deploydocs`; the workflow and credential model did not change.

Validation must simulate, without creating or pushing a tag, pull-request,
`main`, `v0.2.0-rc1`, `v0.2.0-rc2`, final `v0.2.0`, and malformed-tag
contexts. It must:

- inspect the exact deployment classification for each context;
- build locally under simulated `main` and `v0.2.0-rc2` environments, retain
  exact RC1 classification, and inspect the generated canonical links;
- demonstrate with an isolated folder fixture that the explicit RC1 and RC2
  self-mappings list both prereleases with `dev`, create no self-symlink, and
  do not select `stable`, while a later final release may;
- retain the current job- and step-level credential boundaries;
- pass authority check/self-test, generated-view parity, docs tests,
  manager-log bound, package load, Documenter, YAML/workflow inspection, and
  `git diff --check`;
- pass remote Docs and CI after the implementation push and confirm that the
  ordinary `main` deployment retains `/dev/`, lists RC1 in `versions.js`, and
  leaves `/stable/` absent.

The immutable annotated `v0.2.0-rc1` tag was pushed successfully at tag object
`a4284f0bf448fb9d717de26ccbe1e9fc16db5ed2`, peeling to
`1546c18d3058cce2b5051b50788cda3c12585e51`. Tag-triggered Docs run
`32295705338` published the exact canonical RC1 folder without changing
`/dev/` or creating `/stable/`. Its default version expansion omitted the
prerelease from `versions.js`; the explicit self-mapping above is the sole
authorized correction. The immutable tag must not be moved, replaced,
deleted, or retriggered as a repair mechanism.

Main-deployment Docs run `32302304167` and CI run `32302304185` passed at
the repair commit. The live `versions.js` lists exactly `v0.2.0-rc1` and
`dev`; both folders retain their exact canonical URLs, and `/stable/`
remains absent. With no final release present, Documenter sets its internal
`DOCUMENTER_STABLE` JavaScript fallback to the first listed version, RC1.
That fallback is not a `stable` selector entry, symlink, or published path
and grants no final-release status.

Commit `abee269eed7028c864fa18ae44b4b946af63dfcf` implemented the classifier,
canonical paths, workflow trigger, and focused tests within the authorized
files. Its deltas were `3/1` workflow lines, `32/7` `docs/make.jl` lines, and
`30/0` existing-test lines. Docs run `32264133694` and general CI run
`32264133755` passed; `gh-pages` advanced from `e164c91b5` to `8a43fc8f3`,
built from that commit, and the live `/dev/` contract page reflected the new
behavior. No package version, tag, release, source API, credential model, or
documentation schema changed.

Versioning, changelog and citation work, candidate tagging, GitHub release
creation, and final-release authorization remain separate decisions.

Commit `62a1a4821` restored this standard same-repository Documenter deployment
path. Docs run `32072728238` and general CI run `32072728319` passed; the
`gh-pages` branch advanced from `a9b74566e` to `255ea4ed4`, built from that
commit. The live producer current-status page returned successfully with the
corresponding current-main content.

No custom credential, repository-setting change, alternate deployment
framework, source/API change, or documentation-navigation redesign is part of
this contract. GitHub's repository browser remains a raw-file view; the
published site provides the intended rendered navigation.

## v0.2.1 Candidate Preparation

Pass 625 authorizes candidate preparation only under
`HP-PQS-PUBLIC-V021-FN-01` and `HP-PQS-PUBLIC-V021-TEST-01`.
Reviewed main is `57e61403262042923605f94a28e910dd2b543550`; released
v0.2.0 is `adfcaba32d4db06d9d796d947276433717bd2d89`.
Repo-manager's conversation-delivered readiness review is audit evidence, not
an independent grant. Comparison confirms unchanged Project declarations and
examples, the seven-line compatibility removal, and the local-tag fetch
collision. No additional RC is required absent a concrete candidate blocker.

### Compatibility restoration

Restore exactly the original three root exports:
`QiuWhiteResidualGaussianOperators`, `OneCenterAtomicNestedLayerStructure`,
and `TimedNestedFixedBlockBuild`. Preserve the existing alias and diagnostic
type unchanged; restore the original four-line parametric timing carrier from
v0.2.0, including `fixed_block::F` and `timings::TimeG.TimingReport`.
This supersedes the absence condition in export-integrity maintenance, not
the historical evidence. No renewed development, warning, shim, replacement
interface, field, method, or numerical behavior is authorized.

Add concise compatibility-only documentation at existing definitions/alias,
curate the three names in the existing export reference, and extend existing
mechanical public-surface and core owners. Verify exports, alias identity,
original type/field/constructor behavior, and documentation resolution.
Keep the exact five reserved undocumented names unchanged. Compare the
candidate export inventory to released v0.2.0; report any additional
compatibility discrepancy rather than inferring that private layout changes
are public guarantees.

### Tag verification repair

Only the tag-identity step of `.github/workflows/ci.yml` may change.
Fetch `refs/tags/${tag}` with `--no-tags` into a dedicated non-tag reference,
preferably `refs/verification/release-tag`. Resolve that reference's
`^{tag}`, verify its object type, then its `^{commit}` peel, event SHA,
tree against checkout HEAD, canonical version spelling and Project version.
Preserve the existing remote installation/package-load step and all other
CI behavior. A preexisting verification reference must fail closed; no force,
tag replacement, tag deletion, or silently trusted checkout tag is permitted.

Rehearse the actual verification commands using bounded disposable local Git
repositories: annotated remote with a colliding local lightweight tag succeeds
without changing it; missing remote tag, lightweight remote tag, wrong commit,
tree or version, malformed tag, and preexisting verification ref fail.
Use existing focused docs checks for the production wiring; keep temporary
fixtures outside the repository, with no new framework or checked-in owner.

### Candidate and documentation boundaries

Change only Project version to `0.2.1`. Add a concise `v0.2.1 (unreleased)`
changelog section above byte-identical release history: centered/displaced
kernel arithmetic repairs, fifth radial construction attempt with early exit
and opt-out, affected-case increased cost, truthful exhausted-refinement
warnings, and documentation/validation improvements. Relocations are not new
numerical implementations. Neither whole-matrix convergence nor energy-error
control is promised; preserve the near-origin inverse-radius limitation.

Keep the README's available installation pinned to v0.2.0 until publication;
at most add a clearly unreleased v0.2.1 candidate note. Never present an absent
tag as installable. Update docs footer/version checks to distinguish candidate
documentation from released API; retain the frozen refresh context's v0.2.0
label. Preserve `stable => release-0.2.0`, site-root behavior, original
version folders, and all refresh safeguards. No Docs workflow edit is needed.

### Exact implementation surfaces and budgets

No new implementation file. Added-line preferred/hard budgets:
- `src/GaussletBases.jl`: exactly three restored export lines.
- `src/cartesian/cartesian_nested_faces.jl`: original four-line definition
  plus documentation; `src/cartesian/cartesian_nested_atomic.jl` and
  `src/ordinary/ordinary_qw_types_and_bases.jl`: documentation only.
  Combined source documentation: 24/36 lines; executable restoration: seven.
- `.github/workflows/ci.yml`: 6/10 added lines, tag step only.
- `Project.toml`: one version substitution, no dependency/compat change.
- `CHANGELOG.md`: 18/24 added lines; `README.md`: 4/8.
- `docs/make.jl`: 5/8 added lines for truthful version labeling only.
- `docs/src/reference/export.md`: 10/16 added lines.
- `test/core/runtests.jl`: 12/20 added compatibility lines.
- `test/docs/public_surface_runtests.jl`: 4/6 added inventory lines.
- `test/docs/runtests.jl`: 30/50 added mechanical checks/substitutions.
Normal lifecycle evidence and mechanically required generated-view/digest
reconciliation are separate from these implementation budgets.

### Acceptance and handback

Repo-manager waits for this authority commit and its checks. Prefer one
implementation/candidate push so source/workflow/version changes receive one
full existing matrix. Reuse centered CI 34364324843, displaced CI 34377518052
(148 regressions and 64 complete-matrix checks), and radial CI 34392873423.
No new numerical policy, broad angular run, repeated paper example, or
benchmark campaign is authorized.

Require focused compatibility checks, tag-fixture rejection evidence, existing
three numerical gates (including Julia 1.10 Supported floor), package load,
docs_fast/full docs, authority/self-test, deterministic views, Documenter,
manager-log bound, YAML inspection and diff checks. Perform one isolated
Julia 1.12.6 candidate-archive installation/load and small public Example 01.
Freeze commit, tree, archive entry/byte counts and SHA-256; exclude a root
Manifest and both protected handoffs. Record exact deltas and workflow IDs.
Candidate acceptance awaits independent review; no tag identity exists yet.

Stop without an implementation commit if restoration changes original
definitions, another compatibility defect appears, tag checks weaken, stable
moves, budgets are exceeded, or broader machinery is necessary.

Tagging, publication, registration, release assets, stable promotion, and
old-tag reruns are not authorized. A later separate bounded promotion may
change the stable pin and associated selector/root-policy checks only after
the new immutable versioned documentation exists, is verified release-compatible,
and its canonical URLs, source links, README destinations and preserved old
folders pass review. That later transaction must keep future main/tag
deployments from undoing the selected pin. Standard60, represented-Hartree,
arbitrary-position work and other cleanup remain outside this packet.
