# Documentation Deployment

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
one explicit prerelease selector entry:

```julia
"v0.2.0-rc1" => "v0.2.0-rc1"
```

This self-mapping keeps the existing RC1 folder in `versions.js` across both
tag and later `main` deployments; it creates no alias. Documenter still
excludes prerelease folders from the release set used to create or advance
`stable`. A later final `v0.2.0` tag may become `stable` under the standard
policy. No other prerelease entry, custom `stable` alias, or dynamic release
index is authorized. These semantics follow
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

The Pass 495 repair is narrower than those standing limits. It may add only
the explicit RC1 self-mapping to the existing `deploydocs` call and focused
coverage in `test/docs/runtests.jl`: at most `8` added `docs/make.jl` lines and
at most `25` added existing-test lines. `.github/workflows/docs.yml` must not
change. If Documenter cannot retain RC1 without changing `stable`, the tag,
or the credential/deployment model, stop and report.

Validation must simulate, without creating or pushing a tag, pull-request,
`main`, `v0.2.0-rc1`, final `v0.2.0`, and malformed-tag contexts. It must:

- inspect the exact deployment classification for each context;
- build locally under simulated `main` and `v0.2.0-rc1` environments and
  inspect the generated canonical links;
- demonstrate with an isolated folder fixture that the explicit RC1
  self-mapping lists `v0.2.0-rc1`, creates no self-symlink, and does not select
  `stable`, while a later final release may;
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
