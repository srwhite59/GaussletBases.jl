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

Tag-aware deployment is approved under
`HP-PQS-DOCS-TAGDEPLOY-FN-01` and
`HP-PQS-DOCS-TAGDEPLOY-TEST-01`, pending implementation. This amendment does
not authorize a package-version change, a tag, or a release.

The deployment classifier must preserve three distinct contexts:

- a pull request is build-only with `contents: read`, receives no deployment
  credential, and cannot call `deploydocs`;
- a push to `main` builds and deploys the `dev` folder with canonical URL
  `https://srwhite59.github.io/GaussletBases.jl/dev/`;
- a pushed tag of the exact form `vMAJOR.MINOR.PATCH` or
  `vMAJOR.MINOR.PATCH-PRERELEASE`, accepted as a Julia/Documenter semantic
  version and containing no build-metadata suffix, builds and deploys the exact
  tag folder. For example, `v0.2.0-rc1` uses canonical URL
  `https://srwhite59.github.io/GaussletBases.jl/v0.2.0-rc1/`.

The workflow may use a `v*` trigger so that GitHub starts the check for a
candidate version tag, but `docs/make.jl` must independently validate the full
tag before deployment. A non-version tag, a malformed `v` tag, an unsupported
event, or a mismatch between the resolved build context and deployment request
must fail before `deploydocs`; it must never fall back to `dev` or another
release folder.

Deployment continues to use Documenter's standard version policy. A
prerelease folder remains independently addressable but is excluded from the
set used to create or advance `stable`. A later final `v0.2.0` tag may become
`stable`; no custom alias policy is authorized. These semantics follow
[Documenter's deployment criteria](https://documenter.juliadocs.org/stable/lib/public/#Documenter.deploy_folder)
and [versioned deployment policy](https://documenter.juliadocs.org/stable/man/hosting/#Documentation-Versions).

The implementation is limited to `.github/workflows/docs.yml`, `docs/make.jl`,
and focused assertions in `test/docs/runtests.jl`. It may make the existing
deployment note and compact current status truthful, but it may not add a
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

Acceptance must simulate, without creating or pushing a tag, pull-request,
`main`, `v0.2.0-rc1`, final `v0.2.0`, and malformed-tag contexts. It must:

- inspect the exact deployment classification for each context;
- build locally under simulated `main` and `v0.2.0-rc1` environments and
  inspect the generated canonical links;
- demonstrate under Documenter's standard version expansion that the
  prerelease does not select `stable` while a final release may;
- retain the current job- and step-level credential boundaries;
- pass authority check/self-test, generated-view parity, docs tests,
  manager-log bound, package load, Documenter, YAML/workflow inspection, and
  `git diff --check`;
- pass remote Docs and CI after the implementation push and confirm that the
  ordinary `main` deployment still updates the live `/dev/` site.

The implementation commit returns these two records for a separate docs-only
lifecycle closeout. Versioning, changelog and citation work, candidate tagging,
GitHub release creation, and final-release authorization remain separate
decisions.

Commit `62a1a4821` restored this standard same-repository Documenter deployment
path. Docs run `32072728238` and general CI run `32072728319` passed; the
`gh-pages` branch advanced from `a9b74566e` to `255ea4ed4`, built from that
commit. The live producer current-status page returned successfully with the
corresponding current-main content.

No custom credential, repository-setting change, alternate deployment
framework, source/API change, or documentation-navigation redesign is part of
this contract. GitHub's repository browser remains a raw-file view; the
published site provides the intended rendered navigation.
