# Documentation Deployment Note

The current Documenter site structure is already close to the intended
KrylovKit-style package-docs layout:

- Home
- Manual
- Reference
- Developer Notes

The remaining problem is not page organization inside the docs source tree. The
remaining problem is delivery surface.

GitHub's repository browser will always present the `docs/` directory as a file
listing. It cannot display the rendered Documenter sidebar, landing pages, or
top-level navigation in the way a deployed docs site can.

So the correct next step is to deploy the rendered site and treat that as the
primary documentation surface.

That means:

- keep the local docs build working with `julia --project=docs docs/make.jl`
- add standard `deploydocs(...)` deployment in `docs/make.jl`
- update the docs CI workflow so pushes to `main` deploy the site
- serve the published site through the usual GitHub Pages `gh-pages` branch
  target
- point README and repository-side docs links to the rendered docs URL rather
  than expecting users to browse raw markdown through the GitHub tree

This is still a documentation-delivery pass, not a new content or navigation
pass. The intended left-panel experience comes from the deployed Documenter
site, not from GitHub's file browser.
