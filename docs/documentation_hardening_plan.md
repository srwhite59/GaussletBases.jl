# Documentation Hardening Plan

This note records the first docs-hardening pass after the Documenter
transition and the first curated API-reference pass.

## 1. What is in place now

The package now has the main documentation structure it needed:

- a Documenter-based docs site
- curated tutorial, how-to, explanation, and reference sections
- curated API pages built from real docstrings

That is enough documentation structure for now. The next job is not another
large content pass. It is to make the docs build more robust and easier to
maintain.

## 2. What the hardening pass should do

The current priorities are:

- make the docs build run automatically in CI
- add a small doctest slice for the highest-value entry points
- improve cross-links between workflow pages and reference pages
- keep the API reference curated rather than turning it into a giant dump

This is a maintenance and discipline pass, not a structural rewrite.

## 3. Why this is not a `@autodocs` pass

The curated API pages are already the right first reference layer. The package
still has a broad export surface, and many exports are more advanced than the
current public story. A giant undifferentiated `@autodocs` dump would work
against the current teaching flow.

So the hardening pass should keep:

- curated `@docs` pages
- selected docstrings for the real entry points
- narrative pages for workflow and interpretation

## 4. `checkdocs` stance for this pass

This pass should revisit `checkdocs`, but not force an all-at-once strictness
jump. The curated API slice is real, but it still covers only the main public
entry points rather than every exported symbol.

So the right first move is:

- keep `checkdocs` conservative for now
- add CI docs builds
- add a small doctest slice
- keep expanding docstrings intentionally

Once a broader fraction of the export surface has real docstrings, tightening
`checkdocs` toward a more standard `:exports` stance will be much easier.

## 5. Success criterion

After this pass, the docs should:

- build automatically in CI
- have a few real doctested examples
- have better manual/reference cross-links
- feel stable enough to leave mostly alone unless the science changes
