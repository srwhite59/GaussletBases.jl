# Curated API Reference First Pass

This note records the first API-reference pass on top of the new Documenter
site.

## 1. Current site structure

The site now has a standard Julia-package skeleton:

- tutorials
- how-to pages
- explanations
- a small reference section

That is enough structure to start adding API material without forcing the
whole flat note tree into a reference role.

## 2. Why the next step is curated API reference

Once a package has a working Documenter site, the next missing piece is
usually not more narrative text. It is a stable reference layer for the
exported entry points that users actually call.

This package has now reached that point:

- the workflow pages explain what to do
- the branch pages explain what the current scientific lines mean
- the reference section should now explain what the main exported functions and
  objects are

## 3. Why `@docs` is the right first tool

For the first API pass, curated `@docs` blocks are better than a giant
undifferentiated `@autodocs` dump.

The reasons are:

- the export surface is broad enough that a full dump would be noisy
- many exports are still internal-leaning or advanced
- the first reference section should match the package’s actual public story
- the package is still choosing which layers are mature enough to foreground

So the first API pass should present a deliberate slice of the export surface,
grouped by user task rather than by source-file accident.

## 4. Division of labor: docstrings versus narrative pages

The narrative pages and the docstrings should not compete with each other.

Narrative pages should answer questions like:

- where do I start?
- what workflow do I use?
- what is the current interpretation of the atomic or ordinary branch?

Docstrings should answer questions like:

- what is this function or object?
- what does it return?
- what are the important keywords?
- what is the intended interpretation of the current model?

That is why the first API docstrings should be short, plain-language, and
contract-oriented rather than essay-like.

## 5. First curated reference slice

The first reference section should stay small and focus on the main user entry
points:

- basis and mapping constructors
- quadrature and diagnostics
- main operator builders
- atomic IDA entry points
- ordinary-branch user entry points
- export writers

That is enough to make the docs site feel like a recognizably standard Julia
package docs site without pretending that every export already has finished
reference documentation.
