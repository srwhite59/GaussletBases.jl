# Manual

This is the main user-facing manual for `GaussletBases.jl`.

If you are using the package as a scientist, this is the section to read
first. The manual is intentionally small. It tells you:

- what the package can already do reliably
- which workflow is the recommended starting point
- where the mature radial/atomic workflow, the newer ordinary/cartesian
  workflow, and the advanced research line currently stand
- which pages are background material rather than first-read documents

## Who this manual is for

Read this section if you want to:

- build and use gausslet bases
- reproduce the current radial and atomic workflows
- understand whether the ordinary Cartesian branch is ready for your use case
- find the right example or workflow page without digging through note history

If you already know the object or function you want, jump to the
[Reference](../reference/index.md) instead.

## Recommended reading order

For a new reader, the shortest useful path is:

1. [First radial workflow](../tutorials/first_radial_workflow.md)
2. [Recommended atomic setup](../howto/recommended_atomic_setup.md)
3. [Example guide](../howto/example_guide.md)

That path gets you from basis construction to a real hydrogen calculation and
then to the current atomic examples.

## Branch-specific paths

If you want atom-centered radial and atomic work:

1. [First radial workflow](../tutorials/first_radial_workflow.md)
2. [Recommended atomic setup](../howto/recommended_atomic_setup.md)
3. [Current atomic branch](../explanations/current_atomic_branch.md)

If you want the ordinary Cartesian branch:

1. [Current ordinary branch](../explanations/current_ordinary_branch.md)
2. [Example guide](../howto/example_guide.md)

The ordinary branch is worth reading after you understand the radial line. It
is the right place to learn the current mapped and hybrid ordinary workflows,
but it is not the best first entry point into the package.

## If you want more depth later

After the main workflow pages are clear, the next useful documents are:

- [Examples](../examples/index.md)
- [Reference](../reference/index.md)
- [Developer Notes](../developer/index.md)

The Developer Notes section is where the design history, architecture
background, and narrower supporting notes live. It is there when you want more
depth, without crowding the main manual.

## Pages in this manual

- [First radial workflow](../tutorials/first_radial_workflow.md)
- [Recommended atomic setup](../howto/recommended_atomic_setup.md)
- [Current atomic branch](../explanations/current_atomic_branch.md)
- [Current ordinary branch](../explanations/current_ordinary_branch.md)
