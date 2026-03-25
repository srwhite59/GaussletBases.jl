# Manual

This is the main user-facing manual for `GaussletBases.jl`.

If you are using the package as a scientist, this is the section to read
first. The manual is intentionally small. It tells you:

- what the package can already do reliably
- which workflow is the recommended starting point
- where the mature radial/atomic workflow, the newer ordinary/cartesian
  workflow, and the advanced research line currently stand
- which pages are background material rather than first-read documents

Detailed basis-construction and operator-construction recipes live in the
[Algorithms](../algorithms/index.md) section. The Manual tells you what to
read and run; Algorithms records how the main constructions are built.

The scientific motivation throughout is the gausslet one from the papers:
localized, orthonormal Gaussian-built basis functions with systematic spacing
refinement and a two-index or diagonal Coulomb structure. The manual is
organized around the workflows that currently realize that idea most cleanly.

## Who this manual is for

Read this section if you want to:

- build and use gausslet bases
- reproduce the current radial and atomic workflows
- understand whether the ordinary Cartesian branch is ready for your use case
- find the right example or workflow page without digging through note history

If you already know the object or function you want, jump to the
[Reference](../reference/index.md). If you want the exact construction order
for a basis or operator path, jump to [Algorithms](../algorithms/index.md).

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

If you want the active experimental angular line:

1. [Angular research track](../explanations/angular_research_track.md)

That page is intentionally narrow. It is there to mark the current research
boundary, not to present a finished angular workflow.

## If you want more depth later

After the main workflow pages are clear, the next useful documents are:

- [Algorithms](../algorithms/index.md)
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
- [Angular research track](../explanations/angular_research_track.md)
