# Documentation Structure Plan

This note records the purpose of the current documentation-structure pass.

## 1. What the top-level layers are right now

The repository already has strong scientific notes, but they are spread across
several different layers:

- the `README.md` landing page
- `docs/example_guide.md`
- current workflow notes for the radial, atomic, and ordinary branches
- narrower research and development notes that record why later design choices
  were made

Each of those layers is useful on its own.

## 2. Where reader overload is happening

The current overload is mostly architectural rather than scientific.

Three things have drifted together:

- onboarding docs for a new reader
- current workflow docs for the main lines of the package
- supporting notes that explain how the current state was reached

As a result, a new reader has to infer too much:

- which pages are start-here pages
- which pages describe the current recommended interpretation
- which pages are narrower supporting notes

## 3. What this pass should do

This pass should separate those layers more clearly without rewriting the
science:

- keep `README.md` as a front door
- add one small docs map page
- reshape the example guide into an actual reading/running guide
- add one short current-status page for the atomic branch
- add one short current-status page for the ordinary branch
- leave the narrower historical and development notes in place, but make them
  easier to recognize as supporting material

## 4. What this pass should not do

This is not a scientific redesign pass.

It should not:

- change the public math story
- archive large piles of notes
- do a large style-only rewrite
- blur the current wording discipline around the radial numerical path, the
  ordinary experimental PGDG path, or the provisional status of the current
  mapping heuristics

## 5. Intended result

After this pass, the intended reading structure should be:

1. `README.md` for a quick landing page
2. `docs/index.md` for the current docs map
3. `docs/example_guide.md` for example sequencing
4. `docs/current_atomic_branch.md` or `docs/current_ordinary_branch.md`
   depending on branch interest
5. narrower supporting notes only after the current branch-status pages are
   clear
