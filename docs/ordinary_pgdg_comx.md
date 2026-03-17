# Ordinary Mapped PGDG: Minimal Cleanup and COMX Follow-On

This note describes the next narrow experimental step after the first analytic
primitive/contraction-level PGDG-style prototype.

## 1. What the current pre-COMX prototype is doing

The present `MappedPGDGPrototype1D` object does one specific thing:

- it replaces the explicitly distorted primitive layer behind a
  `MappedUniformBasis` by a locally matched plain-Gaussian primitive proxy
- it keeps the current working basis structure through the existing contraction
  matrix
- it then evaluates one-body matrices analytically on that plain-Gaussian
  primitive layer

So the current prototype is a test of the **analytic primitive/operator route**,
not yet a recreation of the full historical PGDG basis-construction path.

## 2. What is still missing before the path is historically faithful

The missing pieces are the basis-construction steps that come after the
primitive layer has been chosen.

At minimum, a historically faithful follow-on needs:

- overlap cleanup / orthogonalization of the working space
- a COMX-style position-localization step
- consistent ordering and sign conventions after localization

Those steps matter because they change the representation of the working basis,
even when the primitive layer and the one-body operators are already fixed.

## 3. Where COMX enters

In the present package language, COMX belongs to the basis layer, not to the
primitive integral layer.

The sequence is:

1. build primitive-space matrices analytically
2. contract them to the current working basis
3. clean up the basis overlap
4. diagonalize the position operator in that cleaned-up space
5. order and phase-fix the localized functions

So COMX is not a different integral formula. It is a controlled change of basis
inside the working space.

## 4. What this first minimal follow-on includes

The present follow-on stays narrow:

- projected overlap matrix
- overlap cleanup / orthogonalization
- COMX-style localization through the one-dimensional position operator
- deterministic ordering by localized center
- a simple sign convention

It still does **not** attempt to port the larger historical PGDG driver logic.

## 5. How comparisons should be read

Three levels are useful now:

1. current mapped numerical ordinary path  
   same current working basis, numerical mapped one-body matrices

2. current pre-COMX analytic prototype  
   same current working basis, analytic plain-Gaussian primitive proxy

3. post-COMX analytic path  
   analytic plain-Gaussian primitive proxy plus the first cleanup/localization
   step

The first comparison tests whether the analytic primitive route is promising at
all.

The third path then tests a different question:

**does the first faithful cleanup/localization step materially improve the
analytic path, or does the remaining discrepancy mainly come from the primitive
proxy itself?**
