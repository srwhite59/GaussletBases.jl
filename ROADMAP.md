# Roadmap

This is a maintenance and research boundary, not a schedule or a commitment
to additional capabilities. For usage, start with the README and manuals;
STATUS.md distinguishes supported interfaces from experimental work.

## Current maintenance

- Preserve radial/atomic construction as the beginner workflow.
- Maintain PQS as the standard Cartesian base construction, with
  White-Lindsey as a matched comparison and alternative.
- Maintain the bounded supplied-field Hartree screening and external
  Cartesian GTO transfer interfaces as separate capabilities.
- Keep basis, operator, and transfer conventions explicit, with validation
  tied to supported public examples and downstream contracts.

## Research boundaries

Supported construction does not establish convergence for arbitrary systems.
General molecular geometries, broad SCF/post-HF workflows, and native DMRG
remain outside the package's public workflow contract. Angular injection,
hierarchy/contraction, and experimental chain/square nested producers retain
their narrower research status.

Future extensions require a named scientific endpoint and independent review;
this roadmap does not authorize promotion, new solvers, or compatibility
changes. Historical prototypes and private diagnostics are not new public
commitments.
