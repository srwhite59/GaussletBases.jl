# High-Order Doside Terminology

This note locks in the current terminology for the two main high-order
three-dimensional basis constructions.

## Standard terms

### Full-block union (FBU)

The **full-block union (FBU)** is the basis obtained by taking the union of the
full transformed local blocks over the selected side ladder, followed by the
usual metric cleanup / orthonormalization.

This is the standard name for what was previously described in various places
as:

- full-cube union
- filled-cube union
- full-block union

Use **full-block union (FBU)** as the preferred name.

### Full-shell basis (FSB)

The **full-shell basis (FSB)** is the basis obtained from:

- the central full transformed block at the starting side, and
- the outer-shell functions taken from larger full transformed blocks,
- with the usual projection against the accumulated inner span and final
  cleanup / orthonormalization.

This is the standard name for what was previously described in various places
as:

- structured shell basis
- doside stack
- shell-stacked basis

Use **full-shell basis (FSB)** as the preferred name.

## Reporting convention

In summaries and notes:

- write the full name on first use
- use the acronym after that

Examples:

- “The **full-block union (FBU)** closes the span gap.”
- “The **full-shell basis (FSB)** is cheaper but may miss some directions.”
- “Later in the same summary, refer to them simply as **FBU** and **FSB**.”

## Scope of the names

The name **full-shell basis (FSB)** should remain in place even if the internal
construction is improved, for example by changing the final orthogonalization
or shell cleanup details, as long as the construction is still fundamentally:

- full transformed block at the center
- shell-only additions from larger full transformed blocks
- hierarchical projection / cleanup

So the names are meant to track the conceptual objects, not every internal
implementation detail.
