# Public Story Status Note

This is a small documentation-story/status pass after the public opening of the
repository and the deployment of the docs site.

The current public documentation already does several important things well:

- it introduces gausslets in plain scientific language
- it gives the radial/atomic line a clear first-path status
- it has a visible rendered docs site with Manual, Examples, Reference, and
  Developer Notes
- it labels the ordinary branch honestly as newer and still experimental

The main remaining imbalance is subtler:

- the ordinary/cartesian line needs to read as a real second workflow, not just
  as an experimental footnote
- the contraction/hierarchy/prototype PGDG material needs to stay clearly
  separated as the advanced/research line
- `STATUS.md` and `ROADMAP.md` should use the same three-layer public story as
  the front-door docs

This pass sharpens that distinction:

- mature radial / atomic workflow
- real newer ordinary / Cartesian workflow
- advanced/developer research line

I also audited README, docs pages, and runnable examples for generalized-eigen
teaching language. No user-facing teaching path was still presenting a
generalized eigensolve as the normal route. The public examples and tutorials
already use ordinary Hermitian eigensolves for the orthonormal teaching path,
with overlap matrices treated as diagnostics or orthogonalization ingredients
rather than as the main solve interface.
