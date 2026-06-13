Accepted pass 178.

The audit did the right thing: it stopped Be2/Cr2 artifact work and re-centered
the loop on an atom-first PQS physics target. The important correction is now
explicit:

- The compact PQS 223 path is a seam/readiness route, not the physics target.
- The Be2 221 path is a midpoint/contact-slab source-box spine, not a sensible
  first Be2 WL/PQS comparison target.
- The old high-value WL reference was the Fig.8-style one-center He route with
  fixed local core and multiple fixed-width shells.
- Current route-owned PQS multi-layer source planning exists, but its source
  realization grows shell source dimensions from each layer's inner box. That
  gives the old q-ladder/growing-q behavior, not the fixed-q shell inventory we
  need to match WL.

The specific implementation seam is narrower than the previous confusion made
it look. `PQSMultilayerShellLayerRegion` already carries
`source_mode_shape`, and terminal lowering records it on PQS lowering
contracts. The source-plan realization currently ignores that fact and uses
`length.(inner_box)` to select `q`/`raw_source_dims`.

Next pass should make source realization consume the route-owned fixed source
shape for shellification/lowering-backed PQS layers, while preserving the
explicit-box bridge's old growing-q behavior. The first target is a focused
one-center He q=5/n_s=5 fixed-shell inventory: core 125 plus three shell
sectors of 98 retained modes each, for gausslet-only retained dimension 419.

No Be2/Cr2 handoff work should resume until this atom path is clarified.

-- repo-manager@macmini
