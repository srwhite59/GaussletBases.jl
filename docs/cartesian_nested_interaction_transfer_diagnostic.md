# Cartesian Nested Interaction-transfer Diagnostic

This note records the first interaction-side diagnostic after the
coverage-fixed complete shell sequence.

The current starting point is:

- shell-plus-core is the last trusted physical nonrecursive state
- the corrected complete four-shell sequence now captures the parent low-energy
  one-body content very well
- `E1` is back in the right regime
- but `⟨Vee⟩` is still far too large

So the next question is whether the remaining error comes from:

- the fixed-fixed interaction packet itself
- or the later residual-sector nearest/GGT interaction extension

## Diagnostic Definition

The comparison was done on the stabilized He count-17 fixed-`a` case used in
the earlier nesting notes.

Two closely related diagnostics were run.

### 1. Fixed-only projected parent-ground test

Build the parent fixed-block one-body ground orbital `u_parent` and evaluate
the fixed-only interaction expectation on:

- the parent fixed space
- the shell-plus-core retained fixed space
- the corrected complete four-shell retained fixed space

For one retained fixed block with interaction matrix `V_fixed`, project the
parent ground orbital into that fixed space and evaluate

```math
\langle V_{ee} \rangle_{\mathrm{fixed}} =
w^T V_{\mathrm{fixed}} w,
\qquad
w_i = \frac{|c_i|^2}{\sum_j |c_j|^2}.
```

This removes the residual-sector nearest/GGT extension entirely.

### 2. Final nearest/GGT sector decomposition

For the actual final `1s^2` orbital from
`ordinary_cartesian_qiu_white_operators(...; interaction_treatment = :ggt_nearest)`,
split the final expectation into:

- fixed-fixed
- fixed-residual
- residual-residual

using the final orbital probability weights.

This tests whether the remaining excess `⟨Vee⟩` is mainly coming from the
residual extension or is already present in the fixed-fixed packet.

The scratch driver is
[nested_interaction_transfer_diagnostic.jl](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/tmp/work/nested_interaction_transfer_diagnostic.jl).

## Test Case

Same stabilized He case:

- `count = 17`
- `a = 1/4`
- `xmax = 10`
- `s = 0.626026121152214`
- base QW-PGDG nearest/GGT path

## Fixed-only Projected Parent-ground Result

Parent fixed-space one-body ground orbital:

- energy `-1.9984765858974547`
- parent fixed-space interaction expectation `1.2486322227019246`

Projected into shell-plus-core fixed block:

- fixed dim `1403`
- projected fixed-space interaction expectation `1.2486612176560734`
- shift from parent `+2.8994954148764407e-5`

Projected into corrected complete four-shell fixed block:

- fixed dim `589`
- projected fixed-space interaction expectation `1.818451883926288`
- shift from parent `+0.5698196612243633`

This is the key number.

Even after removing the residual-sector nearest/GGT extension completely, the
corrected complete shell sequence already gives the wrong interaction
expectation on the projected parent ground orbital.

So the remaining `⟨Vee⟩` failure is not mainly a residual-extension artifact.

## Final nearest/GGT Sector Decomposition

Shell-plus-core final `1s^2`:

- `nfixed = 1403`
- `nresidual = 2`
- orbital energy `-1.9985629071304167`
- total `⟨Vee⟩ = 1.2489197930296585`
- fixed weight mass `0.9999836623020556`
- residual weight mass `1.633769794439465e-5`
- fixed-fixed contribution `1.2488545499498696`
- fixed-residual contribution `6.52387430635677e-5`
- residual-residual contribution `4.3367250344410704e-9`

Corrected complete four-shell final `1s^2`:

- `nfixed = 589`
- `nresidual = 2`
- orbital energy `-1.9981842264017424`
- total `⟨Vee⟩ = 1.816303754935915`
- fixed weight mass `0.9999903655915707`
- residual weight mass `9.6344084292005e-6`
- fixed-fixed contribution `1.8162601086253625`
- fixed-residual contribution `4.364480245189598e-5`
- residual-residual contribution `1.5081004472931866e-9`

So in both cases the final orbital is almost entirely in the fixed sector, and
for the corrected complete sequence almost all of the excessive `⟨Vee⟩` is
already in the fixed-fixed block.

## Interpretation

This gives a concrete interaction-side diagnosis.

What is ruled out:

- the residual nearest/GGT extension is not the main problem
- residual-sector couplings are numerically tiny here

What remains:

- the corrected complete shell sequence is still transferring the interaction
  packet poorly inside the compressed fixed block itself
- this happens even though the same retained space captures the parent
  low-energy one-body orbitals very well

So the failure is now more specific:

- not a general one-body retained-span failure
- not a residual-extension failure
- but a fixed-fixed interaction representation failure

The most likely structural reading is:

- the current shell contraction language preserves the low-energy amplitude
  space much better than it preserves the induced diagonal-density / pair
  content needed by the two-index IDA interaction

In other words:

- one-body orbital capture is already good enough
- but the contracted shell packet is still not preserving the interaction
  functional seen by `|psi|^2`

## Conclusion

This is strong enough to guide the next nonrecursive step.

The next nesting diagnostic should stay on the interaction side and target the
fixed-fixed pair representation directly, for example:

- projected parent-density interaction transfer
- low-rank pair-density capture on the shell pieces
- or a contraction policy that preserves both amplitude and diagonal-density
  content

The roadmap change from this pass is:

- do not spend the next pass on residual nearest/GGT cleanup first
- do not immediately retune shell counts first
- focus next on the fixed-block interaction / pair-content representation
  inside the complete nonrecursive shell language
