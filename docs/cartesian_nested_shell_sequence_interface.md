## Cartesian Nested Shell-Sequence Interface

The shell-plus-core milestone settles the main pre-recursion question:

- the nested fixed-block adapter works
- a shell-only fixed block is not physically adequate
- a shell-plus-core fixed block already is

So the next step should stay in the same fixed-block language and generalize
only from one shell-plus-core block to a short shell sequence.

The intended interface is:

- one direct interior core block
- plus an ordered list of shell layers from inside to outside
- one combined fixed coefficient matrix
- one propagated fixed-block packet

A natural first internal object would be something like
`_CartesianNestedShellSequence3D`, carrying:

- `core_indices`
- `shell_layers::Vector{_CartesianNestedShell3D}`
- `layer_column_ranges`
- `coefficient_matrix`
- `support_indices`
- `support_states`
- `packet`

The consumer model should not change again. The existing nested fixed-block
adapter should still read:

- fixed-fixed data directly from the combined packet
- fixed-Gaussian raw blocks by contraction through the combined coefficient
  matrix
- Gaussian-Gaussian blocks on the existing analytic route

So the shell-sequence step is only a richer source of the same fixed-block
object, not a new consumer path.

### First Nonrecursive Multi-Shell Test Case

The first practical nonrecursive sequence should stay on the stabilized
fixed-`a` He family and just add one more shell layer.

A good target is:

- He, fixed-`a` family with `a = 1/4`, `xmax = 10`
- `count = 17`
- direct core box on `4:14`
- first shell around that core:
  - tangential intervals `4:14`
  - fixed indices `(3, 15)`
- second shell outside that:
  - tangential intervals `3:15`
  - fixed indices `(2, 16)`

This is still nonrecursive in implementation:

- construct two shell layers explicitly
- concatenate them with the direct core block
- propagate one combined packet

The first validation should ask only:

- overlap quality of the combined sequence packet
- disjointness of the shell interiors
- end-to-end nearest/GGT QW-PGDG behavior
- and whether the second shell materially changes the scalar relative to the
  one-shell-plus-core baseline
