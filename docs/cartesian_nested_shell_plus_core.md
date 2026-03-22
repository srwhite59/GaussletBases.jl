## Cartesian Nested Shell-Plus-Core Fixed Block

The shell-only fixed-block adapter was the right plumbing experiment:

- the shell packet itself is structurally valid
- the fixed-Gaussian cross-block contraction pattern works
- and the nearest/GGT QW-PGDG consumer can already read a preassembled 3D fixed packet

But the shell-only scalar result also showed the expected physical limitation:
the shell packet by itself does not represent the interior/core region well
enough to replace the full parent fixed space.

So the next pre-recursion step is not more shell refinement. It is a
shell-plus-core fixed block.

The first core choice is intentionally simple:

- keep one interior rectangular block of parent fixed functions directly
- keep the current shell faces unchanged
- assemble one combined fixed coefficient matrix from
  - direct parent core columns
  - plus the propagated shell-face columns

This keeps the geometry logic clean:

- shell faces still cover disjoint face interiors
- the core fills the missing interior volume
- the combined fixed block is still read by the same nested fixed-block adapter

The first consumer target stays narrow:

- QW-PGDG
- `interaction_treatment = :ggt_nearest`
- no MWG
- no recursion

If this restores a physically meaningful fixed block, the next remaining
pre-recursion step is to decide how the current shell-plus-core object should
grow into a sequence of shell layers, not whether the nested consumer adapter
itself is viable.
