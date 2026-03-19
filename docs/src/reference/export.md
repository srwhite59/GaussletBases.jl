# Export layer

The package’s current export layer is producer-side only. These writers export
the current atomic IDA Hamiltonian model honestly as dense or sliced/block
density-density data for downstream solver consumers.

For the narrative explanation of the current producer-side story, see:

- [Current atomic branch](../explanations/current_atomic_branch.md)
- [Example guide](../howto/example_guide.md)

```@docs
write_fullida_dense_jld2
write_sliced_ham_jld2
```
