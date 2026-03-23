# Export layer

The package’s current export layer is producer-side only. It now exposes both:

- public in-memory payload builders
- JLD2 writers that delegate to those payload builders

The current atomic IDA Hamiltonian model is exported honestly as dense or
sliced/block density-density data for downstream solver consumers.

For the narrative explanation of the current producer-side story, see:

- [Current atomic branch](../explanations/current_atomic_branch.md)
- [Example guide](../howto/example_guide.md)

```@docs
fullida_dense_payload
sliced_ham_payload
write_fullida_dense_jld2
write_sliced_ham_jld2
```
