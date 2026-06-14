Review 225 - accepted

Pass 225 completed the no-edit audit of the H2 PQS MWG/GTO supplement provider
seam.

Accepted findings:

```text
Provider-local reusable blocks:
  CartesianCPBBlockProviders
    cpb_mixed_gto_overlap_block
    cpb_gto_overlap_operator_block
    cpb_gto_nuclear_by_center_block
    cpb_gto_supplement_local_operator_bundle

Route-global combined-GTO machinery:
  route_global_combined_gto_basis_layout
  route_global_mixed_gto_blocks_from_decomposed_units
  route_global_combined_gto_one_electron_matrices
  route_global_combined_gto_final_basis_projection
  route_global_residual_gto_mwg_representation
  route_global_combined_gto_final_basis_density_density_matrix
```

The important architectural conclusion is correct: the provider kernels and
combined-basis algebra are reusable in principle, but the current route-global
mixed-block adapter is WL/decomposed-unit shaped. PQS needs a route-owned
request/representation seam first, not immediate provider-block construction.

Current PQS H2 physical facts already available:

```text
centers: (0,0,-2), (0,0,2)
nuclear charges: (1,1)
parent axis counts: (9,9,15)
physical support counts: (275,578,362)
retained counts: (251,98,114)
gausslet final dimension: 463
source/final basis payloads and no-supplement H1/H1-J/RHF diagnostics
```

The next implementation should therefore add a compact supplement request
payload before provider matrices. It should bind the policy, geometry,
supplement basis identity, residual policy, and representation status. If it
does not create a real supplement representation yet, the blocker should sharpen
to the representation-level missing fact; if it does, the blocker can remain
`:missing_provider_gto_supplement_blocks`.

No Julia validation was needed for this read-only pass. The tree remained clean
and even with origin.

Next direction:

Add `_PQSDiatomicPhysicalGaussletSupplementRequestPayload`, keep it
matrix-free, write compact artifact fields, and keep source/test/bin
net-negative by shrinking repeated supplement-preflight assertions rather than
adding another field-cloud test.

-- repo-manager@macmini
