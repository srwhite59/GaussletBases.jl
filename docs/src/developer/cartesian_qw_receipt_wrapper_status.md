# Historical Cartesian QW Receipt Wrapper Status

Status: historical evidence. The receipt implementation in
`src/cartesian_qw_operator_carried_spaces.jl` is approved for deletion under
`HP-RETIRE-QW-DONOR-*` and grants no current source authority.

## What The Wrapper Was

The receipt line wrapped existing QW builders without owning Hamiltonian
kernels. It normalized a carried space, recorded route/backend/storage facts,
delegated to `ordinary_cartesian_qiu_white_operators(...)`, and returned an
audit record around the unchanged operator payload.

Its main public/internal objects were:

- `CartesianOperatorBuildSource3D`;
- `CartesianQWOperatorConstructionRecord3D`;
- `CartesianQWOperatorConstructionReceipt3D`;
- carried-space sidecars and diagnostics.

The only committed consumer was the experimental chain/square implementation
in `ordinary_qw_experimental_paths.jl`. Current PQS/WL producer composition,
artifacts, and consumers do not use the receipt layer.

## Retirement Interpretation

The wrapper proved that pre-build request facts could be compared with an
existing QW operator payload, but its 24-name submodule surface never became a
needed producer or consumer contract. Keeping it solely because it was exported
would preserve an untested experimental API with no caller.

Delete the receipt module together with its sole consumer. Do not create
stubs, deprecation wrappers, aliases, or a smaller compatibility receipt.

This decision does not delete or prejudge:

- `cartesian_qw_hybrid_representation.jl`;
- `cartesian_carried_spaces.jl`, which requires a separate caller/API audit;
- active QW builders and residual/raw/operator kernels;
- chain/square basis types and geometry diagnostics.

See
[QW and high-order experimental cluster retirement](designs/cartesian_hamiltonian_producer/qw_high_order_experimental_retirement.md)
for the complete deletion and validation contract.
