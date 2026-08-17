# Contracted-Parent And Multilayer Retirement

## Status And Authority

This page records the completed retirement of:

- `HP-RETIRE-CONTRACTED-PARENT-FN-01`;
- `HP-RETIRE-CONTRACTED-PARENT-TEST-01`.

Status: source retirement completed at `7fbe512dc` and was accepted in Pass
467. The function ID is retired and the validation ID is completed; neither
grants further work. Current PGDG parent, source-box-first PQS contraction,
White--Lindsey realization, terminal basis, PRF mechanics, pair
materialization, Gaussian raw blocks, represented Hartree, and IDA assembly
were unchanged by this retirement. The pair ladder now has its own later
retirement contract.

## Decision

The contracted-parent and multilayer cluster was an abandoned alternative
representation, planning, and oracle route. No current numerical consumer
reached it. Pass 467 accepted its deletion rather than retaining internal
APIs, metadata plans, or test pressure through compatibility code.

The approved reduction is exact:

```text
source                 9,541 lines
obsolete test             66 lines
tracked probes            983 lines
total                  10,590 lines
```

The deletion added no source, test, probe, alias, stub, deprecation, adapter,
replacement helper, or moved kernel.

## Retired Source

Commit `7fbe512dc` deleted these files in full:

```text
src/CartesianContractedParents.jl
src/CartesianContractedParentMetrics.jl
src/cartesian_contracted_parent_metrics/core.jl
src/cartesian_contracted_parent_metrics/product_staged_metric_fallbacks.jl
src/cartesian_contracted_parent_metrics/source_box_pair_shadow.jl
src/pqs_multilayer_shell_region_plan.jl
src/pqs_multilayer_support_one_body.jl
src/pqs_multilayer_support_density.jl
src/pqs_multilayer_shell_source_plan.jl
src/pqs_multilayer_complete_core_shell_h1.jl
src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl
```

It removed only the associated seven root includes from
`src/GaussletBases.jl`. In
`src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`,
it removed the complete-core-shell include, its adjacent include comment, and the
stale module-doc claim that this utility serves live migration paths.

It deleted only these orphaned partial sections from otherwise retained files:

- `src/cartesian_gto_probes.jl:58-210`, the private source-box GTO shadow;
- `src/pqs_source_box_low_order_materialization.jl:1-131`, the dead atomic
  multilayer adapter.

The remaining contents and current callers of both partial files were
preserved.

## Retired Test And Probes

The pass deleted `test/ida/runtests.jl:3-68`, the historical PGDG-factor test
that called
the retired metrics submodule. Do not replace it: current parent/terminal IDA
behavior is protected by the unchanged active test groups and public
Cartesian endpoint.

It also deleted these tracked probes:

```text
tmp/work/benchmark_cartesian_contracted_parent_metric_packet.jl
tmp/work/validate_product_doside_source_box_unification.jl
tmp/work/validate_route_shaped_safe_term_consumer.jl
```

Other ignored probes may remain as historical evidence. They are not
compatibility obligations and must not receive adapters.

## Caller Proof

At baseline `6c9cf6032d0a062051109f796d7d3e4734fdcc09`:

- `CartesianContractedParents` is consumed only by
  `CartesianContractedParentMetrics`;
- the metrics module is reached only by the private GTO shadow, the isolated
  multilayer family, and the obsolete IDA testset;
- the apparent one-center atomic adapter has no caller;
- the complete-core-shell final-basis utility is called only by the
  multilayer H1 owner;
- no `bin`, tool, example, public endpoint, artifact, reflection,
  registration, serialization, `__init__`, or load-time consumer exists.

The two submodule APIs and the multilayer function family are internal
experimental vocabulary, not compatibility contracts.

No live committed caller or load-time side effect was found. Reintroduction
of any retired surface or compatibility glue requires new docs-only
authority.

## Authority Reconciliation

`HP-PQS-ASPECTSHELL-FN-01` remains implemented and in maintenance. Its active
owners are the base Hamiltonian, White--Lindsey terminal realizer, and
source-box route helper. The retired multilayer region/source plans no longer
own or carry the matched aspect-shell policy.

`CartesianContractedParents` is not an artifact or representation donor.
Future representation and provenance work starts from the current basis,
overlap, transfer, parent, terminal, and artifact owners.

## Preserved Owners

The accepted deletion did not delete, move, or reinterpret:

- `src/CartesianParentGaussletBases.jl` and
  `src/CartesianParentAxisFactors.jl`;
- current shellification and terminal PQS/White--Lindsey realization;
- the mapped-COMX axis-transform helper under
  `src/cartesian_pair_block_materialization/`;
- Gaussian analytic and raw-block owners;
- source-shell realization, PRF, residual-GTO, represented-Hartree, and IDA
  facilities;
- current basis representation, cross-overlap, transfer, artifacts, drivers,
  or public endpoints.

Historical contracted-parent plans may retain a short retirement banner or
link. They grant no restoration authority.

The no-caller pair-planning/materialization ladder was preserved by Pass 467
and is separately approved for retirement under
`HP-RETIRE-PAIR-LADDER-FN-01`. That later decision does not change this pass's
accepted line accounting or validation.

## Accepted Validation

The source retirement pass:

1. proved all audited module, function, include, test, and tracked-probe
   symbols absent from committed current source;
2. deleted exactly `9,541` source, `66` test, and `983` tracked-probe lines,
   with no replacement lines;
3. loaded the package in `6.58` seconds and passed the unchanged `core`,
   `ida`, and `cartesian` groups at `12,720/12,720` in `90.83` seconds;
4. passed the public Cartesian endpoint at `232/232` with unchanged terminal
   due diligence;
5. passed the anti-bloat scan, `git diff --check`, and GitHub numerical CI.

The pass did not edit current numerical tests merely to preserve retired
names. Reintroduction of any deleted surface requires a new docs-only design
amendment.
