# Contracted-Parent And Multilayer Retirement

## Status And Authority

This page owns the approved retirement contract for:

- `HP-RETIRE-CONTRACTED-PARENT-FN-01`;
- `HP-RETIRE-CONTRACTED-PARENT-TEST-01`.

The records authorize one exact deletion pass. They do not retire the current
PGDG parent, source-box-first PQS contraction, White--Lindsey realization,
terminal basis, PRF mechanics, pair materialization, Gaussian raw blocks,
represented Hartree, or IDA assembly.

## Decision

The contracted-parent and multilayer cluster is an abandoned alternative
representation, planning, and oracle route. No current numerical consumer
reaches it. Delete it rather than retain its internal APIs, metadata plans,
or test pressure through compatibility code.

The approved reduction is exact:

```text
source                 9,541 lines
obsolete test             66 lines
tracked probes            983 lines
total                  10,590 lines
```

Add no source, test, probe, alias, stub, deprecation, adapter, replacement
helper, or moved kernel.

## Source Deletion

Delete these files in full:

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

Remove only the associated seven root includes from `src/GaussletBases.jl`.
In
`src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`,
remove the complete-core-shell include, its adjacent include comment, and the
stale module-doc claim that this utility serves live migration paths.

Delete only these orphaned partial sections from otherwise retained files:

- `src/cartesian_gto_probes.jl:58-210`, the private source-box GTO shadow;
- `src/pqs_source_box_low_order_materialization.jl:1-131`, the dead atomic
  multilayer adapter.

Preserve the remaining contents and current callers of both partial files.

## Test And Probe Deletion

Delete `test/ida/runtests.jl:3-68`, the historical PGDG-factor test that calls
the retired metrics submodule. Do not replace it: current parent/terminal IDA
behavior is protected by the unchanged active test groups and public
Cartesian endpoint.

Delete these tracked probes:

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

If a live committed caller or load-time side effect is found before deletion,
stop without committing and report it. Do not narrow the retirement by adding
glue.

## Authority Reconciliation

`HP-PQS-ASPECTSHELL-FN-01` remains implemented and in maintenance. Its active
owners are the base Hamiltonian, White--Lindsey terminal realizer, and
source-box route helper. The retired multilayer region/source plans no longer
own or carry the matched aspect-shell policy.

`CartesianContractedParents` is not an artifact or representation donor.
Future representation and provenance work starts from the current basis,
overlap, transfer, parent, terminal, and artifact owners.

## Preserved Owners

Do not delete, move, or reinterpret:

- `src/CartesianParentGaussletBases.jl` and
  `src/CartesianParentAxisFactors.jl`;
- current shellification and terminal PQS/White--Lindsey realization;
- `src/cartesian_pair_block_materialization/`;
- Gaussian analytic and raw-block owners;
- source-shell realization, PRF, residual-GTO, represented-Hartree, and IDA
  facilities;
- current basis representation, cross-overlap, transfer, artifacts, drivers,
  or public endpoints.

The historical contracted-parent plans may retain a short retirement banner
or link. They grant no restoration authority.

## Acceptance

The source retirement pass must:

1. prove all audited module, function, include, test, and tracked-probe symbols
   are absent from committed current source;
2. confirm exactly `9,541` source, `66` test, and `983` tracked-probe lines are
   deleted, with no replacement lines;
3. pass package load and the unchanged `core`, `ida`, and `cartesian` groups;
4. inspect the public Cartesian endpoint's terminal due-diligence report;
5. pass the anti-bloat scan and `git diff --check`.

The pass must not edit current numerical tests merely to preserve retired
names. Reintroduction of any deleted surface requires a new docs-only design
amendment.
