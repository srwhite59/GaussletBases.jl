# Complete-Core-Shell RHF Retirement

Status: approved deletion authority under `HP-RETIRE-CCS-RHF-FN-01` and
`HP-RETIRE-CCS-RHF-TEST-01`.

## Decision

The old complete-core-shell RHF payload stack is approved for retirement. The
current Cartesian Hamiltonian producer workflow is the canonical staged driver
and `CartesianIDAHamiltonian` artifact path. CR2-facing atom/diatomic,
base/supplemented, and PQS/WL use now runs through
`bin/cartesian_ham_builder.jl`, staged producer functions, and Hamiltonian
artifacts, not through the old complete-core-shell RHF payload stack.

The deletion target is:

```text
src/pqs_multilayer_complete_core_shell_rhf.jl
```

and its root include in:

```text
src/GaussletBases.jl
```

## Evidence

The source file is currently included by `src/GaussletBases.jl`, but focused
search found no live `src`, `bin`, `test`, or `tool` caller outside the file
itself. Remaining references are the include, the file internals, and
docs/history material.

The file carries stale route-era payload/status vocabulary, including:

- `pqs_multilayer_complete_core_shell_rhf_input_contract`;
- `pqs_multilayer_complete_core_shell_rhf_scf_payload`;
- `pqs_multilayer_complete_core_shell_rhf_one_step_payload`;
- blocked payload/status helpers.

It depends on the older complete-core-shell H1/final-basis workflow rather than
the current canonical driver, base/supplemented producer, and
`CartesianIDAHamiltonian` artifact path.

## Approved IDs

- `HP-RETIRE-CCS-RHF-FN-01` - remove the stale RHF stack and include.
- `HP-RETIRE-CCS-RHF-TEST-01` - validate active producer paths after deletion.

## Approved Source Surface

Approved later source files:

```text
src/GaussletBases.jl
src/pqs_multilayer_complete_core_shell_rhf.jl
```

Approved behavior:

- remove the `pqs_multilayer_complete_core_shell_rhf.jl` include from
  `src/GaussletBases.jl`;
- delete `src/pqs_multilayer_complete_core_shell_rhf.jl`;
- remove only docs/index references that describe this RHF stack as active
  current code, if encountered during the deletion pass;
- do not add replacements, adapters, compatibility wrappers, status objects,
  or tests.

Expected line result:

- net deletion of roughly `1879` source lines;
- only a small include deletion and any minimal stale active-reference cleanup.

## Forbidden

This retirement lane does not approve:

- canonical driver changes;
- source changes outside the approved file/include except minimal stale
  active-reference cleanup if directly required;
- changes to `pqs_multilayer_complete_core_shell_h1.jl`;
- changes to `pqs_complete_core_shell_final_basis.jl`;
- changes to `pqs_source_box_low_order_materialization.jl`;
- ordinary or Qiu-White donor-kernel changes;
- artifact schema, provenance, reader, or manifest changes;
- route, shellification, terminal-lowering, raw-block, Residual Gaussian,
  MWG, IDA, or Hamiltonian assembly changes;
- replacements, adapters, compatibility wrappers, reports, status fields, or
  payload objects;
- committed tests or fixtures;
- Cr2 workflow.

## Validation

`HP-RETIRE-CCS-RHF-TEST-01` approves only:

- `git diff --check`;
- package load;
- focused `rg` showing no remaining live references to
  `pqs_multilayer_complete_core_shell_rhf`;
- canonical small base artifact/readback smoke;
- canonical small supplemented artifact/readback smoke;
- unchanged H2 Residual Gaussian endpoint;
- no Cr2 run.

## Failure Rule

If any live `src`, `bin`, `test`, or `tool` caller depends on the RHF stack,
make no source commit and report the exact caller. Do not preserve the path
through an adapter.
