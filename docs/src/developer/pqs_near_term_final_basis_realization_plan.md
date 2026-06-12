# PQS Near-Term Final-Basis Realization Plan

This note captures the immediate module-boundary, algorithm-shape, and test
discipline decisions after the 2026-06-12 PQS source-box pivot baton loop.
It is a near-term implementation guide, not a new framework layer.

## Current Status

PQS has crossed from metadata/source-box scaffolding into a real, oracle-backed
one-electron final-basis seam for the cubic projected-q-shell fixture:

```text
raw source box
-> PQS boundary source-mode retained rule
-> retained source-mode overlap / kinetic / nuclear blocks
-> shell-realization final basis
-> final overlap / kinetic / by-center nuclear
-> final one-electron Hamiltonian
-> ordinary final-basis H1 eigensolve
```

The final H1 probe showed that the final overlap is identity to roundoff, the
final Hamiltonian agrees with a shell-support oracle to roundoff, and the
ordinary final-basis H1 eigenvalue agrees with the oracle eigensolve to
roundoff.

This is not yet a production PQS route. The shell projection and Lowdin cleanup
inputs are still supplied from the shell-realization/oracle layer. The correct
status is:

```text
source-box retained operators are real;
the final-basis one-electron seam is validated against an oracle;
source-owned shell-realization construction is not yet fully route-owned;
IDA, density-density, RHF, driver adoption, exports, and artifacts are not done.
```

## Module Boundary Decision

The overnight work revealed a likely missing concept module. It is not a
general `PQS` module and not a replacement for `CartesianRawProductSources`.
The missing concept is final-basis realization and final-basis operator
transfer.

Current conceptual split:

```text
CartesianRawProductSources
    raw source CPB, source-mode dimensions/order, axis/source facts,
    and narrow PQS source-mode boundary selector metadata

CartesianPairBlockMaterialization
    raw source operator blocks, retained boundary/source operator blocks,
    direct/direct and White-Lindsey local/pair materialization seams

CartesianFinalBasisRealization       # proposed
    shell/support realization data, Lowdin cleanup, final overlap diagnostics,
    and retained-boundary/operator -> final-basis operator transfer
```

`CartesianRawProductSources` is a raw-source fact/bundle module. It may own the
narrow PQS source-mode boundary selector metadata because that selector is tied
directly to source-mode ordering and raw source facts. It should not own
general retained-rule policy, shell projection, Lowdin cleanup, final retained
units, IDA weights, pair blocks, Hamiltonians, exports, or artifacts.

`CartesianPairBlockMaterialization` should not keep accumulating later-stage
PQS final-basis concepts. The stable later-stage seam is:

```text
boundary retained source-space basis
-> shell/support realization
-> Lowdin final orthonormal basis
-> final-basis operator transform
```

That seam should move into a module such as:

```text
CartesianFinalBasisRealization
```

The first extraction should be mechanical. Do not combine it with result-type
redesign, IDA, density-density, or driver adoption.

Candidate functions to move first:

```text
pqs_source_shell_realization_final_basis
pqs_source_shell_final_one_body_from_boundary_matrix
pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block
pqs_source_shell_projected_one_body_matrix    # oracle/helper
```

`pqs_source_shell_final_one_electron_hamiltonian` may move only if the new
module's docstring explicitly includes narrow final one-body assembly from
already-final operators. It is acceptable to leave Hamiltonian assembly in
CPBM for the first extraction.

## Algorithm-Shape Optimization

The current retained PQS source-mode contraction materializes a full raw
source-space block and then selects boundary rows/columns:

```text
O_retained = O_source[left_columns, right_columns]
```

That is fine for the current `5 x 5 x 5` fixture, where the raw source count is
`125` and the retained boundary count is `98`. It should not become the
scalable algorithmic habit.

The next algorithm-shape optimization is a direct retained-boundary product
kernel:

```text
left retained mode tuples
right retained mode tuples
1D axis factors
-> retained boundary block directly
```

For overlap:

```text
S[a,b] =
    Sx[ix_a, ix_b] * Sy[iy_a, iy_b] * Sz[iz_a, iz_b]
```

For kinetic:

```text
K[a,b] =
    Kx[ix_a, ix_b] * Sy[iy_a, iy_b] * Sz[iz_a, iz_b]
  + Sx[ix_a, ix_b] * Ky[iy_a, iy_b] * Sz[iz_a, iz_b]
  + Sx[ix_a, ix_b] * Sy[iy_a, iy_b] * Kz[iz_a, iz_b]
```

The existing dense raw-source block plus selector path should remain a
small-fixture oracle/reference, not the production scaling shape.

After overlap/kinetic, apply the same direct retained-boundary idea to
centered Gaussian/electron-nuclear factors. Do not start with density-density
or RHF.

## Shell-Support Operator Policy

`pqs_source_shell_projected_one_body_matrix(...)` is useful as an oracle or
comparison helper:

```text
shell-support operator
-> shell projection
-> final operator
```

It should not become the production PQS operator path unless an explicit policy
change promotes it. The active production direction should remain:

```text
source factors / raw source block
-> retained boundary source-mode operator
-> final-basis operator via Lowdin cleanup
```

## Test Discipline

Do not keep growing
`test/nested/cartesian_pair_block_materialization_contract_runtests.jl` as a
development notebook.

Near-term test rule:

```text
Implementation pass:
    one compact module-contract check max

Physics probe:
    tmp/work only unless promoted to an acceptance gate

Acceptance gate:
    one compact endpoint test, and it must replace or shrink older oracle
    pressure
```

New PQS tests should not be added to the large CPBM contract file unless the
edit is a tiny pure CPBM kernel contract. PQS final-basis behavior should move
to either:

```text
a compact FinalBasisRealization module-contract test
or one compact physics/workflow H1 smoke gate
```

If a durable final H1 gate is added, keep only:

```text
source dims/count: 5 x 5 x 5 / 125
retained boundary count: 98
final retained count: 98
final overlap identity
H finite/symmetric
ordinary eigensolve, no generalized overlap solve
H/eigenvalue agreement with shell-support oracle
no _pqs_current_route_safe_term_matrices
no IDA/RHF/driver/export/artifact claim
```

Do not assert every report field, nonclaim flag, timing field, or blocker
vocabulary.

## Validation Policy

The slow nested harness is not a routine baton-loop validation target. In
unattended baton work, poll baton files about once per minute and stop the
unattended loop if the expected response has not appeared after about one hour.

Validation commands may run longer than one minute. If validation progress is
unclear for a long time, prefer a response checkpoint with the last visible
output and validation status over indefinite waiting. A user may set a shorter
temporary limit for a specific validation run; that should be recorded as a
run-specific instruction, not promoted to standing baton policy.

Doer agents must not request UI escalation during unattended baton mode. If a
command needs permission or sandbox escape, write `.agent_handoffs/ATTENTION.md`
with the exact command, reason, and blocker, then stop.

## Near-Term Sequence

### Pass 1: Status And Ownership Corrections

- Fix the PQS framework wording so `CartesianRawProductSources` owns the narrow
  PQS source-mode boundary selector metadata, not general retained-rule policy.
- Add/review `review.031.md` for the incomplete slow validation.
- Record the no-escalation baton rule and the one-hour unattended polling
  stop/checkpoint rule.
- Record that the final H1 seam is oracle-backed and validated, not a
  production-owned full PQS route.

### Pass 2: Final-Basis Realization Module Extraction

- Create `CartesianFinalBasisRealization`.
- Move only stable final-basis realization/operator-transfer functions.
- Keep function names and return shapes mostly unchanged.
- Update existing references.
- Do not redesign result types, add IDA/RHF/density-density, add drivers, or
  grow broad tests.

### Pass 3: Direct Retained-Boundary Overlap/Kinetic

- Add a direct retained-boundary product kernel for PQS overlap and kinetic.
- Compare against existing raw-source block then selector on a tiny synthetic
  fixture.
- Keep the raw-source selector path as oracle/reference.
- Do not touch shell realization or H1 in this pass.

### Pass 4: Direct Retained-Boundary Nuclear

- Build retained by-center nuclear blocks directly from retained mode tuples
  and centered Gaussian source factors.
- Compare against the existing raw-source nuclear block then selector path.
- Preserve the uncharged by-center convention.

### Pass 5: Compact Final H1 Gate And Oracle Shrink

- Promote the final H1 probe only if it replaces or shrinks older oracle
  pressure.
- Shrink or quarantine `_pqs_current_route_safe_term_matrices(...)` callers
  that only preserve old helper vocabulary.

## Non-Goals

Do not do these as part of the near-term cleanup:

```text
create a broad PQS everything-module
introduce one giant route-result struct
rewrite all NamedTuples into structs
add IDA weights and density-density before physical H1 is sensible
add RHF before H1 and self-Coulomb diagnostics are stable
promote shell-support operator projection to production path
expand slow integration tests as routine gates
```
