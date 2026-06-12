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

The durable final H1 gate is now the complete one-center core plus surrounding
shell route, not the earlier boundary-shell-only mechanical fixture:

```text
current_box: 1:7 x 1:7 x 1:7
inner direct core: 2:6 x 2:6 x 2:6
raw source dims: 5 x 5 x 5
core support count: 125
surrounding shell support count: 218
shell final retained count: 98
complete final retained count: 223
final overlap identity
H finite/symmetric
ordinary eigensolve, no generalized overlap solve
H1 ~= -0.48047934800387226
same-geometry fixed-block oracle H1 ~= -0.48047920531279725
electron-nuclear factor source: pgdg_intermediate.gaussian_factor_terms
no _pqs_current_route_safe_term_matrices
no IDA/RHF/driver/export/artifact claim
```

The old `98`-function boundary-shell-only H1 path should be treated as
mechanical/source-box coverage or an oracle aid, not as a physical one-center
acceptance basis. Do not assert every report field, nonclaim flag, timing
field, or blocker vocabulary.

Subsequent complete core/shell PQS He RHF probes validated the final/pre-final
density convention and route scaling, but they are not yet accepted physical
convergence fixtures. The q ladder used one direct core plus exactly one
surrounding shell while keeping the mapped box radius near 8 bohr:

```text
q=5   final dim 223    RHF -2.7213372828531668
q=7   final dim 561    RHF -2.810068050134403
q=9   final dim 1115   RHF -2.8499091618019303
q=11  final dim 1933   RHF -2.8559475204289022
```

Those values are strong route/scaling evidence, especially for the localized
pre-final positive-weight density gauge consumed by final-basis orbitals through
the combined Lowdin cleanup. They do not prove physical convergence because
increasing q alone does not review the number of shell layers,
mapping/distortion parameters, central spacing, or box-size policy together.
Do not promote q=9 or q=11 as a permanent He RHF gate until the fixture is
reviewed. The q=11 point remains the stronger exploratory/reference point; q=9
is the cheaper candidate if a compact gate is later approved.

The next PQS physical probes should be WL-aligned fixture probes rather than
more q-only scaling points. In WL terms, the fixture must choose the parent box
radius, central spacing `d`, mapping/distortion `s`, and shell depth together.
The current complete core/shell PQS helper still uses a projected shell source
with a one-cell raw boundary around the inner box. A side-13, `d = 0.1`,
`s = 1.0`, three-surrounding-shell final-basis smoke is therefore blocked on a
multi-layer PQS shell/source producer before H1 or RHF should be interpreted.
A follow-up seam audit showed that this does not require immediate
generalization of the final-basis helper. Three legal one-cell shell descriptors
for `(3:11)^3/(4:10)^3`, `(2:12)^3/(3:11)^3`, and `(1:13)^3/(2:12)^3` have
disjoint supports, cover the side-13 parent together with the `(4:10)^3` core,
and can be collapsed into one block-diagonal shell sector before calling
`pqs_complete_core_shell_final_basis`. The smallest next implementation target
is therefore a route-owned multi-layer PQS shell source plan that builds and
combines repeated one-cell shell descriptors; old fixed-block matrices remain
oracle/reference material only.

That seam is now implemented as `pqs_multilayer_shell_source_plan(...)`. The
first side-13 final-basis/H1 smoke, using `AsinhMapping(c = 0.1, s = 1.0,
tail_spacing = 10.0)`, core `(4:10)^3`, and outer box `(1:13)^3`, materialized
three shell layers, 1,854 shell support rows, 1,206 shell-retained columns, and
a 1,549-dimensional final basis. The final overlap identity error was about
`5.51e-13`; the Z = 1 H1 energy was about `-0.494223730383033`, and the Z = 2
H1 energy was about `-1.975561823201342`. This is a route smoke/probe only, not
an acceptance gate or RHF fixture. The mapping relation among `Z`, core spacing
`d`, distortion `s`, radius, and shell depth remains provisional for PQS.
The matching pre-RHF density diagnostic materialized the pre-final density
interaction with positive weights, a finite symmetric `(1549, 1549)` pair
matrix, and no signed-final-weight, raw-no-division, or fixed-block pair
authority. The Z = 2 H1 energy was about `-1.975561823201342`, and the
self-Coulomb diagnostic was about `1.216926438886032` versus the hydrogenic
`5Z/8 = 1.25` and the WL side-13 value near `1.215829476773570`. This makes a
side-13 PQS RHF probe reasonable as the next probe, but still not a gate.
The follow-up side-13 PQS RHF probe converged in 8 iterations with one-electron
energy about `-3.846945861201986`, electron-electron energy about
`1.009690214812516`, and total energy about `-2.837255646389471`. The result is
about `+0.024424349222768` Hartree above the He HF reference and about
`-0.000757646688557` Hartree relative to the WL side-13 RHF probe. The final
and pre-final density traces were both 1 to roundoff, giving two electrons
under the restricted closed-shell convention. This remains a non-acceptance
reference probe and should be reviewed as a future compact gate candidate only
after the PQS fixture rule is settled.

## Validation Policy

The slow nested harness is not a routine baton-loop validation target. In
unattended baton work, poll baton files about once per minute and stop the
unattended loop if the expected response has not appeared after about one hour.

Validation commands may run longer than one minute. If validation progress is
unclear for a long time, prefer a response checkpoint with the last visible
output and validation status over indefinite waiting.

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
