# Producer-Wide Coulomb Accuracy Policy

Status: approved narrow source/design authority under
`HP-PQS-COULOMB-ACCURACY-FN-01` and
`HP-PQS-COULOMB-ACCURACY-TEST-01`.

This authority adds one expert accuracy choice to the Cartesian/PQS producer
and requires one resolved `CoulombGaussianExpansion` to govern every
Coulomb-expanded part of that Hamiltonian construction. It does not change the
default, expose custom expansion parameters, or authorize a solver or
Cr2-specific workflow.

## Physics Target

Cr and Cr2 consumers need internally consistent high-accuracy Hamiltonians.
The current producer selects the compact expansion independently in
parent/PGDG construction, base unit-nuclear and IDA assembly, residual-GTO
mixed/self and augmented unit-nuclear construction, and residual
matched-width Gaussian (MWG) interaction assembly.

Choosing high accuracy at only one of those points would not define one
Hamiltonian approximation. Parent factor packets, base `V_GG`, augmented
unit-nuclear blocks, and residual-containing MWG blocks must use the same
expansion.

## Expert Input And Presets

The source-backed base and supplemented producer facades may accept:

```julia
coulomb_accuracy = :compact  # default
coulomb_accuracy = :high
```

The option belongs with producer basis/construction inputs. This authority does
not add it to `bin/cartesian_ham_builder.jl` or create a canonical CLI option.
It is route-family-neutral wherever current PQS and White-Lindsey constructions
share the parent/base/supplemented machinery; neither route may re-resolve a
different expansion.

The only accepted values and exact presets are:

| policy | `doacc` | terms | `del` | `s` | `c` | `maxu` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `:compact` | `false` | 45 | 0.6 | 0.5 | 0.03 | 27.0 |
| `:high` | `true` | 135 | 1.0 | 0.16 | 0.01 | 135.0 |

These are names for the existing deterministic
`coulomb_gaussian_expansion(...)` presets. They are not a new fitting
algorithm. The default of the general expansion utility is outside this lane;
the Cartesian/PQS producer default remains `:compact`.

Do not expose `doacc`, `del`, `s`, `c`, `maxu`, coefficient vectors,
exponent vectors, or custom expansion objects as new user inputs.

## One-Expansion Construction Contract

The producer must resolve the policy exactly once, before parent-axis PGDG
construction:

```text
producer input
-> resolve policy
-> one CoulombGaussianExpansion
-> parent/PGDG factors
-> base unit-nuclear and Vee
-> residual-GTO exact Coulomb-expanded blocks
-> residual MWG interaction
-> artifact summary
```

The existing `CoulombGaussianExpansion` is the construction-stage object.
Carry that object through the working-basis construction; do not copy its seven
summary fields into route or stage records. A small summary is appropriate only
at serialization/reporting boundaries.

Required behavior:

1. Resolve the expansion before `cartesian_parent(...)` or equivalent parent
   construction.
2. Build all parent-axis PGDG factor packets with that expansion's exponent
   vector.
3. Use the same carried expansion for base unit-nuclear attraction and base
   IDA electron-electron assembly.
4. Use it for residual-GTO mixed/self Coulomb blocks, augmented unit-nuclear
   blocks, and residual-containing MWG assembly.
5. Require exact term-count and exponent-order parity between the carried
   expansion and every PGDG axis packet consumed by base or augmented
   assembly.
6. Fail on mismatch. Do not silently rebuild compact packets, truncate the high
   expansion, or mix coefficients from one preset with exponents from another.

The current MWG blanket rejection of an explicit expansion may be removed only
because the producer now supplies Hamiltonian-wide expansion authority. MWG
must validate its expansion against the parent PGDG exponent sequence before
constructing residual-containing interaction blocks.

The old private `_cartesian_base_ida_hamiltonian(...)` helper currently
chooses compact accuracy independently. A focused caller scan should decide its
fate:

- if no live caller remains, delete it;
- if it is still live, require an explicit carried expansion argument.

Do not preserve an uncalled helper through an adapter, and do not leave any
private base helper free to select its own expansion.

## Stable Analytic Gaussian Amendment

The first high-accuracy implementation attempt correctly stopped after the
135-term expansion produced raw pair factors near `1e235` and terminal
non-finite values. Follow-up measurement showed that the failure is
catastrophic cancellation in analytic Gaussian formulas, not a fundamental
PGDG carrier or terminal-contraction limit.

This amendment approves algebraically equivalent stable formulas in exactly
three functions:

```text
GaussianAnalyticIntegrals.gaussian_factor
GaussianAnalyticIntegrals.gaussian_pair_factor
CartesianGaussianRawBlocks._factor_axis_integral
```

For `gaussian_factor`, with `alpha_g = 2g` and
`A = alpha_a + alpha_b + alpha_g`, replace
`sum(alpha*c^2) - A*mean^2` by the pairwise weighted-distance identity:

```text
Q = (
      alpha_a*alpha_b*(c_a - c_b)^2
    + alpha_a*alpha_g*(c_a - c_g)^2
    + alpha_b*alpha_g*(c_b - c_g)^2
    ) / A
```

The value remains `sqrt(2pi/A) * exp(-Q/2)`.

For `gaussian_pair_factor`, do not form the determinant as a subtraction
of two `O(g^2)` terms. Use:

```text
D = alpha_a*alpha_b + 2g*(alpha_a + alpha_b)
Q = 2g*alpha_a*alpha_b*(c_a - c_b)^2 / D
value = 2pi/sqrt(D) * exp(-Q/2)
```

For `_factor_axis_integral`, with
`gamma = alpha_l + alpha_r + alpha_f`, compute the constant as:

```text
constant = (
      alpha_l*alpha_r*(c_l - c_r)^2
    + alpha_l*alpha_f*(c_l - c_f)^2
    + alpha_r*alpha_f*(c_r - c_f)^2
    ) / gamma
```

The polynomial moment and prefactor convention remains unchanged.

These are numerical rewrites of the same integrals. Do not clamp negative
intermediate values, take absolute values, reduce the exponent range, or add a
scaled/log PGDG carrier to mask the cancellation. If another analytic formula
outside these three functions fails the high-range oracle, stop and request a
separate amendment rather than sweeping the file.

Because `GaussianAnalyticIntegrals` is shared, existing ordinary/Qiu-White
callers may inherit the stable evaluation. That is algebraic kernel repair, not
authority for route-specific rewiring, cleanup, default changes, or new
ordinary/Qiu-White validation frameworks.

The manager audit using only these stable identities found:

- H/H2 135-term pair-factor maximum scale below `2.6`;
- tightest-term scale about `1.60e-7`;
- finite base IDA and unit-nuclear matrices with symmetry error near roundoff;
- finite supplemented H2 residual-GTO products and MWG interaction;
- bounded supplemented runtime about `36.1 s`.

No scaled carrier, log carrier, terminal IDA redesign, or broad PGDG change is
needed before resuming the approved producer-wide wiring.

## Numerical And Object Boundaries

This lane does not add expansion fields to `CartesianIDAHamiltonian`. That
object remains the finished matrix/electron-count Hamiltonian. Expansion
authority belongs to the producer working object and artifact provenance, not
to numerical matrix ownership.

The source pass may add one compact expansion field to the current base
working-basis construction. It must not add a flat policy/term/parameter field
cloud to route stages or duplicate the expansion summary across base and
augmentation objects.

The high preset changes the numerical Hamiltonian and is expected to increase
term-dependent time and memory. It does not change shellification, terminal
basis dimensions, residual-selection policy, injection policy, EGOI, or
screened-Hartree equations.

## Hamiltonian Artifact Contract

New base and supplemented Hamiltonian artifacts must write one
Hamiltonian-wide summary as one compact group:

```text
coulomb_expansion/policy
coulomb_expansion/doacc
coulomb_expansion/term_count
coulomb_expansion/del
coulomb_expansion/s
coulomb_expansion/c
coulomb_expansion/maxu
```

The summary describes the expansion used for parent/PGDG, base IDA, exact
augmented Coulomb-expanded blocks, and MWG. A supplemented artifact must not
write separate base and augmentation policies.

The ordinary matrix-only `read_cartesian_ida_hamiltonian(...)` contract may
continue to ignore this provenance. New protected-localized artifacts and
protected ladder members/manifests must preserve and expose the summary on
readback because consumers use those artifacts to resume expensive Cr/Cr2
work.

Legacy artifacts without the group remain readable where they were readable
before. Missing provenance is `nothing`/unavailable; it must never be
inferred as `:high`. A new artifact claiming `:compact` or `:high` must
validate the stored term count and parameters against the named deterministic
preset.

Protected-localized matrices remain in their existing native order, and this
summary does not change protected-localized interaction semantics, sector maps,
or ladder transfer rules.

## Atomic Reference Packet Exception

An atomic HF reference packet is not one producer Hamiltonian. It records
several separately evaluated reference objects, so its expansion provenance is
role-qualified:

- pure-GTO packet RHF uses `:high` (`doacc = true`, 135 terms);
- density-fit/self-energy evaluation currently uses `:compact`;
- the fitted-potential broad-tail scaffold currently uses `:compact`.

The packet writer/readback must record those roles explicitly. It may reuse a
small common expansion-summary serializer, but it must not label the whole
packet with one Hamiltonian-wide policy. Fitted density and potential terms
remain evaluation devices, not protected orbitals or producer Coulomb inputs.

Compact packet-local self-energy and potential fitting remain approved because
the measured Cr screened scalar-constant error is about `0.0402 mHa`. This is
a recorded evaluated approximation, not permission to mix compact and high
expansions inside a produced Hamiltonian.

## Approved Source Surface

Only these files are approved:

```text
src/cartesian_base_hamiltonian.jl
src/pqs_source_box_route_driver_helpers.jl
src/pqs_source_box_low_order_materialization.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/cartesian_residual_gaussians/mwg_interaction.jl
src/cartesian_ida_hamiltonian.jl
src/cartesian_protected_ladder_bundle.jl
src/cartesian_reference_density/atomic_hf_reference_packets.jl
src/GaussianAnalyticIntegrals.jl
src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl
```

`src/cartesian_ida_hamiltonian.jl` is approved only for compact summary
serialization/readback shared by current artifact owners. It must not make
`CartesianIDAHamiltonian` choose or carry a construction expansion.

No new source file, struct, public export, driver input, or general reporting
framework is approved.

## Validation Authority

`HP-PQS-COULOMB-ACCURACY-TEST-01` approves:

- `git diff --check`;
- package load;
- omitted policy versus explicit `:compact` matrix equality for a bounded
  base construction;
- a bounded `:high` base construction with 135-term artifact provenance and
  exact parent/PGDG exponent parity;
- a bounded White-Lindsey base smoke confirming the same policy is carried
  through the shared parent/base machinery;
- a bounded supplemented residual-GTO/MWG construction proving the same
  expansion reaches base, exact augmented, and MWG work;
- protected-localized member and protected ladder manifest write/readback of
  the summary;
- atomic packet roundtrip preserving separate RHF, density/self-energy, and
  potential-tail expansion provenance;
- stage timing and expansion-dependent allocation reporting for at least one
  bounded compact/high comparison;
- terminal due-diligence inspection for every endpoint-style base or
  supplemented probe used to interpret the result.

Existing tests may be updated only where they already own these contracts:

```text
test/driver_public/cartesian_base_hamiltonian_runtests.jl
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
test/nested/cartesian_atomic_hf_reference_packet_runtests.jl
test/core/runtests.jl
```

The core test addition is limited to a small focused analytic-kernel test. It
must compare the stable formulas with a BigFloat oracle across compact and high
exponent ranges, preserve moderate-exponent values within roundoff, and cover
finite nonnegative s-type factors plus translated centers that trigger the old
cancellation.

The high supplemented and protected-ladder checks may remain ignored bounded
probes if adding them to a committed endpoint would materially increase normal
test runtime. No new committed test file or Cr/Cr2 endpoint assertion is
approved.

## Stop Rules

Stop without a source commit and report the exact boundary if:

- parent/PGDG construction cannot consume the carried expansion without a new
  route-stage object or broad carrier redesign;
- base or augmented assembly needs independently generated expansion data;
- MWG cannot prove exponent parity with the parent PGDG packets;
- protected/ladder artifacts cannot distinguish known preset provenance from
  missing legacy provenance without a broader format redesign;
- packet provenance would be mistaken for Hamiltonian-wide policy;
- the source pass needs files outside the approved list.

## Explicit Exclusions

This authority does not approve:

- changing the producer default to `:high`;
- custom expansion parameters or coefficient/exponent inputs;
- canonical CLI changes;
- ordinary Qiu-White, legacy, or experimental path cleanup;
- scaled/log PGDG carriers, new stage objects, or terminal contraction changes;
- shellification, terminal realization, retained selection, mapping, residual
  selection, injection, EGOI, or screened-Hartree formula changes;
- solver/HF/MP2-NO workflow;
- a Cr/Cr2-specific producer branch or committed endpoint assertion;
- changing protected-localized `Vee` or ladder transfer semantics;
- treating atomic packet fitted Gaussians as basis functions.

Target source growth is at most about 300 added lines across the approved
files, offset where practical by deleting independent compact selectors and
the uncalled private base helper. If the change requires a new policy framework
or materially exceeds that size, stop for design review.

The stable-formula amendment itself should remain below about 60 added source
lines and must not introduce a new carrier, cache, status object, or module.
