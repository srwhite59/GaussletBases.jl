# Producer-Wide Coulomb Accuracy Policy

Status: compact/high producer policy implemented; the fixed standard tier and
narrow canonical-driver exposure are approved for implementation under
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

The source-backed base and supplemented producer facades accept:

```julia
coulomb_accuracy = :compact  # default
coulomb_accuracy = :standard
coulomb_accuracy = :high
```

The option belongs with producer basis/construction inputs. The canonical
human-facing driver may expose the same policy name and default; it does not
own a second policy or expansion resolver. The option is route-family-neutral
wherever current PQS and White-Lindsey constructions share the
parent/base/supplemented machinery; neither route may re-resolve a different
expansion.

The only accepted values and exact presets are:

| policy | `doacc` | terms | `del` | `s` | `c` | `maxu` | role |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `:compact` | `false` | 45 | 0.6 | 0.5 | 0.03 | 27.0 | cheapest legacy approximation |
| `:standard` | `false` | 60 | 1.0 | 0.34257593251905827 | 0.042605721927199074 | 60.0 | normal accuracy/cost balance |
| `:high` | `true` | 135 | 1.0 | 0.16 | 0.01 | 135.0 | reference-grade production |

`:compact` and `:high` retain their existing deterministic generator paths.
`:standard` is the fixed analytic K60 quadrature with the existing implicit
midpoint phase `theta = 0.5`. For `k = 0:59`, its canonical coefficient and
exponent construction is:

```text
u_k = s * (k + 0.5)
x_k = c * sinh(u_k)
exponent_k = x_k^2
coefficient_k = (2s/sqrt(pi)) * sqrt(x_k^2 + c^2)
```

The coefficient vector followed by the exponent vector, encoded as canonical
little-endian Float64 bytes, has SHA-256 fingerprint:

```text
2de3ec44fc3d6b11ea26b7551e6b5ddef8bb2de1898fe0702d65f91cbf6c0f3a
```

The fixed operation order above is part of the deterministic preset. Calling
the generic keyword-override utility with algebraically equivalent parameters
can differ in final Float64 bits and therefore does not establish preset
identity. The implementation may add one private fixed-preset constructor in
the existing Coulomb expansion owner; it must not change compact/high bit
patterns or expose a new public custom-expansion interface.

`doacc` is a legacy compatibility field, not preset identity: it is `false`
for both `:compact` and `:standard`. Policy, exact parameters, term count, and
the coefficient/exponent fingerprint define a preset. These are fixed
quadratures, not runtime fits. The default of the general expansion utility is
outside this lane; the Cartesian/PQS producer default remains `:compact`.

The accepted bounded evidence is:

- compact45 misses the Cr closed-shell s/p-only RHF control by about
  `5.75 mHa`;
- standard60 differs from high135 by about `2.43e-10 Ha` for that control;
- H, Be, Ne, and Cr s/p/d integral controls agree closely with high135.

Standard60 does not strictly dominate compact45 at extremely long range. It is
the recommended opt-in for serious atomic and ordinary molecular work, not a
claim of uniform pointwise superiority. The controlled Cr2 screened off/on
calculation remains `:high`; changing it to `:standard` would confound that
comparison. A producer-default change requires separate authority after
bounded molecular/PQS validation.

Do not expose `doacc`, `del`, `s`, `c`, `maxu`, coefficient vectors,
exponent vectors, or custom expansion objects as new user inputs.

## Canonical Driver Exposure

`bin/cartesian_ham_builder.jl` may add exactly one public expert input:

```julia
coulomb_accuracy = :compact  # :compact, :standard, or :high
```

Required driver behavior:

1. Keep the visible editable default `:compact` near the other basis inputs.
2. Include `:coulomb_accuracy` in the existing trusted input-file and
   `key=value` override allowlist.
3. Normalize with `Symbol(...)` and reject values outside
   `(:compact, :standard, :high)` before construction.
4. Add the normalized symbol to the existing `common_basis` `NamedTuple`.
5. Print the resolved policy in the existing basis contract summary.

The driver must pass the symbol through the public basis contract and let the
existing producer resolve the one `CoulombGaussianExpansion`. It must not call
`coulomb_gaussian_expansion(...)`, inspect coefficients/exponents, add a second
default, or create a new parser/configuration object. Omitted driver input must
remain exactly equivalent to explicit `coulomb_accuracy = :compact`.

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
coulomb_expansion/fingerprint
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
before. Existing compact/high summaries without `fingerprint` remain readable
as legacy provenance, but the missing fingerprint is unavailable and must not
be invented. No legacy `:standard` summary exists: a summary claiming
`:standard` without the exact fingerprint is invalid. Every new summary for
any named preset writes and validates the fingerprint together with policy,
term count, and parameters. Missing provenance must never be inferred as
`:standard` or `:high`.

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
src/ordinary_coulomb.jl
src/GaussianAnalyticIntegrals.jl
src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl
bin/cartesian_ham_builder.jl
```

`src/cartesian_ida_hamiltonian.jl` is approved only for compact summary
serialization/readback shared by current artifact owners. It must not make
`CartesianIDAHamiltonian` choose or carry a construction expansion.

No new source file, struct, public export, driver input other than the exact
`coulomb_accuracy` symbol above, or general reporting framework is approved.

## Validation Authority

`HP-PQS-COULOMB-ACCURACY-TEST-01` approves:

- `git diff --check`;
- package load;
- omitted policy versus explicit `:compact` matrix equality for a bounded
  base construction;
- a bounded `:standard` base construction with 60-term artifact provenance,
  exact K60 fingerprint, and parent/PGDG exponent parity;
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
  bounded compact/standard/high comparison;
- terminal due-diligence inspection for every endpoint-style base or
  supplemented probe used to interpret the result.

Existing tests may be updated only where they already own these contracts:

```text
test/driver_public/cartesian_base_hamiltonian_runtests.jl
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
test/nested/cartesian_atomic_hf_reference_packet_runtests.jl
test/core/runtests.jl
test/docs/cartesian_ham_builder_policy_runtests.jl
```

The existing docs policy test may assert only the canonical-driver contract:
visible compact default, public-input allowlisting, value validation,
`common_basis` forwarding, contract printing, and absence of private expansion
resolution. Do not create a new committed driver test or input fixture.

Source validation must also run bounded one-center canonical-driver smokes for:

- omitted policy and explicit `:compact`, with exact matrix/artifact parity;
- explicit `:standard`, with successful input acceptance, finite/symmetric
  output, and exact `:standard`/60-term/fingerprint provenance;
- explicit `:high`, with successful input acceptance, finite/symmetric output,
  and `coulomb_expansion/` provenance reporting `:high` and `135` terms.

These may use ignored temporary input/output paths. Inspect terminal due
diligence, but do not add a high-accuracy endpoint or energy baseline.

The core test addition is limited to a small focused analytic-kernel test. It
must compare the stable formulas with a BigFloat oracle across compact,
standard, and high exponent ranges, preserve moderate-exponent values within
roundoff, reproduce the exact standard K60 coefficient/exponent fingerprint,
and cover finite nonnegative s-type factors plus translated centers that
trigger the old cancellation.

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
- changing the producer default to `:standard`;
- custom expansion parameters or coefficient/exponent inputs;
- canonical driver inputs or CLI behavior beyond the single policy symbol
  approved above;
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

The driver-only portion should target at most 25 added `bin`/test lines. If
that portion needs a parser abstraction, new configuration carrier, or another
committed test file, stop and report the missing seam. The fixed K60 resolver
and fingerprint provenance belong to the producer/artifact implementation
portion approved above, not to the driver.
