# Producer-Wide Coulomb Accuracy Policy

Status: compact/high producer policy implemented and maintenance-only under
`HP-PQS-COULOMB-ACCURACY-FN-01` and
`HP-PQS-COULOMB-ACCURACY-TEST-01`.

Pass 616 places standard60 implementation, fingerprint/provenance additions,
canonical-driver exposure, and their new tests on hold pending assessment of
the Gaussian-kernel determinant cancellation reported September 8. The exact
K60 parameters and fingerprint remain retained design, not execution authority.
This restriction overrides pending-implementation language and budgets below;
those portions are deferred specifications, not permission to resume. Kernel
repair also requires separate review and authority. No preset retuning, default,
source, test, artifact-format, or release change is granted by this packet.

One resolved `CoulombGaussianExpansion` must still govern every Coulomb-expanded
part of an existing producer construction. Compact/high behavior is unchanged.

## Physics Target

Cr and Cr2 consumers need internally consistent high-accuracy Hamiltonians.
The implemented compact/high producer now resolves and carries one expansion
through parent/PGDG construction, base unit-nuclear and IDA assembly,
residual-GTO mixed/self and augmented unit-nuclear construction, and residual
matched-width Gaussian (MWG) interaction assembly. The fixed standard tier is
the remaining producer extension, now on hold.

Choosing high accuracy at only one of those points would not define one
Hamiltonian approximation. Parent factor packets, base `V_GG`, augmented
unit-nuclear blocks, and residual-containing MWG blocks must use the same
expansion.

## Expert Input And Presets

The source-backed base and supplemented producer facades currently accept:

```julia
coulomb_accuracy = :compact  # default
coulomb_accuracy = :high
```

The retained, deferred design names the following option, which committed
source does not accept and this authority no longer permits implementing:

```julia
coulomb_accuracy = :standard
```

The option belongs with producer basis/construction inputs. The canonical
human-facing driver does not yet expose it; a separately reauthorized amendment may
add the same policy name and default, but the driver never owns a second policy
or expansion resolver. The option is route-family-neutral
wherever current PQS and White-Lindsey constructions share the
parent/base/supplemented machinery; neither route may re-resolve a different
expansion.

The existing presets and retained standard60 design are:

| policy | `doacc` | terms | `del` | `s` | `c` | `maxu` | role |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `:compact` | `false` | 45 | 0.6 | 0.5 | 0.03 | 27.0 | cheapest legacy approximation |
| `:standard` | `false` | 60 | 1.0 | 0.34257593251905827 | 0.042605721927199074 | 60.0 | atomic/local Gaussian-integral accuracy-cost option; implementation held |
| `:high` | `true` | 135 | 1.0 | 0.16 | 0.01 | 135.0 | higher-accuracy finite reference approximation |

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
identity. The deferred design allows one private fixed-preset constructor in
the existing Coulomb expansion owner; it must not change compact/high bit
patterns or expose a new public custom-expansion interface.

`doacc` is a legacy compatibility field, not preset identity: it is `false`
for both `:compact` and `:standard`. Policy, exact parameters, term count, and
the coefficient/exponent fingerprint define a preset. These are fixed
quadratures, not runtime fits. The default of the general expansion utility is
outside this lane; the Cartesian/PQS producer default remains `:compact`.

Standard60 is an atomic/local Gaussian-integral accuracy-cost option, not a
monotone intermediate tier for pointwise, diffuse, or long-range accuracy.
Molecular extent, diffuseness, and PQS use require validation before broader
recommendation. The controlled Cr2 screened off/on comparison remains `:high`;
this neither certifies high135 at arbitrary scales nor permits changing that
comparison. Defaults remain owner-specific: the general expansion utility uses
high, while the Cartesian producer and pair-matrix utility use compact.

### September 8 Evidence Reconciliation

Reviewed at `1c6a4a67da78f1963b0f0ccc6a99b21324b43439` using
`tmp/reviews/coulomb-review-2026-09-08/REPORT.md` (SHA-256
`8132de61ad6a08bf6b250a7cea3fb7d09fd6d99b445539b54f98bd426b48eda9`).
Recovered TSV rows, precision logs, and current source were independently
inspected; no RHF or numerical suite was rerun for this policy reconciliation.

The recovered RHF control is **Cr6+, 18 electrons, nine doubly occupied
orbitals, and 33 Cartesian s/p orbitals from cc-pV5Z**, not neutral Cr or Cr2.
Both attraction and repulsion used the same expansion. High135 gave
`-1033.8673200814285 Ha`; standard60 gave `-1033.8673200811857 Ha`, a signed
standard-minus-high difference of `+2.4283508537e-10 Ha`. Both converged with
electron trace 18. Standard started from converged high density, so the SCF
timings are not a fair cost comparison. The recorded PySCF value was
`-1033.8673200862515 Ha`; neither finite expansion is an exact reference.
Compact45's approximately `5.75 mHa` difference is a historical Pass 364
claim: its original numerical row was not recovered. Do not label it freshly
verified or substitute an unrelated extreme-primitive calculation.

Separate nuclear/ERI controls, independent of the single RHF energy, found
maximum standard-minus-high absolute differences (Ha):

| Bounded control | Orbitals | Nuclear | ERI |
| --- | ---: | ---: | ---: |
| Cr full s/p | 33 | 7.28e-9 | 3.07e-10 |
| Cr selected s/p/d | 10 | 7.28e-9 | 1.89e-10 |
| H selected 2s+3p | 5 | 1.16e-10 | 1.33e-10 |
| Be selected 2s+3p | 5 | 3.17e-10 | 4.48e-11 |
| Ne selected 2s+3p | 5 | 8.55e-10 | 2.46e-10 |

These are selected atomic controls against high135, not full-basis surveys or
analytic all-integral certification. Historical paired integral construction
cost was 14.628 s for K60 versus 21.810 s for high135, not a current producer
performance promise. Offline tuning used a weighted log-grid objective, not a
uniform-error fit; the exact retained candidate was not a rejected nonlinear fit.

Fresh independent s-Gaussian controls use normalized densities with exponents
`p,q`, separation `R`, and `beta=p*q/(p+q)`. Their exact interaction is
`erf(sqrt(beta)*R)/R`, with limit `2sqrt(beta/pi)` at zero; the expansion is
`sum(c*(beta/(beta+zeta))^(3/2)*exp(-beta*zeta/(beta+zeta)*R^2))`.
Orbital-product overlaps multiply both expressions; nuclear attraction uses
`beta=p` and one overlap factor. These formulas bypass the live determinant
kernel and distinguish quadrature error from implementation error.

For same-center orbital exponents `a,b=2a`, standard60 signed pair errors at
`a=0.01,1,100` were `+2.45e-13,+8.44e-12,-6.59e-13 Ha`: excellent tested
overlapping-density accuracy, supported by 256-bit accumulation. For `a=b=1`
at `R=30,50,100` bohr, errors were `-9.57e-8,+2.49e-7,-1.40e-4 Ha`.
At `a=1e-4,b=2e-4,R=0`, the independent finite-expansion error was
`-5.13e-4 Ha` (about 3.94%); broad densities sample the deficient tail even
at one center. Sampled pointwise relative maxima on `[0.01,50]` were
`5.229e-7,1.272e-5,1.285e-13` for compact/standard/high respectively, not
continuous interval bounds. No uniform integrated guarantee follows.
High135 is also finite: at the artificial normalized `a=1.3554e8,b=2a`
pair its relative error was `6.23e-7`. This omits contraction weights and is
not the recovered Cr6+ energy. It rules out an arbitrary-scale exact-oracle claim.

### Kernel Assessment Hold

At `a=1e-4,b=2e-4`, the actual standard/high pair builder rejects its quadratic
form, although the independent integral is well-defined. In
`src/foundation/GaussianAnalyticIntegrals.jl`,
`centered_polynomial_gaussian_pair_factor_integral` computes
`(beta_left+zeta)*(beta_right+zeta)-zeta^2`; large terms can cancel although
the exact determinant is positive. This is distinct from the previously
repaired `gaussian_pair_factor`, and is not evidence for retuning K60.
Independent inspection also qualifies the report's rounding-level agreement:
selected `(11|22)` entries agree, but `integrals.tsv` records large complete
pair-matrix discrepancies at R=100/300 for standard/high. Their cause is not
established here; do not certify full matrices from the selected-entry result.
Assess numerical range and these discrepancies separately before reauthorizing
standard60. No kernel fix, clamping, new regression, or implementation is
granted here. Existing compact/high maintenance excludes this new repair.

Do not expose `doacc`, `del`, `s`, `c`, `maxu`, coefficient vectors,
exponent vectors, or custom expansion objects as new user inputs.

## Canonical Driver Exposure

The canonical driver does not currently expose this input. The following
deferred design requires fresh authority after the kernel assessment:

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

The old private `_cartesian_base_ida_hamiltonian(...)` helper was deleted after
its caller scan proved empty. Do not restore it through an adapter or leave any
new private base helper free to select its own expansion.

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

The existing maintenance paths are listed below. Pending additions described
in this document are held by Pass 616; this list does not reopen them.

```text
src/cartesian/cartesian_base_hamiltonian.jl
src/cartesian/pqs_source_box_route_driver_helpers.jl
src/cartesian/pqs_source_box_low_order_materialization.jl
src/cartesian/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/cartesian/cartesian_residual_gaussians/mwg_interaction.jl
src/cartesian/cartesian_ida_hamiltonian.jl
src/cartesian/cartesian_protected_ladder_bundle.jl
src/cartesian/cartesian_reference_density/atomic_hf_reference_packets.jl
src/ordinary/ordinary_coulomb.jl
src/foundation/GaussianAnalyticIntegrals.jl
src/cartesian/cartesian_gaussian_raw_blocks/nuclear_blocks.jl
bin/cartesian_ham_builder.jl
```

`src/cartesian/cartesian_ida_hamiltonian.jl` is approved only for compact summary
serialization/readback shared by current artifact owners. It must not make
`CartesianIDAHamiltonian` choose or carry a construction expansion.

No new source file, struct, public export, driver input other than the exact
`coulomb_accuracy` symbol above, or general reporting framework is approved.

## Validation Authority

`HP-PQS-COULOMB-ACCURACY-TEST-01` maintains existing compact/high coverage only.
Its completed/maintenance state does not certify the unimplemented standard60
or driver extension. The following original acceptance checklist is retained
for future review, not as an execution grant for new tests or smokes:

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
