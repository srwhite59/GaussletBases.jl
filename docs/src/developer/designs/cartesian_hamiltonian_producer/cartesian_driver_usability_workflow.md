# Cartesian Driver Usability Workflow

Status: implemented canonical contract for the human-facing Cartesian
Hamiltonian driver and its non-exported staged producer calls, plus one
implemented private PQS/WL paper-validation driver with an implemented H2+
endpoint, implemented neutral-H2 supplemented one-body preflight, and
implemented five-row fixed-state density-density measurement. The current
PQS/WL paper rows are provisional because they predate the approved matched
shared-shell aspect correction; replay under that correction precedes the
external same-density oracle.
Registry entries own
the lifecycle and source permissions for the corresponding IDs. Driver
validation IDs without committed fixtures are completed evidence or explicitly
named probe records, not continuing test authority.

## Boundary

The canonical driver is:

```text
bin/cartesian_ham_builder.jl
```

It is a trusted local scientific script, not a public parser or a second
library API. It owns editable defaults, trusted input loading, command-line
overrides, construction of visible `system`, `basis`, and `supplement`
records, terminal due-diligence presentation, coarse stage timing, artifact
production, and optional readback.

The exported library interface remains the base-only facade:

```julia
cartesian_base_hamiltonian(system; basis, hamfile = nothing)
```

That facade accepts exact public records and returns
`CartesianIDAHamiltonian{Float64}` directly. The driver instead calls the
non-exported staged functions in `src/cartesian_base_hamiltonian.jl` so its
physics-level timings remain visible and so it can compose the separately
implemented supplemented workflow. Driver variable names are not facade
keywords, and driver defaults are not hidden producer defaults.

The canonical driver does not currently accept or forward
`coulomb_accuracy`. Omission therefore selects the producer's `:compact`
default. The facade's implemented `:high` opt-in and the unimplemented
`:standard` tier do not imply canonical-driver support.

## Invocation And Inputs

```text
julia --project=. bin/cartesian_ham_builder.jl [input.jl] [key=value ...]
```

If the first argument has no `=`, it is included as trusted Julia code. It may
assign the public driver variables or return a `NamedTuple` or `AbstractDict`
of replacements. Later `key=value` arguments are trusted Julia expressions
and take precedence. Unknown driver keys throw `ArgumentError`.

The implemented inputs and checked-in defaults are:

| Input | Default | Meaning |
| --- | --- | --- |
| `Natom` | `2` | `1` for an origin-centered atom; `2` for a homonuclear z-axis diatomic |
| `atom` | `"H"` | per-center label only |
| `Z` | `1.0` | explicit per-center nuclear charge |
| `R` | `4.0` | full diatomic bond length; unused for an atom |
| `nup`, `ndn` | `1`, `1` | explicit spin-sector electron counts |
| `ns` | `5` | public source/cube/nesting size |
| `nesting` | `:pqs` | `:pqs` or `:wl` |
| `source_span` | `:ordinary` | `:ordinary` or PQS-only `:mapped_comx` |
| `core_spacing` | `0.3` | near-nucleus physical spacing |
| `s_factor` | `1.0` | finite positive expert mapping-strength factor |
| `gausslet_family` | `:G10` | producer `parent_axis_family`; only `:G10` is implemented |
| `padding` | `10.0` | atom radius or padding beyond diatomic nuclei |
| `basisname` | `nothing` | `nothing` selects base; a label selects supplementation |
| `lmax` | `1` | supplement angular cutoff |
| `uncontracted` | `false` | supplement contraction choice |
| `supplement_width_max` | `nothing` | optional positive maximum supplement width |
| `basisfile` | `nothing` | optional trusted local BasisSets path |
| `hamfile` | `"cartesian_ida_hamiltonian.jld2"` | nonempty artifact filename |
| `check_file` | `true` | read the written artifact and compare its dimension |
| `print_contract` | `true` | print the public contract and terminal due diligence |
| `print_timing` | `true` | print coarse stage timings |
| `expected_dimension` | `nothing` | optional exact final-dimension check |

The driver has no public `q`, `d`, `reference_spacing`, `tail_spacing`, route
control, raw-provider control, or Coulomb-policy input. The producer derives
route-local `q`, rejects unsupported combinations, and supplies its own
`reference_spacing = 1.0`, `tail_spacing = 10.0`, and compact Coulomb defaults.

## Contract Construction

For `Natom = 1`, the driver constructs one center at the origin. For
`Natom = 2`, it constructs equal labels and charges at
`(0, 0, -R/2)` and `(0, 0, R/2)`. In both cases the public `system` has exactly
`atom_symbols`, `nuclear_charges`, `atom_locations`, `nup`, and `ndn`.
Symbols are labels; explicit charges and electron counts are authority.

The common visible basis fields are:

```text
ns, core_spacing, s_factor, parent_axis_family, nesting, source_span
```

An atom adds `radius = padding`. A diatomic adds:

```text
xmax_parallel = R/2 + padding
xmax_transverse = padding
```

`nesting = :pqs` derives `q = ns`; `nesting = :wl` derives `q = ns - 2`.
Mapped-COMX source spans are PQS-only, and WL diatomics require `ns >= 4`.
All detailed geometry, electron-sector, spacing, and policy checks remain
producer owned. The driver supplies explicit `nup` and `ndn`; it does not
infer a neutral sector or accept a separate `net_charge` authority.

`basisname === nothing` selects the base workflow. Otherwise the driver builds:

```text
basis_by_center = fill(String(basisname), Natom)
lmax
uncontracted
width_filtering = nothing | (; max_width = supplement_width_max)
basisfile
```

and selects the shared residual-GTO/MWG supplemented workflow. Supported
origin-centered atoms and homonuclear z-axis diatomics compose with both PQS
and White-Lindsey nesting through the same terminal-basis boundary.

## Staged Workflow

The implemented top-level sequence is:

```text
construct and optionally print public contract
build base working basis and terminal due diligence
build base product/moment operators
build base unit-nuclear operators
build base localized-IDA interaction
assemble base Hamiltonian
[load supplement and select owner-local residual Gaussians]
[build augmented products, unit-nuclear operators, and residual MWG/IDA]
[assemble supplemented Hamiltonian]
write artifact
check expected dimension and optional readback
print dimension and coarse timings
```

Bracketed stages run only when `basisname !== nothing`. For a base run, base
assembly writes `hamfile`. For a supplemented run, the intermediate base
Hamiltonian stays in memory and final supplemented assembly writes `hamfile`.
The result is always the existing `CartesianIDAHamiltonian{Float64}` artifact,
with sidecars owned by the base or supplemented producer and the artifact
manifest contracts.

When `print_contract = true`, the driver also prints the terminal
due-diligence report: normalized geometry, parent bounds and axes, weight
statistics, dimensions, terminal rows, retained/source counts, slab topology,
and warning flags. This is a bounded human review surface, not a route report.
Consumers must inspect it before interpreting endpoint, residual, injection,
screened-Hartree, EGOI, Be2, or Cr2 results.

`expected_dimension` and `check_file` run after construction and artifact
writing. A dimension mismatch throws `ArgumentError`; `check_file` uses
`read_cartesian_ida_hamiltonian` and checks the readback dimension. It does not
write a separate check record or validate every provenance sidecar.

## Source Ownership

- `bin/cartesian_ham_builder.jl` owns the canonical script, inputs,
  presentation, staged calls, and artifact/readback workflow.
- `src/cartesian_base_hamiltonian.jl` owns the non-exported staged producer
  composition used by the script.
- `src/pqs_source_box_low_order_materialization.jl` and
  `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl` own the
  existing behavior-preserving base operator-class factoring.
- `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` owns the
  existing supplemented operator-class factoring and artifact write.

Artifact format and reader behavior remain with their canonical artifact
owners. This workflow owns no other `bin`, `src`, `tools`, committed test, or
committed input-fixture file.

## Failure Behavior

The driver rejects unknown keys, `Natom` outside `1:2`, unsupported `nesting`
or `source_span`, and an empty `hamfile` before construction where practical.
The producer then rejects malformed records, invalid electron sectors,
unsupported geometry, invalid sizes/spacings/mapping factors, incompatible
source spans, and invalid supplement inputs. Input evaluation, numerical
construction, artifact writing, and readback errors propagate. Failures are
not converted to status payloads or partial-success objects.

## Non-Goals

This contract does not add route-stage controls, stop points, raw-provider
switches, diagnostic dumps, allocation probes, custom Coulomb parameters,
solver/RHF/ECP/EGOI behavior, Cr2-specific branches, public exports, result
wrappers, artifact schemas, committed input fixtures, or a safe untrusted-code
parser. Route and ladder diagnostics remain outside the canonical driver.

## Private PQS/WL Paper H2 Driver

`HP-PQS-PAPER-H2-DRV-FN-01` owns one tracked, non-exported scientific
validation script:

```text
bin/pqs_paper_h2_driver.jl
```

This script is private to the matched full-parent/White-Lindsey/PQS H2+ paper
campaign and the one frozen neutral-H2 supplemented one-body preflight below.
All other H2 endpoints remain outside the current implementation authority.
It is not a public API, a general paper-driver framework, a replacement for
`cartesian_ham_builder.jl`, or a source of producer defaults. It may be removed
after the campaign if no durable use remains.

### Required Construction

The private driver continues to use the staged route because its full-parent
reference and terminal diagnostics are not returned by the public base facade;
charged-sector support does not replace this campaign-specific construction.
For each bare H2+ method, the script must carry the actual two-center system,
with `nup=1`, `ndn=0`, and nuclei at
`(0,0,-R/2)` and `(0,0,R/2)`, through the existing staged multicenter route:

```text
cartesian_system
-> cartesian_recipe
-> cartesian_parent
-> cartesian_shells
-> cartesian_units
-> cartesian_transforms
```

This route must use the existing positive combined inverse-sqrt parent mapping.
The PQS case uses the current source-box-first terminal realization and the
White-Lindsey case uses the corresponding current terminal realizer. The
driver may call these staged owners, terminal products, arbitrary-center
nuclear attraction, terminal inventory, and
`_cartesian_terminal_due_diligence_report` by qualified private names. It must
not add a source helper or route around the public facade for ordinary charged
Hamiltonian construction.

Both methods use the same `ns=5`-derived physical parent. PQS uses route-local
`q=5`; White-Lindsey uses route-local `q=3`. The White-Lindsey `q` must not
define a different parent scale. Independently constructed method instances
must agree in parent centers, weights, mapped axis operators, physical bounds,
and fingerprints. The script must validate route identity and live objects
rather than infer them from labels.

For every eligible shared complete shell, both methods must also consume the
same aspect-selected outer source shape `(ns,ns,L)`. PQS retains its boundary
quotient; White-Lindsey uses axis-specific inner counts
`(ns-2,ns-2,L-2)`. Direct-core, aggregate complete-shell, slab, and bare final
dimensions must match. The driver does not choose this policy or patch a
shape; it reports and validates the source-backed construction.

Physical nuclear attraction must use the actual two centers through the
existing arbitrary-center kernel. No route may use
`pqs_source_shell_realization_final_basis`, CPBM source-shell bridge/readiness
paths, or shell/support-row compatibility materialization. The supplemented
preflight instead uses the existing supplemented working-basis seam exactly as
specified below; it remains the separately frozen neutral-H2 endpoint.

The frozen paper defaults are:

```text
public ns = 5
core_spacing = 1.2 / (Z * (ns - 1)) = 0.30
mapping_s_standard = sqrt(Z * core_spacing) = sqrt(0.30)
mapping_s_effective = mapping_s_standard
s_factor = 1.0
tail_spacing = 2.8 (default)
source_span = :ordinary
coulomb_accuracy = :high
Coulomb expansion terms = 135
PQS route-local q = 5
White-Lindsey route-local q = 3
```

Parent extents follow the physical two-center padding convention:

```text
xmax_parallel = R/2 + padding
xmax_transverse = padding
```

The accepted `R=2` baseline uses `padding=10.0` bohr. The completed bounded box
control used `padding=20.0` bohr. These are the only campaign-authorized
padding values; accepting a numeric `padding` input does not authorize a
general padding scan.

The completed finer-tail control used `tail_spacing=2.0`, only with `R=2.0`
and `padding=10.0`, and favored retaining `2.8`. The exact
campaign-authorized combinations remain:

```text
(padding=10.0, tail_spacing=2.8)
(padding=20.0, tail_spacing=2.8)
(padding=10.0, tail_spacing=2.0)
```

The driver must reject other combinations. Equivalently, both inputs must
belong to their listed two-value sets and
`tail_spacing == 2.8 || padding == 10.0`. In particular, the accepted
padding-20 control is evidence only at `tail_spacing=2.8`; the
`padding=20.0`, `tail_spacing=2.0` combination is not authorized.

The allowed controls remain only:

```text
system = :h2plus | :h2
method = :pqs | :wl | :both
R
output_dir
padding
tail_spacing
optional trusted config include
name=value overrides for these approved inputs
```

`system=:h2` is valid only for the frozen supplemented preflight with
`method=:both`, `R=2`, `padding=10.0`, and `tail_spacing=2.8`. The driver
records every resolved input. It is not a mapping, `q`, source-span,
supplement, Coulomb-policy, padding, tail-spacing, or general parameter
scanner.

### Mandatory Full-Parent Reference Row

At `R=2`, `method=:both` must emit three rows from the same live construction:

```text
method = parent
method = pqs
method = wl
```

`parent` is an internal output-row identity, not a new public input mode. The
accepted padding-10 row used a `21 x 21 x 29` (`12789`-function) PGDG parent;
that dimension is evidence, not a fixed parent-shape contract. Each authorized
padding must use the live constructed parent without forming its dense
three-dimensional matrix.

The live parent axis counts must be positive odd integers, their product must
equal the reported parent dimension, and every parent weight must be finite
and strictly positive. For each authorized padding/tail-spacing combination,
the independently constructed PQS and White-Lindsey routes must retain
identical live parent objects, physical bounds, and fingerprints. The resolved
tail spacing must agree with the live parent setup and every live mapping.
Dimensions and bounds must be reported from those live objects rather than
inferred from the padding-10 evidence.

The driver may consume the live PGDG axis overlap and kinetic matrices,
analytic high135 nuclear Gaussian factors at the two physical centers, and
the existing apply-based `_lanczos_ground_state_apply` routine. It must apply
the separable three-dimensional parent `H1` after axis-wise symmetric overlap
orthogonalization. Only small file-local mode-product functions needed for
that application are approved. No source helper, reusable parent-oracle
framework, generalized eigensolver, quadrature path, or new dependency is
allowed.

After the lowest state is returned, the driver must independently apply the
complete parent `H1` again and report the recomputed residual. It must also
apply and contract the kinetic, left-nuclear, and right-nuclear terms
separately, with

```text
H1_expectation = T_expectation + U_left_expectation + U_right_expectation
```

closing at the existing numerical tolerance. The second parent eigenvalue and
gap may be unavailable because the current apply-based solver returns only the
lowest state.

### Full-Parent State Capture

The implemented capture diagnostic is limited to `R=2`, `padding=10.0`,
`tail_spacing=2.8`, and `method=:both`. It measures how much of the normalized
full-parent ground state is retained by the frozen PQS `q=5` and
White-Lindsey `q=3` terminal spans. It changes neither terminal construction
nor terminal coefficients and introduces no new input mode.

Let `y` be the Euclidean-normalized eigenvector of the overlap-orthogonalized
parent problem, let `S` be the separable parent overlap, let
`X=S^(-1/2)`, and let `C` denote a terminal realization's support-local parent
coefficients. The diagnostic is

```text
c = X y
a = C' S c = C' S^(1/2) y
capture = a'a
r = y - S^(1/2) C a
lost_norm = r'r
```

The driver must evaluate `S^(1/2)` through the existing one-dimensional
factors `S_axis * X_axis`. It must apply `C` blockwise from each
`CartesianTerminalBasisBlock`'s support indices, coefficients, and native
column range. A dense parent-by-terminal coefficient, overlap, or projection
matrix is forbidden.

Before interpreting capture, the driver must establish that terminal block
supports form a disjoint exact partition of `1:parent_dimension`, that
due-diligence rows match blocks in native order and exact column ranges, and
that both methods use the same parent fingerprint and the same parent
eigenvector. It must then require:

- finite capture fraction within `[0,1]` up to numerical tolerance;
- the parent-space residual to be orthogonal to the terminal span;
- captured norm plus lost norm to equal the parent-state norm;
- nonnegative regional lost norms that sum to the global lost norm; and
- regional captured contributions that sum to the global captured norm.

For a block, regional loss is `sum(abs2, r[block.support_indices])`, and
regional capture is the squared norm of `a` over the block's column range.
The readable report aggregates these values by physical `region_key`,
`region_kind`, and `shell_index`; slab entries additionally retain axis, side,
stack index, and stack count. It must not define local loss as source norm
minus captured norm, because symmetric parent orthogonalization does not make
that a unique local decomposition.

The existing TSV and readable report may add exactly:

```text
parent_ground_state_norm
terminal_capture_fraction
terminal_lost_norm
capture_closure_error
```

`parent_ground_state_norm` is reported for all three rows. The terminal-only
fields are not applicable on the parent row. Variable-length regional
breakdowns belong only in the readable report; no additional output file,
payload, or schema family is approved. PQS and White-Lindsey capture are
reported without an ordering gate. The White-Lindsey
`source_mode_shape`-unavailable advisory remains truthful and must not be
replaced with a fabricated PQS-only field.

### Supplemented H2 One-Body Preflight

One private supplemented endpoint is implemented only for:

```text
system = :h2
method = :both
R = 2.0
padding = 10.0
tail_spacing = 2.8
nup = 1
ndn = 1
```

The neutral electron metadata is required by the existing supplemented
producer. The endpoint stops after exact one-body assembly, so its operator is
the same two-proton operator as the accepted H2+ rows. It must be described as
a neutral-metadata supplemented one-body preflight, not as a charged
supplemented H2+ production calculation.

The run preserves the accepted bare `parent`, `pqs`, and `wl` rows and adds
exactly `pqs_supplemented` and `wl_supplemented` rows. Bare and supplemented
rows must use the same physical parent settings and live parent fingerprint.
The supplemented rows must use the existing neutral
`cartesian_base_working_basis(...; supplemented=true)` path and existing
residual/augmented-operator owners. They must not reconstruct the terminal
basis, residual algorithm, or exact operators in the driver.

The supplement is frozen and has no new user controls:

```text
basis_by_center = ["cc-pVTZ", "cc-pVTZ"]
lmax = 1
uncontracted = false
width_filtering = nothing
basisfile = nothing
residual_occupation_cutoff = 1.0e-6
```

Omitting `basisfile` means the bundled trusted legacy `BasisSets` data. The
same ordered raw supplement must be used for both routes. Validate and report
its labels, Cartesian angular powers, centers, owner assignment, primitive
exponents and contraction coefficients, exact self-overlap, conditioning, and
one deterministic fingerprint. Historical candidate or retained counts are
not acceptance constants.

Let `A` be the ordered raw supplement, `S_AA=<A|A>`,
`B=<parent|A>`, and `X=<G|A>` for the bare terminal basis. After validating
positive supplement rank without eigenvalue flooring, report the metric-aware
capture operators

```text
K_parent = S_AA^(-1/2) B' S_parent^(-1) B S_AA^(-1/2)
K_terminal = S_AA^(-1/2) X' X S_AA^(-1/2)
```

and their singular/capture spectra. Apply `S_parent^(-1)` through the existing
separable axis factors. A bounded parent-by-supplement cross overlap from the
existing non-nuclear raw-block owner is allowed; a dense parent overlap,
generalized parent-plus-supplement solve, or new parent augmentation is not.

For each physical owner, independently inspect

```text
M_a = S_AA[a,a] - X[:,a]' X[:,a]
```

and require its retained rank and occupations to match the existing
owner-local residual builder at the production cutoff `1.0e-6`. Report the
complete owner spectra, minimum retained occupation, maximum discarded
occupation, and cutoff margins. A numerically marginal direction is a stop
condition; it is not permission to switch to the numerical-complete
`1.0e-10` policy, floor an eigenvalue, or inject a direction.

Build the exact augmented kinetic and by-center unit-nuclear matrices only
through:

```text
cartesian_residual_gto_supplement_basis
cartesian_residual_gto_augmentation
cartesian_residual_gto_augmented_products
cartesian_residual_gto_augmented_unit_nuclear
```

The trusted bare kinetic and unit-nuclear blocks may be passed through their
existing keywords. Assemble the two-proton `H1` directly from those exact
matrices. Do not call `cartesian_base_vee`,
`cartesian_residual_gto_augmented_vee`, or a Hamiltonian assembly owner.
Require finite symmetric kinetic, both nuclear matrices, and `H1`; bare
`G-G` block parity; terminal-residual orthogonality; augmented-overlap
identity; lowest-state residual; and a variational augmented one-body energy
relative to the corresponding bare terminal row.

The existing private TSV/report may add exactly:

```text
supplement_fingerprint_sha256
supplement_candidate_count
supplement_parent_capture_min_sv
supplement_terminal_capture_min_sv
residual_occupation_cutoff
residual_dimension
residual_min_retained_occupation
residual_max_discarded_occupation
terminal_residual_orthogonality_error
bare_terminal_energy_change_Ha
```

Use the existing `final_dimension`, `overlap_identity_error`, one-body
decomposition, symmetry, residual, timing, allocation, RSS, and path fields.
Nonapplicable rows must say so explicitly. Complete raw-supplement identity,
capture spectra, owner spectra/ranks/margins, base/residual/augmented
dimensions, operator dimensions/fingerprints, and due diligence belong in the
readable report only. Do not dump full matrices or add an output file, payload,
artifact, or schema family.

The preflight passes only if both routes use the same parent and exact raw
supplement; loading and ownership close; capture spectra are finite and
physical; production residual selection is reproduced; all metric, symmetry,
identity, residual, and variational gates pass; and the report cleanly
separates parent capture, terminal contraction, and external residual content.
The accepted one-body rows and diagnostics remain frozen under the bounded
interaction extension below. Any other `Vee`, IDA/MWG,
RHF/SCF, cutoff control, or supplement variation requires later authority.

### Fixed-State H2 Density-Density Interaction

The implemented private endpoint is limited to the same neutral-metadata
five-row construction:

```text
system = :h2
method = :both
R = 2.0
padding = 10.0
tail_spacing = 2.8
```

It introduces no new input. The common physical target is the already accepted
full-parent H2+ ground-state spatial orbital. The driver must project that
single target independently into the full parent, bare PQS, bare
White-Lindsey, supplemented PQS, and supplemented White-Lindsey spaces,
normalize each projection, and hold every orbital fixed. It must not use the
lowest neutral-H2 orbital already reported by the one-body preflight and must
not perform an SCF, orbital update, or variational interaction solve.

Let `y` be the normalized overlap-orthogonalized parent state,
`c_P=S_parent^(-1/2)y` its physical parent coefficients, `C_G` the terminal
coefficient map, `B=<parent|A>`, and `(T_G,T_A)` the existing residual
transform. The unnormalized projection coefficients are:

```text
g = C_G' S_parent^(1/2) y
a = B' c_P
r = T_G' g + T_A' a

u_parent = c_P
u_bare = g / sqrt(g'g)
u_supplemented = [g; r] / sqrt(g'g + r'r)
```

The driver must report `g'g` and `g'g+r'r` as target captures and validate
normalization, finite values, dimensions, native ordering, and nondecreasing
capture under augmentation. It must use the actual accepted `B`, `T_G`, and
`T_A`; reconstructing a second residual basis or appending target directions
is forbidden. Parent normalization uses the physical `S_parent` metric;
terminal and supplemented normalization uses the Euclidean metric of their
orthonormal native bases.

Each normalized state is evaluated with the exact one-body matrix already
owned by its row and the category-appropriate production density-density
interaction:

| Row | Interaction owner |
| --- | --- |
| full parent | Matrix-free ordinary parent IDA using the live high135 separable PGDG pair factors |
| bare PQS / White-Lindsey | `cartesian_base_vee` terminal/site IDA |
| supplemented PQS / White-Lindsey | `cartesian_residual_gto_augmented_vee`, retaining terminal IDA in `G-G` and existing MWG in `G-R` and `R-R` |

The parent contraction must use
`q_P=abs2(c_P)/sum(abs2(c_P))` and apply the existing high135 separable pair
factors term by term. It must not form a dense
`parent_dimension x parent_dimension` matrix. Bare and augmented matrices must
come from the existing base/Hamiltonian and residual-MWG owners. The augmented
`G-G` block must equal the matching freshly built bare terminal IDA block; no
`C' V C`, interaction rotation, fitted-density substitute, or reused
interaction from another basis is allowed.

For a normalized real spatial orbital `u` occupied once in each spin channel,
let `q=abs2(u)` and `J=q'Vq`. Report the current density-density accounting
explicitly:

```text
direct term = 2J
same-spin exchange/self-cancellation contribution = -J
Vee expectation = J
fixed electronic energy = 2<u|H1|u> + J
fixed total energy = fixed electronic energy + 1/R
```

The direct and same-spin terms are an accounting decomposition of the current
density-density Hamiltonian. This endpoint does not claim continuum-exact
exchange or add a transition-density interaction.

The original five one-body rows and every pre-existing field must remain
numerically unchanged. The same private TSV/report may add exactly:

```text
fixed_target_capture
fixed_target_norm_error
fixed_target_fingerprint_sha256
fixed_H1_expectation_Ha
fixed_density_density_direct_Ha
fixed_density_density_exchange_Ha
fixed_density_density_Vee_Ha
fixed_electronic_energy_Ha
fixed_total_energy_Ha
interaction_convention
interaction_symmetry_error_Ha
interaction_elapsed_s
interaction_allocated_bytes
```

The state fingerprint is over the sign-canonicalized normalized coefficient
vector in that row's native basis and binds the later independent replay. Fix
the sign by making the first largest-magnitude coefficient positive before
hashing the native `Float64` bytes.
Complete interaction dimensions, base/augmented block parity, matrix or
factor-resource fingerprints, charge closure, finite/symmetry checks, and
phase timing belong in the readable report. Matrix rows use the matrix
infinity-norm symmetry error; the matrix-free parent uses a deterministic
bilinear apply-symmetry check. No coefficient dump, additional output file,
artifact, payload, or schema family is approved.

The production result is a measurement awaiting independent validation. A
paper-validation script outside this repository must reconstruct each state
from the exact clean commit, match parent/basis/supplement fingerprints,
captures, normalized-state fingerprints, and fixed one-body expectations, and
then converge a direct-density Coulomb oracle on that same reconstructed
density. Compare that oracle scalar with
`fixed_density_density_Vee_Ha=J`, not with the intermediate `2J` direct term.
Numerical integration, density fitting, and the oracle implementation must not
enter `GaussletBases`. No production interaction error or method ordering may
be interpreted until this same-density replay passes.

Implementation must stop without a commit if the exact projected states cannot
be reconstructed from the existing parent/supplement/residual objects, if the
augmented interaction requires a new source seam or state carrier, if a
category block differs from its existing owner, or if any capture,
normalization, charge, symmetry, finiteness, decomposition, fingerprint, or
unchanged-one-body gate fails.

### Endpoints And Output

The implemented endpoints are the three-row lowest-`H1` H2+ comparison, the
five-row bare-plus-supplemented one-body preflight, and the five-row
fixed-state density-density measurement above. They build no RHF/SCF state or
Hamiltonian artifact.

Each method writes one compact TSV row and a readable text report containing
the resolved configuration, git commit and dirty state, dimensions, shell
ledger, energy decomposition, convergence result, phase timings and
allocations, process peak RSS with units, and output paths. The readable report
must include the complete existing structured terminal due-diligence object
for both methods, including full parent-axis, dimension, shell/slab-row, and
warning content. A summary or driver reconstruction is not a substitute. If
the White-Lindsey route cannot supply complete shell/slab rows through the
existing owner, implementation must stop and report that exact source-backed
gap.

The parent row has no terminal ledger or terminal due-diligence rows; those
fields are explicitly not applicable. PQS and White-Lindsey retain complete
terminal due diligence and exact column accounting. All three rows must carry
one clean Git commit and the same live parent fingerprint.

The existing private TSV/report may add:

```text
independent_reference_error
parent_resolution_error
contraction_error
```

using the frozen `R=2` H2+ total-energy reference
`-0.6026342144949465 Ha`. Define

```text
independent_reference_error(row) = E_total(row) - E_reference
parent_resolution_error = E_total(parent) - E_reference
contraction_error(method) = E_total(method) - E_total(parent)
```

The parent contraction error and parent terminal-ledger fields may be marked
not applicable. This is an extension of the existing private report, not a new
schema, payload, metadata family, or report-carried basis.

### Implementation And Acceptance Limits

- The accepted capture-enabled driver was `299` lines. The implemented
  supplemented preflight is `416` lines. Reusing its live construction for the
  fixed-state interaction is preferred over a companion that would duplicate
  parent, supplement, and residual assembly. Prefer a final size at or below
  `500` lines and require at most `525` lines.
- Added `src` lines, committed tests, probes, fixtures, modules, helper files,
  exports, status vocabularies, adapters, and companion drivers are all zero.
- Do not copy the scratch parent oracle, sampled-density oracle,
  generalized-overlap solver, AddNest include chain, or scratch reporting
  framework.
- The H2+ implementation pass changed the exact authority path state from
  `planned` to `existing` and regenerated the checked views in the same commit.
- If the endpoint cannot be implemented in the one approved file and hard
  budget using existing kernels, make no implementation commit and report the
  missing operation.

The initial three-row scientific gate was accepted with:

```text
julia --project=. bin/pqs_paper_h2_driver.jl \
    system=:h2plus method=:both R=2.0 padding=10.0 \
    output_dir='"/tmp/pqs_paper_h2plus_R2"'
```

The accepted padding convergence control was:

```text
julia --project=. bin/pqs_paper_h2_driver.jl \
    system=:h2plus method=:both R=2.0 padding=20.0 \
    output_dir='"/tmp/pqs_paper_h2plus_R2_full_parent_padding20_clean"'
```

It must report the parent/PQS/White-Lindsey rows from one clean commit and
shared live parent fingerprint. PQS and White-Lindsey retain their route-local
`q=5` and `q=3`, finite symmetric `H1`, converged lowest eigenpairs, and
complete due-diligence/output records. The parent row must pass the matrix-free
Lanczos, independent residual, decomposition, and frozen-reference gates.
Compare each padding-20 energy with its preserved padding-10 counterpart using
a provisional absolute shift tolerance of `0.01 mHa`. Report the live
dimensions, topology, bounds, actual padding, residuals, timing, allocations,
and peak RSS for the control. No PQS-versus-White-Lindsey energy ordering is an
acceptance criterion.

The one authorized finer-tail control is:

```text
julia --project=. bin/pqs_paper_h2_driver.jl \
    system=:h2plus method=:both R=2.0 padding=10.0 tail_spacing=2.0 \
    output_dir='"/tmp/pqs_paper_h2plus_R2_tail2_clean"'
```

It must build the same matrix-free parent/PQS/White-Lindsey object with
`ns=5`, `core_spacing=0.30`, `s_factor=1`, `source_span=:ordinary`, high135,
and route-local `q=5/q=3`. Compare its three energies, reference and
contraction errors, axes, dimensions, bounds, actual padding, shell/slab
topology, warnings, residuals, timing, allocations, and peak RSS against the
preserved `tail_spacing=2.8`, padding-10 rows. Report the shifts without
imposing an energy ordering or convergence threshold.

The authorized capture diagnostic uses the accepted baseline command:

```text
julia --project=. bin/pqs_paper_h2_driver.jl \
    system=:h2plus method=:both R=2.0 padding=10.0 tail_spacing=2.8 \
    output_dir='"/tmp/pqs_paper_h2plus_R2_capture_clean"'
```

It must preserve the complete PQS/White-Lindsey due-diligence records and
report the four fixed capture fields plus the readable regional breakdown.
The capture implementation is accepted only if all support, native-order,
orthogonality, global-closure, and regional-closure gates above pass.

The implemented supplemented preflight command is:

```text
julia --project=. bin/pqs_paper_h2_driver.jl \
    system=:h2 method=:both R=2.0 padding=10.0 tail_spacing=2.8 \
    output_dir='"/tmp/pqs_paper_h2_supplemented_one_body_clean"'
```

Run it from a clean worktree. It must preserve the accepted bare rows, append
the two supplemented rows, include complete terminal due diligence for both
routes, and satisfy every supplement, capture, residual-rank, exact one-body,
and variational gate above.

Repo acceptance at `119bc17c7` retained the exact bare dimensions
`12789/1285/1109`, added `18` residual directions to each terminal route, and
produced supplemented dimensions `1303/1127`. The minimum parent supplement
capture singular value was `0.9991415015`; terminal minima were
`0.9781753530/0.7793743583` for PQS/WL. Both owner-local residual ranks were
`9+9`, with no discarded direction at the production cutoff. This is
repository acceptance of the bounded one-body endpoint, not independent paper
acceptance; that commit alone did not establish the interaction stage below.
The `1285/1109` terminal difference was subsequently traced to the old
PQS-only aspect policy: common physical shells retained `960/784` PQS/WL
complete-shell functions. These rows remain valid old-policy measurements but
are provisional for the matched paper comparison.

The implemented fixed-state interaction command is the same frozen H2 command:

```text
julia --project=. bin/pqs_paper_h2_driver.jl \
    system=:h2 method=:both R=2.0 padding=10.0 tail_spacing=2.8 \
    output_dir='"/tmp/pqs_paper_h2_fixed_interaction_clean"'
```

Run it from a clean worktree. The five original one-body rows must match the
accepted replay before the new interaction fields are considered. The repo
gate establishes construction, projection, category ownership, and internal
energy accounting only. Independent same-density oracle convergence is a
separate paper-validation acceptance gate.

Repository acceptance at `7d2b6dc61` retained dimensions
`12789/1285/1109/1303/1127` and target captures
`1.0/0.9999709356/0.9991007021/0.9999972239/0.9999933045`. The corresponding
production density-density `J` values are
`0.6613194326/0.6614182833/0.6632012296/0.6614042003/0.6625621285 Ha`.
These are construction-backed measurements, not an interaction-accuracy
ordering. They must be replayed after the matched shared-shell correction;
the external same-density oracle then uses the corrected WL state
fingerprints. Old WL values are evidence, not parity targets.

This authority changes no public input, canonical-driver behavior, numerical
default, artifact schema, AddNest/exact-W/injection/PRF/external-RG/MWG/
screening/EGOI behavior, or correlated-solver surface. H2 beyond this exact
fixed-state gate, RHF/SCF, orbital relaxation, continuum-exact interaction
claims, transition-density exchange, numerical integration in the repo,
capture at any other campaign point, `q` ladders, general padding or
tail-spacing scans, geometry curves, enrichment, supplement variation, PRFs,
test changes, additional output columns beyond the exact lists above, and new
artifacts remain forbidden.
