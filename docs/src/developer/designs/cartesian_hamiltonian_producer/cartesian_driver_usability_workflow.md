# Cartesian Driver Usability Workflow

Status: implemented canonical contract for the human-facing Cartesian
Hamiltonian driver and its non-exported staged producer calls, plus one
approved-pending private PQS/WL paper-validation driver. Registry entries own
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
All detailed geometry, neutrality, spacing, and policy checks remain producer
owned.

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
The producer then rejects malformed records, nonneutral systems, unsupported
geometry, invalid sizes/spacings/mapping factors, incompatible source spans,
and invalid supplement inputs. Input evaluation, numerical construction,
artifact writing, and readback errors propagate. Failures are not converted
to status payloads or partial-success objects.

## Non-Goals

This contract does not add route-stage controls, stop points, raw-provider
switches, diagnostic dumps, allocation probes, custom Coulomb parameters,
solver/RHF/ECP/EGOI behavior, Cr2-specific branches, public exports, result
wrappers, artifact schemas, committed input fixtures, or a safe untrusted-code
parser. Route and ladder diagnostics remain outside the canonical driver.

## Private PQS/WL Paper H2 Driver

`HP-PQS-PAPER-H2-DRV-FN-01` approves one tracked, non-exported scientific
validation script:

```text
bin/pqs_paper_h2_driver.jl
```

This script is private to the matched White-Lindsey/PQS H2+/H2 paper campaign.
It is not a public API, a general paper-driver framework, a replacement for
`cartesian_ham_builder.jl`, or a source of producer defaults. It may be removed
after the campaign if no durable use remains.

### Required Construction

For each requested method, the script builds the same frozen one-center
neutral-H representative-parent construction as the contraction template.
Independently built PQS and White-Lindsey parent instances must agree in
centers, weights, mapped axis operators, bounds, and fingerprints. The PQS
path must be the current source-box-first chain:

```text
cartesian_base_working_basis
-> cartesian_transforms
-> pqs_terminal_basis_realization
```

The White-Lindsey case must use the corresponding current terminal realizer
through the same staged working-basis boundary. The driver may call these and
the required operator/interaction routines by qualified non-exported names.
It must validate the resulting method, route-local `q`, parent, shell ledger,
and terminal realization rather than infer route identity from labels.

The representative parent is common to the two physical centers. Physical
nuclei are placed at `(0,0,-R/2)` and `(0,0,R/2)` only through the existing
arbitrary-center terminal unit-nuclear kernel. The script must not use the
public multicenter diatomic facade,
`pqs_source_shell_realization_final_basis`, CPBM source-shell
bridge/readiness paths, or shell/support-row compatibility materialization.

The frozen paper defaults are:

```text
core_spacing = 0.15
parent mapping d = core_spacing
s_factor = 1.0
s = sqrt(core_spacing)
tail_spacing = 2.8
ns = 5
source_span = :ordinary
coulomb_accuracy = :high
Coulomb expansion terms = 135
```

The first `R=2` gate uses representative-parent radius `6.0`; a later curve
may use one explicitly recorded campaign-fixed extent. The allowed controls
are only:

```text
system = :h2plus | :h2
method = :pqs | :wl | :both
R
output_dir
parent_radius (campaign-fixed for a curve)
optional trusted config include
name=value overrides for these approved inputs
```

The driver records every resolved input. It is not a mapping, `q`,
source-span, supplement, Coulomb-policy, or parameter scanner.

### Endpoints And Output

The first endpoint is the lowest terminal `H1` eigenpair for H2+. It builds no
`Vee` and runs no RHF. Only after that endpoint is accepted may the same script
build ordinary 135-term terminal IDA for H2, evaluate the fixed state formed
from the H2+ orbital, and run closed-shell density-density RHF through the
existing matrix contractions.

Each method writes one compact TSV row and a readable text report containing
the resolved configuration, git commit and dirty state, dimensions, shell
ledger, energy decomposition, convergence result, phase timings and
allocations, process peak RSS with units, and output paths. H2 may optionally
write and read the existing minimal Hamiltonian artifact. The one-center
template provenance must not be attached as if it described the physical
two-center geometry. No new JLD2 schema, result wrapper, payload, metadata
family, or report-carried basis is approved.

### Implementation And Acceptance Limits

- Preferred implementation size is at most `125` added `bin` lines; the hard
  limit is `150`.
- Added `src` lines, committed tests, probes, fixtures, modules, helper files,
  exports, status vocabularies, and adapters are all zero.
- Do not copy the scratch parent oracle, sampled-density oracle,
  generalized-overlap solver, AddNest include chain, or scratch reporting
  framework.
- The implementation pass must change the exact authority path state from
  `planned` to `existing` and regenerate the checked views in the same commit.
- If the endpoint cannot be implemented in the one approved file and hard
  budget using existing kernels, make no implementation commit and report the
  missing operation.

The sole initial scientific acceptance command is:

```text
julia --project=. bin/pqs_paper_h2_driver.jl \
    system=:h2plus method=:both R=2.0 \
    output_dir='"/tmp/pqs_paper_h2plus_R2"'
```

It must report a `25 x 25 x 25` parent (`15625` functions), PQS `q=5`, WL
`q=3`, ten `98`-column shells, a `125`-column core, final dimension `1105`
for both methods, a 135-term expansion, finite symmetric `H1`, a converged
lowest state, complete due-diligence/output records, and lower PQS than WL
variational error. A material mismatch stops the pass before any H2 endpoint.

This authority changes no public input, canonical-driver behavior, numerical
default, artifact schema, AddNest/exact-W/injection/PRF/external-RG/MWG/
screening/EGOI behavior, or correlated-solver surface.
