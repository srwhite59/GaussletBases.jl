# High-Order Endcap/Panel H2 Chemistry Reproduction

Date: 2026-05-16

This memo records a mainline scratch reproduction of the old-standard-aligned
H2 HF/ED chemistry comparison for the guarded endcap/panel shared-shell path.
The route is still internal, experimental, and disabled by default.

## Question

Can the mainline guarded endcap/panel path reproduce the high-order R = 4.0
old-standard H2 chemistry row after the source has been imported through
ordinary/QW operator construction?

The comparison uses:

- H2 at `R = 4.0` bohr
- `core_spacing = 0.5`
- `xmax_parallel = 6.0`
- `xmax_transverse = 4.0`
- H/cc-pVTZ molecular Gaussian supplement with `lmax = 1`
- `nside = 5`
- finite IDA/QW density-density Hamiltonian
- restricted closed-shell HF
- `Diag2ptle.davidson_ground_sylv` two-electron ED

The chemistry reproduction intentionally uses `interaction_treatment =
:ggt_nearest`, because that is the interaction route used by the high-order
old-standard rows being reproduced. MWG remains the preferred residual-Gaussian
interaction route for current operator preflight, but it was not used for this
old-standard chemistry match.

## Bottom Line

The mainline guarded endcap/panel route reproduces the high-order R = 4.0
old-standard rows to numerical precision.

| route | fixed block size | final dim | residuals | HF total | ED total | BO error (mHa) | ED-HF lowering (mHa) |
|---|---:|---:|---:|---:|---:|---:|---:|
| default complete rectangular | `(1215, 463)` | 481 | 18 | -0.910938264352 | -1.015613837691 | 0.776415 | 104.675573 |
| endcap/panel q=4,L=4 | `(1215, 443)` | 461 | 18 | -0.910977315003 | -1.015663743783 | 0.726509 | 104.686429 |

The endcap/panel route saves 20 final functions and lowers the ED total by
`0.049906` mHa relative to the default complete-rectangular route.

Against the high-order reference rows:

- mainline default minus high-order old S/P ED: `4.9e-11` mHa
- mainline endcap minus high-order all-shared q=4,L=4 ED: `1.25e-10` mHa

This is an exact reproduction for practical purposes.

## Operator and Solver Preflight

Both routes passed the finite-Hamiltonian checks used by the scratch script.

Default complete-rectangular route:

- shared layer types: `_CartesianNestedCompleteShell3D`,
  `_CartesianNestedCompleteShell3D`
- shared layer columns: `(98, 114)`
- fixed overlap error: `1.38e-14`
- operator overlap error: `2.44e-12`
- H symmetry error: `0.0`
- V symmetry error: `0.0`
- residual width status: `all_nan` under `:ggt_nearest`
- residual owner set: `(1, 2)`
- HF converged in 9 iterations
- ED converged by tolerance in 20 iterations with residual `2.99e-11`

Endcap/panel q=4,L=4 route:

- shared layer types: `_CartesianNestedEndcapPanelShellLayer3D`,
  `_CartesianNestedEndcapPanelShellLayer3D`
- shared layer columns: `(96, 96)`
- fixed overlap error: `1.38e-14`
- operator overlap error: `1.03e-12`
- H symmetry error: `0.0`
- V symmetry error: `0.0`
- residual width status: `all_nan` under `:ggt_nearest`
- residual owner set: `(1, 2)`
- HF converged in 9 iterations
- ED converged by tolerance in 19 iterations with residual `4.58e-11`

## Timing and Allocation

The scratch run used the existing `Diag2ptle` module at
`/Users/srwhite/Dropbox/GaussletModules/Diag2ptle.jl`.

Main stage timings from the fresh run:

| stage | time (s) | allocation |
|---|---:|---:|
| Coulomb expansion | 0.064513 | 17.4 MB |
| bond-aligned basis | 0.025996 | 1.7 MB |
| H cc-pVTZ S/P supplement | 0.004835 | 0.4 MB |
| default source and fixed block | 12.362971 | 1.16 GB |
| default ordinary QW operators | 12.587400 | 1.79 GB |
| default HF | 0.378750 | 322.8 MB |
| default ED | 0.310978 | 226.6 MB |
| endcap source and fixed block | 0.459292 | 150.1 MB |
| endcap ordinary QW operators | 7.777075 | 874.0 MB |
| endcap HF | 0.228844 | 297.3 MB |
| endcap ED | 0.241863 | 200.0 MB |

The default source/fixed-block stage is much slower because it pays the first
large compilation/construction cost in this scratch process. The endcap route
still shows the expected smaller operator construction and ED dimensions.

## Reproduction Artifact

The scratch script is:

```sh
julia --project=. tmp/work/h2_endcap_panel_hf_ed_chemistry_reproduction.jl
```

It prints fixed-block dimensions, final dimensions, residual counts, HF totals,
ED totals, BO errors, ED-HF lowering, overlap and symmetry checks, solver
metadata, and timing/allocation rows.

## Boundary

This is a mainline experimental chemistry validation artifact only.

The current mainline frontend exposes the path only as an opt-in experimental
keyword, `shared_shell_layer_policy = :endcap_panel_owned`; the default remains
`:complete_rectangular`. No CR2 scripts were edited. No broad OPCU, FSB, or FBU
machinery was imported. No high-order branch merge was done.

The current evidence says the internal guarded path is now chemistry-ready for
repo-manager review at the old-standard H2 R = 4.0 level. The next decision is
whether repo-manager wants a public/internal frontend control, a broader stretch
sweep, or a separate MWG chemistry comparison.
