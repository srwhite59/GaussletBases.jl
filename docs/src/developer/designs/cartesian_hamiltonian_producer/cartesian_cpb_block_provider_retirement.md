# Cartesian CPB Block Provider Retirement

## Status And Authority

This page owns the retirement contract for:

- `HP-RETIRE-CPB-PROVIDER-FN-01`
- `HP-RETIRE-CPB-PROVIDER-TEST-01`

Status: source retirement is approved. Deletion and validation require a
separate reviewed implementation pass.

The retirement removes an orphaned experimental provider pilot. It does not
retire Cartesian coordinate-product-box geometry, parent-axis factors,
Gaussian raw-block kernels, pair materialization, terminal realization, or any
current Hamiltonian producer.

## Retirement Decision

At baseline `34f4c656a1f64a5285fa0bac3291b2288c81d444`, delete:

- `src/CartesianCPBBlockProviders.jl`, in full;
- only `include("CartesianCPBBlockProviders.jl")` from
  `src/GaussletBases.jl`.

The module contains `6,445` lines and the include contributes one line. The
required source reduction is therefore exactly `6,446` lines, with zero added
source lines.

Do not add an alias, stub, deprecation, compatibility module, replacement
provider, helper, or test. Qualified access to this internal, unadvertised
submodule is not a compatibility obligation.

## Retired Surface

The module's complete 47-name export surface retires with the file:

```text
CPBIntervalPair3D
CPBOverlapAxisBlockSet
CPBAxisProductOperatorBlock
CPBSumOfAxisProductsOperatorBlock
CPBOverlapOperatorBlock
CPBKineticOperatorBlock
CPBOneBodyAxisOperatorBlock
CPBMixedGTOLocalOverlapBlock
CPBMixedGTOCoordinateMomentLocalBlock
CPBMixedGTOKineticLocalBlock
CPBMixedGTOSupplementLocalBlock
CPBGTOSupplementOneBodyBlock
CPBMixedGTOSupplementNuclearByCenterBlock
CPBGTOSupplementNuclearByCenterBlock
CPBGTOSupplementLocalOperatorBundle
CPBElectronElectronLocalBlock
CPBElectronNuclearByCenterLocalBlock
CPBLocalIntegralWeights
CPBOverlapDenseBlock
CPBLocalOverlapBlockRecord
CPBLocalOverlapBlockCollection
cpb_interval_pair
cpb_overlap_axis_blocks
cpb_axis_product_operator_block
cpb_sum_of_axis_products_operator_block
cpb_overlap_operator_block
cpb_kinetic_operator_block
cpb_position_operator_block
cpb_x2_operator_block
cpb_mixed_gto_overlap_block
cpb_mixed_gto_position_operator_block
cpb_mixed_gto_x2_operator_block
cpb_mixed_gto_kinetic_operator_block
cpb_gto_overlap_operator_block
cpb_gto_position_operator_block
cpb_gto_x2_operator_block
cpb_gto_kinetic_operator_block
cpb_mixed_gto_nuclear_by_center_block
cpb_gto_nuclear_by_center_block
cpb_gto_supplement_local_operator_bundle
cpb_electron_electron_local_block
cpb_electron_nuclear_by_center_local_block
cpb_local_integral_weights
cpb_overlap_dense_block
cpb_local_overlap_block_record
cpb_local_overlap_block_collection
summary
```

## Caller And Side-Effect Audit

At the approved baseline:

- the only committed reference outside the file is its root include;
- no committed source, test, driver, tool, example, or docs build path calls
  the module;
- no root export exposes the submodule or its names;
- the module owns no artifact, serialization, registration, or load-time side
  effect.

Three completed July 2026 pure-GTO runners under external `work/cr2` use the
qualified module:

- `2026-07-07_cr2_sto6g_pure_gto_uhf`;
- `2026-07-07_cr2_5zs_sto_pd_pure_gto_uhf`;
- `2026-07-07_cr_5zs_sto_pd_pure_gto_uhf`.

They are historical probes and may remain pinned to a pre-retirement commit.
They do not justify compatibility code in current source.

If a live committed caller or load-time side effect is found before deletion,
stop without committing and report it. Do not narrow the retirement by adding
glue.

## Preserved Owners

The retirement must not delete, move, or reinterpret:

- `src/cartesian_cpb/` coordinate-product-box geometry;
- `src/CartesianParentGaussletBases.jl` and
  `src/CartesianParentAxisFactors.jl`;
- `src/GaussianAnalyticIntegrals.jl`, `src/ordinary_coulomb.jl`, and
  `src/cartesian_gaussian_raw_blocks/`;
- Cartesian terminal, PQS, White--Lindsey, PRF, residual-GTO, and
  represented-Hartree owners;
- `src/cartesian_pair_block_materialization/`;
- contracted-parent and multilayer source surfaces.

The old CPB overlap-placement ladder and CPB-provider destinations in cleanup
notes become historical evidence only. They grant no replacement or deletion
authority over the preserved owners.

## Validation Contract

The source retirement pass must:

1. remove the file, include, module name, and all 47 exported names from
   committed source;
2. confirm package load;
3. run the unchanged `core`, `ida`, and `cartesian` test groups;
4. inspect and report terminal due diligence for the Cartesian endpoint;
5. confirm exactly `6,446` deleted and zero added source lines;
6. run `git diff --check` and the mechanical anti-bloat scan.

Tests are validation inputs, not edit surfaces for this retirement. If any test
change is required, stop and request a separate amendment. No replacement test
is justified for a deleted API with no committed consumer.

## Non-Goals

This authority does not permit:

- pair-materialization, contracted-parent, or multilayer deletion;
- changes to CPB geometry, parent factors, raw Gaussian blocks, terminal
  realization, PQS/WL, PRFs, residual Gaussians, represented Hartree, IDA, MWG,
  artifacts, drivers, or solvers;
- a new provider abstraction or revival of the retired placement ladder;
- changes to external historical runners.
