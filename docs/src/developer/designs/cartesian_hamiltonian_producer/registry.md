# Cartesian Hamiltonian Producer Authority Registry

> **Generated authority view. Do not edit.** The record-level source is
> [authority.toml](authority.toml), SHA-256 `13e7ec2bc3ce4e6909bcc8798c725a8d38d6377ed9663727295e1e16f76d34f9`.

Tracked producer work is authorized only when a unique record has an
execution grant and surface, and the requested change stays within its exact
owned paths, scope, `current.md`, `invariants.md`, and canonical contract.
Lifecycle never grants work by itself. Any missing or conflicting fact fails closed.

## Records

### HP-CGAI-FN-01 - historical Cartesian Gaussian axis-helper proposal

- **Lifecycle:** `superseded`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [cartesian\_gaussian\_raw\_blocks\_nuclear.md](cartesian_gaussian_raw_blocks_nuclear.md); heading `Cartesian Gaussian Raw Blocks - Nuclear Slice`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** `HP-CGRB-FN-02`, `HP-CGRB-NN-FN-01`
- **Scope:** none. Any future cross-owner in-place helper requires a new docs-only amendment rather than reactivation of this ID.

### HP-CGRB-FILE-01 - neutral Cartesian Gaussian raw-block module files

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_gaussian\_raw\_blocks\_nuclear.md](cartesian_gaussian_raw_blocks_nuclear.md); heading `Cartesian Gaussian Raw Blocks - Nuclear Slice`
- **Owned paths:**
  - `source` / `existing`: `src/GaussletBases.jl`
  - `source` / `existing`: `src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl`
  - `source` / `existing`: `src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl`
- **Evidence:**
  - `repo_path`: `test/core/runtests.jl`
  - `repo_path`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Dependencies:** none
- **Scope:** maintain the internal module and nuclear owner; no public export.

### HP-CGRB-FN-01 - exact uncharged Gaussian nuclear raw blocks

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_gaussian\_raw\_blocks\_nuclear.md](cartesian_gaussian_raw_blocks_nuclear.md); heading `Cartesian Gaussian Raw Blocks - Nuclear Slice`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl`
- **Evidence:** none
- **Dependencies:** `HP-CGRB-FILE-01`
- **Scope:** maintain exact uncharged by-center \`G-A\`/\`A-A\` output and the stable pairwise analytic factor formula.

### HP-CGRB-FN-02 - nuclear one-dimensional axis-family reuse

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_gaussian\_raw\_blocks\_nuclear.md](cartesian_gaussian_raw_blocks_nuclear.md); heading `Cartesian Gaussian Raw Blocks - Nuclear Slice`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl`
- **Evidence:** none
- **Dependencies:** `HP-CGRB-FN-01`
- **Scope:** maintain function-local family maps, unique center/family-pair tables, orientation flags, term-first filling, and coupled primitive products.

### HP-CGRB-NN-FILE-01 - non-nuclear raw-block file

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_gaussian\_raw\_blocks\_non\_nuclear.md](cartesian_gaussian_raw_blocks_non_nuclear.md); heading `Cartesian Gaussian Raw Blocks - Non-Nuclear Slice`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl`
  - `source` / `existing`: `src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl`
- **Evidence:**
  - `repo_path`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Dependencies:** none
- **Scope:** maintain the internal non-nuclear owner and module include.

### HP-CGRB-NN-FN-01 - exact non-nuclear Gaussian raw blocks

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_gaussian\_raw\_blocks\_non\_nuclear.md](cartesian_gaussian_raw_blocks_non_nuclear.md); heading `Cartesian Gaussian Raw Blocks - Non-Nuclear Slice`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_gaussian_axis_integrals.jl`
  - `source` / `existing`: `src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl`
- **Evidence:** none
- **Dependencies:** `HP-CGRB-NN-FILE-01`
- **Scope:** maintain exact \`G-A\`/\`A-A\` overlap, kinetic, x/y/z, and x2/y2/z2 blocks, symmetry, family reuse, and the overlap-only path.

### HP-CGRB-NN-TEST-01 - non-nuclear extraction validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_gaussian\_raw\_blocks\_non\_nuclear.md](cartesian_gaussian_raw_blocks_non_nuclear.md); heading `Cartesian Gaussian Raw Blocks - Non-Nuclear Slice`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain the H2 endpoint and existing indirect parity coverage. No dedicated raw-block test file currently exists.

### HP-CGRB-NN-WIRE-01 - Residual Gaussian and Qiu-White non-nuclear rewiring

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_gaussian\_raw\_blocks\_non\_nuclear.md](cartesian_gaussian_raw_blocks_non_nuclear.md); heading `Cartesian Gaussian Raw Blocks - Non-Nuclear Slice`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  - `source` / `existing`: `src/ordinary_qw_operator_assembly.jl`
  - `source` / `existing`: `src/ordinary_qw_raw_blocks.jl`
- **Evidence:**
  - `repo_path`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Dependencies:** `HP-CGRB-NN-FN-01`
- **Scope:** maintain the Residual Gaussian and main diatomic Qiu-White callers of the neutral blocks.

### HP-CGRB-TEST-01 - nuclear extraction validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_gaussian\_raw\_blocks\_nuclear.md](cartesian_gaussian_raw_blocks_nuclear.md); heading `Cartesian Gaussian Raw Blocks - Nuclear Slice`
- **Owned paths:**
  - `test` / `existing`: `test/core/runtests.jl`
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain the H2 endpoint and stable-factor oracle. No dedicated raw-block test file currently exists.

### HP-CGRB-WIRE-01 - Residual Gaussian and Qiu-White rewiring

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_gaussian\_raw\_blocks\_nuclear.md](cartesian_gaussian_raw_blocks_nuclear.md); heading `Cartesian Gaussian Raw Blocks - Nuclear Slice`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  - `source` / `existing`: `src/ordinary_qw_operator_assembly.jl`
  - `source` / `existing`: `src/ordinary_qw_raw_blocks.jl`
- **Evidence:** none
- **Dependencies:** `HP-CGRB-FN-01`
- **Scope:** maintain the two direct neutral-kernel call sites and indirect Qiu-White assembly compatibility.

### HP-CHANGE-01 - return shell overlap from existing shell plan — rejected/deferred

- **Lifecycle:** `rejected`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [terminal\_basis\_and\_base\_assembly.md](terminal_basis_and_base_assembly.md); heading `Terminal Basis And Base Assembly`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** rejected as standalone authority; returning shell overlap may exist only as a private implementation detail under \`HP-FN-00\` and creates no independent source surface.

### HP-COMP-ANGBOX-AUDIT-01 - angular-balanced shellification geometry audit

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [common\_terminal\_shell\_decomposition.md](common_terminal_shell_decomposition.md); heading `Common Terminal Shell Decomposition`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** angular-balanced shellification geometry audit.

### HP-COMP-ANGBOX-FN-01 - angular-balanced diatomic shellification

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [common\_terminal\_shell\_decomposition.md](common_terminal_shell_decomposition.md); heading `Common Terminal Shell Decomposition`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_shellification/terminal_geometry.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** emit native ordered \`:angular\_z\_extension\_slab\` stacks so the ordinary shell body plus planned axial extensions realizes the physical outer-nucleus angular target. It does not change real-shell retained policy or central-gap/contact ownership.

### HP-COMP-ANGBOX-TEST-01 - angular shellification validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [common\_terminal\_shell\_decomposition.md](common_terminal_shell_decomposition.md); heading `Common Terminal Shell Decomposition`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** none. The completed evidence does not authorize test maintenance or a committed Cr2 fixture.

### HP-COMP-ATOMBOX-FN-01 - one-center atom physical extent

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r1\_one\_center\_base\_atoms.md](r1_one_center_base_atoms.md); heading `R1 One-Center Base Atoms`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** derive atom parent counts from physical \`basis.radius\` and the existing mapping/spacing policy; \`ns\` remains resolution/nesting input.

### HP-COMP-ATOMBOX-TEST-01 - atom physical-extent validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r1\_one\_center\_base\_atoms.md](r1_one_center_base_atoms.md); heading `R1 One-Center Base Atoms`
- **Owned paths:**
  - `test` / `existing`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** validation maintenance.

### HP-COMP-BASEDIAT-FN-01 - base homonuclear z-axis diatomics

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [nesting\_supplement\_composition\_plan.md](nesting_supplement_composition_plan.md); heading `Nesting/Supplement Composition`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** validate explicit equal-symbol/equal-charge neutral all-electron diatomics at two finite distinct z-axis centers and send them through the existing PQS/WL base path.

### HP-COMP-BASEDIAT-TEST-01 - base diatomic validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [nesting\_supplement\_composition\_plan.md](nesting_supplement_composition_plan.md); heading `Nesting/Supplement Composition`
- **Owned paths:**
  - `test` / `existing`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** validation maintenance.

### HP-COMP-FACEPROD-FN-01 - neutral terminal face-product helper

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [common\_terminal\_shell\_decomposition.md](common_terminal_shell_decomposition.md); heading `Common Terminal Shell Decomposition`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/terminal_face_product_blocks.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** one route-neutral face/face-stack coefficient assembly over fixed normal-axis indices. It is not a new terminal-basis policy.

### HP-COMP-FACEPROD-TEST-01 - face-product validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [common\_terminal\_shell\_decomposition.md](common_terminal_shell_decomposition.md); heading `Common Terminal Shell Decomposition`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** none. The completed evidence does not authorize test maintenance or a dedicated committed fixture.

### HP-COMP-NS-FN-01 - public ns and derived q

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `driver`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [nesting\_supplement\_composition\_plan.md](nesting_supplement_composition_plan.md); heading `Nesting/Supplement Composition`
- **Owned paths:**
  - `driver` / `existing`: `bin/cartesian_ham_builder.jl`
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** normalize public \`ns\`, derive \`q = ns\` for PQS and \`q = ns - 2\` for WL, validate any legacy \`q\` compatibility, and record the existing compact provenance.

### HP-COMP-NS-TEST-01 - public ns validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [nesting\_supplement\_composition\_plan.md](nesting_supplement_composition_plan.md); heading `Nesting/Supplement Composition`
- **Owned paths:**
  - `test` / `existing`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** validation maintenance.

### HP-COMP-NSCORE-FN-01 - public ns direct-core side

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [public\_ns\_core\_side\_parity.md](public_ns_core_side_parity.md); heading `Public ns Direct-Core Side Parity`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/pqs_source_box_route_driver_helpers.jl`
- **Evidence:** none
- **Dependencies:** `HP-COMP-NS-FN-01`
- **Scope:** preserve \`direct\_core\_side = isodd(ns) ? ns : ns + 1\` for direct nucleus-centered identity blocks only. Boundary retained construction remains route-local.

### HP-COMP-NSCORE-TEST-01 - public ns direct-core validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [public\_ns\_core\_side\_parity.md](public_ns_core_side_parity.md); heading `Public ns Direct-Core Side Parity`
- **Owned paths:**
  - `test` / `existing`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** validation maintenance.

### HP-COMP-OUTERMM-FN-01 - outer-mismatch-only correction

- **Lifecycle:** `superseded`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [common\_terminal\_shell\_decomposition.md](common_terminal_shell_decomposition.md); heading `Common Terminal Shell Decomposition`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** outer-mismatch-only correction.

### HP-COMP-OUTERMM-TEST-01 - outer-mismatch-only validation

- **Lifecycle:** `superseded`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [common\_terminal\_shell\_decomposition.md](common_terminal_shell_decomposition.md); heading `Common Terminal Shell Decomposition`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** outer-mismatch-only validation.

### HP-COMP-SHELLGEOM-DIAT-FN-01 - diatomic common shellifier entry

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [common\_terminal\_shell\_decomposition.md](common_terminal_shell_decomposition.md); heading `Common Terminal Shell Decomposition`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_shellification/terminal_geometry.jl`
  - `source` / `existing`: `src/pqs_source_box_route_driver_helpers.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** feed the same public \`ns\`, direct-core side, centers, bond axis, and parent facts into common z-axis diatomic shellification before family lowering. Central-gap/contact redesign is not approved.

### HP-COMP-SHELLGEOM-DIAT-TEST-01 - diatomic shellifier-entry validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [common\_terminal\_shell\_decomposition.md](common_terminal_shell_decomposition.md); heading `Common Terminal Shell Decomposition`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** none. The completed evidence does not authorize test maintenance or a committed Cr2 gate.

### HP-COMP-SHELLGEOM-FN-01 - common shell decomposition

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [common\_terminal\_shell\_decomposition.md](common_terminal_shell_decomposition.md); heading `Common Terminal Shell Decomposition`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_shellification/terminal_geometry.jl`
  - `source` / `existing`: `src/pqs_source_box_route_driver_helpers.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain route-family-free direct-core/shell regions, ordering, coverage, and owned support before PQS/WL retained construction diverges.

### HP-COMP-SHELLGEOM-TEST-01 - common shell validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [common\_terminal\_shell\_decomposition.md](common_terminal_shell_decomposition.md); heading `Common Terminal Shell Decomposition`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** none. The completed evidence does not authorize test maintenance or a committed fixture.

### HP-COMP-SUPPATOM-FN-01 - supplemented one-center atoms

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `driver`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [nesting\_supplement\_composition\_plan.md](nesting_supplement_composition_plan.md); heading `Nesting/Supplement Composition`
- **Owned paths:**
  - `driver` / `existing`: `bin/cartesian_ham_builder.jl`
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** origin-centered all-electron atoms use the existing atomic supplement loader and the same PQS/WL residual-GTO/MWG path as supported diatomics.

### HP-COMP-SUPPATOM-TEST-01 - supplemented atom validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [nesting\_supplement\_composition\_plan.md](nesting_supplement_composition_plan.md); heading `Nesting/Supplement Composition`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** validation maintenance.

### HP-COMP-SUPPWL-FN-01 - supplemented WL z-axis diatomics

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [nesting\_supplement\_composition\_plan.md](nesting_supplement_composition_plan.md); heading `Nesting/Supplement Composition`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** a WL terminal basis enters the same residual-GTO, exact augmented one-body, residual MWG/IDA, assembly, and artifact path as PQS.

### HP-COMP-SUPPWL-TEST-01 - supplemented WL validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [nesting\_supplement\_composition\_plan.md](nesting_supplement_composition_plan.md); heading `Nesting/Supplement Composition`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** validation maintenance.

### HP-COMP-THINSLAB-FN-01 - common compact thin-slab lowering

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [common\_terminal\_shell\_decomposition.md](common_terminal_shell_decomposition.md); heading `Common Terminal Shell Decomposition`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`
  - `source` / `existing`: `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`
  - `source` / `existing`: `src/cartesian_retained_units/lower_contract_units.jl`
  - `source` / `existing`: `src/cartesian_shellification/terminal_geometry.jl`
  - `source` / `existing`: `src/cartesian_terminal_lowering/region_contracts.jl`
  - `source` / `existing`: `src/cartesian_terminal_lowering/selection.jl`
  - `source` / `existing`: `src/pqs_source_box_diatomic_complete_core_shell.jl`
  - `source` / `existing`: `src/pqs_source_box_route_driver_helpers.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** midpoint, outer-mismatch, and angular-z-extension slabs lower as compact face stacks for both PQS and WL, never as full identity CPBs. Maintenance includes both terminal realizers plus only conditionally required native slab metadata and route-summary caller support in the named shellification/helper files. Real shells remain family-specific.

### HP-COMP-THINSLAB-META-FN-01 - thin-slab inventory metadata

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [common\_terminal\_shell\_decomposition.md](common_terminal_shell_decomposition.md); heading `Common Terminal Shell Decomposition`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_terminal_shellification_geometry.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** describe native compact slab kinds consistently in internal inventory/scaffold summaries. It does not materialize coefficients or create artifact/report payloads.

### HP-COMP-THINSLAB-META-TEST-01 - thin-slab inventory validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [common\_terminal\_shell\_decomposition.md](common_terminal_shell_decomposition.md); heading `Common Terminal Shell Decomposition`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** none. The completed evidence does not authorize test maintenance, coefficient work, or report/artifact payloads.

### HP-COMP-THINSLAB-TEST-01 - compact thin-slab validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [common\_terminal\_shell\_decomposition.md](common_terminal_shell_decomposition.md); heading `Common Terminal Shell Decomposition`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** none. The completed evidence does not authorize test maintenance or a committed Cr2 gate.

### HP-COMP-WLDIAT-FN-01 - WL z-axis diatomic base terminal records

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [white\_lindsey\_terminal\_basis\_realization.md](white_lindsey_terminal_basis_realization.md); heading `White-Lindsey Terminal Basis Realization`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`
  - `source` / `existing`: `src/cartesian_terminal_lowering/region_contracts.jl`
  - `source` / `existing`: `src/cartesian_terminal_lowering/selection.jl`
  - `source` / `existing`: `src/cartesian_terminal_shellification_geometry.jl`
  - `source` / `existing`: `src/pqs_source_box_diatomic_complete_core_shell.jl`
  - `source` / `existing`: `src/pqs_source_box_route_driver_helpers.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain native WL z-axis diatomic terminal records and the shared base Hamiltonian path, including truthful route provenance \`:z\_axis\_diatomic\_wl\_base\`.

### HP-COMP-WLDIAT-TEST-01 - WL z-axis diatomic base validation

- **Lifecycle:** `completed`
- **Grant:** `measurement`
- **Surfaces:** `measurement`
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [white\_lindsey\_terminal\_basis\_realization.md](white_lindsey_terminal_basis_realization.md); heading `White-Lindsey Terminal Basis Realization`
- **Owned paths:**
  - `measurement` / `optional_local`: `tmp/work/wl_diatomic_base_validation.jl`
- **Evidence:**
  - `manager_pass`: `137`
- **Dependencies:** none
- **Scope:** maintain only the exact optional WL diatomic validation probe as completed evidence for \`HP-COMP-WLDIAT-FN-01\`; no committed test or execution-whitelist authority.

### HP-COMP-WLNS-FN-01 - WL diatomic ns guard

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [nesting\_supplement\_composition\_plan.md](nesting_supplement_composition_plan.md); heading `Nesting/Supplement Composition`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
- **Evidence:** none
- **Dependencies:** `HP-COMP-NS-FN-01`
- **Scope:** reject normalized WL z-axis diatomic \`ns \< 4\` before route construction and preserve retained-support saturation as valid behavior.

### HP-COMP-WLNS-TEST-01 - WL diatomic ns validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [nesting\_supplement\_composition\_plan.md](nesting_supplement_composition_plan.md); heading `Nesting/Supplement Composition`
- **Owned paths:**
  - `test` / `existing`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** validation maintenance.

### HP-CONTRACT-VEC-FN-01 - vector-backed contract plans

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [route\_stage\_metadata\_contract.md](route_stage_metadata_contract.md); heading `Route/Stage Metadata Contract`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
  - `source` / `existing`: `src/cartesian_retained_unit_transform_contracts/records.jl`
  - `source` / `existing`: `src/cartesian_retained_unit_transform_contracts/summaries.jl`
  - `source` / `existing`: `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`
  - `source` / `existing`: `src/cartesian_terminal_lowering/contracts.jl`
  - `source` / `existing`: `src/cartesian_terminal_lowering/selection.jl`
  - `source` / `existing`: `src/cartesian_terminal_lowering/summaries.jl`
  - `source` / `existing`: `src/pqs_source_box_route_driver_helpers.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain vector-backed available/selected lowering contracts and retained-unit transform contracts with unchanged accessor and order semantics. Per-contract \`source\_cpbs\` and fixed mathematical tuples remain outside this cleanup.

### HP-CONTRACT-VEC-TEST-01 - contract-plan validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [route\_stage\_metadata\_contract.md](route_stage_metadata_contract.md); heading `Route/Stage Metadata Contract`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** contract-plan validation.

### HP-DRV-ATOM-CLEAN-01 - remove hidden atom \`d\` driver residue

- **Lifecycle:** `implemented`
- **Grant:** `preservation`
- **Surfaces:** `driver`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_driver\_atom\_workflow.md](cartesian_driver_atom_workflow.md); heading `Cartesian Driver Atom Workflow`
- **Owned paths:**
  - `driver` / `existing`: `bin/cartesian_ham_builder.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** preserve the absence of hidden atom \`d\`; visible atom basis uses \`ns\`, \`core\_spacing\`, \`radius\`, and current optional fields. Do not restore the compatibility field or use this ID for new inputs, diagnostics, source algorithms, tools, fixtures, or Cr2 workflow.

### HP-DRV-ATOM-FN-01 - explicit base atom driver workflow

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `driver`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_driver\_atom\_workflow.md](cartesian_driver_atom_workflow.md); heading `Cartesian Driver Atom Workflow`
- **Owned paths:**
  - `driver` / `existing`: `bin/cartesian_ham_builder.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain \`Natom=1\`, \`basisname === nothing\` base-atom selection; explicit origin, charge, spin-sector, neutral-count, \`ns\`, \`core\_spacing\`, \`s\_factor\`, source-span/nesting, and radius-from-padding inputs; and clear unsupported-input failures. There is no public \`mode=:base\` input.

### HP-DRV-ATOM-TEST-01 - base atom driver validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [cartesian\_driver\_atom\_workflow.md](cartesian_driver_atom_workflow.md); heading `Cartesian Driver Atom Workflow`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** base atom driver validation.

### HP-DRV-ATOM-WIRE-01 - driver atom-to-base-facade wiring

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `driver`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_driver\_atom\_workflow.md](cartesian_driver_atom_workflow.md); heading `Cartesian Driver Atom Workflow`
- **Owned paths:**
  - `driver` / `existing`: `bin/cartesian_ham_builder.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** pass the explicit atom contract through the same named base producer stages, artifact writer, provenance, and readback as the base facade. Supplemented atoms remain separately governed by \`HP-COMP-SUPPATOM-\*\`.

### HP-DRV-FILE-01 - canonical driver file

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `driver`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_driver\_usability\_workflow.md](cartesian_driver_usability_workflow.md); heading `Cartesian Driver Usability Workflow`
- **Owned paths:**
  - `driver` / `existing`: `bin/cartesian_ham_builder.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain the canonical trusted local scientific driver. No other \`bin\`, tool, source, test, or committed input fixture is authorized by this ID.

### HP-DRV-FN-01 - compact functional driver workflow

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `driver`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_driver\_usability\_workflow.md](cartesian_driver_usability_workflow.md); heading `Cartesian Driver Usability Workflow`
- **Owned paths:**
  - `driver` / `existing`: `bin/cartesian_ham_builder.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain trusted input-file and \`key=value\` overrides, visible system/basis/supplement construction, \`basisname === nothing\` base selection, coarse physics timing, terminal inventory/due diligence, artifact write, and optional readback. Exact live inputs and defaults are canonical in the linked driver contract.

### HP-DRV-INV-FN-01 - canonical driver terminal-region inventory

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `driver`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_driver\_usability\_workflow.md](cartesian_driver_usability_workflow.md); heading `Cartesian Driver Usability Workflow`
- **Owned paths:**
  - `driver` / `existing`: `bin/cartesian_ham_builder.jl`
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain the bounded human-facing terminal-region inventory with region/lowering/realization kind, shell index, support/retained counts, compression, identity/product status, index/physical bounds, slab facts, and base/supplemented dimensions.

### HP-DRV-INV-TEST-01 - terminal-region inventory validation

- **Lifecycle:** `completed`
- **Grant:** `measurement`
- **Surfaces:** `measurement`
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [cartesian\_driver\_usability\_workflow.md](cartesian_driver_usability_workflow.md); heading `Cartesian Driver Usability Workflow`
- **Owned paths:**
  - `measurement` / `optional_local`: `tmp/work/terminal_inventory_native_shell_index_probe.jl`
- **Evidence:**
  - `manager_pass`: `188`
- **Dependencies:** none
- **Scope:** maintain only the exact optional terminal-inventory probe as completed evidence for \`HP-DRV-INV-FN-01\`; no committed test or driver authority.

### HP-DRV-NEST-FN-01 - construction-family driver input

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `driver`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_driver\_usability\_workflow.md](cartesian_driver_usability_workflow.md); heading `Cartesian Driver Usability Workflow`
- **Owned paths:**
  - `driver` / `existing`: `bin/cartesian_ham_builder.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain visible \`nesting = :pqs \| :wl\`, default \`:pqs\`, in driver contract/summary/readback facts. It is a construction-family choice, not a route diagnostic.

### HP-DRV-NEST-TEST-01 - construction-family validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [cartesian\_driver\_usability\_workflow.md](cartesian_driver_usability_workflow.md); heading `Cartesian Driver Usability Workflow`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** construction-family validation.

### HP-DRV-NEST-WIRE-01 - construction-family route mapping

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `driver`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_driver\_usability\_workflow.md](cartesian_driver_usability_workflow.md); heading `Cartesian Driver Usability Workflow`
- **Owned paths:**
  - `driver` / `existing`: `bin/cartesian_ham_builder.jl`
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** map \`:pqs\` to \`:pqs\_source\_box\` and \`:wl\` to \`:white\_lindsey\_low\_order\`, preserve public stage/artifact behavior, and reject unsupported combinations without exposing internal route vocabulary.

### HP-DRV-SHELLDD-FN-01 - terminal shellification due-diligence report

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `driver`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [terminal\_shellification\_due\_diligence.md](terminal_shellification_due_diligence.md); heading `Terminal Shellification Due Diligence`
- **Owned paths:**
  - `driver` / `existing`: `bin/cartesian_ham_builder.jl`
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain the bounded in-memory/report table joining terminal inventory with retained/support facts; system/geometry, axis/center/weight, dimension, and shell-row diagnostics; actual and expected source shapes; retained/final ranges; slab metadata; and advisory warning flags.

### HP-DRV-SHELLDD-TEST-01 - terminal shellification due-diligence validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `docs`
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [terminal\_shellification\_due\_diligence.md](terminal_shellification_due_diligence.md); heading `Terminal Shellification Due Diligence`
- **Owned paths:**
  - `docs` / `existing`: `docs/src/developer/designs/cartesian_hamiltonian_producer/terminal_shellification_due_diligence.md`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** docs/evidence maintenance only.

### HP-DRV-STAGE-FN-01 - visible physics-stage producer surface

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_driver\_usability\_workflow.md](cartesian_driver_usability_workflow.md); heading `Cartesian Driver Usability Workflow`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  - `source` / `existing`: `src/pqs_source_box_low_order_materialization.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain separate non-exported stages for working basis, product/moment, unit-nuclear, IDA/MWG interaction, residual augmentation, and Hamiltonian assembly so the driver can bind and time physical objects. Facades remain wrappers over the same construction.

### HP-DRV-STAGE-TEST-01 - staged driver validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [cartesian\_driver\_usability\_workflow.md](cartesian_driver_usability_workflow.md); heading `Cartesian Driver Usability Workflow`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** staged driver validation.

### HP-DRV-STAGE-WIRE-01 - canonical driver staged wiring

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `driver`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_driver\_usability\_workflow.md](cartesian_driver_usability_workflow.md); heading `Cartesian Driver Usability Workflow`
- **Owned paths:**
  - `driver` / `existing`: `bin/cartesian_ham_builder.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** call the named producer stages directly and print coarse timings for product/moment, unit-nuclear, interaction, and assembly work. Do not replace them with an opaque wrapper or expose underscored route stages, stop controls, providers, allocation probes, or solver controls.

### HP-DRV-TEST-01 - driver workflow validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [cartesian\_driver\_usability\_workflow.md](cartesian_driver_usability_workflow.md); heading `Cartesian Driver Usability Workflow`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** driver workflow validation.

### HP-FILE-01 - terminal realization file

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [terminal\_basis\_and\_base\_assembly.md](terminal_basis_and_base_assembly.md); heading `Terminal Basis And Base Assembly`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
- **Evidence:** none
- **Dependencies:** `HP-OBJ-01`, `HP-OBJ-02`
- **Scope:** maintain the implemented terminal object and PQS realization owner.

### HP-FN-00 - block-local terminal shell realization

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [terminal\_basis\_and\_base\_assembly.md](terminal_basis_and_base_assembly.md); heading `Terminal Basis And Base Assembly`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain support-local shell seed, Gram/Lowdin realization, and deterministic sign canonicalization.

### HP-FN-01 - terminal basis realizer

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [terminal\_basis\_and\_base\_assembly.md](terminal_basis_and_base_assembly.md); heading `Terminal Basis And Base Assembly`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
- **Evidence:**
  - `repo_path`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
- **Dependencies:** `HP-FN-00`
- **Scope:** maintain the live \`pqs\_terminal\_basis\_realization(...)\` signature and direct/shell/compact-slab dispatch recorded in the canonical contract.

### HP-FN-02 - structural terminal support checks

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [terminal\_basis\_and\_base\_assembly.md](terminal_basis_and_base_assembly.md); heading `Terminal Basis And Base Assembly`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain exact support equality, duplicate-row rejection, pairwise disjointness, and shell-local identity validation.

### HP-FN-03 - blockwise one-body assembly

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [terminal\_basis\_and\_base\_assembly.md](terminal_basis_and_base_assembly.md); heading `Terminal Basis And Base Assembly`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`
- **Evidence:**
  - `repo_path`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
  - `repo_path`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Dependencies:** `HP-OBJ-01`, `HP-OBJ-02`
- **Scope:** maintain block-pair product assembly and the file-local term-first Gaussian-sum accumulator within the current workspace bound.

### HP-FN-04 - localized IDA assembly

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [terminal\_basis\_and\_base\_assembly.md](terminal_basis_and_base_assembly.md); heading `Terminal Basis And Base Assembly`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_ida.jl`
- **Evidence:**
  - `repo_path`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
- **Dependencies:** none
- **Scope:** maintain blockwise final-weight-normalized localized IDA assembly and symmetry validation.

### HP-FN-05 - final Hamiltonian construction

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [terminal\_basis\_and\_base\_assembly.md](terminal_basis_and_base_assembly.md); heading `Terminal Basis And Base Assembly`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_ida_hamiltonian.jl`
- **Evidence:**
  - `repo_path`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
  - `repo_path`: `test/ida/cartesian_ida_hamiltonian_runtests.jl`
- **Dependencies:** `HP-FN-03`, `HP-FN-04`
- **Scope:** construct the existing \`CartesianIDAHamiltonian\` directly and assemble \`H1 = K + sum\_A Z\_A U\_A\` through its current accounting helpers.

### HP-HAM-MANIFEST-FN-01 - compact Hamiltonian artifact manifest

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `artifacts`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_hamiltonian\_artifact\_manifest.md](cartesian_hamiltonian_artifact_manifest.md); heading `Cartesian Hamiltonian Artifact Manifest`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  - `source` / `existing`: `src/cartesian_ida_hamiltonian.jl`
- **Evidence:**
  - `repo_path`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
  - `repo_path`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Dependencies:** none
- **Scope:** maintain matrix-order \`final\_basis\_labels/\` and \`recipe\_provenance/\` on existing facade artifacts.

### HP-HAM-MANIFEST-SRC-FN-01 - source-mode provenance seam

- **Lifecycle:** `implemented`
- **Grant:** `implementation`
- **Surfaces:** `artifacts`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_hamiltonian\_artifact\_manifest.md](cartesian_hamiltonian_artifact_manifest.md); heading `Cartesian Hamiltonian Artifact Manifest`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
  - `source` / `existing`: `src/cartesian_raw_product_sources/CartesianRawProductSources.jl`
  - `source` / `existing`: `src/cartesian_raw_product_sources/records.jl`
  - `source` / `existing`: `src/cartesian_raw_product_sources/source_mode_indices.jl`
  - `source` / `existing`: `src/cartesian_retained_unit_transform_contracts/CartesianRetainedUnitTransformContracts.jl`
  - `source` / `existing`: `src/cartesian_retained_unit_transform_contracts/records.jl`
  - `source` / `existing`: `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`
  - `source` / `existing`: `src/cartesian_retained_units/CartesianRetainedUnits.jl`
  - `source` / `existing`: `src/cartesian_retained_units/lower_contract_units.jl`
  - `source` / `existing`: `src/cartesian_retained_units/records.jl`
  - `source` / `existing`: `src/cartesian_terminal_lowering/contracts.jl`
  - `source` / `existing`: `src/cartesian_terminal_lowering/region_contracts.jl`
- **Evidence:**
  - `repo_path`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
  - `repo_path`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Dependencies:** `HP-HAM-MANIFEST-FN-01`
- **Scope:** maintain the existing \`source\_mode\_provenance\` carrier and add only missing construction-native relations or labels on the listed surfaces.

### HP-HAM-MANIFEST-SRC-TEST-01 - source-mode provenance seam validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_hamiltonian\_artifact\_manifest.md](cartesian_hamiltonian_artifact_manifest.md); heading `Cartesian Hamiltonian Artifact Manifest`
- **Owned paths:**
  - `test` / `existing`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:** none
- **Dependencies:** `HP-HAM-MANIFEST-SRC-FN-01`
- **Scope:** maintain native source-group presence, status, no-inference, and readback checks for the implemented subset.

### HP-HAM-MANIFEST-TEST-01 - artifact manifest validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_hamiltonian\_artifact\_manifest.md](cartesian_hamiltonian_artifact_manifest.md); heading `Cartesian Hamiltonian Artifact Manifest`
- **Owned paths:**
  - `test` / `existing`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:** none
- **Dependencies:** `HP-HAM-MANIFEST-FN-01`
- **Scope:** maintain base/supplemented readback and selected provenance-root checks. The accepted direct-JLD2 evidence owns detailed field/status coverage.

### HP-MCOMX-DRV-FN-01 - canonical source-span selector

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `driver`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [mapped\_comx\_source\_span.md](mapped_comx_source_span.md); heading `Mapped-COMX Source Span`
- **Owned paths:**
  - `driver` / `existing`: `bin/cartesian_ham_builder.jl`
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/pqs_source_box_route_driver_helpers.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** expose \`source\_span = :ordinary \| :mapped\_comx\`, default ordinary, and reject mapped-COMX with White-Lindsey. This is not a diagnostic route switch.

### HP-MCOMX-DRV-TEST-01 - driver source-span validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [mapped\_comx\_source\_span.md](mapped_comx_source_span.md); heading `Mapped-COMX Source Span`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** none. The completed evidence does not authorize test maintenance or a new committed test surface.

### HP-MCOMX-FILE-01 - mapped-COMX source ownership

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [mapped\_comx\_source\_span.md](mapped_comx_source_span.md); heading `Mapped-COMX Source Span`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_nested_faces.jl`
  - `source` / `existing`: `src/cartesian_pair_block_materialization/pqs_source_axis_transforms.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** source maintenance in the existing owners only; no new production file or second COMX owner.

### HP-MCOMX-FN-01 - mapped source-span construction

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [mapped\_comx\_source\_span.md](mapped_comx_source_span.md); heading `Mapped-COMX Source Span`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_nested_faces.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** construct normalized-local-coordinate mapped enrichment before the existing physical-coordinate COMX cleanup. Ordinary behavior remains default.

### HP-MCOMX-OBJ-01 - mapped source specification

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [mapped\_comx\_source\_span.md](mapped_comx_source_span.md); heading `Mapped-COMX Source Span`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_nested_faces.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** fixed protected-P2, mapped-Chebyshev, lambda/no-sqrt-J, and physical-localization facts. No public export or general tuning object.

### HP-MCOMX-TERM-FN-01 - terminal shell-seed consumption

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [mapped\_comx\_source\_span.md](mapped_comx_source_span.md); heading `Mapped-COMX Source Span`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** validate and consume materialized carried axis facts as the PQS shell seed while preserving ordinary fallback, boundary selection, support restriction, Lowdin, and canonicalization.

### HP-MCOMX-TERM-TEST-01 - terminal seam validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [mapped\_comx\_source\_span.md](mapped_comx_source_span.md); heading `Mapped-COMX Source Span`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** none. The completed evidence does not authorize test maintenance or a new committed test surface.

### HP-MCOMX-TEST-01 - mapped source validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [mapped\_comx\_source\_span.md](mapped_comx_source_span.md); heading `Mapped-COMX Source Span`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** none. The completed evidence does not authorize test maintenance, a committed fixture, or default-promotion work.

### HP-MCOMX-WIRE-01 - PQS axis-transform wiring

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [mapped\_comx\_source\_span.md](mapped_comx_source_span.md); heading `Mapped-COMX Source Span`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_pair_block_materialization/pqs_source_axis_transforms.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** pass the internal source-span choice into the existing doside seam and return ordinary carried \`AxisSourceTransformFact\` objects.

### HP-NEST-ART-FN-01 - nesting artifact-truth cleanup

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `artifacts`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_hamiltonian\_artifact\_manifest.md](cartesian_hamiltonian_artifact_manifest.md); heading `Cartesian Hamiltonian Artifact Manifest`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`
- **Evidence:**
  - `repo_path`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
  - `repo_path`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Dependencies:** none
- **Scope:** maintain truthful \`nesting\` and route labels in producer/recipe provenance and the nesting-neutral final-basis module description.

### HP-NEST-ART-TEST-01 - nesting artifact-truth validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_hamiltonian\_artifact\_manifest.md](cartesian_hamiltonian_artifact_manifest.md); heading `Cartesian Hamiltonian Artifact Manifest`
- **Owned paths:**
  - `test` / `existing`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:** none
- **Dependencies:** `HP-COMP-SUPPWL-TEST-01`, `HP-NEST-ART-FN-01`
- **Scope:** maintain small artifact/readback provenance checks for truthful nesting and route labels.

### HP-OBJ-01 - \`CartesianTerminalBasisBlock\`

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [terminal\_basis\_and\_base\_assembly.md](terminal_basis_and_base_assembly.md); heading `Terminal Basis And Base Assembly`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
- **Evidence:**
  - `repo_path`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
  - `repo_path`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Dependencies:** none
- **Scope:** preserve the exact support-local block object and direct-identity versus compact-coefficient semantics.

### HP-OBJ-02 - \`CartesianTerminalBasisRealization\`

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [terminal\_basis\_and\_base\_assembly.md](terminal_basis_and_base_assembly.md); heading `Terminal Basis And Base Assembly`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
- **Evidence:**
  - `repo_path`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
- **Dependencies:** `HP-OBJ-01`
- **Scope:** preserve the exact realization fields, native block order, final dimension, and structural-overlap interpretation. \`max\_cross\_overlap\` remains legacy object shape, not a physical repair signal.

### HP-OBJ-03 - generic build-result wrapper — rejected

- **Lifecycle:** `rejected`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [terminal\_basis\_and\_base\_assembly.md](terminal_basis_and_base_assembly.md); heading `Terminal Basis And Base Assembly`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** generic build-result wrapper — rejected.

### HP-PQS-ASPECTSHELL-FN-01 - PQS aspect-aware source modes

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_complete\_shell\_aspect\_source\_modes.md](pqs_complete_shell_aspect_source_modes.md); heading `PQS Complete-Shell Aspect-Aware Source Modes`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/pqs_multilayer_shell_region_plan.jl`
  - `source` / `existing`: `src/pqs_multilayer_shell_source_plan.jl`
  - `source` / `existing`: `src/pqs_source_box_route_driver_helpers.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** preserve the existing post-shellification angular-band \`L\` selection and one authoritative \`(q,q,L)\` shape through lowering, retention, realization, and reporting.

### HP-PQS-ASPECTSHELL-TEST-01 - aspect-source validation

- **Lifecycle:** `completed`
- **Grant:** `measurement`
- **Surfaces:** `measurement`
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [pqs\_complete\_shell\_aspect\_source\_modes.md](pqs_complete_shell_aspect_source_modes.md); heading `PQS Complete-Shell Aspect-Aware Source Modes`
- **Owned paths:**
  - `measurement` / `optional_local`: `tmp/work/h2plus_aspect_shell_completeness_replay.jl`
  - `measurement` / `optional_local`: `tmp/work/pqs_aspect_shell_validation.jl`
- **Evidence:**
  - `manager_pass`: `247`
  - `manager_pass`: `248`
- **Dependencies:** none
- **Scope:** maintain only the two exact optional aspect-source probes as completed evidence for \`HP-PQS-ASPECTSHELL-FN-01\`; no committed test authority.

### HP-PQS-ATOMREF-PACKET-FN-01 - reusable atomic HF reference packets

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `artifacts`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [atomic\_hf\_reference\_packets.md](atomic_hf_reference_packets.md); heading `Atomic HF Reference Packets`
- **Owned paths:**
  - `source` / `existing`: `src/GaussletBases.jl`
  - `source` / `existing`: `src/cartesian_reference_density/CartesianReferenceDensity.jl`
  - `source` / `existing`: `src/cartesian_reference_density/atomic_hf_reference_packets.jl`
  - `source` / `existing`: `src/cartesian_reference_density/screened_hartree_correction.jl`
- **Evidence:** none
- **Dependencies:** `HP-RHO0-MIXH-GAAA-FN-01`, `HP-RHO0-MIXH-GG-FN-01`
- **Scope:** maintain converged one-center determinant packets, exact packet self-integrity, exact owner/order/placement mapping, numerical owner-local overlap equivalence at \`1e-10\`, ordinary density and radial-potential fits, read/write validation, and explicit fit/provenance diagnostics. Density fits own \`E0\`; potential fits approximate \`J0\`. Polished legacy packets reject.

### HP-PQS-ATOMREF-PACKET-TEST-01 - atomic packet validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [atomic\_hf\_reference\_packets.md](atomic_hf_reference_packets.md); heading `Atomic HF Reference Packets`
- **Owned paths:**
  - `test` / `existing`: `test/misc/runtests.jl`
  - `test` / `existing`: `test/nested/cartesian_atomic_hf_reference_packet_runtests.jl`
  - `test` / `existing`: `test/nested/cartesian_screened_hartree_correction_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain packet roundtrip, convergence, determinant/fit, fingerprint, embedding-equivalence, compact-Coulomb-role, and malformed-input coverage. The scientific body of \`data/legacy/BasisSets\` is not test authority to rewrite that data.

### HP-PQS-ATOMREF-POTMOM-FN-01 - retired determinant-moment polish

- **Lifecycle:** `retired`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [atomic\_hf\_reference\_packets.md](atomic_hf_reference_packets.md); heading `Atomic HF Reference Packets`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** retired determinant-moment polish.

### HP-PQS-ATOMREF-POTMOM-TEST-01 - retired polish validation

- **Lifecycle:** `retired`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [atomic\_hf\_reference\_packets.md](atomic_hf_reference_packets.md); heading `Atomic HF Reference Packets`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** retired polish validation.

### HP-PQS-COULOMB-ACCURACY-FN-01 - producer-wide Coulomb accuracy policy

- **Lifecycle:** `implemented`
- **Grant:** `implementation`
- **Surfaces:** `artifacts`, `driver`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [coulomb\_accuracy\_policy.md](coulomb_accuracy_policy.md); heading `Producer-Wide Coulomb Accuracy Policy`
- **Owned paths:**
  - `driver` / `existing`: `bin/cartesian_ham_builder.jl`
  - `source` / `existing`: `src/GaussianAnalyticIntegrals.jl`
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  - `source` / `existing`: `src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl`
  - `source` / `existing`: `src/cartesian_ida_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_protected_ladder_bundle.jl`
  - `source` / `existing`: `src/cartesian_reference_density/atomic_hf_reference_packets.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/mwg_interaction.jl`
  - `source` / `existing`: `src/ordinary_coulomb.jl`
  - `source` / `existing`: `src/pqs_source_box_low_order_materialization.jl`
  - `source` / `existing`: `src/pqs_source_box_route_driver_helpers.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** source maintenance, including only completion of the already approved Standard60 fingerprint/provenance and canonical-driver exposure.

### HP-PQS-COULOMB-ACCURACY-TEST-01 - Coulomb accuracy validation

- **Lifecycle:** `approved`
- **Grant:** `implementation`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [coulomb\_accuracy\_policy.md](coulomb_accuracy_policy.md); heading `Producer-Wide Coulomb Accuracy Policy`
- **Owned paths:**
  - `test` / `existing`: `test/core/runtests.jl`
  - `test` / `existing`: `test/docs/cartesian_ham_builder_policy_runtests.jl`
  - `test` / `existing`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
  - `test` / `existing`: `test/nested/cartesian_atomic_hf_reference_packet_runtests.jl`
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** validation maintenance for implemented compact/high behavior and the already-approved Standard60/driver completion only.

### HP-PQS-MAP-SFACTOR-FN-01 - expert mapping \`s\_factor\` keyword

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `artifacts`, `driver`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_mapping\_s\_factor.md](pqs_mapping_s_factor.md); heading `` PQS/WL Mapping `s_factor` ``
- **Owned paths:**
  - `driver` / `existing`: `bin/cartesian_ham_builder.jl`
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_protected_ladder_bundle.jl`
  - `source` / `existing`: `src/mappings.jl`
  - `source` / `existing`: `src/pqs_source_box_route_driver_helpers.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain finite positive \`s\_factor\`, default \`1.0\`, one-center \`effective\_s = s\_factor\*sqrt(Z\*core\_spacing)\`, the analogous per-center multicenter combined-inverse-sqrt input, and explicit standard/effective provenance.

### HP-PQS-MAP-SFACTOR-TEST-01 - mapping \`s\_factor\` validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_mapping\_s\_factor.md](pqs_mapping_s_factor.md); heading `` PQS/WL Mapping `s_factor` ``
- **Owned paths:**
  - `test` / `existing`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain only the existing committed \`s\_factor\` assertions in that file. New fixtures or endpoint policy require separate authority.

### HP-PQS-SCREEN-HARTREE-AUDIT-01 - protected-GTO screened Hartree residual-density audit

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [screened\_hartree\_residual\_density.md](screened_hartree_residual_density.md); heading `Screened Hartree Residual-Density Formalism`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** protected-GTO screened Hartree residual-density audit.

### HP-PQS-SCREEN-HARTREE-CORR-FN-01 - internal screened-Hartree correction

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [screened\_hartree\_correction\_assembly.md](screened_hartree_correction_assembly.md); heading `Screened Hartree Correction Assembly`
- **Owned paths:**
  - `source` / `existing`: `src/GaussletBases.jl`
  - `source` / `existing`: `src/cartesian_reference_density/CartesianReferenceDensity.jl`
  - `source` / `existing`: `src/cartesian_reference_density/atomic_hf_reference_packets.jl`
  - `source` / `existing`: `src/cartesian_reference_density/screened_hartree_correction.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** consume represented converged references and same-basis \`V\_IDA\`, \`J0\_G\`, and \`E0\_G\`; return in-memory \`Delta\_J0 = J0\_G - Diagonal(V\_IDA\*q0)\` and \`C = 0.5\*q0'V\_IDA\*q0 - 0.5\*E0\_G\`; preserve strict representation, finiteness, symmetry, convergence, and derivative/algebra failures while reporting ordinary fitted-potential energy inconsistency.

### HP-PQS-SCREEN-HARTREE-CORR-TEST-01 - correction validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [screened\_hartree\_correction\_assembly.md](screened_hartree_correction_assembly.md); heading `Screened Hartree Correction Assembly`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_screened_hartree_correction_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain packet/reference consistency, same-basis, anchor, derivative, symmetry/finiteness, fitted-potential reporting, and malformed input coverage without adding physics endpoint assertions.

### HP-PQS-SCREEN-HARTREE-NE-AUDIT-01 - Ne screened Hartree endpoint measurement

- **Lifecycle:** `superseded`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [screened\_hartree\_residual\_density.md](screened_hartree_residual_density.md); heading `Screened Hartree Residual-Density Formalism`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** Ne screened Hartree endpoint measurement.

### HP-PQS-SCREEN-HARTREE-NE-FITCLOUD-AUDIT-01 - Ne fitted-cloud measurement

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [screened\_hartree\_residual\_density.md](screened_hartree_residual_density.md); heading `Screened Hartree Residual-Density Formalism`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** Ne fitted-cloud measurement.

### HP-PQS-SCREEN-HARTREE-POTFIT-AUDIT-01 - fitted-potential measurement

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [atomic\_hf\_reference\_packets.md](atomic_hf_reference_packets.md); heading `Atomic HF Reference Packets`
  - `evidence` [screened\_hartree\_residual\_density.md](screened_hartree_residual_density.md); heading `Screened Hartree Residual-Density Formalism`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** fitted-potential measurement.

### HP-PQS-SHELLQ-OVERRIDE-FN-01 - semantic per-shell PQS source-q overrides

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_semantic\_shell\_q\_overrides.md](pqs_semantic_shell_q_overrides.md); heading `Semantic Per-Shell PQS Source-q Overrides`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_protected_ladder_bundle.jl`
  - `source` / `existing`: `src/pqs_source_box_route_driver_helpers.jl`
- **Evidence:** none
- **Dependencies:** `HP-PQS-ASPECTSHELL-FN-01`, `HP-RG-NUMCOMP-FN-01`
- **Scope:** Maintain \`owner = :all\` overrides for positive semantic \`:atom\_local\_shell\` or \`:shared\_molecular\_shell\` indices with non-Boolean integer \`source\_q \>= 3\` and \`source\_q \!= route\_q\`. Values above route q refine the shell contraction; values below route q coarsen it while parent axes, support, ownership, cores, slabs, and route metadata remain unchanged. Atom-local shells use \`(source\_q,source\_q,source\_q)\`; shared shells rerun the existing angular-band selector for \`(source\_q,source\_q,L)\`. Selector retention, \`nside\`, and \`selected\_q\` use \`source\_q\`. One authoritative shape must reach existing retained/support/transform/realization/due-diligence consumers unchanged.

### HP-PQS-SHELLQ-OVERRIDE-TEST-01 - semantic source-q override validation

- **Lifecycle:** `approved`
- **Grant:** `implementation`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_semantic\_shell\_q\_overrides.md](pqs_semantic_shell_q_overrides.md); heading `Semantic Per-Shell PQS Source-q Overrides`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** Validate route-q 7 to source-q 6 and 5 coarsening, expected retained-count reduction, unchanged parent/support/ownership/cores/slabs/route metadata, orthonormal contraction columns, omitted/empty parity, finite symmetric full construction, and rejection of malformed, below-3, Boolean, equal-route, asymmetric, and unmatched requests. Preserve existing refinement, residual, packet-capture, \`J0/E0\`, correction, dimension, and due-diligence gates. No new accessor, dense baseline-to-variant overlap API, source-pass HF, endpoint energy assertion, or production claim is approved.

### HP-R1-ART-01 - public base producer artifact provenance

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `artifacts`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r1\_public\_base\_producer.md](r1_public_base_producer.md); heading `R1 Public Base Producer`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_ida_hamiltonian.jl`
- **Evidence:**
  - `repo_path`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
- **Dependencies:** none
- **Scope:** maintain the fixed \`producer\_provenance/\` keys and truthful route/size/mapping/system values listed in the canonical R1 contract.

### HP-R1-ATOM-FN-01 - explicit one-center all-electron base atom facade

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r1\_one\_center\_base\_atoms.md](r1_one_center_base_atoms.md); heading `R1 One-Center Base Atoms`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain explicit origin-centered, neutral, all-electron atom validation in the existing \`cartesian\_base\_hamiltonian(system; basis, hamfile)\` facade. Charge, electron counts, spin sectors, basis, and ECP behavior must never be inferred from the atom label.

### HP-R1-ATOM-TEST-01 - one-center base atom validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r1\_one\_center\_base\_atoms.md](r1_one_center_base_atoms.md); heading `R1 One-Center Base Atoms`
- **Owned paths:**
  - `test` / `existing`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain H regression, malformed atom/basis input rejection, finite/symmetric operator, mapping/provenance, and artifact/readback checks.

### HP-R1-ATOM-WIRE-01 - one-center atom shared workflow wiring

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r1\_one\_center\_base\_atoms.md](r1_one_center_base_atoms.md); heading `R1 One-Center Base Atoms`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** map explicit charge and \`core\_spacing\` into the private atomic mapping, derive physical parent extent from \`radius\`, and use the same terminal, one-body, IDA, Hamiltonian, writer, and provenance machinery as the supported base producer. Atom routes remain \`:one\_center\_pqs\_base\` or \`:one\_center\_wl\_base\` according to nesting.

### HP-R1-CORE-FN-01 - unified core-spacing producer contract

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r1\_public\_base\_producer.md](r1_public_base_producer.md); heading `R1 Public Base Producer`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
- **Evidence:**
  - `repo_path`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
- **Dependencies:** none
- **Scope:** maintain \`core\_spacing\` as the single public near-nucleus physical scale and atom-only compatibility \`d == core\_spacing\`.

### HP-R1-FILE-01 - public base producer source file

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r1\_public\_base\_producer.md](r1_public_base_producer.md); heading `R1 Public Base Producer`
- **Owned paths:**
  - `source` / `existing`: `src/GaussletBases.jl`
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
- **Evidence:**
  - `repo_path`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
- **Dependencies:** none
- **Scope:** maintain the public base producer in its current owner file and the existing \`cartesian\_base\_hamiltonian\` export.

### HP-R1-FN-01 - public base Hamiltonian producer facade

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r1\_public\_base\_producer.md](r1_public_base_producer.md); heading `R1 Public Base Producer`
- **Owned paths:**
  - `source` / `existing`: `src/GaussletBases.jl`
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
- **Evidence:**
  - `repo_path`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
- **Dependencies:** `HP-R1-CORE-FN-01`
- **Scope:** maintain \`cartesian\_base\_hamiltonian(system; basis, hamfile=nothing)\` with plain \`NamedTuple\` inputs and direct \`CartesianIDAHamiltonian{Float64}\` return.

### HP-R1-TEST-01 - public base producer endpoint test/example

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r1\_public\_base\_producer.md](r1_public_base_producer.md); heading `R1 Public Base Producer`
- **Owned paths:**
  - `test` / `existing`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain the standalone atom/H2 endpoint, malformed-input, deprecated-\`d\`, geometry, matrix, Coulomb, artifact, and provenance checks.

### HP-R1-WIRE-01 - report-free base producer wiring

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r1\_public\_base\_producer.md](r1_public_base_producer.md); heading `R1 Public Base Producer`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
- **Evidence:**
  - `repo_path`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
- **Dependencies:** none
- **Scope:** maintain the report-free staged path used by the public facade and human-facing driver.

### HP-R3-ART-01 - compact supplemented artifact provenance

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `artifacts`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_hamiltonian\_artifact\_manifest.md](cartesian_hamiltonian_artifact_manifest.md); heading `Cartesian Hamiltonian Artifact Manifest`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
- **Evidence:**
  - `repo_path`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Dependencies:** `HP-HAM-MANIFEST-FN-01`, `HP-HAM-MANIFEST-SRC-FN-01`
- **Scope:** maintain the compact \`supplement\_provenance/\` group and existing writer composition recorded in the artifact contract.

### HP-R3-FN-01 - residual-basis construction

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [residual\_gaussian\_domain\_module.md](residual_gaussian_domain_module.md); heading `Residual Gaussian Domain Module`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/residual_basis.jl`
- **Evidence:**
  - `repo_path`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Dependencies:** `HP-RG-FN-01`
- **Scope:** maintain delegation to owner-local residual selection and one final inter-owner merge.

### HP-R3-FN-02 - exact augmented one-body and moment assembly

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [residual\_gaussian\_domain\_module.md](residual_gaussian_domain_module.md); heading `Residual Gaussian Domain Module`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/augmented_operators.jl`
- **Evidence:**
  - `repo_path`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Dependencies:** `HP-RG-FN-02`
- **Scope:** maintain exact augmented kinetic, by-center unit-nuclear, coordinate, and second-moment transformations.

### HP-R3-FN-03 - residual MWG/IDA and in-memory Hamiltonian

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r3\_residual\_gto\_mwg\_augmentation.md](r3_residual_gto_mwg_augmentation.md); heading `R3 Residual-GTO/MWG Compatibility History`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/mwg_interaction.jl`
- **Evidence:**
  - `repo_path`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Dependencies:** `HP-RG-FN-03`, `HP-RG-FN-04`
- **Scope:** maintain \`pqs\_terminal\_residual\_gto\_augmented\_hamiltonian(...)\` and direct return of the existing \`CartesianIDAHamiltonian{Float64}\`.

### HP-R3-OBJ-01 - residual-GTO augmentation object

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r3\_residual\_gto\_mwg\_augmentation.md](r3_residual_gto_mwg_augmentation.md); heading `R3 Residual-GTO/MWG Compatibility History`
  - `canonical` [residual\_gaussian\_domain\_module.md](residual_gaussian_domain_module.md); heading `Residual Gaussian Domain Module`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/residual_basis.jl`
- **Evidence:**
  - `repo_path`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Dependencies:** `HP-RG-OBJ-01`
- **Scope:** maintain \`CartesianTerminalResidualGTOAugmentation\` only as the live alias of \`CartesianResidualGaussianBasis\` required by compatibility callers.

### HP-R3-TEST-01 - residual-GTO/MWG compatibility endpoint

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r3\_residual\_gto\_mwg\_augmentation.md](r3_residual_gto_mwg_augmentation.md); heading `R3 Residual-GTO/MWG Compatibility History`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain the standalone H2 residual geometry, exact-operator, independent interaction, in-memory Hamiltonian, and artifact compatibility gate.

### HP-R3BASE-DRV-TEST-01 - canonical-driver reuse validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `docs`
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [r3\_same\_construction\_base\_reuse.md](r3_same_construction_base_reuse.md); heading `R3 Same-Construction Base Reuse`
- **Owned paths:**
  - `docs` / `existing`: `docs/src/developer/designs/cartesian_hamiltonian_producer/r3_same_construction_base_reuse.md`
- **Evidence:**
  - `repo_path`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Dependencies:** `HP-R3BASE-DRV-WIRE-01`
- **Scope:** maintain the accepted call-site validation record.

### HP-R3BASE-DRV-WIRE-01 - canonical-driver K/U reuse wiring

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `driver`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r3\_same\_construction\_base\_reuse.md](r3_same_construction_base_reuse.md); heading `R3 Same-Construction Base Reuse`
- **Owned paths:**
  - `driver` / `existing`: `bin/cartesian_ham_builder.jl`
- **Evidence:** none
- **Dependencies:** `HP-R3BASE-FN-01`
- **Scope:** maintain passing the already-built base kinetic and by-center unit nuclear matrices to supplemented exact-operator construction.

### HP-R3BASE-FN-01 - same-construction base K/U reuse

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r3\_same\_construction\_base\_reuse.md](r3_same_construction_base_reuse.md); heading `R3 Same-Construction Base Reuse`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
- **Evidence:** none
- **Dependencies:** `HP-R3GG-FN-01`, `HP-R3UN-FN-01`
- **Scope:** maintain trusted \`base\_kinetic\` and \`base\_unit\_nuclear\` handoff, dimension/center checks, and live exact recomputation fallbacks.

### HP-R3BASE-TEST-01 - same-construction reuse validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r3\_same\_construction\_base\_reuse.md](r3_same_construction_base_reuse.md); heading `R3 Same-Construction Base Reuse`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain focused fallback/reuse parity, exact-operator, endpoint, and artifact-readback coverage.

### HP-R3GG-FN-01 - terminal G-G product matrices

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r3\_terminal\_gg\_product\_matrices.md](r3_terminal_gg_product_matrices.md); heading `R3 Terminal G-G Product Matrices`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain exact kinetic and first/second moment \`G-G\` assembly with function-local scratch reuse.

### HP-R3GG-TEST-01 - terminal G-G validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r3\_terminal\_gg\_product\_matrices.md](r3_terminal_gg_product_matrices.md); heading `R3 Terminal G-G Product Matrices`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:** none
- **Dependencies:** `HP-R3GG-FN-01`
- **Scope:** maintain the focused endpoint, parity, finiteness, and symmetry coverage. No separate product-framework fixture is authorized.

### HP-R3REM-AUDIT-01 - remaining exact-operator allocation audit

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [r3\_remaining\_exact\_operator\_allocation\_audit.md](r3_remaining_exact_operator_allocation_audit.md); heading `R3 Remaining Exact-Operator Allocation Audit`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** remaining exact-operator allocation audit.

### HP-R3U-FILE-01 - supplemented workflow source and validation files

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`, `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r3\_usability\_supplemented\_workflow.md](r3_usability_supplemented_workflow.md); heading `R3 Usability Supplemented Workflow`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain the non-exported facade in the existing source and test files.

### HP-R3U-FN-01 - non-exported supplemented Hamiltonian facade

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r3\_usability\_supplemented\_workflow.md](r3_usability_supplemented_workflow.md); heading `R3 Usability Supplemented Workflow`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
- **Evidence:**
  - `repo_path`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Dependencies:** none
- **Scope:** maintain \`cartesian\_residual\_gto\_mwg\_hamiltonian(system; basis, supplement, hamfile)\` as a non-exported direct-\`CartesianIDAHamiltonian{Float64}\` facade.

### HP-R3U-TEST-01 - supplemented facade endpoint

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r3\_usability\_supplemented\_workflow.md](r3_usability_supplemented_workflow.md); heading `R3 Usability Supplemented Workflow`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:** none
- **Dependencies:** `HP-R3-TEST-01`
- **Scope:** maintain H2 facade type, endpoint, malformed-input, artifact/readback, and provenance checks.

### HP-R3U-WIRE-01 - base-to-RG same-construction workflow

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r3\_usability\_supplemented\_workflow.md](r3_usability_supplemented_workflow.md); heading `R3 Usability Supplemented Workflow`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
- **Evidence:**
  - `repo_path`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Dependencies:** none
- **Scope:** maintain one same-construction path from validated input through the direct in-memory Hamiltonian and optional existing artifact.

### HP-R3U-ZDI-FN-01 - homonuclear z-axis diatomic supplemented facade

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r3\_homonuclear\_diatomic\_supplemented\_workflow.md](r3_homonuclear_diatomic_supplemented_workflow.md); heading `R3 Homonuclear Z-Axis Diatomic Supplemented Workflow`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
- **Evidence:** none
- **Dependencies:** `HP-R3U-FN-01`
- **Scope:** maintain explicit neutral all-electron homonuclear two-center z-axis validation and optional trusted supplement \`basisfile\`.

### HP-R3U-ZDI-TEST-01 - homonuclear diatomic validation

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r3\_homonuclear\_diatomic\_supplemented\_workflow.md](r3_homonuclear_diatomic_supplemented_workflow.md); heading `R3 Homonuclear Z-Axis Diatomic Supplemented Workflow`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:** none
- **Dependencies:** `HP-R3U-TEST-01`
- **Scope:** maintain only the existing homonuclear H2 assertions in that committed file. New fixtures or broader geometry require separate authority.

### HP-R3U-ZDI-WIRE-01 - canonical driver supplemented-mode wiring

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `driver`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r3\_homonuclear\_diatomic\_supplemented\_workflow.md](r3_homonuclear_diatomic_supplemented_workflow.md); heading `R3 Homonuclear Z-Axis Diatomic Supplemented Workflow`
- **Owned paths:**
  - `driver` / `existing`: `bin/cartesian_ham_builder.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain generic supplemented-mode wiring for explicit homonuclear z-axis inputs through the same producer construction.

### HP-R3UN-FN-01 - terminal unit-nuclear U\_GG Gaussian sum

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r3\_unit\_nuclear\_ugg\_gaussian\_sum.md](r3_unit_nuclear_ugg_gaussian_sum.md); heading `R3 Unit-Nuclear U_GG Gaussian Sum`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain exact uncharged by-center \`U\_GG\`, term-first assembly, and function-local scratch reuse.

### HP-R3UN-TEST-01 - terminal unit-nuclear validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r3\_unit\_nuclear\_ugg\_gaussian\_sum.md](r3_unit_nuclear_ugg_gaussian_sum.md); heading `R3 Unit-Nuclear U_GG Gaussian Sum`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:** none
- **Dependencies:** `HP-R3UN-FN-01`
- **Scope:** maintain the focused endpoint, fallback, parity, finiteness, and symmetry coverage. No separate Gaussian-sum fixture is authorized.

### HP-RAW-SRCMODE-FN-01 - raw product source-mode inventory

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [route\_stage\_metadata\_contract.md](route_stage_metadata_contract.md); heading `Route/Stage Metadata Contract`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
  - `source` / `existing`: `src/cartesian_raw_product_sources/records.jl`
  - `source` / `existing`: `src/cartesian_raw_product_sources/source_mode_indices.jl`
  - `source` / `existing`: `src/cartesian_raw_product_sources/summaries.jl`
  - `source` / `existing`: `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain vector-backed mode/column inventories while preserving fixed \`NTuple{3,Int}\` coordinates, deterministic mode order, retained-rule association, and manifest source provenance.

### HP-RAW-SRCMODE-TEST-01 - raw source-mode validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [route\_stage\_metadata\_contract.md](route_stage_metadata_contract.md); heading `Route/Stage Metadata Contract`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** raw source-mode validation.

### HP-REP-XGTO-IMPORT-FN-01 - external GTO orbital import

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [external\_gto\_orbital\_import.md](external_gto_orbital_import.md); heading `External GTO Orbital Import`
- **Owned paths:**
  - `source` / `existing`: `src/GaussletBases.jl`
  - `source` / `existing`: `src/cartesian_external_gto_import.jl`
- **Evidence:**
  - `repo_path`: `test/nested/cartesian_external_gto_import_runtests.jl`
- **Dependencies:** none
- **Scope:** maintain explicit packet validation, \`S\_FG\*C\_G\` import, spin-resolved capture, and direct unorthonormalized protected \`S\_LG\*C\_G\` composition.

### HP-REP-XGTO-IMPORT-TEST-01 - external GTO orbital import validation

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [external\_gto\_orbital\_import.md](external_gto_orbital_import.md); heading `External GTO Orbital Import`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_external_gto_import_runtests.jl`
- **Evidence:** none
- **Dependencies:** `HP-REP-XGTO-IMPORT-FN-01`
- **Scope:** maintain packet identity/order/\`S\_GG\`, restricted/spin-resolved import, capture, rotation-invariance, and malformed-input checks.

### HP-REP-XGTO-PROTECT-SIDECAR-FN-01 - protected external-GTO representation sidecar

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `artifacts`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [external\_gto\_orbital\_import.md](external_gto_orbital_import.md); heading `External GTO Orbital Import`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_external_gto_import.jl`
- **Evidence:**
  - `repo_path`: `test/nested/cartesian_external_gto_import_runtests.jl`
- **Dependencies:** none
- **Scope:** maintain the standalone native-order, final-by-external v1 sidecar, exact \`S\_LG\`, direct spin imports, packet/member fingerprints, and metric-aware capture diagnostics.

### HP-REP-XGTO-PROTECT-SIDECAR-TEST-01 - protected external-GTO sidecar validation

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [external\_gto\_orbital\_import.md](external_gto_orbital_import.md); heading `External GTO Orbital Import`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_external_gto_import_runtests.jl`
- **Evidence:** none
- **Dependencies:** `HP-REP-XGTO-PROTECT-SIDECAR-FN-01`
- **Scope:** maintain exact key/identity, roundtrip, saved-overlap reimport, rectangular capture, packet/member/artifact, and tamper-rejection checks.

### HP-RES-01 - terminal basis build result — rejected

- **Lifecycle:** `rejected`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [terminal\_basis\_and\_base\_assembly.md](terminal_basis_and_base_assembly.md); heading `Terminal Basis And Base Assembly`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** terminal basis build result — rejected.

### HP-RETIRE-CARRIED-SPACE-FN-01 - retire orphaned Cartesian carried-space adapter

- **Lifecycle:** `retired`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [qw\_high\_order\_experimental\_retirement.md](qw_high_order_experimental_retirement.md); heading `Cartesian Carried-Space Adapter Retirement`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** Delete the orphaned internal CartesianCarriedSpaces adapter and its sole GaussletBases include without aliases, stubs, deprecations, or replacement machinery.

### HP-RETIRE-CARRIED-SPACE-TEST-01 - validate Cartesian carried-space adapter retirement

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [qw\_high\_order\_experimental\_retirement.md](qw_high_order_experimental_retirement.md); heading `Cartesian Carried-Space Adapter Retirement`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** Use the unchanged core geometry gate plus package-load and deleted-symbol scans to validate adapter retirement; add no replacement tests.

### HP-RETIRE-CCS-RHF-FN-01 - remove stale RHF payload stack

- **Lifecycle:** `retired`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [complete\_core\_shell\_rhf\_retirement.md](complete_core_shell_rhf_retirement.md); heading `Complete-Core-Shell RHF Retirement`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** remove stale RHF payload stack.

### HP-RETIRE-CCS-RHF-TEST-01 - retirement validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [complete\_core\_shell\_rhf\_retirement.md](complete_core_shell_rhf_retirement.md); heading `Complete-Core-Shell RHF Retirement`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** retirement validation.

### HP-RETIRE-DRV-MAT-DOC-01 - active docs cleanup

- **Lifecycle:** `retired`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [route\_driver\_materialization\_retirement.md](route_driver_materialization_retirement.md); heading `Route-Driver Materialization Workflow Retirement`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** active docs cleanup.

### HP-RETIRE-DRV-MAT-FN-01 - remove old materialization/report/save wrappers

- **Lifecycle:** `retired`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [route\_driver\_materialization\_retirement.md](route_driver_materialization_retirement.md); heading `Route-Driver Materialization Workflow Retirement`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** remove old materialization/report/save wrappers.

### HP-RETIRE-DRV-MAT-TEST-01 - retirement validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [route\_driver\_materialization\_retirement.md](route_driver_materialization_retirement.md); heading `Route-Driver Materialization Workflow Retirement`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** retirement validation.

### HP-RETIRE-DRV-MAT-TOOL-01 - old wrapper-tool quarantine

- **Lifecycle:** `retired`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [route\_driver\_materialization\_retirement.md](route_driver_materialization_retirement.md); heading `Route-Driver Materialization Workflow Retirement`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** old wrapper-tool quarantine.

### HP-RETIRE-LADDER-RUNNERS-FN-01 - delete dangling ladder runners

- **Lifecycle:** `retired`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [route\_driver\_materialization\_retirement.md](route_driver_materialization_retirement.md); heading `Route-Driver Materialization Workflow Retirement`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** delete dangling ladder runners.

### HP-RETIRE-LADDER-RUNNERS-TEST-01 - ladder runner deletion validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [route\_driver\_materialization\_retirement.md](route_driver_materialization_retirement.md); heading `Route-Driver Materialization Workflow Retirement`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** ladder runner deletion validation.

### HP-RETIRE-QW-DONOR-FN-01 - retire obsolete QW and high-order experimental source cluster

- **Lifecycle:** `retired`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [qw\_high\_order\_experimental\_retirement.md](qw_high_order_experimental_retirement.md); heading `QW And High-Order Experimental Cluster Retirement`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** Delete the obsolete four-file QW/high-order experimental cluster and its GaussletBases includes, exports, generics, and carried-space submodule surface without compatibility glue.

### HP-RETIRE-QW-DONOR-TEST-01 - validate obsolete QW and high-order cluster retirement

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [qw\_high\_order\_experimental\_retirement.md](qw_high_order_experimental_retirement.md); heading `QW And High-Order Experimental Cluster Retirement`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** Remove only stale cluster expectations if encountered and validate surviving chain/square constructors plus the current WL producer; add no replacement tests.

### HP-RG-CUTOFF-FN-01 - residual occupation cutoff and identity tolerance defaults

- **Lifecycle:** `superseded`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [residual\_gaussian\_orthogonality\_robustness.md](residual_gaussian_orthogonality_robustness.md); heading `Residual Gaussian Orthogonality And Cutoff Policy`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** residual occupation cutoff and identity tolerance defaults.

### HP-RG-CUTOFF-FN-02 - production residual cutoff tightening

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [residual\_gaussian\_orthogonality\_robustness.md](residual_gaussian_orthogonality_robustness.md); heading `Residual Gaussian Orthogonality And Cutoff Policy`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/residual_basis.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** source maintenance.

### HP-RG-CUTOFF-TEST-01 - residual cutoff/tolerance validation

- **Lifecycle:** `superseded`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [residual\_gaussian\_orthogonality\_robustness.md](residual_gaussian_orthogonality_robustness.md); heading `Residual Gaussian Orthogonality And Cutoff Policy`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** `HP-RG-CUTOFF-FN-01`
- **Scope:** residual cutoff/tolerance validation.

### HP-RG-CUTOFF-TEST-02 - production residual cutoff validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [residual\_gaussian\_orthogonality\_robustness.md](residual_gaussian_orthogonality_robustness.md); heading `Residual Gaussian Orthogonality And Cutoff Policy`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:** none
- **Dependencies:** `HP-RG-CUTOFF-FN-02`
- **Scope:** test maintenance, not source authority.

### HP-RG-FILE-01 - Residual Gaussian module files

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [residual\_gaussian\_domain\_module.md](residual_gaussian_domain_module.md); heading `Residual Gaussian Domain Module`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_residual_gaussians/CartesianResidualGaussians.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/augmented_operators.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/mwg_interaction.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/residual_basis.jl`
- **Evidence:**
  - `repo_path`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Dependencies:** none
- **Scope:** source maintenance.

### HP-RG-FN-01 - residual Gaussian basis construction

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [residual\_gaussian\_domain\_module.md](residual_gaussian_domain_module.md); heading `Residual Gaussian Domain Module`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_residual_gaussians/residual_basis.jl`
- **Evidence:**
  - `repo_path`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Dependencies:** none
- **Scope:** source maintenance.

### HP-RG-FN-02 - exact augmented operator transformation

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [residual\_gaussian\_domain\_module.md](residual_gaussian_domain_module.md); heading `Residual Gaussian Domain Module`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_residual_gaussians/augmented_operators.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** source maintenance.

### HP-RG-FN-03 - moment-matched Gaussian descriptors

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [residual\_gaussian\_domain\_module.md](residual_gaussian_domain_module.md); heading `Residual Gaussian Domain Module`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_residual_gaussians/mwg_interaction.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** source maintenance.

### HP-RG-FN-04 - residual IDA interaction assembly

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [residual\_gaussian\_domain\_module.md](residual_gaussian_domain_module.md); heading `Residual Gaussian Domain Module`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_residual_gaussians/mwg_interaction.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** source maintenance.

### HP-RG-IDTOL-FN-01 - residual final-identity tolerance default

- **Lifecycle:** `superseded`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [residual\_gaussian\_orthogonality\_robustness.md](residual_gaussian_orthogonality_robustness.md); heading `Residual Gaussian Orthogonality And Cutoff Policy`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** `HP-RG-ORTHO-FN-01`
- **Scope:** historical \`identity\_atol = 1e-8\` transition, superseded by \`HP-RG-CUTOFF-FN-01\`; no execution authority.

### HP-RG-IDTOL-TEST-01 - residual final-identity tolerance validation

- **Lifecycle:** `superseded`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [residual\_gaussian\_orthogonality\_robustness.md](residual_gaussian_orthogonality_robustness.md); heading `Residual Gaussian Orthogonality And Cutoff Policy`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** `HP-RG-IDTOL-FN-01`
- **Scope:** residual final-identity tolerance validation.

### HP-RG-INJECT-AUDIT-01 - direct-G injection measurement audit

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [residual\_gaussian\_injection\_hybrid.md](residual_gaussian_injection_hybrid.md); heading `Default-Off Direct-G Residual Injection`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** direct-G injection measurement audit.

### HP-RG-INJECT-FN-01 - default-off direct-G injection compatibility

- **Lifecycle:** `implemented`
- **Grant:** `preservation`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [residual\_gaussian\_injection\_hybrid.md](residual_gaussian_injection_hybrid.md); heading `Default-Off Direct-G Residual Injection`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/augmented_operators.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/mwg_interaction.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/residual_basis.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** source maintenance for existing default-off behavior only; no feature expansion or new caller.

### HP-RG-NUMCOMP-FN-01 - numerical-complete residual basis and additive consumer

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [numerical\_complete\_residual\_basis.md](numerical_complete_residual_basis.md); heading `Numerical-Complete Residual Gaussian Basis`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_protected_ladder_bundle.jl`
  - `source` / `existing`: `src/cartesian_reference_density/atomic_hf_reference_packets.jl`
  - `source` / `existing`: `src/cartesian_reference_density/screened_hartree_correction.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/augmented_operators.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/residual_basis.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** source maintenance for the fixed \`eta\_num = 1e-10\` numerical- complete composition and existing private additive consumer only.

### HP-RG-NUMCOMP-TEST-01 - numerical-complete residual validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [numerical\_complete\_residual\_basis.md](numerical_complete_residual_basis.md); heading `Numerical-Complete Residual Gaussian Basis`
- **Owned paths:**
  - `test` / `existing`: `test/misc/runtests.jl`
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** validation maintenance for the existing compact rank/malformed- metric and bounded H2 surfaces; ignored H2/Be2 and one gated Cr2 comparison remain measurement-only.

### HP-RG-OBJ-01 - residual Gaussian basis object

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [residual\_gaussian\_domain\_module.md](residual_gaussian_domain_module.md); heading `Residual Gaussian Domain Module`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_residual_gaussians/residual_basis.jl`
- **Evidence:** none
- **Dependencies:** `HP-RG-FILE-01`
- **Scope:** source maintenance.

### HP-RG-OCC-FIRST-INJECT-AUDIT-01 - occupied-first global injection measurement audit

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [occupied\_first\_injection.md](occupied_first_injection.md); heading `Occupied-First Injection Geometry`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** occupied-first global injection measurement audit.

### HP-RG-OCC-FIRST-INJECT-FN-01 - occupied-first injection geometry

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [occupied\_first\_injection.md](occupied_first_injection.md); heading `Occupied-First Injection Geometry`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_residual_gaussians/residual_basis.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** validate physical capture geometry, make supplied \`Y\_occ\` mandatory, keep pre-inclusion capture distinct from post-inclusion recovery, and capture-select optional supplement directions. Weak rejected directions never become MWG residual channels.

### HP-RG-OCC-FIRST-INJECT-TEST-01 - occupied-first validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [occupied\_first\_injection.md](occupied_first_injection.md); heading `Occupied-First Injection Geometry`
- **Owned paths:**
  - `test` / `existing`: `test/misc/runtests.jl`
  - `test` / `existing`: `test/nested/cartesian_occupied_first_injection_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain synthetic pre/post and malformed-capture checks plus the bounded real packet-driven Be/Ne PQS gate and terminal due diligence.

### HP-RG-ORTHO-FN-01 - residual final-orthogonality robustness

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [residual\_gaussian\_orthogonality\_robustness.md](residual_gaussian_orthogonality_robustness.md); heading `Residual Gaussian Orthogonality And Cutoff Policy`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/residual_basis.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** source maintenance.

### HP-RG-ORTHO-TEST-01 - residual final-orthogonality validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [residual\_gaussian\_orthogonality\_robustness.md](residual_gaussian_orthogonality_robustness.md); heading `Residual Gaussian Orthogonality And Cutoff Policy`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:** none
- **Dependencies:** `HP-RG-ORTHO-FN-01`
- **Scope:** test maintenance, not source authority.

### HP-RG-PROTECT-ADDREF-FN-01 - protected additive atomic reference correction

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [protected\_additive\_reference\_correction.md](protected_additive_reference_correction.md); heading `Protected Additive Atomic Reference Correction`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl`
  - `source` / `existing`: `src/cartesian_protected_ladder_bundle.jl`
  - `source` / `existing`: `src/cartesian_reference_density/atomic_hf_reference_packets.jl`
  - `source` / `existing`: `src/cartesian_reference_density/screened_hartree_correction.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/augmented_operators.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/residual_basis.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** build compact \`R\` once; use staged protected geometry with the full-rank occupied union mandatory for basis protection; preserve original per-packet occupied blocks for additive \`P0\`; build placed fitted-potential \`GG/GA/AA\`; include all self and twice-cross \`E0\` terms; transform \`J0\` through native protected/localized one-body operators; and return the existing in-memory \`ScreenedHartreeCorrection\` plus reference diagnostics. The private seam returns \`(member, correction, reference)\`; the no-reference path remains unchanged.

### HP-RG-PROTECT-ADDREF-TEST-01 - additive-reference validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [protected\_additive\_reference\_correction.md](protected_additive_reference_correction.md); heading `Protected Additive Atomic Reference Correction`
- **Owned paths:**
  - `test` / `existing`: `test/misc/runtests.jl`
  - `test` / `existing`: `test/nested/cartesian_screened_hartree_correction_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain mandatory-union recovery, packet embedding/failure, per-packet trace, additive \`P0/q0\`, self/cross \`E0\`, placed raw-block, protected/localized \`J0\`, correction-anchor, no-reference parity, and ordinary-packet rejection/diagnostic checks. The accepted padded Be2 smoke is structural evidence only; its retired polish-assisted energy value is not a current endpoint gate.

### HP-RG-PROTECT-ART-FN-01 - protected-localized Hamiltonian artifact variant

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `artifacts`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [protected\_localized\_artifact.md](protected_localized_artifact.md); heading `Protected-Localized Artifact Contract`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_ida_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/augmented_operators.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** write and read the recognized, versioned, opt-in protected-localized artifact without changing its native matrix order.

### HP-RG-PROTECT-ART-TEST-01 - protected artifact validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [protected\_localized\_artifact.md](protected_localized_artifact.md); heading `Protected-Localized Artifact Contract`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** `HP-RG-PROTECT-ART-FN-01`
- **Scope:** protected artifact validation.

### HP-RG-PROTECT-ARTLOC-FN-01 - protected artifact row-locality metadata

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [protected\_localized\_artifact.md](protected_localized_artifact.md); heading `Protected-Localized Artifact Contract`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_ida_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/augmented_operators.jl`
- **Evidence:** none
- **Dependencies:** `HP-RG-PROTECT-ART-FN-01`
- **Scope:** attach validated native-order center, sector, inverse-permutation, and optional all-or-none spread metadata to the protected artifact.

### HP-RG-PROTECT-ARTLOC-TEST-01 - row-locality validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [protected\_localized\_artifact.md](protected_localized_artifact.md); heading `Protected-Localized Artifact Contract`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** `HP-RG-PROTECT-ART-TEST-01`, `HP-RG-PROTECT-ARTLOC-FN-01`
- **Scope:** row-locality validation.

### HP-RG-PROTECT-EGOI-AUDIT-01 - protected-localized EGOI measurement audit

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [retained\_gto\_egoi.md](retained_gto_egoi.md); heading `Retained-GTO Local-Product EGOI`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** protected-localized EGOI measurement audit.

### HP-RG-PROTECT-EGOI-FN-01 - retained-GTO local-product EGOI helper

- **Lifecycle:** `approved`
- **Grant:** `implementation`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [retained\_gto\_egoi.md](retained_gto_egoi.md); heading `Retained-GTO Local-Product EGOI`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_residual_gaussians/augmented_operators.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/residual_basis.jl`
  - `source` / `existing`: `src/hamiltonian_corrections.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** build the owner-balanced retained-original \`s1+s2\` target, native \`Qtarget\`, local symmetric products, \`AA-AA\` / \`BB-BB\` / \`AA-BB\` acceptance blocks, exactly local \`M2\` \`DeltaV\`, and compact diagnostics in memory.

### HP-RG-PROTECT-EGOI-TEST-01 - retained-GTO EGOI validation

- **Lifecycle:** `approved`
- **Grant:** `implementation`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [retained\_gto\_egoi.md](retained_gto_egoi.md); heading `Retained-GTO Local-Product EGOI`
- **Owned paths:**
  - `test` / `planned`: `test/nested/cartesian_retained_gto_egoi_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** test implementation on one exact planned path.

### HP-RG-PROTECT-INJECT-DESIGN-01 - protected-original compact-main injection design

- **Lifecycle:** `completed`
- **Grant:** `design`
- **Surfaces:** `docs`
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [protected\_localized\_basis.md](protected_localized_basis.md); heading `Protected-Localized Basis Convention`
- **Owned paths:**
  - `docs` / `existing`: `docs/src/developer/designs/cartesian_hamiltonian_producer/protected_localized_basis.md`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** design-only rationale maintenance; no source, test, artifact, or workflow execution authority.

### HP-RG-PROTECT-INJECT-FN-01 - staged protected-original geometry

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [protected\_localized\_basis.md](protected_localized_basis.md); heading `Protected-Localized Basis Convention`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_residual_gaussians/residual_basis.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** consume the already-built compact residual object; build protected and broad original subspaces; keep Gaussian Gram, representability, and fake-RDM gates distinct; and return transform-ready \`Z\`, \`B\`, \`Q\_perp\`, \`F\`, and diagnostics. Rejected broad directions never become MWG channels.

### HP-RG-PROTECT-INJECT-TEST-01 - protected geometry validation

- **Lifecycle:** `completed`
- **Grant:** `measurement`
- **Surfaces:** `docs`, `measurement`
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [protected\_localized\_basis.md](protected_localized_basis.md); heading `Protected-Localized Basis Convention`
- **Owned paths:**
  - `docs` / `existing`: `docs/src/developer/reports/cr2_staged_subspace_filter_870498b54/README.md`
  - `measurement` / `optional_local`: `tmp/work/cr2_source_backed_staged_protected_geometry_probe.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain only the exact optional geometry probe and tracked report as completed evidence for \`HP-RG-PROTECT-INJECT-FN-01\`; no committed test or execution-whitelist authority.

### HP-RG-PROTECT-LADDER-BUNDLE-FN-01 - protected-localized ladder bundle facility

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [protected\_localized\_ladder.md](protected_localized_ladder.md); heading `Protected-Localized Ladder Bundles`
- **Owned paths:**
  - `source` / `existing`: `src/GaussletBases.jl`
  - `source` / `existing`: `src/cartesian_ida_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_protected_ladder_bundle.jl`
  - `source` / `existing`: `src/cartesian_representation_transfer.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** write a versioned directory manifest, protected member artifacts, exact adjacent \`S\_BA\` sidecars, optional native-order restart sidecars, and bounded summaries; transfer only as \`C\_B = S\_BA \* C\_A\` and evaluate with target \`H1\_L\` / \`Vee\_L\`.

### HP-RG-PROTECT-LADDER-BUNDLE-TEST-01 - ladder bundle validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [protected\_localized\_ladder.md](protected_localized_ladder.md); heading `Protected-Localized Ladder Bundles`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** ladder bundle validation.

### HP-RG-PROTECT-LADDER-XFER-AUDIT-01 - same-parent ladder transfer audit

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [protected\_localized\_ladder.md](protected_localized_ladder.md); heading `Protected-Localized Ladder Bundles`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** same-parent ladder transfer audit.

### HP-RG-PROTECT-ONEBODY-AUDIT-01 - protected fixed-sector one-body audit

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [protected\_localized\_basis.md](protected_localized_basis.md); heading `Protected-Localized Basis Convention`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** protected fixed-sector one-body audit.

### HP-RG-PROTECT-ONEBODY-FN-01 - exact protected one-body transform

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [protected\_localized\_basis.md](protected_localized_basis.md); heading `Protected-Localized Basis Convention`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_residual_gaussians/augmented_operators.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** construct dense exact fixed-sector kinetic, per-center unit nuclear, and assembled \`H1\_F\` matrices through the actual protected/localized one-body transform, with orthogonality and symmetry diagnostics.

### HP-RG-PROTECT-ONEBODY-TEST-01 - protected one-body validation

- **Lifecycle:** `completed`
- **Grant:** `measurement`
- **Surfaces:** `docs`, `measurement`
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [protected\_localized\_basis.md](protected_localized_basis.md); heading `Protected-Localized Basis Convention`
- **Owned paths:**
  - `docs` / `existing`: `docs/src/developer/reports/cr2_protected_onebody_audit_eaf05a38c/README.md`
  - `measurement` / `optional_local`: `tmp/work/cr2_protected_onebody_dense_source_replay.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain only the exact optional dense replay and tracked report as completed evidence for \`HP-RG-PROTECT-ONEBODY-FN-01\`; no committed test or execution-whitelist authority.

### HP-RG-PROTECT-VEE-AUDIT-01 - protected interaction decision audit

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [protected\_localized\_basis.md](protected_localized_basis.md); heading `Protected-Localized Basis Convention`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** protected interaction decision audit.

### HP-RG-RHO0-GAL-AUDIT-01 - row-gauge rho0 audit

- **Lifecycle:** `superseded`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [rho0\_reference\_density\_matrix.md](rho0_reference_density_matrix.md); heading `Rho0 And Reference-Density Correction History`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** row-gauge rho0 audit.

### HP-RG-SPECTRAL-AUDIT-01 - residual-sector spectral audit

- **Lifecycle:** `approved`
- **Grant:** `measurement`
- **Surfaces:** `measurement`
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [residual\_gaussian\_orthogonality\_robustness.md](residual_gaussian_orthogonality_robustness.md); heading `Residual Gaussian Orthogonality And Cutoff Policy`
- **Owned paths:**
  - `measurement` / `optional_local`: `tmp/work/rg_spectral_cutoff1e6_audit.jl`
- **Evidence:**
  - `external_path`: `/Users/srw/dmrgtmp/cr2_r1p68_ns7_lmax2_d0p00847_fixed95fec2b8/rg_spectral_cutoff1e6_modes.tsv`
  - `external_path`: `/Users/srw/dmrgtmp/cr2_r1p68_ns7_lmax2_d0p00847_fixed95fec2b8/rg_spectral_cutoff1e6_owner_metrics.tsv`
  - `external_path`: `/Users/srw/dmrgtmp/cr2_r1p68_ns7_lmax2_d0p00847_fixed95fec2b8/rg_spectral_cutoff1e6_spectra.tsv`
  - `external_path`: `/Users/srw/dmrgtmp/cr2_r1p68_ns7_lmax2_d0p00847_fixed95fec2b8/rg_spectral_cutoff1e6_stage_timings.tsv`
  - `external_path`: `/Users/srw/dmrgtmp/cr2_r1p68_ns7_lmax2_d0p00847_fixed95fec2b8/rg_spectral_cutoff1e6_summary.txt`
- **Dependencies:** none
- **Scope:** run only the exact ignored residual-sector spectral probe and external text/TSV reporting governed by the canonical contract; no tracked source, test, artifact, driver, or workflow changes.

### HP-RG-TEST-01 - Residual Gaussian endpoint validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [residual\_gaussian\_domain\_module.md](residual_gaussian_domain_module.md); heading `Residual Gaussian Domain Module`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** test maintenance, not source authority.

### HP-RG-WIRE-01 - terminal and facade compatibility wiring

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `driver`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [residual\_gaussian\_domain\_module.md](residual_gaussian_domain_module.md); heading `Residual Gaussian Domain Module`
- **Owned paths:**
  - `driver` / `existing`: `bin/cartesian_ham_builder.jl`
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
- **Evidence:** none
- **Dependencies:** `HP-RG-FN-01`, `HP-RG-FN-02`, `HP-RG-FN-03`, `HP-RG-FN-04`
- **Scope:** compatibility/caller maintenance.

### HP-RHO0-ANCHOR-FN-01 - old full-interaction anchor

- **Lifecycle:** `superseded`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [rho0\_reference\_density\_matrix.md](rho0_reference_density_matrix.md); heading `Rho0 And Reference-Density Correction History`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** old full-interaction anchor.

### HP-RHO0-ANCHOR-TEST-01 - old anchor validation

- **Lifecycle:** `superseded`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [rho0\_reference\_density\_matrix.md](rho0_reference_density_matrix.md); heading `Rho0 And Reference-Density Correction History`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** old anchor validation.

### HP-RHO0-CORR-AUDIT-01 - corrected-Hamiltonian audit

- **Lifecycle:** `superseded`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [rho0\_reference\_density\_matrix.md](rho0_reference_density_matrix.md); heading `Rho0 And Reference-Density Correction History`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** corrected-Hamiltonian audit.

### HP-RHO0-FAPP-AUDIT-01 - approximate fixed-P0 Fock audit

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [rho0\_reference\_density\_matrix.md](rho0_reference_density_matrix.md); heading `Rho0 And Reference-Density Correction History`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** approximate fixed-P0 Fock audit.

### HP-RHO0-FAPP-FN-01 - approximate IDA energy/Fock seam

- **Lifecycle:** `implemented`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [rho0\_reference\_density\_matrix.md](rho0_reference_density_matrix.md); heading `Rho0 And Reference-Density Correction History`
- **Owned paths:** none
- **Evidence:**
  - `repo_path`: `src/cartesian_ida_hamiltonian.jl`
- **Dependencies:** none
- **Scope:** approximate IDA energy/Fock seam.

### HP-RHO0-FAPP-TEST-01 - approximate IDA derivative validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [rho0\_reference\_density\_matrix.md](rho0_reference_density_matrix.md); heading `Rho0 And Reference-Density Correction History`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** approximate IDA derivative validation.

### HP-RHO0-JANCHOR-FN-01 - direct-Hartree anchor helper

- **Lifecycle:** `superseded`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [rho0\_reference\_density\_matrix.md](rho0_reference_density_matrix.md); heading `Rho0 And Reference-Density Correction History`
- **Owned paths:** none
- **Evidence:**
  - `repo_path`: `src/cartesian_ida_hamiltonian.jl`
  - `repo_path`: `src/cartesian_residual_gaussians/augmented_operators.jl`
- **Dependencies:** none
- **Scope:** direct-Hartree anchor helper.

### HP-RHO0-JANCHOR-TEST-01 - direct-Hartree anchor validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [rho0\_reference\_density\_matrix.md](rho0_reference_density_matrix.md); heading `Rho0 And Reference-Density Correction History`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** direct-Hartree anchor validation.

### HP-RHO0-MIXH-FEXACT-FN-01 - protected exact-Hartree transform

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [reference\_hartree\_numerics.md](reference_hartree_numerics.md); heading `Reference Hartree Numerics`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_residual_gaussians/augmented_operators.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** source maintenance for the exact protected fixed-sector transform, localized \`sym(W' \* J0\_F \* W)\`, and existing convenience composition only.

### HP-RHO0-MIXH-FEXACT-TEST-01 - protected exact-Hartree validation

- **Lifecycle:** `completed`
- **Grant:** `measurement`
- **Surfaces:** `measurement`
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [reference\_hartree\_numerics.md](reference_hartree_numerics.md); heading `Reference Hartree Numerics`
- **Owned paths:**
  - `measurement` / `optional_local`: `tmp/work/rho0_mixh_fexact_source_validation.jl`
- **Evidence:**
  - `git_commit`: `40a6f7e99`
  - `manager_pass`: `288`
- **Dependencies:** none
- **Scope:** maintain only the exact optional bounded H/Be/Be2 replay as completed evidence for \`HP-RHO0-MIXH-FEXACT-FN-01\`; no committed test authority.

### HP-RHO0-MIXH-GAAA-FN-01 - exact mixed-Hartree GA/AA blocks

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [reference\_hartree\_numerics.md](reference_hartree_numerics.md); heading `Reference Hartree Numerics`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl`
  - `source` / `existing`: `src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** source maintenance for exact \`GA\`, symmetric \`AA\`, and compact diagnostics only.

### HP-RHO0-MIXH-GAAA-TEST-01 - exact mixed-Hartree GA/AA validation

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [reference\_hartree\_numerics.md](reference_hartree_numerics.md); heading `Reference Hartree Numerics`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_screened_hartree_correction_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** inherited committed-test maintenance for the placed-potential \`A-A\` kernel checks only.

### HP-RHO0-MIXH-GG-FN-01 - exact mixed-Hartree GG block

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [reference\_hartree\_numerics.md](reference_hartree_numerics.md); heading `Reference Hartree Numerics`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl`
  - `source` / `existing`: `src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl`
  - `source` / `existing`: `src/gaussian_coulomb_reference.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** source maintenance for exact one-center finite symmetric-\`P\_A\` \`GG\` construction and compact diagnostics only.

### HP-RHO0-MIXH-GG-TEST-01 - exact mixed-Hartree GG validation

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [reference\_hartree\_numerics.md](reference_hartree_numerics.md); heading `Reference Hartree Numerics`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_screened_hartree_correction_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** validation maintenance for the stated compact gate and existing consumer test only.

### HP-RHO0-REFDENS-AUDIT-01 - fixed-P0 audit

- **Lifecycle:** `superseded`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [rho0\_reference\_density\_matrix.md](rho0_reference_density_matrix.md); heading `Rho0 And Reference-Density Correction History`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** fixed-P0 audit.

### HP-RHO0-REFDENS-ERI-01 - historical candidate mixed-ERI owner

- **Lifecycle:** `superseded`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [rho0\_reference\_density\_matrix.md](rho0_reference_density_matrix.md); heading `Rho0 And Reference-Density Correction History`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** historical candidate mixed-ERI owner.

### HP-RHO0-REFDENS-FN-01 - historical candidate correction owner

- **Lifecycle:** `superseded`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [rho0\_reference\_density\_matrix.md](rho0_reference_density_matrix.md); heading `Rho0 And Reference-Density Correction History`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** historical candidate correction owner.

### HP-RHO0-REFDENS-MIXH-AUDIT-01 - exact mixed-Hartree seam audit

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [rho0\_reference\_density\_matrix.md](rho0_reference_density_matrix.md); heading `Rho0 And Reference-Density Correction History`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** exact mixed-Hartree seam audit.

### HP-RHO0-XPAIR-AUDIT-01 - exchange/direct pairing question

- **Lifecycle:** `deferred`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [rho0\_reference\_density\_matrix.md](rho0_reference_density_matrix.md); heading `Rho0 And Reference-Density Correction History`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** measurement-only after explicit design-manager reactivation, for ignored H/Be/Be2 diagnostics only.

### HP-ROUTE-INV-FN-01 - retained-unit route inventory

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [route\_stage\_metadata\_contract.md](route_stage_metadata_contract.md); heading `Route/Stage Metadata Contract`
- **Owned paths:**
  - `source` / `existing`: `src/pqs_source_box_route_driver_helpers.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain vector-backed ordered retained-unit and pair-family rows with label lookup; labels remain data rather than concrete type parameters.

### HP-ROUTE-INV-TEST-01 - route inventory validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [route\_stage\_metadata\_contract.md](route_stage_metadata_contract.md); heading `Route/Stage Metadata Contract`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** route inventory validation.

### HP-ROUTE-RECIPE-FN-01 - family-selective route recipes

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [nesting\_supplement\_composition\_plan.md](nesting_supplement_composition_plan.md); heading `Nesting/Supplement Composition`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/pqs_source_box_route_driver_helpers.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** build only the selected \`:pqs\_source\_box\` or \`:white\_lindsey\_low\_order\` subrecipe and leave inactive family vocabulary absent or \`nothing\` without merging the algorithms.

### HP-ROUTE-RECIPE-TEST-01 - route recipe validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [nesting\_supplement\_composition\_plan.md](nesting_supplement_composition_plan.md); heading `Nesting/Supplement Composition`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** route recipe validation.

### HP-ROUTE-STAGE-CARRIER-FN-01 - route/stage carriers

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [route\_stage\_metadata\_contract.md](route_stage_metadata_contract.md); heading `Route/Stage Metadata Contract`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
  - `source` / `existing`: `src/pqs_source_box_diatomic_complete_core_shell.jl`
  - `source` / `existing`: `src/pqs_source_box_route_driver_helpers.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** keep only live compact plans/realizations/summaries across stage boundaries. This does not retire route skeletons or pair/assembly/report stages.

### HP-ROUTE-STAGE-CARRIER-TEST-01 - carrier validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [route\_stage\_metadata\_contract.md](route_stage_metadata_contract.md); heading `Route/Stage Metadata Contract`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** carrier validation.

### HP-ROUTE-STAGE-TYPE-FN-01 - route/stage type surface

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [route\_stage\_metadata\_contract.md](route_stage_metadata_contract.md); heading `Route/Stage Metadata Contract`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_terminal_shellification_geometry.jl`
  - `source` / `existing`: `src/pqs_source_box_route_driver_helpers.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** preserve compact vector-backed route/shellification summaries and narrow stage returns without duplicate lowering-plan ownership.

### HP-ROUTE-STAGE-TYPE-TEST-01 - type-surface validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [route\_stage\_metadata\_contract.md](route_stage_metadata_contract.md); heading `Route/Stage Metadata Contract`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** type-surface validation.

### HP-TEST-01 - new committed terminal smoke — rejected

- **Lifecycle:** `rejected`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [terminal\_basis\_and\_base\_assembly.md](terminal_basis_and_base_assembly.md); heading `Terminal Basis And Base Assembly`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** new committed terminal smoke — rejected.

### HP-WIRE-01 - terminal-basis stage integration

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [terminal\_basis\_and\_base\_assembly.md](terminal_basis_and_base_assembly.md); heading `Terminal Basis And Base Assembly`
- **Owned paths:**
  - `source` / `existing`: `src/pqs_source_box_route_driver_helpers.jl`
- **Evidence:**
  - `repo_path`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
- **Dependencies:** `HP-FN-01`
- **Scope:** pass current PQS support, retained, transform, and bundle objects directly into the terminal realizer.

### HP-WIRE-02 - historical direct materialization Hamiltonian handoff

- **Lifecycle:** `retired`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [route\_driver\_materialization\_retirement.md](route_driver_materialization_retirement.md); heading `Route-Driver Materialization Workflow Retirement`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** historical direct materialization Hamiltonian handoff.

### HP-WLDIAT-COMPACT-FN-01 - WL diatomic compact retained basis

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [white\_lindsey\_terminal\_basis\_realization.md](white_lindsey_terminal_basis_realization.md); heading `White-Lindsey Terminal Basis Realization`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`
  - `source` / `existing`: `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`
  - `source` / `existing`: `src/cartesian_retained_units/lower_contract_units.jl`
  - `source` / `existing`: `src/cartesian_shellification/terminal_geometry.jl`
  - `source` / `existing`: `src/cartesian_terminal_lowering/region_contracts.jl`
  - `source` / `existing`: `src/pqs_source_box_route_driver_helpers.jl`
- **Evidence:** none
- **Dependencies:** `HP-COMP-FACEPROD-FN-01`
- **Scope:** maintain compact products of one-dimensional contractions for WL boundary units. Identity realization remains valid only for true direct/core units; support rows are not themselves retained functions.

### HP-WLDIAT-COMPACT-TEST-01 - WL compact-basis validation

- **Lifecycle:** `completed`
- **Grant:** `measurement`
- **Surfaces:** `measurement`
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [white\_lindsey\_terminal\_basis\_realization.md](white_lindsey_terminal_basis_realization.md); heading `White-Lindsey Terminal Basis Realization`
- **Owned paths:**
  - `measurement` / `optional_local`: `tmp/work/wl_compact_block_size_audit.jl`
  - `measurement` / `optional_local`: `tmp/work/wl_diatomic_base_validation.jl`
- **Evidence:**
  - `manager_pass`: `152`
- **Dependencies:** none
- **Scope:** maintain only the two exact optional WL compact-basis probes as completed evidence for \`HP-WLDIAT-COMPACT-FN-01\`; no committed test authority.

### HP-WLDIAT-PARITY-FN-01 - WL boundary retained-count parity

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [white\_lindsey\_terminal\_basis\_realization.md](white_lindsey_terminal_basis_realization.md); heading `White-Lindsey Terminal Basis Realization`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** keep odd-side enforcement only for direct nucleus-centered cores; WL boundary products retain the requested count. Canonical examples remain \`ns=4 -\> 56\` and \`ns=5 -\> 98\` boundary columns.

### HP-WLDIAT-PARITY-TEST-01 - WL parity validation

- **Lifecycle:** `completed`
- **Grant:** `measurement`
- **Surfaces:** `measurement`
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [white\_lindsey\_terminal\_basis\_realization.md](white_lindsey_terminal_basis_realization.md); heading `White-Lindsey Terminal Basis Realization`
- **Owned paths:**
  - `measurement` / `optional_local`: `tmp/work/wl_compact_block_size_audit.jl`
  - `measurement` / `optional_local`: `tmp/work/wl_diatomic_base_validation.jl`
- **Evidence:**
  - `manager_pass`: `154`
- **Dependencies:** none
- **Scope:** maintain only the two exact optional WL parity probes as completed evidence for \`HP-WLDIAT-PARITY-FN-01\`; no committed test authority.

### HP-WLTERM-FILE-01 - WL terminal realization file

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [white\_lindsey\_terminal\_basis\_realization.md](white_lindsey_terminal_basis_realization.md); heading `White-Lindsey Terminal Basis Realization`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** maintain the existing WL-specific sibling that returns the shared \`CartesianTerminalBasisRealization\`. No new module, basis object, route result, artifact, report, or export is authorized.

### HP-WLTERM-FN-01 - WL low-order terminal basis realization

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [white\_lindsey\_terminal\_basis\_realization.md](white_lindsey_terminal_basis_realization.md); heading `White-Lindsey Terminal Basis Realization`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`
- **Evidence:** none
- **Dependencies:** `HP-COMP-FACEPROD-FN-01`, `HP-WLTERM-FILE-01`
- **Scope:** realize direct identity blocks and compact WL facet/edge/corner, boundary-stratum, and thin-slab products on authoritative owned supports while preserving retained/transform order and block-local identity checks.

### HP-WLTERM-TEST-01 - WL terminal-basis validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [white\_lindsey\_terminal\_basis\_realization.md](white_lindsey_terminal_basis_realization.md); heading `White-Lindsey Terminal Basis Realization`
- **Owned paths:**
  - `test` / `existing`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** inherited committed-test maintenance only.

### HP-WLTERM-WIRE-01 - WL route helper terminal-basis wiring

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [white\_lindsey\_terminal\_basis\_realization.md](white_lindsey_terminal_basis_realization.md); heading `White-Lindsey Terminal Basis Realization`
- **Owned paths:**
  - `source` / `existing`: `src/pqs_source_box_route_driver_helpers.jl`
- **Evidence:** none
- **Dependencies:** none
- **Scope:** pass native \`:white\_lindsey\_low\_order\` support, retained, and transform records into the WL realizer without changing PQS behavior or restoring old WL materialization.
