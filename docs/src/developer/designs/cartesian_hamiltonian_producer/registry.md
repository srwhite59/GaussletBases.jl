# Cartesian Hamiltonian Producer Authority Registry

> **Generated authority view. Do not edit.** The record-level source is
> [authority.toml](authority.toml), SHA-256 `da28f6f9f15fdb79bac8ce0b4120351d8b3aa54d33fb458d6c3b11e7c37710ee`.

Tracked producer work is authorized only when a unique record has an
execution grant and surface, and the requested change stays within its exact
owned paths, scope, `current.md`, `invariants.md`, and canonical contract.
Lifecycle never grants work by itself. Any missing or conflicting fact fails closed.

## Records

### HP-CARTESIAN-INTERNAL-MAINTENANCE-CI-FN-01 - scheduled Cartesian internal numerical-maintenance workflow

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `tools`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [test\_suite\_reorganization\_plan.md](../../test_suite_reorganization_plan.md); heading `Scheduled Cartesian Internal Maintenance Gate`
- **Owned paths:**
  - `tool` / `existing`: `.github/workflows/cartesian-internal-maintenance.yml`
- **Evidence:**
  - `git_commit`: `1f550899ffac6b682081595b482e1bdfb2e65da2`
  - `git_commit`: `b8e89dae08840665847628e2c8b2e8bf42730387`
  - `manager_pass`: `510`
  - `manager_pass`: `512`
  - `manager_pass`: `514`
- **Dependencies:** `HP-PQS-ATOMREF-PACKET-TEST-01`, `HP-PQS-SCREEN-HARTREE-CORR-TEST-01`, `HP-R3-TEST-01`, `HP-REP-XGTO-PROTECT-SIDECAR-TEST-01`
- **Scope:** Maintain exactly the implemented \`.github/workflows/cartesian-internal-maintenance.yml\`: Julia 1.12, weekly cron \`17 10 \* \* 3\` plus \`workflow\_dispatch\`, one sequential job with job-level \`contents: read\`, one instantiate and package-load check, then atomic packet, protected sidecar, R3A, and screened Hartree as four separate processes in that order, with \`timeout-minutes: 15\`. Commit \`b8e89dae08840665847628e2c8b2e8bf42730387\` changed only the timeout, and manual run \`32657600781\` passed in 9m59s. Maintenance is limited to semantics-preserving workflow syntax or action-version upkeep. No command, test composition, ordering, cadence, trigger, permission, Julia version, timeout, public gate, source, test, RC1, or release-state change is authorized without a new amendment.

### HP-CARTESIAN-INTERNAL-MAINTENANCE-CI-TEST-01 - scheduled internal gate and screening deduplication validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [test\_suite\_reorganization\_plan.md](../../test_suite_reorganization_plan.md); heading `Scheduled Cartesian Internal Maintenance Gate`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_screened_hartree_correction_runtests.jl`
- **Evidence:**
  - `git_commit`: `1f550899ffac6b682081595b482e1bdfb2e65da2`
  - `manager_pass`: `510`
  - `manager_pass`: `512`
- **Dependencies:** `HP-CARTESIAN-INTERNAL-MAINTENANCE-CI-FN-01`, `HP-PQS-ATOMREF-PACKET-TEST-01`, `HP-PQS-PUBLIC-SCREEN-TEST-01`, `HP-PQS-SCREEN-HARTREE-CORR-TEST-01`
- **Scope:** Maintain the completed screening deduplication from commit \`1f550899ffac6b682081595b482e1bdfb2e65da2\`: 47 duplicate lines were removed while the remaining screened-Hartree owner passes 85 checks. Generic supplied-field malformed cases remain owned by the 22-check public screening test, and packet readback/schema failures remain owned by the 117-check atomic-packet test. Preserve packet correction, occupied embedding, additive P0/q0 and self/cross accounting, translated density fields, GG/GA/AA checks, fitted-field validation, and consumer-specific rejection. The timeout repair authorizes no test edit, assertion, fixture, file, runner row, or other nested-suite change. Occupied-first injection and represented molecular Hartree remain unscheduled and unclassified.

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
- **Scope:** maintain explicit equal-symbol/equal-charge diatomics at two finite distinct z-axis centers and send them through the existing PQS/WL base path. Electron-sector independence and charged-sector acceptance are owned separately by HP-R1-ESECTOR-FN-01.

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
- **Scope:** maintain \`Natom=1\`, \`basisname === nothing\` base-atom selection; explicit origin, nuclear charge, \`nup\`/\`ndn\`, \`ns\`, \`core\_spacing\`, \`s\_factor\`, source-span/nesting, and radius-from-padding inputs; and clear unsupported-input failures. The driver infers neither neutrality nor spin. There is no public \`mode=:base\` input.

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
- **Evidence:**
  - `external_path`: `/Users/srw/dmrgtmp/pqs_perf_audit_f9ca5a7b_20260825/compile_and_nuclear_followup.md`
  - `git_commit`: `3419da6132810d8c4454f5b013c6302ef7842cb3`
  - `manager_pass`: `532`
  - `manager_pass`: `533`
- **Dependencies:** `HP-FN-03`
- **Scope:** Maintain the separate non-exported working-basis, product/moment, unit-nuclear, IDA/MWG interaction, residual-augmentation, and Hamiltonian-assembly stages. The base kinetic and unit-nuclear stages must reuse their accepted pre-sized lexical terminal action/tile/block buffers across the three Cartesian product terms and physical centers while preserving exact matrices, term/block/center order, 64 MiB tiling, stage return shapes, and all endpoint facts. Facades remain wrappers over the same construction. Do not add a stage/API/result/metadata change, public control, workspace type, global/task-local/persistent cache, escaping scratch state, new source owner, or arithmetic reorder.

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
- **Evidence:**
  - `git_commit`: `e1d5ca2ddb3a39134fddb476d00029ec590c431f`
  - `manager_pass`: `526`
- **Dependencies:** none
- **Scope:** Maintain the implemented direct support-local shell-seed assembly in \`pqs\_terminal\_basis\_realization.jl\`. \`\_shell\_seed\` must assemble exactly \`support.support\_states\` by retained boundary modes in existing order for ordinary and carried-axis-fact inputs, with the established Float64 conversion points, literal \`0.0\` for skipped zero products, and exact \`vx \* vy \* vz\` arithmetic. Preserve every support/index, retained mode, column range, Gram/Lowdin input, topology, provenance, due-diligence fact, and public result. Do not restore \`\_shell\_seed\_full\_coefficients\_from\_axis\_facts\` or \`\_shell\_seed\_full\_coefficients\`, add a second construction path, or add a fallback, workspace, cache, type, API, metadata, dependency, or validation-policy change. Cold-compilation specialization, scalar-loop optimization, Gram-policy changes, and compatibility cleanup require separate authority.

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
  - `external_path`: `/Users/srw/dmrgtmp/pqs_perf_audit_f9ca5a7b_20260825/compile_and_nuclear_followup.md`
  - `git_commit`: `3419da6132810d8c4454f5b013c6302ef7842cb3`
  - `git_commit`: `94ec277d954b5435a04b0ad68ae352c95b0434c7`
  - `manager_pass`: `532`
  - `manager_pass`: `533`
  - `manager_pass`: `534`
  - `manager_pass`: `535`
  - `repo_path`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
  - `repo_path`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Dependencies:** `HP-OBJ-01`, `HP-OBJ-02`
- **Scope:** Maintain the accepted exact-order four-element reduction in \`\_fill\_terminal\_gaussian\_sum\_action\!\` from commit \`94ec277d954b5435a04b0ad68ae352c95b0434c7\`. One call-local nonescaping Float64 \`coefficients\[term\] \* fx\[term, ix, jx\]\` table is reused across existing block pairs; each complete group contains exactly four independent support-pair outputs with four explicit scalar accumulators and ordinary three-index access; the zero-to-three-element remainder stays scalar. Preserve each output element's original 135-term order, addition order, left-associated \`((coefficients \* fx) \* fy) \* fz\` Float64 arithmetic, conversion points, and bitwise PQS/White-Lindsey matrices for both nuclei. Preserve all dimensions, fingerprints, energies, captures, eigen-residuals, symmetry, topology, warnings, due diligence, release tolerances, call-local ownership, pre-sized buffers, and 64 MiB tiling. The accepted one-file source delta is \`+50/-14\`, exactly at the authorized hard added-line limit; the old scalar loop is deleted. Accepted measurements are PQS \`5.021 -\> 2.687 s\`, White-Lindsey \`4.125 -\> 2.182 s\`, warmed complete comparison \`32.286 -\> 24.020--24.102 s\`, fresh complete comparison \`60.770 -\> 51.911--52.961 s\`, maximum compilation increase \`0.134 s\`, allocation-free scalar accumulation, and \`491,600\` bytes of call-local preweighting. CI run \`33023091521\` and Docs run \`33023091519\` passed. Do not add a fallback, helper, tuple machinery, eight-lane batching, further loop restructuring, explicit full unrolling, term-major traversal, layout-dependent offsets, type, file, API, cache, persistent state, dependency, metadata, test, cold reporting/compilation change, compatibility cleanup, or release work.

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
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
- **Evidence:**
  - `git_commit`: `e1d5ca2ddb3a39134fddb476d00029ec590c431f`
  - `manager_pass`: `526`
- **Dependencies:** none
- **Scope:** Maintain consumption of the three validated materialized carried axis matrices through the same direct support-local shell-seed loop as ordinary source matrices. Preserve carried and ordinary byte-exact coefficients, interval/shape validation, retained boundary order, support restriction, Gram/Lowdin, canonicalization, provenance, and all mapped-COMX behavior. This grants only \`pqs\_terminal\_basis\_realization.jl\`; do not restore global parent-row/source-column materialization or add mapped-source construction, module wiring, a second branch, fallback, test edit, or source-span policy change.

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

### HP-PQS-ASPECTSHELL-FN-01 - matched PQS/WL aspect-aware shell modes

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_complete\_shell\_aspect\_source\_modes.md](pqs_complete_shell_aspect_source_modes.md); heading `Matched PQS/White-Lindsey Complete-Shell Aspect Source Modes`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`
  - `source` / `existing`: `src/pqs_source_box_route_driver_helpers.jl`
- **Evidence:**
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/Papers/PQS/validation/pqs_wl_shared_shell_policy_audit_2026-07-29.md`
  - `git_commit`: `6e2c97704`
  - `manager_pass`: `441`
- **Dependencies:** `HP-COMP-SHELLGEOM-DIAT-FN-01`, `HP-WLDIAT-COMPACT-FN-01`, `HP-WLDIAT-PARITY-FN-01`, `HP-WLTERM-FN-01`
- **Scope:** Correct only eligible shared complete shells so PQS and White-Lindsey consume one post-shellification outer shape \`D=(ns,ns,L)\`: preserve PQS's boundary quotient, make WL derive axis-specific inner counts \`D.-2\`, and require equal aggregate shell and bare terminal dimensions without changing the parent, common shell supports, public \`ns\`, route-local \`q\`, direct cores, slabs, PQS numerics, or ineligible one-center/nonshared WL behavior. Reuse the existing \`source\_mode\_shape\`; the sole anti-bloat exception is one WL-realizer read of that existing metadata key. Add no shape field, route object, generalized metadata algorithm, driver branch, artifact, or interaction change. Limit the correction to 60 added source lines in the approved existing files. Stop without committing if the current strata cannot consume noncubic counts through existing products.

### HP-PQS-ASPECTSHELL-TEST-01 - matched aspect-source validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_complete\_shell\_aspect\_source\_modes.md](pqs_complete_shell_aspect_source_modes.md); heading `Matched PQS/White-Lindsey Complete-Shell Aspect Source Modes`
- **Owned paths:**
  - `test` / `existing`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
- **Evidence:**
  - `manager_pass`: `247`
  - `manager_pass`: `248`
  - `git_commit`: `6e2c97704`
  - `manager_pass`: `441`
- **Dependencies:** `HP-PQS-ASPECTSHELL-FN-01`
- **Scope:** Add at most 30 lines for one bounded matched-H2 regression in the existing public Cartesian base test. Require common parent/shell/support parity, one outer shape per eligible shell, equal PQS boundary and aggregated WL stratum counts, equal bare terminal dimensions, exact PQS parity, WL Gram/column checks, and finite symmetric operators. Complete acceptance with the unchanged clean private paper-driver replay; add no test file, fixture, probe, or committed paper output.

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

### HP-PQS-DOCS-TAGDEPLOY-FN-01 - tag-aware versioned documentation deployment

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `docs`, `tools`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [documentation\_deployment.md](documentation_deployment.md); heading `Tag-Aware Version Deployment`
- **Owned paths:**
  - `docs` / `existing`: `docs/src/developer/designs/cartesian_hamiltonian_producer/current.md`
  - `docs` / `existing`: `docs/src/developer/designs/cartesian_hamiltonian_producer/documentation_deployment.md`
  - `tool` / `existing`: `.github/workflows/docs.yml`
  - `tool` / `existing`: `docs/make.jl`
- **Evidence:**
  - `git_commit`: `2b3c23970144aa030ae52b875a5cf01b32886b6e`
  - `git_commit`: `31caa87d3b83599de7f7295678ee599209113552`
  - `git_commit`: `abee269eed7028c864fa18ae44b4b946af63dfcf`
  - `manager_pass`: `490`
  - `manager_pass`: `491`
  - `manager_pass`: `495`
  - `manager_pass`: `496`
  - `manager_pass`: `519`
  - `manager_pass`: `520`
  - `manager_pass`: `522`
- **Dependencies:** `HP-PQS-PUBLIC-DOC-01`, `HP-PQS-PUBLIC-DOC-PARITY-FN-01`
- **Scope:** Maintain only the implemented tag-aware deployment in the existing Docs workflow and docs/make.jl. Preserve PR build-only contents:read behavior, main-to-dev deployment, exact canonical full vMAJOR.MINOR.PATCH or vMAJOR.MINOR.PATCH-PRERELEASE tag folders without build metadata, fail-closed rejection before deploydocs, and Documenter's standard exclusion of prereleases from the real stable alias. Retain the existing v0.2.0-rc2 and v0.2.0-rc1 folders through their exact self-mappings; neither creates an alias. No other prerelease entry or custom stable policy is authorized. GITHUB\_TOKEN and GAUSSLETBASES\_DOCS\_DEPLOY remain deployment-step-only under the contents:write deployment job. No new file, workflow change, tag mutation, release, source/API/example/scientific-doc/dependency/citation change, custom credential, alternate host, manifest, artifact, arbitrary-tag deployment, dynamic version index, or release framework.

### HP-PQS-DOCS-TAGDEPLOY-TEST-01 - tag-aware documentation deployment validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [documentation\_deployment.md](documentation_deployment.md); heading `Tag-Aware Version Deployment`
- **Owned paths:**
  - `test` / `existing`: `test/docs/runtests.jl`
- **Evidence:**
  - `git_commit`: `2b3c23970144aa030ae52b875a5cf01b32886b6e`
  - `git_commit`: `31caa87d3b83599de7f7295678ee599209113552`
  - `git_commit`: `abee269eed7028c864fa18ae44b4b946af63dfcf`
  - `manager_pass`: `490`
  - `manager_pass`: `491`
  - `manager_pass`: `495`
  - `manager_pass`: `496`
  - `manager_pass`: `519`
  - `manager_pass`: `520`
  - `manager_pass`: `522`
- **Dependencies:** `HP-PQS-DOCS-TAGDEPLOY-FN-01`
- **Scope:** Maintain only the focused existing docs-test assertions for PR/main/prerelease/final/malformed-tag classification, exact main and tag canonical paths, semantic-version rejection, unchanged least-privilege workflow boundaries, exact v0.2.0-rc1 and v0.2.0-rc2 selector retention, and standard prerelease exclusion from stable. Use an isolated Documenter version-expansion fixture to require RC2, RC1, and dev entries, no self-symlink, and no stable alias without a final release. Add no test file, parser framework, actual tag/deployment, numerical gate, source test, or release behavior.

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

### HP-PQS-PAPER-H2-DRV-FN-01 - private matched PQS/WL H2 paper driver

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `driver`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_driver\_usability\_workflow.md](cartesian_driver_usability_workflow.md); heading `Private PQS/WL Paper H2 Driver`
- **Owned paths:**
  - `driver` / `existing`: `bin/pqs_paper_h2_driver.jl`
- **Evidence:**
  - `manager_pass`: `423`
  - `manager_pass`: `424`
  - `manager_pass`: `426`
  - `manager_pass`: `427`
  - `git_commit`: `171e2f368`
  - `manager_pass`: `428`
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/Papers/PQS/provenance/h2plus_R2_full_parent_clean_2026-07-25/`
  - `manager_pass`: `429`
  - `git_commit`: `67351aca6`
  - `manager_pass`: `430`
  - `external_path`: `/Users/srw/Dropbox/Papers/PQS/provenance/h2plus_R2_padding20_clean_2026-07-25/`
  - `manager_pass`: `431`
  - `git_commit`: `33357da4f`
  - `manager_pass`: `432`
  - `external_path`: `/Users/srw/Dropbox/Papers/PQS/provenance/h2plus_R2_tail2_clean_2026-07-26/`
  - `manager_pass`: `433`
  - `git_commit`: `1e373d0c1`
  - `manager_pass`: `434`
  - `external_path`: `/Users/srw/Dropbox/Papers/PQS/validation/h2_supplemented_one_body_preflight_contract_2026-07-26.md`
  - `manager_pass`: `435`
  - `git_commit`: `119bc17c7`
  - `manager_pass`: `436`
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/Papers/PQS/validation/h2_supplemented_one_body_2026-07-26.md`
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/Papers/PQS/validation/h2plus_h2_acceptance_contract.md`
  - `manager_pass`: `437`
  - `git_commit`: `7d2b6dc61`
  - `manager_pass`: `438`
- **Dependencies:** `HP-CGRB-NN-FN-01`, `HP-COMP-SUPPWL-FN-01`, `HP-FN-04`, `HP-FN-05`, `HP-PQS-COULOMB-ACCURACY-FN-01`, `HP-R3-FN-01`, `HP-R3-FN-02`, `HP-R3-FN-03`, `HP-R3U-ZDI-FN-01`, `HP-RG-CUTOFF-FN-02`
- **Scope:** Maintain only the existing private \`bin/pqs\_paper\_h2\_driver.jl\` H2+ controls, neutral-H2 supplemented one-body preflight, and frozen five-row fixed-state density-density measurement at R=2, padding=10.0, tail\_spacing=2.8, and method=:both. Preserve every accepted one-body field, the common normalized full-parent H2+ target, exact parent-metric and native terminal/residual projections, blockwise terminal-coordinate identity, sign-canonical state fingerprints, matrix-free high135 parent ordinary IDA, bare terminal/site IDA, augmented terminal-IDA plus residual-MWG, fresh G-G parity, thirteen fixed interaction fields, and readable category/fingerprint/charge/symmetry/timing diagnostics. Keep the driver \<=525 lines with zero source/test/tool/helper additions. Independent same-density direct-Coulomb validation remains external. No new input, output file, artifact, payload, state carrier, source helper, test, companion driver, RHF/SCF, orbital relaxation, numerical integration or density oracle in the repo, continuum-exact or transition-exchange claim, scan, other H2 endpoint, supplement/cutoff change, enrichment, PRF, screening, EGOI, public API, canonical-driver change, or framework is maintenance-authorized.

### HP-PQS-PARENT-GDIRECT-FN-01 - parent-backed Gaussian direct interaction

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [parent\_residual\_functions.md](parent_residual_functions.md); heading `Parent-Backed Gaussian Direct Interaction`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_ida.jl`
  - `source` / `existing`: `src/ordinary_coulomb.jl`
- **Evidence:**
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/cr2/reports/cr2_prf_final_g_interaction_audit_2026-07-14.md`
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/cr2/reports/cr2_prf_gaussian_distance_interaction_audit_2026-07-14.md`
  - `git_commit`: `5b46ae073`
  - `manager_pass`: `415`
  - `manager_pass`: `416`
- **Dependencies:** `HP-OBJ-02`, `HP-PQS-COULOMB-ACCURACY-FN-01`, `HP-PQS-PRF-FN-01`
- **Scope:** Maintain only the implemented internal onsite-calibrated Gaussian direct resource for explicit current-parent coefficients, including mapped centers, same-expansion positive parent-IDA onsite values, unrenormalized parent charges, tiled PRF-PRF/PRF-G direct blocks, and the bounded parent-IDA direct comparator. Existing G-G Vee, exchange, GTO residuals, public inputs, artifacts, and solvers remain unchanged.

### HP-PQS-PARENT-GDIRECT-TEST-01 - parent-backed Gaussian direct validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [parent\_residual\_functions.md](parent_residual_functions.md); heading `Parent-Backed Gaussian Direct Interaction`
- **Owned paths:**
  - `test` / `existing`: `test/core/runtests.jl`
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:**
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/cr2/reports/cr2_prf_final_g_interaction_audit_2026-07-14.md`
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/cr2/reports/cr2_prf_gaussian_distance_interaction_audit_2026-07-14.md`
  - `git_commit`: `5b46ae073`
  - `manager_pass`: `415`
  - `manager_pass`: `416`
- **Dependencies:** `HP-PQS-PARENT-GDIRECT-FN-01`
- **Scope:** Maintain only the committed analytic onsite, symmetry, finiteness, positive-semidefinite, far-tail, expansion-identity, charge, tiled-parity, bounded parent-IDA, and unchanged-G-G checks. No transition-exchange, PRF-GTO-residual, Cr2, HF, endpoint, or complete-Hamiltonian assertion.

### HP-PQS-PRF-CONSUMER-FN-01 - private PRF consumer diagnostics

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [parent\_backed\_injected\_composition.md](parent_backed_injected_composition.md); heading `Consumer API Contract`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
- **Evidence:**
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/hooke/reports/hooke_be_ns5_prf_completion_status_2026-08-06.md`
  - `manager_pass`: `449`
  - `git_commit`: `25d670e15`
  - `manager_pass`: `450`
  - `git_commit`: `be59362b8`
  - `manager_pass`: `451`
  - `manager_pass`: `475`
  - `git_commit`: `820c2fc8bc9fa25ca699994e3eff0805b5b5eb98`
  - `manager_pass`: `476`
- **Dependencies:** `HP-PQS-PRF-FN-01`, `HP-PQS-PRF-INJECT-COMP-FN-01`, `HP-PQS-PRF-INJECT-INTERACT-FN-01`, `HP-R3U-FN-01`
- **Scope:** Maintain only the six unexported PRF diagnostic/provenance definitions in src/cartesian\_base\_hamiltonian.jl: CartesianParentResidualRegion, CartesianParentBackedHamiltonianResult, cartesian\_parent\_residual\_regions, cartesian\_parent\_residual\_block, cartesian\_parent\_backed\_composition, and cartesian\_parent\_backed\_hamiltonian. Preserve exact descriptor-to-PRF source binding, delegation to the independently owned PRF/injection/operator/interaction numerics, and unchanged historical Cr2 reproduction capability. Do not re-export them or add an alias, shim, deprecation, helper, replacement, type, field, signature, metadata, default, artifact, driver, solver, numerical change, broader PRF cleanup, or public compatibility promise.

### HP-PQS-PRF-CONSUMER-TEST-01 - qualify internal PRF consumer validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [parent\_backed\_injected\_composition.md](parent_backed_injected_composition.md); heading `Consumer API Contract`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:**
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/hooke/reports/hooke_be_ns5_prf_completion_status_2026-08-06.md`
  - `manager_pass`: `449`
  - `git_commit`: `25d670e15`
  - `manager_pass`: `450`
  - `git_commit`: `be59362b8`
  - `manager_pass`: `451`
  - `manager_pass`: `475`
  - `git_commit`: `820c2fc8bc9fa25ca699994e3eff0805b5b5eb98`
  - `manager_pass`: `476`
- **Dependencies:** `HP-PQS-PRF-CONSUMER-FN-01`, `HP-PQS-PRF-INJECT-COMP-TEST-01`, `HP-PQS-PRF-INJECT-INTERACT-TEST-01`
- **Scope:** Maintain only the existing qualified internal PRF validation in test/nested/cartesian\_r3a\_h2\_augmented\_one\_body\_runtests.jl: semantic-region uniqueness, stale-basis and cross-region rejection, additive/fixed-span composition parity, exact final-Hamiltonian parity, and omitted-path and malformed-input gates. Keep every PRF call and type assertion module-qualified. Add or remove no assertion, fixture, owner, file, probe, endpoint policy, public API promise, or numerical behavior.

### HP-PQS-PRF-FN-01 - parent residual function mechanics

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [parent\_residual\_functions.md](parent_residual_functions.md); heading `Parent Residual Function Construction`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`
- **Evidence:**
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/cr2/reports/cr2_fixed_parent_shell_q5_q9_completion_ladder_2026-07-14.md`
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/cr2/reports/cr2_shell_local_parent_completion_occupancy_audit_2026-07-14.md`
  - `git_commit`: `5b46ae073`
  - `manager_pass`: `415`
  - `manager_pass`: `416`
- **Dependencies:** `HP-FN-03`, `HP-OBJ-02`
- **Scope:** Maintain only the implemented internal support-local PRF mechanics for consumer-selected parent target columns, including strict support/rank/metric failures and exact native PRF-G/PRF-PRF one-body blocks. No selection, RDM, shell-q, artifact, Hamiltonian, or solver policy.

### HP-PQS-PRF-INJECT-COMP-FN-01 - parent-backed injected basis composition

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [parent\_backed\_injected\_composition.md](parent_backed_injected_composition.md); heading `Fixed-Span Injection`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  - `source` / `existing`: `src/cartesian_protected_ladder_bundle.jl`
- **Evidence:**
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/cr2/reports/cr2_q7_injected_numcomplete_screened_uhf_2026-07-15.md`
  - `git_commit`: `cdd2c27af`
  - `manager_pass`: `418`
  - `manager_pass`: `419`
- **Dependencies:** `HP-PQS-PRF-FN-01`, `HP-RG-NUMCOMP-FN-01`
- **Scope:** Maintain only the implemented private fixed-span angular-style composition: project every old localized terminal seed before symmetric Lowdin, retain the exact parent-backed complement, preserve native B = \[Ginj,Rnew,RGexternal\] span/dimension, residualize the supplement at the numerical-complete 1e-10 policy, and rebuild exact native one-body operators. No consumer selection, RDM, shell, q, cutoff, count, orientation, public, artifact, solver, or Cr2 policy.

### HP-PQS-PRF-INJECT-COMP-TEST-01 - parent-backed injected basis validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [parent\_backed\_injected\_composition.md](parent_backed_injected_composition.md); heading `Fixed-Span Injection`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:**
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/cr2/reports/cr2_q7_injected_numcomplete_screened_uhf_2026-07-15.md`
  - `git_commit`: `cdd2c27af`
  - `manager_pass`: `418`
  - `manager_pass`: `419`
- **Dependencies:** `HP-PQS-PRF-INJECT-COMP-FN-01`
- **Scope:** Maintain only the committed bounded H2 construction, complete old-seed projection, malformed/rank/support/weight/span failure, stale augmentation, exact one-body oracle, numerical-complete residual, native-order, and omitted-path parity checks. Padded Be2 packet/correction validation is external manager-reviewed evidence, not a committed fixture or owned test surface. No Cr2, HF, endpoint, selection, or localization assertion.

### HP-PQS-PRF-INJECT-INTERACT-FN-01 - parent-backed injected mixed interaction

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [parent\_backed\_injected\_composition.md](parent_backed_injected_composition.md); heading `Interaction Block Contract`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_ida.jl`
  - `source` / `existing`: `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  - `source` / `existing`: `src/cartesian_protected_ladder_bundle.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/mwg_interaction.jl`
- **Evidence:**
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/cr2/reports/cr2_q7_injected_numcomplete_screened_uhf_2026-07-15.md`
  - `git_commit`: `006432e9d`
  - `manager_pass`: `418`
  - `manager_pass`: `420`
- **Dependencies:** `HP-PQS-PARENT-GDIRECT-FN-01`, `HP-PQS-PRF-INJECT-COMP-FN-01`, `HP-PQS-SCREEN-HARTREE-CORR-FN-01`
- **Scope:** Maintain only the implemented native \[Ginj,Rnew,RGexternal\] assembly: freshly rebuilt terminal IDA, onsite-calibrated parent-Gaussian direct blocks involving Rnew, established MWG blocks involving RGexternal, fitted atomic fields represented in the complete basis, and unchanged screened-Hartree delegation after H and V are fixed. No C'VC, old-interaction reuse, exact PRF-GTO, transition exchange, public, artifact, sidecar, solver, HF, or Cr2 policy.

### HP-PQS-PRF-INJECT-INTERACT-TEST-01 - parent-backed injected interaction validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [parent\_backed\_injected\_composition.md](parent_backed_injected_composition.md); heading `Interaction Block Contract`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:**
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/cr2/reports/cr2_q7_injected_numcomplete_screened_uhf_2026-07-15.md`
  - `git_commit`: `006432e9d`
  - `manager_pass`: `418`
  - `manager_pass`: `420`
- **Dependencies:** `HP-PQS-PRF-INJECT-INTERACT-FN-01`
- **Scope:** Maintain only the committed bounded H2 freshly rebuilt terminal, separated IDA/Gaussian/MWG block, finite positive moment, stale-input, correction-input, unchanged-path, and parity checks. Padded Be2 packet accounting, correction recomputation, derivative anchor, and unchanged H/V are external manager-reviewed evidence, not a committed fixture or owned test surface. No Cr2, HF, endpoint, exact-exchange, or exact PRF-GTO assertion.

### HP-PQS-PRF-TEST-01 - parent residual function validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [parent\_residual\_functions.md](parent_residual_functions.md); heading `Parent Residual Function Construction`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:**
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/cr2/reports/cr2_fixed_parent_shell_q5_q9_completion_ladder_2026-07-14.md`
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/cr2/reports/cr2_shell_local_parent_completion_occupancy_audit_2026-07-14.md`
  - `git_commit`: `5b46ae073`
  - `manager_pass`: `415`
  - `manager_pass`: `416`
- **Dependencies:** `HP-PQS-PRF-FN-01`
- **Scope:** Maintain only the committed synthetic and bounded H2 checks for consumer-selected targets, strict malformed-input failures, deterministic phases, terminal/PRF metric identities, exact one-body oracle parity, and unchanged G construction. No Cr2, HF, or selection-policy assertion.

### HP-PQS-PUBLIC-COMPAT-FN-01 - PQS v0.2 dependency compatibility declaration

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `Dependency Compatibility Lifecycle`
- **Owned paths:**
  - `source` / `existing`: `Project.toml`
- **Evidence:**
  - `git_commit`: `873b6e27a0f3e049fda70b4ff16dc0682354efce`
  - `git_commit`: `8bb5b88009cad222fcd5fddbfcd8ad65555e811a`
  - `manager_pass`: `485`
  - `manager_pass`: `486`
- **Dependencies:** `HP-PQS-PUBLIC-DOC-01`, `HP-PQS-PUBLIC-MATCHED-FN-01`, `HP-PQS-PUBLIC-SCREEN-FN-01`
- **Scope:** Maintain exactly JLD2 = 0.6.4, LinearAlgebra = 1.10, SHA = 0.7, SparseArrays = 1.10, SpecialFunctions = 2.8, and TOML = 1 in the root Project.toml \[compat\] table while retaining julia = 1.10. Package version is owned separately by an explicit candidate/release record. Any added key, widened range, lowered floor, tag, release, registration, or root-manifest change requires separate authority. No src code, export, test, example, docs-content, workflow, or manifest-policy change.

### HP-PQS-PUBLIC-COMPAT-TEST-01 - PQS v0.2 dependency compatibility validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `Dependency Compatibility Lifecycle`
- **Owned paths:** none
- **Evidence:**
  - `git_commit`: `873b6e27a0f3e049fda70b4ff16dc0682354efce`
  - `git_commit`: `8bb5b88009cad222fcd5fddbfcd8ad65555e811a`
  - `manager_pass`: `485`
  - `manager_pass`: `486`
- **Dependencies:** `HP-PQS-PUBLIC-COMPAT-FN-01`, `HP-PQS-PUBLIC-MATCHED-TEST-01`, `HP-PQS-PUBLIC-SCREEN-TEST-01`
- **Scope:** Completed validation evidence only; this record grants no test or workflow edit. The accepted declaration passed exact six-bound Project.toml parsing, fresh isolated Julia 1.12.6 resolution with archived manifest and direct-version proof, package load, unchanged bounded public examples and focused release gate, docs resolution/build, authority and docs checks, bounded numerical groups, successful Julia 1.10 CI, clean diff, and confirmation that no root Manifest.toml is tracked. Stop on any future incompatible floor rather than widening or lowering a bound.

### HP-PQS-PUBLIC-DOC-01 - PQS v0.2 public documentation

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `docs`
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `User Documentation And Layering`
- **Owned paths:**
  - `docs` / `existing`: `docs/src/developer/designs/cartesian_hamiltonian_producer/README.md`
  - `docs` / `existing`: `docs/src/developer/designs/cartesian_hamiltonian_producer/pqs_public_surface.md`
  - `docs` / `existing`: `docs/src/examples/index.md`
  - `docs` / `existing`: `docs/src/howto/example_guide.md`
  - `docs` / `existing`: `docs/src/manual/index.md`
  - `docs` / `existing`: `docs/src/manual/projected_q_shells.md`
  - `docs` / `existing`: `docs/src/manual/reference_density_hartree_screening.md`
- **Evidence:**
  - `external_path`: `/Users/srw/Dropbox/Papers/PQS/notes/assignments/repo_design_manager_pqs_v0p2_public_surface_2026-08-17.md`
  - `git_commit`: `058ee54f45c759949f70b54a699ccc318476f8ac`
  - `manager_pass`: `481`
  - `manager_pass`: `482`
- **Dependencies:** `HP-PQS-PUBLIC-MATCHED-FN-01`, `HP-PQS-PUBLIC-SCREEN-FN-01`
- **Scope:** Maintain only the implemented v0.2 reader contract separating stable PQS/WL construction, supplied-field screened-Hartree assembly, and provenance-only Ximg/XHF. Keep the two bounded examples and rendered navigation aligned with the exported API without duplicating subsystem algorithms, exposing private packets/stages, or creating a documentation framework.

### HP-PQS-PUBLIC-DOC-PARITY-FN-01 - PQS v0.2 reader-document parity repair

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `docs`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `Reader-Facing Documentation Parity`
- **Owned paths:**
  - `docs` / `existing`: `docs/example_guide.md`
  - `docs` / `existing`: `docs/src/algorithms/cartesian_nested_diatomic_box_policy.md`
  - `docs` / `existing`: `docs/src/algorithms/cartesian_nested_diatomic_coordinate_distortion.md`
  - `docs` / `existing`: `docs/src/howto/example_guide.md`
  - `docs` / `existing`: `docs/src/reference/export.md`
  - `source` / `existing`: `src/cartesian_reference_density/screened_hartree_correction.jl`
  - `source` / `existing`: `src/pqs_matched_h2plus.jl`
- **Evidence:**
  - `git_commit`: `c78defcc9b299ee5f32cf42910812f5581657d93`
  - `manager_pass`: `487`
  - `manager_pass`: `488`
  - `manager_pass`: `489`
- **Dependencies:** `HP-PQS-PUBLIC-DOC-01`, `HP-PQS-PUBLIC-MATCHED-FN-01`, `HP-PQS-PUBLIC-SCREEN-FN-01`, `HP-QW-NESTED-DIAT-FN-01`
- **Scope:** Maintain only the implemented v0.2 reader parity: concise public docstrings on the exact ten accepted bindings added by c78defcc9, the unchanged existing pqs\_h2plus\_comparison docstring, all eleven existing root exports in the curated @docs reference, path-free historical provenance in the two owned nested-diatomic pages, and example 28's internal/non-contractual classification. Preserve the typed field semantics, common-basis and signed-consistency conventions, and accessor-only ScreenedHartreeCorrection contract. No executable source, definition, signature, field, export, dispatch, numerical, example-code, dependency, workflow, version, tag, release, citation, paper claim, API/example promotion, or broad rewrite.

### HP-PQS-PUBLIC-DOC-PARITY-TEST-01 - PQS v0.2 reader-document parity validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `Reader-Facing Documentation Parity`
- **Owned paths:**
  - `test` / `existing`: `test/docs/runtests.jl`
- **Evidence:**
  - `git_commit`: `c78defcc9b299ee5f32cf42910812f5581657d93`
  - `manager_pass`: `487`
  - `manager_pass`: `488`
  - `manager_pass`: `489`
- **Dependencies:** `HP-PQS-PUBLIC-DOC-PARITY-FN-01`, `HP-PQS-READER-TEST-01`
- **Scope:** Maintain only the fixed eleven-name root-export/curated-reference parity check, absence of personal absolute or Dropbox paths in the two owned rendered algorithm pages, and explicit internal/non-contractual example-28 classification in both reader guides. Preserve Documenter as the executable docstring-resolution gate. Add no test file, parser, generalized export inventory, example execution, workflow, numerical gate, or broad documentation scan.

### HP-PQS-PUBLIC-MATCHED-FN-01 - matched H2+ public comparison

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `Matched H2+ Public Comparison`
- **Owned paths:**
  - `source` / `existing`: `src/pqs_matched_h2plus.jl`
- **Evidence:**
  - `external_path`: `/Users/srw/Dropbox/Papers/PQS/notes/assignments/repo_design_manager_pqs_v0p2_public_surface_2026-08-17.md`
  - `external_path`: `/Users/srw/Dropbox/Papers/PQS/reproduction/pqs_release_candidate_2026-08-17/data/table1_h2plus/table1_h2plus_matched_composite.tsv`
  - `external_path`: `/Users/srw/dmrgtmp/pqs_perf_audit_f9ca5a7b_20260825/report.md`
  - `external_path`: `/Users/srw/dmrgtmp/pass530_parent_workspace_20260825/report.md`
  - `git_commit`: `058ee54f45c759949f70b54a699ccc318476f8ac`
  - `git_commit`: `fbef95e60ca3aadafe2082b4a63f8522a724e7be`
  - `manager_pass`: `481`
  - `manager_pass`: `482`
  - `manager_pass`: `530`
  - `manager_pass`: `531`
- **Dependencies:** `HP-PQS-ASPECTSHELL-FN-01`, `HP-PQS-COULOMB-ACCURACY-FN-01`, `HP-R1-FN-01`
- **Scope:** Maintain only the accepted call-local workspace implementation in \`src/pqs\_matched\_h2plus.jl\`: each \`\_pqs\_h2plus\_parent\_solution\` invocation owns disjoint output and scratch storage used by the mutating three-axis product action and \`raw\_nuclear\` accumulation. Preserve exact deterministic action parity, axis/product/Gaussian-term arithmetic order, no parent-sized per-action allocation, independent concurrent-call ownership, and all dimensions/fingerprints/topology/energies/captures/residuals/symmetry/warnings and release tolerances. Do not restore the allocating helpers or add a mutable global, task-local store, pool, cache, lock, hidden persistent state, escaping scratch view, public workspace/API/type, second source owner, test edit, dependency, metadata, workflow, example, fixture, or release change.

### HP-PQS-PUBLIC-MATCHED-TEST-01 - matched H2+ public release validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `Validation Contract`
- **Owned paths:**
  - `test` / `existing`: `test/pqs_h2plus_table1_release_runtests.jl`
  - `test` / `existing`: `test/runtests.jl`
- **Evidence:**
  - `external_path`: `/Users/srw/Dropbox/Papers/PQS/reproduction/pqs_release_candidate_2026-08-17/data/table1_h2plus/table1_h2plus_matched_composite.tsv`
  - `git_commit`: `058ee54f45c759949f70b54a699ccc318476f8ac`
  - `git_commit`: `72c46f9ea0dd6b2da7a6a302d34ea1c501d18647`
  - `git_commit`: `b0dbd9ea37317590334a24883ef0667bdb0195a5`
  - `manager_pass`: `481`
  - `manager_pass`: `482`
  - `manager_pass`: `483`
  - `manager_pass`: `484`
  - `manager_pass`: `536`
  - `manager_pass`: `537`
- **Dependencies:** `HP-PQS-PUBLIC-MATCHED-FN-01`
- **Scope:** Maintain the accepted single-execution matched-H2+ release arrangement from commit b0dbd9ea37317590334a24883ef0667bdb0195a5. Example 41 remains independently executable, writes the same three-row, eight-column TSV and concise summary, and returns its already-computed PQSH2PlusComparison. The release owner includes that example once and applies all existing 18 assertions to the returned comparison; test/runtests.jl must not restore the duplicate Example 41 subprocess. Preserve result-type checks, exact common topology and 12789/1285/1285 dimensions, capture atol 2e-8, energy/error atol 1e-7 hartree, accounting checks, public API coverage, nonfinite-reference rejection, required-check identity, and documentation links. The successful pqs\_release path constructs the complete comparison exactly once. Accepted fresh cost is 54.91 seconds and 9.716 GiB versus the former 110.56 seconds and 19.342 GiB. Preserve the +4/-5 three-file implementation without adding a cache, fixture, artifact, helper, framework, file, API, numerical change, workflow, docs, or release action. The clean candidate replay rules remain unchanged.

### HP-PQS-PUBLIC-RC1-FN-01 - v0.2.0-rc1 candidate preparation

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `docs`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `v0.2.0-rc1 Candidate Preparation`
- **Owned paths:**
  - `docs` / `existing`: `CHANGELOG.md`
  - `source` / `existing`: `Project.toml`
- **Evidence:**
  - `git_commit`: `41fa897ae919510b179a425027a8ce2d4a2167b3`
  - `manager_pass`: `492`
  - `manager_pass`: `493`
- **Dependencies:** `HP-PQS-DOCS-TAGDEPLOY-FN-01`, `HP-PQS-PRF-CONSUMER-FN-01`, `HP-PQS-PUBLIC-COMPAT-FN-01`, `HP-PQS-PUBLIC-DOC-PARITY-FN-01`, `HP-PQS-PUBLIC-MATCHED-FN-01`, `HP-PQS-PUBLIC-SCREEN-FN-01`, `HP-PQS-READER-DOC-01`, `HP-QW-NESTED-DIAT-FN-01`
- **Scope:** Maintain the prepared v0.2.0-rc1 candidate identity implemented by commit 41fa897ae919510b179a425027a8ce2d4a2167b3: root Project.toml version 0.2.0-rc1 and one concise root CHANGELOG.md grouped as Added, Changed, Fixed, Public-surface reduction, and Scope. Preserve traceability to accepted public APIs, examples 39-41, dependency bounds, documentation deployment, ordinary-QW fixes, and six-name PRF de-promotion without pass history, private numerical evidence, manuscript claims, dates, DOI/citation invention, or internal deletion inventories. This record itself authorizes no tag or release; the completed annotated tag is owned by HP-PQS-PUBLIC-RC1-TAG-FN-01 and the exact GitHub prerelease operation is separately owned by HP-PQS-PUBLIC-RC1-RELEASE-FN-01. No registration, CITATION.cff, source/export/test/example/dependency/workflow/API/numerical/scientific-doc change, stable link, badge, homepage change, tracked root manifest, or paper-workspace import.

### HP-PQS-PUBLIC-RC1-RELEASE-FN-01 - v0.2.0-rc1 GitHub prerelease publication

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `v0.2.0-rc1 GitHub Prerelease Lifecycle`
- **Owned paths:** none
- **Evidence:**
  - `external_path`: `/Users/srw/dmrgtmp/release_replays/gaussletbases_v0p2_rc1_release_2026-08-19`
  - `git_commit`: `1546c18d3058cce2b5051b50788cda3c12585e51`
  - `git_commit`: `d653b267625e19059c72ae6925040b91b77f85fa`
  - `git_commit`: `e63cdf6359293e7c274a6556357dc831b2f7eb02`
  - `manager_pass`: `497`
  - `manager_pass`: `498`
- **Dependencies:** `HP-PQS-DOCS-TAGDEPLOY-FN-01`, `HP-PQS-PUBLIC-RC1-FN-01`
- **Scope:** Completed GitHub prerelease publication only; this record grants no further release, tag, file, workflow, repository-metadata, or asset mutation. Release 373460389 is visible at the exact v0.2.0-rc1 tag with title GaussletBases v0.2.0-rc1, prerelease=true, draft=false, no uploaded assets, no latest final release, and the exact 1455-byte package-centered body separating Projected q-shells (PQS) from reference-density Hartree screening. Annotated tag object a4284f0bf448fb9d717de26ccbe1e9fc16db5ed2 still peels to frozen commit 1546c18d3058cce2b5051b50788cda3c12585e51. Automatic source archives and a fresh isolated Julia 1.12.6 installation validate tree 1b53a9eb51d11cfc31b8b0356349c62f0de8915f; RC1 and dev documentation remain live, versions.js lists RC1 and dev, and /stable/ remains absent. Preserve the published release without editing, deleting, recreating, retargeting, adding assets, or changing its narrative. Final v0.2.0, General registration, CITATION.cff, paper titles/status/citations, stable links, repository metadata, and source/API/test/example/dependency/workflow/numerical changes require separate authority.

### HP-PQS-PUBLIC-RC1-RELEASE-TEST-01 - v0.2.0-rc1 GitHub prerelease validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `v0.2.0-rc1 GitHub Prerelease Lifecycle`
- **Owned paths:** none
- **Evidence:**
  - `external_path`: `/Users/srw/dmrgtmp/release_replays/gaussletbases_v0p2_rc1_release_2026-08-19`
  - `git_commit`: `1546c18d3058cce2b5051b50788cda3c12585e51`
  - `git_commit`: `d653b267625e19059c72ae6925040b91b77f85fa`
  - `git_commit`: `e63cdf6359293e7c274a6556357dc831b2f7eb02`
  - `manager_pass`: `497`
  - `manager_pass`: `498`
- **Dependencies:** `HP-PQS-PUBLIC-RC1-RELEASE-FN-01`
- **Scope:** Completed post-publication evidence only; this record grants no file, release, tag, workflow, repository-metadata, or asset mutation. GitHub API and rendered-page checks confirm release 373460389 has exact tag v0.2.0-rc1, title and 1455-byte canonical body, draft=false, prerelease=true, no uploaded assets, and no latest final release. Both automatic source archives validated; a fresh isolated Julia 1.12.6 installation loaded GaussletBases v0.2.0-rc1 with manifest tree 1b53a9eb51d11cfc31b8b0356349c62f0de8915f. RC1 and dev documentation retain exact canonical URLs, versions.js lists RC1 and dev, and /stable/ remains absent. Preserve the accepted release and evidence without edit, deletion, recreation, retargeting, asset upload, or silent retry.

### HP-PQS-PUBLIC-RC1-TAG-FN-01 - v0.2.0-rc1 annotated tag publication

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `v0.2.0-rc1 Annotated Tag Lifecycle`
- **Owned paths:** none
- **Evidence:**
  - `git_commit`: `1546c18d3058cce2b5051b50788cda3c12585e51`
  - `manager_pass`: `494`
  - `manager_pass`: `495`
- **Dependencies:** `HP-PQS-DOCS-TAGDEPLOY-FN-01`, `HP-PQS-PUBLIC-RC1-FN-01`
- **Scope:** Completed immutable annotated-tag publication only; this record grants no further tag, file, workflow, or release mutation. Local and remote tag object a4284f0bf448fb9d717de26ccbe1e9fc16db5ed2 peels to frozen target 1546c18d3058cce2b5051b50788cda3c12585e51 with tree 1b53a9eb51d11cfc31b8b0356349c62f0de8915f and target git-archive SHA-256 2a0b6938d3b341900d73668e7f0644c34b8a851e1b823356c53c2866fd19522a. Tag-triggered Docs run 32295705338 published the canonical RC1 folder without creating stable or changing dev. Never move, replace, delete, or recreate the tag. Selector-index repair is separately maintained by HP-PQS-DOCS-TAGDEPLOY-FN-01. This completed record grants no GitHub release; the exact RC1 prerelease operation is separately owned by HP-PQS-PUBLIC-RC1-RELEASE-FN-01. No registration, CITATION.cff, final v0.2.0, homepage/stable-link edit, or tracked source/API/dependency/workflow/numerical/manuscript change.

### HP-PQS-PUBLIC-RC1-TAG-TEST-01 - v0.2.0-rc1 annotated tag validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `v0.2.0-rc1 Annotated Tag Lifecycle`
- **Owned paths:** none
- **Evidence:**
  - `git_commit`: `1546c18d3058cce2b5051b50788cda3c12585e51`
  - `git_commit`: `31caa87d3b83599de7f7295678ee599209113552`
  - `manager_pass`: `494`
  - `manager_pass`: `495`
  - `manager_pass`: `496`
- **Dependencies:** `HP-PQS-PUBLIC-RC1-TAG-FN-01`
- **Scope:** Completed tag-acceptance evidence only; this record grants no file, workflow, tag-repair, or release edit. Annotated local and remote tag object a4284f0bf448fb9d717de26ccbe1e9fc16db5ed2 peels to frozen target 1546c18d3058cce2b5051b50788cda3c12585e51; tag-triggered Docs run 32295705338 and repair main-deployment Docs run 32302304167 passed; /v0.2.0-rc1/ and /dev/ retain exact canonical URLs; versions.js lists RC1 and dev; /stable/ remains absent. Documenter's internal DOCUMENTER\_STABLE fallback naming RC1 creates no stable selector entry, symlink, path, or final-release status. Preserve the immutable tag; do not retarget, replace, delete, or recreate it. This completed record grants no release publication; the exact RC1 GitHub prerelease is separately owned by HP-PQS-PUBLIC-RC1-RELEASE-FN-01.

### HP-PQS-PUBLIC-RC1-TEST-01 - v0.2.0-rc1 candidate-preparation validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `v0.2.0-rc1 Candidate Preparation`
- **Owned paths:** none
- **Evidence:**
  - `git_commit`: `41fa897ae919510b179a425027a8ce2d4a2167b3`
  - `manager_pass`: `492`
  - `manager_pass`: `493`
- **Dependencies:** `HP-PQS-PUBLIC-RC1-FN-01`
- **Scope:** Completed acceptance evidence only; this record grants no test, workflow, source, docs, tag, or release edit. Commit 41fa897ae919510b179a425027a8ce2d4a2167b3 passed exact Project/changelog review, fresh Julia 1.12.6 GitHub installation, package load, public examples 01/39/40/41, focused H2+ 18/18, Julia 1.10 CI, authority/docs/manager-log checks, exact prerelease documentation classification without stable advancement, and a clean git archive excluding a root manifest and untracked handoffs. The completed annotated tag is owned by HP-PQS-PUBLIC-RC1-TAG-FN-01; the exact GitHub prerelease operation is separately owned by HP-PQS-PUBLIC-RC1-RELEASE-FN-01.

### HP-PQS-PUBLIC-RC2-FN-01 - v0.2.0-rc2 reader front door and candidate preparation

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `docs`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `v0.2.0-rc2 Reader Front Door And Candidate Preparation`
- **Owned paths:**
  - `docs` / `existing`: `CHANGELOG.md`
  - `docs` / `existing`: `README.md`
  - `source` / `existing`: `Project.toml`
- **Evidence:**
  - `git_commit`: `2b3c23970144aa030ae52b875a5cf01b32886b6e`
  - `git_commit`: `de8991a75d5565e65dcabb9d80c5626a5b86905d`
  - `manager_pass`: `519`
  - `manager_pass`: `520`
- **Dependencies:** `HP-PQS-DOCS-TAGDEPLOY-FN-01`, `HP-PQS-PUBLIC-RC1-FN-01`, `HP-PUBLIC-EXPORT-INTEGRITY-FN-01`, `HP-PUBLIC-PAPER-CI-FN-01`, `HP-REP-PQS-RG-WORKING-FN-01`, `HP-REP-XGTO-CLOSESTDET-FN-01`, `HP-REP-XGTO-INTERCHANGE-FN-01`, `HP-REP-XGTO-PYSCF-EXPORT-FN-01`, `HP-REP-XGTO-READER-DOC-01`
- **Scope:** Maintain the exact v0.2.0-rc2 candidate identity implemented by commit 2b3c23970144aa030ae52b875a5cf01b32886b6e: root Project.toml version 0.2.0-rc2, the concise post-RC1 CHANGELOG section above a byte-unchanged RC1 section, and the three distinct root-README capability links with the radial beginner route retained and the external exporter described as checkpoint-only. PySCF and NumPy remain optional external-command dependencies. Exact RC2/RC1/dev selector maintenance belongs to HP-PQS-DOCS-TAGDEPLOY-FN-01. This record authorizes no tag or release; the completed immutable annotated tag is recorded by HP-PQS-PUBLIC-RC2-TAG-FN-01 and the exact GitHub prerelease operation is separately owned by HP-PQS-PUBLIC-RC2-RELEASE-FN-01. No production source/API/export/numerical/dependency/example/workflow/fixture-format/manifest-policy/test-owner change, tutorial duplication, manuscript or benchmark claim, paper status, DOI, stable alias, registration, citation, or final-v0.2 action is authorized.

### HP-PQS-PUBLIC-RC2-RELEASE-FN-01 - v0.2.0-rc2 GitHub prerelease publication

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `v0.2.0-rc2 GitHub Prerelease Lifecycle`
- **Owned paths:** none
- **Evidence:**
  - `git_commit`: `2b3c23970144aa030ae52b875a5cf01b32886b6e`
  - `git_commit`: `69e11242809be85f4d892d082d25ac66467ae373`
  - `git_commit`: `b341be8814e654eb6039e18d1033c3a71936019b`
  - `manager_pass`: `523`
  - `manager_pass`: `524`
- **Dependencies:** `HP-PQS-DOCS-TAGDEPLOY-FN-01`, `HP-PQS-PUBLIC-RC2-FN-01`
- **Scope:** Completed GitHub prerelease publication only; this record grants no further release, tag, file, workflow, repository-metadata, or asset mutation. Release 376503169 is visible at the exact v0.2.0-rc2 tag with title GaussletBases v0.2.0-rc2, prerelease=true, draft=false, no uploaded assets, no latest final release, and the exact 2200-byte ASCII body with SHA-256 a2cbaaa2a349857e897d6d58fb728c6ccd9d731c7371bb39997bfc0360f3653a. Annotated tag object 7c8a21b998a838d245e0b5a7f4915910e2a091bc still peels to frozen commit 2b3c23970144aa030ae52b875a5cf01b32886b6e and tree 7a4b51aec25f62436620f4ff938262d0f6b2fd62. The 2564018-byte tarball and 2957570-byte zipball reconstruct that tree; a fresh isolated Julia 1.12.6 installation loaded GaussletBases 0.2.0-rc2 with the same manifest tree. RC2, RC1, and dev documentation remain live, and /stable/ remains absent. Preserve the published release without editing, deleting, recreating, retargeting, adding assets, or changing its narrative. Final v0.2.0, registration, citation metadata, stable documentation, repository metadata, and source/API/test/example/dependency/workflow/numerical changes require separate authority.

### HP-PQS-PUBLIC-RC2-RELEASE-TEST-01 - v0.2.0-rc2 GitHub prerelease validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `v0.2.0-rc2 GitHub Prerelease Lifecycle`
- **Owned paths:** none
- **Evidence:**
  - `git_commit`: `2b3c23970144aa030ae52b875a5cf01b32886b6e`
  - `git_commit`: `69e11242809be85f4d892d082d25ac66467ae373`
  - `git_commit`: `b341be8814e654eb6039e18d1033c3a71936019b`
  - `manager_pass`: `523`
  - `manager_pass`: `524`
- **Dependencies:** `HP-PQS-PUBLIC-RC2-RELEASE-FN-01`
- **Scope:** Completed post-publication evidence only; this record grants no file, release, tag, workflow, repository-metadata, or asset mutation. GitHub API and rendered-page checks confirm release 376503169 has exact tag v0.2.0-rc2, title GaussletBases v0.2.0-rc2, 2200-byte canonical body with SHA-256 a2cbaaa2a349857e897d6d58fb728c6ccd9d731c7371bb39997bfc0360f3653a, draft=false, prerelease=true, no uploaded assets, and no latest final release. The automatic tarball has 2564018 bytes and SHA-256 5e92245b92350865415facf1350ba186b58052f0befb6380a5397d1eadaef445; the zipball has 2957570 bytes and SHA-256 84a5619ab679650ceb6e7c0ab0e728d078cb0785d1b26d4f83f6e81471b24d30; both reconstruct tree 7a4b51aec25f62436620f4ff938262d0f6b2fd62. A fresh isolated Julia 1.12.6 tag installation loaded GaussletBases 0.2.0-rc2 with that manifest tree. RC2, RC1, and dev documentation remain live, versions.js lists all three, and /stable/ remains absent. Preserve the accepted release and immutable tag without edit, deletion, recreation, retargeting, asset upload, or silent retry.

### HP-PQS-PUBLIC-RC2-TAG-FN-01 - v0.2.0-rc2 annotated tag publication

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `v0.2.0-rc2 Annotated Tag Lifecycle`
- **Owned paths:** none
- **Evidence:**
  - `git_commit`: `2b3c23970144aa030ae52b875a5cf01b32886b6e`
  - `manager_pass`: `521`
  - `manager_pass`: `522`
- **Dependencies:** `HP-PQS-DOCS-TAGDEPLOY-FN-01`, `HP-PQS-PUBLIC-RC2-FN-01`, `HP-PUBLIC-PAPER-CI-FN-01`
- **Scope:** Completed immutable annotated-tag publication only; this record grants no further tag, file, workflow, or release mutation. Local and remote tag object 7c8a21b998a838d245e0b5a7f4915910e2a091bc has message GaussletBases v0.2.0-rc2 and peels to frozen target 2b3c23970144aa030ae52b875a5cf01b32886b6e with tree 7a4b51aec25f62436620f4ff938262d0f6b2fd62. The accepted archive remains 676 entries and 9994240 bytes with SHA-256 6728b80c1397f13b367c2d898fbdda3176c6cb87c39597817c590a0f41c1e2ac. Tag-triggered Docs run 32798625043 and CI run 32798625045 passed; Pages deployment 32798719038 passed. /v0.2.0-rc2/ and /dev/ are live, versions.js lists RC2, RC1, and dev, and /stable/ remains absent. Never move, replace, delete, or recreate the tag. This completed record itself grants no GitHub release; the exact RC2 prerelease operation is separately owned by HP-PQS-PUBLIC-RC2-RELEASE-FN-01. No asset upload, registration, citation, stable alias, final v0.2.0, or tracked mutation is authorized here.

### HP-PQS-PUBLIC-RC2-TAG-TEST-01 - v0.2.0-rc2 annotated tag validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `v0.2.0-rc2 Annotated Tag Lifecycle`
- **Owned paths:** none
- **Evidence:**
  - `git_commit`: `2b3c23970144aa030ae52b875a5cf01b32886b6e`
  - `manager_pass`: `521`
  - `manager_pass`: `522`
- **Dependencies:** `HP-PQS-DOCS-TAGDEPLOY-TEST-01`, `HP-PQS-PUBLIC-RC2-TAG-FN-01`, `HP-PUBLIC-PAPER-CI-TEST-01`
- **Scope:** Completed tag-acceptance evidence only; this record grants no file, workflow, tag-repair, or release edit. Local and remote annotated tag object 7c8a21b998a838d245e0b5a7f4915910e2a091bc peels to frozen target 2b3c23970144aa030ae52b875a5cf01b32886b6e and tree 7a4b51aec25f62436620f4ff938262d0f6b2fd62. Tag-triggered Docs run 32798625043 passed; CI run 32798625045 passed PQS, Supported floor, and Screening; Pages deployment 32798719038 passed. /v0.2.0-rc2/ has its exact canonical URL, versions.js lists RC2, RC1, and dev, /stable/ remains absent, and /dev/ remains intact. Only the tag ref changed; main and origin/main remained at 707e5bcbb35f375a09c663cfd949f40377d37e19 with both handoffs untouched. Preserve the immutable tag without movement, replacement, deletion, recreation, or silent retry. This completed record grants no release publication; the exact RC2 GitHub prerelease is separately owned by HP-PQS-PUBLIC-RC2-RELEASE-FN-01.

### HP-PQS-PUBLIC-RC2-TEST-01 - v0.2.0-rc2 candidate validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `v0.2.0-rc2 Reader Front Door And Candidate Preparation`
- **Owned paths:** none
- **Evidence:**
  - `git_commit`: `2b3c23970144aa030ae52b875a5cf01b32886b6e`
  - `git_commit`: `de8991a75d5565e65dcabb9d80c5626a5b86905d`
  - `manager_pass`: `519`
  - `manager_pass`: `520`
- **Dependencies:** `HP-PQS-DOCS-TAGDEPLOY-TEST-01`, `HP-PQS-PUBLIC-RC2-FN-01`, `HP-PUBLIC-EXPORT-INTEGRITY-TEST-01`, `HP-PUBLIC-PAPER-CI-TEST-01`, `HP-REP-PQS-RG-WORKING-CI-TEST-01`, `HP-REP-XGTO-INTERCHANGE-TEST-01`, `HP-REP-XGTO-READER-DOC-TEST-01`
- **Scope:** Completed v0.2.0-rc2 candidate-acceptance evidence only; this record grants no test, workflow, source, docs, tag, or release edit. Commit 2b3c23970144aa030ae52b875a5cf01b32886b6e passed clean Julia 1.10.12 and 1.12.6 package loads, the three public CI gates, export integrity, public residual-GTO/external-transfer and frozen H2/cc-pVTZ coverage, public examples 01/39/40/41, focused H2+ 18/18, screening 22/22, docs 99/99 and 10/10, and an archive/install replay excluding a root Manifest and both handoffs. The candidate archive contains 676 entries, is 9994240 bytes, and has SHA-256 6728b80c1397f13b367c2d898fbdda3176c6cb87c39597817c590a0f41c1e2ac. RC2/RC1/dev selector simulation retained no stable alias or path. The accepted C2 replay remains canonical and was not repeated because implementation was unchanged. The completed immutable annotated tag is recorded by HP-PQS-PUBLIC-RC2-TAG-FN-01; the exact GitHub prerelease operation is separately owned by HP-PQS-PUBLIC-RC2-RELEASE-FN-01. Registration, citation, stable documentation, and final v0.2 remain unauthorized.

### HP-PQS-PUBLIC-SCREEN-FN-01 - public supplied-field screened-Hartree assembly

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `Screened-Hartree Public Assembly`
- **Owned paths:**
  - `source` / `existing`: `src/GaussletBases.jl`
  - `source` / `existing`: `src/cartesian_reference_density/CartesianReferenceDensity.jl`
  - `source` / `existing`: `src/cartesian_reference_density/screened_hartree_correction.jl`
- **Evidence:**
  - `external_path`: `/Users/srw/Dropbox/Papers/PQS/notes/assignments/repo_design_manager_pqs_v0p2_public_surface_2026-08-17.md`
  - `git_commit`: `058ee54f45c759949f70b54a699ccc318476f8ac`
  - `manager_pass`: `481`
  - `manager_pass`: `482`
- **Dependencies:** `HP-PQS-SCREEN-HARTREE-CORR-FN-01`
- **Scope:** Maintain only the typed supplied-field public assembly: \`ExactRepresentedHartreeField\`, \`FittedReferenceHartreeField\`, existing \`ScreenedHartreeCorrection\`, \`screened\_hartree\_correction\`, and four correction accessors. Preserve one orthonormal basis/order, spin-summed P0 with fractional occupations allowed, public Coulomb-self-integral terminology, signed Tr(P0\*J0)-self-integral consistency, and density\_nonnegativity\_atol semantics. Atomic fitting, general terminal-plus-supplement field evaluation, packet placement/additivity, exchange, artifacts, solvers, and Ximg/XHF remain internal or forbidden.

### HP-PQS-PUBLIC-SCREEN-TEST-01 - public screened-Hartree assembly validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `Validation Contract`
- **Owned paths:**
  - `test` / `existing`: `test/ida/runtests.jl`
  - `test` / `existing`: `test/runtests.jl`
- **Evidence:**
  - `git_commit`: `058ee54f45c759949f70b54a699ccc318476f8ac`
  - `git_commit`: `72c46f9ea0dd6b2da7a6a302d34ea1c501d18647`
  - `manager_pass`: `481`
  - `manager_pass`: `482`
  - `manager_pass`: `483`
  - `manager_pass`: `484`
- **Dependencies:** `HP-PQS-PUBLIC-SCREEN-FN-01`, `HP-PQS-SCREEN-HARTREE-CORR-TEST-01`
- **Scope:** Maintain bounded public exact/fitted constructor, accessor, malformed-input, signed-consistency, energy/action-closure, and unchanged-physical-H1 tests plus the small one-center two-electron example. Obtain V\_IDA from a CartesianIDAHamiltonian and the accurate represented field from the existing two-orbital pure-GTO oracle. The candidate replay shares the archived fresh-resolution manifest and tracked-clean provenance required by HP-PQS-PUBLIC-MATCHED-TEST-01. No historical-He reproduction, packet fixture, external SCF, general four-index engine, exchange, solver, or mixed-basis field-construction claim.

### HP-PQS-PUBLIC-V020-FN-01 - v0.2.0 final candidate maintenance

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `docs`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `v0.2.0 Final Candidate And Conditional Publication Process`
- **Owned paths:**
  - `docs` / `existing`: `CHANGELOG.md`
  - `docs` / `existing`: `README.md`
  - `source` / `existing`: `Project.toml`
- **Evidence:**
  - `git_commit`: `1ab36c04753e0adaa03dafbaea4850d99b05ecba`
  - `git_commit`: `2b3c23970144aa030ae52b875a5cf01b32886b6e`
  - `git_commit`: `3cce96d40f9e4f06f23a190f782834271fca884b`
  - `git_commit`: `9ddc689c1bc806c7ec899cac7a39d77cb7fad3bf`
  - `git_commit`: `e1d5ca2ddb3a39134fddb476d00029ec590c431f`
  - `git_commit`: `fbef95e60ca3aadafe2082b4a63f8522a724e7be`
  - `git_commit`: `3419da6132810d8c4454f5b013c6302ef7842cb3`
  - `git_commit`: `94ec277d954b5435a04b0ad68ae352c95b0434c7`
  - `git_commit`: `b0dbd9ea37317590334a24883ef0667bdb0195a5`
  - `git_commit`: `adfcaba32d4db06d9d796d947276433717bd2d89`
  - `manager_pass`: `538`
  - `manager_pass`: `539`
- **Dependencies:** `HP-PQS-DOCS-TAGDEPLOY-FN-01`, `HP-PQS-PUBLIC-MATCHED-FN-01`, `HP-PQS-PUBLIC-RC2-FN-01`, `HP-PQS-PUBLIC-SCREEN-FN-01`, `HP-PUBLIC-EXPORT-INTEGRITY-FN-01`, `HP-PUBLIC-PAPER-CI-FN-01`, `HP-REP-PQS-RG-WORKING-FN-01`, `HP-REP-XGTO-INTERCHANGE-FN-01`
- **Scope:** Maintain only the exact accepted v0.2.0 candidate at commit adfcaba32d4db06d9d796d947276433717bd2d89 and tree f64ba21e06ff57e2b5e78d91214398115afbe8de: root version 0.2.0, concise post-RC2 CHANGELOG section above byte-identical RC2/RC1 history, README rev = "v0.2.0" installation and stable links with radial-first onboarding, and separate PQS, reference-density Hartree-screening, and external Cartesian-GTO stories. Preserve the modest package-nearest-software claim and accepted public numerics. The clean archive has 677 entries, 10137600 bytes, and SHA-256 df09cc6fd7dc144daa168c9feb4a41be9b974ef450e1e81bf586787318ad1566. The candidate itself authorizes no further edit, tag, release, registration, citation, or stable deployment. The exact conditional final transaction is separately owned by HP-PQS-PUBLIC-V020-RELEASE-FN-01. Preserve source, exports, APIs, dependencies/compat, examples, numerical behavior, workflows, docs/make.jl, fixture formats, manifest policy, and old changelog sections.

### HP-PQS-PUBLIC-V020-RELEASE-FN-01 - v0.2.0 final GitHub release publication

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `Final Publication Closeout`
- **Owned paths:** none
- **Evidence:**
  - `git_commit`: `adfcaba32d4db06d9d796d947276433717bd2d89`
  - `git_commit`: `d9e2189931039e039de0caafff5b18c6c696cec5`
  - `git_commit`: `4323ed4dbc9199bfec92d522ddc870f336eeb6e8`
  - `manager_pass`: `539`
  - `manager_pass`: `540`
  - `manager_pass`: `541`
- **Dependencies:** `HP-PQS-DOCS-TAGDEPLOY-FN-01`, `HP-PQS-PUBLIC-V020-FN-01`, `HP-PUBLIC-PAPER-CI-FN-01`
- **Scope:** Completed final GitHub release publication only; this record grants no further release, tag, file, workflow, repository-metadata, or asset mutation. Release 378216554 is visible at the immutable v0.2.0 tag with title GaussletBases v0.2.0, draft=false, prerelease=false, latest-final status, zero uploaded assets, and the exact 2278-byte ASCII body including final newline with SHA-256 e9ae9bcdad74b33bb66fb3e7e6a149d26285cb9bcc2f4c9555ac713be8bc90d2. Annotated tag object 722e8e8752a9d23f45e95d2f88e1749f9f3002e4 has message GaussletBases v0.2.0 and still peels to frozen commit adfcaba32d4db06d9d796d947276433717bd2d89 and tree f64ba21e06ff57e2b5e78d91214398115afbe8de. The fresh-scratch namespaced-ref lane proved tag/message/peel/tree/version and ls-remote/API identity after tag CI 33130411193 encountered only its local lightweight-tag collision. Docs 33130411176 and Pages 33130489319 passed; /v0.2.0/ and the real /stable/ alias are identical and canonicalize to /v0.2.0/, versions.js retains stable/v0.2/RC2/RC1/dev, and prior folders remain intact. Both automatic archives reconstruct the frozen tree, and a fresh isolated Julia 1.12.6 installation from rev = "v0.2.0" loaded GaussletBases 0.2.0 with that tree. No numerical gate was rerun. Preserve the published release and tag without editing, deleting, recreating, retargeting, moving, adding assets, or silent retry. Main-branch tag-lane repair, registration, citation metadata, later releases, and any source/API/test/example/dependency/workflow/numerical change require separate authority.

### HP-PQS-PUBLIC-V020-RELEASE-TEST-01 - v0.2.0 final publication validation

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `Final Publication Closeout`
- **Owned paths:** none
- **Evidence:**
  - `git_commit`: `adfcaba32d4db06d9d796d947276433717bd2d89`
  - `git_commit`: `4323ed4dbc9199bfec92d522ddc870f336eeb6e8`
  - `manager_pass`: `539`
  - `manager_pass`: `540`
  - `manager_pass`: `541`
- **Dependencies:** `HP-PQS-DOCS-TAGDEPLOY-TEST-01`, `HP-PQS-PUBLIC-V020-RELEASE-FN-01`, `HP-PQS-PUBLIC-V020-TEST-01`, `HP-PUBLIC-PAPER-CI-TEST-01`
- **Scope:** Completed final-publication evidence only; this record grants no file, release, tag, workflow, repository-metadata, or asset mutation. GitHub API and rendered-page checks confirm release 378216554 has exact tag v0.2.0, title GaussletBases v0.2.0, 2278-byte canonical body with final newline and SHA-256 e9ae9bcdad74b33bb66fb3e7e6a149d26285cb9bcc2f4c9555ac713be8bc90d2, draft=false, prerelease=false, latest-final status, and zero uploaded assets. The fresh namespaced verification ref, git ls-remote, GitHub API, and local tag agree on annotated object 722e8e8752a9d23f45e95d2f88e1749f9f3002e4, message GaussletBases v0.2.0, peel adfcaba32d4db06d9d796d947276433717bd2d89, tree f64ba21e06ff57e2b5e78d91214398115afbe8de, and version 0.2.0. Both automatic source archives reconstruct that tree; a fresh isolated Julia 1.12.6 tag installation loaded GaussletBases 0.2.0 with the same tree. /v0.2.0/, /stable/, selector entries, /dev/, RC2, and RC1 remain correct. No numerical gate was rerun. Preserve the accepted release and immutable tag without edit, deletion, recreation, retargeting, movement, asset upload, or silent retry.

### HP-PQS-PUBLIC-V020-TEST-01 - v0.2.0 final candidate acceptance

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `v0.2.0 Final Candidate And Conditional Publication Process`
- **Owned paths:** none
- **Evidence:**
  - `git_commit`: `1ab36c04753e0adaa03dafbaea4850d99b05ecba`
  - `git_commit`: `2b3c23970144aa030ae52b875a5cf01b32886b6e`
  - `git_commit`: `adfcaba32d4db06d9d796d947276433717bd2d89`
  - `manager_pass`: `538`
  - `manager_pass`: `539`
- **Dependencies:** `HP-PQS-DOCS-TAGDEPLOY-TEST-01`, `HP-PQS-PUBLIC-MATCHED-TEST-01`, `HP-PQS-PUBLIC-SCREEN-TEST-01`, `HP-PQS-PUBLIC-V020-FN-01`, `HP-PUBLIC-EXPORT-INTEGRITY-TEST-01`, `HP-PUBLIC-PAPER-CI-TEST-01`, `HP-REP-PQS-RG-WORKING-CI-TEST-01`, `HP-REP-XGTO-INTERCHANGE-TEST-01`, `HP-REP-XGTO-READER-DOC-TEST-01`
- **Scope:** Completed final-candidate acceptance evidence only; this record grants no test, workflow, source, docs, tag, or release edit. Candidate adfcaba32d4db06d9d796d947276433717bd2d89 and tree f64ba21e06ff57e2b5e78d91214398115afbe8de changed only Project.toml, CHANGELOG.md, README.md, and test/docs/runtests.jl. Direct GitHub install/load passed on Julia 1.10.12 and 1.12.6; fresh archive installation passed on 1.12.6; examples 01/39/40/41, H2+ 18/18, screening, export integrity, residual-GTO/external-transfer frozen fixture, authority, docs 114/114 + 10/10, Documenter, and diff checks passed. CI 33126022579 passed all three numerical gates and Docs 33126022531 passed. The clean archive has 677 entries, 10137600 bytes, SHA-256 df09cc6fd7dc144daa168c9feb4a41be9b974ef450e1e81bf586787318ad1566, no root Manifest, and neither handoff. The accepted exact ASCII release body has 2278 bytes including its final newline and SHA-256 e9ae9bcdad74b33bb66fb3e7e6a149d26285cb9bcc2f4c9555ac713be8bc90d2. The docs-test diff was +38/-9: nine one-for-one replacements and 29 substantive additions within the hard 30-line bound. The exact conditional final transaction is separately owned by HP-PQS-PUBLIC-V020-RELEASE-FN-01.

### HP-PQS-READER-DOC-01 - reader-facing PQS documentation entrance

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `docs`
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [r1\_public\_base\_producer.md](r1_public_base_producer.md); heading `Reader-Facing PQS Entrance`
- **Owned paths:**
  - `docs` / `existing`: `docs/src/algorithms/cartesian_ida_overview.md`
  - `docs` / `existing`: `docs/src/algorithms/ida_hamiltonian_and_counterpoise.md`
  - `docs` / `existing`: `docs/src/algorithms/pqs_shell_construction.md`
  - `docs` / `existing`: `docs/src/developer/designs/cartesian_hamiltonian_producer/r1_public_base_producer.md`
  - `docs` / `existing`: `docs/src/examples/index.md`
  - `docs` / `existing`: `docs/src/howto/example_guide.md`
  - `docs` / `existing`: `docs/src/index.md`
  - `docs` / `existing`: `docs/src/manual/index.md`
  - `docs` / `existing`: `docs/src/manual/projected_q_shells.md`
- **Evidence:**
  - `git_commit`: `01f0fe002`
  - `manager_pass`: `479`
  - `manager_pass`: `480`
- **Dependencies:** `HP-PQS-ASPECTSHELL-FN-01`, `HP-R1-FN-01`
- **Scope:** Maintain only the implemented reader-facing PQS entrance, its visible links to the three owned PQS/IDA algorithm pages, and radial-first onboarding. Preserve the bounded one-center and bond-aligned homonuclear geometry statement and construction-smoke interpretation. Add no source/API, parser, documentation framework, duplicated algorithm account, artifact, solver, screening, supplement, PRF, paper-driver, manuscript result, or broad molecular claim.

### HP-PQS-READER-TEST-01 - public PQS H2+ example validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r1\_public\_base\_producer.md](r1_public_base_producer.md); heading `Reader-Facing PQS Entrance`
- **Owned paths:**
  - `test` / `existing`: `test/docs/runtests.jl`
  - `test` / `existing`: `test/runtests.jl`
- **Evidence:**
  - `git_commit`: `01f0fe002`
  - `manager_pass`: `479`
  - `manager_pass`: `480`
- **Dependencies:** `HP-PQS-ASPECTSHELL-FN-01`, `HP-PQS-READER-DOC-01`, `HP-R1-FN-01`, `HP-R1-TEST-01`
- **Scope:** Maintain only the committed examples/39\_pqs\_h2plus.jl public-only fixture, its one quick-example invocation, focused documentation link/public-name checks, the existing Julia-1.10 examples-group selection, and the canonical 1e-10 symmetry/eigen-residual gates. Add no test file, source/API/default, private call, parser, artifact, solver, supplement, PRF, screening, numerical-order gate, workflow framework, or other CI behavior.

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
- **Evidence:**
  - `git_commit`: `4de885c90`
  - `manager_pass`: `414`
- **Dependencies:** `HP-PQS-ASPECTSHELL-FN-01`, `HP-RG-NUMCOMP-FN-01`
- **Scope:** Maintain \`owner = :all\` overrides for positive semantic \`:atom\_local\_shell\` or \`:shared\_molecular\_shell\` indices with non-Boolean integer \`source\_q \>= 3\` and \`source\_q \!= route\_q\`. Values above route q refine the shell contraction; values below route q coarsen it while parent axes, support, ownership, cores, slabs, and route metadata remain unchanged. Atom-local shells use \`(source\_q,source\_q,source\_q)\`; shared shells rerun the existing angular-band selector for \`(source\_q,source\_q,L)\`. Selector retention, \`nside\`, and \`selected\_q\` use \`source\_q\`. One authoritative shape must reach existing retained/support/transform/realization/due-diligence consumers unchanged.

### HP-PQS-SHELLQ-OVERRIDE-TEST-01 - semantic source-q override validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_semantic\_shell\_q\_overrides.md](pqs_semantic_shell_q_overrides.md); heading `Semantic Per-Shell PQS Source-q Overrides`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:**
  - `git_commit`: `4de885c90`
  - `manager_pass`: `414`
- **Dependencies:** `HP-PQS-SHELLQ-OVERRIDE-FN-01`
- **Scope:** Validate route-q 7 to source-q 6 and 5 coarsening, expected retained-count reduction, unchanged parent/support/ownership/cores/slabs/route metadata, orthonormal contraction columns, omitted/empty parity, finite symmetric full construction, and rejection of malformed, below-3, Boolean, equal-route, asymmetric, and unmatched requests. Preserve existing refinement, residual, packet-capture, \`J0/E0\`, correction, dimension, and due-diligence gates. No new accessor, dense baseline-to-variant overlap API, source-pass HF, endpoint energy assertion, or production claim is approved.

### HP-PUBLIC-EXPORT-INTEGRITY-FN-01 - reduce unsupported package exports

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `Package Export Integrity`
- **Owned paths:**
  - `source` / `existing`: `src/GaussletBases.jl`
  - `source` / `existing`: `src/cartesian_nested_faces.jl`
- **Evidence:**
  - `git_commit`: `41ab3e13121f5af1e145775500e91f9ac61c9760`
  - `git_commit`: `b72500f7e619db5875918e3290ed2b306be51f43`
  - `git_commit`: `6e7bcbb7dae4e865dbdc0362b8f39ffd23f0a468`
  - `manager_pass`: `499`
  - `manager_pass`: `500`
  - `manager_pass`: `543`
  - `manager_pass`: `545`
- **Dependencies:** none
- **Scope:** Maintain the post-v0.2 export cleanup implemented by commit 6e7bcbb7dae4e865dbdc0362b8f39ffd23f0a468. TimedNestedFixedBlockBuild and its root export remain absent. OneCenterAtomicNestedLayerStructure and QiuWhiteResidualGaussianOperators remain defined but unexported; retain the one-line QiuWhiteResidualGaussianOperators alias to OrdinaryCartesianOperators3D for qualified compatibility. The accepted production delta is +0/-7 across src/GaussletBases.jl and src/cartesian\_nested\_faces.jl, and the undocumented exported-binding backlog is 71. Preserve CartesianBasisBundle3D and nested\_fixed\_block\_timing\_report absence, final\_units and unit\_keys internal-export absence, lowering\_recipe, ShellLocalAngularProfileKey, @timeg, CuratedSpherePointSet, LegacySGaussianData, QiuWhiteHybridOrbital3D, every ordinary/QW/nested numerical result, and immutable v0.2.0 tags/releases. Add no shim, warning, replacement API, alias deletion, docstring, test, helper, metadata, file, dependency, numerical change, workflow, version, tag, release, registration, or citation change. Ignored historical probes and conflicted copies remain untouched.

### HP-PUBLIC-EXPORT-INTEGRITY-TEST-01 - validate package export integrity

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `Package Export Integrity`
- **Owned paths:**
  - `test` / `existing`: `test/core/runtests.jl`
- **Evidence:**
  - `git_commit`: `41ab3e13121f5af1e145775500e91f9ac61c9760`
  - `git_commit`: `b72500f7e619db5875918e3290ed2b306be51f43`
  - `git_commit`: `6e7bcbb7dae4e865dbdc0362b8f39ffd23f0a468`
  - `manager_pass`: `499`
  - `manager_pass`: `500`
  - `manager_pass`: `543`
  - `manager_pass`: `545`
- **Dependencies:** `HP-PUBLIC-EXPORT-INTEGRITY-FN-01`
- **Scope:** Maintain the compact dynamic regression in test/core/runtests.jl unchanged. It audits GaussletBases and every defined direct package-owned child module discovered with names(GaussletBases; all=true, imported=false) and parentmodule(child) === GaussletBases. For each name returned by names(module; all=false, imported=false), require the binding to be defined; when the value is a Function, require methods(value) to be nonempty. Commit 6e7bcbb7dae4e865dbdc0362b8f39ffd23f0a468 passed this owner while runtime inspection proved TimedNestedFixedBlockBuild absent, OneCenterAtomicNestedLayerStructure and QiuWhiteResidualGaussianOperators defined but unexported, and the undocumented export backlog reduced from 74 to 71. Preserve this invariant without editing tests, snapshotting an export list, classifying broader use, traversing external/imported modules, or adding numerical assertions. Add no test file, fixture, dependency, helper, compatibility assertion, version/release check, or broad architecture audit under this record.

### HP-PUBLIC-PAPER-CI-FN-01 - paper-aligned PQS and screening CI workflow

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `tools`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `Paper-Aligned CI Boundary`
- **Owned paths:**
  - `tool` / `existing`: `.github/workflows/ci.yml`
- **Evidence:**
  - `git_commit`: `409bc41b46afced92e3e711b8a55760869ac5d3d`
  - `git_commit`: `a4e85e820fd4056e985a18e20da87180f370ef66`
  - `git_commit`: `e65764377bd4640916e80342071da754d80aca32`
  - `git_commit`: `3cce96d40f9e4f06f23a190f782834271fca884b`
  - `git_commit`: `72446880603c7e554f6ae71b2de2dc6edf28b31b`
  - `git_commit`: `9ddc689c1bc806c7ec899cac7a39d77cb7fad3bf`
  - `git_commit`: `15676153aec1569f5224ffa6ff5ed67b054c837f`
  - `manager_pass`: `503`
  - `manager_pass`: `504`
  - `manager_pass`: `505`
  - `manager_pass`: `525`
  - `manager_pass`: `527`
  - `manager_pass`: `528`
  - `manager_pass`: `529`
  - `manager_pass`: `536`
  - `manager_pass`: `537`
  - `manager_pass`: `542`
  - `manager_pass`: `544`
- **Dependencies:** `HP-PQS-DOCS-TAGDEPLOY-FN-01`, `HP-PQS-PUBLIC-COMPAT-FN-01`, `HP-PQS-PUBLIC-MATCHED-TEST-01`, `HP-PQS-PUBLIC-SCREEN-TEST-01`
- **Scope:** Maintain the implemented three-gate workflow with the exact Julia 1.10 Supported-floor selection core,ida,cartesian,examples,radial,misc. Commit 15676153aec1569f5224ffa6ff5ed67b054c837f changed only that existing group list; preserve the workflow file, three job names and rows, Julia versions, 30-minute timeout, path classifier and four-path allowlist, documentation-only markers and lightweight lane, commands, triggers, permissions, pull-request behavior, PQS and Screening groups, annotated-tag identity/install lane, disabled slow tests, and fail-closed behavior. Angular remains explicitly excluded pending a separate fast-versus-acceptance audit after its first package-owned one-body fixture ran 13m49s without reaching an assertion. Add no row, job, workflow, action dependency, helper, command, trigger, permission, source, numerical policy, dependency, manifest, version, tag, release, registration, citation, or other test-owner change.

### HP-PUBLIC-PAPER-CI-TEST-01 - paper-aligned PQS and screening release validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_public\_surface.md](pqs_public_surface.md); heading `Paper-Aligned CI Boundary`
- **Owned paths:**
  - `test` / `existing`: `test/docs/runtests.jl`
  - `test` / `existing`: `test/pqs_h2plus_table1_release_runtests.jl`
  - `test` / `existing`: `test/runtests.jl`
- **Evidence:**
  - `git_commit`: `409bc41b46afced92e3e711b8a55760869ac5d3d`
  - `git_commit`: `a4e85e820fd4056e985a18e20da87180f370ef66`
  - `git_commit`: `e65764377bd4640916e80342071da754d80aca32`
  - `git_commit`: `3cce96d40f9e4f06f23a190f782834271fca884b`
  - `git_commit`: `72446880603c7e554f6ae71b2de2dc6edf28b31b`
  - `git_commit`: `9ddc689c1bc806c7ec899cac7a39d77cb7fad3bf`
  - `git_commit`: `b0dbd9ea37317590334a24883ef0667bdb0195a5`
  - `git_commit`: `15676153aec1569f5224ffa6ff5ed67b054c837f`
  - `manager_pass`: `503`
  - `manager_pass`: `504`
  - `manager_pass`: `505`
  - `manager_pass`: `525`
  - `manager_pass`: `527`
  - `manager_pass`: `528`
  - `manager_pass`: `529`
  - `manager_pass`: `536`
  - `manager_pass`: `537`
  - `manager_pass`: `542`
  - `manager_pass`: `544`
- **Dependencies:** `HP-PQS-DOCS-TAGDEPLOY-TEST-01`, `HP-PQS-PUBLIC-MATCHED-TEST-01`, `HP-PQS-PUBLIC-SCREEN-TEST-01`, `HP-PUBLIC-PAPER-CI-FN-01`
- **Scope:** Maintain only the focused documentation-policy assertion requiring the exact Supported-floor selection core,ida,cartesian,examples,radial,misc. Commit 15676153aec1569f5224ffa6ff5ed67b054c837f added that assertion without changing numerical coverage. Preserve every numerical assertion, tolerance, fixture, pqs\_release single-execution contract, screening\_release contract, required-check identity, all three CI names/rows, path-aware routing, workflow command, and documentation link. Radial passed 322/322, misc passed 59/59, their combined selection passed 381/381, and remote CI run 33141930944 passed the expanded Julia 1.10 Supported-floor gate plus unchanged PQS and Screening gates. Add no numerical assertion, fixture, helper, cache, artifact, file, framework, group, row, workflow, source, API, dependency, docs, version, tag, or release change.

### HP-QW-NESTED-DIAT-FN-01 - repair exported ordinary-QW nested diatomic front doors

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_nested\_diatomic\_box\_policy.md](../../../algorithms/cartesian_nested_diatomic_box_policy.md); heading `Exported Front-Door Maintenance Contract`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_nested_diatomic.jl`
  - `source` / `existing`: `src/ordinary_qw_nested_frontends.jl`
- **Evidence:**
  - `git_commit`: `6a365699102ebe30c31a9499f338cb192334ef1c`
  - `manager_pass`: `472`
  - `manager_pass`: `474`
- **Dependencies:** none
- **Scope:** Maintain the three existing root-exported bond-aligned ordinary-QW nested diatomic front doors: a legal unsplit/no-shared-layer source owns the packet required by fixed-block construction; shared shells retain their actual rectangular or endcap/panel provenance without invalid concrete-type coercion; and basis-level geometry diagnostics forward the existing shared-shell policy, q, and L controls. Preserve source reuse, signatures, exports, defaults, coefficients, dimensions, geometry, shell policy, and numerical operators. Add no helper file, type family, field, metadata, status, keyword, export, adapter, fallback, PQS/WL change, or broader ordinary-QW work.

### HP-QW-NESTED-DIAT-TEST-01 - validate ordinary-QW nested diatomic front-door repairs

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [cartesian\_nested\_diatomic\_box\_policy.md](../../../algorithms/cartesian_nested_diatomic_box_policy.md); heading `Exported Front-Door Maintenance Contract`
- **Owned paths:**
  - `test` / `existing`: `test/core/runtests.jl`
- **Evidence:**
  - `git_commit`: `6a365699102ebe30c31a9499f338cb192334ef1c`
  - `manager_pass`: `472`
  - `manager_pass`: `474`
- **Dependencies:** `HP-QW-NESTED-DIAT-FN-01`
- **Scope:** Maintain exactly two bounded regressions in the existing core owner: a legal unsplit H2 case proving packet availability and source/fixed-block reuse, and an endcap/panel diagnostics case proving actual provenance plus nondefault existing policy/q/L forwarding. Preserve the default baseline checks for unchanged dimensions, coefficients, geometry, overlap, and kinetic operators. Add no test file, historical nested-suite wiring, metadata/status assertion, PQS/WL fixture, or numerical-policy change.

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
- **Scope:** maintain explicit origin-centered all-electron atom geometry and input validation in the existing \`cartesian\_base\_hamiltonian(system; basis, hamfile)\` facade. Charge, electron counts, spin sectors, basis, and ECP behavior must never be inferred from the atom label. Charged-sector acceptance is owned separately by HP-R1-ESECTOR-FN-01.

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

### HP-R1-ESECTOR-FN-01 - explicit charged electron sectors

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r1\_public\_base\_producer.md](r1_public_base_producer.md); heading `Electron-Sector Independence`
- **Owned paths:**
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
- **Evidence:**
  - `git_commit`: `540cef00b`
  - `manager_pass`: `446`
- **Dependencies:** `HP-COMP-BASEDIAT-FN-01`, `HP-R1-ATOM-FN-01`, `HP-R1-FN-01`, `HP-R3U-ZDI-FN-01`
- **Scope:** Remove neutrality-derived and integer-charge validation from the existing base and supplemented atom/homonuclear-z-diatomic normalization in src/cartesian\_base\_hamiltonian.jl. Preserve explicit required nup/ndn, positive total electron count, orbital-dimension validation, geometry restrictions, neutral-call behavior, and all basis/operator algorithms. At fixed nuclei and basis, changing only nup/ndn must leave parent-axis, terminal-coefficient/support, supplement/residual, kinetic, unit-nuclear, H1, and Vee numerical arrays plus nuclear repulsion exact; containers and provenance may differ only in explicit sector-derived fields. Delete the obsolete integer-charge helper if it has no live caller. Add no API, type, field, artifact key/schema, cache key, correction framework, driver branch, source helper/file, or zero-electron-sector support. Maximum 20 added source lines; net source should decrease.

### HP-R1-ESECTOR-TEST-01 - charged-sector independence validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [r1\_public\_base\_producer.md](r1_public_base_producer.md); heading `Electron-Sector Independence`
- **Owned paths:**
  - `test` / `existing`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:**
  - `git_commit`: `540cef00b`
  - `manager_pass`: `446`
- **Dependencies:** `HP-R1-ESECTOR-FN-01`
- **Scope:** In existing test owners only, prove neutral/charged He and He2 basis/operator exact parity, supplemented He2/He2(2+) exact operator parity, charged artifact/readback sector preservation, unchanged neutral endpoints, and rejection of malformed or invalid sectors. Require bitwise equality where arrays share the same deterministic construction; report any non-bitwise difference and stop rather than adding tolerances or sector-dependent branches. Maximum 60 added test lines and no new test/probe/fixture file.

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
  - `manager_pass`: `475`
  - `repo_path`: `test/driver_public/cartesian_base_hamiltonian_runtests.jl`
- **Dependencies:** none
- **Scope:** maintain the public base producer in its current owner file, the existing \`cartesian\_base\_hamiltonian\` export, and the existing expert/unstable \`cartesian\_base\_working\_basis\` export. The working-basis return fields are not a stable public schema, and this ownership does not export PRF-specific wrappers.

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
- **Scope:** maintain explicit all-electron homonuclear two-center z-axis geometry, system validation, and optional trusted supplement \`basisfile\`. Charged-sector acceptance and sector-independent operator parity are owned separately by HP-R1-ESECTOR-FN-01.

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

### HP-REP-MIXDENS-HARTREE-FN-01 - represented multicenter mixed-density direct Hartree producer

- **Lifecycle:** `approved`
- **Grant:** `implementation`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [represented\_mixed\_density\_hartree.md](represented_mixed_density_hartree.md); heading `Represented Mixed-Density Hartree Producer`
- **Owned paths:**
  - `source` / `existing`: `src/GaussianAnalyticIntegrals.jl`
  - `source` / `existing`: `src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl`
  - `source` / `existing`: `src/cartesian_reference_density/CartesianReferenceDensity.jl`
  - `source` / `existing`: `src/cartesian_reference_density/represented_molecular_hartree.jl`
  - `source` / `planned`: `src/cartesian_reference_density/represented_hartree_contractions.jl`
  - `source` / `existing`: `src/cartesian_residual_gaussians/residual_basis.jl`
  - `source` / `existing`: `src/gaussian_coulomb_reference.jl`
- **Evidence:**
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/Papers/PQS/validation/cr2_req084_molecular_full_hartree_capability_2026-08-13.md`
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/Papers/PQS/validation/work/REQ-084/repo_design_manager_assignment.md`
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/cr2/reports/req084_cr2_molecular_full_hartree_capability_2026-08-13.md`
  - `manager_pass`: `459`
  - `manager_pass`: `460`
  - `git_commit`: `a77ceed5d1b1f074c04a2b11cda2b962be2d47a7`
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/Papers/PQS/validation/work/REQ-084/repo_design_manager_scaling_assignment.md`
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/cr2/reports/req084_cr2_molecular_full_hartree_continuation_preflight_2026-08-14.md`
  - `manager_pass`: `461`
  - `manager_pass`: `462`
- **Dependencies:** `HP-RG-FN-01`, `HP-RG-FN-02`, `HP-RHO0-MIXH-GAAA-FN-01`, `HP-RHO0-MIXH-GG-FN-01`
- **Scope:** Replace the bounded global component-pair production loop with the exact occupied-contracted block/separable action for complete GG, GA/AG, and AA source and target sectors under every expansion term. Revalidate residual cross/identity with the authoritative RG 1e-10/scale-aware 5e-8 contract while keeping state charge/Gram at 1e-10. The correction may add 650 preferred/800 hard source lines but must delete or hard-bound the old unbounded path in the same commit. No fit, screening, truncation, IDA substitute, residual mutation, public API, Cr2 branch, or partial scaffolding.

### HP-REP-MIXDENS-HARTREE-TEST-01 - represented multicenter mixed-density Hartree validation

- **Lifecycle:** `approved`
- **Grant:** `implementation`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [represented\_mixed\_density\_hartree.md](represented_mixed_density_hartree.md); heading `Represented Mixed-Density Hartree Producer`
- **Owned paths:**
  - `test` / `existing`: `test/nested/cartesian_represented_molecular_hartree_runtests.jl`
- **Evidence:** none
- **Dependencies:** `HP-REP-MIXDENS-HARTREE-FN-01`
- **Scope:** Extend only the existing bounded mixed-basis test owner by 120 preferred/180 hard lines, with final length at most 360: compare occupied-contracted block/separable GG, GA/AG, and AA sectors and complete native fields with the bounded component and independent Gaussian oracles; distinguish residual 1e-10 cross/scale-aware 5e-8 identity validity from strict 1e-10 state recovery; and validate deterministic two-size contraction/resource accounting plus unchanged atomic and screened owners. No Cr2 fixture, fit, solver, or metadata/status test.

### HP-REP-PQS-RG-WORKING-CI-FN-01 - PQS residual-GTO public Cartesian CI wiring

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_residual\_gto\_working\_basis.md](pqs_residual_gto_working_basis.md); heading `Public Cartesian CI Extraction`
- **Owned paths:**
  - `test` / `existing`: `test/runtests.jl`
- **Evidence:**
  - `git_commit`: `346589a6dafd8fbfe1794947dc8cb7da4d07606b`
  - `git_commit`: `0adac132c213efc712d27347b043dd57c204d745`
  - `git_commit`: `50b5ee9a1fd77bbcf7b8e8e95b3f0a79a53d11bc`
  - `manager_pass`: `506`
  - `manager_pass`: `507`
- **Dependencies:** `HP-REP-PQS-RG-WORKING-FN-01`
- **Scope:** Maintain exactly one include for test/driver\_public/cartesian\_residual\_gto\_mwg\_system\_runtests.jl in the existing :cartesian group in test/runtests.jl, as accepted in 50b5ee9a1fd77bbcf7b8e8e95b3f0a79a53d11bc. The Julia 1.10 Supported-floor row owns the regression through the unchanged cartesian selection. Add no CI row, workflow edit, group, example smoke, nested-suite wiring, standalone fixture/helper file, or source/API behavior; inline public test construction is owned only by HP-REP-PQS-RG-WORKING-CI-TEST-01. Do not alter the separately named PQS and screening paper gates.

### HP-REP-PQS-RG-WORKING-CI-TEST-01 - PQS residual-GTO public Cartesian validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_residual\_gto\_working\_basis.md](pqs_residual_gto_working_basis.md); heading `Public Cartesian CI Extraction`
- **Owned paths:**
  - `test` / `existing`: `test/driver_public/cartesian_residual_gto_mwg_system_runtests.jl`
  - `test` / `existing`: `test/nested/cartesian_external_gto_import_runtests.jl`
- **Evidence:**
  - `git_commit`: `346589a6dafd8fbfe1794947dc8cb7da4d07606b`
  - `git_commit`: `0adac132c213efc712d27347b043dd57c204d745`
  - `git_commit`: `50b5ee9a1fd77bbcf7b8e8e95b3f0a79a53d11bc`
  - `manager_pass`: `506`
  - `manager_pass`: `507`
  - `manager_pass`: `508`
  - `git_commit`: `c633d5db0e51c1795f4a4ea194929eca3fc53b69`
  - `manager_pass`: `509`
- **Dependencies:** `HP-REP-PQS-RG-WORKING-CI-FN-01`, `HP-REP-PQS-RG-WORKING-TEST-01`, `HP-REP-XGTO-IMPORT-TEST-01`
- **Scope:** Maintain the accepted 136-line, 47-check public residual-GTO/import owner from c633d5db0e51c1795f4a4ea194929eca3fc53b69: one paid H2 construction; normalized matched \`px\`/\`py\` probes selected through public representation APIs with \`S\_GG = I\`; basic and equal-occupation rotated imports; explicit alpha/beta imports; ordering/metric/fingerprint identity and rejection; invalid spin combinations; and nonorthonormal source-coefficient rejection. Keep the public owner free of private names and duplicate restricted-import clusters. Preserve the complete 49-check protected-sidecar testset and its required helpers as direct-run internal evidence. The accepted extraction added/deleted 75/7 public-owner lines and deleted 111 nested lines, for total test delta +75/-118 and net -43. Add no source, API, fixture file, runner, workflow, CI row, dependency, numerical policy, duplicate nested coverage, or other nested-suite wiring under maintenance.

### HP-REP-PQS-RG-WORKING-FN-01 - PQS residual-GTO same-construction working basis

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `docs`, `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_residual\_gto\_working\_basis.md](pqs_residual_gto_working_basis.md); heading `PQS Residual-GTO Same-Construction Working Basis`
- **Owned paths:**
  - `source` / `existing`: `src/GaussletBases.jl`
  - `source` / `existing`: `src/cartesian_base_hamiltonian.jl`
  - `source` / `existing`: `src/cartesian_gto_probes.jl`
  - `docs` / `existing`: `docs/src/reference/export.md`
- **Evidence:**
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/Papers/PQS/validation/requests/REQ-101_amendment_public_external_gto_import.md`
  - `external_path`: `/Users/srw/Library/CloudStorage/Dropbox/Papers/PQS/validation/work/REQ-101/public_import_resume/preflight.md`
  - `git_commit`: `346589a6dafd8fbfe1794947dc8cb7da4d07606b`
  - `manager_pass`: `501`
  - `manager_pass`: `502`
- **Dependencies:** `HP-R3U-FN-01`, `HP-R3U-WIRE-01`, `HP-REP-XGTO-IMPORT-FN-01`
- **Scope:** Maintain the expert \`cartesian\_residual\_gto\_mwg\_system\` constructor and its unexported concrete same-construction result with stable \`.hamiltonian\`, \`gto\_overlap\_matrix(result, probes)\`, and unchanged external-import behavior. Retain only terminal/factorized-parent/supplement/residual data and preserve native \`S\_RX = T\_G' \* S\_GX + T\_A' \* S\_AX\`, exact direct-facade/artifact parity, block selection, finiteness and dimension failures, and the opaque overlap-only boundary accepted in 346589a6dafd8fbfe1794947dc8cb7da4d07606b. Add no new file, exported type, generalized representation, persistence, parser, solver, screening automation, dense parent-to-terminal/raw-to-final map, or paper-specific branch under this record.

### HP-REP-PQS-RG-WORKING-TEST-01 - PQS residual-GTO working-basis validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [pqs\_residual\_gto\_working\_basis.md](pqs_residual_gto_working_basis.md); heading `PQS Residual-GTO Same-Construction Working Basis`
- **Owned paths:**
  - `test` / `existing`: `test/docs/runtests.jl`
  - `test` / `existing`: `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`
- **Evidence:**
  - `git_commit`: `346589a6dafd8fbfe1794947dc8cb7da4d07606b`
  - `manager_pass`: `501`
  - `manager_pass`: `502`
- **Dependencies:** `HP-REP-PQS-RG-WORKING-FN-01`, `HP-REP-XGTO-IMPORT-TEST-01`
- **Scope:** Maintain the accepted private-oracle and documentation checks in the existing supplemented-H2 and docs owners: exact direct-Hamiltonian/artifact parity, native \`S\_GX\`/\`S\_AX\`/\`S\_RX\` formula and private assembly-boundary validation, artifact readback/provenance, terminal due diligence, and curated-reference parity. Public-only result, overlap, packet-import, and malformed-input coverage transfers solely under HP-REP-PQS-RG-WORKING-CI-\*. Preserve unique private evidence while removing duplicate public assertions. REQ-101 remains external acceptance; add no C2 fixture, generalized-representation assertion, solver, or screening policy under this record.

### HP-REP-XGTO-CLOSESTDET-FN-01 - external Cartesian GTO closest determinant

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [external\_cartesian\_gto\_interchange.md](external_cartesian_gto_interchange.md); heading `External Cartesian GTO Interchange`
- **Owned paths:**
  - `source` / `existing`: `src/GaussletBases.jl`
  - `source` / `existing`: `src/cartesian_external_gto_interchange.jl`
- **Evidence:**
  - `git_commit`: `bd1d11cec71c2257b311054bb064e361bb533042`
  - `manager_pass`: `515`
  - `manager_pass`: `516`
- **Dependencies:** `HP-REP-PQS-RG-WORKING-FN-01`, `HP-REP-XGTO-INTERCHANGE-FN-01`
- **Scope:** Maintain the accepted \`closest\_external\_gto\_determinant(working, packet; minimum\_gram\_eigenvalue)\` operation from bd1d11cec71c2257b311054bb064e361bb533042. Preserve the existing raw \`S\_FG\*C\_G\` import unchanged, the separate full-rank symmetric-Lowdin determinant, Gram/principal-angle/orthonormality diagnostics, explicit caller threshold, and failures on determinant occupation, geometry, electron sector, supported Hamiltonian, or rank mismatch. Accept only a same-construction working object exposing its Hamiltonian; add no solver state, detached-coefficient overload, separately supplied Hamiltonian, public type, flooring, dropping, appending, fractional determinant semantics, or automatic low-capture acceptance.

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
  - `test` / `existing`: `test/driver_public/cartesian_residual_gto_mwg_system_runtests.jl`
  - `test` / `existing`: `test/nested/cartesian_external_gto_import_runtests.jl`
- **Evidence:**
  - `git_commit`: `c633d5db0e51c1795f4a4ea194929eca3fc53b69`
  - `manager_pass`: `509`
- **Dependencies:** `HP-REP-XGTO-IMPORT-FN-01`
- **Scope:** Maintain the accepted split of packet identity/order/\`S\_GG\`, restricted/spin-resolved import, capture, rotation invariance, and malformed-input checks from c633d5db0e51c1795f4a4ea194929eca3fc53b69. Public importer assertions live only in the wired residual-GTO owner under HP-REP-PQS-RG-WORKING-CI-TEST-01; the protected-sidecar testset remains direct-run under HP-REP-XGTO-PROTECT-SIDECAR-TEST-01. Restore no duplicate general-import nested testset.

### HP-REP-XGTO-INTERCHANGE-FN-01 - external Cartesian GTO interchange reader

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [external\_cartesian\_gto\_interchange.md](external_cartesian_gto_interchange.md); heading `External Cartesian GTO Interchange`
- **Owned paths:**
  - `source` / `existing`: `src/GaussletBases.jl`
  - `source` / `existing`: `src/cartesian_external_gto_interchange.jl`
- **Evidence:**
  - `git_commit`: `bd1d11cec71c2257b311054bb064e361bb533042`
  - `manager_pass`: `515`
  - `manager_pass`: `516`
- **Dependencies:** `HP-REP-XGTO-IMPORT-FN-01`
- **Scope:** Maintain the strict version-1 \`read\_external\_cartesian\_gto\_packet(path; overlap\_atol, overlap\_rtol)\` reader accepted in bd1d11cec71c2257b311054bb064e361bb533042. Preserve the same-stem TOML/little-endian column-major Float64 bundle, ordered explicit Cartesian probes, positive diagonal AO normalization reconciliation, complete source-overlap parity through the existing overlap kernel, unchanged MO coefficients, and reuse of \`ExternalGTOOrbitalPacket\`. Add no second overlap implementation, metric solve, runtime permutation ledger, label sorting, basis lookup, public type, mandatory dependency, alternate packet/import path, Hamiltonian payload, checkpoint/HDF5/Molden/QCSchema parser, or artifact schema.

### HP-REP-XGTO-INTERCHANGE-TEST-01 - external Cartesian GTO interchange validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [external\_cartesian\_gto\_interchange.md](external_cartesian_gto_interchange.md); heading `External Cartesian GTO Interchange`
- **Owned paths:**
  - `test` / `existing`: `test/driver_public/cartesian_residual_gto_mwg_system_runtests.jl`
  - `test` / `existing`: `test/driver_public/external_cartesian_gto_h2_ccpvtz_v1.f64`
  - `test` / `existing`: `test/driver_public/external_cartesian_gto_h2_ccpvtz_v1.toml`
- **Evidence:**
  - `git_commit`: `0ca4e18c7ff550ec6f254f63b9f075193094ccd0`
  - `git_commit`: `bd1d11cec71c2257b311054bb064e361bb533042`
  - `external_path`: `/Users/srw/Dropbox/Papers/PQS/validation/work/REQ-101/hfdmrg_seed_viability_resume/run_single_pass_case.jl`
  - `external_path`: `/Users/srw/Dropbox/Papers/PQS/validation/work/REQ-101/req046_public_gate/audit_bridge.py`
  - `external_path`: `/Users/srw/Dropbox/Papers/PQS/validation/work/REQ-101/schema_adapter_resume_v2/replay_and_packet.py`
  - `manager_pass`: `515`
  - `manager_pass`: `516`
- **Dependencies:** `HP-REP-XGTO-CLOSESTDET-FN-01`, `HP-REP-XGTO-INTERCHANGE-FN-01`, `HP-REP-XGTO-PYSCF-EXPORT-FN-01`
- **Scope:** Maintain the accepted frozen two-center H2/cc-pVTZ Cartesian d-shell bundle and its wired public validation from bd1d11cec71c2257b311054bb064e361bb533042 and 0ca4e18c7ff550ec6f254f63b9f075193094ccd0. Preserve payload/AO integrity, full source-overlap parity, source MO metric identity, exact existing-import reconstruction, occupied capture, closest-determinant orthogonality with unchanged raw projection, malformed bundle/metric/order/occupation/geometry/sector/threshold rejection, multi-column rotation covariance, and exporter-version identity. The accepted read-only C2/aug-cc-pV6Z replay remains external evidence summarized in the canonical contract; add no large packet, NWChem snapshot, permutation ledger, new test owner, workflow, scheduled regeneration, private-helper call, paper identifier, or production/release claim.

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

### HP-REP-XGTO-PYSCF-EXPORT-FN-01 - PySCF external Cartesian GTO exporter

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `driver`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [external\_cartesian\_gto\_interchange.md](external_cartesian_gto_interchange.md); heading `External Cartesian GTO Interchange`
- **Owned paths:**
  - `driver` / `existing`: `bin/export_pyscf_cartesian_gto.py`
- **Evidence:**
  - `git_commit`: `0ca4e18c7ff550ec6f254f63b9f075193094ccd0`
  - `git_commit`: `bd1d11cec71c2257b311054bb064e361bb533042`
  - `manager_pass`: `515`
  - `manager_pass`: `516`
- **Dependencies:** `HP-REP-XGTO-INTERCHANGE-FN-01`
- **Scope:** Maintain the optional checkpoint-to-v1-bundle command \`bin/export\_pyscf\_cartesian\_gto.py\` accepted in bd1d11cec71c2257b311054bb064e361bb533042 and fixture-regenerated in 0ca4e18c7ff550ec6f254f63b9f075193094ccd0. Preserve exact PySCF Cartesian AO order and explicit primitive data; for real spherical RHF/UHF states use only \`mol.cart2sph\_coeff(normalized="sp")\` and \`C\_cart = X\*C\_sph\`, use \`mol.bas\_ctr\_coeff\` or its exact normalization-equivalent, validate the Cartesian source overlap/MO metric, and emit deterministic hashes. PySCF and NumPy remain external command dependencies. Add no PySCF calculation driver, basis-only or live-mean-field mode, global metric solve, Julia dependency, general format framework, basis lookup, GHF/spinor/complex/periodic/ECP/ROHF claim, Hamiltonian/Fock/ERI payload, or extra output file.

### HP-REP-XGTO-READER-DOC-01 - external Cartesian GTO reader manual

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `docs`
- **Execution whitelist:** `false`
- **Documents:**
  - `canonical` [external\_cartesian\_gto\_interchange.md](external_cartesian_gto_interchange.md); heading `Reader-Facing Manual Authority`
- **Owned paths:**
  - `docs` / `existing`: `docs/src/developer/designs/cartesian_hamiltonian_producer/external_cartesian_gto_interchange.md`
  - `docs` / `existing`: `docs/src/manual/external_cartesian_gto_transfer.md`
  - `docs` / `existing`: `docs/src/manual/index.md`
  - `docs` / `existing`: `docs/src/reference/export.md`
- **Evidence:**
  - `git_commit`: `a746038f0528c791314f273313e2f7142f3a03b0`
  - `manager_pass`: `517`
  - `manager_pass`: `518`
- **Dependencies:** `HP-REP-XGTO-CLOSESTDET-FN-01`, `HP-REP-XGTO-IMPORT-FN-01`, `HP-REP-XGTO-INTERCHANGE-FN-01`, `HP-REP-XGTO-PYSCF-EXPORT-FN-01`
- **Scope:** Maintain only the accepted reader-facing external Cartesian-GTO transfer manual, its Manual navigation and index link, and its curated seven-binding API section from a746038f0528c791314f273313e2f7142f3a03b0. Preserve checkpoint-only PySCF export, same-geometry working construction, final-by-source overlap inspection, unchanged raw projection, separately thresholded closest-determinant preparation, explicit capture diagnostics, Hamiltonian attestation, optional external PySCF/NumPy dependencies, and all unsupported-case boundaries. Add no README promotion, heavy example, source/API/fixture/format/workflow/dependency/version/tag/release/registration/citation/paper change, basis-only/live-mean-field exporter, universal capture threshold, solver state, private name, paper identifier, C2 production number, or duplicate gto\_overlap\_matrix reference entry.

### HP-REP-XGTO-READER-DOC-TEST-01 - external Cartesian GTO reader documentation validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [external\_cartesian\_gto\_interchange.md](external_cartesian_gto_interchange.md); heading `Reader-Facing Manual Authority`
- **Owned paths:**
  - `test` / `existing`: `test/docs/runtests.jl`
- **Evidence:**
  - `git_commit`: `a746038f0528c791314f273313e2f7142f3a03b0`
  - `manager_pass`: `517`
  - `manager_pass`: `518`
- **Dependencies:** `HP-REP-XGTO-READER-DOC-01`
- **Scope:** Maintain only the accepted focused documentation checks from a746038f0528c791314f273313e2f7142f3a03b0: manual navigation and index linkage, resolution and curated-reference presence of the seven owned bindings, raw projection versus determinant cleanup, checkpoint-only export, explicit minimum-Gram threshold, Hamiltonian attestation, and absence of personal/Dropbox paths, request identifiers, C2 production numbers, and private names. Add no test file, executable example, numerical assertion, PySCF installation, broad docs scan, runner/workflow change, or duplicate gto\_overlap\_matrix ownership.

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

### HP-RETIRE-CONTRACTED-PARENT-FN-01 - retire obsolete contracted-parent and multilayer route

- **Lifecycle:** `retired`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [contracted\_parent\_multilayer\_retirement.md](contracted_parent_multilayer_retirement.md); heading `Contracted-Parent And Multilayer Retirement`
- **Owned paths:** none
- **Evidence:**
  - `git_commit`: `7fbe512dc06df0be6c12a99f23aafe392c32e8a3`
  - `manager_pass`: `467`
- **Dependencies:** none
- **Scope:** Delete exactly the no-caller contracted-parent/metrics submodules, five-file PQS multilayer family, complete-core-shell final-basis utility, seven root includes, final-basis include/doc claim, private GTO source-box shadow, dead atomic multilayer adapter, and three tracked probes, for 9,541 source plus 983 probe lines. Add no replacement, alias, shim, helper, moved kernel, or compatibility surface; preserve every current parent, terminal, PQS/WL, PRF, pair-materialization, raw-block, represented-Hartree, IDA, artifact, and driver owner. Stop without committing if a live caller or side effect is found.

### HP-RETIRE-CONTRACTED-PARENT-TEST-01 - validate contracted-parent and multilayer retirement

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [contracted\_parent\_multilayer\_retirement.md](contracted_parent_multilayer_retirement.md); heading `Contracted-Parent And Multilayer Retirement`
- **Owned paths:** none
- **Evidence:**
  - `git_commit`: `7fbe512dc06df0be6c12a99f23aafe392c32e8a3`
  - `manager_pass`: `467`
- **Dependencies:** none
- **Scope:** Delete only the obsolete 66-line contracted-parent PGDG-factor testset from test/ida/runtests.jl, add or edit no replacement test, then validate package load, unchanged core/ida/cartesian groups, public Cartesian terminal due diligence, exact line accounting, the anti-bloat scan, and git diff --check.

### HP-RETIRE-CPB-PROVIDER-FN-01 - retire orphaned Cartesian CPB block-provider pilot

- **Lifecycle:** `retired`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [cartesian\_cpb\_block\_provider\_retirement.md](cartesian_cpb_block_provider_retirement.md); heading `Cartesian CPB Block Provider Retirement`
- **Owned paths:** none
- **Evidence:**
  - `git_commit`: `181ed6968dbe0845d1d8fb2fdca0c6597a96dec0`
  - `manager_pass`: `465`
- **Dependencies:** none
- **Scope:** Delete only the orphaned CartesianCPBBlockProviders module and its sole GaussletBases include, retiring its 47-name qualified internal surface without replacement, compatibility glue, aliases, helpers, tests, or changes to preserved CPB, parent, raw-block, terminal, pair-materialization, residual, or represented-Hartree owners.

### HP-RETIRE-CPB-PROVIDER-TEST-01 - validate Cartesian CPB block-provider retirement

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [cartesian\_cpb\_block\_provider\_retirement.md](cartesian_cpb_block_provider_retirement.md); heading `Cartesian CPB Block Provider Retirement`
- **Owned paths:** none
- **Evidence:**
  - `git_commit`: `181ed6968dbe0845d1d8fb2fdca0c6597a96dec0`
  - `manager_pass`: `465`
- **Dependencies:** none
- **Scope:** Run the existing core, ida, and Cartesian endpoint gates unchanged after exact provider/include deletion; add or edit no test, and stop for a separate amendment if validation would require test changes.

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

### HP-RETIRE-PAIR-LADDER-FN-01 - retire orphaned Cartesian pair planning and materialization ladder

- **Lifecycle:** `retired`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [cartesian\_pair\_planning\_materialization\_retirement.md](cartesian_pair_planning_materialization_retirement.md); heading `Cartesian Pair Planning And Materialization Retirement`
- **Owned paths:** none
- **Evidence:**
  - `git_commit`: `32bb3e2a7c58cf935353e34217f4226f84a6557a`
  - `manager_pass`: `470`
- **Dependencies:** none
- **Scope:** Delete the orphaned unit-pair, pair-operator-plan, pair-block materialization, source-shell final-basis, standalone module-wrapper, and ladder-tool surfaces exactly as specified, while retaining pqs\_source\_axis\_transforms.jl in place behind a minimal CPBM owner. Add no alias, shim, replacement framework, helper, test, metadata, or status vocabulary; require at least approximately 9,900 net source lines removed and stop if a live caller, side effect, or mapped-COMX output change is found.

### HP-RETIRE-PAIR-LADDER-TEST-01 - validate Cartesian pair-ladder retirement

- **Lifecycle:** `completed`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `evidence` [cartesian\_pair\_planning\_materialization\_retirement.md](cartesian_pair_planning_materialization_retirement.md); heading `Cartesian Pair Planning And Materialization Retirement`
- **Owned paths:** none
- **Evidence:**
  - `git_commit`: `32bb3e2a7c58cf935353e34217f4226f84a6557a`
  - `manager_pass`: `470`
- **Dependencies:** none
- **Scope:** Run the existing core, ida, and Cartesian endpoint tests unchanged plus one transient before/after mapped-COMX terminal parity construction; add or edit no committed test or probe, and stop for a separate amendment if validation requires test changes.

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

### HP-RG-PROTECT-EGOI-FN-01 - deferred retained-GTO local-product EGOI helper

- **Lifecycle:** `deferred`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [retained\_gto\_egoi.md](retained_gto_egoi.md); heading `Retained-GTO Local-Product EGOI`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** deferred with no source authority. The unaccepted retained-original \`s1+s2\` / local-product / \`M2\` adapter is preserved only on archive branch \`archive/retained-gto-egoi-wip-2026-08-03\`; any implementation requires a new docs-only reactivation tied to a current physics target.

### HP-RG-PROTECT-EGOI-TEST-01 - deferred retained-GTO EGOI validation

- **Lifecycle:** `deferred`
- **Grant:** `none`
- **Surfaces:** none
- **Execution whitelist:** `false`
- **Documents:**
  - `history` [retained\_gto\_egoi.md](retained_gto_egoi.md); heading `Retained-GTO Local-Product EGOI`
- **Owned paths:** none
- **Evidence:** none
- **Dependencies:** none
- **Scope:** deferred with no test authority. The planned focused test was never created; future validation requires the same docs-only reactivation as the specialized helper.

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

### HP-SLICE-HCHAIN-FN-01 - minimal sliced hydrogen-chain producer

- **Lifecycle:** `implemented`
- **Grant:** `maintenance`
- **Surfaces:** `source`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [sliced\_hydrogen\_chain.md](sliced_hydrogen_chain.md); heading `Sliced Hydrogen-Chain Producer`
- **Owned paths:**
  - `source` / `existing`: `src/GaussletBases.jl`
  - `source` / `existing`: `src/sliced_hydrogen_chain.jl`
- **Evidence:**
  - `git_commit`: `c895e353d`
  - `git_commit`: `47023f190`
  - `git_commit`: `6bbb01ade`
  - `manager_pass`: `455`
  - `manager_pass`: `457`
- **Dependencies:** none
- **Scope:** Maintain the accepted six-export expert interface and compact numerical owner. Preserve sliced\_h1\_bandwidth(chain)::Int as the O(1), zero-allocation query over existing construction state for the represented H1 structural half-bandwidth, including finite, periodic-template, and singleton behavior and exact represented zeros outside the reported band. Preserve bounded G10 donor parity, deterministic H/STO-6G transverse construction, full analytic Galerkin H1, longitudinal IntegralDiagonal Vee, common high135 near evaluation with checked long-range continuations, bare-operator sector independence, and the existing stable consumer fields. Add no further export, field, threshold, dense matrix, framework, driver, artifact, solver, MPO, multi-transverse support, Cartesian/PQS coupling, or publication policy under this grant.

### HP-SLICE-HCHAIN-TEST-01 - sliced hydrogen-chain producer validation

- **Lifecycle:** `completed`
- **Grant:** `maintenance`
- **Surfaces:** `tests`
- **Execution whitelist:** `true`
- **Documents:**
  - `canonical` [sliced\_hydrogen\_chain.md](sliced_hydrogen_chain.md); heading `Sliced Hydrogen-Chain Producer`
- **Owned paths:**
  - `test` / `existing`: `test/ida/runtests.jl`
- **Evidence:**
  - `git_commit`: `c895e353d`
  - `git_commit`: `47023f190`
  - `git_commit`: `6bbb01ade`
  - `manager_pass`: `455`
  - `manager_pass`: `457`
- **Dependencies:** `HP-SLICE-HCHAIN-FN-01`
- **Scope:** Maintain the committed IDA-owner checks, including finite, periodic-template, and singleton structural H1-bandwidth bounds, exact represented-operator zeros outside the reported band, and zero-allocation constant-time query access. Preserve deterministic H/STO-6G transverse construction, bounded G10 donor parity, analytic H1 and direct-quadrature Vee oracles, scalar/row/action parity, finite tail diagnostics, malformed-input rejection, the H10 dense fixture, and H1000 compact construction. H10000 remains transient evidence. Add no HF, DMRG, MPO, publication-energy, fixture file, or test owner under this grant.

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
- **Scope:** Keep odd-side enforcement only for direct nucleus-centered cores. Cubic WL boundary products retain the requested route-local count, with canonical examples \`ns=4 -\> 56\` and \`ns=5 -\> 98\`; eligible matched aspect shells are separately governed by \`HP-PQS-ASPECTSHELL-FN-01\` and use axis-specific inner counts from the shared outer shape.

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
