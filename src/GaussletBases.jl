module GaussletBases

using JLD2
using LinearAlgebra
using SHA
using SparseArrays
using SpecialFunctions
using TOML
const _PACKAGE_ROOT = normpath(joinpath(@__DIR__, ".."))
_package_data_path(parts...) = normpath(joinpath(_PACKAGE_ROOT, "data", parts...))

export AbstractFunction1D,
       AbstractPrimitiveFunction1D,
       AbstractBasisFunction1D,
       AbstractCoordinateMapping,
       AbstractBasisSpec,
       AbstractDiagonalApproximation,
       PrimitiveSet1D,
       BasisMetadata1D,
       BasisRepresentation1D,
       CartesianBasisMetadata3D,
       CartesianBasisRepresentation3D,
       CartesianGaussianShellOrbitalRepresentation3D,
       CartesianGaussianShellSupplementRepresentation3D,
       CartesianIDAHamiltonian,
       PQSH2PlusRow,
       PQSH2PlusComparison,
       pqs_h2plus_comparison,
       ExactRepresentedHartreeField,
       FittedReferenceHartreeField,
       ScreenedHartreeCorrection,
       screened_hartree_correction,
       screened_hartree_delta_one_body,
       screened_hartree_energy_constant,
       screened_hartree_consistency_error,
       screened_hartree_field_kind,
       cartesian_base_hamiltonian,
       cartesian_base_working_basis,
       cartesian_residual_gto_mwg_system,
       CartesianBasisTransferDiagnostics,
       CartesianBasisProjector3D,
       CartesianOrbitalTransferResult,
       @timeg,
       timing_enabled,
       timing_live_enabled,
       set_timing!,
       set_timing_live!,
       set_timing_thresholds!,
       reset_timing_report!,
       current_timing_report,
       timing_report,
       cross_overlap,
       gto_overlap_matrix,
       gto_occupancy_matrix,
       ExternalGTOOrbitalSpinBlock,
       ExternalGTOOrbitalPacket,
       ExternalGTOOrbitalSpinImport,
       ExternalGTOOrbitalImportResult,
       external_gto_overlap_fingerprint,
       external_gto_ordering_fingerprint,
       import_external_gto_orbitals,
       read_external_cartesian_gto_packet,
       closest_external_gto_determinant,
       basis_projector,
       transfer_orbitals,
       gaussian_coulomb_pair_index,
       gaussian_coulomb_pair_matrix,
       EGOIDensityDensityCorrectionResult,
       StationaryFockCorrectionResult,
       HamiltonianCorrectionResult,
       OrdinaryProjectedHamiltonianCorrectionTarget,
       egoi_target_product_matrix,
       egoi_target_coulomb_matrix,
       egoi_density_density_correction,
       projected_orbital_density,
       density_density_restricted_fock,
       occupied_virtual_fock_residual,
       stationary_fock_one_body_correction,
       egoi_stationary_hamiltonian_correction,
       ordinary_cartesian_projected_gaussian_target,
       ordinary_cartesian_egoi_stationary_correction,
       CartesianGTOSubspaceProjectionResult,
       RadialYlmSolidHarmonicGTOFit,
       RadialYlmCartesianGTOAdapter,
       RadialYlmCartesianProjectionResult,
       fit_radial_ylm_to_solid_harmonic_gto,
       evaluate_radial_ylm_gto_fit,
       radial_ylm_fit_cartesian_gto_adapter,
       project_radial_ylm_gto_adapter_to_cartesian,
       project_cartesian_gto_to_supplement_subspace,
       BasisBox1D,
       BasisPartition1D,
       HierarchicalBasisBox1D,
       HierarchicalBasisPartition1D,
       LeafLocalPGDG1D,
       LeafGaussianSpec1D,
       GlobalMappedPrimitiveLayer1D,
       LeafBoxContraction1D,
       LeafBoxContractionLayer1D,
       UniformBasisSpec,
       MappedUniformBasisSpec,
       LegacyAtomicGaussianShell,
       LegacyAtomicGaussianSupplement,
       LegacyBondAlignedDiatomicGaussianSupplement,
       LegacyBondAlignedHeteronuclearGaussianSupplement,
       LegacySGaussianData,
       BondAlignedDiatomicQWBasis3D,
       BondAlignedHomonuclearChainQWBasis3D,
       AxisAlignedHomonuclearSquareLatticeQWBasis3D,
       BondAlignedDiatomicGeometryPoint3D,
       BondAlignedDiatomicGeometryNucleus3D,
       BondAlignedDiatomicGeometryBox3D,
       BondAlignedDiatomicGeometryPayload3D,
       BondAlignedDiatomicGeometryPlaneSlice3D,
       MappedOrdinaryOneBody1D,
       CartesianProductOrbital3D,
       OrdinaryCartesianOrbital3D,
       QiuWhiteHybridOrbital3D,
       QWRGResidualSpaceDiagnostics,
       OrdinaryCartesianIDAOperators,
       OrdinaryCartesianOperators3D,
       HalfLineBasisSpec,
       RadialBasisSpec,
       RadialBoundaryPrototype,
       Gaussian,
       HalfLineGaussian,
       XGaussian,
       recommended_xgaussians,
       Distorted,
       GaussletFamily,
       Gausslet,
       UniformBasis,
       MappedUniformBasis,
       MappedPGDGPrototype1D,
       MappedPGDGLocalized1D,
       HalfLineBasis,
       RadialBasis,
       MappedGausslet,
       BoundaryGausslet,
       RadialGausslet,
       StencilTerm,
       FunctionStencil,
       value,
       direct_value,
       derivative,
       center,
       reference_center,
       moment_center,
       integral_weight,
       stencil,
       stencil_matrix,
       build_basis,
       radial_boundary_prototype_names,
       radial_boundary_prototype,
       build_paper_parity_radial_basis,
       basis_metadata,
       basis_representation,
       basis_partition,
       hierarchical_partition,
       build_leaf_pgdg,
       augment_leaf_pgdg,
       build_global_mapped_primitive_layer,
       contract_leaf_boxes,
       refine_partition,
        primitive_set,
        boxes,
        leaf_boxes,
        leaf_contractions,
        leaf_primitive_indices,
        primitive_origins,
        primitive_leaf_boxes,
        box_indices,
        box_level,
        box_parent,
       box_children,
       box_block,
       box_coupling,
       basis_spec,
       family,
       mapping,
       centers,
       reference_centers,
       integral_weights,
       contract_primitive_vector,
       contract_primitive_diagonal,
       contract_primitive_matrix,
       CoulombGaussianExpansion,
       coulomb_gaussian_expansion,
       SlicedHydrogenChain,
       sliced_hydrogen_chain,
       sliced_h1,
       sliced_h1_bandwidth,
       sliced_vee,
       sliced_row!,
       mapped_pgdg_prototype,
       mapped_pgdg_localized,
       legacy_atomic_gaussian_supplement,
       legacy_bond_aligned_diatomic_gaussian_supplement,
       legacy_bond_aligned_heteronuclear_gaussian_supplement,
       legacy_s_gaussian_data,
       bond_aligned_homonuclear_qw_basis,
       bond_aligned_homonuclear_chain_qw_basis,
       axis_aligned_homonuclear_square_lattice_qw_basis,
       bond_aligned_heteronuclear_qw_basis,
       bond_aligned_diatomic_nested_fixed_source,
       bond_aligned_diatomic_nested_fixed_block,
       bond_aligned_diatomic_nested_geometry_diagnostics,
       bond_aligned_homonuclear_chain_geometry_diagnostics,
       axis_aligned_homonuclear_square_lattice_geometry_diagnostics,
       build_one_center_atomic_full_parent_shell_sequence,
       one_center_atomic_full_parent_fixed_block,
       build_one_center_atomic_legacy_profile_shell_sequence,
       one_center_atomic_legacy_profile_fixed_block,
       OneCenterAtomicNestedStructureDiagnostics,
       one_center_atomic_nested_structure_diagnostics,
       one_center_atomic_nested_structure_report,
       bond_aligned_diatomic_geometry_payload,
       bond_aligned_diatomic_source_geometry_payload,
       bond_aligned_diatomic_plane_slice,
       SpherePointSetProvenance,
       SpherePointSet,
       CuratedSpherePointSet,
       ShellLocalInjectedAngularBasis,
       ShellLocalAngularProfileKey,
       ShellLocalAngularProfile,
       ShellLocalAngularProfileOverlap,
       AtomicFixedRadialAngularSequenceLevel,
       AtomicFixedRadialAngularSequenceOverlapSidecar,
       AtomicFixedRadialAngularSequence,
       AtomicShellLocalInjectedAngularAssembly,
       AtomicInjectedAngularOneBodyBenchmark,
       AtomicInjectedAngularCartesianMomentBundle,
       AtomicInjectedAngularHFStyleBenchmark,
       AtomicInjectedAngularHFDMRGHFAdapter,
       AtomicInjectedAngularSmallEDBenchmark,
       fibonacci_sphere_point_set,
       optimize_sphere_point_set,
       sphere_point_set,
       sphere_point_set_orders,
       curated_sphere_point_set,
       curated_sphere_point_set_orders,
       shell_local_angular_profile,
       adjacent_shell_local_angular_profile_overlap,
       build_atomic_fixed_radial_angular_sequence,
       atomic_fixed_radial_angular_level_dense_payload,
       atomic_fixed_radial_angular_overlap_sidecar_payload,
       atomic_fixed_radial_legacy_dmrgatom_payload,
       build_shell_local_injected_angular_basis,
       assign_atomic_angular_shell_orders,
       build_atomic_shell_local_angular_assembly,
       build_atomic_injected_angular_one_body_benchmark,
       build_atomic_injected_angular_cartesian_moments,
       build_atomic_injected_angular_hf_style_benchmark,
       build_atomic_injected_angular_hfdmrg_payload,
       build_atomic_injected_angular_hfdmrg_hf_adapter,
       build_atomic_injected_angular_hfdmrg_hf_seeds,
       build_atomic_injected_angular_small_ed_benchmark,
       shell_local_injected_angular_diagnostics,
       atomic_shell_local_angular_diagnostics,
       atomic_injected_angular_one_body_diagnostics,
       atomic_injected_angular_hf_style_diagnostics,
       atomic_injected_angular_hfdmrg_hf_adapter_diagnostics,
       atomic_injected_angular_small_ed_diagnostics,
       run_atomic_injected_angular_hfdmrg_hf,
       mapped_ordinary_one_body_operators,
       mapped_cartesian_hydrogen_energy,
       ordinary_sho_hamiltonian,
       ordinary_sho_spectrum,
       ordinary_cartesian_ida_operators,
       ordinary_cartesian_product_operators,
       ordinary_cartesian_qiu_white_operators,
       HydrogenicCoreProjectorCorrectionSpec,
       HydrogenicCoreBranchCorrectionSpec,
       OrdinaryCartesianCorrectionResult,
       OrdinaryCartesianBranchCorrectionResult,
       apply_ordinary_cartesian_corrections,
       ordinary_cartesian_corrected_branch,
       assembled_one_body_hamiltonian,
       nested_cartesian_operators,
       diagnose_qwrg_residual_space,
       ordinary_cartesian_vee_expectation,
       angular_benchmark_exact_hamv6_payload,
       atomic_hamv6_payload,
       atomic_ida_density_interaction_matrix,
       fullida_dense_payload,
       sliced_ham_payload,
       write_angular_benchmark_exact_hamv6_jld2,
       write_atomic_fixed_radial_angular_level_jld2,
       write_atomic_fixed_radial_angular_overlap_sidecar_jld2,
       write_atomic_fixed_radial_legacy_dmrgatom_jld2,
       write_atomic_hamv6_jld2,
       write_cartesian_ida_hamiltonian,
       write_fullida_dense_jld2,
       write_sliced_ham_jld2,
       read_cartesian_ida_hamiltonian,
       one_body_hamiltonian,
       nuclear_repulsion,
       gaussian_factor_matrix,
       gaussian_factor_matrices,
       RadialQuadratureGrid,
       radial_quadrature,
       quadrature_points,
       quadrature_weights,
       basis_diagnostics,
       IntegralDiagonal,
       overlap_matrix,
       position_matrix,
       kinetic_matrix,
       nuclear_matrix,
       centrifugal_matrix,
       multipole_matrix,
       YlmChannel,
       YlmChannelSet,
       AtomicOrbital,
       RadialAtomicOperators,
       AtomicOneBodyOperators,
       AtomicIDAOperators,
       AtomicIDATwoElectronState,
       AtomicIDATwoElectronProblem,
       ylm_channels,
       atomic_one_body_operators,
       atomic_ida_operators,
       atomic_ida_two_electron_problem,
       channel_range,
       channel_hamiltonian,
       channel_overlap,
       orbitals,
       two_electron_states,
       radial_multipole,
       direct_matrix,
       exchange_matrix,
       fock_matrix,
       fock_matrix_alpha,
       fock_matrix_beta,
       density_matrix,
       uhf_energy,
       uhf_step,
       uhf_scf,
       gaunt_tensor,
       gaunt_coefficient,
       angular_kernel,
       apply_overlap,
       apply_hamiltonian,
       ground_state_energy,
       lanczos_ground_state,
       atomic_operators,
       centrifugal,
       multipole,
       coefficients,
       primitives,
       terms,
       IdentityMapping,
       AsinhMapping,
       white_lindsey_atomic_mapping,
       CombinedInvsqrtMapping,
       fit_asinh_mapping_for_extent,
       fit_asinh_mapping_for_strength,
       fit_combined_invsqrt_mapping,
       uofx,
       xofu,
       dudx,
       du2dx2

"""
    AbstractFunction1D

Abstract supertype for callable one-dimensional function objects.
"""
abstract type AbstractFunction1D end

"""
    AbstractPrimitiveFunction1D <: AbstractFunction1D

Abstract supertype for lowest-level primitive function objects.
"""
abstract type AbstractPrimitiveFunction1D <: AbstractFunction1D end

"""
    AbstractBasisFunction1D <: AbstractFunction1D

Abstract supertype for higher-level callable basis functions.
"""
abstract type AbstractBasisFunction1D <: AbstractFunction1D end

"""
    AbstractCoordinateMapping

Abstract supertype for coordinate maps between physical `x` and reference `u`.
"""
abstract type AbstractCoordinateMapping end

"""
    AbstractBasisSpec

Abstract supertype for public basis-construction recipes.
"""
abstract type AbstractBasisSpec end

"""
    AbstractDiagonalApproximation

Abstract supertype for supported diagonal-approximation choices in the radial
electron-electron operator layer.
"""
abstract type AbstractDiagonalApproximation end

"""
    value(f, x)

Evaluate `f` at physical coordinate `x` through its normal evaluation route.
Calling an `AbstractFunction1D` as `f(x)` delegates here.
"""
function value end

"""
    direct_value(f, x)

Evaluate `f` at physical `x` directly from `stencil(f)`. This is exact for
that representation, not an independent analytic oracle.
"""
function direct_value end

"""
    derivative(f, x; order=1)

Evaluate the requested nonnegative derivative order of `f` with respect to
physical coordinate `x`.
"""
function derivative end

"""
    center(f)

Return the physical-coordinate center of `f`.
"""
function center end

"""
    reference_center(f)

Return the pre-mapping, reference-coordinate center corresponding to `f`.
"""
function reference_center end
function moment_center end

"""
    integral_weight(f)

Return the integral of one function `f`. This is distinct from basis-level
`integral_weights` and from quadrature weights.
"""
function integral_weight end

"""
    stencil(f)

Return the ordered primitive expansion representing `f` as a `FunctionStencil`.
"""
function stencil end
function stencil_matrix end
function build_basis end

"""
    basis_metadata(object)

Return metadata associated with a supported basis, layer, representation, or
supplement. Its concrete type, fields, ownership, and copy behavior are defined
by the owning object; no universal metadata schema is promised.
"""
function basis_metadata end
function basis_representation end
function cross_overlap end
function gto_overlap_matrix end
function gto_occupancy_matrix end
function basis_projector end
function transfer_orbitals end
function basis_partition end
function hierarchical_partition end
function build_leaf_pgdg end
function augment_leaf_pgdg end
function build_global_mapped_primitive_layer end
function contract_leaf_boxes end
function refine_partition end
function primitive_set end

"""
    boxes(partition)

Return stored box records in partition order; treat the vector and records as read-only.
"""
function boxes end

"""
    leaf_boxes(hierarchy)

Select hierarchy nodes with no children, preserving stored order. Treat the
returned records and their contained vectors as read-only.
"""
function leaf_boxes end

"""
    leaf_contractions(layer)

Return stored `LeafBoxContraction1D` records in order; treat them as read-only.
"""
function leaf_contractions end

"""
    leaf_primitive_indices(generator, leaf_box_index)

Return the global primitive-index range owned by the selected leaf box.
"""
function leaf_primitive_indices end

"""
    primitive_origins(generator)

Return the stored primitive-aligned provenance labels. Treat the returned
vector as read-only.
"""
function primitive_origins end

"""
    primitive_leaf_boxes(generator)

Return the stored primitive-aligned leaf-box provenance. Treat the returned
vector as read-only.
"""
function primitive_leaf_boxes end

"""
    box_indices(partition, box_index)

Return the selected box's stored basis-index vector; treat it as read-only.
"""
function box_indices end

"""
    box_level(hierarchy, box_index)

Return the selected hierarchy node's level.
"""
function box_level end

"""
    box_parent(hierarchy, box_index)

Return the parent identifier, or `nothing` when the selected node is a root.
"""
function box_parent end

"""
    box_children(hierarchy, box_index)

Return the selected node's stored child identifiers. Treat the returned vector
as read-only.
"""
function box_children end

"""
    box_block(matrix, partition, box_index)
    box_block(representation, partition, operator_name, box_index)

Materialize a fresh `Matrix{Float64}` copy of one diagonal box block. A
`BasisRepresentation1D` call selects the named operator.
"""
function box_block end

"""
    box_coupling(matrix, partition, box_i, box_j)
    box_coupling(representation, partition, operator_name, box_i, box_j)

Materialize a fresh `Matrix{Float64}` copy of the rectangular block between
two boxes. A `BasisRepresentation1D` call selects the named operator.
"""
function box_coupling end

"""
    basis_spec(object)

Return the construction specification associated with a supported basis-like object.
"""
function basis_spec end

"""
    family(object)

Identify the gausslet family associated with a supported basis-like object.
"""
function family end

"""
    mapping(object)

Return the map connecting the object's reference and physical coordinates.
"""
function mapping end

"""
    centers(object)

Return the object's basis centers in physical coordinates.
"""
function centers end

"""
    reference_centers(object)

Return the object's basis centers in reference coordinates.
"""
function reference_centers end

"""
    integral_weights(object)

Return the integrals of the object's basis functions. These are not quadrature weights.
"""
function integral_weights end
function contract_primitive_vector end
function contract_primitive_diagonal end
function contract_primitive_matrix end
function coulomb_gaussian_expansion end
function mapped_pgdg_prototype end
function mapped_pgdg_localized end
function hybrid_mapped_ordinary_basis end
function bond_aligned_diatomic_nested_fixed_source end
function bond_aligned_diatomic_nested_fixed_block end
function bond_aligned_diatomic_nested_geometry_diagnostics end
function bond_aligned_homonuclear_chain_geometry_diagnostics end
function axis_aligned_homonuclear_square_lattice_qw_basis end
function axis_aligned_homonuclear_square_lattice_geometry_diagnostics end
function build_one_center_atomic_full_parent_shell_sequence end
function one_center_atomic_full_parent_fixed_block end
function build_one_center_atomic_legacy_profile_shell_sequence end
function one_center_atomic_legacy_profile_fixed_block end
function one_center_atomic_nested_structure_diagnostics end
function one_center_atomic_nested_structure_report end
"""
    bond_aligned_diatomic_geometry_payload(object...)

Construct bond-aligned geometry records for inspection and visualization.

Overloads cover supported bond-aligned QW bases, ordinary operators, nested
sources and fixed blocks, and residual-augmented operators. Multi-object
overloads require the same basis geometry where the corresponding method
enforces it. Nested payloads contain compressed fixed centers; use
[`bond_aligned_diatomic_source_geometry_payload`](@ref) to inspect raw
source-region points separately.

This operation neither constructs nor mutates a basis or Hamiltonian. Treat
the returned records and their vectors as read-only inspection data. The
interface follows the current bond-aligned QW/nested geometry and provenance
conventions; it is not a general molecular-geometry API, plotting framework,
or stable serialization format.
"""
function bond_aligned_diatomic_geometry_payload end
function bond_aligned_diatomic_source_geometry_payload end
function bond_aligned_diatomic_plane_slice end
function _qwrg_bond_aligned_axis_bundles end
function mapped_ordinary_one_body_operators end
function mapped_cartesian_hydrogen_energy end
function ordinary_sho_hamiltonian end
function ordinary_sho_spectrum end
function ordinary_cartesian_ida_operators end
function ordinary_cartesian_product_operators end
function ordinary_cartesian_qiu_white_operators end
function apply_ordinary_cartesian_corrections end
function ordinary_cartesian_corrected_branch end
function assembled_one_body_hamiltonian end
function diagnose_qwrg_residual_space end
function ordinary_cartesian_1s2_check end
function ordinary_cartesian_vee_expectation end
function nested_cartesian_operators end
function angular_benchmark_exact_hamv6_payload end
function build_atomic_fixed_radial_angular_sequence end
function atomic_fixed_radial_angular_level_dense_payload end
function atomic_fixed_radial_angular_overlap_sidecar_payload end
function build_atomic_injected_angular_cartesian_moments end
function fullida_dense_payload end
function sliced_ham_payload end
function write_angular_benchmark_exact_hamv6_jld2 end
function write_atomic_fixed_radial_angular_level_jld2 end
function write_atomic_fixed_radial_angular_overlap_sidecar_jld2 end
function write_fullida_dense_jld2 end
function write_sliced_ham_jld2 end
function gaussian_factor_matrix end
function gaussian_factor_matrices end
function egoi_target_product_matrix end
function egoi_target_coulomb_matrix end
function egoi_density_density_correction end
function projected_orbital_density end
function density_density_restricted_fock end
function occupied_virtual_fock_residual end
function stationary_fock_one_body_correction end
function egoi_stationary_hamiltonian_correction end
function ordinary_cartesian_projected_gaussian_target end
function ordinary_cartesian_egoi_stationary_correction end
function fit_radial_ylm_to_solid_harmonic_gto end
function evaluate_radial_ylm_gto_fit end
function radial_ylm_fit_cartesian_gto_adapter end
function project_radial_ylm_gto_adapter_to_cartesian end
function project_cartesian_gto_to_supplement_subspace end
function radial_quadrature end
function quadrature_points end
function quadrature_weights end
function basis_diagnostics end
function overlap_matrix end
function position_matrix end
function kinetic_matrix end
function nuclear_matrix end
function centrifugal_matrix end
function multipole_matrix end
function ylm_channels end
function atomic_one_body_operators end
function atomic_ida_operators end
function atomic_ida_two_electron_problem end
function channel_range end
function channel_hamiltonian end
function channel_overlap end

"""
    orbitals(owner)

Return the owner's stored orbital records in its construction-defined order.
Treat the returned storage and records as read-only; no universal ordering is implied.
"""
function orbitals end

"""
    two_electron_states(problem)

Return the tiny problem's stored state records in its defined order. Treat the
returned storage and records as read-only.
"""
function two_electron_states end

"""
    radial_multipole(operators, L)

Return the requested stored radial multipole matrix from `AtomicIDAOperators`.
Treat the matrix as read-only.
"""
function radial_multipole end
function direct_matrix end
function exchange_matrix end
function fock_matrix end
function fock_matrix_alpha end
function fock_matrix_beta end
function density_matrix end
function uhf_energy end
function uhf_step end
function uhf_scf end

"""
    gaunt_tensor(operators, L)

Reconstruct the requested dense Gaunt tensor on demand. The result is not a
retained cache or a scalable large-angular-momentum representation.
"""
function gaunt_tensor end

"""
    gaunt_coefficient(operators, L, M, left, right)

Return one complex spherical-harmonic Gaunt coefficient.
"""
function gaunt_coefficient end

"""
    angular_kernel(operators, L)

Reconstruct the requested dense angular kernel on demand. The result is not a
retained cache or a scalable large-angular-momentum representation.
"""
function angular_kernel end

"""
    apply_overlap(problem, coefficients)

Validate the coefficient dimension and apply the tiny problem's stored dense
overlap matrix.
"""
function apply_overlap end

"""
    apply_hamiltonian(problem, coefficients)

Validate the coefficient dimension and apply the tiny problem's stored dense
Hamiltonian matrix.
"""
function apply_hamiltonian end

"""
    ground_state_energy(problem)

Return the lowest eigenvalue of the tiny problem's stored dense Hermitian
Hamiltonian.
"""
function ground_state_energy end

"""
    lanczos_ground_state(problem; krylovdim=200, maxiter=200, tol=1e-10, v0=nothing)

Run the small, fully reorthogonalized reference Lanczos routine. It returns
`value`, `vector`, `residual`, `iterations`, and `converged` and is not a
general many-electron or production eigensolver.
"""
function lanczos_ground_state end
function atomic_operators end
function centrifugal end
function multipole end

"""
    coefficients(object)

Return stored coefficient data. For a `FunctionStencil`, its order matches
`primitives(stencil)`.
"""
function coefficients end
function primitives end

"""
    terms(stencil)

Materialize the ordered coefficient/primitive pairs as `StencilTerm` values.
"""
function terms end

"""
    uofx(mapping, x)

Map physical coordinate `x` to reference coordinate `u`.
"""
function uofx end

"""
    xofu(mapping, u)

Map reference coordinate `u` to physical coordinate `x`, the inverse of `uofx`.
"""
function xofu end

"""
    dudx(mapping, x)

Evaluate the first derivative `du/dx` of the reference map `u(x)` at physical `x`.
"""
function dudx end

"""
    du2dx2(mapping, x)

Evaluate the second derivative `d^2u/dx^2` of the reference map `u(x)` at physical `x`.
"""
function du2dx2 end
function fit_asinh_mapping_for_extent end

(f::AbstractFunction1D)(x::Real) = value(f, x)
(mapping::AbstractCoordinateMapping)(x::Real) = uofx(mapping, x)

value(f::AbstractFunction1D, x::Real) = direct_value(f, x)
direct_value(f::AbstractFunction1D, x::Real) = stencil(f)(x)

function derivative(f::AbstractFunction1D, x::Real; order::Int = 1)
    order >= 0 || throw(ArgumentError("derivative order must be nonnegative"))
    order == 0 && return value(f, x)
    st = stencil(f)
    total = 0.0
    for i in eachindex(coefficients(st))
        total += coefficients(st)[i] * derivative(primitives(st)[i], x; order = order)
    end
    return total
end

reference_center(f::AbstractFunction1D) = center(f)

function integral_weight(f::AbstractFunction1D)
    st = stencil(f)
    total = 0.0
    for term in terms(st)
        total += term.coefficient * integral_weight(term.primitive)
    end
    return total
end

include("mappings.jl")
include("stencils.jl")
include("functions.jl")
include("GaussianAnalyticIntegrals.jl")
include("timing.jl")
using .TimeG: @timeg,
             timing_enabled,
             timing_live_enabled,
             set_timing!,
             set_timing_live!,
             set_timing_thresholds!,
             reset_timing_report!,
             current_timing_report,
             timing_report
include("internal/wavelet_filters.jl")
include("families.jl")
include("bases.jl")
include("quadrature.jl")
include("primitive_sets.jl")
include("ordinary_coulomb.jl")
include("ordinary_pgdg.jl")
include("ordinary_pgdg_refinement_masks.jl")
include("ordinary_hybrid.jl")
include("legacy_basis_adapter.jl")
include("ordinary_mapped_backends.jl")
include("ordinary_sho.jl")
include("ordinary_cartesian_ida.jl")
include("cartesian_nested_faces.jl")
include("cartesian_nested_owned_units.jl")
include("cartesian_nested_atomic.jl")
include("cartesian_nested_diatomic.jl")
include("cartesian_cpb/CartesianCPB.jl")
include("cartesian_raw_product_sources/CartesianRawProductSources.jl")
include("cartesian_residual_gaussians/CartesianResidualGaussians.jl")
include("cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl")
include("cartesian_final_basis_realization/CartesianFinalBasisRealization.jl")
include("cartesian_route_core/CartesianRouteCore.jl")
include("cartesian_terminal_shellification_geometry.jl")
include("cartesian_shellification/CartesianShellification.jl")
include("cartesian_terminal_lowering/CartesianTerminalLowering.jl")
include("cartesian_retained_units/CartesianRetainedUnits.jl")
include("cartesian_retained_unit_transform_contracts/CartesianRetainedUnitTransformContracts.jl")
include("cartesian_pair_block_materialization/CartesianPairBlockMaterialization.jl")
include("cartesian_nested_experimental_geometries.jl")
include("cartesian_gaussian_axis_integrals.jl")
include("cartesian_ida_hamiltonian.jl")
include("ordinary_qw_types_and_bases.jl")
include("CartesianParentGaussletBases.jl")
include("pqs_source_box_route_driver_helpers.jl")
include("pqs_source_box_low_order_materialization.jl")
include("cartesian_base_hamiltonian.jl")
include("pqs_matched_h2plus.jl")
include("cartesian_reference_density/CartesianReferenceDensity.jl")
using .CartesianReferenceDensity: ExactRepresentedHartreeField,
                                  FittedReferenceHartreeField,
                                  ScreenedHartreeCorrection,
                                  screened_hartree_correction,
                                  screened_hartree_delta_one_body,
                                  screened_hartree_energy_constant,
                                  screened_hartree_consistency_error,
                                  screened_hartree_field_kind
include("pqs_source_box_diatomic_complete_core_shell.jl")
include("pqs_source_box_route_driver_skeletons.jl")
include("ordinary_qw_nested_frontends.jl")
include("ordinary_qw_residuals.jl")
include("ordinary_qw_raw_blocks.jl")
include("ordinary_qw_operator_assembly.jl")
include("ordinary_qw_corrections.jl")
include("cartesian_basis_representation.jl")
include("cartesian_cross_overlap.jl")
include("cartesian_representation_constructors.jl")
include("cartesian_representation_transfer.jl")
include("cartesian_protected_ladder_bundle.jl")
include("cartesian_gto_probes.jl")
include("cartesian_external_gto_import.jl")
include("cartesian_external_gto_interchange.jl")
include("cartesian_qw_hybrid_representation.jl")
include("gaussian_coulomb_reference.jl")
include("hamiltonian_corrections.jl")
include("bond_aligned_diatomic_geometry.jl")
include("partitions.jl")
include("hierarchical_partitions.jl")
include("leaf_pgdg.jl")
include("global_leaf_contraction.jl")
include("diagnostics.jl")
include("radial_boundary_prototypes.jl")
include("operators.jl")
include("sliced_hydrogen_chain.jl")
include("atomic_ylm.jl")
include("radial_ylm_gto_bridge.jl")
include("gaunt_tables.jl")
include("angular_point_sets.jl")
include("angular_shell_basis.jl")
include("angular_shell_assembly.jl")
include("angular_atomic_benchmark.jl")
include("atomic_angular_sectors.jl")
include("atomic_ida.jl")
include("atomic_ida_direct.jl")
include("atomic_ida_exchange.jl")
include("atomic_ida_fock.jl")
include("atomic_ida_uhf.jl")
include("atomic_ida_two_electron.jl")
include("fullida_dense_export.jl")
include("angular_sequence_export.jl")
include("sliced_ham_export.jl")

end
