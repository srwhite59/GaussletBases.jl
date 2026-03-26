module GaussletBases

using JLD2
using LinearAlgebra
using SpecialFunctions
using TOML

export AbstractFunction1D,
       AbstractPrimitiveFunction1D,
       AbstractBasisFunction1D,
       AbstractCoordinateMapping,
       AbstractBasisSpec,
       AbstractDiagonalApproximation,
       PrimitiveSet1D,
       BasisMetadata1D,
       BasisRepresentation1D,
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
       HybridMappedOrdinaryBasis1D,
       LegacyAtomicGaussianShell,
       LegacyAtomicGaussianSupplement,
       LegacyBondAlignedDiatomicGaussianSupplement,
       LegacyBondAlignedHeteronuclearGaussianSupplement,
       LegacySGaussianData,
       BondAlignedDiatomicQWBasis3D,
       BondAlignedDiatomicGeometryPoint3D,
       BondAlignedDiatomicGeometryNucleus3D,
       BondAlignedDiatomicGeometryBox3D,
       BondAlignedDiatomicGeometryPayload3D,
       BondAlignedDiatomicGeometryPlaneSlice3D,
       MappedOrdinaryOneBody1D,
       CartesianProductOrbital3D,
       QiuWhiteHybridOrbital3D,
       OrdinaryCartesianIDAOperators,
       QiuWhiteResidualGaussianOperators,
       HalfLineBasisSpec,
       RadialBasisSpec,
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
       mapped_pgdg_prototype,
       mapped_pgdg_localized,
       hybrid_mapped_ordinary_basis,
       legacy_atomic_gaussian_supplement,
       legacy_bond_aligned_diatomic_gaussian_supplement,
       legacy_bond_aligned_heteronuclear_gaussian_supplement,
       legacy_s_gaussian_data,
       bond_aligned_homonuclear_qw_basis,
       bond_aligned_heteronuclear_qw_basis,
       bond_aligned_diatomic_geometry_payload,
       bond_aligned_diatomic_source_geometry_payload,
       bond_aligned_diatomic_plane_slice,
       write_bond_aligned_diatomic_points3d,
       write_bond_aligned_diatomic_plane_projection,
       SpherePointSetProvenance,
       CuratedSpherePointSet,
       ShellLocalInjectedAngularBasis,
       AtomicShellLocalInjectedAngularAssembly,
       AtomicInjectedAngularOneBodyBenchmark,
       AtomicInjectedAngularHFStyleBenchmark,
       AtomicInjectedAngularSmallEDBenchmark,
       curated_sphere_point_set,
       curated_sphere_point_set_orders,
       build_shell_local_injected_angular_basis,
       assign_atomic_angular_shell_orders,
       build_atomic_shell_local_angular_assembly,
       build_atomic_injected_angular_one_body_benchmark,
       build_atomic_injected_angular_hf_style_benchmark,
       build_atomic_injected_angular_small_ed_benchmark,
       shell_local_injected_angular_diagnostics,
       atomic_shell_local_angular_diagnostics,
       atomic_injected_angular_one_body_diagnostics,
       atomic_injected_angular_hf_style_diagnostics,
       atomic_injected_angular_small_ed_diagnostics,
       mapped_ordinary_one_body_operators,
       mapped_cartesian_hydrogen_energy,
       ordinary_sho_hamiltonian,
       ordinary_sho_spectrum,
       ordinary_cartesian_ida_operators,
       ordinary_cartesian_qiu_white_operators,
       ordinary_cartesian_vee_expectation,
       angular_benchmark_exact_hamv6_payload,
       atomic_hamv6_payload,
       atomic_ida_density_interaction_matrix,
       fullida_dense_payload,
       sliced_ham_payload,
       write_angular_benchmark_exact_hamv6_jld2,
       write_atomic_hamv6_jld2,
       write_fullida_dense_jld2,
       write_sliced_ham_jld2,
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

function value end
function direct_value end
function derivative end
function center end
function reference_center end
function moment_center end
function integral_weight end
function stencil end
function stencil_matrix end
function build_basis end
function basis_metadata end
function basis_representation end
function basis_partition end
function hierarchical_partition end
function build_leaf_pgdg end
function augment_leaf_pgdg end
function build_global_mapped_primitive_layer end
function contract_leaf_boxes end
function refine_partition end
function primitive_set end
function boxes end
function leaf_boxes end
function leaf_contractions end
function leaf_primitive_indices end
function primitive_origins end
function primitive_leaf_boxes end
function box_indices end
function box_level end
function box_parent end
function box_children end
function box_block end
function box_coupling end

function basis_spec end
function family end
function mapping end
function centers end
function reference_centers end
function integral_weights end
function contract_primitive_vector end
function contract_primitive_diagonal end
function contract_primitive_matrix end
function coulomb_gaussian_expansion end
function mapped_pgdg_prototype end
function mapped_pgdg_localized end
function hybrid_mapped_ordinary_basis end
function bond_aligned_diatomic_geometry_payload end
function bond_aligned_diatomic_source_geometry_payload end
function bond_aligned_diatomic_plane_slice end
function write_bond_aligned_diatomic_points3d end
function write_bond_aligned_diatomic_plane_projection end
function mapped_ordinary_one_body_operators end
function mapped_cartesian_hydrogen_energy end
function ordinary_sho_hamiltonian end
function ordinary_sho_spectrum end
function ordinary_cartesian_ida_operators end
function ordinary_cartesian_qiu_white_operators end
function ordinary_cartesian_1s2_check end
function ordinary_cartesian_vee_expectation end
function angular_benchmark_exact_hamv6_payload end
function fullida_dense_payload end
function sliced_ham_payload end
function write_angular_benchmark_exact_hamv6_jld2 end
function write_fullida_dense_jld2 end
function write_sliced_ham_jld2 end
function gaussian_factor_matrix end
function gaussian_factor_matrices end
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
function orbitals end
function two_electron_states end
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
function gaunt_tensor end
function gaunt_coefficient end
function angular_kernel end
function apply_overlap end
function apply_hamiltonian end
function ground_state_energy end
function lanczos_ground_state end
function atomic_operators end
function centrifugal end
function multipole end

function coefficients end
function primitives end
function terms end

function uofx end
function xofu end
function dudx end
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
include("ordinary_qiu_white_rg.jl")
include("bond_aligned_diatomic_geometry.jl")
include("bond_aligned_diatomic_geometry_export.jl")
include("partitions.jl")
include("hierarchical_partitions.jl")
include("leaf_pgdg.jl")
include("global_leaf_contraction.jl")
include("diagnostics.jl")
include("operators.jl")
include("atomic_ylm.jl")
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
include("sliced_ham_export.jl")

end
