using GaussletBases, Test
const _DOCUMENTED_PUBLIC_SURFACE = (
    v0_2 = (
        :PQSH2PlusRow, :PQSH2PlusComparison, :pqs_h2plus_comparison,
        :ExactRepresentedHartreeField, :FittedReferenceHartreeField,
        :ScreenedHartreeCorrection, :screened_hartree_correction,
        :screened_hartree_delta_one_body, :screened_hartree_energy_constant,
        :screened_hartree_consistency_error, :screened_hartree_field_kind,
        :cartesian_residual_gto_mwg_system),
    external_gto = (
        :ExternalGTOOrbitalSpinBlock, :ExternalGTOOrbitalPacket,
        :ExternalGTOOrbitalSpinImport, :ExternalGTOOrbitalImportResult,
        :import_external_gto_orbitals, :read_external_cartesian_gto_packet,
        :closest_external_gto_determinant),
    foundation = (
        :uofx, :xofu, :dudx, :du2dx2, :basis_spec,
        :family, :mapping, :centers, :reference_centers, :integral_weights),
    function_stencil = (
        :value, :direct_value, :derivative, :center, :reference_center,
        :integral_weight, :stencil, :coefficients, :terms),
    partition_leaf = (
        :boxes, :leaf_boxes, :box_indices, :box_level, :box_parent, :box_children,
        :box_block, :box_coupling, :leaf_primitive_indices, :primitive_origins,
        :primitive_leaf_boxes, :leaf_contractions),
    atomic_ida = (
        :orbitals, :two_electron_states, :radial_multipole, :gaunt_coefficient,
        :gaunt_tensor, :angular_kernel, :apply_overlap, :apply_hamiltonian,
        :ground_state_energy, :lanczos_ground_state),
    supported = (
        :BondAlignedDiatomicQWBasis3D, :CoulombGaussianExpansion, :basis_metadata,
        :cartesian_base_hamiltonian, :external_gto_ordering_fingerprint,
        :external_gto_overlap_fingerprint),
    geometry = (
        :bond_aligned_diatomic_geometry_payload, :BondAlignedDiatomicGeometryPoint3D,
        :BondAlignedDiatomicGeometryNucleus3D, :BondAlignedDiatomicGeometryBox3D,
        :BondAlignedDiatomicGeometryPayload3D, :BondAlignedDiatomicGeometryPlaneSlice3D,
        :bond_aligned_diatomic_source_geometry_payload, :bond_aligned_diatomic_plane_slice),
    angular = (
        :ShellLocalAngularProfile, :ShellLocalAngularProfileOverlap,
        :shell_local_angular_profile, :adjacent_shell_local_angular_profile_overlap,
        :AtomicFixedRadialAngularSequenceLevel,
        :AtomicFixedRadialAngularSequenceOverlapSidecar, :AtomicFixedRadialAngularSequence),
    radial_parity = (
        :RadialBoundaryPrototype, :radial_boundary_prototype,
        :radial_boundary_prototype_names, :build_paper_parity_radial_basis),
    qw_geometry = (
        :bond_aligned_diatomic_nested_geometry_diagnostics, :bond_aligned_homonuclear_chain_geometry_diagnostics,
        :axis_aligned_homonuclear_square_lattice_geometry_diagnostics),
    sliced_chain = (
        :SlicedHydrogenChain, :sliced_hydrogen_chain, :sliced_h1,
        :sliced_h1_bandwidth, :sliced_vee, Symbol("sliced_row!")),
    working_basis = (:cartesian_base_working_basis,),
)

const _RESERVED_UNDOCUMENTED_PUBLIC_SURFACE = Set((
    :OneCenterAtomicNestedStructureDiagnostics, :one_center_atomic_nested_structure_diagnostics,
    :one_center_atomic_nested_structure_report, :diagnose_qwrg_residual_space, :ShellLocalAngularProfileKey,
))
# Mirror Julia 1.12's documentation lookup for the supported Julia 1.10 floor.
function _public_surface_hasdoc(binding::Base.Docs.Binding)
    docs = Base.Docs
    if docs.defined(binding) && !isnothing(docs.getdoc(docs.resolve(binding), Union{}))
        return true
    end
    for owner in docs.modules
        metadata = docs.meta(owner; autoinit = false)
        metadata !== nothing && haskey(metadata, binding) && return true
    end
    alias = docs.aliasof(binding)
    return alias == binding ? false : _public_surface_hasdoc(alias)
end

_public_surface_hasdoc(owner::Module, name::Symbol) = _public_surface_hasdoc(Base.Docs.Binding(owner, name))

function _public_surface_undocumented_names(owner::Module)
    return Set(name for name in names(owner) if name != nameof(owner) && !_public_surface_hasdoc(owner, name))
end

@testset "Public surface documentation" begin
    public_names = names(GaussletBases)
    documented_names = collect(Iterators.flatten(values(_DOCUMENTED_PUBLIC_SURFACE)))
    @test length(documented_names) == 95
    @test allunique(documented_names)
    @test all(name -> name in public_names, documented_names)
    @test all(name -> _public_surface_hasdoc(GaussletBases, name), documented_names)
    alias_binding = Base.Docs.Binding(GaussletBases, :CuratedSpherePointSet)
    @test Base.Docs.aliasof(alias_binding) != alias_binding && _public_surface_hasdoc(alias_binding)
    undocumented = _public_surface_undocumented_names(GaussletBases)
    @test undocumented == _RESERVED_UNDOCUMENTED_PUBLIC_SURFACE
    if isdefined(Base.Docs, :hasdoc)
        native_hasdoc, native_undocumented_names =
            getfield(Base.Docs, :hasdoc), getfield(Base.Docs, :undocumented_names)
        @test all(name -> _public_surface_hasdoc(GaussletBases, name) ==
                          native_hasdoc(GaussletBases, name), public_names)
        @test undocumented == Set(filter(!=(nameof(GaussletBases)),
            native_undocumented_names(GaussletBases)))
    end
end
