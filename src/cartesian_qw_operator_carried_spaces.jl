module CartesianQWOperatorCarriedSpaces

import ..GaussletBases:
    AbstractBondAlignedOrdinaryQWBasis3D,
    BondAlignedDiatomicQWBasis3D,
    CartesianBasisRepresentation3D,
    LegacyAtomicGaussianSupplement,
    LegacyBondAlignedDiatomicGaussianSupplement,
    LegacyBondAlignedHeteronuclearGaussianSupplement,
    MappedUniformBasis,
    OrdinaryCartesianOperators3D,
    _NestedFixedBlock3D,
    _normalized_atomic_build_context,
    _normalized_bond_aligned_build_context,
    _resolve_atomic_qw_gausslet_backend,
    _resolved_nuclear_term_storage,
    _validate_atomic_qw_nested_backend_consistency,
    _validate_nuclear_term_storage,
    _validate_operator_route_backend,
    _validate_operator_route_interaction_treatment,
    basis_representation
import ..GaussletBases.CartesianCarriedSpaces:
    CartesianCarriedSpace3D,
    carried_space_diagnostics,
    carried_space_representation,
    cartesian_carried_space

export CartesianQWOperatorCarriedSpaceSidecar,
       CartesianOperatorBuildSource3D,
       CartesianQWOperatorConstructionRecord3D,
       cartesian_qw_operator_carried_space_sidecar,
       cartesian_qw_operator_build_source,
       cartesian_qw_operator_construction_record,
       qw_operator_carried_space,
       qw_operator_basis_representation,
       qw_operator_carried_space_diagnostics,
       qw_operator_carried_space_provenance,
       operator_build_source_carried_space,
       operator_build_source_diagnostics,
       operator_build_source_provenance,
       qw_operator_construction_record_sidecar,
       qw_operator_construction_record_diagnostics,
       qw_operator_construction_record_provenance

"""
    CartesianQWOperatorCarriedSpaceSidecar

Internal sidecar that normalizes the carried Cartesian space behind an existing
`OrdinaryCartesianOperators3D` payload.

This is a bridge object only. It does not replace QW operator assembly, does
not own Hamiltonian matrices, and does not build metric or interaction packets.
"""
struct CartesianQWOperatorCarriedSpaceSidecar{C<:CartesianCarriedSpace3D,R,D,V}
    carried_space::C
    operator_representation::R
    diagnostics::D
    provenance::V
end

"""
    CartesianOperatorBuildSource3D

Internal pre-build summary for an existing Cartesian QW operator route.

This records the normalized carried space plus route/backend/storage metadata
before operator construction. It is not a builder and deliberately stores no
Hamiltonian matrices or residual packets.
"""
struct CartesianOperatorBuildSource3D{C<:CartesianCarriedSpace3D,CM,D,V}
    carried_space::C
    basis_family::Symbol
    carried_space_kind::Symbol
    nuclei::Vector{NTuple{3,Float64}}
    nuclear_charges::Vector{Float64}
    gausslet_backend::Symbol
    requested_gausslet_backend::Symbol
    interaction_treatment::Symbol
    nuclear_term_storage::Symbol
    requested_nuclear_term_storage::Symbol
    capabilities::CM
    diagnostics::D
    provenance::V
end

"""
    CartesianQWOperatorConstructionRecord3D

Internal audit record tying a pre-build `CartesianOperatorBuildSource3D` to an
already-built `OrdinaryCartesianOperators3D` payload.

The record derives the existing post-build carried-space sidecar and compares
route metadata only. It never builds operators, does not own Hamiltonian
matrices, and does not change numerical payloads.
"""
struct CartesianQWOperatorConstructionRecord3D{
    S<:CartesianOperatorBuildSource3D,
    C<:CartesianQWOperatorCarriedSpaceSidecar,
    D,
    V,
}
    build_source::S
    sidecar::C
    diagnostics::D
    provenance::V
end

qw_operator_carried_space(sidecar::CartesianQWOperatorCarriedSpaceSidecar) =
    sidecar.carried_space
qw_operator_basis_representation(sidecar::CartesianQWOperatorCarriedSpaceSidecar) =
    sidecar.operator_representation
qw_operator_carried_space_diagnostics(sidecar::CartesianQWOperatorCarriedSpaceSidecar) =
    sidecar.diagnostics
qw_operator_carried_space_provenance(sidecar::CartesianQWOperatorCarriedSpaceSidecar) =
    sidecar.provenance

operator_build_source_carried_space(source::CartesianOperatorBuildSource3D) =
    source.carried_space
operator_build_source_diagnostics(source::CartesianOperatorBuildSource3D) =
    source.diagnostics
operator_build_source_provenance(source::CartesianOperatorBuildSource3D) =
    source.provenance

qw_operator_construction_record_sidecar(record::CartesianQWOperatorConstructionRecord3D) =
    record.sidecar
qw_operator_construction_record_diagnostics(record::CartesianQWOperatorConstructionRecord3D) =
    record.diagnostics
qw_operator_construction_record_provenance(record::CartesianQWOperatorConstructionRecord3D) =
    record.provenance

function _qw_operator_carried_source(operators::OrdinaryCartesianOperators3D)
    basis = operators.basis
    basis isa _NestedFixedBlock3D && return basis, :nested_fixed_block_operator
    basis isa MappedUniformBasis && return basis, :mapped_uniform_direct_product_operator
    basis isa AbstractBondAlignedOrdinaryQWBasis3D &&
        return basis, :bond_aligned_direct_product_operator
    throw(
        ArgumentError(
            "cannot derive a Cartesian carried-space sidecar for QW operator basis $(typeof(basis))",
        ),
    )
end

function _qw_operator_sidecar_diagnostics(
    operators::OrdinaryCartesianOperators3D,
    carried_space::CartesianCarriedSpace3D,
    operator_representation::CartesianBasisRepresentation3D,
)
    carried_diagnostics = carried_space_diagnostics(carried_space)
    operator_dimension = size(operators.overlap, 1)
    return (
        operator_dimension = operator_dimension,
        operator_gausslet_count = operators.gausslet_count,
        operator_residual_count = operators.residual_count,
        raw_parent_dimension = size(operators.raw_to_final, 1),
        raw_to_final_shape = size(operators.raw_to_final),
        gausslet_backend = operators.gausslet_backend,
        interaction_treatment = operators.interaction_treatment,
        nuclear_term_storage = operators.nuclear_term_storage,
        carried_parent_dimension = carried_diagnostics.parent_dimension,
        carried_dimension = carried_diagnostics.representation_final_dimension,
        carried_has_contracted_parent = carried_diagnostics.has_contracted_parent,
        carried_has_staged_sidecar = carried_diagnostics.has_staged_sidecar,
        carried_staged_by_center_path = carried_diagnostics.staged_by_center_path,
        carried_dimension_matches_operator_gausslet_count =
            carried_diagnostics.representation_final_dimension == operators.gausslet_count,
        operator_representation_basis_kind =
            operator_representation.metadata.basis_kind,
        operator_representation_parent_kind =
            operator_representation.metadata.parent_kind,
        operator_representation_final_dimension =
            operator_representation.metadata.final_dimension,
        operator_representation_matches_operator_dimension =
            operator_representation.metadata.final_dimension == operator_dimension,
        dense_parent_matrix_used = false,
        heavy_metric_packet_built = false,
        numerical_outputs_changed = false,
    )
end

function _qw_operator_sidecar_provenance(
    operators::OrdinaryCartesianOperators3D,
    input_kind::Symbol,
)
    return (
        source = :cartesian_qw_operator_carried_space_sidecar,
        input_kind = input_kind,
        operator_type = nameof(typeof(operators)),
        basis_type = nameof(typeof(operators.basis)),
        gaussian_data_type =
            isnothing(operators.gaussian_data) ? nothing : nameof(typeof(operators.gaussian_data)),
    )
end

_maybe_property(object, name::Symbol, default = nothing) =
    hasproperty(object, name) ? getproperty(object, name) : default

function _float_nuclear_charges(charges, expected_count::Int)
    values = Float64[Float64(value) for value in charges]
    length(values) == expected_count || throw(
        ArgumentError("Cartesian QW operator build source requires one nuclear charge per nucleus"),
    )
    return values
end

function _capability_metadata(capabilities)
    return (
        route_label = _maybe_property(capabilities, :route_label),
        allowed_gausslet_backends =
            _maybe_property(capabilities, :allowed_gausslet_backends),
        backend_support_scope =
            _maybe_property(capabilities, :backend_support_scope),
        allowed_interaction_treatments =
            _maybe_property(capabilities, :allowed_interaction_treatments),
        localized_parent_kind =
            _maybe_property(capabilities, :localized_parent_kind),
        localized_parent_kind_backend =
            _maybe_property(capabilities, :localized_parent_kind_backend),
        timing_label = _maybe_property(capabilities, :timing_label),
    )
end

function _build_source_diagnostics(
    context,
    carried_space::CartesianCarriedSpace3D,
    basis_family::Symbol,
    nuclei::Vector{NTuple{3,Float64}},
    charges::Vector{Float64};
    gausslet_backend::Symbol,
    requested_gausslet_backend::Symbol,
    interaction_treatment::Symbol,
    nuclear_term_storage::Symbol,
    requested_nuclear_term_storage::Symbol,
)
    carried_diagnostics = carried_space_diagnostics(carried_space)
    capabilities = _capability_metadata(context.capabilities)
    return (
        basis_family = basis_family,
        carried_space_kind = context.carried_space_kind,
        parent_dimension = carried_diagnostics.parent_dimension,
        carried_dimension = carried_diagnostics.representation_final_dimension,
        carried_has_contracted_parent = carried_diagnostics.has_contracted_parent,
        carried_has_staged_sidecar = carried_diagnostics.has_staged_sidecar,
        carried_staged_by_center_path = carried_diagnostics.staged_by_center_path,
        nucleus_count = length(nuclei),
        nuclear_charge_count = length(charges),
        gausslet_backend = gausslet_backend,
        requested_gausslet_backend = requested_gausslet_backend,
        interaction_treatment = interaction_treatment,
        nuclear_term_storage = nuclear_term_storage,
        requested_nuclear_term_storage = requested_nuclear_term_storage,
        capability_metadata = capabilities,
        route_label = capabilities.route_label,
        dense_parent_matrix_used = false,
        heavy_metric_packet_built = false,
        operator_built = false,
    )
end

function _build_source_provenance(
    input_kind::Symbol,
    carried_source,
    gaussian_data,
    context,
)
    return (
        source = :cartesian_qw_operator_build_source,
        input_kind = input_kind,
        carried_source_type = nameof(typeof(carried_source)),
        gaussian_data_type = isnothing(gaussian_data) ? nothing : nameof(typeof(gaussian_data)),
        route_label = _maybe_property(context.capabilities, :route_label),
        carried_route_metadata = _maybe_property(context, :route_metadata),
        parent_route_metadata = _maybe_property(context, :parent_route_metadata),
    )
end

function _cartesian_operator_build_source(
    context,
    carried_source,
    gaussian_data,
    input_kind::Symbol;
    basis_family::Symbol,
    nuclei::Vector{NTuple{3,Float64}},
    nuclear_charges::Vector{Float64},
    gausslet_backend::Symbol,
    requested_gausslet_backend::Symbol,
    interaction_treatment::Symbol,
    nuclear_term_storage::Symbol,
    requested_nuclear_term_storage::Symbol,
)
    carried_space = cartesian_carried_space(carried_source)
    diagnostics = _build_source_diagnostics(
        context,
        carried_space,
        basis_family,
        nuclei,
        nuclear_charges;
        gausslet_backend,
        requested_gausslet_backend,
        interaction_treatment,
        nuclear_term_storage,
        requested_nuclear_term_storage,
    )
    provenance = _build_source_provenance(input_kind, carried_source, gaussian_data, context)
    return CartesianOperatorBuildSource3D(
        carried_space,
        basis_family,
        context.carried_space_kind,
        nuclei,
        nuclear_charges,
        gausslet_backend,
        requested_gausslet_backend,
        interaction_treatment,
        nuclear_term_storage,
        requested_nuclear_term_storage,
        _capability_metadata(context.capabilities),
        diagnostics,
        provenance,
    )
end

function _bond_aligned_build_source(
    context,
    carried_source,
    gaussian_data,
    input_kind::Symbol;
    nuclear_charges,
    nuclear_term_storage::Symbol,
    interaction_treatment::Symbol,
    gausslet_backend::Symbol,
)
    _validate_operator_route_backend(context, gausslet_backend)
    _validate_operator_route_interaction_treatment(context, interaction_treatment)
    requested_nuclear_term_storage = nuclear_term_storage
    resolved_nuclear_term_storage = _resolved_nuclear_term_storage(
        _validate_nuclear_term_storage(nuclear_term_storage),
        context.carried,
    )
    charges = _float_nuclear_charges(nuclear_charges, length(context.nuclei))
    return _cartesian_operator_build_source(
        context,
        carried_source,
        gaussian_data,
        input_kind;
        basis_family = context.basis_family,
        nuclei = NTuple{3,Float64}[context.nuclei...],
        nuclear_charges = charges,
        gausslet_backend,
        requested_gausslet_backend = gausslet_backend,
        interaction_treatment,
        nuclear_term_storage = resolved_nuclear_term_storage,
        requested_nuclear_term_storage,
    )
end

"""
    cartesian_qw_operator_build_source(...)

Return the internal pre-build summary for an existing Cartesian QW route.

The supported routes mirror existing QW entry points for direct-product and
nested fixed-block Cartesian inputs. This helper validates the same backend and
interaction capability metadata as the current builders, but it does not build
operators or numerical packets.
"""
function cartesian_qw_operator_build_source(
    basis::AbstractBondAlignedOrdinaryQWBasis3D;
    nuclear_charges::AbstractVector{<:Real} = basis.nuclear_charges,
    nuclear_term_storage::Symbol = :auto,
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
)
    context = _normalized_bond_aligned_build_context(basis)
    return _bond_aligned_build_source(
        context,
        basis,
        nothing,
        :bond_aligned_direct_product_input;
        nuclear_charges,
        nuclear_term_storage,
        interaction_treatment,
        gausslet_backend,
    )
end

function cartesian_qw_operator_build_source(
    fixed_block::_NestedFixedBlock3D{<:AbstractBondAlignedOrdinaryQWBasis3D};
    nuclear_charges::AbstractVector{<:Real} = fixed_block.parent_basis.nuclear_charges,
    nuclear_term_storage::Symbol = :auto,
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
)
    context = _normalized_bond_aligned_build_context(fixed_block)
    return _bond_aligned_build_source(
        context,
        fixed_block,
        nothing,
        :bond_aligned_nested_fixed_block_input;
        nuclear_charges,
        nuclear_term_storage,
        interaction_treatment,
        gausslet_backend,
    )
end

function cartesian_qw_operator_build_source(
    basis::BondAlignedDiatomicQWBasis3D,
    gaussian_data::Union{
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    };
    nuclear_charges::AbstractVector{<:Real} = basis.nuclear_charges,
    nuclear_term_storage::Symbol = :auto,
    interaction_treatment::Symbol = :mwg,
    gausslet_backend::Symbol = :numerical_reference,
)
    context = _normalized_bond_aligned_build_context(basis, gaussian_data)
    return _bond_aligned_build_source(
        context,
        basis,
        gaussian_data,
        :bond_aligned_molecular_direct_product_input;
        nuclear_charges,
        nuclear_term_storage,
        interaction_treatment,
        gausslet_backend,
    )
end

function cartesian_qw_operator_build_source(
    fixed_block::_NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D},
    gaussian_data::Union{
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    };
    nuclear_charges::AbstractVector{<:Real} = fixed_block.parent_basis.nuclear_charges,
    nuclear_term_storage::Symbol = :auto,
    interaction_treatment::Symbol = :mwg,
    gausslet_backend::Symbol = :numerical_reference,
)
    context = _normalized_bond_aligned_build_context(fixed_block, gaussian_data)
    return _bond_aligned_build_source(
        context,
        fixed_block,
        gaussian_data,
        :bond_aligned_molecular_nested_fixed_block_input;
        nuclear_charges,
        nuclear_term_storage,
        interaction_treatment,
        gausslet_backend,
    )
end

function _atomic_build_source(
    context,
    carried_source,
    gaussian_data::LegacyAtomicGaussianSupplement,
    input_kind::Symbol;
    Z::Real,
    interaction_treatment::Symbol,
    gausslet_backend::Symbol,
    nuclear_term_storage::Symbol,
)
    nuclear_term_storage == :total_only || throw(
        ArgumentError("atomic QW operator build source currently supports nuclear_term_storage = :total_only"),
    )
    resolved_gausslet_backend = _resolve_atomic_qw_gausslet_backend(context, gausslet_backend)
    _validate_operator_route_backend(context, resolved_gausslet_backend)
    _validate_atomic_qw_nested_backend_consistency(context, resolved_gausslet_backend)
    _validate_operator_route_interaction_treatment(context, interaction_treatment)
    return _cartesian_operator_build_source(
        context,
        carried_source,
        gaussian_data,
        input_kind;
        basis_family = :one_center_atomic,
        nuclei = NTuple{3,Float64}[(0.0, 0.0, 0.0)],
        nuclear_charges = Float64[Float64(Z)],
        gausslet_backend = resolved_gausslet_backend,
        requested_gausslet_backend = gausslet_backend,
        interaction_treatment,
        nuclear_term_storage = :total_only,
        requested_nuclear_term_storage = nuclear_term_storage,
    )
end

function cartesian_qw_operator_build_source(
    basis::MappedUniformBasis,
    gaussian_data::LegacyAtomicGaussianSupplement;
    Z::Real = 2.0,
    interaction_treatment::Symbol = :mwg,
    gausslet_backend::Symbol = :auto,
    nuclear_term_storage::Symbol = :total_only,
)
    context = _normalized_atomic_build_context(basis, gaussian_data)
    return _atomic_build_source(
        context,
        basis,
        gaussian_data,
        :atomic_direct_product_input;
        Z,
        interaction_treatment,
        gausslet_backend,
        nuclear_term_storage,
    )
end

function cartesian_qw_operator_build_source(
    fixed_block::_NestedFixedBlock3D,
    gaussian_data::LegacyAtomicGaussianSupplement;
    Z::Real = 2.0,
    interaction_treatment::Symbol = :mwg,
    gausslet_backend::Symbol = :auto,
    nuclear_term_storage::Symbol = :total_only,
)
    context = _normalized_atomic_build_context(fixed_block, gaussian_data)
    return _atomic_build_source(
        context,
        fixed_block,
        gaussian_data,
        :atomic_nested_fixed_block_input;
        Z,
        interaction_treatment,
        gausslet_backend,
        nuclear_term_storage,
    )
end

function _qw_operator_basis_representation(
    operators::OrdinaryCartesianOperators3D,
    carried_space::CartesianCarriedSpace3D,
)
    if operators.residual_count == 0 &&
       size(operators.raw_to_final, 1) == operators.gausslet_count &&
       size(operators.raw_to_final, 2) == operators.gausslet_count
        return carried_space_representation(carried_space)
    end
    return basis_representation(operators)
end

"""
    cartesian_qw_operator_carried_space_sidecar(operators)

Derive the normalized carried-space sidecar for an existing Cartesian QW
operator payload.

Supported carried sources are direct-product mapped/bond-aligned bases and
nested fixed blocks. Unsupported routes throw instead of guessing a parent or
falling back to numerical quadrature.

Pure carried-only operators reuse the carried-space representation. Hybrid or
residual operators derive the existing full operator representation.
"""
function cartesian_qw_operator_carried_space_sidecar(
    operators::OrdinaryCartesianOperators3D,
)
    carried_source, input_kind = _qw_operator_carried_source(operators)
    carried_space = cartesian_carried_space(carried_source)
    operator_representation = _qw_operator_basis_representation(operators, carried_space)
    diagnostics = _qw_operator_sidecar_diagnostics(
        operators,
        carried_space,
        operator_representation,
    )
    provenance = _qw_operator_sidecar_provenance(operators, input_kind)
    return CartesianQWOperatorCarriedSpaceSidecar(
        carried_space,
        operator_representation,
        diagnostics,
        provenance,
    )
end

_record_field_comparison(field::Symbol, source_value, operator_value) = (
    field = field,
    source_value = source_value,
    operator_value = operator_value,
    matches = source_value == operator_value,
    ambiguous = false,
)

_record_nuclear_charge_comparison(source::CartesianOperatorBuildSource3D, operators) = (
    field = :nuclear_charges,
    source_value = source.nuclear_charges,
    operator_value = operators.nuclear_charges,
    matches = operators.nuclear_charges !== nothing &&
              source.nuclear_charges == operators.nuclear_charges,
    ambiguous = operators.nuclear_charges === nothing,
)

function _expected_operator_input_kind(source_input_kind::Symbol)
    source_input_kind == :atomic_direct_product_input &&
        return :mapped_uniform_direct_product_operator
    source_input_kind == :bond_aligned_direct_product_input &&
        return :bond_aligned_direct_product_operator
    source_input_kind == :bond_aligned_molecular_direct_product_input &&
        return :bond_aligned_direct_product_operator
    source_input_kind == :atomic_nested_fixed_block_input &&
        return :nested_fixed_block_operator
    source_input_kind == :bond_aligned_nested_fixed_block_input &&
        return :nested_fixed_block_operator
    source_input_kind == :bond_aligned_molecular_nested_fixed_block_input &&
        return :nested_fixed_block_operator
    return :unknown
end

function _construction_record_comparisons(
    source::CartesianOperatorBuildSource3D,
    operators::OrdinaryCartesianOperators3D,
    sidecar::CartesianQWOperatorCarriedSpaceSidecar,
)
    source_diagnostics = operator_build_source_diagnostics(source)
    sidecar_diagnostics = qw_operator_carried_space_diagnostics(sidecar)
    source_provenance = operator_build_source_provenance(source)
    sidecar_provenance = qw_operator_carried_space_provenance(sidecar)
    return (
        operator_input_kind = _record_field_comparison(
            :operator_input_kind,
            _expected_operator_input_kind(source_provenance.input_kind),
            sidecar_provenance.input_kind,
        ),
        gausslet_backend = _record_field_comparison(
            :gausslet_backend,
            source.gausslet_backend,
            operators.gausslet_backend,
        ),
        interaction_treatment = _record_field_comparison(
            :interaction_treatment,
            source.interaction_treatment,
            operators.interaction_treatment,
        ),
        nuclear_term_storage = _record_field_comparison(
            :nuclear_term_storage,
            source.nuclear_term_storage,
            operators.nuclear_term_storage,
        ),
        nuclear_charges = _record_nuclear_charge_comparison(source, operators),
        carried_dimension = _record_field_comparison(
            :carried_dimension,
            source_diagnostics.carried_dimension,
            sidecar_diagnostics.carried_dimension,
        ),
        carried_has_contracted_parent = _record_field_comparison(
            :carried_has_contracted_parent,
            source_diagnostics.carried_has_contracted_parent,
            sidecar_diagnostics.carried_has_contracted_parent,
        ),
        carried_has_staged_sidecar = _record_field_comparison(
            :carried_has_staged_sidecar,
            source_diagnostics.carried_has_staged_sidecar,
            sidecar_diagnostics.carried_has_staged_sidecar,
        ),
        carried_staged_by_center_path = _record_field_comparison(
            :carried_staged_by_center_path,
            source_diagnostics.carried_staged_by_center_path,
            sidecar_diagnostics.carried_staged_by_center_path,
        ),
    )
end

function _construction_record_mismatch_fields(comparisons)
    return Symbol[
        name for name in keys(comparisons) if !getfield(comparisons, name).matches
    ]
end

function _construction_record_ambiguous_fields(comparisons)
    return Symbol[
        name for name in keys(comparisons) if getfield(comparisons, name).ambiguous
    ]
end

function _construction_record_diagnostics(
    source::CartesianOperatorBuildSource3D,
    operators::OrdinaryCartesianOperators3D,
    sidecar::CartesianQWOperatorCarriedSpaceSidecar,
)
    source_diagnostics = operator_build_source_diagnostics(source)
    sidecar_diagnostics = qw_operator_carried_space_diagnostics(sidecar)
    comparisons = _construction_record_comparisons(source, operators, sidecar)
    mismatch_fields = _construction_record_mismatch_fields(comparisons)
    ambiguous_fields = _construction_record_ambiguous_fields(comparisons)
    return (
        source_basis_family = source.basis_family,
        source_carried_space_kind = source.carried_space_kind,
        source_input_kind = operator_build_source_provenance(source).input_kind,
        sidecar_input_kind = qw_operator_carried_space_provenance(sidecar).input_kind,
        operator_dimension = sidecar_diagnostics.operator_dimension,
        operator_gausslet_count = sidecar_diagnostics.operator_gausslet_count,
        operator_residual_count = sidecar_diagnostics.operator_residual_count,
        compared_fields = keys(comparisons),
        comparisons = comparisons,
        mismatch_fields = mismatch_fields,
        ambiguous_mismatch_fields = ambiguous_fields,
        source_sidecar_agree = isempty(mismatch_fields),
        source_carried_dimension = source_diagnostics.carried_dimension,
        sidecar_carried_dimension = sidecar_diagnostics.carried_dimension,
        dense_parent_matrix_used = false,
        heavy_metric_packet_built = false,
        operator_built = false,
        numerical_outputs_changed = false,
    )
end

function _construction_record_provenance(
    source::CartesianOperatorBuildSource3D,
    operators::OrdinaryCartesianOperators3D,
    sidecar::CartesianQWOperatorCarriedSpaceSidecar,
)
    return (
        source = :cartesian_qw_operator_construction_record,
        build_source_provenance = operator_build_source_provenance(source),
        sidecar_provenance = qw_operator_carried_space_provenance(sidecar),
        operator_type = nameof(typeof(operators)),
        basis_type = nameof(typeof(operators.basis)),
        gaussian_data_type =
            isnothing(operators.gaussian_data) ? nothing : nameof(typeof(operators.gaussian_data)),
    )
end

"""
    cartesian_qw_operator_construction_record(source, operators; throw_on_mismatch=true)

Return an internal audit record tying an existing
`CartesianOperatorBuildSource3D` to an already-built
`OrdinaryCartesianOperators3D` payload.

This helper derives the post-build carried-space sidecar and compares route
metadata. It does not build operators, dispatch kernels, construct metric
packets, or mutate numerical matrices.
"""
function cartesian_qw_operator_construction_record(
    source::CartesianOperatorBuildSource3D,
    operators::OrdinaryCartesianOperators3D;
    throw_on_mismatch::Bool = true,
)
    sidecar = cartesian_qw_operator_carried_space_sidecar(operators)
    diagnostics = _construction_record_diagnostics(source, operators, sidecar)
    if throw_on_mismatch && !diagnostics.source_sidecar_agree
        throw(
            ArgumentError(
                "Cartesian QW operator construction record mismatch for fields " *
                join(string.(diagnostics.mismatch_fields), ", "),
            ),
        )
    end
    provenance = _construction_record_provenance(source, operators, sidecar)
    return CartesianQWOperatorConstructionRecord3D(
        source,
        sidecar,
        diagnostics,
        provenance,
    )
end

end
