module CartesianQWOperatorCarriedSpaces

using LinearAlgebra: I, Symmetric, diag, norm

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
    _nested_bond_aligned_diatomic_high_order_recipe_policy,
    _nested_bond_aligned_diatomic_high_order_recipe_policy_diagnostics,
    _nested_bond_aligned_diatomic_high_order_recipe_realization_audit,
    _nested_bond_aligned_diatomic_high_order_recipe_realization_diagnostics,
    _nested_bond_aligned_diatomic_high_order_recipe_source_construction,
    _nested_bond_aligned_diatomic_high_order_recipe_source_construction_diagnostics,
    _nested_bond_aligned_diatomic_high_order_recipe_source_fixed_block,
    _nested_bond_aligned_diatomic_high_order_recipe_source_readiness,
    _nested_bond_aligned_diatomic_high_order_recipe_source_readiness_diagnostics,
    _qwrg_bond_aligned_axis_bundles,
    _normalized_atomic_build_context,
    _normalized_bond_aligned_build_context,
    _resolve_atomic_qw_gausslet_backend,
    _resolved_nuclear_term_storage,
    _validate_atomic_qw_nested_backend_consistency,
    _validate_nuclear_term_storage,
    _validate_operator_route_backend,
    _validate_operator_route_interaction_treatment,
    basis_representation,
    bond_aligned_homonuclear_qw_basis,
    coulomb_gaussian_expansion,
    gto_overlap_matrix,
    legacy_bond_aligned_diatomic_gaussian_supplement,
    ordinary_cartesian_qiu_white_operators
import ..GaussletBases.CartesianCarriedSpaces:
    CartesianCarriedSpace3D,
    carried_space_diagnostics,
    carried_space_parent,
    carried_space_provenance,
    carried_space_representation,
    cartesian_carried_space

export CartesianQWOperatorCarriedSpaceSidecar,
       CartesianOperatorBuildSource3D,
       CartesianQWOperatorConstructionRecord3D,
       CartesianQWOperatorConstructionReceipt3D,
       cartesian_qw_operator_carried_space_sidecar,
       cartesian_qw_operator_build_source,
       cartesian_qw_operator_construction_record,
       cartesian_qw_operator_construction_receipt,
       cartesian_qw_operator_receipt_coverage,
       qw_operator_carried_space,
       qw_operator_basis_representation,
       qw_operator_carried_space_diagnostics,
       qw_operator_carried_space_provenance,
       operator_build_source_carried_space,
       operator_build_source_diagnostics,
       operator_build_source_provenance,
       qw_operator_construction_record_sidecar,
       qw_operator_construction_record_diagnostics,
       qw_operator_construction_record_provenance,
       qw_operator_construction_receipt_source,
       qw_operator_construction_receipt_operators,
       qw_operator_construction_receipt_record,
       qw_operator_construction_receipt_diagnostics,
       qw_operator_construction_receipt_provenance

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
route metadata plus lightweight carried-space identity metadata. It never
builds operators, does not own Hamiltonian matrices, and does not change
numerical payloads.
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

"""
    CartesianQWOperatorConstructionReceipt3D

Internal audited construction receipt for existing Cartesian QW builders.

The receipt builds a pre-build `CartesianOperatorBuildSource3D`, delegates all
numerical work to the matching `ordinary_cartesian_qiu_white_operators` method,
and then records the existing source/operator audit. It is deliberately a thin
wrapper: it does not implement Hamiltonian kernels, construct metric packets,
or change backend/route semantics.
"""
struct CartesianQWOperatorConstructionReceipt3D{
    S<:CartesianOperatorBuildSource3D,
    O<:OrdinaryCartesianOperators3D,
    R<:CartesianQWOperatorConstructionRecord3D,
    D,
    V,
}
    build_source::S
    operators::O
    construction_record::R
    diagnostics::D
    provenance::V
end

"""
    _BondAlignedDiatomicHighOrderQRowRouteReceipt3D

Internal q-row route receipt for the experimental high-order diatomic
atom-growth/endcap-panel path. It records existing policy, source, fixed-block,
and QW receipt objects for diagnostics only; all numerical work is delegated to
the existing nested source and QW builders.
"""
struct _BondAlignedDiatomicHighOrderQRowRouteReceipt3D{
    B,
    E,
    A,
    P,
    PD,
    RA,
    RD,
    C,
    CD,
    R,
    RDD,
    F,
    QR,
    QD,
    QCD,
    D,
    PR,
}
    basis::B
    expansion::E
    axis_bundles::A
    policy::P
    policy_diagnostics::PD
    realization_audit::RA
    realization_diagnostics::RD
    source_construction::C
    source_construction_diagnostics::CD
    readiness::R
    readiness_diagnostics::RDD
    fixed_block::F
    qw_receipt::QR
    qw_receipt_diagnostics::QD
    qw_record_diagnostics::QCD
    diagnostics::D
    provenance::PR
end

"""
    _BondAlignedHomonuclearHighOrderQRowFixtureReceipt3D

Internal homonuclear fixture receipt for the experimental q-row route. It
constructs the existing `BondAlignedDiatomicQWBasis3D` from explicit geometry,
delegates numerical work to the basis-object q-row route receipt, and records
fixture provenance only.
"""
struct _BondAlignedHomonuclearHighOrderQRowFixtureReceipt3D{B,QR,D,V}
    basis::B
    q_row_route_receipt::QR
    diagnostics::D
    provenance::V
end

"""
    _BondAlignedHomonuclearHighOrderQRowFixtureSupplementReceipt3D

Internal supplement-aware homonuclear q-row fixture receipt. It composes the
private homonuclear q-row fixture wrapper with the existing nested molecular
supplement QW receipt; it does not implement Hamiltonian kernels.
"""
struct _BondAlignedHomonuclearHighOrderQRowFixtureSupplementReceipt3D{FR,S,QR,D,V}
    fixture_receipt::FR
    supplement::S
    supplement_qw_receipt::QR
    diagnostics::D
    provenance::V
end

"""
    _BondAlignedHomonuclearHighOrderQRowFixtureSupplementCaptureH1Diagnostic3D

Private parent-grid target capture/H1 diagnostic for the homonuclear q-row
fixture supplement route. It reuses the existing q-row/supplement receipt and
QW operators; it does not implement Hamiltonian kernels.
"""
struct _BondAlignedHomonuclearHighOrderQRowFixtureSupplementCaptureH1Diagnostic3D{
    R,
    P,
    T,
    F,
    H,
    ROWS,
    S,
    D,
    V,
}
    route_receipt::R
    parent_qw_receipt::P
    target_coefficients::T
    fixed_projected_coefficients::F
    final_projected_coefficients::H
    target_rows::ROWS
    summary::S
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

qw_operator_construction_receipt_source(
    receipt::CartesianQWOperatorConstructionReceipt3D,
) = receipt.build_source
qw_operator_construction_receipt_operators(
    receipt::CartesianQWOperatorConstructionReceipt3D,
) = receipt.operators
qw_operator_construction_receipt_record(
    receipt::CartesianQWOperatorConstructionReceipt3D,
) = receipt.construction_record
qw_operator_construction_receipt_diagnostics(
    receipt::CartesianQWOperatorConstructionReceipt3D,
) = receipt.diagnostics
qw_operator_construction_receipt_provenance(
    receipt::CartesianQWOperatorConstructionReceipt3D,
) = receipt.provenance

_nested_bond_aligned_diatomic_high_order_q_row_route_diagnostics(
    receipt::_BondAlignedDiatomicHighOrderQRowRouteReceipt3D,
) = receipt.diagnostics

_nested_bond_aligned_diatomic_high_order_q_row_route_provenance(
    receipt::_BondAlignedDiatomicHighOrderQRowRouteReceipt3D,
) = receipt.provenance

_nested_bond_aligned_homonuclear_high_order_q_row_fixture_diagnostics(
    receipt::_BondAlignedHomonuclearHighOrderQRowFixtureReceipt3D,
) = receipt.diagnostics

_nested_bond_aligned_homonuclear_high_order_q_row_fixture_provenance(
    receipt::_BondAlignedHomonuclearHighOrderQRowFixtureReceipt3D,
) = receipt.provenance

_nested_bond_aligned_homonuclear_high_order_q_row_fixture_supplement_diagnostics(
    receipt::_BondAlignedHomonuclearHighOrderQRowFixtureSupplementReceipt3D,
) = receipt.diagnostics

_nested_bond_aligned_homonuclear_high_order_q_row_fixture_supplement_provenance(
    receipt::_BondAlignedHomonuclearHighOrderQRowFixtureSupplementReceipt3D,
) = receipt.provenance

_nested_bond_aligned_homonuclear_high_order_q_row_fixture_supplement_capture_h1_diagnostics(
    diagnostic::_BondAlignedHomonuclearHighOrderQRowFixtureSupplementCaptureH1Diagnostic3D,
) = diagnostic.diagnostics

_nested_bond_aligned_homonuclear_high_order_q_row_fixture_supplement_capture_h1_rows(
    diagnostic::_BondAlignedHomonuclearHighOrderQRowFixtureSupplementCaptureH1Diagnostic3D,
) = diagnostic.target_rows

_nested_bond_aligned_homonuclear_high_order_q_row_fixture_supplement_capture_h1_summary(
    diagnostic::_BondAlignedHomonuclearHighOrderQRowFixtureSupplementCaptureH1Diagnostic3D,
) = diagnostic.summary

_nested_bond_aligned_homonuclear_high_order_q_row_fixture_supplement_capture_h1_provenance(
    diagnostic::_BondAlignedHomonuclearHighOrderQRowFixtureSupplementCaptureH1Diagnostic3D,
) = diagnostic.provenance

const _CARTESIAN_QW_OPERATOR_RECEIPT_COVERAGE = (
    contract = :delegate_plus_audit_only,
    covered_route_families = (
        (
            route_family = :bond_aligned_direct_product,
            input_kind = :bond_aligned_direct_product_input,
            inputs = "AbstractBondAlignedOrdinaryQWBasis3D",
            supplement = nothing,
            builder = :ordinary_cartesian_qiu_white_operators,
            status = :covered,
        ),
        (
            route_family = :bond_aligned_nested_fixed_block,
            input_kind = :bond_aligned_nested_fixed_block_input,
            inputs = "_NestedFixedBlock3D{<:AbstractBondAlignedOrdinaryQWBasis3D}",
            supplement = nothing,
            builder = :ordinary_cartesian_qiu_white_operators,
            status = :covered,
        ),
        (
            route_family = :atomic_direct_product_supplement,
            input_kind = :atomic_direct_product_input,
            inputs = "MappedUniformBasis + LegacyAtomicGaussianSupplement",
            supplement = :legacy_atomic_gaussian,
            builder = :ordinary_cartesian_qiu_white_operators,
            status = :covered,
        ),
        (
            route_family = :atomic_nested_fixed_block_supplement,
            input_kind = :atomic_nested_fixed_block_input,
            inputs = "_NestedFixedBlock3D + LegacyAtomicGaussianSupplement",
            supplement = :legacy_atomic_gaussian,
            builder = :ordinary_cartesian_qiu_white_operators,
            status = :covered,
        ),
        (
            route_family = :bond_aligned_molecular_direct_product_supplement,
            input_kind = :bond_aligned_molecular_direct_product_input,
            inputs = "BondAlignedDiatomicQWBasis3D + molecular LegacyGaussianSupplement",
            supplement = :legacy_bond_aligned_molecular_gaussian,
            builder = :ordinary_cartesian_qiu_white_operators,
            status = :covered,
        ),
        (
            route_family = :bond_aligned_molecular_nested_fixed_block_supplement,
            input_kind = :bond_aligned_molecular_nested_fixed_block_input,
            inputs = "_NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D} + molecular LegacyGaussianSupplement",
            supplement = :legacy_bond_aligned_molecular_gaussian,
            builder = :ordinary_cartesian_qiu_white_operators,
            status = :covered,
        ),
    ),
    intentionally_uncovered_route_families = (
        (
            route_family = :alias_frontends,
            inputs = "ordinary_cartesian_product_operators / nested_cartesian_operators",
            reason = "aliases remain direct-builder-only; call the receipt helper with the underlying canonical basis or fixed-block inputs",
        ),
        (
            route_family = :already_built_operators,
            inputs = "OrdinaryCartesianOperators3D",
            reason = "receipts are pre-build wrappers; use cartesian_qw_operator_construction_record or cartesian_qw_operator_carried_space_sidecar for already-built operators",
        ),
        (
            route_family = :non_diatomic_molecular_supplement,
            inputs = "chain/square or arbitrary fixed-block inputs with molecular LegacyGaussianSupplement",
            reason = "there is no existing direct molecular supplement builder route to delegate to",
        ),
        (
            route_family = :atomic_without_supplement,
            inputs = "MappedUniformBasis without LegacyAtomicGaussianSupplement",
            reason = "the receipt layer currently covers QW supplement or bond-aligned product routes, not mapped one-body-only helpers",
        ),
    ),
)

cartesian_qw_operator_receipt_coverage() = _CARTESIAN_QW_OPERATOR_RECEIPT_COVERAGE

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
    matches = isequal(source_value, operator_value),
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

function _carried_space_identity_metadata(space::CartesianCarriedSpace3D)
    diagnostics = carried_space_diagnostics(space)
    representation = carried_space_representation(space)
    provenance = carried_space_provenance(space)
    parent = carried_space_parent(space)
    return (
        parent_axis_counts = diagnostics.parent_axis_counts,
        parent_dimension = diagnostics.parent_dimension,
        representation_basis_kind = diagnostics.representation_basis_kind,
        representation_parent_kind = representation.metadata.parent_kind,
        representation_final_dimension = diagnostics.representation_final_dimension,
        axis_sharing = parent.axis_sharing,
        has_contracted_parent = diagnostics.has_contracted_parent,
        has_staged_sidecar = diagnostics.has_staged_sidecar,
        staged_by_center_path = diagnostics.staged_by_center_path,
        provenance_input_kind = provenance.input_kind,
        provenance_source_type = provenance.source_type,
        provenance_representation_basis_kind = provenance.representation_basis_kind,
        provenance_route_metadata = provenance.route_metadata,
    )
end

function _construction_record_comparisons(
    source::CartesianOperatorBuildSource3D,
    operators::OrdinaryCartesianOperators3D,
    sidecar::CartesianQWOperatorCarriedSpaceSidecar,
)
    source_carried_metadata = _carried_space_identity_metadata(source.carried_space)
    sidecar_carried_metadata = _carried_space_identity_metadata(sidecar.carried_space)
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
            source_carried_metadata.representation_final_dimension,
            sidecar_carried_metadata.representation_final_dimension,
        ),
        carried_parent_axis_counts = _record_field_comparison(
            :carried_parent_axis_counts,
            source_carried_metadata.parent_axis_counts,
            sidecar_carried_metadata.parent_axis_counts,
        ),
        carried_parent_dimension = _record_field_comparison(
            :carried_parent_dimension,
            source_carried_metadata.parent_dimension,
            sidecar_carried_metadata.parent_dimension,
        ),
        carried_representation_basis_kind = _record_field_comparison(
            :carried_representation_basis_kind,
            source_carried_metadata.representation_basis_kind,
            sidecar_carried_metadata.representation_basis_kind,
        ),
        carried_representation_parent_kind = _record_field_comparison(
            :carried_representation_parent_kind,
            source_carried_metadata.representation_parent_kind,
            sidecar_carried_metadata.representation_parent_kind,
        ),
        carried_representation_final_dimension = _record_field_comparison(
            :carried_representation_final_dimension,
            source_carried_metadata.representation_final_dimension,
            sidecar_carried_metadata.representation_final_dimension,
        ),
        carried_axis_sharing = _record_field_comparison(
            :carried_axis_sharing,
            source_carried_metadata.axis_sharing,
            sidecar_carried_metadata.axis_sharing,
        ),
        carried_has_contracted_parent = _record_field_comparison(
            :carried_has_contracted_parent,
            source_carried_metadata.has_contracted_parent,
            sidecar_carried_metadata.has_contracted_parent,
        ),
        carried_has_staged_sidecar = _record_field_comparison(
            :carried_has_staged_sidecar,
            source_carried_metadata.has_staged_sidecar,
            sidecar_carried_metadata.has_staged_sidecar,
        ),
        carried_staged_by_center_path = _record_field_comparison(
            :carried_staged_by_center_path,
            source_carried_metadata.staged_by_center_path,
            sidecar_carried_metadata.staged_by_center_path,
        ),
        carried_provenance_input_kind = _record_field_comparison(
            :carried_provenance_input_kind,
            source_carried_metadata.provenance_input_kind,
            sidecar_carried_metadata.provenance_input_kind,
        ),
        carried_provenance_source_type = _record_field_comparison(
            :carried_provenance_source_type,
            source_carried_metadata.provenance_source_type,
            sidecar_carried_metadata.provenance_source_type,
        ),
        carried_provenance_representation_basis_kind = _record_field_comparison(
            :carried_provenance_representation_basis_kind,
            source_carried_metadata.provenance_representation_basis_kind,
            sidecar_carried_metadata.provenance_representation_basis_kind,
        ),
        carried_provenance_route_metadata = _record_field_comparison(
            :carried_provenance_route_metadata,
            source_carried_metadata.provenance_route_metadata,
            sidecar_carried_metadata.provenance_route_metadata,
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
    source_carried_metadata = _carried_space_identity_metadata(source.carried_space)
    sidecar_carried_metadata = _carried_space_identity_metadata(sidecar.carried_space)
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
        source_carried_dimension = source_carried_metadata.representation_final_dimension,
        sidecar_carried_dimension = sidecar_carried_metadata.representation_final_dimension,
        source_parent_axis_counts = source_carried_metadata.parent_axis_counts,
        sidecar_parent_axis_counts = sidecar_carried_metadata.parent_axis_counts,
        source_parent_dimension = source_carried_metadata.parent_dimension,
        sidecar_parent_dimension = sidecar_carried_metadata.parent_dimension,
        intentionally_not_compared = (
            :parent_axis_object_identity,
            :parent_axis_tables,
            :coefficient_matrix_values,
            :support_index_values,
            :overlap_matrix_values,
            :one_body_matrix_values,
            :interaction_matrix_values,
            :residual_packet_values,
            :metric_packet_values,
        ),
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
metadata plus cheap carried-space identity metadata. It does not build
operators, dispatch kernels, construct metric packets, or mutate numerical
matrices.
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

function _receipt_diagnostics(
    source::CartesianOperatorBuildSource3D,
    operators::OrdinaryCartesianOperators3D,
    record::CartesianQWOperatorConstructionRecord3D,
    input_kind::Symbol,
    forwarded_keyword_names,
)
    record_diagnostics = qw_operator_construction_record_diagnostics(record)
    return (
        input_kind = input_kind,
        source_input_kind = operator_build_source_provenance(source).input_kind,
        builder = :ordinary_cartesian_qiu_white_operators,
        delegated_to_existing_builder = true,
        forwarded_keyword_names = forwarded_keyword_names,
        source_sidecar_agree = record_diagnostics.source_sidecar_agree,
        mismatch_fields = record_diagnostics.mismatch_fields,
        ambiguous_mismatch_fields = record_diagnostics.ambiguous_mismatch_fields,
        compared_fields = record_diagnostics.compared_fields,
        operator_dimension = size(operators.overlap, 1),
        operator_gausslet_count = operators.gausslet_count,
        operator_residual_count = operators.residual_count,
        gausslet_backend = operators.gausslet_backend,
        interaction_treatment = operators.interaction_treatment,
        nuclear_term_storage = operators.nuclear_term_storage,
        dense_parent_matrix_used = false,
        heavy_metric_packet_built = false,
        metric_packet_built = false,
        new_hamiltonian_kernel_used = false,
        numerical_outputs_changed = false,
        operator_built = true,
    )
end

function _receipt_provenance(
    source::CartesianOperatorBuildSource3D,
    operators::OrdinaryCartesianOperators3D,
    record::CartesianQWOperatorConstructionRecord3D,
    input_kind::Symbol,
    forwarded_keyword_names,
)
    return (
        source = :cartesian_qw_operator_construction_receipt,
        input_kind = input_kind,
        builder = :ordinary_cartesian_qiu_white_operators,
        build_source_provenance = operator_build_source_provenance(source),
        construction_record_provenance = qw_operator_construction_record_provenance(record),
        operator_type = nameof(typeof(operators)),
        basis_type = nameof(typeof(operators.basis)),
        gaussian_data_type =
            isnothing(operators.gaussian_data) ? nothing : nameof(typeof(operators.gaussian_data)),
        forwarded_keyword_names = forwarded_keyword_names,
    )
end

function _cartesian_qw_operator_construction_receipt(
    source::CartesianOperatorBuildSource3D,
    operators::OrdinaryCartesianOperators3D,
    input_kind::Symbol,
    forwarded_keyword_names,
)
    record = cartesian_qw_operator_construction_record(source, operators)
    diagnostics = _receipt_diagnostics(
        source,
        operators,
        record,
        input_kind,
        forwarded_keyword_names,
    )
    provenance = _receipt_provenance(
        source,
        operators,
        record,
        input_kind,
        forwarded_keyword_names,
    )
    return CartesianQWOperatorConstructionReceipt3D(
        source,
        operators,
        record,
        diagnostics,
        provenance,
    )
end

"""
    cartesian_qw_operator_construction_receipt(input; kwargs...)

Build an internal audited receipt for selected existing Cartesian QW routes.

This helper is a construction wrapper only: it first builds the normalized
pre-build source, delegates operator assembly to the matching existing
`ordinary_cartesian_qiu_white_operators` method, and records the source/operator
audit. It is not a new Hamiltonian builder and does not change backend routing
or numerical kernels.
"""
function cartesian_qw_operator_construction_receipt(
    basis::AbstractBondAlignedOrdinaryQWBasis3D;
    nuclear_charges::AbstractVector{<:Real} = basis.nuclear_charges,
    nuclear_term_storage::Symbol = :auto,
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    kwargs...,
)
    forwarded_keyword_names = (:nuclear_charges, :nuclear_term_storage,
        :interaction_treatment, :gausslet_backend, keys((; kwargs...))...)
    source = cartesian_qw_operator_build_source(
        basis;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
    )
    operators = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        kwargs...,
    )
    return _cartesian_qw_operator_construction_receipt(
        source,
        operators,
        :bond_aligned_direct_product_input,
        forwarded_keyword_names,
    )
end

function cartesian_qw_operator_construction_receipt(
    fixed_block::_NestedFixedBlock3D{<:AbstractBondAlignedOrdinaryQWBasis3D};
    nuclear_charges::AbstractVector{<:Real} = fixed_block.parent_basis.nuclear_charges,
    nuclear_term_storage::Symbol = :auto,
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    kwargs...,
)
    forwarded_keyword_names = (:nuclear_charges, :nuclear_term_storage,
        :interaction_treatment, :gausslet_backend, keys((; kwargs...))...)
    source = cartesian_qw_operator_build_source(
        fixed_block;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
    )
    operators = ordinary_cartesian_qiu_white_operators(
        fixed_block;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        kwargs...,
    )
    return _cartesian_qw_operator_construction_receipt(
        source,
        operators,
        :bond_aligned_nested_fixed_block_input,
        forwarded_keyword_names,
    )
end

function cartesian_qw_operator_construction_receipt(
    basis::MappedUniformBasis,
    gaussian_data::LegacyAtomicGaussianSupplement;
    Z::Real = 2.0,
    nuclear_term_storage::Symbol = :total_only,
    interaction_treatment::Symbol = :mwg,
    gausslet_backend::Symbol = :auto,
    residual_keep_policy::Symbol = :near_null_only,
    kwargs...,
)
    forwarded_keyword_names =
        (:Z, :nuclear_term_storage, :interaction_treatment, :gausslet_backend, :residual_keep_policy,
            keys((; kwargs...))...)
    source = cartesian_qw_operator_build_source(
        basis,
        gaussian_data;
        Z = Z,
        nuclear_term_storage = nuclear_term_storage,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
    )
    operators = ordinary_cartesian_qiu_white_operators(
        basis,
        gaussian_data;
        Z = Z,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        residual_keep_policy = residual_keep_policy,
        kwargs...,
    )
    return _cartesian_qw_operator_construction_receipt(
        source,
        operators,
        :atomic_direct_product_input,
        forwarded_keyword_names,
    )
end

function cartesian_qw_operator_construction_receipt(
    fixed_block::_NestedFixedBlock3D,
    gaussian_data::LegacyAtomicGaussianSupplement;
    Z::Real = 2.0,
    nuclear_term_storage::Symbol = :total_only,
    interaction_treatment::Symbol = :mwg,
    gausslet_backend::Symbol = :auto,
    residual_keep_policy::Symbol = :near_null_only,
    kwargs...,
)
    forwarded_keyword_names =
        (:Z, :nuclear_term_storage, :interaction_treatment, :gausslet_backend, :residual_keep_policy,
            keys((; kwargs...))...)
    source = cartesian_qw_operator_build_source(
        fixed_block,
        gaussian_data;
        Z = Z,
        nuclear_term_storage = nuclear_term_storage,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
    )
    operators = ordinary_cartesian_qiu_white_operators(
        fixed_block,
        gaussian_data;
        Z = Z,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        residual_keep_policy = residual_keep_policy,
        kwargs...,
    )
    return _cartesian_qw_operator_construction_receipt(
        source,
        operators,
        :atomic_nested_fixed_block_input,
        forwarded_keyword_names,
    )
end

function cartesian_qw_operator_construction_receipt(
    basis::BondAlignedDiatomicQWBasis3D,
    gaussian_data::Union{
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    };
    nuclear_charges::AbstractVector{<:Real} = basis.nuclear_charges,
    nuclear_term_storage::Symbol = :auto,
    interaction_treatment::Symbol = :mwg,
    gausslet_backend::Symbol = :numerical_reference,
    kwargs...,
)
    forwarded_keyword_names = (:nuclear_charges, :nuclear_term_storage,
        :interaction_treatment, :gausslet_backend, keys((; kwargs...))...)
    source = cartesian_qw_operator_build_source(
        basis,
        gaussian_data;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
    )
    operators = ordinary_cartesian_qiu_white_operators(
        basis,
        gaussian_data;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        kwargs...,
    )
    return _cartesian_qw_operator_construction_receipt(
        source,
        operators,
        :bond_aligned_molecular_direct_product_input,
        forwarded_keyword_names,
    )
end

function cartesian_qw_operator_construction_receipt(
    fixed_block::_NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D},
    gaussian_data::Union{
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    };
    nuclear_charges::AbstractVector{<:Real} = fixed_block.parent_basis.nuclear_charges,
    nuclear_term_storage::Symbol = :auto,
    interaction_treatment::Symbol = :mwg,
    gausslet_backend::Symbol = :numerical_reference,
    kwargs...,
)
    forwarded_keyword_names = (:nuclear_charges, :nuclear_term_storage,
        :interaction_treatment, :gausslet_backend, keys((; kwargs...))...)
    source = cartesian_qw_operator_build_source(
        fixed_block,
        gaussian_data;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
    )
    operators = ordinary_cartesian_qiu_white_operators(
        fixed_block,
        gaussian_data;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        kwargs...,
    )
    return _cartesian_qw_operator_construction_receipt(
        source,
        operators,
        :bond_aligned_molecular_nested_fixed_block_input,
        forwarded_keyword_names,
    )
end

function _nested_q_row_retained_counts(region_builds)
    return Tuple(build.retained_count for build in region_builds)
end

function _nested_q_row_shared_retained_counts(region_builds)
    return Tuple(
        build.retained_count for build in region_builds
        if build.role == :regular_shared_molecular_shell
    )
end

function _nested_q_row_shared_retained_count(region_builds)
    shared_counts = _nested_q_row_shared_retained_counts(region_builds)
    return length(shared_counts) == 1 ? only(shared_counts) : nothing
end

function _nested_q_row_region_roles(region_builds)
    return Tuple(build.role for build in region_builds)
end

"""
    _nested_bond_aligned_diatomic_high_order_q_row_route_receipt(basis; shared_q, ...)

Internal q-row receipt for the experimental high-order diatomic
atom-growth/endcap-panel path. The helper only composes existing recipe,
fixed-block, and QW receipt helpers; it does not add a public frontend, change
defaults, or implement Hamiltonian kernels.
"""
function _nested_bond_aligned_diatomic_high_order_q_row_route_receipt(
    basis::BondAlignedDiatomicQWBasis3D;
    shared_q::Integer,
    shared_order::Integer = shared_q,
    protected_atom_side_count::Integer = 5,
    q_min::Integer = 4,
    nside::Integer = 5,
    expansion = coulomb_gaussian_expansion(doacc = false),
    term_coefficients = nothing,
    packet_kernel::Symbol = :factorized_direct,
    nuclear_charges::AbstractVector{<:Real} = basis.nuclear_charges,
    nuclear_term_storage::Symbol = :total_only,
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :pgdg_localized_experimental,
)
    gausslet_backend == :pgdg_localized_experimental || throw(
        ArgumentError(
            "experimental high-order q-row route receipt requires gausslet_backend = :pgdg_localized_experimental",
        ),
    )
    q_min_int = Int(q_min)
    q_min_int == 4 || throw(
        ArgumentError("experimental high-order q-row route receipt currently requires q_min = 4"),
    )
    shared_q_int = Int(shared_q)
    shared_order_int = Int(shared_order)
    nside_int = Int(nside)
    protected_atom_side_count_int = Int(protected_atom_side_count)
    coefficients =
        isnothing(term_coefficients) ? Float64.(expansion.coefficients) :
        Float64.(term_coefficients)

    axis_bundles = _qwrg_bond_aligned_axis_bundles(basis, expansion)
    policy = _nested_bond_aligned_diatomic_high_order_recipe_policy(
        basis,
        axis_bundles;
        protected_atom_side_count = protected_atom_side_count_int,
        q_min = q_min_int,
        atom_q = 4,
        atom_order = 4,
        shared_q = shared_q_int,
        shared_order = shared_order_int,
        contact_q = 4,
        contact_order = 4,
        outer_mismatch_q = 4,
        outer_mismatch_order = 4,
    )
    policy_diagnostics = _nested_bond_aligned_diatomic_high_order_recipe_policy_diagnostics(
        policy,
    )
    realization_audit = _nested_bond_aligned_diatomic_high_order_recipe_realization_audit(
        policy,
    )
    realization_diagnostics =
        _nested_bond_aligned_diatomic_high_order_recipe_realization_diagnostics(
            realization_audit,
        )
    source_construction =
        _nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            axis_bundles,
            policy;
            nside = nside_int,
            term_coefficients = coefficients,
            packet_kernel = packet_kernel,
            build_sequence_packet = true,
        )
    source_diagnostics =
        _nested_bond_aligned_diatomic_high_order_recipe_source_construction_diagnostics(
            source_construction,
        )
    readiness = _nested_bond_aligned_diatomic_high_order_recipe_source_readiness(
        source_construction;
        build_fixed_block = true,
    )
    readiness_diagnostics =
        _nested_bond_aligned_diatomic_high_order_recipe_source_readiness_diagnostics(
            readiness,
        )
    fixed_block = readiness.fixed_block
    isnothing(fixed_block) && throw(
        ArgumentError("experimental high-order q-row route receipt did not produce a fixed block"),
    )

    # Keep the existing helper as the fixed-block contract oracle without
    # rebuilding the block: readiness used it when `build_fixed_block=true`.
    readiness_diagnostics.can_produce_fixed_block || throw(
        ArgumentError(
            "experimental high-order q-row route receipt is not fixed-block ready",
        ),
    )

    qw_receipt = cartesian_qw_operator_construction_receipt(
        fixed_block;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        expansion = expansion,
    )
    qw_receipt_diagnostics = qw_operator_construction_receipt_diagnostics(qw_receipt)
    qw_record_diagnostics = qw_operator_construction_record_diagnostics(
        qw_operator_construction_receipt_record(qw_receipt),
    )
    operators = qw_operator_construction_receipt_operators(qw_receipt)
    region_builds = source_diagnostics.region_builds
    diagnostics = (
        route_label = :bond_aligned_diatomic_high_order_q_row_route,
        receipt_contract = :delegate_existing_recipe_fixed_block_and_qw_receipt,
        shared_q = shared_q_int,
        shared_order = shared_order_int,
        q_min = q_min_int,
        protected_atom_side_count = protected_atom_side_count_int,
        nside = nside_int,
        non_shared_q_policy = get(source_diagnostics.metadata, :non_shared_q_policy, :unknown),
        parent_dimension = readiness_diagnostics.parent_dimension,
        fixed_dimension = readiness_diagnostics.fixed_dimension,
        retained_counts_by_region = _nested_q_row_retained_counts(region_builds),
        shared_retained_counts = _nested_q_row_shared_retained_counts(region_builds),
        shared_retained_count = _nested_q_row_shared_retained_count(region_builds),
        region_roles = _nested_q_row_region_roles(region_builds),
        support_coverage = source_diagnostics.support_coverage,
        overlap_error = readiness_diagnostics.overlap_error,
        staged_sidecar_available =
            readiness_diagnostics.staged_by_center_sidecar_available,
        backend = qw_receipt_diagnostics.gausslet_backend,
        interaction_treatment = qw_receipt_diagnostics.interaction_treatment,
        nuclear_term_storage = qw_receipt_diagnostics.nuclear_term_storage,
        residual_count = operators.residual_count,
        gausslet_count = operators.gausslet_count,
        source_sidecar_agree = qw_receipt_diagnostics.source_sidecar_agree &&
                               qw_record_diagnostics.source_sidecar_agree,
        mismatch_fields = Tuple(qw_receipt_diagnostics.mismatch_fields),
        dense_parent_matrix_used = qw_receipt_diagnostics.dense_parent_matrix_used,
        heavy_metric_packet_built = qw_receipt_diagnostics.heavy_metric_packet_built,
        new_hamiltonian_kernel_used = qw_receipt_diagnostics.new_hamiltonian_kernel_used,
        numerical_outputs_changed = qw_receipt_diagnostics.numerical_outputs_changed,
        default_source_builder_changed =
            get(source_diagnostics.metadata, :default_source_builder_changed, false),
        active_builder_consumes = source_diagnostics.active_builder_consumes,
    )
    provenance = (
        source = :_nested_bond_aligned_diatomic_high_order_q_row_route_receipt,
        basis_type = nameof(typeof(basis)),
        q_row_scope = :shared_endcap_panel_only,
        non_shared_q_policy = :fixed_q4_order4,
        source_builder =
            :_nested_bond_aligned_diatomic_high_order_recipe_source_construction,
        fixed_block_builder =
            :_nested_bond_aligned_diatomic_high_order_recipe_source_fixed_block,
        qw_builder = :cartesian_qw_operator_construction_receipt,
        public_api = false,
        science_validation = false,
    )

    return _BondAlignedDiatomicHighOrderQRowRouteReceipt3D(
        basis,
        expansion,
        axis_bundles,
        policy,
        policy_diagnostics,
        realization_audit,
        realization_diagnostics,
        source_construction,
        source_diagnostics,
        readiness,
        readiness_diagnostics,
        fixed_block,
        qw_receipt,
        qw_receipt_diagnostics,
        qw_record_diagnostics,
        diagnostics,
        provenance,
    )
end

"""
    _nested_bond_aligned_homonuclear_high_order_q_row_fixture_receipt(; ...)

Private homonuclear fixture wrapper for the experimental q-row route. The
helper requires geometry-defining inputs, constructs the existing
`BondAlignedDiatomicQWBasis3D` with `bond_aligned_homonuclear_qw_basis`, and
delegates unchanged to
`_nested_bond_aligned_diatomic_high_order_q_row_route_receipt`.
"""
function _nested_bond_aligned_homonuclear_high_order_q_row_fixture_receipt(;
    bond_length,
    core_spacing,
    xmax_parallel,
    xmax_transverse,
    shared_q::Integer,
    family = :G10,
    bond_axis::Symbol = :z,
    nuclear_charge::Real = 1.0,
    reference_spacing::Real = 1.0,
    tail_spacing::Real = 10.0,
    shared_order::Integer = shared_q,
    protected_atom_side_count::Integer = 5,
    q_min::Integer = 4,
    nside::Integer = 5,
    expansion = coulomb_gaussian_expansion(doacc = false),
    packet_kernel::Symbol = :factorized_direct,
    nuclear_term_storage::Symbol = :total_only,
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :pgdg_localized_experimental,
)
    basis = bond_aligned_homonuclear_qw_basis(
        family = family,
        bond_length = bond_length,
        core_spacing = core_spacing,
        xmax_parallel = xmax_parallel,
        xmax_transverse = xmax_transverse,
        bond_axis = bond_axis,
        nuclear_charge = nuclear_charge,
        reference_spacing = reference_spacing,
        tail_spacing = tail_spacing,
    )
    route_receipt = _nested_bond_aligned_diatomic_high_order_q_row_route_receipt(
        basis;
        shared_q = shared_q,
        shared_order = shared_order,
        protected_atom_side_count = protected_atom_side_count,
        q_min = q_min,
        nside = nside,
        expansion = expansion,
        packet_kernel = packet_kernel,
        nuclear_charges = basis.nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
    )
    route_diagnostics =
        _nested_bond_aligned_diatomic_high_order_q_row_route_diagnostics(
            route_receipt,
        )
    route_provenance =
        _nested_bond_aligned_diatomic_high_order_q_row_route_provenance(
            route_receipt,
        )
    parent_axis_counts = (
        length(basis.basis_x),
        length(basis.basis_y),
        length(basis.basis_z),
    )
    parent_dimension = prod(parent_axis_counts)
    flat_index_convention = (
        one_based = true,
        formula = "flat = (ix - 1) * ny * nz + (iy - 1) * nz + iz",
        fastest_axis = :z,
        slowest_axis = :x,
    )
    diagnostics = (
        route_label = :bond_aligned_homonuclear_high_order_q_row_fixture,
        receipt_contract = :construct_homonuclear_basis_then_delegate_q_row_route,
        basis_constructor = :bond_aligned_homonuclear_qw_basis,
        family = family,
        bond_length = Float64(bond_length),
        core_spacing = Float64(core_spacing),
        xmax_parallel = Float64(xmax_parallel),
        xmax_transverse = Float64(xmax_transverse),
        bond_axis = bond_axis,
        reference_spacing = Float64(reference_spacing),
        tail_spacing = Float64(tail_spacing),
        nuclear_charge = Float64(nuclear_charge),
        basis_nuclear_charges = Tuple(Float64.(basis.nuclear_charges)),
        basis_nuclei = Tuple(basis.nuclei),
        parent_axis_counts = parent_axis_counts,
        parent_dimension = parent_dimension,
        flat_index_convention = flat_index_convention,
        shared_q = route_diagnostics.shared_q,
        shared_order = route_diagnostics.shared_order,
        q_min = route_diagnostics.q_min,
        protected_atom_side_count = route_diagnostics.protected_atom_side_count,
        nside = route_diagnostics.nside,
        fixed_dimension = route_diagnostics.fixed_dimension,
        retained_counts_by_region = route_diagnostics.retained_counts_by_region,
        shared_retained_counts = route_diagnostics.shared_retained_counts,
        shared_retained_count = route_diagnostics.shared_retained_count,
        q_row_route_diagnostics = route_diagnostics,
    )
    provenance = (
        source = :_nested_bond_aligned_homonuclear_high_order_q_row_fixture_receipt,
        basis_constructor = :bond_aligned_homonuclear_qw_basis,
        delegated_route = route_provenance.source,
        charge_policy = :basis_nuclear_charges_only,
        homonuclear_only = true,
        heteronuclear_support = false,
        element_label_provenance = false,
        public_api = false,
        science_validation = false,
    )
    return _BondAlignedHomonuclearHighOrderQRowFixtureReceipt3D(
        basis,
        route_receipt,
        diagnostics,
        provenance,
    )
end

"""
    _nested_bond_aligned_homonuclear_high_order_q_row_fixture_supplement_receipt(; ...)

Private supplement-aware homonuclear q-row fixture wrapper. It constructs the
q-row fixed block through
`_nested_bond_aligned_homonuclear_high_order_q_row_fixture_receipt`, constructs
a legacy homonuclear bond-aligned molecular supplement on the same nuclei, and
delegates QW construction to the existing nested molecular supplement receipt.
"""
function _nested_bond_aligned_homonuclear_high_order_q_row_fixture_supplement_receipt(;
    bond_length,
    core_spacing,
    xmax_parallel,
    xmax_transverse,
    shared_q::Integer,
    atom::AbstractString,
    basis_name::AbstractString,
    family = :G10,
    bond_axis::Symbol = :z,
    nuclear_charge::Real = 1.0,
    reference_spacing::Real = 1.0,
    tail_spacing::Real = 10.0,
    shared_order::Integer = shared_q,
    protected_atom_side_count::Integer = 5,
    q_min::Integer = 4,
    nside::Integer = 5,
    expansion = coulomb_gaussian_expansion(doacc = false),
    packet_kernel::Symbol = :factorized_direct,
    lmax::Integer = 0,
    basisfile::Union{Nothing,AbstractString} = nothing,
    max_width::Union{Nothing,Real} = nothing,
    nuclear_term_storage::Symbol = :by_center,
    interaction_treatment::Symbol = :mwg,
    gausslet_backend::Symbol = :pgdg_localized_experimental,
)
    gausslet_backend == :pgdg_localized_experimental || throw(
        ArgumentError(
            "experimental high-order q-row fixture supplement receipt requires gausslet_backend = :pgdg_localized_experimental",
        ),
    )
    fixture_receipt =
        _nested_bond_aligned_homonuclear_high_order_q_row_fixture_receipt(
            bond_length = bond_length,
            core_spacing = core_spacing,
            xmax_parallel = xmax_parallel,
            xmax_transverse = xmax_transverse,
            shared_q = shared_q,
            family = family,
            bond_axis = bond_axis,
            nuclear_charge = nuclear_charge,
            reference_spacing = reference_spacing,
            tail_spacing = tail_spacing,
            shared_order = shared_order,
            protected_atom_side_count = protected_atom_side_count,
            q_min = q_min,
            nside = nside,
            expansion = expansion,
            packet_kernel = packet_kernel,
            gausslet_backend = gausslet_backend,
        )
    basis = fixture_receipt.basis
    fixed_block = fixture_receipt.q_row_route_receipt.fixed_block
    supplement = legacy_bond_aligned_diatomic_gaussian_supplement(
        atom,
        basis_name,
        basis.nuclei;
        lmax = lmax,
        basisfile = basisfile,
        max_width = max_width,
    )
    supplement_qw_receipt = cartesian_qw_operator_construction_receipt(
        fixed_block,
        supplement;
        nuclear_charges = fixed_block.parent_basis.nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        expansion = expansion,
    )
    fixture_diagnostics =
        _nested_bond_aligned_homonuclear_high_order_q_row_fixture_diagnostics(
            fixture_receipt,
        )
    fixture_provenance =
        _nested_bond_aligned_homonuclear_high_order_q_row_fixture_provenance(
            fixture_receipt,
        )
    supplement_diagnostics = qw_operator_construction_receipt_diagnostics(
        supplement_qw_receipt,
    )
    supplement_record_diagnostics = qw_operator_construction_record_diagnostics(
        qw_operator_construction_receipt_record(supplement_qw_receipt),
    )
    operators = qw_operator_construction_receipt_operators(supplement_qw_receipt)
    diagnostics = (
        route_label = :bond_aligned_homonuclear_high_order_q_row_fixture_supplement,
        receipt_contract =
            :construct_homonuclear_q_row_fixture_then_delegate_nested_molecular_supplement_receipt,
        fixture_diagnostics = fixture_diagnostics,
        supplement_constructor = :legacy_bond_aligned_diatomic_gaussian_supplement,
        supplement_atom = String(atom),
        supplement_basis_name = String(basis_name),
        supplement_lmax = Int(lmax),
        supplement_basisfile = basisfile,
        supplement_max_width = max_width === nothing ? nothing : Float64(max_width),
        supplement_nuclei = Tuple(supplement.nuclei),
        nuclear_charges = Tuple(Float64.(fixed_block.parent_basis.nuclear_charges)),
        backend = supplement_diagnostics.gausslet_backend,
        interaction_treatment = supplement_diagnostics.interaction_treatment,
        nuclear_term_storage = supplement_diagnostics.nuclear_term_storage,
        gausslet_count = operators.gausslet_count,
        residual_count = operators.residual_count,
        final_dimension = size(operators.overlap, 1),
        source_sidecar_agree = supplement_diagnostics.source_sidecar_agree &&
                               supplement_record_diagnostics.source_sidecar_agree,
        mismatch_fields = Tuple(supplement_diagnostics.mismatch_fields),
        dense_parent_matrix_used = supplement_diagnostics.dense_parent_matrix_used,
        heavy_metric_packet_built = supplement_diagnostics.heavy_metric_packet_built,
        new_hamiltonian_kernel_used =
            supplement_diagnostics.new_hamiltonian_kernel_used,
        numerical_outputs_changed = supplement_diagnostics.numerical_outputs_changed,
        supplement_qw_receipt_diagnostics = supplement_diagnostics,
    )
    provenance = (
        source =
            :_nested_bond_aligned_homonuclear_high_order_q_row_fixture_supplement_receipt,
        fixture_source = fixture_provenance.source,
        supplement_constructor = :legacy_bond_aligned_diatomic_gaussian_supplement,
        qw_builder = :cartesian_qw_operator_construction_receipt,
        charge_policy = :fixed_block_parent_basis_nuclear_charges,
        homonuclear_only = true,
        heteronuclear_support = false,
        public_api = false,
        science_validation = false,
    )
    return _BondAlignedHomonuclearHighOrderQRowFixtureSupplementReceipt3D(
        fixture_receipt,
        supplement,
        supplement_qw_receipt,
        diagnostics,
        provenance,
    )
end

function _nested_capture_target_matrix(
    target_coefficients::AbstractMatrix{<:Real},
)
    matrix = Matrix{Float64}(target_coefficients)
    all(isfinite, matrix) || throw(
        ArgumentError("q-row fixture supplement capture/H1 target coefficients must be finite"),
    )
    return matrix
end

function _nested_capture_target_labels(labels, ntarget::Int)
    labels === nothing && return ["target_$(index)" for index in 1:ntarget]
    length(labels) == ntarget || throw(
        DimensionMismatch(
            "target label count $(length(labels)) does not match target column count $(ntarget)",
        ),
    )
    return String[string(label) for label in labels]
end

function _nested_capture_target_occupations(occupations, ntarget::Int)
    occupations === nothing && return nothing
    length(occupations) == ntarget || throw(
        DimensionMismatch(
            "target occupation count $(length(occupations)) does not match target column count $(ntarget)",
        ),
    )
    values = Float64[Float64(occupation) for occupation in occupations]
    all(isfinite, values) || throw(ArgumentError("target occupations must be finite"))
    all(>=(0.0), values) || throw(ArgumentError("target occupations must be nonnegative"))
    return values
end

function _nested_capture_column_gram_values(
    gram::AbstractMatrix{<:Real},
    label::AbstractString,
)
    values = Float64[Float64(value) for value in real.(diag(gram))]
    all(isfinite, values) || throw(ArgumentError("$(label) diagonal values must be finite"))
    return values
end

function _nested_capture_positive_source_norms(source_norms::AbstractVector{<:Real})
    all(>(0.0), source_norms) || throw(
        ArgumentError(
            "q-row fixture supplement capture/H1 requires every target column to have positive parent-overlap norm",
        ),
    )
    return source_norms
end

function _nested_capture_expectations(
    coefficients::AbstractMatrix{<:Real},
    matrix::AbstractMatrix{<:Real},
    norms::AbstractVector{<:Real},
)
    projected = transpose(coefficients) * Matrix{Float64}(matrix) * coefficients
    values = _nested_capture_column_gram_values(projected, "projected expectation")
    return Float64[values[index] / Float64(norms[index]) for index in eachindex(values)]
end

function _nested_capture_rows(
    labels::AbstractVector{<:AbstractString},
    occupations,
    source_norms::AbstractVector{<:Real},
    fixed_norms::AbstractVector{<:Real},
    final_norms::AbstractVector{<:Real},
    parent_h1::AbstractVector{<:Real},
    fixed_h1::AbstractVector{<:Real},
    final_h1::AbstractVector{<:Real},
)
    rows = NamedTuple[]
    for index in eachindex(labels)
        source_norm = Float64(source_norms[index])
        fixed_norm = Float64(fixed_norms[index])
        final_norm = Float64(final_norms[index])
        fixed_fraction = fixed_norm / source_norm
        final_fraction = final_norm / source_norm
        parent_value = Float64(parent_h1[index])
        fixed_value = Float64(fixed_h1[index])
        final_value = Float64(final_h1[index])
        push!(
            rows,
            (
                column = index,
                label = labels[index],
                occupation = occupations === nothing ? nothing : occupations[index],
                source_norm = source_norm,
                fixed_captured_norm = fixed_norm,
                final_captured_norm = final_norm,
                fixed_capture_fraction = fixed_fraction,
                final_capture_fraction = final_fraction,
                final_minus_fixed_capture = final_fraction - fixed_fraction,
                parent_h1_expectation = parent_value,
                fixed_projected_h1_expectation = fixed_value,
                final_projected_h1_expectation = final_value,
                fixed_h1_delta = fixed_value - parent_value,
                final_h1_delta = final_value - parent_value,
            ),
        )
    end
    return rows
end

function _nested_capture_summary(rows::AbstractVector, occupations)
    weights = occupations === nothing ? ones(Float64, length(rows)) : Float64.(occupations)
    source_total = sum(weights[index] * rows[index].source_norm for index in eachindex(rows))
    fixed_total = sum(weights[index] * rows[index].fixed_captured_norm for index in eachindex(rows))
    final_total = sum(weights[index] * rows[index].final_captured_norm for index in eachindex(rows))
    source_total > 0.0 || throw(
        ArgumentError("q-row fixture supplement capture/H1 summary requires positive total source norm"),
    )
    fixed_fractions = Float64[row.fixed_capture_fraction for row in rows]
    final_fractions = Float64[row.final_capture_fraction for row in rows]
    worst_fixed_index = argmin(fixed_fractions)
    worst_final_index = argmin(final_fractions)
    max_fixed_h1_delta = maximum(abs(row.fixed_h1_delta) for row in rows)
    max_final_h1_delta = maximum(abs(row.final_h1_delta) for row in rows)
    return (
        target_count = length(rows),
        occupation_policy = occupations === nothing ? :unit_weights : :explicit_occupations,
        source_norm_total = source_total,
        fixed_captured_norm_total = fixed_total,
        final_captured_norm_total = final_total,
        fixed_capture_fraction_total = fixed_total / source_total,
        final_capture_fraction_total = final_total / source_total,
        final_minus_fixed_capture_total = (final_total - fixed_total) / source_total,
        worst_fixed_label = rows[worst_fixed_index].label,
        worst_fixed_capture = rows[worst_fixed_index].fixed_capture_fraction,
        worst_final_label = rows[worst_final_index].label,
        worst_final_capture = rows[worst_final_index].final_capture_fraction,
        max_abs_fixed_h1_delta = max_fixed_h1_delta,
        max_abs_final_h1_delta = max_final_h1_delta,
    )
end

"""
    _nested_bond_aligned_homonuclear_high_order_q_row_fixture_supplement_capture_h1(target_coefficients; ...)

Private parent-grid target capture/H1 diagnostic for the homonuclear q-row
fixture supplement route. The target columns must be expressed in the exact
parent Cartesian grid of the fixture. Final-basis self-overlap is reported as a
trust gate only; the final projection uses the existing orthonormal final-basis
contract.
"""
function _nested_bond_aligned_homonuclear_high_order_q_row_fixture_supplement_capture_h1(
    target_coefficients::AbstractMatrix{<:Real};
    target_labels = nothing,
    target_occupations = nothing,
    bond_length,
    core_spacing,
    xmax_parallel,
    xmax_transverse,
    shared_q::Integer,
    atom::AbstractString,
    basis_name::AbstractString,
    family = :G10,
    bond_axis::Symbol = :z,
    nuclear_charge::Real = 1.0,
    reference_spacing::Real = 1.0,
    tail_spacing::Real = 10.0,
    shared_order::Integer = shared_q,
    protected_atom_side_count::Integer = 5,
    q_min::Integer = 4,
    nside::Integer = 5,
    expansion = coulomb_gaussian_expansion(doacc = false),
    packet_kernel::Symbol = :factorized_direct,
    lmax::Integer = 0,
    basisfile::Union{Nothing,AbstractString} = nothing,
    max_width::Union{Nothing,Real} = nothing,
    nuclear_term_storage::Symbol = :by_center,
    interaction_treatment::Symbol = :mwg,
    gausslet_backend::Symbol = :pgdg_localized_experimental,
)
    target_matrix = _nested_capture_target_matrix(target_coefficients)
    gausslet_backend == :pgdg_localized_experimental || throw(
        ArgumentError(
            "experimental q-row fixture supplement capture/H1 requires gausslet_backend = :pgdg_localized_experimental",
        ),
    )
    route_receipt =
        _nested_bond_aligned_homonuclear_high_order_q_row_fixture_supplement_receipt(
            bond_length = bond_length,
            core_spacing = core_spacing,
            xmax_parallel = xmax_parallel,
            xmax_transverse = xmax_transverse,
            shared_q = shared_q,
            atom = atom,
            basis_name = basis_name,
            family = family,
            bond_axis = bond_axis,
            nuclear_charge = nuclear_charge,
            reference_spacing = reference_spacing,
            tail_spacing = tail_spacing,
            shared_order = shared_order,
            protected_atom_side_count = protected_atom_side_count,
            q_min = q_min,
            nside = nside,
            expansion = expansion,
            packet_kernel = packet_kernel,
            lmax = lmax,
            basisfile = basisfile,
            max_width = max_width,
            nuclear_term_storage = nuclear_term_storage,
            interaction_treatment = interaction_treatment,
            gausslet_backend = gausslet_backend,
        )
    fixture_receipt = route_receipt.fixture_receipt
    fixed_block = fixture_receipt.q_row_route_receipt.fixed_block
    basis = fixture_receipt.basis
    fixture_diagnostics =
        _nested_bond_aligned_homonuclear_high_order_q_row_fixture_diagnostics(
            fixture_receipt,
        )
    parent_dimension = fixture_diagnostics.parent_dimension
    size(target_matrix, 1) == parent_dimension || throw(
        DimensionMismatch(
            "q-row fixture supplement capture/H1 target row count $(size(target_matrix, 1)) does not match parent dimension $(parent_dimension)",
        ),
    )
    ntarget = size(target_matrix, 2)
    ntarget > 0 || throw(
        ArgumentError("q-row fixture supplement capture/H1 requires at least one target column"),
    )
    labels = _nested_capture_target_labels(target_labels, ntarget)
    occupations = _nested_capture_target_occupations(target_occupations, ntarget)

    parent_qw_receipt = cartesian_qw_operator_construction_receipt(
        basis;
        nuclear_charges = fixed_block.parent_basis.nuclear_charges,
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = gausslet_backend,
        expansion = expansion,
    )
    parent_ops = qw_operator_construction_receipt_operators(parent_qw_receipt)
    fixed_ops = qw_operator_construction_receipt_operators(
        fixture_receipt.q_row_route_receipt.qw_receipt,
    )
    final_ops = qw_operator_construction_receipt_operators(
        route_receipt.supplement_qw_receipt,
    )

    parent_overlap = Matrix{Float64}(parent_ops.overlap)
    parent_hamiltonian = Matrix{Float64}(parent_ops.one_body_hamiltonian)
    fixed_overlap = Matrix{Float64}(fixed_block.overlap)
    fixed_hamiltonian = Matrix{Float64}(fixed_ops.one_body_hamiltonian)
    final_overlap = Matrix{Float64}(final_ops.overlap)
    final_hamiltonian = Matrix{Float64}(final_ops.one_body_hamiltonian)
    fixed_coefficients = Matrix{Float64}(fixed_block.coefficient_matrix)
    parent_target_overlap = parent_overlap * target_matrix
    source_gram = transpose(target_matrix) * parent_target_overlap
    source_norms = _nested_capture_positive_source_norms(
        _nested_capture_column_gram_values(source_gram, "parent target source norm"),
    )

    fixed_target_overlap = transpose(fixed_coefficients) * parent_target_overlap
    fixed_projected = Symmetric(fixed_overlap) \ fixed_target_overlap
    fixed_capture_gram =
        transpose(fixed_projected) * fixed_overlap * fixed_projected
    fixed_norms = _nested_capture_column_gram_values(
        fixed_capture_gram,
        "fixed capture norm",
    )

    supplement_parent_overlap = gto_overlap_matrix(basis, route_receipt.supplement)
    supplement_target_overlap = transpose(supplement_parent_overlap) * target_matrix
    raw_target_overlap = vcat(fixed_target_overlap, supplement_target_overlap)
    size(raw_target_overlap, 1) == size(final_ops.raw_to_final, 1) || throw(
        DimensionMismatch(
            "q-row fixture supplement capture/H1 raw target overlap rows do not match raw_to_final rows",
        ),
    )
    final_projected =
        Matrix{Float64}(transpose(final_ops.raw_to_final) * raw_target_overlap)
    final_capture_gram = transpose(final_projected) * final_projected
    final_norms = _nested_capture_column_gram_values(
        final_capture_gram,
        "final capture norm",
    )

    parent_h1 = _nested_capture_expectations(
        target_matrix,
        parent_hamiltonian,
        source_norms,
    )
    fixed_h1 = _nested_capture_expectations(
        fixed_projected,
        fixed_hamiltonian,
        fixed_norms,
    )
    final_h1 = _nested_capture_expectations(
        final_projected,
        final_hamiltonian,
        final_norms,
    )
    rows = _nested_capture_rows(
        labels,
        occupations,
        source_norms,
        fixed_norms,
        final_norms,
        parent_h1,
        fixed_h1,
        final_h1,
    )
    summary = _nested_capture_summary(rows, occupations)
    route_diagnostics =
        _nested_bond_aligned_homonuclear_high_order_q_row_fixture_supplement_diagnostics(
            route_receipt,
        )
    parent_receipt_diagnostics = qw_operator_construction_receipt_diagnostics(
        parent_qw_receipt,
    )
    diagnostics = (
        route_label = :bond_aligned_homonuclear_high_order_q_row_fixture_supplement_capture_h1,
        diagnostic_contract = :parent_grid_target_capture_h1,
        target_space = :parent_grid_coefficients,
        target_shape = size(target_matrix),
        parent_dimension = parent_dimension,
        fixed_dimension = route_diagnostics.fixture_diagnostics.fixed_dimension,
        final_dimension = route_diagnostics.final_dimension,
        residual_count = route_diagnostics.residual_count,
        source_gram = Matrix{Float64}(source_gram),
        fixed_capture_gram = Matrix{Float64}(fixed_capture_gram),
        final_capture_gram = Matrix{Float64}(final_capture_gram),
        fixed_overlap_error = norm(fixed_overlap - I, Inf),
        final_overlap_error = norm(final_overlap - I, Inf),
        final_basis_policy = :assume_orthonormal_final_basis,
        final_overlap_usage = :diagnostic_trust_gate_only,
        fixed_projection_formula = :solve_fixed_overlap_against_parent_target_overlap,
        final_projection_formula = :raw_target_overlap_times_raw_to_final,
        route_diagnostics = route_diagnostics,
        parent_receipt_diagnostics = parent_receipt_diagnostics,
        backend = route_diagnostics.backend,
        parent_backend = parent_ops.gausslet_backend,
        source_sidecar_agree = route_diagnostics.source_sidecar_agree &&
                               parent_receipt_diagnostics.source_sidecar_agree,
        mismatch_fields = Tuple(route_diagnostics.mismatch_fields),
        dense_parent_matrix_used = route_diagnostics.dense_parent_matrix_used ||
                                   parent_receipt_diagnostics.dense_parent_matrix_used,
        heavy_metric_packet_built = route_diagnostics.heavy_metric_packet_built ||
                                    parent_receipt_diagnostics.heavy_metric_packet_built,
        new_hamiltonian_kernel_used =
            route_diagnostics.new_hamiltonian_kernel_used ||
            parent_receipt_diagnostics.new_hamiltonian_kernel_used,
        numerical_outputs_changed =
            route_diagnostics.numerical_outputs_changed ||
            parent_receipt_diagnostics.numerical_outputs_changed,
        omitted_diagnostics = (
            spin_summary = :not_in_first_pass,
            tsv_parsing = :not_in_first_pass,
            normalized_capture_singular_values = :not_in_first_pass,
        ),
    )
    provenance = (
        source =
            :_nested_bond_aligned_homonuclear_high_order_q_row_fixture_supplement_capture_h1,
        route_receipt = route_receipt.provenance.source,
        target_contract = :exact_parent_grid_coefficients,
        public_api = false,
        science_validation = false,
        cr2_file_parser = false,
    )
    return _BondAlignedHomonuclearHighOrderQRowFixtureSupplementCaptureH1Diagnostic3D(
        route_receipt,
        parent_qw_receipt,
        target_matrix,
        Matrix{Float64}(fixed_projected),
        final_projected,
        rows,
        summary,
        diagnostics,
        provenance,
    )
end

function _cartesian_qw_operator_receipt_supported_input_summary()
    coverage = cartesian_qw_operator_receipt_coverage()
    return join((row.inputs for row in coverage.covered_route_families), "; ")
end

function _cartesian_qw_operator_receipt_unsupported_detail(args)
    isempty(args) && return "no receipt input was provided"
    first_arg = args[1]
    if first_arg isa Function
        return "receipt wrappers do not accept frontend/alias function objects; call cartesian_qw_operator_construction_receipt with the underlying basis, fixed block, and optional supplement inputs"
    elseif first_arg isa OrdinaryCartesianOperators3D
        return "received an already-built OrdinaryCartesianOperators3D payload; use cartesian_qw_operator_construction_record or cartesian_qw_operator_carried_space_sidecar for post-build audits"
    else
        return "first input type $(typeof(first_arg)) is not in the receipt coverage table"
    end
end

function _cartesian_qw_operator_receipt_unsupported_error(args)
    input_types = isempty(args) ? "()" : join(string.(typeof.(args)), ", ")
    detail = _cartesian_qw_operator_receipt_unsupported_detail(args)
    supported = _cartesian_qw_operator_receipt_supported_input_summary()
    return ArgumentError(
        "unsupported Cartesian QW operator construction receipt route for input types ($input_types): $detail. Supported receipt input families: $supported. The receipt helper is a delegate/audit wrapper only and will not guess aliases or synthesize new Hamiltonian routes.",
    )
end

function cartesian_qw_operator_construction_receipt(args...; kwargs...)
    throw(_cartesian_qw_operator_receipt_unsupported_error(args))
end

end
