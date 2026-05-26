module CartesianQWOperatorCarriedSpaces

import ..GaussletBases:
    AbstractBondAlignedOrdinaryQWBasis3D,
    CartesianBasisRepresentation3D,
    MappedUniformBasis,
    OrdinaryCartesianOperators3D,
    _NestedFixedBlock3D,
    basis_representation
import ..GaussletBases.CartesianCarriedSpaces:
    CartesianCarriedSpace3D,
    carried_space_diagnostics,
    carried_space_representation,
    cartesian_carried_space

export CartesianQWOperatorCarriedSpaceSidecar,
       cartesian_qw_operator_carried_space_sidecar,
       qw_operator_carried_space,
       qw_operator_basis_representation,
       qw_operator_carried_space_diagnostics,
       qw_operator_carried_space_provenance

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

qw_operator_carried_space(sidecar::CartesianQWOperatorCarriedSpaceSidecar) =
    sidecar.carried_space
qw_operator_basis_representation(sidecar::CartesianQWOperatorCarriedSpaceSidecar) =
    sidecar.operator_representation
qw_operator_carried_space_diagnostics(sidecar::CartesianQWOperatorCarriedSpaceSidecar) =
    sidecar.diagnostics
qw_operator_carried_space_provenance(sidecar::CartesianQWOperatorCarriedSpaceSidecar) =
    sidecar.provenance

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

end
