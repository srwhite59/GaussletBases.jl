"""
    HybridMappedOrdinaryBasis1D

Hybrid ordinary one-dimensional basis built from a mapped ordinary backbone
plus a small explicit set of added core Gaussians.

This is the current practical ordinary-branch object for the friendlier
hybrid/core-supported regime. It keeps the mapped gausslet backbone visible,
but augments it with centered Gaussian support near the origin.
"""
struct HybridMappedOrdinaryBasis1D
    source_basis::MappedUniformBasis
    backend::Symbol
    core_gaussians::Vector{Gaussian}
    primitive_layer::PrimitiveSet1D
    coefficient_matrix::Matrix{Float64}
    center_data::Vector{Float64}
    integral_weight_data::Vector{Float64}
end

function Base.show(io::IO, basis::HybridMappedOrdinaryBasis1D)
    print(
        io,
        "HybridMappedOrdinaryBasis1D(backend=:",
        basis.backend,
        ", nbasis=",
        length(basis),
        ", ncore=",
        length(basis.core_gaussians),
    )
    if basis.backend != :numerical_reference
        print(io, ", experimental=true")
    end
    print(io, ")")
end

Base.length(basis::HybridMappedOrdinaryBasis1D) = size(basis.coefficient_matrix, 2)
primitive_set(basis::HybridMappedOrdinaryBasis1D) = basis.primitive_layer
primitives(basis::HybridMappedOrdinaryBasis1D) = primitives(basis.primitive_layer)
stencil_matrix(basis::HybridMappedOrdinaryBasis1D) = basis.coefficient_matrix
basis_spec(basis::HybridMappedOrdinaryBasis1D) = basis_spec(basis.source_basis)
family(basis::HybridMappedOrdinaryBasis1D) = family(basis.source_basis)
mapping(basis::HybridMappedOrdinaryBasis1D) = mapping(basis.source_basis)
centers(basis::HybridMappedOrdinaryBasis1D) = basis.center_data
reference_centers(basis::HybridMappedOrdinaryBasis1D) = basis.center_data
integral_weights(basis::HybridMappedOrdinaryBasis1D) = basis.integral_weight_data

function overlap_matrix(basis::HybridMappedOrdinaryBasis1D)
    return contract_primitive_matrix(basis, overlap_matrix(primitive_set(basis)))
end

function position_matrix(basis::HybridMappedOrdinaryBasis1D)
    return contract_primitive_matrix(basis, position_matrix(primitive_set(basis)))
end

function kinetic_matrix(basis::HybridMappedOrdinaryBasis1D)
    return contract_primitive_matrix(basis, kinetic_matrix(primitive_set(basis)))
end

function gaussian_factor_matrix(
    basis::HybridMappedOrdinaryBasis1D;
    exponent::Real,
    center::Real = 0.0,
)
    primitive_matrix = gaussian_factor_matrix(
        primitive_set(basis);
        exponent = exponent,
        center = center,
    )
    return contract_primitive_matrix(basis, primitive_matrix)
end

function gaussian_factor_matrices(
    basis::HybridMappedOrdinaryBasis1D;
    exponents::AbstractVector{<:Real},
    center::Real = 0.0,
)
    primitive_matrices = gaussian_factor_matrices(
        primitive_set(basis);
        exponents = exponents,
        center = center,
    )
    return [contract_primitive_matrix(basis, matrix) for matrix in primitive_matrices]
end

"""
    hybrid_mapped_ordinary_basis(
        basis::MappedUniformBasis;
        core_gaussians,
        backend = :pgdg_localized_experimental,
    )

Build a hybrid one-dimensional ordinary basis from a mapped ordinary backbone
plus a small explicit set of core Gaussians.

The backbone comes from the chosen mapped ordinary backend, while the added
`core_gaussians` are kept as explicit primitives before the final overlap
cleanup and COMX-style localization step.

Typical usage is:

```julia
hybrid = hybrid_mapped_ordinary_basis(
    basis;
    core_gaussians = [Gaussian(center = 0.0, width = 0.2)],
)
```
"""
function hybrid_mapped_ordinary_basis(
    basis::MappedUniformBasis;
    core_gaussians::AbstractVector{<:Gaussian},
    backend::Symbol = :pgdg_localized_experimental,
)
    layer = _mapped_ordinary_backend_layer(basis, backend)
    backbone_primitives = primitives(primitive_set(layer))
    core_primitive_list = Gaussian[gaussian for gaussian in core_gaussians]
    primitive_layer = PrimitiveSet1D(
        AbstractPrimitiveFunction1D[vcat(backbone_primitives, core_primitive_list)...];
        name = :hybrid_mapped_ordinary_primitives,
    )

    backbone_coefficients = Matrix{Float64}(stencil_matrix(layer))
    nbackbone_primitives = length(backbone_primitives)
    ncore = length(core_primitive_list)
    nbackbone_basis = size(backbone_coefficients, 2)
    coefficient_matrix = zeros(Float64, nbackbone_primitives + ncore, nbackbone_basis + ncore)
    coefficient_matrix[1:nbackbone_primitives, 1:nbackbone_basis] .= backbone_coefficients
    for index in 1:ncore
        coefficient_matrix[nbackbone_primitives + index, nbackbone_basis + index] = 1.0
    end

    overlap_seed = transpose(coefficient_matrix) * overlap_matrix(primitive_layer) * coefficient_matrix
    position_seed = transpose(coefficient_matrix) * position_matrix(primitive_layer) * coefficient_matrix
    primitive_weights = Float64[integral_weight(primitive) for primitive in primitives(primitive_layer)]
    sign_vector = vec(transpose(coefficient_matrix) * primitive_weights)
    transform, center_values = _cleanup_comx_transform(overlap_seed, position_seed, sign_vector)
    final_coefficients = coefficient_matrix * transform
    integral_weight_data = vec(transpose(final_coefficients) * primitive_weights)

    return HybridMappedOrdinaryBasis1D(
        basis,
        backend,
        core_primitive_list,
        primitive_layer,
        final_coefficients,
        center_values,
        integral_weight_data,
    )
end
