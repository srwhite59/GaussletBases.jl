function _cartesian_hybrid_parent_basis(
    operators::OrdinaryCartesianOperators3D,
)
    return operators.basis isa _NestedFixedBlock3D ? operators.basis.parent_basis : operators.basis
end

function _cartesian_atomic_axis_bundles(
    operators::OrdinaryCartesianOperators3D,
    parent_basis::MappedUniformBasis,
)
    gausslet_bundle = _mapped_ordinary_gausslet_1d_bundle(
        parent_basis;
        exponents = operators.expansion.exponents,
        center = 0.0,
        backend = operators.gausslet_backend,
    )
    return (x = gausslet_bundle, y = gausslet_bundle, z = gausslet_bundle)
end

function _cartesian_hybrid_supplement_axis_tables(
    operators::OrdinaryCartesianOperators3D,
    factorized_cartesian_parent_basis::_CartesianNestedFactorizedBasis3D,
    supplement3d,
    axis_bundles::NamedTuple{(:x, :y, :z)},
    route_label::AbstractString,
)
    proxy_x = _qwrg_mapped_supplement_proxy_layer(axis_bundles.x.basis, axis_bundles.x)
    proxy_y = _qwrg_mapped_supplement_proxy_layer(axis_bundles.y.basis, axis_bundles.y)
    proxy_z = _qwrg_mapped_supplement_proxy_layer(axis_bundles.z.basis, axis_bundles.z)
    norbitals = length(supplement3d.orbitals)
    x_table = zeros(Float64, size(factorized_cartesian_parent_basis.x_functions, 2), norbitals)
    y_table = zeros(Float64, size(factorized_cartesian_parent_basis.y_functions, 2), norbitals)
    z_table = zeros(Float64, size(factorized_cartesian_parent_basis.z_functions, 2), norbitals)
    for (orbital_index, orbital) in pairs(supplement3d.orbitals)
        x_data = _qwrg_atomic_axis_cross_data(proxy_x, orbital, :x, operators.expansion)
        y_data = _qwrg_atomic_axis_cross_data(proxy_y, orbital, :y, operators.expansion)
        z_data = _qwrg_atomic_axis_cross_data(proxy_z, orbital, :z, operators.expansion)
        coefficients = Vector{Float64}(orbital.coefficients)
        x_table[:, orbital_index] .=
            transpose(factorized_cartesian_parent_basis.x_functions) * (x_data.overlap * coefficients)
        y_table[:, orbital_index] .=
            transpose(factorized_cartesian_parent_basis.y_functions) * (y_data.overlap * coefficients)
        z_table[:, orbital_index] .=
            transpose(factorized_cartesian_parent_basis.z_functions) * (z_data.overlap * coefficients)
    end
    return (x = x_table, y = y_table, z = z_table)
end

function _cartesian_atomic_hybrid_overlap_sidecars(
    operators::OrdinaryCartesianOperators3D,
    cartesian_parent::CartesianBasisRepresentation3D,
)
    factorized_cartesian_parent_basis = _cartesian_optional_factorized_parent_basis(cartesian_parent)
    parent_basis = _cartesian_hybrid_parent_basis(operators)
    parent_basis isa MappedUniformBasis || throw(
        ArgumentError(
            "atomic hybrid overlap sidecars currently require a MappedUniformBasis or nested fixed block built from one",
        ),
    )
    gausslet_bundle = _mapped_ordinary_gausslet_1d_bundle(
        parent_basis;
        exponents = operators.expansion.exponents,
        center = 0.0,
        backend = operators.gausslet_backend,
    )
    supplement3d = _atomic_cartesian_shell_supplement_3d(operators.gaussian_data)
    raw_blocks = _qwrg_atomic_cartesian_blocks_3d(
        gausslet_bundle,
        supplement3d,
        operators.expansion,
    )
    parent_to_fixed_coefficients =
        cartesian_parent.coefficient_matrix === nothing ?
        Matrix{Float64}(I, cartesian_parent.metadata.final_dimension, cartesian_parent.metadata.final_dimension) :
        Matrix{Float64}(cartesian_parent.coefficient_matrix)
    base_sidecars = (
        hybrid_overlap_kind =
            isnothing(factorized_cartesian_parent_basis) ?
            :dense_atomic_mixed_raw :
            :factorized_atomic_mixed_raw,
        cartesian_probe_overlap_kind = :atomic_qw_dense_parent,
        cartesian_probe_parent_basis = parent_basis,
        cartesian_probe_expansion = operators.expansion,
        cartesian_probe_gausslet_backend = operators.gausslet_backend,
        hybrid_overlap_audit_contract =
            :block_support_sparse_parent_to_fixed_map_materialized_densely_for_exact_parent_to_supplement_overlap_audit,
        hybrid_overlap_audit_coefficient_scope = :parent_to_fixed,
        hybrid_overlap_audit_cross_overlap = :parent_to_supplement,
        hybrid_overlap_audit_parent_parent_operator = false,
        exact_cartesian_supplement_overlap =
            Matrix{Float64}(transpose(parent_to_fixed_coefficients) * raw_blocks.overlap_ga),
        exact_supplement_overlap = Matrix{Float64}(raw_blocks.overlap_aa),
    )
    isnothing(factorized_cartesian_parent_basis) && return base_sidecars
    return merge(
        base_sidecars,
        (
            factorized_cartesian_parent_basis = factorized_cartesian_parent_basis,
            cartesian_supplement_axis_tables = _cartesian_hybrid_supplement_axis_tables(
                operators,
                factorized_cartesian_parent_basis,
                supplement3d,
                _cartesian_atomic_axis_bundles(operators, parent_basis),
                "atomic",
            ),
        ),
    )
end

function _cartesian_diatomic_hybrid_overlap_sidecars(
    operators::OrdinaryCartesianOperators3D,
    cartesian_parent::CartesianBasisRepresentation3D,
)
    factorized_cartesian_parent_basis = _cartesian_optional_factorized_parent_basis(cartesian_parent)
    parent_basis = _cartesian_hybrid_parent_basis(operators)
    parent_basis isa BondAlignedDiatomicQWBasis3D || throw(
        ArgumentError(
            "bond-aligned diatomic hybrid overlap sidecars currently require a BondAlignedDiatomicQWBasis3D or nested fixed block built from one",
        ),
    )
    bundles = _qwrg_bond_aligned_axis_bundles(
        parent_basis,
        operators.expansion;
        gausslet_backend = operators.gausslet_backend,
    )
    supplement3d = _bond_aligned_diatomic_cartesian_shell_supplement_3d(operators.gaussian_data)
    overlap_blocks = _qwrg_diatomic_cartesian_shell_overlap_blocks_3d(
        bundles,
        supplement3d,
        parent_basis,
        operators.expansion,
    )
    parent_to_fixed_coefficients =
        cartesian_parent.coefficient_matrix === nothing ?
        Matrix{Float64}(I, cartesian_parent.metadata.final_dimension, cartesian_parent.metadata.final_dimension) :
        Matrix{Float64}(cartesian_parent.coefficient_matrix)
    base_sidecars = (
        hybrid_overlap_kind =
            isnothing(factorized_cartesian_parent_basis) ?
            :dense_bond_aligned_diatomic_mixed_raw :
            :factorized_bond_aligned_diatomic_mixed_raw,
        hybrid_overlap_audit_contract =
            :block_support_sparse_parent_to_fixed_map_materialized_densely_for_exact_parent_to_supplement_overlap_audit,
        hybrid_overlap_audit_coefficient_scope = :parent_to_fixed,
        hybrid_overlap_audit_cross_overlap = :parent_to_supplement,
        hybrid_overlap_audit_parent_parent_operator = false,
        exact_cartesian_supplement_overlap =
            Matrix{Float64}(transpose(parent_to_fixed_coefficients) * overlap_blocks.overlap_ga),
        exact_supplement_overlap = Matrix{Float64}(overlap_blocks.overlap_aa),
    )
    isnothing(factorized_cartesian_parent_basis) && return base_sidecars
    return merge(
        base_sidecars,
        (
            factorized_cartesian_parent_basis = factorized_cartesian_parent_basis,
            cartesian_supplement_axis_tables = _cartesian_hybrid_supplement_axis_tables(
                operators,
                factorized_cartesian_parent_basis,
                supplement3d,
                (x = bundles.bundle_x, y = bundles.bundle_y, z = bundles.bundle_z),
                "bond-aligned diatomic",
            ),
        ),
    )
end

function _cartesian_hybrid_overlap_sidecars(
    operators::OrdinaryCartesianOperators3D,
    cartesian_parent::CartesianBasisRepresentation3D,
)
    if operators.gaussian_data isa LegacyAtomicGaussianSupplement
        return _cartesian_atomic_hybrid_overlap_sidecars(operators, cartesian_parent)
    end
    parent_basis = _cartesian_hybrid_parent_basis(operators)
    if parent_basis isa BondAlignedDiatomicQWBasis3D &&
       operators.gaussian_data isa Union{
           LegacyBondAlignedDiatomicGaussianSupplement,
           LegacyBondAlignedHeteronuclearGaussianSupplement,
       }
        return _cartesian_diatomic_hybrid_overlap_sidecars(operators, cartesian_parent)
    end
    throw(
        ArgumentError(
            "exact hybrid overlap sidecars currently support atomic and bond-aligned diatomic routes; got parent basis $(typeof(parent_basis)) with supplement $(typeof(operators.gaussian_data))",
        ),
    )
end

function _qwrg_cartesian_parent_representation(
    operators::OrdinaryCartesianOperators3D,
)
    basis = operators.basis
    basis isa _NestedFixedBlock3D && return basis_representation(basis)
    basis isa AbstractBondAlignedOrdinaryQWBasis3D && return basis_representation(basis)
    basis isa MappedUniformBasis && return _cartesian_direct_product_representation(basis)
    throw(
        ArgumentError(
            "Cartesian basis representation does not yet support QW parent basis $(typeof(basis))",
        ),
    )
end

function basis_representation(operators::OrdinaryCartesianOperators3D)
    cartesian_parent = _qwrg_cartesian_parent_representation(operators)
    axis_representations = cartesian_parent.axis_representations
    axis_metadata = _cartesian_axis_metadata(axis_representations)
    size(cartesian_parent.metadata.basis_centers, 1) == operators.gausslet_count || throw(
        ArgumentError(
            "QW Cartesian basis representation requires the parent Cartesian representation to match gausslet_count",
        ),
    )
    supplement_representation = _cartesian_supplement_representation(operators.gaussian_data)
    parent_labels = vcat(
        cartesian_parent.metadata.basis_labels,
        [orbital.label for orbital in supplement_representation.orbitals],
    )
    parent_centers = vcat(
        cartesian_parent.metadata.basis_centers,
        _cartesian_supplement_center_matrix(supplement_representation),
    )
    size(operators.raw_to_final, 1) == length(parent_labels) || throw(
        ArgumentError(
            "QW Cartesian basis representation requires raw_to_final rows to match the explicit raw parent basis dimension",
        ),
    )
    final_labels = _cartesian_labels(operators.orbital_data)
    final_centers = _cartesian_center_matrix(operators.orbital_data)
    metadata = CartesianBasisMetadata3D(
        :hybrid_residual,
        _cartesian_axis_sharing(axis_representations),
        axis_metadata,
        :cartesian_plus_supplement_raw,
        cartesian_parent.metadata.parent_axis_counts,
        size(operators.raw_to_final, 1),
        size(operators.raw_to_final, 2),
        cartesian_parent.metadata.working_box,
        final_labels,
        final_centers,
        (
            gausslet_count = operators.gausslet_count,
            residual_count = operators.residual_count,
            supplement_kind = supplement_representation.supplement_kind,
            supplement_lmax = getfield(supplement_representation.metadata, :lmax),
        ),
    )
    overlap_sidecars = _cartesian_hybrid_overlap_sidecars(operators, cartesian_parent)
    return CartesianBasisRepresentation3D(
        metadata,
        axis_representations,
        :dense,
        Matrix{Float64}(operators.raw_to_final),
        parent_labels,
        parent_centers,
        nothing,
        nothing,
        (;
            cartesian_parent_representation = cartesian_parent,
            supplement_representation = supplement_representation,
            overlap_sidecars...,
        ),
    )
end

basis_metadata(operators::OrdinaryCartesianOperators3D) = basis_representation(operators).metadata
