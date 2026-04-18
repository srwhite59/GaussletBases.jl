function _qwrg_bond_aligned_uses_midpoint_slab(
    basis::BondAlignedDiatomicQWBasis3D;
    atol::Float64 = 1.0e-10,
    rtol::Float64 = 1.0e-8,
)
    order, coordinates = _qwrg_bond_axis_order(basis)
    length(order) == 2 || return false
    midpoint = 0.5 * (coordinates[order[1]] + coordinates[order[2]])
    distances = abs.([coordinates[index] - midpoint for index in order])
    spacings = _qwrg_bond_axis_local_spacings(basis)[order]
    return isapprox(distances[1], distances[2]; atol = atol, rtol = rtol) &&
        isapprox(spacings[1], spacings[2]; atol = atol, rtol = rtol)
end

function _qwrg_bond_aligned_preferred_split_side(
    basis::BondAlignedDiatomicQWBasis3D;
    atol::Float64 = 1.0e-10,
    rtol::Float64 = 1.0e-8,
)
    order, _coordinates = _qwrg_bond_axis_order(basis)
    length(order) == 2 || return :left
    left_index, right_index = order
    spacings = _qwrg_bond_axis_local_spacings(basis)
    left_spacing = spacings[left_index]
    right_spacing = spacings[right_index]
    if left_spacing < right_spacing && !isapprox(left_spacing, right_spacing; atol = atol, rtol = rtol)
        return :left
    elseif right_spacing < left_spacing && !isapprox(left_spacing, right_spacing; atol = atol, rtol = rtol)
        return :right
    else
        return :left
    end
end

function _qwrg_bond_aligned_axis_bundles(
    basis::AbstractBondAlignedOrdinaryQWBasis3D,
    expansion::CoulombGaussianExpansion;
    gausslet_backend::Symbol = :numerical_reference,
)
    bundle_x = _mapped_ordinary_gausslet_1d_bundle(
        basis.basis_x;
        exponents = expansion.exponents,
        center = 0.0,
        backend = gausslet_backend,
    )
    bundle_y = _mapped_ordinary_gausslet_1d_bundle(
        basis.basis_y;
        exponents = expansion.exponents,
        center = 0.0,
        backend = gausslet_backend,
    )
    bundle_z = _mapped_ordinary_gausslet_1d_bundle(
        basis.basis_z;
        exponents = expansion.exponents,
        center = 0.0,
        backend = gausslet_backend,
    )
    return _CartesianNestedAxisBundles3D(bundle_x, bundle_y, bundle_z)
end

function _require_reference_only_gausslet_backend(
    route_label::AbstractString,
    gausslet_backend::Symbol,
)
    gausslet_backend == :numerical_reference || throw(
        ArgumentError(
            "$(route_label) is currently a numerical-reference-only route; PGDG production-contract support is not yet implemented here (got gausslet_backend = :$(gausslet_backend))",
        ),
    )
    return gausslet_backend
end

function _resolved_nested_term_coefficients(
    expansion::CoulombGaussianExpansion,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}},
)
    if isnothing(term_coefficients)
        return expansion.coefficients
    elseif term_coefficients isa Vector{Float64}
        return term_coefficients
    else
        return Float64[Float64(value) for value in term_coefficients]
    end
end

function bond_aligned_diatomic_nested_fixed_source(
    basis::BondAlignedDiatomicQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_unsplit_parallel_to_transverse_ratio_for_split::Float64 = 3.0,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    reference_fudge_factor::Float64 = 1.2,
    core_near_nucleus_protect_rows::Union{Symbol,Integer} = :auto,
    shared_shell_retain_xy::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_xz::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_yz::Union{Nothing,Tuple{Int,Int}} = nothing,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    return @timeg "diatomic.fixed_source.total" begin
        _require_reference_only_gausslet_backend(
            "bond-aligned diatomic nested fixed source",
            gausslet_backend,
        )
        resolved_term_coefficients =
            _resolved_nested_term_coefficients(expansion, term_coefficients)
        midpoint =
            sum(_qwrg_axis_coordinate(nucleus, basis.bond_axis) for nucleus in basis.nuclei) /
            length(basis.nuclei)
        bundles = @timeg "diatomic.fixed_source.axis_bundles" begin
            _qwrg_bond_aligned_axis_bundles(
                basis,
                expansion;
                gausslet_backend = gausslet_backend,
            )
        end
        use_midpoint_slab = _qwrg_bond_aligned_uses_midpoint_slab(basis)
        preferred_split_side = _qwrg_bond_aligned_preferred_split_side(basis)
        @timeg "diatomic.fixed_source.source_assembly" begin
            _nested_bond_aligned_diatomic_source(
                basis,
                bundles;
                bond_axis = basis.bond_axis,
                midpoint = midpoint,
                nside = nside,
                min_unsplit_parallel_to_transverse_ratio_for_split =
                    min_unsplit_parallel_to_transverse_ratio_for_split,
                min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
                reference_fudge_factor = reference_fudge_factor,
                core_near_nucleus_protect_rows = core_near_nucleus_protect_rows,
                use_midpoint_slab = use_midpoint_slab,
                prefer_midpoint_tie_side = preferred_split_side,
                shared_shell_retain_xy = shared_shell_retain_xy,
                shared_shell_retain_xz = shared_shell_retain_xz,
                shared_shell_retain_yz = shared_shell_retain_yz,
                term_coefficients = resolved_term_coefficients,
            )
        end
    end
end

function bond_aligned_diatomic_nested_fixed_block(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    return (
        source = source,
        fixed_block = _nested_fixed_block(source),
    )
end

function bond_aligned_diatomic_nested_fixed_block(
    basis::BondAlignedDiatomicQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_unsplit_parallel_to_transverse_ratio_for_split::Float64 = 3.0,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    reference_fudge_factor::Float64 = 1.2,
    core_near_nucleus_protect_rows::Union{Symbol,Integer} = :auto,
    shared_shell_retain_xy::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_xz::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_yz::Union{Nothing,Tuple{Int,Int}} = nothing,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    source = bond_aligned_diatomic_nested_fixed_source(
        basis;
        expansion = expansion,
        gausslet_backend = gausslet_backend,
        nside = nside,
        min_unsplit_parallel_to_transverse_ratio_for_split =
            min_unsplit_parallel_to_transverse_ratio_for_split,
        min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        reference_fudge_factor = reference_fudge_factor,
        core_near_nucleus_protect_rows = core_near_nucleus_protect_rows,
        shared_shell_retain_xy = shared_shell_retain_xy,
        shared_shell_retain_xz = shared_shell_retain_xz,
        shared_shell_retain_yz = shared_shell_retain_yz,
        term_coefficients = term_coefficients,
    )
    return bond_aligned_diatomic_nested_fixed_block(source)
end

function _bond_aligned_diatomic_nested_geometry_diagnostics(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    return @timeg "diatomic.geometry_diagnostics" begin
        contract_audit = _nested_source_contract_audit(source)
        shared_shell_dimensions = Int[size(shell.coefficient_matrix, 2) for shell in source.shared_shell_layers]
        child_sequence_dimensions = Int[
            size(sequence.coefficient_matrix, 2) for sequence in source.child_sequences
        ]
        (
            source = source,
            geometry = source.geometry,
            nside = source.nside,
            child_shell_retention_contract = source.child_shell_retention_contract,
            shared_shell_retention_contract = source.shared_shell_retention_contract,
            shared_shell_count = length(source.shared_shell_layers),
            shared_shell_dimensions = shared_shell_dimensions,
            shared_shells_match_contract =
                all(==(source.shared_shell_retention_contract.shell_increment), shared_shell_dimensions),
            child_sequence_count = length(source.child_sequences),
            child_sequence_dimensions = child_sequence_dimensions,
            fixed_dimension = size(source.sequence.coefficient_matrix, 2),
            contract_audit = contract_audit,
        )
    end
end

function bond_aligned_diatomic_nested_geometry_diagnostics(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    return _bond_aligned_diatomic_nested_geometry_diagnostics(source)
end

function bond_aligned_diatomic_nested_geometry_diagnostics(
    basis::BondAlignedDiatomicQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_unsplit_parallel_to_transverse_ratio_for_split::Float64 = 3.0,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    reference_fudge_factor::Float64 = 1.2,
    core_near_nucleus_protect_rows::Union{Symbol,Integer} = :auto,
    shared_shell_retain_xy::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_xz::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_yz::Union{Nothing,Tuple{Int,Int}} = nothing,
)
    source = bond_aligned_diatomic_nested_fixed_source(
        basis;
        expansion = expansion,
        gausslet_backend = gausslet_backend,
        nside = nside,
        min_unsplit_parallel_to_transverse_ratio_for_split =
            min_unsplit_parallel_to_transverse_ratio_for_split,
        min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        reference_fudge_factor = reference_fudge_factor,
        core_near_nucleus_protect_rows = core_near_nucleus_protect_rows,
        shared_shell_retain_xy = shared_shell_retain_xy,
        shared_shell_retain_xz = shared_shell_retain_xz,
        shared_shell_retain_yz = shared_shell_retain_yz,
    )
    return _bond_aligned_diatomic_nested_geometry_diagnostics(source)
end
