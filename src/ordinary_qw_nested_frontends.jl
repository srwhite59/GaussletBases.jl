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

struct _QWRGNestedSourceFrontendContext{B,O,C}
    basis::B
    expansion::CoulombGaussianExpansion
    gausslet_backend::Symbol
    build_options::O
    capabilities::C
end

_qwrg_optional_timeg(label::Nothing, builder::F) where {F<:Function} = builder()

function _qwrg_optional_timeg(label::AbstractString, builder::F) where {F<:Function}
    return @timeg label begin
        builder()
    end
end

_qwrg_optional_timeg(builder::F, label::Nothing) where {F<:Function} = builder()

function _qwrg_optional_timeg(builder::F, label::AbstractString) where {F<:Function}
    return @timeg label begin
        builder()
    end
end

function _normalized_nested_source_frontend_context(
    basis::BondAlignedDiatomicQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_unsplit_parallel_to_transverse_ratio_for_split::Float64 = 3.0,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    reference_fudge_factor::Float64 = 1.2,
    shared_shell_angular_resolution_scale::Float64 = 1.4,
    core_near_nucleus_protect_rows::Union{Symbol,Integer} = :auto,
    shared_shell_retain_xy::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_xz::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_yz::Union{Nothing,Tuple{Int,Int}} = nothing,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    midpoint =
        sum(_qwrg_axis_coordinate(nucleus, basis.bond_axis) for nucleus in basis.nuclei) /
        length(basis.nuclei)
    return _QWRGNestedSourceFrontendContext(
        basis,
        expansion,
        gausslet_backend,
        (
            nside = nside,
            midpoint = midpoint,
            min_unsplit_parallel_to_transverse_ratio_for_split =
                min_unsplit_parallel_to_transverse_ratio_for_split,
            min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
            reference_fudge_factor = reference_fudge_factor,
            shared_shell_angular_resolution_scale = shared_shell_angular_resolution_scale,
            core_near_nucleus_protect_rows = core_near_nucleus_protect_rows,
            use_midpoint_slab = _qwrg_bond_aligned_uses_midpoint_slab(basis),
            preferred_split_side = _qwrg_bond_aligned_preferred_split_side(basis),
            shared_shell_retain_xy = shared_shell_retain_xy,
            shared_shell_retain_xz = shared_shell_retain_xz,
            shared_shell_retain_yz = shared_shell_retain_yz,
            term_coefficients = term_coefficients,
        ),
        (
            route_label = "bond-aligned diatomic nested fixed source",
            total_timing_label = "diatomic.fixed_source.total",
            axis_bundles_timing_label = "diatomic.fixed_source.axis_bundles",
            source_assembly_timing_label = "diatomic.fixed_source.source_assembly",
        ),
    )
end

function _nested_source_from_frontend_context(
    context::_QWRGNestedSourceFrontendContext{<:BondAlignedDiatomicQWBasis3D},
    bundles::_CartesianNestedAxisBundles3D,
    term_coefficients::AbstractVector{Float64},
)
    basis = context.basis
    options = context.build_options
    return _nested_bond_aligned_diatomic_source(
        basis,
        bundles;
        bond_axis = basis.bond_axis,
        midpoint = options.midpoint,
        nside = options.nside,
        min_unsplit_parallel_to_transverse_ratio_for_split =
            options.min_unsplit_parallel_to_transverse_ratio_for_split,
        min_parallel_to_transverse_ratio = options.min_parallel_to_transverse_ratio,
        reference_fudge_factor = options.reference_fudge_factor,
        shared_shell_angular_resolution_scale = options.shared_shell_angular_resolution_scale,
        core_near_nucleus_protect_rows = options.core_near_nucleus_protect_rows,
        use_midpoint_slab = options.use_midpoint_slab,
        prefer_midpoint_tie_side = options.preferred_split_side,
        shared_shell_retain_xy = options.shared_shell_retain_xy,
        shared_shell_retain_xz = options.shared_shell_retain_xz,
        shared_shell_retain_yz = options.shared_shell_retain_yz,
        term_coefficients = term_coefficients,
    )
end

function _nested_source_frontend_source(context::_QWRGNestedSourceFrontendContext)
    return _qwrg_optional_timeg(context.capabilities.total_timing_label) do
        _require_reference_only_gausslet_backend(
            context.capabilities.route_label,
            context.gausslet_backend,
        )
        bundles = _qwrg_optional_timeg(context.capabilities.axis_bundles_timing_label) do
            _qwrg_bond_aligned_axis_bundles(
                context.basis,
                context.expansion;
                gausslet_backend = context.gausslet_backend,
            )
        end
        resolved_term_coefficients = _resolved_nested_term_coefficients(
            context.expansion,
            context.build_options.term_coefficients,
        )
        _qwrg_optional_timeg(context.capabilities.source_assembly_timing_label) do
            _nested_source_from_frontend_context(
                context,
                bundles,
                resolved_term_coefficients,
            )
        end
    end
end

function _nested_source_fixed_block(
    source,
)
    return (
        source = source,
        fixed_block = _nested_fixed_block(source),
    )
end

function _nested_source_frontend_fixed_block(
    context::_QWRGNestedSourceFrontendContext,
)
    return _nested_source_fixed_block(_nested_source_frontend_source(context))
end

struct _CartesianNestedSourceGlassBoxContract{A}
    fixed_dimension::Int
    contract_audit::A
    shared_shell_dimensions::Vector{Int}
    shared_shell_provenance::Vector{_CartesianNestedShellLayerProvenance3D}
    leaf_count::Int
end

struct _CartesianNestedGlassBoxContract{A}
    fixed_dimension::Int
    contract_audit::A
    layer_dimensions::Vector{Int}
    layer_provenance::Vector{_CartesianNestedShellLayerProvenance3D}
    leaf_count::Union{Nothing,Int}
end

function _nested_glass_box_contract_audit end

function _nested_glass_box_layer_dimensions end

function _nested_glass_box_layer_provenance end

_nested_glass_box_fixed_dimension(subject) = size(subject.sequence.coefficient_matrix, 2)

function _nested_glass_box_leaf_count(source)
    return _nested_source_leaf_count(source)
end

function _nested_glass_box_contract(subject)
    return _CartesianNestedGlassBoxContract(
        _nested_glass_box_fixed_dimension(subject),
        _nested_glass_box_contract_audit(subject),
        _nested_glass_box_layer_dimensions(subject),
        _nested_glass_box_layer_provenance(subject),
        _nested_glass_box_leaf_count(subject),
    )
end

function _nested_source_leaf_count end

function _nested_source_shared_shell_dimensions end

function _nested_source_shared_shell_provenance end

_nested_source_fixed_dimension(source) = _nested_glass_box_fixed_dimension(source)

_nested_glass_box_contract_audit(source) = _nested_source_contract_audit(source)

_nested_glass_box_layer_dimensions(source) = _nested_source_shared_shell_dimensions(source)

_nested_glass_box_layer_provenance(source) = _nested_source_shared_shell_provenance(source)

function _nested_source_common_contract(source)
    common_contract = _nested_glass_box_contract(source)
    isnothing(common_contract.leaf_count) && throw(
        ArgumentError(
            "split nested source contract requires a meaningful leaf_count on the source-backed route",
        ),
    )
    return _CartesianNestedSourceGlassBoxContract(
        common_contract.fixed_dimension,
        common_contract.contract_audit,
        common_contract.layer_dimensions,
        common_contract.layer_provenance,
        something(common_contract.leaf_count),
    )
end

function _nested_source_leaf_count(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    return length(source.child_sequences)
end

function _nested_source_shared_shell_dimensions(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    return Int[size(shell.coefficient_matrix, 2) for shell in source.shared_shell_layers]
end

function _nested_source_shared_shell_provenance(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    return _CartesianNestedShellLayerProvenance3D[
        shell.provenance for shell in source.shared_shell_layers
    ]
end

function _nested_glass_box_contract_audit(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    return _nested_source_contract_audit(source)
end

function _nested_glass_box_layer_dimensions(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    return _nested_source_shared_shell_dimensions(source)
end

function _nested_glass_box_layer_provenance(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    return _nested_source_shared_shell_provenance(source)
end

function _nested_glass_box_leaf_count(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    return _nested_source_leaf_count(source)
end

function _nested_glass_box_fixed_dimension(
    diagnostics::OneCenterAtomicNestedStructureDiagnostics,
)
    return diagnostics.total_actual_gausslet_count
end

function _nested_glass_box_contract_audit(
    diagnostics::OneCenterAtomicNestedStructureDiagnostics,
)
    return diagnostics
end

function _nested_glass_box_layer_dimensions(
    diagnostics::OneCenterAtomicNestedStructureDiagnostics,
)
    return Int[layer.retained_dimension for layer in diagnostics.layer_structures]
end

function _nested_glass_box_layer_provenance(
    diagnostics::OneCenterAtomicNestedStructureDiagnostics,
)
    return _CartesianNestedShellLayerProvenance3D[
        layer.provenance for layer in diagnostics.layer_structures
    ]
end

function _nested_glass_box_leaf_count(
    ::OneCenterAtomicNestedStructureDiagnostics,
)
    return nothing
end

function bond_aligned_diatomic_nested_fixed_source(
    basis::BondAlignedDiatomicQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_unsplit_parallel_to_transverse_ratio_for_split::Float64 = 3.0,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    reference_fudge_factor::Float64 = 1.2,
    shared_shell_angular_resolution_scale::Float64 = 1.4,
    core_near_nucleus_protect_rows::Union{Symbol,Integer} = :auto,
    shared_shell_retain_xy::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_xz::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_yz::Union{Nothing,Tuple{Int,Int}} = nothing,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    context = _normalized_nested_source_frontend_context(
        basis;
        expansion = expansion,
        gausslet_backend = gausslet_backend,
        nside = nside,
        min_unsplit_parallel_to_transverse_ratio_for_split =
            min_unsplit_parallel_to_transverse_ratio_for_split,
        min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        reference_fudge_factor = reference_fudge_factor,
        shared_shell_angular_resolution_scale = shared_shell_angular_resolution_scale,
        core_near_nucleus_protect_rows = core_near_nucleus_protect_rows,
        shared_shell_retain_xy = shared_shell_retain_xy,
        shared_shell_retain_xz = shared_shell_retain_xz,
        shared_shell_retain_yz = shared_shell_retain_yz,
        term_coefficients = term_coefficients,
    )
    return _nested_source_frontend_source(context)
end

function bond_aligned_diatomic_nested_fixed_block(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    return _nested_source_fixed_block(source)
end

function bond_aligned_diatomic_nested_fixed_block(
    basis::BondAlignedDiatomicQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_unsplit_parallel_to_transverse_ratio_for_split::Float64 = 3.0,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    reference_fudge_factor::Float64 = 1.2,
    shared_shell_angular_resolution_scale::Float64 = 1.4,
    core_near_nucleus_protect_rows::Union{Symbol,Integer} = :auto,
    shared_shell_retain_xy::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_xz::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_yz::Union{Nothing,Tuple{Int,Int}} = nothing,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    context = _normalized_nested_source_frontend_context(
        basis;
        expansion = expansion,
        gausslet_backend = gausslet_backend,
        nside = nside,
        min_unsplit_parallel_to_transverse_ratio_for_split =
            min_unsplit_parallel_to_transverse_ratio_for_split,
        min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        reference_fudge_factor = reference_fudge_factor,
        shared_shell_angular_resolution_scale = shared_shell_angular_resolution_scale,
        core_near_nucleus_protect_rows = core_near_nucleus_protect_rows,
        shared_shell_retain_xy = shared_shell_retain_xy,
        shared_shell_retain_xz = shared_shell_retain_xz,
        shared_shell_retain_yz = shared_shell_retain_yz,
        term_coefficients = term_coefficients,
    )
    return _nested_source_frontend_fixed_block(context)
end

function _bond_aligned_diatomic_nested_geometry_diagnostics(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    return @timeg "diatomic.geometry_diagnostics" begin
        common_contract = _nested_source_common_contract(source)
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
            shared_shell_dimensions = common_contract.shared_shell_dimensions,
            shared_shell_provenance = common_contract.shared_shell_provenance,
            shared_shells_match_contract =
                all(
                    ==(source.shared_shell_retention_contract.shell_increment),
                    common_contract.shared_shell_dimensions,
                ),
            child_sequence_count = common_contract.leaf_count,
            child_sequence_dimensions = child_sequence_dimensions,
            fixed_dimension = common_contract.fixed_dimension,
            contract_audit = common_contract.contract_audit,
        )
    end
end

function _nested_source_geometry_diagnostics(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    return _bond_aligned_diatomic_nested_geometry_diagnostics(source)
end

function _nested_source_frontend_geometry_diagnostics(
    context::_QWRGNestedSourceFrontendContext,
)
    return _nested_source_geometry_diagnostics(_nested_source_frontend_source(context))
end

function bond_aligned_diatomic_nested_geometry_diagnostics(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    return _nested_source_geometry_diagnostics(source)
end

function bond_aligned_diatomic_nested_geometry_diagnostics(
    basis::BondAlignedDiatomicQWBasis3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    nside::Int = 5,
    min_unsplit_parallel_to_transverse_ratio_for_split::Float64 = 3.0,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    reference_fudge_factor::Float64 = 1.2,
    shared_shell_angular_resolution_scale::Float64 = 1.4,
    core_near_nucleus_protect_rows::Union{Symbol,Integer} = :auto,
    shared_shell_retain_xy::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_xz::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_yz::Union{Nothing,Tuple{Int,Int}} = nothing,
)
    context = _normalized_nested_source_frontend_context(
        basis;
        expansion = expansion,
        gausslet_backend = gausslet_backend,
        nside = nside,
        min_unsplit_parallel_to_transverse_ratio_for_split =
            min_unsplit_parallel_to_transverse_ratio_for_split,
        min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        reference_fudge_factor = reference_fudge_factor,
        shared_shell_angular_resolution_scale = shared_shell_angular_resolution_scale,
        core_near_nucleus_protect_rows = core_near_nucleus_protect_rows,
        shared_shell_retain_xy = shared_shell_retain_xy,
        shared_shell_retain_xz = shared_shell_retain_xz,
        shared_shell_retain_yz = shared_shell_retain_yz,
    )
    return _nested_source_frontend_geometry_diagnostics(context)
end
