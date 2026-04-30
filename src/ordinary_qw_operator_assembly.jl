function _qwrg_capture_timeg_report(
    build::F,
    label::AbstractString,
    timing::Bool;
    timing_io::IO = stdout,
) where {F<:Function}
    if !timing
        return @timeg label build()
    end

    old_config = TimeG._TIMING_CONFIG[]
    state = TimeG._timing_state()
    parent = isempty(state.stack) ? nothing : state.stack[end]
    root_count = length(state.roots)
    child_count = isnothing(parent) ? 0 : length(parent.children)
    try
        set_timing!(true)
        set_timing_live!(false)
        set_timing_thresholds!(expand = 0.0, drop = 0.0)
        result = @timeg label build()
        node = isnothing(parent) ? state.roots[root_count + 1] : parent.children[child_count + 1]
        isnothing(node) && throw(
            ArgumentError("QW TimeG capture expected exactly one root timing node"),
        )
        timing_report(timing_io, TimeG.TimingReport(TimeG.TimingNode[deepcopy(node)]))
        return result
    finally
        TimeG._TIMING_CONFIG[] = old_config
    end
end

function _qwrg_print_basis_counts(
    io::IO,
    gausslet_count::Integer,
    raw_to_final::AbstractMatrix{<:Real},
)
    return _qwrg_print_basis_counts(io, "gausslet_count", gausslet_count, raw_to_final)
end

function _qwrg_print_basis_counts(
    io::IO,
    carried_label::AbstractString,
    carried_count::Integer,
    raw_to_final::AbstractMatrix{<:Real},
)
    total_basis_dim = size(raw_to_final, 2)
    residual_count = total_basis_dim - carried_count
    println(
        io,
        "QW-RG basis counts  ",
        carried_label,
        "=",
        carried_count,
        "  residual_count=",
        residual_count,
        "  total_basis_dim=",
        total_basis_dim,
    )
    flush(io)
    return nothing
end

function _validate_nuclear_term_storage(storage::Symbol)
    storage in (:auto, :total_only, :by_center) || throw(
        ArgumentError(
            "nuclear_term_storage must be :auto, :total_only, or :by_center",
        ),
    )
    return storage
end

struct _QWRGBondAlignedBuildContext{
    PB,
    C,
    CR,
    PR,
    GD,
    RM,
    PM,
    CM,
}
    basis_family::Symbol
    carried_space_kind::Symbol
    parent_basis::PB
    carried::C
    carried_representation::CR
    parent_representation::PR
    nuclei::Vector{NTuple{3,Float64}}
    default_nuclear_charges::Vector{Float64}
    contraction::Union{Nothing,_CartesianCoefficientMap}
    gaussian_data::GD
    route_metadata::RM
    parent_route_metadata::PM
    capabilities::CM
end

function _pure_bond_aligned_nested_route_label(
    basis_family::Symbol,
)
    basis_family == :bond_aligned_diatomic &&
        return "bond-aligned diatomic nested ordinary_cartesian_qiu_white_operators"
    basis_family == :bond_aligned_homonuclear_chain &&
        return "experimental bond-aligned homonuclear chain nested ordinary_cartesian_qiu_white_operators"
    basis_family == :axis_aligned_homonuclear_square_lattice &&
        return "experimental axis-aligned homonuclear square-lattice nested ordinary_cartesian_qiu_white_operators"
    throw(ArgumentError("unsupported pure bond-aligned basis_family = :$(basis_family)"))
end

function _normalized_bond_aligned_build_context(
    basis::AbstractBondAlignedOrdinaryQWBasis3D,
)
    metadata = _cartesian_bond_aligned_build_metadata(basis)
    return _QWRGBondAlignedBuildContext(
        metadata.basis_family,
        metadata.carried_metadata.basis_kind,
        basis,
        basis,
        metadata.carried_representation,
        metadata.parent_representation,
        NTuple{3,Float64}[Tuple(Float64.(nucleus)) for nucleus in basis.nuclei],
        Float64[1.0 for _ in basis.nuclei],
        nothing,
        nothing,
        metadata.carried_route_metadata,
        metadata.parent_route_metadata,
        (
            route_label = "bond-aligned ordinary_cartesian_qiu_white_operators",
            allowed_gausslet_backends = (:numerical_reference, :pgdg_localized_experimental),
            backend_support_scope = nothing,
            allowed_interaction_treatments = (:ggt_nearest, :mwg),
            localized_parent_kind = nothing,
            localized_parent_kind_backend = nothing,
            localized_parent_kind_requirement = nothing,
            interaction_treatment_error_message = nothing,
            timing_label = "qwrg.bond_aligned_ordinary.total",
        ),
    )
end

function _normalized_bond_aligned_build_context(
    fixed_block::_NestedFixedBlock3D{<:AbstractBondAlignedOrdinaryQWBasis3D},
)
    metadata = _cartesian_bond_aligned_build_metadata(fixed_block)
    parent_basis = fixed_block.parent_basis
    return _QWRGBondAlignedBuildContext(
        metadata.basis_family,
        metadata.carried_metadata.basis_kind,
        parent_basis,
        fixed_block,
        metadata.carried_representation,
        metadata.parent_representation,
        NTuple{3,Float64}[Tuple(Float64.(nucleus)) for nucleus in parent_basis.nuclei],
        Float64[1.0 for _ in parent_basis.nuclei],
        fixed_block.coefficient_matrix,
        nothing,
        metadata.carried_route_metadata,
        metadata.parent_route_metadata,
        (
            route_label = _pure_bond_aligned_nested_route_label(metadata.basis_family),
            allowed_gausslet_backends = (:numerical_reference, :pgdg_localized_experimental),
            backend_support_scope = " on pure Cartesian-parent nested fixed blocks",
            allowed_interaction_treatments = (:ggt_nearest,),
            localized_parent_kind = :cartesian_product_basis,
            localized_parent_kind_backend = :pgdg_localized_experimental,
            localized_parent_kind_requirement =
                "only for nested fixed blocks with parent_kind = {parent_kind}; supplement-bearing or transformed parent spaces remain numerical-reference-only",
            interaction_treatment_error_message = nothing,
            timing_label = "qwrg.bond_aligned_nested_fixed.total",
        ),
    )
end

function _normalized_bond_aligned_build_context(
    basis::BondAlignedDiatomicQWBasis3D,
    gaussian_data::Union{
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    },
)
    metadata = _cartesian_bond_aligned_build_metadata(basis)
    return _QWRGBondAlignedBuildContext(
        metadata.basis_family,
        metadata.carried_metadata.basis_kind,
        basis,
        basis,
        metadata.carried_representation,
        metadata.parent_representation,
        NTuple{3,Float64}[Tuple(Float64.(nucleus)) for nucleus in basis.nuclei],
        Float64[1.0 for _ in basis.nuclei],
        nothing,
        gaussian_data,
        metadata.carried_route_metadata,
        metadata.parent_route_metadata,
        (
            route_label = "bond-aligned diatomic molecular QW path",
            allowed_gausslet_backends = (:numerical_reference, :pgdg_localized_experimental),
            backend_support_scope = " on the direct-contracted molecular one-body backbone",
            allowed_interaction_treatments = (:ggt_nearest,),
            localized_parent_kind = nothing,
            localized_parent_kind_backend = nothing,
            localized_parent_kind_requirement = nothing,
            interaction_treatment_error_message = nothing,
            timing_label = "qwrg.diatomic_shell.total",
            timing_prefix = "qwrg.diatomic_shell",
            bundles_label = "qwrg.diatomic_shell.shared_bundles",
            residual_center_label = "bond-aligned diatomic residual-center extraction",
        ),
    )
end

function _normalized_bond_aligned_build_context(
    fixed_block::_NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D},
    gaussian_data::Union{
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    },
)
    metadata = _cartesian_bond_aligned_build_metadata(fixed_block)
    parent_basis = fixed_block.parent_basis
    return _QWRGBondAlignedBuildContext(
        metadata.basis_family,
        metadata.carried_metadata.basis_kind,
        parent_basis,
        fixed_block,
        metadata.carried_representation,
        metadata.parent_representation,
        NTuple{3,Float64}[Tuple(Float64.(nucleus)) for nucleus in parent_basis.nuclei],
        Float64[1.0 for _ in parent_basis.nuclei],
        fixed_block.coefficient_matrix,
        gaussian_data,
        metadata.carried_route_metadata,
        metadata.parent_route_metadata,
        (
            route_label = "bond-aligned diatomic nested molecular QW path",
            allowed_gausslet_backends = (:numerical_reference, :pgdg_localized_experimental),
            backend_support_scope = " on the direct-contracted molecular one-body backbone",
            allowed_interaction_treatments = (:ggt_nearest,),
            localized_parent_kind = :cartesian_product_basis,
            localized_parent_kind_backend = :pgdg_localized_experimental,
            localized_parent_kind_requirement =
                "only when the carried nested fixed block still has parent_kind = {parent_kind}; transformed parent spaces remain numerical-reference-only",
            interaction_treatment_error_message = nothing,
            timing_label = "qwrg.nested_diatomic_shell.total",
            timing_prefix = "qwrg.nested_diatomic_shell",
            bundles_label = "qwrg.nested_diatomic_shell.parent_bundles",
            residual_center_label = "bond-aligned diatomic nested residual-center extraction",
        ),
    )
end

_qwrg_capability_value(capabilities, name::Symbol, default = nothing) =
    hasproperty(capabilities, name) ? getproperty(capabilities, name) : default

function _qwrg_supported_symbol_list(symbols::Tuple{Vararg{Symbol}})
    return join((string(":", value) for value in symbols), " or ")
end

function _qwrg_backend_validation_error(capabilities, gausslet_backend::Symbol)
    route_label = capabilities.route_label
    allowed_gausslet_backends = _qwrg_capability_value(
        capabilities,
        :allowed_gausslet_backends,
        (:numerical_reference,),
    )
    if allowed_gausslet_backends == (:numerical_reference,)
        return "$(route_label) is currently a numerical-reference-only route; PGDG production-contract support is not yet implemented here (got gausslet_backend = :$(gausslet_backend))"
    end
    backend_support_scope = something(
        _qwrg_capability_value(capabilities, :backend_support_scope, nothing),
        "",
    )
    return "$(route_label) currently supports gausslet_backend = $(_qwrg_supported_symbol_list(allowed_gausslet_backends))$(backend_support_scope); broader PGDG production-contract support is not yet implemented here (got gausslet_backend = :$(gausslet_backend))"
end

function _validate_operator_route_backend(context, gausslet_backend::Symbol)
    allowed_gausslet_backends = _qwrg_capability_value(
        context.capabilities,
        :allowed_gausslet_backends,
        (:numerical_reference,),
    )
    gausslet_backend in allowed_gausslet_backends || throw(
        ArgumentError(
            _qwrg_backend_validation_error(context.capabilities, gausslet_backend),
        ),
    )

    localized_parent_kind = _qwrg_capability_value(
        context.capabilities,
        :localized_parent_kind,
        nothing,
    )
    localized_parent_kind_backend = _qwrg_capability_value(
        context.capabilities,
        :localized_parent_kind_backend,
        nothing,
    )
    if !isnothing(localized_parent_kind) &&
       !isnothing(localized_parent_kind_backend) &&
       gausslet_backend == localized_parent_kind_backend
        parent_kind = context.carried_representation.metadata.parent_kind
        parent_kind == localized_parent_kind || throw(
            ArgumentError(
                "$(context.capabilities.route_label) currently supports gausslet_backend = :$(localized_parent_kind_backend) " *
                replace(
                    something(
                        _qwrg_capability_value(
                            context.capabilities,
                            :localized_parent_kind_requirement,
                            nothing,
                        ),
                        "only when the carried nested fixed block still has parent_kind = {parent_kind}",
                    ),
                    "{parent_kind}" => ":" * string(localized_parent_kind),
                ) *
                " (got parent_kind = :$(parent_kind))",
            ),
        )
    end
    return gausslet_backend
end

function _validate_operator_route_interaction_treatment(
    context,
    interaction_treatment::Symbol,
)
    allowed_interaction_treatments = context.capabilities.allowed_interaction_treatments
    interaction_treatment in allowed_interaction_treatments && return interaction_treatment
    throw(
        ArgumentError(
            something(
                _qwrg_capability_value(
                    context.capabilities,
                    :interaction_treatment_error_message,
                    nothing,
                ),
                "$(context.capabilities.route_label) requires interaction_treatment = $(_qwrg_supported_symbol_list(allowed_interaction_treatments))",
            ),
        ),
    )
end

function _validate_bond_aligned_molecular_inputs(
    context::_QWRGBondAlignedBuildContext,
    nuclear_charges::AbstractVector{<:Real},
)
    length(nuclear_charges) == length(context.nuclei) || throw(
        ArgumentError(
            "$(context.capabilities.route_label) requires one nuclear charge per nucleus",
        ),
    )
    gaussian_data = something(context.gaussian_data)
    _qwrg_same_nuclei(gaussian_data.nuclei, context.nuclei) || throw(
        ArgumentError("bond-aligned diatomic molecular supplement nuclei must match the bond-aligned basis nuclei"),
    )
    return nuclear_charges
end

_resolved_nuclear_term_storage(storage::Symbol, ::BondAlignedDiatomicQWBasis3D) =
    storage == :auto ? :by_center : storage
_resolved_nuclear_term_storage(
    storage::Symbol,
    ::Union{BondAlignedHomonuclearChainQWBasis3D,AxisAlignedHomonuclearSquareLatticeQWBasis3D},
) = storage == :auto ? :total_only : storage
_resolved_nuclear_term_storage(
    storage::Symbol,
    fixed_block::_NestedFixedBlock3D{<:AbstractBondAlignedOrdinaryQWBasis3D},
) = _resolved_nuclear_term_storage(storage, fixed_block.parent_basis)

struct _QWRGAtomicBuildContext{PB,C,GD,CM}
    carried_space_kind::Symbol
    parent_basis::PB
    carried::C
    gaussian_data::GD
    contraction::Union{Nothing,_CartesianCoefficientMap}
    capabilities::CM
end

function _normalized_atomic_build_context(
    basis::MappedUniformBasis,
    gaussian_data::LegacyAtomicGaussianSupplement,
)
    return _QWRGAtomicBuildContext(
        :direct_product,
        basis,
        basis,
        gaussian_data,
        nothing,
        (
            route_label = "ordinary_cartesian_qiu_white_operators",
            allowed_gausslet_backends = (:numerical_reference,),
            allowed_interaction_treatments = (:ggt_nearest, :mwg),
            interaction_treatment_error_message =
                "Qiu-White interaction_treatment must be :ggt_nearest or :mwg",
            timing_label = "qwrg.atomic_shell.total",
            bundle_label = "qwrg.atomic_shell.shared_bundle",
            raw_blocks_label = "qwrg.atomic_shell.raw_blocks",
            residual_space_label = "qwrg.atomic_shell.residual_space",
            one_body_label = "qwrg.atomic_shell.one_body",
            centers_label = "qwrg.atomic_shell.centers",
            interaction_label = "qwrg.atomic_shell.interaction",
            carried_count_label = "gausslet_count",
            residual_center_label = "Qiu-White residual-center extraction",
        ),
    )
end

function _normalized_atomic_build_context(
    fixed_block::_NestedFixedBlock3D,
    gaussian_data::LegacyAtomicGaussianSupplement,
)
    return _QWRGAtomicBuildContext(
        :nested_fixed_block,
        fixed_block.parent_basis,
        fixed_block,
        gaussian_data,
        fixed_block.coefficient_matrix,
        (
            route_label = "nested ordinary_cartesian_qiu_white_operators",
            allowed_gausslet_backends = (:numerical_reference,),
            allowed_interaction_treatments = (:ggt_nearest,),
            interaction_treatment_error_message =
                "nested ordinary_cartesian_qiu_white_operators currently supports only interaction_treatment = :ggt_nearest",
            timing_label = "qwrg.nested_atomic_shell.total",
            bundle_label = "qwrg.nested_atomic_shell.parent_bundle",
            raw_blocks_label = "qwrg.nested_atomic_shell.raw_blocks",
            residual_space_label = "qwrg.nested_atomic_shell.residual_space",
            one_body_label = "qwrg.nested_atomic_shell.one_body",
            centers_label = "qwrg.nested_atomic_shell.centers",
            interaction_label = "qwrg.nested_atomic_shell.interaction",
            carried_count_label = "fixed_count",
            residual_center_label = "nested QW-PGDG residual-center extraction",
        ),
    )
end

function _same_nuclear_charge_configuration(
    left::Union{Nothing,AbstractVector{<:Real}},
    right::Union{Nothing,AbstractVector{<:Real}};
    tol::Real = 1.0e-12,
)
    (left === nothing || right === nothing) && return left === right
    length(left) == length(right) || return false
    for index in eachindex(left, right)
        abs(Float64(left[index]) - Float64(right[index])) <= tol || return false
    end
    return true
end

function _assemble_one_body_hamiltonian(
    kinetic_one_body::AbstractMatrix{<:Real},
    nuclear_one_body_by_center::AbstractVector{<:AbstractMatrix{<:Real}},
    nuclear_charges::AbstractVector{<:Real},
)
    length(nuclear_charges) == length(nuclear_one_body_by_center) || throw(
        ArgumentError("one-body reassembly requires one nuclear charge per stored nuclear term"),
    )
    matrix = Matrix{Float64}(kinetic_one_body)
    for (charge, nuclear_matrix) in zip(nuclear_charges, nuclear_one_body_by_center)
        matrix .+= Float64(charge) .* nuclear_matrix
    end
    return Matrix{Float64}(0.5 .* (matrix .+ transpose(matrix)))
end

"""
    assembled_one_body_hamiltonian(operators::OrdinaryCartesianOperators3D;
                                   nuclear_charges = operators.nuclear_charges)

Reassemble the final-basis one-body Hamiltonian from stored kinetic and
per-center nuclear-attraction terms when available.

This reassembly uses the kinetic matrix carried by the operator payload itself.
For nested fixed-block routes, that kinetic follows the fixed-block packet
contract and should not be assumed to equal a later contraction of some
separate parent-space ordinary one-body path.

If the operator payload does not retain per-center nuclear terms, the stored
`one_body_hamiltonian` can only be returned for the original nuclear charges;
requesting modified charges then throws an `ArgumentError`.
"""
function assembled_one_body_hamiltonian(
    operators::OrdinaryCartesianOperators3D;
    nuclear_charges = operators.nuclear_charges,
)
    if isnothing(operators.kinetic_one_body) || isnothing(operators.nuclear_one_body_by_center)
        if nuclear_charges === operators.nuclear_charges ||
           _same_nuclear_charge_configuration(nuclear_charges, operators.nuclear_charges)
            return Matrix{Float64}(operators.one_body_hamiltonian)
        end
        throw(
            ArgumentError(
                "this operator payload does not retain per-center nuclear one-body terms; rebuild with nuclear_term_storage = :by_center",
            ),
        )
    end
    charges = isnothing(nuclear_charges) ? operators.nuclear_charges : nuclear_charges
    isnothing(charges) && return Matrix{Float64}(operators.one_body_hamiltonian)
    return _assemble_one_body_hamiltonian(
        operators.kinetic_one_body,
        operators.nuclear_one_body_by_center,
        charges,
    )
end

function ordinary_cartesian_vee_expectation(
    operators::OrdinaryCartesianOperators3D,
    orbital::AbstractVector;
    overlap_tol::Real = 1.0e-8,
)
    length(orbital) == length(operators.orbital_data) ||
        throw(ArgumentError("orbital length must match the Qiu-White ordinary Cartesian orbital dimension"))
    overlap_error = norm(operators.overlap - I, Inf)
    overlap_error <= Float64(overlap_tol) || throw(
        ArgumentError(
            "ordinary_cartesian_vee_expectation currently requires an orthonormal final basis; got overlap error $(overlap_error)",
        ),
    )
    weights = Float64[abs2(coefficient) for coefficient in orbital]
    norm2 = sum(weights)
    norm2 > 0.0 || throw(ArgumentError("orbital must have nonzero norm"))
    weights ./= norm2
    return Float64(real(dot(weights, operators.interaction_matrix * weights)))
end

function ordinary_cartesian_1s2_check(
    operators::OrdinaryCartesianOperators3D;
    overlap_tol::Real = 1.0e-8,
)
    decomposition = eigen(Hermitian(operators.one_body_hamiltonian))
    orbital = decomposition.vectors[:, 1]
    return (
        orbital_energy = Float64(decomposition.values[1]),
        orbital = orbital,
        vee_expectation = ordinary_cartesian_vee_expectation(
            operators,
            orbital;
            overlap_tol = overlap_tol,
        ),
        overlap_error = norm(operators.overlap - I, Inf),
    )
end

function _qwrg_carried_plus_residual_split(
    raw_to_final::AbstractMatrix{<:Real},
    carried_count::Integer,
    added_count::Integer;
    structure_tol::Real = 1.0e-12,
)
    size(raw_to_final, 1) == carried_count + added_count || return nothing
    size(raw_to_final, 2) >= carried_count || return nothing

    carried_block = Matrix{Float64}(raw_to_final[:, 1:carried_count])
    carried_identity = Matrix{Float64}(I, carried_count, carried_count)
    isapprox(
        carried_block[1:carried_count, :],
        carried_identity;
        atol = structure_tol,
        rtol = 0.0,
    ) || return nothing
    if added_count > 0 &&
       maximum(abs.(carried_block[(carried_count + 1):end, :])) > structure_tol
        return nothing
    end

    return (
        residual_carried_coefficients = Matrix{Float64}(raw_to_final[1:carried_count, (carried_count + 1):end]),
        residual_added_coefficients = Matrix{Float64}(raw_to_final[(carried_count + 1):end, (carried_count + 1):end]),
    )
end

function _qwrg_structured_final_one_body_matrices(
    carried_one_body::AbstractMatrix{<:Real},
    one_body_ga::AbstractMatrix{<:Real},
    one_body_aa::AbstractMatrix{<:Real},
    raw_to_final::AbstractMatrix{<:Real};
    structure_tol::Real = 1.0e-12,
)
    carried_count = size(carried_one_body, 1)
    size(carried_one_body, 2) == carried_count || throw(
        ArgumentError("structured final one-body mix requires a square carried block"),
    )
    size(one_body_ga, 1) == carried_count || throw(
        ArgumentError("structured final one-body mix requires a GA block with carried-row count"),
    )

    added_count = size(one_body_ga, 2)
    size(one_body_aa, 1) == added_count || throw(
        ArgumentError("structured final one-body mix requires a square added-added block"),
    )
    size(one_body_aa, 2) == added_count || throw(
        ArgumentError("structured final one-body mix requires a square added-added block"),
    )

    split = _qwrg_carried_plus_residual_split(
        raw_to_final,
        carried_count,
        added_count;
        structure_tol = structure_tol,
    )
    isnothing(split) && return nothing

    carried = Matrix{Float64}(carried_one_body)
    residual_count = size(raw_to_final, 2) - carried_count
    residual_count >= 0 || throw(
        ArgumentError("structured final one-body mix requires final dimension >= carried dimension"),
    )
    if residual_count == 0
        return nothing, 0.5 .* (carried .+ transpose(carried))
    end

    coupling = Matrix{Float64}(one_body_ga)
    added = Matrix{Float64}(one_body_aa)
    residual_carried_coefficients = split.residual_carried_coefficients
    residual_added_coefficients = split.residual_added_coefficients

    carried_residual =
        carried * residual_carried_coefficients +
        coupling * residual_added_coefficients
    residual_carried =
        transpose(residual_carried_coefficients) * carried +
        transpose(residual_added_coefficients) * transpose(coupling)
    residual_residual =
        transpose(residual_carried_coefficients) * carried_residual +
        transpose(residual_added_coefficients) *
        (
            transpose(coupling) * residual_carried_coefficients +
            added * residual_added_coefficients
        )

    final_one_body = Matrix{Float64}(undef, size(raw_to_final, 2), size(raw_to_final, 2))
    final_one_body[1:carried_count, 1:carried_count] = carried
    final_one_body[1:carried_count, (carried_count + 1):end] = carried_residual
    final_one_body[(carried_count + 1):end, 1:carried_count] = residual_carried
    final_one_body[(carried_count + 1):end, (carried_count + 1):end] = residual_residual
    return nothing, 0.5 .* (final_one_body .+ transpose(final_one_body))
end

function _qwrg_final_one_body_matrix_from_blocks(
    carried_one_body::AbstractMatrix{<:Real},
    one_body_ga::AbstractMatrix{<:Real},
    one_body_aa::AbstractMatrix{<:Real},
    raw_to_final::AbstractMatrix{<:Real},
)
    structured = _qwrg_structured_final_one_body_matrices(
        carried_one_body,
        one_body_ga,
        one_body_aa,
        raw_to_final,
    )
    !isnothing(structured) && return structured

    raw_one_body = [
        Matrix{Float64}(carried_one_body) Matrix{Float64}(one_body_ga)
        Matrix{Float64}(transpose(one_body_ga)) Matrix{Float64}(one_body_aa)
    ]
    final_one_body = Matrix{Float64}(transpose(raw_to_final) * raw_one_body * raw_to_final)
    return raw_one_body, 0.5 .* (final_one_body .+ transpose(final_one_body))
end

function _qwrg_one_body_matrices(
    gausslet_one_body::AbstractMatrix{<:Real},
    one_body_ga::AbstractMatrix{<:Real},
    one_body_aa::AbstractMatrix{<:Real},
    raw_to_final::AbstractMatrix{<:Real},
)
    return _qwrg_final_one_body_matrix_from_blocks(
        gausslet_one_body,
        one_body_ga,
        one_body_aa,
        raw_to_final,
    )
end

function _qwrg_final_nuclear_one_body_by_center(
    carried_nuclear_one_body_by_center::AbstractVector{<:AbstractMatrix{<:Real}},
    nuclear_ga_by_center::AbstractVector{<:AbstractMatrix{<:Real}},
    nuclear_aa_by_center::AbstractVector{<:AbstractMatrix{<:Real}},
    raw_to_final::AbstractMatrix{<:Real},
)
    length(carried_nuclear_one_body_by_center) == length(nuclear_ga_by_center) || throw(
        ArgumentError("final nuclear one-body assembly requires one GA block per center"),
    )
    length(carried_nuclear_one_body_by_center) == length(nuclear_aa_by_center) || throw(
        ArgumentError("final nuclear one-body assembly requires one AA block per center"),
    )
    return [
        _qwrg_one_body_matrices(
            carried_nuclear_one_body_by_center[index],
            nuclear_ga_by_center[index],
            nuclear_aa_by_center[index],
            raw_to_final,
        )[2] for index in eachindex(carried_nuclear_one_body_by_center)
    ]
end

function _qwrg_nearest_indices(
    gausslet_orbitals::AbstractVector{<:CartesianProductOrbital3D},
    residual_centers::AbstractMatrix{<:Real},
)
    return Int[
        argmin(
            [
                (orbital.x - residual_centers[index, 1])^2 +
                (orbital.y - residual_centers[index, 2])^2 +
                (orbital.z - residual_centers[index, 3])^2 for orbital in gausslet_orbitals
            ],
        ) for index in axes(residual_centers, 1)
    ]
end

function _qwrg_nearest_indices(
    fixed_centers::AbstractMatrix{<:Real},
    residual_centers::AbstractMatrix{<:Real},
)
    size(fixed_centers, 2) == 3 || throw(
        ArgumentError("nearest-center selection for a nested fixed block requires centers as an n×3 matrix"),
    )
    return Int[
        argmin(
            [
                (fixed_centers[fixed, 1] - residual_centers[index, 1])^2 +
                (fixed_centers[fixed, 2] - residual_centers[index, 2])^2 +
                (fixed_centers[fixed, 3] - residual_centers[index, 3])^2 for fixed in axes(fixed_centers, 1)
            ],
        ) for index in axes(residual_centers, 1)
    ]
end

function _qwrg_effective_gaussians(
    residual_centers::AbstractMatrix{<:Real},
    residual_widths::AbstractMatrix{<:Real},
)
    nresidual = size(residual_centers, 1)
    x_gaussians = Gaussian[]
    y_gaussians = Gaussian[]
    z_gaussians = Gaussian[]
    for index in 1:nresidual
        push!(
            x_gaussians,
            Gaussian(center = residual_centers[index, 1], width = residual_widths[index, 1]),
        )
        push!(
            y_gaussians,
            Gaussian(center = residual_centers[index, 2], width = residual_widths[index, 2]),
        )
        push!(
            z_gaussians,
            Gaussian(center = residual_centers[index, 3], width = residual_widths[index, 3]),
        )
    end
    return x_gaussians, y_gaussians, z_gaussians
end

function _qwrg_same_gaussians(
    left::AbstractVector{<:Gaussian},
    right::AbstractVector{<:Gaussian};
    tol::Real = 1.0e-12,
)
    length(left) == length(right) || return false
    for index in eachindex(left, right)
        abs(left[index].center_value - right[index].center_value) <= tol || return false
        abs(left[index].width - right[index].width) <= tol || return false
    end
    return true
end

function _qwrg_interaction_matrix_nearest(
    gausslet_interaction::AbstractMatrix{<:Real},
    gausslet_orbitals::AbstractVector{<:CartesianProductOrbital3D},
    residual_centers::AbstractMatrix{<:Real},
)
    ngausslet = size(gausslet_interaction, 1)
    nresidual = size(residual_centers, 1)
    interaction = zeros(Float64, ngausslet + nresidual, ngausslet + nresidual)
    interaction[1:ngausslet, 1:ngausslet] .= Matrix{Float64}(gausslet_interaction)
    nearest = _qwrg_nearest_indices(gausslet_orbitals, residual_centers)
    for residual in 1:nresidual
        index = nearest[residual]
        interaction[1:ngausslet, ngausslet + residual] .= gausslet_interaction[:, index]
        interaction[ngausslet + residual, 1:ngausslet] .= gausslet_interaction[:, index]
    end
    for i in 1:nresidual, j in i:nresidual
        value = gausslet_interaction[nearest[i], nearest[j]]
        interaction[ngausslet + i, ngausslet + j] = value
        interaction[ngausslet + j, ngausslet + i] = value
    end
    return interaction
end

function _qwrg_interaction_matrix_nearest(
    fixed_interaction::AbstractMatrix{<:Real},
    fixed_centers::AbstractMatrix{<:Real},
    residual_centers::AbstractMatrix{<:Real},
)
    nfixed = size(fixed_interaction, 1)
    nresidual = size(residual_centers, 1)
    interaction = zeros(Float64, nfixed + nresidual, nfixed + nresidual)
    interaction[1:nfixed, 1:nfixed] .= Matrix{Float64}(fixed_interaction)
    nearest = _qwrg_nearest_indices(fixed_centers, residual_centers)
    for residual in 1:nresidual
        index = nearest[residual]
        interaction[1:nfixed, nfixed + residual] .= fixed_interaction[:, index]
        interaction[nfixed + residual, 1:nfixed] .= fixed_interaction[:, index]
    end
    for i in 1:nresidual, j in i:nresidual
        value = fixed_interaction[nearest[i], nearest[j]]
        interaction[nfixed + i, nfixed + j] = value
        interaction[nfixed + j, nfixed + i] = value
    end
    return interaction
end

# Alg QW-RG step 8b, 8c, and 8d: Match exact residual moments by effective
# Gaussians and keep RG terms in the same two-index IDA representation.
# See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
function _qwrg_interaction_matrix_mwg(
    gausslet_bundle::_MappedOrdinaryGausslet1DBundle,
    gausslet_interaction::AbstractMatrix{<:Real},
    expansion::CoulombGaussianExpansion,
    residual_centers::AbstractMatrix{<:Real},
    residual_widths::AbstractMatrix{<:Real},
)
    ngausslet = size(gausslet_interaction, 1)
    nresidual = size(residual_centers, 1)
    interaction = zeros(Float64, ngausslet + nresidual, ngausslet + nresidual)
    interaction[1:ngausslet, 1:ngausslet] .= Matrix{Float64}(gausslet_interaction)

    x_gaussians, y_gaussians, z_gaussians = _qwrg_effective_gaussians(residual_centers, residual_widths)
    # Alg QW-RG step 8d: Evaluate RG-gausslet terms from exact contracted
    # 1D raw-space pair-factor blocks, then assemble the 3D IDA interaction.
    # See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
    pair_x = _qwrg_split_block_matrices(
        gausslet_bundle,
        x_gaussians,
        expansion,
    )
    pair_y = _qwrg_same_gaussians(x_gaussians, y_gaussians) ? pair_x :
        _qwrg_split_block_matrices(gausslet_bundle, y_gaussians, expansion)
    pair_z = if _qwrg_same_gaussians(x_gaussians, z_gaussians)
        pair_x
    elseif _qwrg_same_gaussians(y_gaussians, z_gaussians)
        pair_y
    else
        _qwrg_split_block_matrices(gausslet_bundle, z_gaussians, expansion)
    end
    analytic_x = _qwrg_gaussian_analytic_blocks(x_gaussians, expansion)
    analytic_y = _qwrg_same_gaussians(x_gaussians, y_gaussians) ? analytic_x : _qwrg_gaussian_analytic_blocks(y_gaussians, expansion)
    analytic_z = if _qwrg_same_gaussians(x_gaussians, z_gaussians)
        analytic_x
    elseif _qwrg_same_gaussians(y_gaussians, z_gaussians)
        analytic_y
    else
        _qwrg_gaussian_analytic_blocks(z_gaussians, expansion)
    end

    for residual in 1:nresidual
        column = zeros(Float64, ngausslet)
        scratch = zeros(Float64, ngausslet)
        for term in eachindex(expansion.coefficients)
            fx = view(pair_x.pair_ga[term], :, residual)
            fy = view(pair_y.pair_ga[term], :, residual)
            fz = view(pair_z.pair_ga[term], :, residual)
            _qwrg_fill_product_column!(scratch, fx, fy, fz)
            column .+= expansion.coefficients[term] .* scratch
        end
        interaction[1:ngausslet, ngausslet + residual] .= column
        interaction[ngausslet + residual, 1:ngausslet] .= column
    end

    for i in 1:nresidual, j in i:nresidual
        value = 0.0
        for term in eachindex(expansion.coefficients)
            value += expansion.coefficients[term] *
                analytic_x.pair_aa[term][i, j] *
                analytic_y.pair_aa[term][i, j] *
                analytic_z.pair_aa[term][i, j]
        end
        interaction[ngausslet + i, ngausslet + j] = value
        interaction[ngausslet + j, ngausslet + i] = value
    end

    return interaction
end

function _qwrg_residual_label(
    residual_index::Integer,
    owner_nucleus_index::Integer,
    owner_counts::Dict{Int,Int},
)
    owner = Int(owner_nucleus_index)
    if owner == 0
        return string("rg", residual_index)
    end
    owner_counts[owner] = get(owner_counts, owner, 0) + 1
    owner_label = 1 <= owner <= 26 ? string(Char(Int('A') + owner - 1)) : string("N", owner)
    return string("rg", owner_label, owner_counts[owner])
end

function _qwrg_orbital_data(
    gausslet_orbitals::AbstractVector{<:CartesianProductOrbital3D},
    residual_centers::AbstractMatrix{<:Real},
    residual_widths::AbstractMatrix{<:Real},
    residual_nucleus_indices::AbstractVector{<:Integer} = zeros(Int, size(residual_centers, 1)),
)
    length(residual_nucleus_indices) == size(residual_centers, 1) || throw(
        DimensionMismatch("residual orbital data requires one owner nucleus index per residual center"),
    )
    orbitals_out = OrdinaryCartesianOrbital3D[]
    for orbital in gausslet_orbitals
        push!(
            orbitals_out,
            OrdinaryCartesianOrbital3D(
                orbital.index,
                :gausslet,
                "g($(orbital.ix),$(orbital.iy),$(orbital.iz))",
                orbital.x,
                orbital.y,
                orbital.z,
                NaN,
                NaN,
                NaN,
            ),
        )
    end
    base_index = length(gausslet_orbitals)
    residual_owner_counts = Dict{Int,Int}()
    for index in axes(residual_centers, 1)
        owner_nucleus_index = Int(residual_nucleus_indices[index])
        push!(
            orbitals_out,
            OrdinaryCartesianOrbital3D(
                base_index + index,
                :residual_gaussian,
                _qwrg_residual_label(index, owner_nucleus_index, residual_owner_counts),
                residual_centers[index, 1],
                residual_centers[index, 2],
                residual_centers[index, 3],
                residual_widths[index, 1],
                residual_widths[index, 2],
                residual_widths[index, 3],
                owner_nucleus_index,
            ),
        )
    end
    return orbitals_out
end

function _qwrg_orbital_data(
    fixed_centers::AbstractMatrix{<:Real},
    residual_centers::AbstractMatrix{<:Real},
    residual_widths::AbstractMatrix{<:Real};
    fixed_kind::Symbol = :nested_fixed,
    fixed_label_prefix::AbstractString = "nf",
    residual_nucleus_indices::AbstractVector{<:Integer} = zeros(Int, size(residual_centers, 1)),
)
    size(fixed_centers, 2) == 3 || throw(
        ArgumentError("nested fixed-block orbital data requires an n×3 center matrix"),
    )
    length(residual_nucleus_indices) == size(residual_centers, 1) || throw(
        DimensionMismatch("residual orbital data requires one owner nucleus index per residual center"),
    )
    orbitals_out = OrdinaryCartesianOrbital3D[]
    for index in axes(fixed_centers, 1)
        push!(
            orbitals_out,
            OrdinaryCartesianOrbital3D(
                index,
                fixed_kind,
                string(fixed_label_prefix, index),
                fixed_centers[index, 1],
                fixed_centers[index, 2],
                fixed_centers[index, 3],
                NaN,
                NaN,
                NaN,
            ),
        )
    end
    base_index = size(fixed_centers, 1)
    residual_owner_counts = Dict{Int,Int}()
    for index in axes(residual_centers, 1)
        owner_nucleus_index = Int(residual_nucleus_indices[index])
        push!(
            orbitals_out,
            OrdinaryCartesianOrbital3D(
                base_index + index,
                :residual_gaussian,
                _qwrg_residual_label(index, owner_nucleus_index, residual_owner_counts),
                residual_centers[index, 1],
                residual_centers[index, 2],
                residual_centers[index, 3],
                residual_widths[index, 1],
                residual_widths[index, 2],
                residual_widths[index, 3],
                owner_nucleus_index,
            ),
        )
    end
    return orbitals_out
end

function _qwrg_diatomic_one_body_matrix(
    basis::AbstractBondAlignedOrdinaryQWBasis3D,
    bundle_x::_MappedOrdinaryGausslet1DBundle,
    bundle_y::_MappedOrdinaryGausslet1DBundle,
    bundle_z::_MappedOrdinaryGausslet1DBundle,
    expansion::CoulombGaussianExpansion,
    nuclear_charges::AbstractVector{<:Real},
)
    length(nuclear_charges) == length(basis.nuclei) || throw(
        ArgumentError("bond-aligned ordinary QW path requires one nuclear charge per nucleus"),
    )
    all(Float64(value) >= 0.0 for value in nuclear_charges) || throw(
        ArgumentError("bond-aligned ordinary QW path requires nonnegative nuclear charges"),
    )
    kinetic_one_body = _qwrg_diatomic_kinetic_matrix(bundle_x, bundle_y, bundle_z)
    nuclear_one_body_by_center = _qwrg_diatomic_nuclear_one_body_by_center(
        basis,
        bundle_x,
        bundle_y,
        bundle_z,
        expansion,
    )
    return _assemble_one_body_hamiltonian(
        kinetic_one_body,
        nuclear_one_body_by_center,
        nuclear_charges,
    )
end

function _ordinary_cartesian_qiu_white_operators_pure_bond_aligned_direct(
    context::_QWRGBondAlignedBuildContext;
    nuclear_charges::AbstractVector{<:Real},
    nuclear_term_storage::Symbol,
    expansion::CoulombGaussianExpansion,
    interaction_treatment::Symbol,
    gausslet_backend::Symbol,
)
    basis = context.parent_basis
    resolved_nuclear_term_storage = _resolved_nuclear_term_storage(
        _validate_nuclear_term_storage(nuclear_term_storage),
        basis,
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

    overlap = _qwrg_diatomic_overlap_matrix(bundle_x, bundle_y, bundle_z)
    kinetic_one_body = _qwrg_diatomic_kinetic_matrix(bundle_x, bundle_y, bundle_z)
    direct_contracted_nuclear =
        resolved_nuclear_term_storage == :by_center ?
        _qwrg_bond_aligned_direct_contracted_nuclear_one_body_by_center(
            basis,
            _qwrg_full_cartesian_product_factorized_basis((
                size(bundle_x.pgdg_intermediate.overlap, 1),
                size(bundle_y.pgdg_intermediate.overlap, 1),
                size(bundle_z.pgdg_intermediate.overlap, 1),
            )),
            bundle_x,
            bundle_y,
            bundle_z,
            expansion,
        ) :
        nothing
    nuclear_one_body_by_center =
        resolved_nuclear_term_storage == :by_center ?
        direct_contracted_nuclear :
        nothing
    one_body_hamiltonian =
        !isnothing(direct_contracted_nuclear) ?
        _assemble_one_body_hamiltonian(
            kinetic_one_body,
            direct_contracted_nuclear,
            nuclear_charges,
        ) :
        (
            isnothing(nuclear_one_body_by_center) ?
            _qwrg_diatomic_one_body_matrix(
                basis,
                bundle_x,
                bundle_y,
                bundle_z,
                expansion,
                nuclear_charges,
            ) :
            _assemble_one_body_hamiltonian(
                kinetic_one_body,
                nuclear_one_body_by_center,
                nuclear_charges,
            )
        )
    interaction_matrix = _qwrg_diatomic_interaction_matrix(
        bundle_x,
        bundle_y,
        bundle_z,
        expansion,
    )

    gausslet_orbitals = _mapped_cartesian_orbitals(
        centers(basis.basis_x),
        centers(basis.basis_y),
        centers(basis.basis_z),
    )
    gausslet_count = length(gausslet_orbitals)
    zero_residual_centers = zeros(Float64, 0, 3)
    zero_residual_widths = zeros(Float64, 0, 3)

    return OrdinaryCartesianOperators3D(
        basis,
        nothing,
        gausslet_backend,
        interaction_treatment,
        expansion,
        overlap,
        one_body_hamiltonian,
        interaction_matrix,
        _qwrg_orbital_data(
            gausslet_orbitals,
            zero_residual_centers,
            zero_residual_widths,
        ),
        gausslet_count,
        0,
        Matrix{Float64}(I, gausslet_count, gausslet_count),
        zero_residual_centers,
        zero_residual_widths,
        Float64[Float64(value) for value in nuclear_charges],
        resolved_nuclear_term_storage == :by_center ? Matrix{Float64}(kinetic_one_body) : nothing,
        resolved_nuclear_term_storage == :by_center ? nuclear_one_body_by_center : nothing,
        resolved_nuclear_term_storage,
    )
end

"""
    ordinary_cartesian_qiu_white_operators(
        basis::BondAlignedDiatomicQWBasis3D;
        nuclear_charges = fill(1.0, length(basis.nuclei)),
        expansion = coulomb_gaussian_expansion(doacc = false),
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        timing = false,
    )

Build the first bond-aligned diatomic ordinary QW reference Hamiltonian on the
diatomic distortion path.

This first molecular pass is intentionally narrow:

- one bond-aligned homonuclear diatomic basis object
- no nested fixed block yet
- no molecular Gaussian supplement yet
- the final basis is the distorted 3D gausslet product basis itself, so the
  residual-Gaussian sector is empty
- this pure direct-product route accepts `gausslet_backend = :pgdg_localized_experimental`
"""
function ordinary_cartesian_qiu_white_operators(
    basis::BondAlignedDiatomicQWBasis3D;
    nuclear_charges::AbstractVector{<:Real} = fill(1.0, length(basis.nuclei)),
    nuclear_term_storage::Symbol = :auto,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    timing::Bool = false,
)
    context = _normalized_bond_aligned_build_context(basis)
    return _ordinary_cartesian_qiu_white_operators_pure_bond_aligned(
        context;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        timing = timing,
    )
end

"""
    ordinary_cartesian_qiu_white_operators(
        basis::BondAlignedHomonuclearChainQWBasis3D;
        nuclear_charges = fill(1.0, length(basis.nuclei)),
        expansion = coulomb_gaussian_expansion(doacc = false),
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        timing = false,
    )

Build the first experimental ordinary QW reference Hamiltonian on the
homonuclear chain distortion path.

This chain milestone is intentionally narrow:

- homonuclear linear chains only
- ordinary QW product basis only
- no nested chain splitting
- no molecular-shell supplement route
- this pure direct-product route accepts `gausslet_backend = :pgdg_localized_experimental`
"""
function ordinary_cartesian_qiu_white_operators(
    basis::BondAlignedHomonuclearChainQWBasis3D;
    nuclear_charges::AbstractVector{<:Real} = fill(1.0, length(basis.nuclei)),
    nuclear_term_storage::Symbol = :auto,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    timing::Bool = false,
)
    context = _normalized_bond_aligned_build_context(basis)
    return _ordinary_cartesian_qiu_white_operators_pure_bond_aligned(
        context;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        timing = timing,
    )
end

"""
    ordinary_cartesian_qiu_white_operators(
        basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D;
        nuclear_charges = fill(1.0, length(basis.nuclei)),
        expansion = coulomb_gaussian_expansion(doacc = false),
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        timing = false,
    )

Build the first experimental ordinary QW reference Hamiltonian on the
homonuclear square-lattice distortion path.

This lattice milestone is intentionally narrow:

- homonuclear `n × n` square lattices only
- ordinary QW product basis only
- no nested planar splitting
- no molecular-shell supplement route
- this pure direct-product route accepts `gausslet_backend = :pgdg_localized_experimental`
"""
function ordinary_cartesian_qiu_white_operators(
    basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D;
    nuclear_charges::AbstractVector{<:Real} = fill(1.0, length(basis.nuclei)),
    nuclear_term_storage::Symbol = :auto,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    timing::Bool = false,
)
    context = _normalized_bond_aligned_build_context(basis)
    return _ordinary_cartesian_qiu_white_operators_pure_bond_aligned(
        context;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        timing = timing,
    )
end

function _qwrg_bond_aligned_molecular_carried_data(
    context::_QWRGBondAlignedBuildContext,
    bundles,
)
    if context.carried_space_kind == :direct_product
        basis = context.parent_basis
        orbitals = _mapped_cartesian_orbitals(
            centers(basis.basis_x),
            centers(basis.basis_y),
            centers(basis.basis_z),
        )
        return (
            count = length(orbitals),
            orbitals = orbitals,
            overlap = _qwrg_diatomic_overlap_matrix(
                bundles.bundle_x,
                bundles.bundle_y,
                bundles.bundle_z,
            ),
        )
    end

    fixed_block = context.carried
    return (
        count = size(fixed_block.overlap, 1),
        fixed_centers = fixed_block.fixed_centers,
        overlap = fixed_block.overlap,
    )
end

function _qwrg_bond_aligned_molecular_residual_space(
    context::_QWRGBondAlignedBuildContext,
    carried_data,
    blocks,
    supplement::_BondAlignedDiatomicCartesianShellSupplement3D,
)
    supplement_owner_indices = Int[orbital.owner_nucleus_index for orbital in supplement.orbitals]
    if context.carried_space_kind == :direct_product
        return _qwrg_residual_space_by_owner(
            carried_data.overlap,
            blocks.overlap_ga,
            blocks.overlap_aa,
            supplement_owner_indices,
        )
    end

    contraction = something(context.contraction)
    overlap_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.overlap_ga)
    return _qwrg_residual_space_by_owner(
        carried_data.overlap,
        overlap_fg,
        blocks.overlap_aa,
        supplement_owner_indices,
    )
end

function _qwrg_bond_aligned_molecular_carried_one_body_blocks(
    context::_QWRGBondAlignedBuildContext,
    bundles,
    expansion::CoulombGaussianExpansion,
    nuclear_charges::AbstractVector{<:Real},
    resolved_nuclear_term_storage::Symbol,
    gausslet_backend::Symbol,
    use_by_center_final_mix::Bool,
)
    basis = context.parent_basis
    if context.carried_space_kind == :direct_product
        gausslet_kinetic = _qwrg_diatomic_kinetic_matrix(
            bundles.bundle_x,
            bundles.bundle_y,
            bundles.bundle_z,
        )
        carried_nuclear_one_body_by_center =
            resolved_nuclear_term_storage == :by_center ?
            _qwrg_bond_aligned_direct_contracted_nuclear_one_body_by_center(
                basis,
                _qwrg_full_cartesian_product_factorized_basis(
                    (
                        size(bundles.bundle_x.pgdg_intermediate.overlap, 1),
                        size(bundles.bundle_y.pgdg_intermediate.overlap, 1),
                        size(bundles.bundle_z.pgdg_intermediate.overlap, 1),
                    ),
                ),
                bundles.bundle_x,
                bundles.bundle_y,
                bundles.bundle_z,
                expansion,
            ) :
            nothing
        carried_one_body =
            use_by_center_final_mix ?
            nothing :
            _qwrg_diatomic_one_body_matrix(
                basis,
                bundles.bundle_x,
                bundles.bundle_y,
                bundles.bundle_z,
                expansion,
                nuclear_charges,
            )
        return (
            kinetic = gausslet_kinetic,
            nuclear_by_center = carried_nuclear_one_body_by_center,
            one_body = carried_one_body,
        )
    end

    fixed_block = context.carried
    contraction = something(context.contraction)
    timing_prefix = context.capabilities.timing_prefix
    fixed_kinetic = @timeg "$(timing_prefix).one_body.carried.kinetic" begin
        Matrix{Float64}(fixed_block.kinetic)
    end
    fixed_nuclear_one_body_by_center =
        resolved_nuclear_term_storage == :by_center ?
        let
            factorized_basis = @timeg "$(timing_prefix).one_body.carried.factorized_basis" begin
                _nested_factorized_parent_basis(fixed_block)
            end
            _qwrg_bond_aligned_direct_contracted_nuclear_one_body_by_center(
                basis,
                factorized_basis,
                bundles.bundle_x,
                bundles.bundle_y,
                bundles.bundle_z,
                expansion;
                timing_setup_label = "$(timing_prefix).one_body.carried.nuclear_setup",
                timing_contract_label = "$(timing_prefix).one_body.carried.nuclear_contract",
            )
        end :
        nothing
    fixed_one_body =
        use_by_center_final_mix ?
        nothing :
        @timeg "$(timing_prefix).one_body.carried.reassembly" begin
            parent_one_body = _qwrg_diatomic_one_body_matrix(
                basis,
                bundles.bundle_x,
                bundles.bundle_y,
                bundles.bundle_z,
                expansion,
                nuclear_charges,
            )
            matrix = Matrix{Float64}(transpose(contraction) * parent_one_body * contraction)
            0.5 .* (matrix .+ transpose(matrix))
        end
    return (
        kinetic = fixed_kinetic,
        nuclear_by_center = fixed_nuclear_one_body_by_center,
        one_body = fixed_one_body,
    )
end

function _qwrg_bond_aligned_molecular_coupling_one_body_blocks(
    context::_QWRGBondAlignedBuildContext,
    blocks,
    use_by_center_final_mix::Bool,
)
    if context.carried_space_kind == :direct_product
        return (
            kinetic = Matrix{Float64}(blocks.kinetic_ga),
            nuclear_by_center =
                use_by_center_final_mix ?
                [Matrix{Float64}(matrix) for matrix in blocks.nuclear_ga_by_center] :
                nothing,
            one_body =
                use_by_center_final_mix ?
                nothing :
                Matrix{Float64}(blocks.one_body_ga),
        )
    end

    contraction = something(context.contraction)
    return (
        kinetic = _qwrg_contract_parent_ga_matrix(contraction, blocks.kinetic_ga),
        nuclear_by_center =
            use_by_center_final_mix ?
            [
                _qwrg_contract_parent_ga_matrix(contraction, matrix) for
                matrix in blocks.nuclear_ga_by_center
            ] :
            nothing,
        one_body =
            use_by_center_final_mix ?
            nothing :
            _qwrg_contract_parent_ga_matrix(contraction, blocks.one_body_ga),
    )
end

function _qwrg_bond_aligned_molecular_supplement_one_body_blocks(
    blocks,
    use_by_center_final_mix::Bool,
)
    return (
        kinetic = Matrix{Float64}(blocks.kinetic_aa),
        nuclear_by_center =
            use_by_center_final_mix ?
            [Matrix{Float64}(matrix) for matrix in blocks.nuclear_aa_by_center] :
            nothing,
        one_body =
            use_by_center_final_mix ?
            nothing :
            Matrix{Float64}(blocks.one_body_aa),
    )
end

function _qwrg_bond_aligned_molecular_residual_centers(
    context::_QWRGBondAlignedBuildContext,
    carried_data,
    bundles,
    blocks,
    residual_data,
)
    if context.carried_space_kind == :direct_product
        x_carried = _qwrg_diatomic_gausslet_axis_matrix(
            bundles.bundle_x,
            bundles.bundle_y,
            bundles.bundle_z,
            :x,
        )
        y_carried = _qwrg_diatomic_gausslet_axis_matrix(
            bundles.bundle_x,
            bundles.bundle_y,
            bundles.bundle_z,
            :y,
        )
        z_carried = _qwrg_diatomic_gausslet_axis_matrix(
            bundles.bundle_x,
            bundles.bundle_y,
            bundles.bundle_z,
            :z,
        )
    else
        fixed_block = context.carried
        contraction = something(context.contraction)
        x_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.position_x_ga)
        y_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.position_y_ga)
        z_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.position_z_ga)
        x_raw = [Matrix{Float64}(fixed_block.position_x) x_fg; transpose(x_fg) blocks.position_x_aa]
        y_raw = [Matrix{Float64}(fixed_block.position_y) y_fg; transpose(y_fg) blocks.position_y_aa]
        z_raw = [Matrix{Float64}(fixed_block.position_z) z_fg; transpose(z_fg) blocks.position_z_aa]
        center_data = _qwrg_residual_center_data(
            residual_data.raw_overlap,
            x_raw,
            y_raw,
            z_raw,
            residual_data.raw_to_final,
            carried_data.count,
        )
        residual_centers = center_data.centers
        residual_widths = fill(NaN, size(residual_centers, 1), 3)
        center_data.overlap_error <= 1.0e-8 || throw(
            ArgumentError("$(context.capabilities.residual_center_label) requires an orthonormal residual block"),
        )
        return residual_centers, residual_widths
    end

    x_raw = [x_carried blocks.position_x_ga; transpose(blocks.position_x_ga) blocks.position_x_aa]
    y_raw = [y_carried blocks.position_y_ga; transpose(blocks.position_y_ga) blocks.position_y_aa]
    z_raw = [z_carried blocks.position_z_ga; transpose(blocks.position_z_ga) blocks.position_z_aa]
    center_data = _qwrg_residual_center_data(
        residual_data.raw_overlap,
        x_raw,
        y_raw,
        z_raw,
        residual_data.raw_to_final,
        carried_data.count,
    )
    residual_centers = center_data.centers
    residual_widths = fill(NaN, size(residual_centers, 1), 3)
    center_data.overlap_error <= 1.0e-8 || throw(
        ArgumentError("$(context.capabilities.residual_center_label) requires an orthonormal residual block"),
    )
    return residual_centers, residual_widths
end

function _qwrg_bond_aligned_molecular_interaction_matrix(
    context::_QWRGBondAlignedBuildContext,
    carried_data,
    bundles,
    residual_centers::AbstractMatrix{<:Real},
    expansion::CoulombGaussianExpansion,
)
    if context.carried_space_kind == :direct_product
        gausslet_interaction = _qwrg_diatomic_interaction_matrix(
            bundles.bundle_x,
            bundles.bundle_y,
            bundles.bundle_z,
            expansion,
        )
        return _qwrg_interaction_matrix_nearest(
            gausslet_interaction,
            carried_data.orbitals,
            residual_centers,
        )
    end

    fixed_block = context.carried
    fixed_interaction = _qwrg_fixed_block_interaction_matrix(fixed_block, expansion)
    return _qwrg_interaction_matrix_nearest(
        fixed_interaction,
        carried_data.fixed_centers,
        residual_centers,
    )
end

function _qwrg_bond_aligned_molecular_operator_orbitals(
    context::_QWRGBondAlignedBuildContext,
    carried_data,
    residual_centers::AbstractMatrix{<:Real},
    residual_widths::AbstractMatrix{<:Real},
    residual_nucleus_indices::AbstractVector{<:Integer},
)
    if context.carried_space_kind == :direct_product
        return _qwrg_orbital_data(
            carried_data.orbitals,
            residual_centers,
            residual_widths,
            residual_nucleus_indices,
        )
    end

    return _qwrg_orbital_data(
        carried_data.fixed_centers,
        residual_centers,
        residual_widths;
        fixed_kind = :nested_fixed,
        fixed_label_prefix = "nf",
        residual_nucleus_indices = residual_nucleus_indices,
    )
end

function _ordinary_cartesian_qiu_white_operators_bond_aligned_molecular(
    context::_QWRGBondAlignedBuildContext;
    nuclear_charges::AbstractVector{<:Real},
    nuclear_term_storage::Symbol,
    expansion::CoulombGaussianExpansion,
    interaction_treatment::Symbol,
    gausslet_backend::Symbol,
    timing::Bool,
)
    _validate_operator_route_backend(context, gausslet_backend)
    _validate_operator_route_interaction_treatment(context, interaction_treatment)
    _validate_bond_aligned_molecular_inputs(context, nuclear_charges)
    gaussian_data = something(context.gaussian_data)
    basis = context.parent_basis
    carried = context.carried
    timing_prefix = context.capabilities.timing_prefix
    resolved_nuclear_term_storage = _resolved_nuclear_term_storage(
        _validate_nuclear_term_storage(nuclear_term_storage),
        carried,
    )

    return _qwrg_capture_timeg_report(context.capabilities.timing_label, timing) do
        bundles = @timeg context.capabilities.bundles_label begin
            _qwrg_bond_aligned_axis_bundles(
                basis,
                expansion;
                gausslet_backend = gausslet_backend,
            )
        end

        supplement3d = _bond_aligned_diatomic_cartesian_shell_supplement_3d(gaussian_data)
        blocks = @timeg "$(timing_prefix).raw_blocks" begin
            _qwrg_diatomic_cartesian_shell_blocks_3d(
                bundles,
                supplement3d,
                basis,
                expansion,
                nuclear_charges,
            )
        end

        carried_data = _qwrg_bond_aligned_molecular_carried_data(context, bundles)
        residual_data = @timeg "$(timing_prefix).residual_space" begin
            _qwrg_bond_aligned_molecular_residual_space(
                context,
                carried_data,
                blocks,
                supplement3d,
            )
        end

        final_kinetic, final_nuclear_one_body_by_center, final_one_body = @timeg "$(timing_prefix).one_body" begin
            use_by_center_final_mix =
                gausslet_backend == :pgdg_localized_experimental ||
                resolved_nuclear_term_storage == :by_center

            carried_blocks = @timeg "$(timing_prefix).one_body.carried" begin
                _qwrg_bond_aligned_molecular_carried_one_body_blocks(
                    context,
                    bundles,
                    expansion,
                    nuclear_charges,
                    resolved_nuclear_term_storage,
                    gausslet_backend,
                    use_by_center_final_mix,
                )
            end

            coupling_blocks = @timeg "$(timing_prefix).one_body.coupling" begin
                _qwrg_bond_aligned_molecular_coupling_one_body_blocks(
                    context,
                    blocks,
                    use_by_center_final_mix,
                )
            end

            supplement_blocks = @timeg "$(timing_prefix).one_body.supplement" begin
                _qwrg_bond_aligned_molecular_supplement_one_body_blocks(
                    blocks,
                    use_by_center_final_mix,
                )
            end

            final_kinetic_local = if use_by_center_final_mix
                @timeg "$(timing_prefix).one_body.final_mix" begin
                    _qwrg_final_one_body_matrix_from_blocks(
                        carried_blocks.kinetic,
                        coupling_blocks.kinetic,
                        supplement_blocks.kinetic,
                        residual_data.raw_to_final,
                    )[2]
                end
            else
                _qwrg_final_one_body_matrix_from_blocks(
                    carried_blocks.kinetic,
                    coupling_blocks.kinetic,
                    supplement_blocks.kinetic,
                    residual_data.raw_to_final,
                )[2]
            end

            assembled_nuclear_one_body_by_center_local = if use_by_center_final_mix
                @timeg "$(timing_prefix).one_body.by_center_final_mix" begin
                    _qwrg_final_nuclear_one_body_by_center(
                        carried_blocks.nuclear_by_center,
                        coupling_blocks.nuclear_by_center,
                        supplement_blocks.nuclear_by_center,
                        residual_data.raw_to_final,
                    )
                end
            else
                nothing
            end

            final_one_body_local = if isnothing(assembled_nuclear_one_body_by_center_local)
                @timeg "$(timing_prefix).one_body.final_mix" begin
                    _qwrg_final_one_body_matrix_from_blocks(
                        carried_blocks.one_body,
                        coupling_blocks.one_body,
                        supplement_blocks.one_body,
                        residual_data.raw_to_final,
                    )[2]
                end
            else
                _assemble_one_body_hamiltonian(
                    final_kinetic_local,
                    assembled_nuclear_one_body_by_center_local,
                    nuclear_charges,
                )
            end

            final_nuclear_one_body_by_center_local =
                resolved_nuclear_term_storage == :by_center ?
                assembled_nuclear_one_body_by_center_local :
                nothing
            (
                final_kinetic_local,
                final_nuclear_one_body_by_center_local,
                final_one_body_local,
            )
        end

        residual_centers, residual_widths = @timeg "$(timing_prefix).centers" begin
            _qwrg_bond_aligned_molecular_residual_centers(
                context,
                carried_data,
                bundles,
                blocks,
                residual_data,
            )
        end

        interaction_matrix = @timeg "$(timing_prefix).interaction" begin
            _qwrg_bond_aligned_molecular_interaction_matrix(
                context,
                carried_data,
                bundles,
                residual_centers,
                expansion,
            )
        end

        return OrdinaryCartesianOperators3D(
            carried,
            gaussian_data,
            gausslet_backend,
            interaction_treatment,
            expansion,
            Matrix{Float64}(residual_data.final_overlap),
            final_one_body,
            Matrix{Float64}(0.5 .* (interaction_matrix .+ transpose(interaction_matrix))),
            _qwrg_bond_aligned_molecular_operator_orbitals(
                context,
                carried_data,
                residual_centers,
                residual_widths,
                residual_data.residual_nucleus_indices,
            ),
            carried_data.count,
            size(residual_centers, 1),
            Matrix{Float64}(residual_data.raw_to_final),
            Matrix{Float64}(residual_centers),
            Matrix{Float64}(residual_widths),
            Int[residual_data.residual_nucleus_indices...],
            Float64[Float64(value) for value in nuclear_charges],
            resolved_nuclear_term_storage == :by_center ? Matrix{Float64}(final_kinetic) : nothing,
            resolved_nuclear_term_storage == :by_center ? final_nuclear_one_body_by_center : nothing,
            resolved_nuclear_term_storage,
        )
    end
end

function _ordinary_cartesian_qiu_white_operators_pure_bond_aligned_nested(
    context::_QWRGBondAlignedBuildContext;
    nuclear_charges::AbstractVector{<:Real},
    nuclear_term_storage::Symbol,
    expansion::CoulombGaussianExpansion,
    interaction_treatment::Symbol,
    gausslet_backend::Symbol,
)
    fixed_block = context.carried
    basis = context.parent_basis
    contraction = something(context.contraction)
    resolved_nuclear_term_storage = _resolved_nuclear_term_storage(
        _validate_nuclear_term_storage(nuclear_term_storage),
        fixed_block,
    )
    length(nuclear_charges) == length(basis.nuclei) || throw(
        ArgumentError("$(context.capabilities.route_label) requires one nuclear charge per nucleus"),
    )

    bundles = _qwrg_bond_aligned_axis_bundles(
        basis,
        expansion;
        gausslet_backend = gausslet_backend,
    )

    # Nested fixed blocks carry packet-level kinetic directly. This path should
    # not reinterpret that payload as a late contraction of a separate parent
    # ordinary one-body build.
    kinetic_one_body = Matrix{Float64}(fixed_block.kinetic)
    direct_contracted_nuclear =
        resolved_nuclear_term_storage == :by_center ?
        _qwrg_bond_aligned_direct_contracted_nuclear_one_body_by_center(
            basis,
            _nested_factorized_parent_basis(fixed_block),
            bundles.bundle_x,
            bundles.bundle_y,
            bundles.bundle_z,
            expansion,
        ) :
        nothing
    nuclear_one_body_by_center =
        resolved_nuclear_term_storage == :by_center ?
        direct_contracted_nuclear :
        nothing
    one_body_hamiltonian =
        !isnothing(direct_contracted_nuclear) ?
        _assemble_one_body_hamiltonian(
            kinetic_one_body,
            direct_contracted_nuclear,
            nuclear_charges,
        ) :
        (
            isnothing(nuclear_one_body_by_center) ?
            _qwrg_diatomic_one_body_matrix(
                basis,
                bundles.bundle_x,
                bundles.bundle_y,
                bundles.bundle_z,
                expansion,
                nuclear_charges,
            ) :
            _assemble_one_body_hamiltonian(
                kinetic_one_body,
                nuclear_one_body_by_center,
                nuclear_charges,
            )
        )
    if isnothing(nuclear_one_body_by_center) && isnothing(direct_contracted_nuclear)
        one_body_hamiltonian = Matrix{Float64}(transpose(contraction) * one_body_hamiltonian * contraction)
        one_body_hamiltonian = 0.5 .* (one_body_hamiltonian .+ transpose(one_body_hamiltonian))
    end
    interaction_matrix = _qwrg_fixed_block_interaction_matrix(fixed_block, expansion)

    fixed_count = size(fixed_block.overlap, 1)
    zero_residual_centers = zeros(Float64, 0, 3)
    zero_residual_widths = zeros(Float64, 0, 3)
    return OrdinaryCartesianOperators3D(
        fixed_block,
        nothing,
        gausslet_backend,
        interaction_treatment,
        expansion,
        Matrix{Float64}(fixed_block.overlap),
        one_body_hamiltonian,
        Matrix{Float64}(0.5 .* (interaction_matrix .+ transpose(interaction_matrix))),
        _qwrg_orbital_data(
            fixed_block.fixed_centers,
            zero_residual_centers,
            zero_residual_widths;
            fixed_kind = :nested_fixed,
            fixed_label_prefix = "nf",
        ),
        fixed_count,
        0,
        Matrix{Float64}(I, fixed_count, fixed_count),
        zero_residual_centers,
        zero_residual_widths,
        Float64[Float64(value) for value in nuclear_charges],
        resolved_nuclear_term_storage == :by_center ? kinetic_one_body : nothing,
        resolved_nuclear_term_storage == :by_center ? nuclear_one_body_by_center : nothing,
        resolved_nuclear_term_storage,
    )
end

function _ordinary_cartesian_qiu_white_operators_pure_bond_aligned(
    context::_QWRGBondAlignedBuildContext;
    nuclear_charges::AbstractVector{<:Real},
    nuclear_term_storage::Symbol,
    expansion::CoulombGaussianExpansion,
    interaction_treatment::Symbol,
    gausslet_backend::Symbol,
    timing::Bool,
)
    _validate_operator_route_backend(context, gausslet_backend)
    _validate_operator_route_interaction_treatment(context, interaction_treatment)
    return _qwrg_capture_timeg_report(context.capabilities.timing_label, timing) do
        context.carried_space_kind == :direct_product &&
            return _ordinary_cartesian_qiu_white_operators_pure_bond_aligned_direct(
                context;
                nuclear_charges = nuclear_charges,
                nuclear_term_storage = nuclear_term_storage,
                expansion = expansion,
                interaction_treatment = interaction_treatment,
                gausslet_backend = gausslet_backend,
            )
        context.carried_space_kind == :nested_fixed_block &&
            return _ordinary_cartesian_qiu_white_operators_pure_bond_aligned_nested(
                context;
                nuclear_charges = nuclear_charges,
                nuclear_term_storage = nuclear_term_storage,
                expansion = expansion,
                interaction_treatment = interaction_treatment,
                gausslet_backend = gausslet_backend,
            )
        throw(
            ArgumentError(
                "pure bond-aligned ordinary_cartesian_qiu_white_operators does not support carried_space_kind = :$(context.carried_space_kind)",
            ),
        )
    end
end

function _qwrg_atomic_carried_data(
    context::_QWRGAtomicBuildContext,
    gausslet_bundle,
)
    if context.carried_space_kind == :direct_product
        gg_blocks = _qwrg_gausslet_1d_blocks(gausslet_bundle)
        orbitals = _mapped_cartesian_orbitals(gausslet_bundle.pgdg_intermediate.centers)
        count = length(orbitals)
        overlap = zeros(Float64, count, count)
        _qwrg_fill_product_matrix!(
            overlap,
            gg_blocks.overlap_gg,
            gg_blocks.overlap_gg,
            gg_blocks.overlap_gg,
        )
        return (
            count = count,
            orbitals = orbitals,
            overlap = overlap,
            gg_blocks = gg_blocks,
        )
    end

    fixed_block = context.carried
    return (
        count = size(fixed_block.overlap, 1),
        fixed_centers = fixed_block.fixed_centers,
        overlap = fixed_block.overlap,
    )
end

function _qwrg_atomic_residual_space(
    context::_QWRGAtomicBuildContext,
    carried_data,
    blocks,
    residual_keep_policy::Symbol,
    residual_keep_tol::Float64,
    residual_accept_tol::Float64,
)
    if context.carried_space_kind == :direct_product
        return _qwrg_residual_space(
            carried_data.overlap,
            blocks.overlap_ga,
            blocks.overlap_aa;
            keep_policy = residual_keep_policy,
            keep_abs_tol = residual_keep_tol,
            accept_tol = residual_accept_tol,
        )
    end

    contraction = something(context.contraction)
    overlap_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.overlap_ga)
    return _qwrg_residual_space(
        carried_data.overlap,
        overlap_fg,
        blocks.overlap_aa;
        keep_policy = residual_keep_policy,
        keep_abs_tol = residual_keep_tol,
        accept_tol = residual_accept_tol,
    )
end

function _qwrg_atomic_final_one_body(
    context::_QWRGAtomicBuildContext,
    carried_data,
    gausslet_bundle,
    blocks,
    expansion::CoulombGaussianExpansion,
    Z::Real,
    raw_to_final::AbstractMatrix{<:Real},
)
    one_body_aa = _qwrg_atomic_cartesian_one_body_aa(blocks, expansion; Z = Z)
    if context.carried_space_kind == :direct_product
        gg_blocks = carried_data.gg_blocks
        gausslet_one_body = _qwrg_gausslet_one_body_matrix(gg_blocks, expansion; Z = Z)
        one_body_ga = Matrix{Float64}(blocks.kinetic_ga)
        for term in eachindex(expansion.coefficients)
            one_body_ga .-= Float64(Z) * expansion.coefficients[term] .* blocks.factor_ga[term]
        end
        return _qwrg_one_body_matrices(
            gausslet_one_body,
            one_body_ga,
            one_body_aa,
            raw_to_final,
        )[2]
    end

    fixed_block = context.carried
    contraction = something(context.contraction)
    kinetic_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.kinetic_ga)
    factor_fg = _qwrg_contract_parent_ga_terms(contraction, blocks.factor_ga)
    one_body_fg = Matrix{Float64}(kinetic_fg)
    for term in eachindex(expansion.coefficients)
        one_body_fg .-= Float64(Z) * expansion.coefficients[term] .* factor_fg[term]
    end
    fixed_one_body = _qwrg_fixed_block_one_body_matrix(fixed_block, expansion; Z = Z)
    return _qwrg_one_body_matrices(
        fixed_one_body,
        one_body_fg,
        one_body_aa,
        raw_to_final,
    )[2]
end

function _qwrg_atomic_residual_centers_and_widths(
    context::_QWRGAtomicBuildContext,
    carried_data,
    gausslet_bundle,
    blocks,
    residual_data,
    residual_accept_tol::Float64,
    interaction_treatment::Symbol,
)
    if context.carried_space_kind == :direct_product
        gg_blocks = carried_data.gg_blocks
        x_gg = _qwrg_gausslet_axis_matrix(gg_blocks, :x)
        y_gg = _qwrg_gausslet_axis_matrix(gg_blocks, :y)
        z_gg = _qwrg_gausslet_axis_matrix(gg_blocks, :z)
        x_raw = [x_gg blocks.position_x_ga; transpose(blocks.position_x_ga) blocks.position_x_aa]
        y_raw = [y_gg blocks.position_y_ga; transpose(blocks.position_y_ga) blocks.position_y_aa]
        z_raw = [z_gg blocks.position_z_ga; transpose(blocks.position_z_ga) blocks.position_z_aa]
        center_data = _qwrg_residual_center_data(
            residual_data.raw_overlap,
            x_raw,
            y_raw,
            z_raw,
            residual_data.raw_to_final,
            carried_data.count,
        )
        residual_centers = center_data.centers
        residual_widths = fill(NaN, size(residual_centers, 1), 3)
        center_data.overlap_error <= residual_accept_tol || throw(
            ArgumentError("$(context.capabilities.residual_center_label) requires an orthonormal residual block"),
        )

        if interaction_treatment == :mwg
            gg_x2 = (
                overlap_gg = gg_blocks.overlap_gg,
                position_gg = gg_blocks.position_gg,
                x2_gg = Matrix{Float64}(_x2_matrix(gausslet_bundle.basis)),
            )
            x2_x_gg = _qwrg_gausslet_axis_matrix(gg_x2, :x; squared = true)
            x2_y_gg = _qwrg_gausslet_axis_matrix(gg_x2, :y; squared = true)
            x2_z_gg = _qwrg_gausslet_axis_matrix(gg_x2, :z; squared = true)
            x2_raw = [x2_x_gg blocks.x2_x_ga; transpose(blocks.x2_x_ga) blocks.x2_x_aa]
            y2_raw = [x2_y_gg blocks.x2_y_ga; transpose(blocks.x2_y_ga) blocks.x2_y_aa]
            z2_raw = [x2_z_gg blocks.x2_z_ga; transpose(blocks.x2_z_ga) blocks.x2_z_aa]
            moment_data = _qwrg_residual_moment_data(
                residual_data.raw_overlap,
                x_raw,
                x2_raw,
                y2_raw,
                z2_raw,
                center_data,
            )
            residual_centers = moment_data.centers
            residual_widths = moment_data.widths
        end
        return residual_centers, residual_widths
    end

    fixed_block = context.carried
    contraction = something(context.contraction)
    x_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.position_x_ga)
    y_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.position_y_ga)
    z_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.position_z_ga)
    x_raw = [Matrix{Float64}(fixed_block.position_x) x_fg; transpose(x_fg) blocks.position_x_aa]
    y_raw = [Matrix{Float64}(fixed_block.position_y) y_fg; transpose(y_fg) blocks.position_y_aa]
    z_raw = [Matrix{Float64}(fixed_block.position_z) z_fg; transpose(z_fg) blocks.position_z_aa]
    center_data = _qwrg_residual_center_data(
        residual_data.raw_overlap,
        x_raw,
        y_raw,
        z_raw,
        residual_data.raw_to_final,
        carried_data.count,
    )
    residual_centers = center_data.centers
    residual_widths = fill(NaN, size(residual_centers, 1), 3)
    center_data.overlap_error <= residual_accept_tol || throw(
        ArgumentError("$(context.capabilities.residual_center_label) requires an orthonormal residual block"),
    )
    return residual_centers, residual_widths
end

function _qwrg_atomic_interaction_matrix(
    context::_QWRGAtomicBuildContext,
    carried_data,
    gausslet_bundle,
    expansion::CoulombGaussianExpansion,
    residual_centers::AbstractMatrix{<:Real},
    residual_widths::AbstractMatrix{<:Real},
    interaction_treatment::Symbol,
)
    if context.carried_space_kind == :direct_product
        gausslet_interaction = _qwrg_gausslet_interaction_matrix(
            carried_data.gg_blocks,
            expansion,
        )
        if interaction_treatment == :ggt_nearest
            return _qwrg_interaction_matrix_nearest(
                gausslet_interaction,
                carried_data.orbitals,
                residual_centers,
            )
        end
        return _qwrg_interaction_matrix_mwg(
            gausslet_bundle,
            gausslet_interaction,
            expansion,
            residual_centers,
            residual_widths,
        )
    end

    fixed_block = context.carried
    fixed_interaction = _qwrg_fixed_block_interaction_matrix(fixed_block, expansion)
    return _qwrg_interaction_matrix_nearest(
        fixed_interaction,
        fixed_block.fixed_centers,
        residual_centers,
    )
end

function _qwrg_atomic_operator_orbitals(
    context::_QWRGAtomicBuildContext,
    carried_data,
    residual_centers::AbstractMatrix{<:Real},
    residual_widths::AbstractMatrix{<:Real},
)
    if context.carried_space_kind == :direct_product
        return _qwrg_orbital_data(
            carried_data.orbitals,
            residual_centers,
            residual_widths,
        )
    end

    fixed_block = context.carried
    return _qwrg_orbital_data(
        fixed_block.fixed_centers,
        residual_centers,
        residual_widths;
        fixed_kind = :nested_fixed,
        fixed_label_prefix = "nf",
    )
end

function _ordinary_cartesian_qiu_white_operators_atomic(
    context::_QWRGAtomicBuildContext;
    expansion::CoulombGaussianExpansion,
    Z::Real,
    interaction_treatment::Symbol,
    gausslet_backend::Symbol,
    residual_keep_policy::Symbol = :near_null_only,
    timing::Bool,
)
    _validate_operator_route_backend(context, gausslet_backend)
    _validate_operator_route_interaction_treatment(context, interaction_treatment)
    residual_keep_policy_value = _qwrg_atomic_residual_keep_policy(residual_keep_policy)
    residual_keep_tol = _qwrg_atomic_residual_keep_tol()
    residual_accept_tol = _qwrg_atomic_residual_accept_tol()
    parent_basis = context.parent_basis
    carried = context.carried

    return _qwrg_capture_timeg_report(context.capabilities.timing_label, timing) do
        gausslet_bundle = @timeg context.capabilities.bundle_label begin
            _mapped_ordinary_gausslet_1d_bundle(
                parent_basis;
                exponents = expansion.exponents,
                center = 0.0,
                backend = gausslet_backend,
            )
        end

        supplement3d = _atomic_cartesian_shell_supplement_3d(context.gaussian_data)
        blocks = @timeg context.capabilities.raw_blocks_label begin
            _qwrg_atomic_cartesian_blocks_3d(
                gausslet_bundle,
                supplement3d,
                expansion,
            )
        end

        carried_data = _qwrg_atomic_carried_data(context, gausslet_bundle)
        residual_data = @timeg context.capabilities.residual_space_label begin
            _qwrg_atomic_residual_space(
                context,
                carried_data,
                blocks,
                residual_keep_policy_value,
                residual_keep_tol,
                residual_accept_tol,
            )
        end
        timing && _qwrg_print_basis_counts(
            stdout,
            context.capabilities.carried_count_label,
            carried_data.count,
            residual_data.raw_to_final,
        )

        final_one_body = @timeg context.capabilities.one_body_label begin
            _qwrg_atomic_final_one_body(
                context,
                carried_data,
                gausslet_bundle,
                blocks,
                expansion,
                Z,
                residual_data.raw_to_final,
            )
        end

        residual_centers, residual_widths = @timeg context.capabilities.centers_label begin
            _qwrg_atomic_residual_centers_and_widths(
                context,
                carried_data,
                gausslet_bundle,
                blocks,
                residual_data,
                residual_accept_tol,
                interaction_treatment,
            )
        end

        interaction_matrix = @timeg context.capabilities.interaction_label begin
            _qwrg_atomic_interaction_matrix(
                context,
                carried_data,
                gausslet_bundle,
                expansion,
                residual_centers,
                residual_widths,
                interaction_treatment,
            )
        end

        return OrdinaryCartesianOperators3D(
            carried,
            context.gaussian_data,
            gausslet_backend,
            interaction_treatment,
            expansion,
            Matrix{Float64}(residual_data.final_overlap),
            final_one_body,
            Matrix{Float64}(0.5 .* (interaction_matrix .+ transpose(interaction_matrix))),
            _qwrg_atomic_operator_orbitals(
                context,
                carried_data,
                residual_centers,
                residual_widths,
            ),
            carried_data.count,
            size(residual_centers, 1),
            Matrix{Float64}(residual_data.raw_to_final),
            Matrix{Float64}(residual_centers),
            Matrix{Float64}(residual_widths),
        )
    end
end

"""
    ordinary_cartesian_qiu_white_operators(
        basis::MappedUniformBasis,
        gaussian_data::LegacyAtomicGaussianSupplement;
        expansion = coulomb_gaussian_expansion(doacc = false),
        Z = 2.0,
        interaction_treatment = :mwg,
        gausslet_backend = :numerical_reference,
        residual_keep_policy = :near_null_only,
        timing = false,
    )

Build the paper-faithful Qiu-White residual-Gaussian ordinary Cartesian
reference Hamiltonian.

This path is intentionally separate from the current COMX/localized hybrid
route. It:

- keeps the full 3D gausslet product basis fixed
- orthogonalizes the added 3D Gaussian supplement against that full 3D space
- builds the one-body matrices exactly in the raw gausslet-plus-GTO space
- keeps the two-electron interaction in the same two-index integral-diagonal
  approximation (IDA) representation used for the gausslet channel

Allowed `interaction_treatment` values are:

- `:ggt_nearest`
- `:mwg`

Set `timing = true` to print a TimeG phase tree for debugging the reference
implementation. This is intentionally narrow and does not add timing noise to
the broader library.

This is a reference path for validating the Qiu-White formulation, not yet a
claim that the ordinary branch is solver-ready.

This path now uses the explicit atomic-centered 3D Cartesian shell route for
all active atomic supplement content up to `lmax <= 2`, including pure `s`
shells.

`nested_profile` selects the shell/working-box geometry contract. It is
independent of supplement residual retention.

`residual_keep_policy = :near_null_only` is now the canonical atomic
supplement keep contract. It keeps orthogonalized residual directions with
residual-overlap eigenvalue `> 1e-7`.

`residual_keep_policy = :legacy_profile` remains accepted only as a
compatibility alias for `:near_null_only`.

After residual directions are selected, the retained residual block is
explicitly stabilized back to overlap-orthonormality before any downstream
center/width extraction or raw-to-final transforms are used. The atomic
residual-center acceptance gate on that stabilized block is `1e-7`.
"""
function ordinary_cartesian_qiu_white_operators(
    basis::MappedUniformBasis,
    gaussian_data::LegacyAtomicGaussianSupplement;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
    interaction_treatment::Symbol = :mwg,
    gausslet_backend::Symbol = :numerical_reference,
    residual_keep_policy::Symbol = :near_null_only,
    timing::Bool = false,
    )
    context = _normalized_atomic_build_context(basis, gaussian_data)
    return _ordinary_cartesian_qiu_white_operators_atomic(
        context;
        expansion = expansion,
        Z = Z,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        residual_keep_policy = residual_keep_policy,
        timing = timing,
    )
end

"""
    ordinary_cartesian_qiu_white_operators(
        fixed_block::_NestedFixedBlock3D,
        gaussian_data::LegacyAtomicGaussianSupplement;
        expansion = coulomb_gaussian_expansion(doacc = false),
        Z = 2.0,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        residual_keep_policy = :near_null_only,
        timing = false,
    )

Build the first nested fixed-block QW-PGDG consumer path.

This overload reuses the stabilized parent PGDG raw fixed-to-Gaussian blocks,
contracts them through the supplied shell-level fixed map, and then runs the
same downstream residual-space / one-body / nearest-GGT algebra as the
unnested QW-PGDG route.

This first adapter is intentionally narrow:

- it consumes an already-assembled nonseparable 3D fixed packet
- it keeps the shell packet as the fixed-fixed block directly
- it supports only `interaction_treatment = :ggt_nearest`

It now uses the explicit atomic-centered 3D Cartesian shell route for all
active atomic supplement content up to `lmax <= 2`, including pure `s`
shells.

`nested_profile` selects the shell/working-box geometry contract. It is
independent of supplement residual retention.

`residual_keep_policy = :near_null_only` is now the canonical atomic
supplement keep contract. It keeps orthogonalized residual directions with
residual-overlap eigenvalue `> 1e-7`.

`residual_keep_policy = :legacy_profile` remains accepted only as a
compatibility alias for `:near_null_only`.

After residual directions are selected, the retained residual block is
explicitly stabilized back to overlap-orthonormality before any downstream
center/width extraction or raw-to-final transforms are used. The atomic
residual-center acceptance gate on that stabilized block is `1e-7`.
"""
function ordinary_cartesian_qiu_white_operators(
    fixed_block::_NestedFixedBlock3D,
    gaussian_data::LegacyAtomicGaussianSupplement;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    residual_keep_policy::Symbol = :near_null_only,
    timing::Bool = false,
)
    context = _normalized_atomic_build_context(fixed_block, gaussian_data)
    return _ordinary_cartesian_qiu_white_operators_atomic(
        context;
        expansion = expansion,
        Z = Z,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        residual_keep_policy = residual_keep_policy,
        timing = timing,
    )
end

"""
    ordinary_cartesian_qiu_white_operators(
        fixed_block::_NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D};
        nuclear_charges = fill(1.0, length(fixed_block.parent_basis.nuclei)),
        expansion = coulomb_gaussian_expansion(doacc = false),
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        timing = false,
    )

Build the first bond-aligned diatomic nested fixed-block ordinary QW path.

This is intentionally narrower than the atomic hybrid route:

- it consumes the already-assembled bond-aligned diatomic `_NestedFixedBlock3D`
- it keeps the residual-Gaussian sector empty
- it supports only `interaction_treatment = :ggt_nearest`
- this pure prebuilt nested fixed-block route accepts
  `gausslet_backend = :pgdg_localized_experimental` when the nested
  representation remains a pure Cartesian parent space

The point of this first pass is to validate the diatomic fixed-block geometry
and transferred packet cleanly before molecular Gaussian completion is added.
"""
function ordinary_cartesian_qiu_white_operators(
    fixed_block::_NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D};
    nuclear_charges::AbstractVector{<:Real} = fill(1.0, length(fixed_block.parent_basis.nuclei)),
    nuclear_term_storage::Symbol = :auto,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    timing::Bool = false,
)
    context = _normalized_bond_aligned_build_context(fixed_block)
    return _ordinary_cartesian_qiu_white_operators_pure_bond_aligned(
        context;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        timing = timing,
    )
end

"""
    ordinary_cartesian_qiu_white_operators(
        fixed_block::_NestedFixedBlock3D{<:BondAlignedHomonuclearChainQWBasis3D};
        nuclear_charges = fill(1.0, length(fixed_block.parent_basis.nuclei)),
        expansion = coulomb_gaussian_expansion(doacc = false),
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        timing = false,
    )

Build the first experimental bond-aligned homonuclear-chain nested fixed-block
ordinary QW path.

This path is intentionally narrow:

- homonuclear chains only
- nested fixed-block route only
- no supplement route
- no residual-Gaussian sector
- no claim that the odd-chain split policy is settled globally
- this pure prebuilt nested fixed-block route accepts
  `gausslet_backend = :pgdg_localized_experimental` when the nested
  representation remains a pure Cartesian parent space
"""
function ordinary_cartesian_qiu_white_operators(
    fixed_block::_NestedFixedBlock3D{<:BondAlignedHomonuclearChainQWBasis3D};
    nuclear_charges::AbstractVector{<:Real} = fill(1.0, length(fixed_block.parent_basis.nuclei)),
    nuclear_term_storage::Symbol = :auto,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    timing::Bool = false,
)
    context = _normalized_bond_aligned_build_context(fixed_block)
    return _ordinary_cartesian_qiu_white_operators_pure_bond_aligned(
        context;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        timing = timing,
    )
end

"""
    ordinary_cartesian_qiu_white_operators(
        fixed_block::_NestedFixedBlock3D{<:AxisAlignedHomonuclearSquareLatticeQWBasis3D};
        nuclear_charges = fill(1.0, length(fixed_block.parent_basis.nuclei)),
        expansion = coulomb_gaussian_expansion(doacc = false),
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        timing = false,
    )

Build the first experimental nested fixed-block ordinary-QW square-lattice
Hamiltonian.

This path is intentionally narrow:

- homonuclear square lattices only
- nested fixed-block route only
- no supplement route
- no residual-Gaussian sector
- no claim that the planar split policy is settled globally
- this pure prebuilt nested fixed-block route accepts
  `gausslet_backend = :pgdg_localized_experimental` when the nested
  representation remains a pure Cartesian parent space
"""
function ordinary_cartesian_qiu_white_operators(
    fixed_block::_NestedFixedBlock3D{<:AxisAlignedHomonuclearSquareLatticeQWBasis3D};
    nuclear_charges::AbstractVector{<:Real} = fill(1.0, length(fixed_block.parent_basis.nuclei)),
    nuclear_term_storage::Symbol = :auto,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    timing::Bool = false,
)
    context = _normalized_bond_aligned_build_context(fixed_block)
    return _ordinary_cartesian_qiu_white_operators_pure_bond_aligned(
        context;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        timing = timing,
    )
end


"""
    ordinary_cartesian_qiu_white_operators(
        basis::BondAlignedDiatomicQWBasis3D,
        gaussian_data::LegacyBondAlignedDiatomicGaussianSupplement;
        nuclear_charges = fill(1.0, length(basis.nuclei)),
        expansion = coulomb_gaussian_expansion(doacc = false),
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        timing = false,
    )

Build the first bond-aligned diatomic ordinary QW route with a true molecular
supplement and residual-Gaussian completion.

This first molecular supplement pass is intentionally narrow:

- one bond-aligned diatomic basis
- one explicit two-center molecular shell supplement built from a named atomic
  basis
- only `interaction_treatment = :ggt_nearest`
- this direct-product molecular supplement route now accepts
  `gausslet_backend = :pgdg_localized_experimental` on the carried GG
  one-body backbone only; GA/AA supplement closure remains unchanged
"""
function ordinary_cartesian_qiu_white_operators(
    basis::BondAlignedDiatomicQWBasis3D,
    gaussian_data::Union{
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    };
    nuclear_charges::AbstractVector{<:Real} = fill(1.0, length(basis.nuclei)),
    nuclear_term_storage::Symbol = :auto,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    timing::Bool = false,
)
    context = _normalized_bond_aligned_build_context(basis, gaussian_data)
    return _ordinary_cartesian_qiu_white_operators_bond_aligned_molecular(
        context;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        timing = timing,
    )
end

"""
    ordinary_cartesian_qiu_white_operators(
        fixed_block::_NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D},
        gaussian_data::Union{
            LegacyBondAlignedDiatomicGaussianSupplement,
            LegacyBondAlignedHeteronuclearGaussianSupplement,
        };
        nuclear_charges = fill(1.0, length(fixed_block.parent_basis.nuclei)),
        expansion = coulomb_gaussian_expansion(doacc = false),
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        timing = false,
    )

Build the first bond-aligned diatomic nested fixed-block QW route with a true
molecular supplement and residual-Gaussian completion.

This nested molecular supplement route now accepts
`gausslet_backend = :pgdg_localized_experimental` only when the carried nested
fixed block still represents a pure Cartesian parent space; the widened scope
is the fixed-space GG one-body backbone, while GA/AA supplement closure
remains unchanged.
"""
function ordinary_cartesian_qiu_white_operators(
    fixed_block::_NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D},
    gaussian_data::Union{
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    };
    nuclear_charges::AbstractVector{<:Real} = fill(1.0, length(fixed_block.parent_basis.nuclei)),
    nuclear_term_storage::Symbol = :auto,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    timing::Bool = false,
)
    context = _normalized_bond_aligned_build_context(fixed_block, gaussian_data)
    return _ordinary_cartesian_qiu_white_operators_bond_aligned_molecular(
        context;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        timing = timing,
    )
end

"""
    ordinary_cartesian_product_operators(...)

Compatibility-preserving clearer name for the direct-product Cartesian routes
within the ordinary Cartesian operator family. This includes both pure
Cartesian product routes and direct-product residual-Gaussian hybrid routes.
"""
function ordinary_cartesian_product_operators(
    basis::AbstractBondAlignedOrdinaryQWBasis3D;
    kwargs...,
)
    return ordinary_cartesian_qiu_white_operators(basis; kwargs...)
end

function ordinary_cartesian_product_operators(
    basis::MappedUniformBasis,
    gaussian_data::LegacyAtomicGaussianSupplement;
    kwargs...,
)
    return ordinary_cartesian_qiu_white_operators(basis, gaussian_data; kwargs...)
end

function ordinary_cartesian_product_operators(
    basis::BondAlignedDiatomicQWBasis3D,
    gaussian_data::Union{
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    };
    kwargs...,
)
    return ordinary_cartesian_qiu_white_operators(basis, gaussian_data; kwargs...)
end

"""
    nested_cartesian_operators(...)

Compatibility-preserving clearer name for the nested fixed-block Cartesian
routes within the ordinary Cartesian operator family. This includes both pure
nested fixed-block routes and nested residual-Gaussian hybrid routes.
"""
function nested_cartesian_operators(
    fixed_block::Union{
        _NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D},
        _NestedFixedBlock3D{<:BondAlignedHomonuclearChainQWBasis3D},
        _NestedFixedBlock3D{<:AxisAlignedHomonuclearSquareLatticeQWBasis3D},
    };
    kwargs...,
)
    return ordinary_cartesian_qiu_white_operators(fixed_block; kwargs...)
end

function nested_cartesian_operators(
    fixed_block::_NestedFixedBlock3D,
    gaussian_data::LegacyAtomicGaussianSupplement;
    kwargs...,
)
    return ordinary_cartesian_qiu_white_operators(fixed_block, gaussian_data; kwargs...)
end

function nested_cartesian_operators(
    fixed_block::_NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D},
    gaussian_data::Union{
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    };
    kwargs...,
)
    return ordinary_cartesian_qiu_white_operators(fixed_block, gaussian_data; kwargs...)
end
