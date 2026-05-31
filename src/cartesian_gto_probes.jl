function _gto_working_representation(working)
    representation =
        working isa CartesianBasisRepresentation3D ? working :
        working isa MappedUniformBasis ? _cartesian_direct_product_representation(working) :
        basis_representation(working)
    representation isa CartesianBasisRepresentation3D || throw(
        ArgumentError(
            "gto_overlap_matrix requires a Cartesian working-basis representation; got $(typeof(representation)) from $(typeof(working))",
        ),
    )
    _cartesian_supports_gto_overlap(representation) || throw(
        ArgumentError(
            "gto_overlap_matrix supports explicit Cartesian product, nested fixed-block, and hybrid residual-Gaussian Cartesian working bases with exact factorized sidecars or dense parent analytic GTO sidecars; got parent_kind :$(representation.metadata.parent_kind)",
        ),
    )
    return representation
end

function _cartesian_supports_gto_overlap(
    representation::CartesianBasisRepresentation3D,
)
    representation.metadata.parent_kind == :cartesian_product_basis && return true
    representation.metadata.parent_kind == :cartesian_plus_supplement_raw || return false
    hasproperty(representation.parent_data, :cartesian_parent_representation) || return false
    hasproperty(representation.parent_data, :supplement_representation) || return false
    parent_representation = representation.parent_data.cartesian_parent_representation
    supplement_representation = representation.parent_data.supplement_representation
    parent_representation isa CartesianBasisRepresentation3D || return false
    supplement_representation isa CartesianGaussianShellSupplementRepresentation3D || return false
    parent_representation.metadata.parent_kind == :cartesian_product_basis || return false
    supplement_representation.supplement_kind in (
        :atomic_cartesian_shell,
        :bond_aligned_diatomic_cartesian_shell,
        :bond_aligned_heteronuclear_cartesian_shell,
    ) || return false
    all(
        orbital -> orbital.primitive_normalization == :axiswise_normalized_cartesian_gaussian,
        supplement_representation.orbitals,
    ) || return false
    return (
        _cartesian_supports_exact_hybrid_overlap(representation) ||
        hasproperty(representation.parent_data, :cartesian_probe_overlap_kind)
    )
end

function _gto_probe_representation(probes)
    representation =
        probes isa CartesianGaussianShellSupplementRepresentation3D ? probes :
        basis_representation(probes)
    representation isa CartesianGaussianShellSupplementRepresentation3D || throw(
        ArgumentError(
            "gto_overlap_matrix requires a legacy Gaussian probe family or CartesianGaussianShellSupplementRepresentation3D; got $(typeof(representation)) from $(typeof(probes))",
        ),
    )
    return representation
end

const _PQS_SOURCE_BOX_GTO_CROSS_OVERLAP_UNSUPPORTED_TERMS = (
    :shell_projection,
    :lowdin_cleanup,
    :support_coefficient_matrix,
    :retained_pqs_weights,
    :ida_weight_division,
)

function _pqs_source_box_gto_axis_projection(
    axis_representation::BasisRepresentation1D,
    orbital::CartesianGaussianShellOrbitalRepresentation3D,
    axis::Symbol,
    interval::UnitRange{Int},
    coefficients::AbstractMatrix{<:Real},
)
    axis_cross = _cartesian_basis_supplement_axis_primitive_cross(
        axis_representation,
        orbital,
        axis,
    )
    length(interval) == size(coefficients, 1) || throw(
        DimensionMismatch(
            "PQS source-box GTO projection interval length $(length(interval)) does not match local coefficient row count $(size(coefficients, 1)) on axis $(axis)",
        ),
    )
    first(interval) >= 1 && last(interval) <= size(axis_cross, 1) || throw(
        BoundsError(axes(axis_cross, 1), interval),
    )
    return Matrix{Float64}(transpose(Matrix{Float64}(coefficients)) * axis_cross[interval, :])
end

function _pqs_source_box_gto_cross_overlap_shadow(
    raw_product_box_plan,
    parent_representation::CartesianBasisRepresentation3D,
    probes;
    provenance = :pqs_source_box_gto_cross_overlap_shadow,
)
    raw_product_box_plan.representation == :orthogonal_raw_product_box || throw(
        ArgumentError("PQS source-box GTO shadow requires a raw product-box plan"),
    )
    parent_representation.metadata.parent_kind == :cartesian_product_basis || throw(
        ArgumentError(
            "PQS source-box GTO shadow requires a Cartesian-product parent representation; got :$(parent_representation.metadata.parent_kind)",
        ),
    )
    probe_representation = _gto_probe_representation(probes)
    source_mode_dims = raw_product_box_plan.source_mode_dims
    axis_interval_lengths =
        ntuple(axis -> length(raw_product_box_plan.axis_intervals[axis]), 3)
    axis_coefficient_shapes = ntuple(
        axis -> size(raw_product_box_plan.axis_local_coefficients[axis]),
        3,
    )
    boundary_mode_count =
        length(raw_product_box_plan.boundary_selector.mode_indices)
    gto_count = length(probe_representation.orbitals)
    cross_overlap = zeros(Float64, boundary_mode_count, gto_count)
    axis_names = (:x, :y, :z)
    axis_representations = parent_representation.axis_representations

    for (column, orbital) in pairs(probe_representation.orbitals)
        axis_projections = ntuple(
            axis -> _pqs_source_box_gto_axis_projection(
                getproperty(axis_representations, axis_names[axis]),
                orbital,
                axis_names[axis],
                raw_product_box_plan.axis_intervals[axis],
                raw_product_box_plan.axis_local_coefficients[axis],
            ),
            3,
        )
        @inbounds for (row, mode) in pairs(raw_product_box_plan.boundary_selector.mode_indices)
            value = 0.0
            for primitive in eachindex(orbital.coefficients)
                value +=
                    Float64(orbital.coefficients[primitive]) *
                    axis_projections[1][mode[1], primitive] *
                    axis_projections[2][mode[2], primitive] *
                    axis_projections[3][mode[3], primitive]
            end
            cross_overlap[row, column] = value
        end
    end

    diagnostics = (
        path = :pqs_source_box_gto_cross_overlap_shadow,
        raw_plan_first_path = true,
        descriptor_adapter = false,
        private_shadow_only = true,
        pqs_representation = :mode_selected_raw_product_box,
        source_box_mode_span = true,
        raw_product_box_plan_used = true,
        shared_raw_product_box_plan_available =
            raw_product_box_plan.shared_raw_product_box_plan_available,
        shared_raw_product_box_plan_used =
            raw_product_box_plan.shared_raw_product_box_plan_used,
        source_mode_ordering = raw_product_box_plan.source_mode_ordering,
        source_mode_dims = source_mode_dims,
        axis_interval_lengths = axis_interval_lengths,
        axis_coefficient_shapes = axis_coefficient_shapes,
        boundary_mode_count = boundary_mode_count,
        gto_count = gto_count,
        primitive_axis_overlap_source =
            :_cartesian_basis_supplement_axis_primitive_cross,
        projected_1d_axis_tables_used = true,
        product_box_column_selection_reference = true,
        numerical_reference_fallback = false,
        final_basis_handoff_adoption = false,
        shell_projection_used = false,
        lowdin_cleanup_used = false,
        support_coefficient_matrix_used = false,
        support_local_pqs_oracle_used = false,
        retained_weight_semantics = :not_positive_quadrature_weights,
        retained_pqs_weights_used = false,
        retained_pqs_weights_positive_checked = false,
        ida_weight_division_allowed = false,
        packet_adoption = false,
        fixed_block_routing = false,
        qwhamiltonian_consumes = false,
        public_default_consumes = false,
        cr2_science_status_changed = false,
        local_ecp_gaussian_mwg_implemented = false,
        generic_retained_unit_framework = false,
        unsupported_terms = _PQS_SOURCE_BOX_GTO_CROSS_OVERLAP_UNSUPPORTED_TERMS,
        output_finite = all(isfinite, cross_overlap),
        max_abs_cross_overlap =
            isempty(cross_overlap) ? 0.0 : maximum(abs, cross_overlap),
        provenance = provenance,
    )
    return (
        path = :pqs_source_box_gto_cross_overlap_shadow,
        cross_overlap = cross_overlap,
        probe_representation = probe_representation,
        diagnostics = diagnostics,
    )
end

function _pqs_source_box_gto_cross_overlap_shadow(
    descriptor::_CartesianNestedProjectedQShellStagedUnitDescriptor3D,
    parent_representation::CartesianBasisRepresentation3D,
    probes;
    provenance = :pqs_source_box_gto_cross_overlap_shadow,
)
    raw_plan = CartesianContractedParentMetrics._pqs_raw_product_box_plan(
        descriptor,
    )
    return _pqs_source_box_gto_cross_overlap_shadow(
        raw_plan,
        parent_representation,
        probes;
        provenance = provenance,
    )
end

function _gto_block_indices(
    final_dimension::Int,
    block_indices::Nothing,
)
    return nothing
end

function _gto_block_indices(
    final_dimension::Int,
    block_indices::AbstractVector{<:Integer},
)
    indices = Int[index for index in block_indices]
    all(index -> 1 <= index <= final_dimension, indices) || throw(
        BoundsError(1:final_dimension, indices),
    )
    return indices
end

function _gto_block_indices(final_dimension::Int, block_indices)
    throw(
        ArgumentError(
            "block_indices must be nothing or an explicit vector of basis indices; got $(typeof(block_indices))",
        ),
    )
end

function _gto_weight_vector(weights::Symbol, norbitals::Int)
    weights == :uniform && return ones(Float64, norbitals)
    weights == :shell_equalized && throw(
        ArgumentError(
            "gto_occupancy_matrix weights = :shell_equalized is not implemented yet; pass :uniform or an explicit weight vector",
        ),
    )
    throw(
        ArgumentError(
            "unsupported gto_occupancy_matrix weights :$(weights); supported first-pass modes are :uniform and explicit vectors",
        ),
    )
end

function _gto_weight_vector(weights::AbstractVector{<:Real}, norbitals::Int)
    length(weights) == norbitals || throw(
        DimensionMismatch(
            "explicit GTO weights length $(length(weights)) does not match probe dimension $(norbitals)",
        ),
    )
    values = Float64[Float64(weight) for weight in weights]
    all(isfinite, values) || throw(ArgumentError("explicit GTO weights must be finite"))
    return values
end

function _gto_weight_vector(weights, norbitals::Int)
    throw(
        ArgumentError(
            "gto_occupancy_matrix weights must be :uniform or an explicit real vector; got $(typeof(weights))",
        ),
    )
end

function _cartesian_exact_cartesian_probe_cross(
    raw,
    probes::CartesianGaussianShellSupplementRepresentation3D,
)
    hasproperty(raw, :exact_cartesian_supplement_overlap) || return nothing
    exact_cross = raw.exact_cartesian_supplement_overlap
    exact_cross === nothing && return nothing
    _cartesian_same_supplement_raw_identity(raw.supplement_representation, probes) || return nothing
    return Matrix{Float64}(exact_cross)
end

function _cartesian_exact_supplement_probe_cross(
    raw,
    probes::CartesianGaussianShellSupplementRepresentation3D,
)
    hasproperty(raw, :exact_supplement_overlap) || return nothing
    exact_cross = raw.exact_supplement_overlap
    exact_cross === nothing && return nothing
    _cartesian_same_supplement_raw_identity(raw.supplement_representation, probes) || return nothing
    return Matrix{Float64}(exact_cross)
end

function _cartesian_gto_factorized_parent_basis(raw)
    hasproperty(raw, :factorized_cartesian_parent_basis) || return nothing
    factorized = raw.factorized_cartesian_parent_basis
    factorized isa _CartesianNestedFactorizedBasis3D && return factorized
    return nothing
end

function _cartesian_atomic_qw_dense_parent_probe_overlap(
    raw,
    probes::CartesianGaussianShellSupplementRepresentation3D,
)
    hasproperty(raw, :cartesian_probe_overlap_kind) || return nothing
    raw.cartesian_probe_overlap_kind == :atomic_qw_dense_parent || return nothing
    parent_basis = raw.cartesian_probe_parent_basis
    parent_basis isa MappedUniformBasis || throw(
        ArgumentError(
            "atomic QW dense-parent GTO overlap fallback requires a MappedUniformBasis parent; got $(typeof(parent_basis))",
        ),
    )
    gausslet_bundle = _mapped_ordinary_gausslet_1d_bundle(
        parent_basis;
        exponents = raw.cartesian_probe_expansion.exponents,
        center = 0.0,
        backend = raw.cartesian_probe_gausslet_backend,
    )
    raw_blocks = _qwrg_atomic_cartesian_blocks_3d(
        gausslet_bundle,
        probes,
        raw.cartesian_probe_expansion,
    )
    parent_to_fixed_coefficients =
        raw.cartesian_representation.coefficient_matrix === nothing ?
        Matrix{Float64}(
            I,
            raw.cartesian_representation.metadata.final_dimension,
            raw.cartesian_representation.metadata.final_dimension,
        ) :
        Matrix{Float64}(raw.cartesian_representation.coefficient_matrix)
    return Matrix{Float64}(transpose(parent_to_fixed_coefficients) * raw_blocks.overlap_ga)
end

function _cartesian_cartesian_probe_overlap(
    raw,
    probes::CartesianGaussianShellSupplementRepresentation3D,
)
    isempty(probes.orbitals) &&
        return zeros(Float64, raw.cartesian_representation.metadata.final_dimension, 0)
    exact_cross = _cartesian_exact_cartesian_probe_cross(raw, probes)
    exact_cross !== nothing && return exact_cross
    try
        return _cartesian_basis_supplement_cross(raw.cartesian_representation, probes)
    catch error
        if error isa ArgumentError
            dense_parent_cross = _cartesian_atomic_qw_dense_parent_probe_overlap(raw, probes)
            dense_parent_cross !== nothing && return dense_parent_cross
            factorized = _cartesian_gto_factorized_parent_basis(raw)
            factorized === nothing && throw(
                ArgumentError(
                    "gto_overlap_matrix cannot build analytic Cartesian-GTO overlap for this working/probe pair: dense parent overlap failed with $(sprint(showerror, error)), and no exact factorized sidecar or dense QW parent sidecar is available",
                ),
            )
            return _cartesian_factorized_basis_supplement_cross(
                factorized,
                raw.cartesian_representation,
                probes,
                raw,
            )
        end
        rethrow()
    end
end

function _gto_raw_components(
    representation::CartesianBasisRepresentation3D,
)
    if representation.metadata.parent_kind == :cartesian_product_basis
        nraw = representation.metadata.final_dimension
        factorized = _cartesian_optional_factorized_parent_basis(representation)
        return merge(
            (
                cartesian_representation = representation,
                supplement_representation = _cartesian_empty_supplement_representation(),
                raw_to_final = Matrix{Float64}(I, nraw, nraw),
                cartesian_supplement_axis_tables = nothing,
                exact_cartesian_supplement_overlap = nothing,
                exact_supplement_overlap = nothing,
            ),
            isnothing(factorized) ? (;) : (factorized_cartesian_parent_basis = factorized,),
        )
    elseif representation.metadata.parent_kind == :cartesian_plus_supplement_raw
        _cartesian_supports_gto_overlap(representation) || throw(
            ArgumentError(
                "gto_overlap_matrix requires hybrid residual-Gaussian representations to expose an exact factorized sidecar or dense parent analytic GTO sidecar",
            ),
        )
        raw_to_final = representation.coefficient_matrix === nothing ?
            throw(
                ArgumentError(
                    "hybrid Cartesian GTO overlap requires an explicit raw_to_final coefficient matrix",
                ),
            ) :
            Matrix{Float64}(representation.coefficient_matrix)
        return (
            cartesian_representation = representation.parent_data.cartesian_parent_representation,
            supplement_representation = representation.parent_data.supplement_representation,
            raw_to_final = raw_to_final,
            factorized_cartesian_parent_basis =
                hasproperty(representation.parent_data, :factorized_cartesian_parent_basis) ?
                representation.parent_data.factorized_cartesian_parent_basis : nothing,
            cartesian_supplement_axis_tables =
                hasproperty(representation.parent_data, :cartesian_supplement_axis_tables) ?
                representation.parent_data.cartesian_supplement_axis_tables : nothing,
            exact_cartesian_supplement_overlap =
                hasproperty(representation.parent_data, :exact_cartesian_supplement_overlap) ?
                representation.parent_data.exact_cartesian_supplement_overlap : nothing,
            exact_supplement_overlap =
                hasproperty(representation.parent_data, :exact_supplement_overlap) ?
                representation.parent_data.exact_supplement_overlap : nothing,
            cartesian_probe_overlap_kind =
                hasproperty(representation.parent_data, :cartesian_probe_overlap_kind) ?
                representation.parent_data.cartesian_probe_overlap_kind : nothing,
            cartesian_probe_parent_basis =
                hasproperty(representation.parent_data, :cartesian_probe_parent_basis) ?
                representation.parent_data.cartesian_probe_parent_basis : nothing,
            cartesian_probe_expansion =
                hasproperty(representation.parent_data, :cartesian_probe_expansion) ?
                representation.parent_data.cartesian_probe_expansion : nothing,
            cartesian_probe_gausslet_backend =
                hasproperty(representation.parent_data, :cartesian_probe_gausslet_backend) ?
                representation.parent_data.cartesian_probe_gausslet_backend : nothing,
        )
    end
    throw(
        ArgumentError(
            "gto_overlap_matrix does not support Cartesian parent kind :$(representation.metadata.parent_kind)",
        ),
    )
end

function _cartesian_supplement_probe_overlap(
    raw,
    probes::CartesianGaussianShellSupplementRepresentation3D,
)
    isempty(raw.supplement_representation.orbitals) &&
        return zeros(Float64, 0, length(probes.orbitals))
    exact_cross = _cartesian_exact_supplement_probe_cross(raw, probes)
    exact_cross !== nothing && return exact_cross
    return _cartesian_supplement_cross_overlap(raw.supplement_representation, probes)
end

function _cartesian_gto_hybrid_audit_diagnostics(
    representation::CartesianBasisRepresentation3D,
)
    parent_data = representation.parent_data
    if parent_data !== nothing &&
            hasproperty(parent_data, :hybrid_overlap_audit_contract)
        return (
            dense_audit_materialization_used = true,
            dense_audit_materialization_contract =
                parent_data.hybrid_overlap_audit_contract,
            dense_audit_coefficient_scope =
                hasproperty(parent_data, :hybrid_overlap_audit_coefficient_scope) ?
                parent_data.hybrid_overlap_audit_coefficient_scope : nothing,
            dense_audit_cross_overlap =
                hasproperty(parent_data, :hybrid_overlap_audit_cross_overlap) ?
                parent_data.hybrid_overlap_audit_cross_overlap : nothing,
            dense_audit_parent_parent_operator =
                hasproperty(parent_data, :hybrid_overlap_audit_parent_parent_operator) ?
                parent_data.hybrid_overlap_audit_parent_parent_operator : false,
        )
    end
    return (
        dense_audit_materialization_used = false,
        dense_audit_materialization_contract = nothing,
        dense_audit_coefficient_scope = nothing,
        dense_audit_cross_overlap = nothing,
        dense_audit_parent_parent_operator = false,
    )
end

function _cartesian_final_gto_cross_overlap_handoff(
    fixed_working,
    final_supplement,
    raw_to_final::AbstractMatrix{<:Real},
    probes;
    provenance::Symbol = :private_final_gto_cross_overlap_handoff,
)
    fixed_representation = _gto_working_representation(fixed_working)
    supplement_representation = _gto_probe_representation(final_supplement)
    probe_representation = _gto_probe_representation(probes)
    raw_to_final_matrix = Matrix{Float64}(raw_to_final)
    fixed_gto_overlap = Matrix{Float64}(gto_overlap_matrix(fixed_representation, probe_representation))
    supplement_gto_overlap = Matrix{Float64}(
        _cartesian_supplement_cross_overlap(
            supplement_representation,
            probe_representation,
        ),
    )
    raw_dimension = size(fixed_gto_overlap, 1) + size(supplement_gto_overlap, 1)
    size(raw_to_final_matrix, 1) == raw_dimension || throw(
        DimensionMismatch(
            "final-basis/GTO handoff raw_to_final rows $(size(raw_to_final_matrix, 1)) must match fixed/GTO rows $(size(fixed_gto_overlap, 1)) plus supplement/GTO rows $(size(supplement_gto_overlap, 1))",
        ),
    )
    size(fixed_gto_overlap, 2) == size(supplement_gto_overlap, 2) || throw(
        DimensionMismatch(
            "fixed/GTO and supplement/GTO overlaps must have the same GTO dimension; got $(size(fixed_gto_overlap, 2)) and $(size(supplement_gto_overlap, 2))",
        ),
    )
    raw_gto_overlap = [fixed_gto_overlap; supplement_gto_overlap]
    final_gto_overlap = Matrix{Float64}(transpose(raw_to_final_matrix) * raw_gto_overlap)
    audit = _cartesian_gto_hybrid_audit_diagnostics(fixed_representation)
    diagnostics = merge(
        (
            cross_overlap_contract = :final_basis_uses_cross_overlap_only,
            self_overlaps_downstream_data = false,
            fixed_dimension = size(fixed_gto_overlap, 1),
            supplement_dimension = size(supplement_gto_overlap, 1),
            residual_count = size(supplement_gto_overlap, 1),
            raw_dimension = raw_dimension,
            final_dimension = size(final_gto_overlap, 1),
            gto_dimension = size(final_gto_overlap, 2),
            orientation = :final_by_gto,
            cross_overlap_label = :S_final_gto,
            transfer_convention = :C_final_equals_S_final_gto_times_C_gto,
            raw_piece_order = (:fixed_gto_overlap, :supplement_gto_overlap),
            final_transform = :transpose_raw_to_final_times_raw_gto_overlap,
            fixed_overlap_source = :gto_overlap_matrix,
            supplement_overlap_source = :_cartesian_supplement_cross_overlap,
            final_self_overlap_used = false,
            gto_self_overlap_used = false,
            output_finite = all(isfinite, final_gto_overlap),
            max_abs_cross_overlap =
                isempty(final_gto_overlap) ? 0.0 : maximum(abs, final_gto_overlap),
            provenance = provenance,
        ),
        audit,
    )
    return (
        cross_overlap = final_gto_overlap,
        final_gto_overlap = final_gto_overlap,
        raw_gto_overlap = raw_gto_overlap,
        fixed_gto_overlap = fixed_gto_overlap,
        supplement_gto_overlap = supplement_gto_overlap,
        fixed_representation = fixed_representation,
        supplement_representation = supplement_representation,
        probe_representation = probe_representation,
        diagnostics = diagnostics,
    )
end

function _gto_overlap_matrix(
    working::CartesianBasisRepresentation3D,
    probes::CartesianGaussianShellSupplementRepresentation3D,
)
    raw = _gto_raw_components(working)
    cg = _cartesian_cartesian_probe_overlap(raw, probes)
    gg = _cartesian_supplement_probe_overlap(raw, probes)
    raw_overlap = [cg; gg]
    return Matrix{Float64}(transpose(raw.raw_to_final) * raw_overlap)
end

"""
    gto_overlap_matrix(working, probes; block_indices = nothing)

Return the exact overlap matrix `S_BG = <B|G>` between a supported Cartesian
working basis `B` and a legacy Gaussian probe family `G`.

`working` may be a `CartesianBasisRepresentation3D`, a supported public
Cartesian basis object such as a bond-aligned ordinary basis or nested fixed
block, a `MappedUniformBasis` interpreted as its 3D Cartesian direct product, or
an `OrdinaryCartesianOperators3D` whose public representation carries the exact
hybrid Cartesian/supplement sidecars. `probes` may be a legacy Gaussian
supplement or a `CartesianGaussianShellSupplementRepresentation3D`.

When `block_indices` is supplied, only those working-basis rows are returned.
Block policy is intentionally external: this helper does not choose octants,
shells, DG regions, or any other geometry partition.
"""
function gto_overlap_matrix(working, probes; block_indices = nothing)
    working_representation = _gto_working_representation(working)
    probe_representation = _gto_probe_representation(probes)
    overlap = _gto_overlap_matrix(working_representation, probe_representation)
    indices = _gto_block_indices(working_representation.metadata.final_dimension, block_indices)
    indices === nothing && return overlap
    return Matrix{Float64}(overlap[indices, :])
end

function gto_overlap_matrix(
    working,
    probes,
    block_indices::AbstractVector{<:Integer},
)
    return gto_overlap_matrix(working, probes; block_indices = block_indices)
end

"""
    gto_occupancy_matrix(working, probes; weights = :uniform, block_indices = nothing)

Build a density-like GTO importance operator for basis-design work:

`M_I = S_{I,G} W S_{G,I}`

where `S_{I,G}` is the exact block-restricted overlap from
`gto_overlap_matrix(...)`. The result is an RDM-like design/occupancy matrix,
not a physical one-particle density matrix and not an SCF object.

First-pass weighting supports `weights = :uniform` and explicit real weight
vectors with one entry per probe orbital. `weights = :shell_equalized` is a
reserved future option and is intentionally not implemented yet.
"""
function gto_occupancy_matrix(working, probes; weights = :uniform, block_indices = nothing)
    overlap = gto_overlap_matrix(working, probes; block_indices = block_indices)
    weight_vector = _gto_weight_vector(weights, size(overlap, 2))
    weighted_overlap = overlap .* reshape(weight_vector, 1, :)
    return Matrix{Float64}(weighted_overlap * transpose(overlap))
end

function gto_occupancy_matrix(
    working,
    probes,
    block_indices::AbstractVector{<:Integer};
    weights = :uniform,
)
    return gto_occupancy_matrix(
        working,
        probes;
        weights = weights,
        block_indices = block_indices,
    )
end

function _cartesian_empty_centers()
    return zeros(Float64, 0, 3)
end

function _cartesian_supplement_center_matrix(
    supplement::CartesianGaussianShellSupplementRepresentation3D,
)
    isempty(supplement.orbitals) && return _cartesian_empty_centers()
    matrix = Matrix{Float64}(undef, length(supplement.orbitals), 3)
    for (row, orbital) in pairs(supplement.orbitals)
        matrix[row, 1] = orbital.center[1]
        matrix[row, 2] = orbital.center[2]
        matrix[row, 3] = orbital.center[3]
    end
    return matrix
end

function _cartesian_supplement_orbital_representation(
    orbital::_AtomicCartesianShellOrbital3D,
)
    return CartesianGaussianShellOrbitalRepresentation3D(
        String(orbital.label),
        (Int(orbital.lx), Int(orbital.ly), Int(orbital.lz)),
        (
            Float64(orbital.center[1]),
            Float64(orbital.center[2]),
            Float64(orbital.center[3]),
        ),
        Float64[Float64(value) for value in orbital.exponents],
        Float64[Float64(value) for value in orbital.coefficients],
        :axiswise_normalized_cartesian_gaussian,
    )
end

function _cartesian_supplement_kind(::Nothing)
    return :none
end

function _cartesian_supplement_kind(::LegacyAtomicGaussianSupplement)
    return :atomic_cartesian_shell
end

function _cartesian_supplement_kind(::LegacyBondAlignedDiatomicGaussianSupplement)
    return :bond_aligned_diatomic_cartesian_shell
end

function _cartesian_supplement_kind(::LegacyBondAlignedHeteronuclearGaussianSupplement)
    return :bond_aligned_heteronuclear_cartesian_shell
end

function _cartesian_supplement_lmax(::Nothing)
    return nothing
end

function _cartesian_supplement_lmax(data::LegacyAtomicGaussianSupplement)
    return data.lmax
end

function _cartesian_supplement_lmax(data::LegacyBondAlignedDiatomicGaussianSupplement)
    return data.atomic_source.lmax
end

function _cartesian_supplement_lmax(data::LegacyBondAlignedHeteronuclearGaussianSupplement)
    return maximum(source.lmax for source in data.atomic_sources)
end

function _cartesian_supplement_metadata(::Nothing)
    return (
        source_kind = :none,
        atom = nothing,
        basis_name = nothing,
        basisfile = nothing,
        lmax = nothing,
        nuclei = NTuple{3,Float64}[],
        uncontracted = nothing,
        max_width = nothing,
    )
end

function _cartesian_supplement_metadata(data::LegacyAtomicGaussianSupplement)
    return (
        source_kind = :legacy_atomic_gaussian_supplement,
        atom = String(data.atom),
        basis_name = String(data.basis_name),
        basisfile = String(data.basisfile),
        lmax = data.lmax,
        nuclei = NTuple{3,Float64}[(0.0, 0.0, 0.0)],
        uncontracted = data.uncontracted,
        max_width = data.max_width,
    )
end

function _cartesian_supplement_metadata(data::LegacyBondAlignedDiatomicGaussianSupplement)
    return (
        source_kind = :legacy_bond_aligned_diatomic_gaussian_supplement,
        atom = String(data.atomic_source.atom),
        basis_name = String(data.atomic_source.basis_name),
        basisfile = String(data.atomic_source.basisfile),
        lmax = data.atomic_source.lmax,
        nuclei = NTuple{3,Float64}[
            (
                Float64(nucleus[1]),
                Float64(nucleus[2]),
                Float64(nucleus[3]),
            ) for nucleus in data.nuclei
        ],
        uncontracted = data.atomic_source.uncontracted,
        max_width = data.max_width,
    )
end

function _cartesian_supplement_metadata(data::LegacyBondAlignedHeteronuclearGaussianSupplement)
    return (
        source_kind = :legacy_bond_aligned_heteronuclear_gaussian_supplement,
        atom = nothing,
        basis_name = nothing,
        basisfile = nothing,
        lmax = maximum(source.lmax for source in data.atomic_sources),
        nuclei = NTuple{3,Float64}[
            (
                Float64(nucleus[1]),
                Float64(nucleus[2]),
                Float64(nucleus[3]),
            ) for nucleus in data.nuclei
        ],
        uncontracted = all(source.uncontracted for source in data.atomic_sources),
        max_width = data.max_width,
    )
end

function _cartesian_empty_supplement_representation()
    return CartesianGaussianShellSupplementRepresentation3D(
        :none,
        CartesianGaussianShellOrbitalRepresentation3D[],
        _cartesian_supplement_metadata(nothing),
    )
end

function _cartesian_supplement_representation(::Nothing)
    return _cartesian_empty_supplement_representation()
end

function _cartesian_supplement_representation(
    data::Union{
        LegacyAtomicGaussianSupplement,
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    },
)
    supplement3d =
        data isa LegacyAtomicGaussianSupplement ? _atomic_cartesian_shell_supplement_3d(data) :
        _bond_aligned_diatomic_cartesian_shell_supplement_3d(data)
    orbitals = CartesianGaussianShellOrbitalRepresentation3D[
        _cartesian_supplement_orbital_representation(orbital) for orbital in supplement3d.orbitals
    ]
    return CartesianGaussianShellSupplementRepresentation3D(
        _cartesian_supplement_kind(data),
        orbitals,
        _cartesian_supplement_metadata(data),
    )
end

function basis_representation(
    data::Union{
        LegacyAtomicGaussianSupplement,
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    },
)
    return _cartesian_supplement_representation(data)
end

function basis_metadata(
    data::Union{
        LegacyAtomicGaussianSupplement,
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    },
)
    return basis_representation(data).metadata
end

function _cartesian_supplement_orbital_signature(
    orbital::CartesianGaussianShellOrbitalRepresentation3D,
)
    return (
        label = orbital.label,
        angular_powers = orbital.angular_powers,
        center = orbital.center,
        exponents = Tuple(Float64[Float64(value) for value in orbital.exponents]),
        coefficients = Tuple(Float64[Float64(value) for value in orbital.coefficients]),
        primitive_normalization = orbital.primitive_normalization,
    )
end

function _cartesian_same_supplement_raw_identity(
    left::CartesianGaussianShellSupplementRepresentation3D,
    right::CartesianGaussianShellSupplementRepresentation3D,
)
    left.supplement_kind == right.supplement_kind || return false
    isequal(left.metadata, right.metadata) || return false
    length(left.orbitals) == length(right.orbitals) || return false
    for index in eachindex(left.orbitals)
        _cartesian_supplement_orbital_signature(left.orbitals[index]) ==
        _cartesian_supplement_orbital_signature(right.orbitals[index]) || return false
    end
    return true
end
