# Route-owned source planning for multi-layer PQS shell realizations.

function _pqs_multilayer_axis_metrics(bundles::_CartesianNestedAxisBundles3D)
    pgdg_x = _nested_axis_pgdg(bundles, :x)
    pgdg_y = _nested_axis_pgdg(bundles, :y)
    pgdg_z = _nested_axis_pgdg(bundles, :z)
    return (;
        x = (overlap = pgdg_x.overlap, kinetic = pgdg_x.kinetic),
        y = (overlap = pgdg_y.overlap, kinetic = pgdg_y.kinetic),
        z = (overlap = pgdg_z.overlap, kinetic = pgdg_z.kinetic),
    )
end

function _pqs_multilayer_box_depth(
    outer_box::NTuple{3,UnitRange{Int}},
    core_box::NTuple{3,UnitRange{Int}},
)
    lower_depths = ntuple(axis -> first(core_box[axis]) - first(outer_box[axis]), 3)
    upper_depths = ntuple(axis -> last(outer_box[axis]) - last(core_box[axis]), 3)
    all(depth -> depth >= 1, lower_depths) ||
        throw(ArgumentError("multi-layer PQS shell plan requires the core box to be strictly inside the outer box"))
    lower_depths == upper_depths ||
        throw(ArgumentError("multi-layer PQS shell plan currently requires symmetric shell depth on every axis"))
    all(depth -> depth == lower_depths[1], lower_depths) ||
        throw(ArgumentError("multi-layer PQS shell plan currently requires the same shell depth on every axis"))
    return lower_depths[1]
end

function _pqs_multilayer_core_box_at_depth(
    core_box::NTuple{3,UnitRange{Int}},
    depth::Int,
)
    return ntuple(axis -> (first(core_box[axis]) - depth):(last(core_box[axis]) + depth), 3)
end

function _pqs_multilayer_block_concatenate_shell_coefficients(shell_records)
    total_support = sum(record -> record.support_count, shell_records)
    total_retained = sum(record -> record.retained_count, shell_records)
    coefficients = zeros(Float64, total_support, total_retained)
    support_offset = 0
    retained_offset = 0
    for record in shell_records
        support_range = (support_offset + 1):(support_offset + record.support_count)
        retained_range = (retained_offset + 1):(retained_offset + record.retained_count)
        coefficients[support_range, retained_range] .= record.shell_final_coefficients
        support_offset += record.support_count
        retained_offset += record.retained_count
    end
    return coefficients
end

function _pqs_multilayer_duplicate_count(values)
    return length(values) - length(unique(values))
end

function _pqs_multilayer_property(plan, key::Symbol, default = nothing)
    return hasproperty(plan, key) ? getproperty(plan, key) : default
end

function _pqs_multilayer_support_product_matrix(
    left_states::AbstractVector{<:NTuple{3,Int}},
    right_states::AbstractVector{<:NTuple{3,Int}},
    operator_x::AbstractMatrix{<:Real},
    operator_y::AbstractMatrix{<:Real},
    operator_z::AbstractMatrix{<:Real},
)
    result = Matrix{Float64}(undef, length(left_states), length(right_states))
    @inbounds for (left_index, (ix, iy, iz)) in pairs(left_states)
        for (right_index, (jx, jy, jz)) in pairs(right_states)
            result[left_index, right_index] =
                Float64(operator_x[ix, jx]) *
                Float64(operator_y[iy, jy]) *
                Float64(operator_z[iz, jz])
        end
    end
    return result
end

function _pqs_multilayer_axis_tuple(axis_data, name::AbstractString)
    if axis_data isa NamedTuple
        all(key -> haskey(axis_data, key), (:x, :y, :z)) ||
            throw(ArgumentError("$name must provide x, y, and z entries"))
        return (axis_data.x, axis_data.y, axis_data.z)
    end
    length(axis_data) == 3 ||
        throw(ArgumentError("$name must provide exactly three axis entries"))
    return (axis_data[1], axis_data[2], axis_data[3])
end

function _pqs_multilayer_center_records(center_records)
    isnothing(center_records) &&
        throw(ArgumentError("PQS multi-layer support electron-nuclear requires center records"))
    center_records isa NamedTuple && return [center_records]
    return collect(center_records)
end

function _pqs_multilayer_center_property(center_record, key::Symbol, default = nothing)
    center_record isa NamedTuple && haskey(center_record, key) &&
        return getfield(center_record, key)
    hasproperty(center_record, key) && return getproperty(center_record, key)
    return default
end

function _pqs_multilayer_center_summary(center_record)
    charge = _pqs_multilayer_center_property(center_record, :charge)
    isnothing(charge) &&
        (charge = _pqs_multilayer_center_property(center_record, :nuclear_charge))
    location = _pqs_multilayer_center_property(center_record, :location)
    isnothing(charge) &&
        throw(ArgumentError("PQS multi-layer support electron-nuclear requires recorded nuclear charge"))
    isnothing(location) &&
        throw(ArgumentError("PQS multi-layer support electron-nuclear requires center location"))
    length(location) == 3 ||
        throw(ArgumentError("PQS multi-layer support electron-nuclear requires a 3D center location"))
    return (;
        center_key = _pqs_multilayer_center_property(center_record, :center_key, :unknown),
        center_index = _pqs_multilayer_center_property(center_record, :center_index, :unknown),
        location = (Float64(location[1]), Float64(location[2]), Float64(location[3])),
        nuclear_charge = Float64(charge),
    )
end

function _pqs_multilayer_explicit_factor_terms(
    gaussian_factor_terms_by_center,
    center_index::Int,
    center_count::Int,
)
    if gaussian_factor_terms_by_center isa AbstractArray
        center_count == 1 ||
            throw(ArgumentError("single explicit Gaussian factor array is only valid for one center"))
        return (
            gaussian_factor_terms_by_center,
            gaussian_factor_terms_by_center,
            gaussian_factor_terms_by_center,
            :explicit_single_axis_factor_reused_for_xyz,
        )
    end
    if gaussian_factor_terms_by_center isa NamedTuple &&
       all(key -> haskey(gaussian_factor_terms_by_center, key), (:x, :y, :z))
        center_count == 1 ||
            throw(ArgumentError("single explicit x/y/z Gaussian factor set is only valid for one center"))
        return (
            gaussian_factor_terms_by_center.x,
            gaussian_factor_terms_by_center.y,
            gaussian_factor_terms_by_center.z,
            :explicit_xyz_axis_factors,
        )
    end
    factor_sets = collect(gaussian_factor_terms_by_center)
    length(factor_sets) == center_count ||
        throw(ArgumentError("explicit Gaussian factor terms must match center record count"))
    factors = factor_sets[center_index]
    if factors isa AbstractArray
        return (
            factors,
            factors,
            factors,
            :explicit_single_axis_factor_reused_for_xyz_by_center,
        )
    end
    if factors isa NamedTuple && all(key -> haskey(factors, key), (:x, :y, :z))
        return (factors.x, factors.y, factors.z, :explicit_xyz_axis_factors_by_center)
    end
    axes = _pqs_multilayer_axis_tuple(factors, "explicit Gaussian factor terms")
    return (axes[1], axes[2], axes[3], :explicit_xyz_axis_factors_by_center)
end

function _pqs_multilayer_centered_factor_terms(axis_layers, coulomb_expansion, center)
    isnothing(axis_layers) &&
        throw(ArgumentError("off-origin support electron-nuclear requires axis_layers or explicit Gaussian factor terms"))
    exponents = Float64.(coulomb_expansion.exponents)
    layers = _pqs_multilayer_axis_tuple(axis_layers, "axis_layers")
    factors = ntuple(
        axis -> gaussian_factor_matrices(
            layers[axis];
            exponents,
            center = center.location[axis],
        ),
        3,
    )
    return (
        _pqs_multilayer_term_first_factor_array(factors[1]),
        _pqs_multilayer_term_first_factor_array(factors[2]),
        _pqs_multilayer_term_first_factor_array(factors[3]),
        :centered_axis_layers,
    )
end

function _pqs_multilayer_term_first_factor_array(factors::AbstractArray{<:Real,3})
    return Array{Float64,3}(factors)
end

function _pqs_multilayer_term_first_factor_array(factors)
    matrices = collect(factors)
    isempty(matrices) &&
        throw(ArgumentError("PQS multi-layer centered Gaussian factors cannot be empty"))
    first_matrix = Matrix{Float64}(first(matrices))
    result = Array{Float64,3}(
        undef,
        length(matrices),
        size(first_matrix, 1),
        size(first_matrix, 2),
    )
    result[1, :, :] .= first_matrix
    for term_index in 2:length(matrices)
        matrix = Matrix{Float64}(matrices[term_index])
        size(matrix) == size(first_matrix) ||
            throw(ArgumentError("PQS multi-layer centered Gaussian factor matrices must share one shape"))
        result[term_index, :, :] .= matrix
    end
    return result
end

function _pqs_multilayer_validate_factor_terms(axis_terms, term_count::Int)
    for axis in 1:3
        ndims(axis_terms[axis]) == 3 ||
            throw(ArgumentError("PQS multi-layer Gaussian factor terms must be term-first 3D arrays"))
        size(axis_terms[axis], 1) == term_count ||
            throw(ArgumentError("PQS multi-layer Gaussian factor term count mismatch"))
    end
    return axis_terms
end

function _pqs_multilayer_support_electron_nuclear_matrix(
    states,
    axis_terms,
    coefficients,
)
    result = zeros(Float64, length(states), length(states))
    for term_index in eachindex(coefficients)
        result .+=
            -Float64(coefficients[term_index]) *
            _pqs_multilayer_support_product_matrix(
                states,
                states,
                @view(axis_terms[1][term_index, :, :]),
                @view(axis_terms[2][term_index, :, :]),
                @view(axis_terms[3][term_index, :, :]),
            )
    end
    return result
end

"""
    pqs_multilayer_shell_source_plan(bundles, core_box, outer_box; ...)

Build a compact route-owned source plan for a direct PQS core plus repeated
one-cell surrounding shell layers. Each layer is materialized with the existing
projected-q-shell descriptor and shell-realization plan, then the disjoint shell
supports and shell isometries are collapsed into one shell sector suitable for
`CartesianFinalBasisRealization.pqs_complete_core_shell_final_basis`.

This helper plans shell/source data only. It does not build final-basis overlap,
one-body operators, H1, IDA, RHF, driver wiring, exports, or artifacts.
"""
function pqs_multilayer_shell_source_plan(
    bundles::_CartesianNestedAxisBundles3D,
    core_box::NTuple{3,UnitRange{Int}},
    outer_box::NTuple{3,UnitRange{Int}};
    bond_axis::Symbol = :z,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
    metadata = (;),
)
    dims = _nested_axis_lengths(bundles)
    for axis in 1:3
        first(outer_box[axis]) >= 1 && last(outer_box[axis]) <= dims[axis] ||
            throw(ArgumentError("multi-layer PQS outer box must lie inside parent dimensions"))
    end
    layer_count = _pqs_multilayer_box_depth(outer_box, core_box)
    metrics = _pqs_multilayer_axis_metrics(bundles)
    core_support_indices =
        _nested_box_support_indices(core_box[1], core_box[2], core_box[3], dims)
    core_support_states =
        [_cartesian_unflat_index(index, dims) for index in core_support_indices]

    shell_records = Vector{NamedTuple}(undef, layer_count)
    for layer_index in 1:layer_count
        inner_box =
            layer_index == 1 ?
            core_box :
            _pqs_multilayer_core_box_at_depth(core_box, layer_index - 1)
        current_box = _pqs_multilayer_core_box_at_depth(core_box, layer_index)
        q_values = length.(inner_box)
        all(q -> q == q_values[1], q_values) ||
            throw(ArgumentError("multi-layer PQS shell layers currently require cubic raw source dimensions"))
        q = q_values[1]
        layer = _nested_projected_q_shell_layer(
            bundles,
            current_box,
            inner_box;
            bond_axis,
            q,
            L = q,
            raw_source_dims = (q, q, q),
            selected_q = q,
            term_coefficients,
        )
        descriptor = _nested_projected_q_shell_staged_unit_descriptor(layer)
        shell_plan =
            CartesianContractedParentMetrics._pqs_shell_realization_plan(
                descriptor,
                metrics,
            )
        shell_final_coefficients =
            shell_plan.shell_projection_matrix * shell_plan.lowdin_cleanup
        shell_records[layer_index] = (;
            layer_index,
            current_box,
            inner_box,
            raw_source_dims = (q, q, q),
            q,
            descriptor,
            shell_plan,
            shell_support_indices = descriptor.support_indices,
            shell_support_states = descriptor.support_states,
            shell_final_coefficients,
            support_count = descriptor.support_count,
            mode_count = descriptor.mode_count,
            retained_count = descriptor.retained_count,
            isometry_error = shell_plan.isometry_error,
        )
    end

    shell_support_indices =
        reduce(vcat, (record.shell_support_indices for record in shell_records); init = Int[])
    shell_support_states =
        reduce(vcat, (record.shell_support_states for record in shell_records); init = Tuple{Int,Int,Int}[])
    shell_final_coefficients =
        _pqs_multilayer_block_concatenate_shell_coefficients(shell_records)
    core_shell_duplicate_count = length(intersect(core_support_indices, shell_support_indices))
    shell_duplicate_count = _pqs_multilayer_duplicate_count(shell_support_indices)
    combined_support_indices = vcat(core_support_indices, shell_support_indices)
    intended_support_indices =
        _nested_box_support_indices(outer_box[1], outer_box[2], outer_box[3], dims)
    combined_unique = sort!(unique(combined_support_indices))
    intended_sorted = sort!(collect(intended_support_indices))
    covers_outer_box = combined_unique == intended_sorted
    blocker =
        shell_duplicate_count != 0 ? :duplicate_shell_support :
        core_shell_duplicate_count != 0 ? :core_shell_support_overlap :
        !covers_outer_box ? :combined_support_does_not_cover_outer_box :
        nothing
    status = isnothing(blocker) ?
        :available_pqs_multilayer_shell_source_plan :
        :blocked_pqs_multilayer_shell_source_plan

    summary = (;
        status,
        blocker,
        layer_count,
        core_support_count = length(core_support_indices),
        shell_support_count = length(shell_support_indices),
        shell_final_retained_count = size(shell_final_coefficients, 2),
        combined_support_count = length(combined_support_indices),
        intended_support_count = length(intended_sorted),
        shell_duplicate_count,
        core_shell_duplicate_count,
        covers_outer_box,
        final_basis_helper = :pqs_complete_core_shell_final_basis,
        collapsed_shell_sector = true,
        fixed_block_matrix_authority_used = false,
        h1_materialized = false,
        rhf_materialized = false,
    )

    return (;
        object_kind = :pqs_multilayer_shell_source_plan,
        status,
        blocker,
        source_kind = :repeated_one_cell_projected_q_shell_layers,
        bundles,
        metrics,
        core_box,
        outer_box,
        bond_axis,
        layer_count,
        shell_records,
        core_support_indices,
        core_support_states,
        shell_support_indices,
        shell_support_states,
        shell_final_coefficients,
        summary,
        metadata = merge(
            NamedTuple(metadata),
            (;
                source = :pqs_multilayer_shell_source_plan,
                collapsed_shell_sector = true,
                old_fixed_block_matrix_authority_used = false,
            ),
        ),
    )
end

"""
    pqs_multilayer_complete_core_shell_final_basis(plan; metadata = (;), ...)

Build the complete core/shell final-basis realization for a route-owned
multi-layer PQS shell source plan. This helper performs only the overlap
assembly needed by `CartesianFinalBasisRealization`; it does not materialize
one-body operators, H1, IDA, density-density, RHF, driver wiring, exports, or
artifacts.
"""
function pqs_multilayer_complete_core_shell_final_basis(
    plan;
    identity_atol::Real = 1.0e-8,
    rank_atol::Real = 1.0e-10,
    metadata = (;),
)
    _pqs_multilayer_property(plan, :object_kind) ===
        :pqs_multilayer_shell_source_plan ||
        throw(ArgumentError("PQS multi-layer final basis requires a pqs_multilayer_shell_source_plan"))
    if plan.status !== :available_pqs_multilayer_shell_source_plan
        return _blocked_pqs_multilayer_complete_core_shell_final_basis(
            plan;
            blocker = isnothing(plan.blocker) ?
                :pqs_multilayer_shell_source_plan_not_available :
                plan.blocker,
            metadata,
        )
    end

    metrics = plan.metrics
    core_overlap = _pqs_multilayer_support_product_matrix(
        plan.core_support_states,
        plan.core_support_states,
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    core_shell_overlap = _pqs_multilayer_support_product_matrix(
        plan.core_support_states,
        plan.shell_support_states,
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    shell_overlap = _pqs_multilayer_support_product_matrix(
        plan.shell_support_states,
        plan.shell_support_states,
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    merged_metadata = merge(
        NamedTuple(plan.metadata),
        NamedTuple(metadata),
        (;
            source = :pqs_multilayer_complete_core_shell_final_basis,
            input_source_plan = :pqs_multilayer_shell_source_plan,
            overlap_blocks_built_from_plan_support_states = true,
            h1_materialized = false,
            ida_data_materialized = false,
            density_density_materialized = false,
            rhf_materialized = false,
            driver_route_materialized = false,
            exports_materialized = false,
            artifacts_materialized = false,
        ),
    )
    final_basis =
        CartesianFinalBasisRealization.pqs_complete_core_shell_final_basis(
            core_support_indices = plan.core_support_indices,
            shell_support_indices = plan.shell_support_indices,
            core_overlap = core_overlap,
            core_shell_overlap = core_shell_overlap,
            shell_overlap = shell_overlap,
            shell_final_coefficients = plan.shell_final_coefficients,
            identity_atol = identity_atol,
            rank_atol = rank_atol,
            metadata = merged_metadata,
        )
    return merge(
        final_basis,
        (;
            source_plan_object_kind = plan.object_kind,
            source_plan_status = plan.status,
            source_plan_layer_count = plan.layer_count,
            source_plan_core_support_count = length(plan.core_support_indices),
            source_plan_shell_support_count = length(plan.shell_support_indices),
            source_plan_final_basis_helper =
                :pqs_multilayer_complete_core_shell_final_basis,
            h1_materialized = false,
            ida_data_materialized = false,
            density_density_materialized = false,
            rhf_materialized = false,
            driver_route_materialized = false,
            exports_materialized = false,
            artifacts_materialized = false,
        ),
    )
end

"""
    pqs_multilayer_support_kinetic_matrix(plan)

Build the support-space Cartesian kinetic matrix for a route-owned multi-layer
PQS shell source plan. The matrix is ordered over the direct core support rows
followed by the collapsed shell support rows. This helper does not perform
final-basis transfer, H1 assembly, nuclear assembly, IDA, RHF, driver wiring,
exports, or artifacts.
"""
function pqs_multilayer_support_kinetic_matrix(plan)
    _pqs_multilayer_property(plan, :object_kind) ===
        :pqs_multilayer_shell_source_plan ||
        throw(ArgumentError("PQS multi-layer support kinetic requires a pqs_multilayer_shell_source_plan"))
    plan.status === :available_pqs_multilayer_shell_source_plan ||
        throw(ArgumentError("PQS multi-layer support kinetic requires an available source plan"))

    states = vcat(plan.core_support_states, plan.shell_support_states)
    metrics = plan.metrics
    return (
        _pqs_multilayer_support_product_matrix(
            states,
            states,
            metrics.x.kinetic,
            metrics.y.overlap,
            metrics.z.overlap,
        ) +
        _pqs_multilayer_support_product_matrix(
            states,
            states,
            metrics.x.overlap,
            metrics.y.kinetic,
            metrics.z.overlap,
        ) +
        _pqs_multilayer_support_product_matrix(
            states,
            states,
            metrics.x.overlap,
            metrics.y.overlap,
            metrics.z.kinetic,
        )
    )
end

"""
    pqs_multilayer_support_electron_nuclear_by_center_matrices(plan; ...)

Build separated uncharged support-space electron-nuclear by-center matrices for
a route-owned multi-layer PQS shell source plan. Each center matrix represents
negative unit-charge attraction, `-1/r_center`; nuclear charges are recorded
but not applied, and centers are not summed. This helper does not perform
final-basis transfer, H1 assembly, IDA, density-density, RHF, driver wiring,
exports, or artifacts.
"""
function pqs_multilayer_support_electron_nuclear_by_center_matrices(
    plan;
    coulomb_expansion,
    center_records,
    axis_layers = nothing,
    gaussian_factor_terms_by_center = nothing,
)
    _pqs_multilayer_property(plan, :object_kind) ===
        :pqs_multilayer_shell_source_plan ||
        throw(ArgumentError("PQS multi-layer support electron-nuclear requires a pqs_multilayer_shell_source_plan"))
    plan.status === :available_pqs_multilayer_shell_source_plan ||
        throw(ArgumentError("PQS multi-layer support electron-nuclear requires an available source plan"))

    centers = _pqs_multilayer_center_records(center_records)
    isempty(centers) &&
        throw(ArgumentError("PQS multi-layer support electron-nuclear requires at least one center"))
    coefficients = Float64.(coulomb_expansion.coefficients)
    states = vcat(plan.core_support_states, plan.shell_support_states)
    records = map(enumerate(centers)) do (center_index, center_record)
        center = _pqs_multilayer_center_summary(center_record)
        factor_x, factor_y, factor_z, factor_source =
            isnothing(gaussian_factor_terms_by_center) ?
            _pqs_multilayer_centered_factor_terms(axis_layers, coulomb_expansion, center) :
            _pqs_multilayer_explicit_factor_terms(
                gaussian_factor_terms_by_center,
                center_index,
                length(centers),
            )
        axis_terms =
            _pqs_multilayer_validate_factor_terms(
                (factor_x, factor_y, factor_z),
                length(coefficients),
            )
        support_operator =
            _pqs_multilayer_support_electron_nuclear_matrix(
                states,
                axis_terms,
                coefficients,
            )
        (;
            object_kind = :pqs_multilayer_support_electron_nuclear_by_center_matrix,
            status = :materialized_pqs_multilayer_support_electron_nuclear_by_center_matrix,
            blocker = nothing,
            term = :electron_nuclear_by_center,
            support_operator,
            support_operator_shape = size(support_operator),
            support_operator_finite = all(isfinite, support_operator),
            support_state_count = length(states),
            support_ordering = :core_support_states_then_shell_support_states,
            by_center = true,
            center_key = center.center_key,
            center_index = center.center_index,
            center_location = center.location,
            nuclear_charge = center.nuclear_charge,
            nuclear_charge_recorded = true,
            nuclear_charge_applied = false,
            centers_summed = false,
            center_summation = false,
            uncharged_by_center_convention = true,
            gaussian_factor_terms_source = factor_source,
            gaussian_term_count = length(coefficients),
            final_basis_transfer_materialized = false,
            h1_materialized = false,
            ida_data_materialized = false,
            density_density_materialized = false,
            rhf_materialized = false,
            driver_route_materialized = false,
            exports_materialized = false,
            artifacts_materialized = false,
            metadata = (;
                source = :pqs_multilayer_support_electron_nuclear_by_center_matrices,
                physical_operator = :electron_nuclear_attraction,
                negative_unit_charge_attraction = true,
                nuclear_charge_recorded = true,
                nuclear_charge_applied = false,
                centers_summed = false,
                uncharged_by_center_convention = true,
                charge_application_stage = :hamiltonian_assembly,
                gaussian_factor_terms_source = factor_source,
                old_fixed_block_matrix_authority_used = false,
                wl_matrix_authority_used = false,
            ),
        )
    end
    return (;
        object_kind = :pqs_multilayer_support_electron_nuclear_by_center_matrix_set,
        status = :materialized_pqs_multilayer_support_electron_nuclear_by_center_matrix_set,
        blocker = nothing,
        term = :electron_nuclear_by_center,
        records,
        center_count = length(records),
        support_state_count = length(states),
        support_ordering = :core_support_states_then_shell_support_states,
        by_center = true,
        nuclear_charge_recorded = all(record -> record.nuclear_charge_recorded, records),
        nuclear_charge_applied = false,
        centers_summed = false,
        uncharged_by_center_convention = true,
        final_basis_transfer_materialized = false,
        h1_materialized = false,
        ida_data_materialized = false,
        density_density_materialized = false,
        rhf_materialized = false,
        driver_route_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
    )
end

function _blocked_pqs_multilayer_complete_core_shell_final_basis(
    plan;
    blocker,
    metadata = (;),
)
    return (;
        object_kind = :pqs_multilayer_complete_core_shell_final_basis,
        status = :blocked_pqs_multilayer_complete_core_shell_final_basis,
        blocker,
        source_plan_object_kind = _pqs_multilayer_property(plan, :object_kind),
        source_plan_status = _pqs_multilayer_property(plan, :status),
        source_plan_blocker = _pqs_multilayer_property(plan, :blocker),
        source_plan_layer_count = _pqs_multilayer_property(plan, :layer_count, 0),
        source_plan_core_support_count =
            hasproperty(plan, :core_support_indices) ?
            length(plan.core_support_indices) : 0,
        source_plan_shell_support_count =
            hasproperty(plan, :shell_support_indices) ?
            length(plan.shell_support_indices) : 0,
        final_basis_materialized = false,
        one_body_operator_materialized = false,
        h1_materialized = false,
        ida_data_materialized = false,
        density_density_materialized = false,
        rhf_materialized = false,
        driver_route_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
        metadata = merge(
            hasproperty(plan, :metadata) ? NamedTuple(plan.metadata) : (;),
            NamedTuple(metadata),
            (;
                source = :pqs_multilayer_complete_core_shell_final_basis,
                input_source_plan = :pqs_multilayer_shell_source_plan,
                final_basis_helper = :pqs_complete_core_shell_final_basis,
            ),
        ),
    )
end
