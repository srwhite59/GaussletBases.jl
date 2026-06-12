# Dense support-space one-body helpers for the multi-layer PQS H1 seam.

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
