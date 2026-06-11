# By-center electron-nuclear White--Lindsey pair-block pilot.
#
# This path is decomposed-boundary-unit local: it builds one center's support
# block from parent-axis Gaussian factor terms and contracts it with prepared
# White--Lindsey unit coefficient maps. It does not build a full parent-product
# matrix, sum centers, apply physical nuclear charges, assemble Hamiltonians,
# export artifacts, or use IDA/MWG/PQS semantics.

function white_lindsey_boundary_stratum_electron_nuclear_by_center_block(
    pair_unit_coefficients;
    parent_axis_counts,
    parent_axis_bundle_object,
    coulomb_expansion,
    center_record,
    electron_nuclear_axis_context = nothing,
)
    _assert_white_lindsey_pair_unit_coefficients_ready(pair_unit_coefficients)
    axis_context = _white_lindsey_electron_nuclear_axis_context(
        electron_nuclear_axis_context;
        parent_axis_counts,
        parent_axis_bundle_object,
        coulomb_expansion,
        center_record,
    )
    axis_counts = axis_context.axis_counts
    center_summary = axis_context.center_summary
    coefficients = axis_context.coefficients
    axis_terms = axis_context.axis_terms

    left_support = pair_unit_coefficients.left_support_indices
    right_support = pair_unit_coefficients.right_support_indices
    _assert_white_lindsey_overlap_support(left_support, :left)
    _assert_white_lindsey_overlap_support(right_support, :right)

    support_block = _white_lindsey_electron_nuclear_support_block(
        left_support,
        right_support,
        axis_counts,
        axis_terms,
        coefficients,
    )
    left_support_coefficients = _white_lindsey_support_coefficient_matrix(
        pair_unit_coefficients.left_coefficient_matrix,
        left_support,
        pair_unit_coefficients.left_coefficient_space,
        :left,
    )
    right_support_coefficients = _white_lindsey_support_coefficient_matrix(
        pair_unit_coefficients.right_coefficient_matrix,
        right_support,
        pair_unit_coefficients.right_coefficient_space,
        :right,
    )
    block = Matrix{Float64}(
        transpose(left_support_coefficients) *
        support_block *
        right_support_coefficients,
    )

    return PairBlockMaterializationResult(
        :electron_nuclear_by_center,
        pair_unit_coefficients.pair_key,
        block,
        true,
        true,
        true,
        false,
        false,
        false,
        (;
            materialization_path =
                :white_lindsey_boundary_stratum_electron_nuclear_by_center_adapter,
            pair_family = pair_unit_coefficients.pair_family,
            left_stratum_kind = pair_unit_coefficients.left_stratum_kind,
            right_stratum_kind = pair_unit_coefficients.right_stratum_kind,
            left_support_count = length(left_support),
            right_support_count = length(right_support),
            left_retained_column_count =
                pair_unit_coefficients.left_retained_column_count,
            right_retained_column_count =
                pair_unit_coefficients.right_retained_column_count,
            parent_axis_counts = axis_counts,
            support_block_shape = size(support_block),
            left_support_coefficient_shape =
                size(left_support_coefficients),
            right_support_coefficient_shape =
                size(right_support_coefficients),
            coefficient_source =
                :white_lindsey_boundary_stratum_pair_unit_coefficients,
            center_key = center_summary.center_key,
            center_index = center_summary.center_index,
            center_location = center_summary.location,
            nuclear_charge = center_summary.charge,
            nuclear_charge_recorded = true,
            nuclear_charge_applied = false,
            charge_application_stage =
                :acceptance_or_hamiltonian_assembly,
            by_center = true,
            centers_summed = false,
            uncharged_by_center_convention = true,
            factor_source_path =
                :axis_pgdg_intermediate_gaussian_factor_terms,
            axis_term_cache_scope = axis_context.cache_scope,
            axis_term_cache_status = axis_context.cache_status,
            gaussian_expansion_loop = :inner_support_contraction,
            gaussian_term_count = length(coefficients),
            axis_factor_term_shapes =
                _white_lindsey_electron_nuclear_axis_term_shapes(axis_terms),
            local_pair_block_materialized = true,
            source_operator_blocks_materialized = true,
            final_pair_blocks_materialized = true,
            operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
            artifacts_materialized = false,
            ida_mwg_data_materialized = false,
            dense_parent_product_matrix_materialized = false,
            full_parent_window_cpb_used = false,
        ),
    )
end

function white_lindsey_boundary_stratum_electron_nuclear_by_center_block(
    unit_pair::CUP.UnitPairRecord;
    parent_axis_counts,
    parent_axis_bundle_object,
    coulomb_expansion,
    center_record,
    electron_nuclear_axis_context = nothing,
)
    return white_lindsey_boundary_stratum_electron_nuclear_by_center_block(
        white_lindsey_boundary_stratum_pair_unit_coefficients(unit_pair);
        parent_axis_counts,
        parent_axis_bundle_object,
        coulomb_expansion,
        center_record,
        electron_nuclear_axis_context,
    )
end

function _white_lindsey_electron_nuclear_axis_context(
    context;
    parent_axis_counts,
    parent_axis_bundle_object,
    coulomb_expansion,
    center_record,
)
    if !isnothing(context)
        _assert_white_lindsey_electron_nuclear_axis_context(context)
        return context
    end
    return _white_lindsey_electron_nuclear_axis_context(
        parent_axis_counts,
        parent_axis_bundle_object,
        coulomb_expansion,
        center_record,
    )
end

function _white_lindsey_electron_nuclear_axis_context(
    parent_axis_counts,
    parent_axis_bundle_object,
    coulomb_expansion,
    center_record,
)
    axis_counts = _axis_counts_tuple(parent_axis_counts)
    center_summary = _white_lindsey_electron_nuclear_center_summary(center_record)
    center_summary.status === :available_electron_nuclear_center_record ||
        throw(
            ArgumentError(
                "White--Lindsey electron-nuclear by-center block requires an available center record",
            ),
        )
    coefficients =
        _white_lindsey_electron_nuclear_expansion_coefficients(coulomb_expansion)
    exponents =
        _white_lindsey_electron_nuclear_expansion_exponents(coulomb_expansion)
    axis_terms = _white_lindsey_electron_nuclear_axis_terms(
        parent_axis_bundle_object,
        exponents,
        center_summary,
    )
    _assert_white_lindsey_electron_nuclear_axis_terms(
        axis_terms,
        coefficients,
        axis_counts,
    )
    return (;
        object_kind = :white_lindsey_electron_nuclear_axis_context,
        status = :available_white_lindsey_electron_nuclear_axis_context,
        axis_counts,
        center_summary,
        coefficients,
        exponents,
        axis_terms,
        cache_scope = :center_axis_expansion_parent_context,
        cache_status = :centered_axis_terms_reused_across_unit_pairs,
    )
end

function _assert_white_lindsey_electron_nuclear_axis_context(context)
    _white_lindsey_descriptor_property(context, :object_kind) ===
    :white_lindsey_electron_nuclear_axis_context || throw(
        ArgumentError(
            "White--Lindsey electron-nuclear block requires an axis context built by the electron-nuclear axis context helper",
        ),
    )
    _white_lindsey_descriptor_property(context, :status) ===
    :available_white_lindsey_electron_nuclear_axis_context || throw(
        ArgumentError(
            "White--Lindsey electron-nuclear block requires an available axis context",
        ),
    )
    return nothing
end

function _white_lindsey_electron_nuclear_center_summary(center_record)
    isnothing(center_record) && return (;
        status = :blocked_electron_nuclear_center_record,
        blocker = :missing_electron_nuclear_center_record,
        center_key = :unavailable,
        center_index = :unavailable,
        charge = :unavailable,
        location = :unavailable,
        x = :unavailable,
        y = :unavailable,
        z = :unavailable,
    )
    charge = _white_lindsey_descriptor_property(center_record, :charge)
    isnothing(charge) &&
        (charge = _white_lindsey_descriptor_property(center_record, :nuclear_charge))
    isnothing(charge) &&
        (charge = _white_lindsey_descriptor_property(center_record, :Z))
    location = _white_lindsey_descriptor_property(center_record, :location)
    if isnothing(location)
        x = _white_lindsey_descriptor_property(center_record, :x)
        y = _white_lindsey_descriptor_property(center_record, :y)
        z = _white_lindsey_descriptor_property(center_record, :z)
        if !isnothing(x) && !isnothing(y) && !isnothing(z)
            location = (x, y, z)
        end
    end
    location_tuple =
        _white_lindsey_electron_nuclear_location_tuple(location)
    blocker =
        isnothing(charge) ?
        :missing_electron_nuclear_center_charge :
        (
            isnothing(location_tuple) ?
            :missing_electron_nuclear_center_location :
            nothing
        )
    return (;
        status =
            isnothing(blocker) ?
            :available_electron_nuclear_center_record :
            :blocked_electron_nuclear_center_record,
        blocker,
        center_key = something(
            _white_lindsey_descriptor_property(center_record, :center_key),
            :unavailable,
        ),
        center_index = something(
            _white_lindsey_descriptor_property(center_record, :center_index),
            :unavailable,
        ),
        charge = isnothing(charge) ? :unavailable : Float64(charge),
        location = isnothing(location_tuple) ? :unavailable : location_tuple,
        x = isnothing(location_tuple) ? :unavailable : location_tuple[1],
        y = isnothing(location_tuple) ? :unavailable : location_tuple[2],
        z = isnothing(location_tuple) ? :unavailable : location_tuple[3],
    )
end

function _white_lindsey_electron_nuclear_location_tuple(location)
    isnothing(location) && return nothing
    length(location) == 3 || return nothing
    return (Float64(location[1]), Float64(location[2]), Float64(location[3]))
end

function _white_lindsey_electron_nuclear_expansion_coefficients(expansion)
    coefficients = _white_lindsey_descriptor_property(expansion, :coefficients)
    coefficients isa AbstractVector ||
        throw(
            ArgumentError(
                "White--Lindsey electron-nuclear by-center block requires expansion coefficients",
            ),
        )
    return Float64[-Float64(value) for value in coefficients]
end

function _white_lindsey_electron_nuclear_expansion_exponents(expansion)
    exponents = _white_lindsey_descriptor_property(expansion, :exponents)
    exponents isa AbstractVector ||
        throw(
            ArgumentError(
                "White--Lindsey electron-nuclear by-center block requires expansion exponents",
            ),
        )
    return exponents
end

function _white_lindsey_electron_nuclear_axis_terms(
    parent_axis_bundle_object,
    exponents,
    center_summary,
)
    return (;
        x = _white_lindsey_electron_nuclear_axis_terms(
            parent_axis_bundle_object,
            :x,
            exponents,
            center_summary.x,
        ),
        y = _white_lindsey_electron_nuclear_axis_terms(
            parent_axis_bundle_object,
            :y,
            exponents,
            center_summary.y,
        ),
        z = _white_lindsey_electron_nuclear_axis_terms(
            parent_axis_bundle_object,
            :z,
            exponents,
            center_summary.z,
        ),
    )
end

function _white_lindsey_electron_nuclear_axis_terms(
    parent_axis_bundle_object,
    axis::Symbol,
    exponents,
    center_coordinate,
)
    axis_bundle = _white_lindsey_axis_bundle(parent_axis_bundle_object, axis)
    basis = _white_lindsey_descriptor_property(axis_bundle, :basis)
    backend = _white_lindsey_descriptor_property(axis_bundle, :backend)
    isnothing(basis) &&
        throw(ArgumentError("missing $(axis) axis basis for electron-nuclear block"))
    isnothing(backend) &&
        throw(ArgumentError("missing $(axis) axis backend for electron-nuclear block"))
    centered_bundle = ParentGaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents,
        center = center_coordinate,
        backend,
    )
    pgdg_intermediate =
        _white_lindsey_descriptor_property(centered_bundle, :pgdg_intermediate)
    terms =
        _white_lindsey_descriptor_property(
            pgdg_intermediate,
            :gaussian_factor_terms,
        )
    terms isa AbstractArray{<:Real,3} ||
        throw(
            ArgumentError(
                "missing $(axis) electron-nuclear Gaussian factor terms",
            ),
        )
    return terms
end

function _white_lindsey_axis_bundle(parent_axis_bundle_object, axis::Symbol)
    value = _white_lindsey_descriptor_property(parent_axis_bundle_object, axis)
    isnothing(value) &&
        throw(ArgumentError("missing $(axis) axis bundle for electron-nuclear block"))
    return value
end

function _assert_white_lindsey_electron_nuclear_axis_terms(
    axis_terms,
    coefficients,
    axis_counts::NTuple{3,Int},
)
    for (axis_index, axis) in enumerate((:x, :y, :z))
        terms = getproperty(axis_terms, axis)
        size(terms, 1) == length(coefficients) ||
            throw(
                ArgumentError(
                    "electron-nuclear Gaussian term count mismatch on $(axis)",
                ),
            )
        size(terms, 2) == axis_counts[axis_index] &&
            size(terms, 3) == axis_counts[axis_index] ||
            throw(
                ArgumentError(
                    "electron-nuclear axis factor shape mismatch on $(axis)",
                ),
            )
    end
    return nothing
end

function _white_lindsey_electron_nuclear_support_block(
    left_support,
    right_support,
    axis_counts::NTuple{3,Int},
    axis_terms,
    coefficients,
)
    block = zeros(Float64, length(left_support), length(right_support))
    left_states =
        Tuple(
            ParentGaussletBases._cartesian_unflat_index(index, axis_counts)
            for index in left_support
        )
    right_states =
        Tuple(
            ParentGaussletBases._cartesian_unflat_index(index, axis_counts)
            for index in right_support
        )
    for (row, left_state) in pairs(left_states)
        ix_left, iy_left, iz_left = left_state
        for (column, right_state) in pairs(right_states)
            ix_right, iy_right, iz_right = right_state
            value = 0.0
            @inbounds for term in eachindex(coefficients)
                value +=
                    coefficients[term] *
                    axis_terms.x[term, ix_left, ix_right] *
                    axis_terms.y[term, iy_left, iy_right] *
                    axis_terms.z[term, iz_left, iz_right]
            end
            block[row, column] = value
        end
    end
    return block
end

function _white_lindsey_electron_nuclear_axis_term_shapes(axis_terms)
    return (x = size(axis_terms.x), y = size(axis_terms.y), z = size(axis_terms.z))
end
