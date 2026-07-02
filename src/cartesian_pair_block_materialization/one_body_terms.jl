# Shared one-body term descriptors for future mixed pair-block consumers.
#
# This file only describes term vocabulary and selector surfaces. It does not
# construct 1D factors, pair blocks, Hamiltonian data, exports, artifacts, or
# PQS Lowdin realization.

const _ONE_BODY_TERMS = (
    :overlap,
    :position_x,
    :position_y,
    :position_z,
    :x2_x,
    :x2_y,
    :x2_z,
    :kinetic,
)

# Block-set consumers use this descriptor to keep one-body results term-separated.
# Callers own factor supply; this metadata layer does not build factors or
# assemble operator/Hamiltonian sums.
function _one_body_term_set_descriptor()
    return _one_body_term_set_descriptor(_ONE_BODY_TERMS)
end

function _one_body_term_set_descriptor(term::Symbol)
    return _one_body_term_set_descriptor((term,))
end

function _one_body_term_set_descriptor(terms::AbstractVector)
    return _one_body_term_set_descriptor(Tuple(terms))
end

function _one_body_term_set_descriptor(terms::Tuple)
    isempty(terms) &&
        throw(ArgumentError("one-body term set must contain at least one term"))
    descriptors = Tuple(_one_body_term_descriptor(term) for term in terms)
    return (;
        object_kind = :cartesian_pair_block_one_body_term_set_descriptor,
        status = :available_one_body_term_set_descriptor,
        requested_terms = terms,
        terms,
        term_count = length(terms),
        term_descriptors = descriptors,
        term_families = Tuple(descriptor.family for descriptor in descriptors),
        required_factor_names =
            _one_body_ordered_unique_required_factor_names(descriptors),
        result_terms_remain_separated = true,
        block_set_results_summed = false,
        factor_provider_scope = :caller_supplied_or_family_provider,
    )
end

function _one_body_term_set_descriptor(terms)
    throw(
        ArgumentError(
            "one-body term set must be a Symbol, Tuple, or AbstractVector, got $(typeof(terms))",
        ),
    )
end

function _one_body_term_descriptor(term)
    throw(ArgumentError("one-body term must be a Symbol, got $(typeof(term))"))
end

function _one_body_term_descriptor(term::Symbol)
    term === :overlap && return _one_body_term_descriptor_result(
        term,
        :overlap,
        :one_body_overlap,
        nothing,
        (:overlap,),
        (:overlap_1d,),
    )

    if term in (:position_x, :position_y, :position_z)
        axis = _one_body_axis_for_prefixed_term(term, "position_")
        return _one_body_term_descriptor_result(
            term,
            :position,
            :one_body_axis_position,
            axis,
            (:overlap, :position),
            (:overlap_1d, :position_1d),
        )
    end

    if term in (:x2_x, :x2_y, :x2_z)
        axis = _one_body_axis_for_prefixed_term(term, "x2_")
        return _one_body_term_descriptor_result(
            term,
            :x2,
            :one_body_axis_x2,
            axis,
            (:overlap, :x2),
            (:overlap_1d, :x2_1d),
        )
    end

    term === :kinetic && return _one_body_term_descriptor_result(
        term,
        :kinetic,
        :one_body_cartesian_kinetic_sum,
        nothing,
        (:overlap, :kinetic),
        (:overlap_1d, :kinetic_1d),
    )

    throw(ArgumentError("unsupported one-body term: $(term)"))
end

function _one_body_term_descriptor_result(
    requested_term::Symbol,
    family::Symbol,
    term_kind::Symbol,
    axis,
    required_factor_roles::Tuple,
    required_factor_names::Tuple,
)
    return (;
        object_kind = :cartesian_pair_block_one_body_term_descriptor,
        status = :available_one_body_term_descriptor,
        requested_term,
        term = requested_term,
        family,
        term_family = family,
        term_kind,
        axis,
        axis_index = _one_body_axis_index(axis),
        required_factor_roles,
        required_factor_names,
        factor_provider_scope = :caller_supplied_or_family_provider,
    )
end

function _one_body_ordered_unique_required_factor_names(descriptors::Tuple)
    names = Symbol[]
    for descriptor in descriptors
        for name in descriptor.required_factor_names
            name in names && continue
            push!(names, name)
        end
    end
    return Tuple(names)
end

function _one_body_axis_for_prefixed_term(term::Symbol, prefix::AbstractString)
    suffix = replace(String(term), prefix => ""; count = 1)
    suffix == "x" && return :x
    suffix == "y" && return :y
    suffix == "z" && return :z
    throw(ArgumentError("unsupported axis one-body term: $(term)"))
end

function _one_body_axis_index(axis)
    isnothing(axis) && return nothing
    axis === :x && return 1
    axis === :y && return 2
    axis === :z && return 3
    throw(ArgumentError("unsupported one-body axis: $(axis)"))
end
