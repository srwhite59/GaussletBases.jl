"""
    MappedPGDGPrototype1D

Prototype one-dimensional pure-Gaussian approximation to a
`MappedUniformBasis`.

The idea is the same as in the older PGDG work: replace each explicitly
distorted primitive

```text
g(u(x)) * sqrt(du/dx)
```

by a locally matched plain Gaussian in physical space. The local amplitude
`sqrt(du/dx)` is folded into the contraction matrix, so the resulting primitive
layer can reuse the analytic Gaussian matrix builders already present in the
package.
"""
struct MappedPGDGPrototype1D
    source_basis::MappedUniformBasis
    primitive_layer::PrimitiveSet1D
    coefficient_matrix::Matrix{Float64}
    primitive_width_data::Vector{Float64}
    primitive_amplitude_data::Vector{Float64}
end

"""
    MappedPGDGLocalized1D

Experimental post-cleanup / COMX-localized follow-on to
`MappedPGDGPrototype1D`.

This is still not the full historical PGDG path, but it adds the first
basis-cleanup and localization stage after the analytic primitive proxy has
been built.
"""
struct MappedPGDGLocalized1D
    source_prototype::MappedPGDGPrototype1D
    coefficient_matrix::Matrix{Float64}
    center_data::Vector{Float64}
    integral_weight_data::Vector{Float64}
end

"""
    MappedPGDGLogFitPrototype1D

Internal diagnostic proxy for `MappedUniformBasis`, using a sampled short-window
weighted log-quadratic Gaussian fit for each distorted primitive. Public PGDG
production backends do not call this helper.
"""
struct MappedPGDGLogFitPrototype1D
    source_basis::MappedUniformBasis
    primitive_layer::PrimitiveSet1D
    coefficient_matrix::Matrix{Float64}
    primitive_width_data::Vector{Float64}
    primitive_amplitude_data::Vector{Float64}
    primitive_center_data::Vector{Float64}
end

function Base.show(io::IO, prototype::MappedPGDGPrototype1D)
    print(
        io,
        "MappedPGDGPrototype1D(experimental, nbasis=",
        length(prototype),
        ", nprimitive=",
        length(prototype.primitive_layer),
        ", family=:",
        family(prototype.source_basis).name,
        ")",
    )
end

function Base.show(io::IO, localized::MappedPGDGLocalized1D)
    print(
        io,
        "MappedPGDGLocalized1D(experimental, nbasis=",
        length(localized),
        ", nprimitive=",
        length(primitive_set(localized)),
        ", family=:",
        family(localized.source_prototype.source_basis).name,
        ")",
    )
end

function Base.show(io::IO, prototype::MappedPGDGLogFitPrototype1D)
    print(
        io,
        "MappedPGDGLogFitPrototype1D(experimental, nbasis=",
        length(prototype),
        ", nprimitive=",
        length(prototype.primitive_layer),
        ", family=:",
        family(prototype.source_basis).name,
        ")",
    )
end

Base.length(prototype::MappedPGDGPrototype1D) = size(prototype.coefficient_matrix, 2)
Base.length(localized::MappedPGDGLocalized1D) = size(localized.coefficient_matrix, 2)
Base.length(prototype::MappedPGDGLogFitPrototype1D) = size(prototype.coefficient_matrix, 2)

primitive_set(prototype::MappedPGDGPrototype1D) = prototype.primitive_layer
primitives(prototype::MappedPGDGPrototype1D) = primitives(prototype.primitive_layer)
stencil_matrix(prototype::MappedPGDGPrototype1D) = prototype.coefficient_matrix
primitive_set(localized::MappedPGDGLocalized1D) = primitive_set(localized.source_prototype)
primitives(localized::MappedPGDGLocalized1D) = primitives(primitive_set(localized))
stencil_matrix(localized::MappedPGDGLocalized1D) = localized.coefficient_matrix
primitive_set(prototype::MappedPGDGLogFitPrototype1D) = prototype.primitive_layer
primitives(prototype::MappedPGDGLogFitPrototype1D) = primitives(prototype.primitive_layer)
stencil_matrix(prototype::MappedPGDGLogFitPrototype1D) = prototype.coefficient_matrix

basis_spec(prototype::MappedPGDGPrototype1D) = basis_spec(prototype.source_basis)
family(prototype::MappedPGDGPrototype1D) = family(prototype.source_basis)
mapping(prototype::MappedPGDGPrototype1D) = mapping(prototype.source_basis)
centers(prototype::MappedPGDGPrototype1D) = centers(prototype.source_basis)
reference_centers(prototype::MappedPGDGPrototype1D) = reference_centers(prototype.source_basis)
basis_spec(localized::MappedPGDGLocalized1D) = basis_spec(localized.source_prototype.source_basis)
family(localized::MappedPGDGLocalized1D) = family(localized.source_prototype.source_basis)
mapping(localized::MappedPGDGLocalized1D) = mapping(localized.source_prototype.source_basis)
centers(localized::MappedPGDGLocalized1D) = localized.center_data
reference_centers(localized::MappedPGDGLocalized1D) = localized.center_data
basis_spec(prototype::MappedPGDGLogFitPrototype1D) = basis_spec(prototype.source_basis)
family(prototype::MappedPGDGLogFitPrototype1D) = family(prototype.source_basis)
mapping(prototype::MappedPGDGLogFitPrototype1D) = mapping(prototype.source_basis)
centers(prototype::MappedPGDGLogFitPrototype1D) = centers(prototype.source_basis)
reference_centers(prototype::MappedPGDGLogFitPrototype1D) = reference_centers(prototype.source_basis)

function integral_weights(prototype::MappedPGDGPrototype1D)
    primitive_weights = Float64[integral_weight(primitive) for primitive in primitives(prototype)]
    return vec(transpose(stencil_matrix(prototype)) * primitive_weights)
end

integral_weights(localized::MappedPGDGLocalized1D) = localized.integral_weight_data

function integral_weights(prototype::MappedPGDGLogFitPrototype1D)
    primitive_weights = Float64[integral_weight(primitive) for primitive in primitives(prototype)]
    return vec(transpose(stencil_matrix(prototype)) * primitive_weights)
end

function _mapped_pgdg_primitive_data(basis::MappedUniformBasis)
    gaussian_primitives = Gaussian[]
    primitive_widths = Float64[]
    primitive_amplitudes = Float64[]

    for primitive in primitives(basis)
        if primitive isa Distorted
            primitive.primitive isa Gaussian || throw(ArgumentError("mapped_pgdg_prototype currently requires distorted full-line Gaussian primitives"))

            mapped_center = center(primitive)
            local_density = dudx(primitive.mapping, mapped_center)
            local_density > 0.0 || throw(ArgumentError("mapped_pgdg_prototype requires a positive local mapping density"))

            width_value = primitive.primitive.width / local_density
            amplitude_value = sqrt(local_density)
            push!(gaussian_primitives, Gaussian(center = mapped_center, width = width_value))
            push!(primitive_widths, width_value)
            push!(primitive_amplitudes, amplitude_value)
        elseif primitive isa Gaussian
            push!(gaussian_primitives, primitive)
            push!(primitive_widths, primitive.width)
            push!(primitive_amplitudes, 1.0)
        else
            throw(ArgumentError("mapped_pgdg_prototype currently requires full-line Gaussian primitives"))
        end
    end

    return gaussian_primitives, primitive_widths, primitive_amplitudes
end

"""
    mapped_pgdg_prototype(basis::MappedUniformBasis)

Build the narrow PGDG-style analytic prototype associated with `basis`.

This is not a new production basis family. It is a locally linearized
plain-Gaussian proxy used to test whether the mapped ordinary one-body path can
be rebuilt on analytic primitive integrals and the existing contraction layer.
"""
function mapped_pgdg_prototype(basis::MappedUniformBasis)
    gaussian_primitives, primitive_widths, primitive_amplitudes = _mapped_pgdg_primitive_data(basis)
    primitive_layer = PrimitiveSet1D(
        AbstractPrimitiveFunction1D[primitive for primitive in gaussian_primitives];
        name = :mapped_pgdg_primitives,
    )
    amplitude_matrix = Diagonal(primitive_amplitudes)
    coefficient_matrix = amplitude_matrix * Matrix{Float64}(stencil_matrix(basis))
    return MappedPGDGPrototype1D(
        basis,
        primitive_layer,
        coefficient_matrix,
        primitive_widths,
        primitive_amplitudes,
    )
end

function _logfit_gaussian_proxy(primitive::Gaussian)
    return primitive, primitive.center_value, primitive.width, 1.0
end

_derivativefit_gaussian_proxy(primitive::Gaussian) = _logfit_gaussian_proxy(primitive)

function _logfit_gaussian_proxy(primitive::Distorted{<:Gaussian})
    x_center = center(primitive)
    local_density = dudx(primitive.mapping, x_center)
    local_density > 0.0 || throw(ArgumentError("mapped PGDG log-fit prototype requires a positive local mapping density"))

    sigma_linear = primitive.primitive.width / local_density
    ulo, uhi = _reference_bounds(primitive.primitive)
    xlo = xofu(primitive.mapping, ulo)
    xhi = xofu(primitive.mapping, uhi)
    window_sigma = 2.0
    left = max(xlo, x_center - window_sigma * sigma_linear)
    right = min(xhi, x_center + window_sigma * sigma_linear)
    right > left || return Gaussian(center = x_center, width = sigma_linear), x_center, sigma_linear, sqrt(local_density)

    sample_count = 121
    xs = collect(range(left, right; length = sample_count))
    values = Float64[value(primitive, x) for x in xs]
    weights = values .^ 2
    all(iszero, weights) && return Gaussian(center = x_center, width = sigma_linear), x_center, sigma_linear, sqrt(local_density)

    deltas = xs .- x_center
    design = hcat(ones(length(xs)), deltas, deltas .^ 2)
    logs = log.(max.(values, 1.0e-300))
    scaled_design = design .* sqrt.(weights)
    scaled_logs = logs .* sqrt.(weights)
    coefficients = scaled_design \ scaled_logs
    _, linear_term, quadratic_term = coefficients

    quadratic_term < 0.0 || return Gaussian(center = x_center, width = sigma_linear), x_center, sigma_linear, sqrt(local_density)

    width_value = sqrt(-1.0 / (2.0 * quadratic_term))
    shift_value = linear_term * width_value^2
    center_value = x_center + shift_value
    amplitude_value = exp(coefficients[1] + shift_value^2 / (2.0 * width_value^2))
    gaussian = Gaussian(center = center_value, width = width_value)
    return gaussian, center_value, width_value, amplitude_value
end

function _derivativefit_gaussian_proxy(primitive::Distorted{<:Gaussian})
    x_center = center(primitive)
    local_density = dudx(primitive.mapping, x_center)
    local_density > 0.0 || throw(ArgumentError("mapped PGDG derivative-fit prototype requires a positive local mapping density"))

    sigma_linear = primitive.primitive.width / local_density
    ulo, uhi = _reference_bounds(primitive.primitive)
    xlo = xofu(primitive.mapping, ulo)
    xhi = xofu(primitive.mapping, uhi)
    window_sigma = 2.0
    left = max(xlo, x_center - window_sigma * sigma_linear)
    right = min(xhi, x_center + window_sigma * sigma_linear)
    right > left || return _logfit_gaussian_proxy(primitive)

    sample_count = 121
    xs = collect(range(left, right; length = sample_count))
    values = Float64[value(primitive, x) for x in xs]
    derivatives = Float64[derivative(primitive, x; order = 1) for x in xs]
    weights = values .^ 2
    maximum(weights) <= 0.0 && return _logfit_gaussian_proxy(primitive)

    keep = weights .> (maximum(weights) * 1.0e-10)
    count(keep) ≥ 5 || return _logfit_gaussian_proxy(primitive)

    xs_fit = xs[keep]
    deltas = xs_fit .- x_center
    weights_fit = weights[keep]
    ratios = derivatives[keep] ./ values[keep]
    design = hcat(ones(length(xs_fit)), deltas)
    scaled_design = design .* sqrt.(weights_fit)
    scaled_ratios = ratios .* sqrt.(weights_fit)
    coefficients = scaled_design \ scaled_ratios
    intercept, slope = coefficients

    slope < 0.0 || return _logfit_gaussian_proxy(primitive)

    width_value = sqrt(-1.0 / slope)
    isfinite(width_value) || return _logfit_gaussian_proxy(primitive)
    width_value > 0.0 || return _logfit_gaussian_proxy(primitive)

    center_value = x_center + intercept * width_value^2
    center_value = clamp(center_value, left, right)

    shift = ((xs_fit .- center_value) .^ 2) ./ (2.0 * width_value^2)
    amplitude_log = sum(weights_fit .* (log.(max.(values[keep], 1.0e-300)) .+ shift)) / sum(weights_fit)
    amplitude_value = exp(amplitude_log)
    isfinite(amplitude_value) || return _logfit_gaussian_proxy(primitive)

    gaussian = Gaussian(center = center_value, width = width_value)
    return gaussian, center_value, width_value, amplitude_value
end

function _mapped_pgdg_logfit_primitive_data(basis::MappedUniformBasis)
    gaussian_primitives = Gaussian[]
    primitive_centers = Float64[]
    primitive_widths = Float64[]
    primitive_amplitudes = Float64[]

    for primitive in primitives(basis)
        if primitive isa Distorted
            primitive.primitive isa Gaussian || throw(ArgumentError("mapped_pgdg_logfit_prototype currently requires distorted full-line Gaussian primitives"))
        elseif !(primitive isa Gaussian)
            throw(ArgumentError("mapped_pgdg_logfit_prototype currently requires full-line Gaussian primitives"))
        end

        gaussian, center_value, width_value, amplitude_value = _logfit_gaussian_proxy(primitive)
        push!(gaussian_primitives, gaussian)
        push!(primitive_centers, center_value)
        push!(primitive_widths, width_value)
        push!(primitive_amplitudes, amplitude_value)
    end

    return gaussian_primitives, primitive_centers, primitive_widths, primitive_amplitudes
end

function _mapped_pgdg_derivativefit_primitive_data(basis::MappedUniformBasis)
    gaussian_primitives = Gaussian[]
    primitive_centers = Float64[]
    primitive_widths = Float64[]
    primitive_amplitudes = Float64[]

    for primitive in primitives(basis)
        if primitive isa Distorted
            primitive.primitive isa Gaussian || throw(ArgumentError("mapped_pgdg_derivativefit_prototype currently requires distorted full-line Gaussian primitives"))
        elseif !(primitive isa Gaussian)
            throw(ArgumentError("mapped_pgdg_derivativefit_prototype currently requires full-line Gaussian primitives"))
        end

        gaussian, center_value, width_value, amplitude_value = _derivativefit_gaussian_proxy(primitive)
        push!(gaussian_primitives, gaussian)
        push!(primitive_centers, center_value)
        push!(primitive_widths, width_value)
        push!(primitive_amplitudes, amplitude_value)
    end

    return gaussian_primitives, primitive_centers, primitive_widths, primitive_amplitudes
end

"""
    mapped_pgdg_logfit_prototype(basis::MappedUniformBasis)

Internal diagnostic proxy for `basis`, using sampled short-window weighted
log-quadratic Gaussian replacements for the distorted primitives.
"""
function mapped_pgdg_logfit_prototype(basis::MappedUniformBasis)
    gaussian_primitives, primitive_centers, primitive_widths, primitive_amplitudes =
        _mapped_pgdg_logfit_primitive_data(basis)
    primitive_layer = PrimitiveSet1D(
        AbstractPrimitiveFunction1D[primitive for primitive in gaussian_primitives];
        name = :mapped_pgdg_logfit_primitives,
    )
    amplitude_matrix = Diagonal(primitive_amplitudes)
    coefficient_matrix = amplitude_matrix * Matrix{Float64}(stencil_matrix(basis))
    return MappedPGDGLogFitPrototype1D(
        basis,
        primitive_layer,
        coefficient_matrix,
        primitive_widths,
        primitive_amplitudes,
        primitive_centers,
    )
end

"""
    mapped_pgdg_derivativefit_prototype(basis::MappedUniformBasis)

Internal diagnostic proxy for `basis`, using sampled short-window value and
curvature matching for the distorted primitives. Public PGDG production
backends do not call this helper.
"""
function mapped_pgdg_derivativefit_prototype(basis::MappedUniformBasis)
    gaussian_primitives, primitive_centers, primitive_widths, primitive_amplitudes =
        _mapped_pgdg_derivativefit_primitive_data(basis)
    primitive_layer = PrimitiveSet1D(
        AbstractPrimitiveFunction1D[primitive for primitive in gaussian_primitives];
        name = :mapped_pgdg_derivativefit_primitives,
    )
    amplitude_matrix = Diagonal(primitive_amplitudes)
    coefficient_matrix = amplitude_matrix * Matrix{Float64}(stencil_matrix(basis))
    return MappedPGDGLogFitPrototype1D(
        basis,
        primitive_layer,
        coefficient_matrix,
        primitive_widths,
        primitive_amplitudes,
        primitive_centers,
    )
end

function _sign_fix_transform_columns!(transform::Matrix{Float64}, sign_vector::AbstractVector{<:Real})
    sign_values = Float64[Float64(value) for value in sign_vector]
    for column in 1:size(transform, 2)
        sign_probe = dot(sign_values, view(transform, :, column))
        if abs(sign_probe) <= 1.0e-12
            pivot = argmax(abs.(view(transform, :, column)))
            sign_probe = transform[pivot, column]
        end
        if sign_probe < 0.0
            transform[:, column] .*= -1.0
        end
    end
    return transform
end

function _cleanup_comx_transform(
    overlap::AbstractMatrix{<:Real},
    position::AbstractMatrix{<:Real},
    sign_vector::AbstractVector{<:Real},
)
    vectors, invhalf = _s_invsqrt_reduced(overlap)
    localizer, center_values = _comx_reduced(vectors, invhalf, position)
    transform = Matrix{Float64}(vectors * (invhalf * localizer))
    ordering = sortperm(center_values)
    transform = transform[:, ordering]
    centers_ordered = Float64[Float64(center_values[index]) for index in ordering]
    _sign_fix_transform_columns!(transform, sign_vector)
    return transform, centers_ordered
end

"""
    mapped_pgdg_localized(basis::MappedUniformBasis)
    mapped_pgdg_localized(prototype::MappedPGDGPrototype1D)

Experimental post-cleanup / COMX-localized follow-on to the pre-COMX analytic
PGDG-style prototype.

This adds the first overlap cleanup and position-localization step on top of
the analytic primitive proxy, while still stopping well short of the larger
historical PGDG driver logic.
"""
function mapped_pgdg_localized(prototype::MappedPGDGPrototype1D)
    overlap = overlap_matrix(prototype)
    position = position_matrix(prototype)
    transform, center_values = _cleanup_comx_transform(overlap, position, integral_weights(prototype))
    coefficient_matrix = stencil_matrix(prototype) * transform
    primitive_weights = Float64[integral_weight(primitive) for primitive in primitives(prototype)]
    integral_weight_data = vec(transpose(coefficient_matrix) * primitive_weights)
    return MappedPGDGLocalized1D(
        prototype,
        coefficient_matrix,
        center_values,
        integral_weight_data,
    )
end

mapped_pgdg_localized(basis::MappedUniformBasis) = mapped_pgdg_localized(mapped_pgdg_prototype(basis))

function mapped_pgdg_localized(prototype::MappedPGDGLogFitPrototype1D)
    overlap = overlap_matrix(prototype)
    position = position_matrix(prototype)
    transform, center_values = _cleanup_comx_transform(overlap, position, integral_weights(prototype))
    coefficient_matrix = stencil_matrix(prototype) * transform
    primitive_weights = Float64[integral_weight(primitive) for primitive in primitives(prototype)]
    integral_weight_data = vec(transpose(coefficient_matrix) * primitive_weights)
    return MappedPGDGLocalized1D(
        MappedPGDGPrototype1D(
            prototype.source_basis,
            prototype.primitive_layer,
            prototype.coefficient_matrix,
            prototype.primitive_width_data,
            prototype.primitive_amplitude_data,
        ),
        coefficient_matrix,
        center_values,
        integral_weight_data,
    )
end

function overlap_matrix(prototype::MappedPGDGPrototype1D)
    return contract_primitive_matrix(prototype, overlap_matrix(primitive_set(prototype)))
end

function overlap_matrix(localized::MappedPGDGLocalized1D)
    return contract_primitive_matrix(localized, overlap_matrix(primitive_set(localized)))
end

function overlap_matrix(prototype::MappedPGDGLogFitPrototype1D)
    return contract_primitive_matrix(prototype, overlap_matrix(primitive_set(prototype)))
end

function position_matrix(prototype::MappedPGDGPrototype1D)
    return contract_primitive_matrix(prototype, position_matrix(primitive_set(prototype)))
end

function position_matrix(localized::MappedPGDGLocalized1D)
    return contract_primitive_matrix(localized, position_matrix(primitive_set(localized)))
end

function position_matrix(prototype::MappedPGDGLogFitPrototype1D)
    return contract_primitive_matrix(prototype, position_matrix(primitive_set(prototype)))
end

function kinetic_matrix(prototype::MappedPGDGPrototype1D)
    return contract_primitive_matrix(prototype, kinetic_matrix(primitive_set(prototype)))
end

function kinetic_matrix(localized::MappedPGDGLocalized1D)
    return contract_primitive_matrix(localized, kinetic_matrix(primitive_set(localized)))
end

function kinetic_matrix(prototype::MappedPGDGLogFitPrototype1D)
    return contract_primitive_matrix(prototype, kinetic_matrix(primitive_set(prototype)))
end

function gaussian_factor_matrix(
    prototype::MappedPGDGPrototype1D;
    exponent::Real,
    center::Real = 0.0,
)
    primitive_matrix = gaussian_factor_matrix(
        primitive_set(prototype);
        exponent = exponent,
        center = center,
    )
    return contract_primitive_matrix(prototype, primitive_matrix)
end

function gaussian_factor_matrix(
    localized::MappedPGDGLocalized1D;
    exponent::Real,
    center::Real = 0.0,
)
    primitive_matrix = gaussian_factor_matrix(
        primitive_set(localized);
        exponent = exponent,
        center = center,
    )
    return contract_primitive_matrix(localized, primitive_matrix)
end

function gaussian_factor_matrix(
    prototype::MappedPGDGLogFitPrototype1D;
    exponent::Real,
    center::Real = 0.0,
)
    primitive_matrix = gaussian_factor_matrix(
        primitive_set(prototype);
        exponent = exponent,
        center = center,
    )
    return contract_primitive_matrix(prototype, primitive_matrix)
end

function gaussian_factor_matrices(
    prototype::MappedPGDGPrototype1D;
    exponents::AbstractVector{<:Real},
    center::Real = 0.0,
)
    primitive_matrices = gaussian_factor_matrices(
        primitive_set(prototype);
        exponents = exponents,
        center = center,
    )
    return [contract_primitive_matrix(prototype, matrix) for matrix in primitive_matrices]
end

function gaussian_factor_matrices(
    localized::MappedPGDGLocalized1D;
    exponents::AbstractVector{<:Real},
    center::Real = 0.0,
)
    primitive_matrices = gaussian_factor_matrices(
        primitive_set(localized);
        exponents = exponents,
        center = center,
    )
    return [contract_primitive_matrix(localized, matrix) for matrix in primitive_matrices]
end

function gaussian_factor_matrices(
    prototype::MappedPGDGLogFitPrototype1D;
    exponents::AbstractVector{<:Real},
    center::Real = 0.0,
)
    primitive_matrices = gaussian_factor_matrices(
        primitive_set(prototype);
        exponents = exponents,
        center = center,
    )
    return [contract_primitive_matrix(prototype, matrix) for matrix in primitive_matrices]
end
