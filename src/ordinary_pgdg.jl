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

Base.length(prototype::MappedPGDGPrototype1D) = size(prototype.coefficient_matrix, 2)
Base.length(localized::MappedPGDGLocalized1D) = size(localized.coefficient_matrix, 2)

primitive_set(prototype::MappedPGDGPrototype1D) = prototype.primitive_layer
primitives(prototype::MappedPGDGPrototype1D) = primitives(prototype.primitive_layer)
stencil_matrix(prototype::MappedPGDGPrototype1D) = prototype.coefficient_matrix
primitive_set(localized::MappedPGDGLocalized1D) = primitive_set(localized.source_prototype)
primitives(localized::MappedPGDGLocalized1D) = primitives(primitive_set(localized))
stencil_matrix(localized::MappedPGDGLocalized1D) = localized.coefficient_matrix

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

function integral_weights(prototype::MappedPGDGPrototype1D)
    primitive_weights = Float64[integral_weight(primitive) for primitive in primitives(prototype)]
    return vec(transpose(stencil_matrix(prototype)) * primitive_weights)
end

integral_weights(localized::MappedPGDGLocalized1D) = localized.integral_weight_data

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

function overlap_matrix(prototype::MappedPGDGPrototype1D)
    return contract_primitive_matrix(prototype, overlap_matrix(primitive_set(prototype)))
end

function overlap_matrix(localized::MappedPGDGLocalized1D)
    return contract_primitive_matrix(localized, overlap_matrix(primitive_set(localized)))
end

function position_matrix(prototype::MappedPGDGPrototype1D)
    return contract_primitive_matrix(prototype, position_matrix(primitive_set(prototype)))
end

function position_matrix(localized::MappedPGDGLocalized1D)
    return contract_primitive_matrix(localized, position_matrix(primitive_set(localized)))
end

function kinetic_matrix(prototype::MappedPGDGPrototype1D)
    return contract_primitive_matrix(prototype, kinetic_matrix(primitive_set(prototype)))
end

function kinetic_matrix(localized::MappedPGDGLocalized1D)
    return contract_primitive_matrix(localized, kinetic_matrix(primitive_set(localized)))
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
