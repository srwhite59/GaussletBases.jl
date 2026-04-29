"""
    CartesianBasisMetadata3D

Minimal public metadata bundle identifying a three-dimensional Cartesian basis.

This is the Cartesian analogue of `BasisMetadata1D`. It records:

- the high-level basis kind, such as `:direct_product`, `:nested_fixed_block`,
  or `:hybrid_residual`
- the axis-level one-dimensional metadata reused from the existing 1D contract
- the defining parent-space contract and final basis size
- light basis-level geometric data such as labels, centers, and working-box
  metadata when that construction metadata already exists

The object is basis-defining metadata only. It is not an operator cache.
"""
struct CartesianBasisMetadata3D{AT <: NamedTuple, RT <: NamedTuple}
    basis_kind::Symbol
    axis_sharing::Symbol
    axis_metadata::AT
    parent_kind::Symbol
    parent_axis_counts::NTuple{3,Int}
    parent_dimension::Int
    final_dimension::Int
    working_box::Union{Nothing,NTuple{3,UnitRange{Int}}}
    basis_labels::Vector{String}
    basis_centers::Matrix{Float64}
    route_metadata::RT
end

function Base.show(io::IO, metadata::CartesianBasisMetadata3D)
    print(
        io,
        "CartesianBasisMetadata3D(kind=:",
        metadata.basis_kind,
        ", nfinal=",
        metadata.final_dimension,
        ", parent_kind=:",
        metadata.parent_kind,
        ", axis_sharing=:",
        metadata.axis_sharing,
        ")",
    )
end

"""
    CartesianBasisRepresentation3D

Public in-memory three-dimensional Cartesian basis representation.

The representation stores:

- `metadata`
- explicit 1D axis representations reused from the existing 1D contract
- the contraction from the stored parent/raw basis to the final basis
- parent-space labels and centers
- support-space indexing data when the construction is support-local
- optional parent-data sidecars needed to make mixed raw-space contracts
  explicit, for example the Cartesian parent representation and supplement
  orbital metadata in the QW residual-Gaussian route

This first-pass contract is representation-only. It deliberately avoids any
operator cache beyond the basis-defining contraction data itself.
"""
struct CartesianBasisRepresentation3D{
    MT <: CartesianBasisMetadata3D,
    AT <: NamedTuple,
    PT <: NamedTuple,
}
    metadata::MT
    axis_representations::AT
    contraction_kind::Symbol
    coefficient_matrix::Union{Nothing,_CartesianCoefficientMap}
    parent_labels::Vector{String}
    parent_centers::Matrix{Float64}
    support_indices::Union{Nothing,Vector{Int}}
    support_states::Union{Nothing,Vector{NTuple{3,Int}}}
    parent_data::PT
end

function Base.show(io::IO, representation::CartesianBasisRepresentation3D)
    print(
        io,
        "CartesianBasisRepresentation3D(kind=:",
        representation.metadata.basis_kind,
        ", nfinal=",
        representation.metadata.final_dimension,
        ", parent_kind=:",
        representation.metadata.parent_kind,
        ", contraction=:",
        representation.contraction_kind,
        ")",
    )
end

basis_metadata(representation::CartesianBasisRepresentation3D) = representation.metadata
basis_representation(representation::CartesianBasisRepresentation3D) = representation

"""
    CartesianGaussianShellOrbitalRepresentation3D

Public basis-defining representation of one contracted 3D Cartesian Gaussian
shell orbital used on the supplement side of the hybrid residual-Gaussian
routes.

The stored data is exact enough for overlap:

- `label`
- `angular_powers`
- `center`
- `exponents`
- `coefficients`
- `primitive_normalization`

The current normalization contract is
`:axiswise_normalized_cartesian_gaussian`, meaning each primitive uses the same
axiswise normalized Cartesian-Gaussian prefactor already used in the live
QW/nested residual construction.
"""
struct CartesianGaussianShellOrbitalRepresentation3D
    label::String
    angular_powers::NTuple{3,Int}
    center::NTuple{3,Float64}
    exponents::Vector{Float64}
    coefficients::Vector{Float64}
    primitive_normalization::Symbol
end

function Base.show(io::IO, orbital::CartesianGaussianShellOrbitalRepresentation3D)
    print(
        io,
        "CartesianGaussianShellOrbitalRepresentation3D(label=",
        repr(orbital.label),
        ", l=(",
        orbital.angular_powers[1],
        ",",
        orbital.angular_powers[2],
        ",",
        orbital.angular_powers[3],
        "), nprimitive=",
        length(orbital.exponents),
        ")",
    )
end

"""
    CartesianGaussianShellSupplementRepresentation3D

Public representation of the explicit supplement-orbital sector used in the
hybrid residual-Gaussian Cartesian routes.

It stores:

- `supplement_kind`
- explicit orbital representations
- light source metadata such as basis name, `lmax`, and nuclei when that data
  already exists on the legacy supplement object
"""
struct CartesianGaussianShellSupplementRepresentation3D{MT <: NamedTuple}
    supplement_kind::Symbol
    orbitals::Vector{CartesianGaussianShellOrbitalRepresentation3D}
    metadata::MT
end

function Base.show(io::IO, supplement::CartesianGaussianShellSupplementRepresentation3D)
    print(
        io,
        "CartesianGaussianShellSupplementRepresentation3D(kind=:",
        supplement.supplement_kind,
        ", norbitals=",
        length(supplement.orbitals),
        ")",
    )
end

basis_representation(supplement::CartesianGaussianShellSupplementRepresentation3D) = supplement
basis_metadata(supplement::CartesianGaussianShellSupplementRepresentation3D) = supplement.metadata

function _cartesian_product_states(axis_counts::NTuple{3,Int})
    nx, ny, nz = axis_counts
    states = Vector{NTuple{3,Int}}(undef, nx * ny * nz)
    index = 0
    for ix in 1:nx, iy in 1:ny, iz in 1:nz
        index += 1
        states[index] = (ix, iy, iz)
    end
    return states
end

function _cartesian_product_cross_overlap(
    left_states::AbstractVector{<:NTuple{3,Int}},
    right_states::AbstractVector{<:NTuple{3,Int}},
    overlap_x::AbstractMatrix{<:Real},
    overlap_y::AbstractMatrix{<:Real},
    overlap_z::AbstractMatrix{<:Real},
)
    matrix = zeros(Float64, length(left_states), length(right_states))
    for (row, (ix_left, iy_left, iz_left)) in pairs(left_states)
        for (column, (ix_right, iy_right, iz_right)) in pairs(right_states)
            matrix[row, column] =
                overlap_x[ix_left, ix_right] *
                overlap_y[iy_left, iy_right] *
                overlap_z[iz_left, iz_right]
        end
    end
    return matrix
end

function _cartesian_identity_factorized_basis(
    axis_counts::NTuple{3,Int},
)
    nx, ny, nz = axis_counts
    return _CartesianNestedFactorizedBasis3D(
        axis_counts,
        Matrix{Float64}(LinearAlgebra.I, nx, nx),
        Matrix{Float64}(LinearAlgebra.I, ny, ny),
        Matrix{Float64}(LinearAlgebra.I, nz, nz),
        _cartesian_product_states(axis_counts),
        ones(Float64, nx * ny * nz),
        0.0,
    )
end

function _cartesian_factorized_parent_basis(
    representation::CartesianBasisRepresentation3D,
)
    representation.metadata.parent_kind == :cartesian_product_basis || throw(
        ArgumentError(
            "Cartesian factorized parent extraction requires a pure Cartesian-product parent space, got :$(representation.metadata.parent_kind)",
        ),
    )
    if hasproperty(representation.parent_data, :factorized_cartesian_parent_basis)
        cached = representation.parent_data.factorized_cartesian_parent_basis
        cached isa _CartesianNestedFactorizedBasis3D && return cached
    end
    if representation.contraction_kind == :identity
        return _cartesian_identity_factorized_basis(representation.metadata.parent_axis_counts)
    elseif representation.contraction_kind == :dense
        representation.coefficient_matrix === nothing && throw(
            ArgumentError(
                "Cartesian factorized parent extraction requires an explicit coefficient matrix",
            ),
        )
        return _nested_extract_factorized_basis(
            representation.coefficient_matrix,
            representation.metadata.parent_axis_counts,
        )
    end
    throw(
        ArgumentError(
            "Cartesian factorized parent extraction does not support contraction kind :$(representation.contraction_kind)",
        ),
    )
end

function _cartesian_factorized_axis_cross_table(
    left_functions::AbstractMatrix{<:Real},
    left_axis::BasisRepresentation1D,
    right_functions::AbstractMatrix{<:Real},
    right_axis::BasisRepresentation1D,
)
    basis_cross = _basis_cross_overlap_1d(left_axis, right_axis)
    scratch = Matrix{Float64}(undef, size(left_functions, 2), size(basis_cross, 2))
    table = Matrix{Float64}(undef, size(left_functions, 2), size(right_functions, 2))
    mul!(scratch, transpose(left_functions), basis_cross)
    mul!(table, scratch, right_functions)
    return table
end

Base.@kwdef mutable struct _CartesianCrossOverlapStageTimer
    factorized_parent_ns::Int = 0
    axis_cross_table_ns::Int = 0
    cartesian_cartesian_block_ns::Int = 0
    cartesian_supplement_block_ns::Int = 0
    supplement_supplement_block_ns::Int = 0
    final_contraction_ns::Int = 0
end

function _cartesian_time_stage!(
    f::F,
    timer::_CartesianCrossOverlapStageTimer,
    field::Symbol,
) where {F}
    start_ns = time_ns()
    value = f()
    elapsed_ns = Int(time_ns() - start_ns)
    setfield!(timer, field, getfield(timer, field) + elapsed_ns)
    return value
end

function _cartesian_cross_overlap_stage_timings(
    timer::_CartesianCrossOverlapStageTimer,
    axis_cross_tables = nothing,
)
    return (
        factorized_parent_seconds = timer.factorized_parent_ns / 1.0e9,
        axis_cross_table_seconds = timer.axis_cross_table_ns / 1.0e9,
        axis_cross_tables = axis_cross_tables,
        cartesian_cartesian_block_seconds = timer.cartesian_cartesian_block_ns / 1.0e9,
        cartesian_supplement_block_seconds = timer.cartesian_supplement_block_ns / 1.0e9,
        supplement_supplement_block_seconds = timer.supplement_supplement_block_ns / 1.0e9,
        final_contraction_seconds = timer.final_contraction_ns / 1.0e9,
        total_seconds =
            (
                timer.factorized_parent_ns +
                timer.axis_cross_table_ns +
                timer.cartesian_cartesian_block_ns +
                timer.cartesian_supplement_block_ns +
                timer.supplement_supplement_block_ns +
                timer.final_contraction_ns
            ) / 1.0e9,
    )
end

function _cartesian_factorized_axis_cross_table_with_stage_timings(
    left_functions::AbstractMatrix{<:Real},
    left_axis::BasisRepresentation1D,
    right_functions::AbstractMatrix{<:Real},
    right_axis::BasisRepresentation1D,
)
    basis_cross_start_ns = time_ns()
    basis_cross = _basis_cross_overlap_1d(left_axis, right_axis)
    basis_cross_seconds = (time_ns() - basis_cross_start_ns) / 1.0e9

    left_contract_start_ns = time_ns()
    scratch = Matrix{Float64}(undef, size(left_functions, 2), size(basis_cross, 2))
    mul!(scratch, transpose(left_functions), basis_cross)
    left_contract_seconds = (time_ns() - left_contract_start_ns) / 1.0e9

    right_contract_start_ns = time_ns()
    table = Matrix{Float64}(undef, size(left_functions, 2), size(right_functions, 2))
    mul!(table, scratch, right_functions)
    right_contract_seconds = (time_ns() - right_contract_start_ns) / 1.0e9

    return (
        table = table,
        timings = (
            basis_cross_overlap_seconds = basis_cross_seconds,
            left_dictionary_contraction_seconds = left_contract_seconds,
            right_dictionary_contraction_seconds = right_contract_seconds,
            reused_from = nothing,
            total_seconds =
                basis_cross_seconds + left_contract_seconds + right_contract_seconds,
        ),
    )
end

function _cartesian_same_axis_cross_table_work(
    left_axis_a::BasisRepresentation1D,
    left_functions_a::AbstractMatrix{<:Real},
    right_axis_a::BasisRepresentation1D,
    right_functions_a::AbstractMatrix{<:Real},
    left_axis_b::BasisRepresentation1D,
    left_functions_b::AbstractMatrix{<:Real},
    right_axis_b::BasisRepresentation1D,
    right_functions_b::AbstractMatrix{<:Real},
)
    return _structural_isequal(left_axis_a.metadata, left_axis_b.metadata) &&
           _structural_isequal(right_axis_a.metadata, right_axis_b.metadata) &&
           isequal(left_functions_a, left_functions_b) &&
           isequal(right_functions_a, right_functions_b)
end

function _cartesian_factorized_axis_cross_tables(
    left::CartesianBasisRepresentation3D,
    right::CartesianBasisRepresentation3D,
    left_factorized::_CartesianNestedFactorizedBasis3D,
    right_factorized::_CartesianNestedFactorizedBasis3D,
)
    function _reused_axis_result(table, axis::Symbol)
        return (
            table = table,
            timings = (
                basis_cross_overlap_seconds = 0.0,
                left_dictionary_contraction_seconds = 0.0,
                right_dictionary_contraction_seconds = 0.0,
                reused_from = axis,
                total_seconds = 0.0,
            ),
        )
    end

    x = _cartesian_factorized_axis_cross_table_with_stage_timings(
        left_factorized.x_functions,
        left.axis_representations.x,
        right_factorized.x_functions,
        right.axis_representations.x,
    )
    y =
        if _cartesian_same_axis_cross_table_work(
            left.axis_representations.x,
            left_factorized.x_functions,
            right.axis_representations.x,
            right_factorized.x_functions,
            left.axis_representations.y,
            left_factorized.y_functions,
            right.axis_representations.y,
            right_factorized.y_functions,
        )
            _reused_axis_result(x.table, :x)
        else
            _cartesian_factorized_axis_cross_table_with_stage_timings(
                left_factorized.y_functions,
                left.axis_representations.y,
                right_factorized.y_functions,
                right.axis_representations.y,
            )
        end
    z =
        if _cartesian_same_axis_cross_table_work(
            left.axis_representations.x,
            left_factorized.x_functions,
            right.axis_representations.x,
            right_factorized.x_functions,
            left.axis_representations.z,
            left_factorized.z_functions,
            right.axis_representations.z,
            right_factorized.z_functions,
        )
            _reused_axis_result(x.table, :x)
        elseif _cartesian_same_axis_cross_table_work(
            left.axis_representations.y,
            left_factorized.y_functions,
            right.axis_representations.y,
            right_factorized.y_functions,
            left.axis_representations.z,
            left_factorized.z_functions,
            right.axis_representations.z,
            right_factorized.z_functions,
        )
            _reused_axis_result(y.table, :y)
        else
            _cartesian_factorized_axis_cross_table_with_stage_timings(
                left_factorized.z_functions,
                left.axis_representations.z,
                right_factorized.z_functions,
                right.axis_representations.z,
            )
        end
    return (
        tables = (x = x.table, y = y.table, z = z.table),
        timings = (x = x.timings, y = y.timings, z = z.timings),
    )
end

function _cartesian_factorized_basis_supplement_cross_from_axis_tables(
    factorized::_CartesianNestedFactorizedBasis3D,
    axis_tables::NamedTuple{(:x, :y, :z)},
    norbitals::Int,
)
    matrix = zeros(Float64, length(factorized.basis_triplets), norbitals)
    @inbounds for basis_index in eachindex(factorized.basis_triplets)
        ix, iy, iz = factorized.basis_triplets[basis_index]
        amplitude = factorized.basis_amplitudes[basis_index]
        for orbital_index in 1:norbitals
            matrix[basis_index, orbital_index] =
                amplitude *
                axis_tables.x[ix, orbital_index] *
                axis_tables.y[iy, orbital_index] *
                axis_tables.z[iz, orbital_index]
        end
    end
    return matrix
end

function _cartesian_factorized_supplement_axis_tables(
    basis::CartesianBasisRepresentation3D,
    supplement::CartesianGaussianShellSupplementRepresentation3D,
    candidates...,
)
    for candidate in candidates
        candidate === nothing && continue
        hasproperty(candidate, :cartesian_supplement_axis_tables) || continue
        axis_tables = candidate.cartesian_supplement_axis_tables
        axis_tables === nothing && continue
        _cartesian_same_basis_identity(basis, candidate.cartesian_representation) || continue
        _cartesian_same_supplement_raw_identity(
            supplement,
            candidate.supplement_representation,
        ) || continue
        return axis_tables
    end
    throw(
        ArgumentError(
            "exact hybrid cartesian-supplement overlap requires stored factorized axis tables on a hybrid raw side with matching Cartesian parent and supplement raw identities",
        ),
    )
end

function _cartesian_factorized_contract_axis_x(
    overlap_x::AbstractMatrix{<:Real},
    tensor::Array{Float64,4},
)
    contracted = overlap_x * reshape(tensor, size(tensor, 1), :)
    return reshape(
        contracted,
        size(overlap_x, 1),
        size(tensor, 2),
        size(tensor, 3),
        size(tensor, 4),
    )
end

function _cartesian_factorized_contract_axis_y(
    overlap_y::AbstractMatrix{<:Real},
    tensor::Array{Float64,4},
)
    permuted = PermutedDimsArray(tensor, (2, 1, 3, 4))
    contracted = overlap_y * reshape(permuted, size(tensor, 2), :)
    return permutedims(
        reshape(
            contracted,
            size(overlap_y, 1),
            size(tensor, 1),
            size(tensor, 3),
            size(tensor, 4),
        ),
        (2, 1, 3, 4),
    )
end

function _cartesian_factorized_contract_axis_z(
    overlap_z::AbstractMatrix{<:Real},
    tensor::Array{Float64,4},
)
    permuted = PermutedDimsArray(tensor, (3, 1, 2, 4))
    contracted = overlap_z * reshape(permuted, size(tensor, 3), :)
    return permutedims(
        reshape(
            contracted,
            size(overlap_z, 1),
            size(tensor, 1),
            size(tensor, 2),
            size(tensor, 4),
        ),
        (2, 3, 1, 4),
    )
end

function _cartesian_factorized_cartesian_cross_overlap(
    left_factorized::_CartesianNestedFactorizedBasis3D,
    right_factorized::_CartesianNestedFactorizedBasis3D,
    axis_cross_tables::NamedTuple{(:x, :y, :z)};
    block_columns::Int = 32,
)
    ncols = length(right_factorized.basis_triplets)
    left_dimension = length(left_factorized.basis_triplets)
    result = zeros(Float64, left_dimension, ncols)
    right_triplets = right_factorized.basis_triplets
    right_amplitudes = right_factorized.basis_amplitudes
    left_triplets = left_factorized.basis_triplets
    left_amplitudes = left_factorized.basis_amplitudes
    nx_right = size(right_factorized.x_functions, 2)
    ny_right = size(right_factorized.y_functions, 2)
    nz_right = size(right_factorized.z_functions, 2)

    for first_column in 1:block_columns:ncols
        last_column = min(ncols, first_column + block_columns - 1)
        block_range = first_column:last_column
        block_size = length(block_range)
        aggregated = zeros(Float64, nx_right, ny_right, nz_right, block_size)
        @inbounds for (local_column, global_column) in enumerate(block_range)
            ix, iy, iz = right_triplets[global_column]
            aggregated[ix, iy, iz, local_column] = right_amplitudes[global_column]
        end

        aggregated = _cartesian_factorized_contract_axis_x(axis_cross_tables.x, aggregated)
        aggregated = _cartesian_factorized_contract_axis_y(axis_cross_tables.y, aggregated)
        aggregated = _cartesian_factorized_contract_axis_z(axis_cross_tables.z, aggregated)

        @inbounds for basis_index in eachindex(left_triplets)
            ix, iy, iz = left_triplets[basis_index]
            amplitude = left_amplitudes[basis_index]
            for (local_column, global_column) in enumerate(block_range)
                result[basis_index, global_column] =
                    amplitude * aggregated[ix, iy, iz, local_column]
            end
        end
    end
    return result
end

function _cartesian_factorized_cartesian_cross_apply(
    left_factorized::_CartesianNestedFactorizedBasis3D,
    right_factorized::_CartesianNestedFactorizedBasis3D,
    axis_cross_tables::NamedTuple{(:x, :y, :z)},
    right_coefficients::AbstractMatrix{<:Real},
    ;
    block_columns::Int = 8,
)
    nright = length(right_factorized.basis_triplets)
    size(right_coefficients, 1) == nright || throw(
        DimensionMismatch(
            "factorized Cartesian cross apply expected $nright source Cartesian rows, got $(size(right_coefficients, 1))",
        ),
    )
    ncols = size(right_coefficients, 2)
    result = zeros(Float64, length(left_factorized.basis_triplets), ncols)
    right_triplets = right_factorized.basis_triplets
    right_amplitudes = right_factorized.basis_amplitudes
    left_triplets = left_factorized.basis_triplets
    left_amplitudes = left_factorized.basis_amplitudes
    nx_right = size(right_factorized.x_functions, 2)
    ny_right = size(right_factorized.y_functions, 2)
    nz_right = size(right_factorized.z_functions, 2)

    for first_column in 1:block_columns:ncols
        last_column = min(ncols, first_column + block_columns - 1)
        block_range = first_column:last_column
        block_size = length(block_range)
        aggregated = zeros(Float64, nx_right, ny_right, nz_right, block_size)
        @inbounds for source_index in eachindex(right_triplets)
            ix, iy, iz = right_triplets[source_index]
            amplitude = right_amplitudes[source_index]
            for (local_column, global_column) in enumerate(block_range)
                aggregated[ix, iy, iz, local_column] +=
                    amplitude * Float64(right_coefficients[source_index, global_column])
            end
        end

        aggregated = _cartesian_factorized_contract_axis_x(axis_cross_tables.x, aggregated)
        aggregated = _cartesian_factorized_contract_axis_y(axis_cross_tables.y, aggregated)
        aggregated = _cartesian_factorized_contract_axis_z(axis_cross_tables.z, aggregated)

        @inbounds for target_index in eachindex(left_triplets)
            ix, iy, iz = left_triplets[target_index]
            amplitude = left_amplitudes[target_index]
            for (local_column, global_column) in enumerate(block_range)
                result[target_index, global_column] =
                    amplitude * aggregated[ix, iy, iz, local_column]
            end
        end
    end
    return result
end

function _cartesian_factorized_basis_supplement_cross(
    factorized::_CartesianNestedFactorizedBasis3D,
    basis::CartesianBasisRepresentation3D,
    supplement::CartesianGaussianShellSupplementRepresentation3D,
    candidates...,
)
    isempty(supplement.orbitals) && return zeros(Float64, basis.metadata.final_dimension, 0)
    axis_tables =
        _cartesian_factorized_supplement_axis_tables(basis, supplement, candidates...)
    return _cartesian_factorized_basis_supplement_cross_from_axis_tables(
        factorized,
        axis_tables,
        length(supplement.orbitals),
    )
end

function _cartesian_factorized_basis_supplement_cross(
    basis::CartesianBasisRepresentation3D,
    supplement::CartesianGaussianShellSupplementRepresentation3D,
    candidates...,
)
    return _cartesian_factorized_basis_supplement_cross(
        _cartesian_factorized_parent_basis(basis),
        basis,
        supplement,
        candidates...,
    )
end

function _cartesian_same_parent_raw_identity(
    left::CartesianBasisRepresentation3D,
    right::CartesianBasisRepresentation3D,
)
    left.metadata.parent_kind == :cartesian_product_basis || return false
    right.metadata.parent_kind == :cartesian_product_basis || return false
    left.metadata.parent_axis_counts == right.metadata.parent_axis_counts || return false
    left.metadata.parent_dimension == right.metadata.parent_dimension || return false
    isequal(left.parent_labels, right.parent_labels) || return false
    isequal(left.parent_centers, right.parent_centers) || return false
    return true
end

function _cartesian_same_basis_identity(
    left::CartesianBasisRepresentation3D,
    right::CartesianBasisRepresentation3D,
)
    _cartesian_same_parent_raw_identity(left, right) || return false
    left.metadata.basis_kind == right.metadata.basis_kind || return false
    left.metadata.axis_sharing == right.metadata.axis_sharing || return false
    left.metadata.final_dimension == right.metadata.final_dimension || return false
    isequal(left.metadata.working_box, right.metadata.working_box) || return false
    isequal(left.metadata.basis_labels, right.metadata.basis_labels) || return false
    isequal(left.metadata.basis_centers, right.metadata.basis_centers) || return false
    isequal(left.metadata.route_metadata, right.metadata.route_metadata) || return false
    left.contraction_kind == right.contraction_kind || return false
    isequal(left.coefficient_matrix, right.coefficient_matrix) || return false
    isequal(left.support_indices, right.support_indices) || return false
    isequal(left.support_states, right.support_states) || return false
    return true
end

function _cartesian_supports_exact_hybrid_overlap(
    representation::CartesianBasisRepresentation3D,
)
    if representation.metadata.parent_kind == :cartesian_product_basis
        return true
    elseif representation.metadata.parent_kind == :cartesian_plus_supplement_raw
        hasproperty(representation.parent_data, :cartesian_parent_representation) || return false
        hasproperty(representation.parent_data, :supplement_representation) || return false
        hasproperty(representation.parent_data, :factorized_cartesian_parent_basis) || return false
        hasproperty(representation.parent_data, :cartesian_supplement_axis_tables) || return false
        parent_representation = representation.parent_data.cartesian_parent_representation
        supplement_representation = representation.parent_data.supplement_representation
        axis_tables = representation.parent_data.cartesian_supplement_axis_tables
        parent_representation isa CartesianBasisRepresentation3D || return false
        supplement_representation isa CartesianGaussianShellSupplementRepresentation3D || return false
        parent_representation.metadata.parent_kind == :cartesian_product_basis || return false
        supplement_representation.supplement_kind in (
            :atomic_cartesian_shell,
            :bond_aligned_diatomic_cartesian_shell,
            :bond_aligned_heteronuclear_cartesian_shell,
        ) || return false
        hasproperty(axis_tables, :x) || return false
        hasproperty(axis_tables, :y) || return false
        hasproperty(axis_tables, :z) || return false
        all(
            orbital -> orbital.primitive_normalization == :axiswise_normalized_cartesian_gaussian,
            supplement_representation.orbitals,
        ) || return false
        return true
    end
    return false
end

function _cartesian_cross_overlap_error(
    left::CartesianBasisRepresentation3D,
    right::CartesianBasisRepresentation3D,
)
    if left.metadata.parent_kind == :cartesian_plus_supplement_raw ||
       right.metadata.parent_kind == :cartesian_plus_supplement_raw
        return ArgumentError(
            "exact Cartesian cross overlap on hybrid residual-Gaussian representations is currently supported only for the atomic Cartesian-plus-supplement lane with explicit Cartesian parent and supplement orbital representations",
        )
    end
    return ArgumentError(
        "exact Cartesian cross overlap is not yet implemented for parent kinds :$(left.metadata.parent_kind) and :$(right.metadata.parent_kind)",
    )
end

function _cartesian_parent_state_basis(
    representation::CartesianBasisRepresentation3D,
)
    if representation.contraction_kind == :dense &&
       representation.support_indices !== nothing &&
       representation.support_states !== nothing
        return (
            states = representation.support_states,
            coefficients = Matrix{Float64}(representation.coefficient_matrix[representation.support_indices, :]),
        )
    elseif representation.contraction_kind == :dense
        return (
            states = _cartesian_product_states(representation.metadata.parent_axis_counts),
            coefficients = Matrix{Float64}(representation.coefficient_matrix),
        )
    elseif representation.contraction_kind == :identity
        return (
            states = _cartesian_product_states(representation.metadata.parent_axis_counts),
            coefficients = nothing,
        )
    end
    throw(
        ArgumentError(
            "unsupported Cartesian representation contraction kind :$(representation.contraction_kind)",
        ),
    )
end

function _cartesian_parent_cross_overlap(
    left::CartesianBasisRepresentation3D,
    right::CartesianBasisRepresentation3D,
)
    overlap_x = _basis_cross_overlap_1d(left.axis_representations.x, right.axis_representations.x)
    overlap_y = _basis_cross_overlap_1d(left.axis_representations.y, right.axis_representations.y)
    overlap_z = _basis_cross_overlap_1d(left.axis_representations.z, right.axis_representations.z)

    if left.contraction_kind == :identity && right.contraction_kind == :identity
        return Matrix{Float64}(kron(overlap_x, kron(overlap_y, overlap_z)))
    end

    left_parent = _cartesian_parent_state_basis(left)
    right_parent = _cartesian_parent_state_basis(right)
    raw_cross = _cartesian_product_cross_overlap(
        left_parent.states,
        right_parent.states,
        overlap_x,
        overlap_y,
        overlap_z,
    )

    if left_parent.coefficients !== nothing
        raw_cross = Matrix{Float64}(transpose(left_parent.coefficients) * raw_cross)
    end
    if right_parent.coefficients !== nothing
        raw_cross = Matrix{Float64}(raw_cross * right_parent.coefficients)
    end
    return raw_cross
end

function _cartesian_basis_supplement_axis_cross(
    basis::BasisRepresentation1D,
    supplement::CartesianGaussianShellSupplementRepresentation3D,
    axis::Symbol,
)
    isempty(supplement.orbitals) && return zeros(Float64, size(basis.coefficient_matrix, 2), 0)
    basis_primitives = collect(primitives(primitive_set(basis)))
    all(primitive -> primitive isa Gaussian, basis_primitives) || throw(
        ArgumentError(
            "exact hybrid cartesian-supplement overlap currently requires Gaussian 1D primitives on the Cartesian axis representation",
        ),
    )
    primitive_cross = zeros(Float64, length(basis_primitives), length(supplement.orbitals))
    for (column, orbital) in pairs(supplement.orbitals)
        center_value =
            axis == :x ? orbital.center[1] :
            axis == :y ? orbital.center[2] :
            axis == :z ? orbital.center[3] :
            throw(ArgumentError("axis must be :x, :y, or :z"))
        power =
            axis == :x ? orbital.angular_powers[1] :
            axis == :y ? orbital.angular_powers[2] :
            axis == :z ? orbital.angular_powers[3] :
            throw(ArgumentError("axis must be :x, :y, or :z"))
        for (row, primitive) in pairs(basis_primitives)
            alpha_basis = _qwrg_gaussian_exponent(primitive::Gaussian)
            value = 0.0
            for index in eachindex(orbital.exponents)
                exponent = Float64(orbital.exponents[index])
                coefficient = Float64(orbital.coefficients[index])
                prefactor = _qwrg_atomic_shell_prefactor(exponent, power)
                value += coefficient * _qwrg_atomic_basic_integral(
                    alpha_basis,
                    primitive.center_value,
                    0,
                    1.0,
                    exponent,
                    center_value,
                    power,
                    prefactor,
                )
            end
            primitive_cross[row, column] = value
        end
    end
    return Matrix{Float64}(transpose(basis.coefficient_matrix) * primitive_cross)
end

function _cartesian_basis_supplement_cross(
    basis::CartesianBasisRepresentation3D,
    supplement::CartesianGaussianShellSupplementRepresentation3D,
)
    isempty(supplement.orbitals) && return zeros(Float64, basis.metadata.final_dimension, 0)
    overlap_x =
        _cartesian_basis_supplement_axis_cross(basis.axis_representations.x, supplement, :x)
    overlap_y =
        _cartesian_basis_supplement_axis_cross(basis.axis_representations.y, supplement, :y)
    overlap_z =
        _cartesian_basis_supplement_axis_cross(basis.axis_representations.z, supplement, :z)
    parent = _cartesian_parent_state_basis(basis)
    raw_cross = zeros(Float64, length(parent.states), length(supplement.orbitals))
    for (row, (ix, iy, iz)) in pairs(parent.states)
        @inbounds for column in eachindex(supplement.orbitals)
            raw_cross[row, column] =
                overlap_x[ix, column] *
                overlap_y[iy, column] *
                overlap_z[iz, column]
        end
    end
    if parent.coefficients !== nothing
        return Matrix{Float64}(transpose(parent.coefficients) * raw_cross)
    end
    return raw_cross
end

function _cartesian_supplement_orbital_axis_overlap(
    left::CartesianGaussianShellOrbitalRepresentation3D,
    right::CartesianGaussianShellOrbitalRepresentation3D,
    axis::Symbol,
)
    left.primitive_normalization == :axiswise_normalized_cartesian_gaussian || throw(
        ArgumentError("unsupported left supplement primitive normalization :$(left.primitive_normalization)"),
    )
    right.primitive_normalization == :axiswise_normalized_cartesian_gaussian || throw(
        ArgumentError("unsupported right supplement primitive normalization :$(right.primitive_normalization)"),
    )
    center_left =
        axis == :x ? left.center[1] :
        axis == :y ? left.center[2] :
        axis == :z ? left.center[3] :
        throw(ArgumentError("axis must be :x, :y, or :z"))
    center_right =
        axis == :x ? right.center[1] :
        axis == :y ? right.center[2] :
        right.center[3]
    power_left =
        axis == :x ? left.angular_powers[1] :
        axis == :y ? left.angular_powers[2] :
        left.angular_powers[3]
    power_right =
        axis == :x ? right.angular_powers[1] :
        axis == :y ? right.angular_powers[2] :
        right.angular_powers[3]
    matrix = zeros(Float64, length(left.exponents), length(right.exponents))
    for j in eachindex(right.exponents)
        exponent_right = Float64(right.exponents[j])
        prefactor_right = _qwrg_atomic_shell_prefactor(exponent_right, power_right)
        for i in eachindex(left.exponents)
            exponent_left = Float64(left.exponents[i])
            prefactor_left = _qwrg_atomic_shell_prefactor(exponent_left, power_left)
            matrix[i, j] = _qwrg_atomic_basic_integral(
                exponent_left,
                center_left,
                power_left,
                prefactor_left,
                exponent_right,
                center_right,
                power_right,
                prefactor_right,
            )
        end
    end
    return Float64(dot(left.coefficients, matrix * right.coefficients))
end

function _cartesian_supplement_cross_overlap(
    left::CartesianGaussianShellSupplementRepresentation3D,
    right::CartesianGaussianShellSupplementRepresentation3D,
)
    matrix = zeros(Float64, length(left.orbitals), length(right.orbitals))
    for (row, left_orbital) in pairs(left.orbitals), (column, right_orbital) in pairs(right.orbitals)
        matrix[row, column] =
            _cartesian_supplement_orbital_axis_overlap(left_orbital, right_orbital, :x) *
            _cartesian_supplement_orbital_axis_overlap(left_orbital, right_orbital, :y) *
            _cartesian_supplement_orbital_axis_overlap(left_orbital, right_orbital, :z)
    end
    return matrix
end

function _cartesian_raw_components(
    representation::CartesianBasisRepresentation3D,
)
    if representation.metadata.parent_kind == :cartesian_product_basis
        nraw = representation.metadata.final_dimension
        return (
            cartesian_representation = representation,
            supplement_representation = _cartesian_empty_supplement_representation(),
            raw_to_final = Matrix{Float64}(I, nraw, nraw),
            factorized_cartesian_parent_basis = _cartesian_factorized_parent_basis(representation),
            cartesian_supplement_axis_tables = nothing,
        )
    elseif representation.metadata.parent_kind == :cartesian_plus_supplement_raw
        _cartesian_supports_exact_hybrid_overlap(representation) || throw(
            _cartesian_cross_overlap_error(representation, representation),
        )
        raw_to_final = representation.coefficient_matrix === nothing ?
            throw(
                ArgumentError(
                    "hybrid Cartesian basis representation requires an explicit raw_to_final coefficient matrix",
                ),
            ) :
            Matrix{Float64}(representation.coefficient_matrix)
        return (
            cartesian_representation = representation.parent_data.cartesian_parent_representation,
            supplement_representation = representation.parent_data.supplement_representation,
            raw_to_final = raw_to_final,
            factorized_cartesian_parent_basis =
                representation.parent_data.factorized_cartesian_parent_basis,
            cartesian_supplement_axis_tables =
                representation.parent_data.cartesian_supplement_axis_tables,
            exact_cartesian_supplement_overlap =
                hasproperty(representation.parent_data, :exact_cartesian_supplement_overlap) ?
                representation.parent_data.exact_cartesian_supplement_overlap : nothing,
            exact_supplement_overlap =
                hasproperty(representation.parent_data, :exact_supplement_overlap) ?
                representation.parent_data.exact_supplement_overlap : nothing,
        )
    end
    throw(_cartesian_cross_overlap_error(representation, representation))
end

function _cartesian_exact_cartesian_supplement_cross(
    left_raw,
    right_raw,
)
    hasproperty(left_raw, :exact_cartesian_supplement_overlap) || return nothing
    exact_cross = left_raw.exact_cartesian_supplement_overlap
    exact_cross === nothing && return nothing
    _cartesian_same_supplement_raw_identity(
        left_raw.supplement_representation,
        right_raw.supplement_representation,
    ) || return nothing
    return Matrix{Float64}(exact_cross)
end

function _cartesian_exact_supplement_cross_overlap(
    left_raw,
    right_raw,
)
    hasproperty(left_raw, :exact_supplement_overlap) || return nothing
    hasproperty(right_raw, :exact_supplement_overlap) || return nothing
    left_exact = left_raw.exact_supplement_overlap
    right_exact = right_raw.exact_supplement_overlap
    left_exact === nothing && return nothing
    right_exact === nothing && return nothing
    _cartesian_same_supplement_raw_identity(
        left_raw.supplement_representation,
        right_raw.supplement_representation,
    ) || return nothing
    return Matrix{Float64}(0.5 .* (Matrix{Float64}(left_exact) .+ Matrix{Float64}(right_exact)))
end

function _cartesian_hybrid_cartesian_supplement_cross_dense_reference(
    left_raw,
    right_raw,
)
    isempty(right_raw.supplement_representation.orbitals) &&
        return zeros(Float64, left_raw.cartesian_representation.metadata.final_dimension, 0)
    exact_cross = _cartesian_exact_cartesian_supplement_cross(left_raw, right_raw)
    exact_cross !== nothing && return exact_cross
    try
        return _cartesian_basis_supplement_cross(
            left_raw.cartesian_representation,
            right_raw.supplement_representation,
        )
    catch error
        if error isa ArgumentError
            return _cartesian_factorized_basis_supplement_cross(
                left_raw.cartesian_representation,
                right_raw.supplement_representation,
                left_raw,
                right_raw,
            )
        end
        rethrow()
    end
end

function _cartesian_mixed_raw_cross_overlap_dense_reference(
    left::CartesianBasisRepresentation3D,
    right::CartesianBasisRepresentation3D,
)
    left_raw = _cartesian_raw_components(left)
    right_raw = _cartesian_raw_components(right)

    cc = cross_overlap(
        left_raw.cartesian_representation,
        right_raw.cartesian_representation,
    )
    cg =
        isempty(right_raw.supplement_representation.orbitals) ?
        zeros(Float64, size(cc, 1), 0) :
        _cartesian_hybrid_cartesian_supplement_cross_dense_reference(left_raw, right_raw)
    gc =
        isempty(left_raw.supplement_representation.orbitals) ?
        zeros(Float64, 0, size(cc, 2)) :
        transpose(_cartesian_hybrid_cartesian_supplement_cross_dense_reference(right_raw, left_raw))
    gg =
        isempty(left_raw.supplement_representation.orbitals) ||
        isempty(right_raw.supplement_representation.orbitals) ?
        zeros(
            Float64,
            length(left_raw.supplement_representation.orbitals),
            length(right_raw.supplement_representation.orbitals),
        ) : begin
            exact_gg = _cartesian_exact_supplement_cross_overlap(left_raw, right_raw)
            exact_gg === nothing ?
            _cartesian_supplement_cross_overlap(
                left_raw.supplement_representation,
                right_raw.supplement_representation,
            ) :
            exact_gg
        end
    raw_cross = [cc cg; gc gg]
    return Matrix{Float64}(
        transpose(left_raw.raw_to_final) * raw_cross * right_raw.raw_to_final,
    )
end

function _cartesian_mixed_raw_cross_overlap_with_stage_timings(
    left::CartesianBasisRepresentation3D,
    right::CartesianBasisRepresentation3D,
)
    timer = _CartesianCrossOverlapStageTimer()
    left_raw = nothing
    right_raw = nothing
    _cartesian_time_stage!(timer, :factorized_parent_ns) do
        left_raw = _cartesian_raw_components(left)
        right_raw = _cartesian_raw_components(right)
        nothing
    end

    axis_cross_tables = _cartesian_time_stage!(timer, :axis_cross_table_ns) do
        _cartesian_factorized_axis_cross_tables(
            left_raw.cartesian_representation,
            right_raw.cartesian_representation,
            left_raw.factorized_cartesian_parent_basis,
            right_raw.factorized_cartesian_parent_basis,
        )
    end

    cc = _cartesian_time_stage!(timer, :cartesian_cartesian_block_ns) do
        _cartesian_factorized_cartesian_cross_overlap(
            left_raw.factorized_cartesian_parent_basis,
            right_raw.factorized_cartesian_parent_basis,
            axis_cross_tables.tables,
        )
    end

    cg = _cartesian_time_stage!(timer, :cartesian_supplement_block_ns) do
        isempty(right_raw.supplement_representation.orbitals) ? zeros(Float64, size(cc, 1), 0) :
        begin
            exact_cross = _cartesian_exact_cartesian_supplement_cross(left_raw, right_raw)
            exact_cross === nothing ?
            _cartesian_factorized_basis_supplement_cross(
                left_raw.factorized_cartesian_parent_basis,
                left_raw.cartesian_representation,
                right_raw.supplement_representation,
                left_raw,
                right_raw,
            ) :
            exact_cross
        end
    end

    gc = _cartesian_time_stage!(timer, :cartesian_supplement_block_ns) do
        isempty(left_raw.supplement_representation.orbitals) ? zeros(Float64, 0, size(cc, 2)) :
        begin
            exact_cross = _cartesian_exact_cartesian_supplement_cross(right_raw, left_raw)
            exact_cross === nothing ?
            transpose(
                _cartesian_factorized_basis_supplement_cross(
                    right_raw.factorized_cartesian_parent_basis,
                    right_raw.cartesian_representation,
                    left_raw.supplement_representation,
                    right_raw,
                    left_raw,
                ),
            ) :
            transpose(exact_cross)
        end
    end

    gg = _cartesian_time_stage!(timer, :supplement_supplement_block_ns) do
        isempty(left_raw.supplement_representation.orbitals) ||
        isempty(right_raw.supplement_representation.orbitals) ?
        zeros(
            Float64,
            length(left_raw.supplement_representation.orbitals),
            length(right_raw.supplement_representation.orbitals),
        ) : begin
            exact_gg = _cartesian_exact_supplement_cross_overlap(left_raw, right_raw)
            exact_gg === nothing ?
            _cartesian_supplement_cross_overlap(
                left_raw.supplement_representation,
                right_raw.supplement_representation,
            ) :
            exact_gg
        end
    end

    matrix = _cartesian_time_stage!(timer, :final_contraction_ns) do
        raw_cross = [cc cg; gc gg]
        Matrix{Float64}(transpose(left_raw.raw_to_final) * raw_cross * right_raw.raw_to_final)
    end

    return (
        matrix = matrix,
        stage_timings = _cartesian_cross_overlap_stage_timings(timer, axis_cross_tables.timings),
    )
end

function _cartesian_mixed_raw_cross_apply_with_stage_timings(
    left::CartesianBasisRepresentation3D,
    right::CartesianBasisRepresentation3D,
    right_coefficients::AbstractMatrix{<:Real},
)
    size(right_coefficients, 1) == right.metadata.final_dimension || throw(
        DimensionMismatch(
            "mixed raw transfer apply expected $(right.metadata.final_dimension) source final rows, got $(size(right_coefficients, 1))",
        ),
    )
    timer = _CartesianCrossOverlapStageTimer()
    left_raw = nothing
    right_raw = nothing
    _cartesian_time_stage!(timer, :factorized_parent_ns) do
        left_raw = _cartesian_raw_components(left)
        right_raw = _cartesian_raw_components(right)
        nothing
    end

    axis_cross_tables = _cartesian_time_stage!(timer, :axis_cross_table_ns) do
        _cartesian_factorized_axis_cross_tables(
            left_raw.cartesian_representation,
            right_raw.cartesian_representation,
            left_raw.factorized_cartesian_parent_basis,
            right_raw.factorized_cartesian_parent_basis,
        )
    end

    right_raw_coefficients = _cartesian_time_stage!(timer, :final_contraction_ns) do
        Matrix{Float64}(right_raw.raw_to_final * right_coefficients)
    end
    right_cartesian_dimension = right_raw.cartesian_representation.metadata.final_dimension
    right_supplement_dimension = length(right_raw.supplement_representation.orbitals)
    right_cartesian_coefficients =
        @view right_raw_coefficients[1:right_cartesian_dimension, :]
    right_supplement_coefficients =
        @view right_raw_coefficients[(right_cartesian_dimension + 1):(right_cartesian_dimension + right_supplement_dimension), :]

    target_cartesian = _cartesian_time_stage!(timer, :cartesian_cartesian_block_ns) do
        _cartesian_factorized_cartesian_cross_apply(
            left_raw.factorized_cartesian_parent_basis,
            right_raw.factorized_cartesian_parent_basis,
            axis_cross_tables.tables,
            right_cartesian_coefficients,
        )
    end

    _cartesian_time_stage!(timer, :cartesian_supplement_block_ns) do
        if right_supplement_dimension > 0
            cg = _cartesian_exact_cartesian_supplement_cross(left_raw, right_raw)
            if cg === nothing
                cg = _cartesian_factorized_basis_supplement_cross(
                    left_raw.factorized_cartesian_parent_basis,
                    left_raw.cartesian_representation,
                    right_raw.supplement_representation,
                    left_raw,
                    right_raw,
                )
            end
            target_cartesian .+= cg * right_supplement_coefficients
        end
        nothing
    end

    left_supplement_dimension = length(left_raw.supplement_representation.orbitals)
    target_supplement = _cartesian_time_stage!(timer, :cartesian_supplement_block_ns) do
        if left_supplement_dimension == 0
            zeros(Float64, 0, size(right_coefficients, 2))
        else
            gc = _cartesian_exact_cartesian_supplement_cross(right_raw, left_raw)
            if gc === nothing
                gc = _cartesian_factorized_basis_supplement_cross(
                    right_raw.factorized_cartesian_parent_basis,
                    right_raw.cartesian_representation,
                    left_raw.supplement_representation,
                    right_raw,
                    left_raw,
                )
            end
            Matrix{Float64}(transpose(gc) * right_cartesian_coefficients)
        end
    end

    _cartesian_time_stage!(timer, :supplement_supplement_block_ns) do
        if left_supplement_dimension > 0 && right_supplement_dimension > 0
            exact_gg = _cartesian_exact_supplement_cross_overlap(left_raw, right_raw)
            gg =
                exact_gg === nothing ?
                _cartesian_supplement_cross_overlap(
                    left_raw.supplement_representation,
                    right_raw.supplement_representation,
                ) :
                exact_gg
            target_supplement .+= gg * right_supplement_coefficients
        end
        nothing
    end

    coefficients = _cartesian_time_stage!(timer, :final_contraction_ns) do
        raw_target_coefficients =
            left_supplement_dimension == 0 ?
            target_cartesian :
            [target_cartesian; target_supplement]
        Matrix{Float64}(transpose(left_raw.raw_to_final) * raw_target_coefficients)
    end

    return (
        coefficients = coefficients,
        stage_timings = _cartesian_cross_overlap_stage_timings(timer, axis_cross_tables.timings),
    )
end

function _cartesian_mixed_raw_cross_overlap(
    left::CartesianBasisRepresentation3D,
    right::CartesianBasisRepresentation3D,
)
    return _cartesian_mixed_raw_cross_overlap_with_stage_timings(left, right).matrix
end

function _cartesian_uses_exact_mixed_raw_cross_overlap(
    left::CartesianBasisRepresentation3D,
    right::CartesianBasisRepresentation3D,
)
    return (
        left.metadata.parent_kind == :cartesian_plus_supplement_raw ||
        right.metadata.parent_kind == :cartesian_plus_supplement_raw
    ) &&
           _cartesian_supports_exact_hybrid_overlap(left) &&
           _cartesian_supports_exact_hybrid_overlap(right)
end

"""
    cross_overlap(left::CartesianBasisRepresentation3D, right::CartesianBasisRepresentation3D)

Return the exact Cartesian basis cross-overlap matrix between `left` and
`right`, with rows ordered by the left basis and columns ordered by the right
basis.

This first pass supports two exact representation classes:

- direct-product Cartesian representations
- nested fixed-block Cartesian representations
- hybrid residual-Gaussian final bases whose public representation includes the
  explicit Cartesian parent representation and supplement-orbital
  representation needed to build the exact mixed raw-space overlap

The implementation is algebraic:

- if both sides expose an explicit Cartesian product parent, build the parent
  cross overlap from the stored 1D representation layers
- then contract by the stored 3D coefficient matrices
- if either side is hybrid residual-Gaussian, build the exact mixed raw-space
  block from:
  - Cartesian-Cartesian overlap on the public Cartesian parent representations
  - Cartesian-supplement overlap against the public supplement-orbital
    representation
  - supplement-supplement overlap on the public supplement representation
  - then contract by each side's `raw_to_final`
"""
function cross_overlap(
    left::CartesianBasisRepresentation3D,
    right::CartesianBasisRepresentation3D,
)
    if left.metadata.parent_kind == :cartesian_product_basis &&
       right.metadata.parent_kind == :cartesian_product_basis
        return _cartesian_parent_cross_overlap(left, right)
    elseif _cartesian_uses_exact_mixed_raw_cross_overlap(left, right)
        return _cartesian_mixed_raw_cross_overlap(left, right)
    end
    throw(_cartesian_cross_overlap_error(left, right))
end

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
    (
        representation.metadata.parent_kind == :cartesian_product_basis ||
        _cartesian_supports_exact_hybrid_overlap(representation)
    ) || throw(
        ArgumentError(
            "gto_overlap_matrix supports explicit Cartesian product, nested fixed-block, and exact hybrid residual-Gaussian Cartesian working bases; got parent_kind :$(representation.metadata.parent_kind)",
        ),
    )
    return representation
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
            return _cartesian_factorized_basis_supplement_cross(
                raw.factorized_cartesian_parent_basis,
                raw.cartesian_representation,
                probes,
                raw,
            )
        end
        rethrow()
    end
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

function _gto_overlap_matrix(
    working::CartesianBasisRepresentation3D,
    probes::CartesianGaussianShellSupplementRepresentation3D,
)
    raw = _cartesian_raw_components(working)
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

"""
    CartesianBasisTransferDiagnostics

Small public diagnostic bundle for exact Cartesian basis-to-basis transfer on
the supported Cartesian representation families.

The diagnostics record:

- the source and target basis dimensions
- whether the exact transfer used the same-parent pure-Cartesian,
  general pure-Cartesian, or hybrid mixed-raw path
- the defect of the materialized transfer operator relative to the exact final
  basis cross overlap used by the working transfer contract
- optional self-overlap diagnostics when that data was explicitly constructed
- optional source/target Euclidean norm traces and transfer residuals when
  specific orbital coefficients were transferred
"""
struct CartesianBasisTransferDiagnostics
    source_dimension::Int
    target_dimension::Int
    transfer_path::Symbol
    projector_residual_inf::Float64
    source_self_overlap_error::Float64
    target_self_overlap_error::Float64
    source_metric_trace::Union{Nothing,Float64}
    target_metric_trace::Union{Nothing,Float64}
    transferred_residual_inf::Union{Nothing,Float64}
end

function Base.show(io::IO, diagnostics::CartesianBasisTransferDiagnostics)
    print(
        io,
        "CartesianBasisTransferDiagnostics(source_dim=",
        diagnostics.source_dimension,
        ", target_dim=",
        diagnostics.target_dimension,
        ", path=:",
        diagnostics.transfer_path,
        ", projector_residual_inf=",
        diagnostics.projector_residual_inf,
        ")",
    )
end

"""
    CartesianBasisProjector3D

Exact metric-consistent Cartesian basis transfer matrix together with light
diagnostics describing how it was built.
"""
struct CartesianBasisProjector3D
    matrix::Matrix{Float64}
    diagnostics::CartesianBasisTransferDiagnostics
end

function Base.show(io::IO, projector::CartesianBasisProjector3D)
    print(
        io,
        "CartesianBasisProjector3D(size=",
        size(projector.matrix),
        ", path=:",
        projector.diagnostics.transfer_path,
        ")",
    )
end

"""
    CartesianOrbitalTransferResult

Returned by `transfer_orbitals(...)`. Bundles the transferred coefficients,
the exact basis projector when one was materialized, and transfer diagnostics.
"""
struct CartesianOrbitalTransferResult{CT <: Union{Vector{Float64},Matrix{Float64}}}
    coefficients::CT
    projector::Union{Nothing,CartesianBasisProjector3D}
    diagnostics::CartesianBasisTransferDiagnostics
end

function Base.show(io::IO, result::CartesianOrbitalTransferResult)
    print(
        io,
        "CartesianOrbitalTransferResult(size=",
        size(result.coefficients),
        ", path=:",
        result.diagnostics.transfer_path,
        ", projector_size=",
        result.projector === nothing ? "none" : string(size(result.projector.matrix)),
        ")",
    )
end

function _cartesian_transfer_path(
    source::CartesianBasisRepresentation3D,
    target::CartesianBasisRepresentation3D,
)
    if source.metadata.parent_kind == :cartesian_plus_supplement_raw ||
       target.metadata.parent_kind == :cartesian_plus_supplement_raw
        return :hybrid_mixed_raw_cross_overlap_transfer
    end
    return _cartesian_same_parent_raw_identity(source, target) ?
           :same_parent_cross_overlap_transfer :
           :pure_cartesian_cross_overlap_transfer
end

function _cartesian_columns_are_nearly_orthonormal(
    coefficients::AbstractMatrix{<:Real};
    atol::Float64 = 1.0e-8,
)
    size(coefficients, 2) > 1 || return false
    gram = Matrix{Float64}(transpose(coefficients) * coefficients)
    identity = Matrix{Float64}(LinearAlgebra.I, size(gram, 1), size(gram, 2))
    return norm(gram - identity, Inf) <= atol
end

function _cartesian_euclidean_orthonormalize_columns(
    coefficients::AbstractMatrix{<:Real},
)
    ncols = size(coefficients, 2)
    factorization = qr(Matrix{Float64}(coefficients))
    orthonormal = Matrix{Float64}(factorization.Q[:, 1:ncols])
    rblock = Matrix{Float64}(factorization.R[1:ncols, 1:ncols])
    for column in 1:ncols
        if rblock[column, column] < 0.0
            orthonormal[:, column] .*= -1.0
        end
    end
    return orthonormal
end

function _cartesian_finalize_transferred_coefficients(
    source_coefficients::AbstractVector{<:Real},
    transferred_coefficients_matrix::AbstractMatrix{<:Real},
)
    return (
        vec(Matrix{Float64}(transferred_coefficients_matrix)),
        false,
    )
end

function _cartesian_finalize_transferred_coefficients(
    source_coefficients::AbstractMatrix{<:Real},
    transferred_coefficients_matrix::AbstractMatrix{<:Real},
)
    source_coefficients_matrix = _cartesian_coefficients_matrix(source_coefficients)
    transferred_matrix = Matrix{Float64}(transferred_coefficients_matrix)
    if _cartesian_columns_are_nearly_orthonormal(source_coefficients_matrix) &&
       !_cartesian_columns_are_nearly_orthonormal(transferred_matrix)
        return (
            _cartesian_euclidean_orthonormalize_columns(transferred_matrix),
            true,
        )
    end
    return (transferred_matrix, false)
end

function _cartesian_transfer_diagnostics(
    source::CartesianBasisRepresentation3D,
    target::CartesianBasisRepresentation3D,
    target_source_overlap::AbstractMatrix{<:Real},
    projector::AbstractMatrix{<:Real};
    source_coefficients::Union{Nothing,AbstractVecOrMat{<:Real}} = nothing,
    transferred_coefficients::Union{Nothing,AbstractVecOrMat{<:Real}} = nothing,
    orthonormalized::Bool = false,
)
    source_dimension = source.metadata.final_dimension
    target_dimension = target.metadata.final_dimension
    transfer_path = _cartesian_transfer_path(source, target)
    projector_residual_inf = norm(projector - target_source_overlap, Inf)

    if source_coefficients === nothing || transferred_coefficients === nothing
        return CartesianBasisTransferDiagnostics(
            source_dimension,
            target_dimension,
            transfer_path,
            projector_residual_inf,
            NaN,
            NaN,
            nothing,
            nothing,
            nothing,
        )
    end

    source_coefficients_matrix = _cartesian_coefficients_matrix(source_coefficients)
    transferred_coefficients_matrix = _cartesian_coefficients_matrix(transferred_coefficients)
    source_metric =
        Matrix{Float64}(transpose(source_coefficients_matrix) * source_coefficients_matrix)
    target_metric =
        Matrix{Float64}(transpose(transferred_coefficients_matrix) * transferred_coefficients_matrix)
    transferred_residual_inf =
        orthonormalized ? NaN : norm(
            transferred_coefficients_matrix - projector * source_coefficients_matrix,
            Inf,
        )

    return CartesianBasisTransferDiagnostics(
        source_dimension,
        target_dimension,
        transfer_path,
        projector_residual_inf,
        NaN,
        NaN,
        tr(source_metric),
        tr(target_metric),
        transferred_residual_inf,
    )
end

function _cartesian_coefficients_matrix(source_coefficients::AbstractVector{<:Real})
    return reshape(Float64[Float64(value) for value in source_coefficients], :, 1)
end

function _cartesian_coefficients_matrix(source_coefficients::AbstractMatrix{<:Real})
    return Matrix{Float64}(source_coefficients)
end

function _cartesian_transfer_diagnostics_without_projector(
    source::CartesianBasisRepresentation3D,
    target::CartesianBasisRepresentation3D,
    source_coefficients::AbstractVecOrMat{<:Real},
    transferred_coefficients::AbstractVecOrMat{<:Real},
)
    source_coefficients_matrix = _cartesian_coefficients_matrix(source_coefficients)
    transferred_coefficients_matrix = _cartesian_coefficients_matrix(transferred_coefficients)
    return CartesianBasisTransferDiagnostics(
        source.metadata.final_dimension,
        target.metadata.final_dimension,
        _cartesian_transfer_path(source, target),
        NaN,
        NaN,
        NaN,
        tr(transpose(source_coefficients_matrix) * source_coefficients_matrix),
        tr(transpose(transferred_coefficients_matrix) * transferred_coefficients_matrix),
        NaN,
    )
end

function _cartesian_basis_projector_with_stage_timings(
    source::CartesianBasisRepresentation3D,
    target::CartesianBasisRepresentation3D,
)
    overlap_result =
        _cartesian_uses_exact_mixed_raw_cross_overlap(target, source) ?
        _cartesian_mixed_raw_cross_overlap_with_stage_timings(target, source) :
        (matrix = Matrix{Float64}(cross_overlap(target, source)), stage_timings = nothing)
    diagnostics = _cartesian_transfer_diagnostics(
        source,
        target,
        overlap_result.matrix,
        overlap_result.matrix,
    )
    return (
        projector = CartesianBasisProjector3D(overlap_result.matrix, diagnostics),
        stage_timings = overlap_result.stage_timings,
    )
end

"""
    basis_projector(
        source::CartesianBasisRepresentation3D,
        target::CartesianBasisRepresentation3D,
    )

Return the exact metric-consistent basis-to-basis transfer operator from
`source` into `target` for the currently supported Cartesian representations.

For the final working-basis contract used in this repo, the transfer matrix is
the exact final-basis cross overlap

`T_{B <- A} = S_{BA} = <B|A>`.

This first pass supports the same representation families as `cross_overlap(...)`,
including hybrid residual-Gaussian final bases whose public representation
exposes the exact mixed raw-space identity.

Final-basis self-overlaps are construction/debug diagnostics only. Normal
transfer on these orthonormal working bases uses only the exact cross overlap;
raw/nonorthogonal generalized-overlap algebra remains outside this helper.
"""
function basis_projector(
    source::CartesianBasisRepresentation3D,
    target::CartesianBasisRepresentation3D,
)
    return _cartesian_basis_projector_with_stage_timings(source, target).projector
end

function _cartesian_transfer_coefficients(
    source_coefficients::AbstractVector{<:Real},
    projector::CartesianBasisProjector3D,
)
    source_dimension = size(projector.matrix, 2)
    length(source_coefficients) == source_dimension || throw(
        DimensionMismatch(
            "source coefficient vector length $(length(source_coefficients)) does not match source basis dimension $source_dimension",
        ),
    )
    return Vector{Float64}(projector.matrix * Vector{Float64}(source_coefficients))
end

function _cartesian_transfer_coefficients(
    source_coefficients::AbstractMatrix{<:Real},
    projector::CartesianBasisProjector3D,
)
    source_dimension = size(projector.matrix, 2)
    size(source_coefficients, 1) == source_dimension || throw(
        DimensionMismatch(
            "source coefficient matrix row count $(size(source_coefficients, 1)) does not match source basis dimension $source_dimension",
        ),
    )
    return Matrix{Float64}(projector.matrix * Matrix{Float64}(source_coefficients))
end

function transfer_orbitals(
    source_coefficients::AbstractVecOrMat{<:Real},
    projector::CartesianBasisProjector3D,
)
    transferred_coefficients_matrix = _cartesian_transfer_coefficients(source_coefficients, projector)
    transferred_coefficients, orthonormalized =
        _cartesian_finalize_transferred_coefficients(
            source_coefficients,
            _cartesian_coefficients_matrix(transferred_coefficients_matrix),
        )
    source_coefficients_matrix = _cartesian_coefficients_matrix(source_coefficients)
    transferred_coefficients_matrix = _cartesian_coefficients_matrix(transferred_coefficients)
    diagnostics = CartesianBasisTransferDiagnostics(
        projector.diagnostics.source_dimension,
        projector.diagnostics.target_dimension,
        projector.diagnostics.transfer_path,
        projector.diagnostics.projector_residual_inf,
        projector.diagnostics.source_self_overlap_error,
        projector.diagnostics.target_self_overlap_error,
        tr(transpose(source_coefficients_matrix) * source_coefficients_matrix),
        tr(transpose(transferred_coefficients_matrix) * transferred_coefficients_matrix),
        orthonormalized ? NaN : norm(
            transferred_coefficients_matrix -
            projector.matrix * source_coefficients_matrix,
            Inf,
        ),
    )
    return CartesianOrbitalTransferResult(
        transferred_coefficients,
        projector,
        diagnostics,
    )
end

"""
    transfer_orbitals(
        source_coefficients::AbstractVecOrMat,
        source::CartesianBasisRepresentation3D,
        target::CartesianBasisRepresentation3D,
    )

Transfer orbital coefficients from `source` into `target` using the exact final
basis cross overlap returned by `basis_projector(...)`. Pass
`materialize_projector = false` on supported mixed-raw hybrid transfers to apply
the same transfer to the supplied coefficient block without building the full
dense projector.

For the final orthonormal working-basis contract used in this repo, the transfer
formula is

`C_B = S_{BA} * C_A`.

Final-basis self-overlaps are diagnostic only and are not part of this normal
transfer path. When the source columns already represent an orthonormal occupied
block, the transferred target block may be Euclidean-orthonormalized as one
final cleanup step.
"""
function transfer_orbitals(
    source_coefficients::AbstractVecOrMat{<:Real},
    source::CartesianBasisRepresentation3D,
    target::CartesianBasisRepresentation3D;
    materialize_projector::Bool = true,
)
    if !materialize_projector
        _cartesian_uses_exact_mixed_raw_cross_overlap(target, source) || throw(
            ArgumentError(
                "materialize_projector=false currently requires an exact mixed-raw Cartesian transfer path",
            ),
        )
        source_coefficients_matrix = _cartesian_coefficients_matrix(source_coefficients)
        transfer_result = _cartesian_mixed_raw_cross_apply_with_stage_timings(
            target,
            source,
            source_coefficients_matrix,
        )
        transferred_coefficients, _ =
            _cartesian_finalize_transferred_coefficients(
                source_coefficients,
                transfer_result.coefficients,
            )
        return CartesianOrbitalTransferResult(
            transferred_coefficients,
            nothing,
            _cartesian_transfer_diagnostics_without_projector(
                source,
                target,
                source_coefficients,
                transferred_coefficients,
            ),
        )
    end
    projector = basis_projector(source, target)
    return transfer_orbitals(source_coefficients, projector)
end

function _cartesian_basis_labels(prefix::AbstractString, count::Int)
    return [string(prefix, index) for index in 1:count]
end

function _cartesian_product_labels(
    orbitals::AbstractVector{<:CartesianProductOrbital3D},
)
    return [string("g(", orbital.ix, ",", orbital.iy, ",", orbital.iz, ")") for orbital in orbitals]
end

function _cartesian_labels(orbitals::AbstractVector{<:OrdinaryCartesianOrbital3D})
    return String[String(orbital.label) for orbital in orbitals]
end

function _cartesian_labels(orbitals::AbstractVector{<:_AtomicCartesianShellOrbital3D})
    return String[String(orbital.label) for orbital in orbitals]
end

function _cartesian_center_matrix(
    orbitals::AbstractVector{<:CartesianProductOrbital3D},
)
    centers_out = Matrix{Float64}(undef, length(orbitals), 3)
    for (row, orbital) in pairs(orbitals)
        centers_out[row, 1] = orbital.x
        centers_out[row, 2] = orbital.y
        centers_out[row, 3] = orbital.z
    end
    return centers_out
end

function _cartesian_center_matrix(
    orbitals::AbstractVector{<:OrdinaryCartesianOrbital3D},
)
    centers_out = Matrix{Float64}(undef, length(orbitals), 3)
    for (row, orbital) in pairs(orbitals)
        centers_out[row, 1] = orbital.x
        centers_out[row, 2] = orbital.y
        centers_out[row, 3] = orbital.z
    end
    return centers_out
end

function _cartesian_center_matrix(
    orbitals::AbstractVector{<:_AtomicCartesianShellOrbital3D},
)
    centers_out = Matrix{Float64}(undef, length(orbitals), 3)
    for (row, orbital) in pairs(orbitals)
        centers_out[row, 1] = orbital.center[1]
        centers_out[row, 2] = orbital.center[2]
        centers_out[row, 3] = orbital.center[3]
    end
    return centers_out
end

function _cartesian_axis_representation(basis::MappedUniformBasis)
    # The mapped ordinary Cartesian lane should expose the post-proxy,
    # post-cleanup working layer rather than the upstream distorted primitives.
    return _mapped_ordinary_working_basis_representation(basis)
end

function _cartesian_axis_representation(basis)
    return basis_representation(basis)
end

function _cartesian_axis_representations(basis::MappedUniformBasis)
    representation = _cartesian_axis_representation(basis)
    return (x = representation, y = representation, z = representation)
end

function _cartesian_axis_representations(basis::AbstractBondAlignedOrdinaryQWBasis3D)
    return (
        x = _cartesian_axis_representation(basis.basis_x),
        y = _cartesian_axis_representation(basis.basis_y),
        z = _cartesian_axis_representation(basis.basis_z),
    )
end

function _cartesian_axis_metadata(axis_representations::NamedTuple)
    return (
        x = axis_representations.x.metadata,
        y = axis_representations.y.metadata,
        z = axis_representations.z.metadata,
    )
end

function _cartesian_axis_signature(representation::BasisRepresentation1D)
    metadata = representation.metadata
    return (
        metadata.basis_kind,
        metadata.family_name,
        metadata.center_data,
        metadata.reference_center_data,
        metadata.integral_weight_data,
        metadata.basis_labels,
        metadata.coefficient_matrix,
    )
end

function _cartesian_axis_sharing(axis_representations::NamedTuple)
    x_signature = _cartesian_axis_signature(axis_representations.x)
    y_signature = _cartesian_axis_signature(axis_representations.y)
    z_signature = _cartesian_axis_signature(axis_representations.z)
    xy = isequal(x_signature, y_signature)
    xz = isequal(x_signature, z_signature)
    yz = isequal(y_signature, z_signature)
    xy && yz && return :shared_xyz
    xy && return :shared_xy
    xz && return :shared_xz
    yz && return :shared_yz
    return :separate_axes
end

function _cartesian_axis_counts(axis_representations::NamedTuple)
    return (
        length(axis_representations.x.metadata.center_data),
        length(axis_representations.y.metadata.center_data),
        length(axis_representations.z.metadata.center_data),
    )
end

function _cartesian_axis_centers(axis_representations::NamedTuple)
    return (
        axis_representations.x.metadata.center_data,
        axis_representations.y.metadata.center_data,
        axis_representations.z.metadata.center_data,
    )
end

function _cartesian_direct_product_representation(
    axis_representations::NamedTuple;
    route_metadata::NamedTuple = (;),
)
    axis_metadata = _cartesian_axis_metadata(axis_representations)
    axis_counts = _cartesian_axis_counts(axis_representations)
    x_centers, y_centers, z_centers = _cartesian_axis_centers(axis_representations)
    orbitals = _mapped_cartesian_orbitals(x_centers, y_centers, z_centers)
    labels = _cartesian_product_labels(orbitals)
    centers_matrix = _cartesian_center_matrix(orbitals)
    metadata = CartesianBasisMetadata3D(
        :direct_product,
        _cartesian_axis_sharing(axis_representations),
        axis_metadata,
        :cartesian_product_basis,
        axis_counts,
        length(orbitals),
        length(orbitals),
        nothing,
        labels,
        centers_matrix,
        route_metadata,
    )
    return CartesianBasisRepresentation3D(
        metadata,
        axis_representations,
        :identity,
        nothing,
        copy(labels),
        copy(centers_matrix),
        nothing,
        nothing,
        (;),
    )
end

function _cartesian_direct_product_representation(
    basis::MappedUniformBasis;
    route_metadata::NamedTuple = (;),
)
    return _cartesian_direct_product_representation(
        _cartesian_axis_representations(basis);
        route_metadata = route_metadata,
    )
end

function _cartesian_direct_product_representation(
    basis::AbstractBondAlignedOrdinaryQWBasis3D;
    route_metadata::NamedTuple = (;),
)
    return _cartesian_direct_product_representation(
        _cartesian_axis_representations(basis);
        route_metadata = route_metadata,
    )
end

function _cartesian_basis_route_metadata(basis::BondAlignedDiatomicQWBasis3D)
    return (
        basis_family = :bond_aligned_diatomic,
        bond_axis = basis.bond_axis,
        nuclei = copy(basis.nuclei),
        target_core_spacing = basis.target_core_spacing,
    )
end

function _cartesian_basis_route_metadata(basis::BondAlignedHomonuclearChainQWBasis3D)
    return (
        basis_family = :bond_aligned_homonuclear_chain,
        chain_axis = basis.chain_axis,
        chain_coordinates = Float64[Float64(value) for value in basis.chain_coordinates],
        nuclei = copy(basis.nuclei),
        target_core_spacing = basis.target_core_spacing,
    )
end

function _cartesian_basis_route_metadata(basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D)
    return (
        basis_family = :axis_aligned_homonuclear_square_lattice,
        lattice_size = basis.lattice_size,
        x_coordinates = Float64[Float64(value) for value in basis.x_coordinates],
        y_coordinates = Float64[Float64(value) for value in basis.y_coordinates],
        nuclei = copy(basis.nuclei),
        target_core_spacing = basis.target_core_spacing,
    )
end

function basis_representation(basis::BondAlignedDiatomicQWBasis3D)
    return _cartesian_direct_product_representation(
        basis;
        route_metadata = _cartesian_basis_route_metadata(basis),
    )
end

function basis_representation(basis::BondAlignedHomonuclearChainQWBasis3D)
    return _cartesian_direct_product_representation(
        basis;
        route_metadata = _cartesian_basis_route_metadata(basis),
    )
end

function basis_representation(basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D)
    return _cartesian_direct_product_representation(
        basis;
        route_metadata = _cartesian_basis_route_metadata(basis),
    )
end

basis_metadata(basis::AbstractBondAlignedOrdinaryQWBasis3D) = basis_representation(basis).metadata

function _cartesian_bond_aligned_build_metadata(
    basis::AbstractBondAlignedOrdinaryQWBasis3D,
)
    representation = basis_representation(basis)
    route_metadata = representation.metadata.route_metadata
    return (
        carried_representation = representation,
        parent_representation = representation,
        carried_metadata = representation.metadata,
        parent_metadata = representation.metadata,
        carried_route_metadata = route_metadata,
        parent_route_metadata = route_metadata,
        basis_family = route_metadata.basis_family,
    )
end

function _cartesian_parent_axis_representations(parent_basis)
    parent_basis isa MappedUniformBasis && return _cartesian_axis_representations(parent_basis)
    parent_basis isa AbstractBondAlignedOrdinaryQWBasis3D &&
        return _cartesian_axis_representations(parent_basis)
    throw(
        ArgumentError(
            "Cartesian basis representation does not yet support parent basis $(typeof(parent_basis))",
        ),
    )
end

function _cartesian_full_parent_box(axis_counts::NTuple{3,Int})
    return (1:axis_counts[1], 1:axis_counts[2], 1:axis_counts[3])
end

function _cartesian_shell_kind(::Any)
    return :unknown_shell
end

function _cartesian_shell_kind(::_CartesianNestedShell3D)
    return :shell
end

function _cartesian_shell_kind(::_CartesianNestedCompleteShell3D)
    return :complete_shell
end

function _cartesian_shell_kind(::_CartesianNestedShellPlusCore3D)
    return :shell_plus_core
end

function _cartesian_shell_kind(::_CartesianNestedShellSequence3D)
    return :shell_sequence
end

function _cartesian_shell_working_box(shell)
    return hasproperty(shell, :working_box) ? shell.working_box : nothing
end

function _cartesian_working_box_profile(
    working_box::Union{Nothing,NTuple{3,UnitRange{Int}}},
    axis_counts::NTuple{3,Int},
)
    working_box === nothing && return nothing
    return working_box == _cartesian_full_parent_box(axis_counts) ? :full_parent : :explicit_inner_box
end

function _cartesian_maybe_complete_shell_nside(shell)
    shell isa _CartesianNestedShellSequence3D || return nothing
    working_box = shell.working_box
    length(working_box[1]) == length(working_box[2]) == length(working_box[3]) || return nothing
    core_side = round(Int, cbrt(length(shell.core_indices)))
    core_side^3 == length(shell.core_indices) || return nothing
    expected_increment = _one_center_atomic_shell_increment(core_side)
    all(length(range) == expected_increment for range in shell.layer_column_ranges) || return nothing
    return core_side
end

function _cartesian_support_states(
    fixed_block::_NestedFixedBlock3D,
    axis_counts::NTuple{3,Int},
)
    if hasproperty(fixed_block.shell, :support_states)
        return NTuple{3,Int}[state for state in fixed_block.shell.support_states]
    end
    return NTuple{3,Int}[
        _cartesian_unflat_index(index, axis_counts) for index in fixed_block.support_indices
    ]
end

function _cartesian_fixed_block_route_metadata(
    fixed_block::_NestedFixedBlock3D,
    axis_counts::NTuple{3,Int},
)
    working_box = _cartesian_shell_working_box(fixed_block.shell)
    return (
        shell_kind = _cartesian_shell_kind(fixed_block.shell),
        working_box_profile = _cartesian_working_box_profile(working_box, axis_counts),
        nside = _cartesian_maybe_complete_shell_nside(fixed_block.shell),
        support_count = length(fixed_block.support_indices),
        term_storage = fixed_block.term_storage,
    )
end

function basis_representation(fixed_block::_NestedFixedBlock3D)
    axis_representations = _cartesian_parent_axis_representations(fixed_block.parent_basis)
    axis_metadata = _cartesian_axis_metadata(axis_representations)
    axis_counts = _cartesian_axis_counts(axis_representations)
    x_centers, y_centers, z_centers = _cartesian_axis_centers(axis_representations)
    parent_orbitals = _mapped_cartesian_orbitals(x_centers, y_centers, z_centers)
    parent_labels = _cartesian_product_labels(parent_orbitals)
    parent_centers = _cartesian_center_matrix(parent_orbitals)
    support_states = _cartesian_support_states(fixed_block, axis_counts)
    final_labels = _cartesian_basis_labels("nf", size(fixed_block.coefficient_matrix, 2))
    metadata = CartesianBasisMetadata3D(
        :nested_fixed_block,
        _cartesian_axis_sharing(axis_representations),
        axis_metadata,
        :cartesian_product_basis,
        axis_counts,
        length(parent_orbitals),
        size(fixed_block.coefficient_matrix, 2),
        _cartesian_shell_working_box(fixed_block.shell),
        final_labels,
        Matrix{Float64}(fixed_block.fixed_centers),
        _cartesian_fixed_block_route_metadata(fixed_block, axis_counts),
    )
    parent_data =
        isnothing(fixed_block.factorized_cartesian_parent_basis[]) ?
        (;) :
        (;
            factorized_cartesian_parent_basis =
                fixed_block.factorized_cartesian_parent_basis[],
        )
    return CartesianBasisRepresentation3D(
        metadata,
        axis_representations,
        :dense,
        _cartesian_coefficient_map_storage(fixed_block.coefficient_matrix),
        parent_labels,
        parent_centers,
        Vector{Int}(fixed_block.support_indices),
        support_states,
        parent_data,
    )
end

basis_metadata(fixed_block::_NestedFixedBlock3D) = basis_representation(fixed_block).metadata

function _cartesian_bond_aligned_build_metadata(
    fixed_block::_NestedFixedBlock3D{<:AbstractBondAlignedOrdinaryQWBasis3D},
)
    carried_representation = basis_representation(fixed_block)
    parent_representation = basis_representation(fixed_block.parent_basis)
    return (
        carried_representation = carried_representation,
        parent_representation = parent_representation,
        carried_metadata = carried_representation.metadata,
        parent_metadata = parent_representation.metadata,
        carried_route_metadata = carried_representation.metadata.route_metadata,
        parent_route_metadata = parent_representation.metadata.route_metadata,
        basis_family = parent_representation.metadata.route_metadata.basis_family,
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

function _cartesian_hybrid_parent_basis(
    operators::OrdinaryCartesianOperators3D,
)
    return operators.basis isa _NestedFixedBlock3D ? operators.basis.parent_basis : operators.basis
end

function _cartesian_atomic_axis_bundles(
    operators::OrdinaryCartesianOperators3D,
    parent_basis::MappedUniformBasis,
)
    gausslet_bundle = _mapped_ordinary_gausslet_1d_bundle(
        parent_basis;
        exponents = operators.expansion.exponents,
        center = 0.0,
        backend = operators.gausslet_backend,
    )
    return (x = gausslet_bundle, y = gausslet_bundle, z = gausslet_bundle)
end

function _cartesian_hybrid_supplement_axis_tables(
    operators::OrdinaryCartesianOperators3D,
    factorized_cartesian_parent_basis::_CartesianNestedFactorizedBasis3D,
    supplement3d,
    axis_bundles::NamedTuple{(:x, :y, :z)},
    route_label::AbstractString,
)
    proxy_x = axis_bundles.x.pgdg_intermediate.auxiliary_layer
    proxy_y = axis_bundles.y.pgdg_intermediate.auxiliary_layer
    proxy_z = axis_bundles.z.pgdg_intermediate.auxiliary_layer
    proxy_x isa _MappedLegacyProxyLayer1D || throw(
        ArgumentError(
            "$(route_label) factorized hybrid overlap sidecars require the refinement_levels = 0 legacy proxy line on the x axis",
        ),
    )
    proxy_y isa _MappedLegacyProxyLayer1D || throw(
        ArgumentError(
            "$(route_label) factorized hybrid overlap sidecars require the refinement_levels = 0 legacy proxy line on the y axis",
        ),
    )
    proxy_z isa _MappedLegacyProxyLayer1D || throw(
        ArgumentError(
            "$(route_label) factorized hybrid overlap sidecars require the refinement_levels = 0 legacy proxy line on the z axis",
        ),
    )
    norbitals = length(supplement3d.orbitals)
    x_table = zeros(Float64, size(factorized_cartesian_parent_basis.x_functions, 2), norbitals)
    y_table = zeros(Float64, size(factorized_cartesian_parent_basis.y_functions, 2), norbitals)
    z_table = zeros(Float64, size(factorized_cartesian_parent_basis.z_functions, 2), norbitals)
    for (orbital_index, orbital) in pairs(supplement3d.orbitals)
        x_data = _qwrg_atomic_axis_cross_data(proxy_x, orbital, :x, operators.expansion)
        y_data = _qwrg_atomic_axis_cross_data(proxy_y, orbital, :y, operators.expansion)
        z_data = _qwrg_atomic_axis_cross_data(proxy_z, orbital, :z, operators.expansion)
        coefficients = Vector{Float64}(orbital.coefficients)
        x_table[:, orbital_index] .=
            transpose(factorized_cartesian_parent_basis.x_functions) * (x_data.overlap * coefficients)
        y_table[:, orbital_index] .=
            transpose(factorized_cartesian_parent_basis.y_functions) * (y_data.overlap * coefficients)
        z_table[:, orbital_index] .=
            transpose(factorized_cartesian_parent_basis.z_functions) * (z_data.overlap * coefficients)
    end
    return (x = x_table, y = y_table, z = z_table)
end

function _cartesian_atomic_hybrid_overlap_sidecars(
    operators::OrdinaryCartesianOperators3D,
    cartesian_parent::CartesianBasisRepresentation3D,
)
    factorized_cartesian_parent_basis = _cartesian_factorized_parent_basis(cartesian_parent)
    parent_basis = _cartesian_hybrid_parent_basis(operators)
    parent_basis isa MappedUniformBasis || throw(
        ArgumentError(
            "factorized atomic hybrid overlap sidecars currently require a MappedUniformBasis or nested fixed block built from one",
        ),
    )
    gausslet_bundle = _mapped_ordinary_gausslet_1d_bundle(
        parent_basis;
        exponents = operators.expansion.exponents,
        center = 0.0,
        backend = operators.gausslet_backend,
    )
    supplement3d = _atomic_cartesian_shell_supplement_3d(operators.gaussian_data)
    raw_blocks = _qwrg_atomic_cartesian_blocks_3d(
        gausslet_bundle,
        supplement3d,
        operators.expansion,
    )
    parent_coefficients =
        cartesian_parent.coefficient_matrix === nothing ?
        Matrix{Float64}(I, cartesian_parent.metadata.final_dimension, cartesian_parent.metadata.final_dimension) :
        Matrix{Float64}(cartesian_parent.coefficient_matrix)
    return (
        hybrid_overlap_kind = :factorized_atomic_mixed_raw,
        factorized_cartesian_parent_basis = factorized_cartesian_parent_basis,
        cartesian_supplement_axis_tables = _cartesian_hybrid_supplement_axis_tables(
            operators,
            factorized_cartesian_parent_basis,
            supplement3d,
            _cartesian_atomic_axis_bundles(operators, parent_basis),
            "atomic",
        ),
        exact_cartesian_supplement_overlap =
            Matrix{Float64}(transpose(parent_coefficients) * raw_blocks.overlap_ga),
        exact_supplement_overlap = Matrix{Float64}(raw_blocks.overlap_aa),
    )
end

function _cartesian_diatomic_hybrid_overlap_sidecars(
    operators::OrdinaryCartesianOperators3D,
    cartesian_parent::CartesianBasisRepresentation3D,
)
    factorized_cartesian_parent_basis = _cartesian_factorized_parent_basis(cartesian_parent)
    parent_basis = _cartesian_hybrid_parent_basis(operators)
    parent_basis isa BondAlignedDiatomicQWBasis3D || throw(
        ArgumentError(
            "factorized bond-aligned diatomic hybrid overlap sidecars currently require a BondAlignedDiatomicQWBasis3D or nested fixed block built from one",
        ),
    )
    bundles = _qwrg_bond_aligned_axis_bundles(
        parent_basis,
        operators.expansion;
        gausslet_backend = operators.gausslet_backend,
    )
    supplement3d = _bond_aligned_diatomic_cartesian_shell_supplement_3d(operators.gaussian_data)
    overlap_blocks = _qwrg_diatomic_cartesian_shell_overlap_blocks_3d(
        bundles,
        supplement3d,
        parent_basis,
        operators.expansion,
    )
    parent_coefficients =
        cartesian_parent.coefficient_matrix === nothing ?
        Matrix{Float64}(I, cartesian_parent.metadata.final_dimension, cartesian_parent.metadata.final_dimension) :
        Matrix{Float64}(cartesian_parent.coefficient_matrix)
    return (
        hybrid_overlap_kind = :factorized_bond_aligned_diatomic_mixed_raw,
        factorized_cartesian_parent_basis = factorized_cartesian_parent_basis,
        cartesian_supplement_axis_tables = _cartesian_hybrid_supplement_axis_tables(
            operators,
            factorized_cartesian_parent_basis,
            supplement3d,
            (x = bundles.bundle_x, y = bundles.bundle_y, z = bundles.bundle_z),
            "bond-aligned diatomic",
        ),
        exact_cartesian_supplement_overlap =
            Matrix{Float64}(transpose(parent_coefficients) * overlap_blocks.overlap_ga),
        exact_supplement_overlap = Matrix{Float64}(overlap_blocks.overlap_aa),
    )
end

function _cartesian_hybrid_overlap_sidecars(
    operators::OrdinaryCartesianOperators3D,
    cartesian_parent::CartesianBasisRepresentation3D,
)
    if operators.gaussian_data isa LegacyAtomicGaussianSupplement
        return _cartesian_atomic_hybrid_overlap_sidecars(operators, cartesian_parent)
    end
    parent_basis = _cartesian_hybrid_parent_basis(operators)
    if parent_basis isa BondAlignedDiatomicQWBasis3D &&
       operators.gaussian_data isa Union{
           LegacyBondAlignedDiatomicGaussianSupplement,
           LegacyBondAlignedHeteronuclearGaussianSupplement,
       }
        return _cartesian_diatomic_hybrid_overlap_sidecars(operators, cartesian_parent)
    end
    throw(
        ArgumentError(
            "factorized exact hybrid overlap sidecars currently support atomic and bond-aligned diatomic routes; got parent basis $(typeof(parent_basis)) with supplement $(typeof(operators.gaussian_data))",
        ),
    )
end

function _qwrg_cartesian_parent_representation(
    operators::OrdinaryCartesianOperators3D,
)
    basis = operators.basis
    basis isa _NestedFixedBlock3D && return basis_representation(basis)
    basis isa AbstractBondAlignedOrdinaryQWBasis3D && return basis_representation(basis)
    basis isa MappedUniformBasis && return _cartesian_direct_product_representation(basis)
    throw(
        ArgumentError(
            "Cartesian basis representation does not yet support QW parent basis $(typeof(basis))",
        ),
    )
end

function basis_representation(operators::OrdinaryCartesianOperators3D)
    cartesian_parent = _qwrg_cartesian_parent_representation(operators)
    axis_representations = cartesian_parent.axis_representations
    axis_metadata = _cartesian_axis_metadata(axis_representations)
    size(cartesian_parent.metadata.basis_centers, 1) == operators.gausslet_count || throw(
        ArgumentError(
            "QW Cartesian basis representation requires the parent Cartesian representation to match gausslet_count",
        ),
    )
    supplement_representation = _cartesian_supplement_representation(operators.gaussian_data)
    parent_labels = vcat(
        cartesian_parent.metadata.basis_labels,
        [orbital.label for orbital in supplement_representation.orbitals],
    )
    parent_centers = vcat(
        cartesian_parent.metadata.basis_centers,
        _cartesian_supplement_center_matrix(supplement_representation),
    )
    size(operators.raw_to_final, 1) == length(parent_labels) || throw(
        ArgumentError(
            "QW Cartesian basis representation requires raw_to_final rows to match the explicit raw parent basis dimension",
        ),
    )
    final_labels = _cartesian_labels(operators.orbital_data)
    final_centers = _cartesian_center_matrix(operators.orbital_data)
    metadata = CartesianBasisMetadata3D(
        :hybrid_residual,
        _cartesian_axis_sharing(axis_representations),
        axis_metadata,
        :cartesian_plus_supplement_raw,
        cartesian_parent.metadata.parent_axis_counts,
        size(operators.raw_to_final, 1),
        size(operators.raw_to_final, 2),
        cartesian_parent.metadata.working_box,
        final_labels,
        final_centers,
        (
            gausslet_count = operators.gausslet_count,
            residual_count = operators.residual_count,
            supplement_kind = supplement_representation.supplement_kind,
            supplement_lmax = getfield(supplement_representation.metadata, :lmax),
        ),
    )
    overlap_sidecars = _cartesian_hybrid_overlap_sidecars(operators, cartesian_parent)
    return CartesianBasisRepresentation3D(
        metadata,
        axis_representations,
        :dense,
        Matrix{Float64}(operators.raw_to_final),
        parent_labels,
        parent_centers,
        nothing,
        nothing,
        (;
            cartesian_parent_representation = cartesian_parent,
            supplement_representation = supplement_representation,
            overlap_sidecars...,
        ),
    )
end

basis_metadata(operators::OrdinaryCartesianOperators3D) = basis_representation(operators).metadata
