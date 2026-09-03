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

function _cartesian_optional_factorized_parent_basis(
    representation::CartesianBasisRepresentation3D,
)
    try
        return _cartesian_factorized_parent_basis(representation)
    catch err
        _nested_factorized_basis_optional_failure(err) || rethrow()
        return nothing
    end
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

function _cartesian_basis_supplement_axis_primitive_cross(
    basis::BasisRepresentation1D,
    orbital::CartesianGaussianShellOrbitalRepresentation3D,
    axis::Symbol,
)
    basis_primitives = collect(primitives(primitive_set(basis)))
    all(primitive -> primitive isa Gaussian, basis_primitives) || throw(
        ArgumentError(
            "exact hybrid cartesian-supplement overlap currently requires Gaussian 1D primitives on the Cartesian axis representation",
        ),
    )
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
    primitive_cross = zeros(Float64, length(basis_primitives), length(orbital.exponents))
    for (row, primitive) in pairs(basis_primitives)
        alpha_basis = _qwrg_gaussian_exponent(primitive::Gaussian)
        for column in eachindex(orbital.exponents)
            exponent = Float64(orbital.exponents[column])
            prefactor = _qwrg_atomic_shell_prefactor(exponent, power)
            primitive_cross[row, column] = _qwrg_atomic_basic_integral(
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
    end
    return Matrix{Float64}(transpose(basis.coefficient_matrix) * primitive_cross)
end

function _cartesian_basis_supplement_cross(
    basis::CartesianBasisRepresentation3D,
    supplement::CartesianGaussianShellSupplementRepresentation3D,
)
    isempty(supplement.orbitals) && return zeros(Float64, basis.metadata.final_dimension, 0)
    parent = _cartesian_parent_state_basis(basis)
    raw_cross = zeros(Float64, length(parent.states), length(supplement.orbitals))
    for (column, orbital) in pairs(supplement.orbitals)
        overlap_x = _cartesian_basis_supplement_axis_primitive_cross(
            basis.axis_representations.x,
            orbital,
            :x,
        )
        overlap_y = _cartesian_basis_supplement_axis_primitive_cross(
            basis.axis_representations.y,
            orbital,
            :y,
        )
        overlap_z = _cartesian_basis_supplement_axis_primitive_cross(
            basis.axis_representations.z,
            orbital,
            :z,
        )
        @inbounds for (row, (ix, iy, iz)) in pairs(parent.states)
            value = 0.0
            for primitive in eachindex(orbital.coefficients)
                value +=
                    Float64(orbital.coefficients[primitive]) *
                    overlap_x[ix, primitive] *
                    overlap_y[iy, primitive] *
                    overlap_z[iz, primitive]
            end
            raw_cross[row, column] = value
        end
    end
    if parent.coefficients !== nothing
        return Matrix{Float64}(transpose(parent.coefficients) * raw_cross)
    end
    return raw_cross
end

function _cartesian_supplement_orbital_axis_overlap_matrix(
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
    return matrix
end

function _cartesian_supplement_cross_overlap(
    left::CartesianGaussianShellSupplementRepresentation3D,
    right::CartesianGaussianShellSupplementRepresentation3D,
)
    matrix = zeros(Float64, length(left.orbitals), length(right.orbitals))
    for (row, left_orbital) in pairs(left.orbitals), (column, right_orbital) in pairs(right.orbitals)
        overlap_x = _cartesian_supplement_orbital_axis_overlap_matrix(
            left_orbital,
            right_orbital,
            :x,
        )
        overlap_y = _cartesian_supplement_orbital_axis_overlap_matrix(
            left_orbital,
            right_orbital,
            :y,
        )
        overlap_z = _cartesian_supplement_orbital_axis_overlap_matrix(
            left_orbital,
            right_orbital,
            :z,
        )
        primitive_overlap =
            Matrix{Float64}(overlap_x) .* Matrix{Float64}(overlap_y) .* Matrix{Float64}(overlap_z)
        matrix[row, column] =
            Float64(dot(left_orbital.coefficients, primitive_overlap * right_orbital.coefficients))
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
