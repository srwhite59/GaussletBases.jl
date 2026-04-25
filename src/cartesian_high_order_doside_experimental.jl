const _ExperimentalHighOrderCoefficientMap =
    Union{Matrix{Float64},SparseArrays.SparseMatrixCSC{Float64,Int}}

struct _ExperimentalHighOrderAxisData1D
    basis::MappedUniformBasis
    pgdg_intermediate
    backend::Symbol
    centers::Vector{Float64}
    weights::Vector{Float64}
    overlap::Matrix{Float64}
    position::Matrix{Float64}
    kinetic::Matrix{Float64}
    gaussian_factors::Vector{Matrix{Float64}}
    pair_factors_1d::Vector{Matrix{Float64}}
end

struct _ExperimentalHighOrderBlock1D
    side::Int
    interval::UnitRange{Int}
    coefficients::_ExperimentalHighOrderCoefficientMap
    local_coefficients::Matrix{Float64}
    local_overlap::Matrix{Float64}
    local_position::Matrix{Float64}
    local_weights::Vector{Float64}
    local_centers::Vector{Float64}
    localized_centers::Vector{Float64}
    inner_axis_indices::Vector{Int}
    shell_axis_indices::Vector{Int}
end

struct _ExperimentalHighOrderTensorShell3D
    side::Int
    interval::NTuple{3,UnitRange{Int}}
    full_block_coefficients::_ExperimentalHighOrderCoefficientMap
    shell_coefficients::_ExperimentalHighOrderCoefficientMap
    shell_labels::Vector{NTuple{3,Int}}
    shell_kind_counts::NamedTuple{(:faces, :edges, :corners),NTuple{3,Int}}
    column_range::UnitRange{Int}
end

struct ExperimentalHighOrderDosideStack3D
    parent_basis::MappedUniformBasis
    backend::Symbol
    doside::Int
    sides::Vector{Int}
    parent_side::Int
    coefficient_matrix::_ExperimentalHighOrderCoefficientMap
    block_column_ranges::Vector{UnitRange{Int}}
    block_labels::Vector{Symbol}
    shell_layers::Vector{_ExperimentalHighOrderTensorShell3D}
    contracted_weights::Vector{Float64}
    diagnostics::NamedTuple
end

function _experimental_high_order_centered_interval(
    n1d::Int,
    side::Int,
)
    side >= 1 || throw(ArgumentError("experimental high-order doside side lengths must be positive"))
    isodd(side) || throw(ArgumentError("experimental high-order doside lane currently requires odd side lengths"))
    side <= n1d || throw(ArgumentError("experimental high-order doside side length must lie inside the parent basis"))
    start = (n1d - side) ÷ 2 + 1
    return start:(start + side - 1)
end

function _experimental_high_order_validate_request(
    basis::MappedUniformBasis,
    sides::AbstractVector{<:Integer},
    doside::Int,
)
    mapping(basis) isa IdentityMapping || throw(
        ArgumentError("experimental high-order doside stack currently requires an undistorted IdentityMapping parent basis"),
    )
    doside == 5 || throw(ArgumentError("experimental high-order doside stack currently requires doside = 5"))
    n1d = length(basis)
    isodd(n1d) || throw(ArgumentError("experimental high-order doside stack currently requires an odd parent side length"))
    isempty(sides) && throw(ArgumentError("experimental high-order doside stack requires at least one side length"))
    side_values = Int[Int(side) for side in sides]
    first(side_values) == doside || throw(ArgumentError("experimental high-order doside side ladders must start at doside"))
    all(isodd, side_values) || throw(ArgumentError("experimental high-order doside side ladders must be odd"))
    all(side_values[index + 1] == side_values[index] + 2 for index in 1:(length(side_values) - 1)) || throw(
        ArgumentError("experimental high-order doside side ladders must increase by 2"),
    )
    maximum(side_values) == n1d || throw(
        ArgumentError("experimental high-order doside stack currently requires maximum(sides) == length(parent_basis)"),
    )
    return side_values
end

function _experimental_high_order_axis_data_1d(
    basis::MappedUniformBasis;
    backend::Symbol = :numerical_reference,
)
    if backend == :numerical_reference
        representation = basis_representation(basis; operators = (:overlap, :position, :kinetic))
        return _ExperimentalHighOrderAxisData1D(
            basis,
            nothing,
            backend,
            Float64[Float64(value) for value in centers(basis)],
            Float64[Float64(value) for value in integral_weights(basis)],
            Matrix{Float64}(representation.basis_matrices.overlap),
            Matrix{Float64}(representation.basis_matrices.position),
            Matrix{Float64}(representation.basis_matrices.kinetic),
            Matrix{Float64}[],
            Matrix{Float64}[],
        )
    end

    pgdg = _mapped_ordinary_pgdg_intermediate_1d(
        basis;
        exponents = [1.0],
        backend = backend,
    )
    return _ExperimentalHighOrderAxisData1D(
        basis,
        pgdg,
        backend,
        Float64[Float64(value) for value in pgdg.centers],
        Float64[Float64(value) for value in pgdg.weights],
        Matrix{Float64}(pgdg.overlap),
        Matrix{Float64}(pgdg.position),
        Matrix{Float64}(pgdg.kinetic),
        Matrix{Float64}[Matrix{Float64}(factor) for factor in pgdg.gaussian_factors],
        Matrix{Float64}[Matrix{Float64}(factor) for factor in pgdg.pair_factors],
    )
end

function _experimental_high_order_interval_data(
    axis_data::_ExperimentalHighOrderAxisData1D,
    interval::UnitRange{Int},
)
    if !isnothing(axis_data.pgdg_intermediate)
        return _nested_interval_data(axis_data.pgdg_intermediate, interval)
    end
    return (
        overlap = Matrix{Float64}(axis_data.overlap[interval, interval]),
        position = Matrix{Float64}(axis_data.position[interval, interval]),
        weights = Float64[axis_data.weights[index] for index in interval],
        centers = Float64[axis_data.centers[index] for index in interval],
        n1d = length(axis_data.centers),
    )
end

function _experimental_high_order_embed_local_coefficients(
    local_coefficients::AbstractMatrix{<:Real},
    interval::UnitRange{Int},
    n1d::Int,
)
    row_indices = Int[]
    col_indices = Int[]
    values = Float64[]
    for column in axes(local_coefficients, 2), (local_row, full_row) in enumerate(interval)
        value = Float64(local_coefficients[local_row, column])
        iszero(value) && continue
        push!(row_indices, full_row)
        push!(col_indices, column)
        push!(values, value)
    end
    return _nested_sparse_coefficient_map(
        row_indices,
        col_indices,
        values,
        n1d,
        size(local_coefficients, 2),
    )
end

function _experimental_high_order_block_1d(
    axis_data::_ExperimentalHighOrderAxisData1D,
    side::Int;
    doside::Int = 5,
)
    interval = _experimental_high_order_centered_interval(length(axis_data.centers), side)
    interval_data = _experimental_high_order_interval_data(axis_data, interval)
    raw_basis = _nested_retained_span(
        interval_data.weights,
        interval_data.centers,
        interval_data.position,
        interval_data.overlap,
        doside,
    )
    overlap_seed = Matrix{Float64}(transpose(raw_basis) * interval_data.overlap * raw_basis)
    position_seed = Matrix{Float64}(transpose(raw_basis) * interval_data.position * raw_basis)
    sign_vector = vec(transpose(raw_basis) * interval_data.weights)
    transform, localized_centers = _cleanup_comx_transform(overlap_seed, position_seed, sign_vector)
    local_coefficients = Matrix{Float64}(raw_basis * transform)
    local_weights = vec(transpose(interval_data.weights) * local_coefficients)
    coefficients = _experimental_high_order_embed_local_coefficients(
        local_coefficients,
        interval,
        interval_data.n1d,
    )
    return _ExperimentalHighOrderBlock1D(
        side,
        interval,
        coefficients,
        local_coefficients,
        interval_data.overlap,
        interval_data.position,
        local_weights,
        interval_data.centers,
        localized_centers,
        collect(2:(doside - 1)),
        Int[1, doside],
    )
end

function _experimental_high_order_product_coefficients(
    x_block::_ExperimentalHighOrderBlock1D,
    y_block::_ExperimentalHighOrderBlock1D,
    z_block::_ExperimentalHighOrderBlock1D,
    dims::NTuple{3,Int},
)
    nparent = prod(dims)
    row_indices = Int[]
    col_indices = Int[]
    values = Float64[]
    labels = NTuple{3,Int}[]
    column = 0
    for ixcol in axes(x_block.local_coefficients, 2),
        iycol in axes(y_block.local_coefficients, 2),
        izcol in axes(z_block.local_coefficients, 2)
        column += 1
        push!(labels, (ixcol, iycol, izcol))
        for (ix_local, ix) in enumerate(x_block.interval)
            vx = Float64(x_block.local_coefficients[ix_local, ixcol])
            iszero(vx) && continue
            for (iy_local, iy) in enumerate(y_block.interval)
                vy = Float64(y_block.local_coefficients[iy_local, iycol])
                iszero(vy) && continue
                for (iz_local, iz) in enumerate(z_block.interval)
                    vz = Float64(z_block.local_coefficients[iz_local, izcol])
                    iszero(vz) && continue
                    push!(row_indices, _cartesian_flat_index(ix, iy, iz, dims))
                    push!(col_indices, column)
                    push!(values, vx * vy * vz)
                end
            end
        end
    end
    coefficients = _nested_sparse_coefficient_map(
        row_indices,
        col_indices,
        values,
        nparent,
        length(labels),
    )
    return coefficients, labels
end

function _experimental_high_order_shell_kind_counts(
    labels::AbstractVector{<:NTuple{3,Int}};
    doside::Int = 5,
)
    faces = 0
    edges = 0
    corners = 0
    for label in labels
        boundary_axes = count(index -> index == 1 || index == doside, label)
        if boundary_axes == 1
            faces += 1
        elseif boundary_axes == 2
            edges += 1
        elseif boundary_axes == 3
            corners += 1
        end
    end
    return (faces = faces, edges = edges, corners = corners)
end

function _experimental_high_order_tensor_shell_3d(
    axis_data::_ExperimentalHighOrderAxisData1D,
    side::Int;
    doside::Int = 5,
    column_range::UnitRange{Int} = 1:0,
)
    x_block = _experimental_high_order_block_1d(axis_data, side; doside = doside)
    y_block = _experimental_high_order_block_1d(axis_data, side; doside = doside)
    z_block = _experimental_high_order_block_1d(axis_data, side; doside = doside)
    full_block_coefficients, labels = _experimental_high_order_product_coefficients(
        x_block,
        y_block,
        z_block,
        (length(axis_data.centers), length(axis_data.centers), length(axis_data.centers)),
    )

    shell_indices = Int[]
    shell_labels = NTuple{3,Int}[]
    for (column, label) in enumerate(labels)
        if any(index -> index == 1 || index == doside, label)
            push!(shell_indices, column)
            push!(shell_labels, label)
        end
    end

    shell_coefficients = full_block_coefficients[:, shell_indices]
    return _ExperimentalHighOrderTensorShell3D(
        side,
        (x_block.interval, y_block.interval, z_block.interval),
        full_block_coefficients,
        shell_coefficients,
        shell_labels,
        _experimental_high_order_shell_kind_counts(shell_labels; doside = doside),
        column_range,
    )
end

function _experimental_high_order_parent_overlap_3d(
    axis_data::_ExperimentalHighOrderAxisData1D,
)
    overlap = axis_data.overlap
    return Matrix{Float64}(kron(overlap, kron(overlap, overlap)))
end

function _experimental_high_order_parent_weights_3d(
    axis_data::_ExperimentalHighOrderAxisData1D,
)
    return _mapped_cartesian_weights(axis_data.weights)
end

function _experimental_high_order_sign_fix_columns!(
    coefficients::Matrix{Float64},
    sign_vector::AbstractVector{<:Real},
)
    for column in axes(coefficients, 2)
        marker = dot(sign_vector, view(coefficients, :, column))
        if abs(marker) <= 1.0e-12
            local_view = @view coefficients[:, column]
            pivot = argmax(abs.(local_view))
            marker = local_view[pivot]
        end
        if marker < 0.0
            coefficients[:, column] .*= -1.0
        end
    end
    return coefficients
end

function _experimental_high_order_metric_project_out(
    candidates::AbstractMatrix{<:Real},
    basis::AbstractMatrix{<:Real},
    overlap::AbstractMatrix{<:Real},
)
    basis_columns = size(basis, 2)
    basis_columns == 0 && return Matrix{Float64}(candidates)
    basis_value = Matrix{Float64}(basis)
    candidate_value = Matrix{Float64}(candidates)
    couplings = Matrix{Float64}(transpose(basis_value) * overlap * candidate_value)
    return Matrix{Float64}(candidate_value - basis_value * couplings)
end

function _experimental_high_order_lowdin_cleanup(
    coefficients::AbstractMatrix{<:Real},
    overlap::AbstractMatrix{<:Real};
    sign_vector::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    overlap_seed = Matrix{Float64}(transpose(coefficients) * overlap * coefficients)
    vectors, invhalf = _s_invsqrt_reduced(overlap_seed)
    cleaned = Matrix{Float64}(coefficients * (vectors * invhalf))
    if !isnothing(sign_vector)
        _experimental_high_order_sign_fix_columns!(cleaned, sign_vector)
    end
    return cleaned
end

function _experimental_high_order_overlap_error(
    coefficients::AbstractMatrix{<:Real},
    overlap::AbstractMatrix{<:Real},
)
    gram = Matrix{Float64}(transpose(coefficients) * overlap * coefficients)
    return norm(gram - I, Inf)
end

function _experimental_high_order_full_block_union_coefficients(
    axis_data::_ExperimentalHighOrderAxisData1D,
    sides::AbstractVector{<:Integer};
    doside::Int = 5,
)
    blocks = _ExperimentalHighOrderCoefficientMap[]
    for side in sides
        shell = _experimental_high_order_tensor_shell_3d(axis_data, Int(side); doside = doside)
        push!(blocks, shell.full_block_coefficients)
    end
    return _nested_hcat_coefficient_maps(blocks)
end

function _experimental_high_order_doside_stack_3d(
    basis::MappedUniformBasis;
    backend::Symbol = :numerical_reference,
    doside::Int = 5,
    sides::AbstractVector{<:Integer} = [5, 7, 9, 11],
)
    side_values = _experimental_high_order_validate_request(basis, sides, doside)
    axis_data = _experimental_high_order_axis_data_1d(basis; backend = backend)
    parent_overlap = _experimental_high_order_parent_overlap_3d(axis_data)
    parent_weights = _experimental_high_order_parent_weights_3d(axis_data)

    first_block = _experimental_high_order_tensor_shell_3d(axis_data, first(side_values); doside = doside)
    accumulated = Matrix{Float64}(first_block.full_block_coefficients)
    _experimental_high_order_sign_fix_columns!(accumulated, parent_weights)

    block_column_ranges = UnitRange{Int}[1:size(accumulated, 2)]
    block_labels = Symbol[Symbol("side$(first(side_values))_full")]
    shell_layers = _ExperimentalHighOrderTensorShell3D[]
    next_column = size(accumulated, 2) + 1

    for side in side_values[2:end]
        shell = _experimental_high_order_tensor_shell_3d(axis_data, side; doside = doside)
        shell_residual = _experimental_high_order_metric_project_out(
            shell.shell_coefficients,
            accumulated,
            parent_overlap,
        )
        shell_clean = _experimental_high_order_lowdin_cleanup(
            shell_residual,
            parent_overlap;
            sign_vector = parent_weights,
        )
        shell_range = next_column:(next_column + size(shell_clean, 2) - 1)
        next_column = last(shell_range) + 1

        push!(
            shell_layers,
            _ExperimentalHighOrderTensorShell3D(
                shell.side,
                shell.interval,
                shell.full_block_coefficients,
                shell.shell_coefficients,
                shell.shell_labels,
                shell.shell_kind_counts,
                shell_range,
            ),
        )
        accumulated = Matrix{Float64}(hcat(accumulated, shell_clean))
        push!(block_column_ranges, shell_range)
        push!(block_labels, Symbol("side$(side)_shell"))
    end

    contracted_weights = vec(transpose(parent_weights) * accumulated)
    diagnostics = (
        parent_dimension = size(parent_overlap, 1),
        stack_dimension = size(accumulated, 2),
        shell_dimensions = Int[length(shell.shell_labels) for shell in shell_layers],
        shell_kind_counts = [shell.shell_kind_counts for shell in shell_layers],
        overlap_error = _experimental_high_order_overlap_error(accumulated, parent_overlap),
        contracted_weights_finite = all(isfinite, contracted_weights),
    )

    return ExperimentalHighOrderDosideStack3D(
        basis,
        backend,
        doside,
        side_values,
        length(basis),
        accumulated,
        block_column_ranges,
        block_labels,
        shell_layers,
        contracted_weights,
        diagnostics,
    )
end
