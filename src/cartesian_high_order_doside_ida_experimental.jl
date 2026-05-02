struct _ExperimentalHighOrderParentOneBodyData3D{B,O,P,H,A}
    basis::B
    backend::Symbol
    expansion::CoulombGaussianExpansion
    Z::Float64
    axis_data::A
    one_body::O
    parent_overlap::P
    parent_hamiltonian::H
    reference_energy::Float64
end

struct _ExperimentalHighOrderProjectedOneBodyData3D{B,O}
    basis::B
    backend::Symbol
    expansion::CoulombGaussianExpansion
    Z::Float64
    one_body::O
    coefficient_matrix::Matrix{Float64}
    parent_overlap::Matrix{Float64}
    parent_hamiltonian::Matrix{Float64}
    projected_overlap::Matrix{Float64}
    projected_hamiltonian::Matrix{Float64}
end

struct _ExperimentalHighOrderHePlusData3D{P}
    projected_data::P
    orbital_energies::Vector{Float64}
    ground_energy::Float64
    overlap_error::Float64
end

struct _ExperimentalHighOrderHeSingletProblem3D{P}
    projected_data::P
    parent_interaction::Matrix{Float64}
end

struct _ExperimentalHighOrderHeSingletData3D{P,D}
    problem::P
    ground_matrix::Matrix{Float64}
    ground_energy::Float64
    residual::Float64
    iterations::Int
    converged::Bool
    diagnostics::D
end

function Base.getproperty(data::_ExperimentalHighOrderHePlusData3D, name::Symbol)
    if name === :projected_data || name === :orbital_energies || name === :ground_energy || name === :overlap_error
        return getfield(data, name)
    end
    return getproperty(getfield(data, :projected_data), name)
end

function Base.propertynames(data::_ExperimentalHighOrderHePlusData3D, private::Bool = false)
    names = (:projected_data, :orbital_energies, :ground_energy, :overlap_error)
    return (names..., propertynames(getfield(data, :projected_data), private)...)
end

function Base.getproperty(problem::_ExperimentalHighOrderHeSingletProblem3D, name::Symbol)
    if name === :projected_data || name === :parent_interaction
        return getfield(problem, name)
    end
    return getproperty(getfield(problem, :projected_data), name)
end

function Base.propertynames(problem::_ExperimentalHighOrderHeSingletProblem3D, private::Bool = false)
    names = (:projected_data, :parent_interaction)
    return (names..., propertynames(getfield(problem, :projected_data), private)...)
end

function Base.getproperty(data::_ExperimentalHighOrderHeSingletData3D, name::Symbol)
    if name === :problem || name === :ground_matrix || name === :ground_energy ||
       name === :residual || name === :iterations || name === :converged || name === :diagnostics
        return getfield(data, name)
    end
    return getproperty(getfield(data, :problem), name)
end

function Base.propertynames(data::_ExperimentalHighOrderHeSingletData3D, private::Bool = false)
    names = (:problem, :ground_matrix, :ground_energy, :residual, :iterations, :converged, :diagnostics)
    return (names..., propertynames(getfield(data, :problem), private)...)
end

function _experimental_high_order_parent_one_body_data(
    basis::MappedUniformBasis;
    axis_data::Union{Nothing,_ExperimentalHighOrderAxisData1D} = nothing,
    backend::Symbol = :numerical_reference,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
    include_parent_projection_data::Bool = true,
    include_reference_energy::Bool = true,
)
    return @timeg "high_order.parent_one_body.total" begin
        axis_data_value =
            isnothing(axis_data) && backend != :numerical_reference ?
            _experimental_high_order_axis_data_1d(
                basis;
                backend = backend,
                one_body_exponents = expansion.exponents,
                one_body_center = 0.0,
            ) :
            axis_data
        one_body = @timeg "high_order.parent_one_body.build_1d_factors" begin
            if isnothing(axis_data_value)
                mapped_ordinary_one_body_operators(
                    basis;
                    exponents = expansion.exponents,
                    backend = backend,
                )
            else
                _experimental_high_order_axis_one_body_1d(
                    axis_data_value;
                    exponents = expansion.exponents,
                    center = 0.0,
                )
            end
        end
        parent_overlap, parent_hamiltonian = if include_parent_projection_data
            @timeg "high_order.parent_one_body.build_dense_parent_projection" begin
                _mapped_cartesian_one_body_matrix(one_body, expansion; Z = Z)
            end
        else
            (nothing, nothing)
        end
        reference_energy = if include_reference_energy
            @timeg "high_order.parent_one_body.reference_energy" begin
                _mapped_cartesian_hydrogen_energy(one_body, expansion; Z = Z)
            end
        else
            NaN
        end
        _ExperimentalHighOrderParentOneBodyData3D(
            basis,
            backend,
            expansion,
            Float64(Z),
            axis_data_value,
            one_body,
            parent_overlap,
            parent_hamiltonian,
            Float64(reference_energy),
        )
    end
end

function _experimental_high_order_projected_one_body_data(
    parent_data::_ExperimentalHighOrderParentOneBodyData3D,
    coefficient_matrix::AbstractMatrix{<:Real},
)
    (isnothing(parent_data.parent_overlap) || isnothing(parent_data.parent_hamiltonian)) && throw(
        ArgumentError(
            "parent projected one-body data were not built; rerun _experimental_high_order_parent_one_body_data(...; include_parent_projection_data = true)",
        ),
    )
    coefficients = Matrix{Float64}(coefficient_matrix)
    projected_overlap = _symmetrize_ida_matrix(transpose(coefficients) * parent_data.parent_overlap * coefficients)
    projected_hamiltonian = _symmetrize_ida_matrix(transpose(coefficients) * parent_data.parent_hamiltonian * coefficients)
    return _ExperimentalHighOrderProjectedOneBodyData3D(
        parent_data.basis,
        parent_data.backend,
        parent_data.expansion,
        parent_data.Z,
        parent_data.one_body,
        coefficients,
        parent_data.parent_overlap,
        parent_data.parent_hamiltonian,
        projected_overlap,
        projected_hamiltonian,
    )
end

function _experimental_high_order_contract_one_body_1d(
    one_body::MappedOrdinaryOneBody1D,
    coefficient_matrix::AbstractMatrix{<:Real},
)
    coefficients = Matrix{Float64}(coefficient_matrix)
    reduced_overlap = _symmetrize_ida_matrix(transpose(coefficients) * one_body.overlap * coefficients)
    reduced_kinetic = _symmetrize_ida_matrix(transpose(coefficients) * one_body.kinetic * coefficients)
    reduced_gaussian_factors = Matrix{Float64}[
        _symmetrize_ida_matrix(transpose(coefficients) * factor * coefficients) for factor in one_body.gaussian_factors
    ]
    return MappedOrdinaryOneBody1D(
        one_body.basis,
        one_body.backend,
        reduced_overlap,
        reduced_kinetic,
        reduced_gaussian_factors,
        copy(one_body.exponents),
        one_body.center,
    )
end

function _experimental_high_order_physical_reduced_one_body_data(
    basis::MappedUniformBasis,
    side::Int;
    axis_data::Union{Nothing,_ExperimentalHighOrderAxisData1D} = nothing,
    doside::Int = 5,
    backend::Symbol = :numerical_reference,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
    direct_comparison::Symbol = :auto,
)
    return @timeg "high_order.reduced_one_body.total" begin
        perform_direct_comparison =
            direct_comparison === :auto ? _identity_mapping(mapping(basis)) :
            direct_comparison === :always ? true :
            direct_comparison === :never ? false :
            throw(ArgumentError("direct_comparison must be :auto, :always, or :never"))
        axis_data_value = @timeg "high_order.reduced_one_body.build_axis_data" begin
            isnothing(axis_data) ?
            _experimental_high_order_axis_data_1d(
                basis;
                backend = backend,
                one_body_exponents = expansion.exponents,
                one_body_center = 0.0,
            ) :
            axis_data
        end
        physical_shell = @timeg "high_order.reduced_one_body.build_transformed_shell" begin
            _experimental_high_order_physical_shell_3d(axis_data_value, side; doside = doside)
        end
        parent_data = @timeg "high_order.reduced_one_body.build_parent_one_body" begin
            _experimental_high_order_parent_one_body_data(
                basis;
                axis_data = axis_data_value,
                backend = backend,
                expansion = expansion,
                Z = Z,
                include_parent_projection_data = perform_direct_comparison,
                include_reference_energy = false,
            )
        end
        reduced_one_body = @timeg "high_order.reduced_one_body.contract_1d_operators" begin
            _experimental_high_order_contract_one_body_1d(
                parent_data.one_body,
                physical_shell.full_block.block_1d.block.coefficients,
            )
        end
        reduced_full_overlap, reduced_full_hamiltonian = @timeg "high_order.reduced_one_body.assemble_reduced_3d_operator" begin
            _mapped_cartesian_one_body_matrix(
                reduced_one_body,
                expansion;
                Z = Z,
            )
        end
        shell_indices, _ = _experimental_high_order_tensor_shell_indices_and_labels(doside)
        reduced_shell_overlap, reduced_shell_hamiltonian = @timeg "high_order.reduced_one_body.extract_shell_subblock" begin
            (
                _symmetrize_ida_matrix(reduced_full_overlap[shell_indices, shell_indices]),
                _symmetrize_ida_matrix(reduced_full_hamiltonian[shell_indices, shell_indices]),
            )
        end
        direct_full, direct_shell = if perform_direct_comparison
            @timeg "high_order.reduced_one_body.build_dense_parent_reference" begin
                (
                    _experimental_high_order_projected_one_body_data(
                        parent_data,
                        Matrix{Float64}(physical_shell.full_block.shell.full_block_coefficients),
                    ),
                    _experimental_high_order_projected_one_body_data(
                        parent_data,
                        Matrix{Float64}(physical_shell.shell.shell_coefficients),
                    ),
                )
            end
        else
            (nothing, nothing)
        end
        full_summary, shell_summary, direct_full_summary, direct_shell_summary = @timeg "high_order.reduced_one_body.heplus_summaries" begin
            reduced_full_summary = _experimental_high_order_heplus_summary(
                reduced_full_overlap,
                reduced_full_hamiltonian,
            )
            reduced_shell_summary = _experimental_high_order_heplus_summary(
                reduced_shell_overlap,
                reduced_shell_hamiltonian,
            )
            reference_full_summary = if perform_direct_comparison
                _experimental_high_order_heplus_summary(
                    direct_full.projected_overlap,
                    direct_full.projected_hamiltonian,
                )
            else
                nothing
            end
            reference_shell_summary = if perform_direct_comparison
                _experimental_high_order_heplus_summary(
                    direct_shell.projected_overlap,
                    direct_shell.projected_hamiltonian,
                )
            else
                nothing
            end
            (
                reduced_full_summary,
                reduced_shell_summary,
                reference_full_summary,
                reference_shell_summary,
            )
        end
        diagnostics = (
            direct_comparison_performed = perform_direct_comparison,
            full_overlap_error = perform_direct_comparison ? norm(reduced_full_overlap - direct_full.projected_overlap, Inf) : NaN,
            full_hamiltonian_error = perform_direct_comparison ? norm(reduced_full_hamiltonian - direct_full.projected_hamiltonian, Inf) : NaN,
            shell_overlap_error = perform_direct_comparison ? norm(reduced_shell_overlap - direct_shell.projected_overlap, Inf) : NaN,
            shell_hamiltonian_error = perform_direct_comparison ? norm(reduced_shell_hamiltonian - direct_shell.projected_hamiltonian, Inf) : NaN,
            full_summary = full_summary,
            shell_summary = shell_summary,
            direct_full_summary = direct_full_summary,
            direct_shell_summary = direct_shell_summary,
        )
        (
            basis = basis,
            side = side,
            doside = doside,
            axis_data = axis_data_value,
            parent_data = parent_data,
            physical_shell = physical_shell,
            reduced_one_body = reduced_one_body,
            reduced_full_overlap = reduced_full_overlap,
            reduced_full_hamiltonian = reduced_full_hamiltonian,
            reduced_shell_overlap = reduced_shell_overlap,
            reduced_shell_hamiltonian = reduced_shell_hamiltonian,
            direct_full = direct_full,
            direct_shell = direct_shell,
            diagnostics = diagnostics,
        )
    end
end

function _experimental_high_order_projected_one_body_append(
    projected::_ExperimentalHighOrderProjectedOneBodyData3D,
    new_columns::AbstractMatrix{<:Real},
)
    appended_columns = Matrix{Float64}(new_columns)
    overlap_cross = Matrix{Float64}(transpose(projected.coefficient_matrix) * projected.parent_overlap * appended_columns)
    overlap_self = Matrix{Float64}(transpose(appended_columns) * projected.parent_overlap * appended_columns)
    hamiltonian_cross = Matrix{Float64}(transpose(projected.coefficient_matrix) * projected.parent_hamiltonian * appended_columns)
    hamiltonian_self = Matrix{Float64}(transpose(appended_columns) * projected.parent_hamiltonian * appended_columns)
    return _ExperimentalHighOrderProjectedOneBodyData3D(
        projected.basis,
        projected.backend,
        projected.expansion,
        projected.Z,
        projected.one_body,
        Matrix{Float64}(hcat(projected.coefficient_matrix, appended_columns)),
        projected.parent_overlap,
        projected.parent_hamiltonian,
        _symmetrize_ida_matrix([
            projected.projected_overlap overlap_cross
            transpose(overlap_cross) overlap_self
        ]),
        _symmetrize_ida_matrix([
            projected.projected_hamiltonian hamiltonian_cross
            transpose(hamiltonian_cross) hamiltonian_self
        ]),
    )
end

function _experimental_high_order_projected_one_body_data(
    basis::MappedUniformBasis,
    coefficient_matrix::AbstractMatrix{<:Real};
    backend::Symbol = :numerical_reference,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
)
    parent_data = _experimental_high_order_parent_one_body_data(
        basis;
        backend = backend,
        expansion = expansion,
        Z = Z,
    )
    return _experimental_high_order_projected_one_body_data(parent_data, coefficient_matrix)
end

function _experimental_high_order_projected_one_body_data(
    stack::ExperimentalHighOrderDosideStack3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
)
    return _experimental_high_order_projected_one_body_data(
        stack.parent_basis,
        stack.coefficient_matrix;
        backend = stack.backend,
        expansion = expansion,
        Z = Z,
    )
end

function _experimental_high_order_orthonormalized_full_block_union_coefficients(
    axis_data::_ExperimentalHighOrderAxisData1D,
    sides::AbstractVector{<:Integer};
    doside::Int = 5,
)
    parent_overlap = _experimental_high_order_parent_overlap_3d(axis_data)
    parent_weights = _experimental_high_order_parent_weights_3d(axis_data)
    union_coefficients = Matrix{Float64}(
        _experimental_high_order_full_block_union_coefficients(
            axis_data,
            sides;
            doside = doside,
        ),
    )
    return _experimental_high_order_lowdin_cleanup(
        union_coefficients,
        parent_overlap;
        sign_vector = parent_weights,
    )
end

function _experimental_high_order_doside_heplus_data(
    parent_data::_ExperimentalHighOrderParentOneBodyData3D,
    coefficient_matrix::AbstractMatrix{<:Real},
)
    projected = _experimental_high_order_projected_one_body_data(parent_data, coefficient_matrix)
    decomposition = eigen(Symmetric(projected.projected_hamiltonian))
    orbital_energies = Float64[Float64(value) for value in decomposition.values]
    return _ExperimentalHighOrderHePlusData3D(
        projected,
        orbital_energies,
        Float64(decomposition.values[1]),
        norm(projected.projected_overlap - I, Inf),
    )
end

function _experimental_high_order_doside_heplus_data(
    basis::MappedUniformBasis,
    coefficient_matrix::AbstractMatrix{<:Real};
    backend::Symbol = :numerical_reference,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
)
    parent_data = _experimental_high_order_parent_one_body_data(
        basis,
        backend = backend,
        expansion = expansion,
        Z = Z,
    )
    return _experimental_high_order_doside_heplus_data(parent_data, coefficient_matrix)
end

function _experimental_high_order_doside_heplus_data(
    stack::ExperimentalHighOrderDosideStack3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
)
    return _experimental_high_order_doside_heplus_data(
        stack.parent_basis,
        stack.coefficient_matrix;
        backend = stack.backend,
        expansion = expansion,
        Z = Z,
    )
end

function _experimental_high_order_doside_heplus_energy(
    basis::MappedUniformBasis,
    coefficient_matrix::AbstractMatrix{<:Real};
    backend::Symbol = :numerical_reference,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
)
    return _experimental_high_order_doside_heplus_data(
        basis,
        coefficient_matrix;
        backend = backend,
        expansion = expansion,
        Z = Z,
    ).ground_energy
end

function _experimental_high_order_doside_heplus_energy(
    stack::ExperimentalHighOrderDosideStack3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
)
    return _experimental_high_order_doside_heplus_data(
        stack;
        expansion = expansion,
        Z = Z,
    ).ground_energy
end

function _experimental_high_order_heplus_benchmark_row(
    route::Symbol,
    outer_side::Int,
    function_count::Int,
    energy::Float64,
    reference_energy::Float64,
    overlap_error::Float64;
    support_overlap_error::Float64,
    support_overlap_minimum_eigenvalue::Float64 = NaN,
    support_overlap_maximum_eigenvalue::Float64 = NaN,
)
    return (
        route = route,
        outer_side = outer_side,
        function_count = function_count,
        energy = energy,
        error = abs(energy - reference_energy),
        overlap_error = overlap_error,
        support_overlap_error = support_overlap_error,
        support_overlap_minimum_eigenvalue = support_overlap_minimum_eigenvalue,
        support_overlap_maximum_eigenvalue = support_overlap_maximum_eigenvalue,
    )
end

function _experimental_high_order_heplus_summary(
    projected_overlap::AbstractMatrix{<:Real},
    projected_hamiltonian::AbstractMatrix{<:Real},
)
    hamiltonian = _symmetrize_ida_matrix(Matrix{Float64}(projected_hamiltonian))
    overlap = _symmetrize_ida_matrix(Matrix{Float64}(projected_overlap))
    values = eigvals(Symmetric(hamiltonian))
    return (
        ground_energy = Float64(values[1]),
        overlap_error = norm(overlap - I, Inf),
        overlap_spectrum = _experimental_high_order_positive_spectrum(overlap; tol = 1.0e-10),
    )
end

function _experimental_high_order_parent_metric_transfer(
    parent_overlap::AbstractMatrix{<:Real},
    big_coefficients::AbstractMatrix{<:Real},
    small_coefficients::AbstractMatrix{<:Real},
)
    big = Matrix{Float64}(big_coefficients)
    small = Matrix{Float64}(small_coefficients)
    big_overlap = _symmetrize_ida_matrix(transpose(big) * parent_overlap * big)
    small_overlap = _symmetrize_ida_matrix(transpose(small) * parent_overlap * small)
    transfer = Matrix{Float64}(transpose(big) * parent_overlap * small)
    residual = small - big * transfer
    return (
        transfer = transfer,
        coefficient_residual = norm(residual, Inf),
        metric_residual = norm(
            _symmetrize_ida_matrix(transpose(residual) * parent_overlap * residual),
            Inf,
        ),
        overlap_reconstruction_error = norm(
            _symmetrize_ida_matrix(transpose(transfer) * big_overlap * transfer) - small_overlap,
            Inf,
        ),
    )
end

function _experimental_high_order_transfer_heplus_data(
    big_data::_ExperimentalHighOrderHePlusData3D,
    coefficient_matrix::AbstractMatrix{<:Real},
    transfer::AbstractMatrix{<:Real},
)
    projected = _ExperimentalHighOrderProjectedOneBodyData3D(
        big_data.basis,
        big_data.backend,
        big_data.expansion,
        big_data.Z,
        big_data.one_body,
        Matrix{Float64}(coefficient_matrix),
        big_data.parent_overlap,
        big_data.parent_hamiltonian,
        _symmetrize_ida_matrix(transpose(transfer) * big_data.projected_overlap * transfer),
        _symmetrize_ida_matrix(transpose(transfer) * big_data.projected_hamiltonian * transfer),
    )
    decomposition = eigen(Symmetric(projected.projected_hamiltonian))
    orbital_energies = Float64[Float64(value) for value in decomposition.values]
    return _ExperimentalHighOrderHePlusData3D(
        projected,
        orbital_energies,
        Float64(decomposition.values[1]),
        norm(projected.projected_overlap - I, Inf),
    )
end

function _experimental_high_order_route_transfer_rows(
    parent_overlap::AbstractMatrix{<:Real},
    route_data::AbstractVector{<:NamedTuple};
    transfer_tolerance::Real = 1.0e-8,
)
    largest = last(route_data)
    rows = NamedTuple[]
    for row in route_data
        if row.outer_side == largest.outer_side
            push!(
                rows,
                (
                    outer_side = row.outer_side,
                    function_count = row.function_count,
                    construction = :identity,
                    admitted = true,
                    used_for_reuse = true,
                    coefficient_residual = 0.0,
                    metric_residual = 0.0,
                    overlap_reconstruction_error = 0.0,
                    transfer = Matrix{Float64}(I, row.function_count, row.function_count),
                ),
            )
            continue
        end

        transfer_data = _experimental_high_order_parent_metric_transfer(
            parent_overlap,
            largest.coefficient_matrix,
            row.coefficient_matrix,
        )
        admitted =
            transfer_data.metric_residual <= Float64(transfer_tolerance) &&
            transfer_data.overlap_reconstruction_error <= Float64(transfer_tolerance)
        push!(
            rows,
            (
                outer_side = row.outer_side,
                function_count = row.function_count,
                construction = :cross_overlap_parent_metric,
                admitted = admitted,
                used_for_reuse = admitted,
                coefficient_residual = transfer_data.coefficient_residual,
                metric_residual = transfer_data.metric_residual,
                overlap_reconstruction_error = transfer_data.overlap_reconstruction_error,
                transfer = transfer_data.transfer,
            ),
        )
    end
    return rows
end

function _experimental_high_order_route_benchmark_rows(
    parent_data::_ExperimentalHighOrderParentOneBodyData3D,
    route::Symbol,
    route_data::AbstractVector{<:NamedTuple};
    reference_energy::Float64,
    reuse_ladder_transfers::Bool = true,
    transfer_tolerance::Real = 1.0e-8,
)
    transfer_rows = _experimental_high_order_route_transfer_rows(
        parent_data.parent_overlap,
        route_data;
        transfer_tolerance = transfer_tolerance,
    )
    largest = last(route_data)
    largest_data = _experimental_high_order_doside_heplus_data(
        parent_data,
        largest.coefficient_matrix,
    )
    rows = NamedTuple[]
    for (row, transfer_row) in zip(route_data, transfer_rows)
        used_for_reuse = reuse_ladder_transfers && transfer_row.used_for_reuse
        heplus_data =
            row.outer_side == largest.outer_side ? largest_data :
            used_for_reuse ? _experimental_high_order_transfer_heplus_data(
                largest_data,
                row.coefficient_matrix,
                transfer_row.transfer,
            ) :
            _experimental_high_order_doside_heplus_data(
                parent_data,
                row.coefficient_matrix,
            )
        push!(
            rows,
            _experimental_high_order_heplus_benchmark_row(
                route,
                row.outer_side,
                row.function_count,
                heplus_data.ground_energy,
                reference_energy,
                heplus_data.overlap_error;
                support_overlap_error = row.support_overlap_error,
                support_overlap_minimum_eigenvalue = row.support_overlap_minimum_eigenvalue,
                support_overlap_maximum_eigenvalue = row.support_overlap_maximum_eigenvalue,
            ),
        )
    end
    transfer_rows = [
        merge(
            transfer_row,
            (
                used_for_reuse = reuse_ladder_transfers && transfer_row.used_for_reuse,
            ),
        ) for transfer_row in transfer_rows
    ]
    return rows, transfer_rows
end

function _experimental_high_order_high_route_data(
    parent_data::_ExperimentalHighOrderParentOneBodyData3D;
    doside::Int = 5,
    outer_sides::AbstractVector{<:Integer} = [5, 7, 9, 11],
)
    basis = parent_data.basis
    side_values, mapping_family = _experimental_high_order_validate_request(basis, outer_sides, doside)
    axis_data = isnothing(parent_data.axis_data) ?
        _experimental_high_order_axis_data_1d(basis; backend = parent_data.backend) :
        parent_data.axis_data
    parent_overlap = _experimental_high_order_parent_overlap_3d(axis_data)
    parent_weights = _experimental_high_order_parent_weights_3d(axis_data)
    parent_side = length(basis)

    tensor_shells = Dict{Int,_ExperimentalHighOrderTensorShell3D}()
    for side in side_values
        block = _experimental_high_order_physical_block_1d(axis_data, side; doside = doside).block
        tensor_shells[side] = _experimental_high_order_tensor_shell_3d(
            block,
            parent_side;
            doside = doside,
        )
    end

    first_shell = tensor_shells[first(side_values)]
    accumulated = Matrix{Float64}(first_shell.full_block_coefficients)
    _experimental_high_order_sign_fix_columns!(accumulated, parent_weights)

    route_data = NamedTuple[]
    overlap = _symmetrize_ida_matrix(transpose(accumulated) * parent_overlap * accumulated)
    spectrum = _experimental_high_order_positive_spectrum(overlap; tol = 1.0e-10)
    push!(
        route_data,
        (
            outer_side = first(side_values),
            coefficient_matrix = Matrix{Float64}(accumulated),
            function_count = size(accumulated, 2),
            support_overlap_error = norm(overlap - I, Inf),
            support_overlap_minimum_eigenvalue = spectrum.minimum_eigenvalue,
            support_overlap_maximum_eigenvalue = spectrum.maximum_eigenvalue,
        ),
    )

    for side in side_values[2:end]
        shell = tensor_shells[side]
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
        accumulated = Matrix{Float64}(hcat(accumulated, shell_clean))
        overlap = _symmetrize_ida_matrix(transpose(accumulated) * parent_overlap * accumulated)
        spectrum = _experimental_high_order_positive_spectrum(overlap; tol = 1.0e-10)
        push!(
            route_data,
            (
                outer_side = side,
                coefficient_matrix = Matrix{Float64}(accumulated),
                function_count = size(accumulated, 2),
                support_overlap_error = norm(overlap - I, Inf),
                support_overlap_minimum_eigenvalue = spectrum.minimum_eigenvalue,
                support_overlap_maximum_eigenvalue = spectrum.maximum_eigenvalue,
            ),
        )
    end

    return mapping_family, route_data
end

function _experimental_high_order_lower_route_data(
    gausslet_bundle::_MappedOrdinaryGausslet1DBundle;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    outer_sides::AbstractVector{<:Integer} = [5, 7, 9, 11],
    comparator_nside::Int = 5,
)
    route_data = NamedTuple[]
    for outer_side in outer_sides
        working_box = _experimental_high_order_centered_working_box(length(gausslet_bundle.basis), outer_side)
        lower_fixed = one_center_atomic_legacy_profile_fixed_block(
            gausslet_bundle;
            expansion = expansion,
            working_box = working_box,
            nside = comparator_nside,
        )
        overlap = _symmetrize_ida_matrix(Matrix{Float64}(lower_fixed.overlap))
        spectrum = _experimental_high_order_positive_spectrum(overlap; tol = 1.0e-10)
        push!(
            route_data,
            (
                outer_side = outer_side,
                coefficient_matrix = Matrix{Float64}(lower_fixed.coefficient_matrix),
                function_count = size(lower_fixed.coefficient_matrix, 2),
                support_overlap_error = norm(overlap - I, Inf),
                support_overlap_minimum_eigenvalue = spectrum.minimum_eigenvalue,
                support_overlap_maximum_eigenvalue = spectrum.maximum_eigenvalue,
            ),
        )
    end
    return route_data
end

function _experimental_high_order_distorted_parent_heplus_benchmark(
    basis::MappedUniformBasis;
    backend::Symbol = :numerical_reference,
    doside::Int = 5,
    outer_sides::AbstractVector{<:Integer} = [5, 7, 9, 11],
    comparator_nside::Int = 5,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
    reuse_ladder_transfers::Bool = true,
    transfer_tolerance::Real = 1.0e-8,
)
    mapping_family = _experimental_high_order_mapping_family(basis)
    mapping_family == :white_lindsey_atomic_he_d0p2 || throw(
        ArgumentError(
            "distorted-parent experimental He+ benchmark currently requires the explicit White-Lindsey He mapping family",
        ),
    )

    side_values = Int[Int(side) for side in outer_sides]
    gausslet_bundle = _mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = backend,
    )
    axis_data = _experimental_high_order_axis_data_1d(
        basis;
        backend = backend,
        prepared_bundle = gausslet_bundle,
        one_body_exponents = expansion.exponents,
        one_body_center = 0.0,
    )
    parent_data = _experimental_high_order_parent_one_body_data(
        basis;
        axis_data = axis_data,
        backend = backend,
        expansion = expansion,
        Z = Z,
    )
    reference_energy = parent_data.reference_energy

    lower_route_data = _experimental_high_order_lower_route_data(
        gausslet_bundle;
        expansion = expansion,
        outer_sides = side_values,
        comparator_nside = comparator_nside,
    )
    mapping_family, high_route_data = _experimental_high_order_high_route_data(
        parent_data;
        doside = doside,
        outer_sides = side_values,
    )

    lower_order_rows, lower_order_transfer_rows = _experimental_high_order_route_benchmark_rows(
        parent_data,
        :lower_order_atomic_legacy_profile,
        lower_route_data;
        reference_energy = reference_energy,
        reuse_ladder_transfers = reuse_ladder_transfers,
        transfer_tolerance = transfer_tolerance,
    )
    high_order_rows, high_order_transfer_rows = _experimental_high_order_route_benchmark_rows(
        parent_data,
        :high_order_doside,
        high_route_data;
        reference_energy = reference_energy,
        reuse_ladder_transfers = reuse_ladder_transfers,
        transfer_tolerance = transfer_tolerance,
    )
    comparison_rows = NamedTuple[]

    for index in eachindex(side_values)
        lower_row = lower_order_rows[index]
        high_row = high_order_rows[index]
        lower_row.function_count == high_row.function_count || throw(
            ArgumentError(
                "distorted-parent He+ benchmark requires aligned retained-function counts; got lower-order $(lower_row.function_count) and high-order $(high_row.function_count) at outer side $(side_values[index])",
            ),
        )
        lower_error = lower_row.error
        high_error = high_row.error
        push!(
            comparison_rows,
            (
                outer_side = side_values[index],
                function_count = high_row.function_count,
                lower_order_error = lower_error,
                high_order_error = high_error,
                winner =
                    high_error < lower_error ? :high_order_doside :
                    lower_error < high_error ? :lower_order_atomic_legacy_profile :
                    :tie,
                error_ratio = lower_error > 0.0 ? high_error / lower_error : Inf,
            ),
        )
    end

    return (
        mapping_family = mapping_family,
        parent_side = length(basis),
        parent_reference_dimension = length(basis)^3,
        parent_reference_energy = reference_energy,
        lower_order_comparator = :one_center_atomic_legacy_profile_fixed_block,
        comparator_nside = comparator_nside,
        lower_order_rows = lower_order_rows,
        high_order_rows = high_order_rows,
        comparison_rows = comparison_rows,
        lower_order_transfer_rows = lower_order_transfer_rows,
        high_order_transfer_rows = high_order_transfer_rows,
    )
end

function _experimental_high_order_he_singlet_problem(
    basis::MappedUniformBasis,
    coefficient_matrix::AbstractMatrix{<:Real};
    backend::Symbol = :numerical_reference,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
)
    projected = _experimental_high_order_projected_one_body_data(
        basis,
        coefficient_matrix;
        backend = backend,
        expansion = expansion,
        Z = Z,
    )
    ida_parent = ordinary_cartesian_ida_operators(
        basis;
        expansion = expansion,
        Z = Z,
        backend = backend,
    )
    return _ExperimentalHighOrderHeSingletProblem3D(
        projected,
        Matrix{Float64}(ida_parent.interaction_matrix),
    )
end

function _experimental_high_order_he_singlet_problem(
    stack::ExperimentalHighOrderDosideStack3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
)
    return _experimental_high_order_he_singlet_problem(
        stack.parent_basis,
        stack.coefficient_matrix;
        backend = stack.backend,
        expansion = expansion,
        Z = Z,
    )
end

function _experimental_high_order_frobenius_dot(
    left::AbstractMatrix{<:Real},
    right::AbstractMatrix{<:Real},
)
    size(left) == size(right) || throw(DimensionMismatch("Frobenius inner product requires matching matrix sizes"))
    return Float64(sum(Float64.(left) .* Float64.(right)))
end

function _experimental_high_order_frobenius_norm(matrix::AbstractMatrix{<:Real})
    return sqrt(max(_experimental_high_order_frobenius_dot(matrix, matrix), 0.0))
end

function _experimental_high_order_he_column_participation(
    ground_matrix::AbstractMatrix{<:Real},
)
    size(ground_matrix, 1) == size(ground_matrix, 2) || throw(
        ArgumentError("He column participation requires a square ground-state matrix"),
    )
    density = _symmetrize_ida_matrix(2.0 .* Matrix{Float64}(ground_matrix * transpose(ground_matrix)))
    trace_value = Float64(sum(density[index, index] for index in axes(density, 1)))
    trace_value > 0.0 || throw(ArgumentError("He column participation requires a nonzero one-particle density"))
    return (
        density = density,
        trace_value = trace_value,
        column_fractions = Float64[density[index, index] / trace_value for index in axes(density, 1)],
    )
end

function _experimental_high_order_he_shell_participation(
    ground_matrix::AbstractMatrix{<:Real},
    block_labels::AbstractVector{<:Symbol},
    block_column_ranges::AbstractVector{<:UnitRange{Int}},
)
    size(ground_matrix, 1) == size(ground_matrix, 2) || throw(
        ArgumentError("He shell participation requires a square ground-state matrix"),
    )
    length(block_labels) == length(block_column_ranges) || throw(
        ArgumentError("He shell participation requires matching block labels and column ranges"),
    )
    column_participation = _experimental_high_order_he_column_participation(ground_matrix)
    density = column_participation.density
    trace_value = column_participation.trace_value

    occupation_fractions = Float64[]
    peak_orbital_fractions = Float64[]
    for column_range in block_column_ranges
        block_trace = Float64(sum(density[index, index] for index in column_range))
        push!(occupation_fractions, max(block_trace / trace_value, 0.0))
        push!(
            peak_orbital_fractions,
            maximum(Float64[density[index, index] / trace_value for index in column_range]),
        )
    end

    outer_fractions = occupation_fractions[2:end]
    outer_peak_fractions = peak_orbital_fractions[2:end]
    return (
        block_labels = Symbol[block_labels...],
        occupation_fractions = occupation_fractions,
        peak_orbital_fractions = peak_orbital_fractions,
        total_core_fraction = first(occupation_fractions),
        total_outer_shell_fraction = isempty(outer_fractions) ? 0.0 : sum(outer_fractions),
        outermost_shell_fraction = isempty(outer_fractions) ? 0.0 : last(outer_fractions),
        largest_outer_shell_fraction = isempty(outer_fractions) ? 0.0 : maximum(outer_fractions),
        largest_outer_orbital_fraction = isempty(outer_peak_fractions) ? 0.0 : maximum(outer_peak_fractions),
    )
end

function _experimental_high_order_he_shell_participation(
    stack::ExperimentalHighOrderDosideStack3D,
    ground_matrix::AbstractMatrix{<:Real},
)
    return _experimental_high_order_he_shell_participation(
        ground_matrix,
        stack.block_labels,
        stack.block_column_ranges,
    )
end

function _experimental_high_order_he_moment_risk_audit(
    stack::ExperimentalHighOrderDosideStack3D,
    ground_matrix::AbstractMatrix{<:Real};
    top_ks::AbstractVector{<:Integer} = [1, 5, 10],
)
    column_participation = _experimental_high_order_he_column_participation(ground_matrix)
    column_fractions = column_participation.column_fractions
    column_diagnostics = stack.diagnostics.moment_risk.column_diagnostics
    outer_ranked = sort(
        NamedTuple[
            merge(
                diagnostic,
                (
                    state_weight = column_fractions[diagnostic.column_index],
                ),
            ) for diagnostic in column_diagnostics if diagnostic.is_outer_shell
        ];
        by = diagnostic -> (-diagnostic.risk_score, -diagnostic.max_normalized_mu4, -diagnostic.max_normalized_mu3, diagnostic.column_index),
    )
    worst_outer_direction = isempty(outer_ranked) ? nothing : first(outer_ranked)
    top_k_weights = NamedTuple[]
    for k in top_ks
        k >= 1 || throw(ArgumentError("He moment-risk audit requires positive top-k sizes"))
        push!(
            top_k_weights,
            (
                k = Int(k),
                total_state_weight = isempty(outer_ranked) ? 0.0 : sum(diagnostic.state_weight for diagnostic in first(outer_ranked, min(Int(k), length(outer_ranked)))),
            ),
        )
    end
    shell_block_weights = NamedTuple[
        (
            block_label = stack.block_labels[block_index],
            total_state_weight = sum(column_fractions[index] for index in column_range),
        ) for (block_index, column_range) in enumerate(stack.block_column_ranges)
    ]
    return (
        outer_column_count = length(outer_ranked),
        worst_outer_direction = worst_outer_direction,
        top_k_weights = top_k_weights,
        shell_block_weights = shell_block_weights,
        top_ranked_outer_directions = first(outer_ranked, min(10, length(outer_ranked))),
    )
end

function _experimental_high_order_he_singlet_action(
    problem::_ExperimentalHighOrderHeSingletProblem3D,
    coefficients::AbstractMatrix{<:Real},
)
    nstack = size(problem.coefficient_matrix, 2)
    size(coefficients) == (nstack, nstack) || throw(
        DimensionMismatch("He singlet action requires an $(nstack)x$(nstack) symmetric coefficient matrix"),
    )
    symmetric_coefficients = _symmetrize_ida_matrix(coefficients)
    parent_density = Matrix{Float64}(problem.coefficient_matrix * symmetric_coefficients * transpose(problem.coefficient_matrix))
    interaction = Matrix{Float64}(
        transpose(problem.coefficient_matrix) * (problem.parent_interaction .* parent_density) * problem.coefficient_matrix
    )
    one_body = Matrix{Float64}(
        problem.projected_hamiltonian * symmetric_coefficients +
        symmetric_coefficients * problem.projected_hamiltonian
    )
    return _symmetrize_ida_matrix(one_body + interaction)
end

function _experimental_high_order_he_singlet_lanczos(
    problem::_ExperimentalHighOrderHeSingletProblem3D;
    krylovdim::Int = 32,
    maxiter::Int = 32,
    tol::Real = 1.0e-8,
    A0::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
)
    nstack = size(problem.projected_hamiltonian, 1)
    nstack >= 1 || throw(ArgumentError("experimental high-order He Lanczos requires a nonempty stack"))
    krylovdim >= 2 || throw(ArgumentError("experimental high-order He Lanczos requires krylovdim >= 2"))
    maxiter >= 1 || throw(ArgumentError("experimental high-order He Lanczos requires maxiter >= 1"))
    tol > 0 || throw(ArgumentError("experimental high-order He Lanczos requires tol > 0"))

    start_matrix = if A0 === nothing
        decomposition = eigen(Symmetric(problem.projected_hamiltonian))
        lowest = Vector{Float64}(decomposition.vectors[:, 1])
        Matrix{Float64}(lowest * transpose(lowest))
    else
        size(A0) == (nstack, nstack) || throw(
            DimensionMismatch("initial He Lanczos matrix must match the projected stack dimension"),
        )
        Matrix{Float64}(A0)
    end
    start_matrix = _symmetrize_ida_matrix(start_matrix)
    start_norm = _experimental_high_order_frobenius_norm(start_matrix)
    start_norm > 0.0 || throw(ArgumentError("experimental high-order He Lanczos requires a nonzero initial matrix"))
    current = start_matrix ./ start_norm

    vectors = Matrix{Float64}[copy(current)]
    alpha = Float64[]
    beta = Float64[]
    previous = zeros(Float64, nstack, nstack)

    converged = false
    residual = Inf
    iterations = 0
    best_small_vector = ones(Float64, 1)
    best_value = NaN

    maxsteps = min(krylovdim, maxiter)
    for step in 1:maxsteps
        iterations = step
        w = _experimental_high_order_he_singlet_action(problem, current)
        step > 1 && (w .-= beta[end] .* previous)

        a = _experimental_high_order_frobenius_dot(current, w)
        push!(alpha, a)
        w .-= a .* current

        for basis_vector in vectors
            w .-= _experimental_high_order_frobenius_dot(basis_vector, w) .* basis_vector
        end

        w = _symmetrize_ida_matrix(w)
        b = _experimental_high_order_frobenius_norm(w)
        small_eig = eigen(SymTridiagonal(alpha, beta))
        best_value = Float64(real(small_eig.values[1]))
        best_small_vector = Vector{Float64}(small_eig.vectors[:, 1])
        residual = abs(b * best_small_vector[end])

        if residual <= tol || step == maxsteps || b <= sqrt(eps(Float64))
            converged = residual <= tol || b <= sqrt(eps(Float64))
            break
        end

        push!(beta, b)
        previous = current
        current = w ./ b
        push!(vectors, copy(current))
    end

    ground_matrix = zeros(Float64, nstack, nstack)
    for index in eachindex(best_small_vector)
        ground_matrix .+= best_small_vector[index] .* vectors[index]
    end
    ground_matrix = _symmetrize_ida_matrix(ground_matrix)
    ground_matrix ./= _experimental_high_order_frobenius_norm(ground_matrix)

    return (
        value = best_value,
        matrix = ground_matrix,
        residual = residual,
        iterations = iterations,
        converged = converged,
    )
end

function _experimental_high_order_doside_he_singlet_data(
    basis::MappedUniformBasis,
    coefficient_matrix::AbstractMatrix{<:Real};
    backend::Symbol = :numerical_reference,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
    krylovdim::Int = 32,
    maxiter::Int = 32,
    tol::Real = 1.0e-8,
)
    problem = _experimental_high_order_he_singlet_problem(
        basis,
        coefficient_matrix;
        backend = backend,
        expansion = expansion,
        Z = Z,
    )
    result = _experimental_high_order_he_singlet_lanczos(
        problem;
        krylovdim = krylovdim,
        maxiter = maxiter,
        tol = tol,
    )
    diagnostics = (
        projected_overlap_spectrum = _experimental_high_order_positive_spectrum(
            problem.projected_overlap;
            tol = 1.0e-10,
        ),
        shell_participation = _experimental_high_order_he_shell_participation(
            result.matrix,
            Symbol[:all],
            UnitRange{Int}[1:size(result.matrix, 1)],
        ),
        moment_risk_audit = nothing,
    )
    return _ExperimentalHighOrderHeSingletData3D(
        problem,
        result.matrix,
        result.value,
        result.residual,
        result.iterations,
        result.converged,
        diagnostics,
    )
end

function _experimental_high_order_doside_he_singlet_data(
    stack::ExperimentalHighOrderDosideStack3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
    krylovdim::Int = 32,
    maxiter::Int = 32,
    tol::Real = 1.0e-8,
)
    problem = _experimental_high_order_he_singlet_problem(
        stack;
        expansion = expansion,
        Z = Z,
    )
    result = _experimental_high_order_he_singlet_lanczos(
        problem;
        krylovdim = krylovdim,
        maxiter = maxiter,
        tol = tol,
    )
    diagnostics = (
        projected_overlap_spectrum = _experimental_high_order_positive_spectrum(
            problem.projected_overlap;
            tol = 1.0e-10,
        ),
        shell_participation = _experimental_high_order_he_shell_participation(stack, result.matrix),
        moment_risk_audit = _experimental_high_order_he_moment_risk_audit(stack, result.matrix),
    )
    return _ExperimentalHighOrderHeSingletData3D(
        problem,
        result.matrix,
        result.value,
        result.residual,
        result.iterations,
        result.converged,
        diagnostics,
    )
end

function _experimental_high_order_doside_he_singlet_energy(
    basis::MappedUniformBasis,
    coefficient_matrix::AbstractMatrix{<:Real};
    backend::Symbol = :numerical_reference,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
    krylovdim::Int = 32,
    maxiter::Int = 32,
    tol::Real = 1.0e-8,
)
    return _experimental_high_order_doside_he_singlet_data(
        basis,
        coefficient_matrix;
        backend = backend,
        expansion = expansion,
        Z = Z,
        krylovdim = krylovdim,
        maxiter = maxiter,
        tol = tol,
    ).ground_energy
end

function _experimental_high_order_doside_he_singlet_energy(
    stack::ExperimentalHighOrderDosideStack3D;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
    krylovdim::Int = 32,
    maxiter::Int = 32,
    tol::Real = 1.0e-8,
)
    return _experimental_high_order_doside_he_singlet_data(
        stack;
        expansion = expansion,
        Z = Z,
        krylovdim = krylovdim,
        maxiter = maxiter,
        tol = tol,
    ).ground_energy
end
