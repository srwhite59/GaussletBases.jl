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

struct _ExperimentalHighOrderHeSingletData3D{P}
    problem::P
    ground_matrix::Matrix{Float64}
    ground_energy::Float64
    residual::Float64
    iterations::Int
    converged::Bool
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
       name === :residual || name === :iterations || name === :converged
        return getfield(data, name)
    end
    return getproperty(getfield(data, :problem), name)
end

function Base.propertynames(data::_ExperimentalHighOrderHeSingletData3D, private::Bool = false)
    names = (:problem, :ground_matrix, :ground_energy, :residual, :iterations, :converged)
    return (names..., propertynames(getfield(data, :problem), private)...)
end

function _experimental_high_order_projected_one_body_data(
    basis::MappedUniformBasis,
    coefficient_matrix::AbstractMatrix{<:Real};
    backend::Symbol = :numerical_reference,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
)
    coefficients = Matrix{Float64}(coefficient_matrix)
    one_body = mapped_ordinary_one_body_operators(
        basis;
        exponents = expansion.exponents,
        backend = backend,
    )
    parent_overlap, parent_hamiltonian = _mapped_cartesian_one_body_matrix(one_body, expansion; Z = Z)
    projected_overlap = _symmetrize_ida_matrix(transpose(coefficients) * parent_overlap * coefficients)
    projected_hamiltonian = _symmetrize_ida_matrix(transpose(coefficients) * parent_hamiltonian * coefficients)
    return _ExperimentalHighOrderProjectedOneBodyData3D(
        basis,
        backend,
        expansion,
        Float64(Z),
        one_body,
        coefficients,
        parent_overlap,
        parent_hamiltonian,
        projected_overlap,
        projected_hamiltonian,
    )
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
    return _ExperimentalHighOrderHeSingletData3D(
        problem,
        result.matrix,
        result.value,
        result.residual,
        result.iterations,
        result.converged,
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
    return _experimental_high_order_doside_he_singlet_data(
        stack.parent_basis,
        stack.coefficient_matrix;
        backend = stack.backend,
        expansion = expansion,
        Z = Z,
        krylovdim = krylovdim,
        maxiter = maxiter,
        tol = tol,
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
