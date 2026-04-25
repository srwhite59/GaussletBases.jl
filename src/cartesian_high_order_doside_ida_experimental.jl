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
