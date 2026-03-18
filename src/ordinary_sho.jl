function _gaussian_x2(a::Gaussian, b::Gaussian)
    sigma2 = a.width^2 + b.width^2
    weighted_center =
        (a.center_value * b.width^2 + b.center_value * a.width^2) / sigma2
    variance = a.width^2 * b.width^2 / sigma2
    return (weighted_center^2 + variance) * _gaussian_overlap(a, b)
end

function _primitive_x2_matrix(
    set::PrimitiveSet1D,
    ::_NumericalPrimitiveMatrixBackend;
    h = nothing,
)
    xlo, xhi = _primitive_set_bounds(set)
    h_try = h === nothing ? _primitive_matrix_start_h(set) : Float64(h)
    h_try > 0.0 || throw(ArgumentError("numerical primitive x2 matrix requires h > 0"))

    previous = nothing
    current = nothing
    for _ in 1:_PRIMITIVE_MATRIX_MAXITER
        points, weights = _make_midpoint_grid(xlo, xhi, h_try)
        values = _primitive_sample_matrix(set, points)
        current = _symmetrize_primitive_matrix(
            transpose(values) * (((weights .* (points .^ 2))) .* values),
        )
        if previous !== nothing && norm(current - previous, Inf) <= _PRIMITIVE_MATRIX_TOL
            return current
        end
        previous = current
        h_try /= 2.0
    end
    return current
end

function _primitive_x2_matrix(set::PrimitiveSet1D, ::_AnalyticPrimitiveMatrixBackend)
    matrix = zeros(Float64, length(set), length(set))
    for a in 1:length(set)
        pa = primitives(set)[a]
        for b in a:length(set)
            pb = primitives(set)[b]
            value_ab = _gaussian_x2(pa, pb)
            matrix[a, b] = value_ab
            matrix[b, a] = value_ab
        end
    end
    return matrix
end

function _x2_matrix(set::PrimitiveSet1D)
    return _primitive_x2_matrix(set, _select_primitive_matrix_backend(set))
end

function _x2_matrix(basis_like)
    return contract_primitive_matrix(basis_like, _x2_matrix(primitive_set(basis_like)))
end

function _ordinary_position_matrix(basis::MappedUniformBasis)
    representation = basis_representation(basis; operators = (:position,))
    return Matrix{Float64}(representation.basis_matrices.position)
end

function _ordinary_position_matrix(basis_like)
    return Matrix{Float64}(position_matrix(basis_like))
end

function _ordinary_sho_layer(basis::MappedUniformBasis, backend::Symbol)
    return _mapped_ordinary_backend_layer(basis, backend)
end

function _ordinary_sho_layer(basis::HybridMappedOrdinaryBasis1D, backend::Symbol)
    backend == basis.backend ||
        throw(ArgumentError("hybrid mapped ordinary basis backend does not match the requested SHO backend"))
    return basis
end

function _ordinary_sho_core(operators::MappedOrdinaryOneBody1D)
    basis = operators.basis
    layer = _ordinary_sho_layer(basis, operators.backend)
    overlap = Matrix{Float64}(operators.overlap)
    kinetic = Matrix{Float64}(operators.kinetic)
    position = _ordinary_position_matrix(layer)
    x2 = Matrix{Float64}(_x2_matrix(layer))
    return (
        overlap = overlap,
        kinetic = kinetic,
        position = position,
        x2 = x2,
        backend = operators.backend,
    )
end

function _sho_exact_energies(omega::Real, nev::Integer)
    omega_value = Float64(omega)
    return Float64[omega_value * (n + 0.5) for n in 0:(nev - 1)]
end

"""
    ordinary_sho_hamiltonian(
        basis::MappedUniformBasis;
        omega = 1.0,
        center = 0.0,
        backend = :pgdg_localized_experimental,
    )
    ordinary_sho_hamiltonian(
        basis::HybridMappedOrdinaryBasis1D;
        omega = 1.0,
        center = 0.0,
    )
    ordinary_sho_hamiltonian(
        operators::MappedOrdinaryOneBody1D;
        omega = 1.0,
        center = 0.0,
    )

Build the one-dimensional harmonic-oscillator Hamiltonian

`H = T + 0.5 * omega^2 * (x - center)^2`

for the ordinary mapped branch. The returned named tuple contains the overlap,
kinetic, position, displacement-squared, potential, and total Hamiltonian
matrices in the current basis representation.
"""
function ordinary_sho_hamiltonian(
    operators::MappedOrdinaryOneBody1D;
    omega::Real = 1.0,
    center::Real = 0.0,
)
    core = _ordinary_sho_core(operators)
    overlap = core.overlap
    kinetic = core.kinetic
    position = core.position
    x2 = core.x2
    displacement2 = x2 .- 2.0 * Float64(center) .* position .+ Float64(center)^2 .* overlap
    omega_value = Float64(omega)
    potential = 0.5 * omega_value^2 .* displacement2
    hamiltonian = kinetic + potential
    return (
        overlap = overlap,
        kinetic = kinetic,
        position = position,
        x2 = x2,
        displacement2 = displacement2,
        potential = potential,
        hamiltonian = hamiltonian,
        backend = core.backend,
        center = Float64(center),
        omega = omega_value,
    )
end

function ordinary_sho_hamiltonian(
    basis::MappedUniformBasis;
    omega::Real = 1.0,
    center::Real = 0.0,
    backend::Symbol = :pgdg_localized_experimental,
)
    operators = mapped_ordinary_one_body_operators(basis; backend = backend)
    return ordinary_sho_hamiltonian(operators; omega = omega, center = center)
end

function ordinary_sho_hamiltonian(
    basis::HybridMappedOrdinaryBasis1D;
    omega::Real = 1.0,
    center::Real = 0.0,
)
    operators = mapped_ordinary_one_body_operators(basis)
    return ordinary_sho_hamiltonian(operators; omega = omega, center = center)
end

"""
    ordinary_sho_spectrum(args...; omega = 1.0, center = 0.0, nev = 3, backend = ...)

Return a small spectral comparison package for the one-dimensional harmonic
oscillator on the ordinary mapped branch. The result includes the lowest
eigenvalues, the exact SHO values, and the ground-state expectations of
`T` and `(x - center)^2` in the orthonormalized basis.
"""
function ordinary_sho_spectrum(
    operators::MappedOrdinaryOneBody1D;
    omega::Real = 1.0,
    center::Real = 0.0,
    nev::Integer = 3,
)
    return _ordinary_sho_spectrum_from_core(
        _ordinary_sho_core(operators);
        omega = omega,
        center = center,
        nev = nev,
    )
end

function _ordinary_sho_spectrum_from_core(
    core::NamedTuple;
    omega::Real,
    center::Real,
    nev::Integer,
)
    nev > 0 || throw(ArgumentError("ordinary_sho_spectrum requires nev > 0"))
    overlap = core.overlap
    kinetic = core.kinetic
    position = core.position
    x2 = core.x2
    displacement2 = x2 .- 2.0 * Float64(center) .* position .+ Float64(center)^2 .* overlap
    potential = 0.5 * Float64(omega)^2 .* displacement2
    hamiltonian = kinetic + potential
    transformed = _orthonormalize_cartesian_1d(
        overlap,
        [hamiltonian, kinetic, displacement2],
    )
    hamiltonian_orth = transformed[1]
    kinetic_orth = transformed[2]
    displacement2_orth = transformed[3]
    decomposition = eigen(Hermitian(hamiltonian_orth))
    values = Float64[Float64(value) for value in decomposition.values]
    vectors = decomposition.vectors
    ground = view(vectors, :, 1)
    return (
        backend = core.backend,
        omega = Float64(omega),
        center = Float64(center),
        eigenvalues = values[1:min(nev, length(values))],
        exact = _sho_exact_energies(omega, min(nev, length(values))),
        kinetic_expectation = Float64(real(dot(ground, kinetic_orth * ground))),
        displacement2_expectation = Float64(real(dot(ground, displacement2_orth * ground))),
    )
end

function ordinary_sho_spectrum(
    basis::MappedUniformBasis;
    omega::Real = 1.0,
    center::Real = 0.0,
    nev::Integer = 3,
    backend::Symbol = :pgdg_localized_experimental,
)
    operators = mapped_ordinary_one_body_operators(basis; backend = backend)
    return ordinary_sho_spectrum(operators; omega = omega, center = center, nev = nev)
end

function ordinary_sho_spectrum(
    basis::HybridMappedOrdinaryBasis1D;
    omega::Real = 1.0,
    center::Real = 0.0,
    nev::Integer = 3,
)
    operators = mapped_ordinary_one_body_operators(basis)
    return ordinary_sho_spectrum(operators; omega = omega, center = center, nev = nev)
end
