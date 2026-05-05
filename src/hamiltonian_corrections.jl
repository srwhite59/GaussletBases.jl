"""
    EGOIDensityDensityCorrectionResult

Matrix-level result returned by [`egoi_density_density_correction`](@ref).
`interaction_matrix` is the corrected density-density interaction matrix,
`interaction_delta` is the symmetric correction added to the input matrix, and
`diagnostics` records target-space residuals, correction sizes, and product
matrix conditioning.
"""
struct EGOIDensityDensityCorrectionResult
    interaction_matrix::Matrix{Float64}
    interaction_delta::Matrix{Float64}
    diagnostics::NamedTuple
end

"""
    StationaryFockCorrectionResult

Matrix-level result returned by [`stationary_fock_one_body_correction`](@ref).
`fock_matrix` is `F + one_body_delta`; the same symmetric `one_body_delta` is
intended to be added to the one-body Hamiltonian that generated `F`.
"""
struct StationaryFockCorrectionResult
    fock_matrix::Matrix{Float64}
    one_body_delta::Matrix{Float64}
    diagnostics::NamedTuple
end

"""
    HamiltonianCorrectionResult

Combined matrix-level result returned by
[`egoi_stationary_hamiltonian_correction`](@ref) and the ordinary convenience
adapter. `one_body_hamiltonian` and `interaction_matrix` are the final matrices;
`one_body_delta` and `interaction_delta` are final-minus-initial matrix deltas.
"""
struct HamiltonianCorrectionResult
    one_body_hamiltonian::Matrix{Float64}
    interaction_matrix::Matrix{Float64}
    one_body_delta::Matrix{Float64}
    interaction_delta::Matrix{Float64}
    diagnostics::NamedTuple
end

"""
    OrdinaryProjectedHamiltonianCorrectionTarget

Convenience target object produced by
[`ordinary_cartesian_projected_gaussian_target`](@ref). It carries Gaussian
target coefficients, their projection into an orthonormal ordinary working
basis, the dense Gaussian exact target Coulomb matrix, occupations, selected
Gaussian column indices, and projection diagnostics.
"""
struct OrdinaryProjectedHamiltonianCorrectionTarget
    projected_orbitals::Matrix{Float64}
    gaussian_coefficients::Matrix{Float64}
    exact_target::Matrix{Float64}
    occupations::Vector{Float64}
    indices::Vector{Int}
    diagnostics::NamedTuple
end

function _hc_symmetrize(matrix::AbstractMatrix{<:Real})
    dense = Matrix{Float64}(matrix)
    return Matrix{Float64}(0.5 .* (dense .+ transpose(dense)))
end

function _hc_max_abs(values)
    length(values) == 0 && return 0.0
    return Float64(maximum(abs, values))
end

function _hc_rms_abs(values)
    length(values) == 0 && return 0.0
    return Float64(sqrt(sum(abs2, values) / length(values)))
end

function _hc_fro_norm(values)
    return Float64(norm(values))
end

function _hc_relative(numerator::Real, denominator_values)
    denominator = _hc_fro_norm(denominator_values)
    denominator == 0.0 && return (Float64(numerator) == 0.0 ? 0.0 : Inf)
    return Float64(numerator) / denominator
end

function _hc_matrix_delta_diagnostics(delta::AbstractMatrix{<:Real}, baseline::AbstractMatrix{<:Real})
    delta_matrix = Matrix{Float64}(delta)
    baseline_matrix = Matrix{Float64}(baseline)
    delta_max = _hc_max_abs(delta_matrix)
    delta_fro = _hc_fro_norm(delta_matrix)
    return (
        max = delta_max,
        rms = _hc_rms_abs(delta_matrix),
        fro = delta_fro,
        relative_max = _hc_max_abs(baseline_matrix) == 0.0 ?
                       (delta_max == 0.0 ? 0.0 : Inf) :
                       delta_max / _hc_max_abs(baseline_matrix),
        relative_fro = _hc_relative(delta_fro, baseline_matrix),
    )
end

function _hc_residual_diagnostics(residual::AbstractMatrix{<:Real})
    residual_matrix = Matrix{Float64}(residual)
    return (
        fro = _hc_fro_norm(residual_matrix),
        max = _hc_max_abs(residual_matrix),
    )
end

function _hc_regularized_pseudoinverse(
    matrix::AbstractMatrix{<:Real};
    regularization::Real,
)
    regularization_value = Float64(regularization)
    regularization_value >= 0.0 ||
        throw(ArgumentError("regularization must be nonnegative"))
    decomposition = svd(Matrix{Float64}(matrix))
    singular_values = Float64[Float64(value) for value in decomposition.S]
    scale_cutoff =
        isempty(singular_values) ? 0.0 : max(size(matrix)...) * eps(Float64) * maximum(singular_values)
    cutoff = max(regularization_value, scale_cutoff)
    rank = count(value -> value > cutoff, singular_values)
    inverse_values = [
        value > cutoff ? inv(value) : 0.0 for value in singular_values
    ]
    inverse_matrix =
        decomposition.V * Diagonal(inverse_values) * transpose(decomposition.U)
    return inverse_matrix, singular_values, rank
end

function _hc_qr_orthonormalize_columns(matrix::AbstractMatrix{<:Real})
    dense = Matrix{Float64}(matrix)
    column_count = size(dense, 2)
    column_count > 0 || throw(ArgumentError("target orbital matrix must have at least one column"))
    size(dense, 1) >= column_count || throw(
        ArgumentError("target orbital matrix must have at least as many rows as columns"),
    )
    singular_values = Float64[Float64(value) for value in svdvals(dense)]
    rank_cutoff = max(size(dense)...) * eps(Float64) * maximum(singular_values)
    rank = count(value -> value > rank_cutoff, singular_values)
    rank == column_count || throw(
        ArgumentError("target orbital columns must be linearly independent for QR orthonormalization"),
    )
    qr_factor = qr(dense)
    q = Matrix(qr_factor.Q)[:, 1:column_count]
    return q, singular_values, rank
end

"""
    egoi_target_product_matrix(Qtarget)

Return the density-product matrix used by the EGOI correction. `Qtarget`
columns are target orbitals represented in an orthonormal working basis. Output
columns are pointwise products `Qtarget[:, i] .* Qtarget[:, j]`, with target
pair index `(j - 1) * ntarget + i`.
"""
function egoi_target_product_matrix(Qtarget::AbstractMatrix{<:Real})
    q = Matrix{Float64}(Qtarget)
    basis_count, target_count = size(q)
    product = zeros(Float64, basis_count, target_count^2)
    for j in 1:target_count, i in 1:target_count
        product[:, (j - 1) * target_count + i] .= q[:, i] .* q[:, j]
    end
    return product
end

function _egoi_gaussian_pair_transform(Ctarget::AbstractMatrix{<:Real})
    coefficients = Matrix{Float64}(Ctarget)
    gaussian_count, target_count = size(coefficients)
    transform = zeros(Float64, gaussian_count^2, target_count^2)
    for j in 1:target_count, i in 1:target_count
        column = (j - 1) * target_count + i
        for q in 1:gaussian_count, p in 1:gaussian_count
            row = gaussian_coulomb_pair_index(p, q, gaussian_count)
            transform[row, column] = coefficients[p, i] * coefficients[q, j]
        end
    end
    return transform
end

"""
    egoi_target_coulomb_matrix(pair_coulomb, Ctarget)

Transform a dense Gaussian pair Coulomb matrix into target-orbital pair space.
`pair_coulomb` uses the [`gaussian_coulomb_pair_index`](@ref) convention, and
`Ctarget` columns are target orbitals in that Gaussian basis.
"""
function egoi_target_coulomb_matrix(
    pair_coulomb::AbstractMatrix{<:Real},
    Ctarget::AbstractMatrix{<:Real},
)
    pair_matrix = _hc_symmetrize(pair_coulomb)
    pair_count = size(pair_matrix, 1)
    size(pair_matrix, 2) == pair_count ||
        throw(DimensionMismatch("pair_coulomb must be square"))
    gaussian_count = round(Int, sqrt(pair_count))
    gaussian_count^2 == pair_count || throw(
        DimensionMismatch("pair_coulomb dimension must be a perfect square Gaussian pair count"),
    )
    size(Ctarget, 1) == gaussian_count || throw(
        DimensionMismatch("Ctarget row count must match the Gaussian orbital count"),
    )
    transform = _egoi_gaussian_pair_transform(Ctarget)
    return _hc_symmetrize(transpose(transform) * pair_matrix * transform)
end

"""
    egoi_density_density_correction(V, Qtarget, exact_target; regularization=1e-18)

Compute a symmetric density-density correction `ΔV` so the target-orbital pair
space reproduces `exact_target` as closely as possible:
`P' * (V + ΔV) * P ≈ exact_target`, where `P =
egoi_target_product_matrix(Qtarget)`.

The correction is confined to the span of `P` using a regularized pseudoinverse,
so rank-deficient product spaces are handled explicitly and reported in
diagnostics.
"""
function egoi_density_density_correction(
    V::AbstractMatrix{<:Real},
    Qtarget::AbstractMatrix{<:Real},
    exact_target::AbstractMatrix{<:Real};
    regularization::Real = 1.0e-18,
)
    interaction = _hc_symmetrize(V)
    size(interaction, 1) == size(Qtarget, 1) || throw(
        DimensionMismatch("V dimension must match the target working-basis row count"),
    )
    product = egoi_target_product_matrix(Qtarget)
    target = _hc_symmetrize(exact_target)
    size(target) == (size(product, 2), size(product, 2)) || throw(
        DimensionMismatch("exact_target must have one row/column per target orbital pair"),
    )
    initial_target = _hc_symmetrize(transpose(product) * interaction * product)
    target_residual_before = initial_target - target
    product_inverse, singular_values, product_rank = _hc_regularized_pseudoinverse(
        product;
        regularization,
    )
    delta_target = target - initial_target
    delta_v = _hc_symmetrize(transpose(product_inverse) * delta_target * product_inverse)
    corrected = _hc_symmetrize(interaction + delta_v)
    corrected_target = _hc_symmetrize(transpose(product) * corrected * product)
    target_residual_after = corrected_target - target
    delta_diagnostics = _hc_matrix_delta_diagnostics(delta_v, interaction)
    diagnostics = (
        target_residual_fro_before = _hc_fro_norm(target_residual_before),
        target_residual_max_before = _hc_max_abs(target_residual_before),
        target_residual_fro_after = _hc_fro_norm(target_residual_after),
        target_residual_max_after = _hc_max_abs(target_residual_after),
        delta_v_max = delta_diagnostics.max,
        delta_v_rms = delta_diagnostics.rms,
        delta_v_fro = delta_diagnostics.fro,
        delta_v_relative_max = delta_diagnostics.relative_max,
        delta_v_relative_fro = delta_diagnostics.relative_fro,
        product_rank = product_rank,
        product_singular_values = singular_values,
        product_min_singular_value = isempty(singular_values) ? 0.0 : minimum(singular_values),
        product_max_singular_value = isempty(singular_values) ? 0.0 : maximum(singular_values),
        regularization = Float64(regularization),
        effective_singular_cutoff = max(
            Float64(regularization),
            isempty(singular_values) ? 0.0 : max(size(product)...) * eps(Float64) * maximum(singular_values),
        ),
    )
    return EGOIDensityDensityCorrectionResult(corrected, delta_v, diagnostics)
end

"""
    projected_orbital_density(Qtarget, occupations)

Return the density vector `sum_a occupations[a] * abs2(Qtarget[:, a])` in an
orthonormal working basis.
"""
function projected_orbital_density(
    Qtarget::AbstractMatrix{<:Real},
    occupations::AbstractVector{<:Real},
)
    q = Matrix{Float64}(Qtarget)
    occupation_values = Float64[Float64(value) for value in occupations]
    size(q, 2) == length(occupation_values) || throw(
        DimensionMismatch("occupation count must match target orbital count"),
    )
    return Float64[
        sum(occupation_values[column] * abs2(q[row, column]) for column in axes(q, 2))
        for row in axes(q, 1)
    ]
end

"""
    density_density_restricted_fock(H, V, D)

Build the density-density restricted Fock matrix `H + Diagonal(V * D)` in the
repo two-index interaction convention. `D` is a working-basis density vector,
typically from [`projected_orbital_density`](@ref).
"""
function density_density_restricted_fock(
    H::AbstractMatrix{<:Real},
    V::AbstractMatrix{<:Real},
    D::AbstractVector{<:Real},
)
    h = _hc_symmetrize(H)
    interaction = _hc_symmetrize(V)
    density = Float64[Float64(value) for value in D]
    size(h, 1) == size(interaction, 1) == length(density) || throw(
        DimensionMismatch("H, V, and D dimensions must match"),
    )
    return _hc_symmetrize(h + Diagonal(interaction * density))
end

"""
    occupied_virtual_fock_residual(F, Qocc)

Return `(I - P) * F * P` in the Euclidean working-basis metric, where columns
of `Qocc` are QR-orthonormalized before forming `P = Q * Q'`.
"""
function occupied_virtual_fock_residual(
    F::AbstractMatrix{<:Real},
    Qocc::AbstractMatrix{<:Real},
)
    fock = _hc_symmetrize(F)
    q, _singular_values, _rank = _hc_qr_orthonormalize_columns(Qocc)
    size(fock, 1) == size(q, 1) || throw(
        DimensionMismatch("F dimension must match occupied orbital row count"),
    )
    projector = q * transpose(q)
    return Matrix{Float64}(fock * projector - projector * fock * projector)
end

"""
    stationary_fock_one_body_correction(F, Qocc)

Return a symmetric one-body correction that cancels the occupied-virtual Fock
residual for the QR-orthonormalized occupied subspace. This routine assumes an
orthonormal Euclidean working basis and intentionally has no general
metric-overlap API.
"""
function stationary_fock_one_body_correction(
    F::AbstractMatrix{<:Real},
    Qocc::AbstractMatrix{<:Real};
    reference_matrix::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
)
    fock = _hc_symmetrize(F)
    q, singular_values, target_rank = _hc_qr_orthonormalize_columns(Qocc)
    size(fock, 1) == size(q, 1) || throw(
        DimensionMismatch("F dimension must match occupied orbital row count"),
    )
    projector = q * transpose(q)
    residual_before = Matrix{Float64}(fock * projector - projector * fock * projector)
    delta_h = _hc_symmetrize(-(residual_before + transpose(residual_before)))
    corrected_fock = _hc_symmetrize(fock + delta_h)
    residual_after = Matrix{Float64}(
        corrected_fock * projector - projector * corrected_fock * projector,
    )
    baseline = reference_matrix === nothing ? fock : Matrix{Float64}(reference_matrix)
    delta_diagnostics = _hc_matrix_delta_diagnostics(delta_h, baseline)
    diagnostics = (
        occupied_virtual_residual_fro_before = _hc_fro_norm(residual_before),
        occupied_virtual_residual_max_before = _hc_max_abs(residual_before),
        occupied_virtual_residual_fro_after = _hc_fro_norm(residual_after),
        occupied_virtual_residual_max_after = _hc_max_abs(residual_after),
        delta_h_max = delta_diagnostics.max,
        delta_h_rms = delta_diagnostics.rms,
        delta_h_fro = delta_diagnostics.fro,
        delta_h_relative_max = delta_diagnostics.relative_max,
        delta_h_relative_fro = delta_diagnostics.relative_fro,
        target_rank = target_rank,
        target_orthonormality_error = norm(transpose(q) * q - I, Inf),
        target_singular_values = singular_values,
        target_min_singular_value = isempty(singular_values) ? 0.0 : minimum(singular_values),
        target_max_singular_value = isempty(singular_values) ? 0.0 : maximum(singular_values),
    )
    return StationaryFockCorrectionResult(corrected_fock, delta_h, diagnostics)
end

function _hc_occupied_target_columns(
    Qtarget::AbstractMatrix{<:Real},
    occupations::AbstractVector{<:Real},
)
    occupation_values = Float64[Float64(value) for value in occupations]
    occupied_indices = Int[
        index for (index, occupation) in pairs(occupation_values) if occupation > 0.0
    ]
    isempty(occupied_indices) &&
        throw(ArgumentError("at least one occupation must be positive for stationary correction"))
    return Matrix{Float64}(Qtarget)[:, occupied_indices]
end

"""
    egoi_stationary_hamiltonian_correction(H, V, Qtarget, exact_target, occupations;
        include_egoi=true, include_stationary=true, regularization=1e-18)

Compose the EGOI density-density interaction correction and the stationary-Fock
one-body correction on dense matrices in an orthonormal working basis.
"""
function egoi_stationary_hamiltonian_correction(
    H::AbstractMatrix{<:Real},
    V::AbstractMatrix{<:Real},
    Qtarget::AbstractMatrix{<:Real},
    exact_target::AbstractMatrix{<:Real},
    occupations::AbstractVector{<:Real};
    include_egoi::Bool = true,
    include_stationary::Bool = true,
    regularization::Real = 1.0e-18,
)
    h = _hc_symmetrize(H)
    interaction = _hc_symmetrize(V)
    size(h) == size(interaction) ||
        throw(DimensionMismatch("H and V must have matching dimensions"))
    size(h, 1) == size(Qtarget, 1) ||
        throw(DimensionMismatch("Qtarget row count must match H/V dimension"))
    length(occupations) == size(Qtarget, 2) ||
        throw(DimensionMismatch("occupation count must match target orbital count"))

    egoi_result =
        include_egoi ?
        egoi_density_density_correction(
            interaction,
            Qtarget,
            exact_target;
            regularization,
        ) :
        nothing
    corrected_v = egoi_result === nothing ? interaction : egoi_result.interaction_matrix
    interaction_delta =
        egoi_result === nothing ? zeros(Float64, size(interaction)) : egoi_result.interaction_delta

    density = projected_orbital_density(Qtarget, occupations)
    fock = density_density_restricted_fock(h, corrected_v, density)
    stationary_result =
        include_stationary ?
        stationary_fock_one_body_correction(
            fock,
            _hc_occupied_target_columns(Qtarget, occupations);
            reference_matrix = h,
        ) :
        nothing
    one_body_delta =
        stationary_result === nothing ? zeros(Float64, size(h)) : stationary_result.one_body_delta
    corrected_h = _hc_symmetrize(h + one_body_delta)
    diagnostics = (
        include_egoi = include_egoi,
        include_stationary = include_stationary,
        density = density,
        egoi = egoi_result === nothing ? nothing : egoi_result.diagnostics,
        stationary = stationary_result === nothing ? nothing : stationary_result.diagnostics,
    )
    return HamiltonianCorrectionResult(
        corrected_h,
        corrected_v,
        one_body_delta,
        interaction_delta,
        diagnostics,
    )
end

"""
    ordinary_cartesian_projected_gaussian_target(operators, supplement, coefficients;
        indices=nothing, occupations=nothing, expansion=coulomb_gaussian_expansion(doacc=false),
        max_orbitals=64)

Project selected Gaussian target orbitals into an ordinary Cartesian working
basis and build their dense Gaussian exact Coulomb target matrix. This is a
convenience adapter for the shared matrix correction layer, not a separate
ordinary-specific correction model.
"""
function ordinary_cartesian_projected_gaussian_target(
    operators::OrdinaryCartesianOperators3D,
    supplement,
    coefficients::AbstractMatrix{<:Real};
    indices = nothing,
    occupations = nothing,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    max_orbitals = 64,
)
    coefficient_matrix = Matrix{Float64}(coefficients)
    selected_indices =
        indices === nothing ? collect(axes(coefficient_matrix, 2)) : Int[Int(index) for index in indices]
    isempty(selected_indices) &&
        throw(ArgumentError("ordinary projected correction target requires at least one selected index"))
    all(index -> 1 <= index <= size(coefficient_matrix, 2), selected_indices) || throw(
        ArgumentError("target coefficient indices must be valid coefficient matrix columns"),
    )
    target_coefficients = coefficient_matrix[:, selected_indices]
    overlap = gto_overlap_matrix(operators, supplement)
    size(overlap, 2) == size(target_coefficients, 1) || throw(
        DimensionMismatch("Gaussian target coefficient row count must match supplement orbital count"),
    )
    projected_orbitals = Matrix{Float64}(overlap * target_coefficients)
    occupation_values =
        occupations === nothing ?
        ones(Float64, length(selected_indices)) :
        Float64[Float64(value) for value in occupations]
    length(occupation_values) == length(selected_indices) || throw(
        DimensionMismatch("occupation count must match selected target count"),
    )
    pair_coulomb = gaussian_coulomb_pair_matrix(
        supplement;
        expansion,
        max_orbitals,
    )
    exact_target = egoi_target_coulomb_matrix(pair_coulomb, target_coefficients)
    projected_overlap = transpose(projected_orbitals) * projected_orbitals
    diagnostics = (
        selected_indices = Tuple(selected_indices),
        gaussian_orbital_count = size(coefficient_matrix, 1),
        target_count = length(selected_indices),
        projected_overlap_error = norm(projected_overlap - I, Inf),
        projected_column_norms = Float64[norm(view(projected_orbitals, :, i)) for i in axes(projected_orbitals, 2)],
        dense_pair_dimension = size(pair_coulomb, 1),
    )
    return OrdinaryProjectedHamiltonianCorrectionTarget(
        projected_orbitals,
        target_coefficients,
        exact_target,
        occupation_values,
        selected_indices,
        diagnostics,
    )
end

"""
    ordinary_cartesian_egoi_stationary_correction(operators; target,
        nuclear_charges=nothing, include_egoi=true, include_stationary=true,
        regularization=1e-18, overlap_tol=1e-8)

Apply the shared matrix-level EGOI/stationary correction machinery to an
ordinary Cartesian QW operator payload. The result is a
[`HamiltonianCorrectionResult`](@ref), not a new operator payload.
"""
function ordinary_cartesian_egoi_stationary_correction(
    operators::OrdinaryCartesianOperators3D;
    target::OrdinaryProjectedHamiltonianCorrectionTarget,
    nuclear_charges = nothing,
    include_egoi::Bool = true,
    include_stationary::Bool = true,
    regularization::Real = 1.0e-18,
    overlap_tol::Real = 1.0e-8,
)
    overlap_error = norm(operators.overlap - I, Inf)
    overlap_error <= Float64(overlap_tol) || throw(
        ArgumentError(
            "ordinary_cartesian_egoi_stationary_correction requires an orthonormal final ordinary basis; got overlap error $(overlap_error)",
        ),
    )
    charges = nuclear_charges === nothing ? operators.nuclear_charges : nuclear_charges
    h = assembled_one_body_hamiltonian(operators; nuclear_charges = charges)
    result = egoi_stationary_hamiltonian_correction(
        h,
        operators.interaction_matrix,
        target.projected_orbitals,
        target.exact_target,
        target.occupations;
        include_egoi,
        include_stationary,
        regularization,
    )
    diagnostics = merge(
        result.diagnostics,
        (
            ordinary_adapter = (
                overlap_error = overlap_error,
                nuclear_charges = charges === nothing ? nothing : Tuple(Float64.(collect(charges))),
            ),
            target = target.diagnostics,
        ),
    )
    return HamiltonianCorrectionResult(
        result.one_body_hamiltonian,
        result.interaction_matrix,
        result.one_body_delta,
        result.interaction_delta,
        diagnostics,
    )
end
