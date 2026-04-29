# Post-assembly Qiu-White-style correction diagnostics for ordinary Cartesian
# operators. The public surfaces below expose only projector-based one-body
# corrections, with ESOI as explicit opt-in. The broad mode-carrying
# HydrogenicCoreCorrectionSpec remains internal; :local_exact is diagnostic
# only, and a first-order local variant remains unsupported.

"""
    HydrogenicCoreProjectorCorrectionSpec(; Z, nucleus=(0.0, 0.0, 0.0),
        include_esoi=false, local_orbital_index=nothing)

Public post-assembly hydrogenic core correction specification for
[`OrdinaryCartesianOperators3D`](@ref). The one-body correction is the
projector reference route: the lowest finite-basis hydrogenic core orbital is
shifted to the exact `-Z^2/2` eigenvalue. Set `include_esoi=true` to also apply
the local ESOI-style two-body correction that calibrates the same core orbital's
`1s^2` Coulomb scalar to `5Z/8`.

The correction is applied after operator assembly. Corrected payloads keep the
total one-body and interaction matrices as authoritative and deliberately drop
decomposition sidecars (`kinetic_one_body === nothing`,
`nuclear_one_body_by_center === nothing`, `nuclear_term_storage == :total_only`).
`local_orbital_index` is only an optional override for the local basis function
used by ESOI; by default the nearest non-residual local orbital to `nucleus` is
selected.
"""
struct HydrogenicCoreProjectorCorrectionSpec
    Z::Float64
    nucleus::NTuple{3,Float64}
    include_esoi::Bool
    local_orbital_index::Union{Nothing,Int}
end

function HydrogenicCoreProjectorCorrectionSpec(;
    Z::Real,
    nucleus = (0.0, 0.0, 0.0),
    include_esoi::Bool = false,
    local_orbital_index::Union{Nothing,Integer} = nothing,
)
    return HydrogenicCoreProjectorCorrectionSpec(
        Float64(Z),
        _ordinary_cartesian_correction_nucleus_tuple(nucleus),
        include_esoi,
        local_orbital_index === nothing ? nothing : Int(local_orbital_index),
    )
end

"""
    HydrogenicCoreBranchCorrectionSpec(; Z, nucleus=(0.0, 0.0, 0.0),
        include_esoi=false, orbital_selector=:localized_lowest,
        local_orbital_index=nothing, reference_nuclear_charges=nothing)

Public branch-level hydrogenic core correction specification. This is the
matrix-level companion to [`HydrogenicCoreProjectorCorrectionSpec`](@ref) for
counterpoise and branch Hamiltonian work: the one-body branch matrix is first
assembled with caller-provided nuclear charges, then the projector reference
correction is applied to the selected core orbital.

This public branch surface supports one correction or a collection of
corrections for one assembled branch matrix. The default
`orbital_selector = :localized_lowest` partitions final-basis orbitals by their
nearest carried nucleus, diagonalizes the selected center-local subspace, embeds
that local core orbital back into the full final basis, and applies the
projector correction to that embedded vector. Use
`orbital_selector = :global_lowest` as an explicit debug/reference selector
matching the original projector behavior; multi-correction calls require
localized selectors. Set `include_esoi=true` to apply the optional local ESOI
two-body scalar calibration to the same selected core orbital. By default,
localized branch corrections calibrate the one-body projector against a
center-isolated reference branch. Fragment application branches should include
corrections only for nuclei present in that branch. Pass
`reference_nuclear_charges` only to override the calibration branch explicitly;
the override must carry `Z` at the corrected center.
"""
struct HydrogenicCoreBranchCorrectionSpec
    Z::Float64
    nucleus::NTuple{3,Float64}
    include_esoi::Bool
    orbital_selector::Symbol
    local_orbital_index::Union{Nothing,Int}
    reference_nuclear_charges::Union{Nothing,Vector{Float64}}
end

function HydrogenicCoreBranchCorrectionSpec(;
    Z::Real,
    nucleus = (0.0, 0.0, 0.0),
    include_esoi::Bool = false,
    orbital_selector::Symbol = :localized_lowest,
    local_orbital_index::Union{Nothing,Integer} = nothing,
    reference_nuclear_charges = nothing,
)
    orbital_selector in (:localized_lowest, :global_lowest) || throw(
        ArgumentError(
            "unsupported orbital_selector :$(orbital_selector); supported branch selectors are :localized_lowest and :global_lowest",
        ),
    )
    return HydrogenicCoreBranchCorrectionSpec(
        Float64(Z),
        _ordinary_cartesian_correction_nucleus_tuple(nucleus),
        include_esoi,
        orbital_selector,
        local_orbital_index === nothing ? nothing : Int(local_orbital_index),
        reference_nuclear_charges === nothing ? nothing : Float64.(collect(reference_nuclear_charges)),
    )
end

abstract type AbstractOrdinaryCartesianCorrectionSpec end

struct HydrogenicCoreCorrectionSpec <: AbstractOrdinaryCartesianCorrectionSpec
    Z::Float64
    nucleus::NTuple{3,Float64}
    one_body_mode::Symbol
    two_body_mode::Symbol
    local_selection::Symbol
    local_orbital_index::Union{Nothing,Int}
end

function HydrogenicCoreCorrectionSpec(;
    Z::Real,
    nucleus = (0.0, 0.0, 0.0),
    one_body_mode::Symbol = :projector,
    two_body_mode::Symbol = :esoi_local,
    local_selection::Symbol = :nearest_nonresidual,
    local_orbital_index::Union{Nothing,Integer} = nothing,
)
    one_body_mode in (:none, :projector, :local_exact) || throw(
        ArgumentError(
            "unsupported one_body_mode :$(one_body_mode); supported internal staging modes are :none, :projector, and :local_exact",
        ),
    )
    two_body_mode in (:none, :esoi_local) || throw(
        ArgumentError("unsupported two_body_mode :$(two_body_mode); supported modes are :none and :esoi_local"),
    )
    local_selection in (:nearest_nonresidual, :explicit) || throw(
        ArgumentError("unsupported local_selection :$(local_selection); supported modes are :nearest_nonresidual and :explicit"),
    )
    if local_selection == :explicit && local_orbital_index === nothing
        throw(ArgumentError("local_selection = :explicit requires local_orbital_index"))
    end
    nucleus_tuple = _ordinary_cartesian_correction_nucleus_tuple(nucleus)
    return HydrogenicCoreCorrectionSpec(
        Float64(Z),
        nucleus_tuple,
        one_body_mode,
        two_body_mode,
        local_selection,
        local_orbital_index === nothing ? nothing : Int(local_orbital_index),
    )
end

"""
    OrdinaryCartesianCorrectionResult

Result returned by [`apply_ordinary_cartesian_corrections`](@ref). `operators`
contains the corrected [`OrdinaryCartesianOperators3D`](@ref), while
`diagnostics` records the selected local orbital, core eigenvalue shifts,
`1s^2` Coulomb scalar before/after values, and closed-shell proxy energies.
"""
struct OrdinaryCartesianCorrectionResult
    operators::OrdinaryCartesianOperators3D
    diagnostics::NamedTuple
end

"""
    OrdinaryCartesianBranchCorrectionResult

Matrix-level result returned by [`ordinary_cartesian_corrected_branch`](@ref).
`one_body_hamiltonian` is the corrected branch one-body matrix assembled from
the requested nuclear charges, `interaction_matrix` is the corresponding
interaction matrix after any optional ESOI corrections, `one_body_delta` and
`interaction_delta` are the final-minus-initial application-branch matrix
deltas, and `diagnostics` is a branch-level `NamedTuple` containing branch
nuclear charges, correction count, overlap error, corrected center indices, and
a tuple of per-correction diagnostic `NamedTuple`s. This result deliberately
does not claim to be a new [`OrdinaryCartesianOperators3D`](@ref) payload.
"""
struct OrdinaryCartesianBranchCorrectionResult
    one_body_hamiltonian::Matrix{Float64}
    interaction_matrix::Matrix{Float64}
    one_body_delta::Matrix{Float64}
    interaction_delta::Matrix{Float64}
    diagnostics::NamedTuple
end

function _ordinary_cartesian_correction_nucleus_tuple(nucleus)
    values = Tuple(nucleus)
    length(values) == 3 ||
        throw(ArgumentError("HydrogenicCoreCorrectionSpec nucleus must have three coordinates"))
    return (Float64(values[1]), Float64(values[2]), Float64(values[3]))
end

function _ordinary_cartesian_lowest_orbital(h::AbstractMatrix{<:Real})
    decomposition = eigen(Hermitian(Matrix{Float64}(h)))
    orbital = Vector{Float64}(decomposition.vectors[:, 1])
    orbital ./= norm(orbital)
    return Float64(decomposition.values[1]), orbital
end

function _ordinary_cartesian_one_body_expectation(
    h::AbstractMatrix{<:Real},
    orbital::AbstractVector{<:Real},
)
    normalized = Vector{Float64}(orbital)
    normalized ./= norm(normalized)
    return Float64(real(dot(normalized, Matrix{Float64}(h) * normalized)))
end

function _ordinary_cartesian_ida_coulomb_scalar(
    v::AbstractMatrix{<:Real},
    orbital::AbstractVector{<:Real},
)
    weights = abs2.(Vector{Float64}(orbital))
    norm2 = sum(weights)
    norm2 > 0.0 || throw(ArgumentError("orbital must have nonzero coefficient norm"))
    weights ./= norm2
    return Float64(real(dot(weights, Matrix{Float64}(v) * weights)))
end

function _ordinary_cartesian_closed_shell_proxy_energy(
    h::AbstractMatrix{<:Real},
    v::AbstractMatrix{<:Real},
    orbital::AbstractVector{<:Real},
)
    return 2.0 * _ordinary_cartesian_one_body_expectation(h, orbital) +
           _ordinary_cartesian_ida_coulomb_scalar(v, orbital)
end

function _ordinary_cartesian_select_local_orbital(
    operators::OrdinaryCartesianOperators3D,
    spec::HydrogenicCoreCorrectionSpec,
)
    if spec.local_orbital_index !== nothing
        index = spec.local_orbital_index
        1 <= index <= length(operators.orbital_data) ||
            throw(ArgumentError("local_orbital_index $(index) is outside the operator dimension"))
        orbital = operators.orbital_data[index]
        orbital.kind == :residual_gaussian && throw(
            ArgumentError("local_orbital_index $(index) selects a residual Gaussian, not a local gausslet/fixed orbital"),
        )
        return index, orbital, 0.0
    end

    best_index = 0
    best_radius2 = Inf
    best_orbital = nothing
    x0, y0, z0 = spec.nucleus
    for (index, orbital) in pairs(operators.orbital_data)
        orbital.kind == :residual_gaussian && continue
        radius2 = (orbital.x - x0)^2 + (orbital.y - y0)^2 + (orbital.z - z0)^2
        if radius2 < best_radius2
            best_index = index
            best_radius2 = radius2
            best_orbital = orbital
        end
    end
    best_index > 0 ||
        throw(ArgumentError("no non-residual local orbital found for hydrogenic/ESOI correction"))
    return best_index, best_orbital, sqrt(best_radius2)
end

function _ordinary_cartesian_projector_corrected_hamiltonian(
    h::AbstractMatrix{<:Real},
    Z::Real,
)
    initial_eigenvalue, initial_orbital = _ordinary_cartesian_lowest_orbital(h)
    target = -0.5 * Float64(Z)^2
    shift = target - initial_eigenvalue
    corrected_h = Matrix{Float64}(h) .+ shift .* (initial_orbital * transpose(initial_orbital))
    corrected_eigenvalue, corrected_orbital = _ordinary_cartesian_lowest_orbital(corrected_h)
    return (
        h = corrected_h,
        target = target,
        initial_eigenvalue = initial_eigenvalue,
        initial_orbital = initial_orbital,
        shift = shift,
        corrected_eigenvalue = corrected_eigenvalue,
        corrected_orbital = corrected_orbital,
    )
end

function _ordinary_cartesian_projector_corrected_hamiltonian(
    h::AbstractMatrix{<:Real},
    Z::Real,
    orbital::AbstractVector{<:Real},
)
    selected_orbital = Vector{Float64}(orbital)
    selected_orbital ./= norm(selected_orbital)
    initial_expectation = _ordinary_cartesian_one_body_expectation(h, selected_orbital)
    target = -0.5 * Float64(Z)^2
    shift = target - initial_expectation
    corrected_h = Matrix{Float64}(h) .+ shift .* (selected_orbital * transpose(selected_orbital))
    corrected_expectation = _ordinary_cartesian_one_body_expectation(corrected_h, selected_orbital)
    corrected_eigenvalue, corrected_global_orbital = _ordinary_cartesian_lowest_orbital(corrected_h)
    return (
        h = corrected_h,
        target = target,
        initial_expectation = initial_expectation,
        selected_orbital = selected_orbital,
        shift = shift,
        corrected_expectation = corrected_expectation,
        corrected_eigenvalue = corrected_eigenvalue,
        corrected_global_orbital = corrected_global_orbital,
    )
end

function _ordinary_cartesian_projector_delta_from_reference(
    application_h::AbstractMatrix{<:Real},
    reference_h::AbstractMatrix{<:Real},
    Z::Real,
    orbital::AbstractVector{<:Real},
)
    selected_orbital = Vector{Float64}(orbital)
    selected_orbital ./= norm(selected_orbital)
    target = -0.5 * Float64(Z)^2
    calibration_expectation =
        _ordinary_cartesian_one_body_expectation(reference_h, selected_orbital)
    application_initial_expectation =
        _ordinary_cartesian_one_body_expectation(application_h, selected_orbital)
    shift = target - calibration_expectation
    projector = selected_orbital * transpose(selected_orbital)
    delta_h = shift .* projector
    corrected_h = Matrix{Float64}(application_h) .+ delta_h
    application_corrected_expectation =
        _ordinary_cartesian_one_body_expectation(corrected_h, selected_orbital)
    corrected_eigenvalue, corrected_global_orbital = _ordinary_cartesian_lowest_orbital(corrected_h)
    return (
        h = corrected_h,
        delta_h = delta_h,
        target = target,
        calibration_expectation = calibration_expectation,
        application_initial_expectation = application_initial_expectation,
        application_corrected_expectation = application_corrected_expectation,
        selected_orbital = selected_orbital,
        shift = shift,
        corrected_eigenvalue = corrected_eigenvalue,
        corrected_global_orbital = corrected_global_orbital,
    )
end

function _ordinary_cartesian_lowest_with_diagonal_shift(
    h::AbstractMatrix{<:Real},
    local_index::Integer,
    shift::Real,
)
    shifted = Matrix{Float64}(h)
    shifted[local_index, local_index] += Float64(shift)
    eigenvalue, orbital = _ordinary_cartesian_lowest_orbital(shifted)
    return eigenvalue, orbital, shifted
end

function _ordinary_cartesian_local_exact_corrected_hamiltonian(
    h::AbstractMatrix{<:Real},
    Z::Real,
    local_index::Integer;
    eig_tol::Real = 1.0e-11,
    maxiter::Integer = 80,
)
    # Diagnostic-only local alternative to the projector reference: solve for a
    # single diagonal shift whose lowest hydrogenic eigenvalue is exact.
    initial_eigenvalue, initial_orbital = _ordinary_cartesian_lowest_orbital(h)
    target = -0.5 * Float64(Z)^2
    initial_residual = initial_eigenvalue - target
    if abs(initial_residual) <= Float64(eig_tol)
        return (
            h = Matrix{Float64}(h),
            target = target,
            initial_eigenvalue = initial_eigenvalue,
            initial_orbital = initial_orbital,
            shift = 0.0,
            corrected_eigenvalue = initial_eigenvalue,
            corrected_orbital = initial_orbital,
            bracket_iterations = 0,
            bisection_iterations = 0,
        )
    end

    local_weight = max(abs2(initial_orbital[local_index]), eps(Float64))
    trial_step = (target - initial_eigenvalue) / local_weight
    trial_step == 0.0 && (trial_step = initial_residual > 0.0 ? -1.0 : 1.0)

    lower_shift = 0.0
    upper_shift = 0.0
    lower_value = initial_residual
    upper_value = initial_residual
    bracket_iterations = 0
    if initial_residual > 0.0
        lower_shift = trial_step < 0.0 ? trial_step : -abs(trial_step)
        while true
            bracket_iterations += 1
            lower_value =
                _ordinary_cartesian_lowest_with_diagonal_shift(h, local_index, lower_shift)[1] -
                target
            lower_value <= 0.0 && break
            bracket_iterations < maxiter || throw(
                ErrorException("failed to bracket local_exact hydrogenic one-body correction below target"),
            )
            lower_shift *= 2.0
        end
    else
        upper_shift = trial_step > 0.0 ? trial_step : abs(trial_step)
        while true
            bracket_iterations += 1
            upper_value =
                _ordinary_cartesian_lowest_with_diagonal_shift(h, local_index, upper_shift)[1] -
                target
            upper_value >= 0.0 && break
            bracket_iterations < maxiter || throw(
                ErrorException("failed to bracket local_exact hydrogenic one-body correction above target"),
            )
            upper_shift *= 2.0
        end
    end

    corrected_shift = 0.0
    corrected_eigenvalue = initial_eigenvalue
    corrected_orbital = initial_orbital
    corrected_h = Matrix{Float64}(h)
    bisection_iterations = 0
    for iteration in 1:maxiter
        bisection_iterations = iteration
        midpoint = 0.5 * (lower_shift + upper_shift)
        eigenvalue, orbital, shifted =
            _ordinary_cartesian_lowest_with_diagonal_shift(h, local_index, midpoint)
        residual = eigenvalue - target
        corrected_shift = midpoint
        corrected_eigenvalue = eigenvalue
        corrected_orbital = orbital
        corrected_h = shifted
        abs(residual) <= Float64(eig_tol) && break
        if residual <= 0.0
            lower_shift = midpoint
            lower_value = residual
        else
            upper_shift = midpoint
            upper_value = residual
        end
        abs(upper_shift - lower_shift) <= max(1.0, abs(midpoint)) * 1.0e-13 && break
    end

    return (
        h = corrected_h,
        target = target,
        initial_eigenvalue = initial_eigenvalue,
        initial_orbital = initial_orbital,
        shift = corrected_shift,
        corrected_eigenvalue = corrected_eigenvalue,
        corrected_orbital = corrected_orbital,
        bracket_iterations = bracket_iterations,
        bisection_iterations = bisection_iterations,
        bracket_residuals = (lower = lower_value, upper = upper_value),
    )
end

function _ordinary_cartesian_esoi_corrected_interaction(
    v::AbstractMatrix{<:Real},
    orbital::AbstractVector{<:Real},
    Z::Real,
    local_index::Integer,
)
    weights = abs2.(Vector{Float64}(orbital))
    weights ./= sum(weights)
    density_weight = weights[local_index]
    denominator = density_weight^2
    denominator > 0.0 ||
        throw(ArgumentError("selected local orbital has zero density weight in the hydrogenic core orbital"))

    target = 5.0 * Float64(Z) / 8.0
    initial_j = _ordinary_cartesian_ida_coulomb_scalar(v, orbital)
    shift = (target - initial_j) / denominator
    corrected_v = Matrix{Float64}(v)
    corrected_v[local_index, local_index] += shift
    corrected_j = _ordinary_cartesian_ida_coulomb_scalar(corrected_v, orbital)
    return (
        v = corrected_v,
        target = target,
        initial_j = initial_j,
        corrected_j = corrected_j,
        shift = shift,
        scalar_delta = target - initial_j,
        density_weight = density_weight,
    )
end

function _ordinary_cartesian_operators_with_corrections(
    operators::OrdinaryCartesianOperators3D,
    one_body_hamiltonian::AbstractMatrix{<:Real},
    interaction_matrix::AbstractMatrix{<:Real},
)
    # Post-assembly corrections invalidate exact one-body decomposition sidecars:
    # the total matrices remain authoritative, while kinetic/nuclear split
    # provenance is deliberately downgraded.
    nuclear_charges =
        operators.nuclear_charges === nothing ? nothing : copy(operators.nuclear_charges)
    return OrdinaryCartesianOperators3D(
        operators.basis,
        operators.gaussian_data,
        operators.gausslet_backend,
        operators.interaction_treatment,
        operators.expansion,
        copy(operators.overlap),
        Matrix{Float64}(one_body_hamiltonian),
        Matrix{Float64}(interaction_matrix),
        copy(operators.orbital_data),
        operators.gausslet_count,
        operators.residual_count,
        copy(operators.raw_to_final),
        copy(operators.residual_centers),
        copy(operators.residual_widths),
        nuclear_charges,
        nothing,
        nothing,
        :total_only,
    )
end

function _ordinary_cartesian_corrected_matrices(
    operators::OrdinaryCartesianOperators3D,
    initial_h::AbstractMatrix{<:Real},
    initial_v::AbstractMatrix{<:Real},
    spec::HydrogenicCoreCorrectionSpec;
    overlap_error::Real,
    extra_diagnostics::NamedTuple = (;),
)
    initial_h = Matrix{Float64}(initial_h)
    initial_v = Matrix{Float64}(initial_v)
    initial_eigenvalue, initial_orbital = _ordinary_cartesian_lowest_orbital(initial_h)
    target_eigenvalue = -0.5 * spec.Z^2
    local_index, local_orbital, local_distance =
        _ordinary_cartesian_select_local_orbital(operators, spec)

    corrected_h = initial_h
    corrected_eigenvalue = initial_eigenvalue
    corrected_orbital = initial_orbital
    one_body_shift = 0.0
    local_exact_bracket_iterations = 0
    local_exact_bisection_iterations = 0
    if spec.one_body_mode == :projector
        one_body = _ordinary_cartesian_projector_corrected_hamiltonian(initial_h, spec.Z)
        corrected_h = one_body.h
        corrected_eigenvalue = one_body.corrected_eigenvalue
        corrected_orbital = one_body.corrected_orbital
        one_body_shift = one_body.shift
    elseif spec.one_body_mode == :none
        # Leave the one-body matrix untouched and use the uncorrected core
        # orbital for any requested two-body scalar calibration.
    elseif spec.one_body_mode == :local_exact
        one_body =
            _ordinary_cartesian_local_exact_corrected_hamiltonian(initial_h, spec.Z, local_index)
        corrected_h = one_body.h
        corrected_eigenvalue = one_body.corrected_eigenvalue
        corrected_orbital = one_body.corrected_orbital
        one_body_shift = one_body.shift
        local_exact_bracket_iterations = one_body.bracket_iterations
        local_exact_bisection_iterations = one_body.bisection_iterations
    end

    corrected_v = initial_v
    two_body_target = 5.0 * spec.Z / 8.0
    initial_j = _ordinary_cartesian_ida_coulomb_scalar(initial_v, corrected_orbital)
    corrected_j = initial_j
    two_body_shift = 0.0
    two_body_scalar_delta = 0.0
    local_density_weight = begin
        weights = abs2.(corrected_orbital)
        weights ./= sum(weights)
        weights[local_index]
    end
    if spec.two_body_mode == :esoi_local
        two_body =
            _ordinary_cartesian_esoi_corrected_interaction(initial_v, corrected_orbital, spec.Z, local_index)
        corrected_v = two_body.v
        two_body_target = two_body.target
        initial_j = two_body.initial_j
        corrected_j = two_body.corrected_j
        two_body_shift = two_body.shift
        two_body_scalar_delta = two_body.scalar_delta
        local_density_weight = two_body.density_weight
    elseif spec.two_body_mode == :none
        # Keep the interaction unchanged.
    end

    proxy_initial_energy =
        _ordinary_cartesian_closed_shell_proxy_energy(initial_h, initial_v, initial_orbital)
    proxy_corrected_energy =
        _ordinary_cartesian_closed_shell_proxy_energy(corrected_h, corrected_v, corrected_orbital)
    proxy_target_energy = -spec.Z^2 + 5.0 * spec.Z / 8.0

    diagnostics = (
        Z = spec.Z,
        nucleus = spec.nucleus,
        one_body_mode = spec.one_body_mode,
        two_body_mode = spec.two_body_mode,
        local_selection = spec.local_selection,
        local_orbital_index = local_index,
        local_orbital_label = local_orbital.label,
        local_orbital_kind = local_orbital.kind,
        local_orbital_center = (local_orbital.x, local_orbital.y, local_orbital.z),
        local_orbital_distance = local_distance,
        local_density_weight = local_density_weight,
        overlap_error = overlap_error,
        exact_lowest_core_eigenvalue = target_eigenvalue,
        initial_lowest_core_eigenvalue = initial_eigenvalue,
        corrected_lowest_core_eigenvalue = corrected_eigenvalue,
        one_body_shift = one_body_shift,
        local_exact_bracket_iterations = local_exact_bracket_iterations,
        local_exact_bisection_iterations = local_exact_bisection_iterations,
        exact_1s_coulomb = two_body_target,
        initial_1s_coulomb = initial_j,
        corrected_1s_coulomb = corrected_j,
        two_body_scalar_delta = two_body_scalar_delta,
        two_body_local_shift = two_body_shift,
        closed_shell_initial_energy = proxy_initial_energy,
        closed_shell_corrected_energy = proxy_corrected_energy,
        closed_shell_target_energy = proxy_target_energy,
    )
    return (
        one_body_hamiltonian = corrected_h,
        interaction_matrix = corrected_v,
        diagnostics = merge(diagnostics, extra_diagnostics),
    )
end

function _apply_ordinary_cartesian_corrections(
    operators::OrdinaryCartesianOperators3D,
    spec::HydrogenicCoreCorrectionSpec;
    overlap_tol::Real = 1.0e-8,
)
    overlap_error = norm(operators.overlap - I, Inf)
    overlap_error <= Float64(overlap_tol) || throw(
        ArgumentError(
            "_apply_ordinary_cartesian_corrections currently requires an orthonormal final ordinary Cartesian basis; got overlap error $(overlap_error)",
        ),
    )

    matrices = _ordinary_cartesian_corrected_matrices(
        operators,
        operators.one_body_hamiltonian,
        operators.interaction_matrix,
        spec;
        overlap_error = overlap_error,
    )
    corrected_operators = _ordinary_cartesian_operators_with_corrections(
        operators,
        matrices.one_body_hamiltonian,
        matrices.interaction_matrix,
    )
    return OrdinaryCartesianCorrectionResult(corrected_operators, matrices.diagnostics)
end

function _ordinary_cartesian_internal_spec(spec::HydrogenicCoreProjectorCorrectionSpec)
    return HydrogenicCoreCorrectionSpec(;
        Z = spec.Z,
        nucleus = spec.nucleus,
        one_body_mode = :projector,
        two_body_mode = spec.include_esoi ? :esoi_local : :none,
        local_selection = spec.local_orbital_index === nothing ? :nearest_nonresidual : :explicit,
        local_orbital_index = spec.local_orbital_index,
    )
end

function _ordinary_cartesian_internal_spec(spec::HydrogenicCoreBranchCorrectionSpec)
    spec.orbital_selector in (:localized_lowest, :global_lowest) || throw(
        ArgumentError(
            "unsupported orbital_selector :$(spec.orbital_selector); supported branch selectors are :localized_lowest and :global_lowest",
        ),
    )
    return HydrogenicCoreCorrectionSpec(;
        Z = spec.Z,
        nucleus = spec.nucleus,
        one_body_mode = :projector,
        two_body_mode = spec.include_esoi ? :esoi_local : :none,
        local_selection = spec.local_orbital_index === nothing ? :nearest_nonresidual : :explicit,
        local_orbital_index = spec.local_orbital_index,
    )
end

function _ordinary_cartesian_branch_corrections(
    corrections::HydrogenicCoreBranchCorrectionSpec,
)
    return (corrections,)
end

function _ordinary_cartesian_branch_corrections(corrections)
    correction_vector = collect(corrections)
    !isempty(correction_vector) || throw(
        ArgumentError(
            "ordinary_cartesian_corrected_branch requires at least one branch correction spec",
        ),
    )
    for correction in correction_vector
        correction isa HydrogenicCoreBranchCorrectionSpec || throw(
            ArgumentError(
                "ordinary_cartesian_corrected_branch requires HydrogenicCoreBranchCorrectionSpec corrections",
            ),
        )
    end
    if length(correction_vector) > 1
        for correction in correction_vector
            correction.orbital_selector == :localized_lowest || throw(
                ArgumentError(
                    "multi-correction ordinary_cartesian_corrected_branch calls require orbital_selector = :localized_lowest for every correction",
                ),
            )
        end
    end
    return Tuple(correction_vector)
end

function _ordinary_cartesian_branch_partition_centers(
    operators::OrdinaryCartesianOperators3D,
    spec::HydrogenicCoreBranchCorrectionSpec,
)
    if hasproperty(operators.basis, :nuclei)
        centers = NTuple{3,Float64}[
            (Float64(nucleus[1]), Float64(nucleus[2]), Float64(nucleus[3])) for
            nucleus in getproperty(operators.basis, :nuclei)
        ]
        !isempty(centers) || throw(
            ArgumentError("localized branch correction requires at least one carried nucleus"),
        )
        return centers, :basis_nuclei
    end

    if operators.nuclear_charges !== nothing && length(operators.nuclear_charges) > 1
        throw(
            ArgumentError(
                "localized branch correction cannot infer a multi-center partition from this operator payload; use a basis carrying nuclei",
            ),
        )
    end
    return NTuple{3,Float64}[spec.nucleus], :correction_nucleus
end

function _ordinary_cartesian_nearest_center_index(
    point::NTuple{3,Float64},
    centers::AbstractVector{<:NTuple{3,Float64}},
)
    best_index = 0
    best_distance2 = Inf
    px, py, pz = point
    for (index, center) in pairs(centers)
        cx, cy, cz = center
        distance2 = (px - cx)^2 + (py - cy)^2 + (pz - cz)^2
        if distance2 < best_distance2
            best_index = index
            best_distance2 = distance2
        end
    end
    best_index > 0 || throw(ArgumentError("cannot choose nearest center from an empty center list"))
    return best_index, sqrt(best_distance2)
end

function _ordinary_cartesian_branch_localized_orbital_selection(
    operators::OrdinaryCartesianOperators3D,
    h::AbstractMatrix{<:Real},
    spec::HydrogenicCoreBranchCorrectionSpec,
)
    centers, center_source = _ordinary_cartesian_branch_partition_centers(operators, spec)
    selected_center_index, selected_center_distance =
        _ordinary_cartesian_nearest_center_index(spec.nucleus, centers)
    subspace_indices = Int[]
    for (index, orbital) in pairs(operators.orbital_data)
        orbital_center = (orbital.x, orbital.y, orbital.z)
        assigned_center_index, _ = _ordinary_cartesian_nearest_center_index(orbital_center, centers)
        assigned_center_index == selected_center_index && push!(subspace_indices, index)
    end
    !isempty(subspace_indices) || throw(
        ArgumentError(
            "localized branch correction selected an empty center-local subspace for center index $(selected_center_index)",
        ),
    )

    local_h = Matrix{Float64}(h)[subspace_indices, subspace_indices]
    local_eigenvalue, local_orbital = _ordinary_cartesian_lowest_orbital(local_h)
    embedded_orbital = zeros(Float64, size(h, 1))
    embedded_orbital[subspace_indices] .= local_orbital
    embedded_orbital ./= norm(embedded_orbital)
    selected_expectation = _ordinary_cartesian_one_body_expectation(h, embedded_orbital)
    return (
        orbital = embedded_orbital,
        selected_center = centers[selected_center_index],
        selected_center_index = selected_center_index,
        selected_center_source = center_source,
        selected_center_distance = selected_center_distance,
        selected_local_subspace_dimension = length(subspace_indices),
        selected_local_subspace_indices = Tuple(subspace_indices),
        selected_local_eigenvalue = local_eigenvalue,
        selected_core_initial_expectation = selected_expectation,
    )
end

function _ordinary_cartesian_branch_global_orbital_selection(
    operators::OrdinaryCartesianOperators3D,
    h::AbstractMatrix{<:Real},
    spec::HydrogenicCoreBranchCorrectionSpec,
)
    eigenvalue, orbital = _ordinary_cartesian_lowest_orbital(h)
    return (
        orbital = orbital,
        selected_center = spec.nucleus,
        selected_center_index = nothing,
        selected_center_source = :global_full_space,
        selected_center_distance = 0.0,
        selected_local_subspace_dimension = length(operators.orbital_data),
        selected_local_subspace_indices = Tuple(collect(eachindex(operators.orbital_data))),
        selected_local_eigenvalue = eigenvalue,
        selected_core_initial_expectation = eigenvalue,
    )
end

function _ordinary_cartesian_branch_orbital_selection(
    operators::OrdinaryCartesianOperators3D,
    h::AbstractMatrix{<:Real},
    spec::HydrogenicCoreBranchCorrectionSpec,
)
    if spec.orbital_selector == :localized_lowest
        return _ordinary_cartesian_branch_localized_orbital_selection(operators, h, spec)
    elseif spec.orbital_selector == :global_lowest
        return _ordinary_cartesian_branch_global_orbital_selection(operators, h, spec)
    end
    throw(
        ArgumentError(
            "unsupported orbital_selector :$(spec.orbital_selector); supported branch selectors are :localized_lowest and :global_lowest",
        ),
    )
end

function _ordinary_cartesian_localized_selected_center_index(
    operators::OrdinaryCartesianOperators3D,
    spec::HydrogenicCoreBranchCorrectionSpec,
)
    spec.orbital_selector == :localized_lowest || return nothing
    centers, _ = _ordinary_cartesian_branch_partition_centers(operators, spec)
    selected_center_index, _ = _ordinary_cartesian_nearest_center_index(spec.nucleus, centers)
    return selected_center_index
end

function _ordinary_cartesian_validate_localized_branch_charge!(
    charges,
    selected_center_index,
    expected_charge::Real,
    label::AbstractString,
)
    charges === nothing && return nothing
    selected_center_index === nothing && return nothing
    selected_center_index <= length(charges) || throw(
        ArgumentError(
            "$(label) must carry a charge for localized correction center index $(selected_center_index); got $(length(charges)) charges",
        ),
    )
    charge = Float64(charges[selected_center_index])
    isapprox(charge, Float64(expected_charge); atol = 1.0e-10, rtol = 1.0e-10) || throw(
        ArgumentError(
            "$(label) has charge $(charge) at localized correction center index $(selected_center_index), expected $(Float64(expected_charge))",
        ),
    )
    return nothing
end

function _ordinary_cartesian_inferred_reference_nuclear_charges(
    operators::OrdinaryCartesianOperators3D,
    spec::HydrogenicCoreBranchCorrectionSpec;
    selected_center_index = nothing,
)
    spec.reference_nuclear_charges !== nothing && return copy(spec.reference_nuclear_charges)
    spec.orbital_selector == :localized_lowest || return nothing
    hasproperty(operators.basis, :nuclei) || return nothing
    centers, _ = _ordinary_cartesian_branch_partition_centers(operators, spec)
    if selected_center_index === nothing
        selected_center_index, _ = _ordinary_cartesian_nearest_center_index(spec.nucleus, centers)
    end
    charges = zeros(Float64, length(centers))
    charges[selected_center_index] = spec.Z
    return charges
end

function _ordinary_cartesian_reference_one_body_hamiltonian(
    operators::OrdinaryCartesianOperators3D,
    spec::HydrogenicCoreBranchCorrectionSpec,
    application_h::AbstractMatrix{<:Real},
    selected_center_index,
)
    reference_charges = _ordinary_cartesian_inferred_reference_nuclear_charges(
        operators,
        spec;
        selected_center_index = selected_center_index,
    )
    if reference_charges === nothing
        return Matrix{Float64}(application_h), nothing
    end
    if spec.reference_nuclear_charges !== nothing
        _ordinary_cartesian_validate_localized_branch_charge!(
            reference_charges,
            selected_center_index,
            spec.Z,
            "reference_nuclear_charges",
        )
    end
    reference_h = assembled_one_body_hamiltonian(operators; nuclear_charges = reference_charges)
    return reference_h, Tuple(reference_charges)
end

function _ordinary_cartesian_branch_corrected_matrices(
    operators::OrdinaryCartesianOperators3D,
    application_h::AbstractMatrix{<:Real},
    reference_h::AbstractMatrix{<:Real},
    initial_v::AbstractMatrix{<:Real},
    correction::HydrogenicCoreBranchCorrectionSpec,
    internal_spec::HydrogenicCoreCorrectionSpec;
    reference_nuclear_charges,
)
    application_h = Matrix{Float64}(application_h)
    reference_h = Matrix{Float64}(reference_h)
    initial_v = Matrix{Float64}(initial_v)
    application_global_eigenvalue, _ = _ordinary_cartesian_lowest_orbital(application_h)
    calibration_global_eigenvalue, _ = _ordinary_cartesian_lowest_orbital(reference_h)
    selection = _ordinary_cartesian_branch_orbital_selection(operators, reference_h, correction)
    local_index, local_orbital, local_distance =
        _ordinary_cartesian_select_local_orbital(operators, internal_spec)

    one_body = _ordinary_cartesian_projector_delta_from_reference(
        application_h,
        reference_h,
        correction.Z,
        selection.orbital,
    )
    corrected_h = one_body.h
    corrected_orbital = one_body.selected_orbital

    corrected_v = initial_v
    two_body_target = 5.0 * correction.Z / 8.0
    initial_j = _ordinary_cartesian_ida_coulomb_scalar(initial_v, corrected_orbital)
    corrected_j = initial_j
    two_body_shift = 0.0
    two_body_scalar_delta = 0.0
    local_density_weight = begin
        weights = abs2.(corrected_orbital)
        weights ./= sum(weights)
        weights[local_index]
    end
    if internal_spec.two_body_mode == :esoi_local
        two_body = _ordinary_cartesian_esoi_corrected_interaction(
            initial_v,
            corrected_orbital,
            correction.Z,
            local_index,
        )
        corrected_v = two_body.v
        two_body_target = two_body.target
        initial_j = two_body.initial_j
        corrected_j = two_body.corrected_j
        two_body_shift = two_body.shift
        two_body_scalar_delta = two_body.scalar_delta
        local_density_weight = two_body.density_weight
    end

    proxy_initial_energy =
        _ordinary_cartesian_closed_shell_proxy_energy(application_h, initial_v, corrected_orbital)
    proxy_corrected_energy =
        _ordinary_cartesian_closed_shell_proxy_energy(corrected_h, corrected_v, corrected_orbital)
    proxy_target_energy = -correction.Z^2 + 5.0 * correction.Z / 8.0

    diagnostics = (
        Z = correction.Z,
        nucleus = correction.nucleus,
        one_body_mode = :projector,
        two_body_mode = internal_spec.two_body_mode,
        local_selection = internal_spec.local_selection,
        local_orbital_index = local_index,
        local_orbital_label = local_orbital.label,
        local_orbital_kind = local_orbital.kind,
        local_orbital_center = (local_orbital.x, local_orbital.y, local_orbital.z),
        local_orbital_distance = local_distance,
        local_density_weight = local_density_weight,
        reference_nuclear_charges = reference_nuclear_charges,
        exact_lowest_core_eigenvalue = one_body.target,
        initial_lowest_core_eigenvalue = calibration_global_eigenvalue,
        corrected_lowest_core_eigenvalue = one_body.corrected_eigenvalue,
        calibration_lowest_core_eigenvalue = calibration_global_eigenvalue,
        application_lowest_core_eigenvalue = application_global_eigenvalue,
        one_body_shift = one_body.shift,
        local_exact_bracket_iterations = 0,
        local_exact_bisection_iterations = 0,
        exact_1s_coulomb = two_body_target,
        initial_1s_coulomb = initial_j,
        corrected_1s_coulomb = corrected_j,
        two_body_scalar_delta = two_body_scalar_delta,
        two_body_local_shift = two_body_shift,
        closed_shell_initial_energy = proxy_initial_energy,
        closed_shell_corrected_energy = proxy_corrected_energy,
        closed_shell_target_energy = proxy_target_energy,
        orbital_selector = correction.orbital_selector,
        selected_center = selection.selected_center,
        selected_center_index = selection.selected_center_index,
        selected_center_source = selection.selected_center_source,
        selected_center_distance = selection.selected_center_distance,
        selected_local_subspace_dimension = selection.selected_local_subspace_dimension,
        selected_local_subspace_count = selection.selected_local_subspace_dimension,
        selected_local_subspace_indices = selection.selected_local_subspace_indices,
        selected_core_initial_expectation = one_body.calibration_expectation,
        selected_core_corrected_expectation = one_body.application_corrected_expectation,
        calibration_core_initial_expectation = one_body.calibration_expectation,
        calibration_core_corrected_expectation = one_body.target,
        application_core_initial_expectation = one_body.application_initial_expectation,
        application_core_corrected_expectation = one_body.application_corrected_expectation,
        selected_local_eigenvalue = selection.selected_local_eigenvalue,
    )
    return (
        one_body_hamiltonian = corrected_h,
        interaction_matrix = corrected_v,
        diagnostics = diagnostics,
    )
end

function _ordinary_cartesian_check_duplicate_branch_center!(
    corrected_center_indices::Vector,
    selected_center_index,
)
    selected_center_index === nothing && return nothing
    selected_center_index in corrected_center_indices && throw(
        ArgumentError(
            "ordinary_cartesian_corrected_branch received duplicate localized corrections for center index $(selected_center_index)",
        ),
    )
    push!(corrected_center_indices, selected_center_index)
    return nothing
end

"""
    apply_ordinary_cartesian_corrections(operators, spec; overlap_tol=1e-8)
    apply_ordinary_cartesian_corrections(operators; Z, nucleus=(0,0,0),
        include_esoi=false, local_orbital_index=nothing, overlap_tol=1e-8)

Apply the public projector-based hydrogenic core correction to an assembled
[`OrdinaryCartesianOperators3D`](@ref). By default this fixes only the lowest
one-body hydrogenic core eigenvalue to `-Z^2/2`. Pass `include_esoi=true` to
also calibrate the corresponding closed-shell `1s^2` Coulomb scalar to `5Z/8`.

This is a post-assembly transform. The returned operators preserve the total
matrices, but corrected payloads no longer have a clean kinetic/nuclear
decomposition; the invalid sidecars are dropped and `nuclear_term_storage` is
set to `:total_only`.
"""
function apply_ordinary_cartesian_corrections(
    operators::OrdinaryCartesianOperators3D,
    spec::HydrogenicCoreProjectorCorrectionSpec;
    overlap_tol::Real = 1.0e-8,
)
    return _apply_ordinary_cartesian_corrections(
        operators,
        _ordinary_cartesian_internal_spec(spec);
        overlap_tol = overlap_tol,
    )
end

function apply_ordinary_cartesian_corrections(
    operators::OrdinaryCartesianOperators3D;
    Z::Real,
    nucleus = (0.0, 0.0, 0.0),
    include_esoi::Bool = false,
    local_orbital_index::Union{Nothing,Integer} = nothing,
    overlap_tol::Real = 1.0e-8,
)
    spec = HydrogenicCoreProjectorCorrectionSpec(;
        Z = Z,
        nucleus = nucleus,
        include_esoi = include_esoi,
        local_orbital_index = local_orbital_index,
    )
    return apply_ordinary_cartesian_corrections(operators, spec; overlap_tol = overlap_tol)
end

"""
    ordinary_cartesian_corrected_branch(operators; nuclear_charges, corrections,
        overlap_tol=1e-8)

Apply a narrow branch-level hydrogenic correction to an
[`OrdinaryCartesianOperators3D`](@ref) payload and return corrected matrices.
The one-body branch matrix is assembled with `assembled_one_body_hamiltonian`
using `nuclear_charges`, then one or more
[`HydrogenicCoreBranchCorrectionSpec`](@ref)s are applied sequentially.

This is the intended public entry point for branch/counterpoise correction
work. The default `orbital_selector = :localized_lowest` uses the selected
center-local subspace. `orbital_selector = :global_lowest` remains available
only for single-correction debug/reference calls. Multi-correction calls reject
duplicate localized center indices. The result is a matrix-level
[`OrdinaryCartesianBranchCorrectionResult`](@ref), not a transformed
`OrdinaryCartesianOperators3D`.
"""
function ordinary_cartesian_corrected_branch(
    operators::OrdinaryCartesianOperators3D;
    nuclear_charges,
    corrections,
    overlap_tol::Real = 1.0e-8,
)
    overlap_error = norm(operators.overlap - I, Inf)
    overlap_error <= Float64(overlap_tol) || throw(
        ArgumentError(
            "ordinary_cartesian_corrected_branch currently requires an orthonormal final ordinary Cartesian basis; got overlap error $(overlap_error)",
        ),
    )

    correction_specs = _ordinary_cartesian_branch_corrections(corrections)
    branch_charges = nuclear_charges === nothing ? nothing : Float64.(collect(nuclear_charges))
    initial_application_h = assembled_one_body_hamiltonian(operators; nuclear_charges = branch_charges)
    initial_application_v = Matrix{Float64}(operators.interaction_matrix)
    corrected_h = Matrix{Float64}(initial_application_h)
    corrected_v = Matrix{Float64}(initial_application_v)
    per_correction_diagnostics = NamedTuple[]
    corrected_center_indices = Any[]
    for correction in correction_specs
        internal_spec = _ordinary_cartesian_internal_spec(correction)
        selected_center_index = _ordinary_cartesian_localized_selected_center_index(operators, correction)
        _ordinary_cartesian_validate_localized_branch_charge!(
            branch_charges,
            selected_center_index,
            correction.Z,
            "application nuclear_charges",
        )
        reference_h, reference_charges = _ordinary_cartesian_reference_one_body_hamiltonian(
            operators,
            correction,
            corrected_h,
            selected_center_index,
        )
        matrices = _ordinary_cartesian_branch_corrected_matrices(
            operators,
            corrected_h,
            reference_h,
            corrected_v,
            correction,
            internal_spec;
            reference_nuclear_charges = reference_charges,
        )
        diagnostic = matrices.diagnostics
        if length(correction_specs) > 1
            _ordinary_cartesian_check_duplicate_branch_center!(
                corrected_center_indices,
                diagnostic.selected_center_index,
            )
        elseif diagnostic.selected_center_index !== nothing
            push!(corrected_center_indices, diagnostic.selected_center_index)
        end
        corrected_h = matrices.one_body_hamiltonian
        corrected_v = matrices.interaction_matrix
        push!(per_correction_diagnostics, diagnostic)
    end

    diagnostics = (
        branch_nuclear_charges = branch_charges === nothing ? nothing : Tuple(branch_charges),
        correction_count = length(correction_specs),
        corrections = Tuple(per_correction_diagnostics),
        overlap_error = overlap_error,
        corrected_center_indices = Tuple(corrected_center_indices),
    )
    one_body_delta = corrected_h .- initial_application_h
    interaction_delta = corrected_v .- initial_application_v
    return OrdinaryCartesianBranchCorrectionResult(
        corrected_h,
        corrected_v,
        one_body_delta,
        interaction_delta,
        diagnostics,
    )
end
