const _QWRG_RESIDUAL_KEEP_ABS_TOL = 1.0e-8
const _QWRG_RESIDUAL_KEEP_REL_TOL = 1.0e-1
const _QWRG_RESIDUAL_ACCEPT_TOL = 1.0e-8
const _QWRG_ATOMIC_RESIDUAL_KEEP_ABS_TOL = 1.0e-7
const _QWRG_ATOMIC_RESIDUAL_ACCEPT_TOL = 1.0e-7
const _QWRG_RESIDUAL_NULL_ABS_TOL = 1.0e-12
const _QWRG_RESIDUAL_NULL_REL_TOL = 1.0e-12
const _QWRG_RESIDUAL_STABILIZATION_TARGET_TOL = 1.0e-10
const _QWRG_RESIDUAL_STABILIZATION_MAX_PASSES = 3

function _qwrg_residual_null_rank_tol(values::AbstractVector{<:Real})
    isempty(values) && return _QWRG_RESIDUAL_NULL_ABS_TOL
    max_value = maximum(abs.(Float64.(values)))
    return max(_QWRG_RESIDUAL_NULL_ABS_TOL, _QWRG_RESIDUAL_NULL_REL_TOL * max_value)
end

function _qwrg_residual_keep_policy(
    keep_policy::Symbol;
    allow_relative_case_scale::Bool = true,
)
    keep_policy == :legacy_profile && return :near_null_only
    keep_policy == :near_null_only && return keep_policy
    keep_policy == :relative_case_scale && allow_relative_case_scale && return keep_policy
    allowed = allow_relative_case_scale ?
        ":near_null_only, :relative_case_scale, or the compatibility alias :legacy_profile" :
        ":near_null_only or the compatibility alias :legacy_profile"
    throw(
        ArgumentError(
            "QW residual keep policy must be $(allowed)",
        ),
    )
end

function _qwrg_residual_keep_tol(
    values::AbstractVector{<:Real};
    keep_policy::Symbol = :near_null_only,
    abs_tol::Real = _QWRG_RESIDUAL_KEEP_ABS_TOL,
)
    keep_policy_value = _qwrg_residual_keep_policy(keep_policy)
    abs_tol_value = Float64(abs_tol)
    if keep_policy_value == :relative_case_scale
        isempty(values) && return abs_tol_value
        return max(abs_tol_value, _QWRG_RESIDUAL_KEEP_REL_TOL * maximum(values))
    end
    return abs_tol_value
end

function _qwrg_atomic_residual_keep_policy(keep_policy::Symbol)
    return _qwrg_residual_keep_policy(
        keep_policy;
        allow_relative_case_scale = false,
    )
end

_qwrg_atomic_residual_keep_tol() = _QWRG_ATOMIC_RESIDUAL_KEEP_ABS_TOL
_qwrg_atomic_residual_accept_tol() = _QWRG_ATOMIC_RESIDUAL_ACCEPT_TOL

function _qwrg_stabilize_residual_coefficients(
    raw_overlap::AbstractMatrix{<:Real},
    residual_coefficients::AbstractMatrix{<:Real},
)
    function summarize_overlap_block(
        overlap_block::AbstractMatrix{<:Real},
        null_tol::Float64,
    )
        overlap_value = Matrix{Float64}(overlap_block)
        overlap_sym = Matrix{Float64}(0.5 .* (overlap_value .+ transpose(overlap_value)))
        decomposition = eigen(Symmetric(overlap_sym))
        values = Float64[decomposition.values...]
        return (
            matrix = overlap_value,
            symmetrized = overlap_sym,
            vectors = Matrix{Float64}(decomposition.vectors),
            values = values,
            overlap_error = norm(
                overlap_value - Matrix{Float64}(I, size(overlap_value, 1), size(overlap_value, 2)),
                Inf,
            ),
            symmetry_defect = norm(overlap_value - transpose(overlap_value), Inf),
            min_eigenvalue = isempty(values) ? 1.0 : minimum(values),
            max_eigenvalue = isempty(values) ? 1.0 : maximum(values),
            negative_count = count(<(-null_tol), values),
            near_null_count = count(<=(null_tol), values),
        )
    end

    nresidual = size(residual_coefficients, 2)
    nresidual == 0 && return (
        coefficients = Matrix{Float64}(residual_coefficients),
        pre_overlap = Matrix{Float64}(I, 0, 0),
        post_overlap = Matrix{Float64}(I, 0, 0),
        pre_error = 0.0,
        post_error = 0.0,
        pre_symmetry_defect = 0.0,
        post_symmetry_defect = 0.0,
        pre_min_eigenvalue = 1.0,
        pre_max_eigenvalue = 1.0,
        post_min_eigenvalue = 1.0,
        post_max_eigenvalue = 1.0,
        pre_negative_count = 0,
        post_negative_count = 0,
        pre_near_null_count = 0,
        post_near_null_count = 0,
        null_tol = _QWRG_RESIDUAL_NULL_ABS_TOL,
        clipped_count = 0,
        dropped_count = 0,
        correction_passes = 0,
    )
    raw_overlap_value = Matrix{Float64}(raw_overlap)
    stabilized_coefficients = Matrix{Float64}(residual_coefficients)

    pre_overlap = Matrix{Float64}(transpose(stabilized_coefficients) * raw_overlap_value * stabilized_coefficients)
    pre_overlap_sym = Matrix{Float64}(0.5 .* (pre_overlap .+ transpose(pre_overlap)))
    null_tol = _qwrg_residual_null_rank_tol(eigvals(Symmetric(pre_overlap_sym)))
    pre_summary = summarize_overlap_block(pre_overlap, null_tol)
    clipped_count = count(<(null_tol), pre_summary.values)
    pre_summary.negative_count == 0 || throw(
        ArgumentError(
            "QW residual-space stabilization requires the kept residual overlap block to stay positive semidefinite",
        ),
    )

    stabilized_values = max.(pre_summary.values, null_tol)
    inverse_sqrt = Matrix{Float64}(
        pre_summary.vectors *
        Diagonal(1.0 ./ sqrt.(stabilized_values)) *
        transpose(pre_summary.vectors),
    )
    stabilized_coefficients = Matrix{Float64}(stabilized_coefficients * inverse_sqrt)
    correction_passes = 1
    current_overlap = Matrix{Float64}(transpose(stabilized_coefficients) * raw_overlap_value * stabilized_coefficients)
    current_summary = summarize_overlap_block(current_overlap, null_tol)
    current_summary.negative_count == 0 || throw(
        ArgumentError(
            "QW residual-space stabilization produced a numerically indefinite kept residual overlap block",
        ),
    )
    while correction_passes < _QWRG_RESIDUAL_STABILIZATION_MAX_PASSES &&
        (current_summary.overlap_error > _QWRG_RESIDUAL_STABILIZATION_TARGET_TOL ||
         current_summary.symmetry_defect > _QWRG_RESIDUAL_STABILIZATION_TARGET_TOL)
        cholesky_factor = cholesky(Symmetric(current_summary.symmetrized))
        stabilized_coefficients = Matrix{Float64}(stabilized_coefficients / cholesky_factor.U)
        correction_passes += 1
        current_overlap = Matrix{Float64}(transpose(stabilized_coefficients) * raw_overlap_value * stabilized_coefficients)
        current_summary = summarize_overlap_block(current_overlap, null_tol)
        current_summary.negative_count == 0 || throw(
            ArgumentError(
                "QW residual-space stabilization produced a numerically indefinite kept residual overlap block",
            ),
        )
    end

    return (
        coefficients = stabilized_coefficients,
        pre_overlap = pre_summary.symmetrized,
        post_overlap = current_summary.symmetrized,
        pre_error = pre_summary.overlap_error,
        post_error = current_summary.overlap_error,
        pre_symmetry_defect = pre_summary.symmetry_defect,
        post_symmetry_defect = current_summary.symmetry_defect,
        pre_min_eigenvalue = pre_summary.min_eigenvalue,
        pre_max_eigenvalue = pre_summary.max_eigenvalue,
        post_min_eigenvalue = current_summary.min_eigenvalue,
        post_max_eigenvalue = current_summary.max_eigenvalue,
        pre_negative_count = pre_summary.negative_count,
        post_negative_count = current_summary.negative_count,
        pre_near_null_count = pre_summary.near_null_count,
        post_near_null_count = current_summary.near_null_count,
        null_tol = null_tol,
        clipped_count = clipped_count,
        dropped_count = 0,
        correction_passes = correction_passes,
    )
end

# Alg QW-RG step 4: Define residual Gaussians by orthogonalizing 3D GTOs
# to the full 3D fixed working space. For the active PGDG-mediated route,
# the carried 1D PGDG auxiliary line includes the COMX cleanup step inside the
# bundle, so this fixed block should already be orthonormal to numerical
# precision before the residual construction.
# See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
function _qwrg_residual_space_analysis(
    gausslet_overlap::AbstractMatrix{<:Real},
    overlap_ga::AbstractMatrix{<:Real},
    overlap_aa::AbstractMatrix{<:Real},
    ;
    keep_policy::Symbol = :near_null_only,
    keep_abs_tol::Real = _QWRG_RESIDUAL_KEEP_ABS_TOL,
    accept_tol::Real = _QWRG_RESIDUAL_ACCEPT_TOL,
)
    keep_policy_value = _qwrg_residual_keep_policy(keep_policy)
    accept_tol_value = Float64(accept_tol)
    gausslet_overlap_value = Matrix{Float64}(gausslet_overlap)
    overlap_error = norm(
        gausslet_overlap_value - Matrix{Float64}(I, size(gausslet_overlap_value, 1), size(gausslet_overlap_value, 2)),
        Inf,
    )
    overlap_error <= 1.0e-8 || throw(
        ArgumentError(
            "Qiu-White / QW-PGDG residual construction requires an orthonormal fixed 3D overlap block",
        ),
    )

    ngausslet = size(gausslet_overlap_value, 1)
    ngaussian = size(overlap_aa, 1)
    raw_overlap = [
        gausslet_overlap_value Matrix{Float64}(overlap_ga)
        Matrix{Float64}(transpose(overlap_ga)) Matrix{Float64}(overlap_aa)
    ]
    seed_projector = vcat(
        -(gausslet_overlap_value \ Matrix{Float64}(overlap_ga)),
        Matrix{Float64}(I, ngaussian, ngaussian),
    )
    residual_overlap = Matrix{Float64}(transpose(seed_projector) * raw_overlap * seed_projector)
    supplement_decomposition = eigen(Symmetric(Matrix{Float64}(overlap_aa)))
    residual_decomposition = eigen(Symmetric(residual_overlap))

    supplement_null_rank_tol = _qwrg_residual_null_rank_tol(supplement_decomposition.values)
    supplement_numerical_rank = count(>(supplement_null_rank_tol), supplement_decomposition.values)
    residual_null_rank_tol = _qwrg_residual_null_rank_tol(residual_decomposition.values)
    residual_numerical_rank = count(>(residual_null_rank_tol), residual_decomposition.values)
    keep_tol = _qwrg_residual_keep_tol(
        residual_decomposition.values;
        keep_policy = keep_policy_value,
        abs_tol = keep_abs_tol,
    )
    keep = findall(>(keep_tol), residual_decomposition.values)
    discarded = setdiff(collect(1:length(residual_decomposition.values)), keep)
    kept_residual_coefficients = isempty(keep) ? zeros(Float64, size(seed_projector, 1), 0) :
        Matrix{Float64}(
            seed_projector *
            residual_decomposition.vectors[:, keep] *
            Diagonal(1.0 ./ sqrt.(residual_decomposition.values[keep])),
        )
    stabilization = _qwrg_stabilize_residual_coefficients(raw_overlap, kept_residual_coefficients)

    diagnostics = QWRGResidualSpaceDiagnostics(
        ngausslet,
        ngaussian,
        ngausslet + ngaussian,
        size(overlap_ga),
        size(overlap_aa),
        overlap_error,
        supplement_null_rank_tol,
        supplement_numerical_rank,
        residual_null_rank_tol,
        residual_numerical_rank,
        keep_policy_value,
        keep_tol,
        accept_tol_value,
        length(keep),
        ngaussian - length(keep),
        stabilization.null_tol,
        stabilization.correction_passes,
        stabilization.clipped_count,
        stabilization.dropped_count,
        stabilization.pre_error,
        stabilization.post_error,
        stabilization.pre_symmetry_defect,
        stabilization.post_symmetry_defect,
        stabilization.pre_min_eigenvalue,
        stabilization.pre_max_eigenvalue,
        stabilization.post_min_eigenvalue,
        stabilization.post_max_eigenvalue,
        stabilization.pre_negative_count,
        stabilization.post_negative_count,
        stabilization.pre_near_null_count,
        stabilization.post_near_null_count,
        Float64[supplement_decomposition.values...],
        Float64[residual_decomposition.values...],
        keep,
        discarded,
        Float64[residual_decomposition.values[index] for index in keep],
        Float64[residual_decomposition.values[index] for index in discarded],
    )
    return (
        raw_overlap = raw_overlap,
        seed_projector = seed_projector,
        residual_overlap = residual_overlap,
        decomposition = residual_decomposition,
        keep = keep,
        residual_coefficients_pre_stabilization = kept_residual_coefficients,
        residual_coefficients = stabilization.coefficients,
        diagnostics = diagnostics,
    )
end

function diagnose_qwrg_residual_space(
    gausslet_overlap::AbstractMatrix{<:Real},
    overlap_ga::AbstractMatrix{<:Real},
    overlap_aa::AbstractMatrix{<:Real},
    ;
    keep_policy::Symbol = :near_null_only,
    keep_abs_tol::Real = _QWRG_RESIDUAL_KEEP_ABS_TOL,
    accept_tol::Real = _QWRG_RESIDUAL_ACCEPT_TOL,
)
    return _qwrg_residual_space_analysis(
        gausslet_overlap,
        overlap_ga,
        overlap_aa;
        keep_policy = keep_policy,
        keep_abs_tol = keep_abs_tol,
        accept_tol = accept_tol,
    ).diagnostics
end

function _qwrg_residual_space(
    gausslet_overlap::AbstractMatrix{<:Real},
    overlap_ga::AbstractMatrix{<:Real},
    overlap_aa::AbstractMatrix{<:Real},
    ;
    keep_policy::Symbol = :near_null_only,
    keep_abs_tol::Real = _QWRG_RESIDUAL_KEEP_ABS_TOL,
    accept_tol::Real = _QWRG_RESIDUAL_ACCEPT_TOL,
)
    analysis = _qwrg_residual_space_analysis(
        gausslet_overlap,
        overlap_ga,
        overlap_aa;
        keep_policy = keep_policy,
        keep_abs_tol = keep_abs_tol,
        accept_tol = accept_tol,
    )
    ngausslet = analysis.diagnostics.gausslet_count
    keep = analysis.keep
    isempty(keep) && throw(
        ArgumentError("Qiu-White residual-Gaussian construction produced no nontrivial 3D residual directions"),
    )
    residual_coefficients = analysis.residual_coefficients
    gausslet_coefficients = vcat(
        Matrix{Float64}(I, ngausslet, ngausslet),
        zeros(Float64, analysis.diagnostics.gaussian_count, ngausslet),
    )
    raw_to_final = hcat(gausslet_coefficients, residual_coefficients)
    final_overlap = Matrix{Float64}(transpose(raw_to_final) * analysis.raw_overlap * raw_to_final)
    final_overlap = Matrix{Float64}(0.5 .* (final_overlap .+ transpose(final_overlap)))
    return (
        raw_overlap = analysis.raw_overlap,
        raw_to_final = raw_to_final,
        residual_coefficients = residual_coefficients,
        residual_coefficients_pre_stabilization = analysis.residual_coefficients_pre_stabilization,
        final_overlap = final_overlap,
        diagnostics = analysis.diagnostics,
    )
end

function _qwrg_residual_center_data(
    raw_overlap::AbstractMatrix{<:Real},
    x_raw::AbstractMatrix{<:Real},
    y_raw::AbstractMatrix{<:Real},
    z_raw::AbstractMatrix{<:Real},
    raw_to_final::AbstractMatrix{<:Real},
    ngausslet::Int,
)
    residual_coefficients = Matrix{Float64}(raw_to_final[:, (ngausslet + 1):end])
    nresidual = size(residual_coefficients, 2)
    centers = zeros(Float64, nresidual, 3)
    norms = zeros(Float64, nresidual)

    overlap_residual = Matrix{Float64}(transpose(residual_coefficients) * raw_overlap * residual_coefficients)
    for index in 1:nresidual
        vector = view(residual_coefficients, :, index)
        norm_value = Float64(dot(vector, raw_overlap * vector))
        norm_value > 1.0e-12 || throw(
            ArgumentError("residual center extraction requires nonzero residual norm"),
        )

        norms[index] = norm_value
        centers[index, 1] = Float64(dot(vector, x_raw * vector) / norm_value)
        centers[index, 2] = Float64(dot(vector, y_raw * vector) / norm_value)
        centers[index, 3] = Float64(dot(vector, z_raw * vector) / norm_value)
    end

    overlap_error = norm(overlap_residual - I, Inf)
    return (
        centers = centers,
        overlap_error = overlap_error,
        residual_coefficients = residual_coefficients,
        norms = norms,
    )
end

function _qwrg_residual_width_data(
    raw_overlap::AbstractMatrix{<:Real},
    x2_raw::AbstractMatrix{<:Real},
    y2_raw::AbstractMatrix{<:Real},
    z2_raw::AbstractMatrix{<:Real},
    center_data,
)
    residual_coefficients = center_data.residual_coefficients
    centers = center_data.centers
    norms = center_data.norms
    nresidual = size(residual_coefficients, 2)
    widths = zeros(Float64, nresidual, 3)

    for index in 1:nresidual
        vector = view(residual_coefficients, :, index)
        norm_value = norms[index]
        x1 = centers[index, 1]
        y1 = centers[index, 2]
        z1 = centers[index, 3]
        x2 = Float64(dot(vector, x2_raw * vector) / norm_value)
        y2 = Float64(dot(vector, y2_raw * vector) / norm_value)
        z2 = Float64(dot(vector, z2_raw * vector) / norm_value)

        varx = x2 - x1^2
        vary = y2 - y1^2
        varz = z2 - z1^2
        min(varx, vary, varz) > 1.0e-12 || throw(
            ArgumentError("MWG residual moment extraction requires positive residual variances"),
        )

        widths[index, 1] = sqrt(2.0 * varx)
        widths[index, 2] = sqrt(2.0 * vary)
        widths[index, 3] = sqrt(2.0 * varz)
    end

    overlap_residual = Matrix{Float64}(transpose(residual_coefficients) * raw_overlap * residual_coefficients)
    overlap_error = norm(overlap_residual - I, Inf)
    return (
        centers = centers,
        widths = widths,
        overlap_error = overlap_error,
    )
end

function _qwrg_residual_moment_data(
    raw_overlap::AbstractMatrix{<:Real},
    x_raw::AbstractMatrix{<:Real},
    x2_raw::AbstractMatrix{<:Real},
    y2_raw::AbstractMatrix{<:Real},
    z2_raw::AbstractMatrix{<:Real},
    center_data,
)
    width_data = _qwrg_residual_width_data(
        raw_overlap,
        x2_raw,
        y2_raw,
        z2_raw,
        center_data,
    )
    return (
        centers = center_data.centers,
        widths = width_data.widths,
        overlap_error = width_data.overlap_error,
    )
end

function _qwrg_residual_moment_data(
    raw_overlap::AbstractMatrix{<:Real},
    x_raw::AbstractMatrix{<:Real},
    x2_raw::AbstractMatrix{<:Real},
    y_raw::AbstractMatrix{<:Real},
    y2_raw::AbstractMatrix{<:Real},
    z_raw::AbstractMatrix{<:Real},
    z2_raw::AbstractMatrix{<:Real},
    raw_to_final::AbstractMatrix{<:Real},
    ngausslet::Int,
)
    center_data = _qwrg_residual_center_data(
        raw_overlap,
        x_raw,
        y_raw,
        z_raw,
        raw_to_final,
        ngausslet,
    )
    return _qwrg_residual_moment_data(
        raw_overlap,
        x_raw,
        x2_raw,
        y_raw,
        y2_raw,
        z_raw,
        z2_raw,
        center_data,
    )
end
