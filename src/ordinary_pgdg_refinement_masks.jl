"""
    _TernaryGaussianRefinementMask1D

Internal centered `1 -> 1/3` Gaussian quasi-refinement mask for the first
practical PGDG refinement hierarchy step.

The stored coefficients are ordered by integer fine-grid offsets
`-support_radius:support_radius` and are built directly from the analytic
Gaussian quasi-refinement formula described in:

- `docs/gaussian_refinement_analytic_mask_note.md`
- `docs/gaussian_refinement_analytic_mask_comparison.md`

This internal object is intentionally narrow. It provides the first practical
stored mask and local application helpers, but not the full hierarchy
machinery.
"""
struct _TernaryGaussianRefinementMask1D
    factor::Int
    rho::Float64
    support_radius::Int
    half_window::Float64
    coefficients::Vector{Float64}
end

Base.length(mask::_TernaryGaussianRefinementMask1D) = length(mask.coefficients)

function Base.show(io::IO, mask::_TernaryGaussianRefinementMask1D)
    print(
        io,
        "_TernaryGaussianRefinementMask1D(factor=",
        mask.factor,
        ", rho=",
        mask.rho,
        ", support_radius=",
        mask.support_radius,
        ", half_window=",
        mask.half_window,
        ")",
    )
end

_refinement_mask_offsets(mask::_TernaryGaussianRefinementMask1D) = collect(-mask.support_radius:mask.support_radius)

function _default_ternary_refinement_half_window(rho::Real)
    rho_value = Float64(rho)
    if isapprox(rho_value, 1.2; atol = 1.0e-12, rtol = 0.0)
        return 8.0
    elseif isapprox(rho_value, 4.0 / 3.0; atol = 1.0e-12, rtol = 0.0)
        return 10.0
    elseif isapprox(rho_value, sqrt(2.0); atol = 1.0e-12, rtol = 0.0)
        return 10.0
    elseif isapprox(rho_value, 1.1; atol = 1.0e-12, rtol = 0.0)
        return 7.0
    elseif isapprox(rho_value, 1.0; atol = 1.0e-12, rtol = 0.0)
        return 6.0
    else
        # The first integrated default is rho = 1.2. Later alternatives can
        # either use their studied windows or override this explicitly.
        return 8.0
    end
end

function _analytic_ternary_refinement_coefficient(offset::Integer, rho::Real)
    rho_value = Float64(rho)
    return (3.0 / (4.0 * sqrt(pi) * rho_value)) * exp(-(Float64(offset)^2) / (16.0 * rho_value^2))
end

"""
    _analytic_ternary_gaussian_refinement_mask(; rho = 1.2, half_window = nothing)

Build the first practical centered analytic ternary refinement mask for the 1D
PGDG hierarchy.

Conventions:

- coarse spacing `H = 1`
- fine spacing `h = 1/3`
- coarse width `Sigma = rho`
- fine width `sigma = rho / 3`
- `tau^2 = 8 rho^2 / 9`
- coefficients `c_k = phi_tau(k/3)`

The default `rho = 1.2` and default half-window `±8` match the successful
analytic-vs-fitted comparison study.
"""
function _analytic_ternary_gaussian_refinement_mask(; rho::Real = 1.2, half_window = nothing)
    factor = 3
    half_window_value = half_window === nothing ? _default_ternary_refinement_half_window(rho) : Float64(half_window)
    support_radius = round(Int, factor * half_window_value)
    support_radius >= 0 || throw(ArgumentError("analytic ternary refinement mask requires a nonnegative support radius"))
    offsets = -support_radius:support_radius
    coefficients = Float64[_analytic_ternary_refinement_coefficient(offset, rho) for offset in offsets]
    return _TernaryGaussianRefinementMask1D(
        factor,
        Float64(rho),
        support_radius,
        support_radius / factor,
        coefficients,
    )
end

_default_ternary_gaussian_refinement_mask() = _analytic_ternary_gaussian_refinement_mask()

function _refinement_mask_residue_sums(mask::_TernaryGaussianRefinementMask1D)
    sums = zeros(Float64, mask.factor)
    for (offset, coefficient) in zip(_refinement_mask_offsets(mask), mask.coefficients)
        residue = mod(offset, mask.factor) + 1
        sums[residue] += coefficient
    end
    return sums
end

function _refinement_mask_residue_ripple(mask::_TernaryGaussianRefinementMask1D)
    sums = _refinement_mask_residue_sums(mask)
    mean_value = sum(sums) / length(sums)
    return maximum(abs(value / mean_value - 1.0) for value in sums)
end

function _apply_gaussian_refinement_mask(
    mask::_TernaryGaussianRefinementMask1D,
    coefficients::AbstractVector{<:Real};
    offset::Integer = 0,
)
    T = promote_type(Float64, eltype(coefficients))
    ncoarse = length(coefficients)
    fine_offset = mask.factor * Int(offset) - mask.support_radius
    nfine = ncoarse == 0 ? 0 : mask.factor * (ncoarse - 1) + 1 + 2 * mask.support_radius
    fine_coefficients = zeros(T, nfine)
    if ncoarse == 0
        return (offset = fine_offset, coefficients = fine_coefficients)
    end

    for coarse_index in eachindex(coefficients)
        coarse_site = Int(offset) + coarse_index - 1
        center_site = mask.factor * coarse_site
        amplitude = coefficients[coarse_index]
        for (local_index, local_offset) in enumerate(_refinement_mask_offsets(mask))
            fine_index = center_site + local_offset - fine_offset + 1
            fine_coefficients[fine_index] += amplitude * mask.coefficients[local_index]
        end
    end
    return (offset = fine_offset, coefficients = fine_coefficients)
end

function _apply_gaussian_refinement_mask(
    mask::_TernaryGaussianRefinementMask1D,
    line::NamedTuple{(:offset, :coefficients)},
)
    return _apply_gaussian_refinement_mask(mask, line.coefficients; offset = line.offset)
end

function _apply_gaussian_refinement_mask_repeated(
    mask::_TernaryGaussianRefinementMask1D,
    coefficients::AbstractVector{<:Real};
    levels::Integer,
    offset::Integer = 0,
)
    levels >= 0 || throw(ArgumentError("repeated refinement requires levels >= 0"))
    line = (offset = Int(offset), coefficients = collect(coefficients))
    for _ in 1:levels
        line = _apply_gaussian_refinement_mask(mask, line)
    end
    return line
end

function _apply_gaussian_refinement_mask_repeated(
    mask::_TernaryGaussianRefinementMask1D,
    line::NamedTuple{(:offset, :coefficients)},
    ;
    levels::Integer,
)
    return _apply_gaussian_refinement_mask_repeated(mask, line.coefficients; levels = levels, offset = line.offset)
end
