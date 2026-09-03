function _cartesian_weighted_hadamard3(
    left_coefficients::AbstractVector{<:Real},
    right_coefficients::AbstractVector{<:Real},
    x::AbstractMatrix{<:Real},
    y::AbstractMatrix{<:Real},
    z::AbstractMatrix{<:Real},
)
    matrix = Matrix{Float64}(x) .* Matrix{Float64}(y) .* Matrix{Float64}(z)
    return Float64(dot(left_coefficients, matrix * right_coefficients))
end

function _cartesian_gaussian_axis_prefactors(
    exponents::AbstractVector{<:Real},
    powers::AbstractVector{<:Integer},
)
    length(exponents) == length(powers) ||
        throw(DimensionMismatch("Gaussian axis exponents and powers differ"))
    return Float64[
        GaussianAnalyticIntegrals.polynomial_gaussian_shell_prefactor(
            Float64(exponents[index]),
            Int(powers[index]),
        ) for index in eachindex(exponents)
    ]
end

function _cartesian_gaussian_axis_integral(
    exponent_left::Real,
    center_left::Real,
    power_left::Integer,
    prefactor_left::Real,
    exponent_right::Real,
    center_right::Real,
    power_right::Integer,
    prefactor_right::Real,
    term::Symbol;
    factor_exponent::Real = 0.0,
    factor_center::Real = 0.0,
)
    alpha_left = Float64(exponent_left)
    alpha_right = Float64(exponent_right)
    left_center = Float64(center_left)
    right_center = Float64(center_right)
    left_power = Int(power_left)
    right_power = Int(power_right)
    left_prefactor = Float64(prefactor_left)
    right_prefactor = Float64(prefactor_right)
    term === :kinetic && return GaussianAnalyticIntegrals.polynomial_gaussian_kinetic_integral(
        alpha_left,
        left_center,
        left_power,
        left_prefactor,
        alpha_right,
        right_center,
        right_power,
        right_prefactor,
    )
    term === :overlap && return GaussianAnalyticIntegrals.polynomial_gaussian_basic_integral(
        alpha_left,
        left_center,
        left_power,
        left_prefactor,
        alpha_right,
        right_center,
        right_power,
        right_prefactor,
    )
    xpower =
        term === :position ? 1 :
        term === :x2 ? 2 :
        term === :factor ? 0 :
        throw(ArgumentError("unsupported Cartesian Gaussian axis term :$(term)"))
    return GaussianAnalyticIntegrals.polynomial_gaussian_basic_integral(
        alpha_left,
        left_center,
        left_power,
        left_prefactor,
        alpha_right,
        right_center,
        right_power,
        right_prefactor;
        xpower,
        extra_exponent = term === :factor ? Float64(factor_exponent) : 0.0,
        extra_center = term === :factor ? Float64(factor_center) : 0.0,
    )
end

function _cartesian_gaussian_axis_integral_table(
    left_exponents::AbstractVector{<:Real},
    left_centers::AbstractVector{<:Real},
    left_powers::AbstractVector{<:Integer},
    left_prefactors::AbstractVector{<:Real},
    right_exponents::AbstractVector{<:Real},
    right_centers::AbstractVector{<:Real},
    right_powers::AbstractVector{<:Integer},
    right_prefactors::AbstractVector{<:Real},
    term::Symbol;
    factor_exponent::Real = 0.0,
    factor_center::Real = 0.0,
)
    left_count = length(left_exponents)
    right_count = length(right_exponents)
    length(left_centers) == left_count && length(left_powers) == left_count &&
        length(left_prefactors) == left_count ||
        throw(DimensionMismatch("left Gaussian axis arrays have inconsistent lengths"))
    length(right_centers) == right_count && length(right_powers) == right_count &&
        length(right_prefactors) == right_count ||
        throw(DimensionMismatch("right Gaussian axis arrays have inconsistent lengths"))
    matrix = zeros(Float64, left_count, right_count)
    for right in 1:right_count, left in 1:left_count
        matrix[left, right] = _cartesian_gaussian_axis_integral(
            left_exponents[left],
            left_centers[left],
            left_powers[left],
            left_prefactors[left],
            right_exponents[right],
            right_centers[right],
            right_powers[right],
            right_prefactors[right],
            term;
            factor_exponent,
            factor_center,
        )
    end
    return matrix
end
