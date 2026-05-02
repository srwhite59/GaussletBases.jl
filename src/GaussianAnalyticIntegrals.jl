module GaussianAnalyticIntegrals

function gaussian_overlap(a, b)
    sigma2 = a.width^2 + b.width^2
    prefactor = sqrt(2.0 * pi) * a.width * b.width / sqrt(sigma2)
    return prefactor * exp(-0.5 * (a.center_value - b.center_value)^2 / sigma2)
end

function gaussian_kinetic(a, b)
    sigma2 = a.width^2 + b.width^2
    overlap_value = gaussian_overlap(a, b)
    delta = a.center_value - b.center_value
    return 0.5 * overlap_value * (1.0 / sigma2) * (1.0 - delta^2 / sigma2)
end

function gaussian_position(a, b)
    sigma2 = a.width^2 + b.width^2
    weighted_center =
        (a.center_value * b.width^2 + b.center_value * a.width^2) / sigma2
    return weighted_center * gaussian_overlap(a, b)
end

function gaussian_x2(a, b)
    sigma2 = a.width^2 + b.width^2
    weighted_center =
        (a.center_value * b.width^2 + b.center_value * a.width^2) / sigma2
    variance = a.width^2 * b.width^2 / sigma2
    return (weighted_center^2 + variance) * gaussian_overlap(a, b)
end

function gaussian_factor(a, b, exponent::Float64, center_value::Float64)
    exponent >= 0.0 || throw(ArgumentError("gaussian factor integral requires exponent >= 0"))
    exponent == 0.0 && return gaussian_overlap(a, b)

    alpha_a = inv(a.width^2)
    alpha_b = inv(b.width^2)
    alpha_g = 2.0 * exponent
    total_alpha = alpha_a + alpha_b + alpha_g
    weighted_center =
        (
            alpha_a * a.center_value +
            alpha_b * b.center_value +
            alpha_g * center_value
        ) / total_alpha
    constant_term =
        alpha_a * a.center_value^2 +
        alpha_b * b.center_value^2 +
        alpha_g * center_value^2
    return sqrt(2.0 * pi / total_alpha) *
           exp(-0.5 * (constant_term - total_alpha * weighted_center^2))
end

function gaussian_pair_factor(a, b, exponent::Float64)
    exponent >= 0.0 || throw(ArgumentError("pair-factor exponent must be >= 0"))

    alpha_a = inv(a.width^2)
    alpha_b = inv(b.width^2)
    a11 = alpha_a + 2.0 * exponent
    a22 = alpha_b + 2.0 * exponent
    a12 = -2.0 * exponent
    determinant = a11 * a22 - a12^2
    determinant > 0.0 ||
        throw(ArgumentError("pair-factor quadratic form must be positive definite"))

    d1 = alpha_a * a.center_value
    d2 = alpha_b * b.center_value
    constant_term = alpha_a * a.center_value^2 + alpha_b * b.center_value^2
    quadratic_term = (a22 * d1^2 - 2.0 * a12 * d1 * d2 + a11 * d2^2) / determinant
    return (2.0 * pi / sqrt(determinant)) * exp(-0.5 * (constant_term - quadratic_term))
end

gaussian_exponent(gaussian) = 1.0 / (2.0 * gaussian.width^2)

function doublefactorial(n::Int)
    n <= 0 && return 1.0
    value = 1.0
    for k in n:-2:1
        value *= k
    end
    return value
end

function polynomial_gaussian_shell_prefactor(exponent::Float64, power::Int)
    power >= 0 || throw(ArgumentError("polynomial Gaussian shell prefactor requires power >= 0"))
    numerator = 2.0^(2 * power + 0.5) * exponent^(power + 0.5)
    denominator = sqrt(pi) * doublefactorial(2 * power - 1)
    return sqrt(numerator / denominator)
end

function shifted_gaussian_moment(gamma::Float64, power::Int)
    isodd(power) && return 0.0
    k = div(power, 2)
    return sqrt(pi) * doublefactorial(2 * k - 1) / (2.0^k * gamma^(k + 0.5))
end

function polynomial_shift_multiply(
    coefficients::Vector{Float64},
    shift::Float64,
    power::Int,
)
    power >= 0 || throw(ArgumentError("polynomial shift multiply requires power >= 0"))
    result = copy(coefficients)
    for _ in 1:power
        next = zeros(Float64, length(result) + 1)
        for degree in eachindex(result)
            value = result[degree]
            next[degree] += shift * value
            next[degree + 1] += value
        end
        result = next
    end
    return result
end

function polynomial_gaussian_basic_integral(
    alpha_left::Float64,
    center_left::Float64,
    power_left::Int,
    prefactor_left::Float64,
    alpha_right::Float64,
    center_right::Float64,
    power_right::Int,
    prefactor_right::Float64;
    xpower::Int = 0,
    extra_exponent::Float64 = 0.0,
    extra_center::Float64 = 0.0,
)
    gamma = alpha_left + alpha_right + extra_exponent
    gamma > 0.0 ||
        throw(ArgumentError("polynomial Gaussian integral requires positive total exponent"))
    weighted_center =
        (
            alpha_left * center_left +
            alpha_right * center_right +
            extra_exponent * extra_center
        ) / gamma
    constant =
        alpha_left * center_left^2 + alpha_right * center_right^2 +
        extra_exponent * extra_center^2 - gamma * weighted_center^2

    polynomial = Float64[1.0]
    polynomial = polynomial_shift_multiply(polynomial, weighted_center - center_left, power_left)
    polynomial = polynomial_shift_multiply(polynomial, weighted_center - center_right, power_right)
    xpower > 0 && (polynomial = polynomial_shift_multiply(polynomial, weighted_center, xpower))

    value = 0.0
    for degree in eachindex(polynomial)
        value += polynomial[degree] * shifted_gaussian_moment(gamma, degree - 1)
    end
    return prefactor_left * prefactor_right * exp(-constant) * value
end

function polynomial_gaussian_derivative_terms(power::Int, exponent::Float64)
    power == 0 && return ((1, -2.0 * exponent),)
    power == 1 && return ((0, 1.0), (2, -2.0 * exponent))
    power == 2 && return ((1, 2.0), (3, -2.0 * exponent))
    throw(
        ArgumentError(
            "polynomial Gaussian derivative terms currently support only powers 0, 1, and 2",
        ),
    )
end

function polynomial_gaussian_kinetic_integral(
    alpha_left::Float64,
    center_left::Float64,
    power_left::Int,
    prefactor_left::Float64,
    alpha_right::Float64,
    center_right::Float64,
    power_right::Int,
    prefactor_right::Float64,
)
    value = 0.0
    left_terms = polynomial_gaussian_derivative_terms(power_left, alpha_left)
    right_terms = polynomial_gaussian_derivative_terms(power_right, alpha_right)
    for (derived_left_power, left_scale) in left_terms
        for (derived_right_power, right_scale) in right_terms
            value += 0.5 * left_scale * right_scale * polynomial_gaussian_basic_integral(
                alpha_left,
                center_left,
                derived_left_power,
                prefactor_left,
                alpha_right,
                center_right,
                derived_right_power,
                prefactor_right,
            )
        end
    end
    return value
end

end
