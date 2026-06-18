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
