function _lanczos_ground_state_apply(
    apply_hamiltonian!::Function,
    n::Int;
    krylovdim::Int = 200,
    maxiter::Int = 200,
    tol::Real = 1.0e-10,
    v0::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    n >= 1 || throw(ArgumentError("_lanczos_ground_state_apply requires a nonempty problem"))
    krylovdim >= 2 || throw(ArgumentError("_lanczos_ground_state_apply requires krylovdim >= 2"))
    maxiter >= 1 || throw(ArgumentError("_lanczos_ground_state_apply requires maxiter >= 1"))
    tol > 0 || throw(ArgumentError("_lanczos_ground_state_apply requires tol > 0"))

    if v0 === nothing
        v = ones(Float64, n)
    else
        length(v0) == n || throw(DimensionMismatch("initial Lanczos vector length must match the Hamiltonian dimension"))
        v = Vector{Float64}(v0)
    end

    norm(v) > 0 || throw(ArgumentError("initial Lanczos vector must be nonzero"))
    v ./= norm(v)

    vectors = Vector{Vector{Float64}}()
    push!(vectors, copy(v))
    alpha = Float64[]
    beta = Float64[]
    converged = false
    residual = Inf
    iterations = 0
    best_small_vector = ones(Float64, 1)
    best_value = NaN
    previous = zeros(Float64, n)
    scratch = zeros(Float64, n)

    maxsteps = min(maxiter, krylovdim, n)
    for step in 1:maxsteps
        iterations = step
        apply_hamiltonian!(scratch, v)
        w = copy(scratch)
        step > 1 && (w .-= beta[end] .* previous)

        a = dot(v, w)
        push!(alpha, a)
        w .-= a .* v

        for basis_vector in vectors
            w .-= dot(basis_vector, w) .* basis_vector
        end

        b = norm(w)
        small_eig = eigen(SymTridiagonal(alpha, beta))
        best_value = real(small_eig.values[1])
        best_small_vector = Vector{Float64}(small_eig.vectors[:, 1])
        residual = abs(b * best_small_vector[end])

        if residual <= tol || step == maxsteps || b <= sqrt(eps(Float64))
            converged = residual <= tol || b <= sqrt(eps(Float64))
            break
        end

        push!(beta, b)
        previous = v
        v = w ./ b
        push!(vectors, copy(v))
    end

    lanczos_vector = zeros(Float64, n)
    for j in eachindex(best_small_vector)
        lanczos_vector .+= best_small_vector[j] .* vectors[j]
    end
    lanczos_vector ./= norm(lanczos_vector)

    return (
        value = best_value,
        vector = lanczos_vector,
        residual = residual,
        iterations = iterations,
        converged = converged,
    )
end
