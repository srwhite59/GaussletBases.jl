"""
    AtomicIDATwoElectronState

One product-basis state in the first tiny two-electron IDA application layer.

The ordering is up-major:

- the spin-up orbital index runs slowest
- the spin-down orbital index runs fastest

This is the natural product ordering for a `1 up, 1 down` problem built on the
channel-major spatial orbital ordering of `AtomicIDAOperators`.
"""
struct AtomicIDATwoElectronState
    index::Int
    up_orbital::AtomicOrbital
    down_orbital::AtomicOrbital

    function AtomicIDATwoElectronState(index::Int, up_orbital::AtomicOrbital, down_orbital::AtomicOrbital)
        index >= 1 || throw(ArgumentError("AtomicIDATwoElectronState requires index >= 1"))
        return new(index, up_orbital, down_orbital)
    end
end

function Base.show(io::IO, state::AtomicIDATwoElectronState)
    print(
        io,
        "AtomicIDATwoElectronState(index=",
        state.index,
        ", up=",
        state.up_orbital,
        ", down=",
        state.down_orbital,
        ")",
    )
end

"""
    AtomicIDATwoElectronProblem

    Tiny interacting `1 up, 1 down` consumer built on top of `AtomicIDAOperators`.

The object keeps the first many-electron step intentionally explicit:

- product-basis state indexing
- orbital and two-electron overlap diagnostics
- the one-body contribution
- the IDA-style density-density two-body contribution
- the total dense Hamiltonian

It is meant to validate the interacting Hamiltonian assembly, not to provide a
general many-electron framework.
"""
struct AtomicIDATwoElectronProblem{A <: AbstractDiagonalApproximation}
    operators::AtomicIDAOperators{A}
    orbital_data::Vector{AtomicOrbital}
    orbital_overlap::Matrix{Float64}
    state_data::Vector{AtomicIDATwoElectronState}
    overlap::Matrix{Float64}
    one_body::Matrix{Float64}
    two_body::Matrix{Float64}
    hamiltonian::Matrix{Float64}
end

function Base.show(io::IO, problem::AtomicIDATwoElectronProblem)
    print(
        io,
        "AtomicIDATwoElectronProblem(norbitals=",
        length(problem.orbital_data),
        ", nstates=",
        length(problem.state_data),
        ", Lmax=",
        length(problem.operators.angular_kernel_data) - 1,
        ")",
    )
end

orbitals(problem::AtomicIDATwoElectronProblem) = problem.orbital_data
two_electron_states(problem::AtomicIDATwoElectronProblem) = problem.state_data

function apply_overlap(problem::AtomicIDATwoElectronProblem, coefficients::AbstractVector{<:Real})
    length(coefficients) == size(problem.overlap, 1) ||
        throw(DimensionMismatch("coefficient vector length must match the two-electron overlap dimension"))
    return problem.overlap * Vector{Float64}(coefficients)
end

function apply_hamiltonian(problem::AtomicIDATwoElectronProblem, coefficients::AbstractVector{<:Real})
    length(coefficients) == size(problem.hamiltonian, 1) ||
        throw(DimensionMismatch("coefficient vector length must match the two-electron Hamiltonian dimension"))
    return problem.hamiltonian * Vector{Float64}(coefficients)
end

function ground_state_energy(problem::AtomicIDATwoElectronProblem)
    eig = eigen(Hermitian(problem.hamiltonian))
    return minimum(real(eig.values))
end

function lanczos_ground_state(
    problem::AtomicIDATwoElectronProblem;
    krylovdim::Int = 200,
    maxiter::Int = 200,
    tol::Real = 1.0e-10,
    v0::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    n = size(problem.hamiltonian, 1)
    n >= 1 || throw(ArgumentError("lanczos_ground_state requires a nonempty problem"))
    krylovdim >= 2 || throw(ArgumentError("lanczos_ground_state requires krylovdim >= 2"))
    maxiter >= 1 || throw(ArgumentError("lanczos_ground_state requires maxiter >= 1"))
    tol > 0 || throw(ArgumentError("lanczos_ground_state requires tol > 0"))

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

    maxsteps = min(maxiter, krylovdim, n)
    for step in 1:maxsteps
        iterations = step
        w = apply_hamiltonian(problem, v)
        step > 1 && (w .-= beta[end] .* previous)

        a = dot(v, w)
        push!(alpha, a)
        w .-= a .* v

        # Full reorthogonalization keeps this small reference implementation robust.
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

function _atomic_ida_two_electron_states(orbital_data::AbstractVector{AtomicOrbital})
    states = AtomicIDATwoElectronState[]
    index = 1
    for up_orbital in orbital_data
        for down_orbital in orbital_data
            push!(states, AtomicIDATwoElectronState(index, up_orbital, down_orbital))
            index += 1
        end
    end
    return states
end

function _orbital_channel_indices(ops::AtomicIDAOperators, orbital_data::AbstractVector{AtomicOrbital})
    channels = ops.one_body.channels
    return [_channel_index(channels, orbital.channel) for orbital in orbital_data]
end

function _ida_density_interaction_matrix(ops::AtomicIDAOperators, orbital_data::AbstractVector{AtomicOrbital})
    norb = length(orbital_data)
    channel_indices = _orbital_channel_indices(ops, orbital_data)
    vmax = length(ops.angular_kernel_data) - 1
    matrix = zeros(Float64, norb, norb)

    for left in 1:norb
        left_orbital = orbital_data[left]
        alpha = channel_indices[left]
        p = left_orbital.radial_index

        for right in 1:norb
            right_orbital = orbital_data[right]
            beta = channel_indices[right]
            q = right_orbital.radial_index

            total = 0.0
            for L in 0:vmax
                total += radial_multipole(ops, L)[p, q] *
                         angular_kernel(ops, L)[alpha, alpha, beta, beta]
            end
            matrix[left, right] = total
        end
    end

    return 0.5 .* (matrix .+ transpose(matrix))
end

function _ida_two_body_matrix(ops::AtomicIDAOperators, orbital_data::AbstractVector{AtomicOrbital})
    density_interaction = _ida_density_interaction_matrix(ops, orbital_data)
    norb = length(orbital_data)
    matrix = zeros(Float64, norb^2, norb^2)

    for up in 1:norb
        for down in 1:norb
            row = (up - 1) * norb + down
            matrix[row, row] = density_interaction[up, down]
        end
    end

    return matrix
end

"""
    atomic_ida_two_electron_problem(ops::AtomicIDAOperators; orbital_indices=nothing)

Build the first tiny interacting `1 up, 1 down` IDA problem on top of
`AtomicIDAOperators`.

The spatial-orbital basis is the channel-major orbital basis already exposed by
`ops`, optionally restricted to an explicit active subset through
`orbital_indices`. The two-electron product basis is ordered with spin-up
orbital index running slowest and spin-down orbital index running fastest.

The result exposes:

- `problem.orbital_overlap`
- `problem.overlap`
- `problem.one_body`
- `problem.two_body`
- `problem.hamiltonian`

and can be used directly in an ordinary Hermitian solve. The overlap matrices
are retained as diagnostics.
"""
function atomic_ida_two_electron_problem(ops::AtomicIDAOperators; orbital_indices = nothing)
    all_orbitals = orbitals(ops)
    if orbital_indices === nothing
        selected_indices = collect(eachindex(all_orbitals))
    else
        selected_indices = Int[Int(index) for index in orbital_indices]
        isempty(selected_indices) && throw(ArgumentError("atomic_ida_two_electron_problem requires at least one selected orbital"))
        all(index in eachindex(all_orbitals) for index in selected_indices) ||
            throw(BoundsError(all_orbitals, selected_indices))
        length(unique(selected_indices)) == length(selected_indices) ||
            throw(ArgumentError("atomic_ida_two_electron_problem requires selected orbital indices to be unique"))
    end

    orbital_data = all_orbitals[selected_indices]
    norb = length(orbital_data)
    norb >= 1 || throw(ArgumentError("atomic_ida_two_electron_problem requires at least one orbital"))

    overlap_orbital = ops.one_body.overlap[selected_indices, selected_indices]
    one_body_orbital = ops.one_body.hamiltonian[selected_indices, selected_indices]

    identity_orbital = Matrix{Float64}(I, norb, norb)
    overlap = kron(overlap_orbital, overlap_orbital)
    one_body = kron(one_body_orbital, identity_orbital) +
               kron(identity_orbital, one_body_orbital)
    two_body = _ida_two_body_matrix(ops, orbital_data)
    hamiltonian = one_body + two_body
    state_data = _atomic_ida_two_electron_states(orbital_data)
    return AtomicIDATwoElectronProblem(ops, orbital_data, overlap_orbital, state_data, overlap, one_body, two_body, hamiltonian)
end
