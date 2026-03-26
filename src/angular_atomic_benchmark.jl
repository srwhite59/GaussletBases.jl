"""
    AtomicInjectedAngularOneBodyBenchmark

Experimental atom-side one-electron angular benchmark built from a radial
one-body substrate and an `AtomicShellLocalInjectedAngularAssembly`.

This is the first benchmark-layer import on the angular research track. It
keeps the scope narrow: exact shell-pair Galerkin assembly for the one-electron
central atomic Hamiltonian, plus the exact injected low-`l` comparison surface.
"""
struct AtomicInjectedAngularOneBodyBenchmark{R <: RadialAtomicOperators}
    radial_operators::R
    angular_assembly::AtomicShellLocalInjectedAngularAssembly
    shell_ranges::Vector{UnitRange{Int}}
    overlap::Matrix{Float64}
    hamiltonian::Matrix{Float64}
    inverse_r2::Matrix{Float64}
    exact_common_lmax::Int
    exact_channels::YlmChannelSet
    exact_channel_ranges::Vector{UnitRange{Int}}
    exact_transform::Matrix{Float64}
    exact_overlap::Matrix{Float64}
    exact_hamiltonian::Matrix{Float64}
end

"""
    AtomicInjectedAngularHFStyleBenchmark

Experimental closed-shell density-density SCF benchmark on top of the injected
angular assembly layer.

This stays intentionally narrow. It uses the assembled one-electron operators
from `AtomicInjectedAngularOneBodyBenchmark`, a two-index density-density
interaction assembled from shell moment blocks, and a simple closed-shell SCF
iteration. It is an HF-style benchmark path, not a general HF framework.
"""
struct AtomicInjectedAngularHFStyleBenchmark{
    B <: AtomicInjectedAngularOneBodyBenchmark,
    I,
}
    one_body::B
    interaction::Matrix{Float64}
    projected_exact_interaction::Matrix{Float64}
    exact_ida_reference::I
    exact_interaction::Matrix{Float64}
    scf_result::NamedTuple
    exact_scf_result::NamedTuple
end

"""
    AtomicInjectedAngularHFDMRGHFAdapter

Experimental in-memory HFDMRG-facing HF adapter for the current angular
benchmark line.

This adapter is intentionally narrow. It reuses the shell-assembled one-body
and interaction matrices already present in
`AtomicInjectedAngularHFStyleBenchmark` and packages them into the dense
`H, V, psiup0, psidn0` handshake expected by `HFDMRG.solve_hfdmrg(...)`.

It is an HF adapter surface, not a general many-body or file-export contract.
"""
struct AtomicInjectedAngularHFDMRGHFAdapter{
    H <: AtomicInjectedAngularHFStyleBenchmark,
}
    hf_style::H
    route::Symbol
    hamiltonian::Matrix{Float64}
    interaction::Matrix{Float64}
    psiup0::Matrix{Float64}
    psidn0::Matrix{Float64}
    psiup0_source::Symbol
    psidn0_source::Symbol
    nup::Int
    ndn::Int
    overlap_identity_error::Float64
end

"""
    AtomicInjectedAngularSmallEDBenchmark

Experimental tiny `1 up, 1 down` density-density benchmark on top of the
injected angular assembly line.

This is the next narrow ladder step after the first angular HF-style
benchmark. It reuses the same shell-assembled one-body and interaction
ingredients and solves the closed-shell He reference-sized product problem with
a matrix-free Lanczos kernel. It is a benchmark path, not a general ED
framework.
"""
struct AtomicInjectedAngularSmallEDBenchmark{
    H <: AtomicInjectedAngularHFStyleBenchmark,
    P,
}
    hf_style::H
    orbital_count::Int
    state_count::Int
    orbital_overlap::Matrix{Float64}
    orbital_one_body::Matrix{Float64}
    orbital_interaction::Matrix{Float64}
    exact_reference_problem::P
    lanczos_result::NamedTuple
end

function Base.show(io::IO, benchmark::AtomicInjectedAngularSmallEDBenchmark)
    print(
        io,
        "AtomicInjectedAngularSmallEDBenchmark(norbitals=",
        benchmark.orbital_count,
        ", nstates=",
        benchmark.state_count,
        ", exact_common_lmax=",
        benchmark.hf_style.one_body.exact_common_lmax,
        ")",
    )
end

function Base.show(io::IO, adapter::AtomicInjectedAngularHFDMRGHFAdapter)
    print(
        io,
        "AtomicInjectedAngularHFDMRGHFAdapter(route=",
        adapter.route,
        ", dim=",
        size(adapter.hamiltonian, 1),
        ", nup=",
        adapter.nup,
        ", ndn=",
        adapter.ndn,
        ", psiup0_source=",
        adapter.psiup0_source,
        ", psidn0_source=",
        adapter.psidn0_source,
        ", exact_common_lmax=",
        adapter.hf_style.one_body.exact_common_lmax,
        ")",
    )
end

function Base.show(io::IO, benchmark::AtomicInjectedAngularHFStyleBenchmark)
    print(
        io,
        "AtomicInjectedAngularHFStyleBenchmark(nshells=",
        length(benchmark.one_body.angular_assembly.shells),
        ", total_dim=",
        size(benchmark.one_body.overlap, 1),
        ", exact_common_lmax=",
        benchmark.one_body.exact_common_lmax,
        ", converged=",
        benchmark.scf_result.converged,
        ")",
    )
end

function Base.show(io::IO, benchmark::AtomicInjectedAngularOneBodyBenchmark)
    print(
        io,
        "AtomicInjectedAngularOneBodyBenchmark(nshells=",
        length(benchmark.angular_assembly.shells),
        ", total_dim=",
        size(benchmark.overlap, 1),
        ", exact_common_lmax=",
        benchmark.exact_common_lmax,
        ", shell_orders=",
        repr(benchmark.angular_assembly.shell_orders),
        ")",
    )
end

function _orthonormalize_columns(columns::AbstractMatrix{<:Real})
    matrix = Matrix{Float64}(columns)
    ncols = size(matrix, 2)
    q = Matrix(qr(matrix).Q)
    return q[:, 1:ncols]
end

function _validated_hfdmrg_seed(
    seed::AbstractMatrix{<:Real},
    norb::Int,
    nspin::Int,
    name::AbstractString;
    tol::Real = 1.0e-10,
)
    size(seed) == (norb, nspin) ||
        throw(
            DimensionMismatch(
                "$name must have size ($norb, $nspin), got $(size(seed))",
            ),
        )
    matrix = Matrix{Float64}(seed)
    nspin == 0 && return matrix
    orthogonality_error =
        opnorm(transpose(matrix) * matrix - Matrix{Float64}(I, nspin, nspin), Inf)
    orthogonality_error <= tol ||
        throw(
            ArgumentError(
                "$name must have orthonormal columns for HFDMRG; got orthogonality error $orthogonality_error",
            ),
        )
    return matrix
end

function _validate_hfdmrg_spin_count(nspin::Int, norb::Int, name::AbstractString)
    0 <= nspin <= norb ||
        throw(ArgumentError("$name must satisfy 0 <= $name <= $norb"))
    return nspin
end

function _resolve_hfdmrg_module(hfmod)
    hfmod !== nothing && return hfmod
    if isdefined(Main, :HFDMRG)
        return getfield(Main, :HFDMRG)
    end
    error(
        "HFDMRG module not loaded. Add it to LOAD_PATH and run `using HFDMRG` before calling `run_atomic_injected_angular_hfdmrg_hf`.",
    )
end

_default_angular_hfdmrg_spin_count(benchmark::AtomicInjectedAngularHFStyleBenchmark) =
    size(benchmark.scf_result.occupied_coefficients, 2)

_default_angular_hfdmrg_spin_count(benchmark::AtomicInjectedAngularSmallEDBenchmark) =
    _default_angular_hfdmrg_spin_count(benchmark.hf_style)

function _matrix_ranges(offsets::Vector{Int}, dims::Vector{Int})
    return [offsets[i]:(offsets[i] + dims[i] - 1) for i in eachindex(dims)]
end

function _generalized_spectrum(
    hamiltonian::AbstractMatrix{<:Real},
    overlap::AbstractMatrix{<:Real};
    tau::Real = 1.0e-12,
)
    overlap_eigen = eigen(Symmetric(Matrix{Float64}(overlap)))
    max_eval = max(1.0, maximum(abs.(overlap_eigen.values)))
    cutoff = Float64(tau) * max_eval
    any(value -> value < -cutoff, overlap_eigen.values) &&
        throw(ArgumentError("overlap has a negative eigenvalue below the stability cutoff"))

    repaired = max.(overlap_eigen.values, cutoff)
    overlap_inv_sqrt =
        overlap_eigen.vectors * Diagonal(1.0 ./ sqrt.(repaired)) * transpose(overlap_eigen.vectors)
    orthogonal_hamiltonian =
        overlap_inv_sqrt * Matrix{Float64}(hamiltonian) * overlap_inv_sqrt
    return sort(eigvals(Symmetric(orthogonal_hamiltonian)))
end

function _generalized_eigensystem(
    hamiltonian::AbstractMatrix{<:Real},
    overlap::AbstractMatrix{<:Real};
    tau::Real = 1.0e-12,
)
    overlap_eigen = eigen(Symmetric(Matrix{Float64}(overlap)))
    max_eval = max(1.0, maximum(abs.(overlap_eigen.values)))
    cutoff = Float64(tau) * max_eval
    any(value -> value < -cutoff, overlap_eigen.values) &&
        throw(ArgumentError("overlap has a negative eigenvalue below the stability cutoff"))

    repaired = max.(overlap_eigen.values, cutoff)
    overlap_inv_sqrt =
        overlap_eigen.vectors * Diagonal(1.0 ./ sqrt.(repaired)) * transpose(overlap_eigen.vectors)
    orthogonal_hamiltonian =
        overlap_inv_sqrt * Matrix{Float64}(hamiltonian) * overlap_inv_sqrt
    eig = eigen(Symmetric(orthogonal_hamiltonian))
    coefficients = overlap_inv_sqrt * eig.vectors
    return (
        values = Vector{Float64}(eig.values),
        coefficients = Matrix{Float64}(coefficients),
    )
end

function _exact_channel_ranges(nshells::Int, channels::YlmChannelSet)
    nchannels = length(channels)
    return [((i - 1) * nchannels + 1):(i * nchannels) for i in 1:nshells]
end

function _shell_major_orbital_permutation(radial_dim::Int, channels::YlmChannelSet)
    nchannels = length(channels)
    perm = Int[]
    for radial_index in 1:radial_dim
        for channel_index in 1:nchannels
            push!(perm, (channel_index - 1) * radial_dim + radial_index)
        end
    end
    return perm
end

function _assembly_exact_channel_transform(
    assembly::AtomicShellLocalInjectedAngularAssembly,
    exact_channels::YlmChannelSet,
)
    nshells = length(assembly.shells)
    nchannels = length(exact_channels)
    transform = zeros(Float64, nshells * nchannels, size(assembly.overlap, 1))
    shell_ranges = _matrix_ranges(assembly.shell_offsets, assembly.shell_dimensions)

    for shell_index in 1:nshells
        shell = assembly.shells[shell_index]
        shell_range = shell_ranges[shell_index]
        lookup = Dict{YlmChannel,Int}(shell.injected_channels.channel_data[i] => i for i in eachindex(shell.injected_channels.channel_data))
        for channel_index in eachindex(exact_channels.channel_data)
            injected_row = get(lookup, exact_channels[channel_index], 0)
            injected_row == 0 &&
                throw(
                    ArgumentError(
                        "shell $(shell_index) does not contain exact injected channel $(exact_channels[channel_index])",
                    ),
                )
            row = (shell_index - 1) * nchannels + channel_index
            transform[row, shell_range] .= shell.injected_overlap[injected_row, :]
        end
    end

    return transform
end

function _exact_channel_overlap_matrix(
    radial_ops::RadialAtomicOperators,
    exact_channels::YlmChannelSet,
)
    nchannels = length(exact_channels)
    return kron(radial_ops.overlap, Matrix{Float64}(I, nchannels, nchannels))
end

function _exact_channel_hamiltonian_matrix(
    radial_ops::RadialAtomicOperators,
    exact_channels::YlmChannelSet,
    inverse_r2::AbstractMatrix{<:Real},
)
    nchannels = length(exact_channels)
    radial_one_body = radial_ops.kinetic + radial_ops.nuclear
    angular_kinetic =
        Diagonal([_kinetic_eigenvalue(channel) for channel in exact_channels.channel_data])
    return kron(radial_one_body, Matrix{Float64}(I, nchannels, nchannels)) +
           kron(Matrix{Float64}(inverse_r2), Matrix{Float64}(angular_kinetic))
end

function _shell_angular_weights(moment_blocks::Dict{Int,Matrix{Float64}})
    haskey(moment_blocks, 0) ||
        throw(ArgumentError("shell-local moment blocks must contain the monopole L=0 block"))
    weights = sqrt(4 * pi) .* vec(moment_blocks[0][1, :])
    any(weight -> abs(weight) <= 1.0e-12, weights) &&
        throw(ArgumentError("shell-local angular weights are too small for density-density assembly"))
    return weights
end

function _scaled_shell_moment_block(moment_blocks::Dict{Int,Matrix{Float64}}, L::Int)
    haskey(moment_blocks, L) || return nothing
    weights = _shell_angular_weights(moment_blocks)
    return moment_blocks[L] ./ reshape(weights, 1, :)
end

function _assemble_atomic_injected_angular_interaction(
    radial_ops::RadialAtomicOperators,
    assembly::AtomicShellLocalInjectedAngularAssembly,
)
    nshells = length(assembly.shells)
    shell_ranges = _matrix_ranges(assembly.shell_offsets, assembly.shell_dimensions)
    interaction = zeros(Float64, size(assembly.overlap))
    Lmax = length(radial_ops.multipole_data) - 1

    for L in 0:Lmax
        prefactor = 4 * pi / (2 * L + 1)
        radial_block = multipole(radial_ops, L)
        for a in 1:nshells
            ia = shell_ranges[a]
            left = _scaled_shell_moment_block(assembly.shell_moment_blocks[a], L)
            left === nothing && continue
            for b in 1:nshells
                ib = shell_ranges[b]
                right = _scaled_shell_moment_block(assembly.shell_moment_blocks[b], L)
                right === nothing && continue
                interaction[ia, ib] .+= prefactor .* radial_block[a, b] .* (transpose(left) * right)
            end
        end
    end

    return 0.5 .* (interaction .+ transpose(interaction))
end

function _closed_shell_density_density_scf(
    one_body::AbstractMatrix{<:Real},
    overlap::AbstractMatrix{<:Real},
    interaction::AbstractMatrix{<:Real};
    nelec::Int = 2,
    maxiter::Int = 50,
    damping::Real = 0.25,
    tol::Real = 1.0e-8,
)
    iseven(nelec) || throw(ArgumentError("_closed_shell_density_density_scf requires an even electron count"))
    maxiter >= 1 || throw(ArgumentError("_closed_shell_density_density_scf requires maxiter >= 1"))
    0 <= damping < 1 || throw(ArgumentError("_closed_shell_density_density_scf requires 0 <= damping < 1"))
    tol > 0 || throw(ArgumentError("_closed_shell_density_density_scf requires tol > 0"))

    norb = size(one_body, 1)
    size(one_body) == (norb, norb) || throw(DimensionMismatch("one_body must be square"))
    size(overlap) == (norb, norb) || throw(DimensionMismatch("overlap must have the same size as one_body"))
    size(interaction) == (norb, norb) || throw(DimensionMismatch("interaction must have the same size as one_body"))

    nocc = nelec ÷ 2
    nocc <= norb || throw(ArgumentError("occupied closed-shell orbital count cannot exceed the basis dimension"))

    initial = _generalized_eigensystem(one_body, overlap)
    density = 2.0 .* (initial.coefficients[:, 1:nocc] * transpose(initial.coefficients[:, 1:nocc]))
    density = 0.5 .* (density .+ transpose(density))

    energies = Float64[]
    residuals = Float64[]
    energy_change = Inf
    converged = false
    iterations = 0

    result = nothing
    for iteration in 1:maxiter
        iterations = iteration
        occupations = vec(diag(density * overlap))
        fock = Matrix{Float64}(one_body) + Diagonal(Matrix{Float64}(interaction) * occupations)
        eig = _generalized_eigensystem(fock, overlap)
        occupied = eig.coefficients[:, 1:nocc]
        density_new = 2.0 .* (occupied * transpose(occupied))
        density_new = 0.5 .* (density_new .+ transpose(density_new))
        mixed_density = (1 - damping) .* density_new .+ damping .* density
        mixed_density = 0.5 .* (mixed_density .+ transpose(mixed_density))

        mixed_occupations = vec(diag(mixed_density * overlap))
        energy = sum(Matrix{Float64}(one_body) .* mixed_density) +
                 0.5 * dot(mixed_occupations, Matrix{Float64}(interaction) * mixed_occupations)
        residual = norm(mixed_density - density, Inf)
        energy_change = isempty(energies) ? Inf : abs(energy - energies[end])
        push!(energies, energy)
        push!(residuals, residual)

        density = mixed_density
        result = (
            density = mixed_density,
            occupations = mixed_occupations,
            fock = 0.5 .* (fock .+ transpose(fock)),
            coefficients = eig.coefficients,
            occupied_coefficients = occupied,
            orbital_energies = eig.values,
            energy = energy,
        )

        if residual <= tol && energy_change <= tol
            converged = true
            break
        end
    end

    result === nothing && error("_closed_shell_density_density_scf failed to produce an SCF iterate")
    electron_count = sum(result.occupations)
    return merge(
        result,
        (
            energies = energies,
            residuals = residuals,
            iterations = iterations,
            converged = converged,
            energy_change = energy_change,
            electron_count = electron_count,
            electron_count_error = abs(electron_count - nelec),
            max_occupation = maximum(result.occupations),
            min_occupation = minimum(result.occupations),
        ),
    )
end

"""
    build_atomic_injected_angular_one_body_benchmark(radial_ops::RadialAtomicOperators;
                                                     shell_orders=nothing,
                                                     beta=2.0,
                                                     l_inject=:auto,
                                                     tau=1e-12,
                                                     whiten=:svd,
                                                     ord_min=minimum(curated_sphere_point_set_orders()),
                                                     ord_max=maximum(curated_sphere_point_set_orders()),
                                                     r_lo=0.15,
                                                     r_hi=7.0,
                                                     w_lo=0.2,
                                                     w_hi=0.7)

Build the first atom-side angular one-electron benchmark on top of the
shell-local injected angular assembly layer.

The benchmark uses the shell centers already present in `radial_ops` as the
radial-shell schedule and assembles the one-electron central Hamiltonian

    T_r + V_nuc + (1 / r^2) L^2

with exact shell-pair angular overlap and kinetic blocks.
"""
function build_atomic_injected_angular_one_body_benchmark(
    radial_ops::RadialAtomicOperators;
    shell_orders::Union{Nothing,AbstractVector{<:Integer}} = nothing,
    beta::Real = 2.0,
    l_inject::Union{Int,Symbol} = :auto,
    tau::Real = 1.0e-12,
    whiten::Symbol = :svd,
    ord_min::Int = minimum(curated_sphere_point_set_orders()),
    ord_max::Int = maximum(curated_sphere_point_set_orders()),
    r_lo::Real = 0.15,
    r_hi::Real = 7.0,
    w_lo::Real = 0.2,
    w_hi::Real = 0.7,
)
    assembly = build_atomic_shell_local_angular_assembly(
        radial_ops.shell_centers_r;
        shell_orders = shell_orders,
        beta = beta,
        l_inject = l_inject,
        tau = tau,
        whiten = whiten,
        ord_min = ord_min,
        ord_max = ord_max,
        r_lo = r_lo,
        r_hi = r_hi,
        w_lo = w_lo,
        w_hi = w_hi,
    )
    return build_atomic_injected_angular_one_body_benchmark(radial_ops, assembly)
end

"""
    build_atomic_injected_angular_one_body_benchmark(radial_ops::RadialAtomicOperators,
                                                     assembly::AtomicShellLocalInjectedAngularAssembly)

Assemble the first atom-side angular one-electron benchmark from an explicit
radial one-body substrate and a prebuilt shell-local angular assembly.
"""
function build_atomic_injected_angular_one_body_benchmark(
    radial_ops::RadialAtomicOperators,
    assembly::AtomicShellLocalInjectedAngularAssembly,
)
    radial_dim = size(radial_ops.overlap, 1)
    nshells = length(assembly.shells)
    radial_dim == nshells ||
        throw(
            ArgumentError(
                "radial dimension $radial_dim must match the number of angular shells $nshells",
            ),
        )

    shell_ranges = _matrix_ranges(assembly.shell_offsets, assembly.shell_dimensions)
    overlap = zeros(Float64, size(assembly.overlap))
    hamiltonian = zeros(Float64, size(assembly.overlap))
    inverse_r2 = centrifugal(radial_ops, 1)
    radial_one_body = radial_ops.kinetic + radial_ops.nuclear

    for a in 1:nshells
        ia = shell_ranges[a]
        for b in 1:nshells
            ib = shell_ranges[b]
            overlap[ia, ib] .= radial_ops.overlap[a, b] .* assembly.overlap_blocks[a, b]
            hamiltonian[ia, ib] .=
                radial_one_body[a, b] .* assembly.overlap_blocks[a, b] .+
                inverse_r2[a, b] .* assembly.kinetic_blocks[a, b]
        end
    end

    exact_common_lmax = minimum(assembly.shell_exact_lmax)
    exact_channels = ylm_channels(exact_common_lmax)
    exact_channel_ranges = _exact_channel_ranges(nshells, exact_channels)
    exact_transform = _assembly_exact_channel_transform(assembly, exact_channels)
    exact_overlap = _exact_channel_overlap_matrix(radial_ops, exact_channels)
    exact_hamiltonian = _exact_channel_hamiltonian_matrix(radial_ops, exact_channels, inverse_r2)

    return AtomicInjectedAngularOneBodyBenchmark(
        radial_ops,
        assembly,
        shell_ranges,
        overlap,
        hamiltonian,
        inverse_r2,
        exact_common_lmax,
        exact_channels,
        exact_channel_ranges,
        exact_transform,
        exact_overlap,
        exact_hamiltonian,
    )
end

"""
    atomic_injected_angular_one_body_diagnostics(benchmark::AtomicInjectedAngularOneBodyBenchmark)

Return compact diagnostics for the first atom-side angular one-electron
benchmark.
"""
function atomic_injected_angular_one_body_diagnostics(
    benchmark::AtomicInjectedAngularOneBodyBenchmark,
)
    projected_exact_overlap =
        benchmark.exact_transform * benchmark.overlap * transpose(benchmark.exact_transform)
    projected_exact_hamiltonian =
        benchmark.exact_transform * benchmark.hamiltonian * transpose(benchmark.exact_transform)
    benchmark_eigenvalues = _generalized_spectrum(benchmark.hamiltonian, benchmark.overlap)
    exact_eigenvalues = _generalized_spectrum(benchmark.exact_hamiltonian, benchmark.exact_overlap)
    projected_exact_eigenvalues =
        _generalized_spectrum(projected_exact_hamiltonian, projected_exact_overlap)
    ncompare = min(4, length(projected_exact_eigenvalues), length(exact_eigenvalues))

    return (
        nshells = length(benchmark.angular_assembly.shells),
        total_dim = size(benchmark.overlap, 1),
        shell_orders = copy(benchmark.angular_assembly.shell_orders),
        shell_dimensions = copy(benchmark.angular_assembly.shell_dimensions),
        exact_common_lmax = benchmark.exact_common_lmax,
        exact_channel_count = length(benchmark.exact_channels),
        overlap_symmetry_error = opnorm(benchmark.overlap - transpose(benchmark.overlap), Inf),
        hamiltonian_symmetry_error = opnorm(benchmark.hamiltonian - transpose(benchmark.hamiltonian), Inf),
        min_overlap_eigenvalue = minimum(eigvals(Symmetric(benchmark.overlap))),
        projected_exact_overlap_error = opnorm(projected_exact_overlap - benchmark.exact_overlap, Inf),
        projected_exact_hamiltonian_error = opnorm(projected_exact_hamiltonian - benchmark.exact_hamiltonian, Inf),
        benchmark_ground_state_energy = benchmark_eigenvalues[1],
        exact_ground_state_energy = exact_eigenvalues[1],
        benchmark_ground_state_error = abs(benchmark_eigenvalues[1] - exact_eigenvalues[1]),
        projected_exact_low_eigenvalue_count = ncompare,
        projected_exact_low_eigenvalues = projected_exact_eigenvalues[1:ncompare],
        exact_low_eigenvalues = exact_eigenvalues[1:ncompare],
        projected_exact_low_eigenvalue_error =
            maximum(abs.(projected_exact_eigenvalues[1:ncompare] .- exact_eigenvalues[1:ncompare])),
    )
end

function _apply_density_density_two_electron_hamiltonian!(
    out::AbstractVector{Float64},
    coefficients::AbstractVector{<:Real},
    one_body_orbital::AbstractMatrix{<:Real},
    interaction::AbstractMatrix{<:Real},
)
    norb = size(one_body_orbital, 1)
    length(coefficients) == norb^2 ||
        throw(DimensionMismatch("coefficient vector length must match the product-space dimension"))
    X = reshape(Vector{Float64}(coefficients), norb, norb)
    H = Matrix{Float64}(one_body_orbital)
    V = transpose(Matrix{Float64}(interaction))
    Y = H * X + X * transpose(H) + V .* X
    out .= vec(Y)
    return out
end

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

"""
    build_atomic_injected_angular_hf_style_benchmark(radial_ops::RadialAtomicOperators;
                                                     shell_orders=nothing,
                                                     beta=2.0,
                                                     l_inject=:auto,
                                                     tau=1e-12,
                                                     whiten=:svd,
                                                     nelec::Int=2,
                                                     maxiter::Int=50,
                                                     damping::Real=0.25,
                                                     tol::Real=1e-8,
                                                     ord_min=minimum(curated_sphere_point_set_orders()),
                                                     ord_max=maximum(curated_sphere_point_set_orders()),
                                                     r_lo=0.15,
                                                     r_hi=7.0,
                                                     w_lo=0.2,
                                                     w_hi=0.7)

Build the first narrow angular HF-style benchmark on top of the shell-local
angular assembly layer. The mean-field model is the current density-density
closed-shell SCF benchmark, not a full general Hartree-Fock framework.
"""
function build_atomic_injected_angular_hf_style_benchmark(
    radial_ops::RadialAtomicOperators;
    shell_orders::Union{Nothing,AbstractVector{<:Integer}} = nothing,
    beta::Real = 2.0,
    l_inject::Union{Int,Symbol} = :auto,
    tau::Real = 1.0e-12,
    whiten::Symbol = :svd,
    nelec::Int = 2,
    maxiter::Int = 50,
    damping::Real = 0.25,
    tol::Real = 1.0e-8,
    ord_min::Int = minimum(curated_sphere_point_set_orders()),
    ord_max::Int = maximum(curated_sphere_point_set_orders()),
    r_lo::Real = 0.15,
    r_hi::Real = 7.0,
    w_lo::Real = 0.2,
    w_hi::Real = 0.7,
)
    one_body = build_atomic_injected_angular_one_body_benchmark(
        radial_ops;
        shell_orders = shell_orders,
        beta = beta,
        l_inject = l_inject,
        tau = tau,
        whiten = whiten,
        ord_min = ord_min,
        ord_max = ord_max,
        r_lo = r_lo,
        r_hi = r_hi,
        w_lo = w_lo,
        w_hi = w_hi,
    )
    return build_atomic_injected_angular_hf_style_benchmark(
        radial_ops,
        one_body;
        nelec = nelec,
        maxiter = maxiter,
        damping = damping,
        tol = tol,
    )
end

function build_atomic_injected_angular_hf_style_benchmark(
    radial_ops::RadialAtomicOperators,
    one_body::AtomicInjectedAngularOneBodyBenchmark;
    nelec::Int = 2,
    maxiter::Int = 50,
    damping::Real = 0.25,
    tol::Real = 1.0e-8,
)
    interaction =
        _assemble_atomic_injected_angular_interaction(radial_ops, one_body.angular_assembly)
    projected_exact_interaction =
        one_body.exact_transform * interaction * transpose(one_body.exact_transform)
    exact_ida_reference = atomic_ida_operators(radial_ops; lmax = one_body.exact_common_lmax)
    exact_interaction =
        atomic_ida_density_interaction_matrix(exact_ida_reference; ordering = :shell_major)

    scf_result = _closed_shell_density_density_scf(
        one_body.hamiltonian,
        one_body.overlap,
        interaction;
        nelec = nelec,
        maxiter = maxiter,
        damping = damping,
        tol = tol,
    )
    exact_scf_result = _closed_shell_density_density_scf(
        one_body.exact_hamiltonian,
        one_body.exact_overlap,
        exact_interaction;
        nelec = nelec,
        maxiter = maxiter,
        damping = damping,
        tol = tol,
    )

    return AtomicInjectedAngularHFStyleBenchmark(
        one_body,
        interaction,
        projected_exact_interaction,
        exact_ida_reference,
        exact_interaction,
        scf_result,
        exact_scf_result,
    )
end

"""
    atomic_injected_angular_hf_style_diagnostics(benchmark::AtomicInjectedAngularHFStyleBenchmark)

Return compact diagnostics for the first angular HF-style benchmark.
"""
function atomic_injected_angular_hf_style_diagnostics(
    benchmark::AtomicInjectedAngularHFStyleBenchmark,
)
    full = benchmark.scf_result
    exact = benchmark.exact_scf_result
    ncompare = min(4, length(full.orbital_energies), length(exact.orbital_energies))

    return (
        nshells = length(benchmark.one_body.angular_assembly.shells),
        total_dim = size(benchmark.one_body.overlap, 1),
        shell_orders = copy(benchmark.one_body.angular_assembly.shell_orders),
        exact_common_lmax = benchmark.one_body.exact_common_lmax,
        interaction_symmetry_error = opnorm(benchmark.interaction - transpose(benchmark.interaction), Inf),
        full_converged = full.converged,
        full_iterations = full.iterations,
        full_energy = full.energy,
        full_residual = isempty(full.residuals) ? Inf : full.residuals[end],
        full_electron_count_error = full.electron_count_error,
        exact_converged = exact.converged,
        exact_iterations = exact.iterations,
        exact_energy = exact.energy,
        exact_residual = isempty(exact.residuals) ? Inf : exact.residuals[end],
        exact_electron_count_error = exact.electron_count_error,
        energy_difference_to_exact_reference = full.energy - exact.energy,
        ground_orbital_energy_error = abs(full.orbital_energies[1] - exact.orbital_energies[1]),
        orbital_energy_count = ncompare,
        full_low_orbital_energies = full.orbital_energies[1:ncompare],
        exact_low_orbital_energies = exact.orbital_energies[1:ncompare],
        exact_reference_low_orbital_error =
            maximum(abs.(full.orbital_energies[1:ncompare] .- exact.orbital_energies[1:ncompare])),
    )
end

"""
    build_atomic_injected_angular_hfdmrg_hf_adapter(
        radial_ops::RadialAtomicOperators;
        shell_orders=nothing,
        beta=2.0,
        l_inject=:auto,
        tau=1e-12,
        whiten=:svd,
        nelec::Int=2,
        maxiter::Int=50,
        damping::Real=0.25,
        tol::Real=1e-8,
        ord_min=minimum(curated_sphere_point_set_orders()),
        ord_max=maximum(curated_sphere_point_set_orders()),
        r_lo=0.15,
        r_hi=7.0,
        w_lo=0.2,
        w_hi=0.7,
        nup=nothing,
        ndn=nothing,
        psiup0=nothing,
        psidn0=nothing,
    )

Build the first in-memory HFDMRG-facing HF adapter on top of the angular
benchmark line.

The current adapter uses the dense density-density route expected by
`HFDMRG.solve_hfdmrg(H, V, psiup0, psidn0; ...)`. It deliberately avoids the
separate mixed-basis file-export question.
"""
function build_atomic_injected_angular_hfdmrg_hf_adapter(
    radial_ops::RadialAtomicOperators;
    shell_orders::Union{Nothing,AbstractVector{<:Integer}} = nothing,
    beta::Real = 2.0,
    l_inject::Union{Int,Symbol} = :auto,
    tau::Real = 1.0e-12,
    whiten::Symbol = :svd,
    nelec::Int = 2,
    maxiter::Int = 50,
    damping::Real = 0.25,
    tol::Real = 1.0e-8,
    ord_min::Int = minimum(curated_sphere_point_set_orders()),
    ord_max::Int = maximum(curated_sphere_point_set_orders()),
    r_lo::Real = 0.15,
    r_hi::Real = 7.0,
    w_lo::Real = 0.2,
    w_hi::Real = 0.7,
    nup::Union{Nothing,Int} = nothing,
    ndn::Union{Nothing,Int} = nothing,
    psiup0::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
    psidn0::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
)
    hf_style = build_atomic_injected_angular_hf_style_benchmark(
        radial_ops;
        shell_orders = shell_orders,
        beta = beta,
        l_inject = l_inject,
        tau = tau,
        whiten = whiten,
        nelec = nelec,
        maxiter = maxiter,
        damping = damping,
        tol = tol,
        ord_min = ord_min,
        ord_max = ord_max,
        r_lo = r_lo,
        r_hi = r_hi,
        w_lo = w_lo,
        w_hi = w_hi,
    )
    default_nocc = _default_angular_hfdmrg_spin_count(hf_style)
    resolved_nup = something(nup, psiup0 === nothing ? default_nocc : size(psiup0, 2))
    resolved_ndn = something(ndn, psidn0 === nothing ? resolved_nup : size(psidn0, 2))
    return build_atomic_injected_angular_hfdmrg_hf_adapter(
        hf_style;
        nup = resolved_nup,
        ndn = resolved_ndn,
        psiup0 = psiup0,
        psidn0 = psidn0,
    )
end

"""
    build_atomic_injected_angular_hfdmrg_hf_seeds(
        benchmark::AtomicInjectedAngularHFStyleBenchmark;
        nup=size(benchmark.scf_result.occupied_coefficients, 2),
        ndn=nup,
    )

Build the default open-shell-capable HFDMRG seed orbitals from the current
angular HF-style benchmark object.
"""
function build_atomic_injected_angular_hfdmrg_hf_seeds(
    benchmark::AtomicInjectedAngularHFStyleBenchmark;
    nup::Int = size(benchmark.scf_result.occupied_coefficients, 2),
    ndn::Int = nup,
)
    coefficients = benchmark.scf_result.coefficients
    norb = size(coefficients, 1)
    _validate_hfdmrg_spin_count(nup, norb, "nup")
    _validate_hfdmrg_spin_count(ndn, norb, "ndn")
    psiup0 = nup == 0 ? zeros(Float64, norb, 0) : _orthonormalize_columns(coefficients[:, 1:nup])
    psidn0 = ndn == 0 ? zeros(Float64, norb, 0) : _orthonormalize_columns(coefficients[:, 1:ndn])
    return (
        psiup0 = psiup0,
        psidn0 = psidn0,
        psiup0_source = :default_benchmark_orbitals,
        psidn0_source = :default_benchmark_orbitals,
    )
end

"""
    build_atomic_injected_angular_hfdmrg_hf_adapter(
        benchmark::AtomicInjectedAngularHFStyleBenchmark;
        nup=size(benchmark.scf_result.occupied_coefficients, 2),
        ndn=nup,
        psiup0=nothing,
        psidn0=nothing,
    )

Build the dense in-memory HFDMRG-facing HF adapter directly from the current
angular HF-style benchmark object, with explicit open-shell occupation and seed
control.
"""
function build_atomic_injected_angular_hfdmrg_hf_adapter(
    benchmark::AtomicInjectedAngularHFStyleBenchmark;
    nup::Int = size(benchmark.scf_result.occupied_coefficients, 2),
    ndn::Int = nup,
    psiup0::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
    psidn0::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
)
    norb = size(benchmark.one_body.hamiltonian, 1)
    _validate_hfdmrg_spin_count(nup, norb, "nup")
    _validate_hfdmrg_spin_count(ndn, norb, "ndn")

    default_seeds =
        psiup0 === nothing || psidn0 === nothing ?
        build_atomic_injected_angular_hfdmrg_hf_seeds(benchmark; nup = nup, ndn = ndn) :
        nothing

    overlap = benchmark.one_body.overlap
    overlap_identity_error =
        opnorm(overlap - Matrix{Float64}(I, size(overlap, 1), size(overlap, 2)), Inf)
    resolved_psiup0 =
        psiup0 === nothing ?
        default_seeds.psiup0 :
        _validated_hfdmrg_seed(psiup0, norb, nup, "psiup0")
    resolved_psidn0 =
        psidn0 === nothing ?
        default_seeds.psidn0 :
        _validated_hfdmrg_seed(psidn0, norb, ndn, "psidn0")
    psiup0_source = psiup0 === nothing ? default_seeds.psiup0_source : :explicit_seed
    psidn0_source = psidn0 === nothing ? default_seeds.psidn0_source : :explicit_seed

    return AtomicInjectedAngularHFDMRGHFAdapter(
        benchmark,
        :dense_density_density,
        Matrix{Float64}(benchmark.one_body.hamiltonian),
        Matrix{Float64}(benchmark.interaction),
        resolved_psiup0,
        resolved_psidn0,
        psiup0_source,
        psidn0_source,
        nup,
        ndn,
        overlap_identity_error,
    )
end

"""
    build_atomic_injected_angular_hfdmrg_hf_adapter(
        benchmark::AtomicInjectedAngularSmallEDBenchmark;
        kwargs...
    )

Build the same adapter starting from the small-ED benchmark wrapper.
"""
function build_atomic_injected_angular_hfdmrg_hf_adapter(
    benchmark::AtomicInjectedAngularSmallEDBenchmark;
    kwargs...,
)
    return build_atomic_injected_angular_hfdmrg_hf_adapter(benchmark.hf_style; kwargs...)
end

"""
    atomic_injected_angular_hfdmrg_hf_adapter_diagnostics(
        adapter::AtomicInjectedAngularHFDMRGHFAdapter
    )

Return compact diagnostics for the in-memory angular HFDMRG-facing HF adapter.
"""
function atomic_injected_angular_hfdmrg_hf_adapter_diagnostics(
    adapter::AtomicInjectedAngularHFDMRGHFAdapter,
)
    return (
        route = adapter.route,
        basis_dim = size(adapter.hamiltonian, 1),
        nup = adapter.nup,
        ndn = adapter.ndn,
        psiup0_source = adapter.psiup0_source,
        psidn0_source = adapter.psidn0_source,
        shell_orders = copy(adapter.hf_style.one_body.angular_assembly.shell_orders),
        exact_common_lmax = adapter.hf_style.one_body.exact_common_lmax,
        overlap_identity_error = adapter.overlap_identity_error,
        hamiltonian_symmetry_error =
            opnorm(adapter.hamiltonian - transpose(adapter.hamiltonian), Inf),
        interaction_symmetry_error =
            opnorm(adapter.interaction - transpose(adapter.interaction), Inf),
        psiup0_orthogonality_error =
            opnorm(transpose(adapter.psiup0) * adapter.psiup0 - Matrix{Float64}(I, adapter.nup, adapter.nup), Inf),
        psidn0_orthogonality_error =
            opnorm(transpose(adapter.psidn0) * adapter.psidn0 - Matrix{Float64}(I, adapter.ndn, adapter.ndn), Inf),
        benchmark_full_energy = adapter.hf_style.scf_result.energy,
        benchmark_exact_energy = adapter.hf_style.exact_scf_result.energy,
    )
end

"""
    run_atomic_injected_angular_hfdmrg_hf(
        adapter::AtomicInjectedAngularHFDMRGHFAdapter;
        hfmod=nothing,
        maxiter=40,
        blocksize=16,
        cutoff=1e-10,
        scf_cutoff=cutoff/10,
        verbose=false,
    )

Run the first in-memory HFDMRG-facing HF handshake on top of the angular
benchmark line.

This keeps the scope narrow: it delegates directly to
`HFDMRG.solve_hfdmrg(H, V, psiup0, psidn0; ...)` using the adapter's dense
density-density data. It is not a mixed-basis file-export solution and not a
true many-body DMRG adapter.
"""
function run_atomic_injected_angular_hfdmrg_hf(
    adapter::AtomicInjectedAngularHFDMRGHFAdapter;
    hfmod = nothing,
    maxiter::Int = 40,
    blocksize::Int = 16,
    cutoff::Real = 1.0e-10,
    scf_cutoff::Real = cutoff / 10,
    verbose::Bool = false,
)
    hf = _resolve_hfdmrg_module(hfmod)
    psiup, psidn, energy = hf.solve_hfdmrg(
        adapter.hamiltonian,
        adapter.interaction,
        adapter.psiup0,
        adapter.psidn0;
        maxiter = maxiter,
        blocksize = blocksize,
        cutoff = cutoff,
        scf_cutoff = scf_cutoff,
        verbose = verbose,
    )
    return (
        psiup = psiup,
        psidn = psidn,
        energy = energy,
        route = adapter.route,
        maxiter = maxiter,
        blocksize = blocksize,
        cutoff = cutoff,
        scf_cutoff = scf_cutoff,
    )
end

function run_atomic_injected_angular_hfdmrg_hf(
    benchmark::Union{
        AtomicInjectedAngularHFStyleBenchmark,
        AtomicInjectedAngularSmallEDBenchmark,
    };
    nup::Int = _default_angular_hfdmrg_spin_count(benchmark),
    ndn::Int = nup,
    psiup0::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
    psidn0::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
    kwargs...,
)
    adapter = build_atomic_injected_angular_hfdmrg_hf_adapter(
        benchmark;
        nup = nup,
        ndn = ndn,
        psiup0 = psiup0,
        psidn0 = psidn0,
    )
    return run_atomic_injected_angular_hfdmrg_hf(adapter; kwargs...)
end

"""
    build_atomic_injected_angular_small_ed_benchmark(radial_ops::RadialAtomicOperators;
                                                     shell_orders=nothing,
                                                     beta=2.0,
                                                     l_inject=:auto,
                                                     tau=1e-12,
                                                     whiten=:svd,
                                                     nelec::Int=2,
                                                     maxiter::Int=50,
                                                     damping::Real=0.25,
                                                     tol::Real=1e-8,
                                                     ord_min=minimum(curated_sphere_point_set_orders()),
                                                     ord_max=maximum(curated_sphere_point_set_orders()),
                                                     r_lo=0.15,
                                                     r_hi=7.0,
                                                     w_lo=0.2,
                                                     w_hi=0.7)

Build the first narrow interacting small-ED angular benchmark on top of the
existing one-electron and HF-style benchmark layer.
"""
function build_atomic_injected_angular_small_ed_benchmark(
    radial_ops::RadialAtomicOperators;
    shell_orders::Union{Nothing,AbstractVector{<:Integer}} = nothing,
    beta::Real = 2.0,
    l_inject::Union{Int,Symbol} = :auto,
    tau::Real = 1.0e-12,
    whiten::Symbol = :svd,
    nelec::Int = 2,
    maxiter::Int = 50,
    damping::Real = 0.25,
    tol::Real = 1.0e-8,
    ord_min::Int = minimum(curated_sphere_point_set_orders()),
    ord_max::Int = maximum(curated_sphere_point_set_orders()),
    r_lo::Real = 0.15,
    r_hi::Real = 7.0,
    w_lo::Real = 0.2,
    w_hi::Real = 0.7,
)
    hf_style = build_atomic_injected_angular_hf_style_benchmark(
        radial_ops;
        shell_orders = shell_orders,
        beta = beta,
        l_inject = l_inject,
        tau = tau,
        whiten = whiten,
        nelec = nelec,
        maxiter = maxiter,
        damping = damping,
        tol = tol,
        ord_min = ord_min,
        ord_max = ord_max,
        r_lo = r_lo,
        r_hi = r_hi,
        w_lo = w_lo,
        w_hi = w_hi,
    )
    return build_atomic_injected_angular_small_ed_benchmark(radial_ops, hf_style)
end

function build_atomic_injected_angular_small_ed_benchmark(
    radial_ops::RadialAtomicOperators,
    hf_style::AtomicInjectedAngularHFStyleBenchmark,
)
    one_body = hf_style.one_body
    exact_perm = _shell_major_orbital_permutation(
        size(radial_ops.overlap, 1),
        hf_style.exact_ida_reference.one_body.channels,
    )
    exact_reference_problem = atomic_ida_two_electron_problem(
        hf_style.exact_ida_reference;
        orbital_indices = exact_perm,
    )
    initial_orbital = hf_style.scf_result.occupied_coefficients[:, 1]
    initial_product = vec(initial_orbital * transpose(initial_orbital))
    lanczos = _lanczos_ground_state_apply(
        (out, v) -> _apply_density_density_two_electron_hamiltonian!(
            out,
            v,
            one_body.hamiltonian,
            hf_style.interaction,
        ),
        size(one_body.overlap, 1)^2;
        krylovdim = 160,
        maxiter = 160,
        tol = 1.0e-8,
        v0 = initial_product,
    )

    norb = size(one_body.overlap, 1)
    return AtomicInjectedAngularSmallEDBenchmark(
        hf_style,
        norb,
        norb^2,
        Matrix{Float64}(one_body.overlap),
        Matrix{Float64}(one_body.hamiltonian),
        Matrix{Float64}(hf_style.interaction),
        exact_reference_problem,
        lanczos,
    )
end

"""
    atomic_injected_angular_small_ed_diagnostics(benchmark::AtomicInjectedAngularSmallEDBenchmark)

Return compact diagnostics for the first narrow interacting angular small-ED
benchmark.
"""
function atomic_injected_angular_small_ed_diagnostics(
    benchmark::AtomicInjectedAngularSmallEDBenchmark,
)
    exact_problem = benchmark.exact_reference_problem
    full_energy = benchmark.lanczos_result.value
    exact_energy = ground_state_energy(exact_problem)
    norb = benchmark.orbital_count
    state_interaction_diagonal = vec(benchmark.orbital_interaction)

    return (
        nshells = length(benchmark.hf_style.one_body.angular_assembly.shells),
        orbital_count = benchmark.orbital_count,
        state_count = benchmark.state_count,
        shell_orders = copy(benchmark.hf_style.one_body.angular_assembly.shell_orders),
        exact_common_lmax = benchmark.hf_style.one_body.exact_common_lmax,
        orbital_overlap_symmetry_error = opnorm(benchmark.orbital_overlap - transpose(benchmark.orbital_overlap), Inf),
        orbital_one_body_symmetry_error = opnorm(benchmark.orbital_one_body - transpose(benchmark.orbital_one_body), Inf),
        orbital_interaction_symmetry_error = opnorm(benchmark.orbital_interaction - transpose(benchmark.orbital_interaction), Inf),
        orbital_overlap_identity_error =
            opnorm(benchmark.orbital_overlap - Matrix{Float64}(I, benchmark.orbital_count, benchmark.orbital_count), Inf),
        min_orbital_overlap_eigenvalue = minimum(eigvals(Symmetric(benchmark.orbital_overlap))),
        state_overlap_identity_error_estimate =
            opnorm(benchmark.orbital_overlap - Matrix{Float64}(I, norb, norb), Inf) * (1 + opnorm(benchmark.orbital_overlap, Inf)),
        min_state_overlap_eigenvalue_estimate = minimum(eigvals(Symmetric(benchmark.orbital_overlap)))^2,
        state_interaction_diagonal_min = minimum(state_interaction_diagonal),
        state_interaction_diagonal_max = maximum(state_interaction_diagonal),
        full_energy = full_energy,
        full_residual = benchmark.lanczos_result.residual,
        full_iterations = benchmark.lanczos_result.iterations,
        full_converged = benchmark.lanczos_result.converged,
        exact_reference_energy = exact_energy,
        energy_difference_to_exact_reference = full_energy - exact_energy,
    )
end
