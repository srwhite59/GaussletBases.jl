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
