struct ScreenedHartreeCorrection
    delta_one_body::Matrix{Float64}
    energy_constant::Float64
    q0::Vector{Float64}
    P0::Matrix{Float64}
    J0_G::Matrix{Float64}
    f_app_direct::Matrix{Float64}
    Vq0::Vector{Float64}
    E0_G::Float64
    q0Vq0::Float64
    energy_accounting::String
    diagnostics::NamedTuple
    packet_summary::NamedTuple
end

function represented_reference_p0_q0(
    coefficients::AbstractMatrix{<:Real},
    occupations::AbstractVector{<:Real},
)
    C = Matrix{Float64}(coefficients)
    occ = Vector{Float64}(occupations)
    size(C, 2) == length(occ) || throw(
        DimensionMismatch(
            "represented reference has $(size(C, 2)) orbitals but $(length(occ)) occupations",
        ),
    )
    all(isfinite, C) ||
        throw(ArgumentError("represented reference coefficients must be finite"))
    all(isfinite, occ) ||
        throw(ArgumentError("represented reference occupations must be finite"))
    all(>=(0.0), occ) ||
        throw(ArgumentError("represented reference occupations must be nonnegative"))
    P = _sym(C * Diagonal(occ) * transpose(C))
    q = Vector{Float64}(diag(P))
    gram = transpose(C) * C
    return (;
        P0 = P,
        q0 = q,
        trace = tr(P),
        occupation_sum = sum(occ),
        trace_loss = sum(occ) - tr(P),
        occupied_orthogonality_error =
            size(gram, 1) == 0 ? 0.0 :
            norm(gram - Matrix{Float64}(I, size(gram, 1), size(gram, 2)), Inf),
        coefficient_count = size(C, 2),
        dimension = size(C, 1),
    )
end
function represented_additive_reference_p0_q0(coefficient_blocks, occupation_blocks)
    length(coefficient_blocks) == length(occupation_blocks) &&
        !isempty(coefficient_blocks) || throw(ArgumentError(
            "additive represented reference needs matching nonempty blocks"))
    references = [represented_reference_p0_q0(C, occ)
        for (C, occ) in zip(coefficient_blocks, occupation_blocks)]
    all(reference -> reference.dimension == first(references).dimension, references) ||
        throw(DimensionMismatch("additive represented reference dimensions differ"))
    combined = represented_reference_p0_q0(
        hcat(coefficient_blocks...), vcat(occupation_blocks...))
    cross_overlap_max = maximum((maximum(abs, transpose(Matrix{Float64}(coefficient_blocks[a])) *
        Matrix{Float64}(coefficient_blocks[b]))
        for a in eachindex(coefficient_blocks) for b in (a + 1):length(coefficient_blocks));
        init = 0.0)
    return merge(combined, (;
        trace_loss = sum(reference.trace_loss for reference in references),
        occupied_orthogonality_error =
            maximum(reference.occupied_orthogonality_error for reference in references),
        cross_overlap_max,
        block_traces = [reference.trace for reference in references]))
end

function _screened_hartree_matrix(name, matrix)
    M = Matrix{Float64}(matrix)
    size(M, 1) == size(M, 2) ||
        throw(DimensionMismatch("$(name) must be square"))
    all(isfinite, M) || throw(ArgumentError("$(name) must be finite"))
    return _sym(M), norm(M - transpose(M), Inf)
end

function _screened_hartree_trace_product(A, B)
    size(A) == size(B) || throw(DimensionMismatch("trace-product size mismatch"))
    return Float64(sum(Matrix{Float64}(A) .* Matrix{Float64}(B)))
end

function _screened_hartree_check(condition, diagnostic_only, message)
    condition && return nothing
    diagnostic_only && return nothing
    throw(ArgumentError(message))
end

function _screened_hartree_validate_reference(
    reference;
    representation_atol::Real,
    q0_atol::Real,
    diagnostic_only::Bool,
)
    _screened_hartree_check(
        abs(reference.trace_loss) <= representation_atol,
        diagnostic_only,
        "screened-Hartree reference trace loss $(reference.trace_loss) exceeds $(representation_atol)",
    )
    _screened_hartree_check(
        reference.occupied_orthogonality_error <= representation_atol,
        diagnostic_only,
        "screened-Hartree represented occupied orbitals are not orthonormal; error $(reference.occupied_orthogonality_error) exceeds $(representation_atol)",
    )
    _screened_hartree_check(
        all(isfinite, reference.q0),
        diagnostic_only,
        "screened-Hartree q0 contains non-finite entries",
    )
    _screened_hartree_check(
        minimum(reference.q0) >= -q0_atol,
        diagnostic_only,
        "screened-Hartree q0 has negative entries below tolerance $(q0_atol)",
    )
    return nothing
end

function _packet_self_energy_nohalf(packet)
    hasproperty(packet, :density_fit) ||
        throw(ArgumentError("packet does not expose density_fit self-energy data"))
    return Float64(packet.density_fit.row.fit_self_energy)
end

function _packet_summary(packet, reference, source)
    packet === nothing && return (;
        source,
        represented_dimension = reference.dimension,
        reference_orbital_count = reference.coefficient_count,
        represented_trace = reference.trace,
        represented_trace_loss = reference.trace_loss,
        represented_orthogonality_error = reference.occupied_orthogonality_error,
    )
    atom = hasproperty(packet, :atom) ? packet.atom :
        hasproperty(packet, :spec) ? packet.spec.atom : "unknown"
    electron_count = hasproperty(packet, :electron_count) ? packet.electron_count :
        hasproperty(packet, :spec) ? packet.spec.electron_count : NaN
    basis_name = hasproperty(packet, :basis_name) ? packet.basis_name :
        hasproperty(packet, :spec) ? packet.spec.basis_name : "unknown"
    lmax = hasproperty(packet, :lmax) ? packet.lmax :
        hasproperty(packet, :spec) ? packet.spec.lmax : -1
    fill = hasproperty(packet, :fill_shell_convention) ? packet.fill_shell_convention :
        hasproperty(packet, :spec) ? packet.spec.fill_shell_convention : "unknown"
    path = hasproperty(packet, :path) ? packet.path : ""
    density_row = packet.density_fit.row
    potential_row = packet.potential_fit.row
    return (;
        source,
        packet_path = path,
        atom,
        electron_count,
        basis_name,
        lmax,
        fill_shell_convention = fill,
        represented_dimension = reference.dimension,
        reference_orbital_count = reference.coefficient_count,
        represented_trace = reference.trace,
        represented_trace_loss = reference.trace_loss,
        represented_orthogonality_error = reference.occupied_orthogonality_error,
        density_fit_charge = density_row.charge,
        density_fit_charge_error = density_row.charge_error,
        density_fit_self_energy_nohalf = density_row.fit_self_energy,
        density_fit_exact_self_energy_nohalf = density_row.exact_self_energy,
        density_fit_self_energy_relative_error =
            density_row.self_energy_relative_error,
        potential_fit_term_count = length(packet.potential_fit.coefficients),
        potential_fit_radial_relmax = potential_row.relmax,
        potential_fit_tail_charge_error = potential_row.tail_charge_error,
    )
end

function _build_screened_hartree_correction(
    V_IDA,
    J0_G,
    E0_G::Real,
    reference,
    additive_reference;
    packet = nothing,
    J0_G_exact = nothing,
    source::Symbol = :explicit_same_basis_inputs,
    representation_atol::Real = 1.0e-8,
    q0_atol::Real = 1.0e-12,
    input_symmetry_atol::Real = 1.0e-10,
    anchor_atol::Real = 1.0e-8,
    diagnostic_only::Bool = false,
)
    packet === nothing || _require_atomic_reference_converged(
        packet, "screened-Hartree packet consumption")
    V, V_input_symmetry_error = _screened_hartree_matrix("V_IDA", V_IDA)
    J, J_input_symmetry_error = _screened_hartree_matrix("J0_G", J0_G)
    _screened_hartree_check(
        V_input_symmetry_error <= input_symmetry_atol,
        diagnostic_only,
        "screened-Hartree V_IDA input symmetry error $(V_input_symmetry_error) exceeds $(input_symmetry_atol)",
    )
    _screened_hartree_check(
        J_input_symmetry_error <= input_symmetry_atol,
        diagnostic_only,
        "screened-Hartree J0_G input symmetry error $(J_input_symmetry_error) exceeds $(input_symmetry_atol)",
    )
    size(V) == size(J) ||
        throw(DimensionMismatch("V_IDA and J0_G dimensions differ"))
    isfinite(E0_G) || throw(ArgumentError("E0_G must be finite"))
    reference.dimension == size(V, 1) ||
        throw(DimensionMismatch(
            "represented reference dimension $(reference.dimension) does not match matrix dimension $(size(V, 1))",
        ))
    _screened_hartree_validate_reference(
        reference;
        representation_atol,
        q0_atol,
        diagnostic_only,
    )
    q0 = reference.q0
    Vq0 = Vector{Float64}(V * q0)
    f_app_direct = Matrix{Float64}(Diagonal(Vq0))
    delta_raw = Matrix{Float64}(J0_G) - f_app_direct
    delta_input_symmetry_error = norm(delta_raw - transpose(delta_raw), Inf)
    delta = _sym(J - f_app_direct)
    q0Vq0 = Float64(dot(q0, Vq0))
    constant = Float64(0.5 * q0Vq0 - 0.5 * Float64(E0_G))
    current_direct = Float64(0.5 * q0Vq0)
    delta_expectation = _screened_hartree_trace_product(reference.P0, delta)
    corrected_direct = current_direct + delta_expectation + constant
    exact_direct = Float64(0.5 * Float64(E0_G))
    derivative_error = norm(f_app_direct + delta - J, Inf)
    exact_diff = nothing
    exact_relative_fro = NaN
    exact_max_abs = NaN
    if J0_G_exact !== nothing
        exact, exact_input_symmetry_error = _screened_hartree_matrix(
            "J0_G_exact", J0_G_exact)
        _screened_hartree_check(
            exact_input_symmetry_error <= input_symmetry_atol,
            diagnostic_only,
            "screened-Hartree J0_G_exact input symmetry error $(exact_input_symmetry_error) exceeds $(input_symmetry_atol)",
        )
        size(exact) == size(J) ||
            throw(DimensionMismatch("J0_G_exact dimension differs from J0_G"))
        exact_diff = _sym(J - exact)
        exact_relative_fro = norm(exact_diff) / max(norm(exact), eps(Float64))
        exact_max_abs = maximum(abs, exact_diff)
    end
    _screened_hartree_check(
        abs(corrected_direct - exact_direct) <= anchor_atol,
        diagnostic_only,
        "screened-Hartree direct energy anchor error $(corrected_direct - exact_direct) exceeds $(anchor_atol)",
    )
    _screened_hartree_check(
        derivative_error <= anchor_atol,
        diagnostic_only,
        "screened-Hartree derivative anchor error $(derivative_error) exceeds $(anchor_atol)",
    )
    diagnostics = merge(isnothing(additive_reference) ? (;) : (; additive_reference), (;
        accounting =
            :screened_direct_electron_electron_interaction_not_physical_h1,
        q0_charge = sum(q0),
        q0_min = minimum(q0),
        q0_max = maximum(q0),
        P0_trace = reference.trace,
        P0_trace_loss = reference.trace_loss,
        P0_occupation_sum = reference.occupation_sum,
        represented_orthogonality_error = reference.occupied_orthogonality_error,
        J0_G_finite = all(isfinite, J),
        J0_G_input_symmetry_error = J_input_symmetry_error,
        J0_G_symmetry_error = norm(J - transpose(J), Inf),
        Delta_J0_finite = all(isfinite, delta),
        Delta_J0_input_symmetry_error = delta_input_symmetry_error,
        Delta_J0_symmetry_error = norm(delta - transpose(delta), Inf),
        E0_G_nohalf = Float64(E0_G),
        E_exact_direct = exact_direct,
        E_current_direct = current_direct,
        E_delta_expectation = delta_expectation,
        energy_constant = constant,
        E_corrected_direct_at_P0 = corrected_direct,
        direct_hartree_energy_anchor_error = corrected_direct - exact_direct,
        derivative_anchor_error = derivative_error,
        potential_fit_exact_relative_fro = exact_relative_fro,
        potential_fit_exact_max_abs = exact_max_abs,
        V_IDA_input_symmetry_error = V_input_symmetry_error,
        V_IDA_symmetry_error = norm(V - transpose(V), Inf),
        V_IDA_finite = all(isfinite, V),
        diagnostic_only,
    ))
    return ScreenedHartreeCorrection(
        delta,
        constant,
        q0,
        reference.P0,
        J,
        f_app_direct,
        Vq0,
        Float64(E0_G),
        q0Vq0,
        "Delta_J0 + C is accounted as screened direct electron-electron/Hartree interaction, not as a physical kinetic/nuclear one-body change.",
        diagnostics,
        _packet_summary(packet, reference, source),
    )
end
function build_screened_hartree_correction(
    V_IDA, J0_G, E0_G::Real, reference_coefficients::AbstractMatrix{<:Real},
    occupations::AbstractVector{<:Real}; kwargs...)
    return _build_screened_hartree_correction(V_IDA, J0_G, E0_G,
        represented_reference_p0_q0(reference_coefficients, occupations), nothing; kwargs...)
end
function build_additive_screened_hartree_correction(
    V_IDA, J0_G, E0_G::Real, coefficient_blocks, occupation_blocks;
    packets = nothing, kwargs...)
    if !isnothing(packets)
        length(packets) == length(coefficient_blocks) ||
            throw(DimensionMismatch("additive packet and coefficient block counts differ"))
    end
    foreach(packet -> _require_atomic_reference_converged(
        packet, "additive screened-Hartree packet consumption"), something(packets, ()))
    reference = represented_additive_reference_p0_q0(
        coefficient_blocks, occupation_blocks)
    additive_reference = (;
        block_count = length(coefficient_blocks),
        block_traces = reference.block_traces,
        interpacket_occupied_overlap_max = reference.cross_overlap_max)
    return _build_screened_hartree_correction(V_IDA, J0_G, E0_G,
        reference, additive_reference; source = :additive_atomic_packets, kwargs...)
end

function _screened_hartree_vee(ham_or_matrix)
    if hasproperty(ham_or_matrix, :electron_electron_ida)
        return ham_or_matrix.electron_electron_ida
    end
    return ham_or_matrix
end

function build_atomic_packet_screened_hartree_correction(
    base,
    ham_or_V_IDA,
    packet;
    reference_coefficients::AbstractMatrix{<:Real},
    occupations::AbstractVector{<:Real} =
        hasproperty(packet, :occupations) ? packet.occupations :
        hasproperty(packet, :C_occ) ? fill(2.0, size(packet.C_occ, 2)) :
        throw(ArgumentError("packet occupations are required")),
    source::Symbol = :potential_fit,
    check_density_fit_matrix::Bool = false,
    center = nothing,
    kwargs...,
)
    _require_atomic_reference_converged(
        packet, "screened-Hartree atomic packet consumption")
    source in (:potential_fit, :density_fit) ||
        throw(ArgumentError("source must be :potential_fit or :density_fit"))
    J0 = isnothing(center) ?
        atomic_reference_packet_terminal_hartree_gg(base, packet; source) :
        atomic_reference_packet_terminal_hartree_gg(base, packet; source, center)
    J0_exact =
        check_density_fit_matrix && source == :potential_fit ?
        isnothing(center) ?
        atomic_reference_packet_terminal_hartree_gg(base, packet; source = :density_fit) :
        atomic_reference_packet_terminal_hartree_gg(base, packet; source = :density_fit, center) :
        nothing
    return build_screened_hartree_correction(
        _screened_hartree_vee(ham_or_V_IDA),
        J0,
        _packet_self_energy_nohalf(packet),
        reference_coefficients,
        occupations;
        packet,
        J0_G_exact = J0_exact,
        source,
        kwargs...,
    )
end
