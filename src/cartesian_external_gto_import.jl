"""
    ExternalGTOOrbitalSpinBlock(spin, coefficients, occupations)

Spin block from an explicit external GTO orbital packet. `coefficients` are
columns of molecular orbitals in the packet AO order, and `occupations` are the
electron occupations assigned to those columns.
"""
struct ExternalGTOOrbitalSpinBlock
    spin::Symbol
    coefficients::Matrix{Float64}
    occupations::Vector{Float64}
end

function ExternalGTOOrbitalSpinBlock(
    spin::Symbol,
    coefficients::AbstractMatrix{<:Real},
    occupations::AbstractVector{<:Real},
)
    C = Matrix{Float64}(coefficients)
    occ = Vector{Float64}(occupations)
    size(C, 2) == length(occ) || throw(
        DimensionMismatch(
            "external GTO spin block $(spin) has $(size(C, 2)) orbitals but $(length(occ)) occupations",
        ),
    )
    spin in (:restricted, :alpha, :beta) || throw(
        ArgumentError("external GTO spin block must be :restricted, :alpha, or :beta; got $(spin)"),
    )
    all(isfinite, C) || throw(ArgumentError("external GTO MO coefficients for $(spin) are not finite"))
    all(isfinite, occ) || throw(ArgumentError("external GTO occupations for $(spin) are not finite"))
    return ExternalGTOOrbitalSpinBlock(spin, C, occ)
end

ExternalGTOOrbitalSpinBlock(
    coefficients::AbstractMatrix{<:Real},
    occupations::AbstractVector{<:Real};
    spin::Symbol = :restricted,
) = ExternalGTOOrbitalSpinBlock(spin, coefficients, occupations)

"""
    ExternalGTOOrbitalPacket(...)

Resolved external GTO orbital packet. The packet stores the explicit probe basis
used for `gto_overlap_matrix`, packet-owned AO ordering metadata, the source
GTO overlap `S_GG`, spin orbital blocks, and provenance. `S_GG` is used only for
packet validation and source-orbital orthogonality checks.
"""
struct ExternalGTOOrbitalPacket{P, PR}
    probes::P
    centers::Vector{NTuple{3,Float64}}
    angular_powers::Vector{NTuple{3,Int}}
    exponents::Vector{Vector{Float64}}
    contraction_coefficients::Vector{Vector{Float64}}
    primitive_normalizations::Vector{Symbol}
    ao_labels::Vector{String}
    ordering_fingerprint::String
    S_GG::Matrix{Float64}
    S_GG_fingerprint::String
    alpha::ExternalGTOOrbitalSpinBlock
    beta::Union{Nothing,ExternalGTOOrbitalSpinBlock}
    provenance::PR
end

function ExternalGTOOrbitalPacket(
    probes::CartesianGaussianShellSupplementRepresentation3D,
    S_GG::AbstractMatrix{<:Real},
    alpha::ExternalGTOOrbitalSpinBlock;
    beta::Union{Nothing,ExternalGTOOrbitalSpinBlock} = nothing,
    ao_labels::AbstractVector{<:AbstractString} = [orbital.label for orbital in probes.orbitals],
    ordering_fingerprint::AbstractString = external_gto_ordering_fingerprint(
        probes; ao_labels = ao_labels),
    S_GG_fingerprint::AbstractString = external_gto_overlap_fingerprint(S_GG),
    provenance = (;),
)
    n = length(probes.orbitals)
    labels = String[String(label) for label in ao_labels]
    length(labels) == n || throw(
        DimensionMismatch("external GTO packet has $(length(labels)) AO labels for $(n) probe orbitals"),
    )
    S = Matrix{Float64}(S_GG)
    size(S) == (n, n) || throw(
        DimensionMismatch("external GTO S_GG has size $(size(S)); expected ($(n), $(n))"),
    )
    size(alpha.coefficients, 1) == n || throw(
        DimensionMismatch("alpha/restricted coefficient row count $(size(alpha.coefficients, 1)) does not match AO count $(n)"),
    )
    if beta !== nothing
        size(beta.coefficients, 1) == n || throw(
            DimensionMismatch("beta coefficient row count $(size(beta.coefficients, 1)) does not match AO count $(n)"),
        )
    end
    if beta === nothing
        alpha.spin in (:restricted, :alpha) || throw(
            ArgumentError("external GTO packet without beta block requires alpha/restricted spin block; got $(alpha.spin)"),
        )
    else
        alpha.spin == :alpha || throw(
            ArgumentError("spin-resolved external GTO packet requires alpha.spin == :alpha; got $(alpha.spin)"),
        )
        beta.spin == :beta || throw(
            ArgumentError("spin-resolved external GTO packet requires beta.spin == :beta; got $(beta.spin)"),
        )
    end
    all(isfinite, S) || throw(ArgumentError("external GTO S_GG is not finite"))
    centers = NTuple{3,Float64}[orbital.center for orbital in probes.orbitals]
    powers = NTuple{3,Int}[orbital.angular_powers for orbital in probes.orbitals]
    exponents = [Vector{Float64}(orbital.exponents) for orbital in probes.orbitals]
    coefficients = [Vector{Float64}(orbital.coefficients) for orbital in probes.orbitals]
    normalizations = Symbol[orbital.primitive_normalization for orbital in probes.orbitals]
    return ExternalGTOOrbitalPacket(
        probes,
        centers,
        powers,
        exponents,
        coefficients,
        normalizations,
        labels,
        String(ordering_fingerprint),
        S,
        String(S_GG_fingerprint),
        alpha,
        beta,
        provenance,
    )
end

"""
    ExternalGTOOrbitalSpinImport

Per-spin imported orbital coefficients and capture diagnostics.
"""
struct ExternalGTOOrbitalSpinImport
    spin::Symbol
    imported_coefficients::Matrix{Float64}
    source_orthogonality_error::Float64
    capture_matrix::Matrix{Float64}
    orbital_captures::Vector{Float64}
    density_trace_source::Float64
    density_trace_capture::Float64
    density_trace_loss::Float64
    worst_orbital_capture::Float64
end

"""
    ExternalGTOOrbitalImportResult

Result from `import_external_gto_orbitals`. Final-basis orbitals are imported
with `C_F = S_FG * C_G`; the result stores imported coefficients and compact
diagnostics, not Hamiltonian or interaction transforms.
"""
struct ExternalGTOOrbitalImportResult
    cross_overlap_size::Tuple{Int,Int}
    cross_overlap_finite::Bool
    ordering_fingerprint_valid::Bool
    S_GG_fingerprint_valid::Bool
    S_GG_symmetry_error::Float64
    S_GG_expected_error::Float64
    alpha::ExternalGTOOrbitalSpinImport
    beta::Union{Nothing,ExternalGTOOrbitalSpinImport}
end

function external_gto_overlap_fingerprint(S::AbstractMatrix{<:Real})
    data = Matrix{Float64}(S)
    return bytes2hex(sha256(reinterpret(UInt8, vec(data))))
end

function external_gto_ordering_fingerprint(
    probes::CartesianGaussianShellSupplementRepresentation3D;
    ao_labels::AbstractVector{<:AbstractString} = [orbital.label for orbital in probes.orbitals],
)
    signatures = [
        (
            label = String(ao_labels[index]),
            angular_powers = orbital.angular_powers,
            center = orbital.center,
            exponents = Tuple(Float64.(orbital.exponents)),
            coefficients = Tuple(Float64.(orbital.coefficients)),
            primitive_normalization = orbital.primitive_normalization,
        )
        for (index, orbital) in pairs(probes.orbitals)
    ]
    return bytes2hex(sha256(codeunits(repr(signatures))))
end

function _external_gto_validate_packet!(
    packet::ExternalGTOOrbitalPacket;
    validate_fingerprints::Bool,
    S_GG_atol::Real,
    S_GG_rtol::Real,
    S_GG_symmetry_atol::Real,
)
    expected_ordering = external_gto_ordering_fingerprint(
        packet.probes; ao_labels = packet.ao_labels)
    ordering_ok = packet.ordering_fingerprint == expected_ordering
    expected_overlap = external_gto_overlap_fingerprint(packet.S_GG)
    overlap_ok = packet.S_GG_fingerprint == expected_overlap
    if packet.beta === nothing
        packet.alpha.spin in (:restricted, :alpha) || throw(
            ArgumentError("external GTO packet without beta block requires alpha/restricted spin block; got $(packet.alpha.spin)"),
        )
    else
        packet.alpha.spin == :alpha || throw(
            ArgumentError("spin-resolved external GTO packet requires alpha.spin == :alpha; got $(packet.alpha.spin)"),
        )
        packet.beta.spin == :beta || throw(
            ArgumentError("spin-resolved external GTO packet requires beta.spin == :beta; got $(packet.beta.spin)"),
        )
    end
    if validate_fingerprints
        ordering_ok || throw(ArgumentError("external GTO packet AO ordering fingerprint mismatch"))
        overlap_ok || throw(ArgumentError("external GTO packet S_GG fingerprint mismatch"))
    end
    symmetry_error = norm(packet.S_GG - transpose(packet.S_GG), Inf)
    symmetry_error <= S_GG_symmetry_atol || throw(
        ArgumentError(
            "external GTO packet S_GG is not symmetric; error $(symmetry_error) exceeds $(S_GG_symmetry_atol)",
        ),
    )
    S_GG_expected = Matrix{Float64}(_cartesian_supplement_cross_overlap(packet.probes, packet.probes))
    expected_error = norm(packet.S_GG - S_GG_expected, Inf)
    expected_scale = max(norm(S_GG_expected, Inf), 1.0)
    expected_error <= S_GG_atol + S_GG_rtol * expected_scale || throw(
        ArgumentError(
            "external GTO packet S_GG does not match the overlap recomputed from packet.probes; error $(expected_error) exceeds $(S_GG_atol + S_GG_rtol * expected_scale)",
        ),
    )
    return ordering_ok, overlap_ok, Float64(symmetry_error), Float64(expected_error)
end

function _external_gto_spin_import(
    S_FG::AbstractMatrix{<:Real},
    S_GG::AbstractMatrix{<:Real},
    block::ExternalGTOOrbitalSpinBlock;
    orthogonality_atol::Real,
)
    C_G = block.coefficients
    C_F = Matrix{Float64}(S_FG * C_G)
    source_gram = transpose(C_G) * S_GG * C_G
    source_error =
        isempty(source_gram) ? 0.0 :
        norm(source_gram - I(size(source_gram, 1)), Inf)
    source_error <= orthogonality_atol || throw(
        ArgumentError(
            "external GTO $(block.spin) source orbitals are not orthonormal in S_GG; error $(source_error) exceeds $(orthogonality_atol)",
        ),
    )
    capture = Matrix{Float64}(transpose(C_F) * C_F)
    captures = isempty(capture) ? Float64[] : Vector{Float64}(diag(capture))
    density_source = sum(block.occupations)
    density_capture = isempty(captures) ? 0.0 : dot(block.occupations, captures)
    worst_capture = isempty(captures) ? 1.0 : minimum(captures)
    return ExternalGTOOrbitalSpinImport(
        block.spin,
        C_F,
        Float64(source_error),
        capture,
        captures,
        Float64(density_source),
        Float64(density_capture),
        Float64(density_source - density_capture),
        Float64(worst_capture),
    )
end

"""
    import_external_gto_orbitals(working, packet; validate_fingerprints = true,
        orthogonality_atol = 1e-8)

Import external GTO molecular orbitals into an orthonormal Cartesian final basis
with `C_F = S_FG * C_G`, where `S_FG = gto_overlap_matrix(working,
packet.probes)`. The source overlap `S_GG` is validation-only; this helper does
not build generalized final-basis overlap, Hamiltonian transforms, or
interaction transforms.
"""
function import_external_gto_orbitals(
    working,
    packet::ExternalGTOOrbitalPacket;
    validate_fingerprints::Bool = true,
    orthogonality_atol::Real = 1.0e-8,
    S_GG_atol::Real = 1.0e-10,
    S_GG_rtol::Real = 1.0e-10,
    S_GG_symmetry_atol::Real = 1.0e-12,
)
    ordering_ok, overlap_ok, symmetry_error, expected_error = _external_gto_validate_packet!(
        packet;
        validate_fingerprints = validate_fingerprints,
        S_GG_atol = S_GG_atol,
        S_GG_rtol = S_GG_rtol,
        S_GG_symmetry_atol = S_GG_symmetry_atol,
    )
    S_FG = Matrix{Float64}(gto_overlap_matrix(working, packet.probes))
    cross_finite = all(isfinite, S_FG)
    cross_finite || throw(ArgumentError("external GTO final/GTO overlap S_FG is not finite"))
    size(S_FG, 2) == size(packet.S_GG, 1) || throw(
        DimensionMismatch("S_FG column count $(size(S_FG, 2)) does not match S_GG size $(size(packet.S_GG, 1))"),
    )
    alpha = _external_gto_spin_import(
        S_FG,
        packet.S_GG,
        packet.alpha;
        orthogonality_atol = orthogonality_atol,
    )
    beta =
        packet.beta === nothing ? nothing :
        _external_gto_spin_import(
            S_FG,
            packet.S_GG,
            packet.beta;
            orthogonality_atol = orthogonality_atol,
        )
    return ExternalGTOOrbitalImportResult(
        size(S_FG),
        cross_finite,
        ordering_ok,
        overlap_ok,
        symmetry_error,
        expected_error,
        alpha,
        beta,
    )
end
