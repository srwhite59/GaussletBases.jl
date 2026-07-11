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

function _external_gto_import_from_validated_cross_overlap(S_FG::AbstractMatrix{<:Real},
    packet::ExternalGTOOrbitalPacket, validation; orthogonality_atol::Real = 1.0e-8)
    ordering_ok, overlap_ok, symmetry_error, expected_error = validation
    cross = Matrix{Float64}(S_FG)
    all(isfinite, cross) || throw(ArgumentError("external GTO final/GTO overlap is not finite"))
    size(cross, 2) == size(packet.S_GG, 1) || throw(DimensionMismatch(
        "cross-overlap column count does not match packet AO count"))
    alpha = _external_gto_spin_import(
        cross, packet.S_GG, packet.alpha; orthogonality_atol)
    beta = isnothing(packet.beta) ? nothing : _external_gto_spin_import(
        cross, packet.S_GG, packet.beta; orthogonality_atol)
    return ExternalGTOOrbitalImportResult(
        size(cross), true, ordering_ok, overlap_ok, symmetry_error,
        expected_error, alpha, beta)
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
    validation = _external_gto_validate_packet!(packet; validate_fingerprints,
        S_GG_atol, S_GG_rtol, S_GG_symmetry_atol)
    S_FG = Matrix{Float64}(gto_overlap_matrix(working, packet.probes))
    return _external_gto_import_from_validated_cross_overlap(
        S_FG, packet, validation; orthogonality_atol)
end
const _PROTECTED_XGTO_SIDECAR_KIND = :protected_localized_external_gto_representation
const _PROTECTED_XGTO_SIDECAR_CONVENTION = :protected_localized_external_gto_native_v1
struct ProtectedLocalizedExternalGTORepresentation{S,C,I}
    cross_overlap::Matrix{Float64}
    imported::ExternalGTOOrbitalImportResult
    source_metric::S
    occupied_capture::C
    identity::I
end
_external_gto_matrix_fingerprint(A) = external_gto_overlap_fingerprint(A)
_external_gto_recipe_fingerprint(basis, geometry) = bytes2hex(sha256(codeunits(repr((; basis, geometry)))))
_external_gto_prop(x, name::Symbol, default) = hasproperty(x, name) ? getproperty(x, name) : default
_external_gto_require(ok, message) = ok || throw(ArgumentError(message))
_external_gto_validate_packet_strict!(packet) = _external_gto_validate_packet!(packet; validate_fingerprints = true, S_GG_atol = 1.0e-10, S_GG_rtol = 1.0e-10, S_GG_symmetry_atol = 1.0e-12)
function _external_gto_capture_diagnostics(matrix; atol::Real = 1.0e-8)
    gram = transpose(matrix) * matrix
    raw = eigvals(Symmetric(0.5 .* (gram .+ transpose(gram))))
    all(isfinite, raw) || throw(ArgumentError("external GTO capture is not finite"))
    _external_gto_require(isempty(raw) || (minimum(raw) >= -atol && maximum(raw) <= 1 + atol), "external GTO capture lies outside [-$(atol), 1+$(atol)]")
    reported = clamp.(raw, 0.0, 1.0)
    return (; projected_gram_eigenvalues = reported,
        principal_singular_values = svdvals(matrix))
end
function _external_gto_source_metric_capture(S_LG, S_GG)
    metric = eigen(Symmetric(0.5 .* (S_GG .+ transpose(S_GG))))
    scale = isempty(metric.values) ? 0.0 : maximum(metric.values)
    tau = max(1.0e-12, 1.0e-10 * scale)
    _external_gto_require(isempty(metric.values) || minimum(metric.values) >= -tau, "external GTO source metric has a materially negative eigenvalue")
    keep = findall(>(tau), metric.values)
    Q = metric.vectors[:, keep] * Diagonal(inv.(sqrt.(metric.values[keep])))
    T = S_LG * Q
    capture = _external_gto_capture_diagnostics(T)
    limits = isempty(metric.values) ? (0.0, 0.0) : (minimum(metric.values), maximum(metric.values))
    return merge((; retained_rank = length(keep),
        discarded_count = length(metric.values) - length(keep), tau,
        metric_eigenvalue_min = limits[1], metric_eigenvalue_max = limits[2]), capture)
end
_external_gto_spin_capture(imported::ExternalGTOOrbitalSpinImport) = _external_gto_capture_diagnostics(imported.imported_coefficients)
function _external_gto_member_identity(member; member_artifact = nothing)
    recipe = member.recipe
    basis = (; nesting = _external_gto_prop(recipe, :nesting, :pqs), ns = Int(recipe.ns),
        core_spacing = Float64(recipe.core_spacing),
        s_factor = Float64(_external_gto_prop(recipe, :s_factor, 1.0)),
        basisname = String(recipe.basisname), lmax = Int(recipe.lmax))
    locations = reduce(vcat, [reshape(collect(Float64.(x)), 1, :) for x in recipe.atom_locations])
    geometry = (; atom_symbols = String.(recipe.atom_symbols), nuclear_charges = Float64.(recipe.nuclear_charges),
        atom_locations = locations, nup = Int(recipe.nup), ndn = Int(recipe.ndn))
    expansion = hasproperty(member, :coulomb_expansion) ? _cartesian_validate_coulomb_expansion_summary(member.coulomb_expansion) :
        _cartesian_coulomb_expansion_summary(member.inputs.base.input.coulomb_accuracy, member.inputs.base.coulomb_expansion)
    return (; final_dimension = size(member.H, 1),
        artifact_kind = _PROTECTED_LOCALIZED_ARTIFACT_KIND,
        convention_id = _PROTECTED_LOCALIZED_CONVENTION_ID,
        recipe_fingerprint = _external_gto_recipe_fingerprint(basis, geometry),
        H1_L_fingerprint = _external_gto_matrix_fingerprint(member.H),
        Vee_L_fingerprint = _external_gto_matrix_fingerprint(member.V),
        source_artifact = String(_external_gto_prop(recipe, :source_artifact, "")),
        member_artifact = isnothing(member_artifact) ? nothing : String(member_artifact),
        source_commit = String(_external_gto_prop(recipe, :source_commit, "unknown")),
        current_commit = _plb_current_commit(), basis, geometry, expansion)
end
function _external_gto_sidecar_identity(member, packet, S_LG; member_artifact = nothing)
    beta_fingerprint = isnothing(packet.beta) ? nothing : _external_gto_matrix_fingerprint(packet.beta.coefficients)
    return (; artifact_kind = _PROTECTED_XGTO_SIDECAR_KIND,
        format_version = 1, convention_id = _PROTECTED_XGTO_SIDECAR_CONVENTION,
        convention_version = 1, site_order_kind = :native,
        orientation = :final_by_external,
        cross_overlap_fingerprint = _external_gto_matrix_fingerprint(S_LG),
        external = (; ao_count = length(packet.ao_labels),
            ao_labels = copy(packet.ao_labels),
            ordering_fingerprint = packet.ordering_fingerprint,
            S_GG_fingerprint = packet.S_GG_fingerprint,
            alpha_coefficients_fingerprint =
                _external_gto_matrix_fingerprint(packet.alpha.coefficients),
            beta_coefficients_fingerprint = beta_fingerprint,
            alpha_occupations = copy(packet.alpha.occupations),
            beta_occupations = isnothing(packet.beta) ? nothing : copy(packet.beta.occupations),
            provenance = repr(packet.provenance)),
        protected = _external_gto_member_identity(member; member_artifact))
end
function protected_localized_external_gto_import(member,
    packet::ExternalGTOOrbitalPacket; member_artifact = nothing)
    validation = _external_gto_validate_packet_strict!(packet)
    raw_to_L = [member.raw.G_L; member.raw.A_L]
    handoff = _cartesian_final_gto_cross_overlap_handoff(
        member.inputs.base, member.inputs.supplement, raw_to_L, packet.probes;
        provenance = :protected_localized_external_gto_import)
    _external_gto_require(handoff.diagnostics.orientation === :final_by_gto, "protected external GTO handoff orientation must be :final_by_gto")
    S_LG = Matrix{Float64}(handoff.cross_overlap)
    size(S_LG, 1) == size(member.H, 1) || throw(
        DimensionMismatch("protected external GTO final dimension mismatch"))
    imported = _external_gto_import_from_validated_cross_overlap(S_LG, packet, validation)
    source_metric = _external_gto_source_metric_capture(S_LG, packet.S_GG)
    occupied = (; alpha = _external_gto_spin_capture(imported.alpha),
        beta = isnothing(imported.beta) ? nothing : _external_gto_spin_capture(imported.beta))
    identity = _external_gto_sidecar_identity(member, packet, S_LG; member_artifact)
    value = ProtectedLocalizedExternalGTORepresentation(S_LG, imported, source_metric, occupied, identity)
    !isnothing(member_artifact) && _external_gto_validate_saved_artifact(value, member_artifact)
    return value
end
function _external_gto_write_spin!(file, name, imported, capture, occupations)
    prefix = String(name)
    file["imported/$(prefix)/coefficients"] = imported.imported_coefficients
    file["imported/$(prefix)/occupations"] = Float64.(occupations)
    diagnostics = merge((; spin = imported.spin,
        source_orthogonality_error = imported.source_orthogonality_error,
        capture_matrix = imported.capture_matrix,
        orbital_captures = imported.orbital_captures,
        density_trace_source = imported.density_trace_source,
        density_trace_capture = imported.density_trace_capture,
        density_trace_loss = imported.density_trace_loss,
        worst_orbital_capture = imported.worst_orbital_capture), capture)
    _protected_localized_write_simple_group(file, "diagnostics/$(prefix)", diagnostics)
end
function write_protected_localized_external_gto_representation(path, value::ProtectedLocalizedExternalGTORepresentation)
    id, S = value.identity, value.cross_overlap
    _external_gto_require(_external_gto_matrix_fingerprint(S) == id.cross_overlap_fingerprint, "protected external GTO S_LG changed after construction")
    jldopen(String(path), "w") do file
        for name in (:artifact_kind, :format_version, :convention_id,
                :convention_version, :site_order_kind, :orientation)
            file[String(name)] = getproperty(id, name)
        end
        _protected_localized_write_simple_group(file, "cross_overlap", (;
            S_LG = S, fingerprint_sha256 = id.cross_overlap_fingerprint,
            final_dimension = size(S, 1), external_dimension = size(S, 2)))
        ext = id.external
        _protected_localized_write_simple_group(file, "external", (;
            ao_count = ext.ao_count, ao_labels = ext.ao_labels,
            ordering_fingerprint_sha256 = ext.ordering_fingerprint, S_GG_fingerprint_sha256 = ext.S_GG_fingerprint,
            alpha_coefficients_fingerprint_sha256 = ext.alpha_coefficients_fingerprint))
        file["external/provenance/repr"] = ext.provenance
        _external_gto_write_spin!(file, :alpha, value.imported.alpha,
            value.occupied_capture.alpha, ext.alpha_occupations)
        if !isnothing(value.imported.beta)
            file["external/beta_coefficients_fingerprint_sha256"] = ext.beta_coefficients_fingerprint
            _external_gto_write_spin!(file, :beta, value.imported.beta,
                value.occupied_capture.beta, ext.beta_occupations)
        end
        protected = id.protected
        _protected_localized_write_simple_group(file, "protected", (;
            final_dimension = protected.final_dimension, artifact_kind = protected.artifact_kind,
            convention_id = protected.convention_id, recipe_fingerprint_sha256 = protected.recipe_fingerprint,
            H1_L_fingerprint_sha256 = protected.H1_L_fingerprint, Vee_L_fingerprint_sha256 = protected.Vee_L_fingerprint,
            source_artifact = protected.source_artifact, source_commit = protected.source_commit,
            current_commit = protected.current_commit))
        !isnothing(protected.member_artifact) && (file["protected/member_artifact"] = protected.member_artifact)
        _protected_localized_write_simple_group(file, "protected/basis_controls", protected.basis)
        _protected_localized_write_simple_group(file, "protected/geometry_inputs", protected.geometry)
        _cartesian_write_coulomb_expansion_summary!(file, protected.expansion)
        _protected_localized_write_simple_group(file, "diagnostics/cross_overlap", (;
            final_dimension = size(S, 1), external_dimension = size(S, 2),
            finite = all(isfinite, S), max_abs = isempty(S) ? 0.0 : maximum(abs, S),
            fingerprint_sha256 = _external_gto_matrix_fingerprint(S)))
        _protected_localized_write_simple_group(file, "diagnostics/source_metric", value.source_metric)
        file["diagnostics/source_metric/S_GG_symmetry_error"] = value.imported.S_GG_symmetry_error
        file["diagnostics/source_metric/S_GG_expected_error"] = value.imported.S_GG_expected_error
    end
    return path
end
_external_gto_sidecar_key(file, key) = haskey(file, key) ? file[key] : throw(ArgumentError("missing protected external GTO sidecar key $(key)"))
function _external_gto_read_spin(file, name, final_dimension)
    prefix = String(name)
    diagkey(name) = "diagnostics/$(prefix)/$(name)"
    C = Matrix{Float64}(_external_gto_sidecar_key(file, "imported/$(prefix)/coefficients"))
    occupations = Float64.(_external_gto_sidecar_key(file, "imported/$(prefix)/occupations"))
    size(C, 1) == final_dimension || throw(DimensionMismatch("stored $(prefix) coefficient final dimension mismatch"))
    size(C, 2) == length(occupations) || throw(DimensionMismatch("stored $(prefix) coefficient/occupation count mismatch"))
    all(isfinite, C) && all(isfinite, occupations) || throw(ArgumentError("stored $(prefix) import is not finite"))
    capture_matrix = Matrix{Float64}(transpose(C) * C)
    stored_capture = Matrix{Float64}(_external_gto_sidecar_key(file, diagkey(:capture_matrix)))
    norm(capture_matrix - stored_capture, Inf) <= 1.0e-10 || throw(ArgumentError("stored $(prefix) capture matrix mismatch"))
    capture = _external_gto_capture_diagnostics(C)
    for field in (:principal_singular_values, :projected_gram_eigenvalues)
        stored = Float64.(_external_gto_sidecar_key(file, diagkey(field)))
        expected_length = field === :principal_singular_values ? min(size(C)...) : size(C, 2)
        length(stored) == expected_length || throw(DimensionMismatch("stored $(prefix) $(field) length mismatch"))
        norm(stored - getproperty(capture, field), Inf) <= 1.0e-10 || throw(ArgumentError("stored $(prefix) $(field) mismatch"))
    end
    orbital_captures = Vector{Float64}(diag(capture_matrix))
    stored_orbital = Float64.(_external_gto_sidecar_key(file, diagkey(:orbital_captures)))
    norm(stored_orbital - orbital_captures, Inf) <= 1.0e-10 || throw(ArgumentError("stored $(prefix) orbital captures mismatch"))
    source_trace = sum(occupations)
    captured_trace = dot(occupations, orbital_captures)
    source_error = Float64(_external_gto_sidecar_key(file, diagkey(:source_orthogonality_error)))
    _external_gto_require(isfinite(source_error) && 0 <= source_error <= 1.0e-8, "stored $(prefix) source orthogonality error is invalid")
    spin = Symbol(_external_gto_sidecar_key(file, diagkey(:spin)))
    imported = ExternalGTOOrbitalSpinImport(spin, C, source_error, capture_matrix,
        orbital_captures, source_trace, captured_trace, source_trace - captured_trace,
        isempty(orbital_captures) ? 1.0 : minimum(orbital_captures))
    for field in (:density_trace_source, :density_trace_capture, :density_trace_loss,
            :worst_orbital_capture)
        abs(Float64(_external_gto_sidecar_key(file, diagkey(field))) - getproperty(imported, field)) <= 1.0e-10 || throw(ArgumentError("stored $(prefix) $(field) mismatch"))
    end
    return (; imported, capture, occupations)
end
function _external_gto_read_source_metric(file, ao_count, final_dimension)
    prefix = "diagnostics/source_metric"
    rank = Int(_external_gto_sidecar_key(file, "$(prefix)/retained_rank"))
    discarded = Int(_external_gto_sidecar_key(file, "$(prefix)/discarded_count"))
    rank >= 0 && discarded >= 0 && rank + discarded == ao_count || throw(ArgumentError("stored source-metric rank accounting mismatch"))
    reported = Float64.(_external_gto_sidecar_key(file, "$(prefix)/projected_gram_eigenvalues"))
    singular = Float64.(_external_gto_sidecar_key(file, "$(prefix)/principal_singular_values"))
    length(reported) == rank && length(singular) == min(final_dimension, rank) || throw(DimensionMismatch("stored source-metric diagnostic length mismatch"))
    (isempty(reported) || (minimum(reported) >= 0.0 && maximum(reported) <= 1.0)) || throw(ArgumentError("stored source-metric capture is outside its physical range"))
    all(isfinite, reported) && all(isfinite, singular) && all(>=(0.0), singular) &&
        issorted(reported) && issorted(singular; rev = true) ||
        throw(ArgumentError("stored source-metric diagnostics are invalid"))
    squared_singular = sort!(vcat(zeros(rank - length(singular)), singular .^ 2))
    norm(squared_singular - sort(reported), Inf) <= 1.0e-10 || throw(ArgumentError("stored source-metric SVD/Gram diagnostics disagree"))
    tau = Float64(_external_gto_sidecar_key(file, "$(prefix)/tau"))
    extrema = Float64[_external_gto_sidecar_key(file, "$(prefix)/metric_eigenvalue_min"), _external_gto_sidecar_key(file, "$(prefix)/metric_eigenvalue_max")]
    _external_gto_require(isfinite(tau) && tau >= 0 && all(isfinite, extrema) && extrema[1] >= -tau && extrema[2] >= extrema[1] && tau == max(1.0e-12, 1.0e-10 * extrema[2]), "stored source-metric scalar diagnostics are invalid")
    return (; retained_rank = rank, discarded_count = discarded,
        tau, metric_eigenvalue_min = extrema[1], metric_eigenvalue_max = extrema[2],
        projected_gram_eigenvalues = reported, principal_singular_values = singular)
end
function _external_gto_validate_saved_packet(value, packet; source_state::Bool)
    ext = value.identity.external
    length(packet.ao_labels) == ext.ao_count && packet.ao_labels == ext.ao_labels || throw(ArgumentError("external GTO sidecar AO labels/order mismatch"))
    packet.ordering_fingerprint == ext.ordering_fingerprint || throw(ArgumentError("external GTO sidecar ordering fingerprint mismatch"))
    packet.S_GG_fingerprint == ext.S_GG_fingerprint || throw(ArgumentError("external GTO sidecar S_GG fingerprint mismatch"))
    if source_state
        _external_gto_matrix_fingerprint(packet.alpha.coefficients) == ext.alpha_coefficients_fingerprint || throw(ArgumentError("external GTO sidecar alpha source-state mismatch"))
        packet.alpha.occupations == ext.alpha_occupations || throw(ArgumentError("external GTO sidecar alpha occupations mismatch"))
        isnothing(packet.beta) == isnothing(ext.beta_coefficients_fingerprint) || throw(ArgumentError("external GTO sidecar beta source-state mismatch"))
        !isnothing(packet.beta) &&
            _external_gto_matrix_fingerprint(packet.beta.coefficients) !=
                ext.beta_coefficients_fingerprint &&
            throw(ArgumentError("external GTO sidecar beta source-state mismatch"))
        !isnothing(packet.beta) && packet.beta.occupations != ext.beta_occupations && throw(ArgumentError("external GTO sidecar beta occupations mismatch"))
    end
    return nothing
end
function _external_gto_validate_saved_member(value, member)
    stored = value.identity.protected
    current = _external_gto_member_identity(member; member_artifact = stored.member_artifact)
    for name in (:final_dimension, :artifact_kind, :convention_id, :recipe_fingerprint,
            :H1_L_fingerprint, :Vee_L_fingerprint, :source_artifact, :source_commit,
            :basis, :geometry, :expansion)
        getproperty(current, name) == getproperty(stored, name) || throw(ArgumentError("protected external GTO sidecar member $(name) mismatch"))
    end
    return nothing
end
function _external_gto_validate_saved_artifact(value, path)
    stored = value.identity.protected
    artifact_path = String(path)
    artifact = read_protected_localized_ida_hamiltonian(artifact_path)
    s_factor = jldopen(artifact_path, "r") do file
        Float64(_external_gto_sidecar_key(file, "basis_controls/s_factor"))
    end
    basis = (; nesting = artifact.basis_controls.nesting, ns = artifact.basis_controls.ns,
        core_spacing = artifact.basis_controls.core_spacing, s_factor,
        basisname = artifact.basis_controls.basisname, lmax = artifact.basis_controls.lmax)
    geometry = (; atom_symbols = artifact.geometry_inputs.atom_symbols, nuclear_charges = artifact.nuclear_charges,
        atom_locations = artifact.nuclear_positions, nup = artifact.nup, ndn = artifact.ndn)
    current = (; final_dimension = artifact.final_dimension, artifact_kind = artifact.artifact_kind,
        convention_id = artifact.convention_id,
        recipe_fingerprint = _external_gto_recipe_fingerprint(basis, geometry),
        H1_L_fingerprint = _external_gto_matrix_fingerprint(artifact.H1_L),
        Vee_L_fingerprint = _external_gto_matrix_fingerprint(artifact.Vee_L),
        source_artifact = artifact.provenance.source_artifact, source_commit = artifact.provenance.source_commit,
        current_commit = artifact.provenance.current_commit,
        basis, geometry, expansion = artifact.coulomb_expansion)
    _external_gto_require(!isnothing(current.expansion), "protected member artifact requires Coulomb provenance")
    for name in propertynames(current)
        _external_gto_require(getproperty(current, name) == getproperty(stored, name), "protected external GTO sidecar artifact $(name) mismatch")
    end
    return nothing
end
function read_protected_localized_external_gto_representation(path;
    packet = nothing, member = nothing, protected_artifact = nothing)
    value = jldopen(String(path), "r") do file
        Symbol(_external_gto_sidecar_key(file, "artifact_kind")) === _PROTECTED_XGTO_SIDECAR_KIND || throw(ArgumentError("unrecognized protected external GTO sidecar artifact_kind"))
        Int(_external_gto_sidecar_key(file, "format_version")) == 1 || throw(ArgumentError("unsupported protected external GTO sidecar format version"))
        Symbol(_external_gto_sidecar_key(file, "convention_id")) === _PROTECTED_XGTO_SIDECAR_CONVENTION || throw(ArgumentError("unrecognized protected external GTO sidecar convention"))
        Int(_external_gto_sidecar_key(file, "convention_version")) == 1 || throw(ArgumentError("unsupported protected external GTO convention version"))
        Symbol(_external_gto_sidecar_key(file, "site_order_kind")) === :native || throw(ArgumentError("protected external GTO sidecar must use native order"))
        Symbol(_external_gto_sidecar_key(file, "orientation")) === :final_by_external || throw(ArgumentError("protected external GTO sidecar orientation must be :final_by_external"))
        final_dimension = Int(_external_gto_sidecar_key(file, "cross_overlap/final_dimension"))
        external_dimension = Int(_external_gto_sidecar_key(file, "cross_overlap/external_dimension"))
        S = Matrix{Float64}(_external_gto_sidecar_key(file, "cross_overlap/S_LG"))
        size(S) == (final_dimension, external_dimension) || throw(DimensionMismatch("protected external GTO cross-overlap size mismatch"))
        all(isfinite, S) || throw(ArgumentError("stored S_LG is not finite"))
        fingerprint = String(_external_gto_sidecar_key(file, "cross_overlap/fingerprint_sha256"))
        fingerprint == _external_gto_matrix_fingerprint(S) || throw(ArgumentError("stored S_LG fingerprint mismatch"))
        for (key, expected) in (("final_dimension", final_dimension),
                ("external_dimension", external_dimension), ("finite", true),
                ("max_abs", isempty(S) ? 0.0 : maximum(abs, S)),
                ("fingerprint_sha256", fingerprint))
            _external_gto_sidecar_key(file, "diagnostics/cross_overlap/$(key)") == expected || throw(ArgumentError("stored cross-overlap diagnostic $(key) mismatch"))
        end
        ao_count = Int(_external_gto_sidecar_key(file, "external/ao_count"))
        ao_count == external_dimension || throw(DimensionMismatch("protected external GTO AO/cross-overlap dimension mismatch"))
        alpha = _external_gto_read_spin(file, :alpha, final_dimension)
        beta_keys = String["external/beta_coefficients_fingerprint_sha256",
            "imported/beta/coefficients", "imported/beta/occupations"]
        append!(beta_keys, ["diagnostics/beta/$(name)" for name in
            (:spin, :source_orthogonality_error, :capture_matrix, :orbital_captures,
                :density_trace_source, :density_trace_capture, :density_trace_loss,
                :worst_orbital_capture, :projected_gram_eigenvalues,
                :principal_singular_values)])
        beta_present = map(key -> haskey(file, key), beta_keys)
        any(beta_present) && !all(beta_present) && throw(ArgumentError("protected external GTO beta groups must be all present or absent"))
        beta = all(beta_present) ? _external_gto_read_spin(file, :beta, final_dimension) : nothing
        (isnothing(beta) ? alpha.imported.spin in (:restricted, :alpha) : alpha.imported.spin === :alpha && beta.imported.spin === :beta) || throw(ArgumentError("stored external GTO spin roles are inconsistent"))
        source_metric = _external_gto_read_source_metric(file, ao_count, final_dimension)
        metric_errors = Float64[_external_gto_sidecar_key(file, "diagnostics/source_metric/S_GG_symmetry_error"),
            _external_gto_sidecar_key(file, "diagnostics/source_metric/S_GG_expected_error")]
        _external_gto_require(all(isfinite, metric_errors) && all(>=(0.0), metric_errors) && metric_errors[1] <= 1.0e-12 && metric_errors[2] <= 1.0e-8, "stored source-metric overlap errors are invalid")
        imported = ExternalGTOOrbitalImportResult(size(S), true, true, true,
            metric_errors[1], metric_errors[2],
            alpha.imported, isnothing(beta) ? nothing : beta.imported)
        basis = (; nesting = Symbol(file["protected/basis_controls/nesting"]), ns = Int(file["protected/basis_controls/ns"]),
            core_spacing = Float64(file["protected/basis_controls/core_spacing"]), s_factor = Float64(file["protected/basis_controls/s_factor"]),
            basisname = String(file["protected/basis_controls/basisname"]), lmax = Int(file["protected/basis_controls/lmax"]))
        geometry = (; atom_symbols = String.(file["protected/geometry_inputs/atom_symbols"]),
            nuclear_charges = Float64.(file["protected/geometry_inputs/nuclear_charges"]), atom_locations = Matrix{Float64}(file["protected/geometry_inputs/atom_locations"]),
            nup = Int(file["protected/geometry_inputs/nup"]), ndn = Int(file["protected/geometry_inputs/ndn"]))
        member_artifact = haskey(file, "protected/member_artifact") ? String(file["protected/member_artifact"]) : nothing
        expansion = _cartesian_read_coulomb_expansion_summary(file)
        isnothing(expansion) && throw(ArgumentError("protected external GTO sidecar requires Coulomb provenance"))
        protected = (; final_dimension = Int(file["protected/final_dimension"]), artifact_kind = Symbol(file["protected/artifact_kind"]),
            convention_id = Symbol(file["protected/convention_id"]), recipe_fingerprint = String(file["protected/recipe_fingerprint_sha256"]),
            H1_L_fingerprint = String(file["protected/H1_L_fingerprint_sha256"]), Vee_L_fingerprint = String(file["protected/Vee_L_fingerprint_sha256"]),
            source_artifact = String(file["protected/source_artifact"]), member_artifact, source_commit = String(file["protected/source_commit"]),
            current_commit = String(file["protected/current_commit"]), basis, geometry, expansion)
        protected.artifact_kind === _PROTECTED_LOCALIZED_ARTIFACT_KIND && protected.convention_id === _PROTECTED_LOCALIZED_CONVENTION_ID || throw(ArgumentError("stored protected member identity is unrecognized"))
        external = (; ao_count, ao_labels = String.(file["external/ao_labels"]), ordering_fingerprint = String(file["external/ordering_fingerprint_sha256"]),
            S_GG_fingerprint = String(file["external/S_GG_fingerprint_sha256"]), alpha_coefficients_fingerprint = String(file["external/alpha_coefficients_fingerprint_sha256"]),
            beta_coefficients_fingerprint = isnothing(beta) ? nothing : String(file["external/beta_coefficients_fingerprint_sha256"]),
            alpha_occupations = alpha.occupations, beta_occupations = isnothing(beta) ? nothing : beta.occupations,
            provenance = String(file["external/provenance/repr"]))
        length(external.ao_labels) == ao_count || throw(DimensionMismatch("stored external AO label count mismatch"))
        protected.final_dimension == final_dimension || throw(DimensionMismatch("stored protected final dimension mismatch"))
        identity = (; artifact_kind = _PROTECTED_XGTO_SIDECAR_KIND, format_version = 1,
            convention_id = _PROTECTED_XGTO_SIDECAR_CONVENTION, convention_version = 1,
            site_order_kind = :native, orientation = :final_by_external,
            cross_overlap_fingerprint = fingerprint, external, protected)
        capture = (; alpha = alpha.capture, beta = isnothing(beta) ? nothing : beta.capture)
        ProtectedLocalizedExternalGTORepresentation(S, imported, source_metric, capture, identity)
    end
    if !isnothing(packet)
        _external_gto_validate_saved_packet(value, packet; source_state = true)
        validation = _external_gto_validate_packet_strict!(packet)
        fresh = _external_gto_import_from_validated_cross_overlap(
            value.cross_overlap, packet, validation)
        metric = _external_gto_source_metric_capture(value.cross_overlap, packet.S_GG)
        norm(metric.projected_gram_eigenvalues - value.source_metric.projected_gram_eigenvalues, Inf) <= 1.0e-10 || throw(ArgumentError("stored source-metric capture mismatch"))
        norm(metric.principal_singular_values - value.source_metric.principal_singular_values, Inf) <= 1.0e-10 || throw(ArgumentError("stored source-metric singular values mismatch"))
        norm(fresh.alpha.imported_coefficients - value.imported.alpha.imported_coefficients, Inf) <= 1.0e-10 || throw(ArgumentError("stored alpha import does not match S_LG*C_G"))
        abs(fresh.alpha.source_orthogonality_error - value.imported.alpha.source_orthogonality_error) <= 1.0e-10 || throw(ArgumentError("stored alpha source orthogonality mismatch"))
        if !isnothing(fresh.beta)
            norm(fresh.beta.imported_coefficients - value.imported.beta.imported_coefficients, Inf) <= 1.0e-10 || throw(ArgumentError("stored beta import does not match S_LG*C_G"))
            abs(fresh.beta.source_orthogonality_error - value.imported.beta.source_orthogonality_error) <= 1.0e-10 || throw(ArgumentError("stored beta source orthogonality mismatch"))
        end
    end
    !isnothing(member) && _external_gto_validate_saved_member(value, member)
    !isnothing(protected_artifact) && _external_gto_validate_saved_artifact(value, protected_artifact)
    return value
end
function external_gto_import_from_saved_representation(
    value::ProtectedLocalizedExternalGTORepresentation,
    packet::ExternalGTOOrbitalPacket)
    _external_gto_require(_external_gto_matrix_fingerprint(value.cross_overlap) == value.identity.cross_overlap_fingerprint, "protected external GTO S_LG fingerprint mismatch")
    _external_gto_validate_saved_packet(value, packet; source_state = false)
    validation = _external_gto_validate_packet_strict!(packet)
    return _external_gto_import_from_validated_cross_overlap(
        value.cross_overlap, packet, validation)
end
