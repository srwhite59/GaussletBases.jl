const _XGTO_KIND = "gaussletbases.external_cartesian_gto"
const _XGTO_HASH_PATTERN = r"^[0-9a-f]{64}$"

_xgto_fail(message) = throw(ArgumentError("external Cartesian GTO bundle: " * message))

function _xgto_exact_keys(table, name, expected)
    table isa AbstractDict || _xgto_fail("$(name) must be a TOML table")
    actual = sort!(String[String(key) for key in keys(table)])
    wanted = sort!(String[String(key) for key in expected])
    actual == wanted || _xgto_fail("$(name) keys $(actual) do not match $(wanted)")
    return table
end

function _xgto_string(value, name)
    value isa String && !isempty(value) || _xgto_fail("$(name) must be a nonempty string")
    return value
end

function _xgto_int(value, name; minimum = nothing)
    value isa Integer && !(value isa Bool) || _xgto_fail("$(name) must be an integer")
    result = Int(value)
    isnothing(minimum) || result >= minimum || _xgto_fail("$(name) is below $(minimum)")
    return result
end

function _xgto_float(value, name)
    value isa Real && !(value isa Bool) || _xgto_fail("$(name) must be real")
    result = Float64(value)
    isfinite(result) || _xgto_fail("$(name) must be finite")
    return result
end

function _xgto_float_vector(value, name)
    value isa AbstractVector || _xgto_fail("$(name) must be an array")
    return Float64[_xgto_float(item, "$(name) entry") for item in value]
end

function _xgto_int_vector(value, name; minimum = nothing)
    value isa AbstractVector || _xgto_fail("$(name) must be an array")
    return Int[_xgto_int(item, "$(name) entry"; minimum) for item in value]
end

function _xgto_triplet(value, name; integer = false)
    values = integer ? _xgto_int_vector(value, name; minimum = 0) :
             _xgto_float_vector(value, name)
    length(values) == 3 || _xgto_fail("$(name) must contain three entries")
    return integer ? (values[1], values[2], values[3]) :
                     (values[1], values[2], values[3])
end

function _xgto_hash(value, name)
    text = _xgto_string(value, name)
    occursin(_XGTO_HASH_PATTERN, text) || _xgto_fail("$(name) must be lowercase SHA-256")
    return text
end

function _xgto_array_data(payload, records)
    arrays = Dict{String,Any}()
    expected_offset = 0
    for (index, raw) in pairs(records)
        record = _xgto_exact_keys(raw, "arrays[$(index)]",
            ("name", "shape", "byte_offset", "element_count", "sha256"))
        name = _xgto_string(record["name"], "arrays[$(index)].name")
        haskey(arrays, name) && _xgto_fail("duplicate payload array $(name)")
        shape = _xgto_int_vector(record["shape"], "arrays[$(index)].shape"; minimum = 1)
        length(shape) in (1, 2) || _xgto_fail("array $(name) must be one- or two-dimensional")
        count = _xgto_int(record["element_count"], "arrays[$(index)].element_count"; minimum = 1)
        count == prod(shape) || _xgto_fail("array $(name) shape/count mismatch")
        offset = _xgto_int(record["byte_offset"], "arrays[$(index)].byte_offset"; minimum = 0)
        offset == expected_offset || _xgto_fail("array $(name) is not contiguous in manifest order")
        final_offset = offset + 8 * count
        final_offset <= length(payload) || _xgto_fail("array $(name) exceeds payload")
        bytes = payload[(offset + 1):final_offset]
        bytes2hex(sha256(bytes)) == _xgto_hash(record["sha256"], "arrays[$(index)].sha256") ||
            _xgto_fail("array $(name) SHA-256 mismatch")
        words = reinterpret(UInt64, bytes)
        values = Float64[reinterpret(Float64, ltoh(word)) for word in words]
        arrays[name] = length(shape) == 1 ? values : reshape(values, shape[1], shape[2])
        expected_offset = final_offset
    end
    expected_offset == length(payload) || _xgto_fail("payload has trailing bytes")
    return arrays
end

function _xgto_nuclei(records)
    records isa AbstractVector || _xgto_fail("nuclei must be an array of tables")
    charges = Float64[]
    positions = zeros(Float64, length(records), 3)
    symbols = String[]
    for (index, raw) in pairs(records)
        row = _xgto_exact_keys(raw, "nuclei[$(index)]",
            ("atom_index_1based", "symbol", "nuclear_charge", "center_bohr"))
        _xgto_int(row["atom_index_1based"], "nuclei[$(index)].atom_index_1based") == index ||
            _xgto_fail("nuclear indices must be contiguous and ordered")
        push!(symbols, _xgto_string(row["symbol"], "nuclei[$(index)].symbol"))
        charge = _xgto_float(row["nuclear_charge"], "nuclei[$(index)].nuclear_charge")
        charge > 0 || _xgto_fail("nuclear charges must be positive")
        push!(charges, charge)
        positions[index, :] .= _xgto_triplet(row["center_bohr"], "nuclei[$(index)].center_bohr")
    end
    isempty(charges) && _xgto_fail("at least one nucleus is required")
    return symbols, charges, positions
end

function _xgto_raw_probes(records, symbols, nuclear_positions, ao_count)
    records isa AbstractVector || _xgto_fail("aos must be an array of tables")
    length(records) == ao_count || _xgto_fail("AO record count does not match state.ao_count")
    orbitals = CartesianGaussianShellOrbitalRepresentation3D[]
    group_keys = Tuple{Int,Int}[]
    for (index, raw) in pairs(records)
        row = _xgto_exact_keys(raw, "aos[$(index)]", ("ao_index_1based", "label",
            "atom_index_1based", "shell_index_1based", "contraction_index_1based",
            "angular_powers", "center_bohr", "exponents_bohr_inverse_square",
            "contraction_coefficients"))
        _xgto_int(row["ao_index_1based"], "aos[$(index)].ao_index_1based") == index ||
            _xgto_fail("AO indices must be contiguous and ordered")
        atom = _xgto_int(row["atom_index_1based"], "aos[$(index)].atom_index_1based"; minimum = 1)
        atom <= length(symbols) || _xgto_fail("AO $(index) atom index is out of range")
        shell = _xgto_int(row["shell_index_1based"], "aos[$(index)].shell_index_1based"; minimum = 1)
        contraction = _xgto_int(row["contraction_index_1based"], "aos[$(index)].contraction_index_1based"; minimum = 1)
        push!(group_keys, (shell, contraction))
        powers = _xgto_triplet(row["angular_powers"], "aos[$(index)].angular_powers"; integer = true)
        center = _xgto_triplet(row["center_bohr"], "aos[$(index)].center_bohr")
        norm(collect(center) - nuclear_positions[atom, :], Inf) <= 1.0e-10 ||
            _xgto_fail("AO $(index) center does not match its nucleus")
        exponents = _xgto_float_vector(row["exponents_bohr_inverse_square"], "aos[$(index)].exponents")
        coefficients = _xgto_float_vector(row["contraction_coefficients"], "aos[$(index)].coefficients")
        !isempty(exponents) && length(exponents) == length(coefficients) ||
            _xgto_fail("AO $(index) primitive dimensions are invalid")
        all(>(0), exponents) || _xgto_fail("AO $(index) exponents must be positive")
        push!(orbitals, CartesianGaussianShellOrbitalRepresentation3D(
            _xgto_string(row["label"], "aos[$(index)].label"), powers, center,
            exponents, coefficients, :axiswise_normalized_cartesian_gaussian))
    end
    issorted(group_keys) || _xgto_fail("AO shell/contraction groups are not ordered")
    for key in unique(group_keys)
        indices = findall(==(key), group_keys)
        first(indices):last(indices) == indices || _xgto_fail("AO component group $(key) is not contiguous")
        angular_momentum = sum(orbitals[first(indices)].angular_powers)
        expected = NTuple{3,Int}[(lx, ly, angular_momentum - lx - ly)
            for lx in angular_momentum:-1:0 for ly in (angular_momentum - lx):-1:0]
        getproperty.(orbitals[indices], :angular_powers) == expected ||
            _xgto_fail("AO component order mismatch in group $(key)")
    end
    metadata = (; source_kind = :external_cartesian_gto, format_version = 1,
        atom_symbols = symbols, nuclei = [Tuple(nuclear_positions[row, :]) for row in axes(nuclear_positions, 1)])
    return CartesianGaussianShellSupplementRepresentation3D(
        :external_cartesian_gto, orbitals, metadata)
end

function _xgto_orbital_blocks(records, arrays, S_GG, state_kind)
    records isa AbstractVector || _xgto_fail("orbital_blocks must be an array of tables")
    expected_spins = state_kind == "rhf" ? ["restricted"] :
        state_kind == "alpha_only" ? ["alpha"] : state_kind == "uhf" ? ["alpha", "beta"] :
        _xgto_fail("unsupported state.kind $(state_kind)")
    length(records) == length(expected_spins) || _xgto_fail("orbital block count mismatch")
    blocks = ExternalGTOOrbitalSpinBlock[]
    mo_indices = Vector{Vector{Int}}()
    for (index, expected_spin) in pairs(expected_spins)
        row = _xgto_exact_keys(records[index], "orbital_blocks[$(index)]",
            ("spin", "mo_indices_1based", "metric_error_inf"))
        spin_text = _xgto_string(row["spin"], "orbital_blocks[$(index)].spin")
        spin_text == expected_spin || _xgto_fail("orbital block spin/order mismatch")
        indices = _xgto_int_vector(row["mo_indices_1based"], "orbital_blocks[$(index)].mo_indices"; minimum = 1)
        !isempty(indices) && issorted(indices) && allunique(indices) ||
            _xgto_fail("orbital MO indices must be nonempty, sorted, and unique")
        prefix = expected_spin == "beta" ? "beta" : "alpha"
        coefficient_name = prefix * "_coefficients"
        occupation_name = prefix * "_occupations"
        haskey(arrays, coefficient_name) && haskey(arrays, occupation_name) ||
            _xgto_fail("missing arrays for $(expected_spin) block")
        coefficients = arrays[coefficient_name]
        occupations = vec(arrays[occupation_name])
        coefficients isa Matrix{Float64} || _xgto_fail("$(coefficient_name) must be a matrix")
        size(coefficients, 2) == length(indices) == length(occupations) ||
            _xgto_fail("$(expected_spin) orbital dimensions disagree")
        all(>(0), occupations) || _xgto_fail("exported occupations must be positive")
        metric_error = norm(transpose(coefficients) * S_GG * coefficients - I, Inf)
        metric_error <= 1.0e-8 || _xgto_fail("$(expected_spin) source metric error $(metric_error)")
        abs(metric_error - _xgto_float(row["metric_error_inf"], "orbital block metric error")) <= 1.0e-8 ||
            _xgto_fail("$(expected_spin) metric diagnostic mismatch")
        push!(blocks, ExternalGTOOrbitalSpinBlock(Symbol(expected_spin), coefficients, occupations))
        push!(mo_indices, indices)
    end
    expected_arrays = expected_spins == ["alpha", "beta"] ?
        ["source_overlap", "alpha_coefficients", "alpha_occupations", "beta_coefficients", "beta_occupations"] :
        ["source_overlap", "alpha_coefficients", "alpha_occupations"]
    sort!(collect(keys(arrays))) == sort!(expected_arrays) || _xgto_fail("unexpected payload array inventory")
    return blocks, mo_indices
end

"""
    read_external_cartesian_gto_packet(path; overlap_atol=1e-10, overlap_rtol=1e-10)

Read and validate a version-1 external Cartesian-GTO `.toml`/`.f64` bundle.
The returned packet preserves the exported AO order and raw molecular-orbital
coefficients. Its declared Hamiltonian kind is producer attestation, not a
claim derivable from the source checkpoint.
"""
function read_external_cartesian_gto_packet(path::AbstractString;
    overlap_atol::Real = 1.0e-10, overlap_rtol::Real = 1.0e-10)
    atol = _xgto_float(overlap_atol, "overlap_atol")
    rtol = _xgto_float(overlap_rtol, "overlap_rtol")
    atol >= 0 && rtol >= 0 || _xgto_fail("overlap tolerances must be nonnegative")
    manifest_path = abspath(String(path))
    splitext(manifest_path)[2] == ".toml" || _xgto_fail("manifest path must end in .toml")
    isfile(manifest_path) || _xgto_fail("manifest does not exist")
    document = _xgto_exact_keys(TOML.parsefile(manifest_path), "manifest",
        ("format", "encoding", "units", "producer", "state", "source", "payload",
         "nuclei", "aos", "orbital_blocks", "arrays"))
    format = _xgto_exact_keys(document["format"], "format", ("kind", "version"))
    _xgto_string(format["kind"], "format.kind") == _XGTO_KIND || _xgto_fail("unsupported format kind")
    _xgto_int(format["version"], "format.version") == 1 || _xgto_fail("unsupported format version")
    encoding = _xgto_exact_keys(document["encoding"], "encoding",
        ("scalar_type", "byte_order", "array_order"))
    (_xgto_string(encoding["scalar_type"], "encoding.scalar_type"),
     _xgto_string(encoding["byte_order"], "encoding.byte_order"),
     _xgto_string(encoding["array_order"], "encoding.array_order")) ==
        ("ieee754_float64", "little_endian", "column_major") || _xgto_fail("unsupported encoding")
    units = _xgto_exact_keys(document["units"], "units", ("length", "exponent"))
    (_xgto_string(units["length"], "units.length"), _xgto_string(units["exponent"], "units.exponent")) ==
        ("bohr", "bohr_inverse_square") || _xgto_fail("unsupported units")
    producer = _xgto_exact_keys(document["producer"], "producer",
        ("name", "version", "python_version", "numpy_version", "exporter", "exporter_version"))
    for key in keys(producer); _xgto_string(producer[key], "producer.$(key)"); end
    state = _xgto_exact_keys(document["state"], "state", ("kind", "source_ao_representation",
        "exported_ao_representation", "cartesian_normalization", "spherical_conversion",
        "hamiltonian_kind", "molecular_charge", "spin_2s", "electron_count", "nalpha", "nbeta", "ao_count"))
    state_kind = _xgto_string(state["kind"], "state.kind")
    source_kind = _xgto_string(state["source_ao_representation"], "state.source_ao_representation")
    source_kind in ("spherical", "cartesian") || _xgto_fail("unsupported source AO representation")
    _xgto_string(state["exported_ao_representation"], "state.exported_ao_representation") == "cartesian" ||
        _xgto_fail("exported AO representation must be cartesian")
    _xgto_string(state["cartesian_normalization"], "state.cartesian_normalization") == "pyscf_libcint_sp" ||
        _xgto_fail("unsupported Cartesian normalization")
    conversion = _xgto_string(state["spherical_conversion"], "state.spherical_conversion")
    conversion == (source_kind == "spherical" ? "mol.cart2sph_coeff(normalized=sp)" : "not_applied") ||
        _xgto_fail("source representation/conversion mismatch")
    hamiltonian_kind = _xgto_string(state["hamiltonian_kind"], "state.hamiltonian_kind")
    hamiltonian_kind == "all_electron_nuclear_coulomb" || _xgto_fail("unsupported Hamiltonian declaration")
    molecular_charge = _xgto_int(state["molecular_charge"], "state.molecular_charge")
    spin_2s = _xgto_int(state["spin_2s"], "state.spin_2s")
    electron_count = _xgto_int(state["electron_count"], "state.electron_count"; minimum = 1)
    nalpha = _xgto_int(state["nalpha"], "state.nalpha"; minimum = 0)
    nbeta = _xgto_int(state["nbeta"], "state.nbeta"; minimum = 0)
    ao_count = _xgto_int(state["ao_count"], "state.ao_count"; minimum = 1)
    electron_count == nalpha + nbeta && spin_2s == nalpha - nbeta || _xgto_fail("electron-sector fields disagree")
    symbols, nuclear_charges, nuclear_positions = _xgto_nuclei(document["nuclei"])
    abs(sum(nuclear_charges) - electron_count - molecular_charge) <= 1.0e-10 ||
        _xgto_fail("molecular charge disagrees with nuclei and electron count")
    probes_raw = _xgto_raw_probes(document["aos"], symbols, nuclear_positions, ao_count)
    source = _xgto_exact_keys(document["source"], "source", ("checkpoint_sha256",))
    checkpoint_sha256 = _xgto_hash(source["checkpoint_sha256"], "source.checkpoint_sha256")
    payload_record = _xgto_exact_keys(document["payload"], "payload", ("basename", "byte_count", "sha256"))
    payload_name = _xgto_string(payload_record["basename"], "payload.basename")
    basename(payload_name) == payload_name && splitext(payload_name)[2] == ".f64" &&
        splitext(payload_name)[1] == splitext(basename(manifest_path))[1] ||
        _xgto_fail("payload must be the same-stem .f64 basename")
    payload_path = joinpath(dirname(manifest_path), payload_name)
    isfile(payload_path) || _xgto_fail("payload does not exist")
    payload = read(payload_path)
    length(payload) == _xgto_int(payload_record["byte_count"], "payload.byte_count"; minimum = 1) ||
        _xgto_fail("payload byte count mismatch")
    payload_sha256 = _xgto_hash(payload_record["sha256"], "payload.sha256")
    bytes2hex(sha256(payload)) == payload_sha256 || _xgto_fail("payload SHA-256 mismatch")
    document["arrays"] isa AbstractVector || _xgto_fail("arrays must be an array of tables")
    arrays = _xgto_array_data(payload, document["arrays"])
    haskey(arrays, "source_overlap") || _xgto_fail("source_overlap array is missing")
    S_GG = arrays["source_overlap"]
    S_GG isa Matrix{Float64} && size(S_GG) == (ao_count, ao_count) ||
        _xgto_fail("source_overlap dimensions are invalid")
    all(isfinite, S_GG) && norm(S_GG - transpose(S_GG), Inf) <= 1.0e-12 ||
        _xgto_fail("source_overlap must be finite and symmetric")
    raw_overlap = _cartesian_supplement_cross_overlap(probes_raw, probes_raw)
    raw_diagonal = diag(raw_overlap)
    source_diagonal = diag(S_GG)
    all(isfinite, raw_diagonal) && all(>(0), raw_diagonal) &&
        all(isfinite, source_diagonal) && all(>(0), source_diagonal) ||
        _xgto_fail("source/raw overlap diagonals must be finite and positive")
    scales = sqrt.(source_diagonal ./ raw_diagonal)
    all(isfinite, scales) && all(>(0), scales) || _xgto_fail("AO normalization factors must be finite and positive")
    scaled_orbitals = CartesianGaussianShellOrbitalRepresentation3D[
        CartesianGaussianShellOrbitalRepresentation3D(orbital.label, orbital.angular_powers,
            orbital.center, copy(orbital.exponents), scales[index] .* orbital.coefficients,
            orbital.primitive_normalization) for (index, orbital) in pairs(probes_raw.orbitals)]
    probes = CartesianGaussianShellSupplementRepresentation3D(
        :external_cartesian_gto, scaled_orbitals, probes_raw.metadata)
    probe_overlap = _cartesian_supplement_cross_overlap(probes, probes)
    overlap_error = norm(probe_overlap - S_GG, Inf)
    overlap_error <= atol + rtol * max(norm(S_GG, Inf), 1.0) ||
        _xgto_fail("full source-overlap parity failed with error $(overlap_error)")
    blocks, mo_indices = _xgto_orbital_blocks(document["orbital_blocks"], arrays, S_GG, state_kind)
    if state_kind == "rhf"
        nalpha == nbeta && abs(sum(blocks[1].occupations) - electron_count) <= 1.0e-10 ||
            _xgto_fail("restricted occupations disagree with the electron sector")
    else
        (state_kind == "alpha_only" && nbeta == 0) || state_kind == "uhf" ||
            _xgto_fail("state kind/electron sector mismatch")
        abs(sum(blocks[1].occupations) - nalpha) <= 1.0e-10 ||
            _xgto_fail("alpha occupations disagree with the electron sector")
        beta_occupation = length(blocks) == 2 ? sum(blocks[2].occupations) : 0.0
        abs(beta_occupation - nbeta) <= 1.0e-10 ||
            _xgto_fail("beta occupations disagree with the electron sector")
    end
    identity = (; format_version = 1, state_kind = Symbol(state_kind),
        hamiltonian_kind = Symbol(hamiltonian_kind), molecular_charge, spin_2s,
        electron_count, nalpha, nbeta, nuclear_charges, nuclear_positions,
        payload_sha256, checkpoint_sha256, mo_indices,
        producer = (name = String(producer["name"]), version = String(producer["version"]),
            exporter = String(producer["exporter"]), exporter_version = String(producer["exporter_version"])))
    alpha, beta = blocks[1], length(blocks) == 2 ? blocks[2] : nothing
    return ExternalGTOOrbitalPacket(probes, S_GG, alpha; beta,
        ao_labels = [orbital.label for orbital in probes.orbitals],
        provenance = (; external_cartesian_gto = identity))
end

function _xgto_closest_spin(imported, occupations, expected_occupation, threshold, name)
    !isempty(occupations) && all(abs(value - expected_occupation) <= 1.0e-10 for value in occupations) ||
        _xgto_fail("$(name) occupations do not define the required determinant")
    coefficients = imported.imported_coefficients
    gram = Symmetric(transpose(coefficients) * coefficients)
    values = eigvals(gram)
    all(isfinite, values) && minimum(values) >= -1.0e-10 && maximum(values) <= 1.0 + 1.0e-10 ||
        _xgto_fail("$(name) imported Gram spectrum is invalid")
    minimum(values) >= threshold ||
        _xgto_fail("$(name) minimum Gram eigenvalue $(minimum(values)) is below $(threshold)")
    cleaned = coefficients * inv(sqrt(gram))
    orthonormality_error = norm(transpose(cleaned) * cleaned - I, Inf)
    orthonormality_error <= 1.0e-10 || _xgto_fail("$(name) closest determinant is not orthonormal")
    angles = acos.(sqrt.(clamp.(values, 0.0, 1.0)))
    return (; coefficients = cleaned, gram_eigenvalues = values,
        minimum_gram_eigenvalue = minimum(values), principal_angles_radians = angles,
        maximum_principal_angle_radians = maximum(angles), orthonormality_error)
end

"""
    closest_external_gto_determinant(working, packet; minimum_gram_eigenvalue)

Import an external Cartesian-GTO packet and construct the closest full-rank
orthonormal determinant in the same working basis. The raw projected import is
returned unchanged beside the cleaned spin blocks. The caller must choose the
minimum acceptable occupied-Gram eigenvalue.
"""
function closest_external_gto_determinant(working, packet::ExternalGTOOrbitalPacket;
    minimum_gram_eigenvalue)
    threshold = _xgto_float(minimum_gram_eigenvalue, "minimum_gram_eigenvalue")
    0 < threshold <= 1 || _xgto_fail("minimum_gram_eigenvalue must lie in (0,1]")
    hasproperty(working, :hamiltonian) || _xgto_fail("working object must expose its same-construction Hamiltonian")
    hamiltonian = working.hamiltonian
    hamiltonian isa CartesianIDAHamiltonian || _xgto_fail("working Hamiltonian must be CartesianIDAHamiltonian")
    hasproperty(packet.provenance, :external_cartesian_gto) ||
        _xgto_fail("packet lacks versioned interchange identity")
    identity = packet.provenance.external_cartesian_gto
    identity.format_version == 1 && identity.hamiltonian_kind == :all_electron_nuclear_coulomb ||
        _xgto_fail("packet Hamiltonian identity is unsupported")
    hamiltonian.nuclear_charges == identity.nuclear_charges || _xgto_fail("nuclear charges do not match")
    size(hamiltonian.nuclear_positions) == size(identity.nuclear_positions) &&
        norm(hamiltonian.nuclear_positions - identity.nuclear_positions, Inf) <= 1.0e-10 ||
        _xgto_fail("nuclear positions do not match")
    imported = import_external_gto_orbitals(working, packet)
    size(imported.alpha.imported_coefficients, 1) == size(hamiltonian.kinetic, 1) ||
        _xgto_fail("working Hamiltonian/final-basis dimensions disagree")
    if packet.alpha.spin == :restricted
        packet.beta === nothing && identity.state_kind == :rhf || _xgto_fail("restricted packet state mismatch")
        norbitals = length(packet.alpha.occupations)
        hamiltonian.nup == hamiltonian.ndn == identity.nalpha == identity.nbeta == norbitals ||
            _xgto_fail("restricted packet/working electron sectors disagree")
        alpha = _xgto_closest_spin(imported.alpha, packet.alpha.occupations, 2.0, threshold, "restricted")
        return (; imported, alpha, beta = nothing)
    end
    packet.alpha.spin == :alpha || _xgto_fail("spin-resolved packet requires alpha block")
    nalpha = length(packet.alpha.occupations)
    nbeta = isnothing(packet.beta) ? 0 : length(packet.beta.occupations)
    hamiltonian.nup == identity.nalpha == nalpha && hamiltonian.ndn == identity.nbeta == nbeta ||
        _xgto_fail("spin-resolved packet/working electron sectors disagree")
    alpha = _xgto_closest_spin(imported.alpha, packet.alpha.occupations, 1.0, threshold, "alpha")
    beta = isnothing(packet.beta) ? nothing :
        _xgto_closest_spin(imported.beta, packet.beta.occupations, 1.0, threshold, "beta")
    return (; imported, alpha, beta)
end
