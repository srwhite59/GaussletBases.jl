"""
    CartesianIDAHamiltonian(
        kinetic,
        nuclear_attraction_unit_by_center,
        electron_electron_ida,
        nup,
        ndn;
        nuclear_charges,
        nuclear_positions,
    )

One-basis Cartesian IDA Hamiltonian data. The basis is fixed and localized:
`kinetic`, each uncharged unit-nuclear attraction matrix, and
`electron_electron_ida` all share the same `n x n` basis. Dense operator
matrices are treated as owned/read-only by the Hamiltonian object. Nuclear
positions are stored as `ncenter x 3`, and nuclear repulsion is derived from the
stored physical charges and positions.
"""
struct CartesianIDAHamiltonian{T}
    kinetic::Matrix{T}
    nuclear_attraction_unit_by_center::Vector{Matrix{T}}
    electron_electron_ida::Matrix{T}
    nup::Int
    ndn::Int
    nuclear_charges::Vector{T}
    nuclear_positions::Matrix{T}
    nuclear_repulsion::T
end

function _cartesian_nuclear_position_matrix(positions)
    rows = collect(positions)
    matrix = zeros(Float64, length(rows), 3)
    for (row, position) in pairs(rows)
        length(position) == 3 ||
            throw(DimensionMismatch("nuclear positions must have three coordinates"))
        matrix[row, :] .= Float64.(position)
    end
    return matrix
end

_cartesian_dense_float_matrix(matrix::Matrix{Float64}) = matrix
_cartesian_dense_float_matrix(matrix) = Matrix{Float64}(matrix)

_cartesian_float_vector(values::Vector{Float64}) = copy(values)

function _cartesian_float_vector(values)
    return Float64[Float64(value) for value in values]
end

function _cartesian_nuclear_position_matrix(positions::Matrix{Float64})
    size(positions, 2) == 3 ||
        throw(DimensionMismatch("nuclear position matrix must have three columns"))
    return copy(positions)
end

function _cartesian_nuclear_position_matrix(positions::AbstractMatrix{<:Real})
    size(positions, 2) == 3 ||
        throw(DimensionMismatch("nuclear position matrix must have three columns"))
    return _cartesian_dense_float_matrix(positions)
end

function _cartesian_check_symmetric_finite_matrix(
    name::AbstractString,
    matrix::AbstractMatrix{<:Real},
)
    dense = _cartesian_dense_float_matrix(matrix)
    size(dense, 1) == size(dense, 2) ||
        throw(DimensionMismatch("$(name) must be square"))
    all(isfinite, dense) ||
        throw(ArgumentError("$(name) contains non-finite entries"))
    norm(dense - transpose(dense), Inf) <= 1.0e-8 ||
        throw(ArgumentError("$(name) must be symmetric"))
    return dense
end

function CartesianIDAHamiltonian(
    kinetic,
    nuclear_attraction_unit_by_center,
    electron_electron_ida,
    nup::Integer,
    ndn::Integer;
    nuclear_charges,
    nuclear_positions,
)
    kinetic_matrix =
        _cartesian_check_symmetric_finite_matrix("kinetic energy", kinetic)
    electron_electron_matrix =
        _cartesian_check_symmetric_finite_matrix(
            "IDA electron-electron interaction",
            electron_electron_ida,
        )
    size(electron_electron_matrix) == size(kinetic_matrix) ||
        throw(DimensionMismatch("IDA one-body and electron-electron dimensions differ"))
    center_matrices =
        _cartesian_ida_center_matrices(nuclear_attraction_unit_by_center)
    all(size(matrix) == size(kinetic_matrix) for matrix in center_matrices) ||
        throw(DimensionMismatch("unit nuclear attraction dimensions differ"))
    nup_value = Int(nup)
    ndn_value = Int(ndn)
    nup_value >= 0 && ndn_value >= 0 && nup_value + ndn_value > 0 ||
        throw(ArgumentError("spin-sector electron counts must be non-negative and nonzero"))
    norb = size(kinetic_matrix, 1)
    nup_value <= norb && ndn_value <= norb ||
        throw(ArgumentError("spin-sector electron counts must not exceed orbital dimension"))
    charges = _cartesian_float_vector(nuclear_charges)
    all(isfinite, charges) ||
        throw(ArgumentError("nuclear charges contain non-finite entries"))
    length(center_matrices) == length(charges) ||
        throw(DimensionMismatch("unit nuclear attraction count must match charges"))
    positions = _cartesian_nuclear_position_matrix(nuclear_positions)
    size(positions, 1) == length(charges) ||
        throw(DimensionMismatch("nuclear position row count must match charges"))
    all(isfinite, positions) ||
        throw(ArgumentError("nuclear positions contain non-finite entries"))
    repulsion = _cartesian_ida_nuclear_repulsion(charges, positions)
    return CartesianIDAHamiltonian{Float64}(
        kinetic_matrix,
        center_matrices,
        electron_electron_matrix,
        nup_value,
        ndn_value,
        charges,
        positions,
        repulsion,
    )
end

"""
    one_body_hamiltonian(ham::CartesianIDAHamiltonian; center_weights = ones(ncenter))

Assemble `K + sum_A center_weights[A] * Z_A * U_A` in the Hamiltonian's
localized IDA basis.
"""
function one_body_hamiltonian(
    ham::CartesianIDAHamiltonian;
    center_weights = ones(length(ham.nuclear_charges)),
)
    weights = _cartesian_ida_center_weights(center_weights, ham.nuclear_charges)
    matrix = copy(ham.kinetic)
    for (weight, charge, nuclear) in
        zip(weights, ham.nuclear_charges, ham.nuclear_attraction_unit_by_center)
        matrix .+= weight * charge .* nuclear
    end
    return matrix
end

"""
    nuclear_repulsion(ham::CartesianIDAHamiltonian; center_weights = ones(ncenter))

Return `sum_{A<B} (w_A Z_A)(w_B Z_B)/|R_A - R_B|` using stored physical
charges and positions.
"""
function nuclear_repulsion(
    ham::CartesianIDAHamiltonian;
    center_weights = ones(length(ham.nuclear_charges)),
)
    weights = _cartesian_ida_center_weights(center_weights, ham.nuclear_charges)
    return _cartesian_ida_nuclear_repulsion(
        weights .* ham.nuclear_charges,
        ham.nuclear_positions,
    )
end

function _cartesian_ida_density_matrix(
    ham::CartesianIDAHamiltonian,
    density,
    name::AbstractString,
)
    matrix = _cartesian_dense_float_matrix(density)
    size(matrix) == size(ham.kinetic) ||
        throw(DimensionMismatch("$(name) must have size $(size(ham.kinetic))"))
    all(isfinite, matrix) ||
        throw(ArgumentError("$(name) contains non-finite entries"))
    norm(matrix - transpose(matrix), Inf) <= 1.0e-8 ||
        throw(ArgumentError("$(name) must be symmetric"))
    return 0.5 .* (matrix .+ transpose(matrix))
end

function _cartesian_ida_spin_densities(
    ham::CartesianIDAHamiltonian,
    density_alpha,
    density_beta,
)
    alpha = _cartesian_ida_density_matrix(ham, density_alpha, "density_alpha")
    beta = _cartesian_ida_density_matrix(ham, density_beta, "density_beta")
    return alpha, beta
end

function _cartesian_ida_approximate_interaction_energy(
    ham::CartesianIDAHamiltonian,
    density_alpha,
    density_beta,
)
    alpha, beta = _cartesian_ida_spin_densities(ham, density_alpha, density_beta)
    interaction = ham.electron_electron_ida
    total_row_density = diag(alpha) + diag(beta)
    direct = 0.5 * dot(total_row_density, interaction * total_row_density)
    exchange_alpha = 0.5 * sum(interaction .* alpha .* transpose(alpha))
    exchange_beta = 0.5 * sum(interaction .* beta .* transpose(beta))
    return Float64(direct - exchange_alpha - exchange_beta)
end

function _cartesian_ida_approximate_direct_interaction_energy(
    ham::CartesianIDAHamiltonian,
    density_alpha,
    density_beta,
)
    alpha, beta = _cartesian_ida_spin_densities(ham, density_alpha, density_beta)
    row_density = diag(alpha) + diag(beta)
    return Float64(0.5 * dot(row_density, ham.electron_electron_ida * row_density))
end

function _cartesian_ida_approximate_direct_interaction_fock(
    ham::CartesianIDAHamiltonian,
    density_alpha,
    density_beta,
)
    alpha, beta = _cartesian_ida_spin_densities(ham, density_alpha, density_beta)
    row_potential = ham.electron_electron_ida * (diag(alpha) + diag(beta))
    return Matrix(Diagonal(row_potential))
end

function _cartesian_ida_approximate_interaction_fock(
    ham::CartesianIDAHamiltonian,
    same_spin_density::AbstractMatrix{Float64},
    total_row_density::AbstractVector{Float64},
)
    interaction = ham.electron_electron_ida
    fock = Matrix(Diagonal(interaction * total_row_density)) .-
        interaction .* same_spin_density
    return 0.5 .* (fock .+ transpose(fock))
end

function _cartesian_ida_approximate_interaction_fock_alpha(
    ham::CartesianIDAHamiltonian,
    density_alpha,
    density_beta,
)
    alpha, beta = _cartesian_ida_spin_densities(ham, density_alpha, density_beta)
    return _cartesian_ida_approximate_interaction_fock(
        ham,
        alpha,
        diag(alpha) + diag(beta),
    )
end

function _cartesian_ida_approximate_interaction_fock_beta(
    ham::CartesianIDAHamiltonian,
    density_alpha,
    density_beta,
)
    alpha, beta = _cartesian_ida_spin_densities(ham, density_alpha, density_beta)
    return _cartesian_ida_approximate_interaction_fock(
        ham,
        beta,
        diag(alpha) + diag(beta),
    )
end

function _cartesian_ida_approximate_electronic_energy(
    ham::CartesianIDAHamiltonian,
    density_alpha,
    density_beta,
)
    alpha, beta = _cartesian_ida_spin_densities(ham, density_alpha, density_beta)
    one_body = one_body_hamiltonian(ham)
    return Float64(
        sum(one_body .* (alpha + beta)) +
        _cartesian_ida_approximate_interaction_energy(ham, alpha, beta)
    )
end

function _cartesian_ida_approximate_electronic_fock_alpha(
    ham::CartesianIDAHamiltonian,
    density_alpha,
    density_beta,
)
    return one_body_hamiltonian(ham) +
        _cartesian_ida_approximate_interaction_fock_alpha(
            ham,
            density_alpha,
            density_beta,
        )
end

function _cartesian_ida_approximate_electronic_fock_beta(
    ham::CartesianIDAHamiltonian,
    density_alpha,
    density_beta,
)
    return one_body_hamiltonian(ham) +
        _cartesian_ida_approximate_interaction_fock_beta(
            ham,
            density_alpha,
            density_beta,
        )
end

function _cartesian_ida_center_weights(center_weights, charges)
    weights = _cartesian_float_vector(center_weights)
    length(weights) == length(charges) ||
        throw(DimensionMismatch("center weight count must match centers"))
    all(isfinite, weights) ||
        throw(ArgumentError("center weights contain non-finite entries"))
    return weights
end

function _cartesian_ida_nuclear_repulsion(charges, positions)
    repulsion = 0.0
    for right in 2:length(charges), left in 1:(right - 1)
        distance = norm(positions[right, :] .- positions[left, :])
        distance > 0.0 ||
            throw(ArgumentError("nuclear repulsion requires distinct centers"))
        repulsion += charges[left] * charges[right] / distance
    end
    return repulsion
end

function _cartesian_ida_center_matrices(matrices)
    return Matrix{Float64}[
        _cartesian_check_symmetric_finite_matrix(
            "unit nuclear attraction",
            matrix,
        ) for matrix in matrices
    ]
end

function _cartesian_ida_center_matrices(tensor::AbstractArray{<:Real,3})
    return Matrix{Float64}[
        _cartesian_check_symmetric_finite_matrix(
            "unit nuclear attraction",
            view(tensor, :, :, center),
        ) for center in axes(tensor, 3)
    ]
end

function _cartesian_ida_center_tensor(ham::CartesianIDAHamiltonian)
    norb = size(ham.kinetic, 1)
    center_count = length(ham.nuclear_attraction_unit_by_center)
    tensor = Array{Float64}(undef, norb, norb, center_count)
    for center in 1:center_count
        tensor[:, :, center] .= ham.nuclear_attraction_unit_by_center[center]
    end
    return tensor
end

const _CARTESIAN_COULOMB_EXPANSION_SUMMARY_KEYS =
    (:policy, :doacc, :term_count, :del, :s, :c, :maxu)

function _cartesian_coulomb_expansion_preset_summary(policy)
    (policy isa Symbol || policy isa AbstractString) ||
        throw(ArgumentError("Coulomb expansion policy must be :compact or :high"))
    value = Symbol(policy)
    value === :compact && return (;
        policy = value, doacc = false, term_count = 45,
        del = 0.6, s = 0.5, c = 0.03, maxu = 27.0)
    value === :high && return (;
        policy = value, doacc = true, term_count = 135,
        del = 1.0, s = 0.16, c = 0.01, maxu = 135.0)
    throw(ArgumentError("Coulomb expansion policy must be :compact or :high"))
end

function _cartesian_validate_coulomb_expansion_summary(summary)
    all(name -> hasproperty(summary, name), _CARTESIAN_COULOMB_EXPANSION_SUMMARY_KEYS) ||
        throw(ArgumentError("Coulomb expansion summary is incomplete"))
    expected = _cartesian_coulomb_expansion_preset_summary(summary.policy)
    actual = (;
        policy = Symbol(summary.policy),
        doacc = summary.doacc,
        term_count = Int(summary.term_count),
        del = Float64(summary.del),
        s = Float64(summary.s),
        c = Float64(summary.c),
        maxu = Float64(summary.maxu))
    actual.doacc isa Bool || throw(ArgumentError("Coulomb expansion doacc must be Bool"))
    actual == expected || throw(ArgumentError(
        "Coulomb expansion summary does not match the named deterministic preset"))
    return actual
end

function _cartesian_coulomb_expansion_summary(policy, expansion::CoulombGaussianExpansion)
    summary = (;
        policy = Symbol(policy),
        doacc = Symbol(policy) === :high,
        term_count = length(expansion),
        del = expansion.del,
        s = expansion.s,
        c = expansion.c,
        maxu = expansion.maxu)
    return _cartesian_validate_coulomb_expansion_summary(summary)
end

function _cartesian_write_coulomb_expansion_summary!(file, summary;
    prefix::AbstractString = "coulomb_expansion")
    value = _cartesian_validate_coulomb_expansion_summary(summary)
    for name in _CARTESIAN_COULOMB_EXPANSION_SUMMARY_KEYS
        file["$(prefix)/$(name)"] = getproperty(value, name)
    end
    return value
end

function _cartesian_read_coulomb_expansion_summary(file;
    prefix::AbstractString = "coulomb_expansion")
    present = map(name -> haskey(file, "$(prefix)/$(name)"),
        _CARTESIAN_COULOMB_EXPANSION_SUMMARY_KEYS)
    any(present) || return nothing
    all(present) || throw(ArgumentError("Coulomb expansion summary is incomplete"))
    summary = (; (name => file["$(prefix)/$(name)"] for
        name in _CARTESIAN_COULOMB_EXPANSION_SUMMARY_KEYS)...)
    return _cartesian_validate_coulomb_expansion_summary(summary)
end

"""
    write_cartesian_ida_hamiltonian(path, ham::CartesianIDAHamiltonian)

Write the minimal versioned Cartesian IDA Hamiltonian JLD2 artifact.
"""
function write_cartesian_ida_hamiltonian(path, ham::CartesianIDAHamiltonian)
    jldopen(String(path), "w") do file
        file["artifact_kind"] = :cartesian_ida_hamiltonian
        file["format_version"] = 1
        file["kinetic"] = ham.kinetic
        file["nuclear_attraction_unit_by_center"] =
            _cartesian_ida_center_tensor(ham)
        file["electron_electron_ida"] = ham.electron_electron_ida
        file["nuclear_charges"] = ham.nuclear_charges
        file["nuclear_positions"] = ham.nuclear_positions
        file["nup"] = ham.nup
        file["ndn"] = ham.ndn
    end
    return path
end

"""
    read_cartesian_ida_hamiltonian(path)

Read a minimal Cartesian IDA Hamiltonian JLD2 artifact.
"""
function read_cartesian_ida_hamiltonian(path)
    return jldopen(String(path), "r") do file
        file["artifact_kind"] === :cartesian_ida_hamiltonian ||
            throw(ArgumentError("artifact_kind must be :cartesian_ida_hamiltonian"))
        Int(file["format_version"]) == 1 ||
            throw(ArgumentError("unsupported Cartesian IDA Hamiltonian format version"))
        CartesianIDAHamiltonian(
            file["kinetic"],
            file["nuclear_attraction_unit_by_center"],
            file["electron_electron_ida"],
            Int(file["nup"]),
            Int(file["ndn"]);
            nuclear_charges = file["nuclear_charges"],
            nuclear_positions = file["nuclear_positions"],
        )
    end
end

const _PROTECTED_LOCALIZED_ARTIFACT_KIND = :protected_localized_inherited_site_ida_hamiltonian
const _PROTECTED_LOCALIZED_CONVENTION_ID = :protected_localized_inherited_site_ida_v1

function _protected_localized_int(value, name::AbstractString)
    out = Int(value)
    out >= 0 || throw(ArgumentError("$(name) must be non-negative"))
    return out
end

function _protected_localized_prop(values, name::Symbol, context::AbstractString)
    hasproperty(values, name) || throw(ArgumentError("$(context) is missing $(name)"))
    return getproperty(values, name)
end

function _protected_localized_write_simple_group(file, prefix::AbstractString, values)
    for name in propertynames(values)
        file["$(prefix)/$(name)"] = getproperty(values, name)
    end
    return nothing
end

function _protected_localized_sector_counts(sector_counts, final_dimension::Int)
    value(name) = _protected_localized_int(
        _protected_localized_prop(sector_counts, name, "sector_counts"),
        "sector_counts.$(name)",
    )
    base, compact_R = value(:base), value(:compact_R)
    protected_Z, broad_Z = value(:protected_Z), value(:broad_Z)
    M = base + compact_R
    Z = protected_Z + broad_Z
    Qperp = final_dimension - Z
    M == final_dimension ||
        throw(DimensionMismatch("protected-localized final dimension must equal base + compact_R"))
    Qperp >= 0 ||
        throw(DimensionMismatch("protected-localized Z dimension exceeds final dimension"))
    return (; base, compact_R, M, protected_Z, broad_Z, Z, Qperp,
        final_dimension)
end

function _protected_localized_write_sector_metadata(file, counts)
    ranges = (;
        base_G = (1, counts.base),
        compact_R = (counts.base + 1, counts.M),
        protected_Z = (1, counts.protected_Z),
        broad_Z = (counts.protected_Z + 1, counts.Z),
        Z = (1, counts.Z),
        Qperp = (counts.Z + 1, counts.final_dimension),
        L = (1, counts.final_dimension))
    for (name, (first, last)) in pairs(ranges)
        file["sector_ranges/$(name)"] = Int[first, last]
    end
    return nothing
end

function _protected_localized_ordering(localized_ordering)
    return merge((;
            order = :inherited_M_site_order,
            transform = :angular_style_projected_localization,
            parent_M_order = :base_G_then_compact_R,
            fixed_F_order = :protected_Z_then_broad_Z_then_Qperp,
            interaction = :inherited_pre_injection_site_order,
            exact_one_body = :protected_localized_dense_transform),
        localized_ordering)
end

function _protected_localized_read_key(file, key::AbstractString)
    try
        return file[key]
    catch err
        throw(ArgumentError("missing protected-localized artifact key $(key)"))
    end
end
_protected_localized_read_symbol(file, key::AbstractString) =
    Symbol(_protected_localized_read_key(file, key))

function _protected_localized_read_matrix(file, key::AbstractString, dimension::Int)
    matrix = _cartesian_check_symmetric_finite_matrix(key, file[key])
    size(matrix) == (dimension, dimension) ||
        throw(DimensionMismatch("$(key) must have size $((dimension, dimension))"))
    return matrix
end

function _protected_localized_read_counts(file, dimension::Int)
    counts = (; (name => _protected_localized_int(
        file["sector_counts/$(name)"], String(name)) for
        name in (:base, :compact_R, :protected_Z, :broad_Z))...)
    return _protected_localized_sector_counts(counts, dimension)
end

function _protected_localized_read_diagnostics(file)
    float_keys = (:B_min, :B_median, :B_max, :F_S_F_identity_error,
        :Z_S_M_Qperp_error, :Qperp_identity_error, :protected_span_min_sv,
        :L_identity_error, :M_L_diag_delta_max, :M_L_offdiag_max,
        :M_L_fro_delta, :H1_L_symmetry_error, :Vee_L_symmetry_error)
    int_keys = (:B_lt_0p999, :B_lt_0p99, :B_lt_0p98, :B_lt_0p95, :B_lt_0p9)
    return merge(
        NamedTuple{float_keys}(Tuple(Float64(file["diagnostics/$(key)"]) for key in float_keys)),
        NamedTuple{int_keys}(Tuple(Int(file["diagnostics/$(key)"]) for key in int_keys)),
    )
end

function _protected_localized_read_ranges(file)
    range_keys = (:base_G, :compact_R, :protected_Z, :broad_Z, :Z, :Qperp, :L)
    read_range(key) = begin
        values = Int.(file["sector_ranges/$(key)"])
        length(values) == 2 ||
            throw(ArgumentError("sector range $(key) must have two endpoints"))
        (; first = values[1], last = values[2])
    end
    return (; (key => read_range(key) for key in range_keys)...)
end

function _protected_localized_read_metadata(file)
    basis_controls = (;
        nesting = Symbol(file["basis_controls/nesting"]),
        ns = Int(file["basis_controls/ns"]),
        core_spacing = Float64(file["basis_controls/core_spacing"]),
        basisname = String(file["basis_controls/basisname"]),
        lmax = Int(file["basis_controls/lmax"]))
    geometry_inputs = (;
        atom_symbols = String.(file["geometry_inputs/atom_symbols"]),
        nuclear_charges = Float64.(file["geometry_inputs/nuclear_charges"]),
        atom_locations = file["geometry_inputs/atom_locations"])
    localized_ordering = (;
        order = Symbol(file["localized_ordering/order"]),
        transform = Symbol(file["localized_ordering/transform"]),
        parent_M_order = Symbol(file["localized_ordering/parent_M_order"]),
        fixed_F_order = Symbol(file["localized_ordering/fixed_F_order"]),
        interaction = Symbol(file["localized_ordering/interaction"]),
        exact_one_body = Symbol(file["localized_ordering/exact_one_body"]))
    return (; basis_controls, geometry_inputs, localized_ordering)
end

function _protected_localized_float_vector(values, name::AbstractString, dimension::Int)
    vector = Float64.(collect(values))
    length(vector) == dimension ||
        throw(DimensionMismatch("$(name) length must equal final dimension"))
    all(isfinite, vector) || throw(ArgumentError("$(name) must be finite"))
    return vector
end

function _protected_localized_permutation(values, name::AbstractString, dimension::Int)
    vector = Int.(collect(values))
    length(vector) == dimension ||
        throw(DimensionMismatch("$(name) length must equal final dimension"))
    sort(vector) == collect(1:dimension) ||
        throw(ArgumentError("$(name) must be a permutation of 1:final_dimension"))
    return vector
end

function _protected_localized_optional_prop(values, name::Symbol)
    return hasproperty(values, name) ? getproperty(values, name) : nothing
end

function _protected_localized_validate_row_locality(row_locality, counts)
    dimension = counts.final_dimension
    center_x = _protected_localized_float_vector(row_locality.center_x,
        "row_locality.center_x", dimension)
    center_y = _protected_localized_float_vector(row_locality.center_y,
        "row_locality.center_y", dimension)
    center_z = _protected_localized_float_vector(row_locality.center_z,
        "row_locality.center_z", dimension)
    z_order_to_native = _protected_localized_permutation(
        row_locality.z_order_to_native, "row_locality.z_order_to_native", dimension)
    native_to_z_order = _protected_localized_permutation(
        row_locality.native_to_z_order, "row_locality.native_to_z_order", dimension)
    for (z_index, native_index) in pairs(z_order_to_native)
        native_to_z_order[native_index] == z_index ||
            throw(ArgumentError("row_locality permutations are not inverses"))
        if z_index > 1
            previous = z_order_to_native[z_index - 1]
            (center_z[previous], previous) <= (center_z[native_index], native_index) ||
                throw(ArgumentError("row_locality z order is not monotone"))
        end
    end
    sector_label = String.(collect(row_locality.sector_label))
    native_sector_index = Int.(collect(row_locality.native_sector_index))
    length(sector_label) == dimension && length(native_sector_index) == dimension ||
        throw(DimensionMismatch("row_locality sector metadata length mismatch"))
    for row in 1:dimension
        expected_label = row <= counts.base ? "base_G" : "compact_R"
        expected_index = row <= counts.base ? row : row - counts.base
        sector_label[row] == expected_label &&
            native_sector_index[row] == expected_index ||
            throw(ArgumentError("row_locality sector metadata does not match native order"))
    end
    data = (; center_x, center_y, center_z, native_to_z_order,
        z_order_to_native, sector_label, native_sector_index)
    spread_names = (:spread_x, :spread_y, :spread_z)
    spread_values = map(name -> _protected_localized_optional_prop(row_locality, name),
        spread_names)
    all(isnothing, spread_values) && return data
    any(isnothing, spread_values) &&
        throw(ArgumentError("row_locality spreads must be all present or all absent"))
    spreads = NamedTuple{spread_names}(Tuple(_protected_localized_float_vector(
        value, "row_locality.$(name)", dimension) for
        (name, value) in zip(spread_names, spread_values)))
    all(all(>=(0.0), getproperty(spreads, name)) for name in spread_names) ||
        throw(ArgumentError("row_locality spreads must be nonnegative"))
    return merge(data, spreads)
end

function _protected_localized_read_row_locality(file, counts)
    center_x = try
        file["row_locality/center_x"]
    catch
        return nothing
    end
    row_locality = (;
        center_x,
        center_y = file["row_locality/center_y"],
        center_z = file["row_locality/center_z"],
        native_to_z_order = file["row_locality/native_to_z_order"],
        z_order_to_native = file["row_locality/z_order_to_native"],
        sector_label = file["row_locality/sector_label"],
        native_sector_index = file["row_locality/native_sector_index"])
    spread_names = (:spread_x, :spread_y, :spread_z)
    spread_values = map(spread_names) do name
        try
            file["row_locality/$(name)"]
        catch
            nothing
        end
    end
    if any(!isnothing, spread_values)
        all(!isnothing, spread_values) ||
            throw(ArgumentError("row_locality spreads must be all present or all absent"))
        row_locality = merge(row_locality,
            NamedTuple{spread_names}(Tuple(spread_values)))
    end
    return _protected_localized_validate_row_locality(row_locality, counts)
end

function write_protected_localized_ida_hamiltonian(
    path;
    H1_L,
    Vee_L,
    nup::Integer,
    ndn::Integer,
    nuclear_charges,
    nuclear_positions,
    sector_counts,
    diagnostics,
    provenance,
    basis_controls,
    geometry_inputs,
    coulomb_expansion,
    localized_ordering = (;),
    row_locality = nothing,
)
    H1 = _cartesian_check_symmetric_finite_matrix(
        "protected-localized H1_L", H1_L)
    Vee = _cartesian_check_symmetric_finite_matrix(
        "protected-localized inherited-site Vee_L", Vee_L)
    size(H1) == size(Vee) ||
        throw(DimensionMismatch("protected-localized H1_L/Vee_L size mismatch"))
    dimension = size(H1, 1)
    nup_value = _protected_localized_int(nup, "nup")
    ndn_value = _protected_localized_int(ndn, "ndn")
    nup_value <= dimension && ndn_value <= dimension ||
        throw(ArgumentError("protected-localized electron counts exceed dimension"))
    charges = _cartesian_float_vector(nuclear_charges)
    positions = _cartesian_nuclear_position_matrix(nuclear_positions)
    size(positions, 1) == length(charges) ||
        throw(DimensionMismatch("protected-localized center count mismatch"))
    counts = _protected_localized_sector_counts(sector_counts, dimension)
    expansion_summary =
        _cartesian_validate_coulomb_expansion_summary(coulomb_expansion)
    locality = isnothing(row_locality) ? nothing :
        _protected_localized_validate_row_locality(row_locality, counts)
    foreach(key -> _protected_localized_prop(provenance, key, "provenance"),
        (:source_artifact, :source_commit, :current_commit))
    foreach(key -> _protected_localized_prop(basis_controls, key, "basis_controls"),
        (:nesting, :ns, :core_spacing, :basisname, :lmax))
    foreach(key -> _protected_localized_prop(geometry_inputs, key, "geometry_inputs"),
        (:atom_symbols, :nuclear_charges, :atom_locations))
    jldopen(String(path), "w") do file
        file["artifact_kind"] = _PROTECTED_LOCALIZED_ARTIFACT_KIND
        file["format_version"] = 1
        file["convention_id"] = _PROTECTED_LOCALIZED_CONVENTION_ID
        file["convention_version"] = 1
        file["H1_L"] = H1
        file["Vee_L"] = Vee
        file["nup"] = nup_value
        file["ndn"] = ndn_value
        file["final_dimension"] = dimension
        file["nuclear_charges"] = charges
        file["nuclear_positions"] = positions
        file["nuclear_repulsion"] = _cartesian_ida_nuclear_repulsion(charges, positions)
        _protected_localized_write_simple_group(file, "sector_counts", counts)
        _protected_localized_write_sector_metadata(file, counts)
        _protected_localized_write_simple_group(file, "diagnostics", diagnostics)
        _protected_localized_write_simple_group(file, "provenance", provenance)
        _protected_localized_write_simple_group(file, "basis_controls", basis_controls)
        _protected_localized_write_simple_group(file, "geometry_inputs", geometry_inputs)
        _cartesian_write_coulomb_expansion_summary!(file, expansion_summary)
        _protected_localized_write_simple_group(
            file, "localized_ordering", _protected_localized_ordering(localized_ordering))
        if !isnothing(locality)
            _protected_localized_write_simple_group(file, "row_locality", locality)
        end
    end
    return path
end

function read_protected_localized_ida_hamiltonian(path)
    return jldopen(String(path), "r") do file
        _protected_localized_read_symbol(file, "artifact_kind") ===
            _PROTECTED_LOCALIZED_ARTIFACT_KIND ||
            throw(ArgumentError("artifact_kind must be $(_PROTECTED_LOCALIZED_ARTIFACT_KIND)"))
        Int(_protected_localized_read_key(file, "format_version")) == 1 ||
            throw(ArgumentError("unsupported protected-localized format version"))
        _protected_localized_read_symbol(file, "convention_id") ===
            _PROTECTED_LOCALIZED_CONVENTION_ID ||
            throw(ArgumentError("unrecognized protected-localized convention_id"))
        Int(_protected_localized_read_key(file, "convention_version")) == 1 ||
            throw(ArgumentError("unsupported protected-localized convention version"))
        dimension = Int(_protected_localized_read_key(file, "final_dimension"))
        H1 = _protected_localized_read_matrix(file, "H1_L", dimension)
        Vee = _protected_localized_read_matrix(file, "Vee_L", dimension)
        counts = _protected_localized_read_counts(file, dimension)
        nup_value = _protected_localized_int(file["nup"], "nup")
        ndn_value = _protected_localized_int(file["ndn"], "ndn")
        nup_value <= dimension && ndn_value <= dimension ||
            throw(ArgumentError("protected-localized electron counts exceed dimension"))
        charges = _cartesian_float_vector(file["nuclear_charges"])
        positions = _cartesian_nuclear_position_matrix(file["nuclear_positions"])
        diagnostics = _protected_localized_read_diagnostics(file)
        metadata = _protected_localized_read_metadata(file)
        provenance = (;
            source_artifact = String(file["provenance/source_artifact"]),
            source_commit = String(file["provenance/source_commit"]),
            current_commit = String(file["provenance/current_commit"]))
        return (;
            artifact_kind = _PROTECTED_LOCALIZED_ARTIFACT_KIND,
            format_version = 1,
            convention_id = _PROTECTED_LOCALIZED_CONVENTION_ID,
            convention_version = 1,
            H1_L = H1,
            Vee_L = Vee,
            nup = nup_value,
            ndn = ndn_value,
            final_dimension = dimension,
            nuclear_charges = charges,
            nuclear_positions = positions,
            nuclear_repulsion = _cartesian_ida_nuclear_repulsion(charges, positions),
            sector_counts = counts,
            sector_ranges = _protected_localized_read_ranges(file),
            diagnostics,
            provenance,
            coulomb_expansion = _cartesian_read_coulomb_expansion_summary(file),
            row_locality = _protected_localized_read_row_locality(file, counts),
            metadata...)
    end
end
