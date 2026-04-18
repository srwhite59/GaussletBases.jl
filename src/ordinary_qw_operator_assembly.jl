function _qwrg_capture_timeg_report(
    build::F,
    label::AbstractString,
    timing::Bool;
    timing_io::IO = stdout,
) where {F<:Function}
    if !timing
        return @timeg label build()
    end

    old_config = TimeG._TIMING_CONFIG[]
    state = TimeG._timing_state()
    old_children = length(state.root.children)
    try
        set_timing!(true)
        set_timing_live!(false)
        set_timing_thresholds!(expand = 0.0, drop = 0.0)
        result = @timeg label build()
        new_children = state.root.children[(old_children + 1):end]
        length(new_children) == 1 || throw(
            ArgumentError("QW TimeG capture expected exactly one root timing node"),
        )
        timing_report(timing_io, TimeG.TimingReport(TimeG.TimingNode[deepcopy(only(new_children))]))
        return result
    finally
        TimeG._TIMING_CONFIG[] = old_config
    end
end

function _qwrg_print_basis_counts(
    io::IO,
    gausslet_count::Integer,
    raw_to_final::AbstractMatrix{<:Real},
)
    return _qwrg_print_basis_counts(io, "gausslet_count", gausslet_count, raw_to_final)
end

function _qwrg_print_basis_counts(
    io::IO,
    carried_label::AbstractString,
    carried_count::Integer,
    raw_to_final::AbstractMatrix{<:Real},
)
    total_basis_dim = size(raw_to_final, 2)
    residual_count = total_basis_dim - carried_count
    println(
        io,
        "QW-RG basis counts  ",
        carried_label,
        "=",
        carried_count,
        "  residual_count=",
        residual_count,
        "  total_basis_dim=",
        total_basis_dim,
    )
    flush(io)
    return nothing
end

function _validate_nuclear_term_storage(storage::Symbol)
    storage in (:auto, :total_only, :by_center) || throw(
        ArgumentError(
            "nuclear_term_storage must be :auto, :total_only, or :by_center",
        ),
    )
    return storage
end

_resolved_nuclear_term_storage(storage::Symbol, ::BondAlignedDiatomicQWBasis3D) =
    storage == :auto ? :by_center : storage
_resolved_nuclear_term_storage(
    storage::Symbol,
    ::Union{BondAlignedHomonuclearChainQWBasis3D,AxisAlignedHomonuclearSquareLatticeQWBasis3D},
) = storage == :auto ? :total_only : storage
_resolved_nuclear_term_storage(
    storage::Symbol,
    fixed_block::_NestedFixedBlock3D{<:AbstractBondAlignedOrdinaryQWBasis3D},
) = _resolved_nuclear_term_storage(storage, fixed_block.parent_basis)

function _same_nuclear_charge_configuration(
    left::Union{Nothing,AbstractVector{<:Real}},
    right::Union{Nothing,AbstractVector{<:Real}};
    tol::Real = 1.0e-12,
)
    (left === nothing || right === nothing) && return left === right
    length(left) == length(right) || return false
    for index in eachindex(left, right)
        abs(Float64(left[index]) - Float64(right[index])) <= tol || return false
    end
    return true
end

function _assemble_one_body_hamiltonian(
    kinetic_one_body::AbstractMatrix{<:Real},
    nuclear_one_body_by_center::AbstractVector{<:AbstractMatrix{<:Real}},
    nuclear_charges::AbstractVector{<:Real},
)
    length(nuclear_charges) == length(nuclear_one_body_by_center) || throw(
        ArgumentError("one-body reassembly requires one nuclear charge per stored nuclear term"),
    )
    matrix = Matrix{Float64}(kinetic_one_body)
    for (charge, nuclear_matrix) in zip(nuclear_charges, nuclear_one_body_by_center)
        matrix .+= Float64(charge) .* nuclear_matrix
    end
    return Matrix{Float64}(0.5 .* (matrix .+ transpose(matrix)))
end

"""
    assembled_one_body_hamiltonian(operators::OrdinaryCartesianOperators3D;
                                   nuclear_charges = operators.nuclear_charges)

Reassemble the final-basis one-body Hamiltonian from stored kinetic and
per-center nuclear-attraction terms when available.

If the operator payload does not retain per-center nuclear terms, the stored
`one_body_hamiltonian` can only be returned for the original nuclear charges;
requesting modified charges then throws an `ArgumentError`.
"""
function assembled_one_body_hamiltonian(
    operators::OrdinaryCartesianOperators3D;
    nuclear_charges = operators.nuclear_charges,
)
    if isnothing(operators.kinetic_one_body) || isnothing(operators.nuclear_one_body_by_center)
        if nuclear_charges === operators.nuclear_charges ||
           _same_nuclear_charge_configuration(nuclear_charges, operators.nuclear_charges)
            return Matrix{Float64}(operators.one_body_hamiltonian)
        end
        throw(
            ArgumentError(
                "this operator payload does not retain per-center nuclear one-body terms; rebuild with nuclear_term_storage = :by_center",
            ),
        )
    end
    charges = isnothing(nuclear_charges) ? operators.nuclear_charges : nuclear_charges
    isnothing(charges) && return Matrix{Float64}(operators.one_body_hamiltonian)
    return _assemble_one_body_hamiltonian(
        operators.kinetic_one_body,
        operators.nuclear_one_body_by_center,
        charges,
    )
end

function ordinary_cartesian_vee_expectation(
    operators::OrdinaryCartesianOperators3D,
    orbital::AbstractVector;
    overlap_tol::Real = 1.0e-8,
)
    length(orbital) == length(operators.orbital_data) ||
        throw(ArgumentError("orbital length must match the Qiu-White ordinary Cartesian orbital dimension"))
    overlap_error = norm(operators.overlap - I, Inf)
    overlap_error <= Float64(overlap_tol) || throw(
        ArgumentError(
            "ordinary_cartesian_vee_expectation currently requires an orthonormal final basis; got overlap error $(overlap_error)",
        ),
    )
    weights = Float64[abs2(coefficient) for coefficient in orbital]
    norm2 = sum(weights)
    norm2 > 0.0 || throw(ArgumentError("orbital must have nonzero norm"))
    weights ./= norm2
    return Float64(real(dot(weights, operators.interaction_matrix * weights)))
end

function ordinary_cartesian_1s2_check(
    operators::OrdinaryCartesianOperators3D;
    overlap_tol::Real = 1.0e-8,
)
    decomposition = eigen(Hermitian(operators.one_body_hamiltonian))
    orbital = decomposition.vectors[:, 1]
    return (
        orbital_energy = Float64(decomposition.values[1]),
        orbital = orbital,
        vee_expectation = ordinary_cartesian_vee_expectation(
            operators,
            orbital;
            overlap_tol = overlap_tol,
        ),
        overlap_error = norm(operators.overlap - I, Inf),
    )
end

function _qwrg_one_body_matrices(
    gausslet_one_body::AbstractMatrix{<:Real},
    one_body_ga::AbstractMatrix{<:Real},
    one_body_aa::AbstractMatrix{<:Real},
    raw_to_final::AbstractMatrix{<:Real},
)
    raw_one_body = [
        Matrix{Float64}(gausslet_one_body) Matrix{Float64}(one_body_ga)
        Matrix{Float64}(transpose(one_body_ga)) Matrix{Float64}(one_body_aa)
    ]
    final_one_body = Matrix{Float64}(transpose(raw_to_final) * raw_one_body * raw_to_final)
    return raw_one_body, 0.5 .* (final_one_body .+ transpose(final_one_body))
end

function _qwrg_nearest_indices(
    gausslet_orbitals::AbstractVector{<:CartesianProductOrbital3D},
    residual_centers::AbstractMatrix{<:Real},
)
    return Int[
        argmin(
            [
                (orbital.x - residual_centers[index, 1])^2 +
                (orbital.y - residual_centers[index, 2])^2 +
                (orbital.z - residual_centers[index, 3])^2 for orbital in gausslet_orbitals
            ],
        ) for index in axes(residual_centers, 1)
    ]
end

function _qwrg_nearest_indices(
    fixed_centers::AbstractMatrix{<:Real},
    residual_centers::AbstractMatrix{<:Real},
)
    size(fixed_centers, 2) == 3 || throw(
        ArgumentError("nearest-center selection for a nested fixed block requires centers as an n×3 matrix"),
    )
    return Int[
        argmin(
            [
                (fixed_centers[fixed, 1] - residual_centers[index, 1])^2 +
                (fixed_centers[fixed, 2] - residual_centers[index, 2])^2 +
                (fixed_centers[fixed, 3] - residual_centers[index, 3])^2 for fixed in axes(fixed_centers, 1)
            ],
        ) for index in axes(residual_centers, 1)
    ]
end

function _qwrg_effective_gaussians(
    residual_centers::AbstractMatrix{<:Real},
    residual_widths::AbstractMatrix{<:Real},
)
    nresidual = size(residual_centers, 1)
    x_gaussians = Gaussian[]
    y_gaussians = Gaussian[]
    z_gaussians = Gaussian[]
    for index in 1:nresidual
        push!(
            x_gaussians,
            Gaussian(center = residual_centers[index, 1], width = residual_widths[index, 1]),
        )
        push!(
            y_gaussians,
            Gaussian(center = residual_centers[index, 2], width = residual_widths[index, 2]),
        )
        push!(
            z_gaussians,
            Gaussian(center = residual_centers[index, 3], width = residual_widths[index, 3]),
        )
    end
    return x_gaussians, y_gaussians, z_gaussians
end

function _qwrg_same_gaussians(
    left::AbstractVector{<:Gaussian},
    right::AbstractVector{<:Gaussian};
    tol::Real = 1.0e-12,
)
    length(left) == length(right) || return false
    for index in eachindex(left, right)
        abs(left[index].center_value - right[index].center_value) <= tol || return false
        abs(left[index].width - right[index].width) <= tol || return false
    end
    return true
end

function _qwrg_interaction_matrix_nearest(
    gausslet_interaction::AbstractMatrix{<:Real},
    gausslet_orbitals::AbstractVector{<:CartesianProductOrbital3D},
    residual_centers::AbstractMatrix{<:Real},
)
    ngausslet = size(gausslet_interaction, 1)
    nresidual = size(residual_centers, 1)
    interaction = zeros(Float64, ngausslet + nresidual, ngausslet + nresidual)
    interaction[1:ngausslet, 1:ngausslet] .= Matrix{Float64}(gausslet_interaction)
    nearest = _qwrg_nearest_indices(gausslet_orbitals, residual_centers)
    for residual in 1:nresidual
        index = nearest[residual]
        interaction[1:ngausslet, ngausslet + residual] .= gausslet_interaction[:, index]
        interaction[ngausslet + residual, 1:ngausslet] .= gausslet_interaction[:, index]
    end
    for i in 1:nresidual, j in i:nresidual
        value = gausslet_interaction[nearest[i], nearest[j]]
        interaction[ngausslet + i, ngausslet + j] = value
        interaction[ngausslet + j, ngausslet + i] = value
    end
    return interaction
end

function _qwrg_interaction_matrix_nearest(
    fixed_interaction::AbstractMatrix{<:Real},
    fixed_centers::AbstractMatrix{<:Real},
    residual_centers::AbstractMatrix{<:Real},
)
    nfixed = size(fixed_interaction, 1)
    nresidual = size(residual_centers, 1)
    interaction = zeros(Float64, nfixed + nresidual, nfixed + nresidual)
    interaction[1:nfixed, 1:nfixed] .= Matrix{Float64}(fixed_interaction)
    nearest = _qwrg_nearest_indices(fixed_centers, residual_centers)
    for residual in 1:nresidual
        index = nearest[residual]
        interaction[1:nfixed, nfixed + residual] .= fixed_interaction[:, index]
        interaction[nfixed + residual, 1:nfixed] .= fixed_interaction[:, index]
    end
    for i in 1:nresidual, j in i:nresidual
        value = fixed_interaction[nearest[i], nearest[j]]
        interaction[nfixed + i, nfixed + j] = value
        interaction[nfixed + j, nfixed + i] = value
    end
    return interaction
end

# Alg QW-RG step 8b, 8c, and 8d: Match exact residual moments by effective
# Gaussians and keep RG terms in the same two-index IDA representation.
# See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
function _qwrg_interaction_matrix_mwg(
    gausslet_bundle::_MappedOrdinaryGausslet1DBundle,
    gausslet_interaction::AbstractMatrix{<:Real},
    expansion::CoulombGaussianExpansion,
    residual_centers::AbstractMatrix{<:Real},
    residual_widths::AbstractMatrix{<:Real},
)
    ngausslet = size(gausslet_interaction, 1)
    nresidual = size(residual_centers, 1)
    interaction = zeros(Float64, ngausslet + nresidual, ngausslet + nresidual)
    interaction[1:ngausslet, 1:ngausslet] .= Matrix{Float64}(gausslet_interaction)

    x_gaussians, y_gaussians, z_gaussians = _qwrg_effective_gaussians(residual_centers, residual_widths)
    # Alg QW-RG step 8d: Evaluate RG-gausslet terms from exact contracted
    # 1D raw-space pair-factor blocks, then assemble the 3D IDA interaction.
    # See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
    pair_x = _qwrg_split_block_matrices(
        gausslet_bundle,
        x_gaussians,
        expansion,
    )
    pair_y = _qwrg_same_gaussians(x_gaussians, y_gaussians) ? pair_x :
        _qwrg_split_block_matrices(gausslet_bundle, y_gaussians, expansion)
    pair_z = if _qwrg_same_gaussians(x_gaussians, z_gaussians)
        pair_x
    elseif _qwrg_same_gaussians(y_gaussians, z_gaussians)
        pair_y
    else
        _qwrg_split_block_matrices(gausslet_bundle, z_gaussians, expansion)
    end
    analytic_x = _qwrg_gaussian_analytic_blocks(x_gaussians, expansion)
    analytic_y = _qwrg_same_gaussians(x_gaussians, y_gaussians) ? analytic_x : _qwrg_gaussian_analytic_blocks(y_gaussians, expansion)
    analytic_z = if _qwrg_same_gaussians(x_gaussians, z_gaussians)
        analytic_x
    elseif _qwrg_same_gaussians(y_gaussians, z_gaussians)
        analytic_y
    else
        _qwrg_gaussian_analytic_blocks(z_gaussians, expansion)
    end

    for residual in 1:nresidual
        column = zeros(Float64, ngausslet)
        scratch = zeros(Float64, ngausslet)
        for term in eachindex(expansion.coefficients)
            fx = view(pair_x.pair_ga[term], :, residual)
            fy = view(pair_y.pair_ga[term], :, residual)
            fz = view(pair_z.pair_ga[term], :, residual)
            _qwrg_fill_product_column!(scratch, fx, fy, fz)
            column .+= expansion.coefficients[term] .* scratch
        end
        interaction[1:ngausslet, ngausslet + residual] .= column
        interaction[ngausslet + residual, 1:ngausslet] .= column
    end

    for i in 1:nresidual, j in i:nresidual
        value = 0.0
        for term in eachindex(expansion.coefficients)
            value += expansion.coefficients[term] *
                analytic_x.pair_aa[term][i, j] *
                analytic_y.pair_aa[term][i, j] *
                analytic_z.pair_aa[term][i, j]
        end
        interaction[ngausslet + i, ngausslet + j] = value
        interaction[ngausslet + j, ngausslet + i] = value
    end

    return interaction
end

function _qwrg_orbital_data(
    gausslet_orbitals::AbstractVector{<:CartesianProductOrbital3D},
    residual_centers::AbstractMatrix{<:Real},
    residual_widths::AbstractMatrix{<:Real},
)
    orbitals_out = OrdinaryCartesianOrbital3D[]
    for orbital in gausslet_orbitals
        push!(
            orbitals_out,
            OrdinaryCartesianOrbital3D(
                orbital.index,
                :gausslet,
                "g($(orbital.ix),$(orbital.iy),$(orbital.iz))",
                orbital.x,
                orbital.y,
                orbital.z,
                NaN,
                NaN,
                NaN,
            ),
        )
    end
    base_index = length(gausslet_orbitals)
    for index in axes(residual_centers, 1)
        push!(
            orbitals_out,
            OrdinaryCartesianOrbital3D(
                base_index + index,
                :residual_gaussian,
                "rg$index",
                residual_centers[index, 1],
                residual_centers[index, 2],
                residual_centers[index, 3],
                residual_widths[index, 1],
                residual_widths[index, 2],
                residual_widths[index, 3],
            ),
        )
    end
    return orbitals_out
end

function _qwrg_orbital_data(
    fixed_centers::AbstractMatrix{<:Real},
    residual_centers::AbstractMatrix{<:Real},
    residual_widths::AbstractMatrix{<:Real};
    fixed_kind::Symbol = :nested_fixed,
    fixed_label_prefix::AbstractString = "nf",
)
    size(fixed_centers, 2) == 3 || throw(
        ArgumentError("nested fixed-block orbital data requires an n×3 center matrix"),
    )
    orbitals_out = OrdinaryCartesianOrbital3D[]
    for index in axes(fixed_centers, 1)
        push!(
            orbitals_out,
            OrdinaryCartesianOrbital3D(
                index,
                fixed_kind,
                string(fixed_label_prefix, index),
                fixed_centers[index, 1],
                fixed_centers[index, 2],
                fixed_centers[index, 3],
                NaN,
                NaN,
                NaN,
            ),
        )
    end
    base_index = size(fixed_centers, 1)
    for index in axes(residual_centers, 1)
        push!(
            orbitals_out,
            OrdinaryCartesianOrbital3D(
                base_index + index,
                :residual_gaussian,
                "rg$index",
                residual_centers[index, 1],
                residual_centers[index, 2],
                residual_centers[index, 3],
                residual_widths[index, 1],
                residual_widths[index, 2],
                residual_widths[index, 3],
            ),
        )
    end
    return orbitals_out
end

function _qwrg_diatomic_one_body_matrix(
    basis::AbstractBondAlignedOrdinaryQWBasis3D,
    bundle_x::_MappedOrdinaryGausslet1DBundle,
    bundle_y::_MappedOrdinaryGausslet1DBundle,
    bundle_z::_MappedOrdinaryGausslet1DBundle,
    expansion::CoulombGaussianExpansion,
    nuclear_charges::AbstractVector{<:Real},
)
    length(nuclear_charges) == length(basis.nuclei) || throw(
        ArgumentError("bond-aligned ordinary QW path requires one nuclear charge per nucleus"),
    )
    all(Float64(value) >= 0.0 for value in nuclear_charges) || throw(
        ArgumentError("bond-aligned ordinary QW path requires nonnegative nuclear charges"),
    )
    kinetic_one_body = _qwrg_diatomic_kinetic_matrix(bundle_x, bundle_y, bundle_z)
    nuclear_one_body_by_center = _qwrg_diatomic_nuclear_one_body_by_center(
        basis,
        bundle_x,
        bundle_y,
        bundle_z,
        expansion,
    )
    return _assemble_one_body_hamiltonian(
        kinetic_one_body,
        nuclear_one_body_by_center,
        nuclear_charges,
    )
end

function _ordinary_cartesian_qiu_white_operators_bond_aligned_ordinary(
    basis::AbstractBondAlignedOrdinaryQWBasis3D;
    nuclear_charges::AbstractVector{<:Real},
    nuclear_term_storage::Symbol,
    expansion::CoulombGaussianExpansion,
    interaction_treatment::Symbol,
    gausslet_backend::Symbol,
    timing::Bool,
)
    gausslet_backend == :numerical_reference || throw(
        ArgumentError("bond-aligned ordinary_cartesian_qiu_white_operators currently supports only gausslet_backend = :numerical_reference"),
    )
    interaction_treatment in (:ggt_nearest, :mwg) || throw(
        ArgumentError("bond-aligned ordinary_cartesian_qiu_white_operators requires interaction_treatment = :ggt_nearest or :mwg"),
    )
    return _qwrg_capture_timeg_report("qwrg.bond_aligned_ordinary.total", timing) do
        resolved_nuclear_term_storage = _resolved_nuclear_term_storage(
            _validate_nuclear_term_storage(nuclear_term_storage),
            basis,
        )

        bundle_x = _mapped_ordinary_gausslet_1d_bundle(
            basis.basis_x;
            exponents = expansion.exponents,
            center = 0.0,
            backend = gausslet_backend,
        )
        bundle_y = _mapped_ordinary_gausslet_1d_bundle(
            basis.basis_y;
            exponents = expansion.exponents,
            center = 0.0,
            backend = gausslet_backend,
        )
        bundle_z = _mapped_ordinary_gausslet_1d_bundle(
            basis.basis_z;
            exponents = expansion.exponents,
            center = 0.0,
            backend = gausslet_backend,
        )

        overlap = _qwrg_diatomic_overlap_matrix(bundle_x, bundle_y, bundle_z)
        kinetic_one_body = _qwrg_diatomic_kinetic_matrix(bundle_x, bundle_y, bundle_z)
        nuclear_one_body_by_center =
            resolved_nuclear_term_storage == :by_center ?
            _qwrg_diatomic_nuclear_one_body_by_center(
                basis,
                bundle_x,
                bundle_y,
                bundle_z,
                expansion,
            ) :
            nothing
        one_body_hamiltonian =
            isnothing(nuclear_one_body_by_center) ?
            _qwrg_diatomic_one_body_matrix(
                basis,
                bundle_x,
                bundle_y,
                bundle_z,
                expansion,
                nuclear_charges,
            ) :
            _assemble_one_body_hamiltonian(
                kinetic_one_body,
                nuclear_one_body_by_center,
                nuclear_charges,
            )
        interaction_matrix = _qwrg_diatomic_interaction_matrix(
            bundle_x,
            bundle_y,
            bundle_z,
            expansion,
        )

        gausslet_orbitals = _mapped_cartesian_orbitals(
            centers(basis.basis_x),
            centers(basis.basis_y),
            centers(basis.basis_z),
        )
        gausslet_count = length(gausslet_orbitals)
        zero_residual_centers = zeros(Float64, 0, 3)
        zero_residual_widths = zeros(Float64, 0, 3)

        return OrdinaryCartesianOperators3D(
            basis,
            nothing,
            gausslet_backend,
            interaction_treatment,
            expansion,
            overlap,
            one_body_hamiltonian,
            interaction_matrix,
            _qwrg_orbital_data(
                gausslet_orbitals,
                zero_residual_centers,
                zero_residual_widths,
            ),
            gausslet_count,
            0,
            Matrix{Float64}(I, gausslet_count, gausslet_count),
            zero_residual_centers,
            zero_residual_widths,
            Float64[Float64(value) for value in nuclear_charges],
            resolved_nuclear_term_storage == :by_center ? Matrix{Float64}(kinetic_one_body) : nothing,
            resolved_nuclear_term_storage == :by_center ? nuclear_one_body_by_center : nothing,
            resolved_nuclear_term_storage,
        )
    end
end

"""
    ordinary_cartesian_qiu_white_operators(
        basis::BondAlignedDiatomicQWBasis3D;
        nuclear_charges = fill(1.0, length(basis.nuclei)),
        expansion = coulomb_gaussian_expansion(doacc = false),
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        timing = false,
    )

Build the first bond-aligned diatomic ordinary QW reference Hamiltonian on the
diatomic distortion path.

This first molecular pass is intentionally narrow:

- one bond-aligned homonuclear diatomic basis object
- no nested fixed block yet
- no molecular Gaussian supplement yet
- the final basis is the distorted 3D gausslet product basis itself, so the
  residual-Gaussian sector is empty
"""
function ordinary_cartesian_qiu_white_operators(
    basis::BondAlignedDiatomicQWBasis3D;
    nuclear_charges::AbstractVector{<:Real} = fill(1.0, length(basis.nuclei)),
    nuclear_term_storage::Symbol = :auto,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    timing::Bool = false,
)
    return _ordinary_cartesian_qiu_white_operators_bond_aligned_ordinary(
        basis;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        timing = timing,
    )
end

"""
    ordinary_cartesian_qiu_white_operators(
        basis::BondAlignedHomonuclearChainQWBasis3D;
        nuclear_charges = fill(1.0, length(basis.nuclei)),
        expansion = coulomb_gaussian_expansion(doacc = false),
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        timing = false,
    )

Build the first experimental ordinary QW reference Hamiltonian on the
homonuclear chain distortion path.

This chain milestone is intentionally narrow:

- homonuclear linear chains only
- ordinary QW product basis only
- no nested chain splitting
- no molecular-shell supplement route
"""
function ordinary_cartesian_qiu_white_operators(
    basis::BondAlignedHomonuclearChainQWBasis3D;
    nuclear_charges::AbstractVector{<:Real} = fill(1.0, length(basis.nuclei)),
    nuclear_term_storage::Symbol = :auto,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    timing::Bool = false,
)
    return _ordinary_cartesian_qiu_white_operators_bond_aligned_ordinary(
        basis;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        timing = timing,
    )
end

"""
    ordinary_cartesian_qiu_white_operators(
        basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D;
        nuclear_charges = fill(1.0, length(basis.nuclei)),
        expansion = coulomb_gaussian_expansion(doacc = false),
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        timing = false,
    )

Build the first experimental ordinary QW reference Hamiltonian on the
homonuclear square-lattice distortion path.

This lattice milestone is intentionally narrow:

- homonuclear `n × n` square lattices only
- ordinary QW product basis only
- no nested planar splitting
- no molecular-shell supplement route
"""
function ordinary_cartesian_qiu_white_operators(
    basis::AxisAlignedHomonuclearSquareLatticeQWBasis3D;
    nuclear_charges::AbstractVector{<:Real} = fill(1.0, length(basis.nuclei)),
    nuclear_term_storage::Symbol = :auto,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    timing::Bool = false,
)
    return _ordinary_cartesian_qiu_white_operators_bond_aligned_ordinary(
        basis;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        timing = timing,
    )
end

function _ordinary_cartesian_qiu_white_operators_diatomic_shell_3d(
    basis::BondAlignedDiatomicQWBasis3D,
    gaussian_data::Union{
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    };
    nuclear_charges::AbstractVector{<:Real},
    nuclear_term_storage::Symbol,
    expansion::CoulombGaussianExpansion,
    interaction_treatment::Symbol,
    gausslet_backend::Symbol,
    timing::Bool,
)
    gausslet_backend == :numerical_reference || throw(
        ArgumentError("bond-aligned diatomic molecular QW path currently supports only gausslet_backend = :numerical_reference"),
    )
    interaction_treatment == :ggt_nearest || throw(
        ArgumentError("bond-aligned diatomic molecular QW path currently supports only interaction_treatment = :ggt_nearest"),
    )
    length(nuclear_charges) == length(basis.nuclei) || throw(
        ArgumentError("bond-aligned diatomic molecular QW path requires one nuclear charge per nucleus"),
    )
    _qwrg_same_nuclei(gaussian_data.nuclei, basis.nuclei) || throw(
        ArgumentError("bond-aligned diatomic molecular supplement nuclei must match the bond-aligned basis nuclei"),
    )
    return _qwrg_capture_timeg_report("qwrg.diatomic_shell.total", timing) do
        resolved_nuclear_term_storage = _resolved_nuclear_term_storage(
            _validate_nuclear_term_storage(nuclear_term_storage),
            basis,
        )

        bundles = @timeg "qwrg.diatomic_shell.shared_bundles" begin
            _qwrg_bond_aligned_axis_bundles(
                basis,
                expansion;
                gausslet_backend = gausslet_backend,
            )
        end

        supplement3d = _bond_aligned_diatomic_cartesian_shell_supplement_3d(gaussian_data)

        blocks = @timeg "qwrg.diatomic_shell.raw_blocks" begin
            _qwrg_diatomic_cartesian_shell_blocks_3d(
                bundles,
                supplement3d,
                basis,
                expansion,
                nuclear_charges,
            )
        end

        gausslet_orbitals = _mapped_cartesian_orbitals(
            centers(basis.basis_x),
            centers(basis.basis_y),
            centers(basis.basis_z),
        )
        gausslet_count = length(gausslet_orbitals)

        residual_data = @timeg "qwrg.diatomic_shell.residual_space" begin
            gausslet_overlap_3d = _qwrg_diatomic_overlap_matrix(
                bundles.bundle_x,
                bundles.bundle_y,
                bundles.bundle_z,
            )
            _qwrg_residual_space(gausslet_overlap_3d, blocks.overlap_ga, blocks.overlap_aa)
        end

        final_kinetic, final_nuclear_one_body_by_center, final_one_body = @timeg "qwrg.diatomic_shell.one_body" begin
            gausslet_kinetic = _qwrg_diatomic_kinetic_matrix(
                bundles.bundle_x,
                bundles.bundle_y,
                bundles.bundle_z,
            )
            _, final_kinetic_local = _qwrg_one_body_matrices(
                gausslet_kinetic,
                blocks.kinetic_ga,
                blocks.kinetic_aa,
                residual_data.raw_to_final,
            )
            final_nuclear_one_body_by_center_local =
                resolved_nuclear_term_storage == :by_center ?
                [
                    _qwrg_one_body_matrices(
                        gausslet_nuclear,
                        blocks.nuclear_ga_by_center[index],
                        blocks.nuclear_aa_by_center[index],
                        residual_data.raw_to_final,
                    )[2] for (index, gausslet_nuclear) in pairs(
                        _qwrg_diatomic_nuclear_one_body_by_center(
                            basis,
                            bundles.bundle_x,
                            bundles.bundle_y,
                            bundles.bundle_z,
                            expansion,
                        ),
                    )
                ] :
                nothing
            final_one_body_local =
                isnothing(final_nuclear_one_body_by_center_local) ?
                _qwrg_one_body_matrices(
                    _qwrg_diatomic_one_body_matrix(
                        basis,
                        bundles.bundle_x,
                        bundles.bundle_y,
                        bundles.bundle_z,
                        expansion,
                        nuclear_charges,
                    ),
                    blocks.one_body_ga,
                    blocks.one_body_aa,
                    residual_data.raw_to_final,
                )[2] :
                _assemble_one_body_hamiltonian(
                    final_kinetic_local,
                    final_nuclear_one_body_by_center_local,
                    nuclear_charges,
                )
            (
                final_kinetic_local,
                final_nuclear_one_body_by_center_local,
                final_one_body_local,
            )
        end

        residual_centers, residual_widths = @timeg "qwrg.diatomic_shell.centers" begin
            x_gg = _qwrg_diatomic_gausslet_axis_matrix(bundles.bundle_x, bundles.bundle_y, bundles.bundle_z, :x)
            y_gg = _qwrg_diatomic_gausslet_axis_matrix(bundles.bundle_x, bundles.bundle_y, bundles.bundle_z, :y)
            z_gg = _qwrg_diatomic_gausslet_axis_matrix(bundles.bundle_x, bundles.bundle_y, bundles.bundle_z, :z)
            x_raw = [x_gg blocks.position_x_ga; transpose(blocks.position_x_ga) blocks.position_x_aa]
            y_raw = [y_gg blocks.position_y_ga; transpose(blocks.position_y_ga) blocks.position_y_aa]
            z_raw = [z_gg blocks.position_z_ga; transpose(blocks.position_z_ga) blocks.position_z_aa]
            center_data = _qwrg_residual_center_data(
                residual_data.raw_overlap,
                x_raw,
                y_raw,
                z_raw,
                residual_data.raw_to_final,
                gausslet_count,
            )
            residual_centers_local = center_data.centers
            residual_widths_local = fill(NaN, size(residual_centers_local, 1), 3)
            center_data.overlap_error <= 1.0e-8 || throw(
                ArgumentError("bond-aligned diatomic residual-center extraction requires an orthonormal residual block"),
            )
            (residual_centers_local, residual_widths_local)
        end

        interaction_matrix = @timeg "qwrg.diatomic_shell.interaction" begin
            gausslet_interaction = _qwrg_diatomic_interaction_matrix(
                bundles.bundle_x,
                bundles.bundle_y,
                bundles.bundle_z,
                expansion,
            )
            _qwrg_interaction_matrix_nearest(
                gausslet_interaction,
                gausslet_orbitals,
                residual_centers,
            )
        end

        return OrdinaryCartesianOperators3D(
            basis,
            gaussian_data,
            gausslet_backend,
            :ggt_nearest,
            expansion,
            Matrix{Float64}(residual_data.final_overlap),
            final_one_body,
            Matrix{Float64}(0.5 .* (interaction_matrix .+ transpose(interaction_matrix))),
            _qwrg_orbital_data(
                gausslet_orbitals,
                residual_centers,
                residual_widths,
            ),
            gausslet_count,
            size(residual_centers, 1),
            Matrix{Float64}(residual_data.raw_to_final),
            Matrix{Float64}(residual_centers),
            Matrix{Float64}(residual_widths),
            Float64[Float64(value) for value in nuclear_charges],
            resolved_nuclear_term_storage == :by_center ? Matrix{Float64}(final_kinetic) : nothing,
            resolved_nuclear_term_storage == :by_center ? final_nuclear_one_body_by_center : nothing,
            resolved_nuclear_term_storage,
        )
    end
end

function _ordinary_cartesian_qiu_white_operators_bond_aligned_nested_fixed_block(
    fixed_block::_NestedFixedBlock3D{<:AbstractBondAlignedOrdinaryQWBasis3D};
    nuclear_charges::AbstractVector{<:Real},
    nuclear_term_storage::Symbol,
    expansion::CoulombGaussianExpansion,
    interaction_treatment::Symbol,
    gausslet_backend::Symbol,
    timing::Bool,
    context_label::AbstractString,
)
    gausslet_backend == :numerical_reference || throw(
        ArgumentError("$(context_label) currently supports only gausslet_backend = :numerical_reference"),
    )
    interaction_treatment == :ggt_nearest || throw(
        ArgumentError("$(context_label) currently supports only interaction_treatment = :ggt_nearest"),
    )
    return _qwrg_capture_timeg_report("qwrg.bond_aligned_nested_fixed.total", timing) do
        basis = fixed_block.parent_basis
        resolved_nuclear_term_storage = _resolved_nuclear_term_storage(
            _validate_nuclear_term_storage(nuclear_term_storage),
            fixed_block,
        )
        length(nuclear_charges) == length(basis.nuclei) || throw(
            ArgumentError("$(context_label) requires one nuclear charge per nucleus"),
        )

        bundles = _qwrg_bond_aligned_axis_bundles(
            basis,
            expansion;
            gausslet_backend = gausslet_backend,
        )

        contraction = fixed_block.coefficient_matrix
        kinetic_one_body = Matrix{Float64}(fixed_block.kinetic)
        nuclear_one_body_by_center =
            resolved_nuclear_term_storage == :by_center ?
            [
                _qwrg_contract_parent_symmetric_matrix(
                    contraction,
                    matrix,
                ) for matrix in _qwrg_diatomic_nuclear_one_body_by_center(
                    basis,
                    bundles.bundle_x,
                    bundles.bundle_y,
                    bundles.bundle_z,
                    expansion,
                )
            ] :
            nothing
        one_body_hamiltonian =
            isnothing(nuclear_one_body_by_center) ?
            _qwrg_diatomic_one_body_matrix(
                basis,
                bundles.bundle_x,
                bundles.bundle_y,
                bundles.bundle_z,
                expansion,
                nuclear_charges,
            ) :
            _assemble_one_body_hamiltonian(
                kinetic_one_body,
                nuclear_one_body_by_center,
                nuclear_charges,
            )
        if isnothing(nuclear_one_body_by_center)
            one_body_hamiltonian = Matrix{Float64}(transpose(contraction) * one_body_hamiltonian * contraction)
            one_body_hamiltonian = 0.5 .* (one_body_hamiltonian .+ transpose(one_body_hamiltonian))
        end
        interaction_matrix = _qwrg_fixed_block_interaction_matrix(fixed_block, expansion)

        fixed_count = size(fixed_block.overlap, 1)
        zero_residual_centers = zeros(Float64, 0, 3)
        zero_residual_widths = zeros(Float64, 0, 3)
        return OrdinaryCartesianOperators3D(
            fixed_block,
            nothing,
            gausslet_backend,
            interaction_treatment,
            expansion,
            Matrix{Float64}(fixed_block.overlap),
            one_body_hamiltonian,
            Matrix{Float64}(0.5 .* (interaction_matrix .+ transpose(interaction_matrix))),
            _qwrg_orbital_data(
                fixed_block.fixed_centers,
                zero_residual_centers,
                zero_residual_widths;
                fixed_kind = :nested_fixed,
                fixed_label_prefix = "nf",
            ),
            fixed_count,
            0,
            Matrix{Float64}(I, fixed_count, fixed_count),
            zero_residual_centers,
            zero_residual_widths,
            Float64[Float64(value) for value in nuclear_charges],
            resolved_nuclear_term_storage == :by_center ? kinetic_one_body : nothing,
            resolved_nuclear_term_storage == :by_center ? nuclear_one_body_by_center : nothing,
            resolved_nuclear_term_storage,
        )
    end
end

function _ordinary_cartesian_qiu_white_operators_diatomic_fixed_block(
    fixed_block::_NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D};
    nuclear_charges::AbstractVector{<:Real},
    nuclear_term_storage::Symbol,
    expansion::CoulombGaussianExpansion,
    interaction_treatment::Symbol,
    gausslet_backend::Symbol,
    timing::Bool,
)
    return _ordinary_cartesian_qiu_white_operators_bond_aligned_nested_fixed_block(
        fixed_block;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        timing = timing,
        context_label = "bond-aligned diatomic nested ordinary_cartesian_qiu_white_operators",
    )
end

function _ordinary_cartesian_qiu_white_operators_homonuclear_chain_fixed_block(
    fixed_block::_NestedFixedBlock3D{<:BondAlignedHomonuclearChainQWBasis3D};
    nuclear_charges::AbstractVector{<:Real},
    nuclear_term_storage::Symbol,
    expansion::CoulombGaussianExpansion,
    interaction_treatment::Symbol,
    gausslet_backend::Symbol,
    timing::Bool,
)
    return _ordinary_cartesian_qiu_white_operators_bond_aligned_nested_fixed_block(
        fixed_block;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        timing = timing,
        context_label = "experimental bond-aligned homonuclear chain nested ordinary_cartesian_qiu_white_operators",
    )
end

function _ordinary_cartesian_qiu_white_operators_square_lattice_fixed_block(
    fixed_block::_NestedFixedBlock3D{<:AxisAlignedHomonuclearSquareLatticeQWBasis3D};
    nuclear_charges::AbstractVector{<:Real},
    nuclear_term_storage::Symbol,
    expansion::CoulombGaussianExpansion,
    interaction_treatment::Symbol,
    gausslet_backend::Symbol,
    timing::Bool,
)
    return _ordinary_cartesian_qiu_white_operators_bond_aligned_nested_fixed_block(
        fixed_block;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        timing = timing,
        context_label = "experimental axis-aligned homonuclear square-lattice nested ordinary_cartesian_qiu_white_operators",
    )
end

function _ordinary_cartesian_qiu_white_operators_nested_diatomic_shell_3d(
    fixed_block::_NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D},
    gaussian_data::Union{
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    };
    nuclear_charges::AbstractVector{<:Real},
    nuclear_term_storage::Symbol,
    expansion::CoulombGaussianExpansion,
    gausslet_backend::Symbol,
    timing::Bool,
)
    gausslet_backend == :numerical_reference || throw(
        ArgumentError("bond-aligned diatomic nested molecular QW path currently supports only gausslet_backend = :numerical_reference"),
    )

    basis = fixed_block.parent_basis
    length(nuclear_charges) == length(basis.nuclei) || throw(
        ArgumentError("bond-aligned diatomic nested molecular QW path requires one nuclear charge per nucleus"),
    )
    _qwrg_same_nuclei(gaussian_data.nuclei, basis.nuclei) || throw(
        ArgumentError("bond-aligned diatomic molecular supplement nuclei must match the bond-aligned basis nuclei"),
    )
    resolved_nuclear_term_storage = _resolved_nuclear_term_storage(
        _validate_nuclear_term_storage(nuclear_term_storage),
        fixed_block,
    )

    return _qwrg_capture_timeg_report("qwrg.nested_diatomic_shell.total", timing) do
        bundles = @timeg "qwrg.nested_diatomic_shell.parent_bundles" begin
            _qwrg_bond_aligned_axis_bundles(
                basis,
                expansion;
                gausslet_backend = gausslet_backend,
            )
        end

        supplement3d = _bond_aligned_diatomic_cartesian_shell_supplement_3d(gaussian_data)

        blocks = @timeg "qwrg.nested_diatomic_shell.raw_blocks" begin
            _qwrg_diatomic_cartesian_shell_blocks_3d(
                bundles,
                supplement3d,
                basis,
                expansion,
                nuclear_charges,
            )
        end

        contraction = fixed_block.coefficient_matrix
        fixed_count = size(fixed_block.overlap, 1)

        residual_data = @timeg "qwrg.nested_diatomic_shell.residual_space" begin
            overlap_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.overlap_ga)
            _qwrg_residual_space(fixed_block.overlap, overlap_fg, blocks.overlap_aa)
        end

        final_kinetic, nuclear_one_body_by_center, final_one_body = @timeg "qwrg.nested_diatomic_shell.one_body" begin
            fixed_kinetic = Matrix{Float64}(fixed_block.kinetic)
            nuclear_one_body_by_center_local =
                resolved_nuclear_term_storage == :by_center ?
                let
                    parent_nuclear = _qwrg_diatomic_nuclear_one_body_by_center(
                        basis,
                        bundles.bundle_x,
                        bundles.bundle_y,
                        bundles.bundle_z,
                        expansion,
                    )
                    [
                        _qwrg_one_body_matrices(
                            _qwrg_contract_parent_symmetric_matrix(
                                contraction,
                                parent_nuclear[index],
                            ),
                            _qwrg_contract_parent_ga_matrix(contraction, blocks.nuclear_ga_by_center[index]),
                            blocks.nuclear_aa_by_center[index],
                            residual_data.raw_to_final,
                        )[2] for index in eachindex(parent_nuclear)
                    ]
                end :
                nothing
            final_kinetic_local = _qwrg_one_body_matrices(
                fixed_kinetic,
                _qwrg_contract_parent_ga_matrix(contraction, blocks.kinetic_ga),
                blocks.kinetic_aa,
                residual_data.raw_to_final,
            )[2]
            final_one_body_local =
                isnothing(nuclear_one_body_by_center_local) ?
                let
                    parent_one_body = _qwrg_diatomic_one_body_matrix(
                        basis,
                        bundles.bundle_x,
                        bundles.bundle_y,
                        bundles.bundle_z,
                        expansion,
                        nuclear_charges,
                    )
                    fixed_one_body = Matrix{Float64}(transpose(contraction) * parent_one_body * contraction)
                    fixed_one_body = 0.5 .* (fixed_one_body .+ transpose(fixed_one_body))
                    one_body_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.one_body_ga)
                    _qwrg_one_body_matrices(
                        fixed_one_body,
                        one_body_fg,
                        blocks.one_body_aa,
                        residual_data.raw_to_final,
                    )[2]
                end :
                _assemble_one_body_hamiltonian(
                    final_kinetic_local,
                    nuclear_one_body_by_center_local,
                    nuclear_charges,
                )
            (final_kinetic_local, nuclear_one_body_by_center_local, final_one_body_local)
        end

        residual_centers, residual_widths = @timeg "qwrg.nested_diatomic_shell.centers" begin
            x_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.position_x_ga)
            y_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.position_y_ga)
            z_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.position_z_ga)
            x_raw = [Matrix{Float64}(fixed_block.position_x) x_fg; transpose(x_fg) blocks.position_x_aa]
            y_raw = [Matrix{Float64}(fixed_block.position_y) y_fg; transpose(y_fg) blocks.position_y_aa]
            z_raw = [Matrix{Float64}(fixed_block.position_z) z_fg; transpose(z_fg) blocks.position_z_aa]
            center_data = _qwrg_residual_center_data(
                residual_data.raw_overlap,
                x_raw,
                y_raw,
                z_raw,
                residual_data.raw_to_final,
                fixed_count,
            )
            residual_centers_local = center_data.centers
            residual_widths_local = fill(NaN, size(residual_centers_local, 1), 3)
            center_data.overlap_error <= 1.0e-8 || throw(
                ArgumentError("bond-aligned diatomic nested residual-center extraction requires an orthonormal residual block"),
            )
            (residual_centers_local, residual_widths_local)
        end

        interaction_matrix = @timeg "qwrg.nested_diatomic_shell.interaction" begin
            fixed_interaction = _qwrg_fixed_block_interaction_matrix(fixed_block, expansion)
            _qwrg_interaction_matrix_nearest(
                fixed_interaction,
                fixed_block.fixed_centers,
                residual_centers,
            )
        end

        return OrdinaryCartesianOperators3D(
            fixed_block,
            gaussian_data,
            gausslet_backend,
            :ggt_nearest,
            expansion,
            Matrix{Float64}(residual_data.final_overlap),
            final_one_body,
            Matrix{Float64}(0.5 .* (interaction_matrix .+ transpose(interaction_matrix))),
            _qwrg_orbital_data(
                fixed_block.fixed_centers,
                residual_centers,
                residual_widths;
                fixed_kind = :nested_fixed,
                fixed_label_prefix = "nf",
            ),
            fixed_count,
            size(residual_centers, 1),
            Matrix{Float64}(residual_data.raw_to_final),
            Matrix{Float64}(residual_centers),
            Matrix{Float64}(residual_widths),
            Float64[Float64(value) for value in nuclear_charges],
            resolved_nuclear_term_storage == :by_center ? Matrix{Float64}(final_kinetic) : nothing,
            resolved_nuclear_term_storage == :by_center ? nuclear_one_body_by_center : nothing,
            resolved_nuclear_term_storage,
        )
    end
end

function _ordinary_cartesian_qiu_white_operators_atomic_shell_3d(
    basis::MappedUniformBasis,
    gaussian_data::LegacyAtomicGaussianSupplement;
    expansion::CoulombGaussianExpansion,
    Z::Real,
    interaction_treatment::Symbol,
    gausslet_backend::Symbol,
    residual_keep_policy::Symbol = :near_null_only,
    timing::Bool,
)
    residual_keep_tol = _qwrg_atomic_residual_keep_tol()
    residual_accept_tol = _qwrg_atomic_residual_accept_tol()
    return _qwrg_capture_timeg_report("qwrg.atomic_shell.total", timing) do
        gausslet_bundle = @timeg "qwrg.atomic_shell.shared_bundle" begin
            _mapped_ordinary_gausslet_1d_bundle(
                basis,
                exponents = expansion.exponents,
                center = 0.0,
                backend = gausslet_backend,
            )
        end

        supplement3d = _atomic_cartesian_shell_supplement_3d(gaussian_data)
        gg_blocks = _qwrg_gausslet_1d_blocks(gausslet_bundle)

        blocks = @timeg "qwrg.atomic_shell.raw_blocks" begin
            _qwrg_atomic_cartesian_blocks_3d(
                gausslet_bundle,
                supplement3d,
                expansion,
            )
        end

        gausslet_orbitals = _mapped_cartesian_orbitals(gausslet_bundle.pgdg_intermediate.centers)
        gausslet_count = length(gausslet_orbitals)

        residual_data = @timeg "qwrg.atomic_shell.residual_space" begin
            gausslet_overlap_3d = zeros(Float64, gausslet_count, gausslet_count)
            _qwrg_fill_product_matrix!(
                gausslet_overlap_3d,
                gg_blocks.overlap_gg,
                gg_blocks.overlap_gg,
                gg_blocks.overlap_gg,
            )
            _qwrg_residual_space(
                gausslet_overlap_3d,
                blocks.overlap_ga,
                blocks.overlap_aa;
                keep_policy = residual_keep_policy,
                keep_abs_tol = residual_keep_tol,
                accept_tol = residual_accept_tol,
            )
        end
        timing && _qwrg_print_basis_counts(stdout, gausslet_count, residual_data.raw_to_final)

        final_one_body = @timeg "qwrg.atomic_shell.one_body" begin
            gausslet_one_body = _qwrg_gausslet_one_body_matrix(gg_blocks, expansion; Z = Z)
            one_body_ga = Matrix{Float64}(blocks.kinetic_ga)
            for term in eachindex(expansion.coefficients)
                one_body_ga .-= Float64(Z) * expansion.coefficients[term] .* blocks.factor_ga[term]
            end
            one_body_aa = _qwrg_atomic_cartesian_one_body_aa(blocks, expansion; Z = Z)
            _qwrg_one_body_matrices(
                gausslet_one_body,
                one_body_ga,
                one_body_aa,
                residual_data.raw_to_final,
            )[2]
        end

        residual_centers, residual_widths = @timeg "qwrg.atomic_shell.centers" begin
            x_gg = _qwrg_gausslet_axis_matrix(gg_blocks, :x)
            y_gg = _qwrg_gausslet_axis_matrix(gg_blocks, :y)
            z_gg = _qwrg_gausslet_axis_matrix(gg_blocks, :z)
            x_raw = [x_gg blocks.position_x_ga; transpose(blocks.position_x_ga) blocks.position_x_aa]
            y_raw = [y_gg blocks.position_y_ga; transpose(blocks.position_y_ga) blocks.position_y_aa]
            z_raw = [z_gg blocks.position_z_ga; transpose(blocks.position_z_ga) blocks.position_z_aa]
            center_data = _qwrg_residual_center_data(
                residual_data.raw_overlap,
                x_raw,
                y_raw,
                z_raw,
                residual_data.raw_to_final,
                gausslet_count,
            )
            residual_centers_local = center_data.centers
            residual_widths_local = fill(NaN, size(residual_centers_local, 1), 3)
            center_data.overlap_error <= residual_accept_tol || throw(
                ArgumentError("Qiu-White residual-center extraction requires an orthonormal residual block"),
            )

            if interaction_treatment == :mwg
                gg_x2 = (
                    overlap_gg = gg_blocks.overlap_gg,
                    position_gg = gg_blocks.position_gg,
                    x2_gg = Matrix{Float64}(_x2_matrix(gausslet_bundle.basis)),
                )
                x2_x_gg = _qwrg_gausslet_axis_matrix(gg_x2, :x; squared = true)
                x2_y_gg = _qwrg_gausslet_axis_matrix(gg_x2, :y; squared = true)
                x2_z_gg = _qwrg_gausslet_axis_matrix(gg_x2, :z; squared = true)
                x2_raw = [x2_x_gg blocks.x2_x_ga; transpose(blocks.x2_x_ga) blocks.x2_x_aa]
                y2_raw = [x2_y_gg blocks.x2_y_ga; transpose(blocks.x2_y_ga) blocks.x2_y_aa]
                z2_raw = [x2_z_gg blocks.x2_z_ga; transpose(blocks.x2_z_ga) blocks.x2_z_aa]
                moment_data = _qwrg_residual_moment_data(
                    residual_data.raw_overlap,
                    x_raw,
                    x2_raw,
                    y_raw,
                    y2_raw,
                    z_raw,
                    z2_raw,
                    center_data,
                )
                residual_centers_local = moment_data.centers
                residual_widths_local = moment_data.widths
            end
            (residual_centers_local, residual_widths_local)
        end

        interaction_matrix = @timeg "qwrg.atomic_shell.interaction" begin
            gausslet_interaction = _qwrg_gausslet_interaction_matrix(gg_blocks, expansion)
            if interaction_treatment == :ggt_nearest
                _qwrg_interaction_matrix_nearest(
                    gausslet_interaction,
                    gausslet_orbitals,
                    residual_centers,
                )
            elseif interaction_treatment == :mwg
                _qwrg_interaction_matrix_mwg(
                    gausslet_bundle,
                    gausslet_interaction,
                    expansion,
                    residual_centers,
                    residual_widths,
                )
            else
                throw(ArgumentError("Qiu-White interaction_treatment must be :ggt_nearest or :mwg"))
            end
        end

        return OrdinaryCartesianOperators3D(
            basis,
            gaussian_data,
            gausslet_backend,
            interaction_treatment,
            expansion,
            Matrix{Float64}(residual_data.final_overlap),
            final_one_body,
            Matrix{Float64}(0.5 .* (interaction_matrix .+ transpose(interaction_matrix))),
            _qwrg_orbital_data(
                gausslet_orbitals,
                residual_centers,
                residual_widths,
            ),
            gausslet_count,
            size(residual_centers, 1),
            Matrix{Float64}(residual_data.raw_to_final),
            Matrix{Float64}(residual_centers),
            Matrix{Float64}(residual_widths),
        )
    end
end

function _ordinary_cartesian_qiu_white_operators_nested_atomic_shell_3d(
    fixed_block::_NestedFixedBlock3D,
    gaussian_data::LegacyAtomicGaussianSupplement;
    expansion::CoulombGaussianExpansion,
    Z::Real,
    gausslet_backend::Symbol,
    residual_keep_policy::Symbol = :near_null_only,
    timing::Bool,
)
    residual_keep_tol = _qwrg_atomic_residual_keep_tol()
    residual_accept_tol = _qwrg_atomic_residual_accept_tol()
    return _qwrg_capture_timeg_report("qwrg.nested_atomic_shell.total", timing) do
        gausslet_bundle = @timeg "qwrg.nested_atomic_shell.parent_bundle" begin
            _mapped_ordinary_gausslet_1d_bundle(
                fixed_block.parent_basis;
                exponents = expansion.exponents,
                center = 0.0,
                backend = gausslet_backend,
            )
        end

        supplement3d = _atomic_cartesian_shell_supplement_3d(gaussian_data)

        blocks = @timeg "qwrg.nested_atomic_shell.raw_blocks" begin
            _qwrg_atomic_cartesian_blocks_3d(
                gausslet_bundle,
                supplement3d,
                expansion,
            )
        end

        contraction = fixed_block.coefficient_matrix
        fixed_count = size(fixed_block.overlap, 1)

        residual_data = @timeg "qwrg.nested_atomic_shell.residual_space" begin
            overlap_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.overlap_ga)
            _qwrg_residual_space(
                fixed_block.overlap,
                overlap_fg,
                blocks.overlap_aa;
                keep_policy = residual_keep_policy,
                keep_abs_tol = residual_keep_tol,
                accept_tol = residual_accept_tol,
            )
        end
        timing && _qwrg_print_basis_counts(stdout, "fixed_count", fixed_count, residual_data.raw_to_final)

        final_one_body = @timeg "qwrg.nested_atomic_shell.one_body" begin
            kinetic_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.kinetic_ga)
            factor_fg = _qwrg_contract_parent_ga_terms(contraction, blocks.factor_ga)
            one_body_fg = Matrix{Float64}(kinetic_fg)
            for term in eachindex(expansion.coefficients)
                one_body_fg .-= Float64(Z) * expansion.coefficients[term] .* factor_fg[term]
            end
            one_body_aa = _qwrg_atomic_cartesian_one_body_aa(blocks, expansion; Z = Z)
            fixed_one_body = _qwrg_fixed_block_one_body_matrix(fixed_block, expansion; Z = Z)
            _qwrg_one_body_matrices(
                fixed_one_body,
                one_body_fg,
                one_body_aa,
                residual_data.raw_to_final,
            )[2]
        end

        residual_centers, residual_widths = @timeg "qwrg.nested_atomic_shell.centers" begin
            x_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.position_x_ga)
            y_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.position_y_ga)
            z_fg = _qwrg_contract_parent_ga_matrix(contraction, blocks.position_z_ga)
            x_raw = [Matrix{Float64}(fixed_block.position_x) x_fg; transpose(x_fg) blocks.position_x_aa]
            y_raw = [Matrix{Float64}(fixed_block.position_y) y_fg; transpose(y_fg) blocks.position_y_aa]
            z_raw = [Matrix{Float64}(fixed_block.position_z) z_fg; transpose(z_fg) blocks.position_z_aa]
            center_data = _qwrg_residual_center_data(
                residual_data.raw_overlap,
                x_raw,
                y_raw,
                z_raw,
                residual_data.raw_to_final,
                fixed_count,
            )
            residual_centers_local = center_data.centers
            residual_widths_local = fill(NaN, size(residual_centers_local, 1), 3)
            center_data.overlap_error <= residual_accept_tol || throw(
                ArgumentError("nested QW-PGDG residual-center extraction requires an orthonormal residual block"),
            )
            (residual_centers_local, residual_widths_local)
        end

        interaction_matrix = @timeg "qwrg.nested_atomic_shell.interaction" begin
            fixed_interaction = _qwrg_fixed_block_interaction_matrix(fixed_block, expansion)
            _qwrg_interaction_matrix_nearest(
                fixed_interaction,
                fixed_block.fixed_centers,
                residual_centers,
            )
        end

        return OrdinaryCartesianOperators3D(
            fixed_block,
            gaussian_data,
            gausslet_backend,
            :ggt_nearest,
            expansion,
            Matrix{Float64}(residual_data.final_overlap),
            final_one_body,
            Matrix{Float64}(0.5 .* (interaction_matrix .+ transpose(interaction_matrix))),
            _qwrg_orbital_data(
                fixed_block.fixed_centers,
                residual_centers,
                residual_widths;
                fixed_kind = :nested_fixed,
                fixed_label_prefix = "nf",
            ),
            fixed_count,
            size(residual_centers, 1),
            Matrix{Float64}(residual_data.raw_to_final),
            Matrix{Float64}(residual_centers),
            Matrix{Float64}(residual_widths),
        )
    end
end

"""
    ordinary_cartesian_qiu_white_operators(
        basis::MappedUniformBasis,
        gaussian_data::LegacyAtomicGaussianSupplement;
        expansion = coulomb_gaussian_expansion(doacc = false),
        Z = 2.0,
        interaction_treatment = :mwg,
        gausslet_backend = :numerical_reference,
        residual_keep_policy = :near_null_only,
        timing = false,
    )

Build the paper-faithful Qiu-White residual-Gaussian ordinary Cartesian
reference Hamiltonian.

This path is intentionally separate from the current COMX/localized hybrid
route. It:

- keeps the full 3D gausslet product basis fixed
- orthogonalizes the added 3D Gaussian supplement against that full 3D space
- builds the one-body matrices exactly in the raw gausslet-plus-GTO space
- keeps the two-electron interaction in the same two-index integral-diagonal
  approximation (IDA) representation used for the gausslet channel

Allowed `interaction_treatment` values are:

- `:ggt_nearest`
- `:mwg`

Set `timing = true` to print a TimeG phase tree for debugging the reference
implementation. This is intentionally narrow and does not add timing noise to
the broader library.

This is a reference path for validating the Qiu-White formulation, not yet a
claim that the ordinary branch is solver-ready.

This path now uses the explicit atomic-centered 3D Cartesian shell route for
all active atomic supplement content up to `lmax <= 2`, including pure `s`
shells.

`nested_profile` selects the shell/working-box geometry contract. It is
independent of supplement residual retention.

`residual_keep_policy = :near_null_only` is now the canonical atomic
supplement keep contract. It keeps orthogonalized residual directions with
residual-overlap eigenvalue `> 1e-7`.

`residual_keep_policy = :legacy_profile` remains accepted only as a
compatibility alias for `:near_null_only`.

After residual directions are selected, the retained residual block is
explicitly stabilized back to overlap-orthonormality before any downstream
center/width extraction or raw-to-final transforms are used. The atomic
residual-center acceptance gate on that stabilized block is `1e-7`.
"""
function ordinary_cartesian_qiu_white_operators(
    basis::MappedUniformBasis,
    gaussian_data::LegacyAtomicGaussianSupplement;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
    interaction_treatment::Symbol = :mwg,
    gausslet_backend::Symbol = :numerical_reference,
    residual_keep_policy::Symbol = :near_null_only,
    timing::Bool = false,
    )
    gausslet_backend == :numerical_reference || throw(
        ArgumentError("ordinary_cartesian_qiu_white_operators currently supports only gausslet_backend = :numerical_reference"),
    )
    residual_keep_policy_value = _qwrg_atomic_residual_keep_policy(residual_keep_policy)
    return _ordinary_cartesian_qiu_white_operators_atomic_shell_3d(
        basis,
        gaussian_data;
        expansion = expansion,
        Z = Z,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        residual_keep_policy = residual_keep_policy_value,
        timing = timing,
    )
end

"""
    ordinary_cartesian_qiu_white_operators(
        fixed_block::_NestedFixedBlock3D,
        gaussian_data::LegacyAtomicGaussianSupplement;
        expansion = coulomb_gaussian_expansion(doacc = false),
        Z = 2.0,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        residual_keep_policy = :near_null_only,
        timing = false,
    )

Build the first nested fixed-block QW-PGDG consumer path.

This overload reuses the stabilized parent PGDG raw fixed-to-Gaussian blocks,
contracts them through the supplied shell-level fixed map, and then runs the
same downstream residual-space / one-body / nearest-GGT algebra as the
unnested QW-PGDG route.

This first adapter is intentionally narrow:

- it consumes an already-assembled nonseparable 3D fixed packet
- it keeps the shell packet as the fixed-fixed block directly
- it supports only `interaction_treatment = :ggt_nearest`

It now uses the explicit atomic-centered 3D Cartesian shell route for all
active atomic supplement content up to `lmax <= 2`, including pure `s`
shells.

`nested_profile` selects the shell/working-box geometry contract. It is
independent of supplement residual retention.

`residual_keep_policy = :near_null_only` is now the canonical atomic
supplement keep contract. It keeps orthogonalized residual directions with
residual-overlap eigenvalue `> 1e-7`.

`residual_keep_policy = :legacy_profile` remains accepted only as a
compatibility alias for `:near_null_only`.

After residual directions are selected, the retained residual block is
explicitly stabilized back to overlap-orthonormality before any downstream
center/width extraction or raw-to-final transforms are used. The atomic
residual-center acceptance gate on that stabilized block is `1e-7`.
"""
function ordinary_cartesian_qiu_white_operators(
    fixed_block::_NestedFixedBlock3D,
    gaussian_data::LegacyAtomicGaussianSupplement;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    residual_keep_policy::Symbol = :near_null_only,
    timing::Bool = false,
)
    gausslet_backend == :numerical_reference || throw(
        ArgumentError("nested ordinary_cartesian_qiu_white_operators currently supports only gausslet_backend = :numerical_reference"),
    )
    interaction_treatment == :ggt_nearest || throw(
        ArgumentError("nested ordinary_cartesian_qiu_white_operators currently supports only interaction_treatment = :ggt_nearest"),
    )
    residual_keep_policy_value = _qwrg_atomic_residual_keep_policy(residual_keep_policy)
    return _ordinary_cartesian_qiu_white_operators_nested_atomic_shell_3d(
        fixed_block,
        gaussian_data;
        expansion = expansion,
        Z = Z,
        gausslet_backend = gausslet_backend,
        residual_keep_policy = residual_keep_policy_value,
        timing = timing,
    )
end

"""
    ordinary_cartesian_qiu_white_operators(
        fixed_block::_NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D};
        nuclear_charges = fill(1.0, length(fixed_block.parent_basis.nuclei)),
        expansion = coulomb_gaussian_expansion(doacc = false),
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        timing = false,
    )

Build the first bond-aligned diatomic nested fixed-block ordinary QW path.

This is intentionally narrower than the atomic hybrid route:

- it consumes the already-assembled bond-aligned diatomic `_NestedFixedBlock3D`
- it keeps the residual-Gaussian sector empty
- it supports only `interaction_treatment = :ggt_nearest`

The point of this first pass is to validate the diatomic fixed-block geometry
and transferred packet cleanly before molecular Gaussian completion is added.
"""
function ordinary_cartesian_qiu_white_operators(
    fixed_block::_NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D};
    nuclear_charges::AbstractVector{<:Real} = fill(1.0, length(fixed_block.parent_basis.nuclei)),
    nuclear_term_storage::Symbol = :auto,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    timing::Bool = false,
)
    return _ordinary_cartesian_qiu_white_operators_diatomic_fixed_block(
        fixed_block;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        timing = timing,
    )
end

"""
    ordinary_cartesian_qiu_white_operators(
        fixed_block::_NestedFixedBlock3D{<:BondAlignedHomonuclearChainQWBasis3D};
        nuclear_charges = fill(1.0, length(fixed_block.parent_basis.nuclei)),
        expansion = coulomb_gaussian_expansion(doacc = false),
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        timing = false,
    )

Build the first experimental bond-aligned homonuclear-chain nested fixed-block
ordinary QW path.

This path is intentionally narrow:

- homonuclear chains only
- nested fixed-block route only
- no supplement route
- no residual-Gaussian sector
- no claim that the odd-chain split policy is settled globally
"""
function ordinary_cartesian_qiu_white_operators(
    fixed_block::_NestedFixedBlock3D{<:BondAlignedHomonuclearChainQWBasis3D};
    nuclear_charges::AbstractVector{<:Real} = fill(1.0, length(fixed_block.parent_basis.nuclei)),
    nuclear_term_storage::Symbol = :auto,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    timing::Bool = false,
)
    return _ordinary_cartesian_qiu_white_operators_homonuclear_chain_fixed_block(
        fixed_block;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        timing = timing,
    )
end

"""
    ordinary_cartesian_qiu_white_operators(
        fixed_block::_NestedFixedBlock3D{<:AxisAlignedHomonuclearSquareLatticeQWBasis3D};
        nuclear_charges = fill(1.0, length(fixed_block.parent_basis.nuclei)),
        expansion = coulomb_gaussian_expansion(doacc = false),
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        timing = false,
    )

Build the first experimental nested fixed-block ordinary-QW square-lattice
Hamiltonian.

This path is intentionally narrow:

- homonuclear square lattices only
- nested fixed-block route only
- no supplement route
- no residual-Gaussian sector
- no claim that the planar split policy is settled globally
"""
function ordinary_cartesian_qiu_white_operators(
    fixed_block::_NestedFixedBlock3D{<:AxisAlignedHomonuclearSquareLatticeQWBasis3D};
    nuclear_charges::AbstractVector{<:Real} = fill(1.0, length(fixed_block.parent_basis.nuclei)),
    nuclear_term_storage::Symbol = :auto,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    timing::Bool = false,
)
    return _ordinary_cartesian_qiu_white_operators_square_lattice_fixed_block(
        fixed_block;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        timing = timing,
    )
end


"""
    ordinary_cartesian_qiu_white_operators(
        basis::BondAlignedDiatomicQWBasis3D,
        gaussian_data::LegacyBondAlignedDiatomicGaussianSupplement;
        nuclear_charges = fill(1.0, length(basis.nuclei)),
        expansion = coulomb_gaussian_expansion(doacc = false),
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        timing = false,
    )

Build the first bond-aligned diatomic ordinary QW route with a true molecular
supplement and residual-Gaussian completion.

This first molecular supplement pass is intentionally narrow:

- one bond-aligned diatomic basis
- one explicit two-center molecular shell supplement built from a named atomic
  basis
- only `interaction_treatment = :ggt_nearest`
"""
function ordinary_cartesian_qiu_white_operators(
    basis::BondAlignedDiatomicQWBasis3D,
    gaussian_data::Union{
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    };
    nuclear_charges::AbstractVector{<:Real} = fill(1.0, length(basis.nuclei)),
    nuclear_term_storage::Symbol = :auto,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    timing::Bool = false,
)
    return _ordinary_cartesian_qiu_white_operators_diatomic_shell_3d(
        basis,
        gaussian_data;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        interaction_treatment = interaction_treatment,
        gausslet_backend = gausslet_backend,
        timing = timing,
    )
end

"""
    ordinary_cartesian_qiu_white_operators(
        fixed_block::_NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D},
        gaussian_data::Union{
            LegacyBondAlignedDiatomicGaussianSupplement,
            LegacyBondAlignedHeteronuclearGaussianSupplement,
        };
        nuclear_charges = fill(1.0, length(fixed_block.parent_basis.nuclei)),
        expansion = coulomb_gaussian_expansion(doacc = false),
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :numerical_reference,
        timing = false,
    )

Build the first bond-aligned diatomic nested fixed-block QW route with a true
molecular supplement and residual-Gaussian completion.
"""
function ordinary_cartesian_qiu_white_operators(
    fixed_block::_NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D},
    gaussian_data::Union{
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    };
    nuclear_charges::AbstractVector{<:Real} = fill(1.0, length(fixed_block.parent_basis.nuclei)),
    nuclear_term_storage::Symbol = :auto,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    interaction_treatment::Symbol = :ggt_nearest,
    gausslet_backend::Symbol = :numerical_reference,
    timing::Bool = false,
)
    interaction_treatment == :ggt_nearest || throw(
        ArgumentError("bond-aligned diatomic nested molecular QW path currently supports only interaction_treatment = :ggt_nearest"),
    )
    return _ordinary_cartesian_qiu_white_operators_nested_diatomic_shell_3d(
        fixed_block,
        gaussian_data;
        nuclear_charges = nuclear_charges,
        nuclear_term_storage = nuclear_term_storage,
        expansion = expansion,
        gausslet_backend = gausslet_backend,
        timing = timing,
    )
end

"""
    ordinary_cartesian_product_operators(...)

Compatibility-preserving clearer name for the direct-product Cartesian routes
within the ordinary Cartesian operator family. This includes both pure
Cartesian product routes and direct-product residual-Gaussian hybrid routes.
"""
function ordinary_cartesian_product_operators(
    basis::AbstractBondAlignedOrdinaryQWBasis3D;
    kwargs...,
)
    return ordinary_cartesian_qiu_white_operators(basis; kwargs...)
end

function ordinary_cartesian_product_operators(
    basis::MappedUniformBasis,
    gaussian_data::LegacyAtomicGaussianSupplement;
    kwargs...,
)
    return ordinary_cartesian_qiu_white_operators(basis, gaussian_data; kwargs...)
end

function ordinary_cartesian_product_operators(
    basis::BondAlignedDiatomicQWBasis3D,
    gaussian_data::Union{
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    };
    kwargs...,
)
    return ordinary_cartesian_qiu_white_operators(basis, gaussian_data; kwargs...)
end

"""
    nested_cartesian_operators(...)

Compatibility-preserving clearer name for the nested fixed-block Cartesian
routes within the ordinary Cartesian operator family. This includes both pure
nested fixed-block routes and nested residual-Gaussian hybrid routes.
"""
function nested_cartesian_operators(
    fixed_block::Union{
        _NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D},
        _NestedFixedBlock3D{<:BondAlignedHomonuclearChainQWBasis3D},
        _NestedFixedBlock3D{<:AxisAlignedHomonuclearSquareLatticeQWBasis3D},
    };
    kwargs...,
)
    return ordinary_cartesian_qiu_white_operators(fixed_block; kwargs...)
end

function nested_cartesian_operators(
    fixed_block::_NestedFixedBlock3D,
    gaussian_data::LegacyAtomicGaussianSupplement;
    kwargs...,
)
    return ordinary_cartesian_qiu_white_operators(fixed_block, gaussian_data; kwargs...)
end

function nested_cartesian_operators(
    fixed_block::_NestedFixedBlock3D{<:BondAlignedDiatomicQWBasis3D},
    gaussian_data::Union{
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    };
    kwargs...,
)
    return ordinary_cartesian_qiu_white_operators(fixed_block, gaussian_data; kwargs...)
end
