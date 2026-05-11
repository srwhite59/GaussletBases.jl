function _gto_working_representation(working)
    representation =
        working isa CartesianBasisRepresentation3D ? working :
        working isa MappedUniformBasis ? _cartesian_direct_product_representation(working) :
        basis_representation(working)
    representation isa CartesianBasisRepresentation3D || throw(
        ArgumentError(
            "gto_overlap_matrix requires a Cartesian working-basis representation; got $(typeof(representation)) from $(typeof(working))",
        ),
    )
    (
        representation.metadata.parent_kind == :cartesian_product_basis ||
        _cartesian_supports_exact_hybrid_overlap(representation)
    ) || throw(
        ArgumentError(
            "gto_overlap_matrix supports explicit Cartesian product, nested fixed-block, and exact hybrid residual-Gaussian Cartesian working bases; got parent_kind :$(representation.metadata.parent_kind)",
        ),
    )
    return representation
end

function _gto_probe_representation(probes)
    representation =
        probes isa CartesianGaussianShellSupplementRepresentation3D ? probes :
        basis_representation(probes)
    representation isa CartesianGaussianShellSupplementRepresentation3D || throw(
        ArgumentError(
            "gto_overlap_matrix requires a legacy Gaussian probe family or CartesianGaussianShellSupplementRepresentation3D; got $(typeof(representation)) from $(typeof(probes))",
        ),
    )
    return representation
end

function _gto_block_indices(
    final_dimension::Int,
    block_indices::Nothing,
)
    return nothing
end

function _gto_block_indices(
    final_dimension::Int,
    block_indices::AbstractVector{<:Integer},
)
    indices = Int[index for index in block_indices]
    all(index -> 1 <= index <= final_dimension, indices) || throw(
        BoundsError(1:final_dimension, indices),
    )
    return indices
end

function _gto_block_indices(final_dimension::Int, block_indices)
    throw(
        ArgumentError(
            "block_indices must be nothing or an explicit vector of basis indices; got $(typeof(block_indices))",
        ),
    )
end

function _gto_weight_vector(weights::Symbol, norbitals::Int)
    weights == :uniform && return ones(Float64, norbitals)
    weights == :shell_equalized && throw(
        ArgumentError(
            "gto_occupancy_matrix weights = :shell_equalized is not implemented yet; pass :uniform or an explicit weight vector",
        ),
    )
    throw(
        ArgumentError(
            "unsupported gto_occupancy_matrix weights :$(weights); supported first-pass modes are :uniform and explicit vectors",
        ),
    )
end

function _gto_weight_vector(weights::AbstractVector{<:Real}, norbitals::Int)
    length(weights) == norbitals || throw(
        DimensionMismatch(
            "explicit GTO weights length $(length(weights)) does not match probe dimension $(norbitals)",
        ),
    )
    values = Float64[Float64(weight) for weight in weights]
    all(isfinite, values) || throw(ArgumentError("explicit GTO weights must be finite"))
    return values
end

function _gto_weight_vector(weights, norbitals::Int)
    throw(
        ArgumentError(
            "gto_occupancy_matrix weights must be :uniform or an explicit real vector; got $(typeof(weights))",
        ),
    )
end

function _cartesian_exact_cartesian_probe_cross(
    raw,
    probes::CartesianGaussianShellSupplementRepresentation3D,
)
    hasproperty(raw, :exact_cartesian_supplement_overlap) || return nothing
    exact_cross = raw.exact_cartesian_supplement_overlap
    exact_cross === nothing && return nothing
    _cartesian_same_supplement_raw_identity(raw.supplement_representation, probes) || return nothing
    return Matrix{Float64}(exact_cross)
end

function _cartesian_exact_supplement_probe_cross(
    raw,
    probes::CartesianGaussianShellSupplementRepresentation3D,
)
    hasproperty(raw, :exact_supplement_overlap) || return nothing
    exact_cross = raw.exact_supplement_overlap
    exact_cross === nothing && return nothing
    _cartesian_same_supplement_raw_identity(raw.supplement_representation, probes) || return nothing
    return Matrix{Float64}(exact_cross)
end

function _cartesian_cartesian_probe_overlap(
    raw,
    probes::CartesianGaussianShellSupplementRepresentation3D,
)
    isempty(probes.orbitals) &&
        return zeros(Float64, raw.cartesian_representation.metadata.final_dimension, 0)
    exact_cross = _cartesian_exact_cartesian_probe_cross(raw, probes)
    exact_cross !== nothing && return exact_cross
    try
        return _cartesian_basis_supplement_cross(raw.cartesian_representation, probes)
    catch error
        if error isa ArgumentError
            return _cartesian_factorized_basis_supplement_cross(
                raw.factorized_cartesian_parent_basis,
                raw.cartesian_representation,
                probes,
                raw,
            )
        end
        rethrow()
    end
end

function _cartesian_supplement_probe_overlap(
    raw,
    probes::CartesianGaussianShellSupplementRepresentation3D,
)
    isempty(raw.supplement_representation.orbitals) &&
        return zeros(Float64, 0, length(probes.orbitals))
    exact_cross = _cartesian_exact_supplement_probe_cross(raw, probes)
    exact_cross !== nothing && return exact_cross
    return _cartesian_supplement_cross_overlap(raw.supplement_representation, probes)
end

function _gto_overlap_matrix(
    working::CartesianBasisRepresentation3D,
    probes::CartesianGaussianShellSupplementRepresentation3D,
)
    raw = _cartesian_raw_components(working)
    cg = _cartesian_cartesian_probe_overlap(raw, probes)
    gg = _cartesian_supplement_probe_overlap(raw, probes)
    raw_overlap = [cg; gg]
    return Matrix{Float64}(transpose(raw.raw_to_final) * raw_overlap)
end

"""
    gto_overlap_matrix(working, probes; block_indices = nothing)

Return the exact overlap matrix `S_BG = <B|G>` between a supported Cartesian
working basis `B` and a legacy Gaussian probe family `G`.

`working` may be a `CartesianBasisRepresentation3D`, a supported public
Cartesian basis object such as a bond-aligned ordinary basis or nested fixed
block, a `MappedUniformBasis` interpreted as its 3D Cartesian direct product, or
an `OrdinaryCartesianOperators3D` whose public representation carries the exact
hybrid Cartesian/supplement sidecars. `probes` may be a legacy Gaussian
supplement or a `CartesianGaussianShellSupplementRepresentation3D`.

When `block_indices` is supplied, only those working-basis rows are returned.
Block policy is intentionally external: this helper does not choose octants,
shells, DG regions, or any other geometry partition.
"""
function gto_overlap_matrix(working, probes; block_indices = nothing)
    working_representation = _gto_working_representation(working)
    probe_representation = _gto_probe_representation(probes)
    overlap = _gto_overlap_matrix(working_representation, probe_representation)
    indices = _gto_block_indices(working_representation.metadata.final_dimension, block_indices)
    indices === nothing && return overlap
    return Matrix{Float64}(overlap[indices, :])
end

function gto_overlap_matrix(
    working,
    probes,
    block_indices::AbstractVector{<:Integer},
)
    return gto_overlap_matrix(working, probes; block_indices = block_indices)
end

"""
    gto_occupancy_matrix(working, probes; weights = :uniform, block_indices = nothing)

Build a density-like GTO importance operator for basis-design work:

`M_I = S_{I,G} W S_{G,I}`

where `S_{I,G}` is the exact block-restricted overlap from
`gto_overlap_matrix(...)`. The result is an RDM-like design/occupancy matrix,
not a physical one-particle density matrix and not an SCF object.

First-pass weighting supports `weights = :uniform` and explicit real weight
vectors with one entry per probe orbital. `weights = :shell_equalized` is a
reserved future option and is intentionally not implemented yet.
"""
function gto_occupancy_matrix(working, probes; weights = :uniform, block_indices = nothing)
    overlap = gto_overlap_matrix(working, probes; block_indices = block_indices)
    weight_vector = _gto_weight_vector(weights, size(overlap, 2))
    weighted_overlap = overlap .* reshape(weight_vector, 1, :)
    return Matrix{Float64}(weighted_overlap * transpose(overlap))
end

function gto_occupancy_matrix(
    working,
    probes,
    block_indices::AbstractVector{<:Integer};
    weights = :uniform,
)
    return gto_occupancy_matrix(
        working,
        probes;
        weights = weights,
        block_indices = block_indices,
    )
end

function _cartesian_empty_centers()
    return zeros(Float64, 0, 3)
end

function _cartesian_supplement_center_matrix(
    supplement::CartesianGaussianShellSupplementRepresentation3D,
)
    isempty(supplement.orbitals) && return _cartesian_empty_centers()
    matrix = Matrix{Float64}(undef, length(supplement.orbitals), 3)
    for (row, orbital) in pairs(supplement.orbitals)
        matrix[row, 1] = orbital.center[1]
        matrix[row, 2] = orbital.center[2]
        matrix[row, 3] = orbital.center[3]
    end
    return matrix
end

function _cartesian_supplement_orbital_representation(
    orbital::_AtomicCartesianShellOrbital3D,
)
    return CartesianGaussianShellOrbitalRepresentation3D(
        String(orbital.label),
        (Int(orbital.lx), Int(orbital.ly), Int(orbital.lz)),
        (
            Float64(orbital.center[1]),
            Float64(orbital.center[2]),
            Float64(orbital.center[3]),
        ),
        Float64[Float64(value) for value in orbital.exponents],
        Float64[Float64(value) for value in orbital.coefficients],
        :axiswise_normalized_cartesian_gaussian,
    )
end

function _cartesian_supplement_kind(::Nothing)
    return :none
end

function _cartesian_supplement_kind(::LegacyAtomicGaussianSupplement)
    return :atomic_cartesian_shell
end

function _cartesian_supplement_kind(::LegacyBondAlignedDiatomicGaussianSupplement)
    return :bond_aligned_diatomic_cartesian_shell
end

function _cartesian_supplement_kind(::LegacyBondAlignedHeteronuclearGaussianSupplement)
    return :bond_aligned_heteronuclear_cartesian_shell
end

function _cartesian_supplement_lmax(::Nothing)
    return nothing
end

function _cartesian_supplement_lmax(data::LegacyAtomicGaussianSupplement)
    return data.lmax
end

function _cartesian_supplement_lmax(data::LegacyBondAlignedDiatomicGaussianSupplement)
    return data.atomic_source.lmax
end

function _cartesian_supplement_lmax(data::LegacyBondAlignedHeteronuclearGaussianSupplement)
    return maximum(source.lmax for source in data.atomic_sources)
end

function _cartesian_supplement_metadata(::Nothing)
    return (
        source_kind = :none,
        atom = nothing,
        basis_name = nothing,
        basisfile = nothing,
        lmax = nothing,
        nuclei = NTuple{3,Float64}[],
        uncontracted = nothing,
        max_width = nothing,
    )
end

function _cartesian_supplement_metadata(data::LegacyAtomicGaussianSupplement)
    return (
        source_kind = :legacy_atomic_gaussian_supplement,
        atom = String(data.atom),
        basis_name = String(data.basis_name),
        basisfile = String(data.basisfile),
        lmax = data.lmax,
        nuclei = NTuple{3,Float64}[(0.0, 0.0, 0.0)],
        uncontracted = data.uncontracted,
        max_width = data.max_width,
    )
end

function _cartesian_supplement_metadata(data::LegacyBondAlignedDiatomicGaussianSupplement)
    return (
        source_kind = :legacy_bond_aligned_diatomic_gaussian_supplement,
        atom = String(data.atomic_source.atom),
        basis_name = String(data.atomic_source.basis_name),
        basisfile = String(data.atomic_source.basisfile),
        lmax = data.atomic_source.lmax,
        nuclei = NTuple{3,Float64}[
            (
                Float64(nucleus[1]),
                Float64(nucleus[2]),
                Float64(nucleus[3]),
            ) for nucleus in data.nuclei
        ],
        uncontracted = data.atomic_source.uncontracted,
        max_width = data.max_width,
    )
end

function _cartesian_supplement_metadata(data::LegacyBondAlignedHeteronuclearGaussianSupplement)
    return (
        source_kind = :legacy_bond_aligned_heteronuclear_gaussian_supplement,
        atom = nothing,
        basis_name = nothing,
        basisfile = nothing,
        lmax = maximum(source.lmax for source in data.atomic_sources),
        nuclei = NTuple{3,Float64}[
            (
                Float64(nucleus[1]),
                Float64(nucleus[2]),
                Float64(nucleus[3]),
            ) for nucleus in data.nuclei
        ],
        uncontracted = all(source.uncontracted for source in data.atomic_sources),
        max_width = data.max_width,
    )
end

function _cartesian_empty_supplement_representation()
    return CartesianGaussianShellSupplementRepresentation3D(
        :none,
        CartesianGaussianShellOrbitalRepresentation3D[],
        _cartesian_supplement_metadata(nothing),
    )
end

function _cartesian_supplement_representation(::Nothing)
    return _cartesian_empty_supplement_representation()
end

function _cartesian_supplement_representation(
    data::Union{
        LegacyAtomicGaussianSupplement,
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    },
)
    supplement3d =
        data isa LegacyAtomicGaussianSupplement ? _atomic_cartesian_shell_supplement_3d(data) :
        _bond_aligned_diatomic_cartesian_shell_supplement_3d(data)
    orbitals = CartesianGaussianShellOrbitalRepresentation3D[
        _cartesian_supplement_orbital_representation(orbital) for orbital in supplement3d.orbitals
    ]
    return CartesianGaussianShellSupplementRepresentation3D(
        _cartesian_supplement_kind(data),
        orbitals,
        _cartesian_supplement_metadata(data),
    )
end

function basis_representation(
    data::Union{
        LegacyAtomicGaussianSupplement,
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    },
)
    return _cartesian_supplement_representation(data)
end

function basis_metadata(
    data::Union{
        LegacyAtomicGaussianSupplement,
        LegacyBondAlignedDiatomicGaussianSupplement,
        LegacyBondAlignedHeteronuclearGaussianSupplement,
    },
)
    return basis_representation(data).metadata
end

function _cartesian_supplement_orbital_signature(
    orbital::CartesianGaussianShellOrbitalRepresentation3D,
)
    return (
        label = orbital.label,
        angular_powers = orbital.angular_powers,
        center = orbital.center,
        exponents = Tuple(Float64[Float64(value) for value in orbital.exponents]),
        coefficients = Tuple(Float64[Float64(value) for value in orbital.coefficients]),
        primitive_normalization = orbital.primitive_normalization,
    )
end

function _cartesian_same_supplement_raw_identity(
    left::CartesianGaussianShellSupplementRepresentation3D,
    right::CartesianGaussianShellSupplementRepresentation3D,
)
    left.supplement_kind == right.supplement_kind || return false
    isequal(left.metadata, right.metadata) || return false
    length(left.orbitals) == length(right.orbitals) || return false
    for index in eachindex(left.orbitals)
        _cartesian_supplement_orbital_signature(left.orbitals[index]) ==
        _cartesian_supplement_orbital_signature(right.orbitals[index]) || return false
    end
    return true
end
